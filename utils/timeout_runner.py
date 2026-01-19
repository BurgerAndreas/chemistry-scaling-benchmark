"""
Timeout runner with memory tracking for calculations.
"""

import multiprocessing as mp
import psutil
import time
from typing import Callable, Dict, Any, Tuple
import traceback


def _memory_monitor(
    pid: int, memory_queue: mp.Queue, stop_event: mp.Event, sample_interval: float = 0.1
):
    """
    Monitor memory usage of a process.

    Args:
        pid: Process ID to monitor
        memory_queue: Queue to send peak memory usage
        stop_event: Event to signal when to stop monitoring
        sample_interval: How often to sample memory (seconds)
    """
    peak_memory_mb = 0.0
    try:
        process = psutil.Process(pid)
        while not stop_event.is_set():
            try:
                mem_info = process.memory_info()
                memory_mb = mem_info.rss / 1024 / 1024  # Convert to MB
                peak_memory_mb = max(peak_memory_mb, memory_mb)
                time.sleep(sample_interval)
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                break
    except Exception:
        pass
    finally:
        memory_queue.put(peak_memory_mb)


def _run_calculation_worker(
    calc_func: Callable, args: Tuple, kwargs: Dict, result_queue: mp.Queue
):
    """
    Worker function to run calculation in separate process.

    Args:
        calc_func: Function to execute
        args: Positional arguments for calc_func
        kwargs: Keyword arguments for calc_func
        result_queue: Queue to send result
    """
    try:
        result = calc_func(*args, **kwargs)
        result_queue.put(("success", result))
    except Exception as e:
        error_msg = f"{type(e).__name__}: {str(e)}"
        tb = traceback.format_exc()
        result_queue.put(("error", error_msg, tb))


def run_with_timeout(
    calc_func: Callable, timeout_sec: float, *args, **kwargs
) -> Dict[str, Any]:
    """
    Run a calculation function with timeout and memory tracking.

    Args:
        calc_func: Function to execute
        timeout_sec: Maximum time to allow in seconds
        *args: Positional arguments for calc_func
        **kwargs: Keyword arguments for calc_func

    Returns:
        Dictionary with keys:
        - success: bool
        - time_seconds: float (wall-clock time)
        - peak_memory_mb: float
        - result: Any (if success=True)
        - error_message: str (if success=False)
        - traceback: str (if success=False and verbose error available)
    """
    # Create queues and events for IPC
    result_queue = mp.Queue()
    memory_queue = mp.Queue()
    stop_event = mp.Event()

    # Start calculation process
    calc_process = mp.Process(
        target=_run_calculation_worker, args=(calc_func, args, kwargs, result_queue)
    )

    start_time = time.time()
    calc_process.start()

    # Start memory monitor
    monitor_process = mp.Process(
        target=_memory_monitor, args=(calc_process.pid, memory_queue, stop_event)
    )
    monitor_process.start()

    # Wait for calculation to complete or timeout
    calc_process.join(timeout=timeout_sec)
    elapsed_time = time.time() - start_time

    # Check if process completed
    if calc_process.is_alive():
        # Timeout occurred
        calc_process.terminate()
        calc_process.join(timeout=5)
        if calc_process.is_alive():
            calc_process.kill()
            calc_process.join()

        stop_event.set()
        monitor_process.join(timeout=2)
        if monitor_process.is_alive():
            monitor_process.terminate()
            monitor_process.join()

        # Get peak memory
        peak_memory_mb = 0.0
        if not memory_queue.empty():
            peak_memory_mb = memory_queue.get()

        return {
            "success": False,
            "time_seconds": timeout_sec,
            "peak_memory_mb": peak_memory_mb,
            "error_message": "TIMEOUT",
        }

    # Process completed, stop monitor
    stop_event.set()
    monitor_process.join(timeout=2)
    if monitor_process.is_alive():
        monitor_process.terminate()
        monitor_process.join()

    # Get peak memory
    peak_memory_mb = 0.0
    if not memory_queue.empty():
        peak_memory_mb = memory_queue.get()

    # Get result
    if not result_queue.empty():
        result_data = result_queue.get()
        if result_data[0] == "success":
            return {
                "success": True,
                "time_seconds": elapsed_time,
                "peak_memory_mb": peak_memory_mb,
                "result": result_data[1],
            }
        else:
            # Error occurred
            error_msg = result_data[1]
            tb = result_data[2] if len(result_data) > 2 else ""
            return {
                "success": False,
                "time_seconds": elapsed_time,
                "peak_memory_mb": peak_memory_mb,
                "error_message": error_msg,
                "traceback": tb,
            }
    else:
        # No result available (unexpected)
        return {
            "success": False,
            "time_seconds": elapsed_time,
            "peak_memory_mb": peak_memory_mb,
            "error_message": "UNKNOWN_ERROR",
        }
