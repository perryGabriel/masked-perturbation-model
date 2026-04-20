from __future__ import annotations

import concurrent.futures
from dataclasses import dataclass
from typing import Any, Callable


@dataclass
class TimeoutOutcome:
    timed_out: bool
    value: Any = None
    error: str | None = None


def run_with_timeout(fn: Callable[..., Any], timeout_sec: float, *args, **kwargs) -> TimeoutOutcome:
    with concurrent.futures.ProcessPoolExecutor(max_workers=1) as ex:
        fut = ex.submit(fn, *args, **kwargs)
        try:
            return TimeoutOutcome(timed_out=False, value=fut.result(timeout=timeout_sec), error=None)
        except concurrent.futures.TimeoutError:
            fut.cancel()
            ex.shutdown(wait=False, cancel_futures=True)
            return TimeoutOutcome(timed_out=True, value=None, error="timeout")
        except Exception as exc:
            return TimeoutOutcome(timed_out=False, value=None, error=f"{type(exc).__name__}: {exc}")
