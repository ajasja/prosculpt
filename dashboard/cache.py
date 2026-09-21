"""Small in-process cache for the dashboard's expensive filesystem scans.

Cache only values whose staleness is harmless, such as progress counts. Stage,
crash, cancellation and terminal verdicts must never be cached.

Entries expire on a time-to-live measured with time.monotonic(); there is no
mtime-based invalidation.
"""

import os
import threading
import time
from collections import OrderedDict

DEFAULT_TTL_SECONDS = 15.0

# A scan is cached for this multiple of the time it took to produce, so that a
# slow scan stays valid across a whole get_job_status pass, capped at
# MAX_TTL_SECONDS.
SLOW_SCAN_TTL_FACTOR = 20.0
MAX_TTL_SECONDS = 900.0
MAX_ENTRIES = 512

# key -> (stored_at, value, ttl_for_this_entry)
_entries: "OrderedDict[tuple, tuple[float, object, float]]" = OrderedDict()
_dict_lock = threading.Lock()
_key_locks: dict = {}


def _norm(path: str) -> str:
    """Normalised cache key for a path, so that spellings of the same
    directory differing in case or separators share one entry."""
    try:
        return os.path.normcase(os.path.realpath(path))
    except OSError:
        return os.path.normcase(os.path.abspath(path))


def get_or_compute(namespace: str, path: str, compute, ttl: float = DEFAULT_TTL_SECONDS):
    """Cached `compute()` for (namespace, path).

    Concurrent callers for the same key are single-flighted: the second one
    waits for the first rather than launching its own scan of the same
    directory.
    """
    key = (namespace, _norm(path))
    with _dict_lock:
        hit = _entries.get(key)
        if hit is not None and time.monotonic() - hit[0] < hit[2]:
            _entries.move_to_end(key)
            return hit[1]
        lock = _key_locks.get(key)
        if lock is None:
            lock = _key_locks[key] = threading.Lock()

    # compute() runs outside _dict_lock; only callers waiting on this same
    # key are blocked.
    with lock:
        with _dict_lock:
            hit = _entries.get(key)
            if hit is not None and time.monotonic() - hit[0] < hit[2]:
                _entries.move_to_end(key)
                return hit[1]
        started = time.monotonic()
        value = compute()
        elapsed = time.monotonic() - started
        effective_ttl = max(ttl, min(elapsed * SLOW_SCAN_TTL_FACTOR, MAX_TTL_SECONDS))
        with _dict_lock:
            _entries[key] = (time.monotonic(), value, effective_ttl)
            _entries.move_to_end(key)
            while len(_entries) > MAX_ENTRIES:
                evicted, _ = _entries.popitem(last=False)
                _key_locks.pop(evicted, None)
        return value


def clear() -> None:
    """Drop every entry."""
    with _dict_lock:
        _entries.clear()
        _key_locks.clear()


def stats() -> dict:
    with _dict_lock:
        return {
            "entries": len(_entries),
            "max_entries": MAX_ENTRIES,
            "default_ttl_seconds": DEFAULT_TTL_SECONDS,
        }
