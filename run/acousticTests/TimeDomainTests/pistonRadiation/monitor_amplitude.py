#!/usr/bin/env python3
import argparse
import os
import re
import signal
import subprocess
import sys
import time
from typing import List, Tuple


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Monitor probe pressure amplitude and stop when stable.")
    p.add_argument("--case", default=".", help="Case directory")
    p.add_argument("--field", default="p", help="Field to read (p or pr)")
    p.add_argument("--frequency", type=float, required=True, help="Forcing frequency [Hz]")
    p.add_argument("--cycles", type=float, default=5.0, help="Cycles per window")
    p.add_argument("--tol", type=float, default=1e-3, help="Relative amplitude tolerance")
    p.add_argument("--min-samples", type=int, default=50, help="Minimum samples per window")
    p.add_argument("--check-interval", type=float, default=2.0, help="Seconds between checks")
    p.add_argument("--pid", type=int, default=None, help="Solver PID to watch")
    p.add_argument("--endtime-buffer", type=float, default=0.0, help="Extra time added when stopping")
    return p.parse_args()


def latest_probe_file(case_dir: str, field: str) -> str:
    probes_dir = os.path.join(case_dir, "postProcessing", "probes")
    if not os.path.isdir(probes_dir):
        return ""
    subdirs = [d for d in os.listdir(probes_dir) if os.path.isdir(os.path.join(probes_dir, d))]
    if not subdirs:
        return ""
    def to_float(x: str) -> float:
        try:
            return float(x)
        except ValueError:
            return -1.0
    subdirs.sort(key=to_float)
    for d in reversed(subdirs):
        fpath = os.path.join(probes_dir, d, field)
        if os.path.isfile(fpath):
            return fpath
    return ""


def read_probe_timeseries(path: str) -> Tuple[List[float], List[float]]:
    times: List[float] = []
    values: List[float] = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) < 2:
                continue
            try:
                t = float(parts[0])
                v = float(parts[1])
            except ValueError:
                continue
            times.append(t)
            values.append(v)
    return times, values


def window_indices(times: List[float], t_start: float, t_end: float) -> List[int]:
    idx = []
    for i, t in enumerate(times):
        if t_start <= t <= t_end:
            idx.append(i)
    return idx


def amplitude(values: List[float], idx: List[int]) -> float:
    if not idx:
        return 0.0
    vmin = min(values[i] for i in idx)
    vmax = max(values[i] for i in idx)
    return 0.5 * (vmax - vmin)


def solver_alive(pid: int) -> bool:
    if pid is None:
        return True
    try:
        os.kill(pid, 0)
    except OSError:
        return False
    return True


def stop_solver(case_dir: str, stop_time: float) -> None:
    cmd = ["foamDictionary", "-entry", "endTime", "-set", f"{stop_time}", os.path.join(case_dir, "system", "controlDict")]
    subprocess.run(cmd, check=False, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def main() -> int:
    args = parse_args()
    case_dir = os.path.abspath(args.case)

    period = 1.0 / args.frequency
    window = args.cycles * period

    last_amp = None
    while solver_alive(args.pid):
        path = latest_probe_file(case_dir, args.field)
        if not path:
            time.sleep(args.check_interval)
            continue

        times, values = read_probe_timeseries(path)
        if not times:
            time.sleep(args.check_interval)
            continue

        t_end = times[-1]
        t_mid = t_end - window
        t_start = t_end - 2.0 * window

        if t_start < times[0]:
            time.sleep(args.check_interval)
            continue

        idx_prev = window_indices(times, t_start, t_mid)
        idx_curr = window_indices(times, t_mid, t_end)

        if len(idx_prev) < args.min_samples or len(idx_curr) < args.min_samples:
            time.sleep(args.check_interval)
            continue

        amp_prev = amplitude(values, idx_prev)
        amp_curr = amplitude(values, idx_curr)
        denom = amp_prev if amp_prev != 0.0 else max(amp_curr, 1.0)
        rel = abs(amp_curr - amp_prev) / denom

        if last_amp is None:
            last_amp = amp_curr
        else:
            last_amp = amp_curr

        if rel <= args.tol:
            stop_time = t_end + args.endtime_buffer
            stop_solver(case_dir, stop_time)
            return 0

        time.sleep(args.check_interval)

    return 0


if __name__ == "__main__":
    sys.exit(main())
