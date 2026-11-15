#!/usr/bin/env python3
"""Plot DropAlpha strong-scaling and problem-size benchmark CSV files."""

import argparse
import csv
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter


def read_csv(path: Path, required_columns: set[str]) -> list[dict[str, float]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        columns = set(reader.fieldnames or [])
        missing = required_columns - columns
        if missing:
            raise ValueError(f"{path} is missing columns: {', '.join(sorted(missing))}")

        rows = []
        for line_number, row in enumerate(reader, start=2):
            try:
                rows.append({key: float(value) for key, value in row.items()})
            except (TypeError, ValueError) as exc:
                raise ValueError(f"Invalid numeric value in {path}:{line_number}") from exc
    return rows


def newest_results_dir(root: Path) -> Path:
    candidates = [
        path
        for path in root.iterdir()
        if path.is_dir()
        and (
            (path / "strong_scaling.csv").is_file()
            or (path / "problem_size_scaling.csv").is_file()
        )
    ]
    if not candidates:
        raise FileNotFoundError(f"No benchmark result directories found under {root}")
    return max(candidates, key=lambda path: path.stat().st_mtime)


def save_figure(fig: plt.Figure, output_base: Path) -> None:
    output_base.parent.mkdir(parents=True, exist_ok=True)
    png_path = output_base.with_suffix(".png")
    pdf_path = output_base.with_suffix(".pdf")
    fig.savefig(png_path, dpi=220, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {png_path}")
    print(f"Wrote {pdf_path}")


def format_unknowns_millions(value: float, _position: int) -> str:
    label = f"{value / 1.0e6:.2f}".rstrip("0").rstrip(".")
    return label


def plot_strong_scaling(rows: list[dict[str, float]], output_dir: Path) -> None:
    rows.sort(key=lambda row: row["ranks"])
    ranks = [row["ranks"] for row in rows]
    elapsed = [row["elapsed_s"] for row in rows]
    speedup = [row["speedup"] for row in rows]
    efficiency = [100.0 * row["efficiency"] for row in rows]
    baseline_rank = ranks[0]
    ideal_speedup = [rank / baseline_rank for rank in ranks]

    fig, axes = plt.subplots(1, 3, figsize=(12.0, 3.6))

    axes[0].plot(ranks, elapsed, "o-", color="#2f6df6", linewidth=1.8)
    axes[0].set_xlabel("MPI ranks")
    axes[0].set_ylabel("Elapsed time [s]")
    axes[0].set_title("Runtime")

    axes[1].plot(ranks, speedup, "o-", color="#44a35f", linewidth=1.8, label="Measured")
    axes[1].plot(ranks, ideal_speedup, "--", color="black", linewidth=1.2, label="Ideal")
    axes[1].set_xlabel("MPI ranks")
    axes[1].set_ylabel("Speedup")
    axes[1].set_title("Strong scaling")
    axes[1].legend(frameon=False)

    axes[2].plot(ranks, efficiency, "o-", color="#d95f02", linewidth=1.8)
    axes[2].axhline(100.0, color="black", linestyle="--", linewidth=1.2)
    axes[2].set_xlabel("MPI ranks")
    axes[2].set_ylabel("Parallel efficiency [%]")
    axes[2].set_title("Efficiency")

    for axis in axes:
        axis.grid(True, alpha=0.25, linewidth=0.6)
        axis.set_xticks(ranks)

    cells = int(rows[0]["cells"])
    unknowns = int(rows[0]["unknowns"])
    fig.suptitle(f"DropAlpha strong scaling: {cells:,} cells, {unknowns:,} unknowns")
    fig.tight_layout()
    save_figure(fig, output_dir / "strong_scaling")


def plot_problem_size(rows: list[dict[str, float]], output_dir: Path) -> None:
    rows.sort(key=lambda row: row["unknowns"])
    unknowns = [row["unknowns"] for row in rows]
    elapsed = [row["elapsed_s"] for row in rows]
    ranks = sorted({int(row["ranks"]) for row in rows})

    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.6))

    axes[0].loglog(unknowns, elapsed, "o-", color="#2f6df6", linewidth=1.8)
    axes[0].set_xlabel("Block-system unknowns [$10^6$]")
    axes[0].set_ylabel("Elapsed time [s]")
    axes[0].set_title("Problem-size scaling")

    time_per_unknown = [time / count for time, count in zip(elapsed, unknowns)]
    axes[1].loglog(unknowns, time_per_unknown, "o-", color="#d95f02", linewidth=1.8)
    axes[1].set_xlabel("Block-system unknowns [$10^6$]")
    axes[1].set_ylabel("Elapsed time / unknown [s]")
    axes[1].set_title("Normalized cost")

    for axis in axes:
        axis.grid(True, which="both", alpha=0.25, linewidth=0.6)
        axis.xaxis.set_major_locator(FixedLocator(unknowns))
        axis.xaxis.set_major_formatter(FuncFormatter(format_unknowns_millions))
        axis.xaxis.set_minor_formatter(NullFormatter())
        axis.tick_params(axis="x", labelsize=9)

    rank_text = ", ".join(str(rank) for rank in ranks)
    fig.suptitle(f"DropAlpha problem-size scaling: MPI ranks = {rank_text}")
    fig.tight_layout()
    save_figure(fig, output_dir / "problem_size_scaling")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot CSV output generated by runPerformanceScaling.sh."
    )
    parser.add_argument(
        "results_dir",
        nargs="?",
        type=Path,
        help="Benchmark directory containing the CSV files. Defaults to the newest run.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="Output directory. Defaults to <results_dir>/figures.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        results_dir = args.results_dir
        if results_dir is None:
            results_dir = newest_results_dir(Path("performanceResults"))
        results_dir = results_dir.resolve()
        output_dir = (args.output_dir or results_dir / "figures").resolve()

        generated = False
        strong_csv = results_dir / "strong_scaling.csv"
        if strong_csv.is_file():
            rows = read_csv(
                strong_csv,
                {"cells", "unknowns", "ranks", "elapsed_s", "speedup", "efficiency"},
            )
            if rows:
                plot_strong_scaling(rows, output_dir)
                generated = True

        size_csv = results_dir / "problem_size_scaling.csv"
        if size_csv.is_file():
            rows = read_csv(size_csv, {"unknowns", "ranks", "elapsed_s"})
            if rows:
                plot_problem_size(rows, output_dir)
                generated = True

        if not generated:
            raise FileNotFoundError(f"No non-empty performance CSV files found in {results_dir}")
    except (FileNotFoundError, ValueError, OSError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
