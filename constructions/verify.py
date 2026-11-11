#!/usr/bin/env python3
"""Independent, standard-library verification of the final square text files.

This module imports no construction code and uses no mathematical certificates.
Format: first non-comment line is ``n k``; then k blocks of n lines of n
space-separated integers in 0,...,n-1. Empty lines and # comments are ignored.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from time import perf_counter


def read_squares(path: str | Path) -> tuple[int, list[list[list[int]]]]:
    """Read the exact text contract, rejecting missing or additional data."""
    with Path(path).open("r", encoding="ascii") as stream:
        records = ((line_number, line.strip())
                   for line_number, line in enumerate(stream, 1)
                   if line.strip() and not line.lstrip().startswith("#"))
        try:
            line_number, header = next(records)
        except StopIteration as exc:
            raise ValueError("Empty square file") from exc
        try:
            n, k = map(int, header.split())
        except ValueError as exc:
            raise ValueError("Header must contain exactly two integers: n k") from exc
        if n < 2 or not 1 <= k < n:
            raise ValueError("Require n >= 2 and 1 <= k < n")
        squares = []
        for square_index in range(k):
            square = []
            for row_index in range(n):
                try:
                    line_number, line = next(records)
                except StopIteration as exc:
                    raise ValueError(f"Missing square {square_index}, row {row_index}") from exc
                try:
                    row = list(map(int, line.split()))
                except ValueError as exc:
                    raise ValueError(f"Non-integer at line {line_number}") from exc
                if len(row) != n:
                    raise ValueError(f"Line {line_number}: expected {n} entries")
                if any(value < 0 or value >= n for value in row):
                    raise ValueError(f"Line {line_number}: symbol outside 0..{n-1}")
                square.append(row)
            squares.append(square)
        if next(records, None) is not None:
            raise ValueError("Unexpected trailing data")
    return n, squares


def verify_squares(squares: list[list[list[int]]]) -> dict:
    """Exhaustively check dimensions, symbols, Latinness, and every square pair."""
    if not squares or not squares[0]:
        raise ValueError("At least one nonempty square is required")
    n, k = len(squares[0]), len(squares)
    if n < 2 or not 1 <= k < n:
        raise ValueError("Require n >= 2 and 1 <= k < n")
    for j, square in enumerate(squares):
        if len(square) != n or any(len(row) != n for row in square):
            raise ValueError(f"Square {j}: inconsistent dimensions")
        for row in square:
            if any(type(value) is not int or not 0 <= value < n for value in row):
                raise ValueError(f"Square {j}: invalid symbol")
        for r, row in enumerate(square):
            if len(set(row)) != n:
                raise ValueError(f"Square {j}, row {r}: repeated symbol")
        for c in range(n):
            if len({square[r][c] for r in range(n)}) != n:
                raise ValueError(f"Square {j}, column {c}: repeated symbol")
    for i in range(k):
        for j in range(i+1, k):
            seen = bytearray(n*n)
            for r in range(n):
                for c in range(n):
                    pair = n*squares[i][r][c]+squares[j][r][c]
                    if seen[pair]:
                        raise ValueError(f"Squares {i},{j}: repeated pair at ({r},{c})")
                    seen[pair] = 1
            # n^2 distinct in-range pairs exhaust the entire pair alphabet.
    return {"order": n, "mols": k, "rows_checked": n*k,
            "columns_checked": n*k, "square_pairs_checked": k*(k-1)//2,
            "ordered_pairs_checked": n*n*k*(k-1)//2, "valid": True}


def verify_file(path: str | Path) -> dict:
    path = Path(path)
    started = perf_counter()
    _, squares = read_squares(path)
    report = verify_squares(squares)
    report.update(file=path.name,
                  sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                  seconds=round(perf_counter()-started, 6))
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", nargs="*", type=Path,
                        help="Square files or directories (default: squares/)")
    parser.add_argument("--report", type=Path, help="Write the full JSON report")
    args = parser.parse_args()
    paths = args.paths or [Path(__file__).resolve().parent / "squares"]
    files = []
    for path in paths:
        files.extend(sorted(path.glob("mols_*.txt")) if path.is_dir() else [path])
    if not files:
        parser.error("No square files found")
    reports = []
    for path in files:
        report = verify_file(path)
        reports.append(report)
        print(json.dumps(report), flush=True)
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps({"valid": True, "files": reports}, indent=2)+"\n")


if __name__ == "__main__":
    main()
