#!/usr/bin/env python3
"""Construct the catalogued MOLS from the supplied formulas and certificates.

Examples:
    python construct.py all --output regenerated
    python construct.py 22 90 144 --output regenerated
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from time import perf_counter

ROOT = Path(__file__).resolve().parent


def catalog() -> list[dict]:
    return json.loads((ROOT / "catalog.json").read_text(encoding="utf-8"))


def construct(case: str | int) -> list[list[list[int]]]:
    """Return k nested n-by-n squares with labels 0,...,n-1."""
    by_id = {row["id"]: row for row in catalog()}
    try:
        row = by_id[str(case)]
    except KeyError as exc:
        raise ValueError(f"Unknown construction {case}; choices: {', '.join(by_id)}") from exc
    n = row["order"]
    if row["module"] == "base":
        from mols_base.construct import construct as build
        squares = build(n)
    elif row["module"] == "finite":
        from construction import construct as build
        squares = build(n)
    else:
        from mols_geometry import construct as build
        squares = build(n)
    if len(squares) != row["mols"] or any(len(square) != n for square in squares):
        raise ValueError(f"Constructor returned incorrect dimensions for {case}")
    return squares


def write_squares(path: Path, squares: list[list[list[int]]]) -> str:
    """Write one deterministic ASCII file containing the whole set."""
    path.parent.mkdir(parents=True, exist_ok=True)
    n, k = len(squares[0]), len(squares)
    with path.open("w", encoding="ascii", newline="\n") as stream:
        stream.write(f"{n} {k}\n")
        for j, square in enumerate(squares):
            stream.write(f"# Square {j}\n")
            for row in square:
                stream.write(" ".join(map(str, row))+"\n")
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cases", nargs="+", help="Order(s), or all")
    parser.add_argument("--output", type=Path, default=ROOT / "regenerated")
    parser.add_argument("--report", type=Path, help="Write a generation JSON report")
    args = parser.parse_args()
    cases = catalog()
    if args.cases != ["all"]:
        chosen = set(args.cases)
        unknown = chosen - {row["id"] for row in cases}
        if unknown:
            parser.error(f"Unknown construction(s): {', '.join(sorted(unknown))}")
        cases = [row for row in cases if row["id"] in chosen]
    reports = []
    for row in cases:
        started = perf_counter()
        squares = construct(row["id"])
        path = args.output / row["filename"]
        digest = write_squares(path, squares)
        report = {"id": row["id"], "order": row["order"], "mols": len(squares),
                  "file": path.name, "sha256": digest,
                  "seconds": round(perf_counter()-started, 6)}
        reports.append(report)
        print(json.dumps(report), flush=True)
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps({"constructed": reports}, indent=2)+"\n")


if __name__ == "__main__":
    main()
