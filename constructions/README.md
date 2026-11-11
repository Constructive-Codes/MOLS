# MOLS constructions

Use **Python 3.10 or later**. The Python code requires only the standard library.
Run the following commands from this directory:

```bash
python construct.py all --output regenerated
python verify.py regenerated --report verification.json
```

The first command constructs every catalogued set from its formulas and compact
numerical inputs, in the directory `regenerated/`. It does not read the completed squares.
It writes one file per order and prints a JSON summary containing its dimensions, filename,
SHA-256 digest, and construction time. The second command checks the generated
files independently: dimensions and symbol ranges, every row and column, and
every pair of squares.

## Results

`catalog.json` identifies each order, its number of squares, its construction
family, its constructor, and its output filename. Its `benchmark_2024` field
records the comparison value from Miller, Abel, Valkov, and Fraser,
*Implementing the MOLS Table for n Up to 500*, Symmetry 16(12), 1678 (2024),
https://doi.org/10.3390/sym16121678.

| Order | MOLS | Order | MOLS | Order | MOLS | Order | MOLS |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 22 | 4 | 85 | 7 | 183 | 7 | 430 | 11 |
| 45 | 7 | 90 | 8 | 192 | 9 | 442 | 11 |
| 50 | 7 | 92 | 7 | 207 | 9 | 456 | 11 |
| 54 | 9 | 96 | 11 | 288 | 17 | 469 | 11 |
| 72 | 8 | 98 | 7 | 324 | 11 | 482 | 7 |
| 74 | 6 | 144 | 14 | 328 | 9 | 486 | 9 |
| 78 | 7 | 162 | 9 | 335 | 9 | 494 | 11 |
| 80 | 10 | 175 | 7 | 427 | 12 |  |  |

## Package contents

| Path | Purpose |
|---|---|
| `construct.py` | Construction entry point for all 31 orders |
| `catalog.json` | Construction catalogue and output filenames |
| `construction/` | Shared affine construction engine and its entry point |
| `certificates/affine/` | Numerical inputs for six affine constructions |
| `mols_base/` | DM, QDM, relative-DM, and affine-array constructions and their inputs |
| `mols_geometry/` | Finite fields, plane incidence, point selections, and block filling |
| `verify.py` | Independent verifier; it imports no construction modules |

## Square-file format

Each file contains one complete set of squares. The first line is `n k`, where
`n` is the order and `k` is the number of squares. For each square
`j=0,...,k-1`, a comment `# Square j` is followed by exactly `n` lines, each
containing `n` space-separated integers in `0,...,n-1`. Files use ASCII and LF
line endings.

The verifier ignores blank lines and lines beginning with `#`; comments do not
serve as row or symbol labels. It rejects missing or extra data, invalid symbols,
repeated symbols in a row or column, and repeated ordered pairs in the
superposition of any two squares.

