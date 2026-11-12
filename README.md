# Data and code for "Improved Lower Bounds for Mutually Orthogonal Latin Squares (MOLS) at 33 Orders"

This repository has the data and code for:

C.D. Rosin. (November 2026) “Improved Lower Bounds for Mutually Orthogonal Latin Squares (MOLS) at 33 Orders.”

The directory **mols/** has the new sets of MOLS.  Decompress the files first (e.g. by running `gunzip *` in the directory).  Files are named mols_n_m.txt, giving m MOLS of order n.  The first line of the file has n and m.  The remainder of the file has the Latin squares.  Ignore comment lines which start with \#.

To verify the MOLS using the (author-created) Python verification script, cd into the **mols/** directory and run:
`python3 ../checkmols.py`

The construction code (created by GPT-6 Astra) is in the directory **constructions/**.  See the **README.md** file there.

**srg_484_126_40_30.txt** is the strongly regular graph SRG(484,126,40,30) generated from 4 MOLS of order 22.

