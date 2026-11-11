"""Shared affine lift engine and the six larger finite constructions."""
from pathlib import Path
from .affine import load, verify_compact, squares

ROOT = Path(__file__).resolve().parents[1]

def construct(order):
    cert = load(order, ROOT / "certificates" / "affine" / f"n{order}.json")
    verify_compact(cert)
    return squares(cert)[0]
