#!/usr/bin/env python3
"""Regenerate a complete ICOOR_INTERNAL block for a proline-like ring imino-acid
param whose generator left the ICOOR truncated (ring-closure virtual shadow
collapsed onto a ring atom -> d=0). Geometry is computed from the generator's
own 3D PDB; the shadow atom (XD/NV) is repositioned onto the backbone N.

Specific to the azetidine-2-carboxylic-acid (02A) atom set:
  backbone  N CA C O ; ring CB CG (CG closes to N) ; XD = virtual shadow of N
  H: HA 1HB 2HB 1HG 2HG ; connects LOWER(on N) UPPER(on C)
Usage: fix_ring_icoor.py <params_in> <pdb_in> <params_out>
"""
import sys, math

def vsub(a, b): return [a[i]-b[i] for i in range(3)]
def vlen(a):    return math.sqrt(sum(x*x for x in a))
def dot(a, b):  return sum(a[i]*b[i] for i in range(3))
def cross(a, b):
    return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]

def bond_angle(p1, p2, p3):
    """angle at p2 in degrees"""
    v1, v2 = vsub(p1, p2), vsub(p3, p2)
    c = dot(v1, v2)/(vlen(v1)*vlen(v2))
    c = max(-1.0, min(1.0, c))
    return math.degrees(math.acos(c))

def dihedral(p1, p2, p3, p4):
    b1, b2, b3 = vsub(p2, p1), vsub(p3, p2), vsub(p4, p3)
    n1, n2 = cross(b1, b2), cross(b2, b3)
    m = cross(n1, [x/vlen(b2) for x in b2])
    x, y = dot(n1, n2), dot(m, n2)
    return math.degrees(math.atan2(y, x))

params_in, pdb_in, params_out = sys.argv[1], sys.argv[2], sys.argv[3]

# --- read coords from generator PDB (atom-name -> xyz); LOWE/UPPE are stubs ---
co = {}
for ln in open(pdb_in):
    if ln.startswith(("ATOM", "HETATM")):
        nm = ln[12:16].strip()
        xyz = [float(ln[30:38]), float(ln[38:46]), float(ln[46:54])]
        if nm == "LOWE": nm = "LOWER"
        if nm == "UPPE": nm = "UPPER"
        if nm not in co:           # keep first occurrence (residue atoms)
            co[nm] = xyz
# shadow atom XD must coincide with backbone N
co["XD"] = list(co["N"])

# --- ICOOR tree: (atom, parent, grandparent, greatgrandparent) ---
# order guarantees every reference is defined before use
tree = [
    ("N",     "N",  "CA", "C"),    # root (self)
    ("CA",    "N",  "CA", "C"),    # phi 0 / theta 180 convention
    ("C",     "CA", "N",  "C"),
    ("UPPER", "C",  "CA", "N"),
    ("O",     "C",  "CA", "UPPER"),
    ("CB",    "CA", "N",  "C"),
    ("CG",    "CB", "CA", "N"),
    ("XD",    "CG", "CB", "N"),     # shadow of N closing the ring
    ("1HG",   "CG", "CB", "XD"),
    ("2HG",   "CG", "CB", "1HG"),
    ("1HB",   "CB", "CA", "CG"),
    ("2HB",   "CB", "CA", "1HB"),
    ("HA",    "CA", "N",  "CB"),
    ("LOWER", "N",  "CA", "C"),
]

def icoor_line(atom, p, gp, ggp):
    if atom == "N":
        return ("ICOOR_INTERNAL    %-4s %11.6f %11.6f %11.6f   %-4s %-4s %-4s"
                % ("N", 0.0, 0.0, 0.0, "N", "CA", "C"))
    if atom == "CA":
        d = vlen(vsub(co["CA"], co["N"]))
        return ("ICOOR_INTERNAL    %-4s %11.6f %11.6f %11.6f   %-4s %-4s %-4s"
                % ("CA", 0.0, 180.0, d, "N", "CA", "C"))
    d = vlen(vsub(co[atom], co[p]))
    theta = 180.0 - bond_angle(co[atom], co[p], co[gp])
    if atom == "XD":
        phi = 0.0                                  # shadow convention (PRO/NV)
    elif atom == "C":
        phi = 0.0                                  # reference dihedral
    else:
        phi = dihedral(co[atom], co[p], co[gp], co[ggp])
    return ("ICOOR_INTERNAL    %-4s %11.6f %11.6f %11.6f   %-4s %-4s %-4s"
            % (atom, phi, theta, d, p, gp, ggp))

new_icoor = [icoor_line(*t) for t in tree]

# --- splice into params: drop all existing ICOOR_INTERNAL, append fresh block ---
out = []
for ln in open(params_in):
    if ln.startswith("ICOOR_INTERNAL"):
        continue
    out.append(ln.rstrip("\n"))
# ensure block sits after MAINCHAIN_ATOMS (Rosetta is order-tolerant, append at end)
out += new_icoor
open(params_out, "w").write("\n".join(out) + "\n")
print("wrote", params_out, "with", len(new_icoor), "ICOOR lines")
