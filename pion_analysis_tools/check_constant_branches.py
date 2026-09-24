#!/usr/bin/env python3
"""
check_constant_branches.py  --  C4 production sanity check.

For every branch in a *_full_history.npz file, report whether it takes MORE
THAN ONE distinct value across the whole file.  A branch that is declared,
correctly named, and silently constant (all zeros, all -1, one repeated code)
is the failure class that an all-zero scan misses -- e.g. ph_muAncestorID
(constant -1) or ap_par (one constant value).

Branches that are *legitimately* constant are listed in ALLOW_CONSTANT and do
not cause a failure.  Anything else that is constant (or entirely empty) is
reported and the script exits non-zero, so it can be wired into a production
pipeline and fail loudly.

Usage:
    python3 check_constant_branches.py FILE.npz [FILE2.npz ...]
    python3 check_constant_branches.py /path/to/dir            # scans *full_history*.npz

Exit code 0 = all good, 1 = at least one unexpected constant/empty branch.
"""

import sys
import glob
import os
import numpy as np

# --------------------------------------------------------------------------
# Branches that are constant *by design*.  Keep this list short and justified;
# every entry is a promise that a single value is physically correct here.
# --------------------------------------------------------------------------
ALLOW_CONSTANT = {
    "pp_creator",       # every pion_photons row is a Cherenkov photon
    "pp_end_process",   # ... that reached a PMT (Transportation)
    "pp_made_pe",       # pion_photons rows are only written when a PE was made
    "root_file",        # provenance string, one per file
    "trigger_type",     # constant IFF the run is genuinely untriggered -- see note
    # Beam-configuration branches: constant for a single-species, mono-energetic
    # run (the normal case here).  REMOVE these from the allow-list for a
    # mixed-species or multi-energy production, where they must vary.
    "pid",              # primary PDG  (e.g. 211 for a pure pi+ beam)
    "energy",           # primary energy (fixed for a mono-energetic beam)
}

# Substrings marking branches whose constancy is a KNOWN open issue (reported
# as WARN, not FAIL, so the check stays green while the sim side is fixed).
KNOWN_ISSUE_SUBSTR = (
    "digi_hit_trigger",   # C3: DAQ/trigger not configured in current runs
    "trigger_time",       # C3
)


def flatten(entry):
    """Return a flat 1-D array of all values in one npz branch entry.

    Branches are stored as an object array (one element per event), each
    element itself a variable-length array; or as a plain numeric/str array.
    """
    if isinstance(entry, np.ndarray) and entry.dtype == object:
        parts = []
        for e in entry:
            a = np.asarray(e).ravel()
            if a.size:
                parts.append(a)
        if not parts:
            return np.array([])
        try:
            return np.concatenate(parts)
        except ValueError:
            # mixed dtypes (e.g. strings) -> concatenate as python objects
            out = []
            for p in parts:
                out.extend(list(p))
            return np.array(out, dtype=object)
    return np.asarray(entry).ravel()


def n_distinct(flat):
    """Number of distinct values, robust to strings and NaNs."""
    if flat.size == 0:
        return 0
    if flat.dtype == object or flat.dtype.kind in ("U", "S"):
        return len(set(map(str, flat.tolist())))
    # numeric: treat all NaNs as one value
    finite = flat[np.isfinite(flat)] if flat.dtype.kind == "f" else flat
    nd = np.unique(finite).size
    if flat.dtype.kind == "f" and finite.size != flat.size:
        nd += 1  # NaN counts as one extra distinct value
    return nd


def check_file(path):
    d = np.load(path, allow_pickle=True)
    keys = list(d.keys())

    empty, constant, known = [], [], []
    for k in keys:
        flat = flatten(d[k])
        nd = n_distinct(flat)
        if nd == 0:
            empty.append(k)
        elif nd == 1:
            if k in ALLOW_CONSTANT:
                continue
            if any(s in k for s in KNOWN_ISSUE_SUBSTR):
                known.append((k, flat))
            else:
                constant.append((k, flat))

    print(f"\n=== {os.path.basename(path)}  ({len(keys)} branches) ===")

    def one_val(flat):
        try:
            return repr(flat[0]) if flat.size else "n/a"
        except Exception:
            return "?"

    if known:
        print("  KNOWN-ISSUE constant (WARN, not failing):")
        for k, flat in known:
            print(f"    {k:24s} = {one_val(flat)}")
    if empty:
        print("  EMPTY branches (never filled):")
        for k in empty:
            print(f"    {k}")
    if constant:
        print("  UNEXPECTED constant branches (FAIL):")
        for k, flat in constant:
            print(f"    {k:24s} = {one_val(flat)}")
    if not (empty or constant or known):
        print("  OK -- every branch varies (or is an allowed constant).")

    return len(empty) + len(constant)  # known-issues do not fail the check


def collect(args):
    files = []
    for a in args:
        if os.path.isdir(a):
            files += sorted(glob.glob(os.path.join(a, "*full_history*.npz")))
        else:
            files.append(a)
    return files


def main(argv):
    if len(argv) < 2:
        print(__doc__)
        return 2
    files = collect(argv[1:])
    if not files:
        print("No .npz files found.")
        return 2
    fails = sum(check_file(f) for f in files)
    print(f"\n{'FAIL' if fails else 'PASS'}: "
          f"{fails} unexpected empty/constant branch(es) across {len(files)} file(s).")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
