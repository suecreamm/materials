"""
link_ph_outputs.py (safe + flexible for EPW)

Creates symlinks in EPW's dvscf_dir so EPW can open:
  - PREFIX.dyn_qN
  - PREFIX.dvscf_qN  (some EPW builds)
  - PREFIX.dvscfN_1  (other naming conventions)

Safety rules:
  - Never overwrite an existing regular file.
  - Never create self-referential symlinks (src == dst).
  - Exclude already-normalized destination filenames from being re-used as sources.
"""

from __future__ import annotations

from pathlib import Path
import argparse
import os
import re
import sys


def safe_symlink(src: Path, dst: Path) -> None:
    """Create/replace a symlink dst -> src safely."""
    src = src.resolve()
    dst_parent = dst.parent
    dst_parent.mkdir(parents=True, exist_ok=True)

    # If dst exists as a regular file, do NOT touch it.
    if dst.exists() and not dst.is_symlink():
        print(f"[skip] dst is a regular file (keep): {dst}")
        return

    # If dst is a symlink pointing to the same target, skip.
    if dst.is_symlink():
        try:
            if dst.resolve() == src:
                print(f"[skip] exists (same target): {dst}")
                return
        except FileNotFoundError:
            # broken link; we'll recreate
            pass

    # Prevent self-link (dst path equals src path)
    try:
        if dst.resolve() == src:
            print(f"[skip] would create self-link: {dst}")
            return
    except Exception:
        # dst may not exist yet; compare absolute paths
        if dst_parent.resolve() / dst.name == src:
            print(f"[skip] would create self-link: {dst}")
            return

    if dst.is_symlink() or dst.exists():
        dst.unlink()

    os.symlink(src, dst)
    print(f"[link] {dst.name} -> {src}")


def main(prefix: str, dvscf_dir: Path) -> int:
    cwd = Path.cwd()
    dvscf_dir = dvscf_dir.resolve()

    print(f"[info] cwd       : {cwd}")
    print(f"[info] dvscf_dir : {dvscf_dir}")

    # -------------------------
    # dynN -> dyn_qN
    # -------------------------
    dyn_files = sorted(cwd.glob(f"{prefix}.dyn*"))
    if not dyn_files:
        print("[error] no dyn files found in current directory")
        return 1

    dyn_pat = re.compile(rf"^{re.escape(prefix)}\.dyn(\d+)$")
    print(f"[info] found {len(dyn_files)} dyn files")

    for f in dyn_files:
        m = dyn_pat.match(f.name)
        if not m:
            continue
        iq = m.group(1)
        safe_symlink(f, dvscf_dir / f"{prefix}.dyn_q{iq}")

    # -------------------------
    # dvscf sources (flexible)
    # -------------------------
    p = re.escape(prefix)

    # Accept:
    #   PREFIX.dvscfN_M
    #   PREFIX.PREFIX.dvscfN_M
    #   PREFIX.dvscfN
    #   PREFIX.PREFIX.dvscfN
    dvscf_patterns = [
        re.compile(rf"^{p}\.dvscf(\d+)_(\d+)$"),
        re.compile(rf"^{p}\.{p}\.dvscf(\d+)_(\d+)$"),
        re.compile(rf"^{p}\.dvscf(\d+)$"),
        re.compile(rf"^{p}\.{p}\.dvscf(\d+)$"),
    ]

    # Exclude already-normalized destination names from source candidates
    exclude_src = re.compile(rf"^{p}\.dvscf(\d+)_1$|^{p}\.dvscf_q(\d+)$")

    def parse(name: str):
        if exclude_src.match(name):
            return None
        for pat in dvscf_patterns:
            m = pat.match(name)
            if m:
                iq = int(m.group(1))
                ipert = int(m.group(2)) if m.lastindex and m.lastindex >= 2 else 1
                return iq, ipert
        return None

    roots = [
        dvscf_dir,              # ok to search, but we exclude dst-names above
        cwd / "tmp" / "_ph0",
        cwd / "out" / "_ph0",
    ]

    hits: dict[int, tuple[int, Path]] = {}
    for root in roots:
        if not root.exists():
            continue
        for f in root.rglob("*dvscf*"):
            if not f.is_file():
                continue
            parsed = parse(f.name)
            if not parsed:
                continue
            iq, ipert = parsed
            if iq not in hits:
                hits[iq] = (ipert, f)
            else:
                best_ipert, _ = hits[iq]
                if best_ipert != 1 and ipert == 1:
                    hits[iq] = (ipert, f)
                elif best_ipert != 1 and ipert < best_ipert:
                    hits[iq] = (ipert, f)

    if not hits:
        print("[error] no dvscf source files found")
        return 2

    print(f"[info] found dvscf q-points: {sorted(hits.keys())}")

    # Create BOTH naming conventions EPW might request
    for iq in sorted(hits.keys()):
        _ipert, src = hits[iq]
        safe_symlink(src, dvscf_dir / f"{prefix}.dvscf{iq}_1")
        safe_symlink(src, dvscf_dir / f"{prefix}.dvscf_q{iq}")

    print("[done] links created safely")
    return 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("prefix")
    ap.add_argument("--dvscf_dir", default="./tmp/_ph0")
    args = ap.parse_args()
    sys.exit(main(args.prefix, Path(args.dvscf_dir)))

