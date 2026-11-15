"""
Build the design members from the member table and print them, so the structure can be
confirmed before any checking. Read-only.

    python tools\\list_members.py [project file]      (default: Test\\test_model2.toml)
"""

import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from strand7tools.project import load_project  # noqa: E402
from strand7tools.strand7.api import model_session, total_beams  # noqa: E402
from strand7tools.strand7.members import load_members  # noqa: E402
from strand7tools.member_table import read_rows  # noqa: E402

DEFAULT_PROJECT = os.path.join(ROOT, "Test", "test_model2.toml")


def print_member(m):
    d = m.definition
    eff = m.effective_lengths()
    print("=" * 72)
    print(f"Member {d.member_id}" + (f"   ({d.comment})" if d.comment else ""))
    print(f"  Section : {m.section_name}")
    print(f"  Length  : {m.length:.4f}   direction ({m.direction[0]:.3f}, "
          f"{m.direction[1]:.3f}, {m.direction[2]:.3f})")
    a1, a2 = m.ref_axis1, m.ref_axis2
    print(f"  Axes    : axis 1 = ({a1[0]:.3f}, {a1[1]:.3f}, {a1[2]:.3f})  "
          f"axis 2 = ({a2[0]:.3f}, {a2[1]:.3f}, {a2[2]:.3f})   [beam {m.chain[0].beam}]")
    print("            major axis (y-y) = bending about axis 1 -> Moment 2, Shear 2, SA2 "
          "(axis 2 = web direction)")
    print("            minor axis (z-z) = bending about axis 2 -> Moment 1, Shear 1, SA1")
    print("  Beams   : " + " ".join(f"{c.beam}{'' if c.forward else '*'}" for c in m.chain)
          + "     (* = beam runs against the member direction)")
    print("  Sign factors vs member axes (+1 = same direction, -1 = opposite):")
    for c in m.chain:
        print(f"            beam {c.beam:>4}: axis3 {c.f3:+d}   axis1 {c.f1:+d}   axis2 {c.f2:+d}")
    print("  Nodes   : " + " ".join(str(n) for n in m.nodes))
    last = len(m.nodes) - 1
    print("            " + "--".join("X" if i in (0, last) else "o" for i in range(len(m.nodes))))
    print()
    print("  Effective lengths:")
    rows = (
        ("Flexural buckling about y", eff["Lcr_y"], d.kcr_y, d.lcr_y, "Lcr,y"),
        ("Flexural buckling about z", eff["Lcr_z"], d.kcr_z, d.lcr_z, "Lcr,z"),
        ("Lateral-torsional buckling", eff["Lw"], d.kz, d.lw, "Lw"),
    )
    for label, value, k, override, symbol in rows:
        if value is None:
            print(f"    {label:<27}: not checked")
        else:
            src = "typed value" if override is not None else f"k={k:g} x {m.length:.4f}"
            print(f"    {label:<27}: {symbol:<6} = {value:.4f}   ({src})")
    if m.warnings:
        print()
        for w in m.warnings:
            print(f"  WARNING: {w}")
    print()


def main(project_file=DEFAULT_PROJECT):
    project = load_project(project_file)
    rows = read_rows(project.members)
    if any(r.get("restraints") for r in rows):
        print("NOTE: the 'restraints' column is no longer used - members now end at "
              "restraints. Split such members in the table.\n")
    with model_session(project.model) as uid:
        n_beams = total_beams(uid)
        members, errors = load_members(uid, project.members)
    print(f"Model: {project.model}")
    print(f"Beams in model: {n_beams}")
    print(f"Member table: {project.members}\n")
    for mid, msg in errors:
        print(f"Member {mid}: ERROR - {msg}\n")
    for m in members:
        print_member(m)
    print(f"{len(members)} of {len(rows)} members built.")
    owned = {c.beam for m in members for c in m.chain}
    unassigned = [b for b in range(1, n_beams + 1) if b not in owned]
    if unassigned:
        shown = ", ".join(str(b) for b in unassigned[:30])
        more = " ..." if len(unassigned) > 30 else ""
        print(f"Beams not in any member ({len(unassigned)}): {shown}{more}")


if __name__ == "__main__":
    main(*sys.argv[1:2])
