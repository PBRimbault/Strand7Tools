"""
Independent check of the app's member force diagrams against test_model3.st7 (portal frame).
Needs Strand7 (read-only) and the solved results test_model3.LSA.
Run:  python tools\\check_portal.py [model file]      (default: Test\\test_model3.st7)

The model: pitched portal, 10 m span, eaves 4.5 m, apex 6.5 m, fixed bases, UB 305x165x40.
Eaves and apex restrained out of plane (Y). Case 1: gravity self weight + 1.3 kN/m on the
rafters + 1.0 kN/m in Y on the left column; case 2: 0.8 kN/m on the rafters; case 3: Strand7
combination 1.2G + 1.5Q. The left rafter is three beams, one of them reversed; the right
rafter is a single beam running apex to eaves.

1. Reads the model's geometry, section/material data, restraints and beam loads (API).
2. Solves the frame with its own 3D Timoshenko frame stiffness method (no Strand7 solver).
3. Computes N, V, T, M at every station the app reports, by statics from the beam end forces,
   in a member frame defined here from geometry alone.
4. Compares with the app's assembled diagrams (member_forces.assemble), for each primary case
   and for Strand7's combination case.

N is compared with its sign (tension positive). For V, T and M the app's sign depends on the
member's reference axes, so the check is: same magnitude, and ONE sign relation per member and
quantity (the same at every station, including on reversed beams).
"""
import ctypes
import math
import os
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from strand7tools.strand7 import forces as mf  # noqa: E402
from strand7tools.strand7.api import S, check, get_str, model_session  # noqa: E402
from strand7tools.strand7.members import beam_axes, load_members  # noqa: E402

MODEL = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "Test", "test_model3.st7")

MEMBERS = """member,beams
COL_L,1
COL_R,2
RAF_L,6;5;3
RAF_R,4
"""

# ---------------------------------------------------------------- small linear algebra
def solve(A, b):
    n = len(b)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for c in range(n):
        p = max(range(c, n), key=lambda r: abs(M[r][c]))
        M[c], M[p] = M[p], M[c]
        if abs(M[c][c]) < 1e-14:
            raise ValueError("singular stiffness matrix (mechanism?)")
        for r in range(c + 1, n):
            f = M[r][c] / M[c][c]
            if f:
                for k in range(c, n + 1):
                    M[r][k] -= f * M[c][k]
    x = [0.0] * n
    for r in range(n - 1, -1, -1):
        x[r] = (M[r][n] - sum(M[r][k] * x[k] for k in range(r + 1, n))) / M[r][r]
    return x


def matvec(A, v):
    return [sum(a * b for a, b in zip(row, v)) for row in A]


def sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def scl(a, s):
    return tuple(x * s for x in a)


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


# ---------------------------------------------------------------- element
def k_local(E, G, A, Iy, Iz, J, Asy, Asz, L):
    """12x12, DOFs per node: ux uy uz rx ry rz in local x (axial), y, z.
    Iy: about local y (bending in x-z, shear along z, area Asz); Iz: about local z (shear along y, Asy)."""
    k = [[0.0] * 12 for _ in range(12)]

    def put(idx, sub_k):
        for i, a in enumerate(idx):
            for j, b in enumerate(idx):
                k[a][b] += sub_k[i][j]

    ea = E * A / L
    put([0, 6], [[ea, -ea], [-ea, ea]])
    gj = G * J / L
    put([3, 9], [[gj, -gj], [-gj, gj]])
    # x-y plane: v (1, 7), theta_z (5, 11)
    ph = 12 * E * Iz / (G * Asy * L ** 2) if Asy > 0 else 0.0
    a, b = 12 * E * Iz / (L ** 3 * (1 + ph)), 6 * E * Iz / (L ** 2 * (1 + ph))
    c, d = (4 + ph) * E * Iz / (L * (1 + ph)), (2 - ph) * E * Iz / (L * (1 + ph))
    put([1, 5, 7, 11], [[a, b, -a, b], [b, c, -b, d], [-a, -b, a, -b], [b, d, -b, c]])
    # x-z plane: w (2, 8), theta_y (4, 10); theta_y = -dw/dx
    ph = 12 * E * Iy / (G * Asz * L ** 2) if Asz > 0 else 0.0
    a, b = 12 * E * Iy / (L ** 3 * (1 + ph)), 6 * E * Iy / (L ** 2 * (1 + ph))
    c, d = (4 + ph) * E * Iy / (L * (1 + ph)), (2 - ph) * E * Iy / (L * (1 + ph))
    put([2, 4, 8, 10], [[a, -b, -a, -b], [-b, c, b, d], [-a, b, a, b], [-b, d, b, c]])
    return k


def f_equiv(q, L):
    """Equivalent nodal loads (local) for a uniform local load q = (qx, qy, qz) per length."""
    qx, qy, qz = q
    return [qx * L / 2, qy * L / 2, qz * L / 2, 0, -qz * L ** 2 / 12, qy * L ** 2 / 12,
            qx * L / 2, qy * L / 2, qz * L / 2, 0, qz * L ** 2 / 12, -qy * L ** 2 / 12]


# ---------------------------------------------------------------- read the model
tmp = os.path.join(tempfile.gettempdir(), "members_verify.csv")
with open(tmp, "w") as f:
    f.write(MEMBERS)

with model_session(MODEL) as uid:
    scale = mf.force_scale(uid)
    if scale is None:
        raise SystemExit("the model must be in m with N, kN or MN")

    def tot(ent):
        v = ctypes.c_long()
        check(S.St7GetTotal(uid, ent, ctypes.byref(v)))
        return v.value

    nn, nb = tot(S.tyNODE), tot(S.tyBEAM)
    xyz, fixed = {}, {}
    for n in range(1, nn + 1):
        a = (ctypes.c_double * 3)()
        check(S.St7GetNodeXYZ(uid, n, a))
        xyz[n] = (a[0], a[1], a[2])
        st, val, cs = (ctypes.c_long * 6)(), (ctypes.c_double * 6)(), ctypes.c_long()
        err = S.St7GetNodeRestraint6(uid, n, 1, ctypes.byref(cs), st, val)
        fixed[n] = list(st) if err == S.ERR7_NoError else [0] * 6     # "data not found" = free

    beams = {}
    for b in range(1, nb + 1):
        conn = (ctypes.c_long * (S.kMaxElementNode + 1))()
        check(S.St7GetElementConnection(uid, S.tyBEAM, b, conn))
        p = ctypes.c_long()
        check(S.St7GetElementProperty(uid, S.ptBEAMPROP, b, ctypes.byref(p)))
        ints = (ctypes.c_long * 5)()
        sec = (ctypes.c_double * S.kNumBeamSectionData)()
        mat = (ctypes.c_double * S.kNumMaterialData)()
        check(S.St7GetBeamPropertyData(uid, p.value, ints, sec, mat))
        m9 = (ctypes.c_double * 9)()
        check(S.St7GetBeamMaterialData(uid, p.value, m9))
        a1, a2, a3 = beam_axes(uid, b)
        beams[b] = dict(n1=conn[1], n2=conn[2], R=[a3, a1, a2],   # local x, y, z
                        E=m9[S.ipBeamModulus], G=m9[S.ipBeamShear], rho=m9[S.ipBeamDensity],
                        A=sec[S.ipAREA], I11=sec[S.ipI11], I22=sec[S.ipI22], J=sec[S.ipJ],
                        SA1=sec[S.ipSA1], SA2=sec[S.ipSA2])

    nlc = ctypes.c_long()
    check(S.St7GetNumLoadCase(uid, ctypes.byref(nlc)))
    loads = {}       # (case, beam) -> global q (per length)
    for lc in range(1, nlc.value + 1):
        g = (ctypes.c_double * 3)()
        gdir = ctypes.c_long()
        if (S.St7GetLoadCaseGravity(uid, lc, g) != S.ERR7_NoError
                or S.St7GetLoadCaseGravityDir(uid, lc, ctypes.byref(gdir)) != S.ERR7_NoError):
            g[0], gdir.value = 0.0, 0          # not a gravity load case
        ltype = ctypes.c_long()
        check(S.St7GetLoadCaseType(uid, lc, ctypes.byref(ltype)))
        print(f"load case {lc} {get_str(S.St7GetLoadCaseName, uid, lc)!r}: type {ltype.value}, "
              f"gravity {g[0]:g} dir {gdir.value}")
        for b, d in beams.items():
            q = [0.0, 0.0, 0.0]
            for direction in (1, 2, 3):
                for ID in range(1, 6):
                    l1, l2, dd = ctypes.c_long(), ctypes.c_long(), (ctypes.c_double * 6)()
                    if S.St7GetBeamDistributedForceGlobal6ID(uid, b, direction, lc, ID, ctypes.byref(l1),
                                                             ctypes.byref(l2), dd) == S.ERR7_NoError:
                        if abs(dd[0] - dd[1]) > 1e-9:
                            raise SystemExit(f"beam {b}: non-uniform load not handled by this check")
                        q[direction - 1] += dd[0]
            loads[(lc, b)] = q
        loads[("gravity", lc)] = (g[0], gdir.value)

    members, member_errors = load_members(uid, tmp)
    if member_errors:
        raise SystemExit(f"members could not be built: {member_errors}")
    app = {}
    with mf.result_file(uid, MODEL[:-4] + ".LSA", "existing") as (npr, nsec):
        names = {c: get_str(S.St7GetResultCaseName, uid, c) for c in range(1, npr + nsec + 1)}
        for c in names:
            for m in members:
                app[(c, m.definition.member_id)] = mf.assemble(uid, m, c, scale)

# ---------------------------------------------------------------- self weight
print(f"\nforce scale to kN: {scale:g}")
for b, d in beams.items():
    d["L"] = math.dist(xyz[d["n1"]], xyz[d["n2"]])


def case_loads(lc):
    out = {}
    g, gdir = loads[("gravity", lc)]
    for b, d in beams.items():
        q = list(loads[(lc, b)])
        if g and gdir:
            # Strand7 gravity: signed acceleration g along global direction gdir (1=X, 2=Y, 3=Z)
            q[abs(gdir) - 1] += g * d["rho"] * d["A"]
        out[b] = tuple(q)
    return out


# ---------------------------------------------------------------- solve
nodes = sorted(xyz)
dof = {n: 6 * i for i, n in enumerate(nodes)}
ndof = 6 * len(nodes)


def T12(R):
    t = [[0.0] * 12 for _ in range(12)]
    for blk in range(4):
        for i in range(3):
            for j in range(3):
                t[3 * blk + i][3 * blk + j] = R[i][j]
    return t


def transpose(M):
    return [list(r) for r in zip(*M)]


def mmul(A, B):
    Bt = transpose(B)
    return [[dot(r, c) for c in Bt] for r in A]


for d in beams.values():
    d["k"] = k_local(d["E"], d["G"], d["A"], d["I11"], d["I22"], d["J"], d["SA1"], d["SA2"], d["L"])
    d["T"] = T12(d["R"])
    d["kg"] = mmul(mmul(transpose(d["T"]), d["k"]), d["T"])


def analyse(lc):
    q_by_beam = case_loads(lc)
    K = [[0.0] * ndof for _ in range(ndof)]
    F = [0.0] * ndof
    for b, d in beams.items():
        idx = [dof[d["n1"]] + i for i in range(6)] + [dof[d["n2"]] + i for i in range(6)]
        for i, a in enumerate(idx):
            for j, c in enumerate(idx):
                K[a][c] += d["kg"][i][j]
        ql = tuple(dot(r, q_by_beam[b]) for r in d["R"])
        fe = matvec(transpose(d["T"]), f_equiv(ql, d["L"]))
        for i, a in enumerate(idx):
            F[a] += fe[i]
    free = [dof[n] + i for n in nodes for i in range(6) if not fixed[n][i]]
    u = [0.0] * ndof
    sol = solve([[K[a][c] for c in free] for a in free], [F[a] for a in free])
    for a, v in zip(free, sol):
        u[a] = v
    ends = {}
    for b, d in beams.items():
        idx = [dof[d["n1"]] + i for i in range(6)] + [dof[d["n2"]] + i for i in range(6)]
        ul = matvec(d["T"], [u[a] for a in idx])
        ql = tuple(dot(r, q_by_beam[b]) for r in d["R"])
        fl = [x - y for x, y in zip(matvec(d["k"], ul), f_equiv(ql, d["L"]))]
        fg = matvec(transpose(d["T"]), fl)          # global end forces (node on element)
        ends[b] = (fg[0:3], fg[3:6], fg[6:9], fg[9:12])
    return q_by_beam, ends


def internal(member, chain_beam, x, q_by_beam, ends):
    """Member-frame internal forces at member chainage x (independent of the app's axes)."""
    d = beams[chain_beam.beam]
    behind = d["n1"] if chain_beam.forward else d["n2"]
    f0, m0 = (ends[chain_beam.beam][0], ends[chain_beam.beam][1]) if chain_beam.forward else \
        (ends[chain_beam.beam][2], ends[chain_beam.beam][3])
    t = member.direction
    P0 = xyz[behind]
    s = x - chain_beam.x0                         # distance from the behind end
    P = add(xyz[member.nodes[0]], scl(t, x))
    q = q_by_beam[chain_beam.beam]
    Fc = scl(add(f0, scl(q, s)), -1)
    Pm = add(P0, scl(t, s / 2))
    Mc = scl(add(add(m0, cross(sub(P0, P), f0)), cross(sub(Pm, P), scl(q, s))), -1)
    e1 = (0.0, 1.0, 0.0)                          # out of plane (global Y)
    e2 = cross(t, e1)                             # in plane, normal to the member
    return {"N": dot(Fc, t), "V_out": dot(Fc, e1), "V_in": dot(Fc, e2), "T": dot(Mc, t),
            "M_major": dot(Mc, e1), "M_minor": dot(Mc, e2)}


# app quantity -> my quantity (major axis bending = about axis 1 = global Y here)
PAIRS = [("N", "N"), ("Vy", "V_out"), ("Vz", "V_in"), ("T", "T"), ("My", "M_major"), ("Mz", "M_minor")]

results = {lc: analyse(lc) for lc in range(1, nlc.value + 1)}

print("\nReactions check (sum of applied vertical load per case):")
for lc in results:
    qb = results[lc][0]
    print(f"  case {lc}: total load (X, Y, Z) = "
          + ", ".join(f"{sum(qb[b][i] * beams[b]['L'] for b in beams):8.3f}" for i in range(3)))

worst = 0.0
problems = 0
for c, cname in names.items():
    if c <= nlc.value:
        q_by_beam, ends = results[c]
        combo = None
    else:
        combo = {1: 1.2, 2: 1.5}            # Strand7 combination case as named: 1.2G + 1.5Q
    print(f"\n=== result case {c}: {cname}" + ("   (compared with 1.2 x case 1 + 1.5 x case 2)" if combo else ""))
    print(f"  {'member':<6} {'qty':<3} {'peak app':>9} {'peak here':>10} {'max |diff|':>11} "
          f"{'sign rel.':>9}")
    for m in members:
        mid = m.definition.member_id
        rows = []
        for cb, pts in app[(c, mid)]:
            for x, v in pts:
                if combo:
                    mine = {}
                    for lc, fac in combo.items():
                        part = internal(m, cb, x, *results[lc])
                        for k_, val in part.items():
                            mine[k_] = mine.get(k_, 0.0) + fac * val
                else:
                    mine = internal(m, cb, x, q_by_beam, ends)
                rows.append((v, mine))
        for qa, qm in PAIRS:
            peak_a = max(abs(v[qa]) for v, _ in rows)
            peak_m = max(abs(mm[qm]) for _, mm in rows)
            if qa == "N":
                diff = max(abs(v[qa] - mm[qm]) for v, mm in rows)
                rel = "+1"
            else:
                signs = {round(math.copysign(1, v[qa] * mm[qm])) for v, mm in rows
                         if abs(v[qa]) > 0.02 * max(peak_a, 1e-9) and abs(mm[qm]) > 0.02 * max(peak_m, 1e-9)}
                rel = "n/a" if not signs else ("MIXED" if len(signs) > 1 else f"{signs.pop():+d}")
                diff = max(abs(abs(v[qa]) - abs(mm[qm])) for v, mm in rows)
            tol = max(0.01, 0.01 * peak_m)
            flag = "" if diff <= tol and rel != "MIXED" else "   <-- CHECK"
            problems += bool(flag)
            worst = max(worst, diff / max(peak_m, 1.0))
            print(f"  {mid:<6} {qa:<3} {peak_a:>9.3f} {peak_m:>10.3f} {diff:>11.4f} {rel:>9}{flag}")
print(f"\nworst difference relative to the peak: {worst:.2%}")
for c in names:
    for m in members:
        for w in mf.continuity_warnings(m, app[(c, m.definition.member_id)]):
            print(f"CONTINUITY (case {c}, {m.definition.member_id}): {w}")
            problems += 1
print("All checks passed." if not problems else f"{problems} PROBLEM(S) - see CHECK lines above.")

if os.environ.get("DUMP"):          # e.g. set DUMP=RAF_L: station-by-station listing, case 1
    q_by_beam, ends = results[1]
    for m in members:
        if m.definition.member_id != os.environ["DUMP"]:
            continue
        for cb, pts in app[(1, m.definition.member_id)]:
            for x, v in pts:
                mine = internal(m, cb, x, q_by_beam, ends)
                print(f"b{cb.beam} fwd={cb.forward} x={x:6.3f}  T app {v['T']*1e3:9.4f}  here {mine['T']*1e3:9.4f}   "
                      f"Mz app {v['Mz']*1e3:9.4f} here {mine['M_minor']*1e3:9.4f}  Vy app {v['Vy']*1e3:8.4f} here {mine['V_out']*1e3:8.4f}  (x1e-3)")

sys.exit(1 if problems else 0)
