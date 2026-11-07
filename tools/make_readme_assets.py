"""Generate the README figures in assets/.

All quantities are computed here from the discretization the solver uses:
cubic (degree 3) B-splines, periodic in x and z, clamped (open knot vector) in
the wall-normal direction y, on the 80 x 120 x 80 channel configured in
src/lattice.f90 (77 x 120 x 77 coefficients).

Usage (from the repository root):
    python tools/make_readme_assets.py
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch, Polygon
from scipy.interpolate import BSpline

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "assets"
OUT.mkdir(exist_ok=True)
STYLE = ROOT.parent / "_style" / "lab.mplstyle"
if STYLE.exists():
    plt.style.use(str(STYLE))

INK, TEAL, AMBER, SLATE, CORAL, GRAY = (
    "#1B2A41", "#0F9D8A", "#F2A541", "#3D5A80", "#E07A5F", "#8D99AE")

DEG = 3
NX, NY, NZ = 80, 120, 80                     # src/lattice.f90
NXC, NZC = NX - DEG, NZ - DEG                # periodic coefficients
LX, LY, LZ = 4 * np.pi, 2.0, 4 * np.pi / 3


# --------------------------------------------------------------------------
# 1D Galerkin mass matrices (exact Gauss quadrature, deg + 2 points per span)
# --------------------------------------------------------------------------
def gram(knots, nbasis, a, b):
    """Mass matrix M_ij = int_a^b B_i B_j dx for the given knot vector."""
    xg, wg = np.polynomial.legendre.leggauss(DEG + 2)
    brk = np.unique(knots[(knots >= a) & (knots <= b)])
    M = np.zeros((nbasis, nbasis))
    basis = [BSpline.basis_element(knots[i:i + DEG + 2], extrapolate=False)
             for i in range(nbasis)]
    for lo, hi in zip(brk[:-1], brk[1:]):
        x = 0.5 * (hi - lo) * xg + 0.5 * (hi + lo)
        w = 0.5 * (hi - lo) * wg
        V = np.nan_to_num(np.array([bi(x) for bi in basis]))
        M += (V * w) @ V.T
    return M


def periodic_mass(n, L):
    """Circulant mass matrix of n periodic cubic B-splines on [0, L)."""
    h = L / n
    knots = h * np.arange(-DEG, n + DEG + 1)
    M = gram(knots, n + DEG, 0.0, L)
    # fold the DEG wrap-around functions onto their periodic images
    P = np.zeros((n + DEG, n))
    for i in range(n + DEG):
        P[i, (i - DEG) % n] = 1.0
    return P.T @ M @ P


def clamped_mass(nbasis, L):
    """Mass matrix of clamped cubic B-splines on [0, L] (wall-normal direction)."""
    nint = nbasis - DEG
    inner = np.linspace(0.0, L, nint + 1)
    knots = np.r_[[0.0] * DEG, inner, [L] * DEG]
    return gram(knots, nbasis, 0.0, L)


def periodic_symbol(theta):
    """Eigenvalues of the uniform periodic cubic Galerkin mass matrix (h = 1)."""
    c = np.array([2416.0, 1191.0, 120.0, 1.0]) / 5040.0
    return c[0] + 2 * sum(c[k] * np.cos(k * theta) for k in (1, 2, 3))


# --------------------------------------------------------------------------
# Hero banner
# --------------------------------------------------------------------------
def iso(x, y, z):
    """Oblique projection of channel coordinates onto the page."""
    return np.array([x + 0.55 * z, y + 0.32 * z])


def face(ax, pts, **kw):
    ax.add_patch(Polygon([iso(*p) for p in pts], closed=True, **kw))


def hero():
    fig = plt.figure(figsize=(8, 2.8), dpi=200)
    fig.patch.set_facecolor("white")

    # ---- text block
    fig.text(0.035, 0.80, "petsc-bspline-compressible", fontsize=17,
             fontweight="bold", color=INK, ha="left", va="center")
    fig.text(0.035, 0.665,
             "Parallel Fortran + PETSc solver for the 3D compressible\n"
             "Navier-Stokes equations with Galerkin cubic B-splines",
             fontsize=9.6, color=INK, ha="left", va="top", linespacing=1.35)
    items = [
        (TEAL, "DMDA domain decomposition, MPI line communicators"),
        (AMBER, "KSP GMRES / CG wall-normal solves, SuperLU_DIST option"),
        (SLATE, "FFTW diagonalization of periodic mass matrices"),
        (CORAL, "Sum-factorized 3D operators, MMS-verified kernels"),
    ]
    for k, (c, t) in enumerate(items):
        y = 0.40 - k * 0.088
        fig.patches.append(FancyBboxPatch((0.037, y - 0.018), 0.011, 0.036,
                                          boxstyle="round,pad=0.002", fc=c,
                                          ec="none", transform=fig.transFigure))
        fig.text(0.058, y, t, fontsize=8.1, color=INK, va="center")

    # ---- partitioned channel
    ax = fig.add_axes([0.55, 0.25, 0.385, 0.72])
    ax.set_aspect("equal")
    ax.axis("off")
    Lx, Ly, Lz = 4.0, 1.25, 1.9            # drawing proportions
    px, pz = 4, 2                          # ranks in x and z, one rank in y
    tints = ["#E3F4F1", "#F3F6FA"]
    for i in range(px):
        for k in range(pz):
            x0, x1 = Lx * i / px, Lx * (i + 1) / px
            z0, z1 = Lz * k / pz, Lz * (k + 1) / pz
            c = tints[(i + k) % 2]
            # top face of each subdomain
            face(ax, [(x0, Ly, z0), (x1, Ly, z0), (x1, Ly, z1), (x0, Ly, z1)],
                 fc=c, ec=INK, lw=0.6)
    for i in range(px):                    # front faces (z = 0)
        x0, x1 = Lx * i / px, Lx * (i + 1) / px
        face(ax, [(x0, 0, 0), (x1, 0, 0), (x1, Ly, 0), (x0, Ly, 0)],
             fc=tints[i % 2], ec=INK, lw=0.6)
    for k in range(pz):                    # right faces (x = Lx)
        z0, z1 = Lz * k / pz, Lz * (k + 1) / pz
        face(ax, [(Lx, 0, z0), (Lx, 0, z1), (Lx, Ly, z1), (Lx, Ly, z0)],
             fc=tints[(px - 1 + k) % 2], ec=INK, lw=0.6)

    # line communicators: an x-line, a y-line and a z-line through one rank
    xl = np.array([iso(x, 0.62 * Ly, 0.0) for x in (0, Lx)])
    ax.plot(*xl.T, color=TEAL, lw=2.6, solid_capstyle="round")
    yl = np.array([iso(0.62 * Lx, y, 0.0) for y in (0, Ly)])
    ax.plot(*yl.T, color=AMBER, lw=2.6, solid_capstyle="round")
    zl = np.array([iso(Lx, 0.62 * Ly, z) for z in (0, Lz)])
    ax.plot(*zl.T, color=SLATE, lw=2.6, solid_capstyle="round")
    ax.text(*(xl[0] + [0.05, 0.08]), "x: FFTW on commx", fontsize=7, color=TEAL,
            fontweight="bold")
    ax.text(*(yl[1] + [0.10, -0.26]), "y: KSP on commy", fontsize=7,
            color="#B57411", fontweight="bold")
    ax.text(*(iso(Lx, 0, Lz) + [0.04, -0.05]), "z: FFTW\non commz", fontsize=7,
            color=SLATE, fontweight="bold")
    ax.text(*iso(0.03, Ly + 0.07, Lz), "8 ranks (4 x 1 x 2 DMDA)", fontsize=7,
            color=INK, va="bottom")
    ax.text(*(iso(0, 0, 0) + [0.0, -0.22]), "no-slip walls at y = 0 and y = 2, periodic in x and z",
            fontsize=6.5, color=GRAY, ha="left")
    ax.set_xlim(-0.15, Lx + 0.55 * Lz + 0.75)
    ax.set_ylim(-0.3, Ly + 0.32 * Lz + 0.3)

    # ---- B-spline basis inset
    ins = fig.add_axes([0.765, 0.05, 0.17, 0.24])
    knots = np.arange(-3, 11)
    xs = np.linspace(0, 7, 600)
    for i in range(len(knots) - 4):
        b = BSpline.basis_element(knots[i:i + 5], extrapolate=False)
        v = np.nan_to_num(b(xs))
        hl = i == 5
        ins.plot(xs, v, color=TEAL if hl else GRAY, lw=1.8 if hl else 0.8)
    ins.set_xlim(0, 7)
    ins.set_ylim(0, 0.8)
    ins.set_xticks([])
    ins.set_yticks([])
    ins.grid(False)
    for s in ins.spines.values():
        s.set_visible(False)
    ins.set_title("cubic B-spline basis", fontsize=6.5, color=INK, pad=1,
                  fontweight="normal")

    with plt.rc_context({"savefig.bbox": "standard"}):
        fig.savefig(OUT / "hero.png", dpi=200,
                    facecolor="white")
    plt.close(fig)


# --------------------------------------------------------------------------
# Card 1: solver stack
# --------------------------------------------------------------------------
def stack_card():
    fig, ax = plt.subplots(figsize=(8, 3.5))
    ax.axis("off")
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 4.0)

    def box(x, y, w, h, title, body, color):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.12",
                                    fc="white", ec=color, lw=1.8))
        ax.add_patch(FancyBboxPatch((x, y + h - 0.36), w, 0.36,
                                    boxstyle="round,pad=0.02,rounding_size=0.12",
                                    fc=color, ec=color, lw=1.8))
        ax.text(x + w / 2, y + h - 0.18, title, ha="center", va="center",
                fontsize=9.2, color="white", fontweight="bold")
        ax.text(x + 0.12, y + h - 0.48, body, ha="left", va="top",
                fontsize=7.6, color=INK, linespacing=1.32, family="DejaVu Sans")

    def arrow(x0, y0, x1, y1):
        ax.annotate("", (x1, y1), (x0, y0),
                    arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.1))

    box(0.1, 2.15, 3.1, 1.8, "Galerkin assembly",
        "Gauss quadrature, deg + 2 points\n"
        "1D mass / derivative / stiffness\n"
        "  BB, BD, BL per direction\n"
        "1D triple and quadruple product\n"
        "  tensors for nonlinear terms", TEAL)
    box(3.45, 2.15, 3.1, 1.8, "PETSc DM, Mat, Vec",
        "DMDACreate3d: 77 x 120 x 77\n"
        "  periodic / ghosted / periodic\n"
        "DMDACreate1d per line (Comm_split)\n"
        "DMCreateMatrix, MATAIJ\n"
        "DMGlobalToLocal ghost exchange", SLATE)
    box(6.8, 2.15, 3.1, 1.8, "Operator application",
        "3D operators never assembled:\n"
        "  (A_x kron A_y kron A_z) alpha\n"
        "  as MatMult along DMDA lines\n"
        "Nonlinear Galerkin products from\n"
        "  1D triple-product tensors", CORAL)
    box(0.1, 0.05, 4.75, 1.75, "Mass solve, periodic x and z",
        "Circulant mass matrix diagonalized by FFT\n"
        "MatCreateFFT(MATFFTW) on commx / commz\n"
        "VecScatterPetscToFFTW, MatMult (forward)\n"
        "pointwise divide by the mass-matrix symbol\n"
        "MatMultTranspose (inverse), scale 1/N", SLATE)
    box(5.15, 0.05, 4.75, 1.75, "Mass solve, wall-normal y",
        "KSPCreate on commy, KSPSetFromOptions\n"
        "KSPGMRES: mass matrix with Dirichlet rows\n"
        "KSPCG: symmetric positive definite mass\n"
        "rtol = atol = 1e-14, default PC: ILU(0)\n"
        "runtime: preonly + LU via SuperLU_DIST", AMBER)
    arrow(3.2, 3.05, 3.45, 3.05)
    arrow(6.55, 3.05, 6.8, 3.05)
    arrow(2.4, 2.15, 2.4, 1.82)
    arrow(7.5, 2.15, 7.5, 1.82)
    fig.suptitle("Solver stack: Galerkin B-splines on PETSc DMDA, KSP and FFTW",
                 x=0.02, ha="left", fontsize=14, fontweight="bold", color=INK)
    fig.text(0.02, 0.005,
             "Every RK3 stage: Galerkin residuals for mass, momentum and energy, then "
             "one tensor-product mass solve per variable and nonlinear product.",
             fontsize=8, color=GRAY)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(OUT / "card_solver_stack.png")
    plt.close(fig)


# --------------------------------------------------------------------------
# Card 2: operator storage, factored vs assembled
# --------------------------------------------------------------------------
def storage_card():
    Mx = periodic_mass(NXC, LX)
    My = clamped_mass(NY, LY)
    Mz = periodic_mass(NZC, LZ)
    tol = 1e-15
    nnz = [int((np.abs(M) > tol).sum()) for M in (Mx, My, Mz)]
    n3 = NXC * NY * NZC
    nnz3 = nnz[0] * nnz[1] * nnz[2]
    bytes_per = 16 + 4                     # complex double value + int32 column
    fac = sum(nnz) * bytes_per
    full = nnz3 * bytes_per

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(8, 3.7),
                                 gridspec_kw=dict(width_ratios=[1.0, 1.25]))
    a1.spy(np.abs(My[:40, :40]) > tol, markersize=2.6, color=AMBER)
    a1.set_title("wall-normal mass matrix\n(first 40 of 120 rows)", fontsize=10)
    a1.set_xticks([0, 39])
    a1.set_yticks([0, 39])
    a1.grid(False)

    labels = ["1D factors\nx + y + z", "assembled\n3D matrix"]
    vals = [sum(nnz), nnz3]
    bars = a2.bar([0, 1], vals, color=[TEAL, AMBER], width=0.6)
    a2.set_yscale("log")
    a2.set_ylim(1e2, 1e10)
    a2.set_xticks([0, 1], labels)
    a2.set_ylabel("stored nonzeros")
    a2.grid(axis="x", visible=False)
    a2.text(0, vals[0] * 2.2, f"{vals[0]:,}\n{fac / 1e3:.0f} kB", ha="center",
            fontsize=10, fontweight="bold", color=INK)
    a2.text(1, vals[1] * 2.2, f"{vals[1] / 1e6:,.0f} million\n{full / 1e9:.1f} GB",
            ha="center", fontsize=10, fontweight="bold", color=INK)
    a2.set_title(f"one mass operator, {n3:,} unknowns", fontsize=10)
    fig.suptitle("3D Galerkin operators applied from 1D sparse factors",
                 x=0.02, ha="left", fontsize=14, fontweight="bold", color=INK)
    fig.text(0.02, 0.0,
             f"Cubic B-splines on the 80 x 120 x 80 channel ({NXC} x {NY} x {NZC} coefficients), "
             "complex AIJ storage (16 B value + 4 B index). Kronecker nonzeros = product of 1D nonzeros.",
             fontsize=7.4, color=GRAY)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(OUT / "card_operator_storage.png")
    plt.close(fig)
    return nnz, nnz3, fac, full


# --------------------------------------------------------------------------
# Card 3: mass-matrix spectrum and conditioning
# --------------------------------------------------------------------------
def conditioning_card():
    th = np.linspace(0, np.pi, 400)
    lam = periodic_symbol(th)
    ns = [24, 48, 96, 192, 384]
    kp, kc = [], []
    for n in ns:
        ev = np.linalg.eigvalsh(periodic_mass(n, 1.0))
        kp.append(ev.max() / ev.min())
        ev = np.linalg.eigvalsh(clamped_mass(n, 1.0))
        kc.append(ev.max() / ev.min())

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(8, 3.6))
    a1.plot(th, lam, color=SLATE, lw=2.4)
    a1.set_xlabel(r"wavenumber $\theta = k h$")
    a1.set_ylabel(r"eigenvalue / $h$")
    a1.set_xticks([0, np.pi / 2, np.pi], ["0", r"$\pi/2$", r"$\pi$"])
    a1.set_ylim(0, 1.08)
    a1.set_title("periodic symbol (FFT divisor)", fontsize=10.5)
    a1.annotate(f"min {lam[-1]:.4f}", (np.pi, lam[-1]), xytext=(-70, 22),
                textcoords="offset points", fontsize=9, color=INK,
                arrowprops=dict(arrowstyle="-", color=GRAY))

    a2.semilogx(ns, kp, "o-", color=SLATE, lw=2.2, label="periodic (x, z)")
    a2.semilogx(ns, kc, "s-", color=AMBER, lw=2.2, label="clamped (y)")
    a2.set_xticks(ns, [str(n) for n in ns])
    a2.minorticks_off()
    a2.set_xlabel("basis functions per direction")
    a2.set_ylabel(r"$\kappa_2(M)$")
    a2.set_ylim(0, max(kc) * 1.35)
    a2.legend(loc="upper left", fontsize=9)
    a2.set_title("condition number under refinement", fontsize=10.5)
    for n, v in zip(ns[-1:], kp[-1:]):
        a2.text(n, v - 3.2, f"{v:.1f}", ha="center", fontsize=9, color=SLATE)
    for n, v in zip(ns[-1:], kc[-1:]):
        a2.text(n, v + 1.5, f"{v:.1f}", ha="center", fontsize=9, color="#9A6A1E")
    fig.suptitle(f"Mass-matrix condition number is mesh independent ({kp[-1]:.1f} periodic)",
                 x=0.02, ha="left", fontsize=14, fontweight="bold", color=INK)
    fig.text(0.02, 0.0,
             "Cubic B-spline Galerkin mass matrices built with the solver's quadrature; "
             "the bounded condition number bounds the Krylov iteration count under grid refinement.",
             fontsize=7.6, color=GRAY)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(OUT / "card_mass_conditioning.png")
    plt.close(fig)
    return kp, kc


if __name__ == "__main__":
    hero()
    stack_card()
    nnz, nnz3, fac, full = storage_card()
    kp, kc = conditioning_card()
    print("1D nnz (x, y, z):", nnz, "sum", sum(nnz))
    print("3D Kronecker nnz:", nnz3, f"factored {fac} B, assembled {full / 1e9:.2f} GB")
    print("kappa periodic:", [round(float(v), 3) for v in kp])
    print("kappa clamped: ", [round(float(v), 3) for v in kc])
