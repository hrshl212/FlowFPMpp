"""Plot a FlowFPM++ example output (demo_cylinder_vortices or flapping_foil_fpm):
vorticity and the x/y vortex-induced force (F_VIF) densities phi_i div[(u.grad)u],
with influence-potential
contour lines (dashed where negative) on the force panels.
Usage: python plot_fields.py IN.bin OUT.png [x0 x1 y0 y1 [vorticity title [lamb]]]"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

path, out = sys.argv[1], sys.argv[2]
window = [float(v) for v in sys.argv[3:7]] if len(sys.argv) >= 7 else [-1.5, 5.5, -1.6, 1.6]
vorticity_title = sys.argv[7] if len(sys.argv) >= 8 else "Vorticity (synthetic vortex street)"
# Optional 8th argument "lamb": the density files hold -(omega x u).grad(phi_i) instead.
lamb = len(sys.argv) >= 9 and sys.argv[8] == "lamb"
with open(path, "rb") as f:
    nx, ny = np.fromfile(f, np.int32, 2)
    xc = np.fromfile(f, float, nx)
    yc = np.fromfile(f, float, ny)
    oz, fx, fy, phix, phiy, solid = [np.fromfile(f, float, nx * ny).reshape(ny, nx) for _ in range(6)]

ink = "#3b3b38"


def colour_limit(values, floor=0.02, pct=80):
    """Symmetric colour limit set by the flow structures rather than by a few
    near-wall peaks: the pct-th percentile of |value| over the cells whose
    magnitude exceeds `floor` times the 99.9th percentile. Values beyond the
    limit saturate."""
    a = np.abs(values[np.isfinite(values)])
    top = np.percentile(a, 99.9)
    return np.percentile(a[a > floor * top], pct)

panels = [
    (oz, None, r"$\omega_z$", vorticity_title),
    (fx, phix, r"$-(\omega\times u)\cdot\nabla\phi_x$" if lamb else r"$\phi_x\,\nabla\cdot[(u\cdot\nabla)u]$",
     (r"Lamb-vector force density, x (drag), with $\phi_x$" if lamb else r"$F_{VIF}$ density, x (drag), with $\phi_x$")),
    (fy, phiy, r"$-(\omega\times u)\cdot\nabla\phi_y$" if lamb else r"$\phi_y\,\nabla\cdot[(u\cdot\nabla)u]$",
     (r"Lamb-vector force density, y (lift), with $\phi_y$" if lamb else r"$F_{VIF}$ density, y (lift), with $\phi_y$")),
]
fig, axes = plt.subplots(3, 1, figsize=(7.5, 8.4), sharex=True, constrained_layout=True)
for ax, (field, potential, label, title) in zip(axes, panels):
    fluid = np.where(solid > 0, np.nan, field)
    sel = np.ix_((yc >= window[2]) & (yc <= window[3]), (xc >= window[0]) & (xc <= window[1]))
    lim = colour_limit(fluid[sel])
    levels = np.linspace(-lim, lim, 34)  # 33 bins, the middle one centred on zero
    cs = ax.contourf(xc, yc, np.clip(fluid, -lim, lim), levels=levels, cmap="RdBu_r")
    if potential is not None:
        p = np.where(solid > 0, np.nan, potential)
        inside = (xc >= window[0]) & (xc <= window[1])
        rows = (yc >= window[2]) & (yc <= window[3])
        plim = np.nanpercentile(np.abs(p[np.ix_(rows, inside)]), 99)
        ax.contour(xc, yc, p, levels=np.linspace(-plim, plim, 10), colors="#6f6f6b", linewidths=0.5)
    ax.contourf(xc, yc, solid, levels=[0.5, 1.5], colors=["#8a8a86"])
    ax.set_aspect("equal")
    ax.set_xlim(window[0], window[1])
    ax.set_ylim(window[2], window[3])
    ax.set_title(title, loc="left", fontsize=10.5, color=ink)
    ax.set_ylabel("y", color=ink)
    ax.tick_params(colors=ink, labelsize=8)
    for s in ax.spines.values():
        s.set_color("#b5b5b0")
    cb = fig.colorbar(cs, ax=ax, shrink=0.92, ticks=np.linspace(-lim, lim, 5), format="%.2g")
    cb.set_label(label, color=ink, fontsize=9)
    cb.ax.tick_params(labelsize=7, colors=ink)
    cb.outline.set_visible(False)
axes[-1].set_xlabel("x", color=ink)
fig.savefig(out, dpi=150, facecolor="white")
