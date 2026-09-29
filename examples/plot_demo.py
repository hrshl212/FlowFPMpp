"""Plot the demo_cylinder_vortices output: vorticity and the x/y vortex-induced
force densities, with influence-potential contour lines (dashed where
negative) on the force panels. Usage: python plot_demo.py DEMO.bin OUT.png"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

path, out = sys.argv[1], sys.argv[2]
with open(path, "rb") as f:
    nx, ny = np.fromfile(f, np.int32, 2)
    xc = np.fromfile(f, float, nx)
    yc = np.fromfile(f, float, ny)
    oz, fx, fy, phix, phiy, solid = [np.fromfile(f, float, nx * ny).reshape(ny, nx) for _ in range(6)]

ink = "#3b3b38"
panels = [
    (oz, None, r"$\omega_z$", "Vorticity (synthetic vortex street)"),
    (fx, phix, r"$-(\omega\times u)\cdot\nabla\phi_x$", r"Vortex-induced force density, x (drag), with $\phi_x$"),
    (fy, phiy, r"$-(\omega\times u)\cdot\nabla\phi_y$", r"Vortex-induced force density, y (lift), with $\phi_y$"),
]
fig, axes = plt.subplots(3, 1, figsize=(7.5, 8.4), sharex=True, constrained_layout=True)
for ax, (field, potential, label, title) in zip(axes, panels):
    fluid = np.where(solid > 0, np.nan, field)
    lim = np.nanpercentile(np.abs(fluid), 99.7)
    levels = np.linspace(-lim, lim, 34)  # 33 bins, the middle one centred on zero
    cs = ax.contourf(xc, yc, np.clip(fluid, -lim, lim), levels=levels, cmap="RdBu_r")
    if potential is not None:
        p = np.where(solid > 0, np.nan, potential)
        ax.contour(xc, yc, p, levels=np.linspace(-0.45, 0.45, 10), colors="#6f6f6b", linewidths=0.5)
    ax.contourf(xc, yc, solid, levels=[0.5, 1.5], colors=["#8a8a86"])
    ax.set_aspect("equal")
    ax.set_xlim(-1.5, 5.5)
    ax.set_ylim(-1.6, 1.6)
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
