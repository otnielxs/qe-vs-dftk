import matplotlib.pyplot as plt


def apply_style():
    plt.rcParams.update({
        "font.size": 12,
        "axes.titlesize": 13,
        "axes.labelsize": 12,
        "xtick.labelsize": 11,
        "ytick.labelsize": 11,
        "legend.fontsize": 10,
        "lines.linewidth": 1.8,
        "lines.markersize": 6,
        "savefig.dpi": 200,
        "savefig.bbox": "tight",
    })
