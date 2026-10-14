import matplotlib.pyplot as plt


def apply_style():
    plt.rcParams.update({
       "font.family": ["DejaVu Serif"],
        "font.serif": ["Times New Roman"],
        "font.size": 12,
        "axes.titlesize": 15,
        "axes.labelsize": 12,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 10,
        "lines.linewidth": 2.0,
        "lines.markersize": 5,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
    })
