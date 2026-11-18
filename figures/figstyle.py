"""Figure style for the Figure 5 panels.

Trimmed from the style module shared by every figure in the manuscript, so
panels rendered here match the published ones in font, weight and size. Panels
render inline; nothing is written. Panel letters are added at assembly.
"""

import matplotlib as mpl
import matplotlib.pyplot as plt  # noqa: F401  (re-exported for the notebook)

FONT_SCALE = 1.25


def fs(size):
    return size * FONT_SCALE


# One colour and one marker per antibody, constant across every panel. The
# colours are shades within that antibody's antigenic-region hue (FRLY blue,
# betaSW pink, BKLY wheat), so a panel says which antibody a point belongs to
# and which region it targets at the same time.
AB_COLOR = {
    "HEPC108": "#A8C0DE", "hcab55": "#6C8EBF", "hcab64": "#3F5F94",   # FRLY
    "hcab17": "#F0C3DA", "hcab44": "#D97BA9",                         # betaSW
    "hcab23": "#F0D6A6", "hcab40": "#C99A45",                         # BKLY
}
AB_MARKER = {"HEPC108": "o", "hcab55": "s", "hcab64": "^",
             "hcab17": "D", "hcab44": "v",
             "hcab23": "P", "hcab40": "X"}

# print sizes (inches)
W_PANEL, H_PANEL = 2.45, 1.80

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Liberation Sans", "Helvetica", "DejaVu Sans"],
    "font.weight": "bold",
    "axes.labelweight": "bold",
    "axes.titleweight": "bold",
    "font.size": fs(7),
    "axes.labelsize": fs(7.5),
    "axes.titlesize": fs(8),
    "xtick.labelsize": fs(7),
    "ytick.labelsize": fs(7),
    "legend.fontsize": fs(6.5),
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 1.3,
    "xtick.major.width": 1.3,
    "ytick.major.width": 1.3,
    "xtick.major.size": 3.0,
    "ytick.major.size": 3.0,
    "lines.linewidth": 1.6,
    "lines.markeredgewidth": 0.9,
    "figure.dpi": 150,
})
