import os, json, h5py
import csv
import base64
import numpy as np
import pandas as pd
import pickle
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

output = "../ZeroNuFit-dev/output_correct_inputs_uncorrelated_nuisances/v3/fit_9_l200_uniform_1BI_CorrEff"#config["output"]
fit_results = json.load(open(os.path.join(output, "mcmc_files/fit_results.json")))["quantile90"]
limit_signal = fit_results["S"]
limit_half_life = round(10/(limit_signal),1)

with h5py.File('ovbb_plot_entries.h5', 'r') as f:
    bkg_shape = f['bkg_shape'][()].decode('utf-8')
    exp_tot = round(f['exposure'][()],1)
    events = f['events'][:]
    energies = f['energies'][:]
    energies_all_window = f['energies_all_window'][:]
    b_68_left = f['b_68_left'][:]
    b_68_right = f['b_68_right'][:]
    b_mode = f['b_mode'][:]
    b_mode_all_window = f['b_mode_all_window'][:]
    s_90_left = f['s_90_left'][:]
    s_90_right = f['s_90_right'][:]

print(len(b_68_left))
print(len(s_90_left))
print((b_68_left))
print((b_68_right))
print((s_90_left))
print((s_90_right))

fig, ax = plt.subplots(figsize = (5, 2.5))

# ...background 68%
if bkg_shape == "flat":
    band = plt.fill_between(
        [1930, 2190],
        np.full(2, b_68_left[0]),
        np.full(2, b_68_right[0]),
        linewidth = 0,
        alpha = 0.5,
        color = "#228833",
        zorder = 0,
        label="Background best fit and 68% C.I. interval",
    )
else:
    band = plt.fill_between(energies_all_window, b_68_left, b_68_right, color="#228833", alpha=0.5, linewidth=0, zorder=0, label="Background best fit and 68% C.I. interval")

# ...background best fit
line = plt.plot([],[],color = "#228833")
if bkg_shape == "flat":
    line, = plt.plot([1930, 2190], [b_mode[0], b_mode[0]], color = "#228833", zorder = 1)
if bkg_shape != "flat" and b_mode_all_window!=[0]:
    line, = plt.plot(energies_all_window, b_mode_all_window, color = "#228833", zorder = 1)

# ...signal 90%
limit_area = plt.fill_between(
    energies,
    s_90_left,
    s_90_right,
    linewidth = 0,
    color = "#4477AA",
    label = r"$T^{0\nu}_{1/2} > $" + str(limit_half_life) +r"$\times 10^{26}$ yr [90% C.I.]",
)

# plot the data
unbinned = True

# get energies of surviving events
x = events

if unbinned:
    data_art, stemlines, baseline = ax.stem(x, np.full(len(x), 1e-3), "#CC3311", basefmt = "none", label=f"LEGEND-200 data [{exp_tot} kg yr]")
    data_art.set_markersize(1)
    # NOTE: comment the following if the stem has a length!
    data_art.set_zorder(99)
    data_art.set_clip_on(False)

    stemlines.set_linewidth(0.6)
else:
    w = 1 / exp_tot
    bin_width = 1
    xmin = np.min(1930)
    xmax = np.max(2190)
    nbins = np.arange(xmin, xmax + bin_width, bin_width)
    _, _, data_art = plt.hist(x, weights = np.full(len(x), w), bins = nbins, color = "#CC3311", label=f"LEGEND-200 data [{exp_tot} kg yr]")

ax.set_ylim(1E-4, 2)
ax.set_yscale("log")
plt.axvspan(2099, 2109, facecolor="black", alpha=0.4)
plt.axvspan(2114, 2124, facecolor="black", alpha=0.4)
ax.set_xlim(1930, 2190)
ax.set_xlabel("Energy [keV]", loc = "right")
ax.set_ylabel("Counts / (keV kg yr)")
ax.legend(loc='upper left', fontsize=8)
plt.tight_layout()
plt.savefig("l200_result.pdf", bbox_inches = "tight")