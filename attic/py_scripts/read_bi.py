### read_bi.py
#
# Authors: Sofia Calgaro, Toby Dixon
# 
###
import os,json,h5py,math
import numpy as np
import matplotlib
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def smallest_68_ci(bkg, nbins):
    counts, bin_edges = np.histogram(bkg, bins=nbins, density=True)
    max_bin_index = np.argmax(counts)
    max_bin_index = np.argmax(counts)
    max_bin_count = counts[max_bin_index]
    max_bin_start = bin_edges[max_bin_index]
    max_bin_end = bin_edges[max_bin_index + 1]
    marg_mode = (max_bin_end+max_bin_start)/2
    total_area = np.sum(counts * np.diff(bin_edges)) 
    
    target_area = 0.68 * total_area
    
    accumulated_area = 0.0
    left_index = max_bin_index
    right_index = max_bin_index
    bin_widths = np.diff(bin_edges)
    
    while accumulated_area < target_area:
        left_area = counts[left_index] * bin_widths[left_index] if left_index > 0 else 0
        right_area = counts[right_index] * bin_widths[right_index] if right_index < len(counts) - 1 else 0

        if left_area >= right_area and left_index > 0:
            accumulated_area += left_area
            left_index -= 1
        elif right_index < len(counts) - 1:
            accumulated_area += right_area
            right_index += 1
        else:
            break
    
    ci_low = bin_edges[left_index]
    ci_high = bin_edges[right_index + 1]
    
    return marg_mode, ci_low, ci_high

bin_width = 0.000025
xmin = np.min(0)
xmax = np.max(0.002)
nbins = np.arange(xmin, xmax + bin_width, bin_width)

with h5py.File('../bi_tot_all.h5', 'r') as f:
    bkg = f['background'][:]

marg_mode, ci_low, ci_high = smallest_68_ci(bkg, nbins)

bkg_name ="B_l200a_all"
fig, ax = plt.subplots(figsize=(5,3.3))
plt.hist(bkg, bins=nbins, histtype='stepfilled', density=True, color='forestgreen', alpha=0.15, zorder=-12)
plt.hist(bkg, bins=nbins, histtype='step', density=True, color='forestgreen', lw=1, zorder=-12)
#plt.axvline(gms[bkg_name], color='navy', linestyle="-", label="Global mode")
counts, bin_edges = np.histogram(bkg, bins=nbins, density=True)
max_bin_index = np.argmax(counts)
max_bin_count = counts[max_bin_index]
max_bin_start = bin_edges[max_bin_index]
max_bin_end = bin_edges[max_bin_index + 1]
marginalized_mode = (max_bin_end+max_bin_start)/2*1e4
plt.axvline((max_bin_end+max_bin_start)/2, color='r', linestyle="--", label="Marginalized mode")
plt.axvline(ci_low, color='r', linestyle=":", label="Smallest 68% CI")
plt.axvline(ci_high, color='r', linestyle=":")
ax.set_xticks([0,0.0005,0.001,0.0015,0.002])
plt.xlabel(r'BI$_{\rm tot}$ (ckky)')
plt.ylabel('Probability density')
plt.legend(loc="upper right", frameon=False)
plt.savefig("bi_tot.pdf",bbox_inches='tight')