### plot_marginalized_posteriors.py
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

### some colours
c_pistachio = (0.58, 0.87, 0.45)
c_columbiablue = (0.61, 0.87, 1.0)
c_dark_columbiablue = (0.61, 0.87*0.85, 1)
c_frenchblue = (0.0, 0.45, 0.73)
c_dark_gray = (0.33, 0.33, 0.33)
c_dark_pistachio = (0.58*0.85, 0.87*0.85, 0.45*0.85)
c_lava = (0.94, 0.01, 0.05)

### function to read output h5 files with saved posteriors
def read_samples(path,par):
    out=None
    if os.path.isfile(path):
        with h5py.File(path, 'r') as hdf:
            
            v = hdf[f"v/{par}"][:]
            w = hdf[f"weight"][:]
            
            out= np.repeat(v, w)
        return out
    else:
        print(f"File {path} not found!")
        exit()
    

bin_width = 0.000025
xmin = np.min(0)
xmax = np.max(0.1)
nbins = np.arange(xmin, xmax + bin_width, bin_width)
bin_width_qbb = 0.0005
xmax = np.max(0.2)
nbins_qbb = np.arange(xmin, xmax + bin_width_qbb, bin_width_qbb)
bin_width_S = 0.5
xmin_S = np.min(0)
xmax_S = np.max(110)
nbins_S = np.arange(xmin_S, xmax_S + bin_width_S, bin_width_S)
bin_width_S_low = 0.25
nbins_S_low = np.arange(xmin_S, xmax_S + bin_width_S_low, bin_width_S_low)
bin_width_slope = 0.1
xmin_slope = np.min(-30)
xmax_slope = np.max(50)
nbins_slope = np.arange(xmin_slope, xmax_slope + bin_width_slope, bin_width_slope)
nbins_slope_more = np.arange(xmin_slope, xmax_slope + bin_width_slope, 0.05)
bin_width_slope_low = 0.1
nbins_slope_low = np.arange(xmin_slope, xmax_slope + bin_width_slope_low, bin_width_slope_low)

gm_type = "marginalized_modes"

center = 1930
expo = 42.25
fit_range = [[1930, 2099], [2109, 2114], [2124, 2190]]
Qbb = 2039.06

def normalize():
    range_l = [arr[0] for arr in fit_range] 
    range_h = [arr[1] for arr in fit_range] 
    norm = sum([h - l for h, l in zip(range_h, range_l)])

    return norm

def compute_deltaE():
    deltaE = sum([arr[2] - arr[1] for arr in fit_range])
    return deltaE

def compute_window():
    window = fit_range[-1][-1] - fit_range[0][0]
    return window

def normalize_sq():
    range_l = [arr[0] for arr in fit_range] 
    range_h = [arr[1] for arr in fit_range] 
    window_sq = sum([h**2 - l**2 for h, l in zip(range_h, range_l)])
    return window_sq

def exp_stable(x):
    if abs(x) < 1e-6:
        return 1 + x + (x**2) / 2 + (x**3) / 6
    else:
        return math.exp(x)

norm = normalize()
window = compute_window()
window_sq = normalize_sq()
sum_range = normalize()
sum_range_sq = normalize_sq()
print("norm:", norm)
print("window:", window)

def global_mode_at_qbb_linear(global_mode, slope_list_gm):

    bkg_at_qbb_list = 0
    b = global_mode
    s = slope_list_gm
    
    term_1 = (b * expo * norm)
    
    norm_all = sum_range * (1 - s * center / window) + s * sum_range_sq / (2 * window)
    term_2 = (1 + s*(Qbb-center)/window) / norm_all
    
    bkg_at_qbb_list = term_1 * term_2
        
    return bkg_at_qbb_list 

def global_mode_at_qbb(global_mode):
    bkg_at_qbb_list = 0
    b = global_mode
    term_1 = (b * expo * norm)
    term_2 =  1 / norm
    
    bkg_at_qbb_list = term_1 * term_2
        
    return bkg_at_qbb_list 

def global_mode_at_qbb_exponential(global_mode, slope_list_gm):
    range_l = [arr[0] for arr in fit_range] 
    range_h = [arr[1] for arr in fit_range]  
    centers = [center, center, center]
    
    bkg_at_qbb_list = 0
    b = global_mode
    R = slope_list_gm
    Rt = R / window

    term_1 = (b * expo * norm)
    
    if abs(Rt) > 1e-6:
        norm_all = (-sum([exp_stable(-Rt * (center - l)) for l in range_l]) +
                sum([exp_stable(-Rt * (center - h)) for h in range_h])) / Rt
    else:
        norm_all = normalize()
    term_2 = exp_stable((Qbb-center)*Rt) / norm_all
    
    bkg_at_qbb_list = term_1 * term_2
    
    return bkg_at_qbb_list 

def bkg_at_qbb_linear(bkg_list, slope_list):

    bkg_at_qbb_list = []
    for idx,b in enumerate(bkg_list):
        s = slope_list[idx]
        
        term_1 = (b * expo * norm)
        
        norm_all = sum_range * (1 - s * center / window) + s * sum_range_sq / (2 * window)
        term_2 = (1 + s*(Qbb-center)/window) / norm_all
        
        bkg_at_qbb_list.append(term_1 * term_2)
        
    return bkg_at_qbb_list 



def bkg_at_qbb_uniform(bkg_list):
    bkg_at_qbb_list = []
    for idx,b in enumerate(bkg_list):
        term_1 = (b * expo * norm)
        term_2 =  1 / norm
        bkg_at_qbb_list.append(term_1 * term_2)
    return bkg_at_qbb_list 


def bkg_at_qbb_exponential(bkg_list, slope_list):
    range_l = [arr[0] for arr in fit_range] 
    range_h = [arr[1] for arr in fit_range]  
    centers = [center, center, center]
    
    bkg_at_qbb_list = []
    for idx,b in enumerate(bkg_list):
        R = slope_list[idx]
        Rt = R / window

        term_1 = (b * expo * norm)
        
        if abs(Rt) > 1e-6:
            norm_all = (-sum([exp_stable(-Rt * (center - l)) for l in range_l]) +
                    sum([exp_stable(-Rt * (center - h)) for h in range_h])) / Rt
        else:
            norm_all = normalize()
        term_2 = exp_stable((Qbb-center)*Rt) / norm_all
        
        bkg_at_qbb_list.append(term_1 * term_2)
    
    return bkg_at_qbb_list 


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




def different_bkg_models():
    with PdfPages(f"comparison_of_different_bkg_shapes.pdf") as pdf:
        """
        # l200
        flat_file = f"../output_correct_inputs_uncorrelated_nuisances/v2/fit_9_l200_uniform_1BI_CorrEff/mcmc_files/samples.h5"
        linear_file = f"../output_correct_inputs_uncorrelated_nuisances/v2/fit_9_l200_uniform_1BI_CorrEff_linear/mcmc_files/samples.h5"
        expo_file = f"../output_correct_inputs_uncorrelated_nuisances/v2/fit_9_l200_uniform_1BI_CorrEff_exponential/mcmc_files/samples.h5"
        flat_json = json.load(open(f"../output_correct_inputs_uncorrelated_nuisances/v2/fit_9_l200_uniform_1BI_CorrEff/mcmc_files/fit_results.json"))
        linear_json = json.load(open(f"../output_correct_inputs_uncorrelated_nuisances/v2/fit_9_l200_uniform_1BI_CorrEff_linear/mcmc_files/fit_results.json"))
        expo_json = json.load(open(f"../output_correct_inputs_uncorrelated_nuisances/v2/fit_9_l200_uniform_1BI_CorrEff_exponential/mcmc_files/fit_results.json"))
        bkg_name = "B_l200a_all"
        # l200
        flat_file = f"../output_correct_inputs_uncorrelated_nuisances/v3/fit_9_l200_uniform_1BI_CorrEff/mcmc_files/samples.h5"
        linear_file = f"../output_correct_inputs_uncorrelated_nuisances/v3/fit_9_l200_uniform_1BI_CorrEff_linear/mcmc_files/samples.h5"
        expo_file = f"../output_correct_inputs_uncorrelated_nuisances/v3/fit_9_l200_uniform_1BI_CorrEff_exponential/mcmc_files/samples.h5"
        flat_json = json.load(open(f"../output_correct_inputs_uncorrelated_nuisances/v3/fit_9_l200_uniform_1BI_CorrEff/mcmc_files/fit_results.json"))
        linear_json = json.load(open(f"../output_correct_inputs_uncorrelated_nuisances/v3/fit_9_l200_uniform_1BI_CorrEff_linear/mcmc_files/fit_results.json"))
        expo_json = json.load(open(f"../output_correct_inputs_uncorrelated_nuisances/v3/fit_9_l200_uniform_1BI_CorrEff_exponential/mcmc_files/fit_results.json"))
        bkg_name = "B_l200a_Nu24"
        # gerda
        flat_file = f"../../ZeroNuFit.jl/output_1/output_old/wrong_expo_l200/fit_1_gerda_phII_uniform_1BI_CorrEff/mcmc_files/samples.h5"
        linear_file = f"../../ZeroNuFit.jl/output/fit_3_gerda_linear_CorrEff/mcmc_files/samples.h5"
        expo_file = f"../../ZeroNuFit.jl/output/fit_4_gerda_expo_CorrEff/mcmc_files/samples.h5"
        flat_json = json.load(open(f"../../ZeroNuFit.jl/output_1/output_old/wrong_expo_l200/fit_1_gerda_phII_uniform_1BI_CorrEff/mcmc_files/fit_results.json"))
        linear_json = json.load(open(f"../../ZeroNuFit.jl/output/fit_3_gerda_linear_CorrEff/mcmc_files/fit_results.json"))
        expo_json = json.load(open(f"../../ZeroNuFit.jl/output/fit_4_gerda_expo_CorrEff/mcmc_files/fit_results.json"))
        bkg_name = "B_gerda_all_pII"
        """
        
        # combined
        flat_file = f"../output_correct_inputs_uncorrelated_nuisances/v3/fit_11_gerdaIandII_5BI_l200_1BI_mjd_uniform_1BI_CorrEff/mcmc_files/samples.h5"
        linear_file = f"../output_correct_inputs_uncorrelated_nuisances/v3/fit_11_gerdaIandII_5BI_l200_1BI_mjd_uniform_1BI_CorrEff_linear/mcmc_files/samples.h5"
        expo_file = f"../../ZeroNuFit-dev-ovbb/output_correct_inputs_uncorrelated_nuisances/v3/fit_11_gerdaIandII_5BI_l200_1BI_mjd_uniform_1BI_CorrEff_exponential/mcmc_files/samples.h5"
        flat_json = json.load(open(f"../output_correct_inputs_uncorrelated_nuisances/v3/fit_11_gerdaIandII_5BI_l200_1BI_mjd_uniform_1BI_CorrEff/mcmc_files/fit_results.json"))
        linear_json = json.load(open(f"../output_correct_inputs_uncorrelated_nuisances/v3/fit_11_gerdaIandII_5BI_l200_1BI_mjd_uniform_1BI_CorrEff_linear/mcmc_files/fit_results.json"))
        expo_json = json.load(open(f"../../ZeroNuFit-dev-ovbb/output_correct_inputs_uncorrelated_nuisances/v3/fit_11_gerdaIandII_5BI_l200_1BI_mjd_uniform_1BI_CorrEff_exponential/mcmc_files/fit_results.json"))
        custom_colors = [
            "deepskyblue",
            "dodgerblue", 
            "blue",  
            "navy",  
            "darkviolet", 
            "teal",  
            "#228b22", 
            "#ff7f0e", 
            "chocolate",  
            "maroon",  
        ]

        for shape in ["linear","exponential"]:
            fig, ax = plt.subplots(figsize=(8,6))
            for idx,bkg_name in enumerate(["ph1_golden", "ph1_bege", "ph1_silver", "ph1_extra", "B_gerda_all_pII", "B_l200a_Nu24", "B_l200a_extra", "mjd-DS0", "mjd-mod1", "mjd-mod2"]):

                if shape == "linear":
                    slope_to_plot = read_samples(linear_file, f"{bkg_name}_slope")
                    mybins=nbins_slope_more
                else:
                    slope_to_plot =read_samples(expo_file, f"{bkg_name}_slope")
                    mybins=nbins_slope
                #plt.hist(slope_to_plot, bins=nbins_slope, histtype='stepfilled', density=True, color=c_frenchblue, alpha=0.15, zorder=-10)
                plt.hist(slope_to_plot, bins=mybins, histtype='step', density=True, lw=1.2, label=bkg_name, color=custom_colors[idx])

            if shape=="linear":plt.yscale("log")
            if shape=="linear":ax.legend(loc='upper right', ncol=2)
            if shape=="exponential":ax.legend(loc='upper left', ncol=2)
            if shape=="linear": plt.xlabel(r"Slope m$_{\rm lin}$")
            if shape=="exponential": plt.xlabel(r"Slope m$_{\rm exp}$")
            plt.ylabel('Probability density')
            #plt.title(f"Comparison of bkg slope posteriors for {shape} shape")
            if shape=="linear": plt.xlim(-1,10)
            if shape=="exponential": plt.xlim(-15,5)
            pdf.savefig(bbox_inches='tight')
            plt.close()
            print("plotted comparison of posterior pdfs for B SLOPE")

        count_below_zero = sum(x < 0 for x in slope_linear)
        percentage_below_zero = (count_below_zero / len(slope_linear)) * 100
        percentage_above_zero = 100-percentage_below_zero
        print("% of events that have a slope in the linear fit <0: ", percentage_below_zero)
        print("% of events >0: ", percentage_above_zero)
        count_below_zero = sum(x < 0 for x in slope_expo)
        percentage_below_zero = (count_below_zero / len(slope_expo)) * 100
        percentage_above_zero = 100-percentage_below_zero
        print("% of events that have a slope in the exponential fit <0: ", percentage_below_zero)
        print("% of events >0: ", percentage_above_zero)


        fig, ax = plt.subplots(figsize=(5,3.3))
        signal_flat = read_samples(flat_file, 'S')
        plt.hist(signal_flat, bins=nbins_S, histtype='stepfilled', density=True, color='darkorange', alpha=0.15, zorder=-10)
        plt.hist(signal_flat, bins=nbins_S, histtype='step', density=True, color='darkorange', lw=1, zorder=-10)
        
        signal_linear = read_samples(linear_file, 'S')
        plt.hist(signal_linear, bins=nbins_S, histtype='stepfilled', density=True, color=c_frenchblue, alpha=0.15, zorder=-10)
        plt.hist(signal_linear, bins=nbins_S, histtype='step', density=True, color=c_frenchblue, lw=1, zorder=-10)
        
        signal_expo = read_samples(expo_file, 'S')
        plt.hist(signal_expo, bins=nbins_S, histtype='stepfilled', density=True, color=c_dark_pistachio, alpha=0.15, zorder=-8)
        plt.hist(signal_expo, bins=nbins_S, histtype='step', density=True, color=c_dark_pistachio, lw=1, zorder=-8)
        
        flt_patch = mpatches.Patch(color='darkorange', label='Flat bkg', alpha=0.4)
        lin_patch = mpatches.Patch(color=c_frenchblue, label='Linear bkg', alpha=0.4)
        exp_patch = mpatches.Patch(color=c_dark_pistachio, label='Exponential bkg', alpha=0.4)
        ax.legend(handles=[flt_patch, lin_patch, exp_patch], loc='upper right', ncol=1, frameon=False)
        plt.xlabel(r'S ($10^{-27}$ yr$^{-1}$)')
        plt.ylabel('Probability density')
        plt.title("Comparison of signal posteriors")
        plt.xlim(0,100)
        #plt.yscale("log")
        pdf.savefig(bbox_inches='tight')
        plt.close()
        print("plotted comparison of posterior pdfs for S")
        
        fig, ax = plt.subplots(figsize=(5,3.3))
        plt.hist(signal_flat, bins=nbins_S_low, histtype='stepfilled', density=True, color='darkorange', alpha=0.15, zorder=-10)
        plt.hist(signal_flat, bins=nbins_S_low, histtype='step', density=True, color='darkorange', lw=1, zorder=-10)
        plt.hist(signal_linear, bins=nbins_S_low, histtype='stepfilled', density=True, color=c_frenchblue, alpha=0.15, zorder=-10)
        plt.hist(signal_linear, bins=nbins_S_low, histtype='step', density=True, color=c_frenchblue, lw=1, zorder=-10)
        plt.hist(signal_expo, bins=nbins_S_low, histtype='stepfilled', density=True, color=c_dark_pistachio, alpha=0.15, zorder=-8)
        plt.hist(signal_expo, bins=nbins_S_low, histtype='step', density=True, color=c_dark_pistachio, lw=1, zorder=-8)
        
        flt_patch = mpatches.Patch(color='darkorange', label='Flat bkg', alpha=0.4)
        lin_patch = mpatches.Patch(color=c_frenchblue, label='Linear bkg', alpha=0.4)
        exp_patch = mpatches.Patch(color=c_dark_pistachio, label='Exponential bkg', alpha=0.4)
        ax.legend(handles=[flt_patch, lin_patch, exp_patch], loc='upper right', ncol=1, frameon=False)
        plt.xlabel(r'S ($10^{-27}$ yr$^{-1}$)')
        plt.ylabel('Probability density')
        #plt.title("Comparison of signal posteriors - zoom")
        plt.xlim(0,20)
        #plt.yscale("log")
        pdf.savefig(bbox_inches='tight')
        plt.close()
        print("plotted comparison of posterior pdfs for S (zoom)")

        fig, ax = plt.subplots(figsize=(5,3.3))
        bkg_flat = read_samples(flat_file, f'{bkg_name}')
        plt.hist(bkg_flat, bins=nbins, histtype='stepfilled', density=True, color='darkorange', alpha=0.15, zorder=-10)
        plt.hist(bkg_flat, bins=nbins, histtype='step', density=True, color='darkorange', lw=1, zorder=-10)
        
        bkg_linear = read_samples(linear_file, f'{bkg_name}')
        plt.hist(bkg_linear, bins=nbins, histtype='stepfilled', density=True, color=c_frenchblue, alpha=0.15, zorder=-10)
        plt.hist(bkg_linear, bins=nbins, histtype='step', density=True, color=c_frenchblue, lw=1, zorder=-10)
        
        bkg_expo = read_samples(expo_file, f'{bkg_name}')
        plt.hist(bkg_expo, bins=nbins, histtype='stepfilled', density=True, color=c_dark_pistachio, alpha=0.15, zorder=-8)
        plt.hist(bkg_expo, bins=nbins, histtype='step', density=True, color=c_dark_pistachio, lw=1, zorder=-8)
        
        flt_patch = mpatches.Patch(color='darkorange', label='Flat bkg', alpha=0.4)
        lin_patch = mpatches.Patch(color=c_frenchblue, label='Linear bkg', alpha=0.4)
        exp_patch = mpatches.Patch(color=c_dark_pistachio, label='Exponential bkg', alpha=0.4)
        ax.legend(handles=[flt_patch, lin_patch, exp_patch], loc='upper right', ncol=1, frameon=False)
        plt.xlim(0,40e-4)
        x_ticks = ax.get_xticks()
        ax.set_xticks(x_ticks)
        ax.set_xticklabels((x_ticks * 1e4).astype(int))
        ax.set_xticks(np.arange(0, 30e-4, 10e-4))
        ax.set_xticklabels(np.arange(0, 30e-4, 10e-4).astype(int))
        ax.set_xticks(x_ticks)
        ax.set_xticklabels((x_ticks * 1e4).astype(int))
        if "extra" in bkg_name: plt.xlabel(r'BI$_{\rm extra}$ ($10^{-4}$ ckky)')
        else: plt.xlabel(r'BI$_{\rm Nu24}$ ($10^{-4}$ ckky)')
        #plt.title("Comparison of BI posteriors")
        plt.ylabel('Probability density')
        #plt.yscale("log")
        pdf.savefig( bbox_inches='tight')
        plt.close()
        print("plotted comparison of posterior pdfs for BI")

        fig, ax = plt.subplots(figsize=(5,3.3))
        bkg_flat = read_samples(flat_file, f'{bkg_name}')
        bkg_flat = bkg_at_qbb_uniform(bkg_flat)
        bkg_flat_gm = global_mode_at_qbb(flat_json['refined_global_modes'][bkg_name])
        print("Global mode (flat bkg):", bkg_flat_gm)
        plt.hist(bkg_flat, bins=nbins_qbb, histtype='stepfilled', density=True, color='darkorange', alpha=0.15, zorder=-10)
        plt.hist(bkg_flat, bins=nbins_qbb, histtype='step', density=True, color='darkorange', lw=1, zorder=-10)
        marg_mode, ci_low, ci_high = smallest_68_ci(bkg_flat, nbins_qbb)
        print("Marginalized mode (flat bkg):", marg_mode)
        print(f"68% CI: [{ci_low}, {ci_high}]")
        
        bkg_linear = read_samples(linear_file, f'{bkg_name}')
        slope_linear = read_samples(linear_file, f'{bkg_name}_slope')
        bkg_linear = bkg_at_qbb_linear(bkg_linear, slope_linear)
        bkg_linear_gm = global_mode_at_qbb_linear(linear_json['refined_global_modes'][bkg_name], linear_json['refined_global_modes'][f'{bkg_name}_slope'])
        print("Global mode (linear bkg):", bkg_linear_gm)
        plt.hist(bkg_linear, bins=nbins_qbb, histtype='stepfilled', density=True, color=c_frenchblue, alpha=0.15, zorder=-10)
        plt.hist(bkg_linear, bins=nbins_qbb, histtype='step', density=True, color=c_frenchblue, lw=1, zorder=-10)
        marg_mode, ci_low, ci_high = smallest_68_ci(bkg_linear, nbins_qbb)
        print("Marginalized mode (linear bkg):", marg_mode)
        print(f"68% CI: [{ci_low}, {ci_high}]")

        bkg_expo = read_samples(expo_file, f'{bkg_name}')
        slope_expo = read_samples(expo_file, f'{bkg_name}_slope')
        bkg_expo = bkg_at_qbb_exponential(bkg_expo, slope_expo)
        bkg_expo_gm = global_mode_at_qbb_exponential(expo_json['refined_global_modes'][bkg_name], expo_json['refined_global_modes'][f'{bkg_name}_slope'])
        print("Global mode (expo bkg):", bkg_expo_gm)
        plt.hist(bkg_expo, bins=nbins_qbb, histtype='stepfilled', density=True, color=c_dark_pistachio, alpha=0.15, zorder=-8)
        plt.hist(bkg_expo, bins=nbins_qbb, histtype='step', density=True, color=c_dark_pistachio, lw=1, zorder=-8)
        marg_mode, ci_low, ci_high = smallest_68_ci(bkg_expo, nbins_qbb)
        print("Marginalized mode (expo bkg):", marg_mode)
        print(f"68% CI: [{ci_low}, {ci_high}]")
        
        flt_patch = mpatches.Patch(color='darkorange', label='Flat bkg', alpha=0.4)
        lin_patch = mpatches.Patch(color=c_frenchblue, label='Linear bkg', alpha=0.4)
        exp_patch = mpatches.Patch(color=c_dark_pistachio, label='Exponential bkg', alpha=0.4)
        ax.legend(handles=[flt_patch, lin_patch, exp_patch], loc='upper right', ncol=1, frameon=False)
        if "extra" in bkg_name: plt.xlabel(r'Extra background at Q$_{\beta\beta}$ (cts/keV)')
        else: plt.xlabel(r'Nu24 background at Q$_{\beta\beta}$ (cts/keV)')
        #plt.title("Comparison of bkg posteriors")
        plt.ylabel('Probability density')
        #plt.yscale("log")
        plt.xlim(0,0.2)
        pdf.savefig( bbox_inches='tight')
        plt.close()
        print("plotted comparison of posterior pdfs for bkg at Qbb")

bin_width = 0.000025
xmin = np.min(0)
xmax = np.max(0.002)
nbins = np.arange(xmin, xmax + bin_width, bin_width)
def sig_and_bkg():
    exp = "output_correct_inputs_uncorrelated_nuisances/v3/wrong_background_index/fit_9_l200_uniform_1BI_CorrEff" # change me!!
    
    with PdfPages(f"signal_and_bkg_posteriors.pdf") as pdf:

        bkg_name ="B_l200a_all"
        fig, ax = plt.subplots(figsize=(5,3.3))
        json_file = json.load(open(f"../{exp}/mcmc_files/fit_results.json"))
        gms = json_file['refined_global_modes']
        c68 = json_file['ci_68']
        bi_all = read_samples(f"../{exp}/mcmc_files/samples.h5", bkg_name)
        plt.hist(bi_all, bins=nbins, histtype='stepfilled', density=True, color='dodgerblue', alpha=0.15, zorder=-12)
        plt.hist(bi_all, bins=nbins, histtype='step', density=True, color='dodgerblue', lw=1, zorder=-12)
        plt.axvline(gms[bkg_name], color='navy', linestyle="-", label="Global mode")
        counts, bin_edges = np.histogram(bi_all, bins=nbins, density=True)
        max_bin_index = np.argmax(counts)
        max_bin_count = counts[max_bin_index]
        max_bin_start = bin_edges[max_bin_index]
        max_bin_end = bin_edges[max_bin_index + 1]
        print("Marginalized mode:", (max_bin_end+max_bin_start)/2*1e4)
        plt.axvline((max_bin_end+max_bin_start)/2, color='r', linestyle="--", label="Marginalized mode")
        plt.axvline(c68[bkg_name][0]['left'], color='r', linestyle=":", label="Smallest 68% CI")
        plt.axvline(c68[bkg_name][0]['right'], color='r', linestyle=":")
        plt.xlim(0,20e-4)
        x_ticks = ax.get_xticks()
        ax.set_xticks(x_ticks)
        ax.set_xticklabels((x_ticks * 1e4).astype(int))
        ax.set_xticks(np.arange(0, 20e-4, 5e-4))
        ax.set_xticklabels(np.arange(0, 20e-4, 5e-4).astype(int))
        ax.set_xticks(x_ticks)
        ax.set_xticklabels((x_ticks * 1e4).astype(int))
        if "all" in bkg_name: plt.xlabel(r'BI$_{\rm tot}$ (ckky)')
        else: plt.xlabel(r'BI ($10^{-4}$ ckky)')
        ax.set_xticks([0,0.0005,0.001,0.0015,0.002])
        plt.ylabel('Probability density')
        plt.legend(loc="upper right", frameon=False)
        pdf.savefig(bbox_inches='tight')
        plt.close()

        bkg_name = "B_l200a_extra"
        fig, ax = plt.subplots(figsize=(5,3.3))
        json_file = json.load(open(f"../{exp}/mcmc_files/fit_results.json"))
        gms = json_file['refined_global_modes']
        c68 = json_file['ci_68']
        bi_all = read_samples(f"../{exp}/mcmc_files/samples.h5", bkg_name)
        plt.hist(bi_all, bins=nbins, histtype='stepfilled', density=True, color='dodgerblue', alpha=0.15, zorder=-12)
        plt.hist(bi_all, bins=nbins, histtype='step', density=True, color='dodgerblue', lw=1, zorder=-12)
        plt.axvline(gms[bkg_name], color='navy', linestyle="-", label="Global mode")
        counts, bin_edges = np.histogram(bi_all, bins=nbins, density=True)
        max_bin_index = np.argmax(counts)
        max_bin_count = counts[max_bin_index]
        max_bin_start = bin_edges[max_bin_index]
        max_bin_end = bin_edges[max_bin_index + 1]
        print("Marginalized mode:", (max_bin_end+max_bin_start)/2*1e4)
        plt.axvline((max_bin_end+max_bin_start)/2, color='r', linestyle="--", label="Marginalized mode")
        plt.axvline(c68[bkg_name][0]['left'], color='r', linestyle=":", label="Smallest 68% CI")
        plt.axvline(c68[bkg_name][0]['right'], color='r', linestyle=":")
        plt.xlim(0,30e-4)
        x_ticks = ax.get_xticks()
        ax.set_xticks(x_ticks)
        ax.set_xticklabels((x_ticks * 1e4).astype(int))
        ax.set_xticks(np.arange(0, 30e-4, 5e-4))
        ax.set_xticklabels(np.arange(0, 30e-4, 5e-4).astype(int))
        ax.set_xticks(x_ticks)
        ax.set_xticklabels((x_ticks * 1e4).astype(int))
        if "all" in bkg_name: plt.xlabel(r'BI$_{\rm tot}$ (ckky)')
        else: plt.xlabel(r'BI ($10^{-4}$ ckky)')
        plt.ylabel('Probability density')
        plt.legend(loc="upper right", frameon=False)
        pdf.savefig(bbox_inches='tight')
        plt.close()

        fig, ax = plt.subplots(figsize=(5,3.3))
        gms = json_file['refined_global_modes']
        c68 = json_file['ci_68']
        s_all = read_samples(f"../{exp}/mcmc_files/samples.h5", 'S')
        plt.hist(s_all, bins=nbins_S, histtype='stepfilled', density=True, color='dodgerblue', alpha=0.15, zorder=-12)
        plt.hist(s_all, bins=nbins_S, histtype='step', density=True, color='dodgerblue', lw=1, zorder=-12)
        plt.axvline(gms['S'], color='navy', linestyle="-", label="Global mode")
        counts, bin_edges = np.histogram(s_all, bins=nbins_S, density=True)
        max_bin_index = np.argmax(counts)
        max_bin_count = counts[max_bin_index]
        max_bin_start = bin_edges[max_bin_index]
        max_bin_end = bin_edges[max_bin_index + 1]
        plt.axvline((max_bin_end+max_bin_start)/2, color='r', linestyle="--", label="Marginalized mode")
        plt.axvline(c68['S'][0]['left'], color='r', linestyle=":", label="Smallest 68% CI")
        plt.axvline(c68['S'][0]['right'], color='r', linestyle=":")
        plt.xlim(0,50)
        plt.xlabel(r'S ($10^{-27}$ yr$^{-1}$)')
        plt.ylabel('Probability density')
        plt.legend(loc="upper right", frameon=False)
        pdf.savefig(bbox_inches='tight')
        plt.close()


def comparison_corr_uncorr_bkg(path):
    
    with PdfPages(f"comparison_of_different_BI_models.pdf") as pdf:
        
        #################### 3 UNCORR  BI ####################
        exp = "fit_10_l200_uniform_3BI_CorrEff" # change me!!
        fig, ax = plt.subplots(figsize=(5,3.3))
        json_file = json.load(open(f"{path}/ZeroNuFit.jl/output/{exp}/mcmc_files/fit_results.json"))
        gms = json_file[gm_type]
        bi_bege = read_samples(f"{path}/ZeroNuFit.jl/output/{exp}/mcmc_files/samples.h5", 'B_l200a_BG')
        plt.hist(bi_bege, bins=nbins, histtype='stepfilled', density=True, color=c_frenchblue, alpha=0.15, zorder=-10)
        plt.hist(bi_bege, bins=nbins, histtype='step', density=True, color=c_frenchblue, lw=1, zorder=-10)
        bi_icpc = read_samples(f"{path}/ZeroNuFit.jl/output/{exp}/mcmc_files/samples.h5", 'B_l200a_IC')
        plt.hist(bi_icpc, bins=nbins, histtype='stepfilled', density=True, color=c_dark_pistachio, alpha=0.15, zorder=-8)
        plt.hist(bi_icpc, bins=nbins, histtype='step', density=True, color=c_dark_pistachio, lw=1, zorder=-8)
        bi_ppc = read_samples(f"{path}/ZeroNuFit.jl/output/{exp}/mcmc_files/samples.h5", 'B_l200a_PC')
        plt.hist(bi_ppc, bins=nbins, histtype='stepfilled', density=True, color=c_lava, alpha=0.15, zorder=-6)
        plt.hist(bi_ppc, bins=nbins, histtype='step', density=True, color=c_lava, lw=1, zorder=-6)
        #################### TOTAL ####################
        exp = "fit_9_l200_uniform_1BI_CorrEff" # change me!!
        json_file = json.load(open(f"{path}/ZeroNuFit.jl/output/{exp}/mcmc_files/fit_results.json"))
        gms = json_file[gm_type]
        bi_bege = read_samples(f"{path}/ZeroNuFit.jl/output/{exp}/mcmc_files/samples.h5", f'{bkg_name}')
        plt.hist(bi_bege, bins=nbins, histtype='stepfilled', density=True, color='darkorange', alpha=0.15, zorder=-12)
        plt.hist(bi_bege, bins=nbins, histtype='step', density=True, color='darkorange', lw=1, zorder=-12)
        icpc_patch = mpatches.Patch(color=c_dark_pistachio, label='Uncorr. BI, ICPC (5 evts)', alpha=0.2)
        #################### style ####################
        bege_patch = mpatches.Patch(color=c_frenchblue, label='Uncorr. BI, BEGe (1 evts)', alpha=0.2)
        ppc_patch = mpatches.Patch(color=c_lava, label='Uncorr. BI, PPC (1 evts)', alpha=0.2)
        all_patch = mpatches.Patch(color='darkorange', label='Single BI (7 evts)', alpha=0.2)
        ax.legend(handles=[all_patch, bege_patch, icpc_patch, ppc_patch], loc='upper right', ncol=1, frameon=False)
        plt.xlim(0,50e-4)
        x_ticks = ax.get_xticks()
        ax.set_xticks(x_ticks)
        ax.set_xticklabels((x_ticks * 1e4).astype(int))
        ax.set_xticks(np.arange(0, 50e-4, 10e-4))
        ax.set_xticklabels(np.arange(0, 50e-4, 10e-4).astype(int))
        ax.set_xticks(x_ticks)
        ax.set_xticklabels((x_ticks * 1e4).astype(int))
        plt.xlabel(r'BI ($10^{-4}$ ckky)')
        plt.ylabel('Probability density')
        pdf.savefig(bbox_inches='tight')
        plt.close()


        #################### 3 CORR BI ####################
        exp = "fit_14_l200_uniform_3BIhier_CorrEff" # change me!!
        fig, ax = plt.subplots(figsize=(5,3.3))
        json_file = json.load(open(f"{path}/ZeroNuFit.jl/output/{exp}/mcmc_files/fit_results.json"))
        gms = json_file[gm_type]
        bi_bege = read_samples(f"{path}/ZeroNuFit.jl/output/{exp}/mcmc_files/samples.h5", 'B_l200a_BG')
        plt.hist(bi_bege, bins=nbins, histtype='stepfilled', density=True, color=c_frenchblue, alpha=0.15, zorder=-10)
        plt.hist(bi_bege, bins=nbins, histtype='step', density=True, color=c_frenchblue, lw=1, zorder=-10)

        bi_icpc = read_samples(f"{path}/ZeroNuFit.jl/output/{exp}/mcmc_files/samples.h5", 'B_l200a_IC')
        plt.hist(bi_icpc, bins=nbins, histtype='stepfilled', density=True, color=c_dark_pistachio, alpha=0.15, zorder=-8)
        plt.hist(bi_icpc, bins=nbins, histtype='step', density=True, color=c_dark_pistachio, lw=1, zorder=-8)

        bi_ppc = read_samples(f"{path}/ZeroNuFit.jl/output/{exp}/mcmc_files/samples.h5", 'B_l200a_PC')
        plt.hist(bi_ppc, bins=nbins, histtype='stepfilled', density=True, color=c_lava, alpha=0.15, zorder=-6)
        plt.hist(bi_ppc, bins=nbins, histtype='step', density=True, color=c_lava, lw=1, zorder=-6)
        #################### TOTAL ####################
        exp = "fit_9_l200_uniform_1BI_CorrEff" # change me!!
        json_file = json.load(open(f"{path}/ZeroNuFit.jl/output/{exp}/mcmc_files/fit_results.json"))
        gms = json_file[gm_type]
        bi_bege = read_samples(f"{path}/ZeroNuFit.jl/output/{exp}/mcmc_files/samples.h5", f'{bkg_name}')
        plt.hist(bi_bege, bins=nbins, histtype='stepfilled', density=True, color='darkorange', alpha=0.15, zorder=-12)
        plt.hist(bi_bege, bins=nbins, histtype='step', density=True, color='darkorange', lw=1, zorder=-12)
        icpc_patch = mpatches.Patch(color=c_dark_pistachio, label='Corr. BI, ICPC (5 evts)', alpha=0.2)
        bege_patch = mpatches.Patch(color=c_frenchblue, label='Corr. BI, BEGe (1 evts)', alpha=0.2)
        ppc_patch = mpatches.Patch(color=c_lava, label='Corr. BI, PPC (1 evts)', alpha=0.2)
        all_patch = mpatches.Patch(color='darkorange', label='Single BI (7 evts)', alpha=0.2)
        ax.legend(handles=[all_patch, bege_patch, icpc_patch, ppc_patch], loc='upper right', ncol=1, frameon=False)
        #################### style ####################
        plt.xlim(0,50e-4)
        x_ticks = ax.get_xticks()
        ax.set_xticks(x_ticks)
        ax.set_xticklabels((x_ticks * 1e4).astype(int))
        ax.set_xticks(np.arange(0, 50e-4, 10e-4))
        ax.set_xticklabels(np.arange(0, 50e-4, 10e-4).astype(int))
        ax.set_xticks(x_ticks)
        ax.set_xticklabels((x_ticks * 1e4).astype(int))
        plt.xlabel(r'BI ($10^{-4}$ ckky)')
        plt.ylabel(f'Counts / ({bin_width*1e4}' + r'$\times 10^{-4}$ ckky)')
        plt.ylabel('Probability density')
        pdf.savefig(bbox_inches='tight')
        plt.close()
        
from matplotlib.ticker import FuncFormatter
def comparison_newexpo_bkg(path):
    
    with PdfPages(f"comparison_for_new_expo_05122024.pdf") as pdf:
        fig, ax = plt.subplots(figsize=(5,3.3))
        exp = "fit_9_l200_3BI_new_data_05122024" #"fit_9_l200_3BI_new_data"
        json_file = json.load(open(f"{path}/output/{exp}/mcmc_files/fit_results.json"))
        gms = json_file[gm_type]
        bi_all_nu24_2 = read_samples(f"{path}/output/{exp}/mcmc_files/samples.h5", 'B_l200a_all')
        bi_batch5 = read_samples(f"{path}/output/{exp}/mcmc_files/samples.h5", 'B_l200a_V05')
        bi_coax = read_samples(f"{path}/output/{exp}/mcmc_files/samples.h5", 'B_l200a_SC')
        bi_ring = read_samples(f"{path}/output/{exp}/mcmc_files/samples.h5", 'B_l200a_ring')
        plt.hist(bi_all_nu24_2, bins=nbins, histtype='step', density=True, color='forestgreen', zorder=-10)
        plt.hist(bi_batch5, bins=nbins, histtype='step', density=True, color=c_frenchblue, zorder=-10)
        plt.hist(bi_coax, bins=nbins, histtype='step', density=True, color=c_lava, zorder=-10)
        plt.hist(bi_ring, bins=nbins, histtype='step', density=True, color='darkorange', zorder=-10)
        plt.plot([],[], color=c_frenchblue, label='Blinded Batch 05 (0 evts)')
        plt.plot([],[], color=c_lava, label='Blinded COAX (3 evts)')
        plt.plot([],[], color='darkorange', label='Blinded Ringing (1 evt)')
        plt.plot([],[], color='forestgreen', label='Blinded Nu24 detectors (6 evts)')
        ax.legend(loc='upper right', ncol=1, frameon=False)
        plt.xlim(0,0.01)
        def scale_x_labels(x, pos):
            return f'{x * 1e3:.0f}' 
        ax.xaxis.set_major_formatter(FuncFormatter(scale_x_labels))
        plt.xlabel(r'BI ($10^{-3}$ ckky)')
        plt.ylabel('Probability density')
        pdf.savefig(bbox_inches='tight')
        plt.close()

bin_width_S = 0.2
xmin_S = np.min(0)
xmax_S = np.max(110)
nbins_S = np.arange(xmin_S, xmax_S + bin_width_S, bin_width_S)
bin_width_S = 0.1
xmin_S = np.min(0)
xmax_S = np.max(110)
nbins_S_sqrt = np.arange(xmin_S, xmax_S + bin_width_S, bin_width_S)
def signal_priors():
    exp = "output_correct_inputs_uncorrelated_nuisances/v3/fit_9_l200_uniform_1BI_CorrEff" # change me!!
    exp_sqrt = "output_correct_inputs_uncorrelated_nuisances/v3/fit_16_l200_sqrt_1BI_CorrEff" # change me!!
    
    with PdfPages(f"signal_priors_comparison.pdf") as pdf:
      
        fig, ax = plt.subplots(figsize=(5,3.3))
        s_all = read_samples(f"../{exp}/mcmc_files/samples.h5", 'S')
        plt.hist(s_all, bins=nbins_S, histtype='stepfilled', density=True, color='dodgerblue', alpha=0.15, zorder=-12,label="Uniform prior")
        plt.hist(s_all, bins=nbins_S, histtype='step', density=True, color='dodgerblue', lw=1, zorder=-12)
        s_all = read_samples(f"../{exp_sqrt}/mcmc_files/samples.h5", 'S')
        plt.hist(s_all, bins=nbins_S_sqrt, histtype='stepfilled', density=True, color=c_dark_pistachio, alpha=0.15, zorder=-12,label=r"$1/\sqrt{\Gamma_{\rm 1/2}}$ prior")
        plt.hist(s_all, bins=nbins_S_sqrt, histtype='step', density=True, color=c_dark_pistachio, lw=1, zorder=-12)
        plt.xlim(0,20)
        plt.xlabel(r'S ($10^{-27}$ yr$^{-1}$)')
        plt.ylabel('Probability density')
        plt.legend(loc="upper right", frameon=False)
        pdf.savefig(bbox_inches='tight')
        plt.close()




comparison_newexpo_bkg("../")
different_bkg_models()
sig_and_bkg()
different_bkg_models()
signal_priors()