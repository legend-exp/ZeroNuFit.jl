### plot_toys_results.py
#
# Authors: Sofia Calgaro, Toby Dixon
# 
###
import sys
import json
import argparse
import os
import csv
import base64
import numpy as np
import pandas as pd
import pickle
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

## some colours
c_pistachio = (0.58, 0.87, 0.45)
c_columbiablue = (0.61, 0.87, 1.0)
c_dark_columbiablue = (0.61, 0.87*0.85, 1)
c_frenchblue = (0.0, 0.45, 0.73)
c_darl_gray = (0.33, 0.33, 0.33)
c_dark_pistachio = (0.58*0.85, 0.87*0.85, 0.45*0.85)
c_lava = (0.94, 0.01, 0.05)

def main():
    parser = argparse.ArgumentParser(
        description="Script to load fit results"
    )
    parser.add_argument(
        "--fit_fake",
        type=str,
        help="Name of the fit output folder for the NO SIGNAL (=bkg only) fit case",
    )
    parser.add_argument(
        "--fit_real",
        type=str,
        default=None,
        help="Name of the fit output folder for the SIGNAL+BKG fit case",
    )
    parser.add_argument(
        "--path", "--p",
        type=str,
        default="output",
        help="Output path: 'output' is default",
    )
    parser.add_argument(
        "--global_mode", "--gm",
        default="refined_global_modes",
        type=str,
        help="Global mode type: 'refined_global_modes' (default) or 'global_modes'",
    )
    parser.add_argument(
        "--sensitivity", "--s",
        default="sensitivity",
        type=str,
        help="Name of the folder containing sensitivity study outputs: 'sensitivity' (default)",
    )
    parser.add_argument(
        "--width", "--w",
        default="fwhm",
        type=str,
        help="Name of the width paramter: 'fwhm' (default) or 'width'",
    )
    parser.add_argument(
        "--qbb",
        default="False",
        type=str,
        help="'False' if you do not want to break the full distribution into smaller distributions based on the number of events that were generated close to Qbb",
    )
    parser.add_argument(
        "--bkg_name", "--bkg",
        default="B_l200a_all",
        type=str,
        help="Name of the BI in the fit: 'B_l200_all' (default)",
    )
    
    args = parser.parse_args()
    gm_type = args.global_mode
    output = args.path
    fit_name = args.fit_fake
    real_fit = args.fit_real
    sensitivity_folder = args.sensitivity
    width = args.width
    qbb = True if args.qbb == "False" else True
    bkg_name = args.bkg_name
    
    # fit without signal
    json_folder = os.path.join(output, fit_name, sensitivity_folder, "mcmc_files")
    json_files = os.listdir(json_folder)
    json_files = [os.path.join(json_folder, f) for f in json_files if ".json" in f]
    fake_data = [f.replace("mcmc_files", "fake_data").replace("fit_results_", "fake_data") for f in json_files]
    partitions_list = json.load(open(json_files[0]))['config']['partitions'] #it's all the same

    # fit with signal
    if real_fit is not None:
        real_fit_json = json.load(open(os.path.join(output, real_fit, "mcmc_files/fit_results.json"))) 
        s90_real = real_fit_json['quantile90']['S']
        sgm_real = real_fit_json['refined_global_modes']['S']
    else: 
        s90_real = 0
        sgm_real = 0
        
    below_Sreal = []
    tot_S = []

    bgm_values = []
    s90_values = []
    sgm_values = []
    bkg_events = []
    s90_values_0evt = []
    s90_values_1evt = []
    s90_values_2evt = []
    s90_values_other = []
    sgm_values_0evt = []
    sgm_values_1evt = []
    sgm_values_2evt = []
    sgm_values_other = []

    qbb = 2039.06
    gerda_two_phases = False
    s90_values_phII = []
    s90_values_0evt_phI = []
    s90_values_1evt_phI = []
    s90_values_2evt_phI = []
    s90_values_other_phI = []
    s90_values_0evt_phII = []
    s90_values_1evt_phII = []
    s90_values_2evt_phII = []
    s90_values_other_phII = []

    for idx,file in enumerate(json_files):

        json_file = json.load(open(file))
        fake_file = json.load(open(fake_data[idx]))
        energies = []
        ct = 0
        phI = False
        tot_no_fake_data = len(fake_file["events"])
        for entry in fake_file["events"]: # this is a list of # of events
            en = entry["energy"]
            timestamp = entry["timestamp"]
            detector = entry["detector"]
            experiment = entry["experiment"]
            fwhm = 0
            # find the fwhm of this fake event
            for partition_file in partitions_list:
                partition_file = json.load(open(os.path.join("../", partition_file)))
                parts = partition_file["partitions"] # list of partitions
                for p in parts: # list of events for a given partition 
                    p = partition_file["partitions"][p]
                    for evt in p: # 1 event
                        start_ts = evt["start_ts"]
                        end_ts = evt["end_ts"]
                        det = evt["detector"]
                        if timestamp >= start_ts-1 and timestamp<=end_ts+1 and detector == det:
                            fwhm = evt[width] if width=="fwhm" else evt[width]*2.355
                            #if experiment != "MJD" or (experiment=="MJD" and "partitions_mjd_new.json" in partition_file):
                            #    fwhm = evt[width] if width=="fwhm" else evt[width]*2.355
                            #else:
                            #    fwhm = evt[width] if width=="fwhm" else evt[width]*2.355*evt["sigma"]
                            break        

            roi_low = qbb - fwhm
            roi_upp = qbb + fwhm
            if en >= roi_low and en <= roi_upp:
                ct += 1
                # find the experiment to which they belong to
                if gerda_two_phases is True and phI is False:
                    # for sure this is not phI
                    if detector not in ["coax", "bege"]:
                        continue
                    else:
                        phI = True

        quantiles = json_file['quantile90']
        s90 = quantiles['S']
        gms = json_file['refined_global_modes']
        sgm = gms['S']

        tot_S.append(s90)
        if s90<=s90_real and real_fit is not None:
            below_Sreal.append(s90)

        if phI == True: print(en, phI, idx+1)
        s90_values.append(s90)
        sgm_values.append(sgm)
        bgm_values.append(gms[bkg_name]*1e4)
        if gerda_two_phases == True:
            if phI == False:
                s90_values_phII.append(s90)

        bkg_events.append(tot_no_fake_data-ct)        
        
        if ct == 0:
            s90_values_0evt.append(s90)
            sgm_values_0evt.append(sgm)
            if gerda_two_phases == True:
                if phI == False:
                    s90_values_0evt_phII.append(s90)

        if ct == 1:
            s90_values_1evt.append(s90)
            sgm_values_1evt.append(sgm)
            if gerda_two_phases == True:
                if phI == False:
                    s90_values_1evt_phII.append(s90)

        if ct == 2:
            s90_values_2evt.append(s90)
            sgm_values_2evt.append(sgm)
            if gerda_two_phases == True:
                if phI == False:
                    s90_values_2evt_phII.append(s90)

        if ct > 2:
            s90_values_other.append(s90)
            sgm_values_other.append(sgm)
            if gerda_two_phases == True:
                if phI == False:
                    s90_values_other_phII.append(s90)

    print("Number of generated toys:", len(tot_S))
    if qbb is True:
        print("..with 0 cts only in Qbb+-FWHM:", len(s90_values_0evt)/len(tot_S)*100)
        print("..with 1 cts only in Qbb+-FWHM:", len(s90_values_1evt)/len(tot_S)*100)
        print("..with 2 cts only in Qbb+-FWHM:", len(s90_values_2evt)/len(tot_S)*100)
        print("..with >2 cts in Qbb+-FWHM:", len(s90_values_other)/len(tot_S)*100)

    median = np.median(s90_values)
    lower_bound = np.percentile(s90_values, 16)
    upper_bound = np.percentile(s90_values, 84)
    print("*** S ***")
    print("Mean of toys:", np.mean(s90_values))
    print("Median of toys:", median)
    print("Lower bound of the smallest 68% CI:", lower_bound)
    print("Upper bound of the smallest 68% CI:", upper_bound)
    print("*** T ***")
    print("Mean of toys:", 10/np.mean(s90_values))
    print("Median of toys:", 10/median)
    print("Lower bound of the smallest 68% CI:", 10/upper_bound)
    print("Upper bound of the smallest 68% CI:", 10/lower_bound)
    """
    print("*** BI ***")
    print("Mean of toys:", np.mean(bgm_values))
    print("Median of toys:", np.median(bgm_values))
    print("Lower bound of the smallest 68% CI:", np.percentile(bgm_values, 16))
    print("Upper bound of the smallest 68% CI:", np.percentile(bgm_values, 84))
    """

    if real_fit is not None:
        print("Number of toys below S_observed:", len(below_Sreal))
        print("Percentage over the total number of toys:", len(below_Sreal) / len(tot_S), "%") 


    # you can change me (it's for plotting only)
    bin_width = 0.25
    bin_width_gm = 0.25
    min_s90 = np.min(s90_values)-1
    max_s90 = np.max(s90_values)+1
    nbins = np.arange(min_s90, max_s90 + bin_width, bin_width)
    min_sgm = np.min(sgm_values)-1
    max_sgm = np.max(sgm_values)+1
    nbins_gm = np.arange(min_sgm, max_sgm + bin_width_gm, bin_width_gm)

    # let's plot
    with PdfPages(f"{fit_name}_distr.pdf") as pdf:
        
        #### SIGNAL

        ######### log scale ############
        fig, ax = plt.subplots(figsize=(8,5))
        ax.hist(s90_values, bins=nbins, histtype='stepfilled', color='darkblue', alpha=0.1)
        ax.hist(s90_values, histtype="step", bins=nbins,color='darkblue', label="All toys")
        ax.set_ylabel(f"Counts / {bin_width}")
        ax.set_xlabel(r"Signal upper limit 90% CI ($10^{-27}$ yr$^{-1}$)")
        ax.set_yscale("log")
        if real_fit is not None: ax.axvline(s90_real, linestyle='--', color='darkorange', label="Observed limit")
        ax.axvline(np.mean(s90_values), linestyle='--', color='saddlebrown', label="Mean of the distribution")
        ax.axvline(median, linestyle='--', color='red', label="Median of the distribution")
        ax.set_xlim(min_s90, 40)#max_s90)
        ax.legend(loc='upper right', frameon=False)
        pdf.savefig(bbox_inches='tight')
        plt.close()

        fig, ax = plt.subplots(figsize=(8,5))
        ax.hist(s90_values_0evt, bins=nbins, histtype='stepfilled', color='darkblue', alpha=0.1)
        ax.hist(s90_values_0evt, bins=nbins, histtype='step', color='darkblue', lw=1, label=r"Toys w. 0 cts in $Q_{\beta\beta} \pm$FWHM")
        ax.hist(s90_values_1evt, bins=nbins, histtype='stepfilled', color=c_frenchblue, alpha=0.1)
        ax.hist(s90_values_1evt, bins=nbins, histtype='step', color=c_frenchblue, lw=1, label=r"Toys w. 1 cts in $Q_{\beta\beta} \pm$FWHM")
        ax.hist(s90_values_2evt, bins=nbins, histtype='stepfilled', color=c_dark_pistachio, alpha=0.1)
        ax.hist(s90_values_2evt, bins=nbins, histtype='step', color=c_dark_pistachio, lw=1, label=r"Toys w. 2 cts in $Q_{\beta\beta} \pm$FWHM")
        ax.hist(s90_values_other, bins=nbins, histtype='stepfilled', color='gray', alpha=0.3)
        ax.hist(s90_values_other, bins=nbins, histtype='step', color='gray', lw=1, label=r"Toys w. >2 cts in $Q_{\beta\beta} \pm$FWHM")
        if real_fit is not None: ax.axvline(s90_real, linestyle='--', color='darkorange', label="Observed limit")
        ax.axvline(np.mean(s90_values), linestyle='--', color='saddlebrown', label="Mean of the distribution")
        ax.axvline(median, linestyle='--', color='red', label="Median of the distribution")
        ax.set_ylabel(f"Counts / {bin_width}")
        ax.set_xlabel(r"Signal upper limit 90% CI ($10^{-27}$ yr$^{-1}$)")
        ax.set_yscale("log")
        ax.set_xlim(min_s90, 40)#max_s90)
        ax.legend(loc='upper right', frameon=False)
        pdf.savefig(bbox_inches='tight')
        plt.close()

        ######### linear scale ############
        fig, ax = plt.subplots(figsize=(8,5))
        ax.hist(s90_values, bins=nbins, histtype='stepfilled', color='darkblue', alpha=0.1)
        ax.hist(s90_values, histtype="step", bins=nbins,color='darkblue', label="All toys")
        ax.set_ylabel(f"Counts / {bin_width}")
        ax.set_xlabel(r"Signal upper limit 90% CI ($10^{-27}$ yr$^{-1}$)")
        ax.set_xlim(min_s90, 40)#max_s90)
        if real_fit is not None: ax.axvline(s90_real, linestyle='--', color='darkorange', label="Observed limit")
        ax.axvline(np.mean(s90_values), linestyle='--', color='saddlebrown', label="Mean of the distribution")
        ax.axvline(median, linestyle='--', color='red', label="Median of the distribution")
        ax.legend(loc='upper right', frameon=False)
        pdf.savefig(bbox_inches='tight')
        plt.close()

        fig, ax = plt.subplots(figsize=(8,5))
        ax.hist(s90_values_0evt, bins=nbins, histtype='stepfilled', color='darkblue', alpha=0.1)
        ax.hist(s90_values_0evt, bins=nbins, histtype='step', color='darkblue', lw=1, label=r"Toys w. 0 cts in $Q_{\beta\beta} \pm$FWHM")
        ax.hist(s90_values_1evt, bins=nbins, histtype='stepfilled', color=c_frenchblue, alpha=0.1)
        ax.hist(s90_values_1evt, bins=nbins, histtype='step', color=c_frenchblue, lw=1, label=r"Toys w. 1 cts in $Q_{\beta\beta} \pm$FWHM")
        ax.hist(s90_values_2evt, bins=nbins, histtype='stepfilled', color=c_dark_pistachio, alpha=0.1)
        ax.hist(s90_values_2evt, bins=nbins, histtype='step', color=c_dark_pistachio, lw=1, label=r"Toys w. 2 cts in $Q_{\beta\beta} \pm$FWHM")
        ax.hist(s90_values_other, bins=nbins, histtype='stepfilled', color='gray', alpha=0.3)
        ax.hist(s90_values_other, bins=nbins, histtype='step', color='gray', lw=1, label=r"Toys w. >2 cts in $Q_{\beta\beta} \pm$FWHM")
        ax.set_ylabel(f"Counts / {bin_width}")
        ax.set_xlabel(r"Signal upper limit 90% CI ($10^{-27}$ yr$^{-1}$)")
        ax.set_xlim(min_s90, 40)#max_s90)
        if real_fit is not None: ax.axvline(s90_real, linestyle='--', color='darkorange', label="Observed limit")
        ax.axvline(np.mean(s90_values), linestyle='--', color='saddlebrown', label="Mean of the distribution")
        ax.axvline(median, linestyle='--', color='red', label="Median of the distribution")
        
        ax.legend(loc='upper right', frameon=False)
        pdf.savefig(bbox_inches='tight')
        plt.close()
        
        
        #### BACKGROUND 
        fig, ax = plt.subplots(figsize=(8,5))
        ax.scatter(bgm_values, bkg_events, color='darkblue')
        ax.set_ylabel(f"No. of fake bkg events (cts)")
        ax.set_xlabel(r"BI best fit ($10^{-4}$ ckky)")
        if real_fit is not None: ax.axvline(real_fit_json['refined_global_modes'][bkg_name]*1e4, linestyle='--', color='darkorange', label="Observed limit")
        ax.legend(loc='upper right', frameon=False)
        pdf.savefig(bbox_inches='tight')
        plt.close()
        
        min_bgm = 0
        max_bgm = np.max(bgm_values)+1
        nbins_gm = np.arange(min_bgm, max_bgm + 0.1, 0.1)
        fig, ax = plt.subplots(figsize=(8, 5))
        heatmap = ax.hist2d(bgm_values, bkg_events, bins=(np.linspace(-2.5, 31.5, 34), nbins_gm), cmap='viridis')
        cbar = fig.colorbar(heatmap[3], ax=ax)
        cbar.set_label('Number of counts')
        if real_fit:
            ax.axvline(
                real_fit_json['refined_global_modes'][bkg_name] * 1e4, 
                linestyle='--', color='darkorange', label="Observed limit"
            )
        ax.set_ylabel("No. of fake bkg events (cts)")
        ax.set_xlabel(r"BI best fit ($10^{-4}$ ckky)")
        ax.legend(loc='upper right', frameon=False)
        pdf.savefig(bbox_inches='tight')
        plt.close()  

        min_bgm = 0
        max_bgm = np.max(bgm_values)+1
        bin_width = 0.1
        nbins_gm = np.arange(min_bgm, max_bgm + 0.1, 0.1)
        fig, ax = plt.subplots(figsize=(8,5))
        ax.hist(bgm_values, bins=nbins_gm, histtype='stepfilled', color='darkblue', alpha=0.1)
        ax.hist(bgm_values, bins=nbins_gm, histtype="step", color='darkblue', label="All toys")
        ax.set_ylabel(f"Counts / 0.1")
        ax.set_xlabel(r"BI best fit ($10^{-4}$ ckky)")
        if real_fit is not None: ax.axvline(real_fit_json['refined_global_modes'][bkg_name]*1e4, linestyle='--', color='darkorange', label="Observed limit")
        ax.axvline(np.mean(bgm_values), linestyle='--', color='saddlebrown', label="Mean of the distribution")
        ax.axvline(np.median(bgm_values), linestyle='--', color='red', label="Median of the distribution")
        ax.legend(loc='upper right', frameon=False)
        pdf.savefig(bbox_inches='tight')
        plt.close()
        
        
        min_bgm = 0
        max_bgm = np.max(bgm_values)+1
        nbins_gm = np.arange(min_bgm, max_bgm + 0.5, 0.5)
        fig, ax = plt.subplots(figsize=(8,5))
        ax.hist(bgm_values, bins=nbins_gm, histtype='stepfilled', color='darkblue', alpha=0.1)
        ax.hist(bgm_values, bins=nbins_gm, histtype="step", color='darkblue', label="All toys")
        ax.set_ylabel(f"Counts / 0.5")
        ax.set_xlabel(r"BI best fit ($10^{-4}$ ckky)")
        if real_fit is not None: ax.axvline(real_fit_json['refined_global_modes'][bkg_name]*1e4, linestyle='--', color='darkorange', label="Observed limit")
        ax.axvline(np.mean(bgm_values), linestyle='--', color='saddlebrown', label="Mean of the distribution")
        ax.axvline(np.median(bgm_values), linestyle='--', color='red', label="Median of the distribution")
        ax.legend(loc='upper right', frameon=False)
        pdf.savefig(bbox_inches='tight')
        plt.close()       
        
        
        min_bgm = 0
        max_bgm = np.max(bgm_values)+1
        nbins_gm = np.arange(min_bgm, max_bgm + 1, 1)
        fig, ax = plt.subplots(figsize=(8,5))
        ax.hist(bgm_values, bins=nbins_gm, histtype='stepfilled', color='darkblue', alpha=0.1)
        ax.hist(bgm_values, bins=nbins_gm, histtype="step", color='darkblue', label="All toys")
        ax.set_ylabel(f"Counts / 1")
        ax.set_xlabel(r"BI best fit ($10^{-4}$ ckky)")
        if real_fit is not None: ax.axvline(real_fit_json['refined_global_modes'][bkg_name]*1e4, linestyle='--', color='darkorange', label="Observed limit")
        ax.axvline(np.mean(bgm_values), linestyle='--', color='saddlebrown', label="Mean of the distribution")
        ax.axvline(np.median(bgm_values), linestyle='--', color='red', label="Median of the distribution")
        ax.legend(loc='upper right', frameon=False)
        pdf.savefig(bbox_inches='tight')
        plt.close()      


if __name__ == "__main__":
    main()