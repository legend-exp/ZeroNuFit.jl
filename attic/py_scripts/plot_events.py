"""
Script to plot energy events for different experiments given entry JSON files.
Main Authors: Sofia Calgaro
"""
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
    

bin_width = 1
xmin = np.min(1900)
xmax = np.max(3000)
nbins = np.arange(xmin, xmax + bin_width, bin_width)
gm_type = "marginalized_modes"
center = 1930
expo = 42.25
fit_range = [[1930, 2099], [2109, 2114], [2124, 2190]]
Qbb = 2039.06


with PdfPages(f"distribution_of_events.pdf") as pdf:
    gerdaI_events_json = json.load(open("../legend-0vbb-config/gerda/events_gerda_pI.json"))
    gerdaII_events_json = json.load(open("../legend-0vbb-config/gerda/events_gerda.json"))
    l200_events_json = json.load(open("../legend-0vbb-config/legend/events_l200_all_05022025.json"))
    mjd_events_json = json.load(open("../legend-0vbb-config/mjd/events_mjd_og.json"))
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

    json_dict = {
        "ph1_golden": gerdaI_events_json, 
        "ph1_bege": gerdaI_events_json, 
        "ph1_silver": gerdaI_events_json, 
        "ph1_extra": gerdaI_events_json, 
        "B_gerda_all_pII": gerdaII_events_json, 
        "B_l200a_Nu24": l200_events_json, 
        "B_l200a_extra": l200_events_json, 
        "mjd-DS0": mjd_events_json, 
        "mjd-mod1": mjd_events_json, 
        "mjd-mod2": mjd_events_json,
    }

    title_names = {
        "ph1_golden": "GERDA PhI golden", 
        "ph1_bege": "GERDA PhI BEGe", 
        "ph1_silver": "GERDA PhI silver", 
        "ph1_extra": "GERDA PhI extra", 
        "B_gerda_all_pII": "GERDA PhII", 
        "B_l200a_Nu24": "LEGEND Nu24", 
        "B_l200a_extra": "LEGEND extra", 
        "mjd-DS0": "MJD DS0", 
        "mjd-mod1": "MJD module 1", 
        "mjd-mod2": "MJD module 2",
    }
        
    for idx,bkg_name in enumerate(["ph1_golden", "ph1_bege", "ph1_silver", "ph1_extra", "B_gerda_all_pII", "B_l200a_Nu24", "B_l200a_extra", "mjd-DS0", "mjd-mod1", "mjd-mod2"]):
        fig, ax = plt.subplots(figsize=(5,2))
        evt_dict = json_dict[bkg_name]

        events = []
        for entry in evt_dict["events"]:
            if bkg_name=="ph1_golden" and entry["experiment"]=="GERDA_pI" and entry["detector"]=="coax" and (entry["timestamp"] > 1320849782-1 and entry["timestamp"] <1369143782+1): events.append(entry["energy"])
            if bkg_name=="ph1_bege" and entry["experiment"]=="GERDA_pI" and entry["detector"]=="bege" and (entry["timestamp"] > 1320849782-1 and entry["timestamp"] <1369143782+1): events.append(entry["energy"])
            if bkg_name=="ph1_silver" and entry["experiment"]=="GERDA_pI" and (entry["timestamp"] > 339945251-1 and entry["timestamp"] <342589688+1): events.append(entry["energy"])
            if bkg_name=="ph1_extra" and entry["experiment"]=="GERDA_pI" and (entry["timestamp"] > 1370007782-1 and entry["timestamp"] <1380548582+1): events.append(entry["energy"])
                
            if bkg_name=="B_gerda_all_pII" and entry["experiment"]=="GERDA": events.append(entry["energy"])
                
            if bkg_name=="B_l200a_Nu24" and entry["experiment"]=="L200" and entry["detector"][0]!="C" and entry["detector"][0:2]!="V05": events.append(entry["energy"])
            if bkg_name=="B_l200a_extra" and entry["experiment"]=="L200" and (entry["detector"][0]=="C" or entry["detector"][0:2]=="V05"): events.append(entry["energy"])
                
            if bkg_name=="mjd-DS0" and entry["experiment"]=="MJD" and "DS0" in entry["detector"]: events.append(entry["energy"])
            if bkg_name=="mjd-mod1" and entry["experiment"]=="MJD":
                if entry["detector"] in ["DS1","DS2", "DS3", "DS5a", "DS5b", "DS5c", "DS6a", "DS6b", "DS6c","DS7", "DS8P"] and entry["timestamp"] < 140: events.append(entry["energy"])
            if bkg_name=="mjd-mod2" and entry["experiment"]=="MJD":
                if entry["detector"] in ["DS4", "DS5a", "DS5b", "DS5c", "DS6a", "DS6b", "DS6c", "DS8P", "DS8I"] and (entry["timestamp"] >= 140 and entry["timestamp"] <280): events.append(entry["energy"])
                
        plt.hist(events, bins=nbins, histtype='stepfilled', density=False, label=F"{title_names[bkg_name]}: {len(events)} events", color=custom_colors[idx])
        plt.xlabel("Energy (keV)")    
        plt.ylabel('Counts / keV')
        plt.axvspan(2099, 2109, facecolor="black", alpha=0.4)
        plt.axvspan(2114, 2124, facecolor="black", alpha=0.4)
        if "mjd" in bkg_name:
            plt.xlim(1950,2350)
            plt.axvspan(2199.1, 2209.1, facecolor="black", alpha=0.4)
            plt.xticks(np.linspace(1950,2350,4))
        else: 
            plt.xlim(1930,2190)
            plt.xticks(np.linspace(1930,2190,4))
        plt.ylim(0,4.1)
        plt.yticks([0,1,2,3,4])
        ax.legend(loc=(0,1), frameon=False, ncol=3)
        pdf.savefig(bbox_inches='tight')
        plt.close()



