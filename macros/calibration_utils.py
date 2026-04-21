import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import os

plt.rcParams.update({'font.size': 24})

def plot_LO_calibration(csv_file, outdir="/eos/home-c/cgiraldi/www/MTD/MTDTB_CERN_Sep25/energy_intercalibration/LO_calibrations/"):
    df = pd.read_csv(csv_file)
    df_L = df[df['side'] == 'L'].sort_values('bar')
    df_R = df[df['side'] == 'R'].sort_values('bar')
    rms_L = df_L['calib'].std()
    rms_R = df_R['calib'].std()    
    plt.figure(figsize=(12,10))
    #plt.plot(df_L['bar'], df_L['calib'], marker='o', linestyle='-', label=f'L (RMS={rms_L:.3f})', color="red")
    #plt.plot(df_R['bar'], df_R['calib'], marker='s', linestyle='-', label=f'R (RMS={rms_R:.3f})', color="blue")
    plt.plot(df_L['bar'].values, df_L['calib'].values, marker='o', linestyle='-', label=f'L (RMS={rms_L:.3f})', color="red")
    plt.plot(df_R['bar'].values, df_R['calib'].values, marker='s', linestyle='-', label=f'R (RMS={rms_R:.3f})', color="blue")
    plt.axhline(1.0, linestyle=':', linewidth=1)
    plt.ylim(0.65,1.35)
    plt.xlabel("Bar Index")
    plt.ylabel("LO calibration coefficients")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plot_name =	os.path.splitext(os.path.basename(csv_file))[0] + ".png"
    output = os.path.join(outdir, plot_name)
    plt.savefig(output, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved LO calibrations plot in: {output}")

def plot_TOFHIR_calibration(csv_file):
    calib_min = 0.5
    calib_max = 1.5
    csv_path = Path(csv_file)
    outdir = csv_path.parent
    df = pd.read_csv(csv_file)
    sides = ["L", "R"]
    colors = {"L": "red", "R": "blue"}
    # we expect one vov only per moduleCharacterization now, but leaving the possibility of having more than one value
    vovs = sorted(df["vov"].unique())
    thresholds = sorted(df["th"].unique())
    for vov in vovs:
        # calibration vs bar
        df_vov = df[df["vov"] == vov]
        for thr in thresholds:
            subset = df_vov[df_vov["th"] == thr]
            plt.figure(figsize=(12,10))
            for side in sides:
                s = subset[subset["side"] == side].sort_values("bar")
                rms = s["calib"].std()
                label_plot = f"{side} (RMS={rms:.3f})"
                plt.plot(s["bar"].values, s["calib"].values, marker="o", label=label_plot, color=colors[side])
            plt.xlabel("Bar")
            plt.ylabel("TOFHIR calibration")
            plt.title(f"Vov={vov}, th={thr}")
            plt.ylim(calib_min,calib_max)
            plt.legend()
            plt.grid(True)
            filename = outdir / f"TOFHIR_calibration_vs_bar_vov{vov}_th{thr}.png"
            plt.savefig(filename, dpi=300, bbox_inches="tight")
            plt.close()
        # heatmap
        for side in sides:
            subset = df_vov[df_vov["side"] == side]
            heatmap_data = subset.pivot(index="th", columns="bar", values="calib")
            plt.figure(figsize=(12,10))
            plt.imshow(heatmap_data, aspect="auto", vmin=calib_min, vmax=calib_max, cmap="turbo")
            plt.colorbar(label="Calibration")
            plt.xticks(range(len(heatmap_data.columns)), heatmap_data.columns)
            plt.yticks(range(len(heatmap_data.index)), heatmap_data.index)
            plt.title(f"{side} - Vov={vov}")
            plt.xlabel("Bar")
            plt.ylabel("Threshold")
            filename = outdir / f"TOFHIR_calibration_heatmap_{side}_vov{vov}.png"
            plt.savefig(filename, dpi=300, bbox_inches="tight")
            plt.close()
    print(f"Saved TOFHIR calibrations plot in: {outdir}")

def plot_TOFHIR_calibration_multi(csv_files):
    dfs = []
    for f in csv_files:
        dfs.append(pd.read_csv(f))
    df = pd.concat(dfs, ignore_index=True)
    outdir = Path(csv_files[0]).parent
    sides = ["L","R"]
    thresholds = sorted(df["th"].unique())
    vovs = sorted(df["vov"].unique())
    colors = ["green","orange","purple","dodgerblue","magenta","cyan"]
    for thr in thresholds:
        for side in sides:
            subset = df[(df["th"] == thr) & (df["side"] == side)]
            plt.figure(figsize=(12,10))
            for i,vov in enumerate(vovs):
                s = subset[subset["vov"] == vov].sort_values("bar")
                rms = s["calib"].std()
                label = f"Vov={vov} (RMS={rms:.3f})"
                plt.plot(s["bar"].values, s["calib"].values, marker="o", color=colors[i % len(colors)], label=label)
            plt.xlabel("Bar", fontsize=18)
            plt.ylabel("TOFHIR calibration", fontsize=18)
            plt.ylim(0.5,1.5)
            plt.legend(fontsize=14)
            plt.grid(True)
            filename = outdir / f"TOFHIR_calibration_vs_bar_side{side}_th{thr}.png"
            plt.savefig(filename, dpi=300, bbox_inches="tight")
            plt.close()
    print(f"Saved TOFHIR multi-Vov plots in: {outdir}")


def get_rebin_factor(integral):
    # if integral > 15000:
    #     return 1
    # elif integral > 5000:
    #     return 2
    # elif integral > 500:
    #     return 4
    # elif integral > 100:
    #     return 8
    # elif integral > 50:
    #     return 16
    # else:
    #     return 32
    if integral > 15000:
        return 1
    elif integral > 10000:
        return 2
    elif integral > 6000:
        return 4
    elif integral > 4000:
        return 8
    elif integral > 2000:
        return 16
    elif integral > 1000:
        return 32
    else:
        return 64
    
def read_min_energy(txt_file):
    min_energy = {}
    with open(txt_file, "r") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) != 3:
                continue
            bar = int(parts[0])
            vov = float(parts[1])
            energy = float(parts[2])
            min_energy[(bar, vov)] = energy
    return min_energy
