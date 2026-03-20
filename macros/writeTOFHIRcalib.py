#!/usr/bin/env python3
from utils import *
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)

# -------------------------------------------------------
# YOUR PATHS
eos_path = "/eos/home-s/spalluot/MTD/TB_CERN_Sep25/Lab5015Analysis"
outdir_base = "/eos/home-s/spalluot/www/MTD/MTDTB_CERN_Sep25/energy_intercalibration/"
# ------------------------------------------------------- 

# ------- parser -------
parser = argparse.ArgumentParser(description="Energy spectra calibration")
parser.add_argument("-i", "--inputFile", required=True, type=str, help="Input ROOT file")
parser.add_argument("-o", "--outFolder", required=True, type=str, help="Output folder")
parser.add_argument("-sm", "--sensorModuleID", required=True, type=str, help="Sensor module ID for QAQC calibrations")
parser.add_argument("--minEnergy", required=True, type=str, help="minEnergy config file name")
parser.add_argument("--drawCalib", required=False, action="store_true", help="Produce calibration plots vs bar")
parser.add_argument("--writeOutFile", required=False, action="store_true", help="Write output file with corrected energy spectra")
parser.add_argument("--drawComparison", required=False, action="store_true", help="Draw comparison: raw vs LO vs LO+TOFHIR")
parser.add_argument("--fitCheck", required=False, action="store_true", help="Check energy spectrum fit")
parser.add_argument("--drawMPVvsBar", required=False, action="store_true", help="Draw MPV vs bar at different calibration stages")
args = parser.parse_args()
input_file = args.inputFile
sensor_module_id = args.sensorModuleID

# -- directories
inputdir = f"{eos_path}/plots"
outdir = f"{outdir_base}/{args.outFolder}"
os.makedirs(outdir, exist_ok=True)
if args.fitCheck:
    outdir_check = f"{outdir}/fit_check"
    os.makedirs(outdir_check, exist_ok=True)

# -- csv files
LO_csv = f"{inputdir}/module_{sensor_module_id}_LO_calibration_factors.csv"
TOFHIR_LO_csv = f"{outdir}/TOFHIR_LO_calibration_factors.csv"
TOFHIR_csv = f"{outdir}/TOFHIR_calibration_factors.csv"
# -- min energy values for the fit range 
minEnergy_txt = f"{eos_path}/cfg/minEnergies_{args.minEnergy}.txt"
min_energy_dict = read_min_energy(minEnergy_txt)

# ------------------------------------------------------- 
# STEP 1: get LO calibration from csv 
# -------------------------------------------------------
if not os.path.exists(LO_csv):
    print("LO csv file not found:", LO_csv)
    sys.exit(0)
LO_df = pd.read_csv(LO_csv)
calib_LO = LO_df.set_index(["bar","side"])["calib"].to_dict()
print(f"LO calibration factors extracted")

# -------------------------------------------------------
# STEP 2: fit energy calibrated spectra with a Landau fit
# ------------------------------------------------------- 
# -- read energy spectra from moduleChar_step1 file
# -- apply LO calibration factors 
h_energy_raw = {}
h_energy = {}
bars, sides, Vovs, thresholds = set(), set(), set(), set()
f = ROOT.TFile.Open(f"{inputdir}/moduleCharacterization_step1_{input_file}.root")
for key_obj in f.GetListOfKeys():
    name = key_obj.GetName()
    if not name.startswith("h1_energy_bar"):
        continue
    parts = name.split("_")
    bar = int(parts[2][3:5])
    side = parts[2][5:]  # L or R
    vov = float(parts[3][3:])
    thr = int(parts[4][2:])
    bars.add(bar)
    sides.add(side)
    Vovs.add(vov)
    thresholds.add(thr)
    h = f.Get(name).Clone()
    h.SetDirectory(0)
    h_energy_raw[(bar, side, vov, thr)] = h.Clone()
    h_energy_raw[(bar, side, vov, thr)].SetDirectory(0)
    h_calib = h.Clone()
    if "L-R" != side:
        h_calib.Reset()
        h_calib = correct_histogram_x(h, p0=0, p1=calib_LO[(bar,side)], new_name=f"{h.GetName()}_LOcalib")
    h_energy[(bar, side, vov, thr)] = h_calib
    h_energy[(bar, side, vov, thr)].SetDirectory(0)
f.Close()
print("Energy calibrated per LO variations")
bars = sorted(bars)
sides = sorted(sides)
Vovs = sorted(Vovs)
thresholds = sorted(thresholds)

# -- fit Landau MPV for channel intercalibrations
mpv_l = {}
for key,h in h_energy.items():
    bar, side, vov, thr = key
    h_clone = h.Clone()
    f_landau, result = fit_landau_langaus(h_clone, min_energy_dict[(bar,vov)], 850, landau_only=True)
    mpv_l[key] = result["landau_mpv"]    
    # ------------------------------------------------------- 
    # optional : check fit
    # -------------------------------------------------------
    if args.fitCheck:
        c = ROOT.TCanvas("c","c",800,600)
        h_clone.SetLineColor(ROOT.kBlue)
        h_clone.GetXaxis().SetTitle("Energy")
        h_clone.GetYaxis().SetTitle("Entries")
        h_clone.Draw("hist")
        f_landau.SetLineColor(2)
        f_landau.Draw("same")
        c.SaveAs(f"{outdir_check}/energy_fit_bar{bar}{side}_vov{vov}_th{thr}.png")
        del c
    del f_landau

# -------------------------------------------------------
# STEP 3: derive TOFHIR+LO calibrations with MPV values
# ------------------------------------------------------- 
# -- derive intercalibration over MPV mean    
calib_tofhir_lo = {}
for vov in Vovs:
    for thr in thresholds:
        mpv_vals = []
        # mean computed on all channels
        for bar in bars:
            for side in ["L", "R"]:
                key = (bar, side, vov, thr)
                if key in mpv_l:
                    mpv_vals.append(mpv_l[key])        
        if not mpv_vals:
            print("[ERROR] MPV list not found")
            sys.exit()        
        mean_mpv = sum(mpv_vals) / len(mpv_vals)
        for bar in bars:
            for side in ["L", "R"]:
                key = (bar, side, vov, thr)
                if key in mpv_l:
                    calib_tofhir_lo[key] = mean_mpv / mpv_l[key]
                
# -- save tofhir+lo calibration factors to csv file
with open(TOFHIR_LO_csv,"w",newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["bar","side","vov","th","calib"])
    for key,val in calib_tofhir_lo.items():
        bar, side, vov, thr = key
        writer.writerow([bar, side, vov, thr, val])
print("Energy calibrated per TOFHIR+LO variations")

# -------------------------------------------------------
# STEP 4: isolate TOFHIR differen response (LO put back)
# -------------------------------------------------------
calib_tofhir = {}
for key, val in calib_tofhir_lo.items():
    bar, side, vov, thr = key
    # average cannot be calibrated, it will be computed per event in step1
    if side == "L-R":
        continue
    if (bar, side) not in calib_LO:
        print("[ERROR] Missing ", bar, " -  ", side, " in LO calibration factors.")
        sys.exit(0)    
    calib_tofhir[key] = val / calib_LO[(bar, side)]

# -- save to csv
with open(TOFHIR_csv, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["bar","side","vov","th","calib"])
    for key, val in calib_tofhir.items():
        bar, side, vov, thr = key
        writer.writerow([bar, side, vov, thr, val])
print("TOFHIR only calibration factors saved")

# -------------------------------------------------------
# optional : apply TOFHIR+LO calib and save to output file
# ------------------------------------------------------- 
if args.writeOutFile:
    f_out = ROOT.TFile(f"{inputdir}/energy_calib_TOFHIR_LO_{args.inputFile}.root","RECREATE")
    for key,h in h_energy.items():
        bar, side, vov, thr = key
        h_calib = h.Clone()
        if "L-R" != side:
            h_calib.Reset()
            h_calib = correct_histogram_x(h, p0=0, p1=calib_tofhir_lo[key], new_name=f"{h.GetName()}_TOFcalib")
        h_calib.Write()
    f_out.Close()
    print("Calibrated energy spectra saved in:", f_out.GetName())

# -------------------------------------------------------
# optional : draw TOFHIR calibrations
# ------------------------------------------------------- 
if args.drawCalib:
    plot_TOFHIR_calibration(TOFHIR_csv)

# ------------------------------------------------------- 
# optional : draw comparison energy spectra
# -------------------------------------------------------
if args.drawComparison:
    outdir_check = f"{outdir}/calibration_check"
    os.makedirs(outdir_check, exist_ok=True)
    for key,h_lo in h_energy.items():
        bar, side, vov, thr = key
        if side == "L-R":
            continue
        h_raw = h_energy_raw[key].Clone()
        h_localib = h_lo.Clone()
        h_tof = correct_histogram_x(h_lo, p0=0, p1=calib_tofhir_lo[key], new_name=f"{h_lo.GetName()}_TOFcalib")
        h_tofcalib = h_tof.Clone()
        integral = h_raw.Integral()
        rebin = get_rebin_factor(integral)
        h_raw.Rebin(rebin)
        h_localib.Rebin(rebin)
        h_tofcalib.Rebin(rebin)
        c = ROOT.TCanvas("c","c",800,600)
        h_raw.SetLineColor(ROOT.kBlack)
        h_localib.SetLineColor(ROOT.kBlue)
        h_tofcalib.SetLineColor(ROOT.kRed)
        h_raw.GetXaxis().SetTitle("Energy")
        h_raw.GetYaxis().SetTitle("Entries")
        h_raw.Draw("hist")
        h_localib.Draw("hist same")
        h_tofcalib.Draw("hist same")
        leg = ROOT.TLegend(0.6,0.7,0.88,0.88)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.AddEntry(h_raw,"Raw","l")
        leg.AddEntry(h_localib,"LO calib","l")
        leg.AddEntry(h_tofcalib,"LO + TOFHIR","l")
        leg.Draw()
        c.SaveAs(f"{outdir_check}/energy_bar{bar}{side}_vov{vov}_th{thr}.png")
        del c

# -------------------------------------------------------
# optional : draw MPV vs bar at different calibration stages
# -------------------------------------------------------
# -- MPV vs bar plots
if args.drawMPVvsBar:
    outdir_mpv = f"{outdir}/mpv_vs_bar"
    os.makedirs(outdir_mpv, exist_ok=True)
    mpv_raw = {}
    for key, h in h_energy_raw.items():
        bar, side, vov, thr = key
        if side == "L-R":
            continue
        h_clone = h.Clone()
        _,result = fit_landau_langaus(h_clone, min_energy_dict[(bar,vov)], 850, landau_only=True)
        mpv_raw[key] = result["landau_mpv"]

    mpv_lo_tof = {}
    for key, h_ in h_energy.items():
        bar, side, vov, thr = key
        if side == "L-R":
            continue
        h = h_.Clone()
        h_corr = correct_histogram_x(h, p0=0, p1=calib_tofhir_lo[key])
        _,result = fit_landau_langaus(h_corr, min_energy_dict[(bar,vov)], 850, landau_only=True)
        mpv_lo_tof[key] = result["landau_mpv"]

    mpv_tof_only = {}
    for key, h in h_energy_raw.items():
        bar, side, vov, thr = key
        if side == "L-R":
            continue
        if key not in calib_tofhir:
            continue
        h_corr = correct_histogram_x(h, p0=0, p1=calib_tofhir[key])
        _, result = fit_landau_langaus(h_corr, min_energy_dict[(bar,vov)], 850, landau_only=True)
        mpv_tof_only[key] = result["landau_mpv"]
        
    for vov in Vovs:
        for thr in thresholds:
            for side in ["L","R"]:
                c = ROOT.TCanvas("c","c",600,500)
                g_raw = ROOT.TGraph()
                g_lo = ROOT.TGraph()
                g_lo_tof = ROOT.TGraph()
                g_tof = ROOT.TGraph()
                for bar in bars:
                    key = (bar, side, vov, thr)
                    if key not in mpv_raw or key not in mpv_l or key not in mpv_lo_tof or key not in mpv_tof_only:
                        print("[WARNING] something brutto")
                        continue
                    g_raw.SetPoint(g_raw.GetN(), bar, mpv_raw[key])
                    g_lo.SetPoint(g_lo.GetN(), bar, mpv_l[key])
                    g_lo_tof.SetPoint(g_lo_tof.GetN(), bar, mpv_lo_tof[key])
                    g_tof.SetPoint(g_tof.GetN(), bar, mpv_tof_only[key])                    
                # style
                g_raw.SetLineColor(cms_colors[0])
                g_lo.SetLineColor(cms_colors[1])
                g_lo_tof.SetLineColor(cms_colors[2])
                g_tof.SetLineColor(cms_colors[3])
                g_raw.SetMarkerColor(cms_colors[0])
                g_lo.SetMarkerColor(cms_colors[1])
                g_lo_tof.SetMarkerColor(cms_colors[2])
                g_tof.SetMarkerColor(cms_colors[3])
                g_raw.SetMarkerStyle(20)
                g_lo.SetMarkerStyle(21)
                g_lo_tof.SetMarkerStyle(22)
                g_tof.SetMarkerStyle(23)                
                g_raw.SetTitle(f"{side} - Vov={vov} thr={thr};bar;MPV")
                g_raw.GetYaxis().SetRangeUser(150,550)
                g_raw.Draw("APL")
                g_lo.Draw("PL same")
                g_lo_tof.Draw("PL same")
                g_tof.Draw("PL same")
                leg = ROOT.TLegend(0.6,0.7,0.88,0.88)
                leg.SetBorderSize(0)
                leg.SetFillStyle(0)
                leg.AddEntry(g_raw,"Raw","lp")
                leg.AddEntry(g_lo,"LO calib","lp")
                leg.AddEntry(g_lo_tof,"LO + TOFHIR","lp")
                leg.AddEntry(g_tof,"TOFHIR only","lp")
                leg.Draw()
                c.SaveAs(f"{outdir_mpv}/mpv_vs_bar_{side}_vov{vov}_th{thr}.png")
                del c
