#!/usr/bin/env python3
from utils import *
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)

# -------------------------------------------------------
# YOUR PATHS
eos_path = "/eos/home-c/cgiraldi/MTD/TB_CERN_Sep25/Lab5015Analysis"
# ------------------------------------------------------- 

# ------- parser -------
parser = argparse.ArgumentParser(description="Energy spectra calibration")
parser.add_argument("-sm", "--sensorModuleID", required=True, type=str, help="Sensor module ID for QAQC ROOT file")
parser.add_argument("--extraLabel", required=False, type=str, default=None, help="Input file has a name like module_{sensor_module_id}_analysis_{extraLabel}. default=None")
parser.add_argument("--draw", required=False, action="store_true", help="Produce calibration plots vs bar")
#addeed:
parser.add_argument("--mode", required=True, type=str, default="global",choices=["global", "separated_side"],help="Calibration mode: 'global' (default) or 'separated_side'")

args = parser.parse_args()
sensor_module_id = args.sensorModuleID
mode = args.mode
LO_dir = "/eos/cms/store/group/dpg_mtd/comm_mtd/TB/MTDTB_H8_Sep2025/SMs_QAQC/"
LO_csv = f"{eos_path}/plots/intercalibration/LO_calib/module_{sensor_module_id}_LO_calibration_factors_{mode}.csv"

if args.extraLabel is not None:
    extra_label = f"_{args.extraLabel}"
else:
    extra_label = ""
# ------------------------------------------------------- 
# STEP 1: LO calibration
#  produces a LO_calibration csv
# -------------------------------------------------------
f_qaqc = ROOT.TFile.Open(f"{LO_dir}/module_{sensor_module_id}_analysis{extra_label}.root")
g_LO_l = f_qaqc.Get("g_L_light_yield_vs_bar")
g_LO_r = f_qaqc.Get("g_R_light_yield_vs_bar")
LO_L = {}
LO_R = {}
for i in range(g_LO_l.GetN()):
    x, y = c_double(0.0), c_double(0.0)
    g_LO_l.GetPoint(i,x,y)
    bar = int(x.value)
    LO_L[bar] = float(y.value)
for i in range(g_LO_r.GetN()):
    x, y = c_double(0.0), c_double(0.0)
    g_LO_r.GetPoint(i,x,y)
    bar = int(x.value)
    LO_R[bar] = float(y.value)

#LO_LR = list(LO_L.values()) + list(LO_R.values())
#mean_LO = sum(LO_LR) / len(LO_LR)

calib_LO = {}

if mode == "global":
    print("Using GLOBAL calibration (L+R combined)")
    
    LO_LR = list(LO_L.values()) + list(LO_R.values())
    mean_LO = sum(LO_LR) / len(LO_LR)
    
    for bar in LO_L:
        calib_LO[(bar,"L")] = mean_LO / LO_L[bar]
        calib_LO[(bar,"R")] = mean_LO / LO_R[bar]

elif mode == "separated_side":
    print("Using SEPARATED SIDE calibration (L and R independent)")
    
    mean_LO_L = sum(LO_L.values()) / len(LO_L)
    mean_LO_R = sum(LO_R.values()) / len(LO_R)
    
    for bar in LO_L:
        calib_LO[(bar,"L")] = mean_LO_L / LO_L[bar]
        calib_LO[(bar,"R")] = mean_LO_R / LO_R[bar]


#----------------------------------------------
#for bar in LO_L:
#    calib_LO[(bar,"L")] = mean_LO / LO_L[bar]
#    calib_LO[(bar,"R")] = mean_LO / LO_R[bar]
#----------------------------------------------


# save LO calibration factor in a csv file
with open(LO_csv, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["bar","side","calib"])
    for key,val in calib_LO.items():
        bar, side = key
        writer.writerow([bar, side, val])
print(f"LO calibration factors saved in {LO_csv}")

# -------------------------------------------------------
# optional : draw LO calibrations
# ------------------------------------------------------- 
if args.draw == 1:
    plot_LO_calibration(LO_csv)
