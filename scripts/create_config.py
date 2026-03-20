#!/usr/bin/env python
import shutil
from pathlib import Path
import argparse
import sys
from channelMapping import *

# ------ MODIFY YOUR PATHS HERE -------
Lab5015_path = "/eos/home-s/spalluot/MTD/TB_CERN_Sep25/Lab5015Analysis"
plot_path = "/eos/home-s/spalluot/www/MTD/MTDTB_CERN_Sep25"
# -------------------------------------


# parser
# ----------------------------
parser = argparse.ArgumentParser(description='This script creates moduleCharacterization , drawPulseShape, and minEnergy cfg')
parser.add_argument("-ml", "--modulelabel", required=True, type=str, help="module label")
parser.add_argument("-r", "--runs", required=True, type=str, help="comma-separated list of runs to be processed")
parser.add_argument("-t", "--temperature", required=True, type=str, help="temperature")
parser.add_argument("-ov", "--Vov", required=True, type=str, help="overvoltage")
parser.add_argument("-c", "--config", required=True, type=str, help="config number")
parser.add_argument("-th", "--threshold", required=True, type=str, help="threshold used for the scan")
parser.add_argument("-vth1", "--thresholdT1", required=False, default="-1", type=str, help="specific threshold for T1")
parser.add_argument("-vth2", "--thresholdT2", required=False, default="-1", type=str, help="specific threshold for T2")
parser.add_argument("-vthe", "--thresholdE",  required=False, default="-1", type=str, help="specific threshold for E")
parser.add_argument("-e", "--extraLabel", required=False,  type=str, help="eg: angle or check or whatever")
parser.add_argument("--whichEnergyIntercalib", required=False,  default="", type=str, help="i.e. TOFHIR, TOFHIR_LO or empty for None")
parser.add_argument("--dutASIC", required=False, default=7, type=int, help="the DUT ASIC position (default set to 7)")
parser.add_argument("--refASIC", required=False, default=4, type=int, help="the REF ASIC position (default set to 4)")
parser.add_argument("--refBar",  required=False, type=int, help="the REF bar on which ask for coincidence (default set to 7 in the code)")
parser.add_argument("--saveRefInfoFlag", required=False, type=int, default=0, help="0 or 1: flag to set if reference info should be saved")
args = parser.parse_args()

# -- changing the options on REF and DUT ASICs might not be necessary for all the studies
# -   default values are set in the argparse
if args.refBar is None:
    args.refBar = 7
else:
    # - if the reference bar is not the default one, add an extra label
    if args.extraLabel:
        args.extraLabel += "refBar"+str(args.refBar)
    else:
        args.extraLabel = "refBar"+str(args.refBar)
chL = args.refASIC*32 + map_bar_LR[args.refBar][0]
chR = args.refASIC*32 + map_bar_LR[args.refBar][1]
print("channel left ", chL, "   channel right ", chR)

# if calibration factors are to be included, add an extralabel
if args.whichEnergyIntercalib != "":
    args.extraLabel = f"{args.whichEnergyIntercalib}calib"
    intercalib_path = f"Plot_repo_path/energy_intercalibration/baseGeneralLabel/whichCalibration_calibration_factors.csv"
else:
    intercalib_path = "0"
    
# config label
# ----------------------------
Vov_float = float(args.Vov)
base_label = f"{args.modulelabel}_Vov{Vov_float:.2f}_T{args.temperature}C"
if args.extraLabel:
    label = f"{base_label}_{args.extraLabel}"
else:
    label = base_label

# path
# ----------------------------
cfg_path = f"{Lab5015_path}/cfg"
cfgFolder = Path(cfg_path)
temp_min = cfgFolder / f"minEnergies_{args.modulelabel}.txt"


# copy minEnergy file if it does not exist
# ----------------------------
if not temp_min.exists():
    base_min = cfgFolder / "minEnergies_base.txt"
    shutil.copy(base_min, temp_min)

    
# function to write config file
# -------------------------------
def write_cfg(base_path: Path, out_path: Path, replacements: dict, check_Vov: bool = False):
    with open(base_path) as baseCfg, open(out_path, 'w') as newCfg:
        for line in baseCfg:
            original_line = line
            if check_Vov and line.startswith('Vov') and args.Vov not in line:
                print(f"ERROR: missing ov in {base_path.name} file")
                sys.exit(1)

            for key, val in replacements.items():
                if key in line:
                    line = line.replace(key, str(val))
            newCfg.write(line)


# ModuleCharacterization cfg
# ----------------------------
base_module_cfg = cfgFolder / "moduleCharacterization_base.cfg"
out_module_cfg = cfgFolder / f"moduleCharacterization_{label}.cfg"
print(f"writing \t {out_module_cfg.name}")

# replace the words find in the cfg_base (key) with the item of this dict
replacements_module = {
    "interCalibPath" : intercalib_path,
    "Lab5015_repo_path": Lab5015_path,
    "Plot_repo_path": plot_path,
    "runNumbers": args.runs,
    "generalLabel": label,
    "moduleLabel": args.modulelabel,
    "confNumber": args.config,
    "vovLabel": args.Vov,
    "scanVth" : args.threshold,
    "setVth1" : args.thresholdT1,
    "setVth2" : args.thresholdT2,
    "setVthe" : args.thresholdE,
    "channelLeft" : chL,
    "channelRight": chR,
    "refASIC_ch" : args.refASIC,
    "dutASIC_ch" : args.dutASIC,
    "saveReferenceModuleInfoFlag" : args.saveRefInfoFlag,
    "channelMapping_list_from_bar0_to_bar15" : get_TOFHIR_channel_mapping(),
    "whichCalibration" : args.whichEnergyIntercalib,
    "baseGeneralLabel" : base_label
}

write_cfg(base_module_cfg, out_module_cfg, replacements_module, check_Vov=True)
print(f"config : {args.config}")


# drawPulseShapeTB.cpp needs few adjustments for the new structure with DUT and REF ASIC settings

# # drawPulseShapeTB cfg
# # ----------------------------
# base_pulse_cfg = cfgFolder / "drawPulseShapeTB_base.cfg"
# out_pulse_cfg = cfgFolder / f"drawPulseShapeTB_{label}.cfg"
# print(f"writing \t {out_pulse_cfg.name}")

# replacements_pulse = {
#     "runNumbers": args.runs,
#     "generalLabel": label,
#     "moduleLabel": args.modulelabel,
#     "confNumber": args.config
