import sys
import json

LYSO_PROPERTIES = {
    "LYSO818": {"sipm": "HPK_nonIrr", "cell": "HPK 25#mum", "thickness": 3.75, "light": 2410, "stoch_ref": 0, "ov_ref": 3.5, "color": 1, "tau_decay": 41, "tau_rise": 14, "SR_model": 7.5},
    "LYSO813": {"sipm": "HPK_nonIrr", "cell": "HPK 25#mum", "thickness": 3,    "light": 2250, "stoch_ref": 0, "ov_ref": 3.5, "color": 4, "tau_decay": 41, "tau_rise": 14, "SR_model": 7},
    "LYSO816": {"sipm": "HPK_nonIrr", "cell": "HPK 25#mum", "thickness": 2.4,  "light": 2190, "stoch_ref": 0, "ov_ref": 3.5, "color": 2, "tau_decay": 41, "tau_rise": 14, "SR_model": 6.5},
    "LYSO820": {"sipm": "HPK_nonIrr", "cell": "HPK 30#mum", "thickness": 3,    "light": 2400, "stoch_ref": 0,  "ov_ref": 1, "color": 42},

    "LYSO815": {"sipm": "HPK_2E14",   "cell": "HPK 25#mum", "thickness": 3,    "light": 2250, "stoch_ref": 30, "ov_ref": 1, "color": 6},
    "LYSO825": {"sipm": "HPK_2E14",   "cell": "HPK 20#mum", "thickness": 3,    "light": 2050, "stoch_ref": 35, "ov_ref": 1, "color": 8},
    "LYSO819": {"sipm": "HPK_1E14",   "cell": "HPK 25#mum", "thickness": 3.75, "light": 2410, "stoch_ref": 25, "ov_ref": 1, "color": 9},
    "LYSO817": {"sipm": "HPK_1E14",   "cell": "HPK 25#mum", "thickness": 2.4,  "light": 2190, "stoch_ref": 35, "ov_ref": 1, "color": 40},
    "LYSO829": {"sipm": "HPK_1E13",   "cell": "HPK 25#mum", "thickness": 3.75, "light": 2410, "stoch_ref": 35, "ov_ref": 1, "color": 44},
    "LYSO100056":{"sipm":"HPK_2E14",  "cell": "HPK 25#mum", "thickness":3.75,  "light": 2410, "stoch_ref": 25, "ov_ref":1, "color":46},
    "LYSO300032":{"sipm":"HPK_2E14",  "cell": "HPK 25#mum", "thickness":2.4,   "light": 2190, "stoch_ref": 35, "ov_ref":1, "color":48},
    "LYSO200104":{"sipm":"HPK_2E14",  "cell": "HPK 30#mum", "thickness":3,     "light": 2300, "stoch_ref": 30, "ov_ref":1, "color":50},
    
    "LYSO814": {"sipm": "HPK_nonIrr", "cell": "HPK 20#mum", "thickness":3,     "light": 2400, "stoch_ref": 0, "ov_ref": 3.5, "color":52},
    "LYSO528": {"sipm": "HPK_nonIrr", "cell": "HPK 15#mum", "thickness":3,     "light": 1500, "stoch_ref": 0, "ov_ref": 3.5, "color":54},
    "DM":      {"sipm": "HPK_nonIrr", "cell": "HPK 25#mum", "thickness":3.75, "light": 2410, "stoch_ref": 0, "ov_ref":3.5, "color":1}
}

# ----- Note-------
# - LYSO ID is used as SM identifier for pre-prod modules
# - For production modules, all specifications are the same, so only "DM" is used
# -----------------

# -----------------
# The following functions are used to return the above dictionary properties, given the LYSO ID
# -----------------
def lyso_(module):
    for key in LYSO_PROPERTIES:
        if key in module:
            return key
    return None

def get_prop(module, prop):
    key = lyso_(module)
    if key is None:
        return None
    return LYSO_PROPERTIES.get(key, {}).get(prop, None)

# getters
def sipm_(module): return get_prop(module, "sipm")
def sipm_cell_size(module): return get_prop(module, "cell")
def thickness(module): return get_prop(module, "thickness")
def light_output(module): return get_prop(module, "light")
def stoch_reference(module): return get_prop(module, "stoch_ref")
def ov_reference(module): return get_prop(module, "ov_ref")
def color_(module): return get_prop(module, "color")
def tau_decay(module): return get_prop(module, "tau_decay")
def tau_rise(module): return get_prop(module, "tau_rise")
def SR_model(module): return get_prop(module, "SR_model")

# labels
def irradiation(module):
    s = sipm_(module)
    if 'nonIrr' in s:
        return 'non-irr'
    elif 'E1' in s:
        return 'irr ' + s.split('_')[1]
    else:
        print("[ERROR]: which irradiation? ", module)
        return 0

def temperature_(module):
    if '_T' in module:
        return module.split("_T")[1].split("C")[0]
    print("[ERROR]: CANNOT FIND TEMPERATURE in ", module)
    return 'ERROR'

def angle_(module):
    if '_angle' in module:
        return module.split("_angle")[1]
    print("[ERROR]: CANNOT FIND ANGLE ", module)
    return 'ERROR'

def type_(module):
    t = thickness(module)
    if t == 3.75: return 1
    if t == 3:    return 2
    if t == 2.4:  return 3
    print("[ERROR]: TYPE NOT DEFINED ", module)
    return 'ERROR'

def Npe_frac(module):
    t = thickness(module)
    if t:
        t=t/3
    return t

def label_(module):
    if type_(module) and sipm_cell_size(module) and irradiation(module) and temperature_(module):
        return f"T{type_(module)} {sipm_cell_size(module)} {irradiation(module)} T{temperature_(module)}C"
    else:
        return None

# ----------------
# SiPM current was measured and corresponding info stored in json files
# the following functions allow to get the current, VovEff and DCR given
# - the module label
# - the overvoltage
# ----------------
def getVovEffDCR(module, ov) :
    if '829' in module or '819' in module:
        # Hard-coded: Sep23 currents were not reasonable so May23 data are taken instead
        #print("[WARNING] Taking SiPM currents from May TB data\n")
        currents = '/eos/cms/store/group/dpg_mtd/comm_mtd/TB/MTDTB_H8_May2023/VovsEff_v2.json'
    else:
        currents = '/eos/cms/store/group/dpg_mtd/comm_mtd/TB/MTDTB_H8_Sep2023/VovsEff_TOFHIR2C.json'
    ov_temp = round(5*round(float(ov)/5,2),2)
    ov_set = '%.2f'%float(ov_temp)
    irrad = irradiation(module)
    temp = temperature_(module)
    if 'non-irr' in irrad:
        return([ov_set,0])
    elif 'irr' in irrad:
        # Import file with VovEff and DCR
        with open(currents, 'r') as f:
            data = json.load(f)
        if not data[sipm_(module)+'_'+lyso_(module)+'_T'+temp+'C_A'][ov_set][0]:
            print('[ERROR]: module not found in the current json file ---  ',module)
            return 'ERROR'
        else:
            ov_eff_A = float(data[sipm_(module)+'_'+lyso_(module)+'_T'+temp+'C_A'][ov_set][0])
            dcr_A    = float(data[sipm_(module)+'_'+lyso_(module)+'_T'+temp+'C_A'][ov_set][1])
            ov_eff_B = float(data[sipm_(module)+'_'+lyso_(module)+'_T'+temp+'C_B'][ov_set][0])
            dcr_B    = float(data[sipm_(module)+'_'+lyso_(module)+'_T'+temp+'C_B'][ov_set][1])
            current_A= float(data[sipm_(module)+'_'+lyso_(module)+'_T'+temp+'C_A'][ov_set][2])
            current_B= float(data[sipm_(module)+'_'+lyso_(module)+'_T'+temp+'C_B'][ov_set][2])
            current = 0.5*(current_A+current_B)
            ov_eff =  0.5*(ov_eff_A+ov_eff_B)
            dcr    =  0.5*(dcr_A+dcr_B)
            return ([ov_eff, dcr, current]) 
    else:
        print('[ERROR]: CANNOT FIND WHICH IRRADIATION')
        return 'ERROR'

def Vovs_eff(module, ov):
    ov_set_ = '%.2f'%float(ov)
    if 'non-irr' in irradiation(module):
        VovEff_ = float(ov_set_)
    else:    
        VovEff_ = getVovEffDCR(module, ov_set_)[0]
    return VovEff_

def DCR(module, ov):
    ov_set_ = '%.2f'%float(ov)
    if 'non-irr' in irradiation(module):
        DCR_ = 0
    else:    
        DCR_ = getVovEffDCR(module,ov_set_)[1]
    return DCR_

def current_(module, ov):
    ov_set_ = '%.2f'%float(ov)
    if 'non-irr' in irradiation(module):
        cur_ = 0
    else:    
        cur_ = getVovEffDCR(module,ov_set_)[2]
    return cur_

def stoch_ref_fromFit(module):
    values = [0, 0]
    STOCH_REF_FIT = {
        # Type 1 + 25 um 2E14
        ("LYSO100056", "LYSO829", "LYSO819"): {
            "angle32": [41.55, 2.36],
            "angle52": [36.72, 2.17],
            "angle64": [34.07, 1.95],
        },
        # Type 1 + 15 um 2E14
        ("LYSO815",): {
            "angle64": [34.72, 0.42],
            "angle52": [38.06, 2.27],
        },
        # Type 3 + 24 um 2E14
        ("LYSO300032", "LYSO817"): {
            "angle64": [38.71, 0.37],
        },
        # Type 2 + 20 um 2E14
        ("LYSO825",): {
            "angle52": [39.13, 3.46],
        },
        # Type 2 + 30 um 2E14
        ("LYSO200104",): {
            "angle52": [34.85, 2.01],
        }
    }
    for modules, angle_dict in STOCH_REF_FIT.items():
        if any(m in module for m in modules):
            for angle_key, values in angle_dict.items():
                if angle_key in module:
                    return values
            print("[ERROR] : moduleDict - cannot find stoch ref for ", module)
    print("[ERROR] : moduleDict - cannot find stoch ref for ", module)
    return values

# -----------------
# The following dictionary provides the list of bars for which reasonable amount of events
# were recorded at TB. It was clearly a sub-optimal strategy but...
# -----------------

def good_bars(module, ovs, bars):
    good_bars_ = {}
    if '528' in module:
        good_bars_[3.50] = [1,2,3,4,5,6,7,8,9,10,11,12]
        good_bars_[2.00] = [1,2,3,4,5,6,7,8,9,10,11,12]
        good_bars_[1.50] = [1,2,3,4,5,6,7,8,9,10,11,12]
        good_bars_[1.25] = [1,2,3,4,5,6,7,8,9,10,11,12]
        good_bars_[1.00] = [1,2,3,4,5,6,7,8,9,10,11,12]
        good_bars_[0.80] = [1,2,3,4,5,6,7,8,9,10,11,12]
        good_bars_[0.60] = [1,2,6,7,8,9,10,11]
        good_bars_[0.50] = [1,2,6,7,8,9,10,11]
    elif '813' in module:
        good_bars_[3.50] = [7,8,9,10,11,12,13]
        good_bars_[2.00] = [7,8,9,10,11]
        good_bars_[1.50] = [7,8,9,10,11]
        good_bars_[1.25] = [7,8,9,10,11]
        good_bars_[1.00] = [7,8,9,10,11]
        good_bars_[0.80] = [7,8,9,10,11]
        good_bars_[0.60] = [8,9,10,11]
        good_bars_[0.50] = [8,9,10,11]
    elif '818' in module:
        good_bars_[3.50] = [8,9,10,11]
        good_bars_[2.00] = [6,7,8,9,10,11]
        good_bars_[1.50] = [6,7,8,9,10,11]
        good_bars_[1.25] = [6,7,8,9,10,11]
        good_bars_[1.00] = [6,7,8,9,10,11]
        good_bars_[0.80] = [7,8,9,10,11]
        good_bars_[0.60] = [7,8,9,10,11]
        good_bars_[0.50] = [7,8,9,10,11]
    elif '816' in module:
        good_bars_[3.50] = [3,4,5,6,7,8,9,10,11,12,13]
        good_bars_[2.00] = [2,3,4,5,6,7,8,9,10,11,12,13]
        good_bars_[1.50] = [2,3,4,5,6,7,8,9,10,11,12,13]
        good_bars_[1.25] = [2,3,4,5,6,7,8,9,10,11,12,13]
        good_bars_[1.00] = [2,3,4,5,6,7,8,9,10,11,12,13]
        good_bars_[0.80] = [2,3,4,5,6,7,8,9,10]
        good_bars_[0.60] = [2,3,4,5,6,7,8,9,10]
    elif '814' in module:
        good_bars_[3.50] = [2,3,4,5,6,7,8,9,10,11,12,13,14,15]
        good_bars_[2.00] = [2,3,4,5,6,7,8,9,10,11,12,13,14,15]
        good_bars_[1.50] = [2,3,4,5,6,7,8,9,10,11,12,13,14,15]
        good_bars_[1.25] = [2,3,4,5,6,7,8,9,10,11,12,13,14,15]
        good_bars_[1.00] = [2,3,4,5,6,7,8,9,10,11,12,13,14,15]
        good_bars_[0.80] = [2,3,4,5,6,7,8,9,10,11,12,13,14,15]
        good_bars_[0.60] = [2,3,4,5,6,7,8,9,10,11,12,13,14,15]
        good_bars_[0.50] = [2,3,4,5,6,7,8,9,10,11,12,13,14,15]
    elif '815' in module:
        good_bars_[2.00] = [6,7,8,9,10,11,12,13]
        good_bars_[1.50] = [6,7,8,9,10,11,12,13]
        good_bars_[1.25] = [6,7,8,9,10,11,12,13]
        good_bars_[1.00] = [6,7,8,9,10,11,12,13]
        good_bars_[0.80] = [6,7,8,9,10,11,12,13]
        good_bars_[0.60] = [6,7,8,9,10,11,12,13]
    elif '825' in module:
        good_bars_[2.00] = [7,8,9,10]
        good_bars_[1.50] = [7,8,9,10] 
        good_bars_[1.25] = [7,8,9,10]
        good_bars_[1.00] = [7,8,9,10]
        good_bars_[0.80] = [7,8,9,10]
        good_bars_[0.60] = [7,8,9,10]
    elif '819' in module:
        good_bars_[2.50] = [6,7,8,9,10,11,12,13]
        good_bars_[2.00] = [6,7,8,9,10,11,12,13]
        good_bars_[1.50] = [6,7,8,9,10,11,12,13]
        good_bars_[1.25] = [6,7,8,9,10,11,12,13]
        good_bars_[1.00] = [6,7,8,9,10,11,12,13]
        good_bars_[0.80] = [6,7,8,9,10,11,12,13]
        good_bars_[0.60] = [6,7,8,9,10,11,12,13]
    elif '829' in module:
        good_bars_[3.00] = [1,2,3,4,5,6,7,10,11,13,14,15]
        good_bars_[2.50] = [0,1,2,3,4,5,7,8,9,10,11,12]
        good_bars_[2.00] = [1,2,3,4,5,7,8,9,10,11,12]
        good_bars_[1.50] = [3,4,5,7,8,9,10,11,12,13,14,15]
        good_bars_[1.25] = [3,4,5,7,8,9,10,11,12]
        good_bars_[1.00] = [3,4,5,7,8,9,10,11,12,]
        good_bars_[0.80] = [3,4,5,7,8,9,11,12,13,15]
        good_bars_[0.60] = [3,4,5,7,8,9,11,12]
    elif '817' in module:
        good_bars_[3.50] = [2,3,4,5,6,7,8,9,10,11,12,13,14,15]
        good_bars_[2.00] = [2,3,4,5,6,7,8,9,10,11,12,13,14,15]
        good_bars_[1.50] = [2,3,4,5,6,7,8,9,10,11,12,13,14,15]
        good_bars_[1.25] = [2,3,4,5,6,7,8,9,10,11,12,13,14,15]
        good_bars_[1.00] = [2,3,4,5,6,7,8,9,10,11,12,13,14,15]
        good_bars_[0.80] = [2,3,4,5,6,7,8,9,10,11,12,13,14,15]
        good_bars_[0.60] = [2,3,4,5,6,7,8,9,10,11,12,13,14,15]
    elif '820' in module:
        good_bars_[3.50] = [7,8,9,10,11,12]
        good_bars_[2.50] = [7,8,9,10,11,12]
        good_bars_[2.00] = [7,8,9,10,11,12]
        good_bars_[1.50] = [6,7,8,9,10,11,12]
        good_bars_[1.25] = [6,7,8,9,10,11,12]
        good_bars_[1.00] = [6,7,8,9,10,11,12,13]
        good_bars_[0.80] = [6,7,8,9,10,11,12,13]
        good_bars_[0.60] = [6,7,8,9,10,11,12]
        good_bars_[0.50] = [6,7,8,9,10,11,12]
    elif '100056' in module:
        good_bars_[2.00] = [2,4,5,6,8,9,10,11,12,13]
        good_bars_[1.50] = [2,3,4,5,6,7,8,9,10,11,12,13]
        good_bars_[1.25] = [2,3,4,5,6,7,8,9,10,11,12,13]
        good_bars_[1.00] = [2,3,4,5,6,7,8,9,10,11,12,13]
        good_bars_[0.80] = [2,3,4,5,6,7,8,9,10,11,12,13]
        good_bars_[0.60] = [2,3,4,5,6,7,8,9,10,11,12,13]
    elif '300032' in module:
        good_bars_[2.00] = [9,10,11]
        good_bars_[1.50] = [6,7,8,9,10,11,12]
        good_bars_[1.25] = [6,7,8,9,10,11,12]
        good_bars_[1.00] = [6,7,8,9,10,11]
        good_bars_[0.80] = [6,7,8,9,10,11]
        good_bars_[0.60] = [6,7,8,9,10]
    elif '200104' in module:
        good_bars_[2.00] = [2,3,4,5,6,7,8,9,10,11,12,13,14]
        good_bars_[1.50] = [2,3,4,5,6,7,8,9,10,11,12,13,14]
        good_bars_[1.25] = [2,3,4,5,6,7,8,9,10,11,12,13,14]
        good_bars_[1.00] = [2,3,4,5,6,7,8,9,10,11,12,13,14]
        good_bars_[0.80] = [2,3,4,5,6,7,8,9,10,11,12,13,14]
        good_bars_[0.60] = [7,8,9,10,11,12]
    elif 'DM_348' in module:
        good_bars_[3.00] = [6,7,8,9,10]
        good_bars_[2.00] = [6,7,8,9,10]
        good_bars_[1.20] = [6,7,8,9,10]
        good_bars_[0.90] = [6,7,8,9,10]
    else:
        for vov in ovs:
            good_bars_[vov] = bars 
    return good_bars_
