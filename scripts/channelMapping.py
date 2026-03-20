###############################################################
# ----  map dictionary --------
# key  : bar ID
# item : [TOFHIR channel left ID, TOFHIR channel right ID]
# -----------------------------
# ---- notes for reference ----
# barID = x
#    SiPM channel left ID = x
#    SiPM channel right ID = 31 - x
# the TOFHIR channels ID mapping is reported below
# Please note that you should account for the ASIC ID as a rigid shift in the overall channelID : 32*ASIC_ID
# e.g.
#      for ASIC 7 as REF
#       bar 7 ---> 32*7 + channelMapping[7]
#        - chL 233
#        - chR 246
#      for ASIC 4 as REF
#       bar 7
#        - chL 137
#        - chR 150
# -----------------------------
###############################################################
map_bar_LR = {
    0  : [4, 29] ,
    1  : [1, 25] ,
    2  : [0, 24] ,
    3  : [3, 30] ,
    4  : [2, 31] ,
    5  : [6, 28] ,
    6  : [7, 27] ,
    7  : [9, 22] ,
    8  : [5, 26] ,
    9  : [11, 21],
    10 : [8, 23] ,
    11 : [12, 20],
    12 : [10, 19],
    13 : [14, 18],
    14 : [15, 16],
    15 : [13, 17]    
}


def get_TOFHIR_channel_mapping(mapping=map_bar_LR):
    result = []
    # sorted ensures that channels are specified from bar 0 to 15
    for k in sorted(mapping):
        pair = mapping[k]
        if len(pair) != 2:
            raise ValueError(f"Channel mapping for bar {k} does not have exactly 2 elements")
        result.extend(pair)
    return " ".join(map(str, result))
