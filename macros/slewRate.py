#! /usr/bin/env python
import math
import ROOT

# ----- hard-coded values -----
ithmode = 0.313   # equivalence threshold amplitude in uA
max_xErr = 0.015  # max error accepted for pulse shape points
min_x = -5.       # min x value accepted for the pulse shape point
thMin = 1         # min index for low OV + irradiated modules (empirical fix for unstable omino rising edge)
vov_hc = 0.6      # ov below which we fix a threshold index
# --------------------------------


# --- Slew Rate extraction from pulse shape ----
# ----------------------------------------------
def getSlewRateFromPulseShape(g1, timingThreshold, npoints, gtemp, vov, name, canvas=None):
    """
    Fit the slope of the pulse's rising edge in a TGraph:
      1. Find the point corresponding to the chosen timing threshold
      2. Select a window of points around it
      3. Fit a linear function (pol1)
      4. Extract slope and its uncertainty
    Args:
        g1 (TGraph): The input TGraph representing the pulse.
        timingThreshold (float): The timing threshold value.
        npoints (int): Number of points to include in the fit.
        gtemp (TGraph): Temporary TGraph for fitting.
        canvas (TCanvas, optional): Canvas to display the fit (default: None).
    Returns:
        Tuple[float, float]: Slew rate and its error after excluding points with negative slopes.
    """
    # - Safety check 
    if g1.GetN() < npoints:
        print("[WARNING] Not enough points in the pulse shape")
        return (-1, -1)
    # - Find index at the timing threshold
    itiming = 0
    for i in range(0,g1.GetN()):
        if (round(g1.GetY()[i]/ithmode) == timingThreshold):
            itiming = i
            break
    # - NB: this check is actually ineffective (itiming starts at 0)
    #       kept for backward compatibility
    if itiming == -1:
        print("[WARNING] Timing threshold not found")
        return (-1, -1)
    # - Find min time index (start of the pulse)
    ifirst = ROOT.TMath.LocMin(g1.GetN(), g1.GetX())
    # - Define start index for fit 
    imin = max(0, itiming-2)
    if ( imin >= 0 and g1.GetX()[imin+1] < g1.GetX()[imin] ): imin = ifirst
    # - WARNING: hard-code for irradiated omini
    if float(vov)<=vov_hc and 'E1' in name and imin <thMin:
        imin = thMin
    tmin = g1.GetX()[imin]
    tmax = 4
    if tmax<=2:
        print("[WARNING] Insufficient number of points")
        return (-1,-1)
    # - Temporary graph for fitting 
    gtemp = ROOT.TGraphErrors()
    tmax = tmax-1
    if ((imin+npoints) < g1.GetN()): 
        tmax = min(g1.GetX()[imin+npoints],max(g1.GetX()))
        nmax = imin+npoints+1
    else:
        tmax = max(g1.GetX())
        nmax = g1.GetN()-1
    for i in range(imin, nmax):
        if g1.GetErrorX(i) < max_xErr and g1.GetX()[i]>min_x:
            gtemp.SetPoint(gtemp.GetN(), g1.GetX()[i], g1.GetY()[i])
            gtemp.SetPointError(gtemp.GetN()-1, g1.GetErrorX(i), g1.GetErrorY(i))
    # - Fit pulse shape with a linear func
    fitSR = ROOT.TF1('fitSR', 'pol1', tmin, tmax)
    fitSR.SetLineColor(g1.GetMarkerColor()+1)
    fitSR.SetRange(tmin,tmax)
    fitSR.SetParameters(0, 10)
    fitStatus = int(gtemp.Fit(fitSR, 'QRS+'))
    # - Extract slope of the pulse at the timing threshold
    sr = fitSR.Derivative( g1.GetX()[itiming])
    err_sr = fitSR.GetParError(1)
    if (err_sr == 0):
        err_sr = 0.05
    # - (optional) Draw canvas
    if (canvas!=None):
        canvas.cd()
        gtemp.SetMarkerStyle(g1.GetMarkerStyle())
        gtemp.SetMarkerColor(g1.GetMarkerColor()+2)
        gtemp.Draw('psames')
        g1.Draw('psames')
        fitSR.Draw('same')
        canvas.Update()
        ps = gtemp.FindObject("stats")
        # -- Check if 'ps' is not None and is of the TPaveStats class
        if ps and ps.IsA() == ROOT.TPaveStats.Class():
            # --- Cast 'ps' to TPaveStats
            ps.SetTextColor(g1.GetMarkerColor())
            if ('L' in g1.GetName()):
                ps.SetY1NDC(0.85) # new y start position
                ps.SetY2NDC(0.95)# new y end position
            if ('R' in g1.GetName()):
                ps.SetY1NDC(0.73) # new y start position
                ps.SetY2NDC(0.83)# new y end position
        else:
            print("[WARNING] Cannot find stats in the ps")
    return(sr,err_sr)


# -- Find optimal timing threshold (minimum deltaT) --
# ----------------------------------------------------
def findTimingThreshold(g2):
    """
    Scan deltaT vs threshold and return the threshold
    corresponding to the minimum timing resolution
    """
    xmin = 0
    ymin = 9999
    for i in range(0, g2.GetN()):
        y = g2.GetY()[i]
        x = g2.GetX()[i]
        if ( y < ymin ):
            ymin = y
            xmin = x
    return xmin


# --- Noise parametrization ---
# -----------------------------
def get_noise_pars(tofVersion):
    """
    Given a TOFHIR version, it returns the corresponding noise model parameters
    """
    if '2x' in tofVersion:
        n = 420.
        const = 16
        err_n = 35
        err_const = 2
    elif '2c' in tofVersion:
        n = 420.
        const = 16
        err_n = 35
        err_const = 2
    else:
        return None
    return n,err_n,const,err_const


def sigma_noise(sr, tofVersion, err_sr):
    """
    Compute the contribution to time resolution due to electronics noise, given a slew rate and TOFHIR version
    """
    if sr < 0 :
        return 0
    n, err_n, const, err_const = get_noise_pars(tofVersion)
    # Check if the parameters are valid (in case of unknown tofVersion)
    if n is None or const is None:
        return None    
    # Calculate the single noise value
    noise_single = math.sqrt(pow(n / sr, 2) + const ** 2)
    # Tot noise
    noise = noise_single / math.sqrt(2)
    err_noise = (1 / math.sqrt(2)) * (1 / noise) * (n / (sr**2)) * math.sqrt(math.pow(err_n, 2) + math.pow((n / sr) * err_sr, 2) + math.pow(const * err_const, 2))
    return noise, err_noise


# --- Time resolution extraction from histo ----
# ----------------------------------------------
def getTimeResolution(h1_deltaT):
    """
    Fit the deltaT histogram with a gaussian fit to extract sigma
      1. wide initial range 
      2. narrow down to the peak
      3. final fit
    """
    tRes = [-1,-1]
    h1_deltaT.GetXaxis().SetRangeUser(h1_deltaT.GetMean() - 5*h1_deltaT.GetRMS(), h1_deltaT.GetMean() + 5*h1_deltaT.GetRMS())   
    fitFunc = ROOT.TF1('fitFunc','gaus',-10000, 10000)
    fitFunc.SetLineColor(ROOT.kGreen+3)
    fitFunc.SetLineWidth(2)
    fitFunc.SetParameters(h1_deltaT.GetMaximum(),h1_deltaT.GetMean(), h1_deltaT.GetRMS())
    fitXMin = h1_deltaT.GetBinCenter(h1_deltaT.GetMaximumBin()) - 200
    fitXMax = h1_deltaT.GetBinCenter(h1_deltaT.GetMaximumBin()) + 200.
    fitFunc.SetRange(fitXMin, fitXMax)
    h1_deltaT.Fit('fitFunc','QNRL','', fitXMin, fitXMax)
    fitFunc.SetRange(fitFunc.GetParameter(1) - 1.0*fitFunc.GetParameter(2), fitFunc.GetParameter(1) + 1.0*fitFunc.GetParameter(2))
    h1_deltaT.Fit('fitFunc','QNRL')
    fitFunc.SetRange(fitFunc.GetParameter(1) - 2.5*fitFunc.GetParameter(2), fitFunc.GetParameter(1) + 2.5*fitFunc.GetParameter(2))
    h1_deltaT.Fit('fitFunc','QRSL+')
    tRes = [ fitFunc.GetParameter(2),fitFunc.GetParError(2)]
    return tRes
