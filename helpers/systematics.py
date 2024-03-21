
from math import sqrt

CCnueInc_recoNuE_noDet = {
  1:  (0.20419348795344358, 0,    200),
  2:  (0.16270166588632198, 200,  400),
  3:  (0.13507890615510365, 400,  600),
  4:  (0.1387272640038089,  600,  800),
  5:  (0.1424787152580893,  800,  1000),
  6:  (0.14760382763631938, 1000, 1200),
  7:  (0.15605260148049288, 1200, 1400),
  8:  (0.15154633722007868, 1400, 1600),
  9:  (0.15893113706387335, 1600, 1800),
  10: (0.16363994629672654, 1800, 2000),
  11: (0.16377913426655064, 2000, 2200),
  12: (0.17522802582593777, 2200, 2400),
  13: (0.19594277255742099, 2400, 2600),
  14: (0.17968501940462772, 2600, 6000)
}

CCnumuInc_recoNuE_noDet = {
  1:  (0.2006382369967962,  0,    100),
  2:  (0.1688331256214487,  100,  200),
  3:  (0.14190381253341172, 200,  300),
  4:  (0.1410299292207514,  300,  400),
  5:  (0.13453670592879932, 400,  500),
  6:  (0.14831409448047111, 500,  600),
  7:  (0.14703440410424276, 600,  700),
  8:  (0.16493348337160468, 700,  800),
  9:  (0.16394398476878252, 800,  900),
  10: (0.17824939038004237, 900,  1000),
  11: (0.18255422054531917, 1000, 1100),
  12: (0.18593315102626304, 1100, 1200),
  13: (0.2019789337513537,  1200, 1300),
  14: (0.1938856019446503,  1300, 1400),
  15: (0.19289081192305627, 1400, 1500),
  16: (0.21769079989568502, 1500, 1600),
  17: (0.20061425353248324, 1600, 1700),
  18: (0.1917243622936701,  1700, 1800),
  19: (0.1963899364894271,  1800, 1900),
  20: (0.21682219177642362, 1900, 2000),
  21: (0.17564528692156936, 2000, 6000)
}


def SetUncertainties(h_pred, var, sel):

  if sel not in ["CCnue","CCnumu"]:
    print("%s not a valid selection, only \"CCnue\" and \"CCnumu\" supported, plotting stats only error!"%sel)
    return h_pred
  if var == "recoNuE":
    if sel == "CCnue":
      sys = CCnueInc_recoNuE_noDet
    if sel == "CCnumu":
      sys = CCnumuInc_recoNuE_noDet
  else:
    print("no systematics available for %s, plotting stats only error!"%var)
    return h_pred

  #TODO: add more detailed bin range check
  if h_pred.GetNbinsX() != len(sys):
    print("for selection %s, variable %s number of bins in systematics calculation not equal to histo bin count, plotting stats only error!"%(sel, var))
    return h_pred

  print("adding systematics uncertainties for selection %s, variable %s"%(sel, var))

  for i in range(1, len(sys)+1, 1):
    stats_err = h_pred.GetBinError(i)
    sys_err = h_pred.GetBinContent(i)*sys[i][0]
    tot_err = sqrt(stats_err**2 + sys_err**2)
    print("stats, sys, total errors for bin %i: %.2f, %.2f, %.2f"%(i,stats_err, sys_err, tot_err))
    h_pred.SetBinError(i, tot_err)
    
  return h_pred


