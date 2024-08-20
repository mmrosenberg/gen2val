
import ROOT as rt
from math import sqrt


#Public Note files:
#f_CCnue = rt.TFile("systematics/CCnueInc_nuE/SBNfit_covariance_plots_CCnueInc.root")
#f_CCnumu = rt.TFile("systematics/CCnumuInc_nuE/SBNfit_covariance_plots_CCnumuInc.root")

#Final version files:
f_CCnue = rt.TFile("systematics/CCnueInc_nuE_chi2/SBNfit_covariance_plots_CCnueInc_nuE_chi2.root")
f_CCnumu = rt.TFile("systematics/CCnumuInc_nuE_rebinned/SBNfit_covariance_plots_CCnumuInc_nuE_rebinned.root")


flux_vars = ["expskin_FluxUnisim", "horncurrent_FluxUnisim", "nucleoninexsec_FluxUnisim", "nucleonqexsec_FluxUnisim", "nucleontotxsec_FluxUnisim", "pioninexsec_FluxUnisim", "pionqexsec_FluxUnisim", "piontotxsec_FluxUnisim", "kminus_PrimaryHadronNormalization", "kplus_PrimaryHadronFeynmanScaling", "kzero_PrimaryHadronSanfordWang", "piminus_PrimaryHadronSWCentralSplineVariation", "piplus_PrimaryHadronSWCentralSplineVariation"]

xsec_vars = ["All_UBGenie", "AxFFCCQEshape_UBGenie", "DecayAngMEC_UBGenie", "NormCCCOH_UBGenie", "NormNCCOH_UBGenie", "RPA_CCQE_UBGenie", "Theta_Delta2Npi_UBGenie", "ThetaDelta2NRad_UBGenie", "VecFFCCQEshape_UBGenie", "XSecShape_CCMEC_UBGenie", "xsr_scc_Fa3_SCC", "xsr_scc_Fv3_SCC"]

reint_vars = ["reinteractions_piminus_Geant4", "reinteractions_piplus_Geant4", "reinteractions_proton_Geant4"]


#Public note binning:
#CCnue_bins = [[0, 200], [200, 400], [400, 600], [600, 800], [800, 1000], [1000, 1200], [1200, 1400], [1400, 1600], [1600, 1800], [1800, 2000], [2000, 2200], [2200, 2400], [2400, 2600], [2600, 6000]]
#CCnumu_bins = [[0, 100], [100, 200], [200, 300], [300, 400], [400, 500], [500, 600], [600, 700], [700, 800], [800, 900], [900, 1000], [1000, 1100], [1100, 1200], [1200, 1300], [1300, 1400], [1400, 1500], [1500, 1600], [1600, 1700], [1700, 1800], [1800, 1900], [1900, 2000], [2000, 6000]]

#Final version binning:
CCnue_bins = [[0, 400], [400, 600], [600, 800], [800, 1000], [1000, 1200], [1200, 1400], [1400, 1600], [1600, 1800], [1800, 6000]]
CCnumu_bins = [[0, 100], [100, 200], [200, 300], [300, 400], [400, 500], [500, 600], [600, 700], [700, 800], [800, 900], [900, 1000], [1000, 1100], [1100, 1200], [1200, 1300], [1300, 1400], [1400, 1500], [1500, 1600], [1600, 1700], [1700, 1800], [1800, 6000]]


CCnue_flux = []
CCnue_xsec = []
CCnue_reint = []
CCnue_tot = []

CCnumu_flux = []
CCnumu_xsec = []
CCnumu_reint = []
CCnumu_tot = []


for i in range(1,len(CCnue_bins)+1,1):
  flux_err = 0.
  for flux_var in flux_vars:
    c = f_CCnue.Get("individualDir/collapsed fractional covariance matrix%s"%flux_var)
    h = c.GetPrimitive("coll_frac%s"%flux_var)
    flux_err += h.GetBinContent(i,i)
  CCnue_flux.append(sqrt(flux_err))
  xsec_err = 0.
  for xsec_var in xsec_vars:
    c = f_CCnue.Get("individualDir/collapsed fractional covariance matrix%s"%xsec_var)
    h = c.GetPrimitive("coll_frac%s"%xsec_var)
    xsec_err += h.GetBinContent(i,i)
  CCnue_xsec.append(sqrt(xsec_err))
  reint_err = 0.
  for reint_var in reint_vars:
    c = f_CCnue.Get("individualDir/collapsed fractional covariance matrix%s"%reint_var)
    h = c.GetPrimitive("coll_frac%s"%reint_var)
    reint_err += h.GetBinContent(i,i)
  CCnue_reint.append(sqrt(reint_err))
  CCnue_tot.append(sqrt(flux_err+xsec_err+reint_err))

for i in range(1,len(CCnumu_bins)+1,1):
  flux_err = 0.
  for flux_var in flux_vars:
    c = f_CCnumu.Get("individualDir/collapsed fractional covariance matrix%s"%flux_var)
    h = c.GetPrimitive("coll_frac%s"%flux_var)
    flux_err += h.GetBinContent(i,i)
  CCnumu_flux.append(sqrt(flux_err))
  xsec_err = 0.
  for xsec_var in xsec_vars:
    c = f_CCnumu.Get("individualDir/collapsed fractional covariance matrix%s"%xsec_var)
    h = c.GetPrimitive("coll_frac%s"%xsec_var)
    xsec_err += h.GetBinContent(i,i)
  CCnumu_xsec.append(sqrt(xsec_err))
  reint_err = 0.
  for reint_var in reint_vars:
    c = f_CCnumu.Get("individualDir/collapsed fractional covariance matrix%s"%reint_var)
    h = c.GetPrimitive("coll_frac%s"%reint_var)
    reint_err += h.GetBinContent(i,i)
  CCnumu_reint.append(sqrt(reint_err))
  CCnumu_tot.append(sqrt(flux_err+xsec_err+reint_err))


print()
print("CCnue flux =", CCnue_flux)
print()
print("CCnue xsec =", CCnue_xsec)
print()
print("CCnue reint =", CCnue_reint)
print()
print("CCnue tot =", CCnue_tot)
print()

print()
print("CCnumu flux =", CCnumu_flux)
print()
print("CCnumu xsec =", CCnumu_xsec)
print()
print("CCnumu reint =", CCnumu_reint)
print()
print("CCnumu tot =", CCnumu_tot)
print()

print()
print("CCnue reco neutrino energy")
print()
print("{\setlength{\\tabcolsep}{0.5em} \\renewcommand{\\arraystretch}{1.2}")
print("\\begin{table}[ht]")
print("    \centering")
#print("    \\resizebox{\\textwidth}{!}{")
print("    \\begin{tabular}{ |p{1.0in}||p{1.0in}|p{1.0in}|p{1.0in}||p{1.0in}| }")
#print("    \\begin{tabular}{ |c||c|c|c||c| }")
print("    \hline")
titleRow = "    {\\footnotesize \Centering{Reco $E_\\nu$ Bin (MeV)}} & "
#titleRow = "    {\\footnotesize Reco $E_\\nu$ Bin (MeV)} & "
for var in ["Flux","Cross Section","Hadron Reinteraction","Combined"]:
  if var == "Combined":
    titleRow += "\Centering{%s} \\\\"%var
    #titleRow += "%s \\\\"%var
  else:
    titleRow += "\Centering{%s} & "%var
    #titleRow += "%s & "%var
print(titleRow)
print("    \hline\hline")
for iB in range(len(CCnue_bins)):
  row = "    \Centering{%i - %i} & \Centering{%.3f} & \Centering{%.3f} & \Centering{%.3f} & \Centering{%.3f} \\\\"%(CCnue_bins[iB][0],CCnue_bins[iB][1],CCnue_flux[iB],CCnue_xsec[iB],CCnue_reint[iB],CCnue_tot[iB])
  #row = "    %i - %i & %.3f & %.3f & %.3f & %.3f \\\\"%(CCnue_bins[iB][0],CCnue_bins[iB][1],CCnue_flux[iB],CCnue_xsec[iB],CCnue_reint[iB],CCnue_tot[iB])
  print(row)
  print("    \hline")
print("    \end{tabular}")
#print("    }")
print("\end{table}")
print("}")

print()
print("CCnumu reco neutrino energy")
print()
print("{\setlength{\\tabcolsep}{0.5em} \\renewcommand{\\arraystretch}{1.2}")
print("\\begin{table}[ht]")
print("    \centering")
#print("    \\resizebox{\\textwidth}{!}{")
print("    \\begin{tabular}{ |p{1.0in}||p{1.0in}|p{1.0in}|p{1.0in}||p{1.0in}| }")
#print("    \\begin{tabular}{ |c||c|c|c||c| }")
print("    \hline")
titleRow = "    {\\footnotesize \Centering{Reco $E_\\nu$ Bin (MeV)}} & "
#titleRow = "    {\\footnotesize Reco $E_\\nu$ Bin (MeV)} & "
for var in ["Flux","Cross Section","Hadron Reinteraction","Combined"]:
  if var == "Combined":
    titleRow += "\Centering{%s} \\\\"%var
    #titleRow += "%s \\\\"%var
  else:
    titleRow += "\Centering{%s} & "%var
    #titleRow += "%s & "%var
print(titleRow)
print("    \hline\hline")
for iB in range(len(CCnumu_bins)):
  row = "    \Centering{%i - %i} & \Centering{%.3f} & \Centering{%.3f} & \Centering{%.3f} & \Centering{%.3f} \\\\"%(CCnumu_bins[iB][0],CCnumu_bins[iB][1],CCnumu_flux[iB],CCnumu_xsec[iB],CCnumu_reint[iB],CCnumu_tot[iB])
  #row = "    %i - %i & %.3f & %.3f & %.3f & %.3f \\\\"%(CCnumu_bins[iB][0],CCnumu_bins[iB][1],CCnumu_flux[iB],CCnumu_xsec[iB],CCnumu_reint[iB],CCnumu_tot[iB])
  print(row)
  print("    \hline")
print("    \end{tabular}")
#print("    }")
print("\end{table}")
print("}")


