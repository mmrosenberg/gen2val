
import argparse
import ROOT as rt
from math import sqrt
from helpers.plotting_functions import sortHists

parser = argparse.ArgumentParser("combine individual detvar covariance matrices to get full det systematics")
parser.add_argument("-t", "--intag", required=True, type=str, help="tag for input files")
parser.add_argument("--ccnumuOverflow", action="store_true", help="use ccnumu overflow binning")
parser.add_argument("--ccnueChi2", action="store_true", help="use ccnue chi2 binning")
args = parser.parse_args()

def getOverflowLabel(h):
  x = h.GetBinCenter(h.GetNbinsX()) - 0.1*h.GetBinWidth(h.GetNbinsX())
  rt.gPad.Update()
  y = -rt.gPad.GetFrame().GetY2()/27.
  label = rt.TText(x,y,"overflow")
  label.SetTextSize(0.028)
  label.SetTextAngle(-40)
  return label

#detvars = ["LYAtt-bugFix", "LYDown", "LYRayleigh", "recomb2", "SCE", "wiremodThetaXZ", "wiremodThetaYZ", "wiremodX", "wiremodYZ"]
detvars = ["LYAtt", "LYDown", "LYRayleigh", "recomb2", "SCE", "wiremodThetaXZ", "wiremodThetaYZ", "wiremodX", "wiremodYZ"]
detvarsTable = ["LY Att.", "LY Down", "LY Rayleigh", "recomb2", "SCE", "Wiremod $\\theta_{xz}$", "Wiremod $\\theta_{yz}$", "Wiremod $x$", "Wiremod $y,z$"]

rt.gStyle.SetOptStat(0)

if args.ccnueChi2:
  CCnue_frac_covar_nuE = rt.TH2F("CCnue_frac_covar_nuE","Fractional Covariance Matrix",9,0,9,9,0,9)
  CCnue_frac_covar_lepP = rt.TH2F("CCnue_frac_covar_lepP","Fractional Covariance Matrix",8,0,8,8,0,8)
  CCnue_frac_covar_lepCosTheta = rt.TH2F("CCnue_frac_covar_lepCosTheta","Fractional Covariance Matrix",6,0,6,6,0,6)
else:
  CCnue_frac_covar_nuE = rt.TH2F("CCnue_frac_covar_nuE","Fractional Covariance Matrix",14,0,14,14,0,14)
  CCnue_frac_covar_lepP = rt.TH2F("CCnue_frac_covar_lepP","Fractional Covariance Matrix",14,0,14,14,0,14)
  CCnue_frac_covar_lepCosTheta = rt.TH2F("CCnue_frac_covar_lepCosTheta","Fractional Covariance Matrix",16,0,16,16,0,16)
if args.ccnumuOverflow:
  CCnumu_frac_covar_nuE = rt.TH2F("CCnumu_frac_covar_nuE","Fractional Covariance Matrix",19,0,19,19,0,19)
  CCnumu_frac_covar_lepP = rt.TH2F("CCnumu_frac_covar_lepP","Fractional Covariance Matrix",16,0,16,16,0,16)
  #CCnumu_frac_covar_nuE = rt.TH2F("CCnumu_frac_covar_nuE","Fractional Covariance Matrix",15,0,15,15,0,15)
  #CCnumu_frac_covar_lepP = rt.TH2F("CCnumu_frac_covar_lepP","Fractional Covariance Matrix",13,0,13,13,0,13)
  #CCnumu_frac_covar_lepCosTheta = rt.TH2F("CCnumu_frac_covar_lepCosTheta","Fractional Covariance Matrix",14,0,14,14,0,14)
else:
  CCnumu_frac_covar_nuE = rt.TH2F("CCnumu_frac_covar_nuE","Fractional Covariance Matrix",21,0,21,21,0,21)
  CCnumu_frac_covar_lepP = rt.TH2F("CCnumu_frac_covar_lepP","Fractional Covariance Matrix",21,0,21,21,0,21)
CCnumu_frac_covar_lepCosTheta = rt.TH2F("CCnumu_frac_covar_lepCosTheta","Fractional Covariance Matrix",16,0,16,16,0,16)

outfile = rt.TFile("systematics/detvar/detvar_spectra_and_combined_fractional_covariance_matrices_%s.root"%args.intag,"RECREATE")

CCnue_nuE_errors = {}
CCnumu_nuE_errors = {}
CCnue_lepP_errors = {}
CCnumu_lepP_errors = {}
CCnue_lepCosTheta_errors = {}
CCnumu_lepCosTheta_errors = {}

def histAbs(h):
  hnew = h.Clone("hnew")
  for i in range(1,h.GetNbinsX()+1,1):
    for j in range(1,h.GetNbinsY()+1,1):
      hnew.SetBinContent(i,j,abs(h.GetBinContent(i,j)))
  hnew.SetTitle("Collapsed fractional covariance matrix")
  hnew.GetXaxis().SetTitle("Reco Bin i")
  hnew.GetYaxis().SetTitle("Reco Bin j")
  return hnew

for detvar in detvars:
  f = rt.TFile("systematics/detvar/calculate_detvar_systematics_%s_%s_output.root"%(detvar,args.intag), "READ")
  CCnue_nuE = f.Get("CCnue_frac_covar_recoNuE")
  CCnumu_nuE = f.Get("CCnumu_frac_covar_recoNuE")
  CCnue_lepP = f.Get("CCnue_frac_covar_recoLepP")
  CCnumu_lepP = f.Get("CCnumu_frac_covar_recoLepP")
  CCnue_lepCosTheta = f.Get("CCnue_frac_covar_recoLepCosTheta")
  CCnumu_lepCosTheta = f.Get("CCnumu_frac_covar_recoLepCosTheta")
  CCnue_nuE_errors[detvar] = []
  CCnumu_nuE_errors[detvar] = []
  for i in range(1,CCnue_nuE.GetNbinsX()+1,1):
    CCnue_nuE_errors[detvar].append(sqrt(CCnue_nuE.GetBinContent(i,i)))
  for i in range(1,CCnumu_nuE.GetNbinsX()+1,1):
    CCnumu_nuE_errors[detvar].append(sqrt(CCnumu_nuE.GetBinContent(i,i)))
  CCnue_lepP_errors[detvar] = []
  CCnumu_lepP_errors[detvar] = []
  for i in range(1,CCnue_lepP.GetNbinsX()+1,1):
    CCnue_lepP_errors[detvar].append(sqrt(CCnue_lepP.GetBinContent(i,i)))
  for i in range(1,CCnumu_lepP.GetNbinsX()+1,1):
    CCnumu_lepP_errors[detvar].append(sqrt(CCnumu_lepP.GetBinContent(i,i)))
  CCnue_lepCosTheta_errors[detvar] = []
  CCnumu_lepCosTheta_errors[detvar] = []
  for i in range(1,CCnue_lepCosTheta.GetNbinsX()+1,1):
    CCnue_lepCosTheta_errors[detvar].append(sqrt(CCnue_lepCosTheta.GetBinContent(i,i)))
  for i in range(1,CCnumu_lepCosTheta.GetNbinsX()+1,1):
    CCnumu_lepCosTheta_errors[detvar].append(sqrt(CCnumu_lepCosTheta.GetBinContent(i,i)))
  #CCnue_frac_covar_nuE.Add(CCnue_nuE)
  #CCnumu_frac_covar_nuE.Add(CCnumu_nuE)
  #CCnue_frac_covar_lepP.Add(CCnue_lepP)
  #CCnumu_frac_covar_lepP.Add(CCnumu_lepP)
  #CCnue_frac_covar_lepCosTheta.Add(CCnue_lepCosTheta)
  #CCnumu_frac_covar_lepCosTheta.Add(CCnumu_lepCosTheta)
  print("Adding CCnue nuE covar for %s"%detvar)
  CCnue_frac_covar_nuE.Add(CCnue_nuE)
  print("Adding CCnumu nuE covar for %s"%detvar)
  CCnumu_frac_covar_nuE.Add(CCnumu_nuE)
  print("Adding CCnue lepP covar for %s"%detvar)
  CCnue_frac_covar_lepP.Add(CCnue_lepP)
  print("Adding CCnumu lepP covar for %s"%detvar)
  CCnumu_frac_covar_lepP.Add(CCnumu_lepP)
  print("Adding CCnue lepCosTheta covar for %s"%detvar)
  CCnue_frac_covar_lepCosTheta.Add(CCnue_lepCosTheta)
  print("Adding CCnumu lepCosTheta covar for %s"%detvar)
  CCnumu_frac_covar_lepCosTheta.Add(CCnumu_lepCosTheta)
  CCnue_nuE_spec_CV = f.Get("h_CCnue_recoNuE_CV")
  CCnumu_nuE_spec_CV = f.Get("h_CCnumu_recoNuE_CV")
  CCnue_lepP_spec_CV = f.Get("h_CCnue_recoLepP_CV")
  CCnumu_lepP_spec_CV = f.Get("h_CCnumu_recoLepP_CV")
  CCnue_lepCosTheta_spec_CV = f.Get("h_CCnue_recoLepCosTheta_CV")
  CCnumu_lepCosTheta_spec_CV = f.Get("h_CCnumu_recoLepCosTheta_CV")
  CCnue_nuE_spec_CV.SetLineColor(rt.kBlue)
  CCnue_nuE_spec_CV.SetLineWidth(2)
  CCnumu_nuE_spec_CV.SetLineColor(rt.kBlue)
  CCnumu_nuE_spec_CV.SetLineWidth(2)
  CCnue_lepP_spec_CV.SetLineColor(rt.kBlue)
  CCnue_lepP_spec_CV.SetLineWidth(2)
  CCnumu_lepP_spec_CV.SetLineColor(rt.kBlue)
  CCnumu_lepP_spec_CV.SetLineWidth(2)
  CCnue_lepCosTheta_spec_CV.SetLineColor(rt.kBlue)
  CCnue_lepCosTheta_spec_CV.SetLineWidth(2)
  CCnumu_lepCosTheta_spec_CV.SetLineColor(rt.kBlue)
  CCnumu_lepCosTheta_spec_CV.SetLineWidth(2)
  CCnue_nuE_spec_var = f.Get("h_CCnue_recoNuE_var")
  CCnumu_nuE_spec_var = f.Get("h_CCnumu_recoNuE_var")
  CCnue_lepP_spec_var = f.Get("h_CCnue_recoLepP_var")
  CCnumu_lepP_spec_var = f.Get("h_CCnumu_recoLepP_var")
  CCnue_lepCosTheta_spec_var = f.Get("h_CCnue_recoLepCosTheta_var")
  CCnumu_lepCosTheta_spec_var = f.Get("h_CCnumu_recoLepCosTheta_var")
  CCnue_nuE_spec_var.SetLineColor(rt.kRed)
  CCnue_nuE_spec_var.SetLineWidth(2)
  CCnumu_nuE_spec_var.SetLineColor(rt.kRed)
  CCnumu_nuE_spec_var.SetLineWidth(2)
  CCnue_lepP_spec_var.SetLineColor(rt.kRed)
  CCnue_lepP_spec_var.SetLineWidth(2)
  CCnumu_lepP_spec_var.SetLineColor(rt.kRed)
  CCnumu_lepP_spec_var.SetLineWidth(2)
  CCnue_lepCosTheta_spec_var.SetLineColor(rt.kRed)
  CCnue_lepCosTheta_spec_var.SetLineWidth(2)
  CCnumu_lepCosTheta_spec_var.SetLineColor(rt.kRed)
  CCnumu_lepCosTheta_spec_var.SetLineWidth(2)
  CCnue_nuE_hists = sortHists([CCnue_nuE_spec_CV, CCnue_nuE_spec_var])
  CCnue_lepP_hists = sortHists([CCnue_lepP_spec_CV, CCnue_lepP_spec_var])
  CCnue_lepCosTheta_hists = sortHists([CCnue_lepCosTheta_spec_CV, CCnue_lepCosTheta_spec_var])
  CCnumu_nuE_hists = sortHists([CCnumu_nuE_spec_CV, CCnumu_nuE_spec_var])
  CCnumu_lepP_hists = sortHists([CCnumu_lepP_spec_CV, CCnumu_lepP_spec_var])
  CCnumu_lepCosTheta_hists = sortHists([CCnumu_lepCosTheta_spec_CV, CCnumu_lepCosTheta_spec_var])
  outfile.cd()
  cnv_spec_CCnue_nuE = rt.TCanvas("cnv_spec_CCnue_%s_nuE"%detvar,"cnv_spec_CCnue_%s_nuE"%detvar)
  CCnue_nuE_hists[0].Draw("EHIST")
  CCnue_nuE_hists[1].Draw("EHISTSAME")
  leg_CCnue_nuE = rt.TLegend(0.7,0.7,0.9,0.9)
  leg_CCnue_nuE.AddEntry(CCnue_nuE_spec_CV, "CV", "l")
  leg_CCnue_nuE.AddEntry(CCnue_nuE_spec_var, ("%s"%detvar).replace("Att-BugFix","Attenuation"), "l")
  leg_CCnue_nuE.Draw()
  cnv_spec_CCnue_nuE.Write()
  cnv_covar_CCnue_nuE = rt.TCanvas("cnv_covar_CCnue_%s_nuE"%detvar,"cnv_covar_CCnue_%s_nuE"%detvar)
  CCnue_nuE.Draw("COLZ")
  cnv_covar_CCnue_nuE.Write()
  cnv_spec_CCnue_lepP = rt.TCanvas("cnv_spec_CCnue_%s_lepP"%detvar,"cnv_spec_CCnue_%s_lepP"%detvar)
  CCnue_lepP_hists[0].Draw("EHIST")
  CCnue_lepP_hists[1].Draw("EHISTSAME")
  leg_CCnue_lepP = rt.TLegend(0.7,0.7,0.9,0.9)
  leg_CCnue_lepP.AddEntry(CCnue_lepP_spec_CV, "CV", "l")
  leg_CCnue_lepP.AddEntry(CCnue_lepP_spec_var, ("%s"%detvar).replace("Att-BugFix","Attenuation"), "l")
  leg_CCnue_lepP.Draw()
  cnv_spec_CCnue_lepP.Write()
  cnv_covar_CCnue_lepP = rt.TCanvas("cnv_covar_CCnue_%s_lepP"%detvar,"cnv_covar_CCnue_%s_lepP"%detvar)
  CCnue_lepP.Draw("COLZ")
  cnv_covar_CCnue_lepP.Write()
  cnv_spec_CCnue_lepCosTheta = rt.TCanvas("cnv_spec_CCnue_%s_lepCosTheta"%detvar,"cnv_spec_CCnue_%s_lepCosTheta"%detvar)
  CCnue_lepCosTheta_hists[0].Draw("EHIST")
  CCnue_lepCosTheta_hists[1].Draw("EHISTSAME")
  leg_CCnue_lepCosTheta = rt.TLegend(0.7,0.7,0.9,0.9)
  leg_CCnue_lepCosTheta.AddEntry(CCnue_lepCosTheta_spec_CV, "CV", "l")
  leg_CCnue_lepCosTheta.AddEntry(CCnue_lepCosTheta_spec_var, ("%s"%detvar).replace("Att-BugFix","Attenuation"), "l")
  leg_CCnue_lepCosTheta.Draw()
  cnv_spec_CCnue_lepCosTheta.Write()
  cnv_covar_CCnue_lepCosTheta = rt.TCanvas("cnv_covar_CCnue_%s_lepCosTheta"%detvar,"cnv_covar_CCnue_%s_lepCosTheta"%detvar)
  CCnue_lepCosTheta.Draw("COLZ")
  cnv_covar_CCnue_lepCosTheta.Write()
  cnv_spec_CCnumu_nuE = rt.TCanvas("cnv_spec_CCnumu_%s_nuE"%detvar,"cnv_spec_CCnumu_%s_nuE"%detvar)
  CCnumu_nuE_hists[0].Draw("EHIST")
  CCnumu_nuE_hists[1].Draw("EHISTSAME")
  leg_CCnumu_nuE = rt.TLegend(0.7,0.7,0.9,0.9)
  leg_CCnumu_nuE.AddEntry(CCnumu_nuE_spec_CV, "CV", "l")
  leg_CCnumu_nuE.AddEntry(CCnumu_nuE_spec_var, ("%s"%detvar).replace("Att-BugFix","Attenuation"), "l")
  leg_CCnumu_nuE.Draw()
  cnv_spec_CCnumu_nuE.Write()
  cnv_covar_CCnumu_nuE = rt.TCanvas("cnv_covar_CCnumu_%s_nuE"%detvar,"cnv_covar_CCnumu_%s_nuE"%detvar)
  CCnumu_nuE.Draw("COLZ")
  cnv_covar_CCnumu_nuE.Write()
  cnv_spec_CCnumu_lepP = rt.TCanvas("cnv_spec_CCnumu_%s_lepP"%detvar,"cnv_spec_CCnumu_%s_lepP"%detvar)
  CCnumu_lepP_hists[0].Draw("EHIST")
  CCnumu_lepP_hists[1].Draw("EHISTSAME")
  leg_CCnumu_lepP = rt.TLegend(0.7,0.7,0.9,0.9)
  leg_CCnumu_lepP.AddEntry(CCnumu_lepP_spec_CV, "CV", "l")
  leg_CCnumu_lepP.AddEntry(CCnumu_lepP_spec_var, ("%s"%detvar).replace("Att-BugFix","Attenuation"), "l")
  leg_CCnumu_lepP.Draw()
  cnv_spec_CCnumu_lepP.Write()
  cnv_covar_CCnumu_lepP = rt.TCanvas("cnv_covar_CCnumu_%s_lepP"%detvar,"cnv_covar_CCnumu_%s_lepP"%detvar)
  CCnumu_lepP.Draw("COLZ")
  cnv_covar_CCnumu_lepP.Write()
  cnv_spec_CCnumu_lepCosTheta = rt.TCanvas("cnv_spec_CCnumu_%s_lepCosTheta"%detvar,"cnv_spec_CCnumu_%s_lepCosTheta"%detvar)
  CCnumu_lepCosTheta_hists[0].Draw("EHIST")
  CCnumu_lepCosTheta_hists[1].Draw("EHISTSAME")
  leg_CCnumu_lepCosTheta = rt.TLegend(0.7,0.7,0.9,0.9)
  leg_CCnumu_lepCosTheta.AddEntry(CCnumu_lepCosTheta_spec_CV, "CV", "l")
  leg_CCnumu_lepCosTheta.AddEntry(CCnumu_lepCosTheta_spec_var, ("%s"%detvar).replace("Att-BugFix","Attenuation"), "l")
  leg_CCnumu_lepCosTheta.Draw()
  cnv_spec_CCnumu_lepCosTheta.Write()
  cnv_covar_CCnumu_lepCosTheta = rt.TCanvas("cnv_covar_CCnumu_%s_lepCosTheta"%detvar,"cnv_covar_CCnumu_%s_lepCosTheta"%detvar)
  CCnumu_lepCosTheta.Draw("COLZ")
  cnv_covar_CCnumu_lepCosTheta.Write()


def getBins(tag):

  bins = []

  if "CCnue" in tag and ("nuE" in tag or "lepP" in tag) and not args.ccnueChi2:
    bins = [[0, 200], [200, 400], [400, 600], [600, 800], [800, 1000], [1000, 1200], [1200, 1400], [1400, 1600], [1600, 1800], [1800, 2000], [2000, 2200], [2200, 2400], [2400, 2600], [2600, 6000]]

  if "CCnue" in tag and "nuE" in tag and args.ccnueChi2:
    bins = [[0, 400], [400, 600], [600, 800], [800, 1000], [1000, 1200], [1200, 1400], [1400, 1600], [1600, 1800], [1800, 6000]]

  if "CCnue" in tag and "lepP" in tag and args.ccnueChi2:
    bins = [[0, 200], [200, 400], [400, 600], [600, 800], [800, 1000], [1000, 1200], [1200, 1400], [1400, 6000]]

  if "CCnumu" in tag and ("nuE" in tag or "lepP" in tag) and not args.ccnumuOverflow:
    bins = [[0, 100], [100, 200], [200, 300], [300, 400], [400, 500], [500, 600], [600, 700], [700, 800], [800, 900], [900, 1000], [1000, 1100], [1100, 1200], [1200, 1300], [1300, 1400], [1400, 1500], [1500, 1600], [1600, 1700], [1700, 1800], [1800, 1900], [1900, 2000], [2000, 6000]]

  if "CCnumu" in tag and "nuE" in tag and args.ccnumuOverflow:
    bins = [[0, 100], [100, 200], [200, 300], [300, 400], [400, 500], [500, 600], [600, 700], [700, 800], [800, 900], [900, 1000], [1000, 1100], [1100, 1200], [1200, 1300], [1300, 1400], [1400, 1500], [1500, 1600], [1600, 1700], [1700, 1800], [1800, 6000]]
    #bins = [[0, 100], [100, 200], [200, 300], [300, 400], [400, 500], [500, 600], [600, 700], [700, 800], [800, 900], [900, 1000], [1000, 1100], [1100, 1200], [1200, 1300], [1300, 1400], [1400, 6000]]

  if "CCnumu" in tag and "lepP" in tag and args.ccnumuOverflow:
    bins = [[0, 100], [100, 200], [200, 300], [300, 400], [400, 500], [500, 600], [600, 700], [700, 800], [800, 900], [900, 1000], [1000, 1100], [1100, 1200], [1200, 1300], [1300, 1400], [1400, 1500], [1500, 6000]]
    #bins = [[0, 100], [100, 200], [200, 300], [300, 400], [400, 500], [500, 600], [600, 700], [700, 800], [800, 900], [900, 1000], [1000, 1100], [1100, 1200], [1200, 6000]]

  if "CCnue" in tag and "CosTheta" in tag:
    if args.ccnueChi2:
      bins = [[-1.0, 0.375], [0.375, 0.5], [0.5, 0.625], [0.625, 0.75], [0.75, 0.875], [0.875, 1.0]]
    else:
      bins = [[-1.0, -0.875], [-0.875, -0.75], [-0.75, -0.625], [-0.625, -0.5], [-0.5, -0.375], [-0.375, -0.25], [-0.25, -0.125], [-0.125, 0.0], [0.0, 0.125], [0.125, 0.25], [0.25, 0.375], [0.375, 0.5], [0.5, 0.625], [0.625, 0.75], [0.75, 0.875], [0.875, 1.0]]

  if "CCnumu" in tag and "CosTheta" in tag:
    #if args.ccnumuOverflow:
    #  bins = [[-1.0, -0.625], [-0.625, -0.5], [-0.5, -0.375], [-0.375, -0.25], [-0.25, -0.125], [-0.125, 0.0], [0.0, 0.125], [0.125, 0.25], [0.25, 0.375], [0.375, 0.5], [0.5, 0.625], [0.625, 0.75], [0.75, 0.875], [0.875, 1.0]]
    #else:
    bins = [[-1.0, -0.875], [-0.875, -0.75], [-0.75, -0.625], [-0.625, -0.5], [-0.5, -0.375], [-0.375, -0.25], [-0.25, -0.125], [-0.125, 0.0], [0.0, 0.125], [0.125, 0.25], [0.25, 0.375], [0.375, 0.5], [0.5, 0.625], [0.625, 0.75], [0.75, 0.875], [0.875, 1.0]]

  return bins


print()
print("CCnue reco neutrino energy")
print()
print("{\setlength{\\tabcolsep}{0.5em} \\renewcommand{\\arraystretch}{1.2}")
print("\\begin{table}[ht]")
print("    \centering")
print("    \\resizebox{\\textwidth}{!}{")
print("    \\begin{tabular}{ |p{0.9in}||p{0.6in}|p{0.6in}|p{0.6in}|p{0.6in}|p{0.6in}|p{0.6in}|p{0.6in}|p{0.6in}|p{0.6in}||p{0.6in}| }")
print("    \hline")
titleRow = "    {\\footnotesize \Centering{Reco $E_\\nu$ Bin (MeV)}} & "
for dv in detvarsTable:
  titleRow += "\Centering{%s} & "%dv
titleRow += "\Centering{Total} \\\\"
print(titleRow)
print("    \hline\hline")
bins = getBins("CCnue_nuE")
for iB in range(len(bins)):
  row = "    \Centering{%i - %i} & "%(bins[iB][0],bins[iB][1])
  tot = 0.
  for iV, dv in enumerate(detvars):
    row += "\Centering{%.3f} & "%CCnue_nuE_errors[dv][iB]
    tot += CCnue_nuE_errors[dv][iB]**2
  row += "\Centering{%.3f} \\\\"%sqrt(tot)
  print(row)
  print("    \hline")
print("    \end{tabular}")
print("    }")
print("\end{table}")
print("}")

print()
print("CCnumu reco neutrino energy")
print()
print("{\setlength{\\tabcolsep}{0.5em} \\renewcommand{\\arraystretch}{1.2}")
print("\\begin{table}[ht]")
print("    \centering")
print("    \\resizebox{\\textwidth}{!}{")
print("    \\begin{tabular}{ |p{0.9in}||p{0.6in}|p{0.6in}|p{0.6in}|p{0.6in}|p{0.6in}|p{0.6in}|p{0.6in}|p{0.6in}|p{0.6in}||p{0.6in}| }")
print("    \hline")
titleRow = "    {\\footnotesize \Centering{Reco $E_\\nu$ Bin (MeV)}} & "
for dv in detvarsTable:
  titleRow += "\Centering{%s} & "%dv
titleRow += "\Centering{Total} \\\\"
print(titleRow)
print("    \hline\hline")
bins = getBins("CCnumu_nuE")
for iB in range(len(bins)):
  row = "    \Centering{%i - %i} & "%(bins[iB][0],bins[iB][1])
  tot = 0.
  for iV, dv in enumerate(detvars):
    row += "\Centering{%.3f} & "%CCnumu_nuE_errors[dv][iB]
    tot += CCnumu_nuE_errors[dv][iB]**2
  row += "\Centering{%.3f} \\\\"%sqrt(tot)
  print(row)
  print("    \hline")
print("    \end{tabular}")
print("    }")
print("\end{table}")
print("}")
print()

print()
print("CCnue reco electron momentum")
print()
print("{\setlength{\\tabcolsep}{0.5em} \\renewcommand{\\arraystretch}{1.2}")
print("\\begin{table}[ht]")
print("    \centering")
print("    \\resizebox{\\textwidth}{!}{")
print("    \\begin{tabular}{ |p{0.9in}||p{0.6in}|p{0.6in}|p{0.6in}|p{0.6in}|p{0.6in}|p{0.6in}|p{0.6in}|p{0.6in}|p{0.6in}||p{0.6in}| }")
print("    \hline")
titleRow = "    {\\footnotesize \Centering{Reco $p_{e^-}$ Bin (MeV/c)}} & "
for dv in detvarsTable:
  titleRow += "\Centering{%s} & "%dv
titleRow += "\Centering{Total} \\\\"
print(titleRow)
print("    \hline\hline")
bins = getBins("CCnue_lepP")
for iB in range(len(bins)):
  row = "    \Centering{%i - %i} & "%(bins[iB][0],bins[iB][1])
  tot = 0.
  for iV, dv in enumerate(detvars):
    row += "\Centering{%.3f} & "%CCnue_lepP_errors[dv][iB]
    tot += CCnue_lepP_errors[dv][iB]**2
  row += "\Centering{%.3f} \\\\"%sqrt(tot)
  print(row)
  print("    \hline")
print("    \end{tabular}")
print("    }")
print("\end{table}")
print("}")

print()
print("CCnumu reco muon momentum")
print()
print("{\setlength{\\tabcolsep}{0.5em} \\renewcommand{\\arraystretch}{1.2}")
print("\\begin{table}[ht]")
print("    \centering")
print("    \\resizebox{\\textwidth}{!}{")
print("    \\begin{tabular}{ |p{0.9in}||p{0.6in}|p{0.6in}|p{0.6in}|p{0.6in}|p{0.6in}|p{0.6in}|p{0.6in}|p{0.6in}|p{0.6in}||p{0.6in}| }")
print("    \hline")
titleRow = "    {\\footnotesize \Centering{Reco $p_\\mu$ Bin (MeV/c)}} & "
for dv in detvarsTable:
  titleRow += "\Centering{%s} & "%dv
titleRow += "\Centering{Total} \\\\"
print(titleRow)
print("    \hline\hline")
bins = getBins("CCnumu_lepP")
for iB in range(len(bins)):
  row = "    \Centering{%i - %i} & "%(bins[iB][0],bins[iB][1])
  tot = 0.
  for iV, dv in enumerate(detvars):
    row += "\Centering{%.3f} & "%CCnumu_lepP_errors[dv][iB]
    tot += CCnumu_lepP_errors[dv][iB]**2
  row += "\Centering{%.3f} \\\\"%sqrt(tot)
  print(row)
  print("    \hline")
print("    \end{tabular}")
print("    }")
print("\end{table}")
print("}")
print()

print()
print("CCnue reco electron cos(theta)")
print()
print("{\setlength{\\tabcolsep}{0.5em} \\renewcommand{\\arraystretch}{1.2}")
print("\\begin{table}[ht]")
print("    \centering")
print("    \\resizebox{\\textwidth}{!}{")
print("    \\begin{tabular}{ |p{0.9in}||p{0.6in}|p{0.6in}|p{0.6in}|p{0.6in}|p{0.6in}|p{0.6in}|p{0.6in}|p{0.6in}|p{0.6in}||p{0.6in}| }")
print("    \hline")
titleRow = "    {\\footnotesize \Centering{Reco cos($\\theta$) Bin}} & "
for dv in detvarsTable:
  titleRow += "\Centering{%s} & "%dv
titleRow += "\Centering{Total} \\\\"
print(titleRow)
print("    \hline\hline")
bins = getBins("CCnue_lepCosTheta")
for iB in range(len(bins)):
  row = "    \Centering{%i - %i} & "%(bins[iB][0],bins[iB][1])
  tot = 0.
  for iV, dv in enumerate(detvars):
    row += "\Centering{%.3f} & "%CCnue_lepCosTheta_errors[dv][iB]
    tot += CCnue_lepCosTheta_errors[dv][iB]**2
  row += "\Centering{%.3f} \\\\"%sqrt(tot)
  print(row)
  print("    \hline")
print("    \end{tabular}")
print("    }")
print("\end{table}")
print("}")

print()
print("CCnumu reco muon cos(theta)")
print()
print("{\setlength{\\tabcolsep}{0.5em} \\renewcommand{\\arraystretch}{1.2}")
print("\\begin{table}[ht]")
print("    \centering")
print("    \\resizebox{\\textwidth}{!}{")
print("    \\begin{tabular}{ |p{0.9in}||p{0.6in}|p{0.6in}|p{0.6in}|p{0.6in}|p{0.6in}|p{0.6in}|p{0.6in}|p{0.6in}|p{0.6in}||p{0.6in}| }")
print("    \hline")
titleRow = "    {\\footnotesize \Centering{Reco cos($\\theta$) Bin}} & "
for dv in detvarsTable:
  titleRow += "\Centering{%s} & "%dv
titleRow += "\Centering{Total} \\\\"
print(titleRow)
print("    \hline\hline")
bins = getBins("CCnumu_lepCosTheta")
for iB in range(len(bins)):
  row = "    \Centering{%i - %i} & "%(bins[iB][0],bins[iB][1])
  tot = 0.
  for iV, dv in enumerate(detvars):
    row += "\Centering{%.3f} & "%CCnumu_lepCosTheta_errors[dv][iB]
    tot += CCnumu_lepCosTheta_errors[dv][iB]**2
  row += "\Centering{%.3f} \\\\"%sqrt(tot)
  print(row)
  print("    \hline")
print("    \end{tabular}")
print("    }")
print("\end{table}")
print("}")
print()



calcList = [(CCnue_frac_covar_nuE, "CCnue_frac_covar_nuE", getBins("CCnue_nuE")), (CCnumu_frac_covar_nuE, "CCnumu_frac_covar_nuE", getBins("CCnumu_nuE")), (CCnue_frac_covar_lepP, "CCnue_frac_covar_lepP", getBins("CCnue_lepP")), (CCnumu_frac_covar_lepP, "CCnumu_frac_covar_lepP", getBins("CCnumu_lepP")), (CCnue_frac_covar_lepCosTheta, "CCnue_frac_covar_lepCosTheta", getBins("CCnue_lepCosTheta")), (CCnumu_frac_covar_lepCosTheta, "CCnumu_frac_covar_lepCosTheta", getBins("CCnumu_lepCosTheta"))]


for itm in calcList:
  tag = itm[1].replace("CCnue","CCnueInc").replace("CCnumu","CCnumuInc").replace("nuE","recoNuE").replace("lepCosTheta","cosTheta").replace("_frac_covar_","_")
  print("%s_detvar = {"%tag)
  hCovar_frac = itm[0]
  bins = itm[2]
  nbins = len(bins)
  if tag == "CCnueInc_recoNuE":
    print(nbins,(bins[1][0] - (bins[1][1] - bins[1][0]))/1000.,(bins[-1][0]+bins[1][1]-bins[1][0])/1000.)
    h_CCnue_nuE_totErr = rt.TH1F("h_CCnue_nuE_totErr","Total Detector Variation Fractional Uncertainty for CCnue Selection",nbins,(bins[1][0] - (bins[1][1] - bins[1][0]))/1000.,(bins[-1][0]+bins[1][1]-bins[1][0])/1000.)
    h_CCnue_nuE_totErr.GetXaxis().SetTitle("reconstructed neutrino energy (GeV)")
    h_CCnue_nuE_totErr.GetXaxis().CenterTitle(True)
    h_CCnue_nuE_totErr.SetLineWidth(2)
  if tag == "CCnumuInc_recoNuE":
    h_CCnumu_nuE_totErr = rt.TH1F("h_CCnumu_nuE_totErr","Total Detector Variation Fractional Uncertainty for CCnumu Selection",nbins,bins[0][0]/1000.,(bins[-1][0]+bins[0][1]-bins[0][0])/1000.)
    h_CCnumu_nuE_totErr.GetXaxis().SetTitle("reconstructed neutrino energy (GeV)")
    h_CCnumu_nuE_totErr.GetXaxis().CenterTitle(True)
    h_CCnumu_nuE_totErr.SetLineWidth(2)
  if tag == "CCnueInc_lepP":
    h_CCnue_lepP_totErr = rt.TH1F("h_CCnue_lepP_totErr","Total Detector Variation Fractional Uncertainty for CCnue Selection",nbins,bins[0][0]/1000.,(bins[-1][0]+bins[0][1]-bins[0][0])/1000.)
    h_CCnue_lepP_totErr.GetXaxis().SetTitle("reconstructed electron momentum (GeV/c)")
    h_CCnue_lepP_totErr.GetXaxis().CenterTitle(True)
    h_CCnue_lepP_totErr.SetLineWidth(2)
  if tag == "CCnumuInc_lepP":
    h_CCnumu_lepP_totErr = rt.TH1F("h_CCnumu_lepP_totErr","Total Detector Variation Fractional Uncertainty for CCnumu Selection",nbins,bins[0][0]/1000.,(bins[-1][0]+bins[0][1]-bins[0][0])/1000.)
    h_CCnumu_lepP_totErr.GetXaxis().SetTitle("reconstructed muon momentum (GeV/c)")
    h_CCnumu_lepP_totErr.GetXaxis().CenterTitle(True)
    h_CCnumu_lepP_totErr.SetLineWidth(2)
  if tag == "CCnueInc_cosTheta":
    h_CCnue_lepCosTheta_totErr = rt.TH1F("h_CCnue_lepCosTheta_totErr","Total Detector Variation Fractional Uncertainty for CCnue Selection",nbins,bins[0][0]/1000.,(bins[-1][0]+bins[0][1]-bins[0][0])/1000.)
    h_CCnue_lepCosTheta_totErr.GetXaxis().SetTitle("reconstructed electron cos(theta)")
    h_CCnue_lepCosTheta_totErr.GetXaxis().CenterTitle(True)
    h_CCnue_lepCosTheta_totErr.SetLineWidth(2)
  if tag == "CCnumuInc_cosTheta":
    h_CCnumu_lepCosTheta_totErr = rt.TH1F("h_CCnumu_lepCosTheta_totErr","Total Detector Variation Fractional Uncertainty for CCnumu Selection",nbins,bins[0][0]/1000.,(bins[-1][0]+bins[0][1]-bins[0][0])/1000.)
    h_CCnumu_lepCosTheta_totErr.GetXaxis().SetTitle("reconstructed muon cos(theta)")
    h_CCnumu_lepCosTheta_totErr.GetXaxis().CenterTitle(True)
    h_CCnumu_lepCosTheta_totErr.SetLineWidth(2)
  for i in range(1,nbins+1,1):
    uncert = sqrt(hCovar_frac.GetBinContent(i,i))
    if tag == "CCnueInc_recoNuE":
      h_CCnue_nuE_totErr.SetBinContent(i, uncert)
    if tag == "CCnumuInc_recoNuE":
      h_CCnumu_nuE_totErr.SetBinContent(i, uncert)
    if tag == "CCnueInc_lepP":
      h_CCnue_lepP_totErr.SetBinContent(i, uncert)
    if tag == "CCnumuInc_lepP":
      h_CCnumu_lepP_totErr.SetBinContent(i, uncert)
    if tag == "CCnueInc_cosTheta":
      h_CCnue_lepCosTheta_totErr.SetBinContent(i, uncert)
    if tag == "CCnumuInc_cosTheta":
      h_CCnumu_lepCosTheta_totErr.SetBinContent(i, uncert)
    if "cosTheta" in tag:
      if i < nbins:
        print("  %i:  (%f, %f, %f),"%(i,uncert,bins[i-1][0],bins[i-1][1]))
      else:
        print("  %i:  (%f, %f, %f)"%(i,uncert,bins[i-1][0],bins[i-1][1]))
    else:
      if i < nbins:
        print("  %i:  (%f, %i, %i),"%(i,uncert,bins[i-1][0],bins[i-1][1]))
      else:
        print("  %i:  (%f, %i, %i)"%(i,uncert,bins[i-1][0],bins[i-1][1]))
  print("}")
  if tag == "CCnueInc_recoNuE":
    outfile.cd()
    cnv_CCnue_nuE_totErr = rt.TCanvas("cnv_CCnue_nuE_totErr","cnv_CCnue_nuE_totErr")
    h_CCnue_nuE_totErr.Draw("HIST")
    label_CCnue_nuE = getOverflowLabel(h_CCnue_nuE_totErr)
    label_CCnue_nuE.Draw()
    cnv_CCnue_nuE_totErr.Write()
  if tag == "CCnumuInc_recoNuE":
    outfile.cd()
    cnv_CCnumu_nuE_totErr = rt.TCanvas("cnv_CCnumu_nuE_totErr","cnv_CCnumu_nuE_totErr")
    h_CCnumu_nuE_totErr.Draw("HIST")
    label_CCnumu_nuE = getOverflowLabel(h_CCnumu_nuE_totErr)
    label_CCnumu_nuE.Draw()
    cnv_CCnumu_nuE_totErr.Write()
  if tag == "CCnueInc_lepP":
    outfile.cd()
    cnv_CCnue_lepP_totErr = rt.TCanvas("cnv_CCnue_lepP_totErr","cnv_CCnue_lepP_totErr")
    h_CCnue_lepP_totErr.Draw("HIST")
    label_CCnue_lepP = getOverflowLabel(h_CCnue_lepP_totErr)
    label_CCnue_lepP.Draw()
    cnv_CCnue_lepP_totErr.Write()
  if tag == "CCnumuInc_lepP":
    outfile.cd()
    cnv_CCnumu_lepP_totErr = rt.TCanvas("cnv_CCnumu_lepP_totErr","cnv_CCnumu_lepP_totErr")
    h_CCnumu_lepP_totErr.Draw("HIST")
    label_CCnumu_lepP = getOverflowLabel(h_CCnumu_lepP_totErr)
    label_CCnumu_lepP.Draw()
    cnv_CCnumu_lepP_totErr.Write()
  if tag == "CCnueInc_cosTheta":
    outfile.cd()
    cnv_CCnue_lepCosTheta_totErr = rt.TCanvas("cnv_CCnue_lepCosTheta_totErr","cnv_CCnue_lepCosTheta_totErr")
    h_CCnue_lepCosTheta_totErr.Draw("HIST")
    cnv_CCnue_lepCosTheta_totErr.Write()
  if tag == "CCnumuInc_cosTheta":
    outfile.cd()
    cnv_CCnumu_lepCosTheta_totErr = rt.TCanvas("cnv_CCnumu_lepCosTheta_totErr","cnv_CCnumu_lepCosTheta_totErr")
    h_CCnumu_lepCosTheta_totErr.Draw("HIST")
    cnv_CCnumu_lepCosTheta_totErr.Write()




outfile.cd()
cnv_CCnue_nuE = rt.TCanvas("cnv_CCnue_nuE","cnv_CCnue_nuE")
#CCnue_frac_covar_nuE = histAbs(CCnue_frac_covar_nuE)
CCnue_frac_covar_nuE.Draw("COLZ")
cnv_CCnue_nuE.Write()
CCnue_frac_covar_nuE.Write()
cnv_CCnumu_nuE = rt.TCanvas("cnv_CCnumu_nuE","cnv_CCnumu_nuE")
#CCnumu_frac_covar_nuE = histAbs(CCnumu_frac_covar_nuE)
CCnumu_frac_covar_nuE.Draw("COLZ")
cnv_CCnumu_nuE.Write()
CCnumu_frac_covar_nuE.Write()
cnv_CCnue_lepP = rt.TCanvas("cnv_CCnue_lepP","cnv_CCnue_lepP")
#CCnue_frac_covar_lepP = histAbs(CCnue_frac_covar_lepP)
CCnue_frac_covar_lepP.Draw("COLZ")
cnv_CCnue_lepP.Write()
CCnue_frac_covar_lepP.Write()
cnv_CCnumu_lepP = rt.TCanvas("cnv_CCnumu_lepP","cnv_CCnumu_lepP")
#CCnumu_frac_covar_lepP = histAbs(CCnumu_frac_covar_lepP)
CCnumu_frac_covar_lepP.Draw("COLZ")
cnv_CCnumu_lepP.Write()
CCnumu_frac_covar_lepP.Write()
cnv_CCnue_lepCosTheta = rt.TCanvas("cnv_CCnue_lepCosTheta","cnv_CCnue_lepCosTheta")
#CCnue_frac_covar_lepCosTheta = histAbs(CCnue_frac_covar_lepCosTheta)
CCnue_frac_covar_lepCosTheta.Draw("COLZ")
cnv_CCnue_lepCosTheta.Write()
CCnue_frac_covar_lepCosTheta.Write()
cnv_CCnumu_lepCosTheta = rt.TCanvas("cnv_CCnumu_lepCosTheta","cnv_CCnumu_lepCosTheta")
#CCnumu_frac_covar_lepCosTheta = histAbs(CCnumu_frac_covar_lepCosTheta)
CCnumu_frac_covar_lepCosTheta.Draw("COLZ")
cnv_CCnumu_lepCosTheta.Write()
CCnumu_frac_covar_lepCosTheta.Write()



