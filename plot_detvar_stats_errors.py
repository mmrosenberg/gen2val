
import ROOT as rt

def getBinning(h):
  n = h.GetNbinsX()
  l = h.GetBinLowEdge(1)
  h = h.GetBinLowEdge(n)+h.GetBinWidth(n)
  return n, l, h

def getOverflowLabel(h):
  x = h.GetBinCenter(h.GetNbinsX()) - 0.1*h.GetBinWidth(h.GetNbinsX())
  rt.gPad.Update()
  y = -rt.gPad.GetFrame().GetY2()/27.
  label = rt.TText(x,y,"overflow")
  label.SetTextSize(0.028)
  label.SetTextAngle(-40)
  return label

infile = "systematics/detvar/calculate_detvar_systematics_LYAtt_final_origBinning_output.root"
outfile = "systematics/detvar/CCnumu_CV_stats_error_plots_final_origBinning_from_LYAtt.root"


f = rt.TFile(infile)

h_nuE = f.Get("h_CCnumu_recoNuE_CV")
h_lepP = f.Get("h_CCnumu_recoLepP_CV")
h_cosTheta = f.Get("h_CCnumu_recoLepCosTheta_CV")

nuE_n, nuE_l, nuE_h = getBinning(h_nuE)
lepP_n, lepP_l, lepP_h = getBinning(h_lepP)
cosTheta_n, cosTheta_l, cosTheta_h = getBinning(h_cosTheta)

h_nuE_err = rt.TH1F("h_nuE_err", "MC Fractional Stats Error for CCnumu CV Sample", nuE_n, nuE_l, nuE_h)
h_lepP_err = rt.TH1F("h_lepP_err", "MC Fractional Stats Error for CCnumu CV Sample", lepP_n, lepP_l, lepP_h)
h_cosTheta_err = rt.TH1F("h_cosTheta_err", "MC Fractional Stats Error for CCnumu CV Sample", cosTheta_n, cosTheta_l, cosTheta_h)

for i in range(1, nuE_n+1, 1):
  h_nuE_err.SetBinContent(i, h_nuE.GetBinError(i)/h_nuE.GetBinContent(i))
for i in range(1, lepP_n+1, 1):
  h_lepP_err.SetBinContent(i, h_lepP.GetBinError(i)/h_lepP.GetBinContent(i))
for i in range(1, cosTheta_n+1, 1):
  h_cosTheta_err.SetBinContent(i, h_cosTheta.GetBinError(i)/h_cosTheta.GetBinContent(i))


h_nuE_err.GetXaxis().SetTitle("reconstructed neutrino energy (GeV)")
h_nuE_err.GetXaxis().CenterTitle(True)
h_nuE_err.SetLineWidth(2)
h_lepP_err.GetXaxis().SetTitle("reconstructed muon momentum (GeV/c)")
h_lepP_err.GetXaxis().CenterTitle(True)
h_lepP_err.SetLineWidth(2)
h_cosTheta_err.GetXaxis().SetTitle("reconstructed muon cos(theta)")
h_cosTheta_err.SetLineWidth(2)

rt.gStyle.SetOptStat(0)
outFile = rt.TFile(outfile, "RECREATE")

cnv_nuE = rt.TCanvas("cnv_nuE","cnv_nuE")
h_nuE_err.Draw("HIST")
label_nuE = getOverflowLabel(h_nuE_err)
label_nuE.Draw()
cnv_nuE.Write()
h_nuE_err.Write()
cnv_lepP = rt.TCanvas("cnv_lepP","cnv_lepP")
h_lepP_err.Draw("HIST")
label_lepP = getOverflowLabel(h_lepP_err)
label_lepP.Draw()
cnv_lepP.Write()
h_lepP_err.Write()
cnv_cosTheta = rt.TCanvas("cnv_cosTheta","cnv_cosTheta")
h_cosTheta_err.Draw("HIST")
label_cosTheta = getOverflowLabel(h_cosTheta_err)
label_cosTheta.Draw()
cnv_cosTheta.Write()
h_cosTheta_err.Write()

