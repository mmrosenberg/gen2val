
import ROOT as rt
from helpers.plotting_functions import sortHists


rt.TH1.SetDefaultSumw2(rt.kTRUE)
rt.gStyle.SetOptStat(0)

detvars = ["LYAtt", "LYDown", "LYRayleigh", "recomb2", "SCE", "wiremodThetaXZ", "wiremodThetaYZ", "wiremodX", "wiremodYZ"]

outfile = rt.TFile("systematics/detvar/CCnumu_vertex_efficiency_test/plots.root","RECREATE")

pad_x1 = 0.005
pad_x2 = 0.995
uPad_y1 = 0.2525
uPad_y2 = 0.995
lPad_y1 = 0.005
lPad_y2 = 0.2475
bMargin = 0.3

hRat_TS = 0.13
hRat_xTS = 0.13
hRat_yTS = 0.12
hRat_xLS = 0.1
hRat_yLS = 0.1
hRat_xTO = 1.
hRat_yTO = 0.3
hRat_y1 = -0.05
hRat_y2 = 0.05
hRat_nD = 205

h_max = rt.TH1F("h_max","True CCnumu Events with Reconstructed Vertex",28,0.2,3.0)
h_min = rt.TH1F("h_min","True CCnumu Events with Reconstructed Vertex",28,0.2,3.0)
for i in range(1,29):
  h_max.SetBinContent(i,-9.)
  h_min.SetBinContent(i,9.)

hists_fracDiffs = []

for detvar in detvars:
  f = rt.TFile("systematics/detvar/CCnumu_vertex_efficiency_test/calculate_detvar_systematics_CCnumu_vertex_efficiency_%s_output.root"%detvar, "READ")
  h_CV = f.Get("h_trueNuE_CV")
  h_var = f.Get("h_trueNuE_var")
  h_CV.SetLineColor(rt.kBlue)
  h_CV.SetLineWidth(2)
  h_var.SetLineColor(rt.kRed)
  h_var.SetLineWidth(2)
  h_CV.GetXaxis().SetTitle("true neutrino energy (GeV)")
  h_CV.GetYaxis().SetTitle("events per 4.4e+19 POT")
  h_var.GetXaxis().SetTitle("true neutrino energy (GeV)")
  h_var.GetYaxis().SetTitle("events per 4.4e+19 POT")
  hists = sortHists([h_CV,h_var])
  h_diff = h_CV.Clone("h_diff_%s"%detvar)
  h_diff.Add(h_var, -1.)
  h_diff.Divide(h_CV)
  h_diff.SetTitle("")
  h_diff.GetYaxis().SetTitle("(CV - var)/CV")
  h_diff.GetXaxis().SetTitle("true neutrino energy (GeV)")
  h_diff_clone = h_diff.Clone("h_diff_%s_clone"%detvar)
  h_diff_clone.SetTitle("True CCnumu Events with Reconstructed Vertex")
  h_diff_clone.SetDirectory(0)
  hists_fracDiffs.append(h_diff_clone)
  for i in range(1,29):
    this_val = abs(h_diff.GetBinContent(i))
    if this_val > h_max.GetBinContent(i):
      h_max.SetBinContent(i, this_val)
    if this_val < h_min.GetBinContent(i):
      h_min.SetBinContent(i, this_val)
  h_diff.SetLineColor(rt.kBlack)
  h_diff.SetLineWidth(2)
  #h_diff.SetTitle("(CV - %s)/CV"%detvar)
  #h_diff.GetYaxis().SetTitle("")
  h_diff.SetTitleSize(hRat_TS)
  h_diff.GetXaxis().SetTitleSize(hRat_xTS)
  h_diff.GetYaxis().SetTitleSize(hRat_yTS)
  h_diff.GetXaxis().SetLabelSize(hRat_xLS)
  h_diff.GetYaxis().SetLabelSize(hRat_yLS)
  h_diff.GetXaxis().SetTitleOffset(hRat_xTO)
  h_diff.GetYaxis().SetTitleOffset(hRat_yTO)
  #h_diff.GetYaxis().SetRangeUser(hRat_y1,hRat_y2)
  h_diff.GetYaxis().SetNdivisions(hRat_nD)
  outfile.cd()
  cnv = rt.TCanvas("cnv_%s"%detvar,"cnv_%s"%detvar)
  uPad = rt.TPad("uPad_%s"%detvar,"uPad_%s"%detvar,pad_x1,uPad_y1,pad_x2,uPad_y2)
  lPad = rt.TPad("lPad_%s"%detvar,"lPad_%s"%detvar,pad_x1,lPad_y1,pad_x2,lPad_y2)
  lPad.SetBottomMargin(bMargin)
  uPad.Draw()
  lPad.Draw()
  lPad.SetGridy()
  uPad.cd()
  hists[0].Draw("EHIST")
  hists[1].Draw("EHISTSAME")
  leg = rt.TLegend(0.7,0.7,0.9,0.9)
  leg.AddEntry(h_CV, "CV", "l")
  leg.AddEntry(h_var, detvar, "l")
  leg.Draw()
  lPad.cd()
  h_diff.Draw("HIST")
  cnv.Write()

h_max.SetLineColor(13)
h_max.SetLineWidth(2)
h_max.GetXaxis().SetTitle("true neutrino energy (GeV)")
h_max.GetYaxis().SetTitle("|CV - var|/CV")
h_max.GetYaxis().SetTitleSize(0.05)
h_max.GetYaxis().SetRangeUser(0,0.05)
h_min.SetLineColor(8)
h_min.SetLineWidth(2)
outfile.cd()
cnv_bounds = rt.TCanvas("cnv_bounds","cnv_bounds")
h_max.Draw("HIST")
h_min.Draw("HISTSAME")
leg_bounds = rt.TLegend(0.7,0.7,0.9,0.9)
leg_bounds.AddEntry(h_max,"max in bin for any variation","l")
leg_bounds.AddEntry(h_min,"min in bin for any variation","l")
leg_bounds.Draw()
cnv_bounds.Write()


hists_fracDiffs = sortHists(hists_fracDiffs)
hists_fracDiffs_colors = [1,2,3,4,6,8,13,40,46]
for i in range(len(detvars)):
  hists_fracDiffs[i].SetLineColor(hists_fracDiffs_colors[i])
  hists_fracDiffs[i].SetLineWidth(2)

cnv_fracDiffs = rt.TCanvas("cnv_fracDiffs","cnv_fracDiffs")
cnv_fracDiffs.SetGridy()
hists_fracDiffs[0].GetYaxis().SetRangeUser(-0.045,0.085)
hists_fracDiffs[0].GetYaxis().SetTitleOffset(1.4)
hists_fracDiffs[0].GetYaxis().SetNdivisions(715)
hists_fracDiffs[0].Draw("HIST")
for i in range(1,len(hists_fracDiffs)):
  hists_fracDiffs[i].Draw("HISTSAME")
leg_fracDiffs = rt.TLegend(0.7,0.7,0.9,0.9)
for i in range(len(hists_fracDiffs)):
  var = hists_fracDiffs[i].GetName().replace("h_diff_","").replace("_clone","")
  leg_fracDiffs.AddEntry(hists_fracDiffs[i], var, "l")
leg_fracDiffs.Draw()
cnv_fracDiffs.Write()


