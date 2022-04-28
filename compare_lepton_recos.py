
import sys
import argparse
import ROOT as rt
from helpers.plotting_functions import sortHists

parser = argparse.ArgumentParser("Compare Lepton Reco Versions")
parser.add_argument("-f1", "--file1", required=True, type=str, help="plot_lepton_reco_val.py output 1")
parser.add_argument("-f2", "--file2", required=True, type=str, help="plot_lepton_reco_val.py output 2")
parser.add_argument("-f3", "--file3", default="none", type=str, help="plot_lepton_reco_val.py output 3")
parser.add_argument("-v1", "--version1", required=True, type=str, help="version tag for file 1")
parser.add_argument("-v2", "--version2", required=True, type=str, help="version tag for file 2")
parser.add_argument("-v3", "--version3", default="none", type=str, help="version tag for file 3")
parser.add_argument("-o", "--outdir", default="./", type=str, help="output file directory")
parser.add_argument("--nue", help="analyze nue events", action="store_true")
parser.add_argument("--numu", help="analyze numu events", action="store_true")
args = parser.parse_args()

if not args.nue ^ args.numu:
  sys.exit("Need to specify --nue or --numu. Exiting...")

lepton = "electron"
if args.numu:
  lepton = "muon"

rt.TH1.SetDefaultSumw2(rt.kTRUE)
rt.gStyle.SetOptStat(0)

f1 = rt.TFile(args.file1)
f2 = rt.TFile(args.file2)
if args.file3 != "none":
  f3 = rt.TFile(args.file3)

h_recoStatus_1 = f1.Get("h_recoStatus")
h_vtxRes_1 = f1.Get("h_vtxRes")
h_completeness_1 = f1.Get("h_completeness")
h_purity_1 = f1.Get("h_purity")
h_startPosErr_1 = f1.Get("h_startPosErr")
h_startDirErr_1 = f1.Get("h_startDirErr")
g_vtxRes_CDF_1 = f1.Get("g_vtxRes_CDF")
g_startPosErr_CDF_1 = f1.Get("g_startPosErr_CDF")
g_vtxRes_CDF_1.SetLineColor(rt.kBlue)
g_vtxRes_CDF_1.SetMarkerColor(rt.kBlue)
g_startPosErr_CDF_1.SetLineColor(rt.kBlue)
g_startPosErr_CDF_1.SetMarkerColor(rt.kBlue)

h_recoStatus_2 = f2.Get("h_recoStatus")
h_vtxRes_2 = f2.Get("h_vtxRes")
h_completeness_2 = f2.Get("h_completeness")
h_purity_2 = f2.Get("h_purity")
h_startPosErr_2 = f2.Get("h_startPosErr")
h_startDirErr_2 = f2.Get("h_startDirErr")
g_vtxRes_CDF_2 = f2.Get("g_vtxRes_CDF")
g_startPosErr_CDF_2 = f2.Get("g_startPosErr_CDF")
g_vtxRes_CDF_2.SetLineColor(rt.kRed)
g_vtxRes_CDF_2.SetMarkerColor(rt.kRed)
g_startPosErr_CDF_2.SetLineColor(rt.kRed)
g_startPosErr_CDF_2.SetMarkerColor(rt.kRed)

if args.file3 != "none":
  h_recoStatus_3 = f3.Get("h_recoStatus")
  h_vtxRes_3 = f3.Get("h_vtxRes")
  h_completeness_3 = f3.Get("h_completeness")
  h_purity_3 = f3.Get("h_purity")
  h_startPosErr_3 = f3.Get("h_startPosErr")
  h_startDirErr_3 = f3.Get("h_startDirErr")
  g_vtxRes_CDF_3 = f3.Get("g_vtxRes_CDF")
  g_startPosErr_CDF_3 = f3.Get("g_startPosErr_CDF")
  g_vtxRes_CDF_3.SetLineColor(8)
  g_vtxRes_CDF_3.SetMarkerColor(8)
  g_startPosErr_CDF_3.SetLineColor(8)
  g_startPosErr_CDF_3.SetMarkerColor(8)

h_recoStatus_1.Scale(1.0/h_recoStatus_1.Integral())
h_recoStatus_1.GetYaxis().SetTitle("area normalized event count")
h_recoStatus_1.SetLineWidth(2)
h_recoStatus_1.SetLineColor(rt.kBlue)
h_vtxRes_1.Scale(1.0/h_vtxRes_1.Integral())
h_vtxRes_1.GetYaxis().SetTitle("area normalized event count")
h_vtxRes_1.SetLineWidth(2)
h_vtxRes_1.SetLineColor(rt.kBlue)
h_completeness_1.Scale(1.0/h_completeness_1.Integral())
h_completeness_1.GetYaxis().SetTitle("area normalized event count")
h_completeness_1.SetLineWidth(2)
h_completeness_1.SetLineColor(rt.kBlue)
h_purity_1.Scale(1.0/h_purity_1.Integral())
h_purity_1.GetYaxis().SetTitle("area normalized event count")
h_purity_1.SetLineWidth(2)
h_purity_1.SetLineColor(rt.kBlue)
h_startPosErr_1.Scale(1.0/h_startPosErr_1.Integral())
h_startPosErr_1.GetYaxis().SetTitle("area normalized event count")
h_startPosErr_1.SetLineWidth(2)
h_startPosErr_1.SetLineColor(rt.kBlue)
h_startDirErr_1.Scale(1.0/h_startDirErr_1.Integral())
h_startDirErr_1.GetYaxis().SetTitle("area normalized event count")
h_startDirErr_1.SetLineWidth(2)
h_startDirErr_1.SetLineColor(rt.kBlue)

h_recoStatus_2.Scale(1.0/h_recoStatus_2.Integral())
h_recoStatus_2.GetYaxis().SetTitle("area normalized event count")
h_recoStatus_2.SetLineWidth(2)
h_recoStatus_2.SetLineColor(rt.kRed)
h_vtxRes_2.Scale(1.0/h_vtxRes_2.Integral())
h_vtxRes_2.GetYaxis().SetTitle("area normalized event count")
h_vtxRes_2.SetLineWidth(2)
h_vtxRes_2.SetLineColor(rt.kRed)
h_completeness_2.Scale(1.0/h_completeness_2.Integral())
h_completeness_2.GetYaxis().SetTitle("area normalized event count")
h_completeness_2.SetLineWidth(2)
h_completeness_2.SetLineColor(rt.kRed)
h_purity_2.Scale(1.0/h_purity_2.Integral())
h_purity_2.GetYaxis().SetTitle("area normalized event count")
h_purity_2.SetLineWidth(2)
h_purity_2.SetLineColor(rt.kRed)
h_startPosErr_2.Scale(1.0/h_startPosErr_2.Integral())
h_startPosErr_2.GetYaxis().SetTitle("area normalized event count")
h_startPosErr_2.SetLineWidth(2)
h_startPosErr_2.SetLineColor(rt.kRed)
h_startDirErr_2.Scale(1.0/h_startDirErr_2.Integral())
h_startDirErr_2.GetYaxis().SetTitle("area normalized event count")
h_startDirErr_2.SetLineWidth(2)
h_startDirErr_2.SetLineColor(rt.kRed)

if args.file3 != "none":
  h_recoStatus_3.Scale(1.0/h_recoStatus_3.Integral())
  h_recoStatus_3.GetYaxis().SetTitle("area normalized event count")
  h_recoStatus_3.SetLineWidth(2)
  h_recoStatus_3.SetLineColor(8)
  h_vtxRes_3.Scale(1.0/h_vtxRes_3.Integral())
  h_vtxRes_3.GetYaxis().SetTitle("area normalized event count")
  h_vtxRes_3.SetLineWidth(2)
  h_vtxRes_3.SetLineColor(8)
  h_completeness_3.Scale(1.0/h_completeness_3.Integral())
  h_completeness_3.GetYaxis().SetTitle("area normalized event count")
  h_completeness_3.SetLineWidth(2)
  h_completeness_3.SetLineColor(8)
  h_purity_3.Scale(1.0/h_purity_3.Integral())
  h_purity_3.GetYaxis().SetTitle("area normalized event count")
  h_purity_3.SetLineWidth(2)
  h_purity_3.SetLineColor(8)
  h_startPosErr_3.Scale(1.0/h_startPosErr_3.Integral())
  h_startPosErr_3.GetYaxis().SetTitle("area normalized event count")
  h_startPosErr_3.SetLineWidth(2)
  h_startPosErr_3.SetLineColor(8)
  h_startDirErr_3.Scale(1.0/h_startDirErr_3.Integral())
  h_startDirErr_3.GetYaxis().SetTitle("area normalized event count")
  h_startDirErr_3.SetLineWidth(2)
  h_startDirErr_3.SetLineColor(8)


outFileStem = args.outdir+"primary_"+lepton+"_validation_comparison_"+args.version1+"_vs_"+args.version2+"_"
if args.file3 != "none":
  outFileStem = outFileStem+"vs_"+args.version3+"_"

leg = rt.TLegend(0.5,0.7,0.7,0.9)
leg.AddEntry(h_recoStatus_1, args.version1, "l")
leg.AddEntry(h_recoStatus_2, args.version2, "l")
if args.file3 != "none":
  leg.AddEntry(h_recoStatus_3, args.version3, "l")

leg2 = rt.TLegend(0.6,0.4,0.8,0.6)
leg2.AddEntry(h_recoStatus_1, args.version1, "l")
leg2.AddEntry(h_recoStatus_2, args.version2, "l")
if args.file3 != "none":
  leg2.AddEntry(h_recoStatus_3, args.version3, "l")

leg3 = rt.TLegend(0.3,0.7,0.5,0.9)
leg3.AddEntry(h_recoStatus_1, args.version1, "l")
leg3.AddEntry(h_recoStatus_2, args.version2, "l")
if args.file3 != "none":
  leg3.AddEntry(h_recoStatus_3, args.version3, "l")

cnv = rt.TCanvas("cnv","cnv")

h_recoStatus = [h_recoStatus_1, h_recoStatus_2]
if args.file3 != "none":
  h_recoStatus.append(h_recoStatus_3)
h_recoStatus = sortHists(h_recoStatus)
h_recoStatus[0].Draw("EHIST")
h_recoStatus[1].Draw("EHISTSAME")
if args.file3 != "none":
  h_recoStatus[2].Draw("EHISTSAME")
leg.Draw()
cnv.SaveAs(outFileStem+"recoStatus.png")

h_vtxRes = [h_vtxRes_1, h_vtxRes_2]
if args.file3 != "none":
  h_vtxRes.append(h_vtxRes_3)
h_vtxRes = sortHists(h_vtxRes)
h_vtxRes[0].Draw("EHIST")
h_vtxRes[1].Draw("EHISTSAME")
if args.file3 != "none":
  h_vtxRes[2].Draw("EHISTSAME")
leg.Draw()
cnv.SaveAs(outFileStem+"vtxRes.png")

g_vtxRes_CDF = rt.TMultiGraph()
g_vtxRes_CDF.SetTitle(g_vtxRes_CDF_1.GetTitle()+"; "+g_vtxRes_CDF_1.GetXaxis().GetTitle()+"; "+g_vtxRes_CDF_1.GetYaxis().GetTitle())
g_vtxRes_CDF.Add(g_vtxRes_CDF_1)
g_vtxRes_CDF.Add(g_vtxRes_CDF_2)
if args.file3 != "none":
  g_vtxRes_CDF.Add(g_vtxRes_CDF_3)
g_vtxRes_CDF.Draw("ALP")
leg2.Draw()
cnv.SaveAs(outFileStem+"vtxResCDF.png")

h_completeness = [h_completeness_1, h_completeness_2]
if args.file3 != "none":
  h_completeness.append(h_completeness_3)
h_completeness = sortHists(h_completeness)
h_completeness[0].Draw("EHIST")
h_completeness[1].Draw("EHISTSAME")
if args.file3 != "none":
  h_completeness[2].Draw("EHISTSAME")
if args.numu:
  leg.Draw()
if args.nue:
  leg3.Draw()
cnv.SaveAs(outFileStem+"completeness.png")

h_purity = [h_purity_1, h_purity_2]
if args.file3 != "none":
  h_purity.append(h_purity_3)
h_purity = sortHists(h_purity)
h_purity[0].Draw("EHIST")
h_purity[1].Draw("EHISTSAME")
if args.file3 != "none":
  h_purity[2].Draw("EHISTSAME")
leg.Draw()
cnv.SaveAs(outFileStem+"purity.png")

h_startPosErr = [h_startPosErr_1, h_startPosErr_2]
if args.file3 != "none":
  h_startPosErr.append(h_startPosErr_3)
h_startPosErr = sortHists(h_startPosErr)
h_startPosErr[0].Draw("EHIST")
h_startPosErr[1].Draw("EHISTSAME")
if args.file3 != "none":
  h_startPosErr[2].Draw("EHISTSAME")
leg.Draw()
cnv.SaveAs(outFileStem+"startPosErr.png")

g_startPosErr_CDF = rt.TMultiGraph()
g_startPosErr_CDF.SetTitle(g_startPosErr_CDF_1.GetTitle()+"; "+g_startPosErr_CDF_1.GetXaxis().GetTitle()+"; "+g_startPosErr_CDF_1.GetYaxis().GetTitle())
g_startPosErr_CDF.Add(g_startPosErr_CDF_1)
g_startPosErr_CDF.Add(g_startPosErr_CDF_2)
if args.file3 != "none":
  g_startPosErr_CDF.Add(g_startPosErr_CDF_3)
g_startPosErr_CDF.Draw("ALP")
leg2.Draw()
cnv.SaveAs(outFileStem+"startPosErrCDF.png")

h_startDirErr = [h_startDirErr_1, h_startDirErr_2]
if args.file3 != "none":
  h_startDirErr.append(h_startDirErr_3)
h_startDirErr = sortHists(h_startDirErr)
h_startDirErr[0].Draw("EHIST")
h_startDirErr[1].Draw("EHISTSAME")
if args.file3 != "none":
  h_startDirErr[2].Draw("EHISTSAME")
leg.Draw()
cnv.SaveAs(outFileStem+"startDirErr.png")

cnv.SetLogy()

h_startDirErr[0].Draw("EHIST")
h_startDirErr[1].Draw("EHISTSAME")
if args.file3 != "none":
  h_startDirErr[2].Draw("EHISTSAME")
leg.Draw()
cnv.SaveAs(outFileStem+"startDirErr-logY.png")

h_purity[0].Draw("EHIST")
h_purity[1].Draw("EHISTSAME")
if args.file3 != "none":
  h_purity[2].Draw("EHISTSAME")
leg.Draw()
cnv.SaveAs(outFileStem+"purity-logY.png")

