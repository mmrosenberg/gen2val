
import sys
import argparse
import ROOT as rt

def sortHists(hlist):
  iMax = -1
  maxBin = -99.
  i = 0 
  for h in hlist:
    for b in range(h.GetNbinsX()):
      if h.GetBinContent(b+1) > maxBin:
        maxBin = h.GetBinContent(b+1)
        iMax = i 
    i = i+1 
  sortedList = [hlist[iMax]]
  for i in range(len(hlist)):
    if i != iMax:
      sortedList.append(hlist[i])
  return sortedList

parser = argparse.ArgumentParser("Compare Lepton Reco Versions")
parser.add_argument("-f1", "--file1", required=True, type=str, help="plot_lepton_reco_val.py output 1")
parser.add_argument("-f2", "--file2", required=True, type=str, help="plot_lepton_reco_val.py output 2")
parser.add_argument("-v1", "--version1", required=True, type=str, help="version tag for file 1")
parser.add_argument("-v2", "--version2", required=True, type=str, help="version tag for file 2")
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

h_recoStatus_1 = f1.Get("h_recoStatus")
h_vtxRes_1 = f1.Get("h_vtxRes")
h_completeness_1 = f1.Get("h_completeness")
h_purity_1 = f1.Get("h_purity")
h_startPosErr_1 = f1.Get("h_startPosErr")
h_startDirErr_1 = f1.Get("h_startDirErr")

h_recoStatus_2 = f2.Get("h_recoStatus")
h_vtxRes_2 = f2.Get("h_vtxRes")
h_completeness_2 = f2.Get("h_completeness")
h_purity_2 = f2.Get("h_purity")
h_startPosErr_2 = f2.Get("h_startPosErr")
h_startDirErr_2 = f2.Get("h_startDirErr")

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


leg = rt.TLegend(0.5,0.7,0.7,0.9)
leg.AddEntry(h_recoStatus_1, args.version1, "l")
leg.AddEntry(h_recoStatus_2, args.version2, "l")

cnv = rt.TCanvas("cnv","cnv")

h_recoStatus = [h_recoStatus_1, h_recoStatus_2]
h_recoStatus = sortHists(h_recoStatus)
h_recoStatus[0].Draw("EHIST")
h_recoStatus[1].Draw("EHISTSAME")
leg.Draw()
cnv.SaveAs(args.outdir+"primary_"+lepton+"_validation_comparison_"+args.version1+"_vs_"+args.version2+"_recoStatus.png")

h_vtxRes = [h_vtxRes_1, h_vtxRes_2]
h_vtxRes = sortHists(h_vtxRes)
h_vtxRes[0].Draw("EHIST")
h_vtxRes[1].Draw("EHISTSAME")
leg.Draw()
cnv.SaveAs(args.outdir+"primary_"+lepton+"_validation_comparison_"+args.version1+"_vs_"+args.version2+"_vtxRes.png")

h_completeness = [h_completeness_1, h_completeness_2]
h_completeness = sortHists(h_completeness)
h_completeness[0].Draw("EHIST")
h_completeness[1].Draw("EHISTSAME")
leg.Draw()
cnv.SaveAs(args.outdir+"primary_"+lepton+"_validation_comparison_"+args.version1+"_vs_"+args.version2+"_completeness.png")

h_purity = [h_purity_1, h_purity_2]
h_purity = sortHists(h_purity)
h_purity[0].Draw("EHIST")
h_purity[1].Draw("EHISTSAME")
leg.Draw()
cnv.SaveAs(args.outdir+"primary_"+lepton+"_validation_comparison_"+args.version1+"_vs_"+args.version2+"_purity.png")

h_startPosErr = [h_startPosErr_1, h_startPosErr_2]
h_startPosErr = sortHists(h_startPosErr)
h_startPosErr[0].Draw("EHIST")
h_startPosErr[1].Draw("EHISTSAME")
leg.Draw()
cnv.SaveAs(args.outdir+"primary_"+lepton+"_validation_comparison_"+args.version1+"_vs_"+args.version2+"_startPosErr.png")

h_startDirErr = [h_startDirErr_1, h_startDirErr_2]
h_startDirErr = sortHists(h_startDirErr)
h_startDirErr[0].Draw("EHIST")
h_startDirErr[1].Draw("EHISTSAME")
leg.Draw()
cnv.SaveAs(args.outdir+"primary_"+lepton+"_validation_comparison_"+args.version1+"_vs_"+args.version2+"_startDirErr.png")

cnv.SetLogy()

h_startDirErr[0].Draw("EHIST")
h_startDirErr[1].Draw("EHISTSAME")
leg.Draw()
cnv.SaveAs(args.outdir+"primary_"+lepton+"_validation_comparison_"+args.version1+"_vs_"+args.version2+"_startDirErr-logY.png")

h_purity[0].Draw("EHIST")
h_purity[1].Draw("EHISTSAME")
leg.Draw()
cnv.SaveAs(args.outdir+"primary_"+lepton+"_validation_comparison_"+args.version1+"_vs_"+args.version2+"_purity-logY.png")

