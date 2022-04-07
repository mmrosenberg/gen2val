
import sys
import argparse
import ROOT as rt
from array import array
from helpers.plotting_functions import sortHists

parser = argparse.ArgumentParser("Make Lepton Reco Validation Plots")
parser.add_argument("-f", "--file", required=True, type=str, help="input lepton validation file")
parser.add_argument("-o", "--outdir", default="./", type=str, help="output file directory")
parser.add_argument("-v", "--versionTag", required=True, type=str, help="version tag for output files")
parser.add_argument("--nue", help="analyze nue events", action="store_true")
parser.add_argument("--numu", help="analyze numu events", action="store_true")
parser.add_argument("--save1d", help="save 1d dists as png files", action="store_true")
parser.add_argument("--save2d", help="save 2d dists as png files", action="store_true")
parser.add_argument("--prong_mult", help="break down start pos error dists by prong multiplicity", action="store_true")
args = parser.parse_args()

if not args.nue ^ args.numu:
  sys.exit("Need to specify --nue or --numu. Exiting...")

def getCDF(hist, title):
  x, y = array('f'), array('f')
  n = hist.GetNbinsX()
  integral = hist.Integral()
  for i in range(1,n+1):
    x.append(hist.GetBinLowEdge(i)+hist.GetBinWidth(i))
    y.append(hist.Integral(1,i)/integral)
  g = rt.TGraph(n,x,y)
  g.SetTitle(title)
  g.GetXaxis().SetTitle(hist.GetXaxis().GetTitle())
  g.GetYaxis().SetTitle("CDF")
  g.SetLineColor(rt.kBlue)
  g.SetLineWidth(2)
  g.SetMarkerStyle(20)
  g.SetMarkerColor(rt.kBlue)
  g.SetName(hist.GetName().replace("h_","g_")+"_CDF")
  return g

rt.gStyle.SetOptStat(0)

f = rt.TFile(args.file)
t = f.Get("NuIntTree")

lepton = "electron"
title = "CC nue events"
if args.numu:
  lepton = "muon"
  title = "CC numu events"

h_recoStatus = rt.TH1F("h_recoStatus",title,4,0,4)
h_recoStatus.GetXaxis().SetTitle("reco status")
#h_recoStatus.SetLineWidth(2)

h_vtxDist = rt.TH1F("h_vtxRes",title,50,0,3.1)
h_vtxDist.GetXaxis().SetTitle("vertex resolution (cm)")
#h_vtxDist.SetLineWidth(2)

nBinsC = 50
nBinsP = 70

h_trueLepCompletenessI = rt.TH1F("h_completeness",title,nBinsC,0.,1.0+1.0/(nBinsC-1))
h_trueLepCompletenessI.GetXaxis().SetTitle("completeness")
#h_trueLepCompletenessI.SetLineWidth(2)

h_recoLepPixIAnscPurity = rt.TH1F("h_purity",title,nBinsP,0.,1.0+1.0/(nBinsP-1))
h_recoLepPixIAnscPurity.GetXaxis().SetTitle("purity")
#h_recoLepPixIAnscPurity.SetLineWidth(2)

startPosHigh = 55.
lepEHigh = 4.5
if args.numu:
  startPosHigh = 12.
  lepEHigh = 2.

h_recoLepStartPosErr = rt.TH1F("h_startPosErr",title,60,0.,startPosHigh)
h_recoLepStartPosErr.GetXaxis().SetTitle(lepton+" start position error (cm)")
#h_recoLepStartPosErr.SetLineWidth(2)

if args.prong_mult:

  h_recoLepStartPosErr_mult1 = rt.TH1F("h_startPosErr_mult1",title,60,0.,startPosHigh)
  h_recoLepStartPosErr_mult1.GetXaxis().SetTitle(lepton+" start position error (cm)")
  h_recoLepStartPosErr_mult1.SetLineWidth(2)
  h_recoLepStartPosErr_mult1.SetLineColor(rt.kBlue)
  
  h_recoLepStartPosErr_mult2 = rt.TH1F("h_startPosErr_mult2",title,60,0.,startPosHigh)
  h_recoLepStartPosErr_mult2.GetXaxis().SetTitle(lepton+" start position error (cm)")
  h_recoLepStartPosErr_mult2.SetLineWidth(2)
  h_recoLepStartPosErr_mult2.SetLineColor(rt.kRed)
  
  h_recoLepStartPosErr_mult3p = rt.TH1F("h_startPosErr_mult3p",title,60,0.,startPosHigh)
  h_recoLepStartPosErr_mult3p.GetXaxis().SetTitle(lepton+" start position error (cm)")
  h_recoLepStartPosErr_mult3p.SetLineWidth(2)
  h_recoLepStartPosErr_mult3p.SetLineColor(rt.kGreen)

h_recoLepStartDirErr = rt.TH1F("h_startDirErr",title,70,0.,180.)
h_recoLepStartDirErr.GetXaxis().SetTitle(lepton+" start direction error (degrees)")
#h_recoLepStartDirErr.SetLineWidth(2)

h_purVcomp = rt.TH2F("h_purVcomp",title,60,0.,1.+1./59.,60,0.,1.+1./59.)
h_purVcomp.GetXaxis().SetTitle("completeness")
h_purVcomp.GetYaxis().SetTitle("purity")
h_dirVpos = rt.TH2F("h_dirVpos",title,60,0.,startPosHigh,70,0.,180.)
h_dirVpos.GetXaxis().SetTitle(lepton+" start position error (cm)")
h_dirVpos.GetYaxis().SetTitle(lepton+" start direction error (degrees)")
h_posVe = rt.TH2F("h_posVe",title,60,0.,lepEHigh,60,0.,startPosHigh)
h_posVe.GetXaxis().SetTitle(lepton+" true energy (GeV)")
h_posVe.GetYaxis().SetTitle(lepton+" start position error (cm)")
if args.numu:
  h_dirVposER = rt.TH2F("h_dirVposER",title,60,0.,55.,70,0.,180.)
  h_dirVposER.GetXaxis().SetTitle(lepton+" start position error (cm)")
  h_dirVposER.GetYaxis().SetTitle(lepton+" start direction error (cm)")

cnv = rt.TCanvas("cnv","cnv")
t.Draw("recoStatus >> h_recoStatus")
t.Draw("pow(pow(vtxDist_x,2)+pow(vtxDist_y,2)+pow(vtxDist_z,2),0.5) >> h_vtxRes","recoStatus == 0")
#t.Draw("vtxDist >> h_vtxRes","recoStatus == 0")
t.Draw("trueLepCompletenessI >> h_completeness","recoStatus == 0")
t.Draw("recoLepPixIAnscPurity >> h_purity","recoStatus == 0")
t.Draw("recoLepStartPosErr >> h_startPosErr","recoStatus == 0")
if args.prong_mult:
  t.Draw("recoLepStartPosErr >> h_startPosErr_mult1","recoStatus == 0 && vtxNProngs == 1")
  t.Draw("recoLepStartPosErr >> h_startPosErr_mult2","recoStatus == 0 && vtxNProngs == 2")
  t.Draw("recoLepStartPosErr >> h_startPosErr_mult3p","recoStatus == 0 && vtxNProngs > 2")
t.Draw("recoLepStartDirErr >> h_startDirErr","recoStatus == 0")
t.Draw("recoLepPixIAnscPurity:trueLepCompletenessI >> h_purVcomp","recoStatus == 0","COLZ")
t.Draw("recoLepStartDirErr:recoLepStartPosErr >> h_dirVpos","recoStatus == 0","COLZ")
t.Draw("recoLepStartPosErr:trueLepE >> h_posVe", "recoStatus == 0", "COLZ")
if args.numu:
  t.Draw("recoLepStartDirErr:recoLepStartPosErr >> h_dirVposER","recoStatus == 0","COLZ")

g_vtxDist_CDF = getCDF(h_vtxDist, title)
g_recoLepStartPosErr_CDF = getCDF(h_recoLepStartPosErr, title)

if args.save2d:
  h_purVcomp.Draw("COLZ")
  cnv.SaveAs(args.outdir+"primary_"+lepton+"_validation_plots_"+args.versionTag+"_purityVcompleteness.png")
  h_dirVpos.Draw("COLZ")
  cnv.SaveAs(args.outdir+"primary_"+lepton+"_validation_plots_"+args.versionTag+"_dirErrVposErr.png")
  h_posVe.Draw("COLZ")
  cnv.SaveAs(args.outdir+"primary_"+lepton+"_validation_plots_"+args.versionTag+"_posErrVlepE.png")
  if args.numu:
    h_dirVposER.Draw("COLZ")
    cnv.SaveAs(args.outdir+"primary_"+lepton+"_validation_plots_"+args.versionTag+"_dirErrVposErr-eRange.png")

if args.save1d:
  h_recoStatus.Draw("EHIST")
  cnv.SaveAs(args.outdir+"primary_"+lepton+"_validation_plots_"+args.versionTag+"_recoStatus.png")
  h_vtxDist.Draw("EHIST")
  cnv.SaveAs(args.outdir+"primary_"+lepton+"_validation_plots_"+args.versionTag+"_vtxRes.png")
  g_vtxDist_CDF.Draw("ALP")
  cnv.SaveAs(args.outdir+"primary_"+lepton+"_validation_plots_"+args.versionTag+"_vtxRes_CDF.png")
  h_trueLepCompletenessI.Draw("EHIST")
  cnv.SaveAs(args.outdir+"primary_"+lepton+"_validation_plots_"+args.versionTag+"_completeness.png")
  h_recoLepPixIAnscPurity.Draw("EHIST")
  cnv.SaveAs(args.outdir+"primary_"+lepton+"_validation_plots_"+args.versionTag+"_purity.png")
  h_recoLepStartPosErr.Draw("EHIST")
  cnv.SaveAs(args.outdir+"primary_"+lepton+"_validation_plots_"+args.versionTag+"_startPosErr.png")
  g_recoLepStartPosErr_CDF.Draw("ALP")
  cnv.SaveAs(args.outdir+"primary_"+lepton+"_validation_plots_"+args.versionTag+"_startPosErr_CDF.png")
  if args.prong_mult:
    h_startPosErr_mult = [h_recoLepStartPosErr_mult1, h_recoLepStartPosErr_mult2, h_recoLepStartPosErr_mult3p]
    h_startPosErr_mult = sortHists(h_startPosErr_mult)
    h_startPosErr_mult[0].Draw("EHIST")
    h_startPosErr_mult[1].Draw("EHISTSAME")
    h_startPosErr_mult[2].Draw("EHISTSAME")
    leg = rt.TLegend(0.7,0.7,0.9,0.9)
    leg.AddEntry(h_recoLepStartPosErr_mult1, "1 Prong Vertices", "l")
    leg.AddEntry(h_recoLepStartPosErr_mult2, "2 Prong Vertices", "l")
    leg.AddEntry(h_recoLepStartPosErr_mult3p, ">2 Prong Vertices", "l")
    leg.Draw()
    cnv.SaveAs(args.outdir+"primary_"+lepton+"_validation_plots_"+args.versionTag+"_startPosErr_byMultiplicity.png")
  h_recoLepStartDirErr.Draw("EHIST")
  cnv.SaveAs(args.outdir+"primary_"+lepton+"_validation_plots_"+args.versionTag+"_startDirErr.png")
  cnv.SetLogy()
  h_recoLepStartDirErr.Draw("EHIST")
  cnv.SaveAs(args.outdir+"primary_"+lepton+"_validation_plots_"+args.versionTag+"_startDirErr-logY.png")
  h_recoLepPixIAnscPurity.Draw("EHIST")
  cnv.SaveAs(args.outdir+"primary_"+lepton+"_validation_plots_"+args.versionTag+"_purity-logY.png")
  

outfile = rt.TFile(args.outdir+"primary_"+lepton+"_validation_plots_"+args.versionTag+".root","RECREATE")
h_recoStatus.Write()
h_vtxDist.Write()
g_vtxDist_CDF.Write()
h_trueLepCompletenessI.Write()
h_recoLepPixIAnscPurity.Write()
h_recoLepStartPosErr.Write()
g_recoLepStartPosErr_CDF.Write()
if args.prong_mult:
  h_recoLepStartPosErr_mult1.Write()
  h_recoLepStartPosErr_mult2.Write()
  h_recoLepStartPosErr_mult3p.Write()
h_recoLepStartDirErr.Write()
h_purVcomp.Write()
h_dirVpos.Write()
h_posVe.Write()
if args.numu:
  h_dirVposER.Write()

outfile.Close()

