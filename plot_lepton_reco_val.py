
import sys
import argparse
import ROOT as rt

parser = argparse.ArgumentParser("Make Lepton Reco Validation Plots")
parser.add_argument("-f", "--file", required=True, type=str, help="input lepton validation file")
parser.add_argument("-o", "--outdir", default="./", type=str, help="output file directory")
parser.add_argument("-v", "--versionTag", required=True, type=str, help="version tag for output files")
parser.add_argument("--nue", help="analyze nue events", action="store_true")
parser.add_argument("--numu", help="analyze numu events", action="store_true")
parser.add_argument("--save1d", help="save 1d dists as png files", action="store_true")
parser.add_argument("--save2d", help="save 2d dists as png files", action="store_true")
args = parser.parse_args()

if not args.nue ^ args.numu:
  sys.exit("Need to specify --nue or --numu. Exiting...")

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
if args.numu:
  startPosHigh = 12.

h_recoLepStartPosErr = rt.TH1F("h_startPosErr",title,60,0.,startPosHigh)
h_recoLepStartPosErr.GetXaxis().SetTitle(lepton+" start position error (cm)")
#h_recoLepStartPosErr.SetLineWidth(2)

h_recoLepStartDirErr = rt.TH1F("h_startDirErr",title,70,0.,180.)
h_recoLepStartDirErr.GetXaxis().SetTitle(lepton+" start direction error (degrees)")
#h_recoLepStartDirErr.SetLineWidth(2)

h_purVcomp = rt.TH2F("h_purVcomp",title,60,0.,1.+1./59.,60,0.,1.+1./59.)
h_purVcomp.GetXaxis().SetTitle("completeness")
h_purVcomp.GetYaxis().SetTitle("purity")
h_dirVpos = rt.TH2F("h_dirVpos",title,60,0.,startPosHigh,70,0.,180.)
h_dirVpos.GetXaxis().SetTitle(lepton+" start position error (cm)")
h_dirVpos.GetYaxis().SetTitle(lepton+" start direction error (degrees)")
if args.numu:
  h_dirVposER = rt.TH2F("h_dirVposER",title,60,0.,55.,70,0.,180.)
  h_dirVposER.GetXaxis().SetTitle(lepton+" start position error (cm)")
  h_dirVposER.GetYaxis().SetTitle(lepton+" start direction error (cm)")

cnv = rt.TCanvas("cnv","cnv")
t.Draw("recoStatus >> h_recoStatus")
t.Draw("vtxDist >> h_vtxRes","recoStatus == 0")
t.Draw("trueLepCompletenessI >> h_completeness","recoStatus == 0")
t.Draw("recoLepPixIAnscPurity >> h_purity","recoStatus == 0")
t.Draw("recoLepStartPosErr >> h_startPosErr","recoStatus == 0")
t.Draw("recoLepStartDirErr >> h_startDirErr","recoStatus == 0")
t.Draw("recoLepPixIAnscPurity:trueLepCompletenessI >> h_purVcomp","recoStatus == 0","COLZ")
t.Draw("recoLepStartDirErr:recoLepStartPosErr >> h_dirVpos","recoStatus == 0","COLZ")
if args.numu:
  t.Draw("recoLepStartDirErr:recoLepStartPosErr >> h_dirVposER","recoStatus == 0","COLZ")

if args.save2d:
  h_purVcomp.Draw("COLZ")
  cnv.SaveAs(args.outdir+"primary_"+lepton+"_validation_plots_"+args.versionTag+"_purityVcompleteness.png")
  h_dirVpos.Draw("COLZ")
  cnv.SaveAs(args.outdir+"primary_"+lepton+"_validation_plots_"+args.versionTag+"_dirErrVposErr.png")
  if args.numu:
    h_dirVposER.Draw("COLZ")
    cnv.SaveAs(args.outdir+"primary_"+lepton+"_validation_plots_"+args.versionTag+"_dirErrVposErr-eRange.png")

if args.save1d:
  h_recoStatus.Draw("EHIST")
  cnv.SaveAs(args.outdir+"primary_"+lepton+"_validation_plots_"+args.versionTag+"_recoStatus.png")
  h_vtxDist.Draw("EHIST")
  cnv.SaveAs(args.outdir+"primary_"+lepton+"_validation_plots_"+args.versionTag+"_vtxRes.png")
  h_trueLepCompletenessI.Draw("EHIST")
  cnv.SaveAs(args.outdir+"primary_"+lepton+"_validation_plots_"+args.versionTag+"_completeness.png")
  h_recoLepPixIAnscPurity.Draw("EHIST")
  cnv.SaveAs(args.outdir+"primary_"+lepton+"_validation_plots_"+args.versionTag+"_purity.png")
  h_recoLepStartPosErr.Draw("EHIST")
  cnv.SaveAs(args.outdir+"primary_"+lepton+"_validation_plots_"+args.versionTag+"_startPosErr.png")
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
h_trueLepCompletenessI.Write()
h_recoLepPixIAnscPurity.Write()
h_recoLepStartPosErr.Write()
h_recoLepStartDirErr.Write()
h_purVcomp.Write()
h_dirVpos.Write()
if args.numu:
  h_dirVposER.Write()

outfile.Close()

