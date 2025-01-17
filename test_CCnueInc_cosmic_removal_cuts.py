
import sys, argparse
import numpy as np
import ROOT as rt

from math import isinf
from helpers.plotting_functions import sortHists
from helpers.larflowreco_ana_funcs import isFiducial, isFiducialWC, getDistance


parser = argparse.ArgumentParser("Test CCnue Inclusive Secondary Shower Cut")
parser.add_argument("-fnu", "--bnbnu_file", type=str, default="flat_ntuples/dlgen2_reco_v2me06_ntuple_v3_mcc9_v28_wctagger_bnboverlay_trimmed_CCnueInc_passed.root", help="bnb nu input file")
parser.add_argument("-fnue", "--bnbnue_file", type=str, default="flat_ntuples/dlgen2_reco_v2me06_ntuple_v3_mcc9_v28_wctagger_nueintrinsics_trimmed_CCnueInc_passed.root", help="bnb nu input file")
parser.add_argument("-fext", "--extbnb_file", type=str, default="flat_ntuples/dlgen2_reco_v2me06_ntuple_v3_mcc9_v29e_dl_run3_G1_extbnb_trimmed_CCnueInc_passed.root", help="extbnb input file")
parser.add_argument("-o", "--outfile", type=str, default="foo.root", help="output root file name")
parser.add_argument("--overwritePOT", help="overwrite POT for trimmed files", action="store_true")
parser.add_argument("-n", "--nShowersCut", type=int, default=0, help="minimum number of secondary showers")
parser.add_argument("-c", "--compCut", type=float, default=0., help="only consider showers above this completeness")
parser.add_argument("-q", "--chargeCut", type=float, default=0., help="only consider showers above this charge")
parser.add_argument("-qs", "--secondaryChargeCut", type=float, default=0., help="secondary/non-primary charge cut")
args = parser.parse_args()


#rt.gROOT.SetBatch(True)
rt.TH1.SetDefaultSumw2(rt.kTRUE)
rt.gStyle.SetOptStat(0)

fnu = rt.TFile(args.bnbnu_file)
tnu = fnu.Get("EventTree")
if args.overwritePOT:
  fnu_orig = rt.TFile("flat_ntuples/dlgen2_reco_v2me06_ntuple_v1_mcc9_v28_wctagger_bnboverlay_partial.root")
  tnuPOT = fnu_orig.Get("potTree")
else:
  tnuPOT = fnu.Get("potTree")

fnue = rt.TFile(args.bnbnue_file)
tnue = fnue.Get("EventTree")
tnuePOT = fnue.Get("potTree")

fext = rt.TFile(args.extbnb_file)
text = fext.Get("EventTree")


run3POT = 4.3e+19 + 1.701e+20 + 2.97e+19 + 1.524e+17
runs1to3POT = 6.67e+20
#targetPOT = 4.43e19
targetPOT = 4.4e+19
targetPOTstring = "4.4e+19"

tnuePOTsum = 0.
for i in range(tnuePOT.GetEntries()):
  tnuePOT.GetEntry(i)
  tnuePOTsum = tnuePOTsum + tnuePOT.totGoodPOT

tnuPOTsum = 0.
for i in range(tnuPOT.GetEntries()):
  tnuPOT.GetEntry(i)
  tnuPOTsum = tnuPOTsum + tnuPOT.totGoodPOT

if args.extbnb_file == "selection_output/prepare_selection_test_output/prepare_selection_test_reco_v2me05_gen2val_v22_extbnb_file.root":
  textPOTsum = 1.3298521464785359e+19
elif "dlgen2_reco_v2me06_ntuple_v1_mcc9_v29e_dl_run3_G1_extbnb_partial" in args.extbnb_file or args.extbnb_file == "flat_ntuples/dlgen2_reco_v2me06_ntuple_v2_mcc9_v29e_dl_run3_G1_extbnb_trimmed_CCnueInc_passed.root" or args.extbnb_file == "flat_ntuples/dlgen2_reco_v2me06_ntuple_v3_mcc9_v29e_dl_run3_G1_extbnb_trimmed_CCnueInc_passed.root":
  textPOTsum = 2.561872704628622e+19
else:
  sys.exit("POT for input extBNB file unknown")



def configureHists(h_CCnumu, h_NCnumu, h_CCnue, h_NCnue, h_ext):
  h_CCnumu.GetYaxis().SetTitle("events per "+targetPOTstring+" POT")
  h_CCnumu.SetLineColor(rt.kBlue)
  h_CCnumu.SetLineWidth(2)
  h_NCnumu.GetYaxis().SetTitle("events per "+targetPOTstring+" POT")
  h_NCnumu.SetLineColor(40)
  h_NCnumu.SetLineWidth(2)
  h_CCnue.GetYaxis().SetTitle("events per "+targetPOTstring+" POT")
  h_CCnue.SetLineColor(rt.kRed)
  h_CCnue.SetLineWidth(2)
  h_NCnue.GetYaxis().SetTitle("events per "+targetPOTstring+" POT")
  h_NCnue.SetLineColor(8)
  h_NCnue.SetLineWidth(2)
  h_ext.GetYaxis().SetTitle("events per "+targetPOTstring+" POT")
  h_ext.SetLineColor(rt.kBlack)
  h_ext.SetLineWidth(2)
  return h_CCnumu, h_NCnumu, h_CCnue, h_NCnue, h_ext


h_nSS_CCnumu = rt.TH1F("h_nSS_CCnumu","Number of Secondary Showers with Completeness > %.2f and Charge > %.2f"%(args.compCut, args.chargeCut),20,0,20)
h_nSS_NCnumu = rt.TH1F("h_nSS_NCnumu","Number of Secondary Showers with Completeness > %.2f and Charge > %.2f"%(args.compCut, args.chargeCut),20,0,20)
h_nSS_CCnue = rt.TH1F("h_nSS_CCnue","Number of Secondary Showers with Completeness > %.2f and Charge > %.2f"%(args.compCut, args.chargeCut),20,0,20)
h_nSS_NCnue = rt.TH1F("h_nSS_NCnue","Number of Secondary Showers with Completeness > %.2f and Charge > %.2f"%(args.compCut, args.chargeCut),20,0,20)
h_nSS_ext = rt.TH1F("h_nSS_ext","Number of Secondary Showers with Completeness > %.2f and Charge > %.2f"%(args.compCut, args.chargeCut),20,0,20)
h_nSS_CCnumu, h_nSS_NCnumu, h_nSS_CCnue, h_nSS_NCnue, h_nSS_ext = configureHists(h_nSS_CCnumu,
 h_nSS_NCnumu, h_nSS_CCnue, h_nSS_NCnue, h_nSS_ext)

h_nS_CCnumu = rt.TH1F("h_nS_CCnumu","Number of Showers with Completeness > %.2f and Charge > %.2f"%(args.compCut, args.chargeCut),30,0,30)
h_nS_NCnumu = rt.TH1F("h_nS_NCnumu","Number of Showers with Completeness > %.2f and Charge > %.2f"%(args.compCut, args.chargeCut),30,0,30)
h_nS_CCnue = rt.TH1F("h_nS_CCnue","Number of Showers with Completeness > %.2f and Charge > %.2f"%(args.compCut, args.chargeCut),30,0,30)
h_nS_NCnue = rt.TH1F("h_nS_NCnue","Number of Showers with Completeness > %.2f and Charge > %.2f"%(args.compCut, args.chargeCut),30,0,30)
h_nS_ext = rt.TH1F("h_nS_ext","Number of Showers with Completeness > %.2f and Charge > %.2f"%(args.compCut, args.chargeCut),30,0,30)
h_nS_CCnumu, h_nS_NCnumu, h_nS_CCnue, h_nS_NCnue, h_nS_ext = configureHists(h_nS_CCnumu,
 h_nS_NCnumu, h_nS_CCnue, h_nS_NCnue, h_nS_ext)

h_shwChg_CCnumu = rt.TH1F("h_shwChg_CCnumu","Total Charge of Showers with Completeness > %.2f and Charge > %.2f"%(args.compCut, args.chargeCut),800,0,800000)
h_shwChg_NCnumu = rt.TH1F("h_shwChg_NCnumu","Total Charge of Showers with Completeness > %.2f and Charge > %.2f"%(args.compCut, args.chargeCut),800,0,800000)
h_shwChg_CCnue = rt.TH1F("h_shwChg_CCnue","Total Charge of Showers with Completeness > %.2f and Charge > %.2f"%(args.compCut, args.chargeCut),800,0,800000)
h_shwChg_NCnue = rt.TH1F("h_shwChg_NCnue","Total Charge of Showers with Completeness > %.2f and Charge > %.2f"%(args.compCut, args.chargeCut),800,0,800000)
h_shwChg_ext = rt.TH1F("h_shwChg_ext","Total Charge of Showers with Completeness > %.2f and Charge > %.2f"%(args.compCut, args.chargeCut),800,0,800000)
h_shwChg_CCnumu, h_shwChg_NCnumu, h_shwChg_CCnue, h_shwChg_NCnue, h_shwChg_ext = configureHists(h_shwChg_CCnumu,
 h_shwChg_NCnumu, h_shwChg_CCnue, h_shwChg_NCnue, h_shwChg_ext)

h_secShwChg_CCnumu = rt.TH1F("h_secShwChg_CCnumu","Total Charge of Secondary Showers with Completeness > %.2f and Charge > %.2f"%(args.compCut, args.chargeCut),800,0,800000)
h_secShwChg_NCnumu = rt.TH1F("h_secShwChg_NCnumu","Total Charge of Secondary Showers with Completeness > %.2f and Charge > %.2f"%(args.compCut, args.chargeCut),800,0,800000)
h_secShwChg_CCnue = rt.TH1F("h_secShwChg_CCnue","Total Charge of Secondary Showers with Completeness > %.2f and Charge > %.2f"%(args.compCut, args.chargeCut),800,0,800000)
h_secShwChg_NCnue = rt.TH1F("h_secShwChg_NCnue","Total Charge of Secondary Showers with Completeness > %.2f and Charge > %.2f"%(args.compCut, args.chargeCut),800,0,800000)
h_secShwChg_ext = rt.TH1F("h_secShwChg_ext","Total Charge of Secondary Showers with Completeness > %.2f and Charge > %.2f"%(args.compCut, args.chargeCut),800,0,800000)
h_secShwChg_CCnumu, h_secShwChg_NCnumu, h_secShwChg_CCnue, h_secShwChg_NCnue, h_secShwChg_ext = configureHists(h_secShwChg_CCnumu,
 h_secShwChg_NCnumu, h_secShwChg_CCnue, h_secShwChg_NCnue, h_secShwChg_ext)

h_nonPrimShwChg_CCnumu = rt.TH1F("h_nonPrimShwChg_CCnumu","Total Charge of Non-Primary Showers with Completeness > %.2f and Charge > %.2f"%(args.compCut, args.chargeCut),800,0,800000)
h_nonPrimShwChg_NCnumu = rt.TH1F("h_nonPrimShwChg_NCnumu","Total Charge of Non-Primary Showers with Completeness > %.2f and Charge > %.2f"%(args.compCut, args.chargeCut),800,0,800000)
h_nonPrimShwChg_CCnue = rt.TH1F("h_nonPrimShwChg_CCnue","Total Charge of Non-Primary Showers with Completeness > %.2f and Charge > %.2f"%(args.compCut, args.chargeCut),800,0,800000)
h_nonPrimShwChg_NCnue = rt.TH1F("h_nonPrimShwChg_NCnue","Total Charge of Non-Primary Showers with Completeness > %.2f and Charge > %.2f"%(args.compCut, args.chargeCut),800,0,800000)
h_nonPrimShwChg_ext = rt.TH1F("h_nonPrimShwChg_ext","Total Charge of Non-Primary Showers with Completeness > %.2f and Charge > %.2f"%(args.compCut, args.chargeCut),800,0,800000)
h_nonPrimShwChg_CCnumu, h_nonPrimShwChg_NCnumu, h_nonPrimShwChg_CCnue, h_nonPrimShwChg_NCnue, h_nonPrimShwChg_ext = configureHists(h_nonPrimShwChg_CCnumu,
 h_nonPrimShwChg_NCnumu, h_nonPrimShwChg_CCnue, h_nonPrimShwChg_NCnue, h_nonPrimShwChg_ext)

h_pcaRatio_CCnumu = rt.TH1F("h_pcaRatio_CCnumu","PCA Eigenvalue Ratio",50,0,1)
h_pcaRatio_NCnumu = rt.TH1F("h_pcaRatio_NCnumu","PCA Eigenvalue Ratio",50,0,1)
h_pcaRatio_CCnue = rt.TH1F("h_pcaRatio_CCnue","PCA Eigenvalue Ratio",50,0,1)
h_pcaRatio_NCnue = rt.TH1F("h_pcaRatio_NCnue","PCA Eigenvalue Ratio",50,0,1)
h_pcaRatio_ext = rt.TH1F("h_pcaRatio_ext","PCA Eigenvalue Ratio",50,0,1)
h_pcaRatio_CCnumu, h_pcaRatio_NCnumu, h_pcaRatio_CCnue, h_pcaRatio_NCnue, h_pcaRatio_ext = configureHists(h_pcaRatio_CCnumu,
 h_pcaRatio_NCnumu, h_pcaRatio_CCnue, h_pcaRatio_NCnue, h_pcaRatio_ext)

h_pca0y_CCnumu = rt.TH1F("h_pca0y_CCnumu","PC Axis0 y",52,-1.02,1.02)
h_pca0y_NCnumu = rt.TH1F("h_pca0y_NCnumu","PC Axis0 y",52,-1.02,1.02)
h_pca0y_CCnue = rt.TH1F("h_pca0y_CCnue","PC Axis0 y",52,-1.02,1.02)
h_pca0y_NCnue = rt.TH1F("h_pca0y_NCnue","PC Axis0 y",52,-1.02,1.02)
h_pca0y_ext = rt.TH1F("h_pca0y_ext","PC Axis0 y",52,-1.02,1.02)
h_pca0y_CCnumu, h_pca0y_NCnumu, h_pca0y_CCnue, h_pca0y_NCnue, h_pca0y_ext = configureHists(h_pca0y_CCnumu,
 h_pca0y_NCnumu, h_pca0y_CCnue, h_pca0y_NCnue, h_pca0y_ext)

h_elShwCompSum_CCnumu = rt.TH1F("h_elShwCompSum_CCnumu","Electron Shower Completeness Sum",60,-1,6)
h_elShwCompSum_NCnumu = rt.TH1F("h_elShwCompSum_NCnumu","Electron Shower Completeness Sum",60,-1,6)
h_elShwCompSum_CCnue = rt.TH1F("h_elShwCompSum_CCnue","Electron Shower Completeness Sum",60,-1,6)
h_elShwCompSum_NCnue = rt.TH1F("h_elShwCompSum_NCnue","Electron Shower Completeness Sum",60,-1,6)
h_elShwCompSum_ext = rt.TH1F("h_elShwCompSum_ext","Electron Shower Completeness Sum",60,-1,6)
h_elShwCompSum_CCnumu, h_elShwCompSum_NCnumu, h_elShwCompSum_CCnue, h_elShwCompSum_NCnue, h_elShwCompSum_ext = configureHists(h_elShwCompSum_CCnumu,
 h_elShwCompSum_NCnumu, h_elShwCompSum_CCnue, h_elShwCompSum_NCnue, h_elShwCompSum_ext)

h_avgShwComp_CCnumu = rt.TH1F("h_avgShwComp_CCnumu","Average Shower Completeness",40,-1,1)
h_avgShwComp_NCnumu = rt.TH1F("h_avgShwComp_NCnumu","Average Shower Completeness",40,-1,1)
h_avgShwComp_CCnue = rt.TH1F("h_avgShwComp_CCnue","Average Shower Completeness",40,-1,1)
h_avgShwComp_NCnue = rt.TH1F("h_avgShwComp_NCnue","Average Shower Completeness",40,-1,1)
h_avgShwComp_ext = rt.TH1F("h_avgShwComp_ext","Average Shower Completeness",40,-1,1)
h_avgShwComp_CCnumu, h_avgShwComp_NCnumu, h_avgShwComp_CCnue, h_avgShwComp_NCnue, h_avgShwComp_ext = configureHists(h_avgShwComp_CCnumu,
 h_avgShwComp_NCnumu, h_avgShwComp_CCnue, h_avgShwComp_NCnue, h_avgShwComp_ext)

h_avgNonPrimShwComp_CCnumu = rt.TH1F("h_avgNonPrimShwComp_CCnumu","Average Shower Completeness Excluding Primary e-",40,-1,1)
h_avgNonPrimShwComp_NCnumu = rt.TH1F("h_avgNonPrimShwComp_NCnumu","Average Shower Completeness Excluding Primary e-",40,-1,1)
h_avgNonPrimShwComp_CCnue = rt.TH1F("h_avgNonPrimShwComp_CCnue","Average Shower Completeness Excluding Primary e-",40,-1,1)
h_avgNonPrimShwComp_NCnue = rt.TH1F("h_avgNonPrimShwComp_NCnue","Average Shower Completeness Excluding Primary e-",40,-1,1)
h_avgNonPrimShwComp_ext = rt.TH1F("h_avgNonPrimShwComp_ext","Average Shower Completeness Excluding Primary e-",40,-1,1)
h_avgNonPrimShwComp_CCnumu, h_avgNonPrimShwComp_NCnumu, h_avgNonPrimShwComp_CCnue, h_avgNonPrimShwComp_NCnue, h_avgNonPrimShwComp_ext = configureHists(h_avgNonPrimShwComp_CCnumu,
 h_avgNonPrimShwComp_NCnumu, h_avgNonPrimShwComp_CCnue, h_avgNonPrimShwComp_NCnue, h_avgNonPrimShwComp_ext)

h_maxEventGap_CCnumu = rt.TH1F("h_maxEventGap_CCnumu","Largest PCA Projection Gap (cm)",55,-10,100)
h_maxEventGap_NCnumu = rt.TH1F("h_maxEventGap_NCnumu","Largest PCA Projection Gap (cm)",55,-10,100)
h_maxEventGap_CCnue = rt.TH1F("h_maxEventGap_CCnue","Largest PCA Projection Gap (cm)",55,-10,100)
h_maxEventGap_NCnue = rt.TH1F("h_maxEventGap_NCnue","Largest PCA Projection Gap (cm)",55,-10,100)
h_maxEventGap_ext = rt.TH1F("h_maxEventGap_ext","Largest PCA Projection Gap (cm)",55,-10,100)
h_maxEventGap_CCnumu, h_maxEventGap_NCnumu, h_maxEventGap_CCnue, h_maxEventGap_NCnue, h_maxEventGap_ext = configureHists(h_maxEventGap_CCnumu,
 h_maxEventGap_NCnumu, h_maxEventGap_CCnue, h_maxEventGap_NCnue, h_maxEventGap_ext)

h_maxEventDist_CCnumu = rt.TH1F("h_maxEventDist_CCnumu","Largest Continuous PCA Projection Distance (cm)",92,-20,350)
h_maxEventDist_NCnumu = rt.TH1F("h_maxEventDist_NCnumu","Largest Continuous PCA Projection Distance (cm)",92,-20,350)
h_maxEventDist_CCnue = rt.TH1F("h_maxEventDist_CCnue","Largest Continuous PCA Projection Distance (cm)",92,-20,350)
h_maxEventDist_NCnue = rt.TH1F("h_maxEventDist_NCnue","Largest Continuous PCA Projection Distance (cm)",92,-20,350)
h_maxEventDist_ext = rt.TH1F("h_maxEventDist_ext","Largest Continuous PCA Projection Distance (cm)",92,-20,350)
h_maxEventDist_CCnumu, h_maxEventDist_NCnumu, h_maxEventDist_CCnue, h_maxEventDist_NCnue, h_maxEventDist_ext = configureHists(h_maxEventDist_CCnumu,
 h_maxEventDist_NCnumu, h_maxEventDist_CCnue, h_maxEventDist_NCnue, h_maxEventDist_ext)

h2d_gapVScomp_CCnue = rt.TH2F("h_gapVScomp_CCnue","PCA Projection Gap vs. Average Completeness",21,0,1.05,20,0,20)
h2d_gapVScomp_ext = rt.TH2F("h_gapVScomp_ext","PCA Projection Gap vs. Average Completeness",21,0,1.05,20,0,20)
h2d_gapVScomp_CCnue.GetXaxis().SetTitle("average shower completeness")
h2d_gapVScomp_CCnue.GetYaxis().SetTitle("largest PCA projection gap (cm)")
h2d_gapVScomp_ext.GetXaxis().SetTitle("average shower completeness")
h2d_gapVScomp_ext.GetYaxis().SetTitle("largest PCA projection gap (cm)")

maxPCAratio = -1.



def FillNuHistos(h_CCnumu, h_NCnumu, h_NCnue, val, weight, eventType):
  if eventType == 0:
    h_CCnumu.Fill(val, weight)
  if eventType == 1:
    h_NCnumu.Fill(val, weight)
  if eventType == 2:
    h_NCnue.Fill(val, weight)
  return h_CCnumu, h_NCnumu, h_NCnue


for i in range(tnu.GetEntries()):

  tnu.GetEntry(i)

  eventType = -1

  if abs(tnu.trueNuPDG) == 14:
    #CC numu
    if tnu.trueNuCCNC == 0:
      eventType = 0
    #NC numu
    else:
      eventType = 1
  #NC nue
  if abs(tnu.trueNuPDG) == 12 and tnu.trueNuCCNC == 1:
    eventType = 2

  nClassShowers = 0
  nShowers = 0
  nSecShowers = 0
  showerChargeSum = 0.
  secShowerChargeSum = 0.
  recoElCharge = -1.
  nonPrimShwCharge = 0.
  nElShowers = 0
  elShwCompSum = 0.
  avgShwComp = 0.
  recoElComp = 0.

  for iS in range(tnu.nShowers):
    nonPrimShwCharge += tnu.showerCharge[iS]
    if tnu.showerClassified[iS] == 1:
      nClassShowers += 1
      avgShwComp += tnu.showerComp[iS]
    if tnu.showerClassified[iS] == 1 and tnu.showerIsSecondary[iS] != 1 and tnu.showerPID[iS] == 11:
      if tnu.showerCharge[iS] > recoElCharge:
        recoElCharge = tnu.showerCharge[iS]
        recoElComp = tnu.showerComp[iS]
    if tnu.showerClassified[iS] == 1 and tnu.showerPID[iS] == 11:
      nElShowers += 1
      elShwCompSum += tnu.showerComp[iS]
    if tnu.showerCharge[iS] < args.chargeCut:
      continue
    if args.compCut > 0. and not (tnu.showerClassified[iS] == 1 and tnu.showerComp[iS] > args.compCut):
      continue
    nShowers += 1
    showerChargeSum += tnu.showerCharge[iS]
    if tnu.showerIsSecondary[iS] == 1:
      nSecShowers += 1
      secShowerChargeSum += tnu.showerCharge[iS]

  avgNonPrimShwComp = (avgShwComp - recoElComp)/(nClassShowers - 1) if (nClassShowers > 1) else 0.
  avgShwComp = avgShwComp/nClassShowers if (nClassShowers > 0) else -0.5
  nonPrimShwCharge -= recoElCharge
  if tnu.eventPCAxis0[1] >= 0.:
    nonPrimShwCharge = 0.
    showerChargeSum = 0.
    secShowerChargeSum = 0.

  pcaEVsum = tnu.eventPCEigenVals[0] + tnu.eventPCEigenVals[1] + tnu.eventPCEigenVals[2]
  pcaEVratio = tnu.eventPCEigenVals[0]/pcaEVsum if (pcaEVsum > 0. and nShowers >= args.nShowersCut and tnu.eventPCAxis0[1] < 0. and nonPrimShwCharge > args.secondaryChargeCut) else 0.
  if pcaEVratio > maxPCAratio:
    maxPCAratio = pcaEVratio

  maxEventGap = tnu.eventPCProjMaxGap[1]
  maxEventDist = tnu.eventPCProjMaxDist[1]
  if nElShowers < 2 or nonPrimShwCharge < args.secondaryChargeCut or tnu.eventPCAxis0[1] > 0. or pcaEVratio < 0.9:
    elShwCompSum = -1.
  if tnu.eventPCAxis0[1] >= 0. or nonPrimShwCharge < args.secondaryChargeCut or pcaEVratio < 0.9:
    maxEventGap = -9.
    maxEventDist = -19.
    avgShwComp = -1.
    avgNonPrimShwComp = -1.

  h_nSS_CCnumu, h_nSS_NCnumu, h_nSS_NCnue = FillNuHistos(h_nSS_CCnumu,
    h_nSS_NCnumu, h_nSS_NCnue, nSecShowers, tnu.xsecWeight, eventType)
  h_nS_CCnumu, h_nS_NCnumu, h_nS_NCnue = FillNuHistos(h_nS_CCnumu,
    h_nS_NCnumu, h_nS_NCnue, nShowers, tnu.xsecWeight, eventType)
  h_shwChg_CCnumu, h_shwChg_NCnumu, h_shwChg_NCnue = FillNuHistos(h_shwChg_CCnumu,
    h_shwChg_NCnumu, h_shwChg_NCnue, showerChargeSum, tnu.xsecWeight, eventType)
  h_secShwChg_CCnumu, h_secShwChg_NCnumu, h_secShwChg_NCnue = FillNuHistos(h_secShwChg_CCnumu,
    h_secShwChg_NCnumu, h_secShwChg_NCnue, secShowerChargeSum, tnu.xsecWeight, eventType)
  h_nonPrimShwChg_CCnumu, h_nonPrimShwChg_NCnumu, h_nonPrimShwChg_NCnue = FillNuHistos(h_nonPrimShwChg_CCnumu,
    h_nonPrimShwChg_NCnumu, h_nonPrimShwChg_NCnue, nonPrimShwCharge, tnu.xsecWeight, eventType)
  h_pcaRatio_CCnumu, h_pcaRatio_NCnumu, h_pcaRatio_NCnue = FillNuHistos(h_pcaRatio_CCnumu,
    h_pcaRatio_NCnumu, h_pcaRatio_NCnue, pcaEVratio, tnu.xsecWeight, eventType)
  h_pca0y_CCnumu, h_pca0y_NCnumu, h_pca0y_NCnue = FillNuHistos(h_pca0y_CCnumu,
    h_pca0y_NCnumu, h_pca0y_NCnue, tnu.eventPCAxis0[1], tnu.xsecWeight, eventType)
  h_elShwCompSum_CCnumu, h_elShwCompSum_NCnumu, h_elShwCompSum_NCnue = FillNuHistos(h_elShwCompSum_CCnumu,
    h_elShwCompSum_NCnumu, h_elShwCompSum_NCnue, elShwCompSum, tnu.xsecWeight, eventType)
  h_avgShwComp_CCnumu, h_avgShwComp_NCnumu, h_avgShwComp_NCnue = FillNuHistos(h_avgShwComp_CCnumu,
    h_avgShwComp_NCnumu, h_avgShwComp_NCnue, avgShwComp, tnu.xsecWeight, eventType)
  h_avgNonPrimShwComp_CCnumu, h_avgNonPrimShwComp_NCnumu, h_avgNonPrimShwComp_NCnue = FillNuHistos(h_avgNonPrimShwComp_CCnumu,
    h_avgNonPrimShwComp_NCnumu, h_avgNonPrimShwComp_NCnue, avgNonPrimShwComp, tnu.xsecWeight, eventType)
  h_maxEventGap_CCnumu, h_maxEventGap_NCnumu, h_maxEventGap_NCnue = FillNuHistos(h_maxEventGap_CCnumu,
    h_maxEventGap_NCnumu, h_maxEventGap_NCnue, maxEventGap, tnu.xsecWeight, eventType)
  h_maxEventDist_CCnumu, h_maxEventDist_NCnumu, h_maxEventDist_NCnue = FillNuHistos(h_maxEventDist_CCnumu,
    h_maxEventDist_NCnumu, h_maxEventDist_NCnue, maxEventDist, tnu.xsecWeight, eventType)



nSignalCut = 0.

for i in range(tnue.GetEntries()):

  tnue.GetEntry(i)

  nClassShowers = 0
  nShowers = 0
  nSecShowers = 0
  showerChargeSum = 0.
  secShowerChargeSum = 0.
  recoElCharge = -1.
  nonPrimShwCharge = 0.
  nElShowers = 0
  elShwCompSum = 0.
  avgShwComp = 0.
  recoElComp = 0.

  for iS in range(tnue.nShowers):
    nonPrimShwCharge += tnue.showerCharge[iS]
    if tnue.showerClassified[iS] == 1:
      nClassShowers += 1
      avgShwComp += tnue.showerComp[iS]
    if tnue.showerClassified[iS] == 1 and tnue.showerIsSecondary[iS] != 1 and tnue.showerPID[iS] == 11:
      if tnue.showerCharge[iS] > recoElCharge:
        recoElCharge = tnue.showerCharge[iS]
        recoElComp = tnue.showerComp[iS]
    if tnue.showerClassified[iS] == 1 and tnue.showerPID[iS] == 11:
      nElShowers += 1
      elShwCompSum += tnue.showerComp[iS]
    if tnue.showerCharge[iS] < args.chargeCut:
      continue
    if args.compCut > 0. and not (tnue.showerClassified[iS] == 1 and tnue.showerComp[iS] > args.compCut):
      continue
    nShowers += 1
    showerChargeSum += tnue.showerCharge[iS]
    if tnue.showerIsSecondary[iS] == 1:
      nSecShowers += 1
      secShowerChargeSum += tnue.showerCharge[iS]

  avgNonPrimShwComp = (avgShwComp - recoElComp)/(nClassShowers - 1) if (nClassShowers > 1) else 0.
  avgShwComp = avgShwComp/nClassShowers if (nClassShowers > 0) else -0.5
  nonPrimShwCharge -= recoElCharge
  if tnue.eventPCAxis0[1] >= 0.:
    nonPrimShwCharge = 0.
    showerChargeSum = 0.
    secShowerChargeSum = 0.

  pcaEVsum = tnue.eventPCEigenVals[0] + tnue.eventPCEigenVals[1] + tnue.eventPCEigenVals[2]
  pcaEVratio = tnue.eventPCEigenVals[0]/pcaEVsum if (pcaEVsum > 0. and nShowers >= args.nShowersCut and tnue.eventPCAxis0[1] < 0. and nonPrimShwCharge > args.secondaryChargeCut) else 0.
  if pcaEVratio > maxPCAratio:
    maxPCAratio = pcaEVratio

  maxEventGap = tnue.eventPCProjMaxGap[1]
  maxEventDist = tnue.eventPCProjMaxDist[1]
  if nElShowers < 2 or nonPrimShwCharge < args.secondaryChargeCut or tnue.eventPCAxis0[1] > 0. or pcaEVratio < 0.9:
    elShwCompSum = -1.
  if tnue.eventPCAxis0[1] >= 0. or nonPrimShwCharge < args.secondaryChargeCut or pcaEVratio < 0.9:
    maxEventGap = -9.
    maxEventDist = -19.
    avgShwComp = -1.
    avgNonPrimShwComp = -1.
  else:
    h2d_gapVScomp_CCnue.Fill(avgShwComp, maxEventGap, tnue.xsecWeight)

  if tnue.eventPCAxis0[1] < 0. and nonPrimShwCharge > args.secondaryChargeCut and pcaEVratio > 0.9 and avgShwComp < 0.5:
    nSignalCut += tnue.xsecWeight

  #if pcaEVratio > 0.98:
  #  print("example signal event removed by current cosmic cut: (fileid, run, subrun, event) = (%i, %i, %i, %i)"%(tnue.fileid, tnue.run, tnue.subrun, tnue.event))

  h_nSS_CCnue.Fill(nSecShowers, tnue.xsecWeight)
  h_nS_CCnue.Fill(nShowers, tnue.xsecWeight)
  h_shwChg_CCnue.Fill(showerChargeSum, tnue.xsecWeight)
  h_secShwChg_CCnue.Fill(secShowerChargeSum, tnue.xsecWeight)
  h_nonPrimShwChg_CCnue.Fill(nonPrimShwCharge, tnue.xsecWeight)
  h_pcaRatio_CCnue.Fill(pcaEVratio, tnue.xsecWeight)
  h_elShwCompSum_CCnue.Fill(elShwCompSum, tnue.xsecWeight)
  h_avgShwComp_CCnue.Fill(avgShwComp, tnue.xsecWeight)
  h_avgNonPrimShwComp_CCnue.Fill(avgNonPrimShwComp, tnue.xsecWeight)
  h_pca0y_CCnue.Fill(tnue.eventPCAxis0[1], tnue.xsecWeight)
  h_maxEventGap_CCnue.Fill(maxEventGap, tnue.xsecWeight)
  h_maxEventDist_CCnue.Fill(maxEventDist, tnue.xsecWeight)



for i in range(text.GetEntries()):

  text.GetEntry(i)

  nClassShowers = 0
  nShowers = 0
  nSecShowers = 0
  showerChargeSum = 0.
  secShowerChargeSum = 0.
  recoElCharge = -1.
  nonPrimShwCharge = 0.
  nElShowers = 0
  elShwCompSum = 0.
  avgShwComp = 0.
  recoElComp = 0.

  for iS in range(text.nShowers):
    nonPrimShwCharge += text.showerCharge[iS]
    if text.showerClassified[iS] == 1:
      nClassShowers += 1
      avgShwComp += text.showerComp[iS]
    if text.showerClassified[iS] == 1 and text.showerIsSecondary[iS] != 1 and text.showerPID[iS] == 11:
      if text.showerCharge[iS] > recoElCharge:
        recoElCharge = text.showerCharge[iS]
        recoElComp = text.showerComp[iS]
    if text.showerClassified[iS] == 1 and text.showerPID[iS] == 11:
      nElShowers += 1
      elShwCompSum += text.showerComp[iS]
    if text.showerCharge[iS] < args.chargeCut:
      continue
    if args.compCut > 0. and not (text.showerClassified[iS] == 1 and text.showerComp[iS] > args.compCut):
      continue
    nShowers += 1
    showerChargeSum += text.showerCharge[iS]
    if text.showerIsSecondary[iS] == 1:
      nSecShowers += 1
      secShowerChargeSum += text.showerCharge[iS]

  avgNonPrimShwComp = (avgShwComp - recoElComp)/(nClassShowers - 1) if (nClassShowers > 1) else 0.
  avgShwComp = avgShwComp/nClassShowers if (nClassShowers > 0) else -0.5
  nonPrimShwCharge -= recoElCharge
  if text.eventPCAxis0[1] >= 0.:
    nonPrimShwCharge = 0.
    showerChargeSum = 0.
    secShowerChargeSum = 0.

  pcaEVsum = text.eventPCEigenVals[0] + text.eventPCEigenVals[1] + text.eventPCEigenVals[2]
  pcaEVratio = text.eventPCEigenVals[0]/pcaEVsum if (pcaEVsum > 0. and nShowers >= args.nShowersCut and text.eventPCAxis0[1] < 0. and nonPrimShwCharge > args.secondaryChargeCut) else 0.
  if pcaEVratio > maxPCAratio:
    maxPCAratio = pcaEVratio

  maxEventGap = text.eventPCProjMaxGap[1]
  maxEventDist = text.eventPCProjMaxDist[1]
  if nElShowers < 2 or nonPrimShwCharge < args.secondaryChargeCut or text.eventPCAxis0[1] > 0. or pcaEVratio < 0.9:
    elShwCompSum = -1.
  if text.eventPCAxis0[1] >= 0. or nonPrimShwCharge < args.secondaryChargeCut or pcaEVratio < 0.9:
    maxEventGap = -9.
    maxEventDist = -19.
    avgShwComp = -1.
    avgNonPrimShwComp = -1.
  else:
    h2d_gapVScomp_ext.Fill(avgShwComp, maxEventGap)

  eventRemoved = False
  if text.eventPCAxis0[1] < 0. and nonPrimShwCharge > args.secondaryChargeCut and pcaEVratio > 0.9 and avgShwComp < 0.5:
    eventRemoved = True

  if eventRemoved:
    print("cosmic background event REMOVED (%i, %i, %i): nShowers = %i, nSecShowers = %i, showerCharge = %f, secShowerCharge = %f, nonPrimShowerCharge = %f"%(text.run, text.subrun, text.event, nShowers, nSecShowers, showerChargeSum, secShowerChargeSum, nonPrimShwCharge))
  else:
    print("cosmic background event RETAINED (%i, %i, %i): nShowers = %i, nSecShowers = %i, showerCharge = %f, secShowerCharge = %f, nonPrimShowerCharge = %f"%(text.run, text.subrun, text.event, nShowers, nSecShowers, showerChargeSum, secShowerChargeSum, nonPrimShwCharge))

  h_nSS_ext.Fill(nSecShowers)
  h_nS_ext.Fill(nShowers)
  h_shwChg_ext.Fill(showerChargeSum)
  h_secShwChg_ext.Fill(secShowerChargeSum)
  h_nonPrimShwChg_ext.Fill(nonPrimShwCharge)
  h_pcaRatio_ext.Fill(pcaEVratio)
  h_elShwCompSum_ext.Fill(elShwCompSum)
  h_avgShwComp_ext.Fill(avgShwComp)
  h_avgNonPrimShwComp_ext.Fill(avgNonPrimShwComp)
  h_pca0y_ext.Fill(text.eventPCAxis0[1])
  h_maxEventGap_ext.Fill(maxEventGap)
  h_maxEventDist_ext.Fill(maxEventDist)



nSignalCut *= targetPOT/tnuePOTsum

print("number of signal events that will be removed by cosmic cuts: ", nSignalCut)
print("max PCEigenVal ratio:", maxPCAratio)
print("extBNB scaling ratio:", targetPOT/textPOTsum)


h_nSS_CCnumu.Scale(targetPOT/tnuPOTsum)
h_nSS_NCnumu.Scale(targetPOT/tnuPOTsum)
h_nSS_CCnue.Scale(targetPOT/tnuePOTsum)
h_nSS_NCnue.Scale(targetPOT/tnuPOTsum)
h_nSS_ext.Scale(targetPOT/textPOTsum)

h_nS_CCnumu.Scale(targetPOT/tnuPOTsum)
h_nS_NCnumu.Scale(targetPOT/tnuPOTsum)
h_nS_CCnue.Scale(targetPOT/tnuePOTsum)
h_nS_NCnue.Scale(targetPOT/tnuPOTsum)
h_nS_ext.Scale(targetPOT/textPOTsum)

h_shwChg_CCnumu.Scale(targetPOT/tnuPOTsum)
h_shwChg_NCnumu.Scale(targetPOT/tnuPOTsum)
h_shwChg_CCnue.Scale(targetPOT/tnuePOTsum)
h_shwChg_NCnue.Scale(targetPOT/tnuPOTsum)
h_shwChg_ext.Scale(targetPOT/textPOTsum)

h_secShwChg_CCnumu.Scale(targetPOT/tnuPOTsum)
h_secShwChg_NCnumu.Scale(targetPOT/tnuPOTsum)
h_secShwChg_CCnue.Scale(targetPOT/tnuePOTsum)
h_secShwChg_NCnue.Scale(targetPOT/tnuPOTsum)
h_secShwChg_ext.Scale(targetPOT/textPOTsum)

h_nonPrimShwChg_CCnumu.Scale(targetPOT/tnuPOTsum)
h_nonPrimShwChg_NCnumu.Scale(targetPOT/tnuPOTsum)
h_nonPrimShwChg_CCnue.Scale(targetPOT/tnuePOTsum)
h_nonPrimShwChg_NCnue.Scale(targetPOT/tnuPOTsum)
h_nonPrimShwChg_ext.Scale(targetPOT/textPOTsum)

h_pcaRatio_CCnumu.Scale(targetPOT/tnuPOTsum)
h_pcaRatio_NCnumu.Scale(targetPOT/tnuPOTsum)
h_pcaRatio_CCnue.Scale(targetPOT/tnuePOTsum)
h_pcaRatio_NCnue.Scale(targetPOT/tnuPOTsum)
h_pcaRatio_ext.Scale(targetPOT/textPOTsum)

h_pca0y_CCnumu.Scale(targetPOT/tnuPOTsum)
h_pca0y_NCnumu.Scale(targetPOT/tnuPOTsum)
h_pca0y_CCnue.Scale(targetPOT/tnuePOTsum)
h_pca0y_NCnue.Scale(targetPOT/tnuPOTsum)
h_pca0y_ext.Scale(targetPOT/textPOTsum)

h_elShwCompSum_CCnumu.Scale(targetPOT/tnuPOTsum)
h_elShwCompSum_NCnumu.Scale(targetPOT/tnuPOTsum)
h_elShwCompSum_CCnue.Scale(targetPOT/tnuePOTsum)
h_elShwCompSum_NCnue.Scale(targetPOT/tnuPOTsum)
h_elShwCompSum_ext.Scale(targetPOT/textPOTsum)

h_avgShwComp_CCnumu.Scale(targetPOT/tnuPOTsum)
h_avgShwComp_NCnumu.Scale(targetPOT/tnuPOTsum)
h_avgShwComp_CCnue.Scale(targetPOT/tnuePOTsum)
h_avgShwComp_NCnue.Scale(targetPOT/tnuPOTsum)
h_avgShwComp_ext.Scale(targetPOT/textPOTsum)

h_avgNonPrimShwComp_CCnumu.Scale(targetPOT/tnuPOTsum)
h_avgNonPrimShwComp_NCnumu.Scale(targetPOT/tnuPOTsum)
h_avgNonPrimShwComp_CCnue.Scale(targetPOT/tnuePOTsum)
h_avgNonPrimShwComp_NCnue.Scale(targetPOT/tnuPOTsum)
h_avgNonPrimShwComp_ext.Scale(targetPOT/textPOTsum)

h_maxEventGap_CCnumu.Scale(targetPOT/tnuPOTsum)
h_maxEventGap_NCnumu.Scale(targetPOT/tnuPOTsum)
h_maxEventGap_CCnue.Scale(targetPOT/tnuePOTsum)
h_maxEventGap_NCnue.Scale(targetPOT/tnuPOTsum)
h_maxEventGap_ext.Scale(targetPOT/textPOTsum)

h_maxEventDist_CCnumu.Scale(targetPOT/tnuPOTsum)
h_maxEventDist_NCnumu.Scale(targetPOT/tnuPOTsum)
h_maxEventDist_CCnue.Scale(targetPOT/tnuePOTsum)
h_maxEventDist_NCnue.Scale(targetPOT/tnuPOTsum)
h_maxEventDist_ext.Scale(targetPOT/textPOTsum)

h2d_gapVScomp_CCnue.Scale(targetPOT/tnuePOTsum)
h2d_gapVScomp_ext.Scale(targetPOT/textPOTsum)


def configureLegend(leg, h_CCnumu, h_NCnumu, h_CCnue, h_NCnue, h_ext):
  leg.AddEntry(h_CCnue,"CC nue (%.2f)"%h_CCnue.Integral(),"l")
  leg.AddEntry(h_CCnumu,"CC numu (%.2f)"%h_CCnumu.Integral(),"l")
  leg.AddEntry(h_NCnue,"NC nue (%.2f)"%h_NCnue.Integral(),"l")
  leg.AddEntry(h_NCnumu,"NC numu (%.2f)"%h_NCnumu.Integral(),"l")
  leg.AddEntry(h_ext,"cosmic background (%.2f)"%h_ext.Integral(),"l")
  return leg
  
outFile = rt.TFile(args.outfile, "RECREATE")

cnv_nSS = rt.TCanvas("cnv_nSS","cnv_nSS")
hists_nSS = sortHists([h_nSS_CCnumu, h_nSS_NCnumu, h_nSS_CCnue, h_nSS_NCnue, h_nSS_ext])
hists_nSS[0].Draw("EHIST")
for i in range(1,len(hists_nSS)):
  hists_nSS[i].Draw("EHISTSAME")
leg_nSS = rt.TLegend(0.7,0.7,0.9,0.9)
leg_nSS = configureLegend(leg_nSS, h_nSS_CCnumu,
  h_nSS_NCnumu, h_nSS_CCnue, h_nSS_NCnue, h_nSS_ext)
leg_nSS.Draw()
cnv_nSS.Write()

cnv_nS = rt.TCanvas("cnv_nS","cnv_nS")
hists_nS = sortHists([h_nS_CCnumu, h_nS_NCnumu, h_nS_CCnue, h_nS_NCnue, h_nS_ext])
hists_nS[0].Draw("EHIST")
for i in range(1,len(hists_nS)):
  hists_nS[i].Draw("EHISTSAME")
leg_nS = rt.TLegend(0.7,0.7,0.9,0.9)
leg_nS = configureLegend(leg_nS, h_nS_CCnumu,
  h_nS_NCnumu, h_nS_CCnue, h_nS_NCnue, h_nS_ext)
leg_nS.Draw()
cnv_nS.Write()

cnv_shwChg = rt.TCanvas("cnv_shwChg","cnv_shwChg")
hists_shwChg = sortHists([h_shwChg_CCnumu, h_shwChg_NCnumu, h_shwChg_CCnue, h_shwChg_NCnue, h_shwChg_ext])
hists_shwChg[0].Draw("EHIST")
for i in range(1,len(hists_shwChg)):
  hists_shwChg[i].Draw("EHISTSAME")
leg_shwChg = rt.TLegend(0.7,0.7,0.9,0.9)
leg_shwChg = configureLegend(leg_shwChg, h_shwChg_CCnumu,
  h_shwChg_NCnumu, h_shwChg_CCnue, h_shwChg_NCnue, h_shwChg_ext)
leg_shwChg.Draw()
cnv_shwChg.Write()

cnv_secShwChg = rt.TCanvas("cnv_secShwChg","cnv_secShwChg")
hists_secShwChg = sortHists([h_secShwChg_CCnumu, h_secShwChg_NCnumu, h_secShwChg_CCnue, h_secShwChg_NCnue, h_secShwChg_ext])
hists_secShwChg[0].Draw("EHIST")
for i in range(1,len(hists_secShwChg)):
  hists_secShwChg[i].Draw("EHISTSAME")
leg_secShwChg = rt.TLegend(0.7,0.7,0.9,0.9)
leg_secShwChg = configureLegend(leg_secShwChg, h_secShwChg_CCnumu,
  h_secShwChg_NCnumu, h_secShwChg_CCnue, h_secShwChg_NCnue, h_secShwChg_ext)
leg_secShwChg.Draw()
cnv_secShwChg.Write()

cnv_nonPrimShwChg = rt.TCanvas("cnv_nonPrimShwChg","cnv_nonPrimShwChg")
hists_nonPrimShwChg = sortHists([h_nonPrimShwChg_CCnumu, h_nonPrimShwChg_NCnumu, h_nonPrimShwChg_CCnue, h_nonPrimShwChg_NCnue, h_nonPrimShwChg_ext])
hists_nonPrimShwChg[0].Draw("EHIST")
for i in range(1,len(hists_nonPrimShwChg)):
  hists_nonPrimShwChg[i].Draw("EHISTSAME")
leg_nonPrimShwChg = rt.TLegend(0.7,0.7,0.9,0.9)
leg_nonPrimShwChg = configureLegend(leg_nonPrimShwChg, h_nonPrimShwChg_CCnumu,
  h_nonPrimShwChg_NCnumu, h_nonPrimShwChg_CCnue, h_nonPrimShwChg_NCnue, h_nonPrimShwChg_ext)
leg_nonPrimShwChg.Draw()
cnv_nonPrimShwChg.Write()

cnv_pcaRatio = rt.TCanvas("cnv_pcaRatio","cnv_pcaRatio")
hists_pcaRatio = sortHists([h_pcaRatio_CCnumu, h_pcaRatio_NCnumu, h_pcaRatio_CCnue, h_pcaRatio_NCnue, h_pcaRatio_ext])
hists_pcaRatio[0].Draw("EHIST")
for i in range(1,len(hists_pcaRatio)):
  hists_pcaRatio[i].Draw("EHISTSAME")
leg_pcaRatio = rt.TLegend(0.7,0.7,0.9,0.9)
leg_pcaRatio = configureLegend(leg_pcaRatio, h_pcaRatio_CCnumu,
  h_pcaRatio_NCnumu, h_pcaRatio_CCnue, h_pcaRatio_NCnue, h_pcaRatio_ext)
leg_pcaRatio.Draw()
cnv_pcaRatio.Write()

cnv_pca0y = rt.TCanvas("cnv_pca0y","cnv_pca0y")
hists_pca0y = sortHists([h_pca0y_CCnumu, h_pca0y_NCnumu, h_pca0y_CCnue, h_pca0y_NCnue, h_pca0y_ext])
hists_pca0y[0].Draw("EHIST")
for i in range(1,len(hists_pca0y)):
  hists_pca0y[i].Draw("EHISTSAME")
leg_pca0y = rt.TLegend(0.7,0.7,0.9,0.9)
leg_pca0y = configureLegend(leg_pca0y, h_pca0y_CCnumu,
  h_pca0y_NCnumu, h_pca0y_CCnue, h_pca0y_NCnue, h_pca0y_ext)
leg_pca0y.Draw()
cnv_pca0y.Write()

cnv_elShwCompSum = rt.TCanvas("cnv_elShwCompSum","cnv_elShwCompSum")
hists_elShwCompSum = sortHists([h_elShwCompSum_CCnumu, h_elShwCompSum_NCnumu, h_elShwCompSum_CCnue, h_elShwCompSum_NCnue, h_elShwCompSum_ext])
hists_elShwCompSum[0].Draw("EHIST")
for i in range(1,len(hists_elShwCompSum)):
  hists_elShwCompSum[i].Draw("EHISTSAME")
leg_elShwCompSum = rt.TLegend(0.7,0.7,0.9,0.9)
leg_elShwCompSum = configureLegend(leg_elShwCompSum, h_elShwCompSum_CCnumu,
  h_elShwCompSum_NCnumu, h_elShwCompSum_CCnue, h_elShwCompSum_NCnue, h_elShwCompSum_ext)
leg_elShwCompSum.Draw()
cnv_elShwCompSum.Write()

cnv_avgShwComp = rt.TCanvas("cnv_avgShwComp","cnv_avgShwComp")
hists_avgShwComp = sortHists([h_avgShwComp_CCnumu, h_avgShwComp_NCnumu, h_avgShwComp_CCnue, h_avgShwComp_NCnue, h_avgShwComp_ext])
hists_avgShwComp[0].Draw("EHIST")
for i in range(1,len(hists_avgShwComp)):
  hists_avgShwComp[i].Draw("EHISTSAME")
leg_avgShwComp = rt.TLegend(0.7,0.7,0.9,0.9)
leg_avgShwComp = configureLegend(leg_avgShwComp, h_avgShwComp_CCnumu,
  h_avgShwComp_NCnumu, h_avgShwComp_CCnue, h_avgShwComp_NCnue, h_avgShwComp_ext)
leg_avgShwComp.Draw()
cnv_avgShwComp.Write()

cnv_avgNonPrimShwComp = rt.TCanvas("cnv_avgNonPrimShwComp","cnv_avgNonPrimShwComp")
hists_avgNonPrimShwComp = sortHists([h_avgNonPrimShwComp_CCnumu, h_avgNonPrimShwComp_NCnumu, h_avgNonPrimShwComp_CCnue, h_avgNonPrimShwComp_NCnue, h_avgNonPrimShwComp_ext])
hists_avgNonPrimShwComp[0].Draw("EHIST")
for i in range(1,len(hists_avgNonPrimShwComp)):
  hists_avgNonPrimShwComp[i].Draw("EHISTSAME")
leg_avgNonPrimShwComp = rt.TLegend(0.7,0.7,0.9,0.9)
leg_avgNonPrimShwComp = configureLegend(leg_avgNonPrimShwComp, h_avgNonPrimShwComp_CCnumu,
  h_avgNonPrimShwComp_NCnumu, h_avgNonPrimShwComp_CCnue, h_avgNonPrimShwComp_NCnue, h_avgNonPrimShwComp_ext)
leg_avgNonPrimShwComp.Draw()
cnv_avgNonPrimShwComp.Write()

cnv_maxEventGap = rt.TCanvas("cnv_maxEventGap","cnv_maxEventGap")
hists_maxEventGap = sortHists([h_maxEventGap_CCnumu, h_maxEventGap_NCnumu, h_maxEventGap_CCnue, h_maxEventGap_NCnue, h_maxEventGap_ext])
hists_maxEventGap[0].Draw("EHIST")
for i in range(1,len(hists_maxEventGap)):
  hists_maxEventGap[i].Draw("EHISTSAME")
leg_maxEventGap = rt.TLegend(0.7,0.7,0.9,0.9)
leg_maxEventGap = configureLegend(leg_maxEventGap, h_maxEventGap_CCnumu,
  h_maxEventGap_NCnumu, h_maxEventGap_CCnue, h_maxEventGap_NCnue, h_maxEventGap_ext)
leg_maxEventGap.Draw()
cnv_maxEventGap.Write()

cnv_maxEventDist = rt.TCanvas("cnv_maxEventDist","cnv_maxEventDist")
hists_maxEventDist = sortHists([h_maxEventDist_CCnumu, h_maxEventDist_NCnumu, h_maxEventDist_CCnue, h_maxEventDist_NCnue, h_maxEventDist_ext])
hists_maxEventDist[0].Draw("EHIST")
for i in range(1,len(hists_maxEventDist)):
  hists_maxEventDist[i].Draw("EHISTSAME")
leg_maxEventDist = rt.TLegend(0.7,0.7,0.9,0.9)
leg_maxEventDist = configureLegend(leg_maxEventDist, h_maxEventDist_CCnumu,
  h_maxEventDist_NCnumu, h_maxEventDist_CCnue, h_maxEventDist_NCnue, h_maxEventDist_ext)
leg_maxEventDist.Draw()
cnv_maxEventDist.Write()

cnv_gapVScomp_CCnue = rt.TCanvas("cnv_gapVScomp_CCnue","cnv_gapVScomp_CCnue")
h2d_gapVScomp_CCnue.Draw("COLZ")
cnv_gapVScomp_CCnue.Write()

cnv_gapVScomp_ext = rt.TCanvas("cnv_gapVScomp_ext","cnv_gapVScomp_ext")
h2d_gapVScomp_ext.Draw("COLZ")
cnv_gapVScomp_ext.Write()


