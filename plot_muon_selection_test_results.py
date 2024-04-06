
import sys, argparse
import numpy as np
import ROOT as rt

from math import isinf, sqrt
from helpers.plotting_functions import sortHists, getOverflowLabel
from helpers.larflowreco_ana_funcs import isFiducial, getDistance, getCosThetaBeamVector, getCosThetaGravVector 
from helpers.systematics import SetUncertainties


parser = argparse.ArgumentParser("Plot Selection Test Results")
parser.add_argument("-fnu", "--bnbnu_file", type=str, default="flat_ntuples/dlgen2_reco_v2me06_ntuple_v5_mcc9_v28_wctagger_bnboverlay.root", help="bnb nu input file")
parser.add_argument("-fnue", "--bnbnue_file", type=str, default="flat_ntuples/dlgen2_reco_v2me06_ntuple_v5_mcc9_v28_wctagger_nueintrinsics.root", help="bnb nu input file")
parser.add_argument("-fext", "--extbnb_file", type=str, default="flat_ntuples/dlgen2_reco_v2me06_ntuple_v5_mcc9_v29e_dl_run1_all_extbnb.root", help="extbnb input file")
parser.add_argument("-fdata", "--data_file", type=str, default="flat_ntuples/dlgen2_reco_v2me05_ntuple_v5_mcc9_v28_wctagger_bnb5e19.root", help="bnb data input file")
parser.add_argument("-o", "--outfile", type=str, default="selection_output/plot_selection_test_results_output/plot_muon_selection_test_results_output.root", help="output root file name")
parser.add_argument("-vsc", "--vertexScoreCut", type=float, default=0., help="minimum neutrino vertex score")
parser.add_argument("-cgc", "--cosThetaGravCut", type=float, default=-9., help="minimum muon cos(theta)")
parser.add_argument("-nsc", "--fromNeutralScoreCut", type=float, default=9., help="maximum muon from neutral parent score")
parser.add_argument("--procCuts", help="apply process classification/score cuts", action="store_true")
parser.add_argument("--smallFV", help="use 20cm fiducial volume", action="store_true")
parser.add_argument("--recoEOverflow", help="plot overflow bin for final recoE selection plot", action="store_true")
parser.add_argument("--plotPurityVsTrueE", help="plot purity vs. true E (cosmic bkg excluded if present)", action="store_true")
parser.add_argument("--write_ntuples", help="write ntuples with selected events to file", action="store_true")
args = parser.parse_args()

#rt.gROOT.SetBatch(True)
rt.TH1.SetDefaultSumw2(rt.kTRUE)
rt.gStyle.SetOptStat(0)

fnu = rt.TFile(args.bnbnu_file)
tnu = fnu.Get("EventTree")
tnuPOT = fnu.Get("potTree")

fnue = rt.TFile(args.bnbnue_file)
tnue = fnue.Get("EventTree")
tnuePOT = fnue.Get("potTree")

fext = rt.TFile(args.extbnb_file)
text = fext.Get("EventTree")

fdata = rt.TFile(args.data_file)
tdata = fdata.Get("EventTree")

if args.write_ntuples:
  fnu_trimmed = rt.TFile(args.bnbnu_file.replace(".root","_trimmed_CCnumuInc_passed.root"),"RECREATE")
  tnu_trimmed = tnu.CloneTree(0)
  fnue_trimmed = rt.TFile(args.bnbnue_file.replace(".root","_trimmed_CCnumuInc_passed.root"),"RECREATE")
  tnue_trimmed = tnue.CloneTree(0)
  fext_trimmed = rt.TFile(args.extbnb_file.replace(".root","_trimmed_CCnumuInc_passed.root"),"RECREATE")
  text_trimmed = text.CloneTree(0)
  fdata_trimmed = rt.TFile(args.data_file.replace(".root","_trimmed_CCnumuInc_passed.root"),"RECREATE")
  tdata_trimmed = tdata.CloneTree(0)

targetPOT = 4.4e+19
targetPOTstring = "4.4e+19"
BNBspills = 9764047.0

tnuePOTsum = 0.
for i in range(tnuePOT.GetEntries()):
  tnuePOT.GetEntry(i)
  tnuePOTsum = tnuePOTsum + tnuePOT.totGoodPOT

tnuPOTsum = 0.
for i in range(tnuPOT.GetEntries()):
  tnuPOT.GetEntry(i)
  tnuPOTsum = tnuPOTsum + tnuPOT.totGoodPOT

if "dlgen2_reco_v2me06_ntuple_v5_mcc9_v29e_dl_run3_G1_extbnb" in args.extbnb_file:
  #textPOTsum = 2.5631019528728764e+19
  textPOTsum = 1.6357377026304895e+20
elif "runs1to3" in args.extbnb_file:
  EXT = 34202767.0+38971237.0+465951.0+59572045.0+22166992.0+36721376.0+14817082.0+39195178.0+58677653.0+19214565.0+18619185.0
  textPOTsum = (EXT/BNBspills)*targetPOT
elif "run1_all" in args.extbnb_file:
  EXT = 34202767.0+38971237.0
  textPOTsum = (EXT/BNBspills)*targetPOT
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


h_cosFrac_CCnumu = rt.TH1F("h_cosFrac_CCnumu","Fraction of Hits On Cosmics",42,-1.05,1.05)
h_cosFrac_NCnumu = rt.TH1F("h_cosFrac_NCnumu","Fraction of Hits On Cosmics",42,-1.05,1.05)
h_cosFrac_CCnue = rt.TH1F("h_cosFrac_CCnue","Fraction of Hits On Cosmics",42,-1.05,1.05)
h_cosFrac_NCnue = rt.TH1F("h_cosFrac_NCnue","Fraction of Hits On Cosmics",42,-1.05,1.05)
h_cosFrac_ext = rt.TH1F("h_cosFrac_ext","Fraction of Hits On Cosmics",42,-1.05,1.05)
h_cosFrac_CCnumu, h_cosFrac_NCnumu, h_cosFrac_CCnue, h_cosFrac_NCnue, h_cosFrac_ext = configureHists(h_cosFrac_CCnumu, h_cosFrac_NCnumu, h_cosFrac_CCnue, h_cosFrac_NCnue, h_cosFrac_ext)

h_nMu_CCnumu = rt.TH1F("h_nMu_CCnumu","Number of Reco Muons",10,0,10)
h_nMu_NCnumu = rt.TH1F("h_nMu_NCnumu","Number of Reco Muons",10,0,10)
h_nMu_CCnue = rt.TH1F("h_nMu_CCnue","Number of Reco Muons",10,0,10)
h_nMu_NCnue = rt.TH1F("h_nMu_NCnue","Number of Reco Muons",10,0,10)
h_nMu_ext = rt.TH1F("h_nMu_ext","Number of Reco Muons",10,0,10)
h_nMu_CCnumu, h_nMu_NCnumu, h_nMu_CCnue, h_nMu_NCnue, h_nMu_ext = configureHists(h_nMu_CCnumu,
 h_nMu_NCnumu, h_nMu_CCnue, h_nMu_NCnue, h_nMu_ext)

h_vtxScore_CCnumu = rt.TH1F("h_vtxScore_CCnumu","Vertex Keypoint Score",51,0,1.02)
h_vtxScore_NCnumu = rt.TH1F("h_vtxScore_NCnumu","Vertex Keypoint Score",51,0,1.02)
h_vtxScore_CCnue = rt.TH1F("h_vtxScore_CCnue","Vertex Keypoint Score",51,0,1.02)
h_vtxScore_NCnue = rt.TH1F("h_vtxScore_NCnue","Vertex Keypoint Score",51,0,1.02)
h_vtxScore_ext = rt.TH1F("h_vtxScore_ext","Vertex Keypoint Score",51,0,1.02)
h_vtxScore_CCnumu, h_vtxScore_NCnumu, h_vtxScore_CCnue, h_vtxScore_NCnue, h_vtxScore_ext = configureHists(h_vtxScore_CCnumu, h_vtxScore_NCnumu, h_vtxScore_CCnue, h_vtxScore_NCnue, h_vtxScore_ext)

h_cosKPDist_CCnumu = rt.TH1F("h_cosKPDist_CCnumu","Distance to Nearest Cosmic KP (cm)",500,0,1000)
h_cosKPDist_NCnumu = rt.TH1F("h_cosKPDist_NCnumu","Distance to Nearest Cosmic KP (cm)",500,0,1000)
h_cosKPDist_CCnue = rt.TH1F("h_cosKPDist_CCnue","Distance to Nearest Cosmic KP (cm)",500,0,1000)
h_cosKPDist_NCnue = rt.TH1F("h_cosKPDist_NCnue","Distance to Nearest Cosmic KP (cm)",500,0,1000)
h_cosKPDist_ext = rt.TH1F("h_cosKPDist_ext","Distance to Nearest Cosmic KP (cm)",500,0,1000)
h_cosKPDist_CCnumu, h_cosKPDist_NCnumu, h_cosKPDist_CCnue, h_cosKPDist_NCnue, h_cosKPDist_ext = configureHists(h_cosKPDist_CCnumu, h_cosKPDist_NCnumu, h_cosKPDist_CCnue, h_cosKPDist_NCnue, h_cosKPDist_ext)

h_pcaEVRatio_CCnumu = rt.TH1F("h_pcaEVRatio_CCnumu","PCA Eigenvalue Ratio",50,0,1)
h_pcaEVRatio_NCnumu = rt.TH1F("h_pcaEVRatio_NCnumu","PCA Eigenvalue Ratio",50,0,1)
h_pcaEVRatio_CCnue = rt.TH1F("h_pcaEVRatio_CCnue","PCA Eigenvalue Ratio",50,0,1)
h_pcaEVRatio_NCnue = rt.TH1F("h_pcaEVRatio_NCnue","PCA Eigenvalue Ratio",50,0,1)
h_pcaEVRatio_ext = rt.TH1F("h_pcaEVRatio_ext","PCA Eigenvalue Ratio",50,0,1)
h_pcaEVRatio_CCnumu, h_pcaEVRatio_NCnumu, h_pcaEVRatio_CCnue, h_pcaEVRatio_NCnue, h_pcaEVRatio_ext = configureHists(h_pcaEVRatio_CCnumu, h_pcaEVRatio_NCnumu, h_pcaEVRatio_CCnue, h_pcaEVRatio_NCnue, h_pcaEVRatio_ext)

h_pcaCosThetaBeam_CCnumu = rt.TH1F("h_pcaCosThetaBeam_CCnumu","cos(angle between first PCA axis and beam)",42,-1.05,1.05)
h_pcaCosThetaBeam_NCnumu = rt.TH1F("h_pcaCosThetaBeam_NCnumu","cos(angle between first PCA axis and beam)",42,-1.05,1.05)
h_pcaCosThetaBeam_CCnue = rt.TH1F("h_pcaCosThetaBeam_CCnue","cos(angle between first PCA axis and beam)",42,-1.05,1.05)
h_pcaCosThetaBeam_NCnue = rt.TH1F("h_pcaCosThetaBeam_NCnue","cos(angle between first PCA axis and beam)",42,-1.05,1.05)
h_pcaCosThetaBeam_ext = rt.TH1F("h_pcaCosThetaBeam_ext","cos(angle between first PCA axis and beam)",42,-1.05,1.05)
h_pcaCosThetaBeam_CCnumu, h_pcaCosThetaBeam_NCnumu, h_pcaCosThetaBeam_CCnue, h_pcaCosThetaBeam_NCnue, h_pcaCosThetaBeam_ext = configureHists(h_pcaCosThetaBeam_CCnumu, h_pcaCosThetaBeam_NCnumu, h_pcaCosThetaBeam_CCnue, h_pcaCosThetaBeam_NCnue, h_pcaCosThetaBeam_ext)

h_pcaCosThetaGrav_CCnumu = rt.TH1F("h_pcaCosThetaGrav_CCnumu","cos(angle between first PCA axis and gravity)",42,-1.05,1.05)
h_pcaCosThetaGrav_NCnumu = rt.TH1F("h_pcaCosThetaGrav_NCnumu","cos(angle between first PCA axis and gravity)",42,-1.05,1.05)
h_pcaCosThetaGrav_CCnue = rt.TH1F("h_pcaCosThetaGrav_CCnue","cos(angle between first PCA axis and gravity)",42,-1.05,1.05)
h_pcaCosThetaGrav_NCnue = rt.TH1F("h_pcaCosThetaGrav_NCnue","cos(angle between first PCA axis and gravity)",42,-1.05,1.05)
h_pcaCosThetaGrav_ext = rt.TH1F("h_pcaCosThetaGrav_ext","cos(angle between first PCA axis and gravity)",42,-1.05,1.05)
h_pcaCosThetaGrav_CCnumu, h_pcaCosThetaGrav_NCnumu, h_pcaCosThetaGrav_CCnue, h_pcaCosThetaGrav_NCnue, h_pcaCosThetaGrav_ext = configureHists(h_pcaCosThetaGrav_CCnumu, h_pcaCosThetaGrav_NCnumu, h_pcaCosThetaGrav_CCnue, h_pcaCosThetaGrav_NCnue, h_pcaCosThetaGrav_ext)

h_selMu_cosThetaBeam_CCnumu = rt.TH1F("h_selMu_cosThetaBeam_CCnumu","cos(angle between reco primary muon and beam)",42,-1.05,1.05)
h_selMu_cosThetaBeam_NCnumu = rt.TH1F("h_selMu_cosThetaBeam_NCnumu","cos(angle between reco primary muon and beam)",42,-1.05,1.05)
h_selMu_cosThetaBeam_CCnue = rt.TH1F("h_selMu_cosThetaBeam_CCnue","cos(angle between reco primary muon and beam)",42,-1.05,1.05)
h_selMu_cosThetaBeam_NCnue = rt.TH1F("h_selMu_cosThetaBeam_NCnue","cos(angle between reco primary muon and beam)",42,-1.05,1.05)
h_selMu_cosThetaBeam_ext = rt.TH1F("h_selMu_cosThetaBeam_ext","cos(angle between reco primary muon and beam)",42,-1.05,1.05)
h_selMu_cosThetaBeam_CCnumu, h_selMu_cosThetaBeam_NCnumu, h_selMu_cosThetaBeam_CCnue, h_selMu_cosThetaBeam_NCnue, h_selMu_cosThetaBeam_ext = configureHists(h_selMu_cosThetaBeam_CCnumu, h_selMu_cosThetaBeam_NCnumu, h_selMu_cosThetaBeam_CCnue, h_selMu_cosThetaBeam_NCnue, h_selMu_cosThetaBeam_ext)

h_selMu_cosThetaGrav_CCnumu = rt.TH1F("h_selMu_cosThetaGrav_CCnumu","cos(angle between reco primary muon and gravity)",42,-1.05,1.05)
h_selMu_cosThetaGrav_NCnumu = rt.TH1F("h_selMu_cosThetaGrav_NCnumu","cos(angle between reco primary muon and gravity)",42,-1.05,1.05)
h_selMu_cosThetaGrav_CCnue = rt.TH1F("h_selMu_cosThetaGrav_CCnue","cos(angle between reco primary muon and gravity)",42,-1.05,1.05)
h_selMu_cosThetaGrav_NCnue = rt.TH1F("h_selMu_cosThetaGrav_NCnue","cos(angle between reco primary muon and gravity)",42,-1.05,1.05)
h_selMu_cosThetaGrav_ext = rt.TH1F("h_selMu_cosThetaGrav_ext","cos(angle between reco primary muon and gravity)",42,-1.05,1.05)
h_selMu_cosThetaGrav_CCnumu, h_selMu_cosThetaGrav_NCnumu, h_selMu_cosThetaGrav_CCnue, h_selMu_cosThetaGrav_NCnue, h_selMu_cosThetaGrav_ext = configureHists(h_selMu_cosThetaGrav_CCnumu, h_selMu_cosThetaGrav_NCnumu, h_selMu_cosThetaGrav_CCnue, h_selMu_cosThetaGrav_NCnue, h_selMu_cosThetaGrav_ext)

h_selMu_muScore_CCnumu = rt.TH1F("h_selMu_muScore_CCnumu","Reco Primary Muon's Muon Score",41,-20,0.5)
h_selMu_muScore_NCnumu = rt.TH1F("h_selMu_muScore_NCnumu","Reco Primary Muon's Muon Score",41,-20,0.5)
h_selMu_muScore_CCnue = rt.TH1F("h_selMu_muScore_CCnue","Reco Primary Muon's Muon Score",41,-20,0.5)
h_selMu_muScore_NCnue = rt.TH1F("h_selMu_muScore_NCnue","Reco Primary Muon's Muon Score",41,-20,0.5)
h_selMu_muScore_ext = rt.TH1F("h_selMu_muScore_ext","Reco Primary Muon's Muon Score",41,-20,0.5)
h_selMu_muScore_CCnumu, h_selMu_muScore_NCnumu, h_selMu_muScore_CCnue, h_selMu_muScore_NCnue, h_selMu_muScore_ext = configureHists(h_selMu_muScore_CCnumu, h_selMu_muScore_NCnumu, h_selMu_muScore_CCnue, h_selMu_muScore_NCnue, h_selMu_muScore_ext)

h_selMu_piScore_CCnumu = rt.TH1F("h_selMu_piScore_CCnumu","Reco Primary Muon's Pion Score",41,-20,0.5)
h_selMu_piScore_NCnumu = rt.TH1F("h_selMu_piScore_NCnumu","Reco Primary Muon's Pion Score",41,-20,0.5)
h_selMu_piScore_CCnue = rt.TH1F("h_selMu_piScore_CCnue","Reco Primary Muon's Pion Score",41,-20,0.5)
h_selMu_piScore_NCnue = rt.TH1F("h_selMu_piScore_NCnue","Reco Primary Muon's Pion Score",41,-20,0.5)
h_selMu_piScore_ext = rt.TH1F("h_selMu_piScore_ext","Reco Primary Muon's Pion Score",41,-20,0.5)
h_selMu_piScore_CCnumu, h_selMu_piScore_NCnumu, h_selMu_piScore_CCnue, h_selMu_piScore_NCnue, h_selMu_piScore_ext = configureHists(h_selMu_piScore_CCnumu, h_selMu_piScore_NCnumu, h_selMu_piScore_CCnue, h_selMu_piScore_NCnue, h_selMu_piScore_ext)

h_selMu_scoreConf_CCnumu = rt.TH1F("h_selMu_scoreConf_CCnumu","Reco Primary Muon's Muon Score - Pion Score",41,-0.5,20)
h_selMu_scoreConf_NCnumu = rt.TH1F("h_selMu_scoreConf_NCnumu","Reco Primary Muon's Muon Score - Pion Score",41,-0.5,20)
h_selMu_scoreConf_CCnue = rt.TH1F("h_selMu_scoreConf_CCnue","Reco Primary Muon's Muon Score - Pion Score",41,-0.5,20)
h_selMu_scoreConf_NCnue = rt.TH1F("h_selMu_scoreConf_NCnue","Reco Primary Muon's Muon Score - Pion Score",41,-0.5,20)
h_selMu_scoreConf_ext = rt.TH1F("h_selMu_scoreConf_ext","Reco Primary Muon's Muon Score - Pion Score",41,-0.5,20)
h_selMu_scoreConf_CCnumu, h_selMu_scoreConf_NCnumu, h_selMu_scoreConf_CCnue, h_selMu_scoreConf_NCnue, h_selMu_scoreConf_ext = configureHists(h_selMu_scoreConf_CCnumu, h_selMu_scoreConf_NCnumu, h_selMu_scoreConf_CCnue, h_selMu_scoreConf_NCnue, h_selMu_scoreConf_ext)

h_selMu_proc_CCnumu = rt.TH1F("h_selMu_proc_CCnumu","Reco Primary Muon's Production Process Class",3,0,3)
h_selMu_proc_NCnumu = rt.TH1F("h_selMu_proc_NCnumu","Reco Primary Muon's Production Process Class",3,0,3)
h_selMu_proc_CCnue = rt.TH1F("h_selMu_proc_CCnue","Reco Primary Muon's Production Process Class",3,0,3)
h_selMu_proc_NCnue = rt.TH1F("h_selMu_proc_NCnue","Reco Primary Muon's Production Process Class",3,0,3)
h_selMu_proc_ext = rt.TH1F("h_selMu_proc_ext","Reco Primary Muon's Production Process Class",3,0,3)
h_selMu_proc_CCnumu, h_selMu_proc_NCnumu, h_selMu_proc_CCnue, h_selMu_proc_NCnue, h_selMu_proc_ext = configureHists(h_selMu_proc_CCnumu, h_selMu_proc_NCnumu, h_selMu_proc_CCnue, h_selMu_proc_NCnue, h_selMu_proc_ext)

h_selMu_primScore_CCnumu = rt.TH1F("h_selMu_primScore_CCnumu","Reco Primary Muon's Primary Process Score",41,-20,0.5)
h_selMu_primScore_NCnumu = rt.TH1F("h_selMu_primScore_NCnumu","Reco Primary Muon's Primary Process Score",41,-20,0.5)
h_selMu_primScore_CCnue = rt.TH1F("h_selMu_primScore_CCnue","Reco Primary Muon's Primary Process Score",41,-20,0.5)
h_selMu_primScore_NCnue = rt.TH1F("h_selMu_primScore_NCnue","Reco Primary Muon's Primary Process Score",41,-20,0.5)
h_selMu_primScore_ext = rt.TH1F("h_selMu_primScore_ext","Reco Primary Muon's Primary Process Score",41,-20,0.5)
h_selMu_primScore_CCnumu, h_selMu_primScore_NCnumu, h_selMu_primScore_CCnue, h_selMu_primScore_NCnue, h_selMu_primScore_ext = configureHists(h_selMu_primScore_CCnumu, h_selMu_primScore_NCnumu, h_selMu_primScore_CCnue, h_selMu_primScore_NCnue, h_selMu_primScore_ext)

h_selMu_ntrlScore_CCnumu = rt.TH1F("h_selMu_ntrlScore_CCnumu","Reco Primary Muon's Neutral Parent Process Score",41,-20,0.5)
h_selMu_ntrlScore_NCnumu = rt.TH1F("h_selMu_ntrlScore_NCnumu","Reco Primary Muon's Neutral Parent Process Score",41,-20,0.5)
h_selMu_ntrlScore_CCnue = rt.TH1F("h_selMu_ntrlScore_CCnue","Reco Primary Muon's Neutral Parent Process Score",41,-20,0.5)
h_selMu_ntrlScore_NCnue = rt.TH1F("h_selMu_ntrlScore_NCnue","Reco Primary Muon's Neutral Parent Process Score",41,-20,0.5)
h_selMu_ntrlScore_ext = rt.TH1F("h_selMu_ntrlScore_ext","Reco Primary Muon's Neutral Parent Process Score",41,-20,0.5)
h_selMu_ntrlScore_CCnumu, h_selMu_ntrlScore_NCnumu, h_selMu_ntrlScore_CCnue, h_selMu_ntrlScore_NCnue, h_selMu_ntrlScore_ext = configureHists(h_selMu_ntrlScore_CCnumu, h_selMu_ntrlScore_NCnumu, h_selMu_ntrlScore_CCnue, h_selMu_ntrlScore_NCnue, h_selMu_ntrlScore_ext)

h_selMu_chgdScore_CCnumu = rt.TH1F("h_selMu_chgdScore_CCnumu","Reco Primary Muon's Charged Parent Process Score",41,-20,0.5)
h_selMu_chgdScore_NCnumu = rt.TH1F("h_selMu_chgdScore_NCnumu","Reco Primary Muon's Charged Parent Process Score",41,-20,0.5)
h_selMu_chgdScore_CCnue = rt.TH1F("h_selMu_chgdScore_CCnue","Reco Primary Muon's Charged Parent Process Score",41,-20,0.5)
h_selMu_chgdScore_NCnue = rt.TH1F("h_selMu_chgdScore_NCnue","Reco Primary Muon's Charged Parent Process Score",41,-20,0.5)
h_selMu_chgdScore_ext = rt.TH1F("h_selMu_chgdScore_ext","Reco Primary Muon's Charged Parent Process Score",41,-20,0.5)
h_selMu_chgdScore_CCnumu, h_selMu_chgdScore_NCnumu, h_selMu_chgdScore_CCnue, h_selMu_chgdScore_NCnue, h_selMu_chgdScore_ext = configureHists(h_selMu_chgdScore_CCnumu, h_selMu_chgdScore_NCnumu, h_selMu_chgdScore_CCnue, h_selMu_chgdScore_NCnue, h_selMu_chgdScore_ext)

h_selMu_procConf_CCnumu = rt.TH1F("h_selMu_procConf_CCnumu","Reco Primary Muon's Primary Process Score Confidence",41,-0.5,20)
h_selMu_procConf_NCnumu = rt.TH1F("h_selMu_procConf_NCnumu","Reco Primary Muon's Primary Process Score Confidence",41,-0.5,20)
h_selMu_procConf_CCnue = rt.TH1F("h_selMu_procConf_CCnue","Reco Primary Muon's Primary Process Score Confidence",41,-0.5,20)
h_selMu_procConf_NCnue = rt.TH1F("h_selMu_procConf_NCnue","Reco Primary Muon's Primary Process Score Confidence",41,-0.5,20)
h_selMu_procConf_ext = rt.TH1F("h_selMu_procConf_ext","Reco Primary Muon's Primary Process Score Confidence",41,-0.5,20)
h_selMu_procConf_CCnumu, h_selMu_procConf_NCnumu, h_selMu_procConf_CCnue, h_selMu_procConf_NCnue, h_selMu_procConf_ext = configureHists(h_selMu_procConf_CCnumu, h_selMu_procConf_NCnumu, h_selMu_procConf_CCnue, h_selMu_procConf_NCnue, h_selMu_procConf_ext)



h_nuE_CCnumu_wCuts = rt.TH1F("h_nuE_CCnumu_wCuts","Neutrino Energy for True CCnumu Events",30,0,3)
h_nuE_CCnumu_nCuts = rt.TH1F("h_nuE_CCnue_nCuts","Neutrino Energy for True CCnue Events",30,0,3)
h_nuE_CCnue_wCuts = rt.TH1F("h_nuE_CCnue_wCuts","Neutrino Energy for True CCnue Events",30,0,3)
h_nuE_NCnumu_wCuts = rt.TH1F("h_nuE_NCnumu_wCuts","Neutrino Energy for True NCnumu Events",30,0,3)
h_nuE_NCnue_wCuts = rt.TH1F("h_nuE_NCnue_wCuts","Neutrino Energy for True NCnue Events",30,0,3)
h_nuE_all_wCuts = rt.TH1F("h_nuE_all_wCuts","Neutrino Energy for all Events",30,0,3)
h_nuE_CCnumu_nCuts.GetYaxis().SetTitle("events per "+targetPOTstring+" POT")
h_nuE_CCnumu_nCuts.SetLineColor(rt.kBlue)
h_nuE_CCnumu_nCuts.SetLineWidth(2)
h_nuE_CCnumu_nCuts.GetXaxis().SetTitle("true neutrino energy (GeV)")
h_nuE_CCnumu_wCuts.GetYaxis().SetTitle("events per "+targetPOTstring+" POT")
h_nuE_CCnumu_wCuts.SetLineColor(rt.kRed)
h_nuE_CCnumu_wCuts.SetLineWidth(2)
h_nuE_CCnumu_wCuts.GetXaxis().SetTitle("true neutrino energy (GeV)")
h_nuE_all_wCuts.GetYaxis().SetTitle("events per "+targetPOTstring+" POT")
h_nuE_all_wCuts.SetLineColor(rt.kRed)
h_nuE_all_wCuts.SetLineWidth(2)
h_nuE_all_wCuts.GetXaxis().SetTitle("true neutrino energy (GeV)")
h_nuE_CCnumu_eff = rt.TH1F("h_nuE_CCnumu_eff","Inclusive CC numu Selection MC Predictions",30,0,3)
h_nuE_CCnumu_eff.GetXaxis().SetTitle("true neutrino energy (GeV)")
h_nuE_CCnumu_eff.SetLineColor(rt.kBlack)
h_nuE_CCnumu_eff.SetLineWidth(2)
h_nuE_CCnumu_pur = rt.TH1F("h_nuE_CCnumu_pur","Inclusive CC numu Selection MC Predictions",30,0,3)
h_nuE_CCnumu_pur.GetXaxis().SetTitle("true neutrino energy (GeV)")
h_nuE_CCnumu_pur.SetLineColor(8)
h_nuE_CCnumu_pur.SetLineWidth(2)


def configureEffHist(hist, color):
  hist.GetXaxis().SetTitle("true neutrino energy (GeV)")
  hist.SetLineColor(color)
  hist.SetLineWidth(2)
  return hist

h_nuE_CCnumu_wCuts1 = rt.TH1F("h_nuE_CCnumu_wCuts1","Neutrino Energy for True CCnumu Events",30,0,3)
h_nuE_CCnumu_wCuts2 = rt.TH1F("h_nuE_CCnumu_wCuts2","Neutrino Energy for True CCnumu Events",30,0,3)
h_nuE_CCnumu_wCutsAll = rt.TH1F("h_nuE_CCnumu_wCutsAll","Neutrino Energy for True CCnumu Events",30,0,3)
h_nuE_CCnumu_eff1 = rt.TH1F("h_nuE_CCnumu_eff1","Inclusive CCnumu Efficiency",30,0,3)
h_nuE_CCnumu_eff2 = rt.TH1F("h_nuE_CCnumu_eff2","Inclusive CCnumu Efficiency",30,0,3)
h_nuE_CCnumu_effAll = rt.TH1F("h_nuE_CCnumu_effAll","Inclusive CCnumu Efficiency",30,0,3)
h_nuE_CCnumu_eff1 = configureEffHist(h_nuE_CCnumu_eff1, rt.kBlue)
h_nuE_CCnumu_eff2 = configureEffHist(h_nuE_CCnumu_eff2, 8)
h_nuE_CCnumu_effAll = configureEffHist(h_nuE_CCnumu_effAll, rt.kRed)

h_nuEr_CCnumu_wCuts = rt.TH1F("h_nuEr_CCnumu_wCuts","Neutrino Energy for True CCnumu Events",30,0,3)
h_nuEr_CCnumu_nCuts = rt.TH1F("h_nuEr_CCnue_nCuts","Neutrino Energy for True CCnue Events",30,0,3)
h_nuEr_CCnue_wCuts = rt.TH1F("h_nuEr_CCnue_wCuts","Neutrino Energy for True CCnue Events",30,0,3)
h_nuEr_NCnumu_wCuts = rt.TH1F("h_nuEr_NCnumu_wCuts","Neutrino Energy for True NCnumu Events",30,0,3)
h_nuEr_NCnue_wCuts = rt.TH1F("h_nuEr_NCnue_wCuts","Neutrino Energy for True NCnue Events",30,0,3)
h_nuEr_ext_wCuts = rt.TH1F("h_nuEr_ext_wCuts","Neutrino Energy for ExtBNB Events",30,0,3)
h_nuEr_all_wCuts = rt.TH1F("h_nuEr_all_wCuts","Neutrino Energy for all Events",30,0,3)
h_nuEr_CCnumu_nCuts.GetYaxis().SetTitle("events per "+targetPOTstring+" POT")
h_nuEr_CCnumu_nCuts.SetLineColor(rt.kBlue)
h_nuEr_CCnumu_nCuts.SetLineWidth(2)
h_nuEr_CCnumu_nCuts.GetXaxis().SetTitle("reco neutrino energy (GeV)")
h_nuEr_CCnumu_wCuts.GetYaxis().SetTitle("events per "+targetPOTstring+" POT")
h_nuEr_CCnumu_wCuts.SetLineColor(rt.kRed)
h_nuEr_CCnumu_wCuts.SetLineWidth(2)
h_nuEr_CCnumu_wCuts.GetXaxis().SetTitle("reco neutrino energy (GeV)")
h_nuEr_all_wCuts.GetYaxis().SetTitle("events per "+targetPOTstring+" POT")
h_nuEr_all_wCuts.SetLineColor(rt.kRed)
h_nuEr_all_wCuts.SetLineWidth(2)
h_nuEr_all_wCuts.GetXaxis().SetTitle("reco neutrino energy (GeV)")
h_nuEr_CCnumu_eff = rt.TH1F("h_nuEr_CCnumu_eff","Inclusive CC numu Selection MC Predictions",30,0,3)
h_nuEr_CCnumu_eff.GetXaxis().SetTitle("reco neutrino energy (GeV)")
h_nuEr_CCnumu_eff.SetLineColor(rt.kBlack)
h_nuEr_CCnumu_eff.SetLineWidth(2)
h_nuEr_CCnumu_pur = rt.TH1F("h_nuEr_CCnumu_pur","Inclusive CC numu Selection MC Predictions",30,0,3)
h_nuEr_CCnumu_pur.GetXaxis().SetTitle("reco neutrino energy (GeV)")
h_nuEr_CCnumu_pur.SetLineColor(8)
h_nuEr_CCnumu_pur.SetLineWidth(2)

visEn = 60
visEl = 0.
visEh = 6.
if args.recoEOverflow:
  visEn = 21
  visEl = 0.
  visEh = 2.1

h_visE_CCnue_wCuts = rt.TH1F("h_visE_CCnue_wCuts","Reco Nu Energy for True CCnue Events",visEn,visEl,visEh)
h_visE_CCnumu_wCuts = rt.TH1F("h_visE_CCnumu_wCuts","Reco Nu Energy for True CCnumu Events",visEn,visEl,visEh)
h_visE_NCnumu_wCuts = rt.TH1F("h_visE_NCnumu_wCuts","Reco Nu Energy for True NCnumu Events",visEn,visEl,visEh)
h_visE_NCnue_wCuts = rt.TH1F("h_visE_NCnue_wCuts","Reco Nu Energy for True NCnue Events",visEn,visEl,visEh)
h_visE_ext_wCuts = rt.TH1F("h_visE_ext_wCuts","Reco Nu Energy for ExtBNB Events",visEn,visEl,visEh)
h_visE_data_wCuts = rt.TH1F("h_visE_data_wCuts","Reco Nu Energy for BNB Data Events",visEn,visEl,visEh)
h_visE_predErr_wCuts = rt.TH1F("h_visE_predErr_wCuts", "Inclusive CCnumu Selected Events",visEn,visEl,visEh)

h_visE_CCnue_wCutsSet1 = rt.TH1F("h_visE_CCnue_wCutsSet1","Reco Nu Energy for True CCnue Events",visEn,visEl,visEh)
h_visE_CCnumu_wCutsSet1 = rt.TH1F("h_visE_CCnumu_wCutsSet1","Reco Nu Energy for True CCnumu Events",visEn,visEl,visEh)
h_visE_NCnumu_wCutsSet1 = rt.TH1F("h_visE_NCnumu_wCutsSet1","Reco Nu Energy for True NCnumu Events",visEn,visEl,visEh)
h_visE_NCnue_wCutsSet1 = rt.TH1F("h_visE_NCnue_wCutsSet1","Reco Nu Energy for True NCnue Events",visEn,visEl,visEh)
h_visE_ext_wCutsSet1 = rt.TH1F("h_visE_ext_wCutsSet1","Reco Nu Energy for ExtBNB Events",visEn,visEl,visEh)
h_visE_data_wCutsSet1 = rt.TH1F("h_visE_data_wCutsSet1","Reco Nu Energy for BNB Data Events",visEn,visEl,visEh)

h_visE_CCnue_wCutsSet2 = rt.TH1F("h_visE_CCnue_wCutsSet2","Reco Nu Energy for True CCnue Events",visEn,visEl,visEh)
h_visE_CCnumu_wCutsSet2 = rt.TH1F("h_visE_CCnumu_wCutsSet2","Reco Nu Energy for True CCnumu Events",visEn,visEl,visEh)
h_visE_NCnumu_wCutsSet2 = rt.TH1F("h_visE_NCnumu_wCutsSet2","Reco Nu Energy for True NCnumu Events",visEn,visEl,visEh)
h_visE_NCnue_wCutsSet2 = rt.TH1F("h_visE_NCnue_wCutsSet2","Reco Nu Energy for True NCnue Events",visEn,visEl,visEh)
h_visE_ext_wCutsSet2 = rt.TH1F("h_visE_ext_wCutsSet2","Reco Nu Energy for ExtBNB Events",visEn,visEl,visEh)
h_visE_data_wCutsSet2 = rt.TH1F("h_visE_data_wCutsSet2","Reco Nu Energy for BNB Data Events",visEn,visEl,visEh)

h_visE_CCnue_wCutsSet3 = rt.TH1F("h_visE_CCnue_wCutsSet3","Reco Nu Energy for True CCnue Events",visEn,visEl,visEh)
h_visE_CCnumu_wCutsSet3 = rt.TH1F("h_visE_CCnumu_wCutsSet3","Reco Nu Energy for True CCnumu Events",visEn,visEl,visEh)
h_visE_NCnumu_wCutsSet3 = rt.TH1F("h_visE_NCnumu_wCutsSet3","Reco Nu Energy for True NCnumu Events",visEn,visEl,visEh)
h_visE_NCnue_wCutsSet3 = rt.TH1F("h_visE_NCnue_wCutsSet3","Reco Nu Energy for True NCnue Events",visEn,visEl,visEh)
h_visE_ext_wCutsSet3 = rt.TH1F("h_visE_ext_wCutsSet3","Reco Nu Energy for ExtBNB Events",visEn,visEl,visEh)
h_visE_data_wCutsSet3 = rt.TH1F("h_visE_data_wCutsSet3","Reco Nu Energy for BNB Data Events",visEn,visEl,visEh)

h_visE_CCnue_wCutsSet4 = rt.TH1F("h_visE_CCnue_wCutsSet4","Reco Nu Energy for True CCnue Events",visEn,visEl,visEh)
h_visE_CCnumu_wCutsSet4 = rt.TH1F("h_visE_CCnumu_wCutsSet4","Reco Nu Energy for True CCnumu Events",visEn,visEl,visEh)
h_visE_NCnumu_wCutsSet4 = rt.TH1F("h_visE_NCnumu_wCutsSet4","Reco Nu Energy for True NCnumu Events",visEn,visEl,visEh)
h_visE_NCnue_wCutsSet4 = rt.TH1F("h_visE_NCnue_wCutsSet4","Reco Nu Energy for True NCnue Events",visEn,visEl,visEh)
h_visE_ext_wCutsSet4 = rt.TH1F("h_visE_ext_wCutsSet4","Reco Nu Energy for ExtBNB Events",visEn,visEl,visEh)
h_visE_data_wCutsSet4 = rt.TH1F("h_visE_data_wCutsSet4","Reco Nu Energy for BNB Data Events",visEn,visEl,visEh)

h_visE_CCnue_wCutsSet5 = rt.TH1F("h_visE_CCnue_wCutsSet5","Reco Nu Energy for True CCnue Events",visEn,visEl,visEh)
h_visE_CCnumu_wCutsSet5 = rt.TH1F("h_visE_CCnumu_wCutsSet5","Reco Nu Energy for True CCnumu Events",visEn,visEl,visEh)
h_visE_NCnumu_wCutsSet5 = rt.TH1F("h_visE_NCnumu_wCutsSet5","Reco Nu Energy for True NCnumu Events",visEn,visEl,visEh)
h_visE_NCnue_wCutsSet5 = rt.TH1F("h_visE_NCnue_wCutsSet5","Reco Nu Energy for True NCnue Events",visEn,visEl,visEh)
h_visE_ext_wCutsSet5 = rt.TH1F("h_visE_ext_wCutsSet5","Reco Nu Energy for ExtBNB Events",visEn,visEl,visEh)
h_visE_data_wCutsSet5 = rt.TH1F("h_visE_data_wCutsSet5","Reco Nu Energy for BNB Data Events",visEn,visEl,visEh)

h_visE_CCnue_wCutsSet6 = rt.TH1F("h_visE_CCnue_wCutsSet6","Reco Nu Energy for True CCnue Events",visEn,visEl,visEh)
h_visE_CCnumu_wCutsSet6 = rt.TH1F("h_visE_CCnumu_wCutsSet6","Reco Nu Energy for True CCnumu Events",visEn,visEl,visEh)
h_visE_NCnumu_wCutsSet6 = rt.TH1F("h_visE_NCnumu_wCutsSet6","Reco Nu Energy for True NCnumu Events",visEn,visEl,visEh)
h_visE_NCnue_wCutsSet6 = rt.TH1F("h_visE_NCnue_wCutsSet6","Reco Nu Energy for True NCnue Events",visEn,visEl,visEh)
h_visE_ext_wCutsSet6 = rt.TH1F("h_visE_ext_wCutsSet6","Reco Nu Energy for ExtBNB Events",visEn,visEl,visEh)
h_visE_data_wCutsSet6 = rt.TH1F("h_visE_data_wCutsSet6","Reco Nu Energy for BNB Data Events",visEn,visEl,visEh)

def configureStackedHists(h_CCnumu, h_NCnumu, h_CCnue, h_NCnue, h_ext, h_data, title, xtitle):
  h_CCnumu.SetFillColor(rt.kRed)
  h_NCnumu.SetFillColor(8)
  h_CCnue.SetFillColor(rt.kBlue)
  h_NCnue.SetFillColor(40)
  h_ext.SetFillColor(12)
  h_data.SetLineColor(rt.kBlack)
  h_data.SetLineWidth(2)
  h_data.SetTitle(title)
  h_data.GetYaxis().SetTitle("events per "+targetPOTstring+" POT")
  h_data.GetXaxis().SetTitle(xtitle)
  if args.recoEOverflow and "GeV" in xtitle:
    h_data.GetXaxis().CenterTitle(True)
  return h_CCnumu, h_NCnumu, h_CCnue, h_NCnue, h_ext, h_data

h_visE_CCnumu_wCuts, h_visE_NCnumu_wCuts, h_visE_CCnue_wCuts, h_visE_NCnue_wCuts, h_visE_ext_wCuts, h_visE_data_wCuts = configureStackedHists(h_visE_CCnumu_wCuts, h_visE_NCnumu_wCuts, h_visE_CCnue_wCuts, h_visE_NCnue_wCuts, h_visE_ext_wCuts, h_visE_data_wCuts, "Inclusive CCnumu Selected Events", "reconstructed neutrino energy (GeV)")
h_visE_all_wCuts = rt.THStack("h_visE_all_wCuts", "Inclusive CCnumu Selected Events")

h_visE_CCnumu_wCutsSet1, h_visE_NCnumu_wCutsSet1, h_visE_CCnue_wCutsSet1, h_visE_NCnue_wCutsSet1, h_visE_ext_wCutsSet1, h_visE_data_wCutsSet1 = configureStackedHists(h_visE_CCnumu_wCutsSet1, h_visE_NCnumu_wCutsSet1, h_visE_CCnue_wCutsSet1, h_visE_NCnue_wCutsSet1, h_visE_ext_wCutsSet1, h_visE_data_wCutsSet1, "Inclusive CCnumu Selected Events (Cut Subset 1)", "reconstructed neutrino energy (GeV)")
h_visE_all_wCutsSet1 = rt.THStack("h_visE_all_wCutsSet1", "Inclusive CCnumu Selected Events (Cut Subset 1)")
h_visE_CCnumu_wCutsSet2, h_visE_NCnumu_wCutsSet2, h_visE_CCnue_wCutsSet2, h_visE_NCnue_wCutsSet2, h_visE_ext_wCutsSet2, h_visE_data_wCutsSet2 = configureStackedHists(h_visE_CCnumu_wCutsSet2, h_visE_NCnumu_wCutsSet2, h_visE_CCnue_wCutsSet2, h_visE_NCnue_wCutsSet2, h_visE_ext_wCutsSet2, h_visE_data_wCutsSet2, "Inclusive CCnumu Selected Events (Cut Subset 2)", "reconstructed neutrino energy (GeV)")
h_visE_all_wCutsSet2 = rt.THStack("h_visE_all_wCutsSet2", "Inclusive CCnumu Selected Events (Cut Subset 2)")
h_visE_CCnumu_wCutsSet3, h_visE_NCnumu_wCutsSet3, h_visE_CCnue_wCutsSet3, h_visE_NCnue_wCutsSet3, h_visE_ext_wCutsSet3, h_visE_data_wCutsSet3 = configureStackedHists(h_visE_CCnumu_wCutsSet3, h_visE_NCnumu_wCutsSet3, h_visE_CCnue_wCutsSet3, h_visE_NCnue_wCutsSet3, h_visE_ext_wCutsSet3, h_visE_data_wCutsSet3, "Inclusive CCnumu Selected Events (Cut Subset 3)", "reconstructed neutrino energy (GeV)")
h_visE_all_wCutsSet3 = rt.THStack("h_visE_all_wCutsSet3", "Inclusive CCnumu Selected Events (Cut Subset 3)")
h_visE_CCnumu_wCutsSet4, h_visE_NCnumu_wCutsSet4, h_visE_CCnue_wCutsSet4, h_visE_NCnue_wCutsSet4, h_visE_ext_wCutsSet4, h_visE_data_wCutsSet4 = configureStackedHists(h_visE_CCnumu_wCutsSet4, h_visE_NCnumu_wCutsSet4, h_visE_CCnue_wCutsSet4, h_visE_NCnue_wCutsSet4, h_visE_ext_wCutsSet4, h_visE_data_wCutsSet4, "Inclusive CCnumu Selected Events (Cut Subset 4)", "reconstructed neutrino energy (GeV)")
h_visE_all_wCutsSet4 = rt.THStack("h_visE_all_wCutsSet4", "Inclusive CCnumu Selected Events (Cut Subset 4)")
h_visE_CCnumu_wCutsSet5, h_visE_NCnumu_wCutsSet5, h_visE_CCnue_wCutsSet5, h_visE_NCnue_wCutsSet5, h_visE_ext_wCutsSet5, h_visE_data_wCutsSet5 = configureStackedHists(h_visE_CCnumu_wCutsSet5, h_visE_NCnumu_wCutsSet5, h_visE_CCnue_wCutsSet5, h_visE_NCnue_wCutsSet5, h_visE_ext_wCutsSet5, h_visE_data_wCutsSet5, "Inclusive CCnumu Selected Events (Cut Subset 5)", "reconstructed neutrino energy (GeV)")
h_visE_all_wCutsSet5 = rt.THStack("h_visE_all_wCutsSet5", "Inclusive CCnumu Selected Events (Cut Subset 5)")
h_visE_CCnumu_wCutsSet6, h_visE_NCnumu_wCutsSet6, h_visE_CCnue_wCutsSet6, h_visE_NCnue_wCutsSet6, h_visE_ext_wCutsSet6, h_visE_data_wCutsSet6 = configureStackedHists(h_visE_CCnumu_wCutsSet6, h_visE_NCnumu_wCutsSet6, h_visE_CCnue_wCutsSet6, h_visE_NCnue_wCutsSet6, h_visE_ext_wCutsSet6, h_visE_data_wCutsSet6, "Inclusive CCnumu Selected Events (Cut Subset 6)", "reconstructed neutrino energy (GeV)")
h_visE_all_wCutsSet6 = rt.THStack("h_visE_all_wCutsSet6", "Inclusive CCnumu Selected Events (Cut Subset 6)")


cosThetan = 18
cosThetal = -1.125
cosThetah = 1.125
h_cosTheta_CCnue_wCuts = rt.TH1F("h_cosTheta_CCnue_wCuts","Reco Muon cos(theta) for True CCnue Events",cosThetan,cosThetal,cosThetah)
h_cosTheta_CCnumu_wCuts = rt.TH1F("h_cosTheta_CCnumu_wCuts","Reco Muon cos(theta) for True CCnumu Events",cosThetan,cosThetal,cosThetah)
h_cosTheta_NCnumu_wCuts = rt.TH1F("h_cosTheta_NCnumu_wCuts","Reco Muon cos(theta) for True NCnumu Events",cosThetan,cosThetal,cosThetah)
h_cosTheta_NCnue_wCuts = rt.TH1F("h_cosTheta_NCnue_wCuts","Reco Muon cos(theta) for True NCnue Events",cosThetan,cosThetal,cosThetah)
h_cosTheta_ext_wCuts = rt.TH1F("h_cosTheta_ext_wCuts","Reco Muon cos(theta) for ExtBNB Events",cosThetan,cosThetal,cosThetah)
h_cosTheta_data_wCuts = rt.TH1F("h_cosTheta_data_wCuts","Reco Muon cos(theta) for BNB Data Events",cosThetan,cosThetal,cosThetah)
h_cosTheta_all_wCuts = rt.THStack("h_cosTheta_all_wCuts", "Inclusive CCnumu Selected Events")
h_cosTheta_predErr_wCuts = rt.TH1F("h_cosTheta_predErr_wCuts", "Inclusive CCnumu Selected Events",cosThetan,cosThetal,cosThetah)
h_cosTheta_CCnumu_wCuts, h_cosTheta_NCnumu_wCuts, h_cosTheta_CCnue_wCuts, h_cosTheta_NCnue_wCuts, h_cosTheta_ext_wCuts, h_cosTheta_data_wCuts = configureStackedHists(h_cosTheta_CCnumu_wCuts, h_cosTheta_NCnumu_wCuts, h_cosTheta_CCnue_wCuts, h_cosTheta_NCnue_wCuts, h_cosTheta_ext_wCuts, h_cosTheta_data_wCuts, "Inclusive CCnumu Selected Events", "reconstructed muon cos(theta)")

lepPn = 60
lepPl = 0.
lepPh = 6.
if args.recoEOverflow:
  lepPn = 21
  lepPl = 0.
  lepPh = 2.1
h_lepP_CCnue_wCuts = rt.TH1F("h_lepP_CCnue_wCuts","Reco Muon Momentum for True CCnue Events",lepPn,lepPl,lepPh)
h_lepP_CCnumu_wCuts = rt.TH1F("h_lepP_CCnumu_wCuts","Reco Muon Momentum for True CCnumu Events",lepPn,lepPl,lepPh)
h_lepP_NCnumu_wCuts = rt.TH1F("h_lepP_NCnumu_wCuts","Reco Muon Momentum for True NCnumu Events",lepPn,lepPl,lepPh)
h_lepP_NCnue_wCuts = rt.TH1F("h_lepP_NCnue_wCuts","Reco Muon Momentum for True NCnue Events",lepPn,lepPl,lepPh)
h_lepP_ext_wCuts = rt.TH1F("h_lepP_ext_wCuts","Reco Muon Momentum for ExtBNB Events",lepPn,lepPl,lepPh)
h_lepP_data_wCuts = rt.TH1F("h_lepP_data_wCuts","Reco Muon Momentum for BNB Data Events",lepPn,lepPl,lepPh)
h_lepP_all_wCuts = rt.THStack("h_lepP_all_wCuts", "Inclusive CCnumu Selected Events")
h_lepP_predErr_wCuts = rt.TH1F("h_lepP_predErr_wCuts", "Inclusive CCnumu Selected Events",lepPn,lepPl,lepPh)
h_lepP_CCnumu_wCuts, h_lepP_NCnumu_wCuts, h_lepP_CCnue_wCuts, h_lepP_NCnue_wCuts, h_lepP_ext_wCuts, h_lepP_data_wCuts = configureStackedHists(h_lepP_CCnumu_wCuts, h_lepP_NCnumu_wCuts, h_lepP_CCnue_wCuts, h_lepP_NCnue_wCuts, h_lepP_ext_wCuts, h_lepP_data_wCuts, "Inclusive CCnumu Selected Events", "reconstructed muon momentum (GeV/c)")

partScore_n = 47
partScore_l = -23.
partScore_h = 0.5

h_elScr_CCnue_wCuts = rt.TH1F("h_elScr_CCnue_wCuts","Reco Muon's Electron Score for True CCnue Events",partScore_n,partScore_l,partScore_h)
h_elScr_CCnumu_wCuts = rt.TH1F("h_elScr_CCnumu_wCuts","Reco Muon's Electron Score for True CCnumu Events",partScore_n,partScore_l,partScore_h)
h_elScr_NCnumu_wCuts = rt.TH1F("h_elScr_NCnumu_wCuts","Reco Muon's Electron Score for True NCnumu Events",partScore_n,partScore_l,partScore_h)
h_elScr_NCnue_wCuts = rt.TH1F("h_elScr_NCnue_wCuts","Reco Muon's Electron Score for True NCnue Events",partScore_n,partScore_l,partScore_h)
h_elScr_ext_wCuts = rt.TH1F("h_elScr_ext_wCuts","Reco Muon's Electron Score for ExtBNB Events",partScore_n,partScore_l,partScore_h)
h_elScr_data_wCuts = rt.TH1F("h_elScr_data_wCuts","Reco Muon's Electron Score for BNB Data Events",partScore_n,partScore_l,partScore_h)
h_elScr_all_wCuts = rt.THStack("h_elScr_all_wCuts", "Inclusive CCnumu Selected Events")
h_elScr_CCnumu_wCuts, h_elScr_NCnumu_wCuts, h_elScr_CCnue_wCuts, h_elScr_NCnue_wCuts, h_elScr_ext_wCuts, h_elScr_data_wCuts = configureStackedHists(h_elScr_CCnumu_wCuts, h_elScr_NCnumu_wCuts, h_elScr_CCnue_wCuts, h_elScr_NCnue_wCuts, h_elScr_ext_wCuts, h_elScr_data_wCuts, "Inclusive CCnumu Selected Events", "reconstructed muon's electron score")

h_phScr_CCnue_wCuts = rt.TH1F("h_phScr_CCnue_wCuts","Reco Muon's Photon Score for True CCnue Events",partScore_n,partScore_l,partScore_h)
h_phScr_CCnumu_wCuts = rt.TH1F("h_phScr_CCnumu_wCuts","Reco Muon's Photon Score for True CCnumu Events",partScore_n,partScore_l,partScore_h)
h_phScr_NCnumu_wCuts = rt.TH1F("h_phScr_NCnumu_wCuts","Reco Muon's Photon Score for True NCnumu Events",partScore_n,partScore_l,partScore_h)
h_phScr_NCnue_wCuts = rt.TH1F("h_phScr_NCnue_wCuts","Reco Muon's Photon Score for True NCnue Events",partScore_n,partScore_l,partScore_h)
h_phScr_ext_wCuts = rt.TH1F("h_phScr_ext_wCuts","Reco Muon's Photon Score for ExtBNB Events",partScore_n,partScore_l,partScore_h)
h_phScr_data_wCuts = rt.TH1F("h_phScr_data_wCuts","Reco Muon's Photon Score for BNB Data Events",partScore_n,partScore_l,partScore_h)
h_phScr_all_wCuts = rt.THStack("h_phScr_all_wCuts", "Inclusive CCnumu Selected Events")
h_phScr_CCnumu_wCuts, h_phScr_NCnumu_wCuts, h_phScr_CCnue_wCuts, h_phScr_NCnue_wCuts, h_phScr_ext_wCuts, h_phScr_data_wCuts = configureStackedHists(h_phScr_CCnumu_wCuts, h_phScr_NCnumu_wCuts, h_phScr_CCnue_wCuts, h_phScr_NCnue_wCuts, h_phScr_ext_wCuts, h_phScr_data_wCuts, "Inclusive CCnumu Selected Events", "reconstructed muon's photon score")

h_piScr_CCnue_wCuts = rt.TH1F("h_piScr_CCnue_wCuts","Reco Muon's Pion Score for True CCnue Events",partScore_n,partScore_l,partScore_h)
h_piScr_CCnumu_wCuts = rt.TH1F("h_piScr_CCnumu_wCuts","Reco Muon's Pion Score for True CCnumu Events",partScore_n,partScore_l,partScore_h)
h_piScr_NCnumu_wCuts = rt.TH1F("h_piScr_NCnumu_wCuts","Reco Muon's Pion Score for True NCnumu Events",partScore_n,partScore_l,partScore_h)
h_piScr_NCnue_wCuts = rt.TH1F("h_piScr_NCnue_wCuts","Reco Muon's Pion Score for True NCnue Events",partScore_n,partScore_l,partScore_h)
h_piScr_ext_wCuts = rt.TH1F("h_piScr_ext_wCuts","Reco Muon's Pion Score for ExtBNB Events",partScore_n,partScore_l,partScore_h)
h_piScr_data_wCuts = rt.TH1F("h_piScr_data_wCuts","Reco Muon's Pion Score for BNB Data Events",partScore_n,partScore_l,partScore_h)
h_piScr_all_wCuts = rt.THStack("h_piScr_all_wCuts", "Inclusive CCnumu Selected Events")
h_piScr_CCnumu_wCuts, h_piScr_NCnumu_wCuts, h_piScr_CCnue_wCuts, h_piScr_NCnue_wCuts, h_piScr_ext_wCuts, h_piScr_data_wCuts = configureStackedHists(h_piScr_CCnumu_wCuts, h_piScr_NCnumu_wCuts, h_piScr_CCnue_wCuts, h_piScr_NCnue_wCuts, h_piScr_ext_wCuts, h_piScr_data_wCuts, "Inclusive CCnumu Selected Events", "reconstructed muon's pion score")

h_muScr_CCnue_wCuts = rt.TH1F("h_muScr_CCnue_wCuts","Reco Muon's Muon Score for True CCnue Events",partScore_n,partScore_l,partScore_h)
h_muScr_CCnumu_wCuts = rt.TH1F("h_muScr_CCnumu_wCuts","Reco Muon's Muon Score for True CCnumu Events",partScore_n,partScore_l,partScore_h)
h_muScr_NCnumu_wCuts = rt.TH1F("h_muScr_NCnumu_wCuts","Reco Muon's Muon Score for True NCnumu Events",partScore_n,partScore_l,partScore_h)
h_muScr_NCnue_wCuts = rt.TH1F("h_muScr_NCnue_wCuts","Reco Muon's Muon Score for True NCnue Events",partScore_n,partScore_l,partScore_h)
h_muScr_ext_wCuts = rt.TH1F("h_muScr_ext_wCuts","Reco Muon's Muon Score for ExtBNB Events",partScore_n,partScore_l,partScore_h)
h_muScr_data_wCuts = rt.TH1F("h_muScr_data_wCuts","Reco Muon's Muon Score for BNB Data Events",partScore_n,partScore_l,partScore_h)
h_muScr_all_wCuts = rt.THStack("h_muScr_all_wCuts", "Inclusive CCnumu Selected Events")
h_muScr_CCnumu_wCuts, h_muScr_NCnumu_wCuts, h_muScr_CCnue_wCuts, h_muScr_NCnue_wCuts, h_muScr_ext_wCuts, h_muScr_data_wCuts = configureStackedHists(h_muScr_CCnumu_wCuts, h_muScr_NCnumu_wCuts, h_muScr_CCnue_wCuts, h_muScr_NCnue_wCuts, h_muScr_ext_wCuts, h_muScr_data_wCuts, "Inclusive CCnumu Selected Events", "reconstructed muon's muon score")

h_prScr_CCnue_wCuts = rt.TH1F("h_prScr_CCnue_wCuts","Reco Muon's Proton Score for True CCnue Events",partScore_n,partScore_l,partScore_h)
h_prScr_CCnumu_wCuts = rt.TH1F("h_prScr_CCnumu_wCuts","Reco Muon's Proton Score for True CCnumu Events",partScore_n,partScore_l,partScore_h)
h_prScr_NCnumu_wCuts = rt.TH1F("h_prScr_NCnumu_wCuts","Reco Muon's Proton Score for True NCnumu Events",partScore_n,partScore_l,partScore_h)
h_prScr_NCnue_wCuts = rt.TH1F("h_prScr_NCnue_wCuts","Reco Muon's Proton Score for True NCnue Events",partScore_n,partScore_l,partScore_h)
h_prScr_ext_wCuts = rt.TH1F("h_prScr_ext_wCuts","Reco Muon's Proton Score for ExtBNB Events",partScore_n,partScore_l,partScore_h)
h_prScr_data_wCuts = rt.TH1F("h_prScr_data_wCuts","Reco Muon's Proton Score for BNB Data Events",partScore_n,partScore_l,partScore_h)
h_prScr_all_wCuts = rt.THStack("h_prScr_all_wCuts", "Inclusive CCnumu Selected Events")
h_prScr_CCnumu_wCuts, h_prScr_NCnumu_wCuts, h_prScr_CCnue_wCuts, h_prScr_NCnue_wCuts, h_prScr_ext_wCuts, h_prScr_data_wCuts = configureStackedHists(h_prScr_CCnumu_wCuts, h_prScr_NCnumu_wCuts, h_prScr_CCnue_wCuts, h_prScr_NCnue_wCuts, h_prScr_ext_wCuts, h_prScr_data_wCuts, "Inclusive CCnumu Selected Events", "reconstructed muon's proton score")

procScore_n = 27
procScore_l = -13.
procScore_h = 0.5

h_pPScr_CCnue_wCuts = rt.TH1F("h_pPScr_CCnue_wCuts","Reco Muon's Primary Score for True CCnue Events",procScore_n,procScore_l,procScore_h)
h_pPScr_CCnumu_wCuts = rt.TH1F("h_pPScr_CCnumu_wCuts","Reco Muon's Primary Score for True CCnumu Events",procScore_n,procScore_l,procScore_h)
h_pPScr_NCnumu_wCuts = rt.TH1F("h_pPScr_NCnumu_wCuts","Reco Muon's Primary Score for True NCnumu Events",procScore_n,procScore_l,procScore_h)
h_pPScr_NCnue_wCuts = rt.TH1F("h_pPScr_NCnue_wCuts","Reco Muon's Primary Score for True NCnue Events",procScore_n,procScore_l,procScore_h)
h_pPScr_ext_wCuts = rt.TH1F("h_pPScr_ext_wCuts","Reco Muon's Primary Score for ExtBNB Events",procScore_n,procScore_l,procScore_h)
h_pPScr_data_wCuts = rt.TH1F("h_pPScr_data_wCuts","Reco Muon's Primary Score for BNB Data Events",procScore_n,procScore_l,procScore_h)
h_pPScr_all_wCuts = rt.THStack("h_pPScr_all_wCuts", "Inclusive CCnumu Selected Events")
h_pPScr_CCnumu_wCuts, h_pPScr_NCnumu_wCuts, h_pPScr_CCnue_wCuts, h_pPScr_NCnue_wCuts, h_pPScr_ext_wCuts, h_pPScr_data_wCuts = configureStackedHists(h_pPScr_CCnumu_wCuts, h_pPScr_NCnumu_wCuts, h_pPScr_CCnue_wCuts, h_pPScr_NCnue_wCuts, h_pPScr_ext_wCuts, h_pPScr_data_wCuts, "Inclusive CCnumu Selected Events", "reconstructed muon's primary score")

h_pNScr_CCnue_wCuts = rt.TH1F("h_pNScr_CCnue_wCuts","Reco Muon's Secondary (Neutral Parent) Score for True CCnue Events",procScore_n,procScore_l,procScore_h)
h_pNScr_CCnumu_wCuts = rt.TH1F("h_pNScr_CCnumu_wCuts","Reco Muon's Secondary (Neutral Parent) Score for True CCnumu Events",procScore_n,procScore_l,procScore_h)
h_pNScr_NCnumu_wCuts = rt.TH1F("h_pNScr_NCnumu_wCuts","Reco Muon's Secondary (Neutral Parent) Score for True NCnumu Events",procScore_n,procScore_l,procScore_h)
h_pNScr_NCnue_wCuts = rt.TH1F("h_pNScr_NCnue_wCuts","Reco Muon's Secondary (Neutral Parent) Score for True NCnue Events",procScore_n,procScore_l,procScore_h)
h_pNScr_ext_wCuts = rt.TH1F("h_pNScr_ext_wCuts","Reco Muon's Secondary (Neutral Parent) Score for ExtBNB Events",procScore_n,procScore_l,procScore_h)
h_pNScr_data_wCuts = rt.TH1F("h_pNScr_data_wCuts","Reco Muon's Secondary (Neutral Parent) Score for BNB Data Events",procScore_n,procScore_l,procScore_h)
h_pNScr_all_wCuts = rt.THStack("h_pNScr_all_wCuts", "Inclusive CCnumu Selected Events")
h_pNScr_CCnumu_wCuts, h_pNScr_NCnumu_wCuts, h_pNScr_CCnue_wCuts, h_pNScr_NCnue_wCuts, h_pNScr_ext_wCuts, h_pNScr_data_wCuts = configureStackedHists(h_pNScr_CCnumu_wCuts, h_pNScr_NCnumu_wCuts, h_pNScr_CCnue_wCuts, h_pNScr_NCnue_wCuts, h_pNScr_ext_wCuts, h_pNScr_data_wCuts, "Inclusive CCnumu Selected Events", "reconstructed muon's secondary (neutral parent) score")

h_pCScr_CCnue_wCuts = rt.TH1F("h_pCScr_CCnue_wCuts","Reco Muon's Secondary (Charged Parent) Score for True CCnue Events",procScore_n,procScore_l,procScore_h)
h_pCScr_CCnumu_wCuts = rt.TH1F("h_pCScr_CCnumu_wCuts","Reco Muon's Secondary (Charged Parent) Score for True CCnumu Events",procScore_n,procScore_l,procScore_h)
h_pCScr_NCnumu_wCuts = rt.TH1F("h_pCScr_NCnumu_wCuts","Reco Muon's Secondary (Charged Parent) Score for True NCnumu Events",procScore_n,procScore_l,procScore_h)
h_pCScr_NCnue_wCuts = rt.TH1F("h_pCScr_NCnue_wCuts","Reco Muon's Secondary (Charged Parent) Score for True NCnue Events",procScore_n,procScore_l,procScore_h)
h_pCScr_ext_wCuts = rt.TH1F("h_pCScr_ext_wCuts","Reco Muon's Secondary (Charged Parent) Score for ExtBNB Events",procScore_n,procScore_l,procScore_h)
h_pCScr_data_wCuts = rt.TH1F("h_pCScr_data_wCuts","Reco Muon's Secondary (Charged Parent) Score for BNB Data Events",procScore_n,procScore_l,procScore_h)
h_pCScr_all_wCuts = rt.THStack("h_pCScr_all_wCuts", "Inclusive CCnumu Selected Events")
h_pCScr_CCnumu_wCuts, h_pCScr_NCnumu_wCuts, h_pCScr_CCnue_wCuts, h_pCScr_NCnue_wCuts, h_pCScr_ext_wCuts, h_pCScr_data_wCuts = configureStackedHists(h_pCScr_CCnumu_wCuts, h_pCScr_NCnumu_wCuts, h_pCScr_CCnue_wCuts, h_pCScr_NCnue_wCuts, h_pCScr_ext_wCuts, h_pCScr_data_wCuts, "Inclusive CCnumu Selected Events", "reconstructed muon's secondary (charged parent) score")



h_lowENVsX_CCnumu_wCutsSet1 = rt.TH1F("h_lowENVsX_CCnumu_wCutsSet1","Low Energy (< 400 MeV) Event Count vs. X",26,0,260)
h_lowENVsX_NCnumu_wCutsSet1 = rt.TH1F("h_lowENVsX_NCnumu_wCutsSet1","Low Energy (< 400 MeV) Event Count vs. X",26,0,260)
h_lowENVsX_CCnue_wCutsSet1 = rt.TH1F("h_lowENVsX_CCnue_wCutsSet1","Low Energy (< 400 MeV) Event Count vs. X",26,0,260)
h_lowENVsX_NCnue_wCutsSet1 = rt.TH1F("h_lowENVsX_NCnue_wCutsSet1","Low Energy (< 400 MeV) Event Count vs. X",26,0,260)
h_lowENVsX_ext_wCutsSet1 = rt.TH1F("h_lowENVsX_ext_wCutsSet1","Low Energy (< 400 MeV) Event Count vs. X",26,0,260)
h_lowENVsX_data_wCutsSet1 = rt.TH1F("h_lowENVsX_data_wCutsSet1","Low Energy (< 400 MeV) Event Count vs. X",26,0,260)

h_lowEExVsX_wCutsSet1 = rt.TH1F("h_lowEExVsX_wCutsSet1", "Inclusive CCnumu Low Energy (< 400 MeV) Excess (Cut Subset 1)",26,0,260)
h_lowEExVsX_wCutsSet1.SetLineColor(rt.kBlack)
h_lowEExVsX_wCutsSet1.SetLineWidth(2)
h_lowEExVsX_wCutsSet1.GetXaxis().SetTitle("reconstructed vertex x position (cm)")
h_lowEExVsX_wCutsSet1.GetYaxis().SetTitle("excess predicted events per "+targetPOTstring+" POT")

h_lowENVsX_CCnumu_wCutsSet2 = rt.TH1F("h_lowENVsX_CCnumu_wCutsSet2","Low Energy (< 400 MeV) Event Count vs. X",26,0,260)
h_lowENVsX_NCnumu_wCutsSet2 = rt.TH1F("h_lowENVsX_NCnumu_wCutsSet2","Low Energy (< 400 MeV) Event Count vs. X",26,0,260)
h_lowENVsX_CCnue_wCutsSet2 = rt.TH1F("h_lowENVsX_CCnue_wCutsSet2","Low Energy (< 400 MeV) Event Count vs. X",26,0,260)
h_lowENVsX_NCnue_wCutsSet2 = rt.TH1F("h_lowENVsX_NCnue_wCutsSet2","Low Energy (< 400 MeV) Event Count vs. X",26,0,260)
h_lowENVsX_ext_wCutsSet2 = rt.TH1F("h_lowENVsX_ext_wCutsSet2","Low Energy (< 400 MeV) Event Count vs. X",26,0,260)
h_lowENVsX_data_wCutsSet2 = rt.TH1F("h_lowENVsX_data_wCutsSet2","Low Energy (< 400 MeV) Event Count vs. X",26,0,260)

h_lowEExVsX_wCutsSet2 = rt.TH1F("h_lowEExVsX_wCutsSet2", "Inclusive CCnumu Low Energy (< 400 MeV) Excess (Cut Subset 2)",26,0,260)
h_lowEExVsX_wCutsSet2.SetLineColor(rt.kBlack)
h_lowEExVsX_wCutsSet2.SetLineWidth(2)
h_lowEExVsX_wCutsSet2.GetXaxis().SetTitle("reconstructed vertex x position (cm)")
h_lowEExVsX_wCutsSet2.GetYaxis().SetTitle("excess predicted events per "+targetPOTstring+" POT")


h_lowENVsX_CCnumu_wCutsSet3 = rt.TH1F("h_lowENVsX_CCnumu_wCutsSet3","Low Energy (< 400 MeV) Event Count vs. X",13,0,260)
h_lowENVsX_NCnumu_wCutsSet3 = rt.TH1F("h_lowENVsX_NCnumu_wCutsSet3","Low Energy (< 400 MeV) Event Count vs. X",13,0,260)
h_lowENVsX_CCnue_wCutsSet3 = rt.TH1F("h_lowENVsX_CCnue_wCutsSet3","Low Energy (< 400 MeV) Event Count vs. X",13,0,260)
h_lowENVsX_NCnue_wCutsSet3 = rt.TH1F("h_lowENVsX_NCnue_wCutsSet3","Low Energy (< 400 MeV) Event Count vs. X",13,0,260)
h_lowENVsX_ext_wCutsSet3 = rt.TH1F("h_lowENVsX_ext_wCutsSet3","Low Energy (< 400 MeV) Event Count vs. X",13,0,260)
h_lowENVsX_data_wCutsSet3 = rt.TH1F("h_lowENVsX_data_wCutsSet3","Low Energy (< 400 MeV) Event Count vs. X",13,0,260)

h_lowEExVsX_wCutsSet3 = rt.TH1F("h_lowEExVsX_wCutsSet3", "Inclusive CCnumu Low Energy (< 400 MeV) Excess (Cut Subset 3)",13,0,260)
h_lowEExVsX_wCutsSet3.SetLineColor(rt.kBlack)
h_lowEExVsX_wCutsSet3.SetLineWidth(2)
h_lowEExVsX_wCutsSet3.GetXaxis().SetTitle("reconstructed vertex x position (cm)")
h_lowEExVsX_wCutsSet3.GetYaxis().SetTitle("excess predicted events per "+targetPOTstring+" POT")


h_lowENVsX_CCnumu_wCutsSet4 = rt.TH1F("h_lowENVsX_CCnumu_wCutsSet4","Low Energy (< 400 MeV) Event Count vs. X",13,0,260)
h_lowENVsX_NCnumu_wCutsSet4 = rt.TH1F("h_lowENVsX_NCnumu_wCutsSet4","Low Energy (< 400 MeV) Event Count vs. X",13,0,260)
h_lowENVsX_CCnue_wCutsSet4 = rt.TH1F("h_lowENVsX_CCnue_wCutsSet4","Low Energy (< 400 MeV) Event Count vs. X",13,0,260)
h_lowENVsX_NCnue_wCutsSet4 = rt.TH1F("h_lowENVsX_NCnue_wCutsSet4","Low Energy (< 400 MeV) Event Count vs. X",13,0,260)
h_lowENVsX_ext_wCutsSet4 = rt.TH1F("h_lowENVsX_ext_wCutsSet4","Low Energy (< 400 MeV) Event Count vs. X",13,0,260)
h_lowENVsX_data_wCutsSet4 = rt.TH1F("h_lowENVsX_data_wCutsSet4","Low Energy (< 400 MeV) Event Count vs. X",13,0,260)

h_lowEExVsX_wCutsSet4 = rt.TH1F("h_lowEExVsX_wCutsSet4", "Inclusive CCnumu Low Energy (< 400 MeV) Excess (Cut Subset 4)",13,0,260)
h_lowEExVsX_wCutsSet4.SetLineColor(rt.kBlack)
h_lowEExVsX_wCutsSet4.SetLineWidth(2)
h_lowEExVsX_wCutsSet4.GetXaxis().SetTitle("reconstructed vertex x position (cm)")
h_lowEExVsX_wCutsSet4.GetYaxis().SetTitle("excess predicted events per "+targetPOTstring+" POT")


h_lowENVsX_CCnumu_wCutsSet5 = rt.TH1F("h_lowENVsX_CCnumu_wCutsSet5","Low Energy (< 400 MeV) Event Count vs. X",13,0,260)
h_lowENVsX_NCnumu_wCutsSet5 = rt.TH1F("h_lowENVsX_NCnumu_wCutsSet5","Low Energy (< 400 MeV) Event Count vs. X",13,0,260)
h_lowENVsX_CCnue_wCutsSet5 = rt.TH1F("h_lowENVsX_CCnue_wCutsSet5","Low Energy (< 400 MeV) Event Count vs. X",13,0,260)
h_lowENVsX_NCnue_wCutsSet5 = rt.TH1F("h_lowENVsX_NCnue_wCutsSet5","Low Energy (< 400 MeV) Event Count vs. X",13,0,260)
h_lowENVsX_ext_wCutsSet5 = rt.TH1F("h_lowENVsX_ext_wCutsSet5","Low Energy (< 400 MeV) Event Count vs. X",13,0,260)
h_lowENVsX_data_wCutsSet5 = rt.TH1F("h_lowENVsX_data_wCutsSet5","Low Energy (< 400 MeV) Event Count vs. X",13,0,260)

h_lowEExVsX_wCutsSet5 = rt.TH1F("h_lowEExVsX_wCutsSet5", "Inclusive CCnumu Low Energy (< 400 MeV) Excess (Cut Subset 5)",13,0,260)
h_lowEExVsX_wCutsSet5.SetLineColor(rt.kBlack)
h_lowEExVsX_wCutsSet5.SetLineWidth(2)
h_lowEExVsX_wCutsSet5.GetXaxis().SetTitle("reconstructed vertex x position (cm)")
h_lowEExVsX_wCutsSet5.GetYaxis().SetTitle("excess predicted events per "+targetPOTstring+" POT")


h_lowENVsX_CCnumu_wCutsSet6 = rt.TH1F("h_lowENVsX_CCnumu_wCutsSet6","Low Energy (< 400 MeV) Event Count vs. X",13,0,260)
h_lowENVsX_NCnumu_wCutsSet6 = rt.TH1F("h_lowENVsX_NCnumu_wCutsSet6","Low Energy (< 400 MeV) Event Count vs. X",13,0,260)
h_lowENVsX_CCnue_wCutsSet6 = rt.TH1F("h_lowENVsX_CCnue_wCutsSet6","Low Energy (< 400 MeV) Event Count vs. X",13,0,260)
h_lowENVsX_NCnue_wCutsSet6 = rt.TH1F("h_lowENVsX_NCnue_wCutsSet6","Low Energy (< 400 MeV) Event Count vs. X",13,0,260)
h_lowENVsX_ext_wCutsSet6 = rt.TH1F("h_lowENVsX_ext_wCutsSet6","Low Energy (< 400 MeV) Event Count vs. X",13,0,260)
h_lowENVsX_data_wCutsSet6 = rt.TH1F("h_lowENVsX_data_wCutsSet6","Low Energy (< 400 MeV) Event Count vs. X",13,0,260)

h_lowEExVsX_wCutsSet6 = rt.TH1F("h_lowEExVsX_wCutsSet6", "Inclusive CCnumu Low Energy (< 400 MeV) Excess (Cut Subset 6)",13,0,260)
h_lowEExVsX_wCutsSet6.SetLineColor(rt.kBlack)
h_lowEExVsX_wCutsSet6.SetLineWidth(2)
h_lowEExVsX_wCutsSet6.GetXaxis().SetTitle("reconstructed vertex x position (cm)")
h_lowEExVsX_wCutsSet6.GetYaxis().SetTitle("excess predicted events per "+targetPOTstring+" POT")


h_lowENVsX_CCnumu_wCuts = rt.TH1F("h_lowENVsX_CCnumu_wCuts","Low Energy (< 400 MeV) Event Count vs. X",13,0,260)
h_lowENVsX_NCnumu_wCuts = rt.TH1F("h_lowENVsX_NCnumu_wCuts","Low Energy (< 400 MeV) Event Count vs. X",13,0,260)
h_lowENVsX_CCnue_wCuts = rt.TH1F("h_lowENVsX_CCnue_wCuts","Low Energy (< 400 MeV) Event Count vs. X",13,0,260)
h_lowENVsX_NCnue_wCuts = rt.TH1F("h_lowENVsX_NCnue_wCuts","Low Energy (< 400 MeV) Event Count vs. X",13,0,260)
h_lowENVsX_ext_wCuts = rt.TH1F("h_lowENVsX_ext_wCuts","Low Energy (< 400 MeV) Event Count vs. X",13,0,260)
h_lowENVsX_data_wCuts = rt.TH1F("h_lowENVsX_data_wCuts","Low Energy (< 400 MeV) Event Count vs. X",13,0,260)

h_lowEExVsX_wCuts = rt.TH1F("h_lowEExVsX_wCuts", "Inclusive CCnumu Low Energy (< 400 MeV) Excess",13,0,260)
h_lowEExVsX_wCuts.SetLineColor(rt.kBlack)
h_lowEExVsX_wCuts.SetLineWidth(2)
h_lowEExVsX_wCuts.GetXaxis().SetTitle("reconstructed vertex x position (cm)")
h_lowEExVsX_wCuts.GetYaxis().SetTitle("excess predicted events per "+targetPOTstring+" POT")




n_CCnumu_nCuts = 0.
n_CCnumu_wCuts = 0.
n_NCnumu_nCuts = 0.
n_NCnumu_wCuts = 0.
n_CCnue_nCuts = 0.
n_CCnue_wCuts = 0.
n_NCnue_nCuts = 0.
n_NCnue_wCuts = 0.
n_ext_wCuts = 0.
n_ext_nCuts = 0.
n_data_wCuts = 0.
n_data_nCuts = 0.




def FillNuHistos(h_CCnumu, h_NCnumu, h_NCnue, val, weight, eventType):
  if eventType == 0:
    h_CCnumu.Fill(val, weight)
  if eventType == 1:
    h_NCnumu.Fill(val, weight)
  if eventType == 2:
    h_NCnue.Fill(val, weight)
  return h_CCnumu, h_NCnumu, h_NCnue


print("starting loop over bnb nu overlay events")

for i in range(tnu.GetEntries()):

  tnu.GetEntry(i)

  if args.smallFV:
    trueVtxPos = rt.TVector3(tnu.trueVtxX, tnu.trueVtxY, tnu.trueVtxZ)
    if not isFiducial(trueVtxPos):
      continue

  eventType = -1

  if abs(tnu.trueNuPDG) == 14:
    #CC numu
    if tnu.trueNuCCNC == 0:
      eventType = 0
      n_CCnumu_nCuts += tnu.xsecWeight
      h_nuE_CCnumu_nCuts.Fill(tnu.trueNuE, tnu.xsecWeight)
      h_nuEr_CCnumu_nCuts.Fill(tnu.recoNuE/1000., tnu.xsecWeight)
    #NC numu
    else:
      eventType = 1
      n_NCnumu_nCuts += tnu.xsecWeight
  #NC nue
  if abs(tnu.trueNuPDG) == 12 and tnu.trueNuCCNC == 1:
    eventType = 2
    n_NCnue_nCuts += tnu.xsecWeight

  if eventType < 0:
    continue

  if args.smallFV:
    vtxPos = rt.TVector3(tnu.vtxX, tnu.vtxY, tnu.vtxZ)
    vtxIsFiducial = isFiducial(vtxPos)
  else:
    vtxIsFiducial = (tnu.vtxIsFiducial == 1)
  if tnu.foundVertex == 0 or not vtxIsFiducial: #tnu.vtxIsFiducial != 1:
    continue

  if eventType == 0:
    h_nuE_CCnumu_wCuts1.Fill(tnu.trueNuE, tnu.xsecWeight)

  overflowRecoEVal = tnu.recoNuE/1000.
  if overflowRecoEVal > 2.0 and args.recoEOverflow:
    overflowRecoEVal = 2.001

  h_visE_CCnumu_wCutsSet1, h_visE_NCnumu_wCutsSet1, h_visE_NCnue_wCutsSet1 = FillNuHistos(h_visE_CCnumu_wCutsSet1, h_visE_NCnumu_wCutsSet1, h_visE_NCnue_wCutsSet1, overflowRecoEVal, tnu.xsecWeight, eventType)

  if tnu.recoNuE < 400.:
    h_lowENVsX_CCnumu_wCutsSet1, h_lowENVsX_NCnumu_wCutsSet1, h_lowENVsX_NCnue_wCutsSet1 = FillNuHistos(h_lowENVsX_CCnumu_wCutsSet1, h_lowENVsX_NCnumu_wCutsSet1, h_lowENVsX_NCnue_wCutsSet1, tnu.vtxX, tnu.xsecWeight, eventType)

  h_cosFrac_CCnumu, h_cosFrac_NCnumu, h_cosFrac_NCnue = FillNuHistos(h_cosFrac_CCnumu,
    h_cosFrac_NCnumu, h_cosFrac_NCnue, tnu.vtxFracHitsOnCosmic, tnu.xsecWeight, eventType)

  if tnu.vtxFracHitsOnCosmic >= 1.:
    continue

  if eventType == 0:
    h_nuE_CCnumu_wCuts2.Fill(tnu.trueNuE, tnu.xsecWeight)

  h_visE_CCnumu_wCutsSet2, h_visE_NCnumu_wCutsSet2, h_visE_NCnue_wCutsSet2 = FillNuHistos(h_visE_CCnumu_wCutsSet2, h_visE_NCnumu_wCutsSet2, h_visE_NCnue_wCutsSet2, overflowRecoEVal, tnu.xsecWeight, eventType)

  if tnu.recoNuE < 400.:
    h_lowENVsX_CCnumu_wCutsSet2, h_lowENVsX_NCnumu_wCutsSet2, h_lowENVsX_NCnue_wCutsSet2 = FillNuHistos(h_lowENVsX_CCnumu_wCutsSet2, h_lowENVsX_NCnumu_wCutsSet2, h_lowENVsX_NCnue_wCutsSet2, tnu.vtxX, tnu.xsecWeight, eventType)

  nMu = 0
  selMu_cosThetaBeam = -9.
  selMu_cosThetaGrav = -9.
  selMu_recoE = -99.
  selMu_muScore = -99.
  selMu_piScore = -99.
  selMu_prScore = -99.
  selMu_elScore = -99.
  selMu_phScore = -99.
  selMu_scoreConf = -99.
  selMu_proc = -9
  selMu_primScore = -99.
  selMu_ntrlScore = -99.
  selMu_chgdScore = -99.
  selMu_procConf = -99.

  for iT in range(tnu.nTracks):
    if tnu.trackIsSecondary[iT] == 1 or tnu.trackClassified[iT] != 1:
      continue
    if tnu.trackPID[iT] == 13:
      nMu += 1
      if tnu.trackMuScore[iT] > selMu_muScore:
        selMu_cosThetaBeam = getCosThetaBeamVector(tnu.trackStartDirX[iT], tnu.trackStartDirY[iT], 
         tnu.trackStartDirZ[iT])
        selMu_cosThetaGrav = getCosThetaGravVector(tnu.trackStartDirX[iT], tnu.trackStartDirY[iT], 
         tnu.trackStartDirZ[iT])
        selMu_recoE = tnu.trackRecoE[iT]
        selMu_muScore = tnu.trackMuScore[iT]
        selMu_piScore = tnu.trackPiScore[iT]
        selMu_prScore = tnu.trackPrScore[iT]
        selMu_elScore = tnu.trackElScore[iT]
        selMu_phScore = tnu.trackPhScore[iT]
        selMu_scoreConf = tnu.trackMuScore[iT] - tnu.trackPiScore[iT]
        selMu_proc = tnu.trackProcess[iT]
        selMu_primScore = tnu.trackPrimaryScore[iT]
        selMu_ntrlScore = tnu.trackFromNeutralScore[iT]
        selMu_chgdScore = tnu.trackFromChargedScore[iT]
        selMu_procConf = tnu.trackPrimaryScore[iT] - (tnu.trackFromNeutralScore[iT] + tnu.trackFromChargedScore[iT])/2.

  h_nMu_CCnumu, h_nMu_NCnumu, h_nMu_NCnue = FillNuHistos(h_nMu_CCnumu,
    h_nMu_NCnumu, h_nMu_NCnue, nMu, tnu.xsecWeight, eventType)

  if nMu < 1:
    continue

  h_visE_CCnumu_wCutsSet3, h_visE_NCnumu_wCutsSet3, h_visE_NCnue_wCutsSet3 = FillNuHistos(h_visE_CCnumu_wCutsSet3, h_visE_NCnumu_wCutsSet3, h_visE_NCnue_wCutsSet3, overflowRecoEVal, tnu.xsecWeight, eventType)

  if tnu.recoNuE < 400.:
    h_lowENVsX_CCnumu_wCutsSet3, h_lowENVsX_NCnumu_wCutsSet3, h_lowENVsX_NCnue_wCutsSet3 = FillNuHistos(h_lowENVsX_CCnumu_wCutsSet3, h_lowENVsX_NCnumu_wCutsSet3, h_lowENVsX_NCnue_wCutsSet3, tnu.vtxX, tnu.xsecWeight, eventType)

  h_vtxScore_CCnumu, h_vtxScore_NCnumu, h_vtxScore_NCnue = FillNuHistos(h_vtxScore_CCnumu,
    h_vtxScore_NCnumu, h_vtxScore_NCnue, tnu.vtxScore, tnu.xsecWeight, eventType)

  if tnu.vtxScore < args.vertexScoreCut:
    continue

  h_visE_CCnumu_wCutsSet4, h_visE_NCnumu_wCutsSet4, h_visE_NCnue_wCutsSet4 = FillNuHistos(h_visE_CCnumu_wCutsSet4, h_visE_NCnumu_wCutsSet4, h_visE_NCnue_wCutsSet4, overflowRecoEVal, tnu.xsecWeight, eventType)

  if tnu.recoNuE < 400.:
    h_lowENVsX_CCnumu_wCutsSet4, h_lowENVsX_NCnumu_wCutsSet4, h_lowENVsX_NCnue_wCutsSet4 = FillNuHistos(h_lowENVsX_CCnumu_wCutsSet4, h_lowENVsX_NCnumu_wCutsSet4, h_lowENVsX_NCnue_wCutsSet4, tnu.vtxX, tnu.xsecWeight, eventType)

  h_selMu_cosThetaGrav_CCnumu, h_selMu_cosThetaGrav_NCnumu, h_selMu_cosThetaGrav_NCnue = FillNuHistos(h_selMu_cosThetaGrav_CCnumu, h_selMu_cosThetaGrav_NCnumu, h_selMu_cosThetaGrav_NCnue, selMu_cosThetaGrav, tnu.xsecWeight, eventType)
  h_selMu_cosThetaBeam_CCnumu, h_selMu_cosThetaBeam_NCnumu, h_selMu_cosThetaBeam_NCnue = FillNuHistos(h_selMu_cosThetaBeam_CCnumu, h_selMu_cosThetaBeam_NCnumu, h_selMu_cosThetaBeam_NCnue, selMu_cosThetaBeam, tnu.xsecWeight, eventType)

  if selMu_cosThetaGrav < args.cosThetaGravCut:
    continue

  h_visE_CCnumu_wCutsSet5, h_visE_NCnumu_wCutsSet5, h_visE_NCnue_wCutsSet5 = FillNuHistos(h_visE_CCnumu_wCutsSet5, h_visE_NCnumu_wCutsSet5, h_visE_NCnue_wCutsSet5, overflowRecoEVal, tnu.xsecWeight, eventType)

  if tnu.recoNuE < 400.:
    h_lowENVsX_CCnumu_wCutsSet5, h_lowENVsX_NCnumu_wCutsSet5, h_lowENVsX_NCnue_wCutsSet5 = FillNuHistos(h_lowENVsX_CCnumu_wCutsSet5, h_lowENVsX_NCnumu_wCutsSet5, h_lowENVsX_NCnue_wCutsSet5, tnu.vtxX, tnu.xsecWeight, eventType)

  h_selMu_proc_CCnumu, h_selMu_proc_NCnumu, h_selMu_proc_NCnue = FillNuHistos(h_selMu_proc_CCnumu,
    h_selMu_proc_NCnumu, h_selMu_proc_NCnue, selMu_proc, tnu.xsecWeight, eventType)

  if selMu_proc != 0 and args.procCuts:
    continue

  h_visE_CCnumu_wCutsSet6, h_visE_NCnumu_wCutsSet6, h_visE_NCnue_wCutsSet6 = FillNuHistos(h_visE_CCnumu_wCutsSet6, h_visE_NCnumu_wCutsSet6, h_visE_NCnue_wCutsSet6, overflowRecoEVal, tnu.xsecWeight, eventType)

  if tnu.recoNuE < 400.:
    h_lowENVsX_CCnumu_wCutsSet6, h_lowENVsX_NCnumu_wCutsSet6, h_lowENVsX_NCnue_wCutsSet6 = FillNuHistos(h_lowENVsX_CCnumu_wCutsSet6, h_lowENVsX_NCnumu_wCutsSet6, h_lowENVsX_NCnue_wCutsSet6, tnu.vtxX, tnu.xsecWeight, eventType)

  h_selMu_primScore_CCnumu, h_selMu_primScore_NCnumu, h_selMu_primScore_NCnue = FillNuHistos(h_selMu_primScore_CCnumu,
    h_selMu_primScore_NCnumu, h_selMu_primScore_NCnue, selMu_primScore, tnu.xsecWeight, eventType)
  h_selMu_ntrlScore_CCnumu, h_selMu_ntrlScore_NCnumu, h_selMu_ntrlScore_NCnue = FillNuHistos(h_selMu_ntrlScore_CCnumu,
    h_selMu_ntrlScore_NCnumu, h_selMu_ntrlScore_NCnue, selMu_ntrlScore, tnu.xsecWeight, eventType)
  h_selMu_chgdScore_CCnumu, h_selMu_chgdScore_NCnumu, h_selMu_chgdScore_NCnue = FillNuHistos(h_selMu_chgdScore_CCnumu,
    h_selMu_chgdScore_NCnumu, h_selMu_chgdScore_NCnue, selMu_chgdScore, tnu.xsecWeight, eventType)
  h_selMu_procConf_CCnumu, h_selMu_procConf_NCnumu, h_selMu_procConf_NCnue = FillNuHistos(h_selMu_procConf_CCnumu,
    h_selMu_procConf_NCnumu, h_selMu_procConf_NCnue, selMu_procConf, tnu.xsecWeight, eventType)

  if selMu_ntrlScore > args.fromNeutralScoreCut and args.procCuts:
    continue

  if tnu.recoNuE < 400.:
    h_lowENVsX_CCnumu_wCuts, h_lowENVsX_NCnumu_wCuts, h_lowENVsX_NCnue_wCuts = FillNuHistos(h_lowENVsX_CCnumu_wCuts, h_lowENVsX_NCnumu_wCuts, h_lowENVsX_NCnue_wCuts, tnu.vtxX, tnu.xsecWeight, eventType)

  vtxPos = rt.TVector3(tnu.vtxX, tnu.vtxY, tnu.vtxZ)
  cosKPDist = 999.
  for iKP in range(tnu.nKeypoints):
    kpPos = rt.TVector3(tnu.kpMaxPosX[iKP], tnu.kpMaxPosY[iKP], tnu.kpMaxPosZ[iKP])
    distToVtx = getDistance(kpPos, vtxPos)
    if tnu.kpFilterType[iKP] == 1 and distToVtx < cosKPDist:
      cosKPDist = distToVtx

  pcaCosThetaBeam = getCosThetaBeamVector(tnu.eventPCAxis0[0], tnu.eventPCAxis0[1], tnu.eventPCAxis0[2])
  pcaCosThetaGrav = getCosThetaGravVector(tnu.eventPCAxis0[0], tnu.eventPCAxis0[1], tnu.eventPCAxis0[2])
  pcaEVsum = tnu.eventPCEigenVals[0] + tnu.eventPCEigenVals[1] + tnu.eventPCEigenVals[2]
  pcaEVRatio = tnu.eventPCEigenVals[0]/pcaEVsum if (pcaEVsum > 0.) else 0.

  h_cosKPDist_CCnumu, h_cosKPDist_NCnumu, h_cosKPDist_NCnue = FillNuHistos(h_cosKPDist_CCnumu,
    h_cosKPDist_NCnumu, h_cosKPDist_NCnue, cosKPDist, tnu.xsecWeight, eventType)
  h_pcaEVRatio_CCnumu, h_pcaEVRatio_NCnumu, h_pcaEVRatio_NCnue = FillNuHistos(h_pcaEVRatio_CCnumu,
    h_pcaEVRatio_NCnumu, h_pcaEVRatio_NCnue, pcaEVRatio, tnu.xsecWeight, eventType)
  h_pcaCosThetaBeam_CCnumu, h_pcaCosThetaBeam_NCnumu, h_pcaCosThetaBeam_NCnue = FillNuHistos(h_pcaCosThetaBeam_CCnumu,
    h_pcaCosThetaBeam_NCnumu, h_pcaCosThetaBeam_NCnue, pcaCosThetaBeam, tnu.xsecWeight, eventType)
  h_pcaCosThetaGrav_CCnumu, h_pcaCosThetaGrav_NCnumu, h_pcaCosThetaGrav_NCnue = FillNuHistos(h_pcaCosThetaGrav_CCnumu,
    h_pcaCosThetaGrav_NCnumu, h_pcaCosThetaGrav_NCnue, pcaCosThetaGrav, tnu.xsecWeight, eventType)

  h_selMu_muScore_CCnumu, h_selMu_muScore_NCnumu, h_selMu_muScore_NCnue = FillNuHistos(h_selMu_muScore_CCnumu,
    h_selMu_muScore_NCnumu, h_selMu_muScore_NCnue, selMu_muScore, tnu.xsecWeight, eventType)
  h_selMu_piScore_CCnumu, h_selMu_piScore_NCnumu, h_selMu_piScore_NCnue = FillNuHistos(h_selMu_piScore_CCnumu,
    h_selMu_piScore_NCnumu, h_selMu_piScore_NCnue, selMu_piScore, tnu.xsecWeight, eventType)
  h_selMu_scoreConf_CCnumu, h_selMu_scoreConf_NCnumu, h_selMu_scoreConf_NCnue = FillNuHistos(h_selMu_scoreConf_CCnumu,
    h_selMu_scoreConf_NCnumu, h_selMu_scoreConf_NCnue, selMu_scoreConf, tnu.xsecWeight, eventType)

  if args.write_ntuples:
    tnu_trimmed.Fill()

  recoMuP = sqrt( (selMu_recoE + 105.66)**2 - 105.66**2 )/1000.
  if eventType == 0:
    n_CCnumu_wCuts += tnu.xsecWeight
    h_nuE_CCnumu_wCuts.Fill(tnu.trueNuE, tnu.xsecWeight)
    h_nuE_CCnumu_wCutsAll.Fill(tnu.trueNuE, tnu.xsecWeight)
    h_nuEr_CCnumu_wCuts.Fill(tnu.recoNuE/1000., tnu.xsecWeight)
    if args.recoEOverflow and tnu.recoNuE/1000. > 2.0:
      h_visE_CCnumu_wCuts.Fill(2.001, tnu.xsecWeight)
    else:
      h_visE_CCnumu_wCuts.Fill(tnu.recoNuE/1000., tnu.xsecWeight)
    if args.recoEOverflow and recoMuP > 2.0:
      h_lepP_CCnumu_wCuts.Fill(2.001, tnu.xsecWeight)
    else:
      h_lepP_CCnumu_wCuts.Fill(recoMuP, tnu.xsecWeight)
    h_cosTheta_CCnumu_wCuts.Fill(selMu_cosThetaBeam, tnu.xsecWeight)
    h_muScr_CCnumu_wCuts.Fill(selMu_muScore, tnu.xsecWeight)
    h_piScr_CCnumu_wCuts.Fill(selMu_piScore, tnu.xsecWeight)
    h_prScr_CCnumu_wCuts.Fill(selMu_prScore, tnu.xsecWeight)
    h_elScr_CCnumu_wCuts.Fill(selMu_elScore, tnu.xsecWeight)
    h_phScr_CCnumu_wCuts.Fill(selMu_phScore, tnu.xsecWeight)
    h_pPScr_CCnumu_wCuts.Fill(selMu_primScore, tnu.xsecWeight)
    h_pNScr_CCnumu_wCuts.Fill(selMu_ntrlScore, tnu.xsecWeight)
    h_pCScr_CCnumu_wCuts.Fill(selMu_chgdScore, tnu.xsecWeight)
  if eventType == 1:
    n_NCnumu_wCuts += tnu.xsecWeight
    h_nuE_NCnumu_wCuts.Fill(tnu.trueNuE, tnu.xsecWeight)
    h_nuEr_NCnumu_wCuts.Fill(tnu.recoNuE/1000., tnu.xsecWeight)
    if args.recoEOverflow and tnu.recoNuE/1000. > 2.0:
      h_visE_NCnumu_wCuts.Fill(2.001, tnu.xsecWeight)
    else:
      h_visE_NCnumu_wCuts.Fill(tnu.recoNuE/1000., tnu.xsecWeight)
    if args.recoEOverflow and recoMuP > 2.0:
      h_lepP_NCnumu_wCuts.Fill(2.001, tnu.xsecWeight)
    else:
      h_lepP_NCnumu_wCuts.Fill(recoMuP, tnu.xsecWeight)
    h_cosTheta_NCnumu_wCuts.Fill(selMu_cosThetaBeam, tnu.xsecWeight)
    h_muScr_NCnumu_wCuts.Fill(selMu_muScore, tnu.xsecWeight)
    h_piScr_NCnumu_wCuts.Fill(selMu_piScore, tnu.xsecWeight)
    h_prScr_NCnumu_wCuts.Fill(selMu_prScore, tnu.xsecWeight)
    h_elScr_NCnumu_wCuts.Fill(selMu_elScore, tnu.xsecWeight)
    h_phScr_NCnumu_wCuts.Fill(selMu_phScore, tnu.xsecWeight)
    h_pPScr_NCnumu_wCuts.Fill(selMu_primScore, tnu.xsecWeight)
    h_pNScr_NCnumu_wCuts.Fill(selMu_ntrlScore, tnu.xsecWeight)
    h_pCScr_NCnumu_wCuts.Fill(selMu_chgdScore, tnu.xsecWeight)
  if eventType == 2:
    n_NCnue_wCuts += tnu.xsecWeight
    h_nuE_NCnue_wCuts.Fill(tnu.trueNuE, tnu.xsecWeight)
    h_nuEr_NCnue_wCuts.Fill(tnu.recoNuE/1000., tnu.xsecWeight)
    if args.recoEOverflow and tnu.recoNuE/1000. > 2.0:
      h_visE_NCnue_wCuts.Fill(2.001, tnu.xsecWeight)
    else:
      h_visE_NCnue_wCuts.Fill(tnu.recoNuE/1000., tnu.xsecWeight)
    if args.recoEOverflow and recoMuP > 2.0:
      h_lepP_NCnue_wCuts.Fill(2.001, tnu.xsecWeight)
    else:
      h_lepP_NCnue_wCuts.Fill(recoMuP, tnu.xsecWeight)
    h_cosTheta_NCnue_wCuts.Fill(selMu_cosThetaBeam, tnu.xsecWeight)
    h_muScr_NCnue_wCuts.Fill(selMu_muScore, tnu.xsecWeight)
    h_piScr_NCnue_wCuts.Fill(selMu_piScore, tnu.xsecWeight)
    h_prScr_NCnue_wCuts.Fill(selMu_prScore, tnu.xsecWeight)
    h_elScr_NCnue_wCuts.Fill(selMu_elScore, tnu.xsecWeight)
    h_phScr_NCnue_wCuts.Fill(selMu_phScore, tnu.xsecWeight)
    h_pPScr_NCnue_wCuts.Fill(selMu_primScore, tnu.xsecWeight)
    h_pNScr_NCnue_wCuts.Fill(selMu_ntrlScore, tnu.xsecWeight)
    h_pCScr_NCnue_wCuts.Fill(selMu_chgdScore, tnu.xsecWeight)



print("starting loop over bnb intrinsic nue overlay events")

for i in range(tnue.GetEntries()):

  tnue.GetEntry(i)

  if abs(tnue.trueNuPDG) != 12 or tnue.trueNuCCNC != 0:
    continue

  if args.smallFV:
    trueVtxPos = rt.TVector3(tnue.trueVtxX, tnue.trueVtxY, tnue.trueVtxZ)
    if not isFiducial(trueVtxPos):
      continue

  n_CCnue_nCuts += tnue.xsecWeight

  if args.smallFV:
    vtxPos = rt.TVector3(tnue.vtxX, tnue.vtxY, tnue.vtxZ)
    vtxIsFiducial = isFiducial(vtxPos)
  else:
    vtxIsFiducial = (tnue.vtxIsFiducial == 1)
  if tnue.foundVertex == 0 or not vtxIsFiducial: #tnue.vtxIsFiducial != 1:
    continue

  overflowRecoEVal = tnue.recoNuE/1000.
  if overflowRecoEVal > 2.0 and args.recoEOverflow:
    overflowRecoEVal = 2.001

  h_visE_CCnue_wCutsSet1.Fill(overflowRecoEVal, tnue.xsecWeight)

  if tnue.recoNuE < 400.:
    h_lowENVsX_CCnue_wCutsSet1.Fill(tnue.vtxX, tnue.xsecWeight)

  h_cosFrac_CCnue.Fill(tnue.vtxFracHitsOnCosmic, tnue.xsecWeight)

  if tnue.vtxFracHitsOnCosmic >= 1.:
    continue

  h_visE_CCnue_wCutsSet2.Fill(overflowRecoEVal, tnue.xsecWeight)

  if tnue.recoNuE < 400.:
    h_lowENVsX_CCnue_wCutsSet2.Fill(tnue.vtxX, tnue.xsecWeight)

  nMu = 0
  selMu_cosThetaBeam = -9.
  selMu_cosThetaGrav = -9.
  selMu_recoE = -99.
  selMu_muScore = -99.
  selMu_piScore = -99.
  selMu_prScore = -99.
  selMu_elScore = -99.
  selMu_phScore = -99.
  selMu_scoreConf = -99.
  selMu_proc = -9
  selMu_primScore = -99.
  selMu_ntrlScore = -99.
  selMu_chgdScore = -99.
  selMu_procConf = -99.

  for iT in range(tnue.nTracks):
    if tnue.trackIsSecondary[iT] == 1 or tnue.trackClassified[iT] != 1:
      continue
    if tnue.trackPID[iT] == 13:
      nMu += 1
      if tnue.trackMuScore[iT] > selMu_muScore:
        selMu_cosThetaBeam = getCosThetaBeamVector(tnue.trackStartDirX[iT], tnue.trackStartDirY[iT], 
         tnue.trackStartDirZ[iT])
        selMu_cosThetaGrav = getCosThetaGravVector(tnue.trackStartDirX[iT], tnue.trackStartDirY[iT], 
         tnue.trackStartDirZ[iT])
        selMu_recoE = tnue.trackRecoE[iT]
        selMu_muScore = tnue.trackMuScore[iT]
        selMu_piScore = tnue.trackPiScore[iT]
        selMu_prScore = tnue.trackPrScore[iT]
        selMu_elScore = tnue.trackElScore[iT]
        selMu_phScore = tnue.trackPhScore[iT]
        selMu_scoreConf = tnue.trackMuScore[iT] - tnue.trackPiScore[iT]
        selMu_proc = tnue.trackProcess[iT]
        selMu_primScore = tnue.trackPrimaryScore[iT]
        selMu_ntrlScore = tnue.trackFromNeutralScore[iT]
        selMu_chgdScore = tnue.trackFromChargedScore[iT]
        selMu_procConf = tnue.trackPrimaryScore[iT] - (tnue.trackFromNeutralScore[iT] + tnue.trackFromChargedScore[iT])/2.

  h_nMu_CCnue.Fill(nMu, tnue.xsecWeight)

  if nMu < 1:
    continue

  h_visE_CCnue_wCutsSet3.Fill(overflowRecoEVal, tnue.xsecWeight)

  if tnue.recoNuE < 400.:
    h_lowENVsX_CCnue_wCutsSet3.Fill(tnue.vtxX, tnue.xsecWeight)

  h_vtxScore_CCnue.Fill(tnue.vtxScore, tnue.xsecWeight)

  if tnue.vtxScore < args.vertexScoreCut:
    continue

  h_visE_CCnue_wCutsSet4.Fill(overflowRecoEVal, tnue.xsecWeight)

  if tnue.recoNuE < 400.:
    h_lowENVsX_CCnue_wCutsSet4.Fill(tnue.vtxX, tnue.xsecWeight)

  h_selMu_cosThetaGrav_CCnue.Fill(selMu_cosThetaGrav, tnue.xsecWeight)
  h_selMu_cosThetaBeam_CCnue.Fill(selMu_cosThetaBeam, tnue.xsecWeight)

  if selMu_cosThetaGrav < args.cosThetaGravCut:
    continue

  h_visE_CCnue_wCutsSet5.Fill(overflowRecoEVal, tnue.xsecWeight)

  if tnue.recoNuE < 400.:
    h_lowENVsX_CCnue_wCutsSet5.Fill(tnue.vtxX, tnue.xsecWeight)

  h_selMu_proc_CCnue.Fill(selMu_proc, tnue.xsecWeight)

  if selMu_proc != 0 and args.procCuts:
    continue

  h_visE_CCnue_wCutsSet6.Fill(overflowRecoEVal, tnue.xsecWeight)

  if tnue.recoNuE < 400.:
    h_lowENVsX_CCnue_wCutsSet6.Fill(tnue.vtxX, tnue.xsecWeight)

  h_selMu_primScore_CCnue.Fill(selMu_primScore, tnue.xsecWeight)
  h_selMu_ntrlScore_CCnue.Fill(selMu_ntrlScore, tnue.xsecWeight)
  h_selMu_chgdScore_CCnue.Fill(selMu_chgdScore, tnue.xsecWeight)
  h_selMu_procConf_CCnue.Fill(selMu_procConf, tnue.xsecWeight)

  if selMu_ntrlScore > args.fromNeutralScoreCut and args.procCuts:
    continue

  if tnue.recoNuE < 400.:
    h_lowENVsX_CCnue_wCuts.Fill(tnue.vtxX, tnue.xsecWeight)

  vtxPos = rt.TVector3(tnue.vtxX, tnue.vtxY, tnue.vtxZ)
  cosKPDist = 999.
  for iKP in range(tnue.nKeypoints):
    kpPos = rt.TVector3(tnue.kpMaxPosX[iKP], tnue.kpMaxPosY[iKP], tnue.kpMaxPosZ[iKP])
    distToVtx = getDistance(kpPos, vtxPos)
    if tnue.kpFilterType[iKP] == 1 and distToVtx < cosKPDist:
      cosKPDist = distToVtx

  pcaCosThetaBeam = getCosThetaBeamVector(tnue.eventPCAxis0[0], tnue.eventPCAxis0[1], tnue.eventPCAxis0[2])
  pcaCosThetaGrav = getCosThetaGravVector(tnue.eventPCAxis0[0], tnue.eventPCAxis0[1], tnue.eventPCAxis0[2])
  pcaEVsum = tnue.eventPCEigenVals[0] + tnue.eventPCEigenVals[1] + tnue.eventPCEigenVals[2]
  pcaEVRatio = tnue.eventPCEigenVals[0]/pcaEVsum if (pcaEVsum > 0.) else 0.

  h_cosKPDist_CCnue.Fill(cosKPDist, tnue.xsecWeight)
  h_pcaEVRatio_CCnue.Fill(pcaEVRatio, tnue.xsecWeight)
  h_pcaCosThetaBeam_CCnue.Fill(pcaCosThetaBeam, tnue.xsecWeight)
  h_pcaCosThetaGrav_CCnue.Fill(pcaCosThetaGrav, tnue.xsecWeight)

  h_selMu_muScore_CCnue.Fill(selMu_muScore, tnue.xsecWeight)
  h_selMu_piScore_CCnue.Fill(selMu_piScore, tnue.xsecWeight)
  h_selMu_scoreConf_CCnue.Fill(selMu_scoreConf, tnue.xsecWeight)

  if args.write_ntuples:
    tnue_trimmed.Fill()

  n_CCnue_wCuts += tnue.xsecWeight
  h_nuE_CCnue_wCuts.Fill(tnue.trueNuE, tnue.xsecWeight)
  h_nuEr_CCnue_wCuts.Fill(tnue.recoNuE/1000., tnue.xsecWeight)
  if args.recoEOverflow and tnue.recoNuE/1000. > 2.0:
    h_visE_CCnue_wCuts.Fill(2.001, tnue.xsecWeight)
  else:
    h_visE_CCnue_wCuts.Fill(tnue.recoNuE/1000., tnue.xsecWeight)
  recoMuP = sqrt( (selMu_recoE + 105.66)**2 - 105.66**2 )/1000.
  if args.recoEOverflow and recoMuP > 2.0:
    h_lepP_CCnue_wCuts.Fill(2.001, tnue.xsecWeight)
  else:
    h_lepP_CCnue_wCuts.Fill(recoMuP, tnue.xsecWeight)
  h_cosTheta_CCnue_wCuts.Fill(selMu_cosThetaBeam, tnue.xsecWeight)
  h_muScr_CCnue_wCuts.Fill(selMu_muScore, tnue.xsecWeight)
  h_piScr_CCnue_wCuts.Fill(selMu_piScore, tnue.xsecWeight)
  h_prScr_CCnue_wCuts.Fill(selMu_prScore, tnue.xsecWeight)
  h_elScr_CCnue_wCuts.Fill(selMu_elScore, tnue.xsecWeight)
  h_phScr_CCnue_wCuts.Fill(selMu_phScore, tnue.xsecWeight)
  h_pPScr_CCnue_wCuts.Fill(selMu_primScore, tnue.xsecWeight)
  h_pNScr_CCnue_wCuts.Fill(selMu_ntrlScore, tnue.xsecWeight)
  h_pCScr_CCnue_wCuts.Fill(selMu_chgdScore, tnue.xsecWeight)



print("starting loop over extBNB events")

for i in range(text.GetEntries()):

  text.GetEntry(i)

  n_ext_nCuts += 1.

  if args.smallFV:
    vtxPos = rt.TVector3(text.vtxX, text.vtxY, text.vtxZ)
    vtxIsFiducial = isFiducial(vtxPos)
  else:
    vtxIsFiducial = (text.vtxIsFiducial == 1)
  if text.foundVertex == 0 or not vtxIsFiducial: #text.vtxIsFiducial != 1:
    continue

  overflowRecoEVal = text.recoNuE/1000.
  if overflowRecoEVal > 2.0 and args.recoEOverflow:
    overflowRecoEVal = 2.001

  h_visE_ext_wCutsSet1.Fill(overflowRecoEVal)

  if text.recoNuE < 400.:
    h_lowENVsX_ext_wCutsSet1.Fill(text.vtxX)

  h_cosFrac_ext.Fill(text.vtxFracHitsOnCosmic)

  if text.vtxFracHitsOnCosmic >= 1.:
    continue

  h_visE_ext_wCutsSet2.Fill(overflowRecoEVal)

  if text.recoNuE < 400.:
    h_lowENVsX_ext_wCutsSet2.Fill(text.vtxX)

  nMu = 0
  selMu_cosThetaBeam = -9.
  selMu_cosThetaGrav = -9.
  selMu_recoE = -99.
  selMu_muScore = -99.
  selMu_piScore = -99.
  selMu_prScore = -99.
  selMu_elScore = -99.
  selMu_phScore = -99.
  selMu_scoreConf = -99.
  selMu_proc = -9
  selMu_primScore = -99.
  selMu_ntrlScore = -99.
  selMu_chgdScore = -99.
  selMu_procConf = -99.
  selMu_idx = -1

  for iT in range(text.nTracks):
    if text.trackIsSecondary[iT] == 1 or text.trackClassified[iT] != 1:
      continue
    if text.trackPID[iT] == 13:
      nMu += 1
      if text.trackMuScore[iT] > selMu_muScore:
        selMu_cosThetaBeam = getCosThetaBeamVector(text.trackStartDirX[iT], text.trackStartDirY[iT], 
         text.trackStartDirZ[iT])
        selMu_cosThetaGrav = getCosThetaGravVector(text.trackStartDirX[iT], text.trackStartDirY[iT], 
         text.trackStartDirZ[iT])
        selMu_recoE = text.trackRecoE[iT]
        selMu_muScore = text.trackMuScore[iT]
        selMu_piScore = text.trackPiScore[iT]
        selMu_prScore = text.trackPrScore[iT]
        selMu_elScore = text.trackElScore[iT]
        selMu_phScore = text.trackPhScore[iT]
        selMu_scoreConf = text.trackMuScore[iT] - text.trackPiScore[iT]
        selMu_proc = text.trackProcess[iT]
        selMu_primScore = text.trackPrimaryScore[iT]
        selMu_ntrlScore = text.trackFromNeutralScore[iT]
        selMu_chgdScore = text.trackFromChargedScore[iT]
        selMu_procConf = text.trackPrimaryScore[iT] - (text.trackFromNeutralScore[iT] + text.trackFromChargedScore[iT])/2.
        selMu_idx = iT

  h_nMu_ext.Fill(nMu)

  if nMu < 1:
    continue

  h_visE_ext_wCutsSet3.Fill(overflowRecoEVal)

  if text.recoNuE < 400.:
    h_lowENVsX_ext_wCutsSet3.Fill(text.vtxX)

  h_vtxScore_ext.Fill(text.vtxScore)

  if text.vtxScore < args.vertexScoreCut:
    continue

  h_visE_ext_wCutsSet4.Fill(overflowRecoEVal)

  if text.recoNuE < 400.:
    h_lowENVsX_ext_wCutsSet4.Fill(text.vtxX)

  h_selMu_cosThetaGrav_ext.Fill(selMu_cosThetaGrav)
  h_selMu_cosThetaBeam_ext.Fill(selMu_cosThetaBeam)

  if selMu_cosThetaGrav < args.cosThetaGravCut:
    continue

  h_visE_ext_wCutsSet5.Fill(overflowRecoEVal)

  if text.recoNuE < 400.:
    h_lowENVsX_ext_wCutsSet5.Fill(text.vtxX)

  h_selMu_proc_ext.Fill(selMu_proc)

  if selMu_proc != 0 and args.procCuts:
    continue

  h_visE_ext_wCutsSet6.Fill(overflowRecoEVal)

  if text.recoNuE < 400.:
    h_lowENVsX_ext_wCutsSet6.Fill(text.vtxX)

  h_selMu_primScore_ext.Fill(selMu_primScore)
  h_selMu_ntrlScore_ext.Fill(selMu_ntrlScore)
  h_selMu_chgdScore_ext.Fill(selMu_chgdScore)
  h_selMu_procConf_ext.Fill(selMu_procConf)

  if selMu_ntrlScore > args.fromNeutralScoreCut and args.procCuts:
    continue

  if text.recoNuE < 400.:
    h_lowENVsX_ext_wCuts.Fill(text.vtxX)

  vtxPos = rt.TVector3(text.vtxX, text.vtxY, text.vtxZ)
  cosKPDist = 999.
  for iKP in range(text.nKeypoints):
    kpPos = rt.TVector3(text.kpMaxPosX[iKP], text.kpMaxPosY[iKP], text.kpMaxPosZ[iKP])
    distToVtx = getDistance(kpPos, vtxPos)
    if text.kpFilterType[iKP] == 1 and distToVtx < cosKPDist:
      cosKPDist = distToVtx

  pcaCosThetaBeam = getCosThetaBeamVector(text.eventPCAxis0[0], text.eventPCAxis0[1], text.eventPCAxis0[2])
  pcaCosThetaGrav = getCosThetaGravVector(text.eventPCAxis0[0], text.eventPCAxis0[1], text.eventPCAxis0[2])
  pcaEVsum = text.eventPCEigenVals[0] + text.eventPCEigenVals[1] + text.eventPCEigenVals[2]
  pcaEVRatio = text.eventPCEigenVals[0]/pcaEVsum if (pcaEVsum > 0.) else 0.

  h_cosKPDist_ext.Fill(cosKPDist)
  h_pcaEVRatio_ext.Fill(pcaEVRatio)
  h_pcaCosThetaBeam_ext.Fill(pcaCosThetaBeam)
  h_pcaCosThetaGrav_ext.Fill(pcaCosThetaGrav)

  h_selMu_muScore_ext.Fill(selMu_muScore)
  h_selMu_piScore_ext.Fill(selMu_piScore)
  h_selMu_scoreConf_ext.Fill(selMu_scoreConf)

  if args.write_ntuples:
    text_trimmed.Fill()

  n_ext_wCuts += 1.
  h_nuEr_ext_wCuts.Fill(text.recoNuE/1000.)
  if args.recoEOverflow and text.recoNuE/1000. > 2.0:
    h_visE_ext_wCuts.Fill(2.001)
  else:
    h_visE_ext_wCuts.Fill(text.recoNuE/1000.)
  recoMuP = sqrt( (selMu_recoE + 105.66)**2 - 105.66**2 )/1000.
  if args.recoEOverflow and recoMuP > 2.0:
    h_lepP_ext_wCuts.Fill(2.001)
  else:
    h_lepP_ext_wCuts.Fill(recoMuP)
  h_cosTheta_ext_wCuts.Fill(selMu_cosThetaBeam)
  h_muScr_ext_wCuts.Fill(selMu_muScore)
  h_piScr_ext_wCuts.Fill(selMu_piScore)
  h_prScr_ext_wCuts.Fill(selMu_prScore)
  h_elScr_ext_wCuts.Fill(selMu_elScore)
  h_phScr_ext_wCuts.Fill(selMu_phScore)
  h_pPScr_ext_wCuts.Fill(selMu_primScore)
  h_pNScr_ext_wCuts.Fill(selMu_ntrlScore)
  h_pCScr_ext_wCuts.Fill(selMu_chgdScore)

  #if selMu_cosThetaGrav < -0.98:
  #  print(text.fileid, text.run, text.subrun, text.event, selMu_idx, text.trackRecoE[selMu_idx])
  #  print(text.trackStartDirX[selMu_idx], text.trackStartDirY[selMu_idx], text.trackStartDirZ[selMu_idx])



print("starting loop over data events")

for i in range(tdata.GetEntries()):

  tdata.GetEntry(i)

  n_data_nCuts += 1.

  if args.smallFV:
    vtxPos = rt.TVector3(tdata.vtxX, tdata.vtxY, tdata.vtxZ)
    vtxIsFiducial = isFiducial(vtxPos)
  else:
    vtxIsFiducial = (tdata.vtxIsFiducial == 1)
  if tdata.foundVertex == 0 or not vtxIsFiducial: #tdata.vtxIsFiducial != 1:
    continue

  overflowRecoEVal = tdata.recoNuE/1000.
  if overflowRecoEVal > 2.0 and args.recoEOverflow:
    overflowRecoEVal = 2.001

  h_visE_data_wCutsSet1.Fill(overflowRecoEVal)

  if tdata.recoNuE < 400.:
    h_lowENVsX_data_wCutsSet1.Fill(tdata.vtxX)

  if tdata.vtxFracHitsOnCosmic >= 1.:
    continue

  h_visE_data_wCutsSet2.Fill(overflowRecoEVal)

  if tdata.recoNuE < 400.:
    h_lowENVsX_data_wCutsSet2.Fill(tdata.vtxX)

  nMu = 0
  selMu_cosThetaBeam = -9.
  selMu_cosThetaGrav = -9.
  selMu_recoE = -99.
  selMu_muScore = -99.
  selMu_piScore = -99.
  selMu_prScore = -99.
  selMu_elScore = -99.
  selMu_phScore = -99.
  selMu_proc = -9
  selMu_primScore = -99.
  selMu_ntrlScore = -99.
  selMu_chgdScore = -99.

  for iT in range(tdata.nTracks):
    if tdata.trackIsSecondary[iT] == 1 or tdata.trackClassified[iT] != 1:
      continue
    if tdata.trackPID[iT] == 13:
      nMu += 1
      if tdata.trackMuScore[iT] > selMu_muScore:
        selMu_cosThetaBeam = getCosThetaBeamVector(tdata.trackStartDirX[iT], tdata.trackStartDirY[iT], 
         tdata.trackStartDirZ[iT])
        selMu_cosThetaGrav = getCosThetaGravVector(tdata.trackStartDirX[iT], tdata.trackStartDirY[iT], 
         tdata.trackStartDirZ[iT])
        selMu_recoE = tdata.trackRecoE[iT]
        selMu_muScore = tdata.trackMuScore[iT]
        selMu_piScore = tdata.trackPiScore[iT]
        selMu_prScore = tdata.trackPrScore[iT]
        selMu_elScore = tdata.trackElScore[iT]
        selMu_phScore = tdata.trackPhScore[iT]
        selMu_proc = tdata.trackProcess[iT]
        selMu_primScore = tdata.trackPrimaryScore[iT]
        selMu_ntrlScore = tdata.trackFromNeutralScore[iT]
        selMu_chgdScore = tdata.trackFromChargedScore[iT]

  if nMu < 1:
    continue

  h_visE_data_wCutsSet3.Fill(overflowRecoEVal)

  if tdata.recoNuE < 400.:
    h_lowENVsX_data_wCutsSet3.Fill(tdata.vtxX)

  if tdata.vtxScore < args.vertexScoreCut:
    continue

  h_visE_data_wCutsSet4.Fill(overflowRecoEVal)

  if tdata.recoNuE < 400.:
    h_lowENVsX_data_wCutsSet4.Fill(tdata.vtxX)

  if selMu_cosThetaGrav < args.cosThetaGravCut:
    continue

  h_visE_data_wCutsSet5.Fill(overflowRecoEVal)

  if tdata.recoNuE < 400.:
    h_lowENVsX_data_wCutsSet5.Fill(tdata.vtxX)

  if selMu_proc != 0 and args.procCuts:
    continue

  h_visE_data_wCutsSet6.Fill(overflowRecoEVal)

  if tdata.recoNuE < 400.:
    h_lowENVsX_data_wCutsSet6.Fill(tdata.vtxX)

  if selMu_ntrlScore > args.fromNeutralScoreCut and args.procCuts:
    continue

  if tdata.recoNuE < 400.:
    h_lowENVsX_data_wCuts.Fill(tdata.vtxX)

  if args.write_ntuples:
    tdata_trimmed.Fill()

  n_data_wCuts += 1.
  if args.recoEOverflow and tdata.recoNuE/1000. > 2.0:
    h_visE_data_wCuts.Fill(2.001)
  else:
    h_visE_data_wCuts.Fill(tdata.recoNuE/1000.)
  recoMuP = sqrt( (selMu_recoE + 105.66)**2 - 105.66**2 )/1000.
  if args.recoEOverflow and recoMuP > 2.0:
    h_lepP_data_wCuts.Fill(2.001)
  else:
    h_lepP_data_wCuts.Fill(recoMuP)
  h_cosTheta_data_wCuts.Fill(selMu_cosThetaBeam)
  h_muScr_data_wCuts.Fill(selMu_muScore)
  h_piScr_data_wCuts.Fill(selMu_piScore)
  h_prScr_data_wCuts.Fill(selMu_prScore)
  h_elScr_data_wCuts.Fill(selMu_elScore)
  h_phScr_data_wCuts.Fill(selMu_phScore)
  h_pPScr_data_wCuts.Fill(selMu_primScore)
  h_pNScr_data_wCuts.Fill(selMu_ntrlScore)
  h_pCScr_data_wCuts.Fill(selMu_chgdScore)



print("finished event loops")

n_CCnumu_nCuts *= targetPOT/tnuPOTsum
n_NCnumu_nCuts *= targetPOT/tnuPOTsum
n_CCnue_nCuts *= targetPOT/tnuePOTsum
n_NCnue_nCuts *= targetPOT/tnuPOTsum
n_ext_nCuts *= targetPOT/textPOTsum
n_CCnumu_wCuts *= targetPOT/tnuPOTsum
n_NCnumu_wCuts *= targetPOT/tnuPOTsum
n_CCnue_wCuts *= targetPOT/tnuePOTsum
n_NCnue_wCuts *= targetPOT/tnuPOTsum
n_ext_wCuts *= targetPOT/textPOTsum
n_all_nCuts = n_CCnue_nCuts + n_NCnue_nCuts + n_CCnumu_nCuts + n_NCnumu_nCuts + n_ext_nCuts
n_all_wCuts = n_CCnue_wCuts + n_NCnue_wCuts + n_CCnumu_wCuts + n_NCnumu_wCuts + n_ext_wCuts

print()
print("n_CCnumu_nCuts:", n_CCnumu_nCuts)
print("n_ext_nCuts:", n_ext_nCuts)
print("n_NCnumu_nCuts:", n_NCnumu_nCuts)
print("n_CCnue_nCuts:", n_CCnue_nCuts)
print("n_NCnue_nCuts:", n_NCnue_nCuts)
print("n_all_nCuts:", n_all_nCuts)
print("n_data_nCuts:", n_data_nCuts)
print()
print("n_CCnumu_wCuts:", n_CCnumu_wCuts)
print("n_ext_wCuts:", n_ext_wCuts)
print("n_NCnumu_wCuts:", n_NCnumu_wCuts)
print("n_CCnue_wCuts:", n_CCnue_wCuts)
print("n_NCnue_wCuts:", n_NCnue_wCuts)
print("n_all_wCuts:", n_all_wCuts)
print("n_data_wCuts:", n_data_wCuts)
print()
print("CC numu cut efficiency: %f"%(n_CCnumu_wCuts/n_CCnumu_nCuts))
print("CC numu cut purity: %f"%(n_CCnumu_wCuts/n_all_wCuts))
print()


h_nuE_CCnumu_nCuts.Scale(targetPOT/tnuPOTsum)
h_nuE_CCnumu_wCuts.Scale(targetPOT/tnuPOTsum)
h_nuE_NCnumu_wCuts.Scale(targetPOT/tnuPOTsum)
h_nuE_CCnue_wCuts.Scale(targetPOT/tnuePOTsum)
h_nuE_NCnue_wCuts.Scale(targetPOT/tnuPOTsum)
h_nuE_CCnumu_eff.Divide(h_nuE_CCnumu_wCuts,h_nuE_CCnumu_nCuts,1,1,"B")
#h_nuE_CCnumu_eff.Divide(h_nuE_CCnumu_wCuts,h_nuE_CCnumu_nCuts,1,1)
h_nuE_all_wCuts.Add(h_nuE_CCnumu_wCuts)
h_nuE_all_wCuts.Add(h_nuE_NCnumu_wCuts)
h_nuE_all_wCuts.Add(h_nuE_CCnue_wCuts)
h_nuE_all_wCuts.Add(h_nuE_NCnue_wCuts)
#h_nuE_CCnumu_pur.Divide(h_nuE_CCnumu_wCuts,h_nuE_all_wCuts,1,1,"B")
h_nuE_CCnumu_pur.Divide(h_nuE_CCnumu_wCuts,h_nuE_all_wCuts)

h_nuE_CCnumu_wCuts1.Scale(targetPOT/tnuPOTsum)
h_nuE_CCnumu_wCuts2.Scale(targetPOT/tnuPOTsum)
h_nuE_CCnumu_wCutsAll.Scale(targetPOT/tnuPOTsum)
h_nuE_CCnumu_eff1.Divide(h_nuE_CCnumu_wCuts1,h_nuE_CCnumu_nCuts,1,1,"B")
h_nuE_CCnumu_eff2.Divide(h_nuE_CCnumu_wCuts2,h_nuE_CCnumu_nCuts,1,1,"B")
h_nuE_CCnumu_effAll.Divide(h_nuE_CCnumu_wCutsAll,h_nuE_CCnumu_nCuts,1,1,"B")

h_nuEr_CCnumu_nCuts.Scale(targetPOT/tnuPOTsum)
h_nuEr_CCnumu_wCuts.Scale(targetPOT/tnuPOTsum)
h_nuEr_NCnumu_wCuts.Scale(targetPOT/tnuPOTsum)
h_nuEr_CCnue_wCuts.Scale(targetPOT/tnuePOTsum)
h_nuEr_NCnue_wCuts.Scale(targetPOT/tnuPOTsum)
h_nuEr_ext_wCuts.Scale(targetPOT/textPOTsum)
h_nuEr_CCnumu_eff.Divide(h_nuEr_CCnumu_wCuts,h_nuEr_CCnumu_nCuts,1,1,"B")
#h_nuEr_CCnumu_eff.Divide(h_nuEr_CCnumu_wCuts,h_nuEr_CCnumu_nCuts,1,1)
h_nuEr_all_wCuts.Add(h_nuEr_CCnumu_wCuts)
h_nuEr_all_wCuts.Add(h_nuEr_NCnumu_wCuts)
h_nuEr_all_wCuts.Add(h_nuEr_CCnue_wCuts)
h_nuEr_all_wCuts.Add(h_nuEr_NCnue_wCuts)
h_nuEr_all_wCuts.Add(h_nuEr_ext_wCuts)
#h_nuEr_CCnumu_pur.Divide(h_nuEr_CCnumu_wCuts,h_nuEr_all_wCuts,1,1,"B")
h_nuEr_CCnumu_pur.Divide(h_nuEr_CCnumu_wCuts,h_nuEr_all_wCuts)

h_visE_CCnumu_wCuts.Scale(targetPOT/tnuPOTsum)
h_visE_NCnumu_wCuts.Scale(targetPOT/tnuPOTsum)
h_visE_CCnue_wCuts.Scale(targetPOT/tnuePOTsum)
h_visE_NCnue_wCuts.Scale(targetPOT/tnuPOTsum)
h_visE_ext_wCuts.Scale(targetPOT/textPOTsum)
h_visE_all_wCuts.Add(h_visE_ext_wCuts)
h_visE_all_wCuts.Add(h_visE_NCnue_wCuts)
h_visE_all_wCuts.Add(h_visE_CCnue_wCuts)
h_visE_all_wCuts.Add(h_visE_NCnumu_wCuts)
h_visE_all_wCuts.Add(h_visE_CCnumu_wCuts)
h_visE_predErr_wCuts.Add(h_visE_ext_wCuts)
h_visE_predErr_wCuts.Add(h_visE_NCnue_wCuts)
h_visE_predErr_wCuts.Add(h_visE_NCnumu_wCuts)
h_visE_predErr_wCuts.Add(h_visE_CCnumu_wCuts)
h_visE_predErr_wCuts.Add(h_visE_CCnue_wCuts)
h_visE_predErr_wCuts = SetUncertainties(h_visE_predErr_wCuts, "recoNuE", "CCnumu")

h_visE_CCnumu_wCutsSet1.Scale(targetPOT/tnuPOTsum)
h_visE_NCnumu_wCutsSet1.Scale(targetPOT/tnuPOTsum)
h_visE_CCnue_wCutsSet1.Scale(targetPOT/tnuePOTsum)
h_visE_NCnue_wCutsSet1.Scale(targetPOT/tnuPOTsum)
h_visE_ext_wCutsSet1.Scale(targetPOT/textPOTsum)
h_visE_all_wCutsSet1.Add(h_visE_ext_wCutsSet1)
h_visE_all_wCutsSet1.Add(h_visE_NCnue_wCutsSet1)
h_visE_all_wCutsSet1.Add(h_visE_CCnue_wCutsSet1)
h_visE_all_wCutsSet1.Add(h_visE_NCnumu_wCutsSet1)
h_visE_all_wCutsSet1.Add(h_visE_CCnumu_wCutsSet1)

h_visE_CCnumu_wCutsSet2.Scale(targetPOT/tnuPOTsum)
h_visE_NCnumu_wCutsSet2.Scale(targetPOT/tnuPOTsum)
h_visE_CCnue_wCutsSet2.Scale(targetPOT/tnuePOTsum)
h_visE_NCnue_wCutsSet2.Scale(targetPOT/tnuPOTsum)
h_visE_ext_wCutsSet2.Scale(targetPOT/textPOTsum)
h_visE_all_wCutsSet2.Add(h_visE_ext_wCutsSet2)
h_visE_all_wCutsSet2.Add(h_visE_NCnue_wCutsSet2)
h_visE_all_wCutsSet2.Add(h_visE_CCnue_wCutsSet2)
h_visE_all_wCutsSet2.Add(h_visE_NCnumu_wCutsSet2)
h_visE_all_wCutsSet2.Add(h_visE_CCnumu_wCutsSet2)

h_visE_CCnumu_wCutsSet3.Scale(targetPOT/tnuPOTsum)
h_visE_NCnumu_wCutsSet3.Scale(targetPOT/tnuPOTsum)
h_visE_CCnue_wCutsSet3.Scale(targetPOT/tnuePOTsum)
h_visE_NCnue_wCutsSet3.Scale(targetPOT/tnuPOTsum)
h_visE_ext_wCutsSet3.Scale(targetPOT/textPOTsum)
h_visE_all_wCutsSet3.Add(h_visE_ext_wCutsSet3)
h_visE_all_wCutsSet3.Add(h_visE_NCnue_wCutsSet3)
h_visE_all_wCutsSet3.Add(h_visE_CCnue_wCutsSet3)
h_visE_all_wCutsSet3.Add(h_visE_NCnumu_wCutsSet3)
h_visE_all_wCutsSet3.Add(h_visE_CCnumu_wCutsSet3)

h_visE_CCnumu_wCutsSet4.Scale(targetPOT/tnuPOTsum)
h_visE_NCnumu_wCutsSet4.Scale(targetPOT/tnuPOTsum)
h_visE_CCnue_wCutsSet4.Scale(targetPOT/tnuePOTsum)
h_visE_NCnue_wCutsSet4.Scale(targetPOT/tnuPOTsum)
h_visE_ext_wCutsSet4.Scale(targetPOT/textPOTsum)
h_visE_all_wCutsSet4.Add(h_visE_ext_wCutsSet4)
h_visE_all_wCutsSet4.Add(h_visE_NCnue_wCutsSet4)
h_visE_all_wCutsSet4.Add(h_visE_CCnue_wCutsSet4)
h_visE_all_wCutsSet4.Add(h_visE_NCnumu_wCutsSet4)
h_visE_all_wCutsSet4.Add(h_visE_CCnumu_wCutsSet4)

h_visE_CCnumu_wCutsSet5.Scale(targetPOT/tnuPOTsum)
h_visE_NCnumu_wCutsSet5.Scale(targetPOT/tnuPOTsum)
h_visE_CCnue_wCutsSet5.Scale(targetPOT/tnuePOTsum)
h_visE_NCnue_wCutsSet5.Scale(targetPOT/tnuPOTsum)
h_visE_ext_wCutsSet5.Scale(targetPOT/textPOTsum)
h_visE_all_wCutsSet5.Add(h_visE_ext_wCutsSet5)
h_visE_all_wCutsSet5.Add(h_visE_NCnue_wCutsSet5)
h_visE_all_wCutsSet5.Add(h_visE_CCnue_wCutsSet5)
h_visE_all_wCutsSet5.Add(h_visE_NCnumu_wCutsSet5)
h_visE_all_wCutsSet5.Add(h_visE_CCnumu_wCutsSet5)

h_visE_CCnumu_wCutsSet6.Scale(targetPOT/tnuPOTsum)
h_visE_NCnumu_wCutsSet6.Scale(targetPOT/tnuPOTsum)
h_visE_CCnue_wCutsSet6.Scale(targetPOT/tnuePOTsum)
h_visE_NCnue_wCutsSet6.Scale(targetPOT/tnuPOTsum)
h_visE_ext_wCutsSet6.Scale(targetPOT/textPOTsum)
h_visE_all_wCutsSet6.Add(h_visE_ext_wCutsSet6)
h_visE_all_wCutsSet6.Add(h_visE_NCnue_wCutsSet6)
h_visE_all_wCutsSet6.Add(h_visE_CCnue_wCutsSet6)
h_visE_all_wCutsSet6.Add(h_visE_NCnumu_wCutsSet6)
h_visE_all_wCutsSet6.Add(h_visE_CCnumu_wCutsSet6)

h_cosTheta_CCnumu_wCuts.Scale(targetPOT/tnuPOTsum)
h_cosTheta_NCnumu_wCuts.Scale(targetPOT/tnuPOTsum)
h_cosTheta_CCnue_wCuts.Scale(targetPOT/tnuePOTsum)
h_cosTheta_NCnue_wCuts.Scale(targetPOT/tnuPOTsum)
h_cosTheta_ext_wCuts.Scale(targetPOT/textPOTsum)
h_cosTheta_all_wCuts.Add(h_cosTheta_ext_wCuts)
h_cosTheta_all_wCuts.Add(h_cosTheta_NCnue_wCuts)
h_cosTheta_all_wCuts.Add(h_cosTheta_CCnue_wCuts)
h_cosTheta_all_wCuts.Add(h_cosTheta_NCnumu_wCuts)
h_cosTheta_all_wCuts.Add(h_cosTheta_CCnumu_wCuts)
h_cosTheta_predErr_wCuts.Add(h_cosTheta_ext_wCuts)
h_cosTheta_predErr_wCuts.Add(h_cosTheta_NCnue_wCuts)
h_cosTheta_predErr_wCuts.Add(h_cosTheta_CCnue_wCuts)
h_cosTheta_predErr_wCuts.Add(h_cosTheta_NCnumu_wCuts)
h_cosTheta_predErr_wCuts.Add(h_cosTheta_CCnumu_wCuts)
h_cosTheta_predErr_wCuts = SetUncertainties(h_cosTheta_predErr_wCuts, "cosTheta", "CCnumu")

h_lepP_CCnumu_wCuts.Scale(targetPOT/tnuPOTsum)
h_lepP_NCnumu_wCuts.Scale(targetPOT/tnuPOTsum)
h_lepP_CCnue_wCuts.Scale(targetPOT/tnuePOTsum)
h_lepP_NCnue_wCuts.Scale(targetPOT/tnuPOTsum)
h_lepP_ext_wCuts.Scale(targetPOT/textPOTsum)
h_lepP_all_wCuts.Add(h_lepP_ext_wCuts)
h_lepP_all_wCuts.Add(h_lepP_NCnue_wCuts)
h_lepP_all_wCuts.Add(h_lepP_CCnue_wCuts)
h_lepP_all_wCuts.Add(h_lepP_NCnumu_wCuts)
h_lepP_all_wCuts.Add(h_lepP_CCnumu_wCuts)
h_lepP_predErr_wCuts.Add(h_lepP_ext_wCuts)
h_lepP_predErr_wCuts.Add(h_lepP_NCnue_wCuts)
h_lepP_predErr_wCuts.Add(h_lepP_CCnue_wCuts)
h_lepP_predErr_wCuts.Add(h_lepP_NCnumu_wCuts)
h_lepP_predErr_wCuts.Add(h_lepP_CCnumu_wCuts)
h_lepP_predErr_wCuts = SetUncertainties(h_lepP_predErr_wCuts, "lepP", "CCnumu")

h_muScr_CCnumu_wCuts.Scale(targetPOT/tnuPOTsum)
h_muScr_NCnumu_wCuts.Scale(targetPOT/tnuPOTsum)
h_muScr_CCnue_wCuts.Scale(targetPOT/tnuePOTsum)
h_muScr_NCnue_wCuts.Scale(targetPOT/tnuPOTsum)
h_muScr_ext_wCuts.Scale(targetPOT/textPOTsum)
h_muScr_all_wCuts.Add(h_muScr_ext_wCuts)
h_muScr_all_wCuts.Add(h_muScr_NCnue_wCuts)
h_muScr_all_wCuts.Add(h_muScr_CCnue_wCuts)
h_muScr_all_wCuts.Add(h_muScr_NCnumu_wCuts)
h_muScr_all_wCuts.Add(h_muScr_CCnumu_wCuts)

h_piScr_CCnumu_wCuts.Scale(targetPOT/tnuPOTsum)
h_piScr_NCnumu_wCuts.Scale(targetPOT/tnuPOTsum)
h_piScr_CCnue_wCuts.Scale(targetPOT/tnuePOTsum)
h_piScr_NCnue_wCuts.Scale(targetPOT/tnuPOTsum)
h_piScr_ext_wCuts.Scale(targetPOT/textPOTsum)
h_piScr_all_wCuts.Add(h_piScr_ext_wCuts)
h_piScr_all_wCuts.Add(h_piScr_NCnue_wCuts)
h_piScr_all_wCuts.Add(h_piScr_CCnue_wCuts)
h_piScr_all_wCuts.Add(h_piScr_NCnumu_wCuts)
h_piScr_all_wCuts.Add(h_piScr_CCnumu_wCuts)

h_prScr_CCnumu_wCuts.Scale(targetPOT/tnuPOTsum)
h_prScr_NCnumu_wCuts.Scale(targetPOT/tnuPOTsum)
h_prScr_CCnue_wCuts.Scale(targetPOT/tnuePOTsum)
h_prScr_NCnue_wCuts.Scale(targetPOT/tnuPOTsum)
h_prScr_ext_wCuts.Scale(targetPOT/textPOTsum)
h_prScr_all_wCuts.Add(h_prScr_ext_wCuts)
h_prScr_all_wCuts.Add(h_prScr_NCnue_wCuts)
h_prScr_all_wCuts.Add(h_prScr_CCnue_wCuts)
h_prScr_all_wCuts.Add(h_prScr_NCnumu_wCuts)
h_prScr_all_wCuts.Add(h_prScr_CCnumu_wCuts)

h_elScr_CCnumu_wCuts.Scale(targetPOT/tnuPOTsum)
h_elScr_NCnumu_wCuts.Scale(targetPOT/tnuPOTsum)
h_elScr_CCnue_wCuts.Scale(targetPOT/tnuePOTsum)
h_elScr_NCnue_wCuts.Scale(targetPOT/tnuPOTsum)
h_elScr_ext_wCuts.Scale(targetPOT/textPOTsum)
h_elScr_all_wCuts.Add(h_elScr_ext_wCuts)
h_elScr_all_wCuts.Add(h_elScr_NCnue_wCuts)
h_elScr_all_wCuts.Add(h_elScr_CCnue_wCuts)
h_elScr_all_wCuts.Add(h_elScr_NCnumu_wCuts)
h_elScr_all_wCuts.Add(h_elScr_CCnumu_wCuts)

h_phScr_CCnumu_wCuts.Scale(targetPOT/tnuPOTsum)
h_phScr_NCnumu_wCuts.Scale(targetPOT/tnuPOTsum)
h_phScr_CCnue_wCuts.Scale(targetPOT/tnuePOTsum)
h_phScr_NCnue_wCuts.Scale(targetPOT/tnuPOTsum)
h_phScr_ext_wCuts.Scale(targetPOT/textPOTsum)
h_phScr_all_wCuts.Add(h_phScr_ext_wCuts)
h_phScr_all_wCuts.Add(h_phScr_NCnue_wCuts)
h_phScr_all_wCuts.Add(h_phScr_CCnue_wCuts)
h_phScr_all_wCuts.Add(h_phScr_NCnumu_wCuts)
h_phScr_all_wCuts.Add(h_phScr_CCnumu_wCuts)

h_pPScr_CCnumu_wCuts.Scale(targetPOT/tnuPOTsum)
h_pPScr_NCnumu_wCuts.Scale(targetPOT/tnuPOTsum)
h_pPScr_CCnue_wCuts.Scale(targetPOT/tnuePOTsum)
h_pPScr_NCnue_wCuts.Scale(targetPOT/tnuPOTsum)
h_pPScr_ext_wCuts.Scale(targetPOT/textPOTsum)
h_pPScr_all_wCuts.Add(h_pPScr_ext_wCuts)
h_pPScr_all_wCuts.Add(h_pPScr_NCnue_wCuts)
h_pPScr_all_wCuts.Add(h_pPScr_CCnue_wCuts)
h_pPScr_all_wCuts.Add(h_pPScr_NCnumu_wCuts)
h_pPScr_all_wCuts.Add(h_pPScr_CCnumu_wCuts)

h_pNScr_CCnumu_wCuts.Scale(targetPOT/tnuPOTsum)
h_pNScr_NCnumu_wCuts.Scale(targetPOT/tnuPOTsum)
h_pNScr_CCnue_wCuts.Scale(targetPOT/tnuePOTsum)
h_pNScr_NCnue_wCuts.Scale(targetPOT/tnuPOTsum)
h_pNScr_ext_wCuts.Scale(targetPOT/textPOTsum)
h_pNScr_all_wCuts.Add(h_pNScr_ext_wCuts)
h_pNScr_all_wCuts.Add(h_pNScr_NCnue_wCuts)
h_pNScr_all_wCuts.Add(h_pNScr_CCnue_wCuts)
h_pNScr_all_wCuts.Add(h_pNScr_NCnumu_wCuts)
h_pNScr_all_wCuts.Add(h_pNScr_CCnumu_wCuts)

h_pCScr_CCnumu_wCuts.Scale(targetPOT/tnuPOTsum)
h_pCScr_NCnumu_wCuts.Scale(targetPOT/tnuPOTsum)
h_pCScr_CCnue_wCuts.Scale(targetPOT/tnuePOTsum)
h_pCScr_NCnue_wCuts.Scale(targetPOT/tnuPOTsum)
h_pCScr_ext_wCuts.Scale(targetPOT/textPOTsum)
h_pCScr_all_wCuts.Add(h_pCScr_ext_wCuts)
h_pCScr_all_wCuts.Add(h_pCScr_NCnue_wCuts)
h_pCScr_all_wCuts.Add(h_pCScr_CCnue_wCuts)
h_pCScr_all_wCuts.Add(h_pCScr_NCnumu_wCuts)
h_pCScr_all_wCuts.Add(h_pCScr_CCnumu_wCuts)

h_lowENVsX_CCnumu_wCutsSet1.Scale(targetPOT/tnuPOTsum)
h_lowENVsX_NCnumu_wCutsSet1.Scale(targetPOT/tnuPOTsum)
h_lowENVsX_CCnue_wCutsSet1.Scale(targetPOT/tnuePOTsum)
h_lowENVsX_NCnue_wCutsSet1.Scale(targetPOT/tnuPOTsum)
h_lowENVsX_ext_wCutsSet1.Scale(targetPOT/textPOTsum)
h_lowEExVsX_wCutsSet1.Add(h_lowENVsX_CCnumu_wCutsSet1)
h_lowEExVsX_wCutsSet1.Add(h_lowENVsX_NCnumu_wCutsSet1)
h_lowEExVsX_wCutsSet1.Add(h_lowENVsX_CCnue_wCutsSet1)
h_lowEExVsX_wCutsSet1.Add(h_lowENVsX_NCnue_wCutsSet1)
h_lowEExVsX_wCutsSet1.Add(h_lowENVsX_ext_wCutsSet1)
h_lowEExVsX_wCutsSet1.Add(h_lowENVsX_data_wCutsSet1, -1)

h_lowENVsX_CCnumu_wCutsSet2.Scale(targetPOT/tnuPOTsum)
h_lowENVsX_NCnumu_wCutsSet2.Scale(targetPOT/tnuPOTsum)
h_lowENVsX_CCnue_wCutsSet2.Scale(targetPOT/tnuePOTsum)
h_lowENVsX_NCnue_wCutsSet2.Scale(targetPOT/tnuPOTsum)
h_lowENVsX_ext_wCutsSet2.Scale(targetPOT/textPOTsum)
h_lowEExVsX_wCutsSet2.Add(h_lowENVsX_CCnumu_wCutsSet2)
h_lowEExVsX_wCutsSet2.Add(h_lowENVsX_NCnumu_wCutsSet2)
h_lowEExVsX_wCutsSet2.Add(h_lowENVsX_CCnue_wCutsSet2)
h_lowEExVsX_wCutsSet2.Add(h_lowENVsX_NCnue_wCutsSet2)
h_lowEExVsX_wCutsSet2.Add(h_lowENVsX_ext_wCutsSet2)
h_lowEExVsX_wCutsSet2.Add(h_lowENVsX_data_wCutsSet2, -1)

h_lowENVsX_CCnumu_wCutsSet3.Scale(targetPOT/tnuPOTsum)
h_lowENVsX_NCnumu_wCutsSet3.Scale(targetPOT/tnuPOTsum)
h_lowENVsX_CCnue_wCutsSet3.Scale(targetPOT/tnuePOTsum)
h_lowENVsX_NCnue_wCutsSet3.Scale(targetPOT/tnuPOTsum)
h_lowENVsX_ext_wCutsSet3.Scale(targetPOT/textPOTsum)
h_lowEExVsX_wCutsSet3.Add(h_lowENVsX_CCnumu_wCutsSet3)
h_lowEExVsX_wCutsSet3.Add(h_lowENVsX_NCnumu_wCutsSet3)
h_lowEExVsX_wCutsSet3.Add(h_lowENVsX_CCnue_wCutsSet3)
h_lowEExVsX_wCutsSet3.Add(h_lowENVsX_NCnue_wCutsSet3)
h_lowEExVsX_wCutsSet3.Add(h_lowENVsX_ext_wCutsSet3)
h_lowEExVsX_wCutsSet3.Add(h_lowENVsX_data_wCutsSet3, -1)

h_lowENVsX_CCnumu_wCutsSet4.Scale(targetPOT/tnuPOTsum)
h_lowENVsX_NCnumu_wCutsSet4.Scale(targetPOT/tnuPOTsum)
h_lowENVsX_CCnue_wCutsSet4.Scale(targetPOT/tnuePOTsum)
h_lowENVsX_NCnue_wCutsSet4.Scale(targetPOT/tnuPOTsum)
h_lowENVsX_ext_wCutsSet4.Scale(targetPOT/textPOTsum)
h_lowEExVsX_wCutsSet4.Add(h_lowENVsX_CCnumu_wCutsSet4)
h_lowEExVsX_wCutsSet4.Add(h_lowENVsX_NCnumu_wCutsSet4)
h_lowEExVsX_wCutsSet4.Add(h_lowENVsX_CCnue_wCutsSet4)
h_lowEExVsX_wCutsSet4.Add(h_lowENVsX_NCnue_wCutsSet4)
h_lowEExVsX_wCutsSet4.Add(h_lowENVsX_ext_wCutsSet4)
h_lowEExVsX_wCutsSet4.Add(h_lowENVsX_data_wCutsSet4, -1)

h_lowENVsX_CCnumu_wCutsSet5.Scale(targetPOT/tnuPOTsum)
h_lowENVsX_NCnumu_wCutsSet5.Scale(targetPOT/tnuPOTsum)
h_lowENVsX_CCnue_wCutsSet5.Scale(targetPOT/tnuePOTsum)
h_lowENVsX_NCnue_wCutsSet5.Scale(targetPOT/tnuPOTsum)
h_lowENVsX_ext_wCutsSet5.Scale(targetPOT/textPOTsum)
h_lowEExVsX_wCutsSet5.Add(h_lowENVsX_CCnumu_wCutsSet5)
h_lowEExVsX_wCutsSet5.Add(h_lowENVsX_NCnumu_wCutsSet5)
h_lowEExVsX_wCutsSet5.Add(h_lowENVsX_CCnue_wCutsSet5)
h_lowEExVsX_wCutsSet5.Add(h_lowENVsX_NCnue_wCutsSet5)
h_lowEExVsX_wCutsSet5.Add(h_lowENVsX_ext_wCutsSet5)
h_lowEExVsX_wCutsSet5.Add(h_lowENVsX_data_wCutsSet5, -1)

h_lowENVsX_CCnumu_wCutsSet6.Scale(targetPOT/tnuPOTsum)
h_lowENVsX_NCnumu_wCutsSet6.Scale(targetPOT/tnuPOTsum)
h_lowENVsX_CCnue_wCutsSet6.Scale(targetPOT/tnuePOTsum)
h_lowENVsX_NCnue_wCutsSet6.Scale(targetPOT/tnuPOTsum)
h_lowENVsX_ext_wCutsSet6.Scale(targetPOT/textPOTsum)
h_lowEExVsX_wCutsSet6.Add(h_lowENVsX_CCnumu_wCutsSet6)
h_lowEExVsX_wCutsSet6.Add(h_lowENVsX_NCnumu_wCutsSet6)
h_lowEExVsX_wCutsSet6.Add(h_lowENVsX_CCnue_wCutsSet6)
h_lowEExVsX_wCutsSet6.Add(h_lowENVsX_NCnue_wCutsSet6)
h_lowEExVsX_wCutsSet6.Add(h_lowENVsX_ext_wCutsSet6)
h_lowEExVsX_wCutsSet6.Add(h_lowENVsX_data_wCutsSet6, -1)

h_lowENVsX_CCnumu_wCuts.Scale(targetPOT/tnuPOTsum)
h_lowENVsX_NCnumu_wCuts.Scale(targetPOT/tnuPOTsum)
h_lowENVsX_CCnue_wCuts.Scale(targetPOT/tnuePOTsum)
h_lowENVsX_NCnue_wCuts.Scale(targetPOT/tnuPOTsum)
h_lowENVsX_ext_wCuts.Scale(targetPOT/textPOTsum)
h_lowEExVsX_wCuts.Add(h_lowENVsX_CCnumu_wCuts)
h_lowEExVsX_wCuts.Add(h_lowENVsX_NCnumu_wCuts)
h_lowEExVsX_wCuts.Add(h_lowENVsX_CCnue_wCuts)
h_lowEExVsX_wCuts.Add(h_lowENVsX_NCnue_wCuts)
h_lowEExVsX_wCuts.Add(h_lowENVsX_ext_wCuts)
h_lowEExVsX_wCuts.Add(h_lowENVsX_data_wCuts, -1)

print("low energy (< 400 MeV) total excess with cut Set1:",h_lowEExVsX_wCutsSet1.Integral())
print("low energy (< 400 MeV) total excess with cut Set2:",h_lowEExVsX_wCutsSet2.Integral())
print("low energy (< 400 MeV) total excess with cut Set3:",h_lowEExVsX_wCutsSet3.Integral())
print("low energy (< 400 MeV) total excess with cut Set4:",h_lowEExVsX_wCutsSet4.Integral())
print("low energy (< 400 MeV) total excess with cut Set5:",h_lowEExVsX_wCutsSet5.Integral())
print("low energy (< 400 MeV) total excess with cut Set6:",h_lowEExVsX_wCutsSet6.Integral())
print("low energy (< 400 MeV) total excess with all cuts:",h_lowEExVsX_wCuts.Integral())
print()

h_cosFrac_CCnumu.Scale(targetPOT/tnuPOTsum)
h_cosFrac_NCnumu.Scale(targetPOT/tnuPOTsum)
h_cosFrac_CCnue.Scale(targetPOT/tnuePOTsum)
h_cosFrac_NCnue.Scale(targetPOT/tnuPOTsum)
h_cosFrac_ext.Scale(targetPOT/textPOTsum)

h_nMu_CCnumu.Scale(targetPOT/tnuPOTsum)
h_nMu_NCnumu.Scale(targetPOT/tnuPOTsum)
h_nMu_CCnue.Scale(targetPOT/tnuePOTsum)
h_nMu_NCnue.Scale(targetPOT/tnuPOTsum)
h_nMu_ext.Scale(targetPOT/textPOTsum)

h_cosKPDist_CCnumu.Scale(targetPOT/tnuPOTsum)
h_cosKPDist_NCnumu.Scale(targetPOT/tnuPOTsum)
h_cosKPDist_CCnue.Scale(targetPOT/tnuePOTsum)
h_cosKPDist_NCnue.Scale(targetPOT/tnuPOTsum)
h_cosKPDist_ext.Scale(targetPOT/textPOTsum)

h_vtxScore_CCnumu.Scale(targetPOT/tnuPOTsum)
h_vtxScore_NCnumu.Scale(targetPOT/tnuPOTsum)
h_vtxScore_CCnue.Scale(targetPOT/tnuePOTsum)
h_vtxScore_NCnue.Scale(targetPOT/tnuPOTsum)
h_vtxScore_ext.Scale(targetPOT/textPOTsum)

h_pcaEVRatio_CCnumu.Scale(targetPOT/tnuPOTsum)
h_pcaEVRatio_NCnumu.Scale(targetPOT/tnuPOTsum)
h_pcaEVRatio_CCnue.Scale(targetPOT/tnuePOTsum)
h_pcaEVRatio_NCnue.Scale(targetPOT/tnuPOTsum)
h_pcaEVRatio_ext.Scale(targetPOT/textPOTsum)

h_pcaCosThetaBeam_CCnumu.Scale(targetPOT/tnuPOTsum)
h_pcaCosThetaBeam_NCnumu.Scale(targetPOT/tnuPOTsum)
h_pcaCosThetaBeam_CCnue.Scale(targetPOT/tnuePOTsum)
h_pcaCosThetaBeam_NCnue.Scale(targetPOT/tnuPOTsum)
h_pcaCosThetaBeam_ext.Scale(targetPOT/textPOTsum)

h_pcaCosThetaGrav_CCnumu.Scale(targetPOT/tnuPOTsum)
h_pcaCosThetaGrav_NCnumu.Scale(targetPOT/tnuPOTsum)
h_pcaCosThetaGrav_CCnue.Scale(targetPOT/tnuePOTsum)
h_pcaCosThetaGrav_NCnue.Scale(targetPOT/tnuPOTsum)
h_pcaCosThetaGrav_ext.Scale(targetPOT/textPOTsum)

h_selMu_cosThetaBeam_CCnumu.Scale(targetPOT/tnuPOTsum)
h_selMu_cosThetaBeam_NCnumu.Scale(targetPOT/tnuPOTsum)
h_selMu_cosThetaBeam_CCnue.Scale(targetPOT/tnuePOTsum)
h_selMu_cosThetaBeam_NCnue.Scale(targetPOT/tnuPOTsum)
h_selMu_cosThetaBeam_ext.Scale(targetPOT/textPOTsum)

h_selMu_cosThetaGrav_CCnumu.Scale(targetPOT/tnuPOTsum)
h_selMu_cosThetaGrav_NCnumu.Scale(targetPOT/tnuPOTsum)
h_selMu_cosThetaGrav_CCnue.Scale(targetPOT/tnuePOTsum)
h_selMu_cosThetaGrav_NCnue.Scale(targetPOT/tnuPOTsum)
h_selMu_cosThetaGrav_ext.Scale(targetPOT/textPOTsum)

h_selMu_muScore_CCnumu.Scale(targetPOT/tnuPOTsum)
h_selMu_muScore_NCnumu.Scale(targetPOT/tnuPOTsum)
h_selMu_muScore_CCnue.Scale(targetPOT/tnuePOTsum)
h_selMu_muScore_NCnue.Scale(targetPOT/tnuPOTsum)
h_selMu_muScore_ext.Scale(targetPOT/textPOTsum)

h_selMu_piScore_CCnumu.Scale(targetPOT/tnuPOTsum)
h_selMu_piScore_NCnumu.Scale(targetPOT/tnuPOTsum)
h_selMu_piScore_CCnue.Scale(targetPOT/tnuePOTsum)
h_selMu_piScore_NCnue.Scale(targetPOT/tnuPOTsum)
h_selMu_piScore_ext.Scale(targetPOT/textPOTsum)

h_selMu_scoreConf_CCnumu.Scale(targetPOT/tnuPOTsum)
h_selMu_scoreConf_NCnumu.Scale(targetPOT/tnuPOTsum)
h_selMu_scoreConf_CCnue.Scale(targetPOT/tnuePOTsum)
h_selMu_scoreConf_NCnue.Scale(targetPOT/tnuPOTsum)
h_selMu_scoreConf_ext.Scale(targetPOT/textPOTsum)

h_selMu_proc_CCnumu.Scale(targetPOT/tnuPOTsum)
h_selMu_proc_NCnumu.Scale(targetPOT/tnuPOTsum)
h_selMu_proc_CCnue.Scale(targetPOT/tnuePOTsum)
h_selMu_proc_NCnue.Scale(targetPOT/tnuPOTsum)
h_selMu_proc_ext.Scale(targetPOT/textPOTsum)

h_selMu_primScore_CCnumu.Scale(targetPOT/tnuPOTsum)
h_selMu_primScore_NCnumu.Scale(targetPOT/tnuPOTsum)
h_selMu_primScore_CCnue.Scale(targetPOT/tnuePOTsum)
h_selMu_primScore_NCnue.Scale(targetPOT/tnuPOTsum)
h_selMu_primScore_ext.Scale(targetPOT/textPOTsum)

h_selMu_ntrlScore_CCnumu.Scale(targetPOT/tnuPOTsum)
h_selMu_ntrlScore_NCnumu.Scale(targetPOT/tnuPOTsum)
h_selMu_ntrlScore_CCnue.Scale(targetPOT/tnuePOTsum)
h_selMu_ntrlScore_NCnue.Scale(targetPOT/tnuPOTsum)
h_selMu_ntrlScore_ext.Scale(targetPOT/textPOTsum)

h_selMu_chgdScore_CCnumu.Scale(targetPOT/tnuPOTsum)
h_selMu_chgdScore_NCnumu.Scale(targetPOT/tnuPOTsum)
h_selMu_chgdScore_CCnue.Scale(targetPOT/tnuePOTsum)
h_selMu_chgdScore_NCnue.Scale(targetPOT/tnuPOTsum)
h_selMu_chgdScore_ext.Scale(targetPOT/textPOTsum)

h_selMu_procConf_CCnumu.Scale(targetPOT/tnuPOTsum)
h_selMu_procConf_NCnumu.Scale(targetPOT/tnuPOTsum)
h_selMu_procConf_CCnue.Scale(targetPOT/tnuePOTsum)
h_selMu_procConf_NCnue.Scale(targetPOT/tnuPOTsum)
h_selMu_procConf_ext.Scale(targetPOT/textPOTsum)



def configureLegend(leg, h_CCnumu, h_NCnumu, h_CCnue, h_NCnue, h_ext):
  leg.AddEntry(h_CCnue,"CC nue (%.2f)"%h_CCnue.Integral(),"l")
  leg.AddEntry(h_CCnumu,"CC numu (%.2f)"%h_CCnumu.Integral(),"l")
  leg.AddEntry(h_NCnue,"NC nue (%.2f)"%h_NCnue.Integral(),"l")
  leg.AddEntry(h_NCnumu,"NC numu (%.2f)"%h_NCnumu.Integral(),"l")
  leg.AddEntry(h_ext,"cosmic background (%.2f)"%h_ext.Integral(),"l")
  return leg
  

outFile = rt.TFile(args.outfile, "RECREATE")

cnv_cosFrac = rt.TCanvas("cnv_cosFrac","cnv_cosFrac")
hists_cosFrac = sortHists([h_cosFrac_CCnumu, h_cosFrac_NCnumu, h_cosFrac_CCnue, h_cosFrac_NCnue, h_cosFrac_ext])
hists_cosFrac[0].Draw("EHIST")
for i in range(1,len(hists_cosFrac)):
  hists_cosFrac[i].Draw("EHISTSAME")
leg_cosFrac = rt.TLegend(0.7,0.7,0.9,0.9)
leg_cosFrac = configureLegend(leg_cosFrac, h_cosFrac_CCnumu,
  h_cosFrac_NCnumu, h_cosFrac_CCnue, h_cosFrac_NCnue, h_cosFrac_ext)
leg_cosFrac.Draw()
#cnv_cosFrac.SaveAs("cosFrac.png")
cnv_cosFrac.Write()

cnv_nMu = rt.TCanvas("cnv_nMu","cnv_nMu")
hists_nMu = sortHists([h_nMu_CCnumu, h_nMu_NCnumu, h_nMu_CCnue, h_nMu_NCnue, h_nMu_ext])
hists_nMu[0].Draw("EHIST")
for i in range(1,len(hists_nMu)):
  hists_nMu[i].Draw("EHISTSAME")
leg_nMu = rt.TLegend(0.7,0.7,0.9,0.9)
leg_nMu = configureLegend(leg_nMu, h_nMu_CCnumu,
  h_nMu_NCnumu, h_nMu_CCnue, h_nMu_NCnue, h_nMu_ext)
leg_nMu.Draw()
#cnv_nMu.SaveAs("nMu.png")
cnv_nMu.Write()

cnv_vtxScore = rt.TCanvas("cnv_vtxScore","cnv_vtxScore")
hists_vtxScore = sortHists([h_vtxScore_CCnumu, h_vtxScore_NCnumu, h_vtxScore_CCnue, h_vtxScore_NCnue, h_vtxScore_ext])
hists_vtxScore[0].Draw("EHIST")
for i in range(1,len(hists_vtxScore)):
  hists_vtxScore[i].Draw("EHISTSAME")
leg_vtxScore = rt.TLegend(0.7,0.7,0.9,0.9)
leg_vtxScore = configureLegend(leg_vtxScore, h_vtxScore_CCnumu,
  h_vtxScore_NCnumu, h_vtxScore_CCnue, h_vtxScore_NCnue, h_vtxScore_ext)
leg_vtxScore.Draw()
#cnv_vtxScore.SaveAs("vtxScore.png")
cnv_vtxScore.Write()

cnv_pcaEVRatio = rt.TCanvas("cnv_pcaEVRatio","cnv_pcaEVRatio")
hists_pcaEVRatio = sortHists([h_pcaEVRatio_CCnumu, h_pcaEVRatio_NCnumu, h_pcaEVRatio_CCnue, h_pcaEVRatio_NCnue, h_pcaEVRatio_ext])
hists_pcaEVRatio[0].Draw("EHIST")
for i in range(1,len(hists_pcaEVRatio)):
  hists_pcaEVRatio[i].Draw("EHISTSAME")
leg_pcaEVRatio = rt.TLegend(0.7,0.7,0.9,0.9)
leg_pcaEVRatio = configureLegend(leg_pcaEVRatio, h_pcaEVRatio_CCnumu,
  h_pcaEVRatio_NCnumu, h_pcaEVRatio_CCnue, h_pcaEVRatio_NCnue, h_pcaEVRatio_ext)
leg_pcaEVRatio.Draw()
#cnv_pcaEVRatio.SaveAs("pcaEVRatio.png")
cnv_pcaEVRatio.Write()

cnv_pcaCosThetaBeam = rt.TCanvas("cnv_pcaCosThetaBeam","cnv_pcaCosThetaBeam")
hists_pcaCosThetaBeam = sortHists([h_pcaCosThetaBeam_CCnumu, h_pcaCosThetaBeam_NCnumu, h_pcaCosThetaBeam_CCnue, h_pcaCosThetaBeam_NCnue, h_pcaCosThetaBeam_ext])
hists_pcaCosThetaBeam[0].Draw("EHIST")
for i in range(1,len(hists_pcaCosThetaBeam)):
  hists_pcaCosThetaBeam[i].Draw("EHISTSAME")
leg_pcaCosThetaBeam = rt.TLegend(0.7,0.7,0.9,0.9)
leg_pcaCosThetaBeam = configureLegend(leg_pcaCosThetaBeam, h_pcaCosThetaBeam_CCnumu,
  h_pcaCosThetaBeam_NCnumu, h_pcaCosThetaBeam_CCnue, h_pcaCosThetaBeam_NCnue, h_pcaCosThetaBeam_ext)
leg_pcaCosThetaBeam.Draw()
#cnv_pcaCosThetaBeam.SaveAs("pcaCosThetaBeam.png")
cnv_pcaCosThetaBeam.Write()

cnv_pcaCosThetaGrav = rt.TCanvas("cnv_pcaCosThetaGrav","cnv_pcaCosThetaGrav")
hists_pcaCosThetaGrav = sortHists([h_pcaCosThetaGrav_CCnumu, h_pcaCosThetaGrav_NCnumu, h_pcaCosThetaGrav_CCnue, h_pcaCosThetaGrav_NCnue, h_pcaCosThetaGrav_ext])
hists_pcaCosThetaGrav[0].Draw("EHIST")
for i in range(1,len(hists_pcaCosThetaGrav)):
  hists_pcaCosThetaGrav[i].Draw("EHISTSAME")
leg_pcaCosThetaGrav = rt.TLegend(0.7,0.7,0.9,0.9)
leg_pcaCosThetaGrav = configureLegend(leg_pcaCosThetaGrav, h_pcaCosThetaGrav_CCnumu,
  h_pcaCosThetaGrav_NCnumu, h_pcaCosThetaGrav_CCnue, h_pcaCosThetaGrav_NCnue, h_pcaCosThetaGrav_ext)
leg_pcaCosThetaGrav.Draw()
#cnv_pcaCosThetaGrav.SaveAs("pcaCosThetaGrav.png")
cnv_pcaCosThetaGrav.Write()

cnv_selMu_cosThetaBeam = rt.TCanvas("cnv_selMu_cosThetaBeam","cnv_selMu_cosThetaBeam")
hists_selMu_cosThetaBeam = sortHists([h_selMu_cosThetaBeam_CCnumu, h_selMu_cosThetaBeam_NCnumu, h_selMu_cosThetaBeam_CCnue, h_selMu_cosThetaBeam_NCnue, h_selMu_cosThetaBeam_ext])
hists_selMu_cosThetaBeam[0].Draw("EHIST")
for i in range(1,len(hists_selMu_cosThetaBeam)):
  hists_selMu_cosThetaBeam[i].Draw("EHISTSAME")
leg_selMu_cosThetaBeam = rt.TLegend(0.7,0.7,0.9,0.9)
leg_selMu_cosThetaBeam = configureLegend(leg_selMu_cosThetaBeam, h_selMu_cosThetaBeam_CCnumu,
  h_selMu_cosThetaBeam_NCnumu, h_selMu_cosThetaBeam_CCnue, h_selMu_cosThetaBeam_NCnue, h_selMu_cosThetaBeam_ext)
leg_selMu_cosThetaBeam.Draw()
#cnv_selMu_cosThetaBeam.SaveAs("selMu_cosThetaBeam.png")
cnv_selMu_cosThetaBeam.Write()

cnv_selMu_cosThetaGrav = rt.TCanvas("cnv_selMu_cosThetaGrav","cnv_selMu_cosThetaGrav")
hists_selMu_cosThetaGrav = sortHists([h_selMu_cosThetaGrav_CCnumu, h_selMu_cosThetaGrav_NCnumu, h_selMu_cosThetaGrav_CCnue, h_selMu_cosThetaGrav_NCnue, h_selMu_cosThetaGrav_ext])
hists_selMu_cosThetaGrav[0].Draw("EHIST")
for i in range(1,len(hists_selMu_cosThetaGrav)):
  hists_selMu_cosThetaGrav[i].Draw("EHISTSAME")
leg_selMu_cosThetaGrav = rt.TLegend(0.7,0.7,0.9,0.9)
leg_selMu_cosThetaGrav = configureLegend(leg_selMu_cosThetaGrav, h_selMu_cosThetaGrav_CCnumu,
  h_selMu_cosThetaGrav_NCnumu, h_selMu_cosThetaGrav_CCnue, h_selMu_cosThetaGrav_NCnue, h_selMu_cosThetaGrav_ext)
leg_selMu_cosThetaGrav.Draw()
#cnv_selMu_cosThetaGrav.SaveAs("selMu_cosThetaGrav.png")
cnv_selMu_cosThetaGrav.Write()

cnv_cosKPDist = rt.TCanvas("cnv_cosKPDist","cnv_cosKPDist")
hists_cosKPDist = sortHists([h_cosKPDist_CCnumu, h_cosKPDist_NCnumu, h_cosKPDist_CCnue, h_cosKPDist_NCnue, h_cosKPDist_ext])
hists_cosKPDist[0].Draw("EHIST")
for i in range(1,len(hists_cosKPDist)):
  hists_cosKPDist[i].Draw("EHISTSAME")
leg_cosKPDist = rt.TLegend(0.7,0.7,0.9,0.9)
leg_cosKPDist = configureLegend(leg_cosKPDist, h_cosKPDist_CCnumu,
  h_cosKPDist_NCnumu, h_cosKPDist_CCnue, h_cosKPDist_NCnue, h_cosKPDist_ext)
leg_cosKPDist.Draw()
#cnv_cosKPDist.SaveAs("cosKPDist.png")
cnv_cosKPDist.Write()

cnv_selMu_muScore = rt.TCanvas("cnv_selMu_muScore","cnv_selMu_muScore")
hists_selMu_muScore = sortHists([h_selMu_muScore_CCnumu, h_selMu_muScore_NCnumu, h_selMu_muScore_CCnue, h_selMu_muScore_NCnue, h_selMu_muScore_ext])
hists_selMu_muScore[0].Draw("EHIST")
for i in range(1,len(hists_selMu_muScore)):
  hists_selMu_muScore[i].Draw("EHISTSAME")
leg_selMu_muScore = rt.TLegend(0.7,0.7,0.9,0.9)
leg_selMu_muScore = configureLegend(leg_selMu_muScore, h_selMu_muScore_CCnumu,
  h_selMu_muScore_NCnumu, h_selMu_muScore_CCnue, h_selMu_muScore_NCnue, h_selMu_muScore_ext)
leg_selMu_muScore.Draw()
#cnv_selMu_muScore.SaveAs("selMu_muScore.png")
cnv_selMu_muScore.Write()

cnv_selMu_piScore = rt.TCanvas("cnv_selMu_piScore","cnv_selMu_piScore")
hists_selMu_piScore = sortHists([h_selMu_piScore_CCnumu, h_selMu_piScore_NCnumu, h_selMu_piScore_CCnue, h_selMu_piScore_NCnue, h_selMu_piScore_ext])
hists_selMu_piScore[0].Draw("EHIST")
for i in range(1,len(hists_selMu_piScore)):
  hists_selMu_piScore[i].Draw("EHISTSAME")
leg_selMu_piScore = rt.TLegend(0.7,0.7,0.9,0.9)
leg_selMu_piScore = configureLegend(leg_selMu_piScore, h_selMu_piScore_CCnumu,
  h_selMu_piScore_NCnumu, h_selMu_piScore_CCnue, h_selMu_piScore_NCnue, h_selMu_piScore_ext)
leg_selMu_piScore.Draw()
#cnv_selMu_piScore.SaveAs("selMu_piScore.png")
cnv_selMu_piScore.Write()

cnv_selMu_scoreConf = rt.TCanvas("cnv_selMu_scoreConf","cnv_selMu_scoreConf")
hists_selMu_scoreConf = sortHists([h_selMu_scoreConf_CCnumu, h_selMu_scoreConf_NCnumu, h_selMu_scoreConf_CCnue, h_selMu_scoreConf_NCnue, h_selMu_scoreConf_ext])
hists_selMu_scoreConf[0].Draw("EHIST")
for i in range(1,len(hists_selMu_scoreConf)):
  hists_selMu_scoreConf[i].Draw("EHISTSAME")
leg_selMu_scoreConf = rt.TLegend(0.7,0.7,0.9,0.9)
leg_selMu_scoreConf = configureLegend(leg_selMu_scoreConf, h_selMu_scoreConf_CCnumu,
  h_selMu_scoreConf_NCnumu, h_selMu_scoreConf_CCnue, h_selMu_scoreConf_NCnue, h_selMu_scoreConf_ext)
leg_selMu_scoreConf.Draw()
#cnv_selMu_scoreConf.SaveAs("selMu_scoreConf.png")
cnv_selMu_scoreConf.Write()

cnv_selMu_proc = rt.TCanvas("cnv_selMu_proc","cnv_selMu_proc")
hists_selMu_proc = sortHists([h_selMu_proc_CCnumu, h_selMu_proc_NCnumu, h_selMu_proc_CCnue, h_selMu_proc_NCnue, h_selMu_proc_ext])
hists_selMu_proc[0].Draw("EHIST")
for i in range(1,len(hists_selMu_proc)):
  hists_selMu_proc[i].Draw("EHISTSAME")
leg_selMu_proc = rt.TLegend(0.7,0.7,0.9,0.9)
leg_selMu_proc = configureLegend(leg_selMu_proc, h_selMu_proc_CCnumu,
  h_selMu_proc_NCnumu, h_selMu_proc_CCnue, h_selMu_proc_NCnue, h_selMu_proc_ext)
leg_selMu_proc.Draw()
#cnv_selMu_proc.SaveAs("selMu_proc.png")
cnv_selMu_proc.Write()

cnv_selMu_primScore = rt.TCanvas("cnv_selMu_primScore","cnv_selMu_primScore")
hists_selMu_primScore = sortHists([h_selMu_primScore_CCnumu, h_selMu_primScore_NCnumu, h_selMu_primScore_CCnue, h_selMu_primScore_NCnue, h_selMu_primScore_ext])
hists_selMu_primScore[0].Draw("EHIST")
for i in range(1,len(hists_selMu_primScore)):
  hists_selMu_primScore[i].Draw("EHISTSAME")
leg_selMu_primScore = rt.TLegend(0.7,0.7,0.9,0.9)
leg_selMu_primScore = configureLegend(leg_selMu_primScore, h_selMu_primScore_CCnumu,
  h_selMu_primScore_NCnumu, h_selMu_primScore_CCnue, h_selMu_primScore_NCnue, h_selMu_primScore_ext)
leg_selMu_primScore.Draw()
#cnv_selMu_primScore.SaveAs("selMu_primScore.png")
cnv_selMu_primScore.Write()

cnv_selMu_ntrlScore = rt.TCanvas("cnv_selMu_ntrlScore","cnv_selMu_ntrlScore")
hists_selMu_ntrlScore = sortHists([h_selMu_ntrlScore_CCnumu, h_selMu_ntrlScore_NCnumu, h_selMu_ntrlScore_CCnue, h_selMu_ntrlScore_NCnue, h_selMu_ntrlScore_ext])
hists_selMu_ntrlScore[0].Draw("EHIST")
for i in range(1,len(hists_selMu_ntrlScore)):
  hists_selMu_ntrlScore[i].Draw("EHISTSAME")
leg_selMu_ntrlScore = rt.TLegend(0.7,0.7,0.9,0.9)
leg_selMu_ntrlScore = configureLegend(leg_selMu_ntrlScore, h_selMu_ntrlScore_CCnumu,
  h_selMu_ntrlScore_NCnumu, h_selMu_ntrlScore_CCnue, h_selMu_ntrlScore_NCnue, h_selMu_ntrlScore_ext)
leg_selMu_ntrlScore.Draw()
#cnv_selMu_ntrlScore.SaveAs("selMu_ntrlScore.png")
cnv_selMu_ntrlScore.Write()

cnv_selMu_chgdScore = rt.TCanvas("cnv_selMu_chgdScore","cnv_selMu_chgdScore")
hists_selMu_chgdScore = sortHists([h_selMu_chgdScore_CCnumu, h_selMu_chgdScore_NCnumu, h_selMu_chgdScore_CCnue, h_selMu_chgdScore_NCnue, h_selMu_chgdScore_ext])
hists_selMu_chgdScore[0].Draw("EHIST")
for i in range(1,len(hists_selMu_chgdScore)):
  hists_selMu_chgdScore[i].Draw("EHISTSAME")
leg_selMu_chgdScore = rt.TLegend(0.7,0.7,0.9,0.9)
leg_selMu_chgdScore = configureLegend(leg_selMu_chgdScore, h_selMu_chgdScore_CCnumu,
  h_selMu_chgdScore_NCnumu, h_selMu_chgdScore_CCnue, h_selMu_chgdScore_NCnue, h_selMu_chgdScore_ext)
leg_selMu_chgdScore.Draw()
#cnv_selMu_chgdScore.SaveAs("selMu_chgdScore.png")
cnv_selMu_chgdScore.Write()

cnv_selMu_procConf = rt.TCanvas("cnv_selMu_procConf","cnv_selMu_procConf")
hists_selMu_procConf = sortHists([h_selMu_procConf_CCnumu, h_selMu_procConf_NCnumu, h_selMu_procConf_CCnue, h_selMu_procConf_NCnue, h_selMu_procConf_ext])
hists_selMu_procConf[0].Draw("EHIST")
for i in range(1,len(hists_selMu_procConf)):
  hists_selMu_procConf[i].Draw("EHISTSAME")
leg_selMu_procConf = rt.TLegend(0.7,0.7,0.9,0.9)
leg_selMu_procConf = configureLegend(leg_selMu_procConf, h_selMu_procConf_CCnumu,
  h_selMu_procConf_NCnumu, h_selMu_procConf_CCnue, h_selMu_procConf_NCnue, h_selMu_procConf_ext)
leg_selMu_procConf.Draw()
#cnv_selMu_procConf.SaveAs("selMu_procConf.png")
cnv_selMu_procConf.Write()



h_nuE_CCnumu_eff.GetYaxis().SetRangeUser(0,1.003)
h_nuE_CCnumu_pur.GetYaxis().SetRangeUser(0,1.003)
#h_nuE_CCnumu_pur.SetLineColor(rt.kRed)
if not args.plotPurityVsTrueE:
  h_nuE_CCnumu_eff.GetYaxis().SetTitle("efficiency")
cnv_CCnumu_sel_trueE = rt.TCanvas("cnv_CCnumu_sel_trueE","cnv_CCnumu_sel_trueE")
cnv_CCnumu_sel_trueE.SetGrid()
h_nuE_CCnumu_eff.Draw("E")
if args.plotPurityVsTrueE:
  h_nuE_CCnumu_pur.Draw("ESAME")
  leg_CCnumu_sel_trueE = rt.TLegend(0.7,0.7,0.9,0.9)
  leg_CCnumu_sel_trueE.AddEntry(h_nuE_CCnumu_eff, "efficiency", "l")
  leg_CCnumu_sel_trueE.AddEntry(h_nuE_CCnumu_pur, "purity", "l")
  leg_CCnumu_sel_trueE.Draw()
cnv_CCnumu_sel_trueE.Write()

h_nuEr_CCnumu_eff.GetYaxis().SetRangeUser(0,1.003)
h_nuEr_CCnumu_pur.GetYaxis().SetRangeUser(0,1.003)
cnv_CCnumu_sel_recoE = rt.TCanvas("cnv_CCnumu_sel_recoE","cnv_CCnumu_sel_recoE")
cnv_CCnumu_sel_recoE.SetGrid()
h_nuEr_CCnumu_eff.Draw("E")
h_nuEr_CCnumu_pur.Draw("ESAME")
leg_CCnumu_sel_recoE = rt.TLegend(0.7,0.7,0.9,0.9)
leg_CCnumu_sel_recoE.AddEntry(h_nuEr_CCnumu_eff, "efficiency", "l")
leg_CCnumu_sel_recoE.AddEntry(h_nuEr_CCnumu_pur, "purity", "l")
leg_CCnumu_sel_recoE.Draw()
cnv_CCnumu_sel_recoE.Write()

cnv_visE_sel = rt.TCanvas("cnv_visE_sel", "cnv_visE_sel")
h_visE_data_wCuts.Draw("E")
h_visE_all_wCuts.Draw("hist same")
h_visE_data_wCuts.Draw("ESAME")
h_visE_predErr_wCuts.SetFillColor(4)
h_visE_predErr_wCuts.SetFillStyle(3002)
h_visE_predErr_wCuts.Draw("E2SAME")
if args.recoEOverflow:
  label_visE_sel = getOverflowLabel(h_visE_data_wCuts)
  label_visE_sel.Draw()
leg_visE_sel = rt.TLegend(0.7,0.7,0.9,0.9)
leg_visE_sel.AddEntry(h_visE_ext_wCuts, "cosmic background (%.2f)"%h_visE_ext_wCuts.Integral(), "f")
leg_visE_sel.AddEntry(h_visE_NCnue_wCuts, "NC nue (%.2f)"%h_visE_NCnue_wCuts.Integral(), "f")
leg_visE_sel.AddEntry(h_visE_CCnue_wCuts, "CC nue (%.2f)"%h_visE_CCnue_wCuts.Integral(), "f")
leg_visE_sel.AddEntry(h_visE_NCnumu_wCuts, "NC numu (%.2f)"%h_visE_NCnumu_wCuts.Integral(), "f")
leg_visE_sel.AddEntry(h_visE_CCnumu_wCuts, "CC numu (%.2f)"%h_visE_CCnumu_wCuts.Integral(), "f")
leg_visE_sel.AddEntry(h_visE_data_wCuts, "5e19 data (%.2f)"%h_visE_data_wCuts.Integral(), "l")
leg_visE_sel.AddEntry(h_visE_predErr_wCuts, "Pred. stats + sys (no detvar) unc.", "f")
leg_visE_sel.Draw()
cnv_visE_sel.Write()
h_visE_data_wCuts.Write()
h_visE_all_wCuts.Write()
h_visE_predErr_wCuts.Write()

cnv_visE_cutSet1 = rt.TCanvas("cnv_visE_cutSet1", "cnv_visE_cutSet1")
h_visE_data_wCutsSet1.Draw("E")
h_visE_all_wCutsSet1.Draw("hist same")
h_visE_data_wCutsSet1.Draw("ESAME")
if args.recoEOverflow:
  label_visE_cutSet1 = getOverflowLabel(h_visE_data_wCutsSet1)
  label_visE_cutSet1.Draw()
leg_visE_cutSet1 = rt.TLegend(0.7,0.7,0.9,0.9)
leg_visE_cutSet1.AddEntry(h_visE_ext_wCutsSet1, "cosmic background (%.2f)"%h_visE_ext_wCutsSet1.Integral(), "f")
leg_visE_cutSet1.AddEntry(h_visE_NCnue_wCutsSet1, "NC nue (%.2f)"%h_visE_NCnue_wCutsSet1.Integral(), "f")
leg_visE_cutSet1.AddEntry(h_visE_CCnue_wCutsSet1, "CC nue (%.2f)"%h_visE_CCnue_wCutsSet1.Integral(), "f")
leg_visE_cutSet1.AddEntry(h_visE_NCnumu_wCutsSet1, "NC numu (%.2f)"%h_visE_NCnumu_wCutsSet1.Integral(), "f")
leg_visE_cutSet1.AddEntry(h_visE_CCnumu_wCutsSet1, "CC numu (%.2f)"%h_visE_CCnumu_wCutsSet1.Integral(), "f")
leg_visE_cutSet1.AddEntry(h_visE_data_wCutsSet1, "5e19 data (%.2f)"%h_visE_data_wCutsSet1.Integral(), "l")
leg_visE_cutSet1.Draw()
cnv_visE_cutSet1.Write()

cnv_visE_cutSet2 = rt.TCanvas("cnv_visE_cutSet2", "cnv_visE_cutSet2")
h_visE_data_wCutsSet2.Draw("E")
h_visE_all_wCutsSet2.Draw("hist same")
h_visE_data_wCutsSet2.Draw("ESAME")
if args.recoEOverflow:
  label_visE_cutSet2 = getOverflowLabel(h_visE_data_wCutsSet2)
  label_visE_cutSet2.Draw()
leg_visE_cutSet2 = rt.TLegend(0.7,0.7,0.9,0.9)
leg_visE_cutSet2.AddEntry(h_visE_ext_wCutsSet2, "cosmic background (%.2f)"%h_visE_ext_wCutsSet2.Integral(), "f")
leg_visE_cutSet2.AddEntry(h_visE_NCnue_wCutsSet2, "NC nue (%.2f)"%h_visE_NCnue_wCutsSet2.Integral(), "f")
leg_visE_cutSet2.AddEntry(h_visE_CCnue_wCutsSet2, "CC nue (%.2f)"%h_visE_CCnue_wCutsSet2.Integral(), "f")
leg_visE_cutSet2.AddEntry(h_visE_NCnumu_wCutsSet2, "NC numu (%.2f)"%h_visE_NCnumu_wCutsSet2.Integral(), "f")
leg_visE_cutSet2.AddEntry(h_visE_CCnumu_wCutsSet2, "CC numu (%.2f)"%h_visE_CCnumu_wCutsSet2.Integral(), "f")
leg_visE_cutSet2.AddEntry(h_visE_data_wCutsSet2, "5e19 data (%.2f)"%h_visE_data_wCutsSet2.Integral(), "l")
leg_visE_cutSet2.Draw()
cnv_visE_cutSet2.Write()

cnv_visE_cutSet3 = rt.TCanvas("cnv_visE_cutSet3", "cnv_visE_cutSet3")
h_visE_data_wCutsSet3.Draw("E")
h_visE_all_wCutsSet3.Draw("hist same")
h_visE_data_wCutsSet3.Draw("ESAME")
if args.recoEOverflow:
  label_visE_cutSet3 = getOverflowLabel(h_visE_data_wCutsSet3)
  label_visE_cutSet3.Draw()
leg_visE_cutSet3 = rt.TLegend(0.7,0.7,0.9,0.9)
leg_visE_cutSet3.AddEntry(h_visE_ext_wCutsSet3, "cosmic background (%.2f)"%h_visE_ext_wCutsSet3.Integral(), "f")
leg_visE_cutSet3.AddEntry(h_visE_NCnue_wCutsSet3, "NC nue (%.2f)"%h_visE_NCnue_wCutsSet3.Integral(), "f")
leg_visE_cutSet3.AddEntry(h_visE_CCnue_wCutsSet3, "CC nue (%.2f)"%h_visE_CCnue_wCutsSet3.Integral(), "f")
leg_visE_cutSet3.AddEntry(h_visE_NCnumu_wCutsSet3, "NC numu (%.2f)"%h_visE_NCnumu_wCutsSet3.Integral(), "f")
leg_visE_cutSet3.AddEntry(h_visE_CCnumu_wCutsSet3, "CC numu (%.2f)"%h_visE_CCnumu_wCutsSet3.Integral(), "f")
leg_visE_cutSet3.AddEntry(h_visE_data_wCutsSet3, "5e19 data (%.2f)"%h_visE_data_wCutsSet3.Integral(), "l")
leg_visE_cutSet3.Draw()
cnv_visE_cutSet3.Write()

cnv_visE_cutSet4 = rt.TCanvas("cnv_visE_cutSet4", "cnv_visE_cutSet4")
h_visE_data_wCutsSet4.Draw("E")
h_visE_all_wCutsSet4.Draw("hist same")
h_visE_data_wCutsSet4.Draw("ESAME")
if args.recoEOverflow:
  label_visE_cutSet4 = getOverflowLabel(h_visE_data_wCutsSet4)
  label_visE_cutSet4.Draw()
leg_visE_cutSet4 = rt.TLegend(0.7,0.7,0.9,0.9)
leg_visE_cutSet4.AddEntry(h_visE_ext_wCutsSet4, "cosmic background (%.2f)"%h_visE_ext_wCutsSet4.Integral(), "f")
leg_visE_cutSet4.AddEntry(h_visE_NCnue_wCutsSet4, "NC nue (%.2f)"%h_visE_NCnue_wCutsSet4.Integral(), "f")
leg_visE_cutSet4.AddEntry(h_visE_CCnue_wCutsSet4, "CC nue (%.2f)"%h_visE_CCnue_wCutsSet4.Integral(), "f")
leg_visE_cutSet4.AddEntry(h_visE_NCnumu_wCutsSet4, "NC numu (%.2f)"%h_visE_NCnumu_wCutsSet4.Integral(), "f")
leg_visE_cutSet4.AddEntry(h_visE_CCnumu_wCutsSet4, "CC numu (%.2f)"%h_visE_CCnumu_wCutsSet4.Integral(), "f")
leg_visE_cutSet4.AddEntry(h_visE_data_wCutsSet4, "5e19 data (%.2f)"%h_visE_data_wCutsSet4.Integral(), "l")
leg_visE_cutSet4.Draw()
cnv_visE_cutSet4.Write()

cnv_visE_cutSet5 = rt.TCanvas("cnv_visE_cutSet5", "cnv_visE_cutSet5")
h_visE_data_wCutsSet5.Draw("E")
h_visE_all_wCutsSet5.Draw("hist same")
h_visE_data_wCutsSet5.Draw("ESAME")
if args.recoEOverflow:
  label_visE_cutSet5 = getOverflowLabel(h_visE_data_wCutsSet5)
  label_visE_cutSet5.Draw()
leg_visE_cutSet5 = rt.TLegend(0.7,0.7,0.9,0.9)
leg_visE_cutSet5.AddEntry(h_visE_ext_wCutsSet5, "cosmic background (%.2f)"%h_visE_ext_wCutsSet5.Integral(), "f")
leg_visE_cutSet5.AddEntry(h_visE_NCnue_wCutsSet5, "NC nue (%.2f)"%h_visE_NCnue_wCutsSet5.Integral(), "f")
leg_visE_cutSet5.AddEntry(h_visE_CCnue_wCutsSet5, "CC nue (%.2f)"%h_visE_CCnue_wCutsSet5.Integral(), "f")
leg_visE_cutSet5.AddEntry(h_visE_NCnumu_wCutsSet5, "NC numu (%.2f)"%h_visE_NCnumu_wCutsSet5.Integral(), "f")
leg_visE_cutSet5.AddEntry(h_visE_CCnumu_wCutsSet5, "CC numu (%.2f)"%h_visE_CCnumu_wCutsSet5.Integral(), "f")
leg_visE_cutSet5.AddEntry(h_visE_data_wCutsSet5, "5e19 data (%.2f)"%h_visE_data_wCutsSet5.Integral(), "l")
leg_visE_cutSet5.Draw()
cnv_visE_cutSet5.Write()

cnv_visE_cutSet6 = rt.TCanvas("cnv_visE_cutSet6", "cnv_visE_cutSet6")
h_visE_data_wCutsSet6.Draw("E")
h_visE_all_wCutsSet6.Draw("hist same")
h_visE_data_wCutsSet6.Draw("ESAME")
if args.recoEOverflow:
  label_visE_cutSet6 = getOverflowLabel(h_visE_data_wCutsSet6)
  label_visE_cutSet6.Draw()
leg_visE_cutSet6 = rt.TLegend(0.7,0.7,0.9,0.9)
leg_visE_cutSet6.AddEntry(h_visE_ext_wCutsSet6, "cosmic background (%.2f)"%h_visE_ext_wCutsSet6.Integral(), "f")
leg_visE_cutSet6.AddEntry(h_visE_NCnue_wCutsSet6, "NC nue (%.2f)"%h_visE_NCnue_wCutsSet6.Integral(), "f")
leg_visE_cutSet6.AddEntry(h_visE_CCnue_wCutsSet6, "CC nue (%.2f)"%h_visE_CCnue_wCutsSet6.Integral(), "f")
leg_visE_cutSet6.AddEntry(h_visE_NCnumu_wCutsSet6, "NC numu (%.2f)"%h_visE_NCnumu_wCutsSet6.Integral(), "f")
leg_visE_cutSet6.AddEntry(h_visE_CCnumu_wCutsSet6, "CC numu (%.2f)"%h_visE_CCnumu_wCutsSet6.Integral(), "f")
leg_visE_cutSet6.AddEntry(h_visE_data_wCutsSet6, "5e19 data (%.2f)"%h_visE_data_wCutsSet6.Integral(), "l")
leg_visE_cutSet6.Draw()
cnv_visE_cutSet6.Write()

cnv_cosTheta_sel = rt.TCanvas("cnv_cosTheta_sel", "cnv_cosTheta_sel")
h_cosTheta_data_wCuts.Draw("E")
h_cosTheta_all_wCuts.Draw("hist same")
h_cosTheta_data_wCuts.Draw("ESAME")
h_cosTheta_predErr_wCuts.SetFillColor(4)
h_cosTheta_predErr_wCuts.SetFillStyle(3002)
h_cosTheta_predErr_wCuts.Draw("E2SAME")
leg_cosTheta_sel = rt.TLegend(0.7,0.7,0.9,0.9)
leg_cosTheta_sel.AddEntry(h_cosTheta_ext_wCuts, "cosmic background (%.2f)"%h_cosTheta_ext_wCuts.Integral(), "f")
leg_cosTheta_sel.AddEntry(h_cosTheta_NCnue_wCuts, "NC nue (%.2f)"%h_cosTheta_NCnue_wCuts.Integral(), "f")
leg_cosTheta_sel.AddEntry(h_cosTheta_CCnue_wCuts, "CC nue (%.2f)"%h_cosTheta_CCnue_wCuts.Integral(), "f")
leg_cosTheta_sel.AddEntry(h_cosTheta_NCnumu_wCuts, "NC numu (%.2f)"%h_cosTheta_NCnumu_wCuts.Integral(), "f")
leg_cosTheta_sel.AddEntry(h_cosTheta_CCnumu_wCuts, "CC numu (%.2f)"%h_cosTheta_CCnumu_wCuts.Integral(), "f")
leg_cosTheta_sel.AddEntry(h_cosTheta_data_wCuts, "5e19 data (%.2f)"%h_cosTheta_data_wCuts.Integral(), "l")
leg_cosTheta_sel.AddEntry(h_cosTheta_predErr_wCuts, "Pred. stats + sys (no detvar) unc.", "f")
leg_cosTheta_sel.Draw()
cnv_cosTheta_sel.Write()

cnv_lepP_sel = rt.TCanvas("cnv_lepP_sel", "cnv_lepP_sel")
h_lepP_data_wCuts.Draw("E")
h_lepP_all_wCuts.Draw("hist same")
h_lepP_data_wCuts.Draw("ESAME")
h_lepP_predErr_wCuts.SetFillColor(4)
h_lepP_predErr_wCuts.SetFillStyle(3002)
h_lepP_predErr_wCuts.Draw("E2SAME")
if args.recoEOverflow:
  label_lepP_sel = getOverflowLabel(h_lepP_data_wCuts)
  label_lepP_sel.Draw()
leg_lepP_sel = rt.TLegend(0.7,0.7,0.9,0.9)
leg_lepP_sel.AddEntry(h_lepP_ext_wCuts, "cosmic background (%.2f)"%h_lepP_ext_wCuts.Integral(), "f")
leg_lepP_sel.AddEntry(h_lepP_NCnue_wCuts, "NC nue (%.2f)"%h_lepP_NCnue_wCuts.Integral(), "f")
leg_lepP_sel.AddEntry(h_lepP_CCnue_wCuts, "CC nue (%.2f)"%h_lepP_CCnue_wCuts.Integral(), "f")
leg_lepP_sel.AddEntry(h_lepP_NCnumu_wCuts, "NC numu (%.2f)"%h_lepP_NCnumu_wCuts.Integral(), "f")
leg_lepP_sel.AddEntry(h_lepP_CCnumu_wCuts, "CC numu (%.2f)"%h_lepP_CCnumu_wCuts.Integral(), "f")
leg_lepP_sel.AddEntry(h_lepP_data_wCuts, "5e19 data (%.2f)"%h_lepP_data_wCuts.Integral(), "l")
leg_lepP_sel.AddEntry(h_lepP_predErr_wCuts, "Pred. stats + sys (no detvar) unc.", "f")
leg_lepP_sel.Draw()
cnv_lepP_sel.Write()

cnv_muScr_sel = rt.TCanvas("cnv_muScr_sel", "cnv_muScr_sel")
h_muScr_data_wCuts.Draw("E")
h_muScr_all_wCuts.Draw("hist same")
h_muScr_data_wCuts.Draw("ESAME")
leg_muScr_sel = rt.TLegend(0.7,0.7,0.9,0.9)
leg_muScr_sel.AddEntry(h_muScr_ext_wCuts, "cosmic background (%.2f)"%h_muScr_ext_wCuts.Integral(), "f")
leg_muScr_sel.AddEntry(h_muScr_NCnue_wCuts, "NC nue (%.2f)"%h_muScr_NCnue_wCuts.Integral(), "f")
leg_muScr_sel.AddEntry(h_muScr_CCnue_wCuts, "CC nue (%.2f)"%h_muScr_CCnue_wCuts.Integral(), "f")
leg_muScr_sel.AddEntry(h_muScr_NCnumu_wCuts, "NC numu (%.2f)"%h_muScr_NCnumu_wCuts.Integral(), "f")
leg_muScr_sel.AddEntry(h_muScr_CCnumu_wCuts, "CC numu (%.2f)"%h_muScr_CCnumu_wCuts.Integral(), "f")
leg_muScr_sel.AddEntry(h_muScr_data_wCuts, "5e19 data (%.2f)"%h_muScr_data_wCuts.Integral(), "l")
leg_muScr_sel.Draw()
cnv_muScr_sel.Write()

cnv_piScr_sel = rt.TCanvas("cnv_piScr_sel", "cnv_piScr_sel")
h_piScr_data_wCuts.Draw("E")
h_piScr_all_wCuts.Draw("hist same")
h_piScr_data_wCuts.Draw("ESAME")
leg_piScr_sel = rt.TLegend(0.7,0.7,0.9,0.9)
leg_piScr_sel.AddEntry(h_piScr_ext_wCuts, "cosmic background (%.2f)"%h_piScr_ext_wCuts.Integral(), "f")
leg_piScr_sel.AddEntry(h_piScr_NCnue_wCuts, "NC nue (%.2f)"%h_piScr_NCnue_wCuts.Integral(), "f")
leg_piScr_sel.AddEntry(h_piScr_CCnue_wCuts, "CC nue (%.2f)"%h_piScr_CCnue_wCuts.Integral(), "f")
leg_piScr_sel.AddEntry(h_piScr_NCnumu_wCuts, "NC numu (%.2f)"%h_piScr_NCnumu_wCuts.Integral(), "f")
leg_piScr_sel.AddEntry(h_piScr_CCnumu_wCuts, "CC numu (%.2f)"%h_piScr_CCnumu_wCuts.Integral(), "f")
leg_piScr_sel.AddEntry(h_piScr_data_wCuts, "5e19 data (%.2f)"%h_piScr_data_wCuts.Integral(), "l")
leg_piScr_sel.Draw()
cnv_piScr_sel.Write()

cnv_prScr_sel = rt.TCanvas("cnv_prScr_sel", "cnv_prScr_sel")
h_prScr_data_wCuts.Draw("E")
h_prScr_all_wCuts.Draw("hist same")
h_prScr_data_wCuts.Draw("ESAME")
leg_prScr_sel = rt.TLegend(0.7,0.7,0.9,0.9)
leg_prScr_sel.AddEntry(h_prScr_ext_wCuts, "cosmic background (%.2f)"%h_prScr_ext_wCuts.Integral(), "f")
leg_prScr_sel.AddEntry(h_prScr_NCnue_wCuts, "NC nue (%.2f)"%h_prScr_NCnue_wCuts.Integral(), "f")
leg_prScr_sel.AddEntry(h_prScr_CCnue_wCuts, "CC nue (%.2f)"%h_prScr_CCnue_wCuts.Integral(), "f")
leg_prScr_sel.AddEntry(h_prScr_NCnumu_wCuts, "NC numu (%.2f)"%h_prScr_NCnumu_wCuts.Integral(), "f")
leg_prScr_sel.AddEntry(h_prScr_CCnumu_wCuts, "CC numu (%.2f)"%h_prScr_CCnumu_wCuts.Integral(), "f")
leg_prScr_sel.AddEntry(h_prScr_data_wCuts, "5e19 data (%.2f)"%h_prScr_data_wCuts.Integral(), "l")
leg_prScr_sel.Draw()
cnv_prScr_sel.Write()

cnv_elScr_sel = rt.TCanvas("cnv_elScr_sel", "cnv_elScr_sel")
h_elScr_data_wCuts.Draw("E")
h_elScr_all_wCuts.Draw("hist same")
h_elScr_data_wCuts.Draw("ESAME")
leg_elScr_sel = rt.TLegend(0.7,0.7,0.9,0.9)
leg_elScr_sel.AddEntry(h_elScr_ext_wCuts, "cosmic background (%.2f)"%h_elScr_ext_wCuts.Integral(), "f")
leg_elScr_sel.AddEntry(h_elScr_NCnue_wCuts, "NC nue (%.2f)"%h_elScr_NCnue_wCuts.Integral(), "f")
leg_elScr_sel.AddEntry(h_elScr_CCnue_wCuts, "CC nue (%.2f)"%h_elScr_CCnue_wCuts.Integral(), "f")
leg_elScr_sel.AddEntry(h_elScr_NCnumu_wCuts, "NC numu (%.2f)"%h_elScr_NCnumu_wCuts.Integral(), "f")
leg_elScr_sel.AddEntry(h_elScr_CCnumu_wCuts, "CC numu (%.2f)"%h_elScr_CCnumu_wCuts.Integral(), "f")
leg_elScr_sel.AddEntry(h_elScr_data_wCuts, "5e19 data (%.2f)"%h_elScr_data_wCuts.Integral(), "l")
leg_elScr_sel.Draw()
cnv_elScr_sel.Write()

cnv_phScr_sel = rt.TCanvas("cnv_phScr_sel", "cnv_phScr_sel")
h_phScr_data_wCuts.Draw("E")
h_phScr_all_wCuts.Draw("hist same")
h_phScr_data_wCuts.Draw("ESAME")
leg_phScr_sel = rt.TLegend(0.7,0.7,0.9,0.9)
leg_phScr_sel.AddEntry(h_phScr_ext_wCuts, "cosmic background (%.2f)"%h_phScr_ext_wCuts.Integral(), "f")
leg_phScr_sel.AddEntry(h_phScr_NCnue_wCuts, "NC nue (%.2f)"%h_phScr_NCnue_wCuts.Integral(), "f")
leg_phScr_sel.AddEntry(h_phScr_CCnue_wCuts, "CC nue (%.2f)"%h_phScr_CCnue_wCuts.Integral(), "f")
leg_phScr_sel.AddEntry(h_phScr_NCnumu_wCuts, "NC numu (%.2f)"%h_phScr_NCnumu_wCuts.Integral(), "f")
leg_phScr_sel.AddEntry(h_phScr_CCnumu_wCuts, "CC numu (%.2f)"%h_phScr_CCnumu_wCuts.Integral(), "f")
leg_phScr_sel.AddEntry(h_phScr_data_wCuts, "5e19 data (%.2f)"%h_phScr_data_wCuts.Integral(), "l")
leg_phScr_sel.Draw()
cnv_phScr_sel.Write()

cnv_pPScr_sel = rt.TCanvas("cnv_pPScr_sel", "cnv_pPScr_sel")
h_pPScr_data_wCuts.Draw("E")
h_pPScr_all_wCuts.Draw("hist same")
h_pPScr_data_wCuts.Draw("ESAME")
leg_pPScr_sel = rt.TLegend(0.7,0.7,0.9,0.9)
leg_pPScr_sel.AddEntry(h_pPScr_ext_wCuts, "cosmic background (%.2f)"%h_pPScr_ext_wCuts.Integral(), "f")
leg_pPScr_sel.AddEntry(h_pPScr_NCnue_wCuts, "NC nue (%.2f)"%h_pPScr_NCnue_wCuts.Integral(), "f")
leg_pPScr_sel.AddEntry(h_pPScr_CCnue_wCuts, "CC nue (%.2f)"%h_pPScr_CCnue_wCuts.Integral(), "f")
leg_pPScr_sel.AddEntry(h_pPScr_NCnumu_wCuts, "NC numu (%.2f)"%h_pPScr_NCnumu_wCuts.Integral(), "f")
leg_pPScr_sel.AddEntry(h_pPScr_CCnumu_wCuts, "CC numu (%.2f)"%h_pPScr_CCnumu_wCuts.Integral(), "f")
leg_pPScr_sel.AddEntry(h_pPScr_data_wCuts, "5e19 data (%.2f)"%h_pPScr_data_wCuts.Integral(), "l")
leg_pPScr_sel.Draw()
cnv_pPScr_sel.Write()

cnv_pNScr_sel = rt.TCanvas("cnv_pNScr_sel", "cnv_pNScr_sel")
h_pNScr_data_wCuts.Draw("E")
h_pNScr_all_wCuts.Draw("hist same")
h_pNScr_data_wCuts.Draw("ESAME")
leg_pNScr_sel = rt.TLegend(0.7,0.7,0.9,0.9)
leg_pNScr_sel.AddEntry(h_pNScr_ext_wCuts, "cosmic background (%.2f)"%h_pNScr_ext_wCuts.Integral(), "f")
leg_pNScr_sel.AddEntry(h_pNScr_NCnue_wCuts, "NC nue (%.2f)"%h_pNScr_NCnue_wCuts.Integral(), "f")
leg_pNScr_sel.AddEntry(h_pNScr_CCnue_wCuts, "CC nue (%.2f)"%h_pNScr_CCnue_wCuts.Integral(), "f")
leg_pNScr_sel.AddEntry(h_pNScr_NCnumu_wCuts, "NC numu (%.2f)"%h_pNScr_NCnumu_wCuts.Integral(), "f")
leg_pNScr_sel.AddEntry(h_pNScr_CCnumu_wCuts, "CC numu (%.2f)"%h_pNScr_CCnumu_wCuts.Integral(), "f")
leg_pNScr_sel.AddEntry(h_pNScr_data_wCuts, "5e19 data (%.2f)"%h_pNScr_data_wCuts.Integral(), "l")
leg_pNScr_sel.Draw()
cnv_pNScr_sel.Write()

cnv_pCScr_sel = rt.TCanvas("cnv_pCScr_sel", "cnv_pCScr_sel")
h_pCScr_data_wCuts.Draw("E")
h_pCScr_all_wCuts.Draw("hist same")
h_pCScr_data_wCuts.Draw("ESAME")
leg_pCScr_sel = rt.TLegend(0.7,0.7,0.9,0.9)
leg_pCScr_sel.AddEntry(h_pCScr_ext_wCuts, "cosmic background (%.2f)"%h_pCScr_ext_wCuts.Integral(), "f")
leg_pCScr_sel.AddEntry(h_pCScr_NCnue_wCuts, "NC nue (%.2f)"%h_pCScr_NCnue_wCuts.Integral(), "f")
leg_pCScr_sel.AddEntry(h_pCScr_CCnue_wCuts, "CC nue (%.2f)"%h_pCScr_CCnue_wCuts.Integral(), "f")
leg_pCScr_sel.AddEntry(h_pCScr_NCnumu_wCuts, "NC numu (%.2f)"%h_pCScr_NCnumu_wCuts.Integral(), "f")
leg_pCScr_sel.AddEntry(h_pCScr_CCnumu_wCuts, "CC numu (%.2f)"%h_pCScr_CCnumu_wCuts.Integral(), "f")
leg_pCScr_sel.AddEntry(h_pCScr_data_wCuts, "5e19 data (%.2f)"%h_pCScr_data_wCuts.Integral(), "l")
leg_pCScr_sel.Draw()
cnv_pCScr_sel.Write()

cnv_lowEExVsX_cutSet1 = rt.TCanvas("cnv_lowEExVsX_cutSet1", "cnv_lowEExVsX_cutSet1")
h_lowEExVsX_wCutsSet1.Draw("E")
cnv_lowEExVsX_cutSet1.Write()

cnv_lowEExVsX_cutSet2 = rt.TCanvas("cnv_lowEExVsX_cutSet2", "cnv_lowEExVsX_cutSet2")
h_lowEExVsX_wCutsSet2.Draw("E")
cnv_lowEExVsX_cutSet2.Write()

cnv_lowEExVsX_cutSet3 = rt.TCanvas("cnv_lowEExVsX_cutSet3", "cnv_lowEExVsX_cutSet3")
h_lowEExVsX_wCutsSet3.Draw("E")
cnv_lowEExVsX_cutSet3.Write()

cnv_lowEExVsX_cutSet4 = rt.TCanvas("cnv_lowEExVsX_cutSet4", "cnv_lowEExVsX_cutSet4")
h_lowEExVsX_wCutsSet4.Draw("E")
cnv_lowEExVsX_cutSet4.Write()

cnv_lowEExVsX_cutSet5 = rt.TCanvas("cnv_lowEExVsX_cutSet5", "cnv_lowEExVsX_cutSet5")
h_lowEExVsX_wCutsSet5.Draw("E")
cnv_lowEExVsX_cutSet5.Write()

cnv_lowEExVsX_cutSet6 = rt.TCanvas("cnv_lowEExVsX_cutSet6", "cnv_lowEExVsX_cutSet6")
h_lowEExVsX_wCutsSet6.Draw("E")
cnv_lowEExVsX_cutSet6.Write()

cnv_lowEExVsX_allCuts = rt.TCanvas("cnv_lowEExVsX_allCuts", "cnv_lowEExVsX_allCuts")
h_lowEExVsX_wCuts.Draw("E")
cnv_lowEExVsX_allCuts.Write()


h_nuE_CCnumu_eff1.GetYaxis().SetRangeUser(0,1.003)
h_nuE_CCnumu_eff2.GetYaxis().SetRangeUser(0,1.003)
h_nuE_CCnumu_effAll.GetYaxis().SetRangeUser(0,1.003)
cnv_cutSets_eff = rt.TCanvas("cnv_cutSets_eff","cnv_cutSets_eff")
cnv_cutSets_eff.SetGrid()
h_nuE_CCnumu_eff1.Draw("E")
h_nuE_CCnumu_eff2.Draw("ESAME")
h_nuE_CCnumu_effAll.Draw("ESAME")
leg_cutSets_eff = rt.TLegend(0.7,0.7,0.9,0.9)
leg_cutSets_eff.AddEntry(h_nuE_CCnumu_eff1, "cut set 1", "l")
leg_cutSets_eff.AddEntry(h_nuE_CCnumu_eff2, "cut set 2", "l")
leg_cutSets_eff.AddEntry(h_nuE_CCnumu_effAll, "all cuts", "l")
leg_cutSets_eff.Draw()
cnv_cutSets_eff.Write()

h_nuE_CCnumu_eff1.Write()
h_nuE_CCnumu_eff2.Write()
h_nuE_CCnumu_effAll.Write()


if args.write_ntuples:

  fnu_trimmed.cd()
  tnu_trimmed.Write("",rt.TObject.kOverwrite)
  tnuPOT.CloneTree().Write()
  fnu_trimmed.Close()
  
  fnue_trimmed.cd()
  tnue_trimmed.Write("",rt.TObject.kOverwrite)
  tnuePOT.CloneTree().Write()
  fnue_trimmed.Close()
  
  fext_trimmed.cd()
  text_trimmed.Write("",rt.TObject.kOverwrite)
  fext_trimmed.Close()
  
  fdata_trimmed.cd()
  tdata_trimmed.Write("",rt.TObject.kOverwrite)
  fdata_trimmed.Close()



