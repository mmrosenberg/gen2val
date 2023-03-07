
import sys, argparse
import numpy as np
import ROOT as rt

from math import isinf
from helpers.plotting_functions import sortHists


parser = argparse.ArgumentParser("Plot Selection Test Results")
parser.add_argument("-fnu", "--bnbnu_file", type=str, default="prepare_selection_test_output/prepare_selection_test_reco_v2me05_gen2val_v16_nu_file.root", help="bnb nu input file")
parser.add_argument("-fnue", "--bnbnue_file", type=str, default="prepare_selection_test_output/prepare_selection_test_reco_v2me05_gen2val_v16_nue_file.root", help="bnb nu input file")
parser.add_argument("-fext", "--extbnb_file", type=str, default="prepare_selection_test_output/prepare_selection_test_reco_v2me05_gen2val_v16_extbnb_file.root", help="bnb nu input file")
parser.add_argument("-c", "--compCut", type=float, default=0.6, help="completeness cut value")
parser.add_argument("-s", "--confCut", type=float, default=9, help="electron class confidence cut value")
parser.add_argument("-q", "--chargeCut", type=float, default=70000, help="electron charge fraction cut value")
parser.add_argument("-qf", "--chargeFracCut", type=float, default=0.7, help="electron charge fraction cut value")
parser.add_argument("-o", "--outfile", type=str, default="plot_selection_test_results_output.root", help="output root file name")
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

run3POT = 4.3e+19 + 1.701e+20 + 2.97e+19 + 1.524e+17
runs1to3POT = 6.67e+20

tnuePOTsum = 0.
for i in range(tnuePOT.GetEntries()):
  tnuePOT.GetEntry(i)
  tnuePOTsum = tnuePOTsum + tnuePOT.totGoodPOT

tnuPOTsum = 0.
for i in range(tnuPOT.GetEntries()):
  tnuPOT.GetEntry(i)
  tnuPOTsum = tnuPOTsum + tnuPOT.totGoodPOT

#970: number of run3 extbnb merged_dlana files in prepare_selection_test_reco_v2me05_gen2val_v16_extbnb_file.root
#89559: number of run3 extbnb files (# in def: prod_extunbiased_swizzle_crt_inclusive_v6_v6a_goodruns_mcc9_run3)
textPOTsum = (970./89559.)*run3POT

def configureHists(h_CCnumu, h_NCnumu, h_CCnue, h_NCnue, h_ext):
  h_CCnumu.GetYaxis().SetTitle("events per 6.67e+20 POT")
  h_CCnumu.SetLineColor(rt.kBlue)
  h_CCnumu.SetLineWidth(2)
  h_NCnumu.GetYaxis().SetTitle("events per 6.67e+20 POT")
  h_NCnumu.SetLineColor(40)
  h_NCnumu.SetLineWidth(2)
  h_CCnue.GetYaxis().SetTitle("events per 6.67e+20 POT")
  h_CCnue.SetLineColor(rt.kRed)
  h_CCnue.SetLineWidth(2)
  h_NCnue.GetYaxis().SetTitle("events per 6.67e+20 POT")
  h_NCnue.SetLineColor(8)
  h_NCnue.SetLineWidth(2)
  h_ext.GetYaxis().SetTitle("events per 6.67e+20 POT")
  h_ext.SetLineColor(rt.kBlack)
  h_ext.SetLineWidth(2)
  return h_CCnumu, h_NCnumu, h_CCnue, h_NCnue, h_ext


h_cosFrac_CCnumu = rt.TH1F("h_cosFrac_CCnumu","Fraction of Hits On Cosmics",42,-1.05,1.05)
h_cosFrac_NCnumu = rt.TH1F("h_cosFrac_NCnumu","Fraction of Hits On Cosmics",42,-1.05,1.05)
h_cosFrac_CCnue = rt.TH1F("h_cosFrac_CCnue","Fraction of Hits On Cosmics",42,-1.05,1.05)
h_cosFrac_NCnue = rt.TH1F("h_cosFrac_NCnue","Fraction of Hits On Cosmics",42,-1.05,1.05)
h_cosFrac_ext = rt.TH1F("h_cosFrac_ext","Fraction of Hits On Cosmics",42,-1.05,1.05)
h_cosFrac_CCnumu, h_cosFrac_NCnumu, h_cosFrac_CCnue, h_cosFrac_NCnue, h_cosFrac_ext = configureHists(h_cosFrac_CCnumu,
 h_cosFrac_NCnumu, h_cosFrac_CCnue, h_cosFrac_NCnue, h_cosFrac_ext)

h_nEl_CCnumu = rt.TH1F("h_nEl_CCnumu","Number of Reco Electrons",10,0,10)
h_nEl_NCnumu = rt.TH1F("h_nEl_NCnumu","Number of Reco Electrons",10,0,10)
h_nEl_CCnue = rt.TH1F("h_nEl_CCnue","Number of Reco Electrons",10,0,10)
h_nEl_NCnue = rt.TH1F("h_nEl_NCnue","Number of Reco Electrons",10,0,10)
h_nEl_ext = rt.TH1F("h_nEl_ext","Number of Reco Electrons",10,0,10)
h_nEl_CCnumu, h_nEl_NCnumu, h_nEl_CCnue, h_nEl_NCnue, h_nEl_ext = configureHists(h_nEl_CCnumu,
 h_nEl_NCnumu, h_nEl_CCnue, h_nEl_NCnue, h_nEl_ext)

h_nCompEl_CCnumu = rt.TH1F("h_nCompEl_CCnumu","Number of Complete (c > %.2f) Reco Electrons"%args.compCut,10,0,10)
h_nCompEl_NCnumu = rt.TH1F("h_nCompEl_NCnumu","Number of Complete (c > %.2f) Reco Electrons"%args.compCut,10,0,10)
h_nCompEl_CCnue = rt.TH1F("h_nCompEl_CCnue","Number of Complete (c > %.2f) Reco Electrons"%args.compCut,10,0,10)
h_nCompEl_NCnue = rt.TH1F("h_nCompEl_NCnue","Number of Complete (c > %.2f) Reco Electrons"%args.compCut,10,0,10)
h_nCompEl_ext = rt.TH1F("h_nCompEl_ext","Number of Complete (c > %.2f) Reco Electrons"%args.compCut,10,0,10)
h_nCompEl_CCnumu, h_nCompEl_NCnumu, h_nCompEl_CCnue, h_nCompEl_NCnue, h_nCompEl_ext = configureHists(h_nCompEl_CCnumu,
 h_nCompEl_NCnumu, h_nCompEl_CCnue, h_nCompEl_NCnue, h_nCompEl_ext)

h_nMu_CCnumu = rt.TH1F("h_nMu_CCnumu","Number of Reco Muons",10,0,10)
h_nMu_NCnumu = rt.TH1F("h_nMu_NCnumu","Number of Reco Muons",10,0,10)
h_nMu_CCnue = rt.TH1F("h_nMu_CCnue","Number of Reco Muons",10,0,10)
h_nMu_NCnue = rt.TH1F("h_nMu_NCnue","Number of Reco Muons",10,0,10)
h_nMu_ext = rt.TH1F("h_nMu_ext","Number of Reco Muons",10,0,10)
h_nMu_CCnumu, h_nMu_NCnumu, h_nMu_CCnue, h_nMu_NCnue, h_nMu_ext = configureHists(h_nMu_CCnumu,
 h_nMu_NCnumu, h_nMu_CCnue, h_nMu_NCnue, h_nMu_ext)

h_nCompMu_CCnumu = rt.TH1F("h_nCompMu_CCnumu","Number of Complete (c > %.2f) Reco Muons"%args.compCut,10,0,10)
h_nCompMu_NCnumu = rt.TH1F("h_nCompMu_NCnumu","Number of Complete (c > %.2f) Reco Muons"%args.compCut,10,0,10)
h_nCompMu_CCnue = rt.TH1F("h_nCompMu_CCnue","Number of Complete (c > %.2f) Reco Muons"%args.compCut,10,0,10)
h_nCompMu_NCnue = rt.TH1F("h_nCompMu_NCnue","Number of Complete (c > %.2f) Reco Muons"%args.compCut,10,0,10)
h_nCompMu_ext = rt.TH1F("h_nCompMu_ext","Number of Complete (c > %.2f) Reco Muons"%args.compCut,10,0,10)
h_nCompMu_CCnumu, h_nCompMu_NCnumu, h_nCompMu_CCnue, h_nCompMu_NCnue, h_nCompMu_ext = configureHists(h_nCompMu_CCnumu,
 h_nCompMu_NCnumu, h_nCompMu_CCnue, h_nCompMu_NCnue, h_nCompMu_ext)

h_elMaxQF_CCnumu = rt.TH1F("h_elMaxQF_CCnumu","Max Electron Charge Fraction",21,0,1.05)
h_elMaxQF_NCnumu = rt.TH1F("h_elMaxQF_NCnumu","Max Electron Charge Fraction",21,0,1.05)
h_elMaxQF_CCnue = rt.TH1F("h_elMaxQF_CCnue","Max Electron Charge Fraction",21,0,1.05)
h_elMaxQF_NCnue = rt.TH1F("h_elMaxQF_NCnue","Max Electron Charge Fraction",21,0,1.05)
h_elMaxQF_ext = rt.TH1F("h_elMaxQF_ext","Max Electron Charge Fraction",21,0,1.05)
h_elMaxQF_CCnumu, h_elMaxQF_NCnumu, h_elMaxQF_CCnue, h_elMaxQF_NCnue, h_elMaxQF_ext = configureHists(h_elMaxQF_CCnumu,
 h_elMaxQF_NCnumu, h_elMaxQF_CCnue, h_elMaxQF_NCnue, h_elMaxQF_ext)

h_muMaxQF_CCnumu = rt.TH1F("h_muMaxQF_CCnumu","Max Muon Charge Fraction",21,0,1.05)
h_muMaxQF_NCnumu = rt.TH1F("h_muMaxQF_NCnumu","Max Muon Charge Fraction",21,0,1.05)
h_muMaxQF_CCnue = rt.TH1F("h_muMaxQF_CCnue","Max Muon Charge Fraction",21,0,1.05)
h_muMaxQF_NCnue = rt.TH1F("h_muMaxQF_NCnue","Max Muon Charge Fraction",21,0,1.05)
h_muMaxQF_ext = rt.TH1F("h_muMaxQF_ext","Max Muon Charge Fraction",21,0,1.05)
h_muMaxQF_CCnumu, h_muMaxQF_NCnumu, h_muMaxQF_CCnue, h_muMaxQF_NCnue, h_muMaxQF_ext = configureHists(h_muMaxQF_CCnumu,
 h_muMaxQF_NCnumu, h_muMaxQF_CCnue, h_muMaxQF_NCnue, h_muMaxQF_ext)

h_elMaxComp_CCnumu = rt.TH1F("h_elMaxComp_CCnumu","Max Electron Completeness",21,0,1.05)
h_elMaxComp_NCnumu = rt.TH1F("h_elMaxComp_NCnumu","Max Electron Completeness",21,0,1.05)
h_elMaxComp_CCnue = rt.TH1F("h_elMaxComp_CCnue","Max Electron Completeness",21,0,1.05)
h_elMaxComp_NCnue = rt.TH1F("h_elMaxComp_NCnue","Max Electron Completeness",21,0,1.05)
h_elMaxComp_ext = rt.TH1F("h_elMaxComp_ext","Max Electron Completeness",21,0,1.05)
h_elMaxComp_CCnumu, h_elMaxComp_NCnumu, h_elMaxComp_CCnue, h_elMaxComp_NCnue, h_elMaxComp_ext = configureHists(h_elMaxComp_CCnumu,
 h_elMaxComp_NCnumu, h_elMaxComp_CCnue, h_elMaxComp_NCnue, h_elMaxComp_ext)

h_muMaxComp_CCnumu = rt.TH1F("h_muMaxComp_CCnumu","Max Muon Completeness",21,0,1.05)
h_muMaxComp_NCnumu = rt.TH1F("h_muMaxComp_NCnumu","Max Muon Completeness",21,0,1.05)
h_muMaxComp_CCnue = rt.TH1F("h_muMaxComp_CCnue","Max Muon Completeness",21,0,1.05)
h_muMaxComp_NCnue = rt.TH1F("h_muMaxComp_NCnue","Max Muon Completeness",21,0,1.05)
h_muMaxComp_ext = rt.TH1F("h_muMaxComp_ext","Max Muon Completeness",21,0,1.05)
h_muMaxComp_CCnumu, h_muMaxComp_NCnumu, h_muMaxComp_CCnue, h_muMaxComp_NCnue, h_muMaxComp_ext = configureHists(h_muMaxComp_CCnumu,
 h_muMaxComp_NCnumu, h_muMaxComp_CCnue, h_muMaxComp_NCnue, h_muMaxComp_ext)

h_elMaxQComp_CCnumu = rt.TH1F("h_elMaxQComp_CCnumu","Completeness for Largest Electron Shower",21,0,1.05)
h_elMaxQComp_NCnumu = rt.TH1F("h_elMaxQComp_NCnumu","Completeness for Largest Electron Shower",21,0,1.05)
h_elMaxQComp_CCnue = rt.TH1F("h_elMaxQComp_CCnue","Completeness for Largest Electron Shower",21,0,1.05)
h_elMaxQComp_NCnue = rt.TH1F("h_elMaxQComp_NCnue","Completeness for Largest Electron Shower",21,0,1.05)
h_elMaxQComp_ext = rt.TH1F("h_elMaxQComp_ext","Completeness for Largest Electron Shower",21,0,1.05)
h_elMaxQComp_CCnumu, h_elMaxQComp_NCnumu, h_elMaxQComp_CCnue, h_elMaxQComp_NCnue, h_elMaxQComp_ext = configureHists(h_elMaxQComp_CCnumu,
 h_elMaxQComp_NCnumu, h_elMaxQComp_CCnue, h_elMaxQComp_NCnue, h_elMaxQComp_ext)

h_muMaxQComp_CCnumu = rt.TH1F("h_muMaxQComp_CCnumu","Completeness for Largest Muon Track",21,0,1.05)
h_muMaxQComp_NCnumu = rt.TH1F("h_muMaxQComp_NCnumu","Completeness for Largest Muon Track",21,0,1.05)
h_muMaxQComp_CCnue = rt.TH1F("h_muMaxQComp_CCnue","Completeness for Largest Muon Track",21,0,1.05)
h_muMaxQComp_NCnue = rt.TH1F("h_muMaxQComp_NCnue","Completeness for Largest Muon Track",21,0,1.05)
h_muMaxQComp_ext = rt.TH1F("h_muMaxQComp_ext","Completeness for Largest Muon Track",21,0,1.05)
h_muMaxQComp_CCnumu, h_muMaxQComp_NCnumu, h_muMaxQComp_CCnue, h_muMaxQComp_NCnue, h_muMaxQComp_ext = configureHists(h_muMaxQComp_CCnumu,
 h_muMaxQComp_NCnumu, h_muMaxQComp_CCnue, h_muMaxQComp_NCnue, h_muMaxQComp_ext)

h_elMaxQCosTheta_CCnumu = rt.TH1F("h_elMaxQCosTheta_CCnumu","cos(Angle To Beam) for Largest Electron Shower",42,-1.05,1.05)
h_elMaxQCosTheta_NCnumu = rt.TH1F("h_elMaxQCosTheta_NCnumu","cos(Angle To Beam) for Largest Electron Shower",42,-1.05,1.05)
h_elMaxQCosTheta_CCnue = rt.TH1F("h_elMaxQCosTheta_CCnue","cos(Angle To Beam) for Largest Electron Shower",42,-1.05,1.05)
h_elMaxQCosTheta_NCnue = rt.TH1F("h_elMaxQCosTheta_NCnue","cos(Angle To Beam) for Largest Electron Shower",42,-1.05,1.05)
h_elMaxQCosTheta_ext = rt.TH1F("h_elMaxQCosTheta_ext","cos(Angle To Beam) for Largest Electron Shower",42,-1.05,1.05)
h_elMaxQCosTheta_CCnumu, h_elMaxQCosTheta_NCnumu, h_elMaxQCosTheta_CCnue, h_elMaxQCosTheta_NCnue, h_elMaxQCosTheta_ext = configureHists(h_elMaxQCosTheta_CCnumu,
 h_elMaxQCosTheta_NCnumu, h_elMaxQCosTheta_CCnue, h_elMaxQCosTheta_NCnue, h_elMaxQCosTheta_ext)

h_muMaxQCosTheta_CCnumu = rt.TH1F("h_muMaxQCosTheta_CCnumu","cos(Angle To Beam) for Largest Muon Track",42,-1.05,1.05)
h_muMaxQCosTheta_NCnumu = rt.TH1F("h_muMaxQCosTheta_NCnumu","cos(Angle To Beam) for Largest Muon Track",42,-1.05,1.05)
h_muMaxQCosTheta_CCnue = rt.TH1F("h_muMaxQCosTheta_CCnue","cos(Angle To Beam) for Largest Muon Track",42,-1.05,1.05)
h_muMaxQCosTheta_NCnue = rt.TH1F("h_muMaxQCosTheta_NCnue","cos(Angle To Beam) for Largest Muon Track",42,-1.05,1.05)
h_muMaxQCosTheta_ext = rt.TH1F("h_muMaxQCosTheta_ext","cos(Angle To Beam) for Largest Muon Track",42,-1.05,1.05)
h_muMaxQCosTheta_CCnumu, h_muMaxQCosTheta_NCnumu, h_muMaxQCosTheta_CCnue, h_muMaxQCosTheta_NCnue, h_muMaxQCosTheta_ext = configureHists(h_muMaxQCosTheta_CCnumu,
 h_muMaxQCosTheta_NCnumu, h_muMaxQCosTheta_CCnue, h_muMaxQCosTheta_NCnue, h_muMaxQCosTheta_ext)

h_elMaxQElScore_CCnumu = rt.TH1F("h_elMaxQElScore_CCnumu","Electron Score for Largest Electron Shower",21,-20,1)
h_elMaxQElScore_NCnumu = rt.TH1F("h_elMaxQElScore_NCnumu","Electron Score for Largest Electron Shower",21,-20,1)
h_elMaxQElScore_CCnue = rt.TH1F("h_elMaxQElScore_CCnue","Electron Score for Largest Electron Shower",21,-20,1)
h_elMaxQElScore_NCnue = rt.TH1F("h_elMaxQElScore_NCnue","Electron Score for Largest Electron Shower",21,-20,1)
h_elMaxQElScore_ext = rt.TH1F("h_elMaxQElScore_ext","Electron Score for Largest Electron Shower",21,-20,1)
h_elMaxQElScore_CCnumu, h_elMaxQElScore_NCnumu, h_elMaxQElScore_CCnue, h_elMaxQElScore_NCnue, h_elMaxQElScore_ext = configureHists(h_elMaxQElScore_CCnumu,
 h_elMaxQElScore_NCnumu, h_elMaxQElScore_CCnue, h_elMaxQElScore_NCnue, h_elMaxQElScore_ext)

h_elMaxQPhScore_CCnumu = rt.TH1F("h_elMaxQPhScore_CCnumu","Photon Score for Largest Electron Shower",21,-20,1)
h_elMaxQPhScore_NCnumu = rt.TH1F("h_elMaxQPhScore_NCnumu","Photon Score for Largest Electron Shower",21,-20,1)
h_elMaxQPhScore_CCnue = rt.TH1F("h_elMaxQPhScore_CCnue","Photon Score for Largest Electron Shower",21,-20,1)
h_elMaxQPhScore_NCnue = rt.TH1F("h_elMaxQPhScore_NCnue","Photon Score for Largest Electron Shower",21,-20,1)
h_elMaxQPhScore_ext = rt.TH1F("h_elMaxQPhScore_ext","Photon Score for Largest Electron Shower",21,-20,1)
h_elMaxQPhScore_CCnumu, h_elMaxQPhScore_NCnumu, h_elMaxQPhScore_CCnue, h_elMaxQPhScore_NCnue, h_elMaxQPhScore_ext = configureHists(h_elMaxQPhScore_CCnumu,
 h_elMaxQPhScore_NCnumu, h_elMaxQPhScore_CCnue, h_elMaxQPhScore_NCnue, h_elMaxQPhScore_ext)

h_elMaxQPiScore_CCnumu = rt.TH1F("h_elMaxQPiScore_CCnumu","Pion Score for Largest Electron Shower",21,-20,1)
h_elMaxQPiScore_NCnumu = rt.TH1F("h_elMaxQPiScore_NCnumu","Pion Score for Largest Electron Shower",21,-20,1)
h_elMaxQPiScore_CCnue = rt.TH1F("h_elMaxQPiScore_CCnue","Pion Score for Largest Electron Shower",21,-20,1)
h_elMaxQPiScore_NCnue = rt.TH1F("h_elMaxQPiScore_NCnue","Pion Score for Largest Electron Shower",21,-20,1)
h_elMaxQPiScore_ext = rt.TH1F("h_elMaxQPiScore_ext","Pion Score for Largest Electron Shower",21,-20,1)
h_elMaxQPiScore_CCnumu, h_elMaxQPiScore_NCnumu, h_elMaxQPiScore_CCnue, h_elMaxQPiScore_NCnue, h_elMaxQPiScore_ext = configureHists(h_elMaxQPiScore_CCnumu,
 h_elMaxQPiScore_NCnumu, h_elMaxQPiScore_CCnue, h_elMaxQPiScore_NCnue, h_elMaxQPiScore_ext)

h_elMaxQConf_CCnumu = rt.TH1F("h_elMaxQConf_CCnumu","\"Electron Class Confidence\" for Largest Electron Shower",21,-1,20)
h_elMaxQConf_NCnumu = rt.TH1F("h_elMaxQConf_NCnumu","\"Electron Class Confidence\" for Largest Electron Shower",21,-1,20)
h_elMaxQConf_CCnue = rt.TH1F("h_elMaxQConf_CCnue","\"Electron Class Confidence\" for Largest Electron Shower",21,-1,20)
h_elMaxQConf_NCnue = rt.TH1F("h_elMaxQConf_NCnue","\"Electron Class Confidence\" for Largest Electron Shower",21,-1,20)
h_elMaxQConf_ext = rt.TH1F("h_elMaxQConf_ext","\"Electron Class Confidence\" for Largest Electron Shower",21,-1,20)
h_elMaxQConf_CCnumu, h_elMaxQConf_NCnumu, h_elMaxQConf_CCnue, h_elMaxQConf_NCnue, h_elMaxQConf_ext = configureHists(h_elMaxQConf_CCnumu,
 h_elMaxQConf_NCnumu, h_elMaxQConf_CCnue, h_elMaxQConf_NCnue, h_elMaxQConf_ext)

h_elMaxQ_CCnumu = rt.TH1F("h_elMaxQ_CCnumu","Max Electron Charge",34,0,450000)
h_elMaxQ_NCnumu = rt.TH1F("h_elMaxQ_NCnumu","Max Electron Charge",34,0,450000)
h_elMaxQ_CCnue = rt.TH1F("h_elMaxQ_CCnue","Max Electron Charge",34,0,450000)
h_elMaxQ_NCnue = rt.TH1F("h_elMaxQ_NCnue","Max Electron Charge",34,0,450000)
h_elMaxQ_ext = rt.TH1F("h_elMaxQ_ext","Max Electron Charge",34,0,450000)
h_elMaxQ_CCnumu, h_elMaxQ_NCnumu, h_elMaxQ_CCnue, h_elMaxQ_NCnue, h_elMaxQ_ext = configureHists(h_elMaxQ_CCnumu,
 h_elMaxQ_NCnumu, h_elMaxQ_CCnue, h_elMaxQ_NCnue, h_elMaxQ_ext)

h_muMaxQ_CCnumu = rt.TH1F("h_muMaxQ_CCnumu","Max Muon Charge",34,0,450000)
h_muMaxQ_NCnumu = rt.TH1F("h_muMaxQ_NCnumu","Max Muon Charge",34,0,450000)
h_muMaxQ_CCnue = rt.TH1F("h_muMaxQ_CCnue","Max Muon Charge",34,0,450000)
h_muMaxQ_NCnue = rt.TH1F("h_muMaxQ_NCnue","Max Muon Charge",34,0,450000)
h_muMaxQ_ext = rt.TH1F("h_muMaxQ_ext","Max Muon Charge",34,0,450000)
h_muMaxQ_CCnumu, h_muMaxQ_NCnumu, h_muMaxQ_CCnue, h_muMaxQ_NCnue, h_muMaxQ_ext = configureHists(h_muMaxQ_CCnumu,
 h_muMaxQ_NCnumu, h_muMaxQ_CCnue, h_muMaxQ_NCnue, h_muMaxQ_ext)

h_elMaxQComp_wConfCut_CCnumu = rt.TH1F("h_elMaxQComp_wConfCut_CCnumu","Completeness for Largest Electron Shower",21,0,1.05)
h_elMaxQComp_wConfCut_NCnumu = rt.TH1F("h_elMaxQComp_wConfCut_NCnumu","Completeness for Largest Electron Shower",21,0,1.05)
h_elMaxQComp_wConfCut_CCnue = rt.TH1F("h_elMaxQComp_wConfCut_CCnue","Completeness for Largest Electron Shower",21,0,1.05)
h_elMaxQComp_wConfCut_NCnue = rt.TH1F("h_elMaxQComp_wConfCut_NCnue","Completeness for Largest Electron Shower",21,0,1.05)
h_elMaxQComp_wConfCut_ext = rt.TH1F("h_elMaxQComp_wConfCut_ext","Completeness for Largest Electron Shower",21,0,1.05)
h_elMaxQComp_wConfCut_CCnumu, h_elMaxQComp_wConfCut_NCnumu, h_elMaxQComp_wConfCut_CCnue, h_elMaxQComp_wConfCut_NCnue, h_elMaxQComp_wConfCut_ext = configureHists(h_elMaxQComp_wConfCut_CCnumu,
 h_elMaxQComp_wConfCut_NCnumu, h_elMaxQComp_wConfCut_CCnue, h_elMaxQComp_wConfCut_NCnue, h_elMaxQComp_wConfCut_ext)

h_elMaxQCosTheta_wConfCut_CCnumu = rt.TH1F("h_elMaxQCosTheta_wConfCut_CCnumu","cos(Angle To Beam) for Largest Electron Shower",42,-1.05,1.05)
h_elMaxQCosTheta_wConfCut_NCnumu = rt.TH1F("h_elMaxQCosTheta_wConfCut_NCnumu","cos(Angle To Beam) for Largest Electron Shower",42,-1.05,1.05)
h_elMaxQCosTheta_wConfCut_CCnue = rt.TH1F("h_elMaxQCosTheta_wConfCut_CCnue","cos(Angle To Beam) for Largest Electron Shower",42,-1.05,1.05)
h_elMaxQCosTheta_wConfCut_NCnue = rt.TH1F("h_elMaxQCosTheta_wConfCut_NCnue","cos(Angle To Beam) for Largest Electron Shower",42,-1.05,1.05)
h_elMaxQCosTheta_wConfCut_ext = rt.TH1F("h_elMaxQCosTheta_wConfCut_ext","cos(Angle To Beam) for Largest Electron Shower",42,-1.05,1.05)
h_elMaxQCosTheta_wConfCut_CCnumu, h_elMaxQCosTheta_wConfCut_NCnumu, h_elMaxQCosTheta_wConfCut_CCnue, h_elMaxQCosTheta_wConfCut_NCnue, h_elMaxQCosTheta_wConfCut_ext = configureHists(h_elMaxQCosTheta_wConfCut_CCnumu,
 h_elMaxQCosTheta_wConfCut_NCnumu, h_elMaxQCosTheta_wConfCut_CCnue, h_elMaxQCosTheta_wConfCut_NCnue, h_elMaxQCosTheta_wConfCut_ext)

h_elMaxQ_wConfCut_CCnumu = rt.TH1F("h_elMaxQ_wConfCut_CCnumu","Max Electron Charge",34,0,450000)
h_elMaxQ_wConfCut_NCnumu = rt.TH1F("h_elMaxQ_wConfCut_NCnumu","Max Electron Charge",34,0,450000)
h_elMaxQ_wConfCut_CCnue = rt.TH1F("h_elMaxQ_wConfCut_CCnue","Max Electron Charge",34,0,450000)
h_elMaxQ_wConfCut_NCnue = rt.TH1F("h_elMaxQ_wConfCut_NCnue","Max Electron Charge",34,0,450000)
h_elMaxQ_wConfCut_ext = rt.TH1F("h_elMaxQ_wConfCut_ext","Max Electron Charge",34,0,450000)
h_elMaxQ_wConfCut_CCnumu, h_elMaxQ_wConfCut_NCnumu, h_elMaxQ_wConfCut_CCnue, h_elMaxQ_wConfCut_NCnue, h_elMaxQ_wConfCut_ext = configureHists(h_elMaxQ_wConfCut_CCnumu,
 h_elMaxQ_wConfCut_NCnumu, h_elMaxQ_wConfCut_CCnue, h_elMaxQ_wConfCut_NCnue, h_elMaxQ_wConfCut_ext)

h_elMaxQF_wConfCut_CCnumu = rt.TH1F("h_elMaxQF_wConfCut_CCnumu","Max Electron Charge Fraction",21,0,1.05)
h_elMaxQF_wConfCut_NCnumu = rt.TH1F("h_elMaxQF_wConfCut_NCnumu","Max Electron Charge Fraction",21,0,1.05)
h_elMaxQF_wConfCut_CCnue = rt.TH1F("h_elMaxQF_wConfCut_CCnue","Max Electron Charge Fraction",21,0,1.05)
h_elMaxQF_wConfCut_NCnue = rt.TH1F("h_elMaxQF_wConfCut_NCnue","Max Electron Charge Fraction",21,0,1.05)
h_elMaxQF_wConfCut_ext = rt.TH1F("h_elMaxQF_wConfCut_ext","Max Electron Charge Fraction",21,0,1.05)
h_elMaxQF_wConfCut_CCnumu, h_elMaxQF_wConfCut_NCnumu, h_elMaxQF_wConfCut_CCnue, h_elMaxQF_wConfCut_NCnue, h_elMaxQF_wConfCut_ext = configureHists(h_elMaxQF_wConfCut_CCnumu,
 h_elMaxQF_wConfCut_NCnumu, h_elMaxQF_wConfCut_CCnue, h_elMaxQF_wConfCut_NCnue, h_elMaxQF_wConfCut_ext)

h_elMaxComp_wQFcut_CCnumu = rt.TH1F("h_elMaxComp_wQFcut_CCnumu","Max Electron Completeness",21,0,1.05)
h_elMaxComp_wQFcut_NCnumu = rt.TH1F("h_elMaxComp_wQFcut_NCnumu","Max Electron Completeness",21,0,1.05)
h_elMaxComp_wQFcut_CCnue = rt.TH1F("h_elMaxComp_wQFcut_CCnue","Max Electron Completeness",21,0,1.05)
h_elMaxComp_wQFcut_NCnue = rt.TH1F("h_elMaxComp_wQFcut_NCnue","Max Electron Completeness",21,0,1.05)
h_elMaxComp_wQFcut_ext = rt.TH1F("h_elMaxComp_wQFcut_ext","Max Electron Completeness",21,0,1.05)
h_elMaxComp_wQFcut_CCnumu, h_elMaxComp_wQFcut_NCnumu, h_elMaxComp_wQFcut_CCnue, h_elMaxComp_wQFcut_NCnue, h_elMaxComp_wQFcut_ext = configureHists(h_elMaxComp_wQFcut_CCnumu,
 h_elMaxComp_wQFcut_NCnumu, h_elMaxComp_wQFcut_CCnue, h_elMaxComp_wQFcut_NCnue, h_elMaxComp_wQFcut_ext)

h_elMaxQ_wQFcut_CCnumu = rt.TH1F("h_elMaxQ_wQFcut_CCnumu","Max Electron Charge",34,0,450000)
h_elMaxQ_wQFcut_NCnumu = rt.TH1F("h_elMaxQ_wQFcut_NCnumu","Max Electron Charge",34,0,450000)
h_elMaxQ_wQFcut_CCnue = rt.TH1F("h_elMaxQ_wQFcut_CCnue","Max Electron Charge",34,0,450000)
h_elMaxQ_wQFcut_NCnue = rt.TH1F("h_elMaxQ_wQFcut_NCnue","Max Electron Charge",34,0,450000)
h_elMaxQ_wQFcut_ext = rt.TH1F("h_elMaxQ_wQFcut_ext","Max Electron Charge",34,0,450000)
h_elMaxQ_wQFcut_CCnumu, h_elMaxQ_wQFcut_NCnumu, h_elMaxQ_wQFcut_CCnue, h_elMaxQ_wQFcut_NCnue, h_elMaxQ_wQFcut_ext = configureHists(h_elMaxQ_wQFcut_CCnumu,
 h_elMaxQ_wQFcut_NCnumu, h_elMaxQ_wQFcut_CCnue, h_elMaxQ_wQFcut_NCnue, h_elMaxQ_wQFcut_ext)

h_elMaxComp_wQcut_CCnumu = rt.TH1F("h_elMaxComp_wQcut_CCnumu","Max Electron Completeness",21,0,1.05)
h_elMaxComp_wQcut_NCnumu = rt.TH1F("h_elMaxComp_wQcut_NCnumu","Max Electron Completeness",21,0,1.05)
h_elMaxComp_wQcut_CCnue = rt.TH1F("h_elMaxComp_wQcut_CCnue","Max Electron Completeness",21,0,1.05)
h_elMaxComp_wQcut_NCnue = rt.TH1F("h_elMaxComp_wQcut_NCnue","Max Electron Completeness",21,0,1.05)
h_elMaxComp_wQcut_ext = rt.TH1F("h_elMaxComp_wQcut_ext","Max Electron Completeness",21,0,1.05)
h_elMaxComp_wQcut_CCnumu, h_elMaxComp_wQcut_NCnumu, h_elMaxComp_wQcut_CCnue, h_elMaxComp_wQcut_NCnue, h_elMaxComp_wQcut_ext = configureHists(h_elMaxComp_wQcut_CCnumu,
 h_elMaxComp_wQcut_NCnumu, h_elMaxComp_wQcut_CCnue, h_elMaxComp_wQcut_NCnue, h_elMaxComp_wQcut_ext)

h_elMaxQF_wQcut_CCnumu = rt.TH1F("h_elMaxQF_wQcut_CCnumu","Max Electron Charge Fraction",21,0,1.05)
h_elMaxQF_wQcut_NCnumu = rt.TH1F("h_elMaxQF_wQcut_NCnumu","Max Electron Charge Fraction",21,0,1.05)
h_elMaxQF_wQcut_CCnue = rt.TH1F("h_elMaxQF_wQcut_CCnue","Max Electron Charge Fraction",21,0,1.05)
h_elMaxQF_wQcut_NCnue = rt.TH1F("h_elMaxQF_wQcut_NCnue","Max Electron Charge Fraction",21,0,1.05)
h_elMaxQF_wQcut_ext = rt.TH1F("h_elMaxQF_wQcut_ext","Max Electron Charge Fraction",21,0,1.05)
h_elMaxQF_wQcut_CCnumu, h_elMaxQF_wQcut_NCnumu, h_elMaxQF_wQcut_CCnue, h_elMaxQF_wQcut_NCnue, h_elMaxQF_wQcut_ext = configureHists(h_elMaxQF_wQcut_CCnumu,
 h_elMaxQF_wQcut_NCnumu, h_elMaxQF_wQcut_CCnue, h_elMaxQF_wQcut_NCnue, h_elMaxQF_wQcut_ext)


n_raw_CCnumu = 0
n_raw_NCnumu = 0
n_raw_CCnue = 0
n_raw_NCnue = 0
n_raw_ext = 0

n_runs1to3_CCnumu = 0.
n_runs1to3_CCnue = 0.
n_runs1to3_CCnumu_pass = 0.
n_runs1to3_NCnumu_pass = 0.
n_runs1to3_CCnue_pass = 0.
n_runs1to3_NCnue_pass = 0.
n_runs1to3_ext_pass = 0.


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
      n_raw_CCnumu += 1
      n_runs1to3_CCnumu += tnu.xsecWeight
    #NC numu
    else:
      eventType = 1
      n_raw_NCnumu += 1
  #NC nue
  if abs(tnu.trueNuPDG) == 12 and tnu.trueNuCCNC == 1:
    eventType = 2
    n_raw_NCnue += 1

  if eventType < 0:
    continue

  if tnu.nVertices < 1 or tnu.vtxIsFiducial != 1:
    continue

  h_cosFrac_CCnumu, h_cosFrac_NCnumu, h_cosFrac_NCnue = FillNuHistos(h_cosFrac_CCnumu,
    h_cosFrac_NCnumu, h_cosFrac_NCnue, tnu.vtxFracHitsOnCosmic, tnu.xsecWeight, eventType)

  if tnu.vtxFracHitsOnCosmic >= 1.:
    continue

  nMuons = 0
  nElectrons = 0
  nCompMuons = 0
  nCompElectrons = 0
  elMaxComp = -1.
  elMaxQComp = -1.
  elMaxQCosTheta = -1.
  elMaxQElScore = -1.
  elMaxQPhScore = -1.
  elMaxQPiScore = -1.
  elMaxQConf = -1.
  elMaxQFrac = -1.
  elMaxQ = -1.
  muMaxQ = -1.
  muMaxComp = -1.
  muMaxQComp = -1.
  muMaxQCosTheta = -1.
  muMaxQFrac = -1.

  for iT in range(tnu.nTracks):
    if tnu.trackIsSecondary[iT] == 1:
      continue
    if tnu.trackClassified[iT] == 1 and tnu.trackPID[iT] == 13:
      nMuons += 1
      if tnu.trackComp[iT] > args.compCut:
        nCompMuons += 1
      if tnu.trackCharge[iT] > muMaxQ:
        muMaxQ = tnu.trackCharge[iT]
        muMaxQComp = tnu.trackComp[iT]
        muMaxQCosTheta = tnu.trackCosTheta[iT]
      if tnu.trackChargeFrac[iT] > muMaxQFrac:
        muMaxQFrac = tnu.trackChargeFrac[iT]
      if tnu.trackComp[iT] > muMaxComp:
        muMaxComp = tnu.trackComp[iT]
  for iS in range(tnu.nShowers):
    if tnu.showerIsSecondary[iS] == 1:
      continue
    if tnu.showerClassified[iS] == 1 and tnu.showerPID[iS] == 11:
      nElectrons += 1
      if tnu.showerComp[iS] > args.compCut:
        nCompElectrons += 1
      if tnu.showerCharge[iS] > elMaxQ:
        elMaxQ = tnu.showerCharge[iS]
        elMaxQComp = tnu.showerComp[iS]
        elMaxQCosTheta = tnu.showerCosTheta[iS]
        elMaxQElScore = tnu.showerElScore[iS]
        elMaxQPhScore = tnu.showerPhScore[iS]
        elMaxQPiScore = tnu.showerPiScore[iS]
        elMaxQConf = tnu.showerElScore[iS] - (tnu.showerPhScore[iS] + tnu.showerPiScore[iS])/2.
      if tnu.showerChargeFrac[iS] > elMaxQFrac:
        elMaxQFrac = tnu.showerChargeFrac[iS]
      if tnu.showerComp[iS] > elMaxComp:
        elMaxComp = tnu.showerComp[iS]

  h_nEl_CCnumu, h_nEl_NCnumu, h_nEl_NCnue = FillNuHistos(h_nEl_CCnumu,
    h_nEl_NCnumu, h_nEl_NCnue, nElectrons, tnu.xsecWeight, eventType)
  h_nMu_CCnumu, h_nMu_NCnumu, h_nMu_NCnue = FillNuHistos(h_nMu_CCnumu,
    h_nMu_NCnumu, h_nMu_NCnue, nMuons, tnu.xsecWeight, eventType)
  h_nCompEl_CCnumu, h_nCompEl_NCnumu, h_nCompEl_NCnue = FillNuHistos(h_nCompEl_CCnumu,
    h_nCompEl_NCnumu, h_nCompEl_NCnue, nCompElectrons, tnu.xsecWeight, eventType)
  h_nCompMu_CCnumu, h_nCompMu_NCnumu, h_nCompMu_NCnue = FillNuHistos(h_nCompMu_CCnumu,
    h_nCompMu_NCnumu, h_nCompMu_NCnue, nCompMuons, tnu.xsecWeight, eventType)

  if nElectrons >= 1:
    h_elMaxComp_CCnumu, h_elMaxComp_NCnumu, h_elMaxComp_NCnue = FillNuHistos(h_elMaxComp_CCnumu,
      h_elMaxComp_NCnumu, h_elMaxComp_NCnue, elMaxComp, tnu.xsecWeight, eventType)
    h_elMaxQComp_CCnumu, h_elMaxQComp_NCnumu, h_elMaxQComp_NCnue = FillNuHistos(h_elMaxQComp_CCnumu,
      h_elMaxQComp_NCnumu, h_elMaxQComp_NCnue, elMaxQComp, tnu.xsecWeight, eventType)
    h_elMaxQCosTheta_CCnumu, h_elMaxQCosTheta_NCnumu, h_elMaxQCosTheta_NCnue = FillNuHistos(h_elMaxQCosTheta_CCnumu,
      h_elMaxQCosTheta_NCnumu, h_elMaxQCosTheta_NCnue, elMaxQCosTheta, tnu.xsecWeight, eventType)
    h_elMaxQElScore_CCnumu, h_elMaxQElScore_NCnumu, h_elMaxQElScore_NCnue = FillNuHistos(h_elMaxQElScore_CCnumu,
      h_elMaxQElScore_NCnumu, h_elMaxQElScore_NCnue, elMaxQElScore, tnu.xsecWeight, eventType)
    h_elMaxQPhScore_CCnumu, h_elMaxQPhScore_NCnumu, h_elMaxQPhScore_NCnue = FillNuHistos(h_elMaxQPhScore_CCnumu,
      h_elMaxQPhScore_NCnumu, h_elMaxQPhScore_NCnue, elMaxQPhScore, tnu.xsecWeight, eventType)
    h_elMaxQPiScore_CCnumu, h_elMaxQPiScore_NCnumu, h_elMaxQPiScore_NCnue = FillNuHistos(h_elMaxQPiScore_CCnumu,
      h_elMaxQPiScore_NCnumu, h_elMaxQPiScore_NCnue, elMaxQPiScore, tnu.xsecWeight, eventType)
    h_elMaxQConf_CCnumu, h_elMaxQConf_NCnumu, h_elMaxQConf_NCnue = FillNuHistos(h_elMaxQConf_CCnumu,
      h_elMaxQConf_NCnumu, h_elMaxQConf_NCnue, elMaxQConf, tnu.xsecWeight, eventType)
    h_elMaxQF_CCnumu, h_elMaxQF_NCnumu, h_elMaxQF_NCnue = FillNuHistos(h_elMaxQF_CCnumu,
      h_elMaxQF_NCnumu, h_elMaxQF_NCnue, elMaxQFrac, tnu.xsecWeight, eventType)
    h_elMaxQ_CCnumu, h_elMaxQ_NCnumu, h_elMaxQ_NCnue = FillNuHistos(h_elMaxQ_CCnumu,
      h_elMaxQ_NCnumu, h_elMaxQ_NCnue, elMaxQ, tnu.xsecWeight, eventType)

    if elMaxQFrac > args.chargeFracCut:
      h_elMaxComp_wQFcut_CCnumu, h_elMaxComp_wQFcut_NCnumu, h_elMaxComp_wQFcut_NCnue = FillNuHistos(h_elMaxComp_wQFcut_CCnumu,
        h_elMaxComp_wQFcut_NCnumu, h_elMaxComp_wQFcut_NCnue, elMaxComp, tnu.xsecWeight, eventType)
      h_elMaxQ_wQFcut_CCnumu, h_elMaxQ_wQFcut_NCnumu, h_elMaxQ_wQFcut_NCnue = FillNuHistos(h_elMaxQ_wQFcut_CCnumu,
        h_elMaxQ_wQFcut_NCnumu, h_elMaxQ_wQFcut_NCnue, elMaxQ, tnu.xsecWeight, eventType)

    if elMaxQ > args.chargeCut:
      h_elMaxComp_wQcut_CCnumu, h_elMaxComp_wQcut_NCnumu, h_elMaxComp_wQcut_NCnue = FillNuHistos(h_elMaxComp_wQcut_CCnumu,
        h_elMaxComp_wQcut_NCnumu, h_elMaxComp_wQcut_NCnue, elMaxComp, tnu.xsecWeight, eventType)
      h_elMaxQF_wQcut_CCnumu, h_elMaxQF_wQcut_NCnumu, h_elMaxQF_wQcut_NCnue = FillNuHistos(h_elMaxQF_wQcut_CCnumu,
        h_elMaxQF_wQcut_NCnumu, h_elMaxQF_wQcut_NCnue, elMaxQFrac, tnu.xsecWeight, eventType)

    if elMaxQConf > args.confCut:
      h_elMaxQComp_wConfCut_CCnumu, h_elMaxQComp_wConfCut_NCnumu, h_elMaxQComp_wConfCut_NCnue = FillNuHistos(h_elMaxQComp_wConfCut_CCnumu,
        h_elMaxQComp_wConfCut_NCnumu, h_elMaxQComp_wConfCut_NCnue, elMaxQComp, tnu.xsecWeight, eventType)
      h_elMaxQCosTheta_wConfCut_CCnumu, h_elMaxQCosTheta_wConfCut_NCnumu, h_elMaxQCosTheta_wConfCut_NCnue = FillNuHistos(h_elMaxQCosTheta_wConfCut_CCnumu,
        h_elMaxQCosTheta_wConfCut_NCnumu, h_elMaxQCosTheta_wConfCut_NCnue, elMaxQCosTheta, tnu.xsecWeight, eventType)
      h_elMaxQF_wConfCut_CCnumu, h_elMaxQF_wConfCut_NCnumu, h_elMaxQF_wConfCut_NCnue = FillNuHistos(h_elMaxQF_wConfCut_CCnumu,
        h_elMaxQF_wConfCut_NCnumu, h_elMaxQF_wConfCut_NCnue, elMaxQFrac, tnu.xsecWeight, eventType)
      h_elMaxQ_wConfCut_CCnumu, h_elMaxQ_wConfCut_NCnumu, h_elMaxQ_wConfCut_NCnue = FillNuHistos(h_elMaxQ_wConfCut_CCnumu,
        h_elMaxQ_wConfCut_NCnumu, h_elMaxQ_wConfCut_NCnue, elMaxQ, tnu.xsecWeight, eventType)
      if elMaxQ > args.chargeCut and elMaxQFrac > args.chargeFracCut:
        if eventType == 0:
          n_runs1to3_CCnumu_pass += tnu.xsecWeight
        if eventType == 1:
          n_runs1to3_NCnumu_pass += tnu.xsecWeight
        if eventType == 2:
          n_runs1to3_NCnue_pass += tnu.xsecWeight

  if nMuons >= 1:
    h_muMaxComp_CCnumu, h_muMaxComp_NCnumu, h_muMaxComp_NCnue = FillNuHistos(h_muMaxComp_CCnumu,
      h_muMaxComp_NCnumu, h_muMaxComp_NCnue, muMaxComp, tnu.xsecWeight, eventType)
    h_muMaxQComp_CCnumu, h_muMaxQComp_NCnumu, h_muMaxQComp_NCnue = FillNuHistos(h_muMaxQComp_CCnumu,
      h_muMaxQComp_NCnumu, h_muMaxQComp_NCnue, muMaxQComp, tnu.xsecWeight, eventType)
    h_muMaxQCosTheta_CCnumu, h_muMaxQCosTheta_NCnumu, h_muMaxQCosTheta_NCnue = FillNuHistos(h_muMaxQCosTheta_CCnumu,
      h_muMaxQCosTheta_NCnumu, h_muMaxQCosTheta_NCnue, muMaxQCosTheta, tnu.xsecWeight, eventType)
    h_muMaxQF_CCnumu, h_muMaxQF_NCnumu, h_muMaxQF_NCnue = FillNuHistos(h_muMaxQF_CCnumu,
      h_muMaxQF_NCnumu, h_muMaxQF_NCnue, muMaxQFrac, tnu.xsecWeight, eventType)
    h_muMaxQ_CCnumu, h_muMaxQ_NCnumu, h_muMaxQ_NCnue = FillNuHistos(h_muMaxQ_CCnumu,
      h_muMaxQ_NCnumu, h_muMaxQ_NCnue, muMaxQ, tnu.xsecWeight, eventType)
    



for i in range(tnue.GetEntries()):

  tnue.GetEntry(i)

  if abs(tnue.trueNuPDG != 12) or tnue.trueNuCCNC != 0 or isinf(tnue.xsecWeight):
    continue

  n_raw_CCnue += 1
  n_runs1to3_CCnue += tnue.xsecWeight

  if tnu.nVertices < 1 or tnu.vtxIsFiducial != 1:
    continue

  h_cosFrac_CCnue.Fill(tnue.vtxFracHitsOnCosmic, tnue.xsecWeight)
  #if tnue.vtxFracHitsOnCosmic < 0 or tnue.vtxFracHitsOnCosmic > 1.:
  if tnue.vtxFracHitsOnCosmic > 1.:
    print(tnue.vtxFracHitsOnCosmic)

  if tnue.vtxFracHitsOnCosmic >= 1.:
    continue

  nMuons = 0
  nElectrons = 0
  nCompMuons = 0
  nCompElectrons = 0
  elMaxComp = -1.
  elMaxQComp = -1.
  elMaxQCosTheta = -1.
  elMaxQElScore = -1.
  elMaxQPhScore = -1.
  elMaxQPiScore = -1.
  elMaxQConf = -1.
  elMaxQFrac = -1.
  elMaxQ = -1.
  muMaxQ = -1.
  muMaxComp = -1.
  muMaxQComp = -1.
  muMaxQCosTheta = -1.
  muMaxQFrac = -1.

  for iT in range(tnue.nTracks):
    if tnue.trackIsSecondary[iT] == 1:
      continue
    if tnue.trackClassified[iT] == 1 and tnue.trackPID[iT] == 13:
      nMuons += 1
      if tnue.trackComp[iT] > args.compCut:
        nCompMuons += 1
      if tnue.trackCharge[iT] > muMaxQ:
        muMaxQ = tnue.trackCharge[iT]
        muMaxQComp = tnue.trackComp[iT]
        muMaxQCosTheta = tnue.trackCosTheta[iT]
      if tnue.trackChargeFrac[iT] > muMaxQFrac:
        muMaxQFrac = tnue.trackChargeFrac[iT]
      if tnue.trackComp[iT] > muMaxComp:
        muMaxComp = tnue.trackComp[iT]
  for iS in range(tnue.nShowers):
    if tnue.showerIsSecondary[iS] == 1:
      continue
    if tnue.showerClassified[iS] == 1 and tnue.showerPID[iS] == 11:
      nElectrons += 1
      if tnue.showerComp[iS] > args.compCut:
        nCompElectrons += 1
      if tnue.showerCharge[iS] > elMaxQ:
        elMaxQ = tnue.showerCharge[iS]
        elMaxQComp = tnue.showerComp[iS]
        elMaxQCosTheta = tnue.showerCosTheta[iS]
        elMaxQElScore = tnue.showerElScore[iS]
        elMaxQPhScore = tnue.showerPhScore[iS]
        elMaxQPiScore = tnue.showerPiScore[iS]
        elMaxQConf = tnue.showerElScore[iS] - (tnue.showerPhScore[iS] + tnue.showerPiScore[iS])/2.
      if tnue.showerChargeFrac[iS] > elMaxQFrac:
        elMaxQFrac = tnue.showerChargeFrac[iS]
      if tnue.showerComp[iS] > elMaxComp:
        elMaxComp = tnue.showerComp[iS]

  h_nEl_CCnue.Fill(nElectrons, tnue.xsecWeight)
  h_nMu_CCnue.Fill(nMuons, tnue.xsecWeight)
  h_nCompEl_CCnue.Fill(nCompElectrons, tnue.xsecWeight)
  h_nCompMu_CCnue.Fill(nCompMuons, tnue.xsecWeight)

  if nElectrons >= 1:
    h_elMaxComp_CCnue.Fill(elMaxComp, tnue.xsecWeight)
    h_elMaxQComp_CCnue.Fill(elMaxQComp, tnue.xsecWeight)
    h_elMaxQCosTheta_CCnue.Fill(elMaxQCosTheta, tnue.xsecWeight)
    h_elMaxQElScore_CCnue.Fill(elMaxQElScore, tnue.xsecWeight)
    h_elMaxQPhScore_CCnue.Fill(elMaxQPhScore, tnue.xsecWeight)
    h_elMaxQPiScore_CCnue.Fill(elMaxQPiScore, tnue.xsecWeight)
    h_elMaxQConf_CCnue.Fill(elMaxQConf, tnue.xsecWeight)
    h_elMaxQF_CCnue.Fill(elMaxQFrac, tnue.xsecWeight)
    h_elMaxQ_CCnue.Fill(elMaxQ, tnue.xsecWeight)
    if elMaxQFrac > args.chargeFracCut:
      h_elMaxComp_wQFcut_CCnue.Fill(elMaxComp, tnue.xsecWeight)
      h_elMaxQ_wQFcut_CCnue.Fill(elMaxQ, tnue.xsecWeight)
    if elMaxQ > args.chargeCut:
      h_elMaxComp_wQcut_CCnue.Fill(elMaxComp, tnue.xsecWeight)
      h_elMaxQF_wQcut_CCnue.Fill(elMaxQFrac, tnue.xsecWeight)
    if elMaxQConf > args.confCut:
      h_elMaxQComp_wConfCut_CCnue.Fill(elMaxQComp, tnue.xsecWeight)
      h_elMaxQCosTheta_wConfCut_CCnue.Fill(elMaxQCosTheta, tnue.xsecWeight)
      h_elMaxQF_wConfCut_CCnue.Fill(elMaxQFrac, tnue.xsecWeight)
      h_elMaxQ_wConfCut_CCnue.Fill(elMaxQ, tnue.xsecWeight)
      if elMaxQ > args.chargeCut and elMaxQFrac > args.chargeFracCut:
        n_runs1to3_CCnue_pass += tnue.xsecWeight

  if nMuons >= 1:
    h_muMaxComp_CCnue.Fill(muMaxComp, tnue.xsecWeight)
    h_muMaxQComp_CCnue.Fill(muMaxQComp, tnue.xsecWeight)
    h_muMaxQCosTheta_CCnue.Fill(muMaxQCosTheta, tnue.xsecWeight)
    h_muMaxQF_CCnue.Fill(muMaxQFrac, tnue.xsecWeight)
    h_muMaxQ_CCnue.Fill(muMaxQ, tnue.xsecWeight)





for i in range(text.GetEntries()):

  text.GetEntry(i)

  n_raw_ext += 1

  if text.nVertices < 1 or text.vtxIsFiducial != 1:
    continue

  h_cosFrac_ext.Fill(text.vtxFracHitsOnCosmic)

  if text.vtxFracHitsOnCosmic >= 1.:
    continue

  nMuons = 0
  nElectrons = 0
  nCompMuons = 0
  nCompElectrons = 0
  elMaxComp = -1.
  elMaxQComp = -1.
  elMaxQCosTheta = -1.
  elMaxQElScore = -1.
  elMaxQPhScore = -1.
  elMaxQPiScore = -1.
  elMaxQConf = -1.
  elMaxQFrac = -1.
  elMaxQ = -1.
  muMaxQ = -1.
  muMaxComp = -1.
  muMaxQComp = -1.
  muMaxQCosTheta = -1.
  muMaxQFrac = -1.
  for iT in range(text.nTracks):
    if text.trackIsSecondary[iT] == 1:
      continue
    if text.trackClassified[iT] == 1 and text.trackPID[iT] == 13:
      nMuons += 1
      if text.trackComp[iT] > args.compCut:
        nCompMuons += 1
      if text.trackCharge[iT] > muMaxQ:
        muMaxQ = text.trackCharge[iT]
        muMaxQComp = text.trackComp[iT]
        muMaxQCosTheta = text.trackCosTheta[iT]
      if text.trackChargeFrac[iT] > muMaxQFrac:
        muMaxQFrac = text.trackChargeFrac[iT]
      if text.trackComp[iT] > muMaxComp:
        muMaxComp = text.trackComp[iT]
  for iS in range(text.nShowers):
    if text.showerIsSecondary[iS] == 1:
      continue
    if text.showerClassified[iS] == 1 and text.showerPID[iS] == 11:
      nElectrons += 1
      if text.showerComp[iS] > args.compCut:
        nCompElectrons += 1
      if text.showerCharge[iS] > elMaxQ:
        elMaxQ = text.showerCharge[iS]
        elMaxQComp = text.showerComp[iS]
        elMaxQCosTheta = text.showerCosTheta[iS]
        elMaxQElScore = text.showerElScore[iS]
        elMaxQPhScore = text.showerPhScore[iS]
        elMaxQPiScore = text.showerPiScore[iS]
        elMaxQConf = text.showerElScore[iS] - (text.showerPhScore[iS] + text.showerPiScore[iS])/2.
      if text.showerChargeFrac[iS] > elMaxQFrac:
        elMaxQFrac = text.showerChargeFrac[iS]
      if text.showerComp[iS] > elMaxComp:
        elMaxComp = text.showerComp[iS]

  h_nEl_ext.Fill(nElectrons)
  h_nMu_ext.Fill(nMuons)
  h_nCompEl_ext.Fill(nCompElectrons)
  h_nCompMu_ext.Fill(nCompMuons)

  if nElectrons >= 1:
    h_elMaxComp_ext.Fill(elMaxComp)
    h_elMaxQComp_ext.Fill(elMaxQComp)
    h_elMaxQCosTheta_ext.Fill(elMaxQCosTheta)
    h_elMaxQElScore_ext.Fill(elMaxQElScore)
    h_elMaxQPhScore_ext.Fill(elMaxQPhScore)
    h_elMaxQPiScore_ext.Fill(elMaxQPiScore)
    h_elMaxQConf_ext.Fill(elMaxQConf)
    h_elMaxQF_ext.Fill(elMaxQFrac)
    h_elMaxQ_ext.Fill(elMaxQ)
    if elMaxQFrac > args.chargeFracCut:
      h_elMaxComp_wQFcut_ext.Fill(elMaxComp)
      h_elMaxQ_wQFcut_ext.Fill(elMaxQ)
    if elMaxQ > args.chargeCut:
      h_elMaxComp_wQcut_ext.Fill(elMaxComp)
      h_elMaxQF_wQcut_ext.Fill(elMaxQFrac)
    if elMaxQConf > args.confCut:
      h_elMaxQComp_wConfCut_ext.Fill(elMaxQComp)
      h_elMaxQCosTheta_wConfCut_ext.Fill(elMaxQCosTheta)
      h_elMaxQF_wConfCut_ext.Fill(elMaxQFrac)
      h_elMaxQ_wConfCut_ext.Fill(elMaxQ)
      if elMaxQ > args.chargeCut and elMaxQFrac > args.chargeFracCut:
        n_runs1to3_ext_pass += 1.

  if nMuons >= 1:
    h_muMaxComp_ext.Fill(muMaxComp)
    h_muMaxQComp_ext.Fill(muMaxQComp)
    h_muMaxQCosTheta_ext.Fill(muMaxQCosTheta)
    h_muMaxQF_ext.Fill(muMaxQFrac)
    h_muMaxQ_ext.Fill(muMaxQ)




n_runs1to3_CCnumu *= (runs1to3POT/tnuPOTsum)
n_runs1to3_CCnue *= (runs1to3POT/tnuePOTsum)
n_runs1to3_CCnumu_pass *= runs1to3POT/tnuPOTsum
n_runs1to3_NCnumu_pass *= runs1to3POT/tnuPOTsum
n_runs1to3_CCnue_pass *= runs1to3POT/tnuePOTsum
n_runs1to3_NCnue_pass *= runs1to3POT/tnuPOTsum
n_runs1to3_ext_pass *= runs1to3POT/textPOTsum

print()
print("n_raw_CCnumu: ", n_raw_CCnumu)
print("n_raw_NCnumu: ", n_raw_NCnumu)
print("n_raw_CCnue: ", n_raw_CCnue)
print("n_raw_NCnue: ", n_raw_NCnue)
print("n_raw_ext: ", n_raw_ext)
print()
print("n_runs1to3_CCnumu: ", n_runs1to3_CCnumu)
print("n_runs1to3_CCnue: ", n_runs1to3_CCnue)
print()
print("%s integral: %f ; scaling by %.2e/%.2e"%("h_cosFrac_CCnue",h_cosFrac_CCnue.Integral(),runs1to3POT,tnuePOTsum))
print("%s integral: %f ; scaling by %.2e/%.2e"%("h_cosFrac_CCnumu",h_cosFrac_CCnumu.Integral(),runs1to3POT,tnuPOTsum))
print("%s integral: %f ; scaling by %.2e/%.2e"%("h_cosFrac_NCnue",h_cosFrac_NCnue.Integral(),runs1to3POT,tnuPOTsum))
print("%s integral: %f ; scaling by %.2e/%.2e"%("h_cosFrac_NCnumu",h_cosFrac_NCnumu.Integral(),runs1to3POT,tnuPOTsum))
print("%s integral: %f ; scaling by %.2e/%.2e"%("h_cosFrac_ext",h_cosFrac_ext.Integral(),runs1to3POT,textPOTsum))
print()
print("CC nue cut efficiency: %f"%(n_runs1to3_CCnue_pass/n_runs1to3_CCnue))
print("CC nue cut purity: %f"%(n_runs1to3_CCnue_pass/(n_runs1to3_CCnue_pass + n_runs1to3_NCnue_pass + n_runs1to3_CCnumu_pass + n_runs1to3_NCnumu_pass + n_runs1to3_ext_pass)))
print()

#print("h_cosFrac_CCnue unscaled integral: ", h_cosFrac_CCnue.Integral())
#print("h_nEl_CCnue unscaled integral: ", h_nEl_CCnue.Integral())

h_cosFrac_CCnumu.Scale(runs1to3POT/tnuPOTsum)
h_cosFrac_NCnumu.Scale(runs1to3POT/tnuPOTsum)
h_cosFrac_CCnue.Scale(runs1to3POT/tnuePOTsum)
h_cosFrac_NCnue.Scale(runs1to3POT/tnuPOTsum)
h_cosFrac_ext.Scale(runs1to3POT/textPOTsum)

h_nEl_CCnumu.Scale(runs1to3POT/tnuPOTsum)
h_nEl_NCnumu.Scale(runs1to3POT/tnuPOTsum)
h_nEl_CCnue.Scale(runs1to3POT/tnuePOTsum)
h_nEl_NCnue.Scale(runs1to3POT/tnuPOTsum)
h_nEl_ext.Scale(runs1to3POT/textPOTsum)

#print("h_cosFrac_CCnue scaled integral: ", h_cosFrac_CCnue.Integral())
#print("h_nEl_CCnue scaled integral: ", h_nEl_CCnue.Integral())

h_nCompEl_CCnumu.Scale(runs1to3POT/tnuPOTsum)
h_nCompEl_NCnumu.Scale(runs1to3POT/tnuPOTsum)
h_nCompEl_CCnue.Scale(runs1to3POT/tnuePOTsum)
h_nCompEl_NCnue.Scale(runs1to3POT/tnuPOTsum)
h_nCompEl_ext.Scale(runs1to3POT/textPOTsum)

h_nMu_CCnumu.Scale(runs1to3POT/tnuPOTsum)
h_nMu_NCnumu.Scale(runs1to3POT/tnuPOTsum)
h_nMu_CCnue.Scale(runs1to3POT/tnuePOTsum)
h_nMu_NCnue.Scale(runs1to3POT/tnuPOTsum)
h_nMu_ext.Scale(runs1to3POT/textPOTsum)

h_nCompMu_CCnumu.Scale(runs1to3POT/tnuPOTsum)
h_nCompMu_NCnumu.Scale(runs1to3POT/tnuPOTsum)
h_nCompMu_CCnue.Scale(runs1to3POT/tnuePOTsum)
h_nCompMu_NCnue.Scale(runs1to3POT/tnuPOTsum)
h_nCompMu_ext.Scale(runs1to3POT/textPOTsum)

h_elMaxQ_CCnumu.Scale(runs1to3POT/tnuPOTsum)
h_elMaxQ_NCnumu.Scale(runs1to3POT/tnuPOTsum)
h_elMaxQ_CCnue.Scale(runs1to3POT/tnuePOTsum)
h_elMaxQ_NCnue.Scale(runs1to3POT/tnuPOTsum)
h_elMaxQ_ext.Scale(runs1to3POT/textPOTsum)

h_muMaxQ_CCnumu.Scale(runs1to3POT/tnuPOTsum)
h_muMaxQ_NCnumu.Scale(runs1to3POT/tnuPOTsum)
h_muMaxQ_CCnue.Scale(runs1to3POT/tnuePOTsum)
h_muMaxQ_NCnue.Scale(runs1to3POT/tnuPOTsum)
h_muMaxQ_ext.Scale(runs1to3POT/textPOTsum)

h_elMaxQF_CCnumu.Scale(runs1to3POT/tnuPOTsum)
h_elMaxQF_NCnumu.Scale(runs1to3POT/tnuPOTsum)
h_elMaxQF_CCnue.Scale(runs1to3POT/tnuePOTsum)
h_elMaxQF_NCnue.Scale(runs1to3POT/tnuPOTsum)
h_elMaxQF_ext.Scale(runs1to3POT/textPOTsum)

h_muMaxQF_CCnumu.Scale(runs1to3POT/tnuPOTsum)
h_muMaxQF_NCnumu.Scale(runs1to3POT/tnuPOTsum)
h_muMaxQF_CCnue.Scale(runs1to3POT/tnuePOTsum)
h_muMaxQF_NCnue.Scale(runs1to3POT/tnuPOTsum)
h_muMaxQF_ext.Scale(runs1to3POT/textPOTsum)

h_elMaxComp_CCnumu.Scale(runs1to3POT/tnuPOTsum)
h_elMaxComp_NCnumu.Scale(runs1to3POT/tnuPOTsum)
h_elMaxComp_CCnue.Scale(runs1to3POT/tnuePOTsum)
h_elMaxComp_NCnue.Scale(runs1to3POT/tnuPOTsum)
h_elMaxComp_ext.Scale(runs1to3POT/textPOTsum)

h_muMaxComp_CCnumu.Scale(runs1to3POT/tnuPOTsum)
h_muMaxComp_NCnumu.Scale(runs1to3POT/tnuPOTsum)
h_muMaxComp_CCnue.Scale(runs1to3POT/tnuePOTsum)
h_muMaxComp_NCnue.Scale(runs1to3POT/tnuPOTsum)
h_muMaxComp_ext.Scale(runs1to3POT/textPOTsum)

h_elMaxQComp_CCnumu.Scale(runs1to3POT/tnuPOTsum)
h_elMaxQComp_NCnumu.Scale(runs1to3POT/tnuPOTsum)
h_elMaxQComp_CCnue.Scale(runs1to3POT/tnuePOTsum)
h_elMaxQComp_NCnue.Scale(runs1to3POT/tnuPOTsum)
h_elMaxQComp_ext.Scale(runs1to3POT/textPOTsum)

h_muMaxQComp_CCnumu.Scale(runs1to3POT/tnuPOTsum)
h_muMaxQComp_NCnumu.Scale(runs1to3POT/tnuPOTsum)
h_muMaxQComp_CCnue.Scale(runs1to3POT/tnuePOTsum)
h_muMaxQComp_NCnue.Scale(runs1to3POT/tnuPOTsum)
h_muMaxQComp_ext.Scale(runs1to3POT/textPOTsum)

h_elMaxQCosTheta_CCnumu.Scale(runs1to3POT/tnuPOTsum)
h_elMaxQCosTheta_NCnumu.Scale(runs1to3POT/tnuPOTsum)
h_elMaxQCosTheta_CCnue.Scale(runs1to3POT/tnuePOTsum)
h_elMaxQCosTheta_NCnue.Scale(runs1to3POT/tnuPOTsum)
h_elMaxQCosTheta_ext.Scale(runs1to3POT/textPOTsum)

h_muMaxQCosTheta_CCnumu.Scale(runs1to3POT/tnuPOTsum)
h_muMaxQCosTheta_NCnumu.Scale(runs1to3POT/tnuPOTsum)
h_muMaxQCosTheta_CCnue.Scale(runs1to3POT/tnuePOTsum)
h_muMaxQCosTheta_NCnue.Scale(runs1to3POT/tnuPOTsum)
h_muMaxQCosTheta_ext.Scale(runs1to3POT/textPOTsum)

h_elMaxQElScore_CCnumu.Scale(runs1to3POT/tnuPOTsum)
h_elMaxQElScore_NCnumu.Scale(runs1to3POT/tnuPOTsum)
h_elMaxQElScore_CCnue.Scale(runs1to3POT/tnuePOTsum)
h_elMaxQElScore_NCnue.Scale(runs1to3POT/tnuPOTsum)
h_elMaxQElScore_ext.Scale(runs1to3POT/textPOTsum)

h_elMaxQPhScore_CCnumu.Scale(runs1to3POT/tnuPOTsum)
h_elMaxQPhScore_NCnumu.Scale(runs1to3POT/tnuPOTsum)
h_elMaxQPhScore_CCnue.Scale(runs1to3POT/tnuePOTsum)
h_elMaxQPhScore_NCnue.Scale(runs1to3POT/tnuPOTsum)
h_elMaxQPhScore_ext.Scale(runs1to3POT/textPOTsum)

h_elMaxQPiScore_CCnumu.Scale(runs1to3POT/tnuPOTsum)
h_elMaxQPiScore_NCnumu.Scale(runs1to3POT/tnuPOTsum)
h_elMaxQPiScore_CCnue.Scale(runs1to3POT/tnuePOTsum)
h_elMaxQPiScore_NCnue.Scale(runs1to3POT/tnuPOTsum)
h_elMaxQPiScore_ext.Scale(runs1to3POT/textPOTsum)

h_elMaxQConf_CCnumu.Scale(runs1to3POT/tnuPOTsum)
h_elMaxQConf_NCnumu.Scale(runs1to3POT/tnuPOTsum)
h_elMaxQConf_CCnue.Scale(runs1to3POT/tnuePOTsum)
h_elMaxQConf_NCnue.Scale(runs1to3POT/tnuPOTsum)
h_elMaxQConf_ext.Scale(runs1to3POT/textPOTsum)

h_elMaxQComp_wConfCut_CCnumu.Scale(runs1to3POT/tnuPOTsum)
h_elMaxQComp_wConfCut_NCnumu.Scale(runs1to3POT/tnuPOTsum)
h_elMaxQComp_wConfCut_CCnue.Scale(runs1to3POT/tnuePOTsum)
h_elMaxQComp_wConfCut_NCnue.Scale(runs1to3POT/tnuPOTsum)
h_elMaxQComp_wConfCut_ext.Scale(runs1to3POT/textPOTsum)

h_elMaxQCosTheta_wConfCut_CCnumu.Scale(runs1to3POT/tnuPOTsum)
h_elMaxQCosTheta_wConfCut_NCnumu.Scale(runs1to3POT/tnuPOTsum)
h_elMaxQCosTheta_wConfCut_CCnue.Scale(runs1to3POT/tnuePOTsum)
h_elMaxQCosTheta_wConfCut_NCnue.Scale(runs1to3POT/tnuPOTsum)
h_elMaxQCosTheta_wConfCut_ext.Scale(runs1to3POT/textPOTsum)

h_elMaxQ_wConfCut_CCnumu.Scale(runs1to3POT/tnuPOTsum)
h_elMaxQ_wConfCut_NCnumu.Scale(runs1to3POT/tnuPOTsum)
h_elMaxQ_wConfCut_CCnue.Scale(runs1to3POT/tnuePOTsum)
h_elMaxQ_wConfCut_NCnue.Scale(runs1to3POT/tnuPOTsum)
h_elMaxQ_wConfCut_ext.Scale(runs1to3POT/textPOTsum)

h_elMaxQF_wConfCut_CCnumu.Scale(runs1to3POT/tnuPOTsum)
h_elMaxQF_wConfCut_NCnumu.Scale(runs1to3POT/tnuPOTsum)
h_elMaxQF_wConfCut_CCnue.Scale(runs1to3POT/tnuePOTsum)
h_elMaxQF_wConfCut_NCnue.Scale(runs1to3POT/tnuPOTsum)
h_elMaxQF_wConfCut_ext.Scale(runs1to3POT/textPOTsum)

h_elMaxQ_wQFcut_CCnumu.Scale(runs1to3POT/tnuPOTsum)
h_elMaxQ_wQFcut_NCnumu.Scale(runs1to3POT/tnuPOTsum)
h_elMaxQ_wQFcut_CCnue.Scale(runs1to3POT/tnuePOTsum)
h_elMaxQ_wQFcut_NCnue.Scale(runs1to3POT/tnuPOTsum)
h_elMaxQ_wQFcut_ext.Scale(runs1to3POT/textPOTsum)

h_elMaxComp_wQFcut_CCnumu.Scale(runs1to3POT/tnuPOTsum)
h_elMaxComp_wQFcut_NCnumu.Scale(runs1to3POT/tnuPOTsum)
h_elMaxComp_wQFcut_CCnue.Scale(runs1to3POT/tnuePOTsum)
h_elMaxComp_wQFcut_NCnue.Scale(runs1to3POT/tnuPOTsum)
h_elMaxComp_wQFcut_ext.Scale(runs1to3POT/textPOTsum)

h_elMaxQF_wQcut_CCnumu.Scale(runs1to3POT/tnuPOTsum)
h_elMaxQF_wQcut_NCnumu.Scale(runs1to3POT/tnuPOTsum)
h_elMaxQF_wQcut_CCnue.Scale(runs1to3POT/tnuePOTsum)
h_elMaxQF_wQcut_NCnue.Scale(runs1to3POT/tnuPOTsum)
h_elMaxQF_wQcut_ext.Scale(runs1to3POT/textPOTsum)

h_elMaxComp_wQcut_CCnumu.Scale(runs1to3POT/tnuPOTsum)
h_elMaxComp_wQcut_NCnumu.Scale(runs1to3POT/tnuPOTsum)
h_elMaxComp_wQcut_CCnue.Scale(runs1to3POT/tnuePOTsum)
h_elMaxComp_wQcut_NCnue.Scale(runs1to3POT/tnuPOTsum)
h_elMaxComp_wQcut_ext.Scale(runs1to3POT/textPOTsum)



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

cnv_nEl = rt.TCanvas("cnv_nEl","cnv_nEl")
hists_nEl = sortHists([h_nEl_CCnumu, h_nEl_NCnumu, h_nEl_CCnue, h_nEl_NCnue, h_nEl_ext])
hists_nEl[0].Draw("EHIST")
for i in range(1,len(hists_nEl)):
  hists_nEl[i].Draw("EHISTSAME")
leg_nEl = rt.TLegend(0.7,0.7,0.9,0.9)
leg_nEl = configureLegend(leg_nEl, h_nEl_CCnumu,
  h_nEl_NCnumu, h_nEl_CCnue, h_nEl_NCnue, h_nEl_ext)
leg_nEl.Draw()
#cnv_nEl.SaveAs("nEl.png")
cnv_nEl.Write()

cnv_nCompEl = rt.TCanvas("cnv_nCompEl","cnv_nCompEl")
hists_nCompEl = sortHists([h_nCompEl_CCnumu, h_nCompEl_NCnumu, h_nCompEl_CCnue, h_nCompEl_NCnue, h_nCompEl_ext])
hists_nCompEl[0].Draw("EHIST")
for i in range(1,len(hists_nCompEl)):
  hists_nCompEl[i].Draw("EHISTSAME")
leg_nCompEl = rt.TLegend(0.7,0.7,0.9,0.9)
leg_nCompEl = configureLegend(leg_nCompEl, h_nCompEl_CCnumu,
  h_nCompEl_NCnumu, h_nCompEl_CCnue, h_nCompEl_NCnue, h_nCompEl_ext)
leg_nCompEl.Draw()
#cnv_nCompEl.SaveAs("nCompEl.png")
cnv_nCompEl.Write()

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

cnv_nCompMu = rt.TCanvas("cnv_nCompMu","cnv_nCompMu")
hists_nCompMu = sortHists([h_nCompMu_CCnumu, h_nCompMu_NCnumu, h_nCompMu_CCnue, h_nCompMu_NCnue, h_nCompMu_ext])
hists_nCompMu[0].Draw("EHIST")
for i in range(1,len(hists_nCompMu)):
  hists_nCompMu[i].Draw("EHISTSAME")
leg_nCompMu = rt.TLegend(0.7,0.7,0.9,0.9)
leg_nCompMu = configureLegend(leg_nCompMu, h_nCompMu_CCnumu,
  h_nCompMu_NCnumu, h_nCompMu_CCnue, h_nCompMu_NCnue, h_nCompMu_ext)
leg_nCompMu.Draw()
#cnv_nCompMu.SaveAs("nCompMu.png")
cnv_nCompMu.Write()

cnv_elMaxComp = rt.TCanvas("cnv_elMaxComp","cnv_elMaxComp")
hists_elMaxComp = sortHists([h_elMaxComp_CCnumu, h_elMaxComp_NCnumu, h_elMaxComp_CCnue, h_elMaxComp_NCnue, h_elMaxComp_ext])
hists_elMaxComp[0].Draw("EHIST")
for i in range(1,len(hists_elMaxComp)):
  hists_elMaxComp[i].Draw("EHISTSAME")
leg_elMaxComp = rt.TLegend(0.7,0.7,0.9,0.9)
leg_elMaxComp = configureLegend(leg_elMaxComp, h_elMaxComp_CCnumu,
  h_elMaxComp_NCnumu, h_elMaxComp_CCnue, h_elMaxComp_NCnue, h_elMaxComp_ext)
leg_elMaxComp.Draw()
#cnv_elMaxComp.SaveAs("elMaxComp.png")
cnv_elMaxComp.Write()

cnv_elMaxQComp = rt.TCanvas("cnv_elMaxQComp","cnv_elMaxQComp")
hists_elMaxQComp = sortHists([h_elMaxQComp_CCnumu, h_elMaxQComp_NCnumu, h_elMaxQComp_CCnue, h_elMaxQComp_NCnue, h_elMaxQComp_ext])
hists_elMaxQComp[0].Draw("EHIST")
for i in range(1,len(hists_elMaxQComp)):
  hists_elMaxQComp[i].Draw("EHISTSAME")
leg_elMaxQComp = rt.TLegend(0.7,0.7,0.9,0.9)
leg_elMaxQComp = configureLegend(leg_elMaxQComp, h_elMaxQComp_CCnumu,
  h_elMaxQComp_NCnumu, h_elMaxQComp_CCnue, h_elMaxQComp_NCnue, h_elMaxQComp_ext)
leg_elMaxQComp.Draw()
#cnv_elMaxQComp.SaveAs("elMaxQComp.png")
cnv_elMaxQComp.Write()

cnv_elMaxQCosTheta = rt.TCanvas("cnv_elMaxQCosTheta","cnv_elMaxQCosTheta")
hists_elMaxQCosTheta = sortHists([h_elMaxQCosTheta_CCnumu, h_elMaxQCosTheta_NCnumu, h_elMaxQCosTheta_CCnue, h_elMaxQCosTheta_NCnue, h_elMaxQCosTheta_ext])
hists_elMaxQCosTheta[0].Draw("EHIST")
for i in range(1,len(hists_elMaxQCosTheta)):
  hists_elMaxQCosTheta[i].Draw("EHISTSAME")
leg_elMaxQCosTheta = rt.TLegend(0.7,0.7,0.9,0.9)
leg_elMaxQCosTheta = configureLegend(leg_elMaxQCosTheta, h_elMaxQCosTheta_CCnumu,
  h_elMaxQCosTheta_NCnumu, h_elMaxQCosTheta_CCnue, h_elMaxQCosTheta_NCnue, h_elMaxQCosTheta_ext)
leg_elMaxQCosTheta.Draw()
#cnv_elMaxQCosTheta.SaveAs("elMaxQCosTheta.png")
cnv_elMaxQCosTheta.Write()

cnv_elMaxQElScore = rt.TCanvas("cnv_elMaxQElScore","cnv_elMaxQElScore")
hists_elMaxQElScore = sortHists([h_elMaxQElScore_CCnumu, h_elMaxQElScore_NCnumu, h_elMaxQElScore_CCnue, h_elMaxQElScore_NCnue, h_elMaxQElScore_ext])
hists_elMaxQElScore[0].Draw("EHIST")
for i in range(1,len(hists_elMaxQElScore)):
  hists_elMaxQElScore[i].Draw("EHISTSAME")
leg_elMaxQElScore = rt.TLegend(0.7,0.7,0.9,0.9)
leg_elMaxQElScore = configureLegend(leg_elMaxQElScore, h_elMaxQElScore_CCnumu,
  h_elMaxQElScore_NCnumu, h_elMaxQElScore_CCnue, h_elMaxQElScore_NCnue, h_elMaxQElScore_ext)
leg_elMaxQElScore.Draw()
#cnv_elMaxQElScore.SaveAs("elMaxQElScore.png")
cnv_elMaxQElScore.Write()

cnv_elMaxQPhScore = rt.TCanvas("cnv_elMaxQPhScore","cnv_elMaxQPhScore")
hists_elMaxQPhScore = sortHists([h_elMaxQPhScore_CCnumu, h_elMaxQPhScore_NCnumu, h_elMaxQPhScore_CCnue, h_elMaxQPhScore_NCnue, h_elMaxQPhScore_ext])
hists_elMaxQPhScore[0].Draw("EHIST")
for i in range(1,len(hists_elMaxQPhScore)):
  hists_elMaxQPhScore[i].Draw("EHISTSAME")
leg_elMaxQPhScore = rt.TLegend(0.7,0.7,0.9,0.9)
leg_elMaxQPhScore = configureLegend(leg_elMaxQPhScore, h_elMaxQPhScore_CCnumu,
  h_elMaxQPhScore_NCnumu, h_elMaxQPhScore_CCnue, h_elMaxQPhScore_NCnue, h_elMaxQPhScore_ext)
leg_elMaxQPhScore.Draw()
#cnv_elMaxQPhScore.SaveAs("elMaxQPhScore.png")
cnv_elMaxQPhScore.Write()

cnv_elMaxQPiScore = rt.TCanvas("cnv_elMaxQPiScore","cnv_elMaxQPiScore")
hists_elMaxQPiScore = sortHists([h_elMaxQPiScore_CCnumu, h_elMaxQPiScore_NCnumu, h_elMaxQPiScore_CCnue, h_elMaxQPiScore_NCnue, h_elMaxQPiScore_ext])
hists_elMaxQPiScore[0].Draw("EHIST")
for i in range(1,len(hists_elMaxQPiScore)):
  hists_elMaxQPiScore[i].Draw("EHISTSAME")
leg_elMaxQPiScore = rt.TLegend(0.7,0.7,0.9,0.9)
leg_elMaxQPiScore = configureLegend(leg_elMaxQPiScore, h_elMaxQPiScore_CCnumu,
  h_elMaxQPiScore_NCnumu, h_elMaxQPiScore_CCnue, h_elMaxQPiScore_NCnue, h_elMaxQPiScore_ext)
leg_elMaxQPiScore.Draw()
#cnv_elMaxQPiScore.SaveAs("elMaxQPiScore.png")
cnv_elMaxQPiScore.Write()

cnv_elMaxQConf = rt.TCanvas("cnv_elMaxQConf","cnv_elMaxQConf")
hists_elMaxQConf = sortHists([h_elMaxQConf_CCnumu, h_elMaxQConf_NCnumu, h_elMaxQConf_CCnue, h_elMaxQConf_NCnue, h_elMaxQConf_ext])
hists_elMaxQConf[0].Draw("EHIST")
for i in range(1,len(hists_elMaxQConf)):
  hists_elMaxQConf[i].Draw("EHISTSAME")
leg_elMaxQConf = rt.TLegend(0.7,0.7,0.9,0.9)
leg_elMaxQConf = configureLegend(leg_elMaxQConf, h_elMaxQConf_CCnumu,
  h_elMaxQConf_NCnumu, h_elMaxQConf_CCnue, h_elMaxQConf_NCnue, h_elMaxQConf_ext)
leg_elMaxQConf.Draw()
#cnv_elMaxQConf.SaveAs("elMaxQConf.png")
cnv_elMaxQConf.Write()

cnv_elMaxQF = rt.TCanvas("cnv_elMaxQF","cnv_elMaxQF")
hists_elMaxQF = sortHists([h_elMaxQF_CCnumu, h_elMaxQF_NCnumu, h_elMaxQF_CCnue, h_elMaxQF_NCnue, h_elMaxQF_ext])
hists_elMaxQF[0].Draw("EHIST")
for i in range(1,len(hists_elMaxQF)):
  hists_elMaxQF[i].Draw("EHISTSAME")
leg_elMaxQF = rt.TLegend(0.7,0.7,0.9,0.9)
leg_elMaxQF = configureLegend(leg_elMaxQF, h_elMaxQF_CCnumu,
  h_elMaxQF_NCnumu, h_elMaxQF_CCnue, h_elMaxQF_NCnue, h_elMaxQF_ext)
leg_elMaxQF.Draw()
#cnv_elMaxQF.SaveAs("elMaxQF.png")
cnv_elMaxQF.Write()

cnv_elMaxQ = rt.TCanvas("cnv_elMaxQ","cnv_elMaxQ")
hists_elMaxQ = sortHists([h_elMaxQ_CCnumu, h_elMaxQ_NCnumu, h_elMaxQ_CCnue, h_elMaxQ_NCnue, h_elMaxQ_ext])
hists_elMaxQ[0].Draw("EHIST")
for i in range(1,len(hists_elMaxQ)):
  hists_elMaxQ[i].Draw("EHISTSAME")
leg_elMaxQ = rt.TLegend(0.7,0.7,0.9,0.9)
leg_elMaxQ = configureLegend(leg_elMaxQ, h_elMaxQ_CCnumu,
  h_elMaxQ_NCnumu, h_elMaxQ_CCnue, h_elMaxQ_NCnue, h_elMaxQ_ext)
leg_elMaxQ.Draw()
#cnv_elMaxQ.SaveAs("elMaxQ.png")
cnv_elMaxQ.Write()

cnv_elMaxQComp_wConfCut = rt.TCanvas("cnv_elMaxQComp_wConfCut","cnv_elMaxQComp_wConfCut")
hists_elMaxQComp_wConfCut = sortHists([h_elMaxQComp_wConfCut_CCnumu, h_elMaxQComp_wConfCut_NCnumu, h_elMaxQComp_wConfCut_CCnue, h_elMaxQComp_wConfCut_NCnue, h_elMaxQComp_wConfCut_ext])
hists_elMaxQComp_wConfCut[0].Draw("EHIST")
for i in range(1,len(hists_elMaxQComp_wConfCut)):
  hists_elMaxQComp_wConfCut[i].Draw("EHISTSAME")
leg_elMaxQComp_wConfCut = rt.TLegend(0.7,0.7,0.9,0.9)
leg_elMaxQComp_wConfCut = configureLegend(leg_elMaxQComp_wConfCut, h_elMaxQComp_wConfCut_CCnumu,
  h_elMaxQComp_wConfCut_NCnumu, h_elMaxQComp_wConfCut_CCnue, h_elMaxQComp_wConfCut_NCnue, h_elMaxQComp_wConfCut_ext)
leg_elMaxQComp_wConfCut.Draw()
#cnv_elMaxQComp_wConfCut.SaveAs("elMaxQComp_wConfCut.png")
cnv_elMaxQComp_wConfCut.Write()

cnv_elMaxQCosTheta_wConfCut = rt.TCanvas("cnv_elMaxQCosTheta_wConfCut","cnv_elMaxQCosTheta_wConfCut")
hists_elMaxQCosTheta_wConfCut = sortHists([h_elMaxQCosTheta_wConfCut_CCnumu, h_elMaxQCosTheta_wConfCut_NCnumu, h_elMaxQCosTheta_wConfCut_CCnue, h_elMaxQCosTheta_wConfCut_NCnue, h_elMaxQCosTheta_wConfCut_ext])
hists_elMaxQCosTheta_wConfCut[0].Draw("EHIST")
for i in range(1,len(hists_elMaxQCosTheta_wConfCut)):
  hists_elMaxQCosTheta_wConfCut[i].Draw("EHISTSAME")
leg_elMaxQCosTheta_wConfCut = rt.TLegend(0.7,0.7,0.9,0.9)
leg_elMaxQCosTheta_wConfCut = configureLegend(leg_elMaxQCosTheta_wConfCut, h_elMaxQCosTheta_wConfCut_CCnumu,
  h_elMaxQCosTheta_wConfCut_NCnumu, h_elMaxQCosTheta_wConfCut_CCnue, h_elMaxQCosTheta_wConfCut_NCnue, h_elMaxQCosTheta_wConfCut_ext)
leg_elMaxQCosTheta_wConfCut.Draw()
#cnv_elMaxQCosTheta_wConfCut.SaveAs("elMaxQCosTheta_wConfCut.png")
cnv_elMaxQCosTheta_wConfCut.Write()

cnv_elMaxQF_wConfCut = rt.TCanvas("cnv_elMaxQF_wConfCut","cnv_elMaxQF_wConfCut")
hists_elMaxQF_wConfCut = sortHists([h_elMaxQF_wConfCut_CCnumu, h_elMaxQF_wConfCut_NCnumu, h_elMaxQF_wConfCut_CCnue, h_elMaxQF_wConfCut_NCnue, h_elMaxQF_wConfCut_ext])
hists_elMaxQF_wConfCut[0].Draw("EHIST")
for i in range(1,len(hists_elMaxQF_wConfCut)):
  hists_elMaxQF_wConfCut[i].Draw("EHISTSAME")
leg_elMaxQF_wConfCut = rt.TLegend(0.7,0.7,0.9,0.9)
leg_elMaxQF_wConfCut = configureLegend(leg_elMaxQF_wConfCut, h_elMaxQF_wConfCut_CCnumu,
  h_elMaxQF_wConfCut_NCnumu, h_elMaxQF_wConfCut_CCnue, h_elMaxQF_wConfCut_NCnue, h_elMaxQF_wConfCut_ext)
leg_elMaxQF_wConfCut.Draw()
#cnv_elMaxQF_wConfCut.SaveAs("elMaxQF_wConfCut.png")
cnv_elMaxQF_wConfCut.Write()

cnv_elMaxQ_wConfCut = rt.TCanvas("cnv_elMaxQ_wConfCut","cnv_elMaxQ_wConfCut")
hists_elMaxQ_wConfCut = sortHists([h_elMaxQ_wConfCut_CCnumu, h_elMaxQ_wConfCut_NCnumu, h_elMaxQ_wConfCut_CCnue, h_elMaxQ_wConfCut_NCnue, h_elMaxQ_wConfCut_ext])
hists_elMaxQ_wConfCut[0].Draw("EHIST")
for i in range(1,len(hists_elMaxQ_wConfCut)):
  hists_elMaxQ_wConfCut[i].Draw("EHISTSAME")
leg_elMaxQ_wConfCut = rt.TLegend(0.7,0.7,0.9,0.9)
leg_elMaxQ_wConfCut = configureLegend(leg_elMaxQ_wConfCut, h_elMaxQ_wConfCut_CCnumu,
  h_elMaxQ_wConfCut_NCnumu, h_elMaxQ_wConfCut_CCnue, h_elMaxQ_wConfCut_NCnue, h_elMaxQ_wConfCut_ext)
leg_elMaxQ_wConfCut.Draw()
#cnv_elMaxQ_wConfCut.SaveAs("elMaxQ_wConfCut.png")
cnv_elMaxQ_wConfCut.Write()

cnv_elMaxComp_wQFcut = rt.TCanvas("cnv_elMaxComp_wQFcut","cnv_elMaxComp_wQFcut")
hists_elMaxComp_wQFcut = sortHists([h_elMaxComp_wQFcut_CCnumu, h_elMaxComp_wQFcut_NCnumu, h_elMaxComp_wQFcut_CCnue, h_elMaxComp_wQFcut_NCnue, h_elMaxComp_wQFcut_ext])
hists_elMaxComp_wQFcut[0].Draw("EHIST")
for i in range(1,len(hists_elMaxComp_wQFcut)):
  hists_elMaxComp_wQFcut[i].Draw("EHISTSAME")
leg_elMaxComp_wQFcut = rt.TLegend(0.7,0.7,0.9,0.9)
leg_elMaxComp_wQFcut = configureLegend(leg_elMaxComp_wQFcut, h_elMaxComp_wQFcut_CCnumu,
  h_elMaxComp_wQFcut_NCnumu, h_elMaxComp_wQFcut_CCnue, h_elMaxComp_wQFcut_NCnue, h_elMaxComp_wQFcut_ext)
leg_elMaxComp_wQFcut.Draw()
#cnv_elMaxComp_wQFcut.SaveAs("elMaxComp_wQFcut.png")
cnv_elMaxComp_wQFcut.Write()

cnv_elMaxQ_wQFcut = rt.TCanvas("cnv_elMaxQ_wQFcut","cnv_elMaxQ_wQFcut")
hists_elMaxQ_wQFcut = sortHists([h_elMaxQ_wQFcut_CCnumu, h_elMaxQ_wQFcut_NCnumu, h_elMaxQ_wQFcut_CCnue, h_elMaxQ_wQFcut_NCnue, h_elMaxQ_wQFcut_ext])
hists_elMaxQ_wQFcut[0].Draw("EHIST")
for i in range(1,len(hists_elMaxQ_wQFcut)):
  hists_elMaxQ_wQFcut[i].Draw("EHISTSAME")
leg_elMaxQ_wQFcut = rt.TLegend(0.7,0.7,0.9,0.9)
leg_elMaxQ_wQFcut = configureLegend(leg_elMaxQ_wQFcut, h_elMaxQ_wQFcut_CCnumu,
  h_elMaxQ_wQFcut_NCnumu, h_elMaxQ_wQFcut_CCnue, h_elMaxQ_wQFcut_NCnue, h_elMaxQ_wQFcut_ext)
leg_elMaxQ_wQFcut.Draw()
#cnv_elMaxQ_wQFcut.SaveAs("elMaxQ_wQFcut.png")
cnv_elMaxQ_wQFcut.Write()

cnv_elMaxComp_wQcut = rt.TCanvas("cnv_elMaxComp_wQcut","cnv_elMaxComp_wQcut")
hists_elMaxComp_wQcut = sortHists([h_elMaxComp_wQcut_CCnumu, h_elMaxComp_wQcut_NCnumu, h_elMaxComp_wQcut_CCnue, h_elMaxComp_wQcut_NCnue, h_elMaxComp_wQcut_ext])
hists_elMaxComp_wQcut[0].Draw("EHIST")
for i in range(1,len(hists_elMaxComp_wQcut)):
  hists_elMaxComp_wQcut[i].Draw("EHISTSAME")
leg_elMaxComp_wQcut = rt.TLegend(0.7,0.7,0.9,0.9)
leg_elMaxComp_wQcut = configureLegend(leg_elMaxComp_wQcut, h_elMaxComp_wQcut_CCnumu,
  h_elMaxComp_wQcut_NCnumu, h_elMaxComp_wQcut_CCnue, h_elMaxComp_wQcut_NCnue, h_elMaxComp_wQcut_ext)
leg_elMaxComp_wQcut.Draw()
#cnv_elMaxComp_wQcut.SaveAs("elMaxComp_wQcut.png")
cnv_elMaxComp_wQcut.Write()

cnv_elMaxQF_wQcut = rt.TCanvas("cnv_elMaxQF_wQcut","cnv_elMaxQF_wQcut")
hists_elMaxQF_wQcut = sortHists([h_elMaxQF_wQcut_CCnumu, h_elMaxQF_wQcut_NCnumu, h_elMaxQF_wQcut_CCnue, h_elMaxQF_wQcut_NCnue, h_elMaxQF_wQcut_ext])
hists_elMaxQF_wQcut[0].Draw("EHIST")
for i in range(1,len(hists_elMaxQF_wQcut)):
  hists_elMaxQF_wQcut[i].Draw("EHISTSAME")
leg_elMaxQF_wQcut = rt.TLegend(0.7,0.7,0.9,0.9)
leg_elMaxQF_wQcut = configureLegend(leg_elMaxQF_wQcut, h_elMaxQF_wQcut_CCnumu,
  h_elMaxQF_wQcut_NCnumu, h_elMaxQF_wQcut_CCnue, h_elMaxQF_wQcut_NCnue, h_elMaxQF_wQcut_ext)
leg_elMaxQF_wQcut.Draw()
#cnv_elMaxQF_wQcut.SaveAs("elMaxQF_wQcut.png")
cnv_elMaxQF_wQcut.Write()

cnv_muMaxComp = rt.TCanvas("cnv_muMaxComp","cnv_muMaxComp")
hists_muMaxComp = sortHists([h_muMaxComp_CCnumu, h_muMaxComp_NCnumu, h_muMaxComp_CCnue, h_muMaxComp_NCnue, h_muMaxComp_ext])
hists_muMaxComp[0].Draw("EHIST")
for i in range(1,len(hists_muMaxComp)):
  hists_muMaxComp[i].Draw("EHISTSAME")
leg_muMaxComp = rt.TLegend(0.7,0.7,0.9,0.9)
leg_muMaxComp = configureLegend(leg_muMaxComp, h_muMaxComp_CCnumu,
  h_muMaxComp_NCnumu, h_muMaxComp_CCnue, h_muMaxComp_NCnue, h_muMaxComp_ext)
leg_muMaxComp.Draw()
#cnv_muMaxComp.SaveAs("muMaxComp.png")
cnv_muMaxComp.Write()

cnv_muMaxQComp = rt.TCanvas("cnv_muMaxQComp","cnv_muMaxQComp")
hists_muMaxQComp = sortHists([h_muMaxQComp_CCnumu, h_muMaxQComp_NCnumu, h_muMaxQComp_CCnue, h_muMaxQComp_NCnue, h_muMaxQComp_ext])
hists_muMaxQComp[0].Draw("EHIST")
for i in range(1,len(hists_muMaxQComp)):
  hists_muMaxQComp[i].Draw("EHISTSAME")
leg_muMaxQComp = rt.TLegend(0.7,0.7,0.9,0.9)
leg_muMaxQComp = configureLegend(leg_muMaxQComp, h_muMaxQComp_CCnumu,
  h_muMaxQComp_NCnumu, h_muMaxQComp_CCnue, h_muMaxQComp_NCnue, h_muMaxQComp_ext)
leg_muMaxQComp.Draw()
#cnv_muMaxQComp.SaveAs("muMaxQComp.png")
cnv_muMaxQComp.Write()

cnv_muMaxQCosTheta = rt.TCanvas("cnv_muMaxQCosTheta","cnv_muMaxQCosTheta")
hists_muMaxQCosTheta = sortHists([h_muMaxQCosTheta_CCnumu, h_muMaxQCosTheta_NCnumu, h_muMaxQCosTheta_CCnue, h_muMaxQCosTheta_NCnue, h_muMaxQCosTheta_ext])
hists_muMaxQCosTheta[0].Draw("EHIST")
for i in range(1,len(hists_muMaxQCosTheta)):
  hists_muMaxQCosTheta[i].Draw("EHISTSAME")
leg_muMaxQCosTheta = rt.TLegend(0.7,0.7,0.9,0.9)
leg_muMaxQCosTheta = configureLegend(leg_muMaxQCosTheta, h_muMaxQCosTheta_CCnumu,
  h_muMaxQCosTheta_NCnumu, h_muMaxQCosTheta_CCnue, h_muMaxQCosTheta_NCnue, h_muMaxQCosTheta_ext)
leg_muMaxQCosTheta.Draw()
#cnv_muMaxQCosTheta.SaveAs("muMaxQCosTheta.png")
cnv_muMaxQCosTheta.Write()

cnv_muMaxQF = rt.TCanvas("cnv_muMaxQF","cnv_muMaxQF")
hists_muMaxQF = sortHists([h_muMaxQF_CCnumu, h_muMaxQF_NCnumu, h_muMaxQF_CCnue, h_muMaxQF_NCnue, h_muMaxQF_ext])
hists_muMaxQF[0].Draw("EHIST")
for i in range(1,len(hists_muMaxQF)):
  hists_muMaxQF[i].Draw("EHISTSAME")
leg_muMaxQF = rt.TLegend(0.7,0.7,0.9,0.9)
leg_muMaxQF = configureLegend(leg_muMaxQF, h_muMaxQF_CCnumu,
  h_muMaxQF_NCnumu, h_muMaxQF_CCnue, h_muMaxQF_NCnue, h_muMaxQF_ext)
leg_muMaxQF.Draw()
#cnv_muMaxQF.SaveAs("muMaxQF.png")
cnv_muMaxQF.Write()

cnv_muMaxQ = rt.TCanvas("cnv_muMaxQ","cnv_muMaxQ")
hists_muMaxQ = sortHists([h_muMaxQ_CCnumu, h_muMaxQ_NCnumu, h_muMaxQ_CCnue, h_muMaxQ_NCnue, h_muMaxQ_ext])
hists_muMaxQ[0].Draw("EHIST")
for i in range(1,len(hists_muMaxQ)):
  hists_muMaxQ[i].Draw("EHISTSAME")
leg_muMaxQ = rt.TLegend(0.7,0.7,0.9,0.9)
leg_muMaxQ = configureLegend(leg_muMaxQ, h_muMaxQ_CCnumu,
  h_muMaxQ_NCnumu, h_muMaxQ_CCnue, h_muMaxQ_NCnue, h_muMaxQ_ext)
leg_muMaxQ.Draw()
#cnv_muMaxQ.SaveAs("muMaxQ.png")
cnv_muMaxQ.Write()


