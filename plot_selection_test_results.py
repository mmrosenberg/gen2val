
import sys, argparse
import numpy as np
import ROOT as rt

from math import isinf
from helpers.plotting_functions import sortHists
from helpers.larflowreco_ana_funcs import isFiducial, isFiducialWC, getDistance


parser = argparse.ArgumentParser("Plot Selection Test Results")
parser.add_argument("-fnu", "--bnbnu_file", type=str, default="selection_output/prepare_selection_test_output/prepare_selection_test_reco_v2me05_gen2val_v23_nu_file.root", help="bnb nu input file")
parser.add_argument("-fnue", "--bnbnue_file", type=str, default="selection_output/prepare_selection_test_output/prepare_selection_test_reco_v2me05_gen2val_v23_nue_file.root", help="bnb nu input file")
parser.add_argument("-fext", "--extbnb_file", type=str, default="flat_ntuples/dlgen2_reco_v2me06_ntuple_v1_mcc9_v29e_dl_run3_G1_extbnb_partial.root", help="extbnb input file")
parser.add_argument("-fdata", "--data_file", type=str, default="selection_output/prepare_selection_test_output/prepare_selection_test_reco_v2me05_gen2val_v22_bnb5e19_file.root", help="bnb data input file")
parser.add_argument("-vfc", "--vertexFracOnCosCut", type=float, default=1., help="vtxFracHitsOnCosmic cut")
parser.add_argument("-d", "--distCut", type=float, default=9999., help="distance to vertex cut value")
parser.add_argument("-c", "--compCut", type=float, default=0., help="completeness cut value")
parser.add_argument("-p", "--purityCut", type=float, default=0., help="purity cut value")
parser.add_argument("-s", "--confCut", type=float, default=7.3, help="electron class confidence cut value")
parser.add_argument("-q", "--chargeCut", type=float, default=0, help="electron charge fraction cut value")
parser.add_argument("-qf", "--chargeFracCut", type=float, default=0., help="electron charge fraction cut value")
parser.add_argument("-t", "--cosThetaCut", type=float, default=0., help="cos(angle to beam) cut value")
parser.add_argument("-o", "--outfile", type=str, default="selection_output/plot_selection_test_results_output.root", help="output root file name")
parser.add_argument("--recoEOverflow", help="plot overflow bin for final recoE selection plot", action="store_true")
parser.add_argument("--oldVertexVar", help="use old nVertices variable for vertex found cut", action="store_true")
parser.add_argument("--makeKPplots", help="make keypoint plots", action="store_true")
parser.add_argument("--smallFV", help="use 20cm fiducial volume", action="store_true")
parser.add_argument("--NPMLCosmicConfig", help="used extBNB file and incorrect scaling from NPML plots", action="store_true")
parser.add_argument("--printEXTinfo", help="print event info for selected cosmic background", action="store_true")
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

if args.NPMLCosmicConfig:
  fext = rt.TFile("selection_output/prepare_selection_test_output/prepare_selection_test_reco_v2me05_gen2val_v22_extbnb_file.root")
else:
  fext = rt.TFile(args.extbnb_file)
text = fext.Get("EventTree")

fdata = rt.TFile(args.data_file)
tdata = fdata.Get("EventTree")

out_data_sel_evts = open("selection_output/selected_bnb5e19_events.txt","w")
out_data_sel_evts.write("fileid run subrun event recoNuE\n")

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

#970: number of run3 extbnb merged_dlana files in prepare_selection_test_reco_v2me05_gen2val_v16_extbnb_file.root
#17661: number of run3 extbnb merged_dlana files in prepare_selection_test_reco_v2me05_gen2val_v17+_extbnb_file.root
#89559: number of run3 extbnb files (# in def: prod_extunbiased_swizzle_crt_inclusive_v6_v6a_goodruns_mcc9_run3)
#textPOTsum = (17661./89559.)*run3POT
#from Zarko's POT counting for selection_output/prepare_selection_test_output/prepare_selection_test_reco_v2me05_gen2val_v22_extbnb_file.root:
#textPOTsum = 1.3298521464785359e+19
#from Zarko's POT counting for flat_ntuples/dlgen2_reco_v2me06_ntuple_v1_mcc9_v29e_dl_run3_G1_extbnb_partial.root
#textPOTsum = 2.561872704628622e+19

if args.NPMLCosmicConfig:
  textPOTsum = (17661./89559.)*run3POT
else:
  if args.extbnb_file == "selection_output/prepare_selection_test_output/prepare_selection_test_reco_v2me05_gen2val_v22_extbnb_file.root":
    textPOTsum = 1.3298521464785359e+19
  elif "dlgen2_reco_v2me06_ntuple_v1_mcc9_v29e_dl_run3_G1_extbnb_partial" in args.extbnb_file:
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

h_nMu_CCnumu = rt.TH1F("h_nMu_CCnumu","Number of Reco Muons",10,0,10)
h_nMu_NCnumu = rt.TH1F("h_nMu_NCnumu","Number of Reco Muons",10,0,10)
h_nMu_CCnue = rt.TH1F("h_nMu_CCnue","Number of Reco Muons",10,0,10)
h_nMu_NCnue = rt.TH1F("h_nMu_NCnue","Number of Reco Muons",10,0,10)
h_nMu_ext = rt.TH1F("h_nMu_ext","Number of Reco Muons",10,0,10)
h_nMu_CCnumu, h_nMu_NCnumu, h_nMu_CCnue, h_nMu_NCnue, h_nMu_ext = configureHists(h_nMu_CCnumu,
 h_nMu_NCnumu, h_nMu_CCnue, h_nMu_NCnue, h_nMu_ext)

if args.makeKPplots:
  h_cosKPDist_CCnumu = rt.TH1F("h_cosKPDist_CCnumu","Distance to Nearest Cosmic KP (cm)",500,0,1000)
  h_cosKPDist_NCnumu = rt.TH1F("h_cosKPDist_NCnumu","Distance to Nearest Cosmic KP (cm)",500,0,1000)
  h_cosKPDist_CCnue = rt.TH1F("h_cosKPDist_CCnue","Distance to Nearest Cosmic KP (cm)",500,0,1000)
  h_cosKPDist_NCnue = rt.TH1F("h_cosKPDist_NCnue","Distance to Nearest Cosmic KP (cm)",500,0,1000)
  h_cosKPDist_ext = rt.TH1F("h_cosKPDist_ext","Distance to Nearest Cosmic KP (cm)",500,0,1000)
  h_cosKPDist_CCnumu, h_cosKPDist_NCnumu, h_cosKPDist_CCnue, h_cosKPDist_NCnue, h_cosKPDist_ext = configureHists(h_cosKPDist_CCnumu,
   h_cosKPDist_NCnumu, h_cosKPDist_CCnue, h_cosKPDist_NCnue, h_cosKPDist_ext)
  
  h_tEndKPDist_CCnumu = rt.TH1F("h_tEndKPDist_CCnumu","Distance to Nearest Track End KP (cm)",500,0,1000)
  h_tEndKPDist_NCnumu = rt.TH1F("h_tEndKPDist_NCnumu","Distance to Nearest Track End KP (cm)",500,0,1000)
  h_tEndKPDist_CCnue = rt.TH1F("h_tEndKPDist_CCnue","Distance to Nearest Track End KP (cm)",500,0,1000)
  h_tEndKPDist_NCnue = rt.TH1F("h_tEndKPDist_NCnue","Distance to Nearest Track End KP (cm)",500,0,1000)
  h_tEndKPDist_ext = rt.TH1F("h_tEndKPDist_ext","Distance to Nearest Track End KP (cm)",500,0,1000)
  h_tEndKPDist_CCnumu, h_tEndKPDist_NCnumu, h_tEndKPDist_CCnue, h_tEndKPDist_NCnue, h_tEndKPDist_ext = configureHists(h_tEndKPDist_CCnumu,
   h_tEndKPDist_NCnumu, h_tEndKPDist_CCnue, h_tEndKPDist_NCnue, h_tEndKPDist_ext)
  
  h_cosKPDist_wConfCut_CCnumu = rt.TH1F("h_cosKPDist_wConfCut_CCnumu","Distance to Nearest Cosmic KP (cm)",500,0,1000)
  h_cosKPDist_wConfCut_NCnumu = rt.TH1F("h_cosKPDist_wConfCut_NCnumu","Distance to Nearest Cosmic KP (cm)",500,0,1000)
  h_cosKPDist_wConfCut_CCnue = rt.TH1F("h_cosKPDist_wConfCut_CCnue","Distance to Nearest Cosmic KP (cm)",500,0,1000)
  h_cosKPDist_wConfCut_NCnue = rt.TH1F("h_cosKPDist_wConfCut_NCnue","Distance to Nearest Cosmic KP (cm)",500,0,1000)
  h_cosKPDist_wConfCut_ext = rt.TH1F("h_cosKPDist_wConfCut_ext","Distance to Nearest Cosmic KP (cm)",500,0,1000)
  h_cosKPDist_wConfCut_CCnumu, h_cosKPDist_wConfCut_NCnumu, h_cosKPDist_wConfCut_CCnue, h_cosKPDist_wConfCut_NCnue, h_cosKPDist_wConfCut_ext = configureHists(h_cosKPDist_wConfCut_CCnumu,
   h_cosKPDist_wConfCut_NCnumu, h_cosKPDist_wConfCut_CCnue, h_cosKPDist_wConfCut_NCnue, h_cosKPDist_wConfCut_ext)
  
  h_tEndKPDist_wConfCut_CCnumu = rt.TH1F("h_tEndKPDist_wConfCut_CCnumu","Distance to Nearest Track End KP (cm)",500,0,1000)
  h_tEndKPDist_wConfCut_NCnumu = rt.TH1F("h_tEndKPDist_wConfCut_NCnumu","Distance to Nearest Track End KP (cm)",500,0,1000)
  h_tEndKPDist_wConfCut_CCnue = rt.TH1F("h_tEndKPDist_wConfCut_CCnue","Distance to Nearest Track End KP (cm)",500,0,1000)
  h_tEndKPDist_wConfCut_NCnue = rt.TH1F("h_tEndKPDist_wConfCut_NCnue","Distance to Nearest Track End KP (cm)",500,0,1000)
  h_tEndKPDist_wConfCut_ext = rt.TH1F("h_tEndKPDist_wConfCut_ext","Distance to Nearest Track End KP (cm)",500,0,1000)
  h_tEndKPDist_wConfCut_CCnumu, h_tEndKPDist_wConfCut_NCnumu, h_tEndKPDist_wConfCut_CCnue, h_tEndKPDist_wConfCut_NCnue, h_tEndKPDist_wConfCut_ext = configureHists(h_tEndKPDist_wConfCut_CCnumu,
   h_tEndKPDist_wConfCut_NCnumu, h_tEndKPDist_wConfCut_CCnue, h_tEndKPDist_wConfCut_NCnue, h_tEndKPDist_wConfCut_ext)

h_maxElConf_CCnumu = rt.TH1F("h_maxElConf_CCnumu","Max \"Electron Confidence Score\"",21,-1,20)
h_maxElConf_NCnumu = rt.TH1F("h_maxElConf_NCnumu","Max \"Electron Confidence Score\"",21,-1,20)
h_maxElConf_CCnue = rt.TH1F("h_maxElConf_CCnue","Max \"Electron Confidence Score\"",21,-1,20)
h_maxElConf_NCnue = rt.TH1F("h_maxElConf_NCnue","Max \"Electron Confidence Score\"",21,-1,20)
h_maxElConf_ext = rt.TH1F("h_maxElConf_ext","Max \"Electron Confidence Score\"",21,-1,20)
h_maxElConf_CCnumu, h_maxElConf_NCnumu, h_maxElConf_CCnue, h_maxElConf_NCnue, h_maxElConf_ext = configureHists(h_maxElConf_CCnumu,
 h_maxElConf_NCnumu, h_maxElConf_CCnue, h_maxElConf_NCnue, h_maxElConf_ext)

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

h_elMaxPur_CCnumu = rt.TH1F("h_elMaxPur_CCnumu","Max Electron Purity",21,0,1.05)
h_elMaxPur_NCnumu = rt.TH1F("h_elMaxPur_NCnumu","Max Electron Purity",21,0,1.05)
h_elMaxPur_CCnue = rt.TH1F("h_elMaxPur_CCnue","Max Electron Purity",21,0,1.05)
h_elMaxPur_NCnue = rt.TH1F("h_elMaxPur_NCnue","Max Electron Purity",21,0,1.05)
h_elMaxPur_ext = rt.TH1F("h_elMaxPur_ext","Max Electron Purity",21,0,1.05)
h_elMaxPur_CCnumu, h_elMaxPur_NCnumu, h_elMaxPur_CCnue, h_elMaxPur_NCnue, h_elMaxPur_ext = configureHists(h_elMaxPur_CCnumu,
 h_elMaxPur_NCnumu, h_elMaxPur_CCnue, h_elMaxPur_NCnue, h_elMaxPur_ext)

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

h_elMaxQPur_CCnumu = rt.TH1F("h_elMaxQPur_CCnumu","Purity for Largest Electron Shower",21,0,1.05)
h_elMaxQPur_NCnumu = rt.TH1F("h_elMaxQPur_NCnumu","Purity for Largest Electron Shower",21,0,1.05)
h_elMaxQPur_CCnue = rt.TH1F("h_elMaxQPur_CCnue","Purity for Largest Electron Shower",21,0,1.05)
h_elMaxQPur_NCnue = rt.TH1F("h_elMaxQPur_NCnue","Purity for Largest Electron Shower",21,0,1.05)
h_elMaxQPur_ext = rt.TH1F("h_elMaxQPur_ext","Purity for Largest Electron Shower",21,0,1.05)
h_elMaxQPur_CCnumu, h_elMaxQPur_NCnumu, h_elMaxQPur_CCnue, h_elMaxQPur_NCnue, h_elMaxQPur_ext = configureHists(h_elMaxQPur_CCnumu,
 h_elMaxQPur_NCnumu, h_elMaxQPur_CCnue, h_elMaxQPur_NCnue, h_elMaxQPur_ext)

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

h_elMaxQVtxDist_wConfCut_CCnumu = rt.TH1F("h_elMaxQVtxDist_wConfCut_CCnumu","Distance to Vertex for Largest Electron Shower",50,0,100)
h_elMaxQVtxDist_wConfCut_NCnumu = rt.TH1F("h_elMaxQVtxDist_wConfCut_NCnumu","Distance to Vertex for Largest Electron Shower",50,0,100)
h_elMaxQVtxDist_wConfCut_CCnue = rt.TH1F("h_elMaxQVtxDist_wConfCut_CCnue","Distance to Vertex for Largest Electron Shower",50,0,100)
h_elMaxQVtxDist_wConfCut_NCnue = rt.TH1F("h_elMaxQVtxDist_wConfCut_NCnue","Distance to Vertex for Largest Electron Shower",50,0,100)
h_elMaxQVtxDist_wConfCut_ext = rt.TH1F("h_elMaxQVtxDist_wConfCut_ext","Distance to Vertex for Largest Electron Shower",50,0,100)
h_elMaxQVtxDist_wConfCut_CCnumu, h_elMaxQVtxDist_wConfCut_NCnumu, h_elMaxQVtxDist_wConfCut_CCnue, h_elMaxQVtxDist_wConfCut_NCnue, h_elMaxQVtxDist_wConfCut_ext = configureHists(h_elMaxQVtxDist_wConfCut_CCnumu,
 h_elMaxQVtxDist_wConfCut_NCnumu, h_elMaxQVtxDist_wConfCut_CCnue, h_elMaxQVtxDist_wConfCut_NCnue, h_elMaxQVtxDist_wConfCut_ext)

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

h_nuE_CCnue_nCuts = rt.TH1F("h_nuE_CCnue_nCuts","Neutrino Energy for True CCnue Events",30,0,3)
h_nuE_CCnue_wCuts = rt.TH1F("h_nuE_CCnue_wCuts","Neutrino Energy for True CCnue Events",30,0,3)
h_nuE_CCnumu_wCuts = rt.TH1F("h_nuE_CCnumu_wCuts","Neutrino Energy for True CCnumu Events",30,0,3)
h_nuE_NCnumu_wCuts = rt.TH1F("h_nuE_NCnumu_wCuts","Neutrino Energy for True NCnumu Events",30,0,3)
h_nuE_NCnue_wCuts = rt.TH1F("h_nuE_NCnue_wCuts","Neutrino Energy for True NCnue Events",30,0,3)
#h_nuE_ext_wCuts = rt.TH1F("h_nuE_ext_wCuts","Neutrino Energy for ExtBNB Events",30,0,3)
h_nuE_all_wCuts = rt.TH1F("h_nuE_all_wCuts","Neutrino Energy for all Events",30,0,3)
h_nuE_CCnue_nCuts.GetYaxis().SetTitle("events per "+targetPOTstring+" POT")
h_nuE_CCnue_nCuts.SetLineColor(rt.kBlue)
h_nuE_CCnue_nCuts.SetLineWidth(2)
h_nuE_CCnue_nCuts.GetXaxis().SetTitle("true neutrino energy (GeV)")
h_nuE_CCnue_wCuts.GetYaxis().SetTitle("events per "+targetPOTstring+" POT")
h_nuE_CCnue_wCuts.SetLineColor(rt.kRed)
h_nuE_CCnue_wCuts.SetLineWidth(2)
h_nuE_CCnue_wCuts.GetXaxis().SetTitle("true neutrino energy (GeV)")
h_nuE_all_wCuts.GetYaxis().SetTitle("events per "+targetPOTstring+" POT")
h_nuE_all_wCuts.SetLineColor(rt.kRed)
h_nuE_all_wCuts.SetLineWidth(2)
h_nuE_all_wCuts.GetXaxis().SetTitle("true neutrino energy (GeV)")
h_nuE_CCnue_eff = rt.TH1F("h_nuE_CCnue_eff","Inclusive CC nue Selection MC Predictions",30,0,3)
h_nuE_CCnue_eff.GetXaxis().SetTitle("true neutrino energy (GeV)")
h_nuE_CCnue_eff.SetLineColor(rt.kBlack)
h_nuE_CCnue_eff.SetLineWidth(2)
h_nuE_CCnue_pur = rt.TH1F("h_nuE_CCnue_pur","Inclusive CC nue Selection MC Predictions",30,0,3)
h_nuE_CCnue_pur.GetXaxis().SetTitle("true neutrino energy (GeV)")
h_nuE_CCnue_pur.SetLineColor(8)
h_nuE_CCnue_pur.SetLineWidth(2)

h_nuEr_CCnue_nCuts = rt.TH1F("h_nuEr_CCnue_nCuts","Neutrino Energy for True CCnue Events",30,0,3)
h_nuEr_CCnue_wCuts = rt.TH1F("h_nuEr_CCnue_wCuts","Neutrino Energy for True CCnue Events",30,0,3)
h_nuEr_CCnumu_wCuts = rt.TH1F("h_nuEr_CCnumu_wCuts","Neutrino Energy for True CCnumu Events",30,0,3)
h_nuEr_NCnumu_wCuts = rt.TH1F("h_nuEr_NCnumu_wCuts","Neutrino Energy for True NCnumu Events",30,0,3)
h_nuEr_NCnue_wCuts = rt.TH1F("h_nuEr_NCnue_wCuts","Neutrino Energy for True NCnue Events",30,0,3)
h_nuEr_ext_wCuts = rt.TH1F("h_nuEr_ext_wCuts","Neutrino Energy for ExtBNB Events",30,0,3)
h_nuEr_all_wCuts = rt.TH1F("h_nuEr_all_wCuts","Neutrino Energy for all Events",30,0,3)
h_nuEr_CCnue_nCuts.GetYaxis().SetTitle("events per "+targetPOTstring+" POT")
h_nuEr_CCnue_nCuts.SetLineColor(rt.kBlue)
h_nuEr_CCnue_nCuts.SetLineWidth(2)
h_nuEr_CCnue_nCuts.GetXaxis().SetTitle("reco neutrino energy (GeV)")
h_nuEr_CCnue_wCuts.GetYaxis().SetTitle("events per "+targetPOTstring+" POT")
h_nuEr_CCnue_wCuts.SetLineColor(rt.kRed)
h_nuEr_CCnue_wCuts.SetLineWidth(2)
h_nuEr_CCnue_wCuts.GetXaxis().SetTitle("reco neutrino energy (GeV)")
h_nuEr_all_wCuts.GetYaxis().SetTitle("events per "+targetPOTstring+" POT")
h_nuEr_all_wCuts.SetLineColor(rt.kRed)
h_nuEr_all_wCuts.SetLineWidth(2)
h_nuEr_all_wCuts.GetXaxis().SetTitle("reco neutrino energy (GeV)")
h_nuEr_CCnue_eff = rt.TH1F("h_nuEr_CCnue_eff","Inclusive CC nue Selection MC Predictions",30,0,3)
h_nuEr_CCnue_eff.GetXaxis().SetTitle("reco neutrino energy (GeV)")
h_nuEr_CCnue_eff.SetLineColor(rt.kBlack)
h_nuEr_CCnue_eff.SetLineWidth(2)
h_nuEr_CCnue_pur = rt.TH1F("h_nuEr_CCnue_pur","Inclusive CC nue Selection MC Predictions",30,0,3)
h_nuEr_CCnue_pur.GetXaxis().SetTitle("reco neutrino energy (GeV)")
h_nuEr_CCnue_pur.SetLineColor(8)
h_nuEr_CCnue_pur.SetLineWidth(2)

if args.recoEOverflow:
  h_visE_CCnue_wCuts = rt.TH1F("h_visE_CCnue_wCuts","Reco Nu Energy for True CCnue Events",14,0,2.8)
  h_visE_CCnumu_wCuts = rt.TH1F("h_visE_CCnumu_wCuts","Reco Nu Energy for True CCnumu Events",14,0,2.8)
  h_visE_NCnumu_wCuts = rt.TH1F("h_visE_NCnumu_wCuts","Reco Nu Energy for True NCnumu Events",14,0,2.8)
  h_visE_NCnue_wCuts = rt.TH1F("h_visE_NCnue_wCuts","Reco Nu Energy for True NCnue Events",14,0,2.8)
  h_visE_ext_wCuts = rt.TH1F("h_visE_ext_wCuts","Reco Nu Energy for ExtBNB Events",14,0,2.8)
  h_visE_data_wCuts = rt.TH1F("h_visE_data_wCuts","Reco Nu Energy for BNB Data Events",14,0,2.8)
else:
  h_visE_CCnue_wCuts = rt.TH1F("h_visE_CCnue_wCuts","Reco Nu Energy for True CCnue Events",30,0,6)
  h_visE_CCnumu_wCuts = rt.TH1F("h_visE_CCnumu_wCuts","Reco Nu Energy for True CCnumu Events",30,0,6)
  h_visE_NCnumu_wCuts = rt.TH1F("h_visE_NCnumu_wCuts","Reco Nu Energy for True NCnumu Events",30,0,6)
  h_visE_NCnue_wCuts = rt.TH1F("h_visE_NCnue_wCuts","Reco Nu Energy for True NCnue Events",30,0,6)
  h_visE_ext_wCuts = rt.TH1F("h_visE_ext_wCuts","Reco Nu Energy for ExtBNB Events",30,0,6)
  h_visE_data_wCuts = rt.TH1F("h_visE_data_wCuts","Reco Nu Energy for BNB Data Events",30,0,6)
h_visE_CCnue_wCuts.SetFillColor(rt.kRed)
h_visE_NCnue_wCuts.SetFillColor(8)
h_visE_CCnumu_wCuts.SetFillColor(rt.kBlue)
h_visE_NCnumu_wCuts.SetFillColor(40)
h_visE_ext_wCuts.SetFillColor(12)
h_visE_data_wCuts.SetLineColor(rt.kBlack)
h_visE_data_wCuts.SetLineWidth(2)
h_visE_data_wCuts.SetTitle("Inclusive CCnue Selected Events")
h_visE_data_wCuts.GetYaxis().SetTitle("events per "+targetPOTstring+" POT")
h_visE_data_wCuts.GetXaxis().SetTitle("reconstructed neutrino energy (GeV)")
h_visE_all_wCuts = rt.THStack("h_visE_all_wCuts", "Inclusive CCnue Selected Events")

h_phScVsPrSum_allMC_bkg = rt.TH2F("h_phScVsPrSum_allMC_bkg","Photon Score vs. Purity for Largest e- Shower in Neutrino Background",21,0,1.05,21,-20,1)
h_phScVsPrSum_allMC_bkg.GetXaxis().SetTitle("total simulated particle purity")
h_phScVsPrSum_allMC_bkg.GetYaxis().SetTitle("photon score")

h_piScVsPr_CCnumu_bkg = rt.TH2F("h_piScVsPr_CCnumu_bkg","CCnumu Background",21,0,1.05,21,-20,1)
h_piScVsPr_CCnumu_bkg.GetXaxis().SetTitle("true pion purity")
h_piScVsPr_CCnumu_bkg.GetYaxis().SetTitle("pion score")
h_phScVsPr_CCnumu_bkg = rt.TH2F("h_phScVsPr_CCnumu_bkg","CCnumu Background",21,0,1.05,21,-20,1)
h_phScVsPr_CCnumu_bkg.GetXaxis().SetTitle("true photon purity")
h_phScVsPr_CCnumu_bkg.GetYaxis().SetTitle("photon score")
h_elScVsPr_CCnumu_bkg = rt.TH2F("h_elScVsPr_CCnumu_bkg","CCnumu Background",21,0,1.05,21,-20,1)
h_elScVsPr_CCnumu_bkg.GetXaxis().SetTitle("true electron purity")
h_elScVsPr_CCnumu_bkg.GetYaxis().SetTitle("electron score")

h_piScVsPr_NCnumu_bkg = rt.TH2F("h_piScVsPr_NCnumu_bkg","NCnumu Background",21,0,1.05,21,-20,1)
h_piScVsPr_NCnumu_bkg.GetXaxis().SetTitle("true pion purity")
h_piScVsPr_NCnumu_bkg.GetYaxis().SetTitle("pion score")
h_phScVsPr_NCnumu_bkg = rt.TH2F("h_phScVsPr_NCnumu_bkg","NCnumu Background",21,0,1.05,21,-20,1)
h_phScVsPr_NCnumu_bkg.GetXaxis().SetTitle("true photon purity")
h_phScVsPr_NCnumu_bkg.GetYaxis().SetTitle("photon score")
h_elScVsPr_NCnumu_bkg = rt.TH2F("h_elScVsPr_NCnumu_bkg","NCnumu Background",21,0,1.05,21,-20,1)
h_elScVsPr_NCnumu_bkg.GetXaxis().SetTitle("true electron purity")
h_elScVsPr_NCnumu_bkg.GetYaxis().SetTitle("electron score")

h_piScVsPr_NCnue_bkg = rt.TH2F("h_piScVsPr_NCnue_bkg","NCnue Background",21,0,1.05,21,-20,1)
h_piScVsPr_NCnue_bkg.GetXaxis().SetTitle("true pion purity")
h_piScVsPr_NCnue_bkg.GetYaxis().SetTitle("pion score")
h_phScVsPr_NCnue_bkg = rt.TH2F("h_phScVsPr_NCnue_bkg","NCnue Background",21,0,1.05,21,-20,1)
h_phScVsPr_NCnue_bkg.GetXaxis().SetTitle("true photon purity")
h_phScVsPr_NCnue_bkg.GetYaxis().SetTitle("photon score")
h_elScVsPr_NCnue_bkg = rt.TH2F("h_elScVsPr_NCnue_bkg","NCnue Background",21,0,1.05,21,-20,1)
h_elScVsPr_NCnue_bkg.GetXaxis().SetTitle("true electron purity")
h_elScVsPr_NCnue_bkg.GetYaxis().SetTitle("electron score")

h_piScElPrGr0_allMC_bkg = rt.TH1F("h_piScElPrGr0_allMC_bkg","Particle Scores for Largest e- Shower in Neutrino Background",21,-20,1)
h_piScElPrEq0_allMC_bkg = rt.TH1F("h_piScElPrEq0_allMC_bkg","Particle Scores for Largest e- Shower in Neutrino Background",21,-20,1)
h_piScElPrGr0_allMC_bkg.GetYaxis().SetTitle("events per "+targetPOTstring+" POT")
h_piScElPrGr0_allMC_bkg.SetLineWidth(2)
h_piScElPrGr0_allMC_bkg.SetLineColor(rt.kBlue)
h_piScElPrEq0_allMC_bkg.GetYaxis().SetTitle("events per "+targetPOTstring+" POT")
h_piScElPrEq0_allMC_bkg.SetLineWidth(2)
h_piScElPrEq0_allMC_bkg.SetLineColor(rt.kBlue)
h_piScElPrEq0_allMC_bkg.SetLineStyle(rt.kDashed)

h_muScElPrGr0_allMC_bkg = rt.TH1F("h_muScElPrGr0_allMC_bkg","Particle Scores for Largest e- Shower in Neutrino Background",21,-20,1)
h_muScElPrEq0_allMC_bkg = rt.TH1F("h_muScElPrEq0_allMC_bkg","Particle Scores for Largest e- Shower in Neutrino Background",21,-20,1)
h_muScElPrGr0_allMC_bkg.GetYaxis().SetTitle("events per "+targetPOTstring+" POT")
h_muScElPrGr0_allMC_bkg.SetLineWidth(2)
h_muScElPrGr0_allMC_bkg.SetLineColor(rt.kRed)
h_muScElPrEq0_allMC_bkg.GetYaxis().SetTitle("events per "+targetPOTstring+" POT")
h_muScElPrEq0_allMC_bkg.SetLineWidth(2)
h_muScElPrEq0_allMC_bkg.SetLineColor(rt.kRed)
h_muScElPrEq0_allMC_bkg.SetLineStyle(rt.kDashed)

h_phScPrSumGr0_allMC_bkg = rt.TH1F("h_phScPrSumGr0_allMC_bkg","Particle Scores for Largest e- Shower in Neutrino Background",21,-20,1)
h_phScPrSumEq0_allMC_bkg = rt.TH1F("h_phScPrSumEq0_allMC_bkg","Particle Scores for Largest e- Shower in Neutrino Background",21,-20,1)
h_phScPrSumGr0_allMC_bkg.GetYaxis().SetTitle("events per "+targetPOTstring+" POT")
h_phScPrSumGr0_allMC_bkg.SetLineWidth(2)
h_phScPrSumGr0_allMC_bkg.SetLineColor(8)
h_phScPrSumEq0_allMC_bkg.GetYaxis().SetTitle("events per "+targetPOTstring+" POT")
h_phScPrSumEq0_allMC_bkg.SetLineWidth(2)
h_phScPrSumEq0_allMC_bkg.SetLineColor(8)
h_phScPrSumEq0_allMC_bkg.SetLineStyle(rt.kDashed)

h_piScPrGr0_allMC_bkg = rt.TH1F("h_piScPrGr0_allMC_bkg","Particle Scores for Largest e- Shower in Neutrino Background",21,-20,1)
h_piScPrEq0_allMC_bkg = rt.TH1F("h_piScPrEq0_allMC_bkg","Particle Scores for Largest e- Shower in Neutrino Background",21,-20,1)
h_piScPrGr0_allMC_bkg.GetYaxis().SetTitle("events per "+targetPOTstring+" POT")
h_piScPrGr0_allMC_bkg.SetLineWidth(2)
h_piScPrGr0_allMC_bkg.SetLineColor(rt.kBlue)
h_piScPrEq0_allMC_bkg.GetYaxis().SetTitle("events per "+targetPOTstring+" POT")
h_piScPrEq0_allMC_bkg.SetLineWidth(2)
h_piScPrEq0_allMC_bkg.SetLineColor(rt.kBlue)
h_piScPrEq0_allMC_bkg.SetLineStyle(rt.kDashed)

h_phScPrGr0_allMC_bkg = rt.TH1F("h_phScPrGr0_allMC_bkg","Particle Scores for Largest e- Shower in Neutrino Background",21,-20,1)
h_phScPrEq0_allMC_bkg = rt.TH1F("h_phScPrEq0_allMC_bkg","Particle Scores for Largest e- Shower in Neutrino Background",21,-20,1)
h_phScPrGr0_allMC_bkg.GetYaxis().SetTitle("events per "+targetPOTstring+" POT")
h_phScPrGr0_allMC_bkg.SetLineWidth(2)
h_phScPrGr0_allMC_bkg.SetLineColor(8)
h_phScPrEq0_allMC_bkg.GetYaxis().SetTitle("events per "+targetPOTstring+" POT")
h_phScPrEq0_allMC_bkg.SetLineWidth(2)
h_phScPrEq0_allMC_bkg.SetLineColor(8)
h_phScPrEq0_allMC_bkg.SetLineStyle(rt.kDashed)

h_elScPrGr0_allMC_bkg = rt.TH1F("h_elScPrGr0_allMC_bkg","Particle Scores for Largest e- Shower in Neutrino Background",21,-20,1)
h_elScPrEq0_allMC_bkg = rt.TH1F("h_elScPrEq0_allMC_bkg","Particle Scores for Largest e- Shower in Neutrino Background",21,-20,1)
h_elScPrGr0_allMC_bkg.GetYaxis().SetTitle("events per "+targetPOTstring+" POT")
h_elScPrGr0_allMC_bkg.SetLineWidth(2)
h_elScPrGr0_allMC_bkg.SetLineColor(rt.kRed)
h_elScPrEq0_allMC_bkg.GetYaxis().SetTitle("events per "+targetPOTstring+" POT")
h_elScPrEq0_allMC_bkg.SetLineWidth(2)
h_elScPrEq0_allMC_bkg.SetLineColor(rt.kRed)
h_elScPrEq0_allMC_bkg.SetLineStyle(rt.kDashed)

h_piScPrGr0_CCnumu_bkg = rt.TH1F("h_piScPrGr0_CCnumu_bkg","Pion Score for Largest e- Shower in CCnumu Background",21,-20,1)
h_piScPrEq0_CCnumu_bkg = rt.TH1F("h_piScPrEq0_CCnumu_bkg","Pion Score for Largest e- Shower in CCnumu Background",21,-20,1)
h_piScPrGr0_CCnumu_bkg.GetYaxis().SetTitle("events per "+targetPOTstring+" POT")
h_piScPrGr0_CCnumu_bkg.SetLineWidth(2)
h_piScPrGr0_CCnumu_bkg.SetLineColor(rt.kBlue)
h_piScPrEq0_CCnumu_bkg.GetYaxis().SetTitle("events per "+targetPOTstring+" POT")
h_piScPrEq0_CCnumu_bkg.SetLineWidth(2)
h_piScPrEq0_CCnumu_bkg.SetLineColor(rt.kBlue)
h_piScPrEq0_CCnumu_bkg.SetLineStyle(rt.kDashed)

h_phScPrGr0_CCnumu_bkg = rt.TH1F("h_phScPrGr0_CCnumu_bkg","Photon Score for Largest e- Shower in CCnumu Background",21,-20,1)
h_phScPrEq0_CCnumu_bkg = rt.TH1F("h_phScPrEq0_CCnumu_bkg","Photon Score for Largest e- Shower in CCnumu Background",21,-20,1)
h_phScPrGr0_CCnumu_bkg.GetYaxis().SetTitle("events per "+targetPOTstring+" POT")
h_phScPrGr0_CCnumu_bkg.SetLineWidth(2)
h_phScPrGr0_CCnumu_bkg.SetLineColor(rt.kBlue)
h_phScPrEq0_CCnumu_bkg.GetYaxis().SetTitle("events per "+targetPOTstring+" POT")
h_phScPrEq0_CCnumu_bkg.SetLineWidth(2)
h_phScPrEq0_CCnumu_bkg.SetLineColor(rt.kBlue)
h_phScPrEq0_CCnumu_bkg.SetLineStyle(rt.kDashed)

h_elScPrGr0_CCnumu_bkg = rt.TH1F("h_elScPrGr0_CCnumu_bkg","Electron Score for Largest e- Shower in CCnumu Background",21,-20,1)
h_elScPrEq0_CCnumu_bkg = rt.TH1F("h_elScPrEq0_CCnumu_bkg","Electron Score for Largest e- Shower in CCnumu Background",21,-20,1)
h_elScPrGr0_CCnumu_bkg.GetYaxis().SetTitle("events per "+targetPOTstring+" POT")
h_elScPrGr0_CCnumu_bkg.SetLineWidth(2)
h_elScPrGr0_CCnumu_bkg.SetLineColor(rt.kBlue)
h_elScPrEq0_CCnumu_bkg.GetYaxis().SetTitle("events per "+targetPOTstring+" POT")
h_elScPrEq0_CCnumu_bkg.SetLineWidth(2)
h_elScPrEq0_CCnumu_bkg.SetLineColor(rt.kBlue)
h_elScPrEq0_CCnumu_bkg.SetLineStyle(rt.kDashed)

h_piScPrGr0_NCnumu_bkg = rt.TH1F("h_piScPrGr0_NCnumu_bkg","Pion Score for Largest e- Shower in NCnumu Background",21,-20,1)
h_piScPrEq0_NCnumu_bkg = rt.TH1F("h_piScPrEq0_NCnumu_bkg","Pion Score for Largest e- Shower in NCnumu Background",21,-20,1)
h_piScPrGr0_NCnumu_bkg.GetYaxis().SetTitle("events per "+targetPOTstring+" POT")
h_piScPrGr0_NCnumu_bkg.SetLineWidth(2)
h_piScPrGr0_NCnumu_bkg.SetLineColor(40)
h_piScPrEq0_NCnumu_bkg.GetYaxis().SetTitle("events per "+targetPOTstring+" POT")
h_piScPrEq0_NCnumu_bkg.SetLineWidth(2)
h_piScPrEq0_NCnumu_bkg.SetLineColor(40)
h_piScPrEq0_NCnumu_bkg.SetLineStyle(rt.kDashed)

h_phScPrGr0_NCnumu_bkg = rt.TH1F("h_phScPrGr0_NCnumu_bkg","Photon Score for Largest e- Shower in NCnumu Background",21,-20,1)
h_phScPrEq0_NCnumu_bkg = rt.TH1F("h_phScPrEq0_NCnumu_bkg","Photon Score for Largest e- Shower in NCnumu Background",21,-20,1)
h_phScPrGr0_NCnumu_bkg.GetYaxis().SetTitle("events per "+targetPOTstring+" POT")
h_phScPrGr0_NCnumu_bkg.SetLineWidth(2)
h_phScPrGr0_NCnumu_bkg.SetLineColor(40)
h_phScPrEq0_NCnumu_bkg.GetYaxis().SetTitle("events per "+targetPOTstring+" POT")
h_phScPrEq0_NCnumu_bkg.SetLineWidth(2)
h_phScPrEq0_NCnumu_bkg.SetLineColor(40)
h_phScPrEq0_NCnumu_bkg.SetLineStyle(rt.kDashed)

h_elScPrGr0_NCnumu_bkg = rt.TH1F("h_elScPrGr0_NCnumu_bkg","Electron Score for Largest e- Shower in NCnumu Background",21,-20,1)
h_elScPrEq0_NCnumu_bkg = rt.TH1F("h_elScPrEq0_NCnumu_bkg","Electron Score for Largest e- Shower in NCnumu Background",21,-20,1)
h_elScPrGr0_NCnumu_bkg.GetYaxis().SetTitle("events per "+targetPOTstring+" POT")
h_elScPrGr0_NCnumu_bkg.SetLineWidth(2)
h_elScPrGr0_NCnumu_bkg.SetLineColor(40)
h_elScPrEq0_NCnumu_bkg.GetYaxis().SetTitle("events per "+targetPOTstring+" POT")
h_elScPrEq0_NCnumu_bkg.SetLineWidth(2)
h_elScPrEq0_NCnumu_bkg.SetLineColor(40)
h_elScPrEq0_NCnumu_bkg.SetLineStyle(rt.kDashed)

h_piScPrGr0_NCnue_bkg = rt.TH1F("h_piScPrGr0_NCnue_bkg","Pion Score for Largest e- Shower in NCnue Background",21,-20,1)
h_piScPrEq0_NCnue_bkg = rt.TH1F("h_piScPrEq0_NCnue_bkg","Pion Score for Largest e- Shower in NCnue Background",21,-20,1)
h_piScPrGr0_NCnue_bkg.GetYaxis().SetTitle("events per "+targetPOTstring+" POT")
h_piScPrGr0_NCnue_bkg.SetLineWidth(2)
h_piScPrGr0_NCnue_bkg.SetLineColor(8)
h_piScPrEq0_NCnue_bkg.GetYaxis().SetTitle("events per "+targetPOTstring+" POT")
h_piScPrEq0_NCnue_bkg.SetLineWidth(2)
h_piScPrEq0_NCnue_bkg.SetLineColor(8)
h_piScPrEq0_NCnue_bkg.SetLineStyle(rt.kDashed)

h_phScPrGr0_NCnue_bkg = rt.TH1F("h_phScPrGr0_NCnue_bkg","Photon Score for Largest e- Shower in NCnue Background",21,-20,1)
h_phScPrEq0_NCnue_bkg = rt.TH1F("h_phScPrEq0_NCnue_bkg","Photon Score for Largest e- Shower in NCnue Background",21,-20,1)
h_phScPrGr0_NCnue_bkg.GetYaxis().SetTitle("events per "+targetPOTstring+" POT")
h_phScPrGr0_NCnue_bkg.SetLineWidth(2)
h_phScPrGr0_NCnue_bkg.SetLineColor(8)
h_phScPrEq0_NCnue_bkg.GetYaxis().SetTitle("events per "+targetPOTstring+" POT")
h_phScPrEq0_NCnue_bkg.SetLineWidth(2)
h_phScPrEq0_NCnue_bkg.SetLineColor(8)
h_phScPrEq0_NCnue_bkg.SetLineStyle(rt.kDashed)

h_elScPrGr0_NCnue_bkg = rt.TH1F("h_elScPrGr0_NCnue_bkg","Electron Score for Largest e- Shower in NCnue Background",21,-20,1)
h_elScPrEq0_NCnue_bkg = rt.TH1F("h_elScPrEq0_NCnue_bkg","Electron Score for Largest e- Shower in NCnue Background",21,-20,1)
h_elScPrGr0_NCnue_bkg.GetYaxis().SetTitle("events per "+targetPOTstring+" POT")
h_elScPrGr0_NCnue_bkg.SetLineWidth(2)
h_elScPrGr0_NCnue_bkg.SetLineColor(8)
h_elScPrEq0_NCnue_bkg.GetYaxis().SetTitle("events per "+targetPOTstring+" POT")
h_elScPrEq0_NCnue_bkg.SetLineWidth(2)
h_elScPrEq0_NCnue_bkg.SetLineColor(8)
h_elScPrEq0_NCnue_bkg.SetLineStyle(rt.kDashed)



n_raw_CCnumu = 0
n_raw_NCnumu = 0
n_raw_CCnue = 0
n_raw_NCnue = 0
n_raw_ext = 0
n_raw_data = 0

n_runs1to3_CCnumu = 0.
n_runs1to3_CCnue = 0.
n_runs1to3_CCnumu_pass = 0.
n_runs1to3_NCnumu_pass = 0.
n_runs1to3_CCnue_pass = 0.
n_runs1to3_NCnue_pass = 0.
n_runs1to3_ext_pass = 0.
n_runs1to3_data_pass = 0.


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

  if isinf(tnu.xsecWeight):
    continue

  if args.smallFV:
    trueVtxPos = rt.TVector3(tnu.trueVtxX, tnu.trueVtxY, tnu.trueVtxZ)
    if not isFiducial(trueVtxPos):
      continue

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

  vtxPos = rt.TVector3(tnu.vtxX, tnu.vtxY, tnu.vtxZ)
  if args.smallFV:
    vtxIsFiducial = isFiducial(vtxPos)
  else:
    vtxIsFiducial = isFiducialWC(vtxPos)
  if args.oldVertexVar:
    if tnu.nVertices < 1 or not vtxIsFiducial: #tnu.vtxIsFiducial != 1:
      continue
  else:
    if tnu.foundVertex == 0 or not vtxIsFiducial: #tnu.vtxIsFiducial != 1:
      continue

  h_cosFrac_CCnumu, h_cosFrac_NCnumu, h_cosFrac_NCnue = FillNuHistos(h_cosFrac_CCnumu,
    h_cosFrac_NCnumu, h_cosFrac_NCnue, tnu.vtxFracHitsOnCosmic, tnu.xsecWeight, eventType)

  if tnu.vtxFracHitsOnCosmic >= args.vertexFracOnCosCut:
    continue

  if args.makeKPplots:
    nearestCosKP = 999.
    nearestTEndKP = 999.
    for iKP in range(tnu.nKeypoints):
      kpPos = rt.TVector3(tnu.kpMaxPosX[iKP], tnu.kpMaxPosY[iKP], tnu.kpMaxPosZ[iKP])
      distToVtx = getDistance(kpPos, vtxPos)
      if tnu.kpFilterType[iKP] == 1 and distToVtx < nearestCosKP:
        nearestCosKP = distToVtx
      if tnu.kpFilterType[iKP] == 0 and tnu.kpClusterType[iKP] in [2,4,5,7] and distToVtx < nearestTEndKP:
        nearestTEndKP = distToVtx
  
    h_cosKPDist_CCnumu, h_cosKPDist_NCnumu, h_cosKPDist_NCnue = FillNuHistos(h_cosKPDist_CCnumu,
      h_cosKPDist_NCnumu, h_cosKPDist_NCnue, nearestCosKP, tnu.xsecWeight, eventType)
    h_tEndKPDist_CCnumu, h_tEndKPDist_NCnumu, h_tEndKPDist_NCnue = FillNuHistos(h_tEndKPDist_CCnumu,
      h_tEndKPDist_NCnumu, h_tEndKPDist_NCnue, nearestTEndKP, tnu.xsecWeight, eventType)

  visE = 0.
  nMuons = 0
  nElectrons = 0
  maxElConf = -1
  nCompMuons = 0
  elMaxComp = -1.
  elMaxPur = -1.
  elMaxQComp = -1.
  elMaxQPur = -1.
  elMaxQVtxDist = -1.
  elMaxQCosTheta = -1.
  elMaxQElScore = -1.
  elMaxQPhScore = -1.
  elMaxQPiScore = -1.
  elMaxQMuScore = -1.
  elMaxQElPurity = -1.
  elMaxQPhPurity = -1.
  elMaxQPiPurity = -1.
  elMaxQPuritySum = -1.
  elMaxQConf = -1.
  elMaxQFrac = -1.
  elMaxQ = -1.
  muMaxQ = -1.
  muMaxComp = -1.
  muMaxQComp = -1.
  muMaxQCosTheta = -1.
  muMaxQFrac = -1.

  for iT in range(tnu.nTracks):
    visE += tnu.trackCharge[iT]
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
    visE += tnu.showerCharge[iS]
    if tnu.showerIsSecondary[iS] == 1 or tnu.showerClassified[iS] == 0 or tnu.showerDistToVtx[iS] > args.distCut or tnu.showerComp[iS] < args.compCut or tnu.showerPurity[iS] < args.purityCut:
      continue
    if nMuons == 0 and tnu.showerPID[iS] == 11:
      nElectrons += 1
      elConf = tnu.showerElScore[iS] - (tnu.showerPhScore[iS] + tnu.showerPiScore[iS])/2.
      if elConf > maxElConf:
        maxElConf = elConf
      if tnu.showerCharge[iS] > elMaxQ:
        elMaxQ = tnu.showerCharge[iS]
        elMaxQComp = tnu.showerComp[iS]
        elMaxQPur = tnu.showerPurity[iS]
        elMaxQVtxDist = tnu.showerDistToVtx[iS]
        elMaxQCosTheta = tnu.showerCosTheta[iS]
        elMaxQElScore = tnu.showerElScore[iS]
        elMaxQPhScore = tnu.showerPhScore[iS]
        elMaxQPiScore = tnu.showerPiScore[iS]
        elMaxQMuScore = tnu.showerMuScore[iS]
        elMaxQElPurity = tnu.showerTrueElPurity[iS]
        elMaxQPhPurity = tnu.showerTruePhPurity[iS]
        elMaxQPiPurity = tnu.showerTruePiPurity[iS]
        elMaxQPuritySum = tnu.showerTrueElPurity[iS] + tnu.showerTruePhPurity[iS] + tnu.showerTruePiPurity[iS] + tnu.showerTrueMuPurity[iS] + tnu.showerTruePrPurity[iS]
        elMaxQConf = elConf
      if tnu.showerChargeFrac[iS] > elMaxQFrac:
        elMaxQFrac = tnu.showerChargeFrac[iS]
      if tnu.showerComp[iS] > elMaxComp:
        elMaxComp = tnu.showerComp[iS]
      if tnu.showerPurity[iS] > elMaxPur:
        elMaxPur = tnu.showerPurity[iS]

  if nMuons == 0:
    h_nEl_CCnumu, h_nEl_NCnumu, h_nEl_NCnue = FillNuHistos(h_nEl_CCnumu,
      h_nEl_NCnumu, h_nEl_NCnue, nElectrons, tnu.xsecWeight, eventType)
  h_nMu_CCnumu, h_nMu_NCnumu, h_nMu_NCnue = FillNuHistos(h_nMu_CCnumu,
    h_nMu_NCnumu, h_nMu_NCnue, nMuons, tnu.xsecWeight, eventType)
  h_nCompMu_CCnumu, h_nCompMu_NCnumu, h_nCompMu_NCnue = FillNuHistos(h_nCompMu_CCnumu,
    h_nCompMu_NCnumu, h_nCompMu_NCnue, nCompMuons, tnu.xsecWeight, eventType)

  if nElectrons >= 1:
    h_maxElConf_CCnumu, h_maxElConf_NCnumu, h_maxElConf_NCnue = FillNuHistos(h_maxElConf_CCnumu,
      h_maxElConf_NCnumu, h_maxElConf_NCnue, maxElConf, tnu.xsecWeight, eventType)
    h_elMaxComp_CCnumu, h_elMaxComp_NCnumu, h_elMaxComp_NCnue = FillNuHistos(h_elMaxComp_CCnumu,
      h_elMaxComp_NCnumu, h_elMaxComp_NCnue, elMaxComp, tnu.xsecWeight, eventType)
    h_elMaxPur_CCnumu, h_elMaxPur_NCnumu, h_elMaxPur_NCnue = FillNuHistos(h_elMaxPur_CCnumu,
      h_elMaxPur_NCnumu, h_elMaxPur_NCnue, elMaxPur, tnu.xsecWeight, eventType)
    h_elMaxQComp_CCnumu, h_elMaxQComp_NCnumu, h_elMaxQComp_NCnue = FillNuHistos(h_elMaxQComp_CCnumu,
      h_elMaxQComp_NCnumu, h_elMaxQComp_NCnue, elMaxQComp, tnu.xsecWeight, eventType)
    h_elMaxQPur_CCnumu, h_elMaxQPur_NCnumu, h_elMaxQPur_NCnue = FillNuHistos(h_elMaxQPur_CCnumu,
      h_elMaxQPur_NCnumu, h_elMaxQPur_NCnue, elMaxQPur, tnu.xsecWeight, eventType)
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

    #---- the winner ------------------------------------------
    if elMaxQConf > args.confCut:
      h_elMaxQComp_wConfCut_CCnumu, h_elMaxQComp_wConfCut_NCnumu, h_elMaxQComp_wConfCut_NCnue = FillNuHistos(h_elMaxQComp_wConfCut_CCnumu,
        h_elMaxQComp_wConfCut_NCnumu, h_elMaxQComp_wConfCut_NCnue, elMaxQComp, tnu.xsecWeight, eventType)
      h_elMaxQVtxDist_wConfCut_CCnumu, h_elMaxQVtxDist_wConfCut_NCnumu, h_elMaxQVtxDist_wConfCut_NCnue = FillNuHistos(h_elMaxQVtxDist_wConfCut_CCnumu,
        h_elMaxQVtxDist_wConfCut_NCnumu, h_elMaxQVtxDist_wConfCut_NCnue, elMaxQVtxDist, tnu.xsecWeight, eventType)
      h_elMaxQCosTheta_wConfCut_CCnumu, h_elMaxQCosTheta_wConfCut_NCnumu, h_elMaxQCosTheta_wConfCut_NCnue = FillNuHistos(h_elMaxQCosTheta_wConfCut_CCnumu,
        h_elMaxQCosTheta_wConfCut_NCnumu, h_elMaxQCosTheta_wConfCut_NCnue, elMaxQCosTheta, tnu.xsecWeight, eventType)
      h_elMaxQF_wConfCut_CCnumu, h_elMaxQF_wConfCut_NCnumu, h_elMaxQF_wConfCut_NCnue = FillNuHistos(h_elMaxQF_wConfCut_CCnumu,
        h_elMaxQF_wConfCut_NCnumu, h_elMaxQF_wConfCut_NCnue, elMaxQFrac, tnu.xsecWeight, eventType)
      h_elMaxQ_wConfCut_CCnumu, h_elMaxQ_wConfCut_NCnumu, h_elMaxQ_wConfCut_NCnue = FillNuHistos(h_elMaxQ_wConfCut_CCnumu,
        h_elMaxQ_wConfCut_NCnumu, h_elMaxQ_wConfCut_NCnue, elMaxQ, tnu.xsecWeight, eventType)
      if args.makeKPplots:
        h_cosKPDist_wConfCut_CCnumu, h_cosKPDist_wConfCut_NCnumu, h_cosKPDist_wConfCut_NCnue = FillNuHistos(h_cosKPDist_wConfCut_CCnumu,
          h_cosKPDist_wConfCut_NCnumu, h_cosKPDist_wConfCut_NCnue, nearestCosKP, tnu.xsecWeight, eventType)
        h_tEndKPDist_wConfCut_CCnumu, h_tEndKPDist_wConfCut_NCnumu, h_tEndKPDist_wConfCut_NCnue = FillNuHistos(h_tEndKPDist_wConfCut_CCnumu,
          h_tEndKPDist_wConfCut_NCnumu, h_tEndKPDist_wConfCut_NCnue, nearestTEndKP, tnu.xsecWeight, eventType)

      if elMaxQ > args.chargeCut and elMaxQFrac > args.chargeFracCut and elMaxQCosTheta > args.cosThetaCut:
        if eventType == 0:
          n_runs1to3_CCnumu_pass += tnu.xsecWeight
          h_nuE_CCnumu_wCuts.Fill(tnu.trueNuE, tnu.xsecWeight)
          h_nuEr_CCnumu_wCuts.Fill(tnu.recoNuE/1000., tnu.xsecWeight)
          if args.recoEOverflow and tnu.recoNuE/1000. > 2.6:
            h_visE_CCnumu_wCuts.Fill(2.7, tnu.xsecWeight)
          else:
            h_visE_CCnumu_wCuts.Fill(tnu.recoNuE/1000., tnu.xsecWeight)
        if eventType == 1:
          n_runs1to3_NCnumu_pass += tnu.xsecWeight
          h_nuE_NCnumu_wCuts.Fill(tnu.trueNuE, tnu.xsecWeight)
          h_nuEr_NCnumu_wCuts.Fill(tnu.recoNuE/1000., tnu.xsecWeight)
          if args.recoEOverflow and tnu.recoNuE/1000. > 2.6:
            h_visE_NCnumu_wCuts.Fill(2.7, tnu.xsecWeight)
          else:
            h_visE_NCnumu_wCuts.Fill(tnu.recoNuE/1000., tnu.xsecWeight)
        if eventType == 2:
          n_runs1to3_NCnue_pass += tnu.xsecWeight
          h_nuE_NCnue_wCuts.Fill(tnu.trueNuE, tnu.xsecWeight)
          h_nuEr_NCnue_wCuts.Fill(tnu.recoNuE/1000., tnu.xsecWeight)
          if args.recoEOverflow and tnu.recoNuE/1000. > 2.6:
            h_visE_NCnue_wCuts.Fill(2.7, tnu.xsecWeight)
          else:
            h_visE_NCnue_wCuts.Fill(tnu.recoNuE/1000., tnu.xsecWeight)

    else:
      h_phScVsPrSum_allMC_bkg.Fill(elMaxQPuritySum, elMaxQPhScore, tnu.xsecWeight)
      if elMaxQPuritySum < 1e-3:
        h_phScPrSumEq0_allMC_bkg.Fill(elMaxQPhScore, tnu.xsecWeight)
      elif elMaxQPhPurity < 1e-3:
        h_phScPrSumGr0_allMC_bkg.Fill(elMaxQPhScore, tnu.xsecWeight)
      if elMaxQPiPurity < 1e-3:
        h_piScPrEq0_allMC_bkg.Fill(elMaxQPiScore, tnu.xsecWeight)
      else:
        h_piScPrGr0_allMC_bkg.Fill(elMaxQPiScore, tnu.xsecWeight)
      if elMaxQPhPurity < 1e-3:
        h_phScPrEq0_allMC_bkg.Fill(elMaxQPhScore, tnu.xsecWeight)
      else:
        h_phScPrGr0_allMC_bkg.Fill(elMaxQPhScore, tnu.xsecWeight)
      if elMaxQElPurity < 1e-3:
        h_elScPrEq0_allMC_bkg.Fill(elMaxQElScore, tnu.xsecWeight)
        h_piScElPrEq0_allMC_bkg.Fill(elMaxQPiScore, tnu.xsecWeight)
        h_muScElPrEq0_allMC_bkg.Fill(elMaxQMuScore, tnu.xsecWeight)
      else:
        h_elScPrGr0_allMC_bkg.Fill(elMaxQElScore, tnu.xsecWeight)
        h_piScElPrGr0_allMC_bkg.Fill(elMaxQPiScore, tnu.xsecWeight)
        h_muScElPrGr0_allMC_bkg.Fill(elMaxQMuScore, tnu.xsecWeight)
      if eventType == 0:
        h_piScVsPr_CCnumu_bkg.Fill(elMaxQPiPurity, elMaxQPiScore, tnu.xsecWeight)
        h_phScVsPr_CCnumu_bkg.Fill(elMaxQPhPurity, elMaxQPhScore, tnu.xsecWeight)
        h_elScVsPr_CCnumu_bkg.Fill(elMaxQElPurity, elMaxQElScore, tnu.xsecWeight)
        if elMaxQPiPurity < 1e-3:
          h_piScPrEq0_CCnumu_bkg.Fill(elMaxQPiScore, tnu.xsecWeight)
        else:
          h_piScPrGr0_CCnumu_bkg.Fill(elMaxQPiScore, tnu.xsecWeight)
        if elMaxQPhPurity < 1e-3:
          h_phScPrEq0_CCnumu_bkg.Fill(elMaxQPhScore, tnu.xsecWeight)
        else:
          h_phScPrGr0_CCnumu_bkg.Fill(elMaxQPhScore, tnu.xsecWeight)
        if elMaxQElPurity < 1e-3:
          h_elScPrEq0_CCnumu_bkg.Fill(elMaxQElScore, tnu.xsecWeight)
        else:
          h_elScPrGr0_CCnumu_bkg.Fill(elMaxQElScore, tnu.xsecWeight)
      if eventType == 1:
        h_piScVsPr_NCnumu_bkg.Fill(elMaxQPiPurity, elMaxQPiScore, tnu.xsecWeight)
        h_phScVsPr_NCnumu_bkg.Fill(elMaxQPhPurity, elMaxQPhScore, tnu.xsecWeight)
        h_elScVsPr_NCnumu_bkg.Fill(elMaxQElPurity, elMaxQElScore, tnu.xsecWeight)
        if elMaxQPiPurity < 1e-3:
          h_piScPrEq0_NCnumu_bkg.Fill(elMaxQPiScore, tnu.xsecWeight)
        else:
          h_piScPrGr0_NCnumu_bkg.Fill(elMaxQPiScore, tnu.xsecWeight)
        if elMaxQPhPurity < 1e-3:
          h_phScPrEq0_NCnumu_bkg.Fill(elMaxQPhScore, tnu.xsecWeight)
        else:
          h_phScPrGr0_NCnumu_bkg.Fill(elMaxQPhScore, tnu.xsecWeight)
        if elMaxQElPurity < 1e-3:
          h_elScPrEq0_NCnumu_bkg.Fill(elMaxQElScore, tnu.xsecWeight)
        else:
          h_elScPrGr0_NCnumu_bkg.Fill(elMaxQElScore, tnu.xsecWeight)
      if eventType == 2:
        h_piScVsPr_NCnue_bkg.Fill(elMaxQPiPurity, elMaxQPiScore, tnu.xsecWeight)
        h_phScVsPr_NCnue_bkg.Fill(elMaxQPhPurity, elMaxQPhScore, tnu.xsecWeight)
        h_elScVsPr_NCnue_bkg.Fill(elMaxQElPurity, elMaxQElScore, tnu.xsecWeight)
        if elMaxQPiPurity < 1e-3:
          h_piScPrEq0_NCnue_bkg.Fill(elMaxQPiScore, tnu.xsecWeight)
        else:
          h_piScPrGr0_NCnue_bkg.Fill(elMaxQPiScore, tnu.xsecWeight)
        if elMaxQPhPurity < 1e-3:
          h_phScPrEq0_NCnue_bkg.Fill(elMaxQPhScore, tnu.xsecWeight)
        else:
          h_phScPrGr0_NCnue_bkg.Fill(elMaxQPhScore, tnu.xsecWeight)
        if elMaxQElPurity < 1e-3:
          h_elScPrEq0_NCnue_bkg.Fill(elMaxQElScore, tnu.xsecWeight)
        else:
          h_elScPrGr0_NCnue_bkg.Fill(elMaxQElScore, tnu.xsecWeight)
    #--------------------------------------------------------------
    #if eventType == 0 and elMaxQElScore > -1. and (elMaxQPhScore > -1. or elMaxQPiScore > -1.):
    ##elif eventType == 0:
    #  if elMaxQConf > args.confCut:
    #    print("THIS IS SIGNAL!!!!!")
    #  print("CC numu background (fileid, run, subrun, event): (%i, %i, %i, %i) with (e-,ph,pi) scores: (%f, %f, %f) removed"%(tnu.fileid,tnu.run,tnu.subrun,tnu.event,elMaxQElScore,elMaxQPhScore,elMaxQPiScore))

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

  if abs(tnue.trueNuPDG) != 12 or tnue.trueNuCCNC != 0 or isinf(tnue.xsecWeight):
    continue

  if args.smallFV:
    trueVtxPos = rt.TVector3(tnue.trueVtxX, tnue.trueVtxY, tnue.trueVtxZ)
    if not isFiducial(trueVtxPos):
      continue

  n_raw_CCnue += 1
  n_runs1to3_CCnue += tnue.xsecWeight
  h_nuE_CCnue_nCuts.Fill(tnue.trueNuE, tnue.xsecWeight)
  h_nuEr_CCnue_nCuts.Fill(tnue.recoNuE/1000., tnue.xsecWeight)

  vtxPos = rt.TVector3(tnue.vtxX, tnue.vtxY, tnue.vtxZ)
  if args.smallFV:
    vtxIsFiducial = isFiducial(vtxPos)
  else:
    vtxIsFiducial = isFiducialWC(vtxPos)
  if args.oldVertexVar:
    if tnue.nVertices < 1 or not vtxIsFiducial: #tnue.vtxIsFiducial != 1:
      continue
  else:
    if tnue.foundVertex == 0 or not vtxIsFiducial: #tnue.vtxIsFiducial != 1:
      continue

  h_cosFrac_CCnue.Fill(tnue.vtxFracHitsOnCosmic, tnue.xsecWeight)
  #if tnue.vtxFracHitsOnCosmic < 0 or tnue.vtxFracHitsOnCosmic > 1.:
  if tnue.vtxFracHitsOnCosmic > 1.:
    print(tnue.vtxFracHitsOnCosmic)

  if tnue.vtxFracHitsOnCosmic >= args.vertexFracOnCosCut:
    continue

  if args.makeKPplots:
    nearestCosKP = 999.
    nearestTEndKP = 999.
    for iKP in range(tnue.nKeypoints):
      kpPos = rt.TVector3(tnue.kpMaxPosX[iKP], tnue.kpMaxPosY[iKP], tnue.kpMaxPosZ[iKP])
      distToVtx = getDistance(kpPos, vtxPos)
      if tnue.kpFilterType[iKP] == 1 and distToVtx < nearestCosKP:
        nearestCosKP = distToVtx
      if tnue.kpFilterType[iKP] == 0 and tnue.kpClusterType[iKP] in [2,4,5,7] and distToVtx < nearestTEndKP:
        nearestTEndKP = distToVtx
  
    h_cosKPDist_CCnue.Fill(nearestCosKP, tnue.xsecWeight)
    h_tEndKPDist_CCnue.Fill(nearestTEndKP, tnue.xsecWeight)

  visE = 0.
  nMuons = 0
  nElectrons = 0
  maxElConf = -1
  nCompMuons = 0
  elMaxComp = -1.
  elMaxPur = -1.
  elMaxQComp = -1.
  elMaxQPur = -1.
  elMaxQVtxDist = -1.
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
    visE += tnue.trackCharge[iT]
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
    visE += tnue.showerCharge[iS]
    if tnue.showerIsSecondary[iS] == 1 or tnue.showerClassified[iS] == 0 or tnue.showerDistToVtx[iS] > args.distCut or tnue.showerComp[iS] < args.compCut or tnue.showerPurity[iS] < args.purityCut:
      continue
    if nMuons == 0 and tnue.showerPID[iS] == 11:
      nElectrons += 1
      elConf = tnue.showerElScore[iS] - (tnue.showerPhScore[iS] + tnue.showerPiScore[iS])/2.
      if elConf > maxElConf:
        maxElConf = elConf
      if tnue.showerCharge[iS] > elMaxQ:
        elMaxQ = tnue.showerCharge[iS]
        elMaxQComp = tnue.showerComp[iS]
        elMaxQPur = tnue.showerPurity[iS]
        elMaxQVtxDist = tnue.showerDistToVtx[iS]
        elMaxQCosTheta = tnue.showerCosTheta[iS]
        elMaxQElScore = tnue.showerElScore[iS]
        elMaxQPhScore = tnue.showerPhScore[iS]
        elMaxQPiScore = tnue.showerPiScore[iS]
        elMaxQConf = elConf
      if tnue.showerChargeFrac[iS] > elMaxQFrac:
        elMaxQFrac = tnue.showerChargeFrac[iS]
      if tnue.showerComp[iS] > elMaxComp:
        elMaxComp = tnue.showerComp[iS]
      if tnue.showerPurity[iS] > elMaxPur:
        elMaxPur = tnue.showerPurity[iS]

  if nMuons == 0:
    h_nEl_CCnue.Fill(nElectrons, tnue.xsecWeight)
  h_nMu_CCnue.Fill(nMuons, tnue.xsecWeight)
  h_nCompMu_CCnue.Fill(nCompMuons, tnue.xsecWeight)

  if nElectrons >= 1:
    h_maxElConf_CCnue.Fill(maxElConf, tnue.xsecWeight)
    h_elMaxComp_CCnue.Fill(elMaxComp, tnue.xsecWeight)
    h_elMaxPur_CCnue.Fill(elMaxPur, tnue.xsecWeight)
    h_elMaxQComp_CCnue.Fill(elMaxQComp, tnue.xsecWeight)
    h_elMaxQPur_CCnue.Fill(elMaxQPur, tnue.xsecWeight)
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
      h_elMaxQVtxDist_wConfCut_CCnue.Fill(elMaxQVtxDist, tnue.xsecWeight)
      h_elMaxQCosTheta_wConfCut_CCnue.Fill(elMaxQCosTheta, tnue.xsecWeight)
      h_elMaxQF_wConfCut_CCnue.Fill(elMaxQFrac, tnue.xsecWeight)
      h_elMaxQ_wConfCut_CCnue.Fill(elMaxQ, tnue.xsecWeight)
      if args.makeKPplots:
        h_cosKPDist_wConfCut_CCnue.Fill(nearestCosKP, tnue.xsecWeight)
        h_tEndKPDist_wConfCut_CCnue.Fill(nearestTEndKP, tnue.xsecWeight)
      if elMaxQ > args.chargeCut and elMaxQFrac > args.chargeFracCut and elMaxQCosTheta > args.cosThetaCut:
        n_runs1to3_CCnue_pass += tnue.xsecWeight
        h_nuE_CCnue_wCuts.Fill(tnue.trueNuE, tnue.xsecWeight)
        h_nuEr_CCnue_wCuts.Fill(tnue.recoNuE/1000., tnue.xsecWeight)
        if args.recoEOverflow and tnue.recoNuE/1000. > 2.6:
          h_visE_CCnue_wCuts.Fill(2.7, tnue.xsecWeight)
        else:
          h_visE_CCnue_wCuts.Fill(tnue.recoNuE/1000., tnue.xsecWeight)

  if nMuons >= 1:
    h_muMaxComp_CCnue.Fill(muMaxComp, tnue.xsecWeight)
    h_muMaxQComp_CCnue.Fill(muMaxQComp, tnue.xsecWeight)
    h_muMaxQCosTheta_CCnue.Fill(muMaxQCosTheta, tnue.xsecWeight)
    h_muMaxQF_CCnue.Fill(muMaxQFrac, tnue.xsecWeight)
    h_muMaxQ_CCnue.Fill(muMaxQ, tnue.xsecWeight)





for i in range(text.GetEntries()):

  text.GetEntry(i)

  n_raw_ext += 1

  vtxPos = rt.TVector3(text.vtxX, text.vtxY, text.vtxZ)
  if args.smallFV:
    vtxIsFiducial = isFiducial(vtxPos)
  else:
    vtxIsFiducial = isFiducialWC(vtxPos)
  if args.oldVertexVar or args.NPMLCosmicConfig:
    if text.nVertices < 1 or not vtxIsFiducial: #text.vtxIsFiducial != 1:
      continue
  else:
    if text.foundVertex == 0 or not vtxIsFiducial: #text.vtxIsFiducial != 1:
      continue
  #if text.foundVertex == 0 or not vtxIsFiducial: #text.vtxIsFiducial != 1:
  #  continue

  h_cosFrac_ext.Fill(text.vtxFracHitsOnCosmic)

  if text.vtxFracHitsOnCosmic >= args.vertexFracOnCosCut:
    continue

  if args.makeKPplots:
    nearestCosKP = 999.
    nearestTEndKP = 999.
    for iKP in range(text.nKeypoints):
      kpPos = rt.TVector3(text.kpMaxPosX[iKP], text.kpMaxPosY[iKP], text.kpMaxPosZ[iKP])
      distToVtx = getDistance(kpPos, vtxPos)
      if text.kpFilterType[iKP] == 1 and distToVtx < nearestCosKP:
        nearestCosKP = distToVtx
      if text.kpFilterType[iKP] == 0 and text.kpClusterType[iKP] in [2,4,5,7] and distToVtx < nearestTEndKP:
        nearestTEndKP = distToVtx
  
    h_cosKPDist_ext.Fill(nearestCosKP)
    h_tEndKPDist_ext.Fill(nearestTEndKP)

  visE = 0.
  nMuons = 0
  nElectrons = 0
  maxElConf = -1
  nCompMuons = 0
  elMaxComp = -1.
  elMaxPur = -1.
  elMaxQComp = -1.
  elMaxQPur = -1.
  elMaxQVtxDist = -1.
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
    visE += text.trackCharge[iT]
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
    visE += text.showerCharge[iS]
    if text.showerIsSecondary[iS] == 1 or text.showerClassified[iS] == 0 or text.showerDistToVtx[iS] > args.distCut or text.showerComp[iS] < args.compCut or text.showerPurity[iS] < args.purityCut:
      continue
    if nMuons == 0 and text.showerPID[iS] == 11:
      nElectrons += 1
      elConf = text.showerElScore[iS] - (text.showerPhScore[iS] + text.showerPiScore[iS])/2.
      if elConf > maxElConf:
        maxElConf = elConf
      if text.showerCharge[iS] > elMaxQ:
        elMaxQ = text.showerCharge[iS]
        elMaxQComp = text.showerComp[iS]
        elMaxQPur = text.showerPurity[iS]
        elMaxQVtxDist = text.showerDistToVtx[iS]
        elMaxQCosTheta = text.showerCosTheta[iS]
        elMaxQElScore = text.showerElScore[iS]
        elMaxQPhScore = text.showerPhScore[iS]
        elMaxQPiScore = text.showerPiScore[iS]
        elMaxQConf = elConf
      if text.showerChargeFrac[iS] > elMaxQFrac:
        elMaxQFrac = text.showerChargeFrac[iS]
      if text.showerComp[iS] > elMaxComp:
        elMaxComp = text.showerComp[iS]
      if text.showerPurity[iS] > elMaxPur:
        elMaxPur = text.showerPurity[iS]

  if nMuons == 0:
    h_nEl_ext.Fill(nElectrons)
  h_nMu_ext.Fill(nMuons)
  h_nCompMu_ext.Fill(nCompMuons)

  if nElectrons >= 1:
    h_maxElConf_ext.Fill(maxElConf)
    h_elMaxComp_ext.Fill(elMaxComp)
    h_elMaxPur_ext.Fill(elMaxPur)
    h_elMaxQComp_ext.Fill(elMaxQComp)
    h_elMaxQPur_ext.Fill(elMaxQPur)
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
      h_elMaxQVtxDist_wConfCut_ext.Fill(elMaxQVtxDist)
      h_elMaxQCosTheta_wConfCut_ext.Fill(elMaxQCosTheta)
      h_elMaxQF_wConfCut_ext.Fill(elMaxQFrac)
      h_elMaxQ_wConfCut_ext.Fill(elMaxQ)
      if args.makeKPplots:
        h_cosKPDist_wConfCut_ext.Fill(nearestCosKP)
        h_tEndKPDist_wConfCut_ext.Fill(nearestTEndKP)
      if elMaxQ > args.chargeCut and elMaxQFrac > args.chargeFracCut and elMaxQCosTheta > args.cosThetaCut:
        if args.printEXTinfo:
          print("extBNB background event passed selection: (fileid, run, subrun, event, vtxX, vtxY, vtxZ) = (%i, %i, %i, %i, %f, %f, %f)"%(text.fileid, text.run, text.subrun, text.event, text.vtxX, text.vtxY, text.vtxZ))
        n_runs1to3_ext_pass += 1.
        #h_nuE_ext_wCuts.Fill(text.trueNuE)
        h_nuEr_ext_wCuts.Fill(text.recoNuE/1000.)
        if args.recoEOverflow and text.recoNuE/1000. > 2.6:
          h_visE_ext_wCuts.Fill(2.7)
        else:
          h_visE_ext_wCuts.Fill(text.recoNuE/1000.)

  if nMuons >= 1:
    h_muMaxComp_ext.Fill(muMaxComp)
    h_muMaxQComp_ext.Fill(muMaxQComp)
    h_muMaxQCosTheta_ext.Fill(muMaxQCosTheta)
    h_muMaxQF_ext.Fill(muMaxQFrac)
    h_muMaxQ_ext.Fill(muMaxQ)




for i in range(tdata.GetEntries()):

  tdata.GetEntry(i)

  n_raw_data += 1

  vtxPos = rt.TVector3(tdata.vtxX, tdata.vtxY, tdata.vtxZ)
  if args.smallFV:
    vtxIsFiducial = isFiducial(vtxPos)
  else:
    vtxIsFiducial = isFiducialWC(vtxPos)
  if tdata.nVertices < 1 or not vtxIsFiducial: #tdata.vtxIsFiducial != 1:
    continue

  if tdata.vtxFracHitsOnCosmic >= args.vertexFracOnCosCut:
    continue

  visE = 0.
  nMuons = 0
  nElectrons = 0
  elMaxQCosTheta = -1.
  elMaxQConf = -1.
  elMaxQFrac = -1.
  elMaxQ = -1.
  for iT in range(tdata.nTracks):
    visE += tdata.trackCharge[iT]
    if tdata.trackIsSecondary[iT] == 1:
      continue
    if tdata.trackClassified[iT] == 1 and tdata.trackPID[iT] == 13:
      nMuons += 1
  for iS in range(tdata.nShowers):
    visE += tdata.showerCharge[iS]
    if tdata.showerIsSecondary[iS] == 1 or tdata.showerClassified[iS] == 0 or tdata.showerDistToVtx[iS] > args.distCut or tdata.showerComp[iS] < args.compCut or tdata.showerPurity[iS] < args.purityCut:
      continue
    if nMuons == 0 and tdata.showerPID[iS] == 11:
      nElectrons += 1
      elConf = tdata.showerElScore[iS] - (tdata.showerPhScore[iS] + tdata.showerPiScore[iS])/2.
      if tdata.showerCharge[iS] > elMaxQ:
        elMaxQ = tdata.showerCharge[iS]
        elMaxQCosTheta = tdata.showerCosTheta[iS]
        elMaxQConf = elConf
      if tdata.showerChargeFrac[iS] > elMaxQFrac:
        elMaxQFrac = tdata.showerChargeFrac[iS]

  if nElectrons >= 1:
    if elMaxQConf > args.confCut:
      if elMaxQ > args.chargeCut and elMaxQFrac > args.chargeFracCut and elMaxQCosTheta > args.cosThetaCut:
        n_runs1to3_data_pass += 1.
        if args.recoEOverflow and tdata.recoNuE/1000. > 2.6:
          h_visE_data_wCuts.Fill(2.7)
        else:
          h_visE_data_wCuts.Fill(tdata.recoNuE/1000.)
        out_data_sel_evts.write("%i %i %i %i %f\n"%(tdata.fileid, tdata.run, tdata.subrun, tdata.event, tdata.recoNuE))



out_data_sel_evts.close()

n_runs1to3_CCnumu *= (targetPOT/tnuPOTsum)
n_runs1to3_CCnue *= (targetPOT/tnuePOTsum)
n_runs1to3_CCnumu_pass *= targetPOT/tnuPOTsum
n_runs1to3_NCnumu_pass *= targetPOT/tnuPOTsum
n_runs1to3_CCnue_pass *= targetPOT/tnuePOTsum
n_runs1to3_NCnue_pass *= targetPOT/tnuPOTsum
n_runs1to3_ext_pass *= targetPOT/textPOTsum
n_predicted_tot_pass = n_runs1to3_CCnue_pass + n_runs1to3_NCnue_pass + n_runs1to3_CCnumu_pass + n_runs1to3_NCnumu_pass + n_runs1to3_ext_pass

print()
print("n_raw_CCnumu: ", n_raw_CCnumu)
print("n_raw_NCnumu: ", n_raw_NCnumu)
print("n_raw_CCnue: ", n_raw_CCnue)
print("n_raw_NCnue: ", n_raw_NCnue)
print("n_raw_ext: ", n_raw_ext)
print("n_raw_data: ", n_raw_data)
print()
print("n_runs1to3_CCnumu: ", n_runs1to3_CCnumu)
print("n_runs1to3_CCnue: ", n_runs1to3_CCnue)
print()
print("n_runs1to3_CCnue_pass: ", n_runs1to3_CCnue_pass)
print("n_runs1to3_NCnue_pass: ", n_runs1to3_NCnue_pass)
print("n_runs1to3_CCnumu_pass: ", n_runs1to3_CCnumu_pass)
print("n_runs1to3_NCnumu_pass: ", n_runs1to3_NCnumu_pass)
print("n_runs1to3_ext_pass: ", n_runs1to3_ext_pass)
print("n_runs1to3_data_pass: ", n_runs1to3_data_pass)
print("n_predicted_tot_pass: ", n_predicted_tot_pass)
print()
print("%s integral: %f ; scaling by %.2e/%.2e"%("h_cosFrac_CCnue",h_cosFrac_CCnue.Integral(),targetPOT,tnuePOTsum))
print("%s integral: %f ; scaling by %.2e/%.2e"%("h_cosFrac_CCnumu",h_cosFrac_CCnumu.Integral(),targetPOT,tnuPOTsum))
print("%s integral: %f ; scaling by %.2e/%.2e"%("h_cosFrac_NCnue",h_cosFrac_NCnue.Integral(),targetPOT,tnuPOTsum))
print("%s integral: %f ; scaling by %.2e/%.2e"%("h_cosFrac_NCnumu",h_cosFrac_NCnumu.Integral(),targetPOT,tnuPOTsum))
print("%s integral: %f ; scaling by %.2e/%.2e"%("h_cosFrac_ext",h_cosFrac_ext.Integral(),targetPOT,textPOTsum))
print()
print("CC nue cut efficiency: %f"%(n_runs1to3_CCnue_pass/n_runs1to3_CCnue))
print("CC nue cut purity: %f"%(n_runs1to3_CCnue_pass/n_predicted_tot_pass))
print()

#print("h_cosFrac_CCnue unscaled integral: ", h_cosFrac_CCnue.Integral())
#print("h_nEl_CCnue unscaled integral: ", h_nEl_CCnue.Integral())

h_nuE_CCnue_nCuts.Scale(targetPOT/tnuePOTsum)
h_nuE_CCnumu_wCuts.Scale(targetPOT/tnuPOTsum)
h_nuE_NCnumu_wCuts.Scale(targetPOT/tnuPOTsum)
h_nuE_CCnue_wCuts.Scale(targetPOT/tnuePOTsum)
h_nuE_NCnue_wCuts.Scale(targetPOT/tnuPOTsum)
#h_nuE_ext_wCuts.Scale(targetPOT/textPOTsum)
h_nuE_CCnue_eff.Divide(h_nuE_CCnue_wCuts,h_nuE_CCnue_nCuts,1,1,"B")
#h_nuE_CCnue_eff.Divide(h_nuE_CCnue_wCuts,h_nuE_CCnue_nCuts,1,1)
h_nuE_all_wCuts.Add(h_nuE_CCnumu_wCuts)
h_nuE_all_wCuts.Add(h_nuE_NCnumu_wCuts)
h_nuE_all_wCuts.Add(h_nuE_CCnue_wCuts)
h_nuE_all_wCuts.Add(h_nuE_NCnue_wCuts)
#h_nuE_all_wCuts.Add(h_nuE_ext_wCuts)
#h_nuE_CCnue_pur.Divide(h_nuE_CCnue_wCuts,h_nuE_all_wCuts,1,1,"B")
h_nuE_CCnue_pur.Divide(h_nuE_CCnue_wCuts,h_nuE_all_wCuts)

h_nuEr_CCnue_nCuts.Scale(targetPOT/tnuePOTsum)
h_nuEr_CCnumu_wCuts.Scale(targetPOT/tnuPOTsum)
h_nuEr_NCnumu_wCuts.Scale(targetPOT/tnuPOTsum)
h_nuEr_CCnue_wCuts.Scale(targetPOT/tnuePOTsum)
h_nuEr_NCnue_wCuts.Scale(targetPOT/tnuPOTsum)
h_nuEr_ext_wCuts.Scale(targetPOT/textPOTsum)
h_nuEr_CCnue_eff.Divide(h_nuEr_CCnue_wCuts,h_nuEr_CCnue_nCuts,1,1,"B")
#h_nuEr_CCnue_eff.Divide(h_nuEr_CCnue_wCuts,h_nuEr_CCnue_nCuts,1,1)
h_nuEr_all_wCuts.Add(h_nuEr_CCnumu_wCuts)
h_nuEr_all_wCuts.Add(h_nuEr_NCnumu_wCuts)
h_nuEr_all_wCuts.Add(h_nuEr_CCnue_wCuts)
h_nuEr_all_wCuts.Add(h_nuEr_NCnue_wCuts)
h_nuEr_all_wCuts.Add(h_nuEr_ext_wCuts)
#h_nuEr_CCnue_pur.Divide(h_nuEr_CCnue_wCuts,h_nuEr_all_wCuts,1,1,"B")
h_nuEr_CCnue_pur.Divide(h_nuEr_CCnue_wCuts,h_nuEr_all_wCuts)

h_visE_CCnumu_wCuts.Scale(targetPOT/tnuPOTsum)
h_visE_NCnumu_wCuts.Scale(targetPOT/tnuPOTsum)
h_visE_CCnue_wCuts.Scale(targetPOT/tnuePOTsum)
h_visE_NCnue_wCuts.Scale(targetPOT/tnuPOTsum)
h_visE_ext_wCuts.Scale(targetPOT/textPOTsum)

h_visE_all_wCuts.Add(h_visE_ext_wCuts)
h_visE_all_wCuts.Add(h_visE_NCnue_wCuts)
h_visE_all_wCuts.Add(h_visE_NCnumu_wCuts)
h_visE_all_wCuts.Add(h_visE_CCnumu_wCuts)
h_visE_all_wCuts.Add(h_visE_CCnue_wCuts)

h_cosFrac_CCnumu.Scale(targetPOT/tnuPOTsum)
h_cosFrac_NCnumu.Scale(targetPOT/tnuPOTsum)
h_cosFrac_CCnue.Scale(targetPOT/tnuePOTsum)
h_cosFrac_NCnue.Scale(targetPOT/tnuPOTsum)
h_cosFrac_ext.Scale(targetPOT/textPOTsum)

if args.makeKPplots:
  h_cosKPDist_CCnumu.Scale(targetPOT/tnuPOTsum)
  h_cosKPDist_NCnumu.Scale(targetPOT/tnuPOTsum)
  h_cosKPDist_CCnue.Scale(targetPOT/tnuePOTsum)
  h_cosKPDist_NCnue.Scale(targetPOT/tnuPOTsum)
  h_cosKPDist_ext.Scale(targetPOT/textPOTsum)
  
  h_tEndKPDist_CCnumu.Scale(targetPOT/tnuPOTsum)
  h_tEndKPDist_NCnumu.Scale(targetPOT/tnuPOTsum)
  h_tEndKPDist_CCnue.Scale(targetPOT/tnuePOTsum)
  h_tEndKPDist_NCnue.Scale(targetPOT/tnuPOTsum)
  h_tEndKPDist_ext.Scale(targetPOT/textPOTsum)
  
  h_cosKPDist_wConfCut_CCnumu.Scale(targetPOT/tnuPOTsum)
  h_cosKPDist_wConfCut_NCnumu.Scale(targetPOT/tnuPOTsum)
  h_cosKPDist_wConfCut_CCnue.Scale(targetPOT/tnuePOTsum)
  h_cosKPDist_wConfCut_NCnue.Scale(targetPOT/tnuPOTsum)
  h_cosKPDist_wConfCut_ext.Scale(targetPOT/textPOTsum)
  
  h_tEndKPDist_wConfCut_CCnumu.Scale(targetPOT/tnuPOTsum)
  h_tEndKPDist_wConfCut_NCnumu.Scale(targetPOT/tnuPOTsum)
  h_tEndKPDist_wConfCut_CCnue.Scale(targetPOT/tnuePOTsum)
  h_tEndKPDist_wConfCut_NCnue.Scale(targetPOT/tnuPOTsum)
  h_tEndKPDist_wConfCut_ext.Scale(targetPOT/textPOTsum)

h_nEl_CCnumu.Scale(targetPOT/tnuPOTsum)
h_nEl_NCnumu.Scale(targetPOT/tnuPOTsum)
h_nEl_CCnue.Scale(targetPOT/tnuePOTsum)
h_nEl_NCnue.Scale(targetPOT/tnuPOTsum)
h_nEl_ext.Scale(targetPOT/textPOTsum)

#print("h_cosFrac_CCnue scaled integral: ", h_cosFrac_CCnue.Integral())
#print("h_nEl_CCnue scaled integral: ", h_nEl_CCnue.Integral())

h_maxElConf_CCnumu.Scale(targetPOT/tnuPOTsum)
h_maxElConf_NCnumu.Scale(targetPOT/tnuPOTsum)
h_maxElConf_CCnue.Scale(targetPOT/tnuePOTsum)
h_maxElConf_NCnue.Scale(targetPOT/tnuPOTsum)
h_maxElConf_ext.Scale(targetPOT/textPOTsum)

h_nMu_CCnumu.Scale(targetPOT/tnuPOTsum)
h_nMu_NCnumu.Scale(targetPOT/tnuPOTsum)
h_nMu_CCnue.Scale(targetPOT/tnuePOTsum)
h_nMu_NCnue.Scale(targetPOT/tnuPOTsum)
h_nMu_ext.Scale(targetPOT/textPOTsum)

h_nCompMu_CCnumu.Scale(targetPOT/tnuPOTsum)
h_nCompMu_NCnumu.Scale(targetPOT/tnuPOTsum)
h_nCompMu_CCnue.Scale(targetPOT/tnuePOTsum)
h_nCompMu_NCnue.Scale(targetPOT/tnuPOTsum)
h_nCompMu_ext.Scale(targetPOT/textPOTsum)

h_elMaxQ_CCnumu.Scale(targetPOT/tnuPOTsum)
h_elMaxQ_NCnumu.Scale(targetPOT/tnuPOTsum)
h_elMaxQ_CCnue.Scale(targetPOT/tnuePOTsum)
h_elMaxQ_NCnue.Scale(targetPOT/tnuPOTsum)
h_elMaxQ_ext.Scale(targetPOT/textPOTsum)

h_muMaxQ_CCnumu.Scale(targetPOT/tnuPOTsum)
h_muMaxQ_NCnumu.Scale(targetPOT/tnuPOTsum)
h_muMaxQ_CCnue.Scale(targetPOT/tnuePOTsum)
h_muMaxQ_NCnue.Scale(targetPOT/tnuPOTsum)
h_muMaxQ_ext.Scale(targetPOT/textPOTsum)

h_elMaxQF_CCnumu.Scale(targetPOT/tnuPOTsum)
h_elMaxQF_NCnumu.Scale(targetPOT/tnuPOTsum)
h_elMaxQF_CCnue.Scale(targetPOT/tnuePOTsum)
h_elMaxQF_NCnue.Scale(targetPOT/tnuPOTsum)
h_elMaxQF_ext.Scale(targetPOT/textPOTsum)

h_muMaxQF_CCnumu.Scale(targetPOT/tnuPOTsum)
h_muMaxQF_NCnumu.Scale(targetPOT/tnuPOTsum)
h_muMaxQF_CCnue.Scale(targetPOT/tnuePOTsum)
h_muMaxQF_NCnue.Scale(targetPOT/tnuPOTsum)
h_muMaxQF_ext.Scale(targetPOT/textPOTsum)

h_elMaxComp_CCnumu.Scale(targetPOT/tnuPOTsum)
h_elMaxComp_NCnumu.Scale(targetPOT/tnuPOTsum)
h_elMaxComp_CCnue.Scale(targetPOT/tnuePOTsum)
h_elMaxComp_NCnue.Scale(targetPOT/tnuPOTsum)
h_elMaxComp_ext.Scale(targetPOT/textPOTsum)

h_elMaxPur_CCnumu.Scale(targetPOT/tnuPOTsum)
h_elMaxPur_NCnumu.Scale(targetPOT/tnuPOTsum)
h_elMaxPur_CCnue.Scale(targetPOT/tnuePOTsum)
h_elMaxPur_NCnue.Scale(targetPOT/tnuPOTsum)
h_elMaxPur_ext.Scale(targetPOT/textPOTsum)

h_muMaxComp_CCnumu.Scale(targetPOT/tnuPOTsum)
h_muMaxComp_NCnumu.Scale(targetPOT/tnuPOTsum)
h_muMaxComp_CCnue.Scale(targetPOT/tnuePOTsum)
h_muMaxComp_NCnue.Scale(targetPOT/tnuPOTsum)
h_muMaxComp_ext.Scale(targetPOT/textPOTsum)

h_elMaxQComp_CCnumu.Scale(targetPOT/tnuPOTsum)
h_elMaxQComp_NCnumu.Scale(targetPOT/tnuPOTsum)
h_elMaxQComp_CCnue.Scale(targetPOT/tnuePOTsum)
h_elMaxQComp_NCnue.Scale(targetPOT/tnuPOTsum)
h_elMaxQComp_ext.Scale(targetPOT/textPOTsum)

h_elMaxQPur_CCnumu.Scale(targetPOT/tnuPOTsum)
h_elMaxQPur_NCnumu.Scale(targetPOT/tnuPOTsum)
h_elMaxQPur_CCnue.Scale(targetPOT/tnuePOTsum)
h_elMaxQPur_NCnue.Scale(targetPOT/tnuPOTsum)
h_elMaxQPur_ext.Scale(targetPOT/textPOTsum)

h_muMaxQComp_CCnumu.Scale(targetPOT/tnuPOTsum)
h_muMaxQComp_NCnumu.Scale(targetPOT/tnuPOTsum)
h_muMaxQComp_CCnue.Scale(targetPOT/tnuePOTsum)
h_muMaxQComp_NCnue.Scale(targetPOT/tnuPOTsum)
h_muMaxQComp_ext.Scale(targetPOT/textPOTsum)

h_elMaxQCosTheta_CCnumu.Scale(targetPOT/tnuPOTsum)
h_elMaxQCosTheta_NCnumu.Scale(targetPOT/tnuPOTsum)
h_elMaxQCosTheta_CCnue.Scale(targetPOT/tnuePOTsum)
h_elMaxQCosTheta_NCnue.Scale(targetPOT/tnuPOTsum)
h_elMaxQCosTheta_ext.Scale(targetPOT/textPOTsum)

h_muMaxQCosTheta_CCnumu.Scale(targetPOT/tnuPOTsum)
h_muMaxQCosTheta_NCnumu.Scale(targetPOT/tnuPOTsum)
h_muMaxQCosTheta_CCnue.Scale(targetPOT/tnuePOTsum)
h_muMaxQCosTheta_NCnue.Scale(targetPOT/tnuPOTsum)
h_muMaxQCosTheta_ext.Scale(targetPOT/textPOTsum)

h_elMaxQElScore_CCnumu.Scale(targetPOT/tnuPOTsum)
h_elMaxQElScore_NCnumu.Scale(targetPOT/tnuPOTsum)
h_elMaxQElScore_CCnue.Scale(targetPOT/tnuePOTsum)
h_elMaxQElScore_NCnue.Scale(targetPOT/tnuPOTsum)
h_elMaxQElScore_ext.Scale(targetPOT/textPOTsum)

h_elMaxQPhScore_CCnumu.Scale(targetPOT/tnuPOTsum)
h_elMaxQPhScore_NCnumu.Scale(targetPOT/tnuPOTsum)
h_elMaxQPhScore_CCnue.Scale(targetPOT/tnuePOTsum)
h_elMaxQPhScore_NCnue.Scale(targetPOT/tnuPOTsum)
h_elMaxQPhScore_ext.Scale(targetPOT/textPOTsum)

h_elMaxQPiScore_CCnumu.Scale(targetPOT/tnuPOTsum)
h_elMaxQPiScore_NCnumu.Scale(targetPOT/tnuPOTsum)
h_elMaxQPiScore_CCnue.Scale(targetPOT/tnuePOTsum)
h_elMaxQPiScore_NCnue.Scale(targetPOT/tnuPOTsum)
h_elMaxQPiScore_ext.Scale(targetPOT/textPOTsum)

h_elMaxQConf_CCnumu.Scale(targetPOT/tnuPOTsum)
h_elMaxQConf_NCnumu.Scale(targetPOT/tnuPOTsum)
h_elMaxQConf_CCnue.Scale(targetPOT/tnuePOTsum)
h_elMaxQConf_NCnue.Scale(targetPOT/tnuPOTsum)
h_elMaxQConf_ext.Scale(targetPOT/textPOTsum)

h_elMaxQComp_wConfCut_CCnumu.Scale(targetPOT/tnuPOTsum)
h_elMaxQComp_wConfCut_NCnumu.Scale(targetPOT/tnuPOTsum)
h_elMaxQComp_wConfCut_CCnue.Scale(targetPOT/tnuePOTsum)
h_elMaxQComp_wConfCut_NCnue.Scale(targetPOT/tnuPOTsum)
h_elMaxQComp_wConfCut_ext.Scale(targetPOT/textPOTsum)

h_elMaxQVtxDist_wConfCut_CCnumu.Scale(targetPOT/tnuPOTsum)
h_elMaxQVtxDist_wConfCut_NCnumu.Scale(targetPOT/tnuPOTsum)
h_elMaxQVtxDist_wConfCut_CCnue.Scale(targetPOT/tnuePOTsum)
h_elMaxQVtxDist_wConfCut_NCnue.Scale(targetPOT/tnuPOTsum)
h_elMaxQVtxDist_wConfCut_ext.Scale(targetPOT/textPOTsum)

h_elMaxQCosTheta_wConfCut_CCnumu.Scale(targetPOT/tnuPOTsum)
h_elMaxQCosTheta_wConfCut_NCnumu.Scale(targetPOT/tnuPOTsum)
h_elMaxQCosTheta_wConfCut_CCnue.Scale(targetPOT/tnuePOTsum)
h_elMaxQCosTheta_wConfCut_NCnue.Scale(targetPOT/tnuPOTsum)
h_elMaxQCosTheta_wConfCut_ext.Scale(targetPOT/textPOTsum)

h_elMaxQ_wConfCut_CCnumu.Scale(targetPOT/tnuPOTsum)
h_elMaxQ_wConfCut_NCnumu.Scale(targetPOT/tnuPOTsum)
h_elMaxQ_wConfCut_CCnue.Scale(targetPOT/tnuePOTsum)
h_elMaxQ_wConfCut_NCnue.Scale(targetPOT/tnuPOTsum)
h_elMaxQ_wConfCut_ext.Scale(targetPOT/textPOTsum)

h_elMaxQF_wConfCut_CCnumu.Scale(targetPOT/tnuPOTsum)
h_elMaxQF_wConfCut_NCnumu.Scale(targetPOT/tnuPOTsum)
h_elMaxQF_wConfCut_CCnue.Scale(targetPOT/tnuePOTsum)
h_elMaxQF_wConfCut_NCnue.Scale(targetPOT/tnuPOTsum)
h_elMaxQF_wConfCut_ext.Scale(targetPOT/textPOTsum)

h_elMaxQ_wQFcut_CCnumu.Scale(targetPOT/tnuPOTsum)
h_elMaxQ_wQFcut_NCnumu.Scale(targetPOT/tnuPOTsum)
h_elMaxQ_wQFcut_CCnue.Scale(targetPOT/tnuePOTsum)
h_elMaxQ_wQFcut_NCnue.Scale(targetPOT/tnuPOTsum)
h_elMaxQ_wQFcut_ext.Scale(targetPOT/textPOTsum)

h_elMaxComp_wQFcut_CCnumu.Scale(targetPOT/tnuPOTsum)
h_elMaxComp_wQFcut_NCnumu.Scale(targetPOT/tnuPOTsum)
h_elMaxComp_wQFcut_CCnue.Scale(targetPOT/tnuePOTsum)
h_elMaxComp_wQFcut_NCnue.Scale(targetPOT/tnuPOTsum)
h_elMaxComp_wQFcut_ext.Scale(targetPOT/textPOTsum)

h_elMaxQF_wQcut_CCnumu.Scale(targetPOT/tnuPOTsum)
h_elMaxQF_wQcut_NCnumu.Scale(targetPOT/tnuPOTsum)
h_elMaxQF_wQcut_CCnue.Scale(targetPOT/tnuePOTsum)
h_elMaxQF_wQcut_NCnue.Scale(targetPOT/tnuPOTsum)
h_elMaxQF_wQcut_ext.Scale(targetPOT/textPOTsum)

h_elMaxComp_wQcut_CCnumu.Scale(targetPOT/tnuPOTsum)
h_elMaxComp_wQcut_NCnumu.Scale(targetPOT/tnuPOTsum)
h_elMaxComp_wQcut_CCnue.Scale(targetPOT/tnuePOTsum)
h_elMaxComp_wQcut_NCnue.Scale(targetPOT/tnuPOTsum)
h_elMaxComp_wQcut_ext.Scale(targetPOT/textPOTsum)

h_phScVsPrSum_allMC_bkg.Scale(targetPOT/tnuPOTsum)

h_piScVsPr_CCnumu_bkg.Scale(targetPOT/tnuPOTsum)
h_phScVsPr_CCnumu_bkg.Scale(targetPOT/tnuPOTsum)
h_elScVsPr_CCnumu_bkg.Scale(targetPOT/tnuPOTsum)
h_piScVsPr_NCnumu_bkg.Scale(targetPOT/tnuPOTsum)
h_phScVsPr_NCnumu_bkg.Scale(targetPOT/tnuPOTsum)
h_elScVsPr_NCnumu_bkg.Scale(targetPOT/tnuPOTsum)
h_piScVsPr_NCnue_bkg.Scale(targetPOT/tnuPOTsum)
h_phScVsPr_NCnue_bkg.Scale(targetPOT/tnuPOTsum)
h_elScVsPr_NCnue_bkg.Scale(targetPOT/tnuPOTsum)

h_piScElPrGr0_allMC_bkg.Scale(targetPOT/tnuPOTsum)
h_piScElPrEq0_allMC_bkg.Scale(targetPOT/tnuPOTsum)
h_muScElPrGr0_allMC_bkg.Scale(targetPOT/tnuPOTsum)
h_muScElPrEq0_allMC_bkg.Scale(targetPOT/tnuPOTsum)

h_phScPrSumEq0_allMC_bkg.Scale(targetPOT/tnuPOTsum)
h_phScPrSumGr0_allMC_bkg.Scale(targetPOT/tnuPOTsum)

h_piScPrGr0_allMC_bkg.Scale(targetPOT/tnuPOTsum)
h_piScPrEq0_allMC_bkg.Scale(targetPOT/tnuPOTsum)
h_phScPrGr0_allMC_bkg.Scale(targetPOT/tnuPOTsum)
h_phScPrEq0_allMC_bkg.Scale(targetPOT/tnuPOTsum)
h_elScPrGr0_allMC_bkg.Scale(targetPOT/tnuPOTsum)
h_elScPrEq0_allMC_bkg.Scale(targetPOT/tnuPOTsum)
h_piScPrGr0_CCnumu_bkg.Scale(targetPOT/tnuPOTsum)
h_piScPrEq0_CCnumu_bkg.Scale(targetPOT/tnuPOTsum)
h_phScPrGr0_CCnumu_bkg.Scale(targetPOT/tnuPOTsum)
h_phScPrEq0_CCnumu_bkg.Scale(targetPOT/tnuPOTsum)
h_elScPrGr0_CCnumu_bkg.Scale(targetPOT/tnuPOTsum)
h_elScPrEq0_CCnumu_bkg.Scale(targetPOT/tnuPOTsum)
h_piScPrGr0_NCnumu_bkg.Scale(targetPOT/tnuPOTsum)
h_piScPrEq0_NCnumu_bkg.Scale(targetPOT/tnuPOTsum)
h_phScPrGr0_NCnumu_bkg.Scale(targetPOT/tnuPOTsum)
h_phScPrEq0_NCnumu_bkg.Scale(targetPOT/tnuPOTsum)
h_elScPrGr0_NCnumu_bkg.Scale(targetPOT/tnuPOTsum)
h_elScPrEq0_NCnumu_bkg.Scale(targetPOT/tnuPOTsum)
h_piScPrGr0_NCnue_bkg.Scale(targetPOT/tnuPOTsum)
h_piScPrEq0_NCnue_bkg.Scale(targetPOT/tnuPOTsum)
h_phScPrGr0_NCnue_bkg.Scale(targetPOT/tnuPOTsum)
h_phScPrEq0_NCnue_bkg.Scale(targetPOT/tnuPOTsum)
h_elScPrGr0_NCnue_bkg.Scale(targetPOT/tnuPOTsum)
h_elScPrEq0_NCnue_bkg.Scale(targetPOT/tnuPOTsum)


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

if args.makeKPplots:
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
  
  cnv_tEndKPDist = rt.TCanvas("cnv_tEndKPDist","cnv_tEndKPDist")
  hists_tEndKPDist = sortHists([h_tEndKPDist_CCnumu, h_tEndKPDist_NCnumu, h_tEndKPDist_CCnue, h_tEndKPDist_NCnue, h_tEndKPDist_ext])
  hists_tEndKPDist[0].Draw("EHIST")
  for i in range(1,len(hists_tEndKPDist)):
    hists_tEndKPDist[i].Draw("EHISTSAME")
  leg_tEndKPDist = rt.TLegend(0.7,0.7,0.9,0.9)
  leg_tEndKPDist = configureLegend(leg_tEndKPDist, h_tEndKPDist_CCnumu,
    h_tEndKPDist_NCnumu, h_tEndKPDist_CCnue, h_tEndKPDist_NCnue, h_tEndKPDist_ext)
  leg_tEndKPDist.Draw()
  #cnv_tEndKPDist.SaveAs("tEndKPDist.png")
  cnv_tEndKPDist.Write()

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

cnv_maxElConf = rt.TCanvas("cnv_maxElConf","cnv_maxElConf")
hists_maxElConf = sortHists([h_maxElConf_CCnumu, h_maxElConf_NCnumu, h_maxElConf_CCnue, h_maxElConf_NCnue, h_maxElConf_ext])
hists_maxElConf[0].Draw("EHIST")
for i in range(1,len(hists_maxElConf)):
  hists_maxElConf[i].Draw("EHISTSAME")
leg_maxElConf = rt.TLegend(0.7,0.7,0.9,0.9)
leg_maxElConf = configureLegend(leg_maxElConf, h_maxElConf_CCnumu,
  h_maxElConf_NCnumu, h_maxElConf_CCnue, h_maxElConf_NCnue, h_maxElConf_ext)
leg_maxElConf.Draw()
#cnv_maxElConf.SaveAs("maxElConf.png")
cnv_maxElConf.Write()

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

cnv_elMaxPur = rt.TCanvas("cnv_elMaxPur","cnv_elMaxPur")
hists_elMaxPur = sortHists([h_elMaxPur_CCnumu, h_elMaxPur_NCnumu, h_elMaxPur_CCnue, h_elMaxPur_NCnue, h_elMaxPur_ext])
hists_elMaxPur[0].Draw("EHIST")
for i in range(1,len(hists_elMaxPur)):
  hists_elMaxPur[i].Draw("EHISTSAME")
leg_elMaxPur = rt.TLegend(0.7,0.7,0.9,0.9)
leg_elMaxPur = configureLegend(leg_elMaxPur, h_elMaxPur_CCnumu,
  h_elMaxPur_NCnumu, h_elMaxPur_CCnue, h_elMaxPur_NCnue, h_elMaxPur_ext)
leg_elMaxPur.Draw()
#cnv_elMaxPur.SaveAs("elMaxPur.png")
cnv_elMaxPur.Write()

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

cnv_elMaxQPur = rt.TCanvas("cnv_elMaxQPur","cnv_elMaxQPur")
hists_elMaxQPur = sortHists([h_elMaxQPur_CCnumu, h_elMaxQPur_NCnumu, h_elMaxQPur_CCnue, h_elMaxQPur_NCnue, h_elMaxQPur_ext])
hists_elMaxQPur[0].Draw("EHIST")
for i in range(1,len(hists_elMaxQPur)):
  hists_elMaxQPur[i].Draw("EHISTSAME")
leg_elMaxQPur = rt.TLegend(0.7,0.7,0.9,0.9)
leg_elMaxQPur = configureLegend(leg_elMaxQPur, h_elMaxQPur_CCnumu,
  h_elMaxQPur_NCnumu, h_elMaxQPur_CCnue, h_elMaxQPur_NCnue, h_elMaxQPur_ext)
leg_elMaxQPur.Draw()
#cnv_elMaxQPur.SaveAs("elMaxQPur.png")
cnv_elMaxQPur.Write()

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

cnv_elMaxQVtxDist_wConfCut = rt.TCanvas("cnv_elMaxQVtxDist_wConfCut","cnv_elMaxQVtxDist_wConfCut")
hists_elMaxQVtxDist_wConfCut = sortHists([h_elMaxQVtxDist_wConfCut_CCnumu, h_elMaxQVtxDist_wConfCut_NCnumu, h_elMaxQVtxDist_wConfCut_CCnue, h_elMaxQVtxDist_wConfCut_NCnue, h_elMaxQVtxDist_wConfCut_ext])
hists_elMaxQVtxDist_wConfCut[0].Draw("EHIST")
for i in range(1,len(hists_elMaxQVtxDist_wConfCut)):
  hists_elMaxQVtxDist_wConfCut[i].Draw("EHISTSAME")
leg_elMaxQVtxDist_wConfCut = rt.TLegend(0.7,0.7,0.9,0.9)
leg_elMaxQVtxDist_wConfCut = configureLegend(leg_elMaxQVtxDist_wConfCut, h_elMaxQVtxDist_wConfCut_CCnumu,
  h_elMaxQVtxDist_wConfCut_NCnumu, h_elMaxQVtxDist_wConfCut_CCnue, h_elMaxQVtxDist_wConfCut_NCnue, h_elMaxQVtxDist_wConfCut_ext)
leg_elMaxQVtxDist_wConfCut.Draw()
#cnv_elMaxQVtxDist_wConfCut.SaveAs("elMaxQVtxDist_wConfCut.png")
cnv_elMaxQVtxDist_wConfCut.Write()

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

if args.makeKPplots:
  cnv_cosKPDist_wConfCut = rt.TCanvas("cnv_cosKPDist_wConfCut","cnv_cosKPDist_wConfCut")
  hists_cosKPDist_wConfCut = sortHists([h_cosKPDist_wConfCut_CCnumu, h_cosKPDist_wConfCut_NCnumu, h_cosKPDist_wConfCut_CCnue, h_cosKPDist_wConfCut_NCnue, h_cosKPDist_wConfCut_ext])
  hists_cosKPDist_wConfCut[0].Draw("EHIST")
  for i in range(1,len(hists_cosKPDist_wConfCut)):
    hists_cosKPDist_wConfCut[i].Draw("EHISTSAME")
  leg_cosKPDist_wConfCut = rt.TLegend(0.7,0.7,0.9,0.9)
  leg_cosKPDist_wConfCut = configureLegend(leg_cosKPDist_wConfCut, h_cosKPDist_wConfCut_CCnumu,
    h_cosKPDist_wConfCut_NCnumu, h_cosKPDist_wConfCut_CCnue, h_cosKPDist_wConfCut_NCnue, h_cosKPDist_wConfCut_ext)
  leg_cosKPDist_wConfCut.Draw()
  #cnv_cosKPDist_wConfCut.SaveAs("cosKPDist_wConfCut.png")
  cnv_cosKPDist_wConfCut.Write()
  
  cnv_tEndKPDist_wConfCut = rt.TCanvas("cnv_tEndKPDist_wConfCut","cnv_tEndKPDist_wConfCut")
  hists_tEndKPDist_wConfCut = sortHists([h_tEndKPDist_wConfCut_CCnumu, h_tEndKPDist_wConfCut_NCnumu, h_tEndKPDist_wConfCut_CCnue, h_tEndKPDist_wConfCut_NCnue, h_tEndKPDist_wConfCut_ext])
  hists_tEndKPDist_wConfCut[0].Draw("EHIST")
  for i in range(1,len(hists_tEndKPDist_wConfCut)):
    hists_tEndKPDist_wConfCut[i].Draw("EHISTSAME")
  leg_tEndKPDist_wConfCut = rt.TLegend(0.7,0.7,0.9,0.9)
  leg_tEndKPDist_wConfCut = configureLegend(leg_tEndKPDist_wConfCut, h_tEndKPDist_wConfCut_CCnumu,
    h_tEndKPDist_wConfCut_NCnumu, h_tEndKPDist_wConfCut_CCnue, h_tEndKPDist_wConfCut_NCnue, h_tEndKPDist_wConfCut_ext)
  leg_tEndKPDist_wConfCut.Draw()
  #cnv_tEndKPDist_wConfCut.SaveAs("tEndKPDist_wConfCut.png")
  cnv_tEndKPDist_wConfCut.Write()

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

cnv_piScVsPr_CCnumu_bkg = rt.TCanvas("cnv_piScVsPr_CCnumu_bkg","cnv_piScVsPr_CCnumu_bkg")
h_piScVsPr_CCnumu_bkg.Draw("COLZ")
cnv_piScVsPr_CCnumu_bkg.Write()
cnv_phScVsPr_CCnumu_bkg = rt.TCanvas("cnv_phScVsPr_CCnumu_bkg","cnv_phScVsPr_CCnumu_bkg")
h_phScVsPr_CCnumu_bkg.Draw("COLZ")
cnv_phScVsPr_CCnumu_bkg.Write()
cnv_elScVsPr_CCnumu_bkg = rt.TCanvas("cnv_elScVsPr_CCnumu_bkg","cnv_elScVsPr_CCnumu_bkg")
h_elScVsPr_CCnumu_bkg.Draw("COLZ")
cnv_elScVsPr_CCnumu_bkg.Write()

cnv_piScVsPr_NCnumu_bkg = rt.TCanvas("cnv_piScVsPr_NCnumu_bkg","cnv_piScVsPr_NCnumu_bkg")
h_piScVsPr_NCnumu_bkg.Draw("COLZ")
cnv_piScVsPr_NCnumu_bkg.Write()
cnv_phScVsPr_NCnumu_bkg = rt.TCanvas("cnv_phScVsPr_NCnumu_bkg","cnv_phScVsPr_NCnumu_bkg")
h_phScVsPr_NCnumu_bkg.Draw("COLZ")
cnv_phScVsPr_NCnumu_bkg.Write()
cnv_elScVsPr_NCnumu_bkg = rt.TCanvas("cnv_elScVsPr_NCnumu_bkg","cnv_elScVsPr_NCnumu_bkg")
h_elScVsPr_NCnumu_bkg.Draw("COLZ")
cnv_elScVsPr_NCnumu_bkg.Write()

cnv_piScVsPr_NCnue_bkg = rt.TCanvas("cnv_piScVsPr_NCnue_bkg","cnv_piScVsPr_NCnue_bkg")
h_piScVsPr_NCnue_bkg.Draw("COLZ")
cnv_piScVsPr_NCnue_bkg.Write()
cnv_phScVsPr_NCnue_bkg = rt.TCanvas("cnv_phScVsPr_NCnue_bkg","cnv_phScVsPr_NCnue_bkg")
h_phScVsPr_NCnue_bkg.Draw("COLZ")
cnv_phScVsPr_NCnue_bkg.Write()
cnv_elScVsPr_NCnue_bkg = rt.TCanvas("cnv_elScVsPr_NCnue_bkg","cnv_elScVsPr_NCnue_bkg")
h_elScVsPr_NCnue_bkg.Draw("COLZ")
cnv_elScVsPr_NCnue_bkg.Write()

cnv_phScVsPrSum_allMC_bkg = rt.TCanvas("cnv_phScVsPrSum_allMC_bkg","cnv_phScVsPrSum_allMC_bkg")
h_phScVsPrSum_allMC_bkg.Draw("COLZ")
cnv_phScVsPrSum_allMC_bkg.Write()

cnv_scores_allMC_bkg = rt.TCanvas("cnv_scores_allMC_bkg","cnv_scores_allMC_bkg")
hists_scores_allMC_bkg = sortHists([h_piScPrGr0_allMC_bkg, h_piScPrEq0_allMC_bkg, h_phScPrGr0_allMC_bkg, h_phScPrEq0_allMC_bkg])#, h_elScPrGr0_allMC_bkg, h_elScPrEq0_allMC_bkg])
hists_scores_allMC_bkg[0].Draw("EHIST")
for i in range(1,len(hists_scores_allMC_bkg)):
  hists_scores_allMC_bkg[i].Draw("EHISTSAME")
leg_scores_allMC_bkg = rt.TLegend(0.15,0.7,0.4,0.9)
leg_scores_allMC_bkg.AddEntry(h_piScPrEq0_allMC_bkg, "pion score, pion purity = 0", "l")
leg_scores_allMC_bkg.AddEntry(h_piScPrGr0_allMC_bkg, "pion score, pion purity > 0", "l")
leg_scores_allMC_bkg.AddEntry(h_phScPrEq0_allMC_bkg, "photon score, photon purity = 0", "l")
leg_scores_allMC_bkg.AddEntry(h_phScPrGr0_allMC_bkg, "photon score, photon purity > 0", "l")
#leg_scores_allMC_bkg.AddEntry(h_elScPrEq0_allMC_bkg, "electron score, electron purity = 0", "l")
#leg_scores_allMC_bkg.AddEntry(h_elScPrGr0_allMC_bkg, "electron score, electron purity > 0", "l")
leg_scores_allMC_bkg.Draw()
cnv_scores_allMC_bkg.Write()

cnv_crossScores_allMC_bkg = rt.TCanvas("cnv_crossScores_allMC_bkg","cnv_crossScores_allMC_bkg")
hists_crossScores_allMC_bkg = sortHists([h_piScElPrGr0_allMC_bkg, h_piScElPrEq0_allMC_bkg, h_muScElPrGr0_allMC_bkg, h_muScElPrEq0_allMC_bkg])
hists_crossScores_allMC_bkg[0].Draw("EHIST")
for i in range(1,len(hists_crossScores_allMC_bkg)):
  hists_crossScores_allMC_bkg[i].Draw("EHISTSAME")
leg_crossScores_allMC_bkg = rt.TLegend(0.15,0.7,0.4,0.9)
leg_crossScores_allMC_bkg.AddEntry(h_piScElPrEq0_allMC_bkg, "pion score, electron purity = 0", "l")
leg_crossScores_allMC_bkg.AddEntry(h_piScElPrGr0_allMC_bkg, "pion score, electron purity > 0", "l")
leg_crossScores_allMC_bkg.AddEntry(h_muScElPrEq0_allMC_bkg, "muon score, electron purity = 0", "l")
leg_crossScores_allMC_bkg.AddEntry(h_muScElPrGr0_allMC_bkg, "muon score, electron purity > 0", "l")
leg_crossScores_allMC_bkg.Draw()
cnv_crossScores_allMC_bkg.Write()

cnv_cosmicScores_allMC_bkg = rt.TCanvas("cnv_cosmicScores_allMC_bkg","cnv_cosmicScores_allMC_bkg")
hists_cosmicScores_allMC_bkg = sortHists([h_phScPrSumEq0_allMC_bkg,h_phScPrSumGr0_allMC_bkg])
hists_cosmicScores_allMC_bkg[0].Draw("EHIST")
hists_cosmicScores_allMC_bkg[1].Draw("EHISTSAME")
leg_cosmicScores_allMC_bkg = rt.TLegend(0.15,0.7,0.4,0.9)
leg_cosmicScores_allMC_bkg.AddEntry(h_phScPrSumEq0_allMC_bkg, "photon score, sim particle purity = 0","l")
leg_cosmicScores_allMC_bkg.AddEntry(h_phScPrSumGr0_allMC_bkg, "photon score, sim particle purity > 0 and photon purity = 0","l")
leg_cosmicScores_allMC_bkg.Draw()
cnv_cosmicScores_allMC_bkg.Write()

cnv_piSc_CCnumu_bkg = rt.TCanvas("cnv_piSc_CCnumu_bkg","cnv_piSc_CCnumu_bkg")
hists_piScPr_CCnumu_bkg = sortHists([h_piScPrEq0_CCnumu_bkg, h_piScPrGr0_CCnumu_bkg])
hists_piScPr_CCnumu_bkg[0].Draw("EHIST")
hists_piScPr_CCnumu_bkg[1].Draw("EHISTSAME")
leg_piSc_CCnumu_bkg = rt.TLegend(0.15,0.7,0.4,0.9)
leg_piSc_CCnumu_bkg.AddEntry(h_piScPrEq0_CCnumu_bkg, "pion purity = 0", "l")
leg_piSc_CCnumu_bkg.AddEntry(h_piScPrGr0_CCnumu_bkg, "pion purity > 0", "l")
leg_piSc_CCnumu_bkg.Draw()
cnv_piSc_CCnumu_bkg.Write()

cnv_phSc_CCnumu_bkg = rt.TCanvas("cnv_phSc_CCnumu_bkg","cnv_phSc_CCnumu_bkg")
hists_phScPr_CCnumu_bkg = sortHists([h_phScPrEq0_CCnumu_bkg, h_phScPrGr0_CCnumu_bkg])
hists_phScPr_CCnumu_bkg[0].Draw("EHIST")
hists_phScPr_CCnumu_bkg[1].Draw("EHISTSAME")
leg_phSc_CCnumu_bkg = rt.TLegend(0.15,0.7,0.4,0.9)
leg_phSc_CCnumu_bkg.AddEntry(h_phScPrEq0_CCnumu_bkg, "photon purity = 0", "l")
leg_phSc_CCnumu_bkg.AddEntry(h_phScPrGr0_CCnumu_bkg, "photon purity > 0", "l")
leg_phSc_CCnumu_bkg.Draw()
cnv_phSc_CCnumu_bkg.Write()

cnv_elSc_CCnumu_bkg = rt.TCanvas("cnv_elSc_CCnumu_bkg","cnv_elSc_CCnumu_bkg")
hists_elScPr_CCnumu_bkg = sortHists([h_elScPrEq0_CCnumu_bkg, h_elScPrGr0_CCnumu_bkg])
hists_elScPr_CCnumu_bkg[0].Draw("EHIST")
hists_elScPr_CCnumu_bkg[1].Draw("EHISTSAME")
leg_elSc_CCnumu_bkg = rt.TLegend(0.15,0.7,0.4,0.9)
leg_elSc_CCnumu_bkg.AddEntry(h_elScPrEq0_CCnumu_bkg, "electron purity = 0", "l")
leg_elSc_CCnumu_bkg.AddEntry(h_elScPrGr0_CCnumu_bkg, "electron purity > 0", "l")
leg_elSc_CCnumu_bkg.Draw()
cnv_elSc_CCnumu_bkg.Write()

cnv_piSc_NCnumu_bkg = rt.TCanvas("cnv_piSc_NCnumu_bkg","cnv_piSc_NCnumu_bkg")
hists_piScPr_NCnumu_bkg = sortHists([h_piScPrEq0_NCnumu_bkg, h_piScPrGr0_NCnumu_bkg])
hists_piScPr_NCnumu_bkg[0].Draw("EHIST")
hists_piScPr_NCnumu_bkg[1].Draw("EHISTSAME")
leg_piSc_NCnumu_bkg = rt.TLegend(0.15,0.7,0.4,0.9)
leg_piSc_NCnumu_bkg.AddEntry(h_piScPrEq0_NCnumu_bkg, "pion purity = 0", "l")
leg_piSc_NCnumu_bkg.AddEntry(h_piScPrGr0_NCnumu_bkg, "pion purity > 0", "l")
leg_piSc_NCnumu_bkg.Draw()
cnv_piSc_NCnumu_bkg.Write()

cnv_phSc_NCnumu_bkg = rt.TCanvas("cnv_phSc_NCnumu_bkg","cnv_phSc_NCnumu_bkg")
hists_phScPr_NCnumu_bkg = sortHists([h_phScPrEq0_NCnumu_bkg, h_phScPrGr0_NCnumu_bkg])
hists_phScPr_NCnumu_bkg[0].Draw("EHIST")
hists_phScPr_NCnumu_bkg[1].Draw("EHISTSAME")
leg_phSc_NCnumu_bkg = rt.TLegend(0.15,0.7,0.4,0.9)
leg_phSc_NCnumu_bkg.AddEntry(h_phScPrEq0_NCnumu_bkg, "photon purity = 0", "l")
leg_phSc_NCnumu_bkg.AddEntry(h_phScPrGr0_NCnumu_bkg, "photon purity > 0", "l")
leg_phSc_NCnumu_bkg.Draw()
cnv_phSc_NCnumu_bkg.Write()

cnv_elSc_NCnumu_bkg = rt.TCanvas("cnv_elSc_NCnumu_bkg","cnv_elSc_NCnumu_bkg")
hists_elScPr_NCnumu_bkg = sortHists([h_elScPrEq0_NCnumu_bkg, h_elScPrGr0_NCnumu_bkg])
hists_elScPr_NCnumu_bkg[0].Draw("EHIST")
hists_elScPr_NCnumu_bkg[1].Draw("EHISTSAME")
leg_elSc_NCnumu_bkg = rt.TLegend(0.15,0.7,0.4,0.9)
leg_elSc_NCnumu_bkg.AddEntry(h_elScPrEq0_NCnumu_bkg, "electron purity = 0", "l")
leg_elSc_NCnumu_bkg.AddEntry(h_elScPrGr0_NCnumu_bkg, "electron purity > 0", "l")
leg_elSc_NCnumu_bkg.Draw()
cnv_elSc_NCnumu_bkg.Write()

cnv_piSc_NCnue_bkg = rt.TCanvas("cnv_piSc_NCnue_bkg","cnv_piSc_NCnue_bkg")
hists_piScPr_NCnue_bkg = sortHists([h_piScPrEq0_NCnue_bkg, h_piScPrGr0_NCnue_bkg])
hists_piScPr_NCnue_bkg[0].Draw("EHIST")
hists_piScPr_NCnue_bkg[1].Draw("EHISTSAME")
leg_piSc_NCnue_bkg = rt.TLegend(0.15,0.7,0.4,0.9)
leg_piSc_NCnue_bkg.AddEntry(h_piScPrEq0_NCnue_bkg, "pion purity = 0", "l")
leg_piSc_NCnue_bkg.AddEntry(h_piScPrGr0_NCnue_bkg, "pion purity > 0", "l")
leg_piSc_NCnue_bkg.Draw()
cnv_piSc_NCnue_bkg.Write()

cnv_phSc_NCnue_bkg = rt.TCanvas("cnv_phSc_NCnue_bkg","cnv_phSc_NCnue_bkg")
hists_phScPr_NCnue_bkg = sortHists([h_phScPrEq0_NCnue_bkg, h_phScPrGr0_NCnue_bkg])
hists_phScPr_NCnue_bkg[0].Draw("EHIST")
hists_phScPr_NCnue_bkg[1].Draw("EHISTSAME")
leg_phSc_NCnue_bkg = rt.TLegend(0.15,0.7,0.4,0.9)
leg_phSc_NCnue_bkg.AddEntry(h_phScPrEq0_NCnue_bkg, "photon purity = 0", "l")
leg_phSc_NCnue_bkg.AddEntry(h_phScPrGr0_NCnue_bkg, "photon purity > 0", "l")
leg_phSc_NCnue_bkg.Draw()
cnv_phSc_NCnue_bkg.Write()

cnv_elSc_NCnue_bkg = rt.TCanvas("cnv_elSc_NCnue_bkg","cnv_elSc_NCnue_bkg")
hists_elScPr_NCnue_bkg = sortHists([h_elScPrEq0_NCnue_bkg, h_elScPrGr0_NCnue_bkg])
hists_elScPr_NCnue_bkg[0].Draw("EHIST")
hists_elScPr_NCnue_bkg[1].Draw("EHISTSAME")
leg_elSc_NCnue_bkg = rt.TLegend(0.15,0.7,0.4,0.9)
leg_elSc_NCnue_bkg.AddEntry(h_elScPrEq0_NCnue_bkg, "electron purity = 0", "l")
leg_elSc_NCnue_bkg.AddEntry(h_elScPrGr0_NCnue_bkg, "electron purity > 0", "l")
leg_elSc_NCnue_bkg.Draw()
cnv_elSc_NCnue_bkg.Write()

#h_nuE_CCnue_eff.GetYaxis().SetRangeUser(0,1.003)
#h_nuE_CCnue_pur.GetYaxis().SetRangeUser(0,1.003)
h_nuE_CCnue_eff.GetYaxis().SetTitle("efficiency")
cnv_CCnue_sel_trueE = rt.TCanvas("cnv_CCnue_sel_trueE","cnv_CCnue_sel_trueE")
cnv_CCnue_sel_trueE.SetGrid()
h_nuE_CCnue_eff.Draw("E")
#h_nuE_CCnue_pur.Draw("ESAME")
#leg_CCnue_sel_trueE = rt.TLegend(0.7,0.7,0.9,0.9)
#leg_CCnue_sel_trueE.AddEntry(h_nuE_CCnue_eff, "efficiency", "l")
#leg_CCnue_sel_trueE.AddEntry(h_nuE_CCnue_pur, "purity", "l")
#leg_CCnue_sel_trueE.Draw()
cnv_CCnue_sel_trueE.Write()

h_nuEr_CCnue_eff.GetYaxis().SetRangeUser(0,1.003)
h_nuEr_CCnue_pur.GetYaxis().SetRangeUser(0,1.003)
cnv_CCnue_sel_recoE = rt.TCanvas("cnv_CCnue_sel_recoE","cnv_CCnue_sel_recoE")
cnv_CCnue_sel_recoE.SetGrid()
h_nuEr_CCnue_eff.Draw("E")
h_nuEr_CCnue_pur.Draw("ESAME")
leg_CCnue_sel_recoE = rt.TLegend(0.7,0.7,0.9,0.9)
leg_CCnue_sel_recoE.AddEntry(h_nuEr_CCnue_eff, "efficiency", "l")
leg_CCnue_sel_recoE.AddEntry(h_nuEr_CCnue_pur, "purity", "l")
leg_CCnue_sel_recoE.Draw()
cnv_CCnue_sel_recoE.Write()

cnv_visE_sel = rt.TCanvas("cnv_visE_sel", "cnv_visE_sel")
h_visE_data_wCuts.Draw("E")
h_visE_all_wCuts.Draw("hist same")
h_visE_data_wCuts.Draw("ESAME")
leg_visE_sel = rt.TLegend(0.7,0.7,0.9,0.9)
leg_visE_sel.AddEntry(h_visE_ext_wCuts, "cosmic background (%.2f)"%h_visE_ext_wCuts.Integral(), "f")
leg_visE_sel.AddEntry(h_visE_NCnue_wCuts, "NC nue (%.2f)"%h_visE_NCnue_wCuts.Integral(), "f")
leg_visE_sel.AddEntry(h_visE_NCnumu_wCuts, "NC numu (%.2f)"%h_visE_NCnumu_wCuts.Integral(), "f")
leg_visE_sel.AddEntry(h_visE_CCnumu_wCuts, "CC numu (%.2f)"%h_visE_CCnumu_wCuts.Integral(), "f")
leg_visE_sel.AddEntry(h_visE_CCnue_wCuts, "CC nue (%.2f)"%h_visE_CCnue_wCuts.Integral(), "f")
leg_visE_sel.AddEntry(h_visE_data_wCuts, "5e19 data (%.2f)"%h_visE_data_wCuts.Integral(), "l")
leg_visE_sel.Draw()
cnv_visE_sel.Write()

h_nuE_CCnue_nCuts.Write()
h_nuE_CCnue_wCuts.Write()
h_nuE_all_wCuts.Write()
h_nuE_CCnue_eff.Write()
h_nuE_CCnue_pur.Write()

h_nuEr_CCnue_nCuts.Write()
h_nuEr_CCnue_wCuts.Write()
h_nuEr_all_wCuts.Write()
h_nuEr_CCnue_eff.Write()
h_nuEr_CCnue_pur.Write()

h_visE_all_wCuts.Write()
h_visE_data_wCuts.Write()
h_visE_CCnue_wCuts.Write()
h_visE_CCnumu_wCuts.Write()
h_visE_NCnumu_wCuts.Write()
h_visE_NCnue_wCuts.Write()
h_visE_ext_wCuts.Write()


#hMain_xTS = 0.0425
#hMain_yTS = 0.05
#hMain_yTO = 0.6
#hRat_xTS = 0.13
#hRat_yTS = 0.14
#hRat_xLS = 0.1
#hRat_yLS = 0.1
#hRat_xTO = 1.
#hRat_yTO = 0.2
#hRat_y1 = -1
#hRat_y2 = 3
#pad_x1 = 0.005
#pad_x2 = 0.995
#uPad_y1 = 0.2525
#uPad_y2 = 0.995
#lPad_y1 = 0.005
#lPad_y2 = 0.2475
#bMargin = 0.3
#h_nuE_CCnue_nCuts.GetXaxis().SetTitleSize(hMain_xTS)
#h_nuE_CCnue_nCuts.GetYaxis().SetTitleSize(hMain_yTS)
#h_nuE_CCnue_nCuts.GetYaxis().SetTitleOffset(hMain_yTO)
#h_nuE_CCnue_wCuts.GetXaxis().SetTitleSize(hMain_xTS)
#h_nuE_CCnue_wCuts.GetYaxis().SetTitleSize(hMain_yTS)
#h_nuE_CCnue_wCuts.GetYaxis().SetTitleOffset(hMain_yTO)
#h_nuE_CCnue_eff.SetTitle("")
#h_nuE_CCnue_eff.GetYaxis().SetTitle("efficiency")
#h_nuE_CCnue_eff.GetXaxis().SetTitleSize(hRat_xTS)
#h_nuE_CCnue_eff.GetYaxis().SetTitleSize(hRat_yTS)
#h_nuE_CCnue_eff.GetXaxis().SetLabelSize(hRat_xLS)
#h_nuE_CCnue_eff.GetYaxis().SetLabelSize(hRat_yLS)
#h_nuE_CCnue_eff.GetXaxis().SetTitleOffset(hRat_xTO)
#h_nuE_CCnue_eff.GetYaxis().SetTitleOffset(hRat_yTO)
#h_nuE_CCnue_eff.SetLineColor(rt.kBlack)
#cnv_nuE_CCnue = rt.TCanvas("cnv_nuE_CCnue","cnv_nuE_CCnue")
#uPad_nuE_CCnue = rt.TPad("uPad_nuE_CCnue","uPad_nuE_CCnue",pad_x1,uPad_y1,pad_x2,uPad_y2)
#lPad_nuE_CCnue = rt.TPad("lPad_nuE_CCnue","lPad_nuE_CCnue",pad_x1,lPad_y1,pad_x2,lPad_y2)
#lPad_nuE_CCnue.SetBottomMargin(bMargin)
#uPad_nuE_CCnue.Draw()
#lPad_nuE_CCnue.Draw()
#uPad_nuE_CCnue.cd()
#h_nuE_CCnue_nCuts.Draw("EHIST")
#h_nuE_CCnue_wCuts.Draw("EHISTSAME")
#leg_nuE_CCnue = rt.TLegend(0.7,0.7,0.9,0.9)
#leg_nuE_CCnue.AddEntry(h_nuE_CCnue_nCuts, "Without Cuts", "l")
#leg_nuE_CCnue.AddEntry(h_nuE_CCnue_wCuts, "With Cuts", "l")
#leg_nuE_CCnue.Draw()
#lPad_nuE_CCnue.cd()
#h_nuE_CCnue_eff.Draw("E")
#cnv_nuE_CCnue.Write()


