
import sys, argparse
import numpy as np
import ROOT as rt

from math import isinf
from helpers.plotting_functions import sortHists


parser = argparse.ArgumentParser("Plot Selection Test Results")
parser.add_argument("-fnu", "--bnbnu_file", type=str, default="prepare_selection_test_output/prepare_selection_test_reco_v2me05_gen2val_v15_nu_file.root", help="bnb nu input file")
parser.add_argument("-fnue", "--bnbnue_file", type=str, default="prepare_selection_test_output/prepare_selection_test_reco_v2me05_gen2val_v15_nue_file.root", help="bnb nu input file")
parser.add_argument("-fext", "--extbnb_file", type=str, default="prepare_selection_test_output/prepare_selection_test_reco_v2me05_gen2val_v15_extbnb_file.root", help="bnb nu input file")
parser.add_argument("-c", "--compCut", type=float, default=0.6, help="completeness cut value")
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

#970: number of run3 extbnb merged_dlana files in prepare_selection_test_reco_v2me05_gen2val_v15_extbnb_file.root
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

n_runs1to3_CCnumu = 0
n_runs1to3_CCnue = 0


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
  elMaxQFrac = -1.
  elMaxQ = -1.
  muMaxQ = -1.
  muMaxComp = -1.
  muMaxQFrac = -1.

  for iT in range(tnu.nTracks):
    if tnu.trackClassified[iT] == 1 and tnu.trackRecoPID[iT] == 13:
      nMuons += 1
      if tnu.trackRecoComp[iT] > args.compCut:
        nCompMuons += 1
      if tnu.trackCharge[iT] > muMaxQ:
        muMaxQ = tnu.trackCharge[iT]
      if tnu.trackChargeFrac[iT] > muMaxQFrac:
        muMaxQFrac = tnu.trackChargeFrac[iT]
      if tnu.trackRecoComp[iT] > muMaxComp:
        muMaxComp = tnu.trackRecoComp[iT]
  for iS in range(tnu.nShowers):
    if tnu.showerClassified[iS] == 1 and tnu.showerRecoPID[iS] == 11:
      nElectrons += 1
      if tnu.showerRecoComp[iS] > args.compCut:
        nCompElectrons += 1
      if tnu.showerCharge[iS] > elMaxQ:
        elMaxQ = tnu.showerCharge[iS]
      if tnu.showerChargeFrac[iS] > elMaxQFrac:
        elMaxQFrac = tnu.showerChargeFrac[iS]
      if tnu.showerRecoComp[iS] > elMaxComp:
        elMaxComp = tnu.showerRecoComp[iS]

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

  if nMuons >= 1:
    h_muMaxComp_CCnumu, h_muMaxComp_NCnumu, h_muMaxComp_NCnue = FillNuHistos(h_muMaxComp_CCnumu,
      h_muMaxComp_NCnumu, h_muMaxComp_NCnue, muMaxComp, tnu.xsecWeight, eventType)
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
  elMaxQFrac = -1.
  elMaxQ = -1.
  muMaxQ = -1.
  muMaxComp = -1.
  muMaxQFrac = -1.

  for iT in range(tnue.nTracks):
    if tnue.trackClassified[iT] == 1 and tnue.trackRecoPID[iT] == 13:
      nMuons += 1
      if tnue.trackRecoComp[iT] > args.compCut:
        nCompMuons += 1
      if tnue.trackCharge[iT] > muMaxQ:
        muMaxQ = tnue.trackCharge[iT]
      if tnue.trackChargeFrac[iT] > muMaxQFrac:
        muMaxQFrac = tnue.trackChargeFrac[iT]
      if tnue.trackRecoComp[iT] > muMaxComp:
        muMaxComp = tnue.trackRecoComp[iT]
  for iS in range(tnue.nShowers):
    if tnue.showerClassified[iS] == 1 and tnue.showerRecoPID[iS] == 11:
      nElectrons += 1
      if tnue.showerRecoComp[iS] > args.compCut:
        nCompElectrons += 1
      if tnue.showerCharge[iS] > elMaxQ:
        elMaxQ = tnue.showerCharge[iS]
      if tnue.showerChargeFrac[iS] > elMaxQFrac:
        elMaxQFrac = tnue.showerChargeFrac[iS]
      if tnue.showerRecoComp[iS] > elMaxComp:
        elMaxComp = tnue.showerRecoComp[iS]

  h_nEl_CCnue.Fill(nElectrons, tnue.xsecWeight)
  h_nMu_CCnue.Fill(nMuons, tnue.xsecWeight)
  h_nCompEl_CCnue.Fill(nCompElectrons, tnue.xsecWeight)
  h_nCompMu_CCnue.Fill(nCompMuons, tnue.xsecWeight)

  if nElectrons >= 1:
    h_elMaxComp_CCnue.Fill(elMaxComp, tnue.xsecWeight)
    h_elMaxQF_CCnue.Fill(elMaxQFrac, tnue.xsecWeight)
    h_elMaxQ_CCnue.Fill(elMaxQ, tnue.xsecWeight)
    if elMaxQFrac > args.chargeFracCut:
      h_elMaxComp_wQFcut_CCnue.Fill(elMaxComp, tnue.xsecWeight)
      h_elMaxQ_wQFcut_CCnue.Fill(elMaxQ, tnue.xsecWeight)
    if elMaxQ > args.chargeCut:
      h_elMaxComp_wQcut_CCnue.Fill(elMaxComp, tnue.xsecWeight)
      h_elMaxQF_wQcut_CCnue.Fill(elMaxQFrac, tnue.xsecWeight)

  if nMuons >= 1:
    h_muMaxComp_CCnue.Fill(muMaxComp, tnue.xsecWeight)
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
  elMaxQFrac = -1.
  elMaxQ = -1.
  muMaxQ = -1.
  muMaxComp = -1.
  muMaxQFrac = -1.
  for iT in range(text.nTracks):
    if text.trackClassified[iT] == 1 and text.trackRecoPID[iT] == 13:
      nMuons += 1
      if text.trackRecoComp[iT] > args.compCut:
        nCompMuons += 1
      if text.trackCharge[iT] > muMaxQ:
        muMaxQ = text.trackCharge[iT]
      if text.trackChargeFrac[iT] > muMaxQFrac:
        muMaxQFrac = text.trackChargeFrac[iT]
      if text.trackRecoComp[iT] > muMaxComp:
        muMaxComp = text.trackRecoComp[iT]
  for iS in range(text.nShowers):
    if text.showerClassified[iS] == 1 and text.showerRecoPID[iS] == 11:
      nElectrons += 1
      if text.showerRecoComp[iS] > args.compCut:
        nCompElectrons += 1
      if text.showerCharge[iS] > elMaxQ:
        elMaxQ = text.showerCharge[iS]
      if text.showerChargeFrac[iS] > elMaxQFrac:
        elMaxQFrac = text.showerChargeFrac[iS]
      if text.showerRecoComp[iS] > elMaxComp:
        elMaxComp = text.showerRecoComp[iS]

  h_nEl_ext.Fill(nElectrons)
  h_nMu_ext.Fill(nMuons)
  h_nCompEl_ext.Fill(nCompElectrons)
  h_nCompMu_ext.Fill(nCompMuons)

  if nElectrons >= 1:
    h_elMaxComp_ext.Fill(elMaxComp)
    h_elMaxQF_ext.Fill(elMaxQFrac)
    h_elMaxQ_ext.Fill(elMaxQ)
    if elMaxQFrac > args.chargeFracCut:
      h_elMaxComp_wQFcut_ext.Fill(elMaxComp)
      h_elMaxQ_wQFcut_ext.Fill(elMaxQ)
    if elMaxQ > args.chargeCut:
      h_elMaxComp_wQcut_ext.Fill(elMaxComp)
      h_elMaxQF_wQcut_ext.Fill(elMaxQFrac)

  if nMuons >= 1:
    h_muMaxComp_ext.Fill(muMaxComp)
    h_muMaxQF_ext.Fill(muMaxQFrac)
    h_muMaxQ_ext.Fill(muMaxQ)




n_runs1to3_CCnumu *= (runs1to3POT/tnuPOTsum)
n_runs1to3_CCnue *= (runs1to3POT/tnuePOTsum)

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


print("h_cosFrac_CCnue unscaled integral: ", h_cosFrac_CCnue.Integral())
print("h_nEl_CCnue unscaled integral: ", h_nEl_CCnue.Integral())

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

print("h_cosFrac_CCnue scaled integral: ", h_cosFrac_CCnue.Integral())
print("h_nEl_CCnue scaled integral: ", h_nEl_CCnue.Integral())

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


