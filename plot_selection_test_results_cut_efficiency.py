
import sys, argparse
import numpy as np
import ROOT as rt

from math import sqrt, isinf
from helpers.plotting_functions import sortHists


parser = argparse.ArgumentParser("Plot Selection Test Results")
parser.add_argument("-fnue", "--bnbnue_file", type=str, default="flat_ntuples/dlgen2_reco_v2me06_ntuple_v5_mcc9_v28_wctagger_nueintrinsics.root", help="bnb nue input file")
parser.add_argument("-o", "--outfile", type=str, default="selection_output/plot_selection_test_results_output/plot_selection_test_results_cut_efficiency_output.root", help="output root file name")
args = parser.parse_args()

rt.TH1.SetDefaultSumw2(rt.kTRUE)
rt.gStyle.SetOptStat(0)

fnue = rt.TFile(args.bnbnue_file)
tnue = fnue.Get("EventTree")
tnuePOT = fnue.Get("potTree")

#targetPOT = 6.67e+20
targetPOT = 4.4e+19
targetPOTstr = "4.4e+19"

tnuePOTsum = 0.
for i in range(tnuePOT.GetEntries()):
  tnuePOT.GetEntry(i)
  tnuePOTsum = tnuePOTsum + tnuePOT.totGoodPOT


def configureHist(hist, color, particle, ylabel):
  hist.GetXaxis().SetTitle("true "+particle+" energy (GeV)")
  hist.GetYaxis().SetTitle(ylabel)
  hist.SetLineColor(color)
  hist.SetLineWidth(2)
  return hist

h_el_nCuts = rt.TH1F("h_el_nCuts","Inclusive CCnue Electron Energy",40,0,4)
h_el_wCuts1 = rt.TH1F("h_el_wCuts1","Inclusive CCnue Electron Energy",40,0,4)
h_el_wCuts2 = rt.TH1F("h_el_wCuts2","Inclusive CCnue Electron Energy",40,0,4)
h_el_wCuts3 = rt.TH1F("h_el_wCuts3","Inclusive CCnue Electron Energy",40,0,4)
h_el_wCuts4 = rt.TH1F("h_el_wCuts4","Inclusive CCnue Electron Energy",40,0,4)
h_el_wCuts5 = rt.TH1F("h_el_wCuts5","Inclusive CCnue Electron Energy",40,0,4)
h_el_wCuts6 = rt.TH1F("h_el_wCuts6","Inclusive CCnue Electron Energy",40,0,4)
h_el_wCuts7 = rt.TH1F("h_el_wCuts7","Inclusive CCnue Electron Energy",40,0,4)
#h_el_wCuts8 = rt.TH1F("h_el_wCuts8","Inclusive CCnue Electron Energy",40,0,4)
#h_el_wCuts9 = rt.TH1F("h_el_wCuts9","Inclusive CCnue Electron Energy",40,0,4)
#h_el_wCuts10 = rt.TH1F("h_el_wCuts10","Inclusive CCnue Electron Energy",40,0,4)
h_el_eff1 = rt.TH1F("h_el_eff1","Inclusive CCnue Efficiency",40,0,4)
h_el_eff2 = rt.TH1F("h_el_eff2","Inclusive CCnue Efficiency",40,0,4)
h_el_eff3 = rt.TH1F("h_el_eff3","Inclusive CCnue Efficiency",40,0,4)
h_el_eff4 = rt.TH1F("h_el_eff4","Inclusive CCnue Efficiency",40,0,4)
h_el_eff5 = rt.TH1F("h_el_eff5","Inclusive CCnue Efficiency",40,0,4)
h_el_eff6 = rt.TH1F("h_el_eff6","Inclusive CCnue Efficiency",40,0,4)
h_el_eff7 = rt.TH1F("h_el_eff7","Inclusive CCnue Efficiency",40,0,4)
#h_el_eff8 = rt.TH1F("h_el_eff8","Inclusive CCnue Efficiency",40,0,4)
#h_el_eff9 = rt.TH1F("h_el_eff9","Inclusive CCnue Efficiency",40,0,4)
#h_el_eff10 = rt.TH1F("h_el_eff10","Inclusive CCnue Efficiency",40,0,4)
h_el_nCuts = configureHist(h_el_nCuts, 1, "electron", "events per "+targetPOTstr+" POT")
h_el_wCuts1 = configureHist(h_el_wCuts1, 4, "electron", "events per "+targetPOTstr+" POT")
h_el_wCuts2 = configureHist(h_el_wCuts2, 3, "electron", "events per "+targetPOTstr+" POT")
h_el_wCuts3 = configureHist(h_el_wCuts3, 39, "electron", "events per "+targetPOTstr+" POT")
h_el_wCuts4 = configureHist(h_el_wCuts4, 8, "electron", "events per "+targetPOTstr+" POT")
h_el_wCuts5 = configureHist(h_el_wCuts5, 6, "electron", "events per "+targetPOTstr+" POT")
h_el_wCuts6 = configureHist(h_el_wCuts6, 40, "electron", "events per "+targetPOTstr+" POT")
h_el_wCuts7 = configureHist(h_el_wCuts7, 46, "electron", "events per "+targetPOTstr+" POT")
#h_el_wCuts8 = configureHist(h_el_wCuts8, 9, "electron", "events per "+targetPOTstr+" POT")
#h_el_wCuts9 = configureHist(h_el_wCuts9, 41, "electron", "events per "+targetPOTstr+" POT")
#h_el_wCuts10 = configureHist(h_el_wCuts10, 2, "electron", "events per "+targetPOTstr+" POT")
h_el_eff1 = configureHist(h_el_eff1, 4, "electron", "efficiency")
h_el_eff2 = configureHist(h_el_eff2, 3, "electron", "efficiency")
h_el_eff3 = configureHist(h_el_eff3, 39, "electron", "efficiency")
h_el_eff4 = configureHist(h_el_eff4, 8, "electron", "efficiency")
h_el_eff5 = configureHist(h_el_eff5, 6, "electron", "efficiency")
h_el_eff6 = configureHist(h_el_eff6, 40, "electron", "efficiency")
h_el_eff7 = configureHist(h_el_eff7, 46, "electron", "efficiency")
#h_el_eff8 = configureHist(h_el_eff8, 9, "electron", "efficiency")
#h_el_eff9 = configureHist(h_el_eff9, 41, "electron", "efficiency")
#h_el_eff10 = configureHist(h_el_eff10, 2, "electron", "efficiency")

h_nu_nCuts = rt.TH1F("h_nu_nCuts","Inclusive CCnue Neutrino Energy",40,0,4)
h_nu_wCuts1 = rt.TH1F("h_nu_wCuts1","Inclusive CCnue Neutrino Energy",40,0,4)
h_nu_wCuts2 = rt.TH1F("h_nu_wCuts2","Inclusive CCnue Neutrino Energy",40,0,4)
h_nu_wCuts3 = rt.TH1F("h_nu_wCuts3","Inclusive CCnue Neutrino Energy",40,0,4)
h_nu_wCuts4 = rt.TH1F("h_nu_wCuts4","Inclusive CCnue Neutrino Energy",40,0,4)
h_nu_wCuts5 = rt.TH1F("h_nu_wCuts5","Inclusive CCnue Neutrino Energy",40,0,4)
h_nu_wCuts6 = rt.TH1F("h_nu_wCuts6","Inclusive CCnue Neutrino Energy",40,0,4)
h_nu_wCuts7 = rt.TH1F("h_nu_wCuts7","Inclusive CCnue Neutrino Energy",40,0,4)
#h_nu_wCuts8 = rt.TH1F("h_nu_wCuts8","Inclusive CCnue Neutrino Energy",40,0,4)
#h_nu_wCuts9 = rt.TH1F("h_nu_wCuts9","Inclusive CCnue Neutrino Energy",40,0,4)
#h_nu_wCuts10 = rt.TH1F("h_nu_wCuts10","Inclusive CCnue Neutrino Energy",40,0,4)
h_nu_eff1 = rt.TH1F("h_nu_eff1","Inclusive CCnue Efficiency",40,0,4)
h_nu_eff2 = rt.TH1F("h_nu_eff2","Inclusive CCnue Efficiency",40,0,4)
h_nu_eff3 = rt.TH1F("h_nu_eff3","Inclusive CCnue Efficiency",40,0,4)
h_nu_eff4 = rt.TH1F("h_nu_eff4","Inclusive CCnue Efficiency",40,0,4)
h_nu_eff5 = rt.TH1F("h_nu_eff5","Inclusive CCnue Efficiency",40,0,4)
h_nu_eff6 = rt.TH1F("h_nu_eff6","Inclusive CCnue Efficiency",40,0,4)
h_nu_eff7 = rt.TH1F("h_nu_eff7","Inclusive CCnue Efficiency",40,0,4)
#h_nu_eff8 = rt.TH1F("h_nu_eff8","Inclusive CCnue Efficiency",40,0,4)
#h_nu_eff9 = rt.TH1F("h_nu_eff9","Inclusive CCnue Efficiency",40,0,4)
#h_nu_eff10 = rt.TH1F("h_nu_eff10","Inclusive CCnue Efficiency",40,0,4)
h_nu_nCuts = configureHist(h_nu_nCuts, 1, "neutrino", "events per "+targetPOTstr+" POT")
h_nu_wCuts1 = configureHist(h_nu_wCuts1, 4, "neutrino", "events per "+targetPOTstr+" POT")
h_nu_wCuts2 = configureHist(h_nu_wCuts2, 3, "neutrino", "events per "+targetPOTstr+" POT")
h_nu_wCuts3 = configureHist(h_nu_wCuts3, 39, "neutrino", "events per "+targetPOTstr+" POT")
h_nu_wCuts4 = configureHist(h_nu_wCuts4, 8, "neutrino", "events per "+targetPOTstr+" POT")
h_nu_wCuts5 = configureHist(h_nu_wCuts5, 6, "neutrino", "events per "+targetPOTstr+" POT")
h_nu_wCuts6 = configureHist(h_nu_wCuts6, 40, "neutrino", "events per "+targetPOTstr+" POT")
h_nu_wCuts7 = configureHist(h_nu_wCuts7, 46, "neutrino", "events per "+targetPOTstr+" POT")
#h_nu_wCuts8 = configureHist(h_nu_wCuts8, 9, "neutrino", "events per "+targetPOTstr+" POT")
#h_nu_wCuts9 = configureHist(h_nu_wCuts9, 41, "neutrino", "events per "+targetPOTstr+" POT")
#h_nu_wCuts10 = configureHist(h_nu_wCuts10, 2, "neutrino", "events per "+targetPOTstr+" POT")
h_nu_eff1 = configureHist(h_nu_eff1, 4, "neutrino", "efficiency")
h_nu_eff2 = configureHist(h_nu_eff2, 3, "neutrino", "efficiency")
h_nu_eff3 = configureHist(h_nu_eff3, 39, "neutrino", "efficiency")
h_nu_eff4 = configureHist(h_nu_eff4, 8, "neutrino", "efficiency")
h_nu_eff5 = configureHist(h_nu_eff5, 6, "neutrino", "efficiency")
h_nu_eff6 = configureHist(h_nu_eff6, 40, "neutrino", "efficiency")
h_nu_eff7 = configureHist(h_nu_eff7, 46, "neutrino", "efficiency")
#h_nu_eff8 = configureHist(h_nu_eff8, 9, "neutrino", "efficiency")
#h_nu_eff9 = configureHist(h_nu_eff9, 41, "neutrino", "efficiency")
#h_nu_eff10 = configureHist(h_nu_eff10, 2, "neutrino", "efficiency")


n_raw_all = 0
n_scaled_all = 0.
n_scaled_passCuts1 = 0.
n_scaled_passCuts2 = 0.
n_scaled_passCuts3 = 0.
n_scaled_passCuts4 = 0.
n_scaled_passCuts5 = 0.
n_scaled_passCuts6 = 0.
n_scaled_passCuts7 = 0.
#n_scaled_passCuts8 = 0.
#n_scaled_passCuts9 = 0.
#n_scaled_passCuts10 = 0.


for i in range(tnue.GetEntries()):

  tnue.GetEntry(i)

  if abs(tnue.trueNuPDG) != 12 or abs(tnue.trueLepPDG) != 11 or tnue.trueNuCCNC != 0 or isinf(tnue.xsecWeight):
    continue

  n_raw_all += 1
  n_scaled_all += tnue.xsecWeight
  h_el_nCuts.Fill(tnue.trueLepE, tnue.xsecWeight)
  h_nu_nCuts.Fill(tnue.trueNuE, tnue.xsecWeight)

  if tnue.foundVertex != 1:
    continue

  #h_el_wCuts1.Fill(tnue.trueLepE, tnue.xsecWeight)
  #h_nu_wCuts1.Fill(tnue.trueNuE, tnue.xsecWeight)
  #n_scaled_passCuts1 += tnue.xsecWeight

  if tnue.vtxIsFiducial != 1:
    continue

  h_el_wCuts1.Fill(tnue.trueLepE, tnue.xsecWeight)
  h_nu_wCuts1.Fill(tnue.trueNuE, tnue.xsecWeight)
  n_scaled_passCuts1 += tnue.xsecWeight

  if tnue.vtxFracHitsOnCosmic >= 1.:
    continue

  h_el_wCuts2.Fill(tnue.trueLepE, tnue.xsecWeight)
  h_nu_wCuts2.Fill(tnue.trueNuE, tnue.xsecWeight)
  n_scaled_passCuts2 += tnue.xsecWeight

  maxMuScore = -20
  foundMuon = False
  for iT in range(tnue.nTracks):
    if tnue.trackIsSecondary[iT] == 1 or tnue.trackClassified[iT] != 1:
      continue
    if tnue.trackMuScore[iT] > maxMuScore:
      maxMuScore = tnue.trackMuScore[iT]
    if tnue.trackPID[iT] == 13:
      foundMuon = True

  if foundMuon:
    continue

  h_el_wCuts3.Fill(tnue.trueLepE, tnue.xsecWeight)
  h_nu_wCuts3.Fill(tnue.trueNuE, tnue.xsecWeight)
  n_scaled_passCuts3 += tnue.xsecWeight

  foundElectron = False
  elMaxQ = -1.
  elMaxQProc = -1
  elMaxQConf = -1.
  elMaxQCosTheta = -1.
  elMaxQVtxDist = -1.
  for iS in range(tnue.nShowers):
    if tnue.showerIsSecondary[iS] != 1 and tnue.showerClassified[iS] == 1 and tnue.showerPID[iS] == 11:
      foundElectron = True
      elConf = tnue.showerElScore[iS] - (tnue.showerPhScore[iS] + tnue.showerPiScore[iS])/2.
      if tnue.showerCharge[iS] > elMaxQ:
        elMaxQ = tnue.showerCharge[iS]
        elMaxQProc = tnue.showerProcess[iS]
        elMaxQConf = elConf
        elMaxQCosTheta = tnue.showerCosTheta[iS]
        elMaxQVtxDist = tnue.showerDistToVtx[iS]

  if not foundElectron:
    continue

  h_el_wCuts4.Fill(tnue.trueLepE, tnue.xsecWeight)
  h_nu_wCuts4.Fill(tnue.trueNuE, tnue.xsecWeight)
  n_scaled_passCuts4 += tnue.xsecWeight

  if elMaxQProc != 0:
    continue

  h_el_wCuts5.Fill(tnue.trueLepE, tnue.xsecWeight)
  h_nu_wCuts5.Fill(tnue.trueNuE, tnue.xsecWeight)
  n_scaled_passCuts5 += tnue.xsecWeight

  if maxMuScore >= -3.7:
    continue

  h_el_wCuts6.Fill(tnue.trueLepE, tnue.xsecWeight)
  h_nu_wCuts6.Fill(tnue.trueNuE, tnue.xsecWeight)
  n_scaled_passCuts6 += tnue.xsecWeight

  if elMaxQConf <= 7.1:
    continue

  h_el_wCuts7.Fill(tnue.trueLepE, tnue.xsecWeight)
  h_nu_wCuts7.Fill(tnue.trueNuE, tnue.xsecWeight)
  n_scaled_passCuts7 += tnue.xsecWeight

  #if elMaxQCosTheta <= -0.07:
  #  continue

  #h_el_wCuts9.Fill(tnue.trueLepE, tnue.xsecWeight)
  #h_nu_wCuts9.Fill(tnue.trueNuE, tnue.xsecWeight)
  #n_scaled_passCuts9 += tnue.xsecWeight

  #if elMaxQVtxDist >= 35.0:
  #  continue

  #h_el_wCuts10.Fill(tnue.trueLepE, tnue.xsecWeight)
  #h_nu_wCuts10.Fill(tnue.trueNuE, tnue.xsecWeight)
  #n_scaled_passCuts10 += tnue.xsecWeight



n_scaled_all *= targetPOT/tnuePOTsum
n_scaled_passCuts1 *= targetPOT/tnuePOTsum
n_scaled_passCuts2 *= targetPOT/tnuePOTsum
n_scaled_passCuts3 *= targetPOT/tnuePOTsum
n_scaled_passCuts4 *= targetPOT/tnuePOTsum
n_scaled_passCuts5 *= targetPOT/tnuePOTsum
n_scaled_passCuts6 *= targetPOT/tnuePOTsum
n_scaled_passCuts7 *= targetPOT/tnuePOTsum
#n_scaled_passCuts8 *= targetPOT/tnuePOTsum
#n_scaled_passCuts9 *= targetPOT/tnuePOTsum
#n_scaled_passCuts10 *= targetPOT/tnuePOTsum

print()
print("n_raw_all: ", n_raw_all)
print()
print("n_scaled_all: ", n_scaled_all)
print()
print("n_scaled_passCuts1: ", n_scaled_passCuts1)
print("n_scaled_passCuts2: ", n_scaled_passCuts2)
print("n_scaled_passCuts3: ", n_scaled_passCuts3)
print("n_scaled_passCuts4: ", n_scaled_passCuts4)
print("n_scaled_passCuts5: ", n_scaled_passCuts5)
print("n_scaled_passCuts6: ", n_scaled_passCuts6)
print("n_scaled_passCuts7: ", n_scaled_passCuts7)
#print("n_scaled_passCuts8: ", n_scaled_passCuts8)
#print("n_scaled_passCuts9: ", n_scaled_passCuts9)
#print("n_scaled_passCuts10: ", n_scaled_passCuts10)
print()
print("cut set 1 efficiency: %f"%(n_scaled_passCuts1/n_scaled_all))
print("cut set 2 efficiency: %f"%(n_scaled_passCuts2/n_scaled_all))
print("cut set 3 efficiency: %f"%(n_scaled_passCuts3/n_scaled_all))
print("cut set 4 efficiency: %f"%(n_scaled_passCuts4/n_scaled_all))
print("cut set 5 efficiency: %f"%(n_scaled_passCuts5/n_scaled_all))
print("cut set 6 efficiency: %f"%(n_scaled_passCuts6/n_scaled_all))
print("cut set 7 efficiency: %f"%(n_scaled_passCuts7/n_scaled_all))
#print("cut set 8 efficiency: %f"%(n_scaled_passCuts8/n_scaled_all))
#print("cut set 9 efficiency: %f"%(n_scaled_passCuts9/n_scaled_all))
#print("cut set 10 efficiency: %f"%(n_scaled_passCuts10/n_scaled_all))
print()


h_el_nCuts.Scale(targetPOT/tnuePOTsum)
h_el_wCuts1.Scale(targetPOT/tnuePOTsum)
h_el_wCuts2.Scale(targetPOT/tnuePOTsum)
h_el_wCuts3.Scale(targetPOT/tnuePOTsum)
h_el_wCuts4.Scale(targetPOT/tnuePOTsum)
h_el_wCuts5.Scale(targetPOT/tnuePOTsum)
h_el_wCuts6.Scale(targetPOT/tnuePOTsum)
h_el_wCuts7.Scale(targetPOT/tnuePOTsum)
#h_el_wCuts8.Scale(targetPOT/tnuePOTsum)
#h_el_wCuts9.Scale(targetPOT/tnuePOTsum)
#h_el_wCuts10.Scale(targetPOT/tnuePOTsum)
h_nu_nCuts.Scale(targetPOT/tnuePOTsum)
h_nu_wCuts1.Scale(targetPOT/tnuePOTsum)
h_nu_wCuts2.Scale(targetPOT/tnuePOTsum)
h_nu_wCuts3.Scale(targetPOT/tnuePOTsum)
h_nu_wCuts4.Scale(targetPOT/tnuePOTsum)
h_nu_wCuts5.Scale(targetPOT/tnuePOTsum)
h_nu_wCuts6.Scale(targetPOT/tnuePOTsum)
h_nu_wCuts7.Scale(targetPOT/tnuePOTsum)
#h_nu_wCuts8.Scale(targetPOT/tnuePOTsum)
#h_nu_wCuts9.Scale(targetPOT/tnuePOTsum)
#h_nu_wCuts10.Scale(targetPOT/tnuePOTsum)

h_el_eff1.Divide(h_el_wCuts1,h_el_nCuts,1,1,"B")
h_el_eff2.Divide(h_el_wCuts2,h_el_nCuts,1,1,"B")
h_el_eff3.Divide(h_el_wCuts3,h_el_nCuts,1,1,"B")
h_el_eff4.Divide(h_el_wCuts4,h_el_nCuts,1,1,"B")
h_el_eff5.Divide(h_el_wCuts5,h_el_nCuts,1,1,"B")
h_el_eff6.Divide(h_el_wCuts6,h_el_nCuts,1,1,"B")
h_el_eff7.Divide(h_el_wCuts7,h_el_nCuts,1,1,"B")
#h_el_eff8.Divide(h_el_wCuts8,h_el_nCuts,1,1,"B")
#h_el_eff9.Divide(h_el_wCuts9,h_el_nCuts,1,1,"B")
#h_el_eff10.Divide(h_el_wCuts10,h_el_nCuts,1,1,"B")
h_nu_eff1.Divide(h_nu_wCuts1,h_nu_nCuts,1,1,"B")
h_nu_eff2.Divide(h_nu_wCuts2,h_nu_nCuts,1,1,"B")
h_nu_eff3.Divide(h_nu_wCuts3,h_nu_nCuts,1,1,"B")
h_nu_eff4.Divide(h_nu_wCuts4,h_nu_nCuts,1,1,"B")
h_nu_eff5.Divide(h_nu_wCuts5,h_nu_nCuts,1,1,"B")
h_nu_eff6.Divide(h_nu_wCuts6,h_nu_nCuts,1,1,"B")
h_nu_eff7.Divide(h_nu_wCuts7,h_nu_nCuts,1,1,"B")
#h_nu_eff8.Divide(h_nu_wCuts8,h_nu_nCuts,1,1,"B")
#h_nu_eff9.Divide(h_nu_wCuts9,h_nu_nCuts,1,1,"B")
#h_nu_eff10.Divide(h_nu_wCuts10,h_nu_nCuts,1,1,"B")



outFile = rt.TFile(args.outfile, "RECREATE")

h_el_nCuts.Write()
h_el_wCuts1.Write()
h_el_wCuts2.Write()
h_el_wCuts3.Write()
h_el_wCuts4.Write()
h_el_wCuts5.Write()
h_el_wCuts6.Write()
h_el_wCuts7.Write()
#h_el_wCuts8.Write()
#h_el_wCuts9.Write()
#h_el_wCuts10.Write()
h_el_eff1.Write()
h_el_eff2.Write()
h_el_eff3.Write()
h_el_eff4.Write()
h_el_eff5.Write()
h_el_eff6.Write()
h_el_eff7.Write()
#h_el_eff8.Write()
#h_el_eff9.Write()
#h_el_eff10.Write()
h_nu_nCuts.Write()
h_nu_wCuts1.Write()
h_nu_wCuts2.Write()
h_nu_wCuts3.Write()
h_nu_wCuts4.Write()
h_nu_wCuts5.Write()
h_nu_wCuts6.Write()
h_nu_wCuts7.Write()
#h_nu_wCuts8.Write()
#h_nu_wCuts9.Write()
#h_nu_wCuts10.Write()
h_nu_eff1.Write()
h_nu_eff2.Write()
h_nu_eff3.Write()
h_nu_eff4.Write()
h_nu_eff5.Write()
h_nu_eff6.Write()
h_nu_eff7.Write()
#h_nu_eff8.Write()
#h_nu_eff9.Write()
#h_nu_eff10.Write()


#def getPlotYRange(h1, h2):
#  ymin = 9.
#  ymax = -9.
#  for i in range(1,h1.GetNBinsX()+1):
#    low1 = h1.GetBinContent(i) - h1.GetBinError(i)
#    high1 = h1.GetBinContent(i) + h1.GetBinError(i)
#    low2 = h2.GetBinContent(i) - h2.GetBinError(i)
#    high2 = h2.GetBinContent(i) + h2.GetBinError(i)
#    if low1 < ymin:
#      ymin = low1
#    if low2 < ymin:
#      ymin = low2
#    if high1 > ymax:
#      ymax = high1
#    if high2 > ymax:
#      ymax = high2
#  ymax = ymax + (ymax - ymin)*0.05
#  ymin = ymin - (ymax - ymin)*0.05
#  if ymax > 1:
#    ymax = 1.003
#  if ymin < 0:
#    ymin = 0.
#  return ymin, ymax


h_el_eff1.GetYaxis().SetRangeUser(0,1.003)
h_el_eff2.GetYaxis().SetRangeUser(0,1.003)
h_el_eff3.GetYaxis().SetRangeUser(0,1.003)
h_el_eff4.GetYaxis().SetRangeUser(0,1.003)
h_el_eff5.GetYaxis().SetRangeUser(0,1.003)
h_el_eff6.GetYaxis().SetRangeUser(0,1.003)
h_el_eff7.GetYaxis().SetRangeUser(0,1.003)
#h_el_eff8.GetYaxis().SetRangeUser(0,1.003)
#h_el_eff9.GetYaxis().SetRangeUser(0,1.003)
#h_el_eff10.GetYaxis().SetRangeUser(0,1.003)

h_nu_eff1.GetYaxis().SetRangeUser(0,1.003)
h_nu_eff2.GetYaxis().SetRangeUser(0,1.003)
h_nu_eff3.GetYaxis().SetRangeUser(0,1.003)
h_nu_eff4.GetYaxis().SetRangeUser(0,1.003)
h_nu_eff5.GetYaxis().SetRangeUser(0,1.003)
h_nu_eff6.GetYaxis().SetRangeUser(0,1.003)
h_nu_eff7.GetYaxis().SetRangeUser(0,1.003)
#h_nu_eff8.GetYaxis().SetRangeUser(0,1.003)
#h_nu_eff9.GetYaxis().SetRangeUser(0,1.003)
#h_nu_eff10.GetYaxis().SetRangeUser(0,1.003)


cnv_el_E_1 = rt.TCanvas("cnv_el_E_1","cnv_el_E_1")
cnv_el_E_1.SetGrid()
h_el_nCuts.Draw("EHIST")
h_el_wCuts1.Draw("EHISTSAME")
leg_E_1 = rt.TLegend(0.7,0.7,0.9,0.9)
leg_E_1.AddEntry(h_el_nCuts, "all true signal (%.2f)"%h_el_nCuts.Integral(), "l")
leg_E_1.AddEntry(h_el_wCuts1, "cut set 1 (%.2f)"%h_el_wCuts1.Integral(), "l")
leg_E_1.Draw()
cnv_el_E_1.Write()

cnv_el_eff_1 = rt.TCanvas("cnv_el_eff_1","cnv_el_eff_1")
cnv_el_eff_1.SetGrid()
h_el_eff1.Draw("E")
leg_eff_1 = rt.TLegend(0.7,0.7,0.9,0.9)
leg_eff_1.AddEntry(h_el_eff1, "cut set 1", "l")
leg_eff_1.Draw()
cnv_el_eff_1.Write()


cnv_el_E_2 = rt.TCanvas("cnv_el_E_2","cnv_el_E_2")
cnv_el_E_2.SetGrid()
h_el_nCuts.Draw("EHIST")
h_el_wCuts1.Draw("EHISTSAME")
h_el_wCuts2.Draw("EHISTSAME")
leg_E_2 = rt.TLegend(0.7,0.7,0.9,0.9)
leg_E_2.AddEntry(h_el_nCuts, "all true signal (%.2f)"%h_el_nCuts.Integral(), "l")
leg_E_2.AddEntry(h_el_wCuts1, "cut set 1 (%.2f)"%h_el_wCuts1.Integral(), "l")
leg_E_2.AddEntry(h_el_wCuts2, "cut set 2 (%.2f)"%h_el_wCuts2.Integral(), "l")
leg_E_2.Draw()
cnv_el_E_2.Write()

cnv_el_eff_2 = rt.TCanvas("cnv_el_eff_2","cnv_el_eff_2")
cnv_el_eff_2.SetGrid()
h_el_eff1.Draw("E")
h_el_eff2.Draw("ESAME")
leg_eff_2 = rt.TLegend(0.7,0.7,0.9,0.9)
leg_eff_2.AddEntry(h_el_eff1, "cut set 1", "l")
leg_eff_2.AddEntry(h_el_eff2, "cut set 2", "l")
leg_eff_2.Draw()
cnv_el_eff_2.Write()


cnv_el_E_3 = rt.TCanvas("cnv_el_E_3","cnv_el_E_3")
cnv_el_E_3.SetGrid()
h_el_nCuts.Draw("EHIST")
h_el_wCuts2.Draw("EHISTSAME")
h_el_wCuts3.Draw("EHISTSAME")
leg_E_3 = rt.TLegend(0.7,0.7,0.9,0.9)
leg_E_3.AddEntry(h_el_nCuts, "all true signal (%.2f)"%h_el_nCuts.Integral(), "l")
leg_E_3.AddEntry(h_el_wCuts2, "cut set 2 (%.2f)"%h_el_wCuts2.Integral(), "l")
leg_E_3.AddEntry(h_el_wCuts3, "cut set 3 (%.2f)"%h_el_wCuts3.Integral(), "l")
leg_E_3.Draw()
cnv_el_E_3.Write()

cnv_el_eff_3 = rt.TCanvas("cnv_el_eff_3","cnv_el_eff_3")
cnv_el_eff_3.SetGrid()
h_el_eff2.Draw("E")
h_el_eff3.Draw("ESAME")
leg_eff_3 = rt.TLegend(0.7,0.7,0.9,0.9)
leg_eff_3.AddEntry(h_el_eff2, "cut set 2", "l")
leg_eff_3.AddEntry(h_el_eff3, "cut set 3", "l")
leg_eff_3.Draw()
cnv_el_eff_3.Write()


cnv_el_E_4 = rt.TCanvas("cnv_el_E_4","cnv_el_E_4")
cnv_el_E_4.SetGrid()
h_el_nCuts.Draw("EHIST")
h_el_wCuts3.Draw("EHISTSAME")
h_el_wCuts4.Draw("EHISTSAME")
leg_E_4 = rt.TLegend(0.7,0.7,0.9,0.9)
leg_E_4.AddEntry(h_el_nCuts, "all true signal (%.2f)"%h_el_nCuts.Integral(), "l")
leg_E_4.AddEntry(h_el_wCuts3, "cut set 3 (%.2f)"%h_el_wCuts3.Integral(), "l")
leg_E_4.AddEntry(h_el_wCuts4, "cut set 4 (%.2f)"%h_el_wCuts4.Integral(), "l")
leg_E_4.Draw()
cnv_el_E_4.Write()

cnv_el_eff_4 = rt.TCanvas("cnv_el_eff_4","cnv_el_eff_4")
cnv_el_eff_4.SetGrid()
h_el_eff3.Draw("E")
h_el_eff4.Draw("ESAME")
leg_eff_4 = rt.TLegend(0.7,0.7,0.9,0.9)
leg_eff_4.AddEntry(h_el_eff3, "cut set 3", "l")
leg_eff_4.AddEntry(h_el_eff4, "cut set 4", "l")
leg_eff_4.Draw()
cnv_el_eff_4.Write()


cnv_el_E_5 = rt.TCanvas("cnv_el_E_5","cnv_el_E_5")
cnv_el_E_5.SetGrid()
h_el_nCuts.Draw("EHIST")
h_el_wCuts4.Draw("EHISTSAME")
h_el_wCuts5.Draw("EHISTSAME")
leg_E_5 = rt.TLegend(0.7,0.7,0.9,0.9)
leg_E_5.AddEntry(h_el_nCuts, "all true signal (%.2f)"%h_el_nCuts.Integral(), "l")
leg_E_5.AddEntry(h_el_wCuts4, "cut set 4 (%.2f)"%h_el_wCuts4.Integral(), "l")
leg_E_5.AddEntry(h_el_wCuts5, "cut set 5 (%.2f)"%h_el_wCuts5.Integral(), "l")
leg_E_5.Draw()
cnv_el_E_5.Write()

cnv_el_eff_5 = rt.TCanvas("cnv_el_eff_5","cnv_el_eff_5")
cnv_el_eff_5.SetGrid()
h_el_eff4.Draw("E")
h_el_eff5.Draw("ESAME")
leg_eff_5 = rt.TLegend(0.7,0.7,0.9,0.9)
leg_eff_5.AddEntry(h_el_eff4, "cut set 4", "l")
leg_eff_5.AddEntry(h_el_eff5, "cut set 5", "l")
leg_eff_5.Draw()
cnv_el_eff_5.Write()


cnv_el_E_6 = rt.TCanvas("cnv_el_E_6","cnv_el_E_6")
cnv_el_E_6.SetGrid()
h_el_nCuts.Draw("EHIST")
h_el_wCuts5.Draw("EHISTSAME")
h_el_wCuts6.Draw("EHISTSAME")
leg_E_6 = rt.TLegend(0.7,0.7,0.9,0.9)
leg_E_6.AddEntry(h_el_nCuts, "all true signal (%.2f)"%h_el_nCuts.Integral(), "l")
leg_E_6.AddEntry(h_el_wCuts5, "cut set 5 (%.2f)"%h_el_wCuts5.Integral(), "l")
leg_E_6.AddEntry(h_el_wCuts6, "cut set 6 (%.2f)"%h_el_wCuts6.Integral(), "l")
leg_E_6.Draw()
cnv_el_E_6.Write()

cnv_el_eff_6 = rt.TCanvas("cnv_el_eff_6","cnv_el_eff_6")
cnv_el_eff_6.SetGrid()
h_el_eff5.Draw("E")
h_el_eff6.Draw("ESAME")
leg_eff_6 = rt.TLegend(0.7,0.7,0.9,0.9)
leg_eff_6.AddEntry(h_el_eff5, "cut set 5", "l")
leg_eff_6.AddEntry(h_el_eff6, "cut set 6", "l")
leg_eff_6.Draw()
cnv_el_eff_6.Write()


cnv_el_E_7 = rt.TCanvas("cnv_el_E_7","cnv_el_E_7")
cnv_el_E_7.SetGrid()
h_el_nCuts.Draw("EHIST")
h_el_wCuts6.Draw("EHISTSAME")
h_el_wCuts7.Draw("EHISTSAME")
leg_E_7 = rt.TLegend(0.7,0.7,0.9,0.9)
leg_E_7.AddEntry(h_el_nCuts, "all true signal (%.2f)"%h_el_nCuts.Integral(), "l")
leg_E_7.AddEntry(h_el_wCuts6, "cut set 6 (%.2f)"%h_el_wCuts6.Integral(), "l")
leg_E_7.AddEntry(h_el_wCuts7, "all cuts (%.2f)"%h_el_wCuts7.Integral(), "l")
leg_E_7.Draw()
cnv_el_E_7.Write()

cnv_el_eff_7 = rt.TCanvas("cnv_el_eff_7","cnv_el_eff_7")
cnv_el_eff_7.SetGrid()
h_el_eff6.Draw("E")
h_el_eff7.Draw("ESAME")
leg_eff_7 = rt.TLegend(0.7,0.7,0.9,0.9)
leg_eff_7.AddEntry(h_el_eff6, "cut set 6", "l")
leg_eff_7.AddEntry(h_el_eff7, "all cuts", "l")
leg_eff_7.Draw()
cnv_el_eff_7.Write()


#cnv_el_E_8 = rt.TCanvas("cnv_el_E_8","cnv_el_E_8")
#cnv_el_E_8.SetGrid()
#h_el_nCuts.Draw("EHIST")
#h_el_wCuts7.Draw("EHISTSAME")
#h_el_wCuts8.Draw("EHISTSAME")
#leg_E_8 = rt.TLegend(0.7,0.7,0.9,0.9)
#leg_E_8.AddEntry(h_el_nCuts, "all true signal (%.2f)"%h_el_nCuts.Integral(), "l")
#leg_E_8.AddEntry(h_el_wCuts7, "cut set 7 (%.2f)"%h_el_wCuts7.Integral(), "l")
#leg_E_8.AddEntry(h_el_wCuts8, "cut set 8 (%.2f)"%h_el_wCuts8.Integral(), "l")
#leg_E_8.Draw()
#cnv_el_E_8.Write()
#
#cnv_el_eff_8 = rt.TCanvas("cnv_el_eff_8","cnv_el_eff_8")
#cnv_el_eff_8.SetGrid()
#h_el_eff7.Draw("E")
#h_el_eff8.Draw("ESAME")
#leg_eff_8 = rt.TLegend(0.7,0.7,0.9,0.9)
#leg_eff_8.AddEntry(h_el_eff7, "cut set 7", "l")
#leg_eff_8.AddEntry(h_el_eff8, "cut set 8", "l")
#leg_eff_8.Draw()
#cnv_el_eff_8.Write()
#
#
#cnv_el_E_9 = rt.TCanvas("cnv_el_E_9","cnv_el_E_9")
#cnv_el_E_9.SetGrid()
#h_el_nCuts.Draw("EHIST")
#h_el_wCuts8.Draw("EHISTSAME")
#h_el_wCuts9.Draw("EHISTSAME")
#leg_E_9 = rt.TLegend(0.7,0.7,0.9,0.9)
#leg_E_9.AddEntry(h_el_nCuts, "all true signal (%.2f)"%h_el_nCuts.Integral(), "l")
#leg_E_9.AddEntry(h_el_wCuts8, "cut set 8 (%.2f)"%h_el_wCuts8.Integral(), "l")
#leg_E_9.AddEntry(h_el_wCuts9, "cut set 9 (%.2f)"%h_el_wCuts9.Integral(), "l")
#leg_E_9.Draw()
#cnv_el_E_9.Write()
#
#cnv_el_eff_9 = rt.TCanvas("cnv_el_eff_9","cnv_el_eff_9")
#cnv_el_eff_9.SetGrid()
#h_el_eff8.Draw("E")
#h_el_eff9.Draw("ESAME")
#leg_eff_9 = rt.TLegend(0.7,0.7,0.9,0.9)
#leg_eff_9.AddEntry(h_el_eff8, "cut set 8", "l")
#leg_eff_9.AddEntry(h_el_eff9, "cut set 9", "l")
#leg_eff_9.Draw()
#cnv_el_eff_9.Write()
#
#
#cnv_el_E_10 = rt.TCanvas("cnv_el_E_10","cnv_el_E_10")
#cnv_el_E_10.SetGrid()
#h_el_nCuts.Draw("EHIST")
#h_el_wCuts9.Draw("EHISTSAME")
#h_el_wCuts10.Draw("EHISTSAME")
#leg_E_10 = rt.TLegend(0.7,0.7,0.9,0.9)
#leg_E_10.AddEntry(h_el_nCuts, "all true signal (%.2f)"%h_el_nCuts.Integral(), "l")
#leg_E_10.AddEntry(h_el_wCuts9, "cut set 9 (%.2f)"%h_el_wCuts9.Integral(), "l")
#leg_E_10.AddEntry(h_el_wCuts10, "cut set 10 (%.2f)"%h_el_wCuts10.Integral(), "l")
#leg_E_10.Draw()
#cnv_el_E_10.Write()
#
#cnv_el_eff_10 = rt.TCanvas("cnv_el_eff_10","cnv_el_eff_10")
#cnv_el_eff_10.SetGrid()
#h_el_eff9.Draw("E")
#h_el_eff10.Draw("ESAME")
#leg_eff_10 = rt.TLegend(0.7,0.7,0.9,0.9)
#leg_eff_10.AddEntry(h_el_eff9, "cut set 9", "l")
#leg_eff_10.AddEntry(h_el_eff10, "cut set 10", "l")
#leg_eff_10.Draw()
#cnv_el_eff_10.Write()


cnv_nu_E_1 = rt.TCanvas("cnv_nu_E_1","cnv_nu_E_1")
cnv_nu_E_1.SetGrid()
h_nu_nCuts.Draw("EHIST")
h_nu_wCuts1.Draw("EHISTSAME")
leg_E_1 = rt.TLegend(0.7,0.7,0.9,0.9)
leg_E_1.AddEntry(h_nu_nCuts, "all true signal (%.2f)"%h_nu_nCuts.Integral(), "l")
leg_E_1.AddEntry(h_nu_wCuts1, "cut set 1 (%.2f)"%h_nu_wCuts1.Integral(), "l")
leg_E_1.Draw()
cnv_nu_E_1.Write()

cnv_nu_eff_1 = rt.TCanvas("cnv_nu_eff_1","cnv_nu_eff_1")
cnv_nu_eff_1.SetGrid()
h_nu_eff1.Draw("E")
leg_eff_1 = rt.TLegend(0.7,0.7,0.9,0.9)
leg_eff_1.AddEntry(h_nu_eff1, "cut set 1", "l")
leg_eff_1.Draw()
cnv_nu_eff_1.Write()


cnv_nu_E_2 = rt.TCanvas("cnv_nu_E_2","cnv_nu_E_2")
cnv_nu_E_2.SetGrid()
h_nu_nCuts.Draw("EHIST")
h_nu_wCuts1.Draw("EHISTSAME")
h_nu_wCuts2.Draw("EHISTSAME")
leg_E_2 = rt.TLegend(0.7,0.7,0.9,0.9)
leg_E_2.AddEntry(h_nu_nCuts, "all true signal (%.2f)"%h_nu_nCuts.Integral(), "l")
leg_E_2.AddEntry(h_nu_wCuts1, "cut set 1 (%.2f)"%h_nu_wCuts1.Integral(), "l")
leg_E_2.AddEntry(h_nu_wCuts2, "cut set 2 (%.2f)"%h_nu_wCuts2.Integral(), "l")
leg_E_2.Draw()
cnv_nu_E_2.Write()

cnv_nu_eff_2 = rt.TCanvas("cnv_nu_eff_2","cnv_nu_eff_2")
cnv_nu_eff_2.SetGrid()
h_nu_eff1.Draw("E")
h_nu_eff2.Draw("ESAME")
leg_eff_2 = rt.TLegend(0.7,0.7,0.9,0.9)
leg_eff_2.AddEntry(h_nu_eff1, "cut set 1", "l")
leg_eff_2.AddEntry(h_nu_eff2, "cut set 2", "l")
leg_eff_2.Draw()
cnv_nu_eff_2.Write()


cnv_nu_E_3 = rt.TCanvas("cnv_nu_E_3","cnv_nu_E_3")
cnv_nu_E_3.SetGrid()
h_nu_nCuts.Draw("EHIST")
h_nu_wCuts2.Draw("EHISTSAME")
h_nu_wCuts3.Draw("EHISTSAME")
leg_E_3 = rt.TLegend(0.7,0.7,0.9,0.9)
leg_E_3.AddEntry(h_nu_nCuts, "all true signal (%.2f)"%h_nu_nCuts.Integral(), "l")
leg_E_3.AddEntry(h_nu_wCuts2, "cut set 2 (%.2f)"%h_nu_wCuts2.Integral(), "l")
leg_E_3.AddEntry(h_nu_wCuts3, "cut set 3 (%.2f)"%h_nu_wCuts3.Integral(), "l")
leg_E_3.Draw()
cnv_nu_E_3.Write()

cnv_nu_eff_3 = rt.TCanvas("cnv_nu_eff_3","cnv_nu_eff_3")
cnv_nu_eff_3.SetGrid()
h_nu_eff2.Draw("E")
h_nu_eff3.Draw("ESAME")
leg_eff_3 = rt.TLegend(0.7,0.7,0.9,0.9)
leg_eff_3.AddEntry(h_nu_eff2, "cut set 2", "l")
leg_eff_3.AddEntry(h_nu_eff3, "cut set 3", "l")
leg_eff_3.Draw()
cnv_nu_eff_3.Write()


cnv_nu_E_4 = rt.TCanvas("cnv_nu_E_4","cnv_nu_E_4")
cnv_nu_E_4.SetGrid()
h_nu_nCuts.Draw("EHIST")
h_nu_wCuts3.Draw("EHISTSAME")
h_nu_wCuts4.Draw("EHISTSAME")
leg_E_4 = rt.TLegend(0.7,0.7,0.9,0.9)
leg_E_4.AddEntry(h_nu_nCuts, "all true signal (%.2f)"%h_nu_nCuts.Integral(), "l")
leg_E_4.AddEntry(h_nu_wCuts3, "cut set 3 (%.2f)"%h_nu_wCuts3.Integral(), "l")
leg_E_4.AddEntry(h_nu_wCuts4, "cut set 4 (%.2f)"%h_nu_wCuts4.Integral(), "l")
leg_E_4.Draw()
cnv_nu_E_4.Write()

cnv_nu_eff_4 = rt.TCanvas("cnv_nu_eff_4","cnv_nu_eff_4")
cnv_nu_eff_4.SetGrid()
h_nu_eff3.Draw("E")
h_nu_eff4.Draw("ESAME")
leg_eff_4 = rt.TLegend(0.7,0.7,0.9,0.9)
leg_eff_4.AddEntry(h_nu_eff3, "cut set 3", "l")
leg_eff_4.AddEntry(h_nu_eff4, "cut set 4", "l")
leg_eff_4.Draw()
cnv_nu_eff_4.Write()


cnv_nu_E_5 = rt.TCanvas("cnv_nu_E_5","cnv_nu_E_5")
cnv_nu_E_5.SetGrid()
h_nu_nCuts.Draw("EHIST")
h_nu_wCuts4.Draw("EHISTSAME")
h_nu_wCuts5.Draw("EHISTSAME")
leg_E_5 = rt.TLegend(0.7,0.7,0.9,0.9)
leg_E_5.AddEntry(h_nu_nCuts, "all true signal (%.2f)"%h_nu_nCuts.Integral(), "l")
leg_E_5.AddEntry(h_nu_wCuts4, "cut set 4 (%.2f)"%h_nu_wCuts4.Integral(), "l")
leg_E_5.AddEntry(h_nu_wCuts5, "cut set 5 (%.2f)"%h_nu_wCuts5.Integral(), "l")
leg_E_5.Draw()
cnv_nu_E_5.Write()

cnv_nu_eff_5 = rt.TCanvas("cnv_nu_eff_5","cnv_nu_eff_5")
cnv_nu_eff_5.SetGrid()
h_nu_eff4.Draw("E")
h_nu_eff5.Draw("ESAME")
leg_eff_5 = rt.TLegend(0.7,0.7,0.9,0.9)
leg_eff_5.AddEntry(h_nu_eff4, "cut set 4", "l")
leg_eff_5.AddEntry(h_nu_eff5, "cut set 5", "l")
leg_eff_5.Draw()
cnv_nu_eff_5.Write()


cnv_nu_E_6 = rt.TCanvas("cnv_nu_E_6","cnv_nu_E_6")
cnv_nu_E_6.SetGrid()
h_nu_nCuts.Draw("EHIST")
h_nu_wCuts5.Draw("EHISTSAME")
h_nu_wCuts6.Draw("EHISTSAME")
leg_E_6 = rt.TLegend(0.7,0.7,0.9,0.9)
leg_E_6.AddEntry(h_nu_nCuts, "all true signal (%.2f)"%h_nu_nCuts.Integral(), "l")
leg_E_6.AddEntry(h_nu_wCuts5, "cut set 5 (%.2f)"%h_nu_wCuts5.Integral(), "l")
leg_E_6.AddEntry(h_nu_wCuts6, "cut set 6 (%.2f)"%h_nu_wCuts6.Integral(), "l")
leg_E_6.Draw()
cnv_nu_E_6.Write()

cnv_nu_eff_6 = rt.TCanvas("cnv_nu_eff_6","cnv_nu_eff_6")
cnv_nu_eff_6.SetGrid()
h_nu_eff5.Draw("E")
h_nu_eff6.Draw("ESAME")
leg_eff_6 = rt.TLegend(0.7,0.7,0.9,0.9)
leg_eff_6.AddEntry(h_nu_eff5, "cut set 5", "l")
leg_eff_6.AddEntry(h_nu_eff6, "cut set 6", "l")
leg_eff_6.Draw()
cnv_nu_eff_6.Write()


cnv_nu_E_7 = rt.TCanvas("cnv_nu_E_7","cnv_nu_E_7")
cnv_nu_E_7.SetGrid()
h_nu_nCuts.Draw("EHIST")
h_nu_wCuts6.Draw("EHISTSAME")
h_nu_wCuts7.Draw("EHISTSAME")
leg_E_7 = rt.TLegend(0.7,0.7,0.9,0.9)
leg_E_7.AddEntry(h_nu_nCuts, "all true signal (%.2f)"%h_nu_nCuts.Integral(), "l")
leg_E_7.AddEntry(h_nu_wCuts6, "cut set 6 (%.2f)"%h_nu_wCuts6.Integral(), "l")
leg_E_7.AddEntry(h_nu_wCuts7, "all cuts (%.2f)"%h_nu_wCuts7.Integral(), "l")
leg_E_7.Draw()
cnv_nu_E_7.Write()

cnv_nu_eff_7 = rt.TCanvas("cnv_nu_eff_7","cnv_nu_eff_7")
cnv_nu_eff_7.SetGrid()
h_nu_eff6.Draw("E")
h_nu_eff7.Draw("ESAME")
leg_eff_7 = rt.TLegend(0.7,0.7,0.9,0.9)
leg_eff_7.AddEntry(h_nu_eff6, "cut set 6", "l")
leg_eff_7.AddEntry(h_nu_eff7, "all cuts", "l")
leg_eff_7.Draw()
cnv_nu_eff_7.Write()


#cnv_nu_E_8 = rt.TCanvas("cnv_nu_E_8","cnv_nu_E_8")
#cnv_nu_E_8.SetGrid()
#h_nu_nCuts.Draw("EHIST")
#h_nu_wCuts7.Draw("EHISTSAME")
#h_nu_wCuts8.Draw("EHISTSAME")
#leg_E_8 = rt.TLegend(0.7,0.7,0.9,0.9)
#leg_E_8.AddEntry(h_nu_nCuts, "all true signal (%.2f)"%h_nu_nCuts.Integral(), "l")
#leg_E_8.AddEntry(h_nu_wCuts7, "cut set 7 (%.2f)"%h_nu_wCuts7.Integral(), "l")
#leg_E_8.AddEntry(h_nu_wCuts8, "cut set 8 (%.2f)"%h_nu_wCuts8.Integral(), "l")
#leg_E_8.Draw()
#cnv_nu_E_8.Write()
#
#cnv_nu_eff_8 = rt.TCanvas("cnv_nu_eff_8","cnv_nu_eff_8")
#cnv_nu_eff_8.SetGrid()
#h_nu_eff7.Draw("E")
#h_nu_eff8.Draw("ESAME")
#leg_eff_8 = rt.TLegend(0.7,0.7,0.9,0.9)
#leg_eff_8.AddEntry(h_nu_eff7, "cut set 7", "l")
#leg_eff_8.AddEntry(h_nu_eff8, "cut set 8", "l")
#leg_eff_8.Draw()
#cnv_nu_eff_8.Write()
#
#
#cnv_nu_E_9 = rt.TCanvas("cnv_nu_E_9","cnv_nu_E_9")
#cnv_nu_E_9.SetGrid()
#h_nu_nCuts.Draw("EHIST")
#h_nu_wCuts8.Draw("EHISTSAME")
#h_nu_wCuts9.Draw("EHISTSAME")
#leg_E_9 = rt.TLegend(0.7,0.7,0.9,0.9)
#leg_E_9.AddEntry(h_nu_nCuts, "all true signal (%.2f)"%h_nu_nCuts.Integral(), "l")
#leg_E_9.AddEntry(h_nu_wCuts8, "cut set 8 (%.2f)"%h_nu_wCuts8.Integral(), "l")
#leg_E_9.AddEntry(h_nu_wCuts9, "cut set 9 (%.2f)"%h_nu_wCuts9.Integral(), "l")
#leg_E_9.Draw()
#cnv_nu_E_9.Write()
#
#cnv_nu_eff_9 = rt.TCanvas("cnv_nu_eff_9","cnv_nu_eff_9")
#cnv_nu_eff_9.SetGrid()
#h_nu_eff8.Draw("E")
#h_nu_eff9.Draw("ESAME")
#leg_eff_9 = rt.TLegend(0.7,0.7,0.9,0.9)
#leg_eff_9.AddEntry(h_nu_eff8, "cut set 8", "l")
#leg_eff_9.AddEntry(h_nu_eff9, "cut set 9", "l")
#leg_eff_9.Draw()
#cnv_nu_eff_9.Write()
#
#
#cnv_nu_E_10 = rt.TCanvas("cnv_nu_E_10","cnv_nu_E_10")
#cnv_nu_E_10.SetGrid()
#h_nu_nCuts.Draw("EHIST")
#h_nu_wCuts9.Draw("EHISTSAME")
#h_nu_wCuts10.Draw("EHISTSAME")
#leg_E_10 = rt.TLegend(0.7,0.7,0.9,0.9)
#leg_E_10.AddEntry(h_nu_nCuts, "all true signal (%.2f)"%h_nu_nCuts.Integral(), "l")
#leg_E_10.AddEntry(h_nu_wCuts9, "cut set 9 (%.2f)"%h_nu_wCuts9.Integral(), "l")
#leg_E_10.AddEntry(h_nu_wCuts10, "cut set 10 (%.2f)"%h_nu_wCuts10.Integral(), "l")
#leg_E_10.Draw()
#cnv_nu_E_10.Write()
#
#cnv_nu_eff_10 = rt.TCanvas("cnv_nu_eff_10","cnv_nu_eff_10")
#cnv_nu_eff_10.SetGrid()
#h_nu_eff9.Draw("E")
#h_nu_eff10.Draw("ESAME")
#leg_eff_10 = rt.TLegend(0.7,0.7,0.9,0.9)
#leg_eff_10.AddEntry(h_nu_eff9, "cut set 9", "l")
#leg_eff_10.AddEntry(h_nu_eff10, "cut set 10", "l")
#leg_eff_10.Draw()
#cnv_nu_eff_10.Write()




