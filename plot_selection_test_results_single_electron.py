
import sys, argparse
import numpy as np
import ROOT as rt

from math import sqrt, isinf
from helpers.plotting_functions import sortHists
from helpers.larflowreco_ana_funcs import isFiducialWC


parser = argparse.ArgumentParser("Plot Selection Test Results")
#parser.add_argument("-fnu", "--bnbnu_file", type=str, default="/home/matthew/microboone/tufts/files/data/gen2_ntuples/dlgen2_ntuple_reco_v2me05_gen2ntuple_v0_run3b_bnb_nu_overlay.root", help="bnb nu input file")
parser.add_argument("-fnue", "--bnbnue_file", type=str, default="/home/matthew/microboone/tufts/files/data/gen2_ntuples/dlgen2_ntuple_reco_v2me05_gen2ntuple_v0_run3b_bnb_intrinsic_nue_overlay.root", help="bnb nue input file")
#parser.add_argument("-fext", "--extbnb_file", type=str, default="selection_output/prepare_selection_test_output/prepare_selection_test_reco_v2me05_gen2val_v22_extbnb_file.root", help="extbnb input file")
#parser.add_argument("-fdata", "--data_file", type=str, default="selection_output/prepare_selection_test_output/prepare_selection_test_reco_v2me05_gen2val_v22_bnb5e19_file.root", help="bnb data input file")
parser.add_argument("-o", "--outfile", type=str, default="selection_output/plot_selection_test_results_single_electron_output.root", help="output root file name")
#parser.add_argument("--recoEOverflow", help="plot overflow bin for final recoE selection plot", action="store_true")
args = parser.parse_args()

#rt.gROOT.SetBatch(True)
rt.TH1.SetDefaultSumw2(rt.kTRUE)
rt.gStyle.SetOptStat(0)

fnue = rt.TFile(args.bnbnue_file)
tnue = fnue.Get("EventTree")
tnuePOT = fnue.Get("potTree")

targetPOT = 6.67e+20

tnuePOTsum = 0.
for i in range(tnuePOT.GetEntries()):
  tnuePOT.GetEntry(i)
  tnuePOTsum = tnuePOTsum + tnuePOT.totGoodPOT


def configureHist(hist, color, ylabel):
  hist.GetXaxis().SetTitle("true electron energy (GeV)")
  hist.GetYaxis().SetTitle(ylabel)
  hist.SetLineColor(color)
  hist.SetLineWidth(2)
  return hist

h_nCuts = rt.TH1F("h_nCuts","Single Electron CCnue Electron Energy",40,0,4)
h_wCuts1 = rt.TH1F("h_wCuts1","Single Electron CCnue Electron Energy",40,0,4)
h_wCuts2 = rt.TH1F("h_wCuts2","Single Electron CCnue Electron Energy",40,0,4)
h_wCuts3 = rt.TH1F("h_wCuts3","Single Electron CCnue Electron Energy",40,0,4)
h_wCuts4 = rt.TH1F("h_wCuts4","Single Electron CCnue Electron Energy",40,0,4)
h_wCuts5 = rt.TH1F("h_wCuts5","Single Electron CCnue Electron Energy",40,0,4)
h_wCuts6 = rt.TH1F("h_wCuts6","Single Electron CCnue Electron Energy",40,0,4)
h_wCuts7 = rt.TH1F("h_wCuts7","Single Electron CCnue Electron Energy",40,0,4)
h_eff1 = rt.TH1F("h_eff1","Single Electron CCnue Efficiency",40,0,4)
h_eff2 = rt.TH1F("h_eff2","Single Electron CCnue Efficiency",40,0,4)
h_eff3 = rt.TH1F("h_eff3","Single Electron CCnue Efficiency",40,0,4)
h_eff4 = rt.TH1F("h_eff4","Single Electron CCnue Efficiency",40,0,4)
h_eff5 = rt.TH1F("h_eff5","Single Electron CCnue Efficiency",40,0,4)
h_eff6 = rt.TH1F("h_eff6","Single Electron CCnue Efficiency",40,0,4)
h_eff7 = rt.TH1F("h_eff7","Single Electron CCnue Efficiency",40,0,4)
h_nCuts = configureHist(h_nCuts, 1, "events per 6.67e+20 POT")
h_wCuts1 = configureHist(h_wCuts1, 2, "events per 6.67e+20 POT")
h_wCuts2 = configureHist(h_wCuts2, 3, "events per 6.67e+20 POT")
h_wCuts3 = configureHist(h_wCuts3, 4, "events per 6.67e+20 POT")
h_wCuts4 = configureHist(h_wCuts4, 8, "events per 6.67e+20 POT")
h_wCuts5 = configureHist(h_wCuts5, 6, "events per 6.67e+20 POT")
h_wCuts6 = configureHist(h_wCuts6, 40, "events per 6.67e+20 POT")
h_wCuts7 = configureHist(h_wCuts7, 46, "events per 6.67e+20 POT")
h_eff1 = configureHist(h_eff1, 2, "efficiency")
h_eff2 = configureHist(h_eff2, 3, "efficiency")
h_eff3 = configureHist(h_eff3, 4, "efficiency")
h_eff4 = configureHist(h_eff4, 8, "efficiency")
h_eff5 = configureHist(h_eff5, 6, "efficiency")
h_eff6 = configureHist(h_eff6, 40, "efficiency")
h_eff7 = configureHist(h_eff7, 46, "efficiency")

h2d_nCuts = rt.TH2F("h2d_nCuts","Single Electron CCnue Events (all true)",13,0,4,13,0,4)
h2d_wCuts7 = rt.TH2F("h2d_wCuts7","Single Electron CCnue Events (with reco cuts)",13,0,4,13,0,4)
h2d_nCuts_raw = rt.TH2F("h2d_nCuts_raw","Single Electron CCnue Raw Event Counts (all true)",13,0,4,13,0,4)
h2d_wCuts7_raw = rt.TH2F("h2d_wCuts7_raw","Single Electron CCnue Raw Event Counts (with reco cuts)",13,0,4,13,0,4)
h2d_eff7 = rt.TH2F("h2d_eff7","Single Electron CCnue Efficiency",13,0,4,13,0,4)
h2d_nCuts.GetXaxis().SetTitle("true electron energy (GeV)")
h2d_nCuts.GetYaxis().SetTitle("true neutrino energy (GeV)")
h2d_wCuts7.GetXaxis().SetTitle("true electron energy (GeV)")
h2d_wCuts7.GetYaxis().SetTitle("true neutrino energy (GeV)")
h2d_eff7.GetXaxis().SetTitle("true electron energy (GeV)")
h2d_eff7.GetYaxis().SetTitle("true neutrino energy (GeV)")

n_raw_all = 0
n_runs1to3_all = 0.
n_runs1to3_passCuts1 = 0.
n_runs1to3_passCuts2 = 0.
n_runs1to3_passCuts3 = 0.
n_runs1to3_passCuts4 = 0.
n_runs1to3_passCuts5 = 0.
n_runs1to3_passCuts6 = 0.
n_runs1to3_passCuts7 = 0.


extra_particle_dict = {}

for i in range(tnue.GetEntries()):

  tnue.GetEntry(i)

  if abs(tnue.trueNuPDG) != 12 or abs(tnue.trueLepPDG != 11) or tnue.trueNuCCNC != 0 or isinf(tnue.xsecWeight):
    continue

  singleElectron = True
  for p in range(tnue.nTruePrimParts):
    pdg = abs(tnue.truePrimPartPDG[p])
    if pdg == 11:
      continue
    KE = tnue.truePrimPartE[p]
    if pdg not in [12,14,22]:
      mom2 = tnue.truePrimPartPx[p]**2 + tnue.truePrimPartPy[p]**2 + tnue.truePrimPartPz[p]**2
      try:
        mass = sqrt(tnue.truePrimPartE[p]**2 - mom2)
        KE = tnue.truePrimPartE[p] - mass
      except:
        print("pdg: ", pdg)
        print("E2: ", tnue.truePrimPartE[p]**2)
        print("mom2: ", mom2)
        sys.exit()
    if pdg in [111, 22, 1000180400, 311, 321, 3112, 3122, 3212, 3222, 130]:
      singleElectron = False
    elif pdg == 2212:
      if KE > 0.023:
        singleElectron = False
    elif pdg == 13:
      if KE > 0.009:
        singleElectron = False
    elif pdg == 211:
      if KE > 0.002:
        singleElectron = False
    else:
      if pdg in extra_particle_dict:
        extra_particle_dict[pdg] += 1
      else:
        extra_particle_dict[pdg] = 1

  if not singleElectron:
    continue

  n_raw_all += 1
  n_runs1to3_all += tnue.xsecWeight
  h_nCuts.Fill(tnue.trueLepE, tnue.xsecWeight)
  h2d_nCuts.Fill(tnue.trueLepE, tnue.trueNuE, tnue.xsecWeight)
  h2d_nCuts_raw.Fill(tnue.trueLepE, tnue.trueNuE)

  if tnue.foundVertex != 1:
    continue

  h_wCuts1.Fill(tnue.trueLepE, tnue.xsecWeight)
  n_runs1to3_passCuts1 += tnue.xsecWeight

  if tnue.vtxIsFiducial != 1:
    continue

  h_wCuts2.Fill(tnue.trueLepE, tnue.xsecWeight)
  n_runs1to3_passCuts2 += tnue.xsecWeight

  if tnue.vtxFracHitsOnCosmic >= 1.:
    continue

  h_wCuts3.Fill(tnue.trueLepE, tnue.xsecWeight)
  n_runs1to3_passCuts3 += tnue.xsecWeight

  foundMuon = False
  for iT in range(tnue.nTracks):
    if tnue.trackIsSecondary[iT] != 1 and tnue.trackClassified[iT] == 1 and tnue.trackPID[iT] == 13:
      foundMuon = True
      break

  if foundMuon:
    continue

  h_wCuts4.Fill(tnue.trueLepE, tnue.xsecWeight)
  n_runs1to3_passCuts4 += tnue.xsecWeight

  foundElectron = False
  elMaxQ = -1.
  elMaxQCosTheta = -1.
  elMaxQConf = -1.
  for iS in range(tnue.nShowers):
    if tnue.showerIsSecondary[iS] != 1 and tnue.showerClassified[iS] == 1 and tnue.showerPID[iS] == 11:
      foundElectron = True
      elConf = tnue.showerElScore[iS] - (tnue.showerPhScore[iS] + tnue.showerPiScore[iS])/2.
      if tnue.showerCharge[iS] > elMaxQ:
        elMaxQ = tnue.showerCharge[iS]
        elMaxQCosTheta = tnue.showerCosTheta[iS]
        elMaxQConf = elConf

  if not foundElectron:
    continue

  h_wCuts5.Fill(tnue.trueLepE, tnue.xsecWeight)
  n_runs1to3_passCuts5 += tnue.xsecWeight

  if elMaxQConf < 7.3:
    continue

  h_wCuts6.Fill(tnue.trueLepE, tnue.xsecWeight)
  n_runs1to3_passCuts6 += tnue.xsecWeight

  if elMaxQCosTheta < 0.:
    continue

  h_wCuts7.Fill(tnue.trueLepE, tnue.xsecWeight)
  h2d_wCuts7.Fill(tnue.trueLepE, tnue.trueNuE, tnue.xsecWeight)
  h2d_wCuts7_raw.Fill(tnue.trueLepE, tnue.trueNuE)
  n_runs1to3_passCuts7 += tnue.xsecWeight



n_runs1to3_all *= targetPOT/tnuePOTsum
n_runs1to3_passCuts1 *= targetPOT/tnuePOTsum
n_runs1to3_passCuts2 *= targetPOT/tnuePOTsum
n_runs1to3_passCuts3 *= targetPOT/tnuePOTsum
n_runs1to3_passCuts4 *= targetPOT/tnuePOTsum
n_runs1to3_passCuts5 *= targetPOT/tnuePOTsum
n_runs1to3_passCuts6 *= targetPOT/tnuePOTsum
n_runs1to3_passCuts7 *= targetPOT/tnuePOTsum

print()
print("n_raw_all: ", n_raw_all)
print()
print("n_runs1to3_all: ", n_runs1to3_all)
print()
print("n_runs1to3_passCuts1: ", n_runs1to3_passCuts1)
print("n_runs1to3_passCuts2: ", n_runs1to3_passCuts2)
print("n_runs1to3_passCuts3: ", n_runs1to3_passCuts3)
print("n_runs1to3_passCuts4: ", n_runs1to3_passCuts4)
print("n_runs1to3_passCuts5: ", n_runs1to3_passCuts5)
print("n_runs1to3_passCuts6: ", n_runs1to3_passCuts6)
print("n_runs1to3_passCuts7: ", n_runs1to3_passCuts7)
print()
print("cut set 1 efficiency: %f"%(n_runs1to3_passCuts1/n_runs1to3_all))
print("cut set 2 efficiency: %f"%(n_runs1to3_passCuts2/n_runs1to3_all))
print("cut set 3 efficiency: %f"%(n_runs1to3_passCuts3/n_runs1to3_all))
print("cut set 4 efficiency: %f"%(n_runs1to3_passCuts4/n_runs1to3_all))
print("cut set 5 efficiency: %f"%(n_runs1to3_passCuts5/n_runs1to3_all))
print("cut set 6 efficiency: %f"%(n_runs1to3_passCuts6/n_runs1to3_all))
print("cut set 7 efficiency: %f"%(n_runs1to3_passCuts7/n_runs1to3_all))
print()
print("encountered the following ignored particles (pdg number_found)")
for pdg in extra_particle_dict:
  print(pdg, extra_particle_dict[pdg])
print()

h_nCuts.Scale(targetPOT/tnuePOTsum)
h2d_nCuts.Scale(targetPOT/tnuePOTsum)
h_wCuts1.Scale(targetPOT/tnuePOTsum)
h_wCuts2.Scale(targetPOT/tnuePOTsum)
h_wCuts3.Scale(targetPOT/tnuePOTsum)
h_wCuts4.Scale(targetPOT/tnuePOTsum)
h_wCuts5.Scale(targetPOT/tnuePOTsum)
h_wCuts6.Scale(targetPOT/tnuePOTsum)
h_wCuts7.Scale(targetPOT/tnuePOTsum)
h2d_wCuts7.Scale(targetPOT/tnuePOTsum)

h_eff1.Divide(h_wCuts1,h_nCuts,1,1,"B")
h_eff2.Divide(h_wCuts2,h_nCuts,1,1,"B")
h_eff3.Divide(h_wCuts3,h_nCuts,1,1,"B")
h_eff4.Divide(h_wCuts4,h_nCuts,1,1,"B")
h_eff5.Divide(h_wCuts5,h_nCuts,1,1,"B")
h_eff6.Divide(h_wCuts6,h_nCuts,1,1,"B")
h_eff7.Divide(h_wCuts7,h_nCuts,1,1,"B")
h2d_eff7.Divide(h2d_wCuts7,h2d_nCuts)



outFile = rt.TFile(args.outfile, "RECREATE")

h_nCuts.Write()
h2d_nCuts.Write()
h2d_nCuts_raw.Write()
h_wCuts1.Write()
h_wCuts2.Write()
h_wCuts3.Write()
h_wCuts4.Write()
h_wCuts5.Write()
h_wCuts6.Write()
h_wCuts7.Write()
h2d_wCuts7.Write()
h2d_wCuts7_raw.Write()
h_eff1.Write()
h_eff2.Write()
h_eff3.Write()
h_eff4.Write()
h_eff5.Write()
h_eff6.Write()
h_eff7.Write()
h2d_eff7.Write()


cnv_E_1 = rt.TCanvas("cnv_E_1","cnv_E_1")
h_nCuts.Draw("EHIST")
h_wCuts1.Draw("EHISTSAME")
leg_E_1 = rt.TLegend(0.7,0.7,0.9,0.9)
leg_E_1.AddEntry(h_nCuts, "all true signal (%.2f)"%h_nCuts.Integral(), "l")
leg_E_1.AddEntry(h_wCuts1, "cut set 1 (%.2f)"%h_wCuts1.Integral(), "l")
leg_E_1.Draw()
cnv_E_1.Write()

cnv_eff_1 = rt.TCanvas("cnv_eff_1","cnv_eff_1")
h_eff1.Draw("E")
leg_eff_1 = rt.TLegend(0.7,0.7,0.9,0.9)
leg_eff_1.AddEntry(h_eff1, "cut set 1", "l")
leg_eff_1.Draw()
cnv_eff_1.Write()


cnv_E_2 = rt.TCanvas("cnv_E_2","cnv_E_2")
h_nCuts.Draw("EHIST")
h_wCuts1.Draw("EHISTSAME")
h_wCuts2.Draw("EHISTSAME")
leg_E_2 = rt.TLegend(0.7,0.7,0.9,0.9)
leg_E_2.AddEntry(h_nCuts, "all true signal (%.2f)"%h_nCuts.Integral(), "l")
leg_E_2.AddEntry(h_wCuts1, "cut set 1 (%.2f)"%h_wCuts1.Integral(), "l")
leg_E_2.AddEntry(h_wCuts2, "cut set 2 (%.2f)"%h_wCuts2.Integral(), "l")
leg_E_2.Draw()
cnv_E_2.Write()

cnv_eff_2 = rt.TCanvas("cnv_eff_2","cnv_eff_2")
h_eff1.Draw("E")
h_eff2.Draw("ESAME")
leg_eff_2 = rt.TLegend(0.7,0.7,0.9,0.9)
leg_eff_2.AddEntry(h_eff1, "cut set 1", "l")
leg_eff_2.AddEntry(h_eff2, "cut set 2", "l")
leg_eff_2.Draw()
cnv_eff_2.Write()


cnv_E_3 = rt.TCanvas("cnv_E_3","cnv_E_3")
h_nCuts.Draw("EHIST")
h_wCuts2.Draw("EHISTSAME")
h_wCuts3.Draw("EHISTSAME")
leg_E_3 = rt.TLegend(0.7,0.7,0.9,0.9)
leg_E_3.AddEntry(h_nCuts, "all true signal (%.2f)"%h_nCuts.Integral(), "l")
leg_E_3.AddEntry(h_wCuts2, "cut set 2 (%.2f)"%h_wCuts2.Integral(), "l")
leg_E_3.AddEntry(h_wCuts3, "cut set 3 (%.2f)"%h_wCuts3.Integral(), "l")
leg_E_3.Draw()
cnv_E_3.Write()

cnv_eff_3 = rt.TCanvas("cnv_eff_3","cnv_eff_3")
h_eff2.Draw("E")
h_eff3.Draw("ESAME")
leg_eff_3 = rt.TLegend(0.7,0.7,0.9,0.9)
leg_eff_3.AddEntry(h_eff2, "cut set 2", "l")
leg_eff_3.AddEntry(h_eff3, "cut set 3", "l")
leg_eff_3.Draw()
cnv_eff_3.Write()


cnv_E_4 = rt.TCanvas("cnv_E_4","cnv_E_4")
h_nCuts.Draw("EHIST")
h_wCuts3.Draw("EHISTSAME")
h_wCuts4.Draw("EHISTSAME")
leg_E_4 = rt.TLegend(0.7,0.7,0.9,0.9)
leg_E_4.AddEntry(h_nCuts, "all true signal (%.2f)"%h_nCuts.Integral(), "l")
leg_E_4.AddEntry(h_wCuts3, "cut set 3 (%.2f)"%h_wCuts3.Integral(), "l")
leg_E_4.AddEntry(h_wCuts4, "cut set 4 (%.2f)"%h_wCuts4.Integral(), "l")
leg_E_4.Draw()
cnv_E_4.Write()

cnv_eff_4 = rt.TCanvas("cnv_eff_4","cnv_eff_4")
h_eff3.Draw("E")
h_eff4.Draw("ESAME")
leg_eff_4 = rt.TLegend(0.7,0.7,0.9,0.9)
leg_eff_4.AddEntry(h_eff3, "cut set 3", "l")
leg_eff_4.AddEntry(h_eff4, "cut set 4", "l")
leg_eff_4.Draw()
cnv_eff_4.Write()


cnv_E_5 = rt.TCanvas("cnv_E_5","cnv_E_5")
h_nCuts.Draw("EHIST")
h_wCuts4.Draw("EHISTSAME")
h_wCuts5.Draw("EHISTSAME")
leg_E_5 = rt.TLegend(0.7,0.7,0.9,0.9)
leg_E_5.AddEntry(h_nCuts, "all true signal (%.2f)"%h_nCuts.Integral(), "l")
leg_E_5.AddEntry(h_wCuts4, "cut set 4 (%.2f)"%h_wCuts4.Integral(), "l")
leg_E_5.AddEntry(h_wCuts5, "cut set 5 (%.2f)"%h_wCuts5.Integral(), "l")
leg_E_5.Draw()
cnv_E_5.Write()

cnv_eff_5 = rt.TCanvas("cnv_eff_5","cnv_eff_5")
h_eff4.Draw("E")
h_eff5.Draw("ESAME")
leg_eff_5 = rt.TLegend(0.7,0.7,0.9,0.9)
leg_eff_5.AddEntry(h_eff4, "cut set 4", "l")
leg_eff_5.AddEntry(h_eff5, "cut set 5", "l")
leg_eff_5.Draw()
cnv_eff_5.Write()


cnv_E_6 = rt.TCanvas("cnv_E_6","cnv_E_6")
h_nCuts.Draw("EHIST")
h_wCuts5.Draw("EHISTSAME")
h_wCuts6.Draw("EHISTSAME")
leg_E_6 = rt.TLegend(0.7,0.7,0.9,0.9)
leg_E_6.AddEntry(h_nCuts, "all true signal (%.2f)"%h_nCuts.Integral(), "l")
leg_E_6.AddEntry(h_wCuts5, "cut set 5 (%.2f)"%h_wCuts5.Integral(), "l")
leg_E_6.AddEntry(h_wCuts6, "cut set 6 (%.2f)"%h_wCuts6.Integral(), "l")
leg_E_6.Draw()
cnv_E_6.Write()

cnv_eff_6 = rt.TCanvas("cnv_eff_6","cnv_eff_6")
h_eff5.Draw("E")
h_eff6.Draw("ESAME")
leg_eff_6 = rt.TLegend(0.7,0.7,0.9,0.9)
leg_eff_6.AddEntry(h_eff5, "cut set 5", "l")
leg_eff_6.AddEntry(h_eff6, "cut set 6", "l")
leg_eff_6.Draw()
cnv_eff_6.Write()


cnv_E_7 = rt.TCanvas("cnv_E_7","cnv_E_7")
h_nCuts.Draw("EHIST")
h_wCuts6.Draw("EHISTSAME")
h_wCuts7.Draw("EHISTSAME")
leg_E_7 = rt.TLegend(0.7,0.7,0.9,0.9)
leg_E_7.AddEntry(h_nCuts, "all true signal (%.2f)"%h_nCuts.Integral(), "l")
leg_E_7.AddEntry(h_wCuts6, "cut set 6 (%.2f)"%h_wCuts6.Integral(), "l")
leg_E_7.AddEntry(h_wCuts7, "cut set 7 (%.2f)"%h_wCuts7.Integral(), "l")
leg_E_7.Draw()
cnv_E_7.Write()

cnv_eff_7 = rt.TCanvas("cnv_eff_7","cnv_eff_7")
h_eff6.Draw("E")
h_eff7.Draw("ESAME")
leg_eff_7 = rt.TLegend(0.7,0.7,0.9,0.9)
leg_eff_7.AddEntry(h_eff6, "cut set 6", "l")
leg_eff_7.AddEntry(h_eff7, "cut set 7", "l")
leg_eff_7.Draw()
cnv_eff_7.Write()


cnv_E2d_0 = rt.TCanvas("cnv_E2d_0","cnv_E2d_0")
h2d_nCuts.Draw("COLZ")
cnv_E2d_0.Write()

cnv_E2d_7 = rt.TCanvas("cnv_E2d_7","cnv_E2d_7")
h2d_wCuts7.Draw("COLZ")
cnv_E2d_7.Write()

cnv_E2dRaw_0 = rt.TCanvas("cnv_E2dRaw_0","cnv_E2dRaw_0")
h2d_nCuts_raw.Draw("COLZ")
cnv_E2dRaw_0.Write()

cnv_E2dRaw_7 = rt.TCanvas("cnv_E2dRaw_7","cnv_E2dRaw_7")
h2d_wCuts7_raw.Draw("COLZ")
cnv_E2dRaw_7.Write()

cnv_eff2d_7 = rt.TCanvas("cnv_eff2d_7","cnv_eff2d_7")
h2d_eff7.Draw("COLZ")
cnv_eff2d_7.Write()



