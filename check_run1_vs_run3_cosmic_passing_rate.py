

import sys, argparse
import numpy as np
import ROOT as rt

from math import isinf
from helpers.plotting_functions import sortHists
from helpers.larflowreco_ana_funcs import isFiducial, getDistance, getCosThetaBeamVector, getCosThetaGravVector 


parser = argparse.ArgumentParser("Plot Selection Test Results")
parser.add_argument("-fr3", "--extbnb_run3_file", type=str, default="flat_ntuples/dlgen2_reco_v2me06_ntuple_v5_mcc9_v29e_dl_run3_G1_extbnb.root", help="extbnb input file")
parser.add_argument("-fr1", "--extbnb_run1_file", type=str, default="flat_ntuples/dlgen2_reco_v2me06_ntuple_v5_mcc9_v29e_dl_run1_C1_extbnb.root", help="extbnb input file")
parser.add_argument("-tr3", "--n_run3_triggers", type=float, default=39195178.0, help="number of EXT triggers in run3 file")
parser.add_argument("-tr1", "--n_run1_triggers", type=float, default=34202767.0, help="number of EXT triggers in run1 file")
parser.add_argument("-o", "--outfile", type=str, default="check_run1_vs_run3_cosmic_passing_rate_output.root", help="output root file name")
parser.add_argument("-vsc", "--vertexScoreCut", type=float, default=0.9, help="minimum neutrino vertex score")
parser.add_argument("-cgc", "--cosThetaGravCut", type=float, default=-0.9, help="minimum neutrino vertex score")
parser.add_argument("-nsc", "--fromNeutralScoreCut", type=float, default=-3.5, help="minimum neutrino vertex score")
args = parser.parse_args()

#rt.gROOT.SetBatch(True)
rt.TH1.SetDefaultSumw2(rt.kTRUE)
rt.gStyle.SetOptStat(0)

fextR3 = rt.TFile(args.extbnb_run3_file)
textR3 = fextR3.Get("EventTree")

fextR1 = rt.TFile(args.extbnb_run1_file)
textR1 = fextR1.Get("EventTree")

def configure_hists(hists, xtitle, linecolor):
  for hist in hists:
    hist.GetXaxis().SetTitle(xtitle)
    hist.GetYaxis().SetTitle("events per %.2e EXT triggers"%args.n_run3_triggers)
    hist.SetLineColor(linecolor)
    hist.SetLineWidth(2)
  return hists

h_recoE_wCutsSet1_run1 = rt.TH1F("h_recoE_ext_wCutsSet1_run1","ExtBNB Reco Nu Energy (cut subset 1)",11,0,2.2)
h_recoE_wCutsSet2_run1 = rt.TH1F("h_recoE_ext_wCutsSet2_run1","ExtBNB Reco Nu Energy (cut subset 2)",11,0,2.2)
h_recoE_wCutsSet3_run1 = rt.TH1F("h_recoE_ext_wCutsSet3_run1","ExtBNB Reco Nu Energy (cut subset 3)",11,0,2.2)
h_recoE_wCutsSet4_run1 = rt.TH1F("h_recoE_ext_wCutsSet4_run1","ExtBNB Reco Nu Energy (cut subset 4)",11,0,2.2)
h_recoE_wCutsSet5_run1 = rt.TH1F("h_recoE_ext_wCutsSet5_run1","ExtBNB Reco Nu Energy (cut subset 5)",11,0,2.2)
h_recoE_wCutsSet6_run1 = rt.TH1F("h_recoE_ext_wCutsSet6_run1","ExtBNB Reco Nu Energy (cut subset 6)",11,0,2.2)
h_recoE_wCutsSet7_run1 = rt.TH1F("h_recoE_ext_wCutsSet7_run1","ExtBNB Reco Nu Energy (all cuts)",11,0,2.2)
h_recoE_wCutsSet1_run1, h_recoE_wCutsSet2_run1, h_recoE_wCutsSet3_run1, h_recoE_wCutsSet4_run1, h_recoE_wCutsSet5_run1, h_recoE_wCutsSet6_run1, h_recoE_wCutsSet7_run1 = configure_hists((h_recoE_wCutsSet1_run1, h_recoE_wCutsSet2_run1, h_recoE_wCutsSet3_run1, h_recoE_wCutsSet4_run1, h_recoE_wCutsSet5_run1, h_recoE_wCutsSet6_run1, h_recoE_wCutsSet7_run1), "reco neutrino energy (GeV)", rt.kRed)

h_vtxX_wCutsSet1_run1 = rt.TH1F("h_vtxX_ext_wCutsSet1_run1","ExtBNB Reco Nu Vertex X (cut subset 1)",26,0,260)
h_vtxX_wCutsSet2_run1 = rt.TH1F("h_vtxX_ext_wCutsSet2_run1","ExtBNB Reco Nu Vertex X (cut subset 2)",26,0,260)
h_vtxX_wCutsSet3_run1 = rt.TH1F("h_vtxX_ext_wCutsSet3_run1","ExtBNB Reco Nu Vertex X (cut subset 3)",13,0,260)
h_vtxX_wCutsSet4_run1 = rt.TH1F("h_vtxX_ext_wCutsSet4_run1","ExtBNB Reco Nu Vertex X (cut subset 4)",13,0,260)
h_vtxX_wCutsSet5_run1 = rt.TH1F("h_vtxX_ext_wCutsSet5_run1","ExtBNB Reco Nu Vertex X (cut subset 5)",13,0,260)
h_vtxX_wCutsSet6_run1 = rt.TH1F("h_vtxX_ext_wCutsSet6_run1","ExtBNB Reco Nu Vertex X (cut subset 6)",13,0,260)
h_vtxX_wCutsSet7_run1 = rt.TH1F("h_vtxX_ext_wCutsSet7_run1","ExtBNB Reco Nu Vertex X (all cuts)",13,0,260)
h_vtxX_wCutsSet1_run1, h_vtxX_wCutsSet2_run1, h_vtxX_wCutsSet3_run1, h_vtxX_wCutsSet4_run1, h_vtxX_wCutsSet5_run1, h_vtxX_wCutsSet6_run1, h_vtxX_wCutsSet7_run1 = configure_hists((h_vtxX_wCutsSet1_run1, h_vtxX_wCutsSet2_run1, h_vtxX_wCutsSet3_run1, h_vtxX_wCutsSet4_run1, h_vtxX_wCutsSet5_run1, h_vtxX_wCutsSet6_run1, h_vtxX_wCutsSet7_run1), "reco nu vertex x (cm)", rt.kRed)

h_recoE_wCutsSet1_run3 = rt.TH1F("h_recoE_ext_wCutsSet1_run3","ExtBNB Reco Nu Energy (cut subset 1)",11,0,2.2)
h_recoE_wCutsSet2_run3 = rt.TH1F("h_recoE_ext_wCutsSet2_run3","ExtBNB Reco Nu Energy (cut subset 2)",11,0,2.2)
h_recoE_wCutsSet3_run3 = rt.TH1F("h_recoE_ext_wCutsSet3_run3","ExtBNB Reco Nu Energy (cut subset 3)",11,0,2.2)
h_recoE_wCutsSet4_run3 = rt.TH1F("h_recoE_ext_wCutsSet4_run3","ExtBNB Reco Nu Energy (cut subset 4)",11,0,2.2)
h_recoE_wCutsSet5_run3 = rt.TH1F("h_recoE_ext_wCutsSet5_run3","ExtBNB Reco Nu Energy (cut subset 5)",11,0,2.2)
h_recoE_wCutsSet6_run3 = rt.TH1F("h_recoE_ext_wCutsSet6_run3","ExtBNB Reco Nu Energy (cut subset 6)",11,0,2.2)
h_recoE_wCutsSet7_run3 = rt.TH1F("h_recoE_ext_wCutsSet7_run3","ExtBNB Reco Nu Energy (all cuts)",11,0,2.2)
h_recoE_wCutsSet1_run3, h_recoE_wCutsSet2_run3, h_recoE_wCutsSet3_run3, h_recoE_wCutsSet4_run3, h_recoE_wCutsSet5_run3, h_recoE_wCutsSet6_run3, h_recoE_wCutsSet7_run3 = configure_hists((h_recoE_wCutsSet1_run3, h_recoE_wCutsSet2_run3, h_recoE_wCutsSet3_run3, h_recoE_wCutsSet4_run3, h_recoE_wCutsSet5_run3, h_recoE_wCutsSet6_run3, h_recoE_wCutsSet7_run3), "reco neutrino energy (GeV)", rt.kBlue)

h_vtxX_wCutsSet1_run3 = rt.TH1F("h_vtxX_ext_wCutsSet1_run3","ExtBNB Reco Nu Vertex X (cut subset 1)",26,0,260)
h_vtxX_wCutsSet2_run3 = rt.TH1F("h_vtxX_ext_wCutsSet2_run3","ExtBNB Reco Nu Vertex X (cut subset 2)",26,0,260)
h_vtxX_wCutsSet3_run3 = rt.TH1F("h_vtxX_ext_wCutsSet3_run3","ExtBNB Reco Nu Vertex X (cut subset 3)",13,0,260)
h_vtxX_wCutsSet4_run3 = rt.TH1F("h_vtxX_ext_wCutsSet4_run3","ExtBNB Reco Nu Vertex X (cut subset 4)",13,0,260)
h_vtxX_wCutsSet5_run3 = rt.TH1F("h_vtxX_ext_wCutsSet5_run3","ExtBNB Reco Nu Vertex X (cut subset 5)",13,0,260)
h_vtxX_wCutsSet6_run3 = rt.TH1F("h_vtxX_ext_wCutsSet6_run3","ExtBNB Reco Nu Vertex X (cut subset 6)",13,0,260)
h_vtxX_wCutsSet7_run3 = rt.TH1F("h_vtxX_ext_wCutsSet7_run3","ExtBNB Reco Nu Vertex X (all cuts)",13,0,260)
h_vtxX_wCutsSet1_run3, h_vtxX_wCutsSet2_run3, h_vtxX_wCutsSet3_run3, h_vtxX_wCutsSet4_run3, h_vtxX_wCutsSet5_run3, h_vtxX_wCutsSet6_run3, h_vtxX_wCutsSet7_run3 = configure_hists((h_vtxX_wCutsSet1_run3, h_vtxX_wCutsSet2_run3, h_vtxX_wCutsSet3_run3, h_vtxX_wCutsSet4_run3, h_vtxX_wCutsSet5_run3, h_vtxX_wCutsSet6_run3, h_vtxX_wCutsSet7_run3), "reco nu vertex x (cm)", rt.kBlue)



print("starting loop over run1 extBNB events")

for i in range(textR1.GetEntries()):

  textR1.GetEntry(i)

  if textR1.foundVertex == 0 or textR1.vtxIsFiducial != 1:
    continue

  overflowRecoEVal = textR1.recoNuE/1000.
  if overflowRecoEVal > 2.0:
    overflowRecoEVal = 2.1

  h_recoE_wCutsSet1_run1.Fill(overflowRecoEVal)
  h_vtxX_wCutsSet1_run1.Fill(textR1.vtxX)

  if textR1.vtxFracHitsOnCosmic >= 1.:
    continue

  h_recoE_wCutsSet2_run1.Fill(overflowRecoEVal)
  h_vtxX_wCutsSet2_run1.Fill(textR1.vtxX)

  nMu = 0
  selMu_cosThetaBeam = -9.
  selMu_cosThetaGrav = -9.
  selMu_muScore = -99.
  selMu_proc = -9
  selMu_ntrlScore = -99.

  for iT in range(textR1.nTracks):
    if textR1.trackIsSecondary[iT] == 1 or textR1.trackClassified[iT] != 1:
      continue
    if textR1.trackPID[iT] == 13:
      nMu += 1
      if textR1.trackMuScore[iT] > selMu_muScore:
        selMu_cosThetaBeam = getCosThetaBeamVector(textR1.trackStartDirX[iT], textR1.trackStartDirY[iT], 
         textR1.trackStartDirZ[iT])
        selMu_cosThetaGrav = getCosThetaGravVector(textR1.trackStartDirX[iT], textR1.trackStartDirY[iT], 
         textR1.trackStartDirZ[iT])
        selMu_muScore = textR1.trackMuScore[iT]
        selMu_proc = textR1.trackProcess[iT]
        selMu_ntrlScore = textR1.trackFromNeutralScore[iT]

  if nMu < 1:
    continue

  h_recoE_wCutsSet3_run1.Fill(overflowRecoEVal)
  h_vtxX_wCutsSet3_run1.Fill(textR1.vtxX)

  if textR1.vtxScore < args.vertexScoreCut:
    continue

  h_recoE_wCutsSet4_run1.Fill(overflowRecoEVal)
  h_vtxX_wCutsSet4_run1.Fill(textR1.vtxX)

  if selMu_cosThetaGrav < args.cosThetaGravCut:
    continue

  h_recoE_wCutsSet5_run1.Fill(overflowRecoEVal)
  h_vtxX_wCutsSet5_run1.Fill(textR1.vtxX)

  if selMu_proc != 0:
    continue

  h_recoE_wCutsSet6_run1.Fill(overflowRecoEVal)
  h_vtxX_wCutsSet6_run1.Fill(textR1.vtxX)

  if selMu_ntrlScore > args.fromNeutralScoreCut:
    continue

  h_recoE_wCutsSet7_run1.Fill(overflowRecoEVal)
  h_vtxX_wCutsSet7_run1.Fill(textR1.vtxX)



print("starting loop over run3 extBNB events")

for i in range(textR3.GetEntries()):

  textR3.GetEntry(i)

  if textR3.foundVertex == 0 or textR3.vtxIsFiducial != 1:
    continue

  overflowRecoEVal = textR3.recoNuE/1000.
  if overflowRecoEVal > 2.0:
    overflowRecoEVal = 2.1

  h_recoE_wCutsSet1_run3.Fill(overflowRecoEVal)
  h_vtxX_wCutsSet1_run3.Fill(textR3.vtxX)

  if textR3.vtxFracHitsOnCosmic >= 1.:
    continue

  h_recoE_wCutsSet2_run3.Fill(overflowRecoEVal)
  h_vtxX_wCutsSet2_run3.Fill(textR3.vtxX)

  nMu = 0
  selMu_cosThetaBeam = -9.
  selMu_cosThetaGrav = -9.
  selMu_muScore = -99.
  selMu_proc = -9
  selMu_ntrlScore = -99.

  for iT in range(textR3.nTracks):
    if textR3.trackIsSecondary[iT] == 1 or textR3.trackClassified[iT] != 1:
      continue
    if textR3.trackPID[iT] == 13:
      nMu += 1
      if textR3.trackMuScore[iT] > selMu_muScore:
        selMu_cosThetaBeam = getCosThetaBeamVector(textR3.trackStartDirX[iT], textR3.trackStartDirY[iT], 
         textR3.trackStartDirZ[iT])
        selMu_cosThetaGrav = getCosThetaGravVector(textR3.trackStartDirX[iT], textR3.trackStartDirY[iT], 
         textR3.trackStartDirZ[iT])
        selMu_muScore = textR3.trackMuScore[iT]
        selMu_proc = textR3.trackProcess[iT]
        selMu_ntrlScore = textR3.trackFromNeutralScore[iT]

  if nMu < 1:
    continue

  h_recoE_wCutsSet3_run3.Fill(overflowRecoEVal)
  h_vtxX_wCutsSet3_run3.Fill(textR3.vtxX)

  if textR3.vtxScore < args.vertexScoreCut:
    continue

  h_recoE_wCutsSet4_run3.Fill(overflowRecoEVal)
  h_vtxX_wCutsSet4_run3.Fill(textR3.vtxX)

  if selMu_cosThetaGrav < args.cosThetaGravCut:
    continue

  h_recoE_wCutsSet5_run3.Fill(overflowRecoEVal)
  h_vtxX_wCutsSet5_run3.Fill(textR3.vtxX)

  if selMu_proc != 0:
    continue

  h_recoE_wCutsSet6_run3.Fill(overflowRecoEVal)
  h_vtxX_wCutsSet6_run3.Fill(textR3.vtxX)

  if selMu_ntrlScore > args.fromNeutralScoreCut:
    continue

  h_recoE_wCutsSet7_run3.Fill(overflowRecoEVal)
  h_vtxX_wCutsSet7_run3.Fill(textR3.vtxX)



h_recoE_wCutsSet1_run1.Scale(args.n_run3_triggers/args.n_run1_triggers)
h_recoE_wCutsSet2_run1.Scale(args.n_run3_triggers/args.n_run1_triggers)
h_recoE_wCutsSet3_run1.Scale(args.n_run3_triggers/args.n_run1_triggers)
h_recoE_wCutsSet4_run1.Scale(args.n_run3_triggers/args.n_run1_triggers)
h_recoE_wCutsSet5_run1.Scale(args.n_run3_triggers/args.n_run1_triggers)
h_recoE_wCutsSet6_run1.Scale(args.n_run3_triggers/args.n_run1_triggers)
h_recoE_wCutsSet7_run1.Scale(args.n_run3_triggers/args.n_run1_triggers)
h_vtxX_wCutsSet1_run1.Scale(args.n_run3_triggers/args.n_run1_triggers)
h_vtxX_wCutsSet2_run1.Scale(args.n_run3_triggers/args.n_run1_triggers)
h_vtxX_wCutsSet3_run1.Scale(args.n_run3_triggers/args.n_run1_triggers)
h_vtxX_wCutsSet4_run1.Scale(args.n_run3_triggers/args.n_run1_triggers)
h_vtxX_wCutsSet5_run1.Scale(args.n_run3_triggers/args.n_run1_triggers)
h_vtxX_wCutsSet6_run1.Scale(args.n_run3_triggers/args.n_run1_triggers)
h_vtxX_wCutsSet7_run1.Scale(args.n_run3_triggers/args.n_run1_triggers)



outFile = rt.TFile(args.outfile, "RECREATE")


cnv_recoE_cutSet1 = rt.TCanvas("cnv_recoE_cutSet1", "cnv_recoE_cutSet1")
hists_recoE_wCutsSet1 = sortHists([h_recoE_wCutsSet1_run1, h_recoE_wCutsSet1_run3])
hists_recoE_wCutsSet1[0].Draw("EHIST")
hists_recoE_wCutsSet1[1].Draw("EHISTSAME")
leg_recoE_cutSet1 = rt.TLegend(0.7,0.7,0.9,0.9)
leg_recoE_cutSet1.AddEntry(h_recoE_wCutsSet1_run1, "run1 (%.2f)"%h_recoE_wCutsSet1_run1.Integral(), "l")
leg_recoE_cutSet1.AddEntry(h_recoE_wCutsSet1_run3, "run3 (%.2f)"%h_recoE_wCutsSet1_run3.Integral(), "l")
leg_recoE_cutSet1.Draw()
cnv_recoE_cutSet1.Write()

cnv_vtxX_cutSet1 = rt.TCanvas("cnv_vtxX_cutSet1", "cnv_vtxX_cutSet1")
hists_vtxX_wCutsSet1 = sortHists([h_vtxX_wCutsSet1_run1, h_vtxX_wCutsSet1_run3])
hists_vtxX_wCutsSet1[0].Draw("EHIST")
hists_vtxX_wCutsSet1[1].Draw("EHISTSAME")
leg_vtxX_cutSet1 = rt.TLegend(0.7,0.7,0.9,0.9)
leg_vtxX_cutSet1.AddEntry(h_vtxX_wCutsSet1_run1, "run1 (%.2f)"%h_vtxX_wCutsSet1_run1.Integral(), "l")
leg_vtxX_cutSet1.AddEntry(h_vtxX_wCutsSet1_run3, "run3 (%.2f)"%h_vtxX_wCutsSet1_run3.Integral(), "l")
leg_vtxX_cutSet1.Draw()
cnv_vtxX_cutSet1.Write()


cnv_recoE_cutSet2 = rt.TCanvas("cnv_recoE_cutSet2", "cnv_recoE_cutSet2")
hists_recoE_wCutsSet2 = sortHists([h_recoE_wCutsSet2_run1, h_recoE_wCutsSet2_run3])
hists_recoE_wCutsSet2[0].Draw("EHIST")
hists_recoE_wCutsSet2[1].Draw("EHISTSAME")
leg_recoE_cutSet2 = rt.TLegend(0.7,0.7,0.9,0.9)
leg_recoE_cutSet2.AddEntry(h_recoE_wCutsSet2_run1, "run1 (%.2f)"%h_recoE_wCutsSet2_run1.Integral(), "l")
leg_recoE_cutSet2.AddEntry(h_recoE_wCutsSet2_run3, "run3 (%.2f)"%h_recoE_wCutsSet2_run3.Integral(), "l")
leg_recoE_cutSet2.Draw()
cnv_recoE_cutSet2.Write()

cnv_vtxX_cutSet2 = rt.TCanvas("cnv_vtxX_cutSet2", "cnv_vtxX_cutSet2")
hists_vtxX_wCutsSet2 = sortHists([h_vtxX_wCutsSet2_run1, h_vtxX_wCutsSet2_run3])
hists_vtxX_wCutsSet2[0].Draw("EHIST")
hists_vtxX_wCutsSet2[1].Draw("EHISTSAME")
leg_vtxX_cutSet2 = rt.TLegend(0.7,0.7,0.9,0.9)
leg_vtxX_cutSet2.AddEntry(h_vtxX_wCutsSet2_run1, "run1 (%.2f)"%h_vtxX_wCutsSet2_run1.Integral(), "l")
leg_vtxX_cutSet2.AddEntry(h_vtxX_wCutsSet2_run3, "run3 (%.2f)"%h_vtxX_wCutsSet2_run3.Integral(), "l")
leg_vtxX_cutSet2.Draw()
cnv_vtxX_cutSet2.Write()


cnv_recoE_cutSet3 = rt.TCanvas("cnv_recoE_cutSet3", "cnv_recoE_cutSet3")
hists_recoE_wCutsSet3 = sortHists([h_recoE_wCutsSet3_run1, h_recoE_wCutsSet3_run3])
hists_recoE_wCutsSet3[0].Draw("EHIST")
hists_recoE_wCutsSet3[1].Draw("EHISTSAME")
leg_recoE_cutSet3 = rt.TLegend(0.7,0.7,0.9,0.9)
leg_recoE_cutSet3.AddEntry(h_recoE_wCutsSet3_run1, "run1 (%.2f)"%h_recoE_wCutsSet3_run1.Integral(), "l")
leg_recoE_cutSet3.AddEntry(h_recoE_wCutsSet3_run3, "run3 (%.2f)"%h_recoE_wCutsSet3_run3.Integral(), "l")
leg_recoE_cutSet3.Draw()
cnv_recoE_cutSet3.Write()

cnv_vtxX_cutSet3 = rt.TCanvas("cnv_vtxX_cutSet3", "cnv_vtxX_cutSet3")
hists_vtxX_wCutsSet3 = sortHists([h_vtxX_wCutsSet3_run1, h_vtxX_wCutsSet3_run3])
hists_vtxX_wCutsSet3[0].Draw("EHIST")
hists_vtxX_wCutsSet3[1].Draw("EHISTSAME")
leg_vtxX_cutSet3 = rt.TLegend(0.7,0.7,0.9,0.9)
leg_vtxX_cutSet3.AddEntry(h_vtxX_wCutsSet3_run1, "run1 (%.2f)"%h_vtxX_wCutsSet3_run1.Integral(), "l")
leg_vtxX_cutSet3.AddEntry(h_vtxX_wCutsSet3_run3, "run3 (%.2f)"%h_vtxX_wCutsSet3_run3.Integral(), "l")
leg_vtxX_cutSet3.Draw()
cnv_vtxX_cutSet3.Write()


cnv_recoE_cutSet4 = rt.TCanvas("cnv_recoE_cutSet4", "cnv_recoE_cutSet4")
hists_recoE_wCutsSet4 = sortHists([h_recoE_wCutsSet4_run1, h_recoE_wCutsSet4_run3])
hists_recoE_wCutsSet4[0].Draw("EHIST")
hists_recoE_wCutsSet4[1].Draw("EHISTSAME")
leg_recoE_cutSet4 = rt.TLegend(0.7,0.7,0.9,0.9)
leg_recoE_cutSet4.AddEntry(h_recoE_wCutsSet4_run1, "run1 (%.2f)"%h_recoE_wCutsSet4_run1.Integral(), "l")
leg_recoE_cutSet4.AddEntry(h_recoE_wCutsSet4_run3, "run3 (%.2f)"%h_recoE_wCutsSet4_run3.Integral(), "l")
leg_recoE_cutSet4.Draw()
cnv_recoE_cutSet4.Write()

cnv_vtxX_cutSet4 = rt.TCanvas("cnv_vtxX_cutSet4", "cnv_vtxX_cutSet4")
hists_vtxX_wCutsSet4 = sortHists([h_vtxX_wCutsSet4_run1, h_vtxX_wCutsSet4_run3])
hists_vtxX_wCutsSet4[0].Draw("EHIST")
hists_vtxX_wCutsSet4[1].Draw("EHISTSAME")
leg_vtxX_cutSet4 = rt.TLegend(0.7,0.7,0.9,0.9)
leg_vtxX_cutSet4.AddEntry(h_vtxX_wCutsSet4_run1, "run1 (%.2f)"%h_vtxX_wCutsSet4_run1.Integral(), "l")
leg_vtxX_cutSet4.AddEntry(h_vtxX_wCutsSet4_run3, "run3 (%.2f)"%h_vtxX_wCutsSet4_run3.Integral(), "l")
leg_vtxX_cutSet4.Draw()
cnv_vtxX_cutSet4.Write()


cnv_recoE_cutSet5 = rt.TCanvas("cnv_recoE_cutSet5", "cnv_recoE_cutSet5")
hists_recoE_wCutsSet5 = sortHists([h_recoE_wCutsSet5_run1, h_recoE_wCutsSet5_run3])
hists_recoE_wCutsSet5[0].Draw("EHIST")
hists_recoE_wCutsSet5[1].Draw("EHISTSAME")
leg_recoE_cutSet5 = rt.TLegend(0.7,0.7,0.9,0.9)
leg_recoE_cutSet5.AddEntry(h_recoE_wCutsSet5_run1, "run1 (%.2f)"%h_recoE_wCutsSet5_run1.Integral(), "l")
leg_recoE_cutSet5.AddEntry(h_recoE_wCutsSet5_run3, "run3 (%.2f)"%h_recoE_wCutsSet5_run3.Integral(), "l")
leg_recoE_cutSet5.Draw()
cnv_recoE_cutSet5.Write()

cnv_vtxX_cutSet5 = rt.TCanvas("cnv_vtxX_cutSet5", "cnv_vtxX_cutSet5")
hists_vtxX_wCutsSet5 = sortHists([h_vtxX_wCutsSet5_run1, h_vtxX_wCutsSet5_run3])
hists_vtxX_wCutsSet5[0].Draw("EHIST")
hists_vtxX_wCutsSet5[1].Draw("EHISTSAME")
leg_vtxX_cutSet5 = rt.TLegend(0.7,0.7,0.9,0.9)
leg_vtxX_cutSet5.AddEntry(h_vtxX_wCutsSet5_run1, "run1 (%.2f)"%h_vtxX_wCutsSet5_run1.Integral(), "l")
leg_vtxX_cutSet5.AddEntry(h_vtxX_wCutsSet5_run3, "run3 (%.2f)"%h_vtxX_wCutsSet5_run3.Integral(), "l")
leg_vtxX_cutSet5.Draw()
cnv_vtxX_cutSet5.Write()


cnv_recoE_cutSet6 = rt.TCanvas("cnv_recoE_cutSet6", "cnv_recoE_cutSet6")
hists_recoE_wCutsSet6 = sortHists([h_recoE_wCutsSet6_run1, h_recoE_wCutsSet6_run3])
hists_recoE_wCutsSet6[0].Draw("EHIST")
hists_recoE_wCutsSet6[1].Draw("EHISTSAME")
leg_recoE_cutSet6 = rt.TLegend(0.7,0.7,0.9,0.9)
leg_recoE_cutSet6.AddEntry(h_recoE_wCutsSet6_run1, "run1 (%.2f)"%h_recoE_wCutsSet6_run1.Integral(), "l")
leg_recoE_cutSet6.AddEntry(h_recoE_wCutsSet6_run3, "run3 (%.2f)"%h_recoE_wCutsSet6_run3.Integral(), "l")
leg_recoE_cutSet6.Draw()
cnv_recoE_cutSet6.Write()

cnv_vtxX_cutSet6 = rt.TCanvas("cnv_vtxX_cutSet6", "cnv_vtxX_cutSet6")
hists_vtxX_wCutsSet6 = sortHists([h_vtxX_wCutsSet6_run1, h_vtxX_wCutsSet6_run3])
hists_vtxX_wCutsSet6[0].Draw("EHIST")
hists_vtxX_wCutsSet6[1].Draw("EHISTSAME")
leg_vtxX_cutSet6 = rt.TLegend(0.7,0.7,0.9,0.9)
leg_vtxX_cutSet6.AddEntry(h_vtxX_wCutsSet6_run1, "run1 (%.2f)"%h_vtxX_wCutsSet6_run1.Integral(), "l")
leg_vtxX_cutSet6.AddEntry(h_vtxX_wCutsSet6_run3, "run3 (%.2f)"%h_vtxX_wCutsSet6_run3.Integral(), "l")
leg_vtxX_cutSet6.Draw()
cnv_vtxX_cutSet6.Write()


cnv_recoE_cutSet7 = rt.TCanvas("cnv_recoE_cutSet7", "cnv_recoE_cutSet7")
hists_recoE_wCutsSet7 = sortHists([h_recoE_wCutsSet7_run1, h_recoE_wCutsSet7_run3])
hists_recoE_wCutsSet7[0].Draw("EHIST")
hists_recoE_wCutsSet7[1].Draw("EHISTSAME")
leg_recoE_cutSet7 = rt.TLegend(0.7,0.7,0.9,0.9)
leg_recoE_cutSet7.AddEntry(h_recoE_wCutsSet7_run1, "run1 (%.2f)"%h_recoE_wCutsSet7_run1.Integral(), "l")
leg_recoE_cutSet7.AddEntry(h_recoE_wCutsSet7_run3, "run3 (%.2f)"%h_recoE_wCutsSet7_run3.Integral(), "l")
leg_recoE_cutSet7.Draw()
cnv_recoE_cutSet7.Write()

cnv_vtxX_cutSet7 = rt.TCanvas("cnv_vtxX_cutSet7", "cnv_vtxX_cutSet7")
hists_vtxX_wCutsSet7 = sortHists([h_vtxX_wCutsSet7_run1, h_vtxX_wCutsSet7_run3])
hists_vtxX_wCutsSet7[0].Draw("EHIST")
hists_vtxX_wCutsSet7[1].Draw("EHISTSAME")
leg_vtxX_cutSet7 = rt.TLegend(0.7,0.7,0.9,0.9)
leg_vtxX_cutSet7.AddEntry(h_vtxX_wCutsSet7_run1, "run1 (%.2f)"%h_vtxX_wCutsSet7_run1.Integral(), "l")
leg_vtxX_cutSet7.AddEntry(h_vtxX_wCutsSet7_run3, "run3 (%.2f)"%h_vtxX_wCutsSet7_run3.Integral(), "l")
leg_vtxX_cutSet7.Draw()
cnv_vtxX_cutSet7.Write()


