
import sys, argparse
import numpy as np
import ROOT as rt

from math import isinf
from helpers.larflowreco_ana_funcs import isFiducialWC


parser = argparse.ArgumentParser("Make Vertex Resolution Plot")
parser.add_argument("-fnu", "--bnbnu_file", type=str, default="selection_output/prepare_selection_test_output/prepare_selection_test_reco_v2me05_gen2val_v23_nu_file.root", help="bnb nu input file")
parser.add_argument("-fnue", "--bnbnue_file", type=str, default="selection_output/prepare_selection_test_output/prepare_selection_test_reco_v2me05_gen2val_v23_nue_file.root", help="bnb nu input file")
args = parser.parse_args()

rt.TH1.SetDefaultSumw2(rt.kTRUE)

fnu = rt.TFile(args.bnbnu_file)
tnu = fnu.Get("EventTree")
tnuPOT = fnu.Get("potTree")

fnue = rt.TFile(args.bnbnue_file)
tnue = fnue.Get("EventTree")
tnuePOT = fnue.Get("potTree")

runs1to3POT = 6.67e+20

tnuePOTsum = 0.
for i in range(tnuePOT.GetEntries()):
  tnuePOT.GetEntry(i)
  tnuePOTsum = tnuePOTsum + tnuePOT.totGoodPOT

tnuPOTsum = 0.
for i in range(tnuPOT.GetEntries()):
  tnuPOT.GetEntry(i)
  tnuPOTsum = tnuPOTsum + tnuPOT.totGoodPOT

h_nu = rt.TH1F("h_nu","",300,0,7)
h_nue = rt.TH1F("h_nue","",300,0,7)

for i in range(tnu.GetEntries()):
  tnu.GetEntry(i)
  vtxPos = rt.TVector3(tnu.vtxX, tnu.vtxY, tnu.vtxZ)
  if isinf(tnu.xsecWeight) or tnu.nVertices < 1 or not isFiducialWC(vtxPos) or tnu.vtxFracHitsOnCosmic >= 1.:
    continue
  h_nu.Fill(tnu.vtxDistToTrue, tnu.xsecWeight)

for i in range(tnue.GetEntries()):
  tnue.GetEntry(i)
  vtxPos = rt.TVector3(tnue.vtxX, tnue.vtxY, tnue.vtxZ)
  if isinf(tnue.xsecWeight) or tnue.nVertices < 1 or not isFiducialWC(vtxPos) or tnue.vtxFracHitsOnCosmic >= 1.:
    continue
  h_nue.Fill(tnue.vtxDistToTrue, tnue.xsecWeight)

h_nu.Scale(runs1to3POT/tnuPOTsum)
h_nue.Scale(runs1to3POT/tnuePOTsum)

h_vtxRes = rt.TH1F("h_vtxRes","Vertex Resolution for MC Neutrino Interactions",300,0,7)
h_vtxRes.GetYaxis().SetTitle("area normalized event count")
h_vtxRes.GetXaxis().SetTitle("distance from reconstructed to simulated nuetrino vertex (cm)")
h_vtxRes.SetLineWidth(2)

h_vtxRes.Add(h_nu, h_nue)
h_vtxRes.Scale(1.0/h_vtxRes.Integral())

bin_sum = 0.
found68 = False
reached1cm = False
reached2cm = False
reached3cm = False
for i in range(1, h_vtxRes.GetNbinsX()+1):
  bin_sum += h_vtxRes.GetBinContent(i)
  bin_center = h_vtxRes.GetBinCenter(i)
  if bin_sum >= 0.68 and not found68:
    found68 = True
    print("68 percent containment resolution: %.2f cm"%bin_center)
  if bin_center >= 1. and not reached1cm:
    reached1cm = True
    print("fraction within 1cm: %f"%bin_sum)
  if bin_center >= 2. and not reached2cm:
    reached2cm = True
    print("fraction within 2cm: %f"%bin_sum)
  if bin_center >= 3. and not reached3cm:
    reached3cm = True
    print("fraction within 3cm: %f"%bin_sum)

h_vtxRes.Draw("EHIST")

input("Press Enter to continue...")

