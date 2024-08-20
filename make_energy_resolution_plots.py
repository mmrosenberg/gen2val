
import sys, argparse
import numpy as np
import ROOT as rt

parser = argparse.ArgumentParser("Plot Selection Test Results")
parser.add_argument("-fnu", "--bnbnu_file", type=str, default="flat_ntuples/dlgen2_reco_v2me06_ntuple_v5_mcc9_v28_wctagger_bnboverlay.root", help="bnb nu input file")
parser.add_argument("-fnue", "--bnbnue_file", type=str, default="flat_ntuples/dlgen2_reco_v2me06_ntuple_v5_mcc9_v28_wctagger_nueintrinsics.root", help="bnb nue input file")
parser.add_argument("-o", "--outfile", type=str, default="make_energy_resolution_plots_output.root", help="output root file name")
args = parser.parse_args()

rt.TH1.SetDefaultSumw2(rt.kTRUE)
rt.gStyle.SetOptStat(0)

fnu = rt.TFile(args.bnbnu_file)
tnu = fnu.Get("EventTree")
tnuPOT = fnu.Get("potTree")

fnue = rt.TFile(args.bnbnue_file)
tnue = fnue.Get("EventTree")
tnuePOT = fnue.Get("potTree")

targetPOT = 6.67e+20

tnuePOTsum = 0.
for i in range(tnuePOT.GetEntries()):
  tnuePOT.GetEntry(i)
  tnuePOTsum = tnuePOTsum + tnuePOT.totGoodPOT

tnuPOTsum = 0.
for i in range(tnuPOT.GetEntries()):
  tnuPOT.GetEntry(i)
  tnuPOTsum = tnuPOTsum + tnuPOT.totGoodPOT

h_nue = rt.TH2F("h_nue","Reconstructed Contained CCnue MC Events",50,0,5,50,0,5)
h_nue.GetXaxis().SetTitle("true neutrino energy (GeV)")
h_nue.GetYaxis().SetTitle("reconstructed neutrino energy (GeV)")
h_nue.GetZaxis().SetTitle("events per 6.67e+20 POT")
l_nue = rt.TLine(0,0,5,5)
l_nue.SetLineWidth(2)

h_numu = rt.TH2F("h_numu","Reconstructed Contained CCnumu MC Events",50,0,5,50,0,5)
h_numu.GetXaxis().SetTitle("true neutrino energy (GeV)")
h_numu.GetYaxis().SetTitle("reconstructed neutrino energy (GeV)")
h_numu.GetZaxis().SetTitle("events per 6.67e+20 POT")
l_numu = rt.TLine(0,0,5,5)
l_numu.SetLineWidth(2)

n_numu_preCuts = 0
n_numu_leptonReco = 0
n_nue_preCuts = 0
n_nue_leptonReco = 0

for i in range(tnu.GetEntries()):

  tnu.GetEntry(i)

  if abs(tnu.trueNuPDG) != 14 or tnu.trueNuCCNC != 0:
    continue

  if abs(tnu.trueLepPDG) != 13:
    sys.exit("ERROR: primary lepton pdg != 13 for true CCnumu event!!!")

  if tnu.foundVertex == 0 or tnu.vtxIsFiducial == 0 or tnu.vtxFracHitsOnCosmic >= 1 or tnu.vtxDistToTrue > 3.:
    continue

  trueMuTrackID = -1
  eventContained = True
  for iP in range(tnu.nTrueSimParts):
    if tnu.trueSimPartProcess[iP] == 0:
      if tnu.trueSimPartContained[iP] == 0:
        eventContained = False
        break
      if abs(tnu.trueSimPartPDG[iP]) == 13 and abs(tnu.trueLepE - tnu.trueSimPartE[iP]/1000.) < 1e-3:
        if trueMuTrackID != -1:
          sys.exit("ERROR: multiple true primary muons identified!!!")
        trueMuTrackID = tnu.trueSimPartTID[iP]

  if not eventContained:
    continue

  n_numu_preCuts += 1

  if trueMuTrackID == -1:
    sys.exit("ERROR: couldn't find true primary muon in true CCnumu event!!!")

  nuReconstructed = False
  for iT in range(tnu.nTracks):
    if tnu.trackTrueTID[iT] == trueMuTrackID and tnu.trackTrueComp[iT] > 0.6:
      nuReconstructed = True
      break

  if not nuReconstructed:
    continue

  n_numu_leptonReco += 1

  h_numu.Fill(tnu.trueNuE, tnu.recoNuE/1000., tnu.xsecWeight)



for i in range(tnue.GetEntries()):

  tnue.GetEntry(i)

  if abs(tnue.trueNuPDG) != 12 or tnue.trueNuCCNC != 0:
    continue

  if abs(tnue.trueLepPDG) != 11:
    sys.exit("ERROR: primary lepton pdg != 11 for true CCnue event!!!")

  if tnue.foundVertex == 0 or tnue.vtxIsFiducial == 0 or tnue.vtxFracHitsOnCosmic >= 1 or tnue.vtxDistToTrue > 3.:
    continue

  trueElTrackID = -1
  eventContained = True
  for iP in range(tnue.nTrueSimParts):
    if tnue.trueSimPartProcess[iP] == 0:
      if tnue.trueSimPartContained[iP] == 0:
        eventContained = False
        break
      if abs(tnue.trueSimPartPDG[iP]) == 11 and abs(tnue.trueLepE - tnue.trueSimPartE[iP]/1000.) < 1e-3:
        if trueElTrackID != -1:
          sys.exit("ERROR: multiple true primary electrons identified!!!")
        trueElTrackID = tnue.trueSimPartTID[iP]

  if not eventContained:
    continue

  if trueElTrackID == -1:
    electronInPrimary = False
    for iP in range(tnue.nTruePrimParts):
      if abs(tnue.truePrimPartPDG[iP]) == 11 and abs(tnue.trueLepE - tnue.truePrimPartE[iP]) < 1e-3:
        electronInPrimary = True
    if electronInPrimary:
      print("WARNING: primary electron not tracked in detsim, skipping event...")
      continue
    else:
      sys.exit("ERROR: couldn't find true primary electron in true CCnue event!!! entry, rsr = %i, %i, %i, %i"%(i,tnue.run,tnue.subrun,tnue.event))

  n_nue_preCuts += 1

  nuReconstructed = False
  for iS in range(tnue.nShowers):
    if tnue.showerTrueTID[iS] == trueElTrackID and tnue.showerTrueComp[iS] > 0.6:
      nuReconstructed = True
      break

  if not nuReconstructed:
    continue

  n_nue_leptonReco += 1

  h_nue.Fill(tnue.trueNuE, tnue.recoNuE/1000., tnue.xsecWeight)


n_numu_preCuts_scaled = n_numu_preCuts*(targetPOT/tnuPOTsum)
n_numu_leptonReco_scaled = n_numu_leptonReco*(targetPOT/tnuPOTsum)
n_nue_preCuts_scaled = n_nue_preCuts*(targetPOT/tnuePOTsum)
n_nue_leptonReco_scaled = n_nue_leptonReco*(targetPOT/tnuePOTsum)

print("of raw events with a reco vertex and all true final state particles contained...")
print("%i/%i CCnumu events had reconstructed muon"%(n_numu_leptonReco,n_numu_preCuts))
print("%i/%i CCnue events had reconstructed electron"%(n_nue_leptonReco,n_nue_preCuts))

print("of predicted runs1-3 event counts with a reco vertex and all true final state particles contained...")
print("%i/%i CCnumu events had reconstructed muon"%(n_numu_leptonReco_scaled,n_numu_preCuts_scaled))
print("%i/%i CCnue events had reconstructed electron"%(n_nue_leptonReco_scaled,n_nue_preCuts_scaled))

h_numu.Scale(targetPOT/tnuPOTsum)
h_nue.Scale(targetPOT/tnuePOTsum)


outFile = rt.TFile(args.outfile, "RECREATE")

cnv_numu = rt.TCanvas("cnv_numu","cnv_numu")
h_numu.Draw("COLZ")
l_numu.Draw()
cnv_numu.Write()

cnv_nue = rt.TCanvas("cnv_nue","cnv_nue")
h_nue.Draw("COLZ")
l_nue.Draw()
cnv_nue.Write()

h_numu.Write()
h_nue.Write()
l_numu.Write()
l_nue.Write()


