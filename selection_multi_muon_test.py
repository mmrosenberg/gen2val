
import sys, argparse
import numpy as np
import ROOT as rt

parser = argparse.ArgumentParser("test which muon to select for multi-reco-muon true CCnumu events")
parser.add_argument("-i", "--infile", type=str, default="flat_ntuples/dlgen2_reco_v2me06_ntuple_v5_mcc9_v28_wctagger_bnboverlay.root", help="bnb nu input file")
args = parser.parse_args()

fnu = rt.TFile(args.infile)
tnu = fnu.Get("EventTree")

nMultiMu = 0.
nChargeMatches = 0.
nScoreMatches = 0.

for i in range(tnu.GetEntries()):

  tnu.GetEntry(i)

  if abs(tnu.trueNuPDG) != 14 or tnu.trueNuCCNC != 0:
    continue

  if abs(tnu.trueLepPDG) != 13:
    sys.exit("ERROR: primary lepton pdg != 13 for true CCnumu event!!!")

  if tnu.foundVertex == 0 or tnu.vtxIsFiducial != 1:
    continue

  if tnu.vtxFracHitsOnCosmic >= 1.:
    continue

  recoMuons = []
  for iT in range(tnu.nTracks):
    if tnu.trackIsSecondary[iT] == 1 or tnu.trackClassified[iT] != 1:
      continue
    if tnu.trackPID[iT] == 13:
      recoMuons.append( [ tnu.trackCharge[iT], tnu.trackMuScore[iT], tnu.trackTrueTID[iT] ] )

  if len(recoMuons) < 2:
    continue

  trueMuTrackID = -1
  for iP in range(tnu.nTrueSimParts):
    if abs(tnu.trueSimPartPDG[iP]) == 13 and tnu.trueSimPartProcess[iP] == 0 and abs(tnu.trueLepE - tnu.trueSimPartE[iP]/1000.) < 1e-6:
      if trueMuTrackID != -1:
        sys.exit("ERROR: multiple true primary muons identified!!!")
      trueMuTrackID = tnu.trueSimPartTID[iP]
  if trueMuTrackID == -1:
    sys.exit("ERROR: couldn't find true primary muon in true CCnumu event!!!")

  nMultiMu += tnu.xsecWeight

  maxCharge = -999.
  maxScore = -999.
  maxChargeTID = -1
  maxScoreTID = -1
  for recoMuon in recoMuons:
    if recoMuon[0] > maxCharge:
      maxCharge = recoMuon[0]
      maxChargeTID = recoMuon[2]  
    if recoMuon[1] > maxScore:
      maxScore = recoMuon[1]
      maxScoreTID = recoMuon[2]  

  if maxChargeTID == trueMuTrackID:
    nChargeMatches += tnu.xsecWeight

  if maxScoreTID == trueMuTrackID:
    nScoreMatches += tnu.xsecWeight


chargeMatchFraction = nChargeMatches/nMultiMu
scoreMatchFraction = nScoreMatches/nMultiMu

print("fraction of multi-reco-muon CCnumu events where...")
print(f"max charge mu matches true mu: {chargeMatchFraction}")
print(f"max score mu matches true mu: {scoreMatchFraction}")

