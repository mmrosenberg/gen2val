
import sys, argparse
import numpy as np
import ROOT as rt

from math import isinf
from helpers.plotting_functions import sortHists
from helpers.larflowreco_ana_funcs import isFiducial, isFiducialWC, getDistance


parser = argparse.ArgumentParser("Select CCnue Inclusive Events from Flat Ntuples")
parser.add_argument("-fnu", "--bnbnu_file", type=str, default="flat_ntuples/dlgen2_reco_v2me06_ntuple_v1_mcc9_v28_wctagger_bnboverlay_partial.root", help="bnb nu input file")
parser.add_argument("-fnue", "--bnbnue_file", type=str, default="flat_ntuples/dlgen2_reco_v2me06_ntuple_v1_mcc9_v28_wctagger_nueintrinsics_partial.root", help="bnb nu input file")
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
parser.add_argument("--oldVertexVar", help="use old nVertices variable for vertex found cut", action="store_true")
parser.add_argument("--smallFV", help="use 20cm fiducial volume", action="store_true")
parser.add_argument("--printEXTinfo", help="print event info for selected cosmic background", action="store_true")
args = parser.parse_args()


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

fnu_trimmed = rt.TFile(args.bnbnu_file.replace(".root","_trimmed_CCnueInc_passed.root"),"RECREATE")
tnu_trimmed = tnu.CloneTree(0)

fnue_trimmed = rt.TFile(args.bnbnue_file.replace(".root","_trimmed_CCnueInc_passed.root"),"RECREATE")
tnue_trimmed = tnue.CloneTree(0)

fext_trimmed = rt.TFile(args.extbnb_file.replace(".root","_trimmed_CCnueInc_passed.root"),"RECREATE")
text_trimmed = text.CloneTree(0)

fdata_trimmed = rt.TFile(args.data_file.replace(".root","_trimmed_CCnueInc_passed.root"),"RECREATE")
tdata_trimmed = tdata.CloneTree(0)



for i in range(tnu.GetEntries()):

  tnu.GetEntry(i)

  if isinf(tnu.xsecWeight):
    continue

  if args.smallFV:
    trueVtxPos = rt.TVector3(tnu.trueVtxX, tnu.trueVtxY, tnu.trueVtxZ)
    if not isFiducial(trueVtxPos):
      continue

  if abs(tnu.trueNuPDG) != 14 and not (abs(tnu.trueNuPDG) == 12 and tnu.trueNuCCNC == 1):
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

  if tnu.vtxFracHitsOnCosmic >= args.vertexFracOnCosCut:
    continue

  nMuons = 0
  nElectrons = 0
  elMaxQCosTheta = -1.
  elMaxQConf = -1.
  elMaxQFrac = -1.
  elMaxQ = -1.

  for iT in range(tnu.nTracks):
    if tnu.trackIsSecondary[iT] != 1 and tnu.trackClassified[iT] == 1 and tnu.trackPID[iT] == 13:
      nMuons += 1

  if nMuons != 0:
    continue

  for iS in range(tnu.nShowers):
    if tnu.showerIsSecondary[iS] == 1 or tnu.showerClassified[iS] == 0 or tnu.showerDistToVtx[iS] > args.distCut or tnu.showerComp[iS] < args.compCut or tnu.showerPurity[iS] < args.purityCut:
      continue
    if tnu.showerPID[iS] == 11:
      nElectrons += 1
      elConf = tnu.showerElScore[iS] - (tnu.showerPhScore[iS] + tnu.showerPiScore[iS])/2.
      if tnu.showerCharge[iS] > elMaxQ:
        elMaxQ = tnu.showerCharge[iS]
        elMaxQCosTheta = tnu.showerCosTheta[iS]
        elMaxQConf = elConf
        elMaxQFrac = tnu.showerChargeFrac[iS]

  if nElectrons >= 1:
    if elMaxQConf > args.confCut and elMaxQ > args.chargeCut and elMaxQFrac > args.chargeFracCut and elMaxQCosTheta > args.cosThetaCut:
      tnu_trimmed.Fill()



for i in range(tnue.GetEntries()):

  tnue.GetEntry(i)

  if abs(tnue.trueNuPDG) != 12 or tnue.trueNuCCNC != 0 or isinf(tnue.xsecWeight):
    continue

  if args.smallFV:
    trueVtxPos = rt.TVector3(tnue.trueVtxX, tnue.trueVtxY, tnue.trueVtxZ)
    if not isFiducial(trueVtxPos):
      continue

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

  if tnue.vtxFracHitsOnCosmic >= args.vertexFracOnCosCut:
    continue

  nMuons = 0
  nElectrons = 0
  elMaxQCosTheta = -1.
  elMaxQConf = -1.
  elMaxQFrac = -1.
  elMaxQ = -1.

  for iT in range(tnue.nTracks):
    if tnue.trackIsSecondary[iT] != 1 and tnue.trackClassified[iT] == 1 and tnue.trackPID[iT] == 13:
      nMuons += 1

  if nMuons != 0:
    continue

  for iS in range(tnue.nShowers):
    if tnue.showerIsSecondary[iS] == 1 or tnue.showerClassified[iS] == 0 or tnue.showerDistToVtx[iS] > args.distCut or tnue.showerComp[iS] < args.compCut or tnue.showerPurity[iS] < args.purityCut:
      continue
    if tnue.showerPID[iS] == 11:
      nElectrons += 1
      elConf = tnue.showerElScore[iS] - (tnue.showerPhScore[iS] + tnue.showerPiScore[iS])/2.
      if tnue.showerCharge[iS] > elMaxQ:
        elMaxQ = tnue.showerCharge[iS]
        elMaxQCosTheta = tnue.showerCosTheta[iS]
        elMaxQConf = elConf
        elMaxQFrac = tnue.showerChargeFrac[iS]

  if nElectrons >= 1:
    if elMaxQConf > args.confCut and elMaxQ > args.chargeCut and elMaxQFrac > args.chargeFracCut and elMaxQCosTheta > args.cosThetaCut:
      tnue_trimmed.Fill()



for i in range(text.GetEntries()):

  text.GetEntry(i)

  vtxPos = rt.TVector3(text.vtxX, text.vtxY, text.vtxZ)
  if args.smallFV:
    vtxIsFiducial = isFiducial(vtxPos)
  else:
    vtxIsFiducial = isFiducialWC(vtxPos)
  if args.oldVertexVar:
    if text.nVertices < 1 or not vtxIsFiducial: #text.vtxIsFiducial != 1:
      continue
  else:
    if text.foundVertex == 0 or not vtxIsFiducial: #text.vtxIsFiducial != 1:
      continue

  if text.vtxFracHitsOnCosmic >= args.vertexFracOnCosCut:
    continue

  nMuons = 0
  nElectrons = 0
  elMaxQCosTheta = -1.
  elMaxQConf = -1.
  elMaxQFrac = -1.
  elMaxQ = -1.

  for iT in range(text.nTracks):
    if text.trackIsSecondary[iT] != 1 and text.trackClassified[iT] == 1 and text.trackPID[iT] == 13:
      nMuons += 1

  if nMuons != 0:
    continue

  for iS in range(text.nShowers):
    if text.showerIsSecondary[iS] == 1 or text.showerClassified[iS] == 0 or text.showerDistToVtx[iS] > args.distCut or text.showerComp[iS] < args.compCut or text.showerPurity[iS] < args.purityCut:
      continue
    if text.showerPID[iS] == 11:
      nElectrons += 1
      elConf = text.showerElScore[iS] - (text.showerPhScore[iS] + text.showerPiScore[iS])/2.
      if text.showerCharge[iS] > elMaxQ:
        elMaxQ = text.showerCharge[iS]
        elMaxQCosTheta = text.showerCosTheta[iS]
        elMaxQConf = elConf
        elMaxQFrac = text.showerChargeFrac[iS]

  if nElectrons >= 1:
    if elMaxQConf > args.confCut and elMaxQ > args.chargeCut and elMaxQFrac > args.chargeFracCut and elMaxQCosTheta > args.cosThetaCut:
      text_trimmed.Fill()
      if args.printEXTinfo:
        print("extBNB background event passed selection: (fileid, run, subrun, event, vtxX, vtxY, vtxZ) = (%i, %i, %i, %i, %f, %f, %f)"%(text.fileid, text.run, text.subrun, text.event, text.vtxX, text.vtxY, text.vtxZ))



for i in range(tdata.GetEntries()):

  tdata.GetEntry(i)

  vtxPos = rt.TVector3(tdata.vtxX, tdata.vtxY, tdata.vtxZ)
  if args.smallFV:
    vtxIsFiducial = isFiducial(vtxPos)
  else:
    vtxIsFiducial = isFiducialWC(vtxPos)
  if tdata.nVertices < 1 or not vtxIsFiducial: #tdata.vtxIsFiducial != 1:
    continue

  if tdata.vtxFracHitsOnCosmic >= args.vertexFracOnCosCut:
    continue

  nMuons = 0
  nElectrons = 0
  elMaxQCosTheta = -1.
  elMaxQConf = -1.
  elMaxQFrac = -1.
  elMaxQ = -1.

  for iT in range(tdata.nTracks):
    if tdata.trackIsSecondary[iT] != 1 and tdata.trackClassified[iT] == 1 and tdata.trackPID[iT] == 13:
      nMuons += 1

  if nMuons != 0:
    continue

  for iS in range(tdata.nShowers):
    if tdata.showerIsSecondary[iS] == 1 or tdata.showerClassified[iS] == 0 or tdata.showerDistToVtx[iS] > args.distCut or tdata.showerComp[iS] < args.compCut or tdata.showerPurity[iS] < args.purityCut:
      continue
    if tdata.showerPID[iS] == 11:
      nElectrons += 1
      elConf = tdata.showerElScore[iS] - (tdata.showerPhScore[iS] + tdata.showerPiScore[iS])/2.
      if tdata.showerCharge[iS] > elMaxQ:
        elMaxQ = tdata.showerCharge[iS]
        elMaxQCosTheta = tdata.showerCosTheta[iS]
        elMaxQConf = elConf
        elMaxQFrac = tdata.showerChargeFrac[iS]

  if nElectrons >= 1:
    if elMaxQConf > args.confCut and elMaxQ > args.chargeCut and elMaxQFrac > args.chargeFracCut and elMaxQCosTheta > args.cosThetaCut:
      tdata.Fill()


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



