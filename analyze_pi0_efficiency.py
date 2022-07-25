
import os,sys,argparse
import ROOT as rt
from larlite import larlite
from larlite import larutil
from ublarcvapp import ublarcvapp
from larcv import larcv
from larflow import larflow
from array import array
from event_weighting.event_weight_helper import SumPOT, Weights
from helpers.larflowreco_ana_funcs import *

parser = argparse.ArgumentParser("Test mB Excess Selection")
parser.add_argument("-f", "--files", required=True, type=str, nargs="+", help="input mdl files")
parser.add_argument("-o", "--outfile", type=str, default="analyze_pi0_efficiency_output.root", help="output file name")
parser.add_argument("-w", "--weightfile", required=True, type=str, help="weights file (pickled python dict)")
parser.add_argument("--debug", help="print debugging information", action="store_true")
args = parser.parse_args()


outRootFile = rt.TFile(args.outfile, "RECREATE")

potTree = rt.TTree("potTree","potTree")
totPOT = array('f', [0.])
totGoodPOT = array('f', [0.])
potTree.Branch("totPOT", totPOT, 'totPOT/F')
potTree.Branch("totGoodPOT", totGoodPOT, 'totGoodPOT/F')

eventTree = rt.TTree("EventTree","EventTree")
xsecWeight = array('f', [0.])
trueEnu = array('f', [0.])
trueCCNC = array('i', [0])
trueNuPDG = array('i', [0])
nPi0s_gtruth = array('i', [0])
nPi0s_mctruth = array('i', [0])
nGammaPairs = array('i', [0])
nGammasInFV = array('i', [0])
nPrimDecGammas = array('i', [0])
nPrimDecGammasInFV = array('i', [0])
eventTree.Branch("xsecWeight", xsecWeight, 'xsecWeight/F')
eventTree.Branch("trueEnu", trueEnu, 'trueEnu/F')
eventTree.Branch("trueCCNC", trueCCNC, 'trueCCNC/I')
eventTree.Branch("trueNuPDG", trueNuPDG, 'trueNuPDG/I')
eventTree.Branch("nPi0s_gtruth", nPi0s_gtruth, 'nPi0s_gtruth/I')
eventTree.Branch("nPi0s_mctruth", nPi0s_mctruth, 'nPi0s_mctruth/I')
eventTree.Branch("nGammaPairs", nGammaPairs, 'nGammaPairs/I')
eventTree.Branch("nGammasInFV", nGammasInFV, 'nGammasInFV/I')
eventTree.Branch("nPrimDecGammas", nPrimDecGammas, 'nPrimDecGammas/I')
eventTree.Branch("nPrimDecGammasInFV", nPrimDecGammasInFV, 'nPrimDecGammasInFV/I')

sce = larutil.SpaceChargeMicroBooNE()
mcNuVertexer = ublarcvapp.mctools.NeutrinoVertex()

totPOT_ = 0.
totGoodPOT_ = 0.
weights = Weights(args.weightfile)


#-------- begin file loop -----------------------------------------------------#
for mdlfile in args.files:

  ioll = larlite.storage_manager(larlite.storage_manager.kREAD)
  ioll.add_in_filename(mdlfile)
  ioll.open()

  #iolcv = larcv.IOManager(larcv.IOManager.kREAD, "larcv", larcv.IOManager.kTickBackward)
  #iolcv.add_in_file(mdlfile)
  #iolcv.reverse_all_products()
  #iolcv.initialize()

  potInFile, goodPotInFile = SumPOT(mdlfile)
  totPOT_ = totPOT_ + potInFile
  totGoodPOT_ = totGoodPOT_ + goodPotInFile

  #++++++ begin entry loop ++++++++++++++++++++++++++++++++++++++++++++++++++++=
  for ientry in range(ioll.get_entries()):
  
    ioll.go_to(ientry)
    #iolcv.read_entry(ientry)

    mctruth = ioll.get_data(larlite.data.kMCTruth, "generator")
    nuInt = mctruth.at(0).GetNeutrino()
    gtruth = ioll.get_data(larlite.data.kGTruth, "generator")
    npi0s = gtruth.at(0).fNumPi0
    mcNuPos = mctruth.at(0).GetNeutrino().Nu().Position()
    mcNuVertex = mcNuVertexer.getPos3DwSCE(ioll, sce)
    trueVtxPos = rt.TVector3(mcNuVertex[0], mcNuVertex[1], mcNuVertex[2])

    if npi0s < 1 or not isFiducialBig(trueVtxPos):
      continue

    try:
      xsecWeight[0] = weights.get(ioll.run_id(), ioll.subrun_id(), ioll.event_id())
    except:
      print("Couldn't find weight for run %i, subrun %i, event %i in %s!!!"%(ioll.run_id(), ioll.subrun_id(), ioll.event_id(), args.weightfile))
      continue

    trueEnu[0] = nuInt.Nu().Momentum().E()
    trueCCNC[0] = nuInt.CCNC()
    trueNuPDG[0] = nuInt.Nu().PdgCode()
    nPi0s_gtruth[0] = npi0s
    nPi0s_mctruth[0] = 0
    nGammaPairs[0] = 0
    nGammasInFV[0] = 0
    nPrimDecGammas[0] = 0
    nPrimDecGammasInFV[0] = 0

    for mcpart in mctruth.at(0).GetParticles():
      if abs(mcpart.PdgCode()) == 111 and mcpart.StatusCode() == 1:
        nPi0s_mctruth[0] = nPi0s_mctruth[0] + 1

    gammaDict = {}
    mcshowers = ioll.get_data(larlite.data.kMCShower, "mcreco")

    for mcshower in mcshowers:
      if mcshower.PdgCode() != 22 or mcshower.Process() != "Decay":
        continue
      mtid = mcshower.MotherTrackID()
      if mtid in gammaDict:
        gammaDict[mtid].append(mcshower)
      else:
        gammaDict[mtid] = [mcshower]

    for key, gammas in gammaDict.items():
      if len(gammas) == 2:
        nGammaPairs[0] = nGammaPairs[0] + 1
      for gamma in gammas:
        primaryDecay = getDistance(mcNuPos, gamma.Start()) < 1e-4
        if primaryDecay:
          nPrimDecGammas[0] = nPrimDecGammas[0] + 1
        if isFiducialBig(gamma.DetProfile()):
          nGammasInFV[0] = nGammasInFV[0] + 1
          if primaryDecay:
            nPrimDecGammasInFV[0] = nPrimDecGammasInFV[0] + 1

    if args.debug and nPi0s_gtruth[0] == 1 and nPrimDecGammas[0] > 2:
      print("Found >2 primary decay photons for single pi0 event: file %s, run %i, subrun %i, event %i, entry %i"%(mdlfile, ioll.run_id(), ioll.subrun_id(), ioll.event_id(), ientry))

    if args.debug and nPi0s_gtruth[0] > 1 and nPi0s_mctruth[0] == 0:
      print("Found weird pi0 count mismatch: file %s, run %i, subrun %i, event %i, entry %i"%(mdlfile, ioll.run_id(), ioll.subrun_id(), ioll.event_id(), ientry))

    eventTree.Fill()

  #++++++ end entry loop ++++++++++++++++++++++++++++++++++++++++++++++++++++=

  ioll.close()
  #iolcv.finalize()

#-------- end file loop -----------------------------------------------------#


totPOT[0] = totPOT_
totGoodPOT[0] = totGoodPOT_
potTree.Fill()

outRootFile.cd()
eventTree.Write("",rt.TObject.kOverwrite)
potTree.Write("",rt.TObject.kOverwrite)
outRootFile.Close()

