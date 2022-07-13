
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
parser.add_argument("-f", "--files", required=True, type=str, nargs="+", help="input kpsreco files")
parser.add_argument("-t", "--truth", required=True, type=str, help="text file containing merged_dlreco list")
parser.add_argument("-o", "--outfile", type=str, default="test_mBexcess_selection_output.root", help="output file name")
parser.add_argument("-w", "--weightfile", type=str, default="none", help="weights file (pickled python dict)")
parser.add_argument("-mc", "--isMC", help="running over MC input", action="store_true")
parser.add_argument("--printTrueEvents", help="print out event info for true signal events", action="store_true")
parser.add_argument("--printRecoEvents", help="print out event info for reco signal events", action="store_true")
parser.add_argument("--oldVtxBranch", help="use nufitted_v instead of nuvetoed_v for old reco", action="store_true")
args = parser.parse_args()

if args.isMC and args.weightfile=="none":
  sys.exit("Must supply weight file for MC input. Exiting...")

reco2Tag = "merged_dlana_"
if args.isMC:
  reco2Tag = "merged_dlreco_"
files = getFiles(reco2Tag, args.files, args.truth)

outRootFile = rt.TFile(args.outfile, "RECREATE")

if args.isMC:
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
passSelTruth = array('i', [0])
passSelReco = array('i', [0])
trueNMus = array('i', [0])
trueNThreshMus = array('i', [0])
trueNCPis = array('i', [0])
trueNThreshCPis = array('i', [0])
trueNPrs = array('i', [0])
trueNThreshPrs = array('i', [0])
trueNPi0s = array('i', [0])
trueNPhs = array('i', [0])
trueNEls = array('i', [0])
recoSingleVtx = array('i', [0])
recoNNegLLTrks = array('i', [0])
recoNNegLLThreshTrks = array('i', [0])
recoNPosLLTrks = array('i', [0])
recoNPosLLThreshTrks = array('i', [0])
recoVtxDistToTrue = array('f', [0.])
recoShowerNHits = array('i', [0])
eventTree.Branch("xsecWeight", xsecWeight, 'xsecWeight/F')
eventTree.Branch("trueEnu", trueEnu, 'trueEnu/F')
eventTree.Branch("trueCCNC", trueCCNC, 'trueCCNC/I')
eventTree.Branch("trueNuPDG", trueNuPDG, 'trueNuPDG/I')
eventTree.Branch("passSelTruth", passSelTruth, 'passSelTruth/I')
eventTree.Branch("passSelReco", passSelReco, 'passSelReco/I')
eventTree.Branch("trueNMus", trueNMus, 'trueNMus/I')
eventTree.Branch("trueNThreshMus", trueNThreshMus, 'trueNThreshMus/I')
eventTree.Branch("trueNCPis", trueNCPis, 'trueNCPis/I')
eventTree.Branch("trueNThreshCPis", trueNThreshCPis, 'trueNThreshCPis/I')
eventTree.Branch("trueNPrs", trueNPrs, 'trueNPrs/I')
eventTree.Branch("trueNThreshPrs", trueNThreshPrs, 'trueNThreshPrs/I')
eventTree.Branch("trueNPi0s", trueNPi0s, 'trueNPi0s/I')
eventTree.Branch("trueNPhs", trueNPhs, 'trueNPhs/I')
eventTree.Branch("trueNEls", trueNEls, 'trueNEls/I')
eventTree.Branch("recoSingleVtx", recoSingleVtx, 'recoSingleVtx/I')
eventTree.Branch("recoNNegLLTrks", recoNNegLLTrks, 'recoNNegLLTrks/I')
eventTree.Branch("recoNNegLLThreshTrks", recoNNegLLThreshTrks, 'recoNNegLLThreshTrks/I')
eventTree.Branch("recoNPosLLTrks", recoNPosLLTrks, 'recoNPosLLTrks/I')
eventTree.Branch("recoNPosLLThreshTrks", recoNPosLLThreshTrks, 'recoNPosLLThreshTrks/I')
eventTree.Branch("recoVtxDistToTrue", recoVtxDistToTrue, 'recoVtxDistToTrue/F')
eventTree.Branch("recoShowerNHits", recoShowerNHits, 'recoShowerNHits/I')

sce = larutil.SpaceChargeMicroBooNE()
mcNuVertexer = ublarcvapp.mctools.NeutrinoVertex()

if args.isMC:
  totPOT_ = 0.
  totGoodPOT_ = 0.
  weights = Weights(args.weightfile)


#-------- begin file loop -----------------------------------------------------#
for filepair in files:

  ioll = larlite.storage_manager(larlite.storage_manager.kREAD)
  ioll.add_in_filename(filepair[1])
  ioll.open()

  iolcv = larcv.IOManager(larcv.IOManager.kREAD, "larcv", larcv.IOManager.kTickBackward)
  iolcv.add_in_file(filepair[1])
  iolcv.reverse_all_products()
  iolcv.initialize()

  kpsfile = rt.TFile(filepair[0])
  kpst = kpsfile.Get("KPSRecoManagerTree")

  try:
    nKPSTEntries = kpst.GetEntries()
  except:
    print("%s is empty. skipping..."%(filepair[0]))
    ioll.close()
    iolcv.finalize()
    kpsfile.Close()
    continue

  if args.isMC:
    potInFile, goodPotInFile = SumPOT(filepair[1])
    totPOT_ = totPOT_ + potInFile
    totGoodPOT_ = totGoodPOT_ + goodPotInFile

  #++++++ begin entry loop ++++++++++++++++++++++++++++++++++++++++++++++++++++=
  for ientry in range(ioll.get_entries()):
  
    ioll.go_to(ientry)
    iolcv.read_entry(ientry)
    kpst.GetEntry(ientry)

    vertices = kpst.nuvetoed_v
    if args.oldVtxBranch:
      vertices = kpst.nufitted_v
  
    if kpst.run != ioll.run_id() or kpst.subrun != ioll.subrun_id() or kpst.event != ioll.event_id():
      print("EVENTS DON'T MATCH!!!")
      print("truth run/subrun/event: %i/%i/%i"%(ioll.run_id(),ioll.subrun_id(),ioll.event_id()))
      print("reco run/subrun/event: %i/%i/%i"%(kpst.run,kpst.subrun,kpst.event))
      continue

    xsecWeight[0] = -99.
    trueEnu[0] = -99.
    trueCCNC[0] = -1
    trueNuPDG[0] = -1
    passSelTruth[0] = -1
    trueNMus[0] = -1
    trueNThreshMus[0] = -1
    trueNCPis[0] = -1
    trueNThreshCPis[0] = -1
    trueNPrs[0] = -1
    trueNThreshPrs[0] = -1
    trueNPi0s[0] = -1
    trueNPhs[0] = -1
    trueNEls[0] = -1

    if args.isMC:
      mctruth = ioll.get_data(larlite.data.kMCTruth, "generator")
      nuInt = mctruth.at(0).GetNeutrino()
      #lep = nuInt.Lepton()
      mcNuVertex = mcNuVertexer.getPos3DwSCE(ioll, sce)
      trueVtxPos = rt.TVector3(mcNuVertex[0], mcNuVertex[1], mcNuVertex[2])

      #if nuInt.CCNC() != 0 or lep.PdgCode() != 11 or not isFiducial(trueVtxPos):
      if not isFiducial(trueVtxPos):
        continue

      try:
        xsecWeight[0] = weights.get(kpst.run, kpst.subrun, kpst.event)
      except:
        print("Couldn't find weight for run %i, subrun %i, event %i in %s!!!"%(kpst.run, kpst.subrun, kpst.event, args.weightfile))
        continue

      trueEnu[0] = nuInt.Nu().Momentum().E()
      trueCCNC[0] = nuInt.CCNC()
      trueNuPDG[0] = nuInt.Nu().PdgCode()
      passSelTruth[0] = 0
      trueNMus[0] = 0
      trueNThreshMus[0] = 0
      trueNCPis[0] = 0
      trueNThreshCPis[0] = 0
      trueNPrs[0] = 0
      trueNThreshPrs[0] = 0
      trueNPi0s[0] = 0
      trueNPhs[0] = 0
      trueNEls[0] = 0

      mcpg = ublarcvapp.mctools.MCPixelPGraph()
      mcpg.set_adc_treename("wire")
      mcpg.buildgraph(iolcv, ioll)

      for node in mcpg.node_v:

        if node.tid == node.mtid and node.origin == 1:

          if abs(node.pid) == 13:
            trueNMus[0] = trueNMus[0] + 1
            if node.E_MeV > 38.6:
              trueNThreshMus[0] = trueNThreshMus[0] + 1

          if abs(node.pid) == 211:
            trueNCPis[0] = trueNCPis[0] + 1
            if node.E_MeV > 51.0:
              trueNThreshCPis[0] = trueNThreshCPis[0] + 1

          if abs(node.pid) == 2212:
            trueNPrs[0] = trueNPrs[0] + 1
            if node.E_MeV > 343.0:
              trueNThreshPrs[0] = trueNThreshPrs[0] + 1

          if abs(node.pid) == 111:
            trueNPi0s[0] = trueNPi0s[0] + 1

          if abs(node.pid) == 22:
            trueNPhs[0] = trueNPhs[0] + 1

          if abs(node.pid) == 11:
            trueNEls[0] = trueNEls[0] + 1


      #if trueNPi0s[0] == 0 and trueNPhs[0] == 0 and trueNThreshMus[0] == 0 and trueNThreshCPis[0] == 0 and trueNThreshPrs[0] == 0:
      truthSingleShower = (trueNPhs[0] == 1 and trueNEls[0] == 0) or (trueNPhs[0] == 0 and trueNEls[0] == 1)
      if truthSingleShower and trueNPi0s[0] == 0 and trueNThreshMus[0] == 0 and trueNThreshCPis[0] == 0 and trueNThreshPrs[0] == 0:
        passSelTruth[0] = 1


    passSelReco[0] = 0
    recoSingleVtx[0] = 0
    recoNNegLLTrks[0] = -1
    recoNNegLLThreshTrks[0] = -1
    recoNPosLLTrks[0] = -1
    recoNPosLLThreshTrks[0] = -1
    recoVtxDistToTrue[0] = -99.
    recoShowerNHits[0] = -1

    nSSvertices = 0
    nMSvertices = 0
    for vtx in vertices:
      recoVtxPos = rt.TVector3(vtx.pos[0], vtx.pos[1], vtx.pos[2])
      if vtx.shower_v.size() == 1 and isFiducial(recoVtxPos):
        nSSvertices = nSSvertices + 1
        showerVertex = vtx
      if vtx.shower_v.size() > 1:
        nMSvertices = nMSvertices + 1

    if nSSvertices != 1 or nMSvertices != 0:
      eventTree.Fill()
      continue

    recoSingleVtx[0] = 1
    recoNNegLLTrks[0] = 0
    recoNNegLLThreshTrks[0] = 0
    recoNPosLLTrks[0] = 0
    recoNPosLLThreshTrks[0] = 0
    recoVtxDistToTrue[0] = getVertexDistance(trueVtxPos, showerVertex)
    recoShowerNHits[0] = showerVertex.shower_v[0].size()

    for iT in range(showerVertex.track_v.size()):
      length = showerVertex.track_v[iT].Length()
      llr = showerVertex.track_mu_vs_proton_llratio_v.at(iT)
      if llr < 0.0:
        recoNNegLLTrks[0] = recoNNegLLTrks[0] + 1
        if length > 65.:
          recoNNegLLThreshTrks[0] = recoNNegLLThreshTrks[0] + 1
      else:
        recoNPosLLTrks[0] = recoNPosLLTrks[0] + 1
        if length > 7.:
          recoNPosLLThreshTrks[0] = recoNPosLLThreshTrks[0] + 1

    if recoNNegLLThreshTrks[0] == 0 and recoNPosLLThreshTrks[0] == 0:
      passSelReco[0] = 1

    if args.printTrueEvents and passSelTruth[0] == 1:
      mcpg.printAllNodeInfo()
      printStatement = filepair[0]+" "+filepair[1]+" (entry, run, subrun, event): (%i, %i, %i, %i)"%(ientry, kpst.run, kpst.subrun, kpst.event)
      if passSelReco[0] == 1:
        print("PASSED: "+printStatement)
      else:
        print("FAILED: "+printStatement)

    if args.printRecoEvents and passSelReco[0] == 1:
      mcpg.printAllNodeInfo()
      printStatement = filepair[0]+" "+filepair[1]+" (entry, run, subrun, event): (%i, %i, %i, %i)"%(ientry, kpst.run, kpst.subrun, kpst.event)
      if passSelTruth[0] == 1:
        print("TRUE SIGNAL: "+printStatement)
      else:
        print("TRUE BACKGROUND: "+printStatement)
      print("reco shower vertex pos: %f, %f %f"%(showerVertex.pos[0],showerVertex.pos[1],showerVertex.pos[2]))

    eventTree.Fill()

  #++++++ end entry loop ++++++++++++++++++++++++++++++++++++++++++++++++++++=

  ioll.close()
  iolcv.finalize()
  kpsfile.Close()

#-------- end file loop -----------------------------------------------------#


if args.isMC:
  totPOT[0] = totPOT_
  totGoodPOT[0] = totGoodPOT_
  potTree.Fill()

outRootFile.cd()
eventTree.Write("",rt.TObject.kOverwrite)
if args.isMC:
  potTree.Write("",rt.TObject.kOverwrite)
outRootFile.Close()


