
import os,sys,argparse
import ROOT as rt
from larlite import larlite
from larlite import larutil
from ublarcvapp import ublarcvapp
from larcv import larcv
from larflow import larflow
from array import array
from helpers.larflowreco_ana_funcs import *

parser = argparse.ArgumentParser("Test mB Excess Selection")
parser.add_argument("-f", "--files", required=True, type=str, nargs="+", help="input kpsreco files")
parser.add_argument("-t", "--truth", required=True, type=str, help="text file containing merged_dlreco list")
parser.add_argument("-o", "--outfile", type=str, default="test_mBexcess_selection_output.root", help="output file name")
parser.add_argument("--printEvents", help="print out event info for failed and passed events", action="store_true")
parser.add_argument("--oldVtxBranch", help="use nufitted_v instead of nuvetoed_v for old reco", action="store_true")
args = parser.parse_args()

files = getFiles("merged_dlreco_", args.files, args.truth)

outRootFile = rt.TFile(args.outfile, "RECREATE")
eventTree = rt.TTree("EventTree","EventTree")
trueEnue = array('f', [0.])
passSelTruth = array('i', [0])
passSelReco = array('i', [0])
truthNMus = array('i', [0])
truthNThreshMus = array('i', [0])
truthNCPis = array('i', [0])
truthNThreshCPis = array('i', [0])
truthNPrs = array('i', [0])
truthNThreshPrs = array('i', [0])
truthNPi0s = array('i', [0])
truthNPhs = array('i', [0])
recoSingleVtx = array('i', [0])
recoNNegLLTrks = array('i', [0])
recoNNegLLThreshTrks = array('i', [0])
recoNPosLLTrks = array('i', [0])
recoNPosLLThreshTrks = array('i', [0])
eventTree.Branch("trueEnue", trueEnue, 'trueEnue/F')
eventTree.Branch("passSelTruth", passSelTruth, 'passSelTruth/I')
eventTree.Branch("passSelReco", passSelReco, 'passSelReco/I')
eventTree.Branch("truthNMus", truthNMus, 'truthNMus/I')
eventTree.Branch("truthNThreshMus", truthNThreshMus, 'truthNThreshMus/I')
eventTree.Branch("truthNCPis", truthNCPis, 'truthNCPis/I')
eventTree.Branch("truthNThreshCPis", truthNThreshCPis, 'truthNThreshCPis/I')
eventTree.Branch("truthNPrs", truthNPrs, 'truthNPrs/I')
eventTree.Branch("truthNThreshPrs", truthNThreshPrs, 'truthNThreshPrs/I')
eventTree.Branch("truthNPi0s", truthNPi0s, 'truthNPi0s/I')
eventTree.Branch("truthNPhs", truthNPhs, 'truthNPhs/I')
eventTree.Branch("recoSingleVtx", recoSingleVtx, 'recoSingleVtx/I')
eventTree.Branch("recoNNegLLTrks", recoNNegLLTrks, 'recoNNegLLTrks/I')
eventTree.Branch("recoNNegLLThreshTrks", recoNNegLLThreshTrks, 'recoNNegLLThreshTrks/I')
eventTree.Branch("recoNPosLLTrks", recoNPosLLTrks, 'recoNPosLLTrks/I')
eventTree.Branch("recoNPosLLThreshTrks", recoNPosLLThreshTrks, 'recoNPosLLThreshTrks/I')


sce = larutil.SpaceChargeMicroBooNE()
mcNuVertexer = ublarcvapp.mctools.NeutrinoVertex()


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

    mctruth = ioll.get_data(larlite.data.kMCTruth, "generator")
    nuInt = mctruth.at(0).GetNeutrino()
    lep = nuInt.Lepton()
    mcNuVertex = mcNuVertexer.getPos3DwSCE(ioll, sce)
    trueVtxPos = rt.TVector3(mcNuVertex[0], mcNuVertex[1], mcNuVertex[2])

    if nuInt.CCNC() != 0 or lep.PdgCode() != 11 or not isFiducial(trueVtxPos):
      continue

    trueEnue[0] = nuInt.Nu().Momentum().E()
    passSelTruth[0] = 0
    truthNMus[0] = 0
    truthNThreshMus[0] = 0
    truthNCPis[0] = 0
    truthNThreshCPis[0] = 0
    truthNPrs[0] = 0
    truthNThreshPrs[0] = 0
    truthNPi0s[0] = 0
    truthNPhs[0] = 0

    mcpg = ublarcvapp.mctools.MCPixelPGraph()
    mcpg.set_adc_treename("wire")
    mcpg.buildgraph(iolcv, ioll)

    for node in mcpg.node_v:

      if node.tid == node.mtid and node.origin == 1:

        if abs(node.pid) == 13:
          truthNMus[0] = truthNMus[0] + 1
          if node.E_MeV > 38.6:
            truthNThreshMus[0] = truthNThreshMus[0] + 1

        if abs(node.pid) == 211:
          truthNCPis[0] = truthNCPis[0] + 1
          if node.E_MeV > 51.0:
            truthNThreshCPis[0] = truthNThreshCPis[0] + 1

        if abs(node.pid) == 2212:
          truthNPrs[0] = truthNPrs[0] + 1
          if node.E_MeV > 343.0:
            truthNThreshPrs[0] = truthNThreshPrs[0] + 1

        if abs(node.pid) == 111:
          truthNPi0s[0] = truthNPi0s[0] + 1

        if abs(node.pid) == 22:
          truthNPhs[0] = truthNPhs[0] + 1


    if truthNPi0s[0] == 0 and truthNPhs[0] == 0 and truthNThreshMus[0] == 0 and truthNThreshCPis[0] == 0 and truthNThreshPrs[0] == 0:
      passSelTruth[0] = 1


    passSelReco[0] = 0
    recoSingleVtx[0] = 0
    recoNNegLLTrks[0] = -1
    recoNNegLLThreshTrks[0] = -1
    recoNPosLLTrks[0] = -1
    recoNPosLLThreshTrks[0] = -1

    nSSvertices = 0
    for vtx in vertices:
      recoVtxPos = rt.TVector3(vtx.pos[0], vtx.pos[1], vtx.pos[2])
      if vtx.shower_v.size() == 1 and isFiducial(recoVtxPos):
        nSSvertices = nSSvertices + 1
        showerVertex = vtx

    if nSSvertices != 1:
      eventTree.Fill()
      continue

    recoSingleVtx[0] = 1
    recoNNegLLTrks[0] = 0
    recoNNegLLThreshTrks[0] = 0
    recoNPosLLTrks[0] = 0
    recoNPosLLThreshTrks[0] = 0

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

    if args.printEvents and passSelTruth[0] == 1:
      printStatement = filepair[0]+" "+filepair[1]+" (entry, run, subrun, event): (%i, %i, %i, %i)"%(ientry, kpst.run, kpst.subrun, kpst.event)
      if passSelReco[0] == 1:
        print("PASSED: "+printStatement)
      else:
        print("FAILED: "+printStatement)

    eventTree.Fill()

  #++++++ end entry loop ++++++++++++++++++++++++++++++++++++++++++++++++++++=

  ioll.close()
  iolcv.finalize()
  kpsfile.Close()

#-------- end file loop -----------------------------------------------------#


outRootFile.cd()
eventTree.Write("",rt.TObject.kOverwrite)
outRootFile.Close()


