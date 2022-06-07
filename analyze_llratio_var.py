
import os,sys,argparse
import ROOT as rt
from larlite import larlite
from larlite import larutil
from ublarcvapp import ublarcvapp
from larcv import larcv
from larflow import larflow
from array import array
from helpers.larflowreco_ana_funcs import *

parser = argparse.ArgumentParser("Analyze LL Ratio Track Var.")
parser.add_argument("-f", "--files", required=True, type=str, nargs="+", help="input kpsreco files")
parser.add_argument("-t", "--truth", required=True, type=str, help="text file containing merged_dlreco list")
parser.add_argument("-o", "--outfile", type=str, default="analyze_llratio_var_output.root", help="output file directory")
parser.add_argument("-v", "--vertexRange", type=float, default=3., help="search within this distance in cm of true interaction vertex")
parser.add_argument("-m", "--minMatchComp", type=float, default=0.5, help="minimum track completeness for truth matching")
parser.add_argument("--oldVtxBranch", help="use nufitted_v instead of nuvetoed_v for old reco", action="store_true")
args = parser.parse_args()

files = getFiles("merged_dlreco_", args.files, args.truth)

outRootFile = rt.TFile(args.outfile, "RECREATE")
trackTree = rt.TTree("TrackTree","TrackTree")
trackPDG = array('i',[0])
llRatio = array('f', [0.])
trackTree.Branch("trackPDG", trackPDG, 'trackPDG/I')
trackTree.Branch("llRatio", llRatio, 'llRatio/F')

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
    mcNuVertex = mcNuVertexer.getPos3DwSCE(ioll, sce)
    trueVtxPos = rt.TVector3(mcNuVertex[0], mcNuVertex[1], mcNuVertex[2])

    if nuInt.CCNC() != 0 or not isFiducial(trueVtxPos):
      continue

    closeVertices = []
    for vtx in vertices:
      deltaVertex = getVertexDistance(trueVtxPos, vtx)
      if deltaVertex < args.vertexRange:
        closeVertices.append(vtx)

    nodes = []
    mcpg = ublarcvapp.mctools.MCPixelPGraph()
    mcpg.set_adc_treename("wire")
    mcpg.buildgraph(iolcv, ioll)
    for node in mcpg.node_v:
      if abs(node.pid) in [13, 211, 2212] and node.tid == node.mtid and node.E_MeV > 40.:
        nodes.append(node)

    evtImage2D = iolcv.get_data(larcv.kProductImage2D, "wire")
    image2Dvec = evtImage2D.Image2DArray()


    #++++++ begin node loop ++++++++++++++++++++++++++++++++++++++++++++++++++++=
    for node in nodes:

      totNodePixI = 0.
      nodeTickLists = [ [], [], [] ]
      nodePixelDictList = [ {}, {}, {} ]
      for p in range(3):
        nodePix = node.pix_vv[p]
        for iP in range(nodePix.size()//2):
          simTick = nodePix[2*iP]
          simWire = nodePix[2*iP+1]
          row = (simTick - 2400)//6
          totNodePixI = totNodePixI + image2Dvec[p].pixel(row,simWire)
          if simTick not in nodeTickLists[p]:
            nodeTickLists[p].append(simTick)
            nodePixelDictList[p][simTick] = [simWire]
          else:
            nodePixelDictList[p][simTick].append(simWire)

      if not totNodePixI > 0.:
        continue

      maxComp = -1
      ll = -999.

      #++++++ begin track loop ++++++++++++++++++++++++++++++++++++++++++++++++++=
      for vertex in closeVertices:
        for iT in range(vertex.track_v.size()):
          matchedPixels = []
          matchedUniquePixI = 0.
          #------- begin pixel loop --------------------------------------------------
          for hit in vertex.track_hitcluster_v[iT]:
            for p in range(3):
              row = (hit.tick - 2400)//6
              pixVal = image2Dvec[p].pixel(row, hit.targetwire[p])
              if hit.tick in nodeTickLists[p]:
                if hit.targetwire[p] in nodePixelDictList[p][hit.tick]:
                  matchedPixel = [ p, hit.tick, hit.targetwire[p] ]
                  if matchedPixel not in matchedPixels:
                    matchedUniquePixI = matchedUniquePixI + pixVal
                    matchedPixels.append(matchedPixel)
          #------- end pixel loop ----------------------------------------------------
          cI = matchedUniquePixI/totNodePixI
          if cI > maxComp:
            maxComp = cI
            ll = vertex.track_mu_vs_proton_llratio_v.at(iT)
      #++++++ end track loop ++++++++++++++++++++++++++++++++++++++++++++++++++++=
      if maxComp > args.minMatchComp:
        trackPDG[0] = node.pid
        llRatio[0] = ll
        trackTree.Fill()
    #++++++ end node loop ++++++++++++++++++++++++++++++++++++++++++++++++++++++=

  #++++++ end entry loop ++++++++++++++++++++++++++++++++++++++++++++++++++++=

  ioll.close()
  iolcv.finalize()
  kpsfile.Close()

#-------- end file loop -----------------------------------------------------#


outRootFile.cd()
trackTree.Write("",rt.TObject.kOverwrite)
outRootFile.Close()


