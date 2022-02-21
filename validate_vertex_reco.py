
import os,sys,argparse
import ROOT as rt
from larlite import larlite
from larlite import larutil
from ublarcvapp import ublarcvapp
from larcv import larcv
from larflow import larflow
from math import sqrt as sqrt
from math import acos as acos
from math import pi
from array import array
from event_weighting.event_weight_helper import SumPOT, Weights

parser = argparse.ArgumentParser("Validate Vertex Reco")
parser.add_argument("-f", "--files", required=True, type=str, nargs="+", help="input kpsreco files")
parser.add_argument("-t", "--truth", required=True, type=str, help="text file containing merged_dlreco list")
parser.add_argument("-w", "--weightfile", required=True, type=str, help="weights file (pickled python dict)")
parser.add_argument("-o", "--outfile", type=str, default="validate_vertex_reco_output.root", help="output file name")
parser.add_argument("--contained", help="require primary muon containment for numu events", action="store_true")
args = parser.parse_args()

files = []
for kpsfile in args.files:
  dlrecofilelist = open(args.truth, "r")
  for line in dlrecofilelist:
    dlrecofile = line.replace("\n","")
    samtag = dlrecofile[dlrecofile.find("merged_dlreco_"):].replace("merged_dlreco_","").replace(".root","")
    if samtag in kpsfile:
      files.append([kpsfile, dlrecofile])
      break
  dlrecofilelist.close()


def MCLeptonOkay(nuIntLep, mcObjLep):
  if mcObjLep.PdgCode() != nuIntLep.PdgCode():
    return False
  nuIntStart = nuIntLep.Position()
  mcObjStart = mcObjLep.Start()
  xdiffSq = (nuIntStart.X() - mcObjStart.X())**2
  ydiffSq = (nuIntStart.Y() - mcObjStart.Y())**2
  zdiffSq = (nuIntStart.Z() - mcObjStart.Z())**2
  if sqrt(xdiffSq + ydiffSq + zdiffSq) > 1.:
    return False
  return True

detCrds = [[0., 256.35], [-116.5, 116.5], [0, 1036.8]]
fidCrds = [ [detCrds[0][0] + 20. , detCrds[0][1] - 20.] ]
fidCrds.append( [detCrds[1][0] + 20. , detCrds[1][1] - 20.] )
fidCrds.append( [detCrds[2][0] + 20. , detCrds[2][1] - 60.] )

def inRange(x, bnd):
  return (x >= bnd[0] and x <= bnd[1])

def isFiducial(p):
  return (inRange(p.X(),fidCrds[0]) and inRange(p.Y(),fidCrds[1]) and inRange(p.Z(),fidCrds[2]))

def isInDetector(p):
  return (inRange(p.X(),detCrds[0]) and inRange(p.Y(),detCrds[1]) and inRange(p.Z(),detCrds[2]))

def getVertexDistance(pos3v, recoVtx):
  xdiffSq = (pos3v.X() - recoVtx.pos[0])**2
  ydiffSq = (pos3v.Y() - recoVtx.pos[1])**2
  zdiffSq = (pos3v.Z() - recoVtx.pos[2])**2
  return sqrt( xdiffSq + ydiffSq + zdiffSq )

def getTrackGap(track_v, recoVtx):
  maxLength = -99.
  for track in track_v:
    if track.Length() > maxLength:
      maxLength = track.Length()
      gap = getVertexDistance(track.Vertex(), recoVtx)
  return gap

def getShowerGap(shower_trunk_v, recoVtx):
  minGap = 1e9
  for shower_trunk in shower_trunk_v:
    gap = getVertexDistance(shower_trunk.Vertex(), recoVtx)
    if gap < minGap:
      minGap = gap
  return minGap

def getLeptonPixels(pdg, ioll, iolcv):
  # get lepton node
  mcpg = ublarcvapp.mctools.MCPixelPGraph()
  mcpg.set_adc_treename("wire")
  mcpg.buildgraph(iolcv, ioll)
  foundLepNode = False
  for node in mcpg.node_v:
    if node.pid == pdg and node.tid == node.mtid:
      lepNode = node
      foundLepNode = True
      break
  if not foundLepNode:
    sys.exit("Couldn't find lepton node in MCPixelPGraph!!")
  # get pixel value array
  evtImage2D = iolcv.get_data(larcv.kProductImage2D, "wire")
  image2Dvec = evtImage2D.Image2DArray()
  # calculate pixel info return variables
  totLepPixI = 0.
  lepTickLists = [ [], [], [] ]
  lepPixelDictList = [ {}, {}, {} ]
  for p in range(3):
    lepPix = lepNode.pix_vv[p]
    for iP in range(lepPix.size()//2):
      simTick = lepPix[2*iP]
      simWire = lepPix[2*iP+1]
      row = (simTick - 2400)//6
      totLepPixI = totLepPixI + image2Dvec[p].pixel(row,simWire)
      if simTick not in lepTickLists[p]:
        lepTickLists[p].append(simTick)
        lepPixelDictList[p][simTick] = [simWire]
      else:
        lepPixelDictList[p][simTick].append(simWire)
  return totLepPixI, lepTickLists, lepPixelDictList

def getBestCompleteness(iolcv, vertex, lepPDG, totLepPixI, lepTickLists, lepPixelDictList):
  bestComp = -1.
  evtImage2D = iolcv.get_data(larcv.kProductImage2D, "wire")
  image2Dvec = evtImage2D.Image2DArray()
  if lepPDG == 13:
    clusterVec = vertex.track_hitcluster_v
  if lepPDG == 11:
    clusterVec = vertex.shower_v
  for cluster in clusterVec:
    matchedPixels = []
    matchedUniquePixI = 0.
    totPixI = 0.
    for hit in cluster:
      for p in range(3):
        row = (hit.tick - 2400)//6
        pixVal = image2Dvec[p].pixel(row, hit.targetwire[p])
        totPixI = totPixI + pixVal
        if hit.tick in lepTickLists[p]:
          if hit.targetwire[p] in lepPixelDictList[p][hit.tick]:
            matchedPixel = [ p, hit.tick, hit.targetwire[p] ]
            if matchedPixel not in matchedPixels:
              matchedUniquePixI = matchedUniquePixI + pixVal
              matchedPixels.append(matchedPixel)
    comp = matchedUniquePixI/totLepPixI
    if comp > bestComp:
      bestComp = comp
  return bestComp


#NuSelectionVariables classes
prongvars = larflow.reco.NuSelProngVars()
#vertexvars = larflow.reco.NuSelVertexVars()
showertrunkvars = larflow.reco.NuSelShowerTrunkAna()
wcoverlapvars = larflow.reco.NuSelWCTaggerOverlap()
showergapana2d = larflow.reco.NuSelShowerGapAna2D()
#unrecocharge = larflow.reco.NuSelUnrecoCharge()
cosmictagger = larflow.reco.NuSelCosmicTagger()
#muvsproton = larflow.reco.TrackForwardBackwardLL()
showerTrunkCheck = larflow.reco.NuVertexShowerTrunkCheck()

#other needed classes
sce = larutil.SpaceChargeMicroBooNE()
mcNuVertexer = ublarcvapp.mctools.NeutrinoVertex()
pidAlgo = larflow.reco.LikelihoodProtonMuon()

def CalculateNuSelectionVariables(nuvtx, iolcv, ioll):
  nusel = larflow.reco.NuSelectionVariables()
  showergapana2d.analyze(iolcv, ioll, nuvtx, nusel)
  if nusel.nplanes_connected >= 2:
    showerTrunkCheck.checkNuCandidateProngsForMissingCharge(nuvtx, iolcv, ioll)
  nusel.max_proton_pid = 1e3
  for lltrack in nuvtx.track_v:
    trackvars = larflow.reco.NuSelectionVariables.TrackVar_t()
    try:
      trackvars.proton_ll = pidAlgo.calculateLL(lltrack, nuvtx.pos)
      if trackvars.proton_ll < nusel.max_proton_pid:
        nusel.max_proton_pid = trackvars.proton_ll
    except:
      pass
  prongvars.analyze(nuvtx, nusel)
  showertrunkvars.analyze(nuvtx, nusel, iolcv, ioll)
  #vertexvars.analyze(iolcv, ioll, nuvtx, nusel)
  wcoverlapvars.analyze(nuvtx, nusel, iolcv)
  #unrecocharge.analyze(iolcv, ioll, nuvtx, nusel)
  cosmictagger.analyze(nuvtx, nusel)
  #muvsproton.analyze(nuvtx, nusel)
  return nusel


outRootFile = rt.TFile(args.outfile, "RECREATE")

potTree = rt.TTree("potTree","potTree")
totPOT = array('f', [0.])
totGoodPOT = array('f', [0.])
potTree.Branch("totPOT", totPOT, 'totPOT/F')
potTree.Branch("totGoodPOT", totGoodPOT, 'totGoodPOT/F')

trackTree = rt.TTree("TrackTree","TrackTree")
dQdxStatus = array('i', [0])
trackLength = array('f', [0.])
trackTree.Branch("dQdxStatus", dQdxStatus, 'dQdxStatus/I')
trackTree.Branch("trackLength", trackLength, 'trackLength/F')

vertexTree = rt.TTree("VertexTree","VertexTree")
maxNVtxs = 1000
xsecWeight = array('f', [0.])
trueNuE = array('f', [0.])
trueLepE = array('f', [0.])
trueLepPDG = array('i', [0])
bestRecoComp = array('f', [0.])
nVertices = array('i', [0])
vtxDistToTrue = array('f', maxNVtxs*[0.])
vtxBestComp = array('f', maxNVtxs*[0.])
vtxScore = array('f', maxNVtxs*[0.])
vtxAvgScore = array('f', maxNVtxs*[0.])
vtxMaxScore = array('f', maxNVtxs*[0.])
vtxNetScore = array('f', maxNVtxs*[0.])
vtxNetNuScore = array('f', maxNVtxs*[0.])
vtxGap = array('f', maxNVtxs*[0.])
vtxKeyPtType = array('i', maxNVtxs*[0])
vtxNTracks = array('i', maxNVtxs*[0])
vtxNShowers = array('i', maxNVtxs*[0])
vtxHasReco = array('i', maxNVtxs*[0])
vtxMaxTrkTotll = array('f', maxNVtxs*[0])
vtxMaxTrkTotllmuon = array('f', maxNVtxs*[0])
vtxMaxTrkTotllproton = array('f', maxNVtxs*[0])
vtxLongTrkTotll = array('f', maxNVtxs*[0])
vtxLongTrkTotllmuon = array('f', maxNVtxs*[0])
vtxLongTrkTotllproton = array('f', maxNVtxs*[0])
nusel_max_proton_pid = array('f', maxNVtxs*[0])
nusel_ntracks = array('i', maxNVtxs*[0])
nusel_nshowers = array('i', maxNVtxs*[0])
nusel_max_shower_length = array('f', maxNVtxs*[0])
nusel_max_track_length = array('f', maxNVtxs*[0])
nusel_max_shower_nhits = array('i', maxNVtxs*[0])
nusel_max_track_nhits = array('i', maxNVtxs*[0])
nusel_min_shower_gap = array('f', maxNVtxs*[0])
nusel_max_shower_gap = array('f', maxNVtxs*[0])
nusel_min_track_gap = array('f', maxNVtxs*[0])
nusel_max_track_gap = array('f', maxNVtxs*[0])
nusel_closest_shower_ll = array('f', maxNVtxs*[0]) 
nusel_largest_shower_ll = array('f', maxNVtxs*[0]) 
nusel_closest_shower_avedqdx = array('f', maxNVtxs*[0]) 
nusel_largest_shower_avedqdx = array('f', maxNVtxs*[0])
nusel_nplanes_connected = array('i', maxNVtxs*[0])
#nusel_vertex_hip_fraction = array('f', maxNVtxs*[0])
#nusel_vertex_charge_per_pixel = array('f', maxNVtxs*[0])
#nusel_vertex_type = array('i', maxNVtxs*[0])
nusel_frac_trackhits_on_cosmic = array('f', maxNVtxs*[0])
nusel_frac_showerhits_on_cosmic = array('f', maxNVtxs*[0])
nusel_frac_allhits_on_cosmic = array('f', maxNVtxs*[0])
nusel_nshower_pts_on_cosmic = array('i', maxNVtxs*[0])
nusel_ntrack_pts_on_cosmic = array('i', maxNVtxs*[0])
vertexTree.Branch("xsecWeight", xsecWeight, 'xsecWeight/F')
vertexTree.Branch("trueNuE", trueNuE, 'trueNuE/F')
vertexTree.Branch("trueLepE", trueLepE, 'trueLepE/F')
vertexTree.Branch("trueLepPDG", trueLepPDG, 'trueLepPDG/I')
vertexTree.Branch("bestRecoComp", bestRecoComp, 'bestRecoComp/F')
vertexTree.Branch("nVertices", nVertices, 'nVertices/I')
vertexTree.Branch("vtxDistToTrue", vtxDistToTrue, 'vtxDistToTrue[nVertices]/F')
vertexTree.Branch("vtxBestComp", vtxBestComp, 'vtxBestComp[nVertices]/F')
vertexTree.Branch("vtxScore", vtxScore, 'vtxScore[nVertices]/F')
vertexTree.Branch("vtxAvgScore", vtxAvgScore, 'vtxAvgScore[nVertices]/F')
vertexTree.Branch("vtxMaxScore", vtxMaxScore, 'vtxMaxScore[nVertices]/F')
vertexTree.Branch("vtxNetScore", vtxNetScore, 'vtxNetScore[nVertices]/F')
vertexTree.Branch("vtxNetNuScore", vtxNetNuScore, 'vtxNetNuScore[nVertices]/F')
vertexTree.Branch("vtxGap", vtxGap, 'vtxGap[nVertices]/F')
vertexTree.Branch("vtxKeyPtType", vtxKeyPtType, 'vtxKeyPtType[nVertices]/I')
vertexTree.Branch("vtxNTracks", vtxNTracks, 'vtxNTracks[nVertices]/I')
vertexTree.Branch("vtxNShowers", vtxNShowers, 'vtxNShowers[nVertices]/I')
vertexTree.Branch("vtxHasReco", vtxHasReco, 'vtxHasReco[nVertices]/I')
vertexTree.Branch("vtxMaxTrkTotll", vtxMaxTrkTotll, 'vtxMaxTrkTotll[nVertices]/F')
vertexTree.Branch("vtxMaxTrkTotllmuon", vtxMaxTrkTotllmuon, 'vtxMaxTrkTotllmuon[nVertices]/F')
vertexTree.Branch("vtxMaxTrkTotllproton", vtxMaxTrkTotllproton, 'vtxMaxTrkTotllproton[nVertices]/F')
vertexTree.Branch("vtxLongTrkTotll", vtxLongTrkTotll, 'vtxLongTrkTotll[nVertices]/F')
vertexTree.Branch("vtxLongTrkTotllmuon", vtxLongTrkTotllmuon, 'vtxLongTrkTotllmuon[nVertices]/F')
vertexTree.Branch("vtxLongTrkTotllproton", vtxLongTrkTotllproton, 'vtxLongTrkTotllproton[nVertices]/F')
vertexTree.Branch("nusel_max_proton_pid", nusel_max_proton_pid, 'nusel_max_proton_pid[nVertices]/F')
vertexTree.Branch("nusel_ntracks", nusel_ntracks, 'nusel_ntracks[nVertices]/I')
vertexTree.Branch("nusel_nshowers", nusel_nshowers, 'nusel_nshowers[nVertices]/I')
vertexTree.Branch("nusel_max_shower_length", nusel_max_shower_length, 'nusel_max_shower_length[nVertices]/F')
vertexTree.Branch("nusel_max_track_length", nusel_max_track_length, 'nusel_max_track_length[nVertices]/F')
vertexTree.Branch("nusel_max_shower_nhits", nusel_max_shower_nhits, 'nusel_max_shower_nhits[nVertices]/I')
vertexTree.Branch("nusel_max_track_nhits", nusel_max_track_nhits, 'nusel_max_track_nhits[nVertices]/I')
vertexTree.Branch("nusel_min_shower_gap", nusel_min_shower_gap, 'nusel_min_shower_gap[nVertices]/F')
vertexTree.Branch("nusel_max_shower_gap", nusel_max_shower_gap, 'nusel_max_shower_gap[nVertices]/F')
vertexTree.Branch("nusel_min_track_gap", nusel_min_track_gap, 'nusel_min_track_gap[nVertices]/F')
vertexTree.Branch("nusel_max_track_gap", nusel_max_track_gap, 'nusel_max_track_gap[nVertices]/F')
vertexTree.Branch("nusel_closest_shower_ll", nusel_closest_shower_ll, 'nusel_closest_shower_ll[nVertices]/F')
vertexTree.Branch("nusel_largest_shower_ll", nusel_largest_shower_ll, 'nusel_largest_shower_ll[nVertices]/F')
vertexTree.Branch("nusel_closest_shower_avedqdx", nusel_closest_shower_avedqdx, 'nusel_closest_shower_avedqdx[nVertices]/F')
vertexTree.Branch("nusel_largest_shower_avedqdx", nusel_largest_shower_avedqdx, 'nusel_largest_shower_avedqdx[nVertices]/F')
vertexTree.Branch("nusel_nplanes_connected", nusel_nplanes_connected, 'nusel_nplanes_connected[nVertices]/I')
#vertexTree.Branch("nusel_vertex_hip_fraction", nusel_vertex_hip_fraction, 'nusel_vertex_hip_fraction[nVertices]/F')
#vertexTree.Branch("nusel_vertex_charge_per_pixel", nusel_vertex_charge_per_pixel, 'nusel_vertex_charge_per_pixel[nVertices]/F')
#vertexTree.Branch("nusel_vertex_type", nusel_vertex_type, 'nusel_vertex_type[nVertices]/I')
vertexTree.Branch("nusel_frac_trackhits_on_cosmic", nusel_frac_trackhits_on_cosmic, 'nusel_frac_trackhits_on_cosmic[nVertices]/F')
vertexTree.Branch("nusel_frac_showerhits_on_cosmic", nusel_frac_showerhits_on_cosmic, 'nusel_frac_showerhits_on_cosmic[nVertices]/F')
vertexTree.Branch("nusel_frac_allhits_on_cosmic", nusel_frac_allhits_on_cosmic, 'nusel_frac_allhits_on_cosmic[nVertices]/F')
vertexTree.Branch("nusel_nshower_pts_on_cosmic", nusel_nshower_pts_on_cosmic, 'nusel_nshower_pts_on_cosmic[nVertices]/I')
vertexTree.Branch("nusel_ntrack_pts_on_cosmic", nusel_ntrack_pts_on_cosmic, 'nusel_ntrack_pts_on_cosmic[nVertices]/I')


track_total_c = 0
track_good_c = 0
track_missingdQdx_c = 0

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

  potInFile, goodPotInFile = SumPOT(filepair[1])
  totPOT_ = totPOT_ + potInFile
  totGoodPOT_ = totGoodPOT_ + goodPotInFile

  #++++++ begin entry loop ++++++++++++++++++++++++++++++++++++++++++++++++++++=
  for ientry in range(ioll.get_entries()):
  
    ioll.go_to(ientry)
    iolcv.read_entry(ientry)
    kpst.GetEntry(ientry)
  
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

    if nuInt.CCNC() != 0 or lep.PdgCode() not in [11,13] or not isFiducial(trueVtxPos):
      continue

    try:
      xsecWeight[0] = weights.get(kpst.run, kpst.subrun, kpst.event)
    except:
      print("Couldn't find weight for run %i, subrun %i, event %i in %s!!!"%(kpst.run, kpst.subrun, kpst.event, args.weightfile))
      continue

    lepPDG = lep.PdgCode()

    if lepPDG == 13:
      mcleptons = ioll.get_data(larlite.data.kMCTrack, "mcreco")
    if lepPDG == 11:
      mcleptons = ioll.get_data(larlite.data.kMCShower, "mcreco")
    for mclepton in mcleptons:
      if mclepton.PdgCode() == lepPDG and mclepton.Process() == 'primary':
        mcLeptonUnCorr = mclepton
        break
  
    if not MCLeptonOkay(lep, mcLeptonUnCorr):
      print("Couldn't find MC lepton match!!!")
      continue

    if lepPDG == 13 and args.contained and not isInDetector(mcLeptonUnCorr.End()):
        continue

    #if lepPDG == 11 and not isFiducial(mcLeptonUnCorr.End()):
    #  continue

    trueNuE[0] = nuInt.Nu().Momentum().E()
    trueLepE[0] = lep.Momentum().E()
    trueLepPDG[0] = lepPDG
    nVertices[0] = kpst.nufitted_v.size()

    totLepPixI, lepTickLists, lepPixelDictList = getLeptonPixels(lepPDG, ioll, iolcv)

    iV = 0
    bestRecoComp[0] = -1.
    #------- begin vertex loop ------------------------------------------------#
    for vertex in kpst.nufitted_v:

      c = getBestCompleteness(iolcv, vertex, lepPDG, totLepPixI, lepTickLists, lepPixelDictList)
      if c > bestRecoComp[0]:
        bestRecoComp[0] = c
      vtxBestComp[iV] = c
      vtxDistToTrue[iV] = getVertexDistance(trueVtxPos, vertex)
      vtxScore[iV] = vertex.score
      vtxAvgScore[iV] = vertex.avgScore
      vtxMaxScore[iV] = vertex.maxScore
      vtxNetScore[iV] = vertex.netScore
      vtxNetNuScore[iV] = vertex.netNuScore
      vtxKeyPtType[iV] = vertex.keypoint_type
      vtxNTracks[iV] = vertex.track_v.size()
      vtxNShowers[iV] = vertex.shower_v.size()
      vtxHasReco[iV] = 0
      vtxGap[iV] = -99.
      if lepPDG == 11 and vertex.shower_v.size() > 0:
        vtxHasReco[iV] = 1
        vtxGap[iV] = getShowerGap(vertex.shower_trunk_v, vertex)
      if lepPDG == 13 and vertex.track_v.size() > 0:
        vtxHasReco[iV] = 1
        vtxGap[iV] = getTrackGap(vertex.track_v, vertex)

      vtxMaxTrkTotll[iV] = -1.
      vtxMaxTrkTotllmuon[iV] = -1.
      vtxMaxTrkTotllproton[iV] = -1.
      vtxLongTrkTotll[iV] = -1.
      vtxLongTrkTotllmuon[iV] = -1.
      vtxLongTrkTotllproton[iV] = -1.
      maxTrkLength = -99.
      for track in vertex.track_v:
        track_total_c = track_total_c + 1
        try:
          llpid = pidAlgo.calculateLLseparate(track, vertex.pos)
          if llpid[0] > vtxMaxTrkTotll[iV]:
            vtxMaxTrkTotll[iV] = llpid[0]
          if llpid[1] > vtxMaxTrkTotllproton[iV]:
            vtxMaxTrkTotllproton[iV] = llpid[1]
          if llpid[2] > vtxMaxTrkTotllmuon[iV]:
            vtxMaxTrkTotllmuon[iV] = llpid[2]
          if track.Length() > maxTrkLength:
            maxTrkLength = track.Length()
            vtxLongTrkTotll[iV] = llpid[0]
            vtxLongTrkTotllproton[iV] = llpid[1]
            vtxLongTrkTotllmuon[iV] = llpid[2]
          track_good_c = track_good_c + 1
          dQdxStatus[0] = 0
        except:
          track_missingdQdx_c = track_missingdQdx_c + 1
          dQdxStatus[0] = 1
        trackLength[0] = track.Length()
        trackTree.Fill()

      selVars = CalculateNuSelectionVariables(vertex, iolcv, ioll)
      nusel_max_proton_pid[iV] = selVars.max_proton_pid
      nusel_ntracks[iV] = selVars.ntracks
      nusel_nshowers[iV] = selVars.nshowers
      nusel_max_shower_length[iV] = selVars.max_shower_length
      nusel_max_track_length[iV] = selVars.max_track_length
      nusel_max_shower_nhits[iV] = selVars.max_shower_nhits
      nusel_max_track_nhits[iV] = selVars.max_track_nhits
      nusel_min_shower_gap[iV] = selVars.min_shower_gap
      nusel_max_shower_gap[iV] = selVars.max_shower_gap
      nusel_min_track_gap[iV] = selVars.min_track_gap
      nusel_max_track_gap[iV] = selVars.max_track_gap
      nusel_closest_shower_ll[iV] = selVars.closest_shower_ll
      nusel_largest_shower_ll[iV] = selVars.largest_shower_ll
      nusel_closest_shower_avedqdx[iV] = selVars.closest_shower_avedqdx
      nusel_largest_shower_avedqdx[iV] = selVars.largest_shower_avedqdx
      nusel_nplanes_connected[iV] = selVars.nplanes_connected
      #nusel_vertex_hip_fraction[iV] = selVars.vertex_hip_fraction
      #nusel_vertex_charge_per_pixel[iV] = selVars.vertex_charge_per_pixel
      #nusel_vertex_type[iV] = selVars.vertex_type
      nusel_frac_trackhits_on_cosmic[iV] = selVars.frac_trackhits_on_cosmic
      nusel_frac_showerhits_on_cosmic[iV] = selVars.frac_showerhits_on_cosmic
      nusel_frac_allhits_on_cosmic[iV] = selVars.frac_allhits_on_cosmic
      nusel_nshower_pts_on_cosmic[iV] = selVars.nshower_pts_on_cosmic
      nusel_ntrack_pts_on_cosmic[iV] = selVars.ntrack_pts_on_cosmic
      iV = iV + 1
    #------- end vertex loop --------------------------------------------------#

    vertexTree.Fill()

  #++++++ end entry loop ++++++++++++++++++++++++++++++++++++++++++++++++++++=

  ioll.close()
  iolcv.finalize()
  kpsfile.Close()

#-------- end file loop -----------------------------------------------------#


totPOT[0] = totPOT_
totGoodPOT[0] = totGoodPOT_
potTree.Fill()

outRootFile.cd()
vertexTree.Write("",rt.TObject.kOverwrite)
trackTree.Write("",rt.TObject.kOverwrite)
potTree.Write("",rt.TObject.kOverwrite)
outRootFile.Close()


print("track counters (total, good, missing dQdx): %i, %i, %i"%(track_total_c, track_good_c, track_missingdQdx_c))

