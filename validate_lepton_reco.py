
import os,sys,argparse
import ROOT as rt
from larlite import larlite
from larlite import larutil
from ublarcvapp import ublarcvapp
from larcv import larcv
from larflow import larflow
from math import sqrt, acos, cos, sin, pi
from array import array
from helpers.larflowreco_ana_funcs import *

parser = argparse.ArgumentParser("Validate Lepton Reco")
parser.add_argument("-f", "--files", required=True, type=str, nargs="+", help="input kpsreco files")
parser.add_argument("-t", "--truth", required=True, type=str, help="text file containing merged_dlreco list")
parser.add_argument("-o", "--outfile", type=str, default="validate_lepton_reco_output.root", help="output file directory")
parser.add_argument("-v", "--vertexRange", type=float, default=3., help="search within this distance in cm of true interaction vertex")
parser.add_argument("--nue", help="analyze nue events", action="store_true")
parser.add_argument("--numu", help="analyze numu events", action="store_true")
parser.add_argument("--contained", help="require primary muon containment for numu events", action="store_true")
parser.add_argument("--v1", help="reverse track points for v1 bug", action="store_true")
args = parser.parse_args()

if not args.nue ^ args.numu:
  sys.exit("Need to specify --nue or --numu. Exiting...")

files = getFiles("merged_dlreco_", args.files, args.truth)

def getStartDistance(trueStart, recoTrack, p=-1):
  if p < 0:
    trackStartPt = recoTrack.Vertex()
    if args.v1 and args.numu:
      trackStartPt = recoTrack.End()
    xdiffSq = (trueStart.X() - trackStartPt.X())**2
    ydiffSq = (trueStart.Y() - trackStartPt.Y())**2
    zdiffSq = (trueStart.Z() - trackStartPt.Z())**2
  else:
    xdiffSq = (trueStart.X() - recoTrack.LocationAtPoint(p).X())**2
    ydiffSq = (trueStart.Y() - recoTrack.LocationAtPoint(p).Y())**2
    zdiffSq = (trueStart.Z() - recoTrack.LocationAtPoint(p).Z())**2
  return sqrt( xdiffSq + ydiffSq + zdiffSq )

def getMuEndDistance(trueEnd, recoTrack):
  trackEndPt = recoTrack.End()
  if args.v1:
    trackEndPt = recoTrack.Vertex()
  xdiffSq = (trueEnd.X() - trackEndPt.X())**2
  ydiffSq = (trueEnd.Y() - trackEndPt.Y())**2
  zdiffSq = (trueEnd.Z() - trackEndPt.Z())**2
  return sqrt( xdiffSq + ydiffSq + zdiffSq )

def getHitDistance(trackPoint, hit):
  delXsq = (hit[0] - trackPoint.X())**2
  delYsq = (hit[1] - trackPoint.Y())**2
  delZsq = (hit[2] - trackPoint.Z())**2
  return sqrt(delXsq + delYsq + delZsq)

def getTrackSecFlags(trackStart, trackEnd, hit):
  flags = [0,0,0]
  if getHitDistance(trackStart, hit) < 10.:
    flags[0] = 1
  if getHitDistance(trackEnd, hit) < 10.:
    flags[2] = 1
  if not flags[0] and not flags[2]:
    flags[1] = 1
  return flags

def getAngle(true4p, recoTrack):
  trueMag = sqrt(true4p.Px()**2 + true4p.Py()**2 + true4p.Pz()**2)
  tPx = true4p.Px()/trueMag
  tPy = true4p.Py()/trueMag
  tPz = true4p.Pz()/trueMag
  recoStart = recoTrack.Vertex()
  if args.nue:
    if args.v1:
      anglePoint = recoTrack.End()
    else:
      rPx = recoTrack.VertexDirection().X()
      rPy = recoTrack.VertexDirection().Y()
      rPz = recoTrack.VertexDirection().Z()
  if args.numu:
    for i in range(recoTrack.NumberTrajectoryPoints()):
      anglePoint = recoTrack.LocationAtPoint(i)
      if getStartDistance(recoStart, recoTrack, i) > 5.:
        break
    if args.v1:
      recoStart = recoTrack.End()
      for i in reversed(range(recoTrack.NumberTrajectoryPoints())):
        anglePoint = recoTrack.LocationAtPoint(i)
        if getStartDistance(recoStart, recoTrack, i) > 5.:
          break
  if args.numu or (args.nue and args.v1):
    rPx = anglePoint.X() - recoStart.X()
    rPy = anglePoint.Y() - recoStart.Y()
    rPz = anglePoint.Z() - recoStart.Z()
  recoMag = sqrt(rPx**2 + rPy**2 + rPz**2)
  rPx = rPx/recoMag
  rPy = rPy/recoMag
  rPz = rPz/recoMag
  return acos(tPx*rPx + tPy*rPy + tPz*rPz)

def getAxisAngles(true4p):
  trueMag = sqrt(true4p.Px()**2 + true4p.Py()**2 + true4p.Pz()**2)
  tPx = true4p.Px()/trueMag
  tPy = true4p.Py()/trueMag
  tPz = true4p.Pz()/trueMag
  Uy = cos(pi/3.)
  Uz = sin(pi/3.)
  Vy = -sin(pi/6.)
  Vz = cos(pi/6.)
  thetaX = acos(tPx)
  thetaY = acos(tPy)
  thetaZ = acos(tPz)
  thetaU = acos(tPy*Uy + tPz*Uz)
  thetaV = acos(tPy*Vy + tPz*Vz)
  return thetaX, thetaY, thetaZ, thetaU, thetaV

def getMinAngle(angles):
  minAng = 999.
  for ang in angles:
    if ang > pi/2.:
      ang = pi - ang
    if ang < minAng:
      minAng = ang
  return minAng

def getLastPoint(trackUnCorr, trackCorr):
  iStep = trackUnCorr.NumberTrajectoryPoints() - 1
  while not isInDetector(trackUnCorr.LocationAtPoint(iStep)):
    iStep = iStep - 1
  return trackCorr.LocationAtPoint(iStep)



outRootFile = rt.TFile(args.outfile, "RECREATE")

nuIntTree = rt.TTree("NuIntTree","NuIntTree")
trueNuE = array('f', [0.])
trueLepE = array('f', [0.])
trueLepX = array('f', [0.])
trueLepY = array('f', [0.])
trueLepZ = array('f', [0.])
trueLepPx = array('f', [0.])
trueLepPy = array('f', [0.])
trueLepPz = array('f', [0.])
trueLepThX = array('f', [0.])
trueLepThY = array('f', [0.])
trueLepThZ = array('f', [0.])
trueLepThU = array('f', [0.])
trueLepThV = array('f', [0.])
trueLepMinWireAng = array('f', [0.])
trueLepMinSpecAng = array('f', [0.])
recoStatus = array('i',[0])
closestVtxDist = array('f', [0.])
vtxType = array('i', [0])
vtxDist = array('f', [0.])
vtxDist_x = array('f', [0.])
vtxDist_y = array('f', [0.])
vtxDist_z = array('f', [0.])
nCloseVertices05cm = array('i', [0])
nCloseVertices10cm = array('i', [0])
nCloseVertices15cm = array('i', [0])
trueLepCompletenessN = array('f', [0.])
trueLepCompletenessI = array('f', [0.])
recoLep1PixPurity = array('f', [0.])
recoLep2PixPurity = array('f', [0.])
recoLep3PixPurity = array('f', [0.])
recoLepPixIPurity = array('f', [0.])
recoLepPixIAnscPurity = array('f', [0.])
recoLepStartPosErr = array('f', [0.])
recoLepStartDirErr = array('f', [0.])
if args.numu:
  trueMuContained = array('i', [0])
  recoMuEndPosErr = array('f', [0.])
  recoMuTrackLength = array('f', [0.])
  recoMuEarlyPurity = array('f', [0.])
  recoMuMidPurity = array('f', [0.])
  recoMuLatePurity = array('f', [0.])
  recoMuNoVtxPurity = array('f', [0.])
  recoMuNoVtxAnscPurity = array('f', [0.])
nuIntTree.Branch("trueNuE", trueNuE, 'trueNuE/F')
nuIntTree.Branch("trueLepE", trueLepE, 'trueLepE/F')
nuIntTree.Branch("trueLepX", trueLepX, 'trueLepX/F')
nuIntTree.Branch("trueLepY", trueLepY, 'trueLepY/F')
nuIntTree.Branch("trueLepZ", trueLepZ, 'trueLepZ/F')
nuIntTree.Branch("trueLepPx", trueLepPx, 'trueLepPx/F')
nuIntTree.Branch("trueLepPy", trueLepPy, 'trueLepPy/F')
nuIntTree.Branch("trueLepPz", trueLepPz, 'trueLepPz/F')
nuIntTree.Branch("trueLepThX", trueLepThX, 'trueLepThX/F')
nuIntTree.Branch("trueLepThY", trueLepThY, 'trueLepThY/F')
nuIntTree.Branch("trueLepThZ", trueLepThZ, 'trueLepThZ/F')
nuIntTree.Branch("trueLepThU", trueLepThU, 'trueLepThU/F')
nuIntTree.Branch("trueLepThV", trueLepThV, 'trueLepThV/F')
nuIntTree.Branch("trueLepMinWireAng", trueLepMinWireAng, 'trueLepMinWireAng/F')
nuIntTree.Branch("trueLepMinSpecAng", trueLepMinSpecAng, 'trueLepMinSpecAng/F')
nuIntTree.Branch("recoStatus", recoStatus, 'recoStatus/I')
nuIntTree.Branch("closestVtxDist", closestVtxDist, 'closestVtxDist/F')
nuIntTree.Branch("vtxType", vtxType, 'vtxType/I')
nuIntTree.Branch("vtxDist", vtxDist, 'vtxDist/F')
nuIntTree.Branch("vtxDist_x", vtxDist_x, 'vtxDist_x/F')
nuIntTree.Branch("vtxDist_y", vtxDist_y, 'vtxDist_y/F')
nuIntTree.Branch("vtxDist_z", vtxDist_z, 'vtxDist_z/F')
nuIntTree.Branch("nCloseVertices05cm", nCloseVertices05cm, 'nCloseVertices05cm/I')
nuIntTree.Branch("nCloseVertices10cm", nCloseVertices10cm, 'nCloseVertices10cm/I')
nuIntTree.Branch("nCloseVertices15cm", nCloseVertices15cm, 'nCloseVertices15cm/I')
nuIntTree.Branch("trueLepCompletenessN", trueLepCompletenessN, 'trueLepCompletenessN/F')
nuIntTree.Branch("trueLepCompletenessI", trueLepCompletenessI, 'trueLepCompletenessI/F')
nuIntTree.Branch("recoLep1PixPurity", recoLep1PixPurity, 'recoLep1PixPurity/F')
nuIntTree.Branch("recoLep2PixPurity", recoLep2PixPurity, 'recoLep2PixPurity/F')
nuIntTree.Branch("recoLep3PixPurity", recoLep3PixPurity, 'recoLep3PixPurity/F')
nuIntTree.Branch("recoLepPixIPurity", recoLepPixIPurity, 'recoLepPixIPurity/F')
nuIntTree.Branch("recoLepPixIAnscPurity", recoLepPixIAnscPurity, 'recoLepPixIAnscPurity/F')
nuIntTree.Branch("recoLepStartPosErr", recoLepStartPosErr, 'recoLepStartPosErr/F')
nuIntTree.Branch("recoLepStartDirErr", recoLepStartDirErr, 'recoLepStartDirErr/F')
if args.numu:
  nuIntTree.Branch("trueMuContained", trueMuContained, 'trueMuContained/I')
  nuIntTree.Branch("recoMuEndPosErr", recoMuEndPosErr, 'recoMuEndPosErr/F')
  nuIntTree.Branch("recoMuTrackLength", recoMuTrackLength, 'recoMuTrackLength/F')
  nuIntTree.Branch("recoMuEarlyPurity", recoMuEarlyPurity, 'recoMuEarlyPurity/F')
  nuIntTree.Branch("recoMuMidPurity", recoMuMidPurity, 'recoMuMidPurity/F')
  nuIntTree.Branch("recoMuLatePurity", recoMuLatePurity, 'recoMuLatePurity/F')
  nuIntTree.Branch("recoMuNoVtxPurity", recoMuNoVtxPurity, 'recoMuNoVtxPurity/F')
  nuIntTree.Branch("recoMuNoVtxAnscPurity", recoMuNoVtxAnscPurity, 'recoMuNoVtxAnscPurity/F')


sce = larutil.SpaceChargeMicroBooNE()
mcNuVertexer = ublarcvapp.mctools.NeutrinoVertex()
truthTrackSCE = ublarcvapp.mctools.TruthTrackSCE()
truthShowerTrunkSCE = ublarcvapp.mctools.TruthShowerTrunkSCE()


evtCounter = 0

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
    lepPDG = 11 + 2*args.numu

    if nuInt.CCNC() != 0 or lep.PdgCode() != lepPDG or not isFiducial(trueVtxPos):
      continue

    if args.numu:
      mcleptons = ioll.get_data(larlite.data.kMCTrack, "mcreco")
    if args.nue:
      mcleptons = ioll.get_data(larlite.data.kMCShower, "mcreco")
    for mclepton in mcleptons:
      if mclepton.PdgCode() == lepPDG and mclepton.Process() == 'primary':
        mcLeptonUnCorr = mclepton
        break
  
    if not MCLeptonOkay(lep, mcLeptonUnCorr):
      print("Couldn't find MC lepton match!!!")
      continue

    if args.numu:
      trueMuContained[0] = isInDetector(mcLeptonUnCorr.End())
      if args.contained and not trueMuContained[0]:
        continue

    if args.nue and not isFiducial(mcLeptonUnCorr.End()):
      continue

    trueNuE[0] = nuInt.Nu().Momentum().E()
    trueLepE[0] = lep.Momentum().E()
    trueLepX[0] = lep.Position().X()
    trueLepY[0] = lep.Position().Y()
    trueLepZ[0] = lep.Position().Z()
    trueLepPx[0] = lep.Momentum().Px()
    trueLepPy[0] = lep.Momentum().Py()
    trueLepPz[0] = lep.Momentum().Pz()
    trueLepAxisAngles = getAxisAngles(lep.Momentum())
    wireAngles = [ trueLepAxisAngles[1], trueLepAxisAngles[3], trueLepAxisAngles[4] ]
    specAngles = [ trueLepAxisAngles[0], trueLepAxisAngles[1], trueLepAxisAngles[3], trueLepAxisAngles[4] ]
    trueLepThX[0] = trueLepAxisAngles[0]*(180./pi)
    trueLepThY[0] = trueLepAxisAngles[1]*(180./pi)
    trueLepThZ[0] = trueLepAxisAngles[2]*(180./pi)
    trueLepThU[0] = trueLepAxisAngles[3]*(180./pi)
    trueLepThV[0] = trueLepAxisAngles[4]*(180./pi)
    trueLepMinWireAng[0] = getMinAngle(wireAngles)*(180./pi)
    trueLepMinSpecAng[0] = getMinAngle(specAngles)*(180./pi)

    if args.numu:
      mcLepton = truthTrackSCE.applySCE(mcLeptonUnCorr)
      trueEndPos = mcLepton.End()
      if not trueMuContained[0]:
        trueEndPos = getLastPoint(mcLeptonUnCorr, mcLepton)
    if args.nue:
      mcLepton = truthShowerTrunkSCE.applySCE(mcLeptonUnCorr)
    trueStartPos = mcLepton.Vertex()
    trueStartDir = mcLepton.VertexDirection()

    closeVertices = []
    nCloseVertices = [0, 0, 0]
    bestVtxDist = 1e9

    for vtx in kpst.nufitted_v:
      deltaVertex = getVertexDistance(trueVtxPos, vtx)
      if deltaVertex < args.vertexRange:
        closeVertices.append(vtx)
      if deltaVertex < 15.:
        nCloseVertices[2] = nCloseVertices[2] + 1
        if deltaVertex < 10.:
          nCloseVertices[1] = nCloseVertices[1] + 1
          if deltaVertex < 5.:
            nCloseVertices[0] = nCloseVertices[0] + 1
      if deltaVertex < bestVtxDist:
        #bestVertex = vtx
        bestVtxDist = deltaVertex

    nCloseVertices05cm[0] = nCloseVertices[0]
    nCloseVertices10cm[0] = nCloseVertices[1]
    nCloseVertices15cm[0] = nCloseVertices[2]

    if kpst.nufitted_v.size() > 0:
      closestVtxDist[0] = bestVtxDist
    else:
      closestVtxDist[0] = -9.

    #if kpst.nufitted_v.size() < 1:
    if len(closeVertices) < 1:
      recoStatus[0] = 2
      if kpst.nufitted_v.size() < 1:
        recoStatus[0] = 3
      vtxType[0] = -1
      vtxDist[0] = -9.
      vtxDist_x[0] = -999.
      vtxDist_y[0] = -999.
      vtxDist_z[0] = -999.
      trueLepCompletenessN[0] = -1.
      trueLepCompletenessI[0] = -1.
      recoLep1PixPurity[0] = -1.
      recoLep2PixPurity[0] = -1.
      recoLep3PixPurity[0] = -1.
      recoLepPixIPurity[0] = -1.
      recoLepPixIAnscPurity[0] = -1.
      recoLepStartPosErr[0] = -9.
      recoLepStartDirErr[0] = -9.
      if args.numu:
        recoMuEndPosErr[0] = -9.
        recoMuTrackLength[0] = -9.
        recoMuEarlyPurity[0] = -1.
        recoMuMidPurity[0] = -1.
        recoMuLatePurity[0] = -1.
        recoMuNoVtxPurity[0] = -1.
        recoMuNoVtxAnscPurity[0] = -1.
      nuIntTree.Fill()
      continue
    #else:
      #vtxType[0] = bestVertex.keypoint_type
      #vtxDist[0] = bestVtxDist
      #vtxDist_x[0] = bestVertex.pos[0] - trueVtxPos.X()
      #vtxDist_y[0] = bestVertex.pos[1] - trueVtxPos.Y()
      #vtxDist_z[0] = bestVertex.pos[2] - trueVtxPos.Z()

    mcpg = ublarcvapp.mctools.MCPixelPGraph()
    mcpg.set_adc_treename("wire")
    mcpg.buildgraph(iolcv, ioll)
    lepNode = ublarcvapp.mctools.MCPixelPGraph.Node_t()
    lepNodeTID = -1
    foundLepNode = False
    for node in mcpg.node_v:
      if node.pid == lepPDG and node.tid == node.mtid:
        lepNode = node
        lepNodeTID = node.tid
        foundLepNode = True
        break
    if not foundLepNode:
      sys.exit("Couldn't find lepton node in MCPixelPGraph!!")
    lepAndDescNodes = mcpg.getNodeAndDescendentsFromTrackID(lepNodeTID)

    evtImage2D = iolcv.get_data(larcv.kProductImage2D, "wire")
    image2Dvec = evtImage2D.Image2DArray()

    nLepPixels = 0
    totLepPixI = 0.
    lepTickLists = [ [], [], [] ]
    lepPixelDictList = [ {}, {}, {} ]
    for p in range(3):
      lepPix = lepNode.pix_vv[p]
      nLepPixels = nLepPixels + lepPix.size()//2
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

    #nLepADPixels = 0
    #totLepADPixI = 0.
    lepADTickLists = [ [], [], [] ]
    lepADPixelDictList = [ {}, {}, {} ]
    for lepADNode in lepAndDescNodes:
      for p in range(3):
        lepADPix = lepADNode.pix_vv[p]
        #nLepADPixels = nLepADPixels + lepADPix.size()//2
        for iP in range(lepADPix.size()//2):
          simTick = lepADPix[2*iP]
          simWire = lepADPix[2*iP+1]
          #row = (simTick - 2400)//6
          #totLepADPixI = totLepADPixI + image2Dvec[p].pixel(row,simWire)
          if simTick not in lepADTickLists[p]:
            lepADTickLists[p].append(simTick)
            lepADPixelDictList[p][simTick] = [simWire]
          else:
            lepADPixelDictList[p][simTick].append(simWire)

    iVMatch = -1
    iCMatch = -1
    completenessN = -1.
    completenessI = -1.
    purity = [-1., -1., -1.]
    purityI = -1.
    purityIAD = -1.
    purityINV = -1.
    purityINVAD = -1.
    purityI_sec = [-1., -1., -1.]

    #------- begin vertex loop ---------------------------------------------------
    #if True:
    #  vertex = bestVertex
    for iV in range(len(closeVertices)):
      vertex = closeVertices[iV]
      #------- begin track/shower loop --------------------------------------------
      if args.numu:
        clusterVec = vertex.track_hitcluster_v
      if args.nue:
        clusterVec = vertex.shower_v
      for iC in range(clusterVec.size()):
        leptonCluster = clusterVec[iC]
        nUniquePixelsMatched = 0
        matchedPixels = []
        nHits = leptonCluster.size()
        nMatchedHits = [0, 0, 0]
        matchedPixI = 0.
        lad_matchedPixI = 0.
        matchedUniquePixI = 0.
        totPixI = 0.
        if args.numu:
          matchedPixI_sec = [0., 0., 0.]
          totPixI_sec = [0., 0., 0.]
          matchedPixI_nv = 0.
          matchedPixI_nvad = 0.
          totPixI_nv = 0.
          trackStart = vertex.track_v[iC].Vertex()
          trackEnd = vertex.track_v[iC].End()
          if args.v1:
            trackStart = vertex.track_v[iC].End()
            trackEnd = vertex.track_v[iC].Vertex()
        #------- begin hit loop ---------------------------------------------------
        for hit in leptonCluster:
          hit_nPixelsMatched = 0
          if args.numu:
            hitToStartDist = getHitDistance(trackStart, hit)
            trackSecFlags = getTrackSecFlags(trackStart, trackEnd, hit)
          #------- begin plane loop ---------------------------------------------------
          for p in range(3):
            row = (hit.tick - 2400)//6
            pixVal = image2Dvec[p].pixel(row, hit.targetwire[p])
            totPixI = totPixI + pixVal
            if args.numu:
              if hitToStartDist > 20.:
                totPixI_nv = totPixI_nv + pixVal
              for iTS in range(3):
                if trackSecFlags[iTS]:
                  totPixI_sec[iTS] = totPixI_sec[iTS] + pixVal
            # check if pixel belongs to primary lepton
            if hit.tick in lepTickLists[p]:
              if hit.targetwire[p] in lepPixelDictList[p][hit.tick]:
                matchedPixel = [ p, hit.tick, hit.targetwire[p] ]
                hit_nPixelsMatched = hit_nPixelsMatched + 1
                matchedPixI = matchedPixI + pixVal
                if args.numu:
                  if hitToStartDist > 20.:
                    matchedPixI_nv = matchedPixI_nv + pixVal
                  for iTS in range(3):
                    if trackSecFlags[iTS]:
                      matchedPixI_sec[iTS] = matchedPixI_sec[iTS] + pixVal
                if matchedPixel not in matchedPixels:
                  nUniquePixelsMatched = nUniquePixelsMatched + 1
                  matchedUniquePixI = matchedUniquePixI + pixVal
                  matchedPixels.append(matchedPixel)
            # check if pixel belongs to primary lepton or descendents
            if hit.tick in lepADTickLists[p]:
              if hit.targetwire[p] in lepADPixelDictList[p][hit.tick]:
                lad_matchedPixI = lad_matchedPixI + pixVal
                if args.numu:
                  if hitToStartDist > 20.:
                    matchedPixI_nvad = matchedPixI_nvad + pixVal
          #------- end plane loop ---------------------------------------------------
          if hit_nPixelsMatched >= 1:
            nMatchedHits[0] = nMatchedHits[0] + 1
          if hit_nPixelsMatched >= 2:
            nMatchedHits[1] = nMatchedHits[1] + 1
          if hit_nPixelsMatched == 3:
            nMatchedHits[2] = nMatchedHits[2] + 1
        #------- end hit loop ---------------------------------------------------
        cN = (1.0*nUniquePixelsMatched)/nLepPixels
        cI = matchedUniquePixI/totLepPixI
        if cI > completenessI:
          completenessN = cN
          completenessI = cI
          purityI = matchedPixI/totPixI
          purityIAD = lad_matchedPixI/totPixI
          if totPixI_nv > 0.:
            purityINV = matchedPixI_nv/totPixI_nv
            purityINVAD = matchedPixI_nvad/totPixI_nv
          else:
            purityINV = -1.
            purityINVAD = -1.
          for iT3 in range(3):
            if args.numu:
              if totPixI_sec[iT3] > 0.:
                purityI_sec[iT3] = matchedPixI_sec[iT3]/totPixI_sec[iT3]
              else:
                purityI_sec[iT3] = -1.
            purity[iT3] = (1.0*nMatchedHits[iT3])/nHits
          iVMatch = iV
          iCMatch = iC
      #------- end track/shower loop --------------------------------------------
    #------- end vertex loop ---------------------------------------------------
    

    #if iCMatch < 0:
    if iVMatch < 0 or iCMatch < 0:
      recoStatus[0] = 1
      vtxType[0] = -1
      vtxDist[0] = -9.
      vtxDist_x[0] = -999.
      vtxDist_y[0] = -999.
      vtxDist_z[0] = -999.
      trueLepCompletenessN[0] = -1.
      trueLepCompletenessI[0] = -1.
      recoLep1PixPurity[0] = -1.
      recoLep2PixPurity[0] = -1.
      recoLep3PixPurity[0] = -1.
      recoLepPixIPurity[0] = -1.
      recoLepPixIAnscPurity[0] = -1.
      recoLepStartPosErr[0] = -9.
      recoLepStartDirErr[0] = -9.
      if args.numu:
        recoMuEndPosErr[0] = -9.
        recoMuTrackLength[0] = -9.
        recoMuEarlyPurity[0] = -1.
        recoMuMidPurity[0] = -1.
        recoMuLatePurity[0] = -1.
        recoMuNoVtxPurity[0] = -1.
        recoMuNoVtxAnscPurity[0] = -1.
      nuIntTree.Fill()
      continue

    bestVertex = closeVertices[iVMatch]
    if args.numu:
      recoLepTrack = bestVertex.track_v[iCMatch]
      recoMuEndPosErr[0] = getMuEndDistance(trueEndPos, recoLepTrack)
      recoMuTrackLength[0] = recoLepTrack.Length()
      recoMuEarlyPurity[0] = purityI_sec[0]
      recoMuMidPurity[0] = purityI_sec[1]
      recoMuLatePurity[0] = purityI_sec[2]
      recoMuNoVtxPurity[0] = purityINV
      recoMuNoVtxAnscPurity[0] = purityINVAD
    if args.nue:
      recoLepTrack = bestVertex.shower_trunk_v[iCMatch]
    recoStatus[0] = 0
    vtxType[0] = bestVertex.keypoint_type
    vtxDist[0] = bestVtxDist
    vtxDist_x[0] = bestVertex.pos[0] - trueVtxPos.X()
    vtxDist_y[0] = bestVertex.pos[1] - trueVtxPos.Y()
    vtxDist_z[0] = bestVertex.pos[2] - trueVtxPos.Z()
    trueLepCompletenessN[0] = completenessN
    trueLepCompletenessI[0] = completenessI
    recoLep1PixPurity[0] = purity[0]
    recoLep2PixPurity[0] = purity[1]
    recoLep3PixPurity[0] = purity[2]
    recoLepPixIPurity[0] = purityI
    recoLepPixIAnscPurity[0] = purityIAD
    recoLepStartPosErr[0] = getStartDistance(trueStartPos, recoLepTrack)
    recoLepStartDirErr[0] = getAngle(trueStartDir, recoLepTrack)*(180./pi)
    nuIntTree.Fill()

    if recoStatus[0] == 0 and completenessI < 0.05 and purityIAD < 0.05:
      print("BAD RECO EVENT: %s, %s, %i/%i/%i, %i, %f/%i"%(filepair[0], filepair[1], kpst.run, kpst.subrun, kpst.event, ientry, nuInt.Nu().Momentum().E(), nuInt.Nu().PdgCode()))
    if recoStatus[0] == 0 and completenessI > 0.95 and purityIAD > 0.95:
      print("GREAT RECO EVENT: %s, %s, %i/%i/%i, %i, %f/%i"%(filepair[0], filepair[1], kpst.run, kpst.subrun, kpst.event, ientry, nuInt.Nu().Momentum().E(), nuInt.Nu().PdgCode()))

    evtCounter = evtCounter + 1

  #++++++ end entry loop ++++++++++++++++++++++++++++++++++++++++++++++++++++=

  ioll.close()
  iolcv.finalize()
  kpsfile.Close()

#-------- end file loop -----------------------------------------------------#


outRootFile.cd()
nuIntTree.Write("",rt.TObject.kOverwrite)
outRootFile.Close()
print("analyzed %i CC events"%evtCounter)


