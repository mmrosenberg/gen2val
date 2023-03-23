
import os,sys,argparse

import ROOT as rt
import uproot

from larlite import larlite
from larlite import larutil
from ublarcvapp import ublarcvapp
from larcv import larcv
from larflow import larflow

from math import sqrt as sqrt
from math import acos as acos
from math import pi

from array import array
import numpy as np

from event_weighting.event_weight_helper import SumPOT, Weights
from helpers.larflowreco_ana_funcs import *

import torch
from torch import nn
from torch.utils.data import DataLoader
import torchvision.transforms as transforms

parser = argparse.ArgumentParser("Get Vertex and Prong Info Needed for Selections")
parser.add_argument("-f", "--files", required=True, type=str, nargs="+", help="input kpsreco files")
parser.add_argument("-t", "--truth", required=True, type=str, help="text file containing merged_dlreco list")
parser.add_argument("-w", "--weightfile", type=str, default="none", help="weights file (pickled python dict)")
parser.add_argument("-m", "--model_path", type=str, required=True, help="path to prong CNN checkpoint file")
parser.add_argument("-d", "--device", type=str, default="cpu", help="gpu/cpu device")
parser.add_argument("-mc", "--isMC", help="running over MC input", action="store_true")
parser.add_argument("-o", "--outfile", type=str, default="prepare_selection_test_output.root", help="output file name")
parser.add_argument("--makePlots", action="store_true", help="make image plots for specified event")
parser.add_argument("-r", "--run", type=int, default=0, help="run number for event to plot")
parser.add_argument("-sr", "--subrun", type=int, default=0, help="subrun number for event to plot")
parser.add_argument("-e", "--event", type=int, default=0, help="event number for event to plot")
parser.add_argument("--multiGPU", action="store_true", help="use multiple GPUs")
parser.add_argument("--oldVtxBranch", help="use nufitted_v instead of nuvetoed_v for old reco", action="store_true")
args = parser.parse_args()

sys.path.append(args.model_path[:args.model_path.find("/checkpoints")])
from models_instanceNorm_reco_2chan_multiTask import ResBlock, ResNet34
from datasets_reco_5ClassHardLabel_multiTask import mean, std

if args.isMC and args.makePlots:
  import matplotlib.pyplot as plt
  import matplotlib.backends.backend_pdf
  #plotfilename = "selection_test_images_run%i_subrun%i_event%i.pdf"%(args.run,args.subrun,args.event)
  plotfilename = "selection_test_images_for_score_removed_CCnumu_background.pdf"
  outpdf = matplotlib.backends.backend_pdf.PdfPages(plotfilename)

if args.isMC and args.weightfile=="none":
  sys.exit("Must supply weight file for MC input. Exiting...")

reco2Tag = "merged_dlana_"
if args.isMC:
  reco2Tag = "merged_dlreco_"
files = getFiles(reco2Tag, args.files, args.truth)

plotEvents = [[16960, 184, 9249], [16962, 75, 3769], [16935, 300, 15039], [16944, 110, 5549], [16947, 210, 10533], [18508, 291, 14580], [18527, 20, 1034], [18511, 29, 1461], [18514, 195, 9777], [18577, 320, 16019], [18525, 115, 5767], [18520, 13, 676], [18005, 67, 3367], [18021, 3, 166], [18031, 383, 19184], [18032, 515, 25798], [18860, 124, 6248], [18863, 152, 7650], [18874, 72, 3648], [18862, 415, 20765]]


def addClusterCharge(iolcv, cluster, vertexPixels, vertexCharge, threshold):
  evtImage2D = iolcv.get_data(larcv.kProductImage2D, "wire")
  image2Dvec = evtImage2D.Image2DArray()
  clusterPixels = []
  clusterCharge = 0.
  for hit in cluster:
    for p in range(3):
      row = (hit.tick - 2400)//6
      pixel = [ p, hit.tick, hit.targetwire[p] ]
      pixVal = image2Dvec[p].pixel(row, hit.targetwire[p])
      if pixVal < threshold:
        continue
      if pixel not in clusterPixels:
        clusterPixels.append(pixel)
        clusterCharge += pixVal
      if pixel not in vertexPixels:
        vertexPixels.append(pixel)
        vertexCharge += pixVal
  return clusterCharge, vertexPixels, vertexCharge



def getMCProngParticle(sparseimg_vv, mcpg, mcpm, adc_v):

  particleDict = {}
  trackDict = {}
  totalPixI = 0.

  for p in range(3):
    for pix in sparseimg_vv[p]:
      totalPixI += pix.val
      pixContents = mcpm.getPixContent(p, pix.rawRow, pix.rawCol)
      for part in pixContents.particles:
        if abs(part.pdg) in particleDict:
          particleDict[abs(part.pdg)] += pixContents.pixI
        else:
          particleDict[abs(part.pdg)] = pixContents.pixI
        if part.tid in trackDict:
          trackDict[part.tid][2] += pixContents.pixI
        else:
          trackDict[part.tid] = [part.pdg, part.nodeidx, pixContents.pixI]

  maxPartPDG = 0 
  maxPartNID = -1
  maxPartTID = -1
  maxPartI = 0.
  maxPartComp = 0.
  pdglist = []
  puritylist = []

  for part in particleDict:
    pdglist.append(part)
    puritylist.append(particleDict[part]/totalPixI)

  for track in trackDict:
    if trackDict[track][2] > maxPartI:
      maxPartI = trackDict[track][2]
      maxPartPDG = trackDict[track][0]
      maxPartNID = trackDict[track][1]
      maxPartTID = track

  totNodePixI = 0.
  if maxPartI > 0.:
    maxPartNode = mcpg.node_v[maxPartNID]
    if maxPartNode.tid != maxPartTID:
      sys.exit("ERROR: mismatch between node track id from mcpm and mcpg in getMCProngParticle")
    for p in range(3):
      pixels = maxPartNode.pix_vv[p]
      for iP in range(pixels.size()//2):
        row = (pixels[2*iP] - 2400)//6
        col = pixels[2*iP+1]
        totNodePixI += adc_v[p].pixel(row, col)
    if totNodePixI > 0.:
      maxPartComp = maxPartI/totNodePixI

  if maxPartComp > 1.:
    print("ERROR: prong completeness calculated to be >1")

  #return maxPartPDG, maxPartTID, totNodePixI, maxPartI/totalPixI, maxPartComp, pdglist, puritylist
  return maxPartPDG, maxPartI/totalPixI, maxPartComp, pdglist, puritylist


def makeImage(prong_vv):
  plane0pix_row = np.zeros(prong_vv[0].size(), dtype=int)
  plane0pix_col = np.zeros(prong_vv[0].size(), dtype=int)
  plane0pix_val = np.zeros(prong_vv[0].size(), dtype=float)
  plane1pix_row = np.zeros(prong_vv[1].size(), dtype=int)
  plane1pix_col = np.zeros(prong_vv[1].size(), dtype=int)
  plane1pix_val = np.zeros(prong_vv[1].size(), dtype=float)
  plane2pix_row = np.zeros(prong_vv[2].size(), dtype=int)
  plane2pix_col = np.zeros(prong_vv[2].size(), dtype=int)
  plane2pix_val = np.zeros(prong_vv[2].size(), dtype=float)
  raw_plane0pix_row = np.zeros(prong_vv[3].size(), dtype=int)
  raw_plane0pix_col = np.zeros(prong_vv[3].size(), dtype=int)
  raw_plane0pix_val = np.zeros(prong_vv[3].size(), dtype=float)
  raw_plane1pix_row = np.zeros(prong_vv[4].size(), dtype=int)
  raw_plane1pix_col = np.zeros(prong_vv[4].size(), dtype=int)
  raw_plane1pix_val = np.zeros(prong_vv[4].size(), dtype=float)
  raw_plane2pix_row = np.zeros(prong_vv[5].size(), dtype=int)
  raw_plane2pix_col = np.zeros(prong_vv[5].size(), dtype=int)
  raw_plane2pix_val = np.zeros(prong_vv[5].size(), dtype=float)
  for i, pix in enumerate(prong_vv[0]):
    plane0pix_row[i] = pix.row
    plane0pix_col[i] = pix.col
    plane0pix_val[i] = pix.val
  for i, pix in enumerate(prong_vv[1]):
    plane1pix_row[i] = pix.row
    plane1pix_col[i] = pix.col
    plane1pix_val[i] = pix.val
  for i, pix in enumerate(prong_vv[2]):
    plane2pix_row[i] = pix.row
    plane2pix_col[i] = pix.col
    plane2pix_val[i] = pix.val
  for i, pix in enumerate(prong_vv[3]):
    raw_plane0pix_row[i] = pix.row
    raw_plane0pix_col[i] = pix.col
    raw_plane0pix_val[i] = pix.val
  for i, pix in enumerate(prong_vv[4]):
    raw_plane1pix_row[i] = pix.row
    raw_plane1pix_col[i] = pix.col
    raw_plane1pix_val[i] = pix.val
  for i, pix in enumerate(prong_vv[5]):
    raw_plane2pix_row[i] = pix.row
    raw_plane2pix_col[i] = pix.col
    raw_plane2pix_val[i] = pix.val
  image = np.zeros((6,512,512))
  image[0, plane0pix_row, plane0pix_col] = plane0pix_val
  image[2, plane1pix_row, plane1pix_col] = plane1pix_val
  image[4, plane2pix_row, plane2pix_col] = plane2pix_val
  image[1, raw_plane0pix_row, raw_plane0pix_col] = raw_plane0pix_val
  image[3, raw_plane1pix_row, raw_plane1pix_col] = raw_plane1pix_val
  image[5, raw_plane2pix_row, raw_plane2pix_col] = raw_plane2pix_val
  image = torch.from_numpy(image).float()
  norm = transforms.Normalize(mean, std)
  image = norm(image).reshape(1,6,512,512)
  return torch.clamp(image, max=4.0)


def plotEventTitle(r, sr, e, primInfo):
  fig = plt.figure(0, clear=True)
  plt.suptitle("\n\n\nRun %i Subrun %i Event %i \n\n\n %s"%(r, sr, e, primInfo), fontsize=12)
  outpdf.savefig(fig)


def plotImage(X, r, sr, e, pdg, purity, completeness, visE, compPred, elScore, phScore, muScore, piScore, prScore):
  #pltmin = None
  #pltmax = None
  pltmin = -1.
  pltmax = 2.
  suptitlesize=6
  titlesize=8
  ticksize=6
  X0 = X[:,0:2].reshape(X.shape[0], 2, 512, 512)
  X1 = X[:,2:4].reshape(X.shape[0], 2, 512, 512)
  X2 = X[:,4:].reshape(X.shape[0], 2, 512, 512)
  fig = plt.figure(0, clear=True)
  plt.subplot(2,3,1)
  plt.imshow(X0.numpy()[0][0], vmin=pltmin, vmax=pltmax, cmap='jet')
  plt.title("plane 0 prong", fontsize=titlesize)
  plt.xticks(fontsize=ticksize)
  plt.yticks(fontsize=ticksize)
  plt.subplot(2,3,2)
  plt.imshow(X1.numpy()[0][0], vmin=pltmin, vmax=pltmax, cmap='jet')
  plt.title("plane 1 prong", fontsize=titlesize)
  plt.xticks(fontsize=ticksize)
  plt.yticks(fontsize=ticksize)
  plt.subplot(2,3,3)
  plt.imshow(X2.numpy()[0][0], vmin=pltmin, vmax=pltmax, cmap='jet')
  plt.title("plane 2 prong", fontsize=titlesize)
  plt.xticks(fontsize=ticksize)
  plt.yticks(fontsize=ticksize)
  plt.subplot(2,3,4)
  plt.imshow(X0.numpy()[0][1], vmin=pltmin, vmax=pltmax, cmap='jet')
  plt.title("plane 0 all", fontsize=titlesize)
  plt.xticks(fontsize=ticksize)
  plt.yticks(fontsize=ticksize)
  plt.subplot(2,3,5)
  plt.imshow(X1.numpy()[0][1], vmin=pltmin, vmax=pltmax, cmap='jet')
  plt.title("plane 1 all", fontsize=titlesize)
  plt.xticks(fontsize=ticksize)
  plt.yticks(fontsize=ticksize)
  plt.subplot(2,3,6)
  plt.imshow(X2.numpy()[0][1], vmin=pltmin, vmax=pltmax, cmap='jet')
  plt.title("plane 2 all", fontsize=titlesize)
  plt.xticks(fontsize=ticksize)
  plt.yticks(fontsize=ticksize)
  plt.suptitle("Run %i Subrun %i Event %i  |  Prong: pdg %i, purity %.2f, completeness %.2f, visible energy %.2e \n e- score %.2f, photon score %.2f, mu score: %.2f, pi score %.2f, proton score %.2f, completeness prediction: %.2f"%(r, sr, e, pdg, purity, completeness, visE, elScore, phScore, muScore, piScore, prScore, compPred), fontsize=suptitlesize)
  #plt.show()
  outpdf.savefig(fig)
  return


def getPID(cnnClass):
    if cnnClass == 0:
        return 11
    if cnnClass == 1:
        return 22
    if cnnClass == 2:
        return 13
    if cnnClass == 3:
        return 211
    if cnnClass == 4:
        return 2212
    return 0


#other needed classes
sce = larutil.SpaceChargeMicroBooNE()
mcNuVertexer = ublarcvapp.mctools.NeutrinoVertex()
prongvars = larflow.reco.NuSelProngVars()
wcoverlapvars = larflow.reco.NuSelWCTaggerOverlap()
flowTriples = larflow.prep.FlowTriples()

model = ResNet34(2, ResBlock, outputs=5)
if "cuda" in args.device and args.multiGPU:
  model = nn.DataParallel(model)
if args.device == "cpu":
  checkpoint = torch.load(args.model_path, map_location=torch.device('cpu'))
else:
  checkpoint = torch.load(args.model_path)
try:
  model.load_state_dict(checkpoint['model_state_dict'])
except:
  model.module.load_state_dict(checkpoint['model_state_dict'])
model.to(args.device)
model.eval()

outRootFile = rt.TFile(args.outfile, "RECREATE")

if args.isMC:
  potTree = rt.TTree("potTree","potTree")
  totPOT = array('f', [0.])
  totGoodPOT = array('f', [0.])
  potTree.Branch("totPOT", totPOT, 'totPOT/F')
  potTree.Branch("totGoodPOT", totGoodPOT, 'totGoodPOT/F')

eventTree = rt.TTree("EventTree","EventTree")
maxNTrks = 100
maxNShwrs = 100
fileid = array('i', [0])
run = array('i', [0])
subrun = array('i', [0])
event = array('i', [0])
xsecWeight = array('f', [0.])
trueNuE = array('f', [0.])
trueLepE = array('f', [0.])
trueLepPDG = array('i', [0])
trueMuContained = array('i', [0])
trueNuPDG = array('i', [0])
trueNuCCNC = array('i', [0])
nVertices = array('i', [0])
vtxX = array('f', [0.])
vtxY = array('f', [0.])
vtxZ = array('f', [0.])
vtxIsFiducial = array('i', [0])
vtxDistToTrue = array('f', [0.])
vtxBestComp = array('f', [0.])
vtxScore = array('f', [0.])
vtxFracHitsOnCosmic = array('f', [0.])
nTracks = array('i', [0])
trackIsSecondary = array('i', maxNTrks*[0])
trackNHits = array('i', maxNTrks*[0])
trackHitFrac = array('f', maxNTrks*[0.])
trackCharge = array('f', maxNTrks*[0.])
trackChargeFrac = array('f', maxNTrks*[0.])
trackMuKE = array('f', maxNTrks*[0.])
trackPrKE = array('f', maxNTrks*[0.])
trackCosTheta = array('f', maxNTrks*[0.])
trackDistToVtx = array('f', maxNTrks*[0.])
trackClassified = array('i', maxNTrks*[0])
trackPID = array('i', maxNTrks*[0])
trackElScore = array('f', maxNTrks*[0.])
trackPhScore = array('f', maxNTrks*[0.])
trackMuScore = array('f', maxNTrks*[0.])
trackPiScore = array('f', maxNTrks*[0.])
trackPrScore = array('f', maxNTrks*[0.])
trackComp = array('f', maxNTrks*[0.])
trackTruePID = array('i', maxNTrks*[0])
trackTruePurity = array('f', maxNTrks*[0.])
trackTrueComp = array('f', maxNTrks*[0.])
trackTrueElPurity = array('f', maxNTrks*[0.])
trackTruePhPurity = array('f', maxNTrks*[0.])
trackTrueMuPurity = array('f', maxNTrks*[0.])
trackTruePiPurity = array('f', maxNTrks*[0.])
trackTruePrPurity = array('f', maxNTrks*[0.])
nShowers = array('i', [0])
showerIsSecondary = array('i', maxNShwrs*[0])
showerNHits = array('i', maxNShwrs*[0])
showerHitFrac = array('f', maxNShwrs*[0.])
showerCharge = array('f', maxNShwrs*[0.])
showerChargeFrac = array('f', maxNShwrs*[0.])
showerPl0E = array('f', maxNShwrs*[0.])
showerPl1E = array('f', maxNShwrs*[0.])
showerPl2E = array('f', maxNShwrs*[0.])
showerCosTheta = array('f', maxNShwrs*[0.])
showerDistToVtx = array('f', maxNShwrs*[0.])
showerClassified = array('i', maxNShwrs*[0])
showerPID = array('i', maxNShwrs*[0])
showerElScore = array('f', maxNShwrs*[0.])
showerPhScore = array('f', maxNShwrs*[0.])
showerMuScore = array('f', maxNShwrs*[0.])
showerPiScore = array('f', maxNShwrs*[0.])
showerPrScore = array('f', maxNShwrs*[0.])
showerComp = array('f', maxNShwrs*[0])
showerTruePID = array('i', maxNShwrs*[0])
showerTruePurity = array('f', maxNShwrs*[0.])
showerTrueComp = array('f', maxNShwrs*[0.])
showerTrueElPurity = array('f', maxNShwrs*[0.])
showerTruePhPurity = array('f', maxNShwrs*[0.])
showerTrueMuPurity = array('f', maxNShwrs*[0.])
showerTruePiPurity = array('f', maxNShwrs*[0.])
showerTruePrPurity = array('f', maxNShwrs*[0.])
eventTree.Branch("fileid", fileid, 'fileid/I')
eventTree.Branch("run", run, 'run/I')
eventTree.Branch("subrun", subrun, 'subrun/I')
eventTree.Branch("event", event, 'event/I')
eventTree.Branch("xsecWeight", xsecWeight, 'xsecWeight/F')
eventTree.Branch("trueNuE", trueNuE, 'trueNuE/F')
eventTree.Branch("trueLepE", trueLepE, 'trueLepE/F')
eventTree.Branch("trueLepPDG", trueLepPDG, 'trueLepPDG/I')
eventTree.Branch("trueMuContained", trueMuContained, 'trueMuContained/I')
eventTree.Branch("trueNuPDG", trueNuPDG, 'trueNuPDG/I')
eventTree.Branch("trueNuCCNC", trueNuCCNC, 'trueNuCCNC/I')
eventTree.Branch("nVertices", nVertices, 'nVertices/I')
eventTree.Branch("vtxX", vtxX, 'vtxX/F')
eventTree.Branch("vtxY", vtxY, 'vtxY/F')
eventTree.Branch("vtxZ", vtxZ, 'vtxZ/F')
eventTree.Branch("vtxIsFiducial", vtxIsFiducial, 'vtxIsFiducial/I')
eventTree.Branch("vtxDistToTrue", vtxDistToTrue, 'vtxDistToTrue/F')
eventTree.Branch("vtxBestComp", vtxBestComp, 'vtxBestComp/F')
eventTree.Branch("vtxScore", vtxScore, 'vtxScore/F')
eventTree.Branch("vtxFracHitsOnCosmic", vtxFracHitsOnCosmic, 'vtxFracHitsOnCosmic/F')
eventTree.Branch("nTracks", nTracks, 'nTracks/I')
eventTree.Branch("trackIsSecondary", trackIsSecondary, 'trackIsSecondary[nTracks]/I')
eventTree.Branch("trackNHits", trackNHits, 'trackNHits[nTracks]/I')
eventTree.Branch("trackHitFrac", trackHitFrac, 'trackHitFrac[nTracks]/F')
eventTree.Branch("trackCharge", trackCharge, 'trackCharge[nTracks]/F')
eventTree.Branch("trackChargeFrac", trackChargeFrac, 'trackChargeFrac[nTracks]/F')
eventTree.Branch("trackMuKE", trackMuKE, 'trackMuKE[nTracks]/F')
eventTree.Branch("trackPrKE", trackPrKE, 'trackPrKE[nTracks]/F')
eventTree.Branch("trackCosTheta", trackCosTheta, 'trackCosTheta[nTracks]/F')
eventTree.Branch("trackDistToVtx", trackDistToVtx, 'trackDistToVtx[nTracks]/F')
eventTree.Branch("trackClassified", trackClassified, 'trackClassified[nTracks]/I')
eventTree.Branch("trackPID", trackPID, 'trackPID[nTracks]/I')
eventTree.Branch("trackElScore", trackElScore, 'trackElScore[nTracks]/F')
eventTree.Branch("trackPhScore", trackPhScore, 'trackPhScore[nTracks]/F')
eventTree.Branch("trackMuScore", trackMuScore, 'trackMuScore[nTracks]/F')
eventTree.Branch("trackPiScore", trackPiScore, 'trackPiScore[nTracks]/F')
eventTree.Branch("trackPrScore", trackPrScore, 'trackPrScore[nTracks]/F')
eventTree.Branch("trackComp", trackComp, 'trackComp[nTracks]/F')
eventTree.Branch("trackTruePID", trackTruePID, 'trackTruePID[nTracks]/I')
eventTree.Branch("trackTruePurity", trackTruePurity, 'trackTruePurity[nTracks]/F')
eventTree.Branch("trackTrueComp", trackTrueComp, 'trackTrueComp[nTracks]/F')
eventTree.Branch("trackTrueElPurity", trackTrueElPurity, 'trackTrueElPurity[nTracks]/F')
eventTree.Branch("trackTruePhPurity", trackTruePhPurity, 'trackTruePhPurity[nTracks]/F')
eventTree.Branch("trackTrueMuPurity", trackTrueMuPurity, 'trackTrueMuPurity[nTracks]/F')
eventTree.Branch("trackTruePiPurity", trackTruePiPurity, 'trackTruePiPurity[nTracks]/F')
eventTree.Branch("trackTruePrPurity", trackTruePrPurity, 'trackTruePrPurity[nTracks]/F')
eventTree.Branch("nShowers", nShowers, 'nShowers/I')
eventTree.Branch("showerIsSecondary", showerIsSecondary, 'showerIsSecondary[nShowers]/I')
eventTree.Branch("showerNHits", showerNHits, 'showerNHits[nShowers]/I')
eventTree.Branch("showerHitFrac", showerHitFrac, 'showerHitFrac[nShowers]/F')
eventTree.Branch("showerCharge", showerCharge, 'showerCharge[nShowers]/F')
eventTree.Branch("showerChargeFrac", showerChargeFrac, 'showerChargeFrac[nShowers]/F')
eventTree.Branch("showerPl0E", showerPl0E, 'showerPl0E[nShowers]/F')
eventTree.Branch("showerPl1E", showerPl1E, 'showerPl1E[nShowers]/F')
eventTree.Branch("showerPl2E", showerPl2E, 'showerPl2E[nShowers]/F')
eventTree.Branch("showerCosTheta", showerCosTheta, 'showerCosTheta[nShowers]/F')
eventTree.Branch("showerDistToVtx", showerDistToVtx, 'showerDistToVtx[nShowers]/F')
eventTree.Branch("showerClassified", showerClassified, 'showerClassified[nShowers]/I')
eventTree.Branch("showerPID", showerPID, 'showerPID[nShowers]/I')
eventTree.Branch("showerTruePID", showerTruePID, 'showerTruePID[nShowers]/I')
eventTree.Branch("showerTruePurity", showerTruePurity, 'showerTruePurity[nShowers]/F')
eventTree.Branch("showerTrueComp", showerTrueComp, 'showerTrueComp[nShowers]/F')
eventTree.Branch("showerTrueElPurity", showerTrueElPurity, 'showerTrueElPurity[nShowers]/F')
eventTree.Branch("showerTruePhPurity", showerTruePhPurity, 'showerTruePhPurity[nShowers]/F')
eventTree.Branch("showerTrueMuPurity", showerTrueMuPurity, 'showerTrueMuPurity[nShowers]/F')
eventTree.Branch("showerTruePiPurity", showerTruePiPurity, 'showerTruePiPurity[nShowers]/F')
eventTree.Branch("showerTruePrPurity", showerTruePrPurity, 'showerTruePrPurity[nShowers]/F')
eventTree.Branch("showerElScore", showerElScore, 'showerElScore[nShowers]/F')
eventTree.Branch("showerPhScore", showerPhScore, 'showerPhScore[nShowers]/F')
eventTree.Branch("showerMuScore", showerMuScore, 'showerMuScore[nShowers]/F')
eventTree.Branch("showerPiScore", showerPiScore, 'showerPiScore[nShowers]/F')
eventTree.Branch("showerPrScore", showerPrScore, 'showerPrScore[nShowers]/F')
eventTree.Branch("showerComp", showerComp, 'showerComp[nShowers]/F')


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

    #print("reached entry:", ientry)

    ioll.go_to(ientry)
    iolcv.read_entry(ientry)
    kpst.GetEntry(ientry)
  
    if kpst.run != ioll.run_id() or kpst.subrun != ioll.subrun_id() or kpst.event != ioll.event_id():
      print("EVENTS DON'T MATCH!!!")
      print("truth run/subrun/event: %i/%i/%i"%(ioll.run_id(),ioll.subrun_id(),ioll.event_id()))
      print("reco run/subrun/event: %i/%i/%i"%(kpst.run,kpst.subrun,kpst.event))
      continue


    if args.isMC:

      if args.makePlots and [kpst.run, kpst.subrun, kpst.event] not in plotEvents:
        continue

      mctruth = ioll.get_data(larlite.data.kMCTruth, "generator")
      nuInt = mctruth.at(0).GetNeutrino()
      lep = nuInt.Lepton()
      mcNuVertex = mcNuVertexer.getPos3DwSCE(ioll, sce)
      trueVtxPos = rt.TVector3(mcNuVertex[0], mcNuVertex[1], mcNuVertex[2])

      if not isFiducial(trueVtxPos):
        continue

      try:
        xsecWeight[0] = weights.get(kpst.run, kpst.subrun, kpst.event)
      except:
        print("Couldn't find weight for run %i, subrun %i, event %i in %s!!!"%(kpst.run, kpst.subrun, kpst.event, args.weightfile))
        continue

      if nuInt.CCNC() == 0:

        if lep.PdgCode() not in [11,13]:
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

        totLepPixI, lepTickLists, lepPixelDictList = getLeptonPixels(lepPDG, ioll, iolcv)
        if not totLepPixI > 0.:
          continue

        trueMuContained[0] = -1
        if lepPDG == 13:
          if isInDetector(mcLeptonUnCorr.End()):
            trueMuContained[0] = 1
          else:
            trueMuContained[0] = 0

        trueLepPDG[0] = lepPDG
        trueLepE[0] = lep.Momentum().E()

      else: #from "if nuInt.CCNC"
        trueLepPDG[0] = 0
        trueLepE[0] = -9.

      trueNuPDG[0] = nuInt.Nu().PdgCode()
      trueNuCCNC[0] = nuInt.CCNC()
      trueNuE[0] = nuInt.Nu().Momentum().E()
      

    else: #from "if args.isMC"
      trueLepPDG[0] = 0
      trueMuContained[0] = -1
      trueLepE[0] = -9.
      trueNuPDG[0] = 0
      trueNuCCNC[0] = -1
      trueNuE[0] = -9.

    fileid[0] = int(filepair[0][filepair[0].find("fileid")+6:filepair[0].find("fileid")+10])
    run[0] = kpst.run
    subrun[0] = kpst.subrun
    event[0] = kpst.event

    vertices = kpst.nuvetoed_v
    if args.oldVtxBranch:
      vertices = kpst.nufitted_v
  
    nVertices[0] = 0
    vtxScore[0] = -1.
    foundVertex = False
    for vtx in vertices:
      if vtx.keypoint_type != 0:
        continue
      nVertices[0] += 1
      foundVertex = True
      if vtx.netNuScore > vtxScore[0]:
        vtxScore[0] = vtx.netNuScore
        vertex = vtx

    if not foundVertex:
      vtxX[0] = -9999.
      vtxY[0] = -9999.
      vtxZ[0] = -9999.
      vtxIsFiducial[0] = -1
      vtxDistToTrue[0] = -99.
      vtxBestComp[0] = -1.
      vtxFracHitsOnCosmic[0] = -1.
      nTracks[0] = 0
      nShowers[0] = 0
      eventTree.Fill()
      continue

    if args.isMC:
      vtxDistToTrue[0] = getVertexDistance(trueVtxPos, vertex)
      if nuInt.CCNC() == 0:
        vtxBestComp[0] = getBestCompleteness(iolcv, vertex, lepPDG, totLepPixI, lepTickLists, lepPixelDictList)
      else:
        vtxBestComp[0] = -1.
      mcpg = ublarcvapp.mctools.MCPixelPGraph()
      mcpg.set_adc_treename("wire")
      mcpg.buildgraph(iolcv, ioll)
      mcpm = ublarcvapp.mctools.MCPixelPMap()
      mcpm.set_adc_treename("wire")
      mcpm.buildmap(iolcv, mcpg)
      if args.makePlots:
        primaryPartInfo = "Primary Particles:\n\n"
        for node in mcpg.node_v:
          if (node.tid == node.mtid and node.origin == 1 and node.process == "primary"):
            primaryPartInfo += "  %i (%.2f MeV)\n"%(node.pid, node.E_MeV)
        #primaryPartInfo = primaryPartInfo[:-2]
        plotEventTitle(run[0], subrun[0], event[0], primaryPartInfo)
    else:
      vtxDistToTrue[0] = -99.
      vtxBestComp[0] = -1.

    vtxX[0] = vertex.pos[0]
    vtxY[0] = vertex.pos[1]
    vtxZ[0] = vertex.pos[2]
    vtxTVec3 = rt.TVector3(vertex.pos[0], vertex.pos[1], vertex.pos[2])
    vtxIsFiducial[0] = int(isFiducial(vtxTVec3))

    nusel = larflow.reco.NuSelectionVariables()
    wcoverlapvars.analyze(vertex, nusel, iolcv)
    vtxFracHitsOnCosmic[0] = nusel.frac_allhits_on_cosmic

    nTracks[0] = vertex.track_v.size()
    nShowers[0] = vertex.shower_v.size()

    evtImage2D = iolcv.get_data(larcv.kProductImage2D, "wire")
    csmImage2D = iolcv.get_data(larcv.kProductImage2D, "thrumu")
    adc_v = evtImage2D.Image2DArray()
    thrumu_v = csmImage2D.Image2DArray()

    vertexPixels = []
    vertexCharge = 0.
    vertexNHits = 0


    for iTrk, trackCls in enumerate(vertex.track_hitcluster_v):

      vertexNHits += trackCls.size()
      trackIsSecondary[iTrk] = vertex.track_isSecondary_v[iTrk]
      trackNHits[iTrk] = trackCls.size()
      trackCharge[iTrk], vertexPixels, vertexCharge = addClusterCharge(iolcv,trackCls,vertexPixels,vertexCharge,10.)
      trackMuKE[iTrk] = vertex.track_kemu_v[iTrk]
      trackPrKE[iTrk] = vertex.track_keproton_v[iTrk]
      nTrajPoints = vertex.track_v[iTrk].NumberTrajectoryPoints()
      trackLength = getDistance(vertex.track_v[iTrk].Vertex(),vertex.track_v[iTrk].End()) if (nTrajPoints > 1) else -9.
      goodTrack = nTrajPoints > 1 and trackLength > 1e-6
      trackCosTheta[iTrk] = getCosThetaBeamTrack(vertex.track_v[iTrk]) if goodTrack else -9.
      trackDistToVtx[iTrk] = getVertexDistance(vertex.track_v[iTrk].Vertex(), vertex) if goodTrack else -9.

      skip = True
      if goodTrack:
        skip = False
        cropPt = vertex.track_v[iTrk].End()
        prong_vv = flowTriples.make_cropped_initial_sparse_prong_image_reco(adc_v,thrumu_v,trackCls,cropPt,10.,512,512)
        for p in range(3):
          if prong_vv[p].size() < 10:
            skip = True
            break
      if skip:
        trackClassified[iTrk] = 0
        trackPID[iTrk] = 0
        trackElScore[iTrk] = -99.
        trackPhScore[iTrk] = -99.
        trackMuScore[iTrk] = -99.
        trackPiScore[iTrk] = -99.
        trackPrScore[iTrk] = -99.
        trackComp[iTrk] = -1.
        continue

      prongImage = makeImage(prong_vv).to(args.device)
      prongCNN_out = model(prongImage)
      trackClassified[iTrk] = 1
      trackPID[iTrk] = getPID(prongCNN_out[0].argmax(1).item())
      trackElScore[iTrk] = prongCNN_out[0][0][0].item()
      trackPhScore[iTrk] = prongCNN_out[0][0][1].item()
      trackMuScore[iTrk] = prongCNN_out[0][0][2].item()
      trackPiScore[iTrk] = prongCNN_out[0][0][3].item()
      trackPrScore[iTrk] = prongCNN_out[0][0][4].item()
      trackComp[iTrk] = prongCNN_out[1].item()

      if args.isMC:
        pdg, purity, completeness, allPdgs, allPurities = getMCProngParticle(prong_vv, mcpg, mcpm, adc_v)
        trackTruePID[iTrk] = pdg
        trackTruePurity[iTrk] = purity
        trackTrueComp[iTrk] = completeness
        trackTrueElPurity[iTrk] = 0.
        trackTruePhPurity[iTrk] = 0.
        trackTrueMuPurity[iTrk] = 0.
        trackTruePiPurity[iTrk] = 0.
        trackTruePrPurity[iTrk] = 0.
        for iTru, current_pdg in enumerate(allPdgs):
          if current_pdg == 11:
            trackTrueElPurity[iTrk] = allPurities[iTru]
          if current_pdg == 22:
            trackTruePhPurity[iTrk] = allPurities[iTru]
          if current_pdg == 13:
            trackTrueMuPurity[iTrk] = allPurities[iTru]
          if current_pdg == 211:
            trackTruePiPurity[iTrk] = allPurities[iTru]
          if current_pdg == 2212:
            trackTruePrPurity[iTrk] = allPurities[iTru]
        if args.makePlots:
          plotImage(prongImage, run[0], subrun[0], event[0], pdg, purity, completeness,
                    trackCharge[iTrk], trackComp[iTrk], trackElScore[iTrk], trackPhScore[iTrk],
                    trackMuScore[iTrk], trackPiScore[iTrk], trackPrScore[iTrk])
      else:
        trackTruePID[iTrk] = 0
        trackTruePurity[iTrk] = -1.
        trackTrueComp[iTrk] = -1.
        trackTrueElPurity[iTrk] = -1.
        trackTruePhPurity[iTrk] = -1.
        trackTrueMuPurity[iTrk] = -1.
        trackTruePiPurity[iTrk] = -1.
        trackTruePrPurity[iTrk] = -1.


    for iShw, shower in enumerate(vertex.shower_v):

      vertexNHits += shower.size()
      showerIsSecondary[iShw] = vertex.shower_isSecondary_v[iShw]
      showerNHits[iShw] = shower.size()
      showerCharge[iShw], vertexPixels, vertexCharge = addClusterCharge(iolcv,shower,vertexPixels,vertexCharge, 10.)
      showerPl0E[iShw] = vertex.shower_plane_mom_vv[iShw][0].E()
      showerPl1E[iShw] = vertex.shower_plane_mom_vv[iShw][1].E()
      showerPl2E[iShw] = vertex.shower_plane_mom_vv[iShw][2].E()
      showerCosTheta[iShw] = getCosThetaBeamShower(vertex.shower_trunk_v[iShw])
      showerDistToVtx[iShw] = getVertexDistance(vertex.shower_trunk_v[iShw].Vertex(), vertex)

      cropPt = vertex.shower_trunk_v[iShw].Vertex()
      prong_vv = flowTriples.make_cropped_initial_sparse_prong_image_reco(adc_v,thrumu_v,shower,cropPt,10.,512,512)
      skip = False
      for p in range(3):
        if prong_vv[p].size() < 10:
          skip = True
          break
      if skip:
        showerClassified[iShw] = 0
        showerPID[iShw] = 0
        showerElScore[iShw] = -99.
        showerPhScore[iShw] = -99.
        showerMuScore[iShw] = -99.
        showerPiScore[iShw] = -99.
        showerPrScore[iShw] = -99.
        showerComp[iShw] = -1.
        continue

      prongImage = makeImage(prong_vv).to(args.device)
      prongCNN_out = model(prongImage)
      showerClassified[iShw] = 1
      showerPID[iShw] = getPID(prongCNN_out[0].argmax(1).item())
      showerElScore[iShw] = prongCNN_out[0][0][0].item()
      showerPhScore[iShw] = prongCNN_out[0][0][1].item()
      showerMuScore[iShw] = prongCNN_out[0][0][2].item()
      showerPiScore[iShw] = prongCNN_out[0][0][3].item()
      showerPrScore[iShw] = prongCNN_out[0][0][4].item()
      showerComp[iShw] = prongCNN_out[1].item()

      if args.isMC:
        pdg, purity, completeness, allPdgs, allPurities = getMCProngParticle(prong_vv, mcpg, mcpm, adc_v)
        showerTruePID[iShw] = pdg
        showerTruePurity[iShw] = purity
        showerTrueComp[iShw] = completeness
        showerTrueElPurity[iShw] = 0.
        showerTruePhPurity[iShw] = 0.
        showerTrueMuPurity[iShw] = 0.
        showerTruePiPurity[iShw] = 0.
        showerTruePrPurity[iShw] = 0.
        for iTru, current_pdg in enumerate(allPdgs):
          if current_pdg == 11:
            showerTrueElPurity[iShw] = allPurities[iTru]
          if current_pdg == 22:
            showerTruePhPurity[iShw] = allPurities[iTru]
          if current_pdg == 13:
            showerTrueMuPurity[iShw] = allPurities[iTru]
          if current_pdg == 211:
            showerTruePiPurity[iShw] = allPurities[iTru]
          if current_pdg == 2212:
            showerTruePrPurity[iShw] = allPurities[iTru]
        if args.makePlots:
          plotImage(prongImage, run[0], subrun[0], event[0], pdg, purity, completeness,
                    showerCharge[iShw], showerComp[iShw], showerElScore[iShw], showerPhScore[iShw],
                    showerMuScore[iShw], showerPiScore[iShw], showerPrScore[iShw])
      else:
        showerTruePID[iShw] = 0
        showerTruePurity[iShw] = -1.
        showerTrueComp[iShw] = -1.
        showerTrueElPurity[iShw] = -1.
        showerTruePhPurity[iShw] = -1.
        showerTrueMuPurity[iShw] = -1.
        showerTruePiPurity[iShw] = -1.
        showerTruePrPurity[iShw] = -1.


    for i in range(nTracks[0]):
      trackHitFrac[i] = trackNHits[i] / (1.0*vertexNHits)
      trackChargeFrac[i] = trackCharge[i] / vertexCharge

    for i in range(nShowers[0]):
      showerHitFrac[i] = showerNHits[i] / (1.0*vertexNHits)
      showerChargeFrac[i] = showerCharge[i] / vertexCharge

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

if args.isMC and args.makePlots:
  outpdf.close()

