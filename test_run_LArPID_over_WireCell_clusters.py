
import os,sys,argparse

import ROOT as rt
import uproot

from larlite import larlite
from larlite import larutil
from ublarcvapp import ublarcvapp
from larcv import larcv
from larflow import larflow

import numpy as np

import torch
from torch import nn
from torch.utils.data import DataLoader
import torchvision.transforms as transforms

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf



parser = argparse.ArgumentParser("Run LArPID over WC prongs from nue file")
parser.add_argument("-f", "--nuefile", required=True, type=str, help="input single-event Wire-Cell nue file")
parser.add_argument("-t", "--imgfile", required=True, type=str, help="merged_dlreco file")
parser.add_argument("-o", "--outfile", type=str, default="", help="output pdf file name")
parser.add_argument("-m", "--model_path", type=str, required=True, help="path to prong CNN checkpoint file")
parser.add_argument("-d", "--device", type=str, default="cpu", help="gpu/cpu device")
parser.add_argument("--multiGPU", action="store_true", help="use multiple GPUs")
parser.add_argument("--interactive", action="store_true", help="display plots one by one rather than saving to pdf")
args = parser.parse_args()

sys.path.append(args.model_path[:args.model_path.find("/checkpoints")])
from models_instanceNorm_reco_2chan_quadTask import ResBlock, ResNet34
from normalization_constants import mean, std

mpl.rcParams['figure.dpi'] = 300
if not args.interactive:
  plotfilename = args.outfile
  if plotfilename == "":
    plotfilename = args.nuefile.replace(".root","_images.pdf")
  outpdf = matplotlib.backends.backend_pdf.PdfPages(plotfilename)


def getMCPartE(ioll, tid):
  mctracks = ioll.get_data(larlite.data.kMCTrack, "mcreco")
  mcshowers = ioll.get_data(larlite.data.kMCShower, "mcreco")
  for mcparticles in [mctracks, mcshowers]:
    for mcpart in mcparticles:
      if mcpart.TrackID() == tid:
        return mcpart.Start().E()
  return -1.

def getMCProngParticle(sparseimg_vv, mcpg, mcpm, adc_v, ioll):

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
  maxPartE = -1.
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
    maxPartE = getMCPartE(ioll, maxPartTID)
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
  return maxPartPDG, maxPartTID, maxPartE, maxPartI/totalPixI, maxPartComp, pdglist, puritylist


#def plotImage(X, r, sr, e, pdg, purity, completeness, visE, compPred, purPred, elScore, phScore, muScore, piScore, prScore):
#def plotImage(X, compPred, purPred, elScore, phScore, muScore, piScore, prScore):
def plotImage(X, wc_pdg, dl_pdg, dl_proc, compPred, purPred, true_pdg, true_pur, true_comp, elScore, phScore, muScore, piScore, prScore):
  #pltmin = None
  #pltmax = None
  pltmin = -1. 
  pltmax = 2.
  suptitlesize=6
  titlesize=8
  ticksize=6
  X0 = X[0:2].reshape(2, 512, 512)
  X1 = X[2:4].reshape(2, 512, 512)
  X2 = X[4:].reshape(2, 512, 512)
  fig = plt.figure(0, clear=True)
  plt.subplot(2,3,1)
  plt.imshow(X0.numpy()[0], vmin=pltmin, vmax=pltmax, cmap='jet')
  plt.title("plane 0 prong", fontsize=titlesize)
  plt.xticks(fontsize=ticksize)
  plt.yticks(fontsize=ticksize)
  plt.subplot(2,3,2)
  plt.imshow(X1.numpy()[0], vmin=pltmin, vmax=pltmax, cmap='jet')
  plt.title("plane 1 prong", fontsize=titlesize)
  plt.xticks(fontsize=ticksize)
  plt.yticks(fontsize=ticksize)
  plt.subplot(2,3,3)
  plt.imshow(X2.numpy()[0], vmin=pltmin, vmax=pltmax, cmap='jet')
  plt.title("plane 2 prong", fontsize=titlesize)
  plt.xticks(fontsize=ticksize)
  plt.yticks(fontsize=ticksize)
  plt.subplot(2,3,4)
  plt.imshow(X0.numpy()[1], vmin=pltmin, vmax=pltmax, cmap='jet')
  plt.title("plane 0 all", fontsize=titlesize)
  plt.xticks(fontsize=ticksize)
  plt.yticks(fontsize=ticksize)
  plt.subplot(2,3,5)
  plt.imshow(X1.numpy()[1], vmin=pltmin, vmax=pltmax, cmap='jet')
  plt.title("plane 1 all", fontsize=titlesize)
  plt.xticks(fontsize=ticksize)
  plt.yticks(fontsize=ticksize)
  plt.subplot(2,3,6)
  plt.imshow(X2.numpy()[1], vmin=pltmin, vmax=pltmax, cmap='jet')
  plt.title("plane 2 all", fontsize=titlesize)
  plt.xticks(fontsize=ticksize)
  plt.yticks(fontsize=ticksize)
  #plt.suptitle("Run %i Subrun %i Event %i  |  Prong: pdg %i, purity %.2f, completeness %.2f, visible energy %.2e \n e- score %.2f, photon score %.2f, mu score: %.2f, pi score %.2f, proton score %.2f, completeness prediction: %.2f, purity prediction: %.2f"%(r, sr, e, pdg, purity, completeness, visE, elScore, phScore, muScore, piScore, prScore, compPred, purPred), fontsize=suptitlesize)
  #plt.suptitle("e- score %.2f, photon score %.2f, mu score: %.2f, pi score %.2f, proton score %.2f \n completeness prediction: %.2f, purity prediction: %.2f"%(elScore, phScore, muScore, piScore, prScore, compPred, purPred), fontsize=suptitlesize)
  plt.suptitle("true pdg: %i, true purity: %.2f, true completeness: %.2f, Wire-Cell pdg: %i \n LArPID pdg: %i, LArPID process: %s, LArPID purity: %.2f, LArPID completeness: %.2f \n e- score %.2f, photon score %.2f, mu score: %.2f, pi score %.2f, proton score %.2f"%(true_pdg,true_pur,true_comp,wc_pdg,dl_pdg,dl_proc,purPred,compPred,elScore, phScore, muScore, piScore, prScore), fontsize=suptitlesize)
  if args.interactive:
    plt.show()
    input("Press Enter to continue...")
  else: 
    outpdf.savefig(fig)
  return



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

def getProcString(cnnProcClass):
  if cnnProcClass == 0:
    return "primary"
  if cnnProcClass == 1:
    return "secondary, neutral parent"
  if cnnProcClass == 2:
    return "secondary, charged parent"
  return "unknown class"


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


geo = larutil.GeometryHelper.GetME()
flowTriples = larflow.prep.FlowTriples()


fnue = rt.TFile(args.nuefile)
Trun = fnue.Get("Trun")
TMC = fnue.Get("TMC")
T_rec = fnue.Get("T_rec")

Trun.GetEntry(0)
run = Trun.runNo
subrun = Trun.subRunNo
event = Trun.eventNo

recoIDs = {}
TMC.GetEntry(0)
for i in range(TMC.mc_Ntrack):
  #recoIDs[ TMC.mc_id[i] ] = TMC.mc_pdg[i]
  recoIDs[ TMC.mc_id[i] ] = i

print("recoIDs:", recoIDs)

recoSPs = {}
recoSPkeys = []
recoSPskipped = []


ioll = larlite.storage_manager(larlite.storage_manager.kREAD)
ioll.add_in_filename(args.imgfile)
ioll.open()

iolcv = larcv.IOManager(larcv.IOManager.kREAD, "larcv", larcv.IOManager.kTickBackward)
iolcv.add_in_file(args.imgfile)
iolcv.reverse_all_products()
iolcv.initialize()


for ientry in range(ioll.get_entries()):

  ioll.go_to(ientry)

  if run != ioll.run_id() or subrun != ioll.subrun_id() or event != ioll.event_id():
    continue

  iolcv.read_entry(ientry)

  print("found matching larcv event for rsr = %i,%i,%i"%(run,subrun,event))

  evtImage2D = iolcv.get_data(larcv.kProductImage2D, "wire")
  csmImage2D = iolcv.get_data(larcv.kProductImage2D, "thrumu")
  adc_v = evtImage2D.Image2DArray()
  thrumu_v = csmImage2D.Image2DArray()

  mcpg = ublarcvapp.mctools.MCPixelPGraph()
  mcpg.set_adc_treename("wire")
  mcpg.buildgraph(iolcv, ioll)
  mcpm = ublarcvapp.mctools.MCPixelPMap()
  mcpm.set_adc_treename("wire")
  mcpm.buildmap(iolcv, mcpg)

  mcpg.printAllNodeInfo()

  for iSP in range(T_rec.GetEntries()):
    T_rec.GetEntry(iSP)
    if T_rec.real_cluster_id not in recoIDs:
      if T_rec.real_cluster_id not in recoSPskipped:
        recoSPskipped.append(T_rec.real_cluster_id)
      continue
    if T_rec.real_cluster_id not in recoSPs:
      recoSPs[T_rec.real_cluster_id] = [(T_rec.x,T_rec.y,T_rec.z)]
      recoSPkeys.append(T_rec.real_cluster_id)
    else:
      recoSPs[T_rec.real_cluster_id].append( (T_rec.x,T_rec.y,T_rec.z) )

  print("recoSPkeys:", recoSPkeys)
  print("recoSPskipped:", recoSPskipped)

  for cluster in recoSPs:

    print("analyzing WC cluster with ID", cluster)

    lfcluster = larlite.larflowcluster()
    for sp in recoSPs[cluster]:
      hit = larlite.larflow3dhit()
      for p in range(adc_v.size()):
        center2D = geo.Point_3Dto2D(sp[0],sp[1],sp[2],p)
        if p == 0:
          hit.tick = int(center2D.t/geo.TimeToCm() + 3200.)
        if p < hit.targetwire.size():
          hit.targetwire[p] = int(center2D.w/geo.WireToCm())
        else:
          hit.targetwire.push_back(int(center2D.w/geo.WireToCm()))
      lfcluster.push_back(hit)

    iMC = recoIDs[cluster]
    print("  WC pdg:", TMC.mc_pdg[iMC])
    if TMC.mc_pdg[iMC] in [11, 22]:
      cropPt = rt.TVector3(TMC.mc_startXYZT[iMC*4],TMC.mc_startXYZT[iMC*4+1],TMC.mc_startXYZT[iMC*4+2])
    else:
      cropPt = rt.TVector3(TMC.mc_endXYZT[iMC*4],TMC.mc_endXYZT[iMC*4+1],TMC.mc_endXYZT[iMC*4+2])

    prong_vv = flowTriples.make_cropped_initial_sparse_prong_image_reco(adc_v,thrumu_v,lfcluster,cropPt,10.,512,512)

    skip = False
    for p in range(3):
      if prong_vv[p].size() < 5:
        skip = True
        break
    if skip:
      continue

    with torch.no_grad():
      prongImage = makeImage(prong_vv).to(args.device)
      prongCNN_out = model(prongImage)
      pdg, trackid, trueE, purity, completeness, allPdgs, allPurities = getMCProngParticle(prong_vv, mcpg, mcpm, adc_v, ioll)
      print("  LArPID PID:", getPID(prongCNN_out[0].argmax(1).item()))
      print("  LArPID ElScore:", prongCNN_out[0][0][0].item())
      print("  LArPID PhScore:", prongCNN_out[0][0][1].item())
      print("  LArPID MuScore:", prongCNN_out[0][0][2].item())
      print("  LArPID PiScore:", prongCNN_out[0][0][3].item())
      print("  LArPID PrScore:", prongCNN_out[0][0][4].item())
      print("  LArPID Comp:", prongCNN_out[1].item())
      print("  LArPID Purity:", prongCNN_out[2].item())
      print("  LArPID Process:", prongCNN_out[3].argmax(1).item())
      print("  LArPID PrimaryScore:", prongCNN_out[3][0][0].item())
      print("  LArPID FromNeutralScore:", prongCNN_out[3][0][1].item())
      print("  LArPID FromChargedScore:", prongCNN_out[3][0][2].item())
      print("  true pdg:", pdg)
      print("  true purity:", purity)
      print("  true completeness:", completeness)
      print("  true allPDGs:", allPdgs)
      print("  true allPurities:", allPurities)
      #plotImage(prongImage[0].cpu(), prongCNN_out[1].item(), prongCNN_out[2].item(), prongCNN_out[0][0][0].item(), prongCNN_out[0][0][1].item(), prongCNN_out[0][0][2].item(), prongCNN_out[0][0][3].item(), prongCNN_out[0][0][4].item())
      plotImage(prongImage[0].cpu(), TMC.mc_pdg[iMC], getPID(prongCNN_out[0].argmax(1).item()), getProcString(prongCNN_out[3].argmax(1).item()), prongCNN_out[1].item(), prongCNN_out[2].item(), pdg, purity, completeness, prongCNN_out[0][0][0].item(), prongCNN_out[0][0][1].item(), prongCNN_out[0][0][2].item(), prongCNN_out[0][0][3].item(), prongCNN_out[0][0][4].item())


if not args.interactive:
  outpdf.close()


