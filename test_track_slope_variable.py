
import os,sys,argparse
import ROOT as rt
import uproot
from larlite import larlite
from larflow import larflow
from larcv import larcv
import numpy as np
import matplotlib.pyplot as plt
from math import sqrt, acos, cos, sin, pi
from scipy.stats import linregress


parser = argparse.ArgumentParser("Plot pixSum vs hit distance for input track")
parser.add_argument("-f", "--kpsfile", required=True, type=str, help="input kpsreco file")
parser.add_argument("-ft", "--mdlfile", required=True, type=str, help="input merged_dlreco file")
parser.add_argument("-r", "--run", required=True, type=int, help="run")
parser.add_argument("-sr", "--subrun", required=True, type=int, help="subrun")
parser.add_argument("-e", "--event", required=True, type=int, help="event")
parser.add_argument("-v", "--vertexID", required=True, type=int, help="index for neutrino vertex")
parser.add_argument("-t", "--trackID", required=True, type=int, help="index for track")
args = parser.parse_args()

def getHitDistance(trackPoint, hit):
  delXsq = (hit[0] - trackPoint.X())**2
  delYsq = (hit[1] - trackPoint.Y())**2
  delZsq = (hit[2] - trackPoint.Z())**2
  return sqrt(delXsq + delYsq + delZsq)

def getPixelVal(hit, image2Dvec, p=2, sumPlanes=False):
  row = (hit.tick - 2400)//6
  if sumPlanes:
    pixSum = 0.
    for p in range(3):
      pixSum += image2Dvec[p].pixel(row, hit.targetwire[p])
    return pixSum
  return image2Dvec[p].pixel(row, hit.targetwire[p])

def fitTrackPoints(dists, pixSums, cutoffL, cutoffH):
  dists_fit = []
  pixSums_fit = []
  for iH, dist in enumerate(dists):
    if dist > cutoffL and dist < cutoffH:
      dists_fit.append(dist)
      pixSums_fit.append(pixSums[iH])
  if len(dists_fit) > 5:
    fitRes = linregress(dists_fit, pixSums_fit)
    return fitRes.slope, fitRes.intercept
  return 0., 0.

def sortFunc(vals):
  return vals[0]

f = rt.TFile(args.kpsfile)
kpst = f.Get("KPSRecoManagerTree")

ioll = larlite.storage_manager(larlite.storage_manager.kREAD)
ioll.add_in_filename(args.mdlfile)
ioll.open()

iolcv = larcv.IOManager(larcv.IOManager.kREAD, "larcv", larcv.IOManager.kTickBackward)
iolcv.add_in_file(args.mdlfile)
iolcv.reverse_all_products()
iolcv.initialize()

for i in range(kpst.GetEntries()):

  kpst.GetEntry(i)

  if kpst.run != args.run or kpst.subrun != args.subrun or kpst.event != args.event:
    continue

  ioll.go_to(i)
  iolcv.read_entry(i)

  if kpst.run != ioll.run_id() or kpst.subrun != ioll.subrun_id() or kpst.event != ioll.event_id():
    print("EVENTS DON'T MATCH!!!")
    print("truth run/subrun/event: %i/%i/%i"%(ioll.run_id(),ioll.subrun_id(),ioll.event_id()))
    print("reco run/subrun/event: %i/%i/%i"%(kpst.run,kpst.subrun,kpst.event))
    continue

  print("plotting for run %i subrun %i event %i"%(kpst.run, kpst.subrun, kpst.event))

  
  evtImage2D = iolcv.get_data(larcv.kProductImage2D, "wire")
  image2Dvec = evtImage2D.Image2DArray()

  track = kpst.nuvetoed_v[args.vertexID].track_v[args.trackID]
  trackCluster = kpst.nuvetoed_v[args.vertexID].track_hitcluster_v[args.trackID]

  fitStart = 3.
  fitEnd = 10.
  minSlope = 1e9

  if track.Length() < fitStart:
    sys.exit("can't fit track starting from %fcm from vertex, track length = %f"%(fitStart, track.Length()))

  for p in range(3):

    trackPoints = []
    for hit in trackCluster:
      trackPoints.append([getHitDistance(track.Vertex(), hit), getPixelVal(hit, image2Dvec, p, False)])

    trackPoints.sort(key=sortFunc)
    dists = [ vals[0] for vals in trackPoints ]
    pixVals = [vals[1] for vals in trackPoints ]

    slope, intercept = fitTrackPoints(dists, pixVals, fitStart, fitEnd)
    if slope < minSlope:
      minSlope = slope
      fitIntercept = intercept
      fitPlaneDists = dists
      fitPlanePixVals = pixVals

  print("minSlope: ", minSlope)

  fitDists = np.linspace(fitStart,fitEnd,100)
  fitPoints = [minSlope*fitDists[iF] + fitIntercept for iF in range(100)]
  fitPoints = np.array(fitPoints)

  plt.plot(fitPlaneDists, fitPlanePixVals, 'b-')
  plt.plot(fitDists, fitPoints, 'r-')
  plt.show()




