
import os,sys,argparse
import ROOT as rt
from larlite import larlite
from larlite import larutil
from ublarcvapp import ublarcvapp
from helpers.larflowreco_ana_funcs import *
from math import sqrt
from array import array

parser = argparse.ArgumentParser("Compare Range vs KE for Different Particles")
parser.add_argument("-i", "--input", required=True, type=str, nargs="+", help="input truth files")
parser.add_argument("-o", "--output", type=str, default="analyze_range_mass_scaling_output.root", help="output file path")
args = parser.parse_args()

sce = larutil.SpaceChargeMicroBooNE()
mcNuVertexer = ublarcvapp.mctools.NeutrinoVertex()

outRootFile = rt.TFile(args.output, "RECREATE")

tree = rt.TTree("ParticleTree","ParticleTree")
pdg = array('f', [0.])
length = array('f', [0.])
KE = array('f', [0.])
tree.Branch("pdg", pdg, 'pdg/F')
tree.Branch("length", length, 'length/F')
tree.Branch("KE", KE, 'KE/F')

def getTrackLength(track):
  tracklen = 0.
  for i in range(track.size()-1):
    x1, y1, z1 = track[i].X(), track[i].Y(), track[i].Z()
    x2, y2, z2 = track[i+1].X(), track[i+1].Y(), track[i+1].Z()
    tracklen += sqrt(pow(x2-x1,2) + pow(y2-y1,2) + pow(z2-z1,2))
  return tracklen


for infile in args.input:

  ioll = larlite.storage_manager(larlite.storage_manager.kREAD)
  ioll.add_in_filename(infile)
  ioll.open()

  for ientry in range(ioll.get_entries()):
  
    ioll.go_to(ientry)

    mcNuVertex = mcNuVertexer.getPos3DwSCE(ioll, sce)
    trueVtxPos = rt.TVector3(mcNuVertex[0], mcNuVertex[1], mcNuVertex[2])

    if not isFiducial(trueVtxPos):
      continue

    mctracks = ioll.get_data(larlite.data.kMCTrack, "mcreco")

    for mctrack in mctracks:
      if abs(mctrack.PdgCode()) in [13, 211] and mctrack.Process() == 'primary' and isInDetector(mctrack.End()):
        pdg[0] = abs(mctrack.PdgCode())
        length[0] = getTrackLength(mctrack)
        mass = 105.7
        if abs(mctrack.PdgCode()) == 211:
          mass = 139.6
        KE[0] = mctrack.Start().E() - mass
        tree.Fill()

  ioll.close()


outRootFile.cd()
tree.Write("",rt.TObject.kOverwrite)
outRootFile.Close()

