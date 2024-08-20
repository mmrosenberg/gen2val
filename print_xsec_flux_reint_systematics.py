
import ROOT as rt
import argparse
from math import sqrt

parser = argparse.ArgumentParser("Print Systematics Uncertainties")
parser.add_argument("-t", "--tag", type=str, required=True, help="covariance matrix file tag")
args = parser.parse_args()

bins = []

if "CCnueInc" in args.tag and ("nuE" in args.tag or "lepP" in args.tag) and "chi2" not in args.tag:
  bins = [[0, 200], [200, 400], [400, 600], [600, 800], [800, 1000], [1000, 1200], [1200, 1400], [1400, 1600], [1600, 1800], [1800, 2000], [2000, 2200], [2200, 2400], [2400, 2600], [2600, 6000]]

if "CCnueInc" in args.tag and "nuE" in args.tag and "chi2" in args.tag:
  bins = [[0, 400], [400, 600], [600, 800], [800, 1000], [1000, 1200], [1200, 1400], [1400, 1600], [1600, 1800], [1800, 6000]]

if "CCnueInc" in args.tag and "lepP" in args.tag and "chi2" in args.tag:
  bins = [[0, 200], [200, 400], [400, 600], [600, 800], [800, 1000], [1000, 1200], [1200, 1400], [1400, 6000]]

if "CCnumuInc" in args.tag and ("nuE" in args.tag or "lepP" in args.tag) and "rebinned" not in args.tag:
  bins = [[0, 100], [100, 200], [200, 300], [300, 400], [400, 500], [500, 600], [600, 700], [700, 800], [800, 900], [900, 1000], [1000, 1100], [1100, 1200], [1200, 1300], [1300, 1400], [1400, 1500], [1500, 1600], [1600, 1700], [1700, 1800], [1800, 1900], [1900, 2000], [2000, 6000]]

if "CCnumuInc" in args.tag and "nuE" in args.tag and "rebinned" in args.tag:
  bins = [[0, 100], [100, 200], [200, 300], [300, 400], [400, 500], [500, 600], [600, 700], [700, 800], [800, 900], [900, 1000], [1000, 1100], [1100, 1200], [1200, 1300], [1300, 1400], [1400, 1500], [1500, 1600], [1600, 1700], [1700, 1800], [1800, 6000]]

if "CCnumuInc" in args.tag and "lepP" in args.tag and "rebinned" in args.tag:
  bins = [[0, 100], [100, 200], [200, 300], [300, 400], [400, 500], [500, 600], [600, 700], [700, 800], [800, 900], [900, 1000], [1000, 1100], [1100, 1200], [1200, 1300], [1300, 1400], [1400, 1500], [1500, 6000]]


if "cosTheta" in args.tag and "chi2" not in args.tag:
  bins = [[-1.0, -0.875], [-0.875, -0.75], [-0.75, -0.625], [-0.625, -0.5], [-0.5, -0.375], [-0.375, -0.25], [-0.25, -0.125], [-0.125, 0.0], [0.0, 0.125], [0.125, 0.25], [0.25, 0.375], [0.375, 0.5], [0.5, 0.625], [0.625, 0.75], [0.75, 0.875], [0.875, 1.0]]

if "CCnueInc" in args.tag and "cosTheta" in args.tag and "chi2" in args.tag:
  bins = [[-1.0, 0.375], [0.375, 0.5], [0.5, 0.625], [0.625, 0.75], [0.75, 0.875], [0.875, 1.0]]

if "nuE" in args.tag and "cutSet" not in args.tag and "chi2" not in args.tag and "rebinned" not in args.tag:
  fCovar = rt.TFile("systematics/%s/SBNfit_covariance_plots_%s.root"%(args.tag,args.tag.replace("_nuE","")))
else:
  fCovar = rt.TFile("systematics/%s/SBNfit_covariance_plots_%s.root"%(args.tag,args.tag))
hCovar_frac = fCovar.Get("coll_frac")

nbins = len(bins)

print("%s_noDet = {"%args.tag)

for i in range(1,nbins+1,1):
  uncert = sqrt(hCovar_frac.GetBinContent(i,i))
  if "cosTheta" in args.tag:
    if i < nbins:
      print("  %i:  (%f, %f, %f),"%(i,uncert,bins[i-1][0],bins[i-1][1]))
    else:
      print("  %i:  (%f, %f, %f)"%(i,uncert,bins[i-1][0],bins[i-1][1]))
  else:
    if i < nbins:
      print("  %i:  (%f, %i, %i),"%(i,uncert,bins[i-1][0],bins[i-1][1]))
    else:
      print("  %i:  (%f, %i, %i)"%(i,uncert,bins[i-1][0],bins[i-1][1]))

print("}")
