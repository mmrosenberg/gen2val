
import os,sys,argparse

import ROOT as rt
import uproot

from larlite import larlite
from larcv import larcv

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser("2D larcv event display script")
parser.add_argument("-i", "--larcv_file", type=str, required=True, help="input larcv images file")
parser.add_argument("-w","--wireName", type=str, default="wire", help="wire tree name")
parser.add_argument("-e","--event", type=int, default=0, help="first entry to plot")
parser.add_argument("-m","--min", type=float, default=-9999., help="minimum pixel value (none if < -999)")
parser.add_argument("-M","--max", type=float, default=-9999., help="maximum pixel value (none if < -999)")
args = parser.parse_args()


def plotImage(adc_v):
  if args.min < -999:
    pltmin = None
  else:
    pltmin = args.min
  if args.max < -999:
    pltmax = None
  else:
    pltmax = args.max
  ticksize=6
  fig = plt.figure(0, clear=True)
  fig.subplots_adjust(left=0.07, right=0.97, bottom=0.07, top=0.95, hspace=0.33)
  for p in range(3):
    plt.subplot(3,1,p+1)
    img = np.zeros((adc_v[p].meta().rows(), adc_v[p].meta().cols()))
    for r in range(adc_v[p].meta().rows()):
      for c in range(adc_v[p].meta().cols()):
        img[r][c] = adc_v[p].pixel(r, c)
    plt.imshow(img, vmin=pltmin, vmax=pltmax, cmap='jet')
    plt.xticks(fontsize=ticksize)
    plt.yticks(fontsize=ticksize)
  plt.show()
  input("Press Enter to continue...")
  return

iolcv = larcv.IOManager(larcv.IOManager.kREAD, "larcv", larcv.IOManager.kTickBackward)
iolcv.add_in_file(args.larcv_file)
iolcv.reverse_all_products()
iolcv.initialize()

for i in range(args.event, iolcv.get_n_entries()):
  iolcv.read_entry(i)
  evtImage2D = iolcv.get_data(larcv.kProductImage2D, args.wireName)
  plotImage(evtImage2D.Image2DArray())
