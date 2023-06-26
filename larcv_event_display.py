
import os,sys,argparse

import ROOT as rt
import uproot

from larlite import larlite
from larcv import larcv

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

rootutils = larcv.rootutils()

parser = argparse.ArgumentParser("2D larcv event display script")
parser.add_argument("-i", "--larcv_file", type=str, required=True, help="input larcv images file")
parser.add_argument("-o", "--output_dir", type=str, default="", help="save plot to this output directory (optional)")
parser.add_argument("-w","--wireName", type=str, default="wire", help="wire tree name")
parser.add_argument("-en","--entry", type=int, default=0, help="first entry to plot")
parser.add_argument("-m","--min", type=float, default=-9999., help="minimum pixel value (none if < -999)")
parser.add_argument("-M","--max", type=float, default=-9999., help="maximum pixel value (none if < -999)")
parser.add_argument("-r", "--run", type=int, default=-1, help="plot event with this run (if > 0)")
parser.add_argument("-sr", "--subrun", type=int, default=-1, help="plot event with this subrun (if > 0)")
parser.add_argument("-e", "--event", type=int, default=-1, help="plot event with this event id (if > 0)")
parser.add_argument("--auto_zoom", help="zoom in on in-time hits above threshold", action="store_true")
parser.add_argument("--root", help="plot with root instead of matplotlib", action="store_true")
parser.add_argument("-zr1", "--zoom_r1", type=int, default=-1, help="optional plot zoom coord: min row")
parser.add_argument("-zr2", "--zoom_r2", type=int, default=-1, help="optional plot zoom coord: max row")
parser.add_argument("-z0c1", "--zoom_p0c1", type=int, default=-1, help="optional plot zoom coord: plane 0 min column")
parser.add_argument("-z0c2", "--zoom_p0c2", type=int, default=-1, help="optional plot zoom coord: plane 0 max column")
parser.add_argument("-z1c1", "--zoom_p1c1", type=int, default=-1, help="optional plot zoom coord: plane 1 min column")
parser.add_argument("-z1c2", "--zoom_p1c2", type=int, default=-1, help="optional plot zoom coord: plane 1 max column")
parser.add_argument("-z2c1", "--zoom_p2c1", type=int, default=-1, help="optional plot zoom coord: plane 2 min column")
parser.add_argument("-z2c2", "--zoom_p2c2", type=int, default=-1, help="optional plot zoom coord: plane 2 max column")
parser.add_argument("--look_up_event", help="print additional bnb5e19 event info in title", action="store_true")
args = parser.parse_args()

z_r1, z_r2 = args.zoom_r1, args.zoom_r2
z_p0c1, z_p0c2 = args.zoom_p0c1, args.zoom_p0c2
z_p1c1, z_p1c2 = args.zoom_p1c1, args.zoom_p1c2
z_p2c1, z_p2c2 = args.zoom_p2c1, args.zoom_p2c2

outfile = rt.TFile("event_display_hists.root","RECREATE")
rt.gStyle.SetOptStat(0)
rt.gROOT.SetBatch(rt.kTRUE)

def getBNB5e19EventString(lineBreaks=True):
  selDir="/home/matthew/microboone/tufts/gen2val/selection_output"
  with open("%s/selected_bnb5e19_events_sorted_with_WC_recoE_plotBounds.txt"%selDir,"r") as evtFile:
    for line in evtFile:
      if "fileid" in line:
        continue
      data = line.split()
      if int(data[1]) == args.run and int(data[2]) == args.subrun and int(data[3]) == args.event:
        foundInGen2 = int(data[4])
        foundInWC = int(data[5])
        foundInGen2Str = "Yes" if foundInGen2 else "No"
        foundInWCStr = "Yes" if foundInWC else "No"
        recoNuE = data[6]+" MeV"
        break
  if lineBreaks:
    return "run %i subrun %i event %i \n DL Gen2 reco neutrino energy: %s \n found in DL Gen2: %s \n found in wire cell: %s"%(args.run,args.subrun,args.event,recoNuE,foundInGen2Str,foundInWCStr)
  return "run %i subrun %i event %i / DL Gen2 reco neutrino energy: %s / found in DL Gen2: %s / found in wire cell: %s"%(args.run,args.subrun,args.event,recoNuE,foundInGen2Str,foundInWCStr)

def getBNB5e19EventString_visE(lineBreaks=True):
  selDir="/home/matthew/microboone/tufts/gen2val/selection_output"
  with open("%s/selected_bnb5e19_events_sorted_withWC.txt"%selDir,"r") as evtFile1:
    for line in evtFile1:
      if "fileid" in line:
        continue
      data = line.split()
      if int(data[1]) == args.run and int(data[2]) == args.subrun and int(data[3]) == args.event:
        foundInGen2 = int(data[4])
        foundInWC = int(data[5])
        foundInGen2Str = "Yes" if foundInGen2 else "No"
        foundInWCStr = "Yes" if foundInWC else "No"
        if foundInWC and not foundInGen2:
          with open("%s/selected_bnb5e19_WConly_events.txt"%selDir) as wcEvtFile:
            for wcLine in wcEvtFile:
              if "fileid" in wcLine:
                continue
              wcData = wcLine.split()
              if int(wcData[1]) == args.run and int(wcData[2]) == args.subrun and int(wcData[3]) == args.event:
                visE = float(wcData[4])
                break
        else:
          with open("%s/selected_bnb5e19_events_sorted.txt"%selDir,"r") as gen2EvtFile:
            for gen2Line in gen2EvtFile:
              if "fileid" in gen2Line:
                continue
              gen2Data = gen2Line.split()
              if int(gen2Data[1]) == args.run and int(gen2Data[2]) == args.subrun and int(gen2Data[3]) == args.event:
                visE = float(gen2Data[4])
                break
        break
  if lineBreaks:
    return "run %i subrun %i event %i \n reco gen2 vertex visible energy: %.2f \n found in DL Gen2: %s \n found in wire cell: %s"%(args.run,args.subrun,args.event,visE,foundInGen2Str,foundInWCStr)
  return "run %i subrun %i event %i / reco gen2 vertex visible energy: %.2f / found in DL Gen2: %s / found in wire cell: %s"%(args.run,args.subrun,args.event,visE,foundInGen2Str,foundInWCStr)

def zoom_in(img, p):
  if p == 0:
    img = img[z_r1:z_r2,z_p0c1:z_p0c2]
  if p == 1:
    img = img[z_r1:z_r2,z_p1c1:z_p1c2]
  if p == 2:
    img = img[z_r1:z_r2,z_p2c1:z_p2c2]
  return img

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
  zoom = True if (z_p0c1 >= 0) else False
  if zoom:
    fig.subplots_adjust(left=0.053, right=0.98, wspace=0.28)
  else:
    fig.subplots_adjust(left=0.07, right=0.97, bottom=0.07, top=0.95, hspace=0.33)
  for p in range(3):
    if zoom:
      plt.subplot(1,3,p+1)
    else:
      plt.subplot(3,1,p+1)
    img = np.zeros((adc_v[p].meta().rows(), adc_v[p].meta().cols()))
    for r in range(adc_v[p].meta().rows()):
      for c in range(adc_v[p].meta().cols()):
        img[r][c] = adc_v[p].pixel(r, c)
    if zoom:
      img = zoom_in(img, p)
    plt.imshow(img, vmin=pltmin, vmax=pltmax, cmap='jet')
    plt.xticks(fontsize=ticksize)
    plt.yticks(fontsize=ticksize)
  if zoom:
    if args.look_up_event:
      plt.suptitle(getBNB5e19EventString())
    else:
      plt.suptitle("run %i subrun %i event %i"%(args.run,args.subrun,args.event))
  if args.output_dir != "":
    outpath = "%s/r%i_sr%i_e%i_%s_event_display.png"%(args.output_dir,args.run,args.subrun,args.event,args.wireName)
    if zoom:
      outpath = outpath.replace(".png","_zoom.png")
    plt.savefig(outpath, dpi=300)
  else:
    plt.show()
    input("Press Enter to continue...")
  return

def zoom_in_root(hists):
  for p in range(3):
    zr1, zr2 = z_r1*6 + 2400, z_r2*6 + 2400
    if p == 0:
      zc1, zc2 = z_p0c1, z_p0c2
    if p == 1:
      zc1, zc2 = z_p1c1, z_p1c2
    if p == 2:
      zc1, zc2 = z_p2c1, z_p2c2
    nBinsX, nBinsY = 0, 0
    xBinLow, xBinHigh, yBinLow, yBinHigh = 0, 0, 0, 0
    firstBin = True
    for i in range(1,hists[p].GetNbinsX()+1):
      binLow = hists[p].GetXaxis().GetBinLowEdge(i)
      binHigh = binLow + hists[p].GetXaxis().GetBinWidth(i)
      #print(binLow, binHigh, zc1, zc2)
      if binLow >= zc1 and binHigh < zc2:
        if firstBin:
          xBinLow = i
          firstBin = False
        xBinHigh = i
        nBinsX += 1
    firstBin = True
    for i in range(1,hists[p].GetNbinsY()+1):
      binLow = hists[p].GetYaxis().GetBinLowEdge(i)
      binHigh = binLow + hists[p].GetYaxis().GetBinWidth(i)
      #print(binLow, binHigh, zr1, zr2)
      if binLow >= zr1 and binHigh < zr2:
        if firstBin:
          yBinLow = i
          firstBin = False
        yBinHigh = i
        nBinsY += 1
    xLow = hists[p].GetXaxis().GetBinLowEdge(xBinLow)
    xHigh = hists[p].GetXaxis().GetBinLowEdge(xBinHigh) + hists[p].GetXaxis().GetBinWidth(xBinHigh)
    yLow = hists[p].GetYaxis().GetBinLowEdge(yBinLow)
    yHigh = hists[p].GetYaxis().GetBinLowEdge(yBinHigh) + hists[p].GetYaxis().GetBinWidth(yBinHigh)
    newHist = rt.TH2D("zoomHist%i"%p, "", nBinsX,xLow,xHigh, nBinsY,yLow,yHigh)
    for iX in range(nBinsX):
      for iY in range(nBinsY):
        newHist.SetBinContent(iX+1,iY+1,hists[p].GetBinContent(xBinLow+iX,yBinLow+iY))
    hists[p] = newHist
  return hists

def plotImageRoot(adc_v):
  if args.min > 0:
    for p in range(3):
      adc_v[p].threshold(args.min, 0.)
  hists = rootutils.as_th2d_v(adc_v, "larcv_hist")
  zoom = True if (z_p0c1 >= 0) else False
  cnv = rt.TCanvas("cnv","cnv")
  if zoom:
    cnv.Divide(3,1)
    hists = zoom_in_root(hists)
    cnv.SetCanvasSize(1800,600)
  else:
    cnv.Divide(1,3)
    cnv.SetCanvasSize(1800,1800)
  cnv.cd(1)
  hists[0].Draw("COLZ")
  cnv.cd(2)
  hists[1].Draw("COLZ")
  cnv.cd(3)
  hists[2].Draw("COLZ")
  cnv.cd()
  titlePad = rt.TPad("all","all",0,0,1,1)
  titlePad.SetFillStyle(4000)
  titlePad.Draw()
  titlePad.cd()
  lat = rt.TLatex()
  if args.look_up_event:
    if zoom:
      lat.DrawLatexNDC(.05,.95,getBNB5e19EventString(False))
    else:
      lat.DrawLatexNDC(.4,.95,getBNB5e19EventString(False))
  else:
    lat.DrawLatexNDC(.4,.95,"run %i subrun %i event %i"%(args.run,args.subrun,args.event))
  outfile.cd()
  cnv.Write()
  outdir = args.output_dir if (args.output_dir != "") else "./"
  outpath = "%s/r%i_sr%i_e%i_%s_event_display.png"%(outdir,args.run,args.subrun,args.event,args.wireName)
  cnv.SaveAs(outpath)


def set_zoom_coords(iolcv):
  global z_r1
  global z_r2
  global z_p0c1
  global z_p0c2
  global z_p1c1
  global z_p1c2
  global z_p2c1
  global z_p2c2
  threshold = 50
  adc_v = iolcv.get_data(larcv.kProductImage2D, "wire").Image2DArray()
  adc_v_thrumu = iolcv.get_data(larcv.kProductImage2D, "thrumu").Image2DArray()
  pixBounds = []
  for p in range(3):
    rmin, rmax = 9999, -9999
    cmin, cmax = 9999, -9999
    for r in range(adc_v[p].meta().rows()):
      for c in range(adc_v[p].meta().cols()):
        if (adc_v[p].pixel(r, c) - adc_v_thrumu[p].pixel(r, c)) > threshold:
          rmin = r if (r < rmin) else rmin
          rmax = r if (r > rmax) else rmax
          cmin = c if (c < cmin) else cmin
          cmax = c if (c > cmax) else cmax
    pixBounds.append([ [rmin, rmax], [cmin, cmax] ])
  #imgRange = 0
  #for p in range(3):
  #  rowRange = pixBounds[p][0][1] - pixBounds[p][0][0]
  #  colRange = pixBounds[p][1][1] - pixBounds[p][1][0]
  #  imgRange = rowRange if (rowRange > imgRange) else imgRange
  #  imgRange = colRange if (colRange > imgRange) else imgRange
  row_mid = (pixBounds[2][0][0] + pixBounds[2][0][1])/2
  print("pixBounds:", pixBounds)
  print("row_mid:", row_mid)
  #print("imgRange:", imgRange)
  #z_r1 = int(row_mid - 0.6*imgRange)
  #z_r2 = int(row_mid + 0.6*imgRange)
  z_r1 = int(row_mid - 0.6*(pixBounds[2][0][1] - pixBounds[2][0][0]))
  z_r2 = int(row_mid + 0.6*(pixBounds[2][0][1] - pixBounds[2][0][0]))
  nRows = adc_v[2].meta().rows()
  #if imgRange > nRows:
  #  z_r1, z_r2 = 0, nRows
  if (pixBounds[2][0][1] - pixBounds[2][0][0]) > nRows:
    z_r1, z_r2 = 0, nRows
  else:
    if z_r1 < 0:
      z_r1, z_r2 = 0, (z_r2 - z_r1)
    if z_r2 > nRows:
      z_r1, z_r2 = (nRows - (z_r2 - z_r1)), nRows
  for p in range(3):
    nCols = adc_v[p].meta().cols()
    col_mid = (pixBounds[p][1][0] + pixBounds[p][1][1])/2
    #z_c1 = int(col_mid - 0.6*imgRange)
    #z_c2 = int(col_mid + 0.6*imgRange)
    z_c1 = int(col_mid - 0.6*(pixBounds[p][1][1] - pixBounds[p][1][0]))
    z_c2 = int(col_mid + 0.6*(pixBounds[p][1][1] - pixBounds[p][1][0]))
    if z_c1 < 0:
      z_c1, z_c2 = 0, (z_c2 - z_c1)
    if z_c2 > nCols:
      z_c1, z_c2 = (nCols - (z_c2 - z_c1)), nCols
    if p == 0:
      z_p0c1, z_p0c2 = z_c1, z_c2
    if p == 1:
      z_p1c1, z_p1c2 = z_c1, z_c2
    if p == 2:
      z_p2c1, z_p2c2 = z_c1, z_c2



iolcv = larcv.IOManager(larcv.IOManager.kREAD, "larcv", larcv.IOManager.kTickBackward)
iolcv.add_in_file(args.larcv_file)
iolcv.reverse_all_products()
iolcv.initialize()

have_larlite = True
try:
  ioll = larlite.storage_manager(larlite.storage_manager.kREAD)
  ioll.add_in_filename(args.larcv_file)
  ioll.open()
except:
  have_larlite = False


for i in range(args.entry, iolcv.get_n_entries()):
  iolcv.read_entry(i)
  if have_larlite:
    ioll.go_to(i)
    r = ioll.run_id()
    sr = ioll.subrun_id()
    e = ioll.event_id()
    if args.run >= 0 and args.subrun >= 0 and args.event >= 0:
      if args.run != r or args.subrun != sr or args.event != e:
        print("skipping entry %i (run %i, subrun %i, event %i)"%(i,r,sr,e))
        continue
  print("plotting image for run %i subrun %i event %i"%(r, sr, e))
  if args.auto_zoom:
    set_zoom_coords(iolcv)
  print(z_r1, z_r2, z_p0c1, z_p0c2, z_p1c1, z_p1c2, z_p2c1, z_p2c2)
  if args.wireName == "inTime":
    adc_v = iolcv.get_data(larcv.kProductImage2D, "wire").Image2DArray()
    adc_v_thrumu = iolcv.get_data(larcv.kProductImage2D, "thrumu").Image2DArray()
    for p in range(3):
      adc_v[p].overlay(adc_v_thrumu[p], larcv.Image2D.kSubtract)
  else:
    adc_v = iolcv.get_data(larcv.kProductImage2D, args.wireName).Image2DArray()
  if args.root:
    plotImageRoot(adc_v)
  else:
    plotImage(adc_v)


