
import os,sys,argparse

import ROOT as rt

from larlite import larlite
from larcv import larcv

parser = argparse.ArgumentParser("compare larcv image pixel values")
parser.add_argument("-a", "--larcv_file1", type=str, required=True, help="first input larcv images file")
parser.add_argument("-b", "--larcv_file2", type=str, required=True, help="second input larcv images file")
parser.add_argument("-e", "--entry", type=int, default=0, help="entry in input files to print pixel values for")
parser.add_argument("--tickforward1", help="read first larcv image file with tickforward option", action="store_true")
parser.add_argument("--tickforward2", help="read second larcv image file with tickforward option", action="store_true")
parser.add_argument("--nomc", help="don't read MC truth images", action="store_true")
args = parser.parse_args()

if args.tickforward1:
  iolcv1 = larcv.IOManager(larcv.IOManager.kREAD, "larcv", larcv.IOManager.kTickForward)
  iolcv1.add_in_file(args.larcv_file1)
  iolcv1.initialize()
else:
  iolcv1 = larcv.IOManager(larcv.IOManager.kREAD, "larcv", larcv.IOManager.kTickBackward)
  iolcv1.add_in_file(args.larcv_file1)
  iolcv1.reverse_all_products()
  iolcv1.initialize()

if args.tickforward2:
  iolcv2 = larcv.IOManager(larcv.IOManager.kREAD, "larcv", larcv.IOManager.kTickForward)
  iolcv2.add_in_file(args.larcv_file2)
  iolcv2.initialize()
else:
  iolcv2 = larcv.IOManager(larcv.IOManager.kREAD, "larcv", larcv.IOManager.kTickBackward)
  iolcv2.add_in_file(args.larcv_file2)
  iolcv2.reverse_all_products()
  iolcv2.initialize()

ioll1 = larlite.storage_manager(larlite.storage_manager.kREAD)
ioll1.add_in_filename(args.larcv_file1)
ioll1.open()

ioll2 = larlite.storage_manager(larlite.storage_manager.kREAD)
ioll2.add_in_filename(args.larcv_file2)
ioll2.open()

iolcv1.read_entry(args.entry)
ioll1.go_to(args.entry)
iolcv2.read_entry(args.entry)
ioll2.go_to(args.entry)

if ioll1.run_id() != ioll2.run_id() or ioll1.subrun_id() != ioll2.subrun_id() or ioll1.event_id() != ioll2.event_id():
  print("MISMATCH IN EVENT NUMBERS!!!")
  print("  file 1: read in run %i subrun %i event %i"%(ioll1.run_id(),ioll1.subrun_id(),ioll1.event_id()))
  print("  file 2: read in run %i subrun %i event %i"%(ioll2.run_id(),ioll2.subrun_id(),ioll2.event_id()))
  sys.exit()

images = ["wire","thrumu","ancestor","instance","larflow","segment"]
if args.nomc:
  images = ["wire","thrumu"]

for image_type in images:
  adc_v1 = iolcv1.get_data(larcv.kProductImage2D, image_type).Image2DArray()
  adc_v2 = iolcv2.get_data(larcv.kProductImage2D, image_type).Image2DArray()
  if adc_v1.size() != adc_v2.size():
    print("MISMATCH IN PLANE COUNT FOR %s IMAGE!!!"%image_type)
    print("  file1 adc Image2DArray size: %i"%adc_v1.size())
    print("  file2 adc Image2DArray size: %i"%adc_v2.size())
    sys.exit()
  for p in range(adc_v1.size()):
    if adc_v1[p].meta().rows() != adc_v2[p].meta().rows() or adc_v1[p].meta().cols() != adc_v2[p].meta().cols():
      print("MISMATCH IN %s IMAGE PLANE %i META DATA!!!")
      print("  file1: %i rows, %i columns"%(adc_v1[p].meta().rows(),adc_v1[p].meta().cols()))
      print("  file2: %i rows, %i columns"%(adc_v2[p].meta().rows(),adc_v2[p].meta().cols()))
      continue
    for r in range(adc_v1[p].meta().rows()):
      for c in range(adc_v1[p].meta().cols()):
        if adc_v1[p].pixel(r,c) != adc_v2[p].pixel(r,c):
          print("MISMATCH IN PIXEL VALUE IN %s IMAGE FOR PLANE %i ROW %i COLUMN %i!!!"%(image_type,p,r,c))
          print("  file1 pixel value: %f"%adc_v1[p].pixel(r,c))
          print("  file2 pixel value: %f"%adc_v2[p].pixel(r,c))

