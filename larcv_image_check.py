
import os,sys,argparse

import ROOT as rt

from larlite import larlite
from larcv import larcv

parser = argparse.ArgumentParser("Pring larcv image pixel values")
parser.add_argument("-i", "--larcv_file", type=str, required=True, help="input larcv images file")
parser.add_argument("-o", "--output_file", type=str, default="larcv_image_check_output.txt", help="output text file with image pixel values")
parser.add_argument("-e", "--entry", type=int, default=0, help="entry in input file to print pixel values for")
parser.add_argument("--tickforward", help="read in larcv images with tickforward option", action="store_true")
args = parser.parse_args()

if args.tickforward:
  iolcv = larcv.IOManager(larcv.IOManager.kREAD, "larcv", larcv.IOManager.kTickForward)
  iolcv.add_in_file(args.larcv_file)
  iolcv.initialize()
else:
  iolcv = larcv.IOManager(larcv.IOManager.kREAD, "larcv", larcv.IOManager.kTickBackward)
  iolcv.add_in_file(args.larcv_file)
  iolcv.reverse_all_products()
  iolcv.initialize()

ioll = larlite.storage_manager(larlite.storage_manager.kREAD)
ioll.add_in_filename(args.larcv_file)
ioll.open()

iolcv.read_entry(args.entry)
ioll.go_to(args.entry)

images = ["wire","thrumu","ancestor","instance","larflow","segment"]

with open(args.output_file, "w") as output:
  output.write("Image Values for Run %i Subrun %i Event %i\n"%(ioll.run_id(),ioll.subrun_id(),ioll.event_id()))
  for image_type in images:
    output.write("plane row col pixel_value for %s image:\n"%image_type)
    adc_v = iolcv.get_data(larcv.kProductImage2D, image_type).Image2DArray()
    for p in range(adc_v.size()):
      for r in range(adc_v[p].meta().rows()):
        for c in range(adc_v[p].meta().cols()):
          output.write("  %i %i %i %f\n"%(p, r, c, adc_v[p].pixel(r, c)))

