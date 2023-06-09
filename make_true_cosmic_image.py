
import os,sys,argparse

import ROOT as rt

from larcv import larcv

parser = argparse.ArgumentParser("replace merged_dlreco thrumu cosmic image using truth")
parser.add_argument("-i", "--input", type=str, required=True, help="input merged dlreco file")
parser.add_argument("-od", "--output_dir", type=str, default="./", help="output directory")
parser.add_argument("-op", "--output_path", type=str, default="", help="output file path (overrides output_dir option if supplied")
args = parser.parse_args()


if args.input[-5:] != ".root":
  sys.exit("Expected input ending in .root, ensure you are supplying a merged_dlreco file")

larcv_output_name = args.input.replace(".root","_copy1.root")
if args.output_path != "":
  final_output_name = args.output_path
else:
  final_output_name = args.output_dir+"/"+args.input.replace(".root","_trueThrumu.root")


iolcv = larcv.IOManager(larcv.IOManager.kBOTH, "larcv", larcv.IOManager.kTickBackward)
iolcv.add_in_file(args.input)
iolcv.set_out_file(larcv_output_name)
#iolcv.reverse_all_products()
iolcv.initialize()

for i in range(0, iolcv.get_n_entries()):

  iolcv.read_entry(i)
  
  adc_v_thrumu = iolcv.get_data(larcv.kProductImage2D, "thrumu").Image2DArray()
  adc_v_instance = iolcv.get_data(larcv.kProductImage2D, "instance").Image2DArray()

  for p in range(3):    
    for r in range(adc_v_thrumu[p].meta().rows()):
      for c in range(adc_v_thrumu[p].meta().cols()):
        if adc_v_instance[p].pixel(r, c) > 0:
          adc_v_thrumu[p].set_pixel(r, c, 0)

  iolcv.save_entry()

iolcv.finalize()


forig = rt.TFile(args.input)
flarcv = rt.TFile(larcv_output_name)
fout = rt.TFile(final_output_name,"RECREATE")

larcv_tree_list = []

for key in flarcv.GetListOfKeys():
  if key.GetClassName() == "TTree":
    tname = key.GetName()
    tlarcv = flarcv.Get(tname)
    tcopy = tlarcv.CloneTree()
    fout.cd()
    tcopy.Write("",rt.TObject.kOverwrite)
    larcv_tree_list.append(tname)

for key in forig.GetListOfKeys():
  if key.GetClassName() == "TTree":
    tname = key.GetName()
    if tname in larcv_tree_list:
      continue
    torig = forig.Get(tname)
    tcopy = torig.CloneTree()
    fout.cd()
    tcopy.Write("",rt.TObject.kOverwrite)

forig.Close()
flarcv.Close()
fout.Close()

os.system("rm %s"%larcv_output_name)

