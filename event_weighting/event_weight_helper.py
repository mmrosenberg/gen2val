
import sys
import pickle
import gc
import ROOT as rt


class Weights:

  def __init__(self, infile):
    f = open(infile, "rb")
    self.weightDict = pickle.load(f)
    f.close()

  def get(self, run, subrun, event):
    return self.weightDict[run][subrun][event]

  def clear(self):
    del self.weightDict
    gc.collect()


def SumPOT(infile):
  mdlfile = rt.TFile(infile)
  pottree = mdlfile.Get("potsummary_generator_tree")
  totPOT = 0.
  totGoodPOT = 0.
  for entry in range(pottree.GetEntries()):
    pottree.GetEntry(entry)
    totPOT = totPOT + pottree.potsummary_generator_branch.totpot
    totGoodPOT = totGoodPOT + pottree.potsummary_generator_branch.totgoodpot
  mdlfile.Close()
  return totPOT, totGoodPOT


def WriteWeights(inrootfile):
  f = rt.TFile(inrootfile)
  t = f.Get("eventweight_tree")
  weights = {}
  for i in range(t.GetEntries()):
    t.GetEntry(i)
    if t.run in weights:
      if t.subrun in weights[t.run]:
        weights[t.run][t.subrun][t.event] = t.xsec_corr_weight
      else:
        weights[t.run][t.subrun] = {t.event: t.xsec_corr_weight}
    else:
      weights[t.run] = {t.subrun: {t.event: t.xsec_corr_weight} }
  outfile = open(inrootfile.replace(".root",".pkl"), "wb")
  pickle.dump(weights, outfile)
  outfile.close()

if __name__ == "__main__":
  #print(SumPOT(sys.argv[1]))
  WriteWeights(sys.argv[1])

