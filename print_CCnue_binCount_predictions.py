
import sys, argparse
import ROOT as rt
rt.TH1.SetDefaultSumw2(rt.kTRUE)

parser = argparse.ArgumentParser("Print CCnue Predictions")
parser.add_argument("-i", "--infile", type=str, required=True, help="input file (output of plot_selection_test_results.py)")
args = parser.parse_args()

f_CCnue = rt.TFile(args.infile)

hData_CCnue_nuE = f_CCnue.Get("h_visE_data_wCuts")
hData_CCnue_lepP = f_CCnue.Get("h_lepP_data_wCuts")
hData_CCnue_cosTheta = f_CCnue.Get("h_cosTheta_data_wCuts")

def getPred(hList):
  h0 = hList.At(0)
  h1 = hList.At(1)
  h2 = hList.At(2)
  h3 = hList.At(3)
  h4 = hList.At(4)
  hPred = h0.Clone("hPred")
  hPred.Add(h1)
  hPred.Add(h2)
  hPred.Add(h3)
  hPred.Add(h4)
  return hPred

hStack_CCnue_nuE = f_CCnue.Get("h_visE_all_wCuts")
hStack_CCnue_lepP = f_CCnue.Get("h_lepP_all_wCuts")
hStack_CCnue_cosTheta = f_CCnue.Get("h_cosTheta_all_wCuts")

hList_CCnue_nuE = hStack_CCnue_nuE.GetHists()
hList_CCnue_lepP = hStack_CCnue_lepP.GetHists()
hList_CCnue_cosTheta = hStack_CCnue_cosTheta.GetHists()

hPred_CCnue_nuE = getPred(hList_CCnue_nuE)
hPred_CCnue_lepP = getPred(hList_CCnue_lepP)
hPred_CCnue_cosTheta = getPred(hList_CCnue_cosTheta)

vPred_CCnue_nuE = [hPred_CCnue_nuE.GetBinContent(i) for i in range(1,hPred_CCnue_nuE.GetNbinsX()+1,1)]
vPred_CCnue_lepP = [hPred_CCnue_lepP.GetBinContent(i) for i in range(1,hPred_CCnue_lepP.GetNbinsX()+1,1)]
vPred_CCnue_cosTheta = [hPred_CCnue_cosTheta.GetBinContent(i) for i in range(2,hPred_CCnue_cosTheta.GetNbinsX(),1)]
vData_CCnue_nuE = [hData_CCnue_nuE.GetBinContent(i) for i in range(1,hData_CCnue_nuE.GetNbinsX()+1,1)]
vData_CCnue_lepP = [hData_CCnue_lepP.GetBinContent(i) for i in range(1,hData_CCnue_lepP.GetNbinsX()+1,1)]
vData_CCnue_cosTheta = [hData_CCnue_cosTheta.GetBinContent(i) for i in range(2,hData_CCnue_cosTheta.GetNbinsX(),1)]

print("vPred_CCnue_nuE")
print(vPred_CCnue_nuE)
print("vPred_CCnue_lepP")
print(vPred_CCnue_lepP)
print("vPred_CCnue_cosTheta")
print(vPred_CCnue_cosTheta)
print("vData_CCnue_nuE")
print(vData_CCnue_nuE)
print("vData_CCnue_lepP")
print(vData_CCnue_lepP)
print("vData_CCnue_cosTheta")
print(vData_CCnue_cosTheta)


