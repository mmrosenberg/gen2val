
import ROOT as rt

rt.TH1.SetDefaultSumw2(rt.kTRUE)

#Public note files:
#f_CCnue = rt.TFile("selection_output/plot_selection_test_results_output/plot_selection_test_results_run1MC_withProdProcessAndMuCuts_noCosThetaOrVtxDistCuts_output.root")
#f_CCnumu = rt.TFile("selection_output/plot_selection_test_results_output/plot_muon_selection_test_results_output_justNMuAndWCCosCuts.root")

#Final version files:
f_CCnue = rt.TFile("selection_output/plot_selection_test_results_output/plot_selection_test_results_run1MC_withProdProcessAndMuCuts_noCosThetaOrVtxDistCuts_chi2binning_output.root")
f_CCnumu = rt.TFile("selection_output/plot_selection_test_results_output/plot_muon_selection_test_results_output_justNMuAndWCCosCuts_rebinned.root")

hStack_CCnue = f_CCnue.Get("h_visE_all_wCuts")
hStack_CCnumu = f_CCnumu.Get("h_visE_all_wCuts")

hList_CCnue = hStack_CCnue.GetHists()
hList_CCnumu = hStack_CCnumu.GetHists()

h0_CCnue = hList_CCnue.At(0)
h1_CCnue = hList_CCnue.At(1)
h2_CCnue = hList_CCnue.At(2)
h3_CCnue = hList_CCnue.At(3)
h4_CCnue = hList_CCnue.At(4)
h0_CCnumu = hList_CCnumu.At(0)
h1_CCnumu = hList_CCnumu.At(1)
h2_CCnumu = hList_CCnumu.At(2)
h3_CCnumu = hList_CCnumu.At(3)
h4_CCnumu = hList_CCnumu.At(4)

hPred_CCnue = h0_CCnue.Clone("hPred_CCnue")
hPred_CCnue.Add(h1_CCnue)
hPred_CCnue.Add(h2_CCnue)
hPred_CCnue.Add(h3_CCnue)
hPred_CCnue.Add(h4_CCnue)
hPred_CCnumu = h0_CCnumu.Clone("hPred_CCnumu")
hPred_CCnumu.Add(h1_CCnumu)
hPred_CCnumu.Add(h2_CCnumu)
hPred_CCnumu.Add(h3_CCnumu)
hPred_CCnumu.Add(h4_CCnumu)

fracErr_CCnue = []
countErr_CCnue = []
fracErr_CCnumu = []
countErr_CCnumu = []

for i in range(1,hPred_CCnue.GetNbinsX()+1,1):
  countErr_CCnue.append(hPred_CCnue.GetBinError(i))
  fracErr_CCnue.append(hPred_CCnue.GetBinError(i)/hPred_CCnue.GetBinContent(i))
for i in range(1,hPred_CCnumu.GetNbinsX()+1,1):
  countErr_CCnumu.append(hPred_CCnumu.GetBinError(i))
  fracErr_CCnumu.append(hPred_CCnumu.GetBinError(i)/hPred_CCnumu.GetBinContent(i))

print()
print("CCnue stats error =", countErr_CCnue)
print()
print("CCnue fractional stats error =", fracErr_CCnue)
print()
print("CCnumu stats error =", countErr_CCnumu)
print()
print("CCnumu fractional stats error =", fracErr_CCnumu)
print()


