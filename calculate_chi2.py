
import ROOT as rt
from math import sqrt
import numpy as np
from scipy.stats.distributions import chi2

rt.TH1.SetDefaultSumw2(rt.kTRUE)

#look at print_xsec_flux_reint_systematics.py to get flux+xsec+reint covariance matrix
#look at combine_detvar_matrices.py to get detvar covariance matrix
#look at make_gen2_paper_systematics_plots.py or print_stats_errors.py for diag. stats covariance matrix
#look at print_stats_errors.py to get bin counts for prediction and converting fraction covar to full covar
#chi2_CNP documented in https://arxiv.org/abs/1903.07185
#  add diagonal 3/(1/M_i + 2/P_i) covar matrix to V_stats + V_sys for full covariance matrix V
#  then chi2_CNP = (M - P)^T V^-1 (M - P)


#========= Get Data and MC predictions as TH1F ===========================================

#f_CCnue = rt.TFile("selection_output/plot_selection_test_results_output/plot_selection_test_results_run1MC_withProdProcessAndMuCuts_noCosThetaOrVtxDistCuts_output.root")
#f_CCnumu = rt.TFile("selection_output/plot_selection_test_results_output/plot_muon_selection_test_results_output_justNMuAndWCCosCuts.root")
#f_CCnue = rt.TFile("plot_selection_test_results_output.root")
f_CCnue = rt.TFile("plot_selection_test_results_chi2_output.root")
f_CCnumu = rt.TFile("plot_muon_selection_test_results_output.root")

hData_CCnue_nuE = f_CCnue.Get("h_visE_data_wCuts")
hData_CCnumu_nuE = f_CCnumu.Get("h_visE_data_wCuts")
hData_CCnue_lepP = f_CCnue.Get("h_lepP_data_wCuts")
hData_CCnumu_lepP = f_CCnumu.Get("h_lepP_data_wCuts")
hData_CCnue_cosTheta_raw = f_CCnue.Get("h_cosTheta_data_wCuts")
hData_CCnumu_cosTheta_raw = f_CCnumu.Get("h_cosTheta_data_wCuts")
hData_CCnue_cosTheta = rt.TH1F("hData_CCnue_cosTheta","hData_CCnue_cosTheta",6,0,6)
hData_CCnumu_cosTheta = rt.TH1F("hData_CCnumu_cosTheta","hData_CCnumu_cosTheta",16,0,16)

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
hStack_CCnumu_nuE = f_CCnumu.Get("h_visE_all_wCuts")
hStack_CCnue_lepP = f_CCnue.Get("h_lepP_all_wCuts")
hStack_CCnumu_lepP = f_CCnumu.Get("h_lepP_all_wCuts")
hStack_CCnue_cosTheta = f_CCnue.Get("h_cosTheta_all_wCuts")
hStack_CCnumu_cosTheta = f_CCnumu.Get("h_cosTheta_all_wCuts")

hList_CCnue_nuE = hStack_CCnue_nuE.GetHists()
hList_CCnumu_nuE = hStack_CCnumu_nuE.GetHists()
hList_CCnue_lepP = hStack_CCnue_lepP.GetHists()
hList_CCnumu_lepP = hStack_CCnumu_lepP.GetHists()
hList_CCnue_cosTheta = hStack_CCnue_cosTheta.GetHists()
hList_CCnumu_cosTheta = hStack_CCnumu_cosTheta.GetHists()

hPred_CCnue_nuE = getPred(hList_CCnue_nuE)
hPred_CCnumu_nuE = getPred(hList_CCnumu_nuE)
hPred_CCnue_lepP = getPred(hList_CCnue_lepP)
hPred_CCnumu_lepP = getPred(hList_CCnumu_lepP)
hPred_CCnue_cosTheta_raw = getPred(hList_CCnue_cosTheta)
hPred_CCnumu_cosTheta_raw = getPred(hList_CCnumu_cosTheta)
hPred_CCnue_cosTheta = rt.TH1F("hPred_CCnue_cosTheta","hPred_CCnue_cosTheta",6,0,6)
hPred_CCnumu_cosTheta = rt.TH1F("hPred_CCnumu_cosTheta","hPred_CCnumu_cosTheta",16,0,16)

for i in range(1,17,1):
  if i <= 6:
    hPred_CCnue_cosTheta.SetBinContent(i, hPred_CCnue_cosTheta_raw.GetBinContent(i+1))
    hPred_CCnue_cosTheta.SetBinError(i, hPred_CCnue_cosTheta_raw.GetBinError(i+1))
    hData_CCnue_cosTheta.SetBinContent(i, hData_CCnue_cosTheta_raw.GetBinContent(i+1))
    hData_CCnue_cosTheta.SetBinError(i, hData_CCnue_cosTheta_raw.GetBinError(i+1))
  hPred_CCnumu_cosTheta.SetBinContent(i, hPred_CCnumu_cosTheta_raw.GetBinContent(i+1))
  hPred_CCnumu_cosTheta.SetBinError(i, hPred_CCnumu_cosTheta_raw.GetBinError(i+1))
  hData_CCnumu_cosTheta.SetBinContent(i, hData_CCnumu_cosTheta_raw.GetBinContent(i+1))
  hData_CCnumu_cosTheta.SetBinError(i, hData_CCnumu_cosTheta_raw.GetBinError(i+1))

#-----------------------------------------------------------------------------------------



#============== Make MCstat and CNP Covar as TH2F =======================================

h_stat_CCnue_nuE = rt.TH2F("h_stat_CCnue_nuE","Fractional Covariance Matrix",9,0,9,9,0,9)
h_stat_CCnumu_nuE = rt.TH2F("h_stat_CCnumu_nuE","Fractional Covariance Matrix",21,0,21,21,0,21)
h_stat_CCnue_lepP = rt.TH2F("h_stat_CCnue_lepP","Fractional Covariance Matrix",8,0,8,8,0,8)
h_stat_CCnumu_lepP = rt.TH2F("h_stat_CCnumu_lepP","Fractional Covariance Matrix",21,0,21,21,0,21)
h_stat_CCnue_cosTheta = rt.TH2F("h_stat_CCnue_cosTheta","Fractional Covariance Matrix",6,0,6,6,0,6)
h_stat_CCnumu_cosTheta = rt.TH2F("h_stat_CCnumu_cosTheta","Fractional Covariance Matrix",16,0,16,16,0,16)

h_cnp_CCnue_nuE = rt.TH2F("h_cnp_CCnue_nuE","Fractional Covariance Matrix",9,0,9,9,0,9)
h_cnp_CCnumu_nuE = rt.TH2F("h_cnp_CCnumu_nuE","Fractional Covariance Matrix",21,0,21,21,0,21)
h_cnp_CCnue_lepP = rt.TH2F("h_cnp_CCnue_lepP","Fractional Covariance Matrix",8,0,8,8,0,8)
h_cnp_CCnumu_lepP = rt.TH2F("h_cnp_CCnumu_lepP","Fractional Covariance Matrix",21,0,21,21,0,21)
h_cnp_CCnue_cosTheta = rt.TH2F("h_cnp_CCnue_cosTheta","Fractional Covariance Matrix",6,0,6,6,0,6)
h_cnp_CCnumu_cosTheta = rt.TH2F("h_cnp_CCnumu_cosTheta","Fractional Covariance Matrix",16,0,16,16,0,16)

for i in range(1,22,1):
  for j in range(1,22,1):
    if i <= 6 and j <= 6:
      if i != j:
        h_stat_CCnue_cosTheta.SetBinContent(i,j,0.)
        h_cnp_CCnue_cosTheta.SetBinContent(i,j,0.)
      else:
        cnp_CCnue_cosTheta = 3.0/((1.0/hData_CCnue_cosTheta.GetBinContent(i)) + (2.0/hPred_CCnue_cosTheta.GetBinContent(i))) if hData_CCnue_cosTheta.GetBinContent(i) > 0. else hPred_CCnue_cosTheta.GetBinContent(i)/2.
        h_stat_CCnue_cosTheta.SetBinContent(i,i, hPred_CCnue_cosTheta.GetBinError(i))
        h_cnp_CCnue_cosTheta.SetBinContent(i,i, cnp_CCnue_cosTheta)
        print(i, hData_CCnue_cosTheta.GetBinContent(i), hPred_CCnue_cosTheta.GetBinContent(i), cnp_CCnue_cosTheta)
    if i <= 8 and j <= 8:
      if i != j:
        h_stat_CCnue_lepP.SetBinContent(i,j,0.)
        h_cnp_CCnue_lepP.SetBinContent(i,j,0.)
      else:
        cnp_lepP = 3.0/((1.0/hData_CCnue_lepP.GetBinContent(i)) + (2.0/hPred_CCnue_lepP.GetBinContent(i))) if hData_CCnue_lepP.GetBinContent(i) > 0. else hPred_CCnue_lepP.GetBinContent(i)/2.
        h_stat_CCnue_lepP.SetBinContent(i,i, hPred_CCnue_lepP.GetBinError(i))
        h_cnp_CCnue_lepP.SetBinContent(i,i, cnp_lepP)
    if i <= 9 and j <= 9:
      if i != j:
        h_stat_CCnue_nuE.SetBinContent(i,j,0.)
        h_cnp_CCnue_nuE.SetBinContent(i,j,0.)
      else:
        cnp_nuE = 3.0/((1.0/hData_CCnue_nuE.GetBinContent(i)) + (2.0/hPred_CCnue_nuE.GetBinContent(i))) if hData_CCnue_nuE.GetBinContent(i) > 0. else hPred_CCnue_nuE.GetBinContent(i)/2.
        h_stat_CCnue_nuE.SetBinContent(i,i, hPred_CCnue_nuE.GetBinError(i))
        h_cnp_CCnue_nuE.SetBinContent(i,i, cnp_nuE)
    if i <= 16 and j <= 16:
      if i != j:
        h_stat_CCnumu_cosTheta.SetBinContent(i,j,0.)
        h_cnp_CCnumu_cosTheta.SetBinContent(i,j,0.)
      else:
        cnp_CCnumu_cosTheta = 3.0/((1.0/hData_CCnumu_cosTheta.GetBinContent(i)) + (2.0/hPred_CCnumu_cosTheta.GetBinContent(i))) if hData_CCnumu_cosTheta.GetBinContent(i) > 0. else hPred_CCnumu_cosTheta.GetBinContent(i)/2.
        h_stat_CCnumu_cosTheta.SetBinContent(i,i, hPred_CCnumu_cosTheta.GetBinError(i))
        h_cnp_CCnumu_cosTheta.SetBinContent(i,i, cnp_CCnumu_cosTheta)
    if i != j:
      h_stat_CCnumu_nuE.SetBinContent(i,j,0.)
      h_cnp_CCnumu_nuE.SetBinContent(i,j,0.)
      h_stat_CCnumu_lepP.SetBinContent(i,j,0.)
      h_cnp_CCnumu_lepP.SetBinContent(i,j,0.)
    else:
      cnp_nuE = 3.0/((1.0/hData_CCnumu_nuE.GetBinContent(i)) + (2.0/hPred_CCnumu_nuE.GetBinContent(i))) if hData_CCnumu_nuE.GetBinContent(i) > 0. else hPred_CCnumu_nuE.GetBinContent(i)/2.
      cnp_lepP = 3.0/((1.0/hData_CCnumu_lepP.GetBinContent(i)) + (2.0/hPred_CCnumu_lepP.GetBinContent(i))) if hData_CCnumu_lepP.GetBinContent(i) > 0. else hPred_CCnumu_lepP.GetBinContent(i)/2.
      h_stat_CCnumu_nuE.SetBinContent(i,i, hPred_CCnumu_nuE.GetBinError(i))
      h_cnp_CCnumu_nuE.SetBinContent(i,i, cnp_nuE)
      h_stat_CCnumu_lepP.SetBinContent(i,i, hPred_CCnumu_lepP.GetBinError(i))
      h_cnp_CCnumu_lepP.SetBinContent(i,i, cnp_lepP)

#--------------------------------------------------------------------------------------------------




# =========== Get flux+xsec+reint covar as TH2F ===================================

#f_sbn_CCnue_nuE = rt.TFile("systematics/CCnueInc_nuE/SBNfit_covariance_plots_CCnueInc.root")
f_sbn_CCnue_nuE = rt.TFile("systematics/CCnueInc_nuE_chi2/SBNfit_covariance_plots_CCnueInc_nuE_chi2.root")
#f_sbn_CCnue_lepP = rt.TFile("systematics/CCnueInc_lepP/SBNfit_covariance_plots_CCnueInc_lepP.root")
f_sbn_CCnue_lepP = rt.TFile("systematics/CCnueInc_lepP_chi2/SBNfit_covariance_plots_CCnueInc_lepP_chi2.root")
#f_sbn_CCnue_cosTheta = rt.TFile("systematics/CCnueInc_cosTheta/SBNfit_covariance_plots_CCnueInc_cosTheta.root")
f_sbn_CCnue_cosTheta = rt.TFile("systematics/CCnueInc_cosTheta_chi2/SBNfit_covariance_plots_CCnueInc_cosTheta_chi2.root")
f_sbn_CCnumu_nuE = rt.TFile("systematics/CCnumuInc_nuE/SBNfit_covariance_plots_CCnumuInc.root")
f_sbn_CCnumu_lepP = rt.TFile("systematics/CCnumuInc_lepP/SBNfit_covariance_plots_CCnumuInc_lepP.root")
f_sbn_CCnumu_cosTheta = rt.TFile("systematics/CCnumuInc_cosTheta/SBNfit_covariance_plots_CCnumuInc_cosTheta.root")

h_sbn_CCnue_nuE = f_sbn_CCnue_nuE.Get("coll_frac")
h_sbn_CCnue_lepP = f_sbn_CCnue_lepP.Get("coll_frac")
h_sbn_CCnue_cosTheta = f_sbn_CCnue_cosTheta.Get("coll_frac")
h_sbn_CCnumu_nuE = f_sbn_CCnumu_nuE.Get("coll_frac")
h_sbn_CCnumu_lepP = f_sbn_CCnumu_lepP.Get("coll_frac")
h_sbn_CCnumu_cosTheta = f_sbn_CCnumu_cosTheta.Get("coll_frac")

for i in range(1,22,1):
  for j in range(1,22,1):
    if i <= 6 and j <= 6:
      vij_CCnue_cosTheta = h_sbn_CCnue_cosTheta.GetBinContent(i,j)*hPred_CCnue_cosTheta.GetBinContent(i)*hPred_CCnue_cosTheta.GetBinContent(j)
      h_sbn_CCnue_cosTheta.SetBinContent(i,j, vij_CCnue_cosTheta)
    if i <= 8 and j <= 8:
      vij_lepP = h_sbn_CCnue_lepP.GetBinContent(i,j)*hPred_CCnue_lepP.GetBinContent(i)*hPred_CCnue_lepP.GetBinContent(j)
      h_sbn_CCnue_lepP.SetBinContent(i,j, vij_lepP)
    if i <= 9 and j <= 9:
      vij_nuE = h_sbn_CCnue_nuE.GetBinContent(i,j)*hPred_CCnue_nuE.GetBinContent(i)*hPred_CCnue_nuE.GetBinContent(j)
      h_sbn_CCnue_nuE.SetBinContent(i,j, vij_nuE)
    if i <= 16 and j <= 16:
      vij_CCnumu_cosTheta = h_sbn_CCnumu_cosTheta.GetBinContent(i,j)*hPred_CCnumu_cosTheta.GetBinContent(i)*hPred_CCnumu_cosTheta.GetBinContent(j)
      h_sbn_CCnumu_cosTheta.SetBinContent(i,j, vij_CCnumu_cosTheta)
    vij_nuE = h_sbn_CCnumu_nuE.GetBinContent(i,j)*hPred_CCnumu_nuE.GetBinContent(i)*hPred_CCnumu_nuE.GetBinContent(j)
    vij_lepP = h_sbn_CCnumu_lepP.GetBinContent(i,j)*hPred_CCnumu_lepP.GetBinContent(i)*hPred_CCnumu_lepP.GetBinContent(j)
    h_sbn_CCnumu_nuE.SetBinContent(i,j, vij_nuE)
    h_sbn_CCnumu_lepP.SetBinContent(i,j, vij_lepP)

#-----------------------------------------------------------------------------------


#========= Get detvar covar as TH2F ===============================================

h_det_CCnue_nuE = rt.TH2F("h_det_CCnue_nuE","Fractional Covariance Matrix",9,0,9,9,0,9)
h_det_CCnumu_nuE = rt.TH2F("h_det_CCnumu_nuE","Fractional Covariance Matrix",21,0,21,21,0,21)
h_det_CCnumu_nuE_raw = rt.TH2F("h_det_CCnumu_nuE_raw","Fractional Covariance Matrix",15,0,15,15,0,15)
h_det_CCnue_lepP = rt.TH2F("h_det_CCnue_lepP","Fractional Covariance Matrix",8,0,8,8,0,8)
h_det_CCnumu_lepP = rt.TH2F("h_det_CCnumu_lepP","Fractional Covariance Matrix",21,0,21,21,0,21)
h_det_CCnumu_lepP_raw = rt.TH2F("h_det_CCnumu_lepP_raw","Fractional Covariance Matrix",13,0,13,13,0,13)
h_det_CCnue_cosTheta = rt.TH2F("h_det_CCnue_cosTheta","Fractional Covariance Matrix",6,0,6,6,0,6)
h_det_CCnumu_cosTheta = rt.TH2F("h_det_CCnumu_cosTheta","Fractional Covariance Matrix",16,0,16,16,0,16)

detvars = ["LYAtt-bugFix", "LYDown", "LYRayleigh", "recomb2", "SCE", "wiremodThetaXZ", "wiremodThetaYZ", "wiremodX", "wiremodYZ"]
for detvar in detvars:
  f = rt.TFile("systematics/detvar/calculate_detvar_systematics_%s_bkgHack_chi2_output.root"%detvar, "READ")
  CCnue_nuE = f.Get("CCnue_frac_covar_recoNuE")
  CCnumu_nuE = f.Get("CCnumu_frac_covar_recoNuE")
  CCnue_lepP = f.Get("CCnue_frac_covar_recoLepP")
  CCnumu_lepP = f.Get("CCnumu_frac_covar_recoLepP")
  CCnue_cosTheta = f.Get("CCnue_frac_covar_recoLepCosTheta")
  CCnumu_cosTheta = f.Get("CCnumu_frac_covar_recoLepCosTheta")
  h_det_CCnue_nuE.Add(CCnue_nuE)
  h_det_CCnumu_nuE_raw.Add(CCnumu_nuE)
  h_det_CCnue_lepP.Add(CCnue_lepP)
  h_det_CCnumu_lepP_raw.Add(CCnumu_lepP)
  h_det_CCnue_cosTheta.Add(CCnue_cosTheta)
  h_det_CCnumu_cosTheta.Add(CCnumu_cosTheta)

for i in range(1,22,1):
  for j in range(1,22,1):
    if i <= 15 and j <= 15:
      h_det_CCnumu_nuE.SetBinContent(i,j, h_det_CCnumu_nuE_raw.GetBinContent(i,j))
    elif j > 15 and i <= 15:
      h_det_CCnumu_nuE.SetBinContent(i,j, h_det_CCnumu_nuE_raw.GetBinContent(i,15))
    elif i > 15 and j <= 15:
      h_det_CCnumu_nuE.SetBinContent(i,j, h_det_CCnumu_nuE_raw.GetBinContent(15,j))
    elif i > 15 and j > 15:
      h_det_CCnumu_nuE.SetBinContent(i,j, h_det_CCnumu_nuE_raw.GetBinContent(15,15))
    else:
      print("ERROR IN FILLING EMPTY CCNUMU NU ENERGY MATRIX ENTRIES!!!")
    if i <= 13 and j <= 13:
      h_det_CCnumu_lepP.SetBinContent(i,j, h_det_CCnumu_lepP_raw.GetBinContent(i,j))
    elif j > 13 and i <= 13:
      h_det_CCnumu_lepP.SetBinContent(i,j, h_det_CCnumu_lepP_raw.GetBinContent(i,13))
    elif i > 13 and j <= 13:
      h_det_CCnumu_lepP.SetBinContent(i,j, h_det_CCnumu_lepP_raw.GetBinContent(13,j))
    elif i > 13 and j > 13:
      h_det_CCnumu_lepP.SetBinContent(i,j, h_det_CCnumu_lepP_raw.GetBinContent(13,13))
    else:
      print("ERROR IN FILLING EMPTY CCNUMU MU MOMENTUM MATRIX ENTRIES!!!")

for i in range(1,22,1):
  for j in range(1,22,1):
    if i <= 6 and j <= 6:
      vij_CCnue_cosTheta = h_det_CCnue_cosTheta.GetBinContent(i,j)*hPred_CCnue_cosTheta.GetBinContent(i)*hPred_CCnue_cosTheta.GetBinContent(j)
      h_det_CCnue_cosTheta.SetBinContent(i,j, vij_CCnue_cosTheta)
    if i <= 8 and j <= 8:
      vij_lepP = h_det_CCnue_lepP.GetBinContent(i,j)*hPred_CCnue_lepP.GetBinContent(i)*hPred_CCnue_lepP.GetBinContent(j)
      h_det_CCnue_lepP.SetBinContent(i,j, vij_lepP)
    if i <= 9 and j <= 9:
      vij_nuE = h_det_CCnue_nuE.GetBinContent(i,j)*hPred_CCnue_nuE.GetBinContent(i)*hPred_CCnue_nuE.GetBinContent(j)
      h_det_CCnue_nuE.SetBinContent(i,j, vij_nuE)
    if i <= 16 and j <= 16:
      vij_CCnumu_cosTheta = h_det_CCnumu_cosTheta.GetBinContent(i,j)*hPred_CCnumu_cosTheta.GetBinContent(i)*hPred_CCnumu_cosTheta.GetBinContent(j)
      h_det_CCnumu_cosTheta.SetBinContent(i,j, vij_CCnumu_cosTheta)
    vij_nuE = h_det_CCnumu_nuE.GetBinContent(i,j)*hPred_CCnumu_nuE.GetBinContent(i)*hPred_CCnumu_nuE.GetBinContent(j)
    vij_lepP = h_det_CCnumu_lepP.GetBinContent(i,j)*hPred_CCnumu_lepP.GetBinContent(i)*hPred_CCnumu_lepP.GetBinContent(j)
    h_det_CCnumu_nuE.SetBinContent(i,j, vij_nuE)
    h_det_CCnumu_lepP.SetBinContent(i,j, vij_lepP)

#-----------------------------------------------------------------------------------------



#============= write full fractional covariance matrices to file ===========================

h_fullFrac_CCnue_nuE = rt.TH2F("h_fullFrac_CCnue_nuE","Fractional Covariance Matrix for CC #nu_{e} Selection Binned in E_{#nu}",9,0,9,9,0,9)
h_fullFrac_CCnumu_nuE = rt.TH2F("h_fullFrac_CCnumu_nuE","Fractional Covariance Matrix for CC #nu_{#mu} Selection Binned in E_{#nu}",21,0,21,21,0,21)
h_fullFrac_CCnue_lepP = rt.TH2F("h_fullFrac_CCnue_lepP","Fractional Covariance Matrix for CC #nu_{e} Selection Binned in p_{e^{-}}",8,0,8,8,0,8)
h_fullFrac_CCnumu_lepP = rt.TH2F("h_fullFrac_CCnumu_lepP","Fractional Covariance Matrix for CC #nu_{#mu} Selection Binned in p_{#mu}",21,0,21,21,0,21)
h_fullFrac_CCnue_cosTheta = rt.TH2F("h_fullFrac_CCnue_cosTheta","Fractional Covariance Matrix for CC #nu_{e} Selection Binned in cos(#theta)",6,0,6,6,0,6)
h_fullFrac_CCnumu_cosTheta = rt.TH2F("h_fullFrac_CCnumu_cosTheta","Fractional Covariance Matrix for CC #nu_{#mu} Selection Binned in cos(#theta)",16,0,16,16,0,16)
h_fullFrac_CCnue_nuE.GetXaxis().SetTitle("Reco Bin i")
h_fullFrac_CCnue_nuE.GetYaxis().SetTitle("Reco Bin j")
h_fullFrac_CCnue_lepP.GetXaxis().SetTitle("Reco Bin i")
h_fullFrac_CCnue_lepP.GetYaxis().SetTitle("Reco Bin j")
h_fullFrac_CCnue_cosTheta.GetXaxis().SetTitle("Reco Bin i")
h_fullFrac_CCnue_cosTheta.GetYaxis().SetTitle("Reco Bin j")
h_fullFrac_CCnumu_nuE.GetXaxis().SetTitle("Reco Bin i")
h_fullFrac_CCnumu_nuE.GetYaxis().SetTitle("Reco Bin j")
h_fullFrac_CCnumu_lepP.GetXaxis().SetTitle("Reco Bin i")
h_fullFrac_CCnumu_lepP.GetYaxis().SetTitle("Reco Bin j")
h_fullFrac_CCnumu_cosTheta.GetXaxis().SetTitle("Reco Bin i")
h_fullFrac_CCnumu_cosTheta.GetYaxis().SetTitle("Reco Bin j")

h_fullFrac_CCnue_nuE.Add(h_cnp_CCnue_nuE)
h_fullFrac_CCnue_nuE.Add(h_stat_CCnue_nuE)
h_fullFrac_CCnue_nuE.Add(h_sbn_CCnue_nuE)
h_fullFrac_CCnue_nuE.Add(h_det_CCnue_nuE)
h_fullFrac_CCnue_lepP.Add(h_cnp_CCnue_lepP)
h_fullFrac_CCnue_lepP.Add(h_stat_CCnue_lepP)
h_fullFrac_CCnue_lepP.Add(h_sbn_CCnue_lepP)
h_fullFrac_CCnue_lepP.Add(h_det_CCnue_lepP)
h_fullFrac_CCnue_cosTheta.Add(h_cnp_CCnue_cosTheta)
h_fullFrac_CCnue_cosTheta.Add(h_stat_CCnue_cosTheta)
h_fullFrac_CCnue_cosTheta.Add(h_sbn_CCnue_cosTheta)
h_fullFrac_CCnue_cosTheta.Add(h_det_CCnue_cosTheta)
h_fullFrac_CCnumu_nuE.Add(h_cnp_CCnumu_nuE)
h_fullFrac_CCnumu_nuE.Add(h_stat_CCnumu_nuE)
h_fullFrac_CCnumu_nuE.Add(h_sbn_CCnumu_nuE)
h_fullFrac_CCnumu_nuE.Add(h_det_CCnumu_nuE)
h_fullFrac_CCnumu_lepP.Add(h_cnp_CCnumu_lepP)
h_fullFrac_CCnumu_lepP.Add(h_stat_CCnumu_lepP)
h_fullFrac_CCnumu_lepP.Add(h_sbn_CCnumu_lepP)
h_fullFrac_CCnumu_lepP.Add(h_det_CCnumu_lepP)
h_fullFrac_CCnumu_cosTheta.Add(h_cnp_CCnumu_cosTheta)
h_fullFrac_CCnumu_cosTheta.Add(h_stat_CCnumu_cosTheta)
h_fullFrac_CCnumu_cosTheta.Add(h_sbn_CCnumu_cosTheta)
h_fullFrac_CCnumu_cosTheta.Add(h_det_CCnumu_cosTheta)

for i in range(1,h_fullFrac_CCnue_nuE.GetNbinsX()+1,1):
  for j in range(1,h_fullFrac_CCnue_nuE.GetNbinsY()+1,1):
    h_fullFrac_CCnue_nuE.SetBinContent(i,j, h_fullFrac_CCnue_nuE.GetBinContent(i,j)/(hPred_CCnue_nuE.GetBinContent(i)*hPred_CCnue_nuE.GetBinContent(j)))
for i in range(1,h_fullFrac_CCnue_lepP.GetNbinsX()+1,1):
  for j in range(1,h_fullFrac_CCnue_lepP.GetNbinsY()+1,1):
    h_fullFrac_CCnue_lepP.SetBinContent(i,j, h_fullFrac_CCnue_lepP.GetBinContent(i,j)/(hPred_CCnue_lepP.GetBinContent(i)*hPred_CCnue_lepP.GetBinContent(j)))
for i in range(1,h_fullFrac_CCnue_cosTheta.GetNbinsX()+1,1):
  for j in range(1,h_fullFrac_CCnue_cosTheta.GetNbinsY()+1,1):
    h_fullFrac_CCnue_cosTheta.SetBinContent(i,j, h_fullFrac_CCnue_cosTheta.GetBinContent(i,j)/(hPred_CCnue_cosTheta.GetBinContent(i)*hPred_CCnue_cosTheta.GetBinContent(j)))
for i in range(1,h_fullFrac_CCnumu_nuE.GetNbinsX()+1,1):
  for j in range(1,h_fullFrac_CCnumu_nuE.GetNbinsY()+1,1):
    h_fullFrac_CCnumu_nuE.SetBinContent(i,j, h_fullFrac_CCnumu_nuE.GetBinContent(i,j)/(hPred_CCnumu_nuE.GetBinContent(i)*hPred_CCnumu_nuE.GetBinContent(j)))
for i in range(1,h_fullFrac_CCnumu_lepP.GetNbinsX()+1,1):
  for j in range(1,h_fullFrac_CCnumu_lepP.GetNbinsY()+1,1):
    h_fullFrac_CCnumu_lepP.SetBinContent(i,j, h_fullFrac_CCnumu_lepP.GetBinContent(i,j)/(hPred_CCnumu_lepP.GetBinContent(i)*hPred_CCnumu_lepP.GetBinContent(j)))
for i in range(1,h_fullFrac_CCnumu_cosTheta.GetNbinsX()+1,1):
  for j in range(1,h_fullFrac_CCnumu_cosTheta.GetNbinsY()+1,1):
    h_fullFrac_CCnumu_cosTheta.SetBinContent(i,j, h_fullFrac_CCnumu_cosTheta.GetBinContent(i,j)/(hPred_CCnumu_cosTheta.GetBinContent(i)*hPred_CCnumu_cosTheta.GetBinContent(j)))

outfile = rt.TFile("calculate_chi2_full_coviance_plots.root", "RECREATE")
rt.gStyle.SetOptStat(0)

cnv_CCnue_nuE = rt.TCanvas("cnv_CCnue_nuE","cnv_CCnue_nuE")
h_fullFrac_CCnue_nuE.Draw("COLZ")
cnv_CCnue_nuE.Write()
h_fullFrac_CCnue_nuE.Write()
cnv_CCnumu_nuE = rt.TCanvas("cnv_CCnumu_nuE","cnv_CCnumu_nuE")
h_fullFrac_CCnumu_nuE.Draw("COLZ")
cnv_CCnumu_nuE.Write()
h_fullFrac_CCnumu_nuE.Write()
cnv_CCnue_lepP = rt.TCanvas("cnv_CCnue_lepP","cnv_CCnue_lepP")
h_fullFrac_CCnue_lepP.Draw("COLZ")
cnv_CCnue_lepP.Write()
h_fullFrac_CCnue_lepP.Write()
cnv_CCnumu_lepP = rt.TCanvas("cnv_CCnumu_lepP","cnv_CCnumu_lepP")
h_fullFrac_CCnumu_lepP.Draw("COLZ")
cnv_CCnumu_lepP.Write()
h_fullFrac_CCnumu_lepP.Write()
cnv_CCnue_cosTheta = rt.TCanvas("cnv_CCnue_cosTheta","cnv_CCnue_cosTheta")
h_fullFrac_CCnue_cosTheta.Draw("COLZ")
cnv_CCnue_cosTheta.Write()
h_fullFrac_CCnue_cosTheta.Write()
cnv_CCnumu_cosTheta = rt.TCanvas("cnv_CCnumu_cosTheta","cnv_CCnumu_cosTheta")
h_fullFrac_CCnumu_cosTheta.Draw("COLZ")
cnv_CCnumu_cosTheta.Write()
h_fullFrac_CCnumu_cosTheta.Write()

#-------------------------------------------------------------------------------------------



#============= convert root objects into numpy matrices and vectors ========================

vPred_CCnue_nuE = np.array([hPred_CCnue_nuE.GetBinContent(i) for i in range(1,hPred_CCnue_nuE.GetNbinsX()+1,1)])
vPred_CCnue_lepP = np.array([hPred_CCnue_lepP.GetBinContent(i) for i in range(1,hPred_CCnue_lepP.GetNbinsX()+1,1)])
vPred_CCnue_cosTheta = np.array([hPred_CCnue_cosTheta.GetBinContent(i) for i in range(1,hPred_CCnue_cosTheta.GetNbinsX()+1,1)])
vData_CCnue_nuE = np.array([hData_CCnue_nuE.GetBinContent(i) for i in range(1,hData_CCnue_nuE.GetNbinsX()+1,1)])
vData_CCnue_lepP = np.array([hData_CCnue_lepP.GetBinContent(i) for i in range(1,hData_CCnue_lepP.GetNbinsX()+1,1)])
vData_CCnue_cosTheta = np.array([hData_CCnue_cosTheta.GetBinContent(i) for i in range(1,hData_CCnue_cosTheta.GetNbinsX()+1,1)])
vPred_CCnumu_nuE = np.array([hPred_CCnumu_nuE.GetBinContent(i) for i in range(1,hPred_CCnumu_nuE.GetNbinsX()+1,1)])
vPred_CCnumu_lepP = np.array([hPred_CCnumu_lepP.GetBinContent(i) for i in range(1,hPred_CCnumu_lepP.GetNbinsX()+1,1)])
vPred_CCnumu_cosTheta = np.array([hPred_CCnumu_cosTheta.GetBinContent(i) for i in range(1,hPred_CCnumu_cosTheta.GetNbinsX()+1,1)])
vData_CCnumu_nuE = np.array([hData_CCnumu_nuE.GetBinContent(i) for i in range(1,hData_CCnumu_nuE.GetNbinsX()+1,1)])
vData_CCnumu_lepP = np.array([hData_CCnumu_lepP.GetBinContent(i) for i in range(1,hData_CCnumu_lepP.GetNbinsX()+1,1)])
vData_CCnumu_cosTheta = np.array([hData_CCnumu_cosTheta.GetBinContent(i) for i in range(1,hData_CCnumu_cosTheta.GetNbinsX()+1,1)])

cov_cnp_CCnue_nuE = np.array([ [h_cnp_CCnue_nuE.GetBinContent(i,j) for j in range(1,h_cnp_CCnue_nuE.GetNbinsY()+1,1)] for i in range(1,h_cnp_CCnue_nuE.GetNbinsX()+1,1) ])
cov_cnp_CCnue_lepP = np.array([ [h_cnp_CCnue_lepP.GetBinContent(i,j) for j in range(1,h_cnp_CCnue_lepP.GetNbinsY()+1,1)] for i in range(1,h_cnp_CCnue_lepP.GetNbinsX()+1,1) ])
cov_cnp_CCnue_cosTheta = np.array([ [h_cnp_CCnue_cosTheta.GetBinContent(i,j) for j in range(1,h_cnp_CCnue_cosTheta.GetNbinsY()+1,1)] for i in range(1,h_cnp_CCnue_cosTheta.GetNbinsX()+1,1) ])
cov_cnp_CCnumu_nuE = np.array([ [h_cnp_CCnumu_nuE.GetBinContent(i,j) for j in range(1,h_cnp_CCnumu_nuE.GetNbinsY()+1,1)] for i in range(1,h_cnp_CCnumu_nuE.GetNbinsX()+1,1) ])
cov_cnp_CCnumu_lepP = np.array([ [h_cnp_CCnumu_lepP.GetBinContent(i,j) for j in range(1,h_cnp_CCnumu_lepP.GetNbinsY()+1,1)] for i in range(1,h_cnp_CCnumu_lepP.GetNbinsX()+1,1) ])
cov_cnp_CCnumu_cosTheta = np.array([ [h_cnp_CCnumu_cosTheta.GetBinContent(i,j) for j in range(1,h_cnp_CCnumu_cosTheta.GetNbinsY()+1,1)] for i in range(1,h_cnp_CCnumu_cosTheta.GetNbinsX()+1,1) ])

cov_stat_CCnue_nuE = np.array([ [h_stat_CCnue_nuE.GetBinContent(i,j) for j in range(1,h_stat_CCnue_nuE.GetNbinsY()+1,1)] for i in range(1,h_stat_CCnue_nuE.GetNbinsX()+1,1) ])
cov_stat_CCnue_lepP = np.array([ [h_stat_CCnue_lepP.GetBinContent(i,j) for j in range(1,h_stat_CCnue_lepP.GetNbinsY()+1,1)] for i in range(1,h_stat_CCnue_lepP.GetNbinsX()+1,1) ])
cov_stat_CCnue_cosTheta = np.array([ [h_stat_CCnue_cosTheta.GetBinContent(i,j) for j in range(1,h_stat_CCnue_cosTheta.GetNbinsY()+1,1)] for i in range(1,h_stat_CCnue_cosTheta.GetNbinsX()+1,1) ])
cov_stat_CCnumu_nuE = np.array([ [h_stat_CCnumu_nuE.GetBinContent(i,j) for j in range(1,h_stat_CCnumu_nuE.GetNbinsY()+1,1)] for i in range(1,h_stat_CCnumu_nuE.GetNbinsX()+1,1) ])
cov_stat_CCnumu_lepP = np.array([ [h_stat_CCnumu_lepP.GetBinContent(i,j) for j in range(1,h_stat_CCnumu_lepP.GetNbinsY()+1,1)] for i in range(1,h_stat_CCnumu_lepP.GetNbinsX()+1,1) ])
cov_stat_CCnumu_cosTheta = np.array([ [h_stat_CCnumu_cosTheta.GetBinContent(i,j) for j in range(1,h_stat_CCnumu_cosTheta.GetNbinsY()+1,1)] for i in range(1,h_stat_CCnumu_cosTheta.GetNbinsX()+1,1) ])

cov_sbn_CCnue_nuE = np.array([ [h_sbn_CCnue_nuE.GetBinContent(i,j) for j in range(1,h_sbn_CCnue_nuE.GetNbinsY()+1,1)] for i in range(1,h_sbn_CCnue_nuE.GetNbinsX()+1,1) ])
cov_sbn_CCnue_lepP = np.array([ [h_sbn_CCnue_lepP.GetBinContent(i,j) for j in range(1,h_sbn_CCnue_lepP.GetNbinsY()+1,1)] for i in range(1,h_sbn_CCnue_lepP.GetNbinsX()+1,1) ])
cov_sbn_CCnue_cosTheta = np.array([ [h_sbn_CCnue_cosTheta.GetBinContent(i,j) for j in range(1,h_sbn_CCnue_cosTheta.GetNbinsY()+1,1)] for i in range(1,h_sbn_CCnue_cosTheta.GetNbinsX()+1,1) ])
cov_sbn_CCnumu_nuE = np.array([ [h_sbn_CCnumu_nuE.GetBinContent(i,j) for j in range(1,h_sbn_CCnumu_nuE.GetNbinsY()+1,1)] for i in range(1,h_sbn_CCnumu_nuE.GetNbinsX()+1,1) ])
cov_sbn_CCnumu_lepP = np.array([ [h_sbn_CCnumu_lepP.GetBinContent(i,j) for j in range(1,h_sbn_CCnumu_lepP.GetNbinsY()+1,1)] for i in range(1,h_sbn_CCnumu_lepP.GetNbinsX()+1,1) ])
cov_sbn_CCnumu_cosTheta = np.array([ [h_sbn_CCnumu_cosTheta.GetBinContent(i,j) for j in range(1,h_sbn_CCnumu_cosTheta.GetNbinsY()+1,1)] for i in range(1,h_sbn_CCnumu_cosTheta.GetNbinsX()+1,1) ])

cov_det_CCnue_nuE = np.array([ [h_det_CCnue_nuE.GetBinContent(i,j) for j in range(1,h_det_CCnue_nuE.GetNbinsY()+1,1)] for i in range(1,h_det_CCnue_nuE.GetNbinsX()+1,1) ])
cov_det_CCnue_lepP = np.array([ [h_det_CCnue_lepP.GetBinContent(i,j) for j in range(1,h_det_CCnue_lepP.GetNbinsY()+1,1)] for i in range(1,h_det_CCnue_lepP.GetNbinsX()+1,1) ])
cov_det_CCnue_cosTheta = np.array([ [h_det_CCnue_cosTheta.GetBinContent(i,j) for j in range(1,h_det_CCnue_cosTheta.GetNbinsY()+1,1)] for i in range(1,h_det_CCnue_cosTheta.GetNbinsX()+1,1) ])
cov_det_CCnumu_nuE = np.array([ [h_det_CCnumu_nuE.GetBinContent(i,j) for j in range(1,h_det_CCnumu_nuE.GetNbinsY()+1,1)] for i in range(1,h_det_CCnumu_nuE.GetNbinsX()+1,1) ])
cov_det_CCnumu_lepP = np.array([ [h_det_CCnumu_lepP.GetBinContent(i,j) for j in range(1,h_det_CCnumu_lepP.GetNbinsY()+1,1)] for i in range(1,h_det_CCnumu_lepP.GetNbinsX()+1,1) ])
cov_det_CCnumu_cosTheta = np.array([ [h_det_CCnumu_cosTheta.GetBinContent(i,j) for j in range(1,h_det_CCnumu_cosTheta.GetNbinsY()+1,1)] for i in range(1,h_det_CCnumu_cosTheta.GetNbinsX()+1,1) ])

print()
print("cov_cnp_CCnue_nuE")
print(cov_cnp_CCnue_nuE)
print("cov_cnp_CCnue_lepP")
print(cov_cnp_CCnue_lepP)
print("cov_cnp_CCnue_cosTheta")
print(cov_cnp_CCnue_cosTheta)
print("cov_cnp_CCnumu_nuE")
print(cov_cnp_CCnumu_nuE)
print("cov_cnp_CCnumu_lepP")
print(cov_cnp_CCnumu_lepP)
print("cov_cnp_CCnumu_cosTheta")
print(cov_cnp_CCnumu_cosTheta)
print()
print("cov_stat_CCnue_nuE")
print(cov_stat_CCnue_nuE)
print("cov_stat_CCnue_lepP")
print(cov_stat_CCnue_lepP)
print("cov_stat_CCnue_cosTheta")
print(cov_stat_CCnue_cosTheta)
print("cov_stat_CCnumu_nuE")
print(cov_stat_CCnumu_nuE)
print("cov_stat_CCnumu_lepP")
print(cov_stat_CCnumu_lepP)
print("cov_stat_CCnumu_cosTheta")
print(cov_stat_CCnumu_cosTheta)
print()
print("cov_sbn_CCnue_nuE")
print(cov_sbn_CCnue_nuE)
print("cov_sbn_CCnue_lepP")
print(cov_sbn_CCnue_lepP)
print("cov_sbn_CCnue_cosTheta")
print(cov_sbn_CCnue_cosTheta)
print("cov_sbn_CCnumu_nuE")
print(cov_sbn_CCnumu_nuE)
print("cov_sbn_CCnumu_lepP")
print(cov_sbn_CCnumu_lepP)
print("cov_sbn_CCnumu_cosTheta")
print(cov_sbn_CCnumu_cosTheta)
print()
print("cov_det_CCnue_nuE")
print(cov_det_CCnue_nuE)
print("cov_det_CCnue_lepP")
print(cov_det_CCnue_lepP)
print("cov_det_CCnue_cosTheta")
print(cov_det_CCnue_cosTheta)
print("cov_det_CCnumu_nuE")
print(cov_det_CCnumu_nuE)
print("cov_det_CCnumu_lepP")
print(cov_det_CCnumu_lepP)
print("cov_det_CCnumu_cosTheta")
print(cov_det_CCnumu_cosTheta)
print()
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
print("vPred_CCnumu_nuE")
print(vPred_CCnumu_nuE)
print("vPred_CCnumu_lepP")
print(vPred_CCnumu_lepP)
print("vPred_CCnumu_cosTheta")
print(vPred_CCnumu_cosTheta)
print("vData_CCnumu_nuE")
print(vData_CCnumu_nuE)
print("vData_CCnumu_lepP")
print(vData_CCnumu_lepP)
print("vData_CCnumu_cosTheta")
print(vData_CCnumu_cosTheta)
print()

#-------------------------------------------------------------------------------------



#=============== calculate chi2 and p-values ============================================

cov_full_CCnue_nuE = cov_cnp_CCnue_nuE + cov_stat_CCnue_nuE + cov_sbn_CCnue_nuE + cov_det_CCnue_nuE
cov_full_CCnue_lepP = cov_cnp_CCnue_lepP + cov_stat_CCnue_lepP + cov_sbn_CCnue_lepP + cov_det_CCnue_lepP
cov_full_CCnue_cosTheta = cov_cnp_CCnue_cosTheta + cov_stat_CCnue_cosTheta + cov_sbn_CCnue_cosTheta + cov_det_CCnue_cosTheta
cov_full_CCnumu_nuE = cov_cnp_CCnumu_nuE + cov_stat_CCnumu_nuE + cov_sbn_CCnumu_nuE + cov_det_CCnumu_nuE
cov_full_CCnumu_lepP = cov_cnp_CCnumu_lepP + cov_stat_CCnumu_lepP + cov_sbn_CCnumu_lepP + cov_det_CCnumu_lepP
cov_full_CCnumu_cosTheta = cov_cnp_CCnumu_cosTheta + cov_stat_CCnumu_cosTheta + cov_sbn_CCnumu_cosTheta + cov_det_CCnumu_cosTheta

covInv_full_CCnue_nuE = np.linalg.inv(cov_full_CCnue_nuE)
covInv_full_CCnue_lepP = np.linalg.inv(cov_full_CCnue_lepP)
covInv_full_CCnue_cosTheta = np.linalg.inv(cov_full_CCnue_cosTheta)
covInv_full_CCnumu_nuE = np.linalg.inv(cov_full_CCnumu_nuE)
covInv_full_CCnumu_lepP = np.linalg.inv(cov_full_CCnumu_lepP)
covInv_full_CCnumu_cosTheta = np.linalg.inv(cov_full_CCnumu_cosTheta)

diff_CCnue_nuE = vData_CCnue_nuE - vPred_CCnue_nuE
diff_CCnue_lepP = vData_CCnue_lepP - vPred_CCnue_lepP
diff_CCnue_cosTheta = vData_CCnue_cosTheta - vPred_CCnue_cosTheta
diff_CCnumu_nuE = vData_CCnumu_nuE - vPred_CCnumu_nuE
diff_CCnumu_lepP = vData_CCnumu_lepP - vPred_CCnumu_lepP
diff_CCnumu_cosTheta = vData_CCnumu_cosTheta - vPred_CCnumu_cosTheta

chi2_CCnue_nuE = np.dot(diff_CCnue_nuE, np.dot(covInv_full_CCnue_nuE, diff_CCnue_nuE))
chi2_CCnue_lepP = np.dot(diff_CCnue_lepP, np.dot(covInv_full_CCnue_lepP, diff_CCnue_lepP))
chi2_CCnue_cosTheta = np.dot(diff_CCnue_cosTheta, np.dot(covInv_full_CCnue_cosTheta, diff_CCnue_cosTheta))
chi2_CCnumu_nuE = np.dot(diff_CCnumu_nuE, np.dot(covInv_full_CCnumu_nuE, diff_CCnumu_nuE))
chi2_CCnumu_lepP = np.dot(diff_CCnumu_lepP, np.dot(covInv_full_CCnumu_lepP, diff_CCnumu_lepP))
chi2_CCnumu_cosTheta = np.dot(diff_CCnumu_cosTheta, np.dot(covInv_full_CCnumu_cosTheta, diff_CCnumu_cosTheta))

pval_CCnue_nuE = chi2.sf(chi2_CCnue_nuE, h_stat_CCnue_nuE.GetNbinsX())
pval_CCnue_lepP = chi2.sf(chi2_CCnue_lepP, h_stat_CCnue_lepP.GetNbinsX())
pval_CCnue_cosTheta = chi2.sf(chi2_CCnue_cosTheta, h_stat_CCnue_cosTheta.GetNbinsX())
pval_CCnumu_nuE = chi2.sf(chi2_CCnumu_nuE, h_stat_CCnumu_nuE.GetNbinsX())
pval_CCnumu_lepP = chi2.sf(chi2_CCnumu_lepP, h_stat_CCnumu_lepP.GetNbinsX())
pval_CCnumu_cosTheta = chi2.sf(chi2_CCnumu_cosTheta, h_stat_CCnumu_cosTheta.GetNbinsX())

#print("diff_CCnue_nuE")
#print(diff_CCnue_nuE)
#print("cov_full_CCnue_nuE")
#print(cov_full_CCnue_nuE)
#print("covInv_full_CCnue_nuE")
#print(covInv_full_CCnue_nuE)
#print("np.dot(covInv_full_CCnue_nuE,diff_CCnue_nuE)")
#print(np.dot(covInv_full_CCnue_nuE,diff_CCnue_nuE))
#print("np.dot(diff_CCnue_nuE, np.dot(covInv_full_CCnue_nuE, diff_CCnue_nuE))")
#print(np.dot(diff_CCnue_nuE, np.dot(covInv_full_CCnue_nuE, diff_CCnue_nuE)))
#print("diff_CCnue_cosTheta")
#print(diff_CCnue_cosTheta)
#print("cov_full_CCnue_cosTheta")
#print(cov_full_CCnue_cosTheta)
#print("covInv_full_CCnue_cosTheta")
#print(covInv_full_CCnue_cosTheta)
#print("np.dot(covInv_full_CCnue_cosTheta,diff_CCnue_cosTheta)")
#print(np.dot(covInv_full_CCnue_cosTheta,diff_CCnue_cosTheta))
#print("np.dot(diff_CCnue_cosTheta, np.dot(covInv_full_CCnue_cosTheta, diff_CCnue_cosTheta))")
#print(np.dot(diff_CCnue_cosTheta, np.dot(covInv_full_CCnue_cosTheta, diff_CCnue_cosTheta)))

print()
print()
print("CCnue_nuE chi2/ndf, p-val = %f/%i, %f"%(chi2_CCnue_nuE,h_stat_CCnue_nuE.GetNbinsX(),pval_CCnue_nuE))
print("CCnue_lepP chi2/ndf, p-val = %f/%i, %f"%(chi2_CCnue_lepP,h_stat_CCnue_lepP.GetNbinsX(),pval_CCnue_lepP))
print("CCnue_cosTheta chi2/ndf, p-val = %f/%i, %f"%(chi2_CCnue_cosTheta,h_stat_CCnue_cosTheta.GetNbinsX(),pval_CCnue_cosTheta))
print("CCnumu_nuE chi2/ndf, p-val = %f/%i, %f"%(chi2_CCnumu_nuE,h_stat_CCnumu_nuE.GetNbinsX(),pval_CCnumu_nuE))
print("CCnumu_lepP chi2/ndf, p-val = %f/%i, %f"%(chi2_CCnumu_lepP,h_stat_CCnumu_lepP.GetNbinsX(),pval_CCnumu_lepP))
print("CCnumu_cosTheta chi2/ndf, p-val = %f/%i, %f"%(chi2_CCnumu_cosTheta,h_stat_CCnumu_cosTheta.GetNbinsX(),pval_CCnumu_cosTheta))
print()


