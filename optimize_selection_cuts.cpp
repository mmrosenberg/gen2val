
#include "TFile.h"
#include "TTree.h"
#include "TGraph2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TROOT.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <cmath>


std::vector<float> runSelection(float pCut, float sCut, float piCut, float phCut,
                                float procCut, float chgdCut, float ntrlCut, float maxMuCut){

  float n_runs1to3_CCnue = 0.;
  float n_runs1to3_CCnumu_pass = 0.;
  float n_runs1to3_NCnumu_pass = 0.;
  float n_runs1to3_CCnue_pass = 0.;
  float n_runs1to3_NCnue_pass = 0.;
  float n_runs1to3_ext_pass = 0.;


  std::string infiledir = "flat_ntuples/";

  TFile fnu((infiledir+"dlgen2_reco_v2me06_ntuple_v5_mcc9_v28_wctagger_bnboverlay.root").c_str());
  TTree* tnu = (TTree*)fnu.Get("EventTree");
  TTree* tnuPOT = (TTree*)fnu.Get("potTree");

  TFile fnue((infiledir+"dlgen2_reco_v2me06_ntuple_v5_mcc9_v28_wctagger_nueintrinsics.root").c_str());
  TTree* tnue = (TTree*)fnue.Get("EventTree");
  TTree* tnuePOT = (TTree*)fnue.Get("potTree");

  TFile fext((infiledir+"dlgen2_reco_v2me06_ntuple_v5_mcc9_v29e_dl_runs1to3_extbnb.root").c_str());
  TTree* text = (TTree*)fext.Get("EventTree");


  float targetPOT = 4.4e+19;
  float BNBspills = 9764047.0;
  float tnuPOTsum = 0.;
  float tnuePOTsum = 0.;
  float totGoodPOT;

  tnuPOT -> SetBranchAddress("totGoodPOT", &totGoodPOT);
  for(int i = 0; i < tnuPOT->GetEntries(); ++i){
    tnuPOT -> GetEntry(i);
    tnuPOTsum += totGoodPOT;
  }

  tnuePOT -> SetBranchAddress("totGoodPOT", &totGoodPOT);
  for(int i = 0; i < tnuePOT->GetEntries(); ++i){
    tnuePOT -> GetEntry(i);
    tnuePOTsum += totGoodPOT;
  }

  float nEXTtrigs = 34202767.0+38971237.0+465951.0+59572045.0+22166992.0+36721376.0+14817082.0+39195178.0+58677653.0+19214565.0+18619185.0;


  float xsecWeight;
  int trueNuPDG;
  int trueNuCCNC;
  int foundVertex;
  int vtxIsFiducial;
  float vtxFracHitsOnCosmic;
  int nTracks;
  int trackIsSecondary[100];
  int trackClassified[100];
  int trackPID[100];
  float trackMuScore[100];
  int nShowers;
  int showerIsSecondary[100];
  int showerClassified[100];
  int showerPID[100];
  float showerCharge[100];
  float showerCosTheta[100];
  float showerPurity[100];
  float showerElScore[100];
  float showerPiScore[100];
  float showerPhScore[100];
  int showerProcess[100];
  float showerPrimaryScore[100];
  float showerFromNeutralScore[100];
  float showerFromChargedScore[100];


  tnu -> SetBranchAddress("xsecWeight", &xsecWeight);
  tnu -> SetBranchAddress("trueNuPDG", &trueNuPDG);
  tnu -> SetBranchAddress("trueNuCCNC", &trueNuCCNC);
  tnu -> SetBranchAddress("foundVertex", &foundVertex);
  tnu -> SetBranchAddress("vtxIsFiducial", &vtxIsFiducial);
  tnu -> SetBranchAddress("vtxFracHitsOnCosmic", &vtxFracHitsOnCosmic);
  tnu -> SetBranchAddress("nTracks", &nTracks);
  tnu -> SetBranchAddress("trackIsSecondary", trackIsSecondary);
  tnu -> SetBranchAddress("trackClassified", trackClassified);
  tnu -> SetBranchAddress("trackPID", trackPID);
  tnu -> SetBranchAddress("trackMuScore", trackMuScore);
  tnu -> SetBranchAddress("nShowers", &nShowers);
  tnu -> SetBranchAddress("showerIsSecondary", showerIsSecondary);
  tnu -> SetBranchAddress("showerClassified", showerClassified);
  tnu -> SetBranchAddress("showerPID", showerPID);
  tnu -> SetBranchAddress("showerCharge", showerCharge);
  tnu -> SetBranchAddress("showerCosTheta", showerCosTheta);
  tnu -> SetBranchAddress("showerPurity", showerPurity);
  tnu -> SetBranchAddress("showerElScore", showerElScore);
  tnu -> SetBranchAddress("showerPiScore", showerPiScore);
  tnu -> SetBranchAddress("showerPhScore", showerPhScore);
  tnu -> SetBranchAddress("showerProcess", showerProcess);
  tnu -> SetBranchAddress("showerPrimaryScore", showerPrimaryScore);
  tnu -> SetBranchAddress("showerFromNeutralScore", showerFromNeutralScore);
  tnu -> SetBranchAddress("showerFromChargedScore", showerFromChargedScore);

  for(int i = 0; i < tnu->GetEntries(); ++i){

    tnu -> GetEntry(i);

    if(std::isinf(xsecWeight)) continue;

    int eventType = -1;
    if(std::abs(trueNuPDG) == 14){
      if(trueNuCCNC == 0) eventType = 0; // CC numu
      else eventType = 1; //NC numu
    }
    if(std::abs(trueNuPDG) == 12 && trueNuCCNC == 1) eventType = 2; //NC nue

    if(eventType < 0 || foundVertex == 0 || vtxIsFiducial != 1 || vtxFracHitsOnCosmic >= 1.) continue;

    int nMuons = 0;
    float maxMuScore = -20.;
    for(int iT = 0; iT < nTracks; ++iT){
      if(trackIsSecondary[iT] == 1 || trackClassified[iT] != 1) continue;
      if(trackPID[iT] == 13) ++nMuons;
      if(trackMuScore[iT] > maxMuScore) maxMuScore = trackMuScore[iT];
    }
    if(nMuons > 0) continue;

    int nElectrons = 0;
    float elMaxQ = -1.;
    float elMaxQCosTheta = -1.;
    float elMaxQPur = -1.;
    float elMaxQConf = -1.;
    float elMaxQPiDiff = -1.;
    float elMaxQPhDiff = -1.;
    float elMaxQProc = -1.;
    float elMaxQProcConf = -1.;
    float elMaxQNtrlDiff = -1.;
    float elMaxQChgdDiff = -1.;
 
    for(int iS = 0; iS < nShowers; ++iS){
      if(showerIsSecondary[iS] == 1 || showerClassified[iS] == 0 || showerPID[iS] != 11) continue;
      ++nElectrons;
      if(showerCharge[iS] > elMaxQ){
        elMaxQ = showerCharge[iS];
        elMaxQCosTheta = showerCosTheta[iS];
        elMaxQPur = showerPurity[iS];
        elMaxQConf = showerElScore[iS] - (showerPhScore[iS] + showerPiScore[iS])/2.;
        elMaxQPiDiff = showerElScore[iS] - showerPiScore[iS];
        elMaxQPhDiff = showerElScore[iS] - showerPhScore[iS];
        elMaxQProc = showerProcess[iS];
        elMaxQProcConf = showerPrimaryScore[iS] - (showerFromNeutralScore[iS] + showerFromChargedScore[iS])/2.;
        elMaxQNtrlDiff = showerPrimaryScore[iS] - showerFromNeutralScore[iS];
        elMaxQChgdDiff = showerPrimaryScore[iS] - showerFromChargedScore[iS];
      }
    }

    //if(nElectrons > 0 && elMaxQCosTheta > 0. && elMaxQPur > pCut &&
    // elMaxQConf > sCut && elMaxQPiDiff > piCut && elMaxQPhDiff > phCut){
    if(nElectrons > 0 && elMaxQProc == 0 && maxMuScore < maxMuCut && elMaxQPur > pCut && elMaxQConf > sCut && elMaxQPiDiff > piCut && elMaxQPhDiff > phCut && elMaxQProcConf > procCut && elMaxQNtrlDiff > ntrlCut && elMaxQChgdDiff > chgdCut){
      if(eventType == 0) n_runs1to3_CCnumu_pass += xsecWeight;
      if(eventType == 1) n_runs1to3_NCnumu_pass += xsecWeight;
      if(eventType == 2) n_runs1to3_NCnue_pass += xsecWeight;
    }

  }


  tnue -> SetBranchAddress("xsecWeight", &xsecWeight);
  tnue -> SetBranchAddress("trueNuPDG", &trueNuPDG);
  tnue -> SetBranchAddress("trueNuCCNC", &trueNuCCNC);
  tnue -> SetBranchAddress("foundVertex", &foundVertex);
  tnue -> SetBranchAddress("vtxIsFiducial", &vtxIsFiducial);
  tnue -> SetBranchAddress("vtxFracHitsOnCosmic", &vtxFracHitsOnCosmic);
  tnue -> SetBranchAddress("nTracks", &nTracks);
  tnue -> SetBranchAddress("trackIsSecondary", trackIsSecondary);
  tnue -> SetBranchAddress("trackClassified", trackClassified);
  tnue -> SetBranchAddress("trackPID", trackPID);
  tnue -> SetBranchAddress("trackMuScore", trackMuScore);
  tnue -> SetBranchAddress("nShowers", &nShowers);
  tnue -> SetBranchAddress("showerIsSecondary", showerIsSecondary);
  tnue -> SetBranchAddress("showerClassified", showerClassified);
  tnue -> SetBranchAddress("showerPID", showerPID);
  tnue -> SetBranchAddress("showerCharge", showerCharge);
  tnue -> SetBranchAddress("showerCosTheta", showerCosTheta);
  tnue -> SetBranchAddress("showerPurity", showerPurity);
  tnue -> SetBranchAddress("showerElScore", showerElScore);
  tnue -> SetBranchAddress("showerPiScore", showerPiScore);
  tnue -> SetBranchAddress("showerPhScore", showerPhScore);
  tnue -> SetBranchAddress("showerProcess", showerProcess);
  tnue -> SetBranchAddress("showerPrimaryScore", showerPrimaryScore);
  tnue -> SetBranchAddress("showerFromNeutralScore", showerFromNeutralScore);
  tnue -> SetBranchAddress("showerFromChargedScore", showerFromChargedScore);

  for(int i = 0; i < tnue->GetEntries(); ++i){

    tnue -> GetEntry(i);

    if(std::isinf(xsecWeight) || std::abs(trueNuPDG) != 12 or trueNuCCNC != 0) continue;
    n_runs1to3_CCnue += xsecWeight;
    if(foundVertex == 0 || vtxIsFiducial != 1 || vtxFracHitsOnCosmic >= 1.) continue;

    int nMuons = 0;
    float maxMuScore = -20.;
    for(int iT = 0; iT < nTracks; ++iT){
      if(trackIsSecondary[iT] == 1 || trackClassified[iT] != 1) continue;
      if(trackPID[iT] == 13) ++nMuons;
      if(trackMuScore[iT] > maxMuScore) maxMuScore = trackMuScore[iT];
    }
    if(nMuons > 0) continue;

    int nElectrons = 0;
    float elMaxQ = -1.;
    float elMaxQCosTheta = -1.;
    float elMaxQPur = -1.;
    float elMaxQConf = -1.;
    float elMaxQPiDiff = -1.;
    float elMaxQPhDiff = -1.;
    float elMaxQProc = -1.;
    float elMaxQProcConf = -1.;
    float elMaxQNtrlDiff = -1.;
    float elMaxQChgdDiff = -1.;
 
    for(int iS = 0; iS < nShowers; ++iS){
      if(showerIsSecondary[iS] == 1 || showerClassified[iS] == 0 || showerPID[iS] != 11) continue;
      ++nElectrons;
      if(showerCharge[iS] > elMaxQ){
        elMaxQ = showerCharge[iS];
        elMaxQCosTheta = showerCosTheta[iS];
        elMaxQPur = showerPurity[iS];
        elMaxQConf = showerElScore[iS] - (showerPhScore[iS] + showerPiScore[iS])/2.;
        elMaxQPiDiff = showerElScore[iS] - showerPiScore[iS];
        elMaxQPhDiff = showerElScore[iS] - showerPhScore[iS];
        elMaxQProc = showerProcess[iS];
        elMaxQProcConf = showerPrimaryScore[iS] - (showerFromNeutralScore[iS] + showerFromChargedScore[iS])/2.;
        elMaxQNtrlDiff = showerPrimaryScore[iS] - showerFromNeutralScore[iS];
        elMaxQChgdDiff = showerPrimaryScore[iS] - showerFromChargedScore[iS];
      }
    }

    //if(nElectrons > 0 && elMaxQCosTheta > 0. && elMaxQPur > pCut &&
    // elMaxQConf > sCut && elMaxQPiDiff > piCut && elMaxQPhDiff > phCut){
    if(nElectrons > 0 && elMaxQProc == 0 && maxMuScore < maxMuCut && elMaxQPur > pCut && elMaxQConf > sCut && elMaxQPiDiff > piCut && elMaxQPhDiff > phCut && elMaxQProcConf > procCut && elMaxQNtrlDiff > ntrlCut && elMaxQChgdDiff > chgdCut){
      n_runs1to3_CCnue_pass += xsecWeight;
    }

  }


  text -> SetBranchAddress("foundVertex", &foundVertex);
  text -> SetBranchAddress("vtxIsFiducial", &vtxIsFiducial);
  text -> SetBranchAddress("vtxFracHitsOnCosmic", &vtxFracHitsOnCosmic);
  text -> SetBranchAddress("nTracks", &nTracks);
  text -> SetBranchAddress("trackIsSecondary", trackIsSecondary);
  text -> SetBranchAddress("trackClassified", trackClassified);
  text -> SetBranchAddress("trackPID", trackPID);
  text -> SetBranchAddress("trackMuScore", trackMuScore);
  text -> SetBranchAddress("nShowers", &nShowers);
  text -> SetBranchAddress("showerIsSecondary", showerIsSecondary);
  text -> SetBranchAddress("showerClassified", showerClassified);
  text -> SetBranchAddress("showerPID", showerPID);
  text -> SetBranchAddress("showerCharge", showerCharge);
  text -> SetBranchAddress("showerCosTheta", showerCosTheta);
  text -> SetBranchAddress("showerPurity", showerPurity);
  text -> SetBranchAddress("showerElScore", showerElScore);
  text -> SetBranchAddress("showerPiScore", showerPiScore);
  text -> SetBranchAddress("showerPhScore", showerPhScore);
  text -> SetBranchAddress("showerProcess", showerProcess);
  text -> SetBranchAddress("showerPrimaryScore", showerPrimaryScore);
  text -> SetBranchAddress("showerFromNeutralScore", showerFromNeutralScore);
  text -> SetBranchAddress("showerFromChargedScore", showerFromChargedScore);

  for(int i = 0; i < text->GetEntries(); ++i){

    text -> GetEntry(i);

    if(foundVertex == 0 || vtxIsFiducial != 1 || vtxFracHitsOnCosmic >= 1.) continue;

    int nMuons = 0;
    float maxMuScore = -20.;
    for(int iT = 0; iT < nTracks; ++iT){
      if(trackIsSecondary[iT] == 1 || trackClassified[iT] != 1) continue;
      if(trackPID[iT] == 13) ++nMuons;
      if(trackMuScore[iT] > maxMuScore) maxMuScore = trackMuScore[iT];
    }
    if(nMuons > 0) continue;

    int nElectrons = 0;
    float elMaxQ = -1.;
    float elMaxQCosTheta = -1.;
    float elMaxQPur = -1.;
    float elMaxQConf = -1.;
    float elMaxQPiDiff = -1.;
    float elMaxQPhDiff = -1.;
    float elMaxQProc = -1.;
    float elMaxQProcConf = -1.;
    float elMaxQNtrlDiff = -1.;
    float elMaxQChgdDiff = -1.;
 
    for(int iS = 0; iS < nShowers; ++iS){
      if(showerIsSecondary[iS] == 1 || showerClassified[iS] == 0 || showerPID[iS] != 11) continue;
      ++nElectrons;
      if(showerCharge[iS] > elMaxQ){
        elMaxQ = showerCharge[iS];
        elMaxQCosTheta = showerCosTheta[iS];
        elMaxQPur = showerPurity[iS];
        elMaxQConf = showerElScore[iS] - (showerPhScore[iS] + showerPiScore[iS])/2.;
        elMaxQPiDiff = showerElScore[iS] - showerPiScore[iS];
        elMaxQPhDiff = showerElScore[iS] - showerPhScore[iS];
        elMaxQProc = showerProcess[iS];
        elMaxQProcConf = showerPrimaryScore[iS] - (showerFromNeutralScore[iS] + showerFromChargedScore[iS])/2.;
        elMaxQNtrlDiff = showerPrimaryScore[iS] - showerFromNeutralScore[iS];
        elMaxQChgdDiff = showerPrimaryScore[iS] - showerFromChargedScore[iS];
      }
    }

    //if(nElectrons > 0 && elMaxQCosTheta > 0. && elMaxQPur > pCut &&
    // elMaxQConf > sCut && elMaxQPiDiff > piCut && elMaxQPhDiff > phCut){
    if(nElectrons > 0 && elMaxQProc == 0 && maxMuScore < maxMuCut && elMaxQPur > pCut && elMaxQConf > sCut && elMaxQPiDiff > piCut && elMaxQPhDiff > phCut && elMaxQProcConf > procCut && elMaxQNtrlDiff > ntrlCut && elMaxQChgdDiff > chgdCut){
      n_runs1to3_ext_pass += 1.;
    }

  }


  n_runs1to3_CCnue *= (targetPOT/tnuePOTsum);
  n_runs1to3_CCnumu_pass *= targetPOT/tnuPOTsum;
  n_runs1to3_NCnumu_pass *= targetPOT/tnuPOTsum;
  n_runs1to3_CCnue_pass *= targetPOT/tnuePOTsum;
  n_runs1to3_NCnue_pass *= targetPOT/tnuPOTsum;
  n_runs1to3_ext_pass *= BNBspills/nEXTtrigs;

  float total_passed = n_runs1to3_CCnue_pass + n_runs1to3_NCnue_pass + n_runs1to3_CCnumu_pass +
   n_runs1to3_NCnumu_pass + n_runs1to3_ext_pass;
  float purity = 0.;
  if(total_passed > 0.) purity = n_runs1to3_CCnue_pass / total_passed;
  //float purity = n_runs1to3_CCnue_pass / (n_runs1to3_CCnue_pass + n_runs1to3_NCnue_pass +
  // n_runs1to3_CCnumu_pass + n_runs1to3_NCnumu_pass + n_runs1to3_ext_pass);
  float efficiency = n_runs1to3_CCnue_pass / n_runs1to3_CCnue;

  std::vector<float> results{purity, efficiency, purity*efficiency};
  return results;

}


int main(int argc, char** argv){

  bool combineScores = false;
  bool run2dScoreCuts = false;
  bool run3dScorePurityOpt = false;
  bool singlePIDScore = false;
  bool singleProcScore = false;
  bool run2dProcScoreCuts = false;
  bool run2dConfScoreCuts = false;
  bool run2dConfMaxMuCuts = true;
  bool run4dScoreCuts = false;

  if(combineScores){

    std::vector<float> pCuts; pCuts.reserve(100);
    std::vector<float> sCuts; sCuts.reserve(100);
    std::vector<float> purity; purity.reserve(100);
    std::vector<float> efficiency; efficiency.reserve(100);
    std::vector<float> purXeff; purXeff.reserve(100);

    float pCut_i = 0.5;
    float pCut_f = 0.96;
    float pCut_delta = 0.01;
    float sCut_i = 5.;
    float sCut_f = 9.;
    float sCut_delta = 0.1;

    for(float pCut = pCut_i; pCut < (pCut_f + 0.1*pCut_delta); pCut += pCut_delta){
      for(float sCut = sCut_i; sCut < (sCut_f + 0.1*sCut_delta); sCut += sCut_delta){

        pCuts.push_back(pCut);
        sCuts.push_back(sCut);

        std::vector<float> selectionResults = runSelection(pCut, sCut, 0., 0., 0., 0., 0., -3.5);
        purity.push_back(selectionResults[0]);
        efficiency.push_back(selectionResults[1]);
        purXeff.push_back(selectionResults[2]);

        std::cout << "for pCut = " << pCut << " and sCut = " << sCut << ",  purity = "<<selectionResults[0]
         << ", efficiency = " << selectionResults[1] << ", purity*efficiency = " << selectionResults[2] << std::endl;

      }
    }

    float max_pXe = -1.;
    unsigned int iMax = 0;
    for(unsigned int i = 0; i < purXeff.size(); ++i){
      if(purXeff[i] > max_pXe){
        max_pXe = purXeff[i];
        iMax = i;
      }
    }

    std::cout << std::endl << std::endl;
    std::cout << "maximum purity*efficiency = " << max_pXe << std::endl;
    std::cout << "purity, effciency at max = " << purity[iMax] << ", " << efficiency[iMax] << std::endl;
    std::cout << "purity, score confidence cuts at max = " << pCuts[iMax] << ", " << sCuts[iMax] << std::endl;

    TFile f_out("optimize_selection_cuts_output.root","RECREATE");

    float purXeff_arr[purXeff.size()];
    float pCuts_arr[pCuts.size()];
    float sCuts_arr[sCuts.size()];
    std::copy(purXeff.begin(), purXeff.end(), purXeff_arr);
    std::copy(pCuts.begin(), pCuts.end(), pCuts_arr);
    std::copy(sCuts.begin(), sCuts.end(), sCuts_arr);
    TGraph2D g_opt(purXeff.size(), pCuts_arr, sCuts_arr, purXeff_arr);
    TCanvas cnv_opt("cnv_opt","cnv_opt");
    g_opt.Draw("surf1");
    cnv_opt.Write();

    std::ofstream arrayFile;
    arrayFile.open("optimize_selection_cuts_output_arrays.py");
    arrayFile << "pCuts = [";
    for(unsigned int i = 0; i < pCuts.size(); ++i){
      if(i == pCuts.size()-1) arrayFile << pCuts[i] << "]\n";
      else arrayFile << pCuts[i] << ", "; 
    }
    arrayFile << "sCuts = [";
    for(unsigned int i = 0; i < sCuts.size(); ++i){
      if(i == sCuts.size()-1) arrayFile << sCuts[i] << "]\n";
      else arrayFile << sCuts[i] << ", "; 
    }
    arrayFile << "purity = [";
    for(unsigned int i = 0; i < purity.size(); ++i){
      if(i == purity.size()-1) arrayFile << purity[i] << "]\n";
      else arrayFile << purity[i] << ", "; 
    }
    arrayFile << "efficiency = [";
    for(unsigned int i = 0; i < efficiency.size(); ++i){
      if(i == efficiency.size()-1) arrayFile << efficiency[i] << "]\n";
      else arrayFile << efficiency[i] << ", "; 
    }
    arrayFile << "purXeff = [";
    for(unsigned int i = 0; i < purXeff.size(); ++i){
      if(i == purXeff.size()-1) arrayFile << purXeff[i] << "]\n";
      else arrayFile << purXeff[i] << ", "; 
    }
    arrayFile.close();

  }


  if(run2dScoreCuts){

    std::vector<float> piCuts; piCuts.reserve(100);
    std::vector<float> phCuts; phCuts.reserve(100);
    std::vector<float> purity; purity.reserve(100);
    std::vector<float> efficiency; efficiency.reserve(100);
    std::vector<float> purXeff; purXeff.reserve(100);

    float piCut_i = 8.4;
    float piCut_f = 9.2;
    float piCut_delta = 0.1;
    float phCut_i = 3.2;
    float phCut_f = 4.0;
    float phCut_delta = 0.1;

    for(float piCut = piCut_i; piCut < (piCut_f + 0.1*piCut_delta); piCut += piCut_delta){
      for(float phCut = phCut_i; phCut < (phCut_f + 0.1*phCut_delta); phCut += phCut_delta){

        piCuts.push_back(piCut);
        phCuts.push_back(phCut);

        std::vector<float> selectionResults = runSelection(0., 0., piCut, phCut, 0., 0., 0., -3.5);
        purity.push_back(selectionResults[0]);
        efficiency.push_back(selectionResults[1]);
        purXeff.push_back(selectionResults[2]);

        std::cout << "for piCut = " << piCut << " and phCut = " << phCut << ",  purity = "<<selectionResults[0]
         << ", efficiency = " << selectionResults[1] << ", purity*efficiency = " << selectionResults[2] << std::endl;

      }
    }

    float max_pXe = -1.;
    unsigned int iMax = 0;
    for(unsigned int i = 0; i < purXeff.size(); ++i){
      if(purXeff[i] > max_pXe){
        max_pXe = purXeff[i];
        iMax = i;
      }
    }

    std::cout << std::endl << std::endl;
    std::cout << "maximum purity*efficiency = " << max_pXe << std::endl;
    std::cout << "purity, effciency at max = " << purity[iMax] << ", " << efficiency[iMax] << std::endl;
    std::cout << "pion, photon score cuts at max = " << piCuts[iMax] << ", " << phCuts[iMax] << std::endl;

    TFile f_out("optimize_selection_cuts_output_piPhCuts_fine.root","RECREATE");

    float purXeff_arr[purXeff.size()];
    float piCuts_arr[piCuts.size()];
    float phCuts_arr[phCuts.size()];
    std::copy(purXeff.begin(), purXeff.end(), purXeff_arr);
    std::copy(piCuts.begin(), piCuts.end(), piCuts_arr);
    std::copy(phCuts.begin(), phCuts.end(), phCuts_arr);
    TGraph2D g_opt(purXeff.size(), piCuts_arr, phCuts_arr, purXeff_arr);
    TCanvas cnv_opt("cnv_opt","cnv_opt");
    g_opt.Draw("surf1");
    cnv_opt.Write();

    std::ofstream arrayFile;
    arrayFile.open("optimize_selection_cuts_output_arrays_piPhCuts_fine.py");
    arrayFile << "piCuts = [";
    for(unsigned int i = 0; i < piCuts.size(); ++i){
      if(i == piCuts.size()-1) arrayFile << piCuts[i] << "]\n";
      else arrayFile << piCuts[i] << ", "; 
    }
    arrayFile << "phCuts = [";
    for(unsigned int i = 0; i < phCuts.size(); ++i){
      if(i == phCuts.size()-1) arrayFile << phCuts[i] << "]\n";
      else arrayFile << phCuts[i] << ", "; 
    }
    arrayFile << "purity = [";
    for(unsigned int i = 0; i < purity.size(); ++i){
      if(i == purity.size()-1) arrayFile << purity[i] << "]\n";
      else arrayFile << purity[i] << ", "; 
    }
    arrayFile << "efficiency = [";
    for(unsigned int i = 0; i < efficiency.size(); ++i){
      if(i == efficiency.size()-1) arrayFile << efficiency[i] << "]\n";
      else arrayFile << efficiency[i] << ", "; 
    }
    arrayFile << "purXeff = [";
    for(unsigned int i = 0; i < purXeff.size(); ++i){
      if(i == purXeff.size()-1) arrayFile << purXeff[i] << "]\n";
      else arrayFile << purXeff[i] << ", "; 
    }
    arrayFile.close();

  }


  if(run3dScorePurityOpt){

    std::vector<float> pCuts; pCuts.reserve(100);
    std::vector<float> piCuts; piCuts.reserve(100);
    std::vector<float> phCuts; phCuts.reserve(100);
    std::vector<float> purity; purity.reserve(100);
    std::vector<float> efficiency; efficiency.reserve(100);
    std::vector<float> purXeff; purXeff.reserve(100);

    float pCut_i = 0.7;
    float pCut_f = 0.86;
    float pCut_delta = 0.02;
    float piCut_i = 8.;
    float piCut_f = 10.;
    float piCut_delta = 0.2;
    float phCut_i = 2.;
    float phCut_f = 4.;
    float phCut_delta = 0.2;

    for(float pCut = pCut_i; pCut < (pCut_f + 0.1*pCut_delta); pCut += pCut_delta){
      for(float piCut = piCut_i; piCut < (piCut_f + 0.1*piCut_delta); piCut += piCut_delta){
        for(float phCut = phCut_i; phCut < (phCut_f + 0.1*phCut_delta); phCut += phCut_delta){

          pCuts.push_back(pCut);
          piCuts.push_back(piCut);
          phCuts.push_back(phCut);

          std::vector<float> selectionResults = runSelection(pCut, 0., piCut, phCut, 0., 0., 0., -3.5);
          purity.push_back(selectionResults[0]);
          efficiency.push_back(selectionResults[1]);
          purXeff.push_back(selectionResults[2]);

          std::cout << "for pCut = " << pCut << ", piCut = " << piCut << ", and phCut = " << phCut
           << ",  purity = "<<selectionResults[0] << ", efficiency = " << selectionResults[1]
           << ", purity*efficiency = " << selectionResults[2] << std::endl;

        }
      }
    }

    float max_pXe = -1.;
    unsigned int iMax = 0;
    for(unsigned int i = 0; i < purXeff.size(); ++i){
      if(purXeff[i] > max_pXe){
        max_pXe = purXeff[i];
        iMax = i;
      }
    }

    std::cout << std::endl << std::endl;
    std::cout << "maximum purity*efficiency = " << max_pXe << std::endl;
    std::cout << "purity, effciency at max = " << purity[iMax] << ", " << efficiency[iMax] << std::endl;
    std::cout << "prong purity, pion, photon score cuts at max = "
              << pCuts[iMax] << ", " << piCuts[iMax] << ", " << phCuts[iMax] << std::endl;

    std::ofstream arrayFile;
    arrayFile.open("optimize_selection_cuts_output_arrays_purPiPhCuts.py");
    arrayFile << "pCuts = [";
    for(unsigned int i = 0; i < pCuts.size(); ++i){
      if(i == pCuts.size()-1) arrayFile << pCuts[i] << "]\n";
      else arrayFile << pCuts[i] << ", "; 
    }
    arrayFile << "piCuts = [";
    for(unsigned int i = 0; i < piCuts.size(); ++i){
      if(i == piCuts.size()-1) arrayFile << piCuts[i] << "]\n";
      else arrayFile << piCuts[i] << ", "; 
    }
    arrayFile << "phCuts = [";
    for(unsigned int i = 0; i < phCuts.size(); ++i){
      if(i == phCuts.size()-1) arrayFile << phCuts[i] << "]\n";
      else arrayFile << phCuts[i] << ", "; 
    }
    arrayFile << "purity = [";
    for(unsigned int i = 0; i < purity.size(); ++i){
      if(i == purity.size()-1) arrayFile << purity[i] << "]\n";
      else arrayFile << purity[i] << ", "; 
    }
    arrayFile << "efficiency = [";
    for(unsigned int i = 0; i < efficiency.size(); ++i){
      if(i == efficiency.size()-1) arrayFile << efficiency[i] << "]\n";
      else arrayFile << efficiency[i] << ", "; 
    }
    arrayFile << "purXeff = [";
    for(unsigned int i = 0; i < purXeff.size(); ++i){
      if(i == purXeff.size()-1) arrayFile << purXeff[i] << "]\n";
      else arrayFile << purXeff[i] << ", "; 
    }
    arrayFile.close();

  }


  if(singlePIDScore){

    std::vector<float> sCuts; sCuts.reserve(100);
    std::vector<float> purity; purity.reserve(100);
    std::vector<float> efficiency; efficiency.reserve(100);
    std::vector<float> purXeff; purXeff.reserve(100);

    float sCut_i = 4.;
    float sCut_f = 9.;
    if(argc == 4){
      sCut_i = std::atof(argv[1]);
      sCut_f = std::atof(argv[2]);
    }
    float sCut_delta = 0.1;

    for(float sCut = sCut_i; sCut < (sCut_f + 0.1*sCut_delta); sCut += sCut_delta){

      sCuts.push_back(sCut);

      std::vector<float> selectionResults = runSelection(0., sCut, 0., 0., 0., 0., 0., -3.5);
      purity.push_back(selectionResults[0]);
      efficiency.push_back(selectionResults[1]);
      purXeff.push_back(selectionResults[2]);

      std::cout << "for sCut = " << sCut << ",  purity = "<<selectionResults[0]
       << ", efficiency = " << selectionResults[1] << ", purity*efficiency = " << selectionResults[2] << std::endl;

    }

    float max_pXe = -1.;
    unsigned int iMax = 0;
    for(unsigned int i = 0; i < purXeff.size(); ++i){
      if(purXeff[i] > max_pXe){
        max_pXe = purXeff[i];
        iMax = i;
      }
    }

    std::cout << std::endl << std::endl;
    std::cout << "maximum purity*efficiency = " << max_pXe << std::endl;
    std::cout << "purity, effciency at max = " << purity[iMax] << ", " << efficiency[iMax] << std::endl;
    std::cout << "score confidence cut at max = " << sCuts[iMax] << std::endl;

    std::ofstream arrayFile;
    std::string filename = "optimize_selection_cuts_output_arrays_1dPIDConfCut.py";
    if(argc == 4){
      std::string filetag(argv[3]);
      filename = "optimize_selection_cuts_output_arrays_1dPIDConfCut_"+filetag+".py";
    }
    arrayFile.open(filename.c_str());
    arrayFile << "sCuts = [";
    for(unsigned int i = 0; i < sCuts.size(); ++i){
      if(i == sCuts.size()-1) arrayFile << sCuts[i] << "]\n";
      else arrayFile << sCuts[i] << ", "; 
    }
    arrayFile << "purity = [";
    for(unsigned int i = 0; i < purity.size(); ++i){
      if(i == purity.size()-1) arrayFile << purity[i] << "]\n";
      else arrayFile << purity[i] << ", "; 
    }
    arrayFile << "efficiency = [";
    for(unsigned int i = 0; i < efficiency.size(); ++i){
      if(i == efficiency.size()-1) arrayFile << efficiency[i] << "]\n";
      else arrayFile << efficiency[i] << ", "; 
    }
    arrayFile << "purXeff = [";
    for(unsigned int i = 0; i < purXeff.size(); ++i){
      if(i == purXeff.size()-1) arrayFile << purXeff[i] << "]\n";
      else arrayFile << purXeff[i] << ", "; 
    }
    arrayFile.close();

  }


  if(singleProcScore){

    std::vector<float> sCuts; sCuts.reserve(100);
    std::vector<float> purity; purity.reserve(100);
    std::vector<float> efficiency; efficiency.reserve(100);
    std::vector<float> purXeff; purXeff.reserve(100);

    float sCut_i = 2.;
    float sCut_f = 8.;
    float sCut_delta = 0.1;

    for(float sCut = sCut_i; sCut < (sCut_f + 0.1*sCut_delta); sCut += sCut_delta){

      sCuts.push_back(sCut);

      std::vector<float> selectionResults = runSelection(0., 0., 0., 0., sCut, 0., 0., -3.5);
      purity.push_back(selectionResults[0]);
      efficiency.push_back(selectionResults[1]);
      purXeff.push_back(selectionResults[2]);

      std::cout << "for sCut = " << sCut << ",  purity = "<<selectionResults[0]
       << ", efficiency = " << selectionResults[1] << ", purity*efficiency = " << selectionResults[2] << std::endl;

    }

    float max_pXe = -1.;
    unsigned int iMax = 0;
    for(unsigned int i = 0; i < purXeff.size(); ++i){
      if(purXeff[i] > max_pXe){
        max_pXe = purXeff[i];
        iMax = i;
      }
    }

    std::cout << std::endl << std::endl;
    std::cout << "maximum purity*efficiency = " << max_pXe << std::endl;
    std::cout << "purity, effciency at max = " << purity[iMax] << ", " << efficiency[iMax] << std::endl;
    std::cout << "score confidence cut at max = " << sCuts[iMax] << std::endl;

    std::ofstream arrayFile;
    arrayFile.open("optimize_selection_cuts_output_arrays_1dProcConfCut.py");
    arrayFile << "sCuts = [";
    for(unsigned int i = 0; i < sCuts.size(); ++i){
      if(i == sCuts.size()-1) arrayFile << sCuts[i] << "]\n";
      else arrayFile << sCuts[i] << ", "; 
    }
    arrayFile << "purity = [";
    for(unsigned int i = 0; i < purity.size(); ++i){
      if(i == purity.size()-1) arrayFile << purity[i] << "]\n";
      else arrayFile << purity[i] << ", "; 
    }
    arrayFile << "efficiency = [";
    for(unsigned int i = 0; i < efficiency.size(); ++i){
      if(i == efficiency.size()-1) arrayFile << efficiency[i] << "]\n";
      else arrayFile << efficiency[i] << ", "; 
    }
    arrayFile << "purXeff = [";
    for(unsigned int i = 0; i < purXeff.size(); ++i){
      if(i == purXeff.size()-1) arrayFile << purXeff[i] << "]\n";
      else arrayFile << purXeff[i] << ", "; 
    }
    arrayFile.close();

  }


  if(run2dConfScoreCuts){

    std::vector<float> pCuts; pCuts.reserve(100);
    std::vector<float> sCuts; sCuts.reserve(100);
    std::vector<float> purity; purity.reserve(100);
    std::vector<float> efficiency; efficiency.reserve(100);
    std::vector<float> purXeff; purXeff.reserve(100);

    float pCut_i = 4.4;
    float pCut_f = 6.4;
    float pCut_delta = 0.1;
    float sCut_i = 6.1;
    float sCut_f = 8.1;
    float sCut_delta = 0.1;

    for(float pCut = pCut_i; pCut < (pCut_f + 0.1*pCut_delta); pCut += pCut_delta){
      for(float sCut = sCut_i; sCut < (sCut_f + 0.1*sCut_delta); sCut += sCut_delta){

        pCuts.push_back(pCut);
        sCuts.push_back(sCut);

        std::vector<float> selectionResults = runSelection(0., sCut, 0., 0., pCut, 0., 0., -3.5);
        purity.push_back(selectionResults[0]);
        efficiency.push_back(selectionResults[1]);
        purXeff.push_back(selectionResults[2]);

        std::cout << "for pCut = " << pCut << " and sCut = " << sCut << ",  purity = "<<selectionResults[0]
         << ", efficiency = " << selectionResults[1] << ", purity*efficiency = " << selectionResults[2] << std::endl;

      }
    }

    float max_pXe = -1.;
    unsigned int iMax = 0;
    for(unsigned int i = 0; i < purXeff.size(); ++i){
      if(purXeff[i] > max_pXe){
        max_pXe = purXeff[i];
        iMax = i;
      }
    }

    std::cout << std::endl << std::endl;
    std::cout << "maximum purity*efficiency = " << max_pXe << std::endl;
    std::cout << "purity, effciency at max = " << purity[iMax] << ", " << efficiency[iMax] << std::endl;
    std::cout << "process, PID score confidence cuts at max = " << pCuts[iMax] << ", " << sCuts[iMax] << std::endl;

    TFile f_out("optimize_selection_cuts_output_2dConfScoreCuts.root","RECREATE");

    float purXeff_arr[purXeff.size()];
    float pCuts_arr[pCuts.size()];
    float sCuts_arr[sCuts.size()];
    std::copy(purXeff.begin(), purXeff.end(), purXeff_arr);
    std::copy(pCuts.begin(), pCuts.end(), pCuts_arr);
    std::copy(sCuts.begin(), sCuts.end(), sCuts_arr);
    TGraph2D g_opt(purXeff.size(), pCuts_arr, sCuts_arr, purXeff_arr);
    TCanvas cnv_opt("cnv_opt","cnv_opt");
    g_opt.Draw("surf1");
    cnv_opt.Write();

    std::ofstream arrayFile;
    arrayFile.open("optimize_selection_cuts_output_arrays_2dConfScoreCuts.py");
    arrayFile << "pCuts = [";
    for(unsigned int i = 0; i < pCuts.size(); ++i){
      if(i == pCuts.size()-1) arrayFile << pCuts[i] << "]\n";
      else arrayFile << pCuts[i] << ", "; 
    }
    arrayFile << "sCuts = [";
    for(unsigned int i = 0; i < sCuts.size(); ++i){
      if(i == sCuts.size()-1) arrayFile << sCuts[i] << "]\n";
      else arrayFile << sCuts[i] << ", "; 
    }
    arrayFile << "purity = [";
    for(unsigned int i = 0; i < purity.size(); ++i){
      if(i == purity.size()-1) arrayFile << purity[i] << "]\n";
      else arrayFile << purity[i] << ", "; 
    }
    arrayFile << "efficiency = [";
    for(unsigned int i = 0; i < efficiency.size(); ++i){
      if(i == efficiency.size()-1) arrayFile << efficiency[i] << "]\n";
      else arrayFile << efficiency[i] << ", "; 
    }
    arrayFile << "purXeff = [";
    for(unsigned int i = 0; i < purXeff.size(); ++i){
      if(i == purXeff.size()-1) arrayFile << purXeff[i] << "]\n";
      else arrayFile << purXeff[i] << ", "; 
    }
    arrayFile.close();

  }


  if(run2dProcScoreCuts){

    std::vector<float> chgdCuts; chgdCuts.reserve(100);
    std::vector<float> ntrlCuts; ntrlCuts.reserve(100);
    std::vector<float> purity; purity.reserve(100);
    std::vector<float> efficiency; efficiency.reserve(100);
    std::vector<float> purXeff; purXeff.reserve(100);

    float chgdCut_i = 5.4;
    float chgdCut_f = 6.2;
    float chgdCut_delta = 0.1;
    float ntrlCut_i = 2.8;
    float ntrlCut_f = 3.6;
    float ntrlCut_delta = 0.1;

    for(float chgdCut = chgdCut_i; chgdCut < (chgdCut_f + 0.1*chgdCut_delta); chgdCut += chgdCut_delta){
      for(float ntrlCut = ntrlCut_i; ntrlCut < (ntrlCut_f + 0.1*ntrlCut_delta); ntrlCut += ntrlCut_delta){

        chgdCuts.push_back(chgdCut);
        ntrlCuts.push_back(ntrlCut);

        std::vector<float> selectionResults = runSelection(0., 0., 0., 0., 0., chgdCut, ntrlCut, -3.5);
        purity.push_back(selectionResults[0]);
        efficiency.push_back(selectionResults[1]);
        purXeff.push_back(selectionResults[2]);

        std::cout << "for chgdCut = " << chgdCut << " and ntrlCut = " << ntrlCut << ",  purity = "<<selectionResults[0]
         << ", efficiency = " << selectionResults[1] << ", purity*efficiency = " << selectionResults[2] << std::endl;

      }
    }

    float max_pXe = -1.;
    unsigned int iMax = 0;
    for(unsigned int i = 0; i < purXeff.size(); ++i){
      if(purXeff[i] > max_pXe){
        max_pXe = purXeff[i];
        iMax = i;
      }
    }

    std::cout << std::endl << std::endl;
    std::cout << "maximum purity*efficiency = " << max_pXe << std::endl;
    std::cout << "purity, effciency at max = " << purity[iMax] << ", " << efficiency[iMax] << std::endl;
    std::cout << "charged, neutral parent process score cuts at max = " << chgdCuts[iMax] << ", " << ntrlCuts[iMax] << std::endl;

    TFile f_out("optimize_selection_cuts_output_2dProcScoreCuts_fine.root","RECREATE");

    float purXeff_arr[purXeff.size()];
    float chgdCuts_arr[chgdCuts.size()];
    float ntrlCuts_arr[ntrlCuts.size()];
    std::copy(purXeff.begin(), purXeff.end(), purXeff_arr);
    std::copy(chgdCuts.begin(), chgdCuts.end(), chgdCuts_arr);
    std::copy(ntrlCuts.begin(), ntrlCuts.end(), ntrlCuts_arr);
    TGraph2D g_opt(purXeff.size(), chgdCuts_arr, ntrlCuts_arr, purXeff_arr);
    TCanvas cnv_opt("cnv_opt","cnv_opt");
    g_opt.Draw("surf1");
    cnv_opt.Write();

    std::ofstream arrayFile;
    arrayFile.open("optimize_selection_cuts_output_arrays_2dProcScoreCuts_fine.py");
    arrayFile << "chgdCuts = [";
    for(unsigned int i = 0; i < chgdCuts.size(); ++i){
      if(i == chgdCuts.size()-1) arrayFile << chgdCuts[i] << "]\n";
      else arrayFile << chgdCuts[i] << ", "; 
    }
    arrayFile << "ntrlCuts = [";
    for(unsigned int i = 0; i < ntrlCuts.size(); ++i){
      if(i == ntrlCuts.size()-1) arrayFile << ntrlCuts[i] << "]\n";
      else arrayFile << ntrlCuts[i] << ", "; 
    }
    arrayFile << "purity = [";
    for(unsigned int i = 0; i < purity.size(); ++i){
      if(i == purity.size()-1) arrayFile << purity[i] << "]\n";
      else arrayFile << purity[i] << ", "; 
    }
    arrayFile << "efficiency = [";
    for(unsigned int i = 0; i < efficiency.size(); ++i){
      if(i == efficiency.size()-1) arrayFile << efficiency[i] << "]\n";
      else arrayFile << efficiency[i] << ", "; 
    }
    arrayFile << "purXeff = [";
    for(unsigned int i = 0; i < purXeff.size(); ++i){
      if(i == purXeff.size()-1) arrayFile << purXeff[i] << "]\n";
      else arrayFile << purXeff[i] << ", "; 
    }
    arrayFile.close();

  }


  if(run2dConfMaxMuCuts){

    std::vector<float> mCuts; mCuts.reserve(100);
    std::vector<float> sCuts; sCuts.reserve(100);
    std::vector<float> purity; purity.reserve(100);
    std::vector<float> efficiency; efficiency.reserve(100);
    std::vector<float> purXeff; purXeff.reserve(100);

    float mCut_i = -6.5;
    float mCut_f = -0.5;
    float mCut_delta = 0.1;
    float sCut_i = 3.0;
    float sCut_f = 9.0;
    float sCut_delta = 0.1;

    for(float mCut = mCut_i; mCut < (mCut_f + 0.1*mCut_delta); mCut += mCut_delta){
      for(float sCut = sCut_i; sCut < (sCut_f + 0.1*sCut_delta); sCut += sCut_delta){

        mCuts.push_back(mCut);
        sCuts.push_back(sCut);

        std::vector<float> selectionResults = runSelection(0., sCut, 0., 0., 0., 0., 0., mCut);
        purity.push_back(selectionResults[0]);
        efficiency.push_back(selectionResults[1]);
        purXeff.push_back(selectionResults[2]);

        std::cout << "for mCut = " << mCut << " and sCut = " << sCut << ",  purity = "<<selectionResults[0]
         << ", efficiency = " << selectionResults[1] << ", purity*efficiency = " << selectionResults[2] << std::endl;

      }
    }

    float max_pXe = -1.;
    unsigned int iMax = 0;
    for(unsigned int i = 0; i < purXeff.size(); ++i){
      if(purXeff[i] > max_pXe){
        max_pXe = purXeff[i];
        iMax = i;
      }
    }

    std::cout << std::endl << std::endl;
    std::cout << "maximum purity*efficiency = " << max_pXe << std::endl;
    std::cout << "purity, effciency at max = " << purity[iMax] << ", " << efficiency[iMax] << std::endl;
    std::cout << "max muon, PID score confidence cuts at max = " << mCuts[iMax] << ", " << sCuts[iMax] << std::endl;

    TFile f_out("optimize_selection_cuts_output_2dConfMaxMuCuts.root","RECREATE");

    float purXeff_arr[purXeff.size()];
    float mCuts_arr[mCuts.size()];
    float sCuts_arr[sCuts.size()];
    std::copy(purXeff.begin(), purXeff.end(), purXeff_arr);
    std::copy(mCuts.begin(), mCuts.end(), mCuts_arr);
    std::copy(sCuts.begin(), sCuts.end(), sCuts_arr);
    TGraph2D g_opt(purXeff.size(), mCuts_arr, sCuts_arr, purXeff_arr);
    TCanvas cnv_opt("cnv_opt","cnv_opt");
    g_opt.Draw("surf1");
    cnv_opt.Write();

    std::ofstream arrayFile;
    arrayFile.open("optimize_selection_cuts_output_arrays_2dConfMaxMuCuts.py");
    arrayFile << "mCuts = [";
    for(unsigned int i = 0; i < mCuts.size(); ++i){
      if(i == mCuts.size()-1) arrayFile << mCuts[i] << "]\n";
      else arrayFile << mCuts[i] << ", "; 
    }
    arrayFile << "sCuts = [";
    for(unsigned int i = 0; i < sCuts.size(); ++i){
      if(i == sCuts.size()-1) arrayFile << sCuts[i] << "]\n";
      else arrayFile << sCuts[i] << ", "; 
    }
    arrayFile << "purity = [";
    for(unsigned int i = 0; i < purity.size(); ++i){
      if(i == purity.size()-1) arrayFile << purity[i] << "]\n";
      else arrayFile << purity[i] << ", "; 
    }
    arrayFile << "efficiency = [";
    for(unsigned int i = 0; i < efficiency.size(); ++i){
      if(i == efficiency.size()-1) arrayFile << efficiency[i] << "]\n";
      else arrayFile << efficiency[i] << ", "; 
    }
    arrayFile << "purXeff = [";
    for(unsigned int i = 0; i < purXeff.size(); ++i){
      if(i == purXeff.size()-1) arrayFile << purXeff[i] << "]\n";
      else arrayFile << purXeff[i] << ", "; 
    }
    arrayFile.close();

  }


  if(run4dScoreCuts){

    std::vector<float> chgdCuts; chgdCuts.reserve(100);
    std::vector<float> ntrlCuts; ntrlCuts.reserve(100);
    std::vector<float> piCuts; piCuts.reserve(100);
    std::vector<float> phCuts; phCuts.reserve(100);
    std::vector<float> purity; purity.reserve(100);
    std::vector<float> efficiency; efficiency.reserve(100);
    std::vector<float> purXeff; purXeff.reserve(100);

    float chgdCut_i = 1.;
    float chgdCut_f = 6.;
    float chgdCut_delta = 0.2;
    float ntrlCut_i = 1.;
    float ntrlCut_f = 5.;
    float ntrlCut_delta = 0.2;
    float piCut_i = 7.;
    float piCut_f = 11.;
    float piCut_delta = 0.2;
    float phCut_i = 1.;
    float phCut_f = 5.;
    float phCut_delta = 0.2;

    for(float chgdCut = chgdCut_i; chgdCut < (chgdCut_f + 0.1*chgdCut_delta); chgdCut += chgdCut_delta){
      for(float ntrlCut = ntrlCut_i; ntrlCut < (ntrlCut_f + 0.1*ntrlCut_delta); ntrlCut += ntrlCut_delta){
        for(float piCut = piCut_i; piCut < (piCut_f + 0.1*piCut_delta); piCut += piCut_delta){
          for(float phCut = phCut_i; phCut < (phCut_f + 0.1*phCut_delta); phCut += phCut_delta){

            chgdCuts.push_back(chgdCut);
            ntrlCuts.push_back(ntrlCut);
            piCuts.push_back(piCut);
            phCuts.push_back(phCut);

            std::vector<float> selectionResults = runSelection(0., 0., piCut, phCut, 0., chgdCut, ntrlCut, -3.5);
            purity.push_back(selectionResults[0]);
            efficiency.push_back(selectionResults[1]);
            purXeff.push_back(selectionResults[2]);

            std::cout << "for chgdCut = " << chgdCut << ", ntrlCut = " << ntrlCut
             << ", piCut = " << piCut << ", and phCut = " << phCut
             << ",  purity = "<<selectionResults[0] << ", efficiency = " << selectionResults[1]
             << ", purity*efficiency = " << selectionResults[2] << std::endl;

          }
        }
      }
    }

    float max_pXe = -1.;
    unsigned int iMax = 0;
    for(unsigned int i = 0; i < purXeff.size(); ++i){
      if(purXeff[i] > max_pXe){
        max_pXe = purXeff[i];
        iMax = i;
      }
    }

    std::cout << std::endl << std::endl;
    std::cout << "maximum purity*efficiency = " << max_pXe << std::endl;
    std::cout << "purity, effciency at max = " << purity[iMax] << ", " << efficiency[iMax] << std::endl;
    std::cout << "charged parent, neutral parent, pion, photon score cuts at max = "
     << chgdCuts[iMax] << ", " << ntrlCuts[iMax] << ", " << piCuts[iMax] << ", " << phCuts[iMax] << std::endl;

    std::ofstream arrayFile;
    arrayFile.open("optimize_selection_cuts_output_arrays_4dScoreCuts.py");
    arrayFile << "chgdCuts = [";
    for(unsigned int i = 0; i < chgdCuts.size(); ++i){
      if(i == chgdCuts.size()-1) arrayFile << chgdCuts[i] << "]\n";
      else arrayFile << chgdCuts[i] << ", "; 
    }
    arrayFile << "ntrlCuts = [";
    for(unsigned int i = 0; i < ntrlCuts.size(); ++i){
      if(i == ntrlCuts.size()-1) arrayFile << ntrlCuts[i] << "]\n";
      else arrayFile << ntrlCuts[i] << ", "; 
    }
    arrayFile << "piCuts = [";
    for(unsigned int i = 0; i < piCuts.size(); ++i){
      if(i == piCuts.size()-1) arrayFile << piCuts[i] << "]\n";
      else arrayFile << piCuts[i] << ", "; 
    }
    arrayFile << "phCuts = [";
    for(unsigned int i = 0; i < phCuts.size(); ++i){
      if(i == phCuts.size()-1) arrayFile << phCuts[i] << "]\n";
      else arrayFile << phCuts[i] << ", "; 
    }
    arrayFile << "purity = [";
    for(unsigned int i = 0; i < purity.size(); ++i){
      if(i == purity.size()-1) arrayFile << purity[i] << "]\n";
      else arrayFile << purity[i] << ", "; 
    }
    arrayFile << "efficiency = [";
    for(unsigned int i = 0; i < efficiency.size(); ++i){
      if(i == efficiency.size()-1) arrayFile << efficiency[i] << "]\n";
      else arrayFile << efficiency[i] << ", "; 
    }
    arrayFile << "purXeff = [";
    for(unsigned int i = 0; i < purXeff.size(); ++i){
      if(i == purXeff.size()-1) arrayFile << purXeff[i] << "]\n";
      else arrayFile << purXeff[i] << ", "; 
    }
    arrayFile.close();

  }


}


