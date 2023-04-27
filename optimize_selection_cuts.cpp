
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


std::vector<float> runSelection(float pCut, float sCut, float piCut, float phCut){

  float n_runs1to3_CCnue = 0.;
  float n_runs1to3_CCnumu_pass = 0.;
  float n_runs1to3_NCnumu_pass = 0.;
  float n_runs1to3_CCnue_pass = 0.;
  float n_runs1to3_NCnue_pass = 0.;
  float n_runs1to3_ext_pass = 0.;


  std::string infiledir = "selection_output/prepare_selection_test_output/";

  TFile fnu((infiledir+"prepare_selection_test_reco_v2me05_gen2val_v19_nu_file.root").c_str());
  TTree* tnu = (TTree*)fnu.Get("EventTree");
  TTree* tnuPOT = (TTree*)fnu.Get("potTree");

  TFile fnue((infiledir+"prepare_selection_test_reco_v2me05_gen2val_v19_nue_file.root").c_str());
  TTree* tnue = (TTree*)fnue.Get("EventTree");
  TTree* tnuePOT = (TTree*)fnue.Get("potTree");

  TFile fext((infiledir+"prepare_selection_test_reco_v2me05_gen2val_v19_extbnb_file.root").c_str());
  TTree* text = (TTree*)fext.Get("EventTree");


  float run3POT = 4.3e+19 + 1.701e+20 + 2.97e+19 + 1.524e+17;
  float runs1to3POT = 6.67e+20;
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

  //970: number of run3 extbnb merged_dlana files in prepare_selection_test_reco_v2me05_gen2val_v16_extbnb_file.root
  //17661: number of run3 extbnb merged_dlana files in prepare_selection_test_reco_v2me05_gen2val_v17+_extbnb_file.root
  //89559: number of run3 extbnb files (# in def: prod_extunbiased_swizzle_crt_inclusive_v6_v6a_goodruns_mcc9_run3)
  float textPOTsum = (17661./89559.)*run3POT;


  float xsecWeight;
  int trueNuPDG;
  int trueNuCCNC;
  int nVertices;
  int vtxIsFiducial;
  float vtxFracHitsOnCosmic;
  int nTracks;
  int trackIsSecondary[100];
  int trackClassified[100];
  int trackPID[100];
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


  tnu -> SetBranchAddress("xsecWeight", &xsecWeight);
  tnu -> SetBranchAddress("trueNuPDG", &trueNuPDG);
  tnu -> SetBranchAddress("trueNuCCNC", &trueNuCCNC);
  tnu -> SetBranchAddress("nVertices", &nVertices);
  tnu -> SetBranchAddress("vtxIsFiducial", &vtxIsFiducial);
  tnu -> SetBranchAddress("vtxFracHitsOnCosmic", &vtxFracHitsOnCosmic);
  tnu -> SetBranchAddress("nTracks", &nTracks);
  tnu -> SetBranchAddress("trackIsSecondary", trackIsSecondary);
  tnu -> SetBranchAddress("trackClassified", trackClassified);
  tnu -> SetBranchAddress("trackPID", trackPID);
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

  for(int i = 0; i < tnu->GetEntries(); ++i){

    tnu -> GetEntry(i);

    if(std::isinf(xsecWeight)) continue;

    int eventType = -1;
    if(std::abs(trueNuPDG) == 14){
      if(trueNuCCNC == 0) eventType = 0; // CC numu
      else eventType = 1; //NC numu
    }
    if(std::abs(trueNuPDG) == 12 && trueNuCCNC == 1) eventType = 2; //NC nue

    if(eventType < 0 || nVertices < 1 || vtxIsFiducial != 1 || vtxFracHitsOnCosmic >= 1.) continue;

    int nMuons = 0;
    for(int iT = 0; iT < nTracks; ++iT){
      if(trackIsSecondary[iT] != 1 && trackClassified[iT] == 1 && trackPID[iT] == 13) ++nMuons;
    }
    if(nMuons > 0) continue;

    int nElectrons = 0;
    float elMaxQ = -1.;
    float elMaxQCosTheta = -1.;
    float elMaxQPur = -1.;
    float elMaxQConf = -1.;
    float elMaxQPiDiff = -1.;
    float elMaxQPhDiff = -1.;
 
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
      }
    }

    if(nElectrons > 0 && elMaxQCosTheta > 0. && elMaxQPur > pCut &&
     elMaxQConf > sCut && elMaxQPiDiff > piCut && elMaxQPhDiff > phCut){
      if(eventType == 0) n_runs1to3_CCnumu_pass += xsecWeight;
      if(eventType == 1) n_runs1to3_NCnumu_pass += xsecWeight;
      if(eventType == 2) n_runs1to3_NCnue_pass += xsecWeight;
    }

  }


  tnue -> SetBranchAddress("xsecWeight", &xsecWeight);
  tnue -> SetBranchAddress("trueNuPDG", &trueNuPDG);
  tnue -> SetBranchAddress("trueNuCCNC", &trueNuCCNC);
  tnue -> SetBranchAddress("nVertices", &nVertices);
  tnue -> SetBranchAddress("vtxIsFiducial", &vtxIsFiducial);
  tnue -> SetBranchAddress("vtxFracHitsOnCosmic", &vtxFracHitsOnCosmic);
  tnue -> SetBranchAddress("nTracks", &nTracks);
  tnue -> SetBranchAddress("trackIsSecondary", trackIsSecondary);
  tnue -> SetBranchAddress("trackClassified", trackClassified);
  tnue -> SetBranchAddress("trackPID", trackPID);
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

  for(int i = 0; i < tnue->GetEntries(); ++i){

    tnue -> GetEntry(i);

    if(std::isinf(xsecWeight) || std::abs(trueNuPDG) != 12 or trueNuCCNC != 0) continue;
    n_runs1to3_CCnue += xsecWeight;
    if(nVertices < 1 || vtxIsFiducial != 1 || vtxFracHitsOnCosmic >= 1.) continue;

    int nMuons = 0;
    for(int iT = 0; iT < nTracks; ++iT){
      if(trackIsSecondary[iT] != 1 && trackClassified[iT] == 1 && trackPID[iT] == 13) ++nMuons;
    }
    if(nMuons > 0) continue;

    int nElectrons = 0;
    float elMaxQ = -1.;
    float elMaxQCosTheta = -1.;
    float elMaxQPur = -1.;
    float elMaxQConf = -1.;
    float elMaxQPiDiff = -1.;
    float elMaxQPhDiff = -1.;
 
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
      }
    }

    if(nElectrons > 0 && elMaxQCosTheta > 0. && elMaxQPur > pCut &&
     elMaxQConf > sCut && elMaxQPiDiff > piCut && elMaxQPhDiff > phCut){
      n_runs1to3_CCnue_pass += xsecWeight;
    }

  }


  text -> SetBranchAddress("xsecWeight", &xsecWeight);
  text -> SetBranchAddress("trueNuPDG", &trueNuPDG);
  text -> SetBranchAddress("trueNuCCNC", &trueNuCCNC);
  text -> SetBranchAddress("nVertices", &nVertices);
  text -> SetBranchAddress("vtxIsFiducial", &vtxIsFiducial);
  text -> SetBranchAddress("vtxFracHitsOnCosmic", &vtxFracHitsOnCosmic);
  text -> SetBranchAddress("nTracks", &nTracks);
  text -> SetBranchAddress("trackIsSecondary", trackIsSecondary);
  text -> SetBranchAddress("trackClassified", trackClassified);
  text -> SetBranchAddress("trackPID", trackPID);
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

  for(int i = 0; i < text->GetEntries(); ++i){

    text -> GetEntry(i);

    if(nVertices < 1 || vtxIsFiducial != 1 || vtxFracHitsOnCosmic >= 1.) continue;

    int nMuons = 0;
    for(int iT = 0; iT < nTracks; ++iT){
      if(trackIsSecondary[iT] != 1 && trackClassified[iT] == 1 && trackPID[iT] == 13) ++nMuons;
    }
    if(nMuons > 0) continue;

    int nElectrons = 0;
    float elMaxQ = -1.;
    float elMaxQCosTheta = -1.;
    float elMaxQPur = -1.;
    float elMaxQConf = -1.;
    float elMaxQPiDiff = -1.;
    float elMaxQPhDiff = -1.;
 
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
      }
    }

    if(nElectrons > 0 && elMaxQCosTheta > 0. && elMaxQPur > pCut &&
     elMaxQConf > sCut && elMaxQPiDiff > piCut && elMaxQPhDiff > phCut){
      n_runs1to3_ext_pass += 1.;
    }

  }


  n_runs1to3_CCnue *= (runs1to3POT/tnuePOTsum);
  n_runs1to3_CCnumu_pass *= runs1to3POT/tnuPOTsum;
  n_runs1to3_NCnumu_pass *= runs1to3POT/tnuPOTsum;
  n_runs1to3_CCnue_pass *= runs1to3POT/tnuePOTsum;
  n_runs1to3_NCnue_pass *= runs1to3POT/tnuPOTsum;
  n_runs1to3_ext_pass *= runs1to3POT/textPOTsum;

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
  bool run3dScorePurityOpt = true;

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

        std::vector<float> selectionResults = runSelection(pCut, sCut, 0., 0.);
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

    float piCut_i = 5.;
    float piCut_f = 12.;
    float piCut_delta = 0.1;
    float phCut_i = 0.;
    float phCut_f = 6.;
    float phCut_delta = 0.1;

    for(float piCut = piCut_i; piCut < (piCut_f + 0.1*piCut_delta); piCut += piCut_delta){
      for(float phCut = phCut_i; phCut < (phCut_f + 0.1*phCut_delta); phCut += phCut_delta){

        piCuts.push_back(piCut);
        phCuts.push_back(phCut);

        std::vector<float> selectionResults = runSelection(0., 0., piCut, phCut);
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

    TFile f_out("optimize_selection_cuts_output_piPhCuts.root","RECREATE");

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
    arrayFile.open("optimize_selection_cuts_output_arrays_piPhCuts.py");
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

          std::vector<float> selectionResults = runSelection(pCut, 0., piCut, phCut);
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


}


