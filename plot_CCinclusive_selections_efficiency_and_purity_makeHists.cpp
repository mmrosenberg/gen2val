
#include <iostream>
//#include <fstream>
//#include <vector>
#include <cmath>
//#include <unistd.h>
//#include <getopt.h>
#include <cstdlib>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TObject.h>


struct SelResults {
  bool foundVtx, foundNuVtx;
  bool noMu, noMuLike;
  bool primEl, confPrimEl;
  bool passedCCnue, passedCCnumu;
  SelResults() : foundVtx(false), foundNuVtx(false),
                 noMu(false), noMuLike(false), primEl(false), confPrimEl(false),
                 passedCCnue(false), passedCCnumu(false) {}
};


SelResults checkCuts(int foundVertex, int vtxIsFiducial, float vtxFracHitsOnCosmic,
 int nTracks, int nShowers, int trackIsSecondary[], int trackClassified[],
 int showerIsSecondary[], int showerClassified[], int trackPID[], int showerPID[],
 int showerProcess[], float trackMuScore[], float showerCharge[], float showerElScore[],
 float showerPhScore[], float showerPiScore[]){

  SelResults results;

  if(foundVertex == 0 || vtxIsFiducial != 1) return results;
  results.foundVtx = true;
  if(vtxFracHitsOnCosmic >= (1. - 1e-6)) return results;
  results.foundNuVtx = true;

  int nMuons = 0;
  float maxMuScore = -99.;  
  for(int iT = 0; iT < nTracks; ++iT){
    if(trackIsSecondary[iT] == 1 || trackClassified[iT] != 1) continue;
    if(trackPID[iT] == 13) ++nMuons;
    if(trackMuScore[iT] > maxMuScore) maxMuScore = trackMuScore[iT];
  }

  if(nMuons == 0) results.noMu = true;
  //if(nMuons == 0 && maxMuScore < -3.7) results.noMuLike = true;

  int nElectrons = 0;
  float elMaxQ = -99.;
  float elMaxQConf = -9.;
  int elMaxQProc = -1;
  for(int iS = 0; iS < nShowers; ++iS){
    if(showerIsSecondary[iS] == 1 || showerClassified[iS] != 1) continue;
    if(showerPID[iS] == 11){
      ++nElectrons;
      if(showerCharge[iS] > elMaxQ){
        elMaxQ = showerCharge[iS];
        elMaxQProc = showerProcess[iS];
        elMaxQConf = showerElScore[iS] - (showerPhScore[iS] + showerPiScore[iS])/2.;
      }
    }
  }

  //if(nElectrons >= 1 && elMaxQProc == 0) results.primEl = true;
  //if(nElectrons >= 1 && elMaxQProc == 0 && elMaxQConf > 7.1) results.confPrimEl = true;
  if(nMuons == 0 && nElectrons >= 1 && elMaxQProc == 0) results.primEl = true;
  if(results.primEl && maxMuScore < -3.7) results.noMuLike = true;
  if(results.noMuLike && elMaxQConf > 7.1){ results.confPrimEl = true; results.passedCCnue = true; }

  //results.passedCCnue = (results.foundNuVtx && results.noMuLike && results.confPrimEl);
  results.passedCCnumu = (results.foundNuVtx && !results.noMu);

  return results;

}



class HistFiller {

  private:

    std::string sampleType;
    float samplePOT;
    float targetPOT;

    //ntuple variables
    int trueNuPDG, trueNuCCNC;
    int foundVertex, vtxIsFiducial;
    float trueNuE, recoNuE;
    float vtxFracHitsOnCosmic;
    float xsecWeight;
    int nTracks, nShowers;
    int trackIsSecondary[100];
    int trackClassified[100];
    int showerIsSecondary[100];
    int showerClassified[100];
    int trackPID[100];
    int showerPID[100];
    int showerProcess[100];
    float trackMuScore[100];
    float showerCharge[100];
    float showerElScore[100];
    float showerPhScore[100];
    float showerPiScore[100];

    TFile* ntuple;
    TTree* eventTree;

    float nueH;
    float numuH;

  public:

    TH1F* h_nue_allTrue;
    TH1F* h_nue_pAll_foundVtx;
    TH1F* h_nue_pTrue_foundVtx;
    TH1F* h_nue_pAll_foundNuVtx;
    TH1F* h_nue_pTrue_foundNuVtx;
    TH1F* h_nue_pAll_noMu;
    TH1F* h_nue_pTrue_noMu;
    TH1F* h_nue_pAll_noMuLike;
    TH1F* h_nue_pTrue_noMuLike;
    TH1F* h_nue_pAll_primEl;
    TH1F* h_nue_pTrue_primEl;
    TH1F* h_nue_pAll_confPrimEl;
    TH1F* h_nue_pTrue_confPrimEl;
    TH1F* h_nue_pAll_allCuts;
    TH1F* h_nue_pTrue_allCuts;

    TH1F* h_numu_allTrue;
    TH1F* h_numu_pAll_foundVtx;
    TH1F* h_numu_pTrue_foundVtx;
    TH1F* h_numu_pAll_foundNuVtx;
    TH1F* h_numu_pTrue_foundNuVtx;
    TH1F* h_numu_pAll_allCuts;
    TH1F* h_numu_pTrue_allCuts;

    TH1F* hT_nue_allTrue;
    TH1F* hT_nue_pTrue_foundVtx;
    TH1F* hT_nue_pTrue_foundNuVtx;
    TH1F* hT_nue_pTrue_noMu;
    TH1F* hT_nue_pTrue_noMuLike;
    TH1F* hT_nue_pTrue_primEl;
    TH1F* hT_nue_pTrue_confPrimEl;
    TH1F* hT_nue_pTrue_allCuts;

    TH1F* hT_numu_allTrue;
    TH1F* hT_numu_pTrue_foundVtx;
    TH1F* hT_numu_pTrue_foundNuVtx;
    TH1F* hT_numu_pTrue_allCuts;


    HistFiller(std::string sampleType_, float targetPOT_=4.4e+19){

      TH1::SetDefaultSumw2(kTRUE);

      sampleType = sampleType_;
      targetPOT = targetPOT_;

      if(sampleType != "bnboverlay" && sampleType != "nueoverlay" && 
       sampleType != "allCosmics" && sampleType != "run1Cosmics"){
        std::cout << "invalid sample requested: " << sampleType << std::endl;
        exit(1);
      }

      std::string ntupleFile = "/home/matthew/microboone/tufts/gen2val/flat_ntuples/";
      if(sampleType == "bnboverlay") ntupleFile += "dlgen2_reco_v2me06_ntuple_v5_mcc9_v28_wctagger_bnboverlay.root";
      if(sampleType == "nueoverlay") ntupleFile += "dlgen2_reco_v2me06_ntuple_v5_mcc9_v28_wctagger_nueintrinsics.root";
      if(sampleType == "allCosmics") ntupleFile += "dlgen2_reco_v2me06_ntuple_v5_mcc9_v29e_dl_runs1to3_extbnb.root";
      if(sampleType == "run1Cosmics") ntupleFile += "dlgen2_reco_v2me06_ntuple_v5_mcc9_v29e_dl_run1_all_extbnb.root";

      ntuple = new TFile(ntupleFile.c_str(), "READ");
      eventTree = (TTree*)ntuple->Get("EventTree");

      eventTree->SetBranchAddress("foundVertex",&foundVertex);
      eventTree->SetBranchAddress("vtxIsFiducial",&vtxIsFiducial);
      eventTree->SetBranchAddress("recoNuE",&recoNuE);
      eventTree->SetBranchAddress("vtxFracHitsOnCosmic",&vtxFracHitsOnCosmic);
      eventTree->SetBranchAddress("nTracks",&nTracks);
      eventTree->SetBranchAddress("nShowers",&nShowers);
      eventTree->SetBranchAddress("trackIsSecondary",trackIsSecondary);
      eventTree->SetBranchAddress("trackClassified",trackClassified);
      eventTree->SetBranchAddress("showerIsSecondary",showerIsSecondary);
      eventTree->SetBranchAddress("showerClassified",showerClassified);
      eventTree->SetBranchAddress("trackPID",trackPID);
      eventTree->SetBranchAddress("showerPID",showerPID);
      eventTree->SetBranchAddress("showerProcess",showerProcess);
      eventTree->SetBranchAddress("trackMuScore",trackMuScore);
      eventTree->SetBranchAddress("showerCharge",showerCharge);
      eventTree->SetBranchAddress("showerElScore",showerElScore);
      eventTree->SetBranchAddress("showerPhScore",showerPhScore);
      eventTree->SetBranchAddress("showerPiScore",showerPiScore);

      if(sampleType == "bnboverlay" || sampleType == "nueoverlay"){

        eventTree->SetBranchAddress("xsecWeight",&xsecWeight);
        eventTree->SetBranchAddress("trueNuPDG",&trueNuPDG);
        eventTree->SetBranchAddress("trueNuCCNC",&trueNuCCNC);
        eventTree->SetBranchAddress("trueNuE",&trueNuE);

        samplePOT = 0.;
        float totGoodPOT;
        TTree* potTree = (TTree*)ntuple->Get("potTree");
        potTree -> SetBranchAddress("totGoodPOT",&totGoodPOT);
        for(int i = 0; i < potTree->GetEntries(); ++i){
          potTree -> GetEntry(i);
          samplePOT += totGoodPOT;
        }

      }

      if(sampleType == "allCosmics"){
        float EXT = 34202767.0+38971237.0+465951.0+59572045.0+22166992.0+36721376.0+
              14817082.0+39195178.0+58677653.0+19214565.0+18619185.0;
        samplePOT = (EXT/9764047.0)*4.4e+19;
      }

      if(sampleType == "run1Cosmics"){
        float EXT = 34202767.0+38971237.0;
        samplePOT = (EXT/9764047.0)*4.4e+19;
      }

      int nueN = 10;
      //float nueL = -0.2;
      nueH = 2.0;
      float nue_edges[nueN+1] = {-0.2,0.0,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0};

      int numuN = 20;
      float numuL = -0.1;
      numuH = 1.9;

      int teN = 40;
      float teL = 0.;
      float teH = 4.;

      if(sampleType == "nueoverlay"){
        h_nue_allTrue = new TH1F(("h_nue_allTrue_"+sampleType).c_str(),"",nueN,nue_edges);
        h_nue_pTrue_foundVtx = new TH1F(("h_nue_pTrue_foundVtx_"+sampleType).c_str(),"",nueN,nue_edges);
        h_nue_pTrue_foundNuVtx = new TH1F(("h_nue_pTrue_foundNuVtx_"+sampleType).c_str(),"",nueN,nue_edges);
        h_nue_pTrue_noMu = new TH1F(("h_nue_pTrue_noMu_"+sampleType).c_str(),"",nueN,nue_edges);
        h_nue_pTrue_noMuLike = new TH1F(("h_nue_pTrue_noMuLike_"+sampleType).c_str(),"",nueN,nue_edges);
        h_nue_pTrue_primEl = new TH1F(("h_nue_pTrue_primEl_"+sampleType).c_str(),"",nueN,nue_edges);
        h_nue_pTrue_confPrimEl = new TH1F(("h_nue_pTrue_confPrimEl_"+sampleType).c_str(),"",nueN,nue_edges);
        h_nue_pTrue_allCuts = new TH1F(("h_nue_pTrue_allCuts_"+sampleType).c_str(),"",nueN,nue_edges);
        hT_nue_allTrue = new TH1F(("hT_nue_allTrue_"+sampleType).c_str(),"",teN,teL,teH);
        hT_nue_pTrue_foundVtx = new TH1F(("hT_nue_pTrue_foundVtx_"+sampleType).c_str(),"",teN,teL,teH);
        hT_nue_pTrue_foundNuVtx = new TH1F(("hT_nue_pTrue_foundNuVtx_"+sampleType).c_str(),"",teN,teL,teH);
        hT_nue_pTrue_noMu = new TH1F(("hT_nue_pTrue_noMu_"+sampleType).c_str(),"",teN,teL,teH);
        hT_nue_pTrue_noMuLike = new TH1F(("hT_nue_pTrue_noMuLike_"+sampleType).c_str(),"",teN,teL,teH);
        hT_nue_pTrue_primEl = new TH1F(("hT_nue_pTrue_primEl_"+sampleType).c_str(),"",teN,teL,teH);
        hT_nue_pTrue_confPrimEl = new TH1F(("hT_nue_pTrue_confPrimEl_"+sampleType).c_str(),"",teN,teL,teH);
        hT_nue_pTrue_allCuts = new TH1F(("hT_nue_pTrue_allCuts_"+sampleType).c_str(),"",teN,teL,teH);
      }

      if(sampleType == "bnboverlay"){
        h_numu_allTrue = new TH1F(("h_numu_allTrue_"+sampleType).c_str(),"",numuN,numuL,numuH);
        h_numu_pTrue_foundVtx = new TH1F(("h_numu_pTrue_foundVtx_"+sampleType).c_str(),"",numuN,numuL,numuH);
        h_numu_pTrue_foundNuVtx = new TH1F(("h_numu_pTrue_foundNuVtx_"+sampleType).c_str(),"",numuN,numuL,numuH);
        h_numu_pTrue_allCuts = new TH1F(("h_numu_pTrue_allCuts_"+sampleType).c_str(),"",numuN,numuL,numuH);
        hT_numu_allTrue = new TH1F(("hT_numu_allTrue_"+sampleType).c_str(),"",teN,teL,teH);
        hT_numu_pTrue_foundVtx = new TH1F(("hT_numu_pTrue_foundVtx_"+sampleType).c_str(),"",teN,teL,teH);
        hT_numu_pTrue_foundNuVtx = new TH1F(("hT_numu_pTrue_foundNuVtx_"+sampleType).c_str(),"",teN,teL,teH);
        hT_numu_pTrue_allCuts = new TH1F(("hT_numu_pTrue_allCuts_"+sampleType).c_str(),"",teN,teL,teH);
      }

      h_nue_pAll_foundVtx = new TH1F(("h_nue_pAll_foundVtx_"+sampleType).c_str(),"",nueN,nue_edges);
      h_nue_pAll_foundNuVtx = new TH1F(("h_nue_pAll_foundNuVtx_"+sampleType).c_str(),"",nueN,nue_edges);
      h_nue_pAll_noMu = new TH1F(("h_nue_pAll_noMu_"+sampleType).c_str(),"",nueN,nue_edges);
      h_nue_pAll_noMuLike = new TH1F(("h_nue_pAll_noMuLike_"+sampleType).c_str(),"",nueN,nue_edges);
      h_nue_pAll_primEl = new TH1F(("h_nue_pAll_primEl_"+sampleType).c_str(),"",nueN,nue_edges);
      h_nue_pAll_confPrimEl = new TH1F(("h_nue_pAll_confPrimEl_"+sampleType).c_str(),"",nueN,nue_edges);
      h_nue_pAll_allCuts = new TH1F(("h_nue_pAll_allCuts_"+sampleType).c_str(),"",nueN,nue_edges);

      h_numu_pAll_foundVtx = new TH1F(("h_numu_pAll_foundVtx_"+sampleType).c_str(),"",numuN,numuL,numuH);
      h_numu_pAll_foundNuVtx = new TH1F(("h_numu_pAll_foundNuVtx_"+sampleType).c_str(),"",numuN,numuL,numuH);
      h_numu_pAll_allCuts = new TH1F(("h_numu_pAll_allCuts_"+sampleType).c_str(),"",numuN,numuL,numuH);

    }

    /*
    ~HistFiller(){
      delete ntuple;
      delete h_nue_allTrue;
      delete h_nue_pTrue_foundVtx;
      delete h_nue_pTrue_foundNuVtx;
      delete h_nue_pTrue_noMu;
      delete h_nue_pTrue_noMuLike;
      delete h_nue_pTrue_primEl;
      delete h_nue_pTrue_confPrimEl;
      delete h_nue_pTrue_allCuts;
      delete h_numu_allTrue;
      delete h_numu_pTrue_foundVtx;
      delete h_numu_pTrue_foundNuVtx;
      delete h_numu_pTrue_allCuts;
      delete h_nue_pAll_foundVtx;
      delete h_nue_pAll_foundNuVtx;
      delete h_nue_pAll_noMu;
      delete h_nue_pAll_noMuLike;
      delete h_nue_pAll_primEl;
      delete h_nue_pAll_confPrimEl;
      delete h_nue_pAll_allCuts;
      delete h_numu_pAll_foundVtx;
      delete h_numu_pAll_foundNuVtx;
      delete h_numu_pAll_allCuts;
    }*/


    void fillHists(){

      TH1::SetDefaultSumw2(kTRUE);

      for(int i = 0; i < eventTree->GetEntries(); ++i){

        eventTree -> GetEntry(i);

        float recoNuE_nue = recoNuE/1000.;
        float recoNuE_numu = recoNuE/1000.;
        if(recoNuE/1000. >= nueH) recoNuE_nue = nueH - 1e-6;
        if(recoNuE/1000. >= numuH) recoNuE_numu = numuH - 1e-6;

        float weight = 1.;

        bool trueCCnue = false;
        if(sampleType == "nueoverlay"){
          if(std::abs(trueNuPDG) != 12 || trueNuCCNC != 0 || std::isinf(xsecWeight)) continue;
          trueCCnue = true;
          weight = xsecWeight;
          h_nue_allTrue -> Fill(recoNuE_nue, weight*(targetPOT/samplePOT));
          hT_nue_allTrue -> Fill(trueNuE, weight*(targetPOT/samplePOT));
        }

        bool trueCCnumu = false;
        if(sampleType == "bnboverlay"){
          if((std::abs(trueNuPDG) == 12 && trueNuCCNC == 0) || std::isinf(xsecWeight)) continue;
          weight = xsecWeight;
          if(std::abs(trueNuPDG) == 14 && trueNuCCNC == 0){
            trueCCnumu = true;
            h_numu_allTrue -> Fill(recoNuE_numu, weight*(targetPOT/samplePOT));
            hT_numu_allTrue -> Fill(trueNuE, weight*(targetPOT/samplePOT));
          }
        }

        SelResults results = checkCuts(foundVertex, vtxIsFiducial, vtxFracHitsOnCosmic,
         nTracks, nShowers, trackIsSecondary, trackClassified, showerIsSecondary, showerClassified,
         trackPID, showerPID, showerProcess, trackMuScore, showerCharge, showerElScore,
         showerPhScore, showerPiScore);

        if(results.foundVtx){
          h_nue_pAll_foundVtx -> Fill(recoNuE_nue, weight*(targetPOT/samplePOT));
          h_numu_pAll_foundVtx -> Fill(recoNuE_numu, weight*(targetPOT/samplePOT));
          if(trueCCnue){
            h_nue_pTrue_foundVtx -> Fill(recoNuE_nue, weight*(targetPOT/samplePOT));
            hT_nue_pTrue_foundVtx -> Fill(trueNuE, weight*(targetPOT/samplePOT));
          }
          if(trueCCnumu){
            h_numu_pTrue_foundVtx -> Fill(recoNuE_numu, weight*(targetPOT/samplePOT));
            hT_numu_pTrue_foundVtx -> Fill(trueNuE, weight*(targetPOT/samplePOT));
          }
        }

        if(results.foundNuVtx){
          h_nue_pAll_foundNuVtx -> Fill(recoNuE_nue, weight*(targetPOT/samplePOT));
          h_numu_pAll_foundNuVtx -> Fill(recoNuE_numu, weight*(targetPOT/samplePOT));
          if(trueCCnue){
            h_nue_pTrue_foundNuVtx -> Fill(recoNuE_nue, weight*(targetPOT/samplePOT));
            hT_nue_pTrue_foundNuVtx -> Fill(trueNuE, weight*(targetPOT/samplePOT));
          }
          if(trueCCnumu){
            h_numu_pTrue_foundNuVtx -> Fill(recoNuE_numu, weight*(targetPOT/samplePOT));
            hT_numu_pTrue_foundNuVtx -> Fill(trueNuE, weight*(targetPOT/samplePOT));
          }
        }

        if(results.noMu){
          h_nue_pAll_noMu -> Fill(recoNuE_nue, weight*(targetPOT/samplePOT));
          if(trueCCnue){
            h_nue_pTrue_noMu -> Fill(recoNuE_nue, weight*(targetPOT/samplePOT));
            hT_nue_pTrue_noMu -> Fill(trueNuE, weight*(targetPOT/samplePOT));
          }
        }

        if(results.primEl){
          h_nue_pAll_primEl -> Fill(recoNuE_nue, weight*(targetPOT/samplePOT));
          if(trueCCnue){
            h_nue_pTrue_primEl -> Fill(recoNuE_nue, weight*(targetPOT/samplePOT));
            hT_nue_pTrue_primEl -> Fill(trueNuE, weight*(targetPOT/samplePOT));
          }
        }

        if(results.noMuLike){
          h_nue_pAll_noMuLike -> Fill(recoNuE_nue, weight*(targetPOT/samplePOT));
          if(trueCCnue){
            h_nue_pTrue_noMuLike -> Fill(recoNuE_nue, weight*(targetPOT/samplePOT));
            hT_nue_pTrue_noMuLike -> Fill(trueNuE, weight*(targetPOT/samplePOT));
          }
        }

        if(results.confPrimEl){
          h_nue_pAll_confPrimEl -> Fill(recoNuE_nue, weight*(targetPOT/samplePOT));
          if(trueCCnue){
            h_nue_pTrue_confPrimEl -> Fill(recoNuE_nue, weight*(targetPOT/samplePOT));
            hT_nue_pTrue_confPrimEl -> Fill(trueNuE, weight*(targetPOT/samplePOT));
          }
        }

        if(results.passedCCnue){
          h_nue_pAll_allCuts -> Fill(recoNuE_nue, weight*(targetPOT/samplePOT));
          if(trueCCnue){
            h_nue_pTrue_allCuts -> Fill(recoNuE_nue, weight*(targetPOT/samplePOT));
            hT_nue_pTrue_allCuts -> Fill(trueNuE, weight*(targetPOT/samplePOT));
          }
        }

        if(results.passedCCnumu){
          h_numu_pAll_allCuts -> Fill(recoNuE_numu, weight*(targetPOT/samplePOT));
          if(trueCCnumu){
            h_numu_pTrue_allCuts -> Fill(recoNuE_numu, weight*(targetPOT/samplePOT));
            hT_numu_pTrue_allCuts -> Fill(trueNuE, weight*(targetPOT/samplePOT));
          }
        }

      }

    }

};




int main(int argc, char* argv[]) {

  TH1::SetDefaultSumw2(kTRUE);

  std::cout << "initializing histogram filler classes" << std::endl;

  HistFiller nue_results("nueoverlay");
  HistFiller bnb_results("bnboverlay");
  HistFiller allCos_results("allCosmics");
  HistFiller run1Cos_results("run1Cosmics");

  std::cout << "looping over intrinsic nue overlay file" << std::endl;
  nue_results.fillHists();
  std::cout << "looping over bnb nu overlay file" << std::endl;
  bnb_results.fillHists();
  std::cout << "looping over runs1-3 extbnb file" << std::endl;
  allCos_results.fillHists();
  std::cout << "looping over run1 extbnb file" << std::endl;
  run1Cos_results.fillHists();

  std::cout << "making purity and efficiency histograms" << std::endl;

  TH1F* h_nue_eff_foundVtx = (TH1F*)nue_results.h_nue_pTrue_foundVtx->Clone("h_nue_eff_foundVtx");
  h_nue_eff_foundVtx -> Divide(h_nue_eff_foundVtx,nue_results.h_nue_allTrue,1,1,"B");

  TH1F* h_nue_eff_foundNuVtx = (TH1F*)nue_results.h_nue_pTrue_foundNuVtx->Clone("h_nue_eff_foundNuVtx");
  h_nue_eff_foundNuVtx -> Divide(h_nue_eff_foundNuVtx,nue_results.h_nue_allTrue,1,1,"B");

  TH1F* h_nue_eff_noMu = (TH1F*)nue_results.h_nue_pTrue_noMu->Clone("h_nue_eff_noMu");
  h_nue_eff_noMu -> Divide(h_nue_eff_noMu,nue_results.h_nue_allTrue,1,1,"B");

  TH1F* h_nue_eff_noMuLike = (TH1F*)nue_results.h_nue_pTrue_noMuLike->Clone("h_nue_eff_noMuLike");
  h_nue_eff_noMuLike -> Divide(h_nue_eff_noMuLike,nue_results.h_nue_allTrue,1,1,"B");

  TH1F* h_nue_eff_primEl = (TH1F*)nue_results.h_nue_pTrue_primEl->Clone("h_nue_eff_primEl");
  h_nue_eff_primEl -> Divide(h_nue_eff_primEl,nue_results.h_nue_allTrue,1,1,"B");

  TH1F* h_nue_eff_confPrimEl = (TH1F*)nue_results.h_nue_pTrue_confPrimEl->Clone("h_nue_eff_confPrimEl");
  h_nue_eff_confPrimEl -> Divide(h_nue_eff_confPrimEl,nue_results.h_nue_allTrue,1,1,"B");

  TH1F* h_nue_eff_allCuts = (TH1F*)nue_results.h_nue_pTrue_allCuts->Clone("h_nue_eff_allCuts");
  h_nue_eff_allCuts -> Divide(h_nue_eff_allCuts,nue_results.h_nue_allTrue,1,1,"B");

  TH1F* h_numu_eff_foundVtx = (TH1F*)bnb_results.h_numu_pTrue_foundVtx->Clone("h_numu_eff_foundVtx");
  h_numu_eff_foundVtx -> Divide(h_numu_eff_foundVtx,bnb_results.h_numu_allTrue,1,1,"B");

  TH1F* h_numu_eff_foundNuVtx = (TH1F*)bnb_results.h_numu_pTrue_foundNuVtx->Clone("h_numu_eff_foundNuVtx");
  h_numu_eff_foundNuVtx -> Divide(h_numu_eff_foundNuVtx,bnb_results.h_numu_allTrue,1,1,"B");

  TH1F* h_numu_eff_allCuts = (TH1F*)bnb_results.h_numu_pTrue_allCuts->Clone("h_numu_eff_allCuts");
  h_numu_eff_allCuts -> Divide(h_numu_eff_allCuts,bnb_results.h_numu_allTrue,1,1,"B");

  TH1F* hT_nue_eff_foundVtx = (TH1F*)nue_results.hT_nue_pTrue_foundVtx->Clone("hT_nue_eff_foundVtx");
  hT_nue_eff_foundVtx -> Divide(hT_nue_eff_foundVtx,nue_results.hT_nue_allTrue,1,1,"B");

  TH1F* hT_nue_eff_foundNuVtx = (TH1F*)nue_results.hT_nue_pTrue_foundNuVtx->Clone("hT_nue_eff_foundNuVtx");
  hT_nue_eff_foundNuVtx -> Divide(hT_nue_eff_foundNuVtx,nue_results.hT_nue_allTrue,1,1,"B");

  TH1F* hT_nue_eff_noMu = (TH1F*)nue_results.hT_nue_pTrue_noMu->Clone("hT_nue_eff_noMu");
  hT_nue_eff_noMu -> Divide(hT_nue_eff_noMu,nue_results.hT_nue_allTrue,1,1,"B");

  TH1F* hT_nue_eff_noMuLike = (TH1F*)nue_results.hT_nue_pTrue_noMuLike->Clone("hT_nue_eff_noMuLike");
  hT_nue_eff_noMuLike -> Divide(hT_nue_eff_noMuLike,nue_results.hT_nue_allTrue,1,1,"B");

  TH1F* hT_nue_eff_primEl = (TH1F*)nue_results.hT_nue_pTrue_primEl->Clone("hT_nue_eff_primEl");
  hT_nue_eff_primEl -> Divide(hT_nue_eff_primEl,nue_results.hT_nue_allTrue,1,1,"B");

  TH1F* hT_nue_eff_confPrimEl = (TH1F*)nue_results.hT_nue_pTrue_confPrimEl->Clone("hT_nue_eff_confPrimEl");
  hT_nue_eff_confPrimEl -> Divide(hT_nue_eff_confPrimEl,nue_results.hT_nue_allTrue,1,1,"B");

  TH1F* hT_nue_eff_allCuts = (TH1F*)nue_results.hT_nue_pTrue_allCuts->Clone("hT_nue_eff_allCuts");
  hT_nue_eff_allCuts -> Divide(hT_nue_eff_allCuts,nue_results.hT_nue_allTrue,1,1,"B");

  TH1F* hT_numu_eff_foundVtx = (TH1F*)bnb_results.hT_numu_pTrue_foundVtx->Clone("hT_numu_eff_foundVtx");
  hT_numu_eff_foundVtx -> Divide(hT_numu_eff_foundVtx,bnb_results.hT_numu_allTrue,1,1,"B");

  TH1F* hT_numu_eff_foundNuVtx = (TH1F*)bnb_results.hT_numu_pTrue_foundNuVtx->Clone("hT_numu_eff_foundNuVtx");
  hT_numu_eff_foundNuVtx -> Divide(hT_numu_eff_foundNuVtx,bnb_results.hT_numu_allTrue,1,1,"B");

  TH1F* hT_numu_eff_allCuts = (TH1F*)bnb_results.hT_numu_pTrue_allCuts->Clone("hT_numu_eff_allCuts");
  hT_numu_eff_allCuts -> Divide(hT_numu_eff_allCuts,bnb_results.hT_numu_allTrue,1,1,"B");

  TH1F* h_nue_all_foundVtx = (TH1F*)nue_results.h_nue_pAll_foundVtx->Clone("h_nue_all_foundVtx");
  h_nue_all_foundVtx->Add(bnb_results.h_nue_pAll_foundVtx);
  h_nue_all_foundVtx->Add(allCos_results.h_nue_pAll_foundVtx);
  TH1F* h_nue_pur_foundVtx = (TH1F*)nue_results.h_nue_pTrue_foundVtx->Clone("h_nue_pur_foundVtx");
  h_nue_pur_foundVtx->Divide(h_nue_all_foundVtx);

  TH1F* h_nue_all_foundNuVtx = (TH1F*)nue_results.h_nue_pAll_foundNuVtx->Clone("h_nue_all_foundNuVtx");
  h_nue_all_foundNuVtx->Add(bnb_results.h_nue_pAll_foundNuVtx);
  h_nue_all_foundNuVtx->Add(allCos_results.h_nue_pAll_foundNuVtx);
  TH1F* h_nue_pur_foundNuVtx = (TH1F*)nue_results.h_nue_pTrue_foundNuVtx->Clone("h_nue_pur_foundNuVtx");
  h_nue_pur_foundNuVtx->Divide(h_nue_all_foundNuVtx);

  TH1F* h_nue_all_noMu = (TH1F*)nue_results.h_nue_pAll_noMu->Clone("h_nue_all_noMu");
  h_nue_all_noMu->Add(bnb_results.h_nue_pAll_noMu);
  h_nue_all_noMu->Add(allCos_results.h_nue_pAll_noMu);
  TH1F* h_nue_pur_noMu = (TH1F*)nue_results.h_nue_pTrue_noMu->Clone("h_nue_pur_noMu");
  h_nue_pur_noMu->Divide(h_nue_all_noMu);

  TH1F* h_nue_all_noMuLike = (TH1F*)nue_results.h_nue_pAll_noMuLike->Clone("h_nue_all_noMuLike");
  h_nue_all_noMuLike->Add(bnb_results.h_nue_pAll_noMuLike);
  h_nue_all_noMuLike->Add(allCos_results.h_nue_pAll_noMuLike);
  TH1F* h_nue_pur_noMuLike = (TH1F*)nue_results.h_nue_pTrue_noMuLike->Clone("h_nue_pur_noMuLike");
  h_nue_pur_noMuLike->Divide(h_nue_all_noMuLike);

  TH1F* h_nue_all_primEl = (TH1F*)nue_results.h_nue_pAll_primEl->Clone("h_nue_all_primEl");
  h_nue_all_primEl->Add(bnb_results.h_nue_pAll_primEl);
  h_nue_all_primEl->Add(allCos_results.h_nue_pAll_primEl);
  TH1F* h_nue_pur_primEl = (TH1F*)nue_results.h_nue_pTrue_primEl->Clone("h_nue_pur_primEl");
  h_nue_pur_primEl->Divide(h_nue_all_primEl);

  TH1F* h_nue_all_confPrimEl = (TH1F*)nue_results.h_nue_pAll_confPrimEl->Clone("h_nue_all_confPrimEl");
  h_nue_all_confPrimEl->Add(bnb_results.h_nue_pAll_confPrimEl);
  h_nue_all_confPrimEl->Add(allCos_results.h_nue_pAll_confPrimEl);
  TH1F* h_nue_pur_confPrimEl = (TH1F*)nue_results.h_nue_pTrue_confPrimEl->Clone("h_nue_pur_confPrimEl");
  h_nue_pur_confPrimEl->Divide(h_nue_all_confPrimEl);

  TH1F* h_nue_all_allCuts = (TH1F*)nue_results.h_nue_pAll_allCuts->Clone("h_nue_all_allCuts");
  h_nue_all_allCuts->Add(bnb_results.h_nue_pAll_allCuts);
  h_nue_all_allCuts->Add(allCos_results.h_nue_pAll_allCuts);
  TH1F* h_nue_pur_allCuts = (TH1F*)nue_results.h_nue_pTrue_allCuts->Clone("h_nue_pur_allCuts");
  h_nue_pur_allCuts->Divide(h_nue_all_allCuts);

  TH1F* h_numu_all_foundVtx = (TH1F*)nue_results.h_numu_pAll_foundVtx->Clone("h_numu_all_foundVtx");
  h_numu_all_foundVtx->Add(bnb_results.h_numu_pAll_foundVtx);
  h_numu_all_foundVtx->Add(run1Cos_results.h_numu_pAll_foundVtx);
  TH1F* h_numu_pur_foundVtx = (TH1F*)bnb_results.h_numu_pTrue_foundVtx->Clone("h_numu_pur_foundVtx");
  h_numu_pur_foundVtx->Divide(h_numu_all_foundVtx);

  TH1F* h_numu_all_foundNuVtx = (TH1F*)nue_results.h_numu_pAll_foundNuVtx->Clone("h_numu_all_foundNuVtx");
  h_numu_all_foundNuVtx->Add(bnb_results.h_numu_pAll_foundNuVtx);
  h_numu_all_foundNuVtx->Add(run1Cos_results.h_numu_pAll_foundNuVtx);
  TH1F* h_numu_pur_foundNuVtx = (TH1F*)bnb_results.h_numu_pTrue_foundNuVtx->Clone("h_numu_pur_foundNuVtx");
  h_numu_pur_foundNuVtx->Divide(h_numu_all_foundNuVtx);

  TH1F* h_numu_all_allCuts = (TH1F*)nue_results.h_numu_pAll_allCuts->Clone("h_numu_all_allCuts");
  h_numu_all_allCuts->Add(bnb_results.h_numu_pAll_allCuts);
  h_numu_all_allCuts->Add(run1Cos_results.h_numu_pAll_allCuts);
  TH1F* h_numu_pur_allCuts = (TH1F*)bnb_results.h_numu_pTrue_allCuts->Clone("h_numu_pur_allCuts");
  h_numu_pur_allCuts->Divide(h_numu_all_allCuts);


  TH1F* h_nue_sig_noCuts = (TH1F*)nue_results.h_nue_allTrue->Clone("h_nue_sig_noCuts");
  TH1F* h_nue_sig_foundVtx = (TH1F*)nue_results.h_nue_pTrue_foundVtx->Clone("h_nue_sig_foundVtx");
  TH1F* h_nue_sig_foundNuVtx = (TH1F*)nue_results.h_nue_pTrue_foundNuVtx->Clone("h_nue_sig_foundNuVtx");
  TH1F* h_nue_sig_noMu = (TH1F*)nue_results.h_nue_pTrue_noMu->Clone("h_nue_sig_noMu");
  TH1F* h_nue_sig_primEl = (TH1F*)nue_results.h_nue_pTrue_primEl->Clone("h_nue_sig_primEl");
  TH1F* h_nue_sig_noMuLike = (TH1F*)nue_results.h_nue_pTrue_noMuLike->Clone("h_nue_sig_noMuLike");
  TH1F* h_nue_sig_confPrimEl = (TH1F*)nue_results.h_nue_pTrue_confPrimEl->Clone("h_nue_sig_confPrimEl");
  TH1F* h_nue_sig_allCuts = (TH1F*)nue_results.h_nue_pTrue_allCuts->Clone("h_nue_sig_allCuts");

  TH1F* h_numu_sig_noCuts = (TH1F*)bnb_results.h_numu_allTrue->Clone("h_numu_sig_noCuts");
  TH1F* h_numu_sig_foundVtx = (TH1F*)bnb_results.h_numu_pTrue_foundVtx->Clone("h_numu_sig_foundVtx");
  TH1F* h_numu_sig_foundNuVtx = (TH1F*)bnb_results.h_numu_pTrue_foundNuVtx->Clone("h_numu_sig_foundNuVtx");
  TH1F* h_numu_sig_allCuts = (TH1F*)bnb_results.h_numu_pTrue_allCuts->Clone("h_numu_sig_allCuts");

  TH1F* h_nue_bkg_foundVtx = (TH1F*)h_nue_all_foundVtx->Clone("h_nue_bkg_foundVtx");
  h_nue_bkg_foundVtx -> Add(nue_results.h_nue_pTrue_foundVtx, -1);
  TH1F* h_nue_bkg_foundNuVtx = (TH1F*)h_nue_all_foundNuVtx->Clone("h_nue_bkg_foundNuVtx");
  h_nue_bkg_foundNuVtx -> Add(nue_results.h_nue_pTrue_foundNuVtx, -1);
  TH1F* h_nue_bkg_noMu = (TH1F*)h_nue_all_noMu->Clone("h_nue_bkg_noMu");
  h_nue_bkg_noMu -> Add(nue_results.h_nue_pTrue_noMu, -1);
  TH1F* h_nue_bkg_primEl = (TH1F*)h_nue_all_primEl->Clone("h_nue_bkg_primEl");
  h_nue_bkg_primEl -> Add(nue_results.h_nue_pTrue_primEl, -1);
  TH1F* h_nue_bkg_noMuLike = (TH1F*)h_nue_all_noMuLike->Clone("h_nue_bkg_noMuLike");
  h_nue_bkg_noMuLike -> Add(nue_results.h_nue_pTrue_noMuLike, -1);
  TH1F* h_nue_bkg_confPrimEl = (TH1F*)h_nue_all_confPrimEl->Clone("h_nue_bkg_confPrimEl");
  h_nue_bkg_confPrimEl -> Add(nue_results.h_nue_pTrue_confPrimEl, -1);
  TH1F* h_nue_bkg_allCuts = (TH1F*)h_nue_all_allCuts->Clone("h_nue_bkg_allCuts");
  h_nue_bkg_allCuts -> Add(nue_results.h_nue_pTrue_allCuts, -1);

  TH1F* h_numu_bkg_foundVtx = (TH1F*)h_numu_all_foundVtx->Clone("h_numu_bkg_foundVtx");
  h_numu_bkg_foundVtx -> Add(bnb_results.h_numu_pTrue_foundVtx, -1);
  TH1F* h_numu_bkg_foundNuVtx = (TH1F*)h_numu_all_foundNuVtx->Clone("h_numu_bkg_foundNuVtx");
  h_numu_bkg_foundNuVtx -> Add(bnb_results.h_numu_pTrue_foundNuVtx, -1);
  TH1F* h_numu_bkg_allCuts = (TH1F*)h_numu_all_allCuts->Clone("h_numu_bkg_allCuts");
  h_numu_bkg_allCuts -> Add(bnb_results.h_numu_pTrue_allCuts, -1);


  double CCnue_eff_foundVtx = nue_results.h_nue_pTrue_foundVtx->Integral() / nue_results.h_nue_allTrue->Integral();
  double CCnue_pur_foundVtx = nue_results.h_nue_pTrue_foundVtx->Integral() / h_nue_all_foundVtx->Integral();
  double CCnumu_eff_foundVtx = bnb_results.h_numu_pTrue_foundVtx->Integral() / bnb_results.h_numu_allTrue->Integral();
  double CCnumu_pur_foundVtx = bnb_results.h_numu_pTrue_foundVtx->Integral() / h_numu_all_foundVtx->Integral();

  double CCnue_eff_foundNuVtx = nue_results.h_nue_pTrue_foundNuVtx->Integral() / nue_results.h_nue_allTrue->Integral();
  double CCnue_pur_foundNuVtx = nue_results.h_nue_pTrue_foundNuVtx->Integral() / h_nue_all_foundNuVtx->Integral();
  double CCnumu_eff_foundNuVtx = bnb_results.h_numu_pTrue_foundNuVtx->Integral() / bnb_results.h_numu_allTrue->Integral();
  double CCnumu_pur_foundNuVtx = bnb_results.h_numu_pTrue_foundNuVtx->Integral() / h_numu_all_foundNuVtx->Integral();

  double CCnue_eff_noMu = nue_results.h_nue_pTrue_noMu->Integral() / nue_results.h_nue_allTrue->Integral();
  double CCnue_pur_noMu = nue_results.h_nue_pTrue_noMu->Integral() / h_nue_all_noMu->Integral();

  double CCnue_eff_primEl = nue_results.h_nue_pTrue_primEl->Integral() / nue_results.h_nue_allTrue->Integral();
  double CCnue_pur_primEl = nue_results.h_nue_pTrue_primEl->Integral() / h_nue_all_primEl->Integral();

  double CCnue_eff_noMuLike = nue_results.h_nue_pTrue_noMuLike->Integral() / nue_results.h_nue_allTrue->Integral();
  double CCnue_pur_noMuLike = nue_results.h_nue_pTrue_noMuLike->Integral() / h_nue_all_noMuLike->Integral();

  double CCnue_eff_confPrimEl = nue_results.h_nue_pTrue_confPrimEl->Integral() / nue_results.h_nue_allTrue->Integral();
  double CCnue_pur_confPrimEl = nue_results.h_nue_pTrue_confPrimEl->Integral() / h_nue_all_confPrimEl->Integral();

  double CCnue_eff_allCuts = nue_results.h_nue_pTrue_allCuts->Integral() / nue_results.h_nue_allTrue->Integral();
  double CCnue_pur_allCuts = nue_results.h_nue_pTrue_allCuts->Integral() / h_nue_all_allCuts->Integral();
  double CCnumu_eff_allCuts = bnb_results.h_numu_pTrue_allCuts->Integral() / bnb_results.h_numu_allTrue->Integral();
  double CCnumu_pur_allCuts = bnb_results.h_numu_pTrue_allCuts->Integral() / h_numu_all_allCuts->Integral();

  std::cout << std::endl;
  std::cout << "noCuts CCnue signal count: " << h_nue_sig_noCuts->Integral() << std::endl;
  std::cout << std::endl;

  std::cout << "foundVtx CCnue efficiency: " << CCnue_eff_foundVtx << std::endl;
  std::cout << "foundVtx CCnue purity: " << CCnue_pur_foundVtx << std::endl;
  std::cout << "foundVtx CCnue signal count: " << h_nue_sig_foundVtx->Integral() << std::endl;
  std::cout << "foundVtx CCnue total background count: " << h_nue_bkg_foundVtx->Integral() << std::endl;
  std::cout << "foundVtx CCnue cosmic background count: " <<allCos_results.h_nue_pAll_foundVtx->Integral()<<std::endl;
  std::cout << std::endl;

  std::cout << "foundNuVtx CCnue efficiency: " << CCnue_eff_foundNuVtx << std::endl;
  std::cout << "foundNuVtx CCnue purity: " << CCnue_pur_foundNuVtx << std::endl;
  std::cout << "foundNuVtx CCnue signal count: " << h_nue_sig_foundNuVtx->Integral() << std::endl;
  std::cout << "foundNuVtx CCnue total background count: " << h_nue_bkg_foundNuVtx->Integral() << std::endl;
  std::cout << "foundNuVtx CCnue cosmic background count: " <<allCos_results.h_nue_pAll_foundNuVtx->Integral()<<std::endl;
  std::cout << std::endl;

  std::cout << "noMu CCnue efficiency: " << CCnue_eff_noMu << std::endl;
  std::cout << "noMu CCnue purity: " << CCnue_pur_noMu << std::endl;
  std::cout << "noMu CCnue signal count: " << h_nue_sig_noMu->Integral() << std::endl;
  std::cout << "noMu CCnue total background count: " << h_nue_bkg_noMu->Integral() << std::endl;
  std::cout << "noMu CCnue cosmic background count: " <<allCos_results.h_nue_pAll_noMu->Integral()<<std::endl;
  std::cout << std::endl;

  std::cout << "primEl CCnue efficiency: " << CCnue_eff_primEl << std::endl;
  std::cout << "primEl CCnue purity: " << CCnue_pur_primEl << std::endl;
  std::cout << "primEl CCnue signal count: " << h_nue_sig_primEl->Integral() << std::endl;
  std::cout << "primEl CCnue total background count: " << h_nue_bkg_primEl->Integral() << std::endl;
  std::cout << "primEl CCnue cosmic background count: " <<allCos_results.h_nue_pAll_primEl->Integral()<<std::endl;
  std::cout << std::endl;

  std::cout << "noMuLike CCnue efficiency: " << CCnue_eff_noMuLike << std::endl;
  std::cout << "noMuLike CCnue purity: " << CCnue_pur_noMuLike << std::endl;
  std::cout << "noMuLike CCnue signal count: " << h_nue_sig_noMuLike->Integral() << std::endl;
  std::cout << "noMuLike CCnue total background count: " << h_nue_bkg_noMuLike->Integral() << std::endl;
  std::cout << "noMuLike CCnue cosmic background count: " <<allCos_results.h_nue_pAll_noMuLike->Integral()<<std::endl;
  std::cout << std::endl;

  std::cout << "confPrimEl CCnue efficiency: " << CCnue_eff_confPrimEl << std::endl;
  std::cout << "confPrimEl CCnue purity: " << CCnue_pur_confPrimEl << std::endl;
  std::cout << "confPrimEl CCnue signal count: " << h_nue_sig_confPrimEl->Integral() << std::endl;
  std::cout << "confPrimEl CCnue total background count: " << h_nue_bkg_confPrimEl->Integral() << std::endl;
  std::cout << "confPrimEl CCnue cosmic background count: " <<allCos_results.h_nue_pAll_confPrimEl->Integral()<<std::endl;
  std::cout << std::endl;

  std::cout << "allCuts CCnue efficiency: " << CCnue_eff_allCuts << std::endl;
  std::cout << "allCuts CCnue purity: " << CCnue_pur_allCuts << std::endl;
  std::cout << "allCuts CCnue signal count: " << h_nue_sig_allCuts->Integral() << std::endl;
  std::cout << "allCuts CCnue total background count: " << h_nue_bkg_allCuts->Integral() << std::endl;
  std::cout << "allCuts CCnue cosmic background count: " <<allCos_results.h_nue_pAll_allCuts->Integral()<<std::endl;
  std::cout << std::endl;

  std::cout << "noCuts CCnumu signal count: " << h_numu_sig_noCuts->Integral() << std::endl;
  std::cout << std::endl;

  std::cout << "foundVtx CCnumu efficiency: " << CCnumu_eff_foundVtx << std::endl;
  std::cout << "foundVtx CCnumu purity: " << CCnumu_pur_foundVtx << std::endl;
  std::cout << "foundVtx CCnumu signal count: " << h_numu_sig_foundVtx->Integral() << std::endl;
  std::cout << "foundVtx CCnumu total background count: " << h_numu_bkg_foundVtx->Integral() << std::endl;
  std::cout << "foundVtx CCnumu cosmic background count: " <<run1Cos_results.h_numu_pAll_foundVtx->Integral()<<std::endl;
  std::cout << std::endl;

  std::cout << "foundNuVtx CCnumu efficiency: " << CCnumu_eff_foundNuVtx << std::endl;
  std::cout << "foundNuVtx CCnumu purity: " << CCnumu_pur_foundNuVtx << std::endl;
  std::cout << "foundNuVtx CCnumu signal count: " << h_numu_sig_foundNuVtx->Integral() << std::endl;
  std::cout << "foundNuVtx CCnumu total background count: " << h_numu_bkg_foundNuVtx->Integral() << std::endl;
  std::cout << "foundNuVtx CCnumu cosmic background count: " <<run1Cos_results.h_numu_pAll_foundNuVtx->Integral()<<std::endl;
  std::cout << std::endl;

  std::cout << "allCuts CCnumu efficiency: " << CCnumu_eff_allCuts << std::endl;
  std::cout << "allCuts CCnumu purity: " << CCnumu_pur_allCuts << std::endl;
  std::cout << "allCuts CCnumu signal count: " << h_numu_sig_allCuts->Integral() << std::endl;
  std::cout << "allCuts CCnumu total background count: " << h_numu_bkg_allCuts->Integral() << std::endl;
  std::cout << "allCuts CCnumu cosmic background count: " <<run1Cos_results.h_numu_pAll_allCuts->Integral()<<std::endl;
  std::cout << std::endl;


  TFile* fout = new TFile("plot_CCinclusive_selections_efficiency_and_purity_makeHists_output.root", "RECREATE");
  fout -> cd();
  h_nue_eff_foundVtx->Write("", TObject::kOverwrite);
  hT_nue_eff_foundVtx->Write("", TObject::kOverwrite);
  h_nue_pur_foundVtx->Write("", TObject::kOverwrite);
  h_nue_eff_foundNuVtx->Write("", TObject::kOverwrite);
  hT_nue_eff_foundNuVtx->Write("", TObject::kOverwrite);
  h_nue_pur_foundNuVtx->Write("", TObject::kOverwrite);
  h_nue_eff_noMu->Write("", TObject::kOverwrite);
  hT_nue_eff_noMu->Write("", TObject::kOverwrite);
  h_nue_pur_noMu->Write("", TObject::kOverwrite);
  h_nue_eff_primEl->Write("", TObject::kOverwrite);
  hT_nue_eff_primEl->Write("", TObject::kOverwrite);
  h_nue_pur_primEl->Write("", TObject::kOverwrite);
  h_nue_eff_noMuLike->Write("", TObject::kOverwrite);
  hT_nue_eff_noMuLike->Write("", TObject::kOverwrite);
  h_nue_pur_noMuLike->Write("", TObject::kOverwrite);
  h_nue_eff_confPrimEl->Write("", TObject::kOverwrite);
  hT_nue_eff_confPrimEl->Write("", TObject::kOverwrite);
  h_nue_pur_confPrimEl->Write("", TObject::kOverwrite);
  h_nue_eff_allCuts->Write("", TObject::kOverwrite);
  hT_nue_eff_allCuts->Write("", TObject::kOverwrite);
  h_nue_pur_allCuts->Write("", TObject::kOverwrite);
  h_numu_eff_foundVtx->Write("", TObject::kOverwrite);
  hT_numu_eff_foundVtx->Write("", TObject::kOverwrite);
  h_numu_pur_foundVtx->Write("", TObject::kOverwrite);
  h_numu_eff_foundNuVtx->Write("", TObject::kOverwrite);
  hT_numu_eff_foundNuVtx->Write("", TObject::kOverwrite);
  h_numu_pur_foundNuVtx->Write("", TObject::kOverwrite);
  h_numu_eff_allCuts->Write("", TObject::kOverwrite);
  hT_numu_eff_allCuts->Write("", TObject::kOverwrite);
  h_numu_pur_allCuts->Write("", TObject::kOverwrite);
  h_nue_sig_noCuts->Write("", TObject::kOverwrite);
  h_nue_sig_foundVtx->Write("", TObject::kOverwrite);
  h_nue_sig_foundNuVtx->Write("", TObject::kOverwrite);
  h_nue_sig_noMu->Write("", TObject::kOverwrite);
  h_nue_sig_primEl->Write("", TObject::kOverwrite);
  h_nue_sig_noMuLike->Write("", TObject::kOverwrite);
  h_nue_sig_confPrimEl->Write("", TObject::kOverwrite);
  h_nue_sig_allCuts->Write("", TObject::kOverwrite);
  h_nue_bkg_foundVtx->Write("", TObject::kOverwrite);
  h_nue_bkg_foundNuVtx->Write("", TObject::kOverwrite);
  h_nue_bkg_noMu->Write("", TObject::kOverwrite);
  h_nue_bkg_primEl->Write("", TObject::kOverwrite);
  h_nue_bkg_noMuLike->Write("", TObject::kOverwrite);
  h_nue_bkg_confPrimEl->Write("", TObject::kOverwrite);
  h_nue_bkg_allCuts->Write("", TObject::kOverwrite);
  h_numu_sig_noCuts->Write("", TObject::kOverwrite);
  h_numu_sig_foundVtx->Write("", TObject::kOverwrite);
  h_numu_sig_foundNuVtx->Write("", TObject::kOverwrite);
  h_numu_sig_allCuts->Write("", TObject::kOverwrite);
  h_numu_bkg_foundVtx->Write("", TObject::kOverwrite);
  h_numu_bkg_foundNuVtx->Write("", TObject::kOverwrite);
  h_numu_bkg_allCuts->Write("", TObject::kOverwrite);
  fout -> Close();

}

