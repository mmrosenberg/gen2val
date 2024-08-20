
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <unistd.h>
#include <getopt.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TObject.h>


struct recoLepton {
  float momentum, cosTheta;
  recoLepton(float p, float c) : momentum(p), cosTheta(c) {}
  recoLepton() : momentum(-99.), cosTheta(-9.) {}
};

struct SelEvent {

  int run, subrun, event;
  float recoNuE, trueNuE;
  bool passedCCnue, passedCCnumu;
  recoLepton electron, muon;
  bool matched;
  float xsecWeight;

  SelEvent(int r, int sr, int e, float te, float re, float w, bool pe, bool pm, bool m) : run(r), subrun(sr), event(e), trueNuE(te), recoNuE(re), xsecWeight(w), passedCCnue(pe), passedCCnumu(pm), matched(m) {}
  SelEvent() : run(-1), subrun(-1), event(-1), trueNuE(-99.), recoNuE(-99.), xsecWeight(1.), passedCCnue(false), passedCCnumu(false), matched(false) {}

  bool operator==(const SelEvent& other) const {
    if (run == other.run && subrun == other.subrun && event == other.event && std::abs(trueNuE - other.trueNuE) < 0.01) 
      return true;
    if (run == other.run && subrun == other.subrun && event == other.event)
      std::cout << "true neutrino energy difference for matching events: " << std::abs(trueNuE - other.trueNuE) << std::endl;
    return false;
  }

  bool operator<(const SelEvent& other) const {
    if (run < other.run) return true;
    if (run > other.run) return false;
    if (subrun < other.subrun) return true;
    if (subrun > other.subrun) return false;
    if (event < other.event) return true;
    return false;
  }

};


struct SelResult {
  bool passed;
  float recoLepP, recoLepCosTheta;
  SelResult() : passed(false), recoLepP(-99.), recoLepCosTheta(-9.) {}
};


SelResult isCCnueInc(int foundVertex, int vtxIsFiducial, float vtxFracHitsOnCosmic,
 int nTracks, int nShowers, int trackIsSecondary[], int trackClassified[],
 int showerIsSecondary[], int showerClassified[], int trackPID[], int showerPID[],
 int showerProcess[], float trackMuScore[], float showerCharge[], float showerElScore[],
 float showerPhScore[], float showerPiScore[], float showerRecoE[], float showerCosTheta[]){

  SelResult result;

  if(foundVertex == 0 || vtxIsFiducial != 1) return result;
  if(vtxFracHitsOnCosmic >= (1. - 1e-6)) return result;

  int nMuons = 0;
  float maxMuScore = -99.;  
  for(int iT = 0; iT < nTracks; ++iT){
    if(trackIsSecondary[iT] == 1 || trackClassified[iT] != 1) continue;
    if(trackPID[iT] == 13) ++nMuons;
    if(trackMuScore[iT] > maxMuScore) maxMuScore = trackMuScore[iT];
  }
  if(nMuons > 0 || maxMuScore >= -3.7) return result;

  int nElectrons = 0;
  float elMaxQ = -99.;
  float elMaxQConf = -9.;
  float elMaxQEnergy = -99.;
  float elMaxQCosTheta = -9.;
  int elMaxQProc = -1;
  for(int iS = 0; iS < nShowers; ++iS){
    if(showerIsSecondary[iS] == 1 || showerClassified[iS] != 1) continue;
    if(showerPID[iS] == 11){
      ++nElectrons;
      if(showerCharge[iS] > elMaxQ){
        elMaxQ = showerCharge[iS];
        elMaxQProc = showerProcess[iS];
        elMaxQConf = showerElScore[iS] - (showerPhScore[iS] + showerPiScore[iS])/2.;
        elMaxQEnergy = showerRecoE[iS];
        elMaxQCosTheta = showerCosTheta[iS];
      }
    }
  }
  if(nElectrons < 1 || elMaxQProc != 0 || elMaxQConf <= 7.1) return result;

  result.passed = true;
  result.recoLepP = std::sqrt(std::pow(elMaxQEnergy,2.) - std::pow(0.511,2.));
  result.recoLepCosTheta = elMaxQCosTheta;

  return result;

}


SelResult isCCnumuInc(int foundVertex, int vtxIsFiducial, float vtxFracHitsOnCosmic,
 int nTracks, int trackIsSecondary[], int trackClassified[],
 int trackPID[], float trackMuScore[], float trackRecoE[], float trackCosTheta[]){

  SelResult result;

  if(foundVertex == 0 || vtxIsFiducial != 1) return result;
  if(vtxFracHitsOnCosmic >= (1. - 1e-6)) return result;

  int nMuons = 0;
  float maxMuScore = -99.;
  float muSelEnergy = -99.;
  float muSelCosTheta = -99.;
  for(int iT = 0; iT < nTracks; ++iT){
    if(trackIsSecondary[iT] == 1 || trackClassified[iT] != 1) continue;
    if(trackPID[iT] == 13){
      ++nMuons;
      if(trackMuScore[iT] > maxMuScore){
        maxMuScore = trackMuScore[iT];
        muSelEnergy = trackRecoE[iT];
        muSelCosTheta = trackCosTheta[iT];
      }
    }
  }
  if(nMuons < 1) return result;

  result.passed = true;
  result.recoLepP = std::sqrt(std::pow(muSelEnergy+105.66,2.) - std::pow(105.66,2.));
  result.recoLepCosTheta = muSelCosTheta;

  return result;

}


float overflow(float val, float binH){
  val /= 1000.;
  if(val >= binH) val = binH-0.01;
  return val;
}

float underflow(float val, float binL){
  //val /= 1000.;
  if(val <= binL) val = binL+0.01;
  return val;
}


void getNTuples(const std::string& var, std::vector<std::string>& ntuples){
  std::string ntupleDir = "/home/matthew/microboone/tufts/gen2val/flat_ntuples/detvar/";
  if(var == "SCE" || var == "recomb2"){
    ntuples.push_back(ntupleDir+"dlgen2_reco_v2me06_ntuple_v7_mcc9_v40a_dl_run3b_bnb_nu_overlay_500k_CV_weightsAdded.root");
    ntuples.push_back(ntupleDir+"dlgen2_reco_v2me06_ntuple_v7_mcc9_v40a_dl_run3b_bnb_nu_overlay_500k_"+var+"_weightsAdded.root");
  }
  else{
    ntuples.push_back(ntupleDir+"dlgen2_reco_v2me06_ntuple_v7_altCVnuWeights_mcc9_v29e_dl_run3b_bnb_nu_overlay_1mil_CV.root");
    ntuples.push_back(ntupleDir+"dlgen2_reco_v2me06_ntuple_v7_altCVnuWeights_mcc9_v40a_dl_run3b_bnb_nu_overlay_1mil_"+var+".root");
  }
  ntuples.push_back(ntupleDir+"dlgen2_reco_v2me06_ntuple_v7_mcc9_v29e_dl_run3b_intrinsic_nue_overlay_CV_weightsAdded.root");
  if(var == "LYAtt"){
    ntuples.push_back(ntupleDir+"dlgen2_reco_v2me06_ntuple_v7_mcc9_v40a_dl_run3b_intrinsic_nue_overlay_LYAtt-bugFix_weightsAdded.root");
  }
  else if (var.find("wiremod") != std::string::npos){
    ntuples.push_back(ntupleDir+"dlgen2_reco_v2me06_ntuple_v7_mcc9_v40a_dl_run3b_intrinsic_nue_overlay_"+var+"_weightsAdded.root");
  }
  else{
    ntuples.push_back(ntupleDir+"dlgen2_reco_v2me06_ntuple_v7_mcc9_v29e_dl_run3b_intrinsic_nue_overlay_"+var+"_weightsAdded.root");
  }
}


int main(int argc, char* argv[]) {

  std::string var = "";
  std::string outfile = "calculate_detvar_systematics_output.root";
  bool debug = false;
  bool ccnumuOverflowBinning = false;
  bool ccnueChi2Binning = false;
  bool ccnueBkgHack = false;
  int iarg = 0;

  while(iarg != -1){
    iarg = getopt(argc, argv, "v:o:dmcbh");
    switch(iarg){
      case 'v': var = optarg; break;
      case 'o': outfile = optarg; break;
      case 'd': debug = true; break;
      case 'm': ccnumuOverflowBinning = true; break;
      case 'c': ccnueChi2Binning = true; break;
      case 'b': ccnueBkgHack = true; break;
      case 'h': std::cout << "\t-v\tdetector varation to process" << std::endl;
                std::cout << "\t-o\toutput file name" << std::endl;
                std::cout << "\t-d\tflag, print extra debug info" << std::endl;
                std::cout << "\t-m\tflag, use new CCnumu overflow binning" << std::endl;
                return 0;
      case '?': std::cout << "unkown input option, run with -h to see options" << std::endl;
                return 0;
    }
  }

  if(var == ""){
    std::cout << "must provide detector variation. exiting..." << std::endl;
    return 0;
  }

  TH1::SetDefaultSumw2(kTRUE);

  std::cout << "setting up ntuple TTrees" << std::endl;

  //ntuple file names
  std::vector<std::string> ntuple_filenames;
  getNTuples(var, ntuple_filenames);
  std::string bnb_CV_file = ntuple_filenames[0];
  std::string bnb_var_file = ntuple_filenames[1];
  std::string nue_CV_file = ntuple_filenames[2];
  std::string nue_var_file = ntuple_filenames[3];

  std::cout << "using bnb CV file: " << bnb_CV_file << std::endl;
  std::cout << "using bnb var file: " << bnb_var_file << std::endl;
  std::cout << "using nue CV file: " << nue_CV_file << std::endl;
  std::cout << "using nue var file: " << nue_var_file << std::endl;


  int maxProngs = 100;
  //ntuple variables
  int run, subrun, event;
  int trueNuPDG, trueNuCCNC;
  int foundVertex, vtxIsFiducial;
  float trueNuE, recoNuE;
  float vtxFracHitsOnCosmic;
  float xsecWeight;
  int nTracks, nShowers;
  int trackIsSecondary[maxProngs];
  int trackClassified[maxProngs];
  int showerIsSecondary[maxProngs];
  int showerClassified[maxProngs];
  int trackPID[maxProngs];
  int showerPID[maxProngs];
  int showerProcess[maxProngs];
  float trackMuScore[maxProngs];
  float showerCharge[maxProngs];
  float showerElScore[maxProngs];
  float showerPhScore[maxProngs];
  float showerPiScore[maxProngs];
  float trackRecoE[maxProngs];
  float showerRecoE[maxProngs];
  float trackCosTheta[maxProngs];
  float showerCosTheta[maxProngs];

  //ntuple file/tree declarations
  TFile* ntuple_nue_CV = new TFile(nue_CV_file.c_str(), "READ");
  TTree* tree_nue_CV = (TTree*)ntuple_nue_CV->Get("EventTree");
  tree_nue_CV->SetBranchAddress("run",&run);
  tree_nue_CV->SetBranchAddress("subrun",&subrun);
  tree_nue_CV->SetBranchAddress("event",&event);
  tree_nue_CV->SetBranchAddress("xsecWeight",&xsecWeight);
  tree_nue_CV->SetBranchAddress("trueNuPDG",&trueNuPDG);
  tree_nue_CV->SetBranchAddress("trueNuCCNC",&trueNuCCNC);
  tree_nue_CV->SetBranchAddress("foundVertex",&foundVertex);
  tree_nue_CV->SetBranchAddress("vtxIsFiducial",&vtxIsFiducial);
  tree_nue_CV->SetBranchAddress("trueNuE",&trueNuE);
  tree_nue_CV->SetBranchAddress("recoNuE",&recoNuE);
  tree_nue_CV->SetBranchAddress("vtxFracHitsOnCosmic",&vtxFracHitsOnCosmic);
  tree_nue_CV->SetBranchAddress("nTracks",&nTracks);
  tree_nue_CV->SetBranchAddress("nShowers",&nShowers);
  tree_nue_CV->SetBranchAddress("trackIsSecondary",trackIsSecondary);
  tree_nue_CV->SetBranchAddress("trackClassified",trackClassified);
  tree_nue_CV->SetBranchAddress("showerIsSecondary",showerIsSecondary);
  tree_nue_CV->SetBranchAddress("showerClassified",showerClassified);
  tree_nue_CV->SetBranchAddress("trackPID",trackPID);
  tree_nue_CV->SetBranchAddress("showerPID",showerPID);
  tree_nue_CV->SetBranchAddress("showerProcess",showerProcess);
  tree_nue_CV->SetBranchAddress("trackMuScore",trackMuScore);
  tree_nue_CV->SetBranchAddress("showerCharge",showerCharge);
  tree_nue_CV->SetBranchAddress("showerElScore",showerElScore);
  tree_nue_CV->SetBranchAddress("showerPhScore",showerPhScore);
  tree_nue_CV->SetBranchAddress("showerPiScore",showerPiScore);
  tree_nue_CV->SetBranchAddress("trackRecoE",trackRecoE);
  tree_nue_CV->SetBranchAddress("showerRecoE",showerRecoE);
  tree_nue_CV->SetBranchAddress("trackCosTheta",trackCosTheta);
  tree_nue_CV->SetBranchAddress("showerCosTheta",showerCosTheta);

  TFile* ntuple_nue_var = new TFile(nue_var_file.c_str(), "READ");
  TTree* tree_nue_var = (TTree*)ntuple_nue_var->Get("EventTree");
  tree_nue_var->SetBranchAddress("run",&run);
  tree_nue_var->SetBranchAddress("subrun",&subrun);
  tree_nue_var->SetBranchAddress("event",&event);
  tree_nue_var->SetBranchAddress("xsecWeight",&xsecWeight);
  tree_nue_var->SetBranchAddress("trueNuPDG",&trueNuPDG);
  tree_nue_var->SetBranchAddress("trueNuCCNC",&trueNuCCNC);
  tree_nue_var->SetBranchAddress("foundVertex",&foundVertex);
  tree_nue_var->SetBranchAddress("vtxIsFiducial",&vtxIsFiducial);
  tree_nue_var->SetBranchAddress("trueNuE",&trueNuE);
  tree_nue_var->SetBranchAddress("recoNuE",&recoNuE);
  tree_nue_var->SetBranchAddress("vtxFracHitsOnCosmic",&vtxFracHitsOnCosmic);
  tree_nue_var->SetBranchAddress("nTracks",&nTracks);
  tree_nue_var->SetBranchAddress("nShowers",&nShowers);
  tree_nue_var->SetBranchAddress("trackIsSecondary",trackIsSecondary);
  tree_nue_var->SetBranchAddress("trackClassified",trackClassified);
  tree_nue_var->SetBranchAddress("showerIsSecondary",showerIsSecondary);
  tree_nue_var->SetBranchAddress("showerClassified",showerClassified);
  tree_nue_var->SetBranchAddress("trackPID",trackPID);
  tree_nue_var->SetBranchAddress("showerPID",showerPID);
  tree_nue_var->SetBranchAddress("showerProcess",showerProcess);
  tree_nue_var->SetBranchAddress("trackMuScore",trackMuScore);
  tree_nue_var->SetBranchAddress("showerCharge",showerCharge);
  tree_nue_var->SetBranchAddress("showerElScore",showerElScore);
  tree_nue_var->SetBranchAddress("showerPhScore",showerPhScore);
  tree_nue_var->SetBranchAddress("showerPiScore",showerPiScore);
  tree_nue_var->SetBranchAddress("trackRecoE",trackRecoE);
  tree_nue_var->SetBranchAddress("showerRecoE",showerRecoE);
  tree_nue_var->SetBranchAddress("trackCosTheta",trackCosTheta);
  tree_nue_var->SetBranchAddress("showerCosTheta",showerCosTheta);

  TFile* ntuple_bnb_CV = new TFile(bnb_CV_file.c_str(), "READ");
  TTree* tree_bnb_CV = (TTree*)ntuple_bnb_CV->Get("EventTree");
  tree_bnb_CV->SetBranchAddress("run",&run);
  tree_bnb_CV->SetBranchAddress("subrun",&subrun);
  tree_bnb_CV->SetBranchAddress("event",&event);
  tree_bnb_CV->SetBranchAddress("xsecWeight",&xsecWeight);
  tree_bnb_CV->SetBranchAddress("trueNuPDG",&trueNuPDG);
  tree_bnb_CV->SetBranchAddress("trueNuCCNC",&trueNuCCNC);
  tree_bnb_CV->SetBranchAddress("foundVertex",&foundVertex);
  tree_bnb_CV->SetBranchAddress("vtxIsFiducial",&vtxIsFiducial);
  tree_bnb_CV->SetBranchAddress("trueNuE",&trueNuE);
  tree_bnb_CV->SetBranchAddress("recoNuE",&recoNuE);
  tree_bnb_CV->SetBranchAddress("vtxFracHitsOnCosmic",&vtxFracHitsOnCosmic);
  tree_bnb_CV->SetBranchAddress("nTracks",&nTracks);
  tree_bnb_CV->SetBranchAddress("nShowers",&nShowers);
  tree_bnb_CV->SetBranchAddress("trackIsSecondary",trackIsSecondary);
  tree_bnb_CV->SetBranchAddress("trackClassified",trackClassified);
  tree_bnb_CV->SetBranchAddress("showerIsSecondary",showerIsSecondary);
  tree_bnb_CV->SetBranchAddress("showerClassified",showerClassified);
  tree_bnb_CV->SetBranchAddress("trackPID",trackPID);
  tree_bnb_CV->SetBranchAddress("showerPID",showerPID);
  tree_bnb_CV->SetBranchAddress("showerProcess",showerProcess);
  tree_bnb_CV->SetBranchAddress("trackMuScore",trackMuScore);
  tree_bnb_CV->SetBranchAddress("showerCharge",showerCharge);
  tree_bnb_CV->SetBranchAddress("showerElScore",showerElScore);
  tree_bnb_CV->SetBranchAddress("showerPhScore",showerPhScore);
  tree_bnb_CV->SetBranchAddress("showerPiScore",showerPiScore);
  tree_bnb_CV->SetBranchAddress("trackRecoE",trackRecoE);
  tree_bnb_CV->SetBranchAddress("showerRecoE",showerRecoE);
  tree_bnb_CV->SetBranchAddress("trackCosTheta",trackCosTheta);
  tree_bnb_CV->SetBranchAddress("showerCosTheta",showerCosTheta);

  TFile* ntuple_bnb_var = new TFile(bnb_var_file.c_str(), "READ");
  TTree* tree_bnb_var = (TTree*)ntuple_bnb_var->Get("EventTree");
  tree_bnb_var->SetBranchAddress("run",&run);
  tree_bnb_var->SetBranchAddress("subrun",&subrun);
  tree_bnb_var->SetBranchAddress("event",&event);
  tree_bnb_var->SetBranchAddress("xsecWeight",&xsecWeight);
  tree_bnb_var->SetBranchAddress("trueNuPDG",&trueNuPDG);
  tree_bnb_var->SetBranchAddress("trueNuCCNC",&trueNuCCNC);
  tree_bnb_var->SetBranchAddress("foundVertex",&foundVertex);
  tree_bnb_var->SetBranchAddress("vtxIsFiducial",&vtxIsFiducial);
  tree_bnb_var->SetBranchAddress("trueNuE",&trueNuE);
  tree_bnb_var->SetBranchAddress("recoNuE",&recoNuE);
  tree_bnb_var->SetBranchAddress("vtxFracHitsOnCosmic",&vtxFracHitsOnCosmic);
  tree_bnb_var->SetBranchAddress("nTracks",&nTracks);
  tree_bnb_var->SetBranchAddress("nShowers",&nShowers);
  tree_bnb_var->SetBranchAddress("trackIsSecondary",trackIsSecondary);
  tree_bnb_var->SetBranchAddress("trackClassified",trackClassified);
  tree_bnb_var->SetBranchAddress("showerIsSecondary",showerIsSecondary);
  tree_bnb_var->SetBranchAddress("showerClassified",showerClassified);
  tree_bnb_var->SetBranchAddress("trackPID",trackPID);
  tree_bnb_var->SetBranchAddress("showerPID",showerPID);
  tree_bnb_var->SetBranchAddress("showerProcess",showerProcess);
  tree_bnb_var->SetBranchAddress("trackMuScore",trackMuScore);
  tree_bnb_var->SetBranchAddress("showerCharge",showerCharge);
  tree_bnb_var->SetBranchAddress("showerElScore",showerElScore);
  tree_bnb_var->SetBranchAddress("showerPhScore",showerPhScore);
  tree_bnb_var->SetBranchAddress("showerPiScore",showerPiScore);
  tree_bnb_var->SetBranchAddress("trackRecoE",trackRecoE);
  tree_bnb_var->SetBranchAddress("showerRecoE",showerRecoE);
  tree_bnb_var->SetBranchAddress("trackCosTheta",trackCosTheta);
  tree_bnb_var->SetBranchAddress("showerCosTheta",showerCosTheta);


  std::cout << "analyzing events from intrinsic nue CV ntuple" << std::endl;

  std::vector<SelEvent> nue_CV_events;
  nue_CV_events.reserve(tree_nue_CV->GetEntries());

  for(int i = 0; i < tree_nue_CV->GetEntries(); ++i){
    tree_nue_CV -> GetEntry(i);
    if(std::abs(trueNuPDG) != 12 || trueNuCCNC != 0) continue;
    SelResult CCnueResult = isCCnueInc(foundVertex, vtxIsFiducial, vtxFracHitsOnCosmic,
     nTracks, nShowers, trackIsSecondary, trackClassified, showerIsSecondary, showerClassified,
     trackPID, showerPID, showerProcess, trackMuScore, showerCharge, showerElScore,
     showerPhScore, showerPiScore, showerRecoE, showerCosTheta);
    SelResult CCnumuResult = isCCnumuInc(foundVertex, vtxIsFiducial, vtxFracHitsOnCosmic, nTracks,
     trackIsSecondary, trackClassified, trackPID, trackMuScore, trackRecoE, trackCosTheta);
    SelEvent evt(run, subrun, event, trueNuE*1000., recoNuE, xsecWeight, CCnueResult.passed, CCnumuResult.passed, false);
    if(CCnueResult.passed){
      evt.electron.momentum = CCnueResult.recoLepP;
      evt.electron.cosTheta = CCnueResult.recoLepCosTheta;
    }
    if(CCnumuResult.passed){
      evt.muon.momentum = CCnumuResult.recoLepP;
      evt.muon.cosTheta = CCnumuResult.recoLepCosTheta;
    }
    nue_CV_events.push_back(evt);
  }

  std::sort(nue_CV_events.begin(), nue_CV_events.end());

  std::cout << "looked at "<<nue_CV_events.size()<<" intrinsic nue CV events" << std::endl;

  std::cout << "analyzing events from bnb overlay CV ntuple" << std::endl;

  std::vector<SelEvent> bnb_CV_events;
  bnb_CV_events.reserve(tree_bnb_CV->GetEntries());

  for(int i = 0; i < tree_bnb_CV->GetEntries(); ++i){
    tree_bnb_CV -> GetEntry(i);
    if(std::abs(trueNuPDG) != 14 && trueNuCCNC != 1) continue;
    SelResult CCnueResult = isCCnueInc(foundVertex, vtxIsFiducial, vtxFracHitsOnCosmic,
     nTracks, nShowers, trackIsSecondary, trackClassified, showerIsSecondary, showerClassified,
     trackPID, showerPID, showerProcess, trackMuScore, showerCharge, showerElScore,
     showerPhScore, showerPiScore, showerRecoE, showerCosTheta);
    SelResult CCnumuResult = isCCnumuInc(foundVertex, vtxIsFiducial, vtxFracHitsOnCosmic, nTracks,
     trackIsSecondary, trackClassified, trackPID, trackMuScore, trackRecoE, trackCosTheta);
    SelEvent evt(run, subrun, event, trueNuE*1000., recoNuE, xsecWeight, CCnueResult.passed, CCnumuResult.passed, false);
    if(CCnueResult.passed){
      evt.electron.momentum = CCnueResult.recoLepP;
      evt.electron.cosTheta = CCnueResult.recoLepCosTheta;
    }
    if(CCnumuResult.passed){
      evt.muon.momentum = CCnumuResult.recoLepP;
      evt.muon.cosTheta = CCnumuResult.recoLepCosTheta;
    }
    bnb_CV_events.push_back(evt);
  }

  std::sort(bnb_CV_events.begin(), bnb_CV_events.end());

  std::cout << "looked at "<<bnb_CV_events.size()<<" bnb overlay CV events" << std::endl;

  std::cout << "analyzing matched events from intrinsic nue detvar ntuple" << std::endl;

  std::vector<SelEvent> nue_var_selected_events;
  nue_var_selected_events.reserve(nue_CV_events.size());
  int n_nue_matched = 0;

  for(int i = 0; i < tree_nue_var->GetEntries(); ++i){
    tree_nue_var -> GetEntry(i);
    if(std::abs(trueNuPDG) != 12 || trueNuCCNC != 0) continue;
    SelEvent evt(run, subrun, event, trueNuE*1000., recoNuE, xsecWeight, false, false, false);
    auto it_CV = std::lower_bound(nue_CV_events.begin(), nue_CV_events.end(), evt);
    if(it_CV == nue_CV_events.end() || !(*it_CV == evt)) continue;
    (*it_CV).matched = true; evt.matched = true; ++n_nue_matched;
    SelResult CCnueResult = isCCnueInc(foundVertex, vtxIsFiducial, vtxFracHitsOnCosmic,
     nTracks, nShowers, trackIsSecondary, trackClassified, showerIsSecondary, showerClassified,
     trackPID, showerPID, showerProcess, trackMuScore, showerCharge, showerElScore,
     showerPhScore, showerPiScore, showerRecoE, showerCosTheta);
    SelResult CCnumuResult = isCCnumuInc(foundVertex, vtxIsFiducial, vtxFracHitsOnCosmic, nTracks,
     trackIsSecondary, trackClassified, trackPID, trackMuScore, trackRecoE, trackCosTheta);
    if(CCnueResult.passed || CCnumuResult.passed){
      if(CCnueResult.passed){
        evt.passedCCnue = true;
        evt.electron.momentum = CCnueResult.recoLepP;
        evt.electron.cosTheta = CCnueResult.recoLepCosTheta;
      }
      if(CCnumuResult.passed){
        evt.passedCCnumu = true;
        evt.muon.momentum = CCnumuResult.recoLepP;
        evt.muon.cosTheta = CCnumuResult.recoLepCosTheta;
      }
      nue_var_selected_events.push_back(evt);
    }
  }

  std::cout << "found "<<nue_var_selected_events.size()<<" intrinsic nue detvar events from CV intersection (total of "<<n_nue_matched<<" events) that pass a selection" << std::endl;

  std::cout << "analyzing matched events from bnb overlay detvar ntuple" << std::endl;

  std::vector<SelEvent> bnb_var_selected_events;
  bnb_var_selected_events.reserve(bnb_CV_events.size());
  int n_bnb_matched = 0;

  for(int i = 0; i < tree_bnb_var->GetEntries(); ++i){
    tree_bnb_var -> GetEntry(i);
    if(std::abs(trueNuPDG) != 14 && trueNuCCNC != 1) continue;
    SelEvent evt(run, subrun, event, trueNuE*1000., recoNuE, xsecWeight, false, false, false);
    auto it_CV = std::lower_bound(bnb_CV_events.begin(), bnb_CV_events.end(), evt);
    if(it_CV == bnb_CV_events.end() || !(*it_CV == evt)) continue;
    (*it_CV).matched = true; evt.matched = true; ++n_bnb_matched;
    SelResult CCnueResult = isCCnueInc(foundVertex, vtxIsFiducial, vtxFracHitsOnCosmic,
     nTracks, nShowers, trackIsSecondary, trackClassified, showerIsSecondary, showerClassified,
     trackPID, showerPID, showerProcess, trackMuScore, showerCharge, showerElScore,
     showerPhScore, showerPiScore, showerRecoE, showerCosTheta);
    SelResult CCnumuResult = isCCnumuInc(foundVertex, vtxIsFiducial, vtxFracHitsOnCosmic, nTracks,
     trackIsSecondary, trackClassified, trackPID, trackMuScore, trackRecoE, trackCosTheta);
    if(CCnueResult.passed || CCnumuResult.passed){
      if(CCnueResult.passed){
        evt.passedCCnue = true;
        evt.electron.momentum = CCnueResult.recoLepP;
        evt.electron.cosTheta = CCnueResult.recoLepCosTheta;
      }
      if(CCnumuResult.passed){
        evt.passedCCnumu = true;
        evt.muon.momentum = CCnumuResult.recoLepP;
        evt.muon.cosTheta = CCnumuResult.recoLepCosTheta;
      }
      bnb_var_selected_events.push_back(evt);
    }
  }

  std::cout << "found "<<bnb_var_selected_events.size()<<" bnb overlay detvar events from CV intersection (total of "<<n_bnb_matched<<" events) that pass a selection" << std::endl;

  std::cout << "calculating POT scaling factors" << std::endl;

  float totGoodPOT;
  float nue_POT = 0.;
  float bnb_POT = 0.;
  float targetPOT = 4.4e+19;

  TTree* nue_potTree = (TTree*)ntuple_nue_CV->Get("potTree");
  nue_potTree -> SetBranchAddress("totGoodPOT",&totGoodPOT);

  for(int i = 0; i < nue_potTree->GetEntries(); ++i){
    nue_potTree -> GetEntry(i);
    nue_POT += totGoodPOT;
  }
  nue_POT *= ((1.0*n_nue_matched)/(1.0*nue_CV_events.size()));

  TTree* bnb_potTree = (TTree*)ntuple_bnb_CV->Get("potTree");
  bnb_potTree -> SetBranchAddress("totGoodPOT",&totGoodPOT);

  for(int i = 0; i < bnb_potTree->GetEntries(); ++i){
    bnb_potTree -> GetEntry(i);
    bnb_POT += totGoodPOT;
  }
  bnb_POT *= ((1.0*n_bnb_matched)/(1.0*bnb_CV_events.size()));


  std::cout << "filling histograms" << std::endl;

  //original binning
  int CCnue_recoNuE_nBins = 14;
  float CCnue_recoNuE_binL = 0.;
  float CCnue_recoNuE_binH = 2.8;
  //chi2 binning
  if(ccnueChi2Binning){
    CCnue_recoNuE_nBins = 9;
    CCnue_recoNuE_binL = 0.2;
    CCnue_recoNuE_binH = 2.0;
  }

  //original binning
  int CCnumu_recoNuE_nBins = 21;
  float CCnumu_recoNuE_binL = 0.;
  float CCnumu_recoNuE_binH = 2.1;
  //new CCnumu overflow binning
  if(ccnumuOverflowBinning){
    CCnumu_recoNuE_nBins = 19;
    CCnumu_recoNuE_binL = 0.;
    CCnumu_recoNuE_binH = 1.9;
  }

  //original binning
  int CCnue_recoLepP_nBins = 14;
  float CCnue_recoLepP_binL = 0.;
  float CCnue_recoLepP_binH = 2.8;
  //chi2 binning
  if(ccnueChi2Binning){
    CCnue_recoLepP_nBins = 8;
    CCnue_recoLepP_binL = 0.;
    CCnue_recoLepP_binH = 1.6;
  }

  //original binning
  int CCnumu_recoLepP_nBins = 21;
  float CCnumu_recoLepP_binL = 0.;
  float CCnumu_recoLepP_binH = 2.1;
  //new CCnumu overflow binning
  if(ccnumuOverflowBinning){
    CCnumu_recoLepP_nBins = 16;
    CCnumu_recoLepP_binL = 0.;
    CCnumu_recoLepP_binH = 1.6;
  }

  //original binning
  int CCnue_recoLepCosTheta_nBins = 16;
  float CCnue_recoLepCosTheta_binL = -1.;
  float CCnue_recoLepCosTheta_binH = 1.;
  //chi2 binning
  if(ccnueChi2Binning){
    CCnue_recoLepCosTheta_nBins = 6;
    CCnue_recoLepCosTheta_binL = 0.25;
    CCnue_recoLepCosTheta_binH = 1.;
  }

  //chi2 and original binning
  int CCnumu_recoLepCosTheta_nBins = 16;
  float CCnumu_recoLepCosTheta_binL = -1.;
  float CCnumu_recoLepCosTheta_binH = 1.;


  TH1F h_CCnue_recoNuE_CV("h_CCnue_recoNuE_CV","CCnue Selection, Reco Neutrino Energy (MeV)", CCnue_recoNuE_nBins, CCnue_recoNuE_binL, CCnue_recoNuE_binH);
  TH1F h_CCnue_recoLepP_CV("h_CCnue_recoLepP_CV","CCnue Selection, Reco Electron Momentum (MeV/c)", CCnue_recoLepP_nBins, CCnue_recoLepP_binL, CCnue_recoLepP_binH);
  TH1F h_CCnue_recoLepCosTheta_CV("h_CCnue_recoLepCosTheta_CV","CCnue Selection, Reco Electron cos(theta)", CCnue_recoLepCosTheta_nBins, CCnue_recoLepCosTheta_binL, CCnue_recoLepCosTheta_binH);

  TH1F h_CCnue_recoNuE_var("h_CCnue_recoNuE_var","CCnue Selection, Reco Neutrino Energy (MeV)", CCnue_recoNuE_nBins, CCnue_recoNuE_binL, CCnue_recoNuE_binH);
  TH1F h_CCnue_recoLepP_var("h_CCnue_recoLepP_var","CCnue Selection, Reco Electron Momentum (MeV/c)", CCnue_recoLepP_nBins, CCnue_recoLepP_binL, CCnue_recoLepP_binH);
  TH1F h_CCnue_recoLepCosTheta_var("h_CCnue_recoLepCosTheta_var","CCnue Selection, Reco Electron cos(theta)", CCnue_recoLepCosTheta_nBins, CCnue_recoLepCosTheta_binL, CCnue_recoLepCosTheta_binH);

  TH1F h_CCnumu_recoNuE_CV("h_CCnumu_recoNuE_CV","CCnumu Selection, Reco Neutrino Energy (MeV)", CCnumu_recoNuE_nBins, CCnumu_recoNuE_binL, CCnumu_recoNuE_binH);
  TH1F h_CCnumu_recoNuE_CV_error("h_CCnumu_recoNuE_CV_error","CCnumu Selection, Reco Neutrino Energy (MeV)", CCnumu_recoNuE_nBins, CCnumu_recoNuE_binL, CCnumu_recoNuE_binH);
  TH1F h_CCnumu_recoLepP_CV("h_CCnumu_recoLepP_CV","CCnumu Selection, Reco Muon Momentum (MeV/c)", CCnumu_recoLepP_nBins, CCnumu_recoLepP_binL, CCnumu_recoLepP_binH);
  TH1F h_CCnumu_recoLepCosTheta_CV("h_CCnumu_recoLepCosTheta_CV","CCnumu Selection, Reco Muon cos(theta)", CCnumu_recoLepCosTheta_nBins, CCnumu_recoLepCosTheta_binL, CCnumu_recoLepCosTheta_binH);

  TH1F h_CCnumu_recoNuE_var("h_CCnumu_recoNuE_var","CCnumu Selection, Reco Neutrino Energy (MeV)", CCnumu_recoNuE_nBins, CCnumu_recoNuE_binL, CCnumu_recoNuE_binH);
  TH1F h_CCnumu_recoLepP_var("h_CCnumu_recoLepP_var","CCnumu Selection, Reco Muon Momentum (MeV/c)", CCnumu_recoLepP_nBins, CCnumu_recoLepP_binL, CCnumu_recoLepP_binH);
  TH1F h_CCnumu_recoLepCosTheta_var("h_CCnumu_recoLepCosTheta_var","CCnumu Selection, Reco Muon cos(theta)", CCnumu_recoLepCosTheta_nBins, CCnumu_recoLepCosTheta_binL, CCnumu_recoLepCosTheta_binH);

  TH1F hRawInt_CCnue_recoNuE_CV("hRawInt_CCnue_recoNuE_CV","CCnue Selection, Reco Neutrino Energy (MeV)", CCnue_recoNuE_nBins, CCnue_recoNuE_binL, CCnue_recoNuE_binH);
  TH1F hRawInt_CCnue_recoLepP_CV("hRawInt_CCnue_recoLepP_CV","CCnue Selection, Reco Electron Momentum (MeV/c)", CCnue_recoLepP_nBins, CCnue_recoLepP_binL, CCnue_recoLepP_binH);
  TH1F hRawInt_CCnue_recoLepCosTheta_CV("hRawInt_CCnue_recoLepCosTheta_CV","CCnue Selection, Reco Electron cos(theta)", CCnue_recoLepCosTheta_nBins, CCnue_recoLepCosTheta_binL, CCnue_recoLepCosTheta_binH);

  TH1F hRawInt_CCnue_recoNuE_var("hRawInt_CCnue_recoNuE_var","CCnue Selection, Reco Neutrino Energy (MeV)", CCnue_recoNuE_nBins, CCnue_recoNuE_binL, CCnue_recoNuE_binH);
  TH1F hRawInt_CCnue_recoLepP_var("hRawInt_CCnue_recoLepP_var","CCnue Selection, Reco Electron Momentum (MeV/c)", CCnue_recoLepP_nBins, CCnue_recoLepP_binL, CCnue_recoLepP_binH);
  TH1F hRawInt_CCnue_recoLepCosTheta_var("hRawInt_CCnue_recoLepCosTheta_var","CCnue Selection, Reco Electron cos(theta)", CCnue_recoLepCosTheta_nBins, CCnue_recoLepCosTheta_binL, CCnue_recoLepCosTheta_binH);

  TH1F hRawInt_CCnumu_recoNuE_CV("hRawInt_CCnumu_recoNuE_CV","CCnumu Selection, Reco Neutrino Energy (MeV)", CCnumu_recoNuE_nBins, CCnumu_recoNuE_binL, CCnumu_recoNuE_binH);
  TH1F hRawInt_CCnumu_recoLepP_CV("hRawInt_CCnumu_recoLepP_CV","CCnumu Selection, Reco Muon Momentum (MeV/c)", CCnumu_recoLepP_nBins, CCnumu_recoLepP_binL, CCnumu_recoLepP_binH);
  TH1F hRawInt_CCnumu_recoLepCosTheta_CV("hRawInt_CCnumu_recoLepCosTheta_CV","CCnumu Selection, Reco Muon cos(theta)", CCnumu_recoLepCosTheta_nBins, CCnumu_recoLepCosTheta_binL, CCnumu_recoLepCosTheta_binH);

  TH1F hRawInt_CCnumu_recoNuE_var("hRawInt_CCnumu_recoNuE_var","CCnumu Selection, Reco Neutrino Energy (MeV)", CCnumu_recoNuE_nBins, CCnumu_recoNuE_binL, CCnumu_recoNuE_binH);
  TH1F hRawInt_CCnumu_recoLepP_var("hRawInt_CCnumu_recoLepP_var","CCnumu Selection, Reco Muon Momentum (MeV/c)", CCnumu_recoLepP_nBins, CCnumu_recoLepP_binL, CCnumu_recoLepP_binH);
  TH1F hRawInt_CCnumu_recoLepCosTheta_var("hRawInt_CCnumu_recoLepCosTheta_var","CCnumu Selection, Reco Muon cos(theta)", CCnumu_recoLepCosTheta_nBins, CCnumu_recoLepCosTheta_binL, CCnumu_recoLepCosTheta_binH);

  TH1F hRawBNB_CCnue_recoNuE_CV("hRawBNB_CCnue_recoNuE_CV","CCnue Selection, Reco Neutrino Energy (MeV)", CCnue_recoNuE_nBins, CCnue_recoNuE_binL, CCnue_recoNuE_binH);
  TH1F hRawBNB_CCnue_recoLepP_CV("hRawBNB_CCnue_recoLepP_CV","CCnue Selection, Reco Electron Momentum (MeV/c)", CCnue_recoLepP_nBins, CCnue_recoLepP_binL, CCnue_recoLepP_binH);
  TH1F hRawBNB_CCnue_recoLepCosTheta_CV("hRawBNB_CCnue_recoLepCosTheta_CV","CCnue Selection, Reco Electron cos(theta)", CCnue_recoLepCosTheta_nBins, CCnue_recoLepCosTheta_binL, CCnue_recoLepCosTheta_binH);

  TH1F hRawBNB_CCnue_recoNuE_var("hRawBNB_CCnue_recoNuE_var","CCnue Selection, Reco Neutrino Energy (MeV)", CCnue_recoNuE_nBins, CCnue_recoNuE_binL, CCnue_recoNuE_binH);
  TH1F hRawBNB_CCnue_recoLepP_var("hRawBNB_CCnue_recoLepP_var","CCnue Selection, Reco Electron Momentum (MeV/c)", CCnue_recoLepP_nBins, CCnue_recoLepP_binL, CCnue_recoLepP_binH);
  TH1F hRawBNB_CCnue_recoLepCosTheta_var("hRawBNB_CCnue_recoLepCosTheta_var","CCnue Selection, Reco Electron cos(theta)", CCnue_recoLepCosTheta_nBins, CCnue_recoLepCosTheta_binL, CCnue_recoLepCosTheta_binH);

  TH1F hRawBNB_CCnumu_recoNuE_CV("hRawBNB_CCnumu_recoNuE_CV","CCnumu Selection, Reco Neutrino Energy (MeV)", CCnumu_recoNuE_nBins, CCnumu_recoNuE_binL, CCnumu_recoNuE_binH);
  TH1F hRawBNB_CCnumu_recoLepP_CV("hRawBNB_CCnumu_recoLepP_CV","CCnumu Selection, Reco Muon Momentum (MeV/c)", CCnumu_recoLepP_nBins, CCnumu_recoLepP_binL, CCnumu_recoLepP_binH);
  TH1F hRawBNB_CCnumu_recoLepCosTheta_CV("hRawBNB_CCnumu_recoLepCosTheta_CV","CCnumu Selection, Reco Muon cos(theta)", CCnumu_recoLepCosTheta_nBins, CCnumu_recoLepCosTheta_binL, CCnumu_recoLepCosTheta_binH);

  TH1F hRawBNB_CCnumu_recoNuE_var("hRawBNB_CCnumu_recoNuE_var","CCnumu Selection, Reco Neutrino Energy (MeV)", CCnumu_recoNuE_nBins, CCnumu_recoNuE_binL, CCnumu_recoNuE_binH);
  TH1F hRawBNB_CCnumu_recoLepP_var("hRawBNB_CCnumu_recoLepP_var","CCnumu Selection, Reco Muon Momentum (MeV/c)", CCnumu_recoLepP_nBins, CCnumu_recoLepP_binL, CCnumu_recoLepP_binH);
  TH1F hRawBNB_CCnumu_recoLepCosTheta_var("hRawBNB_CCnumu_recoLepCosTheta_var","CCnumu Selection, Reco Muon cos(theta)", CCnumu_recoLepCosTheta_nBins, CCnumu_recoLepCosTheta_binL, CCnumu_recoLepCosTheta_binH);

  for(const auto& evt : nue_CV_events){
    if(!evt.matched) continue;
    if(evt.passedCCnue){
      h_CCnue_recoNuE_CV.Fill(underflow(overflow(evt.recoNuE,CCnue_recoNuE_binH),CCnue_recoNuE_binL), evt.xsecWeight*(targetPOT/nue_POT));
      hRawInt_CCnue_recoNuE_CV.Fill(underflow(overflow(evt.recoNuE,CCnue_recoNuE_binH),CCnue_recoNuE_binL));
      h_CCnue_recoLepP_CV.Fill(overflow(evt.electron.momentum,CCnue_recoLepP_binH), evt.xsecWeight*(targetPOT/nue_POT));
      hRawInt_CCnue_recoLepP_CV.Fill(overflow(evt.electron.momentum,CCnue_recoLepP_binH));
      h_CCnue_recoLepCosTheta_CV.Fill(underflow(evt.electron.cosTheta,CCnue_recoLepCosTheta_binL), evt.xsecWeight*(targetPOT/nue_POT));
      hRawInt_CCnue_recoLepCosTheta_CV.Fill(underflow(evt.electron.cosTheta,CCnue_recoLepCosTheta_binL));
    }
    if(evt.passedCCnumu){
      h_CCnumu_recoNuE_CV.Fill(overflow(evt.recoNuE,CCnumu_recoNuE_binH), evt.xsecWeight*(targetPOT/nue_POT));
      h_CCnumu_recoNuE_CV_error.Fill(overflow(evt.recoNuE,CCnumu_recoNuE_binH), std::pow((targetPOT/nue_POT)*evt.xsecWeight,2));
      hRawInt_CCnumu_recoNuE_CV.Fill(overflow(evt.recoNuE,CCnumu_recoNuE_binH));
      h_CCnumu_recoLepP_CV.Fill(overflow(evt.muon.momentum,CCnumu_recoLepP_binH), evt.xsecWeight*(targetPOT/nue_POT));
      hRawInt_CCnumu_recoLepP_CV.Fill(overflow(evt.muon.momentum,CCnumu_recoLepP_binH));
      h_CCnumu_recoLepCosTheta_CV.Fill(evt.muon.cosTheta, evt.xsecWeight*(targetPOT/nue_POT));
      hRawInt_CCnumu_recoLepCosTheta_CV.Fill(underflow(evt.muon.cosTheta,CCnumu_recoLepCosTheta_binL));
    }
  }

  float n_bnb_CCnue_CV = 0.0;
  for(const auto& evt : bnb_CV_events){
    if(!evt.matched) continue;
    if(evt.passedCCnue){
      if(debug && evt.recoNuE < 200){
        std::cout << "found bnb CV event with weight "<<evt.xsecWeight<<" passing CCnue selection, run/subrun/event = "<<evt.run<<"/"<<evt.subrun<<"/"<<evt.event << std::endl;
      }
      n_bnb_CCnue_CV += evt.xsecWeight*(targetPOT/bnb_POT);
      if(!ccnueBkgHack){
        h_CCnue_recoNuE_CV.Fill(underflow(overflow(evt.recoNuE,CCnue_recoNuE_binH),CCnue_recoNuE_binL), evt.xsecWeight*(targetPOT/bnb_POT));
        h_CCnue_recoLepP_CV.Fill(overflow(evt.electron.momentum,CCnue_recoLepP_binH), evt.xsecWeight*(targetPOT/bnb_POT));
        h_CCnue_recoLepCosTheta_CV.Fill(underflow(evt.electron.cosTheta,CCnue_recoLepCosTheta_binL), evt.xsecWeight*(targetPOT/bnb_POT));
      }
      hRawBNB_CCnue_recoNuE_CV.Fill(underflow(overflow(evt.recoNuE,CCnue_recoNuE_binH),CCnue_recoNuE_binL));
      hRawBNB_CCnue_recoLepP_CV.Fill(overflow(evt.electron.momentum,CCnue_recoLepP_binH));
      hRawBNB_CCnue_recoLepCosTheta_CV.Fill(underflow(evt.electron.cosTheta,CCnue_recoLepCosTheta_binL));
    }
    if(evt.passedCCnumu){
      h_CCnumu_recoNuE_CV.Fill(overflow(evt.recoNuE,CCnumu_recoNuE_binH), evt.xsecWeight*(targetPOT/bnb_POT));
      h_CCnumu_recoNuE_CV_error.Fill(overflow(evt.recoNuE,CCnumu_recoNuE_binH), std::pow((targetPOT/bnb_POT)*evt.xsecWeight,2));
      hRawBNB_CCnumu_recoNuE_CV.Fill(overflow(evt.recoNuE,CCnumu_recoNuE_binH));
      h_CCnumu_recoLepP_CV.Fill(overflow(evt.muon.momentum,CCnumu_recoLepP_binH), evt.xsecWeight*(targetPOT/bnb_POT));
      hRawBNB_CCnumu_recoLepP_CV.Fill(overflow(evt.muon.momentum,CCnumu_recoLepP_binH));
      h_CCnumu_recoLepCosTheta_CV.Fill(evt.muon.cosTheta, evt.xsecWeight*(targetPOT/bnb_POT));
      hRawBNB_CCnumu_recoLepCosTheta_CV.Fill(underflow(evt.muon.cosTheta,CCnumu_recoLepCosTheta_binL));
    }
  }

  for(int iB = 1; iB <= h_CCnumu_recoNuE_CV_error.GetNbinsX(); ++iB){
    h_CCnumu_recoNuE_CV_error.SetBinContent(iB, std::sqrt(h_CCnumu_recoNuE_CV_error.GetBinContent(iB)));
  }

  for(const auto& evt : nue_var_selected_events){
    if(!evt.matched) continue;
    if(evt.passedCCnue){
      h_CCnue_recoNuE_var.Fill(underflow(overflow(evt.recoNuE,CCnue_recoNuE_binH),CCnue_recoNuE_binL), evt.xsecWeight*(targetPOT/nue_POT));
      hRawInt_CCnue_recoNuE_var.Fill(underflow(overflow(evt.recoNuE,CCnue_recoNuE_binH),CCnue_recoNuE_binL));
      h_CCnue_recoLepP_var.Fill(overflow(evt.electron.momentum,CCnue_recoLepP_binH), evt.xsecWeight*(targetPOT/nue_POT));
      hRawInt_CCnue_recoLepP_var.Fill(overflow(evt.electron.momentum,CCnue_recoLepP_binH));
      h_CCnue_recoLepCosTheta_var.Fill(underflow(evt.electron.cosTheta,CCnue_recoLepCosTheta_binL), evt.xsecWeight*(targetPOT/nue_POT));
      hRawInt_CCnue_recoLepCosTheta_var.Fill(underflow(evt.electron.cosTheta,CCnue_recoLepCosTheta_binL));
    }
    if(evt.passedCCnumu){
      h_CCnumu_recoNuE_var.Fill(overflow(evt.recoNuE,CCnumu_recoNuE_binH), evt.xsecWeight*(targetPOT/nue_POT));
      hRawInt_CCnumu_recoNuE_var.Fill(overflow(evt.recoNuE,CCnumu_recoNuE_binH));
      h_CCnumu_recoLepP_var.Fill(overflow(evt.muon.momentum,CCnumu_recoLepP_binH), evt.xsecWeight*(targetPOT/nue_POT));
      hRawInt_CCnumu_recoLepP_var.Fill(overflow(evt.muon.momentum,CCnumu_recoLepP_binH));
      h_CCnumu_recoLepCosTheta_var.Fill(evt.muon.cosTheta, evt.xsecWeight*(targetPOT/nue_POT));
      hRawInt_CCnumu_recoLepCosTheta_var.Fill(underflow(evt.muon.cosTheta,CCnumu_recoLepCosTheta_binL));
    }
  }

  float n_bnb_CCnue_var = 0.0;
  for(const auto& evt : bnb_var_selected_events){
    if(!evt.matched) continue;
    if(evt.passedCCnue){
      n_bnb_CCnue_var += evt.xsecWeight*(targetPOT/bnb_POT);
      if(!ccnueBkgHack){
        h_CCnue_recoNuE_var.Fill(underflow(overflow(evt.recoNuE,CCnue_recoNuE_binH),CCnue_recoNuE_binL), evt.xsecWeight*(targetPOT/bnb_POT));
        h_CCnue_recoLepP_var.Fill(overflow(evt.electron.momentum,CCnue_recoLepP_binH), evt.xsecWeight*(targetPOT/bnb_POT));
        h_CCnue_recoLepCosTheta_var.Fill(underflow(evt.electron.cosTheta,CCnue_recoLepCosTheta_binL), evt.xsecWeight*(targetPOT/bnb_POT));
      }
      hRawBNB_CCnue_recoNuE_var.Fill(underflow(overflow(evt.recoNuE,CCnue_recoNuE_binH),CCnue_recoNuE_binL));
      hRawBNB_CCnue_recoLepP_var.Fill(overflow(evt.electron.momentum,CCnue_recoLepP_binH));
      hRawBNB_CCnue_recoLepCosTheta_var.Fill(underflow(evt.electron.cosTheta,CCnue_recoLepCosTheta_binL));
    }
    if(evt.passedCCnumu){
      h_CCnumu_recoNuE_var.Fill(overflow(evt.recoNuE,CCnumu_recoNuE_binH), evt.xsecWeight*(targetPOT/bnb_POT));
      hRawBNB_CCnumu_recoNuE_var.Fill(overflow(evt.recoNuE,CCnumu_recoNuE_binH));
      h_CCnumu_recoLepP_var.Fill(overflow(evt.muon.momentum,CCnumu_recoLepP_binH), evt.xsecWeight*(targetPOT/bnb_POT));
      hRawBNB_CCnumu_recoLepP_var.Fill(overflow(evt.muon.momentum,CCnumu_recoLepP_binH));
      h_CCnumu_recoLepCosTheta_var.Fill(evt.muon.cosTheta, evt.xsecWeight*(targetPOT/bnb_POT));
      hRawBNB_CCnumu_recoLepCosTheta_var.Fill(underflow(evt.muon.cosTheta,CCnumu_recoLepCosTheta_binL));
    }
  }

  if(ccnueBkgHack){
    float n_nue_CCnue_CV = h_CCnue_recoNuE_CV.Integral();
    h_CCnue_recoNuE_CV.Scale((n_nue_CCnue_CV+n_bnb_CCnue_CV)/n_nue_CCnue_CV);
    float n_nue_CCnue_var = h_CCnue_recoNuE_var.Integral();
    h_CCnue_recoNuE_var.Scale((n_nue_CCnue_var+n_bnb_CCnue_var)/n_nue_CCnue_var);

    n_nue_CCnue_CV = h_CCnue_recoLepP_CV.Integral();
    h_CCnue_recoLepP_CV.Scale((n_nue_CCnue_CV+n_bnb_CCnue_CV)/n_nue_CCnue_CV);
    n_nue_CCnue_var = h_CCnue_recoLepP_var.Integral();
    h_CCnue_recoLepP_var.Scale((n_nue_CCnue_var+n_bnb_CCnue_var)/n_nue_CCnue_var);

    n_nue_CCnue_CV = h_CCnue_recoLepCosTheta_CV.Integral();
    h_CCnue_recoLepCosTheta_CV.Scale((n_nue_CCnue_CV+n_bnb_CCnue_CV)/n_nue_CCnue_CV);
    n_nue_CCnue_var = h_CCnue_recoLepCosTheta_var.Integral();
    h_CCnue_recoLepCosTheta_var.Scale((n_nue_CCnue_var+n_bnb_CCnue_var)/n_nue_CCnue_var);
  }


  std::cout << "calculating covariance matrices" << std::endl;

  TH2F CCnue_covar_recoNuE("CCnue_covar_recoNuE","Covariance Matrix",CCnue_recoNuE_nBins,0,CCnue_recoNuE_nBins,CCnue_recoNuE_nBins,0,CCnue_recoNuE_nBins);
  TH2F CCnue_frac_covar_recoNuE("CCnue_frac_covar_recoNuE","Fractional Covariance Matrix",CCnue_recoNuE_nBins,0,CCnue_recoNuE_nBins,CCnue_recoNuE_nBins,0,CCnue_recoNuE_nBins);
  TH2F CCnumu_covar_recoNuE("CCnumu_covar_recoNuE","Covariance Matrix",CCnumu_recoNuE_nBins,0,CCnumu_recoNuE_nBins,CCnumu_recoNuE_nBins,0,CCnumu_recoNuE_nBins);
  TH2F CCnumu_frac_covar_recoNuE("CCnumu_frac_covar_recoNuE","Fractional Covariance Matrix",CCnumu_recoNuE_nBins,0,CCnumu_recoNuE_nBins,CCnumu_recoNuE_nBins,0,CCnumu_recoNuE_nBins);
  TH2F CCnue_covar_recoLepP("CCnue_covar_recoLepP","Covariance Matrix",CCnue_recoLepP_nBins,0,CCnue_recoLepP_nBins,CCnue_recoLepP_nBins,0,CCnue_recoLepP_nBins);
  TH2F CCnue_frac_covar_recoLepP("CCnue_frac_covar_recoLepP","Fractional Covariance Matrix",CCnue_recoLepP_nBins,0,CCnue_recoLepP_nBins,CCnue_recoLepP_nBins,0,CCnue_recoLepP_nBins);
  TH2F CCnumu_covar_recoLepP("CCnumu_covar_recoLepP","Covariance Matrix",CCnumu_recoLepP_nBins,0,CCnumu_recoLepP_nBins,CCnumu_recoLepP_nBins,0,CCnumu_recoLepP_nBins);
  TH2F CCnumu_frac_covar_recoLepP("CCnumu_frac_covar_recoLepP","Fractional Covariance Matrix",CCnumu_recoLepP_nBins,0,CCnumu_recoLepP_nBins,CCnumu_recoLepP_nBins,0,CCnumu_recoLepP_nBins);
  TH2F CCnue_covar_recoLepCosTheta("CCnue_covar_recoLepCosTheta","Covariance Matrix",CCnue_recoLepCosTheta_nBins,0,CCnue_recoLepCosTheta_nBins,CCnue_recoLepCosTheta_nBins,0,CCnue_recoLepCosTheta_nBins);
  TH2F CCnue_frac_covar_recoLepCosTheta("CCnue_frac_covar_recoLepCosTheta","Fractional Covariance Matrix",CCnue_recoLepCosTheta_nBins,0,CCnue_recoLepCosTheta_nBins,CCnue_recoLepCosTheta_nBins,0,CCnue_recoLepCosTheta_nBins);
  TH2F CCnumu_covar_recoLepCosTheta("CCnumu_covar_recoLepCosTheta","Covariance Matrix",CCnumu_recoLepCosTheta_nBins,0,CCnumu_recoLepCosTheta_nBins,CCnumu_recoLepCosTheta_nBins,0,CCnumu_recoLepCosTheta_nBins);
  TH2F CCnumu_frac_covar_recoLepCosTheta("CCnumu_frac_covar_recoLepCosTheta","Fractional Covariance Matrix",CCnumu_recoLepCosTheta_nBins,0,CCnumu_recoLepCosTheta_nBins,CCnumu_recoLepCosTheta_nBins,0,CCnumu_recoLepCosTheta_nBins);

  for(int i = 1; i <= h_CCnue_recoNuE_CV.GetNbinsX(); ++i){
    float x_CV_i = h_CCnue_recoNuE_CV.GetBinContent(i);
    float x_var_i = h_CCnue_recoNuE_var.GetBinContent(i);
    for(int j = i; j <= h_CCnue_recoNuE_CV.GetNbinsX(); ++j){
      float x_CV_j = h_CCnue_recoNuE_CV.GetBinContent(j);
      float x_var_j = h_CCnue_recoNuE_var.GetBinContent(j);
      float variance = (x_var_i - x_CV_i)*(x_var_j - x_CV_j);
      float frac_variance = variance/(x_CV_i*x_CV_j);
      CCnue_covar_recoNuE.SetBinContent(i,j,variance);
      CCnue_frac_covar_recoNuE.SetBinContent(i,j,frac_variance);
      if(j != i){
        CCnue_covar_recoNuE.SetBinContent(j,i,variance);
        CCnue_frac_covar_recoNuE.SetBinContent(j,i,frac_variance);
      }
    }
  }

  for(int i = 1; i <= h_CCnue_recoLepP_CV.GetNbinsX(); ++i){
    float x_CV_i = h_CCnue_recoLepP_CV.GetBinContent(i);
    float x_var_i = h_CCnue_recoLepP_var.GetBinContent(i);
    for(int j = i; j <= h_CCnue_recoLepP_CV.GetNbinsX(); ++j){
      float x_CV_j = h_CCnue_recoLepP_CV.GetBinContent(j);
      float x_var_j = h_CCnue_recoLepP_var.GetBinContent(j);
      float variance = (x_var_i - x_CV_i)*(x_var_j - x_CV_j);
      float frac_variance = variance/(x_CV_i*x_CV_j);
      CCnue_covar_recoLepP.SetBinContent(i,j,variance);
      CCnue_frac_covar_recoLepP.SetBinContent(i,j,frac_variance);
      if(j != i){
        CCnue_covar_recoLepP.SetBinContent(j,i,variance);
        CCnue_frac_covar_recoLepP.SetBinContent(j,i,frac_variance);
      }
    }
  }

  for(int i = 1; i <= h_CCnue_recoLepCosTheta_CV.GetNbinsX(); ++i){
    float x_CV_i = h_CCnue_recoLepCosTheta_CV.GetBinContent(i);
    float x_var_i = h_CCnue_recoLepCosTheta_var.GetBinContent(i);
    for(int j = i; j <= h_CCnue_recoLepCosTheta_CV.GetNbinsX(); ++j){
      float x_CV_j = h_CCnue_recoLepCosTheta_CV.GetBinContent(j);
      float x_var_j = h_CCnue_recoLepCosTheta_var.GetBinContent(j);
      float variance = (x_var_i - x_CV_i)*(x_var_j - x_CV_j);
      float frac_variance = variance/(x_CV_i*x_CV_j);
      CCnue_covar_recoLepCosTheta.SetBinContent(i,j,variance);
      CCnue_frac_covar_recoLepCosTheta.SetBinContent(i,j,frac_variance);
      if(j != i){
        CCnue_covar_recoLepCosTheta.SetBinContent(j,i,variance);
        CCnue_frac_covar_recoLepCosTheta.SetBinContent(j,i,frac_variance);
      }
    }
  }

  for(int i = 1; i <= h_CCnumu_recoNuE_CV.GetNbinsX(); ++i){
    float x_CV_i = h_CCnumu_recoNuE_CV.GetBinContent(i);
    float x_var_i = h_CCnumu_recoNuE_var.GetBinContent(i);
    for(int j = i; j <= h_CCnumu_recoNuE_CV.GetNbinsX(); ++j){
      float x_CV_j = h_CCnumu_recoNuE_CV.GetBinContent(j);
      float x_var_j = h_CCnumu_recoNuE_var.GetBinContent(j);
      float variance = (x_var_i - x_CV_i)*(x_var_j - x_CV_j);
      float frac_variance = variance/(x_CV_i*x_CV_j);
      CCnumu_covar_recoNuE.SetBinContent(i,j,variance);
      CCnumu_frac_covar_recoNuE.SetBinContent(i,j,frac_variance);
      if(j != i){
        CCnumu_covar_recoNuE.SetBinContent(j,i,variance);
        CCnumu_frac_covar_recoNuE.SetBinContent(j,i,frac_variance);
      }
    }
  }

  for(int i = 1; i <= h_CCnumu_recoLepP_CV.GetNbinsX(); ++i){
    float x_CV_i = h_CCnumu_recoLepP_CV.GetBinContent(i);
    float x_var_i = h_CCnumu_recoLepP_var.GetBinContent(i);
    for(int j = i; j <= h_CCnumu_recoLepP_CV.GetNbinsX(); ++j){
      float x_CV_j = h_CCnumu_recoLepP_CV.GetBinContent(j);
      float x_var_j = h_CCnumu_recoLepP_var.GetBinContent(j);
      float variance = (x_var_i - x_CV_i)*(x_var_j - x_CV_j);
      float frac_variance = variance/(x_CV_i*x_CV_j);
      CCnumu_covar_recoLepP.SetBinContent(i,j,variance);
      CCnumu_frac_covar_recoLepP.SetBinContent(i,j,frac_variance);
      if(j != i){
        CCnumu_covar_recoLepP.SetBinContent(j,i,variance);
        CCnumu_frac_covar_recoLepP.SetBinContent(j,i,frac_variance);
      }
    }
  }

  for(int i = 1; i <= h_CCnumu_recoLepCosTheta_CV.GetNbinsX(); ++i){
    float x_CV_i = h_CCnumu_recoLepCosTheta_CV.GetBinContent(i);
    float x_var_i = h_CCnumu_recoLepCosTheta_var.GetBinContent(i);
    for(int j = i; j <= h_CCnumu_recoLepCosTheta_CV.GetNbinsX(); ++j){
      float x_CV_j = h_CCnumu_recoLepCosTheta_CV.GetBinContent(j);
      float x_var_j = h_CCnumu_recoLepCosTheta_var.GetBinContent(j);
      float variance = (x_var_i - x_CV_i)*(x_var_j - x_CV_j);
      float frac_variance = variance/(x_CV_i*x_CV_j);
      CCnumu_covar_recoLepCosTheta.SetBinContent(i,j,variance);
      CCnumu_frac_covar_recoLepCosTheta.SetBinContent(i,j,frac_variance);
      if(j != i){
        CCnumu_covar_recoLepCosTheta.SetBinContent(j,i,variance);
        CCnumu_frac_covar_recoLepCosTheta.SetBinContent(j,i,frac_variance);
      }
    }
  }


  std::cout << "writing histograms to output file" << std::endl;

  TFile* fout = new TFile(outfile.c_str(), "RECREATE");
  fout -> cd();
  h_CCnue_recoNuE_CV.Write("", TObject::kOverwrite);
  h_CCnue_recoLepP_CV.Write("", TObject::kOverwrite);
  h_CCnue_recoLepCosTheta_CV.Write("", TObject::kOverwrite);
  h_CCnue_recoNuE_var.Write("", TObject::kOverwrite);
  h_CCnue_recoLepP_var.Write("", TObject::kOverwrite);
  h_CCnue_recoLepCosTheta_var.Write("", TObject::kOverwrite);
  h_CCnumu_recoNuE_CV.Write("", TObject::kOverwrite);
  h_CCnumu_recoNuE_CV_error.Write("", TObject::kOverwrite);
  h_CCnumu_recoLepP_CV.Write("", TObject::kOverwrite);
  h_CCnumu_recoLepCosTheta_CV.Write("", TObject::kOverwrite);
  h_CCnumu_recoNuE_var.Write("", TObject::kOverwrite);
  h_CCnumu_recoLepP_var.Write("", TObject::kOverwrite);
  h_CCnumu_recoLepCosTheta_var.Write("", TObject::kOverwrite);
  hRawInt_CCnue_recoNuE_CV.Write("", TObject::kOverwrite);
  hRawInt_CCnue_recoLepP_CV.Write("", TObject::kOverwrite);
  hRawInt_CCnue_recoLepCosTheta_CV.Write("", TObject::kOverwrite);
  hRawInt_CCnue_recoNuE_var.Write("", TObject::kOverwrite);
  hRawInt_CCnue_recoLepP_var.Write("", TObject::kOverwrite);
  hRawInt_CCnue_recoLepCosTheta_var.Write("", TObject::kOverwrite);
  hRawInt_CCnumu_recoNuE_CV.Write("", TObject::kOverwrite);
  hRawInt_CCnumu_recoLepP_CV.Write("", TObject::kOverwrite);
  hRawInt_CCnumu_recoLepCosTheta_CV.Write("", TObject::kOverwrite);
  hRawInt_CCnumu_recoNuE_var.Write("", TObject::kOverwrite);
  hRawInt_CCnumu_recoLepP_var.Write("", TObject::kOverwrite);
  hRawInt_CCnumu_recoLepCosTheta_var.Write("", TObject::kOverwrite);
  hRawBNB_CCnue_recoNuE_CV.Write("", TObject::kOverwrite);
  hRawBNB_CCnue_recoLepP_CV.Write("", TObject::kOverwrite);
  hRawBNB_CCnue_recoLepCosTheta_CV.Write("", TObject::kOverwrite);
  hRawBNB_CCnue_recoNuE_var.Write("", TObject::kOverwrite);
  hRawBNB_CCnue_recoLepP_var.Write("", TObject::kOverwrite);
  hRawBNB_CCnue_recoLepCosTheta_var.Write("", TObject::kOverwrite);
  hRawBNB_CCnumu_recoNuE_CV.Write("", TObject::kOverwrite);
  hRawBNB_CCnumu_recoLepP_CV.Write("", TObject::kOverwrite);
  hRawBNB_CCnumu_recoLepCosTheta_CV.Write("", TObject::kOverwrite);
  hRawBNB_CCnumu_recoNuE_var.Write("", TObject::kOverwrite);
  hRawBNB_CCnumu_recoLepP_var.Write("", TObject::kOverwrite);
  hRawBNB_CCnumu_recoLepCosTheta_var.Write("", TObject::kOverwrite);
  CCnue_covar_recoNuE.Write("", TObject::kOverwrite);
  CCnue_frac_covar_recoNuE.Write("", TObject::kOverwrite);
  CCnumu_covar_recoNuE.Write("", TObject::kOverwrite);
  CCnumu_frac_covar_recoNuE.Write("", TObject::kOverwrite);
  CCnue_covar_recoLepP.Write("", TObject::kOverwrite);
  CCnue_frac_covar_recoLepP.Write("", TObject::kOverwrite);
  CCnumu_covar_recoLepP.Write("", TObject::kOverwrite);
  CCnumu_frac_covar_recoLepP.Write("", TObject::kOverwrite);
  CCnue_covar_recoLepCosTheta.Write("", TObject::kOverwrite);
  CCnue_frac_covar_recoLepCosTheta.Write("", TObject::kOverwrite);
  CCnumu_covar_recoLepCosTheta.Write("", TObject::kOverwrite);
  CCnumu_frac_covar_recoLepCosTheta.Write("", TObject::kOverwrite);
  fout -> Close();


  std::cout << "done!" << std::endl;

}


