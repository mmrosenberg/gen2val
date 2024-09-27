
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


struct SelEvent {

  int run, subrun, event;
  float recoNuE, trueNuE;
  int recoMu_simTID;
  float recoMu_muScore;
  float recoMu_piScore;
  float recoMu_prScore;
  float recoMu_muScore_var;
  float recoMu_piScore_var;
  float recoMu_prScore_var;
  bool passedCCnumu;
  bool matched;
  float xsecWeight;

  SelEvent(int r, int sr, int e, float te, float re, float w) : run(r), subrun(sr), event(e), trueNuE(te), recoNuE(re), xsecWeight(w), passedCCnumu(false), matched(false), recoMu_simTID(-9), recoMu_muScore(-99.), recoMu_piScore(-99.), recoMu_prScore(-99.), recoMu_muScore_var(-99.), recoMu_piScore_var(-99.), recoMu_prScore_var(-99.) {}
  SelEvent() : run(-1), subrun(-1), event(-1), trueNuE(-99.), recoNuE(-99.), xsecWeight(1.), passedCCnumu(false), matched(false), recoMu_simTID(-9), recoMu_muScore(-99.), recoMu_piScore(-99.), recoMu_prScore(-99.), recoMu_muScore_var(-99.), recoMu_piScore_var(-99.), recoMu_prScore_var(-99.) {}

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


void checkCCnumuInc(int foundVertex, int vtxIsFiducial, float vtxFracHitsOnCosmic,
 int nTracks, int trackIsSecondary[], int trackClassified[], int trackTrueTID[],
 int trackPID[], float trackMuScore[], float trackPiScore[], float trackPrScore[], SelEvent& evt){

  evt.passedCCnumu = false;

  if(foundVertex == 0 || vtxIsFiducial != 1) return;
  if(vtxFracHitsOnCosmic >= (1. - 1e-6)) return;

  int nMuons = 0;
  for(int iT = 0; iT < nTracks; ++iT){
    if(trackIsSecondary[iT] == 1 || trackClassified[iT] != 1) continue;
    if(trackPID[iT] == 13){
      ++nMuons;
      if(trackMuScore[iT] > evt.recoMu_muScore){
        evt.recoMu_muScore = trackMuScore[iT];
        evt.recoMu_piScore = trackPiScore[iT];
        evt.recoMu_prScore = trackPrScore[iT];
        evt.recoMu_simTID = trackTrueTID[iT];
      }
    }
  }

  if(nMuons < 1) return;

  evt.passedCCnumu = true;
  return;

}


int checkMatchingMuon(int foundVertex, int vtxIsFiducial, float vtxFracHitsOnCosmic,
 int nTracks, int trackClassified[], int trackTrueTID[], float trackTruePurity[],
 float trackMuScore[], float trackPiScore[], float trackPrScore[], SelEvent& evt){

  if(foundVertex == 0 || vtxIsFiducial != 1) return 2;
  if(vtxFracHitsOnCosmic >= (1. - 1e-6)) return 2;
  
  int nMatches = 0;
  float highestPurity = -9.;
  for(int iT = 0; iT < nTracks; ++iT){
    if(trackClassified[iT] != 1) continue;
    if(trackTrueTID[iT] == evt.recoMu_simTID){
      ++nMatches;
      if(trackTruePurity[iT] > highestPurity){
        highestPurity = trackTruePurity[iT];
        evt.recoMu_muScore_var = trackMuScore[iT];
        evt.recoMu_piScore_var = trackPiScore[iT];
        evt.recoMu_prScore_var = trackPrScore[iT];
      }
    }
  }

  if(nMatches < 1) return 1;

  evt.matched = true;
  return 0;

}


int main(int argc, char* argv[]) {

  std::string var = "";
  std::string outfile = "calculate_detvar_systematics_CCnumu_particle_score_shifts_output.root";
  int iarg = 0;

  while(iarg != -1){
    iarg = getopt(argc, argv, "v:o:h");
    switch(iarg){
      case 'v': var = optarg; break;
      case 'o': outfile = optarg; break;
      case 'h': std::cout << "\t-v\tdetector varation to process" << std::endl;
                std::cout << "\t-o\toutput file name" << std::endl;
                return 0;
      case '?': std::cout << "unkown input option, run with -h to see options" << std::endl;
                return 0;
    }
  }

  if(var == ""){
    std::cout << "must provide detector variation. exiting..." << std::endl;
    return 0;
  }

  if(outfile == "calculate_detvar_systematics_CCnumu_particle_score_shifts_output.root"){
    outfile = "calculate_detvar_systematics_CCnumu_particle_score_shifts_"+var+"_output.root";
  }

  TH1::SetDefaultSumw2(kTRUE);


  std::cout << "setting up ntuple TTrees" << std::endl;

  //ntuple file names
  std::string bnb_CV_file = "/home/matthew/microboone/tufts/gen2val/flat_ntuples/detvar/";
  std::string bnb_var_file = "/home/matthew/microboone/tufts/gen2val/flat_ntuples/detvar/";
  if(var == "SCE" || var == "recomb2"){
    bnb_CV_file += "dlgen2_reco_v2me06_ntuple_v7_mcc9_v40a_dl_run3b_bnb_nu_overlay_500k_CV_weightsAdded.root";
    bnb_var_file += "dlgen2_reco_v2me06_ntuple_v7_mcc9_v40a_dl_run3b_bnb_nu_overlay_500k_"+var+"_weightsAdded.root";
  }
  else{
    bnb_CV_file += "dlgen2_reco_v2me06_ntuple_v7_altCVnuWeights_mcc9_v29e_dl_run3b_bnb_nu_overlay_1mil_CV.root";
    bnb_var_file += "dlgen2_reco_v2me06_ntuple_v7_altCVnuWeights_mcc9_v40a_dl_run3b_bnb_nu_overlay_1mil_"+var+".root";
  }

  std::cout << "using bnb CV file: " << bnb_CV_file << std::endl;
  std::cout << "using bnb var file: " << bnb_var_file << std::endl;


  int maxProngs = 100;
  //ntuple variables
  int run, subrun, event;
  int trueNuPDG, trueNuCCNC;
  int foundVertex, vtxIsFiducial;
  float trueNuE, recoNuE;
  float vtxFracHitsOnCosmic;
  float xsecWeight;
  int nTracks;
  int trackIsSecondary[maxProngs];
  int trackClassified[maxProngs];
  int trackPID[maxProngs];
  float trackMuScore[maxProngs];
  float trackPiScore[maxProngs];
  float trackPrScore[maxProngs];
  float trackTruePurity[maxProngs];
  int trackTrueTID[maxProngs];

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
  tree_bnb_CV->SetBranchAddress("trackIsSecondary",trackIsSecondary);
  tree_bnb_CV->SetBranchAddress("trackClassified",trackClassified);
  tree_bnb_CV->SetBranchAddress("trackPID",trackPID);
  tree_bnb_CV->SetBranchAddress("trackMuScore",trackMuScore);
  tree_bnb_CV->SetBranchAddress("trackPiScore",trackPiScore);
  tree_bnb_CV->SetBranchAddress("trackPrScore",trackPrScore);
  tree_bnb_CV->SetBranchAddress("trackTruePurity",trackTruePurity);
  tree_bnb_CV->SetBranchAddress("trackTrueTID",trackTrueTID);

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
  tree_bnb_var->SetBranchAddress("trackIsSecondary",trackIsSecondary);
  tree_bnb_var->SetBranchAddress("trackClassified",trackClassified);
  tree_bnb_var->SetBranchAddress("trackPID",trackPID);
  tree_bnb_var->SetBranchAddress("trackMuScore",trackMuScore);
  tree_bnb_var->SetBranchAddress("trackPiScore",trackPiScore);
  tree_bnb_var->SetBranchAddress("trackPrScore",trackPrScore);
  tree_bnb_var->SetBranchAddress("trackTruePurity",trackTruePurity);
  tree_bnb_var->SetBranchAddress("trackTrueTID",trackTrueTID);


  std::cout << "analyzing events from bnb overlay CV ntuple" << std::endl;

  std::vector<SelEvent> bnb_CV_events;
  bnb_CV_events.reserve(tree_bnb_CV->GetEntries());

  for(int i = 0; i < tree_bnb_CV->GetEntries(); ++i){
    tree_bnb_CV -> GetEntry(i);
    if(std::abs(trueNuPDG) != 14 || trueNuCCNC != 0) continue;
    SelEvent evt(run, subrun, event, trueNuE*1000., recoNuE, xsecWeight);
    checkCCnumuInc(foundVertex, vtxIsFiducial, vtxFracHitsOnCosmic, nTracks, trackIsSecondary,
     trackClassified, trackTrueTID, trackPID, trackMuScore, trackPiScore, trackPrScore, evt);
    //if(evt.passedCCnumu && evt.recoNuE > 200 && evt.recoNuE < 800) bnb_CV_events.push_back(evt);
    bnb_CV_events.push_back(evt);
  }

  std::sort(bnb_CV_events.begin(), bnb_CV_events.end());

  //std::cout << "found "<<bnb_CV_events.size()<<" bnb overlay CV events with 200 < recoNuE < 800 that passed selection" << std::endl;
  std::cout << "found "<<bnb_CV_events.size()<<" bnb overlay CV events" << std::endl;


  std::cout << "analyzing matched events from bnb overlay detvar ntuple" << std::endl;

  //std::vector<SelEvent> bnb_var_selected_events;
  //bnb_var_selected_events.reserve(bnb_CV_events.size());
  int n_bnb_matched = 0;
  int n_CV_sel = 0;
  int n_var_withVtx = 0;
  int n_var_sel = 0;
  float n_var_withVtx_xw = 0.;
  float n_var_sel_xw = 0.;

  for(int i = 0; i < tree_bnb_var->GetEntries(); ++i){
    tree_bnb_var -> GetEntry(i);
    if(std::abs(trueNuPDG) != 14 || trueNuCCNC != 0) continue;
    SelEvent evt(run, subrun, event, trueNuE*1000., recoNuE, xsecWeight);
    auto it_CV = std::lower_bound(bnb_CV_events.begin(), bnb_CV_events.end(), evt);
    if(it_CV == bnb_CV_events.end() || !(*it_CV == evt)) continue;
    ++n_bnb_matched;
    if( !((*it_CV).passedCCnumu && (*it_CV).recoNuE > 200 && (*it_CV).recoNuE < 800) ) continue;
    ++n_CV_sel;
    int muonCheck = checkMatchingMuon(foundVertex, vtxIsFiducial, vtxFracHitsOnCosmic,
     nTracks, trackClassified, trackTrueTID, trackTruePurity,
     trackMuScore, trackPiScore, trackPrScore, *it_CV);
    //if(evt.matched) bnb_var_selected_events.push_back(evt);
    if(muonCheck < 2){ ++n_var_withVtx; n_var_withVtx_xw += xsecWeight; }
    if(muonCheck < 1){ ++n_var_sel; n_var_sel_xw += xsecWeight; }
  }

  std::cout << "found "<<n_bnb_matched <<" events from CV intersection. Found matching muon for "<<n_var_sel<<" of "<<n_CV_sel<<" events that passed CV selection and "<<n_var_withVtx<<" that passed CV and had reco vertex in detvar" << std::endl;
  std::cout << "with xsec weights, "<<((n_var_withVtx_xw - n_var_sel_xw)/n_var_withVtx_xw)*100<<"\% of selected CV events had vertex in detvar but original CV muon track not reconstructed" << std::endl;
  //std::cout << "found "<<bnb_var_selected_events.size()<<" bnb overlay detvar events from CV intersection (total of "<<n_bnb_matched<<" events) that have a matching muon" << std::endl;


  std::cout << "calculating POT scaling factors" << std::endl;

  float totGoodPOT;
  float bnb_POT = 0.;
  float targetPOT = 4.4e+19;

  TTree* bnb_potTree = (TTree*)ntuple_bnb_CV->Get("potTree");
  bnb_potTree -> SetBranchAddress("totGoodPOT",&totGoodPOT);

  for(int i = 0; i < bnb_potTree->GetEntries(); ++i){
    bnb_potTree -> GetEntry(i);
    bnb_POT += totGoodPOT;
  }
  bnb_POT *= ((1.0*n_bnb_matched)/(1.0*bnb_CV_events.size()));


  std::cout << "filling histograms" << std::endl;

  int muS_nBins = 40;
  float muS_binL = -1.5;
  float muS_binH = 0.;

  int piS_nBins = 40;
  float piS_binL = -13.;
  float piS_binH = 0.;

  int prS_nBins = 40;
  float prS_binL = -25.;
  float prS_binH = 0.;

  int piD_nBins = 54;
  float piD_binL = -5.;
  float piD_binH = 13.;

  int prD_nBins = 60;
  float prD_binL = -5.;
  float prD_binH = 25.;

  int muEC_nBins = 1000;
  float muEC_binL = -5.;
  float muEC_binH = 5.;

  int piEC_nBins = 1000;
  float piEC_binL = -5.;
  float piEC_binH = 5.;

  int prEC_nBins = 1000;
  float prEC_binL = -5.;
  float prEC_binH = 5.;

  TH1F h_muS_CV("h_muS_CV","Muon Tracks from Selected CC #nu_{#mu} Events with 0.2 GeV < reco E_{#nu} < 0.8 GeV",muS_nBins,muS_binL,muS_binH);
  TH1F h_muS_var("h_muS_var","Muon Tracks from Selected CC #nu_{#mu} Events with 0.2 GeV < reco E_{#nu} < 0.8 GeV",muS_nBins,muS_binL,muS_binH);
  h_muS_CV.GetXaxis()->SetTitle("log(muon score)");
  h_muS_CV.GetYaxis()->SetTitle("events per 4.4e+19 POT");
  h_muS_var.GetXaxis()->SetTitle("log(muon score)");
  h_muS_var.GetYaxis()->SetTitle("events per 4.4e+19 POT");

  TH1F h_piS_CV("h_piS_CV","Muon Tracks from Selected CC #nu_{#mu} Events with 0.2 GeV < reco E_{#nu} < 0.8 GeV",piS_nBins,piS_binL,piS_binH);
  TH1F h_piS_var("h_piS_var","Muon Tracks from Selected CC #nu_{#mu} Events with 0.2 GeV < reco E_{#nu} < 0.8 GeV",piS_nBins,piS_binL,piS_binH);
  h_piS_CV.GetXaxis()->SetTitle("log(pion score)");
  h_piS_CV.GetYaxis()->SetTitle("events per 4.4e+19 POT");
  h_piS_var.GetXaxis()->SetTitle("log(pion score)");
  h_piS_var.GetYaxis()->SetTitle("events per 4.4e+19 POT");

  TH1F h_prS_CV("h_prS_CV","Muon Tracks from Selected CC #nu_{#mu} Events with 0.2 GeV < reco E_{#nu} < 0.8 GeV",prS_nBins,prS_binL,prS_binH);
  TH1F h_prS_var("h_prS_var","Muon Tracks from Selected CC #nu_{#mu} Events with 0.2 GeV < reco E_{#nu} < 0.8 GeV",prS_nBins,prS_binL,prS_binH);
  h_prS_CV.GetXaxis()->SetTitle("log(proton score)");
  h_prS_CV.GetYaxis()->SetTitle("events per 4.4e+19 POT");
  h_prS_var.GetXaxis()->SetTitle("log(proton score)");
  h_prS_var.GetYaxis()->SetTitle("events per 4.4e+19 POT");

  TH1F h_piD_CV("h_piD_CV","Muon Tracks from Selected CC #nu_{#mu} Events with 0.2 GeV < reco E_{#nu} < 0.8 GeV",piD_nBins,piD_binL,piD_binH);
  TH1F h_piD_var("h_piD_var","Muon Tracks from Selected CC #nu_{#mu} Events with 0.2 GeV < reco E_{#nu} < 0.8 GeV",piD_nBins,piD_binL,piD_binH);
  h_piD_CV.GetXaxis()->SetTitle("log(muon score) - log(pion score)");
  h_piD_CV.GetYaxis()->SetTitle("events per 4.4e+19 POT");
  h_piD_var.GetXaxis()->SetTitle("log(muon score) - log(pion score)");
  h_piD_var.GetYaxis()->SetTitle("events per 4.4e+19 POT");

  TH1F h_prD_CV("h_prD_CV","Muon Tracks from Selected CC #nu_{#mu} Events with 0.2 GeV < reco E_{#nu} < 0.8 GeV",prD_nBins,prD_binL,prD_binH);
  TH1F h_prD_var("h_prD_var","Muon Tracks from Selected CC #nu_{#mu} Events with 0.2 GeV < reco E_{#nu} < 0.8 GeV",prD_nBins,prD_binL,prD_binH);
  h_prD_CV.GetXaxis()->SetTitle("log(muon score) - log(proton score)");
  h_prD_CV.GetYaxis()->SetTitle("events per 4.4e+19 POT");
  h_prD_var.GetXaxis()->SetTitle("log(muon score) - log(proton score)");
  h_prD_var.GetYaxis()->SetTitle("events per 4.4e+19 POT");

  TH1F h_muEC("h_muEC","Muon Tracks from Selected CC #nu_{#mu} Events with 0.2 GeV < reco E_{#nu} < 0.8 GeV",muEC_nBins,muEC_binL,muEC_binH);
  h_muEC.GetXaxis()->SetTitle("(log(CV muon score) - log(var muon score))/log(CV muon score)");
  h_muEC.GetYaxis()->SetTitle("events per 4.4e+19 POT");

  TH1F h_piEC("h_piEC","Muon Tracks from Selected CC #nu_{#mu} Events with 0.2 GeV < reco E_{#nu} < 0.8 GeV",piEC_nBins,piEC_binL,piEC_binH);
  h_piEC.GetXaxis()->SetTitle("(log(CV pion score) - log(var pion score))/log(CV pion score)");
  h_piEC.GetYaxis()->SetTitle("events per 4.4e+19 POT");

  TH1F h_prEC("h_prEC","Muon Tracks from Selected CC #nu_{#mu} Events with 0.2 GeV < reco E_{#nu} < 0.8 GeV",prEC_nBins,prEC_binL,prEC_binH);
  h_prEC.GetXaxis()->SetTitle("(log(CV proton score) - log(var proton score))/log(CV proton score)");
  h_prEC.GetYaxis()->SetTitle("events per 4.4e+19 POT");


  float n_var_sel_w = 0.;
  float n_stillMu_w = 0.;
  float n_different_w = 0.;

  for(const auto& evt : bnb_CV_events){
    if(!evt.matched) continue;
    h_muS_CV.Fill(evt.recoMu_muScore, evt.xsecWeight*(targetPOT/bnb_POT));
    h_piS_CV.Fill(evt.recoMu_piScore, evt.xsecWeight*(targetPOT/bnb_POT));
    h_prS_CV.Fill(evt.recoMu_prScore, evt.xsecWeight*(targetPOT/bnb_POT));
    h_piD_CV.Fill(evt.recoMu_muScore - evt.recoMu_piScore, evt.xsecWeight*(targetPOT/bnb_POT));
    h_prD_CV.Fill(evt.recoMu_muScore - evt.recoMu_prScore, evt.xsecWeight*(targetPOT/bnb_POT));
    h_muS_var.Fill(evt.recoMu_muScore_var, evt.xsecWeight*(targetPOT/bnb_POT));
    h_piS_var.Fill(evt.recoMu_piScore_var, evt.xsecWeight*(targetPOT/bnb_POT));
    h_prS_var.Fill(evt.recoMu_prScore_var, evt.xsecWeight*(targetPOT/bnb_POT));
    h_piD_var.Fill(evt.recoMu_muScore_var - evt.recoMu_piScore_var, evt.xsecWeight*(targetPOT/bnb_POT));
    h_prD_var.Fill(evt.recoMu_muScore_var - evt.recoMu_prScore_var, evt.xsecWeight*(targetPOT/bnb_POT));
    h_muEC.Fill((evt.recoMu_muScore - evt.recoMu_muScore_var)/evt.recoMu_muScore, evt.xsecWeight*(targetPOT/bnb_POT));
    h_piEC.Fill((evt.recoMu_piScore - evt.recoMu_piScore_var)/evt.recoMu_piScore, evt.xsecWeight*(targetPOT/bnb_POT));
    h_prEC.Fill((evt.recoMu_prScore - evt.recoMu_prScore_var)/evt.recoMu_prScore, evt.xsecWeight*(targetPOT/bnb_POT));
    n_var_sel_w += evt.xsecWeight*(targetPOT/bnb_POT);
    if(evt.recoMu_muScore_var > evt.recoMu_piScore_var && evt.recoMu_muScore_var > evt.recoMu_prScore_var) n_stillMu_w += evt.xsecWeight*(targetPOT/bnb_POT);
    else n_different_w += evt.xsecWeight*(targetPOT/bnb_POT);
  }

  std::cout << "weighted to 4.4e+19 POT, of "<<n_var_sel_w<<" selected and matched events, original muon is still classified as muon in "<<n_stillMu_w<<" and classified as proton or pion in "<<n_different_w<<" (misclassified in "<<(n_different_w/n_var_sel_w)*100<<"\% of events)" << std::endl;

  /*for(const auto& evt : bnb_var_selected_events){
    if(!evt.matched) continue;
    h_muS_var.Fill(evt.recoMu_muScore, evt.xsecWeight*(targetPOT/bnb_POT));
    h_piS_var.Fill(evt.recoMu_piScore, evt.xsecWeight*(targetPOT/bnb_POT));
    h_prS_var.Fill(evt.recoMu_prScore, evt.xsecWeight*(targetPOT/bnb_POT));
    h_piD_var.Fill(evt.recoMu_muScore - evt.recoMu_piScore, evt.xsecWeight*(targetPOT/bnb_POT));
    h_prD_var.Fill(evt.recoMu_muScore - evt.recoMu_prScore, evt.xsecWeight*(targetPOT/bnb_POT));
  }*/


  std::cout << "writing histograms to output file" << std::endl;

  TFile* fout = new TFile(outfile.c_str(), "RECREATE");
  fout -> cd();
  h_muS_CV.Write();
  h_piS_CV.Write();
  h_prS_CV.Write();
  h_piD_CV.Write();
  h_prD_CV.Write();
  h_muS_var.Write();
  h_piS_var.Write();
  h_prS_var.Write();
  h_piD_var.Write();
  h_prD_var.Write();
  h_muEC.Write();
  h_piEC.Write();
  h_prEC.Write();
  fout -> Close();


  std::cout << "done!" << std::endl;

}


