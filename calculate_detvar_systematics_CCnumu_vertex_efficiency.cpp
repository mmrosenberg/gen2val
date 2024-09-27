
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
  float trueNuE;
  bool passedVertex;
  bool matched;
  float xsecWeight;

  SelEvent(int r, int sr, int e, float te, float w, bool pv, bool m) : run(r), subrun(sr), event(e), trueNuE(te), xsecWeight(w), passedVertex(pv), matched(m) {}
  SelEvent() : run(-1), subrun(-1), event(-1), trueNuE(-99.), xsecWeight(1.), passedVertex(false), matched(false) {}

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


bool vertexReconstructed(int foundVertex, int vtxIsFiducial, float vtxFracHitsOnCosmic, float vtxDistToTrue){

  if(foundVertex == 0 || vtxIsFiducial != 1) return false;
  if(vtxFracHitsOnCosmic >= (1. - 1e-6)) return false;
  if(vtxDistToTrue > 3.) return false;
  return true;

}



int main(int argc, char* argv[]) {

  std::string var = "";
  std::string outfile = "calculate_detvar_systematics_CCnumu_vertex_efficiency_output.root";
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

  if(outfile == "calculate_detvar_systematics_CCnumu_vertex_efficiency_output.root"){
    outfile = "calculate_detvar_systematics_CCnumu_vertex_efficiency_"+var+"_output.root";
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

  //ntuple variables
  int run, subrun, event;
  int trueNuPDG, trueNuCCNC;
  int foundVertex, vtxIsFiducial;
  float trueNuE; //, recoNuE;
  float vtxFracHitsOnCosmic, vtxDistToTrue;
  float xsecWeight;

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
  tree_bnb_CV->SetBranchAddress("vtxFracHitsOnCosmic",&vtxFracHitsOnCosmic);
  tree_bnb_CV->SetBranchAddress("vtxDistToTrue",&vtxDistToTrue);

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
  tree_bnb_var->SetBranchAddress("vtxFracHitsOnCosmic",&vtxFracHitsOnCosmic);
  tree_bnb_var->SetBranchAddress("vtxDistToTrue",&vtxDistToTrue);


  std::cout << "analyzing events from bnb overlay CV ntuple" << std::endl;

  std::vector<SelEvent> bnb_CV_events;
  bnb_CV_events.reserve(tree_bnb_CV->GetEntries());

  for(int i = 0; i < tree_bnb_CV->GetEntries(); ++i){
    tree_bnb_CV -> GetEntry(i);
    if(std::abs(trueNuPDG) != 14 || trueNuCCNC != 0) continue;
    bool passedVertex = vertexReconstructed(foundVertex, vtxIsFiducial, vtxFracHitsOnCosmic, vtxDistToTrue);
    SelEvent evt(run, subrun, event, trueNuE*1000., xsecWeight, passedVertex, false);
    bnb_CV_events.push_back(evt);
  }

  std::sort(bnb_CV_events.begin(), bnb_CV_events.end());

  std::cout << "looked at "<<bnb_CV_events.size()<<" bnb overlay CV events" << std::endl;


  std::cout << "analyzing matched events from bnb overlay detvar ntuple" << std::endl;

  std::vector<SelEvent> bnb_var_selected_events;
  bnb_var_selected_events.reserve(bnb_CV_events.size());
  int n_bnb_matched = 0;

  for(int i = 0; i < tree_bnb_var->GetEntries(); ++i){
    tree_bnb_var -> GetEntry(i);
    if(std::abs(trueNuPDG) != 14 || trueNuCCNC != 0) continue;
    SelEvent evt(run, subrun, event, trueNuE*1000., xsecWeight, false, false);
    auto it_CV = std::lower_bound(bnb_CV_events.begin(), bnb_CV_events.end(), evt);
    if(it_CV == bnb_CV_events.end() || !(*it_CV == evt)) continue;
    (*it_CV).matched = true; evt.matched = true; ++n_bnb_matched;
    bool passedVertex = vertexReconstructed(foundVertex, vtxIsFiducial, vtxFracHitsOnCosmic, vtxDistToTrue);
    if(passedVertex){
      evt.passedVertex = true;
      bnb_var_selected_events.push_back(evt);
    }
  }

  std::cout << "found "<<bnb_var_selected_events.size()<<" bnb overlay detvar events from CV intersection (total of "<<n_bnb_matched<<" events) that pass the selection" << std::endl;


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

  int nBins = 28;
  float binL = 0.2;
  float binH = 3.;

  TH1F h_trueNuE_CV("h_trueNuE_CV","True CCnumu Events with Reconstructed Vertex", nBins, binL, binH);
  TH1F h_trueNuE_var("h_trueNuE_var","True CCnumu Events with Reconstructed Vertex", nBins, binL, binH);

  for(const auto& evt : bnb_CV_events){
    if(!evt.matched || !evt.passedVertex) continue;
    h_trueNuE_CV.Fill(evt.trueNuE/1000., evt.xsecWeight*(targetPOT/bnb_POT));
  }

  for(const auto& evt : bnb_var_selected_events){
    if(!evt.matched || !evt.passedVertex) continue;
    h_trueNuE_var.Fill(evt.trueNuE/1000., evt.xsecWeight*(targetPOT/bnb_POT));
  }


  std::cout << "writing histograms to output file" << std::endl;

  TFile* fout = new TFile(outfile.c_str(), "RECREATE");
  fout -> cd();
  h_trueNuE_CV.Write();
  h_trueNuE_var.Write();
  fout -> Close();


  std::cout << "done!" << std::endl;

}


