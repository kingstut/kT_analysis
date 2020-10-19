/*
#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include "RooUnfoldResponse.h"
#include <RooUnfoldBayes.h>
#include <RooUnfoldSvd.h>
#include <RooUnfoldTUnfold.h>
#endif
*/

#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>

#include <TROOT.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TH2D.h>
#include <THStack.h>
#include <TProfile.h>
#include <iostream>
#include <fstream>
#include <TGraph.h>
#include <TF1.h>
#include <TRandom.h>
#include <TLegend.h>

//gSystem->Load("/global/homes/s/stuti/jets/jet-jet/RooUnfold/libRooUnfold");
//TSystem::Load("/global/homes/s/stuti/jets/jet-jet/RooUnfold/libRooUnfold")
//gSystem->Load("/global/homes/s/stuti/RooUnfold/libRooUnfold");


#include "RooUnfoldResponse.h"
//#include "/global/homes/s/stuti/jets/jet-jet/RooUnfold/src/RooUnfoldBayes.h"
#include "RooUnfoldSvd.h"
#include "RooUnfoldTUnfold.h"


#define NTRACK_MAX (1U << 14)

#include <vector>
#include <math.h>
#include <set>

#include <TSystem.h>
//gSystem->Load("/global/homes/s/stuti/jets/jet-jet/RooUnfold/libRooUnfold");


const int MAX_INPUT_LENGTH = 200;
int probability[19]= {8060000, 877054, 162159, 5230, 563197, 406269, 319087, 231905, 179595, 127286, 92413, 74976, 57540, 40103, 22667, 5230, 5000, 4000};
//int probability[19]= {766,438, 87, 8, 576,433, 336,260,200,157,123,97,75,60,45,37,28,22,18};


float prob_pt( float p){
 float x= (probability[0]+probability[1]+ probability[2]+ probability[3]+ probability[4]+probability[5]+ probability[6]+probability[7]+probability[8]+probability[9]+probability[10]+probability[11]+probability[12]+probability[13]+probability[14]+probability[15]+probability[16]+probability[17]+probability[18])*gRandom->Uniform();
          if (x<probability[0])
              return p;
          else if (x< (probability[0]+probability[1]))
              return p+1;
           else if (x< (probability[0]+probability[1]+probability[2]))
               return p+2;
           else if (x< (probability[0]+probability[1]+probability[2]+probability[3]))
               return p+3;
           else if (x< (probability[0]+probability[1]+ probability[2]+ probability[3]+ probability[4]))
               return p-1;
           else if (x< (probability[0]+probability[1]+ probability[2]+ probability[3]+ probability[4]+probability[5]))
               return p-2;
           else if (x< (probability[0]+probability[1]+ probability[2]+ probability[3]+ probability[4]+probability[5]+ probability[6]))
               return p-3;
           else if (x< (probability[0]+probability[1]+ probability[2]+ probability[3]+ probability[4]+probability[5]+ probability[6]+probability[7]))
               return p-4;
                     else if (x< (probability[0]+probability[1]+ probability[2]+ probability[3]+ probability[4]+probability[5]+ probability[6]+probability[7]+probability[8]))
              return p-5;
                     else if (x< (probability[0]+probability[1]+ probability[2]+ probability[3]+ probability[4]+probability[5]+ probability[6]+probability[7]+probability[8]+probability[9]))
               return p-6;
                     else if (x< (probability[0]+probability[1]+ probability[2]+ probability[3]+ probability[4]+probability[5]+ probability[6]+probability[7]+probability[8]+probability[9]+probability[10]))
              return p-7;
                     else if (x< (probability[0]+probability[1]+ probability[2]+ probability[3]+ probability[4]+probability[5]+ probability[6]+probability[7]+probability[8]+probability[9]+probability[10]+probability[11]))
               return p-probability[3];
                     else if (x< (probability[0]+probability[1]+ probability[2]+ probability[3]+ probability[4]+probability[5]+ probability[6]+probability[7]+probability[8]+probability[9]+probability[10]+probability[11]+probability[12]))
               return p-9;
                               else if (x< (probability[0]+probability[1]+ probability[2]+ probability[3]+ probability[4]+probability[5]+ probability[6]+probability[7]+probability[8]+probability[9]+probability[10]+probability[11]+probability[12]+probability[13]))
               return p-10;
                               else if (x< (probability[0]+probability[1]+ probability[2]+ probability[3]+ probability[4]+probability[5]+ probability[6]+probability[7]+probability[8]+probability[9]+probability[10]+probability[11]+probability[12]+probability[13]+probability[14]))
               return p-11;
                               else if (x< (probability[0]+probability[1]+ probability[2]+ probability[3]+ probability[4]+probability[5]+ probability[6]+probability[7]+probability[8]+probability[9]+probability[10]+probability[11]+probability[12]+probability[13]+probability[14]+probability[15]))
             return p-12;
                               else if (x< (probability[0]+probability[1]+ probability[2]+ probability[3]+ probability[4]+probability[5]+ probability[6]+probability[7]+probability[8]+probability[9]+probability[10]+probability[11]+probability[12]+probability[13]+probability[14]+probability[15]+probability[16]))
             return p-13;
                                         else if (x< (probability[0]+probability[1]+ probability[2]+ probability[3]+ probability[4]+probability[5]+ probability[6]+probability[7]+probability[8]+probability[9]+probability[10]+probability[11]+probability[12]+probability[13]+probability[14]+probability[15]+probability[16]+probability[17]))
               return p-14;
                                         else
               return p-15;


}






int main(int argc, char *argv[]) {
  if (argc < 2) {
    fprintf(stderr,"Format: [command] [root file] \n");
    exit(EXIT_FAILURE);
  }

  //Config File
  FILE* config = fopen("/global/homes/s/stuti/jets/Corr_config.yaml", "r");
  double pT_min = 0;
  double pT_max = 0;
  double Eta_max = 0;

  bool do_pile = false;

  float track_pT_min = 0.0;
  float track_pT_max = 0.0;
  int Track_Cut_Bit = 0;
  int track_chi_max = 0;
  int n_eta_bins = 0;
  int n_phi_bins = 0;

  // zT bins
  //FIXME: Will have to likely set nztbins first, then initialize array
  int nztbins = 7;
  float* ztbins;
  ztbins = new float[nztbins+1];
  ztbins[0] = 0.0; ztbins[1] = 0.1; ztbins[2] = 0.2; ztbins[3] = 0.4; ztbins[4] = 0.6; ztbins[5] = 0.8; ztbins[6] = 1.0; ztbins[7] = 1.2;

  int nptbins = 3;
  float* ptbins;
  ptbins = new float[nptbins+1];
  ptbins[0] = 10.0; ptbins[1] = 11; ptbins[2] = 12.5; ptbins[3] = 16;

  //FIXME: Obviously needs to be put in a header file.
  char line[MAX_INPUT_LENGTH];
  while (fgets(line, MAX_INPUT_LENGTH, config) != NULL) {
      if (line[0] == '#') continue;

      char key[MAX_INPUT_LENGTH];
      char dummy[MAX_INPUT_LENGTH];
      char value[MAX_INPUT_LENGTH];

      // Cap off key[0] and value[0] with null characters and load the key, dummy-characters, and value of the line into their respective arrays
      key[0] = '\0';
      value[0] = '\0';
      sscanf(line, "%[^:]:%[ \t]%100[^\n]", key, dummy, value);


      if (strcmp(key, "pT_min") == 0) {
          pT_min = atof(value);
          std::cout << "pT_min: " << pT_min << std::endl; }

      else if (strcmp(key, "pT_max") == 0) {
          pT_max = atof(value);
          std::cout << "pT_max: " << pT_max << std::endl; }

      else if (strcmp(key, "Eta_max") == 0) {
          Eta_max = atof(value);
          std::cout << "Eta_max: " << Eta_max << std::endl; }

      else if (strcmp(key, "do_pileup_cut") == 0) {
	if (strcmp(value,"true") == 0)
	  do_pile = true;
	std::cout << "do_pileup_cut: " << do_pile << std::endl; }

      else if (strcmp(key, "N_Phi_Bins") == 0) {
	n_phi_bins = atoi(value);
	std::cout << "Number of Phi Bins: " << n_phi_bins << std::endl; }

      else if (strcmp(key, "N_Eta_Bins") == 0) {
        n_eta_bins = atoi(value);
	std::cout << "Number of Eta Bins: " << n_eta_bins << std::endl; }

      else if (strcmp(key, "Track_pT_Min") == 0) {
          track_pT_min = atof(value);
          std::cout << "Track Min pT: " << track_pT_min << std::endl; }

      else if (strcmp(key, "Track_pT_Max") == 0) {
          track_pT_max = atof(value);
          std::cout << "Track Max pT: " << track_pT_max << std::endl; }

      else if (strcmp(key, "Track_Cut_Bit") == 0) {
          Track_Cut_Bit = atoi(value);
          std::cout << "Track Cut Bit: " << Track_Cut_Bit << std::endl; }

      else if (strcmp(key, "Track_Chi_Max") == 0) {
          track_chi_max = atoi(value);
          std::cout << "Track Chi^2 Max: " << track_chi_max << std::endl; }

      else if (strcmp(key, "Zt_bins") == 0) {
          nztbins =-1;
          for (const char *v = value; *v != ']';) {
              while (*v != ']' && !isdigit(*v)) v++;
	      nztbins++;
              while (*v != ']' && (isdigit(*v) || *v == '.')) v++; }

          ztbins = new float[nztbins + 1];
          int i = 0;
          for (const char *v = value; *v != ']' ;) {
              while (*v != ']' && !isdigit(*v)) v++;
              ztbins[i] = atof(v);
              i++;
              while (*v != ']' && (isdigit(*v) || *v == '.')) v++; }

          std::cout << "Number of zT bins: " << nztbins << std::endl << "zT bins: {";
          for (int i = 0; i <= nztbins; i++)
	    std::cout << ztbins[i] << ", ";
          std::cout << "}\n";
      }

      else if (strcmp(key, "Pt_bins") == 0) {
          nptbins =-1;
          for (const char *v = value; *v != ']';) {
            while (*v != ']' && !isdigit(*v)) v++;
	    nptbins++;
	    while (*v != ']' && (isdigit(*v) || *v == '.')) v++; }

          ptbins = new float[nptbins + 1];
	  int i = 0;
	  for (const char *v = value; *v != ']' ;) {
	     while (*v != ']' && !isdigit(*v))  v++;
	     ptbins[i] = atof(v);
	     i++;
	     while (*v != ']' && (isdigit(*v) || *v == '.')) v++; }

	  std::cout << "Number of Pt bins: " << nptbins << std::endl << "Pt bins: {";
          for (int i = 0; i <= nptbins; i++)
	    std::cout << ptbins[i] << ", ";
	  std::cout << "}\n";
      }
  }
  //end Config Loop

  fclose(config);

  for (int i = 0; i <= nztbins; i++)
    std::cout << "zt bound: " << ztbins[i] << std::endl;
  for (int i = 0; i <= nptbins; i++)
    std::cout << "pt bound: " << ptbins[i] << std::endl;


  //HISTOGRAMS
  //TCanvas canvas("canvas", "");
  //TCanvas *canvas = new TCanvas("canvas","multipads",900,700);
    TH1D * h_mainjetpt;
    TH1D * h_mainjetpt_toy;
 THStack * hs_unf;
    
    TH1D * kt_20_40;
    TH1D * kt_40_60;
    TH1D * kt_20_40_truth;
    TH1D * kt_40_60_truth;
    TH1D * corr_20_40;
    TH1D * corr_40_60;
    
    TH2D * kt_pt;
    TH2D * kt_pt_truth;
    THStack * hs;
    TH2D* pt_comp;
    
    
    
    TH1D * h_pt;
    TH1D * h_pt_truth;
    TH1D * h_delta_phi;
    TH1D * h_delta_phi_truth;
    
    
  //INITIALIZE HISTOGRAMS
          h_mainjetpt = new TH1D("mainjetpt","Paired Jet pt distribution", 30,0,60);
          h_mainjetpt_toy = new TH1D("mainjetpt_toy","Paired Jet pt distribution", 60,0,60);
    
    h_delta_phi = new TH1D("delphi","Paired Jet pt distribution",30,0,7);
    h_delta_phi_truth = new TH1D("delphi","Paired Jet pt distribution",30,0,7);
    
    
    hs_unf= new THStack("unf_pT","corrected pT");


          
          

        int  ipt =0;
        int izt=0;

          
     h_pt  = new TH1D("pt","Jet pt distribution", n_phi_bins*2,0,60);
          
      h_pt_truth  = new TH1D("pt_truth","Jet pt_truth distribution", n_phi_bins*2,0,60);

          
            kt_20_40  = new TH1D("20_40_kt","Jet kT distribution 20-40", 8,0,70);
          
                kt_40_60  = new TH1D("40_60_kt","Jet kT distribution 40-60", 8,0,70);
          
    const Int_t NBINS_20 = 5;
   Double_t edges_20[NBINS_20 + 1] = {0.0, 4.0, 8.0, 16.0, 24.0, 32.0};
    
        const Int_t NBINS_40 = 7;
   Double_t edges_40[NBINS_40 + 1] = {0.0, 4.0, 8.0, 16.0, 24.0, 32.0, 40.0, 48.0};
    
                
            kt_20_40_truth  = new TH1D("20_40_kt_truth","Truth Jet kT distribution 20-40", NBINS_20, edges_20);
          
                kt_40_60_truth  = new TH1D("40_60_kt_truth","Truth Jet kT distribution 40-60", NBINS_40, edges_40);    
          
          
          
                      corr_20_40  = new TH1D("20_40correction","Jet kT distribution 20-40", 8,0,70);
          
                corr_40_60  = new TH1D("40_60correction","Jet kT distribution 40-60", n_phi_bins*2,0,70);
          
          
          
      kt_pt  = new TH2D("_KT","kt vs pt", 30,20,70,30, 20, 70);
      kt_pt->Sumw2();
          
      kt_pt_truth  = new TH2D("KT__pT_truth","kt vs pt truth", 30,20,70,30, 20, 70);
      kt_pt_truth->Sumw2();
          
          
      pt_comp  = new TH2D("pT_compare","pt truth vs raw", 15,20,70,15, 20, 70);
      pt_comp->Sumw2();    
          
          
          
      hs = new THStack("stHistogram_5_4","Stacked 2D histograms");
      
          


          kt_20_40->Sumw2();
        kt_40_60->Sumw2();
          
          
          
                    corr_20_40->Sumw2();
        corr_40_60->Sumw2();
          
          
          
                    kt_20_40_truth->Sumw2();
        kt_40_60_truth->Sumw2();
          


      //h_mainjetpt->Sumw2();  
          h_pt->Sumw2(); 
          h_pt_truth->Sumw2(); 
          


    int iarg = 1;
    TString root_file = (TString)argv[iarg];
    std::cout << "Opening: " << (TString)argv[iarg] << std::endl;

    TFile *file = TFile::Open(root_file);

    if (file == NULL) {
      std::cout << " fail" << std::endl;
      exit(EXIT_FAILURE);
    }
    file->Print();

    TTree *_tree_event = dynamic_cast<TTree *>(file->Get("_tree_event"));

    if (_tree_event == NULL) {
      _tree_event = dynamic_cast<TTree *>(file->Get("AliAnalysisTaskNTGJ/_tree_event"));
      if (_tree_event == NULL) {
	std::cout << " fail " << std::endl;
	exit(EXIT_FAILURE);
      }
    }

    //Events
    Bool_t is_pileup_from_spd_5_08;
    Double_t primary_vertex[3];
    Float_t ue_estimate_its_const;
    //Float_t ue_estimate_tpc_const;
    ULong64_t trigger_mask[2];  
    Float_t centrality_v0m;


    // Set the branch addresses of the branches in the TTrees
    _tree_event->SetBranchStatus("*mc*", 0);

    //event Addresses
    _tree_event->SetBranchAddress("centrality_v0m", &centrality_v0m);
    _tree_event->SetBranchAddress("trigger_mask", trigger_mask);
    _tree_event->SetBranchAddress("primary_vertex", primary_vertex);
    _tree_event->SetBranchAddress("is_pileup_from_spd_5_08", &is_pileup_from_spd_5_08);


    //track Addresses
    _tree_event->SetBranchAddress("primary_vertex", primary_vertex);

    //Jet Addresses
    UInt_t njet;
    Float_t jet_e[NTRACK_MAX];
    Float_t jet_pt[NTRACK_MAX];
    Float_t jet_eta[NTRACK_MAX];
    Float_t jet_phi[NTRACK_MAX];
    Float_t jet_multiplicity[NTRACK_MAX];
    Float_t jet_area[NTRACK_MAX];
    
    
    UInt_t njet_truth;
    Float_t jet_e_truth[NTRACK_MAX];
    Float_t jet_pt_truth[NTRACK_MAX];
    Float_t jet_eta_truth[NTRACK_MAX];
    Float_t jet_phi_truth[NTRACK_MAX];


    _tree_event->SetBranchAddress("njet_ak04tpc", &njet);
    _tree_event->SetBranchAddress("jet_ak04tpc_e_raw", jet_e);
    _tree_event->SetBranchAddress("jet_ak04tpc_pt_raw", jet_pt);
    _tree_event->SetBranchAddress("jet_ak04tpc_eta_raw", jet_eta);
    _tree_event->SetBranchAddress("jet_ak04tpc_phi", jet_phi);
    _tree_event->SetBranchAddress("njet_truth_ak04", &njet_truth);
    _tree_event->SetBranchAddress("jet_truth_ak04_pt", jet_pt_truth);
    _tree_event->SetBranchAddress("jet_truth_ak04_eta", jet_eta_truth);
    _tree_event->SetBranchAddress("jet_truth_ak04_phi", jet_phi_truth);
    _tree_event->SetBranchAddress("jet_ak04tpc_area", jet_area);
    
    int num_truth=0;
    int num_raw =0;
    
    //Raw used, Jets not yet callibrated

    Long64_t nentries = _tree_event->GetEntries();
    std::cout << " Total Number of entries in TTree: " << nentries << std::endl;
    Float_t accepted_events = 0;
    
    //ULong64_t myTriggers = (1 << 13) | (1 << 14) | (1 << 15) | (1 << 16);
    
    
    TH1D *density = new TH1D("density", "density", 100, 0, 20);
    TH1D *fluct = new TH1D("fluct", "fluct", 40,-15, 20);
    TH1D *fluct_non_zero = new TH1D("fluct_non_zero", "fluct_non_zero", 40,-15, 20);
    
  TH2D* hRes= new TH2D ("true", "Test Truth",    30, 0, 60, 30, 0, 60);
  TH1D* hMeas= new TH1D ("meas", "Test Measured", 30, 0, 60);
    TH1D* hTrue= new TH1D ("tru", "Test True", 30, 0, 60);

    
    //Unfolding
   // RooUnfoldResponse response (55, 5, 60);
    
    for(int i =0; i<=60; i++){
        for(int j=0; j<probability[0]; j++) //110
        {hRes->Fill (i, i);
         hMeas->Fill (i);
         hTrue->Fill (i);
        //response.Fill(i, i);
        }  
      for(int j=0; j<probability[1]; j++)    //76   
        {hRes->Fill (i, i+1);
        hMeas->Fill (i);
         hTrue->Fill (i+1);
        //response.Fill(i, i+1);
        }  
    for(int j=0; j<probability[2]; j++)    //12
        {hRes->Fill (i, i+2);
         hMeas->Fill (i);
         hTrue->Fill (i+2);
        //response.Fill(i, i+2);
        }   
      for(int j=0; j<probability[3]; j++)     //1
        {hRes->Fill (i, i+3);
                  hMeas->Fill (i);
         hTrue->Fill (i+3);
        //response.Fill(i, i+3);
        }  
    for(int j=0; j<probability[4]; j++) //94
    {hRes->Fill (i, i-1);
                 hMeas->Fill (i);
         hTrue->Fill (i-1);
        //response.Fill(i, i-1);
    }
    for(int j=0; j<probability[5]; j++) //76
    {hRes->Fill (i, i-2);
              hMeas->Fill (i);
         hTrue->Fill (i-2);
       // response.Fill(i, i-2);
    }
    for(int j=0; j<probability[6]; j++)  //4probability[3]
    {hRes->Fill (i, i-3);
              hMeas->Fill (i);
         hTrue->Fill (i-3);
        //response.Fill(i, i-3);
    }
    for(int j=0; j<probability[7]; j++)  //3probability[3]
    {hRes->Fill (i, i-4);
              hMeas->Fill (i);
         hTrue->Fill (i-4);
        //response.Fill(i, i-4);
    }     
    for(int j=0; j<probability[8]; j++)  //30
    {hRes->Fill (i, i-5);
              hMeas->Fill (i);
         hTrue->Fill (i-5);
       // response.Fill(i, i-5);
    }
    for(int j=0; j<probability[9]; j++)  //19
    {hRes->Fill (i, i-6);
              hMeas->Fill (i);
         hTrue->Fill (i-6);
        //response.Fill(i, i-6);
    }
    for(int j=0; j<probability[10]; j++)  //15
    {hRes->Fill (i, i-7);
              hMeas->Fill (i);
         hTrue->Fill (i-7);
       // response.Fill(i, i-7);
    }
            for(int j=0; j<probability[11]; j++)  //13
    {hRes->Fill (i, i-8);
                   hMeas->Fill (i);
         hTrue->Fill (i-8);
       // response.Fill(i, i-probability[3]);
    }
            for(int j=0; j<probability[12]; j++)   //9
    {hRes->Fill (i, i-9);
              hMeas->Fill (i);
         hTrue->Fill (i-9);
        //response.Fill(i, i-9);
    }
            for(int j=0; j<probability[13]; j++)   //7
    {hRes->Fill (i, i-10);
              hMeas->Fill (i);
         hTrue->Fill (i-10);
        //response.Fill(i, i-10);
    }
            for(int j=0; j<probability[14]; j++)   //6
    {hRes->Fill (i, i-11);
              hMeas->Fill (i);
         hTrue->Fill (i-11);
       // response.Fill(i, i-11);
    }
            for(int j=0; j<probability[15]; j++)    //5
    {hRes->Fill (i, i-12);
              hMeas->Fill (i);
         hTrue->Fill (i-12);
       // response.Fill(i, i-12);
    }
            for(int j=0; j<probability[16]; j++)   //4
    {hRes->Fill (i, i-13);
              hMeas->Fill (i);
         hTrue->Fill (i-13);
       // response.Fill(i, i-13);
    }
            for(int j=0; j<probability[17]; j++)    //3
    {hRes->Fill (i, i-14);
              hMeas->Fill (i);
         hTrue->Fill (i-14);
        //response.Fill(i, i-14);
    }
            for(int j=0; j<probability[18]; j++)     //2
    {hRes->Fill (i, i-15);
              hMeas->Fill (i);
         hTrue->Fill (i-15);
        //response.Fill(i, i-15);
    }
}
    
    
     hRes->Scale(1/3776.0);
    hMeas->Scale(1/3776.0);
    hTrue->Scale(1/3776.0);
    
    RooUnfoldResponse response_2(hMeas, hTrue, hRes, "response", "My Response");
    RooUnfoldResponse response_4(hMeas, hTrue, hRes, "response", "My Response");
    RooUnfoldResponse response_10(hMeas, hTrue, hRes, "response", "My Response");
    RooUnfoldResponse response_16(hMeas, hTrue, hRes, "response", "My Response");
    RooUnfoldResponse response_20(hMeas, hTrue, hRes, "response", "My Response");
    RooUnfoldResponse response_25(hMeas, hTrue, hRes, "response", "My Response");
    //RooUnfoldResponse response_50(hMeas, hTrue, hRes, "response", "My Response");
    
    
    Float_t ue_estimate_tpc_const;
    _tree_event->SetBranchAddress("ue_estimate_tpc_const", &ue_estimate_tpc_const);
    
    UInt_t ntrack;
    Float_t track_pt[NTRACK_MAX];
     _tree_event->SetBranchAddress("ntrack", &ntrack);
     _tree_event->SetBranchAddress("track_pt", track_pt);
    
    Float_t track_eta[NTRACK_MAX];
    Float_t track_phi[NTRACK_MAX];
    _tree_event->SetBranchAddress("track_eta", track_eta);
    _tree_event->SetBranchAddress("track_phi", track_phi);
    
    
    //MAIN CORRELATION LOOP
    for(Long64_t ievent = 0; ievent < nentries ; ievent++){

      _tree_event->GetEntry(ievent);
      fprintf(stderr, "\r%s:%d: %llu / %llu", __FILE__, __LINE__, ievent, nentries);

      bool first_jet = true;

     // fprintf(stderr, "event selection done");
        
      //if((myTriggers & trigger_mask[0]) == 0)//MB selection  /////
      //continue;
    //accepted_events++;
        

       //fprintf(stderr,"\n centrality = %f\n",centrality_v0m);
        
      //OUTER JET LOOP   
        
        //std::cout<<"NJET"<<njet<<std::endl;
        
            
    //CAlculate background density 
        /*
         TH1D *h = new TH1D("h", "h", 100, 0, 100);
    
        for (ULong64_t n = 0; n < njet; n++) {
            if (jet_pt[n]<0.15) continue;
            if (jet_area[n]== 0) std::cout<<"ZERO AREA"<<std::endl;
            float p= jet_pt[n]/jet_area[n];
            h->Fill(p);}
    
  Double_t x, q;
  q = 0.5; // 0.5 for "median"
  h->ComputeIntegral(); // just a precaution
  h->GetQuantiles(1, &x, &q);
  //std::cout << "median = " << x << std::endl;
        density->Fill(x);
        
        float rand_eta = 0.1 * gRandom->Uniform();
        if (gRandom->Uniform()<0.5) rand_eta =-rand_eta;
        
        float rand_phi = 2.7 * gRandom->Uniform();
        if (gRandom->Uniform()<0.5) rand_phi =-rand_phi;
        float pt_total =0;
        for (ULong64_t n = 0; n < njet; n++) {
                if (jet_pt[n]<0.15) continue;
                jet_pt[n]-= x* jet_area[n];
                if(jet_eta[n]<rand_eta+ 0.4 && jet_eta[n]>rand_eta- 0.4 && jet_phi[n]<rand_phi+ 0.4 && jet_phi[n]>rand_phi- 0.4 )
                { pt_total += jet_pt[n];
                     passed_jets++;}
        }
        float f= pt_total- (x* 3.14*0.4*0.4);
        fluct->Fill(f);
        */
        float x = ue_estimate_tpc_const;
        density->Fill(x);
        
        float rand_eta = 0.1 * gRandom->Uniform();
        if (gRandom->Uniform()<0.5) rand_eta =-rand_eta;
        float rand_phi = 2.7 * gRandom->Uniform();
        if (gRandom->Uniform()<0.5) rand_phi =-rand_phi;
        
        float pt_total =0;
        for (ULong64_t itrack = 0; itrack < ntrack; itrack++) {   
            if(track_eta[itrack]<rand_eta+ 0.4 && track_eta[itrack]>rand_eta- 0.4 && track_phi[itrack]<rand_phi+ 0.4 && track_phi[itrack]>rand_phi- 0.4 )
            pt_total += track_pt[itrack];
             }
        //std::cout<<pt_total<<std::endl;
        float f= pt_total- (x* 3.14*0.4*0.4);
        fluct->Fill(f);
        if (pt_total!=0)
             fluct_non_zero->Fill(f);
        
          /*      float pt_total_non_zero =0;
           float iterations=0;
        while (pt_total_non_zero ==0){
            //iterations++;
           // if (iterations>50) std::cout<<"non zero not found"<<std::endl;
        for (ULong64_t itrack = 0; itrack < ntrack; itrack++) {   
            if(track_eta[itrack]<rand_eta+ 0.4 && track_eta[itrack]>rand_eta- 0.4 && track_phi[itrack]<rand_phi+ 0.4 && track_phi[itrack]>rand_phi- 0.4 )
            pt_total_non_zero += track_pt[itrack];
             }}
        float f_non_zero= pt_total_non_zero- (x* 3.14*0.4*0.4);
        fluct_non_zero->Fill(f_non_zero);
        
        
        */
        
        
        
        
      for (ULong64_t n = 0; n < njet; n++) {
          jet_pt[n]-= x* jet_area[n];

	if((ievent%10000==0) & n==1)
	  fprintf(stderr,"\n Jet Energy = %f, Jet pT = %f, Jet Eta = %f, Jet Phi = %f\n",jet_e[njet],jet_pt[njet],jet_eta[njet],jet_phi[njet]);
   

        if (jet_pt[n]>20){
         
         h_pt->Fill(jet_pt[n]);
         num_raw++;}
        if (jet_pt[n]>20)
            pt_comp->Fill(jet_pt[n], jet_pt_truth[n]);
        
      
          
	if (jet_pt[n] < 20) continue;
    if (jet_eta[n] <-0.5 or jet_eta[n] > 0.5) continue;

          h_mainjetpt->Fill(jet_pt[n]); 
          
          
          
          
          
    //if (jet_pt[n] < 20) continue;
    float assoc_pt=0;
    float assoc_del;

          
  //INNER JET LOOP 
  for (ULong64_t m = 0; m < njet; m++) {
   jet_pt[m]-= x* jet_area[m];

if (jet_pt[m] < 20 ) continue;
if (jet_eta[m] <-0.5 or jet_eta[m] > 0.5) continue;


	  //Observables:
      Double_t zt = jet_pt[m]/jet_pt[n];
	  Float_t DeltaPhi = TVector2::Phi_mpi_pi(jet_phi[n]- jet_phi[m]);
	  Float_t DeltaEta = jet_eta[n]- jet_eta[m];

      
      
      
	  //if ((TMath::Abs(DeltaPhi) < 0.005) && (TMath::Abs(DeltaEta) < 0.005)) continue; //Match Mixing Cut
      

              

          
        Float_t  del = TMath::Abs(jet_phi[n]-jet_phi[m]);
           
          
            if (2*3.14/3 <= del && del <= 3.14){
                if (jet_pt[m]>assoc_pt){
                    assoc_pt= jet_pt[m];
                    assoc_del= del;
                }
            } 

   




    

         Float_t k_t = jet_pt[n] * TMath::Sin(assoc_del);
         h_delta_phi->Fill(assoc_del);
      
          kt_pt->Fill(jet_pt[n],assoc_pt);


          if (20<jet_pt[n] && jet_pt[n]<40){
              float val= (k_t*0.0152319) +0.868044;
              kt_20_40->Fill(k_t);
            
               }
              
        if (40<jet_pt[n] && jet_pt[n]<60) {
               int val =(k_t*0.01) +0.9;
               kt_40_60->Fill(k_t );
        
        
        }
           



          
          
   

      }//Jet Loop


    } //outer jet loop
        
        
        for (int y=0; y<1000; y++){
        //std::cout<<njet_truth<<std::endl;
          for (ULong64_t n = 0; n < njet; n++) {
              jet_pt[n]-= x* jet_area[n];
       jet_pt[n]= prob_pt(jet_pt[n]);

          
      
        if (jet_pt[n]>20)
         h_pt_truth->Fill(jet_pt[n]);
          num_truth++;
       
    //if (jet_pt_truth[n]!=0) std::cout<<jet_pt_truth[n]<<std::endl;
	if (jet_pt[n] < 20) continue;
    if (jet_eta[n] <-0.5 or jet_eta[n] > 0.5) continue;

     
    float assoc_pt_truth=0;
    float assoc_del_truth;
          
  //INNER JET LOOP 
  for (ULong64_t m = 0; m < njet; m++) {
      jet_pt[m]-= x* jet_area[m];
jet_pt[m]= prob_pt(jet_pt[m]);

if (jet_pt[m]<20) continue;
if (jet_eta[m] <-0.5 or jet_eta[m] > 0.5) continue;
     
       
	  //Observables:
      
	  Float_t DeltaPhi_truth = TVector2::Phi_mpi_pi(jet_phi[n]- jet_phi[m]);
	  Float_t DeltaEta_truth = jet_eta[n]- jet_eta[m];
      
      
      
	  if ((TMath::Abs(DeltaPhi_truth) < 0.005) && (TMath::Abs(DeltaEta_truth) < 0.005)) continue; //Match Mixing Cut
      

          


        Float_t  del_truth = TMath::Abs(jet_phi[n]-jet_phi[m]);
           
          
            if (2*3.14/3 <= del_truth && del_truth <= 3.14){
                if (jet_pt[m]>assoc_pt_truth){
                    assoc_pt_truth= jet_pt[m];
                    assoc_del_truth= del_truth;
                }
            } 
           
  


  
           Float_t k_t_truth = jet_pt[n] * TMath::Sin(assoc_del_truth);
             h_delta_phi_truth->Fill(assoc_del_truth);
          kt_pt_truth->Fill(jet_pt[n],assoc_pt_truth);
          

          if (20<jet_pt[n] && jet_pt[n]<40){
              kt_20_40_truth->Fill(k_t_truth);

            
               }
              
        if (40<jet_pt[n] && jet_pt[n]<60) {
               kt_40_60_truth->Fill(k_t_truth);

        
        
        }
              //kt_pt_truth->Fill(jet_pt_truth[n],k_t_truth);

     
   

      }//Jet Loop


    } //outer jet loop
        }  
        
        
        
        
    }//Event Loop
  /*
   RooUnfoldSvd unfold_2 (&response_2, h_mainjetpt, 2);
    TH1D* hReco_2= (TH1D*) unfold_2.Hreco();
       RooUnfoldSvd unfold_4 (&response_4, h_mainjetpt, 4);
    TH1D* hReco_4= (TH1D*) unfold_4.Hreco();
       RooUnfoldSvd unfold_10 (&response_10, h_mainjetpt, 10);
    TH1D* hReco_10= (TH1D*) unfold_10.Hreco();
       RooUnfoldSvd unfold_16 (&response_16, h_mainjetpt, 16);
    TH1D* hReco_16= (TH1D*) unfold_16.Hreco();
       RooUnfoldSvd unfold_20 (&response_20, h_mainjetpt, 20);
    TH1D* hReco_20= (TH1D*) unfold_20.Hreco();
       RooUnfoldSvd unfold_25 (&response_25, h_mainjetpt, 25);
    TH1D* hReco_25= (TH1D*) unfold_25.Hreco();
     //  RooUnfoldSvd unfold_50 (&response_50, h_mainjetpt, 50);
    //TH1D* hReco_50= (TH1D*) unfold_50.Hreco();
        
       */
            
    //Automatic Output Naming
    size_t lastindex = std::string(root_file).find_last_of(".");
    std::string rawname = std::string(root_file).substr(0, lastindex);
    std::cout<<"Creating new file"<<std::endl;

  TFile* fout;


  fout = new TFile("pPbDATAJune9_regparameter_pt.root","RECREATE");


    

    
        /*
    	  for (int ipt = 0; ipt < nptbins; ipt++){
	    
	      for(int izt = 0; izt<nztbins ; izt++){
             

 

       TH1D * projh2X = Corr->ProjectionX();
      TH1D * projh2Y = Corr->ProjectionY();
                       projh2Y .Sumw2();
          projh2X .Sumw2();
                
                projh2Y .Write();
                projh2X .Write();
                
        }}
            */
    
        

    std::cout<<std::endl;

      
    
    
    
 

        //std::cout<<nptbins<<" "<<nztbins<<"Checking"<<std::endl;
        
        // kt_40_60->Write();
        //kt_20_40->Write();
        
        
        
        
        
        Double_t scale_truth = kt_20_40->GetXaxis()->GetBinWidth(1)*(kt_20_40->GetEntries()); 
        Double_t scalenxt_truth = kt_40_60->GetXaxis()->GetBinWidth(1)*(kt_40_60->GetEntries());
        kt_20_40->Scale(1.0/scale_truth);
        kt_40_60->Scale(1.0/scalenxt_truth);
       
        
                TF1 *func = new TF1("func","0.01*x + 0.9",0,10);
 for (int i=0; i<kt_40_60->GetNbinsX(); i++){
Double_t BinContent = kt_40_60->GetBinContent(i);
Double_t BinCenter =kt_40_60->GetBinCenter(i);
Double_t BinError = kt_40_60->GetBinError(i);
Double_t po = func->Eval(BinCenter);
kt_40_60->SetBinContent(i,BinContent*po);
kt_40_60->SetBinError(i,BinError*po);
}
        
                TF1 *func2 = new TF1("func2","0.0152319*x + 0.868044",0,10);
 for (int i=0; i<kt_20_40->GetNbinsX(); i++){
Double_t BinContent = kt_20_40->GetBinContent(i);
Double_t BinCenter =kt_20_40->GetBinCenter(i);
Double_t BinError = kt_20_40->GetBinError(i);
Double_t po = func2->Eval(BinCenter);
kt_20_40->SetBinContent(i,BinContent*po);
kt_20_40->SetBinError(i,BinError*po);
}
        
        
        
        
        
            Double_t scale = kt_20_40_truth->GetEntries(); 
        Double_t scalenxt = kt_40_60_truth->GetEntries();
        kt_20_40_truth->Scale(1.0/scale);
        kt_40_60_truth->Scale(1.0/scalenxt);
    
    
    for(int i=0;i<5;i++){
Double_t dKty = kt_20_40_truth->GetBinWidth(i+1);
Double_t BinContent = kt_20_40_truth->GetBinContent(i+1);
Double_t BinError = kt_20_40_truth->GetBinError(i+1);
kt_20_40_truth->SetBinContent(i+1,BinContent/dKty);
kt_20_40_truth->SetBinError(i+1,BinError/dKty);
}
        
    for(int i=0;i<7;i++){
Double_t dKty = kt_40_60_truth->GetBinWidth(i+1);
Double_t BinContent = kt_40_60_truth->GetBinContent(i+1);
Double_t BinError = kt_40_60_truth->GetBinError(i+1);
kt_40_60_truth->SetBinContent(i+1,BinContent/dKty);
kt_40_60_truth->SetBinError(i+1,BinError/dKty);
}
    
    

    
    
       
        
                TF1 *funct = new TF1("funct","0.01*x + 0.9",0,10);
 for (int i=0; i<kt_40_60_truth->GetNbinsX(); i++){
Double_t BinContent = kt_40_60_truth->GetBinContent(i);
Double_t BinCenter =kt_40_60_truth->GetBinCenter(i);
Double_t BinError = kt_40_60_truth->GetBinError(i);
Double_t po = funct->Eval(BinCenter);
kt_40_60_truth->SetBinContent(i,BinContent*po);
kt_40_60_truth->SetBinError(i,BinError*po);
}
        
                TF1 *func2t = new TF1("func2t","0.0152319*x + 0.868044",0,10);
 for (int i=0; i<kt_20_40_truth->GetNbinsX(); i++){
Double_t BinContent = kt_20_40_truth->GetBinContent(i);
Double_t BinCenter =kt_20_40_truth->GetBinCenter(i);
Double_t BinError = kt_20_40_truth->GetBinError(i);
Double_t po = func2t->Eval(BinCenter);
kt_20_40_truth->SetBinContent(i,BinContent*po);
kt_20_40_truth->SetBinError(i,BinError*po);
}
        
            Double_t scale_delta_phi = h_delta_phi->GetXaxis()->GetBinWidth(1)*(h_delta_phi->GetEntries()); 
        Double_t scale_delta_phi_truth = h_delta_phi_truth->GetXaxis()->GetBinWidth(1)*(h_delta_phi_truth->GetEntries());
        h_delta_phi->Scale(1.0/scale_delta_phi);
        h_delta_phi_truth->Scale(1.0/scale_delta_phi_truth);
        
        
        
        kt_pt->SetOption("COLZ");
        //kt_pt->Write();
        kt_pt_truth->SetOption("COLZ");
        //kt_pt_truth->Write();
         
           kt_40_60_truth->Write();
        kt_20_40_truth->Write();
    
       
   Double_t kt_20[5] = {0.0758, 0.079, 0.0316, 0.0139, 0.0021};
   Double_t kt_20_err[5] = {0.011, 0.0079, 0.0036, 0.0024, 0.0009};
    
    Double_t kt_40[7] = {0.0489, 0.0473, 0.0345, 0.0172, 0.0109, 0.0099, 0.0044};
   Double_t kt_40_err[7] = {0.0153, 0.0109, 0.0066, 0.0048, 0.0039, 0.0037, 0.0025};
    
     for (int i=1; i<=kt_20_40_truth->GetNbinsX(); i++){
Double_t BinContent = kt_20_40_truth->GetBinContent(i);
Double_t BinError = kt_20_40_truth->GetBinError(i);
Double_t po = kt_20[i-1];
         Double_t qo = kt_20_err[i-1];
         Double_t err = ((BinError/BinContent)+ (qo/ po))* (BinContent/po);
kt_20_40_truth->SetBinContent(i,BinContent/po);
         
kt_20_40_truth->SetBinError(i,err);
}
         for (int i=1; i<=kt_40_60_truth->GetNbinsX(); i++){
Double_t BinContent = kt_40_60_truth->GetBinContent(i);
Double_t BinError = kt_40_60_truth->GetBinError(i);
Double_t po = kt_40[i-1];
         Double_t qo = kt_40_err[i-1];
             Double_t err = ((BinError/BinContent)+ (qo/ po))* (BinContent/po);
kt_40_60_truth->SetBinContent(i,BinContent/po);
kt_40_60_truth->SetBinError(i,err);
}
        
        
        
        
        //kt_40_60->Write();
        //kt_20_40->Write();
        

       
        //std::cout<<kt_40_60_truth->GetEntries()<<" "<<kt_40_60->GetEntries()<<std::endl;
        
  /*
       
        
        corr_40_60-> Fill(kt_40_60_truth->GetEntries()/kt_40_60->GetEntries());
         corr_20_40-> Fill(kt_20_40_truth->GetEntries()/kt_20_40->GetEntries());
        
            corr_40_60->Write();
        corr_20_40->Write();
        
        
        kt_pt->SetOption("COLZ");
        kt_pt->Write();
        kt_pt_truth->SetOption("COLZ");
        kt_pt_truth->Write();
        
        pt_comp->SetOption("COLZ");
        pt_comp->Write();
        
        hs->Add(h_pt );
        hs->Add(h_pt_truth );
        hs->Write();
        

        
        TH1D *locHist_ratio = (TH1D*)kt_20_40_truth->Clone("locHist_ratio");
        locHist_ratio->Sumw2();
        locHist_ratio->Divide(kt_20_40 );
        locHist_ratio->Write();
        
        TH1D *locHist_ratio_40 = (TH1D*)kt_40_60_truth->Clone("locHist_ratio_40");
        locHist_ratio_40->Sumw2();
        locHist_ratio_40->Divide(kt_40_60 );
        locHist_ratio_40->Write();
        
        
        */
        
        
            TCanvas *canvas_kt2040 = new TCanvas("canvas","multipads",900,700);
   // kt_20_40->SetMarkerColor(633);
    kt_20_40_truth->SetMarkerColor(602);
    //kt_20_40->SetLineColor(633);
    kt_20_40_truth->SetLineColor(602);
    //kt_20_40->Draw("ep");
    kt_20_40_truth->Draw("ep");  
    //TLegend *leg204 = new TLegend(0.7,0.3,0.85,0.4,"");
//leg204->AddEntry(kt_20_40_truth ,"#bf{kT with fluctuations}","lp");
//leg204->AddEntry(kt_20_40 ,"#bf{kT}","lp");
    //leg204->Draw();
     canvas_kt2040->Update();
     canvas_kt2040->SaveAs("kt2040_RATIO_bg.pdf");
       
   TCanvas *canvas_kt4060 = new TCanvas("canvas","multipads",900,700);
    //kt_40_60->SetMarkerColor(633);
    kt_40_60_truth->SetMarkerColor(602);
    //kt_40_60->SetLineColor(633);
    kt_40_60_truth->SetLineColor(602);
    //kt_40_60->Draw("ep");
    kt_40_60_truth->Draw("ep");  
    //TLegend *leg40 = new TLegend(0.7,0.3,0.85,0.4,"");
//leg40->AddEntry(kt_40_60_truth ,"#bf{kT with fluctuations}","lp");
//leg40->AddEntry(kt_40_60 ,"#bf{kT}","lp");
  //  leg40->Draw();
     canvas_kt4060->Update();
     canvas_kt4060->SaveAs("kt4060_RATIO_bg.pdf");
        
      
        
        
        


    // h_mainjetpt->Write();
 //density->Write();
  //fluct->Write();
    //fluct_non_zero->Write();
  // hReco->Write();
  hRes->SetOption("COLZ");
    hRes->Scale(1/3776.0);
  //  hRes->Write();
    /*
    
    TCanvas *canvas_phi = new TCanvas("canvas","multipads",900,700);
    h_delta_phi->SetMarkerColor(633);
    h_delta_phi_truth->SetMarkerColor(602);
    h_delta_phi->SetLineColor(633);
    h_delta_phi_truth->SetLineColor(602);
    h_delta_phi->Draw("ep");
    h_delta_phi_truth->Draw("same,ep");  
    TLegend *leg = new TLegend(0.7,0.3,0.85,0.4,"");
leg->AddEntry(h_delta_phi_truth,"#bf{delta phi without bg subtraction}","lp");
leg->AddEntry(h_delta_phi,"#bf{delta phi}","lp");
    leg->Draw();
     canvas_phi->Update();
     //canvas_phi->SaveAs("delta_phi_bg.pdf");
    
    */
    
    /*
     TCanvas *canvas_2 = new TCanvas("canvas","multipads",900,700);
    hReco_2->SetMarkerColor(633);
    h_mainjetpt_toy->SetMarkerColor(602);
        hReco_2->SetLineColor(633);
    h_mainjetpt_toy->SetLineColor(602);
    h_mainjetpt_toy->Draw("ep");
   // hReco_2->Draw("same,ep");  
    TLegend *leg = new TLegend(0.7,0.3,0.85,0.4,"");
//leg->AddEntry(hReco_2,"#bf{pT after unfolding}","lp");
leg->AddEntry(h_mainjetpt_toy,"#bf{pT}","lp");
    leg->Draw();
     canvas_2->Update();
     canvas_2->SaveAs("pt_dist_zeroprob.pdf");
     */
    
    /*
         TCanvas *canvas_4 = new TCanvas("canvas","multipads",900,700);
    hReco_4->SetMarkerColor(633);
    h_mainjetpt->SetMarkerColor(602);
        hReco_4->SetLineColor(633);
    h_mainjetpt->SetLineColor(602);
    h_mainjetpt->Draw("ep");
    hReco_4->Draw("same,ep");  
    TLegend *leg4 = new TLegend(0.7,0.3,0.85,0.4,"");
leg4->AddEntry(hReco_4,"#bf{pT after unfolding}","lp");
leg4->AddEntry(h_mainjetpt,"#bf{pT}","lp");
    leg4->Draw();
     canvas_4->Update();
     canvas_4->SaveAs("unfold4.pdf");
    
         TCanvas *canvas_10 = new TCanvas("canvas","multipads",900,700);
    hReco_10->SetMarkerColor(633);
    h_mainjetpt->SetMarkerColor(602);
        hReco_10->SetLineColor(633);
    h_mainjetpt->SetLineColor(602);
    h_mainjetpt->Draw("ep");
    hReco_10->Draw("same,ep");  
    TLegend *leg10 = new TLegend(0.7,0.3,0.85,0.4,"");
leg10->AddEntry(hReco_10,"#bf{pT after unfolding}","lp");
leg10->AddEntry(h_mainjetpt,"#bf{pT}","lp");
    leg10->Draw();
     canvas_10->Update();
     canvas_10->SaveAs("unfold10.pdf");
    
         TCanvas *canvas_16 = new TCanvas("canvas","multipads",900,700);
    hReco_16->SetMarkerColor(633);
    h_mainjetpt->SetMarkerColor(602);
        hReco_16->SetLineColor(633);
    h_mainjetpt->SetLineColor(602);
    h_mainjetpt->Draw("ep");
    hReco_16->Draw("same,ep");  
    TLegend *leg16 = new TLegend(0.7,0.3,0.85,0.4,"");
leg16->AddEntry(hReco_16,"#bf{pT after unfolding}","lp");
leg16->AddEntry(h_mainjetpt,"#bf{pT}","lp");
    leg16->Draw();
     canvas_16->Update();
     canvas_16->SaveAs("unfold16.pdf");
    
         TCanvas *canvas_20 = new TCanvas("canvas","multipads",900,700);
    hReco_20->SetMarkerColor(633);
    h_mainjetpt->SetMarkerColor(602);
        hReco_20->SetLineColor(633);
    h_mainjetpt->SetLineColor(602);
    h_mainjetpt->Draw("ep");
    hReco_20->Draw("same,ep");  
    TLegend *leg20 = new TLegend(0.7,0.3,0.85,0.4,"");
leg20->AddEntry(hReco_20,"#bf{pT after unfolding}","lp");
leg20->AddEntry(h_mainjetpt,"#bf{pT}","lp");
    leg20->Draw();
     canvas_20->Update();
     canvas_20->SaveAs("unfold20.pdf");
    
         TCanvas *canvas_30 = new TCanvas("canvas","multipads",900,700);
    hReco_25->SetMarkerColor(633);
    h_mainjetpt->SetMarkerColor(602);
        hReco_25->SetLineColor(633);
    h_mainjetpt->SetLineColor(602);
    h_mainjetpt->Draw("ep");
    hReco_25->Draw("same,ep");  
    TLegend *leg30 = new TLegend(0.7,0.3,0.85,0.4,"");
leg30->AddEntry(hReco_25,"#bf{pT after unfolding}","lp");
leg30->AddEntry(h_mainjetpt,"#bf{pT}","lp");
    leg30->Draw();
     canvas_30->Update();
     canvas_30->SaveAs("unfold25.pdf");
    */
    /*
         TCanvas *canvas_50 = new TCanvas("canvas","multipads",900,700);
    hReco_50->SetMarkerColor(633);
    h_mainjetpt->SetMarkerColor(602);
        hReco_50->SetLineColor(633);
    h_mainjetpt->SetLineColor(602);
    h_mainjetpt->Draw("ep");
    hReco_50->Draw("same,ep");  
    TLegend *leg50 = new TLegend(0.7,0.3,0.85,0.4,"");
leg50->AddEntry(hReco_50,"#bf{pT after unfolding}","lp");
leg50->AddEntry(h_mainjetpt,"#bf{pT}","lp");
    leg50->Draw();
     canvas_50->Update();
     canvas_50->SaveAs("unfold50.pdf");
    
    */

  //Seperate zt loops for easier file reading
  fout->Close();
  file->Close();
  std::cout << " ending " << std::endl;
  return EXIT_SUCCESS;
}
