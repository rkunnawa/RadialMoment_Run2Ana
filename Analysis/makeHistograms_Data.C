//! Raghav Kunnawalkam Elayavalli
//! Jan 12th 2017
//! Rutgers 
//! for questions or comments: raghav.k.e at CERN dot CH

//! Macro to read Data hiforest and make simple histograms 

#include <iostream>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <TH1F.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TStyle.h>
#include <TStopwatch.h>
#include <TRandom3.h>
#include <TChain.h>
#include <TProfile.h>
#include <TStopwatch.h>
#include <TEventList.h>
#include <TSystem.h>
#include <TCut.h>
#include <cstdlib>
#include <cmath>
#include "TLegend.h"
#include "TLatex.h"
#include "TMath.h"
#include "TLine.h"

double ptbins_ana[] = {20.0, 30.0, 50.0, 100.0, 150.0, 200.0, 300.0, 600.0};
const int nbins_pt_ana = sizeof(ptbins_ana)/sizeof(double)-1;

double etabins_ana[] = {-3.0, -1.0, 1.0, 3.0};
const int nbins_eta_ana = sizeof(etabins_ana)/sizeof(double)-1; 

static const int trigValue = 4;
static const char trigName [trigValue][256] = {"HLT40","HLT60","HLT80","Combined"};

// divide by bin width
void divideBinWidth(TH1 *h){
  h->Sumw2();
  for (int i=0;i<=h->GetNbinsX();i++){
    Float_t val = h->GetBinContent(i);
    Float_t valErr = h->GetBinError(i);
    val/=h->GetBinWidth(i);
    valErr/=h->GetBinWidth(i);
    h->SetBinContent(i,val);
    h->SetBinError(i,valErr);
  }//binsX loop 
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
}

float deltaphi(float phi1, float phi2)
{
  float pi=TMath::Pi();
 
  float dphi = TMath::Abs(phi1 - phi2);
  if(dphi > pi)dphi -= 2*pi;

  return TMath::Abs(dphi);
}

double eventPrescl(int *trgDec, int *treePrescl, double triggerPt){
  float weight = 0.0;
  if(trgDec[0] && triggerPt>=40.  && triggerPt<60. )
    weight = treePrescl[0];
  if(trgDec[1] && triggerPt>=60.  && triggerPt<80. )
    weight = treePrescl[1];
  if(trgDec[2] && triggerPt>=80. )
    weight = treePrescl[2];
  return weight;
}


using namespace std;

void makeHistograms_Data(int startfile = 0,
			 int endfile = 1,
			 int radius = 4,
			 std::string algo="",
			 std::string jetType= "PF",
			 std::string kfoutname="test.root"){
  
  TStopwatch timer;
  timer.Start();
  
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  gStyle->SetOptStat(0);

  bool printDebug = false;
  if(printDebug)cout<<"radius = "<<radius<<endl;
  
  TDatime date;

  std::string infile_Forest;
  infile_Forest = "list.txt";
  std::ifstream instr_Forest(infile_Forest.c_str(),std::ifstream::in);
  std::string filename_Forest;
  
  if(printDebug)cout<<"reading from "<<startfile<<" to "<<endfile<<endl;
  
  for(int ifile = 0;ifile<startfile;ifile++){
    instr_Forest>>filename_Forest;
  }

  const int N = 7; 

  TChain * jetTree[N];

  string dir[N];
  dir[0] = "hltanalysis";
  dir[1] = "skimanalysis";
  dir[2] = Form("ak%s%d%sJetAnalyzer", algo.c_str(), radius, jetType.c_str());
  dir[3] = "hiEvtAnalyzer";
  dir[4] = "hltobject";
  dir[5] = "hltobject";
  dir[6] = "hltobject";

  string trees[N] = {
    "HltTree",
    "HltTree",
    "t",
    "HiTree",
    "HLT_PAAK4PFJet40_Eta5p1_v",
    "HLT_PAAK4PFJet60_Eta5p1_v",
    "HLT_PAAK4PFJet80_Eta5p1_v"
  };

  for(int t = 0;t<N;t++){
    jetTree[t] = new TChain(string(dir[t]+"/"+trees[t]).data());
  }//tree loop ends

  for(int ifile = startfile; ifile<endfile; ++ifile){

    instr_Forest>>filename_Forest;

    for(int i = 0; i<N; ++i){
      jetTree[i]->Add(filename_Forest.c_str());    
      if(printDebug)cout << "Tree loaded  " << string(dir[i]+"/"+trees[i]).data() << endl;
      if(printDebug)cout << "Entries : " << jetTree[i]->GetEntries() << endl;
    }

  }

  for(int i = 0; i<N; ++i)  
    if(i!=2)
      jetTree[2]->AddFriend(jetTree[i]);

  int jet40_F;
  int jet60_F;
  int jet80_F;
  int jet40_p_F;
  int jet60_p_F;
  int jet80_p_F;
  int jet40_l1_F;
  int jet60_l1_F;
  int jet80_l1_F;
  int jet40_p_l1_F;
  int jet60_p_l1_F;
  int jet80_p_l1_F;
  jetTree[0]->SetBranchAddress("HLT_PAAK4PFJet40_Eta5p1_v3",&jet40_F);
  jetTree[0]->SetBranchAddress("HLT_PAAK4PFJet40_Eta5p1_v3_Prescl",&jet40_p_F);
  jetTree[0]->SetBranchAddress("HLT_PAAK4PFJet60_Eta5p1_v4",&jet60_F);
  jetTree[0]->SetBranchAddress("HLT_PAAK4PFJet60_Eta5p1_v4_Prescl",&jet60_p_F);
  jetTree[0]->SetBranchAddress("HLT_PAAK4PFJet80_Eta5p1_v3",&jet80_F);
  jetTree[0]->SetBranchAddress("HLT_PAAK4PFJet80_Eta5p1_v3_Prescl",&jet80_p_F);
  jetTree[0]->SetBranchAddress("L1_SingleJet16_BptxAND_Final",&jet40_l1_F);
  jetTree[0]->SetBranchAddress("L1_SingleJet16_BptxAND_Prescl",&jet40_p_l1_F);
  jetTree[0]->SetBranchAddress("L1_SingleJet24_BptxAND_Final",&jet60_l1_F);
  jetTree[0]->SetBranchAddress("L1_SingleJet24_BptxAND_Prescl",&jet60_p_l1_F);
  jetTree[0]->SetBranchAddress("L1_SingleJet36_BptxAND_Final",&jet80_l1_F);
  jetTree[0]->SetBranchAddress("L1_SingleJet36_BptxAND_Prescl",&jet80_p_l1_F);
  
  int pprimaryvertexFilter_F;
  int pVertexFilterCutGplus_F;
  int pHBHENoiseFilter_F;
  jetTree[1]->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&pHBHENoiseFilter_F);
  jetTree[1]->SetBranchAddress("pPAprimaryVertexFilter",&pprimaryvertexFilter_F);
  jetTree[1]->SetBranchAddress("pVertexFilterCutGplus",&pVertexFilterCutGplus_F);
  
  int nref_F;
  float pt_F[1000];
  float rawpt_F[1000];
  float eta_F[1000];
  float phi_F[1000];
  float m_F[1000];
  float chMax_F[1000];
  int chN_F[1000];
  int neN_F[1000];
  float trkMax_F[1000];
  float chSum_F[1000];
  float phSum_F[1000];
  float neSum_F[1000];
  float trkSum_F[1000];
  float phMax_F[1000];
  float neMax_F[1000];
  float chHardMax_F[1000];
  float trkHardMax_F[1000];
  float chHardSum_F[1000];
  float phHardSum_F[1000];
  float trkHardSum_F[1000];
  float phHardMax_F[1000];
  float eMax_F[1000];
  float muMax_F[1000];
  float eSum_F[1000];
  float muSum_F[1000];
  float jtpu_F[1000];
  jetTree[2]->SetBranchAddress("nref",&nref_F);
  jetTree[2]->SetBranchAddress("jtpt",&pt_F);
  jetTree[2]->SetBranchAddress("jteta",&eta_F);
  jetTree[2]->SetBranchAddress("jtphi",&phi_F);
  jetTree[2]->SetBranchAddress("rawpt",&rawpt_F);
  jetTree[2]->SetBranchAddress("jtpu",&jtpu_F);
  jetTree[2]->SetBranchAddress("jtm",&m_F);
  jetTree[2]->SetBranchAddress("chargedMax",&chMax_F);
  jetTree[2]->SetBranchAddress("chargedSum",&chSum_F);
  jetTree[2]->SetBranchAddress("chargedN",&chN_F);
  jetTree[2]->SetBranchAddress("neutralN",&neN_F);
  //jetTree[2]->SetBranchAddress("chargedHardMax",&chMax_F);
  jetTree[2]->SetBranchAddress("chargedHardSum",&chSum_F);
  jetTree[2]->SetBranchAddress("trackMax",&trkMax_F);
  jetTree[2]->SetBranchAddress("trackSum",&trkSum_F);
  //jetTree[2]->SetBranchAddress("trackHardMax",&trkMax_F);
  jetTree[2]->SetBranchAddress("trackHardSum",&trkSum_F);
  jetTree[2]->SetBranchAddress("photonMax",&phMax_F);
  jetTree[2]->SetBranchAddress("photonSum",&phSum_F);
  //jetTree[2]->SetBranchAddress("photonHardMax",&phMax_F);
  jetTree[2]->SetBranchAddress("photonHardSum",&phSum_F);
  jetTree[2]->SetBranchAddress("neutralMax",&neMax_F);
  jetTree[2]->SetBranchAddress("neutralSum",&neSum_F);
  jetTree[2]->SetBranchAddress("eSum",&eSum_F);
  jetTree[2]->SetBranchAddress("eMax",&eMax_F);
  jetTree[2]->SetBranchAddress("muSum",&muSum_F);
  jetTree[2]->SetBranchAddress("muMax",&muMax_F);
  //! pf jet variables
  float jtPfCHF_F[1000];
  float jtPfNHF_F[1000];
  float jtPfCEF_F[1000];
  float jtPfNEF_F[1000];
  float jtPfMUF_F[1000];
  float jtPfCHM_F[1000];
  float jtPfNHM_F[1000];
  float jtPfCEM_F[1000];
  float jtPfNEM_F[1000];
  float jtPfMUM_F[1000];
  jetTree[2]->SetBranchAddress("jtPfCHF",&jtPfCHF_F);
  jetTree[2]->SetBranchAddress("jtPfNHF",&jtPfNHF_F);
  jetTree[2]->SetBranchAddress("jtPfCEF",&jtPfCEF_F);
  jetTree[2]->SetBranchAddress("jtPfNEF",&jtPfNEF_F);
  jetTree[2]->SetBranchAddress("jtPfMUF",&jtPfMUF_F);
  jetTree[2]->SetBranchAddress("jtPfCHM",&jtPfCHM_F);
  jetTree[2]->SetBranchAddress("jtPfNHM",&jtPfNHM_F);
  jetTree[2]->SetBranchAddress("jtPfCEM",&jtPfCEM_F);
  jetTree[2]->SetBranchAddress("jtPfNEM",&jtPfNEM_F);
  jetTree[2]->SetBranchAddress("jtPfMUM",&jtPfMUM_F);

  Float_t vz_F;
  ULong64_t evt_F;
  UInt_t run_F;
  UInt_t lumi_F;
  int hiBin_F;
  int hiNpix_F;
  int hiNtracks_F;
  float hiHF_F;
  float hiHFplus_F;
  float hiHFminus_F;
  float hiZDC_F;
  float hiZDCplus_F;
  float hiZDCminus_F;
  float hiHFhit_F;
  float hiHFhitplus_F;
  float hiHFhitminus_F;
  jetTree[3]->SetBranchAddress("evt",&evt_F);
  jetTree[3]->SetBranchAddress("run",&run_F);
  jetTree[3]->SetBranchAddress("lumi",&lumi_F);
  jetTree[3]->SetBranchAddress("vz",&vz_F);
  jetTree[3]->SetBranchAddress("hiBin",&hiBin_F);
  jetTree[3]->SetBranchAddress("hiHF",&hiHF_F);
  jetTree[3]->SetBranchAddress("hiHFplus",&hiHFplus_F);
  jetTree[3]->SetBranchAddress("hiHFminus",&hiHFminus_F);
  jetTree[3]->SetBranchAddress("hiZDC",&hiZDC_F);
  jetTree[3]->SetBranchAddress("hiZDCplus",&hiZDCplus_F);
  jetTree[3]->SetBranchAddress("hiZDCminus",&hiZDCminus_F);
  jetTree[3]->SetBranchAddress("hiHFhit",&hiHFhit_F);
  jetTree[3]->SetBranchAddress("hiHFhitPlus",&hiHFhitplus_F);
  jetTree[3]->SetBranchAddress("hiHFhitMinus",&hiHFhitminus_F);
  jetTree[3]->SetBranchAddress("hiNpix",&hiNpix_F);
  jetTree[3]->SetBranchAddress("hiNtracks",&hiNtracks_F);

  //! vectors for hltobject tree
  std::vector<double>  *trgObjpt_40, *trgObjpt_60, *trgObjpt_80;
  jetTree[4]->SetBranchAddress("pt",&trgObjpt_40);
  jetTree[5]->SetBranchAddress("pt",&trgObjpt_60);  
  jetTree[6]->SetBranchAddress("pt",&trgObjpt_80);  
  
  TFile *fout = new TFile(kfoutname.c_str(),"RECREATE");
  fout->cd();
  //! Declare the output histograms to compare data with MC.
  //! Event histograms 
  TH1F * hVz = new TH1F("hVz","",100, -20, 20);
  TH1F * hhiBin = new TH1F("hhiBin","",200, 0, 200);
  TH1F * hNpix = new TH1F("hNpix","", 100, 0, 1200);
  TH1F * hHF = new TH1F("hHF","",500, 0, 500);
  TH1F * hHFplus = new TH1F("hHFplus","",500, 0, 500);
  TH1F * hHFminus = new TH1F("hHFminus","",500, 0, 500);
  TH1F * hZDC = new TH1F("hZDC","",500, 0, 500);
  TH1F * hZDCplus = new TH1F("hZDCplus","",500, 0, 500);
  TH1F * hZDCminus = new TH1F("hZDCminus","",500, 0, 500);
  TH1F * hHFhit = new TH1F("hHFhit","",200, 0, 8000);
  TH1F * hHFhitplus = new TH1F("hHFhitplus","",200, 0, 8000);
  TH1F * hHFhitminus = new TH1F("hHFhitminus","",200, 0, 8000);
  TH1F * hNtracks = new TH1F("hNtracks","", 200, 0, 200);
  //! Jet histograms
  TH1F * hjtpt = new TH1F("hjtpt","", 500, 0, 1000);
  TH1F * hjteta = new TH1F("hjteta","", 100, -5, 5);
  TH1F * hjtphi = new TH1F("hjtphi","", 100, -4, 4);
  TH1F * hjtm = new TH1F("hjtm","", 50, 0, 50);
  TH1F * hCHF = new TH1F("hCHF","", 50, 0, 1.5);
  TH1F * hNHF = new TH1F("hNHF","", 50, 0, 1.5);
  TH1F * hCEF = new TH1F("hCEF","", 50, 0, 1.5);
  TH1F * hNEF = new TH1F("hNEF","", 50, 0, 1.5);
  TH1F * hMUF = new TH1F("hMUF","", 50, 0, 1.5);
  TH1F * hCHM = new TH1F("hCHM","", 70, 0, 70);
  TH1F * hNHM = new TH1F("hNHM","", 70, 0, 70);
  TH1F * hCEM = new TH1F("hCEM","", 70, 0, 70);
  TH1F * hNEM = new TH1F("hNEM","", 70, 0, 70);
  TH1F * hMUM = new TH1F("hMUM","", 70, 0, 70);
  TH1F * hAj = new TH1F("hAj","",50, 0, 1);
  //! 2D histograms with Aj vs different centrality variables
  TH2F * hAjvshiBin = new TH2F("hAjvshiBin","", 200, 0, 200, 50, 0, 1);
  TH2F * hAjvsNpix = new TH2F("hAjvsNpix","", 100, 0, 1200, 50, 0, 1);
  TH2F * hAjvsHF = new TH2F("hAjvsHF","", 500, 0, 500, 50, 0, 1);
  TH2F * hAjvsHFplus = new TH2F("hAjvsHFplus","", 500, 0, 500, 50, 0, 1);
  TH2F * hAjvsHFminus = new TH2F("hAjvsHFminus","", 500, 0, 500, 50, 0, 1);
  TH2F * hAjvsZDC = new TH2F("hAjvsZDC","", 500, 0, 500, 50, 0, 1);
  TH2F * hAjvsZDCplus = new TH2F("hAjvsZDCplus","", 500, 0, 500, 50, 0, 1);
  TH2F * hAjvsZDCminus = new TH2F("hAjvsZDCminus","", 500, 0, 500, 50, 0, 1);
  TH2F * hAjvsHFhit = new TH2F("hAjvsHFhit","", 200, 0, 8000, 50, 0, 1);
  TH2F * hAjvsHFhitplus = new TH2F("hAjvsHFhitplus","", 200, 0, 8000, 50, 0, 1);
  TH2F * hAjvsHFhitminus = new TH2F("hAjvsHFhitminus","", 200, 0, 8000, 50, 0, 1);
  TH2F * hAjvsNtracks = new TH2F("hAjvsNtracks","", 200, 0, 200, 50, 0, 1);
  
  //! now start the event loop for each file.  
  if(printDebug) cout<<"Running through all the events now"<<endl;
  Long64_t nentries = jetTree[0]->GetEntries();
  if(printDebug) nentries = 100;
  
  for(int nEvt = 0; nEvt < nentries; ++ nEvt) {

    if(nEvt%10000 == 0)cout<<nEvt<<"/"<<nentries<<endl;
    if(printDebug)cout<<"nEvt = "<<nEvt<<endl;

    for(int t = 0;t<N;t++)
      jetTree[t]->GetEntry(nEvt);
    
    //! primary event cleaning
    if(fabs(vz_F)>15 || !pprimaryvertexFilter_F ||
       !pVertexFilterCutGplus_F || !pHBHENoiseFilter_F)
      continue;
    if(!jet40_F && !jet60_F && !jet80_F)
      continue;
    
    double ev_weight = 0.0;
    double maxTrigpT = 0.0;
    if(trgObjpt_40->size()>0){
      for(unsigned it = 0; it<trgObjpt_40->size(); ++it){
	if(maxTrigpT < trgObjpt_40->at(it))
	  maxTrigpT = trgObjpt_40->at(it);
      }
    }
    if(trgObjpt_60->size()>0){
      for(unsigned it = 0; it<trgObjpt_60->size(); ++it){
	if(maxTrigpT < trgObjpt_60->at(it))
	  maxTrigpT = trgObjpt_60->at(it);
      }
    }
    if(trgObjpt_80->size()>0){
      for(unsigned it = 0; it<trgObjpt_80->size(); ++it){
	if(maxTrigpT < trgObjpt_80->at(it))
	  maxTrigpT = trgObjpt_80->at(it);
      }
    }

    int presclDes[] = {jet40_F, jet60_F, jet80_F};
    int presclWgt[] = {jet40_p_F*jet40_p_l1_F,
		       jet60_p_F*jet60_p_l1_F,
		       jet80_p_F*jet80_p_l1_F};
    
    ev_weight = eventPrescl(presclDes, presclWgt, maxTrigpT);
    if(printDebug){
      cout<<"Prescl weight = "<<ev_weight<<endl;
    }

    //! fill in the event histograms
    hVz->Fill(vz_F, ev_weight);
    hhiBin->Fill(hiBin_F, ev_weight);
    hNpix->Fill(hiNpix_F, ev_weight);
    hHF->Fill(hiHF_F, ev_weight);
    hHFplus->Fill(hiHFplus_F, ev_weight);
    hHFminus->Fill(hiHFminus_F, ev_weight);
    hZDC->Fill(hiZDC_F, ev_weight);
    hZDCplus->Fill(hiZDCplus_F, ev_weight);
    hZDCminus->Fill(hiZDCminus_F, ev_weight);
    hHFhit->Fill(hiHFhit_F, ev_weight);
    hHFhitplus->Fill(hiHFhitplus_F, ev_weight);
    hHFhitminus->Fill(hiHFhitminus_F, ev_weight);
    hNtracks->Fill(hiNtracks_F, ev_weight);

    double ptlead = 0.0;
    double ptsublead = 0.0;
    int counter = 0;

    //! Loop over the jets 
    for(int jet = 0; jet<nref_F; ++jet){
      //! no JetID at the moment
      if(fabs(eta_F[jet])>3.0) continue;

      if(counter == 0) ptlead = pt_F[jet];
      if(counter == 1) ptsublead = pt_F[jet];
      counter++;
      
      hjtpt->Fill(pt_F[jet], ev_weight);
      hjteta->Fill(eta_F[jet], ev_weight);
      hjtphi->Fill(phi_F[jet], ev_weight);
      hjtm->Fill(m_F[jet], ev_weight);
      hCHF->Fill(jtPfCHF_F[jet], ev_weight);
      hNHF->Fill(jtPfNHF_F[jet], ev_weight);
      hCEF->Fill(jtPfCEF_F[jet], ev_weight);
      hNEF->Fill(jtPfNEF_F[jet], ev_weight);
      hMUF->Fill(jtPfMUF_F[jet], ev_weight);
      hCHM->Fill(jtPfCHM_F[jet], ev_weight);
      hNHM->Fill(jtPfNHM_F[jet], ev_weight);
      hCEM->Fill(jtPfCEM_F[jet], ev_weight);
      hNEM->Fill(jtPfNEM_F[jet], ev_weight);
      hMUM->Fill(jtPfMUM_F[jet], ev_weight);
    }// jet loop
    
    double Aj = 0.0;
    if(ptlead!=0.0 && ptsublead!=0.0){
      Aj = (double)(ptlead - ptsublead)/(ptlead + ptsublead);
      hAj->Fill(Aj, ev_weight);
      hAjvshiBin->Fill(hiBin_F, Aj, ev_weight);
      hAjvsNpix->Fill(hiNpix_F, Aj, ev_weight);
      hAjvsNtracks->Fill(hiNtracks_F, Aj, ev_weight);
      hAjvsHF->Fill(hiHF_F, Aj, ev_weight);
      hAjvsHFplus->Fill(hiHFplus_F, Aj, ev_weight);
      hAjvsHFminus->Fill(hiHFminus_F, Aj, ev_weight);
      hAjvsZDC->Fill(hiZDC_F, Aj, ev_weight);
      hAjvsZDCplus->Fill(hiZDCplus_F, Aj, ev_weight);
      hAjvsZDCminus->Fill(hiZDCminus_F, Aj, ev_weight);
      hAjvsHFhit->Fill(hiHFhit_F, Aj, ev_weight);
      hAjvsHFhitplus->Fill(hiHFhitplus_F, Aj, ev_weight);
      hAjvsHFhitminus->Fill(hiHFhitminus_F, Aj, ev_weight);
    }
    
    if(printDebug)cout<<endl;

  }// event loop

  fout->Write();
  
  timer.Stop();
  cout<<"Macro finished: "<<endl;
  cout<<"CPU time (min)  = "<<(Float_t)timer.CpuTime()/60<<endl;
  cout<<"Real time (min) = "<<(Float_t)timer.RealTime()/60<<endl;
  
}//macro end
