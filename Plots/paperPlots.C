#include "plot.h"
#include "utilities.h"

using namespace std;

void paperPlots()
{

  gStyle->SetOptStat(0);
  TH1::SetDefaultSumw2();
  TFile *fData, *fPY8, *fPY6;
  fData = TFile::Open("../Analysis/pPb_8p16TeV_Data_ak4PF.root");
  // fPY6 = TFile::Open("rootfiles/radialmoment_PY6_pp8p16TeV.root");
  // fPY8 = TFile::Open("rootfiles/radialmoment_PY8_pp8p16TeV.root");

  //! prepare the histograms!
  // TH1F * hData[4][5][3][3][3];
  // TH1F * hPY6A1D[4][3][3][3][2];  
  // TH1F * hPY6A2D[4][3][3][3][2];  
  // TH1F * hPY6B1D[3][3][3][3];  
  // TH2F * hPY6B2D[4][3][3][3];    
  // TH1F * hPY8[4][5][3][3][3];

  int _centedges[] = {0, 40, 80, 120, 200};
  const int _ncentbins = sizeof(_centedges)/sizeof(int)-1;
  std::string centbins[] = {"0_cent_20", "20_cent_40", "40_cent_60", "60_cent_100"};

  double _ptedges[] = {50., 100., 300., 500, 2000.};
  const int _nptbins = sizeof(_ptedges)/sizeof(double)-1;
  std::string ptbins[] = {"50_pt_100", "100_pt_300", "300_pt_500", "pt_gt500"};
  
  double _etaedges[] = {-3.5, -2.0, -0.5, 0.5, 2.0, 3.5};
  const int _netabins = sizeof(_etaedges)/sizeof(double)-1;
  std::string etabins[] = {"m3p5_eta_m2p0", "m2p0_eta_m0p5",
			   "m0p5_eta_p0p5",
			   "p0p5_eta_p2p0", "p2p0_eta_3p5"};

  double _betaValues[] = {0.5, 1., 2., 3.};
  double _radMomMax[] = {0.8, 0.5, 0.3, 0.05};
  const int _nbeta = sizeof(_betaValues)/sizeof(double);
  std::string moments[] = {"halfMoment", "firstMoment", "secondMoment", "thirdMoment"};
  std::string mass[] = {"Inclusive jet mass", "Inclusive groomed jet mass"};

  const int nRebin = 5;
  
  //! 1D histograms 
  // std::string histoType1D[] = {"JetGroomEffect", "JetConst_ungrm", "JetConst_grm", "hJetMass_ungrm", "hJetMass_grm"};
  // const int _nHisto1D = sizeof(histoType1D)/sizeof(std::string);
  
  TH1F * hData_Moment_grm[_ncentbins][_nptbins][_netabins][_nbeta];
  TH1F * hData_Moment_ungrm[_ncentbins][_nptbins][_netabins][_nbeta];  
  TH1F * hData_GroomEffect[_ncentbins][_nptbins][_netabins];
  TH1F * hData_Mass_grm[_ncentbins][_nptbins][_netabins];
  TH1F * hData_Mass_ungrm[_ncentbins][_nptbins][_netabins];			    
  TH1F * hPY6_Moment_grm[_ncentbins][_nptbins][_netabins][_nbeta];
  TH1F * hPY6_Moment_ungrm[_ncentbins][_nptbins][_netabins][_nbeta];  
  TH1F * hPY6_GroomEffect[_ncentbins][_nptbins][_netabins];
  TH1F * hPY6_Mass_grm[_ncentbins][_nptbins][_netabins];
  TH1F * hPY6_Mass_ungrm[_ncentbins][_nptbins][_netabins];			    
  
  //! 2D histograms 
  // std::string histoTypeA[] = {"JetConst_ungrm", "JetConst_grm", "hJetMass_vs_RadMom_ungrm", "hJetMass_vs_RadMom_grm", "hJetMoment_vs_pt_ungrm", "hJetMoment_vs_pt_grm"};
  // std::string histoTypeB[] = {"JetGroomEffect", "hJetMass_ungrm", "hJetMass_grm", "JetConstDistribution_ungrm", "JetConstDistribution_grm", "hJetMass_vs_pt_ungrm", "hJetMass_vs_pt_ungrm"};

  std::string centLeg[] = {"0-20%", "20-40%", "40-60%", "60-100%"};
  std::string ptLeg[] = {"50 < p_{T} < 100 [GeV/c]", "100 < p_{T} < 300 [GeV/c]", "300 < p_{T} < 500 [GeV/c]", " 500 < p_{T} [GeV/c]"};
  std::string etaLeg[] = {"-3.5 < #eta < -2.0", "-2.0 < #eta < -0.5", "-0.5 < #eta < 0.5", "0.5 < #eta < 2.0", "2.0 < #eta < 3.5"};
  std::string momentsLeg[] = {"0.5^{th} Moment", "1^{st} Moment", "2^{nd} Moment", "3^{rd} Moment"};

  std::string radius = {"R = 0.4"};
  std::string jetdes[] = {"anti k_{T} PF Jets",  "anti k_{T} groomed PF Jets"};  
  
  for(int w = 0; w < _ncentbins; ++w){
    for(int x = 0; x < _nptbins; ++x){
      for(int y = 0; y < _netabins; ++y){
	// for(int v = 0; v < 7; ++v){
	//   if(v < 3){
	//     hPY6B1D[w][x][y] = (TH1F*)fPY6->Get(Form("h%s_%s_%s_%s", histoTypeB[v].c_str(), centbins[w].c_str(), ptbins[x].c_str(), etabins[y].c_str()));
	//   }
	//   else{
	//     hPY6B2D[v-1][w][x][y] = (TH2F*)fPY6->Get(Form("h%s_%s_%s_%s", histoTypeB[v].c_str(), centbins[w].c_str(), ptbins[x].c_str(), etabins[y].c_str()));
	//   }
	// }
	// for(int v = 0; v < 8; ++v){
	//   for(int z = 0; z < _nbeta; ++z){
	//     if(v < 4){
	//       hPY6A1D[v][w][x][y][z] = (TH1F*)fPY6->Get(Form("h%s_%s_%s_%s_%s", histoTypeA[v].c_str(), centbins[w].c_str(), ptbins[x].c_str(), etabins[y].c_str(), moments[z].c_str()));
	//       hPY6A1D[v][w][x][y][z]->Rebin(4);
	//     // divideBinWidth(hPY6[v][w][x][y][z]);
	//     }else {
	//       hPY6A2D[v][w][x][y][z] = (TH2F*)fPY6->Get(Form("h%s_%s_%s_%s_%s", histoTypeA[v].c_str(), centbins[w].c_str(), ptbins[x].c_str(), etabins[y].c_str(), moments[z].c_str()));
	//     }
	//   }	
	// }
	hData_GroomEffect[w][x][y] = (TH1F*)fData->Get(Form("hJetGroomEffect_%s_%s_%s", centbins[w].c_str(), ptbins[x].c_str(), etabins[y].c_str()));
	hData_GroomEffect[w][x][y]->Print("base");
	// hData_GroomEffect[w][x][y]->Rebin(nRebin);
	// divideBinWidth(hData_GroomEffect[w][x][y]);
	if(hData_GroomEffect[w][x][y]->GetEntries()>0)
	  hData_GroomEffect[w][x][y]->Scale(1./hData_GroomEffect[w][x][y]->Integral());
	hData_Mass_grm[w][x][y] = (TH1F*)fData->Get(Form("hJetMass_grm_%s_%s_%s", centbins[w].c_str(), ptbins[x].c_str(), etabins[y].c_str()));
	hData_Mass_grm[w][x][y]->Print("base");
	// hData_Mass_grm[w][x][y]->Rebin(nRebin);
	// divideBinWidth(hData_Mass_grm[w][x][y]);
	if(hData_Mass_grm[w][x][y]->GetEntries()>0)
	  hData_Mass_grm[w][x][y]->Scale(1./hData_Mass_grm[w][x][y]->Integral());
	hData_Mass_ungrm[w][x][y] = (TH1F*)fData->Get(Form("hJetMass_ungrm_%s_%s_%s", centbins[w].c_str(), ptbins[x].c_str(), etabins[y].c_str()));	
	hData_Mass_ungrm[w][x][y]->Print("base");
	// hData_Mass_ungrm[w][x][y]->Rebin(nRebin);
	// divideBinWidth(hData_Mass_ungrm[w][x][y]);
	if(hData_Mass_ungrm[w][x][y]->GetEntries()>0)
	  hData_Mass_ungrm[w][x][y]->Scale(1./hData_Mass_ungrm[w][x][y]->Integral());
	
	// hPY6_GroomEffect[w][x][y] = (TH1F*)fPY6->Get(Form("hJetGroomEffect_%s_%s_%s", centbins[w].c_str(), ptbins[x].c_str(), etabins[y].c_str()));
	// hPY6_GroomEffect[w][x][y]->Print("base");
	// hPY6_GroomEffect[w][x][y]->Rebin(nRebin);
	// divideBinWidth(hPY6_GroomEffect[w][x][y]);
	// hPY6_Mass_grm[w][x][y] = (TH1F*)fPY6->Get(Form("hJetMass_grm_%s_%s_%s", centbins[w].c_str(), ptbins[x].c_str(), etabins[y].c_str()));
	// hPY6_Mass_grm[w][x][y]->Print("base");
	// hPY6_Mass_grm[w][x][y]->Rebin(nRebin);
	// divideBinWidth(hPY6_Mass_grm[w][x][y]);
	// hPY6_Mass_ungrm[w][x][y] = (TH1F*)fPY6->Get(Form("hJetMass_ungrm_%s_%s_%s", centbins[w].c_str(), ptbins[x].c_str(), etabins[y].c_str()));	
        // hPY6_Mass_ungrm[w][x][y]->Print("base");
	for(int z = 0; z < _nbeta; ++z){
	  cout<<Form("hJetConst_grm_%s_%s_%s_%s", centbins[w].c_str(), ptbins[x].c_str(), etabins[y].c_str(), moments[z].c_str())<<endl;
	  hData_Moment_grm[w][x][y][z] = (TH1F*)fData->Get(Form("hJetConst_grm_%s_%s_%s_%s", centbins[w].c_str(), ptbins[x].c_str(), etabins[y].c_str(), moments[z].c_str()));
	  // hData_Moment_grm[w][x][y][z]->Print("base");
	  // hData_Moment_grm[w][x][y][z]->Rebin(nRebin);
	  divideBinWidth(hData_Moment_grm[w][x][y][z]);
	  if(hData_Moment_grm[w][x][y][z]->GetEntries()>0)
	    hData_Moment_grm[w][x][y][z]->Scale(1./hData_Moment_grm[w][x][y][z]->Integral());
	  hData_Moment_ungrm[w][x][y][z] = (TH1F*)fData->Get(Form("hJetConst_ungrm_%s_%s_%s_%s", centbins[w].c_str(), ptbins[x].c_str(), etabins[y].c_str(), moments[z].c_str()));
	  hData_Moment_ungrm[w][x][y][z]->Print("base");
	  // hData_Moment_ungrm[w][x][y][z]->Rebin(nRebin);
	  // divideBinWidth(hData_Moment_ungrm[w][x][y][z]);
	  if(hData_Moment_ungrm[w][x][y][z]->GetEntries()>0)
	    hData_Moment_ungrm[w][x][y][z]->Scale(1./hData_Moment_ungrm[w][x][y][z]->Integral());
	  // hPY6_Moment_grm[w][x][y][z] = (TH1F*)fPY6->Get(Form("hJetConst_grm_%s_%s_%s_%s", centbins[w].c_str(), ptbins[x].c_str(), etabins[y].c_str(), moments[z].c_str()));
	  // hPY6_Moment_grm[w][x][y][z]->Print("base");
	  // hPY6_Moment_grm[w][x][y][z]->Rebin(nRebin);
	  // divideBinWidth(hPY6_Moment_grm[w][x][y][z]);
	  // hPY6_Moment_ungrm[w][x][y][z] = (TH1F*)fPY6->Get(Form("hJetConst_ungrm_%s_%s_%s_%s", centbins[w].c_str(), ptbins[x].c_str(), etabins[y].c_str(), moments[z].c_str()));
	  // hPY6_Moment_ungrm[w][x][y][z]->Print("base");
	  // hPY6_Moment_ungrm[w][x][y][z]->Rebin(nRebin);
	  // divideBinWidth(hPY6_Moment_ungrm[w][x][y][z]);
	}
      }
    }
  }
  
  //! prepare the canvas!
  TCanvas * cMoments_grm[_nbeta];  
  TCanvas * cMoments_ungrm[_nbeta];  
  TH1F * hDummyMoment[_nbeta];
  
  TLegend * legMC = myLegend(0.2, 0.6, 0.3, 0.9);
  TLegend * legData = myLegend(0.3, 0.6, 0.4, 0.9);
  TLegend * legComb = myLegend(0.35, 0.6, 0.9, 0.9);

  for(int z = 0; z < _nbeta; ++z){
    cMoments_grm[z] = new TCanvas(Form("%s_grm", moments[z].c_str()),"",400.*_ncentbins, 400.*_nptbins);
    makeMultiPanelCanvas(cMoments_grm[z], _ncentbins, _nptbins, 0.2, 0.2, 0.2, 0.2, 0.08);
    hDummyMoment[z]= new TH1F(Form("hDummyMoment_%s", moments[z].c_str()),"",100, 0, 1);
    hDummyMoment[z]->SetTitle("");
    hDummyMoment[z]->SetYTitle("1/N^{jets} dN/drm");
    hDummyMoment[z]->SetXTitle(Form("%s", momentsLeg[z].c_str()));
    hDummyMoment[z]->SetAxisRange(0.0001, 0.6, "X");
    hDummyMoment[z]->SetAxisRange(0.001, 1, "y");
    hDummyMoment[z]->GetXaxis()->SetNdivisions(505);
    hDummyMoment[z]->GetXaxis()->SetLabelFont(43);
    hDummyMoment[z]->GetXaxis()->SetLabelSize(25);
    hDummyMoment[z]->GetXaxis()->SetTitleSize(25);
    hDummyMoment[z]->GetXaxis()->SetTitleOffset(3.5);
    hDummyMoment[z]->GetXaxis()->SetTitleFont(43);
    hDummyMoment[z]->GetXaxis()->CenterTitle(true);
    hDummyMoment[z]->GetYaxis()->CenterTitle(true);
    hDummyMoment[z]->GetYaxis()->SetNdivisions(610);
    hDummyMoment[z]->GetYaxis()->SetLabelFont(43);
    hDummyMoment[z]->GetYaxis()->SetLabelSize(25);
    hDummyMoment[z]->GetYaxis()->SetTitleSize(25);
    hDummyMoment[z]->GetYaxis()->SetTitleOffset(4.0);
    hDummyMoment[z]->GetYaxis()->SetTitleFont(43);
    int padCounter = 1;
    for(int x = 0; x < _nptbins; ++x){
      for(int w = 0; w < _ncentbins; ++w){
	cMoments_grm[z]->cd(padCounter);
	gPad->SetTickx();
	gPad->SetTicky();
	gStyle->SetErrorX(0);
	cMoments_grm[z]->cd(padCounter)->SetLogy();
	cMoments_grm[z]->cd(padCounter)->SetLogx();
	hDummyMoment[z]->Draw();
	if(z == 0 && padCounter == 1){
	  legMC->AddEntry("","","");
	  legData->AddEntry("","","");
	  legComb->AddEntry("", "", "");
	}
	for(int y = 0; y < _netabins; ++y){
	  // hPY6_Moment_grm[w][x][y][z]->SetLineStyle(2);
	  // hPY6_Moment_grm[w][x][y][z]->SetLineWidth(1);
	  // hPY6_Moment_grm[w][x][y][z]->SetLineColor(color[y]);
	  // hPY6_Moment_grm[w][x][y][z]->Draw("same l");
	  hData_Moment_grm[w][x][y][z]->SetMarkerStyle(33);
	  hData_Moment_grm[w][x][y][z]->SetMarkerSize(1.1);
	  hData_Moment_grm[w][x][y][z]->SetMarkerColor(color[y]);
	  hData_Moment_grm[w][x][y][z]->Draw("same p");
	  if(z == 0 && padCounter == 1){
	    // legMC->AddEntry(hPY6_Moment_grm[w][x][y][z], "", "l");
	    legData->AddEntry(hData_Moment_grm[w][x][y][z], "", "p");
	    legComb->AddEntry("",Form("%s",etaLeg[y].c_str()), "");
	  }
	}
	if(padCounter <= _ncentbins)
	  drawText(Form("%s",centLeg[w].c_str()), 0.4, 0.8, 16);
	if(padCounter%_ncentbins == 1)
	  drawText(Form("%s",ptLeg[x].c_str()), 0.4, 0.7, 16);
	if(padCounter == _ncentbins){
	  drawText("pPb (8.16 TeV)", 0.6, 0.95, 16);
	  drawText(Form("%s",jetdes[1].c_str()),0.1, 0.7, 16);
	  drawText(Form("%s",radius.c_str()), 0.3, 0.6, 16);
	}
	if(x == 0 && w == 0)
	  drawText("CMS Preliminary", 0.3, 0.95, 16);
	padCounter++;
      }
    }
    // legMC->SetTextSize(0.04);
    // legMC->Draw();
    legData->SetTextSize(0.04);
    legData->Draw();
    legComb->SetTextSize(0.04);
    legComb->Draw();
    // drawText("PY6", 0.2, 0.8, 16);
    drawText("Data", 0.3, 0.8, 16);
    cMoments_grm[z]->SaveAs(Form("plots/Data_Grm_RadMom_%s.pdf", moments[z].c_str()),"RECREATE");
  }

  for(int z = 0; z < _nbeta; ++z){
    cMoments_ungrm[z] = new TCanvas(Form("%s_ungrm", moments[z].c_str()),"",400.*_ncentbins, 400.*_nptbins);
    makeMultiPanelCanvas(cMoments_ungrm[z], _ncentbins, _nptbins, 0.2, 0.2, 0.2, 0.2, 0.08);
    int padCounter = 1;
    for(int x = 0; x < _nptbins; ++x){
      for(int w = 0; w < _ncentbins; ++w){
	cMoments_ungrm[z]->cd(padCounter);
	gPad->SetTickx();
	gPad->SetTicky();
	gStyle->SetErrorX(0);
	cMoments_ungrm[z]->cd(padCounter)->SetLogy();
	cMoments_ungrm[z]->cd(padCounter)->SetLogx();
	hDummyMoment[z]->Draw();
	for(int y = 0; y < _netabins; ++y){
	  // hPY6_Moment_ungrm[w][x][y][z]->SetLineStyle(2);
	  // hPY6_Moment_ungrm[w][x][y][z]->SetLineWidth(1);
	  // hPY6_Moment_ungrm[w][x][y][z]->SetLineColor(color[y]);
	  // hPY6_Moment_ungrm[w][x][y][z]->Draw("same l");
	  hData_Moment_ungrm[w][x][y][z]->SetMarkerStyle(33);
	  hData_Moment_ungrm[w][x][y][z]->SetMarkerSize(1.1);
	  hData_Moment_ungrm[w][x][y][z]->SetMarkerColor(color[y]);
	  hData_Moment_ungrm[w][x][y][z]->Draw("same p");
	}
	if(padCounter <= _ncentbins)
	  drawText(Form("%s",centLeg[w].c_str()), 0.4, 0.8, 16);
	if(padCounter%_ncentbins == 1)
	  drawText(Form("%s",ptLeg[x].c_str()), 0.4, 0.7, 16);
	if(padCounter == _ncentbins){
	  drawText("pPb (8.16 TeV)", 0.6, 0.95, 16);
	  drawText(Form("%s",jetdes[0].c_str()),0.1, 0.7, 16);
	  drawText(Form("%s",radius.c_str()), 0.3, 0.6, 16);
	}
	if(x == 0 && w == 0)
	  drawText("CMS Preliminary", 0.3, 0.95, 16);
	padCounter++;
      }
    }
    // legMC->Draw();
    legData->Draw();
    legComb->Draw();
    cMoments_ungrm[z]->SaveAs(Form("plots/Data_unGrm_RadMom_%s.pdf", moments[z].c_str()),"RECREATE");
  }

  

  
  //! prepare the canvas!
  TCanvas * cGroomEffect;  
  TH1F * hDummy_GroomEffect;  
  std::string xtitle[] = {"p^{grm jet}_{T}/p^{jet}_{T}", "r from jet axis", "r from jet axis"};
  std::string ytitle[] = {"1/N_{jets} dN/d(p^{grm jet}_{T}/p^{jet}_{T})", "p^{constituent}_{T}/p^{jet}_{T}", "p^{constituent}_{T}/p^{grm jet}_{T}"};
  
  cGroomEffect = new TCanvas("cgroomEffect","",400.*_ncentbins, 400.*_nptbins);
  makeMultiPanelCanvas(cGroomEffect, _ncentbins, _nptbins, 0.2, 0.2, 0.2, 0.2, 0.08);
  hDummy_GroomEffect= new TH1F("hDummy_GroomEffect","",100, 0, 1);
  hDummy_GroomEffect->SetTitle("");
  hDummy_GroomEffect->SetYTitle(Form("%s",ytitle[0].c_str()));
  hDummy_GroomEffect->SetXTitle(Form("%s",xtitle[0].c_str()));
  hDummy_GroomEffect->SetAxisRange(0.5, 1, "X");
  hDummy_GroomEffect->SetAxisRange(5e-4, 5e-1, "Y");
  hDummy_GroomEffect->GetXaxis()->SetNdivisions(505);
  hDummy_GroomEffect->GetXaxis()->SetLabelFont(43);
  hDummy_GroomEffect->GetXaxis()->SetLabelSize(25);
  hDummy_GroomEffect->GetXaxis()->SetTitleSize(25);
  hDummy_GroomEffect->GetXaxis()->SetTitleOffset(3.5);
  hDummy_GroomEffect->GetXaxis()->SetTitleFont(43);
  hDummy_GroomEffect->GetXaxis()->CenterTitle(true);
  hDummy_GroomEffect->GetYaxis()->CenterTitle(true);
  hDummy_GroomEffect->GetYaxis()->SetNdivisions(610);
  hDummy_GroomEffect->GetYaxis()->SetLabelFont(43);
  hDummy_GroomEffect->GetYaxis()->SetLabelSize(25);
  hDummy_GroomEffect->GetYaxis()->SetTitleSize(25);
  hDummy_GroomEffect->GetYaxis()->SetTitleOffset(4.0);
  hDummy_GroomEffect->GetYaxis()->SetTitleFont(43);
    
  int padCounter = 1;
  for(int x = 0; x < _nptbins; ++x){
    for(int w = 0; w < _ncentbins; ++w){
      cGroomEffect->cd(padCounter);
      gPad->SetTickx();
      gPad->SetTicky();
      gStyle->SetErrorX(0);
      cGroomEffect->cd(padCounter)->SetLogy();
      cGroomEffect->cd(padCounter)->SetLogx();
      hDummy_GroomEffect->Draw();
      for(int y = 0; y < _netabins; ++y){
	// hPY6_GroomEffect[w][x][y]->SetLineStyle(2);
	// hPY6_GroomEffect[w][x][y]->SetLineWidth(2);
	// hPY6_GroomEffect[w][x][y]->SetLineColor(color[y]);
	// hPY6_GroomEffect[w][x][y]->Draw("same l");
	hData_GroomEffect[w][x][y]->SetMarkerStyle(20);
	hData_GroomEffect[w][x][y]->SetMarkerSize(1.1);
	hData_GroomEffect[w][x][y]->SetMarkerColor(color[y]);
	hData_GroomEffect[w][x][y]->Draw("same p");
      }
      if(padCounter <= _ncentbins)
	drawText(Form("%s",centLeg[w].c_str()), 0.4, 0.8, 16);
      if(padCounter%_ncentbins == 1)
	drawText(Form("%s",ptLeg[x].c_str()), 0.3, 0.25, 16);
      if(padCounter == _ncentbins){
	drawText("(8.16 TeV)", 0.6, 0.95, 16);
	drawText(Form("%s",jetdes[0].c_str()),0.1, 0.7, 16);
	drawText(Form("%s",radius.c_str()), 0.3, 0.6, 16);
      }
      // if(x == 0 && w == 0)
      // 	drawText("CMS Simulation", 0.3, 0.95, 16);
      padCounter++;
    }
  }
  // legMC->Draw();
  legData->Draw();
  legComb->Draw();
  cGroomEffect->SaveAs("plots/Data_GroomEffect.pdf", "RECREATE");


  //! prepare the canvas!
  TCanvas * cMass_grm;  
  TH1F * hDummy_Mass_grm;    
  cMass_grm = new TCanvas("cMass_grm","",400.*_ncentbins, 400.*_nptbins);
  makeMultiPanelCanvas(cMass_grm, _ncentbins, _nptbins, 0.2, 0.2, 0.2, 0.2, 0.08);
  hDummy_Mass_grm= new TH1F("hDummy_Mass_grm","",100, 0, 50);
  hDummy_Mass_grm->SetTitle("");
  hDummy_Mass_grm->SetXTitle("SD Groomed Jet Mass m_{J} [GeV/c^2]");
  hDummy_Mass_grm->SetYTitle("1/N^{jets} dN/dm_{J}");
  hDummy_Mass_grm->SetAxisRange(0, 50, "X");
  hDummy_Mass_grm->SetAxisRange(0, 0.2, "Y");
  hDummy_Mass_grm->GetXaxis()->SetNdivisions(505);
  hDummy_Mass_grm->GetXaxis()->SetLabelFont(43);
  hDummy_Mass_grm->GetXaxis()->SetLabelSize(25);
  hDummy_Mass_grm->GetXaxis()->SetTitleSize(25);
  hDummy_Mass_grm->GetXaxis()->SetTitleOffset(3.5);
  hDummy_Mass_grm->GetXaxis()->SetTitleFont(43);
  hDummy_Mass_grm->GetXaxis()->CenterTitle(true);
  hDummy_Mass_grm->GetYaxis()->CenterTitle(true);
  hDummy_Mass_grm->GetYaxis()->SetNdivisions(610);
  hDummy_Mass_grm->GetYaxis()->SetLabelFont(43);
  hDummy_Mass_grm->GetYaxis()->SetLabelSize(25);
  hDummy_Mass_grm->GetYaxis()->SetTitleSize(25);
  hDummy_Mass_grm->GetYaxis()->SetTitleOffset(4.0);
  hDummy_Mass_grm->GetYaxis()->SetTitleFont(43);
    
  padCounter = 1;
  for(int x = 0; x < _nptbins; ++x){
    for(int w = 0; w < _ncentbins; ++w){
      cMass_grm->cd(padCounter);
      gPad->SetTickx();
      gPad->SetTicky();
      gStyle->SetErrorX(0);
      cMass_grm->cd(padCounter)->SetLogy();
      cMass_grm->cd(padCounter)->SetLogx();
      hDummy_Mass_grm->Draw();
      for(int y = 0; y < _netabins; ++y){
	// hPY6_Mass_grm[w][x][y]->SetLineStyle(2);
	// hPY6_Mass_grm[w][x][y]->SetLineWidth(2);
	// hPY6_Mass_grm[w][x][y]->SetLineColor(color[y]);
	// hPY6_Mass_grm[w][x][y]->Draw("same l");
	hData_Mass_grm[w][x][y]->SetMarkerStyle(20);
	hData_Mass_grm[w][x][y]->SetMarkerSize(1.1);
	hData_Mass_grm[w][x][y]->SetMarkerColor(color[y]);
	hData_Mass_grm[w][x][y]->Draw("same p");
      }
      if(padCounter <= _ncentbins)
	drawText(Form("%s",centLeg[w].c_str()), 0.4, 0.8, 16);
      if(padCounter%_ncentbins == 1)
	drawText(Form("%s",ptLeg[x].c_str()), 0.3, 0.25, 16);
      if(padCounter == _ncentbins){
	drawText("(8.16 TeV)", 0.6, 0.95, 16);
	drawText(Form("%s",jetdes[0].c_str()),0.1, 0.7, 16);
	drawText(Form("%s",radius.c_str()), 0.3, 0.6, 16);
      }
      if(x == 0 && w == 0)
	drawText("CMS Preliminary", 0.3, 0.95, 16);
      padCounter++;
    }
  }
  legMC->Draw();
  legData->Draw();
  legComb->Draw();
  cMass_grm->SaveAs("plots/Data_Mass_grm.pdf", "RECREATE");

  
  //! prepare the canvas!
  TCanvas * cMass_ungrm;  
  TH1F * hDummy_Mass_ungrm;    
  cMass_ungrm = new TCanvas("cMass_ungrm","",400.*_ncentbins, 400.*_nptbins);
  makeMultiPanelCanvas(cMass_ungrm, _ncentbins, _nptbins, 0.2, 0.2, 0.2, 0.2, 0.08);
  hDummy_Mass_ungrm= new TH1F("hDummy_Mass_ungrm","",100, 0, 50);
  hDummy_Mass_ungrm->SetTitle("");
  hDummy_Mass_ungrm->SetXTitle("Jet Mass m_{J} [GeV/c^2]");
  hDummy_Mass_ungrm->SetYTitle("1/N^{jets} dN/dm_{J}");
  hDummy_Mass_ungrm->SetAxisRange(0, 50, "X");
  hDummy_Mass_ungrm->SetAxisRange(0, 0.2, "Y");
  hDummy_Mass_ungrm->GetXaxis()->SetNdivisions(505);
  hDummy_Mass_ungrm->GetXaxis()->SetLabelFont(43);
  hDummy_Mass_ungrm->GetXaxis()->SetLabelSize(25);
  hDummy_Mass_ungrm->GetXaxis()->SetTitleSize(25);
  hDummy_Mass_ungrm->GetXaxis()->SetTitleOffset(3.5);
  hDummy_Mass_ungrm->GetXaxis()->SetTitleFont(43);
  hDummy_Mass_ungrm->GetXaxis()->CenterTitle(true);
  hDummy_Mass_ungrm->GetYaxis()->CenterTitle(true);
  hDummy_Mass_ungrm->GetYaxis()->SetNdivisions(610);
  hDummy_Mass_ungrm->GetYaxis()->SetLabelFont(43);
  hDummy_Mass_ungrm->GetYaxis()->SetLabelSize(25);
  hDummy_Mass_ungrm->GetYaxis()->SetTitleSize(25);
  hDummy_Mass_ungrm->GetYaxis()->SetTitleOffset(4.0);
  hDummy_Mass_ungrm->GetYaxis()->SetTitleFont(43);
    
  padCounter = 1;
  for(int x = 0; x < _nptbins; ++x){
    for(int w = 0; w < _ncentbins; ++w){
      cMass_ungrm->cd(padCounter);
      gPad->SetTickx();
      gPad->SetTicky();
      gStyle->SetErrorX(0);
      cMass_ungrm->cd(padCounter)->SetLogy();
      cMass_ungrm->cd(padCounter)->SetLogx();
      hDummy_Mass_ungrm->Draw();
      for(int y = 0; y < _netabins; ++y){
	// hPY6_Mass_ungrm[w][x][y]->SetLineStyle(2);
	// hPY6_Mass_ungrm[w][x][y]->SetLineWidth(2);
	// hPY6_Mass_ungrm[w][x][y]->SetLineColor(color[y]);
	// hPY6_Mass_ungrm[w][x][y]->Draw("same l");
	hData_Mass_ungrm[w][x][y]->SetMarkerStyle(20);
	hData_Mass_ungrm[w][x][y]->SetMarkerSize(1.1);
	hData_Mass_ungrm[w][x][y]->SetMarkerColor(color[y]);
	hData_Mass_ungrm[w][x][y]->Draw("same p");
      }
      if(padCounter <= _ncentbins)
	drawText(Form("%s",centLeg[w].c_str()), 0.4, 0.8, 16);
      if(padCounter%_ncentbins == 1)
	drawText(Form("%s",ptLeg[x].c_str()), 0.3, 0.25, 16);
      if(padCounter == _ncentbins){
	drawText("(8.16 TeV)", 0.6, 0.95, 16);
	drawText(Form("%s",jetdes[0].c_str()),0.1, 0.7, 16);
	drawText(Form("%s",radius.c_str()), 0.3, 0.6, 16);
      }
      if(x == 0 && w == 0)
	drawText("CMS Preliminary", 0.3, 0.95, 16);
      padCounter++;
    }
  }
  legMC->Draw();
  legData->Draw();
  legComb->Draw();
  cMass_ungrm->SaveAs("plots/Data_Mass_ungrm.pdf", "RECREATE");


  
}
