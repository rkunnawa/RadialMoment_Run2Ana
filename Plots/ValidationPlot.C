//! macro to plot variables data vs MC with validation

{
  TFile * fData = TFile::Open("../Analysis/pPb_5p02TeV_MinBias_ak4PF.root");
  TFile * fMC = TFile::Open("../Analysis/pPb_PYEPOS_5p02TeV_MinBiasMC_ak4PF.root");

  const int N = 39;
  
  std::string var[N] = {
    "Vz", "hiBin", "Npix", "HF", "HFplus",	"HFminus", "ZDC", "ZDCplus", "ZDCminus",
    "HFhit", "HFhitplus", "HFhitminus", "Ntracks", "jtpt", "jteta", "dijeteta", "jtphi", "jtm",
    "ChSumoverRawpT", "ChMaxoverRawpT", "NeSumoverRawpT", "NeMaxoverRawpT",
    "PhSumoverRawpT", "PhMaxoverRawpT", "ElSumoverRawpT", "ElMaxoverRawpT",
    "MuSumoverRawpT", "MuMaxoverRawpT", "CHF", "NHF", "CEF", "NEF", "MUF", "CHM",
    "NHM", "CEM", "NEM", "MUM", "Aj"};
  
  TH1F * hJetQA[2][2][N];
  TH1F * hRatio[2][N];

  TCanvas * cJetQA[2][N];
  TLegend* leg0 = new TLegend(0.33,0.03,0.63,0.25);
  TLegend* leg1 = new TLegend(0.53,0.63,0.83,0.85);

  // jet ID loop, 0- no jetID, 1- yes JetID
  for(int k = 0; k<1; ++k){

    // variable loop
    for(int i = 0; i<N; ++i){

      cout<<"Running for Variable: "<<var[i].c_str()<<endl;
      hJetQA[k][0][i] = (TH1F*)fData->Get(Form("h%s", var[i].c_str()));
      hJetQA[k][0][i]->Rebin(5);
      hJetQA[k][0][i]->SetMarkerStyle(20);
      hJetQA[k][0][i]->SetMarkerColor(kBlack);
      hJetQA[k][0][i]->SetTitle("");
      hJetQA[k][0][i]->SetXTitle("");
      hJetQA[k][0][i]->SetYTitle("Normalized Counts");
      // hJetQA[k][0][i]->GetXaxis()->CenterTitle();
      // hJetQA[k][0][i]->GetYaxis()->CenterTitle();
      // hJetQA[k][0][1]->GetYaxis()->SetTitleSize(20);
      // hJetQA[k][0][1]->GetYaxis()->SetTitleFont(43);
      // hJetQA[k][0][1]->GetYaxis()->SetTitleOffset(1.55);
      hJetQA[k][1][i] = (TH1F*)fMC->Get(Form("h%s",var[i].c_str()));
      hJetQA[k][1][i]->Rebin(5);
      hJetQA[k][1][i]->SetFillColor(4);
      hJetQA[k][1][i]->SetLineColor(4);
      hJetQA[k][1][i]->SetFillStyle(3354);
      hJetQA[k][1][i]->GetXaxis()->CenterTitle();
      hJetQA[k][1][i]->GetYaxis()->CenterTitle();

      hJetQA[k][1][i]->Scale(1./hJetQA[k][1][i]->Integral());
      hJetQA[k][0][i]->Scale(1./hJetQA[k][0][i]->Integral());
          
      cJetQA[k][i] = new TCanvas(Form("cJetQA_%dwJetID_%d",k,i),"",800,800);

      // Upper plot will be in pad1
      TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
      pad1->SetBottomMargin(0); // Upper and lower plot are joined
      pad1->SetGridx();         // Vertical grid
      pad1->SetGridy();         // horizontal grid
      if(i!=2 || i!=3) pad1->SetLogy();          // Vertical Logy
      pad1->Draw();             // Draw the upper pad: pad1
      pad1->cd();               // pad1 becomes the current pad
      hJetQA[k][0][i]->SetStats(0);          // No statistics on upper plot
      hJetQA[k][0][i]->Draw("h");               // Draw h1
      hJetQA[k][1][i]->Draw("sameh");         // Draw h2 on top of h1

      if(i == 0 && k == 0){
	leg0->AddEntry("","CMS Preliminary","");
	leg0->AddEntry("","p-Pb MinBias #sqrt{s_{NN}}=5.02 TeV","");
	leg0->AddEntry("","ak4PF Jets","");
	leg0->AddEntry(hJetQA[k][0][0],"Data","pl");
	leg0->AddEntry(hJetQA[k][1][0],"MC","lpf");
	leg0->SetTextSize(0.025);  
      }
      if(i == 0 && k == 1){
	leg1->AddEntry("","CMS Preliminary","");
	leg1->AddEntry("","PP #sqrt{s_{NN}}=5.02 TeV, 2015","");
	leg1->AddEntry("","with Jet ID, 13 TeV","");
	leg1->AddEntry("","ak3PF Jets","");
	leg1->AddEntry(hJetQA[k][0][0],"Data","pl");
	leg1->AddEntry(hJetQA[k][1][0],"MC","lpf");
	leg1->SetTextSize(0.025);  
      }

      if(k == 0) leg0->Draw();
      if(k == 1) leg1->Draw();
      
      // lower plot will be in pad
      cJetQA[k][i]->cd();          // Go back to the main canvas before defining pad2
      TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
      pad2->SetTopMargin(0);
      pad2->SetBottomMargin(0.2);
      pad2->SetGridx(); // vertical grid
      pad2->SetGridy(); // horizontal grid
      pad2->Draw();
      pad2->cd();       // pad2 becomes the current pad
    
      // Define the ratio plot
      hRatio[k][i] = (TH1F*)hJetQA[k][0][i]->Clone(Form("hRatio_%s",var[i].c_str()));
      hRatio[k][i]->Divide(hJetQA[k][1][i]);
      hRatio[k][i]->SetTitle("");
      hRatio[k][i]->SetYTitle("Ratio: Data/MC");
      hRatio[k][i]->SetXTitle(Form("%s",var[i].c_str()));
      hRatio[k][i]->SetMinimum(0);  // Define Y ..
      hRatio[k][i]->SetMaximum(2); // .. range
      hRatio[k][i]->SetStats(0);      // No statistics on lower plot
      hRatio[k][i]->SetLineColor(kBlack);
      hRatio[k][i]->SetMarkerStyle(21);
      hRatio[k][i]->Draw("h");
      hRatio[k][i]->GetYaxis()->SetNdivisions(505);
      hRatio[k][i]->GetYaxis()->SetTitleSize(20);
      hRatio[k][i]->GetYaxis()->SetTitleFont(43);
      hRatio[k][i]->GetYaxis()->SetTitleOffset(1.55);
      hRatio[k][i]->GetYaxis()->SetLabelFont(43); 
      hRatio[k][i]->GetYaxis()->SetLabelSize(15);
      hRatio[k][i]->GetXaxis()->SetTitleSize(20);
      hRatio[k][i]->GetXaxis()->SetTitleFont(43);
      hRatio[k][i]->GetXaxis()->SetTitleOffset(4.);
      hRatio[k][i]->GetXaxis()->SetLabelFont(43); 
      hRatio[k][i]->GetXaxis()->SetLabelSize(15);

      cJetQA[k][i]->SaveAs(Form("plots/QAPlots_pPb_5TeV_%s.pdf",var[i].c_str()),"RECREATE");
    
    }// var loop
   
  }// jet ID loop
    
}
