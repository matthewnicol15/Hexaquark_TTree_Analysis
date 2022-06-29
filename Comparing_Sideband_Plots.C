{
  TFile *RGB_Data = new TFile("/media/mn688/Elements1/PhD/Analysis_Output/Strangeness_Analysis/Sideband_Subtracted/Sideband_Strangeness_Analysis_RGB_Spring2020_Inbending_dst_Tree_Total_11012021_02.root");
  TFile *RGA_Data = new TFile("/media/mn688/Elements1/PhD/Analysis_Output/Strangeness_Analysis/Sideband_Subtracted/Sideband_Strangeness_Analysis_RGA_Spring2019_Inbending_dst_Tree_Total_11012021_02.root");
  TFile *RGA_Data_Smeared = new TFile("/media/mn688/Elements1/PhD/Analysis_Output/Strangeness_Analysis/Sideband_Subtracted/Sideband_Strangeness_Analysis_RGA_Spring2019_Inbending_dst_Tree_Total__Proton_Smeared_11012021_02.root");

  // Strangeness 1
  TH1F *RGB_S1_missing_mass = (TH1F*)RGB_Data->Get("h_projectionx_S1_sig_1");
  TH1F *RGA_S1_missing_mass = (TH1F*)RGA_Data->Get("h_projectionx_S1_sig_1");
  TH1F *RGA_S1_smeared_missing_mass = (TH1F*)RGA_Data_Smeared->Get("h_projectionx_S1_sig_1");

  // Strangeness 2
  TH1F *RGB_S2_missing_mass = (TH1F*)RGB_Data->Get("h_projectionx_S2_sig_1");
  TH1F *RGA_S2_missing_mass = (TH1F*)RGA_Data->Get("h_projectionx_S2_sig_1");
  TH1F *RGA_S2_smeared_missing_mass = (TH1F*)RGA_Data_Smeared->Get("h_projectionx_S2_sig_1");

  // RGB_S1_missing_mass->Rebin(2);
  // RGA_S1_missing_mass->Rebin(2);
  // RGA_S1_smeared_missing_mass->Rebin(2);

  RGB_S2_missing_mass->Rebin(2);
  RGA_S2_missing_mass->Rebin(2);
  RGA_S2_smeared_missing_mass->Rebin(2);

  Double_t RGA_RGB_S1_Scale = RGA_S1_smeared_missing_mass->Integral() / RGB_S1_missing_mass->Integral();
  Double_t RGA_RGA_S1_smeared = RGA_S1_smeared_missing_mass->Integral() / RGA_S1_missing_mass->Integral();

  Double_t RGA_RGB_S2_Scale = RGA_S2_smeared_missing_mass->Integral() / RGB_S2_missing_mass->Integral();
  Double_t RGA_RGA_S2_smeared = RGA_S2_smeared_missing_mass->Integral() / RGA_S2_missing_mass->Integral();

  RGB_S1_missing_mass->Scale(RGA_RGB_S1_Scale);
  RGA_S1_missing_mass->Scale(RGA_RGA_S1_smeared);


  RGB_S2_missing_mass->Scale(RGA_RGB_S2_Scale);
  RGA_S2_missing_mass->Scale(RGA_RGA_S2_smeared);

  gStyle->SetOptStat(0);
  RGB_S1_missing_mass->SetLineColor(kRed);
  RGB_S1_missing_mass->SetMarkerColor(kRed);
  RGB_S2_missing_mass->SetLineColor(kRed);
  RGB_S2_missing_mass->SetMarkerColor(kRed);

  auto *c1 = new TCanvas("c1","RGB and smeared RGA",800,800);
  c1->cd();
  RGB_S1_missing_mass->SetTitle("RGB Strageness 1");
  RGB_S1_missing_mass->Draw("c");
  // RGA_S1_smeared_missing_mass->Draw("hist,same");

  auto *c2 = new TCanvas("c2","RGA and smeared RGA",800,800);
  c2->cd();
  RGA_S1_missing_mass->SetTitle("RGA Strageness 1");
  RGA_S1_missing_mass->Draw("c");
  // RGA_S1_smeared_missing_mass->SetTitle("RGA and smeared RGA");
  // RGA_S1_smeared_missing_mass->Draw();

  auto *c3 = new TCanvas("c3","RGB and smeared RGA",800,800);
  c3->cd();
  RGB_S2_missing_mass->SetTitle("RGB Strageness 2");
  RGB_S2_missing_mass->Draw("c");
  // RGA_S2_smeared_missing_mass->Draw("same");

  auto *c4 = new TCanvas("c4","RGA and smeared RGA",800,800);
  c4->cd();
  // RGA_S2_smeared_missing_mass->SetTitle("RGA and smeared RGA");
  // RGA_S2_smeared_missing_mass->Draw();
  RGA_S2_missing_mass->SetTitle("RGA Strangeness 2");
  RGA_S2_missing_mass->Draw("c");

}
