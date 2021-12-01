{
  TFile *RGB_Data = new TFile("/media/mn688/Elements1/PhD/Analysis_Output/Strangeness_Analysis/Sideband_Subtracted/Sideband_Strangeness_Analysis_RGB_Spring2020_Inbending_dst_Tree_Total_30112021_01.root");
  TFile *RGA_Data = new TFile("/media/mn688/Elements1/PhD/Analysis_Output/Strangeness_Analysis/Sideband_Subtracted/Sideband_Strangeness_Analysis_RGA_Spring2019_Inbending_dst_Tree_Total_30112021_01.root");
  TFile *RGA_Data_Smeared = new TFile("/media/mn688/Elements1/PhD/Analysis_Output/Strangeness_Analysis/Sideband_Subtracted/Sideband_Strangeness_Analysis_RGA_Spring2019_Inbending_dst_Tree_Total_Proton_Smear_30112021_01.root");

  TH1F *RGB_missing_mass = (TH1F*)RGB_Data->Get("h_projectionx_S2_sig_1");
  TH1F *RGA_missing_mass = (TH1F*)RGA_Data->Get("h_projectionx_S2_sig_1");
TH1F *RGA_smeared_missing_mass = (TH1F*)RGA_Data_Smeared->Get("h_projectionx_S2_sig_1");

RGB_missing_mass->Rebin(2);
RGA_missing_mass->Rebin(2);
RGA_smeared_missing_mass->Rebin(2);

Double_t RGA_RGB_Scale = RGA_smeared_missing_mass->Integral() / RGB_missing_mass->Integral();
Double_t RGA_RGA_smeared = RGA_smeared_missing_mass->Integral() / RGA_missing_mass->Integral();

RGB_missing_mass->Scale(RGA_RGB_Scale);
RGA_missing_mass->Scale(RGA_RGA_smeared);

gStyle->SetOptStat(0);
RGB_missing_mass->SetLineColor(kRed);
RGB_missing_mass->SetMarkerColor(kRed);

auto *c1 = new TCanvas("c1","RGB and smeared RGA",800,800);
c1->cd();
RGB_missing_mass->SetTitle("RGB and smeared RGA");
RGB_missing_mass->Draw();
RGA_smeared_missing_mass->Draw("same");

auto *c2 = new TCanvas("c2","RGA and smeared RGA",800,800);
c2->cd();
RGA_smeared_missing_mass->SetTitle("RGA and smeared RGA");
RGA_smeared_missing_mass->Draw();
RGA_missing_mass->Draw("same, hist, L");
}
