{
  TFile *RGB_Data = new TFile("/media/mn688/Elements1/PhD/Analysis_Output/Strangeness_Analysis/Sideband_Subtracted/Sideband_Strangeness_Analysis_RGB_Spring2020_Inbending_dst_Tree_Total_30112021_01.root");
  TFile *RGA_Data = new TFile("/media/mn688/Elements1/PhD/Analysis_Output/Strangeness_Analysis/Sideband_Subtracted/Sideband_Strangeness_Analysis_RGA_Spring2019_Inbending_dst_Tree_Total_30112021_01.root");
  TFile *RGA_Data_Smeared = new TFile("/media/mn688/Elements1/PhD/Analysis_Output/Strangeness_Analysis/Sideband_Subtracted/Sideband_Strangeness_Analysis_RGA_Spring2019_Inbending_dst_Tree_Total_Proton_Smear_30112021_01.root");

  TH1F *RGB_missing_mass = (TH1F*)RGB_Data->Get("h_projectionx_S2_sig_1");
  TH1F *RGA_missing_mass = (TH1F*)RGA_Data->Get("h_projectionx_S2_sig_1");
TH1F *RGA_smeared_missing_mass = (TH1F*)RGA_Data_Smeared->Get("h_projectionx_S2_sig_1");

RGB_missing_mass->Draw();
RGA_missing_mass->Draw("same");
RGA_smeared_missing_mass->Draw("same");
}
