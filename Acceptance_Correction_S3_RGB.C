{
   TFile *f1 = new TFile("/media/mn688/Elements1/PhD/Analysis_Output/Hexaquark/S3_CascadeSigma_Acceptance_Sim_job4936_MC_RGB_Spring2019_Inbending_at_least_1e1KpFD_KpChi3_Tree_Total_29062022_Total_Scaling_01072022_01.root");
   TFile *f2 = new TFile("/media/mn688/Elements1/PhD/Analysis_Output/Hexaquark/S3_CascadeSigma_Acceptance_Sim_job4936_RGB_Spring2019_Inbending_at_least_1e1KpFD_KpChi3_Tree_Total_29062022_Total_Scaling_01072022_01.root");
   // TFile *f3 = new TFile("/media/mn688/Elements1/PhD/Analysis_Output/Hexaquark/RGB_Spring2019_Inbending_at_least_1e1KpFD_Tree_Total_03032022_Total_Scaling_23032022_01.root");
   TFile *f3 = new TFile("/media/mn688/Elements1/PhD/Macros/Hexaquark_TTree_Analysis/Kaon_Background_Subtraction_RGB_Spring2019_Inbending_01072022_02.root");

   // TFile fileOutput1("/media/mn688/Elements1/PhD/Analysis_Output/Hexaquark/test.root","recreate");


   TH3F *h_S3_Kaon_Momentum__Miss_Mass__Kaon_Mass_1 = (TH3F*)f1->Get("h_S3_Kaon_Momentum__Miss_Mass__Kaon_Mass");
   TH3F *h_S3_Kaon_Momentum__Miss_Mass__Kaon_Mass_2 = (TH3F*)f2->Get("h_S3_Kaon_Momentum__Miss_Mass__Kaon_Mass");
   TH1F *EXP = (TH1F*)f3->Get("S3Result_Adjusted");

   TH1F *MC = (TH1F*) h_S3_Kaon_Momentum__Miss_Mass__Kaon_Mass_1->ProjectionY("",1,h_S3_Kaon_Momentum__Miss_Mass__Kaon_Mass_1->GetNbinsX(),1,h_S3_Kaon_Momentum__Miss_Mass__Kaon_Mass_1->GetNbinsZ(),"")->Clone("MC");
   TH1F *OSG = (TH1F*) h_S3_Kaon_Momentum__Miss_Mass__Kaon_Mass_2->ProjectionY("",1,h_S3_Kaon_Momentum__Miss_Mass__Kaon_Mass_2->GetNbinsX(),1,h_S3_Kaon_Momentum__Miss_Mass__Kaon_Mass_2->GetNbinsZ(),"")->Clone("OSG");
   // TH1F *EXP = (TH1F*) h_S3_Kaon_Momentum__Miss_Mass__Kaon_Mass_3->ProjectionY()->Clone("EXP");

   Int_t Rebin = MC->GetNbinsX() / EXP->GetNbinsX();
   MC->Rebin(Rebin);
   OSG->Rebin(Rebin);

   MC->Divide(OSG);

   EXP->Multiply(MC);

   // for(Int_t i = 1; i < EXP->GetNbinsX() + 1; i ++ ){
   //
   //    EXP->SetBinContent(i,EXP->GetBinContent(i)*MC->GetBinContent(i));
   // }

   EXP->SetTitle("Strangeness 3 background subtracted and acceptance corrected;MM(e' K^{+} K^{+} K^{+}) [GeV]");
   EXP->Draw("e");

   auto *S3_Threshold_Line=new TLine(2.503,0,2.503,2500);
   auto *Hex=new TLine(2.677,0,2.677,2500);
   auto *Mol1=new TLine(2.796,0,2.796,2500);
   auto *Mol2=new TLine(3.022,0,3.022,2500);
   S3_Threshold_Line->SetLineColor(kOrange);
   Hex->SetLineColor(kRed);
   Mol1->SetLineColor(kGreen);
   Mol2->SetLineColor(kGreen);
   S3_Threshold_Line->SetLineStyle(2);
   Hex->SetLineStyle(2);
   Mol1->SetLineStyle(2);
   Mol2->SetLineStyle(2);
   S3_Threshold_Line->SetLineWidth(2);
   Hex->SetLineWidth(2);
   Mol1->SetLineWidth(2);
   Mol2->SetLineWidth(2);
   S3_Threshold_Line->Draw("same");
   Hex->Draw("same");
   Mol1->Draw("same");
   Mol2->Draw("same");

   // fileOutput1.Write();
}
