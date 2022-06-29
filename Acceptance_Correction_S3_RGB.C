{

   TFile *f1 = new TFile("/media/mn688/Elements1/PhD/Analysis_Output/Hexaquark/S3_Sim_MC_RGB_Spring2019_Inbending_at_least_1e1KpFD_KpChi3_Tree_Total_23062022_Total_Scaling_23062022_01.root");
   TFile *f2 = new TFile("/media/mn688/Elements1/PhD/Analysis_Output/Hexaquark/S3_Sim_RGB_Spring2019_Inbending_at_least_1e1KpFD_KpChi3_Tree_Total_23062022_Total_Scaling_23062022_01.root");
   TFile *f3 = new TFile("/media/mn688/Elements1/PhD/Analysis_Output/Hexaquark/RGB_Spring2019_Inbending_at_least_1e1KpFD_Tree_Total_03032022_Total_Scaling_23032022_01.root");
   TFile fileOutput1("/media/mn688/Elements1/PhD/Analysis_Output/Hexaquark/test.root","recreate");


   TH3F *h_S3_Kaon_Momentum__Miss_Mass__Kaon_Mass_1 = (TH3F*)f1->Get("h_S3_Kaon_Momentum__Miss_Mass__Kaon_Mass");
   TH3F *h_S3_Kaon_Momentum__Miss_Mass__Kaon_Mass_2 = (TH3F*)f2->Get("h_S3_Kaon_Momentum__Miss_Mass__Kaon_Mass");
   TH3F *h_S3_Kaon_Momentum__Miss_Mass__Kaon_Mass_3 = (TH3F*)f3->Get("h_S3_Kaon_Momentum__Miss_Mass__Kaon_Mass");

   TH1F *MC = (TH1F*) h_S3_Kaon_Momentum__Miss_Mass__Kaon_Mass_1->ProjectionY()->Clone("MC");
   TH1F *OSG = (TH1F*) h_S3_Kaon_Momentum__Miss_Mass__Kaon_Mass_2->ProjectionY()->Clone("OSG");
   TH1F *EXP = (TH1F*) h_S3_Kaon_Momentum__Miss_Mass__Kaon_Mass_3->ProjectionY()->Clone("EXP");
   TH1F *EXP2 = new TH1F("EXP2","copy",300,2,5);

   MC->Divide(OSG);

   for(Int_t i = 1; i < EXP->GetNbinsX() + 1; i ++ ){

      EXP2->SetBinContent(i,EXP->GetBinContent(EXP->FindBin(EXP2->GetBinCenter(i) - 0.940)));
      // if(OSG->GetBinContent(i)!=0)
      EXP2->SetBinContent(i,EXP2->GetBinContent(i)*MC->GetBinContent(i));
   }

EXP2->Draw();

   // auto *S3_Threshold_line=new TLine(2.752,0,2.752,0.75);
   // S3_Threshold_line->SetLineColor(kRed);
   // // OSG->GetXaxis()->SetRange(OSG->FindBin(2.752),OSG->FindBin(5));
   // OSG->Draw();
   // S3_Threshold_line->Draw("same");

   // fileOutput1.Write();
}
