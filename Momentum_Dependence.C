// {
//
//   TFile *f1 = new TFile("/media/mn688/Elements1/PhD/Analysis_Output/Hexaquark/RGA_Fall2018_Outbending_at_least_1e1KpFD_Tree_Total_24022022_Total_Scaling_09032022_02.root");
//   // TFile *f1 = new TFile("/media/mn688/Elements1/PhD/Analysis_Output/Hexaquark/RGA_Fall2018_Outbending_at_least_1e1KpFD_Tree_Total_24022022_Total_Scaling_09032022_02.root");
//   TH3F* hist=(TH3F*)f1->Get("h_S1_Kaon_Momentum__Miss_Mass__Kaon_Mass");
//
//   TH1F* Missing_Mass_1 = (TH1F*) hist->ProjectionY("",hist->GetXaxis()->FindBin(0.5),hist->GetXaxis()->FindBin(1.0),0,hist->GetNbinsZ())->Clone("Missing_Mass_1");
//   TH1F* Missing_Mass_2 = (TH1F*) hist->ProjectionY("",hist->GetXaxis()->FindBin(1.0),hist->GetXaxis()->FindBin(1.5),0,hist->GetNbinsZ())->Clone("Missing_Mass_2");
//   TH1F* Missing_Mass_3 = (TH1F*) hist->ProjectionY("",hist->GetXaxis()->FindBin(1.5),hist->GetXaxis()->FindBin(2.0),0,hist->GetNbinsZ())->Clone("Missing_Mass_3");
//   TH1F* Missing_Mass_4 = (TH1F*) hist->ProjectionY("",hist->GetXaxis()->FindBin(2.0),hist->GetXaxis()->FindBin(2.5),0,hist->GetNbinsZ())->Clone("Missing_Mass_4");
//   TH1F* Missing_Mass_5 = (TH1F*) hist->ProjectionY("",hist->GetXaxis()->FindBin(2.5),hist->GetXaxis()->FindBin(3.0),0,hist->GetNbinsZ())->Clone("Missing_Mass_5");
//
//   TH1F* Kaon_Mass_1 = (TH1F*) hist->ProjectionZ("",hist->GetXaxis()->FindBin(0.5),hist->GetXaxis()->FindBin(1.0),0,hist->GetNbinsY())->Clone("Kaon_Mass_1");
//   TH1F* Kaon_Mass_2 = (TH1F*) hist->ProjectionZ("",hist->GetXaxis()->FindBin(1.0),hist->GetXaxis()->FindBin(1.5),0,hist->GetNbinsY())->Clone("Kaon_Mass_2");
//   TH1F* Kaon_Mass_3 = (TH1F*) hist->ProjectionZ("",hist->GetXaxis()->FindBin(1.5),hist->GetXaxis()->FindBin(2.0),0,hist->GetNbinsY())->Clone("Kaon_Mass_3");
//   TH1F* Kaon_Mass_4 = (TH1F*) hist->ProjectionZ("",hist->GetXaxis()->FindBin(2.0),hist->GetXaxis()->FindBin(2.5),0,hist->GetNbinsY())->Clone("Kaon_Mass_4");
//   TH1F* Kaon_Mass_5 = (TH1F*) hist->ProjectionZ("",hist->GetXaxis()->FindBin(2.5),hist->GetXaxis()->FindBin(3.0),0,hist->GetNbinsY())->Clone("Kaon_Mass_5");
//
//   Double_t Missing_Mass_1_integral = Missing_Mass_1->Integral();
//   Double_t Missing_Mass_2_integral = Missing_Mass_2->Integral();
//   Double_t Missing_Mass_3_integral = Missing_Mass_3->Integral();
//   Double_t Missing_Mass_4_integral = Missing_Mass_4->Integral();
//   Double_t Missing_Mass_5_integral = Missing_Mass_5->Integral();
//
//   Double_t Kaon_Mass_1_integral = Kaon_Mass_1->Integral();
//   Double_t Kaon_Mass_2_integral = Kaon_Mass_2->Integral();
//   Double_t Kaon_Mass_3_integral = Kaon_Mass_3->Integral();
//   Double_t Kaon_Mass_4_integral = Kaon_Mass_4->Integral();
//   Double_t Kaon_Mass_5_integral = Kaon_Mass_5->Integral();
//
//   Missing_Mass_2->Scale(Missing_Mass_1_integral / Missing_Mass_2_integral);
//   Missing_Mass_3->Scale(Missing_Mass_1_integral / Missing_Mass_3_integral);
//   Missing_Mass_4->Scale(Missing_Mass_1_integral / Missing_Mass_4_integral);
//   Missing_Mass_5->Scale(Missing_Mass_1_integral / Missing_Mass_5_integral);
//
//
//   Kaon_Mass_2->Scale(Kaon_Mass_1_integral / Kaon_Mass_2_integral);
//   Kaon_Mass_3->Scale(Kaon_Mass_1_integral / Kaon_Mass_3_integral);
//   Kaon_Mass_4->Scale(Kaon_Mass_1_integral / Kaon_Mass_4_integral);
//   Kaon_Mass_5->Scale(Kaon_Mass_1_integral / Kaon_Mass_5_integral);
//
//
//   Missing_Mass_2->SetLineColor(kRed);
//   Missing_Mass_3->SetLineColor(kOrange);
//   Missing_Mass_4->SetLineColor(kBlack);
//   Missing_Mass_5->SetLineColor(kGreen);
//
//   Kaon_Mass_2->SetLineColor(kRed);
//   Kaon_Mass_3->SetLineColor(kOrange);
//   Kaon_Mass_4->SetLineColor(kBlack);
//   Kaon_Mass_5->SetLineColor(kGreen);
//
//
//   TF1 *func1 = new TF1("func1","gaus(0)+gaus(3)+pol3(6)",0.363,0.6);
//   TF1 *func2 = new TF1("func2","gaus(0)+gaus(3)",0.363,0.6);
//   TF1 *func3 = new TF1("func3","pol3(0)",0.363,0.6);
//   func1->SetParameter(0,Kaon_Mass_3->GetMaximum()/2);
//   func1->SetParameter(1,0.493);
//   func1->SetParameter(2,0.03);
//   func1->SetParameter(3,Kaon_Mass_3->GetMaximum()/3);
//   func1->SetParameter(4,0.493);
//   func1->SetParameter(5,0.05);
//   // func1->SetParameter(5,0.05);
//   // func1->SetParameter(5,0.05);
//   // func1->SetParameter(5,0.05);
//   // func1->SetParameter(5,0.05);
//
//   Kaon_Mass_3->Fit("func1","RB");
//
//   func2->FixParameter(0,func1->GetParameter(0));
//   func2->FixParameter(1,func1->GetParameter(1));
//   func2->FixParameter(2,func1->GetParameter(2));
//   func2->FixParameter(3,func1->GetParameter(3));
//   func2->FixParameter(4,func1->GetParameter(4));
//   func2->FixParameter(5,func1->GetParameter(5));
//   func3->FixParameter(0,func1->GetParameter(6));
//   func3->FixParameter(1,func1->GetParameter(7));
//   func3->FixParameter(2,func1->GetParameter(8));
//   func3->FixParameter(3,func1->GetParameter(9));
//   // func3->FixParameter(4,func1->GetParameter(10));
//
//   func1->SetLineColor(kBlue);
//   func2->SetLineColor(kBlack);
//   func3->SetLineColor(kGreen);
//
//
//
//   auto *c1 = new TCanvas("c1","",800,800);
//   c1->cd();
//   Missing_Mass_1->Draw("same");
//   Missing_Mass_2->Draw("same,hist");
//   Missing_Mass_3->Draw("same,hist");
//   Missing_Mass_4->Draw("same,hist");
//   Missing_Mass_5->Draw("same,hist");
//
//   auto *c2=new TCanvas("c2","",800,800);
//   c2->cd();
//   Kaon_Mass_1->Draw("same");
//   Kaon_Mass_2->Draw("same,hist");
//   Kaon_Mass_3->Draw("same,hist");
//   Kaon_Mass_4->Draw("same,hist");
//   Kaon_Mass_5->Draw("same,hist");
//
//   auto *c3 = new TCanvas("c3","",800,800);
//   c3->cd();
//   Kaon_Mass_3->Draw("same,hist");
//   func1->Draw("same");
//   func2->Draw("same");
//   func3->Draw("same");
//
// }

// Making lots of projections in momentum

{
  TFile *f1 = new TFile("/media/mn688/Elements1/PhD/Analysis_Output/Hexaquark/RGA_Fall2018_Outbending_at_least_1e1KpFD_Tree_Total_24022022_Total_Scaling_10032022_01.root");
  // TFile *f1 = new TFile("/media/mn688/Elements1/PhD/Analysis_Output/Hexaquark/RGA_Fall2018_Outbending_at_least_1e1KpFD_Tree_Total_216022022_Total_Scaling_09032022_02.root");
  TH3F* hist2=(TH3F*)f1->Get("h_S1_Kaon_Momentum__Miss_Mass__Kaon_Mass");

  TH1F *h_sigma_1 = new TH1F("h_sigma_1","sigma 1 relationship",160,0,160);
  TH1F *h_sigma_2 = new TH1F("h_sigma_2","sigma 2 relationship",160,0,160);
  TH1F *h_chi2 = new TH1F("h_chi2","chi^{2} relationship",160,0,160);

  // Changing loop to momentum bin values
  Double_t bin_1, bin_2;
  // Function range
  Double_t Range_min, Range_max;

  // TF1 *func1 = new TF1("func1","gaus(0) + pol3(3)",0.363,0.6);
  TH1F *Missing_Mass_Projections[160];
  TH1F *Kaon_Mass_Projections[160];

  // Parameter values to be set
  Double_t par0, par1, par2, par3, par4, par5, par6, par7, par8, par9;
  for(Int_t i = 0; i < 160; i++){

     bin_1 = 1 + (i*0.01);
     bin_2 = 1 + (i+1)*0.01;
     if(bin_1 > 1.1 && bin_1 < 1.35) continue;
     if(bin_1 > 3) continue;
     cout<<"bin "<<i<<" 1st bin "<<bin_1<<" 2nd bin "<<bin_2<<endl;

     Missing_Mass_Projections[i] = (TH1F*) hist2->ProjectionY("",hist2->GetXaxis()->FindBin(bin_1),hist2->GetXaxis()->FindBin(bin_2),0,hist2->GetNbinsZ())->Clone();
     Kaon_Mass_Projections[i] = (TH1F*) hist2->ProjectionZ("",hist2->GetXaxis()->FindBin(bin_1),hist2->GetXaxis()->FindBin(bin_2),0,hist2->GetNbinsY())->Clone();

     // Setting range of functions


     TF1 *func1 = new TF1("func1","([0] * exp(-pow(x-[1],2) / (2 * pow([2],2)))) + gaus(3) + pol1(6) + ([8] * exp(-pow(x-[1],2) / (2 * pow([9],2))))",0.363,0.6);
     TF1 *signal_1 = new TF1("signal_1","gaus",0.363,0.6);
     TF1 *signal_2 = new TF1("signal_2","gaus",0.363,0.6);
     TF1 *background_1 = new TF1("background_1","gaus",0.363,0.6);
     TF1 *background_2 = new TF1("background_2","pol1",0.363,0.6);


     // Setting parameters based on functions determined
     par0 = Kaon_Mass_Projections[i]->GetMaximum() / 2; // amplitude for 1st signal gaus
     par1 = 0.493; // mean for both signal gauss
     par2 = 0.0000005182 * i * i + 0.00009958 * i + 0.007042; // sigma for 1st signal gaus
     par4 = 0.365; // mean for back gaus
     par5 = -0.000006631 * i * i + 0.0005252 * i + 0.09216; // sigma for back gaus
     par6 = 20; // constant for back pol1
     par7 = -10; // x component for back pol1
     par8 = Kaon_Mass_Projections[i]->GetMaximum() / 3; // amplitude for 2nd signal gaus
     par9 = -0.000000026978 * i * i + 0.0003172 * i + 0.01033; // sigma for 2nd signal gaus


    // Setting parameters
    func1->SetParameter(0,par0); // amplitude for 1st signal gaus
    func1->SetParameter(1,par1); // mean for both signal gaus
    func1->SetParameter(2,par2); // sigma for 1st signal gaus
    func1->FixParameter(4,par4); // mean for back gaus
    func1->SetParameter(5,par5); // sigma for back gaus
    func1->SetParameter(6,par6); // constant for back pol1
    func1->SetParameter(7,par7); // slope for back pol1
    func1->SetParameter(8,par8); // amplitude for 2nd signal gaus
    func1->SetParameter(9,par9); // sigma for 2nd signal gaus

    // Setting parameter limits
    func1->SetParLimits(0,Kaon_Mass_Projections[i]->GetMaximum() / 20, Kaon_Mass_Projections[i]->GetMaximum()); // amplitude for 1st signal gaus
    func1->SetParLimits(1,par1 * 0.99, par1 * 1.01); // mean for both signal gaus
    func1->SetParLimits(2,par2 * 0.95, par2 * 1.05); // sigma for 1st signal gaus
    func1->SetParLimits(3,0,1.0 * Kaon_Mass_Projections[i]->GetBinContent(Kaon_Mass_Projections[i]->FindBin(0.365))); // amplitude for back gaus
    // func1->SetParLimits(4,par4 * 0.98, par4 * 1.02); // mean for back gaus
    func1->SetParLimits(5,par5 * 0.95, par5 * 1.05); // sigma for back gaus
    func1->SetParLimits(6,4, Kaon_Mass_Projections[i]->GetMaximum()); // constant for back pol1
    func1->SetParLimits(7,-1000, 0); // slope for back pol1
    func1->SetParLimits(8,Kaon_Mass_Projections[i]->GetMaximum() / 20, Kaon_Mass_Projections[i]->GetMaximum()); // amplitude for 2nd signal gaus
    func1->SetParLimits(9,par9 * 0.95, par9 * 1.05); // sigma for 2nd signal gaus


    Kaon_Mass_Projections[i]->Fit("func1","RBQ");

    signal_1->SetParameter(0,func1->GetParameter(0));
    signal_1->SetParameter(1,func1->GetParameter(1));
    signal_1->SetParameter(2,func1->GetParameter(2));
    signal_2->SetParameter(0,func1->GetParameter(8));
    signal_2->SetParameter(1,func1->GetParameter(1));
    signal_2->SetParameter(2,func1->GetParameter(9));
    background_1->SetParameter(0,func1->GetParameter(3));
    background_1->SetParameter(1,func1->GetParameter(4));
    background_1->SetParameter(2,func1->GetParameter(5));
    background_2->SetParameter(0,func1->GetParameter(6));
    background_2->SetParameter(1,func1->GetParameter(7));

    cout<<" chi "<<func1->GetChisquare() / func1->GetNDF()<<endl;
    // signal_1->Draw("same");
    // signal_2->Draw("same");
    // background_1->Draw("same");
    // background_2->Draw("same");

    // func1->SetParameter(0,Kaon_Mass_Projections[i]->GetMaximum()/2);
    // func1->SetParameter(1,0.493);
    // func1->SetParameter(2,0.000018756*pow(i,2) + 0.0005938 * i + 0.0069734);
    // func1->SetParameter(3,Kaon_Mass_Projections[i]->GetMaximum()/3);
    // func1->SetParameter(4,0.493);
    // func1->SetParameter(5,0.0018025*pow(i,2) + 0.001222 * i + 0.0048292);
    //
    // func1->SetParLimits(2,0.6*(0.000018756*pow(i,2) + 0.0005938 * i + 0.0069734),1.3*(0.000018756*pow(i,2) + 0.0005938 * i + 0.0069734));
    // func1->SetParLimits(5,0.6*(0.0018025*pow(i,2) + 0.001222 * i + 0.0048292),1.3*(0.0018025*pow(i,2) + 0.001222 * i + 0.0048292));



    // Kaon_Mass_Projections[i]->Fit("func1","RBQ");

    // // cout<<"sigma "<<func1->GetParameter(2)<<endl;
    h_chi2->Fill(i,func1->GetChisquare() / func1->GetNDF());
    h_sigma_1->Fill(i,func1->GetParameter(2));
    h_sigma_1->SetBinError(i,func1->GetParError(2));
    h_sigma_2->Fill(i,func1->GetParameter(9));
    h_sigma_2->SetBinError(i,func1->GetParError(9));

    cout<<func1->GetParameter(2)/func1->GetParameter(8)<<endl;
  }

  // TF1 *func2 = new TF1("func2","[0]*x^2 + [1]*x + [2]",1,22);
  TF1 *func2 = new TF1("func2","pol2(0)",1,160);
  TF1 *func3 = new TF1("func3","pol2(0)",1,160);

  // h_sigma_1->Fit("func2","RB");
  // h_sigma_2->Fit("func3","RB");
  // cout<<func2->GetParameter(0)<<endl;
  // h_sigma_2->Draw();
  // h_sigma_2->Draw();

  // signal_1->SetLineColor(kBlue);
  // signal_2->SetLineColor(kOrange);
  // background_1->SetLineColor(kBlack);
  // background_2->SetLineColor(kGreen);

  Kaon_Mass_Projections[99]->Draw();
  signal_1->Draw("same");
  signal_2->Draw("same");
  background_1->Draw("same");
  background_2->Draw("same");
  // func1->Draw("same");

}
