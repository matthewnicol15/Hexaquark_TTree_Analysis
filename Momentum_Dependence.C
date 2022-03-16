// Making lots of projections in momentum

{
  TFile *f1 = new TFile("/media/mn688/Elements1/PhD/Analysis_Output/Hexaquark/RGA_Fall2018_Outbending_at_least_1e1KpFD_Tree_Total_24022022_Total_Scaling_16032022_01.root");
  // TFile *f1 = new TFile("/media/mn688/Elements1/PhD/Analysis_Output/Hexaquark/RGA_Fall2018_Outbending_at_least_1e1KpFD_Tree_Total_216022022_Total_Scaling_09032022_02.root");
  TH3F* hist2=(TH3F*)f1->Get("h_S3_Kaon_Momentum__Miss_Mass__Kaon_Mass");

  TH1F *h_sigma_1 = new TH1F("h_sigma_1","sigma 1 relationship",160,1,161);
  TH1F *h_sigma_2 = new TH1F("h_sigma_2","sigma 2 relationship",160,1,161);
  TH1F *h_sigma_3 = new TH1F("h_sigma_3","sigma 3 relationship",160,1,161);
  TH1F *h_chi2 = new TH1F("h_chi2","chi^{2} relationship",160,1,161);

  // Changing loop to momentum bin values
  Double_t bin_1, bin_2;
  // Reduced chi^2
  Double_t Chi2;
  // Function range
  Double_t Range_min, Range_max;

  // TF1 *func1 = new TF1(funcname.str().c_str(),"gaus(0) + pol3(3)",0.365,0.6);
  TH1F *Missing_Mass_Projections[160];
  TH1F *Kaon_Mass_Projections[160];

  TF1 *func1[160];
  TF1 *signal_1[160];
  TF1 *signal_2[160];
  TF1 *background_1[160];
  TF1 *background_2[160];

  Int_t counter = 0;

  auto *c1 = new TCanvas("c1","canvas",800,800);
  // Parameter values to be set
  Double_t par0, par1, par2, par3, par4, par5, par6, par7, par8, par9;
  for(Int_t i = 0; i < 160; i++){
     counter++;
     // Determining bin max and min
     bin_1 = 1 + (i*0.01);
     bin_2 = 1 + (i+1)*0.01;


     // Cutting out pion cross over between 1.1 and 1.35 GeV
     if(bin_1 > 1.1 && bin_1 < 1.35) continue;

     // Kaons only plotted up to 3 GeV
     if(bin_1 > 3) continue;
     // cout<<"bin "<<i<<" 1st bin "<<bin_1<<" 2nd bin "<<bin_2<<endl;

     Missing_Mass_Projections[i] = (TH1F*) hist2->ProjectionY("",hist2->GetXaxis()->FindBin(bin_1),hist2->GetXaxis()->FindBin(bin_2),0,hist2->GetNbinsZ())->Clone();
     Kaon_Mass_Projections[i] = (TH1F*) hist2->ProjectionZ("",hist2->GetXaxis()->FindBin(bin_1),hist2->GetXaxis()->FindBin(bin_2),0,hist2->GetNbinsY())->Clone();

     Kaon_Mass_Projections[i]->Rebin(2);

     // Setting range of functions
     ostringstream funcname, sig1, sig2, back1, back2;
     funcname <<"func_"<<i;
     sig1 <<"sig1_"<<i;
     sig2 <<"sig2_"<<i;
     back1 <<"back1_"<<i;
     back2 <<"back2_"<<i;
     func1[i] = new TF1(funcname.str().c_str(),"[0] * exp(-pow(x-[1],2) / (2 * pow([2],2))) + [3] * exp(-pow(x-[1],2) / (2 * pow([2]*[4],2))) + gaus(5) + [8]*[9] + [9]*x",0.365,0.65);
     signal_1[i] = new TF1(sig1.str().c_str(),"gaus",0.365,0.65);
     signal_2[i] = new TF1(sig2.str().c_str(),"gaus",0.365,0.65);
     background_1[i] = new TF1(back1.str().c_str(),"gaus",0.0,0.65);
     background_2[i] = new TF1(back2.str().c_str(),"pol1",0.365,0.65);



     // Setting parameters based on functions determined
     par0 = Kaon_Mass_Projections[i]->GetMaximum() / 2; // amplitude for 1st signal gaus
     par1 = 0.493; // mean for both signal gauss
     par2 = 0.000000857 * pow(i+1,2) + 0.00001848 * (i+1) + 0.01094; // sigma for 1st signal gaus
     par3 = Kaon_Mass_Projections[i]->GetMaximum() / 3; // amplitude for 2nd signal gaus
     par4 = -0.00001617 * pow(i+1,2) + 0.001323 * (i+1) + 1.928; // ratio between signal gaussian sigmas
     par5 = 2 * Kaon_Mass_Projections[i]->GetBinContent(Kaon_Mass_Projections[i]->FindBin(0.365)); // amplitude for back gaus
     par6 = 0.1396; // mean for back gaus (pion mass)
     par7 = 0.3; // sigma for back gaus
     par8 = -5; // constant for back pol1
     par9 = -10; // slope for back pol1

    // Setting parameters
    func1[i]->SetParameter(0,par0); // amplitude for 1st signal gaus
    func1[i]->SetParameter(1,par1); // mean for both signal gaus
    func1[i]->SetParameter(2,par2); // sigma for 1st signal gaus
    func1[i]->SetParameter(3,par3); // amplitude for 2nd signal gaus
    func1[i]->FixParameter(4,par4); // ratio between signal gaussian sigmas
    func1[i]->SetParameter(5,par5); // amplitude for back gaus
    func1[i]->SetParameter(6,par6); // mean for back gaus (pion mass)
    func1[i]->SetParameter(7,par7); // sigma for back gaus
    func1[i]->SetParameter(8,par8); // constant for back pol1
    func1[i]->SetParameter(9,par9); // slope for back pol1

    // Setting parameter limits
    func1[i]->SetParLimits(0,Kaon_Mass_Projections[i]->GetMaximum() / 20, Kaon_Mass_Projections[i]->GetMaximum()); // amplitude for 1st signal gaus
    func1[i]->SetParLimits(1,par1 * 0.98, par1 * 1.02); // mean for both signal gaus
    func1[i]->SetParLimits(2,0.85 * par2, 1.15 * par2); // sigma for 1st signal gaus
    func1[i]->SetParLimits(3,Kaon_Mass_Projections[i]->GetMaximum() / 20, Kaon_Mass_Projections[i]->GetMaximum()); // amplitude for 2nd signal gaus
    func1[i]->SetParLimits(4,0.85 * par4, 1.15 * par4); // ratio between signal gaussian sigmas
    func1[i]->SetParLimits(5,Kaon_Mass_Projections[i]->GetBinContent(Kaon_Mass_Projections[i]->FindBin(0.365)) / 5, 4 * Kaon_Mass_Projections[i]->GetMaximum()); // amplitude for back gaus
    func1[i]->SetParLimits(6,par6 * 0.98, par6 * 1.02); // mean for back gaus (pion mass)
    func1[i]->SetParLimits(7,0.1, 0.8); // sigma for back gaus
    func1[i]->SetParLimits(8,-1000, -0.7); // constant for back pol1
    func1[i]->SetParLimits(9,-500, 0); // slope for back pol1


    Kaon_Mass_Projections[i]->Fit(funcname.str().c_str(),"RBQ");

    Chi2 = func1[i]->GetChisquare() / func1[i]->GetNDF();


    signal_1[i]->SetParameter(0,func1[i]->GetParameter(0));
    signal_1[i]->SetParameter(1,func1[i]->GetParameter(1));
    signal_1[i]->SetParameter(2,func1[i]->GetParameter(2));
    signal_2[i]->SetParameter(0,func1[i]->GetParameter(3));
    signal_2[i]->SetParameter(1,func1[i]->GetParameter(1));
    signal_2[i]->SetParameter(2,func1[i]->GetParameter(2) * func1[i]->GetParameter(4));
    background_1[i]->SetParameter(0,func1[i]->GetParameter(5));
    background_1[i]->SetParameter(1,func1[i]->GetParameter(6));
    background_1[i]->SetParameter(2,func1[i]->GetParameter(7));
    background_2[i]->SetParameter(0,func1[i]->GetParameter(8) * func1[i]->GetParameter(9));
    background_2[i]->SetParameter(1,func1[i]->GetParameter(9));


    h_chi2->Fill(i+1,func1[i]->GetChisquare() / func1[i]->GetNDF());
    h_sigma_1->Fill(i+1,func1[i]->GetParameter(2));
    h_sigma_1->SetBinError(i+1,func1[i]->GetParError(2));
    h_sigma_2->Fill(i+1,func1[i]->GetParameter(4));
    h_sigma_2->SetBinError(i+1,func1[i]->GetParError(4));
    h_sigma_3->Fill(i+1,func1[i]->GetParameter(7));
    h_sigma_3->SetBinError(i+1,func1[i]->GetParError(7));

  }
  TF1 *func2 = new TF1("func2","pol2(0)",38,161);
  TF1 *func3 = new TF1("func3","pol2(0)",38,141);
  h_sigma_1->Fit("func2","RB");
  h_sigma_2->Fit("func3","RB");


}
