// Making lots of projections in momentum

{


   // Get input file
   // TFile *f1 = new TFile("/media/mn688/Elements1/PhD/Analysis_Output/Hexaquark/RGA_Fall2018_Outbending_at_least_1e1KpFD_Tree_Total_24022022_Total_Scaling_16032022_01.root");
   TFile *f1 = new TFile("/media/mn688/Elements1/PhD/Analysis_Output/Hexaquark/RGB_Spring2019_Inbending_at_least_1e1KpFD_Tree_Total_03032022_Total_Scaling_23032022_01.root");
   // TFile *f1=new TFile("/media/mn688/Elements1/PhD/Analysis_Output/Hexaquark/RGB_Spring2019_Inbending_at_least_1e1KpFD_Tree_Total_03032022_Total_Scaling_23032022_01.root");
   // Get histogram from input file
   TH3F* hist_S1=(TH3F*)f1->Get("h_S1_Kaon_Momentum__Miss_Mass__Kaon_Mass");
   TH3F* hist_S2=(TH3F*)f1->Get("h_S2_Kaon_Momentum__Miss_Mass__Kaon_Mass");

   // Rebin strangeness 2 due to low statistics
   hist_S2->Rebin3D(2,1,2);

   TFile fileOutput1("Kaon_Background_Subtraction_RGA_Fall2018_Inbending_Alternative_19052022_02.root","recreate");

   // histograms looking at fit parameters as a function of momentum
   TH1F *h_sigma_1 = new TH1F("h_sigma_1","sigma 1 relationship",300,1,2.6);
   TH1F *h_mean_signal = new TH1F("h_mean_signal","signal mean relationship",300,1,2.6);
   TH1F *h_amp_1 = new TH1F("h_amp_1","amplitude 1 relationship",300,1,2.6);
   TH1F *h_amp_1_S2 = new TH1F("h_amp_1_S2","amplitude 1 relationship",300,1,2.6);
   TH1F *h_sigma_factor = new TH1F("h_sigma_factor","sigma 2 relationship to sigma 1",300,1,2.6);
   TH1F *h_amp_2 = new TH1F("h_amp_2","amplitude 2 relationship",300,1,2.6);
   TH1F *h_amp_2_S2 = new TH1F("h_amp_2_S2","amplitude 2 relationship",300,1,2.6);
   TH1F *h_sigma_3 = new TH1F("h_sigma_3","sigma 3 relationship",300,1,2.6);
   TH1F *h_chi2 = new TH1F("h_chi2","chi^{2} relationship",300,1,2.6);

   TH1F *h_S1_Integral = new TH1F("h_S1_Integral","S1 integral",300,1,2.6);
   TH1F *h_S2_Integral = new TH1F("h_S2_Integral","S2 integral",300,1,2.6);
   // TH1F *h_S1_S2_Ratio = new TH1F("h_S1_S2_Ratio","Ratio between S1 and S2",300,1,2.6);


   // Changing loop to momentum bin values
   Double_t bin_1[300], bin_2[300], momentum_S1[300], momentum_S2[300], momentum_mid_S1, momentum_mid_S2;
   // Ratio between S1 and S2
   Double_t S1_Integral[300], S2_Integral[300];
   // Reduced chi^2
   Double_t Chi2;
   // Parameter values
   Double_t Sigma_1[300], Amp_1[300], Amp_2[300], amp_ratio[300];
   Double_t Sigma_1_S2[300], Amp_1_S2[300], Amp_2_S2[300], amp_ratio_S2[300];
   // Parameter errors
   Double_t Sigma_1_error[300], Amp_1_error[300], Amp_2_error[300], momentum_error[300], Amp_ratio_error[300];
   Double_t Sigma_1_error_S2[300], Amp_1_error_S2[300], Amp_2_error_S2[300], momentum_error_S2[300], Amp_ratio_error_S2[300];
   // Function range
   Double_t Range_min, Range_max;

   TH1F *Missing_Mass_S1_Projections[300];
   TH1F *Kaon_Mass_S1_Projections[300];
   TH1F *Missing_Mass_S2_Projections[300];
   TH1F *Kaon_Mass_S2_Projections[300];

   // Create arrays of the various functions
   // Strangeness 1
   TF1 *func1[300];
   TF1 *signal_1[300];
   TF1 *signal_2[300];
   TF1 *signal_total[300];
   TF1 *background_1[300];
   TF1 *background_2[300];
   TF1 *background_total[300];
   // Strangeness 2
   TF1 *func1_S2[300];
   TF1 *signal_1_S2[300];
   TF1 *signal_2_S2[300];
   TF1 *signal_total_S2[300];
   TF1 *background_1_S2[300];
   TF1 *background_2_S2[300];
   TF1 *background_total_S2[300];

   Int_t counter_S1 = 0, counter_S2 = 0;

   auto *c1 = new TCanvas("c1","canvas",800,800);
   // Parameter values to be set
   Double_t par0, par1, par2, par3, par4, par5, par6, par7, par8;
   Double_t par0_S2, par1_S2, par2_S2, par3_S2, par4_S2, par5_S2, par6_S2, par7_S2, par8_S2;

   /////////////////////////////////////////////////////////////////////////////////
   ///// Strangeness 1     /////////////////////////////////////////////////////////
   /////////////////////////////////////////////////////////////////////////////////

   for(Int_t i = 1; i < hist_S1->GetNbinsX()+1; i++){

      momentum_mid_S1 = hist_S1->GetXaxis()->GetBinCenter(i);

      // Cutting out pion cross over between 1.1 and 1.35 GeV
      if(momentum_mid_S1 > 1.1 && momentum_mid_S1 < 1.35) continue;

      // Kaons only plotted between 1-3 GeV
      if(momentum_mid_S1 < 0.9 || momentum_mid_S1 > 3) continue;


      momentum_S1[counter_S1] = momentum_mid_S1;


      Missing_Mass_S1_Projections[counter_S1] = (TH1F*) hist_S1->ProjectionY("",i,i,0,hist_S1->GetNbinsZ())->Clone();
      Kaon_Mass_S1_Projections[counter_S1] = (TH1F*) hist_S1->ProjectionZ("",i,i,0,hist_S1->GetNbinsY())->Clone();

      // Calculate integrals and ratio for S1
      S1_Integral[counter_S1] = Kaon_Mass_S1_Projections[counter_S1]->Integral();


      // Setting range of functions
      ostringstream funcname, sig1, sig2, sigtotal, back1, back2, backtotal;
      funcname <<"func_"<<i;
      sig1 <<"sig1_"<<i;
      sig2 <<"sig2_"<<i;
      back1 <<"back1_"<<i;
      back2 <<"back2_"<<i;

      // Defining functions
      func1[counter_S1] = new TF1(funcname.str().c_str(),"[0] * exp(-0.5*pow((x-[1]) / [2],2)) + [0] * [8] * exp(-0.5 * pow((x-[1]) / (2*[2]),2)) + gaus(3) + [6]*[7] + [7]*x",0.365,0.65);
      signal_1[counter_S1] = new TF1(sig1.str().c_str(),"gaus(0)",0.365,0.65);
      signal_2[counter_S1] = new TF1(sig2.str().c_str(),"[0] * (3.63259e-01 * exp(-0.5*pow((x - 1.82697e+00) / 1.82220e-01,2)) + 1.27857e+00 * exp(-0.5 * pow((x - 2.00491e+00) / 6.26903e-01,2))) * exp(-0.5 * pow((x-[1]) / (2*[2]),2))",0.365,0.65);
      signal_total[counter_S1] = new TF1(sigtotal.str().c_str(),"[0] * exp(-0.5*pow((x-[1]) / [2],2)) + [0] * [3] * exp(-0.5 * pow((x-[1]) / (2*[2]),2))",0.365,0.65);
      background_1[counter_S1] = new TF1(back1.str().c_str(),"gaus(0)",0.0,0.65);
      background_2[counter_S1] = new TF1(back2.str().c_str(),"[0]*[1] + [1]*x",0.365,0.65);
      background_total[counter_S1] = new TF1(backtotal.str().c_str(),"gaus(0) + [3]*[4] + [4]*x",0.365,0.65);


      // Defining parameters for fits
      // par0 = 0.5 * Kaon_Mass_S1_Projections[counter_S1]->GetBinContent(Kaon_Mass_S1_Projections[counter_S1]->FindBin(0.5));
      par0 = 1179.21*exp(-0.5*pow((momentum_S1[counter_S1]-0.98)/0.3379,2))+960.2+-128.8*momentum_S1[counter_S1];
      par1 = 4.95277e-01 -2.78269e-04*momentum_S1[counter_S1] - 2.19406e-03 * exp(-0.5*pow((momentum_S1[counter_S1]-1.79144e+00) / 2.70089e-01,2)); // mean for both signal gauss
      par2 = exp(-7.03189e+00+1.34979e+00*momentum_S1[counter_S1])+exp(-5.03014e+00+9.75303e-02*momentum_S1[counter_S1]);
      par3 = 2 * Kaon_Mass_S1_Projections[counter_S1]->GetBinContent(Kaon_Mass_S1_Projections[counter_S1]->FindBin(0.365)); // amplitude for back gaus
      par4 = 0.1396; // mean for back gaus (pion mass)
      par5 = 0.3; // sigma for back gaus
      par6 = -5; // constant for back pol1
      par7 = -10; // slope for back pol1
      par8 = 0.446603*exp(-0.5*pow((momentum_S1[counter_S1]-1.72255)/0.23118,2)) + 1.21202*exp(-0.5*pow((momentum_S1[counter_S1]-1.88177)/0.585851,2));


      // Setting parameters
      func1[counter_S1]->FixParameter(0,par0); // amplitude for 1st signal gaus
      func1[counter_S1]->FixParameter(1,par1); // mean for both signal gaus
      func1[counter_S1]->FixParameter(2,par2); // sigma for 1st signal gaus
      func1[counter_S1]->SetParameter(3,par3); // amplitude for back gaus
      func1[counter_S1]->SetParameter(4,par4); // mean for back gaus (pion mass)
      func1[counter_S1]->SetParameter(5,par5); // sigma for back gaus
      func1[counter_S1]->SetParameter(6,par6); // constant for back pol1
      func1[counter_S1]->SetParameter(7,par7); // slope for back pol1
      func1[counter_S1]->SetParameter(8,par8); // ratio between amplitude 1 and 2 for signal


      // Setting parameter limits
      // func1[counter_S1]->SetParLimits(0, 0.85 * par0, 1.15 * par0); // amplitude for sig 1st gaus
      // func1[counter_S1]->SetParLimits(2, 0.8 * par2, 1.2 * par2); // amplitude for sig 1st gaus
      func1[counter_S1]->SetParLimits(3,Kaon_Mass_S1_Projections[counter_S1]->GetBinContent(Kaon_Mass_S1_Projections[counter_S1]->FindBin(0.365)) / 5, 4 * Kaon_Mass_S1_Projections[counter_S1]->GetMaximum()); // amplitude for back gaus
      func1[counter_S1]->SetParLimits(4,par4 * 0.98, par4 * 1.02); // mean for back gaus (pion mass)
      func1[counter_S1]->SetParLimits(5,0.1, 0.8); // sigma for back gaus
      func1[counter_S1]->SetParLimits(6,-1000, -0.7); // constant for back pol1
      func1[counter_S1]->SetParLimits(7,-500, 0); // slope for back pol1


      Kaon_Mass_S1_Projections[counter_S1]->Fit(funcname.str().c_str(),"RBQ");


      Chi2 = func1[counter_S1]->GetChisquare() / func1[counter_S1]->GetNDF();


      signal_1[counter_S1]->SetParameter(0,func1[counter_S1]->GetParameter(0));
      signal_1[counter_S1]->SetParameter(1,func1[counter_S1]->GetParameter(1));
      signal_1[counter_S1]->SetParameter(2,func1[counter_S1]->GetParameter(2));
      signal_2[counter_S1]->SetParameter(0,func1[counter_S1]->GetParameter(0) * func1[counter_S1]->GetParameter(8));
      signal_2[counter_S1]->SetParameter(1,func1[counter_S1]->GetParameter(1));
      signal_2[counter_S1]->SetParameter(2,2 * func1[counter_S1]->GetParameter(2));
      signal_total[counter_S1]->SetParameter(0,func1[counter_S1]->GetParameter(0));
      signal_total[counter_S1]->SetParameter(1,func1[counter_S1]->GetParameter(1));
      signal_total[counter_S1]->SetParameter(2,func1[counter_S1]->GetParameter(2));
      signal_total[counter_S1]->SetParameter(3,func1[counter_S1]->GetParameter(8));
      background_1[counter_S1]->SetParameter(0,func1[counter_S1]->GetParameter(3));
      background_1[counter_S1]->SetParameter(1,func1[counter_S1]->GetParameter(4));
      background_1[counter_S1]->SetParameter(2,func1[counter_S1]->GetParameter(5));
      background_2[counter_S1]->SetParameter(0,func1[counter_S1]->GetParameter(6));
      background_2[counter_S1]->SetParameter(1,func1[counter_S1]->GetParameter(7));
      background_total[counter_S1]->SetParameter(0,func1[counter_S1]->GetParameter(3));
      background_total[counter_S1]->SetParameter(1,func1[counter_S1]->GetParameter(4));
      background_total[counter_S1]->SetParameter(2,func1[counter_S1]->GetParameter(5));
      background_total[counter_S1]->SetParameter(3,func1[counter_S1]->GetParameter(6));
      background_total[counter_S1]->SetParameter(4,func1[counter_S1]->GetParameter(7));


      // Get the projection fit parameters
      Amp_1[counter_S1] = func1[counter_S1]->GetParameter(0);
      Sigma_1[counter_S1] = func1[counter_S1]->GetParameter(2);
      // Amp_2[counter_S1] = func1[counter_S1]->GetParameter(0) * func1[counter_S1]->GetParameter(8);
      Amp_2[counter_S1] = func1[counter_S1]->GetParameter(8);

      // Get the projection fit errors
      Amp_1_error[counter_S1] = func1[counter_S1]->GetParError(0);
      Sigma_1_error[counter_S1] = func1[counter_S1]->GetParError(2);
      Amp_2_error[counter_S1] = func1[counter_S1]->GetParError(0);



      h_chi2->Fill(momentum_S1[counter_S1], Chi2);
      h_sigma_1->Fill(momentum_S1[counter_S1],Sigma_1[counter_S1]);
      h_sigma_1->SetBinError(i, Sigma_1_error[counter_S1]);
      h_amp_1->Fill(momentum_S1[counter_S1],Amp_1[counter_S1]);
      h_amp_1->SetBinError(i, Amp_1_error[counter_S1]);
      h_amp_2->Fill(momentum_S1[counter_S1],Amp_2[counter_S1]);
      h_amp_2->SetBinError(i, Amp_2_error[counter_S1]);
      h_mean_signal->Fill(momentum_S1[counter_S1],func1[counter_S1]->GetParameter(1));
      h_mean_signal->SetBinError(momentum_S1[counter_S1],func1[counter_S1]->GetParError(1));


      h_S1_Integral->Fill(momentum_S1[counter_S1], S1_Integral[counter_S1]);


      counter_S1++;

   }



   for(Int_t i = 1; i < hist_S2->GetNbinsX()+1; i++){

      momentum_mid_S2 = hist_S2->GetXaxis()->GetBinCenter(i);

      // Cutting out pion cross over between 1.1 and 1.35 GeV
      if(momentum_mid_S2 > 1.1 && momentum_mid_S2 < 1.35) continue;

      // Kaons only plotted between 0.9-3 GeV
      if(momentum_mid_S2 < 0.9 || momentum_mid_S2 > 3) continue;


      momentum_S2[counter_S2] = momentum_mid_S2;


      Missing_Mass_S2_Projections[counter_S2] = (TH1F*) hist_S2->ProjectionY("",i,i,0,hist_S2->GetNbinsZ())->Clone();
      Kaon_Mass_S2_Projections[counter_S2] = (TH1F*) hist_S2->ProjectionZ("",i,i,0,hist_S2->GetNbinsY())->Clone();

      // Calculate integrals and ratio for S2
      S2_Integral[counter_S2] = Kaon_Mass_S2_Projections[counter_S2]->Integral();


      // Setting range of functions
      ostringstream funcname_S2, sig1_S2, sig2_S2, sigtotal_S2, back1_S2, back2_S2, backtotal_S2;
      funcname_S2 <<"func_S2_"<<i;
      sig1_S2 <<"sig1_S2_"<<i;
      sig2_S2 <<"sig2_S2_"<<i;
      back1_S2 <<"back1_S2_"<<i;
      back2_S2 <<"back2_S2_"<<i;

      // Defining functions
      func1_S2[counter_S2] = new TF1(funcname_S2.str().c_str(),"[0] * exp(-0.5*pow((x-[1]) / [2],2)) + [0] * [8] * exp(-0.5 * pow((x-[1]) / (2*[2]),2)) + gaus(3) + [6]*[7] + [7]*x",0.365,0.65);
      signal_1_S2[counter_S2] = new TF1(sig1_S2.str().c_str(),"gaus(0)",0.365,0.65);
      signal_2_S2[counter_S2] = new TF1(sig2_S2.str().c_str(),"[0] * (3.63259e-01 * exp(-0.5*pow((x - 1.82697e+00) / 1.82220e-01,2)) + 1.27857e+00 * exp(-0.5 * pow((x - 2.00491e+00) / 6.26903e-01,2))) * exp(-0.5 * pow((x-[1]) / (2*[2]),2))",0.365,0.65);
      signal_total_S2[counter_S2] = new TF1(sigtotal_S2.str().c_str(),"[0] * exp(-0.5*pow((x-[1]) / [2],2)) + [0] * (3.63259e-01 * exp(-0.5*pow((x - 1.82697e+00) / 1.82220e-01,2)) + 1.27857e+00 * exp(-0.5 * pow((x - 2.00491e+00) / 6.26903e-01,2))) * exp(-0.5 * pow((x-[1]) / (2*[2]),2))",0.365,0.65);
      background_1_S2[counter_S2] = new TF1(back1_S2.str().c_str(),"gaus(0)",0.0,0.65);
      background_2_S2[counter_S2] = new TF1(back2_S2.str().c_str(),"[0]*[1] + [1]*x",0.365,0.65);
      background_total_S2[counter_S2] = new TF1(backtotal_S2.str().c_str(),"gaus(0) + [3]*[4] + [4]*x",0.365,0.65);


      // Defining parameters for fits
      par0_S2 = 10.6785*exp(-0.5*pow((momentum_S2[counter_S2]-1.30712)/2.43071e-01,2))+18.6151+-3.27781*momentum_S2[counter_S2];
      par1_S2 = 4.95277e-01 -2.78269e-04*momentum_S2[counter_S2] - 2.19406e-03 * exp(-0.5*pow((momentum_S2[counter_S2]-1.79144e+00) / 2.70089e-01,2)); // mean for both signal gauss
      par2_S2 = exp(-5.87651e+00+9.50189e-01*momentum_S2[counter_S2])+exp(-4.83642e+00+-9.34420e-01*momentum_S2[counter_S2]);
      par3_S2 = 2 * Kaon_Mass_S2_Projections[counter_S2]->GetBinContent(Kaon_Mass_S2_Projections[counter_S2]->FindBin(0.365)); // amplitude for back gaus
      par4_S2 = 0.1396; // mean for back gaus (pion mass)
      par5_S2 = 0.3; // sigma for back gaus
      par6_S2 = -5; // constant for back pol1
      par7_S2 = -10; // slope for back pol1
      par8_S2 = 3.63259e-01*exp(-0.5*pow((momentum_S2[counter_S2]-1.82697e+00)/1.82220e-01,2)) + 1.27857e+00*exp(-0.5*pow((momentum_S2[counter_S2]-2.00491e+00)/6.26903e-01,2));


      // Setting parameters
      func1_S2[counter_S2]->SetParameter(0,par0_S2); // amplitude for 1st signal gaus
      func1_S2[counter_S2]->FixParameter(1,par1_S2); // mean for both signal gaus
      func1_S2[counter_S2]->FixParameter(2,par2_S2); // sigma for 1st signal gaus
      func1_S2[counter_S2]->SetParameter(3,par3_S2); // amplitude for back gaus
      func1_S2[counter_S2]->SetParameter(4,par4_S2); // mean for back gaus (pion mass_S2)
      func1_S2[counter_S2]->SetParameter(5,par5_S2); // sigma for back gaus
      func1_S2[counter_S2]->SetParameter(6,par6_S2); // constant for back pol1
      func1_S2[counter_S2]->SetParameter(7,par7_S2); // slope for back pol1
      func1_S2[counter_S2]->FixParameter(8,par8_S2); // ratio between amplitude 1 and 2 for signal


      // Setting parameter limits
      func1_S2[counter_S2]->SetParLimits(3,Kaon_Mass_S2_Projections[counter_S2]->GetBinContent(Kaon_Mass_S2_Projections[counter_S2]->FindBin(0.365)) / 5, 4 * Kaon_Mass_S2_Projections[counter_S2]->GetMaximum()); // amplitude for back gaus
      func1_S2[counter_S2]->SetParLimits(4,par4_S2 * 0.98, par4_S2 * 1.02); // mean for back gaus (pion mass)
      func1_S2[counter_S2]->SetParLimits(5,0.1, 0.8); // sigma for back gaus
      func1_S2[counter_S2]->SetParLimits(6,-1000, -0.7); // constant for back pol1
      func1_S2[counter_S2]->SetParLimits(7,-500, 0); // slope for back pol1


      Kaon_Mass_S2_Projections[counter_S2]->Fit(funcname_S2.str().c_str(),"RBQ");


      signal_1_S2[counter_S2]->SetParameter(0,func1_S2[counter_S2]->GetParameter(0));
      signal_1_S2[counter_S2]->SetParameter(1,func1_S2[counter_S2]->GetParameter(1));
      signal_1_S2[counter_S2]->SetParameter(2,func1_S2[counter_S2]->GetParameter(2));
      signal_2_S2[counter_S2]->SetParameter(0,func1_S2[counter_S2]->GetParameter(0));
      signal_2_S2[counter_S2]->SetParameter(1,func1_S2[counter_S2]->GetParameter(1));
      signal_2_S2[counter_S2]->SetParameter(2,func1_S2[counter_S2]->GetParameter(2));
      signal_total_S2[counter_S2]->SetParameter(0,func1_S2[counter_S2]->GetParameter(0));
      signal_total_S2[counter_S2]->SetParameter(1,func1_S2[counter_S2]->GetParameter(1));
      signal_total_S2[counter_S2]->SetParameter(2,func1_S2[counter_S2]->GetParameter(2));
      background_1_S2[counter_S2]->SetParameter(0,func1_S2[counter_S2]->GetParameter(3));
      background_1_S2[counter_S2]->SetParameter(1,func1_S2[counter_S2]->GetParameter(4));
      background_1_S2[counter_S2]->SetParameter(2,func1_S2[counter_S2]->GetParameter(5));
      background_2_S2[counter_S2]->SetParameter(0,func1_S2[counter_S2]->GetParameter(6));
      background_2_S2[counter_S2]->SetParameter(1,func1_S2[counter_S2]->GetParameter(7));
      background_total_S2[counter_S2]->SetParameter(0,func1_S2[counter_S2]->GetParameter(3));
      background_total_S2[counter_S2]->SetParameter(1,func1_S2[counter_S2]->GetParameter(4));
      background_total_S2[counter_S2]->SetParameter(2,func1_S2[counter_S2]->GetParameter(5));
      background_total_S2[counter_S2]->SetParameter(3,func1_S2[counter_S2]->GetParameter(6));
      background_total_S2[counter_S2]->SetParameter(4,func1_S2[counter_S2]->GetParameter(7));


      // Get the projection fit parameters
      Amp_1_S2[counter_S2] = func1_S2[counter_S2]->GetParameter(0);
      Sigma_1_S2[counter_S2] = func1_S2[counter_S2]->GetParameter(2);
      // Amp_2_S2[counter_S2] = func1_S2[counter_S2]->GetParameter(0) * func1_S2[counter_S2]->GetParameter(8);
      Amp_2_S2[counter_S2] = func1_S2[counter_S2]->GetParameter(8);

      // Get the projection fit errors
      Amp_1_error_S2[counter_S2] = func1_S2[counter_S2]->GetParError(0);
      Sigma_1_error_S2[counter_S2] = func1_S2[counter_S2]->GetParError(2);
      Amp_2_error_S2[counter_S2] = func1_S2[counter_S2]->GetParError(0);



      h_amp_1_S2->Fill(momentum_S2[counter_S2],Amp_1_S2[counter_S2]);
      h_amp_1_S2->SetBinError(i, Amp_1_error_S2[counter_S2]);
      h_amp_2_S2->Fill(momentum_S2[counter_S2],Amp_2_S2[counter_S2]);
      h_amp_2_S2->SetBinError(i, Amp_2_error_S2[counter_S2]);


      h_S2_Integral->Fill(momentum_S2[counter_S2], S2_Integral[counter_S2]);

      counter_S2++;

   }

   auto *h_S1_S2_Ratio = new TH1F();
   h_S1_S2_Ratio = (TH1F*)h_S1_Integral->Clone("h_S1_S2_Ratio");
   h_S1_S2_Ratio->Divide(h_S2_Integral);

   // Creating TGraphs for the parameters
   // Strangeness 1
   TGraphErrors* gr1 = new TGraphErrors(counter_S1, momentum_S1, Sigma_1, 0, Sigma_1_error);
   TGraphErrors* gr2 = new TGraphErrors(counter_S1, momentum_S1, Amp_1, 0, Amp_1_error);
   TGraphErrors* gr3 = new TGraphErrors(counter_S1, momentum_S1, Amp_2, 0, Amp_2_error);
   TGraphErrors* gr4 = new TGraphErrors(counter_S1, momentum_S1, amp_ratio, 0, Amp_ratio_error);
   // Strangeness 2
   TGraphErrors* gr5 = new TGraphErrors(counter_S2, momentum_S2, Sigma_1_S2, 0, Sigma_1_error_S2);
   TGraphErrors* gr6 = new TGraphErrors(counter_S2, momentum_S2, Amp_1_S2, 0, Amp_1_error_S2);
   TGraphErrors* gr7 = new TGraphErrors(counter_S2, momentum_S2, Amp_2_S2, 0, Amp_2_error_S2);

   // Creating fits for parameters
   // Strageness 1
   TF1 *func_sig_1_amp = new TF1("func_sig_1_amp","gaus(0) + pol1(3)",0.92,2.6);
   TF1 *func_sig_1_mean = new TF1("func_sig_1_mean","pol1(0) - gaus(2)",0.92,2.6);
   TF1 *func_sig_1_sigma = new TF1("func_sig_1_sigma","expo(0) + expo(2)",0.92,2.6);
   TF1 *func_sig_2_amp = new TF1("func_sig_2_amp","gaus(0) + gaus(3)",0.92,2.6);
   // Strageness 2
   TF1 *func_sig_1_amp_S2 = new TF1("func_sig_1_amp_S2","gaus(0) + pol1(3)",0.92,2.6);
   TF1 *func_sig_1_mean_S2 = new TF1("func_sig_1_mean_S2","pol1(0) - gaus(2)",0.92,2.6);
   TF1 *func_sig_1_sigma_S2 = new TF1("func_sig_1_sigma_S2","expo(0) + expo(2)",0.92,2.6);
   TF1 *func_sig_2_amp_S2 = new TF1("func_sig_2_amp_S2","gaus(0) + gaus(3)",0.92,2.6);


   // Setting parameters for fitting parameters

   // signal 1 amplitude
   // gaus + pol1
   func_sig_1_amp->SetParameter(0,1179.21);
   func_sig_1_amp->SetParLimits(0,0.8 * func_sig_1_amp->GetParameter(0), 1.2 * func_sig_1_amp->GetParameter(0));
   func_sig_1_amp->SetParameter(1,0.98);
   func_sig_1_amp->SetParLimits(1,0.98,1.35);
   func_sig_1_amp->SetParameter(2,0.3379);
   func_sig_1_amp->SetParameter(3,960.2);
   func_sig_1_amp->SetParameter(4,-128.8);


   // signal 1 mean
   // pol1 - gaus
   func_sig_1_mean->FixParameter(0,4.95277e-01);
   func_sig_1_mean->FixParameter(1,-2.78269e-04);
   func_sig_1_mean->FixParameter(2,2.19406e-03);
   func_sig_1_mean->FixParameter(3,1.79144e+00);
   func_sig_1_mean->FixParameter(4,2.70089e-01);

   // signal 1 sigma
   // 2 exponentials
   func_sig_1_sigma->SetParameter(0,-7.03189e+00);
   func_sig_1_sigma->SetParameter(1,1.34979e+00);
   func_sig_1_sigma->SetParameter(2,-5.03014e+00);
   func_sig_1_sigma->SetParameter(3,9.75303e-02);

   // signal 2 amplitude ratio to signal 1 amplitude
   func_sig_2_amp->FixParameter(0,3.63259e-01);
   func_sig_2_amp->FixParameter(1,1.82697e+00);
   func_sig_2_amp->FixParameter(2,1.82220e-01);
   func_sig_2_amp->FixParameter(3,1.27857);
   func_sig_2_amp->FixParameter(4,2.00491e+00);
   func_sig_2_amp->FixParameter(5,6.26903e-01);

   // Strangeness 2
   // signal 1 amplitude
   // gaus + pol1
   func_sig_1_amp_S2->SetParameter(0,10);
   func_sig_1_amp_S2->SetParLimits(0,5,20);
   func_sig_1_amp_S2->SetParameter(1,1.14598e+00);
   func_sig_1_amp_S2->SetParLimits(1,1.1,1.35);
   func_sig_1_amp_S2->SetParameter(2,3.10873e-01);
   func_sig_1_amp_S2->SetParameter(3,3.43886e+02);
   func_sig_1_amp_S2->SetParameter(4,-4.66091e+01);


   // signal 1 mean
   // pol1 - gaus
   func_sig_1_mean_S2->FixParameter(0,4.95277e-01);
   func_sig_1_mean_S2->FixParameter(1,-2.78269e-04);
   func_sig_1_mean_S2->FixParameter(2,2.19406e-03);
   func_sig_1_mean_S2->FixParameter(3,1.79144e+00);
   func_sig_1_mean_S2->FixParameter(4,2.70089e-01);

   // signal 1 sigma
   // 2 exponentials
   func_sig_1_sigma_S2->FixParameter(0,-5.87651e+00);
   func_sig_1_sigma_S2->FixParameter(1,9.50189e-01);
   func_sig_1_sigma_S2->FixParameter(2,-4.83642e+00);
   func_sig_1_sigma_S2->FixParameter(3,-9.34420e-01);

   // signal 2 amplitude
   // signal 1 amp * 2 gaus
   func_sig_2_amp_S2->SetParameter(0,3.63259e-01);
   func_sig_2_amp_S2->SetParameter(1,1.82697e+00);
   func_sig_2_amp_S2->SetParameter(2,1.82220e-01);
   func_sig_2_amp_S2->SetParameter(3,1.27857e+00);
   func_sig_2_amp_S2->SetParameter(4,2.00491e+00);
   func_sig_2_amp_S2->SetParameter(5,6.26903e-01);

   // Fitting parameter functions
   gr1->Fit("func_sig_1_sigma","RBQ");
   gr2->Fit("func_sig_1_amp","RBQ");
   gr3->Fit("func_sig_2_amp","RBQ");
   gr5->Fit("func_sig_1_sigma_S2","RBQ");
   gr6->Fit("func_sig_1_amp_S2","RBQ");
   gr7->Fit("func_sig_2_amp_S2","RBQ");


   // Creating strings for parameter values
   ostringstream signal_1_amp, signal_mean, signal_1_sigma, signal_2_amp, signal_2_sigma;
   ostringstream signal_1_amp_S2, signal_mean_S2, signal_1_sigma_S2, signal_2_amp_S2, signal_2_sigma_S2;

   TH2F *Signal_Function = new TH2F("Signal_Function","",500,0.3,0.8,300,0,3);
   TH2F *Background_Function = new TH2F("Background_Function","",500,0.3,0.8,300,0,3);
   TF1 *signal_function = new TF1("signal_function","gaus(0) + gaus(3)",0.3,0.8);
   TH2F *Signal_Function_S2 = new TH2F("Signal_Function_S2","",250,0.3,0.8,150,0,3);
   TH2F *Background_Function_S2 = new TH2F("Background_Function_S2","",250,0.3,0.8,150,0,3);
   TF1 *signal_function_S2 = new TF1("signal_function_S2","gaus(0) + gaus(3)",0.3,0.8);
   TH2F *Result = new TH2F();
   TH2F *Result_S2 = new TH2F();

   // Take copy of the data kaon momentum vs mass
   TH2F *Kaon_Mass_Momentum = new TH2F();
   Kaon_Mass_Momentum = (TH2F*)hist_S1->Project3D("xz")->Clone();

   TH2F *Kaon_Mass_Momentum_S2 = new TH2F();
   Kaon_Mass_Momentum_S2 = (TH2F*)hist_S2->Project3D("xz")->Clone();


   TH3F *S1_Signal = new TH3F();
   S1_Signal = (TH3F*)hist_S1->Clone();
   TH3F *S1_Background = new TH3F();
   S1_Background = (TH3F*)hist_S1->Clone();
   TH3F *S2_Signal = new TH3F();
   S2_Signal = (TH3F*)hist_S2->Clone();
   TH3F *S2_Background = new TH3F();
   S2_Background = (TH3F*)hist_S2->Clone();

   // S1_Signal->GetZaxis()->SetRange(S1_Signal->GetZaxis()->FindBin(0.36),S1_Signal->GetZaxis()->FindBin(0.7));
   // S1_Background->GetZaxis()->SetRange(S1_Background->GetZaxis()->FindBin(0.36),S1_Background->GetZaxis()->FindBin(0.7));

   // Creating signal function and histogram

   // Strangeness 1
   // Loop over kaon momentum
   for(Int_t x_pos=1; x_pos < hist_S1->GetNbinsX() + 1; x_pos++){

      // Set signal function parameters based on values obtained from the 1D fit
      // parameter functions for current bin in momentum
      signal_function->SetParameter(0,func_sig_1_amp->Eval(Signal_Function->GetYaxis()->GetBinCenter(x_pos)));
      signal_function->SetParameter(1,func_sig_1_mean->Eval(Signal_Function->GetYaxis()->GetBinCenter(x_pos)));
      signal_function->SetParameter(2,func_sig_1_sigma->Eval(Signal_Function->GetYaxis()->GetBinCenter(x_pos)));
      signal_function->SetParameter(3,func_sig_1_amp->Eval(Signal_Function->GetYaxis()->GetBinCenter(x_pos)) * func_sig_2_amp->Eval(Signal_Function->GetYaxis()->GetBinCenter(x_pos)));
      signal_function->SetParameter(4,func_sig_1_mean->Eval(Signal_Function->GetYaxis()->GetBinCenter(x_pos)));
      signal_function->SetParameter(5,2 * func_sig_1_sigma->Eval(Signal_Function->GetYaxis()->GetBinCenter(x_pos)));


      // Loop over kaon mass
      for(Int_t z_pos=1; z_pos < hist_S1->GetNbinsZ() + 1; z_pos++){

         // Set bin content for 2D signal histogram based on current momentum
         // and kaon mass bin
         Signal_Function->SetBinContent(z_pos,x_pos,signal_function->Eval(Signal_Function->GetXaxis()->GetBinCenter(z_pos),Signal_Function->GetYaxis()->GetBinCenter(x_pos)));

         // Make sure there are no negative values for the background
         if(Kaon_Mass_Momentum->GetBinContent(z_pos,x_pos) - Signal_Function->GetBinContent(z_pos,x_pos)>0){
            Background_Function->SetBinContent(z_pos,x_pos,Kaon_Mass_Momentum->GetBinContent(z_pos,x_pos) - Signal_Function->GetBinContent(z_pos,x_pos));
         }
         else
         {
            Background_Function->SetBinContent(z_pos,x_pos,0);
         }

         // Loop over missing mass
         for(Int_t y_pos=1; y_pos < hist_S1->GetNbinsY() + 1; y_pos++){

            if(S1_Signal->GetBinContent(x_pos,y_pos,z_pos) > 1){

               S1_Signal->SetBinContent(x_pos, y_pos, z_pos, Signal_Function->GetBinContent(x_pos,z_pos) * S1_Signal->GetBinContent(x_pos,y_pos,z_pos));

               // Make sure there are no negative values for the background
               // if(Background_Function->GetBinContent(x_pos,z_pos) - Signal_Function->GetBinContent(x_pos,z_pos) > 0){
               //    S1_Background->SetBinContent(x_pos, y_pos, z_pos, (Background_Function->GetBinContent(x_pos,z_pos) - Signal_Function->GetBinContent(x_pos,z_pos)) * S1_Background->GetBinContent(x_pos,y_pos,z_pos));
               // }
               // else{
               //    S1_Background->SetBinContent(x_pos, y_pos, z_pos, 0);
               // }
               S1_Background->SetBinContent(x_pos, y_pos, z_pos, Background_Function->GetBinContent(x_pos,z_pos) * S1_Background->GetBinContent(x_pos,y_pos,z_pos));
            }
         }
      }
   }


   // Strangeness 2
   // Loop over kaon momentum
   for(Int_t x_pos=1; x_pos < hist_S2->GetNbinsX() + 1; x_pos++){

      // Set signal function parameters based on values obtained from the 1D fit
      // parameter functions for current bin in momentum
      signal_function_S2->SetParameter(0,func_sig_1_amp_S2->Eval(Signal_Function_S2->GetYaxis()->GetBinCenter(x_pos)));
      signal_function_S2->SetParameter(1,func_sig_1_mean_S2->Eval(Signal_Function_S2->GetYaxis()->GetBinCenter(x_pos)));
      signal_function_S2->SetParameter(2,func_sig_1_sigma_S2->Eval(Signal_Function_S2->GetYaxis()->GetBinCenter(x_pos)));
      signal_function_S2->SetParameter(3,func_sig_1_amp_S2->Eval(Signal_Function_S2->GetYaxis()->GetBinCenter(x_pos)) * func_sig_2_amp_S2->Eval(Signal_Function_S2->GetYaxis()->GetBinCenter(x_pos)));
      signal_function_S2->SetParameter(4,func_sig_1_mean_S2->Eval(Signal_Function_S2->GetYaxis()->GetBinCenter(x_pos)));
      signal_function_S2->SetParameter(5,2 * func_sig_1_sigma_S2->Eval(Signal_Function_S2->GetYaxis()->GetBinCenter(x_pos)));


      // Loop over kaon mass
      for(Int_t z_pos=1; z_pos < hist_S2->GetNbinsZ() + 1; z_pos++){

         // Set bin content for 2D signal histogram based on current momentum
         // and kaon mass bin
         Signal_Function_S2->SetBinContent(z_pos,x_pos,signal_function_S2->Eval(Signal_Function_S2->GetXaxis()->GetBinCenter(z_pos),Signal_Function_S2->GetYaxis()->GetBinCenter(x_pos)));

         if(Kaon_Mass_Momentum_S2->GetBinContent(z_pos,x_pos) - Signal_Function_S2->GetBinContent(z_pos,x_pos)>0){
            Background_Function_S2->SetBinContent(z_pos,x_pos,Kaon_Mass_Momentum_S2->GetBinContent(z_pos,x_pos) - Signal_Function_S2->GetBinContent(z_pos,x_pos));
         }
         else
         {
            Background_Function_S2->SetBinContent(z_pos,x_pos,0);
         }

         // Loop over missing mass
         for(Int_t y_pos=1; y_pos < hist_S2->GetNbinsY() + 1; y_pos++){

            if(S2_Signal->GetBinContent(x_pos,y_pos,z_pos) > 1){

               S2_Signal->SetBinContent(x_pos, y_pos, z_pos, Signal_Function_S2->GetBinContent(x_pos,z_pos) * S2_Signal->GetBinContent(x_pos,y_pos,z_pos));

               // Make sure there are no negative values for the background
               // if(Background_Function_S2->GetBinContent(x_pos,z_pos) - Signal_Function_S2->GetBinContent(x_pos,z_pos) > 0){
               //    S2_Background->SetBinContent(x_pos, y_pos, z_pos, (Background_Function_S2->GetBinContent(x_pos,z_pos) - Signal_Function_S2->GetBinContent(x_pos,z_pos)) * S2_Background->GetBinContent(x_pos,y_pos,z_pos));
               // }
               // else{
               //    S2_Background->SetBinContent(x_pos, y_pos, z_pos, 0);
               // }
               S2_Background->SetBinContent(x_pos, y_pos, z_pos, Background_Function_S2->GetBinContent(x_pos,z_pos) * S2_Background->GetBinContent(x_pos,y_pos,z_pos));

            }
         }
      }
   }



   // Double_t Scale_Low_MM = Signal_Function->Integral(1,Signal_Function->GetNbinsX(),Signal_Function->GetYaxis()->FindBin(0.0),Signal_Function->GetYaxis()->FindBin(0.9)) / Background_Function->Integral(1,Signal_Function->GetNbinsX(),Background_Function->GetYaxis()->FindBin(0.0),Background_Function->GetYaxis()->FindBin(0.9));
   // Background_Function->Scale(Scale_Low_MM);
   // Double_t Scale_Low_MM_S2 = Signal_Function_S2->Integral(1,Signal_Function_S2->GetNbinsX(),Signal_Function_S2->GetYaxis()->FindBin(0.0),Signal_Function_S2->GetYaxis()->FindBin(0.9)) / Background_Function_S2->Integral(1,Signal_Function_S2->GetNbinsX(),Background_Function_S2->GetYaxis()->FindBin(0.0),Background_Function_S2->GetYaxis()->FindBin(0.9));
   // Background_Function_S2->Scale(Scale_Low_MM_S2);

   Result = (TH2F*)Kaon_Mass_Momentum->Clone("Result");

   Result_S2 = (TH2F*)Kaon_Mass_Momentum_S2->Clone("Result_S2");

   Result->Add(Background_Function,-1);
   Result_S2->Add(Background_Function_S2,-1);

   Signal_Function->Scale(Kaon_Mass_Momentum->Integral() / (Background_Function->Integral() + Result->Integral()));
   Background_Function->Scale(Kaon_Mass_Momentum->Integral() / (Background_Function->Integral() + Result->Integral()));

   Result->Divide(Kaon_Mass_Momentum);
   Result_S2->Divide(Kaon_Mass_Momentum_S2);


   fileOutput1.Write();

}
