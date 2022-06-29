// Making lots of projections in momentum

{

   Int_t S2_x_rebin=2;
   Int_t S2_y_rebin=1;
   Int_t S2_z_rebin=2;
   Int_t S3_x_rebin=10;
   Int_t S3_y_rebin=20;
   Int_t S3_z_rebin=30;


   // Get input file
   // TFile *f1 = new TFile("/media/mn688/Elements1/PhD/Analysis_Output/Hexaquark/RGA_Fall2018_Outbending_at_least_1e1KpFD_Tree_Total_24022022_Total_Scaling_16032022_01.root");
   TFile *f1=new TFile("/media/mn688/Elements1/PhD/Analysis_Output/Hexaquark/RGB_Spring2019_Inbending_at_least_1e1KpFD_Tree_Total_03032022_Total_Scaling_23032022_01.root");
   // Get histogram from input file
   TH3F* hist_S1=(TH3F*)f1->Get("h_S1_Kaon_Momentum__Miss_Mass__Kaon_Mass");
   TH3F* hist_S2=(TH3F*)f1->Get("h_S2_Kaon_Momentum__Miss_Mass__Kaon_Mass");
   TH3F* hist_S3=(TH3F*)f1->Get("h_S3_Kaon_Momentum__Miss_Mass__Kaon_Mass");

   // Rebin strangeness 2 due to low statistics
   hist_S2->Rebin3D(S2_x_rebin,S2_y_rebin,S2_z_rebin);

   // Rebin strangeness 3 due to low statistics
   hist_S3->Rebin3D(S3_x_rebin,S3_y_rebin,S3_z_rebin);

   TFile fileOutput1("Kaon_Background_Subtraction_RGB_Spring2019_Inbending_28062022_01.root","recreate");

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
   Double_t bin_1[300], bin_2[300], bin_3[300], momentum_S1[300], momentum_S2[300], momentum_S3[300], momentum_mid_S1, momentum_mid_S2, momentum_mid_S3;
   // Ratio between S1 and S2
   Double_t S1_Integral[300], S2_Integral[300];
   // Reduced chi^2
   Double_t Chi2;
   // Parameter values
   Double_t Sigma_1[300], Amp_1[300], Amp_2[300], amp_ratio[300];
   Double_t Back_Amp_S1[300], Back_Mean_S1[300], Back_Sigma_S1[300], Back_Constant_S1[300], Back_Gradient_S1[300];
   Double_t Sigma_1_S2[300], Amp_1_S2[300], Amp_2_S2[300], amp_ratio_S2[300];
   Double_t Sigma_1_S3[300], Amp_1_S3[300], Amp_2_S3[300], amp_ratio_S3[300];
   // Parameter errors
   Double_t Sigma_1_error[300], Amp_1_error[300], Amp_2_error[300], momentum_error[300], Amp_ratio_error[300];
   Double_t Back_Amp_S1_error[300], Back_Mean_S1_error[300], Back_Sigma_S1_error[300], Back_Constant_S1_error[300], Back_Gradient_S1_error[300];
   Double_t Sigma_1_error_S2[300], Amp_1_error_S2[300], Amp_2_error_S2[300], momentum_error_S2[300], Amp_ratio_error_S2[300];
   Double_t Sigma_1_error_S3[300], Amp_1_error_S3[300], Amp_2_error_S3[300], momentum_error_S3[300], Amp_ratio_error_S3[300];
   // Function range
   Double_t Range_min, Range_max;

   TH1F *Missing_Mass_S1_Projections[300];
   TH1F *Kaon_Mass_S1_Projections[300];
   TH1F *Missing_Mass_S2_Projections[300];
   TH1F *Kaon_Mass_S2_Projections[300];
   TH1F *Missing_Mass_S3_Projections[300];
   TH1F *Kaon_Mass_S3_Projections[300];

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
   // Strangeness 3
   TF1 *func1_S3[300];
   TF1 *signal_1_S3[300];
   TF1 *signal_2_S3[300];
   TF1 *signal_total_S3[300];
   TF1 *background_1_S3[300];
   TF1 *background_2_S3[300];
   TF1 *background_total_S3[300];

   Int_t counter_S1 = 0, counter_S2 = 0, counter_S3 = 0;

   auto *c1 = new TCanvas("c1","canvas",800,800);
   // Parameter values to be set
   Double_t par0, par1, par2, par3, par4, par5, par6, par7, par8;
   Double_t par0_S2, par1_S2, par2_S2, par3_S2, par4_S2, par5_S2, par6_S2, par7_S2, par8_S2;
   Double_t par0_S3, par1_S3, par2_S3, par3_S3, par4_S3, par5_S3, par6_S3, par7_S3, par8_S3;

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
      par0 = 1090*exp(-0.5*pow((momentum_S1[counter_S1]-1.052)/0.272,2))+1299+-322.3*momentum_S1[counter_S1];
      par1 = 4.95277e-01 -2.78269e-04*momentum_S1[counter_S1] - 2.19406e-03 * exp(-0.5*pow((momentum_S1[counter_S1]-1.79144e+00) / 2.70089e-01,2)); // mean for both signal gauss
      par2 = exp(-5.87651e+00+9.50189e-01*momentum_S1[counter_S1])+exp(-4.83642e+00+-9.34420e-01*momentum_S1[counter_S1]);
      par3 = -2466.73 + 8254.67*momentum_S1[counter_S1] -8574.25*pow(momentum_S1[counter_S1],2)+3580.97*pow(momentum_S1[counter_S1],3) -511.388*pow(momentum_S1[counter_S1],4); // amplitude for back gaus
      par4 = 0.1396; // mean for back gaus (pion mass)
      // par5 = 0.3; // sigma for back gaus
      par5 = 0.393*exp(-0.5*pow((momentum_S1[counter_S1]-1.46)/0.613,2))+-0.534+0.224*momentum_S1[counter_S1]; // sigma for back gaus
      par6 = -5; // constant for back pol1
      par7 = -10; // slope for back pol1
      par8 = 0.578473*exp(-0.5*pow((momentum_S1[counter_S1]-1.66762)/0.229272,2)) + 1.34606*exp(-0.5*pow((momentum_S1[counter_S1]-2.12122)/0.701563,2));


      // Setting parameters
      func1[counter_S1]->FixParameter(0,par0); // amplitude for 1st signal gaus
      func1[counter_S1]->FixParameter(1,par1); // mean for both signal gaus
      func1[counter_S1]->FixParameter(2,par2); // sigma for 1st signal gaus
      func1[counter_S1]->SetParameter(3,par3); // amplitude for back gaus
      func1[counter_S1]->FixParameter(4,par4); // mean for back gaus (pion mass)
      func1[counter_S1]->SetParameter(5,par5); // sigma for back gaus
      func1[counter_S1]->SetParameter(6,par6); // constant for back pol1
      func1[counter_S1]->SetParameter(7,par7); // slope for back pol1
      func1[counter_S1]->FixParameter(8,par8); // ratio between amplitude 1 and 2 for signal


      // Setting parameter limits
      // func1[counter_S1]->SetParLimits(3,0, Kaon_Mass_S1_Projections[counter_S1]->GetBinContent(Kaon_Mass_S1_Projections[counter_S1]->FindBin(0.365)) * 10); // amplitude for back gaus
      // func1[counter_S1]->SetParLimits(3,par3 * 0.75, par3 * 1.25); // amplitude for back gaus
      // func1[counter_S1]->SetParLimits(4,par4 * 0.6, par4 * 1.4); // mean for back gaus (pion mass)
      // func1[counter_S1]->SetParLimits(4,par4 * 0.85, par4 * 1.15); // mean for back gaus (pion mass)
      // func1[counter_S1]->SetParLimits(5,par5*0.9, par5*1.1); // sigma for back gaus
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
      Back_Amp_S1[counter_S1] = func1[counter_S1]->GetParameter(3);
      Back_Mean_S1[counter_S1] = func1[counter_S1]->GetParameter(4);
      Back_Sigma_S1[counter_S1] = func1[counter_S1]->GetParameter(5);
      Back_Constant_S1[counter_S1] = func1[counter_S1]->GetParameter(6);
      Back_Gradient_S1[counter_S1] = func1[counter_S1]->GetParameter(7);

      // Get the projection fit errors
      Amp_1_error[counter_S1] = func1[counter_S1]->GetParError(0);
      Sigma_1_error[counter_S1] = func1[counter_S1]->GetParError(2);
      Amp_2_error[counter_S1] = func1[counter_S1]->GetParError(0);
      Back_Amp_S1_error[counter_S1] = func1[counter_S1]->GetParError(3);
      Back_Mean_S1_error[counter_S1] = func1[counter_S1]->GetParError(4);
      Back_Sigma_S1_error[counter_S1] = func1[counter_S1]->GetParError(5);
      Back_Constant_S1_error[counter_S1] = func1[counter_S1]->GetParError(6);
      Back_Gradient_S1_error[counter_S1] = func1[counter_S1]->GetParError(7);

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


   /////////////////////////////////////////////////////////////////////////////////
   ///// Strangeness 2     /////////////////////////////////////////////////////////
   /////////////////////////////////////////////////////////////////////////////////

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
      par0_S2 = 40*exp(-0.5*pow((momentum_S2[counter_S2]-1.2)/0.18,2))+80+-20*momentum_S2[counter_S2];
      par1_S2 = 4.95277e-01 -2.78269e-04*momentum_S2[counter_S2] - 2.19406e-03 * exp(-0.5*pow((momentum_S2[counter_S2]-1.79144e+00) / 2.70089e-01,2)); // mean for both signal gauss
      par2_S2 = exp(-5.87651e+00+9.50189e-01*momentum_S2[counter_S2])+exp(-4.83642e+00+-9.34420e-01*momentum_S2[counter_S2]);
      par3_S2 = -2466.73 + 8254.67*momentum_S2[counter_S2] -8574.25*pow(momentum_S2[counter_S2],2)+3580.97*pow(momentum_S2[counter_S2],3) -511.388*pow(momentum_S2[counter_S2],4); // amplitude for back gaus
      par4_S2 = 0.1396; // mean for back gaus (pion mass)
      // par5_S2 = 0.3; // sigma for back gaus
      par5_S2 = 0.393*exp(-0.5*pow((momentum_S2[counter_S2]-1.46)/0.613,2))+-0.534+0.224*momentum_S2[counter_S2]; // sigma for back gaus
      par6_S2 = -5; // constant for back pol1
      par7_S2 = -10; // slope for back pol1
      par8_S2 = 0.578473*exp(-0.5*pow((momentum_S2[counter_S2]-1.66762)/0.229272,2)) + 1.34606*exp(-0.5*pow((momentum_S2[counter_S2]-2.12122)/0.701563,2));


      // Setting parameters
      func1_S2[counter_S2]->FixParameter(0,par0_S2); // amplitude for 1st signal gaus
      func1_S2[counter_S2]->FixParameter(1,par1_S2); // mean for both signal gaus
      func1_S2[counter_S2]->FixParameter(2,par2_S2); // sigma for 1st signal gaus
      func1_S2[counter_S2]->SetParameter(3,par3_S2); // amplitude for back gaus
      func1_S2[counter_S2]->FixParameter(4,par4_S2); // mean for back gaus (pion mass)
      func1_S2[counter_S2]->SetParameter(5,par5_S2); // sigma for back gaus
      func1_S2[counter_S2]->SetParameter(6,par6_S2); // constant for back pol1
      func1_S2[counter_S2]->SetParameter(7,par7_S2); // slope for back pol1
      func1_S2[counter_S2]->FixParameter(8,par8_S2); // ratio between amplitude 1 and 2 for signal


      // Setting parameter limits
      // func1_S2[counter_S2]->SetParLimits(3,Kaon_Mass_S2_Projections[counter_S2]->GetBinContent(Kaon_Mass_S2_Projections[counter_S2]->FindBin(0.365)) / 5, 4 * Kaon_Mass_S2_Projections[counter_S2]->GetMaximum()); // amplitude for back gaus
      // func1_S2[counter_S2]->SetParLimits(4,par4_S2 * 0.98, par4_S2 * 1.02); // mean for back gaus (pion mass)
      // func1_S2[counter_S2]->SetParLimits(5,0.1, 0.8); // sigma for back gaus
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

   /////////////////////////////////////////////////////////////////////////////////
   ///// Strangeness 3     /////////////////////////////////////////////////////////
   /////////////////////////////////////////////////////////////////////////////////

   for(Int_t i = 1; i < hist_S3->GetNbinsX()+1; i++){

      momentum_mid_S3 = hist_S3->GetXaxis()->GetBinCenter(i);

      // Cutting out pion cross over between 1.1 and 1.35 GeV
      if(momentum_mid_S3 > 1.15 && momentum_mid_S3 < 1.35) continue;

      // Kaons only plotted between 0.9-3 GeV
      if(momentum_mid_S3 < 0.85 || momentum_mid_S3 > 3) continue;


      momentum_S3[counter_S3] = momentum_mid_S3;


      Missing_Mass_S3_Projections[counter_S3] = (TH1F*) hist_S3->ProjectionY("",i,i,0,hist_S3->GetNbinsZ())->Clone();
      Kaon_Mass_S3_Projections[counter_S3] = (TH1F*) hist_S3->ProjectionZ("",i,i,0,hist_S3->GetNbinsY())->Clone();

      // Calculate integrals and ratio for S3
      // S3_Integral[counter_S3] = Kaon_Mass_S3_Projections[counter_S3]->Integral();


      // Setting range of functions
      ostringstream funcname_S3, sig1_S3, sig2_S3, sigtotal_S3, back1_S3, back2_S3, backtotal_S3;
      funcname_S3 <<"func_S3_"<<i;
      sig1_S3 <<"sig1_S3_"<<i;
      sig2_S3 <<"sig2_S3_"<<i;
      back1_S3 <<"back1_S3_"<<i;
      back2_S3 <<"back2_S3_"<<i;

      // Defining functions
      func1_S3[counter_S3] = new TF1(funcname_S3.str().c_str(),"[0] * exp(-0.5*pow((x-[1]) / [2],2)) + [0] * [8] * exp(-0.5 * pow((x-[1]) / (2*[2]),2)) + gaus(3) + [6]*[7] + [7]*x",0.365,0.65);
      signal_1_S3[counter_S3] = new TF1(sig1_S3.str().c_str(),"gaus(0)",0.365,0.65);
      signal_2_S3[counter_S3] = new TF1(sig2_S3.str().c_str(),"[0] * (3.63259e-01 * exp(-0.5*pow((x - 1.82697e+00) / 1.82220e-01,2)) + 1.27857e+00 * exp(-0.5 * pow((x - 2.00491e+00) / 6.26903e-01,2))) * exp(-0.5 * pow((x-[1]) / (2*[2]),2))",0.365,0.65);
      signal_total_S3[counter_S3] = new TF1(sigtotal_S3.str().c_str(),"[0] * exp(-0.5*pow((x-[1]) / [2],2)) + [0] * (3.63259e-01 * exp(-0.5*pow((x - 1.82697e+00) / 1.82220e-01,2)) + 1.27857e+00 * exp(-0.5 * pow((x - 2.00491e+00) / 6.26903e-01,2))) * exp(-0.5 * pow((x-[1]) / (2*[2]),2))",0.365,0.65);
      background_1_S3[counter_S3] = new TF1(back1_S3.str().c_str(),"gaus(0)",0.0,0.65);
      background_2_S3[counter_S3] = new TF1(back2_S3.str().c_str(),"[0]*[1] + [1]*x",0.365,0.65);
      background_total_S3[counter_S3] = new TF1(backtotal_S3.str().c_str(),"gaus(0) + [3]*[4] + [4]*x",0.365,0.65);

      par0_S3 = 43.357*exp(-0.5*pow((momentum_S3[counter_S3]-1.22)/0.13,2))+32.647+-5.934*momentum_S3[counter_S3];
      par1_S3 = 4.95277e-01 -2.78269e-04*momentum_S3[counter_S3] - 2.19406e-03 * exp(-0.5*pow((momentum_S3[counter_S3]-1.79144e+00) / 2.70089e-01,2)); // mean for both signal gauss
      par2_S3 = exp(-5.87651e+00+9.50189e-01*momentum_S3[counter_S3])+exp(-4.83642e+00+-9.34420e-01*momentum_S3[counter_S3]);
      par3_S3 = -2466.73 + 8254.67*momentum_S3[counter_S3] -8574.25*pow(momentum_S3[counter_S3],2)+3580.97*pow(momentum_S3[counter_S3],3) -511.388*pow(momentum_S3[counter_S3],4); // amplitude for back gaus
      par4_S3 = 0.1396; // mean for back gaus (pion mass)
      // par5_S3 = 0.3; // sigma for back gaus
      par5_S3 = 0.393*exp(-0.5*pow((momentum_S3[counter_S3]-1.46)/0.613,2))+-0.534+0.224*momentum_S3[counter_S3]; // sigma for back gaus
      par6_S3 = -5; // constant for back pol1
      par7_S3 = -10; // slope for back pol1
      par8_S3 = 0.578473*exp(-0.5*pow((momentum_S3[counter_S3]-1.66762)/0.229272,2)) + 1.34606*exp(-0.5*pow((momentum_S3[counter_S3]-2.12122)/0.701563,2));

      // Setting parameters
      func1_S3[counter_S3]->FixParameter(0,par0_S3); // amplitude for 1st signal gaus
      func1_S3[counter_S3]->FixParameter(1,par1_S3); // mean for both signal gaus
      func1_S3[counter_S3]->FixParameter(2,par2_S3); // sigma for 1st signal gaus
      func1_S3[counter_S3]->SetParameter(3,par3_S3); // amplitude for back gaus
      func1_S3[counter_S3]->FixParameter(4,par4_S3); // mean for back gaus (pion mass)
      func1_S3[counter_S3]->SetParameter(5,par5_S3); // sigma for back gaus
      func1_S3[counter_S3]->SetParameter(6,par6_S3); // constant for back pol1
      func1_S3[counter_S3]->SetParameter(7,par7_S3); // slope for back pol1
      func1_S3[counter_S3]->FixParameter(8,par8_S3); // ratio between amplitude 1 and 2 for signal


      // Setting parameter limits
      // func1_S3[counter_S3]->SetParLimits(0, 0, 1.2 * Kaon_Mass_S3_Projections[counter_S3]->GetMaximum()); // amplitude for back gaus
      // func1_S3[counter_S3]->SetParLimits(3,Kaon_Mass_S3_Projections[counter_S3]->GetBinContent(Kaon_Mass_S3_Projections[counter_S3]->FindBin(0.365)) / 5, 4 * Kaon_Mass_S3_Projections[counter_S3]->GetMaximum()); // amplitude for back gaus
      // func1_S3[counter_S3]->SetParLimits(4,par4_S3 * 0.98, par4_S3 * 1.02); // mean for back gaus (pion mass)
      // func1_S3[counter_S3]->SetParLimits(5,0.1, 0.8); // sigma for back gaus
      func1_S3[counter_S3]->SetParLimits(6,-1000, -0.7); // constant for back pol1
      func1_S3[counter_S3]->SetParLimits(7,-500, 0); // slope for back pol1


      Kaon_Mass_S3_Projections[counter_S3]->Fit(funcname_S3.str().c_str(),"RBQ");


      signal_1_S3[counter_S3]->SetParameter(0,func1_S3[counter_S3]->GetParameter(0));
      signal_1_S3[counter_S3]->SetParameter(1,func1_S3[counter_S3]->GetParameter(1));
      signal_1_S3[counter_S3]->SetParameter(2,func1_S3[counter_S3]->GetParameter(2));
      signal_2_S3[counter_S3]->SetParameter(0,func1_S3[counter_S3]->GetParameter(0));
      signal_2_S3[counter_S3]->SetParameter(1,func1_S3[counter_S3]->GetParameter(1));
      signal_2_S3[counter_S3]->SetParameter(2,func1_S3[counter_S3]->GetParameter(2));
      signal_total_S3[counter_S3]->SetParameter(0,func1_S3[counter_S3]->GetParameter(0));
      signal_total_S3[counter_S3]->SetParameter(1,func1_S3[counter_S3]->GetParameter(1));
      signal_total_S3[counter_S3]->SetParameter(2,func1_S3[counter_S3]->GetParameter(2));
      background_1_S3[counter_S3]->SetParameter(0,func1_S3[counter_S3]->GetParameter(3));
      background_1_S3[counter_S3]->SetParameter(1,func1_S3[counter_S3]->GetParameter(4));
      background_1_S3[counter_S3]->SetParameter(2,func1_S3[counter_S3]->GetParameter(5));
      background_2_S3[counter_S3]->SetParameter(0,func1_S3[counter_S3]->GetParameter(6));
      background_2_S3[counter_S3]->SetParameter(1,func1_S3[counter_S3]->GetParameter(7));
      background_total_S3[counter_S3]->SetParameter(0,func1_S3[counter_S3]->GetParameter(3));
      background_total_S3[counter_S3]->SetParameter(1,func1_S3[counter_S3]->GetParameter(4));
      background_total_S3[counter_S3]->SetParameter(2,func1_S3[counter_S3]->GetParameter(5));
      background_total_S3[counter_S3]->SetParameter(3,func1_S3[counter_S3]->GetParameter(6));
      background_total_S3[counter_S3]->SetParameter(4,func1_S3[counter_S3]->GetParameter(7));


      // Get the projection fit parameters
      Amp_1_S3[counter_S3] = func1_S3[counter_S3]->GetParameter(0);
      Sigma_1_S3[counter_S3] = func1_S3[counter_S3]->GetParameter(2);
      // Amp_2_S3[counter_S3] = func1_S3[counter_S3]->GetParameter(0) * func1_S3[counter_S3]->GetParameter(8);
      Amp_2_S3[counter_S3] = func1_S3[counter_S3]->GetParameter(8);

      // Get the projection fit errors
      Amp_1_error_S3[counter_S3] = func1_S3[counter_S3]->GetParError(0);
      Sigma_1_error_S3[counter_S3] = func1_S3[counter_S3]->GetParError(2);
      Amp_2_error_S3[counter_S3] = func1_S3[counter_S3]->GetParError(0);



      // h_amp_1_S3->Fill(momentum_S3[counter_S3],Amp_1_S3[counter_S3]);
      // h_amp_1_S3->SetBinError(i, Amp_1_error_S3[counter_S3]);
      // h_amp_2_S3->Fill(momentum_S3[counter_S3],Amp_2_S3[counter_S3]);
      // h_amp_2_S3->SetBinError(i, Amp_2_error_S3[counter_S3]);


      // h_S3_Integral->Fill(momentum_S3[counter_S3], S3_Integral[counter_S3]);

      counter_S3++;

   }

   // auto *h_S1_S2_Ratio = new TH1F();
   // h_S1_S2_Ratio = (TH1F*)h_S1_Integral->Clone("h_S1_S2_Ratio");
   // h_S1_S2_Ratio->Divide(h_S2_Integral);

   // Creating TGraphs for the parameters
   // Strangeness 1
   TGraphErrors* gr1 = new TGraphErrors(counter_S1, momentum_S1, Sigma_1, 0, Sigma_1_error);
   TGraphErrors* gr2 = new TGraphErrors(counter_S1, momentum_S1, Amp_1, 0, Amp_1_error);
   TGraphErrors* gr3 = new TGraphErrors(counter_S1, momentum_S1, Amp_2, 0, Amp_2_error);
   TGraphErrors* gr4 = new TGraphErrors(counter_S1, momentum_S1, amp_ratio, 0, Amp_ratio_error);
   TGraphErrors* grback_amp = new TGraphErrors(counter_S1, momentum_S1, Back_Amp_S1, 0, Back_Amp_S1_error);
   TGraphErrors* grback_mean = new TGraphErrors(counter_S1, momentum_S1, Back_Mean_S1, 0, Back_Mean_S1_error);
   TGraphErrors* grback_sigma = new TGraphErrors(counter_S1, momentum_S1, Back_Sigma_S1, 0, Back_Sigma_S1_error);
   TGraphErrors* grback_constant = new TGraphErrors(counter_S1, momentum_S1, Back_Constant_S1, 0, Back_Constant_S1_error);
   TGraphErrors* grback_gradient = new TGraphErrors(counter_S1, momentum_S1, Back_Gradient_S1, 0, Back_Gradient_S1_error);
   // Strangeness 2
   TGraphErrors* gr5 = new TGraphErrors(counter_S2, momentum_S2, Sigma_1_S2, 0, Sigma_1_error_S2);
   TGraphErrors* gr6 = new TGraphErrors(counter_S2, momentum_S2, Amp_1_S2, 0, Amp_1_error_S2);
   TGraphErrors* gr7 = new TGraphErrors(counter_S2, momentum_S2, Amp_2_S2, 0, Amp_2_error_S2);
   // // Strangeness 3
   TGraphErrors* gr8 = new TGraphErrors(counter_S3, momentum_S3, Sigma_1_S3, 0, Sigma_1_error_S3);
   TGraphErrors* gr9 = new TGraphErrors(counter_S3, momentum_S3, Amp_1_S3, 0, Amp_1_error_S3);
   TGraphErrors* gr10 = new TGraphErrors(counter_S3, momentum_S3, Amp_2_S3, 0, Amp_2_error_S3);

   // Creating fits for parameters
   // Strageness 1
   TF1 *func_sig_1_amp = new TF1("func_sig_1_amp","gaus(0) + pol1(3)",0.9,2.6);
   TF1 *func_sig_1_mean = new TF1("func_sig_1_mean","pol1(0) - gaus(2)",0.9,2.6);
   TF1 *func_sig_1_sigma = new TF1("func_sig_1_sigma","expo(0) + expo(2)",0.9,2.6);
   TF1 *func_sig_2_amp = new TF1("func_sig_2_amp","gaus(0) + gaus(3)",0.9,2.6);
   // Strageness 2
   TF1 *func_sig_1_amp_S2 = new TF1("func_sig_1_amp_S2","gaus(0) + pol1(3)",0.9,2.6);
   TF1 *func_sig_1_mean_S2 = new TF1("func_sig_1_mean_S2","pol1(0) - gaus(2)",0.9,2.6);
   TF1 *func_sig_1_sigma_S2 = new TF1("func_sig_1_sigma_S2","expo(0) + expo(2)",0.9,2.6);
   TF1 *func_sig_2_amp_S2 = new TF1("func_sig_2_amp_S2","gaus(0) + gaus(3)",0.9,2.6);
   // // Strageness 3
   TF1 *func_sig_1_amp_S3 = new TF1("func_sig_1_amp_S3","gaus(0) + pol1(3)",0.9,2.6);
   TF1 *func_sig_1_mean_S3 = new TF1("func_sig_1_mean_S3","pol1(0) - gaus(2)",0.9,2.6);
   TF1 *func_sig_1_sigma_S3 = new TF1("func_sig_1_sigma_S3","expo(0) + expo(2)",0.9,2.6);
   TF1 *func_sig_2_amp_S3 = new TF1("func_sig_2_amp_S3","gaus(0) + gaus(3)",0.9,2.6);

   // Setting parameters for fitting parameters


   // signal 1 amplitude
   // gaus + pol1
   func_sig_1_amp->FixParameter(0,1090);
   func_sig_1_amp->FixParameter(1,1.052);
   func_sig_1_amp->FixParameter(2,0.272);
   func_sig_1_amp->FixParameter(3,1299);
   func_sig_1_amp->FixParameter(4,-322.3);

   // signal 1 mean
   // pol1 - gaus
   func_sig_1_mean->FixParameter(0,4.95277e-01);
   func_sig_1_mean->FixParameter(1,-2.78269e-04);
   func_sig_1_mean->FixParameter(2,2.19406e-03);
   func_sig_1_mean->FixParameter(3,1.79144e+00);
   func_sig_1_mean->FixParameter(4,2.70089e-01);

   // signal 1 sigma
   // 2 exponentials
   func_sig_1_sigma->FixParameter(0,-5.87651e+00);
   func_sig_1_sigma->FixParameter(1,9.50189e-01);
   func_sig_1_sigma->FixParameter(2,-4.83642e+00);
   func_sig_1_sigma->FixParameter(3,-9.34420e-01);

   // signal 2 amplitude ratio to signal 1 amplitude
   func_sig_2_amp->FixParameter(0,0.578473);
   func_sig_2_amp->FixParameter(1,1.66762);
   func_sig_2_amp->FixParameter(2,0.229272);
   func_sig_2_amp->FixParameter(3,1.34606);
   func_sig_2_amp->FixParameter(4,2.12122);
   func_sig_2_amp->FixParameter(5,0.701563);

   // Strangeness 2
   // signal 1 amplitude
   // gaus + pol1
   func_sig_1_amp_S2->SetParameter(0,40);
   func_sig_1_amp_S2->SetParameter(1,1.2);
   func_sig_1_amp_S2->SetParameter(2,0.18);
   func_sig_1_amp_S2->SetParameter(3,80);
   func_sig_1_amp_S2->SetParameter(4,-20);
   //
   //
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

   // signal 2 amplitude ratio to signal 1 amplitude
   func_sig_2_amp_S2->FixParameter(0,0.578473);
   func_sig_2_amp_S2->FixParameter(1,1.66762);
   func_sig_2_amp_S2->FixParameter(2,0.229272);
   func_sig_2_amp_S2->FixParameter(3,1.34606);
   func_sig_2_amp_S2->FixParameter(4,2.12122);
   func_sig_2_amp_S2->FixParameter(5,0.701563);

   // Strangeness 3
   // signal 1 amplitude
   // gaus + pol1
   func_sig_1_amp_S3->FixParameter(0,43.357);
   func_sig_1_amp_S3->FixParameter(1,1.22);
   func_sig_1_amp_S3->FixParameter(2,0.13);
   func_sig_1_amp_S3->FixParameter(3,32.647);
   func_sig_1_amp_S3->FixParameter(4,-5.934);


   // signal 1 mean
   // pol1 - gaus
   func_sig_1_mean_S3->FixParameter(0,4.95277e-01);
   func_sig_1_mean_S3->FixParameter(1,-2.78269e-04);
   func_sig_1_mean_S3->FixParameter(2,2.19406e-03);
   func_sig_1_mean_S3->FixParameter(3,1.79144e+00);
   func_sig_1_mean_S3->FixParameter(4,2.70089e-01);

   // signal 1 sigma
   // 2 exponentials
   func_sig_1_sigma_S3->FixParameter(0,-5.87651e+00);
   func_sig_1_sigma_S3->FixParameter(1,9.50189e-01);
   func_sig_1_sigma_S3->FixParameter(2,-4.83642e+00);
   func_sig_1_sigma_S3->FixParameter(3,-9.34420e-01);

   // signal 2 amplitude ratio to signal 1 amplitude
   func_sig_2_amp_S3->FixParameter(0,0.578473);
   func_sig_2_amp_S3->FixParameter(1,1.66762);
   func_sig_2_amp_S3->FixParameter(2,0.229272);
   func_sig_2_amp_S3->FixParameter(3,1.34606);
   func_sig_2_amp_S3->FixParameter(4,2.12122);
   func_sig_2_amp_S3->FixParameter(5,0.701563);

   // Fitting parameter functions
   gr1->Fit("func_sig_1_sigma","RBQ");
   gr2->Fit("func_sig_1_amp","RBQ");
   gr3->Fit("func_sig_2_amp","RBQ");
   gr5->Fit("func_sig_1_sigma_S2","RBQ");
   gr6->Fit("func_sig_1_amp_S2","RBQ");
   gr7->Fit("func_sig_2_amp_S2","RBQ");
   gr8->Fit("func_sig_1_sigma_S3","RBQ");
   gr9->Fit("func_sig_1_amp_S3","RBQ");
   gr10->Fit("func_sig_2_amp_S3","RBQ");

   // Creating strings for parameter values
   ostringstream signal_1_amp, signal_mean, signal_1_sigma, signal_2_amp, signal_2_sigma;
   ostringstream signal_1_amp_S2, signal_mean_S2, signal_1_sigma_S2, signal_2_amp_S2, signal_2_sigma_S2;

   TH2F *Signal_Function = new TH2F("Signal_Function","",500,0.3,0.8,300,0,3);
   TH2F *Background_Function = new TH2F("Background_Function","",500,0.3,0.8,300,0,3);
   TF1 *signal_function = new TF1("signal_function","gaus(0) + gaus(3)",0.3,0.8);
   TH2F *Signal_Function_S2 = new TH2F("Signal_Function_S2","",500/S2_z_rebin,0.3,0.8,300/S2_x_rebin,0,3);
   TH2F *Background_Function_S2 = new TH2F("Background_Function_S2","",500/S2_z_rebin,0.3,0.8,300/S2_x_rebin,0,3);
   TF1 *signal_function_S2 = new TF1("signal_function_S2","gaus(0) + gaus(3)",0.3,0.8);
   TH2F *Signal_Function_S3 = new TH2F("Signal_Function_S3","",500/S3_z_rebin,0.3,0.8,300/S3_x_rebin,0,3);
   TH2F *Background_Function_S3 = new TH2F("Background_Function_S3","",500/S3_z_rebin,0.3,0.8,300/S3_x_rebin,0,3);
   TF1 *signal_function_S3 = new TF1("signal_function_S3","gaus(0) + gaus(3)",0.3,0.8);
   TH2F *Result = new TH2F();
   TH2F *Result_S2 = new TH2F();
   TH2F *Result_S3 = new TH2F();

   // Take copy of the data kaon momentum vs mass
   TH2F *Kaon_Mass_Momentum = new TH2F();
   Kaon_Mass_Momentum = (TH2F*)hist_S1->Project3D("xz")->Clone();

   TH2F *Kaon_Mass_Momentum_S2 = new TH2F();
   Kaon_Mass_Momentum_S2 = (TH2F*)hist_S2->Project3D("xz")->Clone();

   TH2F *Kaon_Mass_Momentum_S3 = new TH2F();
   Kaon_Mass_Momentum_S3 = (TH2F*)hist_S3->Project3D("xz")->Clone();


   TH3F *S1_Signal = new TH3F("S1_Signal","Sig",300,0,3,300,0,3,500,0.3,0.8);
   TH3F *S1_Background = new TH3F("S1_Background","Back",300,0,3,300,0,3,500,0.3,0.8);
   TH3F *S2_Signal = new TH3F("S2_Signal","Sig2",300/S2_x_rebin,0,3,300/S2_y_rebin,0,3,500/S2_z_rebin,0.3,0.8);
   TH3F *S2_Background = new TH3F("S2_Background","Back2",300/S2_x_rebin,0,3,300/S2_y_rebin,0,3,500/S2_z_rebin,0.3,0.8);
   TH3F *S3_Signal = new TH3F("S3_Signal","Sig3",300/S3_x_rebin,0,3,300/S3_y_rebin,0,3,500/S3_z_rebin,0.3,0.8);
   TH3F *S3_Background = new TH3F("S3_Background","Back3",300/S3_x_rebin,0,3,300/S3_y_rebin,0,3,500/S3_z_rebin,0.3,0.8);


   // Creating signal function and histogram

   // Strangeness 1
   // Loop over kaon momentum
   for(Int_t x_pos=1; x_pos < hist_S1->GetNbinsX() + 1; x_pos++){
   // for(Int_t x_pos=101; x_pos < 112; x_pos++){

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
         if(Kaon_Mass_Momentum->GetBinContent(z_pos,x_pos) > 0){
            if((signal_function->Eval(Signal_Function->GetXaxis()->GetBinCenter(z_pos))/Kaon_Mass_Momentum->GetBinContent(z_pos,x_pos)) < 1){
               Signal_Function->SetBinContent(z_pos,x_pos,signal_function->Eval(Signal_Function->GetXaxis()->GetBinCenter(z_pos))/Kaon_Mass_Momentum->GetBinContent(z_pos,x_pos));
            }
            else{
               Signal_Function->SetBinContent(z_pos,x_pos,1);

            }
            // Make sure there are no negative values for the background
            Background_Function->SetBinContent(z_pos,x_pos,1 - Signal_Function->GetBinContent(z_pos,x_pos));
         }

         // Loop over missing mass
         for(Int_t y_pos=1; y_pos < hist_S1->GetNbinsY() + 1; y_pos++){

            if(hist_S1->GetBinContent(x_pos,y_pos,z_pos) > 0){

               S1_Signal->SetBinContent(x_pos, y_pos, z_pos, Signal_Function->GetBinContent(z_pos,x_pos) * hist_S1->GetBinContent(x_pos,y_pos,z_pos));
               S1_Background->SetBinContent(x_pos, y_pos, z_pos, Background_Function->GetBinContent(z_pos,x_pos) * hist_S1->GetBinContent(x_pos,y_pos,z_pos));

            }
            else{
               S1_Signal->SetBinContent(x_pos, y_pos, z_pos, 0);
               S1_Background->SetBinContent(x_pos, y_pos, z_pos, 0);


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
         if(Kaon_Mass_Momentum_S2->GetBinContent(z_pos,x_pos) > 0){
            if((signal_function_S2->Eval(Signal_Function_S2->GetXaxis()->GetBinCenter(z_pos))/Kaon_Mass_Momentum_S2->GetBinContent(z_pos,x_pos)) < 1){
               Signal_Function_S2->SetBinContent(z_pos,x_pos,signal_function_S2->Eval(Signal_Function_S2->GetXaxis()->GetBinCenter(z_pos))/Kaon_Mass_Momentum_S2->GetBinContent(z_pos,x_pos));
            }
            else{
               Signal_Function_S2->SetBinContent(z_pos,x_pos,1);

            }
            Background_Function_S2->SetBinContent(z_pos,x_pos,1 - Signal_Function_S2->GetBinContent(z_pos,x_pos));
         }

         // Loop over missing mass
         for(Int_t y_pos=1; y_pos < hist_S2->GetNbinsY() + 1; y_pos++){

            if(hist_S2->GetBinContent(x_pos,y_pos,z_pos) > 0){

               S2_Signal->SetBinContent(x_pos, y_pos, z_pos, Signal_Function_S2->GetBinContent(z_pos,x_pos) * hist_S2->GetBinContent(x_pos,y_pos,z_pos));
               S2_Background->SetBinContent(x_pos, y_pos, z_pos, Background_Function_S2->GetBinContent(z_pos,x_pos) * hist_S2->GetBinContent(x_pos,y_pos,z_pos));
            }
            else{
               S2_Signal->SetBinContent(x_pos, y_pos, z_pos, 0);
               S2_Background->SetBinContent(x_pos, y_pos, z_pos, 0);


            }
         }
      }
   }

   // Strangeness 3
   // Loop over kaon momentum
   for(Int_t x_pos=1; x_pos < hist_S3->GetNbinsX() + 1; x_pos++){

      // Set signal function parameters based on values obtained from the 1D fit
      // parameter functions for current bin in momentum
      signal_function_S3->SetParameter(0,func_sig_1_amp_S3->Eval(Signal_Function_S3->GetYaxis()->GetBinCenter(x_pos)));
      signal_function_S3->SetParameter(1,func_sig_1_mean_S3->Eval(Signal_Function_S3->GetYaxis()->GetBinCenter(x_pos)));
      signal_function_S3->SetParameter(2,func_sig_1_sigma_S3->Eval(Signal_Function_S3->GetYaxis()->GetBinCenter(x_pos)));
      signal_function_S3->SetParameter(3,func_sig_1_amp_S3->Eval(Signal_Function_S3->GetYaxis()->GetBinCenter(x_pos)) * func_sig_2_amp_S3->Eval(Signal_Function_S3->GetYaxis()->GetBinCenter(x_pos)));
      signal_function_S3->SetParameter(4,func_sig_1_mean_S3->Eval(Signal_Function_S3->GetYaxis()->GetBinCenter(x_pos)));
      signal_function_S3->SetParameter(5,2 * func_sig_1_sigma_S3->Eval(Signal_Function_S3->GetYaxis()->GetBinCenter(x_pos)));


      // Loop over kaon mass
      for(Int_t z_pos=1; z_pos < hist_S3->GetNbinsZ() + 1; z_pos++){

         // Set bin content for 2D signal histogram based on current momentum
         // and kaon mass bin
         if(Kaon_Mass_Momentum_S3->GetBinContent(z_pos,x_pos) > 0){
            if((signal_function_S3->Eval(Signal_Function_S3->GetXaxis()->GetBinCenter(z_pos))/Kaon_Mass_Momentum_S3->GetBinContent(z_pos,x_pos)) < 1){
               Signal_Function_S3->SetBinContent(z_pos,x_pos,signal_function_S3->Eval(Signal_Function_S3->GetXaxis()->GetBinCenter(z_pos))/Kaon_Mass_Momentum_S3->GetBinContent(z_pos,x_pos));
            }
            else{
               Signal_Function_S3->SetBinContent(z_pos,x_pos,1);

            }
            Background_Function_S3->SetBinContent(z_pos,x_pos,1 - Signal_Function_S3->GetBinContent(z_pos,x_pos));
         }

         // Loop over missing mass
         for(Int_t y_pos=1; y_pos < hist_S3->GetNbinsY() + 1; y_pos++){

            if(hist_S3->GetBinContent(x_pos,y_pos,z_pos) > 0){

               S3_Signal->SetBinContent(x_pos, y_pos, z_pos, Signal_Function_S3->GetBinContent(z_pos,x_pos) * hist_S3->GetBinContent(x_pos,y_pos,z_pos));
               S3_Background->SetBinContent(x_pos, y_pos, z_pos, Background_Function_S3->GetBinContent(z_pos,x_pos) * hist_S3->GetBinContent(x_pos,y_pos,z_pos));
            }
            else{
               S3_Signal->SetBinContent(x_pos, y_pos, z_pos, 0);
               S3_Background->SetBinContent(x_pos, y_pos, z_pos, 0);


            }
         }
      }
   }

   hist_S1->ProjectionY("S1data",1,hist_S1->GetNbinsX(),1,hist_S1->GetNbinsZ(),"");
   S1_Signal->ProjectionY("S1sig",1,hist_S1->GetNbinsX(),1,hist_S1->GetNbinsZ(),"");
   S1_Background->ProjectionY("S1back",1,hist_S1->GetNbinsX(),1,hist_S1->GetNbinsZ(),"");

   S1back->Scale(S1data->Integral(S1data->FindBin(0),S1data->FindBin(0.85)) / S1back->Integral(S1back->FindBin(0),S1back->FindBin(0.85)));
   TH1F *S1Result=(TH1F*)S1data->Clone("S1Result");
   S1Result->Add(S1back,-1);


   hist_S2->ProjectionY("S2data",1,hist_S2->GetNbinsX(),1,hist_S2->GetNbinsZ(),"");
   S2_Signal->ProjectionY("S2sig",1,hist_S2->GetNbinsX(),1,hist_S2->GetNbinsZ(),"");
   S2_Background->ProjectionY("S2back",1,hist_S2->GetNbinsX(),1,hist_S2->GetNbinsZ(),"");

   S2back->Scale(S2data->Integral(S2data->FindBin(0),S2data->FindBin(0.85)) / S2back->Integral(S2back->FindBin(0),S2back->FindBin(0.85)));
   TH1F *S2Result=(TH1F*)S2data->Clone("S2Result");
   S2Result->Add(S2back,-1);

   hist_S3->ProjectionY("S3data",1,hist_S3->GetNbinsX(),1,hist_S3->GetNbinsZ(),"");
   S3_Signal->ProjectionY("S3sig",1,hist_S3->GetNbinsX(),1,hist_S3->GetNbinsZ(),"");
   S3_Background->ProjectionY("S3back",1,hist_S3->GetNbinsX(),1,hist_S3->GetNbinsZ(),"");

   S3back->Scale(S3data->Integral(S3data->FindBin(0),S3data->FindBin(1.55)) / S3back->Integral(S3back->FindBin(0),S3back->FindBin(1.55)));
   TH1F *S3Result=(TH1F*)S3data->Clone("S3Result");
   S3Result->Add(S3back,-1);


   fileOutput1.Write();

}
