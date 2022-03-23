// Making lots of projections in momentum

{
   // Get input file
   TFile *f1 = new TFile("/media/mn688/Elements1/PhD/Analysis_Output/Hexaquark/RGA_Fall2018_Outbending_at_least_1e1KpFD_Tree_Total_24022022_Total_Scaling_16032022_01.root");
   // Get histogram from input file
   TH3F* hist_S1=(TH3F*)f1->Get("h_S1_Kaon_Momentum__Miss_Mass__Kaon_Mass");
   TH3F* hist_S2=(TH3F*)f1->Get("h_S2_Kaon_Momentum__Miss_Mass__Kaon_Mass");

   // histograms looking at fit parameters as a function of momentum
   TH1F *h_sigma_1 = new TH1F("h_sigma_1","sigma 1 relationship",160,1,2.6);
   TH1F *h_mean_signal = new TH1F("h_mean_signal","signal mean relationship",160,1,2.6);
   TH1F *h_amp_1 = new TH1F("h_amp_1","amplitude 1 relationship",160,1,2.6);
   TH1F *h_sigma_factor = new TH1F("h_sigma_factor","sigma 2 relationship to sigma 1",160,1,2.6);
   TH1F *h_amp_2 = new TH1F("h_amp_2","amplitude 2 relationship",160,1,2.6);
   TH1F *h_sigma_3 = new TH1F("h_sigma_3","sigma 3 relationship",160,1,2.6);
   TH1F *h_chi2 = new TH1F("h_chi2","chi^{2} relationship",160,1,2.6);

   TH1F *h_S1_Integral = new TH1F("h_S1_Integral","S1 integral",160,1,2.6);
   TH1F *h_S2_Integral = new TH1F("h_S2_Integral","S2 integral",160,1,2.6);
   // TH1F *h_S1_S2_Ratio = new TH1F("h_S1_S2_Ratio","Ratio between S1 and S2",160,1,2.6);

   // 2d histograms for the signal subtraction
   TH2F *h_signal_2d = new TH2F("h_signal_2d","2D signal",500,0.3,0.8,300,0,3);


   // Changing loop to momentum bin values
   Double_t bin_1[160], bin_2[160], momentum[135], momentum_mid;
   // Ratio between S1 and S2
   Double_t S1_Integral[135], S2_Integral[135];
   // Reduced chi^2
   Double_t Chi2;
   // Parameter values
   Double_t Sigma_1[135], Amp_1[135], Amp_2[135], amp_ratio[135];
   // Parameter errors
   Double_t Sigma_1_error[135], Amp_1_error[135], Amp_2_error[135], momentum_error[135], Amp_ratio_error[135];
   // Function range
   Double_t Range_min, Range_max;

   // TF1 *func1 = new TF1(funcname.str().c_str(),"gaus(0) + pol3(3)",0.365,0.6);
   TH1F *Missing_Mass_S1_Projections[135];
   TH1F *Kaon_Mass_S1_Projections[135];
   TH1F *Kaon_Mass_S2_Projections[135];

   // Create arrays of the various functions
   TF1 *func1[135];
   TF1 *signal_1[135];
   TF1 *signal_2[135];
   TF1 *signal_total[135];
   TF1 *background_1[135];
   TF1 *background_2[135];
   TF1 *background_total[135];

   Int_t counter = 0;

   auto *c1 = new TCanvas("c1","canvas",800,800);
   // Parameter values to be set
   Double_t par0, par1, par2, par3, par4, par5, par6, par7, par8, par9;
   for(Int_t i = 0; i < 160; i++){

      momentum_mid = 1 + ((i + 0.5) * 0.01);


      // Cutting out pion cross over between 1.1 and 1.35 GeV
      if(momentum_mid > 1.1 && momentum_mid < 1.35) continue;

      // Kaons only plotted up to 3 GeV
      if(momentum_mid > 3) continue;


      momentum[counter] = momentum_mid;



      Missing_Mass_S1_Projections[counter] = (TH1F*) hist_S1->ProjectionY("",hist_S1->GetXaxis()->FindBin(momentum_mid),hist_S1->GetXaxis()->FindBin(momentum_mid),0,hist_S1->GetNbinsZ())->Clone();
      Kaon_Mass_S1_Projections[counter] = (TH1F*) hist_S1->ProjectionZ("",hist_S1->GetXaxis()->FindBin(momentum_mid),hist_S1->GetXaxis()->FindBin(momentum_mid),0,hist_S1->GetNbinsY())->Clone();
      Kaon_Mass_S2_Projections[counter] = (TH1F*) hist_S2->ProjectionZ("",hist_S2->GetXaxis()->FindBin(momentum_mid),hist_S2->GetXaxis()->FindBin(momentum_mid),0,hist_S2->GetNbinsY())->Clone();

      // Calculate integrals and ratio for S1 and S2
      S1_Integral[counter] = Kaon_Mass_S1_Projections[counter]->Integral();
      S2_Integral[counter] = Kaon_Mass_S2_Projections[counter]->Integral();

      // cout<< "Integral ratio " << S1_S2_Ratio[counter] << endl;

      // Rebin projections
      // Strangeness 2 looks good with factor of 2
      // Kaon_Mass_S1_Projections[counter]->Rebin(2);

      // Setting range of functions
      ostringstream funcname, sig1, sig2, sigtotal, back1, back2, backtotal;
      funcname <<"func_"<<i;
      sig1 <<"sig1_"<<i;
      sig2 <<"sig2_"<<i;
      back1 <<"back1_"<<i;
      back2 <<"back2_"<<i;
      func1[counter] = new TF1(funcname.str().c_str(),"[0] * exp(-pow(x-[1],2) / (2 * pow([2],2))) + [3] * exp(-pow(x-[1],2) / (2 * pow([2]*[4],2))) + gaus(5) + [8]*[9] + [9]*x",0.365,0.65);
      signal_1[counter] = new TF1(sig1.str().c_str(),"gaus(0)",0.365,0.65);
      signal_2[counter] = new TF1(sig2.str().c_str(),"gaus(0)",0.365,0.65);
      signal_total[counter] = new TF1(sigtotal.str().c_str(),"gaus(0) + gaus(3)",0.365,0.65);
      background_1[counter] = new TF1(back1.str().c_str(),"gaus(0)",0.0,0.65);
      background_2[counter] = new TF1(back2.str().c_str(),"pol1(0)",0.365,0.65);
      background_total[counter] = new TF1(backtotal.str().c_str(),"gaus(0) + pol1(3)",0.365,0.65);



      // Setting parameters based on functions determined
      // par0 = Kaon_Mass_S1_Projections[counter]->GetMaximum() / 2; // amplitude for 1st signal gaus
      // // par0 = 0.08024 * pow(i,2) - 20.9076 * i + 1827.356; // amplitude for 1st signal gaus
      // par1 = 0.493677; // mean for both signal gauss
      // // par2 = 0.000000857 * pow(i+1,2) + 0.00001848 * (i+1) + 0.01094; // sigma for 1st signal gaus
      // par2 = 0.00000135 * pow(i+1,2) - 0.00005776 * (i+1) + 0.01323; // sigma for 1st signal gaus
      // par3 = Kaon_Mass_S1_Projections[counter]->GetMaximum() / 3; // amplitude for 2nd signal gaus
      // par4 = 2; // ratio between signal gaussian sigmas
      // // par4 = -0.00001617 * pow(i+1,2) + 0.001323 * (i+1) + 1.928; // ratio between signal gaussian sigmas
      // par5 = 2 * Kaon_Mass_S1_Projections[counter]->GetBinContent(Kaon_Mass_S1_Projections[counter]->FindBin(0.365)); // amplitude for back gaus
      // par6 = 0.1396; // mean for back gaus (pion mass)
      // par7 = 0.3; // sigma for back gaus
      // par8 = -5; // constant for back pol1
      // par9 = -10; // slope for back pol1

      par0 = 550*exp(-0.5*((momentum[counter]-1.1418174)/0.32175315)*((momentum[counter]-1.1418174)/0.32175315))+(296.52772+-27.580470*momentum[counter]);
      par1 = 4.95277e-01 -2.78269e-04*momentum[counter] - (2.19406e-03 * exp(-pow(momentum[counter]-1.79144e+00,2) / 2 /pow(2.70089e-01,2))); // mean for both signal gauss
      // par1 = 0.493677; // mean for both signal gauss
      par2 = 	exp(-5.87651e+00+9.50189e-01*momentum[counter])+exp(-4.83642e+00+-9.34420e-01*momentum[counter]);
      par3 = par0 * (3.63259e-01 * exp(-pow(momentum[counter] - 1.82697e+00,2) / 2 / pow(1.82220e-01,2)) + 1.27857e+00 * exp(-pow(momentum[counter] - 2.00491e+00,2) / 2 / pow(6.26903e-01,2))); // amplitude for 2nd signal gaus
      par4 = 2;
      par5 = 2 * Kaon_Mass_S1_Projections[counter]->GetBinContent(Kaon_Mass_S1_Projections[counter]->FindBin(0.365)); // amplitude for back gaus
      par6 = 0.1396; // mean for back gaus (pion mass)
      par7 = 0.3; // sigma for back gaus
      par8 = -5; // constant for back pol1
      par9 = -10; // slope for back pol1


      // Setting parameters
      func1[counter]->SetParameter(0,par0); // amplitude for 1st signal gaus
      func1[counter]->FixParameter(1,par1); // mean for both signal gaus
      func1[counter]->FixParameter(2,par2); // sigma for 1st signal gaus
      func1[counter]->FixParameter(3,par3); // amplitude for 2nd signal gaus
      func1[counter]->FixParameter(4,par4); // ratio between signal gaussian sigmas
      func1[counter]->SetParameter(5,par5); // amplitude for back gaus
      func1[counter]->SetParameter(6,par6); // mean for back gaus (pion mass)
      func1[counter]->SetParameter(7,par7); // sigma for back gaus
      func1[counter]->SetParameter(8,par8); // constant for back pol1
      func1[counter]->SetParameter(9,par9); // slope for back pol1

      // Setting parameter limits
      // func1[counter]->SetParLimits(0,Kaon_Mass_S1_Projections[counter]->GetMaximum() / 30, Kaon_Mass_S1_Projections[counter]->GetMaximum()); // amplitude for 1st signal gaus
      func1[counter]->SetParLimits(0,par0 * 0.7, par0 * 1.3); // amplitude for 1st signal gaus
      // func1[counter]->SetParLimits(1,par1 * 0.98, par1 * 1.02); // mean for both signal gaus
      // func1[counter]->SetParLimits(2,0.8 * par2, 1.2 * par2); // sigma for 1st signal gaus
      func1[counter]->SetParLimits(3,Kaon_Mass_S1_Projections[counter]->GetMaximum() / 30, Kaon_Mass_S1_Projections[counter]->GetMaximum()); // amplitude for 2nd signal gaus
      // func1[counter]->SetParLimits(4,0.85 * par4, 1.15 * par4); // ratio between signal gaussian sigmas
      func1[counter]->SetParLimits(5,Kaon_Mass_S1_Projections[counter]->GetBinContent(Kaon_Mass_S1_Projections[counter]->FindBin(0.365)) / 5, 4 * Kaon_Mass_S1_Projections[counter]->GetMaximum()); // amplitude for back gaus
      func1[counter]->SetParLimits(6,par6 * 0.98, par6 * 1.02); // mean for back gaus (pion mass)
      func1[counter]->SetParLimits(7,0.1, 0.8); // sigma for back gaus
      func1[counter]->SetParLimits(8,-1000, -0.7); // constant for back pol1
      func1[counter]->SetParLimits(9,-500, 0); // slope for back pol1


      Kaon_Mass_S1_Projections[counter]->Fit(funcname.str().c_str(),"RBQ");

      Chi2 = func1[counter]->GetChisquare() / func1[counter]->GetNDF();


      signal_1[counter]->SetParameter(0,func1[counter]->GetParameter(0));
      signal_1[counter]->SetParameter(1,func1[counter]->GetParameter(1));
      signal_1[counter]->SetParameter(2,func1[counter]->GetParameter(2));
      signal_2[counter]->SetParameter(0,func1[counter]->GetParameter(3));
      signal_2[counter]->SetParameter(1,func1[counter]->GetParameter(1));
      signal_2[counter]->SetParameter(2,func1[counter]->GetParameter(2) * func1[counter]->GetParameter(4));
      signal_total[counter]->SetParameter(0,func1[counter]->GetParameter(0));
      signal_total[counter]->SetParameter(1,func1[counter]->GetParameter(1));
      signal_total[counter]->SetParameter(2,func1[counter]->GetParameter(2));
      signal_total[counter]->SetParameter(3,func1[counter]->GetParameter(3));
      signal_total[counter]->SetParameter(4,func1[counter]->GetParameter(1));
      signal_total[counter]->SetParameter(5,func1[counter]->GetParameter(2) * func1[counter]->GetParameter(4));
      background_1[counter]->SetParameter(0,func1[counter]->GetParameter(5));
      background_1[counter]->SetParameter(1,func1[counter]->GetParameter(6));
      background_1[counter]->SetParameter(2,func1[counter]->GetParameter(7));
      background_2[counter]->SetParameter(0,func1[counter]->GetParameter(8) * func1[counter]->GetParameter(9));
      background_2[counter]->SetParameter(1,func1[counter]->GetParameter(9));
      background_total[counter]->SetParameter(0,func1[counter]->GetParameter(5));
      background_total[counter]->SetParameter(1,func1[counter]->GetParameter(6));
      background_total[counter]->SetParameter(2,func1[counter]->GetParameter(7));
      background_total[counter]->SetParameter(3,func1[counter]->GetParameter(8));
      background_total[counter]->SetParameter(4,func1[counter]->GetParameter(9));

      // Get the projection fit parameters
      Amp_1[counter] = func1[counter]->GetParameter(0);
      Sigma_1[counter] = func1[counter]->GetParameter(2);
      Amp_2[counter] = func1[counter]->GetParameter(3);
      amp_ratio[counter] = func1[counter]->GetParameter(3) / func1[counter]->GetParameter(0);

      // Get the projection fit errors
      Sigma_1_error[counter] = func1[counter]->GetParError(2);
      Amp_1_error[counter] = func1[counter]->GetParError(0);
      Amp_2_error[counter] = func1[counter]->GetParError(3);
      Amp_ratio_error[counter] = 0.2;
      // momentum_error[counter] = sqrt(Kaon_Mass_S1_Projections[counter]->Integral()) / Kaon_Mass_S1_Projections[counter]->Integral();

      h_chi2->Fill(momentum[counter], Chi2);
      h_sigma_1->Fill(momentum[counter],Sigma_1[counter]);
      h_sigma_1->SetBinError(i, Sigma_1_error[counter]);
      h_amp_1->Fill(momentum[counter],Amp_1[counter]);
      h_amp_1->SetBinError(i, Amp_1_error[counter]);
      h_amp_2->Fill(momentum[counter],Amp_2[counter]);
      h_amp_2->SetBinError(i, Amp_2_error[counter]);
      h_mean_signal->Fill(momentum[counter],func1[counter]->GetParameter(1));
      h_mean_signal->SetBinError(momentum[counter],func1[counter]->GetParError(1));

      h_S1_Integral->Fill(momentum[counter], S1_Integral[counter]);
      h_S2_Integral->Fill(momentum[counter], S2_Integral[counter]);


      counter++;

   }

   h_S1_Integral->Sumw2();
   // h_S2_Integral->Sumw2();
   auto *h_S1_S2_Ratio = new TH1F();
   h_S1_S2_Ratio = (TH1F*)h_S1_Integral->Clone("h_S1_S2_Ratio");
   h_S1_S2_Ratio->Divide(h_S2_Integral);

   // Creating TGraphs for the parameters
   TGraphErrors* gr1 = new TGraphErrors(135, momentum, Sigma_1, 0, Sigma_1_error);
   TGraphErrors* gr2 = new TGraphErrors(135, momentum, Amp_1, 0, Amp_1_error);
   TGraphErrors* gr3 = new TGraphErrors(135, momentum, Amp_2, 0, Amp_2_error);
   TGraphErrors* gr4 = new TGraphErrors(135, momentum, amp_ratio, 0, Amp_ratio_error);

   // Creating fits for parameters
   TF1 *func_sig_1_amp = new TF1("func_sig_1_amp","gaus(0) + pol1(3)",1.0,2.6);
   TF1 *func_sig_1_mean = new TF1("func_sig_1_mean","pol1(0) - gaus(2)",1.0,2.6);
   TF1 *func_sig_1_sigma = new TF1("func_sig_1_sigma","expo(0) + expo(2)",1.0,2.6);
   TF1 *func_sig_2_amp = new TF1("func_sig_2_amp","gaus(0) + gaus(3)",1.0,2.6);


   // Defining parameters
   par0 = 550*exp(-0.5*pow((momentum[counter]-1.1418174) / 0.32175315,2)) +(296.52772+-27.580470*momentum[counter]);
   par1 = 4.95277e-01 -2.78269e-04*momentum[counter] - (2.19406e-03 * exp(-0.5*pow((momentum[counter]-1.79144e+00) / 2.70089e-01,2))); // mean for both signal gauss
   // par1 = 0.493677; // mean for both signal gauss
   par2 = 	exp(-5.87651e+00+9.50189e-01*momentum[counter]) + exp(-4.83642e+00+-9.34420e-01*momentum[counter]);
   par3 = par0 * (3.63259e-01 * exp(-0.5*pow((momentum[counter] - 1.82697e+00) / 1.82220e-01,2)) + 1.27857e+00 * exp(-0.5 * pow((momentum[counter] - 2.00491e+00) / 6.26903e-01,2))); // amplitude for 2nd signal gaus

   // Setting parameters for fitting parameters

   // signal 1 amplitude
   // gaus + pol1
   func_sig_1_amp->FixParameter(0,550);
   func_sig_1_amp->FixParameter(1,1.1418174);
   func_sig_1_amp->FixParameter(2,0.32175315);
   func_sig_1_amp->FixParameter(3,296.52772);
   func_sig_1_amp->FixParameter(4,-27.580470);

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

   // signal 2 amplitude
   // signal 1 amp * 2 gaus
   func_sig_2_amp->FixParameter(0,3.63259e-01);
   func_sig_2_amp->FixParameter(1,1.82697e+00);
   func_sig_2_amp->FixParameter(2,1.82220e-01);
   func_sig_2_amp->FixParameter(3,1.27857e+00);
   func_sig_2_amp->FixParameter(4,2.00491e+00);
   func_sig_2_amp->FixParameter(5,6.26903e-01);

   // Fitting parameter functions
   gr1->Fit("func_sig_1_sigma","RBQ");
   gr2->Fit("func_sig_1_amp","RBQ");
   gr3->Fit("func_sig_2_amp","RBQ");

   // Creating strings for parameter values
   ostringstream signal_1_amp, signal_mean, signal_1_sigma, signal_2_amp, signal_2_sigma, total_function;

   TH2F *Signal_Function = new TH2F("Signal_Function","",500,0.3,0.8,300,0,3);
   TF1 *signal_function = new TF1("signal_function","gaus(0) + gaus(3)",0.35,0.8);
   TF1 *background_function = new TF1("background_function","gaus(0) + gaus(3)",0.35,0.8);

   // Creating signal function and histogram
   // Loop over kaon momentum
   for(Int_t y_pos=1; y_pos < 300; y_pos++){

      // Set signal function parameters based on values obtained from the 1D fit
      // parameter functions for current bin in momentum
      signal_function->SetParameter(0,func_sig_1_amp->Eval(Signal_Function->GetYaxis()->GetBinCenter(y_pos)));
      signal_function->SetParameter(1,func_sig_1_mean->Eval(Signal_Function->GetYaxis()->GetBinCenter(y_pos)));
      signal_function->SetParameter(2,func_sig_1_sigma->Eval(Signal_Function->GetYaxis()->GetBinCenter(y_pos)));
      signal_function->SetParameter(3,func_sig_1_amp->Eval(Signal_Function->GetYaxis()->GetBinCenter(y_pos))*func_sig_2_amp->Eval(Signal_Function->GetYaxis()->GetBinCenter(y_pos)));
      signal_function->SetParameter(4,func_sig_1_mean->Eval(Signal_Function->GetYaxis()->GetBinCenter(y_pos)));
      signal_function->SetParameter(5,2 * func_sig_1_sigma->Eval(Signal_Function->GetYaxis()->GetBinCenter(y_pos)));

      // Loop over kaon mass
      for(Int_t x_pos=1; x_pos < 500; x_pos++){

         // Set bin content for 2D signal histogram based on current momentum
         // and kaon mass bin
         Signal_Function->SetBinContent(x_pos,y_pos,signal_function->Eval(Signal_Function->GetXaxis()->GetBinCenter(x_pos),Signal_Function->GetYaxis()->GetBinCenter(y_pos)));

      }
   }


   TH2F * Kaon_Mass_Momentum = new TH2F();
   Kaon_Mass_Momentum = (TH2F*)hist_S1->Project3D("xz")->Clone();

   Kaon_Mass_Momentum->Add(Signal_Function,-1);
   Kaon_Mass_Momentum->Draw("colz");

}
