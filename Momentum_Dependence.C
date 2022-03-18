// Making lots of projections in momentum

{
   // Get input file
   TFile *f1 = new TFile("/media/mn688/Elements1/PhD/Analysis_Output/Hexaquark/RGA_Fall2018_Outbending_at_least_1e1KpFD_Tree_Total_24022022_Total_Scaling_16032022_01.root");
   // Get histogram from input file
   TH3F* hist_S1=(TH3F*)f1->Get("h_S1_Kaon_Momentum__Miss_Mass__Kaon_Mass");

   // histograms looking at fit parameters as a function of momentum
   TH1F *h_sigma_1 = new TH1F("h_sigma_1","sigma 1 relationship",160,1,161);
   TH1F *h_mean_signal = new TH1F("h_mean_signal","signal mean relationship",160,1,161);
   TH1F *h_amp_1 = new TH1F("h_amp_1","amplitude 1 relationship",160,1,161);
   TH1F *h_sigma_factor = new TH1F("h_sigma_factor","sigma 2 relationship to sigma 1",160,1,161);
   TH1F *h_amp_2 = new TH1F("h_amp_2","amplitude 2 relationship",160,1,161);
   TH1F *h_sigma_3 = new TH1F("h_sigma_3","sigma 3 relationship",160,1,161);
   TH1F *h_chi2 = new TH1F("h_chi2","chi^{2} relationship",160,1,161);

   // 2d histograms for the signal subtraction
   TH2F *h_signal_2d = new TH2F("h_signal_2d","2D signal",500,0.3,0.8,300,0,3);


   // Changing loop to momentum bin values
   Double_t bin_1, bin_2;
   // Reduced chi^2
   Double_t Chi2;
   // Function range
   Double_t Range_min, Range_max;

   // TF1 *func1 = new TF1(funcname.str().c_str(),"gaus(0) + pol3(3)",0.365,0.6);
   TH1F *Missing_Mass_S1_Projections[160];
   TH1F *Kaon_Mass_S1_Projections[160];

   // Create arrays of the various functions
   TF1 *func1[160];
   TF1 *signal_1[160];
   TF1 *signal_2[160];
   TF1 *signal_total[160];
   TF1 *background_1[160];
   TF1 *background_2[160];
   TF1 *background_total[160];

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
      if(bin_1 < 1.35) continue;
      // if(bin_1 > 1.1 && bin_1 < 1.35) continue;

      // Kaons only plotted up to 3 GeV
      if(bin_1 > 3) continue;
      // cout<<"bin "<<i<<" 1st bin "<<bin_1<<" 2nd bin "<<bin_2<<endl;

      Missing_Mass_S1_Projections[i] = (TH1F*) hist_S1->ProjectionY("",hist_S1->GetXaxis()->FindBin(bin_1),hist_S1->GetXaxis()->FindBin(bin_1),0,hist_S1->GetNbinsZ())->Clone();
      Kaon_Mass_S1_Projections[i] = (TH1F*) hist_S1->ProjectionZ("",hist_S1->GetXaxis()->FindBin(bin_1),hist_S1->GetXaxis()->FindBin(bin_1),0,hist_S1->GetNbinsY())->Clone();

      // Rebin projections
      // Strangeness 2 looks good with factor of 2
      // Kaon_Mass_S1_Projections[i]->Rebin(2);

      // Setting range of functions
      ostringstream funcname, sig1, sig2, sigtotal, back1, back2, backtotal;
      funcname <<"func_"<<i;
      sig1 <<"sig1_"<<i;
      sig2 <<"sig2_"<<i;
      back1 <<"back1_"<<i;
      back2 <<"back2_"<<i;
      func1[i] = new TF1(funcname.str().c_str(),"[0] * exp(-pow(x-[1],2) / (2 * pow([2],2))) + [3] * exp(-pow(x-[1],2) / (2 * pow([2]*[4],2))) + gaus(5) + [8]*[9] + [9]*x",0.365,0.65);
      signal_1[i] = new TF1(sig1.str().c_str(),"gaus(0)",0.365,0.65);
      signal_2[i] = new TF1(sig2.str().c_str(),"gaus(0)",0.365,0.65);
      signal_total[i] = new TF1(sigtotal.str().c_str(),"gaus(0) + gaus(3)",0.365,0.65);
      background_1[i] = new TF1(back1.str().c_str(),"gaus(0)",0.0,0.65);
      background_2[i] = new TF1(back2.str().c_str(),"pol1(0)",0.365,0.65);
      background_total[i] = new TF1(backtotal.str().c_str(),"gaus(0) + pol1(3)",0.365,0.65);



      // Setting parameters based on functions determined
      par0 = Kaon_Mass_S1_Projections[i]->GetMaximum() / 2; // amplitude for 1st signal gaus
      // par0 = 0.08024 * pow(i,2) - 20.9076 * i + 1827.356; // amplitude for 1st signal gaus
      par1 = 0.493; // mean for both signal gauss
      // par2 = 0.000000857 * pow(i+1,2) + 0.00001848 * (i+1) + 0.01094; // sigma for 1st signal gaus
      par2 = 0.00000135 * pow(i+1,2) - 0.00005776 * (i+1) + 0.01323; // sigma for 1st signal gaus
      par3 = Kaon_Mass_S1_Projections[i]->GetMaximum() / 3; // amplitude for 2nd signal gaus
      par4 = -0.00001617 * pow(i+1,2) + 0.001323 * (i+1) + 1.928; // ratio between signal gaussian sigmas
      par5 = 2 * Kaon_Mass_S1_Projections[i]->GetBinContent(Kaon_Mass_S1_Projections[i]->FindBin(0.365)); // amplitude for back gaus
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
      func1[i]->SetParLimits(0,Kaon_Mass_S1_Projections[i]->GetMaximum() / 30, Kaon_Mass_S1_Projections[i]->GetMaximum()); // amplitude for 1st signal gaus
      // func1[i]->SetParLimits(0,par0 * 0.7, par0 * 1.3); // amplitude for 1st signal gaus
      func1[i]->SetParLimits(1,par1 * 0.98, par1 * 1.02); // mean for both signal gaus
      func1[i]->SetParLimits(2,0.8 * par2, 1.2 * par2); // sigma for 1st signal gaus
      func1[i]->SetParLimits(3,Kaon_Mass_S1_Projections[i]->GetMaximum() / 30, Kaon_Mass_S1_Projections[i]->GetMaximum()); // amplitude for 2nd signal gaus
      func1[i]->SetParLimits(4,0.85 * par4, 1.15 * par4); // ratio between signal gaussian sigmas
      func1[i]->SetParLimits(5,Kaon_Mass_S1_Projections[i]->GetBinContent(Kaon_Mass_S1_Projections[i]->FindBin(0.365)) / 5, 4 * Kaon_Mass_S1_Projections[i]->GetMaximum()); // amplitude for back gaus
      func1[i]->SetParLimits(6,par6 * 0.98, par6 * 1.02); // mean for back gaus (pion mass)
      func1[i]->SetParLimits(7,0.1, 0.8); // sigma for back gaus
      func1[i]->SetParLimits(8,-1000, -0.7); // constant for back pol1
      func1[i]->SetParLimits(9,-500, 0); // slope for back pol1


      Kaon_Mass_S1_Projections[i]->Fit(funcname.str().c_str(),"RBQ");

      Chi2 = func1[i]->GetChisquare() / func1[i]->GetNDF();


      signal_1[i]->SetParameter(0,func1[i]->GetParameter(0));
      signal_1[i]->SetParameter(1,func1[i]->GetParameter(1));
      signal_1[i]->SetParameter(2,func1[i]->GetParameter(2));
      signal_2[i]->SetParameter(0,func1[i]->GetParameter(3));
      signal_2[i]->SetParameter(1,func1[i]->GetParameter(1));
      signal_2[i]->SetParameter(2,func1[i]->GetParameter(2) * func1[i]->GetParameter(4));
      signal_total[i]->SetParameter(0,func1[i]->GetParameter(0));
      signal_total[i]->SetParameter(1,func1[i]->GetParameter(1));
      signal_total[i]->SetParameter(2,func1[i]->GetParameter(2));
      signal_total[i]->SetParameter(3,func1[i]->GetParameter(3));
      signal_total[i]->SetParameter(4,func1[i]->GetParameter(1));
      signal_total[i]->SetParameter(5,func1[i]->GetParameter(2) * func1[i]->GetParameter(4));
      background_1[i]->SetParameter(0,func1[i]->GetParameter(5));
      background_1[i]->SetParameter(1,func1[i]->GetParameter(6));
      background_1[i]->SetParameter(2,func1[i]->GetParameter(7));
      background_2[i]->SetParameter(0,func1[i]->GetParameter(8) * func1[i]->GetParameter(9));
      background_2[i]->SetParameter(1,func1[i]->GetParameter(9));
      background_total[i]->SetParameter(0,func1[i]->GetParameter(5));
      background_total[i]->SetParameter(1,func1[i]->GetParameter(6));
      background_total[i]->SetParameter(2,func1[i]->GetParameter(7));
      background_total[i]->SetParameter(3,func1[i]->GetParameter(8));
      background_total[i]->SetParameter(4,func1[i]->GetParameter(9));

      h_chi2->Fill(i+1,func1[i]->GetChisquare() / func1[i]->GetNDF());
      h_amp_1->Fill(i+1,func1[i]->GetParameter(0));
      h_amp_1->SetBinError(i+1,func1[i]->GetParError(0));
      h_mean_signal->Fill(i+1,func1[i]->GetParameter(1));
      h_mean_signal->SetBinError(i+1,func1[i]->GetParError(1));
      h_sigma_1->Fill(i+1,func1[i]->GetParameter(2));
      h_sigma_1->SetBinError(i+1,func1[i]->GetParError(2));
      h_amp_2->Fill(i+1,func1[i]->GetParameter(3));
      h_amp_2->SetBinError(i+1,func1[i]->GetParError(3));
      h_sigma_factor->Fill(i+1,func1[i]->GetParameter(4));
      h_sigma_factor->SetBinError(i+1,func1[i]->GetParError(4));
      h_sigma_3->Fill(i+1,func1[i]->GetParameter(7));
      h_sigma_3->SetBinError(i+1,func1[i]->GetParError(7));

   }
   TF1 *func_sig_1_amp = new TF1("func_sig_1_amp","pol2(0)",40,161);
   TF1 *func_sig_mean = new TF1("func_sig_mean","pol4(0)",40,161);
   TF1 *func_sig_1_sigma = new TF1("func_sig_1_sigma","pol2(0)",40,161);
   TF1 *func_sig_2_sigma = new TF1("func_sig_2_sigma","pol3(0)",40,161);

   // Setting parameters for fitting parameters
   // signal 1 amplitude
   func_sig_1_amp->FixParameter(0,1002);
   func_sig_1_amp->FixParameter(1,-11.92);
   func_sig_1_amp->FixParameter(2,0.04633);

   // signal mean
   func_sig_mean->FixParameter(0,0.483238);
   func_sig_mean->FixParameter(1,0.000674655);
   func_sig_mean->FixParameter(2,-1.45551e-05);
   func_sig_mean->FixParameter(3,1.22688e-07);
   func_sig_mean->FixParameter(4,-3.49508e-10);

   // signal 1 sigma
   func_sig_1_sigma->FixParameter(0,0.01385);
   func_sig_1_sigma->FixParameter(1,-6.264e-05);
   func_sig_1_sigma->FixParameter(2,1.336e-06);

   // signal 2 sigma
   func_sig_1_sigma->FixParameter(0,1.62832);
   func_sig_1_sigma->FixParameter(1,1.54654e-02);
   func_sig_1_sigma->FixParameter(2,- 1.99258e-04);
   func_sig_1_sigma->FixParameter(3,7.19587e-07);

   // Fitting parameter functions
   h_amp_1->Fit("func_sig_1_amp","RB");
   h_mean_signal->Fit("func_sig_mean","RB");
   h_sigma_1->Fit("func_sig_1_sigma","RB");
   h_sigma_factor->Fit("func_sig_2_sigma","RB");

   // [0] * exp(-pow(x-[1],2) / (2 * pow([2],2))) + [3] * exp(-pow(x-[1],2) / (2 * pow([2]*[4],2)))

   ostringstream signal_1_amp, signal_mean, signal_1_sigma, signal_2_amp, signal_2_sigma, total_function;

   signal_1_amp <<"0.04633 * pow(y,2) - 11.92 * y + 1002";
   signal_mean <<"(-3.49508e-10 * pow(y,4) + 1.22688e-07 * pow(y,3) - 1.45551e-05 * pow(y,2) + 0.000674655 * y + 0.483238)";
   signal_1_sigma <<" 1.336e-06 * pow(y,2) - 6.264e-05 * y + 0.01385";
   signal_2_amp <<"5.52306e-04 * pow(y,3) - 1.80445e-01 * pow(y,2) + 1.49715e+01 * y + 1.30435e+02";
   signal_2_sigma <<" 7.19587e-07 * pow(y,3) - 1.99258e-04 * pow(y,2) + 1.54654e-02 * y + 1.62832";

   total_function << signal_1_amp.str().c_str() << " * exp(-pow(x-"
                  << signal_mean.str().c_str() << ",2) / (2 * pow("
                  << signal_1_sigma.str().c_str() << ",2)))"
                  << " + " << signal_2_amp.str().c_str()
                  << " * exp(-pow(x - " << signal_mean.str().c_str()
                  << " ,2) / (2 * pow(" << signal_1_sigma.str().c_str()
                  << " * " << signal_2_sigma.str().c_str()
                  << " ,2)))";

   TF2 *signal_functions = new TF2("signal_functions",total_function.str().c_str(),0.3, 0.8, 0, 3);
   cout << total_function.str().c_str() << endl;
   // h_signal_2d->Draw();
   // signal_functions->Draw("colz,same");
}
