
{

   //////////////////////////////////////////////////////////////////////////////
   //// Setting up files and histograms    ////////////////////////////////////
   //////////////////////////////////////////////////////////////////////////////

   // Input file
   // RGA data
   TFile *f1=new TFile("/media/mn688/Elements1/PhD/Analysis_Output/Hexaquark/RGB_Sping2019_Inbending_at_least_1e1KpFD_Tree_Total_01032022_Total_Scaling_07032022_01.root");

   // Getting the multidimensional histogram plots
   TH3F *h_S1_Photon_Energy__Miss_Mass__Kaon_Mass = (TH3F*)f1->Get("h_S1_Photon_Energy__Miss_Mass__Kaon_Mass");
   TH3F *h_S2_Photon_Energy__Miss_Mass__Kaon_Mass = (TH3F*)f1->Get("h_S2_Photon_Energy__Miss_Mass__Kaon_Mass");
   TH3F *h_S3_Photon_Energy__Miss_Mass__Kaon_Mass = (TH3F*)f1->Get("h_S3_Photon_Energy__Miss_Mass__Kaon_Mass");

   // Make array for all the 1D projections
   TH1D *h_projectionx_S1_sig[500];
   TH1D *h_projectionx_S1_back[500];
   TH1D *h_projectionx_S2_sig[500];
   TH1D *h_projectionx_S2_back[500];
   TH1D *h_projectionx_S3_sig[500];
   TH1D *h_projectionx_S3_back[500];

   //////////////////////////////////////////////////////////////////////////////
   //// Take all the projections    ///////////////////////////////////////
   //////////////////////////////////////////////////////////////////////////////

   // Strangeness 1
   // Photon energy
   TH1D *S1_photon_energy_total = (TH1D*)h_S1_Photon_Energy__Miss_Mass__Kaon_Mass->Project3D("x")->Clone("S1_photon_energy_total");
   // missing mass of electron and 1 kaon
   TH1D *S1_miss_mass_total = (TH1D*)h_S1_Photon_Energy__Miss_Mass__Kaon_Mass->Project3D("y")->Clone("S1_miss_mass_total");
   // Mass of kaon 1
   TH1D *S1_kp1_mass_total = (TH1D*)h_S1_Photon_Energy__Miss_Mass__Kaon_Mass->Project3D("z")->Clone("S1_kp1_mass_total");


   h_S2_Photon_Energy__Miss_Mass__Kaon_Mass->Rebin3D(1,6,2);

   // Strangeness 2
   // Photon energy
   TH1D *S2_photon_energy_total = (TH1D*)h_S2_Photon_Energy__Miss_Mass__Kaon_Mass->Project3D("x")->Clone("S2_photon_energy_total");
   // missing mass of electron and 2 kaons
   TH1D *S2_miss_mass_total = (TH1D*)h_S2_Photon_Energy__Miss_Mass__Kaon_Mass->Project3D("y")->Clone("S2_miss_mass_total");
   // Mass of kaon 1
   TH1D *S2_kp1_mass_total = (TH1D*)h_S2_Photon_Energy__Miss_Mass__Kaon_Mass->Project3D("z")->Clone("S2_kp1_mass_total");


   h_S3_Photon_Energy__Miss_Mass__Kaon_Mass->Rebin3D(1,100,5);

   // Strangeness 3
   // Photon energy
   TH1D *S3_photon_energy_total = (TH1D*)h_S3_Photon_Energy__Miss_Mass__Kaon_Mass->Project3D("x")->Clone("S3_photon_energy_total");
   // missing mass of electron and 3 kaons
   TH1D *S3_miss_mass_total = (TH1D*)h_S3_Photon_Energy__Miss_Mass__Kaon_Mass->Project3D("y")->Clone("S3_miss_mass_total");
   // Mass of kaon 1
   TH1D *S3_kp1_mass_total = (TH1D*)h_S3_Photon_Energy__Miss_Mass__Kaon_Mass->Project3D("z")->Clone("S3_kp1_mass_total");

   //////////////////////////////////////////////////////////////////////////////
   //// Creating and fitting functions    ///////////////////////////////////////
   //////////////////////////////////////////////////////////////////////////////

   // Define functions for fitting kaon calculated mass
   // Functions for strangeness 1 - kaon 1
   // Total function
   TF1 *func1_S1 = new TF1("func1_S1","gaus(0) + gaus(3) + gaus(6)",0.365,0.72);
   // Signal narrow gaussian
   TF1 *func2_S1 = new TF1("func2_S1","gaus(0)",0.365,0.72);
   // Signal wide gaussian
   TF1 *func3_S1 = new TF1("func3_S1","gaus(0)",0.365,0.72);
   // Total signal
   TF1 *func4_S1 = new TF1("func4_S1","gaus(0) + gaus(3)",0.365,0.72);
   // Background
   TF1 *func5_S1 = new TF1("func5_S1","gaus(0)",0.365,0.72);

   // Functions for strangeness 2 - kaon 1
   // Total function
   TF1 *func1_S2 = new TF1("func1_S2","gaus(0) + gaus(3) + gaus(6)",0.365,0.72);
   // Signal narrow gaussian
   TF1 *func2_S2 = new TF1("func2_S2","gaus(0)",0.365,0.72);
   // Signal wide gaussian
   TF1 *func3_S2 = new TF1("func3_S2","gaus(0)",0.365,0.72);
   // Total signal
   TF1 *func4_S2 = new TF1("func4_S2","gaus(0) + gaus(3)",0.365,0.72);
   // Background
   TF1 *func5_S2 = new TF1("func5_S2","gaus(0)",0.365,0.72);

   // Functions for strangeness 3 - kaon 1
   // Total function
   TF1 *func1_S3 = new TF1("func1_S3","gaus(0) + gaus(3) + pol3(6)",0.365,0.72);
   // Signal narrow gaussian
   TF1 *func2_S3 = new TF1("func2_S3","gaus(0)",0.365,0.72);
   // Signal wide gaussian
   TF1 *func3_S3 = new TF1("func3_S3","gaus(0)",0.365,0.72);
   // Total signal
   TF1 *func4_S3 = new TF1("func4_S3","gaus(0) + gaus(3)",0.365,0.72);
   // Background
   TF1 *func5_S3 = new TF1("func5_S3","gaus(0)",0.365,0.72);


   // Setting parameters before fitting
   // Strangeness 1 - kaon 1
   func1_S1->SetParameter(0,S1_kp1_mass_total->GetMaximum() / 2);
   func1_S1->SetParameter(1,0.493);
   func1_S1->SetParameter(2,0.01);
   func1_S1->SetParameter(3,S1_kp1_mass_total->GetMaximum() / 2);
   func1_S1->SetParameter(4,0.493);
   func1_S1->SetParameter(5,0.03);
   func1_S1->SetParameter(6,S1_kp1_mass_total->GetMaximum() / 4);
   func1_S1->SetParameter(7,0.345);
   func1_S1->SetParameter(8,0.5);

   // func5_S1->SetParameter(0,S1_kp1_mass_total->GetMaximum() / 4);
   // func5_S1->SetParLimits(0,S1_kp1_mass_total->GetMaximum() / 9, S1_kp1_mass_total->GetMaximum() / 2.5);
   // func5_S1->SetParameter(1,0.345);
   // func5_S1->SetParameter(2,0.5);

   // Setting parameter limits before fitting
   func1_S1->SetParLimits(0,S1_kp1_mass_total->GetMaximum() / 5,S1_kp1_mass_total->GetMaximum()); // amplitude for first gauss
   func1_S1->SetParLimits(1,0.480,0.505); // mean for first gauss
   func1_S1->SetParLimits(2,0.005,0.02); // sigma for first gauss
   func1_S1->SetParLimits(3,S1_kp1_mass_total->GetMaximum() / 5,S1_kp1_mass_total->GetMaximum()); // amplitude for second gauss
   func1_S1->SetParLimits(4,0.480,0.505); // mean for second gauss
   func1_S1->SetParLimits(5,0.005,0.04); // sigma for second gauss
   func1_S1->SetParLimits(7,0.1,0.38); // mean for background gauss

   // Strangeness 2 - kaon 1
   func1_S2->SetParameter(0,S2_kp1_mass_total->GetMaximum() / 2);
   func1_S2->SetParameter(1,0.493);
   func1_S2->SetParameter(2,0.01);
   func1_S2->SetParameter(3,S2_kp1_mass_total->GetMaximum() / 2);
   func1_S2->SetParameter(4,0.493);
   func1_S2->SetParameter(5,0.03);
   func1_S2->SetParameter(6,S1_kp1_mass_total->GetMaximum() / 4);
   func1_S2->SetParameter(7,0.345);
   func1_S2->SetParameter(8,0.5);


   // Setting parameter limits before fitting
   func1_S2->SetParLimits(0,S2_kp1_mass_total->GetMaximum() / 5,S2_kp1_mass_total->GetMaximum()); // amplitude for first gauss
   func1_S2->SetParLimits(1,0.480,0.505); // mean for first gauss
   func1_S2->SetParLimits(2,0.005,0.02); // sigma for first gauss
   func1_S2->SetParLimits(3,S2_kp1_mass_total->GetMaximum() / 5,S2_kp1_mass_total->GetMaximum()); // amplitude for second gauss
   func1_S2->SetParLimits(4,0.480,0.505); // mean for second gauss
   func1_S2->SetParLimits(5,0.005,0.04); // sigma for second gauss
   func1_S2->SetParLimits(7,0.1,0.38); // mean for background gauss

   // Strangeness 3 - kaon 1
   func1_S3->SetParameter(0,S3_kp1_mass_total->GetMaximum() / 2);
   func1_S3->SetParameter(1,0.493);
   func1_S3->SetParameter(2,0.01);
   func1_S3->SetParameter(3,S3_kp1_mass_total->GetMaximum() / 2);
   func1_S3->SetParameter(4,0.493);
   func1_S3->SetParameter(5,0.03);
   func1_S3->SetParameter(6,S1_kp1_mass_total->GetMaximum() / 4);
   func1_S3->SetParameter(7,0.345);
   func1_S3->SetParameter(8,0.5);


   // Setting parameter limits before fitting
   func1_S3->SetParLimits(0,S3_kp1_mass_total->GetMaximum() / 5,S3_kp1_mass_total->GetMaximum()); // amplitude for first gauss
   func1_S3->SetParLimits(1,0.480,0.505); // mean for first gauss
   func1_S3->SetParLimits(2,0.005,0.02); // sigma for first gauss
   func1_S3->SetParLimits(3,S3_kp1_mass_total->GetMaximum() / 5,S3_kp1_mass_total->GetMaximum()); // amplitude for second gauss
   func1_S3->SetParLimits(4,0.480,0.505); // mean for second gauss
   func1_S3->SetParLimits(5,0.005,0.04); // sigma for second gauss
   func1_S3->SetParLimits(7,0.1,0.38); // mean for background gauss



   // Fitting functions and getting parameters
   // Strangeness 1 - kaon 1
   S1_kp1_mass_total->Fit("func1_S1","RB");
   func2_S1->FixParameter(0, func1_S1->GetParameter(0));
   func2_S1->FixParameter(1, func1_S1->GetParameter(1));
   func2_S1->FixParameter(2, func1_S1->GetParameter(2));
   func3_S1->FixParameter(0, func1_S1->GetParameter(3));
   func3_S1->FixParameter(1, func1_S1->GetParameter(4));
   func3_S1->FixParameter(2, func1_S1->GetParameter(5));
   func4_S1->FixParameter(0, func1_S1->GetParameter(0));
   func4_S1->FixParameter(1, func1_S1->GetParameter(1));
   func4_S1->FixParameter(2, func1_S1->GetParameter(2));
   func4_S1->FixParameter(3, func1_S1->GetParameter(3));
   func4_S1->FixParameter(4, func1_S1->GetParameter(4));
   func4_S1->FixParameter(5, func1_S1->GetParameter(5));
   func5_S1->FixParameter(0, func1_S1->GetParameter(6));
   func5_S1->FixParameter(1, func1_S1->GetParameter(7));
   func5_S1->FixParameter(2, func1_S1->GetParameter(8));
   // func5_S1->FixParameter(3, func1_S1->GetParameter(9));
   // func5_S1->FixParameter(4, func1_S1->GetParameter(10));
   // func5_S1->FixParameter(5, func1_S1->GetParameter(11));


   // Strangeness 2 - kaon 1
   S2_kp1_mass_total->Fit("func1_S2","RB");
   func2_S2->FixParameter(0, func1_S2->GetParameter(0));
   func2_S2->FixParameter(1, func1_S2->GetParameter(1));
   func2_S2->FixParameter(2, func1_S2->GetParameter(2));
   func3_S2->FixParameter(0, func1_S2->GetParameter(3));
   func3_S2->FixParameter(1, func1_S2->GetParameter(4));
   func3_S2->FixParameter(2, func1_S2->GetParameter(5));
   func4_S2->FixParameter(0, func1_S2->GetParameter(0));
   func4_S2->FixParameter(1, func1_S2->GetParameter(1));
   func4_S2->FixParameter(2, func1_S2->GetParameter(2));
   func4_S2->FixParameter(3, func1_S2->GetParameter(3));
   func4_S2->FixParameter(4, func1_S2->GetParameter(4));
   func4_S2->FixParameter(5, func1_S2->GetParameter(5));
   func5_S2->FixParameter(0, func1_S2->GetParameter(6));
   func5_S2->FixParameter(1, func1_S2->GetParameter(7));
   func5_S2->FixParameter(2, func1_S2->GetParameter(8));
   // func5_S2->FixParameter(3, func1_S2->GetParameter(9));
   // func5_S2->FixParameter(4, func1_S2->GetParameter(10));
   // func5_S2->FixParameter(5, func1_S2->GetParameter(11));


   // Strangeness 3 - kaon 1
   S3_kp1_mass_total->Fit("func1_S3","RB");
   func2_S3->FixParameter(0, func1_S3->GetParameter(0));
   func2_S3->FixParameter(1, func1_S3->GetParameter(1));
   func2_S3->FixParameter(2, func1_S3->GetParameter(2));
   func3_S3->FixParameter(0, func1_S3->GetParameter(3));
   func3_S3->FixParameter(1, func1_S3->GetParameter(4));
   func3_S3->FixParameter(2, func1_S3->GetParameter(5));
   func4_S3->FixParameter(0, func1_S3->GetParameter(0));
   func4_S3->FixParameter(1, func1_S3->GetParameter(1));
   func4_S3->FixParameter(2, func1_S3->GetParameter(2));
   func4_S3->FixParameter(3, func1_S3->GetParameter(3));
   func4_S3->FixParameter(4, func1_S3->GetParameter(4));
   func4_S3->FixParameter(5, func1_S3->GetParameter(5));
   func5_S3->FixParameter(0, func1_S3->GetParameter(6));
   func5_S3->FixParameter(1, func1_S3->GetParameter(7));
   func5_S3->FixParameter(2, func1_S3->GetParameter(8));
   func5_S3->FixParameter(3, func1_S3->GetParameter(9));
   // func5_S3->FixParameter(4, func1_S3->GetParameter(10));
   // func5_S3->FixParameter(5, func1_S3->GetParameter(11));

   //////////////////////////////////////////////////////////////////////////////
   //// Background subtraction         //////////////////////////////////////////
   //////////////////////////////////////////////////////////////////////////////

   // Create scaling factors for background and signal for each strangeness
   Double_t S1_Scaling_Factor_sig, S1_Scaling_Factor_back;
   Double_t S2_Scaling_Factor_sig, S2_Scaling_Factor_back;
   Double_t S3_Scaling_Factor_sig, S3_Scaling_Factor_back;

   //////////////////////////////////////////////////////////////////////////////
   // Strangeness 1 - kaon 1

   // Loop over the kaon mass bins (z-axis) and make projection of MM (y-axis)
   for(Int_t bin = 1; bin < h_S1_Photon_Energy__Miss_Mass__Kaon_Mass->GetNbinsZ(); bin++){

      // Create names for the signal histograms
      ostringstream h_projectionx_S1_sig_name;
      h_projectionx_S1_sig_name << "h_projectionx_S1_sig_" << bin;

      // Taking projection of MM from kaon mass bin for sig and back
      h_projectionx_S1_sig[bin] = (TH1D*)h_S1_Photon_Energy__Miss_Mass__Kaon_Mass->ProjectionY("",0,h_S1_Photon_Energy__Miss_Mass__Kaon_Mass->GetNbinsX(),bin,bin)->Clone(h_projectionx_S1_sig_name.str().c_str());
      h_projectionx_S1_back[bin] = (TH1D*)h_S1_Photon_Energy__Miss_Mass__Kaon_Mass->ProjectionY("",0,h_S1_Photon_Energy__Miss_Mass__Kaon_Mass->GetNbinsX(),bin,bin)->Clone("h_projectionx_S1_back");

      // Make sure there are reasonable statistics in the projection
      if(h_projectionx_S1_sig[bin]->Integral() > 100){

         // Determine the scaling factor looking at the signal function
         S1_Scaling_Factor_sig = func5_S1->Eval(S1_kp1_mass_total->GetBinCenter(bin)) / h_projectionx_S1_sig[bin]->Integral() ;

         // Determine the scaling factor looking at the background function
         S1_Scaling_Factor_back = (func5_S1->Eval(S1_kp1_mass_total->GetBinCenter(bin)) -
         func4_S1->Eval(S1_kp1_mass_total->GetBinCenter(bin))) / h_projectionx_S1_back[bin]->Integral();

         // If background scaling factor is less than 0 set to zero
         if(S1_Scaling_Factor_back < 0) S1_Scaling_Factor_back = 0;

         // Scale the signal and background histograms accordingly
         h_projectionx_S1_sig[bin]->Scale(S1_Scaling_Factor_sig);
         h_projectionx_S1_back[bin]->Scale(S1_Scaling_Factor_back);

      }
   }

   // Adding all the histograms together to get total continous sideband subtracted result
   for(Int_t bin = 2; bin < h_S1_Photon_Energy__Miss_Mass__Kaon_Mass->GetNbinsZ(); bin++){

      h_projectionx_S1_sig[1]->Add(h_projectionx_S1_sig[bin]);
      h_projectionx_S1_back[1]->Add(h_projectionx_S1_back[bin]);
   }

   // Determining the integral of signal and background where MM < 1, this is
   // done to ensure below threshold background is removed properly
   Double_t back_integral_S1 = h_projectionx_S1_back[1]->Integral(h_projectionx_S1_back[1]->FindBin(0), h_projectionx_S1_back[1]->FindBin(1));
   Double_t sig_integral_S1 = h_projectionx_S1_sig[1]->Integral(h_projectionx_S1_sig[1]->FindBin(0), h_projectionx_S1_sig[1]->FindBin(1));
   // Scale the background accordingly for low MM
   h_projectionx_S1_back[1]->Scale(sig_integral_S1 / back_integral_S1);

   // Subtract the background histogram from the signal
   h_projectionx_S1_sig[1]->Add(h_projectionx_S1_back[1], -1);


   //////////////////////////////////////////////////////////////////////////////
   // Strangeness 2 - kaon 1

   // Loop over the kaon mass bins (z-axis) and make projection of MM (y-axis)
   for(Int_t bin = 1; bin < S2_kp1_mass_total->GetNbinsX(); bin++){

      // Create names for the signal histograms
      ostringstream h_projectionx_S2_sig_name;
      h_projectionx_S2_sig_name << "h_projectionx_S2_sig_" << bin;

      // Taking projection of MM from kaon mass bin for sig and back
      h_projectionx_S2_sig[bin] = (TH1D*)h_S2_Photon_Energy__Miss_Mass__Kaon_Mass->ProjectionY("",0,h_S2_Photon_Energy__Miss_Mass__Kaon_Mass->GetNbinsX(),bin,bin)->Clone(h_projectionx_S2_sig_name.str().c_str());
      h_projectionx_S2_back[bin] = (TH1D*)h_S2_Photon_Energy__Miss_Mass__Kaon_Mass->ProjectionY("",0,h_S2_Photon_Energy__Miss_Mass__Kaon_Mass->GetNbinsX(),bin,bin)->Clone("h_projectionx_S2_back");


      // Make sure there are reasonable statistics in the projection
      if(h_projectionx_S2_sig[bin]->Integral() > 100){

         // Determine the scaling factor looking at the signal function
         S2_Scaling_Factor_sig = func5_S2->Eval(S2_kp1_mass_total->GetBinCenter(bin)) / h_projectionx_S2_sig[bin]->Integral() ;

         // Determine the scaling factor looking at the background function
         S2_Scaling_Factor_back = (func5_S2->Eval(S2_kp1_mass_total->GetBinCenter(bin)) -
         func4_S2->Eval(S2_kp1_mass_total->GetBinCenter(bin))) / h_projectionx_S2_back[bin]->Integral();

         // If background scaling factor is less than 0 set to zero
         if(S2_Scaling_Factor_back < 0) S2_Scaling_Factor_back = 0;

         // Scale the signal and background histograms accordingly
         h_projectionx_S2_sig[bin]->Scale(S2_Scaling_Factor_sig);
         h_projectionx_S2_back[bin]->Scale(S2_Scaling_Factor_back);

      }
   }

   // Adding all the histograms together to get total continous sideband subtracted result
   for(Int_t bin = 2; bin < S2_kp1_mass_total->GetNbinsX(); bin++){

      h_projectionx_S2_sig[1]->Add(h_projectionx_S2_sig[bin]);
      h_projectionx_S2_back[1]->Add(h_projectionx_S2_back[bin]);
   }

   // Determining the integral of signal and background where MM < 1, this is
   // done to ensure below threshold background is removed properly
   Double_t back_integral_S2 = h_projectionx_S2_back[1]->Integral(h_projectionx_S2_back[1]->FindBin(0), h_projectionx_S2_back[1]->FindBin(1));
   Double_t sig_integral_S2 = h_projectionx_S2_sig[1]->Integral(h_projectionx_S2_sig[1]->FindBin(0), h_projectionx_S2_sig[1]->FindBin(1));
   // Scale the background accordingly for low MM
   h_projectionx_S2_back[1]->Scale(sig_integral_S2 / back_integral_S2);

   // Subtract the background histogram from the signal
   h_projectionx_S2_sig[1]->Add(h_projectionx_S2_back[1], -1);


   //////////////////////////////////////////////////////////////////////////////
   // Strangeness 3 - kaon 1

   // Loop over the kaon mass bins (z-axis) and make projection of MM (y-axis)
   for(Int_t bin = 1; bin < S3_kp1_mass_total->GetNbinsX(); bin++){

      // Create names for the signal histograms
      ostringstream h_projectionx_S3_sig_name;
      h_projectionx_S3_sig_name << "h_projectionx_S3_sig_" << bin;

      // Taking projection of MM from kaon mass bin for sig and back
      h_projectionx_S3_sig[bin] = (TH1D*)h_S3_Photon_Energy__Miss_Mass__Kaon_Mass->ProjectionY("",0,h_S3_Photon_Energy__Miss_Mass__Kaon_Mass->GetNbinsX(),bin,bin)->Clone(h_projectionx_S3_sig_name.str().c_str());
      h_projectionx_S3_back[bin] = (TH1D*)h_S3_Photon_Energy__Miss_Mass__Kaon_Mass->ProjectionY("",0,h_S3_Photon_Energy__Miss_Mass__Kaon_Mass->GetNbinsX(),bin,bin)->Clone("h_projectionx_S3_back");


      // Make sure there are reasonable statistics in the projection
      if(h_projectionx_S3_sig[bin]->Integral() > 20){

         // Determine the scaling factor looking at the signal function
         S3_Scaling_Factor_sig = func5_S3->Eval(S3_kp1_mass_total->GetBinCenter(bin)) / h_projectionx_S3_sig[bin]->Integral() ;

         // Determine the scaling factor looking at the background function
         S3_Scaling_Factor_back = (func5_S3->Eval(S3_kp1_mass_total->GetBinCenter(bin)) -
         func4_S3->Eval(S3_kp1_mass_total->GetBinCenter(bin))) / h_projectionx_S3_back[bin]->Integral();

         // If background scaling factor is less than 0 set to zero
         if(S3_Scaling_Factor_back < 0) S3_Scaling_Factor_back = 0;

         // Scale the signal and background histograms accordingly
         h_projectionx_S3_sig[bin]->Scale(S3_Scaling_Factor_sig);
         h_projectionx_S3_back[bin]->Scale(S3_Scaling_Factor_back);

      }
   }

   // Adding all the histograms together to get total continous sideband subtracted result
   for(Int_t bin = 2; bin < S3_kp1_mass_total->GetNbinsX(); bin++){

      h_projectionx_S3_sig[1]->Add(h_projectionx_S3_sig[bin]);
      h_projectionx_S3_back[1]->Add(h_projectionx_S3_back[bin]);
   }

   // Determining the integral of signal and background where MM < 1, this is
   // done to ensure below threshold background is removed properly
   Double_t back_integral_S3 = h_projectionx_S3_back[1]->Integral(h_projectionx_S3_back[1]->FindBin(0), h_projectionx_S3_back[1]->FindBin(1));
   Double_t sig_integral_S3 = h_projectionx_S3_sig[1]->Integral(h_projectionx_S3_sig[1]->FindBin(0), h_projectionx_S3_sig[1]->FindBin(1));
   // Scale the background accordingly for low MM
   h_projectionx_S3_back[1]->Scale(sig_integral_S3 / back_integral_S3);

   // Subtract the background histogram from the signal
   h_projectionx_S3_sig[1]->Add(h_projectionx_S3_back[1], -1);


   // Strangness 1
   func2_S1->SetLineColor(kBlue);
   func3_S1->SetLineColor(kBlue);
   func4_S1->SetLineColor(kBlack);
   func5_S1->SetLineColor(kGreen);

   // Strangness 2
   func2_S2->SetLineColor(kBlue);
   func3_S2->SetLineColor(kBlue);
   func4_S2->SetLineColor(kBlack);
   func5_S2->SetLineColor(kGreen);

   // Strangness 3
   func2_S3->SetLineColor(kBlue);
   func3_S3->SetLineColor(kBlue);
   func4_S3->SetLineColor(kBlack);
   func5_S3->SetLineColor(kGreen);


   // Creating lines for the canvases
   // Zero lines for after sideband subtraction
   auto *l1 = new TLine(0,0,3,0);
   l1->SetLineColor(kRed);
   l1->SetLineWidth(2);

   // Cascade ground state
   auto *cascade_line = new TLine(1.32171,0,1.32171,350);
   // Cascade 1530 state
   auto *cascade_1530_line = new TLine(1.535,0,1.535,350);

   // Change colours of functions to made it clearer
   // Particle lines
   cascade_line->SetLineColor(kRed);
   cascade_1530_line->SetLineColor(kRed);


   auto *c1 = new TCanvas("c1","Strangeness 1 kaon 1 mass",800,800);
   c1->cd();
   S1_kp1_mass_total->Draw();
   func1_S1->Draw("same");
   func2_S1->Draw("same");
   func3_S1->Draw("same");
   func4_S1->Draw("same");
   func5_S1->Draw("same");

   // auto *c2 = new TCanvas("c2","Strageness 1 Before background Subtraction",800,800);
   // c2->cd();
   // S1_miss_mass_total->Draw();
   //
   //
   // auto *c3 = new TCanvas("c3","Strageness 1 After background Subtraction",800,800);
   // c3->cd();
   // h_projectionx_S1_sig[1]->Draw();
   // l1->Draw("same");

   auto *c4 = new TCanvas("c4","Strangeness 2 kaon 1 mass",800,800);
   c4->cd();
   S2_kp1_mass_total->Draw();
   func1_S2->Draw("same");
   func2_S2->Draw("same");
   func3_S2->Draw("same");
   func4_S2->Draw("same");
   func5_S2->Draw("same");

   // auto *c5 = new TCanvas("c5","Strageness 2 Before background Subtraction",800,800);
   // c5->cd();
   // S2_miss_mass_total->Draw();
   //
   //
   // auto *c6 = new TCanvas("c6","Strageness 2 After background Subtraction",800,800);
   // c6->cd();
   // h_projectionx_S2_sig[1]->Draw();
   // l1->Draw("same");
   // cascade_line->Draw("same");
   // cascade_1530_line->Draw("same");

   // auto *c7 = new TCanvas("c7","Strangeness 3 kaon 1 mass",800,800);
   // c7->cd();
   // S3_kp1_mass_total->Draw();
   // func1_S3->Draw("same");
   // func2_S3->Draw("same");
   // func3_S3->Draw("same");
   // func4_S3->Draw("same");
   // func5_S3->Draw("same");

   // auto *c8 = new TCanvas("c8","Strageness 3 Before background Subtraction",800,800);
   // c8->cd();
   // S3_miss_mass_total->Draw();
   //
   //
   // auto *c9 = new TCanvas("c9","Strageness 3 After background Subtraction",800,800);
   // c9->cd();
   // h_projectionx_S3_sig[1]->Draw();
   // l1->Draw("same");
}
