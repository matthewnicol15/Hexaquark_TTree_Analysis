
{

   //////////////////////////////////////////////////////////////////////////////
   ////Define variables for naming    ///////////////////////////////////////////
   //////////////////////////////////////////////////////////////////////////////

   // Information for canvas and histogram name
   ostringstream File_Path;
   ostringstream Data;
   ostringstream Quantity;
   ostringstream Date;
   ostringstream Version;
   ostringstream Output_File_Name;


   //////////////////////////////////////////////////////////////////////////////
   //// Setting up files and histograms    ////////////////////////////////////
   //////////////////////////////////////////////////////////////////////////////

   // Input file
   // RGA data
   TFile *f1=new TFile("/mnt/c/Users/Nics/Documents/RGA_Fall2018_Inbending_at_least_1e1KpFD_Tree_Total_24022022_1M_Scaling_01032022_02.root");


   // Getting the multidimensional histogram plots
   TH3F *h_S1_Photon_Energy__Miss_Mass__Kaon_Mass = (TH3F*)f1->Get("h_S1_Photon_Energy__Miss_Mass__Kaon_Mass");
   TH3F *h_S2_Photon_Energy__Miss_Mass__Kaon_Mass = (TH3F*)f1->Get("h_S2_Photon_Energy__Miss_Mass__Kaon_Mass");
   TH3F *h_S3_Photon_Energy__Miss_Mass__Kaon_Mass = (TH3F*)f1->Get("h_S3_Photon_Energy__Miss_Mass__Kaon_Mass");
   // Get the number of bins
   Int_t zbinmax = h_S2_Photon_Energy__Miss_Mass__Kaon_Mass->GetNbinsZ();


   //////////////////////////////////////////////////////////////////////////////
   //// Take all the projections    ///////////////////////////////////////
   //////////////////////////////////////////////////////////////////////////////

   // Strangeness 1
   // Photon energy
   TH1D *S1_kp1_mass_total = (TH1D*)h_S1_Photon_Energy__Miss_Mass__Kaon_Mass->Project3D("z")->Clone("S1_kp1_mass_total");

   // Strangeness 2
   // Mass of kaon 1
   TH1D *S2_kp1_mass_total = (TH1D*)h_S2_Photon_Energy__Miss_Mass__Kaon_Mass->Project3D("z")->Clone("S2_kp1_mass_total");

   // Strangeness 3
   // Mass of kaon 1
   TH1D *S3_kp1_mass_total = (TH1D*)h_S3_Photon_Energy__Miss_Mass__Kaon_Mass->Project3D("z")->Clone("S3_kp1_mass_total");

   //////////////////////////////////////////////////////////////////////////////
   //// Creating and fitting functions    ///////////////////////////////////////
   //////////////////////////////////////////////////////////////////////////////


   // Define functions for fitting kaon calculated mass
   // Functions for strangeness 1 - kaon 1
   TF1 *func1 = new TF1("func1","gaus(0) + pol4(3)",0.36,0.7);
   TF1 *func2 = new TF1("func2","gaus(0)",0.36,0.7);
   // TF1 *func3 = new TF1("func3","gaus(0)",0.36,0.7);
   // TF1 *func4 = new TF1("func4","gaus(0) + gaus(3)",0.36,0.7);
   TF1 *func5 = new TF1("func5","pol4(0)",0.36,0.7);

   // // Functions for strangeness 2 - kaon 1
   // TF1 *func1_s2_kp1 = new TF1("func1_s2_kp1","gaus(0) + pol3(3) + gaus(7)",0.36,0.7);
   // TF1 *func2_s2_kp1 = new TF1("func2_s2_kp1","gaus(0)",0.36,0.7);
   // TF1 *func3_s2_kp1 = new TF1("func3_s2_kp1","pol3(0)",0.36,0.7);
   // TF1 *func4_s2_kp1 = new TF1("func4_s2_kp1","gaus(0)",0.36,0.7);
   // TF1 *func5_s2_kp1 = new TF1("func5_s2_kp1","gaus(0) + gaus(3)",0.36,0.7);
   //
   // // Functions for strangeness 3 - kaon 1
   // TF1 *func1_s3_kp1 = new TF1("func1_s3_kp1","gaus(0) + pol3(3) + gaus(7)",0.36,0.7);
   // TF1 *func2_s3_kp1 = new TF1("func2_s3_kp1","gaus(0)",0.36,0.7);
   // TF1 *func3_s3_kp1 = new TF1("func3_s3_kp1","pol3(0)",0.36,0.7);
   // TF1 *func4_s3_kp1 = new TF1("func4_s3_kp1","gaus(0)",0.36,0.7);
   // TF1 *func5_s3_kp1 = new TF1("func5_s3_kp1","gaus(0) + gaus(3)",0.36,0.7);

   // Setting parameters before fitting
   // Strangeness 1 - kaon 1
   func1->SetParameter(0,S1_kp1_mass_total->GetMaximum() / 2);
   func1->SetParameter(1,0.493);
   func1->SetParameter(2,0.02);
   // func1->SetParameter(3,S1_kp1_mass_total->GetMaximum() / 3);
   // func1->SetParameter(4,0.493);
   // func1->SetParameter(5,0.03);

   // Setting parameter limits before fitting
   func1->SetParLimits(0,S1_kp1_mass_total->GetMaximum() / 3,S1_kp1_mass_total->GetMaximum()); // amplitude for first gauss
   func1->SetParLimits(1,0.485,0.505); // mean for first gauss
   func1->SetParLimits(2,0.005,0.06); // sigma for first gauss
   // func1->SetParLimits(3,S1_kp1_mass_total->GetMaximum() / 4,S1_kp1_mass_total->GetMaximum()); // amplitude for second gauss
   // func1->SetParLimits(4,0.485,0.505); // mean for second gauss
   // func1->SetParLimits(5,0.005,0.04); // sigma for second gauss

   // // Strangeness 2 - kaon 1
   // func1_s2_kp1->SetParameters(S2_kp1_mass_total->GetMaximum() / 2,0.493,0.02);
   // func1_s2_kp1->SetParameter(7,S2_kp1_mass_total->GetMaximum() / 2);
   // func1_s2_kp1->SetParameter(8,0.493);
   // func1_s2_kp1->SetParameter(9,0.02);
   // // Setting parameter limits before fitting
   // func1_s2_kp1->SetParLimits(0,S2_kp1_mass_total->GetMaximum() / 3,S2_kp1_mass_total->GetMaximum());
   // func1_s2_kp1->SetParLimits(1,0.480,0.505);
   // func1_s2_kp1->SetParLimits(2,0.005,0.03);
   // func1_s2_kp1->SetParLimits(7,S2_kp1_mass_total->GetMaximum() / 3,S2_kp1_mass_total->GetMaximum());
   // func1_s2_kp1->SetParLimits(8,0.480,0.505);
   // func1_s2_kp1->SetParLimits(9,0.005,0.05);
   //
   // // Strangeness 3 - kaon 1
   // func1_s3_kp1->SetParameters(S3_kp1_mass_total->GetMaximum() / 2,0.493,0.02);
   // func1_s3_kp1->SetParameter(7,S3_kp1_mass_total->GetMaximum() / 2);
   // func1_s3_kp1->SetParameter(8,0.493);
   // func1_s3_kp1->SetParameter(9,0.02);
   // // Setting parameter limits before fitting
   // func1_s3_kp1->SetParLimits(0,S3_kp1_mass_total->GetMaximum() / 3,S3_kp1_mass_total->GetMaximum());
   // func1_s3_kp1->SetParLimits(1,0.480,0.505);
   // func1_s3_kp1->SetParLimits(2,0.005,0.03);
   // func1_s3_kp1->SetParLimits(7,S3_kp1_mass_total->GetMaximum() / 3,S3_kp1_mass_total->GetMaximum());
   // func1_s3_kp1->SetParLimits(8,0.480,0.505);
   // func1_s3_kp1->SetParLimits(9,0.005,0.05);

   // Fitting functions and getting parameters
   // Strangeness 1 - kaon 1
   S1_kp1_mass_total->Fit("func1","RB");
   func2->FixParameter(0, func1->GetParameter(0));
   func2->FixParameter(1, func1->GetParameter(1));
   func2->FixParameter(2, func1->GetParameter(2));
   // func3->FixParameter(0, func1->GetParameter(3));
   // func3->FixParameter(1, func1->GetParameter(4));
   // func3->FixParameter(2, func1->GetParameter(5));
   // func4->FixParameter(0, func1->GetParameter(0));
   // func4->FixParameter(1, func1->GetParameter(1));
   // func4->FixParameter(2, func1->GetParameter(2));
   // func4->FixParameter(3, func1->GetParameter(3));
   // func4->FixParameter(4, func1->GetParameter(4));
   // func4->FixParameter(5, func1->GetParameter(5));
   func5->FixParameter(0, func1->GetParameter(3));
   func5->FixParameter(1, func1->GetParameter(4));
   func5->FixParameter(2, func1->GetParameter(5));
   func5->FixParameter(3, func1->GetParameter(6));
   // func5->FixParameter(4, func1->GetParameter(10));
   // func5->FixParameter(5, func1->GetParameter(11));

   // // Strangeness 2 - kaon 1
   // S2_kp1_mass_total->Fit("func1_s2_kp1","RB");
   // func2_s2_kp1->FixParameter(0, func1_s2_kp1->GetParameter(0));
   // func2_s2_kp1->FixParameter(1, func1_s2_kp1->GetParameter(1));
   // func2_s2_kp1->FixParameter(2, func1_s2_kp1->GetParameter(2));
   // func3_s2_kp1->FixParameter(0, func1_s2_kp1->GetParameter(3));
   // func3_s2_kp1->FixParameter(1, func1_s2_kp1->GetParameter(4));
   // func3_s2_kp1->FixParameter(2, func1_s2_kp1->GetParameter(5));
   // func3_s2_kp1->FixParameter(3, func1_s2_kp1->GetParameter(6));
   // func4_s2_kp1->FixParameter(0, func1_s2_kp1->GetParameter(7));
   // func4_s2_kp1->FixParameter(1, func1_s2_kp1->GetParameter(8));
   // func4_s2_kp1->FixParameter(2, func1_s2_kp1->GetParameter(9));
   // func5_s2_kp1->FixParameter(0, func1_s2_kp1->GetParameter(0));
   // func5_s2_kp1->FixParameter(1, func1_s2_kp1->GetParameter(1));
   // func5_s2_kp1->FixParameter(2, func1_s2_kp1->GetParameter(2));
   // func5_s2_kp1->FixParameter(3, func1_s2_kp1->GetParameter(7));
   // func5_s2_kp1->FixParameter(4, func1_s2_kp1->GetParameter(8));
   // func5_s2_kp1->FixParameter(5, func1_s2_kp1->GetParameter(9));
   //
   // // Strangeness 3 - kaon 1
   // S3_kp1_mass_total->Fit("func1_s3_kp1","RB");
   // func2_s3_kp1->FixParameter(0, func1_s3_kp1->GetParameter(0));
   // func2_s3_kp1->FixParameter(1, func1_s3_kp1->GetParameter(1));
   // func2_s3_kp1->FixParameter(2, func1_s3_kp1->GetParameter(2));
   // func3_s3_kp1->FixParameter(0, func1_s3_kp1->GetParameter(3));
   // func3_s3_kp1->FixParameter(1, func1_s3_kp1->GetParameter(4));
   // func3_s3_kp1->FixParameter(2, func1_s3_kp1->GetParameter(5));
   // func3_s3_kp1->FixParameter(3, func1_s3_kp1->GetParameter(6));
   // func4_s3_kp1->FixParameter(0, func1_s3_kp1->GetParameter(7));
   // func4_s3_kp1->FixParameter(1, func1_s3_kp1->GetParameter(8));
   // func4_s3_kp1->FixParameter(2, func1_s3_kp1->GetParameter(9));
   // func5_s3_kp1->FixParameter(0, func1_s3_kp1->GetParameter(0));
   // func5_s3_kp1->FixParameter(1, func1_s3_kp1->GetParameter(1));
   // func5_s3_kp1->FixParameter(2, func1_s3_kp1->GetParameter(2));
   // func5_s3_kp1->FixParameter(3, func1_s3_kp1->GetParameter(7));
   // func5_s3_kp1->FixParameter(4, func1_s3_kp1->GetParameter(8));
   // func5_s3_kp1->FixParameter(5, func1_s3_kp1->GetParameter(9));

   func5->SetLineColor(kGreen);


   S1_kp1_mass_total->Draw();
   func1->Draw("same");
   func2->Draw("same");
   // func3->Draw("same");
   // func4->Draw("same");
   func5->Draw("same");

}
