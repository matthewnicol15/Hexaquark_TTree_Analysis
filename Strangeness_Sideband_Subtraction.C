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

  // Setting the strings for canvas name
  File_Path<<"/media/mn688/Elements1/PhD/Analysis_Output/Strangeness_Analysis/Sideband_Subtracted/";
  Data<<"RGA_Spring2019_Inbending_dst_Tree_Total_Proton_Smear";
  Quantity<<"";
  Date<<"30112021";
  Version<<"01";

  // Setting the output file name
  Output_File_Name<<File_Path.str().c_str()<<"Sideband_Strangeness_Analysis_"<<Data.str().c_str()<<Quantity.str().c_str()<<"_"<<Date.str().c_str()<<"_"<<Version.str().c_str()<<".root";


  //////////////////////////////////////////////////////////////////////////////
  //// Setting up files and histograms    ////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////

  // Input file
  // RGA data
  // TFile *f1=new TFile("/media/mn688/Elements1/PhD/Analysis_Output/Strangeness_Analysis/Strangeness_Analysis_RGA_Spring2019_Inbending_dst_Tree_Total_22112021_01.root");
  // RGB data
  // TFile *f1=new TFile("/media/mn688/Elements1/PhD/Analysis_Output/Strangeness_Analysis/Strangeness_Analysis_RGB_Spring2020_Inbending_dst_Tree_Total_23112021_01.root");
  // RGA data with proton smearing
  TFile *f1=new TFile("/media/mn688/Elements1/PhD/Analysis_Output/Strangeness_Analysis/Strangeness_Analysis__proton_smear_RGA_Spring2019_Inbending_dst_Tree_Total__25112021_01.root");


  // Getting the multidimensional histogram plots
  TH2D *hmiss_1_a__S1_kp_1 = (TH2D*)f1->Get("hmiss_1_a__S1_kp_1");
  TH3F *hmiss_s2_a__S2_kp_1__S2_kp_2 = (TH3F*)f1->Get("hmiss_s2_a__S2_kp_1__S2_kp_2");
  // Get the number of bins
  Int_t maxbin = hmiss_1_a__S1_kp_1->GetNbinsY();
  Int_t zbinmax = hmiss_s2_a__S2_kp_1__S2_kp_2->GetNbinsZ();


  // Make array for all the 1D projections
  TH1D *h_projectionx_S1_sig[100];
  TH1D *h_projectionx_S1_back[100];
  TH1D *h_projectionx_S2_sig[100];
  TH1D *h_projectionx_S2_back[100];


  // Creating output file
  TFile *output_file=new TFile(Output_File_Name.str().c_str(),"recreate");


  //////////////////////////////////////////////////////////////////////////////
  //// Take all the projections    ///////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////

  // Strangeness 1
  // missing mass of electron and 1 kaon
  TH1D *S1_miss_mass_total = (TH1D*)hmiss_1_a__S1_kp_1->ProjectionX()->Clone("S1_miss_mass_total");
  // Mass of kaon 1
  TH1D *S1_kp1_mass_total = (TH1D*)hmiss_1_a__S1_kp_1->ProjectionY()->Clone("S1_kp1_mass_total");

  // Strangeness 2 file 1
  // missing mass of electron and 2 kaons
  TH1D *S2_miss_mass_total = (TH1D*)hmiss_s2_a__S2_kp_1__S2_kp_2->Project3D("x")->Clone("S2_miss_mass_total");
  // Mass of kaon 1
  TH1D *S2_kp1_mass_total = (TH1D*)hmiss_s2_a__S2_kp_1__S2_kp_2->Project3D("y")->Clone("S2_kp1_mass_total");
  // Mass of kaon 2
  TH1D *S2_kp2_mass_total = (TH1D*)hmiss_s2_a__S2_kp_1__S2_kp_2->Project3D("z")->Clone("S2_kp2_mass_total");


  //////////////////////////////////////////////////////////////////////////////
  //// Creating and fitting functions    ///////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////

  // Create sigma limits for sidebands
  // Strangeness 1 - kaon 1
  Double_t S1_Signal_Sigma_Limits, S1_Background_Sigma_Lower_Limit, S1_Background_Sigma_Upper_Limit;
  // Strangeness 2 - kaon 1
  Double_t S2_Signal_Sigma_Limits, S2_Background_Sigma_Lower_Limit, S2_Background_Sigma_Upper_Limit;

  // Define sigma limits for sidebands
  // Strangeness 1 - kaon 1
  S1_Peak_Sigma_Limits = 2;
  S1_Background_Sigma_Lower_Limit = 4;
  S1_Background_Sigma_Upper_Limit = 6;
  // Strangeness 2 - kaon 1
  S2_Peak_Sigma_Limits = 2;
  S2_Background_Sigma_Lower_Limit = 4;
  S2_Background_Sigma_Upper_Limit = 6;


  // Define functions for fitting kaon calculated mass
  // Functions for strangeness 1 - kaon 1
  TF1 *func1 = new TF1("func1","gaus(0) + pol3(3) + gaus(7)",0.36,0.7);
  TF1 *func2 = new TF1("func2","gaus(0)",0.36,0.7);
  TF1 *func3 = new TF1("func3","pol3(0)",0.36,0.7);
  TF1 *func4 = new TF1("func4","gaus(0)",0.36,0.7);
  TF1 *func5 = new TF1("func5","gaus(0) + gaus(3)",0.36,0.7);

  // Functions for strangeness 2 - kaon 1 file 1
  TF1 *func1_s2_kp1 = new TF1("func1_s2_kp1","gaus(0) + pol3(3) + gaus(7)",0.36,0.7);
  TF1 *func2_s2_kp1 = new TF1("func2_s2_kp1","gaus(0)",0.36,0.7);
  TF1 *func3_s2_kp1 = new TF1("func3_s2_kp1","pol3(0)",0.36,0.7);
  TF1 *func4_s2_kp1 = new TF1("func4_s2_kp1","gaus(0)",0.36,0.7);
  TF1 *func5_s2_kp1 = new TF1("func5_s2_kp1","gaus(0) + gaus(3)",0.36,0.7);


  // Setting parameters before fitting
  // Strangeness 1 - kaon 1
  func1->SetParameters(S1_kp1_mass_total->GetMaximum(),0.493,0.02); // Amplitude, mean, sigma for firs gauss
  func1->SetParameter(7,S1_kp1_mass_total->GetMaximum() / 2); // amplitude for second gauss
  func1->SetParameter(8,0.493); // mean for second gauss
  func1->SetParameter(9,0.02); // sigma for second gauss
  // Setting parameter limits before fitting
  func1->SetParLimits(0,S1_kp1_mass_total->GetMaximum() / 3,S1_kp1_mass_total->GetMaximum()); // amplitude for first gauss
  func1->SetParLimits(1,0.480,0.505); // mean for first gauss
  func1->SetParLimits(2,0.005,0.03); // sigma for first gauss
  func1->SetParLimits(7,S1_kp1_mass_total->GetMaximum() / 3,S1_kp1_mass_total->GetMaximum()); // amplitude for second gauss
  func1->SetParLimits(8,0.480,0.505); // mean for second gauss
  func1->SetParLimits(9,0.005,0.05); // sigma for second gauss

  // Strangeness 2 - kaon 1 file 1
  func1_s2_kp1->SetParameters(S2_kp1_mass_total->GetMaximum() / 2,0.493,0.02);
  func1_s2_kp1->SetParameter(7,S2_kp1_mass_total->GetMaximum() / 2);
  func1_s2_kp1->SetParameter(8,0.493);
  func1_s2_kp1->SetParameter(9,0.02);
  // Setting parameter limits before fitting
  func1_s2_kp1->SetParLimits(0,S2_kp1_mass_total->GetMaximum() / 3,S2_kp1_mass_total->GetMaximum());
  func1_s2_kp1->SetParLimits(1,0.480,0.505);
  func1_s2_kp1->SetParLimits(2,0.005,0.03);
  func1_s2_kp1->SetParLimits(7,S2_kp1_mass_total->GetMaximum() / 3,S2_kp1_mass_total->GetMaximum());
  func1_s2_kp1->SetParLimits(8,0.480,0.505);
  func1_s2_kp1->SetParLimits(9,0.005,0.05);


  // Fitting functions and getting parameters
  // Strangeness 1 - kaon 1
  S1_kp1_mass_total->Fit("func1","RB");
  func2->FixParameter(0, func1->GetParameter(0));
  func2->FixParameter(1, func1->GetParameter(1));
  func2->FixParameter(2, func1->GetParameter(2));
  func3->FixParameter(0, func1->GetParameter(3));
  func3->FixParameter(1, func1->GetParameter(4));
  func3->FixParameter(2, func1->GetParameter(5));
  func3->FixParameter(3, func1->GetParameter(6));
  func4->FixParameter(0, func1->GetParameter(7));
  func4->FixParameter(1, func1->GetParameter(8));
  func4->FixParameter(2, func1->GetParameter(9));
  func5->FixParameter(0, func1->GetParameter(0));
  func5->FixParameter(1, func1->GetParameter(1));
  func5->FixParameter(2, func1->GetParameter(2));
  func5->FixParameter(3, func1->GetParameter(7));
  func5->FixParameter(4, func1->GetParameter(8));
  func5->FixParameter(5, func1->GetParameter(9));

  // Strangeness 2 - kaon 1 file 1
  S2_kp1_mass_total->Fit("func1_s2_kp1","RB");
  func2_s2_kp1->FixParameter(0, func1_s2_kp1->GetParameter(0));
  func2_s2_kp1->FixParameter(1, func1_s2_kp1->GetParameter(1));
  func2_s2_kp1->FixParameter(2, func1_s2_kp1->GetParameter(2));
  func3_s2_kp1->FixParameter(0, func1_s2_kp1->GetParameter(3));
  func3_s2_kp1->FixParameter(1, func1_s2_kp1->GetParameter(4));
  func3_s2_kp1->FixParameter(2, func1_s2_kp1->GetParameter(5));
  func3_s2_kp1->FixParameter(3, func1_s2_kp1->GetParameter(6));
  func4_s2_kp1->FixParameter(0, func1_s2_kp1->GetParameter(7));
  func4_s2_kp1->FixParameter(1, func1_s2_kp1->GetParameter(8));
  func4_s2_kp1->FixParameter(2, func1_s2_kp1->GetParameter(9));
  func5_s2_kp1->FixParameter(0, func1_s2_kp1->GetParameter(0));
  func5_s2_kp1->FixParameter(1, func1_s2_kp1->GetParameter(1));
  func5_s2_kp1->FixParameter(2, func1_s2_kp1->GetParameter(2));
  func5_s2_kp1->FixParameter(3, func1_s2_kp1->GetParameter(7));
  func5_s2_kp1->FixParameter(4, func1_s2_kp1->GetParameter(8));
  func5_s2_kp1->FixParameter(5, func1_s2_kp1->GetParameter(9));


  //////////////////////////////////////////////////////////////////////////////
  //// Creating and define sideband limits    ///////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////

  // Sideband limits
  // Strangeness 1 - kaon 1
  Double_t S1_Peak_Lower_Limit, S1_Peak_Upper_Limit;
  Double_t S1_Background_Left_Lower_Limit, S1_Background_Left_Upper_Limit;
  Double_t S1_Background_Right_Lower_Limit, S1_Background_Right_Upper_Limit;
  // Strangeness 2 - kaon 1 file 1
  Double_t S2_Peak_Lower_Limit, S2_Peak_Upper_Limit;
  Double_t S2_Background_Left_Lower_Limit, S2_Background_Left_Upper_Limit;
  Double_t S2_Background_Right_Lower_Limit, S2_Background_Right_Upper_Limit;


  // Define Sideband Limits
  // Strangeness 1 - kaon 1
  S1_Peak_Lower_Limit = func1->GetParameter(1) - S1_Peak_Sigma_Limits * func1->GetParameter(2);
  S1_Peak_Upper_Limit = func1->GetParameter(1) + S1_Peak_Sigma_Limits * func1->GetParameter(2);
  S1_Background_Left_Lower_Limit = func1->GetParameter(1) - S1_Background_Sigma_Upper_Limit * func1->GetParameter(2);
  S1_Background_Left_Upper_Limit = func1->GetParameter(1) - S1_Background_Sigma_Lower_Limit * func1->GetParameter(2);
  S1_Background_Right_Lower_Limit = func1->GetParameter(1) + S1_Background_Sigma_Lower_Limit * func1->GetParameter(2);
  S1_Background_Right_Upper_Limit = func1->GetParameter(1) + S1_Background_Sigma_Upper_Limit * func1->GetParameter(2);


  // Strangeness 2 - kaon 1 file 1
  S2_Peak_Lower_Limit = func1_s2_kp1->GetParameter(1) - S2_Peak_Sigma_Limits * func1_s2_kp1->GetParameter(2);
  S2_Peak_Upper_Limit = func1_s2_kp1->GetParameter(1) + S2_Peak_Sigma_Limits * func1_s2_kp1->GetParameter(2);
  S2_Background_Left_Lower_Limit = func1_s2_kp1->GetParameter(1) - S2_Background_Sigma_Upper_Limit * func1_s2_kp1->GetParameter(2);
  S2_Background_Left_Upper_Limit = func1_s2_kp1->GetParameter(1) - S2_Background_Sigma_Lower_Limit * func1_s2_kp1->GetParameter(2);
  S2_Background_Right_Lower_Limit = func1_s2_kp1->GetParameter(1) + S2_Background_Sigma_Lower_Limit * func1_s2_kp1->GetParameter(2);
  S2_Background_Right_Upper_Limit = func1_s2_kp1->GetParameter(1) + S2_Background_Sigma_Upper_Limit * func1_s2_kp1->GetParameter(2);

  //////////////////////////////////////////////////////////////////////////////
  //// Calculating scaling factor    //////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////

  Double_t S1_Scaling_Factor_sig, S1_Scaling_Factor_back, S2_Scaling_Factor_sig, S2_Scaling_Factor_back;
  //
  // // Strangeness 1 - kaon 1
  // // Loop over the x bins to determine scaling factors
  for(Int_t bin = 1; bin < 100; bin++){
    ostringstream h_projectionx_S1_sig_name;
     h_projectionx_S1_sig_name << "h_projectionx_S1_sig_" << bin;
    // Get the x projection for the current bin
    h_projectionx_S1_sig[bin] = (TH1D*)hmiss_1_a__S1_kp_1->ProjectionX("",bin,bin)->Clone(h_projectionx_S1_sig_name.str().c_str());
    h_projectionx_S1_back[bin] = (TH1D*)hmiss_1_a__S1_kp_1->ProjectionX("",bin,bin)->Clone("h_projectionx_S1_back");
    cout<<bin<<endl;

    if(h_projectionx_S1_sig[bin]->Integral() > 100){
      cout<<"integral big"<<endl;
      // Determine the scaling factor looking at the signal function
      S1_Scaling_Factor_sig = func5->Eval(S1_kp1_mass_total->GetBinCenter(bin)) / h_projectionx_S1_sig[bin]->Integral() ;


      // Determine the scaling factor looking at the signal function
      S1_Scaling_Factor_back = (func3->Eval(S1_kp1_mass_total->GetBinCenter(bin)) - func5->Eval(S1_kp1_mass_total->GetBinCenter(bin))) /
      h_projectionx_S1_back[bin]->Integral();

      if(S1_Scaling_Factor_back < 0) S1_Scaling_Factor_back = 0;

      // Scale the signal and background histograms accordingly
      h_projectionx_S1_sig[bin]->Scale(S1_Scaling_Factor_sig);
      h_projectionx_S1_back[bin]->Scale(S1_Scaling_Factor_back);
    }
  }
  // Adding all the histograms together to get total continous sideband subtracted result
  for(Int_t bin = 2; bin < 100; bin++){
    cout<<bin<<endl;
    h_projectionx_S1_sig[1]->Add(h_projectionx_S1_sig[bin]);
    h_projectionx_S1_back[1]->Add(h_projectionx_S1_back[bin]);
  }
  Double_t back_integral_S1 = h_projectionx_S1_back[1]->Integral(h_projectionx_S1_back[1]->FindBin(0), h_projectionx_S1_back[1]->FindBin(1));
  Double_t sig_integral_S1 = h_projectionx_S1_sig[1]->Integral(h_projectionx_S1_sig[1]->FindBin(0), h_projectionx_S1_sig[1]->FindBin(1));
  h_projectionx_S1_back[1]->Scale(sig_integral_S1 / back_integral_S1);

  h_projectionx_S1_sig[1]->Add(h_projectionx_S1_back[1], -1);



  // Strangeness 2 - kaon 1 file 1
  // Loop over the x bins to determine scaling factors
  for(Int_t bin = 1; bin < 100; bin++){
    ostringstream h_projectionx_S2_sig_name;
     h_projectionx_S2_sig_name << "h_projectionx_S2_sig_" << bin;

    // Get the x projection for the current bin
    h_projectionx_S2_sig[bin] = (TH1D*)hmiss_s2_a__S2_kp_1__S2_kp_2->ProjectionX("",bin,bin,0,zbinmax)->Clone(h_projectionx_S2_sig_name.str().c_str());
    h_projectionx_S2_back[bin] = (TH1D*)hmiss_s2_a__S2_kp_1__S2_kp_2->ProjectionX("",bin,bin,0,zbinmax)->Clone("h_projectionx_S2_back");
    cout<<bin<<endl;

    if(h_projectionx_S2_sig[bin]->Integral() > 100){
      cout<<"integral big"<<endl;
      // Determine the scaling factor looking at the signal function
      S2_Scaling_Factor_sig = func5_s2_kp1->Eval(S2_kp1_mass_total->GetBinCenter(bin)) / h_projectionx_S2_sig[bin]->Integral() ;


      // Determine the scaling factor looking at the signal function
      S2_Scaling_Factor_back = (func3_s2_kp1->Eval(S2_kp1_mass_total->GetBinCenter(bin)) - func5_s2_kp1->Eval(S2_kp1_mass_total->GetBinCenter(bin)))/h_projectionx_S2_back[bin]->Integral();
      if(S2_Scaling_Factor_back < 0) S2_Scaling_Factor_back = 0;

      // Scale the signal and background histograms accordingly
      h_projectionx_S2_sig[bin]->Scale(S2_Scaling_Factor_sig);
      h_projectionx_S2_back[bin]->Scale(S2_Scaling_Factor_back);
    }
  }
  // Adding all the histograms together to get total continous sideband subtracted result
  for(Int_t bin = 2; bin < 100; bin++){
    cout<<bin<<endl;
    h_projectionx_S2_sig[1]->Add(h_projectionx_S2_sig[bin]);
    h_projectionx_S2_back[1]->Add(h_projectionx_S2_back[bin]);
  }
  Double_t back_integral_S2 = h_projectionx_S2_back[1]->Integral(h_projectionx_S2_back[1]->FindBin(0), h_projectionx_S2_back[1]->FindBin(1));
  Double_t sig_integral_S2 = h_projectionx_S2_sig[1]->Integral(h_projectionx_S2_sig[1]->FindBin(0), h_projectionx_S2_sig[1]->FindBin(1));
  h_projectionx_S2_back[1]->Scale(sig_integral_S2 / back_integral_S2);

  h_projectionx_S2_sig[1]->Add(h_projectionx_S2_back[1], -1);


  ////////////////////////////////
  //////////////////////////////////////////////
  //// Sideband Subtraction       //////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////

  // Defining the limits of the peak and backgroud regions
  // Strangeness 1 - kaon 1
  Int_t Lower_end = S1_kp1_mass_total->FindBin(func1->GetParameter(1) - 2 * func1->GetParameter(2));
  Int_t Upper_end = S1_kp1_mass_total->FindBin(func1->GetParameter(1) + 2 * func1->GetParameter(2));

  // Strangeness 2 - kaon 1
  Int_t Peak_Lower_end_s2 = S2_kp1_mass_total->FindBin(S2_Peak_Lower_Limit);
  Int_t Peak_Upper_end_s2 = S2_kp1_mass_total->FindBin(S2_Peak_Upper_Limit);
  Int_t Back_Left_Lower_end_s2 = S2_kp1_mass_total->FindBin(S2_Background_Left_Lower_Limit);
  Int_t Back_Left_Upper_end_s2 = S2_kp1_mass_total->FindBin(S2_Background_Left_Upper_Limit);
  Int_t Back_Right_Lower_end_s2 = S2_kp1_mass_total->FindBin(S2_Background_Right_Lower_Limit);
  Int_t Back_Right_Upper_end_s2 = S2_kp1_mass_total->FindBin(S2_Background_Right_Upper_Limit);

  // Get the missing mass in the kaon mass peak region
  TH1D *S2_miss_mass_peak = (TH1D*)hmiss_s2_a__S2_kp_1__S2_kp_2->ProjectionX
  ("_px",Peak_Lower_end_s2, Peak_Upper_end_s2,0,zbinmax)->Clone("S2_miss_mass_peak");

  // Get the missing mass in the kaon mass background left region
  TH1D *S2_miss_mass_background_left = (TH1D*)hmiss_s2_a__S2_kp_1__S2_kp_2->ProjectionX
  ("",Back_Left_Lower_end_s2, Back_Left_Upper_end_s2,0,zbinmax)->Clone("S2_miss_mass_background_left");

  // Get the missing mass in the kaon mass background right region
  TH1D *S2_miss_mass_background_Right = (TH1D*)hmiss_s2_a__S2_kp_1__S2_kp_2->ProjectionX
  ("",Back_Right_Lower_end_s2, Back_Right_Upper_end_s2,0,zbinmax)->Clone("S2_miss_mass_background_Right");

  // Combining the two backgroud sidebands
  TH1D *S2_miss_mass_background_Total = (TH1D*)S2_miss_mass_background_left->Clone("S2_miss_mass_background_Total");
  S2_miss_mass_background_Total->Add(S2_miss_mass_background_Right);

  TH1D *S2_miss_mass_background_Total_continous = (TH1D*)S2_miss_mass_background_Total->Clone("S2_miss_mass_background_Total_continous");


  // Creating variables for multiplication factors
  // Strangeness 2 - kaon 1
  Double_t Peak_Background, Left_Background, Right_Background;
  Double_t Multiplication_factor;

  // Determining the multiplication factor
  // Strangeness 2
  Peak_Background = func3_s2_kp1->Integral(S2_miss_mass_total->FindBin(S2_Peak_Lower_Limit),S2_miss_mass_total->FindBin(S2_Peak_Upper_Limit));

  Left_Background = func3_s2_kp1->Integral(S2_miss_mass_total->FindBin(S2_Background_Left_Lower_Limit),S2_miss_mass_total->FindBin(S2_Background_Left_Upper_Limit));

  Right_Background = func3_s2_kp1->Integral(S2_miss_mass_total->FindBin(S2_Background_Right_Lower_Limit),S2_miss_mass_total->FindBin(S2_Background_Right_Upper_Limit));

  Multiplication_factor = Peak_Background / (Left_Background + Right_Background);

  S2_miss_mass_background_Total->Scale(Multiplication_factor);
  // cout<<Peak_Background<<" "<<Left_Background<<" "<<Right_Background<<" "<<Multiplication_factor<<endl;

  // Producing the sideband subtracted result
  TH1D *S2_miss_mass_Result = (TH1D*)S2_miss_mass_peak->Clone("S2_miss_mass_Result");
  S2_miss_mass_Result->Add(S2_miss_mass_background_Total,-1);

  //////////////////////////////////////////////////////////////////////////////
  //// Continous Sideband Subtraction     /////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////

  // Get the missing mass from limits of sidebands
  TH1D *S2_miss_mass_total_continous = (TH1D*)hmiss_s2_a__S2_kp_1__S2_kp_2->ProjectionX
  ("",S2_Background_Left_Lower_Limit, S2_Background_Right_Upper_Limit,0,zbinmax)->Clone("S2_miss_mass_total_continous");
  // Get the missing mass from limitis of sidebands for result plot
  TH1D *S2_miss_mass_total_continous_result = (TH1D*)hmiss_s2_a__S2_kp_1__S2_kp_2->ProjectionX
  ("",S2_Background_Left_Lower_Limit, S2_Background_Right_Upper_Limit,0,zbinmax)->Clone("S2_miss_mass_total_continous_result");

  // Determine integral of sidebands
  Double_t Sideband_Left_Integral = func1_s2_kp1->Integral(S2_Background_Left_Lower_Limit, S2_Background_Left_Upper_Limit);
  Double_t Sideband_Right_Integral = func1_s2_kp1->Integral(S2_Background_Right_Lower_Limit, S2_Background_Right_Upper_Limit);
  Double_t Sideband_Integral = Sideband_Left_Integral + Sideband_Right_Integral;
  // Determine integral of background function from sideband limits
  Double_t Background_Total_Integral = func3_s2_kp1->Integral
  (S2_Background_Left_Lower_Limit,
    S2_Background_Right_Upper_Limit);
    // Double_t Background_Total_Integral = func3_s2_kp1->Integral
    // (S2_miss_mass_total->FindBin(S2_Background_Left_Lower_Limit),
    // S2_miss_mass_total->FindBin(S2_Background_Right_Upper_Limit));


    // Scale background distribution to entire background function
    S2_miss_mass_background_Total_continous->Scale(Background_Total_Integral / Sideband_Integral);
    // S2_miss_mass_background_Total_continous->Scale(1.5);

    // Fit total background histogram
    TF1 *S2_back_total_func1 = new TF1("S2_back_total_func1","pol5(0)");
    S2_miss_mass_background_Total_continous->Fit("S2_back_total_func1");
    TH1F *hbackground_total_function = new TH1F("hbackground_total_function","total background function",300,0,3);

    // Loop over bins and set bin content according to background function
    for(int back_bin = 1; back_bin < 300; back_bin++){
      hbackground_total_function->SetBinContent(back_bin,
        S2_back_total_func1->Eval(hbackground_total_function->GetBinCenter(back_bin)));
    }

    // Take away background from total distribution
    S2_miss_mass_total_continous_result->Add(hbackground_total_function, -1);
    // S2_miss_mass_total_continous_result->Add(S2_miss_mass_background_Total_continous, -1);

    //////////////////////////////////////////////////////////////////////////////
    //// Styling histograms and plots    //////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////

    // Matching gstyle stuff for CAA
    gROOT->Reset();
    TStyle *plain  = new TStyle("Plain","Plain Style (no colors/fill areas)");

    plain->SetCanvasBorderMode(0);
    plain->SetPadBorderMode(0);
    plain->SetPadColor(0);
    plain->SetCanvasColor(0);
    plain->SetTitleColor(0);
    plain->SetStatColor(0);
    plain->SetPalette(1);
    gROOT->SetStyle("Plain");
    gStyle->SetPalette(1);
    gStyle->SetPalette(1);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadBottomMargin(0.15);
    //gStyle->SetTitleFontSize(0.15);
    gStyle->SetOptStat(0);
    gStyle->SetLineStyleString(13,"10 30");  // dotted
    gStyle->SetLineStyleString(12,"40 40"); // dashed
    gStyle->SetLineStyleString(14,"50 25 10 25");

    // Set colour for the total backgroud histogram
    hbackground_total_function->SetLineColor(kRed);

    // Matching style of other plots in CAA for the continous sideband subtraction
    // result with my method
    S2_miss_mass_total_continous_result->SetTitle("Background Subtracted Missing Mass");
    S2_miss_mass_total_continous_result->Rebin(4);
    S2_miss_mass_total_continous_result->SetLineWidth(3);
    // S2_miss_mass_total_continous_result->SetLineColor(4);
    S2_miss_mass_total_continous_result->SetLineColor(kRed);
    S2_miss_mass_total_continous_result->SetMarkerStyle(8);
    S2_miss_mass_total_continous_result->SetMarkerSize(1);
    S2_miss_mass_total_continous_result->SetMarkerColor(4);
    S2_miss_mass_total_continous_result->GetXaxis()->SetTitleOffset(1.2);
    // S2_miss_mass_total_continous_result->GetYaxis()->SetTitle("Counts");
    // S2_miss_mass_total_continous_result->GetXaxis()->SetTitle("P_{n} [GeV/c]");


    h_projectionx_S2_sig[1]->SetTitle("Background Subtracted Missing Mass");
    h_projectionx_S2_sig[1]->Rebin(4);
    h_projectionx_S2_sig[1]->SetLineWidth(3);
    h_projectionx_S2_sig[1]->SetLineColor(4);
    h_projectionx_S2_sig[1]->SetMarkerStyle(8);
    h_projectionx_S2_sig[1]->SetMarkerSize(1);
    h_projectionx_S2_sig[1]->SetMarkerColor(4);
    h_projectionx_S2_sig[1]->GetXaxis()->SetTitleOffset(1.2);
    h_projectionx_S2_sig[1]->GetYaxis()->SetTitle("Counts");

    //////////////////////////////////////////////////////////////////////////////
    //// Creating lines to show sidebands and resonances    /////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////

    // Signal lower end
    auto *l1 = new TLine(S2_Peak_Lower_Limit, 0,S2_Peak_Lower_Limit, S2_kp1_mass_total->GetMaximum());
    // Signal upper end
    auto *l2 = new TLine(S2_Peak_Upper_Limit, 0,S2_Peak_Upper_Limit, S2_kp1_mass_total->GetMaximum());
    // background left lower end
    auto *l3 = new TLine(S2_Background_Left_Lower_Limit, 0,S2_Background_Left_Lower_Limit, S2_kp1_mass_total->GetMaximum());
    // background left upper end
    auto *l4 = new TLine(S2_Background_Left_Upper_Limit, 0,S2_Background_Left_Upper_Limit, S2_kp1_mass_total->GetMaximum());
    // background right lower end
    auto *l5 = new TLine(S2_Background_Right_Lower_Limit, 0,S2_Background_Right_Lower_Limit, S2_kp1_mass_total->GetMaximum());
    // background right upper end
    auto *l6 = new TLine(S2_Background_Right_Upper_Limit, 0,S2_Background_Right_Upper_Limit, S2_kp1_mass_total->GetMaximum());

    // line at 0 to highlight background subtraction
    auto *l7 = new TLine(0, 0, 3, 0);
    l7->SetLineColor(kRed);

    // Lines for strangeness 2
    // cascade ground state
    auto *l8 = new TLine(1.321, 0, 1.321, h_projectionx_S2_sig[1]->GetMaximum());
    auto *l9 = new TLine(1.535, 0, 1.535, h_projectionx_S2_sig[1]->GetMaximum());
    auto *l10 = new TLine(1.614, 0, 1.614, h_projectionx_S2_sig[1]->GetMaximum() * 0.8);

    // l8->SetLineWidth(2);
    // l9->SetLineWidth(2);
    l8->SetLineStyle(9);
    l9->SetLineStyle(9);
    l10->SetLineStyle(9);
    l8->SetLineColor(15);
    l9->SetLineColor(15);
    l10->SetLineColor(15);

    ////////////////////////////////////////////////////////////////////////////
    //// Creating text boxes and arrows    /////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////

    TPaveLabel *Cascade_ground_state = new TPaveLabel(1.25,250,1.35,280,"#Xi^{-}");
    Cascade_ground_state->SetFillColor(0);
    Cascade_ground_state->SetTextSize(0.8);

    TPaveLabel *Cascade_1530 = new TPaveLabel(1.5,280,1.6,310,"#Xi(1530)^{-}");
    Cascade_1530->SetFillColor(0);
    Cascade_1530->SetTextSize(0.8);

    TPaveLabel *Threshold = new TPaveLabel(1.62,10,2.4,40,"K^{+}#Lambda K^{+} K^{-}");
    Threshold->SetFillColor(0);
    Threshold->SetTextSize(0.6);

    TArrow *threshold_arrow = new TArrow(1.614,5,2.2,5,0.05,">");
    threshold_arrow->SetAngle(40);
    threshold_arrow->SetLineWidth(2);

    //////////////////////////////////////////////////////////////////////////////
    //// Creating canvases and drawing plots    //////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////

    // canvas for projections
    auto *c1 = new TCanvas("c1","original and backgroud subtracted",800,800);
    c1->cd();
    S2_miss_mass_total_continous_result->Draw("hist,L");
    h_projectionx_S2_sig[1]->Draw("same,hist,L");

    // S2_miss_mass_total->Draw();
    // hbackground_total_function->Draw("hist,same");
    // S2_kp1_mass_total->Draw();
    // // func1_s2_kp1->Draw("same");
    // func2_s2_kp1->Draw("same");
    // func3_s2_kp1->Draw("same");
    // func4_s2_kp1->Draw("same");
    // // func5_s2_kp1->Draw("same");
    // l1->Draw("same");
    // l2->Draw("same");
    // l3->Draw("same");
    // l4->Draw("same");
    // l5->Draw("same");
    // l6->Draw("same");

    auto *c2 = new TCanvas("c2","original and backgroud subtracted",800,800);
    c2->cd();
    h_projectionx_S1_sig[1]->Draw("same,hist,L");


    auto *c4 = new TCanvas("c4","strangeness 2 peak",800,800);
    c4->cd();
    // h_projectionx_S2_back[1]->Draw("hist");
    h_projectionx_S2_sig[1]->Draw("same,hist,L");
    l7->Draw("same");
    l8->Draw("same");
    l9->Draw("same");
    l10->Draw("same");
    Cascade_ground_state->Draw("same");
    Cascade_1530->Draw("same");
    Threshold->Draw("same");
    threshold_arrow->Draw();

    output_file->Write();
  }
