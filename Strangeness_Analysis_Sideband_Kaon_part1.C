#include <cstdlib>
#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TF1.h>
#include <TH2.h>
#include <TH3.h>
#include <TMath.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <TChain.h>
#include <TBenchmark.h>
#include <vector>

// Macro name
void Strangeness_Analysis_Sideband_Kaon_part1(){

  // Beam energies
  // RGA Fall 2018
  // Double_t beam_energy = 10.604;
  // RGA Spring 2019
  Double_t beam_energy = 10.199;
  // RGB Spring 2020, for runno less than 11394
  // Double_t beam_energy = 10.2129;
  // RGB Spring 2020, for runno 11394 and above
  // Double_t beam_energy = 10.3894;

  //////////////////////////////////////////////////////////////////////////////
  ////Define variables for naming and limits ///////////
  //////////////////////////////////////////////////////////////////////////////

  // Information for canvas and histogram name
  ostringstream Additional_information;
  ostringstream Data;
  ostringstream Quantity;
  ostringstream Date;
  ostringstream Version;
  ostringstream Output_File_Name;

  // Setting the strings for canvas name
  Additional_information<<"";
  Data<<"RGA_Spring2019_Inbending_dst_Tree_Total_201021";
  Quantity<<"";
  Date<<"09122021";
  Version<<"03";

  Output_File_Name<<"/media/mn688/Elements1/PhD/Analysis_Output/Strangeness_Analysis_"<<Additional_information.str().c_str()<<Data.str().c_str()<<"_"<<Quantity.str().c_str()<<"_"<<Date.str().c_str()<<"_"<<Version.str().c_str()<<".root";


  //////////////////////////////////////////////////////////////////////////////
  ////Creating components to read from TTree ///////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////

  // Read input root file and assign it to 'f'
  TFile *f = new TFile("/shared/storage/physhad/JLab/mn688/Trees/Dibaryon/RGA/RGA_Spring2019_Inbending_at_least_1eFD_1Kp_Tree_Total_201021.root");

  // Read TTree within root file and assign it to 't1'
  TTree *t1 = (TTree*)f->Get("RGA_Spring2019_Inbending_201021");

  // Read file with information on vectors
  gROOT->ProcessLine(".L ~/work/Macros/Loader.C+");

  // Event information
  TLorentzVector *readbeam=NULL;  // Information on the beam
  TLorentzVector *readtarget=NULL; // Information on the target
  TLorentzVector target; // Information on the target
  Double_t start_time; // Event start time
  Int_t readrunno;
  Int_t readeventno;
  Int_t readtriggerno;
  // Number of given particle or charge track in each event
  Int_t readchargetracks; // Number of positive or negative charge tracks
  Int_t readelno; // e^-
  Int_t readothertracks; // Number of particles excluding p, e^-, pions or kaons
  Int_t region; // which region the particles go in (FT, FD, CD)

  // Set any vectors to 0 before reading from the TTree
  // Particle information
  vector<TLorentzVector> *v_p4=0;   // 4-vectors for the detected particles
  vector<TLorentzVector> *v_vertex=0;   // Vertex information for particles
  vector<double> *v_path=0;   // Measured path of particles
  vector<double> *v_time=0;   // Measured time of flight of particles
  vector<double> *v_beta=0;   // Beta measured from FTOF
  vector<double> *v_energy=0;   // Energy of particle
  vector<double> *v_charge=0;   // Charge of particle from Drift Chambers
  vector<double> *v_PID=0;   // PID from FTOF
  vector<double> *v_chi2PID=0;   // Chi^2 of the PID
  vector<Int_t> *v_region=0; // region particle goes in

  // Setting the branch addresses to read from
  t1->SetBranchAddress("p4",&v_p4);
  t1->SetBranchAddress("vertex",&v_vertex);
  t1->SetBranchAddress("path",&v_path);
  t1->SetBranchAddress("beta",&v_beta);
  t1->SetBranchAddress("beam",&readbeam);
  t1->SetBranchAddress("target",&readtarget);
  t1->SetBranchAddress("start_time",&start_time);
  t1->SetBranchAddress("energy",&v_energy);
  t1->SetBranchAddress("charge",&v_charge);
  t1->SetBranchAddress("PID",&v_PID);
  t1->SetBranchAddress("chi2PID",&v_chi2PID);
  t1->SetBranchAddress("region",&v_region);
  t1->SetBranchAddress("time",&v_time);
  t1->SetBranchAddress("elno",&readelno);
  t1->SetBranchAddress("eventno",&readeventno);
  t1->SetBranchAddress("runno",&readrunno);
  t1->SetBranchAddress("triggerno",&readtriggerno);

  // Path and name for the output file to save
  TFile fileOutput1(Output_File_Name.str().c_str(),"recreate");


  // Getting particle database to use for masses
  auto db=TDatabasePDG::Instance();

  ////////////////////////////////////////////////////////////////////////////////
  ////Creating functions for kaon mass fit ///////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////

  // Define functions for fitting kaon calculated mass
  // Function for strangeness 1 - kaon 1
  TF1 *func1 = new TF1("func1","gaus(0) + pol3(3) + gaus(7)",0.36,0.7);
  TF1 *func2 = new TF1("func2","gaus(0)",0.36,0.7);
  TF1 *func3 = new TF1("func3","pol3(0)",0.36,0.7);
  TF1 *func4 = new TF1("func4","gaus(0)",0.36,0.7);
  TF1 *func5 = new TF1("func5","gaus(0) + gaus(3)",0.36,0.7);

  // Function for strangeness 2 - kaon 1
  TF1 *func1_s2_kp1 = new TF1("func1_s2_kp1","gaus(0) + pol3(3) + gaus(7)",0.36,0.7);
  TF1 *func2_s2_kp1 = new TF1("func2_s2_kp1","gaus(0)",0.36,0.7);
  TF1 *func3_s2_kp1 = new TF1("func3_s2_kp1","pol3(0)",0.36,0.7);
  TF1 *func4_s2_kp1 = new TF1("func4_s2_kp1","gaus(0)",0.36,0.7);
  TF1 *func5_s2_kp1 = new TF1("func5_s2_kp1","gaus(0) + gaus(3)",0.36,0.7);

  // Function for strangeness 3 - kaon 1

  //////////////////////////////////////////////////////////////////////////////
  ////Creating functions for proton smearing on RGA data   /////////////////////
  //////////////////////////////////////////////////////////////////////////////

  // Random number for phi distribution

  // Sin distribution for theta
  TF1 *theta_function = new TF1("theta_function","sin(x)",0,TMath::Pi());
  // Momentum distribution for implementing fermi motion
  TF1 *fermi_momentum_function = new TF1("fermi_momentum_function",
  "pow((x*(pow(0.26,2)-pow(0.0456,2))/(pow(x,2)+pow(0.0456,2))/(pow(x,2)+pow(0.26,2))),2)",0.0,2.0);

  TRandom *phi_function = new TRandom3(); // Random number generator to produce mass of dsss

  Double_t generated_theta, generated_phi, generated_momentum;
  Double_t generated_Px, generated_Py, generated_Pz;
  Double_t generated_mass;

  //////////////////////////////////////////////////////////////////////////////
  ////Create histograms here ///////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////

  // Histograms for events
  auto* hbeam=new TH1D("hbeam","Beam mass; Beam Mass [GeV];Counts",200,0,11);
  auto* hkaon=new TH1D("hkaon","kaon momentum; kaon momentum [GeV];Counts",200,0,10);
  auto* hkaons=new TH1D("hkaons","kaon numbers; Kaons in event;Counts",6,0,6);
  auto* hproton=new TH1D("hproton","proton momentum; proton momentum [GeV];Counts",200,0,10);
  auto* hregion=new TH1D("hregion","Regions;Region;Counts",3,1,4);
  auto* h_photon_energy_S1 = new TH2F("h_photon_energy_S1","Photon energy for S1 vs calculated kaon mass",100,0.2,0.8,120,0,12);
  auto* h_photon_energy_S2 = new TH2F("h_photon_energy_S2","Photon energy for S2 vs calculated kaon mass",100,0.2,0.8,120,0,12);
  auto* h_photon_energy_S3 = new TH2F("h_photon_energy_S3","Photon energy for S3 vs calculated kaon mass",100,0.2,0.8,120,0,12);

  // Histograms looking at kaon properties
  auto* hmass_S1_kp_1=new TH1F("hmass_S1_kp_1","K^{+} mass;M(K^{+});Counts",100,0.2,0.8);
  auto* hmass_S2_kp_1=new TH1F("hmass_S2_kp_1","K^{+} mass;M(K^{+});Counts",100,0.2,0.8);
  auto* hmass_S2_kp_2=new TH1F("hmass_S2_kp_2","K^{+} mass;M(K^{+});Counts",100,0.2,0.8);
  auto* hmass_S3_kp_1=new TH1F("hmass_S3_kp_1","K^{+} mass;M(K^{+});Counts",100,0.2,0.8);
  auto* hmass_S3_kp_2=new TH1F("hmass_S3_kp_2","K^{+} mass;M(K^{+});Counts",100,0.2,0.8);
  auto* hmass_S3_kp_3=new TH1F("hmass_S3_kp_3","K^{+} mass;M(K^{+});Counts",100,0.2,0.8);
  auto* hmass_S3_kp_1_a=new TH1F("hmass_S3_kp_1_a","K^{+} mass;M(K^{+});Counts",100,0.2,0.8);
  auto* hmass_S3_kp_2_a=new TH1F("hmass_S3_kp_2_a","K^{+} mass;M(K^{+});Counts",100,0.2,0.8);
  auto* hmass_S3_kp_3_a=new TH1F("hmass_S3_kp_3_a","K^{+} mass;M(K^{+});Counts",100,0.2,0.8);
  auto* h_theta_s1_kp = new TH1F("h_theta_s1_kp","K^{+} #Theta distribution",180,0,180);
  auto* h_theta_s2_kp1 = new TH1F("h_theta_s2_kp1","K^{+} #Theta distribution",180,0,180);
  auto* h_theta_s2_kp2 = new TH1F("h_theta_s2_kp2","K^{+} #Theta distribution",180,0,180);
  auto* h_theta_s3_kp1 = new TH1F("h_theta_s3_kp1","K^{+} #Theta distribution",180,0,180);
  auto* h_theta_s3_kp2 = new TH1F("h_theta_s3_kp2","K^{+} #Theta distribution",180,0,180);
  auto* h_theta_s3_kp3 = new TH1F("h_theta_s3_kp3","K^{+} #Theta distribution",180,0,180);

  // Histograms for strangeness 1 channel
  auto* hmiss_mass_all=new TH1D("miss_all","MM^2(e' K^{+} p #pi^{-});MM^2(e' K^{+} p #pi^{-}) [GeV];Counts",200,-1,1);
  auto* hmiss_mass_all_a=new TH1D("hmiss_mass_all_a","MM^2(e' K^{+} p #pi^{-});MM^2(e' K^{+} p #pi^{-}) [GeV];Counts",200,-1,1);
  auto* hmiss_momentum_all=new TH1D("hmiss_momentum_all","P(B + T - e' - K^{+} - p - #pi^{-});P(B + T - e' - K^{+} - p - #pi^{-}) [GeV];Counts",200,0,2);
  auto* hmiss_momentum_all_a=new TH1D("hmiss_momentum_all_a","P(B + T - e' - K^{+} - p - #pi^{-});P(B + T - e' - K^{+} - p - #pi^{-}) [GeV];Counts",200,0,2);
  auto* hinv_lambda=new TH1D("hinv_lambda","Invariant mass of p #pi^{-};M(p #pi^{-}) [GeV];Counts",200,0.5,2.5);
  auto* hinv_lambda_a=new TH1D("hinv_lambda_a","Invariant mass of p #pi^{-};M(p #pi^{-}) [GeV];Counts",200,0.5,2.5);
  auto* hmiss_1=new TH1D("hmiss_1","MM(e' K^{+});MM(e' K^{+}) [GeV];Counts",200,0,4);
  auto* hmiss_1_a__S1_kp_1=new TH2D("hmiss_1_a__S1_kp_1","Kaon mass against missing mass;MM(e' K^{+}) [GeV]; M(K^{+}) [GeV]",200,0,4,100,0.2,0.8);
  auto* hmiss_1_b=new TH1D("hmiss_1_b","MM(e' K^{+});MM(e' K^{+}) [GeV];Counts",200,0,4);
  auto* hmiss_1_c=new TH1D("hmiss_1_c","MM(e' K^{+});MM(e' K^{+}) [GeV];Counts",200,0,4);
  auto* h_delta_beta_kp_s1_1=new TH2D("h_delta_beta_kp_s1_1","#Delta #Beta K^{+}; P [GeV]; #Delta #Beta",200,0,11,200,-1,1);
  auto* h_delta_beta_kp_s1_1FD=new TH2D("h_delta_beta_kp_s1_1FD","#Delta #Beta K^{+}; P [GeV]; #Delta #Beta",200,0,11,200,-1,1);

  // Histograms for strangeness 2 channel
  auto* h_delta_beta_kp_s2_1=new TH2D("h_delta_beta_kp_s2_1","#Delta #Beta K^{+}; P [GeV]; #Delta #Beta",200,0,11,200,-1,1);
  auto* h_delta_beta_kp_s2_2=new TH2D("h_delta_beta_kp_s2_2","#Delta #Beta K^{+}; P [GeV]; #Delta #Beta",200,0,11,200,-1,1);
  auto* h_delta_beta_kp_s2_1FD=new TH2D("h_delta_beta_kp_s2_1FD","#Delta #Beta K^{+}; P [GeV]; #Delta #Beta",200,0,11,200,-1,1);
  auto* h_delta_beta_kp_s2_2FD=new TH2D("h_delta_beta_kp_s2_2FD","#Delta #Beta K^{+}; P [GeV]; #Delta #Beta",200,0,11,200,-1,1);
  auto* hmiss_2=new TH1D("hmiss_2","MM^{2}(e' K^{+} p);MM^{2}(e' K^{+} p) [GeV^{2}];Counts",200,-1,3);
  auto* hmiss_s2=new TH1D("hmiss_s2","MM(e' K^{+} K^{+});MM(e' K^{+} K^{+}) [GeV];Counts",300,0,3);
  auto* hmiss_s2_a__S2_kp_1__S2_kp_2=new TH3F("hmiss_s2_a__S2_kp_1__S2_kp_2",
  "MM against M(K^{+}) (1) against M(K^{+}) (2);MM(e' K^{+} K^{+}) [GeV];M(K^{+}) (1) [GeV]; M(K^{+}) (2) [GeV]",
  300,0,3,100,0.2,0.8,100,0.2,0.8);

  // Histograms for strangeness 3 channel
  auto* hmiss_s3=new TH1D("hmiss_s3","MM(e' K^{+} K^{+} K^{+});MM(e' K^{+} K^{+} K^{+}) [GeV];Counts",300,0,3);
  auto* hmiss_s3_a=new TH1D("hmiss_s3_a","MM(e' K^{+} K^{+} K^{+});MM(e' K^{+} K^{+} K^{+}) [GeV];Counts",300,0,3);
  auto* h_delta_beta_kp_s3_1=new TH2D("h_delta_beta_kp_s3_1","#Delta #Beta K^{+}; P [GeV]; #Delta #Beta",200,0,11,200,-1,1);
  auto* h_delta_beta_kp_s3_2=new TH2D("h_delta_beta_kp_s3_2","#Delta #Beta K^{+}; P [GeV]; #Delta #Beta",200,0,11,200,-1,1);
  auto* h_delta_beta_kp_s3_3=new TH2D("h_delta_beta_kp_s3_3","#Delta #Beta K^{+}; P [GeV]; #Delta #Beta",200,0,11,200,-1,1);
  auto* h_delta_beta_kp_s3_1FD=new TH2D("h_delta_beta_kp_s3_1FD","#Delta #Beta K^{+}; P [GeV]; #Delta #Beta",200,0,11,200,-1,1);
  auto* h_delta_beta_kp_s3_2FD=new TH2D("h_delta_beta_kp_s3_2FD","#Delta #Beta K^{+}; P [GeV]; #Delta #Beta",200,0,11,200,-1,1);
  auto* h_delta_beta_kp_s3_3FD=new TH2D("h_delta_beta_kp_s3_3FD","#Delta #Beta K^{+}; P [GeV]; #Delta #Beta",200,0,11,200,-1,1);
  auto* hmiss_s3_a__S3_kp_1=new TH2F("hmiss_s3_a__S3_kp_1",
  "MM against M(K^{+}) (1) ;MM(e' K^{+} K^{+} K^{+}) [GeV];M(K^{+}) (1) [GeV]",
  300,0,3,100,0.2,0.8);

  // Create vectors of TLorentzVectors to store information of
  // all particles of a given type (important when you have more than 1
  // of any particle in the same event)
  vector<TLorentzVector> v_el;  // e^-
  vector<TLorentzVector> v_pip; // pi^+
  vector<TLorentzVector> v_pim; // pi^-
  vector<TLorentzVector> v_pr; // protons
  vector<TLorentzVector> v_kp; // K^+
  vector<TLorentzVector> v_km; // K^-
  vector<TLorentzVector> v_unidentified_neg; // Particles with a PID of 0
  vector<TLorentzVector> v_othertracks; // Any other particles are assigned to this

  // TLorentzVectors for individual particles
  TLorentzVector el;
  TLorentzVector pip;
  TLorentzVector pim;
  TLorentzVector pr;
  TLorentzVector kp;
  TLorentzVector km;
  TLorentzVector othertracks;
  TLorentzVector unidentified_neg;

  // Vertex information for different particle types
  vector<TLorentzVector> v_vertex_el; // e^-
  vector<TLorentzVector> v_vertex_pip; // pi^+
  vector<TLorentzVector> v_vertex_pim; // pi^-
  vector<TLorentzVector> v_vertex_pr; // protons
  vector<TLorentzVector> v_vertex_kp; // K^+
  vector<TLorentzVector> v_vertex_km; // K^-
  TLorentzVector vertex_el;
  TLorentzVector vertex_pip;
  TLorentzVector vertex_pim;
  TLorentzVector vertex_pr;
  TLorentzVector vertex_kp;
  TLorentzVector vertex_km;

  // These are used to define the missing masses later
  TLorentzVector beam;
  TLorentzVector missall;
  TLorentzVector miss1;
  TLorentzVector miss_s2, miss_s3;
  TLorentzVector miss2;
  TLorentzVector lambda;


  // After information is read from the TTree, particles are identified using
  // PID and assigned to that particle type
  // e^-
  vector<Double_t> v_beta_tof_el;  // Beta measured
  Double_t beta_tof_el;
  vector<Double_t> v_P_el;  // Momentum measured
  Double_t P_el;
  vector<Double_t>v_path_el; // Path measured
  Double_t path_el;
  vector<Double_t>v_TOF_el;  // TOF measured
  Double_t TOF_el;
  vector<Double_t> v_beta_calc_el;  // Beta calculated
  Double_t beta_calc_el;
  vector<Double_t> v_delta_beta_el;  // Beta calculated - beta measured
  Double_t delta_beta_el;
  vector<Double_t> v_vertex_time_el;  // Vertex time calculated
  Double_t vertex_time_el;
  vector<Double_t> v_region_el; // region hit
  Double_t region_el;

  // pi^+
  vector<Double_t> v_beta_tof_pip;  // Beta measured
  Double_t beta_tof_pip;
  vector<Double_t> v_P_pip;  // Momentum measured
  Double_t P_pip;
  vector<Double_t>v_path_pip; // Path measured
  Double_t path_pip;
  vector<Double_t>v_TOF_pip;  // TOF measured
  Double_t TOF_pip;
  vector<Double_t> v_beta_calc_pip;  // Beta calculated
  Double_t beta_calc_pip;
  vector<Double_t> v_delta_beta_pip;  // Beta calculated - beta measured
  Double_t delta_beta_pip;
  vector<Double_t> v_vertex_time_pip;  // Vertex time calculated
  Double_t vertex_time_pip;
  vector<Double_t> v_region_pip; // region hit
  Double_t region_pip;

  // pi^-
  vector<Double_t> v_beta_tof_pim;  // Beta measured
  Double_t beta_tof_pim;
  vector<Double_t> v_P_pim;  // Momentum measured
  Double_t P_pim;
  vector<Double_t>v_path_pim; // Path measured
  Double_t path_pim;
  vector<Double_t>v_TOF_pim;  // TOF measured
  Double_t TOF_pim;
  vector<Double_t> v_beta_calc_pim;  // Beta calculated
  Double_t beta_calc_pim;
  vector<Double_t> v_delta_beta_pim;  // Beta calculated - beta measured
  Double_t delta_beta_pim;
  vector<Double_t> v_vertex_time_pim;  // Vertex time calculated
  Double_t vertex_time_pim;
  vector<Double_t> v_region_pim; // region hit
  Double_t region_pim;

  // protons
  vector<Double_t> v_beta_tof_pr;  // Beta measured
  Double_t beta_tof_pr;
  vector<Double_t> v_P_pr;  // Momentum measured
  Double_t P_pr;
  vector<Double_t>v_path_pr; // Path measured
  Double_t path_pr;
  vector<Double_t>v_TOF_pr;  // TOF measured
  Double_t TOF_pr;
  vector<Double_t> v_beta_calc_pr;  // Beta calculated
  Double_t beta_calc_pr;
  vector<Double_t> v_delta_beta_pr;  // Beta calculated - beta measured
  Double_t delta_beta_pr;
  vector<Double_t> v_vertex_time_pr;  // Vertex time calculated
  Double_t vertex_time_pr;
  vector<Double_t> v_region_pr; // region hit
  Double_t region_pr;

  // K^+
  vector<Double_t> v_beta_tof_kp;  // Beta measured
  Double_t beta_tof_kp;
  vector<Double_t> v_P_kp;  // Momentum measured
  Double_t P_kp;
  vector<Double_t>v_path_kp; // Path measured
  Double_t path_kp;
  vector<Double_t>v_TOF_kp;  // TOF measured
  Double_t TOF_kp;
  vector<Double_t> v_beta_calc_kp;  // Beta calculated
  Double_t beta_calc_kp;
  vector<Double_t> v_delta_beta_kp;  // Beta calculated - beta measured
  Double_t delta_beta_kp;
  vector<Double_t> v_vertex_time_kp;  // Vertex time calculated
  Double_t vertex_time_kp;
  vector<Double_t> v_region_kp; // region hit
  Double_t region_kp;
  vector<Double_t> v_Mass_kp; // region hit
  Double_t mass_kp;

  // K^-
  vector<Double_t> v_beta_tof_km;  // Beta measured
  Double_t beta_tof_km;
  vector<Double_t> v_P_km;  // Momentum measured
  Double_t P_km;
  vector<Double_t>v_path_km; // Path measured
  Double_t path_km;
  vector<Double_t>v_TOF_km;  // TOF measured
  Double_t TOF_km;
  vector<Double_t> v_beta_calc_km;  // Beta calculated
  Double_t beta_calc_km;
  vector<Double_t> v_delta_beta_km;  // Beta calculated - beta measured
  Double_t delta_beta_km;
  vector<Double_t> v_vertex_time_km;  // Vertex time calculated
  Double_t vertex_time_km;
  vector<Double_t> v_region_km; // region hit
  Double_t region_km;

  // Particles with PID = 0
  vector<Double_t> v_beta_tof_unidentified_neg;  // Beta measured
  Double_t beta_tof_unidentified_neg;
  vector<Double_t> v_P_unidentified_neg;  // Momentum measured
  Double_t P_unidentified_neg;
  vector<Double_t>v_path_unidentified_neg; // Path measured
  Double_t path_unidentified_neg;
  vector<Double_t>v_TOF_unidentified_neg;  // TOF measured
  Double_t TOF_unidentified_neg;
  vector<Double_t> v_beta_calc_unidentified_neg;  // Beta calculated
  Double_t beta_calc_unidentified_neg;
  vector<Double_t> v_delta_beta_unidentified_neg;  // Beta calculated - beta measured
  Double_t delta_beta_unidentified_neg;
  vector<Double_t> v_vertex_time_unidentified_neg;  // Vertex time calculated
  Double_t vertex_time_unidentified_neg;
  vector<Double_t> v_region_unidentified_neg; // region hit
  Double_t region_unidentified_neg;

  Double_t c=30;  // Speed of light used for calculating vertex time

  // Reads the total number of entries in the TTree
  Long64_t nentries = t1->GetEntries();
  // You can just run over a set number of events for fast analysis
  // Long64_t nentries = 100000;
  cout<<nentries<<endl; // Printing out the total number of events

  // This is used to print out the percentage of events completed so far
  Int_t Percentage = nentries/100;

  // This loops over all the entries in the TTree
  for(Long64_t i=0; i<nentries;i++){
    t1->GetEntry(i);

    // This prints out the percentage of events completed so far
    if (i % Percentage == 0){
      fprintf (stderr, "%lld\r", i/Percentage);
      fflush (stderr);
    }
    // All the vectors must be cleared at the start of each event entry
    // e^-
    v_el.clear();
    v_beta_tof_el.clear();
    v_P_el.clear();
    v_path_el.clear();
    v_TOF_el.clear();
    v_beta_calc_el.clear();
    v_delta_beta_el.clear();
    v_vertex_time_el.clear();
    v_vertex_el.clear();
    v_region_el.clear();

    // pi^+
    v_pip.clear();
    v_beta_tof_pip.clear();
    v_P_pip.clear();
    v_path_pip.clear();
    v_TOF_pip.clear();
    v_beta_calc_pip.clear();
    v_delta_beta_pip.clear();
    v_vertex_pip.clear();
    v_vertex_time_pip.clear();
    v_vertex_pip.clear();
    v_region_pip.clear();

    // pi^-
    v_pim.clear();
    v_beta_tof_pim.clear();
    v_P_pim.clear();
    v_path_pim.clear();
    v_TOF_pim.clear();
    v_beta_calc_pim.clear();
    v_delta_beta_pim.clear();
    v_vertex_pim.clear();
    v_vertex_time_pim.clear();
    v_vertex_pim.clear();
    v_region_pim.clear();

    // protons
    v_pr.clear();
    v_beta_tof_pr.clear();
    v_P_pr.clear();
    v_path_pr.clear();
    v_TOF_pr.clear();
    v_beta_calc_pr.clear();
    v_delta_beta_pr.clear();
    v_vertex_pr.clear();
    v_vertex_time_pr.clear();
    v_vertex_pr.clear();
    v_region_pr.clear();

    // K^+
    v_kp.clear();
    v_beta_tof_kp.clear();
    v_P_kp.clear();
    v_path_kp.clear();
    v_TOF_kp.clear();
    v_beta_calc_kp.clear();
    v_delta_beta_kp.clear();
    v_vertex_kp.clear();
    v_vertex_time_kp.clear();
    v_vertex_kp.clear();
    v_region_kp.clear();
    v_Mass_kp.clear();

    // K^-
    v_km.clear();
    v_beta_tof_km.clear();
    v_P_km.clear();
    v_path_km.clear();
    v_TOF_km.clear();
    v_beta_calc_km.clear();
    v_delta_beta_km.clear();
    v_vertex_km.clear();
    v_vertex_time_km.clear();
    v_vertex_km.clear();
    v_region_km.clear();

    // Particles with PID = 0
    v_unidentified_neg.clear();
    v_beta_tof_unidentified_neg.clear();
    v_P_unidentified_neg.clear();
    v_path_unidentified_neg.clear();
    v_TOF_unidentified_neg.clear();
    v_beta_calc_unidentified_neg.clear();
    v_delta_beta_unidentified_neg.clear();
    // v_vertex_unidentified_neg.clear();
    v_vertex_time_unidentified_neg.clear();
    // v_vertex_unidentified_neg.clear();
    v_region_unidentified_neg.clear();

    // Other particles
    // v_othertracks.clear();
    // v_beta_tof_othertracks.clear();
    // v_P_othertracks.clear();
    // v_beta_calc_othertracks.clear();
    // v_delta_beta_othertracks.clear();

    // This reads the number of particles in the current entry
    Int_t Nparticles = v_p4->size();

    // This loops over all the particles in the current entry
    for(Int_t j=0; j<Nparticles; j++){

      // Filling histogram to show the different regions hit
      hregion->Fill(v_region->at(j));
      // Can put a selection on which regions particles are hitting

      // Checking PID and assigning particles
      // e^-
      if(v_PID->at(j)==11){
        // Setting the 4-vector and assigning mass from PDG
        el.SetXYZM(v_p4->at(j).Px(),v_p4->at(j).Py(),v_p4->at(j).Pz(),db->GetParticle(11)->Mass());
        TOF_el = v_time->at(j); // Measured time
        path_el = v_path->at(j); // Measured path
        beta_tof_el = v_beta->at(j); // Measured beta from FTOF
        // Calculating beta from momentum and mass
        beta_calc_el = el.Rho()/(sqrt((pow(el.Rho(),2))+(pow(el.M(),2))));
        // Difference between calculated and measured beta
        delta_beta_el = beta_calc_el-beta_tof_el;
        vertex_time_el = TOF_el - path_el / (beta_tof_el*c); // Calculate vertex time
        // Setting the vertex information now vertex time has been calculated
        vertex_el.SetXYZT(v_vertex->at(j).X(), v_vertex->at(j).Y(), v_vertex->at(j).Z(), vertex_time_el);
        region_el = v_region->at(j);

        // Pushing back all that iformation into the vectors
        // Again this is done so you can store information on multiple particles
        // of the same type in one place
        v_el.push_back(el);
        v_beta_tof_el.push_back(beta_tof_el);
        v_P_el.push_back(P_el);
        v_path_el.push_back(path_el);
        v_TOF_el.push_back(TOF_el);
        v_beta_calc_el.push_back(beta_calc_el);
        v_delta_beta_el.push_back(delta_beta_el);
        v_vertex_time_el.push_back(vertex_time_el);
        v_vertex_el.push_back(vertex_el);
        v_region_el.push_back(region_el);
      }

      // pi^+
      else if(v_PID->at(j)==211){
        // Setting the 4-vector and assigning mass from PDG
        pip.SetXYZM(v_p4->at(j).Px(),v_p4->at(j).Py(),v_p4->at(j).Pz(),db->GetParticle(211)->Mass());
        TOF_pip = v_time->at(j); // Measured time
        path_pip = v_path->at(j); // Measured path
        beta_tof_pip = v_beta->at(j); // Measured beta from FTOF
        // Calculating beta from momentum and mass
        beta_calc_pip = pip.Rho()/(sqrt((pow(pip.Rho(),2))+(pow(pip.M(),2))));
        // Difference between calculated and measured beta
        delta_beta_pip = beta_calc_pip-beta_tof_pip;
        vertex_time_pip = TOF_pip - path_pip / (beta_tof_pip*c); // Calculate vertex time
        // Setting the vertex information now vertex time has been calculated
        vertex_pip.SetXYZT(v_vertex->at(j).X(), v_vertex->at(j).Y(), v_vertex->at(j).Z(), vertex_time_pip);
        region_pip = v_region->at(j);

        // Pushing back all that iformation into the vectors
        // Again this is done so you can store information on multiple particles
        // of the same type in one place
        v_pip.push_back(pip);
        v_beta_tof_pip.push_back(beta_tof_pip);
        v_P_pip.push_back(P_pip);
        v_path_pip.push_back(path_pip);
        v_TOF_pip.push_back(TOF_pip);
        v_beta_calc_pip.push_back(beta_calc_pip);
        v_delta_beta_pip.push_back(delta_beta_pip);
        v_vertex_time_pip.push_back(vertex_time_pip);
        v_vertex_pip.push_back(vertex_pip);
        v_region_pip.push_back(region_pip);
      }

      // pi^-
      else if(v_PID->at(j)==-211){
        // Setting the 4-vector and assigning mass from PDG
        pim.SetXYZM(v_p4->at(j).Px(),v_p4->at(j).Py(),v_p4->at(j).Pz(),db->GetParticle(-211)->Mass());
        TOF_pim = v_time->at(j); // Measured time
        path_pim = v_path->at(j); // Measured path
        beta_tof_pim = v_beta->at(j); // Measured beta from FTOF
        // Calculating beta from momentum and mass
        beta_calc_pim = pim.Rho()/(sqrt((pow(pim.Rho(),2))+(pow(pim.M(),2))));
        // Difference between calculated and measured beta
        delta_beta_pim = beta_calc_pim-beta_tof_pim;
        vertex_time_pim = TOF_pim - path_pim / (beta_tof_pim*c); // Calculate vertex time
        // Setting the vertex information now vertex time has been calculated
        vertex_pim.SetXYZT(v_vertex->at(j).X(), v_vertex->at(j).Y(), v_vertex->at(j).Z(), vertex_time_pim);
        region_pim = v_region->at(j);

        // Pushing back all that iformation into the vectors
        // Again this is done so you can store information on multiple particles
        // of the same type in one place
        v_pim.push_back(pim);
        v_beta_tof_pim.push_back(beta_tof_pim);
        v_P_pim.push_back(P_pim);
        v_path_pim.push_back(path_pim);
        v_TOF_pim.push_back(TOF_pim);
        v_beta_calc_pim.push_back(beta_calc_pim);
        v_delta_beta_pim.push_back(delta_beta_pim);
        v_vertex_time_pim.push_back(vertex_time_pim);
        v_vertex_pim.push_back(vertex_pim);
        v_region_pim.push_back(region_pim);
      }
      else if(v_PID->at(j)==0 && v_charge->at(j)<0){
        // Setting the 4-vector and assigning mass from PDG
        unidentified_neg.SetXYZM(v_p4->at(j).Px(),v_p4->at(j).Py(),v_p4->at(j).Pz(),db->GetParticle(-211)->Mass());
        TOF_unidentified_neg = v_time->at(j); // Measured time
        path_unidentified_neg = v_path->at(j); // Measured path
        beta_tof_unidentified_neg = v_beta->at(j); // Measured beta from FTOF
        // Calculating beta from momentum and mass
        beta_calc_unidentified_neg = unidentified_neg.Rho()/(sqrt((pow(unidentified_neg.Rho(),2))+(pow(unidentified_neg.M(),2))));
        // Difference between calculated and measured beta
        delta_beta_unidentified_neg = beta_calc_unidentified_neg-beta_tof_unidentified_neg;
        vertex_time_unidentified_neg = TOF_unidentified_neg - path_unidentified_neg / (beta_tof_unidentified_neg*c); // Calculate vertex time
        // Setting the vertex information now vertex time has been calculated
        // vertex_unidentified_neg.SetXYZT(v_vertex->at(j).X(), v_vertex->at(j).Y(), v_vertex->at(j).Z(), vertex_time_unidentified_neg);
        region_unidentified_neg = v_region->at(j);

        // Pushing back all that iformation into the vectors
        // Again this is done so you can store information on multiple particles
        // of the same type in one place
        v_unidentified_neg.push_back(unidentified_neg);
        v_beta_tof_unidentified_neg.push_back(beta_tof_unidentified_neg);
        v_P_unidentified_neg.push_back(P_unidentified_neg);
        v_path_unidentified_neg.push_back(path_unidentified_neg);
        v_TOF_unidentified_neg.push_back(TOF_unidentified_neg);
        v_beta_calc_unidentified_neg.push_back(beta_calc_unidentified_neg);
        v_delta_beta_unidentified_neg.push_back(delta_beta_unidentified_neg);
        v_vertex_time_unidentified_neg.push_back(vertex_time_unidentified_neg);
        // v_vertex_unidentified_neg.push_back(vertex_unidentified_neg);
        v_region_unidentified_neg.push_back(region_unidentified_neg);
      }
      // protons
      else if(v_PID->at(j)==2212){
        // Setting the 4-vector and assigning mass from PDG
        pr.SetXYZM(v_p4->at(j).Px(),v_p4->at(j).Py(),v_p4->at(j).Pz(),db->GetParticle(2212)->Mass());
        TOF_pr = v_time->at(j); // Measured time
        path_pr = v_path->at(j); // Measured path
        beta_tof_pr = v_beta->at(j); // Measured beta from FTOF
        // Calculating beta from momentum and mass
        beta_calc_pr = pr.Rho()/(sqrt((pow(pr.Rho(),2))+(pow(pr.M(),2))));
        // Difference between calculated and measured beta
        delta_beta_pr = beta_calc_pr-beta_tof_pr;
        vertex_time_pr = TOF_pr - path_pr / (beta_tof_pr*c); // Calculate vertex time
        // Setting the vertex information now vertex time has been calculated
        vertex_pr.SetXYZT(v_vertex->at(j).X(), v_vertex->at(j).Y(), v_vertex->at(j).Z(), vertex_time_pr);
        region_pr = v_region->at(j);

        // Pushing back all that iformation into the vectors
        // Again this is done so you can store information on multiple particles
        // of the same type in one place
        v_pr.push_back(pr);
        v_beta_tof_pr.push_back(beta_tof_pr);
        v_P_pr.push_back(P_pr);
        v_path_pr.push_back(path_pr);
        v_TOF_pr.push_back(TOF_pr);
        v_beta_calc_pr.push_back(beta_calc_pr);
        v_delta_beta_pr.push_back(delta_beta_pr);
        v_vertex_time_pr.push_back(vertex_time_pr);
        v_vertex_pr.push_back(vertex_pr);
        v_region_pr.push_back(region_pr);
      }

      // K^+
      else if(v_PID->at(j)==321){
        // Setting the 4-vector and assigning mass from PDG
        kp.SetXYZM(v_p4->at(j).Px(),v_p4->at(j).Py(),v_p4->at(j).Pz(),db->GetParticle(321)->Mass());
        TOF_kp = v_time->at(j); // Measured time
        path_kp = v_path->at(j); // Measured path
        beta_tof_kp = v_beta->at(j); // Measured beta from FTOF
        // Calculating beta from momentum and mass
        beta_calc_kp = kp.Rho()/(sqrt((pow(kp.Rho(),2))+(pow(kp.M(),2))));
        // Difference between calculated and measured beta
        delta_beta_kp = beta_calc_kp-beta_tof_kp;
        vertex_time_kp = TOF_kp - path_kp / (beta_tof_kp*c); // Calculate vertex time
        // Setting the vertex information now vertex time has been calculated
        vertex_kp.SetXYZT(v_vertex->at(j).X(), v_vertex->at(j).Y(), v_vertex->at(j).Z(), vertex_time_kp);
        region_kp = v_region->at(j);
        mass_kp = sqrt((pow(v_p4->at(j).Rho(),2) / (pow(beta_tof_kp,2))) - pow(v_p4->at(j).Rho(),2));

        // Pushing back all that iformation into the vectors
        // Again this is done so you can store information on multiple particles
        // of the same type in one place
        v_kp.push_back(kp);
        v_beta_tof_kp.push_back(beta_tof_kp);
        v_P_kp.push_back(P_kp);
        v_path_kp.push_back(path_kp);
        v_TOF_kp.push_back(TOF_kp);
        v_beta_calc_kp.push_back(beta_calc_kp);
        v_delta_beta_kp.push_back(delta_beta_kp);
        v_vertex_time_kp.push_back(vertex_time_kp);
        v_vertex_kp.push_back(vertex_kp);
        v_region_kp.push_back(region_kp);
        v_Mass_kp.push_back(mass_kp);
      }

      // K^-
      else if(v_PID->at(j)==-321){
        // Setting the 4-vector and assigning mass from PDG
        km.SetXYZM(v_p4->at(j).Px(),v_p4->at(j).Py(),v_p4->at(j).Pz(),db->GetParticle(-321)->Mass());
        TOF_km = v_time->at(j); // Measured time
        path_km = v_path->at(j); // Measured path
        beta_tof_km = v_beta->at(j); // Measured beta from FTOF
        // Calculating beta from momentum and mass
        beta_calc_km = km.Rho()/(sqrt((pow(km.Rho(),2))+(pow(km.M(),2))));
        // Difference between calculated and measured beta
        delta_beta_km = beta_calc_km-beta_tof_km;
        vertex_time_km = TOF_km - path_km / (beta_tof_km*c); // Calculate vertex time
        // Setting the vertex information now vertex time has been calculated
        vertex_km.SetXYZT(v_vertex->at(j).X(), v_vertex->at(j).Y(), v_vertex->at(j).Z(), vertex_time_km);
        region_km = v_region->at(j);

        // Pushing back all that iformation into the vectors
        // Again this is done so you can store information on multiple particles
        // of the same type in one place
        v_km.push_back(km);
        v_beta_tof_km.push_back(beta_tof_km);
        v_P_km.push_back(P_km);
        v_path_km.push_back(path_km);
        v_TOF_km.push_back(TOF_km);
        v_beta_calc_km.push_back(beta_calc_km);
        v_delta_beta_km.push_back(delta_beta_km);
        v_vertex_time_km.push_back(vertex_time_km);
        v_vertex_km.push_back(vertex_km);
        v_region_km.push_back(region_km);
      }
    }

    // Setting beam energy to value given at top of macro
    beam.SetXYZM(0,0,beam_energy,0);
    hbeam->Fill(beam.Rho());

    // Determine photon four vector
    TLorentzVector Photon;
    Photon = beam - v_el.at(0);


    // Here you can apply conditions on the events you want to analyse
    if(v_kp.size() > 0 && v_el.size() == 1 &&v_region_el.at(0) == 1){
      target = (TLorentzVector)*readtarget;

      miss1 = beam + (TLorentzVector)*readtarget - v_el.at(0) - v_kp.at(0);

      hkaons->Fill(v_kp.size());

      // Looking at stragneness 1 channel
      if(v_kp.size()==1){
        hmiss_1->Fill(miss1.M());

        h_delta_beta_kp_s1_1->Fill(v_kp.at(0).Rho(),v_delta_beta_kp.at(0));

        if(v_region_kp.at(0)!=1) continue;

        // Looking at the angular distribution for forward going kaons
        h_theta_s1_kp->Fill(v_kp.at(0).Theta() * TMath::RadToDeg());

        // Looking at delta beta vs momentum for the kaon
        h_delta_beta_kp_s1_1FD->Fill(v_kp.at(0).Rho(),v_delta_beta_kp.at(0));

        // Put a cut on the delta beta of the K+
        if(fabs(v_delta_beta_kp.at(0))<0.02){

          // Plot the missing mass against calculated kaon mass
          hmiss_1_a__S1_kp_1->Fill(miss1.M(),v_Mass_kp.at(0));
          // Plot the calculated mass of the kaon
          hmass_S1_kp_1->Fill(v_Mass_kp.at(0));

          // Plot the photon energy distribution
          h_photon_energy_S1->Fill(v_Mass_kp.at(0),Photon.E());


          if(v_pr.size()==1){
            // Select which region you want the particles to go in
            if(v_region_pr.at(0) == 1){
              miss2 = beam + (TLorentzVector)*readtarget - v_el.at(0) - v_kp.at(0) - v_pr.at(0);
              hmiss_2->Fill(miss2.M2());
              hkaon->Fill(v_kp.at(0).Rho());
              hproton->Fill(v_pr.at(0).Rho());


              // Selecting events where the pion is also detected
              if(v_pim.size()==1){
                lambda = v_pr.at(0) + v_pim.at(0);
                missall = beam + (TLorentzVector)*readtarget - v_el.at(0) - v_kp.at(0) - v_pr.at(0) - v_pim.at(0);
                hmiss_mass_all->Fill(missall.M2());
                hmiss_momentum_all->Fill(missall.Rho());
                hinv_lambda->Fill(lambda.M());
                hmiss_1_b->Fill(miss1.M());


                // Cut around the missing mass of all detected particles (neutron mass)
                if(fabs(missall.M2()) < 0.1){
                  hinv_lambda_a->Fill(lambda.M());
                  hmiss_1_c->Fill(miss1.M());
                }
                // Cutting around the invariant mass of lambda
                if(lambda.M() <1.14){
                  hmiss_mass_all_a->Fill(missall.M2());
                  hmiss_momentum_all_a->Fill(missall.Rho());
                }
              }
            }
          }
        }
      }

      // Looking at stragneness 2 channel
      if(v_kp.size()==2){

        ////////////////////////////////////////////////////////////////////////
        //// Implementing fermi motion    //////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////

        // generated_theta = theta_function -> GetRandom(0,TMath::Pi());
        // generated_phi = phi_function -> Uniform(-TMath::Pi(),TMath::Pi());
        // generated_momentum = fermi_momentum_function -> GetRandom();
        // // cout<<"phi "<<generated_phi<<" theta "<<generated_theta<<" P "<<generated_momentum<<endl;
        //
        // generated_Px = generated_momentum * TMath::Sin(generated_theta) * TMath::Cos(generated_phi);
        // generated_Py = generated_momentum * TMath::Sin(generated_theta) * TMath::Sin(generated_phi);
        // generated_Pz = generated_momentum * TMath::Cos(generated_theta);
        //
        // generated_mass = sqrt((pow(1.876,2) / 4) - pow(generated_momentum,2));
        //
        // target.SetXYZM(generated_Px, generated_Py, generated_Pz, generated_mass);

        ////////////////////////////////////////////////////////////////////////


        miss_s2 = beam + target - v_el.at(0) - v_kp.at(0)- v_kp.at(1);

        hmiss_s2->Fill(miss_s2.M());
        h_delta_beta_kp_s2_1->Fill(v_kp.at(0).Rho(),v_delta_beta_kp.at(0));
        h_delta_beta_kp_s2_2->Fill(v_kp.at(1).Rho(),v_delta_beta_kp.at(1));

        if(v_region_kp.at(0) != 1 || v_region_kp.at(1) != 1) continue;
        h_delta_beta_kp_s2_1FD->Fill(v_kp.at(0).Rho(),v_delta_beta_kp.at(0));
        h_delta_beta_kp_s2_2FD->Fill(v_kp.at(1).Rho(),v_delta_beta_kp.at(1));

        // Looking at the angular distribution for forward going kaons
        h_theta_s2_kp1->Fill(v_kp.at(0).Theta() * TMath::RadToDeg());
        h_theta_s2_kp2->Fill(v_kp.at(1).Theta() * TMath::RadToDeg());

        if(fabs(v_delta_beta_kp.at(0))<0.02 && fabs(v_delta_beta_kp.at(1))<0.02){
          hmiss_s2_a__S2_kp_1__S2_kp_2->Fill(miss_s2.M(), v_Mass_kp.at(0), v_Mass_kp.at(1));
          hmass_S2_kp_1->Fill(v_Mass_kp.at(0));
          hmass_S2_kp_2->Fill(v_Mass_kp.at(1));

          // Plot the photon energy distribution
          h_photon_energy_S2->Fill(v_Mass_kp.at(0),Photon.E());


        }
      }
      if(v_kp.size() == 3){
        miss_s3 = beam + (TLorentzVector)*readtarget - v_el.at(0) - v_kp.at(0) - v_kp.at(1) - v_kp.at(2);

        hmiss_s3->Fill(miss_s3.M());
        hmass_S3_kp_1->Fill(v_Mass_kp.at(0));
        hmass_S3_kp_2->Fill(v_Mass_kp.at(1));
        hmass_S3_kp_3->Fill(v_Mass_kp.at(2));
        // Looking at delta beta of all kaons
        h_delta_beta_kp_s3_1->Fill(v_kp.at(0).Rho(),v_delta_beta_kp.at(0));
        h_delta_beta_kp_s3_2->Fill(v_kp.at(1).Rho(),v_delta_beta_kp.at(1));
        h_delta_beta_kp_s3_3->Fill(v_kp.at(2).Rho(),v_delta_beta_kp.at(2));

        if(v_region_kp.at(0) != 1 || v_region_kp.at(1) != 1 || v_region_kp.at(2) != 1) continue;
        h_delta_beta_kp_s3_1FD->Fill(v_kp.at(0).Rho(),v_delta_beta_kp.at(0));
        h_delta_beta_kp_s3_2FD->Fill(v_kp.at(1).Rho(),v_delta_beta_kp.at(1));
        h_delta_beta_kp_s3_3FD->Fill(v_kp.at(2).Rho(),v_delta_beta_kp.at(2));

        // Looking at the angular distribution for forward going kaons
        h_theta_s3_kp1->Fill(v_kp.at(0).Theta() * TMath::RadToDeg());
        h_theta_s3_kp2->Fill(v_kp.at(1).Theta() * TMath::RadToDeg());
        h_theta_s3_kp3->Fill(v_kp.at(2).Theta() * TMath::RadToDeg());

        if(fabs(v_delta_beta_kp.at(0))<0.02 && fabs(v_delta_beta_kp.at(1))<0.02 && fabs(v_delta_beta_kp.at(2))<0.02){
          hmiss_s3_a__S3_kp_1->Fill(miss_s3.M(), v_Mass_kp.at(0));
          hmass_S3_kp_1_a->Fill(v_Mass_kp.at(0));
          hmass_S3_kp_2_a->Fill(v_Mass_kp.at(1));
          hmass_S3_kp_3_a->Fill(v_Mass_kp.at(2));

          // Plot the photon energy distribution
          h_photon_energy_S3->Fill(v_Mass_kp.at(0),Photon.E());
        }
      }
    }
  }

  //////////////////////////////////////////////////////////////////////////////////////
  //// Fitting functions to calculated kaon mass  //////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////

  // Set parameteres for kaon mass fit

  // Strangeness 1 - kaon 1
  // Setting parameters before fitting
  func1->SetParameters(hmass_S1_kp_1->GetMaximum(),0.493,0.02); // Amplitude, mean, sigma for firs gauss
  func1->SetParameter(7,hmass_S1_kp_1->GetMaximum() / 2); // amplitude for second gauss
  func1->SetParameter(8,0.493); // mean for second gauss
  func1->SetParameter(9,0.02); // sigma for second gauss
  // Setting parameter limits before fitting
  func1->SetParLimits(0,hmass_S1_kp_1->GetMaximum() / 3,hmass_S1_kp_1->GetMaximum()); // amplitude for first gauss
  func1->SetParLimits(1,0.480,0.505); // mean for first gauss
  func1->SetParLimits(2,0.005,0.03); // sigma for first gauss
  func1->SetParLimits(7,hmass_S1_kp_1->GetMaximum() / 3,hmass_S1_kp_1->GetMaximum()); // amplitude for second gauss
  func1->SetParLimits(8,0.480,0.505); // mean for second gauss
  func1->SetParLimits(9,0.005,0.05); // sigma for second gauss

  // Strangeness 1 - kaon 1
  hmass_S1_kp_1->Fit("func1","RB");
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


  // Strangeness 2 - kaon 1
  // Setting parameters before fitting
  func1_s2_kp1->SetParameters(hmass_S2_kp_1->GetMaximum() / 2,0.493,0.02);
  func1_s2_kp1->SetParameter(7,hmass_S2_kp_1->GetMaximum() / 2);
  func1_s2_kp1->SetParameter(8,0.493);
  func1_s2_kp1->SetParameter(9,0.02);
  // Setting parameter limits before fitting
  func1_s2_kp1->SetParLimits(0,hmass_S2_kp_1->GetMaximum() / 3,hmass_S2_kp_1->GetMaximum());
  func1_s2_kp1->SetParLimits(1,0.480,0.505);
  func1_s2_kp1->SetParLimits(2,0.005,0.03);
  func1_s2_kp1->SetParLimits(7,hmass_S2_kp_1->GetMaximum() / 3,hmass_S2_kp_1->GetMaximum());
  func1_s2_kp1->SetParLimits(8,0.480,0.505);
  func1_s2_kp1->SetParLimits(9,0.005,0.05);

  // Strangeness 2 - kaon 1
  hmass_S2_kp_1->Fit("func1_s2_kp1","RB");
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


  // Strangeness 3 - kaon 1



  // Saving the function for part 2
  func1->Write();
  func2->Write();
  func3->Write();
  func4->Write();
  func5->Write();
  func1_s2_kp1->Write();
  func2_s2_kp1->Write();
  func3_s2_kp1->Write();
  func4_s2_kp1->Write();
  func5_s2_kp1->Write();


  fileOutput1.Write();

}
