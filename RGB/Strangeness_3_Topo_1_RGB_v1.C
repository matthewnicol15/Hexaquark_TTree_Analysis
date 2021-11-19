#include <cstdlib>
#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TH2.h>
#include <TChain.h>
#include <TBenchmark.h>
#include <vector>

// Macro name
void Strangeness_3_Topo_1a_RGB_v1(){

  // Read file with information on vectors
  gROOT->ProcessLine(".L /mnt/f/PhD/Macros/Loader.C+");

  // Read input root file and assign it to 'f'
  TFile *f = new TFile("/mnt/f/PhD/Trees/RGB_Inc_Spring2019_Inbending_Pass1_v0_e_Kp_FD_Tree_220920_01.root");
  // Read TTree within root file and assign it to 't1'
  TTree *t1 = (TTree*)f->Get("RGB_Inc_Tree_220920_01");


  // Creating components to read from TTree
  // Set any vectors to 0 before reading from the TTree
  // Event information
  TLorentzVector *readbeam=NULL;  // Information on the beam
  TLorentzVector *readtarget=NULL; // Information on the target
  Double_t start_time; // Event start time
  Int_t readrunno;
  Int_t readeventno;
  Int_t readtriggerno;
  // Number of given particle or charge track in each event
  Int_t readchargetracks; // Number of positive or negative charge tracks
  Int_t readelno; // e^-
  Int_t pimno,protonno,kaonpno; // Counts number of particles in an event
  Int_t readothertracks; // Number of particles excluding p, e^-, pions or kaons
  Int_t region; // which region the particles go in (FT, FD, CD)


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
  TFile fileOutput1("/mnt/f/PhD/Analysis_Output/RGB/Inclusive/Inbending/Strangeness_3/PID/Strangeness_3_Topo_1a_RGB_Inc_Pass1_Inbending_e_Kp_FD_011020_01.root","recreate");


  // Getting particle database to use for masses
  auto db=TDatabasePDG::Instance();

  // Create histograms here
  auto* hbeam=new TH1F("hbeam","Beam mass; Beam Mass [GeV];Counts",200,0,11);
  auto* hkaon=new TH1F("hkaon","kaon momentum; kaon momentum [GeV];Counts",200,0,10);
  auto* hproton=new TH1F("hproton","proton momentum; proton momentum [GeV];Counts",200,0,10);
  auto* hpimno=new TH1F("hpimno","number of #pi^{-}; no. of #pi^{-};Counts",10,0,10);
  auto* hpimno_1=new TH1F("hpimno_1","number of #pi^{-}; no. of #pi^{-};Counts",10,0,10);
  auto* hprotonno=new TH1F("hprotonno","number of p; no. of p;Counts",10,0,10);
  auto* hregion_kp1=new TH1F("hregion_kp1","Regions hit by K^{+} (1);Regions;Counts",3,0,3);
  auto* hregion_kp2=new TH1F("hregion_kp2","Regions hit by K^{+} (1);Regions;Counts",3,0,3);
  auto* hregion_kp3=new TH1F("hregion_kp3","Regions hit by K^{+} (1);Regions;Counts",3,0,3);
  auto* hmass_kp1=new TH1F("hmass_kp1","Mass of K^{+} (1);Mass [GeV];Counts",400,0,3);
  auto* hmass_kp2=new TH1F("hmass_kp2","Mass of K^{+} (1);Mass [GeV];Counts",400,0,3);
  auto* hmass_kp3=new TH1F("hmass_kp3","Mass of K^{+} (1);Mass [GeV];Counts",400,0,3);
  auto* hangular_distribution_kp1=new TH2D("hangular_distribution_kp1","Theta vs P for K^{+} (1);P [GeV];#theta [rad]",200,0,11,200,0,200);
  auto* hangular_distribution_kp2=new TH2D("hangular_distribution_kp2","Theta vs P for K^{+} (1);P [GeV];#theta [rad]",200,0,11,200,0,200);
  auto* hangular_distribution_kp3=new TH2D("hangular_distribution_kp3","Theta vs P for K^{+} (1);P [GeV];#theta [rad]",200,0,11,200,0,200);

  auto* hinv_lambda=new TH1F("hinv_lambda","Invariant mass of p #pi^{-};M(p #pi^{-}) [GeV];Counts",400,0,2);
  auto* hmiss_5_1=new TH1F("hmiss_5_1","MM(e' K^{+} K^{+} K^{+});MM(e' K^{+}) [GeV];Counts",200,-5,5);
  auto* hmiss_5_2=new TH1F("hmiss_5_2","MM(e' K^{+} K^{+} K^{+});MM(e' K^{+}) [GeV];Counts",200,-5,5);
  auto* hmiss_5_3=new TH1F("hmiss_5_3","MM(e' K^{+} K^{+} K^{+});MM(e' K^{+}) [GeV];Counts",200,-5,5);
  auto* hmiss_5_4=new TH1F("hmiss_5_4","MM(e' K^{+} K^{+} K^{+});MM(e' K^{+}) [GeV];Counts",200,-5,5);
  auto* hmiss_5_5=new TH1F("hmiss_5_5","MM(e' K^{+} K^{+} K^{+});MM(e' K^{+}) [GeV];Counts",200,-5,5);
  auto* hmiss_5_6=new TH1F("hmiss_5_6","MM(e' K^{+} K^{+} K^{+});MM(e' K^{+}) [GeV];Counts",200,-5,5);
  auto* hmiss_5_7=new TH1F("hmiss_5_7","MM(e' K^{+} K^{+} K^{+});MM(e' K^{+}) [GeV];Counts",200,-5,5);
  auto* hmiss_5_8=new TH1F("hmiss_5_8","MM(e' K^{+} K^{+} K^{+});MM(e' K^{+}) [GeV];Counts",200,-5,5);
  auto* hmiss_5_9=new TH1F("hmiss_5_9","MM(e' K^{+} K^{+} K^{+});MM(e' K^{+}) [GeV];Counts",200,-5,5);
  auto* hmiss_5_10=new TH1F("hmiss_5_10","MM(e' K^{+} K^{+} K^{+});MM(e' K^{+}) [GeV];Counts",200,-5,5);
  auto* hmass=new TH1F("hmass","Calculated Mass;Mass [GeV];Counts",200,0,2);
  auto* hregion=new TH1F("hregion","Regions hit;Region;Counts",5,0,4);


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
  TLorentzVector miss5;
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
  vector<Double_t> v_Mass_kp; // Calculated mass
  Double_t Mass_kp;

  // K^+
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
  Double_t Mass; // Calculated mass

  // Reads the total number of entries in the TTree
  // Long64_t nentries = t1->GetEntries();
  // You can just run over a set number of events for fast analysis
  Long64_t nentries = 1000000;

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

    pimno = 0;
    protonno = 0;
    kaonpno = 0;

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


    // This reads the number of particles in the current entry/event
    Int_t Nparticles = v_p4->size();

    // This loops over all the particles in the current entry/event
    for(Int_t j=0; j<Nparticles; j++){

      // Calculating the mass for each particle using their beta and momentum
      Mass = sqrt((pow(v_p4->at(j).Rho(),2) / (pow(beta_tof_kp,2))) - pow(v_p4->at(j).Rho(),2));

      // Filling histogram to show the different regions hit
      hregion->Fill(v_region->at(j));
      // Can put a selection on which regions particles are hitting

      // Filling histogram showing calculated masses of particles
      hmass->Fill(Mass);

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
        Mass_kp = sqrt((pow(v_p4->at(j).Rho(),2) / (pow(beta_tof_kp,2))) - pow(v_p4->at(j).Rho(),2));

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
        v_Mass_kp.push_back(Mass_kp);
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
    beam = (TLorentzVector)*readbeam;
    if(readrunno < 6400) beam.SetXYZM(0,0,10.6,0);
    else beam.SetXYZM(0,0,10.2,0);
    hbeam->Fill(beam.Rho());

    pimno = v_pim.size();
    protonno = v_pr.size();
    kaonpno = v_kp.size();

    // Here you can apply conditions on the events you want to analyse
    if(v_kp.size()==3 && v_el.size()==1){
      // Missing mass of 3 K+, looking for the dsss
      miss5 = beam + (TLorentzVector)*readtarget - v_el.at(0) - v_kp.at(0) - v_kp.at(1) - v_kp.at(2);
      hmass_kp1->Fill(v_Mass_kp.at(0));
      hmass_kp2->Fill(v_Mass_kp.at(1));
      hmass_kp3->Fill(v_Mass_kp.at(2));
      hangular_distribution_kp1->Fill(v_kp.at(0).Theta()*TMath::RadToDeg(),v_kp.at(0).Rho());
      hangular_distribution_kp2->Fill(v_kp.at(1).Theta()*TMath::RadToDeg(),v_kp.at(0).Rho());
      hangular_distribution_kp3->Fill(v_kp.at(2).Theta()*TMath::RadToDeg(),v_kp.at(0).Rho());
      hmiss_5_1->Fill(miss5.M());


      // Selecting events where the K- is also detected
      if(v_km.size() == 1){
        hmiss_5_2->Fill(miss5.M());

        // Selecting events where the proton is also detected
        if(v_pr.size()==1){
          hmiss_5_3->Fill(miss5.M());


          // Selecting events where at least 1 pion is also detected
          if(v_pim.size() == 2){
            hmiss_5_4->Fill(miss5.M());

            // Looking at effect of having the protons kaons in FD
            if(v_region_kp.at(0) > 0.8 && v_region_kp.at(0) < 1.2 &&
            v_region_kp.at(1) > 0.8 && v_region_kp.at(1) < 1.2 &&
            v_region_kp.at(2) > 0.8 && v_region_kp.at(2) < 1.2 &&
            v_region_pr.at(0)> 0.8 && v_region_pr.at(0) < 1.2){

              hmiss_5_5->Fill(miss5.M());
            }


            // Looking at all combinations of protons and pions to identify lambdas
            // Looping over all protons in event
            for(Int_t k=0;k<protonno;k++){
              // Looping over all pi^{-}s in event
              for(Int_t l=0;l<pimno;l++){

                lambda = v_pr.at(k) + v_pim.at(l);
                hinv_lambda->Fill(lambda.M());
                if(lambda.M() < 1.2)hmiss_5_6->Fill(miss5.M());
              }
            }
          }
        }
      }
    }
  }
  fileOutput1.Write();
}
