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
void Tree_Reader_S3_At_Least_1p2pim(){

  //////////////////////////////////////////////////////////////////////////////
  ////Define variables for naming    ///////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////

  // Create strings for naming output file
  ostringstream File_Path;
  ostringstream Data;
  ostringstream Quantity;
  ostringstream Date;
  ostringstream Version;
  ostringstream Output_File_Name;

  // Setting the strings for output file name
  File_Path<<"/media/mn688/Elements1/PhD/Analysis_Output/Hexaquark/";
  Data<<"RGB_Spring2020_Inbending_S3_eFD_At_Least_1p2pim_Tree_030821_01";
  Quantity<<"Total";
  Date<<"03122021";
  Version<<"01";

  Output_File_Name<<File_Path.str().c_str()<<Data.str().c_str()<<"_"<<Quantity.str().c_str()<<"_"<<Date.str().c_str()<<"_"<<Version.str().c_str();

  //////////////////////////////////////////////////////////////////////////////
  //// Setting up input tree and variables    //////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////

  // Read file with information on vectors
  gROOT->ProcessLine(".L /media/mn688/Elements1/PhD/Macros/Loader.C+");

  // Read input root file and assign it to 'f'
  TFile *f = new TFile("/media/mn688/Elements1/PhD/Trees/Dibaryon/RGB/Strangeness_3/RGB_Spring2020_Inbending_1eFD_at_least_1p2pim_Tree_02122021_01.root");
  // Read TTree within root file and assign it to 't1'
  TTree *t1 = (TTree*)f->Get("RGB_Spring2020_Inbending_02122021");


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
  TFile fileOutput1(Output_File_Name.str().c_str(),"recreate");


  // Getting particle database to use for masses
  auto db=TDatabasePDG::Instance();

  //////////////////////////////////////////////////////////////////////////////
  //// Create histograms here    ///////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////

  // Event Information
  auto* hbeam=new TH1F("hbeam","Beam mass; Beam Mass [GeV];Counts",200,0,11);
  auto* hkaonpno=new TH1F("hkaonpno", "Number of K^{+};# of K^{+}; Counts",10,0,10);
  auto* hregion=new TH1F("hregion","Regions hit;Region;Counts",5,0,4);

  // Particle Information
  auto* hmass=new TH1F("hmass","Calculated Mass;Mass [GeV];Counts",200,0,2);

  // Analysis Plots
  auto* hlambdas=new TH2D("hlambdas", "Invariant mass of #Lambda (2) against #Lambda (1);M(p #pi^{-}) [GeV];M(p #pi^{-}) [GeV]",600,1,3,600,1,3);
  auto* hinv_lambda=new TH1F("hinv_lambda","Invariant mass of p #pi^{-};M(p #pi^{-}) [GeV];Counts",400,1,3);
  auto* hinv_lambda_1=new TH1F("hinv_lambda_1","Invariant mass of p #pi^{-};M(p #pi^{-}) [GeV];Counts",400,1,3);
  auto* hinv_lambda_2=new TH1F("hinv_lambda_2","Invariant mass of p #pi^{-};M(p #pi^{-}) [GeV];Counts",400,1,3);
  auto* hinv_lambda_3=new TH1F("hinv_lambda_3","Invariant mass of p #pi^{-};M(p #pi^{-}) [GeV];Counts",400,1,3);
  auto* hinv_cascade=new TH1F("hinv_cascade","Invariant mass of p #pi^{-} #pi^{-};M(p #pi^{-} #pi^{-}) [GeV];Counts",600,1,3);


  // Create vectors of TLorentzVectors to store information of
  // all particles of a given type (important when you have more than 1
  // of any particle in the same event)
  vector<TLorentzVector> v_el;  // e^-
  vector<TLorentzVector> v_pip; // pi^+
  vector<TLorentzVector> v_pim; // pi^-
  vector<TLorentzVector> v_pr; // protons
  vector<TLorentzVector> v_neutron; // neutrons
  vector<TLorentzVector> v_kp; // K^+

  // TLorentzVectors for individual particles
  TLorentzVector el;
  TLorentzVector pip;
  TLorentzVector pim;
  TLorentzVector pr;
  TLorentzVector kp;
  TLorentzVector neutron;
  TLorentzVector othertracks;

  // Vertex information for different particle types
  vector<TLorentzVector> v_vertex_el; // e^-
  vector<TLorentzVector> v_vertex_pip; // pi^+
  vector<TLorentzVector> v_vertex_pim; // pi^-
  vector<TLorentzVector> v_vertex_pr; // protons
  vector<TLorentzVector> v_vertex_kp; // K^+
  vector<TLorentzVector> v_vertex_neutron; // neutron
  TLorentzVector vertex_el;
  TLorentzVector vertex_pip;
  TLorentzVector vertex_pim;
  TLorentzVector vertex_pr;
  TLorentzVector vertex_kp;
  TLorentzVector vertex_neutron;

  // These are used to define the missing masses later
  TLorentzVector beam;
  TLorentzVector proton_pion_1, proton_pion_2, proton_pion_3;
  TLorentzVector cascade_12, cascade_13, cascade_23;


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

  // neutron
  vector<Double_t> v_beta_tof_neutron;  // Beta measured
  Double_t beta_tof_neutron;
  vector<Double_t> v_P_neutron;  // Momentum measured
  Double_t P_neutron;
  vector<Double_t>v_path_neutron; // Path measured
  Double_t path_neutron;
  vector<Double_t>v_TOF_neutron;  // TOF measured
  Double_t TOF_neutron;
  vector<Double_t> v_beta_calc_neutron;  // Beta calculated
  Double_t beta_calc_neutron;
  vector<Double_t> v_delta_beta_neutron;  // Beta calculated - beta measured
  Double_t delta_beta_neutron;
  vector<Double_t> v_vertex_time_neutron;  // Vertex time calculated
  Double_t vertex_time_neutron;
  vector<Double_t> v_region_neutron; // region hit
  Double_t region_neutron;


  Double_t c=30;  // Speed of light used for calculating vertex time
  Double_t Mass; // Calculated mass

  // Reads the total number of entries in the TTree
  Long64_t nentries = t1->GetEntries();
  // You can just run over a set number of events for fast analysis
  // Long64_t nentries = 1000000;

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

    // neutrons
    v_neutron.clear();
    v_beta_tof_neutron.clear();
    v_P_neutron.clear();
    v_path_neutron.clear();
    v_TOF_neutron.clear();
    v_beta_calc_neutron.clear();
    v_delta_beta_neutron.clear();
    v_vertex_neutron.clear();
    v_vertex_time_neutron.clear();
    v_vertex_neutron.clear();
    v_region_neutron.clear();


    // This reads the number of particles in the current entry/event
    Int_t Nparticles = v_p4->size();

    // This loops over all the particles in the current entry/event
    for(Int_t j=0; j<Nparticles; j++){

      // Calculating the mass for each particle using their beta and momentum
      Mass = sqrt((pow(v_p4->at(j).Rho(),2) / (pow(v_beta.at(j),2))) - pow(v_p4->at(j).Rho(),2));

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
        v_Mass_kp.push_back(Mass);
      }


      // neutron
      else if(v_PID->at(j)==2112){
        // Setting the 4-vector and assigning mass from PDG
        neutron.SetXYZM(v_p4->at(j).Px(),v_p4->at(j).Py(),v_p4->at(j).Pz(),db->GetParticle(2112)->Mass());
        TOF_neutron = v_time->at(j); // Measured time
        path_neutron = v_path->at(j); // Measured path
        beta_tof_neutron = v_beta->at(j); // Measured beta from FTOF
        // Calculating beta from momentum and mass
        beta_calc_neutron = neutron.Rho()/(sqrt((pow(neutron.Rho(),2))+(pow(neutron.M(),2))));
        // Difference between calculated and measured beta
        delta_beta_neutron = beta_calc_neutron-beta_tof_neutron;
        vertex_time_neutron = TOF_neutron - path_neutron / (beta_tof_neutron*c); // Calculate vertex time
        // Setting the vertex information now vertex time has been calculated
        vertex_neutron.SetXYZT(v_vertex->at(j).X(), v_vertex->at(j).Y(), v_vertex->at(j).Z(), vertex_time_neutron);
        region_neutron = v_region->at(j);

        // Pushing back all that iformation into the vectors
        // Again this is done so you can store information on multiple particles
        // of the same type in one place
        v_neutron.push_back(neutron);
        v_beta_tof_neutron.push_back(beta_tof_neutron);
        v_P_neutron.push_back(P_neutron);
        v_path_neutron.push_back(path_neutron);
        v_TOF_neutron.push_back(TOF_neutron);
        v_beta_calc_neutron.push_back(beta_calc_neutron);
        v_delta_beta_neutron.push_back(delta_beta_neutron);
        v_vertex_time_neutron.push_back(vertex_time_neutron);
        v_vertex_neutron.push_back(vertex_neutron);
        v_region_neutron.push_back(region_neutron);
      }
    }

    beam = (TLorentzVector)*readbeam;
    // Setting the beam momentum based on run
    beam.SetXYZM(0,0,10.4,0);
    // if(readrunno < 6400) beam.SetXYZM(0,0,10.6,0);
    // else beam.SetXYZM(0,0,10.2,0);
    hbeam->Fill(beam.Rho());

    //////////////////////////////////////////////////////////////////////////////
    //// Selecting Events for specific topolgies    //////////////////////////////
    //////////////////////////////////////////////////////////////////////////////


    // Here you can apply conditions on the events you want to analyse
    // Requiring exact topology 3
    if(v_el.size() != 1 || v_pr.size() < 1 || v_pim.size() < 3 || v_neutron.size() < 1) continue;
    {
      hkaonpno->Fill(v_kp.size());

      //////////////////////////////////////////////////////////////////////////////
      //// Calculating invariant mass combinations    //////////////////////////////
      //////////////////////////////////////////////////////////////////////////////

      // Lambdas or sigmas
      proton_pion_1 = v_pr.at(0) + v_pim.at(0); // First pion
      proton_pion_2 = v_pr.at(0) + v_pim.at(1); // Second pion
      proton_pion_3 = v_pr.at(0) + v_pim.at(2); // Third pion

      // Cascades
      cascade_12 = v_pr.at(0) + v_pim.at(0) + v_pim.at(1); // First and second pion
      cascade_13 = v_pr.at(0) + v_pim.at(0) + v_pim.at(2); // First and third pion
      cascade_23 = v_pr.at(0) + v_pim.at(1) + v_pim.at(2); // Second and third pion

      //////////////////////////////////////////////////////////////////////////////
      //// Checking possible invariant mass combinations    ////////////////////////
      //////////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////////
      // Check if pion 1 is from lambda
      if(proton_pion_1 < 1.18){

        // Check if pion 1 and 2 are from cascade
        if(cascade_12 < 1.4){

          // Check if pion 3 is from sigma
          if(proton_pion_3 > 1.18 && proton_pion_3 < 1.27){


          }

        }

        // Check if pion 1 and 3 are from cascade
        else if(cascade_13 < 1.4){

          // Check if pion 2 is from sigma
          if(proton_pion_2 > 1.18 && proton_pion_2 < 1.27){


          }
        }
      }

      //////////////////////////////////////////////////////////////////////////////
      // Check if pion 2 is from lambda
      if(proton_pion_2 < 1.18){

        // Check if pion 1 and 2 are from cascade
        if(cascade_12 < 1.4){

          // Check if pion 3 is from sigma
          if(proton_pion_3 > 1.18 && proton_pion_3 < 1.27){


          }

        }

        // Check if pion 2 and 3 are from cascade
        else if(cascade_23 < 1.4){

          // Check if pion 1 is from sigma
          if(proton_pion_1 > 1.18 && proton_pion_1 < 1.27){


          }
        }
      }

      //////////////////////////////////////////////////////////////////////////////
      // Check if pion 3 is from lambda
      if(proton_pion_3 < 1.18){

        // Check if pion 1 and 3 are from cascade
        if(cascade_13 < 1.4){

          // Check if pion 2 is from sigma
          if(proton_pion_2 > 1.18 && proton_pion_2 < 1.27){


          }

        }

        // Check if pion 2 and 3 are from cascade
        else if(cascade_23 < 1.4){

          // Check if pion 1 is from sigma
          if(proton_pion_1 > 1.18 && proton_pion_1 < 1.27){


          }
        }
      }


    }

  }
  fileOutput1.Write();
}
