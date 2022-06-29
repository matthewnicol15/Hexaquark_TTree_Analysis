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
#include <TChain.h>
#include <TBenchmark.h>
#include <vector>

// Macro name
void Strangeness_Analysis_Sideband_Kaon_part2(){

  // Read file with information on vectors
  gROOT->ProcessLine(".L ~/work/Macros/Loader.C+");

  // Grabbing part 1 file for the kaon mass fit function
  TFile *fin1 = new TFile("/media/mn688/Elements1/PhD/Analysis_Output/Strangeness_Analysis_RGA_Smeared_Spring2019_Inbending_dst_Tree_Total_221021__10012022_01.root");
  // Getting the kaon mass fit functions from part 1
  // Strangeness 1 - kaon 1
  // Total function
  TF1 *func1 = (TF1*)fin1->Get("func1");
  // Background function
  TF1 *func3 = (TF1*)fin1->Get("func3");


  TF1 *func1_s2_kp1 = (TF1*)fin1->Get("func1_s2_kp1");
  // cout<<func1->GetParameter(0)<<" "<<func1->GetParameter(1)<<" "<<func1->GetParameter(2)<<endl;
  // cout<<func1_s2_kp1->GetParameter(0)<<" "<<func1_s2_kp1->GetParameter(1)<<" "<<func1_s2_kp1->GetParameter(2)<<endl;



  // Read input root file and assign it to 'f'
  TFile *f = new TFile("/media/mn688/Elements1/PhD/Trees/Dibaryon/RGA/RGA_Spring2019_Inbending_at_least_1eFD_1Kp_Tree_Total_201021.root");

  // Read TTree within root file and assign it to 't1'
  TTree *t1 = (TTree*)f->Get("RGA_Spring2019_Inbending_201021");


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
  // TFile fileOutput1("/media/mn688/Elements1/PhD/Analysis_Output/Strangeness_RGA_SPRING_2019_Inbending_eFD_Kp_201021_04_part2_Total_231121_01.root","recreate");


  // Getting particle database to use for masses
  auto db=TDatabasePDG::Instance();

  // Function for strangeness 2 - kaon 2
  TF1 *func1_s2_kp2 = new TF1("func1_s2_kp2","gaus(0) + pol3(3) + gaus(7)",0.36,0.7);
  TF1 *func2_s2_kp2 = new TF1("func2_s2_kp2","gaus(0)",0.36,0.7);
  TF1 *func3_s2_kp2 = new TF1("func3_s2_kp2","pol3(0)",0.36,0.7);
  TF1 *func4_s2_kp2 = new TF1("func4_s2_kp2","gaus(0)",0.36,0.7);
  TF1 *func5_s2_kp2 = new TF1("func5_s2_kp2","gaus(0) + gaus(3)",0.36,0.7);

  // Create histograms here
  // Histograms for events
  auto* hbeam=new TH1D("hbeam","Beam mass; Beam Mass [GeV];Counts",200,0,11);
  auto* hkaon=new TH1D("hkaon","kaon momentum; kaon momentum [GeV];Counts",200,0,10);
  auto* hkaons=new TH1D("hkaons","kaon numbers; Kaons in event;Counts",6,0,6);
  auto* hproton=new TH1D("hproton","proton momentum; proton momentum [GeV];Counts",200,0,10);

  // Histograms for strangeness 1 channel
  auto* hmiss_1_a=new TH1D("hmiss_1_a","MM(e' K^{+});MM(e' K^{+}) [GeV];Counts",200,0,5);
  auto* hmiss_1_a_sig=new TH1D("hmiss_1_a_sig","MM(e' K^{+});MM(e' K^{+}) [GeV];Counts",200,0,5);
  auto* hmiss_1_a_back=new TH1D("hmiss_1_a_back","MM(e' K^{+});MM(e' K^{+}) [GeV];Counts",200,0,5);

  auto* hregion=new TH1D("hregion","Regions;Region;Counts",3,1,4);

  // Histograms for strangeness 2 channel
  auto* hmass_S2_kp_2=new TH1F("hmass_S2_kp_2","K^{+} mass;M(K^{+});Counts",100,0.2,0.8);
  auto* hmass_S2_kp_2_sig=new TH1F("hmass_S2_kp_2_sig","K^{+} mass;M(K^{+});Counts",100,0.2,0.8);
  auto* hmass_S2_kp_2_back=new TH1F("hmass_S2_kp_2_back","K^{+} mass;M(K^{+});Counts",100,0.2,0.8);

  // Histograms for strangeness 3 channel


  // Create vectors of TLorentzVectors to store information of
  // all particles of a given type (important when you have more than 1
  // of any particle in the same event)
  vector<TLorentzVector> v_el;  // e^-
  vector<TLorentzVector> v_pip; // pi^+
  vector<TLorentzVector> v_pim; // pi^-
  vector<TLorentzVector> v_pr; // protons
  vector<TLorentzVector> v_kp; // K^+
  vector<TLorentzVector> v_km; // K^-

  // TLorentzVectors for individual particles
  TLorentzVector el;
  TLorentzVector pip;
  TLorentzVector pim;
  TLorentzVector pr;
  TLorentzVector kp;
  TLorentzVector km;

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


  Double_t c=30;  // Speed of light used for calculating vertex time

  // Reads the total number of entries in the TTree
  Long64_t nentries = t1->GetEntries();
  // You can just run over a set number of events for fast analysis
  // Long64_t nentries = 90000;
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

    // Setting beam energy to 10.2 GeV for RGB Spring 2020
    beam.SetXYZM(0,0,10.2,0);

    // beam = (TLorentzVector)*readbeam;

    hbeam->Fill(beam.Rho());

    ////////////////////////////////////////////////////////////////////////////
    //// Filling histograms for peak and background regions   //////////////////
    ////////////////////////////////////////////////////////////////////////////

    // Here you can apply conditions on the events you want to analyse
    if(v_kp.size() > 0 && v_el.size() == 1 &&v_region_el.at(0) == 1){
      miss1 = beam + (TLorentzVector)*readtarget - v_el.at(0) - v_kp.at(0);

      // Looking at stragneness 1 channel
      if(v_kp.size()==1){

        // Selecting events where kaon is in the FD
        if(v_region_kp.at(0)!=1){
          continue;
        }

        // Cutting on delta beta and momentum of kaons
        if(fabs(v_delta_beta_kp.at(0))<0.02 && (v_kp.at(0).Rho() < 0.55 || v_kp.at(0).Rho() > 0.95)){
          hmiss_1_a->Fill(miss1.M());

          // Grabbing the missing mass of e'K^+ for the signal part of the
          if(v_Mass_kp.at(0) > func1->GetParameter(1) - 2*func1->GetParameter(2) && v_Mass_kp.at(0) < func1->GetParameter(1) + 2*func1->GetParameter(2)){
            hmiss_1_a_sig->Fill(miss1.M());

          }

          // Grabbing the missing mass of e'K^+ for the signal part of the
          else if((v_Mass_kp.at(0) > func1->GetParameter(1) - 7 * func1->GetParameter(2) && v_Mass_kp.at(0) < func1->GetParameter(1) - 5 * func1->GetParameter(2))
          || (v_Mass_kp.at(0) > func1->GetParameter(1) + 5 * func1->GetParameter(2) && v_Mass_kp.at(0) < func1->GetParameter(1) + 7 * func1->GetParameter(2))){

            hmiss_1_a_back->Fill(miss1.M());

          }
        }
      }

      // Looking at stragneness 2 channel
      if(v_kp.size()==2){
        miss_s2 = beam + (TLorentzVector)*readtarget - v_el.at(0) - v_kp.at(0)- v_kp.at(1);

        // Selecting events where kaons are in the FD
        if(v_region_kp.at(0) != 1 || v_region_kp.at(1) != 1) continue;

        // Cutting on delta beta and momentum of kaons
        if(fabs(v_delta_beta_kp.at(0))<0.02 && fabs(v_delta_beta_kp.at(1))<0.02 &&
        (v_kp.at(0).Rho() < 0.55 || v_kp.at(0).Rho() > 0.95) &&
        (v_kp.at(1).Rho() < 0.55 || v_kp.at(1).Rho() > 0.95)){
          hmass_S2_kp_2->Fill(v_Mass_kp.at(1));

          // Grabbing the missing mass of e'K^+ for the signal part of the
          if(v_Mass_kp.at(0) > func1->GetParameter(1) - 2*func1->GetParameter(2) && v_Mass_kp.at(0) < func1->GetParameter(1) + 2*func1->GetParameter(2)){
            hmass_S2_kp_2_sig->Fill(v_Mass_kp.at(1));
          }
          // Grabbing the missing mass of e'K^+ for the signal part of the
          else if((v_Mass_kp.at(0) > func1->GetParameter(1) - 7 * func1->GetParameter(2) && v_Mass_kp.at(0) < func1->GetParameter(1) - 5 * func1->GetParameter(2))
          || (v_Mass_kp.at(0) > func1->GetParameter(1) + 5 * func1->GetParameter(2) && v_Mass_kp.at(0) < func1->GetParameter(1) + 7 * func1->GetParameter(2))){
            hmass_S2_kp_2_back->Fill(v_Mass_kp.at(1));
          }
        }
      }

      // Strangeness 3
      if(v_kp.size() == 3){
        miss_s3 = beam + (TLorentzVector)*readtarget - v_el.at(0) - v_kp.at(0) - v_kp.at(1) - v_kp.at(2);

        // Selecting events where kaons are in the FD
        if(v_region_kp.at(0) != 1 || v_region_kp.at(1) != 1 || v_region_kp.at(2) != 1) continue;

        // Cutting on delta beta and momentum of kaons
        if(fabs(v_delta_beta_kp.at(0))<0.02 && fabs(v_delta_beta_kp.at(1))<0.02 && fabs(v_delta_beta_kp.at(2))<0.02 &&
        (v_kp.at(0).Rho() < 0.55 || v_kp.at(0).Rho() > 0.95) &&
        (v_kp.at(1).Rho() < 0.55 || v_kp.at(1).Rho() > 0.95) &&
        (v_kp.at(2).Rho() < 0.55 || v_kp.at(2).Rho() > 0.95)){
        }
      }
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  //// Perform sideband subtraction          ///////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////

  // Creating variables for multiplication factors
  // Strangeness 1 - kaon 1
  Double_t Peak_Background, Left_Background, Right_Background;
  Double_t Multiplication_factor;

  // Determining the multiplication factor
  // Strangeness 1
  Peak_Background = func1->Integral(hmiss_1_a->FindBin(func1->GetParameter(1) - 2 * func1->GetParameter(2)),
  hmiss_1_a->FindBin(func1->GetParameter(1) + 2 * func1->GetParameter(2)));

  Left_Background = func1->Integral(hmiss_1_a->FindBin(func1->GetParameter(1) - 7 * func1->GetParameter(2)),
  hmiss_1_a->FindBin(func1->GetParameter(1) - 5 * func1->GetParameter(2)));

  Right_Background = func1->Integral(hmiss_1_a->FindBin(func1->GetParameter(1) + 5 * func1->GetParameter(2)),
  hmiss_1_a->FindBin(func1->GetParameter(1) + 7 * func1->GetParameter(2)));
  cout<<Peak_Background<<" "<<Left_Background<<" "<<Right_Background<<endl;

  Multiplication_factor = Peak_Background / (Left_Background + Right_Background);

  // Scaling the background histograms by multiplication factor
  hmiss_1_a_back->Scale(Multiplication_factor);

  // Cloning the signal histograms to use for the sbtracted results
  TH1D *hmiss_1_a_result = (TH1D*)hmiss_1_a_sig->Clone("hmiss_1_a_result");
  TH1D *hmass_S2_kp_2_result = (TH1D*)hmass_S2_kp_2_sig->Clone("hmass_S2_kp_2_result");

  // Subtracting the background histograms from signal
  hmiss_1_a_result ->Add(hmiss_1_a_back,-1);
  hmass_S2_kp_2_result ->Add(hmass_S2_kp_2_back,-1);

  //////////////////////////////////////////////////////////////////////////////
  //// Fitting functions to sideband subtracted result  ////////////////////////
  //////////////////////////////////////////////////////////////////////////////

  // Strangeness 2 - kaon 2
  // Setting parameters before fitting
  func1_s2_kp2->SetParameters(hmass_S2_kp_2_result->GetMaximum() / 2,0.493,0.02);
  func1_s2_kp2->SetParameter(7,hmass_S2_kp_2_result->GetMaximum() / 2);
  func1_s2_kp2->SetParameter(8,0.493);
  func1_s2_kp2->SetParameter(9,0.02);
  // Setting parameter limits before fitting
  func1_s2_kp2->SetParLimits(0,hmass_S2_kp_2_result->GetMaximum() / 3,hmass_S2_kp_2_result->GetMaximum());
  func1_s2_kp2->SetParLimits(1,0.480,0.505);
  func1_s2_kp2->SetParLimits(2,0.005,0.03);
  func1_s2_kp2->SetParLimits(7,hmass_S2_kp_2_result->GetMaximum() / 3,hmass_S2_kp_2_result->GetMaximum());
  func1_s2_kp2->SetParLimits(8,0.480,0.505);
  func1_s2_kp2->SetParLimits(9,0.005,0.05);

  // Fitting function to strangeness 2 - kaon 2
  hmass_S2_kp_2_result->Fit("func1_s2_kp2","RB");

  // Getting individual signal and background functions from the total
  func2_s2_kp2->FixParameter(0, func1_s2_kp2->GetParameter(0));
  func2_s2_kp2->FixParameter(1, func1_s2_kp2->GetParameter(1));
  func2_s2_kp2->FixParameter(2, func1_s2_kp2->GetParameter(2));
  func3_s2_kp2->FixParameter(0, func1_s2_kp2->GetParameter(3));
  func3_s2_kp2->FixParameter(1, func1_s2_kp2->GetParameter(4));
  func3_s2_kp2->FixParameter(2, func1_s2_kp2->GetParameter(5));
  func3_s2_kp2->FixParameter(3, func1_s2_kp2->GetParameter(6));
  func4_s2_kp2->FixParameter(0, func1_s2_kp2->GetParameter(7));
  func4_s2_kp2->FixParameter(1, func1_s2_kp2->GetParameter(8));
  func4_s2_kp2->FixParameter(2, func1_s2_kp2->GetParameter(9));
  func5_s2_kp2->FixParameter(0, func1_s2_kp2->GetParameter(0));
  func5_s2_kp2->FixParameter(1, func1_s2_kp2->GetParameter(1));
  func5_s2_kp2->FixParameter(2, func1_s2_kp2->GetParameter(2));
  func5_s2_kp2->FixParameter(3, func1_s2_kp2->GetParameter(7));
  func5_s2_kp2->FixParameter(4, func1_s2_kp2->GetParameter(8));
  func5_s2_kp2->FixParameter(5, func1_s2_kp2->GetParameter(9));

  // Save functions to extract parameters in part 3
  func1_s2_kp2->Write();
  func2_s2_kp2->Write();
  func3_s2_kp2->Write();
  func4_s2_kp2->Write();
  func5_s2_kp2->Write();

  // Save the output file
  // fileOutput1.Write();

}
