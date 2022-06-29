#include <cstdlib>
#include <iostream>
#include <sstream>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TChain.h>
#include <TBenchmark.h>
#include <vector>

// Macro name
void Tree_Reader_Strangeness_Scaling(){


   //////////////////////////////////////////////////////////////////////////////
   ////Define variables for naming    ///////////////////////////////////////////
   //////////////////////////////////////////////////////////////////////////////

   // Create strings for naming output file
   ostringstream File_Path;
   ostringstream Data;
   ostringstream Quantity;
   ostringstream Date;
   ostringstream Version;
   ostringstream Topology;
   ostringstream Output_File_Name;

   // Setting the strings for output file name
   File_Path<<"/mnt/d/PhD/Analysis_Output/Hexaquark/";
   // File_Path<<"/mnt/c/Users/Nics/Documents/";
   Data<<"S3_Sim_RGB_Spring2019_Inbending_at_least_1e1KpFD_KpChi3_Tree_Total_23062022";
   Quantity<<"Total";
   Topology<<"Scaling";
   Date<<"23062022";
   Version<<"01";


   Output_File_Name<<File_Path.str().c_str()<<Data.str().c_str()<<"_"<<Quantity.str().c_str()<<
   "_"<<Topology.str().c_str()<<"_"<<Date.str().c_str()<<"_"<<Version.str().c_str()<<".root";

   cout<<Output_File_Name.str().c_str()<<endl;

   //////////////////////////////////////////////////////////////////////////////
   //// Setting up input tree and variables    //////////////////////////////////
   //////////////////////////////////////////////////////////////////////////////

   // Read file with information on vectors
   gROOT->ProcessLine(".L /mnt/d/PhD/Macros/Loader.C+");

   // Read input root file and assign it to 'f'
   TFile *f = new TFile("/mnt/d/PhD/Trees/Simulations/Hexaquark/S3_Sim_job4904_Deuteron_Target_210622_01/S3_Sim_job4904_10_6_Beam_Deuteron_Target_230622_01.root");
   // TFile *f = new TFile("/mnt/c/Users/Nics/Documents/RGA_Fall18_Inbending_Scaling_Strangeness_24022022_01.root");
   // Read TTree within root file and assign it to 't1'
   TTree *t1 = (TTree*)f->Get("Tree");


   // Creating components to read from TTree
   // Set any vectors to 0 before reading from the TTree
   // Event information
   TLorentzVector *readbeam=NULL;  // Information on the beam
   TLorentzVector *readtarget=NULL; // Information on the target
   TLorentzVector Target_Sim(0,0,0,1.8756); // Information on the target
   TLorentzVector Beam_Sim;
   Beam_Sim.SetXYZM(0,0,10.6,0.000511);
   // TLorentzVector *el=NULL;   // 4-vectors for the scattered electron
   Double_t start_time; // Event start time
   Int_t readrunno;
   Int_t readeventno;
   Int_t readtriggerno;
   Int_t positive_charge_tracks; // Number of positive charge tracks in event
   Int_t negative_charge_tracks; // Number of negative charge tracks in event
   // Number of given particle or charge track in each event
   Int_t readchargetracks; // Number of positive or negative charge tracks
   Int_t readelno; // e^-
   Int_t region; // which region the particles go in (FT, FD, CD)
   Double_t electronFD; // Number of e in FD
   Double_t Kaonp_FD_Chi3; // Number of K+ in FD with Chi^2 PID < 3

   // Particle information
   vector<TLorentzVector> *v_p4=0;   // 4-vectors for the detected particles
   vector<TLorentzVector> *v_vertex=0;   // Vertex information for particles
   vector<double> *v_path=0;   // Measured path of particles
   vector<double> *v_time=0;   // Measured time of flight of particles
   vector<double> *v_beta=0;   // Beta measured from FTOF
   vector<double> *v_energy=0;   // Energy of particle
   vector<double> *v_status=0;   // Status of particle
   vector<double> *v_charge=0;   // Charge of particle from Drift Chambers
   vector<double> *v_PID=0;   // PID from FTOF
   vector<double> *v_chi2PID=0;   // Chi^2 of the PID
   vector<Int_t> *v_region=0; // region particle goes in

   // Setting the branch addresses to read from
   t1->SetBranchAddress("eventno",&readeventno);
   t1->SetBranchAddress("runno",&readrunno);
   t1->SetBranchAddress("triggerno",&readtriggerno);
   t1->SetBranchAddress("start_time",&start_time);
   t1->SetBranchAddress("p4",&v_p4);
   t1->SetBranchAddress("vertex",&v_vertex);
   t1->SetBranchAddress("beta",&v_beta);
   t1->SetBranchAddress("status",&v_status);
   t1->SetBranchAddress("energy",&v_energy);
   t1->SetBranchAddress("charge",&v_charge);
   t1->SetBranchAddress("PID",&v_PID);
   t1->SetBranchAddress("chi2PID",&v_chi2PID);
   t1->SetBranchAddress("region",&v_region);
   t1->SetBranchAddress("time",&v_time);
   t1->SetBranchAddress("path",&v_path);
   t1->SetBranchAddress("beam",&readbeam);
   t1->SetBranchAddress("target",&readtarget);
   t1->SetBranchAddress("elno",&readelno);
   t1->SetBranchAddress("positive_charge_tracks",&positive_charge_tracks);


   // Path and name for the output file to save
   TFile fileOutput1(Output_File_Name.str().c_str(),"recreate");

   // Getting particle database to use for masses
   auto db=TDatabasePDG::Instance();

   //////////////////////////////////////////////////////////////////////////////
   //// Create histograms here    ///////////////////////////////////////////////
   //////////////////////////////////////////////////////////////////////////////

   // Event Information
   auto* hbeam=new TH1F("hbeam","Beam mass; Beam Mass [GeV];Counts",200,0,11);
   auto* h_S1_Photon_Energy__Miss_Mass__Kaon_Mass = new TH3F("h_S1_Photon_Energy__Miss_Mass__Kaon_Mass","Photon energy, missing mass and kaon mass",120,0,12,300,0,3,500,0.3,0.8);
   auto* h_S2_Photon_Energy__Miss_Mass__Kaon_Mass = new TH3F("h_S2_Photon_Energy__Miss_Mass__Kaon_Mass","Photon energy, missing mass and kaon mass",120,0,12,300,0,3,500,0.3,0.8);
   auto* h_S3_Photon_Energy__Miss_Mass__Kaon_Mass = new TH3F("h_S3_Photon_Energy__Miss_Mass__Kaon_Mass","Photon energy, missing mass and kaon mass",120,0,12,300,2,5,500,0.3,0.8);


   // Missing mass, kaon momentum and kaon mass
   auto* h_S1_Kaon_Momentum__Miss_Mass__Kaon_Mass = new TH3F("h_S1_Kaon_Momentum__Miss_Mass__Kaon_Mass","Kaon momentum, missing mass and kaon mass",300,0,3,300,0,3,500,0.3,0.8);
   auto* h_S2_Kaon_Momentum__Miss_Mass__Kaon_Mass = new TH3F("h_S2_Kaon_Momentum__Miss_Mass__Kaon_Mass","Kaon momentum, missing mass and kaon mass",300,0,3,300,0,3,500,0.3,0.8);
   auto* h_S3_Kaon_Momentum__Miss_Mass__Kaon_Mass = new TH3F("h_S3_Kaon_Momentum__Miss_Mass__Kaon_Mass","Kaon momentum, missing mass and kaon mass",300,0,3,300,2,5,500,0.3,0.8);

   auto* h_S1_Miss_Mass = new TH1F("h_S1_Miss_Mass","Missing Mass",300,0,3);


   // Histograms looking at kaon properties
   auto* h_delta_beta_S1_kp_1 = new TH2F("h_delta_beta_S1_kp_1","#Delta#Beta of K^{+} (1);P [GeV];#Delta#Beta",480,0,12,400,-0.4,0.4);
   auto* h_delta_beta_S2_kp_1 = new TH2F("h_delta_beta_S2_kp_1","#Delta#Beta of K^{+} (1);P [GeV];#Delta#Beta",480,0,12,400,-0.4,0.4);
   auto* h_delta_beta_S2_kp_2 = new TH2F("h_delta_beta_S2_kp_2","#Delta#Beta of K^{+} (2);P [GeV];#Delta#Beta",480,0,12,400,-0.4,0.4);
   auto* h_delta_beta_S3_kp_1 = new TH2F("h_delta_beta_S3_kp_1","#Delta#Beta of K^{+} (1);P [GeV];#Delta#Beta",480,0,12,400,-0.4,0.4);
   auto* h_delta_beta_S3_kp_2 = new TH2F("h_delta_beta_S3_kp_2","#Delta#Beta of K^{+} (2);P [GeV];#Delta#Beta",480,0,12,400,-0.4,0.4);
   auto* h_delta_beta_S3_kp_3 = new TH2F("h_delta_beta_S3_kp_3","#Delta#Beta of K^{+} (3);P [GeV];#Delta#Beta",480,0,12,400,-0.4,0.4);

   // Particle Information
   auto* hmass_pos=new TH1F("hmass_pos","Calculated Mass positives;Mass [GeV];Counts",200,0,2);
   auto* hmass_neg=new TH1F("hmass_neg","Calculated Mass negatives;Mass [GeV];Counts",200,0,2);

   //////////////////////////////////////////////////////////////////////////////
   //// Making variables    ///////////////////////////////////////////////
   //////////////////////////////////////////////////////////////////////////////

   // Create vectors of TLorentzVectors to store information of
   // all particles of a given type (important when you have more than 1
   // of any particle in the same event)
   vector<TLorentzVector> v_el;  // e^-
   vector<TLorentzVector> v_pip; // pi^+
   vector<TLorentzVector> v_pim; // pi^-
   vector<TLorentzVector> v_pr; // protons
   vector<TLorentzVector> v_neutron; // neutrons
   vector<TLorentzVector> v_kp; // K^+
   vector<TLorentzVector> v_km; // K^-

   // TLorentzVectors for individual particles
   TLorentzVector Photon;
   TLorentzVector el;
   TLorentzVector pip;
   TLorentzVector pim;
   TLorentzVector pr;
   TLorentzVector kp;
   TLorentzVector km;
   TLorentzVector neutron;

   // Vertex information for different particle types
   vector<TLorentzVector> v_vertex_el; // e^-
   vector<TLorentzVector> v_vertex_pip; // pi^+
   vector<TLorentzVector> v_vertex_pim; // pi^-
   vector<TLorentzVector> v_vertex_pr; // protons
   vector<TLorentzVector> v_vertex_kp; // K^+
   vector<TLorentzVector> v_vertex_km; // K^-
   vector<TLorentzVector> v_vertex_neutron; // neutron
   TLorentzVector vertex_el;
   TLorentzVector vertex_pip;
   TLorentzVector vertex_pim;
   TLorentzVector vertex_pr;
   TLorentzVector vertex_kp;
   TLorentzVector vertex_km;
   TLorentzVector vertex_neutron;


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
   vector<Double_t> v_Mass_km; // Calculated mass
   Double_t Mass_km;

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

   //////////////////////////////////////////////////////////////////////////////
   // Other variables

   // These are used to define the missing masses later
   TLorentzVector beam;
   TLorentzVector Miss_S1_eKp, Miss_S2_eKpKp, Miss_S3_eKpKpKp;

   Double_t c = 30;  // Speed of light used for calculating vertex time
   Double_t Mass; // Calculated mass

   //////////////////////////////////////////////////////////////////////////////
   //// Looping over events in the tree    //////////////////////////////////////
   //////////////////////////////////////////////////////////////////////////////


   // Reads the total number of entries in the TTree
   Long64_t nentries = t1->GetEntries();
   // You can just run over a set number of events for fast analysis
   // Long64_t nentries = 100;

   cout<<nentries<<endl;

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
      v_Mass_km.clear();

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

      //////////////////////////////////////////////////////////////////////////////
      //// Looping over particles in the current event    //////////////////////////
      //////////////////////////////////////////////////////////////////////////////
      // This reads the number of particles in the current entry/event
      Int_t Nparticles = v_p4->size();

      // This loops over all the particles in the current entry/event
      for(Int_t j=0; j<Nparticles; j++){

         // Calculating the mass for each particle using their beta and momentum
         Mass = sqrt((pow(v_p4->at(j).Rho(),2) / (pow(v_beta->at(j),2))) - pow(v_p4->at(j).Rho(),2));

         // Filling histogram showing calculated masses of particles
         if(v_charge->at(j) > 0) hmass_pos->Fill(Mass);
         else if(v_charge->at(j) < 0) hmass_neg->Fill(Mass);

         // Checking PID and assigning particles
         // e^-
         if(v_PID->at(j)==11){
            // Setting the 4-vector and assigning mass from PDG
            el.SetXYZM(v_p4->at(j).Px(),v_p4->at(j).Py(),v_p4->at(j).Pz(),db->GetParticle(11)->Mass());
            TOF_el = v_time->at(j); // Measured time
            path_el = v_path->at(j); // Measured path
            beta_tof_el = v_beta->at(j); // Measured beta from FTOF
            // Calculating beta from momentum and mass
            beta_calc_el = el.Rho()/(sqrt(pow(el.Rho(),2)+pow(el.M(),2)));
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
            beta_calc_pip = pip.Rho()/(sqrt(pow(pip.Rho(),2)+pow(pip.M(),2)));
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
            beta_calc_pim = pim.Rho()/(sqrt(pow(pim.Rho(),2)+pow(pim.M(),2)));
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
            beta_calc_pr = pr.Rho()/(sqrt(pow(pr.Rho(),2)+pow(pr.M(),2)));
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
            beta_calc_kp = kp.Rho()/(sqrt(pow(kp.Rho(),2)+pow(kp.M(),2)));
            // Difference between calculated and measured beta
            delta_beta_kp = beta_calc_kp-beta_tof_kp;
            vertex_time_kp = TOF_kp - path_kp / (beta_tof_kp*c); // Calculate vertex time
            // Setting the vertex information now vertex time has been calculated
            vertex_kp.SetXYZT(v_vertex->at(j).X(), v_vertex->at(j).Y(), v_vertex->at(j).Z(), vertex_time_kp);
            region_kp = v_region->at(j);
            Mass_kp = sqrt((pow(v_p4->at(j).Rho(),2) / (pow(beta_tof_kp,2))) - pow(v_p4->at(j).Rho(),2));

            // Only keep kaons in FD
            if(region_kp != 1) continue;

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
            km.SetXYZM(v_p4->at(j).Px(),v_p4->at(j).Py(),v_p4->at(j).Pz(),db->GetParticle(321)->Mass());
            TOF_km = v_time->at(j); // Measured time
            path_km = v_path->at(j); // Measured path
            beta_tof_km = v_beta->at(j); // Measured beta from FTOF
            // Calculating beta from momentum and mass
            beta_calc_km = km.Rho()/(sqrt(pow(km.Rho(),2)+pow(km.M(),2)));
            // Difference between calculated and measured beta
            delta_beta_km = beta_calc_km-beta_tof_km;
            vertex_time_km = TOF_km - path_km / (beta_tof_km*c); // Calculate vertex time
            // Setting the vertex information now vertex time has been calculated
            vertex_km.SetXYZT(v_vertex->at(j).X(), v_vertex->at(j).Y(), v_vertex->at(j).Z(), vertex_time_km);
            region_km = v_region->at(j);
            Mass_km = sqrt((pow(v_p4->at(j).Rho(),2) / (pow(beta_tof_km,2))) - pow(v_p4->at(j).Rho(),2));

            // Only keep kaons in FD
            if(region_km != 1) continue;

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
            v_Mass_km.push_back(Mass_km);
         }



         // neutron
         else if(v_PID->at(j)==2112){
            // Setting the 4-vector and assigning mass from PDG
            neutron.SetXYZM(v_p4->at(j).Px(),v_p4->at(j).Py(),v_p4->at(j).Pz(),db->GetParticle(2112)->Mass());
            TOF_neutron = v_time->at(j); // Measured time
            path_neutron = v_path->at(j); // Measured path
            beta_tof_neutron = v_beta->at(j); // Measured beta from FTOF
            // Calculating beta from momentum and mass
            beta_calc_neutron = neutron.Rho()/(sqrt(pow(neutron.Rho(),2)+pow(neutron.M(),2)));
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
      // if(readrunno < 11394)beam.SetXYZM(0,0,beam_energy_low,0);
      // else beam.SetXYZM(0,0,beam_energy_high,06);
      // if(readrunno < 6400)   beam.SetXYZM(0,0,beam_energy,0);
      // else beam.SetXYZM(0,0,beam_energy_low,0);
      // beam.SetXYZM(0,0,10.6,0);
      hbeam->Fill(beam.Rho());

      //////////////////////////////////////////////////////////////////////////////
      //// Selecting Events for specific topolgies    //////////////////////////////
      //////////////////////////////////////////////////////////////////////////////


      // Here you can apply conditions on the events you want to analyse
      // Requiring exactly in electron and at least 1 kaon
      if(v_el.size() != 1 || v_kp.size() < 1) continue;

      // Reconstructing photon
      Photon = beam - v_el.at(0);

      // Strangness 1
      if(v_kp.size() == 1){


         // Missing mass
         Miss_S1_eKp = beam + (TLorentzVector)*readtarget - v_el.at(0) - v_kp.at(0);

         // K+ properties
         h_delta_beta_S1_kp_1->Fill(v_kp.at(0).Rho(),v_delta_beta_kp.at(0));
         h_S1_Kaon_Momentum__Miss_Mass__Kaon_Mass->Fill(v_kp.at(0).Rho(), Miss_S1_eKp.M(),v_Mass_kp.at(0));

         h_S1_Miss_Mass->Fill(Miss_S1_eKp.M());

         // Look at photon energy, missing mass and kaon mass
         h_S1_Photon_Energy__Miss_Mass__Kaon_Mass->Fill(Photon.E(),Miss_S1_eKp.M(),v_Mass_kp.at(0));
      }

      // Strangness 2
      if(v_kp.size() == 2){

         // Missing mass
         Miss_S2_eKpKp = beam + (TLorentzVector)*readtarget - v_el.at(0) - v_kp.at(0) - v_kp.at(1);

         // K+ properties
         h_delta_beta_S2_kp_1->Fill(v_kp.at(0).Rho(),v_delta_beta_kp.at(0));
         h_delta_beta_S2_kp_2->Fill(v_kp.at(1).Rho(),v_delta_beta_kp.at(1));
         h_S2_Kaon_Momentum__Miss_Mass__Kaon_Mass->Fill(v_kp.at(0).Rho(), Miss_S2_eKpKp.M(),v_Mass_kp.at(0));

         // Delta beta cut
         // if(v_kp.at(0).Rho() > 1.33 && v_kp.at(1).Rho() > 0.9 && fabs(v_delta_beta_kp.at(0)) < 0.02 && fabs(v_delta_beta_kp.at(1)) < 0.02){

         // Look at photon energy, missing mass and kaon mass
         h_S2_Photon_Energy__Miss_Mass__Kaon_Mass->Fill(Photon.E(),Miss_S2_eKpKp.M(),v_Mass_kp.at(0));
         // }
      }

      // Strangness 3
      if(v_kp.size() == 3){

         // Missing mass
         // Miss_S3_eKpKpKp = beam + (TLorentzVector)*readtarget - v_el.at(0) - v_kp.at(0) - v_kp.at(1) - v_kp.at(2);
         Miss_S3_eKpKpKp = Beam_Sim + (TLorentzVector)*readtarget - v_el.at(0) - v_kp.at(0) - v_kp.at(1) - v_kp.at(2);

         // K+ properties
         h_delta_beta_S3_kp_1->Fill(v_kp.at(0).Rho(),v_delta_beta_kp.at(0));
         h_delta_beta_S3_kp_2->Fill(v_kp.at(1).Rho(),v_delta_beta_kp.at(1));
         h_delta_beta_S3_kp_3->Fill(v_kp.at(2).Rho(),v_delta_beta_kp.at(2));
         h_S3_Kaon_Momentum__Miss_Mass__Kaon_Mass->Fill(v_kp.at(0).Rho(), Miss_S3_eKpKpKp.M(),v_Mass_kp.at(0));

         // Delta beta cut
         // if(fabs(v_delta_beta_kp.at(0)) < 0.02 && fabs(v_delta_beta_kp.at(1)) < 0.02 &&
         // fabs(v_delta_beta_kp.at(2)) < 0.02){

         // Look at photon energy, missing mass and kaon mass
         h_S3_Photon_Energy__Miss_Mass__Kaon_Mass->Fill(Photon.E(),Miss_S3_eKpKpKp.M(),v_Mass_kp.at(0));
         // }
      }
   }

   fileOutput1.Write();
}
