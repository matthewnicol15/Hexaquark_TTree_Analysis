#include <cstdlib>
#include <iostream>
#include <chrono>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TH2.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TBenchmark.h>
#include "clas12reader.h"
#include <vector>
using namespace clas12;


void Tree_Maker_RGB_Exact(){
  auto start = std::chrono::high_resolution_clock::now();
  gBenchmark->Start("timer");
  int counter=0;

  gROOT->ProcessLine(".L ./Loader.C+"); // Uses Loader.C file, make sure Loader.C is in this file path

  // Filepath to input files, use *.hipo to analyse all hipo files in a directory
  TString inputFile("/volatile/clas12/rg-b/production/recon/spring2019/torus-1/pass1/v0/dst/train/inc/*.hipo");

  // Creating a chain for the data from different files
  TChain fake("hipo");
  fake.Add(inputFile.Data());

  // Shortcut to a list of all the input file names
  auto files=fake.GetListOfFiles();

  // Create root file to save TTree in
  TFile f("/volatile/clas12/matthewn/RGB_Spring2019_Inbending_Pass1_v0_inc_1e_2pos_1neg_Tree_140920_01.root","recreate");
  // Creating TTree object
  TTree RGB_inc_Tree_140920_01("RGB_inc_Tree_140920_01","it's a tree!");

  // Access the PDG database to get information usind PID (e.g db->GetParticle(211)->Mass() for pi^{+} mass)
  auto db=TDatabasePDG::Instance();

  // Creating diagnostic histograms
  auto* h_invariant_lambda=new TH1F("h_invariant_lambda","Invariant mass of lambda",200,-1,3);
  auto* h_missing_mass=new TH1F("h_missing_mass","Missing mass of lambda",200,-1,3);
  auto* h_mass_P=new TH1F("h_mass_P","Missing mass of lambda",200,-1,3);
  auto* h_mass_N=new TH1F("h_mass_N","Missing mass of lambda",200,-1,3);


  // Creating TLorentzVectors of particles for diagnostics
  TLorentzVector Positive_1,Positive_2; // Unidentified
  TLorentzVector e_scattered,Proton,Kaon_Positive,Pion_Minus; // Identified
  Int_t Pos_1, Pos_2;

  // Creating TLorentzVectors of invariant and missing mass for diagnostics
  TLorentzVector Invariant_Lambda,Missing_Mass_Lambda;

  // Information to save to the tree
  // Any information specific to an event
  Int_t eventno; // Records the event number
  Int_t runno; // Records the run number
  Int_t triggerno; // Records the trigger number
  Double_t start_time; // Records the start time for the event

  // Any information specific to an individual particle
  TLorentzVector p4; // TLorentzVector for four vector of each particle (Px,Py,Pz,E)
  TLorentzVector vertex;   // Vertex position and time (Vx,Vy,Vz,Vt)
  Double_t beta; // beta from time of flight (TOF)
  Double_t status; // Gives information on which detectors and how many picked it up
  Double_t energy; // Energy measured by detector
  Double_t charge; // Charge measured by drift chambers
  Double_t PID; // Records PID deterimined by TOF
  Double_t chi2PID; // Chi^2 for PID of TOF
  Double_t time; // Time recorded by by FTCAL,FTOF or CTOF
  Double_t path; // Path measured by FTCAL,FTOF or CTOF
  Double_t vertex_time; // Calculated vertex time from TOF information
  Double_t Region; // Records 0.0 for FT, 1.0 for FD and 2.0 for CD
  Double_t Mass; // Records the calculated mass



  // Vectors of particle measurables for when you have more than one of a type of particle (e.g 2 pi^{-})
  vector<TLorentzVector> v_p4;
  vector<TLorentzVector> v_vertex;
  vector<double> v_beta;
  vector<double> v_status;
  vector<double> v_energy;
  vector<double> v_charge;
  vector<double> v_PID;
  vector<double> v_chi2PID;
  vector<double> v_time;
  vector<double> v_path;
  vector<double> v_region;
  vector<double> v_Mass;


  // Record the number of particles in an event
  Int_t elno; // electrons
  Int_t positronno; // positrons
  Int_t protonno; // protons
  Int_t antiprotonno; // antiprotons
  Int_t neutronno; // neutrons
  Int_t photonno; // photons
  Int_t pipno; // pi^{+}
  Int_t pimno; // pi^{-}
  Int_t pi0no; // pi^{0}
  Int_t kaonpno; // K^{+}
  Int_t kaonmno; // K^{-}
  Int_t positive_charge_tracks; // Count the number of positive charge tracks
  Int_t negative_charge_tracks; // Count the number of negative charge tracks

  // Setting TLorentzVectors for beam and target, for final analysis beam energy
  // will have to be changed depending on which run/runs you are analysing. Here
  // it is just set to 10.6 GeV
  TLorentzVector beam(0,0,10.6,10.6); // Set 4 vector four the beam, all momentum in z-direction
  TLorentzVector target(0,0,0,1.8756); // Set 4 vector for target, stationary so no momentum

  Double_t c = 30; // Speed of light, used to calculate vertex time

  // Assign a branch to each measurable and name it
  RGB_inc_Tree_140920_01.Branch("eventno",&eventno);
  RGB_inc_Tree_140920_01.Branch("runno",&runno,"runno/I");
  RGB_inc_Tree_140920_01.Branch("triggerno",&triggerno,"triggerno/I");
  RGB_inc_Tree_140920_01.Branch("start_time",&start_time);
  RGB_inc_Tree_140920_01.Branch("p4",&v_p4);
  RGB_inc_Tree_140920_01.Branch("vertex",&v_vertex);
  RGB_inc_Tree_140920_01.Branch("beta",&v_beta);
  RGB_inc_Tree_140920_01.Branch("status",&v_status);
  RGB_inc_Tree_140920_01.Branch("energy",&v_energy);
  RGB_inc_Tree_140920_01.Branch("charge",&v_charge);
  RGB_inc_Tree_140920_01.Branch("PID",&v_PID);
  RGB_inc_Tree_140920_01.Branch("chi2PID",&v_chi2PID);
  RGB_inc_Tree_140920_01.Branch("region",&v_region);
  RGB_inc_Tree_140920_01.Branch("Mass",&v_Mass);
  RGB_inc_Tree_140920_01.Branch("time",&v_time);
  RGB_inc_Tree_140920_01.Branch("path",&v_path);
  RGB_inc_Tree_140920_01.Branch("beam",&beam);
  RGB_inc_Tree_140920_01.Branch("target",&target);
  RGB_inc_Tree_140920_01.Branch("elno",&elno);
  RGB_inc_Tree_140920_01.Branch("negative_charge_tracks",&negative_charge_tracks, "neg");
  RGB_inc_Tree_140920_01.Branch("positive_charge_tracks",&positive_charge_tracks);


  // Going over all the input files listed above
  for(Int_t i=0;i<files->GetEntries();i++){

    // Create the CLAS12 event reader
    clas12reader c12(files->At(i)->GetTitle());

    // Prints out the file currently being analysed
    cout<<"file: "<<i<<endl;

    // This loop goes over the events within each file
    while(c12.next()==true){


      counter++;
      // Clear the vectors from the previous event
      v_p4.clear();
      v_vertex.clear();
      v_beta.clear();
      v_status.clear();
      v_energy.clear();
      v_charge.clear();
      v_PID.clear();
      v_chi2PID.clear();
      v_time.clear();
      v_path.clear();
      v_region.clear();
      v_Mass.clear();

      // Setting particles to zero, will be update later once identified
      Proton.SetXYZM(0,0,0,0);
      Kaon_Positive.SetXYZM(0,0,0,0);
      Pion_Minus.SetXYZM(0,0,0,0);

      // Define how to access the information for each event
      eventno = c12.runconfig()->getEvent(); // Getting the event number
      runno = c12.runconfig()->getRun(); // Getting the run number
      triggerno = c12.runconfig()->getTrigger(); // Getting the trigger bit
      start_time = c12.event()->getStartTime(); // Getting start time for each event

      elno = 0;
      positronno = 0;
      protonno = 0;
      antiprotonno = 0;
      neutronno = 0;
      photonno = 0;
      pipno = 0;
      pimno = 0;
      pi0no = 0;
      kaonpno = 0;
      kaonmno = 0;
      positive_charge_tracks = 0;
      negative_charge_tracks = 0;

      // This loop goes over each particle within the current event
      // This loop goes over each particle within the current event
      for(auto& p : c12.getDetParticles()){

        // Define how to access the information for each particle
        p4.SetXYZM(p->par()->getPx(), p->par()->getPy(), p->par()->getPz(),0); // Sets the four vector for each particle, mass will be set once PID determined
        time = p->getTime(); // Gets the TOF
        path = p->getPath(); // Gets the path
        beta = p->par()->getBeta(); // Gets the beta from TOF
        vertex_time = time - path / (beta*c); // Calculated vertex time
        vertex.SetXYZT(p->par()->getVx(), p->par()->getVy(), p->par()->getVz(),vertex_time); // Sets the vertex vector for each particle
        status = p->par()->getStatus(); // Information on which detectors are involved for this particle
        energy =  p->getDetEnergy(); // Energy measured by detector
        charge = p->par()->getCharge(); // Charge of the particle measured (FTB isn't recognised!!!)
        PID = p->par()->getPid();  // PID determined by TOF
        chi2PID = p->par()->getChi2Pid(); // Chi^2 for PID from TOF
        Mass = sqrt((pow(p4.Rho(),2) / (pow(beta,2))) - pow(p4.Rho(),2)); // Calculating Mass from beta and momentum

        // Save the particle information in the vectors by pushing it back
        v_vertex.push_back(vertex);
        v_beta.push_back(beta);
        v_status.push_back(status);
        v_energy.push_back(energy);
        v_charge.push_back(charge);
        v_PID.push_back(PID);
        v_chi2PID.push_back(chi2PID);
        v_time.push_back(time);
        v_path.push_back(path);
        v_p4.push_back(p4); // Recording the 4 vector for each neutron
        v_Mass.push_back(Mass);

        if(PID == 11){
          elno++;
          e_scattered.SetXYZM(p4.Px(),p4.Py(),p4.Pz(),db->GetParticle(11)->Mass());
        }
        else if(charge > 0){
          h_mass_P->Fill(Mass);
          positive_charge_tracks++;
        }
        else if(charge < 0){
          h_mass_N->Fill(Mass);
          negative_charge_tracks++;
        }
        if(p->getRegion()==FT){
          Region=0.0;
        }
        else if(p->getRegion()==FD){
          Region=1.0;
        }
        else if(p->getRegion()==CD){
          Region=2.0;
        }
        v_region.push_back(Region);

      }

      // Here you can apply a basic skim for events you want to save in your tree
      if(positive_charge_tracks == 2 && elno == 1 && negative_charge_tracks == 1){

        Pos_1 = -100;
        Pos_2 = -100;

        // This assigns particles based on charge and mass
        // Loops over the 4 particles (1 e-, 1 other negative and 2 positive
        // assumed to be proton and K+)
        for(Int_t j=0;j<4;j++){
          // Loop over negative particles

          if(v_charge.at(j) < -0.01 && v_PID.at(j) != 11){
            Pion_Minus.SetXYZM(v_p4.at(j).Px(),v_p4.at(j).Py(),v_p4.at(j).Pz(),db->GetParticle(211)->Mass());
          }

          // Loop over positive particles
          else if(v_charge.at(j) > 0.01){
            if(Pos_1<0) Pos_1 = j;
            else Pos_2 = j;
          }
        }
        // Check to see which of the positive particles is heavier and
        // set it to be the proton
        if(v_Mass.at(Pos_1) > v_Mass.at(Pos_2)){
          Proton.SetXYZM(v_p4.at(Pos_1).Px(),v_p4.at(Pos_1).Py(),v_p4.at(Pos_1).Pz(),db->GetParticle(2212)->Mass());
          Kaon_Positive.SetXYZM(v_p4.at(Pos_2).Px(),v_p4.at(Pos_2).Py(),v_p4.at(Pos_2).Pz(),db->GetParticle(321)->Mass());
        }
        else {
          Proton.SetXYZM(v_p4.at(Pos_2).Px(),v_p4.at(Pos_2).Py(),v_p4.at(Pos_2).Pz(),db->GetParticle(2212)->Mass());
          Kaon_Positive.SetXYZM(v_p4.at(Pos_1).Px(),v_p4.at(Pos_1).Py(),v_p4.at(Pos_1).Pz(),db->GetParticle(321)->Mass());
        }


        // Defining the invariant and missing masses
        Invariant_Lambda = Proton + Pion_Minus;
        Missing_Mass_Lambda = beam + target - e_scattered - Kaon_Positive;

        // Filling the diagnostic histograms
        h_invariant_lambda->Fill(Invariant_Lambda.M());
        h_missing_mass->Fill(Missing_Mass_Lambda.M());

        //Fill the TTree
        if(Invariant_Lambda.M() < 1.2 && Missing_Mass_Lambda.M() < 2) RGB_inc_Tree_140920_01.Fill();
      }
    }
  }
  auto finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = finish - start;
  std::cout << "Elapsed time: " << elapsed.count()<< " events = "<<counter<< " s\n";
  RGB_inc_Tree_140920_01.Write();
  f.Write();
  f.Close();
}
