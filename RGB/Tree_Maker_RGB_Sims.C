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


void Tree_Maker_RGB_Sims(){
  auto start = std::chrono::high_resolution_clock::now();
  gBenchmark->Start("timer");
  int counter=0;
  gROOT->ProcessLine(".L /home/matthewn/Documents/Macros/Loader.C+"); // Uses Loader.C file, make sure Loader.C is in this file path

  // Filepath to input files, use *.hipo to analyse all hipo files in a directory
  TString inputFile("/volatile/clas12/matthewn/Simulations/Dibaryon/Cascade_Sigma_10M_RGB_Inben_240221.hipo");
  // TString inputFile("/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass1/v0/dst/train/skim4/skim4_005032.hipo");

  // Creating a chain for the data from different files
  TChain fake("hipo");
  fake.Add(inputFile.Data());

  // Shortcut to a list of all the input file names
  auto files=fake.GetListOfFiles();

  // Create root file to save TTree in
  TFile f("/volatile/clas12/matthewn/Simulations/Dibaryon/Cascade_Sigma_Sim_10M_240221_RGB_Fall2018_Inbending_S3_Tree_240221.root","recreate");
  // Creating TTree object
  TTree Cascade_Sigma_Sim_10M_240221_RGB_Fall2018_Inbending("Cascade_Sigma_Sim_10M_240221_RGB_Fall2018_Inbending","it's a tree!");


  // Access the PDG database to get information usind PID (e.g db->GetParticle(211)->Mass() for pi^{+} mass)
  auto db=TDatabasePDG::Instance();
  // Information to save to the tree
  // Any information specific to an event
  Int_t eventno; // Records the event number
  Int_t runno; // Records the run number
  Int_t triggerno; // Records the trigger number
  Double_t start_time; // Records the start time for the event
  Int_t Nparticles; // Records the number of particles in an event (charged only)
  Int_t kaonpno = 0; // Records the number of K^{+} in each event
  Int_t elno = 0; // Records the number of e^{-} in each event
  Int_t GoodEvents = 0; // Records the number of events that pass the cuts

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
  Int_t Region; // Records 0 for FT, 1 for FD and 2 for CD
  Double_t mass; // calculated mass of positive particles
  int Pos_position;

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
  vector<int> v_region;
  vector<Double_t> v_mass;
  vector<int> v_Pos_position;

  // Record the number of particles in an event
  Int_t positive,negative; // Positive and negative tracks in an event

  // Setting TLorentzVectors for beam and target, for final analysis beam energy
  // will have to be changed depending on which run/runs you are analysing. Here
  // it is just set to 10.6 GeV
  TLorentzVector beam(0,0,10.6,10.6); // Set 4 vector four the beam, all momentum in z-direction
  TLorentzVector target(0,0,0,1.8756); // Set 4 vector for target, stationary so no momentum

  Double_t c=30; // Speed of light, used to calculate vertex time

  // Assign a branch to each measurable and name it
  Cascade_Sigma_Sim_10M_240221_RGB_Fall2018_Inbending.Branch("eventno",&eventno);
  Cascade_Sigma_Sim_10M_240221_RGB_Fall2018_Inbending.Branch("runno",&runno,"runno/I");
  Cascade_Sigma_Sim_10M_240221_RGB_Fall2018_Inbending.Branch("triggerno",&triggerno,"triggerno/I");
  Cascade_Sigma_Sim_10M_240221_RGB_Fall2018_Inbending.Branch("start_time",&start_time);
  Cascade_Sigma_Sim_10M_240221_RGB_Fall2018_Inbending.Branch("p4",&v_p4);
  Cascade_Sigma_Sim_10M_240221_RGB_Fall2018_Inbending.Branch("vertex",&v_vertex);
  Cascade_Sigma_Sim_10M_240221_RGB_Fall2018_Inbending.Branch("beta",&v_beta);
  Cascade_Sigma_Sim_10M_240221_RGB_Fall2018_Inbending.Branch("status",&v_status);
  Cascade_Sigma_Sim_10M_240221_RGB_Fall2018_Inbending.Branch("energy",&v_energy);
  Cascade_Sigma_Sim_10M_240221_RGB_Fall2018_Inbending.Branch("charge",&v_charge);
  Cascade_Sigma_Sim_10M_240221_RGB_Fall2018_Inbending.Branch("PID",&v_PID);
  Cascade_Sigma_Sim_10M_240221_RGB_Fall2018_Inbending.Branch("chi2PID",&v_chi2PID);
  Cascade_Sigma_Sim_10M_240221_RGB_Fall2018_Inbending.Branch("region",&v_region);
  Cascade_Sigma_Sim_10M_240221_RGB_Fall2018_Inbending.Branch("time",&v_time);
  Cascade_Sigma_Sim_10M_240221_RGB_Fall2018_Inbending.Branch("path",&v_path);
  Cascade_Sigma_Sim_10M_240221_RGB_Fall2018_Inbending.Branch("mass",&v_mass);
  Cascade_Sigma_Sim_10M_240221_RGB_Fall2018_Inbending.Branch("Pos_position",&v_Pos_position);
  Cascade_Sigma_Sim_10M_240221_RGB_Fall2018_Inbending.Branch("beam",&beam);
  Cascade_Sigma_Sim_10M_240221_RGB_Fall2018_Inbending.Branch("target",&target);
  Cascade_Sigma_Sim_10M_240221_RGB_Fall2018_Inbending.Branch("positive",&positive);
  Cascade_Sigma_Sim_10M_240221_RGB_Fall2018_Inbending.Branch("negative",&negative);


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
      v_mass.clear();
      v_Pos_position.clear();

      // Define how to access the information for each event
      eventno = c12.runconfig()->getEvent(); // Getting the event number
      runno = c12.runconfig()->getRun(); // Getting the run number
      triggerno = c12.runconfig()->getTrigger(); // Getting the trigger bit
      start_time = c12.event()->getStartTime(); // Getting start time for each event
      Nparticles = 0; // Count the number of particles in an event (charged only)

      kaonpno = 0; // Setting the kaon number to 0 abetat start of each event
      elno = 0; // Setting the electron number to 0 at start of each event
      positive = 0; // Setting the number of positives to 0 at start of each event
      negative = 0; // Setting the number of negatives to 0 at start of each event

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
        if(p->getRegion()==FT){
          Region=0;
        }
        if(p->getRegion()==FD){
          Region=1;
        }
        if(p->getRegion()==CD){
          Region=2;
        }

        // To assign particle mass we check the PID and use PDG database
        // Skip neutral particles
        if(charge==0)continue;

        // Counting the number of charged particles
        else if(charge > 0)positive++;
        else if(charge < 0)negative++;

        // Electrons
        if(PID==11){
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
          p4.SetXYZM(p4.Px(),p4.Py(),p4.Pz(),db->GetParticle(11)->Mass());
          v_p4.push_back(p4); // Recording the 4 vector for each electron
          v_region.push_back(Region);
          elno++; // Increasing the count of electrons
        }

        // protons
        else if(PID==2212){
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
          p4.SetXYZM(p4.Px(),p4.Py(),p4.Pz(),db->GetParticle(2212)->Mass());
          v_p4.push_back(p4); // Recording the 4 vector for each proton
          v_region.push_back(Region);
          // protonno++; // Increasing the count of protons
        }
        // K^{+}
        else if(PID==321){
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
          p4.SetXYZM(p4.Px(),p4.Py(),p4.Pz(),db->GetParticle(321)->Mass());
          v_p4.push_back(p4); // Recording the 4 vector for each K^{+}
          v_region.push_back(Region);
          kaonpno++; // Increasing the count of K^{+}
        }
        // K^{-}
        else if(PID==-321){
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
          p4.SetXYZM(p4.Px(),p4.Py(),p4.Pz(),db->GetParticle(-321)->Mass());
          v_p4.push_back(p4); // Recording the 4 vector for each K^{-}
          v_region.push_back(Region);
          // kaonmno++; // Increasing the count of K^{-}
        }
        else if(PID==211){
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
          p4.SetXYZM(p4.Px(),p4.Py(),p4.Pz(),db->GetParticle(211)->Mass());
          v_p4.push_back(p4); // Recording the 4 vector for each pi^{+}
          v_region.push_back(Region);
          // pipno++; // Increasing the count of pi^{+}
        }
        else if(PID==-211){
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
          p4.SetXYZM(p4.Px(),p4.Py(),p4.Pz(),db->GetParticle(-211)->Mass());
          v_p4.push_back(p4); // Recording the 4 vector for each pi^{-}
          v_region.push_back(Region);
          // pimno++; // Increasing the count of pi^{-}
        }
      }
      // Here you can apply a basic skim for events you want to save in your tree
      // Requires at least 2 kaons and 1 electron to be stored in tree
      if(kaonpno > 2 && elno == 1){
        //Fill the TTree
        Cascade_Sigma_Sim_10M_240221_RGB_Fall2018_Inbending.Fill();
        GoodEvents++;
      }
    }

  }
  auto finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = finish - start;
  std::cout << "Elapsed time: " << elapsed.count()<< " events = "<<counter<< " s\n";
  cout<<"Good Events: "<<GoodEvents<<endl;

  Cascade_Sigma_Sim_10M_240221_RGB_Fall2018_Inbending.Write();
  f.Close();
}
