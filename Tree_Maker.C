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


void Tree_Maker(){
  auto start = std::chrono::high_resolution_clock::now();
  gBenchmark->Start("timer");
  int counter=0;

  // Data files to process
  TString inputFile("/home/matthewn/links/RGB_Spring_2020_Inbending_dst/*.hipo");


  gROOT->ProcessLine(".L /home/matthewn/Documents/Macros/Loader.C+"); // Uses Loader.C file, make sure Loader.C is in this file path


  // Creating a chain for the data from different files
  TChain fake("hipo");
  fake.Add(inputFile.Data());

  // Shortcut to a list of all the input file names
  auto files=fake.GetListOfFiles();

  // Create root file to save TTree in
  TFile f("/lustre19/expphy/volatile/clas12/matthewn/Trees/Dibaryon/RGB_Spring2020_Inbending_1eFD_at_least_1p2pim_Tree_02122021_01.root","recreate");
  // Creating TTree object
  TTree RGB_Spring2020_Inbending_221021("RGB_Spring2020_Inbending_02122021","it's a tree!");

  // Access the PDG database to get information usind PID (e.g db->GetParticle(211)->Mass() for pi^{+} mass)
  auto db=TDatabasePDG::Instance();

  // Creating histograms for negative particles that are not electron
  auto* hnegatives=new TH1F("hnegatives","PID of negatives;PID;Counts",2000,-1000,1000);
  auto* hpositives=new TH1F("hpositives","PID of positives;PID;Counts",2000,-1000,1000);
  auto* hkaonpno=new TH1F("hkaonpno","Number of K^{+};No. of K^{+};Counts",100,0,10);
  auto* hPID=new TH1F("hPID","PID ;PID;Counts",2000,-1000,1000);
  auto* hmass_n=new TH1F("hmass_n","Mass of negatives;Mass [GeV];Counts",1000,-10000,10000);
  auto* hmass_p=new TH1F("hmass_p","Mass of negatives;Mass [GeV];Counts",1000,-10000,10000);
  auto* hbeta_n=new TH2D("hbeta_n","beta against momentum of positive PID 0 particles;Momentum [GeV]; Beta",200,0,11,1000,-100,10);
  auto* hbeta_p=new TH2D("hbeta_n","beta against momentum of positive PID 0 particles;Momentum [GeV]; Beta",200,0,11,1000,-100,10);
  auto* hangular_distribution_n=new TH2D("hangular_distribution_n","Angular distribution of negative PID 0 particles;Momentum [GeV]; Beta",200,0,11,200,0,200);
  auto* hangular_distribution_p=new TH2D("hangular_distribution_p","Angular distribution of positive PID 0 particles;Momentum [GeV]; Beta",200,0,11,200,0,200);

  // Information to save to the tree
  // Any information specific to an event
  Int_t eventno; // Records the event number
  Int_t runno; // Records the run number
  Int_t triggerno; // Records the trigger number
  Double_t start_time; // Records the start time for the event
  Int_t Tree_Events=0; // Records the number of events that are put in the tree

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
  Double_t Mass; // Records the calculated mass from beta and momentum


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
  Int_t kaonpFD; // K^{+} in FD
  Int_t electronFD; // e' in FD
  Int_t kaonmno; // K^{-}
  Int_t positive_charge_tracks; // Count the number of positive charge tracks
  Int_t negative_charge_tracks; // Count the number of negative charge tracks

  // Setting TLorentzVectors for beam and target, for final analysis beam energy
  // will have to be changed depending on which run/runs you are analysing. Here
  // it is just set to 10.6 GeV
  TLorentzVector beam(0,0,0,0); // Set 4 vector four the beam, all momentum in z-direction
  TLorentzVector target(0,0,0,0.93827); // Set 4 vector for target, stationary so no momentum
  // TLorentzVector target(0,0,0,1.8756); // Set 4 vector for target, stationary so no momentum

  Double_t c = 30; // Speed of light, used to calculate vertex time

  // Assign a branch to each measurable and name it
  RGB_Spring2020_Inbending_221021.Branch("eventno",&eventno);
  RGB_Spring2020_Inbending_221021.Branch("runno",&runno,"runno/I");
  RGB_Spring2020_Inbending_221021.Branch("triggerno",&triggerno,"triggerno/I");
  RGB_Spring2020_Inbending_221021.Branch("start_time",&start_time);
  RGB_Spring2020_Inbending_221021.Branch("p4",&v_p4);
  RGB_Spring2020_Inbending_221021.Branch("vertex",&v_vertex);
  RGB_Spring2020_Inbending_221021.Branch("beta",&v_beta);
  RGB_Spring2020_Inbending_221021.Branch("status",&v_status);
  RGB_Spring2020_Inbending_221021.Branch("energy",&v_energy);
  RGB_Spring2020_Inbending_221021.Branch("charge",&v_charge);
  RGB_Spring2020_Inbending_221021.Branch("PID",&v_PID);
  RGB_Spring2020_Inbending_221021.Branch("chi2PID",&v_chi2PID);
  RGB_Spring2020_Inbending_221021.Branch("region",&v_region);
  RGB_Spring2020_Inbending_221021.Branch("time",&v_time);
  RGB_Spring2020_Inbending_221021.Branch("path",&v_path);
  RGB_Spring2020_Inbending_221021.Branch("beam",&beam);
  RGB_Spring2020_Inbending_221021.Branch("target",&target);
  RGB_Spring2020_Inbending_221021.Branch("elno",&elno);
  RGB_Spring2020_Inbending_221021.Branch("negative_charge_tracks",&negative_charge_tracks, "neg");
  RGB_Spring2020_Inbending_221021.Branch("positive_charge_tracks",&positive_charge_tracks);
  RGB_Spring2020_Inbending_221021.Branch("kaonpFD",&kaonpFD);
  RGB_Spring2020_Inbending_221021.Branch("electronFD",&electronFD);


  // Going over all the input files listed above
  for(Int_t i=30000; i<files->GetEntries(); i++){

    // Create the CLAS12 event reader
    clas12reader c12(files->At(i)->GetTitle());


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

      // Define how to access the information for each event
      eventno = c12.runconfig()->getEvent(); // Getting the event number
      runno = c12.runconfig()->getRun(); // Getting the run number
      triggerno = c12.runconfig()->getTrigger(); // Getting the trigger bit
      start_time = c12.event()->getStartTime(); // Getting start time for each event

      // Setting beam energy depending on which run it is
      // if(runno < 6400) beam.SetXYZM(0,0,10.6,0);
      // if(runno > 6400) beam.SetXYZM(0,0,10.2,0);
      beam.SetXYZM(0,0,10.2,0);

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
      kaonpFD = 0;
      electronFD = 0;
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
        Mass = sqrt((pow(p4.Rho(),2) / (pow(beta,2))) - pow(p4.Rho(),2));

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
        // cout<<p4.Rho()<<endl;
        v_p4.push_back(p4); // Recording the 4 vector for each neutron

        hPID->Fill(PID);

        if(PID == 11) elno++; // Count number of electrons

        // Looking at positive particles
        if(charge > 0){
          if(PID==321) kaonpno++; // Count the number of positive kaons
          if(PID==2212)protonno++; // Count the number of protons
          positive_charge_tracks++; // Count number of positive charges
          hpositives->Fill(PID); // Seeing PID of positive particles

          // Calculated mass for positive particles
          hmass_p->Fill(Mass);
        }

        // Looking at negative particles
        else if(charge < 0) {
          if(PID==-211) pimno++; // Count the number of negative pions
          if(PID==-321) kaonmno++; // Count the number of negative pions
          negative_charge_tracks++;
          hnegatives->Fill(PID); // Seeing PID of negative particles

          // Calculated mass for negative particles
          hmass_n->Fill(Mass);
        }


        // Recording the region the particles hit
        if(p->getRegion()==FT){
          Region=0.0;
        }
        else if(p->getRegion()==FD){
          Region=1.0;
          if(PID == 321)kaonpFD++;
          if(PID == 11)electronFD++;
        }
        else if(p->getRegion()==CD){
          Region=2.0;
        }
        v_region.push_back(Region); // pushing back the region of all particles in the event
      }

      // Checking number of kaons per event
      hkaonpno->Fill(kaonpno);


      // Here you can apply a basic skim for events you want to save in your tree
      if(protonno > 0 && pimno > 1 && elno == 1 && electronFD == 1){
        RGB_Spring2020_Inbending_221021.Fill();
        Tree_Events++;
      }
    }
    // Prints out the file currently analysed and how many events are in tree
    cout<<"file: "<<i<<"events in tree "<<Tree_Events<<endl;
  }
  RGB_Spring2020_Inbending_221021.Write(); // Write information to the TTree
  auto finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = finish - start;
  std::cout << "Elapsed time: " << elapsed.count()<< " events = "<<counter<< " s\n";
  cout<<"events in tree "<<Tree_Events<<endl;
  f.Write(); // Write information to the root file
  f.Close(); // Close the root file at the end
}
