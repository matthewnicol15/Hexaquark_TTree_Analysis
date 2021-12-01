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


void Tree_Event_Selection_Test(){
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

  // Access the PDG database to get information usind PID (e.g db->GetParticle(211)->Mass() for pi^{+} mass)
  auto db=TDatabasePDG::Instance();

  // Any information specific to an individual particle
  TLorentzVector p4; // TLorentzVector for four vector of each particle (Px,Py,Pz,E)
  Double_t beta; // beta from time of flight (TOF)
  Double_t status; // Gives information on which detectors and how many picked it up
  Double_t energy; // Energy measured by detector
  Double_t charge; // Charge measured by drift chambers
  Double_t PID; // Records PID deterimined by TOF
  Double_t chi2PID; // Chi^2 for PID of TOF
  Double_t time; // Time recorded by by FTCAL,FTOF or CTOF
  Double_t path; // Path measured by FTCAL,FTOF or CTOF
  Double_t Region; // Records 0.0 for FT, 1.0 for FD and 2.0 for CD
  Double_t Mass; // Records the calculated mass from beta and momentum

  // Record the number of particles in an event
  Int_t elno; // electrons
  Int_t positronno; // positrons
  Int_t protonno; // protons
  Int_t antiprotonno; // antiprotons
  Int_t neutronno; // neutrons
  Int_t photonno; // photons
  Int_t pipno; // pi^{+}
  Int_t pimno; // pi^{-}
  Int_t kaonpno; // K^{+}
  Int_t kaonpFD; // K^{+} in FD
  Int_t electronFD; // e' in FD
  Int_t kaonmno; // K^{-}
  Int_t positive_charge_tracks; // Count the number of positive charge tracks
  Int_t negative_charge_tracks; // Count the number of negative charge tracks
  Int_t eFD_Events = 0; // Count the number of negative charge tracks
  Int_t topology_1 = 0, topology_2 = 0, topology_3 = 0; // Count the number of negative charge tracks


  // Going over all the input files listed above
  for(Int_t i=0; i<10; i++){

    // Create the CLAS12 event reader
    clas12reader c12(files->At(i)->GetTitle());


    // This loop goes over the events within each file
    while(c12.next()==true){
      counter++;

      elno = 0;
      protonno = 0;
      neutronno = 0;
      pimno = 0;
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
        status = p->par()->getStatus(); // Information on which detectors are involved for this particle
        energy =  p->getDetEnergy(); // Energy measured by detector
        charge = p->par()->getCharge(); // Charge of the particle measured (FTB isn't recognised!!!)
        PID = p->par()->getPid();  // PID determined by TOF
        chi2PID = p->par()->getChi2Pid(); // Chi^2 for PID from TOF
        Mass = sqrt((pow(p4.Rho(),2) / (pow(beta,2))) - pow(p4.Rho(),2));

        // Looking at positive particles
        if(charge > 0){
          if(PID==321) kaonpno++; // Count the number of positive kaons
          else if(PID==2212)protonno++; // Count the number of protons
          positive_charge_tracks++; // Count number of positive charges
        }

        // Looking at negative particles
        else if(charge < 0) {
          if(PID == 11) elno++; // Count number of electrons

          if(PID==-211) pimno++; // Count the number of negative pions
          if(PID==-321) kaonmno++; // Count the number of negative pions
          negative_charge_tracks++;
        }


        // Recording the region the particles hit
        if(p->getRegion()==FT){
          Region=0.0;
        }
        else if(p->getRegion()==FD){
          Region=1.0;
          if(PID == 321)kaonpFD++;
          else if(PID == 11)electronFD++;
        }
        else if(p->getRegion()==CD){
          Region=2.0;
        }
      }


      // Count the number of events in each topology
      if(electronFD == 1){
        eFD_Events++;
        if(kaonpno > 2) topology_1++;
        if(protonno > 0 && pimno > 1) topology_2++;
        if(positive_charge_tracks > 0 && negative_charge_tracks > 2) topology_3++;
      }
    }
  }
  auto finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = finish - start;
  std::cout << "Elapsed time: " << elapsed.count()<< " events = "<<counter<< " s\n";
  cout<<"eFD "<<eFD_Events<<" topology 1 "<<topology_1<<" topology 2 "<<topology_2<<" topology 3 "<<topology_3<<endl;
}
