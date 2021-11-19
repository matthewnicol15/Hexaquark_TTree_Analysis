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
#include <vector>

void Mass_PID_RGA(){

  // Read file with information on vectors
  gROOT->ProcessLine(".L ./Loader.C+");

  // Read input root file and assign it to 'f'
  TFile *f = new TFile("/mnt/f/PhD/Trees/RGA_Fall2018_Inbending_Pass1_Skim4_0050_e_2pos_1neg_Tree_100920_01.root");
  // Read TTree within root file and assign it to 't1'
  TTree *t1 = (TTree*)f->Get("skim4_Tree_100920_01");

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
  Int_t readprotonno; // Protons
  Int_t readpipno; // pi^+
  Int_t readpimno; // pi^-
  Int_t readelno; // e^-
  Int_t readkaonpno; // K^+
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
  // t1->SetBranchAddress("protonno",&readprotonno);
  // t1->SetBranchAddress("pipno",&readpipno);
  // t1->SetBranchAddress("pimno",&readpimno);
  t1->SetBranchAddress("elno",&readelno);
  // t1->SetBranchAddress("kaonpno",&readkaonpno);
  t1->SetBranchAddress("eventno",&readeventno);
  t1->SetBranchAddress("runno",&readrunno);
  t1->SetBranchAddress("triggerno",&readtriggerno);

  // Path and name for the output file to save
  // TFile fileOutput1("../New_Data/Mass_Skim4_0050_Pass1_nop_nopi_110920_01.root","recreate");

  // Getting particle database to use for masses
  auto db=TDatabasePDG::Instance();


  //Creating histograms

  // Test histograms
  auto* helno=new TH1F("helno","Number of electrons",5,0,5);
  auto* hprotonno=new TH1F("hprotonno","Number of protons",5,0,5);
  auto* hkaonpno=new TH1F("hkaonpno","Number of K^{+}",5,0,5);
  auto* hpimno=new TH1F("hpimno","Number of #pi^{-}",5,0,5);
  auto* hregion=new TH1F("hregion","Regions;Region;Counts",4,-0.5,3.5);
  auto* hregion_kp=new TH1F("hregion_kp","Regions;Region;Counts",4,-0.5,3.5);

  auto* hmass=new TH1F("hmass","Calculated mass",700,-0.2,1.2);
  auto* hmass_a=new TH1F("hmass_a","Calculated mass in kaon region",700,-0.2,1.2);
  auto* hmass_b=new TH1F("hmass_b","Calculated mass in kaon region",700,-0.2,1.2);
  auto* hmass_c=new TH1F("hmass_c","Calculated mass in kaon region",700,-0.2,1.2);
  auto* hmass_d=new TH1F("hmass_d","Calculated mass in kaon region",700,-0.2,1.2);
  auto* hmass_e=new TH1F("hmass_e","Calculated mass in kaon region",700,-0.2,1.2);

  auto* hbeta=new TH2D("hbeta","Beta measured", 500,0,11,400,-1,1);
  auto* hdelta_beta_el=new TH2D("hdelta_beta_el"," Delta Beta e^{-}", 500,0,11,400,-1,1);
  auto* hdelta_beta_pim=new TH2D("hdelta_beta_pim"," Delta Beta #pi^{-}", 500,0,11,400,-1,1);
  auto* hdelta_beta_km=new TH2D("hdelta_beta_km"," Delta Beta K^{-}", 500,0,11,400,-1,1);
  auto* hdelta_beta_kp=new TH2D("hdelta_beta_kp"," Delta Beta K^{+}", 500,0,11,400,-1,1);
  auto* hdelta_beta_pr=new TH2D("hdelta_beta_pr"," Delta Beta p", 500,0,11,400,-1,1);
  auto* hdelta_beta_pip=new TH2D("hdelta_beta_pip"," Delta Beta #pi^{+}", 500,0,11,400,-1,1);

  //PID only cut
  auto* hmiss_mass_all=new TH1F("miss_all","MM^2(e' K^{+} p #pi^{-});MM^2(e' K^{+} p #pi^{-}) [GeV];Counts",200,-1,3);
  auto* hmiss_momentum_all=new TH1F("hmiss_momentum_all","P(B + T - e' - K^{+} - p - #pi^{-});P(B + T - e' - K^{+} - p - #pi^{-}) [GeV];Counts",200,-1,3);
  auto* hmiss_1=new TH1F("miss_1","MM(e' K^{+});MM(e' K^{+}) [GeV];Counts",200,-1,3);
  auto* hmiss_1_a=new TH1F("hmiss_1_a","MM(e' K^{+});MM(e' K^{+}) [GeV];Counts",200,-1,3);
  auto* hmiss_1_b=new TH1F("hmiss_1_b","MM(e' K^{+});MM(e' K^{+}) [GeV];Counts",200,-1,3);
  auto* hmiss_1_c=new TH1F("hmiss_1_c","MM(e' K^{+});MM(e' K^{+}) [GeV];Counts",200,-1,3);
  auto* hmiss_1_d=new TH1F("hmiss_1_d","MM(e' K^{+});MM(e' K^{+}) [GeV];Counts",200,-1,3);
  auto* hmiss_1_e=new TH1F("hmiss_1_e","MM(e' K^{+});MM(e' K^{+}) [GeV];Counts",200,-1,3);
  auto* hmiss_1_2=new TH1F("hmiss_1_2","MM(e' K^{+});MM(e' K^{+}) [GeV];Counts",200,-1,3);
  auto* hmiss_3=new TH1F("hmiss_3","MM(e' K^{+} K^{+});MM(e' K^{+} K^{+}) [GeV];Counts",200,-1,3);
  auto* hmiss_3_f=new TH1F("hmiss_3_f","MM(e' K^{+} K^{+});MM(e' K^{+} K^{+}) [GeV];Counts",200,-1,3);
  auto* hmiss_3_g=new TH1F("hmiss_3_g","MM(e' K^{+} K^{+});MM(e' K^{+} K^{+}) [GeV];Counts",200,-1,3);
  auto* hmiss_5=new TH1F("hmiss_5","MM(e' K^{+} K^{+} K^{+});MM(e' K^{+} K^{+} K^{+}) [GeV];Counts",600,-3,3);
  auto* hlambda_1=new TH1F("lambda_1","M(p #pi^{-});M(p #pi^{-}) [GeV];Counts",750,0.5,2.0);
  auto* hS1_kaon_beta=new TH2D("hS1_kaon_beta","#Delta#Beta K^{+};P [GeV];#Delta#Beta K^{+}",200,0,11,100,-0.05,0.05);
  auto* hS1_kaon_beta_1=new TH2D("hS1_kaon_beta_1","#Delta#Beta K^{+};P [GeV];#Delta#Beta K^{+}",200,0,11,100,-0.05,0.05);
  auto* hproton_beta=new TH2D("hproton_beta","#Delta#Beta p;P [GeV];#Delta#Beta p",200,0,11,100,-0.05,0.05);
  auto* hproton_beta_1=new TH2D("hproton_beta_1","#Delta#Beta p;P [GeV];#Delta#Beta p",200,0,11,100,-0.05,0.05);
  auto* hpion_beta=new TH2D("hpion_beta","#Delta#Beta #pi^{-};P [GeV];#Delta#Beta #pi^{-}",200,0,11,100,-0.05,0.05);
  auto* hpion_beta_1=new TH2D("hpion_beta_1","#Delta#Beta #pi^{-};P [GeV];#Delta#Beta #pi^{-}",200,0,11,100,-0.05,0.05);
  auto* hS2_kaon1_beta=new TH2D("hS2_kaon1_beta","#Delta#Beta K^{+} (1));P [GeV];#Delta#Beta K^{+} (1)",200,0,11,100,-0.05,0.05);
  auto* hS2_kaon2_beta=new TH2D("hS2_kaon2_beta","#Delta#Beta K^{+} (2));P [GeV],#Delta#Beta K^{+} (2)",200,0,11,100,-0.05,0.05);

  // delta beta and momentum cuts on kaons and protons
  auto* hmiss_mass_allc=new TH1F("miss_allc","MM^2(e' K^{+} p #pi^{-}) after cuts on delta beta and momentum of K^{+} and p;MM^2(e' K^{+} p #pi^{-}) [GeV];Counts",200,-1,3);
  auto* hmiss_momentum_allc=new TH1F("hmiss_momentum_allc","P(B + T - e' - K^{+} - p - #pi^{-}) after cuts on delta beta and momentum of K^{+} and p;P(B + T - e' - K^{+} - p - #pi^{-}) [GeV];Counts",200,-1,3);
  auto* hmiss_1c=new TH1F("miss_1c","MM(e' K^{+}) after cuts on delta beta and momentum of K^{+} and p;MM(e' K^{+}) [GeV];Counts",200,-1,3);
  auto* hmiss_3c=new TH1F("miss_3c","MM(e' K^{+} K^{+}) after cuts on delta beta and momentum of K^{+} and p;MM(e' K^{+} p) [GeV];Counts",200,-1,3);
  auto* hlambda_1c=new TH1F("lambda_1c","M(p #pi^{-}) after cuts on delta beta and momentum of K^{+} and p;M(p #pi^{-}) [GeV];Counts",750,0.5,2.0);

  //Exclusivity cuts as well
  auto* hmiss_1t=new TH1F("miss_1t","MM(e' K^{+}) after cuts on delta beta and momentum of K^{+} and p and MM^2(e' K^{+} p) cut;MM(e' K^{+}) [GeV];Counts",200,-1,3);
  auto* hmiss_3t=new TH1F("miss_3t","MM(e' K^{+} K^{+}) after cuts on delta beta and momentum of K^{+} and p and MM(e' K^{+} p) cut;MM^2(e' K^{+} p) [GeV];Counts",200,-1,3);
  auto* hlambda_1t=new TH1F("lambda_1t","M(p #pi^{-}) after cuts on delta beta and momentum of K^{+} and p and MM^2(e' K^{+} p) cut;M(p #pi^{-}) [GeV];Counts",750,0.5,2.0);


  vector<TLorentzVector> v_el;
  vector<TLorentzVector> v_pip;
  vector<TLorentzVector> v_pim;
  vector<TLorentzVector> v_pr;
  vector<TLorentzVector> v_kp;
  vector<TLorentzVector> v_kp_b;
  vector<TLorentzVector> v_km;
  vector<TLorentzVector> v_kp_pi;
  vector<TLorentzVector> v_othertracks;
  vector<TLorentzVector> v_vertex_el;
  vector<TLorentzVector> v_vertex_pip;
  vector<TLorentzVector> v_vertex_pim;
  vector<TLorentzVector> v_vertex_pr;
  vector<TLorentzVector> v_vertex_kp;
  vector<TLorentzVector> v_vertex_kp_b;
  vector<TLorentzVector> v_vertex_km;


  TLorentzVector missall;
  TLorentzVector missallcasc;
  TLorentzVector miss1;
  TLorentzVector miss1_b;
  TLorentzVector miss3;
  TLorentzVector miss5;
  TLorentzVector lambda;
  TLorentzVector vertex_el;
  TLorentzVector vertex_pip;
  TLorentzVector vertex_pim;
  TLorentzVector vertex_pr;
  TLorentzVector vertex_kp;
  TLorentzVector vertex_kp_b;
  TLorentzVector vertex_km;

  Double_t P_p4;
  Double_t beta_tof;
  Double_t beta_calc_pi;
  Double_t delta_beta_pi;
  Double_t beta_calc_kaon;
  Double_t delta_beta_kaon;
  Double_t Mass;

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
  vector<Double_t> v_region_el;  // Region particle is detected
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
  vector<Double_t> v_region_pip;  // Region particle is detected
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
  vector<Double_t> v_region_pim;  // Region particle is detected
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
  vector<Double_t> v_region_pr;  // Region particle is detected
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
  vector<Double_t> v_region_kp;  // Region particle is detected
  Double_t region_kp;

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
  vector<Double_t> v_region_km;  // Region particle is detected
  Double_t region_km;

  vector<Double_t> v_beta_tof_othertracks;
  Double_t beta_tof_othertracks;
  vector<Double_t> v_P_othertracks;
  Double_t P_othertracks;
  vector<Double_t> v_beta_calc_othertracks;
  Double_t beta_calc_othertracks;
  vector<Double_t> v_delta_beta_othertracks;
  Double_t delta_beta_othertracks;
  vector<Double_t> v_vt_othertracks;
  Double_t vt_othertracks;

  Double_t delta_theta, delta_phi_S1_pr_pim, delta_P, Delta_Phi_Unidentified;

  Double_t Delta_Vertex_Time;

  Double_t c=30;

  TLorentzVector el;
  TLorentzVector pip;
  TLorentzVector pim;
  TLorentzVector pr;
  TLorentzVector kp;
  TLorentzVector km;
  TLorentzVector kp_b;
  TLorentzVector othertracks;


  // Long64_t nentries = 800000;
  Long64_t nentries = t1->GetEntries();
  Int_t Percentage = nentries/100;
  for(Long64_t i=0; i<nentries;i++){
    t1->GetEntry(i);
    if (i % Percentage == 0){
      fprintf (stderr, "%lld\r", i/Percentage);
      fflush (stderr);
    }

    cout<<readrunno<<endl;

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

    // cout<<"test1"<<endl;

    Int_t Nparticles = v_p4->size();
    for(Int_t j=0; j<Nparticles; j++){

      // Filling histogram to show the different regions hit
      hregion->Fill(v_region->at(j));

      // Getting the measured momentum and beta for calculating delta beta
      P_p4 = sqrt((pow(v_p4->at(j).Px(),2))+(pow(v_p4->at(j).Py(),2))+(pow(v_p4->at(j).Pz(),2))); // Measured momentum
      beta_tof = v_beta->at(j); // Measured beta values
      // cout<<"test2"<<endl;

      // Calculated beta for various particles
      beta_calc_el = P_p4/(sqrt((pow(P_p4,2))+(pow(db->GetParticle(11)->Mass(),2))));
      beta_calc_pi = P_p4/(sqrt((pow(P_p4,2))+(pow(db->GetParticle(211)->Mass(),2))));
      beta_calc_pr = P_p4/(sqrt((pow(P_p4,2))+(pow(db->GetParticle(2212)->Mass(),2))));
      beta_calc_kaon = P_p4/(sqrt((pow(P_p4,2))+(pow(db->GetParticle(321)->Mass(),2))));
      // cout<<"test3"<<endl;

      // Calculating delta beta for various particles
      delta_beta_el = beta_calc_el-beta_tof;
      delta_beta_pi = beta_calc_pi-beta_tof;
      delta_beta_pr = beta_calc_pr-beta_tof;
      delta_beta_kaon = beta_calc_kaon-beta_tof;

      hbeta->Fill(P_p4,v_beta->at(j));

      // Calculating Mass from beta and momentum
      Mass = sqrt((pow(P_p4,2) / (pow(beta_tof,2))) - pow(P_p4,2));

      // cout<<"test4"<<endl;
      if(v_charge->at(j) < 0){
        // cout<<"test5"<<endl;

        hdelta_beta_el->Fill(P_p4,delta_beta_el);
        hdelta_beta_pim->Fill(P_p4,delta_beta_pi);
        hdelta_beta_km->Fill(P_p4,delta_beta_kaon);

        if(v_PID->at(j) == 11){
          // cout<<"test6"<<endl;
          el.SetXYZM(v_p4->at(j).Px(),v_p4->at(j).Py(),v_p4->at(j).Pz(),db->GetParticle(11)->Mass());
          P_el = sqrt((pow(v_p4->at(j).Px(),2))+(pow(v_p4->at(j).Py(),2))+(pow(v_p4->at(j).Pz(),2))); // Measured momentum
          TOF_el = v_time->at(j);
          path_el = v_path->at(j);
          beta_tof_el = v_beta->at(j);
          beta_calc_el = P_el/(sqrt((pow(P_el,2))+(pow(db->GetParticle(11)->Mass(),2))));
          delta_beta_el = beta_calc_el-beta_tof_el;
          vertex_time_el = TOF_el - path_el / (beta_tof_el*c);
          vertex_el.SetXYZT(v_vertex->at(j).X(), v_vertex->at(j).Y(), v_vertex->at(j).Z(), vertex_time_el);

          v_el.push_back(el);
          v_beta_tof_el.push_back(beta_tof_el);
          v_P_el.push_back(P_el);
          v_path_el.push_back(path_el);
          v_TOF_el.push_back(TOF_el);
          v_beta_calc_el.push_back(beta_calc_el);
          v_delta_beta_el.push_back(delta_beta_el);
          v_vertex_time_el.push_back(vertex_time_el);
          v_vertex_el.push_back(vertex_el);
        }

        else if(Mass > 0.05 && Mass < 0.25){
          // cout<<"test6"<<endl;
          pim.SetXYZM(v_p4->at(j).Px(),v_p4->at(j).Py(),v_p4->at(j).Pz(),db->GetParticle(-211)->Mass());
          P_pim = sqrt((pow(v_p4->at(j).Px(),2))+(pow(v_p4->at(j).Py(),2))+(pow(v_p4->at(j).Pz(),2))); // Measured momentum
          TOF_pim = v_time->at(j);
          path_pim = v_path->at(j);
          beta_tof_pim = v_beta->at(j);
          beta_calc_pim = P_pim/(sqrt((pow(P_pim,2))+(pow(db->GetParticle(-211)->Mass(),2))));
          delta_beta_pim = beta_calc_pim-beta_tof_pim;
          vertex_time_pim = TOF_pim - path_pim / (beta_tof_pim*c);
          vertex_pim.SetXYZT(v_vertex->at(j).X(), v_vertex->at(j).Y(), v_vertex->at(j).Z(), vertex_time_pim);

          v_pim.push_back(pim);
          v_beta_tof_pim.push_back(beta_tof_pim);
          v_P_pim.push_back(P_pim);
          v_path_pim.push_back(path_pim);
          v_TOF_pim.push_back(TOF_pim);
          v_beta_calc_pim.push_back(beta_calc_pim);
          v_delta_beta_pim.push_back(delta_beta_pim);
          v_vertex_time_pim.push_back(vertex_time_pim);
          v_vertex_pim.push_back(vertex_pim);
        }

        else if(Mass > 0.3 && Mass < 0.7){
          // cout<<"test6"<<endl;
          km.SetXYZM(v_p4->at(j).Px(),v_p4->at(j).Py(),v_p4->at(j).Pz(),db->GetParticle(-321)->Mass());
          P_km = sqrt((pow(v_p4->at(j).Px(),2))+(pow(v_p4->at(j).Py(),2))+(pow(v_p4->at(j).Pz(),2))); // Measured momentum
          TOF_km = v_time->at(j);
          path_km = v_path->at(j);
          beta_tof_km = v_beta->at(j);
          beta_calc_km = P_km/(sqrt((pow(P_km,2))+(pow(db->GetParticle(-321)->Mass(),2))));
          delta_beta_km = beta_calc_km-beta_tof_km;
          vertex_time_km = TOF_km - path_km / (beta_tof_km*c);
          vertex_km.SetXYZT(v_vertex->at(j).X(), v_vertex->at(j).Y(), v_vertex->at(j).Z(), vertex_time_km);

          v_km.push_back(km);
          v_beta_tof_km.push_back(beta_tof_km);
          v_P_km.push_back(P_km);
          v_path_km.push_back(path_km);
          v_TOF_km.push_back(TOF_km);
          v_beta_calc_km.push_back(beta_calc_km);
          v_delta_beta_km.push_back(delta_beta_km);
          v_vertex_time_km.push_back(vertex_time_km);
          v_vertex_km.push_back(vertex_km);
        }
      }

      if(v_charge->at(j) > 0){

        hdelta_beta_pr->Fill(P_p4,delta_beta_pr);
        hdelta_beta_pip->Fill(P_p4,delta_beta_pi);

        if(Mass > 0.3 && Mass < 0.7){
          // cout<<"test6"<<endl;
          kp.SetXYZM(v_p4->at(j).Px(),v_p4->at(j).Py(),v_p4->at(j).Pz(),db->GetParticle(321)->Mass());
          P_kp = sqrt((pow(v_p4->at(j).Px(),2))+(pow(v_p4->at(j).Py(),2))+(pow(v_p4->at(j).Pz(),2))); // Measured momentum
          TOF_kp = v_time->at(j);
          path_kp = v_path->at(j);
          beta_tof_kp = v_beta->at(j);
          beta_calc_kp = P_kp/(sqrt((pow(P_kp,2))+(pow(db->GetParticle(321)->Mass(),2))));
          delta_beta_kp = beta_calc_kp-beta_tof_kp;
          vertex_time_kp = TOF_kp - path_kp / (beta_tof_kp*c);
          vertex_kp.SetXYZT(v_vertex->at(j).X(), v_vertex->at(j).Y(), v_vertex->at(j).Z(), vertex_time_kp);
          region_kp = 1.0*v_region->at(j); // Regions are now FT=0, FD=1.5, CD=3

          hregion_kp->Fill(region_kp);
          hmass_a->Fill(Mass);
          if(region_kp < 0.5){
            hmass_b->Fill(Mass);
          }
          else if(region_kp < 1.5){
            hmass_c->Fill(Mass);
          }
          else if(region_kp < 2.5){
            hmass_d->Fill(Mass);
          }
          else{hmass_e->Fill(Mass);
          }

          hdelta_beta_kp->Fill(P_p4,delta_beta_kp);

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
        }

        else if(Mass > 0.7 && Mass < 1.2){
          // cout<<"test6"<<endl;
          pr.SetXYZM(v_p4->at(j).Px(),v_p4->at(j).Py(),v_p4->at(j).Pz(),db->GetParticle(2212)->Mass());
          P_pr = sqrt((pow(v_p4->at(j).Px(),2))+(pow(v_p4->at(j).Py(),2))+(pow(v_p4->at(j).Pz(),2))); // Measured momentum
          TOF_pr = v_time->at(j);
          path_pr = v_path->at(j);
          beta_tof_pr = v_beta->at(j);
          beta_calc_pr = P_pr/(sqrt((pow(P_pr,2))+(pow(db->GetParticle(2212)->Mass(),2))));
          delta_beta_pr = beta_calc_pr-beta_tof_pr;
          vertex_time_pr = TOF_pr - path_pr / (beta_tof_pr*c);
          vertex_pr.SetXYZT(v_vertex->at(j).X(), v_vertex->at(j).Y(), v_vertex->at(j).Z(), vertex_time_pr);

          v_pr.push_back(pr);
          v_beta_tof_pr.push_back(beta_tof_pr);
          v_P_pr.push_back(P_pr);
          v_path_pr.push_back(path_pr);
          v_TOF_pr.push_back(TOF_pr);
          v_beta_calc_pr.push_back(beta_calc_pr);
          v_delta_beta_pr.push_back(delta_beta_pr);
          v_vertex_time_pr.push_back(vertex_time_pr);
          v_vertex_pr.push_back(vertex_pr);
        }
        else if(Mass > 0.02 && Mass < 0.3){
          // cout<<"test6"<<endl;
          pip.SetXYZM(v_p4->at(j).Px(),v_p4->at(j).Py(),v_p4->at(j).Pz(),db->GetParticle(211)->Mass());
          P_pip = sqrt((pow(v_p4->at(j).Px(),2))+(pow(v_p4->at(j).Py(),2))+(pow(v_p4->at(j).Pz(),2))); // Measured momentum
          TOF_pip = v_time->at(j);
          path_pip = v_path->at(j);
          beta_tof_pip = v_beta->at(j);
          beta_calc_pip = P_pip/(sqrt((pow(P_pip,2))+(pow(db->GetParticle(211)->Mass(),2))));
          delta_beta_pip = beta_calc_pip-beta_tof_pip;
          vertex_time_pip = TOF_pip - path_pip / (beta_tof_pip*c);
          vertex_pip.SetXYZT(v_vertex->at(j).X(), v_vertex->at(j).Y(), v_vertex->at(j).Z(), vertex_time_pip);

          v_pip.push_back(pip);
          v_beta_tof_pip.push_back(beta_tof_pip);
          v_P_pip.push_back(P_pip);
          v_path_pip.push_back(path_pip);
          v_TOF_pip.push_back(TOF_pip);
          v_beta_calc_pip.push_back(beta_calc_pip);
          v_delta_beta_pip.push_back(delta_beta_pip);
          v_vertex_time_pip.push_back(vertex_time_pip);
          v_vertex_pip.push_back(vertex_pip);
        }


      }

      if(v_charge->at(j) != 0)hmass->Fill(Mass);

    }

    helno->Fill(v_el.size());
    hpimno->Fill(v_pim.size());
    hkaonpno->Fill(v_kp.size());
    hprotonno->Fill(v_pr.size());



    // Strangeness 1 channels
    if(v_kp.size()==1 && v_el.size()==1 /*&& v_pr.size()==1 && v_pim.size()==1*/){

      // missall = (TLorentzVector)*readbeam + (TLorentzVector)*readtarget - v_el.at(0) - v_kp.at(0) - v_pr.at(0) - v_pim.at(0);
      miss1 = (TLorentzVector)*readbeam + (TLorentzVector)*readtarget - v_el.at(0) - v_kp.at(0);
      // lambda = v_pr.at(0) + v_pim.at(0);

      // hlambda_1->Fill(lambda.M());
      hmiss_1->Fill(miss1.M());
      // hmiss_mass_all->Fill(missall.M2());
      // hmiss_momentum_all->Fill(missall.Rho());

      if(fabs(v_delta_beta_kp.at(0))<0.01 && v_P_kp.at(0)>0.9 && v_P_kp.at(0)<3 /*&& fabs(v_delta_beta_pr.at(0))<0.01 && v_P_pr.at(0)>0.2 && v_P_pr.at(0)<4.8 && fabs(v_delta_beta_pim.at(0))<0.01 && v_P_pim.at(0)>0.2 && v_P_pim.at(0)<3*/){
        hmiss_1_a->Fill(miss1.M()); // All events

        // Filling histograms depending on where the K^+ goes
        if(v_region_kp.at(0) < 0.5){
          hmiss_1_b->Fill(miss1.M()); // Kaons in FT
        }
        else if(v_region_kp.at(0) < 1.5){
          hmiss_1_c->Fill(miss1.M()); // Kaons in FD
        }
        else if(v_region_kp.at(0) < 2.5){
          hmiss_1_d->Fill(miss1.M()); // Kaons in CD
        }
        else{
          hmiss_1_e->Fill(miss1.M()); // Kaons anywhere else
        }
        // hlambda_1c->Fill(lambda.M());
        // hmiss_mass_allc->Fill(missall.M2());
        // hmiss_momentum_allc->Fill(missall.Rho());

        // if(fabs(missall.M2())<0.05 && fabs(missall.Rho())<0.1){
        //   hmiss_1t->Fill(miss1.M());
        //   hlambda_1t->Fill(lambda.M());
        // }
        }
      }

      // Strangeness 2 channels
      if(v_kp.size()==2  && v_el.size()==1 /*&& v_pr.size()==1 && v_pim.size()==2*/){
        // missallcasc = (TLorentzVector)*readbeam + (TLorentzVector)*readtarget - v_el.at(0) - v_kp.at(0) - v_kp.at(1) - v_pr.at(0) - v_pim.at(0) - v_pim.at(1);
        miss3 = (TLorentzVector)*readbeam + (TLorentzVector)*readtarget - v_el.at(0) - v_kp.at(0) - v_kp.at(1);
        // cascade_lambda_1 = v_pr.at(0) + v_pim.at(0);
        // cascade_lambda_2 = v_pr.at(0) + v_pim.at(1);
        Delta_Vertex_Time = v_vertex_kp.at(0).T() - v_vertex_kp.at(1).T();
        // inv_cascade = v_pr.at(0) + v_pim.at(0) + v_pim.at(1);

        hmiss_3->Fill(miss3.M());
        // Filling histograms depending on where the K^+ goes
        if(v_region_kp.at(0) < 1.5 && v_region_kp.at(0) > 0.5){
          hmiss_3_f->Fill(miss3.M()); // Kaons in FT
          if(v_region_kp.at(1) < 1.5 && v_region_kp.at(1) > 0.5){
            hmiss_3_g->Fill(miss3.M()); // Kaons in FD
          }
        }


      }

      // Strangeness 3 channels with electron
      if(v_kp.size()==3 && v_el.size()==1){
        miss5 = (TLorentzVector)*readbeam + (TLorentzVector)*readtarget - v_el.at(0) - v_kp.at(0) - v_kp.at(1) - v_kp.at(2);
        // Filling histograms with PID cuts only
        hmiss_5->Fill(miss5.M());
        // if(fabs(v_delta_beta_kp.at(0))<0.01 && v_P_kp.at(0)>0.5 && v_P_kp.at(0)<2.6 &&
        // fabs(v_delta_beta_kp.at(1))<0.01 && v_P_kp.at(1)>0.5 && v_P_kp.at(1)<2.6 &&
        // fabs(v_delta_beta_kp.at(2))<0.01 && v_P_kp.at(2)>0.5 && v_P_kp.at(2)<2.6){
        //
        //   // Filling histograms with PID, delta beta and momentum cuts
        //   hmiss_5c->Fill(miss5.M());
        // }
      }
    }
    // fileOutput1.Write();
  }
