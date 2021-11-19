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

void Mass_PID_RGB(){

  gROOT->ProcessLine(".L /mnt/f/PhD/Macros/Loader.C+");
  TFile *f = new TFile("/mnt/f/PhD/Trees/RGB_Inc_Spring2019_Inbending_Pass1_v0_1_Pos_Tree_280920_01.root");
  TTree *t1 = (TTree*)f->Get("RGB_Inc_Tree_280920_01");


  TLorentzVector *readbeam=NULL;
  TLorentzVector *readtarget=NULL;
  TLorentzVector Target(0,0,0,1.8756);

  vector<TLorentzVector> *v_p4=0;
  vector<TLorentzVector> *v_vertex=0;

  vector<double> *v_path=0;
  Double_t path;

  vector<double> *v_time=0;
  Double_t time;

  vector<double> *v_beta=0;

  Double_t start_time;
  vector<double> *v_energy=0;
  vector<double> *v_charge=0;
  vector<double> *v_PID=0;
  vector<double> *v_chi2PID=0;

  Int_t readchargetracks;
  Int_t readprotonno;
  Int_t readpipno;
  Int_t readpimno;
  Int_t readelno;
  Int_t readkaonpno;
  Int_t readothertracks;
  Int_t readeventno;
  Int_t readrunno;
  Int_t readtriggerno;

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
  t1->SetBranchAddress("time",&v_time);
  t1->SetBranchAddress("protonno",&readprotonno);
  t1->SetBranchAddress("pipno",&readpipno);
  t1->SetBranchAddress("pimno",&readpimno);
  t1->SetBranchAddress("elno",&readelno);
  t1->SetBranchAddress("kaonpno",&readkaonpno);
  t1->SetBranchAddress("eventno",&readeventno);
  t1->SetBranchAddress("runno",&readrunno);
  t1->SetBranchAddress("triggerno",&readtriggerno);

  TFile fileOutput1("/mnt/f/PhD/New_Data/Mass_RGB_Spring_Inbending_Inclusive_Pass1_290920_01.root","recreate");

  auto db=TDatabasePDG::Instance();


  //Creating histograms

  // Test histograms
  auto* helno=new TH1F("helno","Number of electrons",5,0,5);
  auto* hprotonno=new TH1F("hprotonno","Number of protons",5,0,5);
  auto* hkaonpno=new TH1F("hkaonpno","Number of K^{+}",5,0,5);
  auto* hpimno=new TH1F("hpimno","Number of #pi^{-}",5,0,5);

  auto* hmass_P=new TH1F("hmass_P","Calculated mass",700,-0.2,1.2);
  auto* hmass_N=new TH1F("hmass_N","Calculated mass",700,-0.2,1.2);

  auto* hbeta=new TH2D("hbeta","Beta measured", 500,0,11,400,-1,1);
  auto* hdelta_beta_el=new TH2D("hdelta_beta_el"," Delta Beta e^{-}", 500,0,11,400,-1,1);
  auto* hdelta_beta_pim=new TH2D("hdelta_beta_pim"," Delta Beta #pi^{-}", 500,0,11,400,-1,1);
  auto* hdelta_beta_km=new TH2D("hdelta_beta_km"," Delta Beta K^{-}", 500,0,11,400,-1,1);
  auto* hdelta_beta_kp=new TH2D("hdelta_beta_kp"," Delta Beta K^{+}", 500,0,11,400,-1,1);
  auto* hdelta_beta_pr=new TH2D("hdelta_beta_pr"," Delta Beta p", 500,0,11,400,-1,1);
  auto* hdelta_beta_pip=new TH2D("hdelta_beta_pip"," Delta Beta #pi^{+}", 500,0,11,400,-1,1);

  auto* hmiss_mass_all=new TH1F("miss_all","MM^2(e' K^{+} p #pi^{-});MM^2(e' K^{+} p #pi^{-}) [GeV];Counts",400,-1,3);
  auto* hmiss_momentum_all=new TH1F("hmiss_momentum_all","P(B + T - e' - K^{+} - p - #pi^{-});P(B + T - e' - K^{+} - p - #pi^{-}) [GeV];Counts",400,-1,3);
  auto* hmiss_1=new TH1F("miss_1","MM(e' K^{+});MM(e' K^{+}) [GeV];Counts",400,-1,3);
  auto* hmiss_1_b=new TH1F("hmiss_1_b","MM(e' K^{+});MM(e' K^{+}) [GeV];Counts",400,-1,3);
  auto* hmiss_1_2=new TH1F("hmiss_1_2","MM(e' K^{+});MM(e' K^{+}) [GeV];Counts",400,-1,3);
  auto* hmiss_3=new TH1F("hmiss_3","MM(e' K^{+} K^{+});MM(e' K^{+} K^{+}) [GeV];Counts",400,-1,3);
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
  auto* hmiss_mass_allc=new TH1F("miss_allc","MM^2(e' K^{+} p #pi^{-}) after cuts on delta beta and momentum of K^{+} and p;MM^2(e' K^{+} p #pi^{-}) [GeV];Counts",400,-1,3);
  auto* hmiss_momentum_allc=new TH1F("hmiss_momentum_allc","P(B + T - e' - K^{+} - p - #pi^{-}) after cuts on delta beta and momentum of K^{+} and p;P(B + T - e' - K^{+} - p - #pi^{-}) [GeV];Counts",400,-1,3);
  auto* hmiss_1c=new TH1F("miss_1c","MM(e' K^{+}) after cuts on delta beta and momentum of K^{+} and p;MM(e' K^{+}) [GeV];Counts",400,-1,3);
  auto* hmiss_3c=new TH1F("miss_3c","MM(e' K^{+} K^{+}) after cuts on delta beta and momentum of K^{+} and p;MM(e' K^{+} p) [GeV];Counts",400,-1,3);
  auto* hlambda_1c=new TH1F("lambda_1c","M(p #pi^{-}) after cuts on delta beta and momentum of K^{+} and p;M(p #pi^{-}) [GeV];Counts",750,0.5,2.0);

  //Exclusivity cuts as well
  auto* hmiss_1t=new TH1F("miss_1t","MM(e' K^{+}) after cuts on delta beta and momentum of K^{+} and p and MM^2(e' K^{+} p) cut;MM(e' K^{+}) [GeV];Counts",400,-1,3);
  auto* hmiss_3t=new TH1F("miss_3t","MM(e' K^{+} K^{+}) after cuts on delta beta and momentum of K^{+} and p and MM(e' K^{+} p) cut;MM^2(e' K^{+} p) [GeV];Counts",400,-1,3);
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

  vector<Double_t> v_beta_tof_el;
  Double_t beta_tof_el;
  vector<Double_t> v_P_el;
  Double_t P_el;
  vector<Double_t>v_path_el;
  Double_t path_el;
  vector<Double_t>v_TOF_el;
  Double_t TOF_el;
  vector<Double_t> v_beta_calc_el;
  Double_t beta_calc_el;
  vector<Double_t> v_delta_beta_el;
  Double_t delta_beta_el;
  vector<Double_t> v_vt_el;
  Double_t vt_el;
  vector<Double_t> v_vertex_time_el;
  Double_t vertex_time_el;


  vector<Double_t> v_beta_tof_pip;
  Double_t beta_tof_pip;
  vector<Double_t> v_P_pip;
  Double_t P_pip;
  vector<Double_t>v_path_pip;
  Double_t path_pip;
  vector<Double_t>v_TOF_pip;
  Double_t TOF_pip;
  vector<Double_t> v_beta_calc_pip;
  Double_t beta_calc_pip;
  vector<Double_t> v_delta_beta_pip;
  Double_t delta_beta_pip;
  vector<Double_t> v_vt_pip;
  Double_t vt_pip;
  vector<Double_t> v_vertex_time_pip;
  Double_t vertex_time_pip;

  vector<Double_t> v_beta_tof_pim;
  Double_t beta_tof_pim;
  vector<Double_t> v_P_pim;
  Double_t P_pim;
  vector<Double_t>v_path_pim;
  Double_t path_pim;
  vector<Double_t>v_TOF_pim;
  Double_t TOF_pim;
  vector<Double_t> v_beta_calc_pim;
  Double_t beta_calc_pim;
  vector<Double_t> v_delta_beta_pim;
  Double_t delta_beta_pim;
  vector<Double_t> v_vt_pim;
  Double_t vt_pim;
  vector<Double_t> v_vertex_time_pim;
  Double_t vertex_time_pim;

  vector<Double_t> v_beta_tof_pr;
  Double_t beta_tof_pr;
  vector<Double_t> v_P_pr;
  Double_t P_pr;
  vector<Double_t>v_path_pr;
  Double_t path_pr;
  vector<Double_t>v_TOF_pr;
  Double_t TOF_pr;
  vector<Double_t> v_beta_calc_pr;
  Double_t beta_calc_pr;
  vector<Double_t> v_delta_beta_pr;
  Double_t delta_beta_pr;
  vector<Double_t> v_vt_pr;
  Double_t vt_pr;
  vector<Double_t> v_vertex_time_pr;
  Double_t vertex_time_pr;

  vector<Double_t> v_beta_tof_kp;
  Double_t beta_tof_kp;
  vector<Double_t> v_P_kp;
  Double_t P_kp;
  vector<Double_t>v_path_kp;
  Double_t path_kp;
  vector<Double_t>v_TOF_kp;
  Double_t TOF_kp;
  vector<Double_t> v_beta_calc_kp;
  Double_t beta_calc_kp;
  vector<Double_t> v_delta_beta_kp;
  Double_t delta_beta_kp;
  vector<Double_t> v_vt_kp;
  Double_t vt_kp;
  vector<Double_t> v_time_kp;
  Double_t time_kp;
  vector<Double_t> v_vertex_time_kp;
  Double_t vertex_time_kp;

  vector<Double_t> v_beta_tof_kp_b;
  Double_t beta_tof_kp_b;
  vector<Double_t> v_P_kp_b;
  Double_t P_kp_b;
  vector<Double_t>v_path_kp_b;
  Double_t path_kp_b;
  vector<Double_t>v_TOF_kp_b;
  Double_t TOF_kp_b;
  vector<Double_t> v_beta_calc_kp_b;
  Double_t beta_calc_kp_b;
  vector<Double_t> v_delta_beta_kp_b;
  Double_t delta_beta_kp_b;
  vector<Double_t> v_vt_kp_b;
  Double_t vt_kp_b;
  vector<Double_t> v_time_kp_b;
  Double_t time_kp_b;
  vector<Double_t> v_vertex_time_kp_b;
  Double_t vertex_time_kp_b;

  vector<Double_t> v_beta_tof_km;
  Double_t beta_tof_km;
  vector<Double_t> v_P_km;
  Double_t P_km;
  vector<Double_t>v_path_km;
  Double_t path_km;
  vector<Double_t>v_TOF_km;
  Double_t TOF_km;
  vector<Double_t> v_beta_calc_km;
  Double_t beta_calc_km;
  vector<Double_t> v_delta_beta_km;
  Double_t delta_beta_km;
  vector<Double_t> v_vt_km;
  Double_t vt_km;
  vector<Double_t> v_vertex_time_km;
  Double_t vertex_time_km;

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


  // Long64_t nentries = 8000000;
  Long64_t nentries = t1->GetEntries();
  Int_t Percentage = nentries/100;
  for(Long64_t i=0; i<nentries;i++){
    t1->GetEntry(i);
    if (i % Percentage == 0){
      fprintf (stderr, "%lld\r", i/Percentage);
      fflush (stderr);
    }

    // Set the beam energy depending on the run number
    if(readrunno > 6400) continue;

    // readbeam.SetXYZM(0,0,10.6,10.6);
    // else readbeam.SetXYZM(0,0,10.2,10.2);

    v_el.clear();
    v_beta_tof_el.clear();
    v_P_el.clear();
    v_path_el.clear();
    v_TOF_el.clear();
    v_beta_calc_el.clear();
    v_delta_beta_el.clear();
    v_vertex_time_el.clear();
    v_vertex_el.clear();

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

    v_kp.clear();
    v_beta_tof_kp.clear();
    v_P_kp.clear();
    v_path_kp.clear();
    v_TOF_kp.clear();
    v_beta_calc_kp.clear();
    v_delta_beta_kp.clear();
    v_time_kp.clear();
    v_vertex_kp.clear();
    v_vertex_time_kp.clear();
    v_vertex_kp.clear();

    v_kp_b.clear();
    v_beta_tof_kp_b.clear();
    v_P_kp_b.clear();
    v_path_kp_b.clear();
    v_TOF_kp_b.clear();
    v_beta_calc_kp_b.clear();
    v_delta_beta_kp_b.clear();
    v_time_kp_b.clear();
    v_vertex_kp_b.clear();
    v_vertex_time_kp_b.clear();
    v_vertex_kp_b.clear();

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

    v_othertracks.clear();
    v_beta_tof_othertracks.clear();
    v_P_othertracks.clear();
    v_beta_calc_othertracks.clear();
    v_delta_beta_othertracks.clear();

    // cout<<"test1"<<endl;

    Int_t Nparticles = v_p4->size();
    for(Int_t j=0; j<Nparticles; j++){
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

        else if(Mass > 0.4 && Mass < 0.6){
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
        }

        // else if((Mass > 0.434225 && Mass < 0.453853) || (Mass > 0.532365 && Mass < 0.551993)){
        //   // cout<<"test6"<<endl;
        //   kp_b.SetXYZM(v_p4->at(j).Px(),v_p4->at(j).Py(),v_p4->at(j).Pz(),db->GetParticle(321)->Mass());
        //   P_kp_b = sqrt((pow(v_p4->at(j).Px(),2))+(pow(v_p4->at(j).Py(),2))+(pow(v_p4->at(j).Pz(),2))); // Measured momentum
        //   TOF_kp_b = v_time->at(j);
        //   path_kp_b = v_path->at(j);
        //   beta_tof_kp_b = v_beta->at(j);
        //   beta_calc_kp_b = P_kp_b/(sqrt((pow(P_kp_b,2))+(pow(db->GetParticle(321)->Mass(),2))));
        //   delta_beta_kp_b = beta_calc_kp_b-beta_tof_kp_b;
        //   vertex_time_kp_b = TOF_kp_b - path_kp_b / (beta_tof_kp_b*c);
        //   vertex_kp_b.SetXYZT(v_vertex->at(j).X(), v_vertex->at(j).Y(), v_vertex->at(j).Z(), vertex_time_kp_b);
        //
        //
        //   v_kp_b.push_back(kp_b);
        //   v_beta_tof_kp_b.push_back(beta_tof_kp_b);
        //   v_P_kp_b.push_back(P_kp_b);
        //   v_path_kp_b.push_back(path_kp_b);
        //   v_TOF_kp_b.push_back(TOF_kp_b);
        //   v_beta_calc_kp_b.push_back(beta_calc_kp_b);
        //   v_delta_beta_kp_b.push_back(delta_beta_kp_b);
        //   v_vertex_time_kp_b.push_back(vertex_time_kp_b);
        //   v_vertex_kp_b.push_back(vertex_kp_b);
        // }

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

      if(v_charge->at(j) < 0)hmass_N->Fill(Mass);
      if(v_charge->at(j) > 0)hmass_P->Fill(Mass);

    }

    helno->Fill(v_el.size());
    hpimno->Fill(v_pim.size());
    hkaonpno->Fill(v_kp.size());
    hprotonno->Fill(v_pr.size());



    // Strangeness 1 channels
    if(v_kp.size()==1 && v_el.size()==1 && v_pr.size()==1 && v_pim.size()==1){


      missall = (TLorentzVector)*readbeam + Target - v_el.at(0) - v_kp.at(0) - v_pr.at(0) - v_pim.at(0);
      miss1 = (TLorentzVector)*readbeam + Target - v_el.at(0) - v_kp.at(0);
      lambda = v_pr.at(0) + v_pim.at(0);

      hlambda_1->Fill(lambda.M());
      hmiss_1->Fill(miss1.M());
      hmiss_mass_all->Fill(missall.M2());
      hmiss_momentum_all->Fill(missall.Rho());

      if(fabs(v_delta_beta_kp.at(0))<0.01 && v_P_kp.at(0)>0.9 && v_P_kp.at(0)<3 && fabs(v_delta_beta_pr.at(0))<0.01 && v_P_pr.at(0)>0.2 && v_P_pr.at(0)<4.8 && fabs(v_delta_beta_pim.at(0))<0.01 && v_P_pim.at(0)>0.2 && v_P_pim.at(0)<3){
        hlambda_1c->Fill(lambda.M());
        hmiss_1c->Fill(miss1.M());
        hmiss_mass_allc->Fill(missall.M2());
        hmiss_momentum_allc->Fill(missall.Rho());

        if(fabs(missall.M2())<0.05 && fabs(missall.Rho())<0.1){
          hmiss_1t->Fill(miss1.M());
          hlambda_1t->Fill(lambda.M());
        }
      }
    }
//     if(v_kp_b.size()==1 && v_el.size()==1 && v_pr.size()==1 && v_pim.size()==1){
//       miss1_b = (TLorentzVector)*readbeam + Target - v_el.at(0) - v_kp_b.at(0);
//       hmiss_1_b->Fill(miss1_b.M());
//       if(fabs(v_delta_beta_kp_b.at(0))<0.01 && v_P_kp_b.at(0)>0.9 && v_P_kp_b.at(0)<3 && fabs(v_delta_beta_pr.at(0))<0.01 && v_P_pr.at(0)>0.2 && v_P_pr.at(0)<4.8 && fabs(v_delta_beta_pim.at(0))<0.01 && v_P_pim.at(0)>0.2 && v_P_pim.at(0)<3){
//         hmiss_1_2->Fill(miss1_b.M());
//
// }
//     }

    // // Strangeness 2 channels
    // if(v_kp.size()==2  && v_el.size()==1 /*&& v_pr.size()==1 && v_pim.size()==2*/){
    //   // missallcasc = (TLorentzVector)*readbeam + Target - v_el.at(0) - v_kp.at(0) - v_kp.at(1) - v_pr.at(0) - v_pim.at(0) - v_pim.at(1);
    //   miss3 = (TLorentzVector)*readbeam + Target - v_el.at(0) - v_kp.at(0) - v_kp.at(1);
    //   // cascade_lambda_1 = v_pr.at(0) + v_pim.at(0);
    //   // cascade_lambda_2 = v_pr.at(0) + v_pim.at(1);
    //   Delta_Vertex_Time = v_vertex_kp.at(0).T() - v_vertex_kp.at(1).T();
    //   // inv_cascade = v_pr.at(0) + v_pim.at(0) + v_pim.at(1);
    //
    //   hmiss_3->Fill(miss3.M());
    //
    //
    // }
    //
    // // Strangeness 3 channels with electron
    // if(v_kp.size()==3 && v_el.size()==1){
    //   miss5 = (TLorentzVector)*readbeam + Target - v_el.at(0) - v_kp.at(0) - v_kp.at(1) - v_kp.at(2);
    //   // Filling histograms with PID cuts only
    //   hmiss_5->Fill(miss5.M());
    //   // if(fabs(v_delta_beta_kp.at(0))<0.01 && v_P_kp.at(0)>0.5 && v_P_kp.at(0)<2.6 &&
    //   // fabs(v_delta_beta_kp.at(1))<0.01 && v_P_kp.at(1)>0.5 && v_P_kp.at(1)<2.6 &&
    //   // fabs(v_delta_beta_kp.at(2))<0.01 && v_P_kp.at(2)>0.5 && v_P_kp.at(2)<2.6){
    //   //
    //   //   // Filling histograms with PID, delta beta and momentum cuts
    //   //   hmiss_5c->Fill(miss5.M());
    //   // }
    // }
  }
  fileOutput1.Write();
}
