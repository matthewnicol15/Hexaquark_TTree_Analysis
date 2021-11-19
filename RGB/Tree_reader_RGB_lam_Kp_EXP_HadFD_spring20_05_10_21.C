#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TH1.h>
#include <TF1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <vector>
#include <TLorentzVector.h>
#include <math.h>
#include <TMath.h>
#include <TLine.h>

double Calc_dtfInterDOCA(const TVector3 &locUnitDir1, const TVector3 &locUnitDir2, const TVector3 &locVertex1, const TVector3 &locVertex2, TVector3 &locInterDOCA1, TVector3 &locInterDOCA2){
    //originated from code by Jörn Langheinrich
    //you can use this function to find the DOCA to a fixed point by calling this function with locUnitDir1 and 2 parallel, and the fixed vertex as locVertex2
    double locUnitDot = locUnitDir1*locUnitDir2;
    double locDenominator = locUnitDot*locUnitDot - 1.0; /// scalar product of directions
    double locDistVertToInterDOCA1 = 0.0, locDistVertToInterDOCA2 = 0.0; //distance from vertex to DOCA point


    if(fabs(locDenominator) < 1.0e-15) //parallel
        locDistVertToInterDOCA1 = (locVertex2 - locVertex1)*locUnitDir2/locUnitDot; //the opposite
    else{
        double locA = (locVertex1 - locVertex2)*locUnitDir1;
        double locB = (locVertex1 - locVertex2)*locUnitDir2;
        locDistVertToInterDOCA1 = (locA - locUnitDot*locB)/locDenominator;
        locDistVertToInterDOCA2 = (locUnitDot*locA - locB)/locDenominator;
    }


    locInterDOCA1 = locVertex1 + locDistVertToInterDOCA1*locUnitDir1; //intersection point of DOCA line and track 1
    locInterDOCA2 = locVertex2 + locDistVertToInterDOCA2*locUnitDir2; //intersection point of DOCA line and track 2
    double locDOCA = (locInterDOCA1 - locInterDOCA2).Mag();
    return ((locVertex2.Z() > locVertex1.Z()) ? locDOCA : -1.0*locDOCA);
}

void Tree_reader_RGB_lam_Kp_EXP_HadFD_spring20_05_10_21(){

  gROOT->ProcessLine(".L /media/mn688/Elements1/PhD/Macros/Loader.C+");
  TFile *f = new TFile("/media/mn688/Elements1/PhD/Trees/Dibaryon/RGB/RGB_kp_pim_p_sprig20_05_10_21.root");
  TTree *t1 = (TTree*)f->Get("RGB_kp_pim_p");

  vector<TLorentzVector> *v_p4=0, *v_neutral_p4=0;

  TLorentzVector *readbeam=NULL;
  TLorentzVector *readtarget=NULL;
  TLorentzVector proton, neutron;

  proton.SetXYZM(0, 0, 0, 0.9383);
  neutron.SetXYZM(0, 0, 0, 0.9396);

  vector<TLorentzVector> *v_vertex=0;

  vector<double> *v_beta=0;
  vector<double> *v_start_time=0;
  vector<double> *v_TOF=0;
  vector<double> *v_path=0;
  vector<int> *v_region=0;
  vector<double> *v_status=0;

  vector<double> *v_energy=0, *v_chi2PID=0;

  vector<int> *v_chargetracks=0, *v_protonno=0, *v_kaonpno=0, *v_kaonmno=0, *v_kpno=0, *v_pipno=0, *v_pimno=0, *v_gamno=0, *v_elno=0, *v_other=0;

  vector<int> *v_event_num=0, *v_charge=0, *v_PID=0;

  Int_t PID;

  t1->SetBranchAddress("p4",&v_p4);

  t1->SetBranchAddress("vertex",&v_vertex);

  t1->SetBranchAddress("beta",&v_beta);
  t1->SetBranchAddress("start_time",&v_start_time);
  t1->SetBranchAddress("TOF",&v_TOF);
  t1->SetBranchAddress("path",&v_path);
  t1->SetBranchAddress("region",&v_region);
  t1->SetBranchAddress("status",&v_status);

  t1->SetBranchAddress("beam",&readbeam);
  t1->SetBranchAddress("target",&readtarget);

  t1->SetBranchAddress("energy",&v_energy);
  t1->SetBranchAddress("charge",&v_charge);
  t1->SetBranchAddress("PID",&v_PID);
  t1->SetBranchAddress("chi2PID",&v_chi2PID);
  t1->SetBranchAddress("chargetracks",&v_chargetracks);
  t1->SetBranchAddress("protonno",&v_protonno);
  t1->SetBranchAddress("kaonpno",&v_kaonpno);
  t1->SetBranchAddress("pipno",&v_pipno);
  t1->SetBranchAddress("pimno",&v_pimno);
  t1->SetBranchAddress("gamno",&v_gamno);
  t1->SetBranchAddress("elno",&v_elno);

  t1->SetBranchAddress("EventNum",&v_event_num);

  TFile fileOutput1("/media/mn688/Elements1/PhD/Analysis_Output/RGB/Inclusive/Inbending/Strangeness_1/PID/Strangeness_1_RGB_SPRING_2020_Inbending_e_Kp_FD_131020_01.root","recreate");

  auto* hbeta_kp1=new TH2F("hbeta_kp1","#Delta#beta for K^{+} (post PID);P [GeV/c];#Delta#beta",600,0,6,200,-0.02,0.02);
  hbeta_kp1->GetXaxis()->SetLabelSize(0.05);
  hbeta_kp1->GetYaxis()->SetLabelSize(0.05);

  auto* hbeta_pim1=new TH2F("hbeta_pim1","#Delta#beta for #pi^{-} (post PID);P [GeV/c];#Delta#beta",600,0,6,200,-0.02,0.02);
  hbeta_pim1->GetXaxis()->SetLabelSize(0.05);
  hbeta_pim1->GetYaxis()->SetLabelSize(0.05);

  auto* hbeta_pr1=new TH2F("hbeta_pr1","#Delta#beta for p (post PID);P [GeV/c];#Delta#beta",600,0,6,200,-0.02,0.02);
  hbeta_pr1->GetXaxis()->SetLabelSize(0.05);
  hbeta_pr1->GetYaxis()->SetLabelSize(0.05);

  auto* hMM_all = new TH1F("hMM_all","Missing Mass of All Detected Particles;MM(e' K^{+} #pi^{-} p) [GeV]",500,-1,1);
  hMM_all->GetXaxis()->SetLabelSize(0.05);
  hMM_all->GetYaxis()->SetLabelSize(0.05);

  auto* hMM_all_LamCut = new TH1F("hMM_all_LamCut","Missing Mass of All Detected Particles;MM(e' K^{+} #pi^{-} p) [GeV]",200,-1,1);
  hMM_all_LamCut->GetXaxis()->SetLabelSize(0.05);
  hMM_all_LamCut->GetYaxis()->SetLabelSize(0.05);

  auto* hMM_all_LamCut_LamMass = new TH1F("hMM_all_LamCut_LamMass","Missing Mass of All Detected Particles;MM(e' K^{+} #pi^{-} p) [GeV]",200,-1,1);
  hMM_all_LamCut_LamMass->GetXaxis()->SetLabelSize(0.05);
  hMM_all_LamCut_LamMass->GetYaxis()->SetLabelSize(0.05);

  auto* hMM_all_momCut = new TH1F("hMM_all_momCut","Missing Mass of All Detected Particles;MM(e' K^{+} #pi^{-} p) [GeV]",500,-1,1);
  hMM_all_momCut->GetXaxis()->SetLabelSize(0.05);
  hMM_all_momCut->GetYaxis()->SetLabelSize(0.05);

  auto* hMM_all_LamCut_momCut = new TH1F("hMM_all_LamCut_momCut","Missing Mass of All Detected Particles;MM(e' K^{+} #pi^{-} p) [GeV]",200,-1,1);
  hMM_all_LamCut_momCut->GetXaxis()->SetLabelSize(0.05);
  hMM_all_LamCut_momCut->GetYaxis()->SetLabelSize(0.05);

  auto* hMM_all_LamCut_momCut_LamMass = new TH1F("hMM_all_LamCut_momCut_LamMass","Missing Mass of All Detected Particles;MM(e' K^{+} #pi^{-} p) [GeV]",200,-1,1);
  hMM_all_LamCut_momCut_LamMass->GetXaxis()->SetLabelSize(0.05);
  hMM_all_LamCut_momCut_LamMass->GetYaxis()->SetLabelSize(0.05);

  auto* hMrho_all = new TH1F("hMrho_all","Missing Momentum of All Detected Particles;ME(e' K^{+} #pi^{-} p) [GeV]",500,0,1);
  hMrho_all->GetXaxis()->SetLabelSize(0.05);
  hMrho_all->GetYaxis()->SetLabelSize(0.05);

  auto* hMrho_all_nCut = new TH1F("hMrho_all_nCut","Missing Momentum of All Detected Particles;ME(e' K^{+} #pi^{-} p) [GeV]",100,0,1);
  hMrho_all_nCut->GetXaxis()->SetLabelSize(0.05);
  hMrho_all_nCut->GetYaxis()->SetLabelSize(0.05);

  auto* hMrho_all_LamCut = new TH1F("hMrho_all_LamCut","Missing Momentum of All Detected Particles;ME(e' K^{+} #pi^{-} p) [GeV]",100,0,1);
  hMrho_all_LamCut->GetXaxis()->SetLabelSize(0.05);
  hMrho_all_LamCut->GetYaxis()->SetLabelSize(0.05);

  auto* hMrho_all_nCut_LamCut = new TH1F("hMrho_all_nCut_LamCut","Missing Momentum of All Detected Particles;ME(e' K^{+} #pi^{-} p) [GeV]",100,0,1);
  hMrho_all_nCut_LamCut->GetXaxis()->SetLabelSize(0.05);
  hMrho_all_nCut_LamCut->GetYaxis()->SetLabelSize(0.05);

  auto* hInv_lam = new TH1F("hInv_lam","Invariant Mass of p and #pi^{-};M(p #pi^{-}) [GeV]",250,1,1.5);
  hInv_lam->GetXaxis()->SetLabelSize(0.05);
  hInv_lam->GetYaxis()->SetLabelSize(0.05);

  auto* hInv_lam_nCut = new TH1F("hInv_lam_nCut","Invariant Mass of p and #pi^{-};M(p #pi^{-}) [GeV]",250,1,1.5);
  hInv_lam_nCut->GetXaxis()->SetLabelSize(0.05);
  hInv_lam_nCut->GetYaxis()->SetLabelSize(0.05);

  auto* hInv_lam_nCut_momCut = new TH1F("hInv_lam_nCut_momCut","Invariant Mass of p and #pi^{-};M(p #pi^{-}) [GeV]",250,1,1.5);
  hInv_lam_nCut_momCut->GetXaxis()->SetLabelSize(0.05);
  hInv_lam_nCut_momCut->GetYaxis()->SetLabelSize(0.05);

  auto* hInvLamVS_MMall=new TH2F("hInvLamVS_MMall","Invariant #Lambda Vs Missing n;M [GeV];MM [GeV]",250,0.5,1.5,125,1,1.5);
  hInvLamVS_MMall->GetXaxis()->SetLabelSize(0.05);
  hInvLamVS_MMall->GetYaxis()->SetLabelSize(0.05);

  auto* hMMds = new TH1F("hMMds","Missing Mass of the e' and K^{+};MM(e' K^{+}) [GeV]",250,1.5,3.5);
  hMMds->GetXaxis()->SetLabelSize(0.05);
  hMMds->GetYaxis()->SetLabelSize(0.05);

  //TLorentzVectors
  vector <TLorentzVector> v_el, v_kp, v_pim, v_pr;
  TLorentzVector el, el1, kp, kp1, kp2, pim, pim1, pim2, pim3, pr, pr1, pr2;

  vector <int> v_region_el, v_region_pr, v_region_kp, v_region_pim;
  Int_t region_el, region_el1, region_pr, region_pr1, region_pr2, region_kp, region_kp1, region_kp2, region_pim, region_pim1, region_pim2, region_pim3;

  vector <TLorentzVector> v_ver_el, v_ver_kp, v_ver_pim, v_ver_pr;
  TLorentzVector ver_el, ver_el1, ver_kp, ver_kp1, ver_kp2, ver_pim, ver_pim1, ver_pim2, ver_pim3, ver_pr, ver_pr1, ver_pr2;

  TLorentzVector miss_all, miss_all_lamMass, miss_Lam, virtual_photon, Inv_lam, Inv_lam_PDGm, pr1_Lam, pim1_Lam, MMds, beam(0, 0, 10.4, 10.4);

  TVector3 Lam_boost;

  Double_t P_ei, P_ef, q_sqr;
  //TLorentzVectors

  //Beta
  vector <double> v_beta_tof_el, v_beta_calc_el;
  Double_t beta_tof_el, beta_tof_el1, beta_calc_el, beta_calc_el1;

  vector <double> v_beta_tof_kp, v_beta_calc_kp;
  Double_t beta_tof_kp, beta_tof_kp1, beta_tof_kp2, beta_calc_kp, beta_calc_kp1, delta_beta_kp1, beta_calc_kp2, delta_beta_kp2;

  vector <double> v_beta_tof_pim, v_beta_calc_pim;
  Double_t beta_tof_pim, beta_tof_pim1, beta_tof_pim2, beta_tof_pim3, beta_calc_pim, beta_calc_pim1, delta_beta_pim1, beta_calc_pim2, delta_beta_pim2, beta_calc_pim3, delta_beta_pim3;

  vector <double> v_beta_tof_pr, v_beta_calc_pr;
  Double_t beta_tof_pr, beta_tof_pr1, beta_tof_pr2, beta_calc_pr, beta_calc_pr1, delta_beta_pr1, beta_calc_pr2, delta_beta_pr2;
  //Beta

  //Flight time and path
  vector <double> v_TOF_el, v_TOF_pr, v_TOF_kp, v_TOF_pim, v_TOF_gam;
  Double_t TOF_el, TOF_el1, TOF_pr, TOF_pr1, TOF_pr2, TOF_kp, TOF_kp1, TOF_kp2, TOF_pim, TOF_pim1, TOF_pim2, TOF_pim3;

  vector <double> v_path_el, v_path_pr, v_path_kp, v_path_pim;
  Double_t path_el, path_el1, path_pr, path_pr1, path_pr2, path_kp, path_kp1, path_kp2, path_pim, path_pim1, path_pim2, path_pim3;

  Double_t c = 29.9792458;
  //Flight time and path

  //Vertex time
  vector <double> v_ver_time_el, v_ver_time_pr, v_ver_time_kp, v_ver_time_pim, v_ver_time_gam;
  Double_t ver_time_el, ver_time_el1, ver_time_pr, ver_time_pr1, ver_time_pr2, ver_time_kp, ver_time_kp1, ver_time_kp2, ver_time_pim, ver_time_pim1, ver_time_pim2, ver_time_pim3;

  Double_t diff_ver_el_pr, diff_ver_el_kp, diff_ver_el_pim, diff_ver_pr_kp, diff_ver_pr_pim, diff_ver_kp_pim;

  Double_t diff_ver_time_el_pr, diff_ver_time_el_kp, diff_ver_time_el_pim, diff_ver_time_pr_kp, diff_ver_time_pr_pim, diff_ver_time_kp_pim;
  //Vertex time

  //chi2PID
  vector <double> v_chi2_kp;
  Double_t chi2_kp, chi2_kp1;

  vector <double> v_chi2_pr;
  Double_t chi2_pr, chi2_pr1;

  vector <double> v_chi2_pim;
  Double_t chi2_pim, chi2_pim1;
  //chi2PID

  //status
  vector <double> v_status_pim;
  Double_t status_pim, status_pim1;

  vector <double> v_status_kp;
  Double_t status_kp, status_kp1;

  vector <double> v_status_pr;
  Double_t status_pr, status_pr1;
  //status

  //DOCA
  double doca_K0, doca_p_kp1, doca_p_pim1;
  TVector3 old_vert_kp1, old_vert_pim1, old_vert_pr1, old_dir_kp1, old_dir_pim1, old_dir_pr1, new_vert_kp1, new_vert_pim1, new_vert_pr1;
  TVector3 loc_DOCA_vertex;
  //DOCA

  //Angles
  Double_t Lam_open_angle, photon_open_angle;
  //Angles

  Long64_t nentries = t1->GetEntries();
  // Long64_t nentries = 1000000;
  Int_t Percentage = nentries/100;
  cout<<nentries<<endl;
  Int_t Np;

  //cout<<nentries<<endl;

  for(Long64_t i=0; i<nentries;i++){
    t1->GetEntry(i);

    //cout<<v_event_num->at(i)<<endl;

    if (i % Percentage == 0){
      fprintf (stderr, "%lld\r", i/Percentage);
      fflush (stderr);
    }

    Np = v_p4->size();

    v_el.clear();
    v_beta_tof_el.clear();
    v_beta_calc_el.clear();
    v_TOF_el.clear();
    v_path_el.clear();
    v_ver_time_el.clear();
    v_ver_el.clear();
    v_region_el.clear();

    v_kp.clear();
    v_beta_tof_kp.clear();
    v_beta_calc_kp.clear();
    v_TOF_kp.clear();
    v_path_kp.clear();
    v_ver_time_kp.clear();
    v_ver_kp.clear();
    v_region_kp.clear();
    v_chi2_kp.clear();
    v_status_kp.clear();

    v_pim.clear();
    v_beta_tof_pim.clear();
    v_beta_calc_pim.clear();
    v_TOF_pim.clear();
    v_path_pim.clear();
    v_ver_time_pim.clear();
    v_ver_pim.clear();
    v_region_pim.clear();
    v_chi2_pim.clear();
    v_status_pim.clear();

    v_pr.clear();
    v_beta_tof_pr.clear();
    v_beta_calc_pr.clear();
    v_TOF_pr.clear();
    v_path_pr.clear();
    v_ver_time_pr.clear();
    v_ver_pr.clear();
    v_region_pr.clear();
    v_chi2_pr.clear();
    v_status_pr.clear();

    for(Int_t j=0; j<Np; j++){

      PID = v_PID->at(j);

      if(v_p4->at(j).M()>0.0001){

        if(PID==11){
          el.SetXYZM(v_p4->at(j).Px(),v_p4->at(j).Py(),v_p4->at(j).Pz(),v_p4->at(j).M());

          beta_tof_el = v_beta->at(j);
          beta_calc_el = el.Rho()/el.E();
          TOF_el = v_TOF->at(j);
          path_el = v_path->at(j);
          ver_time_el = TOF_el - path_el/(beta_calc_el*c);

          ver_el.SetXYZT(v_vertex->at(j).X(),v_vertex->at(j).Y(),v_vertex->at(j).Z(),ver_time_el);

          region_el = v_region->at(j);

          v_el.push_back(el);
          v_beta_tof_el.push_back(beta_tof_el);
          v_beta_calc_el.push_back(beta_calc_el);
          v_TOF_el.push_back(TOF_el);
          v_path_el.push_back(path_el);
          v_ver_time_el.push_back(ver_time_el);
          v_ver_el.push_back(ver_el);
          v_region_el.push_back(region_el);
        }
        else if(PID==321){
          kp.SetXYZM(v_p4->at(j).Px(),v_p4->at(j).Py(),v_p4->at(j).Pz(),v_p4->at(j).M());

          beta_tof_kp = v_beta->at(j);
          beta_calc_kp = kp.Rho()/kp.E();

          TOF_kp = v_TOF->at(j);
          path_kp = v_path->at(j);
          ver_time_kp = TOF_kp - path_kp/(beta_tof_kp*c);

          ver_kp.SetXYZT(v_vertex->at(j).X(),v_vertex->at(j).Y(),v_vertex->at(j).Z(),ver_time_kp);

          region_kp = v_region->at(j);

          chi2_kp = v_chi2PID->at(j);

          status_kp = v_status->at(j);

          v_kp.push_back(kp);
          v_beta_tof_kp.push_back(beta_tof_kp);
          v_beta_calc_kp.push_back(beta_calc_kp);
          v_TOF_kp.push_back(TOF_kp);
          v_path_kp.push_back(path_kp);
          v_ver_time_kp.push_back(ver_time_kp);
          v_ver_kp.push_back(ver_kp);
          v_region_kp.push_back(region_kp);
          v_chi2_kp.push_back(chi2_kp);
          v_status_kp.push_back(status_kp);
        }
        else if(PID==-211){
          pim.SetXYZM(v_p4->at(j).Px(),v_p4->at(j).Py(),v_p4->at(j).Pz(),v_p4->at(j).M());

          beta_tof_pim = v_beta->at(j);
          beta_calc_pim = pim.Rho()/pim.E();

          TOF_pim = v_TOF->at(j);
          path_pim = v_path->at(j);
          ver_time_pim = TOF_pim - path_pim/(beta_calc_pim*c);

          ver_pim.SetXYZT(v_vertex->at(j).X(),v_vertex->at(j).Y(),v_vertex->at(j).Z(),ver_time_pim);

          region_pim = v_region->at(j);

          chi2_pim = v_chi2PID->at(j);

          status_pim = v_status->at(j);

          v_pim.push_back(pim);
          v_beta_tof_pim.push_back(beta_tof_pim);
          v_beta_calc_pim.push_back(beta_calc_pim);
          v_TOF_pim.push_back(TOF_pim);
          v_path_pim.push_back(path_pim);
          v_ver_time_pim.push_back(ver_time_pim);
          v_ver_pim.push_back(ver_pim);
          v_region_pim.push_back(region_pim);
          v_chi2_pim.push_back(chi2_pim);
          v_status_pim.push_back(status_pim);

        }
        else if(PID==2212){
          pr.SetXYZM(v_p4->at(j).Px(),v_p4->at(j).Py(),v_p4->at(j).Pz(),v_p4->at(j).M());

          beta_tof_pr = v_beta->at(j);
          beta_calc_pr = pr.Rho()/pr.E();

          TOF_pr = v_TOF->at(j);
          path_pr = v_path->at(j);
          ver_time_pr = TOF_pr - path_pr/(beta_calc_pr*c);

          ver_pr.SetXYZT(v_vertex->at(j).X(),v_vertex->at(j).Y(),v_vertex->at(j).Z(),ver_time_pr);

          region_pr = v_region->at(j);

          chi2_pr = v_chi2PID->at(j);

          status_pr = v_status->at(j);

          v_pr.push_back(pr);
          v_beta_tof_pr.push_back(beta_tof_pr);
          v_beta_calc_pr.push_back(beta_calc_pr);
          v_TOF_pr.push_back(TOF_pr);
          v_path_pr.push_back(path_pr);
          v_ver_time_pr.push_back(ver_time_pr);
          v_ver_pr.push_back(ver_pr);
          v_region_pr.push_back(region_pr);
          v_chi2_pr.push_back(chi2_pr);
          v_status_pr.push_back(status_pr);

        }
      }

    }

    if(v_el.size()<1 || v_kp.size()<1 || v_pim.size()<1 || v_pr.size()<1) continue;

    // for(Int_t i = 0; i < v_el.size(); i++){
    //
    // }

    el1 = v_el.at(0);
    beta_tof_el1 = v_beta_tof_el.at(0);
    beta_calc_el1 = v_beta_calc_el.at(0);
    TOF_el1 = v_TOF_el.at(0);
    path_el1 = v_path_el.at(0);
    ver_time_el1 = v_ver_time_el.at(0);
    ver_el1 = v_ver_el.at(0);
    region_el1 = v_region_el.at(0);

    kp1 = v_kp.at(0);
    beta_tof_kp1 = v_beta_tof_kp.at(0);
    beta_calc_kp1 = v_beta_calc_kp.at(0);
    TOF_kp1 = v_TOF_kp.at(0);
    path_kp1 = v_path_kp.at(0);
    ver_time_kp1 = v_ver_time_kp.at(0);
    ver_kp1 = v_ver_kp.at(0);
    region_kp1 = v_region_kp.at(0);
    chi2_kp1 = v_chi2_kp.at(0);
    status_kp1 = v_status_kp.at(0);

    pim1 = v_pim.at(0);
    beta_tof_pim1 = v_beta_tof_pim.at(0);
    beta_calc_pim1 = v_beta_calc_pim.at(0);
    TOF_pim1 = v_TOF_pim.at(0);
    path_pim1 = v_path_pim.at(0);
    ver_time_pim1 = v_ver_time_pim.at(0);
    ver_pim1 = v_ver_pim.at(0);
    region_pim1 = v_region_pim.at(0);
    chi2_pim1 = v_chi2_pim.at(0);
    status_pim1 = v_status_pim.at(0);

    pr1 = v_pr.at(0);
    beta_tof_pr1 = v_beta_tof_pr.at(0);
    beta_calc_pr1 = v_beta_calc_pr.at(0);
    TOF_pr1 = v_TOF_pr.at(0);
    path_pr1 = v_path_pr.at(0);
    ver_time_pr1 = v_ver_time_pr.at(0);
    ver_pr1 = v_ver_pr.at(0);
    region_pr1 = v_region_pr.at(0);
    chi2_pr1 = v_chi2_pr.at(0);
    status_pr1 = v_status_pr.at(0);

    if(region_el1 != 1 || region_kp1 != 1 || region_pr1 != 1 || region_pim1 != 1) continue;

      delta_beta_kp1 = beta_calc_kp1-beta_tof_kp1;

      delta_beta_pim1 = beta_calc_pim1-beta_tof_pim1;

      delta_beta_pr1 = beta_calc_pr1-beta_tof_pr1;

    if(abs(delta_beta_kp1) > 0.01 || abs(delta_beta_pim1) > 0.01 || abs(delta_beta_pr1) > 0.01 || (kp1.Rho() > 0.65 && kp1.Rho() < 0.83)) continue;

    hbeta_kp1->Fill(kp1.Rho(),delta_beta_kp1);

    hbeta_pim1->Fill(pim1.Rho(),delta_beta_pim1);

    hbeta_pr1->Fill(pr1.Rho(),delta_beta_pr1);

    Inv_lam = pr1 + pim1;

    miss_all = (TLorentzVector)*readbeam + proton - el1 - Inv_lam - kp1;

    virtual_photon = (TLorentzVector)*readbeam - el1;

    q_sqr = -virtual_photon.M2();

    hMM_all->Fill(miss_all.M2());

    hInv_lam->Fill(Inv_lam.M());

    hInvLamVS_MMall->Fill(miss_all.M2(),Inv_lam.M());



    hMrho_all->Fill(miss_all.Rho());

    //MMn cut
    if(fabs(miss_all.M2()) < 0.02){
      hMrho_all_nCut->Fill(miss_all.Rho());

      hInv_lam_nCut->Fill(Inv_lam.M());

      //hLam_openAng->Fill(Lam_open_angle);
    }
    //MMn cut

    //InvLamCut
    if(Inv_lam.M() < 1.14){
      Inv_lam_PDGm.SetXYZM(Inv_lam.Px(), Inv_lam.Py(), Inv_lam.Pz(), 1.1157);

      miss_all_lamMass = (TLorentzVector)*readbeam + proton - el1 - Inv_lam_PDGm - kp1;

      hMM_all_LamCut->Fill(miss_all.M2());

      hMM_all_LamCut_LamMass->Fill(miss_all_lamMass.M());

      hMrho_all_LamCut->Fill(miss_all.Rho());
    }
    //InvLamCut

    //MMn cut + InvLamCut
    if(fabs(miss_all.M2()) < 0.02 && Inv_lam.M() < 1.14){
      hMrho_all_nCut_LamCut->Fill(miss_all.Rho());
    }
    //MMn cut + InvLamCut

    //InvLamCut + momCut
    if(Inv_lam.M() < 1.14 && miss_all.Rho() > 0.2){
      hMM_all_LamCut_momCut->Fill(miss_all.M2());
      hMM_all_LamCut_momCut_LamMass->Fill(miss_all_lamMass.M());
    }
    //InvLamCut + momCut

    //MMn cut + momCut
    if(fabs(miss_all.M2()) < 0.02 && miss_all.Rho() > 0.2){
      hInv_lam_nCut_momCut->Fill(Inv_lam.M());
    }
    //MMn cut + momCut

    MMds = (TLorentzVector)*readbeam + proton - el1 - kp1;

    if(fabs(miss_all.M2()) < 0.02 && Inv_lam.M() < 1.14 && miss_all.Rho() > 0.2){

      hMMds->Fill(MMds.M());

    }


  }

  fileOutput1.Write();

  f->Close();

}
