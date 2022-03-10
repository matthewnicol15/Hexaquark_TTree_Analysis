{

TFile *f1 = new TFile("/media/mn688/Elements1/PhD/Analysis_Output/Hexaquark/RGA_Fall2018_Outbending_at_least_1e1KpFD_Tree_Total_24022022_Total_Scaling_09032022_02.root");
TH3F* hist=(TH3F*)f1->Get("h_S1_Kaon_Momentum__Miss_Mass__Kaon_Mass");

TH1F* Missing_Mass_1 = (TH1F*) hist->ProjectionY("",hist->GetXaxis()->FindBin(0.5),hist->GetXaxis()->FindBin(1.0),0,hist->GetNbinsZ())->Clone("Missing_Mass_1");
TH1F* Missing_Mass_2 = (TH1F*) hist->ProjectionY("",hist->GetXaxis()->FindBin(1.0),hist->GetXaxis()->FindBin(1.5),0,hist->GetNbinsZ())->Clone("Missing_Mass_2");
TH1F* Missing_Mass_3 = (TH1F*) hist->ProjectionY("",hist->GetXaxis()->FindBin(1.5),hist->GetXaxis()->FindBin(2.0),0,hist->GetNbinsZ())->Clone("Missing_Mass_3");
TH1F* Missing_Mass_4 = (TH1F*) hist->ProjectionY("",hist->GetXaxis()->FindBin(2.0),hist->GetXaxis()->FindBin(2.5),0,hist->GetNbinsZ())->Clone("Missing_Mass_4");
TH1F* Missing_Mass_5 = (TH1F*) hist->ProjectionY("",hist->GetXaxis()->FindBin(2.5),hist->GetXaxis()->FindBin(3.0),0,hist->GetNbinsZ())->Clone("Missing_Mass_5");

TH1F* Kaon_Mass_1 = (TH1F*) hist->ProjectionZ("",hist->GetXaxis()->FindBin(0.5),hist->GetXaxis()->FindBin(1.0),0,hist->GetNbinsY())->Clone("Kaon_Mass_1");
TH1F* Kaon_Mass_2 = (TH1F*) hist->ProjectionZ("",hist->GetXaxis()->FindBin(1.0),hist->GetXaxis()->FindBin(1.5),0,hist->GetNbinsY())->Clone("Kaon_Mass_2");
TH1F* Kaon_Mass_3 = (TH1F*) hist->ProjectionZ("",hist->GetXaxis()->FindBin(1.5),hist->GetXaxis()->FindBin(2.0),0,hist->GetNbinsY())->Clone("Kaon_Mass_3");
TH1F* Kaon_Mass_4 = (TH1F*) hist->ProjectionZ("",hist->GetXaxis()->FindBin(2.0),hist->GetXaxis()->FindBin(2.5),0,hist->GetNbinsY())->Clone("Kaon_Mass_4");
TH1F* Kaon_Mass_5 = (TH1F*) hist->ProjectionZ("",hist->GetXaxis()->FindBin(2.5),hist->GetXaxis()->FindBin(3.0),0,hist->GetNbinsY())->Clone("Kaon_Mass_5");

Double_t Missing_Mass_1_integral = Missing_Mass_1->Integral();
Double_t Missing_Mass_2_integral = Missing_Mass_2->Integral();
Double_t Missing_Mass_3_integral = Missing_Mass_3->Integral();
Double_t Missing_Mass_4_integral = Missing_Mass_4->Integral();
Double_t Missing_Mass_5_integral = Missing_Mass_5->Integral();

Double_t Kaon_Mass_1_integral = Kaon_Mass_1->Integral();
Double_t Kaon_Mass_2_integral = Kaon_Mass_2->Integral();
Double_t Kaon_Mass_3_integral = Kaon_Mass_3->Integral();
Double_t Kaon_Mass_4_integral = Kaon_Mass_4->Integral();
Double_t Kaon_Mass_5_integral = Kaon_Mass_5->Integral();

Missing_Mass_2->Scale(Missing_Mass_1_integral / Missing_Mass_2_integral);
Missing_Mass_3->Scale(Missing_Mass_1_integral / Missing_Mass_3_integral);
Missing_Mass_4->Scale(Missing_Mass_1_integral / Missing_Mass_4_integral);
Missing_Mass_5->Scale(Missing_Mass_1_integral / Missing_Mass_5_integral);


Kaon_Mass_2->Scale(Kaon_Mass_1_integral / Kaon_Mass_2_integral);
Kaon_Mass_3->Scale(Kaon_Mass_1_integral / Kaon_Mass_3_integral);
Kaon_Mass_4->Scale(Kaon_Mass_1_integral / Kaon_Mass_4_integral);
Kaon_Mass_5->Scale(Kaon_Mass_1_integral / Kaon_Mass_5_integral);


Missing_Mass_2->SetLineColor(kRed);
Missing_Mass_3->SetLineColor(kOrange);
Missing_Mass_4->SetLineColor(kBlack);
Missing_Mass_5->SetLineColor(kGreen);

Kaon_Mass_2->SetLineColor(kRed);
Kaon_Mass_3->SetLineColor(kOrange);
Kaon_Mass_4->SetLineColor(kBlack);
Kaon_Mass_5->SetLineColor(kGreen);


TF1 *func1 = new TF1("func1","gaus(0)+gaus(3)+pol3(6)",0.363,0.6);
TF1 *func2 = new TF1("func2","gaus(0)+gaus(3)",0.363,0.6);
TF1 *func3 = new TF1("func3","pol3(0)",0.363,0.6);
func1->SetParameter(0,Kaon_Mass_3->GetMaximum()/2);
func1->SetParameter(1,0.493);
func1->SetParameter(2,0.03);
func1->SetParameter(3,Kaon_Mass_3->GetMaximum()/3);
func1->SetParameter(4,0.493);
func1->SetParameter(5,0.05);
// func1->SetParameter(5,0.05);
// func1->SetParameter(5,0.05);
// func1->SetParameter(5,0.05);
// func1->SetParameter(5,0.05);

Kaon_Mass_3->Fit("func1","RB");

func2->FixParameter(0,func1->GetParameter(0));
func2->FixParameter(1,func1->GetParameter(1));
func2->FixParameter(2,func1->GetParameter(2));
func2->FixParameter(3,func1->GetParameter(3));
func2->FixParameter(4,func1->GetParameter(4));
func2->FixParameter(5,func1->GetParameter(5));
func3->FixParameter(0,func1->GetParameter(6));
func3->FixParameter(1,func1->GetParameter(7));
func3->FixParameter(2,func1->GetParameter(8));
func3->FixParameter(3,func1->GetParameter(9));
// func3->FixParameter(4,func1->GetParameter(10));

func1->SetLineColor(kBlue);
func2->SetLineColor(kBlack);
func3->SetLineColor(kGreen);



auto *c1 = new TCanvas("c1","",800,800);
c1->cd();
Missing_Mass_1->Draw("same");
Missing_Mass_2->Draw("same,hist");
Missing_Mass_3->Draw("same,hist");
Missing_Mass_4->Draw("same,hist");
Missing_Mass_5->Draw("same,hist");

auto *c2=new TCanvas("c2","",800,800);
c2->cd();
Kaon_Mass_1->Draw("same");
Kaon_Mass_2->Draw("same,hist");
Kaon_Mass_3->Draw("same,hist");
Kaon_Mass_4->Draw("same,hist");
Kaon_Mass_5->Draw("same,hist");

auto *c3 = new TCanvas("c3","",800,800);
c3->cd();
Kaon_Mass_3->Draw("same,hist");
func1->Draw("same");
func2->Draw("same");
func3->Draw("same");

}
