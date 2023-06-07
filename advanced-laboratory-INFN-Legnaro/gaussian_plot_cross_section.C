#include "TH1.h"
#include "TMath.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"

{
ofstream RISULTATI;
	RISULTATI.open("RISULTATI.txt");

TF1 *g1 = new TF1("g1","[0]*TMath::Gaus(x,[1],[2])+[3]*TMath::Gaus(x,[4],[5])", -0.3, 0.3 );
Double_t p0=5000;
Double_t p1=0;
Double_t p2=0.0392699082;
Double_t p3=1200;
Double_t p4=0;
Double_t p5=0.0872664626;
g1->SetParameters(p0,p1,p2,p3,p4,p5);
/*
g1->SetParLimit(0, 4000, 9000);
g1->FixParameter(1, 0);
g1->SetParLimit(2, 0.03, 0.06);
g1->SetParLimit(3, 200, 1000);
g1->FixParameter(4,0);
g1->SetParLimit(5, 0.07, 0.095);
*/
TCanvas *c1 = new TCanvas("c1");

TGraph *data47=new TGraph("47rad.dat");
data47->GetYaxis()->SetRangeUser(0,8000);
data47->Fit("g1", "", "", -0.201, 0.201);
data47->Draw("");

data47->SetTitle ("47 MeV");
c1->SaveAs("47MeV_grad.pdf");

TF1 *fit1=  data47->GetFunction("g1");
Double_t q0=fit1->GetParameter(0);
Double_t q1=fit1->GetParameter(1);
Double_t q2=fit1->GetParameter(2);
Double_t q3=fit1->GetParameter(3);
Double_t q4=fit1->GetParameter(4);
Double_t q5=fit1->GetParameter(5);


Double_t j1=g1->Eval(0.0523598776);
RISULTATI << "The differential cross section for 47 MeV a 3 degrees is:  " << j1 << endl;

TF1 *g1fit = new TF1("g1fit","([0]*TMath::Gaus(x,[1],[2])+[3]*TMath::Gaus(x,[4],[5]))*TMath::Sin(x)", -0.3, 0.3 );
g1fit->SetParameters(q0,q1,q2,q3,q4,q5);
g1fit->SetLineColor(kBlue);
//g1fit->Draw("same");

Double_t i1=g1fit->Integral(0.,0.3);
RISULTATI << "The integral is: " << i1 << endl; //è l'integrale tra 0 e max_angolo di dsigma/domega*sen(x) --> va moltiplicato per 2*pi

Double_t k1=(2*TMath::Pi()*i1)/j1;
RISULTATI << "The value of K for 47 MeV is: " << k1 << endl << endl;

//_______________________________________________________________________________________________________________

TF1 *g2 = new TF1("g2","[0]*TMath::Gaus(x,[1],[2])+[3]*TMath::Gaus(x,[4],[5])", -0.3, 0.3 );
p0=350;
p1=0;
p2=0.0492699082;
p3=60;
p4=0;
p5=0.0772664626;
g2->SetParameters(p0,p1,p2,p3,p4,p5);

/*g2->SetParLimit(0, 150, 800);
g2->FixParameter(1, 0);
g2->SetParLimit(2, 0.03, 0.06);
g2->SetParLimit(3, 0, 300);
g2->FixParameter(4,0);
g2->SetParLimit(5, 0.07, 0.095);
*/

TCanvas *c2 = new TCanvas("c2");

TGraph *data40=new TGraph("40rad.dat");
data40->GetYaxis()->SetRangeUser(0,400);
data40->Fit("g2","", "", -0.201, 0.201);
data40->Draw("");

data40->SetTitle ("40 MeV");
c2->SaveAs("40MeV_grad.pdf");

TF1 *fit2=  data40->GetFunction("g2");
q0=fit2->GetParameter(0);
q1=fit2->GetParameter(1);
q2=fit2->GetParameter(2);
q3=fit2->GetParameter(3);
q4=fit2->GetParameter(4);
q5=fit2->GetParameter(5);


Double_t j2=g2->Eval(0.0523598776);
RISULTATI << "The differential cross section for 40 MeV a 3 degrees is:  " << j2 << endl;

TF1 *g2fit = new TF1("g2fit","([0]*TMath::Gaus(x,[1],[2])+[3]*TMath::Gaus(x,[4],[5]))*TMath::Sin(x)", -0.3, 0.3 );
g2fit->SetParameters(q0,q1,q2,q3,q4,q5);
g2fit->SetLineColor(kBlue);
//g2fit->Draw("same");

Double_t i2=g2fit->Integral(0.,0.3);
RISULTATI << "The integral is: " << i2 << endl; //è l'integrale tra 0 e max_angolo di dsigma/domega*sen(x) --> va moltiplicato per 2*pi

Double_t k2=(2*TMath::Pi()*i2)/j2;
RISULTATI << "The value of K for 40 MeV is: " << k2 << endl << endl;

//________________________________________________________________________________________________________________

TF1 *g3 = new TF1("g3","[0]*TMath::Gaus(x,[1],[2])+[3]*TMath::Gaus(x,[4],[5])", -0.3, 0.3 );
p0=50;
p1=0;
p2=0.0392699082;
p3=1;
p4=0;
p5=0.0872664626;
g3->SetParameters(p0,p1,p2,p3,p4,p5);
/*
g2->SetParLimit(0, 4000, 9000);
g2->FixParameter(1, 0);
g2->SetParLimit(2, 0.03, 0.06);
g2->SetParLimit(3, 200, 1000);
g2->FixParameter(4,0);
g2->SetParLimit(5, 0.07, 0.095);
*/

TCanvas *c3 = new TCanvas("c3");

TGraph *data37=new TGraph("37rad.dat");
data37->GetYaxis()->SetRangeUser(0,70);
data37->Fit("g3","", "", -0.201, 0.201);
data37->Draw("");

data37->SetTitle ("37 MeV");
c3->SaveAs("37MeV_grad.pdf");


TF1 *fit3=  data37->GetFunction("g3");
q0=fit3->GetParameter(0);
q1=fit3->GetParameter(1);
q2=fit3->GetParameter(2);
q3=fit3->GetParameter(3);
q4=fit3->GetParameter(4);
q5=fit3->GetParameter(5);


Double_t j3=g3->Eval(0.0523598776);
RISULTATI << "The differential cross section for 37 MeV a 3 degrees is:  " << j3 << endl;

TF1 *g3fit = new TF1("g3fit","([0]*TMath::Gaus(x,[1],[2])+[3]*TMath::Gaus(x,[4],[5]))*TMath::Sin(x)", -0.3, 0.3 );
g3fit->SetParameters(q0,q1,q2,q3,q4,q5);
g3fit->SetLineColor(kBlue);
//g3fit->Draw("same");

Double_t i3=g3fit->Integral(0.,0.3);
RISULTATI << "The integral is: " << i3 << endl; //è l'integrale tra 0 e max_angolo di dsigma/domega*sen(x) --> va moltiplicato per 2*pi

Double_t k3=(2*TMath::Pi()*i3)/j3;
RISULTATI << "The value of K for 37 MeV is: " << k3 << endl << endl;

}
