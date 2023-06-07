{

long long int z = 10000000;//z= 1cm= 10^7 nm

double lunghezzadonda[401];
double assorbanza_exp[401];

double assorbanza_teorica[401];

double R;
double rho;
double epsilonm;
double omega_p=1.37*pow(10.,16); //in mn-1
double gamma_bulk=1.08*pow(10.,14); //in mn-1
double v_F=1.4*pow(10.,15); //in nm/s
double omega;
double gamma_corr=0;
double epsilon1_corr;
double epsilon2_corr;

double epsilon1F= -4.43375; //inseriti manualmente dalla condizione di Frolich
double epsilon2F= 2.38035;
//double epsilonmF = -epsilon1F/2; //epsilonmF vale 2.216875
double epsilonmF = 2.28;
double epsilon1[401];
double epsilon2[401];
double lambda_exp[401];

long long int c =299792458*pow(10.0,9.0);//velocit� della luce in nanometri. Ricorda che la lunghezza d'onda � in nanometri

double sigmaext[401];

/*---------------------------------------------------------------------------------------------------------------
Acquisisco i dati sperimentali*/

ifstream AbsorbanceSpectra;
AbsorbanceSpectra.open("AbsorbanceSpectra.dat");
    for(int i=0; i<401; i++){
        AbsorbanceSpectra >> lunghezzadonda[i] >> assorbanza_exp[i];
    }
    AbsorbanceSpectra.close();

    
ifstream BulkDielectricFunction;
    BulkDielectricFunction.open("BulkDielectricFunction.dat");
    for(int k=0;k<401; k++){
        BulkDielectricFunction >> lambda_exp[k] >> epsilon1[k] >> epsilon2[k];
    }
    BulkDielectricFunction.close();

ofstream DatiSperimentali;
    DatiSperimentali.open("DatiSperimentali.dat");
    for(int j=0; j<401; j++){
        DatiSperimentali << lunghezzadonda[j] << "\t" << assorbanza_exp[j] << "\t" << epsilon1[j] << "\t" << epsilon2[j] << endl;
    }
    DatiSperimentali.close();

/*-----------------------------------------------------------------------------
Ricavo Assorbanza_fit_corr*/


R= 3.65;
rho= 10.98*pow(10,-9);
epsilonm= 2.02;
gamma_corr=gamma_bulk + (TMath::Pi()*v_F)/(4.*R);

ofstream Assorbanza_fit_corr;
Assorbanza_fit_corr.open("Assorbanza_fit_corr.dat");

for(int i=0; i<401; i++){
    omega=2*TMath::Pi()*c/lunghezzadonda[i];
    epsilon1_corr=epsilon1[i]+(pow(omega_p,2.))*(1/(pow(omega,2.)+pow(gamma_bulk,2.))-1/(pow(omega,2.)+pow(gamma_corr,2.)));
    epsilon2_corr=epsilon2[i]-(pow(omega_p,2.))/omega*(gamma_bulk/(pow(omega,2.)+pow(gamma_bulk,2.))-gamma_corr/(pow(omega,2.)+pow(gamma_corr,2.)));
    sigmaext[i]= 9.*(2.*TMath::Pi()/lunghezzadonda[i])*pow(epsilonm,1.5)*(4./3.)*TMath::Pi()*pow(R,3)*epsilon2_corr/(pow((epsilon1_corr+2*epsilonm),2)+pow(epsilon2_corr,2));
    assorbanza_teorica[i]= log10(TMath::E())*z*rho*sigmaext[i];

    Assorbanza_fit_corr << lunghezzadonda[i] << "\t" << assorbanza_teorica[i] << "\t" << endl;
    }

Assorbanza_fit_corr.close();


/*-------------------------------------------------------------------------------------------------------
Faccio il grafico dell'assorbanza sperimentale con quella ricavata*/

TMultiGraph *assorbanze = new TMultiGraph ();

TGraph *assorbanza_exp_graph = new TGraph("AbsorbanceSpectra.dat");
	assorbanza_exp_graph->SetLineColor(4);
    assorbanza_exp_graph->SetTitle("Assorbanza sperimentale");
	gStyle->SetTitleFontSize(0.07);
	assorbanza_exp_graph->GetXaxis()->SetRangeUser(400,800);
	assorbanza_exp_graph->GetXaxis()->SetTitle("#lambda (nm)");
	assorbanza_exp_graph->GetXaxis()->SetTitleSize(0.04);
	assorbanza_exp_graph->GetXaxis()->CenterTitle();
	assorbanza_exp_graph->GetYaxis()->SetTitle("Absorbance");
	assorbanza_exp_graph->GetYaxis()->SetTitleSize(0.04);
	assorbanza_exp_graph->GetYaxis()->CenterTitle();
	

TGraph *assorbanza_fit_graph = new TGraph("Assorbanza_fit_corr.dat");
	assorbanza_fit_graph->SetLineColor(3);
    assorbanza_fit_graph->SetTitle("Assorbanza fit");
	gStyle->SetTitleFontSize(0.07);
	assorbanza_fit_graph->GetXaxis()->SetRangeUser(400,800);
	assorbanza_fit_graph->GetXaxis()->SetTitle("#lambda (nm)");
	assorbanza_fit_graph->GetXaxis()->SetTitleSize(0.04);
	assorbanza_fit_graph->GetXaxis()->CenterTitle();
	assorbanza_fit_graph->GetYaxis()->SetTitle("Absorbance");
	assorbanza_fit_graph->GetYaxis()->SetTitleSize(0.04);
	assorbanza_fit_graph->GetYaxis()->CenterTitle();

    assorbanze->Add(assorbanza_exp_graph);
    assorbanze->Add(assorbanza_fit_graph);
	assorbanze->Draw("AL");

    TLegend *legenda = new TLegend(0.5,0.5,0.5,0.5);
    legenda->AddEntry(assorbanza_exp_graph, "assorbanza sperimentale");
    legenda->AddEntry(assorbanza_fit_graph, "assorbanza chi minimo");
    
    legenda->Draw();

   
	


}

