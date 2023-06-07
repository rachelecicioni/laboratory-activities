{

#define r 300 //dimensione di R
#define p 1000 //dimensione di rho
#define m 150 //dimensione di epsilonm
#define L 401//numero di lunghezze donda

int i_min=59; //99 è 500nm, 79 è 480 nm, 59 è 460 nm
int i_max=189; //179 è 580 nm
long long int z = 10000000;//z= 1cm= 10^7 nm

double lunghezzadonda[L];
double assorbanza_exp[L];

double assorbanza_teorica[L];

double R[r];
double rho[p];

double epsilon1F= -4.43375; //inseriti manualmente dalla condizione di Frolich
double epsilon2F= 2.38035;
double epsilonmF = -epsilon1F/2; //epsilonmF vale 2.216875
//double epsilonmF = 2.02;
double epsilon1[L];
double epsilon2[L];
double lambda_exp[L];

long long int c =299792458*pow(10.0,9.0);//velocit� della luce in nanometri. Ricorda che la lunghezza d'onda � in nanometri

double sigmaext[L];

double chiquadro_R_rho[r][p];
double chiquadro_R_rho_log[r][p];

for(int k=0;k<r; k++){
       for(int l=0;l<m; l++){
           chiquadro_R_rho[k][l]=0;
           chiquadro_R_rho_log[k][l]=0;
       }
    }

/*---------------------------------------------------------------------------------------------------------------
Acquisisco i dati sperimentali*/   //!!!lunghezzadonda e lambda_exp sono identici

for(int k=0;k<r; k++){
       for(int n=0;n<p; n++){
           chiquadro_R_rho[k][n]=0;
           chiquadro_R_rho_log[k][n]=0;
       }
    }

ifstream AbsorbanceSpectra;
AbsorbanceSpectra.open("AbsorbanceSpectra.dat");
    for(int i=0; i<L; i++){
        AbsorbanceSpectra >> lunghezzadonda[i] >> assorbanza_exp[i];
    }
    AbsorbanceSpectra.close();

    
ifstream BulkDielectricFunction;
    BulkDielectricFunction.open("BulkDielectricFunction.dat");
    for(int k=0;k<L; k++){
        BulkDielectricFunction >> lambda_exp[k] >> epsilon1[k] >> epsilon2[k];
    }
    BulkDielectricFunction.close();

ofstream DatiSperimentali;
    DatiSperimentali.open("DatiSperimentali.dat");
    for(int j=0; j<L; j++){
        DatiSperimentali << lunghezzadonda[j] << "\t" << assorbanza_exp[j] << "\t" << epsilon1[j] << "\t" << epsilon2[j] << endl;
    }
    DatiSperimentali.close();

/*------------------------------------------------------------------------------------------------------------------
Faccio variare (R,rho) con le correzioni di taglia*/

//Definisco le costanti necessarie per la correzione di taglia
/*
double omega_p=1.37*pow(10.,16); //in mn-1
double gamma_bulk=1.08*pow(10.,14); //in mn-1
double v_F=1.4*pow(10.,15); //in nm/s
double omega;
double gamma_corr=0;
double epsilon1_corr;
double epsilon2_corr;

ofstream file_chiquadro_R_rho_corr;//genera il file con ChiQuadrorho e rho
file_chiquadro_R_rho_corr.open("file_chiquadro_R_rho_corr.dat");
ofstream file_chiquadro_minimo_R_rho_corr;
file_chiquadro_minimo_R_rho_corr.open("file_chiquadro_minimo_R_rho_corr.dat");
ofstream assorbanza_R_rho_corr;
assorbanza_R_rho_corr.open("assorbanza_R_rho_corr.dat");



R[0]=3.0; //R in nm e va da 0.01 a 10 nm con passo 0.01
rho[0]=100*pow(10,-11); //rho in nm^-3 e va da 1e-11 nm a 1e-8 nm^-3 con step 1*10-11

double chiquadro_R_rho_minimo=1000000000000;
double R_minimo_rho=0;
double rho_minimo=0;


for(int k=0;k<r; k++){
    if (k<r-1){
    R[k+1]=R[k]+0.01;
        }
            
    for(int n=0; n<p; n++){
        if (n<p-1){
        rho[n+1]=rho[n]+pow(10,-11);
            }
        
        for(int i=0; i<401; i++){
            omega=2*TMath::Pi()*c/lunghezzadonda[i];
            gamma_corr=gamma_bulk + (TMath::Pi()*v_F)/(4.*R[k]);
            epsilon1_corr=epsilon1[i]+(pow(omega_p,2.))*(1/(pow(omega,2.)+pow(gamma_bulk,2.))-1/(pow(omega,2.)+pow(gamma_corr,2.)));
            epsilon2_corr=epsilon2[i]-((pow(omega_p,2.))/omega)*(gamma_bulk/(pow(omega,2.)+pow(gamma_bulk,2.))-gamma_corr/(pow(omega,2.)+pow(gamma_corr,2.)));
            sigmaext[i]= 9.*(2.*TMath::Pi()/lunghezzadonda[i])*pow(epsilonmF,1.5)*(4./3.)*TMath::Pi()*pow(R[k],3.)*epsilon2_corr/(pow((epsilon1_corr+2*epsilonmF),2.)+pow(epsilon2_corr,2.));
            assorbanza_teorica[i]= log10(TMath::E())*z*rho[n]*sigmaext[i];
        if(i>=i_min && i <=i_max){
            chiquadro_R_rho[k][n]+=pow((assorbanza_exp[i]-assorbanza_teorica[i]),2)/assorbanza_teorica[i];  
        }
        }
        if (chiquadro_R_rho[k][n]<=chiquadro_R_rho_minimo){
                chiquadro_R_rho_minimo=chiquadro_R_rho[k][n];
                R_minimo_rho=R[k];
                rho_minimo=rho[n];
            }
    file_chiquadro_R_rho_corr <<  R[k] << "\t"<< rho[n] << "\t" << chiquadro_R_rho[k][n] << endl;
    chiquadro_R_rho_log[k][n] = log(chiquadro_R_rho[k][n]);
    }
    cout << "Indice " << k << endl;
}

file_chiquadro_minimo_R_rho_corr << "rho minimo: " << rho_minimo << "\t" << "R minimo: " << R_minimo_rho << "\t" << "ChiQuadro minimo: " << "\t" << chiquadro_R_rho_minimo << endl;

//File con assorbanza_R_rho_corr

for(int i=0; i<401; i++){
    omega=2*TMath::Pi()*c/lunghezzadonda[i];
    epsilon1_corr=epsilon1[i]+(pow(omega_p,2.))*(1/(pow(omega,2.)+pow(gamma_bulk,2.))-1/(pow(omega,2.)+pow(gamma_corr,2.)));
    epsilon2_corr=epsilon2[i]-(pow(omega_p,2.))/omega*(gamma_bulk/(pow(omega,2.)+pow(gamma_bulk,2.))-gamma_corr/(pow(omega,2.)+pow(gamma_corr,2.)));
    sigmaext[i]= 9.*(2.*TMath::Pi()/lunghezzadonda[i])*pow(epsilonmF,1.5)*(4./3.)*TMath::Pi()*pow(R_minimo_rho,3)*epsilon2_corr/(pow((epsilon1_corr+2*epsilonmF),2)+pow(epsilon2_corr,2));
    assorbanza_teorica[i]= log10(TMath::E())*z*rho_minimo*sigmaext[i];
    assorbanza_R_rho_corr << lunghezzadonda[i] << "\t" << assorbanza_teorica[i] << "\t" << endl;
    }

assorbanza_R_rho_corr.close();
*/
/*----------------------------------------------------------------------------------------------------------
grafico l'andamento del chiquadro*/
/*
int R_histo[r];
int rho_histo[p];

R_histo[0]=1;
rho_histo[0]=1;

for (int k=0; k<(r-1); k++){
			R_histo[k+1]=R_histo[k] + 1; 
	for (int n=0; n<(p-1); n++){
			rho_histo[n+1] = rho_histo[n] + 1; 
	}
}

TH2D *histo2D_rho = new TH2D("histo2D_rho","Chi Square",r,R_histo[0],R_histo[r-1],p,rho_histo[0],rho_histo[p-1]);

histo2D_rho->GetXaxis()->SetTitle("R (10^{-1} nm)");
histo2D_rho->GetXaxis()->SetTitleSize(0.045);
histo2D_rho->GetXaxis()->SetTickSize(0.02);
histo2D_rho->GetXaxis()->CenterTitle();

histo2D_rho->GetYaxis()->SetTitle("#rho (10^{-11})");
histo2D_rho->GetYaxis()->SetTitleSize(0.045);
histo2D_rho->GetYaxis()->SetTickSize(0.02);
histo2D_rho->GetYaxis()->CenterTitle();

histo2D_rho->SetMaximum(10);

for(int k=0; k<r; k++) {
    for (int n=0; n<p; n++){
    histo2D_rho->SetBinContent(R_histo[k], rho_histo[n], chiquadro_R_rho_log[k][n]);
    }
}


TCanvas *c2 = new TCanvas("c2");
c2->SetCanvasSize(1047 , 711);
c2->SetWindowSize(1067, 741);
	
gStyle->SetOptStat(0);
histo2D_rho->SetContour(10000);
gStyle->SetPalette(56); //il nome della palette è kInvertedDarkBodyRadiator=56
	
histo2D_rho->Draw("colz");
c2->SaveAs("histo_2d_R_rho_corr.png");

file_chiquadro_R_rho_corr.close();
file_chiquadro_minimo_R_rho_corr.close();
*/
/*--------------------------------------------------------------------------------------------------------------------------
Faccio variare (R,epsilonm) con le correzioni di taglia*/
/*
double epsilonm[m]; //R già ce l'ho e anche rho è il valore minimo
double chiquadro_R_epsilonm[r][m];
double chiquadro_R_epsilonm_log[r][m];

for(int k=0;k<r; k++){
       for(int l=0;l<m; l++){
           chiquadro_R_epsilonm[k][l]=0;
           chiquadro_R_epsilonm_log[k][l]=0;
       }
    }

ofstream file_chiquadro_R_epsilonm_corr;//genera il file con ChiQuadrorho e epsilonm
file_chiquadro_R_epsilonm_corr.open("file_chiquadro_R_epsilonm_corr.dat");
ofstream file_chiquadro_minimo_R_epsilonm_corr;
file_chiquadro_minimo_R_epsilonm_corr.open("file_chiquadro_minimo_R_epsilonm_corr.dat");
ofstream assorbanza_R_epsilonm_corr;
assorbanza_R_epsilonm_corr.open("assorbanza_R_epsilonm_corr.dat");

R[0]=3.0; //va da 0.01 a 10nm con passo 0.01
epsilonm[0]=1.3; //va da 1 a 3 con passo 0.01

double chiquadro_R_epsilonm_minimo=1000000000000;
double R_minimo_epsilonm=0;
double epsilonm_minimo=0;

for(int k=0;k<r; k++){
    if (k<r-1){
        R[k+1]=R[k]+0.01;
            }
    
    for(int l=0; l<m; l++){
        if (l<m-1){
        epsilonm[l+1]=epsilonm[l]+0.01;
            }

        for(int i=0; i<401; i++){
            omega=2*TMath::Pi()*c/lunghezzadonda[i];
            gamma_corr=gamma_bulk + (TMath::Pi()*v_F)/(4.*R[k]);
            epsilon1_corr=epsilon1[i]+(pow(omega_p,2.))*(1/(pow(omega,2.)+pow(gamma_bulk,2.))-1/(pow(omega,2.)+pow(gamma_corr,2.)));
            epsilon2_corr=epsilon2[i]-(pow(omega_p,2.))/omega*(gamma_bulk/(pow(omega,2.)+pow(gamma_bulk,2.))-gamma_corr/(pow(omega,2.)+pow(gamma_corr,2.)));
            sigmaext[i]= 9.*(2.*TMath::Pi()/lunghezzadonda[i])*pow(epsilonm[l],1.5)*(4./3.)*TMath::Pi()*pow(R[k],3)*epsilon2_corr/(pow((epsilon1_corr+2*epsilonm[l]),2)+pow(epsilon2_corr,2));
            assorbanza_teorica[i]= log10(TMath::E())*z*rho_minimo*sigmaext[i];
        if(i>=i_min && i <=i_max){
            chiquadro_R_epsilonm[k][l]+=pow((assorbanza_exp[i]-assorbanza_teorica[i]),2)/assorbanza_teorica[i];
            
        }
        }
        if (chiquadro_R_epsilonm[k][l]<=chiquadro_R_epsilonm_minimo){
                chiquadro_R_epsilonm_minimo=chiquadro_R_epsilonm[k][l];
                R_minimo_epsilonm=R[k];
                epsilonm_minimo=epsilonm[l];
            }
    file_chiquadro_R_epsilonm_corr <<  R[k] << "\t"<< epsilonm[l] << "\t" << chiquadro_R_epsilonm[k][l] << endl;
    chiquadro_R_epsilonm_log[k][l] = log(chiquadro_R_epsilonm[k][l]);
    }
    cout << "Indice " << k << endl;
}


file_chiquadro_minimo_R_epsilonm_corr << "epsilonm minimo: " << epsilonm_minimo << "\t" << "R minimo: " << R_minimo_epsilonm << "\t" << "ChiQuadro minimo: " << "\t" << chiquadro_R_epsilonm_minimo << endl;

//File con assorbanza_R_epsilonm_corr

for(int i=0; i<401; i++){
    omega=2*TMath::Pi()*c/lunghezzadonda[i];
    epsilon1_corr=epsilon1[i]+(pow(omega_p,2.))*(1/(pow(omega,2.)+pow(gamma_bulk,2.))-1/(pow(omega,2.)+pow(gamma_corr,2.)));
    epsilon2_corr=epsilon2[i]-(pow(omega_p,2.))/omega*(gamma_bulk/(pow(omega,2.)+pow(gamma_bulk,2.))-gamma_corr/(pow(omega,2.)+pow(gamma_corr,2.)));
    sigmaext[i]= 9.*(2.*TMath::Pi()/lunghezzadonda[i])*pow(epsilonm_minimo,1.5)*(4./3.)*TMath::Pi()*pow(R_minimo_epsilonm,3)*epsilon2_corr/(pow((epsilon1_corr+2*epsilonm_minimo),2)+pow(epsilon2_corr,2));
    assorbanza_teorica[i]= log10(TMath::E())*z*rho_minimo*sigmaext[i];
    assorbanza_R_epsilonm_corr << lunghezzadonda[i] << "\t" << assorbanza_teorica[i] << "\t" << endl;
    }

assorbanza_R_epsilonm_corr.close();
*/
/*---------------------------------------------------------------------------------------------------------------------
Grafichiamo l'istogramma 2D variando (R,epsilonm) con rho_minimo ottenuto dal predente grafico 2D*/
/*
int epsilonm_histo[m];

R_histo[0]=1;
epsilonm_histo[0]=1;

for (int k=0; k<(r-1); k++){
			R_histo[k+1]=R_histo[k] + 1; }
	for (int l=0; l<(m-1); l++){
			epsilonm_histo[l+1] = epsilonm_histo[l] + 1; 
	}


TH2D *histo2D_epsilonm = new TH2D("histo2D_epsilonm","Chi Square",r,R_histo[0],R_histo[r-1],m,epsilonm_histo[0],epsilonm_histo[m-1]);

histo2D_epsilonm->GetXaxis()->SetTitle("R (10^{-1} nm)");
histo2D_epsilonm->GetXaxis()->SetTitleSize(0.045);
histo2D_epsilonm->GetXaxis()->SetTickSize(0.02);
histo2D_epsilonm->GetXaxis()->CenterTitle();

histo2D_epsilonm->GetYaxis()->SetTitle("#epsilon_{m} (10^{-2})");
histo2D_epsilonm->GetYaxis()->SetTitleSize(0.045);
histo2D_epsilonm->GetYaxis()->SetTickSize(0.02);
histo2D_epsilonm->GetYaxis()->CenterTitle();

histo2D_epsilonm->SetMaximum(4);

for(int k=0; k<r; k++) {
    for (int l=0; l<m; l++){
    histo2D_epsilonm->SetBinContent(R_histo[k], epsilonm_histo[l], chiquadro_R_epsilonm_log[k][l]);
    }
}


TCanvas *c3 = new TCanvas("c3");
c3->SetCanvasSize(1047 , 711);
c3->SetWindowSize(1067, 741);
	
gStyle->SetOptStat(0);
histo2D_epsilonm->SetContour(10000);
gStyle->SetPalette(56); //il nome della palette è kInvertedDarkBodyRadiator=56
	
histo2D_epsilonm->Draw("colz");
c3->SaveAs("histo_2d_R_epsilonm_corr.png");

file_chiquadro_R_epsilonm_corr.close();
file_chiquadro_minimo_R_epsilonm_corr.close();
*/

/*-------------------------------------------------------------------------------------------------------
Faccio il grafico dell'assorbanza sperimentale, assorbanza_R_rho, assorbanza_R_epsilonm*/
/*
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
    
TGraph *assorbanza_R_rho_corr_graph = new TGraph("assorbanza_R_rho_corr.dat");
	assorbanza_R_rho_corr_graph->SetLineColor(3);
	assorbanza_R_rho_corr_graph->SetTitle("Assorbanza sperimentale");
	gStyle->SetTitleFontSize(0.07);
	assorbanza_R_rho_corr_graph->GetXaxis()->SetRangeUser(400,800);
	assorbanza_R_rho_corr_graph->GetXaxis()->SetTitle("#lambda (nm)");
	assorbanza_R_rho_corr_graph->GetXaxis()->SetTitleSize(0.04); 
	assorbanza_R_rho_corr_graph->GetXaxis()->CenterTitle();
	assorbanza_R_rho_corr_graph->GetYaxis()->SetTitle("Absorbance");
	assorbanza_R_rho_corr_graph->GetYaxis()->SetTitleSize(0.04);
	assorbanza_R_rho_corr_graph->GetYaxis()->CenterTitle();

TGraph *assorbanza_R_epsilonm_corr_graph = new TGraph("assorbanza_R_epsilonm_corr.dat");
	assorbanza_R_epsilonm_corr_graph->SetLineColor(2);
	assorbanza_R_epsilonm_corr_graph->SetTitle("Assorbanza sperimentale");
	gStyle->SetTitleFontSize(0.07);
	assorbanza_R_epsilonm_corr_graph->GetXaxis()->SetRangeUser(400,800);
	assorbanza_R_epsilonm_corr_graph->GetXaxis()->SetTitle("#lambda (nm)");
	assorbanza_R_epsilonm_corr_graph->GetXaxis()->SetTitleSize(0.04); 
	assorbanza_R_epsilonm_corr_graph->GetXaxis()->CenterTitle();
	assorbanza_R_epsilonm_corr_graph->GetYaxis()->SetTitle("Absorbance");
	assorbanza_R_epsilonm_corr_graph->GetYaxis()->SetTitleSize(0.04);
	assorbanza_R_epsilonm_corr_graph->GetYaxis()->CenterTitle();

	TCanvas *c1 = new TCanvas("c1");
    c1->SetCanvasSize(1047 , 711);
	c1->SetWindowSize(1067, 741);
	
assorbanze->Add(assorbanza_exp_graph);
assorbanze->Add(assorbanza_R_rho_corr_graph);
assorbanze->Add(assorbanza_R_epsilonm_corr_graph);

assorbanze->Draw("AL");

TLegend *legenda = new TLegend(0.5,0.5,0.5,0.5);
legenda->AddEntry(assorbanza_exp_graph, "assorbanza sperimentale");
legenda->AddEntry(assorbanza_R_rho_corr_graph, "assorbanza (R,rho)");
legenda->AddEntry(assorbanza_R_epsilonm_corr_graph, "assorbanza (R,epsilonm)");
legenda->Draw();

	c1->SaveAs("graph_assorbanze_corr.png");
*/
/*------------------------------------------------------------------------------------------------------------------------------------------------
Vario tutti e tre i parametri*/

#define r3 300
#define p3 1000
#define m3 70

double omega_p=1.37*pow(10.,16); //in mn-1
double gamma_bulk=1.08*pow(10.,14); //in mn-1
double v_F=1.4*pow(10.,15); //in nm/s
double omega;
double gamma_corr=0;
double epsilon1_corr;
double epsilon2_corr;


double R_3[r3];
double rho_3[p3];
double epsilonm_3[m3];
double chiquadro_3[r3][p3][m3];
double chiquadro_3_log[r3][p3][m3];

for(int k=0;k<r3; k++){
       for(int n=0;n<p3; n++){
           for(int l=0; l<m3; l++){
           chiquadro_3[k][n][l]=0;
           chiquadro_3_log[k][n][l]=0;
           }
       }
    }

ofstream file_chiquadro_3_corr;
file_chiquadro_3_corr.open("file_chiquadro_3_corr.dat");
ofstream file_chiquadro_minimo_3_corr;
file_chiquadro_minimo_3_corr.open("file_chiquadro_minimo_3_corr.dat");

R_3[0]=3.00; 
rho_3[0]=1*pow(10,-11); 
epsilonm_3[0]=1.90; 
chiquadro_3[0][0][0]=0;

double chiquadro_3_minimo=10000;
double R_3_minimo=0;
double rho_3_minimo=0;
double epsilonm_3_minimo=0;

for (int k=0; k<r3; k++){
    if(k<(r3-1)){
    R_3[k+1]=R_3[k]+0.01;}

    for(int n=0; n<p3; n++){
        if(n<(p3-1)){
        rho_3[n+1]=rho_3[n]+pow(10,-11);}

        for(int l=0; l<m3; l++){
            if(l<(m3-1)){
            epsilonm_3[l+1]=epsilonm_3[l]+0.01;}
            
        for(int i=0; i<401; i++){
        omega=2*TMath::Pi()*c/lunghezzadonda[i];
        gamma_corr=gamma_bulk + (TMath::Pi()*v_F)/(4.*R_3[k]);
        epsilon1_corr=epsilon1[i]+(pow(omega_p,2.))*(1/(pow(omega,2.)+pow(gamma_bulk,2.))-1/(pow(omega,2.)+pow(gamma_corr,2.)));
        epsilon2_corr=epsilon2[i]-((pow(omega_p,2.))/omega)*(gamma_bulk/(pow(omega,2.)+pow(gamma_bulk,2.))-gamma_corr/(pow(omega,2.)+pow(gamma_corr,2.)));
        sigmaext[i]= 9.*(2.*TMath::Pi()/lunghezzadonda[i])*pow(epsilonm_3[l],1.5)*(4./3.)*TMath::Pi()*pow(R_3[k],3)*epsilon2_corr/(pow((epsilon1_corr+2*epsilonm_3[l]),2)+pow(epsilon2_corr,2));
        assorbanza_teorica[i]= log10(TMath::E())*z*rho_3[n]*sigmaext[i];
        if(i>=i_min && i <=i_max){
            chiquadro_3[k][n][l]+=pow((assorbanza_exp[i]-assorbanza_teorica[i]),2)/(1/assorbanza_teorica[i]);
            }
        }
        if (chiquadro_3[k][n][l]<=chiquadro_3_minimo){
                chiquadro_3_minimo=chiquadro_3[k][n][l];
                R_3_minimo=R_3[k];
                rho_3_minimo=rho_3[n];
                epsilonm_3_minimo=epsilonm_3[l];
            }
    file_chiquadro_3_corr << R_3[k] << "\t"<< rho_3[n] << "\t" << epsilonm_3[l] << "\t" <<chiquadro_3[k][n][l] << endl;
    chiquadro_3_log[k][n][l] = log(chiquadro_3[k][n][l]);

        }

    }
    cout << "Sono all'indice " << k << endl;
}

file_chiquadro_minimo_3_corr << "rho minimo: " << rho_3_minimo << "\t" << "R minimo: " << R_3_minimo << "\t" << "epsilonm minimo: " << epsilonm_3_minimo <<  "\t" << "ChiQuadro minimo: " << "\t" << chiquadro_3_minimo << endl;

file_chiquadro_3_corr.close();
file_chiquadro_minimo_3_corr.close();


}//FINE
