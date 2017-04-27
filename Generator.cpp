
 #include <stdio.h>
 #include <stdlib.h>
 #include <string>
 #include <cstring>
 #include <iostream>
 #include <fstream>
 #include <math.h>
 using namespace std;

 #include "Riostream.h"
 #include "TApplication.h"
 #include "TROOT.h"
 #include "TFile.h"
 #include "TNtuple.h"
 #include "TMath.h"
 #include "TVector3.h"
 #include "TLorentzVector.h"
 #include "TF1.h"
 #include "TRandom.h"
 #include  "TH2.h"

int main(){
//void Generator(){
    
    Double_t a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12;
    Double_t a1p, a2p, a3p, a4p, a5p, a6p, a7p, a8p, a9p, a10p, a11p, a12p;
    
    //Sigmas
    
    const Double_t dist1=5.0; //meter
    const Double_t dist2=4.5; //meter
    const Double_t disttoZDC=14.0; //meter
    
    /*
    const Double_t sigma_phi_P1=0.0; //degrees dphi= 0.005/4.0
    const Double_t sigma_theta_P1=0.0; //degrees for dtheta=0.012/2/4.0
    const Double_t sigma_phi_P2=0.0; //degrees dphi= 0.005/5.0
    const Double_t sigma_theta_P2=0.0; //degrees for dtheta=0.012/2/5.0
    const Double_t sigma_phi_Precoil=0.; //tan(0.10/2/14 m )
    const Double_t sigma_theta_Precoil=0.; //tan(0.10/2/14 m )
    
    const Double_t sigma_timeRPC=0.00; //nsec
    const Double_t sigma_timeZDC=0.00; //ZDC
    
    const Double_t sigma_Pbeam=0.0; //%
    const Double_t sigma_theta_Pbeam=0.; //degrees
    
    const Double_t sigma_phi_Pcm=0.; //tan(0.10/2/14 m )    ||| 0.1 from the tracking before the magnet
    const Double_t sigma_theta_Pcm=0.; //tan(0.10/2/14 m )  ||| 0.1 from the tracking before the magnet
    const Double_t sigma_Pcm= 0.0; // %
    */
    /**/
    const Double_t sigma_phi_P1=0.06; //degrees dphi= 0.005/sqrt(12)/5.0
    const Double_t sigma_theta_P1=0.07; //degrees for dtheta=0.012/sqrt(12)/5.0
    const Double_t sigma_phi_P2=0.07; //degrees dphi= 0.005/sqrt(12)/4.5
    const Double_t sigma_theta_P2=0.08; //degrees for dtheta=0.012/sqrt(12)/4.5
    const Double_t sigma_phi_Precoil=0.1; //tan(0.10/sqrt(12)/14 m )
    const Double_t sigma_theta_Precoil=0.1; //tan(0.10/sqrt(12)/14 m )
    
    const Double_t sigma_timeRPC=0.1; //nsec
    const Double_t sigma_timeZDC=0.15; //ZDC
    
    const Double_t sigma_Pbeam=0;//0.01; //%
    const Double_t sigma_theta_Pbeam=0;//0.05; //degrees
    
    const Double_t sigma_phi_Pcm=0.1; //tan(0.10/2/14 m )
    const Double_t sigma_theta_Pcm=0.1; //tan(0.10/2/14 m )
    const Double_t sigma_Pcm= 0.04; // %
    /**/
    //const Double_t A=10.;
    const Double_t AB11=11.;
    const Double_t Ac12=12.;
    
    //Energy loss calculation constants
    Double_t factor = 0.1535; //MeVcm^2/gr 2*TMath::Pi()*Na*me*c^2
    Double_t me = 0.5109989461*pow(10,6); // MeV/c^2 ;
    //Double_t c = 3*10^8 //m/sec;
    
    Double_t density =  1.032; //silicon gr/cm^3
    Double_t atomic_weight = 13; //atomic weight of Silicon
    Double_t atomic_number = 7.0; //atomic number of Silicon
    Double_t I = 64.7; //mean excitation potential for Silicon
    Double_t charge_B = 5; //charge of incident particle boron;
    Double_t charge_Be = 4; //charge of incident particle Beryllium;
    
    //Double_t density =  2.328; //silicon gr/cm^3
    //Double_t atomic_weight = 28.0855; //atomic weight of Silicon
    //Double_t atomic_number = 14; //atomic number of Silicon
    //Double_t I = 173; //mean excitation potential for Silicon
    
    //Double_t density = 8.96; // copper gr/cm^3
    //Double_t atomic_weight = 63.546; //atomic weight of Copper
    //Double_t atomic_number = 29; //atomic number of Copper
    //Double_t I = 322; // ev mean excitation potential for Copper
    
    // Proton beam energy (momentum)
    const Double_t mp = 0.938272;
    const Double_t mn=  0.939565;
    const Double_t mC12 = 11.1750;
    const Double_t mB10 = 9.3245;
    
    
    Double_t Pbeam = 4.0;
    //Double_t Ek_beam = TMath::Sqrt( pow(Pbeam,2) + pow(mp,2) ) - mp;
    
    //Ek_beam = sqrt(pow(Pbeam,2) + pow(mp,2))-mp;
    //Double_t beta=TMath::Sqrt( pow((Ek_beam+mp+mp)/2,2) - pow(mp,2) );

    cout<<Pbeam<<endl;

    
    
    
    
    
    
    //ub/Sr ub->b b->cm2 Sr Beam Target_thickness transparency neutron_efficiency Time=14 days 90->other
    Double_t time = 14; //days
    Double_t duty_cycle = 2./12; //half time on the accelerator
    Double_t hrstosec = 3600; //hrs to sec
    Double_t hrs = 24; //hrs per day
    Double_t flux = 5e5;
    Double_t transparency = 0.4;
    Double_t neutron_efficiency = 1.0;
    Double_t mbtob = 1e-3;
    Double_t btocmsq = 1e-24;
    Double_t target_thickness = 6.3e23 * 6; // (OR: 15cm LH2 * 6 protons for each 12C beam nucleus)    1/cm^2 = 10%/radius where radius is R=r0*A^1/3 where r0=1.2 fm
    //Double_t factor_solid_angle = 0.08; // 0.25*0.7
    Double_t factor_solid_angle = 1.0;
    
    
    
    
    
    //TF1* F_CrossSection = new TF1("F_CrossSection","expo",3,10); //pp cross section in uB/Sr in the c.m.
    //F_CrossSection->SetParameter(0,8.44564);
    //F_CrossSection->SetParameter(1,-1.25103);
    
    TF1* F_CrossSection = new TF1("F_CrossSection","expo(0)+expo(2)+expo(4)",1.0,14.0); //pp cross section in mB/Sr in the c.m.
    F_CrossSection->SetParameter(0, 6.53798e+00);
    F_CrossSection->SetParameter(1,-4.91538e+00);
    F_CrossSection->SetParameter(2, 1.23210e+00);
    F_CrossSection->SetParameter(3,-1.22331e+00);
    F_CrossSection->SetParameter(4,-6.87842e+00);
    F_CrossSection->SetParameter(5,-3.73291e-01);
    
    
    TF1  *f11 = new TF1("f1","expo",1,2);
    f11->SetParameter(0,4.85781e+00);
    f11->SetParameter(1,-3.22678e+00);
    TF1  *f22 = new TF1("f2","expo",2,8);
    f22->SetParameter(0,1.93109e+00);
    f22->SetParameter(1,-1.32115e+00);
    TF1  *f33 = new TF1("f3","expo",8.0,14.0);
    f33->SetParameter(0,-3.50443e+00);
    f33->SetParameter(1,-6.39896e-01);
    
    //
    // Tree to hold (p,2pn) 'SRC' events and extra kinematical variables.
    // Each Event has two forward scatter protons and a recoil neutron.
    // The reconstructed initial momentum of the proton knockout of the nucleus is lables as "missing momentum"
    //
    TFile* f = new TFile("GSI_4GeVpercperu.root","recreate");
    TTree* T = new TTree("T","T");
    
    Double_t theta_cm;
    Double_t weight;
    Double_t Effective_E;
    Double_t cross_section;
    
    Double_t time_P1, time_P2, time_Precoil; //time_Pcm,
    Double_t time_P1_rec, time_P2_rec, time_Precoil_rec; //time_Pcm_rec,
    
    Double_t theta_Pbeam_smear_P, phi_Pbeam_smear_P, Pbeam_smear_P, P_beam_smear_P[3], Ebeam_smear_P;
    Double_t theta_Pbeam_smear_C12, phi_Pbeam_smear_C12, Pbeam_smear_C12, P_beam_smear_C12[3], Ebeam_smear_C12;
    
    Double_t theta_P1, phi_P1,  P1, P_1[3], E1;
    Double_t theta_P2, phi_P2,  P2, P_2[3], E2;
    Double_t theta_Precoil, phi_Precoil, Precoil, P_recoil[3], Erecoil;
    Double_t theta_Pmiss, phi_Pmiss, Pmiss, P_miss[3], Emiss;
    Double_t Pcm, P_cm[3], Ecm;
    Double_t t,u,s;
    
    Double_t theta_Pbeam_boost, phi_Pbeam_boost, Pbeam_boost, P_beam_boost[3], Ebeam_boost;
    Double_t theta_P1_boost, phi_P1_boost,  P1_boost, P_1_boost[3], E1_boost;
    Double_t theta_P2_boost, phi_P2_boost,  P2_boost, P_2_boost[3], E2_boost;
    Double_t theta_Precoil_boost, phi_Precoil_boost, Precoil_boost, P_recoil_boost[3], Erecoil_boost;
    Double_t theta_Pcm_boost, phi_Pcm_boost, Pcm_boost, P_cm_boost[3], Ecm_boost, Ecm_11_boost;
    Double_t t_boost,u_boost,s_boost;
    Double_t P2_angle_boost, P_2_angle_boost[3], E2_angle_boost, theta_P2_angle_boost, phi_P2_angle_boost;
    Double_t Precoil_angle_boost, P_recoil_angle_boost[3], Erecoil_angle_boost, theta_Precoil_angle_boost, phi_Precoil_angle_boost;
    Double_t Pcm_angle_boost, P_cm_angle_boost[3], Ecm_angle_boost, theta_Pcm_angle_boost, phi_Pcm_angle_boost;
    Double_t beta_cm_B, gamma_cm_B;
    Double_t beta_cm_B11, gamma_cm_B11;
    Double_t beta_cm_Be, gamma_cm_Be;
    Double_t dEdX_B; //MeV*cm^2/gr
    Double_t dEdX_density_B; //MeV/cm
    Double_t dEdX_B11; //MeV*cm^2/gr
    Double_t dEdX_density_B11; //MeV/cm
    Double_t dEdX_Be; //MeV*cm^2/gr
    Double_t dEdX_density_Be; //MeV/cm
    
    Double_t theta_P1_rec_boost, phi_P1_rec_boost, P1_rec_boost, P_1_rec_boost[3], E1_rec_boost;
    Double_t theta_P2_rec_boost, phi_P2_rec_boost, P2_rec_boost, P_2_rec_boost[3], E2_rec_boost;
    Double_t theta_Precoil_rec_boost, phi_Precoil_rec_boost, Precoil_rec_boost, P_recoil_rec_boost[3], Erecoil_rec_boost;
    Double_t theta_Pcm_rec_boost, phi_Pcm_rec_boost, Pcm_rec_boost, P_cm_rec_boost[3], Ecm_rec_boost;
    Double_t t_rec_boost,u_rec_boost,s_rec_boost;
    Double_t t_angle_rec_boost,u_angle_rec_boost,s_angle_rec_boost;
    Double_t P2_angle_rec_boost, P_2_angle_rec_boost[3], E2_angle_rec_boost, theta_P2_angle_rec_boost, phi_P2_angle_rec_boost;
    Double_t Precoil_angle_rec_boost, P_recoil_angle_rec_boost[3], Erecoil_angle_rec_boost, theta_Precoil_angle_rec_boost, phi_Precoil_angle_rec_boost;
    Double_t Pcm_angle_rec_boost, P_cm_angle_rec_boost[3], Ecm_angle_rec_boost, theta_Pcm_angle_rec_boost, phi_Pcm_angle_rec_boost;
    
    Double_t theta_Pbeam_rec, phi_Pbeam_rec, Pbeam_rec, P_beam_rec[3], Ebeam_rec;
    Double_t theta_P1_rec, phi_P1_rec, P1_rec, P_1_rec[3], E1_rec;
    Double_t theta_P2_rec, phi_P2_rec, P2_rec, P_2_rec[3], E2_rec;
    Double_t theta_Precoil_rec, phi_Precoil_rec, Precoil_rec, P_recoil_rec[3], Erecoil_rec;
    Double_t theta_Pmiss_rec, phi_Pmiss_rec, Pmiss_rec, P_miss_rec[3], Emiss_rec;
    Double_t theta_Pcm_rec, phi_Pcm_rec, Pcm_rec, P_cm_rec[3], Ecm_rec;
    Double_t t_rec,u_rec,s_rec;
    Double_t t_angle_rec,u_angle_rec,s_angle_rec;
    Double_t P2_angle_rec, P_2_angle_rec[3], E2_angle_rec, theta_P2_angle_rec, phi_P2_angle_rec;
    Double_t Pmiss_angle_rec, P_miss_angle_rec[3], Emiss_angle_rec, theta_Pmiss_angle_rec, phi_Pmiss_angle_rec;
    Double_t Precoil_angle_rec, P_recoil_angle_rec[3], Erecoil_angle_rec, theta_Precoil_angle_rec, phi_Precoil_angle_rec;
    Double_t Pcm_angle_rec, P_cm_angle_rec[3], Ecm_angle_rec, theta_Pcm_angle_rec, phi_Pcm_angle_rec;
    
    Double_t E_cons_test;
    
    Double_t Recoil_x_NeuLand, Recoil_y_NeuLand, Recoil_z_MagnetHeight;
    Double_t P1_x, P1_y, P1_z, r1;
    Double_t P2_x, P2_y, P2_z, r2;
    
    
    TLorentzVector v1org_copy, v2org_copy, vrecoilorg_copy, vcmorg_copy, vc12_copy;
    
    //////////////////////////////////////////////////////////////////

    T->Branch("weight",&weight,"weight/D");   // Event weight based on cross-section
    T->Branch("theta_cm",&theta_cm,"theta_cm/D");
    T->Branch("Effective_E",&Effective_E,"Effective_E/D");
    T->Branch("cross_section",&cross_section,"cross_section/D");
    T->Branch("Pbeam",&Pbeam,"Pbeam/D");
    
    ////////////////////////////////////////////////////////////////// normal kinematics frame

    T->Branch("P_beam_smear_C12",&P_beam_smear_C12,"P_beam_smear_C12[3]/D");
    T->Branch("P_beam_smear_P"  ,&P_beam_smear_P  ,"P_beam_smear_P[3]/D");
    T->Branch("Pbeam_smear_C12",&Pbeam_smear_C12,"Pbeam_smear_C12/D");
    T->Branch("Pbeam_smear_P  ",&Pbeam_smear_P,"Pbeam_smear_P/D");
    T->Branch("Ebeam_smear_C12",&Ebeam_smear_C12,"Ebeam_smear_C12/D");
    T->Branch("Ebeam_smear_P  ",&Ebeam_smear_P  ,"Ebeam_smear_P/D");
    T->Branch("theta_Pbeam_smear_C12",&theta_Pbeam_smear_C12,"theta_Pbeam_smear_C12/D");
    T->Branch("theta_Pbeam_smear_P  ",&theta_Pbeam_smear_P  ,"theta_Pbeam_smear_P  /D");
    T->Branch("phi_Pbeam_smear_C12",&phi_Pbeam_smear_C12,"phi_Pbeam_smear_C12/D");
    T->Branch("phi_Pbeam_smear_P  ",&phi_Pbeam_smear_P  ,"phi_Pbeam_smear_P  /D");
    
    T->Branch("P_1",&P_1,"P_1[3]/D");
    T->Branch("P1",&P1,"P1/D");
    T->Branch("E1",&E1,"E1/D");
    T->Branch("theta_P1",&theta_P1,"theta_P1/D");
    T->Branch("phi_P1",&phi_P1,"phi_P1/D");
    
    T->Branch("P_2",&P_2,"P_2[3]/D");
    T->Branch("P2",&P2,"P2/D");
    T->Branch("E2",&E2,"E2/D");
    T->Branch("theta_P2",&theta_P2,"theta_P2/D");
    T->Branch("phi_P2",&phi_P2,"phi_P2/D");
    
    T->Branch("P_miss",&P_miss,"P_miss[3]/D");
    T->Branch("Pmiss",&Pmiss,"Pmiss/D");
    T->Branch("Emiss",&Emiss,"Emiss/D");
    T->Branch("theta_Pmiss",&theta_Pmiss,"theta_Pmiss/D");
    T->Branch("phi_Pmiss",&phi_Pmiss,"phi_Pmiss/D");
    
    T->Branch("P_recoil",&P_recoil,"P_recoil[3]/D");
    T->Branch("Precoil",&Precoil,"Precoil/D");
    T->Branch("Erecoil",&Erecoil,"Erecoil/D");
    T->Branch("theta_Precoil",&theta_Precoil,"theta_Precoil/D");   // in-plane angle
    T->Branch("phi_Precoil",&phi_Precoil,"phi_Precoil/D");         // out-of-plane angle
    
    T->Branch("P_cm",&P_cm,"P_cm[3]/D");
    T->Branch("Pcm",&Pcm,"Pcm/D");
    T->Branch("Ecm",&Ecm,"Ecm/D");
    
    T->Branch("t",&t,"t/D");
    T->Branch("u",&u,"u/D");
    T->Branch("s",&s,"s/D");
    
    ////////////////////////////////////////////////////////////////// boosted to the beam frame
    
    T->Branch("P_beam_boost",&P_beam_boost,"P_beam_boost[3]/D");
    T->Branch("Pbeam_boost",&Pbeam_boost,"Pbeam_boost/D");
    T->Branch("Ebeam_boost",&Ebeam_boost,"Ebeam_boost/D");
    T->Branch("theta_Pbeam_boost",&theta_Pbeam_boost,"theta_Pbeam_boost/D");
    T->Branch("phi_Pbeam_boost",&phi_Pbeam_boost,"phi_Pbeam_boost/D");
    
    T->Branch("P_1_boost",&P_1_boost,"P_1_boost[3]/D");
    T->Branch("P1_boost",&P1_boost,"P1_boost/D");
    T->Branch("E1_boost",&E1_boost,"E1_boost/D");
    T->Branch("theta_P1_boost",&theta_P1_boost,"theta_P1_boost/D");
    T->Branch("phi_P1_boost",&phi_P1_boost,"phi_P1_boost/D");
    
    T->Branch("P_2_boost",&P_2_boost,"P_2_boost[3]/D");
    T->Branch("P2_boost",&P2_boost,"P2_boost/D");
    T->Branch("E2_boost",&E2_boost,"E2_boost/D");
    T->Branch("theta_P2_boost",&theta_P2_boost,"theta_P2_boost/D");
    T->Branch("phi_P2_boost",&phi_P2_boost,"phi_P2_boost/D");
    
    T->Branch("P_recoil_boost",&P_recoil_boost,"P_recoil_boost[3]/D");
    T->Branch("Precoil_boost",&Precoil_boost,"Precoil_boost/D");
    T->Branch("theta_Precoil_boost",&theta_Precoil_boost,"theta_Precoil_boost/D");
    T->Branch("phi_Precoil_boost",&phi_Precoil_boost,"phi_Precoil_boost/D");
    
    T->Branch("P_cm_boost",&P_cm_boost,"P_cm_boost[3]/D");
    T->Branch("Pcm_boost",&Pcm_boost,"Pcm_boost/D");
    T->Branch("theta_Pcm_boost",&theta_Pcm_boost,"theta_Pcm_boost/D");
    T->Branch("phi_Pcm_boost",&phi_Pcm_boost,"phi_Pcm_boost/D");
    
    T->Branch("t_boost",&t_boost,"t_boost/D");
    T->Branch("u_boost",&u_boost,"u_boost/D");
    T->Branch("s_boost",&s_boost,"s_boost/D");
    
    T->Branch("P_2_angle_boost",&P_2_angle_boost,"P_2_angle_boost[3]/D");
    T->Branch("P2_angle_boost",&P2_angle_boost,"P2_angle_boost/D");
    T->Branch("E2_angle_boost",&E2_angle_boost,"E2_angle_boost/D");
    T->Branch("theta_P2_angle_boost",&theta_P2_angle_boost,"theta_P2_angle_boost/D");
    T->Branch("phi_P2_angle_boost",&phi_P2_angle_boost,"phi_P2_angle_boost/D");
    
    T->Branch("P_recoil_angle_boost",&P_recoil_angle_boost,"P_recoil_angle_boost[3]/D");
    T->Branch("Precoil_angle_boost",&Precoil_angle_boost,"Precoil_angle_boost/D");
    T->Branch("Erecoil_angle_boost",&Erecoil_angle_boost,"Erecoil_angle_boost/D");
    T->Branch("theta_Precoil_angle_boost",&theta_Precoil_angle_boost,"theta_Precoil_angle_boost/D");
    T->Branch("phi_Precoil_angle_boost",&phi_Precoil_angle_boost,"phi_Precoil_angle_boost/D");
    
    T->Branch("P_cm_angle_boost",&P_cm_angle_boost,"P_cm_angle_boost[3]/D");
    T->Branch("Pcm_angle_boost",&Pcm_angle_boost,"Pcm_angle_boost/D");
    T->Branch("Ecm_angle_boost",&Ecm_angle_boost,"Ecm_angle_boost/D");
    T->Branch("theta_Pcm_angle_boost",&theta_Pcm_angle_boost,"theta_Pcm_angle_boost/D");
    T->Branch("phi_Pcm_angle_boost",&phi_Pcm_angle_boost,"phi_Pcm_angle_boost/D");
    
    T->Branch("dEdX_B",&dEdX_B,"dEdX_B/D");
    T->Branch("dEdX_density_B",&dEdX_density_B,"dEdX_density_B/D");
    T->Branch("beta_cm_B",&beta_cm_B,"beta_cm_B/D");
    T->Branch("gamma_cm_B",&gamma_cm_B,"gamma_cm_B/D");
    
    T->Branch("dEdX_B11",&dEdX_B11,"dEdX_B11/D");
    T->Branch("dEdX_density_B11",&dEdX_density_B11,"dEdX_density_B11/D");
    T->Branch("beta_cm_B11",&beta_cm_B11,"beta_cm_B11/D");
    T->Branch("gamma_cm_B11",&gamma_cm_B11,"gamma_cm_B11/D");
    
    T->Branch("dEdX_Be",&dEdX_Be,"dEdX_Be/D");
    T->Branch("dEdX_density_Be",&dEdX_density_Be,"dEdX_density_Be/D");
    T->Branch("beta_cm_Be",&beta_cm_Be,"beta_cm_Be/D");
    T->Branch("gamma_cm_Be",&gamma_cm_Be,"gamma_cm_Be/D");
    
    T->Branch("a1",&a1,"a1/D");
    T->Branch("a2",&a2,"a2/D");
    T->Branch("a3",&a3,"a3/D");
    T->Branch("a4",&a4,"a4/D");
    T->Branch("a5",&a5,"a5/D");
    T->Branch("a6",&a6,"a6/D");
    T->Branch("a7",&a7,"a7/D");
    T->Branch("a8",&a8,"a8/D");
    T->Branch("a9",&a9,"a9/D");
    T->Branch("a10",&a10,"a10/D");
    T->Branch("a11",&a11,"a11/D");
    T->Branch("a12",&a12,"a12/D");
    
    //////////////////////////////////////////////////////////////////reconstructed after the boost
    
    T->Branch("P1_rec_boost",&P1_rec_boost,"P1_rec_boost/D");
    T->Branch("P_1_rec_boost",&P_1_rec_boost,"P_1_rec_boost[3]/D");
    T->Branch("E1_rec_boost",&E1_rec_boost,"E1_rec_boost/D");
    T->Branch("theta_P1_rec_boost",&theta_P1_rec_boost,"theta_P1_rec_boost/D");
    T->Branch("phi_P1_rec_boost",&phi_P1_rec_boost,"phi_P1_rec_boost/D");
    
    T->Branch("P2_rec_boost",&P2_rec_boost,"P2_rec_boost/D");
    T->Branch("P_2_rec_boost",&P_2_rec_boost,"P_2_rec_boost[3]/D");
    T->Branch("E2_rec_boost",&E2_rec_boost,"E2_rec_boost/D");
    T->Branch("theta_P2_rec_boost",&theta_P2_rec_boost,"theta_P2_rec_boost/D");
    T->Branch("phi_P2_rec_boost",&phi_P2_rec_boost,"phi_P2_rec_boost/D");
    
    T->Branch("P_recoil_rec_boost",&P_recoil_rec_boost,"P_recoil_rec_boost[3]/D");
    T->Branch("Precoil_rec_boost",&Precoil_rec_boost,"Precoil_rec_boost/D");
    T->Branch("theta_Precoil_rec_boost",&theta_Precoil_rec_boost,"theta_Precoil_rec_boost/D");
    T->Branch("phi_Precoil_rec_boost",&phi_Precoil_rec_boost,"phi_Precoil_rec_boost/D");
    
    T->Branch("P_cm_rec_boost",&P_cm_rec_boost,"P_cm_rec_boost[3]/D");
    T->Branch("Pcm_rec_boost",&Pcm_rec_boost,"Pcm_rec_boost/D");
    T->Branch("Ecm_rec_boost",&Ecm_rec_boost,"Ecm_rec_boost/D");
    T->Branch("theta_Pcm_rec_boost",&theta_Pcm_rec_boost,"theta_Pcm_rec_boost/D");
    T->Branch("phi_Pcm_rec_boost",&phi_Pcm_rec_boost,"phi_Pcm_rec_boost/D");
    
    T->Branch("t_rec_boost",&t_rec_boost,"t_rec_boost/D");
    T->Branch("u_rec_boost",&u_rec_boost,"u_rec_boost/D");
    T->Branch("s_rec_boost",&s_rec_boost,"s_rec_boost/D");
    
    T->Branch("t_angle_rec_boost",&t_angle_rec_boost,"t_angle_rec_boost/D");
    T->Branch("u_angle_rec_boost",&u_angle_rec_boost,"u_angle_rec_boost/D");
    T->Branch("s_angle_rec_boost",&s_angle_rec_boost,"s_angle_rec_boost/D");
    
    T->Branch("P_2_angle_rec_boost",&P_2_angle_rec_boost,"P_2_angle_rec_boost[3]/D");
    T->Branch("P2_angle_rec_boost",&P2_angle_rec_boost,"P2_angle_rec_boost/D");
    T->Branch("E2_angle_rec_boost",&E2_angle_rec_boost,"E2_angle_rec_boost/D");
    T->Branch("theta_P2_angle_rec_boost",&theta_P2_angle_rec_boost,"theta_P2_angle_rec_boost/D");
    T->Branch("phi_P2_angle_rec_boost",&phi_P2_angle_rec_boost,"phi_P2_angle_rec_boost/D");
    
    T->Branch("P_recoil_angle_rec_boost",&P_recoil_angle_rec_boost,"P_recoil_angle_rec_boost[3]/D");
    T->Branch("Precoil_angle_rec_boost",&Precoil_angle_rec_boost,"Precoil_angle_rec_boost/D");
    T->Branch("Erecoil_angle_rec_boost",&Erecoil_angle_rec_boost,"Erecoil_angle_rec_boost/D");
    T->Branch("theta_Precoil_angle_rec_boost",&theta_Precoil_angle_rec_boost,"theta_Precoil_angle_rec_boost/D");
    T->Branch("phi_Precoil_angle_rec_boost",&phi_Precoil_angle_rec_boost,"phi_Precoil_angle_rec_boost/D");
    
    T->Branch("P_cm_angle_rec_boost",&P_cm_angle_rec_boost,"P_cm_angle_rec_boost[3]/D");
    T->Branch("Pcm_angle_rec_boost",&Pcm_angle_rec_boost,"Pcm_angle_rec_boost/D");
    T->Branch("Ecm_angle_rec_boost",&Ecm_angle_rec_boost,"Ecm_angle_rec_boost/D");
    T->Branch("theta_Pcm_angle_rec_boost",&theta_Pcm_angle_rec_boost,"theta_Pcm_angle_rec_boost/D");
    T->Branch("phi_Pcm_angle_rec_boost",&phi_Pcm_angle_rec_boost,"phi_Pcm_angle_rec_boost/D");
    
    T->Branch("a1p",&a1p,"a1p/D");
    T->Branch("a2p",&a2p,"a2p/D");
    T->Branch("a3p",&a3p,"a3p/D");
    T->Branch("a4p",&a4p,"a4p/D");
    T->Branch("a5p",&a5p,"a5p/D");
    T->Branch("a6p",&a6p,"a6p/D");
    T->Branch("a7p",&a7p,"a7p/D");
    T->Branch("a8p",&a8p,"a8p/D");
    T->Branch("a9p",&a9p,"a9p/D");
    T->Branch("a10p",&a10p,"a10p/D");
    T->Branch("a11p",&a11p,"a11p/D");
    T->Branch("a12p",&a12p,"a12p/D");
    
    ////////////////////////////////////////////////////////////////// reconstructed back to the normal kinematics
    
    T->Branch("P_beam_rec",&P_beam_rec,"P_beam_rec[3]/D");
    T->Branch("Pbeam_rec",&Pbeam_rec,"Pbeam_rec/D");
    T->Branch("Ebeam_rec",&Ebeam_rec,"Ebeam_rec/D");
    T->Branch("theta_Pbeam_rec",&theta_Pbeam_rec,"theta_Pbeam_rec/D");
    T->Branch("phi_Pbeam_rec",&phi_Pbeam_rec,"phi_Pbeam_rec/D");
    
    T->Branch("P_1_rec",&P_1_rec,"P_1_rec[3]/D");
    T->Branch("P1_rec",&P1_rec,"P1_rec/D");
    T->Branch("E1_rec",&E1_rec,"E1_rec/D");
    T->Branch("theta_P1_rec",&theta_P1_rec,"theta_P1_rec/D");
    T->Branch("phi_P1_rec",&phi_P1_rec,"phi_P1_rec/D");
    
    T->Branch("P_2_rec",&P_2_rec,"P_2_rec[3]/D");
    T->Branch("P2_rec",&P2_rec,"P2_rec/D");
    T->Branch("E2_rec",&E2_rec,"E2_rec/D");
    T->Branch("theta_P2_rec",&theta_P2_rec,"theta_P2_rec/D");
    T->Branch("phi_P2_rec",&phi_P2_rec,"phi_P2_rec/D");
    
    T->Branch("P_miss_rec",&P_miss_rec,"P_miss_rec[3]/D");
    T->Branch("Pmiss_rec",&Pmiss_rec,"Pmiss_rec/D");
    T->Branch("Emiss_rec",&Emiss_rec,"Emiss_rec/D");
    T->Branch("theta_Pmiss_rec",&theta_Pmiss_rec,"theta_Pmiss_rec/D");
    T->Branch("phi_Pmiss_rec",&phi_Pmiss_rec,"phi_Pmiss_rec/D");
    
    T->Branch("P_recoil_rec",&P_recoil_rec,"P_recoil_rec[3]/D");
    T->Branch("Precoil_rec",&Precoil_rec,"Precoil_rec/D");               // momentum vector magnitude
    T->Branch("Erecoil_rec",&Erecoil_rec,"Erecoil_rec/D");
    T->Branch("theta_Precoil_rec",&theta_Precoil_rec,"theta_Precoil_rec/D");   // in-plane angle
    T->Branch("phi_Precoil_rec",&phi_Precoil_rec,"phi_Precoil_rec/D");         // out-of-plane angle
    
    T->Branch("P_cm_rec",&P_cm_rec,"P_cm_rec[3]/D");
    T->Branch("Pcm_rec",&Pcm_rec,"Pcm_rec/D");
    T->Branch("theta_Pcm_rec",&theta_Pcm_rec,"theta_Pcm_rec/D");
    T->Branch("phi_Pcm_rec",&phi_Pcm_rec,"phi_Pcm_rec/D");
    
    T->Branch("t_rec",&t_rec,"t_rec/D");
    T->Branch("u_rec",&u_rec,"u_rec/D");
    T->Branch("s_rec",&s_rec,"s_rec/D");
    
    T->Branch("t_angle_rec",&t_angle_rec,"t_angle_rec/D");
    T->Branch("u_angle_rec",&u_angle_rec,"u_angle_rec/D");
    T->Branch("s_angle_rec",&s_angle_rec,"s_angle_rec/D");
    
    T->Branch("P_2_angle_rec",&P_2_angle_rec,"P_2_angle_angle_rec[3]/D");
    T->Branch("P2_angle_rec",&P2_angle_rec,"P2_angle_rec/D");
    T->Branch("E2_angle_rec",&E2_angle_rec,"E2_angle_rec/D");
    T->Branch("theta_P2_angle_rec",&theta_P2_angle_rec,"theta_P2_angle_rec/D");
    T->Branch("phi_P2_angle_rec",&phi_P2_angle_rec,"phi_P2_angle_rec/D");
    
    T->Branch("P_miss_angle_rec",&P_miss_angle_rec,"P_miss_angle_angle_rec[3]/D");
    T->Branch("Pmiss_angle_rec",&Pmiss_angle_rec,"Pmiss_angle_rec/D");
    T->Branch("Emiss_angle_rec",&Emiss_angle_rec,"Emiss_angle_rec/D");
    T->Branch("theta_Pmiss_angle_rec",&theta_Pmiss_angle_rec,"theta_Pmiss_angle_rec/D");
    T->Branch("phi_Pmiss_angle_rec",&phi_Pmiss_angle_rec,"phi_Pmiss_angle_rec/D");
    
    T->Branch("P_recoil_angle_rec",&P_recoil_angle_rec,"P_recoil_angle_rec[3]/D");
    T->Branch("Precoil_angle_rec",&Precoil_angle_rec,"Precoil_angle_rec/D");
    T->Branch("Erecoil_angle_rec",&Erecoil_angle_rec,"Erecoil_angle_rec/D");
    T->Branch("theta_Precoil_angle_rec",&theta_Precoil_angle_rec,"theta_Precoil_angle_rec/D");
    T->Branch("phi_Precoil_angle_rec",&phi_Precoil_angle_rec,"phi_Precoil_angle_rec/D");
    
    T->Branch("P_cm_angle_rec",&P_cm_angle_rec,"P_cm_angle_rec[3]/D");
    T->Branch("Pcm_angle_rec",&Pcm_angle_rec,"Pcm_angle_rec/D");
    T->Branch("Ecm_angle_rec",&Ecm_angle_rec,"Ecm_angle_rec/D");
    T->Branch("theta_Pcm_angle_rec",&theta_Pcm_angle_rec,"theta_Pcm_angle_rec/D");
    T->Branch("phi_Pcm_angle_rec",&phi_Pcm_angle_rec,"phi_Pcm_angle_rec/D");
    
    
    T->Branch("E_cons_test",&E_cons_test,"E_cons_test/D");

    T->Branch("Recoil_x_NeuLand",&Recoil_x_NeuLand,"Recoil_x_NeuLand/D");
    T->Branch("Recoil_y_NeuLand",&Recoil_y_NeuLand,"Recoil_y_NeuLand/D");
    T->Branch("Recoil_z_MagnetHeight",&Recoil_z_MagnetHeight,"Recoil_z_MagnetHeight/D");
    
    T->Branch("P1_x",&P1_x,"P1_x/D");
    T->Branch("P1_y",&P1_y,"P1_y/D");
    T->Branch("P1_z",&P1_z,"P1_z/D");
    T->Branch("r1",&r1,"r1/D");
    
    T->Branch("P2_x",&P2_x,"P2_x/D");
    T->Branch("P2_y",&P2_y,"P2_y/D");
    T->Branch("P2_z",&P2_z,"P2_z/D");
    T->Branch("r2",&r2,"r2/D");

    //////////////////////////////////////////////////////////////////
    
    // Parametrization of the momentum distribution for high-momentum nucleons (k>k_F)
    //TF1* SRCtail = new TF1("SRCtail","1/(x**4)",0.25,0.7);
    //TF1* n_k     = new TF1("","((x<0.25)?0.8*3/((0.25)**3):0.5*0.2*2.533333/(1./(0.25)-1./5)/(x**4))       ",0,1);
    TF1* n_k_k2  = new TF1("","((x<0.25)?0.8*3/((0.25)**3):0.5*0.2*2.533333/(1./(0.25)-1./5)/(x**4))*(x**2)",0,1);
    TF1* R_2 = new TF1("","x**2",0,0.25);
    
    //
    // Events Generation
    // 1. Raffle two nucleons from a SRC pair.
    //
    const int N = 20000;
    for(int i=0; i<N; i++){
        
        if(i%(N/100)==0) cout << i/(N/100)<<"%" << endl;
        
        //for(int k=0;k<20;k++){
        for(int k=0;k<60;k++){


            //
            // Calculate UN-SMEARED proton beam in nucleus rest frame
            //
            Double_t Ebeam_tmp=TMath::Sqrt( pow(Ac12*Pbeam,2) + pow(mC12,2) );
            TLorentzVector vc12_ForReconstructions           ( TVector3(0,0,Ac12*Pbeam), Ebeam_tmp); // un-smeared Nucleus beam 4-momentum (DUBNA lab frame)
            TLorentzVector v_proton_beam_UnSMEARED( TVector3(0,0,0),                    mp    );  // Standing proton (DUBNA lab frame)
            v_proton_beam_UnSMEARED.Boost(-vc12_ForReconstructions.BoostVector()); // Boost from DUBNA to GSIto GSI frame
            
            // Test to see how to use C beam to boost from GSI to DUBNA
            TLorentzVector V_test_OR( TVector3(0,0,0),                    mC12    );  // Standing proton (DUBNA lab frame)
            V_test_OR.Boost(vc12_ForReconstructions.BoostVector());
            if( TMath::Abs( V_test_OR.Z()-Ac12*Pbeam )>0.0001 || TMath::Abs( V_test_OR.Y())>0.0001 || TMath::Abs( V_test_OR.X())>0.0001 || TMath::Abs( V_test_OR.T()-Ebeam_tmp)>0.0001){
                cout<< "Failed Or test" << endl;
            }
//            if( TMath::Abs( v_proton_beam_UnSMEARED.P() - Pbeam )>0.000001 ){
//                cout << "yaaaaa" << endl;
//            }
            
            
            //
            // Calculate SMEARED proton beam in nucleus rest frame
            //
            Pbeam_smear_C12=Ac12*(Pbeam+gRandom->Gaus(0,sigma_Pbeam*Pbeam));  // Nucleus beam momentum
            phi_Pbeam_smear_C12 = gRandom->Uniform()*2*180.; //deg
            theta_Pbeam_smear_C12 = gRandom->Gaus(0,sigma_theta_Pbeam); //deg
            P_beam_smear_C12[0]=Pbeam_smear_C12*sin(theta_Pbeam_smear_C12*TMath::DegToRad())*cos(phi_Pbeam_smear_C12*TMath::DegToRad());
            P_beam_smear_C12[1]=Pbeam_smear_C12*sin(theta_Pbeam_smear_C12*TMath::DegToRad())*sin(phi_Pbeam_smear_C12*TMath::DegToRad());
            P_beam_smear_C12[2]=Pbeam_smear_C12*cos(theta_Pbeam_smear_C12*TMath::DegToRad());
            Ebeam_smear_C12=TMath::Sqrt( pow(Pbeam_smear_C12,2) + pow(mC12,2) );
            
            TLorentzVector vc12           ( TVector3(P_beam_smear_C12[0],P_beam_smear_C12[1],P_beam_smear_C12[2]), Ebeam_smear_C12); // smeared Nucleus beam 4-momentum (DUBNA lab frame)
            TLorentzVector vc12_tmpForTest( TVector3(P_beam_smear_C12[0],P_beam_smear_C12[1],P_beam_smear_C12[2]), Ebeam_smear_C12); // smeared Nucleus beam 4-momentum (DUBNA lab frame)
            TLorentzVector v_proton_beam_SMEARED( TVector3(0,0,0),                    mp    );  // Standing proton (DUBNA lab frame)
            
            v_proton_beam_SMEARED.Boost(-vc12.BoostVector()); // Boost from DUBNA lab frame to nucleus rest frame
            P_beam_smear_P[0]=v_proton_beam_SMEARED.X();
            P_beam_smear_P[1]=v_proton_beam_SMEARED.Y();
            P_beam_smear_P[2]=v_proton_beam_SMEARED.Z();
            theta_Pbeam_smear_P=v_proton_beam_SMEARED.Theta()*TMath::RadToDeg();
            phi_Pbeam_smear_P=v_proton_beam_SMEARED.Phi()*TMath::RadToDeg();
            Pbeam_smear_P=v_proton_beam_SMEARED.P();
            Ebeam_smear_P=v_proton_beam_SMEARED.T();

        
            vc12_tmpForTest.Boost(-vc12.BoostVector()); // test that the nucleus boos to a standing nucleus
            if( TMath::Abs(vc12_tmpForTest.X())>0.0001 ||TMath::Abs(vc12_tmpForTest.Y())>0.0001 ||TMath::Abs(vc12_tmpForTest.Z())>0.0001 ||TMath::Abs(vc12_tmpForTest.T() - mC12)>0.0001  ){
                cout << "issue in boost!\t" << TMath::Abs(vc12_tmpForTest.X())<<"\t"<< TMath::Abs(vc12_tmpForTest.Y())<<"\t"<< TMath::Abs(vc12_tmpForTest.Z())<<"\t"<< TMath::Abs(vc12_tmpForTest.T()) << endl;
            }

            

            
            //
            //theta_cm = 80.5+k;
            //
            theta_cm = 60.5+k;
            
            
            //
            // Raffle 3-momentum vector for 1st nucleon in the pair.
            //
            Pmiss=n_k_k2->GetRandom();                              // raffle vector magnitude
            gRandom->Sphere(P_miss[0],P_miss[1],P_miss[2], Pmiss);   // smear vector direction
            theta_Pmiss = TMath::ACos(P_miss[2]/Pmiss)*TMath::RadToDeg();      // Calculate theta [degrees]
            phi_Pmiss = TMath::ATan2(P_miss[1],P_miss[0])*TMath::RadToDeg();       // Calculate phi   [degrees]
            
            

            //
            // Calculate 3-momentum vector for the 2nd nucleon in the pair.
            //
            if(Pmiss < 0.25){
                
                Precoil=R_2->GetRandom();                              // raffle vector magnitude
                gRandom->Sphere(P_recoil[0],P_recoil[1],P_recoil[2], Precoil);   // smear vector direction
                theta_Precoil = TMath::ACos(P_recoil[2]/Precoil)*TMath::RadToDeg();      // Calculate theta [degrees]
                phi_Precoil = TMath::ATan2(P_recoil[1],P_recoil[0])*TMath::RadToDeg();       // Calculate phi   [degrees]
                Erecoil = TMath::Sqrt( pow(Precoil,2) + pow(mn,2) );

                P_cm[0] = P_recoil[0] + P_miss[0];
                P_cm[1] = P_recoil[1] + P_miss[1];
                P_cm[2] = P_recoil[2] + P_miss[2];
                Pcm = TMath::Sqrt( pow(P_cm[0],2) + pow(P_cm[1],2) + pow(P_cm[2],2) );
                Ecm=TMath::Sqrt( pow(Pcm,2) + pow(mB10,2) );

            }else{
                //
                // Raffle SRC-pair C.M. momentum.
                //
                P_cm[0] = gRandom->Gaus(0,0.14);
                P_cm[1] = gRandom->Gaus(0,0.14);
                P_cm[2] = gRandom->Gaus(0,0.14);
                Pcm = TMath::Sqrt( pow(P_cm[0],2) + pow(P_cm[1],2) + pow(P_cm[2],2) );
                Ecm=TMath::Sqrt( pow(Pcm,2) + pow(mB10,2) );

                P_recoil[0] = P_cm[0] - P_miss[0];
                P_recoil[1] = P_cm[1] - P_miss[1];
                P_recoil[2] = P_cm[2] - P_miss[2];
                Precoil = TMath::Sqrt(pow(P_recoil[0],2)+pow(P_recoil[1],2)+pow(P_recoil[2],2));
            	theta_Precoil = TMath::ACos(P_recoil[2]/Precoil)*TMath::RadToDeg();
                phi_Precoil = TMath::ATan2(P_recoil[1],P_recoil[0])*TMath::RadToDeg();
                Erecoil = TMath::Sqrt( pow(Precoil,2) + pow(mn,2) );
            }
            
            // NN-SRC pair energy
            Double_t Ed= mC12 - Ecm;
            
            
            
            //
            // 4-momenta for smeared beam, Pmiss, Precoil, Pcm (= P_A-2)
            //
            // beam already defined by:  v_proton_beam_SMEARED
            TLorentzVector v2org     ( TVector3(P_miss  [0],P_miss  [1],P_miss  [2]), Ed-sqrt(pow(mn,2)+pow(Precoil,2)) ); // Missing Momentum ('initial' nucleon in the nucleus)
            TLorentzVector vrecoilorg( TVector3(P_recoil[0],P_recoil[1],P_recoil[2]),    sqrt(pow(mn,2)+pow(Precoil,2)) ); // Recoil Momentum (nucleon in the nucleus)
            TLorentzVector vcmorg( TVector3(-P_cm[0],-P_cm[1],-P_cm[2]), Ecm ); // Recoil Nucleus A-2 initial is 12C final 10B or 10Be
            TLorentzVector vcmorg_TESTTEST( TVector3(-P_cm[0],-P_cm[1],-P_cm[2]), Ecm ); // Recoil Nucleus A-2 initial is 12C final 10B or 10Be
        
            
            

            //
            // calculate cross-section by boosting to the ||nucleon rest frame|| and extracting the 'effective beam Momentum'.
            //
            TLorentzVector m1 = v2org; //v_miss for Lorentz Boost
            m1.SetT(  sqrt(  pow(mp,2)+pow(v2org.P(),2)  )  ); // make Pmiss on-shell for rate claculation (?????)
            v_proton_beam_SMEARED.Boost(-m1.BoostVector());
            Effective_E = v_proton_beam_SMEARED.P(); // Effective momentum
            if (Effective_E<2){// calculate the 'effective' cross section
                cross_section = f11->Eval( Effective_E );
            }else if (Effective_E>=2 && Effective_E<8){
                cross_section = f22->Eval( Effective_E );
            } else if (Effective_E>=8){
                cross_section = f33->Eval( Effective_E );
            }
            cross_section = cross_section * (pow(1-pow(cos(theta_cm*TMath::DegToRad()),2),-4*0.919));
            weight = cross_section * mbtob * btocmsq * (2*TMath::Pi()*(cos((theta_cm-0.5)*TMath::DegToRad())-cos((theta_cm+0.5)*TMath::DegToRad()))) * flux * target_thickness * transparency * transparency *  neutron_efficiency * hrstosec * hrs * time * duty_cycle * factor_solid_angle / N;
            v_proton_beam_SMEARED.Boost(m1.BoostVector());

            
            

            
            //
            // Calculate 90 degrees scattering for a beam proton and the nucleon from the SRC pair (missing momentum)
            //
            
            // boost to the c.m. frame. of the beam + SRC nucleon
            TLorentzVector m = (v_proton_beam_SMEARED + v2org); //C.M. for Lorentz Boost
            
            v_proton_beam_SMEARED.Boost(-m.BoostVector());
            v2org.Boost(-m.BoostVector());
            if( TMath::Abs(v_proton_beam_SMEARED.X() + v2org.X()) > 0.00001 || TMath::Abs(v_proton_beam_SMEARED.Y() + v2org.Y()) > 0.00001 || TMath::Abs(v_proton_beam_SMEARED.Z() + v2org.Z()) > 0.00001 ){
                cout << "In the Pbeam-Pm c.m. frame. momenta not equal!!" << endl;
            }

            // Define the angles of the scattering in the c.m.
            Double_t theta_ll=theta_cm*TMath::DegToRad();              // 90 degrees in the c.m.
            Double_t phi_ll = gRandom->Uniform()*2*TMath::Pi();  // for 90 degrees scattering, phi is arbitraty
            
            
            Double_t v3_Energy_tmp = (v_proton_beam_SMEARED.T() + v2org.T())/2;
            Double_t v3_Momentum_tmp = sqrt(pow(v3_Energy_tmp,2) - pow(mp,2));
            TLorentzVector v3( TVector3(   v3_Momentum_tmp*cos(phi_ll)*sin(theta_ll),   v3_Momentum_tmp*sin(phi_ll)*sin(theta_ll),  v3_Momentum_tmp*cos(theta_ll)  ),                    v3_Energy_tmp    ); // Beam Proton smeared
            TLorentzVector v4( TVector3(  -v3_Momentum_tmp*cos(phi_ll)*sin(theta_ll),  -v3_Momentum_tmp*sin(phi_ll)*sin(theta_ll), -v3_Momentum_tmp*cos(theta_ll)  ),                    v3_Energy_tmp    ); // Beam Proton smeared

            
            v3.RotateY(v_proton_beam_SMEARED.Theta());
            v3.RotateZ(v_proton_beam_SMEARED.Phi());
            v4.RotateY(v_proton_beam_SMEARED.Theta());
            v4.RotateZ(v_proton_beam_SMEARED.Phi());
            if( TMath::Abs(v_proton_beam_SMEARED.Angle(v3.Vect())-theta_ll)>0.000001 ){
                cout << "issue with vector aligment after scattering:\t" << v_proton_beam_SMEARED.Angle(v3.Vect())*TMath::RadToDeg() << "\t" << theta_ll*TMath::RadToDeg()  << endl;
            }
            
            
            // Boost back to the lab
            v3.Boost(m.BoostVector());
            v4.Boost(m.BoostVector());
            v_proton_beam_SMEARED.Boost(m.BoostVector());
            v2org.Boost(m.BoostVector());

            if( TMath::Abs(v3.T() - sqrt( pow(mp,2) + pow(v3.X(),2)+pow(v3.Y(),2)+pow(v3.Z(),2) ) )>0.000001 ){
                cout << "After boost from Pbeam-Pm c.m. frame: protons (v3) don't have E = sqrt(m2+p2)" << endl;
            }
            if( TMath::Abs(v4.T() - sqrt( pow(mp,2) + pow(v4.X(),2)+pow(v4.Y(),2)+pow(v4.Z(),2) ) )>0.000001 ){
                cout << "After boost from Pbeam-Pm c.m. frame: protons (v4) don't have E = sqrt(m2+p2)" << endl;
            }
            
            
            
            
            
            
            
            
            //
            // keep all variables in the GSI frame
            //
            P_1[0] = v3.X();
            P_1[1] = v3.Y();
            P_1[2] = v3.Z();
            P1 = v3.P();
            E1 = v3.T();
            theta_P1 = v3.Theta()*TMath::RadToDeg();
            phi_P1 = v3.Phi()*TMath::RadToDeg();

            P_2[0] = v4.X();
            P_2[1] = v4.Y();
            P_2[2] = v4.Z();
            P2 = v4.P();
            E2 = v4.T();
            theta_P2 = v4.Theta()*TMath::RadToDeg();
            phi_P2 = v4.Phi()*TMath::RadToDeg();

            
            //
            // Calculate t, u, and s
            //
            s= (v3+v4)*(v3+v4);
            u= (v_proton_beam_SMEARED-v4)*(v_proton_beam_SMEARED-v4);
            t= (v_proton_beam_SMEARED-v3)*(v_proton_beam_SMEARED-v3);
            
            
            //////////////////////////////////////////////////////////////////boost to the proton beam rest frame
            
            
            
            //
            // Make copies of original 4-momenta so we use the copis in the boosts (i.e. keep original un touched)
            //
            v1org_copy=v_proton_beam_SMEARED; // beam
            vrecoilorg_copy=vrecoilorg; //recoil
            vcmorg_copy=vcmorg; // A-2
            TLorentzVector v1final = v3; //p1
            TLorentzVector v2final = v4; //p2

        
        
            //
            // boost from GSI frame to DUBNA frame
            //
            v1org_copy.Boost     (vc12.BoostVector());
            vrecoilorg_copy.Boost(vc12.BoostVector());
            vcmorg_copy.Boost    (vc12.BoostVector());
            v1final.Boost        (vc12.BoostVector());
            v2final.Boost        (vc12.BoostVector());
            if( TMath::Abs(v1org_copy.X())>0.000001 || TMath::Abs(v1org_copy.Y())>0.000001 || TMath::Abs(v1org_copy.Z())>0.000001  || TMath::Abs(v1org_copy.T()-mp)>0.00001  ){
                cout << "problem: v1 after boos from GSI frame to DUBNA using v_proton_beam_SMEARED is not standing" << endl;
            }
            TLorentzVector v_StandingC_test( TVector3(0,0,0),                    mC12    );  // Standing proton (DUBNA lab frame)
            v_StandingC_test.Boost     (vc12.BoostVector());
            if( TMath::Abs(v_StandingC_test.X()-P_beam_smear_C12[0])>0.000001 || TMath::Abs(v_StandingC_test.Y()-P_beam_smear_C12[1])>0.000001 || TMath::Abs(v_StandingC_test.Z()-P_beam_smear_C12[2])>0.000001  || TMath::Abs(v_StandingC_test.T()-sqrt(pow(Pbeam_smear_C12,2)+pow(mC12,2)))>0.00001  ){
                cout << "boot issue in GSI -> Dubna. Standing Carbon isn't moving right!\t" << v_StandingC_test.X()<<"\t"<<v_StandingC_test.Y()<<"\t"<<v_StandingC_test.Z()<<"\t"<<P_beam_smear_C12[2]<<"\t"<< v_StandingC_test.T()<<"\t"<< sqrt(pow(Pbeam_smear_C12,2)+pow(mC12,2))<< endl;
            }

            P_beam_boost[0] = v1org_copy.X();
            P_beam_boost[1] = v1org_copy.Y();
            P_beam_boost[2] = v1org_copy.Z();
            Pbeam_boost = v1org_copy.P();
            Ebeam_boost = v1org_copy.T();
            theta_Pbeam_boost = v1org_copy.Theta();
            phi_Pbeam_boost = v1org_copy.Phi();
            // another sanity check that proton beam is standing after boost from GSI to DUBNA
            if( TMath::Abs(P_beam_boost[0])>0.0001 || TMath::Abs(P_beam_boost[1])>0.0001 || TMath::Abs(P_beam_boost[2])>0.0001 || TMath::Abs(Ebeam_boost-mp)>0.0001 ){
                cout << "P_beam_boost not zero" << endl;
            }
            
            P_1_boost[0] = v1final.X();
            P_1_boost[1] = v1final.Y();
            P_1_boost[2] = v1final.Z();
            P1_boost = v1final.P();
            E1_boost= v1final.T();
            theta_P1_boost = v1final.Theta();
            phi_P1_boost = v1final.Phi();
            if( TMath::Abs( theta_P1_boost - TMath::ACos(P_1_boost[2]/P1_boost) )>0.0001  ||
                TMath::Abs( phi_P1_boost   - TMath::ATan2(P_1_boost[1],P_1_boost[0]) )>0.0001 ){
                cout << "angles from root and from our clacualation not consistent! (1)" << endl;
            }
            
            P_2_boost[0] = v2final.X();
            P_2_boost[1] = v2final.Y();
            P_2_boost[2] = v2final.Z();
            P2_boost = TMath::Sqrt(pow(P_2_boost[0],2) + pow(P_2_boost[1],2) + pow(P_2_boost[2],2));
            E2_boost= sqrt( pow(P2_boost,2) + pow(mp,2));
            theta_P2_boost = v2final.Theta();
            phi_P2_boost = v2final.Phi();
            if( TMath::Abs( theta_P2_boost - TMath::ACos(P_2_boost[2]/P2_boost) )>0.0001  ||
                TMath::Abs( phi_P2_boost   - TMath::ATan2(P_2_boost[1],P_2_boost[0]) )>0.0001 ){
                cout << "angles from root and from our clacualation not consistent! (2)" << endl;
            }
            
            P_recoil_boost[0] = vrecoilorg_copy.X();
            P_recoil_boost[1] = vrecoilorg_copy.Y();
            P_recoil_boost[2] = vrecoilorg_copy.Z();
            Precoil_boost = vrecoilorg_copy.P();
            Erecoil_boost= vrecoilorg_copy.T();
            theta_Precoil_boost = vrecoilorg_copy.Theta();
            phi_Precoil_boost = vrecoilorg_copy.Phi();
            if( TMath::Abs( theta_Precoil_boost - TMath::ACos(P_recoil_boost[2]/Precoil_boost) )>0.0001  ||
               TMath::Abs( phi_Precoil_boost   - TMath::ATan2(P_recoil_boost[1],P_recoil_boost[0]) )>0.0001 ){
                cout << "angles from root and from our clacualation not consistent! (3)" << endl;
            }
            
            P_cm_boost[0] = vcmorg_copy.X();
            P_cm_boost[1] = vcmorg_copy.Y();
            P_cm_boost[2] = vcmorg_copy.Z();
            Pcm_boost = vcmorg_copy.P();
            Ecm_boost = vcmorg_copy.T();
            Ecm_11_boost = TMath::Sqrt( pow(Pcm_boost,2) + pow(AB11*mp,2) );// for the special case of 11B
            theta_Pcm_boost = vcmorg_copy.Theta();
            phi_Pcm_boost = vcmorg_copy.Phi();
            if( TMath::Abs( theta_Pcm_boost - TMath::ACos(P_cm_boost[2]/Pcm_boost) )>0.0001  ||
               TMath::Abs( phi_Pcm_boost   - TMath::ATan2(P_cm_boost[1],P_cm_boost[0]) )>0.0001 ){
                cout << "angles from root and from our clacualation not consistent! (4)" << endl;
            }
            if(TMath::Abs( P_cm_boost[0] - Pcm_boost*cos(phi_Pcm_boost)*sin(theta_Pcm_boost) )>0.00001 ||
               TMath::Abs( P_cm_boost[1] - Pcm_boost*sin(phi_Pcm_boost)*sin(theta_Pcm_boost) )>0.00001  ||
               TMath::Abs( P_cm_boost[2] - Pcm_boost*cos(theta_Pcm_boost) )>0.00001 ){
                cout << "cm vector: definition of components using angles incorrect!\t"<< P_cm_boost[0] <<"\t"<< Pcm_boost*cos(phi_Pcm_boost)*sin(theta_Pcm_boost)<<"\t"<<P_cm_boost[1] <<"\t"<< Pcm_boost*sin(phi_Pcm_boost)*sin(theta_Pcm_boost)<<"\t"<<P_cm_boost[2] <<"\t"<< Pcm_boost*cos(theta_Pcm_boost)<< endl;
            }

        
            s_boost= (v1final+v2final)*(v1final+v2final);
            u_boost= (v1org_copy-v2final)*(v1org_copy-v2final);
            t_boost= (v1org_copy-v1final)*(v1org_copy-v1final);
            
            
            
            
            
            
            
            
            
            //
            // Calculate P2, Precoil and Pcm(=P_A-2) in the DUBNA frame using angles and momentum conservation
            //
            a1= TMath::Sin(theta_P1_boost)      * TMath::Cos(phi_P1_boost);
            a2= TMath::Sin(theta_P2_boost)      * TMath::Cos(phi_P2_boost);
            a3= TMath::Sin(theta_Precoil_boost) * TMath::Cos(phi_Precoil_boost);
            a4= TMath::Sin(theta_Pcm_boost)     * TMath::Cos(phi_Pcm_boost);
            
            a5= TMath::Sin(theta_P1_boost)      * TMath::Sin(phi_P1_boost);
            a6= TMath::Sin(theta_P2_boost)      * TMath::Sin(phi_P2_boost);
            a7= TMath::Sin(theta_Precoil_boost) * TMath::Sin(phi_Precoil_boost);
            a8= TMath::Sin(theta_Pcm_boost)     * TMath::Sin(phi_Pcm_boost);
            
            a9 = TMath::Cos(theta_P1_boost);
            a10= TMath::Cos(theta_P2_boost);
            a11= TMath::Cos(theta_Precoil_boost);
            a12= TMath::Cos(theta_Pcm_boost);
            
            P2_angle_boost=
            (a12*(a3*a5*P1_boost - a1*a7*P1_boost + a7*P_beam_smear_C12[0] - a3*P_beam_smear_C12[1]) + a11*(-a4*a5*P1_boost + a1*a8*P1_boost - a8*P_beam_smear_C12[0] + a4*P_beam_smear_C12[1]) +
             (a4*a7 -a3*a8)*(a9*P1_boost - P_beam_smear_C12[2]))/
            (-a12*a3*a6 + a11*a4*a6 + a12*a2*a7 - a10*a4*a7 - a11*a2*a8 + a10*a3*a8);
            
            P_2_angle_boost[0]=P2_angle_boost*a2;
            P_2_angle_boost[1]=P2_angle_boost*a6;
            P_2_angle_boost[2]=P2_angle_boost*a10;
            E2_angle_boost= TMath::Sqrt( pow(P2_angle_boost,2) + pow(mp,2));
            theta_P2_angle_boost = TMath::ACos(P_2_angle_boost[2]/P2_angle_boost) * TMath::RadToDeg();
            phi_P2_angle_boost = TMath::ATan2(P_2_angle_boost[1],P_2_angle_boost[0]) * TMath::RadToDeg();
            TLorentzVector v_p2_angle       ( TVector3(P_2_angle_boost[0],P_2_angle_boost[1],P_2_angle_boost[2]), E2_angle_boost);

            
            Precoil_angle_boost=
            (a12*(a2*a5*P1_boost - a1*a6*P1_boost + a6*P_beam_smear_C12[0] - a2*P_beam_smear_C12[1]) + a10*(-a4*a5*P1_boost + a1*a8*P1_boost - a8*P_beam_smear_C12[0] + a4*P_beam_smear_C12[1]) +
             (a4*a6 - a2*a8)*(a9*P1_boost - P_beam_smear_C12[2]))/
            (a12*a3*a6 - a11*a4*a6 - a12*a2*a7 + a10*a4*a7 + a11*a2*a8 - a10*a3*a8);
            
            P_recoil_angle_boost[0]=Precoil_angle_boost*a3;
            P_recoil_angle_boost[1]=Precoil_angle_boost*a7;
            P_recoil_angle_boost[2]=Precoil_angle_boost*a11;
            Erecoil_angle_boost= sqrt( pow(Precoil_angle_boost,2) + pow(mn,2));
            theta_Precoil_angle_boost = TMath::ACos(P_recoil_angle_boost[2]/Precoil_angle_boost) * TMath::RadToDeg();
            phi_Precoil_angle_boost = TMath::ATan2(P_recoil_angle_boost[1],P_recoil_angle_boost[0]) * TMath::RadToDeg();
            TLorentzVector v_precoil_angle       ( TVector3(P_recoil_angle_boost[0],P_recoil_angle_boost[1],P_recoil_angle_boost[2]), Erecoil_angle_boost);

            
            Pcm_angle_boost=
            (a11*(a2*a5*P1_boost - a1*a6*P1_boost + a6*P_beam_smear_C12[0] - a2*P_beam_smear_C12[1]) + a10*(-a3*a5*P1_boost + a1*a7*P1_boost - a7*P_beam_smear_C12[0] + a3*P_beam_smear_C12[1]) +
             (a3*a6 - a2*a7)*(a9*P1_boost - P_beam_smear_C12[2]))/
            (-a12*a3*a6 + a11*a4*a6 + a12*a2*a7 - a10*a4*a7 - a11*a2*a8 + a10*a3*a8);
            
            P_cm_angle_boost[0]=Pcm_angle_boost*a4;
            P_cm_angle_boost[1]=Pcm_angle_boost*a8;
            P_cm_angle_boost[2]=Pcm_angle_boost*a12;
            Ecm_angle_boost= sqrt( pow(Pcm_angle_boost,2) + pow(mB10,2));
            theta_Pcm_angle_boost = TMath::ACos(P_cm_angle_boost[2]/Pcm_angle_boost) * TMath::RadToDeg();
            phi_Pcm_angle_boost = TMath::ATan2(P_cm_angle_boost[1],P_cm_angle_boost[0]) * TMath::RadToDeg();
            TLorentzVector v_cm_angle       ( TVector3(P_cm_angle_boost[0],P_cm_angle_boost[1],P_cm_angle_boost[2]), Ecm_angle_boost);

            
            E_cons_test = (vc12-v1final-v_p2_angle-v_precoil_angle-v_cm_angle)*(vc12-v1final-v_p2_angle-v_precoil_angle-v_cm_angle) - mp*mp;
            
            
            //
            //  Calculation of betta and gamma for dE/dx calculation for NeuLAND
            //
            beta_cm_B = Pcm_boost/Ecm_boost;
            gamma_cm_B = 1/sqrt(1-pow(beta_cm_B,2));
            dEdX_B= factor*atomic_number/atomic_weight*pow(charge_B/beta_cm_B,2)*(TMath::Log(pow(2*me/I,2)*pow(beta_cm_B*gamma_cm_B,4))-2*pow(beta_cm_B,2));
            dEdX_density_B= factor*density*atomic_number/atomic_weight*pow(charge_B/beta_cm_B,2)*(TMath::Log(pow(2*me/I,2)*pow(beta_cm_B*gamma_cm_B,4))-2*pow(beta_cm_B,2));
            
            beta_cm_B11 = Pcm_boost/Ecm_11_boost;
            gamma_cm_B11 = 1/sqrt(1-pow(beta_cm_B11,2));
            dEdX_B11= factor*atomic_number/atomic_weight*pow(charge_B/beta_cm_B11,2)*(TMath::Log(pow(2*me/I,2)*pow(beta_cm_B11*gamma_cm_B11,4))-2*pow(beta_cm_B11,2));
            dEdX_density_B11= factor*density*atomic_number/atomic_weight*pow(charge_B/beta_cm_B11,2)*(TMath::Log(pow(2*me/I,2)*pow(beta_cm_B11*gamma_cm_B11,4))-2*pow(beta_cm_B11,2));
            
            beta_cm_Be = Pcm_boost/Ecm_boost;
            gamma_cm_Be = 1/sqrt(1-pow(beta_cm_Be,2));
            dEdX_Be= factor*atomic_number/atomic_weight*pow(charge_Be/beta_cm_Be,2)*(TMath::Log(pow(2*me/I,2)*pow(beta_cm_Be*gamma_cm_Be,4))-2*pow(beta_cm_Be,2));
            dEdX_density_Be= factor*density*atomic_number/atomic_weight*pow(charge_Be/beta_cm_Be,2)*(TMath::Log(pow(2*me/I,2)*pow(beta_cm_Be*gamma_cm_Be,4))-2*pow(beta_cm_Be,2));
            
            
            
            
            //
            // Change from radians to degrees
            //
            theta_P1_boost      *= TMath::RadToDeg();
            phi_P1_boost        *= TMath::RadToDeg();
            theta_P2_boost      *= TMath::RadToDeg();
            phi_P2_boost        *= TMath::RadToDeg();
            theta_Precoil_boost *= TMath::RadToDeg();
            phi_Precoil_boost   *= TMath::RadToDeg();
            theta_Pcm_boost     *= TMath::RadToDeg();
            phi_Pcm_boost       *= TMath::RadToDeg();
            

            
            
            
            
            
            //
            //  Smearing DUBNA
            //
            time_P1 = (E1_boost/P1_boost)*3.3*dist1; //1/c in nsec/m
            time_P2 = (E2_boost/P2_boost)*3.3*dist2; //1/c in nsec/m
            time_Precoil = (Erecoil_boost/Precoil_boost)*3.3*disttoZDC; //1/c in nsec/m
            
            time_P1_rec = gRandom->Gaus(time_P1,sigma_timeRPC);
            time_P2_rec = gRandom->Gaus(time_P2,sigma_timeRPC);
            time_Precoil_rec = gRandom->Gaus(time_Precoil,sigma_timeZDC);
            
            theta_P1_rec_boost=theta_P1_boost+gRandom->Gaus(0,sigma_theta_P1);
            phi_P1_rec_boost=phi_P1_boost+gRandom->Gaus(0,sigma_phi_P1);
            P1_rec_boost= sqrt( pow(mp,2) / (pow(time_P1_rec/(3.3*dist1),2) - 1) );
            E1_rec_boost= sqrt( pow(P1_rec_boost,2) + pow(mp,2));
            P_1_rec_boost[0]=P1_rec_boost*cos(phi_P1_rec_boost*TMath::DegToRad())*sin(theta_P1_rec_boost*TMath::DegToRad());
            P_1_rec_boost[1]=P1_rec_boost*sin(phi_P1_rec_boost*TMath::DegToRad())*sin(theta_P1_rec_boost*TMath::DegToRad());
            P_1_rec_boost[2]=P1_rec_boost*cos(theta_P1_rec_boost*TMath::DegToRad());
            
            theta_P2_rec_boost=theta_P2_boost+gRandom->Gaus(0,sigma_theta_P2);
            phi_P2_rec_boost=phi_P2_boost+gRandom->Gaus(0,sigma_phi_P2);
            P2_rec_boost= sqrt( pow(mp,2) / (pow(time_P2_rec/(3.3*dist2),2) - 1) );
            E2_rec_boost= sqrt( pow(P2_rec_boost,2) + pow(mp,2));
            P_2_rec_boost[0]=P2_rec_boost*cos(phi_P2_rec_boost*TMath::DegToRad())*sin(theta_P2_rec_boost*TMath::DegToRad());
            P_2_rec_boost[1]=P2_rec_boost*sin(phi_P2_rec_boost*TMath::DegToRad())*sin(theta_P2_rec_boost*TMath::DegToRad());
            P_2_rec_boost[2]=P2_rec_boost*cos(theta_P2_rec_boost*TMath::DegToRad());
            
            theta_Precoil_rec_boost=theta_Precoil_boost+gRandom->Gaus(0,sigma_theta_Precoil);
            phi_Precoil_rec_boost=phi_Precoil_boost+gRandom->Gaus(0,sigma_phi_Precoil);
            Precoil_rec_boost= sqrt( pow(mn,2) / (pow(time_Precoil_rec/(3.3*disttoZDC),2) - 1) );
            Erecoil_rec_boost= sqrt( pow(Precoil_rec_boost,2) + pow(mn,2));
            P_recoil_rec_boost[0]=Precoil_rec_boost*cos(phi_Precoil_rec_boost*TMath::DegToRad())*sin(theta_Precoil_rec_boost*TMath::DegToRad());
            P_recoil_rec_boost[1]=Precoil_rec_boost*sin(phi_Precoil_rec_boost*TMath::DegToRad())*sin(theta_Precoil_rec_boost*TMath::DegToRad());
            P_recoil_rec_boost[2]=Precoil_rec_boost*cos(theta_Precoil_rec_boost*TMath::DegToRad());
            
            
            theta_Pcm_rec_boost=theta_Pcm_boost+gRandom->Gaus(0,sigma_theta_Pcm);
            phi_Pcm_rec_boost=phi_Pcm_boost+gRandom->Gaus(0,sigma_phi_Pcm);
            Pcm_rec_boost=Pcm_boost+gRandom->Gaus(0,Pcm_boost*sigma_Pcm);
            Ecm_rec_boost=sqrt( pow(Pcm_rec_boost,2) + pow(mB10,2));
            P_cm_rec_boost[0]=Pcm_rec_boost*cos(phi_Pcm_rec_boost*TMath::DegToRad())*sin(theta_Pcm_rec_boost*TMath::DegToRad());
            P_cm_rec_boost[1]=Pcm_rec_boost*sin(phi_Pcm_rec_boost*TMath::DegToRad())*sin(theta_Pcm_rec_boost*TMath::DegToRad());
            P_cm_rec_boost[2]=Pcm_rec_boost*cos(theta_Pcm_rec_boost*TMath::DegToRad());



            //
            // Reconstruction based on smeared angles
            //
            a1p= TMath::Sin(theta_P1_rec_boost*TMath::DegToRad())*TMath::Cos(phi_P1_rec_boost*TMath::DegToRad());
            a2p= TMath::Sin(theta_P2_rec_boost*TMath::DegToRad())*TMath::Cos(phi_P2_rec_boost*TMath::DegToRad());
            a3p= TMath::Sin(theta_Precoil_rec_boost*TMath::DegToRad())*TMath::Cos(phi_Precoil_rec_boost*TMath::DegToRad());
            a4p= TMath::Sin(theta_Pcm_rec_boost*TMath::DegToRad())*TMath::Cos(phi_Pcm_rec_boost*TMath::DegToRad());
            
            a5p= TMath::Sin(theta_P1_rec_boost*TMath::DegToRad())*TMath::Sin(phi_P1_rec_boost*TMath::DegToRad());
            a6p= TMath::Sin(theta_P2_rec_boost*TMath::DegToRad())*TMath::Sin(phi_P2_rec_boost*TMath::DegToRad());
            a7p= TMath::Sin(theta_Precoil_rec_boost*TMath::DegToRad())*TMath::Sin(phi_Precoil_rec_boost*TMath::DegToRad());
            a8p= TMath::Sin(theta_Pcm_rec_boost*TMath::DegToRad())*TMath::Sin(phi_Pcm_rec_boost*TMath::DegToRad());
            
            a9p= TMath::Cos(theta_P1_rec_boost*TMath::DegToRad());
            a10p= TMath::Cos(theta_P2_rec_boost*TMath::DegToRad());
            a11p= TMath::Cos(theta_Precoil_rec_boost*TMath::DegToRad());
            a12p= TMath::Cos(theta_Pcm_rec_boost*TMath::DegToRad());
            
            
            
            
            P2_angle_rec_boost =
            (a12p*(a3p*a5p - a1p*a7p)*P1_rec_boost + a11p*(-a4p*a5p*P1_rec_boost + a1p*a8p*P1_rec_boost) + (a4p*a7p - a3p*a8p)*(a9p*P1_rec_boost - Ac12*Pbeam))/
            (-a12p*a3p*a6p + a11p*a4p*a6p + a12p*a2p*a7p - a10p*a4p*a7p - a11p*a2p*a8p + a10p*a3p*a8p);
            
            P_2_angle_rec_boost[0]=P2_angle_rec_boost*a2p;
            P_2_angle_rec_boost[1]=P2_angle_rec_boost*a6p;
            P_2_angle_rec_boost[2]=P2_angle_rec_boost*a10p;
            E2_angle_rec_boost= TMath::Sqrt( pow(P2_angle_rec_boost,2) + pow(mp,2));
            theta_P2_angle_rec_boost = TMath::ACos(P_2_angle_rec_boost[2]/P2_angle_rec_boost) * TMath::RadToDeg();
            phi_P2_angle_rec_boost = TMath::ATan2(P_2_angle_rec_boost[1],P_2_angle_rec_boost[0]) * TMath::RadToDeg();
            
            if(TMath::Abs(phi_P2_boost-phi_P2_angle_rec_boost) > 30){
                if( phi_P2_angle_rec_boost > 0){
                    phi_P2_angle_rec_boost -= 180;
                }else{
                    phi_P2_angle_rec_boost += 180;
                }
            }
            
            
            
            Precoil_angle_rec_boost=
            (a12p*(a2p*a5p - a1p*a6p)*P1_rec_boost + a10p*(-a4p*a5p*P1_rec_boost + a1p*a8p*P1_rec_boost) + (a4p*a6p - a2p*a8p)*(a9p*P1_rec_boost - Ac12*Pbeam))/
            (a12p*a3p*a6p - a11p*a4p*a6p - a12p*a2p*a7p + a10p*a4p*a7p + a11p*a2p*a8p - a10p*a3p*a8p);
            
            P_recoil_angle_rec_boost[0]=Precoil_angle_rec_boost*a3p;
            P_recoil_angle_rec_boost[1]=Precoil_angle_rec_boost*a7p;
            P_recoil_angle_rec_boost[2]=Precoil_angle_rec_boost*a11p;
            Erecoil_angle_rec_boost= sqrt( pow(Precoil_angle_boost,2) + pow(mn,2));
            theta_Precoil_angle_rec_boost = TMath::ACos(P_recoil_angle_rec_boost[2]/Precoil_angle_rec_boost) * TMath::RadToDeg();
            phi_Precoil_angle_rec_boost = TMath::ATan2(P_recoil_angle_rec_boost[1],P_recoil_angle_rec_boost[0]) * TMath::RadToDeg();
            if(TMath::Abs(phi_Precoil_boost-phi_Precoil_angle_rec_boost) > 30){
                if( phi_Precoil_angle_rec_boost > 0){
                    phi_Precoil_angle_rec_boost -= 180;
                }else{
                    phi_Precoil_angle_rec_boost += 180;
                }
            }

            Pcm_angle_rec_boost=
            (a11p*(a2p*a5p - a1p*a6p)*P1_rec_boost + a10p*(-a3p*a5p*P1_rec_boost + a1p*a7p*P1_rec_boost) + (a3p*a6p - a2p*a7p)*(a9p*P1_rec_boost -Ac12*Pbeam))/
            (-a12p*a3p*a6p + a11p*a4p*a6p + a12p*a2p*a7p - a10p*a4p*a7p - a11p*a2p*a8p + a10p*a3p*a8p);
            
            P_cm_angle_rec_boost[0]=Pcm_angle_rec_boost*a4p;
            P_cm_angle_rec_boost[1]=Pcm_angle_rec_boost*a8p;
            P_cm_angle_rec_boost[2]=Pcm_angle_rec_boost*a12p;
            Ecm_angle_rec_boost= sqrt( pow(Pcm_angle_rec_boost,2) + pow(mB10,2));
            theta_Pcm_angle_rec_boost = TMath::ACos(P_cm_angle_rec_boost[2]/Pcm_angle_rec_boost) * TMath::RadToDeg();
            phi_Pcm_angle_rec_boost = TMath::ATan2(P_cm_angle_rec_boost[1],P_cm_angle_rec_boost[0]) * TMath::RadToDeg();
            if(TMath::Abs(phi_Pcm_boost-phi_Pcm_angle_rec_boost) > 30){
                if( phi_Pcm_angle_rec_boost > 0){
                    phi_Pcm_angle_rec_boost -= 180;
                }else{
                    phi_Pcm_angle_rec_boost += 180;
                }
            }
            
            
            
            // quantities reconstructed from direct measurement
            TLorentzVector v1final_rec( TVector3(P_1_rec_boost[0],P_1_rec_boost[1],P_1_rec_boost[2]),              TMath::Sqrt( pow(P1_rec_boost,2) + pow(mp,2) ) );
            TLorentzVector v2final_rec( TVector3(P_2_rec_boost[0],P_2_rec_boost[1],P_2_rec_boost[2]),              TMath::Sqrt( pow(P2_rec_boost,2) + pow(mp,2) ) );
            TLorentzVector vrecoil_rec( TVector3(P_recoil_rec_boost[0],P_recoil_rec_boost[1],P_recoil_rec_boost[2]),              TMath::Sqrt( pow(Precoil_rec_boost,2) + pow(mn,2) ) );
            TLorentzVector vcm_rec( TVector3(P_cm_rec_boost[0],P_cm_rec_boost[1],P_cm_rec_boost[2]),                              TMath::Sqrt( pow(Pcm_rec_boost,2) + pow(mB10,2) ) );

        
            // quantities reconstructed from angle measurement and conservation equations
            TLorentzVector v2final_angle_rec( TVector3(P_2_angle_rec_boost[0],P_2_angle_rec_boost[1],P_2_angle_rec_boost[2]),TMath::Sqrt( pow(P2_angle_rec_boost,2) + pow(mp,2) ) );
            TLorentzVector vrecoil_angle_rec( TVector3(P_recoil_angle_rec_boost[0],P_recoil_angle_rec_boost[1],P_recoil_angle_rec_boost[2]),TMath::Sqrt( pow(Precoil_angle_rec_boost,2) + pow(mn,2) ) );
            TLorentzVector vcm_angle_rec( TVector3(P_cm_angle_rec_boost[0],P_cm_angle_rec_boost[1],P_cm_angle_rec_boost[2]),TMath::Sqrt( pow(Pcm_angle_rec_boost,2) + pow(mB10,2) ) );
            
            TLorentzVector vstanding_proton( TVector3(0,0,0), mp);
            
            
            //
            // Smeared t, u, s
            //
            t_rec_boost = (vstanding_proton-v1final_rec)*(vstanding_proton-v1final_rec);
            u_rec_boost = (vstanding_proton-v2final_rec)*(vstanding_proton-v2final_rec);
            s_rec_boost = (v1final_rec+v2final_rec)*(v1final_rec+v2final_rec);
            
            //
            // Smeared and angle reconstructed t, u, s
            //
            t_angle_rec_boost = (vstanding_proton-v1final_rec)*(vstanding_proton-v1final_rec);
            u_angle_rec_boost = (vstanding_proton-v2final_angle_rec)*(vstanding_proton-v2final_angle_rec);
            s_angle_rec_boost = (v1final_rec+v2final_angle_rec)*(v1final_rec+v2final_angle_rec);
            
            
            
            
            
            
            
            
            
            
            //
            //  Boost back from smeared DUBNA to smeared GSI
            //
            
            //boost from DUBNA to GSI
            v1final_rec.Boost(-vc12_ForReconstructions.BoostVector());
            v2final_rec.Boost(-vc12_ForReconstructions.BoostVector());
            vrecoil_rec.Boost(-vc12_ForReconstructions.BoostVector());
            vcm_rec.Boost    (-vc12_ForReconstructions.BoostVector());
            
            v2final_angle_rec.Boost(-vc12_ForReconstructions.BoostVector());
            vrecoil_angle_rec.Boost(-vc12_ForReconstructions.BoostVector());
            vcm_angle_rec.Boost    (-vc12_ForReconstructions.BoostVector());
            
            
            
            P_1_rec[0] = v1final_rec.X();
            P_1_rec[1] = v1final_rec.Y();
            P_1_rec[2] = v1final_rec.Z();
            P1_rec = v1final_rec.P();
            E1_rec=v1final_rec.T();
            theta_P1_rec = v1final_rec.Theta()*TMath::RadToDeg();
            phi_P1_rec = v1final_rec.Phi()*TMath::RadToDeg();
            
            P_2_rec[0] = v2final_rec.X();
            P_2_rec[1] = v2final_rec.Y();
            P_2_rec[2] = v2final_rec.Z();
            P2_rec = v2final_rec.P();
            E2_rec = v2final_rec.T();
            theta_P2_rec = v2final_rec.Theta()*TMath::RadToDeg();
            phi_P2_rec = v2final_rec.Phi()*TMath::RadToDeg();
            
            
            P_recoil_rec[0] = vrecoil_rec.X();
            P_recoil_rec[1] = vrecoil_rec.Y();
            P_recoil_rec[2] = vrecoil_rec.Z();
            Precoil_rec = vrecoil_rec.P();
            Erecoil_rec = vrecoil_rec.T();
            theta_Precoil_rec = vrecoil_rec.Theta()*TMath::RadToDeg();
            phi_Precoil_rec = vrecoil_rec.Phi()*TMath::RadToDeg();
            
            P_cm_rec[0] = vcm_rec.X();
            P_cm_rec[1] = vcm_rec.Y();
            P_cm_rec[2] = vcm_rec.Z();
            Pcm_rec = vcm_rec.P();
            Ecm_rec = vcm_rec.T();
            theta_Pcm_rec = vcm_rec.Theta()*TMath::RadToDeg();
            phi_Pcm_rec = vcm_rec.Phi()*TMath::RadToDeg();
            
            P_miss_rec[0] = P_1_rec[0] + P_2_rec[0];
            P_miss_rec[1] = P_1_rec[1] + P_2_rec[1];
            P_miss_rec[2] = P_1_rec[2] + P_2_rec[2] + v_proton_beam_UnSMEARED.P();
            Pmiss_rec = TMath::Sqrt(pow(P_miss_rec[0],2) + pow(P_miss_rec[1],2) + pow(P_miss_rec[2],2));
            Emiss_rec = TMath::Sqrt( pow(Pmiss_rec,2) + pow(mp,2) );
            theta_Pmiss_rec = TMath::ACos(P_miss_rec[2]/Pmiss_rec)*TMath::RadToDeg();      // Calculate theta [degrees]
            phi_Pmiss_rec = TMath::ATan2(P_miss_rec[1],P_miss_rec[0])*TMath::RadToDeg();       // Calculate phi   [degrees]

            
            s_rec= (v1final_rec+v2final_rec)*(v1final_rec+v1final_rec);
            u_rec= (v_proton_beam_UnSMEARED-v2final_rec)*(v_proton_beam_UnSMEARED-v2final_rec);
            t_rec= (v_proton_beam_UnSMEARED-v1final_rec)*(v_proton_beam_UnSMEARED-v1final_rec);
            
            
            P_2_angle_rec[0] = v2final_angle_rec.X();
            P_2_angle_rec[1] = v2final_angle_rec.Y();
            P_2_angle_rec[2] = v2final_angle_rec.Z();
            P2_angle_rec = v2final_angle_rec.P();
            E2_angle_rec = v2final_angle_rec.T();
            theta_P2_angle_rec = v2final_angle_rec.Theta()*TMath::RadToDeg();
            phi_P2_angle_rec = v2final_angle_rec.Phi()*TMath::RadToDeg();
            
            P_recoil_angle_rec[0] = vrecoil_angle_rec.X();
            P_recoil_angle_rec[1] = vrecoil_angle_rec.Y();
            P_recoil_angle_rec[2] = vrecoil_angle_rec.Z();
            Precoil_angle_rec = vrecoil_angle_rec.P();
            Erecoil_angle_rec = vrecoil_angle_rec.T();
            theta_Precoil_angle_rec = vrecoil_angle_rec.Theta()*TMath::RadToDeg();
            phi_Precoil_angle_rec = vrecoil_angle_rec.Phi()*TMath::RadToDeg();
            
            P_cm_angle_rec[0] = vcm_angle_rec.X();
            P_cm_angle_rec[1] = vcm_angle_rec.Y();
            P_cm_angle_rec[2] = vcm_angle_rec.Z();
            Pcm_angle_rec = vcm_angle_rec.P();
            Ecm_angle_rec = vcm_angle_rec.T();
            theta_Pcm_angle_rec = vcm_angle_rec.Theta()*TMath::RadToDeg();
            phi_Pcm_angle_rec = vcm_angle_rec.Phi()*TMath::RadToDeg();
            
            P_miss_angle_rec[0] = P_1_rec[0] + P_2_angle_rec[0];
            P_miss_angle_rec[1] = P_1_rec[1] + P_2_angle_rec[1];
            P_miss_angle_rec[2] = P_1_rec[2] + P_2_angle_rec[2]+v_proton_beam_UnSMEARED.P();
            Pmiss_angle_rec = TMath::Sqrt(pow(P_miss_angle_rec[0],2) + pow(P_miss_angle_rec[1],2) + pow(P_miss_angle_rec[2],2));
            Emiss_angle_rec = TMath::Sqrt( pow(Pmiss_angle_rec,2) + pow(mp,2) );
            theta_Pmiss_angle_rec = TMath::ACos(P_miss_angle_rec[2]/Pmiss_angle_rec)*TMath::RadToDeg();      // Calculate theta [degrees]
            phi_Pmiss_angle_rec = TMath::ATan2(P_miss_angle_rec[1],P_miss_angle_rec[0])*TMath::RadToDeg();       // Calculate phi   [degrees]
            
            
            t_angle_rec = (v_proton_beam_UnSMEARED-v1final_rec)*(v_proton_beam_UnSMEARED-v1final_rec);
            u_angle_rec = (v_proton_beam_UnSMEARED-v2final_angle_rec)*(v_proton_beam_UnSMEARED-v2final_angle_rec);
            s_angle_rec = (v1final_rec+v2final_angle_rec)*(v1final_rec+v2final_angle_rec);
            
            
            
            
            
            //
            // Calculate recoil positions for acceptance cuts
            //
            Recoil_x_NeuLand = (disttoZDC/TMath::Cos(theta_Precoil_boost*TMath::DegToRad())) * TMath::Sin(theta_Precoil_boost*TMath::DegToRad()) * TMath::Cos(phi_Precoil_boost*TMath::DegToRad());
            Recoil_y_NeuLand = (disttoZDC/TMath::Cos(theta_Precoil_boost*TMath::DegToRad())) * TMath::Sin(theta_Precoil_boost*TMath::DegToRad()) * TMath::Sin(phi_Precoil_boost*TMath::DegToRad());
            Recoil_z_MagnetHeight = (0.6/(TMath::Sin(theta_Precoil_boost*TMath::DegToRad())  * TMath::Sin(phi_Precoil_boost*TMath::DegToRad()))) * TMath::Cos(theta_Precoil_boost*TMath::DegToRad());
            
            //
            // Calculate P1 positions for acceptance cuts
            // x= r * cos(phi) * sin(theta)
            // y= r * sin(phi) * sin(theta)
            // z= r * cons(theta)
            // sqrt(x^2 + z^2) = d
            // sqrt(x^2 + y^2 + z^2) = r
            //
            r1 = dist1/TMath::Sqrt( 1-pow(TMath::Sin(phi_P1_boost*TMath::DegToRad())*TMath::Sin(theta_P1_boost*TMath::DegToRad()),2) );
            P1_x = r1 * TMath::Cos(phi_P1_boost*TMath::DegToRad()) * TMath::Sin(theta_P1_boost*TMath::DegToRad());
            P1_y = r1 * TMath::Sin(phi_P1_boost*TMath::DegToRad()) * TMath::Sin(theta_P1_boost*TMath::DegToRad());
            P1_z = r1 * TMath::Cos(theta_P1_boost*TMath::DegToRad());
            
            //
            // Calculate P2 positions for acceptance cuts
            // x= r * cos(phi) * sin(theta)
            // y= r * sin(phi) * sin(theta)
            // z= r * cons(theta)
            // sqrt(x^2 + z^2) = d
            // sqrt(x^2 + y^2 + z^2) = r
            //
            
            r2 = dist2/TMath::Sqrt( 1-pow(TMath::Sin(phi_P2_boost*TMath::DegToRad())*TMath::Sin(theta_P2_boost*TMath::DegToRad()),2) );
            P2_x = r2 * TMath::Cos(phi_P2_boost*TMath::DegToRad()) * TMath::Sin(theta_P2_boost*TMath::DegToRad());
            P2_y = r2 * TMath::Sin(phi_P2_boost*TMath::DegToRad()) * TMath::Sin(theta_P2_boost*TMath::DegToRad());
            P2_z = r2 * TMath::Cos(theta_P2_boost*TMath::DegToRad());
            
            
            
            if(abs(t_rec)>2 && abs(u_rec)>2 && abs(s_rec)>2){ // only consider 'hard' scattering events.
                T->Fill();
            } //end of if statement
            
        } //end of k
        
        
    } //end of i
    
    T->Write();
    f->Write();
    f->Close();
    
    
}






