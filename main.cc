#include "PimEpEm.h"
//#include "PimEpEp.h"
//#include "PimEmEm.h"
#include "data.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include "HFilter.h"


using namespace std;
using namespace PATData;

int main()
{
//******************************************* B E A M   P A R A M E T E R S *****************************/
    /*
    //proj = new TLorentzVector(0,0,1700,1705.71972227); // PION BEAM 1.7 GeV/c momentum
   //targ = new TLorentzVector(0,0,0,938.27231); // PROTON
   beam = new TLorentzVector(0,0,0,0);
   proj = new TLorentzVector(0,0,700,713.7785); // PION BEAM 0.7 GeV/c momentum 
   targ = new TLorentzVector(0,0,0,938.27231); // PROTON
   *beam = *proj + *targ;
   */

   //PION BEAM SCANNER
//   double pion_momentum = 612.0;

  
//     double pion_momentum = 708.0;      // for PWA HADES weights

  double pion_momentum = 701.;             //This is HADES-2020 new momentum for analysis
  //double pion_momentum = 690.;             //This is HADES-2020 new momentum for analysis
  double pion_momentum_B = 690.;             //This is HADES-2020 new momentum for analysis
  
  
  //   double pion_momentum = 687.2;      // for PLUTO generated events
//  double pion_momentum = 691.0;      // for PWA analysis events
//   double pion_momentum = 748.0;
//   double pion_momentum = 800.0;

    //d_pion_mom = TMath::Abs( 801.167 - 788.234 ); // 800

   // ORIGIN: https://hades-wiki.gsi.de/foswiki/bin/view/PionBeam/WebHome

   double pion_energy = sqrt( pion_momentum*pion_momentum + 139.56995*139.56995 );
   double pion_energy_B = sqrt( pion_momentum_B*pion_momentum_B + 139.56995*139.56995 );
   
   proj = new TLorentzVector(0,0, pion_momentum, pion_energy); // PION BEAM momentum as above
   proj_B = new TLorentzVector(0,0, pion_momentum_B, pion_energy_B); // PION BEAM momentum for SIM as above
   
   targ = new TLorentzVector(0,0,0,938.27231); // PROTON
/*******************************************************************************************************/
   beam = new TLorentzVector(0,0,0,0);
   beam_PT = new TLorentzVector(0,0,0,0);
   beam_var = new TLorentzVector(0,0,0,0);
   proj_var = new TLorentzVector(0,0,0,0);
   proj_PT = 0;

   *beam = *proj + *targ;

//*******************************************************************************************************/
/*******************************************************************************************************/
   beam_B = new TLorentzVector(0,0,0,0);
   beam_PT_B = new TLorentzVector(0,0,0,0);
   beam_var_B = new TLorentzVector(0,0,0,0);
   proj_var_B = new TLorentzVector(0,0,0,0);
   proj_PT_B = 0;

   *beam_B = *proj_B + *targ;
   
//*******************************************************************************************************/


   filter = new HFilter;

// *********************** FILES WITH CUTS ****************************
   insideEmS0 = -1;
   insideEmS1 = -1;
   insideEpS0 = -1;
   insideEpS1 = -1;

   insideEm1S0 = -1;
   insideEm1S1 = -1;
   insideEm2S0 = -1;
   insideEm2S1 = -1;

   insideEp1S0 = -1;
   insideEp1S1 = -1;
   insideEp2S0 = -1;
   insideEp2S1 = -1;

   pim = new TLorentzVector(0,0,0,0);
   e1 = new TLorentzVector(0,0,0,0);
   e2 = new TLorentzVector(0,0,0,0);
   gammae1e2 = new TLorentzVector(0,0,0,0);
   e1e2 = new TLorentzVector(0,0,0,0);
   e1_delta = new TLorentzVector(0,0,0,0);
   e2_delta = new TLorentzVector(0,0,0,0);
   e1e2_miss = new TLorentzVector(0,0,0,0);

   pimepem = new TLorentzVector(0,0,0,0);
   pimepem_miss = new TLorentzVector(0,0,0,0);
   ppimepem_miss = new TLorentzVector(0,0,0,0);
   pimepem_miss_PT = new TLorentzVector(0,0,0,0);

 /*******************************************************************************************************/
   pim_B = new TLorentzVector(0,0,0,0);
   e1_B = new TLorentzVector(0,0,0,0);
   e2_B = new TLorentzVector(0,0,0,0);
   gammae1e2_B = new TLorentzVector(0,0,0,0);
   e1e2_B = new TLorentzVector(0,0,0,0);
   e1e2_miss_B = new TLorentzVector(0,0,0,0);
   e1_delta_B = new TLorentzVector(0,0,0,0);
   e2_delta_B = new TLorentzVector(0,0,0,0);
   pimepem_B = new TLorentzVector(0,0,0,0);
   pimepem_miss_B = new TLorentzVector(0,0,0,0);
   pimepem_miss_PT_B = new TLorentzVector(0,0,0,0);
   ppimepem_miss_B = new TLorentzVector(0,0,0,0);


   /*******************************************************************************************************/
   
   
   /*******************************************************************************************************/

   //   outFileData = new TFile("PIMEPEM_PWA_OA9deg_R4_BTMax2_HADESwtgs_P701_Mom_genwtgs_CloseCut0.root","recreate");        //Important Hades weights 
   //outFileData = new TFile("PPIMEPEM_PWA_OA4deg_R4_BT2_MANDLEY_HADESwtgs_P687_2_Mom_NewCut.root","recreate");        //Important Manley-Hades weights   

   //outFileData = new TFile("PIMEPEM_BREM_ANGULAR_OA9deg_P687_2_Mom_CloseCut_FIT4.root","recreate");        //Important Brem Weights
   //outFileData = new TFile("PIMEPEM_BREM_PHASESPACE_OA4deg_P687_2_Mom_CloseCut_FIT4.root","recreate");        //Important Brem Weights

   //outFileData = new TFile("PIMEPEM_PWA_DELTA_OA4deg_HADESwtgs_P701_Mom_genwtgs_CloseCut_FIT4_momTEST.root","recreate");        // PWA PI0-Mom_Correction TEST 

   //---------------------------------------------------------
//---------------------------------------------------------
   //outFileData = new TFile("PIMEPEM_BREMSSTRAHLUNG_OA4deg_HADESwtgs_P690_CloseCut_FIT4_PHASESPACE.root","recreate");        // Bremsstrahlung_FINAL 
   //outFileData = new TFile("PIMEPEM_BREMSSTRAHLUNG_OA4deg_HADESwtgs_P690_CloseCut_FIT4_ANGULAR_MM100MeV_LEVEL3.root","recreate");        // Bremsstrahlung_FINAL 

/////////////////////////////////////////////////////////////////////////////////   ///////////////////////////////////////////////////////////////////////////////// /////////////////////////////////////////////////////////////////////////////////   ////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////   /////////////////////////////////////////////////////////////////////////////////

//outFileData = new TFile("PIMEPEM_BREMSSTRAHLUNG_OA9deg_HADESwtgs_P690_CloseCut_FIT4_ANGULAR_MM100MeV_DEC2020.root","recreate");
// Bremsstrahlung_FINAL 
//outFileData = new TFile("PIMEPEM_CARBON_BREMSSTRAHLUNG_OA4deg_HADESwtgs_P690_CloseCut_FIT4_ANGULAR_MM100MeV_DEC2020.root","recreate");
outFileData = new TFile("PIMEPEM_CARBON_BREM_OA4deg_P690_CloseCut_FIT0_ANGULAR_MM100MeV_DEC2020_EFF_TESTING.root","recreate");
//outFileData = new TFile("PIMEPEM_PROTON_BREM_OA4deg_P690_CloseCut_FIT0_ANGULAR_MM100MeV_DEC2020_EFF_TESTING.root","recreate");
// Bremsstrahlung_FINAL 

//outFileData = new TFile("PIMEPEM_PE_OA4_CC4_TEST.root","recreate");        // Brems
/////////////////////////////////////////////////////////////////////////////////   ///////////////////////////////////////////////////////////////////////////////// /////////////////////////////////////////////////////////////////////////////////   ////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////   /////////////////////////////////////////////////////////////////////////////////

   //outFileData = new TFile("PIMEPEM_CARBON_BREMSSTRAHLUNG_OA4deg_HADESwtgs_P690_CloseCut_FIT4_PHASESPACE_3_3_MM100MeV_LEVEL3.root","recreate");        // Bremsstrahlung_FINAL 
   //outFileData = new TFile("PIMEPEM_CARBON_BREMSSTRAHLUNG_OA9deg_HADESwtgs_P690_CloseCut_FIT4_ANGULAR.root","recreate");        // Bremsstrahlung_FINAL 
//---------------------------------------------------------
//---------------------------------------------------------
   
   //outFileData = new TFile("PIMEPEM_PWA_PI0_OA4deg_HADESwtgs_P701_Mom_genwtgs_CloseCut_FIT4.root","recreate");        // PWA PI0-DALITZ HADES WEIGHTS 
   //outFileData = new TFile("PIMEPEM_PWA_DeltaDalitz_OA9deg_HADESwtgs_P701_Mom_genwtgs_CloseCut_FIT4.root","recreate");        // PWA DELTA-DALITZ HADES WEIGHTS 

   //outFileData = new TFile("PIMEPEM_PWA_DeltaDalitz_OA4deg_HADESwtgs_P701_Mom_genwtgs_CloseCut0.root","recreate");        // PWA DELTA-DALITZ HADES WEIGHTS 
   //outFileData = new TFile("PPIMEPEM_PWA_OA4deg_R4_BT2_MANDLEY_HADESwtgs_P687_2_Mom_NewCut.root","recreate");        //PWA Delta-dalitz Manley-Hades weights   

   //outFileData = new TFile("PIMEPEM_PWA_Bremss_phasespace_OA4deg_HADESwtgs_P687_2_Mom_NewCut_genwtgs.root","recreate");        // PWA Delta-dalitz Hades weights 
   
/************************************** O U T P U T   F I L E ******************************************/
   
   //outFileData = new TFile("PIMEPEM_BREMSSTRAHLUNG_4deg_MAX.root","recreate");
   // outFileData = new TFile("PIMEPEM_BREMSSTRAHLUNG_PHASESPACE_4deg.root","recreate");
   //outFileData = new TFile("PIMEPEM_BREMSSTRAHLUNG_4deg_OLD_MAX.root","recreate");
   //outFileData = new TFile("SIM_Delta1_p_687Mom_OA9deg_8rings.root","recreate");

   //outFileData = new TFile("BREM_Efficiency_PLUTO_SIM_687Mom_9deg_4rings_bt2rings.root","recreate");

   //outFileData = new TFile("PIMEPEM_BREM_9deg_RF4rings_bt2rings_MAX_phasespace22_parametrized.root","recreate");
   //outFileData = new TFile("PIMEPEM_BREM_4deg_RF4rings_bt2rings_MAX_angular22_parametrized.root","recreate");
   //outFileData = new TFile("PIMEPEM_BREM_4deg_RF4rings_bt2rings_MAX_phasespace22.root","recreate");
   //outFileData = new TFile("PIMEPEM_BREM_9deg_RF4rings_bt2rings_MAX_angular22.root","recreate");



   //outFileData = new TFile("PIMEPEM_PWA_1M_691_4deg_RF4rings_bt2rings_MAX.root","recreate");

   //outFileData = new TFile("PIMEPEM_BREMSSTRAHLUNG_9deg_MAX_Neutron_phasespace22.root","recreate");
   
   //outFileData = new TFile("PWA_momentums-ppimpi0_ev_ALL_OA9deg.root","recreate");
   //outFileData = new TFile("PWA_momentums-ppimpi0_2ndSol_dcs_ALL_OA9deg.root","recreate");


/*******************************************************************************************************/

/************************** control ntuple ***************************/
   tlo = new HNtuple("nt","nt");
   tlo->setFile( outFileData );
/*********************************************************************/
   mom_miss = new TH2F("mom_miss","mom_miss",400,680,700,200,7000, 30000); mom_miss->Sumw2();
   
   epem_Mass = new TH1F("epem_Mass","epem_Mass",111,0.0,0.555); epem_Mass->Sumw2();
   epem_Mass_Brem = new TH1F("epem_Mass_Brem","epem_Mass_Brem",111,0.0,0.555); epem_Mass_Brem->Sumw2();
      
   Missing_mass = new TH1F("Missing_mass","Missing_mass",300,0,2);
   Missing_mass->Sumw2();
   Full_signal = new TH1F("Full_signal","Full_signal",300,0,2);
   Full_signal->Sumw2();
   Missing_mass_back1 = new TH1F("Missing_mass_back1","Missing_mass_back1",300,0,2);
   Missing_mass_back1->Sumw2();
   Full_signal_back1 = new TH1F("Full_signal_back1","Full_signal_back1",300,0,2);
   Full_signal_back1->Sumw2();
   Missing_mass_back2 = new TH1F("Missing_mass_back2","Missing_mass_back2",300,0,2);
   Missing_mass_back2->Sumw2();
   Full_signal_back2 = new TH1F("Full_signal_back2","Full_signal_back2",300,0,2);
   Full_signal_back2->Sumw2();

   Miss_Mass = new TH1F("Miss_Mass","Miss_Mass",300,0,2);
   Miss_Mass->Sumw2();
   
   sig_all = new TH1F("sig_all","sig_all",300,0,1);
   sig_all->Sumw2();
   sig_all_back1 = new TH1F("sig_all_back1","sig_all_back1",300,0,1);
   sig_all_back1->Sumw2();
   sig_all_back2 = new TH1F("sig_all_back2","sig_all_back2",300,0,1);
   sig_all_back2->Sumw2();
   sig_OK = 0;

   miss_all = new TH1F("miss_all","miss_all",80,0.6,1.4);
   miss_all->Sumw2();
   miss_all_back1 = new TH1F("miss_all_back1","miss_all_back1",80,0.6,1.4);
   miss_all_back1->Sumw2();
   miss_all_back2 = new TH1F("miss_all_back2","miss_all_back2",80,0.6,1.4);
   miss_all_back2->Sumw2();
   miss_OK = 0;


   EpEm_inv_mass = new TH1F("EpEm_inv_mass","EpEm_inv_mass",300,0,1); EpEm_inv_mass->Sumw2();
   EpEm_inv_mass1 = new TH1F("EpEm_inv_mass_1MeV","EpEm_inv_mass_1MeV",1000,0,1); EpEm_inv_mass1->Sumw2();
   EpEm_inv_mass11 = new TH1F("EpEm_inv_mass_1MeV_more800MeV","EpEm_inv_mass_1MeV_more800MeV",1000,0,1); EpEm_inv_mass11->Sumw2();

   EpEm_inv_mass1_M2 = new TH1F("EpEm_inv_mass_1MeV_M2","EpEm_inv_mass_1MeV_M2",1000,0,1); EpEm_inv_mass1_M2->Sumw2();

   EpEm_inv_mass11_M2 = new TH1F("EpEm_inv_mass_1MeV_more800MeV_M2","EpEm_inv_mass_1MeV_more800MeV_M2",1000,0,1); EpEm_inv_mass11_M2->Sumw2();
   EpEm_inv_mass1_ph = new TH1F("EpEm_inv_mass_1MeV_ph","EpEm_inv_mass_1MeV_ph",1000,0,1); EpEm_inv_mass1_ph->Sumw2();
   EpEm_inv_mass11_ph = new TH1F("EpEm_inv_mass_1MeV_more800MeV_ph","EpEm_inv_mass_1MeV_more800MeV_ph",1000,0,1); EpEm_inv_mass11_ph->Sumw2();
   EpEm_inv_mass1_Brem = new TH1F("EpEm_inv_mass_1MeV_Brem","EpEm_inv_mass_1MeV_Brem",1000,0,1); EpEm_inv_mass1_Brem->Sumw2();
   EpEm_inv_mass11_Brem = new TH1F("EpEm_inv_mass_1MeV_more1400MeV_Brem","EpEm_inv_mass_1MeV_more1400MeV_Brem",1000,0,1); EpEm_inv_mass11_Brem->Sumw2();
   EpEm_inv_mass_Brem_reg = new TH1F("EpEm_inv_mass_region_more1400MeV_Brem","EpEm_inv_mass_region_more1400MeV_Brem",1000,0,1); EpEm_inv_mass_Brem_reg->Sumw2();
   EpEm_inv_mass_Brem_nopi0cut = new TH1F("EpEm_inv_mass_Brem_1MeV_nopi0cut","EpEm_inv_mass_Brem_1MeV_nopi0cut",1000,0,1);
   EpEm_inv_mass_Brem_nopi0cut->Sumw2();
   
   sig_all_MeV = new TH1F("EpEm_5MeV","EpEm_5MeV",200,0,1000); sig_all_MeV->Sumw2();
   
   Missing_mass11 = new TH1F("Three_particle_missing_mass_2MeV","Three_particle_missing_mass_2MeV",1000,0,2); Missing_mass11->Sumw2();
   Missing_mass111 = new TH1F("Three_particle_missing_mass_2MeV_more800MeV","Three_particle_missing_mass_2MeV_more800MeV",1000,0,2); Missing_mass111->Sumw2();
   Missing_mass11_M2 = new TH1F("Three_particle_missing_mass_2MeV_M2","Three_particle_missing_mass_2MeV_M2",1000,0,2); Missing_mass11_M2->Sumw2();
   Missing_mass111_M2 = new TH1F("Three_particle_missing_mass_2MeV_more1400MeV_M2","Three_particle_missing_mass_2MeV_more1400MeV_M2",1000,0,2); Missing_mass111_M2->Sumw2();
   Missing_mass11_ph = new TH1F("Three_particle_missing_mass_2MeV_ph","Three_particle_missing_mass_2MeV_ph",1000,0,2); Missing_mass11_ph->Sumw2();
   Missing_mass111_ph = new TH1F("Three_particle_missing_mass_2MeV_more1400MeV_ph","Three_particle_missing_mass_2MeV_more1400MeV_ph",1000,0,2); Missing_mass111_ph->Sumw2();
   Missing_mass11_Brem = new TH1F("Three_particle_missing_mass_2MeV_Brem","Three_particle_missing_mass_2MeV_Brem",1000,0,2); Missing_mass11_Brem->Sumw2();
   Missing_mass111_Brem = new TH1F("Three_particle_missing_mass_2MeV_more1400MeV_Brem","Three_particle_missing_mass_2MeV_more1400MeV_Brem",1000,0,2); Missing_mass111_Brem->Sumw2();
   
   Missing_mass_Brem_nopi0cut = new TH1F("Missing_mass_Brem_nopi0cut","Missing_mass_Brem_nopi0cut",1000,0,2);
   Missing_mass_Brem_nopi0cut->Sumw2();  
   
   //---------- Angular Distribution -----------------

   
   Miss_ang_dis_CosTheta = new TH1F("Three_particle_missmass_ang_dis_CosTheta","Three_particle_missmass_ang_dis_CosTheta",200,-1,1); Miss_ang_dis_CosTheta->Sumw2();
   pim_ang_dis_CosTheta = new TH1F("Pim_angular_distribution_CosTheta","Pim_angular_distribution_CosTheta",200,-1,1); pim_ang_dis_CosTheta->Sumw2();
   ep_ang_dis_CosTheta = new TH1F("ep_angular_distribution_CosTheta","ep_angular_distribution_CosTheta",200,-1,1); ep_ang_dis_CosTheta->Sumw2();
   em_ang_dis_CosTheta = new TH1F("em_angular_distribution_CosTheta","em_angular_distribution_CosTheta",200,-1,1); em_ang_dis_CosTheta->Sumw2();
   Gamma_ang_dis_CosTheta = new TH1F("Gamma_angular_distribution_CosTheta","Gamma_angular_distribution_CosTheta",200,-1,1); Gamma_ang_dis_CosTheta->Sumw2();


   //---------- Angular Distribution ONLY THETA -----------------

   Miss_ang_dis_Theta = new TH1F("Three_particle_missmass_ang_dis_Theta","Three_particle_missmass_ang_dis_Theta",90,0,90); Miss_ang_dis_Theta->Sumw2();
   pim_ang_dis_Theta = new TH1F("Pim_angular_distribution_Theta","Pim_angular_distribution_Theta",90,0,90); pim_ang_dis_Theta->Sumw2();
   ep_ang_dis_Theta = new TH1F("ep_angular_distribution_Theta","ep_angular_distribution_Theta",90,0,90); ep_ang_dis_Theta->Sumw2();
   em_ang_dis_Theta = new TH1F("em_angular_distribution_Theta","em_angular_distribution_Theta",90,0,90); em_ang_dis_Theta->Sumw2();
   Gamma_ang_dis_Theta = new TH1F("Gamma_angular_distribution_Theta","Gamma_angular_distribution_Theta",90,0,90); Gamma_ang_dis_Theta->Sumw2();

   //----------- Bremsstrahlung reaction ------------
   Missing_mass1 = new TH1F("Three_particle_missing_mass","Three_particle_missing_mass",300,0,2); Missing_mass1->Sumw2();
   Missing_mass2 = new TH1F("Three_particle_missing_mass_2","Three_particle_missing_mass_2",300,0,2); Missing_mass2->Sumw2();
   Missing_momentum = new TH1F("Three_particle_missing_momentum","Three_particle_missing_momentum",150,-1,1500); Missing_momentum->Sumw2();
   Missing_energy = new TH1F("Three_particle_missing_energy","Three_particle_missing_energy",500,-1,1999); Missing_energy->Sumw2();

   //----------- Efficiency histograms ------------------
   EpEm_Pim_P = new TH2F("EpEm_Pim_P","EpEm_Pim_P",20,0,1000,20,0,1000); EpEm_Pim_P->Sumw2();
   EpEm_Pim_Theta = new TH2F("EpEm_Pim_Theta","EpEm_Pim_Theta",36,0,180,36,0,180); EpEm_Pim_Theta->Sumw2();

   pim_P_dis = new TH1F("Pim_P","Pim_P",20,0,1000); pim_P_dis->Sumw2();
   Gamma_P_dis = new TH1F("Gamma_P","Gamma_P",20,0,1000); Gamma_P_dis->Sumw2();

   //----------------------- Momentum Distribution ----------------------------------------------------

   Missing_TOTP = new TH1F("Brems_miss_mom","Brems_miss_mom",100,-500,500); Missing_TOTP->Sumw2();
   Missing_Px = new TH1F("Brems_miss_momX","Brems_miss_momX",100,-500,500); Missing_Px->Sumw2();
   Missing_Py = new TH1F("Brems_miss_momY","Brems_miss_momY",100,-500,500); Missing_Py->Sumw2();
   Missing_Pz = new TH1F("Brems_miss_momZ","Brems_miss_momZ",100,-500,500); Missing_Pz->Sumw2();
Float_t xxbinsang[] =  {0.00,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.93,0.95,0.97,1.00,1.05,1.10,1.15,1.20,1.25,1.30,1.35,1.40,1.45,1.50};
   Int_t nnbinsang = sizeof(xxbinsang)/sizeof(Float_t);

   miss_mass_var = new TH1F("miss_mass_var","miss_mass_var",nnbinsang-1,xxbinsang);
   miss_mass_var->Sumw2(); 

   Float_t xbins[] = {0.000,0.005,0.010,0.015,0.020,0.025,0.030,0.035,0.040,0.045,0.050,0.055,0.070,0.085,0.100,0.120,0.140,0.180,0.230,0.280,0.330,0.400,0.500,0.600};
   Int_t nbins = sizeof(xbins)/sizeof(Float_t);

   sig_all_var = new TH1F("sig_all_var","sig_all_var",nbins-1,xbins);
   sig_all_var->Sumw2();
   sig_all_var_Brem = new TH1F("sig_all_var_Brem","sig_all_var_Brem",nbins-1,xbins);
   sig_all_var_Brem->Sumw2();
   sig_all_var_Brem_nopi0cut = new TH1F("sig_all_var_Brem_nopi0cut","sig_all_var_Brem_nopi0cut",nbins-1,xbins);
   sig_all_var_Brem_nopi0cut->Sumw2();
   
   sig_all_var_back1 = new TH1F("sig_all_var_back1","sig_all_var_back1",nbins-1,xbins);
   sig_all_var_back1->Sumw2();
   sig_all_var_back2 = new TH1F("sig_all_var_back2","sig_all_var_back2",nbins-1,xbins);
   sig_all_var_back2->Sumw2();
   sig_var_OK = 0;

   sig_all_var2 = new TH1F("sig_all_var_no","sig_all_var_no_norm",nbins-1,xbins);
   sig_all_var2->Sumw2();
   sig_all_var2_back1 = new TH1F("sig_all_var_back1_no","sig_all_var_back1_no_norm",nbins-1,xbins);
   sig_all_var2_back1->Sumw2();
   sig_all_var2_back2 = new TH1F("sig_all_var_back2_no","sig_all_var_back2_no_norm",nbins-1,xbins);
   sig_all_var2_back2->Sumw2();
   sig_var2_OK = 0;

   Float_t xbinsang[] =  {-1.0, -0.9, -0.8, -0.7, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.7, 0.8, 0.9, 1.0};
   Int_t nbinsang = sizeof(xbinsang)/sizeof(Float_t);

   cos_ep = new TH1F("cos_ep","cos_ep",10,-1,1); cos_ep->Sumw2();
   cos_em = new TH1F("cos_em","cos_em",10,-1,1); cos_em->Sumw2();
   cos_sum = new TH1F("cos_sum","cos_sum",10,-1,1); cos_sum->Sumw2();
   cos_ep_back1 = new TH1F("cos_ep_back1","cos_ep_back1",10,-1,1); cos_ep_back1->Sumw2();
   cos_em_back1 = new TH1F("cos_em_back1","cos_em_back1",10,-1,1); cos_em_back1->Sumw2();
   cos_sum_back1 = new TH1F("cos_sum_back1","cos_sum_back1",10,-1,1); cos_sum_back1->Sumw2();
   cos_ep_back2 = new TH1F("cos_ep_back2","cos_ep_back2",10,-1,1); cos_ep_back2->Sumw2();
   cos_em_back2 = new TH1F("cos_em_back2","cos_em_back2",10,-1,1); cos_em_back2->Sumw2();
   cos_sum_back2 = new TH1F("cos_sum_back2","cos_sum_back2",10,-1,1); cos_sum_back2->Sumw2();

   cos_ep_cm = new TH1F("cos_ep_cm","cos_ep_cm",10,-1,1); cos_ep_cm->Sumw2();
   cos_back1_cm = new TH1F("cos_back1_cm","cos_back1_cm",10,-1,1); cos_back1_cm->Sumw2();
   cos_back2_cm = new TH1F("cos_back2_cm","cos_back2_cm",10,-1,1); cos_back2_cm->Sumw2();
   cos_ep_cm_OK = new TH1F("cos_ep_cm_OK","cos_ep_cm_OK",10,-1,1); cos_ep_cm_OK->Sumw2();


   rapidity_all = new TH1F("rapidity_all","rapidity_all",50,0,2); rapidity_all->Sumw2();
   rapidity_back1 = new TH1F("rapidity_back1","rapidity_back1",50,0,2); rapidity_back1->Sumw2();
   rapidity_back2 = new TH1F("rapidity_back2","rapidity_back2",50,0,2); rapidity_back2->Sumw2();
   rapidity_OK = 0;

   rapidity_140_all = new TH1F("rapidity_140_all","rapidity_140_all",50,0,2); rapidity_140_all->Sumw2();
   rapidity_140_back1 = new TH1F("rapidity_140_back1","rapidity_140_back1",50,0,2); rapidity_140_back1->Sumw2();
   rapidity_140_back2 = new TH1F("rapidity_140_back2","rapidity_140_back2",50,0,2); rapidity_140_back2->Sumw2();
   rapidity_140_OK = 0;

   pt_all = new TH1F("pt_all","pt_all",50,0,1); pt_all->Sumw2();
   pt_back1 = new TH1F("pt_back1","pt_back1",50,0,1); pt_back1->Sumw2();
   pt_back2 = new TH1F("pt_back2","pt_back2",50,0,1); pt_back2->Sumw2();
   pt_OK = 0;

   pt_140_all = new TH1F("pt_140_all","pt_140_all",50,0,1); pt_140_all->Sumw2();
   pt_140_back1 = new TH1F("pt_140_back1","pt_140_back1",50,0,1); pt_140_back1->Sumw2();
   pt_140_back2 = new TH1F("pt_140_back2","pt_140_back2",50,0,1); pt_140_back2->Sumw2();
   pt_140_OK = 0;



   Four_particle_MM = new TH1F("4_par_Missing_mom","4_par_Missing_mom",300,0,2);Four_particle_MM->Sumw2();




/****************************************************************************************/
/******************* T E S T I N G   P L O T S   S T A R T S ****************************/
/****************************************************************************************/

   Dilepton_OA_0deg = new TH1F("Dilepton_OA_0deg","Dilepton_OA_0deg",201,0,200); Dilepton_OA_0deg->Sumw2();
   epem_Mass_Brem_0deg = new TH1F("epem_Mass_Brem_0deg","epem_Mass_Brem_0deg",111,0.0,0.555); epem_Mass_Brem_0deg->Sumw2();

   Dilepton_OA_1deg = new TH1F("Dilepton_OA_1deg","Dilepton_OA_1deg",201,0,200); Dilepton_OA_1deg->Sumw2();
   epem_Mass_Brem_1deg = new TH1F("epem_Mass_Brem_1deg","epem_Mass_Brem_1deg",111,0.0,0.555); epem_Mass_Brem_1deg->Sumw2();

   Dilepton_OA_2deg = new TH1F("Dilepton_OA_2deg","Dilepton_OA_2deg",201,0,200); Dilepton_OA_2deg->Sumw2();
   epem_Mass_Brem_2deg = new TH1F("epem_Mass_Brem_2deg","epem_Mass_Brem_2deg",111,0.0,0.555); epem_Mass_Brem_2deg->Sumw2();

   Dilepton_OA_3deg = new TH1F("Dilepton_OA_3deg","Dilepton_OA_3deg",201,0,200); Dilepton_OA_3deg->Sumw2();
   epem_Mass_Brem_3deg = new TH1F("epem_Mass_Brem_3deg","epem_Mass_Brem_3deg",111,0.0,0.555); epem_Mass_Brem_3deg->Sumw2();

   Dilepton_OA_4deg = new TH1F("Dilepton_OA_4deg","Dilepton_OA_4deg",201,0,200); Dilepton_OA_4deg->Sumw2();
   epem_Mass_Brem_4deg = new TH1F("epem_Mass_Brem_4deg","epem_Mass_Brem_4deg",111,0.0,0.555); epem_Mass_Brem_4deg->Sumw2();

   Dilepton_OA_5deg = new TH1F("Dilepton_OA_5deg","Dilepton_OA_5deg",201,0,200); Dilepton_OA_5deg->Sumw2();
   epem_Mass_Brem_5deg = new TH1F("epem_Mass_Brem_5deg","epem_Mass_Brem_5deg",111,0.0,0.555); epem_Mass_Brem_5deg->Sumw2();

   Dilepton_OA_6deg = new TH1F("Dilepton_OA_6deg","Dilepton_OA_6deg",201,0,200); Dilepton_OA_6deg->Sumw2();
   epem_Mass_Brem_6deg = new TH1F("epem_Mass_Brem_6deg","epem_Mass_Brem_6deg",111,0.0,0.555); epem_Mass_Brem_6deg->Sumw2();

   Dilepton_OA_7deg = new TH1F("Dilepton_OA_7deg","Dilepton_OA_7deg",201,0,200); Dilepton_OA_7deg->Sumw2();
   epem_Mass_Brem_7deg = new TH1F("epem_Mass_Brem_7deg","epem_Mass_Brem_7deg",111,0.0,0.555); epem_Mass_Brem_7deg->Sumw2();

   Dilepton_OA_8deg = new TH1F("Dilepton_OA_8deg","Dilepton_OA_8deg",201,0,200); Dilepton_OA_8deg->Sumw2();
   epem_Mass_Brem_8deg = new TH1F("epem_Mass_Brem_8deg","epem_Mass_Brem_8deg",111,0.0,0.555); epem_Mass_Brem_8deg->Sumw2();

   Dilepton_OA_9deg = new TH1F("Dilepton_OA_9deg","Dilepton_OA_9deg",201,0,200); Dilepton_OA_9deg->Sumw2();
   epem_Mass_Brem_9deg = new TH1F("epem_Mass_Brem_9deg","epem_Mass_Brem_9deg",111,0.0,0.555); epem_Mass_Brem_9deg->Sumw2();

   Dilepton_OA_10deg = new TH1F("Dilepton_OA_10deg","Dilepton_OA_10deg",201,0,200); Dilepton_OA_10deg->Sumw2();
   epem_Mass_Brem_10deg = new TH1F("epem_Mass_Brem_10deg","epem_Mass_Brem_10deg",111,0.0,0.555); epem_Mass_Brem_10deg->Sumw2();



   Ep_OA = new TH1F("Ep_OA","Ep_OA",201,-50,50); Ep_OA->Sumw2();
   Em_OA = new TH1F("Em_OA","Em_OA",201,-50,50); Em_OA->Sumw2();
   Ep_Mom = new TH1F("Ep_Mom","Ep_Mom",200,0,2000); Ep_Mom->Sumw2();
   Em_Mom = new TH1F("Em_Mom","Em_Mom",200,0,2000); Em_Mom->Sumw2();


/****************************************************************************************/
/********************** T E S T I N G   P A R T   E N D S *******************************/
/****************************************************************************************/

   
/**************************** M A I N   P A R T ****************************************/

   PimEpEm t;
   //PimEpEp t_back1;
   //PimEmEm t_back2;
   cout << "START PIMEPEM!" << endl;
   t.Loop();
   /*cout << "START PIMEPEP!" << endl;
   t_back1.Loop();
   cout << "START PIMEMEM!" << endl;
   t_back2.Loop();*/
   cout << "FINALIZING!" << endl;
   

/***************************** F I N A L     C A L C U L A T I O N S ********************/

   sig_OK = (TH1F*)signal("sig_OK", sig_all, sig_all_back1, sig_all_back2);
   //normalize( sig_all );
   //normalize( sig_all_back1 );
   //normalize( sig_all_back2 );
   //normalize( sig_OK );

   miss_OK = (TH1F*)signal("miss_OK", miss_all, miss_all_back1, miss_all_back2);
   //normalize( miss_all );
   //normalize( miss_all_back1 );
   //normalize( miss_all_back2 );
   //normalize( miss_OK );

   sig_var_OK = (TH1F*)signal("sig_var_OK", sig_all_var, sig_all_var_back1, sig_all_var_back2);
   sig_var2_OK = (TH1F*)signal("sig_var2_OK", sig_all_var2, sig_all_var2_back1, sig_all_var2_back2);
   normalize( sig_all_var );
   normalize( sig_all_var_back1 );
   normalize( sig_all_var_back2 );
   normalize( sig_var_OK );

   SIGNAL = (TH1F*)sig_var2_OK->Clone("SIGNAL");
   SIGNAL->SetNameTitle("SIGNAL","SIGNAL");
   CB = (TH1F*)sig_all_var2_back1->Clone("CB");
   CB->Add(sig_all_var2_back2);
   CB->SetNameTitle("CB","CB");

        SIGNIFICANCE = (TH1F*)CB->Clone("SIGNIFICANCE");
        SIGNIFICANCE->Divide( SIGNAL );
        SIGNIFICANCE->Scale(2.);
        for (Int_t j=1; j<SIGNIFICANCE->GetNbinsX()+1; ++j)
        {
           SIGNIFICANCE->SetBinContent( j, SIGNIFICANCE->GetBinContent(j) + 1. );
        }
        TH1F *temp = (TH1F*)SIGNAL->Clone("temp");
        temp->Divide( SIGNIFICANCE );
        SIGNIFICANCE = temp;
        SIGNIFICANCE->SetNameTitle("SIGNIFICANCE", "SIGNIFICANCE");

   normalize( SIGNIFICANCE );
   normalize( SIGNAL );
   normalize( CB );

   cos_ep_OK = (TH1F*)signal("cos_ep_OK", cos_ep, cos_ep_back1, cos_ep_back2);
//   normalize( cos_ep );
//   normalize( cos_ep_back1 );
//   normalize( cos_ep_back2 );
//   normalize( cos_ep_OK );

   cos_em_OK = (TH1F*)signal("cos_em_OK", cos_em, cos_em_back1, cos_em_back2);
//   normalize( cos_em );
//   normalize( cos_em_back1 );
//   normalize( cos_em_back2 );
//   normalize( cos_em_OK );

   cos_sum_OK = (TH1F*)signal("cos_sum_OK", cos_sum, cos_sum_back1, cos_sum_back2);
//   normalize( cos_sum );
//   normalize( cos_sum_back1 );
//   normalize( cos_sum_back2 );
//   normalize( cos_sum_OK );

   cos_ep_cm_OK = (TH1F*)signal("cos_ep_cm_OK", cos_ep_cm, cos_back1_cm, cos_back2_cm);
//   normalize( cos_ep_cm_OK );
//   normalize( cos_back1_cm );
//   +normalize( cos_back2_cm );
//   normalize( cos_ep_cm );
   
   rapidity_OK = (TH1F*)signal("rapidity_OK", rapidity_all, rapidity_back1, rapidity_back2);
   normalize( rapidity_all );
//   normalize( rapidity_back1 );
//   normalize( rapidity_back2 );
//   normalize( rapidity_OK );

   rapidity_140_OK = (TH1F*)signal("rapidity_140_OK", rapidity_140_all, rapidity_140_back1, rapidity_140_back2);
   normalize( rapidity_140_all );
//  normalize( rapidity_140_back1 );
//   normalize( rapidity_140_back2 );
//   normalize( rapidity_140_OK );

   pt_OK = (TH1F*)signal("pt_OK", pt_all, pt_back1, pt_back2);
   normalize( pt_all );
//   normalize( pt_back1 );
//   normalize( pt_back2 );
//   normalize( pt_OK );

   pt_140_OK = (TH1F*)signal("pt_140_OK", pt_140_all, pt_140_back1, pt_140_back2);
   normalize( pt_140_all );
//   normalize( pt_140_back1 );
//   normalize( pt_140_back2 );
//   normalize( pt_140_OK );

   //normalize( sig_all_var );
     normalize( sig_all);
     normalize( Missing_mass1);
     normalize( Missing_mass);
     normalize( miss_mass_var ); 
     normalize( EpEm_inv_mass );
     //normalize( EpEm_inv_mass1 );
     //normalize( EpEm_inv_mass11 );

     // normalize( EpEm_inv_mass1_ph );
     //normalize( EpEm_inv_mass11_ph );

     //normalize( Missing_mass11 );
     //normalize( Missing_mass111 );

     //normalize( Missing_mass11_ph );
     //normalize( Missing_mass111_ph );
     normalize(epem_Mass_Brem);
     normalize(EpEm_inv_mass1_Brem);
     normalize(Missing_mass11_Brem);
     normalize(EpEm_inv_mass11_Brem);
     normalize(Missing_mass111_Brem);

     normalize(EpEm_inv_mass_Brem_reg);
     normalize(EpEm_inv_mass_Brem_nopi0cut);
     normalize(sig_all_var_Brem_nopi0cut);
     normalize(Missing_mass_Brem_nopi0cut);



/****************************************************************************************/
/********************** T E S T I N G   P A R T   E N D S *******************************/
/****************************************************************************************/
     normalize(Dilepton_OA_0deg);
     normalize(epem_Mass_Brem_0deg);

     normalize(Dilepton_OA_1deg);
     normalize(epem_Mass_Brem_1deg);

     normalize(Dilepton_OA_2deg);
     normalize(epem_Mass_Brem_2deg);

     normalize(Dilepton_OA_3deg);
     normalize(epem_Mass_Brem_3deg);

     normalize(Dilepton_OA_4deg);
     normalize(epem_Mass_Brem_4deg);

     normalize(Dilepton_OA_5deg);
     normalize(epem_Mass_Brem_5deg);

     normalize(Dilepton_OA_6deg);
     normalize(epem_Mass_Brem_6deg);

     normalize(Dilepton_OA_7deg);
     normalize(epem_Mass_Brem_7deg);

     normalize(Dilepton_OA_8deg);
     normalize(epem_Mass_Brem_8deg);

     normalize(Dilepton_OA_9deg);
     normalize(epem_Mass_Brem_9deg);

     normalize(Dilepton_OA_10deg);
     normalize(epem_Mass_Brem_10deg);




     normalize(Ep_OA);
     normalize(Em_OA);
     normalize(Ep_Mom);
     normalize(Em_Mom);

/****************************************************************************************/
/********************** T E S T I N G   P A R T   E N D S *******************************/
/****************************************************************************************/

				     

     /****************************************************************************************/

   outFileData->cd();

   tlo->Write();

   Missing_mass->Write();
   Full_signal->Write();
   Missing_mass_back1->Write();
   Full_signal_back1->Write();
   Missing_mass_back2->Add(Missing_mass_back1,1.);
   Missing_mass_back2->Write();
   Full_signal_back2->Write();

   Miss_Mass->Add(Missing_mass_back2,-1.);
   Miss_Mass->Write();   
   //   Full_signal_back2->Write();

   sig_all->Write();
   sig_all_back1->Write();
   sig_all_back2->Write();
   sig_OK->Write();

   miss_all->Write();
   miss_all_back1->Write();
   miss_all_back2->Write();
   miss_OK->Write();

   sig_all_var->Write();
   sig_all_var_back1->Write();
   sig_all_var_back2->Write();
   sig_var_OK->Write();

   SIGNAL->Write();
   CB->Write();
   SIGNIFICANCE->Write();

   sig_all_var2->Write();
   sig_all_var2_back1->Write();
   sig_all_var2_back2->Write();
   sig_var2_OK->Write();

   cos_ep->Write();
   cos_ep_back1->Write();
   cos_ep_back2->Write();
   cos_ep_OK->Write();
   cos_em->Write();
   cos_em_back1->Write();
   cos_em_back2->Write();
   cos_em_OK->Write();
   cos_sum->Write();
   cos_sum_back1->Write();
   cos_sum_back2->Write();
   cos_sum_OK->Write();
 
   cos_ep_cm_OK->Write();
   cos_ep_cm->Write();
   cos_back1_cm->Write();
   cos_back2_cm->Write();

   rapidity_all->Write();
   rapidity_back1->Write();
   rapidity_back2->Write();
   rapidity_OK->Write();
   
   rapidity_140_all->Write();
   rapidity_140_back1->Write();
   rapidity_140_back2->Write();
   rapidity_140_OK->Write();
   
   pt_all->Write();
   pt_back1->Write();
   pt_back2->Write();
   pt_OK->Write();
   
   pt_140_all->Write();
   pt_140_back1->Write();
   pt_140_back2->Write();
   pt_140_OK->Write();



   Miss_ang_dis_CosTheta->Write();
   pim_ang_dis_CosTheta->Write();
   ep_ang_dis_CosTheta->Write();
   em_ang_dis_CosTheta->Write();
   Gamma_ang_dis_CosTheta->Write();
     
   Miss_ang_dis_Theta->Write();
   pim_ang_dis_Theta->Write();
   ep_ang_dis_Theta->Write();
   em_ang_dis_Theta->Write(); 
   Gamma_ang_dis_Theta->Write();

   Missing_mass1->Write();
   Missing_mass2->Write();
   Missing_momentum->Write();
   Missing_energy->Write();

   pim_P_dis->Write();
   Gamma_P_dis->Write(); 

   EpEm_Pim_P->Write();
   EpEm_Pim_Theta->Write();


   Four_particle_MM->Write();

   Missing_TOTP->Write();
   Missing_Px->Write();
   Missing_Py->Write();
   Missing_Pz->Write();


   EpEm_inv_mass->Write();
     EpEm_inv_mass1->Write();
     EpEm_inv_mass11->Write(); 
     Missing_mass11->Write();
     Missing_mass111->Write();
     

     epem_Mass->Write();
     epem_Mass_Brem->Write();
     sig_all_MeV->Write();


     mom_miss->Write(); //Mom correction test
     

     EpEm_inv_mass1_Brem->Write();
     Missing_mass11_Brem->Write();
     EpEm_inv_mass11_Brem->Write();
     Missing_mass111_Brem->Write();
     sig_all_var_Brem->Write();
     EpEm_inv_mass1_ph->Write();
     EpEm_inv_mass11_ph->Write();
     Missing_mass11_ph->Write();
     Missing_mass111_ph->Write();
     
     EpEm_inv_mass1_M2->Write();
     EpEm_inv_mass11_M2->Write();
     Missing_mass11_M2->Write();
     Missing_mass111_M2->Write();
     miss_mass_var->Write();

     EpEm_inv_mass_Brem_reg->Write();
     EpEm_inv_mass_Brem_nopi0cut->Write();
     sig_all_var_Brem_nopi0cut->Write();
     Missing_mass_Brem_nopi0cut->Write();



/****************************************************************************************/
/*************************** T E S T I N G   P A R T ************************************/
/****************************************************************************************/

     Dilepton_OA_10deg->Write();
     epem_Mass_Brem_10deg->Write();

     Dilepton_OA_1deg->Write();
     epem_Mass_Brem_1deg->Write();

     Dilepton_OA_2deg->Write();
     epem_Mass_Brem_2deg->Write();

     Dilepton_OA_3deg->Write();
     epem_Mass_Brem_3deg->Write();

     Dilepton_OA_4deg->Write();
     epem_Mass_Brem_4deg->Write();

     Dilepton_OA_5deg->Write();
     epem_Mass_Brem_5deg->Write();

     Dilepton_OA_6deg->Write();
     epem_Mass_Brem_6deg->Write();

     Dilepton_OA_7deg->Write();
     epem_Mass_Brem_7deg->Write();

     Dilepton_OA_8deg->Write();
     epem_Mass_Brem_8deg->Write();

     Dilepton_OA_9deg->Write();
     epem_Mass_Brem_9deg->Write();

     Dilepton_OA_10deg->Write();
     epem_Mass_Brem_10deg->Write();
     
     Ep_OA->Write();
     Em_OA->Write();
     Ep_Mom->Write();
     Em_Mom->Write();


     
     outFileData->Close();
}



