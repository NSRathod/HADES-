#ifndef DATA_H
#define DATA_H

#include <TROOT.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCutG.h>
#include <TCut.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include "hntuple.h"
#include "HFilter.h"



/**************************** global histograms repository ***********************************/

namespace PATData {

   extern TFile         *outFileData;

   extern HNtuple       *tlo;

   extern HFilter       *filter;
   extern float          EFF, ACC, EFF_B,EFF_M2;

   extern double         d_pion_mom;
   extern TH2F          *mom_miss;
   extern TH1F          *oa_corr, *Missing_TOTP, *Missing_Px, *Missing_Py, *Missing_Pz, *Missing_mass2, *epem_Mass, *miss_mass_var;

   extern TH2F          *EpEm_Pim_P, *EpEm_Pim_Theta;
   extern TH1F          *Gamma_P_dis, *pim_P_dis, *Four_particle_MM, *epem_Mass_Brem;
   extern TH1F          *Missing_mass, *Full_signal, *Missing_mass_back1, *Full_signal_back1, *Missing_mass_back2, *Full_signal_back2, *Miss_Mass;
   extern TH1F          *Miss_ang_dis_CosTheta, *pim_ang_dis_CosTheta, *ep_ang_dis_CosTheta, *em_ang_dis_CosTheta, *Gamma_ang_dis_CosTheta;
   extern TH1F          *Miss_ang_dis_Theta, *pim_ang_dis_Theta, *ep_ang_dis_Theta, *em_ang_dis_Theta, *Gamma_ang_dis_Theta;
   extern TH1F          *Missing_mass1, *Missing_momentum, *Missing_energy;
   extern TH1F          *EpEm_inv_mass, *EpEm_inv_mass1, *EpEm_inv_mass11, *Missing_mass11, *Missing_mass111,*sig_all_MeV;
   extern TH1F          *EpEm_inv_mass1_Brem, *Missing_mass11_Brem, *EpEm_inv_mass11_Brem, *Missing_mass111_Brem, *sig_all_var_Brem;
   extern TH1F          *EpEm_inv_mass1_ph, *EpEm_inv_mass11_ph, *Missing_mass11_ph, *Missing_mass111_ph;
   extern TH1F          *EpEm_inv_mass1_M2, *EpEm_inv_mass11_M2, *Missing_mass11_M2, *Missing_mass111_M2, *EpEm_inv_mass_Brem_reg;
   extern TH1F          *EpEm_inv_mass_Brem_nopi0cut, *sig_all_var_Brem_nopi0cut, *Missing_mass_Brem_nopi0cut;
   extern TH1F          *Dilepton_OA_0deg, *epem_Mass_Brem_0deg, *Dilepton_OA_1deg, *epem_Mass_Brem_1deg, *Dilepton_OA_2deg, *epem_Mass_Brem_2deg, *Dilepton_OA_3deg, *epem_Mass_Brem_3deg, *Dilepton_OA_4deg, *epem_Mass_Brem_4deg, *Dilepton_OA_5deg, *epem_Mass_Brem_5deg, *Dilepton_OA_6deg, *epem_Mass_Brem_6deg, *Dilepton_OA_7deg, *epem_Mass_Brem_7deg, *Dilepton_OA_8deg, *epem_Mass_Brem_8deg, *Dilepton_OA_9deg, *epem_Mass_Brem_9deg, *Dilepton_OA_10deg, *epem_Mass_Brem_10deg, *Ep_OA, *Em_OA, *Ep_Mom, *Em_Mom;

   
   extern TH1F          *SIGNAL, *CB, *SIGNIFICANCE;
   extern TH1F          *sig_all, *sig_all_back1, *sig_all_back2, *sig_OK;
   extern TH1F          *miss_all, *miss_all_back1, *miss_all_back2, *miss_OK;
   extern TH1F          *sig_all_var, *sig_all_var_back1, *sig_all_var_back2, *sig_var_OK;
   extern TH1F          *sig_all_var2, *sig_all_var2_back1, *sig_all_var2_back2, *sig_var2_OK;
   extern TH1F 		*cos_ep, *cos_em, *cos_sum;
   extern TH1F 		*cos_ep_back1, *cos_em_back1, *cos_sum_back1;
   extern TH1F 		*cos_ep_back2, *cos_em_back2, *cos_sum_back2;
   extern TH1F		*cos_ep_OK, *cos_em_OK, *cos_sum_OK;
   extern TH1F		*cos_ep_cm_OK, *cos_back1_cm, *cos_back2_cm, *cos_ep_cm;
   extern TH1F		*rapidity_all, *rapidity_back1, *rapidity_back2, *rapidity_OK;
   extern TH1F		*rapidity_140_all, *rapidity_140_back1, *rapidity_140_back2, *rapidity_140_OK;
   extern TH1F		*pt_all, *pt_back1, *pt_back2, *pt_OK;
   extern TH1F		*pt_140_all, *pt_140_back1, *pt_140_back2, *pt_140_OK;

   extern TFile *file1_cuts, *file2_cuts;

   extern TCutG *pEpS0, *pEpS1, *pEmS0, *pEmS1;
   extern TCutG *pEm1S0, *pEm1S1, *pEm2S0, *pEm2S1;
   extern TCutG *pEp1S0, *pEp1S1, *pEp2S0, *pEp2S1;
   extern TCutG *pvertex_xy, *pvertex_xz, *pvertex_yz;

   extern Bool_t NoLeptonE1;
   extern Bool_t NoHadronE1;
   extern Bool_t NoLeptonE2;
   extern Bool_t NoHadronE2;

   extern Bool_t Electron;
   extern Bool_t Positron;

   extern Bool_t Electron1;
   extern Bool_t Electron2;
   extern Bool_t Positron1;
   extern Bool_t Positron2;

   extern Bool_t ElectronPositron;
   extern Bool_t ElectronElectron;
   extern Bool_t PositronPositron;

   extern TLorentzVector *pim;
   extern TLorentzVector *e1;
   extern TLorentzVector *e2;
   extern TLorentzVector *beam;
   extern TLorentzVector *beam_PT;
   extern TLorentzVector *proj;
   extern TLorentzVector *proj_PT;
   extern TLorentzVector *targ;
   extern TLorentzVector *gammae1e2;
   extern TLorentzVector *e1e2;
   extern TLorentzVector *e1e2_miss;
   extern TLorentzVector *e1_delta;
   extern TLorentzVector *e2_delta;
   extern TLorentzVector *pimepem;
   extern TLorentzVector *pimepem_miss;
   extern TLorentzVector *pimepem_miss_PT;
   extern TLorentzVector *ppimepem_miss;
   extern TLorentzVector *proj_var;
   extern TLorentzVector *beam_var;
   
   extern TLorentzVector *pim_B;
   extern TLorentzVector *e1_B;
   extern TLorentzVector *e2_B;
   extern TLorentzVector *gammae1e2_B;
   extern TLorentzVector *e1e2_B;
   extern TLorentzVector *e1e2_miss_B;
   extern TLorentzVector *e1_delta_B;
   extern TLorentzVector *e2_delta_B;
   extern TLorentzVector *pimepem_B;
   extern TLorentzVector *pimepem_miss_B;
   extern TLorentzVector *pimepem_miss_PT_B;
   extern TLorentzVector *ppimepem_miss_B;

   extern TLorentzVector *beam_B;
   extern TLorentzVector *proj_B;
   extern TLorentzVector *beam_PT_B;
   extern TLorentzVector *beam_var_B;
   extern TLorentzVector *proj_var_B;
   extern TLorentzVector *proj_PT_B;


   extern Int_t insideTarget;

   extern Int_t insideEmS0;
   extern Int_t insideEmS1;
   extern Int_t insideEpS0;
   extern Int_t insideEpS1;

   extern Int_t insideEm1S0;
   extern Int_t insideEm1S1;
   extern Int_t insideEm2S0;
   extern Int_t insideEm2S1;

   extern Int_t insideEp1S0;
   extern Int_t insideEp1S1;
   extern Int_t insideEp2S0;
   extern Int_t insideEp2S1;

   extern const double D2R;
   extern const double R2D;


   /************************* M E T H O D S *************************************/

   double openingangle(const TLorentzVector& a, const TLorentzVector& b);
   double openingangle(const TVector3& a, const TVector3& b);
   void normalize(TH1* hist);
   TH1* signal(const char* name, TH1* hist, TH1* back1, TH1* back2);

}

/*********************************************************************************************/


#endif // DATA_H
