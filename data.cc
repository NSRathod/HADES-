#include "data.h"

/**************************** global histograms repository ***********************************/

namespace PATData {

   TFile         *outFileData;

   HNtuple       *tlo;

   HFilter       *filter;
  float         EFF, ACC, EFF_B, EFF_M2;

   double        d_pion_mom;
   TH2F          *mom_miss;
  TH1F          *oa_corr, *Missing_TOTP, *Missing_Px, *Missing_Py, *Missing_Pz, *Missing_mass2, *epem_Mass, *miss_mass_var;

   TH2F          *EpEm_Pim_P, *EpEm_Pim_Theta;
  TH1F          *Gamma_P_dis, *pim_P_dis,*Four_particle_MM, *epem_Mass_Brem;
   TH1F          *Missing_mass, *Full_signal, *Missing_mass_back1, *Full_signal_back1, *Missing_mass_back2, *Full_signal_back2, *Miss_Mass;
   TH1F          *Miss_ang_dis_CosTheta, *pim_ang_dis_CosTheta, *ep_ang_dis_CosTheta, *em_ang_dis_CosTheta, *Gamma_ang_dis_CosTheta;
   TH1F          *Miss_ang_dis_Theta, *pim_ang_dis_Theta, *ep_ang_dis_Theta, *em_ang_dis_Theta, *Gamma_ang_dis_Theta;
   TH1F          *Missing_mass1, *Missing_momentum, *Missing_energy;
   TH1F          *EpEm_inv_mass, *EpEm_inv_mass1, *EpEm_inv_mass11, *Missing_mass11, *Missing_mass111,*sig_all_MeV;

   TH1F          *EpEm_inv_mass1_Brem, *Missing_mass11_Brem, *EpEm_inv_mass11_Brem, *Missing_mass111_Brem, *sig_all_var_Brem;
  TH1F          *EpEm_inv_mass1_ph, *EpEm_inv_mass11_ph, *Missing_mass11_ph, *Missing_mass111_ph;
  TH1F          *EpEm_inv_mass1_M2, *EpEm_inv_mass11_M2, *Missing_mass11_M2, *Missing_mass111_M2, *EpEm_inv_mass_Brem_reg;
  TH1F          *EpEm_inv_mass_Brem_nopi0cut, *sig_all_var_Brem_nopi0cut, *Missing_mass_Brem_nopi0cut;

  TH1F          *Dilepton_OA_0deg, *epem_Mass_Brem_0deg, *Dilepton_OA_1deg, *epem_Mass_Brem_1deg, *Dilepton_OA_2deg, *epem_Mass_Brem_2deg, *Dilepton_OA_3deg, *epem_Mass_Brem_3deg, *Dilepton_OA_4deg, *epem_Mass_Brem_4deg, *Dilepton_OA_5deg, *epem_Mass_Brem_5deg, *Dilepton_OA_6deg, *epem_Mass_Brem_6deg, *Dilepton_OA_7deg, *epem_Mass_Brem_7deg, *Dilepton_OA_8deg, *epem_Mass_Brem_8deg, *Dilepton_OA_9deg, *epem_Mass_Brem_9deg, *Dilepton_OA_10deg, *epem_Mass_Brem_10deg, *Ep_OA, *Em_OA, *Ep_Mom, *Em_Mom;
  
   TH1F          *SIGNAL, *CB, *SIGNIFICANCE;
   TH1F          *sig_all, *sig_all_back1, *sig_all_back2, *sig_OK;
   TH1F          *miss_all, *miss_all_back1, *miss_all_back2, *miss_OK;
   TH1F          *sig_all_var, *sig_all_var_back1, *sig_all_var_back2, *sig_var_OK;
   TH1F          *sig_all_var2, *sig_all_var2_back1, *sig_all_var2_back2, *sig_var2_OK;
   TH1F          *cos_ep, *cos_em, *cos_sum;
   TH1F          *cos_ep_cm, *cos_back1_cm, *cos_back2_cm,*cos_ep_cm_OK;
   TH1F          *cos_ep_back1, *cos_em_back1, *cos_sum_back1;
   TH1F          *cos_ep_back2, *cos_em_back2, *cos_sum_back2;
   TH1F          *cos_ep_OK, *cos_em_OK, *cos_sum_OK;
   TH1F          *rapidity_all, *rapidity_back1, *rapidity_back2, *rapidity_OK;
   TH1F          *rapidity_140_all, *rapidity_140_back1, *rapidity_140_back2, *rapidity_140_OK;
   TH1F          *pt_all, *pt_back1, *pt_back2, *pt_OK;
   TH1F          *pt_140_all, *pt_140_back1, *pt_140_back2, *pt_140_OK;


   TFile *file1_cuts, *file2_cuts;

   TCutG *pEpS0, *pEpS1, *pEmS0, *pEmS1;
   TCutG *pEm1S0, *pEm1S1, *pEm2S0, *pEm2S1;
   TCutG *pEp1S0, *pEp1S1, *pEp2S0, *pEp2S1;
   TCutG *pvertex_xy, *pvertex_xz, *pvertex_yz;

   Bool_t NoLeptonE1;
   Bool_t NoHadronE1;
   Bool_t NoLeptonE2;
   Bool_t NoHadronE2;

   Bool_t Electron;
   Bool_t Positron;

   Bool_t Electron1;
   Bool_t Electron2;
   Bool_t Positron1;
   Bool_t Positron2;

   Bool_t ElectronPositron;
   Bool_t ElectronElectron;
   Bool_t PositronPositron;

   TLorentzVector *pim;
   TLorentzVector *e1;
   TLorentzVector *e2;
   TLorentzVector *beam;
   TLorentzVector *beam_PT;
   TLorentzVector *proj;
   TLorentzVector *proj_PT;
   TLorentzVector *targ;
   TLorentzVector *gammae1e2;
   TLorentzVector *e1e2;
   TLorentzVector *e1e2_miss;
   TLorentzVector *e1_delta;
   TLorentzVector *e2_delta;
   TLorentzVector *pimepem;
   TLorentzVector *pimepem_miss;
   TLorentzVector *pimepem_miss_PT;
   TLorentzVector *ppimepem_miss;
   TLorentzVector *proj_var;
   TLorentzVector *beam_var;

   TLorentzVector *pim_B;
   TLorentzVector *e1_B;
   TLorentzVector *e2_B;
   TLorentzVector *gammae1e2_B;
   TLorentzVector *e1e2_B;
   TLorentzVector *e1e2_miss_B;
   TLorentzVector *e1_delta_B;
   TLorentzVector *e2_delta_B;
   TLorentzVector *pimepem_B;
   TLorentzVector *pimepem_miss_B;
   TLorentzVector *pimepem_miss_PT_B;
   TLorentzVector *ppimepem_miss_B;

   TLorentzVector *beam_B;
   TLorentzVector *proj_B;
   TLorentzVector *beam_PT_B;
   TLorentzVector *beam_var_B;
   TLorentzVector *proj_var_B;
   TLorentzVector *proj_PT_B;


  
   Int_t insideTarget;

   Int_t insideEmS0;
   Int_t insideEmS1;
   Int_t insideEpS0;
   Int_t insideEpS1;

   Int_t insideEm1S0;
   Int_t insideEm1S1;
   Int_t insideEm2S0;
   Int_t insideEm2S1;

   Int_t insideEp1S0;
   Int_t insideEp1S1;
   Int_t insideEp2S0;
   Int_t insideEp2S1;

   const double D2R = 1.74532925199432955e-02;
   const double R2D = 57.2957795130823229;


   /************************* M E T H O D S *************************************/

   double openingangle(const TLorentzVector& a, const TLorentzVector& b)
   {
         return TMath::ACos( (a.Px()*b.Px() + a.Py()*b.Py() +  a.Pz()*b.Pz() ) / ( a.Vect().Mag() * b.Vect().Mag() ) );
   }

   double openingangle(const TVector3& a, const TVector3& b)
   {
         return TMath::ACos( (a.Px()*b.Px() + a.Py()*b.Py() +  a.Pz()*b.Pz() ) / ( a.Mag() * b.Mag() ) );
   }


   void normalize(TH1* hist)
   {
      for (Int_t j=1; j<hist->GetNbinsX()+1; ++j)
      {
         hist->SetBinContent( j, hist->GetBinContent(j) / hist->GetBinWidth(j) );
//         hist->SetBinError( j, TMath::Sqrt( hist->GetBinContent(j) ) );
         hist->SetBinError( j, hist->GetBinError(j) / hist->GetBinWidth(j) );
      }
   }

   TH1* signal(const char* name, TH1* hist, TH1* back1, TH1* back2)
   {
      TH1 *ptr = (TH1*)hist->Clone(name);
      for (Int_t j=1; j<hist->GetNbinsX()+1; ++j)
      {
            //ptr->SetBinContent(j, hist->GetBinContent(j) - back1->GetBinContent(j) - back2->GetBinContent(j));
            ptr->SetBinContent(j, hist->GetBinContent(j) - 2*TMath::Sqrt(back1->GetBinContent(j)*back2->GetBinContent(j)));
            ptr->SetBinError(j, TMath::Sqrt( hist->GetBinError(j)*hist->GetBinError(j) + 
                                             back1->GetBinError(j)*back1->GetBinError(j) + 
                                             back2->GetBinError(j)*back2->GetBinError(j) ));
      }

      return ptr;
   }


}

/*********************************************************************************************/

