//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jul 11 12:27:31 2019 by ROOT version 5.34/34
// from TTree PimEpEm_ID/PimEpEm_ID
// found on file: hgeantout_PWA_momentums-ppimpi0_ev_weight_ALL.root
//////////////////////////////////////////////////////////

#ifndef PimEpEm_ID_h
#define PimEpEm_ID_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class PimEpEm_ID {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Float_t         eVertClust_chi2;
   Float_t         eVertClust_x;
   Float_t         eVertClust_y;
   Float_t         eVertClust_z;
   Float_t         eVertReco_chi2;
   Float_t         eVertReco_x;
   Float_t         eVertReco_y;
   Float_t         eVertReco_z;
   Float_t         eVert_chi2;
   Float_t         eVert_x;
   Float_t         eVert_y;
   Float_t         eVert_z;
   Float_t         em_beta;
   Float_t         em_beta_new;
   Float_t         em_btChargeRing;
   Float_t         em_btChargeSum;
   Float_t         em_btChi2;
   Float_t         em_btClusters;
   Float_t         em_btMaxima;
   Float_t         em_btMaximaCharge;
   Float_t         em_btMaximaChargeShared;
   Float_t         em_btMaximaChargeSharedFragment;
   Float_t         em_btMaximaShared;
   Float_t         em_btMaximaSharedFragment;
   Float_t         em_btMeanDist;
   Float_t         em_btNearbyMaxima;
   Float_t         em_btNearbyMaximaShared;
   Float_t         em_btPadsClus;
   Float_t         em_btPadsRing;
   Float_t         em_btRingMatrix;
   Float_t         em_dedx_mdc;
   Float_t         em_dedx_tof;
   Float_t         em_id;
   Float_t         em_isBT;
   Float_t         em_isOffVertexClust;
   Float_t         em_isPrimaryVertex;
   Float_t         em_isUsedVertex;
   Float_t         em_isring;
   Float_t         em_isringmdc;
   Float_t         em_isringnomatch;
   Float_t         em_isringtrack;
   Float_t         em_kIsLepton;
   Float_t         em_kIsUsed;
   Float_t         em_mdcinnerchi2;
   Float_t         em_mdcouterchi2;
   Float_t         em_oa_hadr;
   Float_t         em_oa_lept;
   Float_t         em_p;
   Float_t         em_p_corr_em;
   Float_t         em_p_corr_ep;
   Float_t         em_p_corr_p;
   Float_t         em_p_corr_pim;
   Float_t         em_p_corr_pip;
   Float_t         em_phi;
   Float_t         em_phi_rich;
   Float_t         em_pid;
   Float_t         em_q;
   Float_t         em_r;
   Float_t         em_resolution;
   Float_t         em_resoultion;
   Float_t         em_rich_amp;
   Float_t         em_rich_centr;
   Float_t         em_rich_houtra;
   Float_t         em_rich_padnum;
   Float_t         em_rich_patmat;
   Float_t         em_rkchi2;
   Float_t         em_sector;
   Float_t         em_shw_sum0;
   Float_t         em_shw_sum1;
   Float_t         em_shw_sum2;
   Float_t         em_sim_geninfo;
   Float_t         em_sim_geninfo1;
   Float_t         em_sim_geninfo2;
   Float_t         em_sim_genweight;
   Float_t         em_sim_grandparentid;
   Float_t         em_sim_id;
   Float_t         em_sim_iscommon;
   Float_t         em_sim_isghost;
   Float_t         em_sim_mediumid;
   Float_t         em_sim_p;
   Float_t         em_sim_parentid;
   Float_t         em_sim_processid;
   Float_t         em_sim_px;
   Float_t         em_sim_py;
   Float_t         em_sim_pz;
   Float_t         em_sim_vertexx;
   Float_t         em_sim_vertexy;
   Float_t         em_sim_vertexz;
   Float_t         em_system;
   Float_t         em_theta;
   Float_t         em_theta_rich;
   Float_t         em_tof_mom;
   Float_t         em_tof_new;
   Float_t         em_tof_rec;
   Float_t         em_track_length;
   Float_t         em_tracklength;
   Float_t         em_z;
   Float_t         ep_beta;
   Float_t         ep_beta_new;
   Float_t         ep_btChargeRing;
   Float_t         ep_btChargeSum;
   Float_t         ep_btChi2;
   Float_t         ep_btClusters;
   Float_t         ep_btMaxima;
   Float_t         ep_btMaximaCharge;
   Float_t         ep_btMaximaChargeShared;
   Float_t         ep_btMaximaChargeSharedFragment;
   Float_t         ep_btMaximaShared;
   Float_t         ep_btMaximaSharedFragment;
   Float_t         ep_btMeanDist;
   Float_t         ep_btNearbyMaxima;
   Float_t         ep_btNearbyMaximaShared;
   Float_t         ep_btPadsClus;
   Float_t         ep_btPadsRing;
   Float_t         ep_btRingMatrix;
   Float_t         ep_dedx_mdc;
   Float_t         ep_dedx_tof;
   Float_t         ep_id;
   Float_t         ep_isBT;
   Float_t         ep_isOffVertexClust;
   Float_t         ep_isPrimaryVertex;
   Float_t         ep_isUsedVertex;
   Float_t         ep_isring;
   Float_t         ep_isringmdc;
   Float_t         ep_isringnomatch;
   Float_t         ep_isringtrack;
   Float_t         ep_kIsLepton;
   Float_t         ep_kIsUsed;
   Float_t         ep_mdcinnerchi2;
   Float_t         ep_mdcouterchi2;
   Float_t         ep_oa_hadr;
   Float_t         ep_oa_lept;
   Float_t         ep_p;
   Float_t         ep_p_corr_em;
   Float_t         ep_p_corr_ep;
   Float_t         ep_p_corr_p;
   Float_t         ep_p_corr_pim;
   Float_t         ep_p_corr_pip;
   Float_t         ep_phi;
   Float_t         ep_phi_rich;
   Float_t         ep_pid;
   Float_t         ep_q;
   Float_t         ep_r;
   Float_t         ep_resolution;
   Float_t         ep_resoultion;
   Float_t         ep_rich_amp;
   Float_t         ep_rich_centr;
   Float_t         ep_rich_houtra;
   Float_t         ep_rich_padnum;
   Float_t         ep_rich_patmat;
   Float_t         ep_rkchi2;
   Float_t         ep_sector;
   Float_t         ep_shw_sum0;
   Float_t         ep_shw_sum1;
   Float_t         ep_shw_sum2;
   Float_t         ep_sim_geninfo;
   Float_t         ep_sim_geninfo1;
   Float_t         ep_sim_geninfo2;
   Float_t         ep_sim_genweight;
   Float_t         ep_sim_grandparentid;
   Float_t         ep_sim_id;
   Float_t         ep_sim_iscommon;
   Float_t         ep_sim_isghost;
   Float_t         ep_sim_mediumid;
   Float_t         ep_sim_p;
   Float_t         ep_sim_parentid;
   Float_t         ep_sim_processid;
   Float_t         ep_sim_px;
   Float_t         ep_sim_py;
   Float_t         ep_sim_pz;
   Float_t         ep_sim_vertexx;
   Float_t         ep_sim_vertexy;
   Float_t         ep_sim_vertexz;
   Float_t         ep_system;
   Float_t         ep_theta;
   Float_t         ep_theta_rich;
   Float_t         ep_tof_mom;
   Float_t         ep_tof_new;
   Float_t         ep_tof_rec;
   Float_t         ep_track_length;
   Float_t         ep_tracklength;
   Float_t         ep_z;
   Float_t         event;
   Float_t         evtFlashMDC;
   Float_t         evtGood;
   Float_t         evtGoodTrig;
   Float_t         evtLepMult;
   Float_t         evtMDCMult;
   Float_t         evtPileupMDC;
   Float_t         evtPileupMeta;
   Float_t         evtPileupStart;
   Float_t         evtStart;
   Float_t         evtVertCand;
   Float_t         evtVertClust;
   Float_t         fw_beta_1;
   Float_t         fw_beta_2;
   Float_t         fw_beta_3;
   Float_t         fw_charge_1;
   Float_t         fw_charge_2;
   Float_t         fw_charge_3;
   Float_t         fw_cluster_mult;
   Float_t         fw_distance_1;
   Float_t         fw_distance_2;
   Float_t         fw_distance_3;
   Float_t         fw_mult;
   Float_t         fw_p_1;
   Float_t         fw_p_2;
   Float_t         fw_p_3;
   Float_t         fw_phi_1;
   Float_t         fw_phi_2;
   Float_t         fw_phi_3;
   Float_t         fw_sim_geninfo1_1;
   Float_t         fw_sim_geninfo1_2;
   Float_t         fw_sim_geninfo1_3;
   Float_t         fw_sim_geninfo2_1;
   Float_t         fw_sim_geninfo2_2;
   Float_t         fw_sim_geninfo2_3;
   Float_t         fw_sim_geninfo_1;
   Float_t         fw_sim_geninfo_2;
   Float_t         fw_sim_geninfo_3;
   Float_t         fw_sim_genweight_1;
   Float_t         fw_sim_genweight_2;
   Float_t         fw_sim_genweight_3;
   Float_t         fw_sim_id_1;
   Float_t         fw_sim_id_2;
   Float_t         fw_sim_id_3;
   Float_t         fw_sim_p_1;
   Float_t         fw_sim_p_2;
   Float_t         fw_sim_p_3;
   Float_t         fw_sim_parentid_1;
   Float_t         fw_sim_parentid_2;
   Float_t         fw_sim_parentid_3;
   Float_t         fw_size_1;
   Float_t         fw_size_2;
   Float_t         fw_size_3;
   Float_t         fw_spectator_1;
   Float_t         fw_spectator_2;
   Float_t         fw_spectator_3;
   Float_t         fw_theta_1;
   Float_t         fw_theta_2;
   Float_t         fw_theta_3;
   Float_t         fw_time_1;
   Float_t         fw_time_2;
   Float_t         fw_time_3;
   Float_t         fw_time_min_1;
   Float_t         fw_time_min_2;
   Float_t         fw_time_min_3;
   Float_t         fw_tracknr_1;
   Float_t         fw_tracknr_12;
   Float_t         fw_tracknr_2;
   Float_t         fw_tracknr_22;
   Float_t         fw_tracknr_3;
   Float_t         fw_tracknr_32;
   Float_t         fw_x_lab_1;
   Float_t         fw_x_lab_2;
   Float_t         fw_x_lab_3;
   Float_t         fw_y_lab_1;
   Float_t         fw_y_lab_2;
   Float_t         fw_y_lab_3;
   Float_t         fw_z_lab_1;
   Float_t         fw_z_lab_2;
   Float_t         fw_z_lab_3;
   Float_t         geant_sim_genweight_1;
   Float_t         geant_sim_genweight_2;
   Float_t         geant_sim_id_1;
   Float_t         geant_sim_id_2;
   Float_t         geant_sim_p_1;
   Float_t         geant_sim_p_2;
   Float_t         geant_sim_px_1;
   Float_t         geant_sim_px_2;
   Float_t         geant_sim_py_1;
   Float_t         geant_sim_py_2;
   Float_t         geant_sim_pz_1;
   Float_t         geant_sim_pz_2;
   Float_t         hneg_mult;
   Float_t         hpos_mult;
   Float_t         isBest;
   Float_t         lneg_mult;
   Float_t         lpos_mult;
   Float_t         pim_beta;
   Float_t         pim_beta_new;
   Float_t         pim_btChargeRing;
   Float_t         pim_btChargeSum;
   Float_t         pim_btChi2;
   Float_t         pim_btClusters;
   Float_t         pim_btMaxima;
   Float_t         pim_btMaximaCharge;
   Float_t         pim_btMaximaChargeShared;
   Float_t         pim_btMaximaChargeSharedFragment;
   Float_t         pim_btMaximaShared;
   Float_t         pim_btMaximaSharedFragment;
   Float_t         pim_btMeanDist;
   Float_t         pim_btNearbyMaxima;
   Float_t         pim_btNearbyMaximaShared;
   Float_t         pim_btPadsClus;
   Float_t         pim_btPadsRing;
   Float_t         pim_btRingMatrix;
   Float_t         pim_dedx_mdc;
   Float_t         pim_dedx_tof;
   Float_t         pim_id;
   Float_t         pim_isBT;
   Float_t         pim_isOffVertexClust;
   Float_t         pim_isPrimaryVertex;
   Float_t         pim_isUsedVertex;
   Float_t         pim_isring;
   Float_t         pim_isringmdc;
   Float_t         pim_isringnomatch;
   Float_t         pim_isringtrack;
   Float_t         pim_kIsLepton;
   Float_t         pim_kIsUsed;
   Float_t         pim_mdcinnerchi2;
   Float_t         pim_mdcouterchi2;
   Float_t         pim_oa_hadr;
   Float_t         pim_oa_lept;
   Float_t         pim_p;
   Float_t         pim_p_corr_em;
   Float_t         pim_p_corr_ep;
   Float_t         pim_p_corr_p;
   Float_t         pim_p_corr_pim;
   Float_t         pim_p_corr_pip;
   Float_t         pim_phi;
   Float_t         pim_phi_rich;
   Float_t         pim_pid;
   Float_t         pim_q;
   Float_t         pim_r;
   Float_t         pim_resolution;
   Float_t         pim_resoultion;
   Float_t         pim_rich_amp;
   Float_t         pim_rich_centr;
   Float_t         pim_rich_houtra;
   Float_t         pim_rich_padnum;
   Float_t         pim_rich_patmat;
   Float_t         pim_rkchi2;
   Float_t         pim_sector;
   Float_t         pim_shw_sum0;
   Float_t         pim_shw_sum1;
   Float_t         pim_shw_sum2;
   Float_t         pim_sim_geninfo;
   Float_t         pim_sim_geninfo1;
   Float_t         pim_sim_geninfo2;
   Float_t         pim_sim_genweight;
   Float_t         pim_sim_grandparentid;
   Float_t         pim_sim_id;
   Float_t         pim_sim_iscommon;
   Float_t         pim_sim_isghost;
   Float_t         pim_sim_mediumid;
   Float_t         pim_sim_p;
   Float_t         pim_sim_parentid;
   Float_t         pim_sim_processid;
   Float_t         pim_sim_px;
   Float_t         pim_sim_py;
   Float_t         pim_sim_pz;
   Float_t         pim_sim_vertexx;
   Float_t         pim_sim_vertexy;
   Float_t         pim_sim_vertexz;
   Float_t         pim_system;
   Float_t         pim_theta;
   Float_t         pim_theta_rich;
   Float_t         pim_tof_mom;
   Float_t         pim_tof_new;
   Float_t         pim_tof_rec;
   Float_t         pim_track_length;
   Float_t         pim_tracklength;
   Float_t         pim_z;
   Float_t         runnumber;
   Float_t         totalmult;
   Float_t         trigbit;
   Float_t         trigdec;
   Float_t         trigdownscale;
   Float_t         trigdownscaleflag;

   // List of branches
   TBranch        *b_eVertClust_chi2;   //!
   TBranch        *b_eVertClust_x;   //!
   TBranch        *b_eVertClust_y;   //!
   TBranch        *b_eVertClust_z;   //!
   TBranch        *b_eVertReco_chi2;   //!
   TBranch        *b_eVertReco_x;   //!
   TBranch        *b_eVertReco_y;   //!
   TBranch        *b_eVertReco_z;   //!
   TBranch        *b_eVert_chi2;   //!
   TBranch        *b_eVert_x;   //!
   TBranch        *b_eVert_y;   //!
   TBranch        *b_eVert_z;   //!
   TBranch        *b_em_beta;   //!
   TBranch        *b_em_beta_new;   //!
   TBranch        *b_em_btChargeRing;   //!
   TBranch        *b_em_btChargeSum;   //!
   TBranch        *b_em_btChi2;   //!
   TBranch        *b_em_btClusters;   //!
   TBranch        *b_em_btMaxima;   //!
   TBranch        *b_em_btMaximaCharge;   //!
   TBranch        *b_em_btMaximaChargeShared;   //!
   TBranch        *b_em_btMaximaChargeSharedFragment;   //!
   TBranch        *b_em_btMaximaShared;   //!
   TBranch        *b_em_btMaximaSharedFragment;   //!
   TBranch        *b_em_btMeanDist;   //!
   TBranch        *b_em_btNearbyMaxima;   //!
   TBranch        *b_em_btNearbyMaximaShared;   //!
   TBranch        *b_em_btPadsClus;   //!
   TBranch        *b_em_btPadsRing;   //!
   TBranch        *b_em_btRingMatrix;   //!
   TBranch        *b_em_dedx_mdc;   //!
   TBranch        *b_em_dedx_tof;   //!
   TBranch        *b_em_id;   //!
   TBranch        *b_em_isBT;   //!
   TBranch        *b_em_isOffVertexClust;   //!
   TBranch        *b_em_isPrimaryVertex;   //!
   TBranch        *b_em_isUsedVertex;   //!
   TBranch        *b_em_isring;   //!
   TBranch        *b_em_isringmdc;   //!
   TBranch        *b_em_isringnomatch;   //!
   TBranch        *b_em_isringtrack;   //!
   TBranch        *b_em_kIsLepton;   //!
   TBranch        *b_em_kIsUsed;   //!
   TBranch        *b_em_mdcinnerchi2;   //!
   TBranch        *b_em_mdcouterchi2;   //!
   TBranch        *b_em_oa_hadr;   //!
   TBranch        *b_em_oa_lept;   //!
   TBranch        *b_em_p;   //!
   TBranch        *b_em_p_corr_em;   //!
   TBranch        *b_em_p_corr_ep;   //!
   TBranch        *b_em_p_corr_p;   //!
   TBranch        *b_em_p_corr_pim;   //!
   TBranch        *b_em_p_corr_pip;   //!
   TBranch        *b_em_phi;   //!
   TBranch        *b_em_phi_rich;   //!
   TBranch        *b_em_pid;   //!
   TBranch        *b_em_q;   //!
   TBranch        *b_em_r;   //!
   TBranch        *b_em_resolution;   //!
   TBranch        *b_em_resoultion;   //!
   TBranch        *b_em_rich_amp;   //!
   TBranch        *b_em_rich_centr;   //!
   TBranch        *b_em_rich_houtra;   //!
   TBranch        *b_em_rich_padnum;   //!
   TBranch        *b_em_rich_patmat;   //!
   TBranch        *b_em_rkchi2;   //!
   TBranch        *b_em_sector;   //!
   TBranch        *b_em_shw_sum0;   //!
   TBranch        *b_em_shw_sum1;   //!
   TBranch        *b_em_shw_sum2;   //!
   TBranch        *b_em_sim_geninfo;   //!
   TBranch        *b_em_sim_geninfo1;   //!
   TBranch        *b_em_sim_geninfo2;   //!
   TBranch        *b_em_sim_genweight;   //!
   TBranch        *b_em_sim_grandparentid;   //!
   TBranch        *b_em_sim_id;   //!
   TBranch        *b_em_sim_iscommon;   //!
   TBranch        *b_em_sim_isghost;   //!
   TBranch        *b_em_sim_mediumid;   //!
   TBranch        *b_em_sim_p;   //!
   TBranch        *b_em_sim_parentid;   //!
   TBranch        *b_em_sim_processid;   //!
   TBranch        *b_em_sim_px;   //!
   TBranch        *b_em_sim_py;   //!
   TBranch        *b_em_sim_pz;   //!
   TBranch        *b_em_sim_vertexx;   //!
   TBranch        *b_em_sim_vertexy;   //!
   TBranch        *b_em_sim_vertexz;   //!
   TBranch        *b_em_system;   //!
   TBranch        *b_em_theta;   //!
   TBranch        *b_em_theta_rich;   //!
   TBranch        *b_em_tof_mom;   //!
   TBranch        *b_em_tof_new;   //!
   TBranch        *b_em_tof_rec;   //!
   TBranch        *b_em_track_length;   //!
   TBranch        *b_em_tracklength;   //!
   TBranch        *b_em_z;   //!
   TBranch        *b_ep_beta;   //!
   TBranch        *b_ep_beta_new;   //!
   TBranch        *b_ep_btChargeRing;   //!
   TBranch        *b_ep_btChargeSum;   //!
   TBranch        *b_ep_btChi2;   //!
   TBranch        *b_ep_btClusters;   //!
   TBranch        *b_ep_btMaxima;   //!
   TBranch        *b_ep_btMaximaCharge;   //!
   TBranch        *b_ep_btMaximaChargeShared;   //!
   TBranch        *b_ep_btMaximaChargeSharedFragment;   //!
   TBranch        *b_ep_btMaximaShared;   //!
   TBranch        *b_ep_btMaximaSharedFragment;   //!
   TBranch        *b_ep_btMeanDist;   //!
   TBranch        *b_ep_btNearbyMaxima;   //!
   TBranch        *b_ep_btNearbyMaximaShared;   //!
   TBranch        *b_ep_btPadsClus;   //!
   TBranch        *b_ep_btPadsRing;   //!
   TBranch        *b_ep_btRingMatrix;   //!
   TBranch        *b_ep_dedx_mdc;   //!
   TBranch        *b_ep_dedx_tof;   //!
   TBranch        *b_ep_id;   //!
   TBranch        *b_ep_isBT;   //!
   TBranch        *b_ep_isOffVertexClust;   //!
   TBranch        *b_ep_isPrimaryVertex;   //!
   TBranch        *b_ep_isUsedVertex;   //!
   TBranch        *b_ep_isring;   //!
   TBranch        *b_ep_isringmdc;   //!
   TBranch        *b_ep_isringnomatch;   //!
   TBranch        *b_ep_isringtrack;   //!
   TBranch        *b_ep_kIsLepton;   //!
   TBranch        *b_ep_kIsUsed;   //!
   TBranch        *b_ep_mdcinnerchi2;   //!
   TBranch        *b_ep_mdcouterchi2;   //!
   TBranch        *b_ep_oa_hadr;   //!
   TBranch        *b_ep_oa_lept;   //!
   TBranch        *b_ep_p;   //!
   TBranch        *b_ep_p_corr_em;   //!
   TBranch        *b_ep_p_corr_ep;   //!
   TBranch        *b_ep_p_corr_p;   //!
   TBranch        *b_ep_p_corr_pim;   //!
   TBranch        *b_ep_p_corr_pip;   //!
   TBranch        *b_ep_phi;   //!
   TBranch        *b_ep_phi_rich;   //!
   TBranch        *b_ep_pid;   //!
   TBranch        *b_ep_q;   //!
   TBranch        *b_ep_r;   //!
   TBranch        *b_ep_resolution;   //!
   TBranch        *b_ep_resoultion;   //!
   TBranch        *b_ep_rich_amp;   //!
   TBranch        *b_ep_rich_centr;   //!
   TBranch        *b_ep_rich_houtra;   //!
   TBranch        *b_ep_rich_padnum;   //!
   TBranch        *b_ep_rich_patmat;   //!
   TBranch        *b_ep_rkchi2;   //!
   TBranch        *b_ep_sector;   //!
   TBranch        *b_ep_shw_sum0;   //!
   TBranch        *b_ep_shw_sum1;   //!
   TBranch        *b_ep_shw_sum2;   //!
   TBranch        *b_ep_sim_geninfo;   //!
   TBranch        *b_ep_sim_geninfo1;   //!
   TBranch        *b_ep_sim_geninfo2;   //!
   TBranch        *b_ep_sim_genweight;   //!
   TBranch        *b_ep_sim_grandparentid;   //!
   TBranch        *b_ep_sim_id;   //!
   TBranch        *b_ep_sim_iscommon;   //!
   TBranch        *b_ep_sim_isghost;   //!
   TBranch        *b_ep_sim_mediumid;   //!
   TBranch        *b_ep_sim_p;   //!
   TBranch        *b_ep_sim_parentid;   //!
   TBranch        *b_ep_sim_processid;   //!
   TBranch        *b_ep_sim_px;   //!
   TBranch        *b_ep_sim_py;   //!
   TBranch        *b_ep_sim_pz;   //!
   TBranch        *b_ep_sim_vertexx;   //!
   TBranch        *b_ep_sim_vertexy;   //!
   TBranch        *b_ep_sim_vertexz;   //!
   TBranch        *b_ep_system;   //!
   TBranch        *b_ep_theta;   //!
   TBranch        *b_ep_theta_rich;   //!
   TBranch        *b_ep_tof_mom;   //!
   TBranch        *b_ep_tof_new;   //!
   TBranch        *b_ep_tof_rec;   //!
   TBranch        *b_ep_track_length;   //!
   TBranch        *b_ep_tracklength;   //!
   TBranch        *b_ep_z;   //!
   TBranch        *b_event;   //!
   TBranch        *b_evtFlashMDC;   //!
   TBranch        *b_evtGood;   //!
   TBranch        *b_evtGoodTrig;   //!
   TBranch        *b_evtLepMult;   //!
   TBranch        *b_evtMDCMult;   //!
   TBranch        *b_evtPileupMDC;   //!
   TBranch        *b_evtPileupMeta;   //!
   TBranch        *b_evtPileupStart;   //!
   TBranch        *b_evtStart;   //!
   TBranch        *b_evtVertCand;   //!
   TBranch        *b_evtVertClust;   //!
   TBranch        *b_fw_beta_1;   //!
   TBranch        *b_fw_beta_2;   //!
   TBranch        *b_fw_beta_3;   //!
   TBranch        *b_fw_charge_1;   //!
   TBranch        *b_fw_charge_2;   //!
   TBranch        *b_fw_charge_3;   //!
   TBranch        *b_fw_cluster_mult;   //!
   TBranch        *b_fw_distance_1;   //!
   TBranch        *b_fw_distance_2;   //!
   TBranch        *b_fw_distance_3;   //!
   TBranch        *b_fw_mult;   //!
   TBranch        *b_fw_p_1;   //!
   TBranch        *b_fw_p_2;   //!
   TBranch        *b_fw_p_3;   //!
   TBranch        *b_fw_phi_1;   //!
   TBranch        *b_fw_phi_2;   //!
   TBranch        *b_fw_phi_3;   //!
   TBranch        *b_fw_sim_geninfo1_1;   //!
   TBranch        *b_fw_sim_geninfo1_2;   //!
   TBranch        *b_fw_sim_geninfo1_3;   //!
   TBranch        *b_fw_sim_geninfo2_1;   //!
   TBranch        *b_fw_sim_geninfo2_2;   //!
   TBranch        *b_fw_sim_geninfo2_3;   //!
   TBranch        *b_fw_sim_geninfo_1;   //!
   TBranch        *b_fw_sim_geninfo_2;   //!
   TBranch        *b_fw_sim_geninfo_3;   //!
   TBranch        *b_fw_sim_genweight_1;   //!
   TBranch        *b_fw_sim_genweight_2;   //!
   TBranch        *b_fw_sim_genweight_3;   //!
   TBranch        *b_fw_sim_id_1;   //!
   TBranch        *b_fw_sim_id_2;   //!
   TBranch        *b_fw_sim_id_3;   //!
   TBranch        *b_fw_sim_p_1;   //!
   TBranch        *b_fw_sim_p_2;   //!
   TBranch        *b_fw_sim_p_3;   //!
   TBranch        *b_fw_sim_parentid_1;   //!
   TBranch        *b_fw_sim_parentid_2;   //!
   TBranch        *b_fw_sim_parentid_3;   //!
   TBranch        *b_fw_size_1;   //!
   TBranch        *b_fw_size_2;   //!
   TBranch        *b_fw_size_3;   //!
   TBranch        *b_fw_spectator_1;   //!
   TBranch        *b_fw_spectator_2;   //!
   TBranch        *b_fw_spectator_3;   //!
   TBranch        *b_fw_theta_1;   //!
   TBranch        *b_fw_theta_2;   //!
   TBranch        *b_fw_theta_3;   //!
   TBranch        *b_fw_time_1;   //!
   TBranch        *b_fw_time_2;   //!
   TBranch        *b_fw_time_3;   //!
   TBranch        *b_fw_time_min_1;   //!
   TBranch        *b_fw_time_min_2;   //!
   TBranch        *b_fw_time_min_3;   //!
   TBranch        *b_fw_tracknr_1;   //!
   TBranch        *b_fw_tracknr_12;   //!
   TBranch        *b_fw_tracknr_2;   //!
   TBranch        *b_fw_tracknr_22;   //!
   TBranch        *b_fw_tracknr_3;   //!
   TBranch        *b_fw_tracknr_32;   //!
   TBranch        *b_fw_x_lab_1;   //!
   TBranch        *b_fw_x_lab_2;   //!
   TBranch        *b_fw_x_lab_3;   //!
   TBranch        *b_fw_y_lab_1;   //!
   TBranch        *b_fw_y_lab_2;   //!
   TBranch        *b_fw_y_lab_3;   //!
   TBranch        *b_fw_z_lab_1;   //!
   TBranch        *b_fw_z_lab_2;   //!
   TBranch        *b_fw_z_lab_3;   //!
   TBranch        *b_geant_sim_genweight_1;   //!
   TBranch        *b_geant_sim_genweight_2;   //!
   TBranch        *b_geant_sim_id_1;   //!
   TBranch        *b_geant_sim_id_2;   //!
   TBranch        *b_geant_sim_p_1;   //!
   TBranch        *b_geant_sim_p_2;   //!
   TBranch        *b_geant_sim_px_1;   //!
   TBranch        *b_geant_sim_px_2;   //!
   TBranch        *b_geant_sim_py_1;   //!
   TBranch        *b_geant_sim_py_2;   //!
   TBranch        *b_geant_sim_pz_1;   //!
   TBranch        *b_geant_sim_pz_2;   //!
   TBranch        *b_hneg_mult;   //!
   TBranch        *b_hpos_mult;   //!
   TBranch        *b_isBest;   //!
   TBranch        *b_lneg_mult;   //!
   TBranch        *b_lpos_mult;   //!
   TBranch        *b_pim_beta;   //!
   TBranch        *b_pim_beta_new;   //!
   TBranch        *b_pim_btChargeRing;   //!
   TBranch        *b_pim_btChargeSum;   //!
   TBranch        *b_pim_btChi2;   //!
   TBranch        *b_pim_btClusters;   //!
   TBranch        *b_pim_btMaxima;   //!
   TBranch        *b_pim_btMaximaCharge;   //!
   TBranch        *b_pim_btMaximaChargeShared;   //!
   TBranch        *b_pim_btMaximaChargeSharedFragment;   //!
   TBranch        *b_pim_btMaximaShared;   //!
   TBranch        *b_pim_btMaximaSharedFragment;   //!
   TBranch        *b_pim_btMeanDist;   //!
   TBranch        *b_pim_btNearbyMaxima;   //!
   TBranch        *b_pim_btNearbyMaximaShared;   //!
   TBranch        *b_pim_btPadsClus;   //!
   TBranch        *b_pim_btPadsRing;   //!
   TBranch        *b_pim_btRingMatrix;   //!
   TBranch        *b_pim_dedx_mdc;   //!
   TBranch        *b_pim_dedx_tof;   //!
   TBranch        *b_pim_id;   //!
   TBranch        *b_pim_isBT;   //!
   TBranch        *b_pim_isOffVertexClust;   //!
   TBranch        *b_pim_isPrimaryVertex;   //!
   TBranch        *b_pim_isUsedVertex;   //!
   TBranch        *b_pim_isring;   //!
   TBranch        *b_pim_isringmdc;   //!
   TBranch        *b_pim_isringnomatch;   //!
   TBranch        *b_pim_isringtrack;   //!
   TBranch        *b_pim_kIsLepton;   //!
   TBranch        *b_pim_kIsUsed;   //!
   TBranch        *b_pim_mdcinnerchi2;   //!
   TBranch        *b_pim_mdcouterchi2;   //!
   TBranch        *b_pim_oa_hadr;   //!
   TBranch        *b_pim_oa_lept;   //!
   TBranch        *b_pim_p;   //!
   TBranch        *b_pim_p_corr_em;   //!
   TBranch        *b_pim_p_corr_ep;   //!
   TBranch        *b_pim_p_corr_p;   //!
   TBranch        *b_pim_p_corr_pim;   //!
   TBranch        *b_pim_p_corr_pip;   //!
   TBranch        *b_pim_phi;   //!
   TBranch        *b_pim_phi_rich;   //!
   TBranch        *b_pim_pid;   //!
   TBranch        *b_pim_q;   //!
   TBranch        *b_pim_r;   //!
   TBranch        *b_pim_resolution;   //!
   TBranch        *b_pim_resoultion;   //!
   TBranch        *b_pim_rich_amp;   //!
   TBranch        *b_pim_rich_centr;   //!
   TBranch        *b_pim_rich_houtra;   //!
   TBranch        *b_pim_rich_padnum;   //!
   TBranch        *b_pim_rich_patmat;   //!
   TBranch        *b_pim_rkchi2;   //!
   TBranch        *b_pim_sector;   //!
   TBranch        *b_pim_shw_sum0;   //!
   TBranch        *b_pim_shw_sum1;   //!
   TBranch        *b_pim_shw_sum2;   //!
   TBranch        *b_pim_sim_geninfo;   //!
   TBranch        *b_pim_sim_geninfo1;   //!
   TBranch        *b_pim_sim_geninfo2;   //!
   TBranch        *b_pim_sim_genweight;   //!
   TBranch        *b_pim_sim_grandparentid;   //!
   TBranch        *b_pim_sim_id;   //!
   TBranch        *b_pim_sim_iscommon;   //!
   TBranch        *b_pim_sim_isghost;   //!
   TBranch        *b_pim_sim_mediumid;   //!
   TBranch        *b_pim_sim_p;   //!
   TBranch        *b_pim_sim_parentid;   //!
   TBranch        *b_pim_sim_processid;   //!
   TBranch        *b_pim_sim_px;   //!
   TBranch        *b_pim_sim_py;   //!
   TBranch        *b_pim_sim_pz;   //!
   TBranch        *b_pim_sim_vertexx;   //!
   TBranch        *b_pim_sim_vertexy;   //!
   TBranch        *b_pim_sim_vertexz;   //!
   TBranch        *b_pim_system;   //!
   TBranch        *b_pim_theta;   //!
   TBranch        *b_pim_theta_rich;   //!
   TBranch        *b_pim_tof_mom;   //!
   TBranch        *b_pim_tof_new;   //!
   TBranch        *b_pim_tof_rec;   //!
   TBranch        *b_pim_track_length;   //!
   TBranch        *b_pim_tracklength;   //!
   TBranch        *b_pim_z;   //!
   TBranch        *b_runnumber;   //!
   TBranch        *b_totalmult;   //!
   TBranch        *b_trigbit;   //!
   TBranch        *b_trigdec;   //!
   TBranch        *b_trigdownscale;   //!
   TBranch        *b_trigdownscaleflag;   //!

   PimEpEm_ID(TTree *tree=0);
   virtual ~PimEpEm_ID();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef PimEpEm_ID_cxx
PimEpEm_ID::PimEpEm_ID(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("hgeantout_PWA_momentums-ppimpi0_ev_weight_ALL.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("hgeantout_PWA_momentums-ppimpi0_ev_weight_ALL.root");
      }
      f->GetObject("PimEpEm_ID",tree);

   }
   Init(tree);
}

PimEpEm_ID::~PimEpEm_ID()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t PimEpEm_ID::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t PimEpEm_ID::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void PimEpEm_ID::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("eVertClust_chi2", &eVertClust_chi2, &b_eVertClust_chi2);
   fChain->SetBranchAddress("eVertClust_x", &eVertClust_x, &b_eVertClust_x);
   fChain->SetBranchAddress("eVertClust_y", &eVertClust_y, &b_eVertClust_y);
   fChain->SetBranchAddress("eVertClust_z", &eVertClust_z, &b_eVertClust_z);
   fChain->SetBranchAddress("eVertReco_chi2", &eVertReco_chi2, &b_eVertReco_chi2);
   fChain->SetBranchAddress("eVertReco_x", &eVertReco_x, &b_eVertReco_x);
   fChain->SetBranchAddress("eVertReco_y", &eVertReco_y, &b_eVertReco_y);
   fChain->SetBranchAddress("eVertReco_z", &eVertReco_z, &b_eVertReco_z);
   fChain->SetBranchAddress("eVert_chi2", &eVert_chi2, &b_eVert_chi2);
   fChain->SetBranchAddress("eVert_x", &eVert_x, &b_eVert_x);
   fChain->SetBranchAddress("eVert_y", &eVert_y, &b_eVert_y);
   fChain->SetBranchAddress("eVert_z", &eVert_z, &b_eVert_z);
   fChain->SetBranchAddress("em_beta", &em_beta, &b_em_beta);
   fChain->SetBranchAddress("em_beta_new", &em_beta_new, &b_em_beta_new);
   fChain->SetBranchAddress("em_btChargeRing", &em_btChargeRing, &b_em_btChargeRing);
   fChain->SetBranchAddress("em_btChargeSum", &em_btChargeSum, &b_em_btChargeSum);
   fChain->SetBranchAddress("em_btChi2", &em_btChi2, &b_em_btChi2);
   fChain->SetBranchAddress("em_btClusters", &em_btClusters, &b_em_btClusters);
   fChain->SetBranchAddress("em_btMaxima", &em_btMaxima, &b_em_btMaxima);
   fChain->SetBranchAddress("em_btMaximaCharge", &em_btMaximaCharge, &b_em_btMaximaCharge);
   fChain->SetBranchAddress("em_btMaximaChargeShared", &em_btMaximaChargeShared, &b_em_btMaximaChargeShared);
   fChain->SetBranchAddress("em_btMaximaChargeSharedFragment", &em_btMaximaChargeSharedFragment, &b_em_btMaximaChargeSharedFragment);
   fChain->SetBranchAddress("em_btMaximaShared", &em_btMaximaShared, &b_em_btMaximaShared);
   fChain->SetBranchAddress("em_btMaximaSharedFragment", &em_btMaximaSharedFragment, &b_em_btMaximaSharedFragment);
   fChain->SetBranchAddress("em_btMeanDist", &em_btMeanDist, &b_em_btMeanDist);
   fChain->SetBranchAddress("em_btNearbyMaxima", &em_btNearbyMaxima, &b_em_btNearbyMaxima);
   fChain->SetBranchAddress("em_btNearbyMaximaShared", &em_btNearbyMaximaShared, &b_em_btNearbyMaximaShared);
   fChain->SetBranchAddress("em_btPadsClus", &em_btPadsClus, &b_em_btPadsClus);
   fChain->SetBranchAddress("em_btPadsRing", &em_btPadsRing, &b_em_btPadsRing);
   fChain->SetBranchAddress("em_btRingMatrix", &em_btRingMatrix, &b_em_btRingMatrix);
   fChain->SetBranchAddress("em_dedx_mdc", &em_dedx_mdc, &b_em_dedx_mdc);
   fChain->SetBranchAddress("em_dedx_tof", &em_dedx_tof, &b_em_dedx_tof);
   fChain->SetBranchAddress("em_id", &em_id, &b_em_id);
   fChain->SetBranchAddress("em_isBT", &em_isBT, &b_em_isBT);
   fChain->SetBranchAddress("em_isOffVertexClust", &em_isOffVertexClust, &b_em_isOffVertexClust);
   fChain->SetBranchAddress("em_isPrimaryVertex", &em_isPrimaryVertex, &b_em_isPrimaryVertex);
   fChain->SetBranchAddress("em_isUsedVertex", &em_isUsedVertex, &b_em_isUsedVertex);
   fChain->SetBranchAddress("em_isring", &em_isring, &b_em_isring);
   fChain->SetBranchAddress("em_isringmdc", &em_isringmdc, &b_em_isringmdc);
   fChain->SetBranchAddress("em_isringnomatch", &em_isringnomatch, &b_em_isringnomatch);
   fChain->SetBranchAddress("em_isringtrack", &em_isringtrack, &b_em_isringtrack);
   fChain->SetBranchAddress("em_kIsLepton", &em_kIsLepton, &b_em_kIsLepton);
   fChain->SetBranchAddress("em_kIsUsed", &em_kIsUsed, &b_em_kIsUsed);
   fChain->SetBranchAddress("em_mdcinnerchi2", &em_mdcinnerchi2, &b_em_mdcinnerchi2);
   fChain->SetBranchAddress("em_mdcouterchi2", &em_mdcouterchi2, &b_em_mdcouterchi2);
   fChain->SetBranchAddress("em_oa_hadr", &em_oa_hadr, &b_em_oa_hadr);
   fChain->SetBranchAddress("em_oa_lept", &em_oa_lept, &b_em_oa_lept);
   fChain->SetBranchAddress("em_p", &em_p, &b_em_p);
   fChain->SetBranchAddress("em_p_corr_em", &em_p_corr_em, &b_em_p_corr_em);
   fChain->SetBranchAddress("em_p_corr_ep", &em_p_corr_ep, &b_em_p_corr_ep);
   fChain->SetBranchAddress("em_p_corr_p", &em_p_corr_p, &b_em_p_corr_p);
   fChain->SetBranchAddress("em_p_corr_pim", &em_p_corr_pim, &b_em_p_corr_pim);
   fChain->SetBranchAddress("em_p_corr_pip", &em_p_corr_pip, &b_em_p_corr_pip);
   fChain->SetBranchAddress("em_phi", &em_phi, &b_em_phi);
   fChain->SetBranchAddress("em_phi_rich", &em_phi_rich, &b_em_phi_rich);
   fChain->SetBranchAddress("em_pid", &em_pid, &b_em_pid);
   fChain->SetBranchAddress("em_q", &em_q, &b_em_q);
   fChain->SetBranchAddress("em_r", &em_r, &b_em_r);
   fChain->SetBranchAddress("em_resolution", &em_resolution, &b_em_resolution);
   fChain->SetBranchAddress("em_resoultion", &em_resoultion, &b_em_resoultion);
   fChain->SetBranchAddress("em_rich_amp", &em_rich_amp, &b_em_rich_amp);
   fChain->SetBranchAddress("em_rich_centr", &em_rich_centr, &b_em_rich_centr);
   fChain->SetBranchAddress("em_rich_houtra", &em_rich_houtra, &b_em_rich_houtra);
   fChain->SetBranchAddress("em_rich_padnum", &em_rich_padnum, &b_em_rich_padnum);
   fChain->SetBranchAddress("em_rich_patmat", &em_rich_patmat, &b_em_rich_patmat);
   fChain->SetBranchAddress("em_rkchi2", &em_rkchi2, &b_em_rkchi2);
   fChain->SetBranchAddress("em_sector", &em_sector, &b_em_sector);
   fChain->SetBranchAddress("em_shw_sum0", &em_shw_sum0, &b_em_shw_sum0);
   fChain->SetBranchAddress("em_shw_sum1", &em_shw_sum1, &b_em_shw_sum1);
   fChain->SetBranchAddress("em_shw_sum2", &em_shw_sum2, &b_em_shw_sum2);
   fChain->SetBranchAddress("em_sim_geninfo", &em_sim_geninfo, &b_em_sim_geninfo);
   fChain->SetBranchAddress("em_sim_geninfo1", &em_sim_geninfo1, &b_em_sim_geninfo1);
   fChain->SetBranchAddress("em_sim_geninfo2", &em_sim_geninfo2, &b_em_sim_geninfo2);
   fChain->SetBranchAddress("em_sim_genweight", &em_sim_genweight, &b_em_sim_genweight);
   fChain->SetBranchAddress("em_sim_grandparentid", &em_sim_grandparentid, &b_em_sim_grandparentid);
   fChain->SetBranchAddress("em_sim_id", &em_sim_id, &b_em_sim_id);
   fChain->SetBranchAddress("em_sim_iscommon", &em_sim_iscommon, &b_em_sim_iscommon);
   fChain->SetBranchAddress("em_sim_isghost", &em_sim_isghost, &b_em_sim_isghost);
   fChain->SetBranchAddress("em_sim_mediumid", &em_sim_mediumid, &b_em_sim_mediumid);
   fChain->SetBranchAddress("em_sim_p", &em_sim_p, &b_em_sim_p);
   fChain->SetBranchAddress("em_sim_parentid", &em_sim_parentid, &b_em_sim_parentid);
   fChain->SetBranchAddress("em_sim_processid", &em_sim_processid, &b_em_sim_processid);
   fChain->SetBranchAddress("em_sim_px", &em_sim_px, &b_em_sim_px);
   fChain->SetBranchAddress("em_sim_py", &em_sim_py, &b_em_sim_py);
   fChain->SetBranchAddress("em_sim_pz", &em_sim_pz, &b_em_sim_pz);
   fChain->SetBranchAddress("em_sim_vertexx", &em_sim_vertexx, &b_em_sim_vertexx);
   fChain->SetBranchAddress("em_sim_vertexy", &em_sim_vertexy, &b_em_sim_vertexy);
   fChain->SetBranchAddress("em_sim_vertexz", &em_sim_vertexz, &b_em_sim_vertexz);
   fChain->SetBranchAddress("em_system", &em_system, &b_em_system);
   fChain->SetBranchAddress("em_theta", &em_theta, &b_em_theta);
   fChain->SetBranchAddress("em_theta_rich", &em_theta_rich, &b_em_theta_rich);
   fChain->SetBranchAddress("em_tof_mom", &em_tof_mom, &b_em_tof_mom);
   fChain->SetBranchAddress("em_tof_new", &em_tof_new, &b_em_tof_new);
   fChain->SetBranchAddress("em_tof_rec", &em_tof_rec, &b_em_tof_rec);
   fChain->SetBranchAddress("em_track_length", &em_track_length, &b_em_track_length);
   fChain->SetBranchAddress("em_tracklength", &em_tracklength, &b_em_tracklength);
   fChain->SetBranchAddress("em_z", &em_z, &b_em_z);
   fChain->SetBranchAddress("ep_beta", &ep_beta, &b_ep_beta);
   fChain->SetBranchAddress("ep_beta_new", &ep_beta_new, &b_ep_beta_new);
   fChain->SetBranchAddress("ep_btChargeRing", &ep_btChargeRing, &b_ep_btChargeRing);
   fChain->SetBranchAddress("ep_btChargeSum", &ep_btChargeSum, &b_ep_btChargeSum);
   fChain->SetBranchAddress("ep_btChi2", &ep_btChi2, &b_ep_btChi2);
   fChain->SetBranchAddress("ep_btClusters", &ep_btClusters, &b_ep_btClusters);
   fChain->SetBranchAddress("ep_btMaxima", &ep_btMaxima, &b_ep_btMaxima);
   fChain->SetBranchAddress("ep_btMaximaCharge", &ep_btMaximaCharge, &b_ep_btMaximaCharge);
   fChain->SetBranchAddress("ep_btMaximaChargeShared", &ep_btMaximaChargeShared, &b_ep_btMaximaChargeShared);
   fChain->SetBranchAddress("ep_btMaximaChargeSharedFragment", &ep_btMaximaChargeSharedFragment, &b_ep_btMaximaChargeSharedFragment);
   fChain->SetBranchAddress("ep_btMaximaShared", &ep_btMaximaShared, &b_ep_btMaximaShared);
   fChain->SetBranchAddress("ep_btMaximaSharedFragment", &ep_btMaximaSharedFragment, &b_ep_btMaximaSharedFragment);
   fChain->SetBranchAddress("ep_btMeanDist", &ep_btMeanDist, &b_ep_btMeanDist);
   fChain->SetBranchAddress("ep_btNearbyMaxima", &ep_btNearbyMaxima, &b_ep_btNearbyMaxima);
   fChain->SetBranchAddress("ep_btNearbyMaximaShared", &ep_btNearbyMaximaShared, &b_ep_btNearbyMaximaShared);
   fChain->SetBranchAddress("ep_btPadsClus", &ep_btPadsClus, &b_ep_btPadsClus);
   fChain->SetBranchAddress("ep_btPadsRing", &ep_btPadsRing, &b_ep_btPadsRing);
   fChain->SetBranchAddress("ep_btRingMatrix", &ep_btRingMatrix, &b_ep_btRingMatrix);
   fChain->SetBranchAddress("ep_dedx_mdc", &ep_dedx_mdc, &b_ep_dedx_mdc);
   fChain->SetBranchAddress("ep_dedx_tof", &ep_dedx_tof, &b_ep_dedx_tof);
   fChain->SetBranchAddress("ep_id", &ep_id, &b_ep_id);
   fChain->SetBranchAddress("ep_isBT", &ep_isBT, &b_ep_isBT);
   fChain->SetBranchAddress("ep_isOffVertexClust", &ep_isOffVertexClust, &b_ep_isOffVertexClust);
   fChain->SetBranchAddress("ep_isPrimaryVertex", &ep_isPrimaryVertex, &b_ep_isPrimaryVertex);
   fChain->SetBranchAddress("ep_isUsedVertex", &ep_isUsedVertex, &b_ep_isUsedVertex);
   fChain->SetBranchAddress("ep_isring", &ep_isring, &b_ep_isring);
   fChain->SetBranchAddress("ep_isringmdc", &ep_isringmdc, &b_ep_isringmdc);
   fChain->SetBranchAddress("ep_isringnomatch", &ep_isringnomatch, &b_ep_isringnomatch);
   fChain->SetBranchAddress("ep_isringtrack", &ep_isringtrack, &b_ep_isringtrack);
   fChain->SetBranchAddress("ep_kIsLepton", &ep_kIsLepton, &b_ep_kIsLepton);
   fChain->SetBranchAddress("ep_kIsUsed", &ep_kIsUsed, &b_ep_kIsUsed);
   fChain->SetBranchAddress("ep_mdcinnerchi2", &ep_mdcinnerchi2, &b_ep_mdcinnerchi2);
   fChain->SetBranchAddress("ep_mdcouterchi2", &ep_mdcouterchi2, &b_ep_mdcouterchi2);
   fChain->SetBranchAddress("ep_oa_hadr", &ep_oa_hadr, &b_ep_oa_hadr);
   fChain->SetBranchAddress("ep_oa_lept", &ep_oa_lept, &b_ep_oa_lept);
   fChain->SetBranchAddress("ep_p", &ep_p, &b_ep_p);
   fChain->SetBranchAddress("ep_p_corr_em", &ep_p_corr_em, &b_ep_p_corr_em);
   fChain->SetBranchAddress("ep_p_corr_ep", &ep_p_corr_ep, &b_ep_p_corr_ep);
   fChain->SetBranchAddress("ep_p_corr_p", &ep_p_corr_p, &b_ep_p_corr_p);
   fChain->SetBranchAddress("ep_p_corr_pim", &ep_p_corr_pim, &b_ep_p_corr_pim);
   fChain->SetBranchAddress("ep_p_corr_pip", &ep_p_corr_pip, &b_ep_p_corr_pip);
   fChain->SetBranchAddress("ep_phi", &ep_phi, &b_ep_phi);
   fChain->SetBranchAddress("ep_phi_rich", &ep_phi_rich, &b_ep_phi_rich);
   fChain->SetBranchAddress("ep_pid", &ep_pid, &b_ep_pid);
   fChain->SetBranchAddress("ep_q", &ep_q, &b_ep_q);
   fChain->SetBranchAddress("ep_r", &ep_r, &b_ep_r);
   fChain->SetBranchAddress("ep_resolution", &ep_resolution, &b_ep_resolution);
   fChain->SetBranchAddress("ep_resoultion", &ep_resoultion, &b_ep_resoultion);
   fChain->SetBranchAddress("ep_rich_amp", &ep_rich_amp, &b_ep_rich_amp);
   fChain->SetBranchAddress("ep_rich_centr", &ep_rich_centr, &b_ep_rich_centr);
   fChain->SetBranchAddress("ep_rich_houtra", &ep_rich_houtra, &b_ep_rich_houtra);
   fChain->SetBranchAddress("ep_rich_padnum", &ep_rich_padnum, &b_ep_rich_padnum);
   fChain->SetBranchAddress("ep_rich_patmat", &ep_rich_patmat, &b_ep_rich_patmat);
   fChain->SetBranchAddress("ep_rkchi2", &ep_rkchi2, &b_ep_rkchi2);
   fChain->SetBranchAddress("ep_sector", &ep_sector, &b_ep_sector);
   fChain->SetBranchAddress("ep_shw_sum0", &ep_shw_sum0, &b_ep_shw_sum0);
   fChain->SetBranchAddress("ep_shw_sum1", &ep_shw_sum1, &b_ep_shw_sum1);
   fChain->SetBranchAddress("ep_shw_sum2", &ep_shw_sum2, &b_ep_shw_sum2);
   fChain->SetBranchAddress("ep_sim_geninfo", &ep_sim_geninfo, &b_ep_sim_geninfo);
   fChain->SetBranchAddress("ep_sim_geninfo1", &ep_sim_geninfo1, &b_ep_sim_geninfo1);
   fChain->SetBranchAddress("ep_sim_geninfo2", &ep_sim_geninfo2, &b_ep_sim_geninfo2);
   fChain->SetBranchAddress("ep_sim_genweight", &ep_sim_genweight, &b_ep_sim_genweight);
   fChain->SetBranchAddress("ep_sim_grandparentid", &ep_sim_grandparentid, &b_ep_sim_grandparentid);
   fChain->SetBranchAddress("ep_sim_id", &ep_sim_id, &b_ep_sim_id);
   fChain->SetBranchAddress("ep_sim_iscommon", &ep_sim_iscommon, &b_ep_sim_iscommon);
   fChain->SetBranchAddress("ep_sim_isghost", &ep_sim_isghost, &b_ep_sim_isghost);
   fChain->SetBranchAddress("ep_sim_mediumid", &ep_sim_mediumid, &b_ep_sim_mediumid);
   fChain->SetBranchAddress("ep_sim_p", &ep_sim_p, &b_ep_sim_p);
   fChain->SetBranchAddress("ep_sim_parentid", &ep_sim_parentid, &b_ep_sim_parentid);
   fChain->SetBranchAddress("ep_sim_processid", &ep_sim_processid, &b_ep_sim_processid);
   fChain->SetBranchAddress("ep_sim_px", &ep_sim_px, &b_ep_sim_px);
   fChain->SetBranchAddress("ep_sim_py", &ep_sim_py, &b_ep_sim_py);
   fChain->SetBranchAddress("ep_sim_pz", &ep_sim_pz, &b_ep_sim_pz);
   fChain->SetBranchAddress("ep_sim_vertexx", &ep_sim_vertexx, &b_ep_sim_vertexx);
   fChain->SetBranchAddress("ep_sim_vertexy", &ep_sim_vertexy, &b_ep_sim_vertexy);
   fChain->SetBranchAddress("ep_sim_vertexz", &ep_sim_vertexz, &b_ep_sim_vertexz);
   fChain->SetBranchAddress("ep_system", &ep_system, &b_ep_system);
   fChain->SetBranchAddress("ep_theta", &ep_theta, &b_ep_theta);
   fChain->SetBranchAddress("ep_theta_rich", &ep_theta_rich, &b_ep_theta_rich);
   fChain->SetBranchAddress("ep_tof_mom", &ep_tof_mom, &b_ep_tof_mom);
   fChain->SetBranchAddress("ep_tof_new", &ep_tof_new, &b_ep_tof_new);
   fChain->SetBranchAddress("ep_tof_rec", &ep_tof_rec, &b_ep_tof_rec);
   fChain->SetBranchAddress("ep_track_length", &ep_track_length, &b_ep_track_length);
   fChain->SetBranchAddress("ep_tracklength", &ep_tracklength, &b_ep_tracklength);
   fChain->SetBranchAddress("ep_z", &ep_z, &b_ep_z);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("evtFlashMDC", &evtFlashMDC, &b_evtFlashMDC);
   fChain->SetBranchAddress("evtGood", &evtGood, &b_evtGood);
   fChain->SetBranchAddress("evtGoodTrig", &evtGoodTrig, &b_evtGoodTrig);
   fChain->SetBranchAddress("evtLepMult", &evtLepMult, &b_evtLepMult);
   fChain->SetBranchAddress("evtMDCMult", &evtMDCMult, &b_evtMDCMult);
   fChain->SetBranchAddress("evtPileupMDC", &evtPileupMDC, &b_evtPileupMDC);
   fChain->SetBranchAddress("evtPileupMeta", &evtPileupMeta, &b_evtPileupMeta);
   fChain->SetBranchAddress("evtPileupStart", &evtPileupStart, &b_evtPileupStart);
   fChain->SetBranchAddress("evtStart", &evtStart, &b_evtStart);
   fChain->SetBranchAddress("evtVertCand", &evtVertCand, &b_evtVertCand);
   fChain->SetBranchAddress("evtVertClust", &evtVertClust, &b_evtVertClust);
   fChain->SetBranchAddress("fw_beta_1", &fw_beta_1, &b_fw_beta_1);
   fChain->SetBranchAddress("fw_beta_2", &fw_beta_2, &b_fw_beta_2);
   fChain->SetBranchAddress("fw_beta_3", &fw_beta_3, &b_fw_beta_3);
   fChain->SetBranchAddress("fw_charge_1", &fw_charge_1, &b_fw_charge_1);
   fChain->SetBranchAddress("fw_charge_2", &fw_charge_2, &b_fw_charge_2);
   fChain->SetBranchAddress("fw_charge_3", &fw_charge_3, &b_fw_charge_3);
   fChain->SetBranchAddress("fw_cluster_mult", &fw_cluster_mult, &b_fw_cluster_mult);
   fChain->SetBranchAddress("fw_distance_1", &fw_distance_1, &b_fw_distance_1);
   fChain->SetBranchAddress("fw_distance_2", &fw_distance_2, &b_fw_distance_2);
   fChain->SetBranchAddress("fw_distance_3", &fw_distance_3, &b_fw_distance_3);
   fChain->SetBranchAddress("fw_mult", &fw_mult, &b_fw_mult);
   fChain->SetBranchAddress("fw_p_1", &fw_p_1, &b_fw_p_1);
   fChain->SetBranchAddress("fw_p_2", &fw_p_2, &b_fw_p_2);
   fChain->SetBranchAddress("fw_p_3", &fw_p_3, &b_fw_p_3);
   fChain->SetBranchAddress("fw_phi_1", &fw_phi_1, &b_fw_phi_1);
   fChain->SetBranchAddress("fw_phi_2", &fw_phi_2, &b_fw_phi_2);
   fChain->SetBranchAddress("fw_phi_3", &fw_phi_3, &b_fw_phi_3);
   fChain->SetBranchAddress("fw_sim_geninfo1_1", &fw_sim_geninfo1_1, &b_fw_sim_geninfo1_1);
   fChain->SetBranchAddress("fw_sim_geninfo1_2", &fw_sim_geninfo1_2, &b_fw_sim_geninfo1_2);
   fChain->SetBranchAddress("fw_sim_geninfo1_3", &fw_sim_geninfo1_3, &b_fw_sim_geninfo1_3);
   fChain->SetBranchAddress("fw_sim_geninfo2_1", &fw_sim_geninfo2_1, &b_fw_sim_geninfo2_1);
   fChain->SetBranchAddress("fw_sim_geninfo2_2", &fw_sim_geninfo2_2, &b_fw_sim_geninfo2_2);
   fChain->SetBranchAddress("fw_sim_geninfo2_3", &fw_sim_geninfo2_3, &b_fw_sim_geninfo2_3);
   fChain->SetBranchAddress("fw_sim_geninfo_1", &fw_sim_geninfo_1, &b_fw_sim_geninfo_1);
   fChain->SetBranchAddress("fw_sim_geninfo_2", &fw_sim_geninfo_2, &b_fw_sim_geninfo_2);
   fChain->SetBranchAddress("fw_sim_geninfo_3", &fw_sim_geninfo_3, &b_fw_sim_geninfo_3);
   fChain->SetBranchAddress("fw_sim_genweight_1", &fw_sim_genweight_1, &b_fw_sim_genweight_1);
   fChain->SetBranchAddress("fw_sim_genweight_2", &fw_sim_genweight_2, &b_fw_sim_genweight_2);
   fChain->SetBranchAddress("fw_sim_genweight_3", &fw_sim_genweight_3, &b_fw_sim_genweight_3);
   fChain->SetBranchAddress("fw_sim_id_1", &fw_sim_id_1, &b_fw_sim_id_1);
   fChain->SetBranchAddress("fw_sim_id_2", &fw_sim_id_2, &b_fw_sim_id_2);
   fChain->SetBranchAddress("fw_sim_id_3", &fw_sim_id_3, &b_fw_sim_id_3);
   fChain->SetBranchAddress("fw_sim_p_1", &fw_sim_p_1, &b_fw_sim_p_1);
   fChain->SetBranchAddress("fw_sim_p_2", &fw_sim_p_2, &b_fw_sim_p_2);
   fChain->SetBranchAddress("fw_sim_p_3", &fw_sim_p_3, &b_fw_sim_p_3);
   fChain->SetBranchAddress("fw_sim_parentid_1", &fw_sim_parentid_1, &b_fw_sim_parentid_1);
   fChain->SetBranchAddress("fw_sim_parentid_2", &fw_sim_parentid_2, &b_fw_sim_parentid_2);
   fChain->SetBranchAddress("fw_sim_parentid_3", &fw_sim_parentid_3, &b_fw_sim_parentid_3);
   fChain->SetBranchAddress("fw_size_1", &fw_size_1, &b_fw_size_1);
   fChain->SetBranchAddress("fw_size_2", &fw_size_2, &b_fw_size_2);
   fChain->SetBranchAddress("fw_size_3", &fw_size_3, &b_fw_size_3);
   fChain->SetBranchAddress("fw_spectator_1", &fw_spectator_1, &b_fw_spectator_1);
   fChain->SetBranchAddress("fw_spectator_2", &fw_spectator_2, &b_fw_spectator_2);
   fChain->SetBranchAddress("fw_spectator_3", &fw_spectator_3, &b_fw_spectator_3);
   fChain->SetBranchAddress("fw_theta_1", &fw_theta_1, &b_fw_theta_1);
   fChain->SetBranchAddress("fw_theta_2", &fw_theta_2, &b_fw_theta_2);
   fChain->SetBranchAddress("fw_theta_3", &fw_theta_3, &b_fw_theta_3);
   fChain->SetBranchAddress("fw_time_1", &fw_time_1, &b_fw_time_1);
   fChain->SetBranchAddress("fw_time_2", &fw_time_2, &b_fw_time_2);
   fChain->SetBranchAddress("fw_time_3", &fw_time_3, &b_fw_time_3);
   fChain->SetBranchAddress("fw_time_min_1", &fw_time_min_1, &b_fw_time_min_1);
   fChain->SetBranchAddress("fw_time_min_2", &fw_time_min_2, &b_fw_time_min_2);
   fChain->SetBranchAddress("fw_time_min_3", &fw_time_min_3, &b_fw_time_min_3);
   fChain->SetBranchAddress("fw_tracknr_1", &fw_tracknr_1, &b_fw_tracknr_1);
   fChain->SetBranchAddress("fw_tracknr_12", &fw_tracknr_12, &b_fw_tracknr_12);
   fChain->SetBranchAddress("fw_tracknr_2", &fw_tracknr_2, &b_fw_tracknr_2);
   fChain->SetBranchAddress("fw_tracknr_22", &fw_tracknr_22, &b_fw_tracknr_22);
   fChain->SetBranchAddress("fw_tracknr_3", &fw_tracknr_3, &b_fw_tracknr_3);
   fChain->SetBranchAddress("fw_tracknr_32", &fw_tracknr_32, &b_fw_tracknr_32);
   fChain->SetBranchAddress("fw_x_lab_1", &fw_x_lab_1, &b_fw_x_lab_1);
   fChain->SetBranchAddress("fw_x_lab_2", &fw_x_lab_2, &b_fw_x_lab_2);
   fChain->SetBranchAddress("fw_x_lab_3", &fw_x_lab_3, &b_fw_x_lab_3);
   fChain->SetBranchAddress("fw_y_lab_1", &fw_y_lab_1, &b_fw_y_lab_1);
   fChain->SetBranchAddress("fw_y_lab_2", &fw_y_lab_2, &b_fw_y_lab_2);
   fChain->SetBranchAddress("fw_y_lab_3", &fw_y_lab_3, &b_fw_y_lab_3);
   fChain->SetBranchAddress("fw_z_lab_1", &fw_z_lab_1, &b_fw_z_lab_1);
   fChain->SetBranchAddress("fw_z_lab_2", &fw_z_lab_2, &b_fw_z_lab_2);
   fChain->SetBranchAddress("fw_z_lab_3", &fw_z_lab_3, &b_fw_z_lab_3);
   fChain->SetBranchAddress("geant_sim_genweight_1", &geant_sim_genweight_1, &b_geant_sim_genweight_1);
   fChain->SetBranchAddress("geant_sim_genweight_2", &geant_sim_genweight_2, &b_geant_sim_genweight_2);
   fChain->SetBranchAddress("geant_sim_id_1", &geant_sim_id_1, &b_geant_sim_id_1);
   fChain->SetBranchAddress("geant_sim_id_2", &geant_sim_id_2, &b_geant_sim_id_2);
   fChain->SetBranchAddress("geant_sim_p_1", &geant_sim_p_1, &b_geant_sim_p_1);
   fChain->SetBranchAddress("geant_sim_p_2", &geant_sim_p_2, &b_geant_sim_p_2);
   fChain->SetBranchAddress("geant_sim_px_1", &geant_sim_px_1, &b_geant_sim_px_1);
   fChain->SetBranchAddress("geant_sim_px_2", &geant_sim_px_2, &b_geant_sim_px_2);
   fChain->SetBranchAddress("geant_sim_py_1", &geant_sim_py_1, &b_geant_sim_py_1);
   fChain->SetBranchAddress("geant_sim_py_2", &geant_sim_py_2, &b_geant_sim_py_2);
   fChain->SetBranchAddress("geant_sim_pz_1", &geant_sim_pz_1, &b_geant_sim_pz_1);
   fChain->SetBranchAddress("geant_sim_pz_2", &geant_sim_pz_2, &b_geant_sim_pz_2);
   fChain->SetBranchAddress("hneg_mult", &hneg_mult, &b_hneg_mult);
   fChain->SetBranchAddress("hpos_mult", &hpos_mult, &b_hpos_mult);
   fChain->SetBranchAddress("isBest", &isBest, &b_isBest);
   fChain->SetBranchAddress("lneg_mult", &lneg_mult, &b_lneg_mult);
   fChain->SetBranchAddress("lpos_mult", &lpos_mult, &b_lpos_mult);
   fChain->SetBranchAddress("pim_beta", &pim_beta, &b_pim_beta);
   fChain->SetBranchAddress("pim_beta_new", &pim_beta_new, &b_pim_beta_new);
   fChain->SetBranchAddress("pim_btChargeRing", &pim_btChargeRing, &b_pim_btChargeRing);
   fChain->SetBranchAddress("pim_btChargeSum", &pim_btChargeSum, &b_pim_btChargeSum);
   fChain->SetBranchAddress("pim_btChi2", &pim_btChi2, &b_pim_btChi2);
   fChain->SetBranchAddress("pim_btClusters", &pim_btClusters, &b_pim_btClusters);
   fChain->SetBranchAddress("pim_btMaxima", &pim_btMaxima, &b_pim_btMaxima);
   fChain->SetBranchAddress("pim_btMaximaCharge", &pim_btMaximaCharge, &b_pim_btMaximaCharge);
   fChain->SetBranchAddress("pim_btMaximaChargeShared", &pim_btMaximaChargeShared, &b_pim_btMaximaChargeShared);
   fChain->SetBranchAddress("pim_btMaximaChargeSharedFragment", &pim_btMaximaChargeSharedFragment, &b_pim_btMaximaChargeSharedFragment);
   fChain->SetBranchAddress("pim_btMaximaShared", &pim_btMaximaShared, &b_pim_btMaximaShared);
   fChain->SetBranchAddress("pim_btMaximaSharedFragment", &pim_btMaximaSharedFragment, &b_pim_btMaximaSharedFragment);
   fChain->SetBranchAddress("pim_btMeanDist", &pim_btMeanDist, &b_pim_btMeanDist);
   fChain->SetBranchAddress("pim_btNearbyMaxima", &pim_btNearbyMaxima, &b_pim_btNearbyMaxima);
   fChain->SetBranchAddress("pim_btNearbyMaximaShared", &pim_btNearbyMaximaShared, &b_pim_btNearbyMaximaShared);
   fChain->SetBranchAddress("pim_btPadsClus", &pim_btPadsClus, &b_pim_btPadsClus);
   fChain->SetBranchAddress("pim_btPadsRing", &pim_btPadsRing, &b_pim_btPadsRing);
   fChain->SetBranchAddress("pim_btRingMatrix", &pim_btRingMatrix, &b_pim_btRingMatrix);
   fChain->SetBranchAddress("pim_dedx_mdc", &pim_dedx_mdc, &b_pim_dedx_mdc);
   fChain->SetBranchAddress("pim_dedx_tof", &pim_dedx_tof, &b_pim_dedx_tof);
   fChain->SetBranchAddress("pim_id", &pim_id, &b_pim_id);
   fChain->SetBranchAddress("pim_isBT", &pim_isBT, &b_pim_isBT);
   fChain->SetBranchAddress("pim_isOffVertexClust", &pim_isOffVertexClust, &b_pim_isOffVertexClust);
   fChain->SetBranchAddress("pim_isPrimaryVertex", &pim_isPrimaryVertex, &b_pim_isPrimaryVertex);
   fChain->SetBranchAddress("pim_isUsedVertex", &pim_isUsedVertex, &b_pim_isUsedVertex);
   fChain->SetBranchAddress("pim_isring", &pim_isring, &b_pim_isring);
   fChain->SetBranchAddress("pim_isringmdc", &pim_isringmdc, &b_pim_isringmdc);
   fChain->SetBranchAddress("pim_isringnomatch", &pim_isringnomatch, &b_pim_isringnomatch);
   fChain->SetBranchAddress("pim_isringtrack", &pim_isringtrack, &b_pim_isringtrack);
   fChain->SetBranchAddress("pim_kIsLepton", &pim_kIsLepton, &b_pim_kIsLepton);
   fChain->SetBranchAddress("pim_kIsUsed", &pim_kIsUsed, &b_pim_kIsUsed);
   fChain->SetBranchAddress("pim_mdcinnerchi2", &pim_mdcinnerchi2, &b_pim_mdcinnerchi2);
   fChain->SetBranchAddress("pim_mdcouterchi2", &pim_mdcouterchi2, &b_pim_mdcouterchi2);
   fChain->SetBranchAddress("pim_oa_hadr", &pim_oa_hadr, &b_pim_oa_hadr);
   fChain->SetBranchAddress("pim_oa_lept", &pim_oa_lept, &b_pim_oa_lept);
   fChain->SetBranchAddress("pim_p", &pim_p, &b_pim_p);
   fChain->SetBranchAddress("pim_p_corr_em", &pim_p_corr_em, &b_pim_p_corr_em);
   fChain->SetBranchAddress("pim_p_corr_ep", &pim_p_corr_ep, &b_pim_p_corr_ep);
   fChain->SetBranchAddress("pim_p_corr_p", &pim_p_corr_p, &b_pim_p_corr_p);
   fChain->SetBranchAddress("pim_p_corr_pim", &pim_p_corr_pim, &b_pim_p_corr_pim);
   fChain->SetBranchAddress("pim_p_corr_pip", &pim_p_corr_pip, &b_pim_p_corr_pip);
   fChain->SetBranchAddress("pim_phi", &pim_phi, &b_pim_phi);
   fChain->SetBranchAddress("pim_phi_rich", &pim_phi_rich, &b_pim_phi_rich);
   fChain->SetBranchAddress("pim_pid", &pim_pid, &b_pim_pid);
   fChain->SetBranchAddress("pim_q", &pim_q, &b_pim_q);
   fChain->SetBranchAddress("pim_r", &pim_r, &b_pim_r);
   fChain->SetBranchAddress("pim_resolution", &pim_resolution, &b_pim_resolution);
   fChain->SetBranchAddress("pim_resoultion", &pim_resoultion, &b_pim_resoultion);
   fChain->SetBranchAddress("pim_rich_amp", &pim_rich_amp, &b_pim_rich_amp);
   fChain->SetBranchAddress("pim_rich_centr", &pim_rich_centr, &b_pim_rich_centr);
   fChain->SetBranchAddress("pim_rich_houtra", &pim_rich_houtra, &b_pim_rich_houtra);
   fChain->SetBranchAddress("pim_rich_padnum", &pim_rich_padnum, &b_pim_rich_padnum);
   fChain->SetBranchAddress("pim_rich_patmat", &pim_rich_patmat, &b_pim_rich_patmat);
   fChain->SetBranchAddress("pim_rkchi2", &pim_rkchi2, &b_pim_rkchi2);
   fChain->SetBranchAddress("pim_sector", &pim_sector, &b_pim_sector);
   fChain->SetBranchAddress("pim_shw_sum0", &pim_shw_sum0, &b_pim_shw_sum0);
   fChain->SetBranchAddress("pim_shw_sum1", &pim_shw_sum1, &b_pim_shw_sum1);
   fChain->SetBranchAddress("pim_shw_sum2", &pim_shw_sum2, &b_pim_shw_sum2);
   fChain->SetBranchAddress("pim_sim_geninfo", &pim_sim_geninfo, &b_pim_sim_geninfo);
   fChain->SetBranchAddress("pim_sim_geninfo1", &pim_sim_geninfo1, &b_pim_sim_geninfo1);
   fChain->SetBranchAddress("pim_sim_geninfo2", &pim_sim_geninfo2, &b_pim_sim_geninfo2);
   fChain->SetBranchAddress("pim_sim_genweight", &pim_sim_genweight, &b_pim_sim_genweight);
   fChain->SetBranchAddress("pim_sim_grandparentid", &pim_sim_grandparentid, &b_pim_sim_grandparentid);
   fChain->SetBranchAddress("pim_sim_id", &pim_sim_id, &b_pim_sim_id);
   fChain->SetBranchAddress("pim_sim_iscommon", &pim_sim_iscommon, &b_pim_sim_iscommon);
   fChain->SetBranchAddress("pim_sim_isghost", &pim_sim_isghost, &b_pim_sim_isghost);
   fChain->SetBranchAddress("pim_sim_mediumid", &pim_sim_mediumid, &b_pim_sim_mediumid);
   fChain->SetBranchAddress("pim_sim_p", &pim_sim_p, &b_pim_sim_p);
   fChain->SetBranchAddress("pim_sim_parentid", &pim_sim_parentid, &b_pim_sim_parentid);
   fChain->SetBranchAddress("pim_sim_processid", &pim_sim_processid, &b_pim_sim_processid);
   fChain->SetBranchAddress("pim_sim_px", &pim_sim_px, &b_pim_sim_px);
   fChain->SetBranchAddress("pim_sim_py", &pim_sim_py, &b_pim_sim_py);
   fChain->SetBranchAddress("pim_sim_pz", &pim_sim_pz, &b_pim_sim_pz);
   fChain->SetBranchAddress("pim_sim_vertexx", &pim_sim_vertexx, &b_pim_sim_vertexx);
   fChain->SetBranchAddress("pim_sim_vertexy", &pim_sim_vertexy, &b_pim_sim_vertexy);
   fChain->SetBranchAddress("pim_sim_vertexz", &pim_sim_vertexz, &b_pim_sim_vertexz);
   fChain->SetBranchAddress("pim_system", &pim_system, &b_pim_system);
   fChain->SetBranchAddress("pim_theta", &pim_theta, &b_pim_theta);
   fChain->SetBranchAddress("pim_theta_rich", &pim_theta_rich, &b_pim_theta_rich);
   fChain->SetBranchAddress("pim_tof_mom", &pim_tof_mom, &b_pim_tof_mom);
   fChain->SetBranchAddress("pim_tof_new", &pim_tof_new, &b_pim_tof_new);
   fChain->SetBranchAddress("pim_tof_rec", &pim_tof_rec, &b_pim_tof_rec);
   fChain->SetBranchAddress("pim_track_length", &pim_track_length, &b_pim_track_length);
   fChain->SetBranchAddress("pim_tracklength", &pim_tracklength, &b_pim_tracklength);
   fChain->SetBranchAddress("pim_z", &pim_z, &b_pim_z);
   fChain->SetBranchAddress("runnumber", &runnumber, &b_runnumber);
   fChain->SetBranchAddress("totalmult", &totalmult, &b_totalmult);
   fChain->SetBranchAddress("trigbit", &trigbit, &b_trigbit);
   fChain->SetBranchAddress("trigdec", &trigdec, &b_trigdec);
   fChain->SetBranchAddress("trigdownscale", &trigdownscale, &b_trigdownscale);
   fChain->SetBranchAddress("trigdownscaleflag", &trigdownscaleflag, &b_trigdownscaleflag);
   Notify();
}

Bool_t PimEpEm_ID::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void PimEpEm_ID::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t PimEpEm_ID::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef PimEpEm_ID_cxx
