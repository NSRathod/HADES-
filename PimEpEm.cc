#include "PimEpEm.h"
#include "data.h"
#include <iostream>
#include <TLorentzVector.h>
#include <TVector3.h>
#include "hntuple.h"
#include <TMath.h>


using namespace std;
using namespace PATData;

//TFile *hf = new TFile("/lustre/hades/user/nrathod/PLUTO/TEST/QAplots1.root");
//TH2F * OAvsMass =  (TH2F*)hf->Get("OAvsM");

/////////---------------- BREMSSTRAHLUNG SCALING ------------------///////////////////
//TFile *hf = new TFile("/lustre/nyx/hades/user/nrathod/ANALYSIS_QA/pluto_Brem/PPIM_WORK/FINAL_CAL/Wolf_caln_PimP_Eng_Brem_FINAL_CHECK.root");
//TH1F *Brem_Wgt =  (TH1F*)hf->Get("Wolf_Brem2");
//TFile *hf = new TFile("/lustre/nyx/hades/user/nrathod/ANALYSIS_QA/pluto_Brem/PPIM_WORK/Wolf_caln_PimP_Eng_Brem_FINAL_PS_CHECK.root");
//TH3F *Brem_Wgt =  (TH3F*)hf->Get("epem_4pi_Wolf_wt_Mypt");
//TFile *hf = new TFile("/lustre/nyx/hades/user/nrathod/ANALYSIS_QA/pluto_Brem/PPIM_WORK/Wolf_caln_PimP_Eng_Brem_FINAL_ANGULAR_CHECK1.root");
//TH3F *Brem_Wgt =  (TH3F*)hf->Get("epem_4pi_Wolf_wt_Mypt");
//TFile *hf = new TFile("/lustre/nyx/hades/user/nrathod/ANALYSIS_QA/pluto_Brem/PPIM_WORK/CALn/Wolf_caln_3D_HIST_METHOD_Eng_PimP_Brem_PHASESPACE.root");
//TFile *hf = new TFile("/lustre/nyx/hades/user/nrathod/ANALYSIS_QA/pluto_Brem/PPIM_WORK/CALn/Wolf_caln_3D_HIST_METHOD_Eng_PimP_Brem_ANGULAR.root");

//TFile *hf = new TFile("/lustre/nyx/hades/user/nrathod/ANALYSIS_QA/pluto_Brem/PPIM_WORK/CALn/Wolf_caln_3D_HIST_METHOD_Eng_Pim_CARBON_Brem_PHASESPACE.root");
//TFile *hf = new TFile("/lustre/nyx/hades/user/nrathod/ANALYSIS_QA/pluto_Brem/PPIM_WORK/CALn/Wolf_caln_3D_HIST_METHOD_Eng_Pim_CARBON_Brem_PHASESPACE_3_3.root");
//TH3F *Brem_Wgt =  (TH3F*)hf->Get("epem_4pi_Wolf_wt_Mypt");        

/////////////-LEVEL-2-////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//TFile *hf = new TFile("/lustre/nyx/hades/user/nrathod/ANALYSIS_QA/pluto_Brem/PPIM_WORK/RECHECK_cal/Wolf_caln_3D_HIST_METHOD_Eng_PimP_Brem_PHASESPACE_Kine_Corr_FINAL.root");
//TFile *hf = new TFile("/lustre/nyx/hades/user/nrathod/ANALYSIS_QA/pluto_Brem/PPIM_WORK/Final_Caln/Wolf_caln_3D_HIST_METHOD_Eng_PimP_Brem_ANGULAR_Kine_Corr_FINAL_New_GiBUU.root");
//TFile *hf = new TFile("/lustre/nyx/hades/user/nrathod/ANALYSIS_QA/pluto_Brem/PPIM_WORK/Final_Caln/Wolf_caln_3D_HIST_METHOD_Eng_Pim_CARBON_Brem_PHASESPACE_3_3_Corr_FINAL.root");

///////////////////////////-----LEVEL3-----/////////////////////////////
//TFile *hf = new TFile("/lustre/nyx/hades/user/nrathod/ANALYSIS_QA/pluto_Brem/PPIM_WORK/Final_Caln/Caln_TEST_RUN2.root");
//TFile *hf = new TFile("/lustre/nyx/hades/user/nrathod/ANALYSIS_QA/pluto_Brem/PPIM_WORK/Final_Caln/Caln_CARBON_TEST_RUN2.root");

///////////////////////////-----LEVEL4-----/////////////////////////////
//TFile *hf = new TFile("/lustre/nyx/hades/user/nrathod/ANALYSIS_QA/pluto_Brem/PPIM_WORK/MACRO/RESULTS/Wolf_caln_PimP_BREM_ANGULAR_all_FINAL.root");
TFile *hf = new TFile("/lustre/nyx/hades/user/nrathod/ANALYSIS_QA/pluto_Brem/PPIM_WORK/MACRO/RESULTS/Wolf_caln_PimC_BREM_ANGULAR_all_FINAL.root");
TH3F *Brem_Wgt =  (TH3F*)hf->Get("epem_4pi_Wolf_pt2_wt_Mypt_DIV");

void PimEpEm::Loop()
{
    static long licznik = 0;

   if (fChain == 0) return;


        (*tlo)["ep_mom"] = 0;
        (*tlo)["ep_theta"] = 0;
        (*tlo)["ep_theta_rich"] = 0;
        (*tlo)["ep_phi"] = 0;
        (*tlo)["ep_phi_rich"] = 0;
        (*tlo)["ep_beta"] = 0;
        (*tlo)["em_mom"] = 0;
	(*tlo)["em_theta"] = 0;
        (*tlo)["em_theta_rich"] = 0;
        (*tlo)["em_phi"] = 0;
        (*tlo)["em_phi_rich"] = 0;
	(*tlo)["em_beta"] = 0;
        (*tlo)["oa"] = 0;
        (*tlo)["oa_rich"] = 0;
        (*tlo)["sig"] = 0;
        (*tlo)["ep_m"] = 0;
        (*tlo)["em_m"] = 0;
        (*tlo)["epem_inv_mass"] = 0;
        (*tlo)["epem_inv_mass2"] = 0;
        (*tlo)["epem_miss_mass"] = 0;
        (*tlo)["epem_miss_mass2"] = 0;
        (*tlo)["epem_y"] = 0;
        (*tlo)["epem_pt"] = 0;
        (*tlo)["ep_rich_amp"] = 0;
        (*tlo)["ep_rich_centr"] = 0;
        (*tlo)["ep_rich_padnum"] = 0;
        (*tlo)["ep_rich_patmat"] = 0;
        (*tlo)["ep_rich_houtra"] = 0;
        (*tlo)["em_rich_amp"] = 0;
        (*tlo)["em_rich_centr"] = 0;
        (*tlo)["em_rich_padnum"] = 0;
        (*tlo)["em_rich_patmat"] = 0;
        (*tlo)["em_rich_houtra"] = 0;

	    (*tlo)["eVert_x"] = 0;
	    (*tlo)["eVert_y"] = 0;
	    (*tlo)["eVert_z"] = 0;

        (*tlo)["eVertReco_z"] = -1000.;
        (*tlo)["eVertReco_x"] = -1000.;
        (*tlo)["eVertReco_y"] = -1000.;
        (*tlo)["evtPileupMeta"] = 0.;
        (*tlo)["evtPileupStart"] = 0.;
        (*tlo)["ep_isOffVertexClust"] = 0.;
        (*tlo)["ep_p_corr_ep"] = 0.;
        (*tlo)["em_isOffVertexClust"] = 0.;
        (*tlo)["em_p_corr_em"] = 0.;

        (*tlo)["e1_m"] = 0;
        (*tlo)["e2_m"] = 0;
        (*tlo)["pim_m"] = 0;

        (*tlo)["pimepem_miss_m"] = 0;
        (*tlo)["pimepem_miss_m2"] = 0;
        (*tlo)["pimepem_miss_PT_m"] = 0;
        (*tlo)["pimepem_miss_PT_m2"] = 0;


//        tlo->fill();


   Long64_t nentries = fChain->GetEntries();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      ++licznik;
	  if ((licznik % 100)==0) cout << "Events: " << licznik << endl;


      double F = 1; // 1.006;
      TVector3 v1, v2, v3, v1_B, v2_B, v3_B;
      v1.SetXYZ(F*pim_p_corr_pim*sin(D2R*pim_theta)*cos(D2R*pim_phi),F*pim_p_corr_pim*sin(D2R*pim_theta)*sin(D2R*pim_phi),F*pim_p_corr_pim*cos(D2R*pim_theta));
      v2.SetXYZ(F*ep_p_corr_ep*sin(D2R*ep_theta)*cos(D2R*ep_phi),F*ep_p_corr_ep*sin(D2R*ep_theta)*sin(D2R*ep_phi),F*ep_p_corr_ep*cos(D2R*ep_theta));
      v3.SetXYZ(F*em_p_corr_em*sin(D2R*em_theta)*cos(D2R*em_phi),F*em_p_corr_em*sin(D2R*em_theta)*sin(D2R*em_phi),F*em_p_corr_em*cos(D2R*em_theta));

      v1_B.SetXYZ(pim_sim_px,pim_sim_py,pim_sim_pz);
      v2_B.SetXYZ(ep_sim_px,ep_sim_py,ep_sim_pz);
      v3_B.SetXYZ(em_sim_px,em_sim_py,em_sim_pz);
      
      TVector3 r1, r2;
      r1.SetXYZ(sin(D2R*ep_theta_rich)*cos(D2R*ep_phi_rich),sin(D2R*ep_theta_rich)*sin(D2R*ep_phi_rich),cos(D2R*ep_theta_rich));
      r2.SetXYZ(sin(D2R*em_theta_rich)*cos(D2R*em_phi_rich),sin(D2R*em_theta_rich)*sin(D2R*em_phi_rich),cos(D2R*em_theta_rich));

      pim->SetVectM( v1, 139.56995 );
      e1->SetVectM( v2, 0.51099906 );
      e2->SetVectM( v3, 0.51099906 );
////////////////////////////////////////////////////
      pim_B->SetVectM( v1_B, 139.56995 );
      e1_B->SetVectM( v2_B, 0.51099906 );
      e2_B->SetVectM( v3_B, 0.51099906 );
////////////////////////////////////////////////////  
// PT ----------------------------------------

/*      double FPT = 1;
      double UNIT = 1000.0;
      double pion_px = pt_px_1 * UNIT * FPT - d_pion_mom*pt_px_1/pt_p_1;
      double pion_py = pt_py_1 * UNIT * FPT - d_pion_mom*pt_py_1/pt_p_1;
      double pion_pz = pt_pz_1 * UNIT * FPT - d_pion_mom*pt_pz_1/pt_p_1;
      double pion_p =  pt_p_1 * UNIT * FPT - d_pion_mom;
      double pion_E =  sqrt( pion_p*pion_p + 139.56995*139.56995 );

      proj_PT = new TLorentzVector( pion_px, pion_py, pion_pz, pion_E );
      *beam_PT = *targ + *proj_PT;
      */
// -------------------------------------------

      *gammae1e2 = *e1 + *e2;
      *e1e2 = *e1 + *e2;
      *e1_delta = *e1;
      *e2_delta = *e2;
      *e1e2_miss = *beam - *e1 - *e2;
      //*ppimepem = *p + *pim + *e1 + *e2;
      *pimepem = *pim + *e1 + *e2;
      *pimepem_miss = *beam - *pimepem;
      *pimepem_miss_PT = *beam_PT - *pimepem;
      *ppimepem_miss = *beam_PT - *pimepem;
////////////////////////////////////////////////////
////////////////////////////////////////////////////
      *gammae1e2_B = *e1_B + *e2_B;
      *e1e2_B = *e1_B + *e2_B;
      *e1_delta_B = *e1_B;
      *e2_delta_B = *e2_B;
      *e1e2_miss_B = *beam_B - *e1_B - *e2_B;
      *pimepem_B = *pim_B + *e1_B + *e2_B;
      *pimepem_miss_B = *beam_B - *pimepem_B;
      *pimepem_miss_PT_B = *beam_PT - *pimepem_B;
      *ppimepem_miss_B = *beam_PT - *pimepem_B;
////////////////////////////////////////////////////
      double m2_inv_e1e2_B = gammae1e2_B->M2();
      double m_inv_e1e2_B = gammae1e2_B->M();
      double m_invariant_e1e2_B = gammae1e2_B->M()/1000.;
      double m_inv_e1e2_B_y = gammae1e2_B->Rapidity();
      double m_inv_e1e2_B_pt = gammae1e2_B->Pt()/1000.;
      double m_inv_e1e2_B_pt2 = m_inv_e1e2_B_pt*m_inv_e1e2_B_pt;

      double full_signal_B = pimepem_B->M();
      double missing_mass_B = pimepem_miss_B->M();
      ////////////////////////////////////////////////////  
////////////////////////////////////////////////////

      double em_mom = e1_delta->P();
      double ep_mom = e2_delta->P();
      double m2_inv_e1e2 = gammae1e2->M2();
      double m_inv_e1e2 = gammae1e2->M();
      double m_invariant_e1e2 = gammae1e2->M()/1000.;
      double m_inv_e1e2_y = gammae1e2->Rapidity();
      double m_inv_e1e2_pt = gammae1e2->Pt()/1000.;      
      
      double full_signal = pimepem->M();
      double missing_mass = pimepem_miss->M();
      
      double missing_momentum = pimepem_miss->P();
      double missing_energy = pimepem_miss->E();
      
      double oa = R2D * openingangle(*e1, *e2);
      double oa_rich = R2D * openingangle(r1, r2);
      
      double e1_mass = ep_p*ep_p * (  1. / (ep_beta*ep_beta)  - 1. ) ;
      double e2_mass = em_p*em_p * (  1. / (em_beta*em_beta)  - 1. ) ;
      double pim_mass = pim_p*pim_p * (  1. / (pim_beta*pim_beta)  - 1. ) ;

//---------------------Angular Distribution in RADIAN --------------------------------------------

      double Miss_ang_distribution = pimepem_miss->Theta();
      double pim_ang_distribution = pim->Theta();
      double ep_ang_distribution = e1->Theta();
      double em_ang_distribution = e2->Theta();
      double Gamma_ang_distribution = e1e2->Theta();
      
      //---------------------Angular Distribution in DEGREE --------------------------------------------

      double Miss_ang_dis = pimepem_miss->Theta() * TMath::RadToDeg();
      double pim_ang_dis = pim->Theta() * TMath::RadToDeg();
      double ep_ang_dis = e1->Theta() * TMath::RadToDeg();
      double em_ang_dis = e2->Theta() * TMath::RadToDeg();;
      double Gamma_ang_dis = e1e2->Theta() * TMath::RadToDeg();

      //-------------------------------------------------------------------------------------
      //--------------------- Momentum Distribution -----------------------------------------------------

      double pim_mom_dis = pim->P();
      double Gamma_mom_dis = e1e2->P();
      double brem = ppimepem_miss->P();
      double bremX = ppimepem_miss->Px();
      double bremY = ppimepem_miss->Py();
      double bremZ = ppimepem_miss->Pz();
      

      //-------------------------------------------------------------------------------------


      int binx = Brem_Wgt->GetXaxis()->FindBin(m_invariant_e1e2_B);
      int biny = Brem_Wgt->GetYaxis()->FindBin(m_inv_e1e2_B_y);
      int binz = Brem_Wgt->GetZaxis()->FindBin(m_inv_e1e2_B_pt2);
      //int binx = Brem_Wgt->GetXaxis()->FindBin(m_invariant_e1e2);
      //int biny = ->GetYaxis()->FindBin(oa);
      double bin_content = Brem_Wgt->GetBinContent(binx,biny,binz);
      //cout << "bincontent  :  "<< bin_content <<"   "<<"binX =  "<< binx <<"   "<<"mass =  "<< m_invariant_e1e2_B <<"   "<<pim_sim_genweight<< endl;
      
      ACC = 1.;
      //EFF = bin_content;
      EFF_M2 = 1/TMath::Power(m_inv_e1e2/1000.,2);
      EFF = 1.;
      EFF_B = bin_content;

      //cout << "bincontent  :  "<< bin_content <<"   "<<"binX =  "<< binx <<"   "<<"mass =  "<< m_invariant_e1e2_B <<"   "<<EFF<< endl;
      //-------------------------------------------------------------------------------------
      gammae1e2->Boost(0., 0., -(beam->Beta()));
      e1_delta->Boost(0., 0., -(beam->Beta()));
      e2_delta->Boost(0., 0., -(beam->Beta()));

      e1_delta->Boost( -gammae1e2->Px()/gammae1e2->E(), -gammae1e2->Py()/gammae1e2->E(), -gammae1e2->Pz()/gammae1e2->E());
      e2_delta->Boost( -gammae1e2->Px()/gammae1e2->E(), -gammae1e2->Py()/gammae1e2->E(), -gammae1e2->Pz()/gammae1e2->E());
      //-------------------------------------------------------------------------------------

      gammae1e2_B->Boost(0., 0., -(beam_B->Beta()));
      e1_delta_B->Boost(0., 0., -(beam_B->Beta()));
      e2_delta_B->Boost(0., 0., -(beam_B->Beta()));

      e1_delta_B->Boost( -gammae1e2_B->Px()/gammae1e2_B->E(), -gammae1e2_B->Py()/gammae1e2_B->E(), -gammae1e2_B->Pz()/gammae1e2_B->E());
      e2_delta_B->Boost( -gammae1e2_B->Px()/gammae1e2_B->E(), -gammae1e2_B->Py()/gammae1e2_B->E(), -gammae1e2_B->Pz()/gammae1e2_B->E());

      //-------------------------------------------------------------------------------------
      //cout << "Poczatek obliczen..." << endl;

      double ang_cut = 4.;
      //double ang_cut = 9.;

      double close_cut = 0.;
      //double nonfit_close_cut = -4.;
      //double close_cut = 0.;
      double nonfit_close_cut = 0.;


      insideTarget = 1;

      insideEmS0 = 0;
      insideEmS1 = 0;
      insideEpS0 = 0;
      insideEpS1 = 0;


      NoLeptonE1 = !((ep_oa_lept< close_cut&&ep_oa_lept>0.0) &&ep_oa_lept>nonfit_close_cut );
      NoHadronE1 = 1; // !(ep_oa_hadr< close_cut &&ep_oa_hadr>nonfit_close_cut );
      NoLeptonE2 = !((em_oa_lept< close_cut&&em_oa_lept>0.0) &&em_oa_lept>nonfit_close_cut );
      NoHadronE2 = 1; // !(em_oa_hadr< close_cut &&em_oa_hadr>nonfit_close_cut );
      NoHadronE1 = 1;
      NoHadronE2 = 1;

/*
      NoLeptonE1 = 1;
      NoHadronE1 = 1;
      NoLeptonE2 = 1;
      NoHadronE2 = 1;
*/

      Positron = (((ep_system==0&&insideEpS0==0)||(ep_system==1&&insideEpS1==0)));
      Electron = (((em_system==0&&insideEmS0==0)||(em_system==1&&insideEmS1==0)));

      ElectronPositron = Positron && NoLeptonE1 && NoHadronE1  &&  Electron && NoLeptonE2 && NoHadronE2  &&  insideTarget
	 && ( e1_mass < 5000. && e2_mass < 5000. );
      //ElectronPositron = 1;
      

      if (ElectronPositron && /*(((int)trigbit)&16) &&*/ isBest==1 && oa > 0 && eVertReco_z>-500
          && ep_p > 0 && em_p > 0 && pim_p > 0  ) 
      {
        (*tlo)["ep_mom"] = ep_p;
        (*tlo)["ep_theta"] = ep_theta;
        (*tlo)["ep_theta_rich"] = ep_theta_rich;
        (*tlo)["ep_phi"] = ep_phi;
        (*tlo)["ep_phi_rich"] = ep_phi_rich;
        (*tlo)["ep_beta"] = ep_beta_new;
        (*tlo)["em_mom"] = em_p;
        (*tlo)["em_theta"] = em_theta;
        (*tlo)["em_theta_rich"] = em_theta_rich;
        (*tlo)["em_phi"] = em_phi;
        (*tlo)["em_phi_rich"] = em_phi_rich;
        (*tlo)["em_beta"] = em_beta_new;
        (*tlo)["oa"] = oa;
        (*tlo)["oa_rich"] = oa_rich;
        (*tlo)["sig"] = 1;
        (*tlo)["ep_m"] = e1_mass;
        (*tlo)["em_m"] = e2_mass;
        (*tlo)["epem_inv_mass"] = m_inv_e1e2 / 1000.;
        (*tlo)["epem_inv_mass2"] = m2_inv_e1e2 / 1000000.;
        (*tlo)["epem_miss_mass"] = e1e2_miss->M() / 1000.;
        (*tlo)["epem_miss_mass2"] = e1e2_miss->M2() / 1000000.;
        (*tlo)["epem_y"] = e1e2->Rapidity();
        (*tlo)["epem_pt"] = e1e2->Pt() / 1000.;

        (*tlo)["ep_rich_amp"] = ep_rich_amp;
        (*tlo)["ep_rich_centr"] = ep_rich_centr;
        (*tlo)["ep_rich_padnum"] = ep_rich_padnum;
        (*tlo)["ep_rich_patmat"] = ep_rich_patmat;
        (*tlo)["ep_rich_houtra"] = ep_rich_houtra;
        (*tlo)["em_rich_amp"] = em_rich_amp;
        (*tlo)["em_rich_centr"] = em_rich_centr;
        (*tlo)["em_rich_padnum"] = em_rich_padnum;
        (*tlo)["em_rich_patmat"] = em_rich_patmat;
        (*tlo)["em_rich_houtra"] = em_rich_houtra;

   	    (*tlo)["eVert_x"] = eVert_x;
	    (*tlo)["eVert_y"] = eVert_y;
	    (*tlo)["eVert_z"] = eVert_z;

        (*tlo)["eVertReco_z"] = eVertReco_z;
        (*tlo)["eVertReco_x"] = eVertReco_x;
        (*tlo)["eVertReco_y"] = eVertReco_y;

        (*tlo)["evtPileupMeta"] = evtPileupMeta;
        (*tlo)["evtPileupStart"] = evtPileupStart;
        (*tlo)["ep_isOffVertexClust"] = ep_isOffVertexClust;
        (*tlo)["ep_p_corr_ep"] = ep_p_corr_ep;
        (*tlo)["em_isOffVertexClust"] = em_isOffVertexClust;
        (*tlo)["em_p_corr_em"] = em_p_corr_em;

        (*tlo)["e1_m"] = e1_mass;
        (*tlo)["e2_m"] = e2_mass;
        (*tlo)["pim_m"] = pim_mass;

        (*tlo)["pimepem_miss_m"] = pimepem_miss->M();
        (*tlo)["pimepem_miss_m2"] = pimepem_miss->M2();
        (*tlo)["pimepem_miss_PT_m"] = pimepem_miss_PT->M();
        (*tlo)["pimepem_miss_PT_m2"] = pimepem_miss_PT->M2();


        tlo->fill();

      }
      /*      Int_t btnumpads = 1;
      Int_t richnumpads = 3;

      if (ElectronPositron && isBest==1 && eVertReco_z>-500 && oa > ang_cut && ep_btMaxima >= 2 && em_btMaxima >= 2
	&& (ep_btPadsRing>btnumpads || ep_rich_padnum>richnumpads) && (em_btPadsRing>btnumpads || em_rich_padnum>richnumpads)
	)
     */
      Int_t btnumpads = 1;
      Int_t richnumpads = 3;

//##########################################################################################################################################################
//##########################################################################################################################################################
      
      if (ElectronPositron && isBest==1 && eVertReco_z>-500 && (ep_btMaxima >= 2 || ep_rich_padnum > richnumpads) &&
	  (em_btMaxima >= 2 || em_rich_padnum > richnumpads) && em_beta > 0.8 && ep_beta > 0.8 && ep_sim_id == 2 && em_sim_id == 3 && pim_sim_id == 9 )
	{
	  if(oa > 0) {
	    Dilepton_OA_0deg->Fill(oa,EFF_B);
	    epem_Mass_Brem_0deg->Fill(m_inv_e1e2/1000., EFF_B);}
	  if(oa > 1) {
	    Dilepton_OA_1deg->Fill(oa,EFF_B);
	    epem_Mass_Brem_1deg->Fill(m_inv_e1e2/1000., EFF_B);}
	  if(oa > 2) {
	    Dilepton_OA_2deg->Fill(oa,EFF_B);
	    epem_Mass_Brem_2deg->Fill(m_inv_e1e2/1000., EFF_B);}
	  if(oa > 3) {
	    Dilepton_OA_3deg->Fill(oa,EFF_B);
	    epem_Mass_Brem_3deg->Fill(m_inv_e1e2/1000., EFF_B);}
	  if(oa > 4) {
	    Dilepton_OA_4deg->Fill(oa,EFF_B);
	    epem_Mass_Brem_4deg->Fill(m_inv_e1e2/1000., EFF_B);}
	  if(oa > 5) {
	    Dilepton_OA_5deg->Fill(oa,EFF_B);
	    epem_Mass_Brem_5deg->Fill(m_inv_e1e2/1000., EFF_B);}
	  if(oa > 6) {
	    Dilepton_OA_6deg->Fill(oa,EFF_B);
	    epem_Mass_Brem_6deg->Fill(m_inv_e1e2/1000., EFF_B);}
	  if(oa > 7) {
	    Dilepton_OA_7deg->Fill(oa,EFF_B);
	    epem_Mass_Brem_7deg->Fill(m_inv_e1e2/1000., EFF_B);}
	  if(oa > 8) {
	    Dilepton_OA_8deg->Fill(oa,EFF_B);
	    epem_Mass_Brem_8deg->Fill(m_inv_e1e2/1000., EFF_B);}
	  if(oa > 9) {
	    Dilepton_OA_9deg->Fill(oa,EFF_B);
	    epem_Mass_Brem_9deg->Fill(m_inv_e1e2/1000., EFF_B);}
	  if(oa > 10) {
	    Dilepton_OA_10deg->Fill(oa,EFF_B);
	    epem_Mass_Brem_10deg->Fill(m_inv_e1e2/1000., EFF_B);}



	  Ep_OA->Fill(ep_oa_lept,EFF_B);
	  Em_OA->Fill(em_oa_lept,EFF_B);
	  Ep_Mom->Fill(ep_p,EFF_B);
	  Em_Mom->Fill(em_p,EFF_B);

	}

//##########################################################################################################################################################
//##########################################################################################################################################################



      //if ( ElectronPositron && isBest==1 && eVertReco_z>-500 && oa > ang_cut && ep_btMaxima >= 2 && em_btMaxima >= 2 &&
      //(ep_btPadsRing > btnumpads || ep_rich_padnum > richnumpads) && (em_btPadsRing > btnumpads || em_rich_padnum >  richnumpads) &&
      //em_beta > 0.8 && ep_beta > 0.8 && ep_sim_id == 2 && em_sim_id == 3 && pim_sim_id == 9 )

      //if ( ElectronPositron && isBest==1 && eVertReco_z>-500 && oa > ang_cut
      //	   && em_beta > 0.8 && ep_beta > 0.8 && ep_sim_id == 2 && em_sim_id == 3 && pim_sim_id == 9 ) 

      //if ( ElectronPositron && isBest==1 && eVertReco_z>-500 && oa > ang_cut &&  ep_rich_padnum > richnumpads && em_rich_padnum >  richnumpads &&
      //em_beta > 0.8 && ep_beta > 0.8 && ep_sim_id == 2 && em_sim_id == 3 && pim_sim_id == 9 )

      if (ElectronPositron && isBest==1 && eVertReco_z>-500 && oa > ang_cut /*&& ep_btMaxima >= 2 && em_btMaxima >= 2*/
       && (ep_btMaxima >= 2 || ep_rich_padnum > richnumpads) && (em_btMaxima >= 2 || em_rich_padnum > richnumpads)   
      && em_beta > 0.8 && ep_beta > 0.8 && ep_sim_id == 2 && em_sim_id == 3 && pim_sim_id == 9 )

      //if (ElectronPositron && isBest==1 && eVertReco_z>-500 && oa > ang_cut /*&& ep_btMaxima >= 2 && em_btMaxima >= 2*/
       //  && (ep_btMaxima >= 1 /*|| ep_rich_padnum>richnumpads*/) && (em_btMaxima >= 1 /*|| em_rich_padnum>richnumpads*/) && ep_mom > 100. && em_mom > 100.
      // && ep_sim_id == 2 && em_sim_id == 3 && pim_sim_id == 9 )
	{
//##################################################################################################
   double pi_energy = 0;
   //double mom_start = 790.;
   //double mom_start = 738.;
   double mom_start = 680.;
   //double mom_start = 654.;
   double mom_step = 0.05;
   mom_start -= mom_step;

   
   for (int i=0; i<400; ++i) {

      mom_start += mom_step;
      pi_energy = sqrt( mom_start*mom_start + 139.56995*139.56995 );
      proj_var = new TLorentzVector(0,0, mom_start, pi_energy); // PION BEAM momentum as above

      *beam_var = *proj_var + *targ;
      mom_miss->Fill( mom_start, ( *beam_var - *pim - *e1 - *e2 ).M() );
      delete proj_var;

   }

//##################################################################################################

          sig_all->Fill(m_inv_e1e2/1000., EFF );
	  epem_Mass->Fill(m_inv_e1e2/1000., EFF);
	  epem_Mass_Brem->Fill(m_inv_e1e2/1000., EFF_B);
	  Full_signal->Fill(full_signal/1000., EFF);
	  Missing_mass->Fill(missing_mass/1000., EFF);
          Miss_Mass->Fill(missing_mass/1000., EFF_M2);
	  //cout << "entries = " << jentry<<"   "<< m_inv_e1e2_B/1000. << endl;

	  sig_all_MeV->Fill(m_inv_e1e2, EFF );
          EpEm_inv_mass->Fill(m_inv_e1e2/1000., EFF);
	  EpEm_inv_mass1->Fill(m_inv_e1e2/1000., EFF);

          EpEm_inv_mass1_M2->Fill(m_inv_e1e2/1000., EFF_M2);
	  Missing_mass11_M2->Fill(missing_mass/1000., EFF_M2);
          EpEm_inv_mass1_ph->Fill(m_inv_e1e2/1000., EFF);
	  Missing_mass11_ph->Fill(missing_mass/1000., EFF);
	  EpEm_inv_mass1_Brem->Fill(m_inv_e1e2_B/1000., EFF_B);
	  Missing_mass11_Brem->Fill(missing_mass_B/1000., EFF_B);
	  
	  if(m_inv_e1e2/1000. > 0.100){
	    EpEm_inv_mass11->Fill(m_inv_e1e2/1000., EFF);
	    Missing_mass111->Fill(missing_mass/1000., EFF);

            EpEm_inv_mass11_M2->Fill(m_inv_e1e2/1000., EFF_M2);
	    Missing_mass111_M2->Fill(missing_mass/1000., EFF_M2);
	    EpEm_inv_mass11_ph->Fill(m_inv_e1e2_B/1000., EFF);
	    Missing_mass111_ph->Fill(missing_mass_B/1000., EFF);
	    EpEm_inv_mass11_Brem->Fill(m_inv_e1e2_B/1000., EFF_B);
	    Missing_mass111_Brem->Fill(missing_mass_B/1000., EFF_B);
            miss_mass_var->Fill(missing_mass_B/1000., EFF_B);
          rapidity_140_all->Fill( e1e2->Rapidity(), EFF_B  );
          pt_140_all->Fill( e1e2->Pt() / 1000., EFF_B  );

            if( (missing_mass/1000.) > 0.8 && (missing_mass/1000.) < 1.15 ) EpEm_inv_mass_Brem_reg->Fill(m_inv_e1e2_B/1000., EFF_B);
   
	  }
	  
	  if( (missing_mass/1000.) > 0.8 && (missing_mass/1000.) < 1.20 ) EpEm_inv_mass_Brem_nopi0cut->Fill(m_inv_e1e2/1000., EFF_B);
	  if( (missing_mass/1000.) > 0.8 && (missing_mass/1000.) < 1.20 ) sig_all_var_Brem_nopi0cut->Fill(m_inv_e1e2/1000., EFF_B);
	  if( (missing_mass/1000.) > 0.8 && (missing_mass/1000.) < 1.20 ) Missing_mass_Brem_nopi0cut->Fill(missing_mass/1000., EFF_B);

	  
	  //-----------------------Bremsstrahlung reaction ---------------------------------------------------

	  Missing_energy->Fill(missing_energy, EFF);
	  Missing_mass1->Fill(missing_mass/1000., EFF);
	  Missing_mass11->Fill(missing_mass/1000., EFF);
	  Missing_mass2->Fill(missing_mass/1000., EFF);
	  //Missing_mass2->GetBinContent(jentry)/Missing_mass2->GetBinWidth(jentry);
	  Missing_momentum->Fill(missing_momentum, EFF);

	  //---------------------Angular Distribution in Cos_Theta --------------------------------------------
	  
	  Miss_ang_dis_CosTheta->Fill(cos(Miss_ang_dis), EFF);
	  pim_ang_dis_CosTheta->Fill(cos(pim_ang_dis), EFF); 
	  ep_ang_dis_CosTheta->Fill(cos(ep_ang_dis), EFF);
	  em_ang_dis_CosTheta->Fill(cos(em_ang_dis), EFF);
	  Gamma_ang_dis_CosTheta->Fill(cos(Gamma_ang_dis), EFF);

	  //---------------------Angular Distribution in DEGREE ----------------------------------------------

	  Miss_ang_dis_Theta->Fill(Miss_ang_dis, EFF);
	  pim_ang_dis_Theta->Fill(pim_ang_dis, EFF);
	  ep_ang_dis_Theta->Fill(ep_ang_dis, EFF); 
	  em_ang_dis_Theta->Fill(em_ang_dis, EFF); 
	  Gamma_ang_dis_Theta->Fill(Gamma_ang_dis, EFF);

	  //------------------------------ PLOTS FOR EFFICIENCY ----------------------------------------------

	  pim_P_dis->Fill(pim_mom_dis, EFF);
	  Gamma_P_dis->Fill(Gamma_mom_dis, EFF);

	  EpEm_Pim_P->Fill(pim_mom_dis, Gamma_mom_dis);
	  EpEm_Pim_Theta->Fill(pim_ang_dis, Gamma_ang_dis);
	  
	  //--------------------------------------------------------------------------------------------------

	  Four_particle_MM->Fill(pimepem->P());
	  
	  //----------------------- Momentum Distribution ----------------------------------------------------

	  Missing_TOTP->Fill(brem, EFF);
	  Missing_Px->Fill(bremX, EFF);
	  Missing_Py->Fill(bremY, EFF);
	  Missing_Pz->Fill(bremZ, EFF);
	  
	  //--------------------------------------------------------------------------------------------------
	  
         sig_all_var->Fill(m_inv_e1e2/1000., EFF );
         sig_all_var_Brem->Fill(m_inv_e1e2_B/1000., EFF_B );
         sig_all_var2->Fill(m_inv_e1e2/1000., EFF );
         rapidity_all->Fill( e1e2->Rapidity(), EFF_B  );
         pt_all->Fill( e1e2->Pt() / 1000. , EFF_B );

         if(m_inv_e1e2>140.)  miss_all->Fill(e1e2_miss->M()/1000., EFF );

         if (m_inv_e1e2 > 140. && e1e2_miss->M()>860.&& e1e2_miss->M()<1020.)
         {
          cos_ep->Fill( cos( openingangle( *e1_delta, *gammae1e2 ) ) );
          cos_em->Fill( cos( openingangle( *e2_delta, *gammae1e2 ) ) );
          cos_sum->Fill( cos( openingangle( *e1_delta, *gammae1e2 ) ), 0.5 );
          cos_sum->Fill( cos( openingangle( *e2_delta, *gammae1e2 ) ), 0.5 );
          cos_ep_cm->Fill( cos(gammae1e2->Theta() ));
         }
      }
      //tlo->fill();

      // if (Cut(ientry) < 0) continue;
   } // end of main loop
} // eof Loop 



// -------------------------------------------------------------------------------------------------
PimEpEm::PimEpEm(TTree *tree) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.

   em_acc = em_acc_err = ep_acc = ep_acc_err = 0.;
   em_eff = em_eff_err = ep_eff = ep_eff_err = 0.;

   if (tree == 0) {
	  
      TChain * chain = new TChain("PimEpEm_ID","");

      //chain->Add("/lustre/hades/user/przygoda/PATBT/out/sim/LEPTONS/PWADeltaDalitz_22M/delta_dalitz_H2_22.root/EpEm_ID"); //PWA 2020 HADES DeltaDalitzevents--GOOD
      //chain->Add("/lustre/nyx/hades/user/przygoda/PATBT/out/sim/LEPTONS/PWApi0Dalitz_22M/pi0_dalitz_H2_22.root/EpEm_ID"); //PWA 2020 HADES events--GOOD           

      //----------- NEW PATIRAL WAVE ANALYSIS 2020------------(PLEASE USE THIS)-----------------------------

      //chain->Add("/lustre/hades/user/przygoda/PATBT/out/sim/LEPTONS/PWADeltaDalitz_22M/delta_dalitz_H2_22.root/PimEpEm_ID"); //PWA 2020 HADES events--GOOD
            
      //chain->Add("/lustre/nyx/hades/user/przygoda/PATBT/out/sim/LEPTONS/PWApi0Dalitz_22M/pi0_dalitz_H2_22.root/PimEpEm_ID"); //PWA 2020 HADES events--GOOD

//----------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------

      //chain->Add("/lustre/nyx/hades/user/przygoda/PATBT/out/sim/LEPTONS/PWADeltaDalitz_22M/delta_dalitz_H_22.root/PimEpEm_ID"); //PWA 2020 HADES events

      //chain->Add("/lustre/nyx/hades/user/przygoda/PATBT/out/sim/LEPTONS/PWADeltaDalitz_22M/delta_dalitz_MH_22.root/PimEpEm_ID"); //PWA 2020 MANLEY_HADES events

      //---------------------------------------------------------------
      //---------------------------------------------------------------
      
      //chain->Add("/lustre/nyx/hades/user/nrathod/PAT3_FT/outputs/Events_1760/PimP_690MeV_1760_out/hgeantout_PimP_690_PAT3out.root/PimEpEm_ID");
      
      // -- Full Data --------------------------------------------------------------------------------------------------------
      //chain->Add("/lustre/nyx/hades/user/przygoda/PATBT/out/sim/LEPTONS/PWApi0Dalitz/pwa_PWApi0Dalitz_200k.root/PimEpEm_ID");  //Important PWA OLD corrected---DONT USE THIS FILE

      //chain->Add("/lustre/nyx/hades/user/przygoda/PATBT/out/sim/LEPTONS/PWApi0Dalitz1M/pwa_PWApi0Dalitz_1M.root/PimEpEm_ID");  //Important PWA NEW corrected
      
      ///lustre/nyx/hades/user/przygoda/PATBT/out/sim/LEPTONS/PWApi0Dalitz1M/pwa_PWApi0Dalitz_1M.root/PPimEpEm_ID

      //chain->Add("/lustre/nyx/hades/user/przygoda/PATBT/out/sim/LEPTONS/bremss/brake89.root/PimEpEm_ID");                 // Bremsstrahlung reaction

//----------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------

      //chain->Add("/lustre/nyx/hades/user/przygoda/PATBT/out/sim/LEPTONS/bremss_an/bremss_angular22.root/PimEpEm_ID");
  // Bremsstrahlung angular reaction
      //chain->Add("/lustre/nyx/hades/user/przygoda/PATBT/out/sim/LEPTONS/bremss_ps/bremss_phasespace22.root/PimEpEm_ID");
  // Bremsstrahlung phasespace reaction

      //chain->Add("/lustre/hades/user/przygoda/PATBT/out/sim/LEPTONS/PIMC_ANG/Pim12C_Brem_ANGULAR_22.root/PimEpEm_ID");
  // Brems CARBON angular reaction
      chain->Add("/lustre/hades/user/przygoda/PATBT/out/sim/LEPTONS/PIMC_FLAT/Pim12C_Brem_PHASESPACE_22.root/PimEpEm_ID");
  // Brems CARBON phasespace reaction

//----------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------

      //chain->Add("/lustre/hades/user/przygoda/PATBT/out/sim/LEPTONS/bremss_n/n_phasespace22.root/PimEpEm_ID");

      
      //------------------------------------------ PLUTO SIM PROTON  Target----------------------------------------------------

      //chain->Add("/lustre/nyx/hades/user/przygoda/PATBT/out/sim/LEPTONS/Delta0_large_stat/Delta0_large_stat.root/PimEpEm_ID");         // Important Channel
      //chain->Add("/lustre/nyx/hades/user/przygoda/PATBT/out/sim/LEPTONS/Delta1_large_stat/Delta1_large_stat.root/PimEpEm_ID");         // Important Channel               

      //------------------------------------------ PLUTO SIM CARBON  Target----------------------------------------------------

      //chain->Add(" /lustre/nyx/hades/user/przygoda/PATBT/out/sim/LEPTONS/Delta0_carbon_large_stat/Delta0_carbon_large_stat.root/PimEpEm_ID");
      //chain->Add(" /lustre/nyx/hades/user/przygoda/PATBT/out/sim/LEPTONS/Delta1_carbon_large_stat/Delta1_carbon_large_stat.root/PimEpEm_ID");


      
      //====================================================================== OLD DATA FILES =================================================================================      

      // -- PE 690 ----------------------------------------
      //chain->Add("/lustre/nyx/hades/user/przygoda/PATBT/out/exp/gen2/LEPTONS/leptons690_PE.root/PimEpEm_ID");
      //chain->Add(" /lustre/nyx/hades/user/przygoda/PATBT/out/sim/LEPTONS/Delta0_large_stat/Delta0_large_stat.root/PimEpEm_ID");
      // -- C 690 ----------------------------------------
      //chain->Add("/lustre/nyx/hades/user/przygoda/PATBT/out/exp/gen2/LEPTONS/leptons690_C.root/PPimEpEm_ID");
      //chain->Add(" /lustre/nyx/hades/user/przygoda/PATBT/out/sim/LEPTONS/Delta1_large_stat/Delta1_large_stat.root/PimEpEm_ID");         // Important Channel 
      //chain->Add(" /lustre/nyx/hades/user/przygoda/PATBT/out/sim/LEPTONS/Delta0_large_stat/Delta0_large_stat.root/PimEpEm_ID");         // Important Channel        
//chain->Add("/lustre/nyx/hades/user/przygoda/PATBT/out/sim/LEPTONS/PWApi0/hgeantout_PWA_momentums-ppimpi0_ev_weight_ALL.root/PimEpEm_ID"); //Important PWA events 1
//chain->Add("/lustre/nyx/hades/user/przygoda/PATBT/out/sim/LEPTONS/PWApi0/hgeantout_PWA_momentums-ppimpi0_2ndSol_dcs_ALL.root/PimEpEm_ID");//Important PWA events 2

      //chain->Add("/lustre/nyx/hades/user/przygoda/PATBT/out/sim/LEPTONS/bremss/brake89.root/PimEpEm_ID");                             // Bremsstrahlung reaction

      //=========================================================================================================================================================================
      
      tree = chain; 
   }

   Init(tree);
}


// -------------------------------------------------------------------------------------------------
PimEpEm::~PimEpEm()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

// -------------------------------------------------------------------------------------------------
Int_t PimEpEm::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
// -------------------------------------------------------------------------------------------------
Long64_t PimEpEm::LoadTree(Long64_t entry)
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

// -------------------------------------------------------------------------------------------------
void PimEpEm::Init(TTree *tree)
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

// -------------------------------------------------------------------------------------------------
Bool_t PimEpEm::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

// -------------------------------------------------------------------------------------------------
void PimEpEm::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}

// -------------------------------------------------------------------------------------------------
Int_t PimEpEm::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

