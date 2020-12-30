#include "data.h"
#include "PPim.h"
#include "PPim_buffer.h"
#include <iostream>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <vector>
#include "hntuple.h"

using namespace std;
using namespace PATData;

//#define ELASTICWIDE
//#define ELASTIC
//#define INELASTIC

void PPim::Loop()
{

    static long licznik = 0;
    static long dcounter = 0;
    static long accumulator = 0;
    std::vector< PPim_buffer >  EB;

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntries();
   /*                   CAREFUiL                                nentries = 1000000; */
   cout << "TOTAL: " << nentries << endl;

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      ++licznik; if ((licznik % 1000000)==0) cout << "Events: " << licznik << endl;

      double liczba_plikow = 1;


      // ************************** MAJOR SELECTION CONDITION ************************
      if (isBest > 0) {

      if ( event < accumulator ) { /* cout << " file reset: " << fCurrent << endl; */ accumulator = 0; }
      if ( event == accumulator )
      {
         if ( EB.size() > 0 )
         {
            bool addPion = true;
            for (unsigned i=0; i<EB.size(); ++i)
            {
               if ( isBest < 0 ) addPion = false;
            }
            if ( EB.size() < 2 && isBest > 0 && addPion ) EB.push_back( PPim_buffer( this ) );
         }
         else if ( isBest > 0 )
         {
            EB.push_back( PPim_buffer( this ) );
         }
      }
      else
      {
         if ( EB.size() == 1 )
         {
            double w1 = 1.0;
            filler( EB[0], w1 );
         }
         else  if ( EB.size() > 1 )
         {
            double w1 = 0.5;
            filler( EB[0], w1 );

            double w2 = 0.5;
            filler( EB[1], w2 );
         }

         EB.clear();
         accumulator = static_cast< long >( event );
         if ( isBest > 0 ) EB.push_back( PPim_buffer( this ) );
      }

      } // eof isBest > 0

   } // end of main event loop
}




// ---------------------------------------------------------------------------------------------------------
PPim::PPim(TTree *tree) : fChain(0) 
{
   //p_acc = p_acc_err = pim_acc = pim_acc_err = 0.;
   //p_eff = p_eff_err = pim_eff = pim_eff_err = 0.;

   if (tree == 0) {

      TChain * chain = new TChain("PPim_ID","");
      //TChain * chain = new TChain("PPim","");

      // --- PION EXP  ----------------------------

      // PE 690 July
      //chain->Add("/lustre/nyx/hades/user/przygoda/PAT2/out/exp/gen2/196/hadron196.root/PPim_ID");
      //chain->Add("/lustre/nyx/hades/user/przygoda/PAT2/out/exp/gen2/197/hadron197.root/PPim_ID");
      //chain->Add("/lustre/nyx/hades/user/przygoda/PAT2/out/exp/gen2/198/hadron198.root/PPim_ID");
      // C 690 July
      //chain->Add("/lustre/nyx/hades/user/przygoda/PAT2/out/exp/gen2/195/hadron195.root/PPim_ID");

      // PE 690 Sept
      // A
      //chain->Add("/lustre/nyx/hades/user/przygoda/PAT2/out/exp/gen2/232/hadron232.root/PPim_ID");
      //chain->Add("/lustre/nyx/hades/user/przygoda/PAT2/out/exp/gen2/233/hadron233_filtered.root/PPim_ID");
      //chain->Add("/lustre/nyx/hades/user/przygoda/PAT2/out/exp/gen2/233/hadron233.root/PPim_ID");
      //chain->Add("/lustre/nyx/hades/user/przygoda/PAT2/out/exp/gen2/234/hadron234.root/PPim_ID");
      //chain->Add("/lustre/nyx/hades/user/przygoda/PAT2/out/exp/gen2/235/hadron235.root/PPim_ID");
      // B
      //chain->Add("/lustre/nyx/hades/user/przygoda/PAT2/out/exp/gen2/236/hadron236.root/PPim_ID");
      //chain->Add("/lustre/nyx/hades/user/przygoda/PAT2/out/exp/gen2/237/hadron237.root/PPim_ID");
      //chain->Add("/lustre/nyx/hades/user/przygoda/PAT2/out/exp/gen2/238/hadron238.root/PPim_ID");
      //chain->Add("/lustre/nyx/hades/user/przygoda/PAT2/out/exp/gen2/239/hadron239.root/PPim_ID");
      // C
      //chain->Add("/lustre/nyx/hades/user/przygoda/PAT2/out/exp/gen2/240/hadron240.root/PPim_ID");
      //chain->Add("/lustre/nyx/hades/user/przygoda/PAT2/out/exp/gen2/241/hadron241.root/PPim_ID");
      //chain->Add("/lustre/nyx/hades/user/przygoda/PAT2/out/exp/gen2/242/hadron242.root/PPim_ID");
      ////chain->Add("/lustre/nyx/hades/user/przygoda/PAT2/out/exp/gen2/243/hadron243_filtered.root/PPim_ID");
      //chain->Add("/lustre/nyx/hades/user/przygoda/PAT2/out/exp/gen2/243/hadron243.root/PPim_ID");
      // D
      //chain->Add("/lustre/nyx/hades/user/przygoda/PAT2/out/exp/gen2/244/hadron244.root/PPim_ID");
      //chain->Add("/lustre/nyx/hades/user/przygoda/PAT2/out/exp/gen2/244/hadron244_filtered.root/PPim_ID");
      //chain->Add("/lustre/nyx/hades/user/przygoda/PAT2/out/exp/gen2/245/hadron245.root/PPim_ID");
      //chain->Add("/lustre/nyx/hades/user/przygoda/PAT2/out/exp/gen2/246/hadron246.root/PPim_ID");
      // GOOD
      //chain->Add("/lustre/nyx/hades/user/przygoda/PAT2/out/exp/gen2/244/hadron244evening.root/PPim_ID");
      //chain->Add("/lustre/nyx/hades/user/przygoda/PAT2/out/exp/gen2/245/hadron245.root/PPim_ID");

      // C 690 Sept
      //chain->Add("/lustre/nyx/hades/user/przygoda/PAT2/out/exp/gen2/250/hadron250.root/PPim_ID");
      //chain->Add("/lustre/nyx/hades/user/przygoda/PAT2/out/exp/gen2/251/hadron251.root/PPim_ID");
      //chain->Add("/lustre/nyx/hades/user/przygoda/PAT2/out/exp/gen2/254/hadron254.root/PPim_ID");
      /////chain->Add("/lustre/nyx/hades/user/przygoda/PAT2/out/exp/gen2/254/hadron254_690_filtered.root/PPim_ID");
      chain->Add("/lustre/nyx/hades/user/przygoda/PAT2/out/exp/gen2/255/hadron255_690_filtered.root/PPim_ID");
      //chain->Add("/lustre/nyx/hades/user/przygoda/PAT2/out/exp/gen2/255/hadron255.root/PPim_ID");
      //chain->Add("/lustre/nyx/hades/user/przygoda/PAT2/out/exp/gen2/256/hadron256.root/PPim_ID");


      // PE 656
      ////chain->Add("/lustre/nyx/hades/user/przygoda/PAT2/out/exp/gen2/247/656/hadron247_656_filtered.root/PPim_ID");
      //chain->Add("/lustre/nyx/hades/user/przygoda/PAT2/out/exp/gen2/248/656/hadron248_656_filtered.root/PPim_ID");
      ////chain->Add("/lustre/nyx/hades/user/przygoda/PAT2/out/exp/gen2/249/656/hadron249_656_filtered.root/PPim_ID");
      // C 656
      //chain->Add("/lustre/nyx/hades/user/przygoda/PAT2/out/exp/gen2/252/656/hadron252_656_filtered.root/PPim_ID");
      //chain->Add("/lustre/nyx/hades/user/przygoda/PAT2/out/exp/gen2/253/656/hadron253_656_filtered.root/PPim_ID");
      //chain->Add("/lustre/nyx/hades/user/przygoda/PAT2/out/exp/gen2/254/656/hadron254_656_filtered.root/PPim_ID");

      // PE 748
      //////chain->Add("/lustre/nyx/hades/user/przygoda/PAT2/out/exp/gen2/246/748/hadron246_748_filtered.root/PPim_ID");
      //chain->Add("/lustre/nyx/hades/user/przygoda/PAT2/out/exp/gen2/247/748/hadron247_748_filtered.root/PPim_ID");
      // C 
      //chain->Add("/lustre/nyx/hades/user/przygoda/PAT2/out/exp/gen2/252/748/hadron252_748_filtered.root/PPim_ID");

      // PE 800
      // chain->Add("/lustre/nyx/hades/user/przygoda/PAT2/out/exp/gen2/249/800/hadron249_800_filtered.root/PPim_ID");
      // C
      // chain->Add("/lustre/nyx/hades/user/przygoda/PAT2/out/exp/gen2/251/800/hadron251_800_filtered.root/PPim_ID");
      // chain->Add("/lustre/nyx/hades/user/przygoda/PAT2/out/exp/gen2/252/800/hadron252_800.root/PPim_ID");

      // C 612
      // chain->Add("/lustre/nyx/hades/user/przygoda/PAT2/out/exp/gen2/249/612/hadron249_612.root/PPim_ID");
      // chain->Add("/lustre/nyx/hades/user/przygoda/PAT2/out/exp/gen2/250/612/hadron250_612.root/PPim_ID");


      // ----------------------------------------

      tree = chain;
   }

   Init(tree);
}



// ---------------------------------------------------------------------------------------------------------
 void PPim::filler( const PPim_buffer& s, double WEIGHT )
 {

      static long licz = 0;
      //cout << "TEST: " << endl;

      double F=1.000;//0.9943;//1.0031;
      double FF=1.000;//0.9943;//0.979;
      TVector3 v1, v2;
      // CORR
      v1.SetXYZ(FF*s.p_p_corr_p*sin(D2R*s.p_theta)*cos(D2R*s.p_phi),FF*s.p_p_corr_p*sin(D2R*s.p_theta)*sin(D2R*s.p_phi),FF*s.p_p_corr_p*cos(D2R*s.p_theta));
      v2.SetXYZ(F*s.pim_p_corr_pim*sin(D2R*s.pim_theta)*cos(D2R*s.pim_phi),F*s.pim_p_corr_pim*sin(D2R*s.pim_theta)*sin(D2R*s.pim_phi),F*s.pim_p_corr_pim*cos(D2R*s.pim_theta));

      double P_P_REAL = F*s.p_p_corr_p;
      double PIM_P_REAL = F*s.pim_p_corr_pim;

      // NO corr
      //v1.SetXYZ(F*s.p_p*sin(D2R*s.p_theta)*cos(D2R*s.p_phi),F*s.p_p*sin(D2R*s.p_theta)*sin(D2R*s.p_phi),F*s.p_p*cos(D2R*s.p_theta));
      //v2.SetXYZ(F*s.pim_p*sin(D2R*s.pim_theta)*cos(D2R*s.pim_phi),F*s.pim_p*sin(D2R*s.pim_theta)*sin(D2R*s.pim_phi),F*s.pim_p*cos(D2R*s.pim_theta));


      p->SetVectM( v1, 938.27231 );
      pim->SetVectM( v2, 139.56995 );
      *p_CM = *p;
      *pim_CM = *pim;
      *p_CM_PT = *p;
      *pim_CM_PT = *pim;

      double FPT = 1; // for carbon
      //double FPT = 0.9943; // 0.99864;// let's tune 656
      //double FPT = 0.9930; // 0.9963;// let's tune 690
      //double FPT = 0.9928; // 0.9965;// let's tune 748
      //double FPT = 0.9915; // 0.9961;// this is for 800
      double UNIT = 1000.0;
      double pion_px = s.pt_px_1 * UNIT * FPT - d_pion_mom*s.pt_px_1/s.pt_p_1;
      double pion_py = s.pt_py_1 * UNIT * FPT - d_pion_mom*s.pt_py_1/s.pt_p_1;
      double pion_pz = s.pt_pz_1 * UNIT * FPT - d_pion_mom*s.pt_pz_1/s.pt_p_1;
      double pion_p =  s.pt_p_1 * UNIT * FPT - d_pion_mom;
      double pion_E =  sqrt( pion_p*pion_p + 139.56995*139.56995 );

      proj_PT = new TLorentzVector( pion_px, pion_py, pion_pz, pion_E );
      *beam_PT = *targ + *proj_PT;
//cout << beam_PT->M() << endl;


      *pim_proj = *p + *pim - *targ;
      //cout << "pim proj mom " << pim_proj->P() << " ( " << pim_proj->Px() << ","<< pim_proj->Py()<< "," << pim_proj->Pz() <<" )" << endl;
      *ppim = *p + *pim;
      *ppim_CM = *p + *pim;
      *ppim_CM_PT = *p + *pim;
      *ppim_miss = *beam - (*p + *pim);
      *ppim_miss_PT = *beam_PT - (*p + *pim);
      *ppim_miss_CM = *beam - (*p + *pim);
      *ppim_miss_CM_PT = *beam_PT - (*p + *pim);

      p_CM->Boost( -(*beam).BoostVector() ); 
      pim_CM->Boost( -(*beam).BoostVector() ); 
      p_CM_PT->Boost( -(*beam_PT).BoostVector() ); 
      pim_CM_PT->Boost( -(*beam_PT).BoostVector() ); 
      ppim_CM->Boost( -(*beam).BoostVector() ); 
      ppim_CM_PT->Boost( -(*beam_PT).BoostVector() ); 
      ppim_miss_CM->Boost( -(*beam).BoostVector() ); 
      ppim_miss_CM_PT->Boost( -(*beam_PT).BoostVector() ); 

      ACC = EFF = 1.0;

      double oa = R2D * openingangle(*p, *pim);

      double tantan = TMath::Tan( D2R * s.p_theta ) * TMath::Tan( D2R * s.pim_theta );
      double dphi = TMath::Abs( s.p_phi - s.pim_phi );

      // after smearing I can calculate it myself... ok no smearing so far :-)
      int p_sector_calc = static_cast<int> ( transformPhi( p->Phi()*TMath::RadToDeg() ) / 60 ) - 1;
      int pim_sector_calc = static_cast<int> ( transformPhi( pim->Phi()*TMath::RadToDeg() ) / 60 ) - 1;
      if ( p_sector_calc < 0 ) p_sector_calc = 5;
      if ( pim_sector_calc < 0 ) pim_sector_calc = 5;
      int p_system_calc = p->Theta()*TMath::RadToDeg() < 45. ? 0 : 1;
      int pim_system_calc = pim->Theta()*TMath::RadToDeg() < 45. ? 0 : 1;

      double p_inv_mass2 = P_P_REAL*P_P_REAL * (  1. / (s.p_beta_new*s.p_beta_new)  - 1. ) / 1000000.;
      double pim_inv_mass2 = PIM_P_REAL*PIM_P_REAL * (  1. / (s.pim_beta_new*s.pim_beta_new)  - 1. ) / 1000000.;
      double p_inv_mass = sqrt( P_P_REAL*P_P_REAL * (  1. / (s.p_beta_new*s.p_beta_new)  - 1. ) ) / 1000.;
      double pim_inv_mass = sqrt( PIM_P_REAL*PIM_P_REAL * (  1. / (s.pim_beta_new*s.pim_beta_new)  - 1. ) ) / 1000.;

      double ang_cut = 9.;
      double close_cut = 5.;
      double nonfit_close_cut = -5.;

      NoLeptonP = !( s.p_oa_lept < close_cut && s.p_oa_lept > nonfit_close_cut );
      NoHadronP = !( s.p_oa_hadr < close_cut && s.p_oa_hadr > nonfit_close_cut );
      NoLeptonPim = !( s.pim_oa_lept < close_cut && s.pim_oa_lept > nonfit_close_cut );
      NoHadronPim = !( s.pim_oa_hadr < close_cut && s.pim_oa_hadr > nonfit_close_cut );

      if ( s.isBest > 0 &&
           s.pim_rich_padnum<=0&&s.p_rich_padnum<=0&&oa>2&&s.eVertReco_z>-500 && dphi>0.0 && p_inv_mass>0.6 && p_inv_mass<1.2
           // PWA**** ppim_miss->M2()/1000000. > 0.0 && ppim_miss->M2()/1000000. < 0.1 && s.pim_rich_padnum<=0&&s.p_rich_padnum<=0&&oa>2     
           // && tantan > 1.0 && dphi > 175. && dphi < 185. && dphi>0.0 && p_inv_mass>0.6 && p_inv_mass<1.2 // PWA "elastic && pid"
           // && s.p_z < -10. && s.p_z > -90. && s.p_p > 50. && s.pim_z < -10. && s.pim_z > -90. && s.pim_p > 50.
         )
      {

      (*ptrHNT)["p_inv_mass"] = p_inv_mass;
      (*ptrHNT)["pim_inv_mass"] = pim_inv_mass;
      (*ptrHNT)["p_inv_mass2"] = p_inv_mass2;
      (*ptrHNT)["pim_inv_mass2"] = pim_inv_mass2;
      (*ptrHNT)["ppim_inv_mass"] = ppim_CM->M() / 1000.;
      (*ptrHNT)["ppim_miss_mass2"] = ppim_miss->M2() / 1000000.;
      (*ptrHNT)["ppim_miss_mass2_PT"] = ppim_miss_PT->M2() / 1000000.;
      (*ptrHNT)["ppim_costh_cm"] = ppim_CM->CosTheta();
      (*ptrHNT)["ppim_costh_cm_PT"] = ppim_CM_PT->CosTheta();
      (*ptrHNT)["elastic"] = ( tantan > 1.0 && dphi > 175. && dphi < 185. ) ? 1 : 0;
      (*ptrHNT)["tantan"] = tantan;
      (*ptrHNT)["dphi"] = dphi;
      (*ptrHNT)["oa"] = oa;
      (*ptrHNT)["p_beta"] = s.p_beta_new;
      (*ptrHNT)["p_mom"] = F*s.p_p;
      (*ptrHNT)["p_theta"] = s.p_theta;
      (*ptrHNT)["p_phi"] = s.p_phi;
      (*ptrHNT)["p_sector"] = s.p_sector;
      (*ptrHNT)["p_sector_calc"] = p_sector_calc;
      (*ptrHNT)["pim_beta"] = s.pim_beta_new;
      (*ptrHNT)["pim_mom"] = F*s.pim_p;
      (*ptrHNT)["pim_theta"] = s.pim_theta;
      (*ptrHNT)["pim_phi"] = s.pim_phi;
      (*ptrHNT)["pim_sector"] = s.pim_sector;
      (*ptrHNT)["pim_sector_calc"] = pim_sector_calc;
      (*ptrHNT)["pid"] = dphi>0.0 && p_inv_mass>0.6 && p_inv_mass<1.2;

      (*ptrHNT)["pt_mult"] = s.pt_mult;
      (*ptrHNT)["pt_p"] = s.pt_p_1;
      (*ptrHNT)["pt_theta"] = s.pt_theta_1;
      (*ptrHNT)["pt_phi"] = s.pt_phi_1;
      (*ptrHNT)["pt_costh"] = s.pt_costh_1;
      (*ptrHNT)["pt_x1"] = s.pt_x1_1;
      (*ptrHNT)["pt_y1"] = s.pt_y1_1;
      (*ptrHNT)["pt_x2"] = s.pt_x2_1;
      (*ptrHNT)["pt_y2"] = s.pt_y2_1;
      (*ptrHNT)["pt_xtarg"] = s.pt_xtarg_1;
      (*ptrHNT)["pt_ytarg"] = s.pt_ytarg_1;
      (*ptrHNT)["pt_px"] = s.pt_px_1;
      (*ptrHNT)["pt_py"] = s.pt_py_1;
      (*ptrHNT)["pt_pz"] = s.pt_pz_1;
      (*ptrHNT)["pt_match"] = s.pt_match_1;
      (*ptrHNT)["pt_theta0"] = s.pt_theta0_1;
      (*ptrHNT)["pt_phi0"] = s.pt_phi0_1;

      (*ptrHNT)["pt_pion_p"] = pion_p;
      (*ptrHNT)["pt_pion_px"] = pion_px;
      (*ptrHNT)["pt_pion_py"] = pion_py;
      (*ptrHNT)["pt_pion_pz"] = pion_pz;

      (*ptrHNT)["isBest"] = s.isBest;

      (*ptrHNT)["p_mom_cm"] = p_CM->P();
      (*ptrHNT)["pim_mom_cm"] = pim_CM->P();
      (*ptrHNT)["p_costh_cm"] = p_CM->CosTheta();
      (*ptrHNT)["pim_costh_cm"] = pim_CM->CosTheta();
      (*ptrHNT)["p_theta_cm"] = p_CM->Theta();
      (*ptrHNT)["pim_theta_cm"] = pim_CM->Theta();
      (*ptrHNT)["p_mom_cm_PT"] = p_CM_PT->P();
      (*ptrHNT)["pim_mom_cm_PT"] = pim_CM_PT->P();
      (*ptrHNT)["p_costh_cm_PT"] = p_CM_PT->CosTheta();
      (*ptrHNT)["pim_costh_cm_PT"] = pim_CM_PT->CosTheta();
      (*ptrHNT)["p_theta_cm_PT"] = p_CM_PT->Theta();
      (*ptrHNT)["pim_theta_cm_PT"] = pim_CM_PT->Theta();
      (*ptrHNT)["eVertReco_z"] = s.eVertReco_z;
      (*ptrHNT)["eVertReco_x"] = s.eVertReco_x;
      (*ptrHNT)["eVertReco_y"] = s.eVertReco_y;
      (*ptrHNT)["evtPileupMeta"] = s.evtPileupMeta;
      (*ptrHNT)["evtPileupStart"] = s.evtPileupStart;
      (*ptrHNT)["p_isOffVertexClust"] = s.p_isOffVertexClust;
      (*ptrHNT)["p_p_corr_p"] = P_P_REAL;
      (*ptrHNT)["pim_isOffVertexClust"] = s.pim_isOffVertexClust;
      (*ptrHNT)["pim_p_corr_pim"] = PIM_P_REAL;
      (*ptrHNT)["elastic_cut"] = p_cut->IsInside( P_P_REAL*0.001, s.p_theta*TMath::DegToRad() ) && pim_cut->IsInside( PIM_P_REAL*0.001, s.pim_theta*TMath::DegToRad() );
      (*ptrHNT)["el"] = p_cut->IsInside( P_P_REAL*0.001, s.p_theta*TMath::DegToRad() ) && pim_cut->IsInside( PIM_P_REAL*0.001, s.pim_theta*TMath::DegToRad() ) && ( tantan > 1.0 && dphi > 175. && dphi < 185. ) && (s.eVertReco_z>-500.) && (p_inv_mass>0.6 && p_inv_mass<1.2);
      (*ptrHNT)["elpim"] = pim_cut->IsInside( PIM_P_REAL*0.001, s.pim_theta*TMath::DegToRad() ) && ( tantan > 1.0 && dphi > 175. && dphi < 185. ) && (s.eVertReco_z>-500.) && (p_inv_mass>0.6 && p_inv_mass<1.2);

//cout << " p_cut " << s.p_p*0.001 << "  " << s.p_theta << " = " << p_cut->IsInside( s.p_p*0.001, s.p_theta*TMath::DegToRad() ) << endl;
//cout << " pim_cut " << s.pim_p*0.001 << "  " << s.pim_theta << " = " << pim_cut->IsInside( s.pim_p*0.001, s.pim_theta*TMath::DegToRad() ) << endl;

      (*ptrHNT)["p_px"] = p->Px();
      (*ptrHNT)["p_py"] = p->Py();
      (*ptrHNT)["p_pz"] = p->Pz();
      (*ptrHNT)["pim_px"] = pim->Px();
      (*ptrHNT)["pim_py"] = pim->Py();
      (*ptrHNT)["pim_pz"] = pim->Pz();

      ptrHNT->fill();

         double Qfactor = 1.0;
         if ( s.isBest > 0 && s.pim_rich_padnum<=0&&s.p_rich_padnum<=0&&oa>2 && s.eVertReco_z>-500 && 
              /*s.pt_mult == 1 && */
#if defined(ELASTICWIDE)
             (tantan > 1.0 && dphi>155. && dphi < 205.) && 
#endif
#if defined(ELASTIC)
             (tantan > 1.0 && dphi>175. && dphi < 185.) && 
#endif
#if defined(INELASTIC)
             !(tantan > 1.0 && dphi>175. && dphi < 185.) && 
#endif
             dphi>0.0 && p_inv_mass>0.6 && p_inv_mass<1.2 ) { //

if (true) {
//if (false) {
//if ( ppim_miss_PT->M() > (134.9764 - 50) && ppim_miss_PT->M() < (134.9764 + 50) ) {
//if ( ppim_miss_PT->M() > (134.9764 - 30) && ppim_miss_PT->M() < (134.9764 + 30) ) {

              //p_cut->IsInside( P_P_REAL*0.001, s.p_theta*TMath::DegToRad() ) && pim_cut->IsInside( PIM_P_REAL*0.001, s.pim_theta*TMath::DegToRad() ) ) { // elastic pre-selection

         //if ( eVertReco_z>-500 && !(tantan > 1.0 && dphi>175. && dphi < 185.) && dphi>0.0 && p_inv_mass>0.6 && p_inv_mass<1.2 ) { // inelasttic pre-selection

         /*
         if ( eVertReco_z>-500 && !( tantan > 1.0 && dphi > 175. && dphi < 185. )  && dphi>0.0 && p_inv_mass>0.6 && p_inv_mass<1.2 
              && (ppim_CM->M() / 1000.) < 1.30 
              && !p_cut->IsInside( s.p_theta, s.p_p ) && !pim_cut->IsInside( s.pim_theta, s.pim_p ) ) {
         */
         // PWA output

         //cout << "          " << ++licz << "  " << "1.0  0.0  0.0" << endl;
         cout << "          " << ++licz << "  " << s.totalmult << " " << s.hpos_mult<<" "<<s.hneg_mult<<" "<<s.lpos_mult<<" "<<s.lneg_mult<<endl;
         cout << "      " << pim->Px() << "      " << pim->Py() << "      " << pim->Pz() << "      " << pim->E() << endl;
         cout << "      " << ppim_miss->Px() << "      " << ppim_miss->Py() << "      " << ppim_miss->Pz() << "      " << ppim_miss->E() << endl;
         ////// ********** WARNING ********** PION TRACKER ADDITIONAL **********
//         cout << "      " << ppim_miss_PT->Px() << "      " << ppim_miss_PT->Py() << "      " << ppim_miss_PT->Pz() << "      " << ppim_miss_PT->E() << endl;
         cout << "      " << p->Px() << "      " << p->Py() << "      " << p->Pz() << "      " << p->E() << endl;

 }

         //hbeam_pt->Fill( beam_PT->P() * 0.001 );
         hbeam_pt->Fill( proj_PT->P() * 0.001 );
         hbeam_pt_s->Fill( beam_PT->M() );
         if ( 85 < ppim_miss->M() && 185 > ppim_miss->M() ) {
             hbeam_pt_win->Fill( beam_PT->P() * 0.001 );
             hbeam_pt_s_win->Fill( beam_PT->M() );
         }
         hbeam->Fill( (beam_PT->P() - beam->P()) * 0.001 );
         hbeam_rec->Fill( pim_proj->P() * 0.001 );
         mom_miss_pt->Fill( proj_PT->P(), ( *ppim_miss_PT ).M2() );

         }

         //cout << "P mom theta " << P_P_REAL*0.001 << " = P() = " << p->P() << " " <<  s.p_theta*TMath::DegToRad() << " = Theta() = " << p->Theta() << endl;


         double p_PID = p_inv_mass>0.8 && p_inv_mass< 1.1 && dphi>0 && p_beta>0.1;

          //if (WEIGHT!=1) cout << " masa  " << WEIGHT << endl;
	  //if ( ! ( tantan > 1.0 && dphi > 175. && dphi < 185. ) ) {
  if ( s.isBest > 0 && s.pim_rich_padnum<=0&&s.p_rich_padnum<=0&&oa>2 && s.eVertReco_z>-500 && 
       /*s.pt_mult == 1 && */
#if defined(ELASTICWIDE)
            (tantan > 1.0 && dphi>155. && dphi < 205.) && 
#endif
#if defined(ELASTIC)
            (tantan > 1.0 && dphi>175. && dphi < 185.) && 
#endif
#if defined(INELASTIC)
            !(tantan > 1.0 && dphi>175. && dphi < 185.) && 
#endif
            dphi>0.0 && p_inv_mass>0.6 && p_inv_mass<1.2 ) { 
       
//##################################################################################################
   double pi_energy = 0;
   //double mom_start = 790.;
   //double mom_start = 738.;
   double mom_start = 680.;
   //double mom_start = 654.;
   double mom_step = 0.05;
   mom_start -= mom_step;

/*
   for (int i=0; i<400; ++i) {

      mom_start += mom_step;
      pi_energy = sqrt( mom_start*mom_start + 139.56995*139.56995 );
      proj_var = new TLorentzVector(0,0, mom_start, pi_energy); // PION BEAM momentum as above

      *beam_var = *proj_var + *targ;
      mom_miss->Fill( mom_start, ( *beam_var - *p - *pim ).M2() );
      delete proj_var;

   }
*/
//##################################################################################################
	  
         p_inv_m->Fill( p_inv_mass, EFF*WEIGHT );
         pim_inv_m->Fill( pim_inv_mass, EFF*WEIGHT );

         if ( ppim_miss_CM->M2() / 1000000. > 0 ) ppim_inv_mass->Fill( ppim->M() / 1000., EFF*WEIGHT );
         ppim_inv_mass_IDEAL->Fill( ppim_CM->M() / 1000., EFF*WEIGHT );
         //if ( ppim->M() / 1000. < 1.35 ) ppim_miss_mass->Fill( ppim_miss->M2() / 1000000., EFF*WEIGHT );
         ppim_miss_mass->Fill( ppim_miss->M2() / 1000000., EFF*WEIGHT );
         ppim_miss_mass_PT->Fill( ppim_miss_PT->M2() / 1000000., EFF*WEIGHT );
         ppim_miss_mass_IDEAL->Fill( ppim_miss_CM->M2() / 1000000., EFF*WEIGHT );

         ppim_miss_m->Fill( ppim_miss->M()/1000., EFF*WEIGHT );
         ppim_miss_m2->Fill( ppim_miss->M2(), EFF*WEIGHT );
//if ( ppim_miss_PT->M() > (134.9764 - 30) && ppim_miss_PT->M() < (134.9764 + 30) ) {
         //ppim_miss_m_PT->Fill( ppim_miss_PT->M(), EFF*WEIGHT );
         ppim_miss_m_PT->Fill( ppim_miss_PT->M()/1000. );
//}
         ppim_miss_m2_PT->Fill( ppim_miss_PT->M2(), EFF*WEIGHT );

         ppim_invmass_missmass->Fill( ppim->M() / 1000., ppim_miss->M2() / 1000000., EFF*WEIGHT );
         ppim_costheta_missmass->Fill( ppim_CM->CosTheta(), ppim_miss->M2() / 1000000., EFF*WEIGHT );
  }

         tantan_dphi->Fill( dphi, tantan, EFF*WEIGHT );

      }



 }


// ---------------------------------------------------------------------------------------------------------
PPim::~PPim()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

// ---------------------------------------------------------------------------------------------------------
Int_t PPim::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

// ---------------------------------------------------------------------------------------------------------
Long64_t PPim::LoadTree(Long64_t entry)
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

// ---------------------------------------------------------------------------------------------------------
void PPim::Init(TTree *tree)
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
   fChain->SetBranchAddress("fw_x_lab_1", &fw_x_lab_1, &b_fw_x_lab_1);
   fChain->SetBranchAddress("fw_x_lab_2", &fw_x_lab_2, &b_fw_x_lab_2);
   fChain->SetBranchAddress("fw_x_lab_3", &fw_x_lab_3, &b_fw_x_lab_3);
   fChain->SetBranchAddress("fw_y_lab_1", &fw_y_lab_1, &b_fw_y_lab_1);
   fChain->SetBranchAddress("fw_y_lab_2", &fw_y_lab_2, &b_fw_y_lab_2);
   fChain->SetBranchAddress("fw_y_lab_3", &fw_y_lab_3, &b_fw_y_lab_3);
   fChain->SetBranchAddress("fw_z_lab_1", &fw_z_lab_1, &b_fw_z_lab_1);
   fChain->SetBranchAddress("fw_z_lab_2", &fw_z_lab_2, &b_fw_z_lab_2);
   fChain->SetBranchAddress("fw_z_lab_3", &fw_z_lab_3, &b_fw_z_lab_3);
   fChain->SetBranchAddress("hneg_mult", &hneg_mult, &b_hneg_mult);
   fChain->SetBranchAddress("hpos_mult", &hpos_mult, &b_hpos_mult);
   fChain->SetBranchAddress("isBest", &isBest, &b_isBest);
   fChain->SetBranchAddress("lneg_mult", &lneg_mult, &b_lneg_mult);
   fChain->SetBranchAddress("lpos_mult", &lpos_mult, &b_lpos_mult);
   fChain->SetBranchAddress("p_beta", &p_beta, &b_p_beta);
   fChain->SetBranchAddress("p_beta_new", &p_beta_new, &b_p_beta_new);
   fChain->SetBranchAddress("p_dedx_mdc", &p_dedx_mdc, &b_p_dedx_mdc);
   fChain->SetBranchAddress("p_dedx_tof", &p_dedx_tof, &b_p_dedx_tof);
   fChain->SetBranchAddress("p_id", &p_id, &b_p_id);
   fChain->SetBranchAddress("p_isOffVertexClust", &p_isOffVertexClust, &b_p_isOffVertexClust);
   fChain->SetBranchAddress("p_isPrimaryVertex", &p_isPrimaryVertex, &b_p_isPrimaryVertex);
   fChain->SetBranchAddress("p_isUsedVertex", &p_isUsedVertex, &b_p_isUsedVertex);
   fChain->SetBranchAddress("p_isring", &p_isring, &b_p_isring);
   fChain->SetBranchAddress("p_isringmdc", &p_isringmdc, &b_p_isringmdc);
   fChain->SetBranchAddress("p_isringnomatch", &p_isringnomatch, &b_p_isringnomatch);
   fChain->SetBranchAddress("p_isringtrack", &p_isringtrack, &b_p_isringtrack);
   fChain->SetBranchAddress("p_kIsLepton", &p_kIsLepton, &b_p_kIsLepton);
   fChain->SetBranchAddress("p_kIsUsed", &p_kIsUsed, &b_p_kIsUsed);
   fChain->SetBranchAddress("p_mdcinnerchi2", &p_mdcinnerchi2, &b_p_mdcinnerchi2);
   fChain->SetBranchAddress("p_mdcouterchi2", &p_mdcouterchi2, &b_p_mdcouterchi2);
   fChain->SetBranchAddress("p_oa_hadr", &p_oa_hadr, &b_p_oa_hadr);
   fChain->SetBranchAddress("p_oa_lept", &p_oa_lept, &b_p_oa_lept);
   fChain->SetBranchAddress("p_p", &p_p, &b_p_p);
   fChain->SetBranchAddress("p_p_corr_em", &p_p_corr_em, &b_p_p_corr_em);
   fChain->SetBranchAddress("p_p_corr_ep", &p_p_corr_ep, &b_p_p_corr_ep);
   fChain->SetBranchAddress("p_p_corr_p", &p_p_corr_p, &b_p_p_corr_p);
   fChain->SetBranchAddress("p_p_corr_pim", &p_p_corr_pim, &b_p_p_corr_pim);
   fChain->SetBranchAddress("p_p_corr_pip", &p_p_corr_pip, &b_p_p_corr_pip);
   fChain->SetBranchAddress("p_phi", &p_phi, &b_p_phi);
   fChain->SetBranchAddress("p_phi_rich", &p_phi_rich, &b_p_phi_rich);
   fChain->SetBranchAddress("p_pid", &p_pid, &b_p_pid);
   fChain->SetBranchAddress("p_q", &p_q, &b_p_q);
   fChain->SetBranchAddress("p_r", &p_r, &b_p_r);
   fChain->SetBranchAddress("p_resolution", &p_resolution, &b_p_resolution);
   fChain->SetBranchAddress("p_resoultion", &p_resoultion, &b_p_resoultion);
   fChain->SetBranchAddress("p_rich_amp", &p_rich_amp, &b_p_rich_amp);
   fChain->SetBranchAddress("p_rich_centr", &p_rich_centr, &b_p_rich_centr);
   fChain->SetBranchAddress("p_rich_houtra", &p_rich_houtra, &b_p_rich_houtra);
   fChain->SetBranchAddress("p_rich_padnum", &p_rich_padnum, &b_p_rich_padnum);
   fChain->SetBranchAddress("p_rich_patmat", &p_rich_patmat, &b_p_rich_patmat);
   fChain->SetBranchAddress("p_rkchi2", &p_rkchi2, &b_p_rkchi2);
   fChain->SetBranchAddress("p_sector", &p_sector, &b_p_sector);
   fChain->SetBranchAddress("p_shw_sum0", &p_shw_sum0, &b_p_shw_sum0);
   fChain->SetBranchAddress("p_shw_sum1", &p_shw_sum1, &b_p_shw_sum1);
   fChain->SetBranchAddress("p_shw_sum2", &p_shw_sum2, &b_p_shw_sum2);
   fChain->SetBranchAddress("p_system", &p_system, &b_p_system);
   fChain->SetBranchAddress("p_theta", &p_theta, &b_p_theta);
   fChain->SetBranchAddress("p_theta_rich", &p_theta_rich, &b_p_theta_rich);
   fChain->SetBranchAddress("p_tof_mom", &p_tof_mom, &b_p_tof_mom);
   fChain->SetBranchAddress("p_tof_new", &p_tof_new, &b_p_tof_new);
   fChain->SetBranchAddress("p_tof_rec", &p_tof_rec, &b_p_tof_rec);
   fChain->SetBranchAddress("p_track_length", &p_track_length, &b_p_track_length);
   fChain->SetBranchAddress("p_tracklength", &p_tracklength, &b_p_tracklength);
   fChain->SetBranchAddress("p_z", &p_z, &b_p_z);
   fChain->SetBranchAddress("pim_beta", &pim_beta, &b_pim_beta);
   fChain->SetBranchAddress("pim_beta_new", &pim_beta_new, &b_pim_beta_new);
   fChain->SetBranchAddress("pim_dedx_mdc", &pim_dedx_mdc, &b_pim_dedx_mdc);
   fChain->SetBranchAddress("pim_dedx_tof", &pim_dedx_tof, &b_pim_dedx_tof);
   fChain->SetBranchAddress("pim_id", &pim_id, &b_pim_id);
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
   fChain->SetBranchAddress("pim_system", &pim_system, &b_pim_system);
   fChain->SetBranchAddress("pim_theta", &pim_theta, &b_pim_theta);
   fChain->SetBranchAddress("pim_theta_rich", &pim_theta_rich, &b_pim_theta_rich);
   fChain->SetBranchAddress("pim_tof_mom", &pim_tof_mom, &b_pim_tof_mom);
   fChain->SetBranchAddress("pim_tof_new", &pim_tof_new, &b_pim_tof_new);
   fChain->SetBranchAddress("pim_tof_rec", &pim_tof_rec, &b_pim_tof_rec);
   fChain->SetBranchAddress("pim_track_length", &pim_track_length, &b_pim_track_length);
   fChain->SetBranchAddress("pim_tracklength", &pim_tracklength, &b_pim_tracklength);
   fChain->SetBranchAddress("pim_z", &pim_z, &b_pim_z);
   fChain->SetBranchAddress("pt_costh_1", &pt_costh_1, &b_pt_costh_1);
   fChain->SetBranchAddress("pt_costh_2", &pt_costh_2, &b_pt_costh_2);
   fChain->SetBranchAddress("pt_costh_3", &pt_costh_3, &b_pt_costh_3);
   fChain->SetBranchAddress("pt_costh_4", &pt_costh_4, &b_pt_costh_4);
   fChain->SetBranchAddress("pt_counter", &pt_counter, &b_pt_counter);
   fChain->SetBranchAddress("pt_match_1", &pt_match_1, &b_pt_match_1);
   fChain->SetBranchAddress("pt_match_2", &pt_match_2, &b_pt_match_2);
   fChain->SetBranchAddress("pt_match_3", &pt_match_3, &b_pt_match_3);
   fChain->SetBranchAddress("pt_match_4", &pt_match_4, &b_pt_match_4);
   fChain->SetBranchAddress("pt_mult", &pt_mult, &b_pt_mult);
   fChain->SetBranchAddress("pt_p_1", &pt_p_1, &b_pt_p_1);
   fChain->SetBranchAddress("pt_p_2", &pt_p_2, &b_pt_p_2);
   fChain->SetBranchAddress("pt_p_3", &pt_p_3, &b_pt_p_3);
   fChain->SetBranchAddress("pt_p_4", &pt_p_4, &b_pt_p_4);
   fChain->SetBranchAddress("pt_phi0_1", &pt_phi0_1, &b_pt_phi0_1);
   fChain->SetBranchAddress("pt_phi0_2", &pt_phi0_2, &b_pt_phi0_2);
   fChain->SetBranchAddress("pt_phi0_3", &pt_phi0_3, &b_pt_phi0_3);
   fChain->SetBranchAddress("pt_phi0_4", &pt_phi0_4, &b_pt_phi0_4);
   fChain->SetBranchAddress("pt_phi_1", &pt_phi_1, &b_pt_phi_1);
   fChain->SetBranchAddress("pt_phi_2", &pt_phi_2, &b_pt_phi_2);
   fChain->SetBranchAddress("pt_phi_3", &pt_phi_3, &b_pt_phi_3);
   fChain->SetBranchAddress("pt_phi_4", &pt_phi_4, &b_pt_phi_4);
   fChain->SetBranchAddress("pt_px_1", &pt_px_1, &b_pt_px_1);
   fChain->SetBranchAddress("pt_px_2", &pt_px_2, &b_pt_px_2);
   fChain->SetBranchAddress("pt_px_3", &pt_px_3, &b_pt_px_3);
   fChain->SetBranchAddress("pt_px_4", &pt_px_4, &b_pt_px_4);
   fChain->SetBranchAddress("pt_py_1", &pt_py_1, &b_pt_py_1);
   fChain->SetBranchAddress("pt_py_2", &pt_py_2, &b_pt_py_2);
   fChain->SetBranchAddress("pt_py_3", &pt_py_3, &b_pt_py_3);
   fChain->SetBranchAddress("pt_py_4", &pt_py_4, &b_pt_py_4);
   fChain->SetBranchAddress("pt_pz_1", &pt_pz_1, &b_pt_pz_1);
   fChain->SetBranchAddress("pt_pz_2", &pt_pz_2, &b_pt_pz_2);
   fChain->SetBranchAddress("pt_pz_3", &pt_pz_3, &b_pt_pz_3);
   fChain->SetBranchAddress("pt_pz_4", &pt_pz_4, &b_pt_pz_4);
   fChain->SetBranchAddress("pt_theta0_1", &pt_theta0_1, &b_pt_theta0_1);
   fChain->SetBranchAddress("pt_theta0_2", &pt_theta0_2, &b_pt_theta0_2);
   fChain->SetBranchAddress("pt_theta0_3", &pt_theta0_3, &b_pt_theta0_3);
   fChain->SetBranchAddress("pt_theta0_4", &pt_theta0_4, &b_pt_theta0_4);
   fChain->SetBranchAddress("pt_theta_1", &pt_theta_1, &b_pt_theta_1);
   fChain->SetBranchAddress("pt_theta_2", &pt_theta_2, &b_pt_theta_2);
   fChain->SetBranchAddress("pt_theta_3", &pt_theta_3, &b_pt_theta_3);
   fChain->SetBranchAddress("pt_theta_4", &pt_theta_4, &b_pt_theta_4);
   fChain->SetBranchAddress("pt_x1_1", &pt_x1_1, &b_pt_x1_1);
   fChain->SetBranchAddress("pt_x1_2", &pt_x1_2, &b_pt_x1_2);
   fChain->SetBranchAddress("pt_x1_3", &pt_x1_3, &b_pt_x1_3);
   fChain->SetBranchAddress("pt_x1_4", &pt_x1_4, &b_pt_x1_4);
   fChain->SetBranchAddress("pt_x2_1", &pt_x2_1, &b_pt_x2_1);
   fChain->SetBranchAddress("pt_x2_2", &pt_x2_2, &b_pt_x2_2);
   fChain->SetBranchAddress("pt_x2_3", &pt_x2_3, &b_pt_x2_3);
   fChain->SetBranchAddress("pt_x2_4", &pt_x2_4, &b_pt_x2_4);
   fChain->SetBranchAddress("pt_xtarg_1", &pt_xtarg_1, &b_pt_xtarg_1);
   fChain->SetBranchAddress("pt_xtarg_2", &pt_xtarg_2, &b_pt_xtarg_2);
   fChain->SetBranchAddress("pt_xtarg_3", &pt_xtarg_3, &b_pt_xtarg_3);
   fChain->SetBranchAddress("pt_xtarg_4", &pt_xtarg_4, &b_pt_xtarg_4);
   fChain->SetBranchAddress("pt_y1_1", &pt_y1_1, &b_pt_y1_1);
   fChain->SetBranchAddress("pt_y1_2", &pt_y1_2, &b_pt_y1_2);
   fChain->SetBranchAddress("pt_y1_3", &pt_y1_3, &b_pt_y1_3);
   fChain->SetBranchAddress("pt_y1_4", &pt_y1_4, &b_pt_y1_4);
   fChain->SetBranchAddress("pt_y2_1", &pt_y2_1, &b_pt_y2_1);
   fChain->SetBranchAddress("pt_y2_2", &pt_y2_2, &b_pt_y2_2);
   fChain->SetBranchAddress("pt_y2_3", &pt_y2_3, &b_pt_y2_3);
   fChain->SetBranchAddress("pt_y2_4", &pt_y2_4, &b_pt_y2_4);
   fChain->SetBranchAddress("pt_ytarg_1", &pt_ytarg_1, &b_pt_ytarg_1);
   fChain->SetBranchAddress("pt_ytarg_2", &pt_ytarg_2, &b_pt_ytarg_2);
   fChain->SetBranchAddress("pt_ytarg_3", &pt_ytarg_3, &b_pt_ytarg_3);
   fChain->SetBranchAddress("pt_ytarg_4", &pt_ytarg_4, &b_pt_ytarg_4);
   fChain->SetBranchAddress("runnumber", &runnumber, &b_runnumber);
   fChain->SetBranchAddress("start_cluster_size_1st", &start_cluster_size_1st, &b_start_cluster_size_1st);
   fChain->SetBranchAddress("start_cluster_size_2nd", &start_cluster_size_2nd, &b_start_cluster_size_2nd);
   fChain->SetBranchAddress("start_cluster_size_max", &start_cluster_size_max, &b_start_cluster_size_max);
   fChain->SetBranchAddress("start_corrflag", &start_corrflag, &b_start_corrflag);
   fChain->SetBranchAddress("start_counter", &start_counter, &b_start_counter);
   fChain->SetBranchAddress("start_flag", &start_flag, &b_start_flag);
   fChain->SetBranchAddress("start_module", &start_module, &b_start_module);
   fChain->SetBranchAddress("start_mult", &start_mult, &b_start_mult);
   fChain->SetBranchAddress("start_multipliticy", &start_multipliticy, &b_start_multipliticy);
   fChain->SetBranchAddress("start_resolution", &start_resolution, &b_start_resolution);
   fChain->SetBranchAddress("start_strip", &start_strip, &b_start_strip);
   fChain->SetBranchAddress("start_time", &start_time, &b_start_time);
   fChain->SetBranchAddress("start_time_2nd", &start_time_2nd, &b_start_time_2nd);
   fChain->SetBranchAddress("start_track", &start_track, &b_start_track);
   fChain->SetBranchAddress("start_width", &start_width, &b_start_width);
   fChain->SetBranchAddress("totalmult", &totalmult, &b_totalmult);
   fChain->SetBranchAddress("trigbit", &trigbit, &b_trigbit);
   fChain->SetBranchAddress("trigdec", &trigdec, &b_trigdec);
   fChain->SetBranchAddress("trigdownscale", &trigdownscale, &b_trigdownscale);
   fChain->SetBranchAddress("trigdownscaleflag", &trigdownscaleflag, &b_trigdownscaleflag);
   Notify();
}

// ---------------------------------------------------------------------------------------------------------
Bool_t PPim::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

// ---------------------------------------------------------------------------------------------------------
void PPim::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}

// ---------------------------------------------------------------------------------------------------------
Int_t PPim::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
