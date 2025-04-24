/*////////////////////////////////////////////////////
Original code by J. Estee, MIT
Code updated for comparisons by J. L. Barrow, MIT/TAU; further modified by Caleb Fogler, ODU;
////////////////////////////////////////////////////*/

#include <cstdlib>
#include <iostream>
#include <chrono>
#include <TFile.h>
#include <TTree.h>
#include <THStack.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TF1.h>
#include <TH1.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TBenchmark.h>
#include "clas12reader.h"
#include "clas12writer.h"
#include "skimmer/DC_pion_fiducial.h"
#include "skimmer/DC_pion_fiducial.cpp"
#include "HipoChain.h"
#include <string>
#include <TH2.h>
#include <TStyle.h>
Double_t DeltaThetaPiQBin(Double_t ThetaMin, Double_t ThetaMax);
Double_t ExpFunc(Double_t A, Double_t B, Double_t C, Double_t D, Double_t E, Double_t P);
void EventBinner(bool QorNot, Int_t ichrg, Double_t P_pi, Double_t theta_piq, Double_t q2_epi, Double_t W_epi, Double_t Wmin, Double_t DeltaW, Int_t Wbini, Int_t &iPmom,
                 Int_t &iThPiQ, Int_t &iQ2, Int_t &iW);
Int_t eCutter(Double_t VWcut, Double_t Lv, Double_t Lw, Double_t vz_e, Double_t vz_e_min, Double_t vz_e_max, Double_t Vperp_e, Double_t Vperp_ecut);
Int_t PiCutter(Int_t iCharge, Double_t P_pi, Double_t p5, Double_t p10, Double_t Vperp_pi, Double_t Vdiff_pi,
              Double_t Vperp_pipcut[3], Double_t Vperp_pip_ms[3][2], Double_t Vperp_pimcut[2], Double_t Vperp_pim_ms[2][2], Double_t Nsig = 3);
std::ifstream& GotoLine(std::ifstream& file, unsigned int num);
using namespace clas12;
using namespace std;

void SetLorentzVector(TLorentzVector &p4,clas12::region_part_ptr rp)
{
  p4.SetXYZM(rp->par()->getPx(),rp->par()->getPy(),
	      rp->par()->getPz(),p4.M());

}

//MC comparison added by J. L. Barrow (JLB)
void e4nu_radComparison_4p2GeV(TString inFile_data = "", TString inFile_MC = "", TString inFile_opg = "", TString inFile_opgm = "",
                               TString inFile_tpd = "", TString inFile_opgn = "", TString inFile_opgmn = "",
                               TString outputFile = "", Double_t beamE = 0., bool bQorNot = true, Int_t iGENmod = 0,
                               bool bAllEvents = true, bool bMultiPion = false, bool bRAD = true, bool bAC = true,
                               bool bData = true, bool bOPG = true)
{
  if(bQorNot==true){cout << "Binning in theta_piQ" << endl;}
  else             {cout << "Binning in theta_pi (i.e., NOT Q)" << endl;}
  bool bending = true;
  // Record start time
  auto start = std::chrono::high_resolution_clock::now();
  //Changed by JLB
  clas12root::HipoChain chain_data; clas12root::HipoChain chain_MC; clas12root::HipoChain chain_opg; clas12root::HipoChain chain_opgm;
  clas12root::HipoChain chain_opgn; clas12root::HipoChain chain_opgmn; clas12root::HipoChain chain_MCrad;
  //JLB
  TString path = "/work/clas12/cfogler/e4nu/skimmer/", path_MC = "/cache/clas12/rg-b/MC/e4nu/cfogler/", path_cut = "/volatile/clas12/cfogler/";
  chain_data.Add(path_MC+inFile_data); chain_MC.Add(path_MC+inFile_MC); chain_opg.Add(path_MC+inFile_opg); chain_opgm.Add(path_MC+inFile_opgm);
  chain_opgn.Add(path_MC+inFile_opgn); chain_opgmn.Add(path_cut+inFile_opgmn); chain_MCrad.Add(path_cut+inFile_tpd);
  //JLB (create clas12reader with just tag 0 events)
  chain_data.SetReaderTags({0}); chain_MC.SetReaderTags({0}); chain_opg.SetReaderTags({0}); chain_opgm.SetReaderTags({0});
  chain_opgn.SetReaderTags({0}); chain_opgmn.SetReaderTags({0}); chain_MCrad.SetReaderTags({0});
  auto config_c12_data=chain_data.GetC12Reader(); auto config_c12_MC=chain_MC.GetC12Reader(); auto config_c12_opg=chain_opg.GetC12Reader(); auto config_c12_opgm=chain_opgm.GetC12Reader();
  auto config_c12_opgn=chain_opgn.GetC12Reader(); auto config_c12_opgmn=chain_opgmn.GetC12Reader(); auto config_c12_MCrad=chain_MCrad.GetC12Reader();
  //now get reference to (unique)ptr for accessing data in loop.  This will point to the correct place when file changes
  auto& c12_data=chain_data.C12ref(); auto& c12_MC=chain_MC.C12ref(); auto& c12_opg=chain_opg.C12ref(); auto& c12_opgm=chain_opgm.C12ref();
  auto& c12_opgn=chain_opgn.C12ref(); auto& c12_opgmn=chain_opgmn.C12ref(); auto& c12_gen=chain_MC.C12ref(); auto& c12_MCrad=chain_MCrad.C12ref();
  //struct with all relevent rcdb values
  //auto& rcdbData= config_c12_data->rcdb()->current(); auto& rcdbMC= config_c12_MC->rcdb()->current();
  Double_t weight_data = 1, weight_MC = 1, weight_opg = 1, weight_opgm = 1, weight_opgn = 1, weight_opgmn = 1;

  auto db=TDatabasePDG::Instance();
  chain_data.db()->turnOffQADB(); chain_MC.db()->turnOffQADB(); chain_opg.db()->turnOffQADB();
  chain_opgm.db()->turnOffQADB(); chain_opgn.db()->turnOffQADB(); chain_opgmn.db()->turnOffQADB(); chain_MCrad.db()->turnOffQADB();
  Double_t mass_e = db->GetParticle(11)->Mass();
  Double_t mass_p = db->GetParticle(2212)->Mass();
  Double_t mass_n = db->GetParticle(2112)->Mass();
  Double_t mass_pip = db->GetParticle(211)->Mass();
  Double_t mass_pim = db->GetParticle(-211)->Mass();
  Double_t mass_N = (mass_p + mass_n) / 2.;

  //Make some particles
  TLorentzVector v4_beam_data(0,0,beamE, beamE);
  TLorentzVector v4_beam_MC(0,0,beamE, beamE);   TLorentzVector v4_beam_MCrad(0,0,beamE, beamE);
  TLorentzVector v4_beam_opg(0,0,beamE, beamE);  TLorentzVector v4_beam_opgn(0,0,beamE, beamE);
  TLorentzVector v4_beam_opgm(0,0,beamE, beamE); TLorentzVector v4_beam_opgmn(0,0,beamE, beamE);
  TLorentzVector v4_el(0,0,0,db->GetParticle(11)->Mass());
  TLorentzVector v4_pr(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector v4_nr(0,0,0,db->GetParticle(2112)->Mass());
  TLorentzVector v4_pip(0,0,0,db->GetParticle(211)->Mass());
  TLorentzVector v4_pim(0,0,0,db->GetParticle(-211)->Mass());

  DCPionFiducial dcFid;

  gBenchmark->Start("timer");

  ///Bin max values
  Double_t q2histomax = 3.2, Whistomax = 3., Wthetahistomax = 60., omegahistomax = 5., ThetaQhistomax = 80.;

  //JLB: Stacked versions of 1D histograms (show with "nostack" option)
  //These stacked histograms will store the 1D histograms as they are made
  THStack *q2_stacked = new THStack("q2_stacked","Q^{2} Data to GENIE Comparison;Q^{2};Counts");
  THStack *p_chi2_stacked = new THStack("p_chi2_stacked","Proton #chi^{2} Data to GENIE Comparison;Proton #chi^{2};Counts");
  THStack *el_vz_h_stacked = new THStack("el_vz_h_stacked","Z-Vertex Data to GENIE Comparison;Z-vertex (cm);Counts");
  THStack *W_stacked = new THStack("W_stacked","Invariant Mass, W, Data to GENIE Comparison;W (GeV/c^{2});Counts");
  THStack *omega_stacked = new THStack("omega_stacked","Energy Transfer Data to GENIE Comparison;#omega (GeV);Counts");
  THStack *hs_q2 = new THStack("hs_q2",";Q^{2} (GeV^{2});");
  THStack *hs_Ppi = new THStack("hs_Ppi",";P_{pi} (GeV);");
  THStack *hs_theta_piq = new THStack("hs_theta_piq",";#theta_{piq} (deg);");
  THStack *hs_theta_pi = new THStack("hs_theta_pi",";#theta_{pi} (deg);");
  THStack *hs_q2_data = new THStack("hs_q2_data","Q^{2} Data to GENIE Comparison;Q^{2};Counts");
  THStack *hs_W_data = new THStack("hs_W_data","Invariant Mass, W, Data to GENIE Comparison;W (GeV/c^{2});Counts");
  THStack *hs_Ppi_data = new THStack("hs_Ppi_data",";P_{pi} (GeV);");
  THStack *hs_theta_piq_data = new THStack("hs_theta_piq_data",";#theta_{piq} (deg);");
  THStack *hs_theta_pi_data = new THStack("hs_theta_pi_data",";#theta_{pi} (deg);");

  //Stacked 1D histograms for data comparisons by channel
  THStack *q2_byChannel_stacked = new THStack("q2_byChannel_stacked","Data: Q^{2} Distributions by Reconstructed Channel;Q^{2};Counts");
  THStack *W_byChannel_stacked = new THStack("W_byChannel_stacked","Data: Invariant Mass, W, Distributions by Reconstructed Channel;W (GeV/c^{2});Counts");

  //Histograms for stacking by DC sector
  //W, Invariant mass
  THStack *W_bySector_data_stacked = new THStack("W_bySector_data_stacked","Data: Invariant Mass, W, Distributions by DC Sector;W (GeV/c^{2});Counts");
  TH1D *h_W_bySector_data = new TH1D("h_W_bySector_data","Data, All Sectors;W (GeV/c^{2};Counts",400,0.5,Whistomax);
  TH1D *h_W_S1_data = new TH1D("h_W_S1_data","Datax6, Sector #1;W (GeV/c^{2};Counts",400,0.5,Whistomax);
  TH1D *h_W_S2_data = new TH1D("h_W_S2_data","Datax6, Sector #2;W (GeV/c^{2};Counts",400,0.5,Whistomax);
  TH1D *h_W_S3_data = new TH1D("h_W_S3_data","Datax6, Sector #3;W (GeV/c^{2};Counts",400,0.5,Whistomax);
  TH1D *h_W_S4_data = new TH1D("h_W_S4_data","Datax6, Sector #4;W (GeV/c^{2};Counts",400,0.5,Whistomax);
  TH1D *h_W_S5_data = new TH1D("h_W_S5_data","Datax6, Sector #5;W (GeV/c^{2};Counts",400,0.5,Whistomax);
  TH1D *h_W_S6_data = new TH1D("h_W_S6_data","Datax6, Sector #6;W (GeV/c^{2};Counts",400,0.5,Whistomax);
  //Q^2
  THStack *q2_bySector_data_stacked = new THStack("q2_bySector_data_stacked","Data: Q^{2} Distributions by DC Sector;Q^{2};Counts");
  TH1D *h_q2_bySector_data = new TH1D("h_q2_bySector_data","Data, All Sectors;Q^{2};Counts",400,0.,q2histomax);
  TH1D *h_q2_S1_data = new TH1D("h_q2_S1_data","Datax6, Sector #1;Q^{2};Counts",400,0.,q2histomax);
  TH1D *h_q2_S2_data = new TH1D("h_q2_S2_data","Datax6, Sector #2;Q^{2};Counts",400,0.,q2histomax);
  TH1D *h_q2_S3_data = new TH1D("h_q2_S3_data","Datax6, Sector #3;Q^{2};Counts",400,0.,q2histomax);
  TH1D *h_q2_S4_data = new TH1D("h_q2_S4_data","Datax6, Sector #4;Q^{2};Counts",400,0.,q2histomax);
  TH1D *h_q2_S5_data = new TH1D("h_q2_S5_data","Datax6, Sector #5;Q^{2};Counts",400,0.,q2histomax);
  TH1D *h_q2_S6_data = new TH1D("h_q2_S6_data","Datax6, Sector #6;Q^{2};Counts",400,0.,q2histomax);
  //W vs. Theta
  TH2D *h2_theta_W_S1_data = new TH2D("h2_theta_W_S1_data","Datax6, Sector #1: W vs #theta;W (GeV/c^{2});#theta",400,0.,Whistomax,400,0.,Wthetahistomax);
  TH2D *h2_theta_W_S2_data = new TH2D("h2_theta_W_S2_data","Datax6, Sector #2: W vs #theta;W (GeV/c^{2});#theta",400,0.,Whistomax,400,0.,Wthetahistomax);
  TH2D *h2_theta_W_S3_data = new TH2D("h2_theta_W_S3_data","Datax6, Sector #3: W vs #theta;W (GeV/c^{2});#theta",400,0.,Whistomax,400,0.,Wthetahistomax);
  TH2D *h2_theta_W_S4_data = new TH2D("h2_theta_W_S4_data","Datax6, Sector #4: W vs #theta;W (GeV/c^{2});#theta",400,0.,Whistomax,400,0.,Wthetahistomax);
  TH2D *h2_theta_W_S5_data = new TH2D("h2_theta_W_S5_data","Datax6, Sector #5: W vs #theta;W (GeV/c^{2});#theta",400,0.,Whistomax,400,0.,Wthetahistomax);
  TH2D *h2_theta_W_S6_data = new TH2D("h2_theta_W_S6_data","Datax6, Sector #6: W vs #theta;W (GeV/c^{2});#theta",400,0.,Whistomax,400,0.,Wthetahistomax);
  //Now for GENIE simulation
  //W, Invariant mass
  THStack *W_bySector_MC_stacked = new THStack("W_bySector_MC_stacked","GENIE: Invariant Mass, W, Distributions by DC Sector;W (GeV/c^{2});Counts");
  TH1D *W_bySector_MC = new TH1D("W_bySector_MC","GENIE, All Sectors;W (GeV/c^{2};Counts",400,0.5,Whistomax);
  TH1D *h_W_S1_MC = new TH1D("h_W_S1_MC","GENIEx6, Sector #1;W (GeV/c^{2};Counts",400,0.5,Whistomax);
  TH1D *h_W_S2_MC = new TH1D("h_W_S2_MC","GENIEx6, Sector #2;W (GeV/c^{2};Counts",400,0.5,Whistomax);
  TH1D *h_W_S3_MC = new TH1D("h_W_S3_MC","GENIEx6, Sector #3;W (GeV/c^{2};Counts",400,0.5,Whistomax);
  TH1D *h_W_S4_MC = new TH1D("h_W_S4_MC","GENIEx6, Sector #4;W (GeV/c^{2};Counts",400,0.5,Whistomax);
  TH1D *h_W_S5_MC = new TH1D("h_W_S5_MC","GENIEx6, Sector #5;W (GeV/c^{2};Counts",400,0.5,Whistomax);
  TH1D *h_W_S6_MC = new TH1D("h_W_S6_MC","GENIEx6, Sector #6;W (GeV/c^{2};Counts",400,0.5,Whistomax);
  //Q^2
  THStack *q2_bySector_MC_stacked = new THStack("q2_bySector_MC_stacked","GENIE: Q^{2} Distributions by DC Sector;Q^{2};Counts");
  TH1D *h_q2_bySector_MC = new TH1D("h_q2_bySector_MC","GENIE, All Sectors;Q^{2};Counts",400,0.,q2histomax);
  TH1D *h_q2_S1_MC = new TH1D("h_q2_S1_MC","GENIEx6, Sector #1;Q^{2};Counts",400,0.,q2histomax);
  TH1D *h_q2_S2_MC = new TH1D("h_q2_S2_MC","GENIEx6, Sector #2;Q^{2};Counts",400,0.,q2histomax);
  TH1D *h_q2_S3_MC = new TH1D("h_q2_S3_MC","GENIEx6, Sector #3;Q^{2};Counts",400,0.,q2histomax);
  TH1D *h_q2_S4_MC = new TH1D("h_q2_S4_MC","GENIEx6, Sector #4;Q^{2};Counts",400,0.,q2histomax);
  TH1D *h_q2_S5_MC = new TH1D("h_q2_S5_MC","GENIEx6, Sector #5;Q^{2};Counts",400,0.,q2histomax);
  TH1D *h_q2_S6_MC = new TH1D("h_q2_S6_MC","GENIEx6, Sector #6;Q^{2};Counts",400,0.,q2histomax);
  //W vs. Theta
  TH2D *h2_theta_W_S1_MC = new TH2D("h2_theta_W_S1_MC","GENIEx6, Sector #1: W vs #theta;W (GeV/c^{2});#theta",400,0.,Whistomax,400,0.,Wthetahistomax);
  TH2D *h2_theta_W_S2_MC = new TH2D("h2_theta_W_S2_MC","GENIEx6, Sector #2: W vs #theta;W (GeV/c^{2});#theta",400,0.,Whistomax,400,0.,Wthetahistomax);
  TH2D *h2_theta_W_S3_MC = new TH2D("h2_theta_W_S3_MC","GENIEx6, Sector #3: W vs #theta;W (GeV/c^{2});#theta",400,0.,Whistomax,400,0.,Wthetahistomax);
  TH2D *h2_theta_W_S4_MC = new TH2D("h2_theta_W_S4_MC","GENIEx6, Sector #4: W vs #theta;W (GeV/c^{2});#theta",400,0.,Whistomax,400,0.,Wthetahistomax);
  TH2D *h2_theta_W_S5_MC = new TH2D("h2_theta_W_S5_MC","GENIEx6, Sector #5: W vs #theta;W (GeV/c^{2});#theta",400,0.,Whistomax,400,0.,Wthetahistomax);
  TH2D *h2_theta_W_S6_MC = new TH2D("h2_theta_W_S6_MC","GENIEx6, Sector #6: W vs #theta;W (GeV/c^{2});#theta",400,0.,Whistomax,400,0.,Wthetahistomax);

  //JLB: 1D histograms for data
  //Data
  TH1D *h_q2_data = new TH1D("h_q2_data",";Q^{2};Counts",400,0.,q2histomax);
  TH1D *h_q2pip_data = new TH1D("h_q2pip_data",";Q^{2};Counts",400,0.,q2histomax);
  TH1D *h_q2pim_data = new TH1D("h_q2pim_data",";Q^{2};Counts",400,0.,q2histomax);
  TH1D *h_Ppi_data = new TH1D("h_Ppi_data",";;",400,0.,3.2);
  TH1D *h_Ppip_data = new TH1D("h_Ppip_data",";;",400,0.,3.2);
  TH1D *h_Ppim_data = new TH1D("h_Ppim_data",";;",400,0.,3.2);
  TH1D *h_theta_piq_data = new TH1D("h_theta_piq_data",";;",600,0.,85.);
  TH1D *h_theta_pipq_data = new TH1D("h_theta_pipq_data",";;",600,0.,85.);
  TH1D *h_theta_pimq_data = new TH1D("h_theta_pimq_data",";;",600,0.,85.);
  TH1D *h_theta_pi_data = new TH1D("h_theta_pi_data",";;",600,0.,85.);
  TH1D *h_theta_pip_data = new TH1D("h_theta_pip_data",";;",600,0.,85.);
  TH1D *h_theta_pim_data = new TH1D("h_theta_pim_data",";;",600,0.,85.);
  TH1D *h_energy_transfer_0_data = new TH1D("h_energy_transfer_0_data","Data;#omega (GeV);Counts",1000,0.,omegahistomax);
  TH1D *h_energy_transfer_1_data = new TH1D("h_energy_transfer_1_data","Data;#omega (GeV);Counts",1000,0.,omegahistomax);
  TH1D *h_energy_transfer_2_data = new TH1D("h_energy_transfer_2_data","Data;#omega (GeV);Counts",1000,0.,omegahistomax);
  TH1D *h_energy_transfer_3_data = new TH1D("h_energy_transfer_3_data","Data;#omega (GeV);Counts",1000,0.,omegahistomax);
  TH1D *h_el_vz_h_data = new TH1D("h_el_vz_h_data","Data;Z-vertex (cm);Counts",1000,-10,10);
  TH1D *h_W_data = new TH1D("h_W_data",";W (GeV/c^{2};Counts",400,0.5,Whistomax);
  TH1D *h_Wp_data   = new TH1D("h_Wp_data",";W (GeV/c^{2});",400,0.5,Whistomax);
  TH1D *h_Wpip_data = new TH1D("h_Wpip_data",";W (GeV/c^{2};Counts",400,0.5,Whistomax);
  TH1D *h_Wpim_data = new TH1D("h_Wpim_data",";W (GeV/c^{2};Counts",400,0.5,Whistomax);
  TH1D *h_omega_data = new TH1D("h_omega_data","Data;#omega (GeV);Counts",1000,0.,omegahistomax);
  TH1D *h_p_proton_data = new TH1D("h_p_proton_data","Data;Momentum (GeV/c);Counts",1000,0.,4.);
  TH1D *h_Wpip_data_chi2 = new TH1D("h_Wpip_data_chi2",";W (GeV/c^{2};",400,0.5,Whistomax);
  TH1D *h_Wpip_data_chi3 = new TH1D("h_Wpip_data_chi3",";W (GeV/c^{2};",400,0.5,Whistomax);
  TH1D *h_Wpip_data_ray  = new TH1D("h_Wpip_data_ray",";W (GeV/c^{2};",400,0.5,Whistomax);
  TH1D *h_Betapip_data_ray = new TH1D("h_Betapip_data_ray",";beta (c);",400,0.95,1.02);

  //JLB: 1D histograms for GENIE
  //MC Simulation
  TH1D *h_q2_MC = new TH1D("h_q2_MC","GENIE;Q^{2};Counts",400,0.,q2histomax);
  TH1D *h_q2_qe_MC = new TH1D("h_q2_qe_MC","GENIE QE;Q^{2};Counts",400,0.,q2histomax);
  TH1D *h_q2_mec_MC = new TH1D("h_q2_mec_MC","GENIE MEC;Q^{2};Counts",400,0.,q2histomax);
  TH1D *h_q2_res_MC = new TH1D("h_q2_res_MC","GENIE RES;Q^{2};Counts",400,0.,q2histomax);
  TH1D *h_q2_dis_MC = new TH1D("h_q2_dis_MC","GENIE DIS;Q^{2};Counts",400,0.,q2histomax);
  TH1D *h_Ppi_MC = new TH1D("h_Ppi_MC","GENIE;;",400,0.,4.);
  TH1D *h_theta_piq_MC = new TH1D("h_theta_piq_MC","GENIE;;",600,0.,60.);
  TH1D *h_theta_pi_MC = new TH1D("h_theta_pi_MC","GENIE;;",600,0.,85.);
  TH1D *h_energy_transfer_0_MC = new TH1D("h_energy_transfer_0_MC","GENIE;#omega (GeV);Counts",1000,0.,omegahistomax);
  TH1D *h_energy_transfer_1_MC = new TH1D("h_energy_transfer_1_MC","GENIE;#omega (GeV);Counts",1000,0.,omegahistomax);
  TH1D *h_energy_transfer_2_MC = new TH1D("h_energy_transfer_2_MC","GENIE;#omega (GeV);Counts",1000,0.,omegahistomax);
  TH1D *h_energy_transfer_3_MC = new TH1D("h_energy_transfer_3_MC","GENIE;#omega (GeV);Counts",1000,0.,omegahistomax);
  TH1D *h_el_vz_h_MC = new TH1D("h_el_vz_h_MC","GENIE;Counts",1000,-10,10);
  TH1D *h_W_MC = new TH1D("h_W_MC","GENIE;W (GeV/c^{2};Counts",400,0.5,Whistomax);
  TH1D *h_W_qe_MC = new TH1D("h_W_qe_MC","GENIE QE;W (GeV/c^{2};Counts",400,0.5,Whistomax);
  TH1D *h_W_mec_MC = new TH1D("h_W_mec_MC","GENIE MEC;W (GeV/c^{2};Counts",400,0.5,Whistomax);
  TH1D *h_W_res_MC = new TH1D("h_W_res_MC","GENIE RES;W (GeV/c^{2};Counts",400,0.5,Whistomax);
  TH1D *h_W_dis_MC = new TH1D("h_W_dis_MC","GENIE DIS;W (GeV/c^{2}",400,0.5,Whistomax);
  TH1D *h_omega_MC = new TH1D("h_omega_MC","GENIE;#omega (GeV);Counts",1000,0.,omegahistomax);
  TH1D *h_omega_qe_MC = new TH1D("h_omega_qe_MC","GENIE QE;#omega (GeV);Counts",1000,0.,omegahistomax);
  TH1D *h_omega_mec_MC = new TH1D("h_omega_mec_MC","GENIE MEC;#omega (GeV);Counts",1000,0.,omegahistomax);
  TH1D *h_omega_res_MC = new TH1D("h_omega_res_MC","GENIE RES;#omega (GeV);Counts",1000,0.,omegahistomax);
  TH1D *h_omega_dis_MC = new TH1D("h_omega_dis_MC","GENIE DIS;#omega (GeV);Counts",1000,0.,omegahistomax);
  TH1D *h_p_proton_MC = new TH1D("h_p_proton_MC","GENIE;Momentum (GeV/c);Counts",1000,0.,4.);
  Int_t resnum = 20;
  THStack *hs_W_resid_MC = new THStack("hs_W_resid_MC",";W (GeV/c^{2});");
  TH1D *h_W_resid_MC[resnum];
  //MCrad Simulation
  TH1D *h_q2_MCrad = new TH1D("h_q2_MCrad","GENIE rad;Q^{2};Counts",400,0.,q2histomax);
  TH1D *h_q2_qe_MCrad = new TH1D("h_q2_qe_MCrad","GENIE rad QE;Q^{2};Counts",400,0.,q2histomax);
  TH1D *h_q2_mec_MCrad = new TH1D("h_q2_mec_MCrad","GENIE rad MEC;Q^{2};Counts",400,0.,q2histomax);
  TH1D *h_q2_res_MCrad = new TH1D("h_q2_res_MCrad","GENIE rad RES;Q^{2};Counts",400,0.,q2histomax);
  TH1D *h_q2_dis_MCrad = new TH1D("h_q2_dis_MCrad","GENIE rad DIS;Q^{2};Counts",400,0.,q2histomax);
  TH1D *h_Ppi_MCrad = new TH1D("h_Ppi_MCrad","GENIE rad;;",400,0.,4.);
  TH1D *h_theta_piq_MCrad = new TH1D("h_theta_piq_MCrad","GENIE rad;;",600,0.,60.);
  TH1D *h_theta_pi_MCrad = new TH1D("h_theta_pi_MCrad","GENIE rad;;",600,0.,85.);
  TH1D *h_W_MCrad = new TH1D("h_W_MCrad","GENIE rad;W (GeV/c^{2};Counts",400,0.5,Whistomax);
  TH1D *h_W_qe_MCrad = new TH1D("h_W_qe_MCrad","GENIE rad QE;W (GeV/c^{2};Counts",400,0.5,Whistomax);
  TH1D *h_W_mec_MCrad = new TH1D("h_W_mec_MCrad","GENIE rad MEC;W (GeV/c^{2};Counts",400,0.5,Whistomax);
  TH1D *h_W_res_MCrad = new TH1D("h_W_res_MCrad","GENIE rad RES;W (GeV/c^{2};Counts",400,0.5,Whistomax);
  TH1D *h_W_dis_MCrad = new TH1D("h_W_dis_MCrad","GENIE rad DIS;W (GeV/c^{2}",400,0.5,Whistomax);
  TH1D *h_weight_MCrad = new TH1D("h_weight_MCrad","GENIE Weight;Weight",400,0.,50.);
  TH1D *h_processid_MCrad = new TH1D("h_processid_MCrad","GENIE Weight;Weight",400,0.,50.);

  //2D histograms for specific data/MC samples
  //Data
  TH2D *h2_theta_q2_data = new TH2D("h2_theta_q2_data","Data: Electron Angle vs. Q^{2};Q^{2};#theta",400,0.,4.4,400,0.,Wthetahistomax);//1.2 to 4.4
  TH2D *h2_ecal_data = new TH2D("h2_ecal_data","Data: Electron Momentum Sampling Fraction;Electron Momentum (GeV/c);Sampling Fraction",400,0.,4.4,1000,0.,0.5);//1.2 to 4.4
  TH2D *h2_theta_mom_data = new TH2D("h2_theta_mom_data","Data: Electron Angle vs. Momentum;Momentum (GeV/c);#theta",400,0.,5.,400,0.,Wthetahistomax);
  TH2D *h2_q2_W_data = new TH2D("h2_q2_W_data","Data: Q^{2} vs. W;W (GeV/c^{2});Q^{2}",400,0.,Whistomax,400,0.,4.4);//1.2 to 4.4
  TH2D *h2_theta_W_data = new TH2D("h2_theta_W_data","Data: W vs #theta;W (GeV/c^{2});#theta",400,0.,Whistomax,400,0.,Wthetahistomax);
  TH2D *h2_theta_omega_data = new TH2D("h2_theta_omega_data","Data: Energy Transfer vs. #theta;#omega (GeV);#theta",1000,0.,omegahistomax,400,0.,Wthetahistomax);
  TH2D *h2_mult_p_pis_data = new TH2D("h2_mult_p_pis_data","Data: Multiplicity of Protons vs. Pions #pi^{#pm};Proton Multiplicity;#pi^{#pm} Multiplicity",5,-0.5,4.5,5,-0.5,4.5);
  Double_t Phistomax = 4.;
  TH2D *h2_q2_Wpip_data = new TH2D("h2_q2_Wpip_data","Data: Q^{2} vs. W (Pi+);W (GeV/c^{2});Q^{2}",320,1.,2.6,400,0.,4.4);//1.2 to 4.4
  TH2D *h2_q2_Wpim_data = new TH2D("h2_q2_Wpim_data","Data: Q^{2} vs. W (Pi-);W (GeV/c^{2});Q^{2}",320,1.,2.6,400,0.,4.4);//1.2 to 4.4
  TH2D *h2_Ppip_theta_data = new TH2D("h2_Ppip_theta_data","Data: Pi+ Momentum vs. Electron Angle;#theta_e;Momentum (GeV/c)",400,0.,Wthetahistomax,400,0.,Phistomax);
  TH2D *h2_Ppim_theta_data = new TH2D("h2_Ppim_theta_data","Data: Pi- Momentum vs. Electron Angle;#theta_e;Momentum (GeV/c)",400,0.,Wthetahistomax,400,0.,Phistomax);
  TH2D *h2_Ppip_thetapip_data =  new TH2D("h2_Ppip_thetapip_data","Data: Pi+ Momentum vs. Pi+ Angle;#theta_pip;Momentum (GeV/c)",400,0.,Wthetahistomax,400,0.,Phistomax);
  TH2D *h2_Ppim_thetapim_data =  new TH2D("h2_Ppim_thetapim_data","Data: Pi- Momentum vs. Pi- Angle;#theta_pim;Momentum (GeV/c)",400,0.,Wthetahistomax,400,0.,Phistomax);
  TH2D *h2_Ppip_thetapipq_data = new TH2D("h2_Ppip_thetapipq_data","Data: Pi+ Momentum vs. Pipq Angle;#theta_pipq;Momentum (GeV/c)",400,0.,Wthetahistomax,400,0.,Phistomax);
  TH2D *h2_Ppim_thetapimq_data = new TH2D("h2_Ppim_thetapimq_data","Data: Pi- Momentum vs. Pimq Angle;#theta_pimq;Momentum (GeV/c)",400,0.,Wthetahistomax,400,0.,Phistomax);
  //cut
  TH2D *h2_q2_Wpip_cut = new TH2D("h2_q2_Wpip_cut","Data cut: Q^{2} vs. W (Pi+);W (GeV/c^{2});Q^{2}",320,1.,2.6,400,0.,4.4);//1.2 to 4.4
  TH2D *h2_q2_Wpim_cut = new TH2D("h2_q2_Wpim_cut","Data cut: Q^{2} vs. W (Pi-);W (GeV/c^{2});Q^{2}",320,1.,2.6,400,0.,4.4);//1.2 to 4.4
  TH2D *h2_Ppip_thetapip_cut =  new TH2D("h2_Ppip_thetapip_cut","Data cut: Pi+ Momentum vs. Pi+ Angle;#theta_pip;Momentum (GeV/c)",400,0.,Wthetahistomax,400,0.,Phistomax);
  TH2D *h2_Ppim_thetapim_cut =  new TH2D("h2_Ppim_thetapim_cut","Data cut: Pi- Momentum vs. Pi- Angle;#theta_pim;Momentum (GeV/c)",400,0.,Wthetahistomax,400,0.,Phistomax);
  TH2D *h2_Ppip_thetapipq_cut = new TH2D("h2_Ppip_thetapipq_cut","Data cut: Pi+ Momentum vs. Pipq Angle;#theta_pipq;Momentum (GeV/c)",400,0.,Wthetahistomax,400,0.,Phistomax);
  TH2D *h2_Ppim_thetapimq_cut = new TH2D("h2_Ppim_thetapimq_cut","Data cut: Pi- Momentum vs. Pimq Angle;#theta_pimq;Momentum (GeV/c)",400,0.,Wthetahistomax,400,0.,Phistomax);
  //FTOF panel 2
  TH2D *h2_Ppip_thetapip_scin1 = new TH2D("h2_Ppip_thetapip_scin1","FTOF1: Pi+ Momentum vs. Pi+ Angle;#theta_pip;Momentum (GeV/c)",400,0.,Wthetahistomax,400,0.,Phistomax);
  TH2D *h2_Ppip_thetapip_scin2 = new TH2D("h2_Ppip_thetapip_scin2","FTOF2: Pi+ Momentum vs. Pi+ Angle;#theta_pip;Momentum (GeV/c)",400,0.,Wthetahistomax,400,0.,Phistomax);
  //QE proton analysis
  TH2D *h2_phi_e_phi_pip = new TH2D("h2_phi_e_phi_pip",";phi_pi (deg);phi_e",360,-180.,180.,360,-180.,180.);
  TH2D *h2_E_FTOF_P_pip = new TH2D("h2_E_FTOF_P_pip",";Momentum (GeV/c);E_{FTOF}",400,0.,Phistomax,400,0.,100.);
  TH2D *h2_E_FTOF_P_pipz = new TH2D("h2_E_FTOF_P_pipz",";Momentum (GeV/c);E_{FTOF}",400,0.,Phistomax,800,0.,20.);
  TH2D *h2_E_FTOF_P_pipQE = new TH2D("h2_E_FTOF_P_pipQE",";Momentum (GeV/c);E_{FTOF}",400,0.,Phistomax,400,0.,100.);
  TH2D *h2_phi_e_phi_pim = new TH2D("h2_phi_e_phi_pim",";phi_pi (deg);phi_e",360,-180.,180.,360,-180.,180.);
  TH2D *h2_E_FTOF_P_pim = new TH2D("h2_E_FTOF_P_pim",";Momentum (GeV/c);E_{FTOF} (GeV)",400,0.,Phistomax,400,0.,100.);
  TH2D *h2_beta_t_FTOF = new TH2D("h2_beta_t_FTOF",";t (?);",400,0.,200.,440,-10.,10.1);
  TH2D *h2_beta_t_FTOFa = new TH2D("h2_beta_t_FTOFa",";t (?);",400,0.,200.,440,-10.,10.1);
  TH2D *h2_E_FTOFa_P_pip = new TH2D("h2_E_FTOFa_P_pip",";Momentum (GeV/c);E_{FTOF1A}",400,0.,Phistomax,400,0.,100.);
  TH1D *h_W_data_pipFTOF1B = new TH1D("h_W_data_pipFTOF1B",";W (GeV/c^{2});",400,0.5,Whistomax);
  TH1D *h_W_data_pipFTOF1A = new TH1D("h_W_data_pipFTOF1A",";W (GeV/c^{2});",400,0.5,Whistomax);
  TH1D *h_W_data_pipFTOFnoE = new TH1D("h_W_data_pipFTOFnoE",";W (GeV/c^{2});",400,0.5,Whistomax);
  TH1D *h_Wpip_data_hiTQ = new TH1D("h_Wpip_data_hiTQ",";W (GeV/c^{2});",400,0.5,Whistomax);
  TH1D *h_ThetaQpip_dataQE = new TH1D("h_ThetaQpip_dataQE",";theta_q (deg);",400,0.,80.);
  TH1D *h_ThetaQpip_dataNoQE = new TH1D("h_ThetaQpip_dataNoQE",";theta_q (deg);",400,0.,80.);
  TH1D *h_ThetaQpip_dataFTOF1B = new TH1D("h_ThetaQpip_dataFTOF1B",";theta_q (deg);",400,0.,80.);
  TH2D *h2_q2_Wpip_QE = new TH2D("h2_q2_Wpip_QE",";W (GeV/c^{2});Q^{2} (GeV^{2})",320,0.8,1.1,400,0.,4.4);//1.2 to 4.4
  TH2D *h2_beta_P_pip        = new TH2D("h2_beta_P_pip",";P (GeV/c);",410,0.,4.1,200,0.9,1.1);//440
  TH2D *h2_beta_P_pipzoom    = new TH2D("h2_beta_P_pipzoom",";P (GeV/c);",410,0.,4.1,440,0.8,1.5);//0.8,1.5
  TH2D *h2_beta_P_pim        = new TH2D("h2_beta_P_pim",";P (GeV/c);c",410,0.,4.1,200,0.9,1.1);
  TH2D *h2_beta_P_p          = new TH2D("h2_beta_P_p",";P (GeV/c);c",410,0.,4.1,440,0.,1.1);
  TH2D *h2_beta_P_pipHiTQ    = new TH2D("h2_beta_P_pipHiTQ",";P (GeV/c);",410,0.,4.1,440,0.95,1.02);
  TH2D *h2_beta_P_pipLowW    = new TH2D("h2_beta_P_pipLowW",";P (GeV/c);",410,0.,4.1,440,0.95,1.02);
  TH2D *h2_beta_P_pip_ray    = new TH2D("h2_beta_P_pip_ray",";P (GeV/c);",410,0.,4.1,440,0.95,1.02);
  TH2D *h2_beta_P_pipchi2    = new TH2D("h2_beta_P_pipchi2",";P (GeV/c);",410,0.,4.1,440,0.8,1.5);
  TH2D *h2_beta_P_pipchi3    = new TH2D("h2_beta_P_pipchi3",";P (GeV/c);",410,0.,4.1,440,0.8,1.5);
  TH2D *h2_DC_y_vs_x_data    = new TH2D("h2_DC_y_vs_x_data",";x;y",400,-200.,200.,400,-200.,200.);

  TH1D *h_beta_P8_pip       = new TH1D("h_beta_P8_pip",";;",760,0.,1.9);
  TH1D *h_beta_P8_pim       = new TH1D("h_beta_P8_pim",";;",760,0.,1.9);
  TH1D *h_beta_P8_pipcut      = new TH1D("h_beta_P8_pipcut",";;",200,0.9,1.1);
  TH1D *h_beta_P8_pimcut      = new TH1D("h_beta_P8_pimcut",";;",200,0.9,1.1);
  TH1D *h_beta_P10_pip       = new TH1D("h_beta_P10_pip",";;",760,0.,1.9);
  TH1D *h_beta_P10_pim       = new TH1D("h_beta_P10_pim",";;",760,0.,1.9);
  TH1D *h_beta_P10_pipcut      = new TH1D("h_beta_P10_pipcut",";;",760,0.,1.9);
  TH1D *h_beta_P10_pimcut      = new TH1D("h_beta_P10_pimcut",";;",760,0.,1.9);

  //MC Simulation
  TH2D *h2_theta_q2_MC = new TH2D("h2_theta_q2_MC","GENIE: Electron Angle vs. Q^{2};Q^{2};#theta",400,0.,4.4,400,0.,Wthetahistomax);//1.2 to 4.4
  TH2D *h2_ecal_MC = new TH2D("h2_ecal_MC","GENIE: Electron Momentum Sampling Fraction;Electron Momentum (GeV/c);Sampling Fraction",400,0.,4.4,1000,0.,0.5);//1.2 to 4.4
  TH2D *h2_theta_mom_MC = new TH2D("h2_theta_mom_MC","GENIE: Electron Angle vs. Momentum;Momentum (GeV/c);#theta",400,0.,5.,400,0.,Wthetahistomax);
  TH2D *h2_q2_W_MC = new TH2D("h2_q2_W_MC","GENIE: Q^{2} vs. W;W (GeV/c^{2});Q^{2}",400,0.,Whistomax,400,0.,4.4);//1.2 to 4.4
  TH2D *h2_theta_W_MC = new TH2D("h2_theta_W_MC","GENIE: W vs #theta;W (GeV/c^{2});#theta",400,0.,Whistomax,400,0.,Wthetahistomax);
  TH2D *h2_theta_omega_MC = new TH2D("h2_theta_omega_MC","GENIE: Energy Transfer vs. #theta;#omega (GeV);#theta",1000,0.,omegahistomax,400,0.,Wthetahistomax);
  TH2D *h2_mult_p_pis_MC = new TH2D("h2_mult_p_pis_MC","GENIE: Multiplicity of Protons vs. Pions #pi^{#pm};Proton Multiplicity;#pi^{#pm} Multiplicity",5,-0.5,4.5,5,-0.5,4.5);
  TH2D *h2_Ppip_theta_MC = new TH2D("h2_Ppip_theta_MC","GENIE: Pi+ Momentum vs. Electron Angle;#theta_e;Momentum (GeV/c)",400,0.,Wthetahistomax,400,0.,Phistomax);
  TH2D *h2_Ppim_theta_MC = new TH2D("h2_Ppim_theta_MC","GENIE: Pi- Momentum vs. Electron Angle;#theta_e;Momentum (GeV/c)",400,0.,Wthetahistomax,400,0.,Phistomax);
  //
  TH2D *h2_Ppip_thetapip_MC     =  new TH2D("h2_Ppip_thetapip_MC","GENIE: Pi+ Momentum vs. Pi+ Angle;#theta_pip;Momentum (GeV/c)",400,0.,Wthetahistomax,400,0.,Phistomax);
  TH2D *h2_Ppim_thetapim_MC     =  new TH2D("h2_Ppim_thetapim_MC","GENIE: Pi- Momentum vs. Pi- Angle;#theta_pim;Momentum (GeV/c)",400,0.,Wthetahistomax,400,0.,Phistomax);
  TH2D *h2_q2_W_MCrad = new TH2D("h2_q2_W_MCrad","GENIE rad: Q^{2} vs. W;W (GeV/c^{2});Q^{2}",400,0.,Whistomax,400,0.,4.4);//1.2 to 4.4
  TH2D *h2_theta_W_MCrad = new TH2D("h2_theta_W_MCrad","GENIE rad: W vs #theta;W (GeV/c^{2});#theta",400,0.,Whistomax,400,0.,Wthetahistomax);

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //2D Histograms of particular topologies without cuts on various kinematic variables (what should the cuts be?)
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  //TH2D *h2_theta_mom_opg   = new TH2D("h2_theta_mom_opg","onepigen: Electron Angle vs. Momentum;Momentum (GeV/c);#theta",400,0.,5.,400,0.,Wthetahistomax);
  TH1D *h_q2_opg = new TH1D("h_q2_opg","Onepigen;Q^{2};Counts",400,0.,q2histomax);
  TH1D *h_Ppi_opg = new TH1D("h_Ppi_opg","Onepigen;;",400,0.,4.);
  TH1D *h_theta_piq_opg = new TH1D("h_theta_piq_opg","Onepigen;;",600,0.,60.);
  TH1D *h_theta_pi_opg = new TH1D("h_theta_pi_opg","Onepigen;;",600,0.,60.);
  THStack *hs_twf = new THStack("hs_twf","onepigen test: W < 1.7 GeV and Fermi motion;W (GeV/c^{2});Counts");
  TH1D *h_tw =      new TH1D("h_tw","W < 1.7 GeV;W (GeV/c^{2});Counts",400,1.,Whistomax);
  TH1D *h_tf =      new TH1D("h_tf","Fermi motion;W (GeV/c^{2});Counts",400,1.,Whistomax);
  TH1D *h_twf =     new TH1D("h_twf","W < 1.7 GeV and Fermi motion;W (GeV/c^{2});Counts",400,1.,Whistomax);
  //
  TH2D *h2_Ppip_thetapip_opg     =  new TH2D("h2_Ppip_thetapip_opg","ONEPIGEN: Pi+ Momentum vs. Pi+ Angle;#theta_pip;Momentum (GeV/c)",400,0.,Wthetahistomax,400,0.,Phistomax);
  TH2D *h2_Ppim_thetapim_opg     =  new TH2D("h2_Ppim_thetapim_opg","ONEPIGEN: Pi- Momentum vs. Pi- Angle;#theta_pim;Momentum (GeV/c)",400,0.,Wthetahistomax,400,0.,Phistomax);
  TH2D *h2_Ppip_thetapip_opg2Axx = new TH2D("h2_Ppip_thetapip_opg2Axx","OPG 2AXX: Pi+ Momentum vs. Pi+ Angle;#theta_pip;Momentum (GeV/c)",400,0.,Wthetahistomax,400,0.,Phistomax);
  TH2D *h2_Ppip_thetapip_opg2xAx = new TH2D("h2_Ppip_thetapip_opg2xAx","OPG 2XAX: Pi+ Momentum vs. Pi+ Angle;#theta_pip;Momentum (GeV/c)",400,0.,Wthetahistomax,400,0.,Phistomax);
  TH2D *h2_Ppim_thetapim_opg2Axx = new TH2D("h2_Ppim_thetapim_opg2Axx","OPG 2AXX: Pi- Momentum vs. Pi- Angle;#theta_pim;Momentum (GeV/c)",400,0.,Wthetahistomax,400,0.,Phistomax);
  TH2D *h2_Ppim_thetapim_opg2xAx = new TH2D("h2_Ppim_thetapim_opg2xAx","OPG 2XAX: Pi- Momentum vs. Pi- Angle;#theta_pim;Momentum (GeV/c)",400,0.,Wthetahistomax,400,0.,Phistomax);
  TH2D *h2_thetaQ_Wpip_opg       = new TH2D("h2_thetaQ_Wpip_opg","ONEPIGEN: #theta_q vs. W pip;W (GeV/c^{2});#theta_q",400,0.,Whistomax,400,0.,ThetaQhistomax);
  TH2D *h2_thetaQ_ThetaPip_opg   = new TH2D("h2_thetaQ_ThetaPip_opg","ONEPIGEN: #theta_q vs. #theta_pip;#theta_pip;#theta_q",400,0.,Wthetahistomax,400,0.,ThetaQhistomax);
  TH2D *h2_thetaQ_omegapip_opg   = new TH2D("h2_thetaQ_omegapip_opg","ONEPIGEN: #theta_q vs. #omega_pip;#omega (GeV/c^{2});#theta_q",350,0.,3.5,400,0.,ThetaQhistomax);
  TH2D *h2_thetaQ_vzpip_opg      = new TH2D("h2_thetaQ_vzpip_opg","ONEPIGEN: #theta_q vs. vz_e pip; vz_e (cm);#theta_q",400,-8.,2.,400,0.,ThetaQhistomax);
  TH2D *h2_thetaQ_Wpim_opg       = new TH2D("h2_thetaQ_Wpim_opg","ONEPIGEN: #theta_q vs. W pim;W (GeV/c^{2});#theta_q",400,0.,Whistomax,400,0.,ThetaQhistomax);
  TH2D *h2_thetaQ_ThetaPim_opg   = new TH2D("h2_thetaQ_ThetaPim_opg","ONEPIGEN: #theta_q vs. #theta_pim;#theta_pim;#theta_q",400,0.,Wthetahistomax,400,0.,ThetaQhistomax);
  TH2D *h2_thetaQ_omegapim_opg   = new TH2D("h2_thetaQ_omegapim_opg","ONEPIGEN: #theta_q vs. #omega_pim;#omega (GeV/c^{2});#theta_q",350,0.,3.5,400,0.,ThetaQhistomax);
  TH2D *h2_thetaQ_vzpim_opg      = new TH2D("h2_thetaQ_vzpim_opg","ONEPIGEN: #theta_q vs. vz_e pim; vz_e (cm);#theta_q",400,-8.,2.,400,0.,ThetaQhistomax);
  //Tail analysis
  TH1D *h_opg_tail_vzpip =      new TH1D("h_opg_tail_vzpip","Tail;vzpip (cm);Counts",480,-9.,3.);
  TH1D *h_opg_nota_vzpip =      new TH1D("h_opg_nota_vzpip","No Tail;vzpip (cm);Counts",480,-9.,3.);
  TH1D *h_opg_tail_omega =      new TH1D("h_opg_tail_omega","Tail;omega (GeV);Counts",350,0.,3.5);
  TH1D *h_opg_nota_omega =      new TH1D("h_opg_nota_omega","No Tail;omega (GeV);Counts",350,0.,3.5);
  //
  TH1D *h_phi_e_opg  = new TH1D("h_phi_e_opg","Onepigen;;",360,-180.,180.);
  TH1D *h_phi_pi_opg = new TH1D("h_phi_pi_opg","Onepigen;;",360,-180.,180.);
  TH1D *h_pip_chi2pid= new TH1D("h_pip_chi2pid",";;",200,-5.,5.);
  TH1D *h_pim_chi2pid= new TH1D("h_pim_chi2pid",";;",200,-5.,5.);

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Counters
  Int_t counter_data = 0, counter_MC = 0, counter_gen=0, counter_opg = 0, counter_opgm = 0, counter_opgn = 0, counter_opgmn = 0, counter_MC1 = 0, counter_MCrad = 0;//Good electrons

  ///Cross sections
    Int_t Wbini = 48, Q2bini = 4, models = 6, sectors = 5, Nrec = 2, TPIQi = 4, Pbini = 4;
    Int_t iDATA = 5, iGENIE = 4, iGEN = 3, iGENIEr = 2, iOPGr = 1, iOPGn = 0;
    Int_t i_pi = -1, i_pip = 0, i_pim = 0, canvas = 1, mark = 8, cs_data = 0, cs_MC = 0, cs_opg = 0;
    Int_t cs_opgm = 0, cs_opgn = 0, cs_opgmn = 0; ///files = 7
    Double_t theta_piq = 0.0, theta_piqmax =0.0;
    Double_t P_pi = 0.0, Theta_pi = 0.0, Phi_pi = 0.0, Chi2pid_pi = 0., chi2pid_pip = -0.061, chi2pid_pim = -0.019;
    Double_t Vperp_e = 0.0, Vperp_pi = 0.0, Vdiff_pi = 0.0;
    ///cut values
    Double_t VWcut = 14.0, vz_e_min = -8.0, vz_e_max = 2.0, Vperp_ecut = 3.2, Wcut = 1.08, euler = 2.71828;
    Double_t Vperp_pipcut[3] = {10.5,8.3,5.5}, Vperp_pimcut[2] = {7.4,5.8};
    Double_t Vperp_pip_ms[3][2] = {{0.54,1.8},{0.79,1.6},{0.52,1.3}}, Vperp_pim_ms[2][2] = {{0.040,1.6},{0.073,1.3}}, PiCutSig = 3.0;
    Double_t u_eID = 0.015, u_hID = 0.014, u_Vze = 0.036, u_Vre = 0.012, u_Vrpi = 0.011, u_Vdiff = 0.02, u_TW = 0.01;
    ///CS values
    Double_t CS_Wbin = (Double_t)Wbini, CS_Wmin = 1.075, CS_Wmax = 2.275;// each bin has 0.025 GeV
	Double_t CS_Q2bin = (Double_t)Q2bini, CS_Q2min = 0.6, CS_Q2max = 1.5, CS_Tpiq = (Double_t)TPIQi, CS_Sec = (Double_t)sectors;
	Double_t CS_Pbin = (Double_t)Pbini, CS_Pmin = 0.3, CS_Pmax = 3.3;//990000
	Double_t cm2ub = 1E+30, ub2nb = (6./5.)*1E+3, bin_vol = 1.0; Double_t totalCS = 0.052, gen_evts = 11110000., opg_gen_events[4] = {14.E+6,10.E+6,17.E+6,14.E+6};
	Double_t totalCS_Dipole = 0.4321, totalCS_Rarita = 0.4621, dipole_evts = 1.812E+8, rarita_evts = 1.1111E+8; //ub 1.E+7,1.E+7,5.E+6,5.E+6 Rar=2.616
	Double_t totalCS_DipRad = 0.4321, totalCS_RarRad = 0.4621, diprad_evts = 0.972E+8, rarrad_evts = 1.1111E+8; diprad_evts = 0.576E+8;
	if(bAllEvents == false){dipole_evts = 1.812E+5; diprad_evts = 0.972E+5;}
	Double_t totalCS_MC[2] = {totalCS_Dipole, totalCS_Rarita}, MC_evts[2] = {dipole_evts,rarita_evts}, MC_rad_evts[2] = {diprad_evts,rarrad_evts};
	Double_t tgt_den = 0.169, tgt_len = 5.00, N_A = 6.02214076E+23, tgt_mol = 2.014, q_e = 1.602176634E-19, Q_tot = 6078.43*1E-9;//238300
	Double_t L = 791960*1E-9 * tgt_den * tgt_len * N_A / (q_e * tgt_mol); Double_t opgN_gen_events = 600000;
	/// W and Q2 bin values
	Double_t dW = (CS_Wmax - CS_Wmin)/CS_Wbin, dT = 10.0, dQ2[4] = {0.3, 0.4, 0.5, 0.6}, dP[4] = {0.3, 0.5, 0.7, 0.5};
	Double_t Wi = CS_Wmin + dW/2, w = CS_Wmin + dW, t = 0.0, dt = dT, jpiq = dT;
	Double_t WbinEdge[4][2] = {{1.1,2.275},{1.1,2.175},{1.1,2.075},{1.1,1.9}};
	/// Theta_PiQ bin values
	Double_t ThetaPiQ[5] = {0.0, 12.0, 22.0, 32.0, 46.0};//{0.0, 20., 30.0, 40.0, 60.0};
	Double_t ThetaPi[2][6] = {{7., 13., 19., 25., 31.,32.},{15.,20.,26.,32.,40.,45.}};
	Double_t dPt[2][5] = {{0.75,1.1,1.5,2.2,2.7},{0.6,0.9,1.5,1.9,2.4}};
	/// Sector stuff
  	Int_t sec_i = 0;
  	Double_t y_i[Pbini+1][TPIQi+1][Nrec][sectors][models][Q2bini][Wbini];
  	Double_t tCS[Pbini+1][TPIQi+1][Nrec][sectors][3][Q2bini][Wbini];
  	Double_t y_sum[Pbini+1][TPIQi+1][Nrec][models][Q2bini][Wbini];
  	Double_t y_ave[Pbini+1][TPIQi+1][Nrec][models][Q2bini][Wbini];
	Double_t var_stat[Pbini+1][TPIQi+1][Nrec][sectors][models][Q2bini][Wbini];
  	Double_t var[Pbini+1][TPIQi+1][Nrec][models][Q2bini][Wbini];
    /// Graph Stuff
    Double_t CS_log = 0;
	TVector3 vec_e, vec_q, vec_pi, vec_b;
  	char c_NAME[30];
  	TH1D *hr_CS[Pbini][TPIQi][Nrec][Q2bini]; TH1D *hr_AC[Pbini][TPIQi][Nrec][Q2bini]; TH1D *hr_CSmc[Pbini][TPIQi][Nrec][Q2bini];
  	TH1D *hr_ACq[Nrec][Q2bini]; TH1D *hr_ACp[Pbini][Nrec][Q2bini]; TH1D *hr_ACt[TPIQi][Nrec][Q2bini];
  	TH1D *hr_RCq[Nrec][Q2bini]; TH1D *hr_RCp[Pbini][Nrec][Q2bini]; TH1D *hr_RCt[TPIQi][Nrec][Q2bini];
  	TH1D *hr_ARC[Pbini][TPIQi][Nrec][Q2bini]; TH1D *hr_ARCq[Nrec][Q2bini]; TH1D *hr_ARCp[Pbini][Nrec][Q2bini]; TH1D *hr_ARCt[TPIQi][Nrec][Q2bini];
  	Double_t RadEffCor[Pbini][TPIQi][Nrec][Q2bini];
    TH2D *h2_RC_chi2pip = new TH2D("h2_RC_chi2pip",";chi2/NDF;RC",14,0.,3.5,44,0.,1.1);
    TH2D *h2_RC_chi2pim = new TH2D("h2_RC_chi2pim",";chi2/NDF;RC",14,0.,3.5,44,0.,1.1);
    Double_t RadFitAve[2] = {0.,0.}, RadFitBins[2] = {40.,40.}, RadFitUnc[2] = {0.,0.}, RadFitErr[2] = {0.,0.};
    if(bQorNot == false){RadFitBins[0] = 44.; RadFitBins[1] = 52.;}//number of bins that will have rad corr fits for theta_pi plots
	for(Int_t iPmom = 0; iPmom < Pbini+1; iPmom++)
    {
        for(Int_t iThPiQ = 0; iThPiQ < TPIQi+1; iThPiQ++)
        {
            for(Int_t pm = 0; pm < Nrec; pm++)
            {
              for(Int_t f = 0; f < models; f++)
              {
                for(Int_t iQ2 = 0; iQ2 < Q2bini; iQ2++)
                {
                    if(f==0)
                    {
                        RadEffCor[iPmom][iThPiQ][pm][iQ2] = 0.;
                        sprintf(c_NAME,"hr_CS[%d][%d][%d][%d]",iPmom,iThPiQ,pm,iQ2);
                        hr_CS[iPmom][iThPiQ][pm][iQ2] = new TH1D(c_NAME,"",Wbini,CS_Wmin,CS_Wmax);
                        sprintf(c_NAME,"hr_AC[%d][%d][%d][%d]",iPmom,iThPiQ,pm,iQ2);
                        hr_AC[iPmom][iThPiQ][pm][iQ2] = new TH1D(c_NAME,"",Wbini,CS_Wmin,CS_Wmax);
                        sprintf(c_NAME,"hr_CSmc[%d][%d][%d][%d]",iPmom,iThPiQ,pm,iQ2);
                        hr_CSmc[iPmom][iThPiQ][pm][iQ2] = new TH1D(c_NAME,"",Wbini,CS_Wmin,CS_Wmax);
                        sprintf(c_NAME,"hr_ARC[%d][%d][%d][%d]",iPmom,iThPiQ,pm,iQ2);
                        hr_ARC[iPmom][iThPiQ][pm][iQ2] = new TH1D(c_NAME,"",Wbini,CS_Wmin,CS_Wmax);
                        if(iPmom==0 && iThPiQ<TPIQi)
                        {
                            sprintf(c_NAME,"hr_ACt[%d][%d][%d]",iThPiQ,pm,iQ2);
                            hr_ACt[iThPiQ][pm][iQ2] = new TH1D(c_NAME,"",Wbini,CS_Wmin,CS_Wmax);
                            sprintf(c_NAME,"hr_RCt[%d][%d][%d]",iThPiQ,pm,iQ2);
                            hr_RCt[iThPiQ][pm][iQ2] = new TH1D(c_NAME,"",Wbini,CS_Wmin,CS_Wmax);
                            sprintf(c_NAME,"hr_ARCt[%d][%d][%d]",iThPiQ,pm,iQ2);
                            hr_ARCt[iThPiQ][pm][iQ2] = new TH1D(c_NAME,"",Wbini,CS_Wmin,CS_Wmax);
                            if(iThPiQ==0)
                            {
                                sprintf(c_NAME,"hr_ACq[%d][%d]",pm,iQ2);
                                hr_ACq[pm][iQ2] = new TH1D(c_NAME,"",Wbini,CS_Wmin,CS_Wmax);
                                sprintf(c_NAME,"hr_RCq[%d][%d]",pm,iQ2);
                                hr_RCq[pm][iQ2] = new TH1D(c_NAME,"",Wbini,CS_Wmin,CS_Wmax);
                                sprintf(c_NAME,"hr_ARCq[%d][%d]",pm,iQ2);
                                hr_ARCq[pm][iQ2] = new TH1D(c_NAME,"",Wbini,CS_Wmin,CS_Wmax);
                            }
                        }
                        if(iThPiQ==0 && iPmom<Pbini)
                        {
                            sprintf(c_NAME,"hr_ACp[%d][%d][%d]",iPmom,pm,iQ2);
                            hr_ACp[iPmom][pm][iQ2] = new TH1D(c_NAME,"",Wbini,CS_Wmin,CS_Wmax);
                            sprintf(c_NAME,"hr_RCp[%d][%d][%d]",iPmom,pm,iQ2);
                            hr_RCp[iPmom][pm][iQ2] = new TH1D(c_NAME,"",Wbini,CS_Wmin,CS_Wmax);
                            sprintf(c_NAME,"hr_ARCp[%d][%d][%d]",iPmom,pm,iQ2);
                            hr_ARCp[iPmom][pm][iQ2] = new TH1D(c_NAME,"",Wbini,CS_Wmin,CS_Wmax);
                        }
                    }
                    for(Int_t iW = 0; iW < Wbini; iW++)
                    {
                        for(Int_t s = 0; s < sectors; s++)
                        {
                            y_i[iPmom][iThPiQ][pm][s][f][iQ2][iW] = 0.;
                            tCS[iPmom][iThPiQ][pm][s][f][iQ2][iW] = 0.;
                            var_stat[iPmom][iThPiQ][pm][s][f][iQ2][iW] = 0.;
                        }
                        y_sum[iPmom][iThPiQ][pm][f][iQ2][iW] = 0.;
                        y_ave[iPmom][iThPiQ][pm][f][iQ2][iW] = 0.;
                        var[iPmom][iThPiQ][pm][f][iQ2][iW] = 0.;
                    }
                }
              }
            }
        }
    }
  for(Int_t id=0;id<resnum;id++)
  {
    sprintf(c_NAME,"h_W_resid_MC[%d]",id);
    h_W_resid_MC[id] = new TH1D(c_NAME,"",400,0.5,Whistomax);
  }
  Int_t i_pipr = 0, i_piprF = 0, i_pimr = 0, i_pimrF = 0, i_pipn = 0, i_pipnF = 0, i_pimn = 0, i_pimnF = 0;
  Int_t i_pitot = 0, i_pisum = 0, i_pitom = 0, i_pisut = 0, i_pitop = 0, i_pisup = 0, i_piton = 0, i_pisun = 0;
  Double_t pr01 = 3.72620076E-2, pr02 = 3.47743630E-2, pr03 = 4.93952148E-2, pr04 = 3.87403741E-2, pr05 = 3.55959348E-2, pr06 = 3.82943228E-2, pr07 = 3.39816809E-2;
  Double_t pr08 = 3.57538015E-2, pr09 = 3.72311212E-2, pr10 = 3.81282791E-2, pr11 = 3.89013924E-2, pr12 = 3.67460363E-2, pr13 = 3.94949019E-2, pr14 = 3.81195620E-2;
  Double_t mr01 = 2.56068986E-2, mr02 = 2.96644252E-2, mr03 = 2.59543750E-2, mr04 = 2.99047828E-2, mr05 = 2.88771614E-2, mr06 = 2.83538569E-2;
  Double_t mr07 = 3.01340446E-2, mr08 = 2.81162374E-2, mr09 = 2.90096998E-2, mr10 = 2.96688192E-2, mr11 = 2.79425681E-2, mr12 = 2.92918216E-2;
  Double_t pn01=4.12995666E-2,pn02=4.04881611E-2,pn03=4.56248149E-2,pn04=4.29100879E-2,pn05=4.27147150E-2,pn06=4.06351313E-2,pn07=4.31839041E-2,pn08=4.15647402E-2,pn09=4.25501764E-2;
  Double_t pn10=3.89829502E-2,pn11=4.12300527E-2,pn12=4.53321710E-2,pn13=4.34181802E-2,pn14=4.62041572E-2,pn15=4.42253500E-2,pn16=4.18082289E-2,pn17=3.80347520E-2,pn18=4.04049121E-2;
  Double_t pntw=4.20624651E-2,pntF=4.1958414E-2,pntwF=4.19772863E-2;
  Double_t mn01 = 3.38587388E-2, mn02 = 3.22229117E-2, mn03 = 3.26624736E-2, mn04 = 3.35078128E-2, mn05 = 3.34689766E-2, mn06 = 3.35908160E-2, mn07 = 3.46761085E-2;
  Double_t mn08 = 2.77200956E-2, mn09 = 3.15224752E-2, mn10 = 3.15674283E-2, mn11 = 3.25770155E-2, mn12 = 3.96888182E-2, mn13 = 4.19930257E-2, mn14 = 3.48580852E-2;
  Double_t intCS_pipr[14] = {pr01,pr02,pr03,pr04,pr05,pr06,pr07,pr08,pr09,pr10,pr11,pr12,pr13,pr14};
  Double_t intCS_pimr[10] = {mr01,mr02,mr03,mr04,mr05,mr06,mr07,mr08,mr09,mr12};///12
  Double_t intCS_pipn[17] = {pn01,pn02,pn03,pn05,pn06,pn07,pn08,pn09,pn10,pn11,pn12,pn13,pn14,pn15,pn16,pn17,pn18};
  Double_t intCS_pimn[14] = {mn01,mn02,mn03,mn04,mn05,mn06,mn07,mn08,mn09,mn10,mn11,mn12,mn13,mn14};
  Int_t piprN[14] = {60903,61288,59497,60353,60509,60703,61092,60480,60689,60961,60430,60820,61050,61360};
  Int_t pimrN[10] = {33245,33045,32928,32874,32993,33288,32949,33112,32918,32991};
  Int_t pipnN[17] = {64757,65285,64913,65101,64763,64696,65016,64970,64836,64780,64663,64890,64819,64678,65024,65152,64678};
  Int_t pimnN[14] = {36413,36040,35944,36036,36119,35858,36537,36177,36110,36293,36079,36116,35944,36150};
  i_pitot = 60903+61288+59497+60353+60509+60703+61092+60480+60689+60961+60430+60820+61050+61360;
  i_pitom = 33245+33045+32928+32874+32993+33288+32949+33112+32918+32991;
  i_pitop = 64757+65285+64913+65101+64763+64696+65016+64970+64836+64780+64663+64890+64819+64678+65024+65152+64678;
  i_piton = 36413+36040+35944+36036+36119+35858+36537+36177+36110+36293+36079+36116+35944+36150;
  std::cout << "I'm above the data" << std::endl;
 if(bData==true)
 {
  //////////////////////////////////////////////////////////////////////////////
  //         DATA
  //////////////////////////////////////////////////////////////////////////////
  while(chain_data.Next())
    {
      v4_beam_data.SetE(beamE);
      v4_beam_data.SetPz(beamE);
      if( beamE > 1e-4)
	      {
          //if( counter_data == 0)
          //  cout<<"Beam energy set manually: "<<beamE<<endl;
          v4_beam_data.SetE(beamE);
          v4_beam_data.SetPz(beamE);
	      }

      //TLorentzVector el_mc;
      //c12_data->mcparts()->setEntry(0);

      //el_mc.SetXYZM(c12_data->mcparts()->getPx(), c12_data->mcparts()->getPy(),c12_data->mcparts()->getPz(),db->GetParticle(11)->Mass());

      // get particles by type
      auto electrons=c12_data->getByID(11);
      auto protons=c12_data->getByID(2212);
      auto neutrons=c12_data->getByID(2112);
      auto pips=c12_data->getByID(211);
      auto pims=c12_data->getByID(-211);

      //std::cout << electrons.size() << std::endl;

      if(electrons.size()==1 && electrons[0]->trk(DC)->getSector() != 4 && (pips.size()>=1 || pims.size()>=1))///changed to single pi+ and/or pi-
	      {
          Double_t energy =  (electrons[0]->cal(PCAL)->getEnergy() +  electrons[0]->cal(ECIN)->getEnergy() +  electrons[0]->cal(ECOUT)->getEnergy())/electrons[0]->getP();

          bool ecal = ( (electrons[0]->cal(ECIN)->getLv() >= 14 && electrons[0]->cal(ECIN)->getLw() >= 14) && (energy > 0.18 && energy < 0.28) && electrons[0]->getP() > 0.5 && electrons[0]->getP() < 10 );

          SetLorentzVector(v4_el,electrons[0]);

          TLorentzVector v4_q = v4_beam_data - v4_el; //photon  4-vector

          Double_t q2        = -v4_q.M2();
          Double_t x_b       = q2/(2 * mass_p * v4_q.E() ); //Bjorken-x
          Double_t y         = -v4_q.P() + sqrt( v4_q.E()*v4_q.E() + 2*v4_q.E()*mass_p); //y (Bjorken-y), the scaling variable; from Jourdan, eq. 20
          Double_t vz_e      = electrons[0]->par()->getVz();
          Double_t W         = sqrt(mass_p*mass_p - q2 + 2*v4_q.E()*mass_p); //Invariant mass, assuming electron off of a standing proton
          //std::cout << W << std::endl;
          //std::cout << vz_e << std::endl;
          h_el_vz_h_data->Fill(vz_e);

          if(W < Wcut) continue; /// W min cut (removes QE "pion" [perhaps coincident proton] events with high ThetaQ)
          else if(W > CS_Wmax) continue; /// Save time

          //Electron quality cuts
          bool targ_1 = abs( vz_e - (-1.48)) < 1 * 1.25;
          bool targ_2 = abs( vz_e - (-6.3))  < 1 * 1.25;
          h2_ecal_data->Fill(electrons[0]->getP(),energy);

          ///***********************************************************************************************************************************************************************************************
          //bin: pion momentum, theta_piq, charge (pi+=0,pi-=1), file, Q2, invariant mass
            Int_t iPmom = -1, iThPiQ = -1, iCharge = -1, iMod = -1, iQ2 = -1, iW = -1;
            ///Electron vertex calculation
            vec_q.SetMagThetaPhi(v4_q.P(),v4_q.Theta(),v4_q.Phi());
            Vperp_e = sqrt(pow(electrons[0]->par()->getVx(),2) + pow(electrons[0]->par()->getVy(),2));
            P_pi = 0.0; i_pi = 0; Double_t P_can = 0.0;// pion momentum, pion # for the event, momentum of pion candidate
            Double_t Lv = electrons[0]->cal(PCAL)->getLv(), Lw = electrons[0]->cal(PCAL)->getLw();
            Int_t i_FDstat = 0;
            Double_t beta = 0.;

            ///Electron EC fiducial and vertex cuts
            if(eCutter(VWcut, Lv, Lw, vz_e, vz_e_min, vz_e_max, Vperp_e, Vperp_ecut)!=0) continue;
            //DC fiducial cuts for electron
            bool el_fid_lay_1 = dcFid.DC_pi_fid(electrons[0]->traj(DC,6)->getX(), electrons[0]->traj(DC,6)->getY(),  1, electrons[0]->trk(DC)->getSector(), 1,bending);
            bool el_fid_lay_2 = dcFid.DC_pi_fid(electrons[0]->traj(DC,18)->getX(),electrons[0]->traj(DC,18)->getY(), 1, electrons[0]->trk(DC)->getSector(), 2,bending);
            bool el_fid_lay_3 = dcFid.DC_pi_fid(electrons[0]->traj(DC,36)->getX(),electrons[0]->traj(DC,36)->getY(), 1, electrons[0]->trk(DC)->getSector(), 3,bending);
            if( !(el_fid_lay_1 && el_fid_lay_2 && el_fid_lay_3) ) continue;

            ///Highest momentum Pion Calculator
            if(pips.size()>=1)
            {
                for(Int_t i=0;i<pips.size();i++)
                {
                  if(abs(pips[i]->getChi2Pid()-chi2pid_pip) > 3.){continue;}
                  Double_t beta_pi  = pips[i]->par()->getBeta();
                  P_can = pips[i]->par()->getP();
                  if(P_can > 0.79 && P_can < 0.81)
                    {h_beta_P8_pip->Fill(beta);}
                  else if(P_can > 0.99 && P_can < 1.01)
                    {h_beta_P10_pip->Fill(beta);}
                  ///Remove CD, FTOF2, and sector 4 pions (getSector returns actual sector number)
                  if(pips[i]->par()->getStatus() < 2000 || pips[i]->par()->getStatus() > 3999){continue;}
                  else if(pips[i]->sci(FTOF2)->getEnergy()>0.){continue;}
                  else if(abs(pips[i]->getChi2Pid()-chi2pid_pip) > 3.){continue;}
                  //else if(beta_pi < 0.96+0.015*pips[i]->par()->getP() || beta_pi > ExpFunc(1.008, 0.1, 1, 4.5, euler, P_can)){continue;}
                  Vperp_pi = sqrt(pow(pips[i]->par()->getVx(),2) + pow(pips[i]->par()->getVy(),2));
                  Vdiff_pi = vz_e - pips[i]->par()->getVz();
                  ///Pion perpendicular and electron vertex difference cuts
                  if(PiCutter(0, P_can, 0.5, 1.0, Vperp_pi, Vdiff_pi, Vperp_pipcut, Vperp_pip_ms, Vperp_pimcut, Vperp_pim_ms, PiCutSig)!=0) continue;
                  ///Pion DC fiducial cuts
                  bool pip_fid_lay_1 = dcFid.DC_pi_fid(pips[i]->traj(DC,6)->getX(), pips[i]->traj(DC,6)->getY(), 3, pips[i]->trk(DC)->getSector(), 1,bending);
                  bool pip_fid_lay_2 = dcFid.DC_pi_fid(pips[i]->traj(DC,18)->getX(),pips[i]->traj(DC,18)->getY(),3, pips[i]->trk(DC)->getSector(), 2,bending);
                  bool pip_fid_lay_3 = dcFid.DC_pi_fid(pips[i]->traj(DC,36)->getX(),pips[i]->traj(DC,36)->getY(),3, pips[i]->trk(DC)->getSector(), 3,bending);
                  if( !(pip_fid_lay_1 && pip_fid_lay_2 && pip_fid_lay_3) ) continue;
                  ///Get maximum momentum pion
                  if(P_pi < P_can && pips[i]->par()->getCharge() == 1){P_pi = P_can; i_pi = i; iCharge = 0; beta = beta_pi;}
                }
            }
            if(pims.size()>=1)
            {
                for(Int_t i=0;i<pims.size();i++)
                {
                  if(abs(pims[i]->getChi2Pid()-chi2pid_pim) > 3.){continue;}
                  Double_t beta_pi  = pims[i]->par()->getBeta();
                  P_can = pims[i]->par()->getP();
                  if(P_can > 0.79 && P_can < 0.81)
                    {h_beta_P8_pim->Fill(beta);}
                  else if(P_can > 0.99 && P_can < 1.01)
                    {h_beta_P10_pim->Fill(beta);}
                  ///Remove CD, FTOF2, and sector 4 pions (getSector returns actual sector number)
                  if(pims[i]->par()->getStatus() < 2000 || pims[i]->par()->getStatus() > 3999){continue;}
                  else if(pims[i]->sci(FTOF2)->getEnergy()>0.){continue;}
                  else if(abs(pims[i]->getChi2Pid()-chi2pid_pim) > 3.){continue;}
                  //else if(beta_pi < 0.96+0.015*pims[i]->par()->getP() || beta_pi > 1.){continue;}
                  //else if(abs(pims[i]->getChi2Pid()) > 2.){continue;}
                  Vperp_pi = sqrt(pow(pims[i]->par()->getVx(),2) + pow(pims[i]->par()->getVy(),2));
                  Vdiff_pi = vz_e - pims[i]->par()->getVz();
                  ///Pion perpendicular and electron vertex difference cuts
                  if(PiCutter(1, P_can, 0.5, 1.0, Vperp_pi, Vdiff_pi, Vperp_pipcut, Vperp_pip_ms, Vperp_pimcut, Vperp_pim_ms, PiCutSig)!=0) continue;
                  ///Pion DC fiducial cuts
                  bool pim_fid_lay_1 = dcFid.DC_pi_fid(pims[i]->traj(DC,6)->getX(), pims[i]->traj(DC,6)->getY(), 4, pims[i]->trk(DC)->getSector(), 1,bending);
                  bool pim_fid_lay_2 = dcFid.DC_pi_fid(pims[i]->traj(DC,18)->getX(),pims[i]->traj(DC,18)->getY(),4, pims[i]->trk(DC)->getSector(), 2,bending);
                  bool pim_fid_lay_3 = dcFid.DC_pi_fid(pims[i]->traj(DC,36)->getX(),pims[i]->traj(DC,36)->getY(),4, pims[i]->trk(DC)->getSector(), 3,bending);
                  if( !(pim_fid_lay_1 && pim_fid_lay_2 && pim_fid_lay_3) ) continue;
                  ///Get maximum momentum pion
                  if(P_pi < P_can && pims[i]->par()->getCharge() == -1){P_pi = P_can; i_pi = i; iCharge = 1; beta = beta_pi;}
                }
            }
            if(protons.size()>=1)// && abs(protons[0]->getChi2Pid()+0.061)<3.)
            {
                for(Int_t i=0;i<protons.size();i++)
                {
                    Double_t P_p = protons[i]->par()->getP(), beta_p = protons[i]->par()->getBeta();
                    if(P_p == 0.) continue;
                    h_Wp_data->Fill(W);
                    h2_beta_P_p->Fill(P_p,beta_p);
                }
            }
            ///Remove events with no acceptable pions
            if(P_pi==0.0) continue;
            ///Pion vertex calculation
            if(iCharge==0)
            {
                Vperp_pi = sqrt(pow(pips[i_pi]->par()->getVx(),2) + pow(pips[i_pi]->par()->getVy(),2));
                Vdiff_pi = vz_e - pips[i_pi]->par()->getVz();
                Theta_pi = pips[i_pi]->getTheta();
                Phi_pi   = pips[i_pi]->getPhi();
                Chi2pid_pi=pips[i_pi]->getChi2Pid();
                i_FDstat = pips[i_pi]->par()->getStatus() - 2000;
                  if(P_pi > 0.79 && P_pi < 0.81)
                    {h_beta_P8_pipcut->Fill(beta);}
                  else if(P_pi > 0.99 && P_pi < 1.01)
                    {h_beta_P10_pipcut->Fill(beta);}
            }
            else if(iCharge==1)
            {
                Vperp_pi = sqrt(pow(pims[i_pi]->par()->getVx(),2) + pow(pims[i_pi]->par()->getVy(),2));
                Vdiff_pi = vz_e - pims[i_pi]->par()->getVz();
                Theta_pi = pims[i_pi]->getTheta();
                Phi_pi   = pims[i_pi]->getPhi();
                Chi2pid_pi=pims[i_pi]->getChi2Pid();
                i_FDstat = pims[i_pi]->par()->getStatus() - 2000;
                // Recalculating for neutron mass
                x_b      = q2/(2 * mass_n * v4_q.E() );
                y        = -v4_q.P() + sqrt( v4_q.E()*v4_q.E() + 2*v4_q.E()*mass_n);
                W        = sqrt(mass_n*mass_n - q2 + 2*v4_q.E()*mass_n);
                  if(P_pi > 0.79 && P_pi < 0.81)
                    {h_beta_P8_pimcut->Fill(beta);}
                  else if(P_pi > 0.99 && P_pi < 1.01)
                    {h_beta_P10_pimcut->Fill(beta);}
            }
            Int_t i_SCINstat = i_FDstat/100;
            Int_t i_ECALstat = (i_FDstat-(100*i_SCINstat))/10;
            //if(i_ECALstat == 0) continue;
            vec_pi.SetMagThetaPhi(P_pi,Theta_pi,Phi_pi);
            theta_piq = vec_q.Angle(vec_pi)*TMath::RadToDeg();
            Theta_pi *= TMath::RadToDeg();

          h2_theta_mom_data->Fill(v4_el.P(),v4_el.Theta()*TMath::RadToDeg());
          h_q2_data->Fill(q2);
          if(iCharge==0)      h_q2pip_data->Fill(q2);
          else if(iCharge==1) h_q2pim_data->Fill(q2);
          if(P_pi>0.)
          {
              h_Ppi_data->Fill(P_pi);
              if(iCharge==0)      h_Ppip_data->Fill(P_pi);
              else if(iCharge==1) h_Ppim_data->Fill(P_pi);
          }
          if(theta_piq>0.)
          {
              h_theta_piq_data->Fill(theta_piq);
              if(iCharge==0)      h_theta_pipq_data->Fill(theta_piq);
              else if(iCharge==1) h_theta_pimq_data->Fill(theta_piq);
          }
          if(Theta_pi>0.)
          {
              h_theta_pi_data->Fill(Theta_pi);
              if(iCharge==0)      h_theta_pip_data->Fill(Theta_pi);
              else if(iCharge==1) h_theta_pim_data->Fill(Theta_pi);
          }
          h2_theta_q2_data->Fill(q2,v4_el.Theta()*TMath::RadToDeg());
          h_W_data->Fill(W);
          if(iCharge==0)      h_Wpip_data->Fill(W);
          else if(iCharge==1) h_Wpim_data->Fill(W);
          h2_q2_W_data->Fill(W,q2);
          h_omega_data->Fill(v4_q.E());
          h2_theta_omega_data->Fill(v4_q.E(),v4_el.Theta()*TMath::RadToDeg());
          h2_theta_W_data->Fill(W,v4_el.Theta()*TMath::RadToDeg());
          if(iCharge==0)     {
                h2_DC_y_vs_x_data->Fill(pips[i_pi]->traj(DC,6)->getX(),pips[i_pi]->traj(DC,6)->getY());
                h_pip_chi2pid->Fill(Chi2pid_pi);
                h2_q2_Wpip_data->Fill(W,q2);
                h2_Ppip_thetapip_data->Fill(Theta_pi,P_pi);
                h2_Ppip_thetapipq_data->Fill(theta_piq,P_pi);
                h2_beta_P_pip->Fill(P_pi,beta);
                h2_beta_P_pipzoom->Fill(P_pi,beta);
                if(v4_q.Theta()*TMath::RadToDeg()>50.){
                    h2_beta_P_pipHiTQ->Fill(P_pi,beta);
                }
                if(W > 0.8 && W < 1.1){h2_beta_P_pipLowW->Fill(P_pi,beta);}
                if(pips[i_pi]->sci(FTOF2)->getEnergy()>0.){
                    h2_Ppip_thetapip_scin2->Fill(Theta_pi,P_pi);
                }
                else{
                    h2_Ppip_thetapip_scin1->Fill(Theta_pi,P_pi);
                }
                if(abs(pips[i_pi]->getChi2Pid()+0.061) < 2.)
                {
                    h_Wpip_data_chi2->Fill(W);
                    h2_beta_P_pipchi2->Fill(P_pi,beta);
                }
                else if(abs(pips[i_pi]->getChi2Pid()+0.061) < 3.)
                {
                    h_Wpip_data_chi3->Fill(W);
                    h2_beta_P_pipchi3->Fill(P_pi,beta);
                }
                if(beta < 0.5031*P_pi + 0.4706 + 0.0086 && beta > 0.5031*P_pi + 0.4706 - 0.0086)
                {
                    h_Wpip_data_ray->Fill(W);
                    h_Betapip_data_ray->Fill(beta);
                    h2_beta_P_pip_ray->Fill(P_pi,beta);
                }
                if(pips[i_pi]->sci(FTOF1B)->getEnergy()>0.){
                    if(W > 0.8 && W < 1.1)
                    {
                        h2_phi_e_phi_pip->Fill(Phi_pi*TMath::RadToDeg(),v4_q.Phi()*TMath::RadToDeg());
                        h2_E_FTOF_P_pipQE->Fill(P_pi,pips[i_pi]->sci(FTOF1B)->getEnergy());
                        h2_q2_Wpip_QE->Fill(W,q2);
                        h_ThetaQpip_dataQE->Fill(v4_q.Theta()*TMath::RadToDeg());
                    }
                    else{h_ThetaQpip_dataNoQE->Fill(v4_q.Theta()*TMath::RadToDeg());}
                    h2_beta_t_FTOF->Fill(pips[i_pi]->sci(FTOF1B)->getTime(),beta);
                    h2_E_FTOF_P_pip->Fill(P_pi,pips[i_pi]->sci(FTOF1B)->getEnergy());
                    h2_E_FTOF_P_pipz->Fill(P_pi,pips[i_pi]->sci(FTOF1B)->getEnergy());
                    h_W_data_pipFTOF1B->Fill(W);
                    h_ThetaQpip_dataFTOF1B->Fill(v4_q.Theta()*TMath::RadToDeg());
                }
                else if(pips[i_pi]->sci(FTOF1A)->getEnergy()>0.){
                    h2_beta_t_FTOFa->Fill(pips[i_pi]->sci(FTOF1A)->getTime(),beta);
                    h2_E_FTOFa_P_pip->Fill(P_pi,pips[i_pi]->sci(FTOF1A)->getEnergy());
                    h_W_data_pipFTOF1A->Fill(W);

                }
                else{
                    h_W_data_pipFTOFnoE->Fill(W);
                }
                if(v4_q.Theta()*TMath::RadToDeg()>50.){
                    h_Wpip_data_hiTQ->Fill(W);
                }
          }
          else if(iCharge==1){
                h_pim_chi2pid->Fill(Chi2pid_pi);
                h2_q2_Wpim_data->Fill(W,q2);
                h2_Ppim_thetapim_data->Fill(Theta_pi,P_pi);
                h2_Ppim_thetapimq_data->Fill(theta_piq,P_pi);
                h2_beta_P_pim->Fill(P_pi,beta);
                if(W > 0.9 && W < 1.0)
                {
                    h2_phi_e_phi_pim->Fill(Phi_pi*TMath::RadToDeg(),v4_q.Phi()*TMath::RadToDeg());
                }
                h2_E_FTOF_P_pim->Fill(P_pi,pims[i_pi]->sci(FTOF1B)->getEnergy());
          }


          //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ///Data CS calculations

            if(bQorNot == true){/// theta_piq
                EventBinner(bQorNot, iCharge, P_pi, theta_piq, q2, W, CS_Wmin, dW, Wbini, iPmom, iThPiQ, iQ2, iW);
            }
            else{/// changed theta_piq to Theta_pi
                EventBinner(bQorNot, iCharge, P_pi, Theta_pi, q2, W, CS_Wmin, dW, Wbini, iPmom, iThPiQ, iQ2, iW);
            }
            if(iPmom == -1 || iThPiQ == -1 || iQ2 == -1 || iW == -1 || iCharge == -1) continue;
            if(bMultiPion == true && pips.size() + pims.size() == 1) continue;
            //else if(iCharge == 0){cs_data++;}
          if(iCharge==0)     {
                h2_q2_Wpip_cut->Fill(W,q2);
                h2_Ppip_thetapip_cut->Fill(Theta_pi,P_pi);
                h2_Ppip_thetapipq_cut->Fill(theta_piq,P_pi);
                if(pips[i_pi]->trk(DC)->getSector() == 1)
                    h2_theta_W_S1_data->Fill(W,theta_piq);
                else if(pips[i_pi]->trk(DC)->getSector() == 2)
                    h2_theta_W_S2_data->Fill(W,theta_piq);
                else if(pips[i_pi]->trk(DC)->getSector() == 3)
                    h2_theta_W_S3_data->Fill(W,theta_piq);
                else if(pips[i_pi]->trk(DC)->getSector() == 4)
                    h2_theta_W_S4_data->Fill(W,theta_piq);
                else if(pips[i_pi]->trk(DC)->getSector() == 5)
                    h2_theta_W_S5_data->Fill(W,theta_piq);
                else if(pips[i_pi]->trk(DC)->getSector() == 6)
                    h2_theta_W_S6_data->Fill(W,theta_piq);
          }
          else if(iCharge==1){
                h2_q2_Wpim_cut->Fill(W,q2);
                h2_Ppim_thetapim_cut->Fill(Theta_pi,P_pi);
                h2_Ppim_thetapimq_cut->Fill(theta_piq,P_pi);
                if(pims[i_pi]->trk(DC)->getSector() == 1)
                    h2_theta_W_S1_MC->Fill(W,theta_piq);
                else if(pims[i_pi]->trk(DC)->getSector() == 2)
                    h2_theta_W_S2_MC->Fill(W,theta_piq);
                else if(pims[i_pi]->trk(DC)->getSector() == 3)
                    h2_theta_W_S3_MC->Fill(W,theta_piq);
                else if(pims[i_pi]->trk(DC)->getSector() == 4)
                    h2_theta_W_S4_MC->Fill(W,theta_piq);
                else if(pims[i_pi]->trk(DC)->getSector() == 5)
                    h2_theta_W_S5_MC->Fill(W,theta_piq);
                else if(pims[i_pi]->trk(DC)->getSector() == 6)
                    h2_theta_W_S6_MC->Fill(W,theta_piq);
          }

            //Fill events
            sec_i = electrons[0]->trk(DC)->getSector() - 1;// getSector returns actual sector number
            if(sec_i > 3) sec_i--;
            if(sec_i > -1 && sec_i < sectors)
            {
                y_i[iPmom][iThPiQ][iCharge][sec_i][iDATA][iQ2][iW] += 1.;
                counter_data++;
            }
	      }//if(electrons.size()==1)[line 590ish]
	      //if(bAllEvents==false && counter_data>100000) break;
    }///DATA while loop ************************************************************************************************************************************

  h_q2_data->SetLineColor(kBlack); q2_stacked->Add(h_q2_data); hs_q2_data->Add(h_q2_data);//h_q2_data->Clear();
  h_q2pip_data->SetLineColor(kRed); hs_q2_data->Add(h_q2pip_data);
  h_q2pim_data->SetLineColor(kBlue); hs_q2_data->Add(h_q2pim_data);
  h_el_vz_h_data->SetLineColor(kBlack); el_vz_h_stacked->Add(h_el_vz_h_data); //h_el_vz_h_data->Clear();
  h_W_data->SetLineColor(kBlack); hs_W_data->Add(h_W_data); //h_W_data->Clear();
  h_Wpip_data->SetLineColor(kRed); hs_W_data->Add(h_Wpip_data);
  h_Wpim_data->SetLineColor(kBlue); hs_W_data->Add(h_Wpim_data);
  h_omega_data->SetLineColor(kBlack); omega_stacked->Add(h_omega_data);
  h_p_proton_data->SetLineColor(kBlack);
  h_Ppi_data->SetLineColor(kBlack); hs_Ppi_data->Add(h_Ppi_data);
  h_Ppip_data->SetLineColor(kRed); hs_Ppi_data->Add(h_Ppip_data);
  h_Ppim_data->SetLineColor(kBlue); hs_Ppi_data->Add(h_Ppim_data);
  h_theta_pi_data->SetLineColor(kBlack); hs_theta_pi_data->Add(h_theta_pi_data);
  h_theta_pip_data->SetLineColor(kRed); hs_theta_pi_data->Add(h_theta_pip_data);
  h_theta_pim_data->SetLineColor(kBlue); hs_theta_pi_data->Add(h_theta_pim_data);
  h_theta_piq_data->SetLineColor(kBlack); hs_theta_piq_data->Add(h_theta_piq_data);
  h_theta_pipq_data->SetLineColor(kRed); hs_theta_piq_data->Add(h_theta_pipq_data);
  h_theta_pimq_data->SetLineColor(kBlue); hs_theta_piq_data->Add(h_theta_pimq_data);
 }

  //////////////////////////////////////////////////////////////////////////////
  //         Monte Carlo (GENIE) Reconstructed
  //////////////////////////////////////////////////////////////////////////////
  while(chain_MC.Next())
    {//cout << endl << "Start ";
      v4_beam_MC.SetE(beamE);
      v4_beam_MC.SetPz(beamE);

      if( beamE > 1e-4)
	      {
	        //if( counter_MC == 0)
	        //cout<<"Beam energy set manually: "<<beamE<<endl;
	        v4_beam_MC.SetE(beamE);
	        v4_beam_MC.SetPz(beamE);
	      }
      counter_MC1++;

      TLorentzVector v4_gen_e, v4_gen_pip, v4_gen_pim;
      Int_t Ngen=c12_gen->mcparts()->getRows();
      Int_t iPmom_gen = -1, iThPiQ_gen = -1, iCharge_gen = -1, iFile_gen = -1, iQ2_gen = -1, iW_gen = -1;
      P_pi = 0.0; Double_t P_can = 0.0;
      Double_t mass_gen = mass_p;
      //cout << counter_gen1 << endl;
      for(Int_t i_gen=0; i_gen<Ngen; i_gen++)
      {
          c12_gen->mcparts()->setEntry(i_gen);// Need for gen?
          auto pax=c12_gen->mcparts()->getPx();
          auto pay=c12_gen->mcparts()->getPy();
          auto paz=c12_gen->mcparts()->getPz();
          auto pam=c12_gen->mcparts()->getMass();
          auto pid=c12_gen->mcparts()->getPid();
          //cout << i_gen <<" "<< pid;

          if(pid == 11 && i_gen==0)
          {
              v4_gen_e.SetXYZM(pax, pay, paz, pam);
          }

          ///Highest momentum Pion Calculator
          if(pid == 211)
          {
            P_can = c12_gen->mcparts()->getP();
            if(P_pi < P_can)
            {
                P_pi = P_can;
                iCharge_gen = 0;
                mass_gen = mass_p;
                Theta_pi = c12_gen->mcparts()->getTheta();
                Phi_pi   = c12_gen->mcparts()->getPhi();
                //cout <<" "<< v4_gen_pip.P() <<"="<< P_can;
            }
          }
          if(pid == -211)
          {//sqrt(pax*pax + pay*pay + paz*paz)
            P_can = c12_gen->mcparts()->getP();
            if(P_pi < P_can)
            {
                P_pi = P_can;
                iCharge_gen = 1;
                mass_gen = mass_n;
                Theta_pi = c12_gen->mcparts()->getTheta();
                Phi_pi   = c12_gen->mcparts()->getPhi();
                //cout <<" "<< v4_gen_pim.P() <<"="<< P_can;
            }
          }
          //cout <<" "<< P_pi << "\t";
      }//cout << " Ngen: " << Ngen << " ";

      TLorentzVector v4_gen_q = v4_beam_MC - v4_gen_e; //photon 4-vector

      //Calculate some important kinematic variables
      Double_t q2_gen    = -v4_gen_q.M2();
      Double_t W_gen     = sqrt(mass_gen*mass_gen - q2_gen + 2*v4_gen_q.E()*mass_gen); //Invariant mass, assuming electron off of a standing proton
      ///Theta_piq angle calculation
      vec_q.SetMagThetaPhi(v4_gen_q.P(),v4_gen_q.Theta(),v4_gen_q.Phi());
      vec_pi.SetMagThetaPhi(P_pi,Theta_pi,Phi_pi);
      theta_piq = vec_q.Angle(vec_pi)*TMath::RadToDeg();

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///MC_gen CS calculations

        if(bQorNot == true){/// theta_piq
            EventBinner(bQorNot, iCharge_gen, P_pi, theta_piq, q2_gen, W_gen, CS_Wmin, dW, Wbini, iPmom_gen, iThPiQ_gen, iQ2_gen, iW_gen);
        }
        /*else{
            EventBinner(bQorNot, iCharge_gen, P_pi, Theta_pi, q2_gen, W_gen, CS_Wmin, dW, Wbini, iPmom_gen, iThPiQ_gen, iQ2_gen, iW_gen);
        }*/
        if(iPmom_gen != -1 && iThPiQ_gen != -1 && iQ2_gen != -1 && iW_gen != -1){
            //Fill events
            y_i[iPmom_gen][iThPiQ_gen][iCharge_gen][0][iGEN][iQ2_gen][iW_gen] += 1.;//sector is set to 0
            counter_gen++;
        }
        Theta_pi = 0.0; Phi_pi = 0.0; //cout << iPmom_gen << iThPiQ_gen << iQ2_gen << iW_gen << "\t";

      // get particles by type
      auto electrons=c12_MC->getByID(11);
      auto protons=c12_MC->getByID(2212);
      auto neutrons=c12_MC->getByID(2112);
      auto pips=c12_MC->getByID(211);
      auto pims=c12_MC->getByID(-211);

      //cout << electrons.size() << pips.size() << pims.size();
      //Make sure that we have a good electron and that we are not seeing the weird delta function around zero in MC
      if(electrons.size()==1 && electrons[0]->trk(DC)->getSector() != 4 && (pips.size()>=1 || pims.size()>=1))///changed to single pi+ and/or pi-
	      {//cout << " good event\t";
          Double_t energy =  (electrons[0]->cal(PCAL)->getEnergy() +  electrons[0]->cal(ECIN)->getEnergy() +  electrons[0]->cal(ECOUT)->getEnergy())/electrons[0]->getP();

          bool ecal = ( (electrons[0]->cal(ECIN)->getLv() >= 14 && electrons[0]->cal(ECIN)->getLw() >= 14) && (energy > 0.18 && energy < 0.28) && electrons[0]->getP() > 0.5 && electrons[0]->getP() < 10 );

          SetLorentzVector(v4_el,electrons[0]);

          TLorentzVector v4_q = v4_beam_MC - v4_el; //photon 4-vector

          //Calculate some important kinematic variables
          Double_t q2        = -v4_q.M2();
          Double_t x_b       = q2/(2 * mass_p * v4_q.E() ); //Bjorken-x
          Double_t y         = -v4_q.P() + sqrt( v4_q.E()*v4_q.E() + 2*v4_q.E()*mass_p); //y (Bjorken-y), the scaling variable; from Jourdan, eq. 20
          Double_t vz_e      = electrons[0]->par()->getVz();
          Double_t W         = sqrt(mass_p*mass_p - q2 + 2*v4_q.E()*mass_p); //Invariant mass, assuming electron off of a standing proton

          if(W < Wcut) continue; /// W min cut (removes QE "pion" [perhaps coincident proton] events with high ThetaQ)
          else if(W > CS_Wmax) continue; /// Save time

          Double_t processid = c12_MC->mcevent()->getWeight();//code = 1,2,3,4 = type = qel,mec,res,dis
          Double_t resid = c12_MC->mcevent()->getPtarget();//targP = target polarization, which is not relevant = resid

          //std::cout << vz_e << std::endl;
          h_el_vz_h_MC->Fill(vz_e);

          //Electron quality cuts
          bool targ_1 = abs( vz_e - (-1.48)) < 1 * 1.25;
          bool targ_2 = abs( vz_e - (-6.3))  < 1 * 1.25;

          ///***********************************************************************************************************************************************************************************************
            Int_t iPmom = -1, iThPiQ = -1, iCharge = -1, iMod = -1, iQ2 = -1, iW = -1;
            ///Electron vertex calculation
            vec_q.SetMagThetaPhi(v4_q.P(),v4_q.Theta(),v4_q.Phi());
            Vperp_e = sqrt(pow(electrons[0]->par()->getVx(),2) + pow(electrons[0]->par()->getVy(),2));
            P_pi = 0.0; i_pi = 0; P_can = 0.0;// pion momentum, pion # for the event, momentum of pion candidate
            Double_t Lv = electrons[0]->cal(PCAL)->getLv(), Lw = electrons[0]->cal(PCAL)->getLw();
            Double_t beta = 0.;

            ///Electron EC fiducial and vertex cuts
            if(eCutter(VWcut, Lv, Lw, vz_e, vz_e_min, vz_e_max, Vperp_e, Vperp_ecut)!=0) continue;// cout <<"Passed: EC, ";
            //DC fiducial cuts for electron
            bool el_fid_lay_1 = dcFid.DC_pi_fid(electrons[0]->traj(DC,6)->getX(), electrons[0]->traj(DC,6)->getY(),  1, electrons[0]->trk(DC)->getSector(), 1,bending);
            bool el_fid_lay_2 = dcFid.DC_pi_fid(electrons[0]->traj(DC,18)->getX(),electrons[0]->traj(DC,18)->getY(), 1, electrons[0]->trk(DC)->getSector(), 2,bending);
            bool el_fid_lay_3 = dcFid.DC_pi_fid(electrons[0]->traj(DC,36)->getX(),electrons[0]->traj(DC,36)->getY(), 1, electrons[0]->trk(DC)->getSector(), 3,bending);
            if( !(el_fid_lay_1 && el_fid_lay_2 && el_fid_lay_3) ) continue;// cout << "DC, ";

            ///Highest momentum Pion Calculator
            if(pips.size()>=1)
            {
                for(Int_t i=0;i<pips.size();i++)
                {
                  Double_t beta_pi  = pips[i]->par()->getBeta();
                  P_can = pips[i]->par()->getP();
                  ///Remove CD, FTOF2, and high chi2pid pions (getSector returns actual sector number)
                  if(pips[i]->par()->getStatus() < 2000 || pips[i]->par()->getStatus() > 3999){continue;}
                  else if(pips[i]->sci(FTOF2)->getEnergy()>0.){continue;}
                  else if(abs(pips[i]->getChi2Pid()-chi2pid_pip) > 3.){continue;}
                  Vperp_pi = sqrt(pow(pips[i]->par()->getVx(),2) + pow(pips[i]->par()->getVy(),2));
                  Vdiff_pi = vz_e - pips[i]->par()->getVz();
                  ///Pion perpendicular and electron vertex difference cuts
                  if(PiCutter(0, P_can, 0.5, 1.0, Vperp_pi, Vdiff_pi, Vperp_pipcut, Vperp_pip_ms, Vperp_pimcut, Vperp_pim_ms, PiCutSig)!=0) continue;
                  ///Pion DC fiducial cuts
                  bool pip_fid_lay_1 = dcFid.DC_pi_fid(pips[i]->traj(DC,6)->getX(), pips[i]->traj(DC,6)->getY(), 3, pips[i]->trk(DC)->getSector(), 1,bending);
                  bool pip_fid_lay_2 = dcFid.DC_pi_fid(pips[i]->traj(DC,18)->getX(),pips[i]->traj(DC,18)->getY(),3, pips[i]->trk(DC)->getSector(), 2,bending);
                  bool pip_fid_lay_3 = dcFid.DC_pi_fid(pips[i]->traj(DC,36)->getX(),pips[i]->traj(DC,36)->getY(),3, pips[i]->trk(DC)->getSector(), 3,bending);
                  if( !(pip_fid_lay_1 && pip_fid_lay_2 && pip_fid_lay_3) ) continue;
                  ///Get maximum momentum pion
                  if(P_pi < P_can && pips[i]->par()->getCharge() == 1){P_pi = P_can; i_pi = i; iCharge = 0; beta = beta_pi;}
                }
            }
            if(pims.size()>=1)
            {
                for(Int_t i=0;i<pims.size();i++)
                {
                  Double_t beta_pi  = pims[i]->par()->getBeta();
                  ///Remove CD, FTOF2, and high chi2pid pions (getSector returns actual sector number)
                  if(pims[i]->par()->getStatus() < 2000 || pims[i]->par()->getStatus() > 3999){continue;}
                  else if(pims[i]->sci(FTOF2)->getEnergy()>0.){continue;}
                  else if(abs(pims[i]->getChi2Pid()-chi2pid_pim) > 3.){continue;}
                  Vperp_pi = sqrt(pow(pims[i]->par()->getVx(),2) + pow(pims[i]->par()->getVy(),2));
                  Vdiff_pi = vz_e - pims[i]->par()->getVz();
                  P_can = pims[i]->par()->getP();
                  ///Pion perpendicular and electron vertex difference cuts
                  if(PiCutter(1, P_can, 0.5, 1.0, Vperp_pi, Vdiff_pi, Vperp_pipcut, Vperp_pip_ms, Vperp_pimcut, Vperp_pim_ms, PiCutSig)!=0) continue;
                  ///Pion DC fiducial cuts
                  bool pim_fid_lay_1 = dcFid.DC_pi_fid(pims[i]->traj(DC,6)->getX(), pims[i]->traj(DC,6)->getY(), 4, pims[i]->trk(DC)->getSector(), 1,bending);
                  bool pim_fid_lay_2 = dcFid.DC_pi_fid(pims[i]->traj(DC,18)->getX(),pims[i]->traj(DC,18)->getY(),4, pims[i]->trk(DC)->getSector(), 2,bending);
                  bool pim_fid_lay_3 = dcFid.DC_pi_fid(pims[i]->traj(DC,36)->getX(),pims[i]->traj(DC,36)->getY(),4, pims[i]->trk(DC)->getSector(), 3,bending);
                  if( !(pim_fid_lay_1 && pim_fid_lay_2 && pim_fid_lay_3) ) continue;
                  ///Get maximum momentum pion
                  if(P_pi < P_can && pims[i]->par()->getCharge() == -1){P_pi = P_can; i_pi = i; iCharge = 1; beta = beta_pi;}
                }
            }
            ///Remove events with no acceptable pions
            if(P_pi==0.0) continue;// cout << "accept, ";
            ///Pion vertex calculation
            if(iCharge==0)
            {
                Vperp_pi = sqrt(pow(pips[i_pi]->par()->getVx(),2) + pow(pips[i_pi]->par()->getVy(),2));
                Vdiff_pi = vz_e - pips[i_pi]->par()->getVz();
                Theta_pi = pips[i_pi]->getTheta();
                Phi_pi   = pips[i_pi]->getPhi(); i_pip++;
            }
            else if(iCharge==1)
            {
                Vperp_pi = sqrt(pow(pims[i_pi]->par()->getVx(),2) + pow(pims[i_pi]->par()->getVy(),2));
                Vdiff_pi = vz_e - pims[i_pi]->par()->getVz();
                Theta_pi = pims[i_pi]->getTheta();
                Phi_pi   = pims[i_pi]->getPhi(); i_pim++;
                x_b      = q2/(2 * mass_n * v4_q.E() );
                y        = -v4_q.P() + sqrt( v4_q.E()*v4_q.E() + 2*v4_q.E()*mass_n);
                W        = sqrt(mass_n*mass_n - q2 + 2*v4_q.E()*mass_n);
            }
            vec_pi.SetMagThetaPhi(P_pi,Theta_pi,Phi_pi);
            theta_piq = vec_q.Angle(vec_pi)*TMath::RadToDeg();
            Theta_pi *= TMath::RadToDeg();
            if(iCharge==0){
                h2_Ppip_thetapip_MC->Fill(Theta_pi,P_pi);
            }
            else if(iCharge==1){
                h2_Ppim_thetapim_MC->Fill(Theta_pi,P_pi);
            }

          //Fill various histograms
          h2_theta_mom_MC->Fill(v4_el.P(),v4_el.Theta()*TMath::RadToDeg());
          h_q2_MC->Fill(q2);
          if (processid==1){h_q2_qe_MC->Fill(q2);}
          else if (processid==2){h_q2_mec_MC->Fill(q2);}
          else if (processid==3){h_q2_res_MC->Fill(q2);}
          else if (processid==4){h_q2_dis_MC->Fill(q2);}
          if(P_pi>0.){h_Ppi_MC->Fill(P_pi);}
          if(theta_piq>0.){h_theta_piq_MC->Fill(theta_piq);}
          if(Theta_pi>0.){h_theta_pi_MC->Fill(Theta_pi);}
          h2_theta_q2_MC->Fill(q2,v4_el.Theta()*TMath::RadToDeg());
          h_W_MC->Fill(W);
          if (processid==1){h_W_qe_MC->Fill(W);}
          else if (processid==2){h_W_mec_MC->Fill(W);}
          else if (processid==3){
                h_W_res_MC->Fill(W);
                h_W_resid_MC[Int_t(resid)]->Fill(W);
          }
          else if (processid==4){h_W_dis_MC->Fill(W);}
          h2_q2_W_MC->Fill(W,q2);
          h_omega_MC->Fill(v4_q.E());
          if (processid==1){h_omega_qe_MC->Fill(v4_q.E());}
          else if (processid==2){h_omega_mec_MC->Fill(v4_q.E());}
          else if (processid==3){h_omega_res_MC->Fill(v4_q.E());}
          else if (processid==4){h_omega_dis_MC->Fill(v4_q.E());}
          h2_theta_omega_MC->Fill(v4_q.E(),v4_el.Theta()*TMath::RadToDeg());
          h2_theta_W_MC->Fill(W,v4_el.Theta()*TMath::RadToDeg());
          if(iCharge==0)      h2_Ppip_theta_MC->Fill(v4_el.Theta()*TMath::RadToDeg(),P_pi);
          else if(iCharge==1) h2_Ppim_theta_MC->Fill(v4_el.Theta()*TMath::RadToDeg(),P_pi);


          //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ///MC CS calculations

            if(bQorNot == true){/// theta_piq
                EventBinner(bQorNot, iCharge, P_pi, theta_piq, q2, W, CS_Wmin, dW, Wbini, iPmom, iThPiQ, iQ2, iW);
            }
            else{/// changed theta_piq to Theta_pi
                EventBinner(bQorNot, iCharge, P_pi, Theta_pi, q2, W, CS_Wmin, dW, Wbini, iPmom, iThPiQ, iQ2, iW);
            }
            if(iPmom == -1 || iThPiQ == -1 || iQ2 == -1 || iW == -1) continue;// cout << "binned";
            if(bMultiPion == true && pips.size() + pims.size() == 1) continue;
            //else if(iCharge == 0){cs_MC++; continue;}

            //Fill events
            sec_i = electrons[0]->trk(DC)->getSector() - 1;// getSector returns actual sector number
            if(sec_i > 3) sec_i--;
            if(sec_i > -1 && sec_i < sectors)
            {
                y_i[iPmom][iThPiQ][iCharge][sec_i][iGENIE][iQ2][iW] += 1.;
                counter_MC++;
            }
	        //cout << y_i[iPmom][iThPiQ][iCharge][sec_i][iGENIE][iQ2][iW] << "\t" << counter_MC<<endl;
	      }//if(electrons.size()==1 && (electrons[0]->par()->getVz() < -0.00001 || electrons[0]->par()->getVz() > 0.00001) )
	      if(bAllEvents==false && counter_MC>MC_evts[iGENmod]) break;
    }///GENIE while loop ************************************************************************************************************************************

  //Fudge some scaling
  Double_t MC_to_Data_Scaling = 1.;
  if(bData==true){MC_to_Data_Scaling = ( (Double_t)counter_data ) / ( (Double_t)counter_MC );}
  //Scale, color, and place all histograms
  //Main histograms
  h_q2_MC->SetMarkerColor(kRed); h_q2_MC->Sumw2(); h_q2_MC->Scale(MC_to_Data_Scaling); h_q2_MC->SetLineColor(kRed); q2_stacked->Add(h_q2_MC); //h_q2_MC->Clear();
  h_q2_qe_MC->SetMarkerColor(kGreen); h_q2_qe_MC->Sumw2(); h_q2_qe_MC->Scale(MC_to_Data_Scaling); h_q2_qe_MC->SetLineColor(kGreen); q2_stacked->Add(h_q2_qe_MC);
  h_q2_mec_MC->SetMarkerColor(kOrange); h_q2_mec_MC->Sumw2(); h_q2_mec_MC->Scale(MC_to_Data_Scaling); h_q2_mec_MC->SetLineColor(kOrange); q2_stacked->Add(h_q2_mec_MC);
  h_q2_res_MC->SetMarkerColor(kBlue); h_q2_res_MC->Sumw2(); h_q2_res_MC->Scale(MC_to_Data_Scaling); h_q2_res_MC->SetLineColor(kBlue); q2_stacked->Add(h_q2_res_MC);
  h_q2_dis_MC->SetMarkerColor(kCyan); h_q2_dis_MC->Sumw2(); h_q2_dis_MC->Scale(MC_to_Data_Scaling); h_q2_dis_MC->SetLineColor(kCyan); q2_stacked->Add(h_q2_dis_MC);
  h_el_vz_h_MC->SetMarkerColor(kRed); h_el_vz_h_MC->Sumw2(); h_el_vz_h_MC->Scale(MC_to_Data_Scaling); h_el_vz_h_MC->SetLineColor(kRed); el_vz_h_stacked->Add(h_el_vz_h_MC); //h_el_vz_h_MC->Clear();
  h_W_MC->SetMarkerColor(kRed); h_W_MC->Sumw2(); h_W_MC->Scale(MC_to_Data_Scaling); h_W_MC->SetLineColor(kRed); W_stacked->Add(h_W_MC); //h_W_MC->Clear();
  h_W_qe_MC->SetMarkerColor(kGreen); h_W_qe_MC->Sumw2(); h_W_qe_MC->Scale(MC_to_Data_Scaling); h_W_qe_MC->SetLineColor(kGreen); W_stacked->Add(h_W_qe_MC);
  h_W_mec_MC->SetMarkerColor(kOrange); h_W_mec_MC->Sumw2(); h_W_mec_MC->Scale(MC_to_Data_Scaling); h_W_mec_MC->SetLineColor(kOrange); W_stacked->Add(h_W_mec_MC);
  h_W_res_MC->SetMarkerColor(kBlue); h_W_res_MC->Sumw2(); h_W_res_MC->Scale(MC_to_Data_Scaling); h_W_res_MC->SetLineColor(kBlue); W_stacked->Add(h_W_res_MC);
  h_W_dis_MC->SetMarkerColor(kCyan); h_W_dis_MC->Sumw2(); h_W_dis_MC->Scale(MC_to_Data_Scaling); h_W_dis_MC->SetLineColor(kCyan); W_stacked->Add(h_W_dis_MC);
  h_omega_MC->SetMarkerColor(kRed); h_omega_MC->Sumw2(); h_omega_MC->Scale(MC_to_Data_Scaling); h_omega_MC->SetLineColor(kRed); omega_stacked->Add(h_omega_MC);
  h_omega_qe_MC->SetMarkerColor(kGreen); h_omega_qe_MC->Sumw2(); h_omega_qe_MC->Scale(MC_to_Data_Scaling); h_omega_qe_MC->SetLineColor(kGreen); omega_stacked->Add(h_omega_qe_MC);
  h_omega_mec_MC->SetMarkerColor(kOrange); h_omega_mec_MC->Sumw2(); h_omega_mec_MC->Scale(MC_to_Data_Scaling); h_omega_mec_MC->SetLineColor(kOrange); omega_stacked->Add(h_omega_mec_MC);
  h_omega_res_MC->SetMarkerColor(kBlue); h_omega_res_MC->Sumw2(); h_omega_res_MC->Scale(MC_to_Data_Scaling); h_omega_res_MC->SetLineColor(kBlue); omega_stacked->Add(h_omega_res_MC);
  h_omega_dis_MC->SetMarkerColor(kCyan); h_omega_dis_MC->Sumw2(); h_omega_dis_MC->Scale(MC_to_Data_Scaling); h_omega_dis_MC->SetLineColor(kCyan); omega_stacked->Add(h_omega_dis_MC);
  h_p_proton_MC->SetMarkerColor(kRed); h_p_proton_MC->Sumw2(); h_p_proton_MC->Scale(MC_to_Data_Scaling); h_p_proton_MC->SetLineColor(kRed);
  hs_W_resid_MC->Add(h_W_res_MC);
  for(Int_t id=0;id<20;id++)
  {
    Int_t col = id+1;
    if(id > 3){col = id + 2;}
    else if(id > 9){col = id + 3;}
    h_W_resid_MC[id]->SetMarkerColor(col); h_W_resid_MC[id]->Sumw2(); h_W_resid_MC[id]->Scale(MC_to_Data_Scaling); h_W_resid_MC[id]->SetLineColor(col); hs_W_resid_MC->Add(h_W_resid_MC[id]);
  }

  //////////////////////////////////////////////////////////////////////////////
  //         Monte Carlo (GENIE) Radiative
  //////////////////////////////////////////////////////////////////////////////
  while(chain_MCrad.Next())
    {
      v4_beam_MCrad.SetE(beamE);
      v4_beam_MCrad.SetPz(beamE);

      if( beamE > 1e-4)
	      {
	        //if( counter_MCrad == 0)
	        //cout<<"Beam energy set manually: "<<beamE<<endl;
	        v4_beam_MCrad.SetE(beamE);
	        v4_beam_MCrad.SetPz(beamE);
	      }
      counter_MCrad++;

      // get particles by type
      auto electrons=c12_MCrad->getByID(11);
      auto protons=c12_MCrad->getByID(2212);
      auto neutrons=c12_MCrad->getByID(2112);
      auto pips=c12_MCrad->getByID(211);
      auto pims=c12_MCrad->getByID(-211);

      //Make sure that we have a good electron and that we are not seeing the weird delta function around zero in MC
      if(electrons.size()==1 && electrons[0]->trk(DC)->getSector() != 4 && (pips.size()>=1 || pims.size()>=1))///changed to single pi+ and/or pi-
	      {
          SetLorentzVector(v4_el,electrons[0]);

          TLorentzVector v4_q = v4_beam_MCrad - v4_el; //photon 4-vector

          //Calculate some important kinematic variables
          Double_t q2        = -v4_q.M2();
          Double_t x_b       = q2/(2 * mass_p * v4_q.E() ); //Bjorken-x
          Double_t y         = -v4_q.P() + sqrt( v4_q.E()*v4_q.E() + 2*v4_q.E()*mass_p); //y (Bjorken-y), the scaling variable; from Jourdan, eq. 20
          Double_t vz_e      = electrons[0]->par()->getVz();
          Double_t W         = sqrt(mass_p*mass_p - q2 + 2*v4_q.E()*mass_p); //Invariant mass, assuming electron off of a standing proton

          if(W < Wcut) continue; /// W min cut (removes QE "pion" [perhaps coincident proton] events with high ThetaQ)
          else if(W > CS_Wmax) continue; /// Save time

          Double_t MCrad_weight = c12_MCrad->mcevent()->getWeight();//NOT PROCESSID OF OTHER GENIE FILES
          auto processid_true = c12_MCrad->mcevent()->getProcessid();//code = 1,2,3,4 = type = qel,mec,res,dis???????

          ///***********************************************************************************************************************************************************************************************
            Int_t iPmom = -1, iThPiQ = -1, iCharge = -1, iMod = -1, iQ2 = -1, iW = -1;
            ///Electron vertex calculation
            vec_q.SetMagThetaPhi(v4_q.P(),v4_q.Theta(),v4_q.Phi());
            Vperp_e = sqrt(pow(electrons[0]->par()->getVx(),2) + pow(electrons[0]->par()->getVy(),2));
            P_pi = 0.0; i_pi = 0; Double_t P_can = 0.0;// pion momentum, pion # for the event, momentum of pion candidate
            Double_t Lv = electrons[0]->cal(PCAL)->getLv(), Lw = electrons[0]->cal(PCAL)->getLw();
            Double_t beta = 0.;

            ///Electron EC fiducial and vertex cuts
            if(eCutter(VWcut, Lv, Lw, vz_e, vz_e_min, vz_e_max, Vperp_e, Vperp_ecut)!=0) continue;
            //DC fiducial cuts for electron
            bool el_fid_lay_1 = dcFid.DC_pi_fid(electrons[0]->traj(DC,6)->getX(), electrons[0]->traj(DC,6)->getY(),  1, electrons[0]->trk(DC)->getSector(), 1,bending);
            bool el_fid_lay_2 = dcFid.DC_pi_fid(electrons[0]->traj(DC,18)->getX(),electrons[0]->traj(DC,18)->getY(), 1, electrons[0]->trk(DC)->getSector(), 2,bending);
            bool el_fid_lay_3 = dcFid.DC_pi_fid(electrons[0]->traj(DC,36)->getX(),electrons[0]->traj(DC,36)->getY(), 1, electrons[0]->trk(DC)->getSector(), 3,bending);
            if( !(el_fid_lay_1 && el_fid_lay_2 && el_fid_lay_3) ) continue;

            ///Highest momentum Pion Calculator
            if(pips.size()>=1)
            {
                for(Int_t i=0;i<pips.size();i++)
                {
                  Double_t beta_pi  = pips[i]->par()->getBeta();
                  P_can = pips[i]->par()->getP();
                  ///Remove CD, FTOF2, and high chi2pid pions (getSector returns actual sector number)
                  if(pips[i]->par()->getStatus() < 2000 || pips[i]->par()->getStatus() > 3999){continue;}
                  else if(pips[i]->sci(FTOF2)->getEnergy()>0.){continue;}
                  else if(abs(pips[i]->getChi2Pid()-chi2pid_pip) > 3.){continue;}
                  Vperp_pi = sqrt(pow(pips[i]->par()->getVx(),2) + pow(pips[i]->par()->getVy(),2));
                  Vdiff_pi = vz_e - pips[i]->par()->getVz();
                  ///Pion perpendicular and electron vertex difference cuts
                  if(PiCutter(0, P_can, 0.5, 1.0, Vperp_pi, Vdiff_pi, Vperp_pipcut, Vperp_pip_ms, Vperp_pimcut, Vperp_pim_ms, PiCutSig)!=0) continue;
                  ///Pion DC fiducial cuts
                  bool pip_fid_lay_1 = dcFid.DC_pi_fid(pips[i]->traj(DC,6)->getX(), pips[i]->traj(DC,6)->getY(), 3, pips[i]->trk(DC)->getSector(), 1,bending);
                  bool pip_fid_lay_2 = dcFid.DC_pi_fid(pips[i]->traj(DC,18)->getX(),pips[i]->traj(DC,18)->getY(),3, pips[i]->trk(DC)->getSector(), 2,bending);
                  bool pip_fid_lay_3 = dcFid.DC_pi_fid(pips[i]->traj(DC,36)->getX(),pips[i]->traj(DC,36)->getY(),3, pips[i]->trk(DC)->getSector(), 3,bending);
                  if( !(pip_fid_lay_1 && pip_fid_lay_2 && pip_fid_lay_3) ) continue;
                  ///Get maximum momentum pion
                  if(P_pi < P_can && pips[i]->par()->getCharge() == 1){P_pi = P_can; i_pi = i; iCharge = 0; beta = beta_pi;}
                }
            }
            if(pims.size()>=1)
            {
                for(Int_t i=0;i<pims.size();i++)
                {
                  Double_t beta_pi  = pims[i]->par()->getBeta();
                  ///Remove CD, FTOF2, and high chi2pid pions (getSector returns actual sector number)
                  if(pims[i]->par()->getStatus() < 2000 || pims[i]->par()->getStatus() > 3999){continue;}
                  else if(pims[i]->sci(FTOF2)->getEnergy()>0.){continue;}
                  else if(abs(pims[i]->getChi2Pid()-chi2pid_pim) > 3.){continue;}
                  Vperp_pi = sqrt(pow(pims[i]->par()->getVx(),2) + pow(pims[i]->par()->getVy(),2));
                  Vdiff_pi = vz_e - pims[i]->par()->getVz();
                  P_can = pims[i]->par()->getP();
                  ///Pion perpendicular and electron vertex difference cuts
                  if(PiCutter(1, P_can, 0.5, 1.0, Vperp_pi, Vdiff_pi, Vperp_pipcut, Vperp_pip_ms, Vperp_pimcut, Vperp_pim_ms, PiCutSig)!=0) continue;
                  ///Pion DC fiducial cuts
                  bool pim_fid_lay_1 = dcFid.DC_pi_fid(pims[i]->traj(DC,6)->getX(), pims[i]->traj(DC,6)->getY(), 4, pims[i]->trk(DC)->getSector(), 1,bending);
                  bool pim_fid_lay_2 = dcFid.DC_pi_fid(pims[i]->traj(DC,18)->getX(),pims[i]->traj(DC,18)->getY(),4, pims[i]->trk(DC)->getSector(), 2,bending);
                  bool pim_fid_lay_3 = dcFid.DC_pi_fid(pims[i]->traj(DC,36)->getX(),pims[i]->traj(DC,36)->getY(),4, pims[i]->trk(DC)->getSector(), 3,bending);
                  if( !(pim_fid_lay_1 && pim_fid_lay_2 && pim_fid_lay_3) ) continue;
                  ///Get maximum momentum pion
                  if(P_pi < P_can && pims[i]->par()->getCharge() == -1){P_pi = P_can; i_pi = i; iCharge = 1; beta = beta_pi;}
                }
            }
            ///Remove events with no acceptable pions
            if(P_pi==0.0) continue;
            ///Pion vertex calculation
            if(iCharge==0)
            {
                Vperp_pi = sqrt(pow(pips[i_pi]->par()->getVx(),2) + pow(pips[i_pi]->par()->getVy(),2));
                Vdiff_pi = vz_e - pips[i_pi]->par()->getVz();
                Theta_pi = pips[i_pi]->getTheta();
                Phi_pi   = pips[i_pi]->getPhi(); i_pip++;
            }
            else if(iCharge==1)
            {
                Vperp_pi = sqrt(pow(pims[i_pi]->par()->getVx(),2) + pow(pims[i_pi]->par()->getVy(),2));
                Vdiff_pi = vz_e - pims[i_pi]->par()->getVz();
                Theta_pi = pims[i_pi]->getTheta();
                Phi_pi   = pims[i_pi]->getPhi(); i_pim++;
                x_b      = q2/(2 * mass_n * v4_q.E() );
                y        = -v4_q.P() + sqrt( v4_q.E()*v4_q.E() + 2*v4_q.E()*mass_n);
                W        = sqrt(mass_n*mass_n - q2 + 2*v4_q.E()*mass_n);
            }
            vec_pi.SetMagThetaPhi(P_pi,Theta_pi,Phi_pi);
            theta_piq = vec_q.Angle(vec_pi)*TMath::RadToDeg();
            Theta_pi *= TMath::RadToDeg();
            if(iCharge==0){
                h2_Ppip_thetapip_MC->Fill(Theta_pi,P_pi);
            }
            else if(iCharge==1){
                h2_Ppim_thetapim_MC->Fill(Theta_pi,P_pi);
            }

          //Fill various histograms
          h_q2_MCrad->Fill(q2);
          /*if (processid==1){h_q2_qe_MC->Fill(q2);}
          else if (processid==2){h_q2_mec_MC->Fill(q2);}
          else if (processid==3){h_q2_res_MC->Fill(q2);}
          else if (processid==4){h_q2_dis_MC->Fill(q2);}*/
          if(P_pi>0.){h_Ppi_MCrad->Fill(P_pi);}
          /*if(theta_piq>0.){h_theta_piq_MC->Fill(theta_piq);}
          if(Theta_pi>0.){h_theta_pi_MC->Fill(Theta_pi);}
          h2_theta_q2_MC->Fill(q2,v4_el.Theta()*TMath::RadToDeg());*/
          h_W_MCrad->Fill(W);
          /*if (processid==1){h_W_qe_MC->Fill(W);}
          else if (processid==2){h_W_mec_MC->Fill(W);}
          else if (processid==3){h_W_res_MC->Fill(W);}
          else if (processid==4){h_W_dis_MC->Fill(W);}*/
          h2_q2_W_MCrad->Fill(W,q2);
          h2_theta_W_MCrad->Fill(W,v4_el.Theta()*TMath::RadToDeg());
          //if(iCharge==0)      h2_Ppip_theta_MC->Fill(v4_el.Theta()*TMath::RadToDeg(),P_pi);
          //else if(iCharge==1) h2_Ppim_theta_MC->Fill(v4_el.Theta()*TMath::RadToDeg(),P_pi);
          //h_weight_MCrad->Fill(processid); h_processid_MCrad->Fill(processid_true);


          //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ///MC CS calculations

            if(bQorNot == true){/// theta_piq
                EventBinner(bQorNot, iCharge, P_pi, theta_piq, q2, W, CS_Wmin, dW, Wbini, iPmom, iThPiQ, iQ2, iW);
            }
            else{/// changed theta_piq to Theta_pi
                EventBinner(bQorNot, iCharge, P_pi, Theta_pi, q2, W, CS_Wmin, dW, Wbini, iPmom, iThPiQ, iQ2, iW);
            }
            if(iPmom == -1 || iThPiQ == -1 || iQ2 == -1 || iW == -1 || MCrad_weight<0. || MCrad_weight>100.) continue;
            //else if(iCharge == 0){cs_MC++; continue;}

            //Fill events
            sec_i = electrons[0]->trk(DC)->getSector() - 1;// getSector returns actual sector number
            if(sec_i > 3) sec_i--;
            if(sec_i > -1 && sec_i < sectors)
            {
                y_i[iPmom][iThPiQ][iCharge][sec_i][iGENIEr][iQ2][iW] += 1.;
                tCS[iPmom][iThPiQ][iCharge][sec_i][iGENIEr][iQ2][iW] += MCrad_weight;
                counter_MCrad++;
            }
	        //cout << y_i[iPmom][iThPiQ][iCharge][sec_i][iGENIEr][iQ2][iW] << "\t" << counter_MCrad<<endl;
	      }//if(electrons.size()==1 && (electrons[0]->par()->getVz() < -0.00001 || electrons[0]->par()->getVz() > 0.00001) )
	      if(bAllEvents==false && counter_MCrad>MC_rad_evts[iGENmod]) break;
    }///GENIErad while loop ************************************************************************************************************************************

 if(bOPG==true)
  {
  //////////////////////////////////////////////////////////////////////////////
  //         Monte Carlo (onepigen)
  //////////////////////////////////////////////////////////////////////////////
  while(chain_opg.Next())
    {
      v4_beam_opg.SetE(c12_opg->mcevent()->getEbeam());
      v4_beam_opg.SetPz(c12_opg->mcevent()->getEbeam());

      if( c12_opg->mcevent()->getEbeam() > 1e-4)
	      {
	        //if( counter_opg == 0)
	        //cout<<"Beam energy set manually: "<<beamE<<endl;
	        v4_beam_opg.SetE(c12_opg->mcevent()->getEbeam());
	        v4_beam_opg.SetPz(c12_opg->mcevent()->getEbeam());
	      }

      //intCS
      i_pipr++;
      if(i_pipr > piprN[i_piprF])
      {
          cout << "OLD FILE (opg_pr): "<<piprN[i_piprF]<<"\t"<<i_pisum<<endl;
          i_piprF++; i_pisum += i_pipr-1;
          i_pipr = 1;
          cout << "NEW FILE (opg_pr): "<<piprN[i_piprF]<<"\t"<<i_pisum<<endl;
      }

      // get particles by type
      auto electrons=c12_opg->getByID(11);
      auto protons=c12_opg->getByID(2212);
      auto neutrons=c12_opg->getByID(2112);
      auto pips=c12_opg->getByID(211);
      auto pims=c12_opg->getByID(-211);

      //Make sure that we have a good electron and that we are not seeing the weird delta function around zero in MC && (electrons[0]->par()->getVz() < -0.00001 || electrons[0]->par()->getVz() > 0.00001)
      if(electrons.size()==1 && electrons[0]->trk(DC)->getSector() != 4 && (pips.size()>=1 || pims.size()>=1))// && (pips.size()>=1 || pims.size()>=1)
	      {
          SetLorentzVector(v4_el,electrons[0]);

          TLorentzVector v4_q = v4_beam_opg - v4_el; //photon 4-vector

          //Calculate some important kinematic variables
          Double_t q2        = -v4_q.M2();
          Double_t x_b       = q2/(2 * mass_p * v4_q.E() ); //Bjorken-x
          Double_t y         = -v4_q.P() + sqrt( v4_q.E()*v4_q.E() + 2*v4_q.E()*mass_p); //y (Bjorken-y), the scaling variable; from Jourdan, eq. 20
          Double_t vz_e      = electrons[0]->par()->getVz();
          Double_t W         = sqrt(mass_p*mass_p - q2 + 2*v4_q.E()*mass_p); //Invariant mass, assuming electron off of a standing proton
          weight_opg         = c12_opg->mcevent()->getWeight();

          if(W < Wcut) continue; /// W min cut (removes QE "pion" [perhaps coincident proton] events with high ThetaQ)
          else if(W > CS_Wmax) continue; /// Save time

          ///***********************************************************************************************************************************************************************************************
            Int_t iPmom = -1, iThPiQ = -1, iCharge = -1, iMod = -1, iQ2 = -1, iW = -1;
            ///Electron vertex calculation
            vec_q.SetMagThetaPhi(v4_q.P(),v4_q.Theta(),v4_q.Phi());
            Vperp_e = sqrt(pow(electrons[0]->par()->getVx(),2) + pow(electrons[0]->par()->getVy(),2));
            P_pi = 0.0; i_pi = 0; Double_t P_can = 0.0;// pion momentum, pion # for the event, momentum of pion candidate
            Double_t Lv = electrons[0]->cal(PCAL)->getLv(), Lw = electrons[0]->cal(PCAL)->getLw();
            Double_t beta = 0.;

            ///Electron EC fiducial and vertex cuts
            if(eCutter(VWcut, Lv, Lw, vz_e, vz_e_min, vz_e_max, Vperp_e, Vperp_ecut)!=0) continue;
            //DC fiducial cuts for electron
            bool el_fid_lay_1 = dcFid.DC_pi_fid(electrons[0]->traj(DC,6)->getX(), electrons[0]->traj(DC,6)->getY(),  1, electrons[0]->trk(DC)->getSector(), 1,bending);
            bool el_fid_lay_2 = dcFid.DC_pi_fid(electrons[0]->traj(DC,18)->getX(),electrons[0]->traj(DC,18)->getY(), 1, electrons[0]->trk(DC)->getSector(), 2,bending);
            bool el_fid_lay_3 = dcFid.DC_pi_fid(electrons[0]->traj(DC,36)->getX(),electrons[0]->traj(DC,36)->getY(), 1, electrons[0]->trk(DC)->getSector(), 3,bending);
            if( !(el_fid_lay_1 && el_fid_lay_2 && el_fid_lay_3) ) continue;

            ///Highest momentum Pion Calculator
            if(pips.size()>=1)
            {
                for(Int_t i=0;i<pips.size();i++)
                {
                  Double_t beta_pi  = pips[i]->par()->getBeta();
                  P_can = pips[i]->par()->getP();
                  ///Remove CD, FTOF2, and high chi2pid pions (getSector returns actual sector number)
                  if(pips[i]->par()->getStatus() < 2000 || pips[i]->par()->getStatus() > 3999){continue;}
                  else if(pips[i]->sci(FTOF2)->getEnergy()>0.){continue;}
                  else if(abs(pips[i]->getChi2Pid()-chi2pid_pip) > 3.){continue;}
                  Vperp_pi = sqrt(pow(pips[i]->par()->getVx(),2) + pow(pips[i]->par()->getVy(),2));
                  Vdiff_pi = vz_e - pips[i]->par()->getVz();
                  ///Pion perpendicular and electron vertex difference cuts
                  if(PiCutter(0, P_can, 0.5, 1.0, Vperp_pi, Vdiff_pi, Vperp_pipcut, Vperp_pip_ms, Vperp_pimcut, Vperp_pim_ms, PiCutSig)!=0) continue;
                  ///Pion DC fiducial cuts
                  bool pip_fid_lay_1 = dcFid.DC_pi_fid(pips[i]->traj(DC,6)->getX(), pips[i]->traj(DC,6)->getY(), 3, pips[i]->trk(DC)->getSector(), 1,bending);
                  bool pip_fid_lay_2 = dcFid.DC_pi_fid(pips[i]->traj(DC,18)->getX(),pips[i]->traj(DC,18)->getY(),3, pips[i]->trk(DC)->getSector(), 2,bending);
                  bool pip_fid_lay_3 = dcFid.DC_pi_fid(pips[i]->traj(DC,36)->getX(),pips[i]->traj(DC,36)->getY(),3, pips[i]->trk(DC)->getSector(), 3,bending);
                  if( !(pip_fid_lay_1 && pip_fid_lay_2 && pip_fid_lay_3) ) continue;
                  ///Get maximum momentum pion
                  if(P_pi < P_can && pips[i]->par()->getCharge() == 1){P_pi = P_can; i_pi = i; iCharge = 0; beta = beta_pi;}
                }
            }
            if(pims.size()>=1)
            {
                for(Int_t i=0;i<pims.size();i++)
                {
                  Double_t beta_pi  = pims[i]->par()->getBeta();
                  ///Remove CD, FTOF2, and high chi2pid pions (getSector returns actual sector number)
                  if(pims[i]->par()->getStatus() < 2000 || pims[i]->par()->getStatus() > 3999){continue;}
                  else if(pims[i]->sci(FTOF2)->getEnergy()>0.){continue;}
                  else if(abs(pims[i]->getChi2Pid()-chi2pid_pim) > 3.){continue;}
                  Vperp_pi = sqrt(pow(pims[i]->par()->getVx(),2) + pow(pims[i]->par()->getVy(),2));
                  Vdiff_pi = vz_e - pims[i]->par()->getVz();
                  P_can = pims[i]->par()->getP();
                  ///Pion perpendicular and electron vertex difference cuts
                  if(PiCutter(1, P_can, 0.5, 1.0, Vperp_pi, Vdiff_pi, Vperp_pipcut, Vperp_pip_ms, Vperp_pimcut, Vperp_pim_ms, PiCutSig)!=0) continue;
                  ///Pion DC fiducial cuts
                  bool pim_fid_lay_1 = dcFid.DC_pi_fid(pims[i]->traj(DC,6)->getX(), pims[i]->traj(DC,6)->getY(), 4, pims[i]->trk(DC)->getSector(), 1,bending);
                  bool pim_fid_lay_2 = dcFid.DC_pi_fid(pims[i]->traj(DC,18)->getX(),pims[i]->traj(DC,18)->getY(),4, pims[i]->trk(DC)->getSector(), 2,bending);
                  bool pim_fid_lay_3 = dcFid.DC_pi_fid(pims[i]->traj(DC,36)->getX(),pims[i]->traj(DC,36)->getY(),4, pims[i]->trk(DC)->getSector(), 3,bending);
                  if( !(pim_fid_lay_1 && pim_fid_lay_2 && pim_fid_lay_3) ) continue;
                  ///Get maximum momentum pion
                  if(P_pi < P_can && pims[i]->par()->getCharge() == -1){P_pi = P_can; i_pi = i; iCharge = 1; beta = beta_pi;}
                }
            }
            ///Remove events with no acceptable pions
            if(P_pi==0.0) continue;
            ///Pion vertex calculation
            if(iCharge==0)//pi+
            {
                Vperp_pi = sqrt(pow(pips[i_pi]->par()->getVx(),2) + pow(pips[i_pi]->par()->getVy(),2));
                Vdiff_pi = vz_e - pips[i_pi]->par()->getVz();
                Theta_pi = pips[i_pi]->getTheta();
                Phi_pi   = pips[i_pi]->getPhi(); //i_pip++;
            }
            else if(iCharge==1)//pi-
            {
                Vperp_pi = sqrt(pow(pims[i_pi]->par()->getVx(),2) + pow(pims[i_pi]->par()->getVy(),2));
                Vdiff_pi = vz_e - pims[i_pi]->par()->getVz();
                Theta_pi = pims[i_pi]->getTheta();
                Phi_pi   = pims[i_pi]->getPhi(); //i_pim++;
            }
            vec_pi.SetMagThetaPhi(P_pi,Theta_pi,Phi_pi);
            theta_piq = vec_q.Angle(vec_pi)*TMath::RadToDeg();
            Theta_pi *= TMath::RadToDeg();

            //h2_theta_mom_opg->Fill(el.P(),el.Theta()*TMath::RadToDeg());
            h_q2_opg->Fill(q2);
            if(P_pi>0.) h_Ppi_opg->Fill(P_pi);
            if(theta_piq>0.) h_theta_piq_opg->Fill(theta_piq);
            if(Theta_pi>0.) h_theta_pi_opg->Fill(Theta_pi);/*
            h2_theta_q2_opg->Fill(q2,el.Theta()*TMath::RadToDeg());
            //h_W_opg->Fill(W);
            h2_q2_W_opg->Fill(W,q2);
            //h_omega_opg->Fill(v4_q.E());
            h2_theta_omega_opg->Fill(v4_q.E(),el.Theta()*TMath::RadToDeg());
            h2_theta_W_opg->Fill(W,el.Theta()*TMath::RadToDeg());*/

          //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ///opg CS calculations

            if(bQorNot == true){/// theta_piq
                EventBinner(bQorNot, iCharge, P_pi, theta_piq, q2, W, CS_Wmin, dW, Wbini, iPmom, iThPiQ, iQ2, iW);
            }
            else{/// changed theta_piq to Theta_pi
                EventBinner(bQorNot, iCharge, P_pi, Theta_pi, q2, W, CS_Wmin, dW, Wbini, iPmom, iThPiQ, iQ2, iW);
            }
            if(iPmom == -1 || iThPiQ == -1 || iQ2 == -1 || iW == -1) continue;
            else if(iCharge == 1){cs_opg++; continue;}

            //Fill events
            sec_i = electrons[0]->trk(DC)->getSector() - 1;// getSector returns actual sector number
            if(sec_i > 3) sec_i--;
            if(sec_i > -1 && sec_i < sectors)
            {
                y_i[iPmom][iThPiQ][iCharge][sec_i][iOPGr][iQ2][iW] += 1.;
                tCS[iPmom][iThPiQ][iCharge][sec_i][iOPGr][iQ2][iW] += intCS_pipr[i_piprF];
                counter_opg++; //cout << weight_opg << endl;
            }
	      }//if(electrons.size()==1 && (electrons[0]->par()->getVz() < -0.00001 || electrons[0]->par()->getVz() > 0.00001) )
    }///onepigen while loop ************************************************************************************************************************************
  //////////////////////////////////////////////////////////////////////////////
  //         Monte Carlo (onepigen)
  //////////////////////////////////////////////////////////////////////////////
  while(chain_opgm.Next())
    {
      v4_beam_opgm.SetE(c12_opgm->mcevent()->getEbeam());
      v4_beam_opgm.SetPz(c12_opgm->mcevent()->getEbeam());

      if( c12_opgm->mcevent()->getEbeam() > 1e-4)
	      {
	        //if( counter_opgm == 0)
	        //cout<<"Beam energy set manually: "<<beamE<<endl;
	        v4_beam_opgm.SetE(c12_opgm->mcevent()->getEbeam());
	        v4_beam_opgm.SetPz(c12_opgm->mcevent()->getEbeam());
	      }

      //intCS
      i_pimr++;
      //if(i_pimr < 5) cout << i_pimrF << " " << i_pimr << endl;
      //else if(i_pimr > 32000) cout << i_pimrF << " " << i_pimr << endl;
      if(i_pimr > pimrN[i_pimrF])
      {
          cout << "OLD FILE (opg_mr): "<<pimrN[i_pimrF]<<"\t"<<i_pisut<<endl;
          i_pimrF++; i_pisut += i_pimr-1;
          i_pimr = 1;
          cout << "NEW FILE (opg_mr): "<<pimrN[i_pimrF]<<"\t"<<i_pisut<<endl;
      }

      // get particles by type
      auto electrons=c12_opgm->getByID(11);
      auto protons=c12_opgm->getByID(2212);
      auto neutrons=c12_opgm->getByID(2112);
      auto pips=c12_opgm->getByID(211);
      auto pims=c12_opgm->getByID(-211);

      //Make sure that we have a good electron and that we are not seeing the weird delta function around zero in MC && (electrons[0]->par()->getVz() < -0.00001 || electrons[0]->par()->getVz() > 0.00001)
      if(electrons.size()==1 && electrons[0]->trk(DC)->getSector() != 4 && (pips.size()>=1 || pims.size()>=1))// && (pips.size()>=1 || pims.size()>=1)
	      {
          SetLorentzVector(v4_el,electrons[0]);

          TLorentzVector v4_q = v4_beam_opgm - v4_el; //photon 4-vector

          //Calculate some important kinematic variables
          Double_t q2        = -v4_q.M2();
          Double_t x_b       = q2/(2 * mass_n * v4_q.E() ); //Bjorken-x
          Double_t y         = -v4_q.P() + sqrt( v4_q.E()*v4_q.E() + 2*v4_q.E()*mass_n); //y (Bjorken-y), the scaling variable; from Jourdan, eq. 20
          Double_t vz_e      = electrons[0]->par()->getVz();
          Double_t W         = sqrt(mass_n*mass_n - q2 + 2*v4_q.E()*mass_n); //Invariant mass, assuming electron off of a standing nucleon
          weight_opgm        = c12_opgm->mcevent()->getWeight();

          if(W < Wcut) continue; /// W min cut (removes QE "pion" [perhaps coincident proton] events with high ThetaQ)
          else if(W > CS_Wmax) continue; /// Save time

          ///***********************************************************************************************************************************************************************************************
            Int_t iPmom = -1, iThPiQ = -1, iCharge = -1, iMod = -1, iQ2 = -1, iW = -1;
            ///Electron vertex calculation
            vec_q.SetMagThetaPhi(v4_q.P(),v4_q.Theta(),v4_q.Phi());
            Vperp_e = sqrt(pow(electrons[0]->par()->getVx(),2) + pow(electrons[0]->par()->getVy(),2));
            P_pi = 0.0; i_pi = 0; Double_t P_can = 0.0;// pion momentum, pion # for the event, momentum of pion candidate
            Double_t Lv = electrons[0]->cal(PCAL)->getLv(), Lw = electrons[0]->cal(PCAL)->getLw();
            Double_t beta = 0.;

            ///Electron EC fiducial and vertex cuts
            if(eCutter(VWcut, Lv, Lw, vz_e, vz_e_min, vz_e_max, Vperp_e, Vperp_ecut)!=0) continue;
            //DC fiducial cuts for electron
            bool el_fid_lay_1 = dcFid.DC_pi_fid(electrons[0]->traj(DC,6)->getX(), electrons[0]->traj(DC,6)->getY(),  1, electrons[0]->trk(DC)->getSector(), 1,bending);
            bool el_fid_lay_2 = dcFid.DC_pi_fid(electrons[0]->traj(DC,18)->getX(),electrons[0]->traj(DC,18)->getY(), 1, electrons[0]->trk(DC)->getSector(), 2,bending);
            bool el_fid_lay_3 = dcFid.DC_pi_fid(electrons[0]->traj(DC,36)->getX(),electrons[0]->traj(DC,36)->getY(), 1, electrons[0]->trk(DC)->getSector(), 3,bending);
            if( !(el_fid_lay_1 && el_fid_lay_2 && el_fid_lay_3) ) continue;

            ///Highest momentum Pion Calculator
            if(pips.size()>=1)
            {
                for(Int_t i=0;i<pips.size();i++)
                {
                  Double_t beta_pi  = pips[i]->par()->getBeta();
                  P_can = pips[i]->par()->getP();
                  ///Remove CD, FTOF2, and high chi2pid pions (getSector returns actual sector number)
                  if(pips[i]->par()->getStatus() < 2000 || pips[i]->par()->getStatus() > 3999){continue;}
                  else if(pips[i]->sci(FTOF2)->getEnergy()>0.){continue;}
                  else if(abs(pips[i]->getChi2Pid()-chi2pid_pip) > 3.){continue;}
                  Vperp_pi = sqrt(pow(pips[i]->par()->getVx(),2) + pow(pips[i]->par()->getVy(),2));
                  Vdiff_pi = vz_e - pips[i]->par()->getVz();
                  ///Pion perpendicular and electron vertex difference cuts
                  if(PiCutter(0, P_can, 0.5, 1.0, Vperp_pi, Vdiff_pi, Vperp_pipcut, Vperp_pip_ms, Vperp_pimcut, Vperp_pim_ms, PiCutSig)!=0) continue;
                  ///Pion DC fiducial cuts
                  bool pip_fid_lay_1 = dcFid.DC_pi_fid(pips[i]->traj(DC,6)->getX(), pips[i]->traj(DC,6)->getY(), 3, pips[i]->trk(DC)->getSector(), 1,bending);
                  bool pip_fid_lay_2 = dcFid.DC_pi_fid(pips[i]->traj(DC,18)->getX(),pips[i]->traj(DC,18)->getY(),3, pips[i]->trk(DC)->getSector(), 2,bending);
                  bool pip_fid_lay_3 = dcFid.DC_pi_fid(pips[i]->traj(DC,36)->getX(),pips[i]->traj(DC,36)->getY(),3, pips[i]->trk(DC)->getSector(), 3,bending);
                  if( !(pip_fid_lay_1 && pip_fid_lay_2 && pip_fid_lay_3) ) continue;
                  ///Get maximum momentum pion
                  if(P_pi < P_can && pips[i]->par()->getCharge() == 1){P_pi = P_can; i_pi = i; iCharge = 0; beta = beta_pi;}
                }
            }
            if(pims.size()>=1)
            {
                for(Int_t i=0;i<pims.size();i++)
                {
                  Double_t beta_pi  = pims[i]->par()->getBeta();
                  ///Remove CD, FTOF2, and high chi2pid pions (getSector returns actual sector number)
                  if(pims[i]->par()->getStatus() < 2000 || pims[i]->par()->getStatus() > 3999){continue;}
                  else if(pims[i]->sci(FTOF2)->getEnergy()>0.){continue;}
                  else if(abs(pims[i]->getChi2Pid()-chi2pid_pim) > 3.){continue;}
                  Vperp_pi = sqrt(pow(pims[i]->par()->getVx(),2) + pow(pims[i]->par()->getVy(),2));
                  Vdiff_pi = vz_e - pims[i]->par()->getVz();
                  P_can = pims[i]->par()->getP();
                  ///Pion perpendicular and electron vertex difference cuts
                  if(PiCutter(1, P_can, 0.5, 1.0, Vperp_pi, Vdiff_pi, Vperp_pipcut, Vperp_pip_ms, Vperp_pimcut, Vperp_pim_ms, PiCutSig)!=0) continue;
                  ///Pion DC fiducial cuts
                  bool pim_fid_lay_1 = dcFid.DC_pi_fid(pims[i]->traj(DC,6)->getX(), pims[i]->traj(DC,6)->getY(), 4, pims[i]->trk(DC)->getSector(), 1,bending);
                  bool pim_fid_lay_2 = dcFid.DC_pi_fid(pims[i]->traj(DC,18)->getX(),pims[i]->traj(DC,18)->getY(),4, pims[i]->trk(DC)->getSector(), 2,bending);
                  bool pim_fid_lay_3 = dcFid.DC_pi_fid(pims[i]->traj(DC,36)->getX(),pims[i]->traj(DC,36)->getY(),4, pims[i]->trk(DC)->getSector(), 3,bending);
                  if( !(pim_fid_lay_1 && pim_fid_lay_2 && pim_fid_lay_3) ) continue;
                  ///Get maximum momentum pion
                  if(P_pi < P_can && pims[i]->par()->getCharge() == -1){P_pi = P_can; i_pi = i; iCharge = 1; beta = beta_pi;}
                }
            }
            ///Remove events with no acceptable pions
            if(P_pi==0.0) continue;
            ///Pion vertex calculation
            if(iCharge==0)
            {
                Vperp_pi = sqrt(pow(pips[i_pi]->par()->getVx(),2) + pow(pips[i_pi]->par()->getVy(),2));
                Vdiff_pi = vz_e - pips[i_pi]->par()->getVz();
                Theta_pi = pips[i_pi]->getTheta();
                Phi_pi   = pips[i_pi]->getPhi(); //i_pip++;
            }
            else if(iCharge==1)
            {
                Vperp_pi = sqrt(pow(pims[i_pi]->par()->getVx(),2) + pow(pims[i_pi]->par()->getVy(),2));
                Vdiff_pi = vz_e - pims[i_pi]->par()->getVz();
                Theta_pi = pims[i_pi]->getTheta();
                Phi_pi   = pims[i_pi]->getPhi(); //i_pim++;
            }
            vec_pi.SetMagThetaPhi(P_pi,Theta_pi,Phi_pi);
            theta_piq = vec_q.Angle(vec_pi)*TMath::RadToDeg();
            Theta_pi *= TMath::RadToDeg();

            h_q2_opg->Fill(q2);
            if(P_pi>0.) h_Ppi_opg->Fill(P_pi);
            if(theta_piq>0.) h_theta_piq_opg->Fill(theta_piq);
            if(Theta_pi>0.) h_theta_pi_opg->Fill(Theta_pi);
          //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ///opg CS calculations

            if(bQorNot == true){/// theta_piq
                EventBinner(bQorNot, iCharge, P_pi, theta_piq, q2, W, CS_Wmin, dW, Wbini, iPmom, iThPiQ, iQ2, iW);
            }
            else{/// changed theta_piq to Theta_pi
                EventBinner(bQorNot, iCharge, P_pi, Theta_pi, q2, W, CS_Wmin, dW, Wbini, iPmom, iThPiQ, iQ2, iW);
            }
            if(iPmom == -1 || iThPiQ == -1 || iQ2 == -1 || iW == -1) continue;
            else if(iCharge == 0){cs_opgm++; continue;}

            //Fill events
            sec_i = electrons[0]->trk(DC)->getSector() - 1;// getSector returns actual sector number
            if(sec_i > 3) sec_i--;
            if(sec_i > -1 && sec_i < sectors)
            {
                y_i[iPmom][iThPiQ][iCharge][sec_i][iOPGr][iQ2][iW] += 1.;
                tCS[iPmom][iThPiQ][iCharge][sec_i][iOPGr][iQ2][iW] += intCS_pimr[i_pimrF];
                counter_opgm++; //cout << weight_opg << endl;
            }
	      }//if(electrons.size()==1 && (electrons[0]->par()->getVz() < -0.00001 || electrons[0]->par()->getVz() > 0.00001) )
    }///onepigen while loop ************************************************************************************************************************************

  //////////////////////////////////////////////////////////////////////////////
  //         Monte Carlo (onepigen)
  //////////////////////////////////////////////////////////////////////////////
  while(chain_opgn.Next())
    {
      v4_beam_opgn.SetE(c12_opgn->mcevent()->getEbeam());
      v4_beam_opgn.SetPz(c12_opgn->mcevent()->getEbeam());

      if( c12_opgn->mcevent()->getEbeam() > 1e-4)
	      {
	        //if( counter_opgn == 0)
	        //cout<<"Beam energy set manually: "<<beamE<<endl;
	        v4_beam_opgn.SetE(c12_opgn->mcevent()->getEbeam());
	        v4_beam_opgn.SetPz(c12_opgn->mcevent()->getEbeam());
	      }

      //intCS
      i_pipn++;
      if(i_pipn > pipnN[i_pipnF])
      {
          cout << "OLD FILE (opg_pn): "<<pipnN[i_pipnF]<<"\t"<<i_pisup<<endl;
          i_pipnF++; i_pisup += i_pipn-1;
          i_pipn = 1;
          cout << "NEW FILE (opg_pn): "<<pipnN[i_pipnF]<<"\t"<<i_pisup<<endl;
      }

      // get particles by type
      auto electrons=c12_opgn->getByID(11);
      auto protons=c12_opgn->getByID(2212);
      auto neutrons=c12_opgn->getByID(2112);
      auto pips=c12_opgn->getByID(211);
      auto pims=c12_opgn->getByID(-211);

      //Make sure that we have a good electron and that we are not seeing the weird delta function around zero in MC && (electrons[0]->par()->getVz() < -0.00001 || electrons[0]->par()->getVz() > 0.00001)
      if(electrons.size()==1 && electrons[0]->trk(DC)->getSector() != 4 && (pips.size()>=1 || pims.size()>=1))// && (pips.size()>=1 || pims.size()>=1)
	      {
          SetLorentzVector(v4_el,electrons[0]);

          TLorentzVector v4_q = v4_beam_opgn - v4_el; //photon 4-vector

          //Calculate some important kinematic variables
          Double_t q2        = -v4_q.M2();
          Double_t x_b       = q2/(2 * mass_p * v4_q.E() ); //Bjorken-x
          Double_t y         = -v4_q.P() + sqrt( v4_q.E()*v4_q.E() + 2*v4_q.E()*mass_p); //y (Bjorken-y), the scaling variable; from Jourdan, eq. 20
          Double_t vz_e      = electrons[0]->par()->getVz();
          Double_t W         = sqrt(mass_p*mass_p - q2 + 2*v4_q.E()*mass_p); //Invariant mass, assuming electron off of a standing proton
          //weight_opgn        = c12_opgn->mcevent()->getWeight();

          if(W < Wcut) continue; /// W min cut (removes QE "pion" [perhaps coincident proton] events with high ThetaQ)
          else if(W > CS_Wmax) continue; /// Save time

          ///-**********************************************************************************************************************************************************************************************
            Int_t iPmom = -1, iThPiQ = -1, iCharge = -1, iMod = -1, iQ2 = -1, iW = -1;
            ///Electron vertex calculation
            vec_q.SetMagThetaPhi(v4_q.P(),v4_q.Theta(),v4_q.Phi());
            Vperp_e = sqrt(pow(electrons[0]->par()->getVx(),2) + pow(electrons[0]->par()->getVy(),2));
            P_pi = 0.0; i_pi = 0; Double_t P_can = 0.0;
            Double_t Lv = electrons[0]->cal(PCAL)->getLv(), Lw = electrons[0]->cal(PCAL)->getLw();
            Double_t beta = 0.;

            ///Electron EC fiducial and vertex cuts
            if(eCutter(VWcut, Lv, Lw, vz_e, vz_e_min, vz_e_max, Vperp_e, Vperp_ecut)!=0) continue;
            //DC fiducial cuts for electron
            bool el_fid_lay_1 = dcFid.DC_pi_fid(electrons[0]->traj(DC,6)->getX(), electrons[0]->traj(DC,6)->getY(),  1, electrons[0]->trk(DC)->getSector(), 1,bending);
            bool el_fid_lay_2 = dcFid.DC_pi_fid(electrons[0]->traj(DC,18)->getX(),electrons[0]->traj(DC,18)->getY(), 1, electrons[0]->trk(DC)->getSector(), 2,bending);
            bool el_fid_lay_3 = dcFid.DC_pi_fid(electrons[0]->traj(DC,36)->getX(),electrons[0]->traj(DC,36)->getY(), 1, electrons[0]->trk(DC)->getSector(), 3,bending);
            if( !(el_fid_lay_1 && el_fid_lay_2 && el_fid_lay_3) ) continue;

            ///Highest momentum Pion Calculator
            if(pips.size()>=1)
            {
                for(Int_t i=0;i<pips.size();i++)
                {
                  Double_t beta_pi  = pips[i]->par()->getBeta();
                  P_can = pips[i]->par()->getP();
                  ///Remove CD, FTOF2, and high chi2pid pions (getSector returns actual sector number)
                  if(pips[i]->par()->getStatus() < 2000 || pips[i]->par()->getStatus() > 3999){continue;}
                  else if(pips[i]->sci(FTOF2)->getEnergy()>0.){continue;}
                  else if(abs(pips[i]->getChi2Pid()-chi2pid_pip) > 3.){continue;}
                  Vperp_pi = sqrt(pow(pips[i]->par()->getVx(),2) + pow(pips[i]->par()->getVy(),2));
                  Vdiff_pi = vz_e - pips[i]->par()->getVz();
                  ///Pion perpendicular and electron vertex difference cuts
                  if(PiCutter(0, P_can, 0.5, 1.0, Vperp_pi, Vdiff_pi, Vperp_pipcut, Vperp_pip_ms, Vperp_pimcut, Vperp_pim_ms, PiCutSig)!=0) continue;
                  ///Pion DC fiducial cuts
                  bool pip_fid_lay_1 = dcFid.DC_pi_fid(pips[i]->traj(DC,6)->getX(), pips[i]->traj(DC,6)->getY(), 3, pips[i]->trk(DC)->getSector(), 1,bending);
                  bool pip_fid_lay_2 = dcFid.DC_pi_fid(pips[i]->traj(DC,18)->getX(),pips[i]->traj(DC,18)->getY(),3, pips[i]->trk(DC)->getSector(), 2,bending);
                  bool pip_fid_lay_3 = dcFid.DC_pi_fid(pips[i]->traj(DC,36)->getX(),pips[i]->traj(DC,36)->getY(),3, pips[i]->trk(DC)->getSector(), 3,bending);
                  if( !(pip_fid_lay_1 && pip_fid_lay_2 && pip_fid_lay_3) ) continue;
                  ///Get maximum momentum pion
                  if(P_pi < P_can && pips[i]->par()->getCharge() == 1){P_pi = P_can; i_pi = i; iCharge = 0; beta = beta_pi;}
                }
            }
            if(pims.size()>=1)
            {
                for(Int_t i=0;i<pims.size();i++)
                {
                  Double_t beta_pi  = pims[i]->par()->getBeta();
                  ///Remove CD, FTOF2, and high chi2pid pions (getSector returns actual sector number)
                  if(pims[i]->par()->getStatus() < 2000 || pims[i]->par()->getStatus() > 3999){continue;}
                  else if(pims[i]->sci(FTOF2)->getEnergy()>0.){continue;}
                  else if(abs(pims[i]->getChi2Pid()-chi2pid_pim) > 3.){continue;}
                  Vperp_pi = sqrt(pow(pims[i]->par()->getVx(),2) + pow(pims[i]->par()->getVy(),2));
                  Vdiff_pi = vz_e - pims[i]->par()->getVz();
                  P_can = pims[i]->par()->getP();
                  ///Pion perpendicular and electron vertex difference cuts
                  if(PiCutter(1, P_can, 0.5, 1.0, Vperp_pi, Vdiff_pi, Vperp_pipcut, Vperp_pip_ms, Vperp_pimcut, Vperp_pim_ms, PiCutSig)!=0) continue;
                  ///Pion DC fiducial cuts
                  bool pim_fid_lay_1 = dcFid.DC_pi_fid(pims[i]->traj(DC,6)->getX(), pims[i]->traj(DC,6)->getY(), 4, pims[i]->trk(DC)->getSector(), 1,bending);
                  bool pim_fid_lay_2 = dcFid.DC_pi_fid(pims[i]->traj(DC,18)->getX(),pims[i]->traj(DC,18)->getY(),4, pims[i]->trk(DC)->getSector(), 2,bending);
                  bool pim_fid_lay_3 = dcFid.DC_pi_fid(pims[i]->traj(DC,36)->getX(),pims[i]->traj(DC,36)->getY(),4, pims[i]->trk(DC)->getSector(), 3,bending);
                  if( !(pim_fid_lay_1 && pim_fid_lay_2 && pim_fid_lay_3) ) continue;
                  ///Get maximum momentum pion
                  if(P_pi < P_can && pims[i]->par()->getCharge() == -1){P_pi = P_can; i_pi = i; iCharge = 1; beta = beta_pi;}
                }
            }
            ///Remove events with no acceptable pions
            if(P_pi==0.0) continue;
            ///Pion vertex calculation
            if(iCharge==0)
            {
                Vperp_pi = sqrt(pow(pips[i_pi]->par()->getVx(),2) + pow(pips[i_pi]->par()->getVy(),2));
                Vdiff_pi = vz_e - pips[i_pi]->par()->getVz();
                Theta_pi = pips[i_pi]->getTheta();
                Phi_pi   = pips[i_pi]->getPhi(); //i_pip++;
            }
            else if(iCharge==1)
            {
                Vperp_pi = sqrt(pow(pims[i_pi]->par()->getVx(),2) + pow(pims[i_pi]->par()->getVy(),2));
                Vdiff_pi = vz_e - pims[i_pi]->par()->getVz();
                Theta_pi = pims[i_pi]->getTheta();
                Phi_pi   = pims[i_pi]->getPhi(); //i_pim++;
            }
            vec_pi.SetMagThetaPhi(P_pi,Theta_pi,Phi_pi);
            theta_piq = vec_q.Angle(vec_pi)*TMath::RadToDeg();
            Theta_pi *= TMath::RadToDeg();

          //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ///opg CS calculations

            if(bQorNot == true){/// theta_piq
                EventBinner(bQorNot, iCharge, P_pi, theta_piq, q2, W, CS_Wmin, dW, Wbini, iPmom, iThPiQ, iQ2, iW);
            }
            else{/// changed theta_piq to Theta_pi
                EventBinner(bQorNot, iCharge, P_pi, Theta_pi, q2, W, CS_Wmin, dW, Wbini, iPmom, iThPiQ, iQ2, iW);
            }
            if(iPmom == -1 || iThPiQ == -1 || iQ2 == -1 || iW == -1) continue;
            else if(iCharge == 1){cs_opgn++; continue;}

            /*if(i_pisup<51935)
            {
                h_tf->Fill(W);
            }
            else if(i_pisup<51935+52023)
            {
                h_tw->Fill(W);
            }
            else
            {
                h_twf->Fill(W);
            }*/
            Int_t iLowpg = 0;
            if(iQ2 == 0)
            {
                if(iThPiQ == 0)
                {
                    if(iPmom == 3 && W < 2.05)      iLowpg = 1;
                    else if(iPmom == 2 && W < 1.75) iLowpg = 1;
                    else if(iPmom == 1 && W < 1.45) iLowpg = 1;
                    else if(iPmom == 0 && W < 1.25) iLowpg = 1;
                }
                if(iThPiQ == 1)
                {
                    if(iPmom == 0 && W < 1.2)       iLowpg = 1;
                    else if(iPmom == 1 && W < 1.5)  iLowpg = 1;
                    else if(iPmom == 2 && W < 1.75) iLowpg = 1;
                }
                if(iThPiQ == 2)
                {
                    if(iPmom == 0 && W < 1.3)       iLowpg = 1;
                    else if(iPmom == 1 && W < 1.55) iLowpg = 1;
                }
                if(iThPiQ == 3)
                {
                    if(iPmom == 0 && W < 1.35)      iLowpg = 1;
                }
            }
            if(iLowpg == 1)
            {
                h_opg_tail_vzpip->Fill(vz_e);
                h_opg_tail_omega->Fill(v4_q.E());
            }
            else
            {
                h_opg_nota_vzpip->Fill(vz_e);
                h_opg_nota_omega->Fill(v4_q.E());
            }
            h_phi_e_opg->Fill(v4_q.Phi()*TMath::RadToDeg());
            h_phi_pi_opg->Fill(Phi_pi*TMath::RadToDeg());

            //Fill events
            sec_i = electrons[0]->trk(DC)->getSector() - 1;// getSector returns actual sector number
            if(sec_i > 3) sec_i--;
            if(sec_i > -1 && sec_i < sectors)
            {
                y_i[iPmom][iThPiQ][iCharge][sec_i][iOPGn][iQ2][iW] += 1.;
                tCS[iPmom][iThPiQ][iCharge][sec_i][iOPGn][iQ2][iW] += intCS_pipn[i_pipnF];
                counter_opgn++; //cout << weight_opg << endl;
	        }
	      }//if(electrons.size()==1 && (electrons[0]->par()->getVz() < -0.00001 || electrons[0]->par()->getVz() > 0.00001) )
    }///onepigen while loop ************************************************************************************************************************************
  h_tf->SetMarkerColor(kRed); h_tf->SetLineColor(kRed); hs_twf->Add(h_tf);
  h_tw->SetMarkerColor(kBlue); h_tw->SetLineColor(kBlue); hs_twf->Add(h_tw);
  h_twf->SetMarkerColor(kViolet); h_twf->SetLineColor(kViolet); hs_twf->Add(h_twf);
  while(chain_opgmn.Next())
    {
      v4_beam_opgmn.SetE(c12_opgmn->mcevent()->getEbeam());
      v4_beam_opgmn.SetPz(c12_opgmn->mcevent()->getEbeam());

      if( c12_opgmn->mcevent()->getEbeam() > 1e-4)
	      {
	        //if( counter_opgmn == 0)
	        //cout<<"Beam energy set manually: "<<beamE<<endl;
	        v4_beam_opgmn.SetE(c12_opgmn->mcevent()->getEbeam());
	        v4_beam_opgmn.SetPz(c12_opgmn->mcevent()->getEbeam());
	      }

      //intCS
      i_pimn++;
      if(i_pimn > pimnN[i_pimnF])
      {
          cout << "OLD FILE (opg_mn): "<<pimnN[i_pimnF]<<"\t"<<i_pisun<<endl;
          i_pimnF++; i_pisun += i_pimn-1;
          i_pimn = 1;
          cout << "NEW FILE (opg_mn): "<<pimnN[i_pimnF]<<"\t"<<i_pisun<<endl;
      }

      // get particles by type
      auto electrons=c12_opgmn->getByID(11);
      auto protons=c12_opgmn->getByID(2212);
      auto neutrons=c12_opgmn->getByID(2112);
      auto pips=c12_opgmn->getByID(211);
      auto pims=c12_opgmn->getByID(-211);

      //Make sure that we have a good electron and that we are not seeing the weird delta function around zero in MC && (electrons[0]->par()->getVz() < -0.00001 || electrons[0]->par()->getVz() > 0.00001)
      if(electrons.size()==1 && electrons[0]->trk(DC)->getSector() != 4 && (pips.size()>=1 || pims.size()>=1))// && (pips.size()>=1 || pims.size()>=1)
	      {
          SetLorentzVector(v4_el,electrons[0]);

          TLorentzVector v4_q = v4_beam_opgmn - v4_el; //photon 4-vector

          //Calculate some important kinematic variables
          Double_t q2        = -v4_q.M2();
          Double_t x_b       = q2/(2 * mass_n * v4_q.E() ); //Bjorken-x
          Double_t y         = -v4_q.P() + sqrt( v4_q.E()*v4_q.E() + 2*v4_q.E()*mass_n); //y (Bjorken-y), the scaling variable; from Jourdan, eq. 20
          Double_t vz_e      = electrons[0]->par()->getVz();
          Double_t W         = sqrt(mass_n*mass_n - q2 + 2*v4_q.E()*mass_n); //Invariant mass, assuming electron off of a standing proton
          //weight_opgmn       = c12_opgmn->mcevent()->getWeight();

          if(W < Wcut) continue; /// W min cut (removes QE "pion" [perhaps coincident proton] events with high ThetaQ)
          else if(W > CS_Wmax) continue; /// Save time

          ///-**********************************************************************************************************************************************************************************************
            Int_t iPmom = -1, iThPiQ = -1, iCharge = -1, iMod = -1, iQ2 = -1, iW = -1;
            ///Electron vertex calculation
            vec_q.SetMagThetaPhi(v4_q.P(),v4_q.Theta(),v4_q.Phi());
            Vperp_e = sqrt(pow(electrons[0]->par()->getVx(),2) + pow(electrons[0]->par()->getVy(),2));
            P_pi = 0.0; i_pi = 0; Double_t P_can = 0.0;
            Double_t Lv = electrons[0]->cal(PCAL)->getLv(), Lw = electrons[0]->cal(PCAL)->getLw();
            Double_t beta = 0.;

            ///Electron EC fiducial and vertex cuts
            if(eCutter(VWcut, Lv, Lw, vz_e, vz_e_min, vz_e_max, Vperp_e, Vperp_ecut)!=0) continue;
            //DC fiducial cuts for electron
            bool el_fid_lay_1 = dcFid.DC_pi_fid(electrons[0]->traj(DC,6)->getX(), electrons[0]->traj(DC,6)->getY(),  1, electrons[0]->trk(DC)->getSector(), 1,bending);
            bool el_fid_lay_2 = dcFid.DC_pi_fid(electrons[0]->traj(DC,18)->getX(),electrons[0]->traj(DC,18)->getY(), 1, electrons[0]->trk(DC)->getSector(), 2,bending);
            bool el_fid_lay_3 = dcFid.DC_pi_fid(electrons[0]->traj(DC,36)->getX(),electrons[0]->traj(DC,36)->getY(), 1, electrons[0]->trk(DC)->getSector(), 3,bending);
            if( !(el_fid_lay_1 && el_fid_lay_2 && el_fid_lay_3) ) continue;

            ///Highest momentum Pion Calculator
            if(pips.size()>=1)
            {
                for(Int_t i=0;i<pips.size();i++)
                {
                  Double_t beta_pi  = pips[i]->par()->getBeta();
                  P_can = pips[i]->par()->getP();
                  ///Remove CD, FTOF2, and high chi2pid pions (getSector returns actual sector number)
                  if(pips[i]->par()->getStatus() < 2000 || pips[i]->par()->getStatus() > 3999){continue;}
                  else if(pips[i]->sci(FTOF2)->getEnergy()>0.){continue;}
                  else if(abs(pips[i]->getChi2Pid()-chi2pid_pip) > 3.){continue;}
                  Vperp_pi = sqrt(pow(pips[i]->par()->getVx(),2) + pow(pips[i]->par()->getVy(),2));
                  Vdiff_pi = vz_e - pips[i]->par()->getVz();
                  ///Pion perpendicular and electron vertex difference cuts
                  if(PiCutter(0, P_can, 0.5, 1.0, Vperp_pi, Vdiff_pi, Vperp_pipcut, Vperp_pip_ms, Vperp_pimcut, Vperp_pim_ms, PiCutSig)!=0) continue;
                  ///Pion DC fiducial cuts
                  bool pip_fid_lay_1 = dcFid.DC_pi_fid(pips[i]->traj(DC,6)->getX(), pips[i]->traj(DC,6)->getY(), 3, pips[i]->trk(DC)->getSector(), 1,bending);
                  bool pip_fid_lay_2 = dcFid.DC_pi_fid(pips[i]->traj(DC,18)->getX(),pips[i]->traj(DC,18)->getY(),3, pips[i]->trk(DC)->getSector(), 2,bending);
                  bool pip_fid_lay_3 = dcFid.DC_pi_fid(pips[i]->traj(DC,36)->getX(),pips[i]->traj(DC,36)->getY(),3, pips[i]->trk(DC)->getSector(), 3,bending);
                  if( !(pip_fid_lay_1 && pip_fid_lay_2 && pip_fid_lay_3) ) continue;
                  ///Get maximum momentum pion
                  if(P_pi < P_can && pips[i]->par()->getCharge() == 1){P_pi = P_can; i_pi = i; iCharge = 0; beta = beta_pi;}
                }
            }
            if(pims.size()>=1)
            {
                for(Int_t i=0;i<pims.size();i++)
                {
                  Double_t beta_pi  = pims[i]->par()->getBeta();
                  ///Remove CD, FTOF2, and high chi2pid pions (getSector returns actual sector number)
                  if(pims[i]->par()->getStatus() < 2000 || pims[i]->par()->getStatus() > 3999){continue;}
                  else if(pims[i]->sci(FTOF2)->getEnergy()>0.){continue;}
                  else if(abs(pims[i]->getChi2Pid()-chi2pid_pim) > 3.){continue;}
                  Vperp_pi = sqrt(pow(pims[i]->par()->getVx(),2) + pow(pims[i]->par()->getVy(),2));
                  Vdiff_pi = vz_e - pims[i]->par()->getVz();
                  P_can = pims[i]->par()->getP();
                  ///Pion perpendicular and electron vertex difference cuts
                  if(PiCutter(1, P_can, 0.5, 1.0, Vperp_pi, Vdiff_pi, Vperp_pipcut, Vperp_pip_ms, Vperp_pimcut, Vperp_pim_ms, PiCutSig)!=0) continue;
                  ///Pion DC fiducial cuts
                  bool pim_fid_lay_1 = dcFid.DC_pi_fid(pims[i]->traj(DC,6)->getX(), pims[i]->traj(DC,6)->getY(), 4, pims[i]->trk(DC)->getSector(), 1,bending);
                  bool pim_fid_lay_2 = dcFid.DC_pi_fid(pims[i]->traj(DC,18)->getX(),pims[i]->traj(DC,18)->getY(),4, pims[i]->trk(DC)->getSector(), 2,bending);
                  bool pim_fid_lay_3 = dcFid.DC_pi_fid(pims[i]->traj(DC,36)->getX(),pims[i]->traj(DC,36)->getY(),4, pims[i]->trk(DC)->getSector(), 3,bending);
                  if( !(pim_fid_lay_1 && pim_fid_lay_2 && pim_fid_lay_3) ) continue;
                  ///Get maximum momentum pion
                  if(P_pi < P_can && pims[i]->par()->getCharge() == -1){P_pi = P_can; i_pi = i; iCharge = 1; beta = beta_pi;}
                }
            }
            ///Remove events with no acceptable pions
            if(P_pi==0.0) continue;
            ///Pion vertex calculation
            if(iCharge==0)
            {
                Vperp_pi = sqrt(pow(pips[i_pi]->par()->getVx(),2) + pow(pips[i_pi]->par()->getVy(),2));
                Vdiff_pi = vz_e - pips[i_pi]->par()->getVz();
                Theta_pi = pips[i_pi]->getTheta();
                Phi_pi   = pips[i_pi]->getPhi(); //i_pip++;
            }
            else if(iCharge==1)
            {
                Vperp_pi = sqrt(pow(pims[i_pi]->par()->getVx(),2) + pow(pims[i_pi]->par()->getVy(),2));
                Vdiff_pi = vz_e - pims[i_pi]->par()->getVz();
                Theta_pi = pims[i_pi]->getTheta();
                Phi_pi   = pims[i_pi]->getPhi(); //i_pim++;
            }
            vec_pi.SetMagThetaPhi(P_pi,Theta_pi,Phi_pi);
            theta_piq = vec_q.Angle(vec_pi)*TMath::RadToDeg();
            Theta_pi *= TMath::RadToDeg();

          //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ///opg CS calculations

            if(bQorNot == true){/// theta_piq
                EventBinner(bQorNot, iCharge, P_pi, theta_piq, q2, W, CS_Wmin, dW, Wbini, iPmom, iThPiQ, iQ2, iW);
            }
            else{/// changed theta_piq to Theta_pi
                EventBinner(bQorNot, iCharge, P_pi, Theta_pi, q2, W, CS_Wmin, dW, Wbini, iPmom, iThPiQ, iQ2, iW);
            }
            if(iPmom == -1 || iThPiQ == -1 || iQ2 == -1 || iW == -1) continue;
            else if(iCharge == 0){cs_opgmn++; continue;}

            //Fill events
            sec_i = electrons[0]->trk(DC)->getSector() - 1;// getSector returns actual sector number
            if(sec_i > 3) sec_i--;
            if(sec_i > -1 && sec_i < sectors)
            {
                y_i[iPmom][iThPiQ][iCharge][sec_i][iOPGn][iQ2][iW] += 1.;
                tCS[iPmom][iThPiQ][iCharge][sec_i][iOPGn][iQ2][iW] += intCS_pimn[i_pimnF];
                counter_opgmn++; //cout << weight_opg << endl;
            }
	      }//if(electrons.size()==1 && (electrons[0]->par()->getVz() < -0.00001 || electrons[0]->par()->getVz() > 0.00001) )
    }///onepigen while loop ************************************************************************************************************************************
  }
  Double_t OPG_to_Data_Scaling = 1.;
  if(bData==true && bOPG==true)  OPG_to_Data_Scaling = ( (Double_t)counter_data ) / ( (Double_t)counter_opgn + (Double_t)counter_opgmn );
  hs_q2->Add(h_q2_data);  hs_q2->Add(h_q2_MC);
  h_q2_opg->SetMarkerColor(kBlue); h_q2_opg->Sumw2(); h_q2_opg->Scale(OPG_to_Data_Scaling); h_q2_opg->SetLineColor(kBlue); hs_q2->Add(h_q2_opg);
  //P_pi
  h_Ppi_data->SetLineColor(kBlack); hs_Ppi->Add(h_Ppi_data);
  h_Ppi_MC->SetMarkerColor(kRed); h_Ppi_MC->Sumw2(); h_Ppi_MC->Scale(MC_to_Data_Scaling); h_Ppi_MC->SetLineColor(kRed); hs_Ppi->Add(h_Ppi_MC);
  h_Ppi_opg->SetMarkerColor(kBlue); h_Ppi_opg->Sumw2(); h_Ppi_opg->Scale(OPG_to_Data_Scaling); h_Ppi_opg->SetLineColor(kBlue); hs_Ppi->Add(h_Ppi_opg);
  //theta_piq
  h_theta_piq_data->SetLineColor(kBlack); hs_theta_piq->Add(h_theta_piq_data);
  h_theta_piq_MC->SetMarkerColor(kRed); h_theta_piq_MC->Sumw2(); h_theta_piq_MC->Scale(MC_to_Data_Scaling); h_theta_piq_MC->SetLineColor(kRed); hs_theta_piq->Add(h_theta_piq_MC);
  h_theta_piq_opg->SetMarkerColor(kBlue); h_theta_piq_opg->Sumw2(); h_theta_piq_opg->Scale(OPG_to_Data_Scaling); h_theta_piq_opg->SetLineColor(kBlue); hs_theta_piq->Add(h_theta_piq_opg);
  //theta_pi
  h_theta_pi_data->SetLineColor(kBlack); hs_theta_pi->Add(h_theta_pi_data);
  //h_theta_pip_data->SetLineColor(15); hs_theta_pi->Add(h_theta_pip_data);
  //h_theta_pim_data->SetLineColor(12); hs_theta_pi->Add(h_theta_pim_data);
  h_theta_pi_MC->SetMarkerColor(kRed); h_theta_pi_MC->Sumw2(); h_theta_pi_MC->Scale(MC_to_Data_Scaling); h_theta_pi_MC->SetLineColor(kRed); hs_theta_pi->Add(h_theta_pi_MC);
  h_theta_pi_opg->SetMarkerColor(kBlue); h_theta_pi_opg->Sumw2(); h_theta_pi_opg->Scale(OPG_to_Data_Scaling); h_theta_pi_opg->SetLineColor(kBlue); hs_theta_pi->Add(h_theta_pi_opg);

///--------------------------------------------------------------------OUTPUT--------------------------------------------------------------------///
for(Int_t iCharge = 0; iCharge < Nrec; iCharge++)/// Calculating y_sum
{
  for(Int_t iQ2 = 0; iQ2 < Q2bini; iQ2++)
  {
    for(Int_t iPmom = 0; iPmom < Pbini; iPmom++)
    {
      Int_t iPcanv = (Pbini - 1) - iPmom;//Needed to get plot positions correct for P_pi
      for(Int_t iThPiQ = 0; iThPiQ < TPIQi; iThPiQ++)
      {
        /*dT = DeltaThetaPiQBin(ThetaPiQ[iThPiQ],ThetaPiQ[iThPiQ+1]);
        Double_t dPmom = dP[iPcanv];
        if(bQorNot==false)
        {
          dT = DeltaThetaPiQBin(ThetaPi[iCharge][iThPiQ],ThetaPi[iCharge][iThPiQ+1]);
          dPmom = dPt[iCharge][iPmom+1] - dPt[iCharge][iPmom];
        }
        bin_vol = dPmom * dT * dQ2[iQ2] * dW;*/
        for(Int_t iSector = 0; iSector < sectors; iSector++)
        {
          for(Int_t iMod = 0; iMod < models; iMod++)
          {
            for(Int_t iW = 0; iW < Wbini; iW++)
            {
              y_sum[iPcanv][iThPiQ][iCharge][iMod][iQ2][iW] += y_i[iPcanv][iThPiQ][iCharge][iSector][iMod][iQ2][iW];
            }   // Note: for an array bin # of iW, the ROOT bin # is iW+1
            /*if(iMod==iDATA)//Data
            {
                h_CSS->Scale(cm2ub * ub2nb / (L * bin_vol));
            }
            else if(iMod==iGENIE)//GENIE reconstructed
            {
                h_CSS->Scale(totalCS_MC[iGENmod] * ub2nb / (MC_evts[iGENmod] * bin_vol));
            }
            else if(iMod==iGEN)//GENIE generated
            {
                h_CSS->Scale(totalCS_MC[iGENmod] * ub2nb * (5./6.) / (MC_evts[iGENmod] * bin_vol));
            }
            else if(iMod==iOPGr && iCharge==0)//onepigen
            {
                h_CSS->Scale(ub2nb / (opg_gen_events[0] * bin_vol));
            }
            else if(iMod==iOPGr && iCharge==1)//onepigen pi-
            {
                h_CSS->Scale(ub2nb / (opg_gen_events[1] * bin_vol));/// <-------------------<<
            }
            else if(iMod==iOPGn && iCharge==0)//onepigen pi+ norad
            {
                h_CSS->Scale(ub2nb / (opg_gen_events[2] * bin_vol));
            }
            else if(iMod==iOPGn && iCharge==1)//onepigen pi- norad
            {
                h_CSS->Scale(ub2nb / (opg_gen_events[3] * bin_vol));
            }*/
          }
        }
      }
    }
  }
}
  TFile *outFile = new TFile(outputFile + ".root","RECREATE");
  ofstream txtFile ("BinEntries.txt");
  //ofstream CStxtFile ("allCSvalues.txt");

  TCanvas *c1 = new TCanvas("c1","c1",6000,4000);
  TString fileName = outputFile + ".pdf";
  TString file_CSSq1pip  = "CSSq1pip.C"; // Remember to change nullptr
  TString file_CSSq2pip  = "CSSq2pip.C"; // to NULL before looking at
  TString file_CSSq3pip  = "CSSq3pip.C"; // the file.
  TString file_CSSq1pim  = "CSSq1pim.C";
  TString file_CSSq2pim  = "CSSq2pim.C";
  TString file_CSSq3pim  = "CSSq3pim.C";

  TString file_CSq1pip  = "CSq1pip.C"; // Remember to change nullptr
  TString file_CSq2pip  = "CSq2pip.C"; // to NULL before looking at
  TString file_CSq3pip  = "CSq3pip.C"; // the file.
  TString file_CSq4pip  = "CSq4pip.C";
  TString file_CSq1pim  = "CSq1pim.C";
  TString file_CSq2pim  = "CSq2pim.C";
  TString file_CSq3pim  = "CSq3pim.C";
  TString file_CSq4pim  = "CSq4pim.C";
  //TFile *outFile_CSq1pip = new TFile(outputFile + "CSq1pip.root","RECREATE");
  //TFile *outFile_CSq1pim = new TFile(outputFile + "CSq1pim.root","RECREATE");

  TString file_CSradq1pip  = "CSradq1pip.C"; // Remember to change nullptr
  TString file_CSradq2pip  = "CSradq2pip.C"; // to NULL before looking at
  TString file_CSradq3pip  = "CSradq3pip.C"; // the file.
  TString file_CSradq4pip  = "CSradq4pip.C";
  TString file_CSradq1pim  = "CSradq1pim.C";
  TString file_CSradq2pim  = "CSradq2pim.C";
  TString file_CSradq3pim  = "CSradq3pim.C";
  TString file_CSradq4pim  = "CSradq4pim.C";

  //Set some style options
  gStyle->SetStatY(0.9);
  gStyle->SetStatX(0.9);
  //gStyle->SetOptStat(0);
  gStyle->SetTitleX(0.5);
  gStyle->SetTitleAlign(23);


  //First page of the PDF
  //Divide page into 2x3 grid
  c1->Divide(2,3);
  //Top left
  c1->cd(1);
  el_vz_h_stacked->Draw("hist && nostack");
  c1->cd(1)->BuildLegend(0.7,0.7,0.9,0.9);
  //Top right
  c1->cd(2);
  //c1->cd(2)->SetLogx();
  //c1->cd(2)->SetLogy();
  q2_stacked->Draw("hist && nostack");
  c1->cd(2)->BuildLegend(0.7,0.7,0.9,0.9);
  //Middle left
  c1->cd(3);
  //c1->cd(3)->SetLogx();
  h2_theta_q2_data->Draw("colz"); //h2_theta_q2_data->Write();
  //Middle right
  c1->cd(4);
  //c1->cd(4)->SetLogx();
  h2_theta_q2_MC->Draw("colz"); //h2_theta_q2_MC->Write();
  //Bottom left
  c1->cd(5);
  h2_theta_mom_data->Draw("colz"); //h2_theta_mom_data->Write();
  //Bottom right
  c1->cd(6);
  h2_theta_mom_MC->Draw("colz"); //h2_theta_mom_MC->Write();
  //Start print of page to file
  c1->Print(fileName+"(");
  c1->Clear();

  //Divide page into 2x3 grid
  c1->Divide(2,3);
  //Top left
  c1->cd(1);
  h_q2_MCrad->Draw();
  //Top right
  c1->cd(2);
  h_Ppi_MCrad->Draw();
  //Middle left
  c1->cd(3);
  h_W_MCrad->Draw();
  //Middle right
  c1->cd(4);
  h2_q2_W_MCrad->Draw("colz");
  //Bottom left
  c1->cd(5);
  h2_theta_W_MCrad->Draw("colz");
  //Bottom right
  /*c1->cd(6);
  h2_theta_mom_MC->Draw();*/
  //Start print of page to file
  c1->Print(fileName);
  c1->Clear();

  //Second page
  //Divide page into 2x3 grid
  /*c1->Divide(2,3);
  //Top left
  c1->cd(1);
  W_stacked->Draw("hist && nostack");
  c1->cd(1)->BuildLegend(0.7,0.7,0.9,0.9);
  //Top right
  c1->cd(2);
  //c1->cd(2)->SetLogx();
  c1->cd(2)->SetLogy();
  q2_stacked->Draw("hist && nostack");
  c1->cd(2)->BuildLegend(0.7,0.7,0.9,0.9);
  //Middle right
  c1->cd(3);
  //c1->cd(3)->SetLogy();
  h2_q2_W_data->Draw("colz");
  //Middle left
  c1->cd(4);
  //c1->cd(4)->SetLogy();
  h2_q2_W_MC->Draw("colz");
  //Bottom right
  c1->cd(5);
  h2_theta_W_data->Draw("colz");
  //Bottom left
  c1->cd(6);
  h2_theta_W_MC->Draw("colz");
  //Print page to file
  c1->Print(fileName);
  c1->Clear();


  //Second page
  //Divide page into 2x3 grid
  c1->Divide(2,2);
  //Top left
  c1->cd(1);
  hs_twf->Draw("hist && nostack");
  //Top right
  c1->cd(2);
  h_tf->Draw("colz");
  //Bottom left
  c1->cd(3);
  h_tw->Draw("colz");
  //Bottom right
  c1->cd(4);
  h_twf->Draw("colz");
  //Print page to file
  c1->Print(fileName);
  c1->Clear();*/

  //Second page
  //Divide page into 2x3 grid
  /*c1->Divide(2,1);
  //Top left
  c1->cd(1);
  h_phi_e_opg->Draw();
  //Top right
  c1->cd(2);
  h_phi_pi_opg->Draw();*/
  //Bottom left
  /*c1->cd(3);
  h_opg_tail_omega->Draw("colz");
  //Bottom right
  c1->cd(4);
  h_opg_nota_omega->Draw("colz");*/
  //Print page to file
  //c1->Print(fileName);
  //c1->Clear();

  //Second page
  //Divide page into 2x3 grid
  c1->Divide(2,3);
  //Top left
  c1->cd(1);
  h_W_data->Draw();
  //Top right
  c1->cd(2);
  h_q2_data->Draw();
  //Middle left
  c1->cd(3);
  h_theta_pi_data->Draw();
  //Middle right
  c1->cd(4);
  h_theta_piq_data->Draw();
  //Bottom left
  c1->cd(5);
  h_Ppi_data->Draw();
  //Bottom right
  /*c1->cd(6);
  ->Draw();*/
  //Print page to file
  c1->Print(fileName);
  c1->Clear();

  //Second page
  //Divide page into 2x3 grid
  c1->Divide(2,3);
  //Top left
  c1->cd(1);
  hs_W_data->Draw("hist && nostack");
  //Top right
  c1->cd(2);
  hs_q2_data->Draw("hist && nostack");
  //Middle left
  c1->cd(3);
  hs_theta_pi_data->Draw("hist && nostack");
  //Middle right
  c1->cd(4);
  hs_theta_piq_data->Draw("hist && nostack");
  //Bottom left
  c1->cd(5);
  hs_Ppi_data->Draw("hist && nostack");
  //Bottom right
  /*c1->cd(6);
  ->Draw();*/
  //Print page to file
  c1->Print(fileName);
  c1->Clear();

  //Divide page into 2x2 grid
  c1->Divide(1,2);
  //Top
  c1->cd(1);
  h2_Ppip_thetapip_scin1->Draw("colz");
  //Bottom
  c1->cd(2);
  h2_Ppip_thetapip_scin2->Draw("colz");
  //Print page to file
  c1->Print(fileName);
  c1->Clear();

  //3rd page of the PDF
  //Divide page into 2x4 grid
  c1->Divide(2,4);
  //Top left
  c1->cd(1);
  hs_q2->Draw("hist && nostack");
  c1->cd(1)->BuildLegend(0.7,0.7,0.9,0.9);
  //Top right
  c1->cd(2);
  h_q2_data->Draw(); Double_t x_size = 0.08;
  h_q2_data->GetXaxis()->SetTitleOffset(1.5);
  h_q2_data->GetXaxis()->SetLabelSize(x_size);
  h_q2_data->GetYaxis()->SetLabelSize(x_size);
  //Middle left
  c1->cd(3);
  hs_Ppi->Draw("hist && nostack");
  c1->cd(3)->BuildLegend(0.7,0.7,0.9,0.9);
  //Middle right
  c1->cd(4);
  h_Ppi_data->Draw();
  h_Ppi_data->GetXaxis()->SetTitleOffset(1.5);
  h_Ppi_data->GetXaxis()->SetLabelSize(x_size);
  h_Ppi_data->GetYaxis()->SetLabelSize(x_size);
  //Bottom left
  c1->cd(5);
  //c1->cd(3)->SetLogx();
  hs_theta_piq->Draw("hist && nostack");
  c1->cd(5)->BuildLegend(0.7,0.7,0.9,0.9);
  //Bottom right
  c1->cd(6);
  h_theta_piq_data->Draw();
  h_theta_piq_data->GetXaxis()->SetTitleOffset(1.5);
  h_theta_piq_data->GetXaxis()->SetLabelSize(x_size);
  h_theta_piq_data->GetYaxis()->SetLabelSize(x_size);
  //Bottom left
  c1->cd(7);
  //c1->cd(3)->SetLogx();
  hs_theta_pi->Draw("hist && nostack");
  c1->cd(7)->BuildLegend(0.7,0.7,0.9,0.9);
  //Bottom right
  c1->cd(8);
  h_theta_pi_data->Draw();
  h_theta_pi_data->GetXaxis()->SetTitleOffset(1.5);
  h_theta_pi_data->GetXaxis()->SetLabelSize(x_size);
  h_theta_pi_data->GetYaxis()->SetLabelSize(x_size);
  //Start print of page to file
  c1->Print(fileName);
  c1->Clear();

  //Fourth page
  //Divide page into 2x3 grid
  c1->Divide(2,3);
  //Top left
  c1->cd(1);
  h2_q2_Wpip_data->Draw("colz");
  //Top right
  c1->cd(2);
  h2_q2_Wpim_data->Draw("colz");
  //Middle left
  c1->cd(3);
  h2_Ppip_thetapip_cut->Draw("colz");
  //Middle right
  c1->cd(4);
  h2_Ppim_thetapim_cut->Draw("colz");
  //Bottom left
  c1->cd(5);
  h2_Ppip_thetapipq_cut->Draw("colz");
  //Bottom right
  c1->cd(6);
  h2_Ppim_thetapimq_cut->Draw("colz");
  //Print page to file
  c1->Print(fileName);
  c1->Clear();

  //Fourth page
  //Divide page into 2x3 grid
  //c1->Divide(2,3);
  //Top left
  c1->cd(1);
  h2_q2_Wpip_cut->Draw("colz");
  //Top right
  /*c1->cd(2);
  h2_q2_Wpim_cut->Draw("colz");
  //Middle left
  c1->cd(3);
  h2_Ppip_thetapip_cut->Draw("colz");
  //Middle right
  c1->cd(4);
  h2_Ppim_thetapim_cut->Draw("colz");
  //Bottom left
  c1->cd(5);
  h2_Ppip_thetapipq_cut->Draw("colz");
  //Bottom right
  c1->cd(6);
  h2_Ppim_thetapimq_cut->Draw("colz");*/
  //Print page to file
  c1->Print(fileName);
  c1->Clear();

  //Fifth page
  c1->cd(1);
  h2_q2_Wpim_data->Draw("colz");
  //Print page to file
  c1->Print(fileName);
  c1->Clear();

  //sixth page
  h2_Ppip_thetapip_cut->Draw("colz");
  //Print page to file
  c1->Print(fileName);
  c1->Clear();

  //seventh page
  h2_Ppim_thetapim_cut->Draw("colz");
  //Print page to file
  c1->Print(fileName);
  c1->Clear();

  //eighth page
  h2_Ppip_thetapipq_cut->Draw("colz");
  //Print page to file
  c1->Print(fileName);
  c1->Clear();

  //ninth page
  h2_Ppim_thetapimq_cut->Draw("colz");
  //Print page to file
  c1->Print(fileName);
  c1->Clear();

  //tenth page
  h2_Ppip_thetapipq_data->Draw("colz");
  //Print page to file
  c1->Print(fileName);
  c1->Clear();

  //eleventh page
  h2_Ppim_thetapimq_data->Draw("colz");
  //Print page to file
  c1->Print(fileName);
  c1->Clear();

/*
  //Divide page into 2x2 grid
  c1->Divide(2,2);
  //Top left
  c1->cd(1);
  W_bySector_data_stacked->Draw("hist && nostack");
  c1->cd(1)->BuildLegend(0.7,0.7,0.9,0.9);
  //Top right
  c1->cd(2);
  W_bySector_MC_stacked->Draw("hist && nostack");
  c1->cd(2)->BuildLegend(0.7,0.7,0.9,0.9);
  //Bottom right
  c1->cd(3);
  c1->cd(3)->SetLogx();
  c1->cd(3)->SetLogy();
  q2_bySector_data_stacked->Draw("hist && nostack");
  c1->cd(3)->BuildLegend(0.7,0.7,0.9,0.9);
  //Bottom left
  c1->cd(4);
  c1->cd(4)->SetLogx();
  c1->cd(4)->SetLogy();
  q2_bySector_MC_stacked->Draw("hist && nostack");
  c1->cd(4)->BuildLegend(0.7,0.7,0.9,0.9);
  //Print page to file
  c1->Print(fileName);
  c1->Clear();

*/
  //Divide page into 2x6 grid
  c1->Divide(2,6);
  c1->cd(1);
  h2_theta_W_S1_data->Draw("colz");
  c1->cd(2);
  h2_theta_W_S1_MC->Draw("colz");
  c1->cd(3);
  h2_theta_W_S2_data->Draw("colz");
  c1->cd(4);
  h2_theta_W_S2_MC->Draw("colz");
  c1->cd(5);
  h2_theta_W_S3_data->Draw("colz");
  c1->cd(6);
  h2_theta_W_S3_MC->Draw("colz");
  c1->cd(7);
  h2_theta_W_S4_data->Draw("colz");
  c1->cd(8);
  h2_theta_W_S4_MC->Draw("colz");
  c1->cd(9);
  h2_theta_W_S5_data->Draw("colz");
  c1->cd(10);
  h2_theta_W_S5_MC->Draw("colz");
  c1->cd(11);
  h2_theta_W_S6_data->Draw("colz");
  c1->cd(12);
  h2_theta_W_S6_MC->Draw("colz");
  //Print page to file
  c1->Print(fileName);
  c1->Clear();/*


  //Third page
  //Divide page into 2x2 grid
  c1->Divide(2,2);
  //Top left
  c1->cd(1);
  omega_stacked->Draw("hist && nostack");
  c1->cd(1)->BuildLegend(0.7,0.7,0.9,0.9);
  //mid left
  c1->cd(3);
  h2_theta_omega_data->Draw("colz");
  //mid right
  c1->cd(4);
  h2_theta_omega_MC->Draw("colz");
  //Print page to file
  c1->Print(fileName);
  c1->Clear();

  //Seventh page
  //Divide page into 2x2 grid
  c1->Divide(2,2);
  //Top left
  c1->cd(1);
  q2_stacked->Draw("hist && nostack");
  //c1->cd(1)->SetLogy();
  c1->cd(1)->BuildLegend(0.6,0.6,0.9,0.9);
  //Top right
  c1->cd(2);
  W_stacked->Draw("hist && nostack");
  //c1->cd(2)->SetLogy();
  c1->cd(2)->BuildLegend(0.1,0.6,0.35,0.9);
  //Bottom left
  c1->cd(3);
  q2_stacked->Draw("hist && nostack");
  c1->cd(3)->SetLogy();
  c1->cd(3)->BuildLegend(0.6,0.6,0.9,0.9);
  //Bottom right
  c1->cd(4);
  W_stacked->Draw("hist && nostack");
  c1->cd(4)->SetLogy();
  c1->cd(4)->BuildLegend(0.1,0.6,0.35,0.9);
  //Print page to file
  c1->Print(fileName);
  c1->Clear();*/
  if(bAllEvents == false){MC_evts[0] = 1.812E+7; MC_rad_evts[0] = 0.972E+7;}
  c1->Divide(TPIQi,Pbini);
  //-------------------------------------------------Write to file------------------------------------------//
    Int_t line = 1, size = 3, marc = 2, mark_d = 8;
    Double_t label_size = 0.11;//0.06 0.09
    for(Int_t iCharge = 0; iCharge < Nrec; iCharge++)/// Acceptance Corrections
    {
        for(Int_t iQ2 = 0; iQ2 < Q2bini; iQ2++)
        {
            canvas = 1;
            Double_t dCSq[2][Wbini];
            Double_t dCSt[TPIQi][2][Wbini];
            for(Int_t iW = 0; iW < Wbini; iW++)
            {
                dCSq[0][iW] = 0.; dCSq[1][iW] = 0.;
                for(Int_t iT = 0; iT < TPIQi; iT++)
                {
                    dCSt[iT][0][iW] = 0.; dCSt[iT][1][iW] = 0.;
                }
            }
            sprintf(c_NAME,"h_CSrecQ%d%d",iCharge,iQ2);
            TH1D * h_CSrecQ = new TH1D(c_NAME,"",Wbini,CS_Wmin,CS_Wmax);
            sprintf(c_NAME,"h_CSgenQ%d%d",iCharge,iQ2);
            TH1D * h_CSgenQ = new TH1D(c_NAME,"",Wbini,CS_Wmin,CS_Wmax);
            for(Int_t iPmom = 0; iPmom < Pbini; iPmom++)
            {
                Int_t iPcanv = (Pbini - 1) - iPmom;//Needed to get plot positions correct for P_pi
                Double_t dCSp[2][Wbini];
                for(Int_t iW = 0; iW < Wbini; iW++)
                {
                    dCSp[0][iW] = 0.; dCSp[1][iW] = 0.;
                }
                sprintf(c_NAME,"h_CSrecP%d%d%d",iPcanv,iCharge,iQ2);
                TH1D * h_CSrecP = new TH1D(c_NAME,"",Wbini,CS_Wmin,CS_Wmax);
                sprintf(c_NAME,"h_CSgenP%d%d%d",iPcanv,iCharge,iQ2);
                TH1D * h_CSgenP = new TH1D(c_NAME,"",Wbini,CS_Wmin,CS_Wmax);
                for(Int_t iThPiQ = 0; iThPiQ < TPIQi; iThPiQ++)
                {
                    dT = DeltaThetaPiQBin(ThetaPiQ[iThPiQ],ThetaPiQ[iThPiQ+1]);
                    Double_t dPmom = dP[iPcanv];
                    if(bQorNot==false)
                    {
                        dT = DeltaThetaPiQBin(ThetaPi[iCharge][iThPiQ],ThetaPi[iCharge][iThPiQ+1]);
                        dPmom = dPt[iCharge][iPmom+1] - dPt[iCharge][iPmom];
                    }
                    bin_vol = dPmom * dT * dQ2[iQ2] * dW;
                    c1->cd(canvas); canvas++; //Needed to get plot positions correct
                    sprintf(c_NAME,"hs_CS%d%d%d%d",iPcanv,iThPiQ,iCharge,iQ2);
                    THStack * hs_CS = new THStack(c_NAME,"");
                    sprintf(c_NAME,"h_CSrec%d%d%d%d",iPcanv,iThPiQ,iCharge,iQ2);
                    TH1D * h_CSrec = new TH1D(c_NAME,"",Wbini,CS_Wmin,CS_Wmax);
                    sprintf(c_NAME,"h_CSgen%d%d%d%d",iPcanv,iThPiQ,iCharge,iQ2);
                    TH1D * h_CSgen = new TH1D(c_NAME,"",Wbini,CS_Wmin,CS_Wmax);
                    //gPad->SetLogy();
                    for(Int_t iMod = iGEN; iMod < iDATA; iMod++)
                    {
                        if(iMod == iGENIEr) continue;
                        sprintf(c_NAME,"h_CSrg%d%d%d%d%d",iPcanv,iThPiQ,iCharge,iQ2,iMod);
                        TH1D * h_CSrg = new TH1D(c_NAME,"",Wbini,CS_Wmin,CS_Wmax);
                        //h_CSrg->Sumw2(); //https://root-forum.cern.ch/t/array-of-histograms/13148/2
                        for(Int_t iW = 0; iW < Wbini; iW++)
                        {
                            h_CSrg->SetBinContent(iW+1, y_sum[iPcanv][iThPiQ][iCharge][iMod][iQ2][iW]);
                            dCSq[iMod-3][iW] += y_sum[iPcanv][iThPiQ][iCharge][iMod][iQ2][iW];
                            //cout <<iMod<<" "<<iW <<"\t"<< dCSq[iMod-3][iW] <<" += "<< y_sum[iPcanv][iThPiQ][iCharge][iMod][iQ2][iW]<<endl;
                            dCSp[iMod-3][iW] += y_sum[iPcanv][iThPiQ][iCharge][iMod][iQ2][iW];
                            dCSt[iThPiQ][iMod-3][iW] += y_sum[iPcanv][iThPiQ][iCharge][iMod][iQ2][iW];
                        }   // Note: for an array bin # of iW, the ROOT bin # is iW+1
                        if(iMod==iGENIE)//GENIE reconstructed
                        {
                            h_CSrg->Scale(totalCS_MC[iGENmod] * ub2nb / (MC_evts[iGENmod] * bin_vol));//kRed
                            line = 1; mark = 32; size = 1; marc = 2; h_CSrg->SetLineColor(2);
                            h_CSrec = h_CSrg;
                            //cout << "Rec " << h_CSrg->GetBinContent(25) <<"\t"<< y_sum[iPcanv][iThPiQ][iCharge][iMod][iQ2][24] << endl;
                        }
                        else if(iMod==iGEN)//GENIE generated
                        {
                            h_CSrg->Scale(totalCS_MC[iGENmod] * ub2nb * (5./6.) / (MC_evts[iGENmod] * bin_vol));//kRed
                            line = 1; mark = 32; size = 1; marc = 4; h_CSrg->SetLineColor(3);
                            h_CSgen = h_CSrg; hr_ARC[iPcanv][iThPiQ][iCharge][iQ2] = h_CSrg;
                            //cout << "Gen " << h_CSrg->GetBinContent(25) <<"\t"<< y_sum[iPcanv][iThPiQ][iCharge][iMod][iQ2][24] << endl;
                        }
                        h_CSrg->SetLineStyle(line);
                        h_CSrg->SetMarkerStyle(mark);
                        h_CSrg->SetMarkerSize(size);
                        h_CSrg->SetMarkerColor(marc);
                        hs_CS->Add(h_CSrg);
                    }
                    hr_AC[iPcanv][iThPiQ][iCharge][iQ2]->Sumw2();
                    hr_AC[iPcanv][iThPiQ][iCharge][iQ2]->Divide(h_CSgen,h_CSrec);
                    hs_CS->Draw("nostack,e1p");
                    hs_CS->Write("nostack,e1p");
                    hs_CS->GetYaxis()->SetLabelSize(label_size);//
                    hs_CS->GetYaxis()->SetNdivisions(4);
                    hs_CS->GetXaxis()->SetNdivisions(4);
                    gStyle->SetEndErrorSize(10);
                    gStyle->SetErrorX(0);
                    if(iPcanv == 0 && iThPiQ == 0 && iCharge == 0)
                    {
                        hs_CS->GetXaxis()->SetTitleOffset(1.5);
                        hs_CS->GetXaxis()->SetLabelSize(label_size);
                    }
                    else if(iPcanv == 0 && iThPiQ == 0 && iCharge == 1)
                    {
                        hs_CS->GetXaxis()->SetTitleOffset(1.5);
                        hs_CS->GetXaxis()->SetLabelSize(label_size);
                    }
                    else if(iPcanv == 0)
                    {
                        hs_CS->GetXaxis()->SetTitleOffset(1.5);
                        hs_CS->GetXaxis()->SetLabelSize(label_size);
                    }
                    else if(iThPiQ == 2 && iPcanv == Pbini - 1 && iQ2==0)
                    {
                        if(iCharge == 0){hs_CS->SetTitle("GENIE gen and recon pi+   0.70 < Q2 < 1.0;;");}
                        else if(iCharge == 1){hs_CS->SetTitle("GENIE gen and recon pi-   0.70 < Q2 < 1.0;;");}
                    }
                    else if(iThPiQ == 2 && iPcanv == Pbini - 1 && iQ2==1)
                    {
                        if(iCharge == 0){hs_CS->SetTitle("GENIE gen and recon pi+   1.0 < Q2 < 1.4;;");}
                        else if(iCharge == 1){hs_CS->SetTitle("GENIE gen and recon pi-   1.0 < Q2 < 1.4;;");}
                    }
                    else if(iThPiQ == 2 && iPcanv == Pbini - 1 && iQ2==2)
                    {
                        if(iCharge == 0){hs_CS->SetTitle("GENIE gen and recon pi+   1.4 < Q2 < 1.9;;");}
                        else if(iCharge == 1){hs_CS->SetTitle("GENIE gen and recon pi-   1.4 < Q2 < 1.9;;");}
                    }
                    else if(iThPiQ == 2 && iPcanv == Pbini - 1 && iQ2==3)
                    {
                        if(iCharge == 0){hs_CS->SetTitle("GENIE gen and recon pi+   1.9 < Q2 < 2.5;;");}
                        else if(iCharge == 1){hs_CS->SetTitle("GENIE gen and recon pi-   1.9 < Q2 < 2.5;;");}
                    }
                }
                ///Calculate 3D P_pi AC
                hr_ACp[iPcanv][iCharge][iQ2]->Sumw2();
                for(Int_t iW = 0; iW < Wbini; iW++)
                {
                    /*if(dCSp[iGENIE-3][iW] > 0. && dCSp[iGEN-3][iW] > 0.){hr_ACp[iPcanv][iCharge][iQ2]->Sumw2();
                        hr_ACp[iPcanv][iCharge][iQ2]->SetBinContent(iW+1, dCSp[iGEN-3][iW] * (5./6.) / dCSp[iGENIE-3][iW]);
                    }*/
                    h_CSrecP->SetBinContent(iW+1, dCSp[iGENIE-3][iW] * totalCS_MC[iGENmod] * ub2nb / (MC_evts[iGENmod]));
                    h_CSgenP->SetBinContent(iW+1, dCSp[iGEN-3][iW] * totalCS_MC[iGENmod] * ub2nb * (5./6.) / (MC_evts[iGENmod]));
                    /*if(h_CSrecP->GetBinContent(iW+1) > 0. && h_CSgenP->GetBinContent(iW+1) > 0.){
                        hr_ACp[iPcanv][iCharge][iQ2]->SetBinContent(iW+1, h_CSgenP->GetBinContent(iW+1)/h_CSrecP->GetBinContent(iW+1));
                    }*/
                }   // Note: for an array bin # of iW, the ROOT bin # is iW+1
                hr_ACp[iPcanv][iCharge][iQ2]->Divide(h_CSgenP,h_CSrecP);
                hr_ARCp[iPcanv][iCharge][iQ2] = h_CSgenP;
            }
            ///Calculate 3D Theta_piq AC
            for(Int_t iThPiQ = 0; iThPiQ < TPIQi; iThPiQ++)
            {
                hr_ACt[iThPiQ][iCharge][iQ2]->Sumw2();
                sprintf(c_NAME,"h_CSrecT%d%d%d",iThPiQ,iCharge,iQ2);
                TH1D * h_CSrecT = new TH1D(c_NAME,"",Wbini,CS_Wmin,CS_Wmax);
                sprintf(c_NAME,"h_CSgenT%d%d%d",iThPiQ,iCharge,iQ2);
                TH1D * h_CSgenT = new TH1D(c_NAME,"",Wbini,CS_Wmin,CS_Wmax);
                for(Int_t iW = 0; iW < Wbini; iW++)
                {
                    /*if(dCSt[iThPiQ][iGENIE-3][iW] > 0. && dCSt[iThPiQ][iGEN-3][iW] > 0.){hr_ACt[iThPiQ][iCharge][iQ2]->Sumw2();
                        hr_ACt[iThPiQ][iCharge][iQ2]->SetBinContent(iW+1, dCSt[iThPiQ][iGEN-3][iW] * (5./6.) / dCSt[iThPiQ][iGENIE-3][iW]);
                    }*/
                    h_CSrecT->SetBinContent(iW+1, dCSt[iThPiQ][iGENIE-3][iW] * totalCS_MC[iGENmod] * ub2nb / (MC_evts[iGENmod]));
                    h_CSgenT->SetBinContent(iW+1, dCSt[iThPiQ][iGEN-3][iW] * totalCS_MC[iGENmod] * ub2nb * (5./6.) / (MC_evts[iGENmod]));
                    /*if(h_CSrecT->GetBinContent(iW+1) > 0. && h_CSgenT->GetBinContent(iW+1) > 0.){
                        hr_ACt[iThPiQ][iCharge][iQ2]->SetBinContent(iW+1, h_CSgenT->GetBinContent(iW+1)/h_CSrecT->GetBinContent(iW+1));
                    }*/
                }   // Note: for an array bin # of iW, the ROOT bin # is iW+1
                hr_ACt[iThPiQ][iCharge][iQ2]->Divide(h_CSgenT,h_CSrecT);
                hr_ARCt[iThPiQ][iCharge][iQ2] = h_CSgenT;
            }
            ///Calculate 2D Q2 AC
            hr_ACq[iCharge][iQ2]->Sumw2();
            for(Int_t iW = 0; iW < Wbini; iW++)
            {
                /*if(dCSq[iGENIE-3][iW] > 0. && dCSq[iGEN-3][iW] > 0.){hr_ACq[iCharge][iQ2]->Sumw2();
                    hr_ACq[iCharge][iQ2]->SetBinContent(iW+1, dCSq[iGEN-3][iW] * (5./6.) / dCSq[iGENIE-3][iW]);
                }*/
                h_CSrecQ->SetBinContent(iW+1, dCSq[iGENIE-3][iW] * totalCS_MC[iGENmod] * ub2nb / (MC_evts[iGENmod]));
                h_CSgenQ->SetBinContent(iW+1, dCSq[iGEN-3][iW] * totalCS_MC[iGENmod] * ub2nb * (5./6.) / (MC_evts[iGENmod]));
                /*if(h_CSrecQ->GetBinContent(iW+1) > 0. && h_CSgenQ->GetBinContent(iW+1) > 0.){
                    hr_ACq[iCharge][iQ2]->SetBinContent(iW+1, h_CSgenQ->GetBinContent(iW+1)/h_CSrecQ->GetBinContent(iW+1));
                }*/
            }   // Note: for an array bin # of iW, the ROOT bin # is iW+1
            hr_ACq[iCharge][iQ2]->Divide(h_CSgenQ,h_CSrecQ);
            hr_ARCq[iCharge][iQ2] = h_CSgenQ;
            ///Print 4D GENIE gen and rec cross sections
            c1->Print(fileName);
        }
    }
    for(Int_t iCharge = 0; iCharge < Nrec; iCharge++)
    {
        for(Int_t iQ2 = 0; iQ2 < Q2bini; iQ2++)
        {
            canvas = 1;
            for(Int_t iPmom = 0; iPmom < Pbini; iPmom++)
            {
                Int_t iPcanv = (Pbini - 1) - iPmom;//Needed to get plot positions correct for P_pi
                for(Int_t iThPiQ = 0; iThPiQ < TPIQi; iThPiQ++)
                {
                    c1->cd(canvas); canvas++; //Needed to get plot positions correct
                    gPad->SetLogy(0);
                    hr_AC[iPcanv][iThPiQ][iCharge][iQ2]->Draw("e1p");
                    hr_AC[iPcanv][iThPiQ][iCharge][iQ2]->Write("e1p");
                    hr_AC[iPcanv][iThPiQ][iCharge][iQ2]->GetYaxis()->SetNdivisions(4);
                    hr_AC[iPcanv][iThPiQ][iCharge][iQ2]->GetYaxis()->SetLabelSize(label_size);
                    hr_AC[iPcanv][iThPiQ][iCharge][iQ2]->SetMaximum(20.);
                    hr_AC[iPcanv][iThPiQ][iCharge][iQ2]->SetMinimum(0.);
                    hr_AC[iPcanv][iThPiQ][iCharge][iQ2]->GetXaxis()->SetNdivisions(4);
                    hr_AC[iPcanv][iThPiQ][iCharge][iQ2]->GetXaxis()->SetLabelSize(label_size);
                    if(iThPiQ == 2 && iPcanv == Pbini - 1 && iQ2==0)
                    {
                        if(iCharge == 0){hr_AC[iPcanv][iThPiQ][iCharge][iQ2]->SetTitle("GENIE gen/recon pi+   0.70 < Q2 < 1.0;;");}
                        else if(iCharge == 1){hr_AC[iPcanv][iThPiQ][iCharge][iQ2]->SetTitle("GENIE gen/recon pi-   0.70 < Q2 < 1.0;;");}
                    }
                    else if(iThPiQ == 2 && iPcanv == Pbini - 1 && iQ2==1)
                    {
                        if(iCharge == 0){hr_AC[iPcanv][iThPiQ][iCharge][iQ2]->SetTitle("GENIE gen/recon pi+   1.0 < Q2 < 1.4;;");}
                        else if(iCharge == 1){hr_AC[iPcanv][iThPiQ][iCharge][iQ2]->SetTitle("GENIE gen/recon pi-   1.0 < Q2 < 1.4;;");}
                    }
                    else if(iThPiQ == 2 && iPcanv == Pbini - 1 && iQ2==2)
                    {
                        if(iCharge == 0){hr_AC[iPcanv][iThPiQ][iCharge][iQ2]->SetTitle("GENIE gen/recon pi+   1.4 < Q2 < 1.9;;");}
                        else if(iCharge == 1){hr_AC[iPcanv][iThPiQ][iCharge][iQ2]->SetTitle("GENIE gen/recon pi-   1.4 < Q2 < 1.9;;");}
                    }
                    else if(iThPiQ == 2 && iPcanv == Pbini - 1 && iQ2==3)
                    {
                        if(iCharge == 0){hr_AC[iPcanv][iThPiQ][iCharge][iQ2]->SetTitle("GENIE gen/recon pi+   1.9 < Q2 < 2.5;;");}
                        else if(iCharge == 1){hr_AC[iPcanv][iThPiQ][iCharge][iQ2]->SetTitle("GENIE gen/recon pi-   1.9 < Q2 < 2.5;;");}
                    }
                }
            }
            c1->Print(fileName);
        }
    }
  c1->Clear();
  c1->Divide(Pbini,Q2bini);
    ///Check 3D ratios
    for(Int_t iCharge = 0; iCharge < Nrec; iCharge++)
    {
        canvas = 1;
        for(Int_t iQ2 = 0; iQ2 < Q2bini; iQ2++)
        {
            for(Int_t iPmom = 0; iPmom < Pbini; iPmom++)
            {Int_t iQcanv = (Q2bini - 1) - iQ2;//Needed to get plot positions correct for Q2
                Int_t iPcanv = iPmom;//Needed to get plot positions correct for P_pi
                c1->cd(canvas); canvas++; //Needed to get plot positions correct
                //gPad->SetLogy(0);
                hr_ACp[iPcanv][iCharge][iQcanv]->Draw("e1p");
                hr_ACp[iPcanv][iCharge][iQcanv]->Write("e1p");
                hr_ACp[iPcanv][iCharge][iQcanv]->GetYaxis()->SetNdivisions(4);
                hr_ACp[iPcanv][iCharge][iQcanv]->GetYaxis()->SetLabelSize(label_size);
                hr_ACp[iPcanv][iCharge][iQcanv]->SetMaximum(20.);
                hr_ACp[iPcanv][iCharge][iQcanv]->SetMinimum(0.);
                hr_ACp[iPcanv][iCharge][iQcanv]->GetXaxis()->SetNdivisions(4);
                hr_ACp[iPcanv][iCharge][iQcanv]->GetXaxis()->SetLabelSize(label_size);
                gStyle->SetEndErrorSize(10);
                gStyle->SetErrorX(0);
                if(iQcanv == 3 && iPcanv==2)
                {
                    if(iCharge == 0){hr_ACp[iPcanv][iCharge][iQcanv]->SetTitle("3D GENIE gen/recon pi+   Q2 and P_pi;;");}
                    else if(iCharge == 1){hr_ACp[iPcanv][iCharge][iQcanv]->SetTitle("3D GENIE gen/recon pi-   Q2 and P_pi;;");}
                }
            }
        }
        c1->Print(fileName);
    }
  c1->Clear();
  c1->Divide(TPIQi,Q2bini);
    for(Int_t iCharge = 0; iCharge < Nrec; iCharge++)
    {canvas = 1;
        for(Int_t iQ2 = 0; iQ2 < Q2bini; iQ2++)
        {
            Int_t iQcanv = (Q2bini - 1) - iQ2;//Needed to get plot positions correct for Q2
            for(Int_t iThPiQ = 0; iThPiQ < TPIQi; iThPiQ++)
            {
                c1->cd(canvas); canvas++; //Needed to get plot positions correct
                gPad->SetLogy(0);
                hr_ACt[iThPiQ][iCharge][iQcanv]->Draw("e1p");
                hr_ACt[iThPiQ][iCharge][iQcanv]->Write("e1p");
                hr_ACt[iThPiQ][iCharge][iQcanv]->GetYaxis()->SetNdivisions(4);
                hr_ACt[iThPiQ][iCharge][iQcanv]->GetYaxis()->SetLabelSize(label_size);
                hr_ACt[iThPiQ][iCharge][iQcanv]->SetMaximum(20.);
                hr_ACt[iThPiQ][iCharge][iQcanv]->SetMinimum(0.);
                hr_ACt[iThPiQ][iCharge][iQcanv]->GetXaxis()->SetNdivisions(4);
                hr_ACt[iThPiQ][iCharge][iQcanv]->GetXaxis()->SetLabelSize(label_size);
                gStyle->SetEndErrorSize(10);
                gStyle->SetErrorX(0);
                if(iQcanv == 3 && iThPiQ==2)
                {
                    if(iCharge == 0){hr_ACt[iThPiQ][iCharge][iQcanv]->SetTitle("3D GENIE gen/recon pi+   Q2 and Theta_piq;;");}
                    else if(iCharge == 1){hr_ACt[iThPiQ][iCharge][iQcanv]->SetTitle("3D GENIE gen/recon pi-   Q2 and Theta_piq;;");}
                }
            }
        }
        c1->Print(fileName);
    }
  //Print page to file
  c1->Clear();
    c1->Divide(2,2);
    ///Check 2D ratios
    for(Int_t iCharge = 0; iCharge < Nrec; iCharge++)
    {
        canvas = 1;
        for(Int_t iQ2 = 0; iQ2 < Q2bini; iQ2++)
        {Int_t iQcanv = (Q2bini - 1) - iQ2;//Needed to get plot positions correct for Q2
            c1->cd(canvas); canvas++; //Needed to get plot positions correct
            gPad->SetLogy(0);
            hr_ACq[iCharge][iQcanv]->Draw("e1p");
            hr_ACq[iCharge][iQcanv]->Write("e1p");
            hr_ACq[iCharge][iQcanv]->GetYaxis()->SetNdivisions(4);
            hr_ACq[iCharge][iQcanv]->GetYaxis()->SetLabelSize(label_size);
            hr_ACq[iCharge][iQcanv]->SetMaximum(20.);
            hr_ACq[iCharge][iQcanv]->SetMinimum(0.);
            hr_ACq[iCharge][iQcanv]->GetXaxis()->SetNdivisions(4);
            hr_ACq[iCharge][iQcanv]->GetXaxis()->SetLabelSize(label_size);
            gStyle->SetEndErrorSize(10);
            gStyle->SetErrorX(0);
            if(iCharge == 0 && iQcanv==3){hr_ACq[iCharge][iQcanv]->SetTitle("2D GENIE gen/recon pi+;;");}
            else if(iCharge == 1 && iQcanv==3){hr_ACq[iCharge][iQcanv]->SetTitle("2D GENIE gen/recon pi-;;");}
            //cout <<"2D ACq: "<< hr_ACq[iCharge][iQ2]->GetBinContent(25)<<endl;
        }
        c1->Print(fileName);
    }
  //Print page to file
  c1->Clear();

  //Divide page into 2x2 grid
  c1->Divide(1,2);
  //Top
  c1->cd(1);
  h_pip_chi2pid->Draw();
  //Bottom
  c1->cd(2);
  h_pim_chi2pid->Draw();
  //Print page to file
  c1->Print(fileName);
  c1->Clear();

  /////////////////////////////////////////////////////////////////////////////////
  //2D Histograms of topologies
  /////////////////////////////////////////////////////////////////////////////////

  //2nd-to-Last page (Cross Sections by Sector)
  //Divide page into 5x5 grid
  //c1->Divide(6,5);
  c1->Divide(TPIQi,Pbini);
  line = 1, size = 3, marc = 2;
  label_size = 0.09; Double_t var_sum = 0., var_num = 0.;
  for(Int_t iCharge = 0; iCharge < Nrec; iCharge++)/// SECTOR CROSS SECTIONS
  {
      for(Int_t iQ2 = 0; iQ2 < Q2bini; iQ2++)
      {
          canvas = 1;
          for(Int_t iPmom = 0; iPmom < Pbini; iPmom++)
          {
              Int_t iPcanv = (Pbini - 1) - iPmom;//Needed to get plot positions correct for P_pi
              for(Int_t iThPiQ = 0; iThPiQ < TPIQi; iThPiQ++)
              {
                dT = DeltaThetaPiQBin(ThetaPiQ[iThPiQ],ThetaPiQ[iThPiQ+1]);
                Double_t dPmom = dP[iPcanv];
                if(bQorNot==false)
                {
                    dT = DeltaThetaPiQBin(ThetaPi[iCharge][iThPiQ],ThetaPi[iCharge][iThPiQ+1]);
                    dPmom = dPt[iCharge][iPmom+1] - dPt[iCharge][iPmom];
                }
                bin_vol = dPmom * dT * dQ2[iQ2] * dW;
                c1->cd(canvas); canvas++; //Needed to get plot positions correct
                //Double_t y_ave[iMod][iW] = 0.;
                //Double_t var_stat[iMod][iW] = 0.;
                sprintf(c_NAME,"hs_CSS%d%d%d%d",iPcanv,iThPiQ,iCharge,iQ2);
                THStack * hs_CSS = new THStack(c_NAME,"");
                sprintf(c_NAME,"h_CSdataS%d%d%d%d",iPcanv,iThPiQ,iCharge,iQ2);
                TH1D * h_CSdataS = new TH1D(c_NAME,"",Wbini,CS_Wmin,CS_Wmax);
                h_CSdataS->Sumw2();
                for(Int_t iSector = 0; iSector < sectors; iSector++)
                {
                    ///gPad->SetLogy(); //cout << " s" << iSector;
                    //if(iSector==3) continue;
                    for(Int_t iMod = 0; iMod < models; iMod++)
                    {//cout << "f" << iMod;[iPmom][iThPiQ][pm][s][f][iQ2][iW]
                        sprintf(c_NAME,"h_CSS%d%d%d%d%d%d",iPcanv,iThPiQ,iCharge,iSector,iMod,iQ2);
                        TH1D * h_CSS = new TH1D(c_NAME,"",Wbini,CS_Wmin,CS_Wmax);
                        h_CSS->Sumw2();
                        for(Int_t iW = 0; iW < Wbini; iW++)
                        {
                            /*if(iW == 25)
                                {cout << endl <<"Before "<< iMod <<"\t"<< y_sum[iPcanv][iThPiQ][iCharge][iMod][iQ2][iW];}
                            if(iMod!=iGEN || iMod!=iGENIE)
                            {
                                y_sum[iPcanv][iThPiQ][iCharge][iMod][iQ2][iW] += y_i[iPcanv][iThPiQ][iCharge][iSector][iMod][iQ2][iW];
                            }
                            if(iW == 25)
                                {cout << endl <<"After  "<< iMod <<"\t"<< y_sum[iPcanv][iThPiQ][iCharge][iMod][iQ2][iW];}*/
                            ///Scale h_CSS bins by total CS for onepigen
                            if((iMod==iOPGr || iMod==iOPGn) && y_i[iPcanv][iThPiQ][iCharge][iSector][iMod][iQ2][iW]!=0)
                            {
                                Double_t opg_error = tCS[iPcanv][iThPiQ][iCharge][iSector][iMod][iQ2][iW]/sqrt(y_i[iPcanv][iThPiQ][iCharge][iSector][iMod][iQ2][iW]);
                                h_CSS->SetBinContent(iW+1, tCS[iPcanv][iThPiQ][iCharge][iSector][iMod][iQ2][iW]);
                                h_CSS->SetBinError(iW+1, opg_error);
                            }
                            else/// if(iMod==iDATA || iMod==iGENIE)///Do not scale h_CSS bins by total CS for non-onepigen (or bins have zero events)
                            {
                                h_CSS->SetBinContent(iW+1, y_i[iPcanv][iThPiQ][iCharge][iSector][iMod][iQ2][iW]);
                                h_CSS->SetBinError(iW+1, sqrt(y_i[iPcanv][iThPiQ][iCharge][iSector][iMod][iQ2][iW]));
                                if(iMod==iDATA && iSector==0)
                                {
                                    h_CSdataS->SetBinContent(iW+1, y_i[iPcanv][iThPiQ][iCharge][iSector][iMod][iQ2][iW]);
                                }
                            }
                            /*if(iW == 25 && iPmom ==2 && iThPiQ==2 && iQ2==0)
                                {cout << iMod<<" "<<iSector<<" "<< h_CSS->GetBinContent(iW+1) <<" +- "<< h_CSS->GetBinError(iW+1)<<endl;}*/
                        }   // Note: for an array bin # of iW, the ROOT bin # is iW+1
                        //if (h_CSS->GetSumw2N() == 0) h_CSS->Sumw2(kTRUE);
                        if(iMod==iDATA)//Data
                        {
                            h_CSS->Scale(cm2ub * ub2nb / (L * bin_vol));
                            line = 1; mark = 8;  size = 2; marc = 1 + iSector;
                            /*if(iSector==0)/// For Sector cross section ratios
                            {/// Remember: cannot use this and regular cross sections
                                h_CSdataS->Scale(cm2ub * ub2nb / (L * bin_vol));
                                //h_CSdataS = h_CSS;
                            }*/
                            /*else if(iSector==2) //do not need for ratios
                            {
                                h_CSdataS->Scale(1*Double_t(iSector));
                            }*/
                            //h_CSS->Divide(h_CSS,h_CSdataS);/// end of ratios
                            //hs_CSS->Add(h_CSdataS);
                            hs_CSS->Add(h_CSS);
                        }
                        else if(iMod==iGENIE)//GENIE reconstructed
                        {
                            h_CSS->Scale(totalCS_MC[iGENmod] * ub2nb / (MC_evts[iGENmod] * bin_vol));
                        }
                        else if(iMod==iGEN)//GENIE generated
                        {
                            h_CSS->Scale(totalCS_MC[iGENmod] * ub2nb * (5./6.) / (MC_evts[iGENmod] * bin_vol));
                        }
                        else if(iMod==iOPGr && iCharge==0)//onepigen
                        {
                            h_CSS->Scale(ub2nb / (opg_gen_events[0] * bin_vol));
                        }
                        else if(iMod==iOPGr && iCharge==1)//onepigen pi-
                        {
                            h_CSS->Scale(ub2nb / (opg_gen_events[1] * bin_vol));/// <-------------------<<
                        }
                        else if(iMod==iOPGn && iCharge==0)//onepigen pi+ norad
                        {
                            h_CSS->Scale(ub2nb / (opg_gen_events[2] * bin_vol));
                        }
                        else if(iMod==iOPGn && iCharge==1)//onepigen pi- norad
                        {
                            h_CSS->Scale(ub2nb / (opg_gen_events[3] * bin_vol));
                        }
                        h_CSS->SetLineColor(marc);
                        h_CSS->SetLineStyle(line);
                        h_CSS->SetMarkerStyle(mark);
                        h_CSS->SetMarkerSize(size);
                        h_CSS->SetMarkerColor(marc);
                        gStyle->SetEndErrorSize(10);
                        gStyle->SetErrorX(0);
                        for(Int_t iW = 0; iW < Wbini; iW++)
                        {
                          //Double_t stat_error = cm2ub * ub2nb * sqrt(y_i[iPcanv][iThPiQ][iCharge][iSector][iMod][iQ2][iW]) / (L * bin_vol);
                          //h_CSS->SetBinError(iW+1, stat_error);
                          //if(iW == 25 && iPmom ==2 && iThPiQ==2 && iQ2==0 && iMod==iDATA)
                          //    {cout << iMod<<" "<<iSector<<" "<< h_CSS->GetBinContent(iW+1) <<" +- "<< h_CSS->GetBinError(iW+1)<<endl;}
                          y_i[iPcanv][iThPiQ][iCharge][iSector][iMod][iQ2][iW] = h_CSS->GetBinContent(iW+1);
                          y_ave[iPcanv][iThPiQ][iCharge][iMod][iQ2][iW] += y_i[iPcanv][iThPiQ][iCharge][iSector][iMod][iQ2][iW]/CS_Sec;//# of sectors
                          var_stat[iPcanv][iThPiQ][iCharge][iSector][iMod][iQ2][iW]  += pow(h_CSS->GetBinError(iW+1),2.);
                          /*if(iW == 25 && iPmom ==2 && iThPiQ==2 && iQ2==0 && iMod==iDATA)
                              {cout <<" y_i= "<< y_i[iPcanv][iThPiQ][iCharge][iSector][iMod][iQ2][iW];
                               cout <<"\t y_ave= "<< y_ave[iPcanv][iThPiQ][iCharge][iMod][iQ2][iW];
                               cout <<"\t var_stat= "<< var_stat[iPcanv][iThPiQ][iCharge][iSector][iMod][iQ2][iW] <<endl;}*/
                          if(iSector!=0 && iMod<iGEN)///Set tCS to sum of total CS for all sectors
                          {
                              tCS[iPcanv][iThPiQ][iCharge][0][iMod][iQ2][iW] += tCS[iPcanv][iThPiQ][iCharge][iSector][iMod][iQ2][iW];
                          }
                        } // Note: for an array bin # of iW, the ROOT bin # is iW+1
                    }
                }
                    hs_CSS->Draw("nostack,e1p");
                    hs_CSS->GetYaxis()->SetLabelSize(label_size);
                    hs_CSS->GetXaxis()->SetLabelSize(label_size);
                    hs_CSS->GetYaxis()->SetNdivisions(4);
                    hs_CSS->GetXaxis()->SetNdivisions(5);
                    hs_CSS->SetTitle(";;");//";W (GeV);CS");
                    //hs_CSS->Write();
                for(Int_t iMod = 0; iMod < models; iMod++)
                {
                    for(Int_t iW = 0; iW < Wbini; iW++)
                    {
                        Double_t A = 0., B = 0.;
                        for(Int_t iSector = 0; iSector < sectors; iSector++)
                        {
                          A += pow(y_i[iPcanv][iThPiQ][iCharge][iSector][iMod][iQ2][iW] - (y_ave[iPcanv][iThPiQ][iCharge][iMod][iQ2][iW]),2.);
                          B += var_stat[iPcanv][iThPiQ][iCharge][iSector][iMod][iQ2][iW]/CS_Sec;//# of sectors?
                        }
                        //if(iW == 25 && iPmom ==2 && iThPiQ==2 && iQ2==0 && iMod==iDATA)
                        //    {cout << iMod<<" "<< A <<"/4 - "<< B <<" \t y_ave = "<< y_ave[iPcanv][iThPiQ][iCharge][iMod][iQ2][iW] <<endl;}
                        var[iPcanv][iThPiQ][iCharge][iMod][iQ2][iW] = max(0.,(A/(CS_Sec - 1.)) - B);
                        //if(iW == 25 && iPmom ==2 && iThPiQ==2 && iQ2==0 && iMod==iDATA)
                        //    {cout <<" var = "<< var[iPcanv][iThPiQ][iCharge][iMod][iQ2][iW] <<endl;}
                    }
                }
              }
          }
          c1->Print(fileName);
          /*if(iCharge == 0 && iQ2 == 0)      c1->Print(file_CSSq1pip);
          else if(iCharge == 0 && iQ2 == 1) c1->Print(file_CSSq2pip);
          else if(iCharge == 0 && iQ2 == 2) c1->Print(file_CSSq3pip);
          else if(iCharge == 1 && iQ2 == 0) c1->Print(file_CSSq1pim);
          else if(iCharge == 1 && iQ2 == 1) c1->Print(file_CSSq2pim);
          else if(iCharge == 1 && iQ2 == 2) c1->Print(file_CSSq3pim);*/
      }
  }
  //Print page to file
  c1->Clear();

  //Last page (Cross Sections)
  //Divide page into 5x5 grid
  c1->Divide(TPIQi,Pbini);
  //-------------------------------------------------Write to file------------------------------------------//
    line = 1, size = 3, marc = 2, mark_d = 8;
    label_size = 0.09;//0.06
    for(Int_t iCharge = 0; iCharge < Nrec; iCharge++)/// RADIATIVE cross sections
    {
        for(Int_t iQ2 = 0; iQ2 < Q2bini; iQ2++)
        {
            canvas = 1;
            Double_t dCSq[2][Wbini];
            Double_t dCSt[TPIQi][2][Wbini];
            for(Int_t iW = 0; iW < Wbini; iW++)
            {
                dCSq[0][iW] = 0.; dCSq[1][iW] = 0.;
                for(Int_t iT = 0; iT < TPIQi; iT++)
                {
                    dCSt[iT][0][iW] = 0.; dCSt[iT][1][iW] = 0.;
                }
            }
            sprintf(c_NAME,"h_CSradQ%d%d",iCharge,iQ2);
            TH1D * h_CSradQ = new TH1D(c_NAME,"",Wbini,CS_Wmin,CS_Wmax);
            sprintf(c_NAME,"h_CSnoradQ%d%d",iCharge,iQ2);
            TH1D * h_CSnoradQ = new TH1D(c_NAME,"",Wbini,CS_Wmin,CS_Wmax);
            for(Int_t iPmom = 0; iPmom < Pbini; iPmom++)
            {
                Int_t iPcanv = (Pbini - 1) - iPmom;//Needed to get plot positions correct for P_pi
                Double_t dCSp[2][Wbini];
                for(Int_t iW = 0; iW < Wbini; iW++)
                {
                    dCSp[0][iW] = 0.; dCSp[1][iW] = 0.;
                }
                sprintf(c_NAME,"h_CSradP%d%d%d",iPcanv,iCharge,iQ2);
                TH1D * h_CSradP = new TH1D(c_NAME,"",Wbini,CS_Wmin,CS_Wmax);
                sprintf(c_NAME,"h_CSnoradP%d%d%d",iPcanv,iCharge,iQ2);
                TH1D * h_CSnoradP = new TH1D(c_NAME,"",Wbini,CS_Wmin,CS_Wmax);
                for(Int_t iThPiQ = 0; iThPiQ < TPIQi; iThPiQ++)
                {
                    dT = DeltaThetaPiQBin(ThetaPiQ[iThPiQ],ThetaPiQ[iThPiQ+1]);
                    Double_t dPmom = dP[iPcanv];
                    if(bQorNot==false)
                    {
                        dT = DeltaThetaPiQBin(ThetaPi[iCharge][iThPiQ],ThetaPi[iCharge][iThPiQ+1]);
                        dPmom = dPt[iCharge][iPmom+1] - dPt[iCharge][iPmom];
                    }
                    bin_vol = dPmom * dT * dQ2[iQ2] * dW;
                    c1->cd(canvas); canvas++; //Needed to get plot positions correct
                    //gPad->SetLogy();
                    sprintf(c_NAME,"h_CSopgpr%d%d%d%d",iPcanv,iThPiQ,iCharge,iQ2);
                    TH1D * h_CSopgpr = new TH1D(c_NAME,"",Wbini,CS_Wmin,CS_Wmax);
                    sprintf(c_NAME,"h_CSopgpn%d%d%d%d",iPcanv,iThPiQ,iCharge,iQ2);
                    TH1D * h_CSopgpn = new TH1D(c_NAME,"",Wbini,CS_Wmin,CS_Wmax);
                    sprintf(c_NAME,"h_CSopgmr%d%d%d%d",iPcanv,iThPiQ,iCharge,iQ2);
                    TH1D * h_CSopgmr = new TH1D(c_NAME,"",Wbini,CS_Wmin,CS_Wmax);
                    sprintf(c_NAME,"h_CSopgmn%d%d%d%d",iPcanv,iThPiQ,iCharge,iQ2);
                    TH1D * h_CSopgmn = new TH1D(c_NAME,"",Wbini,CS_Wmin,CS_Wmax);
                    ///GENIE
                    sprintf(c_NAME,"h_CSmcr%d%d%d%d",iPcanv,iThPiQ,iCharge,iQ2);
                    TH1D * h_CSmcr = new TH1D(c_NAME,"",Wbini,CS_Wmin,CS_Wmax);
                    sprintf(c_NAME,"h_CSmcn%d%d%d%d",iPcanv,iThPiQ,iCharge,iQ2);
                    TH1D * h_CSmcn = new TH1D(c_NAME,"",Wbini,CS_Wmin,CS_Wmax);
                    for(Int_t iMod = 0; iMod < iGEN; iMod++)///change iGEN to something more appropriate; also change histogram names.
                    {
                        if(bOPG==false){iMod=iGENIEr;}
                        sprintf(c_NAME,"h_CSopg%d%d%d%d%d",iPcanv,iThPiQ,iCharge,iQ2,iMod);
                        TH1D * h_CSopg = new TH1D(c_NAME,"",Wbini,CS_Wmin,CS_Wmax);
                        for(Int_t iW = 0; iW < Wbini; iW++)
                        {
                            if(iMod==iGENIEr && y_sum[iPcanv][iThPiQ][iCharge][iMod][iQ2][iW]!=0)///Load GENIE norad histograms with y_sum, rad with tCS
                            {
                                h_CSmcr->SetBinContent(iW+1, tCS[iPcanv][iThPiQ][iCharge][0][iMod][iQ2][iW]);
                                h_CSmcn->SetBinContent(iW+1, y_sum[iPcanv][iThPiQ][iCharge][iGENIE][iQ2][iW]);
                                h_CSmcr->SetBinError(iW+1, sqrt(tCS[iPcanv][iThPiQ][iCharge][0][iMod][iQ2][iW]));
                                h_CSmcn->SetBinError(iW+1, sqrt(y_sum[iPcanv][iThPiQ][iCharge][iGENIE][iQ2][iW]));
                                ///Load Rad
                                dCSq[0][iW] += tCS[iPcanv][iThPiQ][iCharge][0][iMod][iQ2][iW];
                                dCSp[0][iW] += tCS[iPcanv][iThPiQ][iCharge][0][iMod][iQ2][iW];
                                dCSt[iThPiQ][0][iW] += tCS[iPcanv][iThPiQ][iCharge][0][iMod][iQ2][iW];
                                ///Load NoRad
                                dCSq[1][iW] += y_sum[iPcanv][iThPiQ][iCharge][iGENIE][iQ2][iW];
                                dCSp[1][iW] += y_sum[iPcanv][iThPiQ][iCharge][iGENIE][iQ2][iW];
                                dCSt[iThPiQ][1][iW] += y_sum[iPcanv][iThPiQ][iCharge][iGENIE][iQ2][iW];
                                if(iCharge==0 && iQ2==0 && iPcanv==2 && iThPiQ==0 && iW==35){
                                    cout << "# of evts: "<< y_sum[iPcanv][iThPiQ][iCharge][iMod][iQ2][iW] << "\t";
                                    cout << "w of evts: "<< tCS[iPcanv][iThPiQ][iCharge][0][iMod][iQ2][iW] << "\t";
                                }
                            }
                            else if(y_sum[iPcanv][iThPiQ][iCharge][iMod][iQ2][iW]!=0)///Scale h_CSopg bins by total CS for onepigen
                            {
                                Double_t opg_error = tCS[iPcanv][iThPiQ][iCharge][0][iMod][iQ2][iW]/sqrt(y_sum[iPcanv][iThPiQ][iCharge][iMod][iQ2][iW]);
                                h_CSopg->SetBinContent(iW+1, tCS[iPcanv][iThPiQ][iCharge][0][iMod][iQ2][iW]);
                                h_CSopg->SetBinError(iW+1, opg_error);
                            }
                        }   // Note: for an array bin # of iW, the ROOT bin # is iW+1
                        if(iMod==iOPGr && iCharge==0)//onepigen pi+
                        {
                            h_CSopg->Scale(ub2nb / (opg_gen_events[0] * bin_vol));//kBlue
                            line = 1; mark = 26; size = 4; marc = 4; h_CSopg->SetLineColor(4);
                            h_CSopgpr = h_CSopg;
                        }
                        else if(iMod==iOPGr && iCharge==1)//onepigen pi-
                        {
                            h_CSopg->Scale(ub2nb / (opg_gen_events[1] * bin_vol));/// <-------------------<<
                            line = 1; mark = 26; size = 4; marc = 4; h_CSopg->SetLineColor(4);//38
                            h_CSopgmr = h_CSopg;
                        }
                        else if(iMod==iOPGn && iCharge==0)//onepigen pi+ norad
                        {
                            h_CSopg->Scale(ub2nb / (opg_gen_events[2] * bin_vol));//kBlue
                            line = 1; mark = 26; size = 4; marc = 3; h_CSopg->SetLineColor(4);
                            h_CSopgpn = h_CSopg;
                        }
                        else if(iMod==iOPGn && iCharge==1)//onepigen pi- norad
                        {
                            h_CSopg->Scale(ub2nb / (opg_gen_events[3] * bin_vol));/// <-------------------<<
                            line = 1; mark = 26; size = 4; marc = 3; h_CSopg->SetLineColor(38);
                            h_CSopgmn = h_CSopg;
                        }
                        else if(iMod==iGENIEr)//GENIE rad
                        {
                            h_CSmcn->Scale(totalCS_MC[iGENmod] * ub2nb / (MC_evts[iGENmod] * bin_vol));/// <-------------------<<
                            h_CSmcr->Scale(totalCS_MC[iGENmod] * ub2nb / (MC_rad_evts[iGENmod] * bin_vol));/// <-------------------<<
                            hr_ARC[iPcanv][iThPiQ][iCharge][iQ2]->Divide(h_CSmcr);
                            //line = 1; mark = 26; size = 4; marc = 3; h_CSopg->SetLineColor(38);
                        }
                        /*else if(iMod==iGENIE)//GENIE norad
                        {
                            h_CSopg->Scale(totalCS_MC[iGENmod] * ub2nb / (MC_evts[iGENmod] * bin_vol));/// <-------------------<<
                            line = 1; mark = 26; size = 4; marc = 3; h_CSopg->SetLineColor(38);
                            h_CSmcn = h_CSopg;
                        }*/
                        gStyle->SetEndErrorSize(10);
                        gStyle->SetErrorX(0);
                    }/// Working on the below:
                    hr_CS[iPcanv][iThPiQ][iCharge][iQ2]->Sumw2(); hr_CSmc[iPcanv][iThPiQ][iCharge][iQ2]->Sumw2();
                    if(iCharge == 0)
                    {
                         hr_CS[iPcanv][iThPiQ][iCharge][iQ2]->Divide(h_CSopgpr,h_CSopgpn);
                    }
                    else if(iCharge == 1)
                    {
                         hr_CS[iPcanv][iThPiQ][iCharge][iQ2]->Divide(h_CSopgmr,h_CSopgmn);
                    }
                    ///GENIE RC
                    hr_CSmc[iPcanv][iThPiQ][iCharge][iQ2]->Divide(h_CSmcr,h_CSmcn);
                    for(Int_t iW = 0; iW < Wbini; iW++)
                    {
                        if(y_sum[iPcanv][iThPiQ][iCharge][iGENIEr][iQ2][iW]!=0)///Calculate GENIE RAD Correction Error
                        {
                            Double_t RadContent = h_CSmcr->GetBinContent(iW+1);
                            Double_t NoRadContent = h_CSmcn->GetBinContent(iW+1);
                            Double_t RadError = h_CSmcr->GetBinError(iW+1);
                            Double_t NoRadError = h_CSmcn->GetBinError(iW+1);
                            Double_t RatioError = RadContent/NoRadContent*sqrt(pow(RadError/RadContent,2)+pow(NoRadError/NoRadContent,2));
                            hr_CSmc[iPcanv][iThPiQ][iCharge][iQ2]->SetBinError(iW+1, RatioError);
                        }
                    }   // Note: for an array bin # of iW, the ROOT bin # is iW+1
                }
                ///Calculate 3D P_pi RC
                hr_RCp[iPcanv][iCharge][iQ2]->Sumw2();
                for(Int_t iW = 0; iW < Wbini; iW++)
                {
                    h_CSradP->SetBinContent(iW+1, dCSp[0][iW]);// * totalCS_MC[iGENmod] * ub2nb / (MC_evts[iGENmod]));
                    h_CSnoradP->SetBinContent(iW+1, dCSp[1][iW]);// * totalCS_MC[iGENmod] * ub2nb / (MC_evts[iGENmod]));
                    if(dCSp[1][iW] > 0. && dCSp[1][iW] < 100. && dCSp[0][iW] > 0. && dCSp[0][iW] < 100.){
                        Double_t ratio_radp = h_CSradP->GetBinContent(iW+1)/h_CSnoradP->GetBinContent(iW+1);
                        hr_RCp[iPcanv][iCharge][iQ2]->SetBinContent(iW+1, h_CSradP->GetBinContent(iW+1)/h_CSnoradP->GetBinContent(iW+1));
                        Double_t error_p = ratio_radp * sqrt(pow(sqrt(dCSp[0][iW])/dCSp[0][iW],2)+pow(sqrt(dCSp[1][iW])/dCSp[1][iW],2));
                        hr_RCp[iPcanv][iCharge][iQ2]->SetBinError(iW+1, error_p);
                        //cout <<"P: "<<iW<<"\t"<<hr_RCp[iPcanv][iCharge][iQ2]->GetBinContent(iW+1)<<endl;
                    }
                }   // Note: for an array bin # of iW, the ROOT bin # is iW+1
                hr_ARCp[iPcanv][iCharge][iQ2]->Divide(h_CSradP);
            }
            ///Calculate 3D Theta_piq RC
            for(Int_t iThPiQ = 0; iThPiQ < TPIQi; iThPiQ++)
            {
                hr_RCt[iThPiQ][iCharge][iQ2]->Sumw2();
                sprintf(c_NAME,"h_CSradT%d%d%d",iThPiQ,iCharge,iQ2);
                TH1D * h_CSradT = new TH1D(c_NAME,"",Wbini,CS_Wmin,CS_Wmax);
                sprintf(c_NAME,"h_CSnoradT%d%d%d",iThPiQ,iCharge,iQ2);
                TH1D * h_CSnoradT = new TH1D(c_NAME,"",Wbini,CS_Wmin,CS_Wmax);
                for(Int_t iW = 0; iW < Wbini; iW++)
                {
                    h_CSradT->SetBinContent(iW+1, dCSt[iThPiQ][0][iW]);// * totalCS_MC[iGENmod] * ub2nb / (MC_evts[iGENmod]));
                    h_CSnoradT->SetBinContent(iW+1, dCSt[iThPiQ][1][iW]);// * totalCS_MC[iGENmod] * ub2nb / (MC_evts[iGENmod]));
                    if(dCSt[iThPiQ][1][iW] > 0. && dCSt[iThPiQ][1][iW] < 100. && dCSt[iThPiQ][0][iW] > 0. && dCSt[iThPiQ][0][iW] < 100.){
                        Double_t ratio_radt = h_CSradT->GetBinContent(iW+1)/h_CSnoradT->GetBinContent(iW+1);
                        hr_RCt[iThPiQ][iCharge][iQ2]->SetBinContent(iW+1, h_CSradT->GetBinContent(iW+1)/h_CSnoradT->GetBinContent(iW+1));
                        Double_t error_t = ratio_radt * sqrt(pow(sqrt(dCSt[iThPiQ][0][iW])/dCSt[iThPiQ][0][iW],2)+pow(sqrt(dCSt[iThPiQ][1][iW])/dCSt[iThPiQ][1][iW],2));
                        hr_RCt[iThPiQ][iCharge][iQ2]->SetBinError(iW+1, error_t);
                        //cout <<"T: "<<iW<<"\t"<<hr_RCt[iThPiQ][iCharge][iQ2]->GetBinContent(iW+1)<<endl;
                    }
                }   // Note: for an array bin # of iW, the ROOT bin # is iW+1
                hr_ARCt[iThPiQ][iCharge][iQ2]->Divide(h_CSradT);
            }
            ///Calculate 2D Q2 RC
            hr_RCq[iCharge][iQ2]->Sumw2();
            for(Int_t iW = 0; iW < Wbini; iW++)
            {
                h_CSradQ->SetBinContent(iW+1, dCSq[0][iW]);// * totalCS_MC[iGENmod] * ub2nb / (MC_evts[iGENmod]));
                h_CSnoradQ->SetBinContent(iW+1, dCSq[1][iW]);// * totalCS_MC[iGENmod] * ub2nb / (MC_evts[iGENmod]));
                if(dCSq[1][iW] > 0. && dCSq[1][iW] < 100. && dCSq[0][iW] > 0. && dCSq[0][iW] < 100.){
                    hr_RCq[iCharge][iQ2]->SetBinContent(iW+1, h_CSradQ->GetBinContent(iW+1)/h_CSnoradQ->GetBinContent(iW+1));
                    Double_t ratio_radq = h_CSradQ->GetBinContent(iW+1)/h_CSnoradQ->GetBinContent(iW+1);
                    Double_t error_q = ratio_radq * sqrt(pow(sqrt(dCSq[0][iW])/dCSq[0][iW],2)+pow(sqrt(dCSq[1][iW])/dCSq[1][iW],2));
                    hr_RCq[iCharge][iQ2]->SetBinError(iW+1, error_q);
                    //cout <<"Q: "<<iW<<"\t"<<hr_RCq[iCharge][iQ2]->GetBinContent(iW+1)<<endl;
                }
            }   // Note: for an array bin # of iW, the ROOT bin # is iW+1
            hr_ARCq[iCharge][iQ2]->Divide(h_CSradQ);
        }
    }///End of RC loops
    for(Int_t iCharge = 0; iCharge < Nrec; iCharge++)///Start of onepigen RC calculator loops
    {if(bOPG==false){continue;}
        for(Int_t iQ2 = 0; iQ2 < Q2bini; iQ2++)
        {
            canvas = 1;
            for(Int_t iPmom = 0; iPmom < Pbini; iPmom++)
            {
                Int_t iPcanv = (Pbini - 1) - iPmom;//Needed to get plot positions correct for P_pi
                for(Int_t iThPiQ = 0; iThPiQ < TPIQi; iThPiQ++)
                {
                    c1->cd(canvas); canvas++; //Needed to get plot positions correct
                    gPad->SetLogy(0);
                    hr_CS[iPcanv][iThPiQ][iCharge][iQ2]->Draw("e1p");
                    //hr_CS[iPcanv][iThPiQ][iCharge][iQ2]->Write("e1p");
                    hr_CS[iPcanv][iThPiQ][iCharge][iQ2]->GetYaxis()->SetNdivisions(4);
                    hr_CS[iPcanv][iThPiQ][iCharge][iQ2]->GetYaxis()->SetLabelSize(0.09);
                    hr_CS[iPcanv][iThPiQ][iCharge][iQ2]->SetMaximum(2.0);
                    hr_CS[iPcanv][iThPiQ][iCharge][iQ2]->SetMinimum(0.0);
                    hr_CS[iPcanv][iThPiQ][iCharge][iQ2]->GetXaxis()->SetNdivisions(4);
                    hr_CS[iPcanv][iThPiQ][iCharge][iQ2]->GetXaxis()->SetLabelSize(0.09);
                    if(iThPiQ == 2 && iPcanv == Pbini - 1 && iQ2==0)
                    {
                        if(iCharge == 0){hr_CS[iPcanv][iThPiQ][iCharge][iQ2]->SetTitle("pi+   0.70 < Q2 < 1.0;;");}
                        else if(iCharge == 1){hr_CS[iPcanv][iThPiQ][iCharge][iQ2]->SetTitle("pi-   0.70 < Q2 < 1.0;;");}
                    }
                    else if(iThPiQ == 2 && iPcanv == Pbini - 1 && iQ2==1)
                    {
                        if(iCharge == 0){hr_CS[iPcanv][iThPiQ][iCharge][iQ2]->SetTitle("pi+   1.0 < Q2 < 1.4;;");}
                        else if(iCharge == 1){hr_CS[iPcanv][iThPiQ][iCharge][iQ2]->SetTitle("pi-   1.0 < Q2 < 1.4;;");}
                    }
                    else if(iThPiQ == 2 && iPcanv == Pbini - 1 && iQ2==2)
                    {
                        if(iCharge == 0){hr_CS[iPcanv][iThPiQ][iCharge][iQ2]->SetTitle("pi+   1.4 < Q2 < 1.9;;");}
                        else if(iCharge == 1){hr_CS[iPcanv][iThPiQ][iCharge][iQ2]->SetTitle("pi-   1.4 < Q2 < 1.9;;");}
                    }
                    else if(iThPiQ == 2 && iPcanv == Pbini - 1 && iQ2==3)
                    {
                        if(iCharge == 0){hr_CS[iPcanv][iThPiQ][iCharge][iQ2]->SetTitle("pi+   1.9 < Q2 < 2.5;;");}
                        else if(iCharge == 1){hr_CS[iPcanv][iThPiQ][iCharge][iQ2]->SetTitle("pi-   1.9 < Q2 < 2.5;;");}
                    }
                    TF1 *RadCorr = new TF1("RadCorr","[0]",WbinEdge[iQ2][0],WbinEdge[iQ2][1]);
                    //RadCorr->SetParameter(0,0.5);
                    hr_CS[iPcanv][iThPiQ][iCharge][iQ2]->Fit("RadCorr");
                    RadFitAve[iCharge] += RadCorr->GetParameter(0)*RadCorr->GetChisquare();
                    RadFitErr[iCharge] += RadCorr->GetChisquare();
                    if(iCharge==0)
                        h2_RC_chi2pip->Fill(RadCorr->GetChisquare()/RadCorr->GetNDF(),RadCorr->GetParameter(0));
                    else if(iCharge==1)
                        h2_RC_chi2pim->Fill(RadCorr->GetChisquare()/RadCorr->GetNDF(),RadCorr->GetParameter(0));
                }
            }
            if(iCharge == 0 && iQ2 == 0)      c1->Print(file_CSradq1pip);
            else if(iCharge == 0 && iQ2 == 1) c1->Print(file_CSradq2pip);
            else if(iCharge == 0 && iQ2 == 2) c1->Print(file_CSradq3pip);
            else if(iCharge == 0 && iQ2 == 3) c1->Print(file_CSradq4pip);
            else if(iCharge == 1 && iQ2 == 0) c1->Print(file_CSradq1pim);
            else if(iCharge == 1 && iQ2 == 1) c1->Print(file_CSradq2pim);
            else if(iCharge == 1 && iQ2 == 2) c1->Print(file_CSradq3pim);
            else if(iCharge == 1 && iQ2 == 3) c1->Print(file_CSradq4pim);
        }
        //RadFitAve[iCharge] /= RadFitBins[iCharge];// divide by number of plots (# of Q2 bins + # of theta_pi(q) bins + # of P_pi bins)
        RadFitAve[iCharge] /= RadFitErr[iCharge];// divide by sum of weights (to get weighted average [or mean])
        RadFitUnc[iCharge] = abs((1./RadFitAve[iCharge]) - 1.) * 0.2;
        for(Int_t iQ2 = 0; iQ2 < Q2bini; iQ2++)
        {
            canvas = 1;
            for(Int_t iPmom = 0; iPmom < Pbini; iPmom++)
            {
                Int_t iPcanv = (Pbini - 1) - iPmom;//Needed to get plot positions correct for P_pi
                for(Int_t iThPiQ = 0; iThPiQ < TPIQi; iThPiQ++)
                {
                    c1->cd(canvas); canvas++; //Needed to get plot positions correct
                    gPad->SetLogy(0);
                    hr_CS[iPcanv][iThPiQ][iCharge][iQ2]->Draw("e1p");
                    //hr_CS[iPcanv][iThPiQ][iCharge][iQ2]->Write("e1p");
                    hr_CS[iPcanv][iThPiQ][iCharge][iQ2]->GetYaxis()->SetNdivisions(4);
                    hr_CS[iPcanv][iThPiQ][iCharge][iQ2]->GetYaxis()->SetLabelSize(0.09);
                    hr_CS[iPcanv][iThPiQ][iCharge][iQ2]->SetMaximum(2.0);
                    hr_CS[iPcanv][iThPiQ][iCharge][iQ2]->SetMinimum(0.0);
                    hr_CS[iPcanv][iThPiQ][iCharge][iQ2]->GetXaxis()->SetNdivisions(4);
                    hr_CS[iPcanv][iThPiQ][iCharge][iQ2]->GetXaxis()->SetLabelSize(0.09);
                    if(iThPiQ == 2 && iPcanv == Pbini - 1 && iQ2==0)
                    {
                        if(iCharge == 0){hr_CS[iPcanv][iThPiQ][iCharge][iQ2]->SetTitle("pi+   0.70 < Q2 < 1.0;;");}
                        else if(iCharge == 1){hr_CS[iPcanv][iThPiQ][iCharge][iQ2]->SetTitle("pi-   0.70 < Q2 < 1.0;;");}
                    }
                    else if(iThPiQ == 2 && iPcanv == Pbini - 1 && iQ2==1)
                    {
                        if(iCharge == 0){hr_CS[iPcanv][iThPiQ][iCharge][iQ2]->SetTitle("pi+   1.0 < Q2 < 1.5;;");}
                        else if(iCharge == 1){hr_CS[iPcanv][iThPiQ][iCharge][iQ2]->SetTitle("pi-   1.0 < Q2 < 1.5;;");}
                    }
                    else if(iThPiQ == 2 && iPcanv == Pbini - 1 && iQ2==2)
                    {
                        if(iCharge == 0){hr_CS[iPcanv][iThPiQ][iCharge][iQ2]->SetTitle("pi+   1.5 < Q2 < 2.5;;");}
                        else if(iCharge == 1){hr_CS[iPcanv][iThPiQ][iCharge][iQ2]->SetTitle("pi-   1.5 < Q2 < 2.5;;");}
                    }
                    TLine *lRadCorrAve = new TLine(WbinEdge[iQ2][0],RadFitAve[iCharge],WbinEdge[iQ2][1],RadFitAve[iCharge]);
                    lRadCorrAve->Draw();
                }
            }
            c1->Print(fileName);
        }
    }
  //c1->Clear();
  //c1->Divide(TPIQi,Pbini);
    for(Int_t iCharge = 0; iCharge < Nrec; iCharge++)///Start of GENIE RC printing loops
    {
        for(Int_t iQ2 = 0; iQ2 < Q2bini; iQ2++)
        {
            canvas = 1;
            for(Int_t iPmom = 0; iPmom < Pbini; iPmom++)
            {
                Int_t iPcanv = (Pbini - 1) - iPmom;//Needed to get plot positions correct for P_pi
                for(Int_t iThPiQ = 0; iThPiQ < TPIQi; iThPiQ++)
                {
                    c1->cd(canvas); canvas++; //Needed to get plot positions correct
                    gPad->SetLogy(0);
                    hr_ARC[iPcanv][iThPiQ][iCharge][iQ2]->Draw("e1p");
                    //hr_CSmc[iPcanv][iThPiQ][iCharge][iQ2]->Write("e1p");
                    hr_ARC[iPcanv][iThPiQ][iCharge][iQ2]->GetYaxis()->SetNdivisions(4);
                    hr_ARC[iPcanv][iThPiQ][iCharge][iQ2]->GetYaxis()->SetLabelSize(0.09);
                    hr_ARC[iPcanv][iThPiQ][iCharge][iQ2]->SetMaximum(2.0);
                    hr_ARC[iPcanv][iThPiQ][iCharge][iQ2]->SetMinimum(0.0);
                    hr_ARC[iPcanv][iThPiQ][iCharge][iQ2]->GetXaxis()->SetNdivisions(4);
                    hr_ARC[iPcanv][iThPiQ][iCharge][iQ2]->GetXaxis()->SetLabelSize(0.09);
                    if(iThPiQ == 2 && iPcanv == Pbini - 1 && iQ2==0)
                    {
                        if(iCharge == 0){hr_ARC[iPcanv][iThPiQ][iCharge][iQ2]->SetTitle("GENIE pi+   0.70 < Q2 < 1.0;;");}
                        else if(iCharge == 1){hr_ARC[iPcanv][iThPiQ][iCharge][iQ2]->SetTitle("GENIE pi-   0.70 < Q2 < 1.0;;");}
                    }
                    else if(iThPiQ == 2 && iPcanv == Pbini - 1 && iQ2==1)
                    {
                        if(iCharge == 0){hr_ARC[iPcanv][iThPiQ][iCharge][iQ2]->SetTitle("GENIE pi+   1.0 < Q2 < 1.4;;");}
                        else if(iCharge == 1){hr_ARC[iPcanv][iThPiQ][iCharge][iQ2]->SetTitle("GENIE pi-   1.0 < Q2 < 1.4;;");}
                    }
                    else if(iThPiQ == 2 && iPcanv == Pbini - 1 && iQ2==2)
                    {
                        if(iCharge == 0){hr_ARC[iPcanv][iThPiQ][iCharge][iQ2]->SetTitle("GENIE pi+   1.4 < Q2 < 1.9;;");}
                        else if(iCharge == 1){hr_ARC[iPcanv][iThPiQ][iCharge][iQ2]->SetTitle("GENIE pi-   1.4 < Q2 < 1.9;;");}
                    }
                    else if(iThPiQ == 2 && iPcanv == Pbini - 1 && iQ2==3)
                    {
                        if(iCharge == 0){hr_ARC[iPcanv][iThPiQ][iCharge][iQ2]->SetTitle("GENIE pi+   1.9 < Q2 < 2.5;;");}
                        else if(iCharge == 1){hr_ARC[iPcanv][iThPiQ][iCharge][iQ2]->SetTitle("GENIE pi-   1.9 < Q2 < 2.5;;");}
                    }
                    TLine *lRadCorrAve = new TLine(WbinEdge[iQ2][0],RadFitAve[iCharge],WbinEdge[iQ2][1],RadFitAve[iCharge]);
                    lRadCorrAve->Draw();
                }
            }
            /*if(iCharge == 0 && iQ2 == 0)      c1->Print(file_CSradq1pip);
            else if(iCharge == 0 && iQ2 == 1) c1->Print(file_CSradq2pip);
            else if(iCharge == 0 && iQ2 == 2) c1->Print(file_CSradq3pip);
            else if(iCharge == 0 && iQ2 == 3) c1->Print(file_CSradq4pip);
            else if(iCharge == 1 && iQ2 == 0) c1->Print(file_CSradq1pim);
            else if(iCharge == 1 && iQ2 == 1) c1->Print(file_CSradq2pim);
            else if(iCharge == 1 && iQ2 == 2) c1->Print(file_CSradq3pim);
            else if(iCharge == 1 && iQ2 == 3) c1->Print(file_CSradq4pim);*/
            c1->Print(fileName);
        }
    }
  c1->Clear();
  c1->Divide(TPIQi,Pbini);
  //-------------------------------------------------Write to file------------------------------------------//
  ///CStxtFile << "*******************************************************************************************" << endl;
  ///if(bQorNot==true){CStxtFile << "4D binning: Q2, P_pi and Theta_piq" << endl;}
  ///else             {CStxtFile << "4D binning: Q2, P_pi and Theta_pi" << endl;}
  ///CStxtFile << "*******************************************************************************************" << endl;
  ///CStxtFile << "iMod\t iChrg\t iQ2\t iPcanv\t iTheta\t W\t CS\t error" << endl;
  ifstream unc_pid("delta_a2_unc_pid.txt");
  ifstream unc_Vze("delta_a2_unc_vze.txt");
  ifstream unc_Vre("delta_a2_unc_vpe.txt");
  ifstream unc_Vrpi("delta_a2_unc_vppi.txt");
  ifstream unc_Vdiff("delta_a2_unc_vd.txt");
  unsigned int iunc_bin = 0;
    line = 1, size = 3, marc = 2, mark_d = 8;
    label_size = 0.11;//0.06 0.09
    Double_t xlab_size = 0.11;//0.06 0.09
    for(Int_t iCharge = 0; iCharge < Nrec; iCharge++)/// cross sections
    {
        //Double_t b2b_error = sqrt(pow(u_eID,2)+pow(u_hID,2)+pow(u_Vze,2)+pow(u_Vre,2)+pow(u_Vrpi,2)+pow(u_Vdiff,2)+pow(u_TW,2)+pow(RadFitUnc[iCharge],2));
        for(Int_t iQ2 = 0; iQ2 < Q2bini; iQ2++)
        {
            canvas = 1;
            for(Int_t iPmom = 0; iPmom < Pbini; iPmom++)
            {
                Int_t iPcanv = (Pbini - 1) - iPmom;//Needed to get plot positions correct for P_pi
                for(Int_t iThPiQ = 0; iThPiQ < TPIQi; iThPiQ++)
                {
                    dT = DeltaThetaPiQBin(ThetaPiQ[iThPiQ],ThetaPiQ[iThPiQ+1]);
                    Double_t dPmom = dP[iPcanv];
                    Double_t CSvalue[models][Wbini];
                    Double_t CSerror[models][Wbini];
                    if(bQorNot==false)
                    {
                        dT = DeltaThetaPiQBin(ThetaPi[iCharge][iThPiQ],ThetaPi[iCharge][iThPiQ+1]);
                        dPmom = dPt[iCharge][iPmom+1] - dPt[iCharge][iPmom];
                    }
                    bin_vol = dPmom * dT * dQ2[iQ2] * dW;
                    c1->cd(canvas); canvas++; //Needed to get plot positions correct
                    sprintf(c_NAME,"hs_CS%d%d%d%d",iPcanv,iThPiQ,iCharge,iQ2);
                    THStack * hs_CS = new THStack(c_NAME,"");
                    //gPad->SetLogy();
                    for(Int_t iMod = 0; iMod < models; iMod++)
                    {
                        //if(iMod==iGENIE) continue;
                        sprintf(c_NAME,"h_CS%d%d%d%d%d",iPcanv,iThPiQ,iCharge,iQ2,iMod);
                        TH1D * h_CS = new TH1D(c_NAME,"",Wbini,CS_Wmin,CS_Wmax);
                        h_CS->Sumw2(); //https://root-forum.cern.ch/t/array-of-histograms/13148/2
                        for(Int_t iW = 0; iW < Wbini; iW++)
                        {
                            CSvalue[iMod][iW] = 0.;
                            CSerror[iMod][iW] = 0.;
                            /// Acceptance and radiative corrections
                            Double_t ACorr = hr_AC[iPcanv][iThPiQ][iCharge][iQ2]->GetBinContent(iW+1);
                            Double_t RCorr = hr_CSmc[iPcanv][iThPiQ][iCharge][iQ2]->GetBinContent(iW+1);
                            if(bAC == false){ACorr=1.;}
                            if(bRAD == false){RCorr=1.;}
                            if(iMod<iGENIEr && y_sum[iPcanv][iThPiQ][iCharge][iMod][iQ2][iW]!=0)///Scale h_CS bins by total CS for onepigen
                            {
                                Double_t opg_error = tCS[iPcanv][iThPiQ][iCharge][0][iMod][iQ2][iW]/sqrt(y_sum[iPcanv][iThPiQ][iCharge][iMod][iQ2][iW]);
                                h_CS->SetBinContent(iW+1, ACorr * tCS[iPcanv][iThPiQ][iCharge][0][iMod][iQ2][iW]);
                                h_CS->SetBinError(iW+1, opg_error);
                                //h_CS->Scale(hr_AC[iPcanv][iThPiQ][iCharge][iQ2]->GetBinContent(iW+1));/// Acceptance corrections
                            }
                            else if(iMod>=iGENIEr && RCorr>0. && RCorr<10.)///Do not scale h_CS bins by total CS for non-onepigen (or bins have zero events)
                            {
                                Double_t y_val = y_sum[iPcanv][iThPiQ][iCharge][iMod][iQ2][iW];
                                /// Set cross section bin values and errors
                                h_CS->SetBinContent(iW+1, y_val * ACorr/RCorr);
                                h_CS->SetBinError(iW+1, sqrt(y_val));
                                //if(RCorr>0.){
                                //    h_CS->Scale(ACorr/RCorr);
                                //}
                                //else{h_CS->Scale(ACorr);}
                                if(iMod==iDATA && txtFile.is_open())///For cut analysis
                                {
                                    Double_t BinEntries = y_sum[iPcanv][iThPiQ][iCharge][iMod][iQ2][iW] *cm2ub*ub2nb/(L*bin_vol*RadFitAve[iCharge]);
                                    BinEntries *= ACorr/RCorr;//hr_AC[iPcanv][iThPiQ][iCharge][iQ2]->GetBinContent(iW+1);
                                    txtFile << BinEntries << endl;
                                }
                            }
                        }   // Note: for an array bin # of iW, the ROOT bin # is iW+1
                        //if (h_CS->GetSumw2N() == 0) h_CS->Sumw2(kTRUE);
                        if(iMod==iDATA)//Data
                        {
                            h_CS->Scale(cm2ub * ub2nb / (L * bin_vol));//kBlack
                            //h_CS->Scale(1. / RadFitAve[iCharge]);/// Radiative corrections
                            line = 1; mark = 8;  size = 5; marc = 1; h_CS->SetLineColor(1);
                        }
                        else if(iMod==iGENIE)//GENIE
                        {
                            h_CS->Scale(totalCS_MC[iGENmod] * ub2nb / (MC_evts[iGENmod] * bin_vol));//kRed
                            //h_CS->Scale(pntF * ub2nb / (MC_evts[iGENmod] * bin_vol));
                            line = 1; mark = 32; size = 5; marc = 2; h_CS->SetLineColor(2);
                        }
                        if(iMod==iGENIEr)//Data
                        {
                            h_CS->Scale(totalCS_MC[iGENmod] * ub2nb / (MC_rad_evts[iGENmod] * bin_vol));
                            line = 1; mark = 32;  size = 5; marc = 6; h_CS->SetLineColor(6);
                        }
                        else if(iMod==iOPGr && iCharge==0)//onepigen pi+
                        {
                            h_CS->Scale(ub2nb / (opg_gen_events[0] * bin_vol));//kBlue
                            //h_CS->Scale(1. / RadFitAve[iCharge]);/// Radiative corrections
                            line = 1; mark = 26; size = 4; marc = 4; h_CS->SetLineColor(4);
                        }
                        else if(iMod==iOPGr && iCharge==1)//onepigen pi-
                        {
                            h_CS->Scale(ub2nb / (opg_gen_events[1] * bin_vol));/// <-------------------<<
                            //h_CS->Scale(1. / RadFitAve[iCharge]);/// Radiative corrections
                            line = 1; mark = 26; size = 4; marc = 4; h_CS->SetLineColor(4);//38
                        }
                        else if(iMod==iOPGn && iCharge==0)//onepigen pi+ norad
                        {
                            h_CS->Scale(ub2nb / (opg_gen_events[2] * bin_vol));//kBlue
                            line = 1; mark = 26; size = 4; marc = 4; h_CS->SetLineColor(4);//marc=3
                        }
                        else if(iMod==iOPGn && iCharge==1)//onepigen pi- norad
                        {
                            h_CS->Scale(ub2nb / (opg_gen_events[3] * bin_vol));/// <-------------------<<
                            line = 1; mark = 26; size = 4; marc = 4; h_CS->SetLineColor(4);//38);
                        }
                        h_CS->SetLineStyle(line);
                        h_CS->SetMarkerStyle(mark);
                        h_CS->SetMarkerSize(2);//0.5 or size
                        h_CS->SetMarkerColor(marc);
                        gStyle->SetEndErrorSize(10);
                        gStyle->SetErrorX(0);
                        if(iMod==iDATA)
                        {
                            for(Int_t iW = 0; iW < Wbini; iW++)
                            {   ///Retrieving point-to-point cut uncertainties
                                Double_t dunc_pid = 0., dunc_Vze = 0., dunc_Vre = 0., dunc_Vrpi = 0., dunc_Vdiff = 0.;
                                if(unc_pid.is_open())///For cut uncertainty
                                {
                                    GotoLine(unc_pid, iunc_bin);
                                    unc_pid >> dunc_pid; //if(iunc_bin<30){cout << dunc_pid << "\t";}
                                }
                                if(unc_Vze.is_open())///For cut uncertainty
                                {
                                    GotoLine(unc_Vze, iunc_bin);
                                    unc_Vze >> dunc_Vze; //if(iunc_bin<30){cout << dunc_Vze << endl;}
                                }
                                if(unc_Vre.is_open())///For cut uncertainty
                                {
                                    GotoLine(unc_Vre, iunc_bin);
                                    unc_Vre >> dunc_Vre;
                                }
                                if(unc_Vrpi.is_open())///For cut uncertainty
                                {
                                    GotoLine(unc_Vrpi, iunc_bin);
                                    unc_Vrpi >> dunc_Vrpi;
                                }
                                if(unc_Vdiff.is_open())///For cut uncertainty
                                {
                                    GotoLine(unc_Vdiff, iunc_bin);
                                    unc_Vdiff >> dunc_Vdiff;
                                }
                                iunc_bin++;
                                Double_t error = pow(h_CS->GetBinError(iW+1),2);                    //statistical [originally h_CSerr]
                                error += var[iPcanv][iThPiQ][iCharge][iMod][iQ2][iW];               //sector-by-sector uncertainty
                                //error += pow(b2b_error * h_CS->GetBinContent(iW+1),2); //bin-to-bin correction uncertainty
                                //error += pow(0.1 * cm2ub * ub2nb / (L * bin_vol),2);                //normalization uncertainty
                                //Double_t radunc = RadFitUnc[iCharge];///radiative correction uncertainty (onepigen)
                                Double_t radunc = hr_CSmc[iPcanv][iThPiQ][iCharge][iQ2]->GetBinError(iW+1);///radiative correction uncertainty (GENIE)
                                error += pow(radunc * h_CS->GetBinContent(iW+1),2); ///radiative correction uncertainty (GENIE)
                                error += pow(u_eID * h_CS->GetBinContent(iW+1),2); //eID uncertainty
                                error += pow(dunc_pid * h_CS->GetBinContent(iW+1),2); //piID uncertainty
                                error += pow(dunc_Vze * h_CS->GetBinContent(iW+1),2); //e- Vertex Z uncertainty
                                error += pow(dunc_Vre * h_CS->GetBinContent(iW+1),2); //e- perpendicular vertex uncertainty
                                error += pow(dunc_Vrpi * h_CS->GetBinContent(iW+1),2); //pi perpendicular vertex uncertainty
                                error += pow(dunc_Vdiff * h_CS->GetBinContent(iW+1),2); //vertex z difference uncertainty
                                error += pow(u_TW * h_CS->GetBinContent(iW+1),2); //target uncertainty
                                if(iMod==iDATA && iQ2==0 && iCharge==0 && iThPiQ<2 && error>0.){
                                cout <<endl<<"***********************************************\n";
                                cout <<endl<<"P: "<<iPcanv<<"   T: "<<iThPiQ<<"   W: "<<iW;
                                cout <<endl<<"Statistical: "<< pow(h_CS->GetBinError(iW+1),2);
                                cout <<endl<<"Sec-by-Sec:  "<< var[iPcanv][iThPiQ][iCharge][iMod][iQ2][iW];
                                cout <<endl<<"Radiative:   "<< pow(radunc * h_CS->GetBinContent(iW+1),2);
                                cout <<endl<<"Electron ID: "<< pow(u_eID * h_CS->GetBinContent(iW+1),2); //eID uncertainty
                                cout <<endl<<"Pion ID:     "<< pow(dunc_pid * h_CS->GetBinContent(iW+1),2); //piID uncertainty
                                cout <<endl<<"e- Vertex Z: "<< pow(dunc_Vze * h_CS->GetBinContent(iW+1),2); //e- Vertex Z uncertainty
                                cout <<endl<<"e- V_perp:   "<< pow(dunc_Vre * h_CS->GetBinContent(iW+1),2); //e- perpendicular vertex uncertainty
                                cout <<endl<<"pi V_perp:   "<< pow(dunc_Vrpi * h_CS->GetBinContent(iW+1),2); //pi perpendicular vertex uncertainty
                                cout <<endl<<"Vdiff:       "<< pow(dunc_Vdiff * h_CS->GetBinContent(iW+1),2); //vertex z difference uncertainty
                                cout <<endl<<"Target:      "<< pow(u_TW * h_CS->GetBinContent(iW+1),2); //target uncertainty
                                cout <<endl<<"Total Uncer: "<< error; //total uncertainty
                                cout <<endl<<"***********************************************\n";
                                }
                                h_CS->SetBinError(iW+1, sqrt(error));
                                var_sum += var[iPcanv][iThPiQ][iCharge][iMod][iQ2][iW];
                                if(var[iPcanv][iThPiQ][iCharge][iMod][iQ2][iW]!=0)
                                    var_num += y_ave[iPcanv][iThPiQ][iCharge][iMod][iQ2][iW];
                                CSvalue[iMod][iW] = h_CS->GetBinContent(iW+1);
                                CSerror[iMod][iW] = h_CS->GetBinError(iW+1);
                            }
                            //hs_CS->Add(h_CSerr);
                        }
                        else
                        {
                            for(Int_t iW = 0; iW < Wbini; iW++)
                            {
                                CSvalue[iMod][iW] = h_CS->GetBinContent(iW+1);
                                CSerror[iMod][iW] = h_CS->GetBinError(iW+1);
                            }
                        }
                        if(iMod==iDATA || iMod==iGENIE || iMod==iGENIEr || (iMod==iOPGn && bMultiPion==false))
                            {
                                hs_CS->Add(h_CS);
                            }
                    }
                    if(iCharge==0)
                    {
                        sprintf(c_NAME,"cs_files/pip_q%d_p%d_t%d.csv",iQ2+1,iPcanv+1,iThPiQ+1);
                    }
                    else if(iCharge==1)
                    {
                        sprintf(c_NAME,"cs_files/pim_q%d_p%d_t%d.csv",iQ2+1,iPcanv+1,iThPiQ+1);
                    }
                    ofstream CStxtFile (c_NAME);
                    if(CStxtFile.is_open())///For cross section values
                    {
                        //CStxtFile << "charge: "<< iCharge <<"\t Q2: "<< iQ2 <<"\t Pcanv: "<< iPcanv <<"\t theta: "<< iThPiQ <<"\n";
                        for(Int_t iW = 0; iW < Wbini; iW++)
                        {
                            Double_t Wvalue  = CS_Wmin + 0.025*(Double_t(iW)+1.);
                            CStxtFile << Wvalue << "," << CSvalue[iDATA][iW] << "," << CSerror[iDATA][iW] << "," << CSvalue[iGENIE][iW];
                            CStxtFile << "," << CSerror[iGENIE][iW] << "," << CSvalue[iOPGn][iW] << "," << CSerror[iOPGn][iW] << "\n";
                        }
                    }
                    CStxtFile.close();
                    hs_CS->Draw("nostack,e1p");
                    //hs_CS->Write("nostack,e1p");
                    hs_CS->GetYaxis()->SetLabelSize(label_size);//0.05 or label_size
                    hs_CS->SetMinimum(0.0);
                    hs_CS->GetYaxis()->SetNdivisions(3);
                    hs_CS->GetXaxis()->SetNdivisions(4);
                    //if(iPcanv == 0 && iThPiQ == 0 && iCharge == 0)
                    //{
                        hs_CS->GetXaxis()->SetTitleOffset(1.5);
                        hs_CS->GetXaxis()->SetLabelSize(xlab_size);
                    /*}
                    else if(iPcanv == 0 && iThPiQ == 0 && iCharge == 1)
                    {
                        hs_CS->GetXaxis()->SetTitleOffset(1.5);
                        hs_CS->GetXaxis()->SetLabelSize(label_size);
                    }
                    else if(iPcanv == 0)
                    {
                        hs_CS->GetXaxis()->SetTitleOffset(1.5);
                        hs_CS->GetXaxis()->SetLabelSize(label_size);
                    }*/
                    if(iThPiQ == 2 && iPcanv == Pbini - 1 && iQ2==0)
                    {
                        if(iCharge == 0){hs_CS->SetTitle("pi+   0.70 < Q2 < 1.0;;");}
                        else if(iCharge == 1){hs_CS->SetTitle("pi-   0.70 < Q2 < 1.0;;");}
                    }
                    else if(iThPiQ == 2 && iPcanv == Pbini - 1 && iQ2==1)
                    {
                        if(iCharge == 0){hs_CS->SetTitle("pi+   1.0 < Q2 < 1.4;;");}
                        else if(iCharge == 1){hs_CS->SetTitle("pi-   1.0 < Q2 < 1.4;;");}
                    }
                    else if(iThPiQ == 2 && iPcanv == Pbini - 1 && iQ2==2)
                    {
                        if(iCharge == 0){hs_CS->SetTitle("pi+   1.4 < Q2 < 1.9;;");}
                        else if(iCharge == 1){hs_CS->SetTitle("pi-   1.4 < Q2 < 1.9;;");}
                    }
                    else if(iThPiQ == 2 && iPcanv == Pbini - 1 && iQ2==3)
                    {
                        if(iCharge == 0){hs_CS->SetTitle("pi+   1.9 < Q2 < 2.5;;");}
                        else if(iCharge == 1){hs_CS->SetTitle("pi-   1.9 < Q2 < 2.5;;");}
                    }
                    /*hs_CS->GetYaxis()->SetLabelSize(0.05);//0.08
                    hs_CS->GetXaxis()->SetLabelSize(0.05);
                    hs_CS->SetTitle(";W (GeV);Cross Section (nb/sr/GeV^{4})");
                    hs_CS->GetYaxis()->SetTitleSize(0.05);
                    hs_CS->GetXaxis()->SetTitleSize(0.05);
                    hs_CS->GetYaxis()->SetTitleOffset(0.8);
                    hs_CS->GetXaxis()->SetTitleOffset(0.5);*/
                    //hs_CS->Write("hist && nostack,e1p");
                    //c1->Print(fileName);
                }
            }
            //if(iCharge < Nrec-1){c1->Print(fileName);}
            //else if(iQ2 < Q2bini-1)
            //{
                c1->Print(fileName);
                  if(iCharge == 0 && iQ2 == 0)      c1->Print(file_CSq1pip);
                  else if(iCharge == 0 && iQ2 == 1) c1->Print(file_CSq2pip);
                  else if(iCharge == 0 && iQ2 == 2) c1->Print(file_CSq3pip);
                  else if(iCharge == 0 && iQ2 == 3) c1->Print(file_CSq4pip);
                  else if(iCharge == 1 && iQ2 == 0) c1->Print(file_CSq1pim);
                  else if(iCharge == 1 && iQ2 == 1) c1->Print(file_CSq2pim);
                  else if(iCharge == 1 && iQ2 == 2) c1->Print(file_CSq3pim);
                  else if(iCharge == 1 && iQ2 == 3) c1->Print(file_CSq4pim);
            //}
        }
    }///End of cross sections

  //Print page to file
  c1->Clear();
  txtFile.close(); unc_pid.close(); unc_Vze.close(); unc_Vre.close(); unc_Vrpi.close(); unc_Vdiff.close(); ///close files

  /// -------------------------------------------------Q2------------------------------------------ ///
  //CStxtFile << "*******************************************************************************************" << endl;
  //CStxtFile << "2D binning: Q2" << endl;
  //CStxtFile << "*******************************************************************************************" << endl;
  //CStxtFile << "iMod\t iChrg\t iQ2\t W\t CS\t error" << endl;
  c1->Divide(2,2);
    line = 1, size = 3, marc = 2, mark_d = 8;
    //label_size = 0.11;//0.06 0.09
    for(Int_t iCharge = 0; iCharge < Nrec; iCharge++)/// cross sections
    {
        Double_t b2b_error = sqrt(pow(u_eID,2)+pow(u_hID,2)+pow(u_Vze,2)+pow(u_Vre,2)+pow(u_Vrpi,2)+pow(u_Vdiff,2)+pow(u_TW,2)+pow(RadFitUnc[iCharge],2));
        canvas = 1;
        for(Int_t iQ2 = 0; iQ2 < Q2bini; iQ2++)
        {
            Double_t CSvalue[models][Wbini];
            Double_t CSerror[models][Wbini];
            Int_t iQcanv = (Q2bini - 1) - iQ2;//Needed to get plot positions correct for Q2
            bin_vol = dQ2[iQcanv] * dW;
            c1->cd(canvas); canvas++; //Needed to get plot positions correct
            sprintf(c_NAME,"hs_CSq%d%d",iCharge,iQcanv);
            THStack * hs_CSq = new THStack(c_NAME,"");
            //gPad->SetLogy();
            for(Int_t iMod = 0; iMod < models; iMod++)
            {
                //if(iMod==iGENIE) continue;
                sprintf(c_NAME,"h_CSq%d%d%d",iCharge,iQcanv,iMod);
                TH1D * h_CSq = new TH1D(c_NAME,"",Wbini,CS_Wmin,CS_Wmax);
                h_CSq->Sumw2(); //https://root-forum.cern.ch/t/array-of-histograms/13148/2
                for(Int_t iW = 0; iW < Wbini; iW++)
                {
                    Double_t bin_content = 0., bin_err = 0.;
                    CSvalue[iMod][iW] = 0.;
                    CSerror[iMod][iW] = 0.;
                    Double_t ACorr = hr_ACq[iCharge][iQcanv]->GetBinContent(iW+1);
                    Double_t RCorr = hr_RCq[iCharge][iQcanv]->GetBinContent(iW+1);
                    if(RCorr<=0. && bRAD==true){cout<<"Negative RC"<<endl; RCorr=1.;}
                    if(bAC == false){ACorr=1.;}
                    if(bRAD == false){RCorr=1.;}
                    for(Int_t iPmom = 0; iPmom < Pbini+1; iPmom++)
                    {
                        Int_t iPcanv = iPmom;//(Pbini - 1) - iPmom;//Needed to get plot positions correct for P_pi
                        for(Int_t iThPiQ = 0; iThPiQ < TPIQi+1; iThPiQ++)
                        {
                            if(iMod<iGENIEr && y_sum[iPcanv][iThPiQ][iCharge][iMod][iQcanv][iW]!=0)///Scale h_CS bins by total CS for onepigen
                            {
                                bin_err += tCS[iPcanv][iThPiQ][iCharge][0][iMod][iQcanv][iW]/sqrt(y_sum[iPcanv][iThPiQ][iCharge][iMod][iQcanv][iW]);
                                bin_content += tCS[iPcanv][iThPiQ][iCharge][0][iMod][iQcanv][iW];
                            }
                            else if(iMod>=iGENIEr)
                            {
                                bin_content += y_sum[iPcanv][iThPiQ][iCharge][iMod][iQcanv][iW];
                            }
                        }
                    }
                    if(iMod<iGENIEr && bin_err!=0)///Scale h_CS bins by total CS for onepigen
                    {
                        h_CSq->SetBinContent(iW+1, bin_content * ACorr);
                        h_CSq->SetBinError(iW+1, bin_err);
                    }
                    else if(iMod>=iGENIEr)///Do not scale h_CS bins by total CS for non-onepigen (or bins have zero events)
                    {
                        h_CSq->SetBinContent(iW+1, bin_content * ACorr / RCorr);
                        h_CSq->SetBinError(iW+1, sqrt(bin_content));
                    }
                }   // Note: for an array bin # of iW, the ROOT bin # is iW+1
                //if (h_CS->GetSumw2N() == 0) h_CS->Sumw2(kTRUE);
                if(iMod==iDATA)//Data
                {
                    h_CSq->Scale(cm2ub * ub2nb / (L * bin_vol));//kBlack
                    //h_CSq->Scale(1. / RadFitAve[iCharge]);/// Radiative corrections
                    line = 1; mark = 8;  size = 5; marc = 1; h_CSq->SetLineColor(1);
                }
                else if(iMod==iGENIE)//GENIE
                {
                    h_CSq->Scale(totalCS_MC[iGENmod] * ub2nb / (MC_evts[iGENmod] * bin_vol));//kRed
                    //h_CS->Scale(pntF * ub2nb / (MC_evts[iGENmod] * bin_vol));
                    line = 1; mark = 32; size = 5; marc = 2; h_CSq->SetLineColor(2);
                }
                if(iMod==iGENIEr)//Data
                {
                    h_CSq->Scale(totalCS_MC[iGENmod] * ub2nb / (MC_rad_evts[iGENmod] * bin_vol));
                    line = 1; mark = 32;  size = 5; marc = 6; h_CSq->SetLineColor(6);
                }
                else if(iMod==iOPGr && iCharge==0)//onepigen pi+
                {
                    h_CSq->Scale(ub2nb / (opg_gen_events[0] * bin_vol));//kBlue
                    h_CSq->Scale(1. / RadFitAve[iCharge]);/// Radiative corrections
                    line = 1; mark = 26; size = 4; marc = 4; h_CSq->SetLineColor(4);
                }
                else if(iMod==iOPGr && iCharge==1)//onepigen pi-
                {
                    h_CSq->Scale(ub2nb / (opg_gen_events[1] * bin_vol));/// <-------------------<<
                    h_CSq->Scale(1. / RadFitAve[iCharge]);/// Radiative corrections
                    line = 1; mark = 26; size = 4; marc = 4; h_CSq->SetLineColor(4);//38
                }
                else if(iMod==iOPGn && iCharge==0)//onepigen pi+ norad
                {
                    h_CSq->Scale(ub2nb / (opg_gen_events[2] * bin_vol));//kBlue
                    line = 1; mark = 26; size = 4; marc = 4; h_CSq->SetLineColor(4);//marc=3
                }
                else if(iMod==iOPGn && iCharge==1)//onepigen pi- norad
                {
                    h_CSq->Scale(ub2nb / (opg_gen_events[3] * bin_vol));/// <-------------------<<
                    line = 1; mark = 26; size = 4; marc = 4; h_CSq->SetLineColor(4);//38);
                }
                h_CSq->SetLineStyle(line);
                h_CSq->SetMarkerStyle(mark);
                h_CSq->SetMarkerSize(4);//0.5 or size
                h_CSq->SetMarkerColor(marc);
                gStyle->SetEndErrorSize(10);
                if(iMod==iDATA)
                {
                    for(Int_t iW = 0; iW < Wbini; iW++)
                    {
                            Double_t error = pow(h_CSq->GetBinError(iW+1),2);       //statistical
                            //error += var[iPcanv][iThPiQ][iCharge][iMod][iQ2][iW];  //sector-by-sector uncertainty
                            error += pow(b2b_error * h_CSq->GetBinContent(iW+1),2); //bin-to-bin correction uncertainty
                            h_CSq->SetBinError(iW+1, sqrt(error));
                            CSvalue[iMod][iW] = h_CSq->GetBinContent(iW+1);
                            CSerror[iMod][iW] = h_CSq->GetBinError(iW+1);
                    }
                }
                else
                {
                    for(Int_t iW = 0; iW < Wbini; iW++)
                    {
                        CSvalue[iMod][iW] = h_CSq->GetBinContent(iW+1);
                        CSerror[iMod][iW] = h_CSq->GetBinError(iW+1);
                    }
                }
                if(iMod==iDATA || iMod==iGENIE || iMod==iGENIEr || (iMod==iOPGn && bMultiPion==false))
                    {
                        hs_CSq->Add(h_CSq);
                    }
                /*for(Int_t iW = 0; iW < Wbini; iW++)
                {
                    if(CStxtFile.is_open())///For cross section values
                    {
                        Double_t CSvalue = h_CSq->GetBinContent(iW+1);
                        Double_t CSerror = h_CSq->GetBinError(iW+1);
                        Double_t Wvalue  = CS_Wmin + 0.025*(Double_t(iW)+1.);
                        CStxtFile << iMod <<"\t "<< iCharge <<"\t "<< iQ2 <<"\t "<< Wvalue <<"\t "<< CSvalue <<"\t "<< CSerror << endl;
                    }
                }*/
            }
            if(iCharge==0)
            {
                sprintf(c_NAME,"cs_files/pip_q%d.csv",iQ2+1);
            }
            else if(iCharge==1)
            {
                sprintf(c_NAME,"cs_files/pim_q%d.csv",iQ2+1);
            }
            ofstream CStxtFile (c_NAME);
            if(CStxtFile.is_open())///For cross section values
            {
                for(Int_t iW = 0; iW < Wbini; iW++)
                {
                    Double_t Wvalue  = CS_Wmin + 0.025*(Double_t(iW)+1.);
                    CStxtFile << Wvalue << "," << CSvalue[iDATA][iW] << "," << CSerror[iDATA][iW] << "," << CSvalue[iGENIE][iW];
                    CStxtFile << "," << CSerror[iGENIE][iW] << "," << CSvalue[iOPGn][iW] << "," << CSerror[iOPGn][iW] << "\n";
                }
            }
            CStxtFile.close();
            hs_CSq->Draw("nostack,e1p");
            hs_CSq->GetYaxis()->SetLabelSize(label_size);//0.05 or label_size
            hs_CSq->GetYaxis()->SetNdivisions(3);
            hs_CSq->GetXaxis()->SetNdivisions(4);
            hs_CSq->GetXaxis()->SetTitleOffset(1.5);
            hs_CSq->GetXaxis()->SetLabelSize(xlab_size);
            if(iQcanv==0)
            {
                if(iCharge == 0){hs_CSq->SetTitle("pi+   0.70 < Q2 < 1.0;;");}
                else if(iCharge == 1){hs_CSq->SetTitle("pi-   0.70 < Q2 < 1.0;;");}
            }
            else if(iQcanv==1)
            {
                if(iCharge == 0){hs_CSq->SetTitle("pi+   1.0 < Q2 < 1.4;;");}
                else if(iCharge == 1){hs_CSq->SetTitle("pi-   1.0 < Q2 < 1.4;;");}
            }
            else if(iQcanv==2)
            {
                if(iCharge == 0){hs_CSq->SetTitle("pi+   1.4 < Q2 < 1.9;;");}
                else if(iCharge == 1){hs_CSq->SetTitle("pi-   1.4 < Q2 < 1.9;;");}
            }
            else if(iQcanv==3)
            {
                if(iCharge == 0){hs_CSq->SetTitle("pi+   1.9 < Q2 < 2.5;;");}
                else if(iCharge == 1){hs_CSq->SetTitle("pi-   1.9 < Q2 < 2.5;;");}
            }
        }
        c1->Print(fileName);
    }///End of cross sections
  //Print page to file
  c1->Clear();
  /// -------------------------------------------------P_PI------------------------------------------ ///
  //CStxtFile << "*******************************************************************************************" << endl;
  //CStxtFile << "3D binning: Q2 and P_pi" << endl;
  //CStxtFile << "*******************************************************************************************" << endl;
  //CStxtFile << "iMod\t iChrg\t iQ2\t iPcanv\t W\t CS\t error" << endl;
  c1->Divide(Pbini,Q2bini);
    line = 1, size = 3, marc = 2, mark_d = 8;
    //label_size = 0.11;//0.06 0.09
    for(Int_t iCharge = 0; iCharge < Nrec; iCharge++)/// cross sections
    {
        Double_t b2b_error = sqrt(pow(u_eID,2)+pow(u_hID,2)+pow(u_Vze,2)+pow(u_Vre,2)+pow(u_Vrpi,2)+pow(u_Vdiff,2)+pow(u_TW,2)+pow(RadFitUnc[iCharge],2));
        canvas = 1;
        for(Int_t iQ2 = 0; iQ2 < Q2bini; iQ2++)
        {
            Int_t iQcanv = (Q2bini - 1) - iQ2;//Needed to get plot positions correct for Q2
            for(Int_t iPmom = 0; iPmom < Pbini; iPmom++)
            {
                Int_t iPcanv = iPmom;//(Pbini - 1) - iPmom;//Needed to get plot positions correct for P_pi
                Double_t dPmom = dP[iPcanv];
                Double_t CSvalue[models][Wbini];
                Double_t CSerror[models][Wbini];
                if(bQorNot==false)
                {
                    dPmom = dPt[iCharge][iPmom+1] - dPt[iCharge][iPmom];
                }
                bin_vol = dPmom * dT * dQ2[iQcanv] * dW;
                c1->cd(canvas); canvas++; //Needed to get plot positions correct
                sprintf(c_NAME,"hs_CSp%d%d%d",iPcanv,iCharge,iQcanv);
                THStack * hs_CSp = new THStack(c_NAME,"");
                //gPad->SetLogy();
                for(Int_t iMod = 0; iMod < models; iMod++)
                {
                    //if(iMod==iGENIE) continue;
                    sprintf(c_NAME,"h_CSp%d%d%d%d",iPcanv,iCharge,iQcanv,iMod);
                    TH1D * h_CSp = new TH1D(c_NAME,"",Wbini,CS_Wmin,CS_Wmax);
                    h_CSp->Sumw2(); //https://root-forum.cern.ch/t/array-of-histograms/13148/2
                    for(Int_t iW = 0; iW < Wbini; iW++)
                    {
                        Double_t bin_content = 0., bin_err = 0.;
                        CSvalue[iMod][iW] = 0.;
                        CSerror[iMod][iW] = 0.;
                        Double_t ACorr = hr_ACp[iPcanv][iCharge][iQcanv]->GetBinContent(iW+1);
                        Double_t RCorr = hr_RCp[iPcanv][iCharge][iQcanv]->GetBinContent(iW+1);
                        if(RCorr<=0. && bRAD==true){cout<<"Negative RC"<<endl; RCorr=1.;}
                        if(bAC == false){ACorr=1.;}
                        if(bRAD == false){RCorr=1.;}
                        for(Int_t iThPiQ = 0; iThPiQ < TPIQi+1; iThPiQ++)
                        {
                            if(iMod<iGENIEr && y_sum[iPcanv][iThPiQ][iCharge][iMod][iQcanv][iW]!=0)///Scale h_CS bins by total CS for onepigen
                            {
                                bin_err += tCS[iPcanv][iThPiQ][iCharge][0][iMod][iQcanv][iW]/sqrt(y_sum[iPcanv][iThPiQ][iCharge][iMod][iQcanv][iW]);
                                bin_content += tCS[iPcanv][iThPiQ][iCharge][0][iMod][iQcanv][iW];
                            }
                            else if(iMod>=iGENIEr)
                            {
                                bin_content += y_sum[iPcanv][iThPiQ][iCharge][iMod][iQcanv][iW];
                            }
                        }
                        if(iMod<iGENIEr && bin_err!=0)///Scale h_CS bins by total CS for onepigen
                        {
                            h_CSp->SetBinContent(iW+1, bin_content * ACorr);
                            h_CSp->SetBinError(iW+1, bin_err);
                        }
                        else if(iMod>=iGENIEr)///Do not scale h_CS bins by total CS for non-onepigen (or bins have zero events)
                        {
                            h_CSp->SetBinContent(iW+1, bin_content * ACorr / RCorr);
                            h_CSp->SetBinError(iW+1, sqrt(bin_content));
                        }
                    }   // Note: for an array bin # of iW, the ROOT bin # is iW+1
                    //if (h_CS->GetSumw2N() == 0) h_CS->Sumw2(kTRUE);
                    if(iMod==iDATA)//Data
                    {
                        h_CSp->Scale(cm2ub * ub2nb / (L * bin_vol));//kBlack
                        //h_CSp->Scale(1. / RadFitAve[iCharge]);/// Radiative corrections
                        line = 1; mark = 8;  size = 5; marc = 1; h_CSp->SetLineColor(1);
                    }
                    else if(iMod==iGENIE)//GENIE
                    {
                        h_CSp->Scale(totalCS_MC[iGENmod] * ub2nb / (MC_evts[iGENmod] * bin_vol));//kRed
                        //h_CS->Scale(pntF * ub2nb / (MC_evts[iGENmod] * bin_vol));
                        line = 1; mark = 32; size = 5; marc = 2; h_CSp->SetLineColor(2);
                    }
                    if(iMod==iGENIEr)//Data
                    {
                        h_CSp->Scale(totalCS_MC[iGENmod] * ub2nb / (MC_rad_evts[iGENmod] * bin_vol));
                        line = 1; mark = 32;  size = 5; marc = 6; h_CSp->SetLineColor(6);
                    }
                    else if(iMod==iOPGr && iCharge==0)//onepigen pi+
                    {
                        h_CSp->Scale(ub2nb / (opg_gen_events[0] * bin_vol));//kBlue
                        //h_CSp->Scale(1. / RadFitAve[iCharge]);/// Radiative corrections
                        line = 1; mark = 26; size = 4; marc = 4; h_CSp->SetLineColor(4);
                    }
                    else if(iMod==iOPGr && iCharge==1)//onepigen pi-
                    {
                        h_CSp->Scale(ub2nb / (opg_gen_events[1] * bin_vol));/// <-------------------<<
                        h_CSp->Scale(1. / RadFitAve[iCharge]);/// Radiative corrections
                        line = 1; mark = 26; size = 4; marc = 4; h_CSp->SetLineColor(4);//38
                    }
                    else if(iMod==iOPGn && iCharge==0)//onepigen pi+ norad
                    {
                        h_CSp->Scale(ub2nb / (opg_gen_events[2] * bin_vol));//kBlue
                        line = 1; mark = 26; size = 4; marc = 4; h_CSp->SetLineColor(4);//marc=3
                    }
                    else if(iMod==iOPGn && iCharge==1)//onepigen pi- norad
                    {
                        h_CSp->Scale(ub2nb / (opg_gen_events[3] * bin_vol));/// <-------------------<<
                        line = 1; mark = 26; size = 4; marc = 4; h_CSp->SetLineColor(4);//38);
                    }
                    h_CSp->SetLineStyle(line);
                    h_CSp->SetMarkerStyle(mark);
                    h_CSp->SetMarkerSize(2);//0.5 or size
                    h_CSp->SetMarkerColor(marc);
                    gStyle->SetEndErrorSize(10);
                    if(iMod==iDATA)
                    {
                        for(Int_t iW = 0; iW < Wbini; iW++)
                        {
                                Double_t error = pow(h_CSp->GetBinError(iW+1),2);       //statistical
                                //error += var[iPcanv][iThPiQ][iCharge][iMod][iQ2][iW];  //sector-by-sector uncertainty
                                error += pow(b2b_error * h_CSp->GetBinContent(iW+1),2); //bin-to-bin correction uncertainty
                                h_CSp->SetBinError(iW+1, sqrt(error));
                                CSvalue[iMod][iW] = h_CSp->GetBinContent(iW+1);
                                CSerror[iMod][iW] = h_CSp->GetBinError(iW+1);
                        }
                    }
                    else
                    {
                        for(Int_t iW = 0; iW < Wbini; iW++)
                        {
                            CSvalue[iMod][iW] = h_CSp->GetBinContent(iW+1);
                            CSerror[iMod][iW] = h_CSp->GetBinError(iW+1);
                        }
                    }
                    if(iMod==iDATA || iMod==iGENIE || iMod==iGENIEr || (iMod==iOPGn && bMultiPion==false))
                        {
                            hs_CSp->Add(h_CSp);
                        }
                    /*for(Int_t iW = 0; iW < Wbini; iW++)
                    {
                        if(CStxtFile.is_open())///For cross section values
                        {
                            Double_t CSvalue = h_CSp->GetBinContent(iW+1);
                            Double_t CSerror = h_CSp->GetBinError(iW+1);
                            Double_t Wvalue  = CS_Wmin + 0.025*(Double_t(iW)+1.);
                            CStxtFile << iMod <<"\t "<< iCharge <<"\t "<< iQ2 <<"\t "<< iPcanv <<"\t "<< Wvalue <<"\t "<< CSvalue <<"\t "<< CSerror << endl;
                        }
                    }*/
                }
                if(iCharge==0)
                {
                    sprintf(c_NAME,"cs_files/pip_q%d_p%d.csv",iQ2+1,iPcanv+1);
                }
                else if(iCharge==1)
                {
                    sprintf(c_NAME,"cs_files/pim_q%d_p%d.csv",iQ2+1,iPcanv+1);
                }
                ofstream CStxtFile (c_NAME);
                if(CStxtFile.is_open())///For cross section values
                {
                    for(Int_t iW = 0; iW < Wbini; iW++)
                    {
                        Double_t Wvalue  = CS_Wmin + 0.025*(Double_t(iW)+1.);
                        CStxtFile << Wvalue << "," << CSvalue[iDATA][iW] << "," << CSerror[iDATA][iW] << "," << CSvalue[iGENIE][iW];
                        CStxtFile << "," << CSerror[iGENIE][iW] << "," << CSvalue[iOPGn][iW] << "," << CSerror[iOPGn][iW] << "\n";
                    }
                }
                CStxtFile.close();
                hs_CSp->Draw("nostack,e1p");
                hs_CSp->GetYaxis()->SetLabelSize(label_size);//0.05 or label_size
                hs_CSp->GetYaxis()->SetNdivisions(3);
                hs_CSp->GetXaxis()->SetNdivisions(4);
                hs_CSp->GetXaxis()->SetTitleOffset(1.5);
                hs_CSp->GetXaxis()->SetLabelSize(xlab_size);
                if(iPcanv == Pbini - 1 && iQcanv==0)
                {
                    if(iCharge == 0){hs_CSp->SetTitle("P_pi pi+   0.70 < Q2 < 1.0;;");}
                    else if(iCharge == 1){hs_CSp->SetTitle("P_pi pi-   0.70 < Q2 < 1.0;;");}
                }
                else if(iPcanv == Pbini - 1 && iQcanv==1)
                {
                    if(iCharge == 0){hs_CSp->SetTitle("pi+   1.0 < Q2 < 1.4;;");}
                    else if(iCharge == 1){hs_CSp->SetTitle("pi-   1.0 < Q2 < 1.4;;");}
                }
                else if(iPcanv == Pbini - 1 && iQcanv==2)
                {
                    if(iCharge == 0){hs_CSp->SetTitle("pi+   1.4 < Q2 < 1.9;;");}
                    else if(iCharge == 1){hs_CSp->SetTitle("pi-   1.4 < Q2 < 1.9;;");}
                }
                else if(iPcanv == Pbini - 1 && iQcanv==3)
                {
                    if(iCharge == 0){hs_CSp->SetTitle("pi+   1.9 < Q2 < 2.5;;");}
                    else if(iCharge == 1){hs_CSp->SetTitle("pi-   1.9 < Q2 < 2.5;;");}
                }
            }
        }
        c1->Print(fileName);
    }///End of cross sections
  //Print page to file
  c1->Clear();
  /// -------------------------------------------------THETA_PIQ------------------------------------------ ///
  //CStxtFile << "*******************************************************************************************" << endl;
  //if(bQorNot==true){CStxtFile << "3D binning: Q2 and Theta_piq" << endl;}
  //else             {CStxtFile << "3D binning: Q2 and Theta_pi" << endl;}
  //CStxtFile << "*******************************************************************************************" << endl;
  //CStxtFile << "iMod\t iChrg\t iQ2\t iTheta\t W\t CS\t error" << endl;
  c1->Divide(TPIQi,Q2bini);
    line = 1, size = 3, marc = 2, mark_d = 8;
    //label_size = 0.11;//0.06 0.09
    for(Int_t iCharge = 0; iCharge < Nrec; iCharge++)/// cross sections
    {
        Double_t b2b_error = sqrt(pow(u_eID,2)+pow(u_hID,2)+pow(u_Vze,2)+pow(u_Vre,2)+pow(u_Vrpi,2)+pow(u_Vdiff,2)+pow(u_TW,2)+pow(RadFitUnc[iCharge],2));
        canvas = 1;
        for(Int_t iQ2 = 0; iQ2 < Q2bini; iQ2++)
        {
            Int_t iQcanv = (Q2bini - 1) - iQ2;//Needed to get plot positions correct for Q2
            for(Int_t iThPiQ = 0; iThPiQ < TPIQi; iThPiQ++)
            {
                Double_t CSvalue[models][Wbini];
                Double_t CSerror[models][Wbini];
                dT = DeltaThetaPiQBin(ThetaPiQ[iThPiQ],ThetaPiQ[iThPiQ+1]);
                if(bQorNot==false)
                {
                    dT = DeltaThetaPiQBin(ThetaPi[iCharge][iThPiQ],ThetaPi[iCharge][iThPiQ+1]);
                }
                bin_vol = dT * dQ2[iQcanv] * dW;
                c1->cd(canvas); canvas++; //Needed to get plot positions correct
                sprintf(c_NAME,"hs_CSt%d%d%d",iThPiQ,iCharge,iQcanv);
                THStack * hs_CSt = new THStack(c_NAME,"");
                //gPad->SetLogy();
                for(Int_t iMod = 0; iMod < models; iMod++)
                {
                    //if(iMod==iGENIE) continue;
                    sprintf(c_NAME,"h_CSt%d%d%d%d",iThPiQ,iCharge,iQcanv,iMod);
                    TH1D * h_CSt = new TH1D(c_NAME,"",Wbini,CS_Wmin,CS_Wmax);
                    h_CSt->Sumw2(); //https://root-forum.cern.ch/t/array-of-histograms/13148/2
                    for(Int_t iW = 0; iW < Wbini; iW++)
                    {
                        Double_t bin_content = 0, bin_err = 0;
                        CSvalue[iMod][iW] = 0.;
                        CSerror[iMod][iW] = 0.;
                        Double_t ACorr = hr_ACt[iThPiQ][iCharge][iQcanv]->GetBinContent(iW+1);
                        Double_t RCorr = hr_RCt[iThPiQ][iCharge][iQcanv]->GetBinContent(iW+1);
                        if(RCorr<=0. && bRAD==true){cout<<"Negative RC"<<endl; RCorr=1.;}
                        if(bAC == false){ACorr=1.;}
                        if(bRAD == false){RCorr=1.;}
                        for(Int_t iPmom = 0; iPmom < Pbini+1; iPmom++)
                        {
                            Int_t iPcanv = iPmom;
                            if(iMod<iGENIEr && y_sum[iPcanv][iThPiQ][iCharge][iMod][iQcanv][iW]!=0)///Scale h_CS bins by total CS for onepigen
                            {
                                bin_err += tCS[iPcanv][iThPiQ][iCharge][0][iMod][iQcanv][iW]/sqrt(y_sum[iPcanv][iThPiQ][iCharge][iMod][iQcanv][iW]);
                                bin_content += tCS[iPcanv][iThPiQ][iCharge][0][iMod][iQcanv][iW];
                            }
                            else if(iMod>=iGENIEr)
                            {
                                bin_content += y_sum[iPcanv][iThPiQ][iCharge][iMod][iQcanv][iW];
                            }
                        }
                        if(iMod<iGENIEr && bin_err!=0)///Scale h_CS bins by total CS for onepigen
                        {
                            h_CSt->SetBinContent(iW+1, bin_content * ACorr);
                            h_CSt->SetBinError(iW+1, bin_err);
                        }
                        else if(iMod>=iGENIEr)///Do not scale h_CS bins by total CS for non-onepigen (or bins have zero events)
                        {
                            h_CSt->SetBinContent(iW+1, bin_content * ACorr / RCorr);
                            h_CSt->SetBinError(iW+1, sqrt(bin_content));
                        }
                    }   // Note: for an array bin # of iW, the ROOT bin # is iW+1
                    //if (h_CS->GetSumw2N() == 0) h_CS->Sumw2(kTRUE);
                    if(iMod==iDATA)//Data
                    {
                        h_CSt->Scale(cm2ub * ub2nb / (L * bin_vol));//kBlack
                        //h_CSt->Scale(1. / RadFitAve[iCharge]);/// Radiative corrections
                        line = 1; mark = 8;  size = 5; marc = 1; h_CSt->SetLineColor(1);
                    }
                    else if(iMod==iGENIE)//GENIE
                    {
                        h_CSt->Scale(totalCS_MC[iGENmod] * ub2nb / (MC_evts[iGENmod] * bin_vol));//kRed
                        //h_CS->Scale(pntF * ub2nb / (MC_evts[iGENmod] * bin_vol));
                        line = 1; mark = 32; size = 5; marc = 2; h_CSt->SetLineColor(2);
                    }
                    if(iMod==iGENIEr)//Data
                    {
                        h_CSt->Scale(totalCS_MC[iGENmod] * ub2nb / (MC_rad_evts[iGENmod] * bin_vol));
                        line = 1; mark = 32;  size = 5; marc = 6; h_CSt->SetLineColor(6);
                    }
                    else if(iMod==iOPGr && iCharge==0)//onepigen pi+
                    {
                        h_CSt->Scale(ub2nb / (opg_gen_events[0] * bin_vol));//kBlue
                        h_CSt->Scale(1. / RadFitAve[iCharge]);/// Radiative corrections
                        line = 1; mark = 26; size = 4; marc = 4; h_CSt->SetLineColor(4);
                    }
                    else if(iMod==iOPGr && iCharge==1)//onepigen pi-
                    {
                        h_CSt->Scale(ub2nb / (opg_gen_events[1] * bin_vol));/// <-------------------<<
                        h_CSt->Scale(1. / RadFitAve[iCharge]);/// Radiative corrections
                        line = 1; mark = 26; size = 4; marc = 4; h_CSt->SetLineColor(4);//38
                    }
                    else if(iMod==iOPGn && iCharge==0)//onepigen pi+ norad
                    {
                        h_CSt->Scale(ub2nb / (opg_gen_events[2] * bin_vol));//kBlue
                        line = 1; mark = 26; size = 4; marc = 4; h_CSt->SetLineColor(4);//marc=3
                    }
                    else if(iMod==iOPGn && iCharge==1)//onepigen pi- norad
                    {
                        h_CSt->Scale(ub2nb / (opg_gen_events[3] * bin_vol));/// <-------------------<<
                        line = 1; mark = 26; size = 4; marc = 4; h_CSt->SetLineColor(4);//38);
                    }
                    h_CSt->SetLineStyle(line);
                    h_CSt->SetMarkerStyle(mark);
                    h_CSt->SetMarkerSize(2);//0.5 or size
                    h_CSt->SetMarkerColor(marc);
                    gStyle->SetEndErrorSize(10);
                    if(iMod==iDATA)
                    {
                        for(Int_t iW = 0; iW < Wbini; iW++)
                        {
                                Double_t error = pow(h_CSt->GetBinError(iW+1),2);       //statistical
                                //error += var[iPcanv][iThPiQ][iCharge][iMod][iQ2][iW];  //sector-by-sector uncertainty
                                error += pow(b2b_error * h_CSt->GetBinContent(iW+1),2); //bin-to-bin correction uncertainty
                                h_CSt->SetBinError(iW+1, sqrt(error));
                                CSvalue[iMod][iW] = h_CSt->GetBinContent(iW+1);
                                CSerror[iMod][iW] = h_CSt->GetBinError(iW+1);
                        }
                    }
                    else
                    {
                        for(Int_t iW = 0; iW < Wbini; iW++)
                        {
                            CSvalue[iMod][iW] = h_CSt->GetBinContent(iW+1);
                            CSerror[iMod][iW] = h_CSt->GetBinError(iW+1);
                        }
                    }
                    if(iMod==iDATA || iMod==iGENIE || iMod==iGENIEr || (iMod==iOPGn && bMultiPion==false))
                        {
                            hs_CSt->Add(h_CSt);
                        }
                    /*for(Int_t iW = 0; iW < Wbini; iW++)
                    {
                        if(CStxtFile.is_open())///For cross section values
                        {
                            Double_t CSvalue = h_CSt->GetBinContent(iW+1);
                            Double_t CSerror = h_CSt->GetBinError(iW+1);
                            Double_t Wvalue  = CS_Wmin + 0.025*(Double_t(iW)+1.);
                            CStxtFile << iMod <<"\t "<< iCharge <<"\t "<< iQ2 <<"\t "<< iThPiQ <<"\t "<< Wvalue <<"\t "<< CSvalue <<"\t "<< CSerror << endl;
                        }
                    }*/
                }
                if(iCharge==0)
                {
                    sprintf(c_NAME,"cs_files/pip_q%d_t%d.csv",iQ2+1,iThPiQ+1);
                }
                else if(iCharge==1)
                {
                    sprintf(c_NAME,"cs_files/pim_q%d_t%d.csv",iQ2+1,iThPiQ+1);
                }
                ofstream CStxtFile (c_NAME);
                if(CStxtFile.is_open())///For cross section values
                {
                    for(Int_t iW = 0; iW < Wbini; iW++)
                    {
                        Double_t Wvalue  = CS_Wmin + 0.025*(Double_t(iW)+1.);
                        CStxtFile << Wvalue << "," << CSvalue[iDATA][iW] << "," << CSerror[iDATA][iW] << "," << CSvalue[iGENIE][iW];
                        CStxtFile << "," << CSerror[iGENIE][iW] << "," << CSvalue[iOPGn][iW] << "," << CSerror[iOPGn][iW] << "\n";
                    }
                }
                CStxtFile.close();
                hs_CSt->Draw("nostack,e1p");
                hs_CSt->GetYaxis()->SetLabelSize(label_size);//0.05 or label_size
                hs_CSt->GetYaxis()->SetNdivisions(3);
                hs_CSt->GetXaxis()->SetNdivisions(4);
                hs_CSt->GetXaxis()->SetTitleOffset(1.5);
                hs_CSt->GetXaxis()->SetLabelSize(xlab_size);
                if(iQcanv==0)
                {
                    if(iCharge == 0){hs_CSt->SetTitle("Theta_piq pi+   0.70 < Q2 < 1.0;;");}
                    else if(iCharge == 1){hs_CSt->SetTitle("Theta_piq pi-   0.70 < Q2 < 1.0;;");}
                }
                else if(iQcanv==1)
                {
                    if(iCharge == 0){hs_CSt->SetTitle("pi+   1.0 < Q2 < 1.4;;");}
                    else if(iCharge == 1){hs_CSt->SetTitle("pi-   1.0 < Q2 < 1.4;;");}
                }
                else if(iQcanv==2)
                {
                    if(iCharge == 0){hs_CSt->SetTitle("pi+   1.4 < Q2 < 1.9;;");}
                    else if(iCharge == 1){hs_CSt->SetTitle("pi-   1.4 < Q2 < 1.9;;");}
                }
                else if(iQcanv==3)
                {
                    if(iCharge == 0){hs_CSt->SetTitle("pi+   1.9 < Q2 < 2.5;;");}
                    else if(iCharge == 1){hs_CSt->SetTitle("pi-   1.9 < Q2 < 2.5;;");}
                }
            }
        }
        c1->Print(fileName);
    }///End of cross sections
  //Print page to file
  c1->Clear();
  //CStxtFile.close();

  /// RC vs. chi2/NDF
  //Print final page to file and close pdf
  c1->Divide(2,1);
  c1->cd(1);
  h2_RC_chi2pip->Draw("colz");
  c1->cd(2);
  h2_RC_chi2pim->Draw("colz");
  c1->Print(fileName+")");
  c1->Clear();

  ///Write out histograms to file
  hs_W_data->Write("nostack");
  hs_W_resid_MC->Write();
  //h_weight_MCrad->Write();
  //h_processid_MCrad->Write();
  hs_q2_data->Write("nostack");
  hs_theta_pi_data->Write("nostack");
  hs_theta_piq_data->Write("nostack");
  hs_Ppi_data->Write("nostack");
  for(Int_t id=0;id<20;id++)
  {
      h_W_resid_MC[id]->Write();
  }
  h2_DC_y_vs_x_data->Write();
  h2_phi_e_phi_pip->Write();
  h2_E_FTOF_P_pip->Write();
  h2_E_FTOF_P_pipQE->Write();
  h2_E_FTOF_P_pipz->Write();
  h_Wpip_data_chi2->Write();
  h_Wpip_data_chi3->Write();
  h2_beta_P_pipchi2->Write();
  h2_beta_P_pipchi3->Write();
  h_Wpip_data_ray->Write();
  h_Betapip_data_ray->Write();
  h2_beta_P_pip_ray->Write();
  h2_beta_P_pip->Write();
  h2_beta_P_pipzoom->Write();
  h2_beta_P_pipHiTQ->Write();
  h2_beta_P_pipLowW->Write();
  h2_beta_P_pim->Write();
  h_beta_P8_pip->Write();
  h_beta_P10_pip->Write();
  h_beta_P8_pipcut->Write();
  h_beta_P10_pipcut->Write();
  h_beta_P8_pim->Write();
  h_beta_P10_pim->Write();
  h_beta_P8_pimcut->Write();
  h_beta_P10_pimcut->Write();
  /*h2_beta_P_p->Write();
  h2_beta_t_FTOF->Write();
  h2_q2_Wpip_QE->Write();
  h2_E_FTOFa_P_pip->Write();
  h2_beta_t_FTOFa->Write();
  h_W_data_pipFTOF1B->Write();
  h_W_data_pipFTOF1A->Write();
  h_W_data_pipFTOFnoE->Write();
  h_Wpip_data_hiTQ->Write();
  h_ThetaQpip_dataQE->Write();
  h_ThetaQpip_dataNoQE->Write();
  h_ThetaQpip_dataFTOF1B->Write();
  h2_phi_e_phi_pim->Write();
  h2_E_FTOF_P_pim->Write();*/
  h_pip_chi2pid->Write();
  h_pim_chi2pid->Write();
  q2_stacked->Write();
  h_q2_data->Write();
  h_q2_MC->Write();
  el_vz_h_stacked->Write();
  h_el_vz_h_data->Write();
  h_el_vz_h_MC->Write();
  W_stacked->Write();
  h_W_data->Write();
  h_Wp_data->Write();
  h_Wpip_data->Write();
  h_W_MC->Write();
  h2_q2_W_data->Write();
  h2_Ppip_thetapipq_data->Write();
  h2_Ppim_thetapimq_data->Write();
  h2_Ppip_thetapip_data->Write();
  h2_Ppim_thetapim_data->Write();


  gBenchmark->Stop("timer");
  gBenchmark->Print("timer");

  auto finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = finish - start;
  std::cout << "Elapsed time: " << elapsed.count() << std:: endl;
  std::cout << "Data events = " << counter_data << ", and MC events = " << counter_MC << " MC1= " << counter_MC1 << std::endl;
  std::cout << "Data CS events = " << cs_data << ", and MC CS events = " << cs_MC << std::endl;
  std::cout << "MC pi+ events =  " << i_pip   << ", and pi- events =   " << i_pim << std::endl;
  std::cout << "opg events (+)=  " << cs_opg   << ", " << counter_opg<< ", (-) " << cs_opgm<< ", " << counter_opgm << std::endl;
  std::cout << "opgpr file events =  " << i_pisum+i_pipr << " = " << i_pitot << "\t Off by:  " << i_pitot-i_pisum-i_pipr << std::endl;
  std::cout << "opgmr file events =  " << i_pisut+i_pimr << " = " << i_pitom << "\t Off by:  " << i_pitom-i_pisut-i_pimr << std::endl;
  std::cout << "opgpn file events =  " << i_pisup+i_pipn << " = " << i_pitop << "\t Off by:  " << i_pitop-i_pisup-i_pipn << std::endl;
  std::cout << "opgmn file events =  " << i_pisun+i_pimn << " = " << i_piton << "\t Off by:  " << i_piton-i_pisun-i_pimn << std::endl;
  std::cout << "RadFitAve (pi+) =  " << RadFitAve[0] <<" +- "<<RadFitUnc[0]<<"\t (pi-) =  " << RadFitAve[1]  <<" +- "<<RadFitUnc[1]<< std::endl;
  std::cout << "Sector uncer =  " << 100*sqrt(var_sum)/var_num <<"%"<< std::endl;
  //std::cout << "tpd events =  " << cs_tpdT   << ", " << cs_tpd2<< ", " << cs_tpd3<< ", " << counter_tpd << std::endl;
}

Int_t eCutter(Double_t VWcut, Double_t Lv, Double_t Lw, Double_t vz_e, Double_t vz_e_min, Double_t vz_e_max, Double_t Vperp_e, Double_t Vperp_ecut)
{
    Int_t cutty = 0;
    if (
        Lv <= VWcut || Lw <= VWcut ||           //e- V,W cuts
        vz_e < vz_e_min || vz_e > vz_e_max ||   //e- vertex z cuts
        Vperp_e > Vperp_ecut                    //e- vertex r cuts
        ) {cutty = 1;}
    return cutty;
}

Int_t PiCutter(Int_t iCharge, Double_t P_pi, Double_t p5, Double_t p10, Double_t Vperp_pi, Double_t Vdiff_pi,
              Double_t Vperp_pipcut[3], Double_t Vperp_pip_ms[3][2], Double_t Vperp_pimcut[2], Double_t Vperp_pim_ms[2][2], Double_t Nsig)
{
    Int_t cutty = 0;
    if (
        iCharge == 0 &&
        (
            (P_pi < p5 &&                 (Vperp_pi>Vperp_pipcut[0] || Vdiff_pi<Vperp_pip_ms[0][0]-Nsig*Vperp_pip_ms[0][1] || Vdiff_pi>Vperp_pip_ms[0][0]+Nsig*Vperp_pip_ms[0][1])) ||
            ((P_pi > p5 && P_pi < p10) && (Vperp_pi>Vperp_pipcut[1] || Vdiff_pi<Vperp_pip_ms[1][0]-Nsig*Vperp_pip_ms[1][1] || Vdiff_pi>Vperp_pip_ms[1][1]+Nsig*Vperp_pip_ms[1][1])) ||
            (P_pi > p10 &&                (Vperp_pi>Vperp_pipcut[2] || Vdiff_pi<Vperp_pip_ms[2][0]-Nsig*Vperp_pip_ms[2][1] || Vdiff_pi>Vperp_pip_ms[2][1]+Nsig*Vperp_pip_ms[2][1]))
        )
       ) {cutty = 2;}
    else if (
        iCharge == 1 &&
        (
            (P_pi < p10 &&                 (Vperp_pi>Vperp_pimcut[0] || Vdiff_pi<Vperp_pim_ms[0][0]-Nsig*Vperp_pim_ms[0][1] || Vdiff_pi>Vperp_pim_ms[0][0]+Nsig*Vperp_pim_ms[0][1])) ||
            (P_pi > p10 &&                 (Vperp_pi>Vperp_pimcut[1] || Vdiff_pi<Vperp_pim_ms[1][0]-Nsig*Vperp_pim_ms[1][1] || Vdiff_pi>Vperp_pim_ms[1][0]+Nsig*Vperp_pim_ms[1][1]))
        )
       ) {cutty = 3;}
    else if (iCharge != 0 && iCharge != 1) {cutty = 4;}
    return cutty;
}

Double_t DeltaThetaPiQBin(Double_t ThetaMin, Double_t ThetaMax)
{
    Double_t answer = 2*3.1415926536*(cos(ThetaMin*3.1415926536/180) - cos(ThetaMax*3.1415926536/180));
    return answer;
}


Double_t ExpFunc(Double_t A, Double_t B, Double_t C, Double_t D, Double_t E, Double_t P)
{
    Double_t answer = A - B*pow(E,C-D*P);
    return answer;
}

void EventBinner(bool QorNot, Int_t ichrg, Double_t P_pi, Double_t theta, Double_t q2_epi, Double_t W_epi, Double_t Wmin, Double_t DeltaW, Int_t Wbini,
                 Int_t &iPmom, Int_t &iThPiQ, Int_t &iQ2, Int_t &iW)
{
    if(QorNot == true)//theta_piq
    {
        Int_t iAllPnTheta = -1;//-1 or 4
        //Selecting pi momentum (P_pi) bin
        Double_t Pbin[2][5] = {{0.7,1.0,1.5,2.2,2.7},{0.7,1.0,1.5,2.2,2.7}};//{0.6,0.9,1.4,2.1,2.6}};
        if(P_pi < Pbin[ichrg][0]){iPmom = iAllPnTheta;}
        else if(P_pi < Pbin[ichrg][1]){iPmom=0;}//0.7 - 1.0 OR 0.6 - 0.9
        else if(P_pi < Pbin[ichrg][2]){iPmom=1;}//1.0 - 1.5 OR 0.9 - 1.4
        else if(P_pi < Pbin[ichrg][3]){iPmom=2;}//1.5 - 2.2 OR 1.4 - 2.1
        else if(P_pi < Pbin[ichrg][4]){iPmom=3;}//2.2 - 2.7 OR 2.1 - 2.6
        else{iPmom = iAllPnTheta;}

        //Selecting angle between pion and momentum transfer (theta_piq) bin
        if(theta<0.){iThPiQ = -1;}
        else if(theta<2.){iThPiQ = 0;}
        else{iThPiQ = Int_t((theta-2.)/10.);}//2 - 12, 12 - 22, 22 - 32, 32 - 46
        if(iThPiQ > 3 && theta < 46.){iThPiQ = 3;}
        else if(iThPiQ > 3){iThPiQ = iAllPnTheta;}
        switch(iPmom)
        {
            case 1:
                if(iThPiQ > 2){iThPiQ = iAllPnTheta;}
                break;
            case 2:
                if(iThPiQ > 1){iThPiQ = iAllPnTheta;}
                break;
            case 3:
                if(iThPiQ > 0){iThPiQ = iAllPnTheta;}
                break;
        }
    }
    else//QorNot == false (theta_pi)
    {
        //Selecting pi momentum (P_pi) bin
        Double_t Pbin[2][5] = {{0.75,1.1,1.5,2.2,2.7},{0.6,0.9,1.5,1.9,2.4}};
        if(P_pi < Pbin[ichrg][0]){iPmom = -1;}
        else if(P_pi < Pbin[ichrg][1]){iPmom=0;}//0.75 - 1.1 OR 0.6 - 0.9
        else if(P_pi < Pbin[ichrg][2]){iPmom=1;}//1.1  - 1.5 OR 0.9 - 1.5
        else if(P_pi < Pbin[ichrg][3]){iPmom=2;}//1.5  - 2.2 OR 1.5 - 1.9
        else if(P_pi < Pbin[ichrg][4]){iPmom=3;}//2.2  - 2.7 OR 1.9 - 2.4
        else{iPmom = -1;}

        //Selecting angle between pion and momentum transfer (theta_pi) bin
        if(ichrg == 0)//(pi+)
        {
            if(theta<7. || theta>31.){iThPiQ = -1;}
            else{iThPiQ = Int_t((theta-7.)/6.);}//7 - 13, 13 - 19, 19 - 25, 25 - 31
            if(P_pi>Pbin[ichrg][3] && iThPiQ > 1){iThPiQ = -1;}
        }
        else
        {
            Double_t Tbin[6] = {15.,20.,26.,32.,40.,45.};
            if(theta < Tbin[0]){iThPiQ = -1;}
            else if(theta < Tbin[1]){iThPiQ=0;}//15 - 20
            else if(theta < Tbin[2]){iThPiQ=1;}//20 - 26
            else if(theta < Tbin[3]){iThPiQ=2;}//26 - 32
            else if(theta < Tbin[4]){iThPiQ=3;}//32 - 40
            else if(theta < Tbin[5]){iThPiQ=4;}//40 - 45
            else{iThPiQ = -1;}
            switch(iPmom)
            {
                case 0:
                    if(iThPiQ < 2){iThPiQ = -1;}
                    break;
                case 1:
                    if(iThPiQ < 1 || iThPiQ > 3){iThPiQ = -1;}
                    break;
                case 2:
                    if(iThPiQ > 2){iThPiQ = -1;}
                    break;
                case 3:
                    if(iThPiQ > 1){iThPiQ = -1;}
                    break;
            }
            if(iThPiQ==4){iThPiQ=0;}//move "40 - 45" bin to bottom left corner
        }
    }

    //Selecting Q2 (Q2) bin
    if(q2_epi < 0.70){iQ2=-1;}
    else if(q2_epi < 1.00){iQ2=0;}//0.75 - 1.00
    else if(q2_epi < 1.40){iQ2=1;}//1.00 - 1.40
    else if(q2_epi < 1.90){iQ2=2;}//1.40 - 1.90
    else if(q2_epi < 2.50){iQ2=3;}//1.90 - 2.50
    else{iQ2=-1;}

    //Selecting invariant mass (W) bin
    iW = Int_t((W_epi - Wmin)/DeltaW);
    Int_t iWmin = 1, iWmax = Wbini;
    switch(iQ2)
    {
        case 0:
            iWmin = 1; iWmax = 48;// 1.1 - 2.275,
            break;
        case 1:
            iWmin = 1; iWmax = 44;// 1.1 - 2.175
            break;
        case 2:
            iWmin = 1; iWmax = 40;// 1.1 - 2.075
            break;
        case 3:
            iWmin = 1; iWmax = 33;// 1.1 - 1.9
            break;
    }
    if(iW < iWmin || iW >= iWmax) {iW = -1;}
}

std::ifstream& GotoLine(std::ifstream& file, unsigned int num)
{
    file.seekg(std::ios::beg);
    for(unsigned int i=0; i < num; ++i)
    {
        file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    }
    return file;
}
