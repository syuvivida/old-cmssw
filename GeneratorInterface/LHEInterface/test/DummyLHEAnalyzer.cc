#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include <iostream>
#include <string>
#include <TH1F.h>
#include <TFile.h>
#include <TLorentzVector.h>
using namespace std;
using namespace edm;
using namespace lhef;



class DummyLHEAnalyzer : public EDAnalyzer {
private: 

  TFile * output;


  TH1F* h_p;
  TH1F* h_pz;
  TH1F* h_n;

  TH1F* h_HT;
  TH1F* h_Xpt; // new particle 
  TH1F* h_Xpz; 
  TH1F* h_Bpt[2];  // two bosons
  TH1F* h_Bpz[2];
  TH1F* h_Dpt[2][2]; // daughter of two bosons
  TH1F* h_Dpz[2][2]; // daughter of two bosons

  TH1F* h_Xm;
  TH1F* h_Xm_extended;
  TH1F* h_XmT;
  TH1F* h_Bm[2];
  TH1F* h_BmT[2];

  TH1F* h_y;
  TH1F* h_Xy;
  TH1F* h_By[2];
  TH1F* h_Dy[2][2];

  TH1F* h_dR;
  TH1F* h_D_dR[2]; // between the daughters

  TH1F* h_dEta;
  TH1F* h_B_dEta;
  TH1F* h_D_dEta[2];
 

  TH1F* h_cos;
  TH1F* h_cosThetaStar; // Boson scattering angle in the X rest frame
  TH1F* h_cosTheta[2]; //Daugher scattering angle in the Boson 1 rest frame
  TH1F* h_cosPhi; // angle between two decay planes
public:
  explicit DummyLHEAnalyzer( const edm::ParameterSet & cfg ) : 
    src_( cfg.getParameter<InputTag>( "src" )),
    fileName_(cfg.getUntrackedParameter<std::string>("histoutputFile")),
    filterProduction_(cfg.getParameter<int>("filterProduction")),
    filterhh_(cfg.getParameter<bool>("filterhh")),
    nTotal_(0),
    nPass_(0)
  {
    switch(filterProduction_){
    case 0: 
      std::cout << "Filtering mode = keeping all process" << std::endl;
      break;
    case 1:
      std::cout << "Filtering mode = keeping only gluon fusion" << std::endl;
      break;
    case 2:
      std::cout << "Filtering mode = keeping only q-qbar annihilation " << std::endl;
      break;
    default:
      std::cout << "Filtering mode = keeping all process" << std::endl;
      break;
    }
    std::cout << "Filtering hh = " << filterhh_ << std::endl;
  }
  void beginJob(){
    output = new TFile(fileName_.data(), "RECREATE");

    h_p = new TH1F("h_p","",125,0,2500);
    h_p->Sumw2();
    h_pz = new TH1F("h_pz","",150,-3000,3000);
    h_pz->Sumw2();
    h_n = new TH1F("h_n","",6,-0.5,5.5);
    h_n->Sumw2();
    h_n->SetXTitle("parton multiplicity");

    h_y = new TH1F("h_y","",300,-3,3);
    h_y-> Sumw2();
    h_dR = new TH1F("h_dR","",200,0,2);
    h_dR->Sumw2();
    h_dEta = new TH1F("h_dEta","",100,0,10);
    h_dEta->Sumw2();

    h_cos = new TH1F("h_cos","",100,-1,1);
    h_cos->Sumw2();

    h_HT = (TH1F*)h_p->Clone("h_HT");
    h_HT->SetXTitle("H_{T} [GeV]");

    h_Xpt = new TH1F("h_Xpt","",100,0,500);
    h_Xpt-> SetXTitle("p_{T}(X) [GeV]");
    h_Xpt-> Sumw2();

    h_Xpz = (TH1F*)h_pz->Clone("h_Xpz");
    h_Xpz-> SetXTitle("p_{z}(X) [GeV]");

    h_Xm  = new TH1F("h_Xm","",1000,0,5000);
    h_Xm -> SetXTitle("M(X) [GeV]");
    h_Xm -> Sumw2();

    h_Xm_extended  = new TH1F("h_Xm_extended","",25000,0,25000);
    h_Xm_extended -> SetXTitle("M(X) [GeV]");
    h_Xm_extended -> Sumw2();    

    h_XmT  = new TH1F("h_XmT","",1000,0,5000);
    h_XmT -> SetXTitle("M_{T}(X) [GeV]");
    h_XmT -> Sumw2();

    h_Xy = (TH1F*)h_y->Clone("h_Xy");
    h_Xy->SetXTitle("Rapidity of X");

    h_cosThetaStar = (TH1F*)h_cos->Clone("h_cosThetaStar");
    h_cosThetaStar -> SetXTitle("cos#theta^{*}");
    h_cosThetaStar -> SetTitle("Cosine of the scattering angle in the rest frame of the X resonance");

    h_cosPhi       = (TH1F*)h_cos->Clone("h_cosPhi");
    h_cosPhi       -> SetXTitle("cos#Phi");
    h_cosPhi       -> SetTitle("Cosine of the angle between the decay planes");

    h_B_dEta       = (TH1F*)h_dEta->Clone("h_B_dEta");
    h_B_dEta       -> SetXTitle("#Delta #eta between the two bosons");

    for(int i=0; i <2; i++)
      {
	h_Bpt[i] = (TH1F*)h_p->Clone(Form("h_Bpt%d",i));
	h_Bpt[i]->SetTitle(Form("Boson %d",i));
	h_Bpt[i]->SetXTitle(Form("p_{T}(Boson %d) [GeV]",i));

	h_Bpz[i] = (TH1F*)h_pz->Clone(Form("h_Bpz%d",i));
	h_Bpz[i]->SetTitle(Form("Boson %d",i));
	h_Bpz[i]->SetXTitle(Form("p_{z}(Boson %d) [GeV]",i));

	h_Bm[i] = new TH1F(Form("h_Bm%d",i), Form("Boson %d",i), 100,50,150);
	h_Bm[i] -> SetXTitle(Form("M(Boson %d) [GeV]",i));
	h_Bm[i]-> Sumw2();

	h_BmT[i] = new TH1F(Form("h_BmT%d",i), Form("Boson %d",i), 100,50,150);
	h_BmT[i] -> SetXTitle(Form("M_{T}(Boson %d) [GeV]",i));
	h_BmT[i]-> Sumw2();
	
	h_By[i] = (TH1F*)h_y->Clone(Form("h_By%d",i));
	h_By[i] -> SetXTitle(Form("Rapidity of Boson %d",i));

	h_D_dR[i] =(TH1F*)h_dR->Clone(Form("h_D_dR%d",i));
	h_D_dR[i]->SetXTitle(Form("#Delta R between the Daughters of Boson %d",i));

	h_D_dEta[i] = (TH1F*)h_dEta->Clone(Form("h_D_dEta%d",i));
	h_D_dEta[i] -> SetXTitle(Form("#Delta #eta between the Daughters of Boson %d",i));

	h_cosTheta[i] =(TH1F*)h_cos->Clone(Form("h_cosTheta%d",i));
	h_cosTheta[i]->SetXTitle(Form("cos#theta_{%d}",i));
	h_cosTheta[i]->SetTitle(Form("Cosine of the scattering angle in the rest frame of Boson %d",i));

	for(int j=0; j<2; j++)
	  {

	    h_Dpt[i][j] =(TH1F*)h_p->Clone(Form("h_Dpt%d_%d",i,j));
	    h_Dpt[i][j]->SetTitle(Form("Daughter %d of Boson %d",j,i));	    
	    h_Dpt[i][j]->SetXTitle(Form("p_{T}(Daughter %d)[GeV]",j));

	    h_Dpz[i][j] =(TH1F*)h_pz->Clone(Form("h_Dpz%d_%d",i,j));
	    h_Dpz[i][j]->SetTitle(Form("Daughter %d of Boson %d",j,i));
	    h_Dpz[i][j]->SetXTitle(Form("p_{z}(Daughter %d)[GeV]",j));

	    h_Dy[i][j]  = (TH1F*)h_y->Clone(Form("h_Dy%d_%d",i,j));
	    h_Dy[i][j] -> SetXTitle(Form("Rapidity of Daughter %d of Boson %d",j,i));
	  }
    
      }
  }

  ~DummyLHEAnalyzer(){
    delete output;
  }

private:
  void analyze( const Event & iEvent, const EventSetup & iSetup ) {
    
    nTotal_++;
    Handle<LHEEventProduct> evt;
    iEvent.getByLabel( src_, evt );

    double weight =1;// evt->originalXWGTUP(); 
    const lhef::HEPEUP hepeup_ = evt->hepeup();

    const int nup_ = hepeup_.NUP; 
    const std::vector<int> idup_ = hepeup_.IDUP;
    const std::vector<lhef::HEPEUP::FiveVector> pup_ = hepeup_.PUP;
    const std::vector<int> istup_ = hepeup_.ISTUP;
    const std::vector<std::pair< int,int > > motup_ = hepeup_.MOTHUP;

    std::vector<TLorentzVector> l4_vector;
    bool gluonFusion=false;
    unsigned int njet=0;
    double ht=0;
    for ( unsigned int icount = 0 ; icount < (unsigned int)nup_; icount++ ) {
      
      int status = istup_[icount];
      double px = (pup_[icount])[0];
      double py = (pup_[icount])[1];
      double pz = (pup_[icount])[2];
      double e  = (pup_[icount])[3];
      // int momIndex = motup_[icount].first-1;
      // int momPID =  momIndex >=0 ? idup_[momIndex]:-1;
      int PID    = idup_[icount];
      if(status==-1 && PID==21)gluonFusion=true;
      if(status!=1)continue;
      if(abs(PID)==21 || ( abs(PID) >=1 && abs(PID)<=5 ))
	{
	  njet++;
	  ht += sqrt(px*px+py*py);
	}
      l4_vector.push_back(TLorentzVector(px,py,pz,e));

    } // end of loop over particles
    h_HT->Fill(ht);
    h_n->Fill(njet);
    if(l4_vector.size()!=4)return;
    if(!gluonFusion && filterProduction_==1)return;
    if( gluonFusion && filterProduction_==2)return;

    
    TLorentzVector l4_d[2][2];
    TLorentzVector l4_X;
    TLorentzVector l4_B[2];
    TVector3 normal_l3[2];
    for(int i=0;i<2;i++)
      {
	for(int j=0;j<2;j++)
	  {
	    l4_d[i][j]=l4_vector[2*i+j];
	    l4_X += l4_d[i][j];
	    l4_B[i] += l4_d[i][j];
	  }
	normal_l3[i] = l4_d[i][0].Vect().Cross( l4_d[i][1].Vect() );
      }

    if(filterhh_ && fabs(l4_B[0].Eta())>2.5)return;
    if(filterhh_ && fabs(l4_B[1].Eta())>2.5)return;

    float deta_B = fabs(l4_B[0].Eta()-l4_B[1].Eta());
    h_B_dEta->Fill(deta_B);
    if(filterhh_ && deta_B>1.3)return;

    double cosphi= TMath::Cos(normal_l3[0].Angle(normal_l3[1]));
    h_cosPhi->Fill(cosphi,weight);

    TLorentzVector l4boost_B[2];
    TLorentzVector l4boost_d[2][2];

    double costheta[2]={-99999,-99999};

    for(int i=0;i<2;i++)
      {
	l4boost_B[i]=l4_B[i];
	l4boost_B[i].Boost(-l4_X.BoostVector());
	for(int j=0;j<2;j++)
	  {
	    l4boost_d[i][j] = l4_d[i][j];
	    l4boost_d[i][j].Boost(-l4_B[i].BoostVector());
	  }
	
	costheta[i] = TMath::Cos(l4boost_d[i][0].Vect().Angle(l4_B[i].Vect()));
	h_cosTheta[i]->Fill(costheta[i],weight);
      }

    double costhetastar = TMath::Cos(l4boost_B[0].Vect().Angle(l4_X.Vect()));
    h_cosThetaStar->Fill(costhetastar,weight);
    


    // these codes apply only to DM model where the last two daughters are dark matters
    // the first boson has a known mass
    TLorentzVector lt_X;
    TLorentzVector lt_B[2];

    lt_B[0].SetPxPyPzE(
		       l4_B[0].Px(),
		       l4_B[0].Py(),
		       0,
		       sqrt(pow(l4_B[0].E(),2)-
			    pow(l4_B[0].Pz(),2))
		       );

    lt_B[1].SetPxPyPzE(
		       l4_B[1].Px(),
		       l4_B[1].Py(),
		       0,
		       l4_B[1].Pt());
    
    lt_X = lt_B[0]+lt_B[1];
    
    //////////////////////////////////////
    h_Xpt->Fill(l4_X.Pt(),weight);
    h_Xpz->Fill(l4_X.Pz(),weight);
    h_Xm->Fill(l4_X.M(),weight);
    h_Xm_extended->Fill(l4_X.M(),weight);
    h_Xy->Fill(l4_X.Rapidity(),weight);
    h_XmT->Fill(lt_X.M(),weight);

    for(int i=0; i<2; i++)
      {
	h_Bpt[i]->Fill(l4_B[i].Pt(),weight);
	h_Bpz[i]->Fill(l4_B[i].Pz(),weight);
	h_Bm[i]->Fill(l4_B[i].M(),weight);
	h_By[i]->Fill(l4_B[i].Rapidity(),weight);
	h_BmT[i]->Fill(lt_B[i].M(),weight);
	h_D_dR[i]->Fill(l4_d[i][0].DeltaR(l4_d[i][1]),weight);

	float deta=fabs(l4_d[i][0].Eta()-l4_d[i][1].Eta());
	h_D_dEta[i]->Fill(deta);

	for(int j=0; j<2; j++){
	  h_Dpt[i][j]->Fill(l4_d[i][j].Pt(),weight);
	  h_Dpz[i][j]->Fill(l4_d[i][j].Pz(),weight);
	  h_Dy[i][j]->Fill(l4_d[i][j].Rapidity(),weight);
	}
      }
    nPass_++;

  }
      
  void endJob(){

    std::cout << "Total number of processed events = " << nTotal_ << std::endl;
    std::cout << "Total number of passed events = " << nPass_ << std::endl;

    output->cd();
    h_n->Write();
    h_HT->Write();
    h_Xpt->Write();
    h_Xpz->Write();
    h_Xm->Write();
    h_Xm_extended->Write();
    h_XmT->Write();
    h_Xy->Write();
    h_cosThetaStar->Write();
    h_cosPhi->Write();

    h_B_dEta->Write();

    for(int i=0; i<2; i++){
      h_Bpt[i]->Write();
      h_Bpz[i]->Write();
      h_Bm[i]->Write();
      h_BmT[i]->Write();
      h_By[i]->Write();
      h_D_dEta[i]->Write();
      h_D_dR[i]->Write();
      h_cosTheta[i]->Write();

      for(int j=0; j<2; j++){
	h_Dpt[i][j]->Write();
	h_Dpz[i][j]->Write();
	h_Dy[i][j]->Write();
      }

    }

    output->Write();
    output->Close();
  }   

  void beginRun(edm::Run const& iRun, edm::EventSetup const& es){


  }

  InputTag src_;
  std::string fileName_;
  int filterProduction_; // 0: all, 1: gluon only, 2: DY only
  bool filterhh_;
  long int nTotal_;
  long int nPass_;
};

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE( DummyLHEAnalyzer );


