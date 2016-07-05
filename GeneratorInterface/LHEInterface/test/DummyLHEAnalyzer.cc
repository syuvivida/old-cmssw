#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"
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

  edm::EDGetTokenT<LHEEventProduct>        lheEventToken;
  edm::EDGetTokenT<reco::GenMETCollection> genMETToken_true;
  edm::EDGetTokenT<reco::GenMETCollection> genMETToken_calo;
  edm::EDGetTokenT<reco::GenMETCollection> genMETToken_caloNonPrompt;

  TFile * output;
  TH1F* h_p;
  TH1F* h_pz;
  TH1F* h_n;

  TH1F* h_HT;
  TH1F* h_MET_true;
  TH1F* h_MET_calo;
  TH1F* h_MET_caloNonPrompt;
  TH1F* h_Xpt; // new particle 
  TH1F* h_Xpz; 
  TH1F* h_Bpt[2];  // two bosons
  TH1F* h_Bpz[2];

  TH1F* h_Xm;
  TH1F* h_Xm_extended;
  TH1F* h_XmT;
  TH1F* h_Bm[2];
  TH1F* h_BmT[2];

  TH1F* h_y;
  TH1F* h_Xy;
  TH1F* h_By[2];

  TH1F* h_dR;

  TH1F* h_dEta;
  TH1F* h_B_dEta;
 

  TH1F* h_cos;
  TH1F* h_cosThetaStar; // Boson scattering angle in the X rest frame

public:
  explicit DummyLHEAnalyzer( const edm::ParameterSet & cfg ) : 
    src_( cfg.getParameter<InputTag>( "src" )),
    fileName_(cfg.getUntrackedParameter<std::string>("histoutputFile")),
    nTotal_(0),
    nPass_(0)
  {
    lheEventToken             = consumes<LHEEventProduct>(src_);
    genMETToken_true          = consumes<reco::GenMETCollection>(edm::InputTag("genMetTrue"));
    genMETToken_calo          = consumes<reco::GenMETCollection>(edm::InputTag("genMetCalo"));
    genMETToken_caloNonPrompt = consumes<reco::GenMETCollection>(edm::InputTag("genMetCaloAndNonPrompt"));
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

    h_MET_true = (TH1F*)h_p->Clone("h_MET_true");
    h_MET_true->SetXTitle("#slash{E}_{T}^{truth} [GeV]");

    h_MET_calo = (TH1F*)h_p->Clone("h_MET_calo");
    h_MET_calo->SetXTitle("#slash{E}_{T}^{calo} [GeV]");

    h_MET_caloNonPrompt = (TH1F*)h_p->Clone("h_MET_caloNonPrompt");
    h_MET_caloNonPrompt->SetXTitle("#slash{E}_{T}^{caloNonPrompt} [GeV]");


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

    
      }
  }

  ~DummyLHEAnalyzer(){
    delete output;
  }

private:
  void analyze( const Event & iEvent, const EventSetup & iSetup ) {
    
    nTotal_++;
    
    Handle<reco::GenMETCollection> metHandle_true;
    if(iEvent.getByToken(genMETToken_true, metHandle_true))
      h_MET_true->Fill(metHandle_true.product()->begin()->pt());
    
    Handle<reco::GenMETCollection> metHandle_calo;
    if(iEvent.getByToken(genMETToken_calo, metHandle_calo))
      h_MET_calo->Fill(metHandle_calo.product()->begin()->pt());
  
    Handle<reco::GenMETCollection> metHandle_caloNonPrompt;
    if(iEvent.getByToken(genMETToken_caloNonPrompt, metHandle_caloNonPrompt))
      h_MET_caloNonPrompt->Fill(metHandle_caloNonPrompt.product()->begin()->pt());

    Handle<LHEEventProduct> evt;
    iEvent.getByToken(lheEventToken,evt);

    double weight =1;// evt->originalXWGTUP(); 
    const lhef::HEPEUP hepeup_ = evt->hepeup();

    const int nup_ = hepeup_.NUP; 
    const std::vector<int> idup_ = hepeup_.IDUP;
    const std::vector<lhef::HEPEUP::FiveVector> pup_ = hepeup_.PUP;
    const std::vector<int> istup_ = hepeup_.ISTUP;
    const std::vector<std::pair< int,int > > motup_ = hepeup_.MOTHUP;

  std::vector<TLorentzVector> l4_vector;
  unsigned int njet=0;
  double ht=0;
  for ( unsigned int icount = 0 ; icount < (unsigned int)nup_; icount++ ) {
      
    int status = istup_[icount];
    double px = (pup_[icount])[0];
    double py = (pup_[icount])[1];
    double pz = (pup_[icount])[2];
    double e  = (pup_[icount])[3];
    int PID    = idup_[icount];
    //      std::cout << PID << "\t" << status << "\t" << px << "\t" << py << "\t" << pz << "\t" << e << std::endl;

    if(status!=1)continue;
    if(PID==25 && (motup_[icount].first !=motup_[icount].second))return;
    if(abs(PID)==21 || ( abs(PID) >=1 && abs(PID)<=5 ))
      {
	njet++;
	ht += sqrt(px*px+py*py);
      }
    l4_vector.push_back(TLorentzVector(px,py,pz,e));

  } // end of loop over particles
    // for(unsigned int i=0; i<3; i++)
    //   l4_vector[i].Print();
  h_HT->Fill(ht);
  h_n->Fill(njet);
  if(l4_vector.size()!=3)return;
    
  TLorentzVector l4_X;
  TLorentzVector l4_B[2];
    
  l4_B[0] = l4_vector[0];
  l4_B[1] = l4_vector[1]+l4_vector[2];

  l4_X = l4_B[0]+l4_B[1];


  float deta_B = fabs(l4_B[0].Eta()-l4_B[1].Eta());
  h_B_dEta->Fill(deta_B);

  TLorentzVector l4boost_B[2];



  for(int i=0;i<2;i++)
    {
      l4boost_B[i]=l4_B[i];
      l4boost_B[i].Boost(-l4_X.BoostVector());
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


    }
  nPass_++;

}
      
  void endJob(){
    
    std::cout << "Total number of processed events = " << nTotal_ << std::endl;
    std::cout << "Total number of passed events = " << nPass_ << std::endl;

    output->cd();
    h_n->Write();
    h_HT->Write();
    h_MET_true->Write();
    h_MET_calo->Write();
    h_MET_caloNonPrompt->Write();
    h_Xpt->Write();
    h_Xpz->Write();
    h_Xm->Write();
    h_Xm_extended->Write();
    h_XmT->Write();
    h_Xy->Write();
    h_cosThetaStar->Write();


    h_B_dEta->Write();

    for(int i=0; i<2; i++){
      h_Bpt[i]->Write();
      h_Bpz[i]->Write();
      h_Bm[i]->Write();
      h_BmT[i]->Write();
      h_By[i]->Write();


    }

    output->Write();
    output->Close();
  }   

void beginRun(edm::Run const& iRun, edm::EventSetup const& es){


}

InputTag src_;
std::string fileName_;
long int nTotal_;
long int nPass_;
};

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE( DummyLHEAnalyzer );


