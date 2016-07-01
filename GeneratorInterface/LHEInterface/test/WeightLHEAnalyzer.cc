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



class WeightLHEAnalyzer : public EDAnalyzer {
private: 

  edm::EDGetTokenT<LHEEventProduct> lheEventToken;
  TFile * output;

  TH1F* h_mass;
  TH1F* h_mZ[3];

public:
  explicit WeightLHEAnalyzer( const edm::ParameterSet & cfg ) : 
    src_( cfg.getParameter<InputTag>( "src" )),
    fileName_(cfg.getUntrackedParameter<std::string>("histoutputFile"))
  {
    lheEventToken=consumes<LHEEventProduct>(src_);
  }
  void beginJob(){
    output = new TFile(fileName_.data(), "RECREATE");

    h_mass = new TH1F("h_mass","",300,60,120);
    h_mass->SetXTitle("M_{Z} [GeV]");
    h_mass->Sumw2();
    
    for(int i=0;i<3;i++)
      h_mZ[i]=(TH1F*)h_mass->Clone(Form("h_mZ%d",i));
    
    h_mZ[0]->SetTitle("Including all weights");
    h_mZ[1]->SetTitle("Only positive weights");
    h_mZ[2]->SetTitle("Only negative weights");

  }
  

  ~WeightLHEAnalyzer(){
    delete output;
  }

private:
  void analyze( const Event & iEvent, const EventSetup & iSetup ) {

    Handle<LHEEventProduct> evt;
    iEvent.getByToken( lheEventToken, evt);

    double weight = evt->originalXWGTUP(); 
    
    //    std::cout << "weight = " << weight << std::endl;

    const lhef::HEPEUP hepeup_ = evt->hepeup();

    const int nup_ = hepeup_.NUP; 
    const std::vector<int> idup_ = hepeup_.IDUP;
    const std::vector<lhef::HEPEUP::FiveVector> pup_ = hepeup_.PUP;
    const std::vector<int> istup_ = hepeup_.ISTUP;
    const std::vector<std::pair< int,int > > motup_ = hepeup_.MOTHUP;


    for ( unsigned int icount = 0 ; icount < (unsigned int)nup_; icount++ ) {
      
      int PID    = idup_[icount];
      // int status = istup_[icount];
      // double px = (pup_[icount])[0];
      // double py = (pup_[icount])[1];
      // double pz = (pup_[icount])[2];
      // double e  = (pup_[icount])[3];
      double m =  (pup_[icount])[4];
      // int momIndex = motup_[icount].first-1;
      // int momPID =  momIndex >=0 ? idup_[momIndex]:-1;
 
      if(PID!=23)continue;
      
      h_mZ[0]->Fill(m,weight);
      if(weight>0)h_mZ[1]->Fill(m,weight);
      else h_mZ[2]->Fill(m,weight);

    } // end of loop over particles


  }
      
  void endJob(){
    output->cd();
    for(int i=0;i<3;i++)
      h_mZ[i]->Write();
    output->Write();
    output->Close();
  }   

  void beginRun(edm::Run const& iRun, edm::EventSetup const& es){


  }

  InputTag src_;
  std::string fileName_;
};

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE( WeightLHEAnalyzer );


