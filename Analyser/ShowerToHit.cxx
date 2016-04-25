#ifndef LARLITE_SHOWERTOHIT_CXX
#define LARLITE_SHOWERTOHIT_CXX

#include "ShowerToHit.h"

namespace larlite {

  bool ShowerToHit::initialize() {

    _evtN = 0;

    // Initialise TTree with all relevant branches
    _t_ch = new TTree("ch_tree","");
    _t_ch->Branch("evt",&_evtN,"evt/i");
    _t_ch->Branch("YTDCvec",&YTDCvec); 
    _t_ch->Branch("YChvec",&YChvec);
    _t_ch->Branch("Times",&Times); 
    _t_ch->Branch("Channels",&Channels);
    _t_ch->Branch("StartX",&_startx);
    _t_ch->SetDirectory(0);


    return true;
  }
  
  bool ShowerToHit::analyze(storage_manager* storage) {
  
bool showerFlag = 1;

    // Get MCShower information
    auto ev_mcs = storage->get_data<event_mcshower>("mcreco");
    // Check if it exists / are there any showers
    if ( (!ev_mcs) || (!ev_mcs->size())){
    	print (msg::kERROR,__FUNCTION__,"MCShower data product not found!");
        // If you still want to do stuff remove line above
    	showerFlag = 0;
    }

    Times.clear();
    Channels.clear();
    _startx.clear();

    float conv = 1e-09/0.5e-06;//ns/tick

    if (showerFlag){
      for (size_t j=0; j < ev_mcs->size(); ++j){
        auto const& shower = (*ev_mcs).at(j);
        float Energy = shower.Start().E();
        float StartX = shower.Start().X(); //cm
        float StartY = shower.Start().Y();
        float StartZ = shower.Start().Z();
        float EndX = shower.End().X();
        float EndY = shower.End().Y();
        float EndZ = shower.End().Z();
        float ST = shower.Start().T(); //ns
        float ET = shower.End().T();
        if(((StartX>0&&StartX<256.04)&&(StartY>-116.25&&StartY<116.25)&&(StartZ>0&&StartZ<1036.8)&&Energy>75)||((EndX>0&&EndX<256.04)&&(EndY>-116.25&&EndY<116.25)&&(EndZ>0&&EndZ<1036.8)&&Energy>75)){
          float sst = (ST*conv)+3200;
          float ssdt = StartX/(1300*100*0.5e-06); //vd*100*tick WHICH SIDE IS WIRE PLANE
          float set = (ET*conv)+3200;
          float sedt = EndX/(1300*100*0.5e-06); //vd*100*tick 
          _startx.push_back(StartX);
          if(StartZ<EndZ){
            int StartCh = (int)(StartZ*10/3) - 20;
            int EndCh = (int)(EndZ*10/3) + 20;
            Channels.push_back(std::make_pair(StartCh, EndCh));
          } else {
            int StartCh = (int)(StartZ*10/3) + 20;
            int EndCh = (int)(EndZ*10/3) - 20;
            Channels.push_back(std::make_pair(StartCh, EndCh));
          }
          if(sst+ssdt<set+sedt){
            float StartT = sst + ssdt - 60;
            float EndT = set + sedt + 60;
            Times.push_back(std::make_pair(StartT, EndT)); //TDC ticks
          } else {
            float StartT = sst + ssdt + 60;
            float EndT = set + sedt - 60;
            Times.push_back(std::make_pair(StartT, EndT)); //TDC ticks
          }
        }
      }
    }

    // Use storage to get larlite::hit object (gaushit, cccluster, pandoraCosmicKHitRemoval)
    auto hitdata = storage->get_data<event_hit>("gaushit");
    // Display error if hit data not present
    if ( (!hitdata) || (!hitdata->size())){
    	print (msg::kERROR,__FUNCTION__,"Hit data product not found!");
    	return false;
    }

    YTDCvec.clear();YChvec.clear();

    // Loop over hits in each event
    for (size_t i=0; i < hitdata->size(); ++i){
      auto const& hit = (*hitdata).at(i);
      // Get the channel number
      int chnum = hit.Channel();
      float tdc = hit.PeakTime();
      float amp = hit.PeakAmplitude();
      if(amp>5){
        // Y plane
        if(chnum>=4800&&chnum<8256){
          YTDCvec.push_back(tdc);
          YChvec.push_back(chnum);
        }
      }
    }

    _t_ch->Fill();

    _evtN += 1;
  
    return true;
  }

  bool ShowerToHit::finalize() {

   // Save histograms to .root file
   if(_fout){
     _fout->cd();
     std::cout << "writing ch tree" << std::endl;
     _t_ch->Write();
   }
  
    return true;
  }

}
#endif
