#ifndef LARLITE_FINDSHOWER_CXX
#define LARLITE_FINDSHOWER_CXX

#include "FindShower.h"

namespace larlite {

  bool FindShower::initialize() {

    _evtN = 0;

    // Initialise TTree with all relevant branches
    _t_ch = new TTree("ch_tree","");
    _t_ch->Branch("TDCvec",&TDCvec);
    _t_ch->Branch("Chvec",&Chvec);
    _t_ch->Branch("Tracks",&Tracks);
    _t_ch->SetDirectory(0);

    return true;
  }
  
  bool FindShower::analyze(storage_manager* storage) {
  
    // Use storage to get larlite::hit object (gaushit, cccluster, pandoraCosmicKHitRemoval)
    auto hitdata = storage->get_data<event_hit>("cccluster");
    if ( (!hitdata) || (!hitdata->size())){
    	print (msg::kERROR,__FUNCTION__,"Hit data product not found!");
    	return false;
    }

    Chvec.clear();
    TDCvec.clear();

    std::vector<std::pair<float,float>> YHits;

    // Get the channel number and hit time for the Y plane
    for (size_t i=0; i < hitdata->size(); ++i){
      auto const& hit = (*hitdata).at(i);
      int chnum = hit.Channel();
      if(chnum>=4800){
        YHits.push_back(std::make_pair(hit.PeakTime(),chnum-4800));
      }
    }

    // Sort hits in time order
    std::sort(YHits.begin(),YHits.end());

    // Vector of bools signifying whether a hit has been checked when finding tracks
    std::vector<bool> isChecked;
    // Vector of temporary tracks
    std::vector<std::pair<float,float>> tempTrack;
    // Track finding tolerance for time
    float T = 5;
    // Track finding tolerance for channel number
    float C = 5;
    // Track finding tolerance for length of track
    float X = 10;

    // Fill the checking vector
    for (size_t j=0; j < YHits.size(); ++j){
      isChecked.push_back(false);
    }

    for (size_t k=0; k < YHits.size()-1; ++k){
      // Skip over points that have been checked
      while(isChecked[k]) ++k;
      // Clear the temporary track
      tempTrack.clear();
      // Get time and channel number
      float time = YHits[k].first;
      float chan = YHits[k].second;
      // Start temporary track
      tempTrack.push_back(std::make_pair(time,chan));
      // Mark hit as checked
      isChecked[k] = true;
      // Look at next point
      int l = k+1;
      float ntime = YHits[l].first;
      float nchan = YHits[l].second;
      // Loop over points until outside of track finding tolerance
      while(ntime-time<T&&l<YHits.size()){
        while(isChecked[l]) ++l;
        ntime = YHits[l].first;
        nchan = YHits[l].second;
        // If hit is within tolerances, increase the temporary track, then look around this new hit
        if(nchan-chan<C&&nchan-chan>-C){
          tempTrack.push_back(std::make_pair(ntime,nchan));
          time = ntime; chan = nchan;
          isChecked[l] = true;
          ++l;
        } else ++l;
      }
      // If the track length is long enough to be a track then add to vector of tracks
      if(tempTrack.size()>X){Tracks.push_back(tempTrack);}
    }
/*
float Dt = 0.3;
    for (size_t k=0; k < YHits.size()-1; ++k){
      while(isChecked[k]) ++k;
      tempTrack.clear();
      float time = YHits[k].first;
      float chan = YHits[k].second;
      tempTrack.push_back(std::make_pair(time,chan));
      isChecked[k] = true;
      int l = k+1;
      float ntime = YHits[l].first;
      float nchan = YHits[l].second;
      //float V1 = std::sqrt((ntime-time)*(ntime-time)+(nchan-chan)*(nchan-chan));
      float theta1 = atan((nchan-chan)/(ntime-time));
      while(ntime-time<T&&l<YHits.size()){
        while(isChecked[l]) ++l;
        ntime = YHits[l].first;
        nchan = YHits[l].second;
        float theta2 = atan((nchan-chan)/(ntime-time));
        if(nchan-chan<C&&nchan-chan>-C&&theta2<theta1+Dt&&theta2>theta1-Dt){
          tempTrack.push_back(std::make_pair(ntime,nchan));
          time = ntime; chan = nchan;
          theta1 = theta2;
          isChecked[l] = true;
          ++l;
        } else ++l;
      }
      if(tempTrack.size()>X){Tracks.push_back(tempTrack);}
    }
*/
    _t_ch->Fill();

    Tracks.clear();

    _evtN += 1;  

    return true;
  }

  bool FindShower::finalize() {

   if(_fout){
     _fout->cd();
     std::cout << "writing ch tree" << std::endl;
     _t_ch->Write();
   }
  
    return true;
  }

}
#endif
