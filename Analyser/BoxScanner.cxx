#ifndef LARLITE_BOXSCANNER_CXX
#define LARLITE_BOXSCANNER_CXX

#include "BoxScanner.h"

// Function to calculate standard deviation from vector of TDC times
float bTDCstd(std::vector<float> const & TDCvec) {

    float stdTDC(0);
    float meanTDC(0);

    for(size_t f=0; f<TDCvec.size(); ++f) {
        meanTDC += (TDCvec[f]);
    }
    meanTDC /= ((float)TDCvec.size());

    for(size_t k=0; k<TDCvec.size(); ++k)
      stdTDC += (TDCvec[k]-meanTDC)*(TDCvec[k]-meanTDC);
    stdTDC = sqrt( stdTDC / ((float)TDCvec.size()));
    return stdTDC;

}

// Function to calculate interquartile range from vector of TDC times
float bTDCiqr(std::vector<float> TDCvec, int hitNo) {

    std::sort(TDCvec.begin(),TDCvec.end());
    int quart = hitNo/4;
    float iqrTDC = (TDCvec[3*quart] - TDCvec[quart]);
    return iqrTDC;

}

float bMaxVal(std::vector<float> ADCvec){

    // Sort in reverse order
    std::sort(ADCvec.rbegin(),ADCvec.rend());
    float avmax = 0;
    if(ADCvec.size()>20){
      for(int i = 0; i < 20; ++i){
        avmax += ADCvec[i];
      }
      return avmax/20;
    } else return avmax;

}

namespace larlite {

  bool BoxScanner::initialize() {
    // Not sure if events work in same way for reco2D
    _evtN = 0;

    // Initialise TTree with all relevant branches
    smallU_tree = new TTree("smallU_tree","");
    smallU_tree->Branch("evt",&_evtN,"evt/i");
    smallU_tree->Branch("box",&_box,"box/i");
    smallU_tree->Branch("hitNoU",&_hitNoU,"hitNoU/i");
    smallU_tree->Branch("TDCstdU",&_TDCstdU,"TDCstdU/F");
    smallU_tree->Branch("TDCiqrU",&_TDCiqrU,"TDCiqrU/F");
    smallU_tree->Branch("ADCampU",&_ADCampU,"ADCampU/F");
    smallU_tree->Branch("WFintU",&_WFintU,"WFintU/F");
    smallU_tree->Branch("MeanampU",&_MeanampU,"MeanampU/F");
    smallU_tree->Branch("MeanintU",&_MeanintU,"MeanintU/F");
    smallU_tree->Branch("LowDenU",&_LowDenU,"LowDenU/i");
    smallU_tree->Branch("HiDenU",&_HiDenU,"HiDenU/i");
    smallU_tree->Branch("MeanRMSU",&_MeanRMSU,"MeanRMSU/F");
    smallU_tree->Branch("MeanMultU",&_MeanMultU,"MeanMultU/F");
    smallU_tree->Branch("WirestdU",&_WirestdU,"WirestdU/F");
    smallU_tree->Branch("WireiqrU",&_WireiqrU,"WireiqrU/F");
    smallU_tree->SetDirectory(0);

    smallU_bkg_tree = smallU_tree->CloneTree();
    smallU_bkg_tree->SetDirectory(0);
    smallU_bkg_tree->SetName("smallU_bkg_tree");

    // Initialise TTree with all relevant branches
    smallV_tree = new TTree("smallV_tree","");
    smallV_tree->Branch("evt",&_evtN,"evt/i");
    smallV_tree->Branch("box",&_box,"box/i");
    smallV_tree->Branch("hitNoV",&_hitNoV,"hitNoV/i");
    smallV_tree->Branch("TDCstdV",&_TDCstdV,"TDCstdV/F");
    smallV_tree->Branch("TDCiqrV",&_TDCiqrV,"TDCiqrV/F");
    smallV_tree->Branch("ADCampV",&_ADCampV,"ADCampV/F");
    smallV_tree->Branch("WFintV",&_WFintV,"WFintV/F");
    smallV_tree->Branch("MeanampV",&_MeanampV,"MeanampV/F");
    smallV_tree->Branch("MeanintV",&_MeanintV,"MeanintV/F");
    smallV_tree->Branch("LowDenV",&_LowDenV,"LowDenV/i");
    smallV_tree->Branch("HiDenV",&_HiDenV,"HiDenV/i");
    smallV_tree->Branch("MeanRMSV",&_MeanRMSV,"MeanRMSV/F");
    smallV_tree->Branch("MeanMultV",&_MeanMultV,"MeanMultV/F");
    smallV_tree->Branch("WirestdV",&_WirestdV,"WirestdV/F");
    smallV_tree->Branch("WireiqrV",&_WireiqrV,"WireiqrV/F");
    smallV_tree->SetDirectory(0);

    return true;
  }
  
  bool BoxScanner::analyze(storage_manager* storage) {

    _box = 0;
  
    bool showerFlag = 1;

    // Get MCShower information
    auto ev_mcs = storage->get_data<event_mcshower>("mcreco");
    // Check if it exists / are there any showers
    if ( (!ev_mcs) || (!ev_mcs->size())){
    	print (msg::kERROR,__FUNCTION__,"MCShower data product not found!");
        // If you still want to do stuff remove line above
    	showerFlag = 0;
    }

    // Re set shower count
    _showerNo = 0;
    _detShowerNo = 0;

    if (showerFlag){
      // Loop over number of showers? not sure this is how MCShower works
      for (size_t j=0; j < ev_mcs->size(); ++j){
        auto const& shower = (*ev_mcs).at(j);
        // Get energy
        float Energy = shower.Start().E();
        float startX = shower.Start().X()/100;
        float startY = shower.Start().Y()/100;
        float startZ = shower.Start().Z()/100;
        // Count number of showers per event
        if((startX>0&&startX<2.5604)&&(startY>-1.2&&startY<1.2)&&(startZ>0&&startZ<10.5)&&Energy>100){
          _detShowerNo += 1;
        }
        _showerNo += 1;
      }
    }

    // Use storage to get larlite::hit object (gaushit, cccluster, pandoraCosmicKHitRemoval)
    auto hitdata = storage->get_data<event_hit>("gaushit");
    // Display error if hit data not present
    if ( (!hitdata) || (!hitdata->size())){
    	print (msg::kERROR,__FUNCTION__,"Hit data product not found!");
    	return false;
    }

    float lowT=0;float highT=960;int lowUVch=0;int highUVch=240;int lowYch=0;int highYch=345;

    for(int ch(0); ch<19; ++ch){
      for(int t(0); t<19; ++t){

        // Initialise hit counters for event
        _hitNoU = 0;_hitNoV = 0;_hitNoY = 0;
        // Initialize integration counters for event
        _WFintU = 0;_WFintV = 0;_WFintY = 0;
        // Clear TDC vectors
        UTDCvec.clear();VTDCvec.clear();YTDCvec.clear();
        // Clear ADC vectors
        UADCvec.clear();VADCvec.clear();YADCvec.clear();
        // Clear Channel number vectors
        UChvec.clear();VChvec.clear();YChvec.clear();
        // Initialize TDC standard deviation counters
        _TDCstdU = 0;_TDCstdV = 0;_TDCstdY = 0;
        // Initialize TDC interquartile range counters
        _TDCiqrU = 0;_TDCiqrV = 0;_TDCiqrY = 0;
        // Initialize max (or total?) ADC amplitude counters
        _ADCampU = 0;_ADCampV = 0;_ADCampY = 0;
        // Initialize mean adc amplitude counters
        _MeanampU = 0;_MeanampV = 0;_MeanampY = 0;
        // Initailize mean integrated waveform counters
        _MeanintU = 0;_MeanintV = 0;_MeanintY = 0;
        // Initialize hit density between 5-35 ADCs
        _LowDenU = 0;_LowDenV = 0;_LowDenY = 0;
        // Initialize hit density > 35
        _HiDenU = 0;_HiDenV = 0;_HiDenY = 0;
        // Initialize mean RMS counters
        _MeanRMSU = 0;_MeanRMSV = 0;_MeanRMSY = 0;
        // Initialize mean multiplicity counters
        _MeanMultU = 0;_MeanMultV = 0;_MeanMultY = 0;
        // Initialize Wire standard deviation counter
        _WirestdU = 0;_WirestdV = 0;_WirestdY = 0;
        // Initialize Wire interquartile range counter
        _WireiqrU = 0;_WireiqrV = 0;_WireiqrY = 0;

        // Loop over hits in each event
        for (size_t i=0; i < hitdata->size(); ++i){
          auto const& hit = (*hitdata).at(i);
          // Get the channel number
          int chnum = hit.Channel();
          float amp = hit.PeakAmplitude();
          float wint = hit.Integral();
          float tdc = hit.PeakTime();
          float rms = hit.RMS();
          float mult = hit.Multiplicity();
          if(amp>10&&tdc>lowT&&tdc<=highT){

            // Do same for U plane 
            if(chnum>lowUVch&&chnum<=highUVch){ 
              _hitNoU += 1;
              _MeanampU += amp;
              _WFintU += wint;
              _MeanRMSU += rms;
              _MeanMultU += mult;
              UTDCvec.push_back(tdc);
              UADCvec.push_back(amp);
              UChvec.push_back(chnum);
              if(amp<=20)_LowDenU += 1;
              if(amp>20)_HiDenU += 1;
            }

            // Do same for V plane CHECK THIS
            if(chnum>2400+lowUVch&&chnum<=2400+highUVch){
              _hitNoV += 1;
              _MeanampV += amp;
              _WFintV += wint;
              _MeanRMSV += rms;
              _MeanMultV += mult;
              VTDCvec.push_back(tdc);
              VADCvec.push_back(amp);
              VChvec.push_back(chnum);
              if(amp<=20)_LowDenV += 1;
              if(amp>20)_HiDenV += 1;
            }

            // Do same for Y plane CHECK THIS
            if(chnum>4800+lowYch&&chnum<=4800+highYch){
              _hitNoY += 1;
              _MeanampY += amp;
              _WFintY += wint;
              _MeanRMSY += rms;
              _MeanMultY += mult;
              YTDCvec.push_back(tdc);
              YADCvec.push_back(amp);
              YChvec.push_back(chnum);
              if(amp<=20)_LowDenY += 1;
              if(amp>20)_HiDenY += 1;
            }
          }
        }

        _ADCampU = bMaxVal(UADCvec);
        _ADCampV = bMaxVal(VADCvec);
        _ADCampY = bMaxVal(YADCvec);

        // Average the ADC amplitudes
        _MeanampU = _MeanampU/_hitNoU;
        _MeanampV = _MeanampV/_hitNoV;
        _MeanampY = _MeanampY/_hitNoY;

        // Average the integrated waveforms
        _MeanintU = _WFintU/_hitNoU;
        _MeanintV = _WFintV/_hitNoV;
        _MeanintY = _WFintY/_hitNoY;

        _MeanRMSU = _MeanRMSU/_hitNoU;
        _MeanRMSV = _MeanRMSV/_hitNoV;
        _MeanRMSY = _MeanRMSY/_hitNoY;

        _MeanMultU = _MeanMultU/_hitNoU;
        _MeanMultV = _MeanMultV/_hitNoV;
        _MeanMultY = _MeanMultY/_hitNoY;

        // U plane
        if(_hitNoU!=0){
          _TDCstdU = bTDCstd(UTDCvec);
          _TDCiqrU = bTDCiqr(UTDCvec,_hitNoU);
          _WirestdU = bTDCstd(UChvec);
          _WireiqrU = bTDCiqr(UChvec,_hitNoU);
        }
        // V plane
        if(_hitNoV!=0){
          _TDCstdV = bTDCstd(VTDCvec);
          _TDCiqrV = bTDCiqr(VTDCvec,_hitNoV);
          _WirestdV = bTDCstd(VChvec);
          _WireiqrV = bTDCiqr(VChvec,_hitNoV);
        }
        // Y plane
        if(_hitNoY!=0){
          _TDCstdY = bTDCstd(YTDCvec);
          _TDCiqrY = bTDCiqr(YTDCvec,_hitNoY);
          _WirestdY = bTDCstd(YChvec);
          _WireiqrY = bTDCiqr(YChvec,_hitNoY);
        }

        // Fill TTree
        if(_hitNoU!=0) smallU_tree->Fill();
        if(_hitNoU>100/*&&_detShowerNo==0*/) smallU_bkg_tree->Fill();

        _box += 1;
        lowT += 480; highT += 480;
      }
      lowUVch += 120; highUVch += 120; lowYch += 173; highYch += 173;
    }

    _evtN += 1;
  
    return true;
  }

  bool BoxScanner::finalize() {

   // Save histograms to .root file
   if(_fout){
     _fout->cd();
     std::cout << "writing ch tree" << std::endl;
     smallU_tree->Write();
     smallU_bkg_tree->Write();
   }
  
    return true;
  }

}
#endif
