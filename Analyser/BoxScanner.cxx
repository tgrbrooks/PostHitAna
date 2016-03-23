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
    } else return ADCvec[0];

}

namespace larlite {

  bool BoxScanner::initialize() {
    // Inithialize event counter
    _evtN = 0;

    // Initialise TTree for U variables and 960x240 boxes
    smallU_tree = new TTree("smallU_tree","");
    smallU_tree->Branch("evt",&_evtN,"evt/i");
    smallU_tree->Branch("box",&_box,"box/i");
    smallU_tree->Branch("hitNo",&_hitNoU,"hitNo/i");
    smallU_tree->Branch("TDCstd",&_TDCstdU,"TDCstd/F");
    smallU_tree->Branch("TDCiqr",&_TDCiqrU,"TDCiqr/F");
    smallU_tree->Branch("ADCamp",&_ADCampU,"ADCamp/F");
    smallU_tree->Branch("WFint",&_WFintU,"WFint/F");
    smallU_tree->Branch("Meanamp",&_MeanampU,"Meanamp/F");
    smallU_tree->Branch("Meanint",&_MeanintU,"Meanint/F");
    smallU_tree->Branch("LowDen",&_LowDenU,"LowDen/i");
    smallU_tree->Branch("HiDen",&_HiDenU,"HiDen/i");
    smallU_tree->Branch("MeanRMS",&_MeanRMSU,"MeanRMS/F");
    smallU_tree->Branch("MeanMult",&_MeanMultU,"MeanMult/F");
    smallU_tree->Branch("Wirestd",&_WirestdU,"Wirestd/F");
    smallU_tree->Branch("Wireiqr",&_WireiqrU,"Wireiqr/F");
    smallU_tree->Branch("detShowerNo",&_detShowerNo,"detShowerNo/i");
    smallU_tree->SetDirectory(0);

    // Initialise background TTree for U variables and 960x240 boxes
    smallU_bkg_tree = smallU_tree->CloneTree();
    smallU_bkg_tree->SetDirectory(0);
    smallU_bkg_tree->SetName("smallU_bkg_tree");

    // Initialise TTree for U variables and 1200x300 boxes
    largeU_tree = smallU_tree->CloneTree();
    largeU_tree->SetDirectory(0);
    largeU_tree->SetName("largeU_tree");

    // Initialise background TTree for U variables and 1200x300 boxes
    largeU_bkg_tree = smallU_tree->CloneTree();
    largeU_bkg_tree->SetDirectory(0);
    largeU_bkg_tree->SetName("largeU_bkg_tree");

    // Initialise TTree for V variables and 960x240 boxes
    smallV_tree = new TTree("smallV_tree","");
    smallV_tree->Branch("evt",&_evtN,"evt/i");
    smallV_tree->Branch("box",&_box,"box/i");
    smallV_tree->Branch("hitNo",&_hitNoV,"hitNo/i");
    smallV_tree->Branch("TDCstd",&_TDCstdV,"TDCstd/F");
    smallV_tree->Branch("TDCiqr",&_TDCiqrV,"TDCiqr/F");
    smallV_tree->Branch("ADCamp",&_ADCampV,"ADCamp/F");
    smallV_tree->Branch("WFint",&_WFintV,"WFint/F");
    smallV_tree->Branch("Meanamp",&_MeanampV,"Meanamp/F");
    smallV_tree->Branch("Meanint",&_MeanintV,"Meanint/F");
    smallV_tree->Branch("LowDen",&_LowDenV,"LowDen/i");
    smallV_tree->Branch("HiDen",&_HiDenV,"HiDen/i");
    smallV_tree->Branch("MeanRMS",&_MeanRMSV,"MeanRMS/F");
    smallV_tree->Branch("MeanMult",&_MeanMultV,"MeanMult/F");
    smallV_tree->Branch("Wirestd",&_WirestdV,"Wirestd/F");
    smallV_tree->Branch("Wireiqr",&_WireiqrV,"Wireiqr/F");
    smallV_tree->SetDirectory(0);

    // Initialise background TTree for V variables and 960x240 boxes
    smallV_bkg_tree = smallV_tree->CloneTree();
    smallV_bkg_tree->SetDirectory(0);
    smallV_bkg_tree->SetName("smallV_bkg_tree");

    // Initialise TTree for V variables and 1200x300 boxes
    largeV_tree = smallV_tree->CloneTree();
    largeV_tree->SetDirectory(0);
    largeV_tree->SetName("largeV_tree");

    // Initialise background TTree for V variables and 1200x300 boxes
    largeV_bkg_tree = smallV_tree->CloneTree();
    largeV_bkg_tree->SetDirectory(0);
    largeV_bkg_tree->SetName("largeV_bkg_tree");

    // Initialise TTree for Y variables and 960x346 boxes
    smallY_tree = new TTree("smallY_tree","");
    smallY_tree->Branch("evt",&_evtN,"evt/i");
    smallY_tree->Branch("box",&_box,"box/i");
    smallY_tree->Branch("hitNo",&_hitNoY,"hitNo/i");
    smallY_tree->Branch("TDCstd",&_TDCstdY,"TDCstd/F");
    smallY_tree->Branch("TDCiqr",&_TDCiqrY,"TDCiqr/F");
    smallY_tree->Branch("ADCamp",&_ADCampY,"ADCamp/F");
    smallY_tree->Branch("WFint",&_WFintY,"WFint/F");
    smallY_tree->Branch("Meanamp",&_MeanampY,"Meanamp/F");
    smallY_tree->Branch("Meanint",&_MeanintY,"Meanint/F");
    smallY_tree->Branch("LowDen",&_LowDenY,"LowDen/i");
    smallY_tree->Branch("HiDen",&_HiDenY,"HiDen/i");
    smallY_tree->Branch("MeanRMS",&_MeanRMSY,"MeanRMS/F");
    smallY_tree->Branch("MeanMult",&_MeanMultY,"MeanMult/F");
    smallY_tree->Branch("Wirestd",&_WirestdY,"Wirestd/F");
    smallY_tree->Branch("Wireiqr",&_WireiqrY,"Wireiqr/F");
    smallY_tree->SetDirectory(0);

    // Initialise background TTree for Y variables and 960x346 boxes
    smallY_bkg_tree = smallY_tree->CloneTree();
    smallY_bkg_tree->SetDirectory(0);
    smallY_bkg_tree->SetName("smallY_bkg_tree");

    // Initialise TTree for Y variables and 1200x432 boxes
    largeY_tree = smallY_tree->CloneTree();
    largeY_tree->SetDirectory(0);
    largeY_tree->SetName("largeY_tree");

    // Initialise background TTree for Y variables and 1200x432 boxes
    largeY_bkg_tree = smallY_tree->CloneTree();
    largeY_bkg_tree->SetDirectory(0);
    largeY_bkg_tree->SetName("largeY_bkg_tree");

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
        // Count number of showers per event with E<100MeV that start in AV
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

    // Set upper and lower limits on box sizes
    float lowT=0;float highT=960;int lowUVch=0;int highUVch=240;int lowYch=0;int highYch=345;

    // Loop over 19 values for time and channel number (boxes increment by 1/2)
    for(int ch(0); ch<19; ++ch){
      lowT = 0; highT = 960;
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
          // Get the relevent data
          int chnum = hit.Channel();
          float amp = hit.PeakAmplitude();
          float wint = hit.Integral();
          float tdc = hit.PeakTime();
          float rms = hit.RMS();
          float mult = hit.Multiplicity();
          // Change ADC tolerance from 5 (gaushit) to 10 and set time range
          if(amp>10&&tdc>lowT&&tdc<=highT){

            // Calculate variables in U channel range
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

            // Calculate variables in V channel range
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

            // Calculate variables in Y channel range
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

        // Get "average" maximum ADC amplitude
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

        // Fill TTrees (background if no high energy showers in AV)
        if(_hitNoU!=0) smallU_tree->Fill();
        if(_hitNoU!=0&&_detShowerNo==0) smallU_bkg_tree->Fill();
        if(_hitNoV!=0) smallV_tree->Fill();
        if(_hitNoV!=0&&_detShowerNo==0) smallV_bkg_tree->Fill();
        if(_hitNoY!=0) smallY_tree->Fill();
        if(_hitNoY!=0&&_detShowerNo==0) smallY_bkg_tree->Fill();

        _box += 1;
        // Increment time limits by 1/2 side of box
        lowT += 480; highT += 480;
      }
      // Increment channel limits by 1/2 side of box
      lowUVch += 120; highUVch += 120; lowYch += 173; highYch += 173;
    }
/*
    // Reset time and channel limits for larger boxes
    lowT=0; highT=1200; lowUVch=0; highUVch=300; lowYch=0; highYch=432;

    // Repeat the process for larger boxes
    for(int ch(0); ch<15; ++ch){
      lowT = 0; highT = 1200;
      for(int t(0); t<15; ++t){

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
          // Get the relevent data
          int chnum = hit.Channel();
          float amp = hit.PeakAmplitude();
          float wint = hit.Integral();
          float tdc = hit.PeakTime();
          float rms = hit.RMS();
          float mult = hit.Multiplicity();
          if(amp>10&&tdc>lowT&&tdc<=highT){

            // Calculate U variables 
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

            // Do same for V plane
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

            // Do same for Y plane
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

        // Fill TTrees
        if(_hitNoU!=0) largeU_tree->Fill();
        if(_hitNoU!=0&&_detShowerNo==0) largeU_bkg_tree->Fill();
        if(_hitNoV!=0) largeV_tree->Fill();
        if(_hitNoV!=0&&_detShowerNo==0) largeV_bkg_tree->Fill();
        if(_hitNoY!=0) largeY_tree->Fill();
        if(_hitNoY!=0&&_detShowerNo==0) largeY_bkg_tree->Fill();

        _box += 1;
        // Increment time limits
        lowT += 600; highT += 600;
      }
      // Increment channel limits
      lowUVch += 150; highUVch += 150; lowYch += 216; highYch += 216;
    }
*/
    _evtN += 1;
  
    return true;
  }

  bool BoxScanner::finalize() {

   // Save TTrees to .root file
   if(_fout){
     _fout->cd();
     std::cout << "writing ch tree" << std::endl;
     smallU_tree->Write();
     smallU_bkg_tree->Write();
/*     smallV_tree->Write();
     smallV_bkg_tree->Write();
     smallY_tree->Write();
     smallY_bkg_tree->Write();
     largeU_tree->Write();
     largeU_bkg_tree->Write();
     largeV_tree->Write();
     largeV_bkg_tree->Write();
     largeY_tree->Write();
     largeY_bkg_tree->Write();*/
   }
  
    return true;
  }

}
#endif
