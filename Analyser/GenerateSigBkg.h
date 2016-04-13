/**
 * \file GenerateSigBkg.h
 *
 * \ingroup Analyser
 * 
 * \brief Class def header for a class GenerateSigBkg
 *
 * @author tom
 */

/** \addtogroup Analyser

    @{*/

#ifndef LARLITE_GENERATESIGBKG_H
#define LARLITE_GENERATESIGBKG_H

#include "Analysis/ana_base.h"
#include <iostream>
#include <vector>
#include <utility>
//#include <cmath>
#include "DataFormat/hit.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/mctrack.h"
#include "DataFormat/mctruth.h"
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <string>
#include <fstream>
#include <math.h>

namespace larlite {
  /**
     \class GenerateSigBkg
     User custom analysis class made by SHELL_USER_NAME
   */
  class GenerateSigBkg : public ana_base{
  
  public:

    /// Default constructor
    GenerateSigBkg(){ _name="GenerateSigBkg"; _fout=0;}

    /// Default destructor
    virtual ~GenerateSigBkg(){}

    /** IMPLEMENT in GenerateSigBkg.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in GenerateSigBkg.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in GenerateSigBkg.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

  protected:

    int _evtN;

    int _showerNo, _detShowerNo;

    int _hitNoY;

    float _TDCstdY;

    float _TDCiqrY;

    float _ADCampY;

    float _MeanampY;

    float _WFintY;

    float _MeanintY;

    int _LowDenY;

    int _HiDenY;

    float _MeanRMSY;

    float _MeanMultY;

    float _WirestdY;

    float _WireiqrY;

    std::vector<float> TDCvec, YTDCvec;
   
    std::vector<float> ADCvec, YADCvec;

    std::vector<float> Chvec, YChvec;

    TTree* sig_tree; TTree* bkg_tree; 
    
  };
}
#endif

//**************************************************************************
// 
// For Analysis framework documentation, read Manual.pdf here:
//
// http://microboone-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=3183
//
//**************************************************************************

/** @} */ // end of doxygen group 
