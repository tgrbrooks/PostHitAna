/**
 * \file BoxScanner.h
 *
 * \ingroup Analyser
 * 
 * \brief Class def header for a class BoxScanner
 *
 * @author tom
 */

/** \addtogroup Analyser

    @{*/

#ifndef LARLITE_BOXSCANNER_H
#define LARLITE_BOXSCANNER_H

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
     \class BoxScanner
     User custom analysis class made by SHELL_USER_NAME
   */
  class BoxScanner : public ana_base{
  
  public:

    /// Default constructor
    BoxScanner(){ _name="BoxScanner"; _fout=0;}

    /// Default destructor
    virtual ~BoxScanner(){}

    /** IMPLEMENT in BoxScanner.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in BoxScanner.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in BoxScanner.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

  protected:

    int _evtN, _box;

    int _showerNo, _detShowerNo;

    int _hitNoU, _hitNoV, _hitNoY;

    float _TDCstdU, _TDCstdV, _TDCstdY;

    float _TDCiqrU, _TDCiqrV, _TDCiqrY;

    float _ADCampU, _ADCampV, _ADCampY;

    float _MeanampU, _MeanampV, _MeanampY;

    float _WFintU, _WFintV, _WFintY;

    float _MeanintU, _MeanintV, _MeanintY;

    int _LowDenU, _LowDenV, _LowDenY;

    int _HiDenU, _HiDenV, _HiDenY;

    float _MeanRMSU, _MeanRMSV, _MeanRMSY;

    float _MeanMultU, _MeanMultV, _MeanMultY;

    float _WirestdU, _WirestdV, _WirestdY;

    float _WireiqrU, _WireiqrV, _WireiqrY;

    std::vector<float> TDCvec, UTDCvec, VTDCvec, YTDCvec;
   
    std::vector<float> ADCvec, UADCvec, VADCvec, YADCvec;

    std::vector<float> Chvec, UChvec, VChvec, YChvec;

    TTree* smallU_tree; TTree* largeU_tree; TTree* largeV_tree; TTree* smallY_tree; TTree* largeY_tree; TTree* smallV_tree;
    TTree* smallU_bkg_tree; TTree* largeU_bkg_tree; TTree* smallV_bkg_tree; TTree* largeV_bkg_tree; TTree* smallY_bkg_tree; TTree* largeY_bkg_tree;
    
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
