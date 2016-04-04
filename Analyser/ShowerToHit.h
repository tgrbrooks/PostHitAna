/**
 * \file ShowerToHit.h
 *
 * \ingroup Analyser
 * 
 * \brief Class def header for a class ShowerToHit
 *
 * @author tom
 */

/** \addtogroup Analyser

    @{*/

#ifndef LARLITE_SHOWERTOHIT_H
#define LARLITE_SHOWERTOHIT_H

#include "Analysis/ana_base.h"
#include "DataFormat/hit.h"
#include "DataFormat/mcshower.h"
#include <iostream>
#include <vector>
#include <utility>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <string>
#include <fstream>
#include <math.h>

namespace larlite {

  class ShowerToHit : public ana_base{
  
  public:

    /// Default constructor
    ShowerToHit(){ _name="ShowerToHit"; _fout=0;}

    /// Default destructor
    virtual ~ShowerToHit(){}

    /** IMPLEMENT in ShowerToHit.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in ShowerToHit.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in ShowerToHit.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

  protected:

    int _evtN;

    std::vector<float> YTDCvec;

    std::vector<float> YChvec;

    std::vector<std::pair<float,float>> Times;

    std::vector<std::pair<int,int>> Channels;

    TTree* _t_ch;
    
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
