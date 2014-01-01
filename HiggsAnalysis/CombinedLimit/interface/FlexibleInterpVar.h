#if ROOT_VERSION_CODE < ROOT_VERSION(5,29,99)
// @(#)root/roostats:$Id: FlexibleInterpVar.h,v 1.1 2011/07/24 20:37:36 gpetrucc Exp $
// Author: Kyle Cranmer, Akira Shibata
/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/
#ifndef ROOSTATS_FLEXIBLEINTERPVAR_COPY
#define ROOSTATS_FLEXIBLEINTERPVAR_COPY

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"

#include <RooRealVar.h>
#include <RooArgList.h>

namespace RooStats{
namespace HistFactory{

  class FlexibleInterpVar : public RooAbsReal {
  public:

    FlexibleInterpVar() ;
    FlexibleInterpVar(const char *name, const char *title,
		      const RooArgList& _paramList, 
		      double nominal, vector<double> low, vector<double> high);

    FlexibleInterpVar(const char *name, const char *title,
		      const RooArgList& _paramList, double nominal, vector<double> low, 
		      vector<double> high,vector<int> code);

    FlexibleInterpVar(const char *name, const char *title);
    FlexibleInterpVar(const FlexibleInterpVar&, const char*);

    void setInterpCode(RooAbsReal& param, int code);
    void setAllInterpCodes(int code);

    void printAllInterpCodes();

    virtual TObject* clone(const char* newname) const { return new FlexibleInterpVar(*this, newname); }
    virtual ~FlexibleInterpVar() ;


  protected:

    RooListProxy _paramList ;
    double _nominal;
    vector<double> _low;
    vector<double> _high;
    vector<int> _interpCode;
    
    TIterator* _paramIter ;  //! do not persist

    Double_t evaluate() const;

    ClassDef(RooStats::HistFactory::FlexibleInterpVar,1) // flexible interpolation
  };
}
}

#endif
#endif
