#ifndef __BackgroundModel_h_included__
#define __BackgroundModel_h_included__

//***********
// BackgroundModel
// Base class for BackgroundModel_MAT, BackgroundModel_MedianGC, BackgroundModel_None
// Any background model can be defined to extend this class provided 
// it gives implementation for all pure virtual functions.  
// Also, can add information to BackgroundDomain::OutputBackgroundModelHeader.  
//***********

#include "stdafx.h"
#include "ProbeBase.h"

class ChipBase;

class BackgroundModel {

public:
	BackgroundModel( ChipBase* i_pChipBase );
	virtual ~BackgroundModel();
	
	virtual const gsl_vector* GetBackgroundCoefficients() const = 0; // this makes the class abstract
	virtual void OutputCoefficients( std::ofstream& i_Outfile ) const = 0;
	virtual double GetBackgroundPrediction( ProbeBase* i_pProbe ) const = 0;

	
private:
	
};



#endif

