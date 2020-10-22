#ifndef __BackgroundModel_None_h_included__
#define __BackgroundModel_None_h_included__


/*
 *  BackgroundModel_None.h
 *  ProbeEffects2.0
 *
 *  Created by Karen Kapur on 7/11/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */



#include "stdafx.h"
#include "BackgroundModel.h"
#include "ChipBase.h"


class BackgroundModel_None : public BackgroundModel {
	
public:
	BackgroundModel_None( ChipBase* i_pChipBase );
	~BackgroundModel_None();
	
	virtual const gsl_vector* GetBackgroundCoefficients() const;
	virtual void OutputCoefficients( std::ofstream& i_Outfile ) const;
	virtual double GetBackgroundPrediction( ProbeBase* i_pProbe ) const;
	
	
private:
		
	void FitModel( ChipBase* i_pChipBase );
	gsl_vector* m_pBackgroundCoefficients;
	
	
};

#endif
