#ifndef __BackgroundModel_MedianGC_h_included__
#define __BackgroundModel_MedianGC_h_included__


/*
 *  BackgroundModel_MedianGC.h
 *  ProbeEffects2.0
 *
 *  Created by Karen Kapur on 7/11/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */


#include "stdafx.h"
#include "BackgroundModel.h"
#include "ChipBase.h"


class BackgroundModel_MedianGC : public BackgroundModel {
	
public:
	BackgroundModel_MedianGC( ChipBase* i_pChipBase );
	~BackgroundModel_MedianGC();
	
	virtual const gsl_vector* GetBackgroundCoefficients() const;
	virtual void OutputCoefficients( std::ofstream& i_Outfile ) const;
	virtual double GetBackgroundPrediction( ProbeBase* i_pProbe ) const;

	
private:
	
	void FitModel( ChipBase* i_pChipBase );
	gsl_vector* m_pBackgroundCoefficients;
	
	
};

#endif
