/*
 *  BackgroundModel_None.cpp
 *  ProbeEffects2.0
 *
 *  Created by Karen Kapur on 7/11/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include "BackgroundModel_None.h"


BackgroundModel_None::BackgroundModel_None( ChipBase* i_pChipBase ) : BackgroundModel( i_pChipBase )
{
	m_pBackgroundCoefficients = gsl_vector_calloc(1); // vector of a single value, set to 0
}


BackgroundModel_None::~BackgroundModel_None() {
	gsl_vector_free( m_pBackgroundCoefficients );
}



const gsl_vector* BackgroundModel_None::GetBackgroundCoefficients() const
{
	return this->m_pBackgroundCoefficients;
}


void BackgroundModel_None::OutputCoefficients( std::ofstream& i_Outfile ) const
{
	i_Outfile << 0 << endl;
}

double BackgroundModel_None::GetBackgroundPrediction( ProbeBase* i_pProbe ) const
{
	return 0;
}


