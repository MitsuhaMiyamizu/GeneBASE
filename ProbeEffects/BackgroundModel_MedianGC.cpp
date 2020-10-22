/*
 *  BackgroundModel_MedianGC.cpp
 *  ProbeEffects2.0
 *
 *  Created by Karen Kapur on 7/11/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include "BackgroundModel_MedianGC.h"
#include "BackgroundDomain.h"
#include <algorithm>


BackgroundModel_MedianGC::BackgroundModel_MedianGC( ChipBase* i_pChipBase ) : BackgroundModel( i_pChipBase )
{
	int maxGCContent = Parameters::PROBELENGTH;
	m_pBackgroundCoefficients = gsl_vector_alloc( maxGCContent + 1 );
	
	this->FitModel( i_pChipBase );
	
}

BackgroundModel_MedianGC::~BackgroundModel_MedianGC() {
	gsl_vector_free( m_pBackgroundCoefficients );
}

void BackgroundModel_MedianGC::FitModel( ChipBase* i_pChipBase ) 
{
	int maxGCContent = Parameters::PROBELENGTH;
	BackgroundDomain* pBackgroundDomain = BackgroundDomain::GetInstance();
	bool* pIsGCSet = new bool[maxGCContent + 1];
	for (int i = 0; i <= maxGCContent; ++i) {
		pIsGCSet[i] = false;
		vector<double> backgroundVals;
		vector<ProbeBase*>::iterator itBegin = pBackgroundDomain->m_Probes.begin();
		vector<ProbeBase*>::iterator itEnd = pBackgroundDomain->m_Probes.end();
		for (; itBegin != itEnd; ++itBegin ) {
			ProbeBase* pProbeBase = *itBegin;
			if ( pProbeBase->GetGcCount() == i ) {
				backgroundVals.push_back( i_pChipBase->GetUnnormalizedIntensity(pProbeBase));
			}
		}
		// get the middle
		int len = backgroundVals.size();
		if ( len <= 1 ) {
			continue; // Later the value will be set appropriately.  
		}
		int middleIndex = len / 2; 
		partial_sort(backgroundVals.begin(), backgroundVals.begin() + middleIndex, backgroundVals.end());
		double median = backgroundVals[ middleIndex - 1 ];
		gsl_vector_set( m_pBackgroundCoefficients, i, median + .5 );
		pIsGCSet[i] = true;
	}
	
	// Now we must go through those gc-content values which do not have any background probes
	// We set their median intensity equal to the median intensity of the closest bin with background probes
	for (int i = 0; i < maxGCContent; ++i ) {
		if ( pIsGCSet[i] == false ) {
			if ( i < (maxGCContent / 2 ) ) { // We search for values larger than i 
				for ( int j = (i + 1); j < maxGCContent; ++j ) {
					if ( pIsGCSet[j] == true ) {
						double medianVal = gsl_vector_get( m_pBackgroundCoefficients, j);
						gsl_vector_set( m_pBackgroundCoefficients, i, medianVal );
						pIsGCSet[i] = true;
						break;
					}
				}
			} else { // We search for values smaller than i
				for ( int j = (i - 1); j >= 0; --j ) {
					if ( pIsGCSet[j] == true ) {
						double medianVal = gsl_vector_get( m_pBackgroundCoefficients, j);
						gsl_vector_set( m_pBackgroundCoefficients, i, medianVal );
						pIsGCSet[i] = true;
						break;
					}
				}
			}
		}
	}
	delete [] pIsGCSet;
}



const gsl_vector* BackgroundModel_MedianGC::GetBackgroundCoefficients() const
{
	return this->m_pBackgroundCoefficients;
}


void BackgroundModel_MedianGC::OutputCoefficients( std::ofstream& i_Outfile ) const
{
	string SEP = "\t";
	int numParameters = m_pBackgroundCoefficients->size;
	for ( int k = 0; k < numParameters; ++k ) {
		i_Outfile << SEP << gsl_vector_get(m_pBackgroundCoefficients, k);
	}
	i_Outfile << endl;
}

double BackgroundModel_MedianGC::GetBackgroundPrediction( ProbeBase* i_pProbe ) const
{
	int gc = i_pProbe->GetGcCount();
	double background = 0;
	if ( gc >= 0 && gc <= this->m_pBackgroundCoefficients->size ) {
		background = gsl_vector_get( this->m_pBackgroundCoefficients, gc );
	}
	return background;
}







