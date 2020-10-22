/*
 *  BackgroundModel_MAT.cpp
 *  ProbeEffects2.0
 *
 *  Created by Karen Kapur on 7/9/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include "BackgroundModel.h"
#include "BackgroundModel_MAT.h"
#include "BackgroundDomain.h"
#include "ChipBase.h"
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_blas.h>


BackgroundModel_MAT::BackgroundModel_MAT( ChipBase* i_pChipBase ) : BackgroundModel( i_pChipBase )
{

	Parameters* pParameters = Parameters::GetInstance();
	int numMATParams = pParameters->GetNumParameters_MAT();
	m_pMATCoefficients = gsl_vector_alloc( numMATParams );
	m_RSquaredLogScale = 0;
	m_RSquaredOriginalScale = 0;
	m_SigmaHat = 1;
	m_pCov = gsl_matrix_alloc( numMATParams, numMATParams );
	
	
	this->FitModel( i_pChipBase );
	
}

BackgroundModel_MAT::~BackgroundModel_MAT() {
	gsl_vector_free( m_pMATCoefficients );
	gsl_matrix_free( m_pCov );
}

void BackgroundModel_MAT::FitModel( ChipBase* i_pChipBase ) {
	Parameters* pParameters = Parameters::GetInstance();
	int numMATParams = pParameters->GetNumParameters_MAT();
	BackgroundDomain* pBackgroundDomain = BackgroundDomain::GetInstance();
	
	int numProbes = pBackgroundDomain->m_Probes.size();
	gsl_matrix*	pDesignMatrix = gsl_matrix_calloc( numProbes, numMATParams );
	gsl_vector*	pLogPM = gsl_vector_calloc( numProbes );
	gsl_vector* pPM = gsl_vector_calloc( numProbes );
	if ( pDesignMatrix == NULL || pLogPM == NULL || pPM == NULL ) {
		Parameters::LOGFILE << "Memory error allocating design matrix or y vector." << endl;
		exit(1);
	}
	
	vector<ProbeBase*>::iterator itProbeBegin = pBackgroundDomain->m_Probes.begin();
	vector<ProbeBase*>::iterator itProbeEnd = pBackgroundDomain->m_Probes.end();
	int selectedProbeCtr = 0;
	for (; itProbeBegin != itProbeEnd; ++itProbeBegin ) {
		ProbeBase* pProbe = *itProbeBegin;
		gsl_vector_view designRowView = gsl_matrix_row (pDesignMatrix, selectedProbeCtr);
		pProbe->SetMATDesignVector( (&designRowView.vector) );
		double pm = static_cast<double>( i_pChipBase->GetUnnormalizedIntensity( pProbe ) );
		double minVal = 1;
		double logPm = log((max(pm, minVal)) );
		gsl_vector_set( pLogPM, selectedProbeCtr, logPm);
		gsl_vector_set( pPM, selectedProbeCtr, pm );
		++selectedProbeCtr;
	}
	gsl_multifit_linear_workspace* pWorkspace = gsl_multifit_linear_alloc (numProbes, numMATParams);
	
	gsl_vector* pLogPredicted = gsl_vector_alloc( numProbes );
	gsl_vector* pPredicted = gsl_vector_alloc( numProbes );
	if ( pWorkspace == NULL || pLogPredicted == NULL || pPredicted == NULL ) {
		Parameters::LOGFILE << "Memory error allocating linear workspace, coefficient vector or predicted vector." << endl;
		exit(1);
	}
	
	double ssy_LogScale = (numProbes - 1) * gsl_stats_variance( pLogPM->data, 1, numProbes );
	double ssy_OriginalScale = (numProbes - 1) * gsl_stats_variance( pPM->data, 1, numProbes );
	double rss;
	gsl_multifit_linear (pDesignMatrix, pLogPM, m_pMATCoefficients, m_pCov, (&rss), pWorkspace); 
	
	gsl_blas_dgemv (CblasNoTrans, 1.0, pDesignMatrix, m_pMATCoefficients, 0.0, pLogPredicted);
	for ( int i = 0; i < numProbes; ++i ) {
		double val = gsl_vector_get( pLogPredicted, i );
		val = exp( val );
		gsl_vector_set( pPredicted, i, val );
	}
	gsl_vector_sub (pLogPredicted, pLogPM);
	double rss_LogScale;
	gsl_blas_ddot (pLogPredicted, pLogPredicted, &rss_LogScale);
	gsl_vector_sub (pPredicted, pPM);
	double rss_OriginalScale;
	gsl_blas_ddot (pPredicted, pPredicted, &rss_OriginalScale);
	m_SigmaHat = sqrt( rss_LogScale / (numProbes - numMATParams) );
	
	// Free memory
	gsl_vector_free( pLogPredicted );
	gsl_vector_free( pPredicted );
	gsl_multifit_linear_free( pWorkspace );
	gsl_matrix_free( pDesignMatrix );
	gsl_vector_free( pLogPM );
	gsl_vector_free( pPM );
	
	m_RSquaredLogScale = 1 - ( rss_LogScale / ssy_LogScale ); 
	m_RSquaredOriginalScale = 1 - ( rss_OriginalScale / ssy_OriginalScale );
	
}


const gsl_vector* BackgroundModel_MAT::GetBackgroundCoefficients() const
{
	return this->m_pMATCoefficients;
}

void BackgroundModel_MAT::OutputCoefficients( std::ofstream& i_Outfile ) const
{
	string SEP = "\t";
	int numParameters = m_pMATCoefficients->size;
	for ( int k = 0; k < numParameters; ++k ) {
		i_Outfile << SEP << gsl_vector_get(m_pMATCoefficients, k);
	}
	i_Outfile << SEP << this->GetRSquaredLogScale()
		<< SEP << this->GetRSquaredOriginalScale() 
		<< SEP << this->GetSigmaHat() << endl;
	
}

double BackgroundModel_MAT::GetBackgroundPrediction( ProbeBase* i_pProbe ) const
{
	gsl_vector* pDesignVec = gsl_vector_alloc( this->m_pMATCoefficients->size );
	i_pProbe->SetMATDesignVector( pDesignVec );
	double logBackground;
	gsl_blas_ddot (pDesignVec, m_pMATCoefficients, &logBackground );
	double background = exp( logBackground ); 
	gsl_vector_free( pDesignVec );
	return background;
}

double BackgroundModel_MAT::GetBackgroundPredictionStandardError( ProbeBase* i_pProbe ) const
{
	double predError;
	gsl_vector* pDesignVec = gsl_vector_alloc( this->m_pMATCoefficients->size );
	i_pProbe->SetMATDesignVector( pDesignVec );
	gsl_vector* pTempVec = gsl_vector_alloc( this->m_pMATCoefficients->size );
	gsl_blas_dgemv (CblasNoTrans, 1, m_pCov, pDesignVec, 0, pTempVec);
	gsl_blas_ddot (pDesignVec, pTempVec, (&predError));
	predError += gsl_pow_2(this->m_SigmaHat); // Prediction error so include sigmaHatSq
	predError = sqrt(predError);
	gsl_vector_free( pTempVec );
	gsl_vector_free( pDesignVec );
	return predError;
}


double BackgroundModel_MAT::GetRSquaredLogScale() const
{
	return this->m_RSquaredLogScale;
}
double BackgroundModel_MAT::GetRSquaredOriginalScale() const
{
	return this->m_RSquaredOriginalScale;
}

double BackgroundModel_MAT::GetSigmaHat() const
{
	return this->m_SigmaHat;
}





