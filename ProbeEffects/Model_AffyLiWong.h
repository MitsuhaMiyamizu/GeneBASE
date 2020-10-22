/*
 *  Model_AffyLiWong.h
 *  ProbeEffects
 *
 *  Created by Karen Kapur on 5/16/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef __Model_AffyLiWong_h_included__
#define __Model_AffyLiWong_h_included__

#include "stdafx.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

class Model_AffyLiWong {

private:
	
	gsl_vector* m_pTheta;
	gsl_vector* m_pPhi;
	int m_NumArrays;
	int m_NumProbes;
	
	// Model fitting parameters
	const double m_NormalArrayQuantile;
	const double m_NormalResidQuantile;
	const int m_LargeThreshold;
	const double m_LargeVariation;
	const double m_OutlierFraction;
	const double m_Delta;
	const int m_Maxit;
	const int m_OuterMaxit;
	
	// Helper functions
	void FitLiWong( const gsl_matrix* i_pDiff );

public:
		
		// ---Model_AffyLiWong---  
		// Input:  I by J matrix of probe intensities, where I = num samples, J = num probes
		// Output:  Stores phi and theta in member varibles m_pPhi and m_pTheta
		// These variables are accessed by calling GetPhi() and GetTheta()
		// Details:  
		// Fits the Li-Wong model described in Li and Wong (2001).  
		// Code is copied from the fit.li.wong function in the R affy package
		// We use the default parameter values from the function.  
		// ------
		Model_AffyLiWong( const gsl_matrix* i_pDiff );
	
	~Model_AffyLiWong(void);
	gsl_vector* GetTheta(); // Vector of expression estimates
	gsl_vector* GetPhi(); // Vector of probe affinities
	
};




#endif

