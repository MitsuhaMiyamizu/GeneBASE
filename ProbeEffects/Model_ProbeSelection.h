/*
 *  Model_ProbeSelection.h
 *  ProbeEffects
 *
 *  Created by Karen Kapur on 2/26/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef __Model_ProbeSelection_h_included__
#define __Model_ProbeSelection_h_included__

#include "stdafx.h"
#include "cluster.h"
#include <gsl/gsl_matrix.h>

class Model_ProbeSelection {

public:
	Model_ProbeSelection( gsl_matrix* i_pInten );
	~Model_ProbeSelection();
	gsl_vector* GetExpression();
private:
	double m_corr_cutoff;
	double m_cut_height;
	gsl_vector* m_pExpression;
	
	void GetInitialEstimates( gsl_matrix* i_pInten, gsl_vector* o_pPrevEstimates, 
							  set<int>& o_SelectedProbeIndexes );
	void UpdateExpressionEstimates( gsl_matrix* i_pInten, gsl_vector* o_pUpdatedEstimates, 
									set<int>& i_SelectedProbeIndexes );

	
};



#endif

