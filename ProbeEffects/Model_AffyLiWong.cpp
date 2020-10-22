/*
 *  Model_AffyLiWong.cpp
 *  ProbeEffects
 *
 *  Created by Karen Kapur on 5/16/08.
 *
 */

#include "Model_AffyLiWong.h"
#include "my_math.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include "DB_VECTOR.h"



Model_AffyLiWong::Model_AffyLiWong( const gsl_matrix* i_pDiff ) : m_NormalArrayQuantile(0.5), m_NormalResidQuantile(0.9), 
m_LargeThreshold(3), m_LargeVariation(0.8), m_OutlierFraction(0.14), m_Delta(1e-06), m_Maxit(50), m_OuterMaxit(50) 
{
	// Set member variables
	m_NumArrays = i_pDiff->size1;
	m_NumProbes = i_pDiff->size2;
	m_pTheta = gsl_vector_alloc( m_NumArrays );
	m_pPhi = gsl_vector_alloc( m_NumProbes );
	
	this->FitLiWong( i_pDiff );
	
}

void Model_AffyLiWong::FitLiWong( const gsl_matrix* i_pDiff ) {
	
	// Setup
	// must check that rows and columns are in correct orientation
	int II = i_pDiff->size1;
	int J = i_pDiff->size2;
	int cI = II;
	int cJ = J;
	
	// Helper variables
	gsl_vector* pOnesII = gsl_vector_alloc( II );
	gsl_vector_set_all( pOnesII, 1 );
	gsl_vector* pOnesJ = gsl_vector_alloc( J );
	gsl_vector_set_all( pOnesJ, 1 );
	
	gsl_vector* pThetaOutliersOld = gsl_vector_calloc( II ); // initialize all to false
	gsl_vector* pPhiOutliersOld = gsl_vector_calloc( J );
	gsl_matrix* pSingleOutliersOld = gsl_matrix_calloc( II, J );
	gsl_vector* pThetaOutliers = gsl_vector_calloc( II );
	gsl_vector* pPhiOutliers = gsl_vector_calloc( J );
	gsl_matrix* pSingleOutliers = gsl_matrix_calloc( II, J );
	
	bool flag1 = true;
	bool flag2 = true;
	gsl_matrix* pDataMatrix = gsl_matrix_alloc( II, J);
	gsl_matrix* pDataMatrixHat = gsl_matrix_alloc( II, J );
	gsl_matrix* pResid = gsl_matrix_alloc( II, J );
	gsl_matrix_memcpy( pDataMatrix, i_pDiff );
	gsl_vector* pPhi = gsl_vector_alloc( J );
	gsl_vector* pTheta = gsl_vector_alloc( II );
	
	// remove.outliers is assumed to be true
	bool changeTheta = true;
	bool changePhi = true;
	bool changeSingle = true;
	int outerIter = 0;
	while ( (flag1 == true) && (flag2 == true) && (changeTheta + changePhi + changeSingle > 0) && (outerIter < m_OuterMaxit) ) {
		++outerIter;
		if ( (outerIter % 3 == 0 && changeTheta == true) || (outerIter % 3 == 1 && changePhi == true) ) {
			// Set phi and theta
			gsl_blas_dgemv (CblasTrans, 1.0, pDataMatrix, pOnesII, 0.0, pPhi);  // phi = colMeans( pDataMatrix )
			gsl_vector_scale( pPhi, (1.0/static_cast<double>(II) ) ); 
			double c = 0;
			for ( int j = 0; j < J; ++j ) {
				if ( gsl_vector_get( pPhiOutliers, j ) == 0 ) {
					c += gsl_pow_2( gsl_vector_get( pPhi, j ) );
				}
			}
			if ( c > 0 ) { // otherwise c = 0
				c = sqrt( cJ / c );  // c <- sqrt(cJ/sum(phi[!phi.outliers]^2))
			} 
			gsl_vector_scale (pPhi, c); // phi <- c * phi			
			for ( int i = 0; i < II; ++i ) { // theta <- (data.matrix[, !phi.outliers, drop = FALSE] %*% phi[!phi.outliers, drop = FALSE])/cJ
				gsl_vector_view currRow = gsl_matrix_row (pDataMatrix, i);
				double valTheta = InnerProduct(&currRow.vector, pPhi, J, pPhiOutliers) / cJ;
				gsl_vector_set( pTheta, i, valTheta );
			}
			int iter = 0;
			bool change = true;
			gsl_vector* pThetaOld = gsl_vector_calloc(II); 
			while ( (change == true) && (iter < m_Maxit) ) {
				++iter;
				// Set phi and theta given the current outlier values
				for ( int j = 0; j < J; ++j ) { // phi <- t(data.matrix[!theta.outliers, , drop = FALSE]) %*% theta[!theta.outliers, drop = FALSE]
					gsl_vector_view currCol = gsl_matrix_column (pDataMatrix, j);
					double valPhi = InnerProduct(&currCol.vector, pTheta, II, pThetaOutliers);
					gsl_vector_set( pPhi, j, valPhi );
				}
				double cInner = 0; // c <- sqrt(cJ/sum(phi[!phi.outliers, drop = FALSE]^2))
				for ( int j = 0; j < J; ++j ) {
					if ( gsl_vector_get( pPhiOutliers, j ) == 0 ) {
						cInner += gsl_pow_2( gsl_vector_get( pPhi, j) );
					}
				}
				if ( cInner > 0 ) { // otherwise cInner = 0
					cInner = sqrt( cJ / cInner );
				}
				gsl_vector_scale( pPhi, cInner ); // phi <- c * phi
				for ( int i = 0; i < II; ++i ) { // theta <- (data.matrix[, !phi.outliers, drop = FALSE] %*% phi[!phi.outliers, drop = FALSE])/cJ
					gsl_vector_view currRow = gsl_matrix_row (pDataMatrix, i);
					double valTheta = InnerProduct(&currRow.vector, pPhi, J, pPhiOutliers) / cJ;
					gsl_vector_set( pTheta, i, valTheta );
				}
				change = !Same(pTheta, pThetaOld, II, pThetaOutliers, m_Delta); // change <- max(abs(theta[!theta.outliers] - theta.old[!theta.outliers]))
				gsl_vector_memcpy( pThetaOld, pTheta );
			}
			gsl_vector_free( pThetaOld );
			if (iter >= m_Maxit) {
				flag1 = false;
			}
			int numPhiNeg = 0;
			int numPhiIncluded = 0;
			for ( int j = 0; j < J; ++j ) {
				if ( gsl_vector_get( pPhiOutliers, j ) == 0 ) {
					if ( gsl_vector_get( pPhi, j ) < 0 ) {
						++numPhiNeg;
					} 
					++numPhiIncluded;
				}
			}
			if ( static_cast<double>(numPhiNeg) / static_cast<double>(numPhiIncluded) > 0.5 ) { // if (mean(phi[!phi.outliers] < 0) > 0.5) ...				
				gsl_vector_scale( pPhi, -1 );
				gsl_vector_scale( pTheta, -1 );
			}
			
			gsl_matrix_set_zero( pDataMatrixHat );
			gsl_blas_dger (1.0, pTheta, pPhi, pDataMatrixHat); //data.matrixhat <- outer(theta, phi)
			gsl_matrix_memcpy( pResid, pDataMatrix );
			gsl_matrix_sub(pResid, pDataMatrixHat); //resid <- data.matrix - data.matrixhat
		}
		if ( outerIter % 3 == 1) {
			for ( int i = 0; i < II; ++i ) { // use pSingleOutliers to sort abs(residuals)
				for ( int j = 0; j < J; ++j ) { 
					double val = gsl_matrix_get( pResid, i, j );
					if ( val < 0 ) { // take absolute value
						val = -val; 
					}
					gsl_matrix_set( pSingleOutliers, i, j, val );
				}
			}			
			double quantVal = partial_sort(pSingleOutliers->data, II*J, m_NormalResidQuantile); 		
			// go through the residual matrix
			for ( int i = 0; i < II; ++i ) {
				for ( int j = 0; j < J; ++j ) {
					if ( gsl_matrix_get( pResid, i, j) > quantVal * m_LargeThreshold ) {
						gsl_matrix_set( pSingleOutliers, i, j, 1 ); // outlier
					} else {
						gsl_matrix_set( pSingleOutliers, i, j, 0 ); // not an outlier
					}
				}
			}
			// if an array has too many single outliers we don't count them
			for ( int i = 0; i < II; ++i ) {
				gsl_vector_view rowView = gsl_matrix_row (pSingleOutliers, i);
				double rowOutlierSum;
				gsl_blas_ddot (&rowView.vector, pOnesJ, &rowOutlierSum);
				if ( rowOutlierSum > ( m_OutlierFraction * cJ) ) {
					gsl_vector_set_all( &rowView.vector, 0 );
				}
			}
			// if a probe has too many single outliers we don't count them
			for ( int j = 0; j < J; ++j ) {
				gsl_vector_view colView = gsl_matrix_column (pSingleOutliers, j);
				double colOutlierSum;
				gsl_blas_ddot (&colView.vector, pOnesII, &colOutlierSum);
				if ( colOutlierSum > ( m_OutlierFraction * cI) ) {
					gsl_vector_set_all( &colView.vector, 0 );
				}
			}
			// go through each element of the matrix, updating pDataMatrix appropriately
			for ( int i = 0; i < II; ++i ) {
				for ( int j = 0; j < J; ++j ) {
					if ( gsl_matrix_get( pSingleOutliers, i, j ) == 0 ) {
						gsl_matrix_set( pDataMatrix, i, j, gsl_matrix_get(i_pDiff, i, j) );
					} else {
						gsl_matrix_set( pDataMatrix, i, j, gsl_matrix_get(pDataMatrixHat, i, j) );
					}
				}
			}
			gsl_matrix_sub( pSingleOutliersOld, pSingleOutliers);
			changeSingle = !gsl_matrix_isnull (pSingleOutliersOld);
			gsl_matrix_memcpy( pSingleOutliersOld, pSingleOutliers );
		} else {
			gsl_vector* pSigmaTheta = gsl_vector_calloc( II );
			gsl_vector* pSigmaPhi = gsl_vector_calloc( J );
			for ( int j = 0; j < J; ++j ) {
				if ( gsl_vector_get( pPhiOutliers, j ) == 0 ) {
					for ( int i = 0; i < II; ++i ) {
						double currval = gsl_vector_get( pSigmaTheta, i );
						gsl_vector_set( pSigmaTheta, i, currval + gsl_pow_2(gsl_matrix_get(pResid, i, j)) );
					}
				}
			}
			for (int i = 0; i < II; ++i ) {
				if ( gsl_vector_get( pThetaOutliers, i ) == 0 ) {
					for ( int j = 0; j < J; ++j ) {
						double currval = gsl_vector_get( pSigmaPhi, j );
						gsl_vector_set( pSigmaPhi, j, currval + gsl_pow_2(gsl_matrix_get(pResid, i, j)) );						
					}
				}
			}
			for ( int i = 0; i < II; ++i ) { //	sigma.theta <- sqrt( rowSums(resid[, !phi.outliers, drop = FALSE]^2)/(cJ - 1) )
				double currval = gsl_vector_get( pSigmaTheta, i );
				gsl_vector_set( pSigmaTheta, i, (sqrt(currval / (cJ - 1) ) ) );
			}
			for ( int j = 0; j < J; ++j ) { // sigma.phi <- sqrt(colSums(resid[!theta.outliers, , drop = FALSE]^2)/(cI - 1))
				double currval = gsl_vector_get( pSigmaPhi, j );
				gsl_vector_set( pSigmaPhi, j, (sqrt(currval / (cI - 1) ) ) );
			}
			
			if (outerIter % 3 == 2) {
				// Changes order of pSigmaTheta
				gsl_vector* pVecII = gsl_vector_alloc( II );
				gsl_vector_memcpy(pVecII, pSigmaTheta);
				double quantVal = partial_sort(pVecII->data, II, m_NormalArrayQuantile); 
				double sumThetaSquared;
				gsl_vector_memcpy( pVecII, pTheta );
				gsl_vector_mul (pVecII, pTheta);
				gsl_blas_ddot (pVecII, pOnesII, &sumThetaSquared);
				gsl_vector_free(pVecII);
				cI = 0;
				for ( int i = 0; i < II; ++i ) { // theta.outliers <- sigma.theta > large.threshold * quantile(sigma.theta, normal.array.quantile) 
					// | theta^2/sum(theta^2) > large.variation
					bool cond1 = gsl_vector_get( pSigmaTheta, i ) > ( m_LargeThreshold * quantVal);
					bool cond2 = ( gsl_pow_2(gsl_vector_get( pTheta, i )) / sumThetaSquared) > m_LargeVariation;
					if ( cond1 == true || cond2 == true  ) {
						gsl_vector_set( pThetaOutliers, i, 1 );
					} else {
						gsl_vector_set( pThetaOutliers, i, 0 );
						++cI; // cI <- sum(!theta.outliers)
					}
				}
				if (cI < 3) {
                    flag2 = false;
				}
				for ( int i = 0; i < II; ++i ) { // single.outliers[theta.outliers, ] <- rep(FALSE, J)
					if ( gsl_vector_get( pThetaOutliers, i ) > 0 ) {
						gsl_vector_view rowVec = gsl_matrix_row (pSingleOutliers, i);
						gsl_vector_set_all( &rowVec.vector, 0 );
					}
					for ( int j = 0; j < J; ++j ) {
						if ( gsl_matrix_get( pSingleOutliers, i, j) > 0 ) { // data.matrix[single.outliers] <- data.matrixhat[single.outliers]
							gsl_matrix_set( pDataMatrix, i, j, gsl_matrix_get(pDataMatrixHat, i, j) );
						} else {
							gsl_matrix_set( pDataMatrix, i, j, gsl_matrix_get(i_pDiff, i, j) ); // data.matrix[!single.outliers] <- original.data.matrix[!single.outliers]
						}
					}
				}
				
				gsl_vector_sub( pThetaOutliersOld, pThetaOutliers);
				changeTheta = !gsl_vector_isnull (pThetaOutliersOld); // change.theta <- sum(abs(theta.outliers.old - theta.outliers))
				gsl_vector_memcpy( pThetaOutliersOld, pThetaOutliers ); // theta.outliers.old <- theta.outliers
				gsl_matrix* pMatIIJ = gsl_matrix_alloc( II, J );
				gsl_matrix_memcpy( pMatIIJ, pSingleOutliersOld );
				gsl_matrix_sub( pMatIIJ, pSingleOutliers);
				changeSingle = !gsl_matrix_isnull (pMatIIJ); // change.single <- sum(abs(single.outliers.old - single.outliers))
				gsl_matrix_free( pMatIIJ );					
			} else {
				
				// Changes order of pSigmaPhi
				gsl_vector* pVecJ = gsl_vector_alloc( J );
				gsl_vector_memcpy(pVecJ, pSigmaPhi);
				double quantVal = partial_sort(pVecJ->data, J, m_NormalArrayQuantile); 
				double sumPhiSquared;
				gsl_vector_memcpy( pVecJ, pPhi );
				gsl_vector_mul (pVecJ, pPhi);
				gsl_blas_ddot (pVecJ, pOnesJ, &sumPhiSquared);
				gsl_vector_free(pVecJ);
				cJ = 0;
				for ( int j = 0; j < J; ++j ) { // phi.outliers <- sigma.phi > large.threshold * quantile(sigma.phi, normal.array.quantile) 
												// | phi^2/sum(phi^2) > large.variation | phi < 0
					bool cond1 = gsl_vector_get( pSigmaPhi, j ) > ( m_LargeThreshold * quantVal);
					bool cond2 = gsl_pow_2(gsl_vector_get( pPhi, j )) / sumPhiSquared > m_LargeVariation;
					bool cond3 = gsl_vector_get( pPhi, j ) < 0;
					if ( cond1 == true || cond2 == true || cond3 == true  ) {
						gsl_vector_set( pPhiOutliers, j, 1 );
					} else {
						gsl_vector_set( pPhiOutliers, j, 0 );
						++cJ; // cJ <- sum(!phi.outliers)
					}
				}
				if (cJ < 3) {
                    flag2 = false;
				}
				
				for ( int j = 0; j < J; ++j ) { // single.outliers[, phi.outliers ] <- rep(FALSE, J)
					if ( gsl_vector_get( pPhiOutliers, j ) > 0 ) {
						gsl_vector_view colVec = gsl_matrix_column (pSingleOutliers, j);
						gsl_vector_set_all( &colVec.vector, 0 );
					}
					for ( int i = 0; i < II; ++i ) {
						if ( gsl_matrix_get( pSingleOutliers, i, j) > 0 ) { // data.matrix[single.outliers] <- data.matrixhat[single.outliers]
							gsl_matrix_set( pDataMatrix, i, j, gsl_matrix_get(pDataMatrixHat, i, j) );
						} else {
							gsl_matrix_set( pDataMatrix, i, j, gsl_matrix_get(i_pDiff, i, j) ); // data.matrix[!single.outliers] <- original.data.matrix[!single.outliers]
						}
					}
				}
				
				gsl_vector_sub( pPhiOutliersOld, pPhiOutliers);
				changePhi = !gsl_vector_isnull (pPhiOutliersOld); // change.phi <- sum(abs(phi.outliers.old - phi.outliers))
				gsl_vector_memcpy( pPhiOutliersOld, pPhiOutliers ); // phi.outliers.old <- phi.outliers
				gsl_matrix* pMatIIJ = gsl_matrix_alloc( II, J );
				gsl_matrix_memcpy( pMatIIJ, pSingleOutliersOld );
				gsl_matrix_sub( pMatIIJ, pSingleOutliers);
				changeSingle = !gsl_matrix_isnull (pMatIIJ); // change.single <- sum(abs(single.outliers.old - single.outliers))
				gsl_matrix_free( pMatIIJ );	
			}
			gsl_vector_free( pSigmaTheta );
			gsl_vector_free( pSigmaPhi );
		}
	}
	if ( outerIter >= m_OuterMaxit) {
		flag2 == false; // no convergence
	}

	// store important stuff
	gsl_vector_memcpy( m_pTheta, pTheta );
	gsl_vector_memcpy( m_pPhi, pPhi );
	// Omit storing outliers, sigma, sigma.theta, sigma.phi
	
	
	// Free memory
	gsl_vector_free( pThetaOutliersOld );
	gsl_vector_free( pPhiOutliersOld );
	gsl_matrix_free( pSingleOutliersOld );
	gsl_vector_free( pThetaOutliers );
	gsl_vector_free( pPhiOutliers );
	gsl_matrix_free( pSingleOutliers );
	gsl_matrix_free( pDataMatrix );
	gsl_matrix_free( pDataMatrixHat );
	gsl_matrix_free( pResid );
		
	gsl_vector_free( pOnesII );
	gsl_vector_free( pOnesJ );

	gsl_vector_free( pPhi );
	gsl_vector_free( pTheta );  
	
}
	
Model_AffyLiWong::~Model_AffyLiWong(void) {
	// free memory
	gsl_vector_free( m_pTheta );
	gsl_vector_free( m_pPhi );
	
}


gsl_vector* Model_AffyLiWong::GetTheta()
{
	return this->m_pTheta; 
}
gsl_vector* Model_AffyLiWong::GetPhi()
{
	return this->m_pPhi; 
}


	