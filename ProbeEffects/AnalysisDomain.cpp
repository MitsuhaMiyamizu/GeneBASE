
/*
 *  AnalysisDomain.cpp
 *  ProbeEffects2.0
 *
 *  Created by Karen Kapur on 7/9/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include <gsl/gsl_vector.h>
#include <gsl/gsl_statistics_double.h>
#include "AnalysisDomain.h"
#include "TranscriptClusterDomain.h"
#include "BackgroundDomain.h"
#include "chip_data.h"
#include "my_math.h"
#include "Model_AffyLiWong.h"
#include "Model_ProbeSelection.h"
#include "PresenceAbsenceModel_MAT.h"
#include "PresenceAbsenceModel_MedianGC.h"


bool AnalysisDomain::m_InstanceFlag = false;
AnalysisDomain* AnalysisDomain::m_pTheOne = NULL;

AnalysisDomain::AnalysisDomain( ) : m_pMeanQuantiles(NULL)
{
}

AnalysisDomain::~AnalysisDomain(void)
{
}

AnalysisDomain* AnalysisDomain::GetInstance( ) {
	if ( !m_InstanceFlag )
	{
		m_pTheOne = new AnalysisDomain( );
		m_InstanceFlag = true;
		return m_pTheOne;
	}
	else {
		return m_pTheOne;
	}
}

bool AnalysisDomain::ModelBackground()
{
	Parameters* pParameters = Parameters::GetInstance();
	TranscriptClusterDomain* pTranscriptClusterDomain = TranscriptClusterDomain::GetInstance();
	BackgroundDomain* pBackgroundDomain = BackgroundDomain::GetInstance();
	
	
	// Output background model header
	ofstream outfile;
	outfile.open( pParameters->GetOutputModelFit().c_str() );
	pBackgroundDomain->OutputBackgroundModelHeader( outfile );
		
	// If quantile normalization
	string normalizationMethod = pParameters->GetNormalizationMethod();
	if ( normalizationMethod == "quantile" ) {
		this->GetMeanQuantiles(); // Reads in each array
								// Calculates the mean quantile values.  
	}
	
	pair< set<string>::iterator, set<string>::iterator > itChipPair = pParameters->GetCelFilesIterators();
	set<string>::iterator itChipBegin = itChipPair.first;
	set<string>::iterator itChipEnd = itChipPair.second;
	for (; itChipBegin != itChipEnd; ++itChipBegin) {
		string filename = pParameters->GetArrayCelFolder() + (*itChipBegin);
		ifstream datafile;
		datafile.open(filename.c_str());
		if (datafile.is_open())
		{
			string::size_type pos = filename.rfind( Parameters::SEPARATOR );
			int len = static_cast<int>(filename.length());
			string smallname = filename.substr( (pos + 1), (len - 1) );
			//!!!
			// CHIP_DATA is for exon arrays 
			CHIP_DATA* pChip = new CHIP_DATA( smallname, pParameters->GetNumCells(), pParameters->GetCellDim() );
			// Read in the chip's data
			bool readOk = pChip->ReadBinaryCel(filename.c_str());
			if (readOk == false) {
				delete pChip;
				continue;
			}
			// If quantile normalization, do it
			// ***Only quantile normalization is done before background correction***
			if ( normalizationMethod == "quantile" ) {
				this->ComputeQuantileNormalizedIntensities( pChip );
			}
			// Model background
			bool modelBkgdOk = pChip->ModelBackground();
			if ( modelBkgdOk == false ) {
				delete pChip;
				continue;
			}
			// Output the parameter estimates and r-squared
			string SEP = "\t";
			outfile << pChip->GetId();
			const BackgroundModel* pBackgroundModel = pChip->GetBackgroundModel();
			pBackgroundModel->OutputCoefficients(outfile );
			// If scalar normalization, do it
			pChip->SetProbeNormalizationConstant();
			// Output background-corrected probe intensities
			if ( pParameters->OutputAllBkgdCorrectNormProbes() == true ) {	
				pChip->OutputAllBkgdCorrectNormProbes();
			}
		} else {
			Parameters::LOGFILE << "ERROR:  " << "Cannot open CEL file " << filename << endl;
		}
	}
	this->DeleteMeanQuantiles(); // Frees space storing m_pMeanQuantiles
	
	return true;
}


bool AnalysisDomain::GetMeanQuantiles()
{
	// Allocate m_pMeanQuantiles
	// Should be deleted after normalizing!
	Parameters* pParameters = Parameters::GetInstance();
	int numcells = pParameters->GetNumCells();
	m_pMeanQuantiles = new double[ numcells ];
	for (int i = 0; i < numcells; ++i) {
		m_pMeanQuantiles[i] = 0;
	}
	
	// for each chip
	// add raw intensities to m_pMeanQuantiles
	// delete the chip	
	int numArrays = 0;
	pair< set<string>::iterator, set<string>::iterator > itChipPair = pParameters->GetCelFilesIterators();
	set<string>::iterator itChipBegin = itChipPair.first;
	set<string>::iterator itChipEnd = itChipPair.second;
	for (; itChipBegin != itChipEnd; ++itChipBegin) {
		string filename = pParameters->GetArrayCelFolder() + (*itChipBegin);
		ifstream datafile;
		datafile.open(filename.c_str());
		if (datafile.is_open()) {
			string::size_type pos = filename.rfind( Parameters::SEPARATOR );
			int len = static_cast<int>(filename.length());
			string smallname = filename.substr( (pos + 1), (len - 1) );
			//!!!
			// CHIP_DATA is for exon arrays 
			CHIP_DATA* pChip = new CHIP_DATA( smallname, pParameters->GetNumCells(), pParameters->GetCellDim() );
			// Read in the chip's data
			bool readOk = pChip->ReadBinaryCel(filename.c_str());
			if (readOk == false) {
				delete pChip;
				continue;
			}
			unsigned short* pSortedVec = new unsigned short[numcells];
			for (int i = 0; i < numcells; ++i) {
				unsigned short val = pChip->m_pInten[i];
				pSortedVec[i] = val; // would like to simply copy the values of pChip->m_pInten to pSortedVec
			}
			std::sort(pSortedVec, pSortedVec + numcells);
			for (int i = 0; i < numcells; ++i) {
				m_pMeanQuantiles[i] += pSortedVec[i];
			}
			++numArrays;
			delete [] pSortedVec;
			delete pChip;
		} else {
			Parameters::LOGFILE << "ERROR:  " << "Unable to open data file " << filename << endl;
		}
	}
	
	if ( numArrays > 0 ) {
		for (int i = 0; i < numcells; ++i) {
			m_pMeanQuantiles[i] = m_pMeanQuantiles[i] / static_cast<double>(numArrays); // double / double
		}
	} else {
		Parameters::LOGFILE << "ERROR:  " << "No CEL files have been read.  Unable to normalize.  "  << endl;
		return false;
	} 
	return true;
}

void AnalysisDomain::ComputeQuantileNormalizedIntensities( ChipBase* i_pChipBase )
{
	int numcells = i_pChipBase->GetNumCells();
	// Compute the normalized intensities
	INT_POINT* seq = new INT_POINT[numcells];
	for (int i = 0; i < numcells; ++i) {
		seq[i].x = i_pChipBase->m_pInten[i];
		seq[i].y = i;
	}
	
	qsort((void *)seq, numcells, sizeof(INT_POINT), compare_int_point);
	for (int i = 0; i < numcells; ++i) {
		int index = seq[i].y;
		i_pChipBase->m_pInten[index] = static_cast<unsigned int>( m_pMeanQuantiles[i] ); // assign the ith quantile to the value with rank i
																					 // the location of the value with rank i is stored
																					 // in seq.  
	}
	delete [] seq;
	
}

void AnalysisDomain::DeleteMeanQuantiles() // Deletes m_pMeanQuantiles
{
	delete [] m_pMeanQuantiles;
}

bool AnalysisDomain::SummarizeExpression()
{
	// Parameters for filtering probes matching off-targets
	Parameters* pParameters = Parameters::GetInstance();
	TranscriptClusterDomain* pTranscriptClusterDomain = TranscriptClusterDomain::GetInstance();
	const double correlationFilterMinSd = pParameters->GetCorrelationFilterMinSd();
	const double minOffTargetCorrelation2Deconvolve = pParameters->GetCorrelationFilterMinCorrelation();
	const bool outputMasked = pParameters->GetCorrelationFilterOutputMasked();
	const int maxEditDistance = pParameters->GetCorrelationFilterEditDistance();
	string summaryMethod = pParameters->GetSummaryMethod(); // liwong, selection, or correlation_filter
	hash_map<unsigned long, pair<bool, gsl_vector*> > deconvolvedEstimates;

	// If quantile normalization
	string normalizationMethod = pParameters->GetNormalizationMethod();
	if ( normalizationMethod == "quantile" ) {
		this->GetMeanQuantiles(); // Reads in each array
								  // Calculates the mean quantile values.  
	}
	
	// For each chip, read intensities, somehow store background-corrected, normalized intensities
	// SummarizeExpression does not need raw probe intensities. 
	vector<CHIP_DATA*> chips;
	pair< set<string>::iterator, set<string>::iterator > itChipPair = pParameters->GetCelFilesIterators();
	set<string>::iterator itChipBegin = itChipPair.first;
	set<string>::iterator itChipEnd = itChipPair.second;
	for (; itChipBegin != itChipEnd; ++itChipBegin) {
		string filename = pParameters->GetArrayCelFolder() + (*itChipBegin);
		ifstream datafile;
		datafile.open(filename.c_str());
		if (datafile.is_open())
		{
			string::size_type pos = filename.rfind( Parameters::SEPARATOR );
			int len = static_cast<int>(filename.length());
			string smallname = filename.substr( (pos + 1), (len - 1) );
			//!!!
			// CHIP_DATA is for exon arrays 
			CHIP_DATA* pChip = new CHIP_DATA( smallname, pParameters->GetNumCells(), pParameters->GetCellDim() );
			// Read in the chip's data
			bool readOk = pChip->ReadBinaryCel(filename.c_str());
			if (readOk == false) {
				delete pChip;
				continue;
			}
			// If quantile normalization, do it
			// ***Only quantile normalization is done before background correction***
			if ( normalizationMethod == "quantile" ) {
				this->ComputeQuantileNormalizedIntensities( pChip );
			}
			// Model background
			bool modelBkgdOk = pChip->ModelBackground();
			if ( modelBkgdOk == false ) {
				delete pChip;
				continue;
			}
			// If scalar normalization, do it
			pChip->SetProbeNormalizationConstant();
			chips.push_back(pChip);
		}
	}
	int numSamples = chips.size();
			
	//********
	// Create output file
	//********
	ofstream outfile;
	outfile.open( pParameters->GetOutputModelFit().c_str() );
	const string SEP = "\t"; 
	outfile << "TranscriptCluster";
	vector<CHIP_DATA*>::iterator itChBegin = chips.begin();
	vector<CHIP_DATA*>::iterator itChEnd = chips.end();
	for (; itChBegin != itChEnd; ++itChBegin ) {
		outfile << SEP << (*itChBegin)->GetId();
	}
	outfile << endl;
	
		
		
	
	// store initial estimates of expression (for li-wong or selection, these are the final estimates)
	pair< hash_map<unsigned long, TranscriptCluster* >::iterator, hash_map<unsigned long, TranscriptCluster*>::iterator > itTcPair = pTranscriptClusterDomain->GetTranscriptClusterIterators();
	hash_map<unsigned long, TranscriptCluster*>::iterator itTcBegin = itTcPair.first;
	hash_map<unsigned long, TranscriptCluster*>::iterator itTcEnd = itTcPair.second;
	for (; itTcBegin != itTcEnd; ++itTcBegin) {
		TranscriptCluster* pTranscriptCluster = itTcBegin->second;
		AnnotationConf annot = pTranscriptCluster->GetHighestAnnotationProbesetType();
		if ( annot != CORE ) {
			continue; // only compute expression estimates for core tcs
		}
		//********
		// Expression estimates generated using all core probes
		//********
		int numProbes = pTranscriptCluster->GetNumProbes(CORE);
		if ( numProbes >= 5 ) {
			gsl_matrix* pInten = gsl_matrix_alloc( numSamples, numProbes );
			this->GetIntensityMatrix( pInten, pTranscriptCluster, chips, false );  // don't exclude cross-hyb probes
			if ( summaryMethod == "liwong" ) {
				Model_AffyLiWong* pModel = new Model_AffyLiWong( pInten );
				gsl_vector* pExpression = pModel->GetTheta();
				outfile << pTranscriptCluster->GetId();
				for ( int i = 0; i < numSamples; ++i ) {
					outfile << SEP << gsl_vector_get( pExpression, i );
				}
				outfile << endl;
				delete pModel;
			} else {
				Model_ProbeSelection* pModel = new Model_ProbeSelection( pInten );
				gsl_vector* pExpression = pModel->GetExpression();
				if ( summaryMethod == "selection" ) { // "selection" or "deconvolve" 						
					outfile << pTranscriptCluster->GetId();
					for ( int i = 0; i < numSamples; ++i ) {
						outfile << SEP << gsl_vector_get( pExpression, i );
					}
					outfile << endl;
				} else if ( summaryMethod == "correlation_filter" ) {
					bool toDeconvolve = pTranscriptCluster->CanDeconvolve(); // Any core probes match a unique off-target?
					gsl_vector* pDeconvolveExpr = gsl_vector_alloc( numSamples );
					gsl_vector_memcpy( pDeconvolveExpr, pExpression);
					deconvolvedEstimates.insert( make_pair(pTranscriptCluster->GetId(), make_pair(toDeconvolve, pDeconvolveExpr) ) );
				}
				delete pModel;  
			}
			gsl_matrix_free( pInten );
		}
	}
		
		
		
	// for the correlation filter, iterate, updating expression estimates using the filter criteria
	if ( summaryMethod == "correlation_filter") {
		bool sometcchanged = true;
		int iter = 0;
		const int maxIter = 10;
		while ( sometcchanged == true && iter < maxIter ) {
			
			// Temporary code:  
			// Create a file listing tc ids, number probes excluded and number of total 
			// core probes.  Used to understand how many genes have expression estimates
			// changing when only a small number of probes are suspected of 
			// cross-hybridizing.  
			// Implementation:  
			// At each iteration, overwrite the previous file
			ofstream tempOutput;
			string tempOutputFilename = "TcsWithExcludedProbes.txt";
			tempOutput.open( tempOutputFilename.c_str() );
			tempOutput << "TranscriptClusterId" << "\t" << "NumberTotalCoreProbes" << "\t" << "NumberExcludedCoreProbes" << endl;
			
			sometcchanged = false;
			hash_map<unsigned long, pair<bool, gsl_vector*> >::iterator itDcvBegin = deconvolvedEstimates.begin();
			hash_map<unsigned long, pair<bool, gsl_vector*> >::iterator itDcvEnd = deconvolvedEstimates.end();
			for (; itDcvBegin != itDcvEnd; ++itDcvBegin) {
				bool toDeconvolve = itDcvBegin->second.first;
				if ( toDeconvolve == false ) { // the tc doesn't have any candidate cross-hyb probes for deconvolving
					continue;
				}
				unsigned long tcId = itDcvBegin->first;
				gsl_vector* pDeconvolvedEstimate = itDcvBegin->second.second;
				TranscriptCluster* pTc = pTranscriptClusterDomain->GetTranscriptCluster(tcId);
				int numCoreProbes = pTc->GetNumProbes(CORE,false);
				int probeColCtr = pTc->GetNumProbes(CORE,true);
				gsl_matrix* pInten = gsl_matrix_calloc( numSamples, numCoreProbes); // initialize to 0
				this->GetIntensityMatrix(pInten, pTc, chips, true); // Expr excluding crosshyb probes
				pair<multimap<AnnotationConf, Probeset*>::iterator, multimap<AnnotationConf, Probeset*>::iterator > itProbesetPair = pTc->GetProbesetIterators( CORE );
				multimap<AnnotationConf, Probeset*>::iterator itProbesetBegin = itProbesetPair.first;
				multimap<AnnotationConf, Probeset*>::iterator itProbesetEnd = itProbesetPair.second;
				for (; itProbesetBegin != itProbesetEnd; ++itProbesetBegin) {
					Probeset* pProbeset = itProbesetBegin->second;
					pair<map<unsigned long, Probe*>::iterator, map<unsigned long, Probe*>::iterator> itProbePair = pProbeset->GetProbeIterators();
					map<unsigned long, Probe*>::iterator itProbeBegin = itProbePair.first;
					map<unsigned long, Probe*>::iterator itProbeEnd = itProbePair.second;
					for (; itProbeBegin != itProbeEnd; ++itProbeBegin) {
						Probe* pProbe = itProbeBegin->second;
						const map< string, pair<int, int> >* pChInfo = pProbe->GetCrossHybInfo();
						if ( pChInfo == NULL || pChInfo->size() == 0 ) { 
							continue; 
							// The probe has no matches to off-targets
							// Gene expression does not change
						}
						gsl_vector* pY = gsl_vector_alloc( numSamples );
						this->GetIntensityVector( pY, pProbe, chips );
						double sdY = gsl_stats_sd( pY->data, 1, numSamples );
						double maxOffTargetCorrelation = -1;
						map<string, pair<int, int> >::const_iterator itChInfoBegin = pChInfo->begin();
						map<string, pair<int, int> >::const_iterator itChInfoEnd = pChInfo->end();
						for (; itChInfoBegin != itChInfoEnd; ++itChInfoBegin) {
							string offTargetTcIdString = itChInfoBegin->first;
							pair<int, int> editDistPair = itChInfoBegin->second;
							int editDist = editDistPair.first + editDistPair.second;
							unsigned long offTargetTcId = strtoul( offTargetTcIdString.c_str(), NULL, 10);
							if ( offTargetTcId == 0 ) {
								continue; // conversion from string to unsigned long failed
							}
							hash_map< unsigned long, pair<bool, gsl_vector*> >::iterator itFind = deconvolvedEstimates.find( offTargetTcId );
							if ( itFind == deconvolvedEstimates.end() ) {
								continue; // no expression estimate available
							}
							gsl_vector* pOffTargetEstimate = itFind->second.second;
							// compute correlation of probe intensity with off-target expression
							//******* Check that a double* can be passed as a double[]
							double corValOffTarget = gsl_stats_covariance (pY->data, 1, pOffTargetEstimate->data, 1, numSamples);
							double sdOT = gsl_stats_sd( pOffTargetEstimate->data, 1, numSamples );
							if ( sdY == 0 || sdOT == 0 ) {
								corValOffTarget = 0;
							} else {
								sdY = max( sdY, correlationFilterMinSd  ); // Guard against high correlations due to
								sdOT = max( sdOT, correlationFilterMinSd );// low standard deviations
									corValOffTarget = corValOffTarget / ( sdY * sdOT);
							}
							if ( editDist <= maxEditDistance ) {
								maxOffTargetCorrelation = max( maxOffTargetCorrelation, corValOffTarget );
							} else {
								// Ignore probes with edit distance greater than maxEditDistance
							}
						}
						if ( maxOffTargetCorrelation > minOffTargetCorrelation2Deconvolve ) { 
							// exclude the probe, don't add probe intensity to pInten
						} else {
							// include the probe
							gsl_matrix_set_col (pInten, probeColCtr, pY);	// add probe intensity info to pInten						
							++probeColCtr;							
						}
						gsl_vector_free( pY );
					}
				}
				// re-compute expression estimate
				if ( probeColCtr < 5 ) {
					gsl_vector_set_all( pDeconvolvedEstimate, -1000); // Don't make null because the # of probes could change in different iterations
				} else {
					gsl_matrix* pNewInten = gsl_matrix_alloc(numSamples, probeColCtr);
					gsl_matrix_const_view intenToCopy = gsl_matrix_const_submatrix (pInten, 0, 0, numSamples, probeColCtr);
					gsl_matrix_memcpy(pNewInten, &intenToCopy.matrix );
					Model_ProbeSelection* pModel = new Model_ProbeSelection( pNewInten );
					gsl_vector* pModelExpr = pModel->GetExpression();
					gsl_vector_sub (pDeconvolvedEstimate, pModelExpr);
					double minDiff, maxDiff;
					gsl_vector_minmax (pDeconvolvedEstimate, &minDiff, &maxDiff);
					if ( minDiff < -.0001 || maxDiff > .0001 ) {
						sometcchanged = true;
					}
					gsl_vector_memcpy( pDeconvolvedEstimate, pModelExpr); // Put the right value in pDeconvolvedEstimate
					delete pModel;
					gsl_matrix_free( pNewInten );
				}
				gsl_matrix_free( pInten );
				
				// Temp code:  Output tc id and number excluded probes, if applicable
				if ( probeColCtr < numCoreProbes) {
					tempOutput << tcId << "\t" << numCoreProbes << "\t" << (numCoreProbes - probeColCtr) << endl;
				}
			}
			tempOutput.close();
			++iter;
		}
		// output results
		// Output the final estimates stored in deconvolvedExpression
		// Estimates will only be stored in deconvolvedEstimates if the summarization method is 
		// correlation_filter.  Otherwise, estimates were output to the output file 
		// after generating expression estimates.  
		hash_map< unsigned long, pair<bool,gsl_vector*> >::iterator itBegin = deconvolvedEstimates.begin();
		hash_map< unsigned long, pair<bool, gsl_vector*> >::iterator itEnd = deconvolvedEstimates.end();
		for (; itBegin != itEnd; ++itBegin ) {
			bool toOutput = false;
			stringstream ss("");
			ss << itBegin->first;		
			gsl_vector* pVec = itBegin->second.second;
			for ( int i = 0; i < numSamples; ++i ) {
				double val = gsl_vector_get( pVec, i );
				if ( val != -1000 ) {
					toOutput = true;
				}
				ss << SEP << val;
			}
			ss << endl;
			if ( toOutput ) {
				outfile << ss.str();
			} else if ( outputMasked == true) {
				TranscriptCluster* pTranscriptCluster = pTranscriptClusterDomain->GetTranscriptCluster( itBegin->first );
				if ( pTranscriptCluster != NULL ) {
					int numProbes = pTranscriptCluster->GetNumProbes(CORE);
					if ( numProbes >= 5 ) {
						gsl_matrix* pInten = gsl_matrix_alloc( numSamples, numProbes );
						this->GetIntensityMatrix( pInten, pTranscriptCluster, chips, false );  // don't exclude cross-hyb probes
						if ( summaryMethod == "liwong" ) {
							Model_AffyLiWong* pModel = new Model_AffyLiWong( pInten );
							gsl_vector* pExpression = pModel->GetTheta();
							outfile << pTranscriptCluster->GetId();
							for ( int i = 0; i < numSamples; ++i ) {
								outfile << SEP << gsl_vector_get( pExpression, i );
							}
							outfile << endl;
							delete pModel;
						} else {
							Model_ProbeSelection* pModel = new Model_ProbeSelection( pInten );
							gsl_vector* pExpression = pModel->GetExpression();
							outfile << pTranscriptCluster->GetId();
							for ( int i = 0; i < numSamples; ++i ) {
								outfile << SEP << gsl_vector_get( pExpression, i );
							}
							outfile << endl;
							delete pModel;  
						}
						gsl_matrix_free( pInten );
					}
				}
			}
			gsl_vector_free( pVec );
		} 
	}
		
	
	outfile.close();		
	// cleanup
	this->DeleteMeanQuantiles(); // harmless if no mean quantiles computed
	itChBegin = chips.begin();
	itChEnd = chips.end();
	for (; itChBegin != itChEnd; ++itChBegin ) {
		CHIP_DATA* pChip = *itChBegin;
		delete pChip;
		pChip = NULL;
	}
	return true;
}

void AnalysisDomain::GetIntensityMatrix( gsl_matrix* o_pInten, TranscriptCluster* i_pTranscriptCluster, 
						 vector<CHIP_DATA*>& i_Chips, const bool& i_ExcludeCrossHyb )
{
	int probeCtr = 0;
	pair<multimap<AnnotationConf, Probeset*>::iterator, multimap<AnnotationConf, Probeset*>::iterator > itProbesetPair = i_pTranscriptCluster->GetProbesetIterators( CORE );
	multimap<AnnotationConf, Probeset*>::iterator itProbesetBegin = itProbesetPair.first;
	multimap<AnnotationConf, Probeset*>::iterator itProbesetEnd = itProbesetPair.second;
	for ( ; itProbesetBegin != itProbesetEnd; ++itProbesetBegin ) {
		Probeset* pProbeset = itProbesetBegin->second;		
		pair<map<unsigned long, Probe*>::iterator, map<unsigned long, Probe*>::iterator> itProbePair = pProbeset->GetProbeIterators();
		map<unsigned long, Probe*>::iterator itProbeBegin = itProbePair.first;
		map<unsigned long, Probe*>::iterator itProbeEnd = itProbePair.second;
		for (; itProbeBegin != itProbeEnd; ++itProbeBegin ) {
			Probe* pProbe = itProbeBegin->second;
			if ( i_ExcludeCrossHyb == true && pProbe->IsCrossHyb() == true ) {
				// do nothing
			} else {
				int chipCtr = 0;
				vector<CHIP_DATA*>::iterator itBegin = i_Chips.begin();
				vector<CHIP_DATA*>::iterator itEnd = i_Chips.end();
				for (; itBegin != itEnd; ++itBegin ) {
					CHIP_DATA* pChip = *itBegin;
					double bkgCorrectedNormalized = pChip->GetBackgroundCorrectedNormalizedIntensity( pProbe );
					bkgCorrectedNormalized = max(0, bkgCorrectedNormalized); // truncate at 0, 05/18/08
					gsl_matrix_set( o_pInten, chipCtr, probeCtr, bkgCorrectedNormalized );
					++chipCtr;
				}
				++probeCtr;
			}
		}
	}
}

void AnalysisDomain::GetIntensityVector(  gsl_vector* o_pInten, Probe* i_pProbe, vector<CHIP_DATA*>& i_Chips)
{
	int chipCtr = 0;
	vector<CHIP_DATA*>::iterator itBegin = i_Chips.begin();
	vector<CHIP_DATA*>::iterator itEnd = i_Chips.end();
	for (; itBegin != itEnd; ++itBegin ) {
		CHIP_DATA* pChip = *itBegin;
		double bkgCorrectedNormalized = pChip->GetBackgroundCorrectedNormalizedIntensity(i_pProbe);		
		bkgCorrectedNormalized = max(0, bkgCorrectedNormalized ); // truncate at 0, 5/18/08
		gsl_vector_set( o_pInten, chipCtr, bkgCorrectedNormalized );
		++chipCtr;
	}	
}



bool AnalysisDomain::OutputCrossHybInfo()
{
	return true;
}

bool AnalysisDomain::ModelPresenceAbsence()
{
	// For each chip, 
	// model the background
	// compute presence/absence calls based on the background model.  
	
	Parameters* pParameters = Parameters::GetInstance();
	BackgroundDomain* pBackgroundDomain = BackgroundDomain::GetInstance();
	string backgroundMethod = pParameters->GetModelBackgroundMethod();
	
	pair< set<string>::iterator, set<string>::iterator > itChipPair = pParameters->GetCelFilesIterators();
	set<string>::iterator itChipBegin = itChipPair.first;
	set<string>::iterator itChipEnd = itChipPair.second;
	for (; itChipBegin != itChipEnd; ++itChipBegin) {
		string filename = pParameters->GetArrayCelFolder() + (*itChipBegin);
		ifstream datafile;
		datafile.open(filename.c_str());
		if (datafile.is_open()) {
			string::size_type pos = filename.rfind( Parameters::SEPARATOR );
			int len = static_cast<int>(filename.length());
			string smallname = filename.substr( (pos + 1), (len - 1) );
			//!!!
			// CHIP_DATA is for exon arrays 
			CHIP_DATA* pChip = new CHIP_DATA( smallname, pParameters->GetNumCells(), pParameters->GetCellDim() );
			// Read in the chip's data
			bool readOk = pChip->ReadBinaryCel(filename.c_str());
			if (readOk == false) {
				delete pChip;
				continue;
			}
			// Model background
			bool modelBkgdOk = pChip->ModelBackground();
			if ( modelBkgdOk == false ) {
				delete pChip;
				continue;
			}

			// Create the appropriate P/A model
			// The P/A calls will be output upon creation of the P/A model
			if ( backgroundMethod == "mat" ) {
				PresenceAbsenceModel_MAT* pPA_MAT = new PresenceAbsenceModel_MAT( pChip );
				delete pPA_MAT;
			} else if ( backgroundMethod == "median_gc" ) {
				PresenceAbsenceModel_MedianGC* pPA_MedianGC = new PresenceAbsenceModel_MedianGC( pChip );
				delete pPA_MedianGC;
			} else {
				Parameters::LOGFILE << "Error:  Cannot model presence/absence using background model " << backgroundMethod << endl;
				return false;
			}	
			delete pChip;
		} else {
			Parameters::LOGFILE << "ERROR:  " << "Cannot open CEL file " << filename << endl;
		}
	}
	
	return true;
}





