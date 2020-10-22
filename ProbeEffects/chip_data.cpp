#include "stdafx.h"
#include "chip_data.h"
#include <math.h>
#include "DB_VECTOR.h"
#include "CELFileData.h"
#include "TranscriptClusterDomain.h"
#include "StringTokenizer.h"
using namespace affxcel;
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <algorithm>

CHIP_DATA::CHIP_DATA( const string& i_Id, const int& i_NumCells, const int& i_CellDim  ) 
: ChipBase( i_Id, i_NumCells, i_CellDim )
{
	m_ProbeNormalizationConstant = 0;
}

CHIP_DATA::~CHIP_DATA() 
{	
	// Calls ChipBase deconstructor
}

double CHIP_DATA::GetBkgSubtractedIntensity( Probe* i_pProbe, const gsl_vector* i_pBeta  )
{
	//debugging!!
	 
	if ( i_pBeta == NULL ) {
		double gcbgd = this->GetGCBgd( i_pProbe );
		ProbeBase* pProbeBase = dynamic_cast< ProbeBase* > ( i_pProbe ); // this should never fail
		if ( pProbeBase == NULL ) {
			Parameters::LOGFILE << "ERROR:  Cannot convert Probe to ProbeBase!" << endl;
			return 0;
		}
		double normIntensity = this->GetUnnormalizedIntensity( pProbeBase );
		return (normIntensity - gcbgd);
	} else {
		// Use MAT background correction
		gsl_vector* pDesignVector = gsl_vector_alloc( 80 );
		i_pProbe->SetMATDesignVector( pDesignVector );
		double fittedVal;
		gsl_blas_ddot (pDesignVector, i_pBeta, &fittedVal);
		fittedVal = exp( fittedVal ); 
		double backgroundSubtracted = this->GetUnnormalizedIntensity( i_pProbe ) - fittedVal;
		gsl_vector_free( pDesignVector );
		return backgroundSubtracted;
		
	}
	return 0;
}

double CHIP_DATA::GetGCBgd( Probe* i_pProbe )
{
	map<unsigned short, unsigned short>::iterator it = this->m_MedianGCBgd.find( i_pProbe->GetGcCount() );
	int gcbgd;
	if (it == this->m_MedianGCBgd.end() ) {
		gcbgd = 0; 
	} else {
		gcbgd = it->second;
	}
		
	return gcbgd;
}


double CHIP_DATA::GetGCBgdPValue( Probe* i_pProbe )
{
	int gcCount = i_pProbe->GetGcCount();
	unsigned short inten = this->GetUnnormalizedIntensity( i_pProbe );
	double pVal;
	int numBackgroundProbes = m_BackgroundData[ gcCount ].size();
	if ( numBackgroundProbes < 10 ) { // This may occur at very low GC counts or very high GC counts
		bool searchHigher;			// We simply search the next up or next down bin.  
		if ( gcCount < 13 ) {		// This has the potential for errors if there are low numbers of background probes.  
			searchHigher = true;
		} else {
			searchHigher = false;
		}
		while ( numBackgroundProbes < 10 ) {
			if (searchHigher ) {
				numBackgroundProbes = m_BackgroundData[ (++gcCount) ].size(); // adjust the value of gcCount
			} else {
				numBackgroundProbes = m_BackgroundData[ (--gcCount) ].size();
			}
		}
	} 

	if ( numBackgroundProbes == 1 ) { // This should never happen
		unsigned short theval = m_BackgroundData[gcCount][0];
		if ( inten < theval) {
			pVal=0.0; // perl dabg returns 1; we invert
		}
		else if (inten > theval) {
			pVal=1.0; // perl dabg returns 0; we invert
		}
		else {
			pVal=0.5;
		}
	} else {
		// l_idx is the index of the value ">=" intensity
		int l_idx=lower_bound(m_BackgroundData[gcCount].begin(),m_BackgroundData[gcCount].end(),inten)
		- m_BackgroundData[gcCount].begin();
		// fix it up to be less than intensity
		if (l_idx == numBackgroundProbes ) {
			l_idx = numBackgroundProbes - 1;
		}
		pVal=( static_cast<double>(l_idx) )/( static_cast<double>(numBackgroundProbes) );
	}

	// now that we have our pval, invert it for our chances of being greater.
	pVal=(1.0-pVal);
  
	return pVal;
}

void CHIP_DATA::SetGCBinBackgroundData()
{
	// load data into m_BackgroundData
	TranscriptClusterDomain* pTranscriptClusterDomain = TranscriptClusterDomain::GetInstance();
	for ( int gcbin = 0; gcbin < 26; ++gcbin ) {
		// First create a vector to store background probe intensities having gc content gcbin
		vector<unsigned short> myVec;
		this->m_BackgroundData.push_back( myVec );

		// genomic background data
		pair< multimap<int, Probe*>::iterator, multimap<int, Probe*>::iterator > itPair = pTranscriptClusterDomain->GetGenomicBgdIterators( gcbin );
		multimap<int, Probe*>::iterator itBegin = itPair.first;
		multimap<int, Probe*>::iterator itEnd = itPair.second;
		for (; itBegin != itEnd; ++itBegin) {
			Probe* pProbe = itBegin->second;
			unsigned short inten = this->GetUnnormalizedIntensity( pProbe );	
			this->m_BackgroundData[gcbin].push_back( inten );
		}
		// antigenomic background data 
		itPair = pTranscriptClusterDomain->GetAntigenomicBgdIterators( gcbin );
		itBegin = itPair.first;
		itEnd = itPair.second;
		for (; itBegin != itEnd; ++itBegin) {
			Probe* pProbe = itBegin->second;
			unsigned short inten = this->GetUnnormalizedIntensity( pProbe );	
			this->m_BackgroundData[gcbin].push_back( inten );
		}
	}

	// For each vector of background probes of a given gc content, sort the values
	int size = static_cast<int>( this->m_BackgroundData.size() );
	for (int gcbin=0; gcbin < size; ++gcbin) {
		sort( this->m_BackgroundData[gcbin].begin(),this->m_BackgroundData[gcbin].end() );
	}

}


void CHIP_DATA::SetProbeNormalizationConstant()
{
	return this->SetProbeNormalizationConstant( m_pBackgroundModel->GetBackgroundCoefficients() );
}

void CHIP_DATA::SetProbeNormalizationConstant( const gsl_vector* i_pBeta )
{
	Parameters* pParameters = Parameters::GetInstance();
	string normMethod = pParameters->GetNormalizationMethod();
	if ( this->m_ProbeNormalizationConstant == 0 ) { // Compute the normalization constant
		if ( normMethod == "core_probe_scaling" ) {
			
			// *****
			// For each core probe belonging to a transcript cluster
			// Compute the background subtracted value
			// Store Unnormalized - background in coreBackgroundCorrectedValues
			// Compute the scaling factor so that median core intensity = 100
			// *****
			const double targetMedianCoreIntensity = 100;
			vector<double> coreBackgroundCorrectedValues;
					
			TranscriptClusterDomain* pTranscriptClusterDomain = TranscriptClusterDomain::GetInstance();
			pair< hash_map<unsigned long, TranscriptCluster*>::iterator, hash_map<unsigned long, TranscriptCluster*>::iterator > itTcPair = pTranscriptClusterDomain->GetTranscriptClusterIterators();
			hash_map<unsigned long, TranscriptCluster*>::iterator itTcBegin = itTcPair.first;
			hash_map<unsigned long, TranscriptCluster*>::iterator itTcEnd = itTcPair.second;
			AnnotationConf coreAnnot = CORE;
			for (; itTcBegin != itTcEnd; ++itTcBegin) {
				TranscriptCluster* pTranscriptCluster = itTcBegin->second;
				pair<multimap<AnnotationConf, Probeset*>::iterator, multimap<AnnotationConf, Probeset*>::iterator > itProbesetPair = pTranscriptCluster->GetProbesetIterators( coreAnnot );
				multimap<AnnotationConf, Probeset*>::iterator itProbesetBegin = itProbesetPair.first;
				multimap<AnnotationConf, Probeset*>::iterator itProbesetEnd = itProbesetPair.second;
				for (; itProbesetBegin != itProbesetEnd; ++itProbesetBegin) {
					Probeset* pProbeset = itProbesetBegin->second;
					pair< map<unsigned long, Probe*>::iterator, map<unsigned long, Probe*>::iterator > itProbePair = pProbeset->GetProbeIterators();
					map<unsigned long, Probe*>::iterator itProbeBegin = itProbePair.first;
					map<unsigned long, Probe*>::iterator itProbeEnd = itProbePair.second;
					for (; itProbeBegin != itProbeEnd; ++itProbeBegin) {
						Probe* pProbe = itProbeBegin->second;
						double unnormalized = this->GetUnnormalizedIntensity( pProbe );
						double background = this->m_pBackgroundModel->GetBackgroundPrediction( pProbe );
						coreBackgroundCorrectedValues.push_back( unnormalized - background );
					}
				}
			}
		
			// Get the median intensity of core probes 'unnormalized - background'
			sort( coreBackgroundCorrectedValues.begin(), coreBackgroundCorrectedValues.end() );
			int medianLocation = coreBackgroundCorrectedValues.size() / 2; // integer division
			double medianCoreIntensity = coreBackgroundCorrectedValues.operator []( medianLocation );
			if ( medianCoreIntensity <= 0 ) {
				Parameters::LOGFILE << "WARNING:  Median core intensity of array " 
					<< this->m_Id << " equals " 
					<< medianCoreIntensity << ".  Normalization scaling factor will be set to 100." 
					<< endl;
				medianCoreIntensity = 1;
			}
			this->m_ProbeNormalizationConstant = targetMedianCoreIntensity / medianCoreIntensity;
		} else {
			this->m_ProbeNormalizationConstant = 1;
		}
	}
}

double CHIP_DATA::GetBackgroundCorrectedNormalizedIntensity( Probe* i_pProbe )
{
	if ( this->m_pBackgroundModel == NULL ) {
		this->ModelBackground();
	}
	if ( this->m_ProbeNormalizationConstant == 0 ) {
		this->SetProbeNormalizationConstant();
	}
	double unnormalized = this->GetUnnormalizedIntensity(i_pProbe); // quantile normalization already taken care of
	double background = this->m_pBackgroundModel->GetBackgroundPrediction(i_pProbe);
	return (unnormalized - background)*this->m_ProbeNormalizationConstant;
	
}

//
// Store the median normalized intensity of background probes.  
// Also finds the variance of the background probe intensities.  
// The intention is to remove the intensity due to non-specific binding (ie cross-hybridization)
// Check whether normalized intensities should be used for 3' arrays...  
// Both genomic and antigenomic background probes are used!  
//
bool CHIP_DATA::SetGCBgd( ) {
	Parameters* pParameters = Parameters::GetInstance();
	int maxGCContent = Parameters::PROBELENGTH;
	TranscriptClusterDomain* pTranscriptClusterDomain = TranscriptClusterDomain::GetInstance();
	pair< multimap<int, Probe*>::iterator, multimap<int, Probe*>::iterator > itPair;
	bool* pIsGCSet = new bool[maxGCContent + 1];
	for (int i = 0; i <= maxGCContent; ++i) {
		pIsGCSet[i] = false;
		int len = 
			pTranscriptClusterDomain->GetNumGenomicBgdProbes( i ) + 
			pTranscriptClusterDomain->GetNumAntigenomicBgdProbes( i );//how big to make the DB_VECTOR
		if (len == 0) {
			continue;										// Later the value will be set appropriately.  
		}
		DB_VECTOR* pVec = new DB_VECTOR(len);
		// Ask ProbesetDomain for iterators of probes with a certain GC content
		itPair = pTranscriptClusterDomain->GetGenomicBgdIterators( i );
		// For each genomic and antigenomic probe, add the probe's intensity to a DB_VECTOR
		multimap<int, Probe*>::iterator itBeg = itPair.first;
		multimap<int, Probe*>::iterator itEnd = itPair.second;
		for (; itBeg != itEnd; ++itBeg) {
			Probe* pProbe = itBeg->second;
			pair<int, int> posPair = pProbe->GetPosition();
			int idx = Parameters::IDX( posPair.first, posPair.second, pParameters->GetCellDim() );
			pVec->Add( this->m_pInten[idx] );
		}
		itPair = pTranscriptClusterDomain->GetAntigenomicBgdIterators( i );
		itBeg = itPair.first;
		itEnd = itPair.second;
		for (; itBeg != itEnd; ++itBeg) {
			Probe* pProbe = itBeg->second;
			pair<int, int> posPair = pProbe->GetPosition();
			int idx = Parameters::IDX( posPair.first, posPair.second, pParameters->GetCellDim() );
			pVec->Add( this->m_pInten[idx] );
		}

		// Ask the DB_VECTOR for its median
		double median = pVec->GetMedian();
		// Ask the DB_VECTOR for its std. deviation.  
		double var = pVec->GetVar();
		// Set this->m_MedianGCBgd[i] = median value
		this->m_MedianGCBgd.insert( make_pair(i, static_cast<unsigned short>(median + .5)) );
		// Set this->m_VarianceGCBgd[i] = var
		//this->m_VarGCBgd.insert( make_pair(i, var ) );
		delete pVec;
		pIsGCSet[i] = true;
	}

	// Now we must go through those gc-content values which do not have any background probes
	// We set their median intensity equal to the median intensity of the closest bin with background probes
	for (int i = 0; i < maxGCContent; ++i ) {
		if ( pIsGCSet[i] == false ) {
			if ( i < (maxGCContent / 2 ) ) { // We search for values larger than i 
				for ( int j = (i + 1); j < maxGCContent; ++j ) {
					if ( pIsGCSet[j] == true ) {
						unsigned short medianVal = (this->m_MedianGCBgd.find( j ))->second;
						this->m_MedianGCBgd.insert( make_pair(i, medianVal) );
						pIsGCSet[i] = true;
						break;
					}
				}
			} else { // We search for values smaller than i
				for ( int j = (i - 1); j >= 0; --j ) {
					if ( pIsGCSet[j] == true ) {
						unsigned short medianVal = (this->m_MedianGCBgd.find( j ))->second;
						this->m_MedianGCBgd.insert( make_pair(i, medianVal) );
						pIsGCSet[i] = true;
						break;
					}
				}
			}
		}
	}
	delete [] pIsGCSet;
	return true;
}

void CHIP_DATA::SetBackgroundSubtractedNormalizedIntensities( gsl_vector* i_pBeta, const double& i_ScalingConstant )
{
	// Sets the m_pInten[] values to the background corrected, normalized values
	Parameters* pParameters = Parameters::GetInstance();
	gsl_vector* pDesignVector = gsl_vector_alloc( pParameters->GetNumParameters_MAT() );
	
	// antigenomic probe iterators
	pair<multimap<int, Probe*>::iterator, multimap<int, Probe*>::iterator > antigenomicItPair = TranscriptClusterDomain::GetInstance()->GetAntigenomicBgdIterators();
	multimap<int, Probe*>::iterator antigenomicItBegin = antigenomicItPair.first;
	multimap<int, Probe*>::iterator antigenomicItEnd = antigenomicItPair.second;
	for (; antigenomicItBegin != antigenomicItEnd; ++antigenomicItBegin) {
		Probe* pProbe = antigenomicItBegin->second;
		pProbe->SetMATDesignVector( pDesignVector );
		pair<int, int> posPair = pProbe->GetPosition();
		int probeIdx = Parameters::IDX(posPair.first, posPair.second, this->m_CellDim );
		double unnormalized = this->m_pInten[ probeIdx ];
		double logBackground;
		gsl_blas_ddot (pDesignVector, i_pBeta, &logBackground);
		double background = exp( logBackground ); // not on log scale
		double backgroundCorrectedIntensity = (unnormalized - background) * i_ScalingConstant;
		this->m_pInten[ probeIdx ] = static_cast<unsigned short>( backgroundCorrectedIntensity );
	}
	// genomic probe iterators
	pair<multimap<int, Probe*>::iterator, multimap<int, Probe*>::iterator > genomicItPair = TranscriptClusterDomain::GetInstance()->GetGenomicBgdIterators();
	multimap<int, Probe*>::iterator genomicItBegin = genomicItPair.first;
	multimap<int, Probe*>::iterator genomicItEnd = genomicItPair.second;
	for (; genomicItBegin != genomicItEnd; ++genomicItBegin) {
		Probe* pProbe = genomicItBegin->second;
		pProbe->SetMATDesignVector( pDesignVector );
		pair<int, int> posPair = pProbe->GetPosition();
		int probeIdx = Parameters::IDX(posPair.first, posPair.second, this->m_CellDim );
		double unnormalized = this->m_pInten[ probeIdx ];
		double logBackground;
		gsl_blas_ddot (pDesignVector, i_pBeta, &logBackground);
		double background = exp( logBackground ); // not on log scale
		double backgroundCorrectedIntensity = (unnormalized - background) * i_ScalingConstant;
		this->m_pInten[ probeIdx ] = static_cast<unsigned short>( backgroundCorrectedIntensity );
	}


	// transcript cluster iterators
	pair< hash_map<unsigned long, TranscriptCluster*>::iterator, hash_map<unsigned long, TranscriptCluster*>::iterator > itTcPair = TranscriptClusterDomain::GetInstance()->GetTranscriptClusterIterators();
	hash_map<unsigned long, TranscriptCluster*>::iterator itTcBegin = itTcPair.first;
	hash_map<unsigned long, TranscriptCluster*>::iterator itTcEnd = itTcPair.second;
	for (; itTcBegin != itTcEnd; ++itTcBegin) {
		TranscriptCluster* pTranscriptCluster = itTcBegin->second;
		pair< multimap<AnnotationConf, Probeset*>::iterator, multimap<AnnotationConf, Probeset*>::iterator > itProbesetPair = pTranscriptCluster->GetProbesetIterators();
		multimap<AnnotationConf, Probeset*>::iterator itProbesetBegin = itProbesetPair.first;
		multimap<AnnotationConf, Probeset*>::iterator itProbesetEnd = itProbesetPair.second;
		for (; itProbesetBegin != itProbesetEnd; ++itProbesetBegin) {
			// For each probeset
			Probeset* pProbeset = itProbesetBegin->second;
			pair<map<unsigned long, Probe*>::iterator, map<unsigned long, Probe*>::iterator> itProbePair = pProbeset->GetProbeIterators();
			map<unsigned long, Probe*>::iterator itProbeBegin = itProbePair.first;
			map<unsigned long, Probe*>::iterator itProbeEnd = itProbePair.second;
			// For each probe
			for (; itProbeBegin != itProbeEnd; ++itProbeBegin) {
				Probe* pProbe = itProbeBegin->second;
				pProbe->SetMATDesignVector( pDesignVector );
				pair<int, int> posPair = pProbe->GetPosition();
				int probeIdx = Parameters::IDX(posPair.first, posPair.second, this->m_CellDim );
				double unnormalized = this->m_pInten[ probeIdx ];
				double logBackground;
				gsl_blas_ddot (pDesignVector, i_pBeta, &logBackground);
				double background = exp( logBackground ); // not on log scale
				double backgroundCorrectedIntensity = (unnormalized - background) * i_ScalingConstant;
				this->m_pInten[ probeIdx ] = static_cast<unsigned short>( backgroundCorrectedIntensity );
			}
		}

	}
	gsl_vector_free( pDesignVector );
}

void CHIP_DATA::OutputBkgdCorrectNormProbesLine( Probe* i_pProbe, gsl_vector* i_pDesignVec, 
												gsl_matrix* i_pCov, 
												gsl_vector* i_pBeta, const double& i_NormalizationConstant, 
												const double& i_SigmaHat, const string& i_ProbesetId, 
												const string& i_ProbesetAnnotationType, const string& i_TcId, 
												ofstream& i_Outfile) 
{	
	int numParameters = i_pBeta->size;
	const string SEP = "\t";
	i_pProbe->SetMATDesignVector( i_pDesignVec );
	double unnormalized = this->GetUnnormalizedIntensity( i_pProbe );
	double logBackground;
	gsl_blas_ddot (i_pDesignVec, i_pBeta, &logBackground);
	double background = exp( logBackground ); // not on log scale
	double backgroundCorrectedIntensity = (unnormalized - background) * i_NormalizationConstant;
	double sqPredError; // Estimate the error from the regression
	gsl_vector* pTempVec = gsl_vector_alloc( numParameters );
	gsl_blas_dgemv (CblasNoTrans, 1, i_pCov, i_pDesignVec, 0, pTempVec);
	gsl_blas_ddot (i_pDesignVec, pTempVec, (&sqPredError));
	sqPredError += gsl_pow_2(i_SigmaHat); // Prediction error of logPM - logBgd.  include sigmaHatSq
	gsl_vector_free( pTempVec );

	double stdError = sqrt( 
		gsl_pow_2(i_NormalizationConstant) 
		* ( exp(sqPredError) - 1) 
		* exp( 2*( logBackground ) + sqPredError ) 
		);
								
	// Output the probe id, the observed value, the fitted value, 
	i_Outfile << i_TcId << SEP
		<< i_ProbesetId << SEP
		<< i_ProbesetAnnotationType << SEP
		<< i_pProbe->GetProbeId() << SEP 
		<< unnormalized << SEP
		<< backgroundCorrectedIntensity << SEP 
		<< stdError << SEP 
		<< i_NormalizationConstant << SEP
		<< logBackground << SEP 
		<< sqPredError << endl; // To determine error models  
}

void CHIP_DATA::OutputBkgdCorrectNormProbesHeader( ofstream& i_Outfile )
{
	const string SEP = "\t";
	i_Outfile 
		<< "TranscriptClusterId" 
		<< SEP 
		<< "ProbesetId" 
		<< SEP 
		<< "ProbesetAnnotationLevel"
		<< SEP
		<< "ProbeId" 
		<< SEP 
		<< "RawProbeIntensity"
		<< SEP
		<< "BackgroundCorrectedMedianScaledIntensity" 
		<< SEP
		<< "StandardError"
		<< SEP 
		<< "NormalizationConstant" 
		<< SEP
		<< "PredictedLogBackground"
		<< SEP
		<< "MATSigmaHatSquared"
		<< SEP
		<< "ProbeSequence" // added 12/07
		<< endl;
}



void CHIP_DATA::OutputAllBkgdCorrectNormProbes( )
{
	this->OutputAllBkgdCorrectNormProbes(this->m_pBackgroundModel->GetBackgroundCoefficients(), this->m_ProbeNormalizationConstant, 0, NULL);
}

void CHIP_DATA::OutputAllBkgdCorrectNormProbes( const gsl_vector* i_pBeta, const double& i_NormalizationConstant, const double& i_SigmaHat, gsl_matrix* i_pCov )
{
	Parameters* pParameters = Parameters::GetInstance();
	TranscriptClusterDomain* pTranscriptClusterDomain = TranscriptClusterDomain::GetInstance();
	string bkgdMethod = pParameters->GetModelBackgroundMethod();	
	bool excludeCrossHybProbes = false;
	string crossHybType = pParameters->GetCrossHybType();
	if ( crossHybType == "filter" ) {
		excludeCrossHybProbes = true;
	}
	
	const string SEP = "\t";
	string fittedFile = pParameters->GetBkgdCorrectNormProbesFile();
	//string::size_type pos = fittedFile.rfind( "." ); // removed 04/21/08
	//if ( pos != string::npos ) { // if there a ".*" at the end, append the chip name and add ".txt"
	//	fittedFile = fittedFile.substr( 0, pos-1 );
	//}
    fittedFile = fittedFile + "_" + this->GetId() + "_" + bkgdMethod + ".txt";
	ofstream fittedOutfile;
	fittedOutfile.open( fittedFile.c_str() );
	this->OutputBkgdCorrectNormProbesHeader( fittedOutfile );
				

	// First output tcs with at least one core probe, at least one extended probe
	// at least one full probe, with free/ambiguous probes
	vector<TranscriptCluster*> coreTcs;
	vector<TranscriptCluster*> extendedTcs;
	vector<TranscriptCluster*> fullTcs;
	vector<TranscriptCluster*> freeAmbTcs;

	// Get transcript cluster iterators from TranscriptClusterDomain
	pair< hash_map<unsigned long, TranscriptCluster* >::iterator, hash_map<unsigned long, TranscriptCluster*>::iterator > itTcPair = pTranscriptClusterDomain->GetTranscriptClusterIterators();
	hash_map<unsigned long, TranscriptCluster*>::iterator itTcBegin = itTcPair.first;
	hash_map<unsigned long, TranscriptCluster*>::iterator itTcEnd = itTcPair.second;
	for (; itTcBegin != itTcEnd; ++itTcBegin) {
		TranscriptCluster* pTranscriptCluster = itTcBegin->second;
		AnnotationConf annot = pTranscriptCluster->GetHighestAnnotationProbesetType();
		if ( annot == CORE ) {
			coreTcs.push_back( pTranscriptCluster );
		} else if ( annot == EXTENDED ) {
			extendedTcs.push_back( pTranscriptCluster );
		} else if ( annot == FULL ) {
			fullTcs.push_back( pTranscriptCluster );
		} else {
			freeAmbTcs.push_back( pTranscriptCluster );
		}
	}
				
	map<AnnotationConf, vector<TranscriptCluster*>*> tcAnnotSets;
	tcAnnotSets.insert( make_pair(CORE, &coreTcs) );
	tcAnnotSets.insert( make_pair(EXTENDED, &extendedTcs) );
	tcAnnotSets.insert( make_pair(FULL, &fullTcs) );
	tcAnnotSets.insert( make_pair(FREE, &freeAmbTcs) );
	map<AnnotationConf, vector<TranscriptCluster*>*>::iterator itTcByAnnotBegin = tcAnnotSets.begin();
	map<AnnotationConf, vector<TranscriptCluster*>*>::iterator itTcByAnnotEnd = tcAnnotSets.end();
	for (; itTcByAnnotBegin != itTcByAnnotEnd; ++itTcByAnnotBegin) {
		AnnotationConf annotSet = itTcByAnnotBegin->first;
		if ( annotSet == CORE ) {
			fittedOutfile << "#" << "CoreTranscriptClusters" << endl;
		} else if ( annotSet == EXTENDED ) {
			fittedOutfile << "#" << "ExtendedTranscriptClusters" << endl;
		} else if ( annotSet == FULL ) {
			fittedOutfile << "#" << "FullTranscriptClusters" << endl;
		} else if ( annotSet == FREE ) {
			fittedOutfile << "#" << "FreeAndAmbiguousTranscriptClusters" << endl;
		}

		vector<TranscriptCluster*>::iterator itTcBegin = itTcByAnnotBegin->second->begin();
		vector<TranscriptCluster*>::iterator itTcEnd = itTcByAnnotBegin->second->end();
		for (; itTcBegin != itTcEnd; ++itTcBegin) {
			TranscriptCluster* pTranscriptCluster = *itTcBegin;
			unsigned long tcId = pTranscriptCluster->GetId();
			ostringstream ostr;
			ostr << tcId;
			string tcIdString = ostr.str();
			string tcIdStringInclCrossHyb = tcIdString;
			ostr.str(""); // ostr will store info from all probes, including cross-hyb probes
			pair< multimap<AnnotationConf, Probeset*>::iterator, multimap<AnnotationConf, Probeset*>::iterator > itProbesetPair = pTranscriptCluster->GetProbesetIterators();
			multimap<AnnotationConf, Probeset*>::iterator itProbesetBegin = itProbesetPair.first;
			multimap<AnnotationConf, Probeset*>::iterator itProbesetEnd = itProbesetPair.second;
			for (; itProbesetBegin != itProbesetEnd; ++itProbesetBegin) {
				Probeset* pProbeset = itProbesetBegin->second;
				unsigned long probesetId = pProbeset->GetProbesetId();
				AnnotationConf annot = itProbesetBegin->first;
				string probesetAnnotationType = "full"; // most frequent type
				if ( annot == CORE ) {
					probesetAnnotationType = "core";
				} else if ( annot == EXTENDED ) {
					probesetAnnotationType = "extended";
				} else if ( annot == FREE ) {
					probesetAnnotationType = "free";
				} else if ( annot == AMBIGUOUS ) {
					probesetAnnotationType = "ambiguous";
				} else if ( annot == NONE ) {
					probesetAnnotationType = "none";
				}
					
				pair<map<unsigned long, Probe*>::iterator, map<unsigned long, Probe*>::iterator> itProbePair = pProbeset->GetProbeIterators();
				map<unsigned long, Probe*>::iterator itProbeBegin = itProbePair.first;
				map<unsigned long, Probe*>::iterator itProbeEnd = itProbePair.second;
				// For each probe
				for (; itProbeBegin != itProbeEnd; ++itProbeBegin) {
					Probe* pProbe = itProbeBegin->second;
					bool isCrossHyb = pProbe->IsCrossHyb();
					if ( excludeCrossHybProbes && (isCrossHyb == true) ) {
						continue; // don't output info for cross-hybridizing probes
					}
					double background = this->m_pBackgroundModel->GetBackgroundPrediction( pProbe );
					double logBackground = log( max(background, 1.0) );
					double sqPredError = 0;
					double stdError = 0;
					double unnormalized = this->GetUnnormalizedIntensity( pProbe );
					double backgroundCorrectedIntensity = (unnormalized - background) * i_NormalizationConstant;
					
					// put all info in the stringstream
					if ( isCrossHyb == false ) {
						// Output the probe id, the observed value, the fitted value, 
						fittedOutfile << tcIdString << SEP
						<< probesetId << SEP
						<< probesetAnnotationType << SEP
						<< pProbe->GetProbeId() << SEP 
						<< unnormalized << SEP
						<< backgroundCorrectedIntensity << SEP 
						<< stdError << SEP 
						<< i_NormalizationConstant << SEP
						<< logBackground << SEP 
						<< sqPredError << SEP // To determine error models  
						<< pProbe->GetSequence() << endl; // added the probe sequence, 12/07
					}
				}
			}
			fittedOutfile << ostr.str(); // ostr is "" if outputTwoEstimates == false
		}
	}

	// Get background probe iterators from Transcript Cluster Domain
	fittedOutfile << "#" << "BackgroundProbes" << endl;
	pair< multimap<int, Probe*>::iterator, multimap<int, Probe*>::iterator > itGenomicBgdPair = pTranscriptClusterDomain->GetGenomicBgdIterators();
	multimap<int, Probe*>::iterator itGenomicBgdBegin = itGenomicBgdPair.first;
	multimap<int, Probe*>::iterator itGenomicBgdEnd = itGenomicBgdPair.second;
	for (; itGenomicBgdBegin != itGenomicBgdEnd; ++itGenomicBgdBegin) {
		Probe* pProbe = itGenomicBgdBegin->second;
		if ( excludeCrossHybProbes && pProbe->IsCrossHyb() ) {
			continue; // don't output info for cross-hybridizing probes
		}	
		double background = this->m_pBackgroundModel->GetBackgroundPrediction( pProbe );
		double logBackground = log( max(background, 1.0) );
		double sqPredError = 0;
		double stdError = 0;
		double unnormalized = this->GetUnnormalizedIntensity( pProbe );
		double backgroundCorrectedIntensity = (unnormalized - background) * i_NormalizationConstant;
					
		// Output the probe id, the observed value, the fitted value, 
		fittedOutfile << "genomic_background" << SEP	
			<< "genomic_background" << SEP
			<< "genomic_background" << SEP
			<< pProbe->GetProbeId() << SEP 
			<< unnormalized << SEP
			<< backgroundCorrectedIntensity << SEP 
			<< stdError << SEP 
			<< i_NormalizationConstant << SEP
			<< logBackground << SEP 
			<< sqPredError << SEP // To determine error models  
			<< pProbe->GetSequence() << endl; // added the probe sequence, 02/08 
	}
	pair< multimap<int, Probe*>::iterator, multimap<int, Probe*>::iterator > itAntigenomicBgdPair = pTranscriptClusterDomain->GetAntigenomicBgdIterators();
	multimap<int, Probe*>::iterator itAntigenomicBgdBegin = itAntigenomicBgdPair.first;
	multimap<int, Probe*>::iterator itAntigenomicBgdEnd = itAntigenomicBgdPair.second;
	for (; itAntigenomicBgdBegin != itAntigenomicBgdEnd; ++itAntigenomicBgdBegin) {
		Probe* pProbe = itAntigenomicBgdBegin->second;
		if ( excludeCrossHybProbes && pProbe->IsCrossHyb() ) {
			continue; // don't output info for cross-hybridizing probes
		}
		double background = this->m_pBackgroundModel->GetBackgroundPrediction( pProbe );
		double logBackground = log( max(background, 1.0) );
		double sqPredError = 0;
		double stdError = 0;
		double unnormalized = this->GetUnnormalizedIntensity( pProbe );
		double backgroundCorrectedIntensity = (unnormalized - background) * i_NormalizationConstant;
					
		// Output the probe id, the observed value, the fitted value, 
		fittedOutfile << "antigenomic_background" << SEP
			<< "antigenomic_background" << SEP
			<< "antigenomic_background" << SEP
			<< pProbe->GetProbeId() << SEP 
			<< unnormalized << SEP
			<< backgroundCorrectedIntensity << SEP 
			<< stdError << SEP 
			<< i_NormalizationConstant << SEP
			<< logBackground << SEP 
			<< sqPredError << SEP // To determine error models 
			<< pProbe->GetSequence() << endl; // added the probe sequence, 02/08 
	}
	fittedOutfile.close();
}


