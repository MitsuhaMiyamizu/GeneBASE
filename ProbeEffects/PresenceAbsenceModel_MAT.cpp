/*
 *  PresenceAbsenceModel_MAT.cpp
 *  ProbeEffects2.0.1
 *
 *  Created by Karen Kapur on 8/11/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include "PresenceAbsenceModel_MAT.h"
#include "stats-distributions.h"
#include "TranscriptClusterDomain.h"
#include "stdafx.h"
#include "BackgroundModel_MAT.h"
#include "BackgroundModel_MedianGC.h"
#include <gsl/gsl_cdf.h>

PresenceAbsenceModel_MAT::PresenceAbsenceModel_MAT( ChipBase* i_pChip ) {
	
	// For each tc
	// for each core probe
	// compute a p-value or score
	// end
	// compute a p-value or score for the tc
	// output the value, df, p-value
	//end	
	
	Parameters* pParameters = Parameters::GetInstance();
	TranscriptClusterDomain* pTranscriptClusterDomain = TranscriptClusterDomain::GetInstance();
	string outputFile = pParameters->GetPresenceAbsenceOutfile() + "_" + i_pChip->GetId() + ".txt";
	ofstream output;
	output.open( outputFile.c_str() );
	const string sep = "\t";
	if ( output.is_open() == true ) {
		output << "TranscriptClusterId" << sep << "Score" << sep << "NumberCoreProbes" << sep << "PValue" << endl;
		pair< hash_map<unsigned long, TranscriptCluster*>::iterator, hash_map<unsigned long, TranscriptCluster*>::iterator> itTcPair = pTranscriptClusterDomain->GetTranscriptClusterIterators();
		hash_map<unsigned long, TranscriptCluster*>::iterator itTcBegin = itTcPair.first;
		hash_map<unsigned long, TranscriptCluster*>::iterator itTcEnd = itTcPair.second;
		for (; itTcBegin != itTcEnd; ++itTcBegin ) {
			unsigned long tcId = itTcBegin->first;
			TranscriptCluster* pTc = itTcBegin->second;
			AnnotationConf annot = CORE;
			vector<double> zScores;
			
			pair<multimap<AnnotationConf, Probeset*>::iterator, multimap<AnnotationConf, Probeset*>::iterator > itProbesetPair = pTc->GetProbesetIterators( annot );
			multimap<AnnotationConf, Probeset*>::iterator itProbesetBegin = itProbesetPair.first;
			multimap<AnnotationConf, Probeset*>::iterator itProbesetEnd = itProbesetPair.second;
			for (; itProbesetBegin != itProbesetEnd; ++itProbesetBegin ) {
				Probeset* pProbeset = itProbesetBegin->second;
				pair<map<unsigned long, Probe*>::iterator, map<unsigned long, Probe*>::iterator> itProbePair = pProbeset->GetProbeIterators();
				map<unsigned long, Probe*>::iterator itProbeBegin = itProbePair.first;
				map<unsigned long, Probe*>::iterator itProbeEnd = itProbePair.second;
				for (; itProbeBegin != itProbeEnd; ++itProbeBegin ) {
					Probe* pProbe = itProbeBegin->second;
					double unnormalized = i_pChip->GetUnnormalizedIntensity( pProbe );
					const BackgroundModel* pBackgroundModel = i_pChip->GetBackgroundModel();
					double backgroundPrediction = pBackgroundModel->GetBackgroundPrediction(pProbe);
					const BackgroundModel_MAT* pBkgMod_MAT = dynamic_cast<const BackgroundModel_MAT*>( pBackgroundModel );
					if ( pBkgMod_MAT == NULL ) {
						Parameters::LOGFILE << "Error:  Cannot apply MAT presence absence modeling.  Unknown background model.  " << endl;
					}
					double predError = pBkgMod_MAT->GetBackgroundPredictionStandardError(pProbe);
					double zScore = ( log( max(unnormalized, 1.0) ) - log( max(backgroundPrediction, 1.0) ) ) / predError;
					zScores.push_back(zScore);
				}
			}
			double probeCount = static_cast<double>( zScores.size() );
			if ( probeCount > 0 ) {
				double test_statistic = 0;
				vector<double>::iterator itScoresBegin = zScores.begin();
				vector<double>::iterator itScoresEnd = zScores.end();
				for (; itScoresBegin != itScoresEnd; ++itScoresBegin) {
					double score = *itScoresBegin;
					test_statistic += score; // add the z-scores
				}
				test_statistic = test_statistic / sqrt( probeCount ); // this ratio should be approximaely N(0,1)
				double pValue = gsl_cdf_ugaussian_Q ( test_statistic );
				output << tcId << sep << test_statistic << sep << probeCount << sep << pValue << endl;			
			}
		}
		output.close();
	
	} else {
		Parameters::LOGFILE << "Error:  Cannot open presence absence output file." << endl;
	}
}

PresenceAbsenceModel_MAT::~PresenceAbsenceModel_MAT() {
	
}





