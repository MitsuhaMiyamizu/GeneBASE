/*
 *  PresenceAbsenceModel_MedianGC.cpp
 *  ProbeEffects2.0.1
 *
 *  Created by Karen Kapur on 8/11/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include "PresenceAbsenceModel_MedianGC.h"
#include "BackgroundDomain.h"
#include "ProbeBase.h"
#include "TranscriptClusterDomain.h"
#include <math.h>
#include "stats-distributions.h"

PresenceAbsenceModel_MedianGC::PresenceAbsenceModel_MedianGC( const ChipBase* i_pChip ) : m_pChip (i_pChip)
{
	
	// First store the background probes from BackgroundDomain
	// by the gc content.  This enables easy determination 
	// of empirical p-values.  
	this->SetGCBinBackgroundData();
		
	
	// For each tc
	// for each core probe
	// compute a p-value or score
	// end
	// compute a p-value or score for the tc
	// output the value, num core probes, p-value
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
			vector<double> gcBinPValues;
			
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
					double pVal = this->GetGCPvalue( pProbe );
					gcBinPValues.push_back( pVal );
				}
			}
			
			int probe_cnt = static_cast<int>( gcBinPValues.size() );
			if ( probe_cnt > 0 ) {
				double test_statistic = 0;	// Changed Affy's code so that the product 
											// of the p-values will not -> 0 for transcript 
											// clusters with large number of probes.  
				vector<double>::iterator gcBinPValItBegin = gcBinPValues.begin();
				vector<double>::iterator gcBinPValItEnd = gcBinPValues.end();
				for (; gcBinPValItBegin != gcBinPValItEnd; ++gcBinPValItBegin) {
					test_statistic += -2 * log((*gcBinPValItBegin));
				}
				double pValueGCBin = affxstat::chisqrprob(2*probe_cnt,(float)test_statistic);
				output << tcId << sep << test_statistic << sep << probe_cnt << sep << pValueGCBin << endl;
			}
		}
	} else {
		Parameters::LOGFILE << "Error:  Cannot open presence absence output file." << endl;
	}
}

PresenceAbsenceModel_MedianGC::~PresenceAbsenceModel_MedianGC()
{	
}

void PresenceAbsenceModel_MedianGC::SetGCBinBackgroundData( ) {
	for ( int gcbin = 0; gcbin < 26; ++gcbin ) {
		// First create a vector to store background probe intensities having gc content gcbin
		vector<unsigned short> myVec;
		this->m_BackgroundData.push_back( myVec );
	}	
	BackgroundDomain* pBackgroundDomain = BackgroundDomain::GetInstance();
	vector<ProbeBase*>::const_iterator itBegin = pBackgroundDomain->m_Probes.begin();
	vector<ProbeBase*>::const_iterator itEnd = pBackgroundDomain->m_Probes.end();
	for (; itBegin != itEnd; ++itBegin) {
		ProbeBase* pProbe = *itBegin;
		unsigned short inten = m_pChip->GetUnnormalizedIntensity( pProbe );	
		int gcbin = pProbe->GetGcCount();
		this->m_BackgroundData[gcbin].push_back( inten );
	}
	
	
	// For each vector of background probes of a given gc content, sort the values
	int size = static_cast<int>( this->m_BackgroundData.size() );
	for (int gcbin=0; gcbin < size; ++gcbin) {
		sort( this->m_BackgroundData[gcbin].begin(),this->m_BackgroundData[gcbin].end() );
	}
	
	return;
}


double PresenceAbsenceModel_MedianGC::GetGCPvalue( ProbeBase* i_pProbe ) {

	int gcCount = i_pProbe->GetGcCount();
	unsigned short inten = m_pChip->GetUnnormalizedIntensity( i_pProbe );
	double pVal;
	int numBackgroundProbes = m_BackgroundData[ gcCount ].size();
	// Handle cases where the number of background probes is less than 10
	if ( numBackgroundProbes < 10 ) { // This may occur at very low GC counts or very high GC counts
		bool searchHigher;			// We simply search the next up or next down bin.  
		if ( gcCount < 13 ) {		// This has the potential for errors if there are low numbers of background probes.  
			searchHigher = true;
		} else {
			searchHigher = false;
		}
		while ( numBackgroundProbes < 10 && gcCount > 0 && gcCount < 26 ) {
			if (searchHigher ) {
				numBackgroundProbes = m_BackgroundData[ (++gcCount) ].size(); // adjust the value of gcCount
			} else {
				numBackgroundProbes = m_BackgroundData[ (--gcCount) ].size();
			}
		}
	}
	

	// If the number of background probes is still too small, handle those cases here.  
	// Should not happen!  
	if ( numBackgroundProbes == 0 ) {
		pVal = 1.0;
	} else if ( numBackgroundProbes < 10 && numBackgroundProbes > 0 ) {
		unsigned short theval = m_BackgroundData[gcCount][0];
		if ( inten < theval) {
			pVal=0.0; // perl dabg returns 1; we invert
		} else if (inten > theval) {
			pVal=1.0; // perl dabg returns 0; we invert
		} else {
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
	// invert to get the chance of being greater.
	pVal=(1.0-pVal);
	return pVal;
}

