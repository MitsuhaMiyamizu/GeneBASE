/*
 *  BackgroundDomain.cpp
 *  ProbeEffects2.0
 *
 *  Created by Karen Kapur on 7/9/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include "BackgroundDomain.h"
#include "TranscriptClusterDomain.h"
#include "Parameters.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


bool BackgroundDomain::m_InstanceFlag = false;
BackgroundDomain* BackgroundDomain::m_pTheOne = NULL;

BackgroundDomain::BackgroundDomain( )
{
	bool success = this->LoadProbes();
	if ( success == false ) {
		Parameters::LOGFILE << "Error:  Cannot load probes in BackgroundDomain" << endl;
		exit(1);
	}

}

BackgroundDomain::~BackgroundDomain(void)
{
	// Clear probes only
	// Does not delete probes because 
	// probes are owned by TranscriptClusterDomain
	this->m_Probes.clear(); 
}

BackgroundDomain* BackgroundDomain::GetInstance( ) {
	if ( !m_InstanceFlag )
	{
		m_pTheOne = new BackgroundDomain( );
		m_InstanceFlag = true;
		return m_pTheOne;
	}
	else {
		return m_pTheOne;
	}
}

bool BackgroundDomain::LoadProbes() {
	Parameters* pParameters = Parameters::GetInstance();
	string trainingProbeType = pParameters->GetMATTrainingProbeType();
	int numProbes = pParameters->GetNumProbesMAT();
	
	try{
		TranscriptClusterDomain* pTranscriptClusterDomain = TranscriptClusterDomain::GetInstance();
		
		bool loadBackground = false;
		AnnotationConf annot = NONE;
		if ( trainingProbeType == "full" ) {
			annot = FULL;
		} else if ( trainingProbeType == "extended" ) {
			annot = EXTENDED;
		} else if ( trainingProbeType == "core" ) {
			annot = CORE;
		} else if ( trainingProbeType == "background" ) {
			loadBackground = true;
		}
		
		if ( loadBackground == true ) {
			pair< multimap<int, Probe*>::iterator, multimap<int, Probe*>::iterator > itBkgProbeItPair = pTranscriptClusterDomain->GetGenomicBgdIterators();
			multimap<int, Probe*>::iterator itBkgProbeItBegin = itBkgProbeItPair.first;
			multimap<int, Probe*>::iterator itBkgProbeItEnd = itBkgProbeItPair.second;
			for (; itBkgProbeItBegin != itBkgProbeItEnd; ++itBkgProbeItBegin ) {
				// Copy info into a new ChipChipProbe object
				// Put the new ChipChipProbe into this->m_BackgroundProbes
				Probe* pProbe = itBkgProbeItBegin->second;
				//ChipChipProbe* pChipChipProbe = new ChipChipProbe( pProbe->GetProbeId(), pProbe->GetSequence() );
				//pChipChipProbe->SetId( pProbe->GetProbeId() );
				//pChipChipProbe->CopyInfo( pProbe );
				this->m_Probes.push_back( pProbe ); //this->m_Probes.push_back( pChipChipProbe );
			}
			
			itBkgProbeItPair = pTranscriptClusterDomain->GetAntigenomicBgdIterators();
			itBkgProbeItBegin = itBkgProbeItPair.first;
			itBkgProbeItEnd = itBkgProbeItPair.second;
			for (; itBkgProbeItBegin != itBkgProbeItEnd; ++itBkgProbeItBegin ) {
				Probe* pProbe = itBkgProbeItBegin->second;
				//ChipChipProbe* pChipChipProbe = new ChipChipProbe( pProbe->GetProbeId(), pProbe->GetSequence());
				//pChipChipProbe->CopyInfo(pProbe);
				this->m_Probes.push_back( pProbe );//this->m_Probes.push_back( pChipChipProbe );
			} 
		} else { // For each transcript cluster, add the probes of the given type into this->m_Probes
			pair< hash_map<unsigned long, TranscriptCluster*>::iterator, hash_map<unsigned long, TranscriptCluster*>::iterator > itTcPair = pTranscriptClusterDomain->GetTranscriptClusterIterators();
			hash_map<unsigned long, TranscriptCluster*>::iterator itTcBegin = itTcPair.first;
			hash_map<unsigned long, TranscriptCluster*>::iterator itTcEnd = itTcPair.second;
			for (; itTcBegin != itTcEnd; ++itTcBegin) {
				TranscriptCluster* pTranscriptCluster = itTcBegin->second;
				pair<multimap<AnnotationConf, Probeset*>::iterator, multimap<AnnotationConf, Probeset*>::iterator > itProbesetPair = 
					pTranscriptCluster->GetProbesetIterators( annot );
				multimap<AnnotationConf, Probeset*>::iterator itProbesetBegin = itProbesetPair.first;
				multimap<AnnotationConf, Probeset*>::iterator itProbesetEnd = itProbesetPair.second;
				for (; itProbesetBegin != itProbesetEnd; ++itProbesetBegin ) {
					Probeset* pProbeset = itProbesetBegin->second;
					// For each Probe
					// Copy info into a new ChipChipProbe object
					// Put the new ChipChipProbe into this->m_Probes
					pair< map<unsigned long, Probe*>::iterator, map<unsigned long, Probe*>::iterator > itProbePair = pProbeset->GetProbeIterators();
					map<unsigned long, Probe*>::iterator itProbeBegin = itProbePair.first;
					map<unsigned long, Probe*>::iterator itProbeEnd = itProbePair.second;
					for (; itProbeBegin != itProbeEnd; ++itProbeBegin ) {
						Probe* pProbe = itProbeBegin->second;
						//ChipChipProbe* pChipChipProbe = new ChipChipProbe( pProbe->GetProbeId(), pProbe->GetSequence() );
						//pair<int, int> posPair = pProbe->GetPosition();
						//pChipChipProbe->SetPosition( posPair.first, posPair.second );
						this->m_Probes.push_back( pProbe ); // m_Probes.push_back( pChipChipProbe );
					}
				}
			}
		}
		
		//***
		// Only keep min( m_Probes.size(), numProbes)
		// Determine which to keep by random permutation
		//***
		if ( this->m_Probes.size() <= numProbes ) {
			// do nothing
		} else {
			//fit the model using i_NumProbes, randomly selected.
			vector<ProbeBase*>* pSelectedProbes = new vector<ProbeBase*>;
			gsl_rng* pRng = pParameters->GetRng();
			int* pSource = new int[ numProbes ];
			for ( int i = 0; i < numProbes; ++i ) {
				pSource[i] = i;
			}
			int* pDest = new int[ numProbes ];
			gsl_ran_choose (pRng, pDest, numProbes, pSource, m_Probes.size(), sizeof(int) );
			// Next select the corresponding elements of this->m_Probes
			for ( int i = 0; i < numProbes; ++i ) {
				int val = pDest[i]; // the ith randomly chosen value
				ProbeBase* pProbe = this->m_Probes.operator[](val);
				pSelectedProbes->push_back( pProbe );
			}
			delete [] pSource;
			delete [] pDest;
			this->m_Probes.clear();
			vector<ProbeBase*>::iterator itBegin = pSelectedProbes->begin();
			vector<ProbeBase*>::iterator itEnd = pSelectedProbes->end();
			for (; itBegin != itEnd; ++itBegin ) {
				ProbeBase* pProbeBase = *itBegin;
				this->m_Probes.push_back( pProbeBase );
			}
			pSelectedProbes->clear();
			delete pSelectedProbes;
		}
	} catch (...) {
		return false;
	}
	return true;
}


void BackgroundDomain::OutputBackgroundModelHeader( ofstream& i_Outfile )
{
	Parameters* pParameters = Parameters::GetInstance();
	string backgroundMethod = pParameters->GetModelBackgroundMethod();
	string SEP = "\t";
	if ( backgroundMethod == "mat" ) {
		char nucs[4] = {'A','C','G','T' };
		i_Outfile << "alpha" << SEP;
		for ( int j = 0; j < 25; ++j ) {
			for ( int k = 0; k < 3; ++k ) {
				i_Outfile << "beta_" << (j+1) << nucs[k] << SEP;
			}
		}
		for ( int k = 0; k < 4; ++k ) {
			i_Outfile << "gamma_" << nucs[k] << SEP;
		}
		i_Outfile << "RSquared_LogScale" << SEP
			<< "RSquared_OriginalScale" << SEP
			<< "SigmaHat" << endl; 
		
	} else if ( backgroundMethod == "median_gc" ) {
		int maxGCContent = Parameters::PROBELENGTH;
		for ( int i = 0; i < maxGCContent; ++i ) {
			i_Outfile << i << SEP;
		}
		i_Outfile << endl;		
	} else if ( backgroundMethod == "none" ) {
		i_Outfile << "none" << endl;
	} else {
		// output nothing
	}
}



