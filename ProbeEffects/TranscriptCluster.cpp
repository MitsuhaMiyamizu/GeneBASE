#include "stdafx.h"
#include "TranscriptCluster.h"
#include "Probeset.h"
#include "Parameters.h"
#include "chip_data.h"
#include <valarray>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_rng.h>

TranscriptCluster::TranscriptCluster(const unsigned long& i_TranscriptClusterId): m_TranscriptClusterId( i_TranscriptClusterId )
{
}

TranscriptCluster::~TranscriptCluster(void)
{
	multimap< AnnotationConf, Probeset* >::iterator itBegin = this->m_ProbesetsByAnnotation.begin();
	multimap< AnnotationConf, Probeset* >::iterator itEnd = this->m_ProbesetsByAnnotation.end();
	for (; itBegin != itEnd; ++itBegin ) {
		Probeset* pProbeset = itBegin->second;
		delete pProbeset;
		pProbeset = NULL;
	}
	this->m_ProbesetsByAnnotation.clear();
}

void TranscriptCluster::AddProbeset(Probeset* i_pProbeset, const AnnotationConf& i_AnnotationConf)  {
	this->m_ProbesetsByAnnotation.insert( make_pair(i_AnnotationConf, i_pProbeset) );
}

pair<multimap<AnnotationConf, Probeset*>::iterator, multimap<AnnotationConf, Probeset*>::iterator > TranscriptCluster::GetProbesetIterators( const AnnotationConf& i_Level )
{
	return make_pair( this->m_ProbesetsByAnnotation.lower_bound( i_Level ), this->m_ProbesetsByAnnotation.upper_bound(i_Level) );
}



int TranscriptCluster::GetNumProbes( const AnnotationConf& i_Annot, const bool& i_RemoveExcludeProbes )
{
	int numProbes = 0;
	for ( AnnotationConf level = i_Annot; level >= CORE; ) {
		pair<multimap<AnnotationConf, Probeset*>::iterator, multimap<AnnotationConf, Probeset*>::iterator> itProbesetPair = this->GetProbesetIterators( level );
		level = static_cast<AnnotationConf>( static_cast<int>(level) - 1 );
		multimap<AnnotationConf, Probeset*>::iterator itProbesetBegin = itProbesetPair.first;
		multimap<AnnotationConf, Probeset*>::iterator itProbesetEnd = itProbesetPair.second;
		for (; itProbesetBegin != itProbesetEnd; ++itProbesetBegin ) {
			Probeset* pProbeset = itProbesetBegin->second;
			numProbes += pProbeset->GetNumProbes( i_RemoveExcludeProbes );
		}
	}
	return numProbes;
}

bool TranscriptCluster::CanDeconvolve()
{
	pair<multimap<AnnotationConf, Probeset*>::iterator, multimap<AnnotationConf, Probeset*>::iterator> itProbesetPair = this->GetProbesetIterators( CORE );
	multimap<AnnotationConf, Probeset*>::iterator itProbesetBegin = itProbesetPair.first;
	multimap<AnnotationConf, Probeset*>::iterator itProbesetEnd = itProbesetPair.second;
	for (; itProbesetBegin != itProbesetEnd; ++itProbesetBegin ) {
		Probeset* pProbeset = itProbesetBegin->second;
		pair<map<unsigned long, Probe*>::iterator, map<unsigned long, Probe*>::iterator> itProbePair = pProbeset->GetProbeIterators();
		map<unsigned long, Probe*>::iterator itProbeBegin = itProbePair.first;
		map<unsigned long, Probe*>::iterator itProbeEnd = itProbePair.second;
		for (; itProbeBegin != itProbeEnd; ++itProbeBegin ) {
			Probe* pProbe = itProbeBegin->second;
			int numOffTargets = pProbe->GetNumChTcIds();
			if ( numOffTargets == 1 ) {
				return true; // at least one probe matches a unique off-target transcript
			}
		}
	}
	return false;
}


pair< multimap<AnnotationConf, Probeset*>::iterator, multimap<AnnotationConf, Probeset*>::iterator > TranscriptCluster::GetProbesetIterators()
{
	return make_pair( this->m_ProbesetsByAnnotation.begin(), this->m_ProbesetsByAnnotation.end() );
}

unsigned long TranscriptCluster::GetId() const
{
	return this->m_TranscriptClusterId;
}

Probe* TranscriptCluster::GetProbe(const string& i_ProbeId)
{
	multimap< AnnotationConf, Probeset* >::iterator itBegin = this->m_ProbesetsByAnnotation.begin();
	multimap< AnnotationConf, Probeset* >::iterator itEnd = this->m_ProbesetsByAnnotation.end();
	for (; itBegin != itEnd; ++itBegin) {
		Probeset* pProbeset = itBegin->second;
		Probe* pProbe = pProbeset->GetProbe( i_ProbeId );
		if ( pProbe != NULL ) {
			return pProbe;
		}
	}

	return NULL; // If the probe has not been found, return NULL
}

AnnotationConf TranscriptCluster::GetHighestAnnotationProbesetType() const
{
	AnnotationConf annot;
	for (annot = CORE; annot <= NONE;  ) {
		multimap< AnnotationConf, Probeset*>::const_iterator itFind = this->m_ProbesetsByAnnotation.find( annot );
		if ( itFind != this->m_ProbesetsByAnnotation.end() ) {
			return annot;
		}
		annot = static_cast<AnnotationConf>(static_cast<int>(annot) + 1);// increment annot
	}
	return annot; // equals NONE here
}





