#ifndef __TranscriptCluster_h_included__
#define __TranscriptCluster_h_included__

#include "Parameters.h"
#include "Probeset.h"
#include <map>
using namespace std;
#include <gsl/gsl_matrix.h>

class TranscriptCluster
{
	friend class TranscriptClusterDomain;
	friend class SimulationDomain;

public:
	TranscriptCluster( const unsigned long& i_TranscriptClusterId );
	~TranscriptCluster(void);
	pair<multimap<AnnotationConf, Probeset*>::iterator, multimap<AnnotationConf, Probeset*>::iterator > GetProbesetIterators( const AnnotationConf& i_Level );
	pair< multimap<AnnotationConf, Probeset*>::iterator, multimap<AnnotationConf, Probeset*>::iterator > GetProbesetIterators();
	int GetNumProbes( const AnnotationConf& i_Annot, const bool& i_RemoveExcludeProbes = false );
	bool CanDeconvolve();
	unsigned long GetId() const;
	Probe* GetProbe(const string& i_ProbeId);

	AnnotationConf GetHighestAnnotationProbesetType() const;
	

private:
	void AddProbeset( Probeset* i_pProbeset, const AnnotationConf& i_Annot );
	unsigned long m_TranscriptClusterId;
	multimap< AnnotationConf, Probeset* > m_ProbesetsByAnnotation;

};

#endif

