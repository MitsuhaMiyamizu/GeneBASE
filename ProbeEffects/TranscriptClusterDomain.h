#ifndef __TranscriptClusterDomain_h_included__
#define __TranscriptClusterDomain_h_included__

#include <string>
#include <map>
#include <vector>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
using namespace std;
#include "TranscriptCluster.h"
#include "ProbeInfo.h"
#include "ClfInfo.h"


class TranscriptClusterDomain
{
public:
	static TranscriptClusterDomain* GetInstance();
	~TranscriptClusterDomain(void);
	pair< hash_map<unsigned long, TranscriptCluster*>::iterator, hash_map<unsigned long, TranscriptCluster*>::iterator> GetTranscriptClusterIterators();
	TranscriptCluster* GetTranscriptCluster( const string& i_TranscriptClusterId );
	TranscriptCluster* GetTranscriptCluster( const unsigned int i_TranscriptClusterId );
	Probeset* GetProbeset( const string& i_ProbesetId );
	bool IsTranscriptCluster( const string& i_TranscriptClusterId ) const;
	Probeset* GetRandomLowConfProbeset( gsl_vector_int* i_pVec );
	int GetNumGenomicBgdProbes( int i_Key ) const;
	int GetNumAntigenomicBgdProbes( int i_Key ) const;
	int GetNumGenomicBgdProbes() const;
	int GetNumAntigenomicBgdProbes() const;
	pair <multimap<int, Probe*>::iterator, multimap<int, Probe*>::iterator> GetGenomicBgdIterators( int i_Key );
	pair <multimap<int, Probe*>::iterator, multimap<int, Probe*>::iterator> GetAntigenomicBgdIterators( int i_Key );
	pair< multimap<int, Probe*>::iterator, multimap<int, Probe*>::iterator > GetGenomicBgdIterators();
	pair< multimap<int, Probe*>::iterator, multimap<int, Probe*>::iterator > GetAntigenomicBgdIterators();

	bool OutputCoreProbes( const string& i_File );
	bool OutputFullProbes( const string& i_File );
	bool OutputBackgroundProbes( const string& i_File ); 
	bool OutputExtendedProbes( const string& i_File );
	bool OutputRandomProbeIntensities( const string& i_ProbeFile, const string& i_Outfile );

private:
	TranscriptClusterDomain(void);
	static bool m_InstanceFlag;
	static TranscriptClusterDomain* m_pTheOne;

	bool LoadProbeInfo();
	bool LoadTranscriptClusters();
	bool ReadProbesets();
	bool AnnotateProbesetsExonArray();
	bool AnnotateProbesetsGeneArray();
	bool ReadClf();
	bool AnnotateCrossHybExonArray();
	hash_map<unsigned long, TranscriptCluster*> m_TranscriptClusters;
	hash_map<unsigned long, Probeset*> m_Probesets; // Probesets must be stored by id to enable look up when assigning to transcript clusters.  
	multimap< std::vector<int>, Probeset*> m_LowConfProbesets; // Stores full and lower annotational confidence probesets.  the vector is the gc content
	map< unsigned long, ClfInfo > m_ClfInfo; // Stores clf info:  probe id, x and y position.  Later erased.  

	multimap<int, Probe*> m_GenomicBgdProbes;
	multimap<int, Probe*> m_AntigenomicBgdProbes;

	hash_multimap<unsigned long, ProbeInfo*> m_ProbeInfo;

};

#endif
