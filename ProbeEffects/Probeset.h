#ifndef __Probeset_h_included__
#define __Probeset_h_included__

#include "Probe.h"
#include "stdafx.h"
#include <string>
#include <map>
using namespace std;

class Probeset
{
	friend class TranscriptClusterDomain;

public:
	Probeset( const unsigned long& i_ProbesetId, const string& i_ProbesetType );
	~Probeset(void);
	void GetGCVec( std::vector<int>& i_GCVec );
	Probe* GetProbe(const string& i_ProbeId);
	bool IsGenomicBackgroundProbeset() const;
	bool IsAntigenomicBackgroundProbeset() const;
	int GetNumProbes( const bool& i_RemoveExcludeProbes = false ) const;
	pair<map<unsigned long, Probe*>::iterator, map<unsigned long, Probe*>::iterator> GetProbeIterators();
	unsigned long GetProbesetId() const;
	

private:
	bool AddProbe( Probe* i_pProbe );
	bool AddProbe(const unsigned long& i_ProbeId, const string& i_Sequence);
	void Annotate(const string& i_GeneAssignment, const string& i_Start, const string& i_Stop, const string& i_Level, const string& i_TranscriptClusterId);


	unsigned long m_ProbesetId;
	bool m_IsGenomicBackgroundProbeset;
	bool m_IsAntigenomicBackgroundProbeset;

	map<unsigned long, Probe*> m_Probes;
};

#endif
