#include "stdafx.h"
#include "Probeset.h"
#include <string>
#include <map>
using namespace std;

Probeset::Probeset( const unsigned long& i_ProbesetId, const string& i_ProbesetType )
{
	m_ProbesetId = i_ProbesetId;
	this->m_IsGenomicBackgroundProbeset = false;
	this->m_IsAntigenomicBackgroundProbeset = false;
	
	if ( i_ProbesetType == "control->bgp->genomic" ) {
		this->m_IsGenomicBackgroundProbeset = true;
	} else if ( i_ProbesetType == "control->bgp->antigenomic" ) {
		this->m_IsAntigenomicBackgroundProbeset = true;
	} 
}

Probeset::~Probeset(void)
{
	// Delete the probes one by one
	for ( map<unsigned long, Probe*>::iterator it = m_Probes.begin(); it != m_Probes.end(); it++)
	{
		Probe* pProbe = it->second;
		delete pProbe;
		pProbe = NULL;
	}
	
}

bool Probeset::AddProbe( Probe* i_pProbe )
{
	unsigned long probeId = i_pProbe->GetProbeId();
	map<unsigned long, Probe*>::iterator it = m_Probes.find( probeId );
	if (it == m_Probes.end() )
	{
		this->m_Probes.insert( make_pair(probeId, i_pProbe) );
		return true;
	} else {
		return false;
	}
}

bool Probeset::AddProbe(const unsigned long& i_ProbeId, 
						const string& i_Sequence)
{
	map<unsigned long, Probe*>::iterator it = m_Probes.find( i_ProbeId );
	if (it == m_Probes.end() )
	{
		string sequence = Parameters::SWAP_ORDER( i_Sequence ); // Exon array probes have sequences 3' to 5'
																// We swap the order to follow the 5' to 3' convention
		Probe* pNewProbe = new Probe(i_ProbeId, sequence);
		this->m_Probes.insert( make_pair(i_ProbeId, pNewProbe) );
		return true;
	} else
	{
		// Error! Trying to insert an already existing probe!
		// This error message should be output somewhere!  
		return false;
	}
}
void Probeset::Annotate(const string& i_GeneAssignment, const string& i_Start, const string& i_Stop, const string& i_Level, const string& i_TranscriptClusterId)
{
	//
	// TBD:  Add code
	// 
}

void Probeset::GetGCVec( std::vector<int>& i_GCVec )
{
	map<unsigned long, Probe*>::iterator probeItBegin = this->m_Probes.begin();
	map<unsigned long, Probe*>::iterator probeItEnd = this->m_Probes.end();
	for(; probeItBegin != probeItEnd; ++probeItBegin ) {
		Probe* pProbe = probeItBegin->second;
		i_GCVec.push_back( pProbe->GetGcCount() );
	}
	return;
}

Probe* Probeset::GetProbe(const string& i_ProbeId)
{
	Probe* pProbe = NULL;
	unsigned long probeIdVal = strtoul( i_ProbeId.c_str(), NULL, 10 );
	map<unsigned long, Probe*>::iterator it = m_Probes.find( probeIdVal );
	if (it != m_Probes.end() )
	{
		pProbe = it->second;
	}
	return pProbe;
}

bool Probeset::IsGenomicBackgroundProbeset() const
{
	return this->m_IsGenomicBackgroundProbeset;
}
bool Probeset::IsAntigenomicBackgroundProbeset() const
{
	return this->m_IsAntigenomicBackgroundProbeset;
}

int Probeset::GetNumProbes( const bool& i_RemoveExcludeProbes ) const
{
	int numProbes = 0;
	if ( i_RemoveExcludeProbes == false ) {
		numProbes = static_cast<int>(this->m_Probes.size());
	} else {
		map<unsigned long, Probe*>::const_iterator itBeg = this->m_Probes.begin();
		map<unsigned long, Probe*>::const_iterator itEnd = this->m_Probes.end();
		for (; itBeg != itEnd; ++itBeg) {
			Probe* pProbe = itBeg->second;
			if ( pProbe->IsCrossHyb() == false) {
				++numProbes;
			}
		}
	}
	return numProbes;
}

pair<map<unsigned long, Probe*>::iterator, map<unsigned long, Probe*>::iterator> Probeset::GetProbeIterators()
{
	return make_pair( this->m_Probes.begin(), this->m_Probes.end() );	
}

unsigned long Probeset::GetProbesetId() const
{
	return this->m_ProbesetId;
}
	




