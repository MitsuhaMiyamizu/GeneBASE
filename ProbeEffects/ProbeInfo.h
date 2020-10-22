#ifndef __ProbeInfo_h_included__
#define __ProbeInfo_h_included__

#include "stdafx.h"

struct ProbeInfo {
public:
	unsigned long m_ProbeId;
	unsigned long m_TranscriptClusterId;
	string m_Level;
	string m_Sequence;

	ProbeInfo (unsigned long i_ProbeId, unsigned long i_TcId, 
		string i_Level, string i_Sequence) : m_ProbeId (i_ProbeId), m_TranscriptClusterId( i_TcId), 
												m_Level( i_Level ), m_Sequence( i_Sequence ) {};
	void Annotate (unsigned long i_TcId, string i_Level) { m_TranscriptClusterId = i_TcId; m_Level = i_Level; };
};




#endif


