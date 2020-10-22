#include "stdafx.h"
#include "Probe.h"

Probe::Probe(const unsigned long& i_ProbeId,  
			 const string& i_Sequence) : ProbeBase(i_ProbeId, i_Sequence )
{
}

Probe::~Probe(void)
{
}

void Probe::SetMATDesignVector( gsl_vector* o_pDesignVector, const bool& i_Binned, const bool& i_UseCopyNumber )
{
	ProbeBase::SetMATDesignVector( o_pDesignVector, i_Binned, i_UseCopyNumber );
}



