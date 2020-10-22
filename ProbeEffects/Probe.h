#ifndef __Probe_h_included__
#define __Probe_h_included__

#include "stdafx.h"
#include "ProbeBase.h"
#include <string>
using namespace std;

class Probe : public ProbeBase
{
public:
	Probe(const unsigned long& i_Id,  
		const string& i_Sequence);
	virtual ~Probe(void);
	virtual void SetMATDesignVector( gsl_vector* o_pDesignVector, const bool& i_Binned = false, const bool& i_UseCopyNumber = false );
	
private:
	
};

#endif
