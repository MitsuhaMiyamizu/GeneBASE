#ifndef __BackgroundDomain_h_included__
#define __BackgroundDomain_h_included__


//***********
// BackgroundDomain
// Loads probes to be used for fitting a background model.  
// The probes are stored here in a vector.  
// Also outputs the background model header.  
//***********

#include "stdafx.h"
#include "ProbeBase.h"

class BackgroundDomain
{
public:
	static BackgroundDomain* GetInstance();
	~BackgroundDomain(void);
	
	void OutputBackgroundModelHeader( ofstream& i_Outfile );
	vector<ProbeBase*> m_Probes; // probes for training a background model
private:
	BackgroundDomain();
	bool LoadProbes();
	
	static bool m_InstanceFlag;
	static BackgroundDomain* m_pTheOne;
};



#endif

