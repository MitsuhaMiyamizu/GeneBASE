/*
 *  PresenceAbsenceModel_MedianGC.h
 *  ProbeEffects2.0.1
 *
 *  Created by Karen Kapur on 8/11/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef __PresenceAbsenceModel_MedianGC_h_included__
#define __PresenceAbsenceModel_MedianGC_h_included__

#include "stdafx.h"
#include "ChipBase.h"

class PresenceAbsenceModel_MedianGC {
public:
	PresenceAbsenceModel_MedianGC( const ChipBase* i_pChip );
	~PresenceAbsenceModel_MedianGC();
	
private:
	void SetGCBinBackgroundData();
	double GetGCPvalue( ProbeBase* i_pProbe );
		
	const ChipBase* m_pChip;
	vector< vector<unsigned short> > m_BackgroundData;
	
	
	
};

#endif