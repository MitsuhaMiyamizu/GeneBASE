
#ifndef __ClfInfo_h_included__
#define __ClfInfo_h_included__

#include "stdafx.h"

class ClfInfo {
public:	
	ClfInfo( int i_X, int i_Y, bool i_IsCrossHyb, string i_ChInfo );
	~ClfInfo(void);
	int GetX() const;
	int GetY() const;
	bool GetIsCrossHyb() const;
	string GetChInfo() const;

private:
	int m_X;
	int m_Y;
	bool m_IsCrossHyb;
	string m_ChInfo;
};



#endif


