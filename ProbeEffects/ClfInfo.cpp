#include "ClfInfo.h"

ClfInfo::ClfInfo( int i_X, int i_Y, bool i_IsCrossHyb, string i_ChInfo ) : m_X( i_X), m_Y( i_Y ), 
						m_IsCrossHyb( i_IsCrossHyb ), m_ChInfo( i_ChInfo )
{								
}

ClfInfo::~ClfInfo() {
}

int ClfInfo::GetX() const
{
	return this->m_X;
}
int ClfInfo::GetY() const
{
	return this->m_Y;
}
bool ClfInfo::GetIsCrossHyb() const
{
	return this->m_IsCrossHyb;
}
string ClfInfo::GetChInfo() const
{
	return this->m_ChInfo;
}




