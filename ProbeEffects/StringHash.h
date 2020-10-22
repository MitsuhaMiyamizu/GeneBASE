#ifndef __StringHash_h_included__
#define __StringHash_h_included__

#include <string>
#include <ext/hash_map>

namespace __gnu_cxx
{
	template<> struct hash< std::string >
	{
		size_t operator()( const std::string& x ) const
	{
		return hash< const char* >()( x.c_str() );
	}
	};
}


#endif

/*
 *  StringHash.h
 *  ProbeEffects
 *
 *  Created by Karen Kapur on 7/27/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

