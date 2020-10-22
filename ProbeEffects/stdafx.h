// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#ifndef __stdafx_h_included__
#define __stdafx_h_included__

#include <iostream>
#include <vector>
#include <set>
#include <string>
#include <map>
#include <list>
#include <fstream>
#include <sstream>

#ifdef WIN32 
#include <hash_map>
#else 
#include <ext/hash_map>
using namespace __gnu_cxx;
#endif

using namespace std;

#include "Parameters.h"
#define _CRT_SECURE_NO_DEPRECATE
#include <gsl/gsl_errno.h>


// TODO: reference additional headers your program requires here

#endif


