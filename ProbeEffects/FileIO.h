/////////////////////////////////////////////////////////////////
//
// Copyright (C) 2004 Affymetrix, Inc.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation; either version 2.1 of the License,
// or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
// for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library; if not, write to the Free Software Foundation, Inc.,
// 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA 
//
/////////////////////////////////////////////////////////////////

#if !defined(AFX_FILEIO_H__8A5CA65E_C2E0_4C36_9B60_30C2020D5DC0__INCLUDED_)
#define AFX_FILEIO_H__8A5CA65E_C2E0_4C36_9B60_30C2020D5DC0__INCLUDED_

//////////////////////////////////////////////////////////////////////

//#ifndef WIN32 // ie if the system uses __PPC__
//#define IS_BIG_ENDIAN ; // Biogibbs has PPC and hence uses big endian
//							// Windows and Macs do not have PPC and hence use little endian.  
//#endif

#ifdef __sparc__ 
#define IS_BIG_ENDIAN 1
#endif
#ifdef __PPC__ // on linux
#define IS_BIG_ENDIAN 1
#endif
#ifdef __POWERPC__ // on darwin
#define IS_BIG_ENDIAN 1
#endif
#ifdef __BIG_ENDIAN__ // we are being told
#undef  IS_BIG_ENDIAN // avoid the warning
#define IS_BIG_ENDIAN 1
#endif


#ifndef IFSTREAM
/*! Input STL stream object */
#define IFSTREAM std::istream
#endif 

#include <fstream>
#include <string>
#include "affy-base-types.h"

//////////////////////////////////////////////////////////////////////

#define FLOAT_SIZE sizeof(float)
#define INT_SIZE sizeof(int)
#define UINT_SIZE sizeof(unsigned int)
#define SHORT_SIZE sizeof(short)
#define USHORT_SIZE sizeof(unsigned short)
#define CHAR_SIZE sizeof(char)
#define UCHAR_SIZE sizeof(unsigned char)

#define ULONG_SIZE sizeof(unsigned long)
#define LONG_SIZE sizeof(long)
/*! The size of a 32 bit integer */
#define INT32_SIZE  sizeof(int32_t )

//////////////////////////////////////////////////////////////////////

/*! Shuffle operator for 32 bit values
 * @param x The value to shuffle
 * @return The shuffled value
 */
#define affy_swap32(x) \
     ((((x) & 0xff000000) >> 24) | (((x) & 0x00ff0000) >>  8) | \
      (((x) & 0x0000ff00) <<  8) | (((x) & 0x000000ff) << 24))

//////////////////////////////////////////////////////////////////////

void SwapBytes(char *src, char *dest, int size);

//////////////////////////////////////////////////////////////////////

void READ_INT(std::ifstream &instr, int &value);
void READ_UINT(std::ifstream &instr, unsigned int &value);
void READ_SHORT(std::ifstream &instr, short &value);
void READ_USHORT(std::ifstream &instr, unsigned short &value);
void READ_CHAR(std::ifstream &instr, char &value);
void READ_UCHAR(std::ifstream &instr, unsigned char &value);
void READ_BOOL(std::ifstream &instr, char &value);
void READ_FLOAT(std::ifstream &instr, float &value);
void READ_STRING(std::ifstream &instr, char * &value);
void READ_FIXED_STRING(std::ifstream &instr, char *value, int len);
void GET_NEXT_LINE(std::ifstream &instr, char *line, int len);
void READ_INT(const char* instr, int &value);

void READ_ULONG(std::ifstream &instr, unsigned long &value);
void READ_LONG(std::ifstream &instr, long &value);

//////////////////////////////////////////////////////////////////////

void ReadInt(std::ifstream &instr, int &ival);
void ReadUChar(std::ifstream &instr, unsigned char &ucval);
void ReadUShort(std::ifstream &instr, unsigned short &usval);
void ReadUInt(std::ifstream &instr, unsigned int &uval);
void ReadFixedString(std::ifstream &instr, std::string &s, int len);
/*! Reads a fixed length string from a file.
 * @param instr The input file stream
 * @param str The returned string value.
 * @param len The length of the string to read
 */
void ReadFixedString(IFSTREAM& instr, std::string& str, uint32_t len);

void ReadFixedUChar(std::ifstream &instr, unsigned char *s, int len);
void ReadUIntLenString(std::ifstream &instr, std::string &s);
void ReadString(std::ifstream &instr, std::string &s);
void ReadFloat(std::ifstream &instr, float &fval);
void ReadFloatFromOldBPMAP(std::ifstream &instr, float &fval);
void ReadUInt(const char *instr, unsigned int &uval);
void ReadFloat(const char *instr, float &fval);
void ReadFloatFromOldBPMAP(const char *instr, float &fval);
void ReadInt(const char *instr, int &ival);



/*! This function is for older BPMAP file with incorrectly written float's
 * @param instr The input file stream.
 * @param fval The returned float
 */
void ReadFloatFromOldBPMAP_N(IFSTREAM &instr, float &fval);

/*! Reads a float value from a big endian file
 * @param instr The input file stream
 * @param val The returned value
 */
void ReadFloat_N(IFSTREAM& instr,float& val);

/*! Reads an unsigned 32 bit value from a big endian file
 * @param instr The input file stream
 * @param val The returned value
 */
void ReadUInt32_N(IFSTREAM& instr,uint32_t& val);

/*! Reads a string frome a big endian file where the string length is stored as an unsigned integer.
 * @param instr The input file stream
 * @param s The returned value
 */
void ReadUIntLenString_N(IFSTREAM &instr, std::string &s);

/*! Reads a fixed length string from a file.
 * @param instr The input file stream
 * @param str The returned string value.
 * @param len The length of the string to read.
 */
void ReadFixedUCString(IFSTREAM& instr, unsigned char* str, uint32_t len);

/*! Reads an unsigned 8 bit value from a file
 * @param instr The input file stream
 * @param val The returned value
 */
void ReadUInt8(IFSTREAM& instr,uint8_t& val);

/*! Gets a float from a big endian data stream.
 * @param ptr A pointer to a big endian data stream
 * @return The host ordered value
 */
float    MmGetFloat_N(float*     ptr);

/*! This function is for older BPMAP file with incorrectly written float's
 * @param ptr A pointer to a little endian data stream
 * @return The host ordered value
 */
float	MmGetFloatFromOldBPMAP_N(float *ptr);

/*! Gets a 32 bit unsigned int from a big endian data stream.
 * @param ptr A pointer to a big endian data stream
 * @return The host ordered value
 */
uint32_t MmGetUInt32_N(uint32_t* ptr);


//////////////////////////////////////////////////////////////////////

#endif // !defined(AFX_FILEIO_H__8A5CA65E_C2E0_4C36_9B60_30C2020D5DC0__INCLUDED_)
