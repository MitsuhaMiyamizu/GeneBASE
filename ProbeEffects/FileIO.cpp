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

#include "stdafx.h"
#include "FileIO.h"
#include "affy-base-types.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>


//////////////////////////////////////////////////////////////////////

#ifdef WIN32
#include <winsock2.h>
#else
#include <sys/types.h>
#include <netinet/in.h>
#include <inttypes.h>
#endif

//////////////////////////////////////////////////////////////////////

void SwapBytes(char *src, char *dest, int size)
{
	int i;
	for (i=0; i<size; i++)
		dest[size - i - 1] = src[i];
}

//////////////////////////////////////////////////////////////////////

void ReadUChar(std::ifstream &instr, unsigned char &ucval)
{
	instr.read((char *)&ucval, UCHAR_SIZE);
}

//////////////////////////////////////////////////////////////////////

void ReadUInt(std::ifstream &instr, unsigned int &uval)
{
	instr.read((char *)&uval, UINT_SIZE);
	uval = ntohl(uval);
}

//////////////////////////////////////////////////////////////////////

void ReadUShort(std::ifstream &instr, unsigned short &usval)
{
	instr.read((char *)&usval, USHORT_SIZE);
	usval = ntohs(usval);
}

//////////////////////////////////////////////////////////////////////

void ReadInt(std::ifstream &instr, int &ival)
{
	instr.read((char *)&ival, INT_SIZE);
	ival = ntohl(ival);
}

//////////////////////////////////////////////////////////////////////

void ReadUIntLenString(std::ifstream &instr, std::string &s)
{
	unsigned int len;
	ReadUInt(instr, len);
	ReadFixedString(instr, s, static_cast<int>(len));
}

//////////////////////////////////////////////////////////////////////

void ReadFixedUChar(std::ifstream &instr, unsigned char *s, int len)
{
	instr.read((char *)s, len);
}

//////////////////////////////////////////////////////////////////////

void ReadFixedString(std::ifstream &instr, std::string &s, int len)
{
	char *buf = new char [len+1];
	instr.read(buf, len);
	buf[len]=0;
	s = buf;
	delete[] buf;
}

//////////////////////////////////////////////////////////////////////

void ReadString(std::ifstream &instr, std::string &s)
{
	int len;
	ReadInt(instr, len);
	ReadFixedString(instr, s, len);
}

//////////////////////////////////////////////////////////////////////

void ReadFloat(std::ifstream &instr, float &fval)
{
	unsigned int uval;
	ReadUInt(instr, uval);
	memcpy(&fval, &uval, INT_SIZE);
}

//////////////////////////////////////////////////////////////////////

void ReadFloatFromOldBPMAP(std::ifstream &instr, float &fval)
{
#ifndef IS_BIG_ENDIAN
	instr.read((char *)&fval, FLOAT_SIZE);
	fval = (float) ntohl(fval);
#else
	float f1=0;
	instr.read((char *)&f1, FLOAT_SIZE);
	float f2=0;
	SwapBytes((char *)&f1, (char *)&f2, FLOAT_SIZE);
	int i2= static_cast<int>(f2);
	int i3;
	SwapBytes((char *)&i2, (char *)&i3, INT_SIZE);
	fval = i3;
#endif
}

//////////////////////////////////////////////////////////////////////

void ReadUInt(const char *instr, unsigned int &uval)
{
	uval = *(unsigned int *)instr;
	uval = ntohl(uval);
}

//////////////////////////////////////////////////////////////////////

void ReadFloat(const char *instr, float &fval)
{
	unsigned int uval;
	ReadUInt(instr, uval);
	memcpy(&fval, &uval, INT_SIZE);
}

//////////////////////////////////////////////////////////////////////

void ReadFloatFromOldBPMAP(const char *instr, float &fval)
{
#ifndef IS_BIG_ENDIAN
	fval = (float) ntohl(*(float *)instr);
#else
	float f1=*(float*) instr;
	float f2=0;
	SwapBytes((char *)&f1, (char *)&f2, FLOAT_SIZE);
	int i2= static_cast<int>(f2);
	int i3;
	SwapBytes((char *)&i2, (char *)&i3, INT_SIZE);
	fval = i3;
#endif
}

//////////////////////////////////////////////////////////////////////

void ReadInt(const char *instr, int &ival)
{
	ival = *(int *)instr;
	ival = ntohl(ival);
}

//////////////////////////////////////////////////////////////////////

void READ_INT(std::ifstream &instr, int &value)
{
	char str[INT_SIZE];
	instr.read(str, INT_SIZE);
#ifdef IS_BIG_ENDIAN
	char str2[INT_SIZE];
	SwapBytes(str, str2, INT_SIZE);
	memcpy(str, str2, INT_SIZE);
#endif
	memcpy(&value, str, INT_SIZE);
}

//////////////////////////////////////////////////////////////////////

void READ_INT(const char* instr, int &value)
{
	char str[INT_SIZE];
	memcpy(str, instr, INT_SIZE);
#ifdef IS_BIG_ENDIAN
	char str2[INT_SIZE];
	SwapBytes(str, str2, INT_SIZE);
	memcpy(str, str2, INT_SIZE);
#endif
	memcpy(&value, str, INT_SIZE);
}

//////////////////////////////////////////////////////////////////////

void READ_UINT(std::ifstream &instr, unsigned int &value)
{
	char str[UINT_SIZE];
	instr.read(str, UINT_SIZE);
#ifdef IS_BIG_ENDIAN
	char str2[UINT_SIZE];
	SwapBytes(str, str2, UINT_SIZE);
	memcpy(str, str2, UINT_SIZE);
#endif
	memcpy(&value, str, UINT_SIZE);
}

//////////////////////////////////////////////////////////////////////

void READ_SHORT(std::ifstream &instr, short &value)
{
	char str[SHORT_SIZE];
	instr.read(str, SHORT_SIZE);
#ifdef IS_BIG_ENDIAN
	char str2[SHORT_SIZE];
	SwapBytes(str, str2, SHORT_SIZE);
	memcpy(str, str2, SHORT_SIZE);
#endif
	memcpy(&value, str, SHORT_SIZE);
}

//////////////////////////////////////////////////////////////////////

void READ_USHORT(std::ifstream &instr, unsigned short &value)
{
	char str[USHORT_SIZE];
	instr.read(str, USHORT_SIZE);
#ifdef IS_BIG_ENDIAN
	char str2[USHORT_SIZE];
	SwapBytes(str, str2, USHORT_SIZE);
	memcpy(str, str2, USHORT_SIZE);
#endif
	memcpy(&value, str, USHORT_SIZE);
}

//////////////////////////////////////////////////////////////////////

void READ_CHAR(std::ifstream &instr, char &value)
{
	instr.read(&value, CHAR_SIZE);
}

//////////////////////////////////////////////////////////////////////

void READ_UCHAR(std::ifstream &instr, unsigned char &value)
{
	instr.read((char *) &value, UCHAR_SIZE);
}

//////////////////////////////////////////////////////////////////////

void READ_BOOL(std::ifstream &instr, char &value)
{
	instr.read(&value, CHAR_SIZE);
}

//////////////////////////////////////////////////////////////////////

void READ_FLOAT(std::ifstream &instr, float &value)
{
	char str[FLOAT_SIZE];
	instr.read(str, FLOAT_SIZE);
#ifdef IS_BIG_ENDIAN
	char str2[FLOAT_SIZE];
	SwapBytes(str, str2, FLOAT_SIZE);
	memcpy(str, str2, FLOAT_SIZE);
#endif
	memcpy(&value, str, FLOAT_SIZE);
}

//////////////////////////////////////////////////////////////////////

void READ_STRING(std::ifstream &instr, char * &value)
{
	int len;
	READ_INT(instr, len);
	value = new char[len+1];
	if (len > 0)
		instr.read(value, len);
	value[len] = '\0';
}

//////////////////////////////////////////////////////////////////////

void READ_FIXED_STRING(std::ifstream &instr, char *value, int len)
{
	instr.read(value, len);
}

//////////////////////////////////////////////////////////////////////

void GET_NEXT_LINE(std::ifstream &instr, char *line, int len)
{
	strcpy(line, "");
	while (!instr.eof())
	{
		instr.getline(line, len);
		if (strlen(line) > 0)
		{
			if (line[strlen(line)-1] == '\r')
				line[strlen(line)-1] = '\0';
			if (strlen(line) > 0)
				return;
		}
	}
}

//////////////////////////////////////////////////////////////////////

void READ_ULONG(std::ifstream &instr, unsigned long &value)
{
	char str[ULONG_SIZE];
	instr.read(str, ULONG_SIZE);
#ifdef IS_BIG_ENDIAN
	char str2[ULONG_SIZE];
	SwapBytes(str, str2, ULONG_SIZE);
	memcpy(str, str2, ULONG_SIZE);
#endif
	memcpy(&value, str, ULONG_SIZE);
}

//////////////////////////////////////////////////////////////////////

void READ_LONG(std::ifstream &instr, long &value)
{
	char str[LONG_SIZE];
	instr.read(str, LONG_SIZE);
#ifdef IS_BIG_ENDIAN
	char str2[LONG_SIZE];
	SwapBytes(str, str2, LONG_SIZE);
	memcpy(str, str2, LONG_SIZE);
#endif
	memcpy(&value, str, LONG_SIZE);
}

//////////////////////////////////////////////////////////////////////

void
ReadFloatFromOldBPMAP_N(IFSTREAM &instr, float &fval)
{
#ifndef IS_BIG_ENDIAN
	instr.read((char *)&fval, FLOAT_SIZE);
	fval = (float) ntohl((uint32_t)fval);
#else
	int i1=0;
	instr.read((char *)&i1, INT32_SIZE);
 	i1=affy_swap32(i1);
	fval=(float)affy_swap32((uint32_t)(*(float *)&i1));
#endif
}

void
ReadFloat_N(IFSTREAM& instr,float& val) 
{
  ReadUInt32_N(instr,(uint32_t&)val);
}

void
ReadUInt32_N(IFSTREAM& instr, uint32_t& val) 
{
  uint32_t v=0;
  instr.read((char*)&v,sizeof(v));
  val=ntohl(v);
}

void
ReadUIntLenString_N(IFSTREAM& instr, std::string &s)
{
	uint32_t len;
	ReadUInt32_N(instr, len);
	ReadFixedString(instr, s, len);
}

// The data file may or may not have the null terminator as
// part of the length. Using a STD::STRING function as the buffer to
// read directly into may result in the string with a length one more
// than is intended. The CHAR * buffer will contain an extra null that
// will get removed when copying to the output STD::STRING.
void
ReadFixedString(IFSTREAM& instr, string& str, uint32_t len)
{
  char *s = new char[len+1];
  instr.read(s,len);
  s[len]=0;
  str = s;
  delete[] s;
}

void
ReadFixedUCString(IFSTREAM& instr, unsigned char* str, uint32_t len)
{
  instr.read((char *)str,len);
}

// No byte swapping needed.
void
ReadUInt8(IFSTREAM& instr, uint8_t& val) 
{
  instr.read((char*)&val,sizeof(val));
}

float
MmGetFloat_N(float* ptr)
{
  uint32_t v=ntohl(*(int*)ptr);
  return *((float*)&v);
}

float
MmGetFloatFromOldBPMAP_N(float *ptr)
{
	float fval;
#ifndef IS_BIG_ENDIAN
	fval = (float) ntohl((uint32_t)*ptr);
#else
	int i1=*(int *)ptr;
 	i1=affy_swap32(i1);
	fval=(float)affy_swap32((uint32_t)(*(float *)&i1));
#endif
	return fval;
}

uint32_t
MmGetUInt32_N(uint32_t* ptr)
{
  uint8_t* cptr=(uint8_t*)ptr;
  uint32_t val=0;
  val|=(*cptr++)<<24;
  val|=(*cptr++)<<16;
  val|=(*cptr++)<< 8;
  val|=(*cptr++);
  return val;
}
