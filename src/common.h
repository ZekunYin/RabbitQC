/*
MIT License

Copyright (c) 2017 OpenGene - Open Source Genetics Toolbox

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

/*
	Last modified by SDU HPC lab: JULY2019.

*/
#ifndef COMMON_H
#define COMMON_H

#define RABBITQC_VER "0.0.1"

#define _DEBUG false

typedef long int64;
typedef unsigned long uint64;

typedef int int32;
typedef unsigned int uint32;

typedef short int16;
typedef unsigned short uint16;

typedef char int8;
typedef unsigned char uint8;

#pragma pack(2) 


#pragma pack() 

// the limit of the queue to store the packs
// error may happen if it generates more packs than this number
static const int PACK_NUM_LIMIT  = 10000000;

// how many reads one pack has
static const int PACK_SIZE = 1000;

// if one pack is produced, but not consumed, it will be kept in the memory
// this number limit the number of in memory packs
// if the number of in memory packs is full, the producer thread should sleep
static const int PACK_IN_MEM_LIMIT = 500;

// if read number is more than this, warn it
static const int WARN_STANDALONE_READ_LIMIT = 10000;

// different filtering results, bigger number means worse
// if r1 and r2 are both failed, then the bigger one of the two results will be recorded
// we reserve some gaps for future types to be added
static const int PASS_FILTER = 0;
static const int FAIL_POLY_X = 4;
static const int FAIL_OVERLAP = 8;
static const int FAIL_N_BASE = 12;
static const int FAIL_LENGTH = 16;
static const int FAIL_TOO_LONG = 17;
static const int FAIL_QUALITY = 20;
static const int FAIL_COMPLEXITY = 24;

// how many types in total we support
static const int FILTER_RESULT_TYPES = 32;




/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/


//#include "../include/dsrc/Globals.h"
#include "Globals.h"

#ifndef NDEBUG
#	define DEBUG 1
#endif

#define BIT(x)							(1 << (x))
#define BIT_ISSET(x, pos)				((x & BIT(pos)) != 0)
#define BIT_SET(x, pos)					(x |= BIT(pos))
#define BIT_UNSET(x, pos)				(x &= ~(BIT(pos)))
#define MIN(x,y)						((x) <= (y) ? (x) : (y))
#define MAX(x,y)						((x) >= (y) ? (x) : (y))
#define ABS(x)							((x) >=  0  ? (x) : -(x))
#define SIGN(x)							((x) >=  0  ?  1  : -1)
#define REC_EXTENSION_FACTOR(size)		( ((size) / 4 > 1024) ? ((size) / 4) : 1024 )
#define MEM_EXTENSION_FACTOR(size)		REC_EXTENSION_FACTOR(size)

#if defined (_WIN32)
#	define _CRT_SECURE_NO_WARNINGS
#	pragma warning(disable : 4996) // D_SCL_SECURE
#	pragma warning(disable : 4244) // conversion uint64 to uint32
#	pragma warning(disable : 4267)
#	pragma warning(disable : 4800) // conversion byte to bool
#endif

// TODO: refactor raw data structs to avoid using <string> as a member
#include <string>



#define COMPILE_TIME_ASSERT(COND, MSG) typedef char static_assertion_##MSG[(!!(COND))*2-1]
#define COMPILE_TIME_ASSERT1(X, L) COMPILE_TIME_ASSERT(X,static_assertion_at_line_##L)
#define COMPILE_TIME_ASSERT2(X, L) COMPILE_TIME_ASSERT1(X,L)
#define STATIC_ASSERT(X)    COMPILE_TIME_ASSERT2(X,__LINE__)


namespace dsrc
{

namespace fq
{

// ********************************************************************************************
struct FastqDatasetType
{
	static const uint32 AutoQualityOffset = 0;
	static const uint32 DefaultQualityOffset = 33;

	uint32	qualityOffset;
	bool	plusRepetition;
	bool	colorSpace;


	FastqDatasetType()
		:	qualityOffset(AutoQualityOffset)
		,	plusRepetition(false)
		,	colorSpace(false)
	{}

	static FastqDatasetType Default()
	{
		FastqDatasetType ds;
		ds.qualityOffset = AutoQualityOffset;
		ds.plusRepetition = false;
		ds.colorSpace = false;
		return ds;
	}
};

struct StreamsInfo
{
	enum StreamName
	{
		MetaStream = 0,
		TagStream,
		DnaStream,
		QualityStream,

		StreamCount = 4
	};

	uint64 sizes[4];

	StreamsInfo()
	{
		Clear();
	}

	void Clear()
	{
		std::fill(sizes, sizes + StreamCount, 0);
	}
};

struct FastqRecord;

} // namespace fq


namespace comp
{

struct CompressionSettings
{
	static const uint32 MaxDnaOrder = 9;
	static const uint32 MaxQualityOrder = 6;
	static const uint32 DefaultDnaOrder = 0;
	static const uint32 DefaultQualityOrder = 0;
	static const uint32 DefaultTagPreserveFlags = 0;		// 0 -- keep all

	uint32	dnaOrder;
	uint32	qualityOrder;
	uint64	tagPreserveFlags;
	bool	lossy;
	bool	calculateCrc32;

	CompressionSettings()
		:	dnaOrder(0)
		,	qualityOrder(0)
		,	tagPreserveFlags(DefaultTagPreserveFlags)
		,	lossy(false)
		,	calculateCrc32(false)
	{}

	static CompressionSettings Default()
	{
		CompressionSettings s;
		s.dnaOrder = DefaultDnaOrder;
		s.qualityOrder = DefaultQualityOrder;
		s.tagPreserveFlags = DefaultTagPreserveFlags;
		s.lossy = false;
		s.calculateCrc32 = false;
		return s;
	}
};

struct InputParameters
{
	static const uint32 DefaultQualityOffset = fq::FastqDatasetType::AutoQualityOffset;
	static const uint32 DefaultDnaCompressionLevel = 0;
	static const uint32 DefaultQualityCompressionLevel = 0;
	static const uint32 DefaultProcessingThreadNum = 2;
	static const uint64 DefaultTagPreserveFlags = 0;
	static const uint32 DefaultFastqBufferSizeMB = 8;

	static const bool DefaultLossyCompressionMode = false;
	static const bool DefaultCalculateCrc32 = false;


	uint32 qualityOffset;
	uint32 dnaCompressionLevel;
	uint32 qualityCompressionLevel;
	uint32 threadNum;
	uint64 tagPreserveFlags;

	uint32 fastqBufferSizeMB;
	bool lossyCompression;
	bool calculateCrc32;
	bool useFastqStdIo;

	std::string inputFilename;
	std::string outputFilename;

	InputParameters()
		:	qualityOffset(DefaultQualityOffset)
		,	dnaCompressionLevel(DefaultDnaCompressionLevel)
		,	qualityCompressionLevel(DefaultQualityCompressionLevel)
		,	threadNum(DefaultProcessingThreadNum)
		,	tagPreserveFlags(DefaultTagPreserveFlags)
		,	fastqBufferSizeMB(DefaultFastqBufferSizeMB)
		,	lossyCompression(DefaultLossyCompressionMode)
		,	calculateCrc32(DefaultCalculateCrc32)
		,	useFastqStdIo(false)
	{}

	static InputParameters Default()
	{
		InputParameters args;
		return args;
	}
};

struct Field;

struct DnaStats;
struct QualityStats;

class BlockCompressor;
class HuffmanEncoder;

struct DsrcDataChunk;

} // namespace comp


namespace core
{

class Buffer;
class BitMemoryReader;
class BitMemoryWriter;
class ErrorHandler;

} // namespace core

} // namespace dsrc


#endif // _COMMON_H
 
