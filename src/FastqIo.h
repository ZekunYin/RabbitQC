/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/

#ifndef H_FASTQREADER
#define H_FASTQREADER

#include "Globals.h"

#include "common.h"
#include "Fastq.h"
#include "DataQueue.h"
#include "DataPool.h"
#include "FastqStream.h"
#include "read.h"
#include <vector>

namespace dsrc
{

namespace fq
{

typedef core::TDataQueue<FastqDataChunk> FastqDataQueue;
typedef core::TDataPool<FastqDataChunk> FastqDataPool;


class FastqReader //: public IFastqIoOperator
{
public:
    FastqReader(FastqFileReader& reader_, FastqDataPool& pool_)
	    :   recordsPool(pool_)
		,	fileReader(reader_)
		,	numParts(0)
	{};

	void readChunk();
	FastqDataChunk* readNextChunk();
	//single pe file
	FastqDataChunk* readNextPairedChunk();
	
	int64 Read(byte* memory_, uint64 size_)
	{
		int64 n = fileReader.Read(memory_, size_);
		return n;
	}

private:

	FastqDataPool&      recordsPool;
	FastqFileReader&	fileReader;
	uint32 numParts;
};

int chunkFormat(FastqDataChunk* &chunk, std::vector<Read*>&,bool);

//single pe file 
int pairedChunkFormat(FastqDataChunk* &chunk, std::vector<ReadPair*>&,bool mHasQuality);
Read* getOnePairedRead(FastqDataChunk* &chunk,int &pos_, bool mHasQuality);
//end single pe file
string getLine(FastqDataChunk* &chunk, int &pos);

} // namespace fq

} // namespace dsrc

#endif
