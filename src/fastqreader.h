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
#ifndef FASTQ_READER_H
#define FASTQ_READER_H

#include <stdio.h>
#include <stdlib.h>
#include "read.h"
#ifdef DYNAMIC_ZLIB
  #include <zlib.h>
#else
  #include "zlib/zlib.h"
#endif
#include "common.h"
#include <iostream>
#include <fstream>

#include "FastqIo.h"
#include "FastqStream.h"

class FastqReader{
public:
	FastqReader(string filename, bool hasQuality = true, bool phred64=false);
	~FastqReader();
	bool isZipped();

	void getBytes(size_t& bytesRead, size_t& bytesTotal);

	//this function is not thread-safe
	//do not call read() of a same FastqReader object from different threads concurrently
	Read* read();
	bool eof();
	bool hasNoLineBreakAtEnd();
public: static bool isZipFastq(string filename);
	static bool isFastq(string filename);
	static bool test();

private:
	void init();
	void close();
	string getLine();
	void clearLineBreaks(char* line);
	void readToBuf();

private:
	string mFilename;
	gzFile mZipFile;
	FILE* mFile;
	bool mZipped;
	bool mHasQuality;
	bool mPhred64;
	char* mBuf;
	int mBufDataLen;
	int mBufUsedLen;
	bool mStdinMode;
	bool mHasNoLineBreakAtEnd;

};

class FastqReaderPair{
public:
	FastqReaderPair(FastqReader* left, FastqReader* right);
	FastqReaderPair(string leftName, string rightName, bool hasQuality = true, bool phred64 = false, bool interleaved = false);
	~FastqReaderPair();
	ReadPair* read();
public:
	FastqReader* mLeft;
	FastqReader* mRight;
	bool mInterleaved;
};

class FastqChunkReaderPair{
public:
    FastqChunkReaderPair(dsrc::fq::FastqReader* left, dsrc::fq::FastqReader* right);
	FastqChunkReaderPair(string leftName, string rightName, bool hasQuality = true, bool phred64 = false, bool interleaved = false);
	~FastqChunkReaderPair();
	void SkipToEol(dsrc::uchar* data_, uint64& pos_, const uint64 size_);
	uint64 GetNextRecordPos(dsrc::uchar* data_, uint64 pos_, const uint64 size_);
	ChunkPair* readNextChunkPair();
	ChunkPair* readNextChunkPair_interleaved();
public:
	dsrc::fq::FastqDataPool* fastqPool_left;
	dsrc::fq::FastqFileReader* fileReader_left;
	dsrc::fq::FastqDataPool* fastqPool_right;
	dsrc::fq::FastqFileReader* fileReader_right;
	dsrc::fq::FastqReader* mLeft;
	dsrc::fq::FastqReader* mRight;
	bool mInterleaved;

	//---mamber variable for readNextchunkpair
	dsrc::core::Buffer swapBuffer_left;
	dsrc::core::Buffer swapBuffer_right;
	bool usesCrlf;
	bool eof;
	uint64 bufferSize_left;
	uint64 bufferSize_right;
	
};

#endif
