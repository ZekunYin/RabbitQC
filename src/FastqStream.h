/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/
/*
  This file is modified by SDU HPC lab for RabbitQC project.

  Last modified: JULY2019 
*/
#ifndef H_FASTQSTREAM
#define H_FASTQSTREAM

#include "Globals.h"

#include "Fastq.h"
#include "Buffer.h"
#include "zlib/zlib.h"
#include "util.h"
#if defined (_WIN32)
#   define _CRT_SECURE_NO_WARNINGS
#   pragma warning(disable : 4996) // D_SCL_SECURE
#   pragma warning(disable : 4244) // conversion uint64 to uint32
//# pragma warning(disable : 4267)
#   define FOPEN    fopen
#   define FSEEK    _fseeki64
#   define FTELL    _ftelli64
#   define FCLOSE   fclose
#elif __APPLE__ // Apple by default suport 64 bit file operations (Darwin 10.5+)
#   define FOPEN    fopen
#   define FSEEK    fseek
#   define FTELL    ftell
#   define FCLOSE   fclose
#else
#   if !defined(_LARGEFILE_SOURCE)
#       define _LARGEFILE_SOURCE
#       if !defined(_LARGEFILE64_SOURCE)
#           define _LARGEFILE64_SOURCE
#       endif
#   endif
#   if defined(_FILE_OFFSET_BITS) && (_FILE_OFFSET_BITS != 64)
#       undef _FILE_OFFSET_BITS
#   endif
#   if !defined(_FILE_OFFSET_BITS)
#       define _FILE_OFFSET_BITS 64
#   endif
#   define FOPEN    fopen64
#   define FSEEK    fseeko64
#   define FTELL    ftello64
#   define FCLOSE   fclose
#endif

namespace dsrc
{

namespace fq
{

class FastqFileReader
{
private:
	static const uint32 SwapBufferSize = 1 << 20;
	//static const uint32 SwapBufferSize = 1 << 13;

public:
	FastqFileReader(const std::string& fileName_)
		:	swapBuffer(SwapBufferSize)
		,	bufferSize(0)
		,	eof(false)
		,	usesCrlf(false)
		,	isZipped(false)
	{	
		if(ends_with(fileName_,".gz")){
			mZipFile = gzopen(fileName_.c_str(),"r");
			isZipped=true;
			gzrewind(mZipFile);

		}else{
			mFile = FOPEN(fileName_.c_str(), "rb");
			if(mFile == NULL){
				throw DsrcException(("Can not open file to read: " + fileName_).c_str()); //--------------need to change----------//
			}
		}
		
			
	}

	~FastqFileReader()
	{
		if(mFile != NULL || mZipFile !=NULL)
			Close();
		delete mFile;
		//delete mZipFile;
	}

	bool Eof() const
	{
		return eof;
	}

	bool ReadNextChunk(FastqDataChunk* chunk_);
	bool ReadNextPairedChunk(FastqDataChunk* chunk_);
	void Close()
	{
		if(mFile != NULL)
			FCLOSE(mFile);
		if(mZipFile !=NULL && isZipped){
			gzclose(mZipFile);
			mZipFile=NULL;
		}

		mFile = NULL;
	}

	int64 Read(byte* memory_, uint64 size_)
	{	if(isZipped){
			int64 n = gzread(mZipFile,memory_,size_);
			if(n == -1)
				cerr<<"Error to read gzip file" <<endl;
			return n;
		}
		else{
			int64 n = fread(memory_, 1, size_, mFile) ;
			return n;
		}
		
		
	}
private:
	core::Buffer	swapBuffer;
	uint64			bufferSize;
	bool			eof;
	bool			usesCrlf;
	bool			isZipped;
	FILE*           mFile;
	gzFile          mZipFile;

	uint64 lastOneReadPos;
	uint64 lastTwoReadPos;

	uint64 GetNextRecordPos(uchar* data_, uint64 pos_, const uint64 size_);
	uint64 GetPreviousRecordPos(uchar* data_, uint64 pos_, const uint64 size_);

	void SkipToEol(uchar* data_, uint64& pos_, const uint64 size_)
	{
		//cout << "SkipToEol " << pos_ << " " << size_ << endl;
		ASSERT(pos_ < size_);

		while (data_[pos_] != '\n' && data_[pos_] != '\r' && pos_ < size_)
			++pos_;

		if (data_[pos_] == '\r' && pos_ < size_)
		{
			if (data_[pos_ + 1] == '\n')
			{
				usesCrlf = true;
				++pos_;
			}
		}
	}
	//跳转到行首
	void SkipToSol(uchar* data_, uint64& pos_, const uint64 size_){
		//std::cout<<"SkipToSol:"<<data_[pos_]<<std::endl;
		ASSERT(pos_ < size_);
		if(data_[pos_] =='\n'){
			--pos_;
		}
		if(data_[pos_] =='\r'){
			usesCrlf = true;
			pos_--;
		}
		//找到换行符
		while(data_[pos_] != '\n' &&  data_[pos_] != '\r'){
			//std::cout<<"pos_;"<<pos_<<std::endl;
			--pos_;
		}
		if(data_[pos_] =='\n'){
			--pos_;
		}
		if(data_[pos_] =='\r'){
			usesCrlf = true;
			pos_--;
		}

	}
	


};


} // namespace fq

} // namespace dsrc

#endif // H_FASTQSTREAM
