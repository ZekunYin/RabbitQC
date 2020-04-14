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
 
#include "FastqStream.h"
#include <iostream>
namespace dsrc
{

namespace fq
{

//bool IFastqStreamReader::ReadNextChunk(FastqDataChunk* chunk_)
bool FastqFileReader::ReadNextChunk(FastqDataChunk* chunk_)
{
	if (Eof())
	{
		chunk_->size = 0;
		return false;
	}

	// flush the data from previous incomplete chunk
	uchar* data = chunk_->data.Pointer();
	const uint64 cbufSize = chunk_->data.Size();
	chunk_->size = 0;
	int64 toRead = cbufSize - bufferSize;
	// buffersize: class IFastqStreamReader de 全局变量 ， 初始化 = 0
	//---------------------------------
	if (bufferSize > 0)
	{
		std::copy(swapBuffer.Pointer(), swapBuffer.Pointer() + bufferSize, data);
		chunk_->size = bufferSize;
		bufferSize = 0;
	}

	// read the next chunk
	int64 r = Read(data + chunk_->size, toRead);
	//std::cout << "r is :" << r << std::endl;
	
	if (r > 0)
	{
		if (r == toRead)	// somewhere before end
		{
		    uint64 chunkEnd = cbufSize - SwapBufferSize; // Swapbuffersize: 1 << 13
			//std::cout << "chunkend  cbufsize Swapbuffersize: " << chunkEnd <<" "<< cbufSize << " " << SwapBufferSize << std::endl;
			chunkEnd = GetNextRecordPos(data, chunkEnd, cbufSize);
			chunk_->size = chunkEnd - 1;
			if (usesCrlf)
				chunk_->size -= 1;

			std::copy(data + chunkEnd, data + cbufSize, swapBuffer.Pointer());
			bufferSize = cbufSize - chunkEnd;
		}
		else				// at the end of file
		{
			chunk_->size += r - 1;	// skip the last EOF symbol
			if (usesCrlf)
				chunk_->size -= 1;

			eof = true;
		}
	}
	else
	{
		eof = true;
	}

	return true;
}
bool FastqFileReader::ReadNextPairedChunk(FastqDataChunk* chunk_)
{	//std::cout<<"ReadNextPairedChunk: enter" <<std::endl;
	if (Eof())
	{
		chunk_->size = 0;
		return false;
	}

	// flush the data from previous incomplete chunk
	uchar* data = chunk_->data.Pointer();
	const uint64 cbufSize = chunk_->data.Size();
	chunk_->size = 0;
	int64 toRead = cbufSize - bufferSize;
	// buffersize: class IFastqStreamReader de 全局变量 ， 初始化 = 0
	//---------------------------------
	if (bufferSize > 0)
	{
		std::copy(swapBuffer.Pointer(), swapBuffer.Pointer() + bufferSize, data);
		chunk_->size = bufferSize;
		bufferSize = 0;
	}

	// read the next chunk
	int64 r = Read(data + chunk_->size, toRead);
	//std::cout << "r is :" << r << std::endl;
	
	if (r > 0)
	{
		if (r == toRead)	// somewhere before end
		{
		    //uint64 chunkEnd = cbufSize - SwapBufferSize; // Swapbuffersize: 1 << 13
			//std::cout << "chunkend  : " << cbufSize-1 <<std::endl;
			lastOneReadPos = GetPreviousRecordPos(data, cbufSize-1,cbufSize);
			//std::cout << "one：" << lastOneReadPos<<std::endl;
			lastTwoReadPos = GetPreviousRecordPos(data, lastOneReadPos-1,cbufSize);
			//std::cout << "two :" << lastTwoReadPos<<std::endl;
			chunk_->size = lastOneReadPos - 1;
			if (usesCrlf) // TODO windows先不管
				chunk_->size -= 1;

			std::copy(data + lastTwoReadPos, data + cbufSize, swapBuffer.Pointer());
			bufferSize = cbufSize - lastTwoReadPos;
		}
		else				// at the end of file
		{
			chunk_->size += r - 1;	// skip the last EOF symbol
			if (usesCrlf)
				chunk_->size -= 1;

			eof = true;
		}
	}
	else
	{
		eof = true;
	}
	//std::cout<<"ReadNextPairedChunk: success@!!!!!!!" <<std::endl;
	return true;
}
//需要向前找到倒数第2条read的开头，不仅仅是@号这么简单
uint64 FastqFileReader::GetPreviousRecordPos(uchar* data_, uint64 pos_,const uint64 size_)
{	int offset =2;

	SkipToSol(data_,pos_,size_);
	if(usesCrlf){
		offset=3;
	}
	while(data_[pos_+offset] !='@'){ //+2
		//std::cout<<"pos_"<<pos_<<std::endl;
		SkipToSol(data_,pos_,size_);
	}
	//标记一下，看是否是质量分
	uint64 pos0=pos_+offset;
	SkipToSol(data_,pos_,size_);
	if(data_[pos_+offset]== '+'){
		//说明上一个@号是质量分
		SkipToSol(data_,pos_,size_);
		SkipToSol(data_,pos_,size_);
		//此时应该是name
		if(data_[pos_+offset]!='@'){
			std::cout << "core dump is " << data_[pos_+offset] << std::endl;
			return pos0; //fix warning
		}else{
			return pos_+offset;
		}
	}
	else{
		return pos0;
	}
}




//uint64 IFastqStreamReader::GetNextRecordPos(uchar* data_, uint64 pos_, const uint64 size_)
uint64 FastqFileReader::GetNextRecordPos(uchar* data_, uint64 pos_, const uint64 size_)
{
	SkipToEol(data_, pos_, size_);
	++pos_;

	// find beginning of the next record
	while (data_[pos_] != '@')
	{
		SkipToEol(data_, pos_, size_);
		++pos_;
	}
	uint64 pos0 = pos_;
		
	SkipToEol(data_, pos_, size_);
	++pos_;

	if (data_[pos_] == '@')			// previous one was a quality field
		return pos_;
	
	SkipToEol(data_, pos_, size_);
	++pos_;
	if(data_[pos_] != '+')
		std::cout << "core dump is pos: " << pos_ << " char: " << data_[pos_] << std::endl;	
	ASSERT(data_[pos_] == '+');	// pos0 was the start of tag
	return pos0;
}

} // namespace fq

} // namespace dsrc

