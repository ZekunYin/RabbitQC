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
  
#include "FastqIo.h"

#include <vector>
#include <map>

#include "Buffer.h"
#include "FastqStream.h"
#include "stdio.h"
#include "read.h"

namespace dsrc
{

namespace fq
{
void FastqReader::readChunk()
{
	FastqDataChunk* part = NULL;

	recordsPool.Acquire(part);
	printf("fastqio: ready to into while\n");

	//while (!errorHandler.IsError() && fileReader.ReadNextChunk(part))
	while(fileReader.ReadNextChunk(part))
	{
		ASSERT(part->size > 0);

		//recordsQueue.Push(numParts, part); //[haoz:] Push:把<numparts, part>组成pair然后放到recordsQueue里面然后notifiy_one
		printf("numParts is %d\n", numParts);
		numParts++;

		recordsPool.Release(part);	
		recordsPool.Acquire(part);
	}

	ASSERT(part->size == 0);
	recordsPool.Release(part);		// the last empty part

	//recordsQueue.SetCompleted();
}

FastqDataChunk* FastqReader::readNextChunk(){
	FastqDataChunk* part = NULL;
	recordsPool.Acquire(part);
	if(fileReader.ReadNextChunk(part))
	{
		return part;
	}
	else
	{
		recordsPool.Release(part);
		return NULL;
	}
}
FastqDataChunk* FastqReader::readNextPairedChunk(){
	FastqDataChunk* part = NULL;
	recordsPool.Acquire(part);
	if(fileReader.ReadNextPairedChunk(part))
	{
		//std::cout <<"read one chunk:   " << std::endl;
		//std::cerr << (char*)part->data.Pointer() << endl;
		//std::cerr << "chunk end" << endl;
		return part;
	}
	else
	{
		//std::cout <<"read one chunk failed :   " << std::endl;
		recordsPool.Release(part);
		return NULL;
	}
}
int chunkFormat(FastqDataChunk* &chunk, std::vector<Read*> &data, bool mHasQuality){
	//format a whole chunk and return number of reads
	int seq_count = 0;
	int line_count = 0;
	int pos_ = 0;
	while(true){
		string name = getLine(chunk, pos_);
		if(name.empty()) break;//dsrc guarantees that read are completed!
		//std::cerr << name << std::endl;

		string sequence = getLine(chunk, pos_);
		//std::cerr<< sequence << std::endl;

		string strand = getLine(chunk, pos_);
		//std::cerr << strand << std::endl;
		if(!mHasQuality){
			string quality = string(sequence.length(), 'K');
			//std::cerr << quality << std::endl;
			data.push_back(new Read(name, sequence, strand, quality));
			seq_count++;

		}else{
			string quality = getLine(chunk, pos_);
			//std::cerr << quality << std::endl;
			data.push_back(new Read(name, sequence, strand, quality));
			seq_count++;
		}
	}

	return seq_count;
}
int pairedChunkFormat(FastqDataChunk* &chunk, std::vector<ReadPair*> &data, bool mHasQuality){
	//format a whole chunk and return number of reads
	int seq_count = 0;
	int pos_ = 0;
	//处理一开始的两条序列
	Read* first = getOnePairedRead(chunk,pos_,mHasQuality);
	Read* second = getOnePairedRead(chunk,pos_,mHasQuality);


	if(first->mName.substr(0, first->mName.find(' ')) == second->mName.substr(0, second->mName.find(' ')) ){
		//cerr << first->mName<<" + " << endl;

		data.push_back(new ReadPair(first,second));
		seq_count++;

	}else{
		//std::cout<<"Right"<<std::endl;
		//cerr << first->mName<<" - "<<endl;
		Read* third = getOnePairedRead(chunk,pos_,mHasQuality);
		data.push_back(new ReadPair(second,third));
		seq_count++;
	}
	while(true){
		Read* t1=getOnePairedRead(chunk,pos_,mHasQuality);
		if(t1 == NULL){
			break;
		}
		Read* t2 =getOnePairedRead(chunk,pos_,mHasQuality);
		if(t2 == NULL){
			//cerr << t1->mName << " with + without -" << endl;
			break;
		}
		data.push_back(new ReadPair(t1,t2));
		seq_count++;
	}

	//for(int i =0;i<seq_count;i++){
	//	cerr << data[i]->mLeft->mName<<endl;
	//	cerr << data[i]->mLeft->mSeq.mStr<<endl;
	//	cerr << data[i]->mLeft->mStrand<<endl;
	//	cerr << data[i]->mLeft->mQuality<<endl;
	//	cerr << data[i]->mRight->mName<<endl;
	//	cerr << data[i]->mRight->mSeq.mStr<<endl;
	//	cerr << data[i]->mRight->mStrand<<endl;
	//	cerr << data[i]->mRight->mQuality<<endl;
	//	
	//}
	//cerr << data[seq_count-1]->mLeft->mName<<endl;
		
	return seq_count;
}
Read* getOnePairedRead(FastqDataChunk* &chunk,int &pos_, bool mHasQuality){
	while(true){
		string name = getLine(chunk, pos_);
		if(name.empty()) return NULL;
		//std::cerr << name << std::endl;
		string sequence = getLine(chunk, pos_);

		string strand = getLine(chunk, pos_);
		if(!mHasQuality){
			string quality = string(sequence.length(), 'K');
			return new Read(name, sequence, strand, quality);
	
		}else{
			string quality = getLine(chunk, pos_);
			//std::cerr << quality << std::endl;
			return new Read(name, sequence, strand, quality);
		}
	}
}


string getLine(FastqDataChunk* &chunk, int &pos){
	int start_pos = pos;
	char* data = (char *)chunk->data.Pointer();

	while(pos <= (chunk->size + 1)){
		if(data[pos] == '\n' || data[pos] == '\r' || pos == (chunk->size + 1)){
			//find a line
			pos++;
			return string(data+start_pos, pos-start_pos - 1);
		}
		else{
			pos++;
		}
	}
	return "";
}

} // namesapce fq

} // namespace dsrc
