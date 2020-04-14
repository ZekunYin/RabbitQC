#include "TGSProcessor.h"
#include "fastqreader.h"
#include <iostream>
#include <unistd.h>
#include <functional>
#include <thread>
#include <memory.h>
#include "util.h"
#include "jsonreporter.h"
#include "htmlreporter.h"
#include "TGSStats.h"


TGSProcessor::TGSProcessor(Options* opt){
    mOptions = opt;
    mProduceFinished = false;
    mFinishedThreads = 0;
    //mFilter = new Filter(opt);
    //mOutStream = NULL;
    //mZipFile = NULL;
    //mUmiProcessor = new UmiProcessor(opt);
    mLeftWriter =  NULL;

    //mDuplicate = NULL;
    //if(mOptions->duplicate.enabled) {
    //    mDuplicate = new Duplicate(mOptions);
    //}

	//dsrc faster reader 
    fastqPool = new dsrc::fq::FastqDataPool(128,1<<22);
}
TGSProcessor::~TGSProcessor(){
	delete fastqPool;
}

/*
void SingleEndProcessor::initOutput() {
    if(mOptions->out1.empty())
        return;
    mLeftWriter = new WriterThread(mOptions, mOptions->out1);
}

void SingleEndProcessor::closeOutput() {
    if(mLeftWriter) {
        delete mLeftWriter;
        mLeftWriter = NULL;
    }
}

void SingleEndProcessor::initConfig(ThreadConfig* config) {
    if(mOptions->out1.empty())
        return;

    if(mOptions->split.enabled) {
        config->initWriterForSplit();
    }
}
*/

bool TGSProcessor::process(){
    //if(!mOptions->split.enabled)
    //    initOutput();

    initPackRepository();
    std::thread producer(std::bind(&TGSProcessor::producerTask, this));

    //TODO: get the correct cycles
    int cycle = 151;
    ThreadConfig** configs = new ThreadConfig*[mOptions->thread];
    for(int t=0; t<mOptions->thread; t++){
		configs[t] = new ThreadConfig(mOptions, t, false);
		//initConfig(configs[t]);
    }

    std::thread** threads = new thread*[mOptions->thread];
    for(int t=0; t<mOptions->thread; t++){
        threads[t] = new std::thread(std::bind(&TGSProcessor::consumerTask, this, configs[t]));
    }

    //std::thread* leftWriterThread = NULL;
    //if(mLeftWriter)
    //    leftWriterThread = new std::thread(std::bind(&SingleEndProcessor::writeTask, this, mLeftWriter));

    producer.join();
    for(int t=0; t<mOptions->thread; t++){
        threads[t]->join();
    }

    //if(!mOptions->split.enabled) {
    //    if(leftWriterThread)
    //        leftWriterThread->join();
    //}

    if(mOptions->verbose)
        loginfo("start to generate reports\n");

    // merge stats and read filter results
    vector<TGSStats*> dataStats;
    //vector<FilterResult*> filterResults;
    for(int t=0; t<mOptions->thread; t++){
        dataStats.push_back(configs[t]->getTGSStats());
    }
    TGSStats* finalStats = TGSStats::merge(dataStats);
	//TGSStats::merge(dataStats, finalStats);
	
    //cerr << "Stats Summary:"<<endl;
    //finalStats->print();
    //cerr << endl;


    //int* dupHist = NULL;
    //double* dupMeanTlen = NULL;
    //double* dupMeanGC = NULL;
    //double dupRate = 0.0;
    //if(mOptions->duplicate.enabled) {
    //    dupHist = new int[mOptions->duplicate.histSize];
    //    memset(dupHist, 0, sizeof(int) * mOptions->duplicate.histSize);
    //    dupMeanGC = new double[mOptions->duplicate.histSize];
    //    memset(dupMeanGC, 0, sizeof(double) * mOptions->duplicate.histSize);
    //    dupRate = mDuplicate->statAll(dupHist, dupMeanGC, mOptions->duplicate.histSize);
    //    cerr << endl;
    //    cerr << "Duplication rate (may be overestimated since this is SE data): " << dupRate * 100.0 << "%" << endl;
    //}

    //// make JSON report
    //JsonReporter jr(mOptions);
    //jr.setDupHist(dupHist, dupMeanGC, dupRate);
    //jr.report(finalFilterResult, finalPreStats, finalPostStats);

    // make HTML report
    HtmlReporter hr(mOptions);
    //hr.setDupHist(dupHist, dupMeanGC, dupRate);
    hr.report3(finalStats);

    // clean up
    for(int t=0; t<mOptions->thread; t++){
        delete threads[t];
        threads[t] = NULL;
        delete configs[t];
        configs[t] = NULL;
    }

    delete finalStats;

    //if(mOptions->duplicate.enabled) {
    //    delete[] dupHist;
    //    delete[] dupMeanGC;
    //}

    delete[] threads;
    delete[] configs;

    //if(leftWriterThread)
    //    delete leftWriterThread;

    //if(!mOptions->split.enabled)
    //    closeOutput();

    return true;
}
bool TGSProcessor::processSingleEnd(ReadPack* pack, ThreadConfig* config){

	//debug for losing last line
		//cerr << pack->data[pack->count - 1]->mName << endl;
		//cerr << pack->data[pack->count - 1]->mQuality << endl;
	//debug
    //string outstr;
    //int readPassed = 0;
	for(int p = 0; p < pack->count; p++){
	  Read* or1 = pack->data[p];
	  config -> getTGSStats() -> tgsStatRead(or1);

	  if(or1 != NULL){
		  delete or1;
	  }
	}
	//1. add TGStats class in project
	//2. add getTGStats function in ThreadConfig file;
	//3. add TGStats menber variable in ThreadConfig class
	//---------------------------------------------------
	
	delete pack;

	return true; //fix warning

}

void TGSProcessor::initPackRepository() {
    //mRepo.packBuffer = new ReadPack*[PACK_NUM_LIMIT];
    mRepo.packBuffer = new dsrc::fq::FastqDataChunk*[PACK_NUM_LIMIT];
    //memset(mRepo.packBuffer, 0, sizeof(ReadPack*)*PACK_NUM_LIMIT);
    memset(mRepo.packBuffer, 0, sizeof(dsrc::fq::FastqDataChunk*)*PACK_NUM_LIMIT);
    mRepo.writePos = 0;
    mRepo.readPos = 0;
    //mRepo.readCounter = 0;
    
}

void TGSProcessor::destroyPackRepository() {
    delete mRepo.packBuffer;
    mRepo.packBuffer = NULL;
}

void TGSProcessor::producePack(dsrc::fq::FastqDataChunk* pack){
    /*while(((mRepo.writePos + 1) % PACK_NUM_LIMIT)
        == mRepo.readPos) {
        //mRepo.repoNotFull.wait(lock);
    }*/

    mRepo.packBuffer[mRepo.writePos] = pack;
	//cerr << "chunk in producePack" << endl;
	//cerr << (char *)pack->data.Pointer() << endl;
    mRepo.writePos++;

}

void TGSProcessor::consumePack(ThreadConfig* config){
	dsrc::fq::FastqDataChunk* chunk;
    ReadPack* data = new ReadPack;
    //std::unique_lock<std::mutex> lock(mRepo.mtx);
    // buffer is empty, just wait here.
    /*while(mRepo.writePos % PACK_NUM_LIMIT == mRepo.readPos % PACK_NUM_LIMIT) {
        if(mProduceFinished){
            //lock.unlock();
            return;
        }
        //mRepo.repoNotEmpty.wait(lock);
    }*/

    mInputMtx.lock();
    while(mRepo.writePos <= mRepo.readPos) {
        usleep(1000);
        if(mProduceFinished) {
            mInputMtx.unlock();
            return;
        }
    }
    //data = mRepo.packBuffer[mRepo.readPos];
    chunk = mRepo.packBuffer[mRepo.readPos];
	//cerr << "read pos is " << mRepo.readPos << endl;
	//cerr << (char*)chunk->data.Pointer() << endl;

    mRepo.readPos++;

    /*if (mRepo.readPos >= PACK_NUM_LIMIT)
        mRepo.readPos = 0;*/
    mInputMtx.unlock();

	//data format for from dsrc to fastp
	data->count = dsrc::fq::chunkFormat(chunk, data->data, true);
	//cerr << (char*)chunk->data.Pointer() << endl;
	fastqPool->Release(chunk);	

    //lock.unlock();
    //mRepo.repoNotFull.notify_all();

    processSingleEnd(data, config);

}

void TGSProcessor::producerTask()
{
    if(mOptions->verbose)
        loginfo("start to load data");
    long lastReported = 0;
    int slept = 0;
    long readNum = 0;
    bool splitSizeReEvaluated = false;
    //Read** data = new Read*[PACK_SIZE];
    //memset(data, 0, sizeof(Read*)*PACK_SIZE);
    //FastqReader reader(mOptions->in1, true, mOptions->phred64);
	//dsrc::fq::FastqDataPool* fastqPool = new dsrc::fq::FastqDataPool(32,1<<22);
	dsrc::fq::FastqFileReader* fileReader = new dsrc::fq::FastqFileReader(mOptions->in1);
	dsrc::fq::FastqReader* dataReader = new dsrc::fq::FastqReader(*fileReader, *fastqPool);


	dsrc::fq::FastqDataChunk* chunk;
	while((chunk = dataReader->readNextChunk())!= NULL){
		//cerr << "chunk content in producer =======================" << endl;
		//cerr << (char*) chunk->data.Pointer() << endl;
        producePack(chunk);
		//cerr << (char*)mRepo.packBuffer[0]->data.Pointer() << endl;
		//usleep(200);
        while(mRepo.writePos - mRepo.readPos > PACK_IN_MEM_LIMIT){
            //cerr<<"sleep"<<endl;
            slept++;
            usleep(100);
        }

	}



    //std::unique_lock<std::mutex> lock(mRepo.readCounterMtx);
    mProduceFinished = true;
    if(mOptions->verbose)
        loginfo("all reads loaded, start to monitor thread status");
    //lock.unlock();
	delete fileReader;
	delete dataReader;
}

void TGSProcessor::consumerTask(ThreadConfig* config)
{
    while(true) {
        if(config->canBeStopped()){
            mFinishedThreads++;
            break;
        }
        while(mRepo.writePos <= mRepo.readPos) {
            if(mProduceFinished)
                break;
            usleep(1000);
        }
        //std::unique_lock<std::mutex> lock(mRepo.readCounterMtx);
        if(mProduceFinished && mRepo.writePos == mRepo.readPos){
            mFinishedThreads++;
            if(mOptions->verbose) {
                string msg = "thread " + to_string(config->getThreadId() + 1) + " data processing completed";
                loginfo(msg);
            }
            //lock.unlock();
            break;
        }
        if(mProduceFinished){
            if(mOptions->verbose) {
                string msg = "thread " + to_string(config->getThreadId() + 1) + " is processing the " + to_string(mRepo.readPos) + " / " + to_string(mRepo.writePos) + " pack";
                loginfo(msg);
            }
            consumePack(config);
            //lock.unlock();
        } else {
            //lock.unlock();
            consumePack(config);
        } 
    }

    if(mFinishedThreads == mOptions->thread) {
        if(mLeftWriter)
            mLeftWriter->setInputCompleted();
    }

    if(mOptions->verbose) {
        string msg = "thread " + to_string(config->getThreadId() + 1) + " finished";
        loginfo(msg);
    }
}

//void TGSProcessor::writeTask(WriterThread* config)
//{
//    while(true) {
//        if(config->isCompleted()){
//            // last check for possible threading related issue
//            config->output();
//            break;
//        }
//        config->output();
//    }
//
//    if(mOptions->verbose) {
//        string msg = config->getFilename() + " writer finished";
//        loginfo(msg);
//    }
//}
