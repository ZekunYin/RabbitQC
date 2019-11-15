#include "peprocessor.h"
#include "seprocessor.h"
#include "fastqreader.h"
#include <iostream>
#include <unistd.h>
#include <functional>
#include <thread>
#include <memory.h>
#include "util.h"
#include "adaptertrimmer.h"
#include "basecorrector.h"
#include "jsonreporter.h"
#include "htmlreporter.h"
#include "polyx.h"

PairEndProcessor::PairEndProcessor(Options* opt){
    mOptions = opt;
    mProduceFinished = false;
    mFinishedThreads = 0;
    mFilter = new Filter(opt);
    mOutStream1 = NULL;
    mZipFile1 = NULL;
    mOutStream2 = NULL;
    mZipFile2 = NULL;
    mUmiProcessor = new UmiProcessor(opt);

    int isizeBufLen = mOptions->insertSizeMax + 1;
    mInsertSizeHist = new long[isizeBufLen];
    memset(mInsertSizeHist, 0, sizeof(long)*isizeBufLen);
    mLeftWriter =  NULL;
    mRightWriter = NULL;

    mDuplicate = NULL;
    if(mOptions->duplicate.enabled) {
        mDuplicate = new Duplicate(mOptions);
    }
	if(mOptions->interleavedInput)
    	fastqPool = new dsrc::fq::FastqDataPool(128,1<<22);
}

PairEndProcessor::~PairEndProcessor() {
    delete mInsertSizeHist;
    if(mDuplicate) {
        delete mDuplicate;
        mDuplicate = NULL;
    }
	if(mOptions->interleavedInput)
    	delete fastqPool;
}

void PairEndProcessor::initOutput() {
    if(mOptions->out1.empty() || mOptions->out2.empty())
        return;
    
    mLeftWriter = new WriterThread(mOptions, mOptions->out1);
    mRightWriter = new WriterThread(mOptions, mOptions->out2);
}

void PairEndProcessor::closeOutput() {
    if(mLeftWriter) {
        delete mLeftWriter;
        mLeftWriter = NULL;
    }
    if(mRightWriter) {
        delete mRightWriter;
        mRightWriter = NULL;
    }
}

void PairEndProcessor::initConfig(ThreadConfig* config) {
    if(mOptions->out1.empty())
        return;
    if(mOptions->split.enabled) {
        config->initWriterForSplit();
    }
}


bool PairEndProcessor::process(){
    if(!mOptions->split.enabled)
        initOutput();

    initPackRepository();
    std::thread producer(std::bind(&PairEndProcessor::producerTask, this));

    //TODO: get the correct cycles
    int cycle = 151;
    ThreadConfig** configs = new ThreadConfig*[mOptions->thread];
    for(int t=0; t<mOptions->thread; t++){
        configs[t] = new ThreadConfig(mOptions, t, true);
        initConfig(configs[t]);
    }

    std::thread** threads = new thread*[mOptions->thread];
    for(int t=0; t<mOptions->thread; t++){
        threads[t] = new std::thread(std::bind(&PairEndProcessor::consumerTask, this, configs[t]));
    }

    std::thread* leftWriterThread = NULL;
    std::thread* rightWriterThread = NULL;
    if(mLeftWriter)
        leftWriterThread = new std::thread(std::bind(&PairEndProcessor::writeTask, this, mLeftWriter));
    if(mRightWriter)
        rightWriterThread = new std::thread(std::bind(&PairEndProcessor::writeTask, this, mRightWriter));

    producer.join();
    for(int t=0; t<mOptions->thread; t++){
        threads[t]->join();
    }

    if(!mOptions->split.enabled) {
        if(leftWriterThread)
            leftWriterThread->join();
        if(rightWriterThread)
            rightWriterThread->join();
    }

    if(mOptions->verbose)
        loginfo("start to generate reports\n");

    // merge stats and filter results
    vector<Stats*> preStats1;
    vector<Stats*> postStats1;
    vector<Stats*> preStats2;
    vector<Stats*> postStats2;
    vector<FilterResult*> filterResults;
    for(int t=0; t<mOptions->thread; t++){
        preStats1.push_back(configs[t]->getPreStats1());
        postStats1.push_back(configs[t]->getPostStats1());
        preStats2.push_back(configs[t]->getPreStats2());
        postStats2.push_back(configs[t]->getPostStats2());
        filterResults.push_back(configs[t]->getFilterResult());
    }
    Stats* finalPreStats1 = Stats::merge(preStats1);
    Stats* finalPostStats1 = Stats::merge(postStats1);
    Stats* finalPreStats2 = Stats::merge(preStats2);
    Stats* finalPostStats2 = Stats::merge(postStats2);
    FilterResult* finalFilterResult = FilterResult::merge(filterResults);

    cerr << "Read1 before filtering:"<<endl;
    finalPreStats1->print();
    cerr << endl;
    cerr << "Read1 after filtering:"<<endl;
    finalPostStats1->print();
    cerr << endl;
    cerr << "Read2 before filtering:"<<endl;
    finalPreStats2->print();
    cerr << endl;
    cerr << "Read2 aftering filtering:"<<endl;
    finalPostStats2->print();

    cerr << endl;
    cerr << "Filtering result:"<<endl;
    finalFilterResult->print();

    int* dupHist = NULL;
    double* dupMeanTlen = NULL;
    double* dupMeanGC = NULL;
    double dupRate = 0.0;
    if(mOptions->duplicate.enabled) {
        dupHist = new int[mOptions->duplicate.histSize];
        memset(dupHist, 0, sizeof(int) * mOptions->duplicate.histSize);
        dupMeanGC = new double[mOptions->duplicate.histSize];
        memset(dupMeanGC, 0, sizeof(double) * mOptions->duplicate.histSize);
        dupRate = mDuplicate->statAll(dupHist, dupMeanGC, mOptions->duplicate.histSize);
        cerr << endl;
        cerr << "Duplication rate: " << dupRate * 100.0 << "%" << endl;
    }

    // insert size distribution
    int peakInsertSize = getPeakInsertSize();
    cerr << endl;
    cerr << "Insert size peak (evaluated by paired-end reads): " << peakInsertSize << endl;

    // make JSON report
    JsonReporter jr(mOptions);
    jr.setDupHist(dupHist, dupMeanGC, dupRate);
    jr.setInsertHist(mInsertSizeHist, peakInsertSize);
    jr.report(finalFilterResult, finalPreStats1, finalPostStats1, finalPreStats2, finalPostStats2);

    // make HTML report
    HtmlReporter hr(mOptions);
    hr.setDupHist(dupHist, dupMeanGC, dupRate);
    hr.setInsertHist(mInsertSizeHist, peakInsertSize);
    hr.report(finalFilterResult, finalPreStats1, finalPostStats1, finalPreStats2, finalPostStats2);

    // clean up
    for(int t=0; t<mOptions->thread; t++){
        delete threads[t];
        threads[t] = NULL;
        delete configs[t];
        configs[t] = NULL;
    }

    delete finalPreStats1;
    delete finalPostStats1;
    delete finalPreStats2;
    delete finalPostStats2;
    delete finalFilterResult;

    if(mOptions->duplicate.enabled) {
        delete[] dupHist;
        delete[] dupMeanGC;
    }

    delete[] threads;
    delete[] configs;

    if(leftWriterThread)
        delete leftWriterThread;
    if(rightWriterThread)
        delete rightWriterThread;

    if(!mOptions->split.enabled)
        closeOutput();

    return true;
}

int PairEndProcessor::getPeakInsertSize() {
    int peak = 0;
    long maxCount = -1;
    for(int i=0; i<mOptions->insertSizeMax; i++) {
        if(mInsertSizeHist[i] > maxCount) {
            peak = i;
            maxCount = mInsertSizeHist[i];
        }
    }
    return peak;
}

bool PairEndProcessor::processPairEnd(ReadPairPack* pack, ThreadConfig* config){
    string outstr1;
    string outstr2;
    string interleaved;
    int readPassed = 0;
	//cerr << "in processPairend, pack->size: " << pack->data.size() << endl;
    for(int p=0;p<pack->count;p++){
        ReadPair* pair = pack->data[p];
        Read* or1 = pair->mLeft;
        Read* or2 = pair->mRight;

        int lowQualNum1 = 0;
        int nBaseNum1 = 0;
        int lowQualNum2 = 0;
        int nBaseNum2 = 0;

        // stats the original read before trimming
        config->getPreStats1()->statRead(or1);
        config->getPreStats2()->statRead(or2);

        // handling the duplication profiling
        if(mDuplicate)
            mDuplicate->statPair(or1, or2);

        // filter by index
        if(mOptions->indexFilter.enabled && mFilter->filterByIndex(or1, or2)) {
            delete pair;
            continue;
        }

        // umi processing
        if(mOptions->umi.enabled)
            mUmiProcessor->process(or1, or2);

        // trim in head and tail, and apply quality cut in sliding window
        Read* r1 = mFilter->trimAndCut(or1, mOptions->trim.front1, mOptions->trim.tail1);
        Read* r2 = mFilter->trimAndCut(or2, mOptions->trim.front2, mOptions->trim.tail2);

        if(r1 != NULL && r2!=NULL) {
            if(mOptions->polyGTrim.enabled)
                PolyX::trimPolyG(r1, r2, config->getFilterResult(), mOptions->polyGTrim.minLen);
            if(mOptions->polyXTrim.enabled)
                PolyX::trimPolyX(r1, r2, config->getFilterResult(), mOptions->polyXTrim.minLen);
        }
        bool isizeEvaluated = false;
        if(r1 != NULL && r2!=NULL && (mOptions->adapter.enabled || mOptions->correction.enabled)){
            OverlapResult ov = OverlapAnalysis::analyze(r1, r2, mOptions->overlapDiffLimit, mOptions->overlapRequire);
            // we only use thread 0 to evaluae ISIZE
            if(config->getThreadId() == 0) {
                statInsertSize(r1, r2, ov);
                isizeEvaluated = true;
            }
            if(mOptions->correction.enabled) {
                BaseCorrector::correctByOverlapAnalysis(r1, r2, config->getFilterResult(), ov);
            }
            if(mOptions->adapter.enabled) {
                bool trimmed = AdapterTrimmer::trimByOverlapAnalysis(r1, r2, config->getFilterResult(), ov);
                if(!trimmed){
                    if(mOptions->adapter.hasSeqR1)
                        AdapterTrimmer::trimBySequence(r1, config->getFilterResult(), mOptions->adapter.sequence, false);
                    if(mOptions->adapter.hasSeqR2)
                        AdapterTrimmer::trimBySequence(r2, config->getFilterResult(), mOptions->adapter.sequenceR2, true);
                }
            }
        }

        if(config->getThreadId() == 0 && !isizeEvaluated && r1 != NULL && r2!=NULL) {
            OverlapResult ov = OverlapAnalysis::analyze(r1, r2, mOptions->overlapDiffLimit, mOptions->overlapRequire);
            statInsertSize(r1, r2, ov);
            isizeEvaluated = true;
        }

        if(r1 != NULL && r2!=NULL) {
            if( mOptions->trim.maxLen1 > 0 && mOptions->trim.maxLen1 < r1->length())
                r1->resize(mOptions->trim.maxLen1);
            if( mOptions->trim.maxLen2 > 0 && mOptions->trim.maxLen2 < r2->length())
                r2->resize(mOptions->trim.maxLen2);
        }

        int result1 = mFilter->passFilter(r1);
        int result2 = mFilter->passFilter(r2);

        config->addFilterResult(max(result1, result2));

        if( r1 != NULL &&  result1 == PASS_FILTER && r2 != NULL && result2 == PASS_FILTER ) {
            
            if(mOptions->outputToSTDOUT) {
                interleaved += r1->toString() + r2->toString();
            } else {
                outstr1 += r1->toString();
                outstr2 += r2->toString();
            }

            // stats the read after filtering
            config->getPostStats1()->statRead(r1);
            config->getPostStats2()->statRead(r2);

            readPassed++;
        }

        delete pair;
        // if no trimming applied, r1 should be identical to or1
        if(r1 != or1 && r1 != NULL)
            delete r1;
        // if no trimming applied, r1 should be identical to or1
        if(r2 != or2 && r2 != NULL)
            delete r2;
    }
    // if splitting output, then no lock is need since different threads write different files
    if(!mOptions->split.enabled)
        mOutputMtx.lock();
    if(mOptions->outputToSTDOUT) {
        // STDOUT output
        fwrite(interleaved.c_str(), 1, interleaved.length(), stdout);
    } else if(mOptions->split.enabled) {
        // split output by each worker thread
        if(!mOptions->out1.empty())
            config->getWriter1()->writeString(outstr1);
        if(!mOptions->out2.empty())
            config->getWriter2()->writeString(outstr2);
    } else {
        // normal output by left/right writer thread
        if(mRightWriter && mLeftWriter) {
            // write PE
            char* ldata = new char[outstr1.size()];
            memcpy(ldata, outstr1.c_str(), outstr1.size());
            mLeftWriter->input(ldata, outstr1.size());

            char* rdata = new char[outstr2.size()];
            memcpy(rdata, outstr2.c_str(), outstr2.size());
            mRightWriter->input(rdata, outstr2.size());
        } else if(mLeftWriter) {
            // write interleaved
            char* ldata = new char[interleaved.size()];
            memcpy(ldata, interleaved.c_str(), interleaved.size());
            mLeftWriter->input(ldata, interleaved.size());
        }
    }
    if(!mOptions->split.enabled)
        mOutputMtx.unlock();

    if(mOptions->split.byFileLines)
        config->markProcessed(readPassed);
    else
        config->markProcessed(pack->count);

    //delete pack->data;
    delete pack;

    return true;
}
    
void PairEndProcessor::statInsertSize(Read* r1, Read* r2, OverlapResult& ov) {
    int isize = mOptions->insertSizeMax;
    if(ov.overlapped) {
        if(ov.offset > 0)
            isize = r1->length() + r2->length() - ov.overlap_len;
        else
            isize = ov.overlap_len;
    }

    if(isize > mOptions->insertSizeMax)
        isize = mOptions->insertSizeMax;

    mInsertSizeHist[isize]++;
}

bool PairEndProcessor::processRead(Read* r, ReadPair* originalPair, bool reversed) {
    // do something here
    return true;
}

void PairEndProcessor::initPackRepository() {
    //mRepo.packBuffer = new ReadPairPack*[PACK_NUM_LIMIT];
	if(mOptions->interleavedInput){
    	mRepo.packBufferInter = new dsrc::fq::FastqDataChunk * [PACK_NUM_LIMIT];
    	memset(mRepo.packBufferInter, 0, sizeof(dsrc::fq::FastqDataChunk*)*PACK_NUM_LIMIT);
	}
	else{
		mRepo.packBuffer = new ChunkPair*[PACK_NUM_LIMIT];
    	memset(mRepo.packBuffer, 0, sizeof(ChunkPair*)*PACK_NUM_LIMIT);
	}
    mRepo.writePos = 0;
    mRepo.readPos = 0;
    
}

void PairEndProcessor::destroyPackRepository() {
	if(mOptions->interleavedInput)
    	delete mRepo.packBufferInter;
	else
    	delete mRepo.packBuffer;

    mRepo.packBuffer = NULL;
}
void PairEndProcessor::producePack_interleaved(dsrc::fq::FastqDataChunk* pack){

    mRepo.packBufferInter[mRepo.writePos] = pack;
    mRepo.writePos++;
}

void PairEndProcessor::producePack(ChunkPair* pack){
    mRepo.packBuffer[mRepo.writePos] = pack;
    mRepo.writePos++;
}

void PairEndProcessor::consumePack(ThreadConfig* config){

	//pe single file
    dsrc::fq::FastqDataChunk* chunk;
	ChunkPair* chunkpair;
	ReadPairPack* data = new ReadPairPack;
	ReadPack* leftPack = new ReadPack;
	ReadPack* rightPack = new ReadPack;
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
	if(mOptions->interleavedInput)
    	chunk = mRepo.packBufferInter[mRepo.readPos];
	else
		chunkpair = mRepo.packBuffer[mRepo.readPos];
	//cerr << "readPos: " << mRepo.readPos << endl;
	//cerr << (char*)chunkpair->leftpart->data.Pointer();
	//cout << "chunk pair lsize:" << chunkpair->leftpart->size << endl;
    mRepo.readPos++;

    /*if (mRepo.readPos >= PACK_NUM_LIMIT)
        mRepo.readPos = 0;*/
    mInputMtx.unlock();
    //mRepo.readPos++;

    //lock.unlock();
    //mRepo.repoNotFull.notify_all();
	//not interleaved case:
	if(mOptions->interleavedInput){
    	data->count = dsrc::fq::pairedChunkFormat(chunk,data->data,true);
    	fastqPool->Release(chunk);
    	processPairEnd(data, config);
	}
	else{
		leftPack->count = dsrc::fq::chunkFormat(chunkpair->leftpart, leftPack->data, true);
		rightPack->count = dsrc::fq::chunkFormat(chunkpair->rightpart, rightPack->data, true);
		//if(leftPack->count != rightPack->count)
		//{
		//	cout << "read pair chunk error: count not equal!!" << leftPack->count << " " <<rightPack->count << endl;
		//	//exit(0); //-----------------TODO:exit()?????????----------------
		//	//std::terminate();
		//	exit(0);
		//}else{
		//	data->count = leftPack->count;
		//	for(int i = 0; i < leftPack->count; ++i){
		//		data->data.push_back(new ReadPair(leftPack->data[i], rightPack->data[i]));
		//	}
		//}

		//ignore the unpaired reads from the file tail
		data->count = leftPack->count < rightPack->count ? leftPack->count : rightPack->count;

		for(int i = 0; i < data->count; ++i){
			data->data.push_back(new ReadPair(leftPack->data[i], rightPack->data[i]));
		}

			
		//if(leftPack->count != rightPack->count)
		//	cerr << "read pair chunk error: count not equal!!" << leftPack->count << " " <<rightPack->count << endl;

		pairReader->fastqPool_left->Release(chunkpair->leftpart);
		pairReader->fastqPool_right->Release(chunkpair->rightpart);
		//cerr << "Read Pair data size is: " << data->data.size() << endl;	
    	processPairEnd(data, config);

		delete leftPack;
		delete rightPack;
	}

}

void PairEndProcessor::producerTask(){
    if(mOptions->verbose)
        loginfo("start to load data");
    long lastReported = 0;
    int slept = 0;
    long readNum = 0;
    bool splitSizeReEvaluated = false;
    if(mOptions->interleavedInput){
        dsrc::fq::FastqFileReader* fileReader = new dsrc::fq::FastqFileReader(mOptions->in1);
        dsrc::fq::FastqReader* dataReader = new dsrc::fq::FastqReader(*fileReader,*fastqPool);
        dsrc::fq::FastqDataChunk* chunk;

        while((chunk = dataReader ->readNextPairedChunk())!= NULL){
            producePack_interleaved(chunk);
            //std::cout<< "chunkSizeL:"<< chunk->size<<std::endl;
            while(mRepo.writePos - mRepo.readPos > PACK_IN_MEM_LIMIT){
                slept++;
                usleep(100);
            }

        }
        delete fileReader;
        delete dataReader;
    }else{
		pairReader = new FastqChunkReaderPair(mOptions->in1, mOptions->in2, true, mOptions->phred64, mOptions->interleavedInput);

		ChunkPair* chunk_pair; 
		//while((chunk_pair = pairReader->readNextChunkPair()) != NULL){
		while((chunk_pair = pairReader->readNextChunkPair()) != NULL){
			//cerr << (char*)chunk_pair->leftpart->data.Pointer();
			producePack(chunk_pair);
			while(mRepo.writePos - mRepo.readPos > PACK_IN_MEM_LIMIT){
				slept++;
				usleep(100);
			}
		}
	}

    mProduceFinished = true;
    if(mOptions->verbose)
        loginfo("all reads loaded, start to monitor thread status");
    //lock.unlock();

}


void PairEndProcessor::consumerTask(ThreadConfig* config)
{
	//std::cout << "in consumerTask " << std::endl;
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
        if(mRightWriter)
            mRightWriter->setInputCompleted();
    }
    
    if(mOptions->verbose) {
        string msg = "thread " + to_string(config->getThreadId() + 1) + " finished";
        loginfo(msg);
    }
}

void PairEndProcessor::writeTask(WriterThread* config)
{
    while(true) {
        if(config->isCompleted()){
            // last check for possible threading related issue
            config->output();
            break;
        }
        config->output();
    }

    if(mOptions->verbose) {
        string msg = config->getFilename() + " writer finished";
        loginfo(msg);
    }
}
