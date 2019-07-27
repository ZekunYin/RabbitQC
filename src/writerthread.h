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
#ifndef WRITER_THREAD_H
#define WRITER_THREAD_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include "writer.h"
#include "options.h"
#include <atomic>
#include <mutex>

using namespace std;

class WriterThread{
public:
    WriterThread(Options* opt, string filename);
    ~WriterThread();

    void initWriter(string filename1);
    void initWriter(ofstream* stream);
    void initWriter(gzFile gzfile);

    void cleanup();

    bool isCompleted();
    void output();
    void input(char* data, size_t size);
    bool setInputCompleted();

    long bufferLength();
    string getFilename() {return mFilename;}

private:
    void deleteWriter();

private:
    Writer* mWriter1;
    Options* mOptions;
    string mFilename;

    // for spliting output
    bool mInputCompleted;
    atomic_long mInputCounter;
    atomic_long mOutputCounter;
    char** mRingBuffer;
    size_t* mRingBufferSizes;

    mutex mtx;

};

#endif
