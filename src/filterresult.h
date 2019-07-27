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
#ifndef FILTER_RESULT_H
#define FILTER_RESULT_H

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include "common.h"
#include "options.h"
#include <fstream>
#include <map>

struct classcomp {
    bool operator() (const string& lhs, const string& rhs) const {
        if (lhs.length() < rhs.length())
            return true;
        else if(lhs.length() == rhs.length()) {
            return lhs < rhs;
        } else
            return false;
    }
};

using namespace std;

class FilterResult{
public:
    FilterResult(Options* opt, bool paired = false);
    ~FilterResult();
    inline long* getFilterReadStats() {return mFilterReadStats;}
    void addFilterResult(int result);
    static FilterResult* merge(vector<FilterResult*>& list);
    void print();
    // for single end
    void addAdapterTrimmed(string adapter, bool isR2 = false);
    // for paired end
    void addAdapterTrimmed(string adapter1, string adapter2);
    // a part of JSON report
    void reportJson(ofstream& ofs, string padding);
    // a part of JSON report for adapters
    void reportAdapterJson(ofstream& ofs, string padding);
    // a part of HTML report
    void reportHtml(ofstream& ofs, long totalReads, long totalBases);
    // a part of HTML report for adapters
    void reportAdapterHtml(ofstream& ofs, long totalBases);
    void outputAdaptersJson(ofstream& ofs, map<string, long, classcomp>& adapterCounts);
    void outputAdaptersHtml(ofstream& ofs, map<string, long, classcomp>& adapterCounts, long totalBases);
    // deal with base correction results
    long* getCorrectionMatrix() {return mCorrectionMatrix;}
    long getTotalCorrectedBases();
    void addCorrection(char from, char to);
    long getCorrectionNum(char from, char to);
    void incCorrectedReads(int count);

public:
    Options* mOptions;
    bool mPaired;
    long mCorrectedReads;
private:
    long mFilterReadStats[FILTER_RESULT_TYPES];
    long mTrimmedAdapterRead;
    long mTrimmedAdapterBases;
    map<string, long, classcomp> mAdapter1;
    map<string, long, classcomp> mAdapter2;
    long* mCorrectionMatrix;
};

#endif
