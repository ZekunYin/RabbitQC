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
#ifndef TGSSTATS_H
#define TGSSTATS_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <map>
#include "read.h"
#include <unordered_map>
//#include "options.h"

using namespace std;

class TGSStats{
public:
	TGSStats(int minLen);
	~TGSStats();
	static TGSStats* merge(vector<TGSStats*>& list);
	void print();
	void tgsStatRead(Read* r);
    
	void reportJson(ofstream& ofs, string padding);
    // a port of HTML report
    void reportHtml(ofstream& ofs, string filteringType, string readName);
    void reportHtmlQuality(ofstream& ofs, string seqFileName, bool isTail, string xAxisName, string yAxisName, double* statsData);
    void reportHtmlContents(ofstream& ofs, string seqFileName, bool isTail, string xAxisName, string yAxisName, double** statsData);
	void reportHistogram(ofstream& ofs);
    bool isLongRead();

    static string list2string(double* list, int size);
    static string list2string(double* list, int size, long* coords);
    static string list2string(long* list, int size);
	static string list2stringReversedOrder(long* list, int size); 
    int base2num(string base);

//private:
	int mMinlen;
	int mHalfMinlen;
	vector<int> mLengths;
	vector<int> mTotalReadsLen; 
	long *head_seq_pos_count[4];
	long *tail_seq_pos_count[4];
	long *head_qual_sum;
	long *tail_qual_sum;
};

#endif
