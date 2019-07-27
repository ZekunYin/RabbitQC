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
#ifndef _WRITER_H
#define _WRITER_H

#include <stdio.h>
#include <stdlib.h>
#ifdef DYNAMIC_ZLIB
  #include <zlib.h>
#else
  #include "zlib/zlib.h"
#endif
#include "common.h"
#include <iostream>
#include <fstream>

using namespace std;

class Writer{
public:
	Writer(string filename, int compression = 3);
	Writer(ofstream* stream);
	Writer(gzFile gzfile);
	~Writer();
	bool isZipped();
	bool writeString(string& s);
	bool writeLine(string& linestr);
	bool write(char* strdata, size_t size);
	string filename();

public:
	static bool test();

private:
	void init();
	void close();

private:
	string mFilename;
	gzFile mZipFile;
	ofstream* mOutStream;
	bool mZipped;
	int mCompression;
	bool haveToClose;
};

#endif
