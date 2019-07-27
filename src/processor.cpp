#include "processor.h"
#include "peprocessor.h"
#include "seprocessor.h"
#include "TGSProcessor.h"

Processor::Processor(Options* opt){
    mOptions = opt;
}


Processor::~Processor(){
}

bool Processor::process() {
    if(mOptions->isPaired()) {
        PairEndProcessor p(mOptions);
        p.process();
    }
	else if(mOptions->isTGS()){
		TGSProcessor p(mOptions);
		p.process();
	}
	else {
        SingleEndProcessor p(mOptions);
        p.process();
    }

    return true;
}
