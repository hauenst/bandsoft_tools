#include "bandhit.h"

ClassImp(bandhit);

bandhit::bandhit(){}	// Empty constructor
bandhit::~bandhit(){}	// Empty destructor

void bandhit::Clear(){
	Sector		= 0;
	Layer		= 0;
	Component	= 0;
	BarID		= 0;
	Edep		= 0;
	Tof		= 0;
	TofFadc		= 0;
	Tdiff		= 0;
	TdiffFadc	= 0;
	X		= 0;
	Y		= 0;
	Z		= 0;
		
	RawLtdc		= 0;
	RawRtdc		= 0;
	RawLtdccorr	= 0;
	RawRtdccorr	= 0;
	RawLtfadc	= 0;
	RawRtfadc	= 0;
	RawLamp		= 0;
	RawRamp		= 0;
	RawLadc		= 0;
	RawRadc		= 0;
		
	PmtLtdc		= 0;
	PmtRtdc		= 0;
	PmtLtfadc	= 0;
	PmtRtfadc	= 0;
	PmtLamp		= 0;
	PmtRamp		= 0;
	PmtLadc		= 0;
	PmtRadc		= 0;
	PmtLped		= 0;
	PmtRped		= 0;
}

