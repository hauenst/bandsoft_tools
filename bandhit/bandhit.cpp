#include "bandhit.h"

ClassImp(bandhit);

bandhit::bandhit(){}	// Empty constructor
bandhit::~bandhit(){}	// Empty destructor

void bandhit::Print(){
	std::cout<< "Sector		" << Sector		<< "\n";
	std::cout<< "Layer		" << Layer		<< "\n";
	std::cout<< "Component	" << Component	<< "\n";
	std::cout<< "BarID		" << BarID		<< "\n";
	std::cout<< "Edep		" << Edep		<< "\n";
	std::cout<< "Tof		" << Tof		<< "\n";
	std::cout<< "TofFadc		" << TofFadc		<< "\n";
	std::cout<< "Tdiff		" << Tdiff		<< "\n";
	std::cout<< "TdiffFadc	" << TdiffFadc	<< "\n";
	std::cout<< "X		" << X		<< "\n";
	std::cout<< "XFadc		" << XFadc		<< "\n";
	std::cout<< "Y		" << Y		<< "\n";
	std::cout<< "Z		" << Z		<< "\n";
	std::cout<< "Status		" << Status		<< "\n";
	std::cout<< "RawLtdc		" << RawLtdc		<< "\n";
	std::cout<< "RawRtdc		" << RawRtdc		<< "\n";
	std::cout<< "RawLtdccorr	" << RawLtdccorr	<< "\n";
	std::cout<< "RawRtdccorr	" << RawRtdccorr	<< "\n";
	std::cout<< "RawLtfadc	" << RawLtfadc	<< "\n";
	std::cout<< "RawRtfadc	" << RawRtfadc	<< "\n";
	std::cout<< "RawLamp		" << RawLamp		<< "\n";
	std::cout<< "RawRamp		" << RawRamp		<< "\n";
	std::cout<< "RawLadc		" << RawLadc		<< "\n";
	std::cout<< "RawRadc		" << RawRadc		<< "\n";
	std::cout<< "PmtLtdc		" << PmtLtdc		<< "\n";
	std::cout<< "PmtRtdc		" << PmtRtdc		<< "\n";
	std::cout<< "PmtLtfadc	" << PmtLtfadc	<< "\n";
	std::cout<< "PmtRtfadc	" << PmtRtfadc	<< "\n";
	std::cout<< "PmtLamp		" << PmtLamp		<< "\n";
	std::cout<< "PmtRamp		" << PmtRamp		<< "\n";
	std::cout<< "PmtLadc		" << PmtLadc		<< "\n";
	std::cout<< "PmtRadc		" << PmtRadc		<< "\n";
	std::cout<< "PmtLped		" << PmtLped		<< "\n";
	std::cout<< "PmtRped		" << PmtRped		<< "\n";

	return;
}

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
	XFadc		= 0;
	Y		= 0;
	Z		= 0;
	Status		= 0;
	DL		.Clear();
		
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

