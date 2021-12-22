#include "particles.h"
#include <iostream>

ClassImp(particles);

particles::particles(){}	// Empty constructor
particles::~particles(){}	// Empty destructor

void particles::Clear(){
	PID		= 0;
	Charge		= 0;
	Status		= 0;

	Time		= 0;
	Beta		= 0;
	Chi2		= 0;
	Pindex		= -1;
	Vtx		= 0;
	Vty		= 0;
	Vtz		= 0;

	Momentum	= 0;
	Theta		= 0;
	Phi		= 0;

	Scint_detector	.clear();
	Scint_status	.clear();
	Scint_sector	.clear();
	Scint_layer	.clear();
	Scint_component	.clear();
	Scint_Edep	.clear();
	Scint_time	.clear();
	Scint_path	.clear();
	Scint_chi2	.clear();
	Scint_x		.clear();
	Scint_y		.clear();
	Scint_z		.clear();
}


void particles::Print(){
	return;
}
