#include "clashit.h"
#include <iostream>

ClassImp(clashit);

clashit::clashit(){}	// Empty constructor
clashit::~clashit(){}	// Empty destructor

void clashit::Clear(){
	Sector		= -1;
	PID		= 0;
	Charge		= 0;
	Status		= 0;

	Time		= 0;
	Beta		= 0;
	Chi2		= 0;
	Etot		= 0;
	Epcal		= 0;
	Eecin		= 0;
	Eecout		= 0;
	EoP		= 0;
	TimeScint	= 0;
	PathScint	= 0;
	U		= 0;
	V		= 0;
	W		= 0;
	PCal_X		= 0;
	PCal_Y		= 0;
	PCal_Z		= 0;
	Vtx		= 0;
	Vty		= 0;
	Vtz		= 0;

	Momentum	= 0;
	Theta		= 0;
	Phi		= 0;

	Q		= 0;
	ThetaQ		= 0;
	PhiQ		= 0;

	Q2		= 0;
	Omega		= 0;
	Xb		= 0;
	W2		= 0;

	DC_chi2         = -999;
	DC_NDF          = -999;
	DC_sector       = -999;

	DC_x1           = -999;
	DC_y1           = -999;
	DC_z1           = -999;

	DC_x2           = -999;
	DC_y2           = -999;
	DC_z2           = -999;

	DC_x3           = -999;
	DC_y3           = -999;
	DC_z3           = -999;

	Nphe		= -999;
	Kov_x		= -999;
	Kov_y		= -999;
	Kov_z		= -999;
	Kov_chi2	= -999;
	Kov_time	= -999;
	Kov_path	= -999;
	Kov_det		= -999;
	Kov_sec		= -999;
	Kov_status	= -999;

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


void clashit::Print(){

	std::cout << "clashit Information REC: PID " << PID << " , Charge " << Charge << " , Status " << Status;
	std::cout << ", Sector(Calo) " << Sector << " , Chi2 " << Chi2 << ", Time " << Time << " , Beta " << Beta << std::endl;
	std::cout << "clashit Information CALO: Etot " << Etot << " , Epcal " << Epcal << ", Eecin " << Eecin;
	std::cout << ", Eecout " << Eecout << ", EoverP " << EoP << " , U(PCal) " << U << " , V(ECal) " << V << " , W(ECal) " << W << std::endl;
	std::cout << "clashit Information Vertex: Vtx " << Vtx << " , Vty " << Vty << ", Vtz " << Vtz;
	std::cout << ", TimeScint(-starttime) " << TimeScint << " , Pathlength Scint " << PathScint << std::endl;
	std::cout << "clashit Information DC: DC_chi2 " << DC_chi2 << " , DC_NDF " << DC_NDF << ", DC_sector " << DC_sector;
	std::cout << ", DC_x1 " << DC_x1 << " , DC_y1 " << DC_y1 << " , DC_z1 " << DC_z1 << std::endl;
	std::cout << "clashit Information DC: DC_x2 " << DC_x2 << " , DC_y2 " << DC_y2<< " , DC_z2 " << DC_z2;
	std::cout << ", DC_x3 " << DC_x3 << " , DC_y3 " << DC_y3 << " , DC_z3 " << DC_z3 << std::endl;
	std::cout << "clashit Information Kinematics: Momentum " << Momentum << " , Theta " << Theta << " ,Phi " << Phi;
	std::cout << ", Q2 " << Q2 << " , Omega/nu " << Omega << " , Xb " << Xb << " , W2 " << W2 << std::endl;
	std::cout << "clashit Information q-vector: Magitude(q) " << Q << " , ThetaQ " << ThetaQ << " ,PhiQ " << PhiQ << std::endl;

}
