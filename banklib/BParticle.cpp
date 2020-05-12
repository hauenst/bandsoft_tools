#include "BParticle.h"
#include "TVector3.h"
// ==============================================================
BParticle::BParticle(hipo::schema __schema) : hipo::bank(__schema){

	pid_order     = __schema.getEntryOrder("pid"    );
	px_order      = __schema.getEntryOrder("px"     );
	py_order      = __schema.getEntryOrder("py"     );
	pz_order      = __schema.getEntryOrder("pz"     );
	vx_order      = __schema.getEntryOrder("vx"     );
	vy_order      = __schema.getEntryOrder("vy"     );
	vz_order      = __schema.getEntryOrder("vz"     );
	vt_order      = __schema.getEntryOrder("vt"	);
	charge_order  = __schema.getEntryOrder("charge" );
	beta_order    = __schema.getEntryOrder("beta"   );
	chi2pid_order = __schema.getEntryOrder("chi2pid");
	status_order  = __schema.getEntryOrder("status" );
}
// ==============================================================
BParticle::~BParticle(){}
// ==============================================================
TVector3 BParticle::getV3v(int index)
{
	TVector3 * V3v = new TVector3();
        float vx = BParticle::getVx(index);
        float vy = BParticle::getVy(index);
        float vz = BParticle::getVz(index);
        V3v -> SetXYZ(vx,vy,vz);
        return * V3v;
}
// ==============================================================
TVector3 BParticle::getV3P(int index)
{
	TVector3 * V3P = new TVector3();
	float px = BParticle::getPx(index);
	float py = BParticle::getPy(index);
	float pz = BParticle::getPz(index);
	V3P -> SetXYZ(px,py,pz);
	return * V3P;
}
// ==============================================================
