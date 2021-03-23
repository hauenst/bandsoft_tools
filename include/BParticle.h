#ifndef BPARTICLE_H
#define BPARTICLE_H

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <vector>
#include <stdint.h>
#include <stdlib.h>
#include "bank.h"

#include "TVector3.h"

class BParticle : public hipo::bank {

	private:

		int pid_order     ;
		int px_order      ;
		int py_order      ;
		int pz_order      ;
		int vx_order      ;
		int vy_order      ;
		int vz_order      ;
		int vt_order      ;
		int charge_order  ;
		int beta_order    ;
		int chi2pid_order ;
		int status_order  ;

	public:

		BParticle(){};

		BParticle(hipo::schema __schema);

		~BParticle();

		int   getPid    (int row) { return getInt  (pid_order    ,row);}
		float getPx     (int row) { return getFloat(px_order     ,row);}
		float getPy     (int row) { return getFloat(py_order     ,row);}
		float getPz     (int row) { return getFloat(pz_order     ,row);}
		float getVx     (int row) { return getFloat(vx_order     ,row);}
		float getVy     (int row) { return getFloat(vy_order     ,row);}
		float getVz     (int row) { return getFloat(vz_order     ,row);}
		float getVt     (int row) { return getFloat(vt_order     ,row);}
		int   getCharge (int row) { return getInt  (charge_order ,row);}
		float getBeta   (int row) { return getFloat(beta_order   ,row);}
		float getChi2pid(int row) { return getFloat(chi2pid_order,row);}
		int   getStatus (int row) { return getInt  (status_order ,row);}

		TVector3 getV3v(int row);
		TVector3 getV3P(int row);
};

#endif
