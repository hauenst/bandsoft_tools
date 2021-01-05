#ifndef __DC_FIDUCIAL_H__
#define __DC_FIDUCIAL_H__

#include "TF1.h"
#include "TVector3.h"

class DCFiducial {
 public:
  DCFiducial();
  virtual ~DCFiducial();

  TVector3 rotate(double x, double y, int sector) const;
  bool DC_e_fid(double x, double y, int sector, int layer, bool bending) const;

  TF1 *get_fmin_in(int sector, int layer) const;
  TF1 *get_fmax_in(int sector, int layer) const;
  TF1 *get_fmin_out(int sector, int layer) const;
  TF1 *get_fmax_out(int sector, int layer) const;
  
 
  
 protected:
  TF1 *fmin_S_L_in[6][3];
  TF1 *fmax_S_L_in[6][3];
  TF1 *fmin_S_L_out[6][3];
  TF1 *fmax_S_L_out[6][3];

  static const double maxparams_in[6][6][3][2];
  static const double minparams_in[6][6][3][2];
  static const double maxparams_out[6][6][3][2];
  static const double minparams_out[6][6][3][2];  
};

#endif
