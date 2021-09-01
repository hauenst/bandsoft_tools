#include "e_pid.h"

#include "TF1.h"
#include "TCanvas.h"
#include "TString.h"
#include "TMath.h"

#include "clashit.h"
#include <iostream>
#include <fstream>
#include <unistd.h>

#define PID_DIR _VAR

using namespace std;

e_pid::e_pid(){

  //Default is 10.6GeV parameters for SF cuts
  setParamsRGB(10.6);
  intervalEpcal = 3.5;
  intervalMom = 3.5;
  color = 2;

}

e_pid::~e_pid(){
}


bool e_pid::isElectron(clashit * eHit){

  //if any of the energy values is 0, cut away the electron
  if (eHit->getEpcal() == 0 || eHit->getEoP() == 0 || eHit->getEecin() == 0) {
    return false;
  }

  bool passSFEpcal = SF_Epcal_Cut(eHit->getSector()-1,eHit->getEpcal(),eHit->getEoP());
  bool passSFMom = SF_Mom_Cut(eHit->getSector()-1,eHit->getMomentum(),eHit->getEoP());
  bool passSFpi = SFpcal_SFecin_Cut(eHit->getMomentum(),eHit->getEpcal(),eHit->getEecin());

  //For Debugging
  //std::cout << "ePID: passSFEpcal = " << passSFEpcal << " , passSFMom = " << passSFMom << " , passSFpi = " << passSFpi << std::endl;


  if(eHit->getPID() != 11){ return false; }
  if(eHit->getCharge() != -1){ return false; }
  if(!passSFEpcal){ return false; }
  if(!passSFMom){ return false; }
  if(!passSFpi){ return false; }
  return true;

}
bool e_pid::isElectronLoose(clashit * eHit){

  bool passSFEpcal = SF_Epcal_Cut(eHit->getSector()-1,eHit->getEpcal(),eHit->getEoP());
  bool passSFMom = SF_Mom_Cut(eHit->getSector()-1,eHit->getMomentum(),eHit->getEoP());
  bool passSFpi = SFpcal_SFecin_Cut(eHit->getMomentum(),eHit->getEpcal(),eHit->getEecin());

  if(eHit->getPID() != 11){ return false; }
  if(eHit->getCharge() != -1){ return false; }
  if(!passSFEpcal){ return false; }
  return true;

}

void e_pid::setParamsRGB(double Ebeam){

  //Determine correct file name
  if(std::fabs(Ebeam-10.6) < 0.05){
    sprintf(paramFileNameEpcal,"%s/SFvEpcal_Params_106_RGB.dat",std::string(PID_DIR).c_str());
    sprintf(paramFileNameMom,"%s/SFvMom_Params_106_RGB.dat",std::string(PID_DIR).c_str());
  }
  else if(std::fabs(Ebeam-10.2) < 0.05){
    sprintf(paramFileNameEpcal,"%s/SFvEpcal_Params_102_RGB.dat",std::string(PID_DIR).c_str());
    sprintf(paramFileNameMom,"%s/SFvMom_Params_102_RGB.dat",std::string(PID_DIR).c_str());
  }
  else if(std::fabs(Ebeam-4.247) < 0.2){//10.2 seems to work well for the low energy run data set
    sprintf(paramFileNameEpcal,"%s/SFvEpcal_Params_102_RGB.dat",std::string(PID_DIR).c_str());
    sprintf(paramFileNameMom,"%s/SFvMom_Params_102_RGB.dat",std::string(PID_DIR).c_str());
  }
  else{
    std::cout<<"Attempting to set a beam energy "<< Ebeam <<" GeV\n"
	     <<"without a defined RGB fiducal cut\n"
	     <<"\t Using E=10.6 GeV as Default. \n\n\n";
    sprintf(paramFileNameEpcal,"%s/SFvEpcal_Params_106_RGB.dat",std::string(PID_DIR).c_str());
    sprintf(paramFileNameMom,"%s/SFvMom_Params_106_RGB.dat",std::string(PID_DIR).c_str());
  }

  //Load file data
  fillParams();
}

void e_pid::fillParams(){

  //Load file data
  ifstream paramFileEpcal(paramFileNameEpcal);
  ifstream paramFileMom(paramFileNameMom);
  if(!paramFileEpcal.is_open() || !paramFileMom.is_open()){
    cout<<"Fiducial cut parameter files failed to load.\n"
	<<"Aborting...\n\n\n";
    exit(-2);
  }

  for(int i = 0; i < 6; i++){
    for(int j = 0; j < 6; j++){
      paramsEpcal[i][j] = 0;
      paramFileEpcal >> paramsEpcal[i][j];
      paramsMom[i][j] = 0;
      paramFileMom >> paramsMom[i][j];
    }
  }

  paramFileEpcal.close();
  paramFileMom.close();

}

void e_pid::setIntervalEpcal(double newInterval){
  intervalEpcal = newInterval;
}

void e_pid::setIntervalMom(double newInterval){
  intervalMom = newInterval;
}

void e_pid::setColor(int newColor){
  color = newColor;
}

void e_pid::drawEpcal(int sector, TCanvas * myCanvas){

  sector--;

  TF1 * meanFunction = new TF1("Mean",[&](double *x, double *p){ return FF(x[0],p[0],p[1],p[2]); },0.06,1.6,3);
  meanFunction->SetLineColor(color);
  meanFunction->SetParameters(paramsEpcal[sector][0],paramsEpcal[sector][1],paramsEpcal[sector][2]);

  TF1 * maxFunction = new TF1("Max",[&](double *x, double *p){ return FF(x[0],p[0],p[1],p[2]) + intervalEpcal * FF(x[0],p[3],p[4],p[5]); },0.06,1.6,6);
  maxFunction->SetLineColor(color);
  maxFunction->SetParameters(paramsEpcal[sector][0],paramsEpcal[sector][1],paramsEpcal[sector][2],paramsEpcal[sector][3],paramsEpcal[sector][4],paramsEpcal[sector][5]);

  TF1 * minFunction = new TF1("Min",[&](double *x, double *p){ return FF(x[0],p[0],p[1],p[2]) - intervalEpcal * FF(x[0],p[3],p[4],p[5]); },0.06,1.6,6);
  minFunction->SetLineColor(color);
  minFunction->SetParameters(paramsEpcal[sector][0],paramsEpcal[sector][1],paramsEpcal[sector][2],paramsEpcal[sector][3],paramsEpcal[sector][4],paramsEpcal[sector][5]);


  myCanvas->cd();
  meanFunction->Draw("SAME");
  maxFunction->Draw("SAME");
  minFunction->Draw("SAME");

}

void e_pid::drawMom(int sector, TCanvas * myCanvas){

  sector--;

  TF1 * meanFunction = new TF1("Mean",[&](double *x, double *p){ return FF(x[0],p[0],p[1],p[2]); },0.5,9.5,3);
  meanFunction->SetLineColor(color);
  meanFunction->SetParameters(paramsMom[sector][0],paramsMom[sector][1],paramsMom[sector][2]);

  TF1 * maxFunction = new TF1("Max",[&](double *x, double *p){ return FF(x[0],p[0],p[1],p[2]) + intervalEpcal * FF(x[0],p[3],p[4],p[5]); },0.5,9.5,6);
  maxFunction->SetLineColor(color);
  maxFunction->SetParameters(paramsMom[sector][0],paramsMom[sector][1],paramsMom[sector][2],paramsMom[sector][3],paramsMom[sector][4],paramsMom[sector][5]);

  TF1 * minFunction = new TF1("Min",[&](double *x, double *p){ return FF(x[0],p[0],p[1],p[2]) - intervalEpcal * FF(x[0],p[3],p[4],p[5]); },0.5,9.5,6);
  minFunction->SetLineColor(color);
  minFunction->SetParameters(paramsMom[sector][0],paramsMom[sector][1],paramsMom[sector][2],paramsMom[sector][3],paramsMom[sector][4],paramsMom[sector][5]);


  myCanvas->cd();
  meanFunction->Draw("SAME");
  maxFunction->Draw("SAME");
  minFunction->Draw("SAME");

}

bool e_pid::SF_Epcal_Cut(int sector, double Epcal, double SF){
  return SF_Cut(sector,Epcal,SF,paramsEpcal,intervalEpcal);
}

bool e_pid::SF_Mom_Cut(int sector, double p, double SF){
  return SF_Cut(sector,p,SF,paramsMom,intervalMom);
}

bool e_pid::SFpcal_SFecin_Cut(double p, double Epcal, double Eecin){

  //This cut is only for high momentum particles
  if(p < 4.5){ return true; }

  double SFpi = (Epcal + Eecin) / p;
  if(SFpi < 0.2){ return false; }
  return true;

}

bool e_pid::SF_Cut(int sector, double x, double SF, double params[6][6], double interval){
  //cout << "in SF_cut: SF =  " << SF << " , minSF = " <<minSF(sector,x,params,interval) << " , maxSF " << maxSF(sector,x,params,interval) << endl;
  if(SF < minSF(sector,x,params,interval)){
    return false;
  }
  if(SF > maxSF(sector,x,params,interval)){
    return false;
  }
  return true;
}

double e_pid::meanSF(int sector, double x, double params[6][6]){
  return FF(x,params[sector][0],params[sector][1],params[sector][2]);
}

double e_pid::sigmaSF(int sector, double x, double params[6][6]){
  return FF(x,params[sector][3],params[sector][4],params[sector][5]);
}

double e_pid::maxSF(int sector, double x, double params[6][6], double interval){
  return meanSF(sector,x,params) + interval * sigmaSF(sector,x,params);
}

double e_pid::minSF(int sector, double x, double params[6][6], double interval){
  return meanSF(sector,x,params) - interval * sigmaSF(sector,x,params);
}

double e_pid::FF(double x, double sf1, double sf2, double sf3){
  return sf1 + (sf2/sqrt(x)) + (sf3/x);
}
