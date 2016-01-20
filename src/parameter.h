#ifndef PARAMETER_H
#define PARAMETER_H

#include <string>
using namespace std;

//file path for results
const string results_path = "results/";

const double e_sq = 1./137*4*M_PI;
const double hbarC = 0.197327053;  // GeV*fm
const double unit_convert = hbarC*hbarC*10e6; // convert 1/(fm^2) to MeV^2

//spacial grid information
const int nx = 11;
const int ny = 11;
const int nz = 1;
const int nt = 20;

const double t_i = 0.0;
const double t_f = 2.0;
const double x_i = -5.0;
const double x_f = 5.0;
const double y_i = -5.0;
const double y_f = 5.0;
const double z_i = 0.0;
const double z_f = 0.0;

//gaussian point for integration
const int n_r = 10;
const int n_phi = 20;

//type of the nucleon
const int Nucleon_type = 2;  //0: point-like. 1: disk-like. 2: gaussian-like
const double disk_radius = 0.4; // disk area for disk-like nucleon
const double norm_disk = 1./(M_PI*disk_radius*disk_radius); //normalization for disk-like nucleon
const double sigma_gaussian = 0.4; // width for gaussian-like nucleon
const double norm_gaussian = 1./(2.*M_PI*sigma_gaussian*sigma_gaussian);  //normalization for gaussian-like nucleon
const int gaussian_Nsigma = 5; //integrated over Nsigma sigma region

#endif
