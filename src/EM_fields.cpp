#include<iostream>
#include<sstream>
#include<fstream>
#include<cmath>
#include<iomanip>

#include "EM_fields.h"
#include "Spectators.h"
#include "gauss_quadrature.h"

using namespace std;

EM_fields::EM_fields()
{

   E_x = new double*** [nt];
   E_y = new double*** [nt];
   E_z = new double*** [nt];
   B_x = new double*** [nt];
   B_y = new double*** [nt];
   B_z = new double*** [nt];
   for(int i=0; i<nt; i++)
   {
      E_x[i] = new double** [nx];
      E_y[i] = new double** [nx];
      E_z[i] = new double** [nx];
      B_x[i] = new double** [nx];
      B_y[i] = new double** [nx];
      B_z[i] = new double** [nx];
      for(int j=0; j<nx; j++)
      {
         E_x[i][j] = new double* [ny];
         E_y[i][j] = new double* [ny];
         E_z[i][j] = new double* [ny];
         B_x[i][j] = new double* [ny];
         B_y[i][j] = new double* [ny];
         B_z[i][j] = new double* [ny];
         for(int k=0; k<ny; k++)
         {
            E_x[i][j][k] = new double [nz];
            E_y[i][j][k] = new double [nz];
            E_z[i][j][k] = new double [nz];
            B_x[i][j][k] = new double [nz];
            B_y[i][j][k] = new double [nz];
            B_z[i][j][k] = new double [nz];
         }
      }
   }
   for(int i=0; i<nt; i++)
      for(int j=0; j<nx; j++)
         for(int k=0; k<ny; k++)
            for(int l=0; l<nz; l++)
            {
               E_x[i][j][k][l] = 0.0e0;
               E_y[i][j][k][l] = 0.0e0;
               E_z[i][j][k][l] = 0.0e0;
               B_x[i][j][k][l] = 0.0e0;
               B_y[i][j][k][l] = 0.0e0;
               B_z[i][j][k][l] = 0.0e0;
            }
   t = new double [nt];
   x = new double [nx];
   y = new double [ny];
   z = new double [nz];
   double eps = 1e-100;
   double dt = (t_f - t_i)/(nt - 1 + eps);
   double dx = (x_f - x_i)/(nx - 1 + eps);
   double dy = (y_f - y_i)/(ny - 1 + eps);
   double dz = (z_f - z_i)/(nz - 1 + eps);
   for(int i=0; i<nt; i++) t[i] = t_i + i*dt;
   for(int i=0; i<nx; i++) x[i] = x_i + i*dx;
   for(int i=0; i<ny; i++) y[i] = y_i + i*dy;
   for(int i=0; i<nz; i++) z[i] = z_i + i*dz;

   r = new double [n_r];
   r_weight = new double [n_r];
   phi = new double [n_phi];
   phi_weight = new double [n_phi];
   cos_phi = new double [n_phi];
   sin_phi = new double [n_phi];
   gauss_quadrature(n_phi, 1, 0.0, 0.0, 0.0, 2*M_PI, phi, phi_weight); 
   for(int i=0; i<n_phi; i++)
   {
      cos_phi[i] = cos(phi[i]);
      sin_phi[i] = sin(phi[i]);
   }
   switch(Nucleon_type)
   {
      case 1:
         gauss_quadrature(n_r, 1, 0.0, 0.0, 0.0, disk_radius, r, r_weight); 
         break;
      case 2:
         gauss_quadrature(n_r, 1, 0.0, 0.0, 0.0, gaussian_Nsigma*sigma_gaussian, r, r_weight); 
         break;
      default:
         break;
   }
   return;
}

EM_fields::~EM_fields()
{
   for(int i=0; i<nt; i++)
   {
      for(int j=0; j<nx; j++)
      {
         for(int k=0; k<ny; k++)
         {
             delete[] E_x[i][j][k];
             delete[] E_y[i][j][k];
             delete[] E_z[i][j][k];
             delete[] B_x[i][j][k];
             delete[] B_y[i][j][k];
             delete[] B_z[i][j][k];
         }
         delete[] E_x[i][j];
         delete[] E_y[i][j];
         delete[] E_z[i][j];
         delete[] B_x[i][j];
         delete[] B_y[i][j];
         delete[] B_z[i][j];
      }
      delete[] E_x[i];
      delete[] E_y[i];
      delete[] E_z[i];
      delete[] B_x[i];
      delete[] B_y[i];
      delete[] B_z[i];
   }
   delete[] E_x;
   delete[] E_y;
   delete[] E_z;
   delete[] B_x;
   delete[] B_y;
   delete[] B_z;
   delete[] t;
   delete[] x;
   delete[] y;
   delete[] z;
   delete[] r;
   delete[] r_weight;
   delete[] phi;
   delete[] phi_weight;
   delete[] cos_phi;
   delete[] sin_phi;

   return;
}

void EM_fields::calculate_EM_fields(Spectators nucleon_list)
{
   int list_length = nucleon_list.get_Spectator_list_length();
   double* nucleon_position_x = new double [list_length];
   double* nucleon_position_y = new double [list_length];
   double* nucleon_rapidity = new double [list_length];
   double* nucleon_cosh_y = new double [list_length];
   double* nucleon_sinh_y = new double [list_length];
   for(int i=0; i<list_length; i++)
   {
      nucleon_position_x[i] = nucleon_list.get_Spectator_position_x(i);
      nucleon_position_y[i] = nucleon_list.get_Spectator_position_y(i);
      nucleon_rapidity[i] = nucleon_list.get_Spectator_rapidity(i);
      nucleon_sinh_y[i] = sinh(nucleon_rapidity[i]);
      nucleon_cosh_y[i] = cosh(nucleon_rapidity[i]);
   }

   for(int i=0; i<nt; i++)
      for(int j=0; j<nx; j++)
         for(int k=0; k<ny; k++)
            for(int l=0; l<nz; l++)
            {
               double temp_sum_Ex = 0.0e0;
               double temp_sum_Ey = 0.0e0;
               double temp_sum_Ez = 0.0e0;
               double temp_sum_Bx = 0.0e0;
               double temp_sum_By = 0.0e0;
               double temp_sum_Bz = 0.0e0;
               for(int list_idx=0; list_idx < list_length; list_idx++)
               {
                  double x_local, y_local, z_local;
                  x_local = x[j] - nucleon_position_x[list_idx];
                  y_local = y[k] - nucleon_position_y[list_idx];
                  z_local = z[l]*nucleon_cosh_y[list_idx] - t[i]*nucleon_sinh_y[list_idx];
                  double denominator, denominator_cubic;
                  double Ex_integrand, Ey_integrand, Ez_integrand, Bx_integrand, By_integrand, Bz_integrand;
                  double Ex_integral, Ey_integral, Ez_integral, Bx_integral, By_integral, Bz_integral;
                  switch(Nucleon_type)
                  {
                     case 0: //point-like nucleons
                       denominator = sqrt(x_local*x_local + y_local*y_local + z_local*z_local);
                       denominator_cubic = denominator*denominator*denominator;
                       Ex_integral = nucleon_cosh_y[list_idx]*x_local/denominator_cubic;
                       Ey_integral = nucleon_cosh_y[list_idx]*y_local/denominator_cubic;
                       Ez_integral = z_local/denominator_cubic;
                       Bx_integral = nucleon_sinh_y[list_idx]*y_local/denominator_cubic;
                       By_integral = - nucleon_sinh_y[list_idx]*x_local/denominator_cubic;
                       Bz_integral = 0.0e0;
                       break;
                     case 1: //disk-like nucleons
                       for(int i_r=0; i_r<n_r; i_r++)
                          for(int i_phi=0; i_phi<n_phi; i_phi++)
                          {
                             denominator = sqrt(x_local*x_local + y_local*y_local + r[i_r]*r[i_r] - 2*r[i_r]*(cos_phi[i_phi]*x_local+sin_phi[i_phi]*y_local) + z_local*z_local);
                             denominator_cubic = denominator*denominator*denominator;
                             Ex_integrand = norm_disk*nucleon_cosh_y[list_idx]*(x_local - r[i_r]*cos_phi[i_phi])/denominator_cubic;
                             Ey_integrand = norm_disk*nucleon_cosh_y[list_idx]*(y_local - r[i_r]*sin_phi[i_phi])/denominator_cubic;
                             Ez_integrand = norm_disk*z_local/denominator_cubic;
                             Bx_integrand = norm_disk*nucleon_sinh_y[list_idx]*(y_local - r[i_r]*sin_phi[i_phi])/denominator_cubic;
                             By_integrand = - norm_disk*nucleon_sinh_y[list_idx]*(x_local - r[i_r]*cos_phi[i_phi])/denominator_cubic;
                             Bz_integrand = 0.0e0;
                             
                             Ex_integral += Ex_integrand*phi_weight[i_phi]*r_weight[i_r]*r[i_r];
                             Ey_integral += Ey_integrand*phi_weight[i_phi]*r_weight[i_r]*r[i_r];
                             Ez_integral += Ez_integrand*phi_weight[i_phi]*r_weight[i_r]*r[i_r];
                             Bx_integral += Bx_integrand*phi_weight[i_phi]*r_weight[i_r]*r[i_r];
                             By_integral += By_integrand*phi_weight[i_phi]*r_weight[i_r]*r[i_r];
                             Bz_integral += Bz_integrand*phi_weight[i_phi]*r_weight[i_r]*r[i_r];
                          }
                       break;
                     case 2: //gaussian-like nucleons
                       for(int i_r=0; i_r<n_r; i_r++)
                          for(int i_phi=0; i_phi<n_phi; i_phi++)
                          {
                             double expon = exp( - r[i_r]*r[i_r]/(2*sigma_gaussian*sigma_gaussian));
                             denominator = sqrt(x_local*x_local + y_local*y_local + r[i_r]*r[i_r] - 2*r[i_r]*(cos_phi[i_phi]*x_local+sin_phi[i_phi]*y_local) + z_local*z_local);
                             denominator_cubic = denominator*denominator*denominator;
                             Ex_integrand = norm_gaussian*nucleon_cosh_y[list_idx]*expon*(x_local - r[i_r]*cos_phi[i_phi])/denominator_cubic;
                             Ey_integrand = norm_gaussian*nucleon_cosh_y[list_idx]*expon*(y_local - r[i_r]*sin_phi[i_phi])/denominator_cubic;
                             Ez_integrand = norm_gaussian*expon*z_local/denominator_cubic;
                             Bx_integrand = norm_gaussian*expon*nucleon_sinh_y[list_idx]*(y_local - r[i_r]*sin_phi[i_phi])/denominator_cubic;
                             By_integrand = - norm_gaussian*expon*nucleon_sinh_y[list_idx]*(x_local - r[i_r]*cos_phi[i_phi])/denominator_cubic;
                             Bz_integrand = 0.0e0;
                             
                             Ex_integral += Ex_integrand*phi_weight[i_phi]*r_weight[i_r]*r[i_r];
                             Ey_integral += Ey_integrand*phi_weight[i_phi]*r_weight[i_r]*r[i_r];
                             Ez_integral += Ez_integrand*phi_weight[i_phi]*r_weight[i_r]*r[i_r];
                             Bx_integral += Bx_integrand*phi_weight[i_phi]*r_weight[i_r]*r[i_r];
                             By_integral += By_integrand*phi_weight[i_phi]*r_weight[i_r]*r[i_r];
                             Bz_integral += Bz_integrand*phi_weight[i_phi]*r_weight[i_r]*r[i_r];
                          }
                       break;
                     default:
                       cout << "EM_fields::calculate_EM_fields error: wrong type of nucleons, Nucleon_type = " << Nucleon_type << endl;
                       exit(1);
                       break;
                  }
                  temp_sum_Ex += e_sq*Ex_integral;
                  temp_sum_Ey += e_sq*Ey_integral;
                  temp_sum_Ez += e_sq*Ez_integral;
                  temp_sum_Bx += e_sq*Bx_integral;
                  temp_sum_By += e_sq*By_integral;
                  temp_sum_Bz += e_sq*Bz_integral;
                  Ex_integral = 0.0;
                  Ey_integral = 0.0;
                  Ez_integral = 0.0;
                  Bx_integral = 0.0;
                  By_integral = 0.0;
                  Bz_integral = 0.0;

               }
               E_x[i][j][k][l] = temp_sum_Ex*unit_convert;   // unit of [e*E_x] is MeV^2
               E_y[i][j][k][l] = temp_sum_Ey*unit_convert;
               E_z[i][j][k][l] = temp_sum_Ez*unit_convert;
               B_x[i][j][k][l] = temp_sum_Bx*unit_convert;
               B_y[i][j][k][l] = temp_sum_By*unit_convert;
               B_z[i][j][k][l] = temp_sum_Bz*unit_convert;
            }
   delete[] nucleon_position_x;
   delete[] nucleon_position_y;
   delete[] nucleon_rapidity;
   delete[] nucleon_sinh_y;
   delete[] nucleon_cosh_y;
   return;
}

void EM_fields::output_EM_fields_transverse_plane(int itime, int iz)
{
   ostringstream filename_string_Ex;
   ostringstream filename_string_Ey;
   ostringstream filename_string_Ez;
   ostringstream filename_string_Bx;
   ostringstream filename_string_By;
   ostringstream filename_string_Bz;
   filename_string_Ex << results_path << "Ex_fields_xy_t=" << setprecision(3) << t[itime] << "_z=" << z[iz] << ".dat";
   filename_string_Ey << results_path << "Ey_fields_xy_t=" << setprecision(3) << t[itime] << "_z=" << z[iz] << ".dat";
   filename_string_Ez << results_path << "Ez_fields_xy_t=" << setprecision(3) << t[itime] << "_z=" << z[iz] << ".dat";
   filename_string_Bx << results_path << "Bx_fields_xy_t=" << setprecision(3) << t[itime] << "_z=" << z[iz] << ".dat";
   filename_string_By << results_path << "By_fields_xy_t=" << setprecision(3) << t[itime] << "_z=" << z[iz] << ".dat";
   filename_string_Bz << results_path << "Bz_fields_xy_t=" << setprecision(3) << t[itime] << "_z=" << z[iz] << ".dat";
   ofstream output_Ex(filename_string_Ex.str().c_str());
   ofstream output_Ey(filename_string_Ey.str().c_str());
   ofstream output_Ez(filename_string_Ez.str().c_str());
   ofstream output_Bx(filename_string_Bx.str().c_str());
   ofstream output_By(filename_string_By.str().c_str());
   ofstream output_Bz(filename_string_Bz.str().c_str());
   for(int i=0; i<nx; i++)
   {
      for(int j=0; j<ny; j++)
      {
         output_Ex << scientific << setprecision(8) << setw(15)  
                   << E_x[itime][i][j][iz] << "   ";
         output_Ey << scientific << setprecision(8) << setw(15)  
                   << E_y[itime][i][j][iz] << "   ";
         output_Ez << scientific << setprecision(8) << setw(15)  
                   << E_z[itime][i][j][iz] << "   ";
         output_Bx << scientific << setprecision(8) << setw(15)  
                   << B_x[itime][i][j][iz] << "   ";
         output_By << scientific << setprecision(8) << setw(15)  
                   << B_y[itime][i][j][iz] << "   ";
         output_Bz << scientific << setprecision(8) << setw(15)  
                   << B_z[itime][i][j][iz] << "   ";
      }
      output_Ex << endl;
      output_Ey << endl;
      output_Ez << endl;
      output_Bx << endl;
      output_By << endl;
      output_Bz << endl;
   }
   output_Ex.close();
   output_Ey.close();
   output_Ez.close();
   output_Bx.close();
   output_By.close();
   output_Bz.close();
   return;
}

void EM_fields::output_EM_fields_time_evolution(int ix, int iy, int iz)
{
   ostringstream filename_string_Ex;
   ostringstream filename_string_Ey;
   ostringstream filename_string_Ez;
   ostringstream filename_string_Bx;
   ostringstream filename_string_By;
   ostringstream filename_string_Bz;
   filename_string_Ex << results_path << "Ex_fields_timeEvo_x=" << setprecision(3) << x[ix] << "_y=" << y[iy] << "_z=" << z[iz] << ".dat";
   filename_string_Ey << results_path << "Ey_fields_timeEvo_x=" << setprecision(3) << x[ix] << "_y=" << y[iy] << "_z=" << z[iz] << ".dat";
   filename_string_Ez << results_path << "Ez_fields_timeEvo_x=" << setprecision(3) << x[ix] << "_y=" << y[iy] << "_z=" << z[iz] << ".dat";
   filename_string_Bx << results_path << "Bx_fields_timeEvo_x=" << setprecision(3) << x[ix] << "_y=" << y[iy] << "_z=" << z[iz] << ".dat";
   filename_string_By << results_path << "By_fields_timeEvo_x=" << setprecision(3) << x[ix] << "_y=" << y[iy] << "_z=" << z[iz] << ".dat";
   filename_string_Bz << results_path << "Bz_fields_timeEvo_x=" << setprecision(3) << x[ix] << "_y=" << y[iy] << "_z=" << z[iz] << ".dat";
   ofstream output_Ex(filename_string_Ex.str().c_str());
   ofstream output_Ey(filename_string_Ey.str().c_str());
   ofstream output_Ez(filename_string_Ez.str().c_str());
   ofstream output_Bx(filename_string_Bx.str().c_str());
   ofstream output_By(filename_string_By.str().c_str());
   ofstream output_Bz(filename_string_Bz.str().c_str());
   for(int i=0; i<nt; i++)
   {
      output_Ex << scientific << setprecision(8) << setw(15)  
                << t[i] << "   " << E_x[i][ix][iy][iz] << endl;
      output_Ey << scientific << setprecision(8) << setw(15)  
                << t[i] << "   " << E_y[i][ix][iy][iz] << endl;
      output_Ez << scientific << setprecision(8) << setw(15)  
                << t[i] << "   " << E_z[i][ix][iy][iz] << endl;
      output_Bx << scientific << setprecision(8) << setw(15)  
                << t[i] << "   " << B_x[i][ix][iy][iz] << endl;
      output_By << scientific << setprecision(8) << setw(15)  
                << t[i] << "   " << B_y[i][ix][iy][iz] << endl;
      output_Bz << scientific << setprecision(8) << setw(15)  
                << t[i] << "   " << B_z[i][ix][iy][iz] << endl;
   }
   output_Ex.close();
   output_Ey.close();
   output_Ez.close();
   output_Bx.close();
   output_By.close();
   output_Bz.close();
   return;
}
