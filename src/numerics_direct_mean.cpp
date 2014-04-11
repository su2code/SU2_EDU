/*!
 * \file numerics_direct_mean.cpp
 * \brief This file contains all the convective term discretization.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 1.1.0
 *
 * SU2, Copyright (C) 2012-2014 Aerospace Design Laboratory (ADL).
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#include "../include/numerics_structure.hpp"
#include <limits>

CCentJST_Flow::CCentJST_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
  /*--- Artifical dissipation part ---*/
  Param_p = 0.3;
  Param_Kappa_2 = config->GetKappa_2nd_Flow();
  Param_Kappa_4 = config->GetKappa_4th_Flow();
  
  /*--- Allocate some structures ---*/
  Diff_U = new double [nVar];
  Diff_Lapl = new double [nVar];
  Velocity_i = new double [nDim];
  Velocity_j = new double [nDim];
  MeanVelocity = new double [nDim];
  Proj_flux_tensor = new double [nVar];
  
}

CCentJST_Flow::~CCentJST_Flow(void) {
  delete [] Diff_U;
  delete [] Diff_Lapl;
  delete [] Velocity_i;
  delete [] Velocity_j;
  delete [] MeanVelocity;
  delete [] Proj_flux_tensor;
}

void CCentJST_Flow::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j,
                                    CConfig *config) {
  
  /*--- Pressure, density, enthalpy, energy, and velocity at points i and j ---*/
  
  Pressure_i = V_i[nDim+1];                       Pressure_j = V_j[nDim+1];
  Density_i = V_i[nDim+2];                        Density_j = V_j[nDim+2];
  Enthalpy_i = V_i[nDim+3];                       Enthalpy_j = V_j[nDim+3];
  SoundSpeed_i = V_i[nDim+4];                     SoundSpeed_j = V_j[nDim+4];
  Energy_i = Enthalpy_i - Pressure_i/Density_i;   Energy_j = Enthalpy_j - Pressure_j/Density_j;
  
  sq_vel_i = 0.0; sq_vel_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim] = V_i[iDim+1];
    Velocity_j[iDim] = V_j[iDim+1];
    sq_vel_i += 0.5*Velocity_i[iDim]*Velocity_i[iDim];
    sq_vel_j += 0.5*Velocity_j[iDim]*Velocity_j[iDim];
  }
  
  /*--- Recompute conservative variables ---*/
  
  double U_i[5], U_j[5];
  U_i[0] = Density_i; U_j[0] = Density_j;
  for (iDim = 0; iDim < nDim; iDim++) {
    U_i[iDim+1] = Density_i*Velocity_i[iDim]; U_j[iDim+1] = Density_j*Velocity_j[iDim];
  }
  U_i[nDim+1] = Density_i*Energy_i; U_j[nDim+1] = Density_j*Energy_j;
  
  /*--- Compute mean values of the variables ---*/
  
  MeanDensity = 0.5*(Density_i+Density_j);
  MeanPressure = 0.5*(Pressure_i+Pressure_j);
  MeanEnthalpy = 0.5*(Enthalpy_i+Enthalpy_j);
  for (iDim = 0; iDim < nDim; iDim++)
    MeanVelocity[iDim] =  0.5*(Velocity_i[iDim]+Velocity_j[iDim]);
  MeanEnergy = 0.5*(Energy_i+Energy_j);
  
  /*--- Get projected flux tensor ---*/
  
  GetInviscidProjFlux(&MeanDensity, MeanVelocity, &MeanPressure, &MeanEnthalpy, Normal, Proj_flux_tensor);
  
  /*--- Residual of the inviscid flux ---*/
  
  for (iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] = Proj_flux_tensor[iVar];
  
  /*--- Jacobians of the inviscid flux, scale = 0.5 because val_residual ~ 0.5*(fc_i+fc_j)*Normal ---*/
  
  if (implicit) {
    GetInviscidProjJac(MeanVelocity, &MeanEnergy, Normal, 0.5, val_Jacobian_i);
    for (iVar = 0; iVar < nVar; iVar++)
      for (jVar = 0; jVar < nVar; jVar++)
        val_Jacobian_j[iVar][jVar] = val_Jacobian_i[iVar][jVar];
  }
  
  /*--- Computes differences btw. Laplacians and conservative variables,
   with a correction for the enthalpy ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    Diff_Lapl[iVar] = Und_Lapl_i[iVar]-Und_Lapl_j[iVar];
    Diff_U[iVar] = U_i[iVar]-U_j[iVar];
  }
  Diff_U[nVar-1] = Density_i*Enthalpy_i-Density_j*Enthalpy_j;
  
  /*--- Compute the local spectral radius and the stretching factor ---*/
  
  ProjVelocity_i = 0.0; ProjVelocity_j = 0.0; Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    ProjVelocity_i += Velocity_i[iDim]*Normal[iDim];
    ProjVelocity_j += Velocity_j[iDim]*Normal[iDim];
    Area += Normal[iDim]*Normal[iDim];
  }
  Area = sqrt(Area);
  
  Local_Lambda_i = (fabs(ProjVelocity_i)+SoundSpeed_i*Area);
  Local_Lambda_j = (fabs(ProjVelocity_j)+SoundSpeed_j*Area);
  MeanLambda = 0.5*(Local_Lambda_i+Local_Lambda_j);
  
  Phi_i = pow(Lambda_i/(4.0*MeanLambda), Param_p);
  Phi_j = pow(Lambda_j/(4.0*MeanLambda), Param_p);
  StretchingFactor = 4.0*Phi_i*Phi_j/(Phi_i+Phi_j);
  
  sc2 = 3.0*(double(Neighbor_i)+double(Neighbor_j))/(double(Neighbor_i)*double(Neighbor_j));
  sc4 = sc2*sc2/4.0;
  
  Epsilon_2 = Param_Kappa_2*0.5*(Sensor_i+Sensor_j)*sc2;
  Epsilon_4 = max(0.0, Param_Kappa_4-Epsilon_2)*sc4;
  
  /*--- Compute viscous part of the residual ---*/
  
  for (iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] += (Epsilon_2*Diff_U[iVar] - Epsilon_4*Diff_Lapl[iVar])*StretchingFactor*MeanLambda;
  
  /*--- Jacobian computation ---*/
  
  if (implicit) {
    
    cte_0 = (Epsilon_2 + Epsilon_4*double(Neighbor_i+1))*StretchingFactor*MeanLambda;
    cte_1 = (Epsilon_2 + Epsilon_4*double(Neighbor_j+1))*StretchingFactor*MeanLambda;
    
    for (iVar = 0; iVar < (nVar-1); iVar++) {
      val_Jacobian_i[iVar][iVar] += cte_0;
      val_Jacobian_j[iVar][iVar] -= cte_1;
    }
    
    /*--- Last row of Jacobian_i ---*/
    
    val_Jacobian_i[nVar-1][0] += cte_0*Gamma_Minus_One*sq_vel_i;
    for (iDim = 0; iDim < nDim; iDim++)
      val_Jacobian_i[nVar-1][iDim+1] -= cte_0*Gamma_Minus_One*Velocity_i[iDim];
    val_Jacobian_i[nVar-1][nVar-1] += cte_0*Gamma;
    
    /*--- Last row of Jacobian_j ---*/
    
    val_Jacobian_j[nVar-1][0] -= cte_1*Gamma_Minus_One*sq_vel_j;
    for (iDim = 0; iDim < nDim; iDim++)
      val_Jacobian_j[nVar-1][iDim+1] += cte_1*Gamma_Minus_One*Velocity_j[iDim];
    val_Jacobian_j[nVar-1][nVar-1] -= cte_1*Gamma;
    
  }
  
}

void CCentJST_Flow::GetInviscidProjFlux(double *val_density,
                                    double *val_velocity,
                                    double *val_pressure,
                                    double *val_enthalpy,
                                    double *val_normal,
                                    double *val_Proj_Flux) {
  double rhou, rhov, rhow;
  
  
  if (nDim == 2) {
    rhou = (*val_density)*val_velocity[0];
    rhov = (*val_density)*val_velocity[1];
    
    val_Proj_Flux[0] = rhou*val_normal[0];
    val_Proj_Flux[1] = (rhou*val_velocity[0]+(*val_pressure))*val_normal[0];
    val_Proj_Flux[2] = rhou*val_velocity[1]*val_normal[0];
    val_Proj_Flux[3] = rhou*(*val_enthalpy)*val_normal[0];
    
    val_Proj_Flux[0] += rhov*val_normal[1];
    val_Proj_Flux[1] += rhov*val_velocity[0]*val_normal[1];
    val_Proj_Flux[2] += (rhov*val_velocity[1]+(*val_pressure))*val_normal[1];
    val_Proj_Flux[3] += rhov*(*val_enthalpy)*val_normal[1];
  }
  else {
    rhou = (*val_density)*val_velocity[0];
    rhov = (*val_density)*val_velocity[1];
    rhow = (*val_density)*val_velocity[2];
    
    val_Proj_Flux[0] = rhou*val_normal[0];
    val_Proj_Flux[1] = (rhou*val_velocity[0]+(*val_pressure))*val_normal[0];
    val_Proj_Flux[2] = rhou*val_velocity[1]*val_normal[0];
    val_Proj_Flux[3] = rhou*val_velocity[2]*val_normal[0];
    val_Proj_Flux[4] = rhou*(*val_enthalpy)*val_normal[0];
    
    val_Proj_Flux[0] += rhov*val_normal[1];
    val_Proj_Flux[1] += rhov*val_velocity[0]*val_normal[1];
    val_Proj_Flux[2] += (rhov*val_velocity[1]+(*val_pressure))*val_normal[1];
    val_Proj_Flux[3] += rhov*val_velocity[2]*val_normal[1];
    val_Proj_Flux[4] += rhov*(*val_enthalpy)*val_normal[1];
    
    val_Proj_Flux[0] += rhow*val_normal[2];
    val_Proj_Flux[1] += rhow*val_velocity[0]*val_normal[2];
    val_Proj_Flux[2] += rhow*val_velocity[1]*val_normal[2];
    val_Proj_Flux[3] += (rhow*val_velocity[2]+(*val_pressure))*val_normal[2];
    val_Proj_Flux[4] += rhow*(*val_enthalpy)*val_normal[2];
  }
  
}

CCentLax_Flow::CCentLax_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  
  /*--- Artifical dissipation part ---*/
  Param_p = 0.3;
  Param_Kappa_0 = config->GetKappa_1st_Flow();
  
  /*--- Allocate some structures ---*/
  Diff_U = new double [nVar];
  Velocity_i = new double [nDim];
  Velocity_j = new double [nDim];
  MeanVelocity = new double [nDim];
  Proj_flux_tensor = new double [nVar];
  
}

CCentLax_Flow::~CCentLax_Flow(void) {
  delete [] Diff_U;
  delete [] Velocity_i;
  delete [] Velocity_j;
  delete [] MeanVelocity;
  delete [] Proj_flux_tensor;
  
}

void CCentLax_Flow::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j,
                                    CConfig *config) {
  
  /*--- Pressure, density, enthalpy, energy, and velocity at points i and j ---*/
  
  Pressure_i = V_i[nDim+1];                       Pressure_j = V_j[nDim+1];
  Density_i = V_i[nDim+2];                        Density_j = V_j[nDim+2];
  Enthalpy_i = V_i[nDim+3];                       Enthalpy_j = V_j[nDim+3];
  SoundSpeed_i = V_i[nDim+4];                     SoundSpeed_j = V_j[nDim+4];
  Energy_i = Enthalpy_i - Pressure_i/Density_i;   Energy_j = Enthalpy_j - Pressure_j/Density_j;
  
  sq_vel_i = 0.0; sq_vel_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim] = V_i[iDim+1];
    Velocity_j[iDim] = V_j[iDim+1];
    sq_vel_i += 0.5*Velocity_i[iDim]*Velocity_i[iDim];
    sq_vel_j += 0.5*Velocity_j[iDim]*Velocity_j[iDim];
  }
  
  /*--- Recompute conservative variables ---*/
  
  double U_i[5], U_j[5];
  U_i[0] = Density_i; U_j[0] = Density_j;
  for (iDim = 0; iDim < nDim; iDim++) {
    U_i[iDim+1] = Density_i*Velocity_i[iDim]; U_j[iDim+1] = Density_j*Velocity_j[iDim];
  }
  U_i[nDim+1] = Density_i*Energy_i; U_j[nDim+1] = Density_j*Energy_j;
  
  /*--- Compute mean values of the variables ---*/
  
  MeanDensity = 0.5*(Density_i+Density_j);
  MeanPressure = 0.5*(Pressure_i+Pressure_j);
  MeanEnthalpy = 0.5*(Enthalpy_i+Enthalpy_j);
  for (iDim = 0; iDim < nDim; iDim++)
    MeanVelocity[iDim] =  0.5*(Velocity_i[iDim]+Velocity_j[iDim]);
  MeanEnergy = 0.5*(Energy_i+Energy_j);
  
  /*--- Get projected flux tensor ---*/
  
  GetInviscidProjFlux(&MeanDensity, MeanVelocity, &MeanPressure, &MeanEnthalpy, Normal, Proj_flux_tensor);
  
  /*--- Residual of the inviscid flux ---*/
  
  for (iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] = Proj_flux_tensor[iVar];
  
  /*--- Jacobians of the inviscid flux, scale = 0.5 because val_residual ~ 0.5*(fc_i+fc_j)*Normal ---*/
  
  if (implicit) {
    GetInviscidProjJac(MeanVelocity, &MeanEnergy, Normal, 0.5, val_Jacobian_i);
    for (iVar = 0; iVar < nVar; iVar++)
      for (jVar = 0; jVar < nVar; jVar++)
        val_Jacobian_j[iVar][jVar] = val_Jacobian_i[iVar][jVar];
  }
  
  /*--- Computes differences btw. conservative variables,
   with a correction for the enthalpy ---*/
  
  for (iVar = 0; iVar < nDim+1; iVar++)
    Diff_U[iVar] = U_i[iVar]-U_j[iVar];
  Diff_U[nDim+1] = Density_i*Enthalpy_i-Density_j*Enthalpy_j;
  
  /*--- Compute the local spectral radius and the stretching factor ---*/
  
  ProjVelocity_i = 0.0; ProjVelocity_j = 0.0; Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    ProjVelocity_i += Velocity_i[iDim]*Normal[iDim];
    ProjVelocity_j += Velocity_j[iDim]*Normal[iDim];
    Area += Normal[iDim]*Normal[iDim];
  }
  Area = sqrt(Area);
  
  Local_Lambda_i = (fabs(ProjVelocity_i)+SoundSpeed_i*Area);
  Local_Lambda_j = (fabs(ProjVelocity_j)+SoundSpeed_j*Area);
  MeanLambda = 0.5*(Local_Lambda_i+Local_Lambda_j);
  
  Phi_i = pow(Lambda_i/(4.0*MeanLambda),Param_p);
  Phi_j = pow(Lambda_j/(4.0*MeanLambda),Param_p);
  StretchingFactor = 4.0*Phi_i*Phi_j/(Phi_i+Phi_j);
  
  sc0 = 3.0*(double(Neighbor_i)+double(Neighbor_j))/(double(Neighbor_i)*double(Neighbor_j));
  Epsilon_0 = Param_Kappa_0*sc0*double(nDim)/3.0;
  
  /*--- Compute viscous part of the residual ---*/
  
  for (iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] += Epsilon_0*Diff_U[iVar]*StretchingFactor*MeanLambda;
  
  /*--- Jacobian computation ---*/
  
  if (implicit) {
    cte = Epsilon_0*StretchingFactor*MeanLambda;
    for (iVar = 0; iVar < (nVar-1); iVar++) {
      val_Jacobian_i[iVar][iVar] += cte;
      val_Jacobian_j[iVar][iVar] -= cte;
    }
    
    /*--- Last row of Jacobian_i ---*/
    
    val_Jacobian_i[nVar-1][0] += cte*Gamma_Minus_One*sq_vel_i;
    for (iDim = 0; iDim < nDim; iDim++)
      val_Jacobian_i[nVar-1][iDim+1] -= cte*Gamma_Minus_One*Velocity_i[iDim];
    val_Jacobian_i[nVar-1][nVar-1] += cte*Gamma;
    
    /*--- Last row of Jacobian_j ---*/
    
    val_Jacobian_j[nVar-1][0] -= cte*Gamma_Minus_One*sq_vel_j;
    for (iDim = 0; iDim < nDim; iDim++)
      val_Jacobian_j[nVar-1][iDim+1] += cte*Gamma_Minus_One*Velocity_j[iDim];
    val_Jacobian_j[nVar-1][nVar-1] -= cte*Gamma;
    
  }
  
}

void CCentLax_Flow::GetInviscidProjFlux(double *val_density,
                                        double *val_velocity,
                                        double *val_pressure,
                                        double *val_enthalpy,
                                        double *val_normal,
                                        double *val_Proj_Flux) {
  double rhou, rhov, rhow;
  
  
  if (nDim == 2) {
    rhou = (*val_density)*val_velocity[0];
    rhov = (*val_density)*val_velocity[1];
    
    val_Proj_Flux[0] = rhou*val_normal[0];
    val_Proj_Flux[1] = (rhou*val_velocity[0]+(*val_pressure))*val_normal[0];
    val_Proj_Flux[2] = rhou*val_velocity[1]*val_normal[0];
    val_Proj_Flux[3] = rhou*(*val_enthalpy)*val_normal[0];
    
    val_Proj_Flux[0] += rhov*val_normal[1];
    val_Proj_Flux[1] += rhov*val_velocity[0]*val_normal[1];
    val_Proj_Flux[2] += (rhov*val_velocity[1]+(*val_pressure))*val_normal[1];
    val_Proj_Flux[3] += rhov*(*val_enthalpy)*val_normal[1];
  }
  else {
    rhou = (*val_density)*val_velocity[0];
    rhov = (*val_density)*val_velocity[1];
    rhow = (*val_density)*val_velocity[2];
    
    val_Proj_Flux[0] = rhou*val_normal[0];
    val_Proj_Flux[1] = (rhou*val_velocity[0]+(*val_pressure))*val_normal[0];
    val_Proj_Flux[2] = rhou*val_velocity[1]*val_normal[0];
    val_Proj_Flux[3] = rhou*val_velocity[2]*val_normal[0];
    val_Proj_Flux[4] = rhou*(*val_enthalpy)*val_normal[0];
    
    val_Proj_Flux[0] += rhov*val_normal[1];
    val_Proj_Flux[1] += rhov*val_velocity[0]*val_normal[1];
    val_Proj_Flux[2] += (rhov*val_velocity[1]+(*val_pressure))*val_normal[1];
    val_Proj_Flux[3] += rhov*val_velocity[2]*val_normal[1];
    val_Proj_Flux[4] += rhov*(*val_enthalpy)*val_normal[1];
    
    val_Proj_Flux[0] += rhow*val_normal[2];
    val_Proj_Flux[1] += rhow*val_velocity[0]*val_normal[2];
    val_Proj_Flux[2] += rhow*val_velocity[1]*val_normal[2];
    val_Proj_Flux[3] += (rhow*val_velocity[2]+(*val_pressure))*val_normal[2];
    val_Proj_Flux[4] += rhow*(*val_enthalpy)*val_normal[2];
  }
  
}

CUpwAUSM_Flow::CUpwAUSM_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
  Diff_U = new double [nVar];
  Velocity_i = new double [nDim];
  Velocity_j = new double [nDim];
  RoeVelocity = new double [nDim];
  delta_vel  = new double [nDim];
  delta_wave = new double [nVar];
  Proj_flux_tensor_i = new double [nVar];
  Proj_flux_tensor_j = new double [nVar];
  Lambda = new double [nVar];
  Epsilon = new double [nVar];
  P_Tensor = new double* [nVar];
  invP_Tensor = new double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    P_Tensor[iVar] = new double [nVar];
    invP_Tensor[iVar] = new double [nVar];
  }
}

CUpwAUSM_Flow::~CUpwAUSM_Flow(void) {
  
  delete [] Diff_U;
  delete [] Velocity_i;
  delete [] Velocity_j;
  delete [] RoeVelocity;
  delete [] delta_vel;
  delete [] delta_wave;
  delete [] Proj_flux_tensor_i;
  delete [] Proj_flux_tensor_j;
  delete [] Lambda;
  delete [] Epsilon;
  for (iVar = 0; iVar < nVar; iVar++) {
    delete [] P_Tensor[iVar];
    delete [] invP_Tensor[iVar];
  }
  delete [] P_Tensor;
  delete [] invP_Tensor;
  
}

void CUpwAUSM_Flow::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {
  
  /*--- Face area (norm or the normal vector) ---*/
  Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    Area += Normal[iDim]*Normal[iDim];
  Area = sqrt(Area);
  
  /*-- Unit Normal ---*/
  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;
  
  /*--- Primitive variables at point i ---*/
  sq_vel = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim] = V_i[iDim+1];
    sq_vel += Velocity_i[iDim]*Velocity_i[iDim];
  }
  Pressure_i = V_i[nDim+1];
  Density_i = V_i[nDim+2];
  Enthalpy_i = V_i[nDim+3];
  Energy_i = Enthalpy_i - Pressure_i/Density_i;
  SoundSpeed_i = sqrt(Gamma*Gamma_Minus_One*(Energy_i-0.5*sq_vel));
  
  /*--- Primitive variables at point j ---*/
  sq_vel = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_j[iDim] = V_j[iDim+1];
    sq_vel += Velocity_j[iDim]*Velocity_j[iDim];
  }
  Pressure_j = V_j[nDim+1];
  Density_j = V_j[nDim+2];
  Enthalpy_j = V_j[nDim+3];
  Energy_j = Enthalpy_j - Pressure_j/Density_j;
  SoundSpeed_j = sqrt(Gamma*Gamma_Minus_One*(Energy_j-0.5*sq_vel));
  
  /*--- Projected velocities ---*/
  ProjVelocity_i = 0.0; ProjVelocity_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    ProjVelocity_i += Velocity_i[iDim]*UnitNormal[iDim];
    ProjVelocity_j += Velocity_j[iDim]*UnitNormal[iDim];
  }
  
  mL	= ProjVelocity_i/SoundSpeed_i;
  mR	= ProjVelocity_j/SoundSpeed_j;
  
  if (fabs(mL) <= 1.0) mLP = 0.25*(mL+1.0)*(mL+1.0);
  else mLP = 0.5*(mL+fabs(mL));
  
  if (fabs(mR) <= 1.0) mRM = -0.25*(mR-1.0)*(mR-1.0);
  else mRM = 0.5*(mR-fabs(mR));
  
  mF = mLP + mRM;
  
  if (fabs(mL) <= 1.0) pLP = 0.25*Pressure_i*(mL+1.0)*(mL+1.0)*(2.0-mL);
  else pLP = 0.5*Pressure_i*(mL+fabs(mL))/mL;
  
  if (fabs(mR) <= 1.0) pRM = 0.25*Pressure_j*(mR-1.0)*(mR-1.0)*(2.0+mR);
  else pRM = 0.5*Pressure_j*(mR-fabs(mR))/mR;
  
  pF = pLP + pRM;
  Phi = fabs(mF);
  
  val_residual[0] = 0.5*(mF*((Density_i*SoundSpeed_i)+(Density_j*SoundSpeed_j))-Phi*((Density_j*SoundSpeed_j)-(Density_i*SoundSpeed_i)));
  for (iDim = 0; iDim < nDim; iDim++)
    val_residual[iDim+1] = 0.5*(mF*((Density_i*SoundSpeed_i*Velocity_i[iDim])+(Density_j*SoundSpeed_j*Velocity_j[iDim]))
                                -Phi*((Density_j*SoundSpeed_j*Velocity_j[iDim])-(Density_i*SoundSpeed_i*Velocity_i[iDim])))+UnitNormal[iDim]*pF;
  val_residual[nVar-1] = 0.5*(mF*((Density_i*SoundSpeed_i*Enthalpy_i)+(Density_j*SoundSpeed_j*Enthalpy_j))-Phi*((Density_j*SoundSpeed_j*Enthalpy_j)-(Density_i*SoundSpeed_i*Enthalpy_i)));
  
  for (iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] *= Area;
  
  /*--- Roe's Jacobian for AUSM (this must be fixed) ---*/
  if (implicit) {
    
    /*--- Promediate Roe variables iPoint and jPoint ---*/
    R = sqrt(Density_j/Density_i);
    RoeDensity = R*Density_i;
    sq_vel = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      RoeVelocity[iDim] = (R*Velocity_j[iDim]+Velocity_i[iDim])/(R+1);
      sq_vel += RoeVelocity[iDim]*RoeVelocity[iDim];
    }
    RoeEnthalpy = (R*Enthalpy_j+Enthalpy_i)/(R+1);
    RoeSoundSpeed = sqrt((Gamma-1)*(RoeEnthalpy-0.5*sq_vel));
    
    /*--- Compute P and Lambda (do it with the Normal) ---*/
    GetPMatrix(&RoeDensity, RoeVelocity, &RoeSoundSpeed, UnitNormal, P_Tensor);
    
    ProjVelocity = 0.0; ProjVelocity_i = 0.0; ProjVelocity_j = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      ProjVelocity   += RoeVelocity[iDim]*UnitNormal[iDim];
      ProjVelocity_i += Velocity_i[iDim]*UnitNormal[iDim];
      ProjVelocity_j += Velocity_j[iDim]*UnitNormal[iDim];
    }
    
    /*--- Flow eigenvalues and Entropy correctors ---*/
    for (iDim = 0; iDim < nDim; iDim++)
      Lambda[iDim] = ProjVelocity;
    Lambda[nVar-2]  = ProjVelocity + RoeSoundSpeed;
    Lambda[nVar-1] = ProjVelocity - RoeSoundSpeed;
    
    /*--- Compute inverse P ---*/
    GetPMatrix_inv(&RoeDensity, RoeVelocity, &RoeSoundSpeed, UnitNormal, invP_Tensor);
    
    /*--- Jacobias of the inviscid flux, scale = 0.5 because val_residual ~ 0.5*(fc_i+fc_j)*Normal ---*/
    GetInviscidProjJac(Velocity_i, &Energy_i, Normal, 0.5, val_Jacobian_i);
    GetInviscidProjJac(Velocity_j, &Energy_j, Normal, 0.5, val_Jacobian_j);
    
    /*--- Roe's Flux approximation ---*/
    for (iVar = 0; iVar < nVar; iVar++) {
      for (jVar = 0; jVar < nVar; jVar++) {
        Proj_ModJac_Tensor_ij = 0.0;
        /*--- Compute |Proj_ModJac_Tensor| = P x |Lambda| x inverse P ---*/
        for (kVar = 0; kVar < nVar; kVar++)
          Proj_ModJac_Tensor_ij += P_Tensor[iVar][kVar]*fabs(Lambda[kVar])*invP_Tensor[kVar][jVar];
        val_Jacobian_i[iVar][jVar] += 0.5*Proj_ModJac_Tensor_ij*Area;
        val_Jacobian_j[iVar][jVar] -= 0.5*Proj_ModJac_Tensor_ij*Area;
      }
    }
  }
}

CUpwHLLC_Flow::CUpwHLLC_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
  Diff_U = new double [nVar];
  Velocity_i = new double [nDim];
  Velocity_j = new double [nDim];
  RoeVelocity = new double [nDim];
  delta_vel  = new double [nDim];
  delta_wave = new double [nVar];
  Proj_flux_tensor_i = new double [nVar];
  Proj_flux_tensor_j = new double [nVar];
  Lambda = new double [nVar];
  Epsilon = new double [nVar];
  P_Tensor = new double* [nVar];
  invP_Tensor = new double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    P_Tensor[iVar] = new double [nVar];
    invP_Tensor[iVar] = new double [nVar];
  }
}

CUpwHLLC_Flow::~CUpwHLLC_Flow(void) {
  
  delete [] Diff_U;
  delete [] Velocity_i;
  delete [] Velocity_j;
  delete [] RoeVelocity;
  delete [] delta_vel;
  delete [] delta_wave;
  delete [] Proj_flux_tensor_i;
  delete [] Proj_flux_tensor_j;
  delete [] Lambda;
  delete [] Epsilon;
  for (iVar = 0; iVar < nVar; iVar++) {
    delete [] P_Tensor[iVar];
    delete [] invP_Tensor[iVar];
  }
  delete [] P_Tensor;
  delete [] invP_Tensor;
  
}

void CUpwHLLC_Flow::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {
  
  /*--- Face area (norm or the normal vector) ---*/
  Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    Area += Normal[iDim]*Normal[iDim];
  Area = sqrt(Area);
  
  /*-- Unit Normal ---*/
  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;
  
  /*--- Primitive variables at point i ---*/
  sq_vel_i = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim] = V_i[iDim+1];
    sq_vel_i += Velocity_i[iDim]*Velocity_i[iDim];
  }
  Pressure_i = V_i[nDim+1];
  Density_i = V_i[nDim+2];
  Enthalpy_i = V_i[nDim+3];
  Energy_i = Enthalpy_i - Pressure_i/Density_i;
  SoundSpeed_i = sqrt(Gamma*Gamma_Minus_One*(Energy_i-0.5*sq_vel_i));
  
  /*--- Primitive variables at point j ---*/
  sq_vel_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_j[iDim] = V_j[iDim+1];
    sq_vel_j += Velocity_j[iDim]*Velocity_j[iDim];
  }
  Pressure_j = V_j[nDim+1];
  Density_j = V_j[nDim+2];
  Enthalpy_j = V_j[nDim+3];
  Energy_j = Enthalpy_j - Pressure_j/Density_j;
  SoundSpeed_j = sqrt(Gamma*Gamma_Minus_One*(Energy_j-0.5*sq_vel_j));
  
  /*--- Projected velocities ---*/
  ProjVelocity_i = 0.0; ProjVelocity_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    ProjVelocity_i += Velocity_i[iDim]*UnitNormal[iDim];
    ProjVelocity_j += Velocity_j[iDim]*UnitNormal[iDim];
  }
  
  /*--- Roe's averaging ---*/
  Rrho = sqrt(Density_j/Density_i);
  tmp = 1.0/(1.0+Rrho);
  for (iDim = 0; iDim < nDim; iDim++)
    velRoe[iDim] = tmp*(Velocity_i[iDim] + Velocity_j[iDim]*Rrho);
  
  uRoe  = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    uRoe += velRoe[iDim]*UnitNormal[iDim];
  
  gamPdivRho = tmp*((Gamma*Pressure_i/Density_i+0.5*(Gamma-1.0)*sq_vel_i) + (Gamma*Pressure_j/Density_j+0.5*(Gamma-1.0)*sq_vel_j)*Rrho);
  sq_velRoe = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    sq_velRoe += velRoe[iDim]*velRoe[iDim];
  
  cRoe  = sqrt(gamPdivRho - ((Gamma+Gamma)*0.5-1.0)*0.5*sq_velRoe);
  
  /*--- Speed of sound at L and R ---*/
  sL = min(uRoe-cRoe, ProjVelocity_i-SoundSpeed_i);
  sR = max(uRoe+cRoe, ProjVelocity_j+SoundSpeed_j);
  
  /*--- speed of contact surface ---*/
  sM = (Pressure_i-Pressure_j
        - Density_i*ProjVelocity_i*(sL-ProjVelocity_i)
        + Density_j*ProjVelocity_j*(sR-ProjVelocity_j))
  /(Density_j*(sR-ProjVelocity_j)-Density_i*(sL-ProjVelocity_i));
  
  /*--- Pressure at right and left (Pressure_j=Pressure_i) side of contact surface ---*/
  pStar = Density_j * (ProjVelocity_j-sR)*(ProjVelocity_j-sM) + Pressure_j;
  
  if (sM >= 0.0) {
    if (sL > 0.0) {
      val_residual[0] = Density_i*ProjVelocity_i;
      for (iDim = 0; iDim < nDim; iDim++)
        val_residual[iDim+1] = Density_i*Velocity_i[iDim]*ProjVelocity_i + Pressure_i*UnitNormal[iDim];
      val_residual[nVar-1] = Energy_i*Density_i*ProjVelocity_i + Pressure_i*ProjVelocity_i;
    }
    else {
      invSLmSs = 1.0/(sL-sM);
      sLmuL = sL-ProjVelocity_i;
      rhoSL = Density_i*sLmuL*invSLmSs;
      for (iDim = 0; iDim < nDim; iDim++)
        rhouSL[iDim] = (Density_i*Velocity_i[iDim]*sLmuL+(pStar-Pressure_i)*UnitNormal[iDim])*invSLmSs;
      eSL = (sLmuL*Energy_i*Density_i-Pressure_i*ProjVelocity_i+pStar*sM)*invSLmSs;
      
      val_residual[0] = rhoSL*sM;
      for (iDim = 0; iDim < nDim; iDim++)
        val_residual[iDim+1] = rhouSL[iDim]*sM + pStar*UnitNormal[iDim];
      val_residual[nVar-1] = (eSL+pStar)*sM;
    }
  }
  else {
    if (sR >= 0.0) {
      invSRmSs = 1.0/(sR-sM);
      sRmuR = sR-ProjVelocity_j;
      rhoSR = Density_j*sRmuR*invSRmSs;
      for (iDim = 0; iDim < nDim; iDim++)
        rhouSR[iDim] = (Density_j*Velocity_j[iDim]*sRmuR+(pStar-Pressure_j)*UnitNormal[iDim])*invSRmSs;
      eSR = (sRmuR*Energy_j*Density_j-Pressure_j*ProjVelocity_j+pStar*sM)*invSRmSs;
      
      val_residual[0] = rhoSR*sM;
      for (iDim = 0; iDim < nDim; iDim++)
        val_residual[iDim+1] = rhouSR[iDim]*sM + pStar*UnitNormal[iDim];
      val_residual[nVar-1] = (eSR+pStar)*sM;
    }
    else {
      val_residual[0] = Density_j*ProjVelocity_j;
      for (iDim = 0; iDim < nDim; iDim++)
        val_residual[iDim+1] = Density_j*Velocity_j[iDim]*ProjVelocity_j + Pressure_j*UnitNormal[iDim];
      val_residual[nVar-1] = Energy_j*Density_j*ProjVelocity_j + Pressure_j*ProjVelocity_j;
    }
  }
  
  for (iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] *= Area;
  
  if (implicit) {
    
    /*--- Average Roe variables iPoint and jPoint ---*/
    R = sqrt(Density_j/Density_i);
    RoeDensity = R*Density_i;
    sq_vel = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      RoeVelocity[iDim] = (R*Velocity_j[iDim]+Velocity_i[iDim])/(R+1);
      sq_vel += RoeVelocity[iDim]*RoeVelocity[iDim];
    }
    RoeEnthalpy = (R*Enthalpy_j+Enthalpy_i)/(R+1);
    RoeSoundSpeed = sqrt((Gamma-1)*(RoeEnthalpy-0.5*sq_vel));
    
    /*--- Compute P and Lambda (do it with the Normal) ---*/
    GetPMatrix(&RoeDensity, RoeVelocity, &RoeSoundSpeed, UnitNormal, P_Tensor);
    
    ProjVelocity = 0.0; ProjVelocity_i = 0.0; ProjVelocity_j = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      ProjVelocity   += RoeVelocity[iDim]*UnitNormal[iDim];
      ProjVelocity_i += Velocity_i[iDim]*UnitNormal[iDim];
      ProjVelocity_j += Velocity_j[iDim]*UnitNormal[iDim];
    }
    
    /*--- Flow eigenvalues and Entropy correctors ---*/
    for (iDim = 0; iDim < nDim; iDim++)
      Lambda[iDim] = ProjVelocity;
    Lambda[nVar-2]  = ProjVelocity + RoeSoundSpeed;
    Lambda[nVar-1] = ProjVelocity - RoeSoundSpeed;
    
    /*--- Compute inverse P ---*/
    GetPMatrix_inv(&RoeDensity, RoeVelocity, &RoeSoundSpeed, UnitNormal, invP_Tensor);
    
    /*--- Jacobias of the inviscid flux, scale = 0.5 because val_residual ~ 0.5*(fc_i+fc_j)*Normal ---*/
    GetInviscidProjJac(Velocity_i, &Energy_i, Normal, 0.5, val_Jacobian_i);
    GetInviscidProjJac(Velocity_j, &Energy_j, Normal, 0.5, val_Jacobian_j);
    
    /*--- Roe's Flux approximation ---*/
    for (iVar = 0; iVar < nVar; iVar++) {
      for (jVar = 0; jVar < nVar; jVar++) {
        Proj_ModJac_Tensor_ij = 0.0;
        /*--- Compute |Proj_ModJac_Tensor| = P x |Lambda| x inverse P ---*/
        for (kVar = 0; kVar < nVar; kVar++)
          Proj_ModJac_Tensor_ij += P_Tensor[iVar][kVar]*fabs(Lambda[kVar])*invP_Tensor[kVar][jVar];
        val_Jacobian_i[iVar][jVar] += 0.5*Proj_ModJac_Tensor_ij*Area;
        val_Jacobian_j[iVar][jVar] -= 0.5*Proj_ModJac_Tensor_ij*Area;
      }
    }
  }
  
}

CUpwRoe_Flow::CUpwRoe_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
  Diff_U = new double [nVar];
  Velocity_i = new double [nDim];
  Velocity_j = new double [nDim];
  RoeVelocity = new double [nDim];
  delta_vel  = new double [nDim];
  delta_wave = new double [nVar];
  Proj_flux_tensor_i = new double [nVar];
  Proj_flux_tensor_j = new double [nVar];
  Lambda = new double [nVar];
  Epsilon = new double [nVar];
  P_Tensor = new double* [nVar];
  invP_Tensor = new double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    P_Tensor[iVar] = new double [nVar];
    invP_Tensor[iVar] = new double [nVar];
  }
}

CUpwRoe_Flow::~CUpwRoe_Flow(void) {
  
  delete [] Diff_U;
  delete [] Velocity_i;
  delete [] Velocity_j;
  delete [] RoeVelocity;
  delete [] delta_vel;
  delete [] delta_wave;
  delete [] Proj_flux_tensor_i;
  delete [] Proj_flux_tensor_j;
  delete [] Lambda;
  delete [] Epsilon;
  for (iVar = 0; iVar < nVar; iVar++) {
    delete [] P_Tensor[iVar];
    delete [] invP_Tensor[iVar];
  }
  delete [] P_Tensor;
  delete [] invP_Tensor;
  
}

void CUpwRoe_Flow::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {
  
  /*--- Face area (norm or the normal vector) ---*/
  
  Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    Area += Normal[iDim]*Normal[iDim];
  Area = sqrt(Area);
  
  /*-- Unit Normal ---*/
  
  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;
  
  /*--- Primitive variables at point i ---*/
  
  for (iDim = 0; iDim < nDim; iDim++)
    Velocity_i[iDim] = V_i[iDim+1];
  Pressure_i = V_i[nDim+1];
  Density_i = V_i[nDim+2];
  Enthalpy_i = V_i[nDim+3];
  Energy_i = Enthalpy_i - Pressure_i/Density_i;
  SoundSpeed_i = sqrt(Pressure_i*Gamma/Density_i);
  
  /*--- Primitive variables at point j ---*/
  
  for (iDim = 0; iDim < nDim; iDim++)
    Velocity_j[iDim] = V_j[iDim+1];
  Pressure_j = V_j[nDim+1];
  Density_j = V_j[nDim+2];
  Enthalpy_j = V_j[nDim+3];
  Energy_j = Enthalpy_j - Pressure_j/Density_j;
  SoundSpeed_j = sqrt(Pressure_j*Gamma/Density_j);
  
  /*--- Recompute conservative variables ---*/
  
  double U_i[5], U_j[5];
  U_i[0] = Density_i; U_j[0] = Density_j;
  for (iDim = 0; iDim < nDim; iDim++) {
    U_i[iDim+1] = Density_i*Velocity_i[iDim]; U_j[iDim+1] = Density_j*Velocity_j[iDim];
  }
  U_i[nDim+1] = Density_i*Energy_i; U_j[nDim+1] = Density_j*Energy_j;
  
  /*--- Roe-averaged variables at interface between i & j ---*/
  
  R = sqrt(fabs(Density_j/Density_i));
  RoeDensity = R*Density_i;
  sq_vel = 0;
  for (iDim = 0; iDim < nDim; iDim++) {
    RoeVelocity[iDim] = (R*Velocity_j[iDim]+Velocity_i[iDim])/(R+1);
    sq_vel += RoeVelocity[iDim]*RoeVelocity[iDim];
  }
  RoeEnthalpy = (R*Enthalpy_j+Enthalpy_i)/(R+1);
  RoeSoundSpeed = sqrt((Gamma-1)*(RoeEnthalpy-0.5*sq_vel));
  
  /*--- Compute Proj_flux_tensor_i ---*/
  GetInviscidProjFlux(&Density_i, Velocity_i, &Pressure_i, &Enthalpy_i, Normal, Proj_flux_tensor_i);
  
  /*--- Compute Proj_flux_tensor_j ---*/
  GetInviscidProjFlux(&Density_j, Velocity_j, &Pressure_j, &Enthalpy_j, Normal, Proj_flux_tensor_j);
  
  /*--- Compute P and Lambda (do it with the Normal) ---*/
  GetPMatrix(&RoeDensity, RoeVelocity, &RoeSoundSpeed, UnitNormal, P_Tensor);
  
  ProjVelocity = 0.0; ProjVelocity_i = 0.0; ProjVelocity_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    ProjVelocity   += RoeVelocity[iDim]*UnitNormal[iDim];
    ProjVelocity_i += Velocity_i[iDim]*UnitNormal[iDim];
    ProjVelocity_j += Velocity_j[iDim]*UnitNormal[iDim];
  }
  
  /*--- Flow eigenvalues and entropy correctors ---*/
  for (iDim = 0; iDim < nDim; iDim++)
    Lambda[iDim] = ProjVelocity;
  
  Lambda[nVar-2] = ProjVelocity + RoeSoundSpeed;
  Lambda[nVar-1] = ProjVelocity - RoeSoundSpeed;
  
  //	/*--- Harten and Hyman (1983) entropy correction ---*/
  //	for (iDim = 0; iDim < nDim; iDim++)
  //		Epsilon[iDim] = 4.0*max(0.0, max(Lambda[iDim]-ProjVelocity_i,ProjVelocity_j-Lambda[iDim]));
  //
  //	Epsilon[nVar-2] = 4.0*max(0.0, max(Lambda[nVar-2]-(ProjVelocity_i+SoundSpeed_i),(ProjVelocity_j+SoundSpeed_j)-Lambda[nVar-2]));
  //	Epsilon[nVar-1] = 4.0*max(0.0, max(Lambda[nVar-1]-(ProjVelocity_i-SoundSpeed_i),(ProjVelocity_j-SoundSpeed_j)-Lambda[nVar-1]));
  //
  //	for (iVar = 0; iVar < nVar; iVar++)
  //		if ( fabs(Lambda[iVar]) < Epsilon[iVar] )
  //			Lambda[iVar] = (Lambda[iVar]*Lambda[iVar] + Epsilon[iVar]*Epsilon[iVar])/(2.0*Epsilon[iVar]);
  //		else
  //			Lambda[iVar] = fabs(Lambda[iVar]);
  
  for (iVar = 0; iVar < nVar; iVar++)
    Lambda[iVar] = fabs(Lambda[iVar]);
  
  if (!implicit) {
    
    /*--- Compute wave amplitudes (characteristics) ---*/
    proj_delta_vel = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      delta_vel[iDim] = Velocity_j[iDim] - Velocity_i[iDim];
      proj_delta_vel += delta_vel[iDim]*Normal[iDim];
    }
    delta_p = Pressure_j - Pressure_i;
    delta_rho = Density_j - Density_i;
    proj_delta_vel = proj_delta_vel/Area;
    
    if (nDim == 2) {
      delta_wave[0] = delta_rho - delta_p/(RoeSoundSpeed*RoeSoundSpeed);
      delta_wave[1] = UnitNormal[1]*delta_vel[0]-UnitNormal[0]*delta_vel[1];
      delta_wave[2] = proj_delta_vel + delta_p/(RoeDensity*RoeSoundSpeed);
      delta_wave[3] = -proj_delta_vel + delta_p/(RoeDensity*RoeSoundSpeed);
    } else {
      delta_wave[0] = delta_rho - delta_p/(RoeSoundSpeed*RoeSoundSpeed);
      delta_wave[1] = UnitNormal[0]*delta_vel[2]-UnitNormal[2]*delta_vel[0];
      delta_wave[2] = UnitNormal[1]*delta_vel[0]-UnitNormal[0]*delta_vel[1];
      delta_wave[3] = proj_delta_vel + delta_p/(RoeDensity*RoeSoundSpeed);
      delta_wave[4] = -proj_delta_vel + delta_p/(RoeDensity*RoeSoundSpeed);
    }
    
    /*--- Roe's Flux approximation ---*/
    for (iVar = 0; iVar < nVar; iVar++) {
      val_residual[iVar] = 0.5*(Proj_flux_tensor_i[iVar]+Proj_flux_tensor_j[iVar]);
      for (jVar = 0; jVar < nVar; jVar++)
        val_residual[iVar] -= 0.5*Lambda[jVar]*delta_wave[jVar]*P_Tensor[iVar][jVar]*Area;
    }
  }
  else {
    
    /*--- Compute inverse P ---*/
    GetPMatrix_inv(&RoeDensity, RoeVelocity, &RoeSoundSpeed, UnitNormal, invP_Tensor);
    
    /*--- Jacobians of the inviscid flux, scaled by
     0.5 because val_resconv ~ 0.5*(fc_i+fc_j)*Normal ---*/
    GetInviscidProjJac(Velocity_i, &Energy_i, Normal, 0.5, val_Jacobian_i);
    GetInviscidProjJac(Velocity_j, &Energy_j, Normal, 0.5, val_Jacobian_j);
    
    /*--- Diference variables iPoint and jPoint ---*/
    for (iVar = 0; iVar < nVar; iVar++)
      Diff_U[iVar] = U_j[iVar]-U_i[iVar];
    
    /*--- Roe's Flux approximation ---*/
    for (iVar = 0; iVar < nVar; iVar++) {
      val_residual[iVar] = 0.5*(Proj_flux_tensor_i[iVar]+Proj_flux_tensor_j[iVar]);
      for (jVar = 0; jVar < nVar; jVar++) {
        Proj_ModJac_Tensor_ij = 0.0;
        /*--- Compute |Proj_ModJac_Tensor| = P x |Lambda| x inverse P ---*/
        for (kVar = 0; kVar < nVar; kVar++)
          Proj_ModJac_Tensor_ij += P_Tensor[iVar][kVar]*Lambda[kVar]*invP_Tensor[kVar][jVar];
        val_residual[iVar] -= 0.5*Proj_ModJac_Tensor_ij*Diff_U[jVar]*Area;
        val_Jacobian_i[iVar][jVar] += 0.5*Proj_ModJac_Tensor_ij*Area;
        val_Jacobian_j[iVar][jVar] -= 0.5*Proj_ModJac_Tensor_ij*Area;
      }
    }
    
  }
  
}

void CUpwRoe_Flow::GetInviscidProjFlux(double *val_density,
                                        double *val_velocity,
                                        double *val_pressure,
                                        double *val_enthalpy,
                                        double *val_normal,
                                        double *val_Proj_Flux) {
  double rhou, rhov, rhow;
  
  
  if (nDim == 2) {
    rhou = (*val_density)*val_velocity[0];
    rhov = (*val_density)*val_velocity[1];
    
    val_Proj_Flux[0] = rhou*val_normal[0];
    val_Proj_Flux[1] = (rhou*val_velocity[0]+(*val_pressure))*val_normal[0];
    val_Proj_Flux[2] = rhou*val_velocity[1]*val_normal[0];
    val_Proj_Flux[3] = rhou*(*val_enthalpy)*val_normal[0];
    
    val_Proj_Flux[0] += rhov*val_normal[1];
    val_Proj_Flux[1] += rhov*val_velocity[0]*val_normal[1];
    val_Proj_Flux[2] += (rhov*val_velocity[1]+(*val_pressure))*val_normal[1];
    val_Proj_Flux[3] += rhov*(*val_enthalpy)*val_normal[1];
  }
  else {
    rhou = (*val_density)*val_velocity[0];
    rhov = (*val_density)*val_velocity[1];
    rhow = (*val_density)*val_velocity[2];
    
    val_Proj_Flux[0] = rhou*val_normal[0];
    val_Proj_Flux[1] = (rhou*val_velocity[0]+(*val_pressure))*val_normal[0];
    val_Proj_Flux[2] = rhou*val_velocity[1]*val_normal[0];
    val_Proj_Flux[3] = rhou*val_velocity[2]*val_normal[0];
    val_Proj_Flux[4] = rhou*(*val_enthalpy)*val_normal[0];
    
    val_Proj_Flux[0] += rhov*val_normal[1];
    val_Proj_Flux[1] += rhov*val_velocity[0]*val_normal[1];
    val_Proj_Flux[2] += (rhov*val_velocity[1]+(*val_pressure))*val_normal[1];
    val_Proj_Flux[3] += rhov*val_velocity[2]*val_normal[1];
    val_Proj_Flux[4] += rhov*(*val_enthalpy)*val_normal[1];
    
    val_Proj_Flux[0] += rhow*val_normal[2];
    val_Proj_Flux[1] += rhow*val_velocity[0]*val_normal[2];
    val_Proj_Flux[2] += rhow*val_velocity[1]*val_normal[2];
    val_Proj_Flux[3] += (rhow*val_velocity[2]+(*val_pressure))*val_normal[2];
    val_Proj_Flux[4] += rhow*(*val_enthalpy)*val_normal[2];
  }
  
}

CUpwRoeTurkel_Flow::CUpwRoeTurkel_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
  Beta_min = config->GetminTurkelBeta();
  Beta_max = config->GetmaxTurkelBeta();
  
  Diff_U = new double [nVar];
  Velocity_i = new double [nDim];
  Velocity_j = new double [nDim];
  RoeVelocity = new double [nDim];
  Proj_flux_tensor_i = new double [nVar];
  Proj_flux_tensor_j = new double [nVar];
  Lambda = new double [nVar];
  Epsilon = new double [nVar];
  absPeJac = new double* [nVar];
  invRinvPe = new double* [nVar];
  R_Tensor  = new double* [nVar];
  Matrix    = new double* [nVar];
  Art_Visc  = new double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    absPeJac[iVar] = new double [nVar];
    invRinvPe[iVar] = new double [nVar];
    Matrix[iVar] = new double [nVar];
    Art_Visc[iVar] = new double [nVar];
    R_Tensor[iVar] = new double [nVar];
  }
}

CUpwRoeTurkel_Flow::~CUpwRoeTurkel_Flow(void) {
  
  delete [] Diff_U;
  delete [] Velocity_i;
  delete [] Velocity_j;
  delete [] RoeVelocity;
  delete [] Proj_flux_tensor_i;
  delete [] Proj_flux_tensor_j;
  delete [] Lambda;
  delete [] Epsilon;
  for (iVar = 0; iVar < nVar; iVar++) {
    delete [] absPeJac[iVar];
    delete [] invRinvPe[iVar];
    delete [] Matrix[iVar];
    delete [] Art_Visc[iVar];
    delete [] R_Tensor[iVar];
  }
  delete [] Matrix;
  delete [] Art_Visc;
  delete [] absPeJac;
  delete [] invRinvPe;
  delete [] R_Tensor;
  
}

void CUpwRoeTurkel_Flow::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {
  
  /*--- Face area (norm or the normal vector) ---*/
  Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    Area += Normal[iDim]*Normal[iDim];
  Area = sqrt(Area);
  
  /*-- Unit Normal ---*/
  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;
  
  /*--- Primitive variables at point i ---*/
  for (iDim = 0; iDim < nDim; iDim++)
    Velocity_i[iDim] = V_i[iDim+1];
  Pressure_i = V_i[nDim+1];
  Density_i = V_i[nDim+2];
  Enthalpy_i = V_i[nDim+3];
  Energy_i = Enthalpy_i - Pressure_i/Density_i;
  
  /*--- Primitive variables at point j ---*/
  for (iDim = 0; iDim < nDim; iDim++)
    Velocity_j[iDim] = V_j[iDim+1];
  Pressure_j = V_j[nDim+1];
  Density_j = V_j[nDim+2];
  Enthalpy_j = V_j[nDim+3];
  Energy_j = Enthalpy_j - Pressure_j/Density_j;
  
  /*--- Recompute conservative variables ---*/
  double U_i[5], U_j[5];
  U_i[0] = Density_i; U_j[0] = Density_j;
  for (iDim = 0; iDim < nDim; iDim++) {
    U_i[iDim+1] = Density_i*Velocity_i[iDim]; U_j[iDim+1] = Density_j*Velocity_j[iDim];
  }
  U_i[nDim+1] = Density_i*Energy_i; U_j[nDim+1] = Density_j*Energy_j;
  
  /*--- Roe-averaged variables at interface between i & j ---*/
  R = sqrt(Density_j/Density_i);
  RoeDensity = R*Density_i;
  sq_vel = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    RoeVelocity[iDim] = (R*Velocity_j[iDim]+Velocity_i[iDim])/(R+1);
    sq_vel += RoeVelocity[iDim]*RoeVelocity[iDim];
  }
  RoeEnthalpy = (R*Enthalpy_j+Enthalpy_i)/(R+1);
  RoeSoundSpeed = sqrt((Gamma-1)*(RoeEnthalpy-0.5*sq_vel));
  RoePressure = RoeDensity/Gamma*RoeSoundSpeed*RoeSoundSpeed;
  
  /*--- Compute Proj_flux_tensor_i ---*/
  GetInviscidProjFlux(&Density_i, Velocity_i, &Pressure_i, &Enthalpy_i, Normal, Proj_flux_tensor_i);
  
  /*--- Compute Proj_flux_tensor_j ---*/
  GetInviscidProjFlux(&Density_j, Velocity_j, &Pressure_j, &Enthalpy_j, Normal, Proj_flux_tensor_j);
  
  ProjVelocity = 0.0; ProjVelocity_i = 0.0; ProjVelocity_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    ProjVelocity   += RoeVelocity[iDim]*UnitNormal[iDim];
    ProjVelocity_i += Velocity_i[iDim]*UnitNormal[iDim];
    ProjVelocity_j += Velocity_j[iDim]*UnitNormal[iDim];
  }
  
  /*--- First few flow eigenvalues of A.Normal with the normal---*/
  for (iDim = 0; iDim < nDim; iDim++)
    Lambda[iDim] = ProjVelocity;
  
  local_Mach = sqrt(sq_vel)/RoeSoundSpeed;
  Beta 	   = max(Beta_min,min(local_Mach,Beta_max));
  Beta2 	   = Beta*Beta;
  
  one_m_Betasqr 		   = 1.0 - Beta2;  // 1-Beta*Beta
  one_p_Betasqr 		   = 1.0 + Beta2;  // 1+Beta*Beta
  sqr_one_m_Betasqr_Lam1 = pow((one_m_Betasqr*Lambda[0]),2); // [(1-Beta^2)*Lambda[0]]^2
  sqr_two_Beta_c_Area    = pow(2.0*Beta*RoeSoundSpeed*Area,2); // [2*Beta*c*Area]^2
  
  /*--- The rest of the flow eigenvalues of preconditioned matrix---*/
  Lambda[nVar-2] = 0.5 * ( one_p_Betasqr*Lambda[0] + sqrt( sqr_one_m_Betasqr_Lam1 + sqr_two_Beta_c_Area));
  Lambda[nVar-1] = 0.5 * ( one_p_Betasqr*Lambda[0] - sqrt( sqr_one_m_Betasqr_Lam1 + sqr_two_Beta_c_Area));
  
  s_hat = 1.0/Area * (Lambda[nVar-1] - Lambda[0]*Beta2);
  r_hat = 1.0/Area * (Lambda[nVar-2] - Lambda[0]*Beta2);
  t_hat = 0.5/Area * (Lambda[nVar-1] - Lambda[nVar-2]);
  rhoB2a2 = RoeDensity*Beta2*RoeSoundSpeed*RoeSoundSpeed;
  
  /*--- Diference variables iPoint and jPoint and absolute value of the eigen values---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    Diff_U[iVar] = U_j[iVar]-U_i[iVar];
    Lambda[iVar] = fabs(Lambda[iVar]);
  }
  
  /*--- Compute the absolute Preconditioned Jacobian in entropic Variables (do it with the Unitary Normal) ---*/
  GetPrecondJacobian(Beta2, r_hat, s_hat, t_hat, rhoB2a2, Lambda, UnitNormal, absPeJac);
  
  /*--- Compute the matrix from entropic variables to conserved variables ---*/
  GetinvRinvPe(Beta2, RoeEnthalpy, RoeSoundSpeed, RoeDensity, RoeVelocity, invRinvPe);
  
  /*--- Compute the matrix from entropic variables to conserved variables ---*/
  GetRMatrix(RoePressure, RoeSoundSpeed, RoeDensity, RoeVelocity, R_Tensor);
  
  if (implicit) {
    /*--- Jacobians of the inviscid flux, scaled by
     0.5 because val_residual ~ 0.5*(fc_i+fc_j)*Normal ---*/
    GetInviscidProjJac(Velocity_i, &Energy_i, Normal, 0.5, val_Jacobian_i);
    GetInviscidProjJac(Velocity_j, &Energy_j, Normal, 0.5, val_Jacobian_j);
  }
  
  for (iVar = 0; iVar < nVar; iVar ++){
    for (jVar = 0; jVar < nVar; jVar ++) {
      Matrix[iVar][jVar] = 0.0;
      for (kVar = 0; kVar < nVar; kVar++)
        Matrix[iVar][jVar]  += absPeJac[iVar][kVar]*R_Tensor[kVar][jVar];
    }
  }
  
  for (iVar = 0; iVar < nVar; iVar ++){
    for (jVar = 0; jVar < nVar; jVar ++) {
      Art_Visc[iVar][jVar] = 0.0;
      for (kVar = 0; kVar < nVar; kVar++)
        Art_Visc[iVar][jVar]  += invRinvPe[iVar][kVar]*Matrix[kVar][jVar];
    }
  }
  
  /*--- Roe's Flux approximation ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    val_residual[iVar] = 0.5*(Proj_flux_tensor_i[iVar]+Proj_flux_tensor_j[iVar]);
    for (jVar = 0; jVar < nVar; jVar++) {
      val_residual[iVar] -= 0.5*Art_Visc[iVar][jVar]*Diff_U[jVar];
      if (implicit) {
        val_Jacobian_i[iVar][jVar] += 0.5*Art_Visc[iVar][jVar];
        val_Jacobian_j[iVar][jVar] -= 0.5*Art_Visc[iVar][jVar];
      }
    }
  }
  
}

void CUpwRoeTurkel_Flow::GetInviscidProjFlux(double *val_density,
                                        double *val_velocity,
                                        double *val_pressure,
                                        double *val_enthalpy,
                                        double *val_normal,
                                        double *val_Proj_Flux) {
  double rhou, rhov, rhow;
  
  
  if (nDim == 2) {
    rhou = (*val_density)*val_velocity[0];
    rhov = (*val_density)*val_velocity[1];
    
    val_Proj_Flux[0] = rhou*val_normal[0];
    val_Proj_Flux[1] = (rhou*val_velocity[0]+(*val_pressure))*val_normal[0];
    val_Proj_Flux[2] = rhou*val_velocity[1]*val_normal[0];
    val_Proj_Flux[3] = rhou*(*val_enthalpy)*val_normal[0];
    
    val_Proj_Flux[0] += rhov*val_normal[1];
    val_Proj_Flux[1] += rhov*val_velocity[0]*val_normal[1];
    val_Proj_Flux[2] += (rhov*val_velocity[1]+(*val_pressure))*val_normal[1];
    val_Proj_Flux[3] += rhov*(*val_enthalpy)*val_normal[1];
  }
  else {
    rhou = (*val_density)*val_velocity[0];
    rhov = (*val_density)*val_velocity[1];
    rhow = (*val_density)*val_velocity[2];
    
    val_Proj_Flux[0] = rhou*val_normal[0];
    val_Proj_Flux[1] = (rhou*val_velocity[0]+(*val_pressure))*val_normal[0];
    val_Proj_Flux[2] = rhou*val_velocity[1]*val_normal[0];
    val_Proj_Flux[3] = rhou*val_velocity[2]*val_normal[0];
    val_Proj_Flux[4] = rhou*(*val_enthalpy)*val_normal[0];
    
    val_Proj_Flux[0] += rhov*val_normal[1];
    val_Proj_Flux[1] += rhov*val_velocity[0]*val_normal[1];
    val_Proj_Flux[2] += (rhov*val_velocity[1]+(*val_pressure))*val_normal[1];
    val_Proj_Flux[3] += rhov*val_velocity[2]*val_normal[1];
    val_Proj_Flux[4] += rhov*(*val_enthalpy)*val_normal[1];
    
    val_Proj_Flux[0] += rhow*val_normal[2];
    val_Proj_Flux[1] += rhow*val_velocity[0]*val_normal[2];
    val_Proj_Flux[2] += rhow*val_velocity[1]*val_normal[2];
    val_Proj_Flux[3] += (rhow*val_velocity[2]+(*val_pressure))*val_normal[2];
    val_Proj_Flux[4] += rhow*(*val_enthalpy)*val_normal[2];
  }
  
}

CAvgGrad_Flow::CAvgGrad_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  
  Flux_Tensor = new double* [nVar];
  for (unsigned short iVar = 0; iVar < (nVar); iVar++)
    Flux_Tensor[iVar] = new double [nDim];
  
  Proj_Flux_Tensor = new double [nVar];

  tau = new double* [nDim];
  delta = new double* [nDim];
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    tau[iDim] = new double [nDim];
    delta[iDim] = new double [nDim];
  }
  
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    for (unsigned short jDim = 0; jDim < nDim; jDim++) {
      if (iDim==jDim) delta[iDim][jDim]=1.0;
      else delta[iDim][jDim]=0.0;
    }
  }
  
  /*--- Compressible flow, primitive variables nDim+3, (T,vx,vy,vz,P,rho) ---*/
  PrimVar_i = new double [nDim+3];
  PrimVar_j = new double [nDim+3];
  Mean_PrimVar = new double [nDim+3];
  
  /*--- Compressible flow, primitive gradient variables nDim+3, (T,vx,vy,vz) ---*/
  Mean_GradPrimVar = new double* [nDim+1];
  for (iVar = 0; iVar < nDim+1; iVar++)
    Mean_GradPrimVar[iVar] = new double [nDim];
}

CAvgGrad_Flow::~CAvgGrad_Flow(void) {
  
  delete [] PrimVar_i;
  delete [] PrimVar_j;
  delete [] Mean_PrimVar;
  
  for (iVar = 0; iVar < nDim+1; iVar++)
    delete [] Mean_GradPrimVar[iVar];
  delete [] Mean_GradPrimVar;
  
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    delete [] tau[iDim];
    delete [] delta[iDim];
  }
  delete [] tau;
  delete [] delta;
  
  delete [] Proj_Flux_Tensor;
  
  for (unsigned short iVar = 0; iVar < nDim+3; iVar++) {
    delete [] Flux_Tensor[iVar];
  }
  delete [] Flux_Tensor;
  
}

void CAvgGrad_Flow::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {
  
  /*--- Normalized normal vector ---*/
  Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    Area += Normal[iDim]*Normal[iDim];
  Area = sqrt(Area);
  
  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;
  
  for (iVar = 0; iVar < nDim+3; iVar++) {
    PrimVar_i[iVar] = V_i[iVar];
    PrimVar_j[iVar] = V_j[iVar];
    Mean_PrimVar[iVar] = 0.5*(PrimVar_i[iVar]+PrimVar_j[iVar]);
  }
  
  /*--- Laminar and Eddy viscosity ---*/
  Laminar_Viscosity_i = V_i[nDim+5];  Laminar_Viscosity_j = V_j[nDim+5];
  Eddy_Viscosity_i = V_i[nDim+6];     Eddy_Viscosity_j = V_j[nDim+6];
  
  /*--- Mean Viscosities and turbulent kinetic energy---*/
  Mean_Laminar_Viscosity = 0.5*(Laminar_Viscosity_i + Laminar_Viscosity_j);
  Mean_Eddy_Viscosity    = 0.5*(Eddy_Viscosity_i + Eddy_Viscosity_j);
  Mean_turb_ke           = 0.5*(turb_ke_i + turb_ke_j);
  
  /*--- Mean gradient approximation ---*/
  for (iVar = 0; iVar < nDim+1; iVar++) {
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradPrimVar[iVar][iDim] = 0.5*(PrimVar_Grad_i[iVar][iDim] + PrimVar_Grad_j[iVar][iDim]);
    }
  }
  
  /*--- Get projected flux tensor ---*/
  GetViscousProjFlux(Mean_PrimVar, Mean_GradPrimVar, Mean_turb_ke, Normal, Mean_Laminar_Viscosity, Mean_Eddy_Viscosity);
  
  /*--- Update viscous residual ---*/
  for (iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] = Proj_Flux_Tensor[iVar];
  
  /*--- Compute the implicit part ---*/
  if (implicit) {
    dist_ij = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      dist_ij += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
    dist_ij = sqrt(dist_ij);
    
    if (dist_ij == 0.0) {
      
      for (iVar = 0; iVar < nVar; iVar++) {
        for (jVar = 0; jVar < nVar; jVar++) {
          val_Jacobian_i[iVar][jVar] = 0.0;
          val_Jacobian_j[iVar][jVar] = 0.0;
        }
      }
    }
    else {
      GetViscousProjJacs(Mean_PrimVar, Mean_Laminar_Viscosity, Mean_Eddy_Viscosity,
                         dist_ij, UnitNormal, Area, Proj_Flux_Tensor, val_Jacobian_i, val_Jacobian_j);
    }
    
  }
  
}


void CAvgGrad_Flow::GetViscousFlux(double *val_primvar, double **val_gradprimvar,
                               double val_laminar_viscosity, double val_eddy_viscosity, double val_mach_inf) {
  
  double total_viscosity = val_laminar_viscosity + val_eddy_viscosity;
  double cp = (Gamma / Gamma_Minus_One) * Gas_Constant;
  double heat_flux_factor = cp * (val_laminar_viscosity/PRANDTL + val_eddy_viscosity/PRANDTL_TURB);
  
  double div_vel = 0.0;
  for (unsigned short iDim = 0 ; iDim < nDim; iDim++)
    div_vel += val_gradprimvar[iDim+1][iDim];
  
  for (unsigned short iDim = 0 ; iDim < nDim; iDim++) {
    for (unsigned short jDim = 0 ; jDim < nDim; jDim++) {
      tau[iDim][jDim] = total_viscosity*( val_gradprimvar[jDim+1][iDim] +
                                         val_gradprimvar[iDim+1][jDim] )
      -TWO3*total_viscosity*div_vel*delta[iDim][jDim];
    }
  }
  
  // Gradient of primitive variables -> [Temp vel_x vel_y vel_z Pressure]
  if (nDim == 3) {
    Flux_Tensor[0][0] = 0.0;
    Flux_Tensor[1][0] = tau[0][0];
    Flux_Tensor[2][0] = tau[0][1];
    Flux_Tensor[3][0] = tau[0][2];
    Flux_Tensor[4][0] = tau[0][0]*val_primvar[1] + tau[0][1]*val_primvar[2] + tau[0][2]*val_primvar[3] +
    heat_flux_factor*val_gradprimvar[0][0];
    
    Flux_Tensor[0][1] = 0.0;
    Flux_Tensor[1][1] = tau[1][0];
    Flux_Tensor[2][1] = tau[1][1];
    Flux_Tensor[3][1] = tau[1][2];
    Flux_Tensor[4][1] = tau[1][0]*val_primvar[1] + tau[1][1]*val_primvar[2] + tau[1][2]*val_primvar[3] +
    heat_flux_factor*val_gradprimvar[0][1];
    
    Flux_Tensor[0][2] = 0.0;
    Flux_Tensor[1][2] = tau[2][0];
    Flux_Tensor[2][2] = tau[2][1];
    Flux_Tensor[3][2] = tau[2][2];
    Flux_Tensor[4][2] = tau[2][0]*val_primvar[1] + tau[2][1]*val_primvar[2] + tau[2][2]*val_primvar[3] +
    heat_flux_factor*val_gradprimvar[0][2];
  }
  if (nDim == 2) {
    Flux_Tensor[0][0] = 0.0;
    Flux_Tensor[1][0] = tau[0][0];
    Flux_Tensor[2][0] = tau[0][1];
    Flux_Tensor[3][0] = tau[0][0]*val_primvar[1] + tau[0][1]*val_primvar[2]+
    heat_flux_factor*val_gradprimvar[0][0];
    
    Flux_Tensor[0][1] = 0.0;
    Flux_Tensor[1][1] = tau[1][0];
    Flux_Tensor[2][1] = tau[1][1];
    Flux_Tensor[3][1] = tau[1][0]*val_primvar[1] + tau[1][1]*val_primvar[2]+
    heat_flux_factor*val_gradprimvar[0][1];
  }
}

void CAvgGrad_Flow::GetViscousProjFlux(double *val_primvar,
                                   double **val_gradprimvar, double val_turb_ke,
                                   double *val_normal,
                                   double val_laminar_viscosity,
                                   double val_eddy_viscosity) {
  
  unsigned short iVar, iDim, jDim;
  double total_viscosity, heat_flux_factor, div_vel, cp, Density;
  Density = val_primvar[nDim+2];
  
  total_viscosity = val_laminar_viscosity + val_eddy_viscosity;
  cp = (Gamma / Gamma_Minus_One) * Gas_Constant;
  heat_flux_factor = cp * (val_laminar_viscosity/PRANDTL + val_eddy_viscosity/PRANDTL_TURB);
  
  div_vel = 0.0;
  for (iDim = 0 ; iDim < nDim; iDim++)
    div_vel += val_gradprimvar[iDim+1][iDim];
  
  for (iDim = 0 ; iDim < nDim; iDim++)
    for (jDim = 0 ; jDim < nDim; jDim++)
      tau[iDim][jDim] = total_viscosity*( val_gradprimvar[jDim+1][iDim] + val_gradprimvar[iDim+1][jDim] )
      - TWO3*total_viscosity*div_vel*delta[iDim][jDim]
      - TWO3*Density*val_turb_ke*delta[iDim][jDim];
  
  
  /*--- Gradient of primitive variables -> [Temp vel_x vel_y vel_z Pressure] ---*/
  if (nDim == 2) {
    Flux_Tensor[0][0] = 0.0;
    Flux_Tensor[1][0] = tau[0][0];
    Flux_Tensor[2][0] = tau[0][1];
    Flux_Tensor[3][0] = tau[0][0]*val_primvar[1] + tau[0][1]*val_primvar[2]+
    heat_flux_factor*val_gradprimvar[0][0];
    
    Flux_Tensor[0][1] = 0.0;
    Flux_Tensor[1][1] = tau[1][0];
    Flux_Tensor[2][1] = tau[1][1];
    Flux_Tensor[3][1] = tau[1][0]*val_primvar[1] + tau[1][1]*val_primvar[2]+
    heat_flux_factor*val_gradprimvar[0][1];
  } else {
    Flux_Tensor[0][0] = 0.0;
    Flux_Tensor[1][0] = tau[0][0];
    Flux_Tensor[2][0] = tau[0][1];
    Flux_Tensor[3][0] = tau[0][2];
    Flux_Tensor[4][0] = tau[0][0]*val_primvar[1] + tau[0][1]*val_primvar[2] + tau[0][2]*val_primvar[3] +
    heat_flux_factor*val_gradprimvar[0][0];
    
    Flux_Tensor[0][1] = 0.0;
    Flux_Tensor[1][1] = tau[1][0];
    Flux_Tensor[2][1] = tau[1][1];
    Flux_Tensor[3][1] = tau[1][2];
    Flux_Tensor[4][1] = tau[1][0]*val_primvar[1] + tau[1][1]*val_primvar[2] + tau[1][2]*val_primvar[3] +
    heat_flux_factor*val_gradprimvar[0][1];
    
    Flux_Tensor[0][2] = 0.0;
    Flux_Tensor[1][2] = tau[2][0];
    Flux_Tensor[2][2] = tau[2][1];
    Flux_Tensor[3][2] = tau[2][2];
    Flux_Tensor[4][2] = tau[2][0]*val_primvar[1] + tau[2][1]*val_primvar[2] + tau[2][2]*val_primvar[3] +
    heat_flux_factor*val_gradprimvar[0][2];
  }
  
  for (iVar = 0; iVar < nVar; iVar++) {
    Proj_Flux_Tensor[iVar] = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      Proj_Flux_Tensor[iVar] += Flux_Tensor[iVar][iDim] * val_normal[iDim];
  }
  
}

void CAvgGrad_Flow::SetLaminarViscosity(double val_lam_viscosity_i, double val_lam_viscosity_j) {
	Laminar_Viscosity_i = val_lam_viscosity_i;
	Laminar_Viscosity_j = val_lam_viscosity_j;
}

CAvgGradCorrected_Flow::CAvgGradCorrected_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

  Flux_Tensor = new double* [nVar];
  for (unsigned short iVar = 0; iVar < (nVar); iVar++)
    Flux_Tensor[iVar] = new double [nDim];
  
  Proj_Flux_Tensor = new double [nVar];
  
  tau = new double* [nDim];
  delta = new double* [nDim];
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    tau[iDim] = new double [nDim];
    delta[iDim] = new double [nDim];
  }
  
  /*--- Compressible flow, primitive variables nDim+3, (T,vx,vy,vz,P,rho) ---*/
  PrimVar_i = new double [nDim+3];
  PrimVar_j = new double [nDim+3];
  Mean_PrimVar = new double [nDim+3];
  
  /*--- Compressible flow, primitive gradient variables nDim+1, (T,vx,vy,vz) ---*/
  Proj_Mean_GradPrimVar_Edge = new double [nDim+1];
  Mean_GradPrimVar = new double* [nDim+1];
  for (iVar = 0; iVar < nDim+1; iVar++)
    Mean_GradPrimVar[iVar] = new double [nDim];
  
  Edge_Vector = new double [nDim];
  
}

CAvgGradCorrected_Flow::~CAvgGradCorrected_Flow(void) {
  
  delete [] PrimVar_i;
  delete [] PrimVar_j;
  delete [] Mean_PrimVar;
  delete [] Proj_Mean_GradPrimVar_Edge;
  delete [] Edge_Vector;
  
  for (iVar = 0; iVar < nDim+1; iVar++)
    delete [] Mean_GradPrimVar[iVar];
  delete [] Mean_GradPrimVar;
  
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    delete [] tau[iDim];
    delete [] delta[iDim];
  }
  delete [] tau;
  delete [] delta;
  
  delete [] Proj_Flux_Tensor;
  
  for (unsigned short iVar = 0; iVar < nDim+3; iVar++) {
    delete [] Flux_Tensor[iVar];
  }
  delete [] Flux_Tensor;
  
}

void CAvgGradCorrected_Flow::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {
  
  /*--- Normalized normal vector ---*/
  Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    Area += Normal[iDim]*Normal[iDim];
  Area = sqrt(Area);
  
  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;
  
  /*--- Compute vector going from iPoint to jPoint ---*/
  dist_ij_2 = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
    dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
  }
  
  /*--- Laminar and Eddy viscosity ---*/
  Laminar_Viscosity_i = V_i[nDim+5];  Laminar_Viscosity_j = V_j[nDim+5];
  Eddy_Viscosity_i = V_i[nDim+6];     Eddy_Viscosity_j = V_j[nDim+6];
  
  for (iVar = 0; iVar < nDim+3; iVar++) {
    PrimVar_i[iVar] = V_i[iVar];
    PrimVar_j[iVar] = V_j[iVar];
    Mean_PrimVar[iVar] = 0.5*(PrimVar_i[iVar]+PrimVar_j[iVar]);
  }
  
  /*--- Mean Viscosities and turbulent kinetic energy ---*/
  Mean_Laminar_Viscosity = 0.5*(Laminar_Viscosity_i + Laminar_Viscosity_j);
  Mean_Eddy_Viscosity    = 0.5*(Eddy_Viscosity_i + Eddy_Viscosity_j);
  Mean_turb_ke           = 0.5*(turb_ke_i + turb_ke_j);
  
  /*--- Projection of the mean gradient in the direction of the edge ---*/
  for (iVar = 0; iVar < nDim+1; iVar++) {
    Proj_Mean_GradPrimVar_Edge[iVar] = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradPrimVar[iVar][iDim] = 0.5*(PrimVar_Grad_i[iVar][iDim] + PrimVar_Grad_j[iVar][iDim]);
      Proj_Mean_GradPrimVar_Edge[iVar] += Mean_GradPrimVar[iVar][iDim]*Edge_Vector[iDim];
    }
    if (dist_ij_2 != 0.0) {
      for (iDim = 0; iDim < nDim; iDim++) {
        Mean_GradPrimVar[iVar][iDim] -= (Proj_Mean_GradPrimVar_Edge[iVar] -
                                         (PrimVar_j[iVar]-PrimVar_i[iVar]))*Edge_Vector[iDim] / dist_ij_2;
      }
    }
  }
  
  /*--- Get projected flux tensor ---*/
  GetViscousProjFlux(Mean_PrimVar, Mean_GradPrimVar, Mean_turb_ke, Normal, Mean_Laminar_Viscosity, Mean_Eddy_Viscosity);
  
  /*--- Save residual value ---*/
  for (iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] = Proj_Flux_Tensor[iVar];
  
  /*--- Compute the implicit part ---*/
  if (implicit) {
    
    if (dist_ij_2 == 0.0) {
      for (iVar = 0; iVar < nVar; iVar++) {
        for (jVar = 0; jVar < nVar; jVar++) {
          val_Jacobian_i[iVar][jVar] = 0.0;
          val_Jacobian_j[iVar][jVar] = 0.0;
        }
      }
    }
    else {
      GetViscousProjJacs(Mean_PrimVar, Mean_Laminar_Viscosity, Mean_Eddy_Viscosity,
                         sqrt(dist_ij_2), UnitNormal, Area, Proj_Flux_Tensor, val_Jacobian_i, val_Jacobian_j);
    }
    
  }
  
}


void CAvgGradCorrected_Flow::GetViscousFlux(double *val_primvar, double **val_gradprimvar,
                               double val_laminar_viscosity, double val_eddy_viscosity, double val_mach_inf) {
  
  double total_viscosity = val_laminar_viscosity + val_eddy_viscosity;
  double cp = (Gamma / Gamma_Minus_One) * Gas_Constant;
  double heat_flux_factor = cp * (val_laminar_viscosity/PRANDTL + val_eddy_viscosity/PRANDTL_TURB);
  
  double div_vel = 0.0;
  for (unsigned short iDim = 0 ; iDim < nDim; iDim++)
    div_vel += val_gradprimvar[iDim+1][iDim];
  
  for (unsigned short iDim = 0 ; iDim < nDim; iDim++) {
    for (unsigned short jDim = 0 ; jDim < nDim; jDim++) {
      tau[iDim][jDim] = total_viscosity*( val_gradprimvar[jDim+1][iDim] +
                                         val_gradprimvar[iDim+1][jDim] )
      -TWO3*total_viscosity*div_vel*delta[iDim][jDim];
    }
  }
  
  // Gradient of primitive variables -> [Temp vel_x vel_y vel_z Pressure]
  if (nDim == 3) {
    Flux_Tensor[0][0] = 0.0;
    Flux_Tensor[1][0] = tau[0][0];
    Flux_Tensor[2][0] = tau[0][1];
    Flux_Tensor[3][0] = tau[0][2];
    Flux_Tensor[4][0] = tau[0][0]*val_primvar[1] + tau[0][1]*val_primvar[2] + tau[0][2]*val_primvar[3] +
    heat_flux_factor*val_gradprimvar[0][0];
    
    Flux_Tensor[0][1] = 0.0;
    Flux_Tensor[1][1] = tau[1][0];
    Flux_Tensor[2][1] = tau[1][1];
    Flux_Tensor[3][1] = tau[1][2];
    Flux_Tensor[4][1] = tau[1][0]*val_primvar[1] + tau[1][1]*val_primvar[2] + tau[1][2]*val_primvar[3] +
    heat_flux_factor*val_gradprimvar[0][1];
    
    Flux_Tensor[0][2] = 0.0;
    Flux_Tensor[1][2] = tau[2][0];
    Flux_Tensor[2][2] = tau[2][1];
    Flux_Tensor[3][2] = tau[2][2];
    Flux_Tensor[4][2] = tau[2][0]*val_primvar[1] + tau[2][1]*val_primvar[2] + tau[2][2]*val_primvar[3] +
    heat_flux_factor*val_gradprimvar[0][2];
  }
  if (nDim == 2) {
    Flux_Tensor[0][0] = 0.0;
    Flux_Tensor[1][0] = tau[0][0];
    Flux_Tensor[2][0] = tau[0][1];
    Flux_Tensor[3][0] = tau[0][0]*val_primvar[1] + tau[0][1]*val_primvar[2]+
    heat_flux_factor*val_gradprimvar[0][0];
    
    Flux_Tensor[0][1] = 0.0;
    Flux_Tensor[1][1] = tau[1][0];
    Flux_Tensor[2][1] = tau[1][1];
    Flux_Tensor[3][1] = tau[1][0]*val_primvar[1] + tau[1][1]*val_primvar[2]+
    heat_flux_factor*val_gradprimvar[0][1];
  }
}

void CAvgGradCorrected_Flow::GetViscousProjFlux(double *val_primvar,
                                   double **val_gradprimvar, double val_turb_ke,
                                   double *val_normal,
                                   double val_laminar_viscosity,
                                   double val_eddy_viscosity) {
  
  unsigned short iVar, iDim, jDim;
  double total_viscosity, heat_flux_factor, div_vel, cp, Density;
  Density = val_primvar[nDim+2];
  
  total_viscosity = val_laminar_viscosity + val_eddy_viscosity;
  cp = (Gamma / Gamma_Minus_One) * Gas_Constant;
  heat_flux_factor = cp * (val_laminar_viscosity/PRANDTL + val_eddy_viscosity/PRANDTL_TURB);
  
  div_vel = 0.0;
  for (iDim = 0 ; iDim < nDim; iDim++)
    div_vel += val_gradprimvar[iDim+1][iDim];
  
  for (iDim = 0 ; iDim < nDim; iDim++)
    for (jDim = 0 ; jDim < nDim; jDim++)
      tau[iDim][jDim] = total_viscosity*( val_gradprimvar[jDim+1][iDim] + val_gradprimvar[iDim+1][jDim] )
      - TWO3*total_viscosity*div_vel*delta[iDim][jDim]
      - TWO3*Density*val_turb_ke*delta[iDim][jDim];
  
  
  /*--- Gradient of primitive variables -> [Temp vel_x vel_y vel_z Pressure] ---*/
  if (nDim == 2) {
    Flux_Tensor[0][0] = 0.0;
    Flux_Tensor[1][0] = tau[0][0];
    Flux_Tensor[2][0] = tau[0][1];
    Flux_Tensor[3][0] = tau[0][0]*val_primvar[1] + tau[0][1]*val_primvar[2]+
    heat_flux_factor*val_gradprimvar[0][0];
    
    Flux_Tensor[0][1] = 0.0;
    Flux_Tensor[1][1] = tau[1][0];
    Flux_Tensor[2][1] = tau[1][1];
    Flux_Tensor[3][1] = tau[1][0]*val_primvar[1] + tau[1][1]*val_primvar[2]+
    heat_flux_factor*val_gradprimvar[0][1];
  } else {
    Flux_Tensor[0][0] = 0.0;
    Flux_Tensor[1][0] = tau[0][0];
    Flux_Tensor[2][0] = tau[0][1];
    Flux_Tensor[3][0] = tau[0][2];
    Flux_Tensor[4][0] = tau[0][0]*val_primvar[1] + tau[0][1]*val_primvar[2] + tau[0][2]*val_primvar[3] +
    heat_flux_factor*val_gradprimvar[0][0];
    
    Flux_Tensor[0][1] = 0.0;
    Flux_Tensor[1][1] = tau[1][0];
    Flux_Tensor[2][1] = tau[1][1];
    Flux_Tensor[3][1] = tau[1][2];
    Flux_Tensor[4][1] = tau[1][0]*val_primvar[1] + tau[1][1]*val_primvar[2] + tau[1][2]*val_primvar[3] +
    heat_flux_factor*val_gradprimvar[0][1];
    
    Flux_Tensor[0][2] = 0.0;
    Flux_Tensor[1][2] = tau[2][0];
    Flux_Tensor[2][2] = tau[2][1];
    Flux_Tensor[3][2] = tau[2][2];
    Flux_Tensor[4][2] = tau[2][0]*val_primvar[1] + tau[2][1]*val_primvar[2] + tau[2][2]*val_primvar[3] +
    heat_flux_factor*val_gradprimvar[0][2];
  }
  
  for (iVar = 0; iVar < nVar; iVar++) {
    Proj_Flux_Tensor[iVar] = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      Proj_Flux_Tensor[iVar] += Flux_Tensor[iVar][iDim] * val_normal[iDim];
  }
  
}

void CAvgGradCorrected_Flow::SetLaminarViscosity(double val_lam_viscosity_i, double val_lam_viscosity_j) {
	Laminar_Viscosity_i = val_lam_viscosity_i;
	Laminar_Viscosity_j = val_lam_viscosity_j;
}
