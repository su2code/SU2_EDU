/*!
 * \file variable_direct_mean.cpp
 * \brief Definition of the solution fields.
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

#include "../include/variable_structure.hpp"

CEulerVariable::CEulerVariable(void) : CVariable() {
  
  /*--- Array initialization ---*/
  
  Primitive          = NULL;
  Gradient_Primitive = NULL;
  Limiter_Primitive  = NULL;
  
}

CEulerVariable::CEulerVariable(double val_density, double *val_velocity,
                               double val_energy, unsigned short val_ndim,
                               unsigned short val_nvar, CConfig *config)
: CVariable(val_ndim, val_nvar, config) {
  
  unsigned short iVar, iDim, iMesh, nMGSmooth = 0;
  
  /*--- Array initialization ---*/
  
  Primitive          = NULL;
  Gradient_Primitive = NULL;
  Limiter_Primitive  = NULL;
  
  /*--- Allocate and initialize the primitive variables and gradients ---*/
  
  nPrimVar = nDim+7; nPrimVarGrad = nDim+4;
  
  /*--- Allocate residual structures ---*/
  
  Res_TruncError = new double [nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Res_TruncError[iVar] = 0.0;
  }
  
  /*--- Residual smoothing (multigrid) ---*/
  
  for (iMesh = 0; iMesh <= config->GetMGLevels(); iMesh++)
    nMGSmooth += config->GetMG_CorrecSmooth(iMesh);
  
  if (nMGSmooth > 0) {
    Residual_Sum = new double [nVar];
    Residual_Old = new double [nVar];
  }
  
  /*--- Allocate undivided laplacian array (centered) ---*/
  
  if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED) {
    Undivided_Laplacian = new double [nVar];
  }
  
  /*--- Allocate arrays for the slope limiter (upwind) ---*/
  
  Limiter_Primitive = new double [nPrimVarGrad];
  for (iVar = 0; iVar < nPrimVarGrad; iVar++) Limiter_Primitive[iVar] = 0.0;
  
  Limiter = new double [nVar];
  for (iVar = 0; iVar < nVar; iVar++) Limiter[iVar] = 0.0;
  
  Solution_Max = new double [nPrimVarGrad];
  Solution_Min = new double [nPrimVarGrad];
  for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
    Solution_Max[iVar] = 0.0;
    Solution_Min[iVar] = 0.0;
  }
  
  /*--- Solution and old solution initialization from input ---*/
  
  Solution[0]     = val_density;
  Solution_Old[0] = val_density;
  for (iDim = 0; iDim < nDim; iDim++) {
    Solution[iDim+1]     = val_density*val_velocity[iDim];
    Solution_Old[iDim+1] = val_density*val_velocity[iDim];
  }
  Solution[nVar-1]     = val_density*val_energy;
  Solution_Old[nVar-1] = val_density*val_energy;
  
  /*--- Incompressible flow, primitive variables nDim+3, (P,vx,vy,vz,rho,beta),
   FreeSurface Incompressible flow, primitive variables nDim+4, (P,vx,vy,vz,rho,beta,dist),
   Compressible flow, primitive variables nDim+5, (T,vx,vy,vz,P,rho,h,c) ---*/
  
  Primitive = new double [nPrimVar];
  for (iVar = 0; iVar < nPrimVar; iVar++) Primitive[iVar] = 0.0;
  
  /*--- Incompressible flow, gradients primitive variables nDim+2, (P,vx,vy,vz,rho),
   FreeSurface Incompressible flow, primitive variables nDim+3, (P,vx,vy,vz,rho,beta,dist),
   Compressible flow, gradients primitive variables nDim+4, (T,vx,vy,vz,P,rho,h)
   We need P, and rho for running the adjoint problem ---*/
  
  Gradient_Primitive = new double* [nPrimVarGrad];
  for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
    Gradient_Primitive[iVar] = new double [nDim];
    for (iDim = 0; iDim < nDim; iDim++)
      Gradient_Primitive[iVar][iDim] = 0.0;
  }
  
}

CEulerVariable::CEulerVariable(double *val_solution, unsigned short val_ndim,
                               unsigned short val_nvar, CConfig *config)
: CVariable(val_ndim, val_nvar, config) {
  
  unsigned short iVar, iDim, iMesh, nMGSmooth = 0;
  
  /*--- Array initialization ---*/
  
  Primitive          = NULL;
  Gradient_Primitive = NULL;
  Limiter_Primitive  = NULL;
  
  /*--- Allocate and initialize the primitive variables and gradients ---*/
  nPrimVar = nDim+7; nPrimVarGrad = nDim+4;
  
  /*--- Allocate residual structures ---*/
  
  Res_TruncError = new double [nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Res_TruncError[iVar] = 0.0;
  }
  
  /*--- Residual smoothing (multigrid) ---*/
  
  for (iMesh = 0; iMesh <= config->GetMGLevels(); iMesh++)
    nMGSmooth += config->GetMG_CorrecSmooth(iMesh);
  
  if (nMGSmooth > 0) {
    Residual_Sum = new double [nVar];
    Residual_Old = new double [nVar];
  }
  
  /*--- Allocate undivided laplacian array (centered)---*/
  
  if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED)
    Undivided_Laplacian = new double [nVar];
  
  /*--- Allocate slope limiter arrays (upwind)---*/
  
  Limiter_Primitive = new double [nPrimVarGrad];
  for (iVar = 0; iVar < nPrimVarGrad; iVar++) Limiter_Primitive[iVar] = 0.0;
  
  Limiter = new double [nVar];
  for (iVar = 0; iVar < nVar; iVar++) Limiter[iVar] = 0.0;
  
  Solution_Max = new double [nPrimVarGrad];
  Solution_Min = new double [nPrimVarGrad];
  for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
    Solution_Max[iVar] = 0.0;
    Solution_Min[iVar] = 0.0;
  }
  
  /*--- Solution initialization from input---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    Solution[iVar]     = val_solution[iVar];
    Solution_Old[iVar] = val_solution[iVar];
  }
  
  /*--- Incompressible flow, primitive variables nDim+3, (P,vx,vy,vz,rho,beta),
   FreeSurface Incompressible flow, primitive variables nDim+4, (P,vx,vy,vz,rho,beta,dist),
   Compressible flow, primitive variables nDim+5, (T,vx,vy,vz,P,rho,h,c) ---*/
  
  Primitive = new double [nPrimVar];
  for (iVar = 0; iVar < nPrimVar; iVar++) Primitive[iVar] = 0.0;
  
  /*--- Incompressible flow, gradients primitive variables nDim+2, (P,vx,vy,vz,rho),
   FreeSurface Incompressible flow, primitive variables nDim+4, (P,vx,vy,vz,rho,beta,dist),
   Compressible flow, gradients primitive variables nDim+4, (T,vx,vy,vz,P,rho,h)
   We need P, and rho for running the adjoint problem ---*/
  
  Gradient_Primitive = new double* [nPrimVarGrad];
  for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
    Gradient_Primitive[iVar] = new double [nDim];
    for (iDim = 0; iDim < nDim; iDim++)
      Gradient_Primitive[iVar][iDim] = 0.0;
  }
  
}

CEulerVariable::~CEulerVariable(void) {
  
  if (Primitive         != NULL) delete [] Primitive;
  if (Limiter_Primitive != NULL) delete [] Limiter_Primitive;
  
  if (Gradient_Primitive != NULL) {
    for (unsigned short iVar = 0; iVar < nPrimVarGrad; iVar++)
      delete Gradient_Primitive[iVar];
    delete [] Gradient_Primitive;
  }
  
}

void CEulerVariable::SetGradient_PrimitiveZero(unsigned short val_primvar) {
  unsigned short iVar, iDim;
  
  for (iVar = 0; iVar < val_primvar; iVar++)
    for (iDim = 0; iDim < nDim; iDim++)
      Gradient_Primitive[iVar][iDim] = 0.0;
}

double CEulerVariable::GetProjVel(double *val_vector) {
  double ProjVel;
  unsigned short iDim;
  
  ProjVel = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    ProjVel += Primitive[iDim+1]*val_vector[iDim];
  
  return ProjVel;
}

bool CEulerVariable::SetPrimVar_Compressible(CConfig *config) {
  
  unsigned short iVar;
  
  bool check_dens = false, check_press = false, check_sos = false;
  bool check_temp = false, physical_solution = true;
  
  double Gas_Constant = config->GetGas_ConstantND();
  double Gamma = config->GetGamma();
  
  /*--- Set the value of the velocity and velocity^2 ---*/
  
  SetVelocity();
  
  /*--- Set the value of the density from the conservative variables ---*/
  
  check_dens = SetDensity();
  
  /*--- Compute pressure (requires calculated velocity^2) ---*/
  
  check_press = SetPressure(Gamma);
  
  /*--- Compute the sound speed (requires calculated pressure) ---*/
  
  check_sos = SetSoundSpeed(Gamma);
  
  /*--- Compute the temperature (requires calculated pressure) ---*/
  
  check_temp = SetTemperature(Gas_Constant);
  
  /*--- Check for and repair any non-physical solutions ---*/
  
  if (check_dens || check_press || check_sos || check_temp) {
    
    /*--- If a variable has become non-physical, use the old solution ---*/
    
    for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = Solution_Old[iVar];
    
    /*--- Recompute the primitive variables using the old solution ---*/
    
    SetVelocity();
    check_dens  = SetDensity();
    check_press = SetPressure(Gamma);
    check_sos   = SetSoundSpeed(Gamma);
    check_temp  = SetTemperature(Gas_Constant);
    
    /*--- Set a flag to keep track of non-physical solutions in the grid ---*/
    
    physical_solution = false;
    
  }
  
  /*--- Compute the enthalpy (requires calculated pressure) ---*/
  
  SetEnthalpy();
  
  return physical_solution;
  
}

CNSVariable::CNSVariable(void) : CEulerVariable() { }

CNSVariable::CNSVariable(double val_density, double *val_velocity,
                         double val_energy, unsigned short val_ndim,
                         unsigned short val_nvar, CConfig *config)
: CEulerVariable(val_density, val_velocity, val_energy, val_ndim, val_nvar, config) {
  
  Temperature_Ref = config->GetTemperature_Ref();
  Viscosity_Ref   = config->GetViscosity_Ref();
  Viscosity_Inf   = config->GetViscosity_FreeStreamND();
  Prandtl_Lam     = config->GetPrandtl_Lam();
  Prandtl_Turb    = config->GetPrandtl_Turb();
  
}

CNSVariable::CNSVariable(double *val_solution, unsigned short val_ndim,
                         unsigned short val_nvar, CConfig *config)
: CEulerVariable(val_solution, val_ndim, val_nvar, config) {
  
  Temperature_Ref = config->GetTemperature_Ref();
  Viscosity_Ref   = config->GetViscosity_Ref();
  Viscosity_Inf   = config->GetViscosity_FreeStreamND();
  Prandtl_Lam     = config->GetPrandtl_Lam();
  Prandtl_Turb    = config->GetPrandtl_Turb();
  
}

CNSVariable::~CNSVariable(void) { }

void CNSVariable::SetVorticity(void) {
  double u_y = Gradient_Primitive[1][1];
  double v_x = Gradient_Primitive[2][0];
  double u_z = 0.0;
  double v_z = 0.0;
  double w_x = 0.0;
  double w_y = 0.0;
  
  if (nDim == 3) {
    u_z = Gradient_Primitive[1][2];
    v_z = Gradient_Primitive[2][2];
    w_x = Gradient_Primitive[3][0];
    w_y = Gradient_Primitive[3][1];
  }
  
  Vorticity[0] = w_y-v_z;
  Vorticity[1] = -(w_x-u_z);
  Vorticity[2] = v_x-u_y;
  
}

void CNSVariable::SetStrainMag(void) {
  
  double div;
  
  if (nDim == 2) {
    div = Gradient_Primitive[1][0] + Gradient_Primitive[2][1];
    StrainMag = 0.0;
    
    // add diagonals
    StrainMag += pow(Gradient_Primitive[1][0] - 1.0/3.0*div, 2.0);
    StrainMag += pow(Gradient_Primitive[2][1] - 1.0/3.0*div, 2.0);
    
    // add off diagonals
    StrainMag += 2.0*pow(0.5*(Gradient_Primitive[1][1] + Gradient_Primitive[2][0]), 2.0);
    
    StrainMag = sqrt(2.0*StrainMag);
    
  }
  else {
    div = Gradient_Primitive[1][0] + Gradient_Primitive[2][1] + Gradient_Primitive[3][2];
    StrainMag = 0.0;
    
    // add diagonals
    StrainMag += pow(Gradient_Primitive[1][0] - 1.0/3.0*div,2.0);
    StrainMag += pow(Gradient_Primitive[2][1] - 1.0/3.0*div,2.0);
    StrainMag += pow(Gradient_Primitive[3][2] - 1.0/3.0*div,2.0);
    
    // add off diagonals
    StrainMag += 2.0*pow(0.5*(Gradient_Primitive[1][1] + Gradient_Primitive[2][0]), 2.0);
    StrainMag += 2.0*pow(0.5*(Gradient_Primitive[1][2] + Gradient_Primitive[3][0]), 2.0);
    StrainMag += 2.0*pow(0.5*(Gradient_Primitive[2][2] + Gradient_Primitive[3][1]), 2.0);
    
    StrainMag = sqrt(2.0*StrainMag);
  }
  
}

bool CNSVariable::SetPrimVar_Compressible(double eddy_visc, double turb_ke, CConfig *config) {
  
  unsigned short iVar;
  
  bool check_dens = false, check_press = false, check_sos = false;
  bool check_temp = false, physical_solution = true;
  
  double Gas_Constant = config->GetGas_ConstantND();
  double Gamma = config->GetGamma();
  
  /*--- Set the value of the velocity and velocity^2 ---*/
  
  SetVelocity();
  
  /*--- Set the value of the density from the conservative variables ---*/
  
  check_dens = SetDensity();
  
  /*--- Compute pressure (requires calculated velocity^2) ---*/
  
  check_press = SetPressure(Gamma, turb_ke);
  
  /*--- Compute the sound speed (requires calculated pressure) ---*/
  
  check_sos = SetSoundSpeed(Gamma);
  
  /*--- Compute the temperature (requires calculated pressure) ---*/
  
  check_temp = SetTemperature(Gas_Constant);
  
  /*--- Check for and repair any non-physical solutions ---*/
  
  if (check_dens || check_press || check_sos || check_temp) {
    
    /*--- If a variable has become non-physical, use the old solution ---*/
    
    for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = Solution_Old[iVar];
    
    /*--- Recompute the primitive variables using the old solution ---*/
    
    SetVelocity();
    check_dens  = SetDensity();
    check_press = SetPressure(Gamma, turb_ke);
    check_sos   = SetSoundSpeed(Gamma);
    check_temp  = SetTemperature(Gas_Constant);
    
    /*--- Set a flag to keep track of non-physical solutions in the grid ---*/
    
    physical_solution = false;
    
  }
  
  /*--- Compute the enthalpy (requires calculated pressure) ---*/
  
  SetEnthalpy();
  
  /*--- Compute the laminar viscosity (requires calculated temperature) ---*/
  
  SetLaminarViscosity(config);
  
  /*--- Set eddy viscosity ---*/
  
  SetEddyViscosity(eddy_visc);
  
  return physical_solution;
  
}
