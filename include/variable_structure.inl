/*!
 * \file variable_structure.inl
 * \brief In-Line subroutines of the <i>variable_structure.hpp</i> file.
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

#pragma once

inline bool CVariable::SetDensity(void) { return 0; }

inline void CVariable::SetSolution(unsigned short val_var, double val_solution) { Solution[val_var] = val_solution; }

inline void CVariable::SetUndivided_Laplacian(unsigned short val_var, double val_undivided_laplacian) { Undivided_Laplacian[val_var] = val_undivided_laplacian; }

inline void CVariable::SetAuxVar(double val_auxvar) { AuxVar = val_auxvar; }

inline void CVariable::SetSolution_Old(unsigned short val_var, double val_solution_old) { Solution_Old[val_var] = val_solution_old; }

inline void CVariable::SetLimiter(unsigned short val_var, double val_limiter) { Limiter[val_var] = val_limiter; }

inline void CVariable::SetSolution_Max(unsigned short val_var, double val_solution) { Solution_Max[val_var] = val_solution; }

inline void CVariable::SetSolution_Min(unsigned short val_var, double val_solution) { Solution_Min[val_var] = val_solution; }

inline void CVariable::SetAuxVarGradient(unsigned short iDim, double val_gradient) { Grad_AuxVar[iDim] = val_gradient; }

inline double *CVariable::GetSolution(void) { return Solution; }

inline double *CVariable::GetSolution_Old(void) { return Solution_Old; }

inline double CVariable::GetAuxVar(void) { return AuxVar; }

inline double *CVariable::GetUndivided_Laplacian(void) { return Undivided_Laplacian; }

inline double CVariable::GetUndivided_Laplacian(unsigned short val_var) { return Undivided_Laplacian[val_var]; }

inline double CVariable::GetSolution(unsigned short val_var) { return Solution[val_var]; }

inline double CVariable::GetSolution_Old(unsigned short val_var) { return Solution_Old[val_var]; }

inline double *CVariable::GetResidual_Sum(void) { return Residual_Sum; }

inline double *CVariable::GetResidual_Old(void) { return Residual_Old; }

inline void CVariable::SetGradient(unsigned short val_var, unsigned short val_dim, double val_value) { Gradient[val_var][val_dim] = val_value; }

inline void CVariable::AddGradient(unsigned short val_var, unsigned short val_dim, double val_value) { Gradient[val_var][val_dim] += val_value; }

inline void CVariable::SubtractGradient(unsigned short val_var, unsigned short val_dim, double val_value) { Gradient[val_var][val_dim] -= val_value; }

inline void CVariable::AddAuxVarGradient(unsigned short val_dim, double val_value) { Grad_AuxVar[val_dim] += val_value; }

inline void CVariable::SubtractAuxVarGradient(unsigned short val_dim, double val_value) { Grad_AuxVar[val_dim] -= val_value; }

inline double CVariable::GetGradient(unsigned short val_var, unsigned short val_dim) { return Gradient[val_var][val_dim]; }

inline double CVariable::GetLimiter(unsigned short val_var) { return Limiter[val_var]; }

inline double CVariable::GetSolution_Max(unsigned short val_var) { return Solution_Max[val_var]; }

inline double CVariable::GetSolution_Min(unsigned short val_var) { return Solution_Min[val_var]; }

inline double CVariable::GetPreconditioner_Beta(void) { return 0; }

inline void CVariable::SetPreconditioner_Beta(double val_Beta) { }

inline double **CVariable::GetGradient(void) { return Gradient; }

inline double *CVariable::GetLimiter(void) { return Limiter; }

inline double *CVariable::GetAuxVarGradient(void) { return Grad_AuxVar; }

inline double CVariable::GetAuxVarGradient(unsigned short val_dim) { return Grad_AuxVar[val_dim]; }

inline double *CVariable::GetResTruncError(void) { return Res_TruncError; }

inline void CVariable::SetDelta_Time(double val_delta_time) { Delta_Time = val_delta_time; }

inline double CVariable::GetDelta_Time(void) { return Delta_Time; }

inline void CVariable::SetMax_Lambda(double val_max_lambda) { Max_Lambda = val_max_lambda; }

inline void CVariable::SetMax_Lambda_Inv(double val_max_lambda) { Max_Lambda_Inv = val_max_lambda; }

inline void CVariable::SetMax_Lambda_Visc(double val_max_lambda) { Max_Lambda_Visc = val_max_lambda; }

inline void CVariable::SetLambda(double val_lambda) { Lambda = val_lambda; }

inline void CVariable::AddMax_Lambda(double val_max_lambda) { Max_Lambda += val_max_lambda; }

inline void CVariable::AddMax_Lambda_Inv(double val_max_lambda) { Max_Lambda_Inv += val_max_lambda; }

inline void CVariable::AddMax_Lambda_Visc(double val_max_lambda) { Max_Lambda_Visc += val_max_lambda; }

inline void CVariable::AddLambda(double val_lambda) { Lambda += val_lambda; }

inline double CVariable::GetMax_Lambda(void) { return Max_Lambda; }

inline double CVariable::GetMax_Lambda_Inv(void) { return Max_Lambda_Inv; }

inline double CVariable::GetMax_Lambda_Visc(void) { return Max_Lambda_Visc; }

inline double CVariable::GetLambda(void) { return Lambda; }

inline double CVariable::GetSensor(void) { return Sensor; }

inline void CVariable::SetSensor(double val_sensor) { Sensor = val_sensor; }

inline double CVariable::GetDensity(void) {	return 0; }

inline double CVariable::GetEnergy(void) { return 0; }

inline double CVariable::GetEddyViscosity(void) { return 0; }

inline double CVariable::GetEnthalpy(void) { return 0; }

inline double CVariable::GetPressure(void) { return 0; }

inline double CVariable::GetProjVel(double *val_vector) { return 0; }

inline double CVariable::GetSoundSpeed(void) { return 0; }

inline double CVariable::GetTemperature(void) { return 0; }

inline double CVariable::GetVelocity(unsigned short val_dim) { return 0; }

inline double CVariable::GetVelocity2(void) { return 0; }

inline double CVariable::GetLaminarViscosity(void) { return 0; }

inline double CVariable::GetVorticity(unsigned short val_dim) { return 0; }

inline double CVariable::GetStrainMag(void) { return 0; }

inline void CVariable::SetEnthalpy(void) { }

inline bool CVariable::SetPrimVar_Compressible(CConfig *config) { return true; }

inline bool CVariable::SetPrimVar_Compressible(double eddy_visc, double turb_ke, CConfig *config) { return true; }

inline double CVariable::GetPrimVar(unsigned short val_var) { return 0; }

inline void CVariable::SetPrimVar(unsigned short val_var, double val_prim) { }

inline void CVariable::SetPrimVar(double *val_prim) { }

inline double *CVariable::GetPrimVar(void) { return NULL; }

inline bool CVariable::SetPressure(double Gamma) { return false; }

inline bool CVariable::SetPressure(CConfig *config) { return false; }

inline bool CVariable::SetPressure(double Gamma, double turb_ke) { return false; }

inline void CVariable::SetPressure(void) { }

inline bool CVariable::SetSoundSpeed(CConfig *config) { return false; }

inline bool CVariable::SetSoundSpeed(void) { return false; }

inline bool CVariable::SetSoundSpeed(double Gamma) { return false; }

inline bool CVariable::SetTemperature(double Gas_Constant) { return false; }

inline bool CVariable::SetTemperature(CConfig *config) { return false; }

inline void CVariable::SetPrimVar(CConfig *config) { }

inline void CVariable::SetPrimVar(CConfig *config, double *Coord) { }

inline void CVariable::SetVelocity(void) { }

inline void CVariable::SetVelocity2(void) { }

inline void CVariable::SetVelocity_Old(double *val_velocity) { }

inline void CVariable::SetLaminarViscosity(CConfig *config) { }

inline void CVariable::SetEddyViscosity(double eddy_visc) { }

inline void CVariable::SetVorticity(void) { }

inline void CVariable::SetStrainMag(void) { }

inline void CVariable::SetGradient_PrimitiveZero(unsigned short val_primvar) { }

inline void CVariable::AddGradient_Primitive(unsigned short val_var, unsigned short val_dim, double val_value) { }

inline void CVariable::SubtractGradient_Primitive(unsigned short val_var, unsigned short val_dim, double val_value) { }

inline double CVariable::GetGradient_Primitive(unsigned short val_var, unsigned short val_dim) { return 0; }

inline double CVariable::GetLimiter_Primitive(unsigned short val_var) { return 0; }

inline void CVariable::SetGradient_Primitive(unsigned short val_var, unsigned short val_dim, double val_value) { }

inline void CVariable::SetLimiter_Primitive(unsigned short val_var, double val_value) { }

inline double **CVariable::GetGradient_Primitive(void) { return NULL; }

inline double *CVariable::GetLimiter_Primitive(void) { return NULL; }

inline void CVariable::SetBlendingFunc(double val_viscosity, double val_dist, double val_density) { }

inline double CVariable::GetF1blending(void) { return 0; }

inline double CVariable::GetF2blending(void) { return 0; }

inline double CVariable::GetmuT() { return 0;}

inline void CVariable::SetmuT(double val_muT) { }

inline double CEulerVariable::GetDensity(void) { return Solution[0]; }

inline double CEulerVariable::GetEnergy(void) { return Solution[nVar-1]/Solution[0]; };

inline double CEulerVariable::GetEnthalpy(void) { return Primitive[nDim+3]; }

inline double CEulerVariable::GetPressure(void) { return Primitive[nDim+1]; }

inline double CEulerVariable::GetSoundSpeed(void) { return Primitive[nDim+4]; }

inline double CEulerVariable::GetTemperature(void) { return Primitive[0]; }

inline double CEulerVariable::GetVelocity(unsigned short val_dim) { return Primitive[val_dim+1]; }

inline double CEulerVariable::GetVelocity2(void) { return Velocity2; }

inline bool CEulerVariable::SetDensity(void) {
  Primitive[nDim+2] = Solution[0];
  if (Primitive[nDim+2] > 0.0) return false;
  else return true;
}

inline bool CEulerVariable::SetPressure(double Gamma) {
  Primitive[nDim+1] = (Gamma-1.0)*Solution[0]*(Solution[nVar-1]/Solution[0]-0.5*Velocity2);
  if (Primitive[nDim+1] > 0.0) return false;
  else return true;
}

inline void CEulerVariable::SetVelocity(void) {
  Velocity2 = 0.0;
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    Primitive[iDim+1] = Solution[iDim+1] / Solution[0];
    Velocity2 += Primitive[iDim+1]*Primitive[iDim+1];
  }
}

inline void CEulerVariable::SetEnthalpy(void) {
  Primitive[nDim+3] = (Solution[nVar-1] + Primitive[nDim+1]) / Solution[0];
}

inline bool CEulerVariable::SetSoundSpeed(double Gamma) {
  double radical = Gamma*Primitive[nDim+1]/Solution[0];
  if (radical < 0.0) return true;
  else {
    Primitive[nDim+4] = sqrt(radical);
    return false;
  }
}

inline bool CEulerVariable::SetTemperature(double Gas_Constant) {
  Primitive[0] = Primitive[nDim+1] / ( Gas_Constant * Solution[0]);
  if (Primitive[0] > 0.0) return false;
  else return true;
}

inline double CEulerVariable::GetPrimVar(unsigned short val_var) { return Primitive[val_var]; }

inline void CEulerVariable::SetPrimVar(unsigned short val_var, double val_prim) { Primitive[val_var] = val_prim; }

inline void CEulerVariable::SetPrimVar(double *val_prim) {
  for (unsigned short iVar = 0; iVar < nPrimVar; iVar++)
  Primitive[iVar] = val_prim[iVar];
}

inline double *CEulerVariable::GetPrimVar(void) { return Primitive; }

inline void CEulerVariable::SetVelocity_Old(double *val_velocity) {
  for (unsigned short iDim = 0; iDim < nDim; iDim++)
  Solution_Old[iDim+1] = val_velocity[iDim]*Solution[0];
}

inline void CEulerVariable::AddGradient_Primitive(unsigned short val_var, unsigned short val_dim, double val_value) { Gradient_Primitive[val_var][val_dim] += val_value; }

inline void CEulerVariable::SubtractGradient_Primitive(unsigned short val_var, unsigned short val_dim, double val_value) { Gradient_Primitive[val_var][val_dim] -= val_value; }

inline double CEulerVariable::GetGradient_Primitive(unsigned short val_var, unsigned short val_dim) { return Gradient_Primitive[val_var][val_dim]; }

inline double CEulerVariable::GetLimiter_Primitive(unsigned short val_var) { return Limiter_Primitive[val_var]; }

inline void CEulerVariable::SetGradient_Primitive(unsigned short val_var, unsigned short val_dim, double val_value) { Gradient_Primitive[val_var][val_dim] = val_value; }

inline void CEulerVariable::SetLimiter_Primitive(unsigned short val_var, double val_value) { Limiter_Primitive[val_var] = val_value; }

inline double **CEulerVariable::GetGradient_Primitive(void) { return Gradient_Primitive; }

inline double *CEulerVariable::GetLimiter_Primitive(void) { return Limiter_Primitive; }

inline double CEulerVariable::GetPreconditioner_Beta() { return Precond_Beta; }

inline void CEulerVariable::SetPreconditioner_Beta(double val_Beta) { Precond_Beta = val_Beta; }

inline double CNSVariable::GetEddyViscosity(void) { return Primitive[nDim+6]; }

inline double CNSVariable::GetLaminarViscosity(void) { return Primitive[nDim+5]; }

inline double CNSVariable::GetVorticity(unsigned short val_dim) { return Vorticity[val_dim]; }

inline double CNSVariable::GetStrainMag(void) { return StrainMag; }

inline void CNSVariable::SetLaminarViscosity(CConfig *config) {
	double Temperature_Dim = Primitive[0]*Temperature_Ref;
	Primitive[nDim+5] = (1.853E-5*(pow(Temperature_Dim/300.0,3.0/2.0) * (300.0+110.3)/(Temperature_Dim+110.3)))/Viscosity_Ref;
}

inline void CNSVariable::SetEddyViscosity(double eddy_visc) { Primitive[nDim+6] = eddy_visc; }

inline bool CNSVariable::SetPressure(double Gamma, double turb_ke) {
  Primitive[nDim+1] = (Gamma-1.0)*Solution[0]*(Solution[nVar-1]/Solution[0]-0.5*Velocity2 - turb_ke);
  if (Primitive[nDim+1] > 0.0) return false;
  else return true;
}

inline double CTurbSSTVariable::GetF1blending(){ return F1; }

inline double CTurbSSTVariable::GetF2blending(){ return F2; }

inline double CTurbSSTVariable::GetCrossDiff(){ return CDkw; }
