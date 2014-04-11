/*!
 * \file numerics_structure.inl
 * \brief In-Line subroutines of the <i>numerics_structure.hpp</i> file.
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

inline void CNumerics::ComputeResidual(double *val_residual, CConfig *config) { }

inline void CNumerics::ComputeResidual(double *val_residual_i, double *val_residual_j) { }

inline void CNumerics::ComputeResidual(double *val_residual_i, double *val_residual_j, CConfig *config) { }

inline void CNumerics::ComputeResidual(double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) { }

inline void CNumerics::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j,
                                       CConfig *config) { }

inline void CNumerics::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j,
                                       double **val_JacobianMeanFlow_i, double **val_JacobianMeanFlow_j, CConfig *config) { }

inline void CNumerics::ComputeResidual(double *val_resconv, double *val_resvisc, double **val_Jacobian_i,
                                       double **val_Jacobian_j, CConfig *config) { }

inline void CNumerics::ComputeResidual(double *val_residual_i, double *val_residual_j,
                                       double **val_Jacobian_ii, double **val_Jacobian_ij,
                                       double **val_Jacobian_ji, double **val_Jacobian_jj, CConfig *config) { }

inline void CNumerics::ComputeResidual(double *val_resconv_i, double *val_resvisc_i, double *val_resconv_j,
                                       double *val_resvisc_j, double **val_Jacobian_ii, double **val_Jacobian_ij,
                                       double **val_Jacobian_ji, double **val_Jacobian_jj, CConfig *config) { }

inline void CNumerics::ComputeResidual(double **val_stiffmatrix_elem, CConfig *config) { }

inline void CNumerics::ComputeResidual(double *val_residual, double **val_Jacobian_i, CConfig *config) { }

inline void CNumerics::ComputeSourceViscous(double *val_residual, CConfig *config) { }

inline double CNumerics::GetPrecond_Beta() { return 0; }

inline void CNumerics::SetUndivided_Laplacian(double *val_und_lapl_i, double *val_und_lapl_j) {
	Und_Lapl_i = val_und_lapl_i;
	Und_Lapl_j = val_und_lapl_j;
}

inline void CNumerics::SetSensor( double val_sensor_i, double val_sensor_j) {
	Sensor_i = val_sensor_i;
	Sensor_j = val_sensor_j;
}

inline void CNumerics::SetConservative(double *val_u_i, double *val_u_j) {
	U_i = val_u_i;
	U_j = val_u_j;
}

inline void CNumerics::SetPrimitive(double *val_v_i, double *val_v_j) {
	V_i = val_v_i;
	V_j = val_v_j;
}

inline void CNumerics::SetTimeStep(double val_timestep) {TimeStep = val_timestep;}

inline void CNumerics::SetEddyViscosity(double val_eddy_viscosity_i, double val_eddy_viscosity_j) {
	Eddy_Viscosity_i=val_eddy_viscosity_i;
	Eddy_Viscosity_j=val_eddy_viscosity_j;
}

inline void CNumerics::SetProduction(double val_production) { }

inline void CNumerics::SetDestruction(double val_destruction) { }

inline void CNumerics::SetCrossProduction(double val_crossproduction) { }

inline double CNumerics::GetProduction(void) { return 0; }

inline double CNumerics::GetDestruction(void) { return 0; }

inline double CNumerics::GetCrossProduction(void) { return 0; }

inline void CNumerics::SetTurbKineticEnergy(double val_turb_ke_i, double val_turb_ke_j) {
	turb_ke_i = val_turb_ke_i;
	turb_ke_j = val_turb_ke_j;
}

inline void CNumerics::SetDistance(double val_dist_i, double val_dist_j) {
	dist_i = val_dist_i;
	dist_j = val_dist_j;
}

inline void CNumerics::SetTurbVar(double *val_turbvar_i, double *val_turbvar_j) {
	TurbVar_i = val_turbvar_i;
	TurbVar_j = val_turbvar_j;
}

inline void CNumerics::SetTurbVarGradient(double **val_turbvar_grad_i, double **val_turbvar_grad_j) {
	TurbVar_Grad_i = val_turbvar_grad_i;
	TurbVar_Grad_j = val_turbvar_grad_j;
}

inline void CNumerics::SetPrimVarGradient(double **val_primvar_grad_i, double **val_primvar_grad_j) {
	PrimVar_Grad_i = val_primvar_grad_i;
	PrimVar_Grad_j = val_primvar_grad_j;
}

inline void CNumerics::SetConsVarGradient(double **val_consvar_grad_i, double **val_consvar_grad_j) {
	ConsVar_Grad_i = val_consvar_grad_i;
	ConsVar_Grad_j = val_consvar_grad_j;
}

inline void CNumerics::SetConsVarGradient(double **val_consvar_grad) {
	ConsVar_Grad = val_consvar_grad;
}

inline void CNumerics::SetCoord(double *val_coord_i, double *val_coord_j) {
	Coord_i = val_coord_i;
	Coord_j = val_coord_j;
}

inline void CNumerics::SetSoundSpeed(double val_soundspeed_i, double val_soundspeed_j) {
	SoundSpeed_i = val_soundspeed_i;
	SoundSpeed_j = val_soundspeed_j;
}

inline void CNumerics::SetEnthalpy(double val_enthalpy_i, double val_enthalpy_j) {
	Enthalpy_i = val_enthalpy_i;
	Enthalpy_j = val_enthalpy_j;
}

inline void CNumerics::SetLambda(double val_lambda_i, double val_lambda_j) {
	Lambda_i = val_lambda_i;
	Lambda_j = val_lambda_j;
}

inline void CNumerics::SetNeighbor(unsigned short val_neighbor_i, unsigned short val_neighbor_j) {
	Neighbor_i = val_neighbor_i;
	Neighbor_j = val_neighbor_j;
}

inline void CNumerics::SetAuxVarGrad(double *val_auxvargrad_i, double *val_auxvargrad_j) {
	AuxVar_Grad_i = val_auxvargrad_i;
	AuxVar_Grad_j = val_auxvargrad_j;
}

inline void CNumerics::SetNormal(double *val_normal) { Normal = val_normal; }

inline void CNumerics::SetVolume(double val_volume) { Volume = val_volume; }

inline void CSourcePieceWise_TurbSST::SetF1blending(double val_F1_i, double val_F1_j){
	F1_i = val_F1_i;
	F1_j = val_F1_j;
}

inline void CSourcePieceWise_TurbSST::SetF2blending(double val_F2_i, double val_F2_j){
	F2_i = val_F2_i;
	F2_j = val_F2_j;
}

inline void CSourcePieceWise_TurbSST::SetStrainMag(double val_StrainMag_i, double val_StrainMag_j){
	StrainMag = val_StrainMag_i;
}

inline void CSourcePieceWise_TurbSST::SetCrossDiff(double val_CDkw_i, double val_CDkw_j){
	CDkw = val_CDkw_i;
}

inline void CSourcePieceWise_TurbSA::SetIntermittency(double intermittency_in) { intermittency = intermittency_in; }

inline void CSourcePieceWise_TurbSA::SetProduction(double val_production) { Production = val_production; }

inline void CSourcePieceWise_TurbSA::SetDestruction(double val_destruction) { Destruction = val_destruction; }

inline void CSourcePieceWise_TurbSA::SetCrossProduction(double val_crossproduction) { CrossProduction = val_crossproduction; }

inline double CSourcePieceWise_TurbSA::GetProduction(void) { return Production; }

inline double CSourcePieceWise_TurbSA::GetDestruction(void) { return Destruction; }

inline double CSourcePieceWise_TurbSA::GetCrossProduction(void) { return CrossProduction; }

inline double CUpwRoeTurkel_Flow::GetPrecond_Beta() { return Beta; }

inline void CNumerics::ComputeResidual(double **val_Jacobian_i, double *val_Jacobian_mui, double ***val_Jacobian_gradi, CConfig *config) { }

inline void CNumerics::ComputeResidual(double **val_Jacobian_i, double *val_Jacobian_mui, double ***val_Jacobian_gradi,
                                       double **val_Jacobian_j, double *val_Jacobian_muj, double ***val_Jacobian_gradj, CConfig *config) { }
