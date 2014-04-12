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


inline void CCentJST_Flow::SetNeighbor(unsigned short val_neighbor_i, unsigned short val_neighbor_j) {
	Neighbor_i = val_neighbor_i;
	Neighbor_j = val_neighbor_j;
}

inline void CCentJST_Flow::SetNormal(double *val_normal) {
  Normal = val_normal;
}

inline void CCentJST_Flow::SetPrimitive(double *val_v_i, double *val_v_j) {
	V_i = val_v_i;
	V_j = val_v_j;
}

inline void CCentJST_Flow::SetSensor( double val_sensor_i, double val_sensor_j) {
	Sensor_i = val_sensor_i;
	Sensor_j = val_sensor_j;
}

inline void CCentJST_Flow::SetLambda(double val_lambda_i, double val_lambda_j) {
	Lambda_i = val_lambda_i;
	Lambda_j = val_lambda_j;
}

inline void CCentJST_Flow::SetUndivided_Laplacian(double *val_und_lapl_i, double *val_und_lapl_j) {
	Und_Lapl_i = val_und_lapl_i;
	Und_Lapl_j = val_und_lapl_j;
}

inline void CCentLax_Flow::SetNeighbor(unsigned short val_neighbor_i, unsigned short val_neighbor_j) {
	Neighbor_i = val_neighbor_i;
	Neighbor_j = val_neighbor_j;
}

inline void CCentLax_Flow::SetNormal(double *val_normal) {
  Normal = val_normal;
}

inline void CCentLax_Flow::SetLambda(double val_lambda_i, double val_lambda_j) {
	Lambda_i = val_lambda_i;
	Lambda_j = val_lambda_j;
}

inline void CCentLax_Flow::SetPrimitive(double *val_v_i, double *val_v_j) {
	V_i = val_v_i;
	V_j = val_v_j;
}

inline void CUpwAUSM_Flow::SetNormal(double *val_normal) {
  Normal = val_normal;
}

inline void CUpwAUSM_Flow::SetPrimitive(double *val_v_i, double *val_v_j) {
	V_i = val_v_i;
	V_j = val_v_j;
}

inline void CUpwHLLC_Flow::SetNormal(double *val_normal) {
  Normal = val_normal;
}

inline void CUpwHLLC_Flow::SetPrimitive(double *val_v_i, double *val_v_j) {
	V_i = val_v_i;
	V_j = val_v_j;
}

inline void CUpwRoe_Flow::SetNormal(double *val_normal) {
  Normal = val_normal;
}

inline void CUpwRoe_Flow::SetPrimitive(double *val_v_i, double *val_v_j) {
	V_i = val_v_i;
	V_j = val_v_j;
}


inline void CUpwRoeTurkel_Flow::SetNormal(double *val_normal) {
  Normal = val_normal;
}

inline void CUpwRoeTurkel_Flow::SetPrimitive(double *val_v_i, double *val_v_j) {
	V_i = val_v_i;
	V_j = val_v_j;
}

inline double CUpwRoeTurkel_Flow::GetPrecond_Beta() { return Beta; }

inline void CAvgGrad_Flow::SetPrimVarGradient(double **val_primvar_grad_i, double **val_primvar_grad_j) {
	PrimVar_Grad_i = val_primvar_grad_i;
	PrimVar_Grad_j = val_primvar_grad_j;
}

inline void CAvgGrad_Flow::SetCoord(double *val_coord_i, double *val_coord_j) {
	Coord_i = val_coord_i;
	Coord_j = val_coord_j;
}

inline void CAvgGrad_Flow::SetNormal(double *val_normal) {
  Normal = val_normal;
}

inline void CAvgGrad_Flow::SetPrimitive(double *val_v_i, double *val_v_j) {
	V_i = val_v_i;
	V_j = val_v_j;
}

inline void CAvgGrad_Flow::SetTurbKineticEnergy(double val_turb_ke_i, double val_turb_ke_j) {
	turb_ke_i = val_turb_ke_i;
	turb_ke_j = val_turb_ke_j;
}

inline void CAvgGrad_Flow::SetEddyViscosity(double val_eddy_viscosity_i, double val_eddy_viscosity_j) {
	Eddy_Viscosity_i=val_eddy_viscosity_i;
	Eddy_Viscosity_j=val_eddy_viscosity_j;
}

inline void CAvgGradCorrected_Flow::SetPrimVarGradient(double **val_primvar_grad_i, double **val_primvar_grad_j) {
	PrimVar_Grad_i = val_primvar_grad_i;
	PrimVar_Grad_j = val_primvar_grad_j;
}

inline void CAvgGradCorrected_Flow::SetCoord(double *val_coord_i, double *val_coord_j) {
	Coord_i = val_coord_i;
	Coord_j = val_coord_j;
}

inline void CAvgGradCorrected_Flow::SetNormal(double *val_normal) {
  Normal = val_normal;
}

inline void CAvgGradCorrected_Flow::SetPrimitive(double *val_v_i, double *val_v_j) {
	V_i = val_v_i;
	V_j = val_v_j;
}

inline void CAvgGradCorrected_Flow::SetTurbKineticEnergy(double val_turb_ke_i, double val_turb_ke_j) {
	turb_ke_i = val_turb_ke_i;
	turb_ke_j = val_turb_ke_j;
}

inline void CAvgGradCorrected_Flow::SetEddyViscosity(double val_eddy_viscosity_i, double val_eddy_viscosity_j) {
	Eddy_Viscosity_i=val_eddy_viscosity_i;
	Eddy_Viscosity_j=val_eddy_viscosity_j;
}