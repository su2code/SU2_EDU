/*!
 * \file solver_structure.inl
 * \brief In-Line subroutines of the <i>solver_structure.hpp</i> file.
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

inline void CSolver::SetIterLinSolver(unsigned short val_iterlinsolver) { IterLinSolver = val_iterlinsolver; }

inline unsigned short CSolver::GetIterLinSolver(void) { return IterLinSolver; }

inline void CSolver::GetSurface_Pressure(CGeometry *geometry, CConfig *config) { }

inline void CSolver::SetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long ExtIter) { }

inline void CSolver::SetPrimVar_Gradient_GG(CGeometry *geometry, CConfig *config) { }

inline void CSolver::SetPrimVar_Gradient_LS(CGeometry *geometry, CConfig *config) { }

inline void CSolver::SetPrimVar_Limiter(CGeometry *geometry, CConfig *config) { }

inline void CSolver::SetPreconditioner(CConfig *config, unsigned short iPoint) { }

inline void CSolver::SetDistance(CGeometry *geometry, CConfig *config) { };

inline double CSolver::GetTotal_CLift() { return 0; }

inline double CSolver::GetTotal_CDrag() { return 0; }

inline double CSolver::GetTotal_CMx() { return 0; }

inline double CSolver::GetTotal_CMy() { return 0; }

inline double CSolver::GetTotal_CMz() { return 0; }

inline double CSolver::GetTotal_CFx() { return 0; }

inline double CSolver::GetTotal_CFy() { return 0; }

inline double CSolver::GetTotal_CFz() { return 0; }

inline double CSolver::GetTotal_CSideForce() { return 0; }

inline double CSolver::GetTotal_CEff() { return 0; }

inline double CSolver::GetTotal_CHeat() { return 0; }

inline void CSolver::SetTotal_CLift(double val_Total_CLift) { }

inline void CSolver::SetTotal_CDrag(double val_Total_CDrag) { }

inline double CSolver::GetCPressure(unsigned short val_marker, unsigned short val_vertex) { return 0; }

inline double *CSolver::GetCharacPrimVar(unsigned short val_marker, unsigned short val_vertex) { return 0; }

inline double CSolver::GetCSkinFriction(unsigned short val_marker, unsigned short val_vertex) { return 0; }

inline double CSolver::GetHeatTransferCoeff(unsigned short val_marker, unsigned short val_vertex) { return 0; }

inline double CSolver::GetYPlus(unsigned short val_marker, unsigned short val_vertex) { return 0; }

inline void CSolver::Viscous_Residual(CGeometry *geometry, CSolver **solver_container,
                                      CNumerics *numerics, CConfig *config, unsigned short iMesh, unsigned short iRKstep) { }

inline void CSolver::Source_Residual(CGeometry *geometry, CSolver **solver_container,
                                     CNumerics *numerics, CNumerics *second_numerics, CConfig *config, unsigned short iMesh) { }

inline double CSolver::GetDensity_Inf(void) { return 0; }

inline double CSolver::GetVelocity_Inf(unsigned short val_dim) { return 0; }

inline double CSolver::GetPressure_Inf(void) { return 0; }

inline double* CSolver::GetConstants() {return NULL;}

inline void CSolver::BC_Euler_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                                   unsigned short val_marker) { }

inline void CSolver::BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                                  CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_Sym_Plane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                                  CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                              CConfig *config, unsigned short val_marker) { }

inline void CSolver::BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                               CConfig *config, unsigned short val_marker) { }


inline void CSolver::BC_Supersonic_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                                         CConfig *config, unsigned short val_marker) { }

inline void CSolver::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                  unsigned short iMesh, unsigned long Iteration) { }

inline void CSolver::Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                    unsigned short iMesh) { }

inline void CSolver::Centered_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                       CConfig *config, unsigned short iMesh, unsigned short iRKStep) { }

inline void CSolver::Upwind_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                     CConfig *config, unsigned short iMesh) { }

inline void CSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem) { }

inline void CSolver::SetDissipation_Switch(CGeometry *geometry, CConfig *config) { }

inline void CSolver::SetUndivided_Laplacian(CGeometry *geometry, CConfig *config) { }

inline void CSolver::SetMax_Eigenvalue(CGeometry *geometry, CConfig *config) { }

inline void CSolver::Inviscid_Forces(CGeometry *geometry, CConfig *config) { }

inline void CSolver::Viscous_Forces(CGeometry *geometry, CConfig *config) { }

inline void CSolver::ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) { }

inline void CSolver::SetRes_RMS(unsigned short val_var, double val_residual) { Residual_RMS[val_var] = val_residual; }

inline void CSolver::AddRes_RMS(unsigned short val_var, double val_residual) { Residual_RMS[val_var] += val_residual; }

inline double CSolver::GetRes_RMS(unsigned short val_var) { return Residual_RMS[val_var]; }

inline void CSolver::SetRes_Max(unsigned short val_var, double val_residual, unsigned long val_point) { Residual_Max[val_var] = val_residual; Point_Max[val_var] = val_point; }

inline void CSolver::AddRes_Max(unsigned short val_var, double val_residual, unsigned long val_point) {
  if (val_residual > Residual_Max[val_var]) {
    Residual_Max[val_var] = val_residual;
    Point_Max[val_var] = val_point; }
}

inline double CSolver::GetRes_Max(unsigned short val_var) { return Residual_Max[val_var]; }

inline unsigned long CSolver::GetPoint_Max(unsigned short val_var) { return Point_Max[val_var]; }

inline void CSolver::Set_OldSolution(CGeometry *geometry) {
	for(unsigned long iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
  node[iPoint]->Set_OldSolution();
}

inline unsigned short CSolver::GetnVar(void) { return nVar; }

inline unsigned short CSolver::GetnOutputVariables(void) { return nOutputVariables; }

inline unsigned short CSolver::GetnPrimVar(void) { return nPrimVar; }

inline unsigned short CSolver::GetnPrimVarGrad(void) { return nPrimVarGrad; }

inline double CSolver::GetMax_Delta_Time(void) { return Max_Delta_Time; }

inline double CSolver::GetMin_Delta_Time(void) { return Min_Delta_Time; }

inline double CEulerSolver::GetDensity_Inf(void) { return Density_Inf; }

inline double CEulerSolver::GetVelocity_Inf(unsigned short val_dim) { return Velocity_Inf[val_dim]; }

inline double CEulerSolver::GetPressure_Inf(void) { return Pressure_Inf; }

inline double CEulerSolver::GetCPressure(unsigned short val_marker, unsigned short val_vertex) { return CPressure[val_marker][val_vertex]; }

inline double *CEulerSolver::GetCharacPrimVar(unsigned short val_marker, unsigned short val_vertex) { return CharacPrimVar[val_marker][val_vertex]; }

inline double CEulerSolver::GetTotal_CLift() { return Total_CLift; }

inline double CEulerSolver::GetTotal_CDrag() { return Total_CDrag; }

inline double CEulerSolver::GetTotal_CMx() { return Total_CMx; }

inline double CEulerSolver::GetTotal_CMy() { return Total_CMy; }

inline double CEulerSolver::GetTotal_CMz() { return Total_CMz; }

inline double CEulerSolver::GetTotal_CFx() { return Total_CFx; }

inline double CEulerSolver::GetTotal_CFy() { return Total_CFy; }

inline double CEulerSolver::GetTotal_CFz() { return Total_CFz; }

inline double CEulerSolver::GetTotal_CSideForce() { return Total_CSideForce; }

inline double CEulerSolver::GetTotal_CEff() { return Total_CEff; }

inline void CEulerSolver::SetTotal_CLift(double val_Total_CLift) { Total_CLift = val_Total_CLift; }

inline void CEulerSolver::SetTotal_CDrag(double val_Total_CDrag) { Total_CDrag = val_Total_CDrag; }

inline double CNSSolver::GetCSkinFriction(unsigned short val_marker, unsigned short val_vertex) { return CSkinFriction[val_marker][val_vertex]; }

inline double CNSSolver::GetHeatTransferCoeff(unsigned short val_marker, unsigned short val_vertex) { return CHeatTransfer[val_marker][val_vertex]; }

inline double CNSSolver::GetYPlus(unsigned short val_marker, unsigned short val_vertex) { return YPlus[val_marker][val_vertex]; }
