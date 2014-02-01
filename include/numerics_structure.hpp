/*!
 * \file numerics_structure.hpp
 * \brief Headers of the main subroutines for the dumerical definition of the problem.
 *        The subroutines and functions are in the <i>numerics_structure.cpp</i>,
 *        <i>numerics_convective.cpp</i>, <i>numerics_viscous.cpp</i>, and
 *        <i>numerics_source.cpp</i> files.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 1.0.0
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

#include <cmath>
#include <iostream>
#include <limits>
#include <cstdlib>

#include "config_structure.hpp"

using namespace std;

/*!
 * \class CNumerics
 * \brief Class for defining the numerical methods.
 * \author F. Palacios.
 * \version 1.0.0
 */
class CNumerics {
protected:
	unsigned short nDim, nVar;	/*!< \brief Number of dimensions and variables. */
	unsigned short nSpecies; 	/*!< \brief No of species present in plasma */
	double Gamma;				/*!< \brief Fluid's Gamma constant (ratio of specific heats). */
	double Gamma_Minus_One;		/*!< \brief Fluids's Gamma - 1.0  . */
	double Gas_Constant;		 		/*!< \brief Gas constant. */
  double *Vector; /*!< \brief Auxiliary vector. */
	unsigned short nDiatomics, nMonatomics;
  
public:
	
  double
  **Flux_Tensor,	/*!< \brief Flux tensor (used for viscous and inviscid purposes. */
	*Proj_Flux_Tensor;		/*!< \brief Flux tensor projected in a direction. */
	
  double
  **tau,		/*!< \brief Viscous stress tensor. */
	**delta;			/*!< \brief Identity matrix. */
  double **dVdU; /*!< \brief Transformation matrix from primitive variables, V, to conserved, U. */
  double
  *Diffusion_Coeff_i, /*!< \brief Species diffusion coefficients at point i. */
  *Diffusion_Coeff_j; /*!< \brief Species diffusion coefficients at point j. */
	double Laminar_Viscosity_i,	/*!< \brief Laminar viscosity at point i. */
	Laminar_Viscosity_j,		/*!< \brief Laminar viscosity at point j. */
	Laminar_Viscosity_id,	/*!< \brief Variation of laminar viscosity at point i. */
	Laminar_Viscosity_jd;		/*!< \brief Variation of laminar viscosity at point j. */
  double Thermal_Conductivity_i, /*!< \brief Thermal conductivity at point i. */
  Thermal_Conductivity_j, /*!< \brief Thermal conductivity at point j. */
  Thermal_Conductivity_ve_i, /*!< \brief Thermal conductivity at point i. */
  Thermal_Conductivity_ve_j; /*!< \brief Thermal conductivity at point j. */
  double *Theta_v; /*!< \brief Characteristic vibrational temperature */
	double Eddy_Viscosity_i,	/*!< \brief Eddy viscosity at point i. */
	Eddy_Viscosity_j;			/*!< \brief Eddy viscosity at point j. */
	double turb_ke_i,	/*!< \brief Turbulent kinetic energy at point i. */
	turb_ke_j;			/*!< \brief Turbulent kinetic energy at point j. */
	double Pressure_i,	/*!< \brief Pressure at point i. */
	Pressure_j;			/*!< \brief Pressure at point j. */
	double GravityForce_i,	/*!< \brief Gravity force at point i. */
	GravityForce_j;			/*!< \brief Gravity force at point j. */
	double Density_i,	/*!< \brief Density at point i. */
	Density_j;			/*!< \brief Density at point j. */
	double Lambda_i,	/*!< \brief Spectral radius at point i. */
	Lambda_j;			/*!< \brief Spectral radius at point j. */
	double LambdaComb_i,	/*!< \brief Spectral radius at point i. */
	LambdaComb_j;			/*!< \brief Spectral radius at point j. */
	double SoundSpeed_i,	/*!< \brief Sound speed at point i. */
	SoundSpeed_j;			/*!< \brief Sound speed at point j. */
	double Enthalpy_i,	/*!< \brief Enthalpy at point i. */
	Enthalpy_j;			/*!< \brief Enthalpy at point j. */
	double dist_i,	/*!< \brief Distance of point i to the nearest wall. */
	dist_j;			/*!< \brief Distance of point j to the nearest wall. */
	double Temp_i,	/*!< \brief Temperature at point i. */
	Temp_j;			/*!< \brief Temperature at point j. */
	double *Temp_tr_i, /*!< \brief Temperature transl-rot at point i. */
	*Temp_tr_j;/*!< \brief Temperature transl-rot at point j. */
	double *Temp_vib_i, /*!< \brief Temperature vibrational at point i. */
	*Temp_vib_j;/*!< \brief Temperature vibrational at point j. */
	double *Und_Lapl_i, /*!< \brief Undivided laplacians at point i. */
	*Und_Lapl_j;		/*!< \brief Undivided laplacians at point j. */
	double Sensor_i,	/*!< \brief Pressure sensor at point i. */
	Sensor_j;			/*!< \brief Pressure sensor at point j. */
	double *GridVel_i,	/*!< \brief Grid velocity at point i. */
	*GridVel_j;			/*!< \brief Grid velocity at point j. */
	double *U_i,		/*!< \brief Vector of conservative variables at point i. */
	*U_id,		/*!< \brief Vector of derivative of conservative variables at point i. */
  *UZeroOrder_i,  /*!< \brief Vector of conservative variables at point i without reconstruction. */
	*U_j,				/*!< \brief Vector of conservative variables at point j. */
  *UZeroOrder_j,  /*!< \brief Vector of conservative variables at point j without reconstruction. */
	*U_jd,				/*!< \brief Vector of derivative of conservative variables at point j. */
	*U_0,				/*!< \brief Vector of conservative variables at node 0. */
	*U_1,				/*!< \brief Vector of conservative variables at node 1. */
	*U_2,				/*!< \brief Vector of conservative variables at node 2. */
	*U_3;				/*!< \brief Vector of conservative variables at node 3. */
	double *V_i,		/*!< \brief Vector of primitive variables at point i. */
	*V_j;				/*!< \brief Vector of primitive variables at point j. */
	double *Psi_i,		/*!< \brief Vector of adjoint variables at point i. */
	*Psi_j;				/*!< \brief Vector of adjoint variables at point j. */
	double *DeltaU_i,	/*!< \brief Vector of linearized variables at point i. */
	*DeltaU_j;			/*!< \brief Vector of linearized variables at point j. */
	double *TurbVar_i,	/*!< \brief Vector of turbulent variables at point i. */
	*TurbVar_id,	/*!< \brief Vector of derivative of turbulent variables at point i. */
	*TurbVar_j,			/*!< \brief Vector of turbulent variables at point j. */
	*TurbVar_jd;	/*!< \brief Vector of derivative of turbulent variables at point j. */
	double *TransVar_i,	/*!< \brief Vector of turbulent variables at point i. */
	*TransVar_j;			/*!< \brief Vector of turbulent variables at point j. */
	double *LevelSetVar_i,	/*!< \brief Vector of turbulent variables at point i. */
	*LevelSetVar_j;			/*!< \brief Vector of turbulent variables at point j. */
	double *TurbPsi_i,	/*!< \brief Vector of adjoint turbulent variables at point i. */
	*TurbPsi_j;			/*!< \brief Vector of adjoint turbulent variables at point j. */
	double **ConsVar_Grad_i,	/*!< \brief Gradient of conservative variables at point i. */
	**ConsVar_Grad_j,			/*!< \brief Gradient of conservative variables at point j. */
	**ConsVar_Grad_0,			/*!< \brief Gradient of conservative variables at point 0. */
	**ConsVar_Grad_1,			/*!< \brief Gradient of conservative variables at point 1. */
	**ConsVar_Grad_2,			/*!< \brief Gradient of conservative variables at point 2. */
	**ConsVar_Grad_3,			/*!< \brief Gradient of conservative variables at point 3. */
	**ConsVar_Grad;				/*!< \brief Gradient of conservative variables which is a scalar. */
	double **PrimVar_Grad_i,	/*!< \brief Gradient of primitive variables at point i. */
	**PrimVar_Grad_j;			/*!< \brief Gradient of primitive variables at point j. */
	double **PsiVar_Grad_i,		/*!< \brief Gradient of adjoint variables at point i. */
	**PsiVar_Grad_j;			/*!< \brief Gradient of adjoint variables at point j. */
	double **TurbVar_Grad_i,	/*!< \brief Gradient of turbulent variables at point i. */
	**TurbVar_Grad_j;			/*!< \brief Gradient of turbulent variables at point j. */
	double **TransVar_Grad_i,	/*!< \brief Gradient of turbulent variables at point i. */
	**TransVar_Grad_j;			/*!< \brief Gradient of turbulent variables at point j. */
	double **LevelSetVar_Grad_i,	/*!< \brief Gradient of level set variables at point i. */
	**LevelSetVar_Grad_j;			/*!< \brief Gradient of level set variables at point j. */
	double **TurbPsi_Grad_i,	/*!< \brief Gradient of adjoint turbulent variables at point i. */
	**TurbPsi_Grad_j;			/*!< \brief Gradient of adjoint turbulent variables at point j. */
	double *AuxVar_Grad_i,		/*!< \brief Gradient of an auxiliary variable at point i. */
	*AuxVar_Grad_j;				/*!< \brief Gradient of an auxiliary variable at point i. */
	double *Coord_i,	/*!< \brief Cartesians coordinates of point i. */
	*Coord_j,			/*!< \brief Cartesians coordinates of point j. */
	*Coord_0,			/*!< \brief Cartesians coordinates of point 0 (Galerkin method, triangle). */
	*Coord_1,			/*!< \brief Cartesians coordinates of point 1 (Galerkin method, tetrahedra). */
	*Coord_2,			/*!< \brief Cartesians coordinates of point 2 (Galerkin method, triangle). */
	*Coord_3;			/*!< \brief Cartesians coordinates of point 3 (Galerkin method, tetrahedra). */
	unsigned short Neighbor_i,	/*!< \brief Number of neighbors of the point i. */
	Neighbor_j;					/*!< \brief Number of neighbors of the point j. */
	double *Normal,	/*!< \brief Normal vector, it norm is the area of the face. */
	*UnitNormal,		/*!< \brief Unitary normal vector. */
	*UnitNormald;		/*!< \brief derivatve of unitary normal vector. */
	double TimeStep,		/*!< \brief Time step useful in dual time method. */
	Area,				/*!< \brief Area of the face i-j. */
	Volume;				/*!< \brief Volume of the control volume around point i. */
	double Volume_n,	/*!< \brief Volume of the control volume at time n. */
	Volume_nM1,		/*!< \brief Volume of the control volume at time n-1. */
	Volume_nP1;		/*!< \brief Volume of the control volume at time n+1. */
	double *U_n,	/*!< \brief Vector of conservative variables at time n. */
	*U_nM1,		/*!< \brief Vector of conservative variables at time n-1. */
	*U_nP1;		/*!< \brief Vector of conservative variables at time n+1. */
	double vel2_inf; /*!< \brief value of the square of freestream speed. */
  
  double *l, *m;
  double *dPdU_i, *dPdU_j;
  double *dTdU_i, *dTdU_j;
  double *dTvedU_i, *dTvedU_j;
  double *Ys, **dFdYj, **dFdYi, *sumdFdYih, *sumdFdYjh, *sumdFdYieve, *sumdFdYjeve;
  unsigned short RHOS_INDEX, T_INDEX, TVE_INDEX, VEL_INDEX, P_INDEX,
  RHO_INDEX, H_INDEX, A_INDEX, RHOCVTR_INDEX, RHOCVVE_INDEX;
  
	/*!
	 * \brief Constructor of the class.
	 */
	CNumerics(void);
  
	/*!
	 * \overload
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CNumerics(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
	/*!
	 * \brief Destructor of the class.
	 */
	virtual ~CNumerics(void);
  
	/*!
	 * \brief Compute the determinant of a 3 by 3 matrix.
	 * \param[in] val_matrix 3 by 3 matrix.
	 * \result Determinant of the matrix
	 */
	double Determinant_3x3(double A00, double A01, double A02,
                         double A10, double A11, double A12,
                         double A20, double A21, double A22);
  
	/*!
	 * \brief Set the solution at different times.
	 * \param[in] val_u_nM1 Conservative solution at time n-1.
	 * \param[in] val_u_n Conservative solution at time n.
	 * \param[in] val_u_nP1 Conservative solution at time n+1.
	 */
	void SetPastSol(double *val_u_nM1, double *val_u_n, double *val_u_nP1);
  
	/*!
	 * \brief Set the control volume at different times.
	 * \param[in] val_volume_nM1 - Control volume at time n-1.
	 * \param[in] val_volume_n - Control volume at time n.
	 * \param[in] val_volume_nP1 - Control volume at time n+1.
	 */
	void SetPastVolume(double val_volume_nM1, double val_volume_n, double val_volume_nP1);
  
	/*!
	 * \brief Set the time step.
	 * \param[in] val_timestep - Value of the time step.
	 */
	void SetTimeStep(double val_timestep);
  
	/*!
	 * \brief Get the Preconditioning Beta.
	 * \return val_Beta - Value of the low Mach Preconditioner.
	 */
	virtual double GetPrecond_Beta();
  
	/*!
	 * \brief Set the freestream velocity square.
	 * \param[in] SetVelocity2_Inf - Value of the square of the freestream velocity.
	 */
	void SetVelocity2_Inf(double val_velocity2);
  
	/*!
	 * \brief Set the value of the conservative variables.
	 * \param[in] val_u_i - Value of the conservative variable at point i.
	 * \param[in] val_u_j - Value of the conservative variable at point j.
	 */
	void SetConservative(double *val_u_i, double *val_u_j);
  
  /*!
	 * \brief Set the value of the conservative variables withour reconstruction.
	 * \param[in] val_u_i - Value of the conservative variable at point i.
	 * \param[in] val_u_j - Value of the conservative variable at point j.
	 */
	void SetConservative_ZeroOrder(double *val_u_i, double *val_u_j);
  
	/*!
	 * \brief Set the value of the primitive variables.
	 * \param[in] val_v_i - Value of the primitive variable at point i.
	 * \param[in] val_v_j - Value of the primitive variable at point j.
	 */
	void SetPrimitive(double *val_v_i, double *val_v_j);
  
	/*!
	 * \brief Set the value of the conservative variables.
	 * \param[in] val_u_0 - Value of the conservative variable at point 0.
	 * \param[in] val_u_1 - Value of the conservative variable at point 1.
	 * \param[in] val_u_2 - Value of the conservative variable at point 2.
	 */
	void SetConservative(double *val_u_0, double *val_u_1, double *val_u_2);
  
	/*!
	 * \brief Set the value of the conservative variables.
	 * \param[in] val_u_0 - Value of the conservative variable at point 0.
	 * \param[in] val_u_1 - Value of the conservative variable at point 1.
	 * \param[in] val_u_2 - Value of the conservative variable at point 2.
	 * \param[in] val_u_3 - Value of the conservative variable at point 3.
	 */
	void SetConservative(double *val_u_0, double *val_u_1, double *val_u_2, double *val_u_3);
  
	/*!
	 * \brief Set the gradient of the conservative variables.
	 * \param[in] val_consvar_grad_i - Gradient of the conservative variable at point i.
	 * \param[in] val_consvar_grad_j - Gradient of the conservative variable at point j.
	 */
	void SetConsVarGradient(double **val_consvar_grad_i, double **val_consvar_grad_j);
  
	/*!
	 * \brief Set the gradient of the conservative variables.
	 * \param[in] val_consvar_grad_0 - Gradient of the conservative variable at point 0.
	 * \param[in] val_consvar_grad_1 - Gradient of the conservative variable at point 1.
	 * \param[in] val_consvar_grad_2 - Gradient of the conservative variable at point 2.
	 */
	void SetConsVarGradient(double **val_consvar_grad_0,
                          double **val_consvar_grad_1,
                          double **val_consvar_grad_2);
  
	/*!
	 * \brief Set the gradient of the conservative variables.
	 * \param[in] val_consvar_grad_0 - Gradient of the conservative variable at point 0.
	 * \param[in] val_consvar_grad_1 - Gradient of the conservative variable at point 1.
	 * \param[in] val_consvar_grad_2 - Gradient of the conservative variable at point 2.
	 * \param[in] val_consvar_grad_3 - Gradient of the conservative variable at point 3.
	 */
	void SetConsVarGradient(double **val_consvar_grad_0,
                          double **val_consvar_grad_1,
                          double **val_consvar_grad_2,
                          double **val_consvar_grad_3);
  
	/*!
	 * \brief Set the gradient of the conservative variables.
	 * \param[in] val_consvar_grad - Gradient of the conservative variable which is a scalar.
	 */
	void SetConsVarGradient(double **val_consvar_grad);
  
	/*!
	 * \brief Set the gradient of the primitive variables.
	 * \param[in] val_primvar_grad_i - Gradient of the primitive variable at point i.
	 * \param[in] val_primvar_grad_j - Gradient of the primitive variable at point j.
	 */
	void SetPrimVarGradient(double **val_primvar_grad_i,
                          double **val_primvar_grad_j);
  
	/*!
	 * \brief Set the value of the adjoint variable.
	 * \param[in] val_psi_i - Value of the adjoint variable at point i.
	 * \param[in] val_psi_j - Value of the adjoint variable at point j.
	 */
	void SetAdjointVar(double *val_psi_i, double *val_psi_j);
  
	/*!
	 * \brief Set the value of the linearized conservative variables.
	 * \param[in] val_deltau_i - Value of the linearized conservative variable at point i.
	 * \param[in] val_deltau_j - Value of the linearized conservative variable at point j.
	 */
	void SetLinearizedVar(double *val_deltau_i, double *val_deltau_j);
  
	/*!
	 * \brief Set the gradient of the adjoint variables.
	 * \param[in] val_psivar_grad_i - Gradient of the adjoint variable at point i.
	 * \param[in] val_psivar_grad_j - Gradient of the adjoint variable at point j.
	 */
	void SetAdjointVarGradient(double **val_psivar_grad_i, double **val_psivar_grad_j);
  
	/*!
	 * \brief Set the value of the turbulent variable.
	 * \param[in] val_turbvar_i - Value of the turbulent variable at point i.
	 * \param[in] val_turbvar_j - Value of the turbulent variable at point j.
	 */
	void SetTurbVar(double *val_turbvar_i, double *val_turbvar_j);
  
	/*!
	 * \brief Set the value of the turbulent variable.
	 * \param[in] val_transvar_i - Value of the turbulent variable at point i.
	 * \param[in] val_transvar_j - Value of the turbulent variable at point j.
	 */
	void SetTransVar(double *val_transvar_i, double *val_transvar_j);
  
	/*!
	 * \brief Set the gradient of the turbulent variables.
	 * \param[in] val_turbvar_grad_i - Gradient of the turbulent variable at point i.
	 * \param[in] val_turbvar_grad_j - Gradient of the turbulent variable at point j.
	 */
	void SetTurbVarGradient(double **val_turbvar_grad_i, double **val_turbvar_grad_j);
  
	/*!
	 * \brief Set the gradient of the turbulent variables.
	 * \param[in] val_turbvar_grad_i - Gradient of the turbulent variable at point i.
	 * \param[in] val_turbvar_grad_j - Gradient of the turbulent variable at point j.
	 */
	void SetTransVarGradient(double **val_transvar_grad_i, double **val_transvar_grad_j);
  
	/*!
	 * \brief Set the value of the level set variable.
	 * \param[in] val_levelsetvar_i - Value of the level set variable at point i.
	 * \param[in] val_levelsetvar_j - Value of the level set variable at point j.
	 */
	void SetLevelSetVar(double *val_levelsetvar_i, double *val_levelsetvar_j);
  
	/*!
	 * \brief Set the gradient of the level set variables.
	 * \param[in] val_levelsetvar_grad_i - Gradient of the level set variable at point i.
	 * \param[in] val_levelsetvar_grad_j - Gradient of the level set variable at point j.
	 */
	void SetLevelSetVarGradient(double **val_levelsetvar_grad_i, double **val_levelsetvar_grad_j);
  
	/*!
	 * \brief Set the value of the adjoint turbulent variable.
	 * \param[in] val_turbpsivar_i - Value of the adjoint turbulent variable at point i.
	 * \param[in] val_turbpsivar_j - Value of the adjoint turbulent variable at point j.
	 */
	void SetTurbAdjointVar(double *val_turbpsivar_i, double *val_turbpsivar_j);
  
	/*!
	 * \brief Set the gradient of the adjoint turbulent variables.
	 * \param[in] val_turbpsivar_grad_i - Gradient of the adjoint turbulent variable at point i.
	 * \param[in] val_turbpsivar_grad_j - Gradient of the adjoint turbulent variable at point j.
	 */
	void SetTurbAdjointGradient (double **val_turbpsivar_grad_i, double **val_turbpsivar_grad_j);
  
	/*!
	 * \brief Set the value of the first blending function.
	 * \param[in] val_F1_i - Value of the first Menter blending function at point i.
	 * \param[in] val_F1_j - Value of the first Menter blending function at point j.
	 */
	virtual void SetF1blending(double val_F1_i, double val_F1_j){/* empty */};
  
	/*!
	 * \brief Set the value of the second blending function.
	 * \param[in] val_F1_i - Value of the second Menter blending function at point i.
	 * \param[in] val_F1_j - Value of the second Menter blending function at point j.
	 */
	virtual void SetF2blending(double val_F1_i, double val_F1_j){/* empty */};
  
	/*!
	 * \brief Set the value of the rate of strain magnitude.
	 * \param[in] val_StrainMag_i - Value of the magnitude of rate of strain at point i.
	 * \param[in] val_StrainMag_j - Value of the magnitude of rate of strain at point j.
	 */
	virtual void SetStrainMag(double val_StrainMag_i, double val_StrainMag_j){/* empty */};
  
	/*!
	 * \brief Set the value of the cross diffusion for the SST model.
	 * \param[in] val_CDkw_i - Value of the cross diffusion at point i.
	 * \param[in] val_CDkw_j - Value of the cross diffusion at point j.
	 */
	virtual void SetCrossDiff(double val_CDkw_i, double val_CDkw_j){/* empty */};
  
	/*!
	 * \brief Set the gradient of the auxiliary variables.
	 * \param[in] val_auxvargrad_i - Gradient of the auxiliary variable at point i.
	 * \param[in] val_auxvargrad_j - Gradient of the auxiliary variable at point j.
	 */
	void SetAuxVarGrad(double *val_auxvargrad_i, double *val_auxvargrad_j);
  
  /*!
	 * \brief Set the diffusion coefficient
	 * \param[in] val_diffusioncoeff_i - Value of the diffusion coefficients at i.
	 * \param[in] val_diffusioncoeff_j - Value of the diffusion coefficients at j
	 */
	void SetDiffusionCoeff(double* val_diffusioncoeff_i,
                         double* val_diffusioncoeff_j);
  
	/*!
	 * \brief Set the laminar viscosity.
	 * \param[in] val_laminar_viscosity_i - Value of the laminar viscosity at point i.
	 * \param[in] val_laminar_viscosity_j - Value of the laminar viscosity at point j.
	 */
	void SetLaminarViscosity(double val_laminar_viscosity_i,
                           double val_laminar_viscosity_j);
  
  /*!
	 * \brief Set the thermal conductivity (translational/rotational)
	 * \param[in] val_thermal_conductivity_i - Value of the thermal conductivity at point i.
	 * \param[in] val_thermal_conductivity_j - Value of the thermal conductivity at point j.
	 * \param[in] iSpecies - Value of the species.
	 */
	void SetThermalConductivity(double val_thermal_conductivity_i,
                              double val_thermal_conductivity_j);
  
	/*!
	 * \brief Set the eddy viscosity.
	 * \param[in] val_eddy_viscosity_i - Value of the eddy viscosity at point i.
	 * \param[in] val_eddy_viscosity_j - Value of the eddy viscosity at point j.
	 */
	void SetEddyViscosity(double val_eddy_viscosity_i,
                        double val_eddy_viscosity_j);
  
	/*!
	 * \brief Set the turbulent kinetic energy.
	 * \param[in] val_turb_ke_i - Value of the turbulent kinetic energy at point i.
	 * \param[in] val_turb_ke_j - Value of the turbulent kinetic energy at point j.
	 */
	void SetTurbKineticEnergy(double val_turb_ke_i, double val_turb_ke_j);
  
	/*!
	 * \brief Set the value of the distance from the nearest wall.
	 * \param[in] val_dist_i - Value of of the distance from point i to the nearest wall.
	 * \param[in] val_dist_j - Value of of the distance from point j to the nearest wall.
	 */
	void SetDistance(double val_dist_i, double val_dist_j);
  
	/*!
	 * \brief Set coordinates of the points.
	 * \param[in] val_coord_i - Coordinates of the point i.
	 * \param[in] val_coord_j - Coordinates of the point j.
	 */
	void SetCoord(double *val_coord_i, double *val_coord_j);
  
	/*!
	 * \overload
	 * \param[in] val_coord_0 - Coordinates of the point 0.
	 * \param[in] val_coord_1 - Coordinates of the point 1.
	 * \param[in] val_coord_2 - Coordinates of the point 2.
	 */
	void SetCoord(double *val_coord_0, double *val_coord_1, double *val_coord_2);
  
	/*!
	 * \overload
	 * \param[in] val_coord_0 - Coordinates of the point 0.
	 * \param[in] val_coord_1 - Coordinates of the point 1.
	 * \param[in] val_coord_2 - Coordinates of the point 2.
	 * \param[in] val_coord_3 - Coordinates of the point 3.
	 */
	void SetCoord(double *val_coord_0, double *val_coord_1, double *val_coord_2,
                double *val_coord_3);
  
	/*!
	 * \brief Set the velocity of the computational grid.
	 * \param[in] val_gridvel_i - Grid velocity of the point i.
	 * \param[in] val_gridvel_j - Grid velocity of the point j.
	 */
	void SetGridVel(double *val_gridvel_i, double *val_gridvel_j);
  
  /*!
	 * \brief Set the value of the pressure.
	 * \param[in] val_pressure_i - Value of the pressure at point i.
	 * \param[in] val_pressure_j - Value of the pressure at point j.
	 */
	void SetPressure(double val_pressure_i, double val_pressure_j);
  
	/*!
	 * \brief Set the value of the sound speed.
	 * \param[in] val_soundspeed_i - Value of the sound speed at point i.
	 * \param[in] val_soundspeed_j - Value of the sound speed at point j.
	 */
	void SetSoundSpeed(double val_soundspeed_i, double val_soundspeed_j);
  
	/*!
	 * \brief Set the value of the temperature.
	 * \param[in] val_temp_i - Value of the temperature at point i.
	 * \param[in] val_temp_j - Value of the temperature at point j.
	 */
	void SetTemperature(double val_temp_i, double val_temp_j);
  
	/*!
	 * \brief Set the value of the species pressures.
	 * \param[in] val_pressure_i - Value of the pressure at point i.
	 * \param[in] val_pressure_j - Value of the pressure at point j.
	 */
	void SetPressure(double* val_pressure_i, double* val_pressure_j);
  
	/*!
	 * \brief Set the value of the enthalpy.
	 * \param[in] val_enthalpy_i - Value of the enthalpy at point i.
	 * \param[in] val_enthalpy_j - Value of the enthalpy at point j.
	 */
	void SetEnthalpy(double val_enthalpy_i, double val_enthalpy_j);
  
	/*!
	 * \brief Set the value of the spectral radius.
	 * \param[in] val_lambda_i - Value of the spectral radius at point i.
	 * \param[in] val_lambda_j - Value of the spectral radius at point j.
	 */
	void SetLambda(double val_lambda_i, double val_lambda_j);
  
	/*!
	 * \brief Set the value of undivided laplacian.
	 * \param[in] val_und_lapl_i Undivided laplacian at point i.
	 * \param[in] val_und_lapl_j Undivided laplacian at point j.
	 */
	void SetUndivided_Laplacian(double *val_und_lapl_i, double *val_und_lapl_j);
  
	/*!
	 * \brief Set the value of the pressure sensor.
	 * \param[in] val_sensor_i Pressure sensor at point i.
	 * \param[in] val_sensor_j Pressure sensor at point j.
	 */
	void SetSensor(double val_sensor_i, double val_sensor_j);
  
	/*!
	 * \brief Set the number of neighbor to a point.
	 * \param[in] val_neighbor_i - Number of neighbor to point i.
	 * \param[in] val_neighbor_j - Number of neighbor to point j.
	 */
	void SetNeighbor(unsigned short val_neighbor_i, unsigned short val_neighbor_j);
  
	/*!
	 * \brief Set the value of the normal vector to the face between two points.
	 * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
	 */
	void SetNormal(double *val_normal);
  
	/*!
	 * \brief Set the value of the volume of the control volume.
	 * \param[in] val_volume Volume of the control volume.
	 */
	void SetVolume(double val_volume);
  
  /*!
	 * \brief Retrieves the value of the species density in the primitive variable vector.
	 * \param[in] iRho_s
	 */
  void SetRhosIndex(unsigned short val_Index);
  
  /*!
	 * \brief Retrieves the value of the species density in the primitive variable vector.
	 * \param[in] iRho_s
	 */
  void SetRhoIndex(unsigned short val_Index);
  
  /*!
	 * \brief Retrieves the value of the species density in the primitive variable vector.
	 * \param[in] iRho_s
	 */
  void SetPIndex(unsigned short val_Index);
  
  /*!
	 * \brief Retrieves the value of the species density in the primitive variable vector.
	 * \param[in] iRho_s
	 */
  void SetTIndex(unsigned short val_Index);
  
  /*!
	 * \brief Retrieves the value of the species density in the primitive variable vector.
	 * \param[in] iRho_s
	 */
  void SetTveIndex(unsigned short val_Index);
  
  /*!
	 * \brief Retrieves the value of the velocity index in the primitive variable vector.
	 * \param[in] i(rho*u)
	 */
  void SetVelIndex(unsigned short val_Index);
  
  /*!
	 * \brief Retrieves the value of the species density in the primitive variable vector.
	 * \param[in] iRho_s
	 */
  void SetHIndex(unsigned short val_Index);
  
  /*!
	 * \brief Retrieves the value of the species density in the primitive variable vector.
	 * \param[in] iRho_s
	 */
  void SetAIndex(unsigned short val_Index);
  
  /*!
	 * \brief Retrieves the value of the species density in the primitive variable vector.
	 * \param[in] iRho_s
	 */
  void SetRhoCvtrIndex(unsigned short val_Index);
  
  /*!
	 * \brief Retrieves the value of the species density in the primitive variable vector.
	 * \param[in] iRho_s
	 */
  void SetRhoCvveIndex(unsigned short val_Index);
  
  /*!
	 * \brief Sets the value of the derivative of pressure w.r.t. species density.
	 * \param[in] iRho_s
	 */
  void SetdPdU(double *val_dPdU_i, double *val_dPdU_j);
  
  /*!
	 * \brief Sets the value of the derivative of temperature w.r.t. species density.
	 * \param[in] iRho_s
	 */
  void SetdTdU(double *val_dTdU_i, double *val_dTdU_j);
  
  /*!
	 * \brief Sets the value of the derivative of vib-el. temperature w.r.t. species density.
	 * \param[in] iRho_s
	 */
  void SetdTvedU(double *val_dTvedU_i, double *val_dTvedU_j);
  
	/*!
	 * \brief Get the inviscid fluxes.
	 * \param[in] val_density - Value of the density.
	 * \param[in] val_velocity - Value of the velocity.
	 * \param[in] val_pressure - Value of the pressure.
	 * \param[in] val_enthalpy - Value of the enthalpy.
	 */
	void GetInviscidFlux(double val_density, double *val_velocity, double val_pressure, double val_enthalpy);
  
	/*!
	 * \brief Get the viscous fluxes.
	 * \param[in] val_primvar - Value of the primitive variables.
	 * \param[in] val_gradprimvar - Gradient of the primitive variables.
	 * \param[in] val_laminar_viscosity - Value of the laminar viscosity.
	 * \param[in] val_eddy_viscosity - Value of the eddy viscosity.
	 * \param[in] val_mach_inf - Value of the Mach number at the infinity.
	 */
	void GetViscousFlux(double *val_primvar, double **val_gradprimvar,
                      double val_laminar_viscosity, double val_eddy_viscosity,
                      double val_mach_inf);
  
	/*!
	 * \brief Compute the projected inviscid flux vector.
	 * \param[in] val_density - Pointer to the density.
	 * \param[in] val_velocity - Pointer to the velocity.
	 * \param[in] val_pressure - Pointer to the pressure.
	 * \param[in] val_enthalpy - Pointer to the enthalpy.
	 * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
	 * \param[out] val_Proj_Flux - Pointer to the projected flux.
	 */
	void GetInviscidProjFlux(double *val_density, double *val_velocity,
                           double *val_pressure, double *val_enthalpy,
                           double *val_normal, double *val_Proj_Flux);
  
  /*!
	 * \brief Compute the projected inviscid flux vector.
	 * \param[in] val_U - Conserved variables
	 * \param[in] val_V - Primitive variables
	 * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
	 * \param[out] val_Proj_Flux - Pointer to the projected flux.
	 */
	void GetInviscidProjFlux(double *val_U, double *val_V, double *val_normal,
                           double *val_Proj_Flux);
  
	/*!
	 * \brief Compute the projection of the viscous fluxes into a direction.
	 * \param[in] val_primvar - Primitive variables.
	 * \param[in] val_gradprimvar - Gradient of the primitive variables.
	 * \param[in] val_turb_ke - Turbulent kinetic energy
	 * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
	 * \param[in] val_laminar_viscosity - Laminar viscosity.
	 * \param[in] val_eddy_viscosity - Eddy viscosity.
	 */
  
	void GetViscousProjFlux(double *val_primvar, double **val_gradprimvar,
                          double val_turb_ke, double *val_normal,
                          double val_laminar_viscosity,
                          double val_eddy_viscosity);
  
  
  /*!
	 * * \brief Compute the projection of the viscous fluxes into a direction.
	 * \brief Overloaded function for multiple species viscous calculations
	 * \param[in] val_primvar - Primitive variables.
	 * \param[in] val_gradprimvar - Gradient of the primitive variables.
	 * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
	 * \param[in] val_laminar_viscosity - Laminar viscosity.
	 * \param[in] val_eddy_viscosity - Eddy viscosity.
	 */
	void GetViscousProjFlux(double *val_primvar,
                          double **val_gradprimvar,
                          double *val_normal,
                          double *val_diffusioncoeff,
                          double val_viscosity,
                          double val_therm_conductivity,
                          double val_therm_conductivity_ve,
                          CConfig *config);
  
	/*!
	 * \brief Compute the projection of the inviscid Jacobian matrices.
	 * \param[in] val_velocity Pointer to the velocity.
	 * \param[in] val_energy Value of the energy.
	 * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
	 * \param[in] val_scale - Scale of the projection.
	 * \param[out] val_Proj_Jac_tensor - Pointer to the projected inviscid Jacobian.
	 */
	void GetInviscidProjJac(double *val_velocity, double *val_energy,
                          double *val_normal, double val_scale,
                          double **val_Proj_Jac_tensor);
  
	/*!
	 * \overload
	 * \brief Compute the projection of the inviscid Jacobian matrices.
	 * \param[in] val_velocity Pointer to the velocity.
	 * \param[in] val_energy Value of the energy.
	 * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
	 * \param[in] val_scale - Scale of the projection.
	 * \param[out] val_Proj_Jac_tensor - Pointer to the projected inviscid Jacobian.
	 */
	void GetInviscidProjJac(double **val_velocity, double *val_energy,
                          double *val_normal, double val_scale,
                          double **val_Proj_Jac_tensor);
  
	/*!
	 * \overload
	 * \brief Compute the projection of the inviscid Jacobian matrices for the two-temperature model.
   * \param[in] val_U - Vector conserved variables.
	 * \param[in] val_V - Vector of primitive variables.
   * \param[in] val_dPdU - Vector of partial derivatives of pressure w.r.t. conserved vars.
	 * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
	 * \param[in] val_scale - Scale of the projection.
	 * \param[out] val_Proj_Jac_tensor - Pointer to the projected inviscid Jacobian.
	 */
  void GetInviscidProjJac(double *val_U, double *val_V, double *val_dPdU,
                          double *val_normal, double val_scale,
                          double **val_Proj_Jac_Tensor);
  
	/*!
	 * \brief TSL-Approximation of Viscous NS Jacobians.
	 * \param[in] val_Mean_PrimVar - Mean value of the primitive variables.
	 * \param[in] val_laminar_viscosity - Value of the laminar viscosity.
	 * \param[in] val_eddy_viscosity - Value of the eddy viscosity.
	 * \param[in] val_dist_ij - Distance between the points.
	 * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
	 * \param[in] val_dS - Area of the face between two nodes.
	 * \param[in] val_Proj_Visc_Flux - Pointer to the projected viscous flux.
	 * \param[out] val_Proj_Jac_Tensor_i - Pointer to the projected viscous Jacobian at point i.
	 * \param[out] val_Proj_Jac_Tensor_j - Pointer to the projected viscous Jacobian at point j.
	 */
	void GetViscousProjJacs(double *val_Mean_PrimVar,
                          double val_laminar_viscosity,
                          double val_eddy_viscosity,
                          double val_dist_ij,
                          double *val_normal, double val_dS,
                          double *val_Proj_Visc_Flux,
                          double **val_Proj_Jac_Tensor_i,
                          double **val_Proj_Jac_Tensor_j);
  
  /*!
	 * \brief TSL-Approximation of Viscous NS Jacobians.
	 * \param[in] val_Mean_PrimVar - Mean value of the primitive variables.
	 * \param[in] val_laminar_viscosity - Value of the laminar viscosity.
	 * \param[in] val_eddy_viscosity - Value of the eddy viscosity.
	 * \param[in] val_dist_ij - Distance between the points.
	 * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
	 * \param[in] val_dS - Area of the face between two nodes.
	 * \param[in] val_Proj_Visc_Flux - Pointer to the projected viscous flux.
	 * \param[out] val_Proj_Jac_Tensor_i - Pointer to the projected viscous Jacobian at point i.
	 * \param[out] val_Proj_Jac_Tensor_j - Pointer to the projected viscous Jacobian at point j.
	 */
	void GetViscousProjJacs(double *val_Mean_PrimVar,
                          double *val_diffusion_coeff,
                          double val_laminar_viscosity,
                          double val_thermal_conductivity,
                          double val_thermal_conductivity_ve,
                          double val_dist_ij,
                          double *val_normal, double val_dS,
                          double *val_Proj_Visc_Flux,
                          double **val_Proj_Jac_Tensor_i,
                          double **val_Proj_Jac_Tensor_j,
                          CConfig *config);
  
	/*!
	 * \brief Computation of the matrix P, this matrix diagonalize the conservative Jacobians in
	 *        the form $P^{-1}(A.Normal)P=Lambda$.
	 * \param[in] val_density - Value of the density.
	 * \param[in] val_velocity - Value of the velocity.
	 * \param[in] val_soundspeed - Value of the sound speed.
	 * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
	 * \param[out] val_p_tensor - Pointer to the P matrix.
	 */
	void GetPMatrix(double *val_density, double *val_velocity,
                  double *val_soundspeed, double *val_normal,
                  double **val_p_tensor);
  
	/*!
	 * \overload
	 * \brief Computation of the matrix P, this matrix diagonalize the conservative Jacobians in
	 *        the form $P^{-1}(A.Normal)P=Lambda$.
	 * \param[in] val_density - Value of the density.
	 * \param[in] val_velocity - Value of the velocity.
	 * \param[in] val_soundspeed - Value of the sound speed.
	 * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
	 * \param[out] val_p_tensor - Pointer to the P matrix.
	 */
	void GetPMatrix(double *val_density, double **val_velocity,
                  double *val_soundspeed, double *val_normal,
                  double **val_p_tensor);
  
  /*!
	 * \overload
	 * \brief Computation of the matrix P, this matrix diagonalizes the conservative Jacobians
	 *        in the form $P^{-1}(A.Normal)P=Lambda$.
	 * \param[in] U - Vector of conserved variables (really only need rhoEve)
	 * \param[in] V - Vector of primitive variables
   * \param[in] val_dPdU - Vector of derivatives of pressure w.r.t. conserved vars.
	 * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[in] l - Tangential vector to face.
   * \param[in] m - Tangential vector to face (mutually orthogonal to val_normal & l).
	 * \param[out] val_invp_tensor - Pointer to inverse of the P matrix.
	 */
  void GetPMatrix(double *U, double *V, double *val_dPdU,
                  double *val_normal, double *l, double *m,
                  double **val_p_tensor) ;
  
	/*!
	 * \brief Computation of the matrix Rinv*Pe.
	 * \param[in] Beta2 - A variable in used to define Pe matrix.
	 * \param[in] val_enthalpy - value of the enthalpy.
	 * \param[in] val_soundspeed - value of the sound speed.
	 * \param[in] val_density - value of the density.
	 * \param[in] val_velocity - value of the velocity.
	 * \param[out] val_invR_invPe - Pointer to the matrix of conversion from entropic to conserved variables.
	 */
	void GetinvRinvPe(double Beta2, double val_enthalpy, double val_soundspeed,
                    double val_density, double* val_velocity,
                    double** val_invR_invPe);
  
	/*!
	 * \brief Computation of the matrix R.
	 * \param[in] val_pressure - value of the pressure.
	 * \param[in] val_soundspeed - value of the sound speed.
	 * \param[in] val_density - value of the density.
	 * \param[in] val_velocity - value of the velocity.
	 * \param[out] val_invR_invPe - Pointer to the matrix of conversion from entropic to conserved variables.
	 */
	void GetRMatrix(double val_pressure, double val_soundspeed,
                  double val_density, double* val_velocity,
                  double** val_invR_invPe);
  
	/*!
	 * \brief Computation of the matrix Td, this matrix diagonalize the preconditioned conservative Jacobians
	 *        in the form $Tg |Lambda| Td = Pc{-1}|Pc (A.Normal)|$.
	 * \param[in] Beta2 - A variable in used to define absPeJacobian matrix.
	 * \param[in] r_hat - A variable in used to define absPeJacobian matrix.
	 * \param[in] s_hat - A variable in used to define absPeJacobian matrix.
	 * \param[in] t_hat - A variable in used to define absPeJacobian matrix.
	 * \param[in] rB2a2 - A variable in used to define absPeJacobian matrix.
	 * \param[in] val_Lambda - Eigenvalues of the Preconditioned Jacobian.
	 * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
	 * \param[out] val_absPeJac - Pointer to the Preconditioned Jacobian matrix.
	 */
	void GetPrecondJacobian(double Beta2, double r_hat, double s_hat, double t_hat, double rB2a2, double* val_Lambda, double* val_normal, double** val_absPeJac);
  
	/*!
	 * \brief Computation of the matrix P^{-1}, this matrix diagonalize the conservative Jacobians
	 *        in the form $P^{-1}(A.Normal)P=Lambda$.
	 * \param[in] val_density - Value of the density.
	 * \param[in] val_velocity - Value of the velocity.
	 * \param[in] val_soundspeed - Value of the sound speed.
	 * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
	 * \param[out] val_invp_tensor - Pointer to inverse of the P matrix.
	 */
	void GetPMatrix_inv(double *val_density, double *val_velocity,
                      double *val_soundspeed, double *val_normal,
                      double **val_invp_tensor);
  
	/*!
	 * \overload
	 * \brief Computation of the matrix P^{-1}, this matrix diagonalize the conservative Jacobians
	 *        in the form $P^{-1}(A.Normal)P=Lambda$.
	 * \param[in] val_density - Value of the density.
	 * \param[in] val_velocity - Value of the velocity.
	 * \param[in] val_soundspeed - Value of the sound speed.
	 * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
	 * \param[out] val_invp_tensor - Pointer to inverse of the P matrix.
	 */
	void GetPMatrix_inv(double *val_density, double **val_velocity,
                      double *val_soundspeed, double *val_normal,
                      double **val_invp_tensor);
  
  /*!
	 * \overload
	 * \brief Computation of the matrix P^{-1}, this matrix diagonalizes the conservative Jacobians
   *        in the form $P^{-1}(A.Normal)P=Lambda$.
   * \param[in] U - Vector of conserved variables.
   * \param[in] V - Vector of primitive variables.
   * \param[in] val_dPdU - Vector of derivatives of pressure w.r.t. conserved variables
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   * \param[in] l - Tangential vector to face.
   * \param[in] m - Tangential vector to face (mutually orthogonal to val_normal & l).
   * \param[out] val_invp_tensor - Pointer to inverse of the P matrix.
	 */
  void GetPMatrix_inv(double *U, double *V, double *val_dPdU,
                      double *val_normal, double *l, double *m,
                      double **val_invp_tensor) ;
  
	/*!
	 * \brief Computation of the projected inviscid lambda (eingenvalues).
	 * \param[in] val_velocity - Value of the velocity.
	 * \param[in] val_soundspeed - Value of the sound speed.
	 * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
	 * \param[in] val_Lambda_Vector - Pointer to Lambda matrix.
	 */
	void GetJacInviscidLambda_fabs(double *val_velocity, double val_soundspeed,
                                 double *val_normal, double *val_Lambda_Vector);
  
	/*!
	 * \brief Compute the numerical residual.
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void ComputeResidual(double *val_residual, CConfig *config);
  
	/*!
	 * \overload
	 * \param[out] val_residual_i - Pointer to the total residual at point i.
	 * \param[out] val_residual_j - Pointer to the total residual at point j.
	 */
	virtual void ComputeResidual(double *val_residual_i, double *val_residual_j);
  
  virtual void ComputeResidual_TransLM(double *val_residual,
                                       double **val_Jacobian_i,
                                       double **val_Jacobian_j, CConfig *config,
                                       double &gamma_sep) ;
  
	/*!
	 * \overload
	 * \param[out] val_residual_i - Pointer to the total residual at point i.
	 * \param[out] val_residual_j - Pointer to the total residual at point j.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void ComputeResidual(double *val_residual_i,
                               double *val_residual_j, CConfig *config);
  
	/*!
	 * \overload
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void ComputeResidual(double *val_residual, double **val_Jacobian_i,
                               double **val_Jacobian_j, CConfig *config);
  
  /*!
	 * \overload
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[out] val_JacobianMeanFlow_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[out] val_JacobianMeanFlow_j - Jacobian of the numerical method at node j (implicit computation).
	 * \param[in] config - Definition of the particular problem.
	 */
  virtual void ComputeResidual(double *val_residual, double **val_Jacobian_i,
                               double **val_Jacobian_j,
                               double **val_JacobianMeanFlow_i,
                               double **val_JacobianMeanFlow_j,
                               CConfig *config);
  
	/*!
	 * \overload
	 * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void ComputeResidual(double **val_Jacobian_i, double **val_Jacobian_j,
                               CConfig *config);
  
	/*!
	 * \overload
	 * \param[out] val_resconv - Pointer to the convective residual.
	 * \param[out] val_resvisc - Pointer to the artificial viscosity residual.
	 * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void ComputeResidual(double *val_resconv, double *val_resvisc,
                               double **val_Jacobian_i, double **val_Jacobian_j,
                               CConfig *config);
  
	/*!
	 * \overload
	 * \param[out] val_residual_i - Pointer to the total residual at point i.
	 * \param[out] val_residual_j - Pointer to the total viscosity residual at point j.
	 * \param[out] val_Jacobian_ii - Jacobian of the numerical method at node i (implicit computation) from node i.
	 * \param[out] val_Jacobian_ij - Jacobian of the numerical method at node i (implicit computation) from node j.
	 * \param[out] val_Jacobian_ji - Jacobian of the numerical method at node j (implicit computation) from node i.
	 * \param[out] val_Jacobian_jj - Jacobian of the numerical method at node j (implicit computation) from node j.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void ComputeResidual(double *val_residual_i, double *val_residual_j,
                               double **val_Jacobian_ii,
                               double **val_Jacobian_ij,
                               double **val_Jacobian_ji,
                               double **val_Jacobian_jj, CConfig *config);
  
	/*!
	 * \overload
	 * \param[out] val_resconv_i - Pointer to the convective residual at point i.
	 * \param[out] val_resvisc_i - Pointer to the artificial viscosity residual at point i.
	 * \param[out] val_resconv_j - Pointer to the convective residual at point j.
	 * \param[out] val_resvisc_j - Pointer to the artificial viscosity residual at point j.
	 * \param[out] val_Jacobian_ii - Jacobian of the numerical method at node i (implicit computation) from node i.
	 * \param[out] val_Jacobian_ij - Jacobian of the numerical method at node i (implicit computation) from node j.
	 * \param[out] val_Jacobian_ji - Jacobian of the numerical method at node j (implicit computation) from node i.
	 * \param[out] val_Jacobian_jj - Jacobian of the numerical method at node j (implicit computation) from node j.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void ComputeResidual(double *val_resconv_i, double *val_resvisc_i,
                               double *val_resconv_j, double *val_resvisc_j,
                               double **val_Jacobian_ii,
                               double **val_Jacobian_ij,
                               double **val_Jacobian_ji,
                               double **val_Jacobian_jj, CConfig *config);
  
	/*!
	 * \overload
	 * \param[out] val_stiffmatrix_elem - Stiffness matrix for Galerkin computation.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void ComputeResidual(double **val_stiffmatrix_elem, CConfig *config);
  
	/*!
	 * \overload
	 * \param[in] config - Definition of the particular problem.
	 * \param[out] val_residual - residual of the source terms
	 * \param[out] val_Jacobian_i - Jacobian of the source terms
	 */
	virtual void ComputeResidual(double *val_residual, double **val_Jacobian_i,
                               CConfig *config);
  
	/*!
	 * \overload
	 * \param[out] - Matrix for storing the constants to be used in the calculation of the equilibrium extent of reaction Keq.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void GetEq_Rxn_Coefficients(double **EqnRxnConstants, CConfig *config);
  
	/*!
	 * \brief Residual for source term integration.
	 * \param[out] val_residual - Pointer to the source residual containing chemistry terms.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void ComputeResidual_Axisymmetric(double *val_residual, CConfig *config);
  
	/*!
	 * \brief Residual for source term integration.
	 * \param[out] val_residual - Pointer to the source residual containing chemistry terms.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void ComputeResidual_Axisymmetric_ad(double *val_residual, double *val_residuald, CConfig *config);
  
	/*!
	 * \brief Calculation of axisymmetric source term Jacobian
	 * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void SetJacobian_Axisymmetric(double **val_Jacobian_i, CConfig *config);
  
  /*!
	 * \brief Calculation of the translational-vibrational energy exchange source term
	 * \param[in] config - Definition of the particular problem.
	 * \param[out] val_residual - residual of the source terms
	 * \param[out] val_Jacobian_i - Jacobian of the source terms
	 */
	virtual void ComputeVibRelaxation(double *val_residual, double **val_Jacobian_i, CConfig *config);
  
  /*!
	 * \brief Calculation of the chemistry source term
	 * \param[in] config - Definition of the particular problem.
	 * \param[out] val_residual - residual of the source terms
	 * \param[out] val_Jacobian_i - Jacobian of the source terms
	 */
	virtual void ComputeChemistry(double *val_residual, double **val_Jacobian_i, CConfig *config);
  
  /*!
	 * \brief Calculates constants used for Keq correlation.
	 * \param[out] A - Pointer to coefficient array.
   * \param[in] val_reaction - Reaction number indicator.
	 * \param[in] config - Definition of the particular problem.
	 */
  virtual void GetKeqConstants(double *A, unsigned short val_reaction, CConfig *config);
  
	/*!
	 * \brief Set intermittency for numerics (used in SA with LM transition model)
	 */
	virtual void SetIntermittency(double intermittency_in);
  
  /*!
	 * \brief Computes the viscous source term for the TNE2 adjoint problem
	 * \param[in] config - Definition of the particular problem.
	 * \param[out] val_residual - residual of the source terms
	 */
  virtual void ComputeSourceViscous(double *val_residual, CConfig *config);
  
  /*!
	 * \brief Residual for source term integration.
	 * \param[in] val_production - Value of the Production.
	 */
  virtual void SetProduction(double val_production);
  
  /*!
	 * \brief Residual for source term integration.
	 * \param[in] val_destruction - Value of the Destruction.
	 */
  virtual void SetDestruction(double val_destruction);
  
  /*!
	 * \brief Residual for source term integration.
	 * \param[in] val_crossproduction - Value of the CrossProduction.
	 */
  virtual void SetCrossProduction(double val_crossproduction);
  
  /*!
	 * \brief Residual for source term integration.
	 * \param[in] val_production - Value of the Production.
	 */
  virtual double GetProduction(void);
  
  /*!
	 * \brief Residual for source term integration.
	 * \param[in] val_destruction - Value of the Destruction.
	 */
  virtual double GetDestruction(void);
  
  /*!
	 * \brief Residual for source term integration.
	 * \param[in] val_crossproduction - Value of the CrossProduction.
	 */
  virtual double GetCrossProduction(void);
  
	/*!
	 * \overload
	 * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void ComputeResidual(double **val_Jacobian_i,
                               double *val_Jacobian_mui,
                               double ***val_Jacobian_gradi, CConfig *config);
  
	/*!
	 * \overload
	 * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void ComputeResidual(double **val_Jacobian_i,
                               double *val_Jacobian_mui,
                               double ***val_Jacobian_gradi,
                               double **val_Jacobian_j,
                               double *val_Jacobian_muj,
                               double ***val_Jacobian_gradj, CConfig *config);
  
  /*!
	 * \brief Computing stiffness matrix of the Galerkin method.
	 * \param[out] val_stiffmatrix_elem - Stiffness matrix for Galerkin computation.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void SetFEA_StiffMatrix2D(double **StiffMatrix_Elem, double CoordCorners[8][3], unsigned short nNodes);
  
  /*!
	 * \brief Computing stiffness matrix of the Galerkin method.
	 * \param[out] val_stiffmatrix_elem - Stiffness matrix for Galerkin computation.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void SetFEA_StiffMatrix3D(double **StiffMatrix_Elem, double CoordCorners[8][3], unsigned short nNodes);
  
  /*!
	 * \brief Computes a basis of orthogonal vectors from a suppled vector
	 * \param[in] config - Normal vector
	 */
  void CreateBasis(double *val_Normal);
  
};

/*!
 * \class CUpwRoe_Flow
 * \brief Class for solving an approximate Riemann solver of Roe for the flow equations.
 * \ingroup ConvDiscr
 * \author A. Bueno (UPM) & F. Palacios (Stanford University).
 * \version 1.0.0
 */
class CUpwRoe_Flow : public CNumerics {
private:
	bool implicit, grid_movement;
	double *Diff_U;
	double *Velocity_i, *Velocity_j, *RoeVelocity;
	double *Proj_flux_tensor_i, *Proj_flux_tensor_j;
	double *delta_wave, *delta_vel;
	double *Lambda, *Epsilon;
	double **P_Tensor, **invP_Tensor;
	double sq_vel, Proj_ModJac_Tensor_ij, Density_i, Energy_i, SoundSpeed_i, Pressure_i, Enthalpy_i,
	Density_j, Energy_j, SoundSpeed_j, Pressure_j, Enthalpy_j, R, RoeDensity, RoeEnthalpy, RoeSoundSpeed,
	ProjVelocity, ProjVelocity_i, ProjVelocity_j, proj_delta_vel, delta_p, delta_rho;
	unsigned short iDim, iVar, jVar, kVar;
  
public:
  
	/*!
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CUpwRoe_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
	/*!
	 * \brief Destructor of the class.
	 */
	~CUpwRoe_Flow(void);
  
	/*!
	 * \brief Compute the Roe's flux between two nodes i and j.
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
	 * \param[in] config - Definition of the particular problem.
	 */
	void ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config);
};

/*!
 * \class CUpwRoeTurkel_Flow
 * \brief Class for solving an approximate Riemann solver of Roe with Turkel Preconditioning for the flow equations.
 * \ingroup ConvDiscr
 * \author A. K. Lonkar (Stanford University)
 * \version 1.0.0
 */
class CUpwRoeTurkel_Flow : public CNumerics {
private:
	bool implicit, grid_movement;
	double *Diff_U;
	double *Velocity_i, *Velocity_j, *RoeVelocity;
	double *Proj_flux_tensor_i, *Proj_flux_tensor_j;
	double *Lambda, *Epsilon;
	double **absPeJac,**invRinvPe,**R_Tensor,**Matrix,**Art_Visc;
	double sq_vel, Proj_ModJac_Tensor_ij, Density_i, Energy_i, SoundSpeed_i, Pressure_i, Enthalpy_i,
	Density_j, Energy_j, SoundSpeed_j, Pressure_j, Enthalpy_j, R, RoePressure, RoeDensity, RoeEnthalpy, RoeSoundSpeed,
	ProjVelocity, ProjVelocity_i, ProjVelocity_j;
	unsigned short iDim, iVar, jVar, kVar;
	double Beta, Beta_min, Beta_max;
  double r_hat, s_hat, t_hat, rhoB2a2, sqr_one_m_Betasqr_Lam1;
	double Beta2, one_m_Betasqr, one_p_Betasqr, sqr_two_Beta_c_Area;
	double local_Mach;
  
public:
  
	/*!
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CUpwRoeTurkel_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
	/*!
	 * \brief Destructor of the class.
	 */
	~CUpwRoeTurkel_Flow(void);
  
	/*!
	 * \brief Compute the Roe's flux between two nodes i and j.
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
	 * \param[in] config - Definition of the particular problem.
	 */
	void ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config);
  
	/*!
	 * \brief Get the Preconditioning Beta.
	 * \return Beta - Value of the low Mach Preconditioner.
	 */
	double GetPrecond_Beta();
};

/*!
 * \class CUpwAUSM_Flow
 * \brief Class for solving an approximate Riemann AUSM.
 * \ingroup ConvDiscr
 * \author F. Palacios
 * \version 1.0.0
 */
class CUpwAUSM_Flow : public CNumerics {
private:
	bool implicit;
	double *Diff_U;
	double *Velocity_i, *Velocity_j, *RoeVelocity;
	double *Proj_flux_tensor_i, *Proj_flux_tensor_j;
	double *delta_wave, *delta_vel;
	double *Lambda, *Epsilon;
	double **P_Tensor, **invP_Tensor;
	double sq_vel, Proj_ModJac_Tensor_ij, Density_i, Energy_i, SoundSpeed_i, Pressure_i, Enthalpy_i,
	Density_j, Energy_j, SoundSpeed_j, Pressure_j, Enthalpy_j, R, RoeDensity, RoeEnthalpy, RoeSoundSpeed,
	ProjVelocity, ProjVelocity_i, ProjVelocity_j, proj_delta_vel, delta_p, delta_rho;
	unsigned short iDim, iVar, jVar, kVar;
  double mL, mR, mLP, mRM, mF, pLP, pRM, pF, Phi;
  
public:
  
	/*!
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CUpwAUSM_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
	/*!
	 * \brief Destructor of the class.
	 */
	~CUpwAUSM_Flow(void);
  
	/*!
	 * \brief Compute the Roe's flux between two nodes i and j.
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
	 * \param[in] config - Definition of the particular problem.
	 */
	void ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config);
};

/*!
 * \class CUpwHLLC_Flow
 * \brief Class for solving an approximate Riemann AUSM.
 * \ingroup ConvDiscr
 * \author F. Palacios, based on the Joe code implementation
 * \version 1.0.0
 */
class CUpwHLLC_Flow : public CNumerics {
private:
	bool implicit;
	double *Diff_U;
	double *Velocity_i, *Velocity_j, *RoeVelocity;
	double *Proj_flux_tensor_i, *Proj_flux_tensor_j;
	double *delta_wave, *delta_vel;
	double *Lambda, *Epsilon;
	double **P_Tensor, **invP_Tensor;
	double sq_vel, sq_vel_i, sq_vel_j, Proj_ModJac_Tensor_ij, Density_i, Energy_i, SoundSpeed_i, Pressure_i, Enthalpy_i,
	Density_j, Energy_j, SoundSpeed_j, Pressure_j, Enthalpy_j, R, RoeDensity, RoeEnthalpy, RoeSoundSpeed,
	ProjVelocity, ProjVelocity_i, ProjVelocity_j, proj_delta_vel, delta_p, delta_rho;
	unsigned short iDim, iVar, jVar, kVar;
  double Rrho, tmp, velRoe[3], uRoe, gamPdivRho, sq_velRoe, cRoe, sL, sR, sM, pStar, invSLmSs, sLmuL, rhoSL, rhouSL[3],
  eSL, invSRmSs, sRmuR, rhoSR, rhouSR[3], eSR;
  
public:
  
	/*!
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CUpwHLLC_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
	/*!
	 * \brief Destructor of the class.
	 */
	~CUpwHLLC_Flow(void);
  
	/*!
	 * \brief Compute the Roe's flux between two nodes i and j.
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
	 * \param[in] config - Definition of the particular problem.
	 */
	void ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config);
};

/*!
 * \class CUpwSca_TurbSA
 * \brief Class for doing a scalar upwind solver for the Spalar-Allmaral turbulence model equations.
 * \ingroup ConvDiscr
 * \author A. Bueno.
 * \version 1.0.0
 */
class CUpwSca_TurbSA : public CNumerics {
private:
	double *Velocity_i, *Velocity_j;
	bool implicit, grid_movement;
	double Density_i, Density_j, q_ij, a0, a1;
	unsigned short iDim;
  
public:
  
	/*!
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CUpwSca_TurbSA(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
	/*!
	 * \brief Destructor of the class.
	 */
	~CUpwSca_TurbSA(void);
  
	/*!
	 * \brief Compute the scalar upwind flux between two nodes i and j.
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
	 * \param[in] config - Definition of the particular problem.
	 */
	void ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config);
};

/*!
 * \class CUpwSca_TurbSST
 * \brief Class for doing a scalar upwind solver for the Menter SST turbulence model equations.
 * \ingroup ConvDiscr
 * \author A. Campos.
 * \version 1.0.0
 */
class CUpwSca_TurbSST : public CNumerics {
private:
	double *Velocity_i, *Velocity_j;
	bool implicit, grid_movement;
	double Density_i, Density_j,
	q_ij,
	a0, a1;
	unsigned short iDim;
  
public:
  
	/*!
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CUpwSca_TurbSST(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
	/*!
	 * \brief Destructor of the class.
	 */
	~CUpwSca_TurbSST(void);
  
	/*!
	 * \brief Compute the scalar upwind flux between two nodes i and j.
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
	 * \param[in] config - Definition of the particular problem.
	 */
	void ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config);
};

/*!
 * \class CCentJST_Flow
 * \brief Class for centered shceme - JST.
 * \ingroup ConvDiscr
 * \author F. Palacios.
 * \version 1.0.0
 */
class CCentJST_Flow : public CNumerics {
  
private:
	unsigned short iDim, iVar, jVar; /*!< \brief Iteration on dimension and variables. */
	double *Diff_U, *Diff_Lapl, /*!< \brief Diference of conservative variables and undivided laplacians. */
	*Velocity_i, *Velocity_j, /*!< \brief Velocity at node 0 and 1. */
	*MeanVelocity, ProjVelocity, ProjVelocity_i, ProjVelocity_j,  /*!< \brief Mean and projected velocities. */
	Density_i, Density_j, Energy_i, Energy_j,  /*!< \brief Mean Density and energies. */
	sq_vel_i, sq_vel_j,   /*!< \brief Modulus of the velocity and the normal vector. */
	MeanDensity, MeanPressure, MeanEnthalpy, MeanEnergy, /*!< \brief Mean values of primitive variables. */
	Param_p, Param_Kappa_2, Param_Kappa_4, /*!< \brief Artificial dissipation parameters. */
	Local_Lambda_i, Local_Lambda_j, MeanLambda, /*!< \brief Local eingenvalues. */
	Phi_i, Phi_j, sc2, sc4, StretchingFactor, /*!< \brief Streching parameters. */
	*Proj_flux_tensor,  /*!< \brief Projected inviscid flux tensor. */
	Epsilon_2, Epsilon_4, cte_0, cte_1, /*!< \brief Artificial dissipation values. */
  ProjGridVel_i, ProjGridVel_j, ProjGridVel;  /*!< \brief Projected grid velocity. */
	bool implicit, /*!< \brief Implicit calculation. */
	grid_movement, /*!< \brief Modification for grid movement. */
	stretching; /*!< \brief Stretching factor. */
  
  
public:
  
	/*!
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimension of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CCentJST_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
	/*!
	 * \brief Destructor of the class.
	 */
	~CCentJST_Flow(void);
  
	/*!
	 * \brief Compute the flow residual using a JST method.
	 * \param[out] val_resconv - Pointer to the convective residual.
	 * \param[out] val_resvisc - Pointer to the artificial viscosity residual.
	 * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
	 * \param[in] config - Definition of the particular problem.
	 */
	void ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j,
                       CConfig *config);
};

/*!
 * \class CCentLax_Flow
 * \brief Class for computing the Lax-Friedrich centered scheme.
 * \ingroup ConvDiscr
 * \author F. Palacios.
 * \version 1.0.0
 */
class CCentLax_Flow : public CNumerics {
private:
	unsigned short iDim, iVar, jVar; /*!< \brief Iteration on dimension and variables. */
	double *Diff_U, /*!< \brief Difference of conservative variables. */
	*Velocity_i, *Velocity_j, /*!< \brief Velocity at node 0 and 1. */
	*MeanVelocity, ProjVelocity, ProjVelocity_i, ProjVelocity_j,  /*!< \brief Mean and projected velocities. */
	*Proj_flux_tensor,  /*!< \brief Projected inviscid flux tensor. */
	Density_i, Density_j, Energy_i, Energy_j,  /*!< \brief Mean Density and energies. */
	sq_vel_i, sq_vel_j,   /*!< \brief Modulus of the velocity and the normal vector. */
	MeanDensity, MeanPressure, MeanEnthalpy, MeanEnergy, /*!< \brief Mean values of primitive variables. */
	Param_p, Param_Kappa_0, /*!< \brief Artificial dissipation parameters. */
	Local_Lambda_i, Local_Lambda_j, MeanLambda, /*!< \brief Local eingenvalues. */
	Phi_i, Phi_j, sc0, StretchingFactor, /*!< \brief Streching parameters. */
	Epsilon_0, cte; /*!< \brief Artificial dissipation values. */
	bool implicit, /*!< \brief Implicit calculation. */
	grid_movement, /*!< \brief Modification for grid movement. */
	stretching, ProjGridVel;
  
public:
  
	/*!
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimension of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CCentLax_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
	/*!
	 * \brief Destructor of the class.
	 */
	~CCentLax_Flow(void);
  
	/*!
	 * \brief Compute the flow residual using a Lax method.
	 * \param[out] val_resconv - Pointer to the convective residual.
	 * \param[out] val_resvisc - Pointer to the artificial viscosity residual.
	 * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
	 * \param[in] config - Definition of the particular problem.
	 */
	void ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j,
                       CConfig *config);
};

/*!
 * \class CAvgGrad_Flow
 * \brief Class for computing viscous term using the average of gradients.
 * \ingroup ViscDiscr
 * \author A. Bueno, and F. Palacios.
 * \version 1.0.0
 */
class CAvgGrad_Flow : public CNumerics {
private:
	unsigned short iDim, iVar, jVar;		/*!< \brief Iterators in dimension an variable. */
	double *Mean_PrimVar,					/*!< \brief Mean primitive variables. */
	*PrimVar_i, *PrimVar_j,				/*!< \brief Primitives variables at point i and 1. */
	**Mean_GradPrimVar,						/*!< \brief Mean value of the gradient. */
	Mean_Laminar_Viscosity, Mean_Eddy_Viscosity, /*!< \brief Mean value of the viscosity. */
	Mean_turb_ke,				/*!< \brief Mean value of the turbulent kinetic energy. */
	*Proj_flux_tensor,	/*!< \brief Projection of the viscous fluxes. */
	dist_ij;						/*!< \brief Length of the edge and face. */
	bool implicit; /*!< \brief Implicit calculus. */
  
public:
  
	/*!
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimension of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CAvgGrad_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
	/*!
	 * \brief Destructor of the class.
	 */
	~CAvgGrad_Flow(void);
  
	/*!
	 * \brief Compute the viscous flow residual using an average of gradients.
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
	 * \param[in] config - Definition of the particular problem.
	 */
	void ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config);
};

/*!
 * \class CAvgGrad_TurbSA
 * \brief Class for computing viscous term using average of gradients (Spalart-Allmaras Turbulence model).
 * \ingroup ViscDiscr
 * \author A. Bueno.
 * \version 1.0.0
 */
class CAvgGrad_TurbSA : public CNumerics {
private:
	double **Mean_GradTurbVar;
	double *Proj_Mean_GradTurbVar_Kappa, *Proj_Mean_GradTurbVar_Edge;
	double *Edge_Vector;
	bool implicit;
	double sigma;
	double nu_i, nu_j, nu_e;
	double dist_ij_2;
	double proj_vector_ij;
	unsigned short iVar, iDim;
	double nu_hat_i;
	double nu_hat_j;
  
public:
  
	/*!
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CAvgGrad_TurbSA(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
	/*!
	 * \brief Destructor of the class.
	 */
	~CAvgGrad_TurbSA(void);
  
	/*!
	 * \brief Compute the viscous turbulence terms residual using an average of gradients.
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[out] Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[out] Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
	 * \param[in] config - Definition of the particular problem.
	 */
	void ComputeResidual(double *val_residual, double **Jacobian_i, double **Jacobian_j, CConfig *config);
};

/*!
 * \class CAvgGradCorrected_Flow
 * \brief Class for computing viscous term using the average of gradients with a correction.
 * \ingroup ViscDiscr
 * \author A. Bueno, and F. Palacios.
 * \version 1.0.0
 */
class CAvgGradCorrected_Flow : public CNumerics {
private:
	unsigned short iDim, iVar, jVar;		/*!< \brief Iterators in dimension an variable. */
	double *Mean_PrimVar,					/*!< \brief Mean primitive variables. */
	*PrimVar_i, *PrimVar_j,				/*!< \brief Primitives variables at point i and 1. */
	*Edge_Vector,									/*!< \brief Vector form point i to point j. */
	**Mean_GradPrimVar, *Proj_Mean_GradPrimVar_Edge,	/*!< \brief Mean value of the gradient. */
	Mean_Laminar_Viscosity, Mean_Eddy_Viscosity,			/*!< \brief Mean value of the viscosity. */
	Mean_turb_ke,				/*!< \brief Mean value of the turbulent kinetic energy. */
	dist_ij_2,					/*!< \brief Length of the edge and face. */
	*Proj_flux_tensor;	/*!< \brief Projection of the viscous fluxes. */
	bool implicit;			/*!< \brief Implicit calculus. */
  
public:
  
	/*!
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimension of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CAvgGradCorrected_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
	/*!
	 * \brief Destructor of the class.
	 */
	~CAvgGradCorrected_Flow(void);
  
	/*!
	 * \brief Compute the viscous flow residual using an average of gradients with correction.
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
	 * \param[in] config - Definition of the particular problem.
	 */
	void ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config);
};

/*!
 * \class CAvgGradCorrected_TurbSA
 * \brief Class for computing viscous term using average of gradients with correction (Spalart-Allmaras turbulence model).
 * \ingroup ViscDiscr
 * \author A. Bueno.
 * \version 1.0.0
 */
class CAvgGradCorrected_TurbSA : public CNumerics {
private:
	double **Mean_GradTurbVar;
	double *Proj_Mean_GradTurbVar_Kappa, *Proj_Mean_GradTurbVar_Edge, *Proj_Mean_GradTurbVar_Corrected;
	double *Edge_Vector;
	bool implicit;
	double sigma, nu_i, nu_j, nu_e, dist_ij_2, proj_vector_ij, nu_hat_i, nu_hat_j;
	unsigned short iVar, iDim;
  
public:
  
	/*!
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CAvgGradCorrected_TurbSA(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
	/*!
	 * \brief Destructor of the class.
	 */
	~CAvgGradCorrected_TurbSA(void);
  
	/*!
	 * \brief Compute the viscous turbulent residual using an average of gradients with correction.
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[out] Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[out] Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
	 * \param[in] config - Definition of the particular problem.
	 */
	void ComputeResidual(double *val_residual, double **Jacobian_i, double **Jacobian_j, CConfig *config);
};

/*!
 * \class CAvgGrad_TurbSST
 * \brief Class for computing viscous term using average of gradient with correction (Menter SST turbulence model).
 * \ingroup ViscDiscr
 * \author A. Bueno.
 * \version 1.0.0
 */
class CAvgGrad_TurbSST : public CNumerics {
private:
	double sigma_k1,                     /*!< \brief Constants for the viscous terms, k-w (1), k-eps (2)*/
	sigma_k2,
	sigma_om1,
	sigma_om2;
  
	double diff_kine,                     /*!< \brief Diffusivity for viscous terms of tke eq */
	diff_omega;                           /*!< \brief Diffusivity for viscous terms of omega eq */
  
	double *Edge_Vector,                  /*!< \brief Vector from node i to node j. */
	dist_ij_2,                            /*!< \brief |Edge_Vector|^2 */
	proj_vector_ij;                       /*!< \brief (Edge_Vector DOT normal)/|Edge_Vector|^2 */
  
	double **Mean_GradTurbVar,            /*!< \brief Average of gradients at cell face */
	*Proj_Mean_GradTurbVar_Normal,        /*!< \brief Mean_gradTurbVar DOT normal */
	*Proj_Mean_GradTurbVar_Edge,          /*!< \brief Mean_gradTurbVar DOT Edge_Vector */
	*Proj_Mean_GradTurbVar_Corrected;
  
	double F1_i, F1_j;                    /*!< \brief Menter's first blending function */
  
	bool implicit;
	unsigned short iVar, iDim;
  
public:
  
	/*!
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CAvgGrad_TurbSST(unsigned short val_nDim, unsigned short val_nVar, double* constants, CConfig *config);
  
	/*!
	 * \brief Destructor of the class.
	 */
	~CAvgGrad_TurbSST(void);
  
	/*!
	 * \brief Sets value of first blending function.
	 */
	void SetF1blending(double val_F1_i, double val_F1_j) { F1_i = val_F1_i; F1_j = val_F1_j;}
  
	/*!
	 * \brief Compute the viscous turbulent residual using an average of gradients wtih correction.
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[out] Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[out] Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
	 * \param[in] config - Definition of the particular problem.
	 */
	void ComputeResidual(double *val_residual, double **Jacobian_i, double **Jacobian_j, CConfig *config);
  
};

/*!
 * \class CAvgGradCorrected_TurbSST
 * \brief Class for computing viscous term using average of gradient with correction (Menter SST turbulence model).
 * \ingroup ViscDiscr
 * \author A. Bueno.
 * \version 1.0.0
 */
class CAvgGradCorrected_TurbSST : public CNumerics {
private:
	double sigma_k1,                     /*!< \brief Constants for the viscous terms, k-w (1), k-eps (2)*/
	sigma_k2,
	sigma_om1,
	sigma_om2;
  
	double diff_kine,                     /*!< \brief Diffusivity for viscous terms of tke eq */
	diff_omega;                           /*!< \brief Diffusivity for viscous terms of omega eq */
  
	double *Edge_Vector,                  /*!< \brief Vector from node i to node j. */
	dist_ij_2,                            /*!< \brief |Edge_Vector|^2 */
	proj_vector_ij;                       /*!< \brief (Edge_Vector DOT normal)/|Edge_Vector|^2 */
  
	double **Mean_GradTurbVar,            /*!< \brief Average of gradients at cell face */
	*Proj_Mean_GradTurbVar_Normal,        /*!< \brief Mean_gradTurbVar DOT normal */
	*Proj_Mean_GradTurbVar_Edge,          /*!< \brief Mean_gradTurbVar DOT Edge_Vector */
	*Proj_Mean_GradTurbVar_Corrected;
  
	double F1_i, F1_j;                    /*!< \brief Menter's first blending function */
  
	bool implicit;
	unsigned short iVar, iDim;
  
public:
  
	/*!
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CAvgGradCorrected_TurbSST(unsigned short val_nDim, unsigned short val_nVar, double* constants, CConfig *config);
  
	/*!
	 * \brief Destructor of the class.
	 */
	~CAvgGradCorrected_TurbSST(void);
  
	/*!
	 * \brief Sets value of first blending function.
	 */
	void SetF1blending(double val_F1_i, double val_F1_j) { F1_i = val_F1_i; F1_j = val_F1_j;}
  
	/*!
	 * \brief Compute the viscous turbulent residual using an average of gradients wtih correction.
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[out] Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[out] Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
	 * \param[in] config - Definition of the particular problem.
	 */
	void ComputeResidual(double *val_residual, double **Jacobian_i, double **Jacobian_j, CConfig *config);
  
};

/*!
 * \class CSourceNothing
 * \brief Dummy class.
 * \ingroup SourceDiscr
 * \author F. Palacios.
 * \version 1.0.0
 */
class CSourceNothing : public CNumerics {
public:
  
	/*!
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CSourceNothing(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
	/*!
	 * \brief Destructor of the class.
	 */
	~CSourceNothing(void);
};

/*!
 * \class CSourcePieceWise_TurbSA
 * \brief Class for integrating the source terms of the Spalart-Allmaras turbulence model equation.
 * \ingroup SourceDiscr
 * \author A. Bueno.
 * \version 1.0.0
 */
class CSourcePieceWise_TurbSA : public CNumerics {
private:
	double cv1_3;
	double k2;
	double cb1;
	double cw2;
	double cw3_6;
  double cb2_sigma;
	double sigma;
	double cb2;
	double cw1;
	double DivVelocity, Vorticity;
	unsigned short iDim;
	double nu, Ji, fv1, fv2, Omega, S, Shat, inv_Shat, dist_i_2, Ji_2, Ji_3, inv_k2_d2;
	double r, g, g_6, glim, fw;
	double norm2_Grad;
	double dfv1, dfv2, dShat;
	double dr, dg, dfw;;
	double nu_hat_i;
	double grad_nu_hat;
	double prod_grads;
  bool transition;
  bool rotating_frame;
  double div, StrainMag;
  double beta, gamma_sep, gamma_eff, intermittency;
  double Freattach, r_t, s1;
  double Production, Destruction, CrossProduction;

public:
  
	/*!
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CSourcePieceWise_TurbSA(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
	/*!
	 * \brief Destructor of the class.
	 */
	~CSourcePieceWise_TurbSA(void);
  
	/*!
	 * \brief Residual for source term integration.
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
	 * \param[in] config - Definition of the particular problem.
	 */
	void ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config);
  
	/*!
	 * \brief Residual for source term integration.
	 * \param[in] intermittency_in - Value of the intermittency.
	 */
  void SetIntermittency(double intermittency_in);
  
  /*!
	 * \brief Residual for source term integration.
	 * \param[in] val_production - Value of the Production.
	 */
  void SetProduction(double val_production);
  
  /*!
	 * \brief Residual for source term integration.
	 * \param[in] val_destruction - Value of the Destruction.
	 */
  void SetDestruction(double val_destruction);
  
  /*!
	 * \brief Residual for source term integration.
	 * \param[in] val_crossproduction - Value of the CrossProduction.
	 */
  void SetCrossProduction(double val_crossproduction);
  
  /*!
	 * \brief Residual for source term integration.
	 * \param[in] val_production - Value of the Production.
	 */
  double GetProduction(void);
  
  /*!
	 * \brief Residual for source term integration.
	 * \param[in] val_destruction - Value of the Destruction.
	 */
  double GetDestruction(void);
  
  /*!
	 * \brief Residual for source term integration.
	 * \param[in] val_crossproduction - Value of the CrossProduction.
	 */
  double GetCrossProduction(void);
};

/*!
 * \class CSourcePieceWise_TurbSST
 * \brief Class for integrating the source terms of the Menter SST turbulence model equations.
 * \ingroup SourceDiscr
 * \author A. Campos.
 * \version 1.0.0
 */
class CSourcePieceWise_TurbSST : public CNumerics {
private:
	double F1_i,
	F1_j,
	F2_i,
	F2_j;
  
	double alfa_1,
	alfa_2,
	beta_1,
	beta_2,
	sigma_omega_1,
	sigma_omega_2,
	beta_star,
	a1;
  
	double StrainMag,
	CDkw,
	norm2_Grad;
    
public:
  
	/*!
	 * \brief Constructor of the class.
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nVar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CSourcePieceWise_TurbSST(unsigned short val_nDim, unsigned short val_nVar, double* constants, CConfig *config);
  
	/*!
	 * \brief Destructor of the class.
	 */
	~CSourcePieceWise_TurbSST(void);
  
	/*!
	 * \brief Set the value of the first blending function.
	 * \param[in] val_F1_i - Value of the first blending function at point i.
	 * \param[in] val_F1_j - Value of the first blending function at point j.
	 */
	void SetF1blending(double val_F1_i, double val_F1_j);
  
	/*!
	 * \brief Set the value of the second blending function.
	 * \param[in] val_F1_i - Value of the second blending function at point i.
	 * \param[in] val_F1_j - Value of the second blending function at point j.
	 */
	void SetF2blending(double val_F2_i, double val_F2_j);
  
	/*!
	 * \brief Set the value of the rate of strain magnitude.
	 * \param[in] val_StrainMag_i - Value of the magnitude of rate of strain at point i.
	 * \param[in] val_StrainMag_j - Value of the magnitude of rate of strain at point j.
	 */
	virtual void SetStrainMag(double val_StrainMag_i, double val_StrainMag_j);
  
	/*!
	 * \brief Set the value of the cross diffusion for the SST model.
	 * \param[in] val_CDkw_i - Value of the cross diffusion at point i.
	 * \param[in] val_CDkw_j - Value of the cross diffusion at point j.
	 */
	virtual void SetCrossDiff(double val_CDkw_i, double val_CDkw_j);
  
	/*!
	 * \brief Residual for source term integration.
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
	 * \param[in] config - Definition of the particular problem.
	 */
	void ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config);
};

#include "numerics_structure.inl"
