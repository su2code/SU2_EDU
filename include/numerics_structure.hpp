/*!
 * \file numerics_structure.hpp
 * \brief Headers of the main subroutines for the dumerical definition of the problem.
 *        The subroutines and functions are in the <i>numerics_structure.cpp</i>,
 *        <i>numerics_convective.cpp</i>, <i>numerics_viscous.cpp</i>, and
 *        <i>numerics_source.cpp</i> files.
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
 * \version 1.1.0
 */
class CNumerics {
  
public:
    
	/*!
	 * \brief Constructor of the class.
	 */
	CNumerics(void);
  
	/*!
	 * \brief Destructor of the class.
	 */
	virtual ~CNumerics(void);
  
  /*!
	 * \brief A pure virtual member.
	 */
  virtual void SetLambda(double val_lambda_i, double val_lambda_j) {/* empty */};

  /*!
	 * \brief A pure virtual member.
	 */
  virtual void SetTurbKineticEnergy(double val_turb_ke_i, double val_turb_ke_j) {/* empty */};
  
  /*!
	 * \brief A pure virtual member.
	 */
  virtual	void SetDistance(double val_dist_i, double val_dist_j) {/* empty */};
  
  /*!
	 * \brief A pure virtual member.
	 */
  virtual	void SetUndivided_Laplacian(double *val_und_lapl_i, double *val_und_lapl_j) {/* empty */};
  
  /*!
	 * \brief A pure virtual member.
	 */
  virtual	void SetSensor(double val_sensor_i, double val_sensor_j) {/* empty */};

  /*!
	 * \brief A pure virtual member.
	 */
  virtual	void SetPrimitive(double *val_v_i, double *val_v_j) {/* empty */};
  
  /*!
	 * \brief A pure virtual member.
	 */
  virtual	void SetNormal(double *val_normal) {/* empty */};

  /*!
	 * \brief A pure virtual member.
	 */
  virtual	void SetCoord(double *val_coord_i, double *val_coord_j) {/* empty */};
  
  /*!
	 * \brief A pure virtual member.
	 */
  virtual	void SetNeighbor(unsigned short val_neighbor_i, unsigned short val_neighbor_j) {/* empty */};
  
  /*!
	 * \brief A pure virtual member.
	 */
  virtual	void SetTurbVar(double *val_turbvar_i, double *val_turbvar_j) {/* empty */};

  /*!
	 * \brief A pure virtual member.
	 */
  virtual	void SetTurbVarGradient(double **val_turbvar_grad_i, double **val_turbvar_grad_j) {/* empty */};
  
  /*!
	 * \brief A pure virtual member.
	 */
  virtual	void SetPrimVarGradient(double **val_primvar_grad_i, double **val_primvar_grad_j) {/* empty */};
  
  /*!
	 * \brief A pure virtual member.
	 */
  virtual	void SetVolume(double val_volume) {/* empty */};
  
  /*!
	 * \brief A pure virtual member.
	 */
	virtual void SetF1blending(double val_F1_i, double val_F1_j) {/* empty */};
  
  /*!
	 * \brief A pure virtual member.
	 */
  virtual void SetF2blending(double val_F2_i, double val_F2_j) {/* empty */};

  /*!
	 * \brief A pure virtual member.
	 */
	virtual void SetStrainMag(double val_StrainMag_i, double val_StrainMag_j) {/* empty */};

  /*!
	 * \brief A pure virtual member.
	 */
	virtual void SetCrossDiff(double val_CDkw_i, double val_CDkw_j){/* empty */};

  /*!
	 * \brief A pure virtual member.
	 */
	virtual void ComputeResidual(double *val_residual, double **val_Jacobian_i,
                               double **val_Jacobian_j, CConfig *config) {/* empty */};
  
  /*!
	 * \brief A pure virtual member.
	 */
  virtual double GetPrecond_Beta(void) { return 0; };
  
};

/*!
 * \class CUpwRoe_Flow
 * \brief Class for solving an approximate Riemann solver of Roe for the flow equations.
 * \ingroup ConvDiscr
 * \author F. Palacios
 * \version 1.1.0
 */
class CUpwRoe_Flow : public CNumerics {
private:
  unsigned short nDim, nVar;	/*!< \brief Number of dimensions and variables. */
	double Gamma;				/*!< \brief Fluid's Gamma constant (ratio of specific heats). */
	double Gamma_Minus_One;		/*!< \brief Fluids's Gamma - 1.0  . */
	double Gas_Constant;		 		/*!< \brief Gas constant. */
  double Area;				/*!< \brief Area of the face i-j. */
  double *Normal;	/*!< \brief Normal vector, it norm is the area of the face. */
  double *UnitNormal;		/*!< \brief Unitary normal vector. */
  double *V_i,		/*!< \brief Vector of primitive variables at point i. */
	*V_j;				/*!< \brief Vector of primitive variables at point j. */
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
	 * \brief Set the value of the normal vector to the face between two points.
	 * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
	 */
	void SetNormal(double *val_normal);
  
  /*!
	 * \brief Set the value of the primitive variables.
	 * \param[in] val_v_i - Value of the primitive variable at point i.
	 * \param[in] val_v_j - Value of the primitive variable at point j.
	 */
	void SetPrimitive(double *val_v_i, double *val_v_j);
  
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
 * \version 1.1.0
 */
class CUpwRoeTurkel_Flow : public CNumerics {
private:
  unsigned short nDim, nVar;	/*!< \brief Number of dimensions and variables. */
	double Gamma;				/*!< \brief Fluid's Gamma constant (ratio of specific heats). */
	double Gamma_Minus_One;		/*!< \brief Fluids's Gamma - 1.0  . */
	double Gas_Constant;		 		/*!< \brief Gas constant. */
  double Area;				/*!< \brief Area of the face i-j. */
  double *UnitNormal;		/*!< \brief Unitary normal vector. */
  double *Normal;	/*!< \brief Normal vector, it norm is the area of the face. */
  double *V_i,		/*!< \brief Vector of primitive variables at point i. */
	*V_j;				/*!< \brief Vector of primitive variables at point j. */
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
	 * \brief Set the value of the normal vector to the face between two points.
	 * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
	 */
	void SetNormal(double *val_normal);
  
  /*!
	 * \brief Set the value of the primitive variables.
	 * \param[in] val_v_i - Value of the primitive variable at point i.
	 * \param[in] val_v_j - Value of the primitive variable at point j.
	 */
	void SetPrimitive(double *val_v_i, double *val_v_j);
  
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
	double GetPrecond_Beta(void);
};

/*!
 * \class CUpwAUSM_Flow
 * \brief Class for solving an approximate Riemann AUSM.
 * \ingroup ConvDiscr
 * \author F. Palacios
 * \version 1.1.0
 */
class CUpwAUSM_Flow : public CNumerics {
private:
  unsigned short nDim, nVar;	/*!< \brief Number of dimensions and variables. */
	double Gamma;				/*!< \brief Fluid's Gamma constant (ratio of specific heats). */
	double Gamma_Minus_One;		/*!< \brief Fluids's Gamma - 1.0  . */
	double Gas_Constant;		 		/*!< \brief Gas constant. */
  double Area;				/*!< \brief Area of the face i-j. */
  double *Normal;	/*!< \brief Normal vector, it norm is the area of the face. */
  double *UnitNormal;		/*!< \brief Unitary normal vector. */
  double *V_i,		/*!< \brief Vector of primitive variables at point i. */
	*V_j;				/*!< \brief Vector of primitive variables at point j. */
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
	 * \brief Set the value of the normal vector to the face between two points.
	 * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
	 */
	void SetNormal(double *val_normal);
  
  /*!
	 * \brief Set the value of the primitive variables.
	 * \param[in] val_v_i - Value of the primitive variable at point i.
	 * \param[in] val_v_j - Value of the primitive variable at point j.
	 */
	void SetPrimitive(double *val_v_i, double *val_v_j);
  
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
 * \version 1.1.0
 */
class CUpwHLLC_Flow : public CNumerics {
private:
  unsigned short nDim, nVar;	/*!< \brief Number of dimensions and variables. */
	double Gamma;				/*!< \brief Fluid's Gamma constant (ratio of specific heats). */
	double Gamma_Minus_One;		/*!< \brief Fluids's Gamma - 1.0  . */
	double Gas_Constant;		 		/*!< \brief Gas constant. */
  double Area;				/*!< \brief Area of the face i-j. */
  double *Normal;	/*!< \brief Normal vector, it norm is the area of the face. */
  double *UnitNormal;		/*!< \brief Unitary normal vector. */
  double *V_i,		/*!< \brief Vector of primitive variables at point i. */
	*V_j;				/*!< \brief Vector of primitive variables at point j. */
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
	 * \brief Set the value of the normal vector to the face between two points.
	 * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
	 */
	void SetNormal(double *val_normal);
  
  /*!
	 * \brief Set the value of the primitive variables.
	 * \param[in] val_v_i - Value of the primitive variable at point i.
	 * \param[in] val_v_j - Value of the primitive variable at point j.
	 */
	void SetPrimitive(double *val_v_i, double *val_v_j);
  
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
 * \author F. Palacios.
 * \version 1.1.0
 */
class CUpwSca_TurbSA : public CNumerics {
private:
  unsigned short nDim, nVar;	/*!< \brief Number of dimensions and variables. */
	double Gamma;				/*!< \brief Fluid's Gamma constant (ratio of specific heats). */
	double Gamma_Minus_One;		/*!< \brief Fluids's Gamma - 1.0  . */
	double Gas_Constant;		 		/*!< \brief Gas constant. */
  double *TurbVar_i,	/*!< \brief Vector of turbulent variables at point i. */
	*TurbVar_j;	/*!< \brief Vector of derivative of turbulent variables at point j. */
  double *Normal;	/*!< \brief Normal vector, it norm is the area of the face. */
  double *V_i,		/*!< \brief Vector of primitive variables at point i. */
	*V_j;				/*!< \brief Vector of primitive variables at point j. */
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
	 * \brief Set the value of the turbulent variable.
	 * \param[in] val_turbvar_i - Value of the turbulent variable at point i.
	 * \param[in] val_turbvar_j - Value of the turbulent variable at point j.
	 */
	void SetTurbVar(double *val_turbvar_i, double *val_turbvar_j);
  
	/*!
	 * \brief Set the value of the normal vector to the face between two points.
	 * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
	 */
	void SetNormal(double *val_normal);
  
  /*!
	 * \brief Set the value of the primitive variables.
	 * \param[in] val_v_i - Value of the primitive variable at point i.
	 * \param[in] val_v_j - Value of the primitive variable at point j.
	 */
	void SetPrimitive(double *val_v_i, double *val_v_j);
  
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
 * \version 1.1.0
 */
class CUpwSca_TurbSST : public CNumerics {
private:
  unsigned short nDim, nVar;	/*!< \brief Number of dimensions and variables. */
	double Gamma;				/*!< \brief Fluid's Gamma constant (ratio of specific heats). */
	double Gamma_Minus_One;		/*!< \brief Fluids's Gamma - 1.0  . */
	double Gas_Constant;		 		/*!< \brief Gas constant. */
  double *TurbVar_i,	/*!< \brief Vector of turbulent variables at point i. */
	*TurbVar_j;	/*!< \brief Vector of derivative of turbulent variables at point j. */
  double *Normal;	/*!< \brief Normal vector, it norm is the area of the face. */
  double *V_i,		/*!< \brief Vector of primitive variables at point i. */
	*V_j;				/*!< \brief Vector of primitive variables at point j. */
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
	 * \brief Set the value of the turbulent variable.
	 * \param[in] val_turbvar_i - Value of the turbulent variable at point i.
	 * \param[in] val_turbvar_j - Value of the turbulent variable at point j.
	 */
	void SetTurbVar(double *val_turbvar_i, double *val_turbvar_j);
  
  /*!
	 * \brief Set the value of the normal vector to the face between two points.
	 * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
	 */
	void SetNormal(double *val_normal);

  /*!
	 * \brief Set the value of the primitive variables.
	 * \param[in] val_v_i - Value of the primitive variable at point i.
	 * \param[in] val_v_j - Value of the primitive variable at point j.
	 */
	void SetPrimitive(double *val_v_i, double *val_v_j);
  
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
 * \version 1.1.0
 */
class CCentJST_Flow : public CNumerics {
  
private:
  unsigned short nDim, nVar;	/*!< \brief Number of dimensions and variables. */
	double Gamma;				/*!< \brief Fluid's Gamma constant (ratio of specific heats). */
	double Gamma_Minus_One;		/*!< \brief Fluids's Gamma - 1.0  . */
	double Gas_Constant;		 		/*!< \brief Gas constant. */
  double Area;				/*!< \brief Area of the face i-j. */
  unsigned short Neighbor_i,	/*!< \brief Number of neighbors of the point i. */
	Neighbor_j;					/*!< \brief Number of neighbors of the point j. */
  double *Normal;	/*!< \brief Normal vector, it norm is the area of the face. */
  double *V_i,		/*!< \brief Vector of primitive variables at point i. */
	*V_j;				/*!< \brief Vector of primitive variables at point j. */
  double Sensor_i,	/*!< \brief Pressure sensor at point i. */
	Sensor_j;			/*!< \brief Pressure sensor at point j. */
  double *Und_Lapl_i, /*!< \brief Undivided laplacians at point i. */
	*Und_Lapl_j;		/*!< \brief Undivided laplacians at point j. */
  double Lambda_i,	/*!< \brief Spectral radius at point i. */
	Lambda_j;			/*!< \brief Spectral radius at point j. */
  double Pressure_i,	/*!< \brief Pressure at point i. */
	Pressure_j;			/*!< \brief Pressure at point j. */
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
  ProjGridVel;  /*!< \brief Projected grid velocity. */
	bool implicit, /*!< \brief Implicit calculation. */
	grid_movement, /*!< \brief Modification for grid movement. */
	stretching; /*!< \brief Stretching factor. */
  double SoundSpeed_i,	/*!< \brief Sound speed at point i. */
	SoundSpeed_j;			/*!< \brief Sound speed at point j. */
	double Enthalpy_i,	/*!< \brief Enthalpy at point i. */
	Enthalpy_j;			/*!< \brief Enthalpy at point j. */
  
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
	 * \brief Set the value of the spectral radius.
	 * \param[in] val_lambda_i - Value of the spectral radius at point i.
	 * \param[in] val_lambda_j - Value of the spectral radius at point j.
	 */
	void SetLambda(double val_lambda_i, double val_lambda_j);
  
  /*!
	 * \brief Set the value of the pressure sensor.
	 * \param[in] val_sensor_i Pressure sensor at point i.
	 * \param[in] val_sensor_j Pressure sensor at point j.
	 */
	void SetSensor(double val_sensor_i, double val_sensor_j);
  
  /*!
	 * \brief Set the value of the primitive variables.
	 * \param[in] val_v_i - Value of the primitive variable at point i.
	 * \param[in] val_v_j - Value of the primitive variable at point j.
	 */
	void SetPrimitive(double *val_v_i, double *val_v_j);
  
  /*!
	 * \brief Set the value of undivided laplacian.
	 * \param[in] val_und_lapl_i Undivided laplacian at point i.
	 * \param[in] val_und_lapl_j Undivided laplacian at point j.
	 */
	void SetUndivided_Laplacian(double *val_und_lapl_i, double *val_und_lapl_j);
  
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
 * \version 1.1.0
 */
class CCentLax_Flow : public CNumerics {
private:
  unsigned short nDim, nVar;	/*!< \brief Number of dimensions and variables. */
	double Gamma;				/*!< \brief Fluid's Gamma constant (ratio of specific heats). */
	double Gamma_Minus_One;		/*!< \brief Fluids's Gamma - 1.0  . */
	double Gas_Constant;		 		/*!< \brief Gas constant. */
  double Area;				/*!< \brief Area of the face i-j. */
  unsigned short Neighbor_i,	/*!< \brief Number of neighbors of the point i. */
	Neighbor_j;					/*!< \brief Number of neighbors of the point j. */
  double *Normal;	/*!< \brief Normal vector, it norm is the area of the face. */
  double *V_i,		/*!< \brief Vector of primitive variables at point i. */
	*V_j;				/*!< \brief Vector of primitive variables at point j. */
  double Lambda_i,	/*!< \brief Spectral radius at point i. */
	Lambda_j;			/*!< \brief Spectral radius at point j. */
  double Pressure_i,	/*!< \brief Pressure at point i. */
	Pressure_j;			/*!< \brief Pressure at point j. */
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
  double SoundSpeed_i,	/*!< \brief Sound speed at point i. */
	SoundSpeed_j;			/*!< \brief Sound speed at point j. */
	double Enthalpy_i,	/*!< \brief Enthalpy at point i. */
	Enthalpy_j;			/*!< \brief Enthalpy at point j. */
  
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
	 * \brief Set the value of the primitive variables.
	 * \param[in] val_v_i - Value of the primitive variable at point i.
	 * \param[in] val_v_j - Value of the primitive variable at point j.
	 */
	void SetPrimitive(double *val_v_i, double *val_v_j);
  
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
	 * \brief Set the value of the spectral radius.
	 * \param[in] val_lambda_i - Value of the spectral radius at point i.
	 * \param[in] val_lambda_j - Value of the spectral radius at point j.
	 */
	void SetLambda(double val_lambda_i, double val_lambda_j);
  
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
 * \author F. Palacios.
 * \version 1.1.0
 */
class CAvgGrad_Flow : public CNumerics {
private:
  unsigned short nDim, nVar;	/*!< \brief Number of dimensions and variables. */
	double Gamma;				/*!< \brief Fluid's Gamma constant (ratio of specific heats). */
	double Gamma_Minus_One;		/*!< \brief Fluids's Gamma - 1.0  . */
	double Gas_Constant;		 		/*!< \brief Gas constant. */
  double Area;				/*!< \brief Area of the face i-j. */
  double **PrimVar_Grad_i,	/*!< \brief Gradient of primitive variables at point i. */
	**PrimVar_Grad_j;			/*!< \brief Gradient of primitive variables at point j. */
  double *Coord_i,	/*!< \brief Cartesians coordinates of point i. */
	*Coord_j;			/*!< \brief Cartesians coordinates of point j. */
  double *Normal;	/*!< \brief Normal vector, it norm is the area of the face. */
  double *UnitNormal;		/*!< \brief Unitary normal vector. */
  double *V_i,		/*!< \brief Vector of primitive variables at point i. */
	*V_j;				/*!< \brief Vector of primitive variables at point j. */
  double turb_ke_i,	/*!< \brief Turbulent kinetic energy at point i. */
	turb_ke_j;			/*!< \brief Turbulent kinetic energy at point j. */
  double Eddy_Viscosity_i,	/*!< \brief Eddy viscosity at point i. */
	Eddy_Viscosity_j;			/*!< \brief Eddy viscosity at point j. */
  double Laminar_Viscosity_i,	/*!< \brief Laminar viscosity at point i. */
	Laminar_Viscosity_j;		/*!< \brief Variation of laminar viscosity at point j. */
  double
  **Flux_Tensor,	/*!< \brief Flux tensor (used for viscous and inviscid purposes. */
	*Proj_Flux_Tensor;		/*!< \brief Flux tensor projected in a direction. */
  double
  **tau,		/*!< \brief Viscous stress tensor. */
	**delta;			/*!< \brief Identity matrix. */
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
	void GetViscousProjJac(double *val_Mean_PrimVar,
                         double val_laminar_viscosity,
                         double val_eddy_viscosity,
                         double val_dist_ij,
                         double *val_normal, double val_dS,
                         double *val_Proj_Visc_Flux,
                         double **val_Proj_Jac_Tensor_i,
                         double **val_Proj_Jac_Tensor_j);
  
  /*!
	 * \brief Set the gradient of the primitive variables.
	 * \param[in] val_primvar_grad_i - Gradient of the primitive variable at point i.
	 * \param[in] val_primvar_grad_j - Gradient of the primitive variable at point j.
	 */
	void SetPrimVarGradient(double **val_primvar_grad_i,
                          double **val_primvar_grad_j);
  
  /*!
	 * \brief Set coordinates of the points.
	 * \param[in] val_coord_i - Coordinates of the point i.
	 * \param[in] val_coord_j - Coordinates of the point j.
	 */
	void SetCoord(double *val_coord_i, double *val_coord_j);
  
  /*!
	 * \brief Set the value of the normal vector to the face between two points.
	 * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
	 */
	void SetNormal(double *val_normal);
  
  /*!
	 * \brief Set the eddy viscosity.
	 * \param[in] val_eddy_viscosity_i - Value of the eddy viscosity at point i.
	 * \param[in] val_eddy_viscosity_j - Value of the eddy viscosity at point j.
	 */
	void SetEddyViscosity(double val_eddy_viscosity_i,
                        double val_eddy_viscosity_j);
  
  /*!
	 * \brief Set the value of the primitive variables.
	 * \param[in] val_v_i - Value of the primitive variable at point i.
	 * \param[in] val_v_j - Value of the primitive variable at point j.
	 */
	void SetPrimitive(double *val_v_i, double *val_v_j);
  
  /*!
	 * \brief Set the turbulent kinetic energy.
	 * \param[in] val_turb_ke_i - Value of the turbulent kinetic energy at point i.
	 * \param[in] val_turb_ke_j - Value of the turbulent kinetic energy at point j.
	 */
	void SetTurbKineticEnergy(double val_turb_ke_i, double val_turb_ke_j);
  
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
	 * \brief Set the laminar viscosity.
	 * \param[in] val_laminar_viscosity_i - Value of the laminar viscosity at point i.
	 * \param[in] val_laminar_viscosity_j - Value of the laminar viscosity at point j.
	 */
	void SetLaminarViscosity(double val_laminar_viscosity_i,
                           double val_laminar_viscosity_j);
  
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
 * \author F. Palacios.
 * \version 1.1.0
 */
class CAvgGrad_TurbSA : public CNumerics {
private:
  unsigned short nDim, nVar;	/*!< \brief Number of dimensions and variables. */
	double Gamma;				/*!< \brief Fluid's Gamma constant (ratio of specific heats). */
	double Gamma_Minus_One;		/*!< \brief Fluids's Gamma - 1.0  . */
	double Gas_Constant;		 		/*!< \brief Gas constant. */
  double **TurbVar_Grad_i,	/*!< \brief Gradient of turbulent variables at point i. */
	**TurbVar_Grad_j;			/*!< \brief Gradient of turbulent variables at point j. */
  double *Coord_i,	/*!< \brief Cartesians coordinates of point i. */
	*Coord_j;			/*!< \brief Cartesians coordinates of point j. */
  double *TurbVar_i,	/*!< \brief Vector of turbulent variables at point i. */
	*TurbVar_j;	/*!< \brief Vector of derivative of turbulent variables at point j. */
  double *Normal;	/*!< \brief Normal vector, it norm is the area of the face. */
  double *V_i,		/*!< \brief Vector of primitive variables at point i. */
	*V_j;				/*!< \brief Vector of primitive variables at point j. */
  double Density_i,	/*!< \brief Density at point i. */
	Density_j;			/*!< \brief Density at point j. */
  double Eddy_Viscosity_i,	/*!< \brief Eddy viscosity at point i. */
	Eddy_Viscosity_j;			/*!< \brief Eddy viscosity at point j. */
  double Laminar_Viscosity_i,	/*!< \brief Laminar viscosity at point i. */
	Laminar_Viscosity_j;		/*!< \brief Variation of laminar viscosity at point j. */
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
	 * \brief Set the gradient of the turbulent variables.
	 * \param[in] val_turbvar_grad_i - Gradient of the turbulent variable at point i.
	 * \param[in] val_turbvar_grad_j - Gradient of the turbulent variable at point j.
	 */
	void SetTurbVarGradient(double **val_turbvar_grad_i, double **val_turbvar_grad_j);
  
  /*!
	 * \brief Set coordinates of the points.
	 * \param[in] val_coord_i - Coordinates of the point i.
	 * \param[in] val_coord_j - Coordinates of the point j.
	 */
	void SetCoord(double *val_coord_i, double *val_coord_j);
  
  /*!
	 * \brief Set the value of the turbulent variable.
	 * \param[in] val_turbvar_i - Value of the turbulent variable at point i.
	 * \param[in] val_turbvar_j - Value of the turbulent variable at point j.
	 */
	void SetTurbVar(double *val_turbvar_i, double *val_turbvar_j);
  
  /*!
	 * \brief Set the value of the normal vector to the face between two points.
	 * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
	 */
	void SetNormal(double *val_normal);
  
  /*!
	 * \brief Set the value of the primitive variables.
	 * \param[in] val_v_i - Value of the primitive variable at point i.
	 * \param[in] val_v_j - Value of the primitive variable at point j.
	 */
	void SetPrimitive(double *val_v_i, double *val_v_j);
  
	/*!
	 * \brief Set the laminar viscosity.
	 * \param[in] val_laminar_viscosity_i - Value of the laminar viscosity at point i.
	 * \param[in] val_laminar_viscosity_j - Value of the laminar viscosity at point j.
	 */
	void SetLaminarViscosity(double val_laminar_viscosity_i,
                           double val_laminar_viscosity_j);
  
  /*!
	 * \brief Set the eddy viscosity.
	 * \param[in] val_eddy_viscosity_i - Value of the eddy viscosity at point i.
	 * \param[in] val_eddy_viscosity_j - Value of the eddy viscosity at point j.
	 */
	void SetEddyViscosity(double val_eddy_viscosity_i,
                        double val_eddy_viscosity_j);
  
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
 * \author F. Palacios.
 * \version 1.1.0
 */
class CAvgGradCorrected_Flow : public CNumerics {
private:
  unsigned short nDim, nVar;	/*!< \brief Number of dimensions and variables. */
	double Gamma;				/*!< \brief Fluid's Gamma constant (ratio of specific heats). */
	double Gamma_Minus_One;		/*!< \brief Fluids's Gamma - 1.0  . */
	double Gas_Constant;		 		/*!< \brief Gas constant. */
  double Area;				/*!< \brief Area of the face i-j. */
  double **PrimVar_Grad_i,	/*!< \brief Gradient of primitive variables at point i. */
	**PrimVar_Grad_j;			/*!< \brief Gradient of primitive variables at point j. */
  double *Coord_i,	/*!< \brief Cartesians coordinates of point i. */
	*Coord_j;			/*!< \brief Cartesians coordinates of point j. */
  double *Normal;	/*!< \brief Normal vector, it norm is the area of the face. */
  double *UnitNormal;		/*!< \brief Unitary normal vector. */
  double *V_i,		/*!< \brief Vector of primitive variables at point i. */
	*V_j;				/*!< \brief Vector of primitive variables at point j. */
  double turb_ke_i,	/*!< \brief Turbulent kinetic energy at point i. */
	turb_ke_j;			/*!< \brief Turbulent kinetic energy at point j. */
  double Eddy_Viscosity_i,	/*!< \brief Eddy viscosity at point i. */
	Eddy_Viscosity_j;			/*!< \brief Eddy viscosity at point j. */
  double Laminar_Viscosity_i,	/*!< \brief Laminar viscosity at point i. */
	Laminar_Viscosity_j;		/*!< \brief Variation of laminar viscosity at point j. */
  double
  **Flux_Tensor,	/*!< \brief Flux tensor (used for viscous and inviscid purposes. */
	*Proj_Flux_Tensor;		/*!< \brief Flux tensor projected in a direction. */
  double
  **tau,		/*!< \brief Viscous stress tensor. */
	**delta;			/*!< \brief Identity matrix. */
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
	void GetViscousProjJac(double *val_Mean_PrimVar,
                         double val_laminar_viscosity,
                         double val_eddy_viscosity,
                         double val_dist_ij,
                         double *val_normal, double val_dS,
                         double *val_Proj_Visc_Flux,
                         double **val_Proj_Jac_Tensor_i,
                         double **val_Proj_Jac_Tensor_j);
  
  /*!
	 * \brief Set the gradient of the primitive variables.
	 * \param[in] val_primvar_grad_i - Gradient of the primitive variable at point i.
	 * \param[in] val_primvar_grad_j - Gradient of the primitive variable at point j.
	 */
	void SetPrimVarGradient(double **val_primvar_grad_i,
                          double **val_primvar_grad_j);
  
  /*!
	 * \brief Set coordinates of the points.
	 * \param[in] val_coord_i - Coordinates of the point i.
	 * \param[in] val_coord_j - Coordinates of the point j.
	 */
	void SetCoord(double *val_coord_i, double *val_coord_j);

  /*!
	 * \brief Set the value of the normal vector to the face between two points.
	 * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
	 */
	void SetNormal(double *val_normal);
  
  /*!
	 * \brief Set the value of the primitive variables.
	 * \param[in] val_v_i - Value of the primitive variable at point i.
	 * \param[in] val_v_j - Value of the primitive variable at point j.
	 */
	void SetPrimitive(double *val_v_i, double *val_v_j);
  
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
	 * \brief Set the laminar viscosity.
	 * \param[in] val_laminar_viscosity_i - Value of the laminar viscosity at point i.
	 * \param[in] val_laminar_viscosity_j - Value of the laminar viscosity at point j.
	 */
	void SetLaminarViscosity(double val_laminar_viscosity_i,
                           double val_laminar_viscosity_j);
  
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
 * \author F. Palacios.
 * \version 1.1.0
 */
class CAvgGradCorrected_TurbSA : public CNumerics {
private:
  unsigned short nDim, nVar;	/*!< \brief Number of dimensions and variables. */
	double Gamma;				/*!< \brief Fluid's Gamma constant (ratio of specific heats). */
	double Gamma_Minus_One;		/*!< \brief Fluids's Gamma - 1.0  . */
	double Gas_Constant;		 		/*!< \brief Gas constant. */
  double **TurbVar_Grad_i,	/*!< \brief Gradient of turbulent variables at point i. */
	**TurbVar_Grad_j;			/*!< \brief Gradient of turbulent variables at point j. */
  double *Coord_i,	/*!< \brief Cartesians coordinates of point i. */
	*Coord_j;			/*!< \brief Cartesians coordinates of point j. */
  double *TurbVar_i,	/*!< \brief Vector of turbulent variables at point i. */
	*TurbVar_j;	/*!< \brief Vector of derivative of turbulent variables at point j. */
  double *Normal;	/*!< \brief Normal vector, it norm is the area of the face. */
  double *V_i,		/*!< \brief Vector of primitive variables at point i. */
	*V_j;				/*!< \brief Vector of primitive variables at point j. */
  double Laminar_Viscosity_i,	/*!< \brief Laminar viscosity at point i. */
	Laminar_Viscosity_j;		/*!< \brief Variation of laminar viscosity at point j. */
	double Eddy_Viscosity_i,	/*!< \brief Eddy viscosity at point i. */
	Eddy_Viscosity_j;			/*!< \brief Eddy viscosity at point j. */
	double **Mean_GradTurbVar;
	double *Proj_Mean_GradTurbVar_Kappa, *Proj_Mean_GradTurbVar_Edge, *Proj_Mean_GradTurbVar_Corrected;
	double *Edge_Vector;
	bool implicit;
	double sigma, nu_i, nu_j, nu_e, dist_ij_2, proj_vector_ij, nu_hat_i, nu_hat_j;
	unsigned short iVar, iDim;
  double Density_i,	/*!< \brief Density at point i. */
	Density_j;			/*!< \brief Density at point j. */
  
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
	 * \brief Set the gradient of the turbulent variables.
	 * \param[in] val_turbvar_grad_i - Gradient of the turbulent variable at point i.
	 * \param[in] val_turbvar_grad_j - Gradient of the turbulent variable at point j.
	 */
	void SetTurbVarGradient(double **val_turbvar_grad_i, double **val_turbvar_grad_j);
  
  /*!
	 * \brief Set the value of the turbulent variable.
	 * \param[in] val_turbvar_i - Value of the turbulent variable at point i.
	 * \param[in] val_turbvar_j - Value of the turbulent variable at point j.
	 */
	void SetTurbVar(double *val_turbvar_i, double *val_turbvar_j);
  
  /*!
	 * \brief Set coordinates of the points.
	 * \param[in] val_coord_i - Coordinates of the point i.
	 * \param[in] val_coord_j - Coordinates of the point j.
	 */
	void SetCoord(double *val_coord_i, double *val_coord_j);
  
  /*!
	 * \brief Set the value of the normal vector to the face between two points.
	 * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
	 */
	void SetNormal(double *val_normal);
  
  /*!
	 * \brief Set the value of the primitive variables.
	 * \param[in] val_v_i - Value of the primitive variable at point i.
	 * \param[in] val_v_j - Value of the primitive variable at point j.
	 */
	void SetPrimitive(double *val_v_i, double *val_v_j);
  
  /*!
	 * \brief Set the eddy viscosity.
	 * \param[in] val_eddy_viscosity_i - Value of the eddy viscosity at point i.
	 * \param[in] val_eddy_viscosity_j - Value of the eddy viscosity at point j.
	 */
	void SetEddyViscosity(double val_eddy_viscosity_i,
                        double val_eddy_viscosity_j);
  
	/*!
	 * \brief Set the laminar viscosity.
	 * \param[in] val_laminar_viscosity_i - Value of the laminar viscosity at point i.
	 * \param[in] val_laminar_viscosity_j - Value of the laminar viscosity at point j.
	 */
	void SetLaminarViscosity(double val_laminar_viscosity_i,
                           double val_laminar_viscosity_j);
  
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
 * \author F. Palacios.
 * \version 1.1.0
 */
class CAvgGrad_TurbSST : public CNumerics {
private:
  unsigned short nDim, nVar;	/*!< \brief Number of dimensions and variables. */
	double Gamma;				/*!< \brief Fluid's Gamma constant (ratio of specific heats). */
	double Gamma_Minus_One;		/*!< \brief Fluids's Gamma - 1.0  . */
	double Gas_Constant;		 		/*!< \brief Gas constant. */
  double **TurbVar_Grad_i,	/*!< \brief Gradient of turbulent variables at point i. */
	**TurbVar_Grad_j;			/*!< \brief Gradient of turbulent variables at point j. */
  double *Coord_i,	/*!< \brief Cartesians coordinates of point i. */
	*Coord_j;			/*!< \brief Cartesians coordinates of point j. */
  double *Normal;	/*!< \brief Normal vector, it norm is the area of the face. */
  double *V_i,		/*!< \brief Vector of primitive variables at point i. */
	*V_j;				/*!< \brief Vector of primitive variables at point j. */
  double Laminar_Viscosity_i,	/*!< \brief Laminar viscosity at point i. */
	Laminar_Viscosity_j;		/*!< \brief Variation of laminar viscosity at point j. */
	double Eddy_Viscosity_i,	/*!< \brief Eddy viscosity at point i. */
	Eddy_Viscosity_j;			/*!< \brief Eddy viscosity at point j. */
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
  double Density_i,	/*!< \brief Density at point i. */
	Density_j;			/*!< \brief Density at point j. */
  
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
	 * \brief Set the gradient of the turbulent variables.
	 * \param[in] val_turbvar_grad_i - Gradient of the turbulent variable at point i.
	 * \param[in] val_turbvar_grad_j - Gradient of the turbulent variable at point j.
	 */
	void SetTurbVarGradient(double **val_turbvar_grad_i, double **val_turbvar_grad_j);
  
  /*!
	 * \brief Set coordinates of the points.
	 * \param[in] val_coord_i - Coordinates of the point i.
	 * \param[in] val_coord_j - Coordinates of the point j.
	 */
	void SetCoord(double *val_coord_i, double *val_coord_j);
  
  /*!
	 * \brief Set the value of the normal vector to the face between two points.
	 * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
	 */
	void SetNormal(double *val_normal);
  
	/*!
	 * \brief Sets value of first blending function.
	 */
	void SetF1blending(double val_F1_i, double val_F1_j) { F1_i = val_F1_i; F1_j = val_F1_j;}
  
  /*!
	 * \brief Set the value of the primitive variables.
	 * \param[in] val_v_i - Value of the primitive variable at point i.
	 * \param[in] val_v_j - Value of the primitive variable at point j.
	 */
	void SetPrimitive(double *val_v_i, double *val_v_j);
  
	/*!
	 * \brief Set the laminar viscosity.
	 * \param[in] val_laminar_viscosity_i - Value of the laminar viscosity at point i.
	 * \param[in] val_laminar_viscosity_j - Value of the laminar viscosity at point j.
	 */
	void SetLaminarViscosity(double val_laminar_viscosity_i,
                           double val_laminar_viscosity_j);
  
  /*!
	 * \brief Set the eddy viscosity.
	 * \param[in] val_eddy_viscosity_i - Value of the eddy viscosity at point i.
	 * \param[in] val_eddy_viscosity_j - Value of the eddy viscosity at point j.
	 */
	void SetEddyViscosity(double val_eddy_viscosity_i,
                        double val_eddy_viscosity_j);
  
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
 * \author F. Palacios.
 * \version 1.1.0
 */
class CAvgGradCorrected_TurbSST : public CNumerics {
private:
  unsigned short nDim, nVar;	/*!< \brief Number of dimensions and variables. */
	double Gamma;				/*!< \brief Fluid's Gamma constant (ratio of specific heats). */
	double Gamma_Minus_One;		/*!< \brief Fluids's Gamma - 1.0  . */
	double Gas_Constant;		 		/*!< \brief Gas constant. */
  double **TurbVar_Grad_i,	/*!< \brief Gradient of turbulent variables at point i. */
	**TurbVar_Grad_j;			/*!< \brief Gradient of turbulent variables at point j. */

  double *TurbVar_i,	/*!< \brief Vector of turbulent variables at point i. */
	*TurbVar_j;	/*!< \brief Vector of derivative of turbulent variables at point j. */

  double *Coord_i,	/*!< \brief Cartesians coordinates of point i. */
	*Coord_j;			/*!< \brief Cartesians coordinates of point j. */
  double *Normal;	/*!< \brief Normal vector, it norm is the area of the face. */
  double *V_i,		/*!< \brief Vector of primitive variables at point i. */
	*V_j;				/*!< \brief Vector of primitive variables at point j. */

  double Laminar_Viscosity_i,	/*!< \brief Laminar viscosity at point i. */
	Laminar_Viscosity_j;		/*!< \brief Variation of laminar viscosity at point j. */
  
	double Eddy_Viscosity_i,	/*!< \brief Eddy viscosity at point i. */
	Eddy_Viscosity_j;			/*!< \brief Eddy viscosity at point j. */
  
  double Density_i,	/*!< \brief Density at point i. */
	Density_j;			/*!< \brief Density at point j. */
  
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
	 * \brief Set the gradient of the turbulent variables.
	 * \param[in] val_turbvar_grad_i - Gradient of the turbulent variable at point i.
	 * \param[in] val_turbvar_grad_j - Gradient of the turbulent variable at point j.
	 */
	void SetTurbVarGradient(double **val_turbvar_grad_i, double **val_turbvar_grad_j);
  
  /*!
	 * \brief Set the value of the turbulent variable.
	 * \param[in] val_turbvar_i - Value of the turbulent variable at point i.
	 * \param[in] val_turbvar_j - Value of the turbulent variable at point j.
	 */
	void SetTurbVar(double *val_turbvar_i, double *val_turbvar_j);
  
  /*!
	 * \brief Set coordinates of the points.
	 * \param[in] val_coord_i - Coordinates of the point i.
	 * \param[in] val_coord_j - Coordinates of the point j.
	 */
	void SetCoord(double *val_coord_i, double *val_coord_j);
  
  /*!
	 * \brief Set the value of the normal vector to the face between two points.
	 * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
	 */
	void SetNormal(double *val_normal);
  
  /*!
	 * \brief Set the value of the primitive variables.
	 * \param[in] val_v_i - Value of the primitive variable at point i.
	 * \param[in] val_v_j - Value of the primitive variable at point j.
	 */
	void SetPrimitive(double *val_v_i, double *val_v_j);
  
	/*!
	 * \brief Sets value of first blending function.
	 */
	void SetF1blending(double val_F1_i, double val_F1_j) { F1_i = val_F1_i; F1_j = val_F1_j;}
  
	/*!
	 * \brief Set the laminar viscosity.
	 * \param[in] val_laminar_viscosity_i - Value of the laminar viscosity at point i.
	 * \param[in] val_laminar_viscosity_j - Value of the laminar viscosity at point j.
	 */
	void SetLaminarViscosity(double val_laminar_viscosity_i,
                           double val_laminar_viscosity_j);
  
  /*!
	 * \brief Set the eddy viscosity.
	 * \param[in] val_eddy_viscosity_i - Value of the eddy viscosity at point i.
	 * \param[in] val_eddy_viscosity_j - Value of the eddy viscosity at point j.
	 */
	void SetEddyViscosity(double val_eddy_viscosity_i,
                        double val_eddy_viscosity_j);
  
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
 * \class CSourcePieceWise_TurbSA
 * \brief Class for integrating the source terms of the Spalart-Allmaras turbulence model equation.
 * \ingroup SourceDiscr
 * \author F. Palacios.
 * \version 1.1.0
 */
class CSourcePieceWise_TurbSA : public CNumerics {
private:
  unsigned short nDim, nVar;	/*!< \brief Number of dimensions and variables. */
	double Gamma;				/*!< \brief Fluid's Gamma constant (ratio of specific heats). */
	double Gamma_Minus_One;		/*!< \brief Fluids's Gamma - 1.0  . */
	double Gas_Constant;		 		/*!< \brief Gas constant. */
	double Volume;				/*!< \brief Volume of the control volume around point i. */
  double **PrimVar_Grad_i,	/*!< \brief Gradient of primitive variables at point i. */
	**PrimVar_Grad_j;			/*!< \brief Gradient of primitive variables at point j. */
  double **TurbVar_Grad_i,	/*!< \brief Gradient of turbulent variables at point i. */
	**TurbVar_Grad_j;			/*!< \brief Gradient of turbulent variables at point j. */
  double *TurbVar_i,	/*!< \brief Vector of turbulent variables at point i. */
	*TurbVar_j;	/*!< \brief Vector of derivative of turbulent variables at point j. */
  double *V_i,		/*!< \brief Vector of primitive variables at point i. */
	*V_j;				/*!< \brief Vector of primitive variables at point j. */
  double dist_i,	/*!< \brief Distance of point i to the nearest wall. */
	dist_j;			/*!< \brief Distance of point j to the nearest wall. */
  double Laminar_Viscosity_i,	/*!< \brief Laminar viscosity at point i. */
	Laminar_Viscosity_j;		/*!< \brief Variation of laminar viscosity at point j. */
	double Density_i,	/*!< \brief Density at point i. */
	Density_j;			/*!< \brief Density at point j. */
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
	 * \brief Set the value of the volume of the control volume.
	 * \param[in] val_volume Volume of the control volume.
	 */
	void SetVolume(double val_volume);
  
  /*!
	 * \brief Set the gradient of the primitive variables.
	 * \param[in] val_primvar_grad_i - Gradient of the primitive variable at point i.
	 * \param[in] val_primvar_grad_j - Gradient of the primitive variable at point j.
	 */
	void SetPrimVarGradient(double **val_primvar_grad_i,
                          double **val_primvar_grad_j);
  
  /*!
	 * \brief Set the gradient of the turbulent variables.
	 * \param[in] val_turbvar_grad_i - Gradient of the turbulent variable at point i.
	 * \param[in] val_turbvar_grad_j - Gradient of the turbulent variable at point j.
	 */
	void SetTurbVarGradient(double **val_turbvar_grad_i, double **val_turbvar_grad_j);
  
  /*!
	 * \brief Set the value of the turbulent variable.
	 * \param[in] val_turbvar_i - Value of the turbulent variable at point i.
	 * \param[in] val_turbvar_j - Value of the turbulent variable at point j.
	 */
	void SetTurbVar(double *val_turbvar_i, double *val_turbvar_j);
  
  /*!
	 * \brief Set the value of the primitive variables.
	 * \param[in] val_v_i - Value of the primitive variable at point i.
	 * \param[in] val_v_j - Value of the primitive variable at point j.
	 */
	void SetPrimitive(double *val_v_i, double *val_v_j);
  
	/*!
	 * \brief Residual for source term integration.
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
	 * \param[in] config - Definition of the particular problem.
	 */
	void ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config);
  
	/*!
	 * \brief Set the laminar viscosity.
	 * \param[in] val_laminar_viscosity_i - Value of the laminar viscosity at point i.
	 * \param[in] val_laminar_viscosity_j - Value of the laminar viscosity at point j.
	 */
	void SetLaminarViscosity(double val_laminar_viscosity_i,
                           double val_laminar_viscosity_j);
  
  /*!
	 * \brief Set the value of the distance from the nearest wall.
	 * \param[in] val_dist_i - Value of of the distance from point i to the nearest wall.
	 * \param[in] val_dist_j - Value of of the distance from point j to the nearest wall.
	 */
	void SetDistance(double val_dist_i, double val_dist_j);
  
	/*!
	 * \brief Residual for source term integration.
	 * \param[in] intermittency_in - Value of the intermittency.
	 */
  void SetIntermittency(double intermittency_in);

};

/*!
 * \class CSourcePieceWise_TurbSST
 * \brief Class for integrating the source terms of the Menter SST turbulence model equations.
 * \ingroup SourceDiscr
 * \author A. Campos.
 * \version 1.1.0
 */
class CSourcePieceWise_TurbSST : public CNumerics {
private:
  unsigned short nDim, nVar;	/*!< \brief Number of dimensions and variables. */
	double Gamma;				/*!< \brief Fluid's Gamma constant (ratio of specific heats). */
	double Gamma_Minus_One;		/*!< \brief Fluids's Gamma - 1.0  . */
	double Gas_Constant;		 		/*!< \brief Gas constant. */
  double Volume;				/*!< \brief Volume of the control volume around point i. */
  double **PrimVar_Grad_i,	/*!< \brief Gradient of primitive variables at point i. */
	**PrimVar_Grad_j;			/*!< \brief Gradient of primitive variables at point j. */
  double *TurbVar_i,	/*!< \brief Vector of turbulent variables at point i. */
	*TurbVar_j;	/*!< \brief Vector of derivative of turbulent variables at point j. */
  double *V_i,		/*!< \brief Vector of primitive variables at point i. */
	*V_j;				/*!< \brief Vector of primitive variables at point j. */
  double dist_i,	/*!< \brief Distance of point i to the nearest wall. */
	dist_j;			/*!< \brief Distance of point j to the nearest wall. */
  double Laminar_Viscosity_i,	/*!< \brief Laminar viscosity at point i. */
	Laminar_Viscosity_j;		/*!< \brief Variation of laminar viscosity at point j. */
	double Eddy_Viscosity_i,	/*!< \brief Eddy viscosity at point i. */
	Eddy_Viscosity_j;			/*!< \brief Eddy viscosity at point j. */
  double Density_i,	/*!< \brief Density at point i. */
	Density_j;			/*!< \brief Density at point j. */
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
	 * \brief Set the value of the volume of the control volume.
	 * \param[in] val_volume Volume of the control volume.
	 */
	void SetVolume(double val_volume);
  
	/*!
	 * \brief Set the gradient of the primitive variables.
	 * \param[in] val_primvar_grad_i - Gradient of the primitive variable at point i.
	 * \param[in] val_primvar_grad_j - Gradient of the primitive variable at point j.
	 */
	void SetPrimVarGradient(double **val_primvar_grad_i,
                          double **val_primvar_grad_j);
  
  /*!
	 * \brief Set the value of the turbulent variable.
	 * \param[in] val_turbvar_i - Value of the turbulent variable at point i.
	 * \param[in] val_turbvar_j - Value of the turbulent variable at point j.
	 */
	void SetTurbVar(double *val_turbvar_i, double *val_turbvar_j);
  
  /*!
	 * \brief Set the value of the primitive variables.
	 * \param[in] val_v_i - Value of the primitive variable at point i.
	 * \param[in] val_v_j - Value of the primitive variable at point j.
	 */
	void SetPrimitive(double *val_v_i, double *val_v_j);
  
  /*!
	 * \brief Set the value of the distance from the nearest wall.
	 * \param[in] val_dist_i - Value of of the distance from point i to the nearest wall.
	 * \param[in] val_dist_j - Value of of the distance from point j to the nearest wall.
	 */
	void SetDistance(double val_dist_i, double val_dist_j);
  
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
	 * \brief Set the laminar viscosity.
	 * \param[in] val_laminar_viscosity_i - Value of the laminar viscosity at point i.
	 * \param[in] val_laminar_viscosity_j - Value of the laminar viscosity at point j.
	 */
	void SetLaminarViscosity(double val_laminar_viscosity_i,
                           double val_laminar_viscosity_j);
  
  /*!
	 * \brief Set the eddy viscosity.
	 * \param[in] val_eddy_viscosity_i - Value of the eddy viscosity at point i.
	 * \param[in] val_eddy_viscosity_j - Value of the eddy viscosity at point j.
	 */
	void SetEddyViscosity(double val_eddy_viscosity_i,
                        double val_eddy_viscosity_j);
	/*!
	 * \brief Residual for source term integration.
	 * \param[out] val_residual - Pointer to the total residual.
	 * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
	 * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
	 * \param[in] config - Definition of the particular problem.
	 */
	void ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config);
};


/*!
 * \class CSourceNothing
 * \brief Dummy class.
 * \ingroup SourceDiscr
 * \author F. Palacios.
 * \version 1.1.0
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

#include "numerics_structure.inl"
