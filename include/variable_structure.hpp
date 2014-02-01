/*!
 * \file variable_structure.hpp
 * \brief Headers of the main subroutines for storing all the variables for
 *        each kind of governing equation (direct, adjoint and linearized).
 *        The subroutines and functions are in the <i>variable_structure.cpp</i> file.
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
#include <cstdlib>

#include "config_structure.hpp"

using namespace std;

/*!
 * \class CVariable
 * \brief Main class for defining the variables.
 * \author F. Palacios.
 * \version 1.0.0
 */
class CVariable {
protected:
  
	double *Solution,		/*!< \brief Solution of the problem. */
	*Solution_Old;			/*!< \brief Old solution of the problem R-K. */
	double *Solution_time_n,	/*!< \brief Solution of the problem at time n for dual-time stepping technique. */
	*Solution_time_n1;			/*!< \brief Solution of the problem at time n-1 for dual-time stepping technique. */
	double **Gradient;		/*!< \brief Gradient of the solution of the problem. */
	double *Limiter;				/*!< \brief Limiter of the solution of the problem. */
	double *Solution_Max;		/*!< \brief Max solution for limiter computation. */
	double *Solution_Min;		/*!< \brief Min solution for limiter computation. */
	double AuxVar;			/*!< \brief Auxiliar variable for gradient computation. */
	double *Grad_AuxVar;	/*!< \brief Gradient of the auxiliar variable. */
	double Delta_Time;	/*!< \brief Time step. */
	double Max_Lambda,	/*!< \brief Maximun eingenvalue. */
	Max_Lambda_Inv,		/*!< \brief Maximun inviscid eingenvalue. */
	Max_Lambda_Visc,	/*!< \brief Maximun viscous eingenvalue. */
	Lambda;				/*!< \brief Value of the eingenvalue. */
	double Sensor;	/*!< \brief Pressure sensor for high order central scheme. */
	double *Undivided_Laplacian;	/*!< \brief Undivided laplacian of the solution. */
	double *Res_TruncError,	/*!< \brief Truncation error for multigrid cycle. */
	*Residual_Old,		/*!< \brief Auxiliar structure for residual smoothing. */
	*Residual_Sum;		/*!< \brief Auxiliar structure for residual smoothing. */
	static unsigned short nDim;		/*!< \brief Number of dimension of the problem. */
	unsigned short nVar;		/*!< \brief Number of variables of the problem,
													 note that this variable cannnot be static, it is possible to
													 have different number of nVar in the same problem. */
  unsigned short nPrimVar, nPrimVarGrad;		/*!< \brief Number of variables of the problem,
                                             note that this variable cannnot be static, it is possible to
                                             have different number of nVar in the same problem. */
  
public:
  
	/*!
	 * \brief Constructor of the class.
	 */
	CVariable(void);
  
  /*!
	 * \overload
	 * \param[in] val_nvar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CVariable(unsigned short val_nvar, CConfig *config);
  
	/*!
	 * \overload
	 * \param[in] val_ndim - Number of dimensions of the problem.
	 * \param[in] val_nvar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CVariable(unsigned short val_ndim, unsigned short val_nvar, CConfig *config);
  
	/*!
	 * \brief Destructor of the class.
	 */
	virtual ~CVariable(void);
  
	/*!
	 * \brief Set the value of the solution.
	 * \param[in] val_solution - Solution of the problem.
	 */
	void SetSolution(double *val_solution);
  
	/*!
	 * \overload
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_solution - Value of the solution for the index <i>val_var</i>.
	 */
	void SetSolution(unsigned short val_var, double val_solution);
  
	/*!
	 * \brief Get the solution.
	 * \param[in] val_var - Index of the variable.
	 * \return Value of the solution for the index <i>val_var</i>.
	 */
	double GetSolution(unsigned short val_var);
  
	/*!
	 * \brief Get the old solution of the problem (Runge-Kutta method)
	 * \param[in] val_var - Index of the variable.
	 * \return Pointer to the old solution vector.
	 */
	double GetSolution_Old(unsigned short val_var);
  
	/*!
	 * \brief Set the value of the old solution.
	 * \param[in] val_solution_old - Pointer to the residual vector.
	 */
	void SetSolution_Old(double *val_solution_old);
  
	/*!
	 * \overload
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_solution_old - Value of the old solution for the index <i>val_var</i>.
	 */
	void SetSolution_Old(unsigned short val_var, double val_solution_old);
  
	/*!
	 * \brief Set old variables to the value of the current variables.
	 */
	void Set_OldSolution(void);
  
	/*!
	 * \brief Set to zero the velocity components of the solution.
	 */
	void SetVelSolutionZero(void);
  
	/*!
	 * \brief Set to zero velocity components of the solution.
	 */
	void SetVelSolutionOldZero(void);
  
	/*!
	 * \brief Set to zero the solution.
	 */
	void SetSolutionZero(void);
  
  /*!
	 * \brief Set to zero a particular solution.
	 */
  void SetSolutionZero(unsigned short val_var);
  
	/*!
	 * \brief Add a value to the solution.
	 * \param[in] val_var - Number of the variable.
	 * \param[in] val_solution - Value that we want to add to the solution.
	 */
	void AddSolution(unsigned short val_var, double val_solution);
  
  /*!
	 * \brief Add a value to the solution, clipping the values.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_solution - Value of the solution change.
   * \param[in] lowerlimit - Lower value.
   * \param[in] upperlimit - Upper value.
	 */
	void AddClippedSolution(unsigned short val_var, double val_solution,
                          double lowerlimit, double upperlimit);
  
	/*!
	 * \brief Update the variables using a conservative format.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_solution - Value of the solution change.
	 * \param[in] val_density - Value of the density.
	 * \param[in] val_density_old - Value of the old density.
	 */
	void AddConservativeSolution(unsigned short val_var, double val_solution,
                               double val_density, double val_density_old, double lowerlimit,
                               double upperlimit);
  
	/*!
	 * \brief Get the solution of the problem.
	 * \return Pointer to the solution vector.
	 */
	double *GetSolution(void);
  
	/*!
	 * \brief Get the old solution of the problem (Runge-Kutta method)
	 * \return Pointer to the old solution vector.
	 */
	double *GetSolution_Old(void);
  
	/*!
	 * \brief Set the value of the old residual.
	 * \param[in] val_residual_old - Pointer to the residual vector.
	 */
	void SetResidual_Old(double *val_residual_old);
  
	/*!
	 * \brief Add a value to the summed residual vector.
	 * \param[in] val_residual - Pointer to the residual vector.
	 */
	void AddResidual_Sum(double *val_residual);
  
	/*!
	 * \brief Set summed residual vector to zero value.
	 */
	void SetResidualSumZero(void);
  
	/*!
	 * \brief Get the value of the summed residual.
	 * \return Pointer to the summed residual.
	 */
	double *GetResidual_Sum(void);
  
	/*!
	 * \brief Get the value of the old residual.
	 * \return Pointer to the old residual.
	 */
	double *GetResidual_Old(void);
  
	/*!
	 * \brief Get the value of the summed residual.
	 * \param[out] val_residual - Pointer to the summed residual.
	 */
	void GetResidual_Sum(double *val_residual);
  
	/*!
	 * \brief Set auxiliar variables, we are looking for the gradient of that variable.
	 * \param[in] val_auxvar - Value of the auxiliar variable.
	 */
	void SetAuxVar(double val_auxvar);
  
	/*!
	 * \brief Get the value of the auxiliary variable.
	 * \return Value of the auxiliary variable.
	 */
	double GetAuxVar(void);
  
	/*!
	 * \brief Set the auxiliary variable gradient to zero value.
	 */
	void SetAuxVarGradientZero(void);
  
	/*!
	 * \brief Set the value of the auxiliary variable gradient.
	 * \param[in] val_dim - Index of the dimension.
	 * \param[in] val_gradient - Value of the gradient for the index <i>val_dim</i>.
	 */
	void SetAuxVarGradient(unsigned short val_dim, double val_gradient);
  
	/*!
	 * \brief Add a value to the auxiliary variable gradient.
	 * \param[in] val_dim - Index of the dimension.
	 * \param[in] val_value - Value of the gradient to be added for the index <i>val_dim</i>.
	 */
	void AddAuxVarGradient(unsigned short val_dim, double val_value);
  
	/*!
	 * \brief Subtract a value to the auxiliary variable gradient.
	 * \param[in] val_dim - Index of the dimension.
	 * \param[in] val_value - Value of the gradient to be subtracted for the index <i>val_dim</i>.
	 */
	void SubtractAuxVarGradient(unsigned short val_dim, double val_value);
  
	/*!
	 * \brief Get the gradient of the auxiliary variable.
	 * \return Value of the gradient of the auxiliary variable.
	 */
	double *GetAuxVarGradient(void);
  
	/*!
	 * \brief Get the gradient of the auxiliary variable.
	 * \param[in] val_dim - Index of the dimension.
	 * \return Value of the gradient of the auxiliary variable for the dimension <i>val_dim</i>.
	 */
	double GetAuxVarGradient(unsigned short val_dim);
  
	/*!
	 * \brief Add a value to the truncation error.
	 * \param[in] val_truncation_error - Value that we want to add to the truncation error.
	 */
	void AddRes_TruncError(double *val_truncation_error);
  
	/*!
	 * \brief Subtract a value to the truncation error.
	 * \param[in] val_truncation_error - Value that we want to subtract to the truncation error.
	 */
	void SubtractRes_TruncError(double *val_truncation_error);
  
	/*!
	 * \brief Set the truncation error to zero.
	 */
	void SetRes_TruncErrorZero(void);
  
  /*!
	 * \brief Set the truncation error to zero.
	 */
	void SetVal_ResTruncError_Zero(unsigned short val_var);
  
	/*!
	 * \brief Set the velocity of the truncation error to zero.
	 */
	void SetVel_ResTruncError_Zero(void);
  
  /*!
	 * \brief Set the velocity of the truncation error to zero.
	 */
	void SetEnergy_ResTruncError_Zero(void);
  
	/*!
	 * \brief Get the truncation error.
	 * \return Pointer to the truncation error.
	 */
	double *GetResTruncError(void);
  
	/*!
	 * \brief Get the truncation error.
	 * \param[out] val_trunc_error - Pointer to the truncation error.
	 */
	void GetResTruncError(double *val_trunc_error);
  
	/*!
	 * \brief Set the gradient of the solution.
	 * \param[in] val_gradient - Gradient of the solution.
	 */
	void SetGradient(double **val_gradient);
  
	/*!
	 * \overload
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_dim - Index of the dimension.
	 * \param[in] val_value - Value of the gradient.
	 */
	void SetGradient(unsigned short val_var, unsigned short val_dim, double val_value);
  
	/*!
	 * \brief Set to zero the gradient of the solution.
	 */
	void SetGradientZero(void);
  
	/*!
	 * \brief Add <i>val_value</i> to the solution gradient.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_dim - Index of the dimension.
	 * \param[in] val_value - Value to add to the solution gradient.
	 */
	void AddGradient(unsigned short val_var, unsigned short val_dim, double val_value);
  
	/*!
	 * \brief Subtract <i>val_value</i> to the solution gradient.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_dim - Index of the dimension.
	 * \param[in] val_value - Value to subtract to the solution gradient.
	 */
	void SubtractGradient(unsigned short val_var, unsigned short val_dim, double val_value);
  
	/*!
	 * \brief Get the value of the solution gradient.
	 * \return Value of the gradient solution.
	 */
	double **GetGradient(void);
  
	/*!
	 * \brief Get the value of the solution gradient.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_dim - Index of the dimension.
	 * \return Value of the solution gradient.
	 */
	double GetGradient(unsigned short val_var, unsigned short val_dim);
  
	/*!
	 * \brief Set the value of the limiter.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_limiter - Value of the limiter for the index <i>val_var</i>.
	 */
	void SetLimiter(unsigned short val_var, double val_limiter);
	
	/*!
	 * \brief Set the value of the max solution.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_limiter - Value of the max solution for the index <i>val_var</i>.
	 */
	void SetSolution_Max(unsigned short val_var, double val_solution);
	
	/*!
	 * \brief Set the value of the min solution.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_limiter - Value of the min solution for the index <i>val_var</i>.
	 */
	void SetSolution_Min(unsigned short val_var, double val_solution);
  
	/*!
	 * \brief Get the value of the slope limiter.
	 * \return Pointer to the limiters vector.
	 */
	double *GetLimiter(void);
  
	/*!
	 * \brief Get the value of the slope limiter.
	 * \param[in] val_var - Index of the variable.
	 * \return Value of the limiter vector for the variable <i>val_var</i>.
	 */
	double GetLimiter(unsigned short val_var);
	
	/*!
	 * \brief Get the value of the min solution.
	 * \param[in] val_var - Index of the variable.
	 * \return Value of the min solution for the variable <i>val_var</i>.
	 */
	double GetSolution_Max(unsigned short val_var);
	
	/*!
	 * \brief Get the value of the min solution.
	 * \param[in] val_var - Index of the variable.
	 * \return Value of the min solution for the variable <i>val_var</i>.
	 */
	double GetSolution_Min(unsigned short val_var);
  
	/*!
	 * \brief Get the value of the preconditioner Beta.
	 * \return Value of the low Mach preconditioner variable Beta
	 */
	virtual double GetPreconditioner_Beta();
  
	/*!
	 * \brief Set the value of the preconditioner Beta.
	 * \param[in] Value of the low Mach preconditioner variable Beta
	 */
	virtual void SetPreconditioner_Beta(double val_Beta);
  
	/*!
	 * \brief Set the value of the time step.
	 * \param[in] val_delta_time - Value of the time step.
	 */
	void SetDelta_Time(double val_delta_time);
  
	/*!
	 * \brief Get the value of the time step.
	 * \return Value of the time step.
	 */
	double GetDelta_Time(void);
  
	/*!
	 * \brief Set the value of the maximum eigenvalue.
	 * \param[in] val_max_lambda - Value of the maximum eigenvalue.
	 */
	void SetMax_Lambda(double val_max_lambda);
  
	/*!
	 * \brief Set the value of the maximum eigenvalue for the inviscid terms of the PDE.
	 * \param[in] val_max_lambda - Value of the maximum eigenvalue for the inviscid terms of the PDE.
	 */
	void SetMax_Lambda_Inv(double val_max_lambda);
  
	/*!
	 * \brief Set the value of the maximum eigenvalue for the viscous terms of the PDE.
	 * \param[in] val_max_lambda - Value of the maximum eigenvalue for the viscous terms of the PDE.
	 */
	void SetMax_Lambda_Visc(double val_max_lambda);
  
	/*!
	 * \brief Add a value to the maximum eigenvalue.
	 * \param[in] val_max_lambda - Value of the maximum eigenvalue.
	 */
	void AddMax_Lambda(double val_max_lambda);
  
	/*!
	 * \brief Add a value to the maximum eigenvalue for the inviscid terms of the PDE.
	 * \param[in] val_max_lambda - Value of the maximum eigenvalue for the inviscid terms of the PDE.
	 */
	void AddMax_Lambda_Inv(double val_max_lambda);
  
	/*!
	 * \brief Add a value to the maximum eigenvalue for the viscous terms of the PDE.
	 * \param[in] val_max_lambda - Value of the maximum eigenvalue for the viscous terms of the PDE.
	 */
	void AddMax_Lambda_Visc(double val_max_lambda);
  
	/*!
	 * \brief Get the value of the maximum eigenvalue.
	 * \return the value of the maximum eigenvalue.
	 */
	double GetMax_Lambda(void);
  
	/*!
	 * \brief Get the value of the maximum eigenvalue for the inviscid terms of the PDE.
	 * \return the value of the maximum eigenvalue for the inviscid terms of the PDE.
	 */
	double GetMax_Lambda_Inv(void);
  
	/*!
	 * \brief Get the value of the maximum eigenvalue for the viscous terms of the PDE.
	 * \return the value of the maximum eigenvalue for the viscous terms of the PDE.
	 */
	double GetMax_Lambda_Visc(void);
  
	/*!
	 * \brief Set the value of the spectral radius.
	 * \param[in] val_lambda - Value of the spectral radius.
	 */
	void SetLambda(double val_lambda);
  
	/*!
	 * \brief Add the value of the spectral radius.
	 * \param[in] val_lambda - Value of the spectral radius.
	 */
	void AddLambda(double val_lambda);
  
	/*!
	 * \brief Get the value of the spectral radius.
	 * \return Value of the spectral radius.
	 */
	double GetLambda(void);
  
	/*!
	 * \brief Set pressure sensor.
	 * \param[in] val_sensor - Value of the pressure sensor.
	 */
	void SetSensor(double val_sensor);
  
	/*!
	 * \brief Get the pressure sensor.
	 * \return Value of the pressure sensor.
	 */
	double GetSensor(void);
  
	/*!
	 * \brief Set the value of the undivided laplacian of the solution.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_undivided_laplacian - Value of the undivided solution for the index <i>val_var</i>.
	 */
	void SetUndivided_Laplacian(unsigned short val_var, double val_undivided_laplacian);
  
	/*!
	 * \brief Add the value of the undivided laplacian of the solution.
	 * \param[in] val_und_lapl - Value of the undivided solution.
	 */
	void AddUnd_Lapl(double *val_und_lapl);
  
	/*!
	 * \brief Subtract the value of the undivided laplacian of the solution.
	 * \param[in] val_und_lapl - Value of the undivided solution.
	 */
	void SubtractUnd_Lapl(double *val_und_lapl);
  
	/*!
	 * \brief Subtract the value of the undivided laplacian of the solution.
	 * \param[in] val_var - Variable of the undivided laplacian.
	 * \param[in] val_und_lapl - Value of the undivided solution.
	 */
	void SubtractUnd_Lapl(unsigned short val_var, double val_und_lapl);
  
	/*!
	 * \brief Set the undivided laplacian of the solution to zero.
	 */
	void SetUnd_LaplZero(void);
  
	/*!
	 * \brief Set a value to the undivided laplacian.
	 * \param[in] val_var - Variable of the undivided laplacian.
	 * \param[in] val_und_lapl - Value of the undivided laplacian.
	 */
	void SetUnd_Lapl(unsigned short val_var, double val_und_lapl);
  
	/*!
	 * \brief Get the undivided laplacian of the solution.
	 * \return Pointer to the undivided laplacian vector.
	 */
	double *GetUndivided_Laplacian(void);
  
	/*!
	 * \brief Get the undivided laplacian of the solution.
	 * \param[in] val_var - Variable of the undivided laplacian.
	 * \return Value of the undivided laplacian vector.
	 */
	double GetUndivided_Laplacian(unsigned short val_var);
  
	/*!
	 * \brief A virtual member.
	 * \return Value of the flow density.
	 */
	virtual double GetDensity(void);
  
	/*!
	 * \brief A virtual member.
	 * \return Value of the flow energy.
	 */
	virtual double GetEnergy(void);
  
	/*!
	 * \brief A virtual member.
	 * \return Value of the eddy viscosity.
	 */
	virtual double GetEddyViscosity(void);
  
	/*!
	 * \brief A virtual member.
	 * \return Value of the flow enthalpy.
	 */
	virtual double GetEnthalpy(void);
  
  /*!
	 * \brief A virtual member.
	 * \return Value of the flow pressure.
	 */
	virtual double GetPressure(void);
  
	/*!
	 * \brief A virtual member.
	 * \param[in] val_vector - Direction of projection.
	 * \return Value of the projected velocity.
	 */
	virtual double GetProjVel(double *val_vector);
  
	/*!
	 * \brief A virtual member.
	 * \return Value of the sound speed.
	 */
	virtual double GetSoundSpeed(void);
  
	/*!
	 * \brief A virtual member.
	 * \return Value of the temperature.
	 */
	virtual double GetTemperature(void);
  
  /*!
	 * \brief A virtual member.
	 * \param[in] val_dim - Index of the dimension.
	 * \return Value of the velocity for the dimension <i>val_dim</i>.
	 */
	virtual double GetVelocity(unsigned short val_dim);
  
	/*!
	 * \brief A virtual member.
	 * \return Norm 2 of the velocity vector.
	 */
	virtual double GetVelocity2(void);
  
	/*!
	 * \brief A virtual member.
	 * \return The laminar viscosity of the flow.
	 */
	virtual double GetLaminarViscosity(void);
  
	/*!
	 * \brief A virtual member.
	 * \param[in] val_dim - Index of the dimension.
	 * \return Value of the vorticity.
	 */
	virtual double GetVorticity(unsigned short val_dim);
  
	/*!
	 * \brief A virtual member.
	 * \return Value of the rate of strain magnitude.
	 */
	virtual double GetStrainMag(void);
  
	/*!
	 * \brief A virtual member.
	 * \param[in] val_eddy_viscosity - Value of the eddy viscosity.
	 */
	virtual void SetEddyViscosity(double eddy_visc);
  
	/*!
	 * \brief A virtual member.
	 */
	virtual void SetEnthalpy(void);
	
	/*!
	 * \brief A virtual member.
	 */
	virtual bool SetPrimVar_Compressible(CConfig *config);
    
	/*!
	 * \brief A virtual member.
	 */
	virtual bool SetPrimVar_Compressible(double eddy_visc, double turb_ke, CConfig *config);
	
	/*!
	 * \brief A virtual member.
	 */
	virtual double GetPrimVar(unsigned short val_var);
  
  /*!
	 * \brief A virtual member.
	 */
  virtual void SetPrimVar(unsigned short val_var, double val_prim);
  
  /*!
	 * \brief A virtual member.
	 */
  virtual void SetPrimVar(double *val_prim);
  
	/*!
	 * \brief A virtual member.
	 */
	virtual double *GetPrimVar(void);
  
	/*!
	 * \brief A virtual member.
	 * \param[in] Gamma - Ratio of Specific heats
	 */
	virtual bool SetPressure(double Gamma);
  
	/*!
	 * \brief A virtual member.
	 * \param[in] Gamma - Ratio of Specific heats
	 */
	virtual bool SetPressure(CConfig *config);
  
	/*!
	 * \brief A virtual member.
	 */
	virtual bool SetPressure(double Gamma, double turb_ke);
  
	/*!
	 * \brief A virtual member.
	 */
	virtual void SetPressure(void);
  
  /*!
	 * \brief A virtual member.
	 */
	virtual bool SetDensity(void);
  
	/*!
	 * \brief A virtual member.
	 * \param[in] Gamma - Ratio of specific heats.
	 */
	virtual bool SetSoundSpeed(double Gamma);
  
	/*!
	 * \brief A virtual member.
	 * \param[in] config - Configuration parameters.
	 */
	virtual bool SetSoundSpeed(CConfig *config);
  
	/*!
	 * \brief A virtual member.
	 */
	virtual bool SetSoundSpeed(void);
  
	/*!
	 * \brief A virtual member.
	 * \param[in] Gas_Constant - Value of the Gas Constant
	 */
	virtual bool SetTemperature(double Gas_Constant);
  
  /*!
	 * \brief A virtual member.
	 * \param[in] config - Configuration parameters.
	 */
	virtual bool SetTemperature(CConfig *config);
  
	/*!
	 * \brief A virtual member.
	 * \param[in] config - Configuration parameters.
	 */
	virtual void SetPrimVar(CConfig *config);
  
  /*!
	 * \brief A virtual member.
	 * \param[in] config - Configuration parameters.
	 */
	virtual void SetPrimVar(CConfig *config, double *Coord);
  
	/*!
	 * \brief A virtual member.
	 */
	virtual void SetVelocity(void);
  
	/*!
	 * \brief A virtual member.
	 */
	virtual void SetVelocity2(void);
  
	/*!
	 * \brief A virtual member.
	 * \param[in] val_velocity - Pointer to the velocity.
	 */
	virtual void SetVelocity_Old(double *val_velocity);
  
	/*!
	 * \brief A virtual member.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void SetLaminarViscosity(CConfig *config);
  
	/*!
	 * \brief A virtual member.
	 */
	virtual void SetVorticity(void);
  
	/*!
	 * \brief A virtual member.
	 */
	virtual void SetStrainMag(void);
  
	/*!
	 * \brief A virtual member.
	 */
	virtual void SetGradient_PrimitiveZero(unsigned short val_primvar);
  
	/*!
	 * \brief A virtual member.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_dim - Index of the dimension.
	 * \param[in] val_value - Value to add to the gradient of the primitive variables.
	 */
	virtual void AddGradient_Primitive(unsigned short val_var, unsigned short val_dim, double val_value);
  
	/*!
	 * \brief A virtual member.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_dim - Index of the dimension.
	 * \param[in] val_value - Value to subtract to the gradient of the primitive variables.
	 */
	virtual void SubtractGradient_Primitive(unsigned short val_var, unsigned short val_dim, double val_value);
  
	/*!
	 * \brief A virtual member.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_dim - Index of the dimension.
	 * \return Value of the primitive variables gradient.
	 */
	virtual double GetGradient_Primitive(unsigned short val_var, unsigned short val_dim);
  
  /*!
	 * \brief A virtual member.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_dim - Index of the dimension.
	 * \return Value of the primitive variables gradient.
	 */
	virtual double GetLimiter_Primitive(unsigned short val_var);
  
	/*!
	 * \brief A virtual member.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_dim - Index of the dimension.
	 * \param[in] val_value - Value of the gradient.
	 */
	virtual void SetGradient_Primitive(unsigned short val_var, unsigned short val_dim, double val_value);
  
  /*!
	 * \brief A virtual member.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_dim - Index of the dimension.
	 * \param[in] val_value - Value of the gradient.
	 */
	virtual void SetLimiter_Primitive(unsigned short val_var, double val_value);
  
	/*!
	 * \brief A virtual member.
	 * \return Value of the primitive variables gradient.
	 */
	virtual double **GetGradient_Primitive(void);
  
  /*!
	 * \brief A virtual member.
	 * \return Value of the primitive variables gradient.
	 */
	virtual double *GetLimiter_Primitive(void);
  
	/*!
	 * \brief Set the blending function for the blending of k-w and k-eps.
	 * \param[in] val_viscosity - Value of the vicosity.
	 * \param[in] val_density - Value of the density.
	 * \param[in] val_dist - Value of the distance to the wall.
	 */
	virtual void SetBlendingFunc(double val_viscosity, double val_dist, double val_density);
  
	/*!
	 * \brief Get the first blending function of the SST model.
	 */
	virtual double GetF1blending(void);
  
	/*!
	 * \brief Get the second blending function of the SST model.
	 */
	virtual double GetF2blending(void);
  
	/*!
	 * \brief Get the value of the cross diffusion of tke and omega.
	 */
	virtual double GetCrossDiff(void){ return 0.0; };
  
	/*!
	 * \brief Get the value of the eddy viscosity.
	 * \return the value of the eddy viscosity.
	 */
	virtual double GetmuT(void);
  
	/*!
	 * \brief Set the value of the eddy viscosity.
	 * \param[in] val_muT
	 */
	virtual void SetmuT(double val_muT);
  
};

/*!
 * \class CEulerVariable
 * \brief Main class for defining the variables of the Euler's solver.
 * \ingroup Euler_Equations
 * \author F. Palacios.
 * \version 1.0.0
 */
class CEulerVariable : public CVariable {
protected:
	double Velocity2;			/*!< \brief Square of the velocity vector. */
	double Precond_Beta;	/*!< \brief Low Mach number preconditioner value, Beta. */
  
	/*--- Primitive variable definition ---*/
	double *Primitive;	/*!< \brief Primitive variables (T,vx,vy,vz,P,rho,h,c) in compressible flows. */
	double **Gradient_Primitive;	/*!< \brief Gradient of the primitive variables (T,vx,vy,vz,P,rho). */
  double *Limiter_Primitive;    /*!< \brief Limiter of the primitive variables (T,vx,vy,vz,P,rho). */
  
public:
  
	/*!
	 * \brief Constructor of the class.
	 */
	CEulerVariable(void);
  
	/*!
	 * \overload
	 * \param[in] val_density - Value of the flow density (initialization value).
	 * \param[in] val_velocity - Value of the flow velocity (initialization value).
	 * \param[in] val_energy - Value of the flow energy (initialization value).
	 * \param[in] val_ndim - Number of dimensions of the problem.
	 * \param[in] val_nvar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CEulerVariable(double val_density, double *val_velocity, double val_energy, unsigned short val_ndim,
                 unsigned short val_nvar, CConfig *config);
  
	/*!
	 * \overload
	 * \param[in] val_solution - Pointer to the flow value (initialization value).
	 * \param[in] val_ndim - Number of dimensions of the problem.
	 * \param[in] val_nvar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CEulerVariable(double *val_solution, unsigned short val_ndim, unsigned short val_nvar, CConfig *config);
  
	/*!
	 * \brief Destructor of the class.
	 */
	virtual ~CEulerVariable(void);
  
	/*!
	 * \brief Set to zero the gradient of the primitive variables.
	 */
	void SetGradient_PrimitiveZero(unsigned short val_primvar);
  
	/*!
	 * \brief Add <i>val_value</i> to the gradient of the primitive variables.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_dim - Index of the dimension.
	 * \param[in] val_value - Value to add to the gradient of the primitive variables.
	 */
	void AddGradient_Primitive(unsigned short val_var, unsigned short val_dim, double val_value);
  
	/*!
	 * \brief Subtract <i>val_value</i> to the gradient of the primitive variables.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_dim - Index of the dimension.
	 * \param[in] val_value - Value to subtract to the gradient of the primitive variables.
	 */
	void SubtractGradient_Primitive(unsigned short val_var, unsigned short val_dim, double val_value);
  
	/*!
	 * \brief Get the value of the primitive variables gradient.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_dim - Index of the dimension.
	 * \return Value of the primitive variables gradient.
	 */
	double GetGradient_Primitive(unsigned short val_var, unsigned short val_dim);
  
  /*!
	 * \brief Get the value of the primitive variables gradient.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_dim - Index of the dimension.
	 * \return Value of the primitive variables gradient.
	 */
	double GetLimiter_Primitive(unsigned short val_var);
  
	/*!
	 * \brief Set the gradient of the primitive variables.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_dim - Index of the dimension.
	 * \param[in] val_value - Value of the gradient.
	 */
	void SetGradient_Primitive(unsigned short val_var, unsigned short val_dim, double val_value);
  
  /*!
	 * \brief Set the gradient of the primitive variables.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_dim - Index of the dimension.
	 * \param[in] val_value - Value of the gradient.
	 */
	void SetLimiter_Primitive(unsigned short val_var, double val_value);
  
	/*!
	 * \brief Get the value of the primitive variables gradient.
	 * \return Value of the primitive variables gradient.
	 */
	double **GetGradient_Primitive(void);
  
  /*!
	 * \brief Get the value of the primitive variables gradient.
	 * \return Value of the primitive variables gradient.
	 */
	double *GetLimiter_Primitive(void);
  
	/*!
	 * \brief Set the value of the pressure.
	 */
	bool SetPressure(double Gamma);
  
	/*!
	 * \brief Set the value of the speed of the sound.
	 * \param[in] Gamma - Value of Gamma.
	 */
	bool SetSoundSpeed(double Gamma);
  
	/*!
	 * \brief Set the value of the enthalpy.
	 */
	void SetEnthalpy(void);
	
	/*!
	 * \brief Set all the primitive variables for compressible flows.
	 */
	bool SetPrimVar_Compressible(CConfig *config);
	
	/*!
	 * \brief Get the primitive variables.
	 * \param[in] val_var - Index of the variable.
	 * \return Value of the primitive variable for the index <i>val_var</i>.
	 */
	double GetPrimVar(unsigned short val_var);
  
  /*!
	 * \brief Set the value of the primitive variables.
	 * \param[in] val_var - Index of the variable.
   * \param[in] val_var - Index of the variable.
	 * \return Set the value of the primitive variable for the index <i>val_var</i>.
	 */
	void SetPrimVar(unsigned short val_var, double val_prim);
  
  /*!
	 * \brief Set the value of the primitive variables.
	 * \param[in] val_prim - Primitive variables.
	 * \return Set the value of the primitive variable for the index <i>val_var</i>.
	 */
	void SetPrimVar(double *val_prim);
  
	/*!
	 * \brief Get the primitive variables of the problem.
	 * \return Pointer to the primitive variable vector.
	 */
	double *GetPrimVar(void);
  
  /*!
	 * \brief Set the value of the density for the incompressible flows.
	 */
	bool SetDensity(void);
  
	/*!
	 * \brief Set the value of the temperature.
	 * \param[in] Gas_Constant - Value of Gas Constant
	 */
	bool SetTemperature(double Gas_Constant);
  
	/*!
	 * \brief Get the norm 2 of the velocity.
	 * \return Norm 2 of the velocity vector.
	 */
	double GetVelocity2(void);
  
	/*!
	 * \brief Get the flow pressure.
	 * \return Value of the flow pressure.
	 */
	double GetPressure(void);

	/*!
	 * \brief Get the speed of the sound.
	 * \return Value of speed of the sound.
	 */
	double GetSoundSpeed(void);
  
	/*!
	 * \brief Get the enthalpy of the flow.
	 * \return Value of the enthalpy of the flow.
	 */
	double GetEnthalpy(void);
  
	/*!
	 * \brief Get the density of the flow.
	 * \return Value of the density of the flow.
	 */
	double GetDensity(void);
  
	/*!
	 * \brief Get the energy of the flow.
	 * \return Value of the energy of the flow.
	 */
	double GetEnergy(void);
  
	/*!
	 * \brief Get the temperature of the flow.
	 * \return Value of the temperature of the flow.
	 */
	double GetTemperature(void);
  
	/*!
	 * \brief Get the velocity of the flow.
	 * \param[in] val_dim - Index of the dimension.
	 * \return Value of the velocity for the dimension <i>val_dim</i>.
	 */
	double GetVelocity(unsigned short val_dim);
  
	/*!
	 * \brief Get the projected velocity in a unitary vector direction (compressible solver).
	 * \param[in] val_vector - Direction of projection.
	 * \return Value of the projected velocity.
	 */
	double GetProjVel(double *val_vector);
  
	/*!
	 * \brief Set the velocity vector from the solution.
	 * \param[in] val_velocity - Pointer to the velocity.
	 */
	void SetVelocity(void);
  
	/*!
	 * \brief Set the velocity vector from the old solution.
	 * \param[in] val_velocity - Pointer to the velocity.
	 */
	void SetVelocity_Old(double *val_velocity);
  
	/*!
	 * \brief Get the value of the preconditioner Beta.
	 * \return Value of the low Mach preconditioner variable Beta
	 */
	double GetPreconditioner_Beta();
  
	/*!
	 * \brief Set the value of the preconditioner Beta.
	 * \param[in] Value of the low Mach preconditioner variable Beta
	 */
	void SetPreconditioner_Beta(double val_Beta);
};

/*!
 * \class CNSVariable
 * \brief Main class for defining the variables of the Navier-Stokes' solver.
 * \ingroup Navier_Stokes_Equations
 * \author F. Palacios.
 * \version 1.0.0
 */
class CNSVariable : public CEulerVariable {
private:
	double Prandtl_Lam;       /*!< \brief Laminar Prandtl number. */
	double Prandtl_Turb;      /*!< \brief Turbulent Prandtl number. */
	double Temperature_Ref;   /*!< \brief Reference temperature of the fluid. */
	double Viscosity_Ref;     /*!< \brief Reference viscosity of the fluid. */
	double Viscosity_Inf;     /*!< \brief Viscosity of the fluid at the infinity. */
	double Vorticity[3];		/*!< \brief Vorticity of the fluid. */
	double StrainMag;           /*!< \brief Magnitude of rate of strain tensor. */
public:
  
	/*!
	 * \brief Constructor of the class.
	 */
	CNSVariable(void);
  
	/*!
	 * \overload
	 * \param[in] val_density - Value of the flow density (initialization value).
	 * \param[in] val_velocity - Value of the flow velocity (initialization value).
	 * \param[in] val_energy - Value of the flow energy (initialization value).
	 * \param[in] val_ndim - Number of dimensions of the problem.
	 * \param[in] val_nvar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CNSVariable(double val_density, double *val_velocity,
              double val_energy, unsigned short val_ndim, unsigned short val_nvar, CConfig *config);
  
	/*!
	 * \overload
	 * \param[in] val_solution - Pointer to the flow value (initialization value).
	 * \param[in] val_ndim - Number of dimensions of the problem.
	 * \param[in] val_nvar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CNSVariable(double *val_solution, unsigned short val_ndim, unsigned short val_nvar, CConfig *config);
  
	/*!
	 * \brief Destructor of the class.
	 */
	~CNSVariable(void);
  
	/*!
	 * \brief Set the laminar viscosity.
	 */
	void SetLaminarViscosity(CConfig *config);
  
	/*!
	 * \brief Set the vorticity value.
	 */
	void SetVorticity(void);
  
	/*!
	 * \brief Set the rate of strain magnitude.
	 */
	void SetStrainMag(void);
  
	/*!
	 * \overload
	 * \param[in] val_eddy_viscosity - Value of the eddy viscosity.
	 */
	void SetEddyViscosity(double eddy_visc);
  
	/*!
	 * \brief Get the laminar viscosity of the flow.
	 * \return Value of the laminar viscosity of the flow.
	 */
	double GetLaminarViscosity(void);
  
	/*!
	 * \brief Get the eddy viscosity of the flow.
	 * \return The eddy viscosity of the flow.
	 */
	double GetEddyViscosity(void);
  
	/*!
	 * \brief Get the value of the vorticity.
	 * \param[in] val_dim - Index of the dimension.
	 * \return Value of the vorticity.
	 */
	double GetVorticity(unsigned short val_dim);
  
	/*!
	 * \brief Get the value of the magnitude of rate of strain.
	 * \return Value of the rate of strain magnitude.
	 */
	double GetStrainMag(void);
  
	/*!
	 * \brief Set the value of pressure.
	 */
	bool SetPressure(double Gamma, double turb_ke);
	
	/*!
	 * \brief Set all the primitive variables for compressible flows
	 */
	bool SetPrimVar_Compressible(double eddy_visc, double turb_ke, CConfig *config);

};

/*!
 * \class CTurbVariable
 * \brief Main class for defining the variables of the turbulence model.
 * \ingroup Turbulence_Model
 * \author A. Bueno.
 * \version 1.0.0
 */
class CTurbVariable : public CVariable {
protected:
	double muT;                /*!< \brief Eddy viscosity. */
  
public:
	/*!
	 * \brief Constructor of the class.
	 */
	CTurbVariable(void);
  
	/*!
	 * \overload
	 * \param[in] val_ndim - Number of dimensions of the problem.
	 * \param[in] val_nvar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CTurbVariable(unsigned short val_ndim, unsigned short val_nvar, CConfig *config);
  
	/*!
	 * \brief Destructor of the class.
	 */
	virtual ~CTurbVariable(void);
  
	/*!
	 * \brief Get the value of the eddy viscosity.
	 * \return the value of the eddy viscosity.
	 */
	double GetmuT();
  
	/*!
	 * \brief Set the value of the eddy viscosity.
	 * \param[in] val_muT - Value of the eddy viscosity.
	 */
	void SetmuT(double val_muT);
};

/*!
 * \class CTurbSAVariable
 * \brief Main class for defining the variables of the turbulence model.
 * \ingroup Turbulence_Model
 * \author A. Bueno.
 * \version 1.0.0
 */

class CTurbSAVariable : public CTurbVariable {
public:
	/*!
	 * \brief Constructor of the class.
	 */
	CTurbSAVariable(void);
  
	/*!
	 * \overload
	 * \param[in] val_nu_tilde - Turbulent variable value (initialization value).
	 * \param[in] val_muT  - The eddy viscosity
	 * \param[in] val_ndim - Number of dimensions of the problem.
	 * \param[in] val_nvar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CTurbSAVariable(double val_nu_tilde, double val_muT, unsigned short val_ndim, unsigned short val_nvar, CConfig *config);
  
	/*!
	 * \brief Destructor of the class.
	 */
	~CTurbSAVariable(void);
  
};

/*!
 * \class CTurbSSTVariable
 * \brief Main class for defining the variables of the turbulence model.
 * \ingroup Turbulence_Model
 * \author A. Bueno.
 * \version 1.0.0
 */

class CTurbSSTVariable : public CTurbVariable {
protected:
	double sigma_om2,
	beta_star;
	double F1,		/*!< \brief Menter blending function for blending of k-w and k-eps. */
	F2,		        /*!< \brief Menter blending function for stress limiter. */
	CDkw;           /*!< \brief Cross-diffusion. */
  
public:
	/*!
	 * \brief Constructor of the class.
	 */
	CTurbSSTVariable(void);
  
	/*!
	 * \overload
	 * \param[in] val_rho_kine - Turbulent variable value (initialization value).
	 * \param[in] val_rho_omega - Turbulent variable value (initialization value).
	 * \param[in] val_ndim - Number of dimensions of the problem.
	 * \param[in] val_nvar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CTurbSSTVariable(double val_rho_kine, double val_rho_omega, double val_muT, unsigned short val_ndim, unsigned short val_nvar,
                   double *constants, CConfig *config);
  
	/*!
	 * \brief Destructor of the class.
	 */
	~CTurbSSTVariable(void);
  
	/*!
	 * \brief Set the blending function for the blending of k-w and k-eps.
	 * \param[in] val_viscosity - Value of the vicosity.
	 * \param[in] val_dist - Value of the distance to the wall.
	 * \param[in] val_density - Value of the density.
	 */
	void SetBlendingFunc(double val_viscosity, double val_dist, double val_density);
  
	/*!
	 * \brief Get the first blending function.
	 */
	double GetF1blending(void);
  
	/*!
	 * \brief Get the second blending function.
	 */
	double GetF2blending(void);
  
	/*!
	 * \brief Get the value of the cross diffusion of tke and omega.
	 */
	double GetCrossDiff(void);
};

#include "variable_structure.inl"
