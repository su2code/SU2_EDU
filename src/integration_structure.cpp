/*!
 * \file integration_structure.cpp
 * \brief This subroutine includes the space and time integration structure.
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

#include "../include/integration_structure.hpp"

CIntegration::CIntegration(CConfig *config) {
	Cauchy_Value = 0;
	Cauchy_Func = 0;
	Old_Func = 0;
	New_Func = 0;
	Cauchy_Counter = 0;
	Convergence = false;
	Convergence_OneShot = false;
	Convergence_FullMG = false;
	Cauchy_Serie = new double [config->GetCauchy_Elems()+1];	
}

CIntegration::~CIntegration(void) {
	delete [] Cauchy_Serie;
}

void CIntegration::Space_Integration(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics,
                                     CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem) {
	unsigned short iMarker;
    
	unsigned short MainSolver = config->GetContainerPosition(RunTime_EqSystem);
    
	/*--- Compute inviscid residuals ---*/
	switch (config->GetKind_ConvNumScheme()) {
		case SPACE_CENTERED:
			solver_container[MainSolver]->Centered_Residual(geometry, solver_container, numerics[CONV_TERM], config, iMesh, iRKStep);
			break;
		case SPACE_UPWIND:
			solver_container[MainSolver]->Upwind_Residual(geometry, solver_container, numerics[CONV_TERM], config, iMesh);
			break;
	}
    
	/*--- Compute viscous residuals ---*/
	switch (config->GetKind_ViscNumScheme()) {
        case AVG_GRAD: case AVG_GRAD_CORRECTED:
            solver_container[MainSolver]->Viscous_Residual(geometry, solver_container, numerics[VISC_TERM], config, iMesh, iRKStep);
            break;
        case GALERKIN:
            solver_container[MainSolver]->Galerkin_Method(geometry, solver_container, numerics[VISC_TERM], config, iMesh);
            break;
	}
    
	/*--- Compute source term residuals ---*/
	switch (config->GetKind_SourNumScheme()) {
        case PIECEWISE_CONSTANT:
            solver_container[MainSolver]->Source_Residual(geometry, solver_container, numerics[SOURCE_FIRST_TERM], numerics[SOURCE_SECOND_TERM], config, iMesh);
            break;
	}
  
	/*--- Weak boundary conditions ---*/
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
		switch (config->GetMarker_All_Boundary(iMarker)) {
			case EULER_WALL:
				solver_container[MainSolver]->BC_Euler_Wall(geometry, solver_container, numerics[CONV_BOUND_TERM], config, iMarker);
				break;
			case INLET_FLOW:
				solver_container[MainSolver]->BC_Inlet(geometry, solver_container, numerics[CONV_BOUND_TERM], numerics[VISC_BOUND_TERM], config, iMarker);
				break;
      case SUPERSONIC_INLET:
				solver_container[MainSolver]->BC_Supersonic_Inlet(geometry, solver_container, numerics[CONV_BOUND_TERM], numerics[VISC_BOUND_TERM], config, iMarker);
        break;
			case OUTLET_FLOW:
				solver_container[MainSolver]->BC_Outlet(geometry, solver_container, numerics[CONV_BOUND_TERM], numerics[VISC_BOUND_TERM], config, iMarker);
				break;
			case FAR_FIELD:
				solver_container[MainSolver]->BC_Far_Field(geometry, solver_container, numerics[CONV_BOUND_TERM], numerics[VISC_BOUND_TERM], config, iMarker);
				break;
			case SYMMETRY_PLANE:
				solver_container[MainSolver]->BC_Sym_Plane(geometry, solver_container, numerics[CONV_BOUND_TERM], numerics[VISC_BOUND_TERM], config, iMarker);
				break;
		}
	}
    
	/*--- Strong boundary conditions (Navier-Stokes and Dirichlet type BCs) ---*/
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		switch (config->GetMarker_All_Boundary(iMarker)) {
      case ISOTHERMAL:
        solver_container[MainSolver]->BC_Isothermal_Wall(geometry, solver_container, numerics[CONV_BOUND_TERM], numerics[VISC_BOUND_TERM], config, iMarker);
				break;
      case HEAT_FLUX:
        solver_container[MainSolver]->BC_HeatFlux_Wall(geometry, solver_container, numerics[CONV_BOUND_TERM], numerics[VISC_BOUND_TERM], config, iMarker);
				break;
		}
  
}

void CIntegration::Time_Integration(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iRKStep,
                                    unsigned short RunTime_EqSystem, unsigned long Iteration) {
	unsigned short MainSolver = config->GetContainerPosition(RunTime_EqSystem);
    
  /*--- Perform the time integration ---*/
  switch (config->GetKind_TimeIntScheme()) {
    case (RUNGE_KUTTA_EXPLICIT):
      solver_container[MainSolver]->ExplicitRK_Iteration(geometry, solver_container, config, iRKStep);
      break;
    case (EULER_EXPLICIT):
      solver_container[MainSolver]->ExplicitEuler_Iteration(geometry, solver_container, config);
      break;
    case (EULER_IMPLICIT):
      solver_container[MainSolver]->ImplicitEuler_Iteration(geometry, solver_container, config);
      break;
  }
  
}

void CIntegration::Convergence_Monitoring(CGeometry *geometry, CConfig *config, unsigned long Iteration, double monitor) {

  unsigned short iCounter;

#ifndef NO_MPI
  int size = MPI::COMM_WORLD.Get_size();
  int rank = MPI::COMM_WORLD.Get_rank();
#endif

	bool Already_Converged = Convergence;
	
  /*--- Cauchi based convergence criteria ---*/
	if (config->GetConvCriteria() == CAUCHY) {
    
    /*--- Initialize at the fist iteration ---*/
		if (Iteration  == 0) {
			Cauchy_Value = 0.0;
			Cauchy_Counter = 0.0;
			for (iCounter = 0; iCounter < config->GetCauchy_Elems(); iCounter++)
				Cauchy_Serie[iCounter] = 0.0;
		}

		Old_Func = New_Func;
		New_Func = monitor;
		Cauchy_Func = fabs(New_Func - Old_Func);

		Cauchy_Serie[Cauchy_Counter] = Cauchy_Func;
		Cauchy_Counter++;

		if (Cauchy_Counter == config->GetCauchy_Elems()) Cauchy_Counter = 0;

		Cauchy_Value = 1;
		if (Iteration  >= config->GetCauchy_Elems()) {
			Cauchy_Value = 0;
			for (iCounter = 0; iCounter < config->GetCauchy_Elems(); iCounter++)
				Cauchy_Value += Cauchy_Serie[iCounter];
		}

		if (Cauchy_Value >= config->GetCauchy_Eps()) Convergence = false;
		else Convergence = true;

		if (Cauchy_Value >= config->GetCauchy_Eps_OneShot()) Convergence_OneShot = false;
		else Convergence_OneShot = true;

		if (Cauchy_Value >= config->GetCauchy_Eps_FullMG()) Convergence_FullMG = false;
		else Convergence_FullMG = true;
	}

  /*--- Residual based convergence criteria ---*/
  if (config->GetConvCriteria() == RESIDUAL) {
    
    /*--- Compute the initial value ---*/
    if (Iteration == config->GetStartConv_Iter() ) InitResidual = monitor;
    if (monitor > InitResidual) InitResidual = monitor;

    /*--- Check the convergence ---*/
    if (((fabs(InitResidual - monitor) >= config->GetOrderMagResidual()) && (monitor < InitResidual))  ||
        (monitor <= config->GetMinLogResidual())) Convergence = true;
    else Convergence = false;
    
  }
  
  /*--- Do not apply any convergence criteria of the number 
   of iterations is less than a particular value ---*/
	if (Iteration < config->GetStartConv_Iter()) {
		Convergence = false;
		Convergence_OneShot = false;
		Convergence_FullMG = false;
	}

	if (Already_Converged) Convergence = true;

  
  /*--- Apply the same convergence criteria to all the processors ---*/
#ifndef NO_MPI

  unsigned short *sbuf_conv = NULL, *rbuf_conv = NULL;
  sbuf_conv = new unsigned short[1]; sbuf_conv[0] = 0;
  rbuf_conv = new unsigned short[1]; rbuf_conv[0] = 0;

  /*--- Convergence criteria ---*/
  sbuf_conv[0] = Convergence;
  MPI::COMM_WORLD.Reduce(sbuf_conv, rbuf_conv, 1, MPI::UNSIGNED_SHORT, MPI::SUM, MASTER_NODE);
  MPI::COMM_WORLD.Barrier();

  /*-- Compute global convergence criteria in the master node --*/
  sbuf_conv[0] = 0;
  if (rank == MASTER_NODE) {
    if (rbuf_conv[0] == size) sbuf_conv[0] = 1;
    else sbuf_conv[0] = 0;
  }
  
  MPI::COMM_WORLD.Bcast(sbuf_conv, 1, MPI::UNSIGNED_SHORT, MASTER_NODE);

  if (sbuf_conv[0] == 1) Convergence = true;
  else Convergence = false;
  
  delete [] sbuf_conv;
  delete [] rbuf_conv;

#endif
  
	/*--- Stop the simulation in case a nan appears, do not save the solution ---*/
	if (monitor != monitor) {

#ifdef NO_MPI
    cout << "\n !!! Error: NaNs detected in solution. Now exiting... !!!" << endl;
		exit(1);
#else
    if (rank == MASTER_NODE)
      cout << "\n !!! Error: NaNs detected in solution. Now exiting... !!!" << endl;
    MPI::COMM_WORLD.Barrier();
    MPI::COMM_WORLD.Abort(1);
#endif
        
	}
  
#ifndef NO_MPI

  MPI::COMM_WORLD.Barrier();

#endif
  
}
