/*!
 * \file iteration_structure.cpp
 * \brief Main subroutines used by SU2_EDU.
 * \author Aerospace Design Laboratory (Stanford University).
 * \version 1.2.0
 *
 * SU2 EDU, Copyright (C) 2014 Aerospace Design Laboratory (Stanford University).
 *
 * SU2 EDU is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 EDU is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2 EDU. If not, see <http://www.gnu.org/licenses/>.
 */

#include "../include/iteration_structure.hpp"

void MeanFlowIteration(COutput *output, CIntegration **integration_container, CGeometry **geometry_container,
                       CSolver ***solver_container, CNumerics ****numerics_container, CConfig *config_container) {
  
	double Physical_dt, Physical_t;
	unsigned short iMesh;
  
  unsigned long IntIter = 0; config_container->SetIntIter(IntIter);
  unsigned long ExtIter = config_container->GetExtIter();
  
  /*--- Set the value of the internal iteration ---*/
  
  IntIter = ExtIter;
  if ((config_container->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
      (config_container->GetUnsteady_Simulation() == DT_STEPPING_2ND)) IntIter = 0;
  
  /*--- Set the initial condition ---*/
  
  solver_container[MESH_0][FLOW_SOL]->SetInitialCondition(geometry_container, solver_container, config_container, ExtIter);
  
  /*--- Update global parameters ---*/
  
  if (config_container->GetKind_Solver() == EULER){
    config_container->SetGlobalParam(EULER, RUNTIME_FLOW_SYS, ExtIter);
  }
  if (config_container->GetKind_Solver() == NAVIER_STOKES){
    config_container->SetGlobalParam(NAVIER_STOKES, RUNTIME_FLOW_SYS, ExtIter);
  }
  if (config_container->GetKind_Solver() == RANS){
    config_container->SetGlobalParam(RANS, RUNTIME_FLOW_SYS, ExtIter);
  }
  
  /*--- Solve the Euler, Navier-Stokes or Reynolds-averaged Navier-Stokes (RANS) equations (one iteration) ---*/
  
  integration_container[FLOW_SOL]->MultiGrid_Iteration(geometry_container, solver_container, numerics_container,
                                                              config_container, RUNTIME_FLOW_SYS, IntIter);
  
  if (config_container->GetKind_Solver() == RANS) {
    
    /*--- Solve the turbulence model ---*/
    
    config_container->SetGlobalParam(RANS, RUNTIME_TURB_SYS, ExtIter);
    integration_container[TURB_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container,
                                                                 config_container, RUNTIME_TURB_SYS, IntIter);
    
  }
  
	/*--- Dual time stepping strategy ---*/
  
	if ((config_container->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
      (config_container->GetUnsteady_Simulation() == DT_STEPPING_2ND)) {
    
		for(IntIter = 1; IntIter < config_container->GetUnst_nIntIter(); IntIter++) {
      
      /*--- Write the convergence history (only screen output) ---*/
      
      output->SetConvergence_History(NULL, geometry_container, solver_container, config_container, integration_container, true, 0.0);
      
      /*--- Set the value of the internal iteration ---*/
      
      config_container->SetIntIter(IntIter);
      
      /*--- Pseudo-timestepping for the Euler, Navier-Stokes or Reynolds-averaged Navier-Stokes equations ---*/
      
      if (config_container->GetKind_Solver() == EULER)
        config_container->SetGlobalParam(EULER, RUNTIME_FLOW_SYS, ExtIter);
      if (config_container->GetKind_Solver() == NAVIER_STOKES)
        config_container->SetGlobalParam(NAVIER_STOKES, RUNTIME_FLOW_SYS, ExtIter);
      if (config_container->GetKind_Solver() == RANS)
        config_container->SetGlobalParam(RANS, RUNTIME_FLOW_SYS, ExtIter);
      
      /*--- Solve the Euler, Navier-Stokes or Reynolds-averaged Navier-Stokes (RANS) equations (one iteration) ---*/
      
      integration_container[FLOW_SOL]->MultiGrid_Iteration(geometry_container, solver_container, numerics_container,
                                                           config_container, RUNTIME_FLOW_SYS, IntIter);
      
      /*--- Pseudo-timestepping the turbulence model ---*/
      
      if (config_container->GetKind_Solver() == RANS) {
        
        /*--- Solve the turbulence model ---*/
        
        config_container->SetGlobalParam(RANS, RUNTIME_TURB_SYS, ExtIter);
        integration_container[TURB_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container,
                                                              config_container, RUNTIME_TURB_SYS, IntIter);
      }
      
      if (integration_container[FLOW_SOL]->GetConvergence()) break;
      
    }
    
    /*--- Update dual time solver on all mesh levels ---*/
    
    for (iMesh = 0; iMesh <= config_container->GetMGLevels(); iMesh++) {
      integration_container[FLOW_SOL]->SetDualTime_Solver(geometry_container[iMesh], solver_container[iMesh][FLOW_SOL], config_container);
      integration_container[FLOW_SOL]->SetConvergence(false);
    }
    
    /*--- Update dual time solver for the turbulence model ---*/
    
    if (config_container->GetKind_Solver() == RANS) {
      integration_container[TURB_SOL]->SetDualTime_Solver(geometry_container[MESH_0], solver_container[MESH_0][TURB_SOL], config_container);
      integration_container[TURB_SOL]->SetConvergence(false);
    }
    
    /*--- Verify convergence criteria (based on total time) ---*/
    
    Physical_dt = config_container->GetDelta_UnstTime();
    Physical_t  = (ExtIter+1)*Physical_dt;
    if (Physical_t >=  config_container->GetTotal_UnstTime())
      integration_container[FLOW_SOL]->SetConvergence(true);
    
  }
    
}

