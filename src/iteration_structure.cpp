/*!
 * \file iteration_structure.cpp
 * \brief Main subroutines used by SU2_CFD.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 1.0.0
 *
 * Stanford University Unstructured (SU2).
 * Copyright (C) 2012-2013 Aerospace Design Laboratory (ADL).
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

#include "../include/iteration_structure.hpp"

void MeanFlowIteration(COutput *output, CIntegration **integration_container, CGeometry **geometry_container,
                       CSolver ***solver_container, CNumerics ****numerics_container, CConfig *config_container) {
  
  unsigned long IntIter = 0; config_container->SetIntIter(IntIter);
  unsigned long ExtIter = config_container->GetExtIter();
  
#ifndef NO_MPI
	int rank = MPI::COMM_WORLD.Get_rank();
#endif
  
  /*--- Set the value of the internal iteration ---*/
  
  IntIter = ExtIter;
  
  /*--- Set the initial condition ---*/
  
  solver_container[MESH_0][FLOW_SOL]->SetInitialCondition(geometry_container, solver_container, config_container, ExtIter);
  
  /*--- Update global parameters ---*/
  
  if (config_container->GetKind_Solver() == EULER) { config_container->SetGlobalParam(EULER, RUNTIME_FLOW_SYS, ExtIter); }
  if (config_container->GetKind_Solver() == NAVIER_STOKES) { config_container->SetGlobalParam(NAVIER_STOKES, RUNTIME_FLOW_SYS, ExtIter); }
  if (config_container->GetKind_Solver() == RANS) { config_container->SetGlobalParam(RANS, RUNTIME_FLOW_SYS, ExtIter); }
  
  /*--- Solve the Euler, Navier-Stokes or Reynolds-averaged Navier-Stokes (RANS) equations (one iteration) ---*/
  
  integration_container[FLOW_SOL]->MultiGrid_Iteration(geometry_container, solver_container, numerics_container,
                                                               config_container, RUNTIME_FLOW_SYS, IntIter);
  
  /*--- Solve the turbulence model ---*/
  
  if (config_container->GetKind_Solver() == RANS) {
    
    config_container->SetGlobalParam(RANS, RUNTIME_TURB_SYS, ExtIter);
    integration_container[TURB_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container,
                                                                  config_container, RUNTIME_TURB_SYS, IntIter);
    
  }
  
}
