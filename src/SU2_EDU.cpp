/*!
 * \file SU2_EDU.cpp
 * \brief Main file of Computational Fluid Dynamics Code (SU2_EDU).
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

#include "../include/SU2_EDU.hpp"

using namespace std;

int main(int argc, char *argv[]) {
  
  bool StopCalc = false;
  unsigned long ExtIter = 0;
  double StartTime = 0.0, StopTime = 0.0, UsedTime = 0.0;
  unsigned short iMesh, iSol, nDim;
  ofstream ConvHist_file;
  int rank = MASTER_NODE;
  int size = SINGLE_NODE;
  
#ifndef NO_MPI
  /*--- MPI initialization, and buffer setting ---*/
  MPI::Init(argc, argv);
  rank = MPI::COMM_WORLD.Get_rank();
  size = MPI::COMM_WORLD.Get_size();
#endif
  
  /*--- Create pointers to all of the classes that may be used throughout
   the SU2_EDU code. In general, the pointers are instantiated down a
   heirarchy over all zones, multigrid levels, equation sets, and equation
   terms as described in the comments below. ---*/
  
  COutput *output                      = NULL;
  CIntegration **integration_container = NULL;
  CGeometry **geometry_container       = NULL;
  CSolver ***solver_container          = NULL;
  CNumerics ****numerics_container     = NULL;
  CConfig *config_container            = NULL;
  CSurfaceMovement *surface_movement   = NULL;
  CVolumetricMovement *grid_movement   = NULL;
  
  /*--- Load in the number of zones and spatial dimensions in the mesh file (If no config
   file is specified, default.cfg is used) ---*/
  
  cout << endl <<"-------------------------------------------------------------------------" << endl;
  cout <<"|    _____   _    _   ___                                               |" << endl;
  cout <<"|   / ____| | |  | | |__ \\    Web: su2.stanford.edu                     |" << endl;
  cout <<"|  | (___   | |  | |    ) |   Twitter: @su2code                         |" << endl;
  cout <<"|   \\___ \\  | |  | |   / /    Forum: www.cfd-online.com/Forums/su2/     |" << endl;
  cout <<"|   ____) | | |__| |  / /_                                              |" << endl;
  cout <<"|  |_____/   \\____/  |____|   Educational Code                          |" << endl;
  cout <<"|                             Release 1.0.0                             |" << endl;
  cout <<"-------------------------------------------------------------------------" << endl;
  cout << "| SU2, Copyright (C) 2012-2014 Aerospace Design Laboratory (ADL).       |" << endl;
  cout << "| SU2 is distributed in the hope that it will be useful,                |" << endl;
  cout << "| but WITHOUT ANY WARRANTY; without even the implied warranty of        |" << endl;
  cout << "| MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU      |" << endl;
  cout << "| Lesser General Public License (version 2.1) for more details.         |" << endl;
  cout <<"-------------------------------------------------------------------------" << endl;
  cout << endl;
  cout <<"-------------------------------------------------------------------------" << endl;
  cout << "| Please check that the SU2_EDU binary, config files, and meshes are    |" << endl;
  cout << "| all located in the same directory (SU2_EDU/bin/ by default).          |" << endl;
  cout << "| For detailed use instructions, see the SU2_EDU/README file.           |" << endl;
  cout <<"-------------------------------------------------------------------------" << endl;
  
  /*--- Get user input for the simulation type (viscous/inviscid) ---*/
  
  char config_file_name[200];
  int SimType = 0;
  string Input = "";
  while(1) {
    
    cout << endl;
    cout << "   [0] Inviscid (Euler)" << endl;
    cout << "   [1] Viscous (RANS)"  << endl;
    cout << "Select simulation type [0]: " ;
    getline(cin, Input);
    stringstream myStream(Input);
    
    /*-- Handle default option --*/
    if (Input.empty())
      myStream << "0";
    
    /*-- Check for valid input --*/
    if (myStream >> SimType) {
      if (SimType == 0) {
        strcpy(config_file_name, "Inviscid_ConfigFile.cfg");
        break;
      } else if(SimType == 1) {
        strcpy(config_file_name, "Viscous_ConfigFile.cfg");
        break;
      }
    }
  }
  
  /*--- The initial version of SU2_EDU assumes 2-D in a single mesh zone. ---*/
  
  nDim  = 2;
  
  /*--- Definition of the configuration option class. In this constructor, 
   the input configuration file is parsed and all options are read and stored. ---*/
  
  config_container = new CConfig(config_file_name, SU2_EDU, VERB_HIGH);
  
#ifndef NO_MPI
  /*--- Change the name of the input/output files for a parallel computation ---*/
  config_container->SetFileNameDomain(rank+1);
#endif
  
  /*--- Perform the non-dimensionalization for the flow equations using the
   specified reference values. ---*/
  
  config_container->SetNondimensionalization(nDim);
  
  /*--- Definition of the geometry class. Within this constructor, the
   mesh file is read and the primal grid is stored (node coords, connectivity,
   & boundary markers). MESH_0 is the index of the original (finest) mesh in
   the multigrid structure which is stored in the CPhysicalGeometry class. ---*/
  
  geometry_container = new CGeometry*[config_container->GetMGLevels()+1];
  geometry_container[MESH_0] = new CPhysicalGeometry(config_container);
  
  if (rank == MASTER_NODE)
    cout << endl <<"------------------------- Geometry Preprocessing ------------------------" << endl;
  
  /*--- Preprocessing of the geometry. In this routine, the edge-
   based data structure is constructed, i.e., node and cell neighbors are
   identified and linked, face areas and volumes of the dual mesh cells are
   computed, and the multigrid levels are created using an agglomeration procedure. ---*/
  
  Geometrical_Preprocessing(geometry_container, config_container);
  
#ifndef NO_MPI
  /*--- Synchronization point after the geometrical definition subroutine ---*/
  MPI::COMM_WORLD.Barrier();
#endif
  
  if (rank == MASTER_NODE)
    cout << endl <<"------------------------- Solver Preprocessing --------------------------" << endl;
  
  /*--- Definition of the solver class: solver_container[#MG_GRIDS][#EQ_SYSTEMS].
   The solver classes are specific to a particular set of governing equations,
   and they contain the subroutines with instructions for computing each spatial
   term of the PDE, i.e., loops over the edges to compute convective and viscous
   fluxes, loops over the nodes to compute source terms, and routines for
   imposing various boundary condition type for the PDE. ---*/
  
  solver_container = new CSolver**[config_container->GetMGLevels()+1];
  for (iMesh = 0; iMesh <= config_container->GetMGLevels(); iMesh++)
    solver_container[iMesh] = NULL;
  
  for (iMesh = 0; iMesh <= config_container->GetMGLevels(); iMesh++) {
    solver_container[iMesh] = new CSolver*[MAX_SOLS];
    for (iSol = 0; iSol < MAX_SOLS; iSol++)
      solver_container[iMesh][iSol] = NULL;
  }
  
  Solver_Preprocessing(solver_container, geometry_container, config_container);
  
#ifndef NO_MPI
  /*--- Synchronization point after the solution preprocessing subroutine ---*/
  MPI::COMM_WORLD.Barrier();
#endif
  
  /*--- Definition of the integration class: integration_container[#EQ_SYSTEMS].
   The integration class orchestrates the execution of the spatial integration
   subroutines contained in the solver class (including multigrid) for computing
   the residual at each node, R(U) and then integrates the equations to a
   steady state or time-accurately. ---*/
  
  integration_container = new CIntegration*[MAX_SOLS];
  Integration_Preprocessing(integration_container, geometry_container, config_container);
  
#ifndef NO_MPI
  /*--- Synchronization point after the integration definition subroutine ---*/
  MPI::COMM_WORLD.Barrier();
#endif
  
  /*--- Definition of the numerical method class:
   numerics_container[#MG_GRIDS][#EQ_SYSTEMS][#EQ_TERMS].
   The numerics class contains the implementation of the numerical methods for
   evaluating convective or viscous fluxes between any two nodes in the edge-based
   data structure (centered, upwind, galerkin), as well as any source terms
   (piecewise constant reconstruction) evaluated in each dual mesh volume. ---*/
  
  numerics_container = new CNumerics***[config_container->GetMGLevels()+1];
  Numerics_Preprocessing(numerics_container, solver_container,geometry_container, config_container);
  
#ifndef NO_MPI
  /*--- Synchronization point after the solver definition subroutine ---*/
  MPI::COMM_WORLD.Barrier();
#endif
  
  /*--- Surface grid deformation using design variables ---*/
  if (rank == MASTER_NODE) cout << endl << "------------------------- Surface grid deformation ----------------------" << endl;
  
  /*--- Definition and initialization of the surface deformation class ---*/
  surface_movement->CopyBoundary(geometry_container[MESH_0], config_container);
  
  /*--- Surface grid deformation ---*/
  if (rank == MASTER_NODE) cout << "Performing the deformation of the surface grid." << endl;
  surface_movement->SetAirfoil(geometry_container[MESH_0], config_container);
  
#ifndef NO_MPI
  /*--- MPI syncronization point ---*/
  MPI::COMM_WORLD.Barrier();
#endif
  
  /*--- Volumetric grid deformation ---*/
  if (rank == MASTER_NODE) cout << endl << "----------------------- Volumetric grid deformation ---------------------" << endl;
  
  /*--- Definition of the Class for grid movement ---*/
  
  if (rank == MASTER_NODE) cout << "Performing the deformation of the volumetric grid." << endl;
  grid_movement = new CVolumetricMovement(geometry_container[MESH_0]);
  grid_movement->SetVolume_Deformation(geometry_container[MESH_0], config_container, true);
  grid_movement->UpdateMultiGrid(geometry_container, config_container);
  
  /*--- Definition of the output class. The output class
   manages the writing of all restart, volume solution, surface solution,
   surface comma-separated value, and convergence history files (both in serial
   and in parallel). ---*/
  
  output = new COutput();
  
  /*--- Open the convergence history file ---*/
  
  if (rank == MASTER_NODE)
    output->SetHistory_Header(&ConvHist_file, config_container);
  
  /*--- Main external loop of the solver. Within this loop, each iteration ---*/
  
  if (rank == MASTER_NODE)
    cout << endl <<"------------------------------ Begin Solver -----------------------------" << endl;
  
#ifdef NO_MPI
  StartTime = double(clock())/double(CLOCKS_PER_SEC);
#else
  MPI::COMM_WORLD.Barrier();
  StartTime = MPI::Wtime();
#endif
  
  while (ExtIter < config_container->GetnExtIter()) {
    
    /*--- Set a timer for each iteration. Store the current iteration and
     update  the value of the CFL number (if there is CFL ramping specified)
     in the config class. ---*/
    
    config_container->SetExtIter(ExtIter);
    config_container->UpdateCFL(ExtIter);
    
    /*--- Perform a single iteration of the Euler/N-S/RANS solver. ---*/
    
    MeanFlowIteration(output, integration_container, geometry_container,
                      solver_container, numerics_container, config_container);
    
    /*--- Synchronization point after a single solver iteration. Compute the
     wall clock time required. ---*/
    
#ifdef NO_MPI
    StopTime = double(clock())/double(CLOCKS_PER_SEC);
#else
    MPI::COMM_WORLD.Barrier();
    StopTime = MPI::Wtime();
#endif
    
    UsedTime = (StopTime - StartTime);
    
    /*--- Update the convergence history file (serial and parallel computations). ---*/
    
    output->SetConvergence_History(&ConvHist_file, geometry_container, solver_container,
                                   config_container, integration_container, false, UsedTime);
    
    /*--- Check whether the current simulation has reached the specified
     convergence criteria, and set StopCalc to true, if so. ---*/
    
    StopCalc = integration_container[FLOW_SOL]->GetConvergence();
    
    /*--- Solution output. Determine whether a solution needs to be written
     after the current iteration, and if so, execute the output file writing
     routines. ---*/
    
    if ((ExtIter+1 == config_container->GetnExtIter()) ||
        ((ExtIter % config_container->GetWrt_Sol_Freq() == 0) && (ExtIter != 0)) ||
        (StopCalc)) {
      
      /*--- Execute the routine for writing restart, volume solution,
       surface solution, and surface comma-separated value files. ---*/
      
      output->SetResult_Files(solver_container, geometry_container, config_container, ExtIter);
      
    }
    
    /*--- If the convergence criteria has been met, terminate the simulation. ---*/
    
    if (StopCalc) break;
    ExtIter++;
    
  }
  
  /*--- Close the convergence history file. ---*/
  
  if (rank == MASTER_NODE) {
    ConvHist_file.close();
    cout << endl <<"History file, closed." << endl;
  }
  
  /*--- Solver class deallocation ---*/
  //    for (iMesh = 0; iMesh <= config_container->GetMGLevels(); iMesh++) {
  //      for (iSol = 0; iSol < MAX_SOLS; iSol++) {
  //        if (solver_container[iMesh][iSol] != NULL) {
  //          delete solver_container[iMesh][iSol];
  //        }
  //      }
  //      delete solver_container[iMesh];
  //    }
  //    delete solver_container;
  //  delete [] solver_container;
  //  if (rank == MASTER_NODE) cout <<"Solution container, deallocated." << endl;
  
  /*--- Geometry class deallocation ---*/
  //    for (iMesh = 0; iMesh <= config_container->GetMGLevels(); iMesh++) {
  //      delete geometry_container[iMesh];
  //    }
  //    delete geometry_container;
  //  delete [] geometry_container;
  //  cout <<"Geometry container, deallocated." << endl;
  
  /*--- Integration class deallocation ---*/
  //  cout <<"Integration container, deallocated." << endl;
  
#ifdef NO_MPI
  StopTime = double(clock())/double(CLOCKS_PER_SEC);
#else
  MPI::COMM_WORLD.Barrier();
  StopTime = MPI::Wtime();
#endif
  
  /*--- Compute/print the total time for performance benchmarking. ---*/
  
  UsedTime = StopTime-StartTime;
  if (rank == MASTER_NODE) {
    cout << "\nCompleted in " << fixed << UsedTime << " seconds on "<< size;
    if (size == 1) cout << " core." << endl; else cout << " cores." << endl;
  }
  
  /*--- Exit the solver cleanly ---*/
  
  if (rank == MASTER_NODE)
    cout << endl <<"------------------------- Exit Success (SU2_EDU) ------------------------" << endl << endl;
  
  /*--- Finalize MPI parallelization ---*/
#ifndef NO_MPI
  MPI::Finalize();
#endif
  
  return EXIT_SUCCESS;
  
}
