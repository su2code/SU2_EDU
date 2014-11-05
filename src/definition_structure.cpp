/*!
 * \file definition_structure.cpp
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

#include "../include/definition_structure.hpp"


unsigned short GetnZone(string val_mesh_filename, unsigned short val_format, CConfig *config) {
  string text_line, Marker_Tag;
  ifstream mesh_file;
  short nZone = 1; // Default value
  unsigned short iLine, nLine = 5;
  char cstr[200];
  string::size_type position;
  
  /*--- Search the mesh file for the 'NZONE' keyword. ---*/
  
  switch (val_format) {
    case SU2:
      
      /*--- Open grid file ---*/
      
      strcpy (cstr, val_mesh_filename.c_str());
      mesh_file.open(cstr, ios::in);
      if (mesh_file.fail()) {
        cout << "cstr=" << cstr << endl;
        cout << "There is no geometry file (GetnZone))!" << endl;
        
        exit(1);

      }
      
      /*--- Read the SU2 mesh file ---*/
      
      for (iLine = 0; iLine < nLine ; iLine++) {
        
        getline (mesh_file,text_line);
        
        /*--- Search for the "NZONE" keyword to see if there are multiple Zones ---*/
        position = text_line.find ("NZONE=",0);
        if (position != string::npos) {
          text_line.erase (0,6); nZone = atoi(text_line.c_str());
        }
      }
      
      break;
      
  }
  
  return (unsigned short) nZone;
}

unsigned short GetnDim(string val_mesh_filename, unsigned short val_format) {
  
  string text_line, Marker_Tag;
  ifstream mesh_file;
  short nDim = 3;
  unsigned short iLine, nLine = 5;
  char cstr[200];
  string::size_type position;
  
  /*--- Open grid file ---*/
  
  strcpy (cstr, val_mesh_filename.c_str());
  mesh_file.open(cstr, ios::in);

  switch (val_format) {
    case SU2:
      
      /*--- Read SU2 mesh file ---*/
      
      for (iLine = 0; iLine < nLine ; iLine++) {
        
        getline (mesh_file,text_line);
        
        /*--- Search for the "NDIM" keyword to see if there are multiple Zones ---*/
        
        position = text_line.find ("NDIME=",0);
        if (position != string::npos) {
          text_line.erase (0,6); nDim = atoi(text_line.c_str());
        }
      }
      break;
      
    case CGNS:
      
#ifdef HAVE_CGNS
      
      /*--- Local variables which are needed when calling the CGNS mid-level API. ---*/
      
      int fn, nbases, nzones, file_type;
      int cell_dim, phys_dim;
      char basename[CGNS_STRING_SIZE];
      
      /*--- Check whether the supplied file is truly a CGNS file. ---*/
      
      if ( cg_is_cgns(val_mesh_filename.c_str(),&file_type) != CG_OK ) {
        printf( "\n\n   !!! Error !!!\n" );
        printf( " %s is not a CGNS file.\n", val_mesh_filename.c_str());
        printf( " Now exiting...\n\n");
        exit(0);
      }
      
      /*--- Open the CGNS file for reading. The value of fn returned
       is the specific index number for this file and will be
       repeatedly used in the function calls. ---*/
      
      if ( cg_open(val_mesh_filename.c_str(),CG_MODE_READ,&fn) ) cg_error_exit();
      
      /*--- Get the number of databases. This is the highest node
       in the CGNS heirarchy. ---*/
      
      if ( cg_nbases(fn, &nbases) ) cg_error_exit();
      
      /*--- Check if there is more than one database. Throw an
       error if there is because this reader can currently
       only handle one database. ---*/
      
      if ( nbases > 1 ) {
        printf("\n\n   !!! Error !!!\n" );
        printf("CGNS reader currently incapable of handling more than 1 database.");
        printf("Now exiting...\n\n");
        exit(0);
      }
      
      /*--- Read the databases. Note that the indexing starts at 1. ---*/
      for ( int i = 1; i <= nbases; i++ ) {
        
        if ( cg_base_read(fn, i, basename, &cell_dim, &phys_dim) ) cg_error_exit();
        
        /*--- Get the number of zones for this base. ---*/
        
        if ( cg_nzones(fn, i, &nzones) ) cg_error_exit();
      
      }
      
      nDim = cell_dim;

#endif
      
      break;
      
  }
  
  mesh_file.close();

  return (unsigned short) nDim;
}

void Geometrical_Preprocessing(CGeometry **geometry, CConfig *config) {
  
  unsigned short iMGlevel;
  unsigned long iPoint;
  int rank = MASTER_NODE;
  
  /*--- Compute elements surrounding points, points surrounding points ---*/
  
  if (rank == MASTER_NODE) cout << "Setting point connectivity." << endl;
  geometry[MESH_0]->SetPoint_Connectivity();
  
  /*--- Renumbering points using Reverse Cuthill McKee ordering ---*/
  
  if (rank == MASTER_NODE) cout << "Renumbering points (Reverse Cuthill McKee Ordering)." << endl;
  geometry[MESH_0]->SetRCM_Ordering(config);
  
  /*--- recompute elements surrounding points, points surrounding points ---*/
  
  if (rank == MASTER_NODE) cout << "Recomputing point connectivity." << endl;
  geometry[MESH_0]->SetPoint_Connectivity();
  
  /*--- Compute elements surrounding elements ---*/
  
  if (rank == MASTER_NODE) cout << "Setting element connectivity." << endl;
  geometry[MESH_0]->SetElement_Connectivity();
  
  /*--- Check the orientation before computing geometrical quantities ---*/
  
  if (rank == MASTER_NODE) cout << "Checking the numerical grid orientation." << endl;
  geometry[MESH_0]->SetBoundVolume();
  geometry[MESH_0]->Check_Orientation(config);
  
  /*--- Create the edge structure ---*/
  
  if (rank == MASTER_NODE) cout << "Identifying edges and vertices." << endl;
  geometry[MESH_0]->SetEdges();
  geometry[MESH_0]->SetVertex(config);
  
  /*--- Compute cell center of gravity ---*/
  
  if (rank == MASTER_NODE) cout << "Computing centers of gravity." << endl;
  geometry[MESH_0]->SetCG();
  
  /*--- Create the control volume structures ---*/
  
  if (rank == MASTER_NODE) cout << "Setting the control volume structure." << endl;
  geometry[MESH_0]->SetControlVolume(config, ALLOCATE);
  geometry[MESH_0]->SetBoundControlVolume(config, ALLOCATE);
  
  /*--- Identify closest normal neighbor ---*/
  
  if (rank == MASTER_NODE) cout << "Searching for the closest normal neighbors to the surfaces." << endl;
  geometry[MESH_0]->FindNormal_Neighbor(config);
  
  if ((config->GetMGLevels() != 0) && (rank == MASTER_NODE))
    cout << "Setting the multigrid structure." <<endl;
  
  /*--- Loop over all the new grid ---*/
  
  for (iMGlevel = 1; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
    
    /*--- Create main agglomeration structure ---*/
    
    geometry[iMGlevel] = new CMultiGridGeometry(geometry, config, iMGlevel);
    
    /*--- Compute points surrounding points. ---*/
    
    geometry[iMGlevel]->SetPoint_Connectivity(geometry[iMGlevel-1]);
    
    /*--- Create the edge structure ---*/
    
    geometry[iMGlevel]->SetEdges();
    geometry[iMGlevel]->SetVertex(geometry[iMGlevel-1], config);
    
    /*--- Create the control volume structures ---*/
    
    geometry[iMGlevel]->SetControlVolume(config,geometry[iMGlevel-1], ALLOCATE);
    geometry[iMGlevel]->SetBoundControlVolume(config,geometry[iMGlevel-1], ALLOCATE);
    geometry[iMGlevel]->SetCoord(geometry[iMGlevel-1]);
    
    /*--- Find closest neighbor to a surface point ---*/
    
    geometry[iMGlevel]->FindNormal_Neighbor(config);
    
  }
  
  /*--- For unsteady simulations, initialize the grid volumes
   and coordinates for previous solutions. Loop over all zones/grids ---*/
  
  if (config->GetUnsteady_Simulation() && config->GetGrid_Movement()) {
    for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
      for (iPoint = 0; iPoint < geometry[iMGlevel]->GetnPoint(); iPoint++) {
        
        /*--- Update cell volume ---*/
        
        geometry[iMGlevel]->node[iPoint]->SetVolume_n();
        geometry[iMGlevel]->node[iPoint]->SetVolume_nM1();
        
        /*--- Update point coordinates ---*/
        geometry[iMGlevel]->node[iPoint]->SetCoord_n();
        geometry[iMGlevel]->node[iPoint]->SetCoord_n1();
        
      }
    }
  }
  
}

void Solver_Preprocessing(CSolver ***solver_container, CGeometry **geometry,
                          CConfig *config) {
  
  unsigned short iMGlevel;
  bool euler, ns, turbulent,
  adj_euler, adj_ns, adj_turb,
  lin_euler, lin_ns, lin_turb,
  tne2_euler, tne2_ns,
  adj_tne2_euler, adj_tne2_ns,
  poisson, wave, fea, heat,
  spalart_allmaras, menter_sst, transition,
  template_solver;
  
  
  /*--- Initialize some useful booleans ---*/
  euler            = false;  ns              = false;  turbulent = false;
  adj_euler        = false;	 adj_ns          = false;	 adj_turb  = false;
  lin_euler        = false;	 lin_ns          = false;  lin_turb  = false;
  tne2_euler       = false;  tne2_ns         = false;
  adj_tne2_euler   = false;  adj_tne2_ns     = false;
  spalart_allmaras = false;  menter_sst      = false;
  poisson          = false;
  wave             = false;
  fea              = false;
  heat             = false;
  transition       = false;
  template_solver  = false;
  
  /*--- Assign booleans ---*/
  switch (config->GetKind_Solver()) {
    case TEMPLATE_SOLVER: template_solver = true; break;
    case EULER : euler = true; break;
    case NAVIER_STOKES: ns = true; break;
    case RANS : ns = true; turbulent = true; if (config->GetKind_Trans_Model() == LM) transition = true; break;
  }
  /*--- Assign turbulence model booleans --- */
  if (turbulent)
    switch (config->GetKind_Turb_Model()){
      case SA: spalart_allmaras = true; break;
      case SST: menter_sst = true; break;
      default: cout << "Specified turbulence model unavailable or none selected" << endl; exit(1); break;
    }
  
  /*--- Definition of the Class for the solution: solver_container[DOMAIN][MESH_LEVEL][EQUATION]. Note that euler, ns
   and potential are incompatible, they use the same position in sol container ---*/
  for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
    
    /*--- Allocate solution for direct problem, and run the preprocessing and postprocessing ---*/
    if (euler) {
      solver_container[iMGlevel][FLOW_SOL] = new CEulerSolver(geometry[iMGlevel], config, iMGlevel);
    }
    if (ns) {
      solver_container[iMGlevel][FLOW_SOL] = new CNSSolver(geometry[iMGlevel], config, iMGlevel);
    }
    if (turbulent) {
      if (spalart_allmaras) {
        solver_container[iMGlevel][TURB_SOL] = new CTurbSASolver(geometry[iMGlevel], config, iMGlevel);
        solver_container[iMGlevel][FLOW_SOL]->Preprocessing(geometry[iMGlevel], solver_container[iMGlevel], config, iMGlevel, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
        solver_container[iMGlevel][TURB_SOL]->Postprocessing(geometry[iMGlevel], solver_container[iMGlevel], config, iMGlevel);
      }
      else if (menter_sst) {
        solver_container[iMGlevel][TURB_SOL] = new CTurbSSTSolver(geometry[iMGlevel], config, iMGlevel);
        solver_container[iMGlevel][FLOW_SOL]->Preprocessing(geometry[iMGlevel], solver_container[iMGlevel], config, iMGlevel, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
        solver_container[iMGlevel][TURB_SOL]->Postprocessing(geometry[iMGlevel], solver_container[iMGlevel], config, iMGlevel);
      }
    }
    
  }
  
}

void Integration_Preprocessing(CIntegration **integration_container,
                               CGeometry **geometry, CConfig *config) {
  
  bool
  euler, adj_euler, lin_euler,
  ns, adj_ns, lin_ns,
  turbulent, adj_turb, lin_turb,
  spalart_allmaras, menter_sst,
  tne2_euler, adj_tne2_euler,
  tne2_ns, adj_tne2_ns,
  poisson, wave, fea, heat, template_solver, transition;
  
  /*--- Initialize some useful booleans ---*/
  euler            = false; adj_euler        = false; lin_euler         = false;
  ns               = false; adj_ns           = false; lin_ns            = false;
  turbulent        = false; adj_turb         = false; lin_turb          = false;
  spalart_allmaras = false; menter_sst       = false;
  tne2_euler       = false; adj_tne2_euler   = false;
  tne2_ns          = false; adj_tne2_ns      = false;
  poisson          = false;
  wave             = false;
  heat             = false;
  fea              = false;
  transition       = false;
  template_solver  = false;
  
  /*--- Assign booleans ---*/
  switch (config->GetKind_Solver()) {
    case TEMPLATE_SOLVER: template_solver = true; break;
    case EULER : euler = true; break;
    case NAVIER_STOKES: ns = true; break;
    case RANS : ns = true; turbulent = true; if (config->GetKind_Trans_Model() == LM) transition = true; break;
  }
  
  /*--- Assign turbulence model booleans --- */
  if (turbulent) {
    switch (config->GetKind_Turb_Model()) {
      case SA: spalart_allmaras = true; break;
      case SST: menter_sst = true; break;
      default: cout << "Specified turbulence model unavailable or none selected" << endl; exit(1); break;
    }
  }
  
  /*--- Allocate solution for a template problem ---*/
  if (template_solver) integration_container[TEMPLATE_SOL] = new CSingleGridIntegration(config);
  
  /*--- Allocate solution for direct problem ---*/
  if (euler) integration_container[FLOW_SOL] = new CMultiGridIntegration(config);
  if (ns) integration_container[FLOW_SOL] = new CMultiGridIntegration(config);
  if (turbulent) integration_container[TURB_SOL] = new CSingleGridIntegration(config);
  
}

void Numerics_Preprocessing(CNumerics ****numerics_container,
                            CSolver ***solver_container, CGeometry **geometry,
                            CConfig *config) {
  
  unsigned short iMGlevel, iSol, nDim,
  
  nVar_Flow             = 0,
  nVar_Turb             = 0;

  
  double *constants = NULL;
  
  bool
  euler, adj_euler, lin_euler,
  ns, adj_ns, lin_ns,
  turbulent, adj_turb, lin_turb,
  tne2_euler, adj_tne2_euler,
  tne2_ns, adj_tne2_ns,
  spalart_allmaras, menter_sst,
  poisson,
  wave,
  fea,
  heat,
  transition,
  template_solver;
  
  /*--- Initialize some useful booleans ---*/
  euler            = false;   ns               = false;   turbulent        = false;
  poisson          = false;
  adj_euler        = false;	  adj_ns           = false;	 adj_turb         = false;
  wave             = false;   heat             = false;   fea              = false;   spalart_allmaras = false;
  tne2_euler       = false;   tne2_ns          = false;
  adj_tne2_euler   = false;	  adj_tne2_ns      = false;
  lin_euler        = false;   lin_ns           = false;   lin_turb         = false;	 menter_sst       = false;
  transition       = false;
  template_solver  = false;
  
  /*--- Assign booleans ---*/
  switch (config->GetKind_Solver()) {
    case EULER : euler = true; break;
    case NAVIER_STOKES: ns = true; break;
    case RANS : ns = true; turbulent = true; if (config->GetKind_Trans_Model() == LM) transition = true; break;
  }
  
  /*--- Assign turbulence model booleans --- */
  if (turbulent)
    switch (config->GetKind_Turb_Model()){
      case SA: spalart_allmaras = true; break;
      case SST: menter_sst = true; constants = solver_container[MESH_0][TURB_SOL]->GetConstants(); break;
      default: cout << "Specified turbulence model unavailable or none selected" << endl; exit(1); break;
    }
  
  /*--- Number of variables for the template ---*/
  if (template_solver) nVar_Flow = solver_container[MESH_0][FLOW_SOL]->GetnVar();
  
  /*--- Number of variables for direct problem ---*/
  if (euler)				nVar_Flow = solver_container[MESH_0][FLOW_SOL]->GetnVar();
  if (ns)	          nVar_Flow = solver_container[MESH_0][FLOW_SOL]->GetnVar();
  if (turbulent)		nVar_Turb = solver_container[MESH_0][TURB_SOL]->GetnVar();
  
  /*--- Number of dimensions ---*/
  nDim = geometry[MESH_0]->GetnDim();
  
  /*--- Definition of the Class for the numerical method: numerics_container[MESH_LEVEL][EQUATION][EQ_TERM] ---*/
  for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
    numerics_container[iMGlevel] = new CNumerics** [MAX_SOLS];
    for (iSol = 0; iSol < MAX_SOLS; iSol++)
      numerics_container[iMGlevel][iSol] = new CNumerics* [MAX_TERMS];
  }
  
  /*--- Solver definition for the Potential, Euler, Navier-Stokes problems ---*/
  if ((euler) || (ns)) {
    
    /*--- Definition of the convective scheme for each equation and mesh level ---*/
    switch (config->GetKind_ConvNumScheme_Flow()) {
      case NO_CONVECTIVE :
        cout << "No convective scheme." << endl; exit(1);
        break;
        
      case SPACE_CENTERED :
          switch (config->GetKind_Centered_Flow()) {
            case NO_CENTERED : cout << "No centered scheme." << endl; break;
            case LAX : numerics_container[MESH_0][FLOW_SOL][CONV_TERM] = new CCentLax_Flow(nDim,nVar_Flow, config); break;
            case JST : numerics_container[MESH_0][FLOW_SOL][CONV_TERM] = new CCentJST_Flow(nDim,nVar_Flow, config); break;
            case JST_KE : numerics_container[MESH_0][FLOW_SOL][CONV_TERM] = new CCentJST_KE_Flow(nDim,nVar_Flow, config); break;
          default : cout << "Centered scheme not implemented." << endl; exit(1); break;
          }
          
          if (!config->GetLowFidelitySim()) {
            for (iMGlevel = 1; iMGlevel <= config->GetMGLevels(); iMGlevel++)
              numerics_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CCentLax_Flow(nDim, nVar_Flow, config);
          }
          else {
            numerics_container[MESH_1][FLOW_SOL][CONV_TERM] = new CCentJST_Flow(nDim, nVar_Flow, config);
            for (iMGlevel = 2; iMGlevel <= config->GetMGLevels(); iMGlevel++)
              numerics_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CCentLax_Flow(nDim, nVar_Flow, config);
          }
          
          /*--- Definition of the boundary condition method ---*/
          for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
            numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwRoe_Flow(nDim, nVar_Flow, config);
          
        break;
      case SPACE_UPWIND :
          switch (config->GetKind_Upwind_Flow()) {
            case NO_UPWIND : cout << "No upwind scheme." << endl; break;
            case ROE:
              for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
                numerics_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwRoe_Flow(nDim, nVar_Flow, config);
                numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwRoe_Flow(nDim, nVar_Flow, config);
              }
              break;
              
            case AUSM:
              for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
                numerics_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwAUSM_Flow(nDim, nVar_Flow, config);
                numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwAUSM_Flow(nDim, nVar_Flow, config);
              }
              break;
              
            case TURKEL:
              for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
                numerics_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwTurkel_Flow(nDim, nVar_Flow, config);
                numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwTurkel_Flow(nDim, nVar_Flow, config);
              }
              break;
              
            case HLLC:
              for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
                numerics_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwHLLC_Flow(nDim, nVar_Flow, config);
                numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwHLLC_Flow(nDim, nVar_Flow, config);
              }
              break;
              
            case MSW:
              for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
                numerics_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwMSW_Flow(nDim, nVar_Flow, config);
                numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwMSW_Flow(nDim, nVar_Flow, config);
              }
              break;
              
            case CUSP:
              for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
                numerics_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwCUSP_Flow(nDim, nVar_Flow, config);
                numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwCUSP_Flow(nDim, nVar_Flow, config);
              }
              break;
              
            default : cout << "Upwind scheme not implemented." << endl; exit(1); break;
          }
        
        break;
        
      default :
        cout << "Convective scheme not implemented (euler and ns)." << endl; exit(1);
        break;
    }
    
    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    switch (config->GetKind_ViscNumScheme_Flow()) {
      case NONE :
        break;
      case AVG_GRAD :
          for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
            numerics_container[iMGlevel][FLOW_SOL][VISC_TERM] = new CAvgGrad_Flow(nDim, nVar_Flow, config);
            numerics_container[iMGlevel][FLOW_SOL][VISC_BOUND_TERM] = new CAvgGrad_Flow(nDim, nVar_Flow, config);
          }
        break;
      case AVG_GRAD_CORRECTED :
          numerics_container[MESH_0][FLOW_SOL][VISC_TERM] = new CAvgGradCorrected_Flow(nDim, nVar_Flow, config);
          for (iMGlevel = 1; iMGlevel <= config->GetMGLevels(); iMGlevel++)
            numerics_container[iMGlevel][FLOW_SOL][VISC_TERM] = new CAvgGrad_Flow(nDim, nVar_Flow, config);
          
          /*--- Definition of the boundary condition method ---*/
          for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
            numerics_container[iMGlevel][FLOW_SOL][VISC_BOUND_TERM] = new CAvgGrad_Flow(nDim, nVar_Flow, config);
        break;
      case GALERKIN :
        cout << "Galerkin viscous scheme not implemented." << endl; exit(1); exit(1);
        break;
      default :
        cout << "Numerical viscous scheme not recognized." << endl; exit(1); exit(1);
        break;
    }
    
    /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
    switch (config->GetKind_SourNumScheme_Flow()) {
      case NONE :
        break;
      case PIECEWISE_CONSTANT :
        
        for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
          numerics_container[iMGlevel][FLOW_SOL][SOURCE_FIRST_TERM] = new CSourceNothing(nDim, nVar_Flow, config);
          numerics_container[iMGlevel][FLOW_SOL][SOURCE_SECOND_TERM] = new CSourceNothing(nDim, nVar_Flow, config);
        }
        
        break;
      default :
        cout << "Source term not implemented." << endl; exit(1);
        break;
    }
    
  }
  
  /*--- Solver definition for the turbulent model problem ---*/
  if (turbulent) {
    
    /*--- Definition of the convective scheme for each equation and mesh level ---*/
    switch (config->GetKind_ConvNumScheme_Turb()) {
      case NONE :
        break;
      case SPACE_UPWIND :
        for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++){
          if (spalart_allmaras) numerics_container[iMGlevel][TURB_SOL][CONV_TERM] = new CUpwSca_TurbSA(nDim, nVar_Turb, config);
          else if (menter_sst) numerics_container[iMGlevel][TURB_SOL][CONV_TERM] = new CUpwSca_TurbSST(nDim, nVar_Turb, config);
        }
        break;
      default :
        cout << "Convective scheme not implemented (turbulent)." << endl; exit(1);
        break;
    }
    
    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    switch (config->GetKind_ViscNumScheme_Turb()) {
      case NONE :
        break;
      case AVG_GRAD :
        for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++){
          if (spalart_allmaras) numerics_container[iMGlevel][TURB_SOL][VISC_TERM] = new CAvgGrad_TurbSA(nDim, nVar_Turb, config);
          else if (menter_sst) numerics_container[iMGlevel][TURB_SOL][VISC_TERM] = new CAvgGrad_TurbSST(nDim, nVar_Turb, constants, config);
        }
        break;
      case AVG_GRAD_CORRECTED :
        for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++){
          if (spalart_allmaras) numerics_container[iMGlevel][TURB_SOL][VISC_TERM] = new CAvgGradCorrected_TurbSA(nDim, nVar_Turb, config);
          else if (menter_sst) numerics_container[iMGlevel][TURB_SOL][VISC_TERM] = new CAvgGradCorrected_TurbSST(nDim, nVar_Turb, constants, config);
        }
        break;
      case GALERKIN :
        cout << "Viscous scheme not implemented." << endl;
        exit(1); break;
      default :
        cout << "Viscous scheme not implemented." << endl; exit(1);
        break;
    }
    
    /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
    switch (config->GetKind_SourNumScheme_Turb()) {
      case NONE :
        break;
      case PIECEWISE_CONSTANT :
        for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
          if (spalart_allmaras) numerics_container[iMGlevel][TURB_SOL][SOURCE_FIRST_TERM] = new CSourcePieceWise_TurbSA(nDim, nVar_Turb, config);
          else if (menter_sst) numerics_container[iMGlevel][TURB_SOL][SOURCE_FIRST_TERM] = new CSourcePieceWise_TurbSST(nDim, nVar_Turb, constants, config);
          numerics_container[iMGlevel][TURB_SOL][SOURCE_SECOND_TERM] = new CSourceNothing(nDim, nVar_Turb, config);
        }
        break;
      default :
        cout << "Source term not implemented." << endl; exit(1);
        break;
    }
    
    /*--- Definition of the boundary condition method ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++){
      if (spalart_allmaras) {
        numerics_container[iMGlevel][TURB_SOL][CONV_BOUND_TERM] = new CUpwSca_TurbSA(nDim, nVar_Turb, config);
        numerics_container[iMGlevel][TURB_SOL][VISC_BOUND_TERM] = new CAvgGrad_TurbSA(nDim, nVar_Turb, config);
      }
      else if (menter_sst) {
        numerics_container[iMGlevel][TURB_SOL][CONV_BOUND_TERM] = new CUpwSca_TurbSST(nDim, nVar_Turb, config);
        numerics_container[iMGlevel][TURB_SOL][VISC_BOUND_TERM] = new CAvgGrad_TurbSST(nDim, nVar_Turb, constants, config);
      }
    }
  }
  
}
