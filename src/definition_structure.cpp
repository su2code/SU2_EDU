/*!
 * \file definition_structure.cpp
 * \brief Main subroutines used by SU2_CFD.
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

#include "../include/definition_structure.hpp"

unsigned short GetnDim(string val_mesh_filename, unsigned short val_format) {
  
  string text_line, Marker_Tag;
  ifstream mesh_file;
  short nDim = 3;
  bool isFound = false;
  char cstr[200];
  string::size_type position;
  
  switch (val_format) {
    case SU2:
      
      /*--- Open grid file ---*/
      strcpy (cstr, val_mesh_filename.c_str());
      mesh_file.open(cstr, ios::in);
      
      /*--- Read SU2 mesh file ---*/
      while (getline (mesh_file,text_line)) {
        /*--- Search for the "NDIM" keyword to see if there are multiple Zones ---*/
        position = text_line.find ("NDIME=",0);
        if (position != string::npos) {
          text_line.erase (0,6); nDim = atoi(text_line.c_str()); isFound = true;
        }
      }
      break;
      
    case CGNS:
      nDim = 3;
      break;
      
    case NETCDF_ASCII:
      nDim = 3;
      break;
  }
  return (unsigned short) nDim;
}

void Geometrical_Preprocessing(CGeometry **geometry, CConfig *config) {
  
  unsigned short iMGlevel;
  
  /*--- Compute elements surrounding points, points surrounding points,
   and elements surrounding elements ---*/
  
  cout << "Setting point connectivity." << endl;
  geometry[MESH_0]->SetEsuP();
  geometry[MESH_0]->SetPsuP();
  
  cout << "Renumbering points." << endl;
  //geometry[MESH_0]->SetRCM(config);
  
  cout << "Recomputing point connectivity." << endl;
  geometry[MESH_0]->SetEsuP();
  geometry[MESH_0]->SetPsuP();
  
  cout << "Setting element connectivity." << endl;
  geometry[MESH_0]->SetEsuE();
  
  /*--- Check the orientation before computing geometrical quantities ---*/
  
  cout << "Checking the numerical grid orientation." << endl;
  geometry[MESH_0]->SetBoundVolume();
  geometry[MESH_0]->Check_Orientation(config);
  
  /*--- Create the edge structure ---*/
  
  cout << "Identifying edges and vertices." << endl;
  geometry[MESH_0]->SetEdges();
  geometry[MESH_0]->SetVertex(config);

  /*--- Color the edges for shared memory parallelism ---*/
  
  geometry[MESH_0]->Color_Edges(config);
  
  cout << "Recomputing point connectivity." << endl;
  geometry[MESH_0]->SetEsuP();
  geometry[MESH_0]->SetPsuP();
  
  cout << "Setting element connectivity." << endl;
  geometry[MESH_0]->SetEsuE();
  
  /*--- Check the orientation before computing geometrical quantities ---*/
  
  cout << "Checking the numerical grid orientation." << endl;
  geometry[MESH_0]->SetBoundVolume();
  geometry[MESH_0]->Check_Orientation(config);
  
  /*--- Create the edge structure ---*/
  
  cout << "Identifying edges and vertices." << endl;
  geometry[MESH_0]->SetEdges();
  geometry[MESH_0]->SetVertex(config);
  
  /*--- Compute cell center of gravity ---*/
  
  cout << "Computing centers of gravity." << endl;
  geometry[MESH_0]->SetCG();
  
  /*--- Create the control volume structures ---*/
  
  cout << "Setting the control volume structure." << endl;
  geometry[MESH_0]->SetControlVolume(config, ALLOCATE);
  geometry[MESH_0]->SetBoundControlVolume(config, ALLOCATE);
  
  /*--- Identify closest normal neighbor ---*/
  
  cout << "Searching for the closest normal neighbors to the surfaces." << endl;
  geometry[MESH_0]->FindNormal_Neighbor(config);
  
  /*--- Computation of wall distances for turbulence modeling ---*/
  
  if (config->GetKind_Solver() == RANS)
    geometry[MESH_0]->ComputeWall_Distance(config);
  
  /*--- Computation of positive surface area in the y-direction (2-D) or the
   z-direction (3-D) which is used to calculate force coefficients (non-dim).
   This value is used if the reference area is set to zero in the config. ---*/
  
  geometry[MESH_0]->ComputeReference_Area(config);
  
  
  if (config->GetMGLevels() != 0)
    cout << "Setting the multigrid structure." <<endl;
  
  /*--- Loop over all the new grid ---*/
  
  for (iMGlevel = 1; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
    
    /*--- Create main agglomeration structure ---*/
    
    geometry[iMGlevel] = new CMultiGridGeometry(geometry, config, iMGlevel);
    
    /*--- Compute points surrounding points. ---*/
    
    geometry[iMGlevel]->SetPsuP(geometry[iMGlevel-1]);
    
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
  
}

void Solver_Preprocessing(CSolver ***solver_container, CGeometry **geometry, CConfig *config) {
  
  unsigned short iMGlevel;
  bool euler, ns, turbulent, spalart_allmaras, menter_sst;
  
  /*--- Initialize some useful booleans ---*/
  euler = false;  ns = false;  turbulent = false;
  menter_sst = false;  spalart_allmaras = false;
  
  /*--- Assign booleans ---*/
  switch (config->GetKind_Solver()) {
    case EULER : euler = true; break;
    case NAVIER_STOKES: ns = true; break;
    case RANS : ns = true; turbulent = true; break;
  }
  /*--- Assign turbulence model booleans --- */
  if (turbulent) {
    switch (config->GetKind_Turb_Model()) {
      case SA: spalart_allmaras = true; break;
      case SST: menter_sst = true; break;
        
      default: cout << "Specified turbulence model unavailable or none selected" << endl; exit(1); break;
    }
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
        solver_container[iMGlevel][FLOW_SOL]->Preprocessing(geometry[iMGlevel], solver_container[iMGlevel], config, iMGlevel, NO_RK_ITER, RUNTIME_FLOW_SYS);
        solver_container[iMGlevel][TURB_SOL]->Postprocessing(geometry[iMGlevel], solver_container[iMGlevel], config, iMGlevel);
      }
      else if (menter_sst) {
        solver_container[iMGlevel][TURB_SOL] = new CTurbSSTSolver(geometry[iMGlevel], config, iMGlevel);
        solver_container[iMGlevel][FLOW_SOL]->Preprocessing(geometry[iMGlevel], solver_container[iMGlevel], config, iMGlevel, NO_RK_ITER, RUNTIME_FLOW_SYS);
        solver_container[iMGlevel][TURB_SOL]->Postprocessing(geometry[iMGlevel], solver_container[iMGlevel], config, iMGlevel);
      }
    }
    
  }
  
}

void Integration_Preprocessing(CIntegration **integration_container,
                               CGeometry **geometry, CConfig *config) {
  
  bool
  euler, ns, turbulent, spalart_allmaras, menter_sst;
  
  /*--- Initialize some useful booleans ---*/
  euler            = false;   ns               = false;  turbulent        = false;
  spalart_allmaras = false;   menter_sst       = false;
  
  /*--- Assign booleans ---*/
  switch (config->GetKind_Solver()) {
    case EULER : euler = true; break;
    case NAVIER_STOKES: ns = true; break;
    case RANS : ns = true; turbulent = true; break;
  }
  
  /*--- Assign turbulence model booleans --- */
  if (turbulent) {
    switch (config->GetKind_Turb_Model()) {
      case SA: spalart_allmaras = true; break;
      case SST: menter_sst = true; break;
      default: cout << "Specified turbulence model unavailable or none selected" << endl; exit(1); break;
    }
  }
  
  /*--- Allocate solution for direct problem ---*/
  if (euler) integration_container[FLOW_SOL] = new CMultiGridIntegration(config);
  if (ns) integration_container[FLOW_SOL] = new CMultiGridIntegration(config);
  if (turbulent) integration_container[TURB_SOL] = new CSingleGridIntegration(config);
  
}

void Numerics_Preprocessing(CNumerics ****numerics_container,
                            CSolver ***solver_container, CGeometry **geometry,
                            CConfig *config) {
  
  unsigned short iMGlevel, iSol, nDim,
  
  nVar_Flow         = 0,
  nVar_Turb         = 0;
  
  double *constants = NULL;
  
  bool
  euler, ns, turbulent,
  spalart_allmaras, menter_sst;
  
  /*--- Initialize some useful booleans ---*/
  euler            = false;   ns               = false;   turbulent        = false;
  spalart_allmaras = false; menter_sst       = false;
  
  /*--- Assign booleans ---*/
  switch (config->GetKind_Solver()) {
    case EULER : euler = true; break;
    case NAVIER_STOKES: ns = true; break;
    case RANS : ns = true; turbulent = true; break;
  }
  
  /*--- Assign turbulence model booleans --- */
  if (turbulent) {
    switch (config->GetKind_Turb_Model()){
      case SA: spalart_allmaras = true; break;
      case SST: menter_sst = true; constants = solver_container[MESH_0][TURB_SOL]->GetConstants(); break;
      default: cout << "Specified turbulence model unavailable or none selected" << endl; exit(1); break;
    }
  }
  
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
        /*--- Compressible flow ---*/
        switch (config->GetKind_Centered_Flow()) {
          case NO_CENTERED : cout << "No centered scheme." << endl; break;
          case LAX : numerics_container[MESH_0][FLOW_SOL][CONV_TERM] = new CCentLax_Flow(nDim,nVar_Flow, config); break;
          case JST : numerics_container[MESH_0][FLOW_SOL][CONV_TERM] = new CCentJST_Flow(nDim,nVar_Flow, config); break;
          default : cout << "Centered scheme not implemented." << endl; exit(1); break;
        }
        
        /*--- Definition of the numerics on the coarse levels ---*/
        
        for (iMGlevel = 1; iMGlevel <= config->GetMGLevels(); iMGlevel++)
          numerics_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CCentLax_Flow(nDim, nVar_Flow, config);
        
        /*--- Definition of the boundary condition method ---*/
        for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
          numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwRoe_Flow(nDim, nVar_Flow, config);
        
        break;
      case SPACE_UPWIND :
        /*--- Compressible flow ---*/
        switch (config->GetKind_Upwind_Flow()) {
          case NO_UPWIND : cout << "No upwind scheme." << endl; break;
          case ROE_1ST : case ROE_2ND :
            for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
              numerics_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwRoe_Flow(nDim, nVar_Flow, config);
              numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwRoe_Flow(nDim, nVar_Flow, config);
            }
            break;
            
          case AUSM_1ST : case AUSM_2ND :
            for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
              numerics_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwAUSM_Flow(nDim, nVar_Flow, config);
              numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwAUSM_Flow(nDim, nVar_Flow, config);
            }
            break;
            
          case ROE_TURKEL_1ST : case ROE_TURKEL_2ND :
            for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
              numerics_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwRoeTurkel_Flow(nDim, nVar_Flow, config);
              numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwRoeTurkel_Flow(nDim, nVar_Flow, config);
            }
            break;
            
          case HLLC_1ST : case HLLC_2ND :
            for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
              numerics_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwHLLC_Flow(nDim, nVar_Flow, config);
              numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwHLLC_Flow(nDim, nVar_Flow, config);
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
        /*--- Compressible flow ---*/
        for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
          numerics_container[iMGlevel][FLOW_SOL][VISC_TERM] = new CAvgGrad_Flow(nDim, nVar_Flow, config);
          numerics_container[iMGlevel][FLOW_SOL][VISC_BOUND_TERM] = new CAvgGrad_Flow(nDim, nVar_Flow, config);
        }
        break;
      case AVG_GRAD_CORRECTED :
        /*--- Compressible flow ---*/
        numerics_container[MESH_0][FLOW_SOL][VISC_TERM] = new CAvgGradCorrected_Flow(nDim, nVar_Flow, config);
        for (iMGlevel = 1; iMGlevel <= config->GetMGLevels(); iMGlevel++)
          numerics_container[iMGlevel][FLOW_SOL][VISC_TERM] = new CAvgGrad_Flow(nDim, nVar_Flow, config);
        
        /*--- Definition of the boundary condition method ---*/
        for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
          numerics_container[iMGlevel][FLOW_SOL][VISC_BOUND_TERM] = new CAvgGrad_Flow(nDim, nVar_Flow, config);
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
