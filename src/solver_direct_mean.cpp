/*!
 * \file solution_direct_mean.cpp
 * \brief Main subrotuines for solving direct problems (Euler, Navier-Stokes, etc.).
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

#include "../include/solver_structure.hpp"

CEulerSolver::CEulerSolver(void) : CSolver() {
  
  /*--- Basic array initialization ---*/
  
  CDrag_Inv = NULL; CLift_Inv = NULL; CSideForce_Inv = NULL;  CEff_Inv = NULL;
  CMx_Inv = NULL; CMy_Inv = NULL; CMz_Inv = NULL;
  CFx_Inv = NULL; CFy_Inv = NULL; CFz_Inv = NULL;
  
  CPressure = NULL; CPressureTarget = NULL; HeatFlux = NULL; HeatFluxTarget = NULL; YPlus = NULL;
  ForceInviscid = NULL; MomentInviscid = NULL;
  
  /*--- Surface based array initialization ---*/
  
  Surface_CLift_Inv = NULL; Surface_CDrag_Inv = NULL;
  Surface_CMx_Inv = NULL; Surface_CMy_Inv = NULL; Surface_CMz_Inv = NULL;
  
  Surface_CLift = NULL; Surface_CDrag = NULL;
  Surface_CMx = NULL; Surface_CMy = NULL; Surface_CMz = NULL;
  
  /*--- Rotorcraft simulation array initialization ---*/
  
  CMerit_Inv = NULL;  CT_Inv = NULL;  CQ_Inv = NULL;
  
  /*--- Supersonic simulation array initialization ---*/
  
  CEquivArea_Inv = NULL;
  CNearFieldOF_Inv = NULL;
  
  /*--- Engine simulation array initialization ---*/
  
  FanFace_MassFlow = NULL;  FanFace_Pressure = NULL;
  FanFace_Mach = NULL;  FanFace_Area = NULL;
  Exhaust_MassFlow = NULL;  Exhaust_Area = NULL;
  
  /*--- Numerical methods array initialization ---*/
  
  iPoint_UndLapl = NULL;
  jPoint_UndLapl = NULL;
  LowMach_Precontioner = NULL;
  Primitive = NULL; Primitive_i = NULL; Primitive_j = NULL;
  CharacPrimVar = NULL;
  
  /*--- Fixed CL mode initialization (cauchy criteria) ---*/
  
  Cauchy_Value = 0;
  Cauchy_Func = 0;
  Old_Func = 0;
  New_Func = 0;
  Cauchy_Counter = 0;
  Cauchy_Serie = NULL;
  
}

CEulerSolver::CEulerSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh) : CSolver() {
  
  unsigned long iPoint, index, counter_local = 0, counter_global = 0, iVertex;
  unsigned short iVar, iDim, iMarker, nLineLets;
  double Density, Velocity2, Pressure, Temperature, dull_val;
  
  unsigned short nZone = geometry->GetnZone();
  bool restart = (config->GetRestart() || config->GetRestart_Flow());
  double Gas_Constant = config->GetGas_ConstantND();
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) || (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  bool roe_turkel = (config->GetKind_Upwind_Flow() == TURKEL);
  bool adjoint = config->GetAdjoint();
  
  int rank = MASTER_NODE;
  
  /*--- Array initialization ---*/
  
  CDrag_Inv = NULL; CLift_Inv = NULL; CSideForce_Inv = NULL;  CEff_Inv = NULL;
  CMx_Inv = NULL; CMy_Inv = NULL; CMz_Inv = NULL;
  CFx_Inv = NULL; CFy_Inv = NULL; CFz_Inv = NULL;
  Surface_CLift_Inv = NULL; Surface_CDrag_Inv = NULL;
  Surface_CMx_Inv = NULL; Surface_CMy_Inv = NULL; Surface_CMz_Inv = NULL;
  Surface_CLift = NULL; Surface_CDrag = NULL;
  Surface_CMx = NULL; Surface_CMy = NULL; Surface_CMz = NULL;
  
  ForceInviscid = NULL; MomentInviscid = NULL;
  CPressure = NULL; CPressureTarget = NULL; HeatFlux = NULL; HeatFluxTarget = NULL; YPlus = NULL;
  
  CMerit_Inv = NULL;  CT_Inv = NULL;  CQ_Inv = NULL;
  
  CEquivArea_Inv = NULL;  CNearFieldOF_Inv = NULL;
  
  FanFace_MassFlow = NULL;  Exhaust_MassFlow = NULL;  Exhaust_Area = NULL;
  FanFace_Pressure = NULL;  FanFace_Mach = NULL;  FanFace_Area = NULL;
  
  iPoint_UndLapl = NULL;  jPoint_UndLapl = NULL;
  LowMach_Precontioner = NULL;
  Primitive = NULL; Primitive_i = NULL; Primitive_j = NULL;
  CharacPrimVar = NULL;
  Cauchy_Serie = NULL;
  
  /*--- Set the gamma value ---*/
  
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
  /*--- Define geometry constants in the solver structure
   Compressible flow, primitive variables nDim+7, (T,vx,vy,vz,P,rho,h,c,lamMu,EddyMu).
   Incompressible flow, primitive variables nDim+5, (P,vx,vy,vz,rho,beta,lamMu,EddyMu).
   FreeSurface Incompressible flow, primitive variables nDim+7, (P,vx,vy,vz,rho,beta,lamMu,EddyMu,LevelSet,Dist).
   ---*/
  
  nDim = geometry->GetnDim();
  nVar = nDim+2; nPrimVar = nDim+7; nPrimVarGrad = nDim+4;
  nMarker      = config->GetnMarker_All();
  nPoint       = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();
  
  /*--- Allocate the node variables (contiguous in memory) ---*/
  
  node = new CVariable*[nPoint];

  /*--- Define some auxiliary vectors related to the residual ---*/
  
  Residual      = new double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual[iVar]      = 0.0;
  Residual_RMS  = new double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar]  = 0.0;
  Residual_Max  = new double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual_Max[iVar]  = 0.0;
  Point_Max     = new unsigned long[nVar];  for (iVar = 0; iVar < nVar; iVar++) Point_Max[iVar]     = 0;
  Residual_i    = new double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual_i[iVar]    = 0.0;
  Residual_j    = new double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual_j[iVar]    = 0.0;
  Res_Conv      = new double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Res_Conv[iVar]      = 0.0;
  Res_Visc      = new double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Res_Visc[iVar]      = 0.0;
  Res_Sour      = new double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Res_Sour[iVar]      = 0.0;
  
  /*--- Define some auxiliary vectors related to the solution ---*/
  
  Solution   = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution[iVar]   = 0.0;
  Solution_i = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution_i[iVar] = 0.0;
  Solution_j = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution_j[iVar] = 0.0;
  
  /*--- Define some auxiliary vectors related to the geometry ---*/
  
  Vector   = new double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector[iDim]   = 0.0;
  Vector_i = new double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector_i[iDim] = 0.0;
  Vector_j = new double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector_j[iDim] = 0.0;
  
  /*--- Define some auxiliary vectors related to the primitive solution ---*/
  
  Primitive   = new double[nPrimVar]; for (iVar = 0; iVar < nPrimVar; iVar++) Primitive[iVar]   = 0.0;
  Primitive_i = new double[nPrimVar]; for (iVar = 0; iVar < nPrimVar; iVar++) Primitive_i[iVar] = 0.0;
  Primitive_j = new double[nPrimVar]; for (iVar = 0; iVar < nPrimVar; iVar++) Primitive_j[iVar] = 0.0;
  
  /*--- Define some auxiliary vectors related to the undivided lapalacian ---*/
  
  if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED) {
    iPoint_UndLapl = new double [nPoint];
    jPoint_UndLapl = new double [nPoint];
  }
  
  /*--- Define some auxiliary vectors related to low-speed preconditioning ---*/
  
  if (roe_turkel) {
    LowMach_Precontioner = new double* [nVar];
    for (iVar = 0; iVar < nVar; iVar ++)
      LowMach_Precontioner[iVar] = new double[nVar];
  }
  
  /*--- Initialize the solution and right hand side vectors for storing
   the residuals and updating the solution (always needed even for
   explicit schemes). ---*/
  
  LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
  
  /*--- Jacobians and vector structures for implicit computations ---*/
  
  if (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT) {
    
    Jacobian_i = new double* [nVar];
    Jacobian_j = new double* [nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Jacobian_i[iVar] = new double [nVar];
      Jacobian_j[iVar] = new double [nVar];
    }
    
    if (rank == MASTER_NODE) cout << "Initialize jacobian structure (Euler). MG level: " << iMesh <<"." << endl;
    Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config);
    
    if ((config->GetKind_Linear_Solver_Prec() == LINELET) ||
        (config->GetKind_Linear_Solver() == SMOOTHER_LINELET)) {
      nLineLets = Jacobian.BuildLineletPreconditioner(geometry, config);
      if (rank == MASTER_NODE) cout << "Compute linelet structure. " << nLineLets << " elements in each line (average)." << endl;
    }
    
  } else {
    if (rank == MASTER_NODE) cout << "Explicit scheme. No jacobian structure (Euler). MG level: " << iMesh <<"." << endl;
  }
  
  /*--- Define some auxiliary vectors for computing flow variable gradients by least squares ---*/
  
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
    
    /*--- S matrix := inv(R)*traspose(inv(R)) ---*/
    
    Smatrix = new double* [nDim];
    for (iDim = 0; iDim < nDim; iDim++)
      Smatrix[iDim] = new double [nDim];
    
    /*--- c vector := transpose(WA)*(Wb) ---*/
    
    cvector = new double* [nPrimVarGrad];
    for (iVar = 0; iVar < nPrimVarGrad; iVar++)
      cvector[iVar] = new double [nDim];
    
  }
  
  /*--- Store the value of the characteristic primitive variables at the boundaries ---*/
  
  CharacPrimVar = new double** [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    CharacPrimVar[iMarker] = new double* [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      CharacPrimVar[iMarker][iVertex] = new double [nPrimVar];
      for (iVar = 0; iVar < nPrimVar; iVar++) {
        CharacPrimVar[iMarker][iVertex][iVar] = 0.0;
      }
    }
  }
  
  /*--- Force definition and coefficient arrays for all of the markers ---*/
  
  CPressure = new double* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    CPressure[iMarker] = new double [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      CPressure[iMarker][iVertex] = 0.0;
    }
  }
  
  /*--- Force definition and coefficient arrays for all of the markers ---*/
  
  CPressureTarget = new double* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    CPressureTarget[iMarker] = new double [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      CPressureTarget[iMarker][iVertex] = 0.0;
    }
  }
  
  /*--- Non-dimensional coefficients ---*/
  
  ForceInviscid     = new double[nDim];
  MomentInviscid    = new double[3];
  CDrag_Inv         = new double[nMarker];
  CLift_Inv         = new double[nMarker];
  CSideForce_Inv    = new double[nMarker];
  CMx_Inv           = new double[nMarker];
  CMy_Inv           = new double[nMarker];
  CMz_Inv           = new double[nMarker];
  CEff_Inv          = new double[nMarker];
  CFx_Inv           = new double[nMarker];
  CFy_Inv           = new double[nMarker];
  CFz_Inv           = new double[nMarker];
  Surface_CLift_Inv = new double[config->GetnMarker_Monitoring()];
  Surface_CDrag_Inv = new double[config->GetnMarker_Monitoring()];
  Surface_CMx_Inv   = new double[config->GetnMarker_Monitoring()];
  Surface_CMy_Inv   = new double[config->GetnMarker_Monitoring()];
  Surface_CMz_Inv   = new double[config->GetnMarker_Monitoring()];
  Surface_CLift     = new double[config->GetnMarker_Monitoring()];
  Surface_CDrag     = new double[config->GetnMarker_Monitoring()];
  Surface_CMx       = new double[config->GetnMarker_Monitoring()];
  Surface_CMy       = new double[config->GetnMarker_Monitoring()];
  Surface_CMz       = new double[config->GetnMarker_Monitoring()];
  
  /*--- Rotorcraft coefficients ---*/
  
  CT_Inv           = new double[nMarker];
  CQ_Inv           = new double[nMarker];
  CMerit_Inv       = new double[nMarker];
  
  /*--- Supersonic coefficients ---*/
  
  CEquivArea_Inv   = new double[nMarker];
  CNearFieldOF_Inv = new double[nMarker];
  
  /*--- Engine simulation ---*/
  
  FanFace_MassFlow  = new double[nMarker];
  Exhaust_MassFlow  = new double[nMarker];
  Exhaust_Area      = new double[nMarker];
  FanFace_Pressure  = new double[nMarker];
  FanFace_Mach      = new double[nMarker];
  FanFace_Area      = new double[nMarker];
  
  /*--- Init total coefficients ---*/
  
  Total_CDrag = 0.0;    Total_CLift = 0.0;      Total_CSideForce = 0.0;
  Total_CMx = 0.0;      Total_CMy = 0.0;        Total_CMz = 0.0;
  Total_CEff = 0.0;     Total_CEquivArea = 0.0; Total_CNearFieldOF = 0.0;
  Total_CFx = 0.0;      Total_CFy = 0.0;        Total_CFz = 0.0;
  Total_CT = 0.0;       Total_CQ = 0.0;         Total_CMerit = 0.0;
  Total_MaxHeat = 0.0;  Total_Heat = 0.0;
  Total_CpDiff = 0.0;   Total_HeatFluxDiff = 0.0;
  
  /*--- Read farfield conditions ---*/
  
  Density_Inf  = config->GetDensity_FreeStreamND();
  Pressure_Inf = config->GetPressure_FreeStreamND();
  Velocity_Inf = config->GetVelocity_FreeStreamND();
  Energy_Inf   = config->GetEnergy_FreeStreamND();
  Mach_Inf     = config->GetMach_FreeStreamND();
  
  /*--- Initializate fan face pressure, fan face mach number, and mass flow rate ---*/
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    FanFace_MassFlow[iMarker] = 0.0;
    Exhaust_MassFlow[iMarker] = 0.0;
    FanFace_Mach[iMarker] = Mach_Inf;
    FanFace_Pressure[iMarker] = Pressure_Inf;
    FanFace_Area[iMarker] = 0.0;
    Exhaust_Area[iMarker] = 0.0;
  }
  
  /*--- Initialize the cauchy critera array for fixed CL mode ---*/
  
  if (config->GetFixed_CL_Mode())
    Cauchy_Serie = new double [config->GetCauchy_Elems()+1];
  
  /*--- Check for a restart and set up the variables at each node
   appropriately. Coarse multigrid levels will be intitially set to
   the farfield values bc the solver will immediately interpolate
   the solution from the finest mesh to the coarser levels. ---*/
  if (!restart || geometry->GetFinestMGLevel() == false || nZone > 1) {
    
    /*--- Restart the solution from the free-stream state ---*/
    for (iPoint = 0; iPoint < nPoint; iPoint++)
      node[iPoint] = new CEulerVariable(Density_Inf, Velocity_Inf, Energy_Inf, nDim, nVar, config);
  }
  
  else {
    
    /*--- Initialize the solution from the restart file information ---*/
    ifstream restart_file;
    string filename = config->GetSolution_FlowFileName();
    
    /*--- Modify file name for an unsteady restart ---*/
    if (dual_time) {
      int Unst_RestartIter;
      if (adjoint) {
        Unst_RestartIter = int(config->GetUnst_AdjointIter()) - 1;
      } else if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
        Unst_RestartIter = int(config->GetUnst_RestartIter())-1;
      else
        Unst_RestartIter = int(config->GetUnst_RestartIter())-2;
      filename = config->GetUnsteady_FileName(filename, Unst_RestartIter);
    }
    
    /*--- Open the restart file, throw an error if this fails. ---*/
    restart_file.open(filename.data(), ios::in);
    if (restart_file.fail()) {
      cout << "There is no flow restart file!! " << filename.data() << "."<< endl;
      exit(1);
    }
    
    /*--- In case this is a parallel simulation, we need to perform the
     Global2Local index transformation first. ---*/
    long *Global2Local = new long[geometry->GetGlobal_nPointDomain()];
    
    /*--- First, set all indices to a negative value by default ---*/
    for(iPoint = 0; iPoint < geometry->GetGlobal_nPointDomain(); iPoint++)
      Global2Local[iPoint] = -1;
    
    /*--- Now fill array with the transform values only for local points ---*/
    for(iPoint = 0; iPoint < nPointDomain; iPoint++)
      Global2Local[geometry->node[iPoint]->GetGlobalIndex()] = iPoint;
    
    /*--- Read all lines in the restart file ---*/
    long iPoint_Local; unsigned long iPoint_Global = 0; string text_line;
    
    /*--- The first line is the header ---*/
    getline (restart_file, text_line);
    
    while (getline (restart_file, text_line)) {
      istringstream point_line(text_line);
      
      /*--- Retrieve local index. If this node from the restart file lives
       on a different processor, the value of iPoint_Local will be -1.
       Otherwise, the local index for this node on the current processor
       will be returned and used to instantiate the vars. ---*/
      iPoint_Local = Global2Local[iPoint_Global];
      
      /*--- Load the solution for this node. Note that the first entry
       on the restart file line is the global index, followed by the
       node coordinates, and then the conservative variables. ---*/
      if (iPoint_Local >= 0) {
        if (nDim == 2) point_line >> index >> dull_val >> dull_val >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3];
        if (nDim == 3) point_line >> index >> dull_val >> dull_val >> dull_val >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3] >> Solution[4];
        node[iPoint_Local] = new CEulerVariable(Solution, nDim, nVar, config);
      }
      iPoint_Global++;
    }
    
    /*--- Instantiate the variable class with an arbitrary solution
     at any halo/periodic nodes. The initial solution can be arbitrary,
     because a send/recv is performed immediately in the solver. ---*/
    for(iPoint = nPointDomain; iPoint < nPoint; iPoint++)
      node[iPoint] = new CEulerVariable(Solution, nDim, nVar, config);
    
    /*--- Close the restart file ---*/
    restart_file.close();
    
    /*--- Free memory needed for the transformation ---*/
    delete [] Global2Local;
  }
  
  /*--- Check that the initial solution is physical, report any non-physical nodes ---*/
  counter_local = 0;
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    Density = node[iPoint]->GetSolution(0);
    Velocity2 = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      Velocity2 += (node[iPoint]->GetSolution(iDim+1)/Density)*(node[iPoint]->GetSolution(iDim+1)/Density);
    Pressure    = Gamma_Minus_One*Density*(node[iPoint]->GetSolution(nDim+1)/Density-0.5*Velocity2);
    Temperature = Pressure / ( Gas_Constant * Density);
    if ((Pressure < 0.0) || (Temperature < 0.0)) {
      Solution[0] = Density_Inf;
      for (iDim = 0; iDim < nDim; iDim++)
        Solution[iDim+1] = Velocity_Inf[iDim]*Density_Inf;
      Solution[nDim+1] = Energy_Inf*Density_Inf;
      node[iPoint]->SetSolution(Solution);
      node[iPoint]->SetSolution_Old(Solution);
      counter_local++;
    }
  }
  counter_global = counter_local;

  if ((rank == MASTER_NODE) && (counter_global != 0))
    cout << "Warning. The original solution contains "<< counter_global << " points that are not physical." << endl;
  
  /*--- Define solver parameters needed for execution of destructor ---*/
  if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED ) space_centered = true;
  else space_centered = false;
  
  if (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT) euler_implicit = true;
  else euler_implicit = false;
  
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) least_squares = true;
  else least_squares = false;
  
}

CEulerSolver::~CEulerSolver(void) {
  unsigned short iVar, iMarker;
  
  /*--- Array deallocation ---*/
  if (CDrag_Inv != NULL)         delete [] CDrag_Inv;
  if (CLift_Inv != NULL)         delete [] CLift_Inv;
  if (CSideForce_Inv != NULL)    delete [] CSideForce_Inv;
  if (CMx_Inv != NULL)           delete [] CMx_Inv;
  if (CMy_Inv != NULL)           delete [] CMy_Inv;
  if (CMz_Inv != NULL)           delete [] CMz_Inv;
  if (CFx_Inv != NULL)           delete [] CFx_Inv;
  if (CFy_Inv != NULL)           delete [] CFy_Inv;
  if (CFz_Inv != NULL)           delete [] CFz_Inv;
  if (Surface_CLift_Inv != NULL) delete[] Surface_CLift_Inv;
  if (Surface_CDrag_Inv != NULL) delete[] Surface_CDrag_Inv;
  if (Surface_CMx_Inv != NULL)  delete [] Surface_CMx_Inv;
  if (Surface_CMy_Inv != NULL)  delete [] Surface_CMy_Inv;
  if (Surface_CMz_Inv != NULL)  delete [] Surface_CMz_Inv;
  if (Surface_CLift != NULL)    delete [] Surface_CLift;
  if (Surface_CDrag != NULL)    delete [] Surface_CDrag;
  if (Surface_CMx != NULL)      delete [] Surface_CMx;
  if (Surface_CMy != NULL)      delete [] Surface_CMy;
  if (Surface_CMz != NULL)      delete [] Surface_CMz;
  if (CEff_Inv != NULL)          delete [] CEff_Inv;
  if (CMerit_Inv != NULL)        delete [] CMerit_Inv;
  if (CT_Inv != NULL)            delete [] CT_Inv;
  if (CQ_Inv != NULL)            delete [] CQ_Inv;
  if (CEquivArea_Inv != NULL)    delete [] CEquivArea_Inv;
  if (CNearFieldOF_Inv != NULL)  delete [] CNearFieldOF_Inv;
  if (ForceInviscid != NULL)     delete [] ForceInviscid;
  if (MomentInviscid != NULL)    delete [] MomentInviscid;
  if (FanFace_MassFlow != NULL)  delete [] FanFace_MassFlow;
  if (Exhaust_MassFlow != NULL)  delete [] Exhaust_MassFlow;
  if (Exhaust_Area != NULL)      delete [] Exhaust_Area;
  if (FanFace_Pressure != NULL)  delete [] FanFace_Pressure;
  if (FanFace_Mach != NULL)      delete [] FanFace_Mach;
  if (FanFace_Area != NULL)      delete [] FanFace_Area;
  if (iPoint_UndLapl != NULL)       delete [] iPoint_UndLapl;
  if (jPoint_UndLapl != NULL)       delete [] jPoint_UndLapl;
  if (Primitive != NULL)        delete [] Primitive;
  if (Primitive_i != NULL)      delete [] Primitive_i;
  if (Primitive_j != NULL)      delete [] Primitive_j;
  
  if (LowMach_Precontioner != NULL) {
    for (iVar = 0; iVar < nVar; iVar ++)
      delete LowMach_Precontioner[iVar];
    delete [] LowMach_Precontioner;
  }
  
  if (CPressure != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      delete CPressure[iMarker];
    delete [] CPressure;
  }
  
  if (CPressureTarget != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      delete CPressureTarget[iMarker];
    delete [] CPressureTarget;
  }
  
  //  if (CharacPrimVar != NULL) {
  //    for (iMarker = 0; iMarker < nMarker; iMarker++) {
  //      for (iVertex = 0; iVertex < nVertex; iVertex++) {
  //        delete CharacPrimVar[iMarker][iVertex];
  //      }
  //    }
  //    delete [] CharacPrimVar;
  //  }
  
  if (HeatFlux != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      delete HeatFlux[iMarker];
    }
    delete [] HeatFlux;
  }
  
  if (HeatFluxTarget != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      delete HeatFluxTarget[iMarker];
    }
    delete [] HeatFluxTarget;
  }
  
  if (YPlus != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      delete YPlus[iMarker];
    }
    delete [] YPlus;
  }
  
  if (Cauchy_Serie != NULL)
    delete [] Cauchy_Serie;
  
}

void CEulerSolver::SetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long ExtIter) {
  unsigned long iPoint, Point_Fine;
  unsigned short iMesh, iChildren, iVar;
  double Area_Children, Area_Parent, *Solution_Fine, *Solution;
  
  bool restart = (config->GetRestart() || config->GetRestart_Flow());
  bool rans = ((config->GetKind_Solver() == RANS) ||
               (config->GetKind_Solver() == ADJ_RANS));
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  
  /*--- If restart solution, then interpolate the flow solution to
   all the multigrid levels, this is important with the dual time strategy ---*/
  
  if (restart && (ExtIter == 0)) {
    Solution = new double[nVar];
    for (iMesh = 1; iMesh <= config->GetMGLevels(); iMesh++) {
      for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
        Area_Parent = geometry[iMesh]->node[iPoint]->GetVolume();
        for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = 0.0;
        for (iChildren = 0; iChildren < geometry[iMesh]->node[iPoint]->GetnChildren_CV(); iChildren++) {
          Point_Fine = geometry[iMesh]->node[iPoint]->GetChildren_CV(iChildren);
          Area_Children = geometry[iMesh-1]->node[Point_Fine]->GetVolume();
          Solution_Fine = solver_container[iMesh-1][FLOW_SOL]->node[Point_Fine]->GetSolution();
          for (iVar = 0; iVar < nVar; iVar++) {
            Solution[iVar] += Solution_Fine[iVar]*Area_Children/Area_Parent;
          }
        }
        solver_container[iMesh][FLOW_SOL]->node[iPoint]->SetSolution(Solution);
      }
    }
    delete [] Solution;
    
    /*--- Interpolate the turblence variable also, if needed ---*/
    if (rans) {
      unsigned short nVar_Turb = solver_container[MESH_0][TURB_SOL]->GetnVar();
      Solution = new double[nVar_Turb];
      for (iMesh = 1; iMesh <= config->GetMGLevels(); iMesh++) {
        for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
          Area_Parent = geometry[iMesh]->node[iPoint]->GetVolume();
          for (iVar = 0; iVar < nVar_Turb; iVar++) Solution[iVar] = 0.0;
          for (iChildren = 0; iChildren < geometry[iMesh]->node[iPoint]->GetnChildren_CV(); iChildren++) {
            Point_Fine = geometry[iMesh]->node[iPoint]->GetChildren_CV(iChildren);
            Area_Children = geometry[iMesh-1]->node[Point_Fine]->GetVolume();
            Solution_Fine = solver_container[iMesh-1][TURB_SOL]->node[Point_Fine]->GetSolution();
            for (iVar = 0; iVar < nVar_Turb; iVar++) {
              Solution[iVar] += Solution_Fine[iVar]*Area_Children/Area_Parent;
            }
          }
          solver_container[iMesh][TURB_SOL]->node[iPoint]->SetSolution(Solution);
        }
        solver_container[iMesh][TURB_SOL]->Postprocessing(geometry[iMesh], solver_container[iMesh], config, iMesh);
      }
      delete [] Solution;
    }
  }
  
  
  /*--- The value of the solution for the first iteration of the dual time ---*/
  
  if (dual_time && (ExtIter == 0 || (restart && ExtIter == config->GetUnst_RestartIter()))) {
    
    /*--- Push back the initial condition to previous solution containers
     for a 1st-order restart or when simply intitializing to freestream. ---*/
    for (iMesh = 0; iMesh <= config->GetMGLevels(); iMesh++) {
      for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
        solver_container[iMesh][FLOW_SOL]->node[iPoint]->Set_Solution_time_n();
        solver_container[iMesh][FLOW_SOL]->node[iPoint]->Set_Solution_time_n1();
        if (rans) {
          solver_container[iMesh][TURB_SOL]->node[iPoint]->Set_Solution_time_n();
          solver_container[iMesh][TURB_SOL]->node[iPoint]->Set_Solution_time_n1();
        }
      }
    }
    
    if ((restart && ExtIter == config->GetUnst_RestartIter()) &&
        (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)) {
      
      /*--- Load an additional restart file for a 2nd-order restart ---*/
      solver_container[MESH_0][FLOW_SOL]->LoadRestart(geometry, solver_container, config, int(config->GetUnst_RestartIter()-1));
      
      /*--- Load an additional restart file for the turbulence model ---*/
      if (rans)
        solver_container[MESH_0][TURB_SOL]->LoadRestart(geometry, solver_container, config, int(config->GetUnst_RestartIter()-1));
      
      /*--- Push back this new solution to time level N. ---*/
      for (iMesh = 0; iMesh <= config->GetMGLevels(); iMesh++) {
        for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
          solver_container[iMesh][FLOW_SOL]->node[iPoint]->Set_Solution_time_n();
          if (rans) {
            solver_container[iMesh][TURB_SOL]->node[iPoint]->Set_Solution_time_n();
          }
        }
      }
    }
  }
  
}

void CEulerSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {
  
  unsigned long iPoint, ErrorCounter = 0;
  bool RightSol;
  int rank;
  
  rank = MASTER_NODE;
  
  unsigned long ExtIter = config->GetExtIter();
  bool adjoint        = config->GetAdjoint();
  bool implicit       = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool second_order   = ((config->GetSpatialOrder_Flow() == SECOND_ORDER) || (config->GetSpatialOrder_Flow() == SECOND_ORDER_LIMITER));
  bool limiter          = ((config->GetSpatialOrder_Flow() == SECOND_ORDER_LIMITER) && (ExtIter <= config->GetLimiterIter()));
  bool center         = (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED) || (adjoint && config->GetKind_ConvNumScheme_AdjFlow() == SPACE_CENTERED);
  bool center_jst     = center && (config->GetKind_Centered_Flow() == JST);
  
  for (iPoint = 0; iPoint < nPoint; iPoint ++) {
    
    /*--- Incompressible flow, primitive variables nDim+3, (P,vx,vy,vz,rho,beta),
     FreeSurface Incompressible flow, primitive variables nDim+5, (P,vx,vy,vz,rho,beta,LevelSet,Dist),
     Compressible flow, primitive variables nDim+5, (T,vx,vy,vz,P,rho,h,c) ---*/
    
    RightSol = node[iPoint]->SetPrimVar_Compressible(config);
    if (!RightSol) ErrorCounter++;
    
  }
  
  for (iPoint = 0; iPoint < nPoint; iPoint ++) {
    
    /*--- Initialize the residual vector (except for output of the residuals) ---*/
    
    if (!Output) LinSysRes.SetBlock_Zero(iPoint);
    
  }
  
  /*--- Upwind second order reconstruction ---*/
  
  if ((second_order && !center) && (iMesh == MESH_0)) {
    
    /*--- Gradient computation ---*/
    
    if (config->GetKind_Gradient_Method() == GREEN_GAUSS) {
      SetPrimVar_Gradient_GG(geometry, config);
    }
    
    if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
      SetPrimVar_Gradient_LS(geometry, config);
    }
    
    /*--- Limiter computation ---*/
    
    if ((limiter) && (iMesh == MESH_0)) {
      SetPrimVar_Limiter(geometry, config);
    }
    
  }
  
  /*--- Artificial dissipation ---*/
  
  if (center) {
    SetMax_Eigenvalue(geometry, config);
    
    if ((center_jst) && (iMesh == MESH_0)) {
      SetDissipation_Switch(geometry, config);
      
      SetUndivided_Laplacian(geometry, config);
    }
  }
  
  /*--- Initialize the jacobian matrices ---*/
  
  if (implicit) {
    Jacobian.SetValZero();
  }
  
  /*--- Error message ---*/
  if (Output && (ErrorCounter >= 10) && (rank == MASTER_NODE) && (iMesh == MESH_0))
    cout <<"The solution contains "<< ErrorCounter << " non-physical points." << endl;
  
}

void CEulerSolver::Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                  unsigned short iMesh) { }

void CEulerSolver::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                unsigned short iMesh, unsigned long Iteration) {
  double *Normal, Area, Vol, Mean_SoundSpeed, Mean_ProjVel, Lambda, Local_Delta_Time,
  Global_Delta_Time = 1E6, Global_Delta_UnstTimeND, ProjVel, ProjVel_i, ProjVel_j;
  unsigned long iEdge, iVertex, iPoint, jPoint;
  unsigned short iDim, iMarker;
  
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool grid_movement = config->GetGrid_Movement();
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  
  Min_Delta_Time = 1.E6; Max_Delta_Time = 0.0;
  
  /*--- Set maximum inviscid eigenvalue to zero, and compute sound speed ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++)
    node[iPoint]->SetMax_Lambda_Inv(0.0);
  
  /*--- Loop interior edges ---*/
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    /*--- Point identification, Normal vector and area ---*/
    
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    
    Normal = geometry->edge[iEdge]->GetNormal();
    Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);
    
    /*--- Mean Values ---*/
    
    Mean_ProjVel = 0.5 * (node[iPoint]->GetProjVel(Normal) + node[jPoint]->GetProjVel(Normal));
    Mean_SoundSpeed = 0.5 * (node[iPoint]->GetSoundSpeed() + node[jPoint]->GetSoundSpeed()) * Area;
    
    /*--- Adjustment for grid movement ---*/
    
    if (grid_movement) {
      double *GridVel_i = geometry->node[iPoint]->GetGridVel();
      double *GridVel_j = geometry->node[jPoint]->GetGridVel();
      ProjVel_i = 0.0; ProjVel_j = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        ProjVel_i += GridVel_i[iDim]*Normal[iDim];
        ProjVel_j += GridVel_j[iDim]*Normal[iDim];
      }
      Mean_ProjVel -= 0.5 * (ProjVel_i + ProjVel_j);
    }
    
    /*--- Inviscid contribution ---*/
    
    Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed;
    if (geometry->node[iPoint]->GetDomain()) node[iPoint]->AddMax_Lambda_Inv(Lambda);
    if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddMax_Lambda_Inv(Lambda);
    
  }
  
  /*--- Loop boundary edges ---*/
  
  for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
    for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
      
      /*--- Point identification, Normal vector and area ---*/
      
      iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
      Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
      Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);
      
      /*--- Mean Values ---*/
      
      Mean_ProjVel = node[iPoint]->GetProjVel(Normal);
      Mean_SoundSpeed = node[iPoint]->GetSoundSpeed() * Area;
      
      /*--- Adjustment for grid movement ---*/
      
      if (grid_movement) {
        double *GridVel = geometry->node[iPoint]->GetGridVel();
        ProjVel = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          ProjVel += GridVel[iDim]*Normal[iDim];
        Mean_ProjVel -= ProjVel;
      }
      
      /*--- Inviscid contribution ---*/
      
      Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed;
      if (geometry->node[iPoint]->GetDomain()) {
        node[iPoint]->AddMax_Lambda_Inv(Lambda);
      }
      
    }
  }
  
  /*--- Each element uses their own speed, steady state simulation ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    
    Vol = geometry->node[iPoint]->GetVolume();
    
    if (Vol != 0.0) {
      Local_Delta_Time = config->GetCFL(iMesh)*Vol / node[iPoint]->GetMax_Lambda_Inv();
      Global_Delta_Time = min(Global_Delta_Time, Local_Delta_Time);
      Min_Delta_Time = min(Min_Delta_Time, Local_Delta_Time);
      Max_Delta_Time = max(Max_Delta_Time, Local_Delta_Time);
      if (Local_Delta_Time > config->GetMax_DeltaTime())
        Local_Delta_Time = config->GetMax_DeltaTime();
      node[iPoint]->SetDelta_Time(Local_Delta_Time);
    }
    else {
      node[iPoint]->SetDelta_Time(0.0);
    }
    
  }
  
  /*--- For exact time solution use the minimum delta time of the whole mesh ---*/
  if (config->GetUnsteady_Simulation() == TIME_STEPPING) {
    for(iPoint = 0; iPoint < nPointDomain; iPoint++)
      node[iPoint]->SetDelta_Time(Global_Delta_Time);
  }
  
  /*--- Recompute the unsteady time step for the dual time strategy
   if the unsteady CFL is diferent from 0 ---*/
  if ((dual_time) && (Iteration == 0) && (config->GetUnst_CFL() != 0.0) && (iMesh == MESH_0)) {
    Global_Delta_UnstTimeND = config->GetUnst_CFL()*Global_Delta_Time/config->GetCFL(iMesh);
    
    config->SetDelta_UnstTimeND(Global_Delta_UnstTimeND);
  }
  
  /*--- The pseudo local time (explicit integration) cannot be greater than the physical time ---*/
  if (dual_time)
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      if (!implicit) {
        Local_Delta_Time = min((2.0/3.0)*config->GetDelta_UnstTimeND(), node[iPoint]->GetDelta_Time());
        node[iPoint]->SetDelta_Time(Local_Delta_Time);
      }
    }
  
}

void CEulerSolver::Centered_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                     CConfig *config, unsigned short iMesh, unsigned short iRKStep) {
  
  unsigned long iEdge, iPoint, jPoint;
  
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool second_order = ((config->GetKind_Centered_Flow() == JST) && (iMesh == MESH_0));
  bool grid_movement = config->GetGrid_Movement();
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    /*--- Points in edge, set normal vectors, and number of neighbors ---*/
    
    iPoint = geometry->edge[iEdge]->GetNode(0); jPoint = geometry->edge[iEdge]->GetNode(1);
    numerics->SetNormal(geometry->edge[iEdge]->GetNormal());
    numerics->SetNeighbor(geometry->node[iPoint]->GetnNeighbor(), geometry->node[jPoint]->GetnNeighbor());
    
    /*--- Set primitive variables w/o reconstruction ---*/
    
    numerics->SetPrimitive(node[iPoint]->GetPrimVar(), node[jPoint]->GetPrimVar());
    
    /*--- Set the largest convective eigenvalue ---*/
    
    numerics->SetLambda(node[iPoint]->GetLambda(), node[jPoint]->GetLambda());
    
    /*--- Set undivided laplacian an pressure based sensor ---*/
    
    if (second_order) {
      numerics->SetUndivided_Laplacian(node[iPoint]->GetUndivided_Laplacian(), node[jPoint]->GetUndivided_Laplacian());
      numerics->SetSensor(node[iPoint]->GetSensor(), node[jPoint]->GetSensor());
    }
    
    /*--- Grid movement ---*/
    
    if (grid_movement) {
      numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[jPoint]->GetGridVel());
    }
    
    /*--- Compute residuals, and jacobians ---*/
    
    numerics->ComputeResidual(Res_Conv, Jacobian_i, Jacobian_j, config);
    
    /*--- Update convective and artificial dissipation residuals ---*/
    
    LinSysRes.AddBlock(iPoint, Res_Conv);
    LinSysRes.SubtractBlock(jPoint, Res_Conv);
    
    /*--- Set implicit computation ---*/
    if (implicit) {
      Jacobian.AddBlock(iPoint,iPoint,Jacobian_i);
      Jacobian.AddBlock(iPoint,jPoint,Jacobian_j);
      Jacobian.SubtractBlock(jPoint,iPoint,Jacobian_i);
      Jacobian.SubtractBlock(jPoint,jPoint,Jacobian_j);
    }
  }
  
}

void CEulerSolver::Upwind_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                   CConfig *config, unsigned short iMesh) {
  double **Gradient_i, **Gradient_j, Project_Grad_i, Project_Grad_j, *V_i, *V_j, *Limiter_i = NULL,
  *Limiter_j = NULL, sqvel;
  unsigned long iEdge, iPoint, jPoint;
  unsigned short iDim, iVar;
  
  bool implicit         = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool second_order     = (((config->GetSpatialOrder_Flow() == SECOND_ORDER) || (config->GetSpatialOrder_Flow() == SECOND_ORDER_LIMITER)) && (iMesh == MESH_0));
  bool limiter          = (config->GetSpatialOrder_Flow() == SECOND_ORDER_LIMITER);
  bool grid_movement    = config->GetGrid_Movement();
  bool roe_turkel       = (config->GetKind_Upwind_Flow() == TURKEL);
  
  for(iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    /*--- Points in edge and normal vectors ---*/
    
    iPoint = geometry->edge[iEdge]->GetNode(0); jPoint = geometry->edge[iEdge]->GetNode(1);
    numerics->SetNormal(geometry->edge[iEdge]->GetNormal());
    
    /*--- Roe Turkel preconditioning ---*/
    
    if (roe_turkel) {
      sqvel = 0.0;
      for (iDim = 0; iDim < nDim; iDim ++)
        sqvel += config->GetVelocity_FreeStream()[iDim]*config->GetVelocity_FreeStream()[iDim];
      numerics->SetVelocity2_Inf(sqvel);
    }
    
    /*--- Grid movement ---*/
    
    if (grid_movement)
      numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[jPoint]->GetGridVel());
    
    /*--- Get primitive variables ---*/
    
    V_i = node[iPoint]->GetPrimVar(); V_j = node[jPoint]->GetPrimVar();
    
    /*--- High order reconstruction using MUSCL strategy ---*/
    
    if (second_order) {
      
      for (iDim = 0; iDim < nDim; iDim++) {
        Vector_i[iDim] = 0.5*(geometry->node[jPoint]->GetCoord(iDim) - geometry->node[iPoint]->GetCoord(iDim));
        Vector_j[iDim] = 0.5*(geometry->node[iPoint]->GetCoord(iDim) - geometry->node[jPoint]->GetCoord(iDim));
      }
      
      Gradient_i = node[iPoint]->GetGradient_Primitive(); Gradient_j = node[jPoint]->GetGradient_Primitive();
      if (limiter) { Limiter_i = node[iPoint]->GetLimiter_Primitive(); Limiter_j = node[jPoint]->GetLimiter_Primitive(); }
      
      for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
        Project_Grad_i = 0.0; Project_Grad_j = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          Project_Grad_i += Vector_i[iDim]*Gradient_i[iVar][iDim];
          Project_Grad_j += Vector_j[iDim]*Gradient_j[iVar][iDim];
        }
        if (limiter) {
          Primitive_i[iVar] = V_i[iVar] + Limiter_i[iVar]*Project_Grad_i;
          Primitive_j[iVar] = V_j[iVar] + Limiter_j[iVar]*Project_Grad_j;
        }
        else {
          Primitive_i[iVar] = V_i[iVar] + Project_Grad_i;
          Primitive_j[iVar] = V_j[iVar] + Project_Grad_j;
        }
      }
      
      /*--- Set conservative variables with reconstruction ---*/
      
      numerics->SetPrimitive(Primitive_i, Primitive_j);
      
    } else {
      
      /*--- Set conservative variables without reconstruction ---*/
      
      numerics->SetPrimitive(V_i, V_j);
      
    }
    
    /*--- Compute the residual ---*/
    
    numerics->ComputeResidual(Res_Conv, Jacobian_i, Jacobian_j, config);
    
    /*--- Update residual value ---*/
    
    LinSysRes.AddBlock(iPoint, Res_Conv);
    LinSysRes.SubtractBlock(jPoint, Res_Conv);
    
    /*--- Set implicit jacobians ---*/
    
    if (implicit) {
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      Jacobian.AddBlock(iPoint, jPoint, Jacobian_j);
      Jacobian.SubtractBlock(jPoint, iPoint, Jacobian_i);
      Jacobian.SubtractBlock(jPoint, jPoint, Jacobian_j);
    }
    
    /*--- Roe Turkel preconditioning, set the value of beta ---*/
    
    if (roe_turkel) {
      node[iPoint]->SetPreconditioner_Beta(numerics->GetPrecond_Beta());
      node[jPoint]->SetPreconditioner_Beta(numerics->GetPrecond_Beta());
    }
    
  }
  
}

void CEulerSolver::Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics,
                                   CConfig *config, unsigned short iMesh) {
  
  unsigned short iVar;
  
  /*--- Initialize the source residual to zero ---*/
  for (iVar = 0; iVar < nVar; iVar++) Residual[iVar] = 0.0;
  
}

void CEulerSolver::SetMax_Eigenvalue(CGeometry *geometry, CConfig *config) {
  double *Normal, Area, Mean_SoundSpeed, Mean_ProjVel, Lambda,
  ProjVel, ProjVel_i, ProjVel_j;
  unsigned long iEdge, iVertex, iPoint, jPoint;
  unsigned short iDim, iMarker;
  
  bool grid_movement = config->GetGrid_Movement();
  
  /*--- Set maximum inviscid eigenvalue to zero, and compute sound speed ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    node[iPoint]->SetLambda(0.0);
  }
  
  /*--- Loop interior edges ---*/
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    /*--- Point identification, Normal vector and area ---*/
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    
    Normal = geometry->edge[iEdge]->GetNormal();
    Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);
    
    /*--- Mean Values ---*/
    Mean_ProjVel = 0.5 * (node[iPoint]->GetProjVel(Normal) + node[jPoint]->GetProjVel(Normal));
    Mean_SoundSpeed = 0.5 * (node[iPoint]->GetSoundSpeed() + node[jPoint]->GetSoundSpeed()) * Area;
    
    /*--- Adjustment for grid movement ---*/
    if (grid_movement) {
      double *GridVel_i = geometry->node[iPoint]->GetGridVel();
      double *GridVel_j = geometry->node[jPoint]->GetGridVel();
      ProjVel_i = 0.0; ProjVel_j =0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        ProjVel_i += GridVel_i[iDim]*Normal[iDim];
        ProjVel_j += GridVel_j[iDim]*Normal[iDim];
      }
      Mean_ProjVel -= 0.5 * (ProjVel_i + ProjVel_j);
    }
    
    /*--- Inviscid contribution ---*/
    Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed;
    if (geometry->node[iPoint]->GetDomain()) node[iPoint]->AddLambda(Lambda);
    if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddLambda(Lambda);
    
  }
  
  /*--- Loop boundary edges ---*/
  for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
    for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
      
      /*--- Point identification, Normal vector and area ---*/
      iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
      Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
      Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);
      
      /*--- Mean Values ---*/
      Mean_ProjVel = node[iPoint]->GetProjVel(Normal);
      Mean_SoundSpeed = node[iPoint]->GetSoundSpeed() * Area;
      
      /*--- Adjustment for grid movement ---*/
      if (grid_movement) {
        double *GridVel = geometry->node[iPoint]->GetGridVel();
        ProjVel = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          ProjVel += GridVel[iDim]*Normal[iDim];
        Mean_ProjVel -= ProjVel;
      }
      
      /*--- Inviscid contribution ---*/
      Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed;
      if (geometry->node[iPoint]->GetDomain()) {
        node[iPoint]->AddLambda(Lambda);
      }
      
    }
  }
  
}

void CEulerSolver::SetUndivided_Laplacian(CGeometry *geometry, CConfig *config) {
  unsigned long iPoint, jPoint, iEdge;
  double Pressure_i = 0, Pressure_j = 0, *Diff;
  unsigned short iVar;
  
  Diff = new double[nVar];
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++)
    node[iPoint]->SetUnd_LaplZero();
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    
    /*--- Solution differences ---*/
    for (iVar = 0; iVar < nVar; iVar++)
      Diff[iVar] = node[iPoint]->GetSolution(iVar) - node[jPoint]->GetSolution(iVar);
    
    /*--- Correction for compressible flows which use the enthalpy ---*/
    Pressure_i = node[iPoint]->GetPressure();
    Pressure_j = node[jPoint]->GetPressure();
    Diff[nVar-1] = (node[iPoint]->GetSolution(nVar-1) + Pressure_i) - (node[jPoint]->GetSolution(nVar-1) + Pressure_j);
    
#ifdef STRUCTURED_GRID
    
    if (geometry->node[iPoint]->GetDomain()) node[iPoint]->SubtractUnd_Lapl(Diff);
    if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddUnd_Lapl(Diff);
    
#else
    
    bool boundary_i = geometry->node[iPoint]->GetPhysicalBoundary();
    bool boundary_j = geometry->node[jPoint]->GetPhysicalBoundary();
    
    /*--- Both points inside the domain, or both in the boundary ---*/
    if ((!boundary_i && !boundary_j) || (boundary_i && boundary_j)) {
      if (geometry->node[iPoint]->GetDomain()) node[iPoint]->SubtractUnd_Lapl(Diff);
      if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddUnd_Lapl(Diff);
    }
    
    /*--- iPoint inside the domain, jPoint on the boundary ---*/
    if (!boundary_i && boundary_j)
      if (geometry->node[iPoint]->GetDomain()) node[iPoint]->SubtractUnd_Lapl(Diff);
    
    /*--- jPoint inside the domain, iPoint on the boundary ---*/
    if (boundary_i && !boundary_j)
      if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddUnd_Lapl(Diff);
    
#endif
    
  }
  
#ifdef STRUCTURED_GRID
  
  unsigned long Point_Normal = 0, iVertex;
  double Pressure_mirror = 0, *U_mirror;
  unsigned short iMarker;
  
  U_mirror = new double[nVar];
  
  /*--- Loop over all boundaries and include an extra contribution
   from a mirror node. Find the nearest normal, interior point
   for a boundary node and make a linear approximation. ---*/
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
    if (config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE &&
        config->GetMarker_All_KindBC(iMarker) != INTERFACE_BOUNDARY &&
        config->GetMarker_All_KindBC(iMarker) != NEARFIELD_BOUNDARY &&
        config->GetMarker_All_KindBC(iMarker) != PERIODIC_BOUNDARY) {
      
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        if (geometry->node[iPoint]->GetDomain()) {
          
          Point_Normal = geometry->vertex[iMarker][iVertex]->GetNormal_Neighbor();
          
          /*--- Interpolate & compute difference in the conserved variables ---*/
          for (iVar = 0; iVar < nVar; iVar++) {
            U_mirror[iVar] = 2.0*node[iPoint]->GetSolution(iVar) - node[Point_Normal]->GetSolution(iVar);
            Diff[iVar]   = node[iPoint]->GetSolution(iVar) - U_mirror[iVar];
          }
          
          /*--- Correction for compressible flows ---*/
          Pressure_mirror = 2.0*node[iPoint]->GetPressure() - node[Point_Normal]->GetPressure();
          Diff[nVar-1] = (node[iPoint]->GetSolution(nVar-1) + node[iPoint]->GetPressure()) - (U_mirror[nVar-1] + Pressure_mirror);
          
          /*--- Subtract contribution at the boundary node only ---*/
          node[iPoint]->SubtractUnd_Lapl(Diff);
        }
      }
    }
  }
  
  delete [] U_mirror;
  
#endif
  
  delete [] Diff;
  
}

void CEulerSolver::SetDissipation_Switch(CGeometry *geometry, CConfig *config) {
  unsigned long iEdge, iPoint, jPoint;
  double Pressure_i, Pressure_j;
  
  /*--- Reset variables to store the undivided pressure ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    iPoint_UndLapl[iPoint] = 0.0;
    jPoint_UndLapl[iPoint] = 0.0;
  }
  
  /*--- Evaluate the pressure sensor ---*/
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    
    /*--- Get the pressure, or density for incompressible solvers ---*/
    Pressure_i = node[iPoint]->GetPressure();
    Pressure_j = node[jPoint]->GetPressure();
    
#ifdef STRUCTURED_GRID
    
    /*--- Compute numerator and denominator ---*/
    if (geometry->node[iPoint]->GetDomain()) {
      iPoint_UndLapl[iPoint] += (Pressure_j - Pressure_i);
      jPoint_UndLapl[iPoint] += (Pressure_i + Pressure_j);
    }
    if (geometry->node[jPoint]->GetDomain()) {
      iPoint_UndLapl[jPoint] += (Pressure_i - Pressure_j);
      jPoint_UndLapl[jPoint] += (Pressure_i + Pressure_j);
    }
    
#else
    
    bool boundary_i = geometry->node[iPoint]->GetPhysicalBoundary();
    bool boundary_j = geometry->node[jPoint]->GetPhysicalBoundary();
    
    /*--- Both points inside the domain, or both on the boundary ---*/
    if ((!boundary_i && !boundary_j) || (boundary_i && boundary_j)){
      if (geometry->node[iPoint]->GetDomain()) { iPoint_UndLapl[iPoint] += (Pressure_j - Pressure_i); jPoint_UndLapl[iPoint] += (Pressure_i + Pressure_j); }
      if (geometry->node[jPoint]->GetDomain()) { iPoint_UndLapl[jPoint] += (Pressure_i - Pressure_j); jPoint_UndLapl[jPoint] += (Pressure_i + Pressure_j); }
    }
    
    /*--- iPoint inside the domain, jPoint on the boundary ---*/
    if (!boundary_i && boundary_j)
      if (geometry->node[iPoint]->GetDomain()) { iPoint_UndLapl[iPoint] += (Pressure_j - Pressure_i); jPoint_UndLapl[iPoint] += (Pressure_i + Pressure_j); }
    
    /*--- jPoint inside the domain, iPoint on the boundary ---*/
    if (boundary_i && !boundary_j)
      if (geometry->node[jPoint]->GetDomain()) { iPoint_UndLapl[jPoint] += (Pressure_i - Pressure_j); jPoint_UndLapl[jPoint] += (Pressure_i + Pressure_j); }
    
#endif
    
  }
  
#ifdef STRUCTURED_GRID
  unsigned short iMarker;
  unsigned long iVertex, Point_Normal;
  double Press_mirror;
  
  /*--- Loop over all boundaries and include an extra contribution
   from a mirror node. Find the nearest normal, interior point
   for a boundary node and make a linear approximation. ---*/
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
    if (config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE &&
        config->GetMarker_All_KindBC(iMarker) != INTERFACE_BOUNDARY &&
        config->GetMarker_All_KindBC(iMarker) != NEARFIELD_BOUNDARY &&
        config->GetMarker_All_KindBC(iMarker) != PERIODIC_BOUNDARY) {
      
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        if (geometry->node[iPoint]->GetDomain()) {
          
          Point_Normal = geometry->vertex[iMarker][iVertex]->GetNormal_Neighbor();
          
          Pressure_i = node[iPoint]->GetPressure();
          
          /*--- Interpolate & compute difference in the conserved variables ---*/
          Pressure_i = node[iPoint]->GetPressure();
          Press_mirror = 2.0*Pressure_i - node[Point_Normal]->GetPressure();
          
          /*--- Add contribution at the boundary node only ---*/
          iPoint_UndLapl[iPoint] += (Press_mirror - Pressure_i);
          jPoint_UndLapl[iPoint] += (Pressure_i + Press_mirror);
          
        }
      }
    }
  }
  
#endif
  
  /*--- Set pressure switch for each point ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++)
    node[iPoint]->SetSensor(fabs(iPoint_UndLapl[iPoint]) / jPoint_UndLapl[iPoint]);
  
}

void CEulerSolver::Inviscid_Forces(CGeometry *geometry, CConfig *config) {
  
  unsigned long iVertex, iPoint;
  unsigned short iDim, iMarker, Boundary, Monitoring, iMarker_Monitoring;
  double Pressure, *Normal = NULL, MomentDist[3], *Coord, Area,
  factor, NFPressOF, RefVel2, RefDensity, RefPressure, Gas_Constant, Mach2Vel, Mach_Motion, UnitNormal[3], Force[3];
  double *Origin = config->GetRefOriginMoment(0);
  string Marker_Tag, Monitoring_Tag;
  
  bool grid_movement      = config->GetGrid_Movement();
  double Alpha            = config->GetAoA()*PI_NUMBER/180.0;
  double Beta             = config->GetAoS()*PI_NUMBER/180.0;
  double RefAreaCoeff     = config->GetRefAreaCoeff();
  double RefLengthMoment  = config->GetRefLengthMoment();
  
  /*--- For dynamic meshes, use the motion Mach number as a reference value
   for computing the force coefficients. Otherwise, use the freestream values,
   which is the standard convention. ---*/
  if (grid_movement) {
    Gas_Constant = config->GetGas_ConstantND();
    Mach2Vel = sqrt(Gamma*Gas_Constant*config->GetTemperature_FreeStreamND());
    Mach_Motion = config->GetMach_Motion();
    RefVel2 = (Mach_Motion*Mach2Vel)*(Mach_Motion*Mach2Vel);
  }
  else {
    RefVel2 = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      RefVel2  += Velocity_Inf[iDim]*Velocity_Inf[iDim];
  }
  
  RefDensity  = Density_Inf;
  RefPressure = Pressure_Inf;
  
  factor = 1.0 / (0.5*RefDensity*RefAreaCoeff*RefVel2);
  
  /*-- Variables initialization ---*/
  
  Total_CDrag = 0.0;  Total_CLift = 0.0;    Total_CSideForce = 0.0;   Total_CEff = 0.0;
  Total_CMx = 0.0;    Total_CMy = 0.0;      Total_CMz = 0.0;
  Total_CFx = 0.0;    Total_CFy = 0.0;      Total_CFz = 0.0;
  Total_CT = 0.0;     Total_CQ = 0.0;       Total_CMerit = 0.0;
  Total_CNearFieldOF = 0.0;
  Total_Heat = 0.0;  Total_MaxHeat = 0.0;
  
  AllBound_CDrag_Inv = 0.0;   AllBound_CLift_Inv = 0.0; AllBound_CSideForce_Inv = 0.0;   AllBound_CEff_Inv = 0.0;
  AllBound_CMx_Inv = 0.0;     AllBound_CMy_Inv = 0.0;   AllBound_CMz_Inv = 0.0;
  AllBound_CFx_Inv = 0.0;     AllBound_CFy_Inv = 0.0;   AllBound_CFz_Inv = 0.0;
  AllBound_CT_Inv = 0.0;      AllBound_CQ_Inv = 0.0;    AllBound_CMerit_Inv = 0.0;
  AllBound_CNearFieldOF_Inv = 0.0;
  
  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    Surface_CLift_Inv[iMarker_Monitoring] = 0.0;
    Surface_CDrag_Inv[iMarker_Monitoring] = 0.0;
    Surface_CMx_Inv[iMarker_Monitoring]   = 0.0;
    Surface_CMy_Inv[iMarker_Monitoring]   = 0.0;
    Surface_CMz_Inv[iMarker_Monitoring]   = 0.0;
    Surface_CLift[iMarker_Monitoring]     = 0.0;
    Surface_CDrag[iMarker_Monitoring]     = 0.0;
    Surface_CMx[iMarker_Monitoring]       = 0.0;
    Surface_CMy[iMarker_Monitoring]       = 0.0;
    Surface_CMz[iMarker_Monitoring]       = 0.0;
  }
  
  /*--- Loop over the Euler and Navier-Stokes markers ---*/
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    Boundary   = config->GetMarker_All_KindBC(iMarker);
    Monitoring = config->GetMarker_All_Monitoring(iMarker);
    
    /*--- Obtain the origin for the moment computation for a particular marker ---*/
    if (Monitoring == YES) {
      for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
        Monitoring_Tag = config->GetMarker_Monitoring(iMarker_Monitoring);
        Marker_Tag = config->GetMarker_All_TagBound(iMarker);
        if (Marker_Tag == Monitoring_Tag)
          Origin = config->GetRefOriginMoment(iMarker_Monitoring);
      }
    }
    
    
    if ((Boundary == EULER_WALL) || (Boundary == HEAT_FLUX) ||
        (Boundary == ISOTHERMAL) || (Boundary == NEARFIELD_BOUNDARY)) {
      
      /*--- Forces initialization at each Marker ---*/
      
      CDrag_Inv[iMarker] = 0.0;         CLift_Inv[iMarker] = 0.0; CSideForce_Inv[iMarker] = 0.0;  CEff_Inv[iMarker] = 0.0;
      CMx_Inv[iMarker] = 0.0;           CMy_Inv[iMarker] = 0.0;   CMz_Inv[iMarker] = 0.0;
      CFx_Inv[iMarker] = 0.0;           CFy_Inv[iMarker] = 0.0;   CFz_Inv[iMarker] = 0.0;
      CT_Inv[iMarker] = 0.0;            CQ_Inv[iMarker] = 0.0;    CMerit_Inv[iMarker] = 0.0;
      CNearFieldOF_Inv[iMarker] = 0.0;
      
      for (iDim = 0; iDim < nDim; iDim++) ForceInviscid[iDim] = 0.0;
      MomentInviscid[0] = 0.0; MomentInviscid[1] = 0.0; MomentInviscid[2] = 0.0;
      NFPressOF = 0.0;
      
      /*--- Loop over the vertices to compute the forces ---*/
      
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
        
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        Pressure = node[iPoint]->GetPressure();
        
        CPressure[iMarker][iVertex] = (Pressure - RefPressure)*factor*RefAreaCoeff;
        
        /*--- Note that the pressure coefficient is computed at the
         halo cells (for visualization purposes), but not the forces ---*/
        
        if ( (geometry->node[iPoint]->GetDomain()) && (Monitoring == YES) ) {
          
          Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
          Coord = geometry->node[iPoint]->GetCoord();
          
          /*--- Quadratic objective function for the near-field.
           This uses the infinity pressure regardless of Mach number. ---*/
          
          NFPressOF += 0.5*(Pressure - Pressure_Inf)*(Pressure - Pressure_Inf)*Normal[nDim-1];
          
          Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);
          for (iDim = 0; iDim < nDim; iDim++) {
            UnitNormal[iDim] = Normal[iDim]/Area;
            MomentDist[iDim] = Coord[iDim] - Origin[iDim];
          }
          
          /*--- Force computation, note the minus sign due to the
           orientation of the normal (outward) ---*/
          
          for (iDim = 0; iDim < nDim; iDim++) {
            Force[iDim] = -(Pressure - Pressure_Inf)*Normal[iDim]*factor;
            ForceInviscid[iDim] += Force[iDim];
          }
          
          /*--- Moment with respect to the reference axis ---*/
          
          if (iDim == 3) {
            MomentInviscid[0] += (Force[2]*MomentDist[1]-Force[1]*MomentDist[2])/RefLengthMoment;
            MomentInviscid[1] += (Force[0]*MomentDist[2]-Force[2]*MomentDist[0])/RefLengthMoment;
          }
          MomentInviscid[2] += (Force[1]*MomentDist[0]-Force[0]*MomentDist[1])/RefLengthMoment;
          
        }
        
      }
      
      /*--- Project forces and store the non-dimensional coefficients ---*/
      
      if  (Monitoring == YES) {
        if (nDim == 2) {
          if (Boundary != NEARFIELD_BOUNDARY) {
            CDrag_Inv[iMarker]  =  ForceInviscid[0]*cos(Alpha) + ForceInviscid[1]*sin(Alpha);
            CLift_Inv[iMarker]  = -ForceInviscid[0]*sin(Alpha) + ForceInviscid[1]*cos(Alpha);
            CEff_Inv[iMarker]   = CLift_Inv[iMarker] / (CDrag_Inv[iMarker]+config->GetCteViscDrag()+EPS);
            CMz_Inv[iMarker]    = MomentInviscid[2];
            CFx_Inv[iMarker]    = ForceInviscid[0];
            CFy_Inv[iMarker]    = ForceInviscid[1];
            CT_Inv[iMarker]     = -CFx_Inv[iMarker];
            CQ_Inv[iMarker]     = -CMz_Inv[iMarker];
            CMerit_Inv[iMarker] = CT_Inv[iMarker]/CQ_Inv[iMarker];
          }
          else { CNearFieldOF_Inv[iMarker] = NFPressOF; }
        }
        if (nDim == 3) {
          if (Boundary != NEARFIELD_BOUNDARY) {
            CDrag_Inv[iMarker]      =  ForceInviscid[0]*cos(Alpha)*cos(Beta) + ForceInviscid[1]*sin(Beta) + ForceInviscid[2]*sin(Alpha)*cos(Beta);
            CLift_Inv[iMarker]      = -ForceInviscid[0]*sin(Alpha) + ForceInviscid[2]*cos(Alpha);
            CSideForce_Inv[iMarker] = -ForceInviscid[0]*sin(Beta)*cos(Alpha) + ForceInviscid[1]*cos(Beta) - ForceInviscid[2]*sin(Beta)*sin(Alpha);
            CEff_Inv[iMarker]       = CLift_Inv[iMarker] / (CDrag_Inv[iMarker]+config->GetCteViscDrag()+EPS);
            CMx_Inv[iMarker]        = MomentInviscid[0];
            CMy_Inv[iMarker]        = MomentInviscid[1];
            CMz_Inv[iMarker]        = MomentInviscid[2];
            CFx_Inv[iMarker]        = ForceInviscid[0];
            CFy_Inv[iMarker]        = ForceInviscid[1];
            CFz_Inv[iMarker]        = ForceInviscid[2];
            CT_Inv[iMarker]         = -CFz_Inv[iMarker];
            CQ_Inv[iMarker]         = -CMz_Inv[iMarker];
            CMerit_Inv[iMarker]     = CT_Inv[iMarker] / CQ_Inv[iMarker];
          }
          else { CNearFieldOF_Inv[iMarker] = NFPressOF; }
        }
        
        AllBound_CDrag_Inv        += CDrag_Inv[iMarker];
        AllBound_CLift_Inv        += CLift_Inv[iMarker];
        AllBound_CSideForce_Inv   += CSideForce_Inv[iMarker];
        AllBound_CMx_Inv          += CMx_Inv[iMarker];
        AllBound_CMy_Inv          += CMy_Inv[iMarker];
        AllBound_CMz_Inv          += CMz_Inv[iMarker];
        AllBound_CFx_Inv          += CFx_Inv[iMarker];
        AllBound_CFy_Inv          += CFy_Inv[iMarker];
        AllBound_CFz_Inv          += CFz_Inv[iMarker];
        AllBound_CT_Inv           += CT_Inv[iMarker];
        AllBound_CQ_Inv           += CQ_Inv[iMarker];
        AllBound_CNearFieldOF_Inv += CNearFieldOF_Inv[iMarker];
        
        AllBound_CEff_Inv = AllBound_CLift_Inv / (AllBound_CDrag_Inv + config->GetCteViscDrag() + EPS);
        AllBound_CMerit_Inv = AllBound_CT_Inv / (AllBound_CQ_Inv + EPS);
        
        /*--- Compute the coefficients per surface ---*/
        
        for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
          Monitoring_Tag = config->GetMarker_Monitoring(iMarker_Monitoring);
          Marker_Tag = config->GetMarker_All_TagBound(iMarker);
          if (Marker_Tag == Monitoring_Tag) {
            Surface_CLift_Inv[iMarker_Monitoring] += CLift_Inv[iMarker];
            Surface_CDrag_Inv[iMarker_Monitoring] += CDrag_Inv[iMarker];
            Surface_CMx_Inv[iMarker_Monitoring]   += CMx_Inv[iMarker];
            Surface_CMy_Inv[iMarker_Monitoring]   += CMy_Inv[iMarker];
            Surface_CMz_Inv[iMarker_Monitoring]   += CMz_Inv[iMarker];
          }
        }
        
      }
      
      
    }
  }
  
  /*--- Update the total coefficients (note that all the nodes have the same value) ---*/
  
  Total_CDrag         = AllBound_CDrag_Inv;
  Total_CLift         = AllBound_CLift_Inv;
  Total_CSideForce    = AllBound_CSideForce_Inv;
  Total_CEff          = Total_CLift / (Total_CDrag + config->GetCteViscDrag() + EPS);
  Total_CMx           = AllBound_CMx_Inv;
  Total_CMy           = AllBound_CMy_Inv;
  Total_CMz           = AllBound_CMz_Inv;
  Total_CFx           = AllBound_CFx_Inv;
  Total_CFy           = AllBound_CFy_Inv;
  Total_CFz           = AllBound_CFz_Inv;
  Total_CT            = AllBound_CT_Inv;
  Total_CQ            = AllBound_CQ_Inv;
  Total_CMerit        = Total_CT / (Total_CQ + EPS);
  Total_CNearFieldOF  = AllBound_CNearFieldOF_Inv;
  
  /*--- Update the total coefficients per surface (note that all the nodes have the same value)---*/
  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    Surface_CLift[iMarker_Monitoring]     = Surface_CLift_Inv[iMarker_Monitoring];
    Surface_CDrag[iMarker_Monitoring]     = Surface_CDrag_Inv[iMarker_Monitoring];
    Surface_CMx[iMarker_Monitoring]       = Surface_CMx_Inv[iMarker_Monitoring];
    Surface_CMy[iMarker_Monitoring]       = Surface_CMy_Inv[iMarker_Monitoring];
    Surface_CMz[iMarker_Monitoring]       = Surface_CMz_Inv[iMarker_Monitoring];
  }
  
}

void CEulerSolver::ExplicitRK_Iteration(CGeometry *geometry, CSolver **solver_container,
                                        CConfig *config, unsigned short iRKStep) {
  double *Residual, *Res_TruncError, Vol, Delta, Res;
  unsigned short iVar;
  unsigned long iPoint;
  
  double RK_AlphaCoeff = config->Get_Alpha_RKStep(iRKStep);
  
  for (iVar = 0; iVar < nVar; iVar++) {
    SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }
  
  /*--- Update the solution ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    Vol = geometry->node[iPoint]->GetVolume();
    Delta = node[iPoint]->GetDelta_Time() / Vol;
    
    Res_TruncError = node[iPoint]->GetResTruncError();
    Residual = LinSysRes.GetBlock(iPoint);
    
    for (iVar = 0; iVar < nVar; iVar++) {
      Res = Residual[iVar] + Res_TruncError[iVar];
      node[iPoint]->AddSolution(iVar, -Res*Delta*RK_AlphaCoeff);
      AddRes_RMS(iVar, Res*Res);
      AddRes_Max(iVar, fabs(Res), geometry->node[iPoint]->GetGlobalIndex());
    }
    
  }
  
  /*--- Compute the root mean square residual ---*/
  SetResidual_RMS(geometry, config);
  
  
}

void CEulerSolver::ExplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
  double *local_Residual, *local_Res_TruncError, Vol, Delta, Res;
  unsigned short iVar;
  unsigned long iPoint;
  
  for (iVar = 0; iVar < nVar; iVar++) {
    SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }
  
  /*--- Update the solution ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    Vol = geometry->node[iPoint]->GetVolume();
    Delta = node[iPoint]->GetDelta_Time() / Vol;
    
    local_Res_TruncError = node[iPoint]->GetResTruncError();
    local_Residual = LinSysRes.GetBlock(iPoint);
    
    for (iVar = 0; iVar < nVar; iVar++) {
      Res = local_Residual[iVar] + local_Res_TruncError[iVar];
      node[iPoint]->AddSolution(iVar, -Res*Delta);
      AddRes_RMS(iVar, Res*Res);
      AddRes_Max(iVar, fabs(Res), geometry->node[iPoint]->GetGlobalIndex());
    }
    
  }
  
  /*--- Compute the root mean square residual ---*/
  
  SetResidual_RMS(geometry, config);
  
}

void CEulerSolver::ImplicitEuler_Iteration(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned short iMesh) {
  
  unsigned short iVar, jVar, iMGlevel;
  unsigned long iPoint, total_index, IterLinSol = 0;
  double Delta, *local_Res_TruncError, Vol;
  
  bool roe_turkel = (config->GetKind_Upwind_Flow() == TURKEL);

  /*--- Set maximum residual to zero ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }
  
  /*--- Build implicit system ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    
    /*--- Read the residual ---*/
    
    local_Res_TruncError = node[iPoint]->GetResTruncError();
    
    /*--- Read the volume ---*/
    
    Vol = geometry[iMesh]->node[iPoint]->GetVolume();
    
    /*--- Modify matrix diagonal to assure diagonal dominance ---*/
    
    if (node[iPoint]->GetDelta_Time() != 0.0) {
      Delta = Vol / node[iPoint]->GetDelta_Time();
      if (roe_turkel) {
        SetPreconditioner(config, iPoint);
        for (iVar = 0; iVar < nVar; iVar ++ )
          for (jVar = 0; jVar < nVar; jVar ++ )
            LowMach_Precontioner[iVar][jVar] = Delta*LowMach_Precontioner[iVar][jVar];
        Jacobian.AddBlock(iPoint, iPoint, LowMach_Precontioner);
      }
      else {
        Jacobian.AddVal2Diag(iPoint, Delta);
      }
    }
    else {
      Jacobian.SetVal2Diag(iPoint, 1.0);
      for (iVar = 0; iVar < nVar; iVar++) {
        total_index = iPoint*nVar + iVar;
        LinSysRes[total_index] = 0.0;
        local_Res_TruncError[iVar] = 0.0;
      }
    }
    
    /*--- Right hand side of the system (-Residual) and initial guess (x = 0) ---*/
    
    for (iVar = 0; iVar < nVar; iVar++) {
      total_index = iPoint*nVar + iVar;
      LinSysRes[total_index] = - (LinSysRes[total_index] + local_Res_TruncError[iVar]);
      LinSysSol[total_index] = 0.0;
      AddRes_RMS(iVar, LinSysRes[total_index]*LinSysRes[total_index]);
      AddRes_Max(iVar, fabs(LinSysRes[total_index]), geometry[iMesh]->node[iPoint]->GetGlobalIndex());
    }
  }
  
  /*--- Initialize residual and solution at the ghost points ---*/
  
  for (iPoint = nPointDomain; iPoint < nPoint; iPoint++) {
    for (iVar = 0; iVar < nVar; iVar++) {
      total_index = iPoint*nVar + iVar;
      LinSysRes[total_index] = 0.0;
      LinSysSol[total_index] = 0.0;
    }
  }

  /*--- Solve or smooth the linear system ---*/
  
  CSysMatrix** Jacobian_Array = new CSysMatrix*[config->GetMGLevels()+1];
  CSysVector** LinSysRes_Array = new CSysVector*[config->GetMGLevels()+1];
  CSysVector** LinSysSol_Array = new CSysVector*[config->GetMGLevels()+1];

  for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
    Jacobian_Array[iMGlevel] = &solver_container[iMGlevel][FLOW_SOL]->Jacobian;
    LinSysRes_Array[iMGlevel] = &solver_container[iMGlevel][FLOW_SOL]->LinSysRes;
    LinSysSol_Array[iMGlevel] = &solver_container[iMGlevel][FLOW_SOL]->LinSysSol;
  }
  
  CSysSolve system;
  
  IterLinSol = system.Solve(Jacobian_Array, LinSysRes_Array, LinSysSol_Array, geometry, config, iMesh);
  
  /*--- The the number of iterations of the linear solver ---*/
  
  SetIterLinSolver(IterLinSol);
  
  /*--- Update solution (system written in terms of increments) ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    for (iVar = 0; iVar < nVar; iVar++) {
      node[iPoint]->AddSolution(iVar, config->GetLinear_Solver_Relax()*LinSysSol[iPoint*nVar+iVar]);
    }
  }
  
  /*--- Compute the root mean square residual ---*/
  
  SetResidual_RMS(geometry[iMesh], config);
  
}

void CEulerSolver::SetPrimVar_Gradient_GG(CGeometry *geometry, CConfig *config) {
  unsigned long iPoint, jPoint, iEdge, iVertex;
  unsigned short iDim, iVar, iMarker;
  double *PrimVar_Vertex, *PrimVar_i, *PrimVar_j, PrimVar_Average,
  Partial_Gradient, Partial_Res, *Normal;
  
  /*--- Gradient primitive variables compressible (temp, vx, vy, vz, P, rho)
   Gradient primitive variables incompressible (rho, vx, vy, vz, beta) ---*/
  PrimVar_Vertex = new double [nPrimVarGrad];
  PrimVar_i = new double [nPrimVarGrad];
  PrimVar_j = new double [nPrimVarGrad];
  
  /*--- Set Gradient_Primitive to zero ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++)
    node[iPoint]->SetGradient_PrimitiveZero(nPrimVarGrad);
  
  /*--- Loop interior edges ---*/
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    
    for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
      PrimVar_i[iVar] = node[iPoint]->GetPrimVar(iVar);
      PrimVar_j[iVar] = node[jPoint]->GetPrimVar(iVar);
    }
    
    Normal = geometry->edge[iEdge]->GetNormal();
    for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
      PrimVar_Average =  0.5 * ( PrimVar_i[iVar] + PrimVar_j[iVar] );
      for (iDim = 0; iDim < nDim; iDim++) {
        Partial_Res = PrimVar_Average*Normal[iDim];
        if (geometry->node[iPoint]->GetDomain())
          node[iPoint]->AddGradient_Primitive(iVar, iDim, Partial_Res);
        if (geometry->node[jPoint]->GetDomain())
          node[jPoint]->SubtractGradient_Primitive(iVar, iDim, Partial_Res);
      }
    }
  }
  
  /*--- Loop boundary edges ---*/
  for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
    for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
      iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
      if (geometry->node[iPoint]->GetDomain()) {
        
        for (iVar = 0; iVar < nPrimVarGrad; iVar++)
          PrimVar_Vertex[iVar] = node[iPoint]->GetPrimVar(iVar);
        
        Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
        for (iVar = 0; iVar < nPrimVarGrad; iVar++)
          for (iDim = 0; iDim < nDim; iDim++) {
            Partial_Res = PrimVar_Vertex[iVar]*Normal[iDim];
            node[iPoint]->SubtractGradient_Primitive(iVar, iDim, Partial_Res);
          }
      }
    }
  }
  
  /*--- Update gradient value ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        Partial_Gradient = node[iPoint]->GetGradient_Primitive(iVar,iDim) / (geometry->node[iPoint]->GetVolume());
        node[iPoint]->SetGradient_Primitive(iVar, iDim, Partial_Gradient);
      }
    }
  }
  
  delete [] PrimVar_Vertex;
  delete [] PrimVar_i;
  delete [] PrimVar_j;
  
}

void CEulerSolver::SetPrimVar_Gradient_LS(CGeometry *geometry, CConfig *config) {
  
  unsigned short iVar, iDim, jDim, iNeigh;
  unsigned long iPoint, jPoint;
  double *PrimVar_i, *PrimVar_j, *Coord_i, *Coord_j, r11, r12, r13, r22, r23, r23_a,
  r23_b, r33, weight, product, z11, z12, z13, z22, z23, z33, detR2;
  bool singular;
  
  /*--- Loop over points of the grid ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    
    /*--- Set the value of the singular ---*/
    singular = false;
    
    /*--- Get coordinates ---*/
    
    Coord_i = geometry->node[iPoint]->GetCoord();
    
    /*--- Get primitives from CVariable ---*/
    
    PrimVar_i = node[iPoint]->GetPrimVar();
    
    /*--- Inizialization of variables ---*/
    
    for (iVar = 0; iVar < nPrimVarGrad; iVar++)
      for (iDim = 0; iDim < nDim; iDim++)
        cvector[iVar][iDim] = 0.0;
    
    r11 = 0.0; r12 = 0.0;   r13 = 0.0;    r22 = 0.0;
    r23 = 0.0; r23_a = 0.0; r23_b = 0.0;  r33 = 0.0; detR2 = 0.0;
    
    for (iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); iNeigh++) {
      jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
      Coord_j = geometry->node[jPoint]->GetCoord();
      
      PrimVar_j = node[jPoint]->GetPrimVar();
      
      weight = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        weight += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
      
      /*--- Sumations for entries of upper triangular matrix R ---*/
      
      if (weight != 0.0) {
        
        r11 += (Coord_j[0]-Coord_i[0])*(Coord_j[0]-Coord_i[0])/weight;
        r12 += (Coord_j[0]-Coord_i[0])*(Coord_j[1]-Coord_i[1])/weight;
        r22 += (Coord_j[1]-Coord_i[1])*(Coord_j[1]-Coord_i[1])/weight;
        
        if (nDim == 3) {
          r13 += (Coord_j[0]-Coord_i[0])*(Coord_j[2]-Coord_i[2])/weight;
          r23_a += (Coord_j[1]-Coord_i[1])*(Coord_j[2]-Coord_i[2])/weight;
          r23_b += (Coord_j[0]-Coord_i[0])*(Coord_j[2]-Coord_i[2])/weight;
          r33 += (Coord_j[2]-Coord_i[2])*(Coord_j[2]-Coord_i[2])/weight;
        }
        
        /*--- Entries of c:= transpose(A)*b ---*/
        
        for (iVar = 0; iVar < nPrimVarGrad; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            cvector[iVar][iDim] += (Coord_j[iDim]-Coord_i[iDim])*(PrimVar_j[iVar]-PrimVar_i[iVar])/weight;
        
      }
      
    }
    
    /*--- Entries of upper triangular matrix R ---*/
    
    if (r11 >= 0.0) r11 = sqrt(r11); else r11 = 0.0;
    if (r11 != 0.0) r12 = r12/r11; else r12 = 0.0;
    if (r22-r12*r12 >= 0.0) r22 = sqrt(r22-r12*r12); else r22 = 0.0;
    
    if (nDim == 3) {
      if (r11 != 0.0) r13 = r13/r11; else r13 = 0.0;
      if ((r22 != 0.0) && (r11*r22 != 0.0)) r23 = r23_a/r22 - r23_b*r12/(r11*r22); else r23 = 0.0;
      if (r33-r23*r23-r13*r13 >= 0.0) r33 = sqrt(r33-r23*r23-r13*r13); else r33 = 0.0;
    }
    
    /*--- Compute determinant ---*/
    
    if (nDim == 2) detR2 = (r11*r22)*(r11*r22);
    else detR2 = (r11*r22*r33)*(r11*r22*r33);
    
    /*--- Detect singular matrices ---*/
    
    if (abs(detR2) <= EPS) { detR2 = 1.0; singular = true; }
    
    /*--- S matrix := inv(R)*traspose(inv(R)) ---*/
    
    if (singular) {
      for (iDim = 0; iDim < nDim; iDim++)
        for (jDim = 0; jDim < nDim; jDim++)
          Smatrix[iDim][jDim] = 0.0;
    }
    else {
      if (nDim == 2) {
        Smatrix[0][0] = (r12*r12+r22*r22)/detR2;
        Smatrix[0][1] = -r11*r12/detR2;
        Smatrix[1][0] = Smatrix[0][1];
        Smatrix[1][1] = r11*r11/detR2;
      }
      else {
        z11 = r22*r33; z12 = -r12*r33; z13 = r12*r23-r13*r22;
        z22 = r11*r33; z23 = -r11*r23; z33 = r11*r22;
        Smatrix[0][0] = (z11*z11+z12*z12+z13*z13)/detR2;
        Smatrix[0][1] = (z12*z22+z13*z23)/detR2;
        Smatrix[0][2] = (z13*z33)/detR2;
        Smatrix[1][0] = Smatrix[0][1];
        Smatrix[1][1] = (z22*z22+z23*z23)/detR2;
        Smatrix[1][2] = (z23*z33)/detR2;
        Smatrix[2][0] = Smatrix[0][2];
        Smatrix[2][1] = Smatrix[1][2];
        Smatrix[2][2] = (z33*z33)/detR2;
      }
    }
    
    /*--- Computation of the gradient: S*c ---*/
    for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        product = 0.0;
        for (jDim = 0; jDim < nDim; jDim++) {
          product += Smatrix[iDim][jDim]*cvector[iVar][jDim];
        }
        
        node[iPoint]->SetGradient_Primitive(iVar, iDim, product);
      }
    }
    
  }
  
}

void CEulerSolver::SetPrimVar_Limiter(CGeometry *geometry, CConfig *config) {
  
  unsigned long iEdge, iPoint, jPoint;
  unsigned short iVar, iDim;
  double **Gradient_i, **Gradient_j, *Coord_i, *Coord_j, *Primitive_i, *Primitive_j,
  dave, LimK, eps1, eps2, dm, dp, du, limiter;
  
  /*--- Initialize solution max and solution min in the entire domain --*/
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
      node[iPoint]->SetSolution_Max(iVar, -EPS);
      node[iPoint]->SetSolution_Min(iVar, EPS);
      node[iPoint]->SetLimiter_Primitive(iVar, 2.0);
    }
  }
  
  /*--- Establish bounds for Spekreijse monotonicity by finding max & min values of neighbor variables --*/
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    /*--- Point identification, Normal vector and area ---*/
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    
    /*--- Get the conserved variables ---*/
    Primitive_i = node[iPoint]->GetPrimVar();
    Primitive_j = node[jPoint]->GetPrimVar();
    
    /*--- Compute the maximum, and minimum values for nodes i & j ---*/
    for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
      du = (Primitive_j[iVar] - Primitive_i[iVar]);
      node[iPoint]->SetSolution_Min(iVar, min(node[iPoint]->GetSolution_Min(iVar), du));
      node[iPoint]->SetSolution_Max(iVar, max(node[iPoint]->GetSolution_Max(iVar), du));
      node[jPoint]->SetSolution_Min(iVar, min(node[jPoint]->GetSolution_Min(iVar), -du));
      node[jPoint]->SetSolution_Max(iVar, max(node[jPoint]->GetSolution_Max(iVar), -du));
    }
  }
  
  /*--- Venkatakrishnan (Venkatakrishnan 1994) limiter ---*/
  
  /*-- Get limiter parameters from the configuration file ---*/
  dave = config->GetRefElemLength();
  LimK = config->GetLimiterCoeff();
  eps1 = LimK*dave;
  eps2 = eps1*eps1*eps1;
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    iPoint     = geometry->edge[iEdge]->GetNode(0);
    jPoint     = geometry->edge[iEdge]->GetNode(1);
    Primitive_i = node[iPoint]->GetPrimVar();
    Primitive_j = node[jPoint]->GetPrimVar();
    Gradient_i = node[iPoint]->GetGradient_Primitive();
    Gradient_j = node[jPoint]->GetGradient_Primitive();
    Coord_i    = geometry->node[iPoint]->GetCoord();
    Coord_j    = geometry->node[jPoint]->GetCoord();
    
    for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
      
      /*--- Calculate the interface left gradient, delta- (dm) ---*/
      
      dm = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        dm += 0.5*(Coord_j[iDim]-Coord_i[iDim])*Gradient_i[iVar][iDim];
      
      /*--- Calculate the interface right gradient, delta+ (dp) ---*/
      
      if ( dm > 0.0 ) dp = node[iPoint]->GetSolution_Max(iVar);
      else dp = node[iPoint]->GetSolution_Min(iVar);
      
      limiter = ( dp*dp + 2.0*dp*dm + eps2 )/( dp*dp + dp*dm + 2.0*dm*dm + eps2);
      
      if (limiter < node[iPoint]->GetLimiter_Primitive(iVar))
        node[iPoint]->SetLimiter_Primitive(iVar, limiter);
      
      /*-- Repeat for point j on the edge ---*/
      
      dm = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        dm += 0.5*(Coord_i[iDim]-Coord_j[iDim])*Gradient_j[iVar][iDim];
      
      if ( dm > 0.0 ) dp = node[jPoint]->GetSolution_Max(iVar);
      else dp = node[jPoint]->GetSolution_Min(iVar);
      
      limiter = ( dp*dp + 2.0*dp*dm + eps2 )/( dp*dp + dp*dm + 2.0*dm*dm + eps2);
      
      if (limiter < node[jPoint]->GetLimiter_Primitive(iVar))
        node[jPoint]->SetLimiter_Primitive(iVar, limiter);
      
    }
  }
  
}

void CEulerSolver::SetPreconditioner(CConfig *config, unsigned short iPoint) {
  unsigned short iDim, jDim, iVar, jVar;
  double Beta, local_Mach, Beta2, rho, enthalpy, soundspeed, sq_vel;
  double *U_i = NULL;
  double Beta_min = config->GetminTurkelBeta();
  double Beta_max = config->GetmaxTurkelBeta();
  
  
  /*--- Variables to calculate the preconditioner parameter Beta ---*/
  local_Mach = sqrt(node[iPoint]->GetVelocity2())/node[iPoint]->GetSoundSpeed();
  Beta 		    = max(Beta_min,min(local_Mach,Beta_max));
  Beta2 		    = Beta*Beta;
  
  U_i = node[iPoint]->GetSolution();
  
  rho = U_i[0];
  enthalpy = node[iPoint]->GetEnthalpy();
  soundspeed = node[iPoint]->GetSoundSpeed();
  sq_vel = node[iPoint]->GetVelocity2();
  
  /*---Calculating the inverse of the preconditioning matrix that multiplies the time derivative  */
  LowMach_Precontioner[0][0] = 0.5*sq_vel;
  LowMach_Precontioner[0][nVar-1] = 1.0;
  for (iDim = 0; iDim < nDim; iDim ++)
    LowMach_Precontioner[0][1+iDim] = -1.0*U_i[iDim+1]/rho;
  
  for (iDim = 0; iDim < nDim; iDim ++) {
    LowMach_Precontioner[iDim+1][0] = 0.5*sq_vel*U_i[iDim+1]/rho;
    LowMach_Precontioner[iDim+1][nVar-1] = U_i[iDim+1]/rho;
    for (jDim = 0; jDim < nDim; jDim ++) {
      LowMach_Precontioner[iDim+1][1+jDim] = -1.0*U_i[jDim+1]/rho*U_i[iDim+1]/rho;
    }
  }
  
  LowMach_Precontioner[nVar-1][0] = 0.5*sq_vel*enthalpy;
  LowMach_Precontioner[nVar-1][nVar-1] = enthalpy;
  for (iDim = 0; iDim < nDim; iDim ++)
    LowMach_Precontioner[nVar-1][1+iDim] = -1.0*U_i[iDim+1]/rho*enthalpy;
  
  
  for (iVar = 0; iVar < nVar; iVar ++ ) {
    for (jVar = 0; jVar < nVar; jVar ++ ) {
      LowMach_Precontioner[iVar][jVar] = (1.0/(Beta2+EPS) - 1.0) * (Gamma-1.0)/(soundspeed*soundspeed)*LowMach_Precontioner[iVar][jVar];
      if (iVar == jVar)
        LowMach_Precontioner[iVar][iVar] += 1.0;
    }
  }
  
}

void CEulerSolver::BC_Euler_Wall(CGeometry *geometry, CSolver **solver_container,
                                 CNumerics *numerics, CConfig *config, unsigned short val_marker) {
  
  unsigned short iDim, iVar, jVar, jDim;
  unsigned long iPoint, iVertex;
  double Pressure, *Normal = NULL, *GridVel = NULL, Area, UnitNormal[3],
  ProjGridVel = 0.0, a2, phi, turb_ke = 0.0;
  
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool grid_movement  = config->GetGrid_Movement();
  bool tkeNeeded = ((config->GetKind_Solver() == RANS) && (config->GetKind_Turb_Model() == SST));
  
  /*--- Loop over all the vertices on this boundary marker ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      
      Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
      
      Area = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
      Area = sqrt (Area);
      
      for (iDim = 0; iDim < nDim; iDim++) UnitNormal[iDim] = -Normal[iDim]/Area;
      
      /*--- Get the pressure ---*/
      
      Pressure = node[iPoint]->GetPressure();
      
      /*--- Add the kinetic energy correction ---*/
      
      if (tkeNeeded) {
        turb_ke = solver_container[TURB_SOL]->node[iPoint]->GetSolution(0);
        Pressure += (2.0/3.0)*node[iPoint]->GetDensity()*turb_ke;
      }
      
      /*--- Compute the residual ---*/
      
      Residual[0] = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        Residual[iDim+1] = Pressure*UnitNormal[iDim]*Area;
      
      Residual[nVar-1] = 0.0;
      
      /*--- Adjustment to energy equation due to grid motion ---*/
      
      if (grid_movement) {
        ProjGridVel = 0.0;
        GridVel = geometry->node[iPoint]->GetGridVel();
        for (iDim = 0; iDim < nDim; iDim++)
          ProjGridVel += GridVel[iDim]*UnitNormal[iDim]*Area;
        Residual[nVar-1] = Pressure*ProjGridVel;
      }
      
      /*--- Add value to the residual ---*/
      
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Form Jacobians for implicit computations ---*/
      
      if (implicit) {
        
        /*--- Initialize jacobian ---*/
        
        for (iVar = 0; iVar < nVar; iVar++) {
          for (jVar = 0; jVar < nVar; jVar++)
            Jacobian_i[iVar][jVar] = 0.0;
        }
        
        a2 = Gamma-1.0;
        phi = 0.5*a2*node[iPoint]->GetVelocity2();
        for (iVar = 0; iVar < nVar; iVar++) {
          Jacobian_i[0][iVar] = 0.0;
          Jacobian_i[nDim+1][iVar] = 0.0;
        }
        for (iDim = 0; iDim < nDim; iDim++) {
          Jacobian_i[iDim+1][0] = -phi*Normal[iDim];
          for (jDim = 0; jDim < nDim; jDim++)
            Jacobian_i[iDim+1][jDim+1] = a2*node[iPoint]->GetVelocity(jDim)*Normal[iDim];
          Jacobian_i[iDim+1][nDim+1] = -a2*Normal[iDim];
        }
        if (grid_movement) {
          ProjGridVel = 0.0;
          GridVel = geometry->node[iPoint]->GetGridVel();
          for (iDim = 0; iDim < nDim; iDim++)
            ProjGridVel += GridVel[iDim]*UnitNormal[iDim]*Area;
          Jacobian_i[nDim+1][0] = phi*ProjGridVel;
          for (jDim = 0; jDim < nDim; jDim++)
            Jacobian_i[nDim+1][jDim+1] = -a2*node[iPoint]->GetVelocity(jDim)*ProjGridVel;
          Jacobian_i[nDim+1][nDim+1] = a2*ProjGridVel;
        }
        Jacobian.AddBlock(iPoint,iPoint,Jacobian_i);
        
      }
    }
  }
  
}

void CEulerSolver::BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  unsigned short iDim;
  unsigned long iVertex, iPoint, Point_Normal;
  
  double *GridVel;
  double Area, UnitNormal[3];
  double Density, Pressure, Velocity[3], Energy;
  double Density_Bound, Pressure_Bound, Vel_Bound[3];
  double Density_Infty, Pressure_Infty, Vel_Infty[3];
  double SoundSpeed, Entropy, Velocity2, Vn;
  double SoundSpeed_Bound, Entropy_Bound, Vel2_Bound, Vn_Bound;
  double SoundSpeed_Infty, Entropy_Infty, Vel2_Infty, Vn_Infty, Qn_Infty;
  double RiemannPlus, RiemannMinus;
  double *V_infty, *V_domain;
  
  double Gas_Constant     = config->GetGas_ConstantND();
  
  bool implicit         = config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT;
  bool grid_movement    = config->GetGrid_Movement();
  bool viscous          = config->GetViscous();
  bool tkeNeeded = ((config->GetKind_Solver() == RANS) && (config->GetKind_Turb_Model() == SST));
  
  double *Normal = new double[nDim];
  
  /*--- Loop over all the vertices on this boundary marker ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Allocate the value at the infinity ---*/
    V_infty = GetCharacPrimVar(val_marker, iVertex);
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Index of the closest interior node ---*/
      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
      
      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);
      
      /*--- Retrieve solution at the farfield boundary node ---*/
      V_domain = node[iPoint]->GetPrimVar();
      
      /*--- Construct solution state at infinity (far-field) ---*/
      
      /*--- Construct solution state at infinity for compressible flow by
       using Riemann invariants, and then impose a weak boundary condition
       by computing the flux using this new state for U. See CFD texts by
       Hirsch or Blazek for more detail. Adapted from an original
       implementation in the Stanford University multi-block (SUmb) solver
       in the routine bcFarfield.f90 written by Edwin van der Weide,
       last modified 06-12-2005. First, compute the unit normal at the
       boundary nodes. ---*/
      
      Area = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
      Area = sqrt(Area);
      
      for (iDim = 0; iDim < nDim; iDim++)
        UnitNormal[iDim] = Normal[iDim]/Area;
      
      /*--- Store primitive variables (density, velocities, velocity squared,
       energy, pressure, and sound speed) at the boundary node, and set some
       other quantities for clarity. Project the current flow velocity vector
       at this boundary node into the local normal direction, i.e. compute
       v_bound.n.  ---*/
      
      Density_Bound = V_domain[nDim+2];
      Vel2_Bound = 0.0; Vn_Bound = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        Vel_Bound[iDim] = V_domain[iDim+1];
        Vel2_Bound     += Vel_Bound[iDim]*Vel_Bound[iDim];
        Vn_Bound       += Vel_Bound[iDim]*UnitNormal[iDim];
      }
      Pressure_Bound   = node[iPoint]->GetPressure();
      SoundSpeed_Bound = sqrt(Gamma*Pressure_Bound/Density_Bound);
      Entropy_Bound    = pow(Density_Bound,Gamma)/Pressure_Bound;
      
      /*--- Store the primitive variable state for the freestream. Project
       the freestream velocity vector into the local normal direction,
       i.e. compute v_infty.n. ---*/
      
      Density_Infty = GetDensity_Inf();
      Vel2_Infty = 0.0; Vn_Infty = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        Vel_Infty[iDim] = GetVelocity_Inf(iDim);
        Vel2_Infty     += Vel_Infty[iDim]*Vel_Infty[iDim];
        Vn_Infty       += Vel_Infty[iDim]*UnitNormal[iDim];
      }
      Pressure_Infty   = GetPressure_Inf();
      SoundSpeed_Infty = sqrt(Gamma*Pressure_Infty/Density_Infty);
      Entropy_Infty    = pow(Density_Infty,Gamma)/Pressure_Infty;
      
      /*--- Adjust the normal freestream velocity for grid movement ---*/
      
      Qn_Infty = Vn_Infty;
      if (grid_movement) {
        GridVel = geometry->node[iPoint]->GetGridVel();
        for (iDim = 0; iDim < nDim; iDim++)
          Qn_Infty -= GridVel[iDim]*UnitNormal[iDim];
      }
      
      /*--- Compute acoustic Riemann invariants: R = u.n +/- 2c/(gamma-1).
       These correspond with the eigenvalues (u+c) and (u-c), respectively,
       which represent the acoustic waves. Positive characteristics are
       incoming, and a physical boundary condition is imposed (freestream
       state). This occurs when either (u.n+c) > 0 or (u.n-c) > 0. Negative
       characteristics are leaving the domain, and numerical boundary
       conditions are required by extrapolating from the interior state
       using the Riemann invariants. This occurs when (u.n+c) < 0 or
       (u.n-c) < 0. Note that grid movement is taken into account when
       checking the sign of the eigenvalue. ---*/
      
      /*--- Check whether (u.n+c) is greater or less than zero ---*/
      
      if (Qn_Infty > -SoundSpeed_Infty) {
        /*--- Subsonic inflow or outflow ---*/
        RiemannPlus = Vn_Bound + 2.0*SoundSpeed_Bound/Gamma_Minus_One;
      } else {
        /*--- Supersonic inflow ---*/
        RiemannPlus = Vn_Infty + 2.0*SoundSpeed_Infty/Gamma_Minus_One;
      }
      
      /*--- Check whether (u.n-c) is greater or less than zero ---*/
      
      if (Qn_Infty > SoundSpeed_Infty) {
        /*--- Supersonic outflow ---*/
        RiemannMinus = Vn_Bound - 2.0*SoundSpeed_Bound/Gamma_Minus_One;
      } else {
        /*--- Subsonic outflow ---*/
        RiemannMinus = Vn_Infty - 2.0*SoundSpeed_Infty/Gamma_Minus_One;
      }
      
      /*--- Compute a new value for the local normal velocity and speed of
       sound from the Riemann invariants. ---*/
      
      Vn = 0.5 * (RiemannPlus + RiemannMinus);
      SoundSpeed = 0.25 * (RiemannPlus - RiemannMinus)*Gamma_Minus_One;
      
      /*--- Construct the primitive variable state at the boundary for
       computing the flux for the weak boundary condition. The values
       that we choose to construct the solution (boundary or freestream)
       depend on whether we are at an inflow or outflow. At an outflow, we
       choose boundary information (at most one characteristic is incoming),
       while at an inflow, we choose infinity values (at most one
       characteristic is outgoing). ---*/
      
      if (Qn_Infty > 0.0)   {
        /*--- Outflow conditions ---*/
        for (iDim = 0; iDim < nDim; iDim++)
          Velocity[iDim] = Vel_Bound[iDim] + (Vn-Vn_Bound)*UnitNormal[iDim];
        Entropy = Entropy_Bound;
      } else  {
        /*--- Inflow conditions ---*/
        for (iDim = 0; iDim < nDim; iDim++)
          Velocity[iDim] = Vel_Infty[iDim] + (Vn-Vn_Infty)*UnitNormal[iDim];
        Entropy = Entropy_Infty;
      }
      
      /*--- Recompute the primitive variables. ---*/
      
      Density = pow(Entropy*SoundSpeed*SoundSpeed/Gamma,1.0/Gamma_Minus_One);
      Velocity2 = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        Velocity2 += Velocity[iDim]*Velocity[iDim];
      }
      Pressure = Density*SoundSpeed*SoundSpeed/Gamma;
      Energy   = Pressure/(Gamma_Minus_One*Density) + 0.5*Velocity2;
      if (tkeNeeded) Energy += GetTke_Inf();
      
      /*--- Store new primitive state for computing the flux. ---*/
      
      V_infty[0] = Pressure/(Gas_Constant*Density);
      for (iDim = 0; iDim < nDim; iDim++)
        V_infty[iDim+1] = Velocity[iDim];
      V_infty[nDim+1] = Pressure;
      V_infty[nDim+2] = Density;
      V_infty[nDim+3] = Energy + Pressure/Density;
      
      
      /*--- Set various quantities in the numerics class ---*/
      conv_numerics->SetPrimitive(V_domain, V_infty);
      
      if (grid_movement) {
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
                                  geometry->node[iPoint]->GetGridVel());
      }
      
      /*--- Compute the convective residual using an upwind scheme ---*/
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      
      /*--- Update residual value ---*/
      
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Convective Jacobian contribution for implicit integration ---*/
      if (implicit)
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
      /*--- Roe Turkel preconditioning, set the value of beta ---*/
      if (config->GetKind_Upwind() == TURKEL)
        node[iPoint]->SetPreconditioner_Beta(conv_numerics->GetPrecond_Beta());
      
      /*--- Viscous residual contribution ---*/
      if (viscous) {
        
        /*--- Set laminar and eddy viscosity at the infinity ---*/
        V_infty[nDim+5] = node[iPoint]->GetLaminarViscosity();
        V_infty[nDim+6] = node[iPoint]->GetEddyViscosity();
        
        /*--- Set the normal vector and the coordinates ---*/
        visc_numerics->SetNormal(Normal);
        visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(),
                                geometry->node[Point_Normal]->GetCoord());
        
        /*--- Primitive variables, and gradient ---*/
        visc_numerics->SetPrimitive(V_domain, V_infty);
        visc_numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(),
                                          node[iPoint]->GetGradient_Primitive());
        
        /*--- Turbulent kinetic energy ---*/
        if (config->GetKind_Turb_Model() == SST)
          visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->node[iPoint]->GetSolution(0),
                                              solver_container[TURB_SOL]->node[iPoint]->GetSolution(0));
        
        /*--- Compute and update viscous residual ---*/
        visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
        LinSysRes.SubtractBlock(iPoint, Residual);
        
        /*--- Viscous Jacobian contribution for implicit integration ---*/
        if (implicit)
          Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
        
      }
    }
  }
  
  /*--- Free locally allocated memory ---*/
  delete [] Normal;
  
}

void CEulerSolver::BC_Inlet(CGeometry *geometry, CSolver **solver_container,
                            CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  unsigned short iDim;
  unsigned long iVertex, iPoint, Point_Normal;
  double P_Total, T_Total, Velocity[3], Velocity2, H_Total, Temperature, Riemann,
  Pressure, Density, Energy, *Flow_Dir, Mach2, SoundSpeed2, SoundSpeed_Total2, Vel_Mag,
  alpha, aa, bb, cc, dd, Area, UnitNormal[3];
  double *V_inlet, *V_domain;
  
  bool implicit             = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool grid_movement        = config->GetGrid_Movement();
  double Two_Gamma_M1       = 2.0/Gamma_Minus_One;
  double Gas_Constant       = config->GetGas_ConstantND();
  unsigned short Kind_Inlet = config->GetKind_Inlet();
  string Marker_Tag         = config->GetMarker_All_TagBound(val_marker);
  bool viscous              = config->GetViscous();
  bool tkeNeeded = ((config->GetKind_Solver() == RANS) && (config->GetKind_Turb_Model() == SST));
  
  double *Normal = new double[nDim];
  
  /*--- Loop over all the vertices on this boundary marker ---*/
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    /*--- Allocate the value at the inlet ---*/
    V_inlet = GetCharacPrimVar(val_marker, iVertex);
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Index of the closest interior node ---*/
      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
      
      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);
      
      Area = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
      Area = sqrt (Area);
      
      for (iDim = 0; iDim < nDim; iDim++)
        UnitNormal[iDim] = Normal[iDim]/Area;
      
      /*--- Retrieve solution at this boundary node ---*/
      V_domain = node[iPoint]->GetPrimVar();
      
      /*--- Build the fictitious intlet state based on characteristics ---*/
      
      /*--- Subsonic inflow: there is one outgoing characteristic (u-c),
       therefore we can specify all but one state variable at the inlet.
       The outgoing Riemann invariant provides the final piece of info.
       Adapted from an original implementation in the Stanford University
       multi-block (SUmb) solver in the routine bcSubsonicInflow.f90
       written by Edwin van der Weide, last modified 04-20-2009. ---*/
      
      switch (Kind_Inlet) {
          
          /*--- Total properties have been specified at the inlet. ---*/
        case TOTAL_CONDITIONS:
          
          /*--- Retrieve the specified total conditions for this inlet. ---*/
          P_Total  = config->GetInlet_Ptotal(Marker_Tag);
          T_Total  = config->GetInlet_Ttotal(Marker_Tag);
          Flow_Dir = config->GetInlet_FlowDir(Marker_Tag);
          
          /*--- Non-dim. the inputs if necessary. ---*/
          P_Total /= config->GetPressure_Ref();
          T_Total /= config->GetTemperature_Ref();
          
          /*--- Store primitives and set some variables for clarity. ---*/
          Density = V_domain[nDim+2];
          Velocity2 = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            Velocity[iDim] = V_domain[iDim+1];
            Velocity2 += Velocity[iDim]*Velocity[iDim];
          }
          Energy      = V_domain[nDim+3] - V_domain[nDim+1]/V_domain[nDim+2];
          Pressure    = V_domain[nDim+1];
          H_Total     = (Gamma*Gas_Constant/Gamma_Minus_One)*T_Total;
          SoundSpeed2 = Gamma*Pressure/Density;
          
          /*--- Compute the acoustic Riemann invariant that is extrapolated
           from the domain interior. ---*/
          Riemann   = 2.0*sqrt(SoundSpeed2)/Gamma_Minus_One;
          for (iDim = 0; iDim < nDim; iDim++)
            Riemann += Velocity[iDim]*UnitNormal[iDim];
          
          /*--- Total speed of sound ---*/
          SoundSpeed_Total2 = Gamma_Minus_One*(H_Total - (Energy + Pressure/Density)+0.5*Velocity2) + SoundSpeed2;
          
          /*--- Dot product of normal and flow direction. This should
           be negative due to outward facing boundary normal convention. ---*/
          alpha = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            alpha += UnitNormal[iDim]*Flow_Dir[iDim];
          
          /*--- Coefficients in the quadratic equation for the velocity ---*/
          aa =  1.0 + 0.5*Gamma_Minus_One*alpha*alpha;
          bb = -1.0*Gamma_Minus_One*alpha*Riemann;
          cc =  0.5*Gamma_Minus_One*Riemann*Riemann
          -2.0*SoundSpeed_Total2/Gamma_Minus_One;
          
          /*--- Solve quadratic equation for velocity magnitude. Value must
           be positive, so the choice of root is clear. ---*/
          dd = bb*bb - 4.0*aa*cc;
          dd = sqrt(max(0.0,dd));
          Vel_Mag   = (-bb + dd)/(2.0*aa);
          Vel_Mag   = max(0.0,Vel_Mag);
          Velocity2 = Vel_Mag*Vel_Mag;
          
          /*--- Compute speed of sound from total speed of sound eqn. ---*/
          SoundSpeed2 = SoundSpeed_Total2 - 0.5*Gamma_Minus_One*Velocity2;
          
          /*--- Mach squared (cut between 0-1), use to adapt velocity ---*/
          Mach2 = Velocity2/SoundSpeed2;
          Mach2 = min(1.0,Mach2);
          Velocity2   = Mach2*SoundSpeed2;
          Vel_Mag     = sqrt(Velocity2);
          SoundSpeed2 = SoundSpeed_Total2 - 0.5*Gamma_Minus_One*Velocity2;
          
          /*--- Compute new velocity vector at the inlet ---*/
          for (iDim = 0; iDim < nDim; iDim++)
            Velocity[iDim] = Vel_Mag*Flow_Dir[iDim];
          
          /*--- Static temperature from the speed of sound relation ---*/
          Temperature = SoundSpeed2/(Gamma*Gas_Constant);
          
          /*--- Static pressure using isentropic relation at a point ---*/
          Pressure = P_Total*pow((Temperature/T_Total),Gamma/Gamma_Minus_One);
          
          /*--- Density at the inlet from the gas law ---*/
          Density = Pressure/(Gas_Constant*Temperature);
          
          /*--- Using pressure, density, & velocity, compute the energy ---*/
          Energy = Pressure/(Density*Gamma_Minus_One) + 0.5*Velocity2;
          if (tkeNeeded) Energy += GetTke_Inf();
          
          /*--- Primitive variables, using the derived quantities ---*/
          V_inlet[0] = Temperature;
          for (iDim = 0; iDim < nDim; iDim++)
            V_inlet[iDim+1] = Velocity[iDim];
          V_inlet[nDim+1] = Pressure;
          V_inlet[nDim+2] = Density;
          V_inlet[nDim+3] = Energy + Pressure/Density;
          
          break;
          
          /*--- Mass flow has been specified at the inlet. ---*/
        case MASS_FLOW:
          
          /*--- Retrieve the specified mass flow for the inlet. ---*/
          Density  = config->GetInlet_Ttotal(Marker_Tag);
          Vel_Mag  = config->GetInlet_Ptotal(Marker_Tag);
          Flow_Dir = config->GetInlet_FlowDir(Marker_Tag);
          
          /*--- Non-dim. the inputs if necessary. ---*/
          Density /= config->GetDensity_Ref();
          Vel_Mag /= config->GetVelocity_Ref();
          
          /*--- Get primitives from current inlet state. ---*/
          for (iDim = 0; iDim < nDim; iDim++)
            Velocity[iDim] = node[iPoint]->GetVelocity(iDim);
          Pressure    = node[iPoint]->GetPressure();
          SoundSpeed2 = Gamma*Pressure/V_domain[nDim+2];
          
          /*--- Compute the acoustic Riemann invariant that is extrapolated
           from the domain interior. ---*/
          Riemann = Two_Gamma_M1*sqrt(SoundSpeed2);
          for (iDim = 0; iDim < nDim; iDim++)
            Riemann += Velocity[iDim]*UnitNormal[iDim];
          
          /*--- Speed of sound squared for fictitious inlet state ---*/
          SoundSpeed2 = Riemann;
          for (iDim = 0; iDim < nDim; iDim++)
            SoundSpeed2 -= Vel_Mag*Flow_Dir[iDim]*UnitNormal[iDim];
          
          SoundSpeed2 = max(0.0,0.5*Gamma_Minus_One*SoundSpeed2);
          SoundSpeed2 = SoundSpeed2*SoundSpeed2;
          
          /*--- Pressure for the fictitious inlet state ---*/
          Pressure = SoundSpeed2*Density/Gamma;
          
          /*--- Energy for the fictitious inlet state ---*/
          Energy = Pressure/(Density*Gamma_Minus_One) + 0.5*Vel_Mag*Vel_Mag;
          if (tkeNeeded) Energy += GetTke_Inf();
          
          /*--- Primitive variables, using the derived quantities ---*/
          V_inlet[0] = Pressure / ( Gas_Constant * Density);
          for (iDim = 0; iDim < nDim; iDim++)
            V_inlet[iDim+1] = Vel_Mag*Flow_Dir[iDim];
          V_inlet[nDim+1] = Pressure;
          V_inlet[nDim+2] = Density;
          V_inlet[nDim+3] = Energy + Pressure/Density;
          
          break;
      }
      
      /*--- Set various quantities in the solver class ---*/
      conv_numerics->SetPrimitive(V_domain, V_inlet);
      
      if (grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());
      
      /*--- Compute the residual using an upwind scheme ---*/
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      
      /*--- Update residual value ---*/
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Jacobian contribution for implicit integration ---*/
      if (implicit)
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
      /*--- Roe Turkel preconditioning, set the value of beta ---*/
      if (config->GetKind_Upwind() == TURKEL)
        node[iPoint]->SetPreconditioner_Beta(conv_numerics->GetPrecond_Beta());
      
      /*--- Viscous contribution ---*/
      if (viscous) {
        
        /*--- Set laminar and eddy viscosity at the infinity ---*/
        V_inlet[nDim+5] = node[iPoint]->GetLaminarViscosity();
        V_inlet[nDim+6] = node[iPoint]->GetEddyViscosity();
        
        /*--- Set the normal vector and the coordinates ---*/
        visc_numerics->SetNormal(Normal);
        visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
        
        /*--- Primitive variables, and gradient ---*/
        visc_numerics->SetPrimitive(V_domain, V_inlet);
        visc_numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(), node[iPoint]->GetGradient_Primitive());
        
        /*--- Turbulent kinetic energy ---*/
        if (config->GetKind_Turb_Model() == SST)
          visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->node[iPoint]->GetSolution(0), solver_container[TURB_SOL]->node[iPoint]->GetSolution(0));
        
        /*--- Compute and update residual ---*/
        visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
        LinSysRes.SubtractBlock(iPoint, Residual);
        
        /*--- Jacobian contribution for implicit integration ---*/
        if (implicit)
          Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
        
      }
      
    }
  }
  
  /*--- Free locally allocated memory ---*/
  delete [] Normal;
  
}

void CEulerSolver::BC_Outlet(CGeometry *geometry, CSolver **solver_container,
                             CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  unsigned short iVar, iDim;
  unsigned long iVertex, iPoint, Point_Normal;
  double Pressure, P_Exit, Velocity[3],
  Velocity2, Entropy, Density, Energy, Riemann, Vn, SoundSpeed, Mach_Exit, Vn_Exit,
  Area, UnitNormal[3];
  double *V_outlet, *V_domain;
  
  bool implicit           = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  double Gas_Constant     = config->GetGas_ConstantND();
  bool grid_movement      = config->GetGrid_Movement();
  string Marker_Tag       = config->GetMarker_All_TagBound(val_marker);
  bool viscous              = config->GetViscous();
  bool tkeNeeded = ((config->GetKind_Solver() == RANS) && (config->GetKind_Turb_Model() == SST));
  
  double *Normal = new double[nDim];
  
  /*--- Loop over all the vertices on this boundary marker ---*/
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    /*--- Allocate the value at the outlet ---*/
    V_outlet = GetCharacPrimVar(val_marker, iVertex);
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Index of the closest interior node ---*/
      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
      
      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);
      
      Area = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
      Area = sqrt (Area);
      
      for (iDim = 0; iDim < nDim; iDim++)
        UnitNormal[iDim] = Normal[iDim]/Area;
      
      /*--- Current solution at this boundary node ---*/
      V_domain = node[iPoint]->GetPrimVar();
      
      /*--- Build the fictitious intlet state based on characteristics ---*/
      
      /*--- Retrieve the specified back pressure for this outlet. ---*/
      P_Exit = config->GetOutlet_Pressure(Marker_Tag);
      
      /*--- Non-dim. the inputs if necessary. ---*/
      P_Exit = P_Exit/config->GetPressure_Ref();
      
      /*--- Check whether the flow is supersonic at the exit. The type
       of boundary update depends on this. ---*/
      Density = V_domain[nDim+2];
      Velocity2 = 0.0; Vn = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        Velocity[iDim] = V_domain[iDim+1];
        Velocity2 += Velocity[iDim]*Velocity[iDim];
        Vn += Velocity[iDim]*UnitNormal[iDim];
      }
      Energy     = V_domain[nDim+3] - V_domain[nDim+1]/V_domain[nDim+2];
      Pressure   = V_domain[nDim+1];
      SoundSpeed = sqrt(Gamma*Pressure/Density);
      Mach_Exit  = sqrt(Velocity2)/SoundSpeed;
      
      if (Mach_Exit >= 1.0) {
        
        /*--- Supersonic exit flow: there are no incoming characteristics,
         so no boundary condition is necessary. Set outlet state to current
         state so that upwinding handles the direction of propagation. ---*/
        for (iVar = 0; iVar < nPrimVar; iVar++) V_outlet[iVar] = V_domain[iVar];
        
      } else {
        
        /*--- Subsonic exit flow: there is one incoming characteristic,
         therefore one variable can be specified (back pressure) and is used
         to update the conservative variables. Compute the entropy and the
         acoustic Riemann variable. These invariants, as well as the
         tangential velocity components, are extrapolated. Adapted from an
         original implementation in the Stanford University multi-block
         (SUmb) solver in the routine bcSubsonicOutflow.f90 by Edwin van
         der Weide, last modified 09-10-2007. ---*/
        
        Entropy = Pressure*pow(1.0/Density,Gamma);
        Riemann = Vn + 2.0*SoundSpeed/Gamma_Minus_One;
        
        /*--- Compute the new fictious state at the outlet ---*/
        Density    = pow(P_Exit/Entropy,1.0/Gamma);
        Pressure   = P_Exit;
        SoundSpeed = sqrt(Gamma*P_Exit/Density);
        Vn_Exit    = Riemann - 2.0*SoundSpeed/Gamma_Minus_One;
        Velocity2  = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          Velocity[iDim] = Velocity[iDim] + (Vn_Exit-Vn)*UnitNormal[iDim];
          Velocity2 += Velocity[iDim]*Velocity[iDim];
        }
        Energy = P_Exit/(Density*Gamma_Minus_One) + 0.5*Velocity2;
        if (tkeNeeded) Energy += GetTke_Inf();
        
        /*--- Conservative variables, using the derived quantities ---*/
        V_outlet[0] = Pressure / ( Gas_Constant * Density);
        for (iDim = 0; iDim < nDim; iDim++)
          V_outlet[iDim+1] = Velocity[iDim];
        V_outlet[nDim+1] = Pressure;
        V_outlet[nDim+2] = Density;
        V_outlet[nDim+3] = Energy + Pressure/Density;
        
      }
      
      /*--- Set various quantities in the solver class ---*/
      conv_numerics->SetPrimitive(V_domain, V_outlet);
      
      if (grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());
      
      /*--- Compute the residual using an upwind scheme ---*/
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      
      /*--- Update residual value ---*/
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Jacobian contribution for implicit integration ---*/
      if (implicit) {
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      }
      
      /*--- Roe Turkel preconditioning, set the value of beta ---*/
      if (config->GetKind_Upwind() == TURKEL)
        node[iPoint]->SetPreconditioner_Beta(conv_numerics->GetPrecond_Beta());
      
      /*--- Viscous contribution ---*/
      if (viscous) {
        
        /*--- Set laminar and eddy viscosity at the infinity ---*/
        V_outlet[nDim+5] = node[iPoint]->GetLaminarViscosity();
        V_outlet[nDim+6] = node[iPoint]->GetEddyViscosity();
        
        /*--- Set the normal vector and the coordinates ---*/
        visc_numerics->SetNormal(Normal);
        visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
        
        /*--- Primitive variables, and gradient ---*/
        visc_numerics->SetPrimitive(V_domain, V_outlet);
        visc_numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(), node[iPoint]->GetGradient_Primitive());
        
        /*--- Turbulent kinetic energy ---*/
        if (config->GetKind_Turb_Model() == SST)
          visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->node[iPoint]->GetSolution(0), solver_container[TURB_SOL]->node[iPoint]->GetSolution(0));
        
        /*--- Compute and update residual ---*/
        visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
        LinSysRes.SubtractBlock(iPoint, Residual);
        
        /*--- Jacobian contribution for implicit integration ---*/
        if (implicit)
          Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
        
      }
      
    }
  }
  
  /*--- Free locally allocated memory ---*/
  delete [] Normal;
  
}

void CEulerSolver::BC_Sym_Plane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                                CConfig *config, unsigned short val_marker) {
  
  /*--- Call the Euler residual --- */
  
  BC_Euler_Wall(geometry, solver_container, conv_numerics, config, val_marker);
  
}

void CEulerSolver::SetResidual_DualTime(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                        unsigned short iRKStep, unsigned short iMesh, unsigned short RunTime_EqSystem) {
  
  /*--- Local variables ---*/
  
  unsigned short iVar, jVar, iMarker, iDim;
  unsigned long iPoint, jPoint, iEdge, iVertex;
  
  double *U_time_nM1, *U_time_n, *U_time_nP1;
  double Volume_nM1, Volume_nP1, TimeStep;
  double *Normal = NULL, *GridVel_i = NULL, *GridVel_j = NULL, Residual_GCL;
  
  bool implicit       = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool grid_movement  = config->GetGrid_Movement();
  
  /*--- Store the physical time step ---*/
  
  TimeStep = config->GetDelta_UnstTimeND();
  
  /*--- Compute the dual time-stepping source term for static meshes ---*/
  
  if (!grid_movement) {
    
    /*--- Loop over all nodes (excluding halos) ---*/
    
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      
      /*--- Retrieve the solution at time levels n-1, n, and n+1. Note that
       we are currently iterating on U^n+1 and that U^n & U^n-1 are fixed,
       previous solutions that are stored in memory. ---*/
      
      U_time_nM1 = node[iPoint]->GetSolution_time_n1();
      U_time_n   = node[iPoint]->GetSolution_time_n();
      U_time_nP1 = node[iPoint]->GetSolution();
      
      /*--- CV volume at time n+1. As we are on a static mesh, the volume
       of the CV will remained fixed for all time steps. ---*/
      
      Volume_nP1 = geometry->node[iPoint]->GetVolume();
      
      /*--- Compute the dual time-stepping source term based on the chosen
       time discretization scheme (1st- or 2nd-order).---*/
      
      for (iVar = 0; iVar < nVar; iVar++) {
        if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
          Residual[iVar] = (U_time_nP1[iVar] - U_time_n[iVar])*Volume_nP1 / TimeStep;
        if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
          Residual[iVar] = ( 3.0*U_time_nP1[iVar] - 4.0*U_time_n[iVar]
                            +1.0*U_time_nM1[iVar])*Volume_nP1 / (2.0*TimeStep);
      }
      
      /*--- Store the residual and compute the Jacobian contribution due
       to the dual time source term. ---*/
      
      LinSysRes.AddBlock(iPoint, Residual);
      if (implicit) {
        for (iVar = 0; iVar < nVar; iVar++) {
          for (jVar = 0; jVar < nVar; jVar++) Jacobian_i[iVar][jVar] = 0.0;
          if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
            Jacobian_i[iVar][iVar] = Volume_nP1 / TimeStep;
          if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
            Jacobian_i[iVar][iVar] = (Volume_nP1*3.0)/(2.0*TimeStep);
        }
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      }
    }
    
  }
  
  else {
    
    /*--- For unsteady flows on dynamic meshes (rigidly transforming or
     dynamically deforming), the Geometric Conservation Law (GCL) should be
     satisfied in conjunction with the ALE formulation of the governing
     equations. The GCL prevents accuracy issues caused by grid motion, i.e.
     a uniform free-stream should be preserved through a moving grid. First,
     we will loop over the edges and boundaries to compute the GCL component
     of the dual time source term that depends on grid velocities. ---*/
    
    for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
      
      /*--- Get indices for nodes i & j plus the face normal ---*/
      
      iPoint = geometry->edge[iEdge]->GetNode(0);
      jPoint = geometry->edge[iEdge]->GetNode(1);
      Normal = geometry->edge[iEdge]->GetNormal();
      
      /*--- Grid velocities stored at nodes i & j ---*/
      
      GridVel_i = geometry->node[iPoint]->GetGridVel();
      GridVel_j = geometry->node[jPoint]->GetGridVel();
      
      /*--- Compute the GCL term by averaging the grid velocities at the
       edge mid-point and dotting with the face normal. ---*/
      
      Residual_GCL = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        Residual_GCL += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
      
      /*--- Compute the GCL component of the source term for node i ---*/
      
      U_time_n = node[iPoint]->GetSolution_time_n();
      for(iVar = 0; iVar < nVar; iVar++)
        Residual[iVar] = U_time_n[iVar]*Residual_GCL;
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Compute the GCL component of the source term for node j ---*/
      
      U_time_n = node[jPoint]->GetSolution_time_n();
      for(iVar = 0; iVar < nVar; iVar++)
        Residual[iVar] = U_time_n[iVar]*Residual_GCL;
      LinSysRes.SubtractBlock(jPoint, Residual);
      
    }
    
    /*---	Loop over the boundary edges ---*/
    
    for(iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
      for(iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
        
        /*--- Get the index for node i plus the boundary face normal ---*/
        
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
        
        /*--- Grid velocities stored at boundary node i ---*/
        
        GridVel_i = geometry->node[iPoint]->GetGridVel();
        
        /*--- Compute the GCL term by dotting the grid velocity with the face
         normal. The normal is negated to match the boundary convention. ---*/
        
        Residual_GCL = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          Residual_GCL -= 0.5*(GridVel_i[iDim]+GridVel_i[iDim])*Normal[iDim];
        
        /*--- Compute the GCL component of the source term for node i ---*/
        
        U_time_n = node[iPoint]->GetSolution_time_n();
        for(iVar = 0; iVar < nVar; iVar++)
          Residual[iVar] = U_time_n[iVar]*Residual_GCL;
        LinSysRes.AddBlock(iPoint, Residual);
        
      }
    }
    
    /*--- Loop over all nodes (excluding halos) to compute the remainder
     of the dual time-stepping source term. ---*/
    
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      
      /*--- Retrieve the solution at time levels n-1, n, and n+1. Note that
       we are currently iterating on U^n+1 and that U^n & U^n-1 are fixed,
       previous solutions that are stored in memory. ---*/
      
      U_time_nM1 = node[iPoint]->GetSolution_time_n1();
      U_time_n   = node[iPoint]->GetSolution_time_n();
      U_time_nP1 = node[iPoint]->GetSolution();
      
      /*--- CV volume at time n-1 and n+1. In the case of dynamically deforming
       grids, the volumes will change. On rigidly transforming grids, the
       volumes will remain constant. ---*/
      
      Volume_nM1 = geometry->node[iPoint]->GetVolume_nM1();
      Volume_nP1 = geometry->node[iPoint]->GetVolume();
      
      /*--- Compute the dual time-stepping source residual. Due to the
       introduction of the GCL term above, the remainder of the source residual
       due to the time discretization has a new form.---*/
      
      for(iVar = 0; iVar < nVar; iVar++) {
        if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
          Residual[iVar] = (U_time_nP1[iVar] - U_time_n[iVar])*(Volume_nP1/TimeStep);
        if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
          Residual[iVar] = (U_time_nP1[iVar] - U_time_n[iVar])*(3.0*Volume_nP1/(2.0*TimeStep))
          + (U_time_nM1[iVar] - U_time_n[iVar])*(Volume_nM1/(2.0*TimeStep));
      }
      
      /*--- Store the residual and compute the Jacobian contribution due
       to the dual time source term. ---*/
      
      LinSysRes.AddBlock(iPoint, Residual);
      if (implicit) {
        for (iVar = 0; iVar < nVar; iVar++) {
          for (jVar = 0; jVar < nVar; jVar++) Jacobian_i[iVar][jVar] = 0.0;
          if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
            Jacobian_i[iVar][iVar] = Volume_nP1/TimeStep;
          if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
            Jacobian_i[iVar][iVar] = (3.0*Volume_nP1)/(2.0*TimeStep);
        }
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      }
    }
  }
  
}

void CEulerSolver::LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter) {
  
  /*--- Restart the solution from file information ---*/
  unsigned short iDim, iVar, iMesh, iMeshFine;
  unsigned long iPoint, index, iChildren, Point_Fine;
  unsigned short turb_model = config->GetKind_Turb_Model();
  double Area_Children, Area_Parent, Coord[3], *Solution_Fine, dull_val;
  bool grid_movement  = config->GetGrid_Movement();
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  string UnstExt, text_line;
  ifstream restart_file;
  
  string restart_filename = config->GetSolution_FlowFileName();
  
  /*--- Modify file name for an unsteady restart ---*/
  if (dual_time)
    restart_filename = config->GetUnsteady_FileName(restart_filename, val_iter);
  
  /*--- Open the restart file, and throw an error if this fails. ---*/
  restart_file.open(restart_filename.data(), ios::in);
  if (restart_file.fail()) {
    cout << "There is no flow restart file!! " << restart_filename.data() << "."<< endl;
    exit(1);
  }
  
  /*--- In case this is a parallel simulation, we need to perform the
   Global2Local index transformation first. ---*/
  long *Global2Local = NULL;
  Global2Local = new long[geometry[MESH_0]->GetGlobal_nPointDomain()];
  /*--- First, set all indices to a negative value by default ---*/
  for(iPoint = 0; iPoint < geometry[MESH_0]->GetGlobal_nPointDomain(); iPoint++) {
    Global2Local[iPoint] = -1;
  }
  
  /*--- Now fill array with the transform values only for local points ---*/
  for(iPoint = 0; iPoint < geometry[MESH_0]->GetnPointDomain(); iPoint++) {
    Global2Local[geometry[MESH_0]->node[iPoint]->GetGlobalIndex()] = iPoint;
  }
  
  /*--- Read all lines in the restart file ---*/
  long iPoint_Local = 0; unsigned long iPoint_Global = 0;
  
  /*--- The first line is the header ---*/
  getline (restart_file, text_line);
  
  while (getline (restart_file,text_line)) {
    istringstream point_line(text_line);
    
    /*--- Retrieve local index. If this node from the restart file lives
     on a different processor, the value of iPoint_Local will be -1, as
     initialized above. Otherwise, the local index for this node on the
     current processor will be returned and used to instantiate the vars. ---*/
    iPoint_Local = Global2Local[iPoint_Global];
    if (iPoint_Local >= 0) {
      
      if (nDim == 2) point_line >> index >> Coord[0] >> Coord[1] >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3];
      if (nDim == 3) point_line >> index >> Coord[0] >> Coord[1] >> Coord[2] >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3] >> Solution[4];
      
      node[iPoint_Local]->SetSolution(Solution);
      
      /*--- For dynamic meshes, read in and store the
       grid coordinates and grid velocities for each node. ---*/
      if (grid_movement) {
        
        /*--- First, remove any variables for the turbulence model that
         appear in the restart file before the grid velocities. ---*/
        if (turb_model == SA || turb_model == ML) {
          point_line >> dull_val;
        } else if (turb_model == SST) {
          point_line >> dull_val >> dull_val;
        }
        
        /*--- Read in the next 2 or 3 variables which are the grid velocities ---*/
        double GridVel[3];
        if (nDim == 2) point_line >> GridVel[0] >> GridVel[1];
        if (nDim == 3) point_line >> GridVel[0] >> GridVel[1] >> GridVel[2];
        for (iDim = 0; iDim < nDim; iDim++) {
          geometry[MESH_0]->node[iPoint_Local]->SetCoord(iDim, Coord[iDim]);
          geometry[MESH_0]->node[iPoint_Local]->SetGridVel(iDim, GridVel[iDim]);
        }
      }
      
    }
    iPoint_Global++;
  }
  
  /*--- Close the restart file ---*/
  restart_file.close();
  
  /*--- Free memory needed for the transformation ---*/
  delete [] Global2Local;
  
  /*--- Interpolate the solution down to the coarse multigrid levels ---*/
  for (iMesh = 1; iMesh <= config->GetMGLevels(); iMesh++) {
    for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
      Area_Parent = geometry[iMesh]->node[iPoint]->GetVolume();
      for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = 0.0;
      for (iChildren = 0; iChildren < geometry[iMesh]->node[iPoint]->GetnChildren_CV(); iChildren++) {
        Point_Fine = geometry[iMesh]->node[iPoint]->GetChildren_CV(iChildren);
        Area_Children = geometry[iMesh-1]->node[Point_Fine]->GetVolume();
        Solution_Fine = solver[iMesh-1][FLOW_SOL]->node[Point_Fine]->GetSolution();
        for (iVar = 0; iVar < nVar; iVar++) {
          Solution[iVar] += Solution_Fine[iVar]*Area_Children/Area_Parent;
        }
      }
      solver[iMesh][FLOW_SOL]->node[iPoint]->SetSolution(Solution);
    }
  }
  
  /*--- Update the geometry for flows on dynamic meshes ---*/
  if (grid_movement) {
        
    /*--- Recompute the edges and  dual mesh control volumes in the
     domain and on the boundaries. ---*/
    geometry[MESH_0]->SetCG();
    geometry[MESH_0]->SetControlVolume(config, UPDATE);
    geometry[MESH_0]->SetBoundControlVolume(config, UPDATE);
    
    /*--- Update the multigrid structure after setting up the finest grid,
     including computing the grid velocities on the coarser levels. ---*/
    
    for (iMesh = 1; iMesh <= config->GetMGLevels(); iMesh++) {
      iMeshFine = iMesh-1;
      geometry[iMesh]->SetControlVolume(config,geometry[iMeshFine], UPDATE);
      geometry[iMesh]->SetBoundControlVolume(config,geometry[iMeshFine],UPDATE);
      geometry[iMesh]->SetCoord(geometry[iMeshFine]);
      geometry[iMesh]->SetRestricted_GridVelocity(geometry[iMeshFine],config);
    }
  }
  
}

CNSSolver::CNSSolver(void) : CEulerSolver() {
  
  /*--- Array initialization ---*/
  CDrag_Visc = NULL;
  CLift_Visc = NULL;
  CSideForce_Visc = NULL;
  CEff_Visc = NULL;
  CMx_Visc = NULL;
  CMy_Visc = NULL;
  CMz_Visc = NULL;
  CFx_Visc = NULL;
  CFy_Visc = NULL;
  CFz_Visc = NULL;
  Surface_CLift_Visc = NULL;
  Surface_CDrag_Visc = NULL;
  Surface_CMx_Visc = NULL;
  Surface_CMy_Visc = NULL;
  Surface_CMz_Visc = NULL;
  CMerit_Visc = NULL;
  CT_Visc = NULL;
  CQ_Visc = NULL;
  ForceViscous = NULL;
  MomentViscous = NULL;
  CSkinFriction = NULL;
  
}

CNSSolver::CNSSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh) : CEulerSolver() {
  
  unsigned long iPoint, index, counter_local = 0, counter_global = 0, iVertex;
  unsigned short iVar, iDim, iMarker, nLineLets;
  double Density, Velocity2, Pressure, Temperature, dull_val;
  
  unsigned short nZone = geometry->GetnZone();
  bool restart = (config->GetRestart() || config->GetRestart_Flow());
  double Gas_Constant = config->GetGas_ConstantND();
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  
  /*--- Array initialization ---*/
  CDrag_Visc = NULL;
  CLift_Visc = NULL;
  CSideForce_Visc = NULL;
  CEff_Visc = NULL;
  CMx_Visc = NULL;
  CMy_Visc = NULL;
  CMz_Visc = NULL;
  CFx_Visc = NULL;
  CFy_Visc = NULL;
  CFz_Visc = NULL;
  Surface_CLift_Visc = NULL;
  Surface_CDrag_Visc = NULL;
  Surface_CMx_Visc = NULL;
  Surface_CMy_Visc = NULL;
  Surface_CMz_Visc = NULL;
  CMerit_Visc = NULL;
  CT_Visc = NULL;
  CQ_Visc = NULL;
  Heat_Visc = NULL;
  MaxHeatFlux_Visc = NULL;
  ForceViscous = NULL;
  MomentViscous = NULL;
  CSkinFriction = NULL;
  
  int rank = MASTER_NODE;
  
  /*--- Set the gamma value ---*/
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
  /*--- Define geometry constants in the solver structure
   Compressible flow, primitive variables nDim+5, (T,vx,vy,vz,P,rho,h,c,lamMu,EddyMu) ---*/
  nDim = geometry->GetnDim();
  nVar = nDim+2; nPrimVar = nDim+7; nPrimVarGrad = nDim+4;
  nMarker      = config->GetnMarker_All();
  nPoint       = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();
  
  /*--- Allocate the node variables ---*/
  
  node = new CVariable*[nPoint];
  
  /*--- Define some auxiliar vector related with the residual ---*/
  Residual      = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual[iVar]      = 0.0;
  Residual_RMS  = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar]  = 0.0;
  Residual_Max  = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_Max[iVar]  = 0.0;
  Residual_i    = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_i[iVar]    = 0.0;
  Residual_j    = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_j[iVar]    = 0.0;
  Res_Conv      = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Res_Conv[iVar]      = 0.0;
  Res_Visc      = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Res_Visc[iVar]      = 0.0;
  Res_Sour      = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Res_Sour[iVar]      = 0.0;
  Point_Max  = new unsigned long[nVar]; for (iVar = 0; iVar < nVar; iVar++) Point_Max[iVar]  = 0;
  
  /*--- Define some auxiliary vectors related to the solution ---*/
  Solution   = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution[iVar]   = 0.0;
  Solution_i = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution_i[iVar] = 0.0;
  Solution_j = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution_j[iVar] = 0.0;
  
  /*--- Define some auxiliary vectors related to the geometry ---*/
  Vector   = new double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector[iDim]   = 0.0;
  Vector_i = new double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector_i[iDim] = 0.0;
  Vector_j = new double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector_j[iDim] = 0.0;
  
  /*--- Define some auxiliary vectors related to the primitive solution ---*/
  Primitive   = new double[nPrimVar]; for (iVar = 0; iVar < nPrimVar; iVar++) Primitive[iVar]   = 0.0;
  Primitive_i = new double[nPrimVar]; for (iVar = 0; iVar < nPrimVar; iVar++) Primitive_i[iVar] = 0.0;
  Primitive_j = new double[nPrimVar]; for (iVar = 0; iVar < nPrimVar; iVar++) Primitive_j[iVar] = 0.0;
  
  /*--- Define some auxiliar vector related with the undivided lapalacian computation ---*/
  if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED) {
    iPoint_UndLapl = new double [nPoint];
    jPoint_UndLapl = new double [nPoint];
  }
  
  /*--- Define some auxiliary vectors related to low-speed preconditioning ---*/
  if (config->GetKind_Upwind_Flow() == TURKEL) {
    LowMach_Precontioner = new double* [nVar];
    for (iVar = 0; iVar < nVar; iVar ++)
      LowMach_Precontioner[iVar] = new double[nVar];
  }
  
  /*--- Initialize the solution and right hand side vectors for storing
   the residuals and updating the solution (always needed even for
   explicit schemes). ---*/
  LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
  
  /*--- Jacobians and vector structures for implicit computations ---*/
  if (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT) {
    
    /*--- Point to point Jacobians ---*/
    Jacobian_i = new double* [nVar];
    Jacobian_j = new double* [nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Jacobian_i[iVar] = new double [nVar];
      Jacobian_j[iVar] = new double [nVar];
    }
    /*--- Initialization of the structure of the whole Jacobian ---*/
    if (rank == MASTER_NODE) cout << "Initialize jacobian structure (Navier-Stokes). MG level: " << iMesh <<"." << endl;
    Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config);
    
    if ((config->GetKind_Linear_Solver_Prec() == LINELET) ||
        (config->GetKind_Linear_Solver() == SMOOTHER_LINELET)) {
      nLineLets = Jacobian.BuildLineletPreconditioner(geometry, config);
      if (rank == MASTER_NODE) cout << "Compute linelet structure. " << nLineLets << " elements in each line (average)." << endl;
    }
    
  } else {
    if (rank == MASTER_NODE)
      cout << "Explicit scheme. No jacobian structure (Navier-Stokes). MG level: " << iMesh <<"." << endl;
  }
  
  /*--- Define some auxiliary vectors for computing flow variable gradients by least squares ---*/
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
    
    /*--- S matrix := inv(R)*traspose(inv(R)) ---*/
    Smatrix = new double* [nDim];
    for (iDim = 0; iDim < nDim; iDim++)
      Smatrix[iDim] = new double [nDim];
    
    /*--- c vector := transpose(WA)*(Wb) ---*/
    cvector = new double* [nPrimVarGrad];
    for (iVar = 0; iVar < nPrimVarGrad; iVar++)
      cvector[iVar] = new double [nDim];
  }
  
  /*--- Store the value of the characteristic primitive variables at the boundaries ---*/
  CharacPrimVar = new double** [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    CharacPrimVar[iMarker] = new double* [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      CharacPrimVar[iMarker][iVertex] = new double [nPrimVar];
      for (iVar = 0; iVar < nPrimVar; iVar++) {
        CharacPrimVar[iMarker][iVertex][iVar] = 0.0;
      }
    }
  }
  
  /*--- Inviscid force definition and coefficient in all the markers ---*/
  CPressure = new double* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    CPressure[iMarker] = new double [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      CPressure[iMarker][iVertex] = 0.0;
    }
  }
  
  /*--- Inviscid force definition and coefficient in all the markers ---*/
  CPressureTarget = new double* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    CPressureTarget[iMarker] = new double [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      CPressureTarget[iMarker][iVertex] = 0.0;
    }
  }
  
  /*--- Heat tranfer in all the markers ---*/
  HeatFlux = new double* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    HeatFlux[iMarker] = new double [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      HeatFlux[iMarker][iVertex] = 0.0;
    }
  }
  
  /*--- Heat tranfer in all the markers ---*/
  HeatFluxTarget = new double* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    HeatFluxTarget[iMarker] = new double [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      HeatFluxTarget[iMarker][iVertex] = 0.0;
    }
  }
  
  /*--- Y plus in all the markers ---*/
  YPlus = new double* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    YPlus[iMarker] = new double [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      YPlus[iMarker][iVertex] = 0.0;
    }
  }
  
  /*--- Skin friction in all the markers ---*/
  CSkinFriction = new double* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    CSkinFriction[iMarker] = new double [geometry->nVertex[iMarker]];
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      CSkinFriction[iMarker][iVertex] = 0.0;
    }
  }
  
  /*--- Non dimensional coefficients ---*/
  ForceInviscid  = new double[3];
  MomentInviscid = new double[3];
  CDrag_Inv      = new double[nMarker];
  CLift_Inv      = new double[nMarker];
  CSideForce_Inv = new double[nMarker];
  CMx_Inv        = new double[nMarker];
  CMy_Inv        = new double[nMarker];
  CMz_Inv        = new double[nMarker];
  CEff_Inv       = new double[nMarker];
  CFx_Inv        = new double[nMarker];
  CFy_Inv        = new double[nMarker];
  CFz_Inv        = new double[nMarker];
  
  Surface_CLift_Inv= new double[config->GetnMarker_Monitoring()];
  Surface_CDrag_Inv= new double[config->GetnMarker_Monitoring()];
  Surface_CMx_Inv  = new double[config->GetnMarker_Monitoring()];
  Surface_CMy_Inv  = new double[config->GetnMarker_Monitoring()];
  Surface_CMz_Inv  = new double[config->GetnMarker_Monitoring()];
  Surface_CLift    = new double[config->GetnMarker_Monitoring()];
  Surface_CDrag    = new double[config->GetnMarker_Monitoring()];
  Surface_CMx      = new double[config->GetnMarker_Monitoring()];
  Surface_CMy      = new double[config->GetnMarker_Monitoring()];
  Surface_CMz      = new double[config->GetnMarker_Monitoring()];
  
  /*--- Rotational coefficients ---*/
  CMerit_Inv = new double[nMarker];
  CT_Inv     = new double[nMarker];
  CQ_Inv     = new double[nMarker];
  
  /*--- Supersonic coefficients ---*/
  CEquivArea_Inv   = new double[nMarker];
  CNearFieldOF_Inv = new double[nMarker];
  
  /*--- Nacelle simulation ---*/
  FanFace_MassFlow  = new double[nMarker];
  Exhaust_MassFlow  = new double[nMarker];
  Exhaust_Area      = new double[nMarker];
  FanFace_Pressure  = new double[nMarker];
  FanFace_Mach      = new double[nMarker];
  FanFace_Area      = new double[nMarker];
  
  /*--- Init total coefficients ---*/
  Total_CDrag   = 0.0;	Total_CLift       = 0.0;  Total_CSideForce   = 0.0;
  Total_CMx     = 0.0;	Total_CMy         = 0.0;  Total_CMz          = 0.0;
  Total_CEff    = 0.0;	Total_CEquivArea  = 0.0;  Total_CNearFieldOF = 0.0;
  Total_CFx     = 0.0;	Total_CFy         = 0.0;  Total_CFz          = 0.0;
  Total_CT      = 0.0;	Total_CQ          = 0.0;  Total_CMerit       = 0.0;
  Total_MaxHeat = 0.0;  Total_Heat        = 0.0;
  Total_CpDiff  = 0.0;  Total_HeatFluxDiff    = 0.0;
  
  ForceViscous    = new double[3];
  MomentViscous   = new double[3];
  CDrag_Visc      = new double[nMarker];
  CLift_Visc      = new double[nMarker];
  CSideForce_Visc = new double[nMarker];
  CMx_Visc        = new double[nMarker];
  CMy_Visc        = new double[nMarker];
  CMz_Visc        = new double[nMarker];
  CEff_Visc       = new double[nMarker];
  CFx_Visc        = new double[nMarker];
  CFy_Visc        = new double[nMarker];
  CFz_Visc        = new double[nMarker];
  CMerit_Visc     = new double[nMarker];
  CT_Visc         = new double[nMarker];
  CQ_Visc         = new double[nMarker];
  Heat_Visc       = new double[nMarker];
  MaxHeatFlux_Visc   = new double[nMarker];
  
  Surface_CLift_Visc = new double[config->GetnMarker_Monitoring()];
  Surface_CDrag_Visc = new double[config->GetnMarker_Monitoring()];
  Surface_CMx_Visc = new double[config->GetnMarker_Monitoring()];
  Surface_CMy_Visc = new double[config->GetnMarker_Monitoring()];
  Surface_CMz_Visc = new double[config->GetnMarker_Monitoring()];
  
  /*--- Read farfield conditions from config ---*/
  Density_Inf   = config->GetDensity_FreeStreamND();
  Pressure_Inf  = config->GetPressure_FreeStreamND();
  Velocity_Inf  = config->GetVelocity_FreeStreamND();
  Energy_Inf    = config->GetEnergy_FreeStreamND();
  Viscosity_Inf = config->GetViscosity_FreeStreamND();
  Mach_Inf      = config->GetMach_FreeStreamND();
  Prandtl_Lam   = config->GetPrandtl_Lam();
  Prandtl_Turb  = config->GetPrandtl_Turb();
  Tke_Inf       = config->GetTke_FreeStreamND();
  
  /*--- Initializate Fan Face Pressure ---*/
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    FanFace_MassFlow[iMarker] = 0.0;
    Exhaust_MassFlow[iMarker] = 0.0;
    Exhaust_Area[iMarker] = 0.0;
    FanFace_Pressure[iMarker] = Pressure_Inf;
    FanFace_Mach[iMarker] = Mach_Inf;
  }
  
  /*--- Check for a restart and set up the variables at each node
   appropriately. Coarse multigrid levels will be intitially set to
   the farfield values bc the solver will immediately interpolate
   the solution from the finest mesh to the coarser levels. ---*/
  
  if (!restart || geometry->GetFinestMGLevel() == false || nZone > 1) {
    
    /*--- Restart the solution from the free-stream state ---*/
    for (iPoint = 0; iPoint < nPoint; iPoint++)
      node[iPoint] = new CNSVariable(Density_Inf, Velocity_Inf, Energy_Inf, nDim, nVar, config);
  }
  
  else {
    
    /*--- Initialize the solution from the restart file information ---*/
    ifstream restart_file;
    string filename = config->GetSolution_FlowFileName();
    
    /*--- Modify file name for an unsteady restart ---*/
    if (dual_time) {
      int Unst_RestartIter;
      if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
        Unst_RestartIter = int(config->GetUnst_RestartIter())-1;
      else
        Unst_RestartIter = int(config->GetUnst_RestartIter())-2;
      filename = config->GetUnsteady_FileName(filename, Unst_RestartIter);
    }
    
    /*--- Open the restart file, throw an error if this fails. ---*/
    restart_file.open(filename.data(), ios::in);
    if (restart_file.fail()) {
      cout << "There is no flow restart file!! " << filename.data() << "."<< endl;
      exit(1);
    }
    
    /*--- In case this is a parallel simulation, we need to perform the
     Global2Local index transformation first. ---*/
    long *Global2Local = new long[geometry->GetGlobal_nPointDomain()];
    
    /*--- First, set all indices to a negative value by default ---*/
    for(iPoint = 0; iPoint < geometry->GetGlobal_nPointDomain(); iPoint++)
      Global2Local[iPoint] = -1;
    
    /*--- Now fill array with the transform values only for local points ---*/
    for(iPoint = 0; iPoint < nPointDomain; iPoint++)
      Global2Local[geometry->node[iPoint]->GetGlobalIndex()] = iPoint;
    
    /*--- Read all lines in the restart file ---*/
    long iPoint_Local; unsigned long iPoint_Global = 0; string text_line;
    
    /*--- The first line is the header ---*/
    getline (restart_file, text_line);
    
    while (getline (restart_file,text_line)) {
      istringstream point_line(text_line);
      
      /*--- Retrieve local index. If this node from the restart file lives
       on a different processor, the value of iPoint_Local will be -1.
       Otherwise, the local index for this node on the current processor
       will be returned and used to instantiate the vars. ---*/
      iPoint_Local = Global2Local[iPoint_Global];
      
      /*--- Load the solution for this node. Note that the first entry
       on the restart file line is the global index, followed by the
       node coordinates, and then the conservative variables. ---*/
      if (iPoint_Local >= 0) {
        if (nDim == 2) point_line >> index >> dull_val >> dull_val >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3];
        if (nDim == 3) point_line >> index >> dull_val >> dull_val >> dull_val >> Solution[0] >> Solution[1] >> Solution[2] >> Solution[3] >> Solution[4];
        node[iPoint_Local] = new CNSVariable(Solution, nDim, nVar, config);
      }
      iPoint_Global++;
    }
    
    /*--- Instantiate the variable class with an arbitrary solution
     at any halo/periodic nodes. The initial solution can be arbitrary,
     because a send/recv is performed immediately in the solver. ---*/
    for(iPoint = nPointDomain; iPoint < nPoint; iPoint++)
      node[iPoint] = new CNSVariable(Solution, nDim, nVar, config);
    
    /*--- Close the restart file ---*/
    restart_file.close();
    
    /*--- Free memory needed for the transformation ---*/
    delete [] Global2Local;
  }
  
  /*--- Check that the initial solution is physical, report any non-physical nodes ---*/
  counter_local = 0;
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    Density = node[iPoint]->GetSolution(0);
    Velocity2 = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      Velocity2 += (node[iPoint]->GetSolution(iDim+1)/Density)*(node[iPoint]->GetSolution(iDim+1)/Density);
    Pressure    = Gamma_Minus_One*Density*(node[iPoint]->GetSolution(nDim+1)/Density-0.5*Velocity2);
    Temperature = Pressure / ( Gas_Constant * Density);
    if ((Pressure < 0.0) || (Temperature < 0.0)) {
      Solution[0] = Density_Inf;
      for (iDim = 0; iDim < nDim; iDim++)
        Solution[iDim+1] = Velocity_Inf[iDim]*Density_Inf;
      Solution[nDim+1] = Energy_Inf*Density_Inf;
      node[iPoint]->SetSolution(Solution);
      node[iPoint]->SetSolution_Old(Solution);
      counter_local++;
    }
  }
  counter_global = counter_local;

  if ((rank == MASTER_NODE) && (counter_global != 0))
    cout << "Warning. The original solution contains "<< counter_global << " points that are not physical." << endl;
  
  /*--- For incompressible solver set the initial values for the density and viscosity,
   unless a freesurface problem, this must be constant during the computation ---*/
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    node[iPoint]->SetDensityInc(Density_Inf);
    node[iPoint]->SetLaminarViscosityInc(Viscosity_Inf);
  }
  
  /*--- Define solver parameters needed for execution of destructor ---*/
  if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED)
    space_centered = true;
  else space_centered = false;
  
  if (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT) euler_implicit = true;
  else euler_implicit = false;
  
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) least_squares = true;
  else least_squares = false;
    
}

CNSSolver::~CNSSolver(void) {
  unsigned short iMarker;
  
  if (CDrag_Visc != NULL)      delete [] CDrag_Visc;
  if (CLift_Visc != NULL)      delete [] CLift_Visc;
  if (CSideForce_Visc != NULL) delete [] CSideForce_Visc;
  if (CMx_Visc != NULL)        delete [] CMx_Visc;
  if (CMy_Visc != NULL)        delete [] CMy_Visc;
  if (CMz_Visc != NULL)        delete [] CMz_Visc;
  if (CFx_Visc != NULL)        delete [] CFx_Visc;
  if (CFy_Visc != NULL)        delete [] CFy_Visc;
  if (CFz_Visc != NULL)        delete [] CFz_Visc;
  if (CEff_Visc != NULL)       delete [] CEff_Visc;
  if (CMerit_Visc != NULL)     delete [] CMerit_Visc;
  if (CT_Visc != NULL)         delete [] CT_Visc;
  if (CQ_Visc != NULL)         delete [] CQ_Visc;
  if (Heat_Visc != NULL)          delete [] Heat_Visc;
  if (MaxHeatFlux_Visc != NULL)       delete [] MaxHeatFlux_Visc;
  if (ForceViscous != NULL)    delete [] ForceViscous;
  if (MomentViscous != NULL)   delete [] MomentViscous;
  
  if (Surface_CLift_Visc != NULL) delete [] Surface_CLift_Visc;
  if (Surface_CDrag_Visc != NULL) delete [] Surface_CDrag_Visc;
  if (Surface_CMx_Visc != NULL)   delete [] Surface_CMx_Visc;
  if (Surface_CMy_Visc != NULL)   delete [] Surface_CMy_Visc;
  if (Surface_CMz_Visc != NULL)   delete [] Surface_CMz_Visc;
  
  if (CSkinFriction != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      delete CSkinFriction[iMarker];
    }
    delete [] CSkinFriction;
  }
  
}

void CNSSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {
  
  unsigned long iPoint, ErrorCounter = 0;
  bool RightSol;
  int rank;
  
  rank = MASTER_NODE;
  
  unsigned long ExtIter     = config->GetExtIter();
  bool adjoint = config->GetAdjoint();
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool center = (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED) || (adjoint && config->GetKind_ConvNumScheme_AdjFlow() == SPACE_CENTERED);
  bool center_jst = center && config->GetKind_Centered_Flow() == JST;
  bool limiter_flow = ((config->GetSpatialOrder_Flow() == SECOND_ORDER_LIMITER) && (ExtIter <= config->GetLimiterIter()));
  bool limiter_turb = ((config->GetSpatialOrder_Turb() == SECOND_ORDER_LIMITER) && (ExtIter <= config->GetLimiterIter()));
  unsigned short turb_model = config->GetKind_Turb_Model();
  bool tkeNeeded = (turb_model == SST);
  double eddy_visc = 0.0, turb_ke = 0.0;
  
  for (iPoint = 0; iPoint < nPoint; iPoint ++) {
    
    if (turb_model != NONE) {
      eddy_visc = solver_container[TURB_SOL]->node[iPoint]->GetmuT();
      if (tkeNeeded) turb_ke = solver_container[TURB_SOL]->node[iPoint]->GetSolution(0);
    }
    
    /*--- Incompressible flow, primitive variables nDim+3, (P,vx,vy,vz,rho,beta),
     FreeSurface Incompressible flow, primitive variables nDim+4, (P,vx,vy,vz,rho,beta,dist),
     Compressible flow, primitive variables nDim+5, (T,vx,vy,vz,P,rho,h,c) ---*/
    RightSol = node[iPoint]->SetPrimVar_Compressible(eddy_visc, turb_ke, config);
    if (!RightSol) ErrorCounter++;
    
  }
  
  for (iPoint = 0; iPoint < nPoint; iPoint ++) {
    
    /*--- Initialize the convective, source and viscous residual vector ---*/
    if (!Output) LinSysRes.SetBlock_Zero(iPoint);
    
  }
  
  /*--- Artificial dissipation ---*/
  if (center) {
    SetMax_Eigenvalue(geometry, config);
    if ((center_jst) && (iMesh == MESH_0)) {
      SetDissipation_Switch(geometry, config);
      SetUndivided_Laplacian(geometry, config);
    }
  }
  
  /*--- Compute gradient of the primitive variables ---*/
  if (config->GetKind_Gradient_Method() == GREEN_GAUSS) {
    SetPrimVar_Gradient_GG(geometry, config);
  }
  
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
    SetPrimVar_Gradient_LS(geometry, config);
  }
  
  /*--- Compute the limiter in case we need it in the turbulence model
   or to limit the viscous terms (check this logic with JST and 2nd order turbulence model) ---*/
  if ((iMesh == MESH_0) && (limiter_flow || limiter_turb)) {
    SetPrimVar_Limiter(geometry, config);
  }
  
  /*--- Initialize the jacobian matrices ---*/
  if (implicit) {
    Jacobian.SetValZero();
  }
  
  /*--- Error message ---*/
  if (Output && (ErrorCounter >= 10) && (rank == MASTER_NODE) && (iMesh == MESH_0))
    cout <<"The solution contains "<< ErrorCounter << " non-physical points." << endl;
  
}

void CNSSolver::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned long Iteration) {
  double *Normal, Area, Vol, Mean_SoundSpeed, Mean_ProjVel, Lambda, Local_Delta_Time, Local_Delta_Time_Visc,
  Global_Delta_Time = 1E6, Mean_LaminarVisc, Mean_EddyVisc, Mean_Density, Lambda_1, Lambda_2, K_v = 0.25, Global_Delta_UnstTimeND;
  unsigned long iEdge, iVertex, iPoint = 0, jPoint = 0;
  unsigned short iDim, iMarker;
  double ProjVel, ProjVel_i, ProjVel_j;
  
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool grid_movement = config->GetGrid_Movement();
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  
  Min_Delta_Time = 1.E6; Max_Delta_Time = 0.0;
  
  /*--- Set maximum inviscid eigenvalue to zero, and compute sound speed and viscosity ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    node[iPoint]->SetMax_Lambda_Inv(0.0);
    node[iPoint]->SetMax_Lambda_Visc(0.0);
  }
  
  /*--- Loop interior edges ---*/
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    /*--- Point identification, Normal vector and area ---*/
    
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    
    Normal = geometry->edge[iEdge]->GetNormal();
    Area = 0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);
    
    /*--- Mean Values ---*/
    
    Mean_ProjVel = 0.5 * (node[iPoint]->GetProjVel(Normal) + node[jPoint]->GetProjVel(Normal));
    Mean_SoundSpeed = 0.5 * (node[iPoint]->GetSoundSpeed() + node[jPoint]->GetSoundSpeed()) * Area;
    
    /*--- Adjustment for grid movement ---*/
    
    if (grid_movement) {
      double *GridVel_i = geometry->node[iPoint]->GetGridVel();
      double *GridVel_j = geometry->node[jPoint]->GetGridVel();
      ProjVel_i = 0.0; ProjVel_j =0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        ProjVel_i += GridVel_i[iDim]*Normal[iDim];
        ProjVel_j += GridVel_j[iDim]*Normal[iDim];
      }
      Mean_ProjVel -= 0.5 * (ProjVel_i + ProjVel_j) ;
    }
    
    /*--- Inviscid contribution ---*/
    
    Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed ;
    if (geometry->node[iPoint]->GetDomain()) node[iPoint]->AddMax_Lambda_Inv(Lambda);
    if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddMax_Lambda_Inv(Lambda);
    
    /*--- Viscous contribution ---*/
    
    Mean_LaminarVisc = 0.5*(node[iPoint]->GetLaminarViscosity() + node[jPoint]->GetLaminarViscosity());
    Mean_EddyVisc    = 0.5*(node[iPoint]->GetEddyViscosity() + node[jPoint]->GetEddyViscosity());
    Mean_Density     = 0.5*(node[iPoint]->GetSolution(0) + node[jPoint]->GetSolution(0));
    
    Lambda_1 = (4.0/3.0)*(Mean_LaminarVisc + Mean_EddyVisc);
    Lambda_2 = (1.0 + (Prandtl_Lam/Prandtl_Turb)*(Mean_EddyVisc/Mean_LaminarVisc))*(Gamma*Mean_LaminarVisc/Prandtl_Lam);
    Lambda = (Lambda_1 + Lambda_2)*Area*Area/Mean_Density;
    
    if (geometry->node[iPoint]->GetDomain()) node[iPoint]->AddMax_Lambda_Visc(Lambda);
    if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddMax_Lambda_Visc(Lambda);
    
  }
  
  /*--- Loop boundary edges ---*/
  
  for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
    for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
      
      /*--- Point identification, Normal vector and area ---*/
      
      iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
      Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
      Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);
      
      /*--- Mean Values ---*/
      
      Mean_ProjVel = node[iPoint]->GetProjVel(Normal);
      Mean_SoundSpeed = node[iPoint]->GetSoundSpeed() * Area;
      
      /*--- Adjustment for grid movement ---*/
      
      if (grid_movement) {
        double *GridVel = geometry->node[iPoint]->GetGridVel();
        ProjVel = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          ProjVel += GridVel[iDim]*Normal[iDim];
        Mean_ProjVel -= ProjVel;
      }
      
      /*--- Inviscid contribution ---*/
      
      Lambda = fabs(Mean_ProjVel) + Mean_SoundSpeed;
      if (geometry->node[iPoint]->GetDomain()) {
        node[iPoint]->AddMax_Lambda_Inv(Lambda);
      }
      
      /*--- Viscous contribution ---*/
      
      Mean_LaminarVisc = node[iPoint]->GetLaminarViscosity();
      Mean_EddyVisc    = node[iPoint]->GetEddyViscosity();
      Mean_Density     = node[iPoint]->GetSolution(0);
      
      Lambda_1 = (4.0/3.0)*(Mean_LaminarVisc + Mean_EddyVisc);
      Lambda_2 = (1.0 + (Prandtl_Lam/Prandtl_Turb)*(Mean_EddyVisc/Mean_LaminarVisc))*(Gamma*Mean_LaminarVisc/Prandtl_Lam);
      Lambda = (Lambda_1 + Lambda_2)*Area*Area/Mean_Density;
      
      if (geometry->node[iPoint]->GetDomain()) node[iPoint]->AddMax_Lambda_Visc(Lambda);
      
    }
  }
  
  /*--- Each element uses their own speed, steady state simulation ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    
    Vol = geometry->node[iPoint]->GetVolume();
    
    if (Vol != 0.0) {
      Local_Delta_Time = config->GetCFL(iMesh)*Vol / node[iPoint]->GetMax_Lambda_Inv();
      Local_Delta_Time_Visc = config->GetCFL(iMesh)*K_v*Vol*Vol/ node[iPoint]->GetMax_Lambda_Visc();
      Local_Delta_Time = min(Local_Delta_Time, Local_Delta_Time_Visc);
      Global_Delta_Time = min(Global_Delta_Time, Local_Delta_Time);
      Min_Delta_Time = min(Min_Delta_Time, Local_Delta_Time);
      Max_Delta_Time = max(Max_Delta_Time, Local_Delta_Time);
      if (Local_Delta_Time > config->GetMax_DeltaTime())
        Local_Delta_Time = config->GetMax_DeltaTime();
      node[iPoint]->SetDelta_Time(Local_Delta_Time);
    }
    else {
      node[iPoint]->SetDelta_Time(0.0);
    }
    
  }
  
  /*--- For exact time solution use the minimum delta time of the whole mesh ---*/
  
  if (config->GetUnsteady_Simulation() == TIME_STEPPING) {
    for(iPoint = 0; iPoint < nPointDomain; iPoint++)
      node[iPoint]->SetDelta_Time(Global_Delta_Time);
  }
  
  /*--- Recompute the unsteady time step for the dual time strategy
   if the unsteady CFL is diferent from 0 ---*/
  
  if ((dual_time) && (Iteration == 0) && (config->GetUnst_CFL() != 0.0) && (iMesh == MESH_0)) {
    Global_Delta_UnstTimeND = config->GetUnst_CFL()*Global_Delta_Time/config->GetCFL(iMesh);
    config->SetDelta_UnstTimeND(Global_Delta_UnstTimeND);
  }
  
  /*--- The pseudo local time (explicit integration) cannot be greater than the physical time ---*/
  
  if (dual_time)
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      if (!implicit) {
        Local_Delta_Time = min((2.0/3.0)*config->GetDelta_UnstTimeND(), node[iPoint]->GetDelta_Time());
        node[iPoint]->SetDelta_Time(Local_Delta_Time);
      }
    }
  
}

void CNSSolver::Viscous_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                 CConfig *config, unsigned short iMesh, unsigned short iRKStep) {
  unsigned long iPoint, jPoint, iEdge;
  
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    /*--- Points, coordinates and normal vector in edge ---*/
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[jPoint]->GetCoord());
    numerics->SetNormal(geometry->edge[iEdge]->GetNormal());
    
    /*--- Primitive variables, and gradient ---*/
    numerics->SetPrimitive(node[iPoint]->GetPrimVar(), node[jPoint]->GetPrimVar());
    numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(), node[jPoint]->GetGradient_Primitive());
    
    /*--- Turbulent kinetic energy ---*/
    if (config->GetKind_Turb_Model() == SST)
      numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->node[iPoint]->GetSolution(0),
                                     solver_container[TURB_SOL]->node[jPoint]->GetSolution(0));
    
    /*--- Compute and update residual ---*/
    numerics->ComputeResidual(Res_Visc, Jacobian_i, Jacobian_j, config);
    
    LinSysRes.SubtractBlock(iPoint, Res_Visc);
    LinSysRes.AddBlock(jPoint, Res_Visc);
    
    /*--- Implicit part ---*/
    if (implicit) {
      Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
      Jacobian.SubtractBlock(iPoint, jPoint, Jacobian_j);
      Jacobian.AddBlock(jPoint, iPoint, Jacobian_i);
      Jacobian.AddBlock(jPoint, jPoint, Jacobian_j);
    }
  }
  
}

void CNSSolver::Viscous_Forces(CGeometry *geometry, CConfig *config) {
  
  unsigned long iVertex, iPoint, iPointNormal;
  unsigned short Boundary, Monitoring, iMarker, iMarker_Monitoring, iDim, jDim;
  double Delta, Viscosity, **Grad_PrimVar, div_vel, *Normal, MomentDist[3], WallDist[3],
  *Coord, *Coord_Normal, Area, WallShearStress, TauNormal, factor, RefVel2,
  RefDensity, GradTemperature, Density, Vel[3], WallDistMod, FrictionVel,
  Mach2Vel, Mach_Motion, *Velocity_Inf, UnitNormal[3], TauElem[3], TauTangent[3], Tau[3][3], Force[3], Cp, thermal_conductivity, MaxNorm = 8.0;
  double *Origin = config->GetRefOriginMoment(0);
  string Marker_Tag, Monitoring_Tag;
  
  double Alpha            = config->GetAoA()*PI_NUMBER/180.0;
  double Beta             = config->GetAoS()*PI_NUMBER/180.0;
  double RefAreaCoeff     = config->GetRefAreaCoeff();
  double RefLengthMoment  = config->GetRefLengthMoment();
  double Gas_Constant     = config->GetGas_ConstantND();
  bool grid_movement      = config->GetGrid_Movement();
  
  /*--- For dynamic meshes, use the motion Mach number as a reference value
   for computing the force coefficients. Otherwise, use the freestream values,
   which is the standard convention. ---*/
  
  if (grid_movement) {
    Mach2Vel = sqrt(Gamma*Gas_Constant*config->GetTemperature_FreeStreamND());
    Mach_Motion = config->GetMach_Motion();
    RefVel2 = (Mach_Motion*Mach2Vel)*(Mach_Motion*Mach2Vel);
  } else {
    Velocity_Inf = config->GetVelocity_FreeStreamND();
    RefVel2 = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      RefVel2  += Velocity_Inf[iDim]*Velocity_Inf[iDim];
  }
  
  RefDensity  = config->GetDensity_FreeStreamND();
  factor = 1.0 / (0.5*RefDensity*RefAreaCoeff*RefVel2);
  
  /*--- Variables initialization ---*/
  
  AllBound_CDrag_Visc = 0.0;  AllBound_CLift_Visc = 0.0;  AllBound_CSideForce_Visc = 0.0;  AllBound_CEff_Visc = 0.0;
  AllBound_CMx_Visc = 0.0;    AllBound_CMy_Visc = 0.0;    AllBound_CMz_Visc = 0.0;
  AllBound_CFx_Visc = 0.0;    AllBound_CFy_Visc = 0.0;    AllBound_CFz_Visc = 0.0;
  AllBound_CT_Visc = 0.0;     AllBound_CQ_Visc = 0.0;     AllBound_CMerit_Visc = 0.0;
  AllBound_HeatFlux_Visc = 0.0;      AllBound_MaxHeatFlux_Visc = 0.0;
  
  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    Surface_CLift_Visc[iMarker_Monitoring] = 0.0;
    Surface_CDrag_Visc[iMarker_Monitoring] = 0.0;
    Surface_CMx_Visc[iMarker_Monitoring]   = 0.0;
    Surface_CMy_Visc[iMarker_Monitoring]   = 0.0;
    Surface_CMz_Visc[iMarker_Monitoring]   = 0.0;
  }
  
  /*--- Loop over the Navier-Stokes markers ---*/
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    Boundary = config->GetMarker_All_KindBC(iMarker);
    Monitoring = config->GetMarker_All_Monitoring(iMarker);
    
    /*--- Obtain the origin for the moment computation for a particular marker ---*/
    if (Monitoring == YES) {
      for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
        Monitoring_Tag = config->GetMarker_Monitoring(iMarker_Monitoring);
        Marker_Tag = config->GetMarker_All_TagBound(iMarker);
        if (Marker_Tag == Monitoring_Tag)
          Origin = config->GetRefOriginMoment(iMarker_Monitoring);
      }
    }
    
    if ((Boundary == HEAT_FLUX) || (Boundary == ISOTHERMAL)) {
      
      /*--- Forces initialization at each Marker ---*/
      CDrag_Visc[iMarker] = 0.0;  CLift_Visc[iMarker] = 0.0; CSideForce_Visc[iMarker] = 0.0;  CEff_Visc[iMarker] = 0.0;
      CMx_Visc[iMarker] = 0.0;    CMy_Visc[iMarker] = 0.0;   CMz_Visc[iMarker] = 0.0;
      CFx_Visc[iMarker] = 0.0;    CFy_Visc[iMarker] = 0.0;   CFz_Visc[iMarker] = 0.0;
      CT_Visc[iMarker] = 0.0;     CQ_Visc[iMarker] = 0.0;    CMerit_Visc[iMarker] = 0.0;
      Heat_Visc[iMarker] = 0.0;      MaxHeatFlux_Visc[iMarker] = 0.0;
      
      for (iDim = 0; iDim < nDim; iDim++) ForceViscous[iDim] = 0.0;
      MomentViscous[0] = 0.0; MomentViscous[1] = 0.0; MomentViscous[2] = 0.0;
      
      /*--- Loop over the vertices to compute the forces ---*/
      
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        iPointNormal = geometry->vertex[iMarker][iVertex]->GetNormal_Neighbor();
        
        Coord = geometry->node[iPoint]->GetCoord();
        Coord_Normal = geometry->node[iPointNormal]->GetCoord();
        
        Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
        Grad_PrimVar = node[iPoint]->GetGradient_Primitive();
        Viscosity = node[iPoint]->GetLaminarViscosity();
        Density = node[iPoint]->GetDensity();
        
        Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);
        for (iDim = 0; iDim < nDim; iDim++) {
          UnitNormal[iDim] = Normal[iDim]/Area;
          MomentDist[iDim] = Coord[iDim] - Origin[iDim];
        }
        
        div_vel = 0.0; for (iDim = 0; iDim < nDim; iDim++) div_vel += Grad_PrimVar[iDim+1][iDim];
        
        for (iDim = 0; iDim < nDim; iDim++) {
          for (jDim = 0 ; jDim < nDim; jDim++) {
            Delta = 0.0; if (iDim == jDim) Delta = 1.0;
            Tau[iDim][jDim] = Viscosity*(Grad_PrimVar[jDim+1][iDim] + Grad_PrimVar[iDim+1][jDim]) -
            TWO3*Viscosity*div_vel*Delta;
          }
          TauElem[iDim] = 0.0;
          for (jDim = 0; jDim < nDim; jDim++)
            TauElem[iDim] += Tau[iDim][jDim]*UnitNormal[jDim];
        }
        
        /*--- Compute wall shear stress (using the stress tensor) ---*/
        
        TauNormal = 0.0; for (iDim = 0; iDim < nDim; iDim++) TauNormal += TauElem[iDim] * UnitNormal[iDim];
        for (iDim = 0; iDim < nDim; iDim++) TauTangent[iDim] = TauElem[iDim] - TauNormal * UnitNormal[iDim];
        WallShearStress = 0.0; for (iDim = 0; iDim < nDim; iDim++) WallShearStress += TauTangent[iDim]*TauTangent[iDim];
        WallShearStress = sqrt(WallShearStress);
        
        for (iDim = 0; iDim < nDim; iDim++)
          Vel[iDim] = node[iPointNormal]->GetVelocity(iDim);
        
        for (iDim = 0; iDim < nDim; iDim++) WallDist[iDim] = (Coord[iDim] - Coord_Normal[iDim]);
        WallDistMod = 0.0; for (iDim = 0; iDim < nDim; iDim++) WallDistMod += WallDist[iDim]*WallDist[iDim]; WallDistMod = sqrt(WallDistMod);
        
        /*--- Compute wall skin friction coefficient, and heat flux on the wall ---*/
        
        CSkinFriction[iMarker][iVertex] = WallShearStress / (0.5*RefDensity*RefVel2);
        
        /*--- Compute y+ and non-dimensional velocity ---*/
        
        FrictionVel = sqrt(fabs(WallShearStress)/Density);
        YPlus[iMarker][iVertex] = WallDistMod*FrictionVel/(Viscosity/Density);
        
        /*--- Compute total and max heat flux on the wall (compressible solver only) ---*/
        
        GradTemperature = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          GradTemperature += Grad_PrimVar[0][iDim]*(-Normal[iDim]);
        
        Cp = (Gamma / Gamma_Minus_One) * Gas_Constant;
        thermal_conductivity = Cp * Viscosity/PRANDTL;
        HeatFlux[iMarker][iVertex] = -thermal_conductivity*GradTemperature;
        Heat_Visc[iMarker] += HeatFlux[iMarker][iVertex]*Area;
        MaxHeatFlux_Visc[iMarker] += pow(HeatFlux[iMarker][iVertex], MaxNorm);
        
        /*--- Note that y+, and heat are computed at the
         halo cells (for visualization purposes), but not the forces ---*/
        
        if ((geometry->node[iPoint]->GetDomain()) && (Monitoring == YES)) {
          
          /*--- Force computation ---*/
          
          for (iDim = 0; iDim < nDim; iDim++) {
            Force[iDim] = TauElem[iDim]*Area*factor;
            ForceViscous[iDim] += Force[iDim];
          }
          
          /*--- Moment with respect to the reference axis ---*/
          
          if (iDim == 3) {
            MomentViscous[0] += (Force[2]*MomentDist[1] - Force[1]*MomentDist[2])/RefLengthMoment;
            MomentViscous[1] += (Force[0]*MomentDist[2] - Force[2]*MomentDist[0])/RefLengthMoment;
          }
          MomentViscous[2] += (Force[1]*MomentDist[0] - Force[0]*MomentDist[1])/RefLengthMoment;
          
        }
        
      }
      
      /*--- Transform ForceViscous and MomentViscous into non-dimensional coefficient ---*/
      
      if  (Monitoring == YES) {
        if (nDim == 2) {
          CDrag_Visc[iMarker]       =  ForceViscous[0]*cos(Alpha) + ForceViscous[1]*sin(Alpha);
          CLift_Visc[iMarker]       = -ForceViscous[0]*sin(Alpha) + ForceViscous[1]*cos(Alpha);
          CEff_Visc[iMarker]        = CLift_Visc[iMarker]/(CDrag_Visc[iMarker]+EPS);
          CMz_Visc[iMarker]         = MomentViscous[2];
          CFx_Visc[iMarker]         = ForceViscous[0];
          CFy_Visc[iMarker]         = ForceViscous[1];
          CT_Visc[iMarker]          = -CFx_Visc[iMarker];
          CQ_Visc[iMarker]          = -CMz_Visc[iMarker];
          CMerit_Visc[iMarker]      = CT_Visc[iMarker]/CQ_Visc[iMarker];
          MaxHeatFlux_Visc[iMarker] = pow(MaxHeatFlux_Visc[iMarker], 1.0/MaxNorm);
        }
        if (nDim == 3) {
          CDrag_Visc[iMarker]       =  ForceViscous[0]*cos(Alpha)*cos(Beta) + ForceViscous[1]*sin(Beta) + ForceViscous[2]*sin(Alpha)*cos(Beta);
          CLift_Visc[iMarker]       = -ForceViscous[0]*sin(Alpha) + ForceViscous[2]*cos(Alpha);
          CSideForce_Visc[iMarker]  = -ForceViscous[0]*sin(Beta)*cos(Alpha) + ForceViscous[1]*cos(Beta) - ForceViscous[2]*sin(Beta)*sin(Alpha);
          CEff_Visc[iMarker]        = CLift_Visc[iMarker]/(CDrag_Visc[iMarker]+EPS);
          CMx_Visc[iMarker]         = MomentViscous[0];
          CMy_Visc[iMarker]         = MomentViscous[1];
          CMz_Visc[iMarker]         = MomentViscous[2];
          CFx_Visc[iMarker]         = ForceViscous[0];
          CFy_Visc[iMarker]         = ForceViscous[1];
          CFz_Visc[iMarker]         = ForceViscous[2];
          CT_Visc[iMarker]          = -CFz_Visc[iMarker];
          CQ_Visc[iMarker]          = -CMz_Visc[iMarker];
          CMerit_Visc[iMarker]      = CT_Visc[iMarker]/CQ_Visc[iMarker];
          MaxHeatFlux_Visc[iMarker] = pow(MaxHeatFlux_Visc[iMarker], 1.0/MaxNorm);
        }
        
        AllBound_CDrag_Visc       += CDrag_Visc[iMarker];
        AllBound_CLift_Visc       += CLift_Visc[iMarker];
        AllBound_CSideForce_Visc  += CSideForce_Visc[iMarker];
        AllBound_CMx_Visc         += CMx_Visc[iMarker];
        AllBound_CMy_Visc         += CMy_Visc[iMarker];
        AllBound_CMz_Visc         += CMz_Visc[iMarker];
        AllBound_CFx_Visc         += CFx_Visc[iMarker];
        AllBound_CFy_Visc         += CFy_Visc[iMarker];
        AllBound_CFz_Visc         += CFz_Visc[iMarker];
        AllBound_CT_Visc          += CT_Visc[iMarker];
        AllBound_CQ_Visc          += CQ_Visc[iMarker];
        AllBound_HeatFlux_Visc    += Heat_Visc[iMarker];
        AllBound_MaxHeatFlux_Visc += pow(MaxHeatFlux_Visc[iMarker], MaxNorm);
        
        /*--- Compute the coefficients per surface ---*/
        
        for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
          Monitoring_Tag = config->GetMarker_Monitoring(iMarker_Monitoring);
          Marker_Tag = config->GetMarker_All_TagBound(iMarker);
          if (Marker_Tag == Monitoring_Tag) {
            Surface_CLift_Visc[iMarker_Monitoring] += CLift_Visc[iMarker];
            Surface_CDrag_Visc[iMarker_Monitoring] += CDrag_Visc[iMarker];
            Surface_CMx_Visc[iMarker_Monitoring]   += CMx_Visc[iMarker];
            Surface_CMy_Visc[iMarker_Monitoring]   += CMy_Visc[iMarker];
            Surface_CMz_Visc[iMarker_Monitoring]   += CMz_Visc[iMarker];
          }
        }
        
      }
      
    }
  }
  
  /*--- Update some global coeffients ---*/
  
  AllBound_CEff_Visc = AllBound_CLift_Visc / (AllBound_CDrag_Visc + EPS);
  AllBound_CMerit_Visc = AllBound_CT_Visc / (AllBound_CQ_Visc + EPS);
  AllBound_MaxHeatFlux_Visc = pow(AllBound_MaxHeatFlux_Visc, 1.0/MaxNorm);
  
  /*--- Update the total coefficients (note that all the nodes have the same value)---*/
  
  Total_CDrag       += AllBound_CDrag_Visc;
  Total_CLift       += AllBound_CLift_Visc;
  Total_CSideForce  += AllBound_CSideForce_Visc;
  Total_CEff        = Total_CLift / (Total_CDrag + EPS);
  Total_CMx         += AllBound_CMx_Visc;
  Total_CMy         += AllBound_CMy_Visc;
  Total_CMz         += AllBound_CMz_Visc;
  Total_CFx         += AllBound_CFx_Visc;
  Total_CFy         += AllBound_CFy_Visc;
  Total_CFz         += AllBound_CFz_Visc;
  Total_CT          += AllBound_CT_Visc;
  Total_CQ          += AllBound_CQ_Visc;
  Total_CMerit      += AllBound_CMerit_Visc;
  Total_Heat        = AllBound_HeatFlux_Visc;
  Total_MaxHeat     = AllBound_MaxHeatFlux_Visc;
  
  /*--- Update the total coefficients per surface (note that all the nodes have the same value)---*/
  
  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    Surface_CLift[iMarker_Monitoring]     += Surface_CLift_Visc[iMarker_Monitoring];
    Surface_CDrag[iMarker_Monitoring]     += Surface_CDrag_Visc[iMarker_Monitoring];
    Surface_CMx[iMarker_Monitoring]       += Surface_CMx_Visc[iMarker_Monitoring];
    Surface_CMy[iMarker_Monitoring]       += Surface_CMy_Visc[iMarker_Monitoring];
    Surface_CMz[iMarker_Monitoring]       += Surface_CMz_Visc[iMarker_Monitoring];
  }
  
}

void CNSSolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  /*--- Local variables ---*/
  unsigned short iDim, jDim, iVar, jVar;
  unsigned long iVertex, iPoint, Point_Normal, total_index;
  
  double Wall_HeatFlux, dist_ij, *Coord_i, *Coord_j, theta2;
  double thetax, thetay, thetaz, etax, etay, etaz, pix, piy, piz, factor;
  double ProjGridVel, *GridVel, GridVel2, *Normal, Area, Pressure;
  double total_viscosity, div_vel, Density, turb_ke, tau_vel[3], UnitNormal[3];
  double laminar_viscosity, eddy_viscosity, **grad_primvar, tau[3][3];
  double delta[3][3] = {{1.0, 0.0, 0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}};
  
  bool implicit       = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool grid_movement  = config->GetGrid_Movement();
  
  /*--- Identify the boundary by string name ---*/
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  
  /*--- Get the specified wall heat flux from config ---*/
  Wall_HeatFlux = config->GetWall_HeatFlux(Marker_Tag);
  
  /*--- Loop over all of the vertices on this boundary marker ---*/
  for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Compute dual-grid area and boundary normal ---*/
      Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
      
      Area = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        Area += Normal[iDim]*Normal[iDim];
      Area = sqrt (Area);
      
      for (iDim = 0; iDim < nDim; iDim++)
        UnitNormal[iDim] = -Normal[iDim]/Area;
      
      /*--- Initialize the convective & viscous residuals to zero ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        Res_Conv[iVar] = 0.0;
        Res_Visc[iVar] = 0.0;
      }
      
      /*--- Store the corrected velocity at the wall which will
       be zero (v = 0), unless there are moving walls (v = u_wall)---*/
      if (grid_movement) {
        GridVel = geometry->node[iPoint]->GetGridVel();
        for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = GridVel[iDim];
      } else {
        for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = 0.0;
      }
      
      /*--- Impose the value of the velocity as a strong boundary
       condition (Dirichlet). Fix the velocity and remove any
       contribution to the residual at this node. ---*/
      node[iPoint]->SetVelocity_Old(Vector);
      
      for (iDim = 0; iDim < nDim; iDim++)
        LinSysRes.SetBlock_Zero(iPoint, iDim+1);
      node[iPoint]->SetVel_ResTruncError_Zero();
      
      /*--- Apply a weak boundary condition for the energy equation.
       Compute the residual due to the prescribed heat flux. ---*/
      Res_Visc[nDim+1] = Wall_HeatFlux * Area;
      
      /*--- If the wall is moving, there are additional residual contributions
       due to pressure (p v_wall.n) and shear stress (tau.v_wall.n). ---*/
      if (grid_movement) {
        
        /*--- Get the grid velocity at the current boundary node ---*/
        GridVel = geometry->node[iPoint]->GetGridVel();
        ProjGridVel = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          ProjGridVel += GridVel[iDim]*UnitNormal[iDim]*Area;
        
        /*--- Retrieve other primitive quantities and viscosities ---*/
        Density  = node[iPoint]->GetSolution(0);
        Pressure = node[iPoint]->GetPressure();
        laminar_viscosity = node[iPoint]->GetLaminarViscosity();
        eddy_viscosity    = node[iPoint]->GetEddyViscosity();
        total_viscosity   = laminar_viscosity + eddy_viscosity;
        grad_primvar      = node[iPoint]->GetGradient_Primitive();
        
        /*--- Turbulent kinetic energy ---*/
        if (config->GetKind_Turb_Model() == SST)
          turb_ke = solver_container[TURB_SOL]->node[iPoint]->GetSolution(0);
        else
          turb_ke = 0.0;
        
        /*--- Divergence of the velocity ---*/
        div_vel = 0.0;
        for (iDim = 0 ; iDim < nDim; iDim++)
          div_vel += grad_primvar[iDim+1][iDim];
        
        /*--- Compute the viscous stress tensor ---*/
        for (iDim = 0; iDim < nDim; iDim++)
          for (jDim = 0; jDim < nDim; jDim++) {
            tau[iDim][jDim] = total_viscosity*( grad_primvar[jDim+1][iDim]
                                               +grad_primvar[iDim+1][jDim] )
            - TWO3*total_viscosity*div_vel*delta[iDim][jDim]
            - TWO3*Density*turb_ke*delta[iDim][jDim];
          }
        
        /*--- Dot product of the stress tensor with the grid velocity ---*/
        for (iDim = 0 ; iDim < nDim; iDim++) {
          tau_vel[iDim] = 0.0;
          for (jDim = 0 ; jDim < nDim; jDim++)
            tau_vel[iDim] += tau[iDim][jDim]*GridVel[jDim];
        }
        
        /*--- Compute the convective and viscous residuals (energy eqn.) ---*/
        Res_Conv[nDim+1] = Pressure*ProjGridVel;
        for (iDim = 0 ; iDim < nDim; iDim++)
          Res_Visc[nDim+1] += tau_vel[iDim]*UnitNormal[iDim]*Area;
        
        /*--- Implicit Jacobian contributions due to moving walls ---*/
        if (implicit) {
          
          /*--- Jacobian contribution related to the pressure term ---*/
          GridVel2 = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            GridVel2 += GridVel[iDim]*GridVel[iDim];
          for (iVar = 0; iVar < nVar; iVar++)
            for (jVar = 0; jVar < nVar; jVar++)
              Jacobian_i[iVar][jVar] = 0.0;
          Jacobian_i[nDim+1][0] = 0.5*(Gamma-1.0)*GridVel2*ProjGridVel;
          for (jDim = 0; jDim < nDim; jDim++)
            Jacobian_i[nDim+1][jDim+1] = -(Gamma-1.0)*GridVel[jDim]*ProjGridVel;
          Jacobian_i[nDim+1][nDim+1] = (Gamma-1.0)*ProjGridVel;
          
          /*--- Add the block to the Global Jacobian structure ---*/
          Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
          
          /*--- Now the Jacobian contribution related to the shear stress ---*/
          for (iVar = 0; iVar < nVar; iVar++)
            for (jVar = 0; jVar < nVar; jVar++)
              Jacobian_i[iVar][jVar] = 0.0;
          
          /*--- Compute closest normal neighbor ---*/
          Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
          
          /*--- Get coordinates of i & nearest normal and compute distance ---*/
          Coord_i = geometry->node[iPoint]->GetCoord();
          Coord_j = geometry->node[Point_Normal]->GetCoord();
          dist_ij = 0;
          for (iDim = 0; iDim < nDim; iDim++)
            dist_ij += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
          dist_ij = sqrt(dist_ij);
          
          theta2 = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            theta2 += UnitNormal[iDim]*UnitNormal[iDim];
          
          factor = total_viscosity*Area/(Density*dist_ij);
          
          if (nDim == 2) {
            thetax = theta2 + UnitNormal[0]*UnitNormal[0]/3.0;
            thetay = theta2 + UnitNormal[1]*UnitNormal[1]/3.0;
            
            etaz   = UnitNormal[0]*UnitNormal[1]/3.0;
            
            pix = GridVel[0]*thetax + GridVel[1]*etaz;
            piy = GridVel[0]*etaz   + GridVel[1]*thetay;
            
            Jacobian_i[nDim+1][0] -= factor*(-pix*GridVel[0]+piy*GridVel[1]);
            Jacobian_i[nDim+1][1] -= factor*pix;
            Jacobian_i[nDim+1][2] -= factor*piy;
          } else {
            thetax = theta2 + UnitNormal[0]*UnitNormal[0]/3.0;
            thetay = theta2 + UnitNormal[1]*UnitNormal[1]/3.0;
            thetaz = theta2 + UnitNormal[2]*UnitNormal[2]/3.0;
            
            etaz = UnitNormal[0]*UnitNormal[1]/3.0;
            etax = UnitNormal[1]*UnitNormal[2]/3.0;
            etay = UnitNormal[0]*UnitNormal[2]/3.0;
            
            pix = GridVel[0]*thetax + GridVel[1]*etaz   + GridVel[2]*etay;
            piy = GridVel[0]*etaz   + GridVel[1]*thetay + GridVel[2]*etax;
            piz = GridVel[0]*etay   + GridVel[1]*etax   + GridVel[2]*thetaz;
            
            Jacobian_i[nDim+1][0] -= factor*(-pix*GridVel[0]+piy*GridVel[1]+piz*GridVel[2]);
            Jacobian_i[nDim+1][1] -= factor*pix;
            Jacobian_i[nDim+1][2] -= factor*piy;
            Jacobian_i[nDim+1][3] -= factor*piz;
          }
          
          /*--- Subtract the block from the Global Jacobian structure ---*/
          Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
        }
      }
      
      /*--- Convective contribution to the residual at the wall ---*/
      LinSysRes.AddBlock(iPoint, Res_Conv);
      
      /*--- Viscous contribution to the residual at the wall ---*/
      LinSysRes.SubtractBlock(iPoint, Res_Visc);
      
      /*--- Enforce the no-slip boundary condition in a strong way by
       modifying the velocity-rows of the Jacobian (1 on the diagonal). ---*/
      if (implicit) {
        for (iVar = 1; iVar <= nDim; iVar++) {
          total_index = iPoint*nVar+iVar;
          Jacobian.DeleteValsRowi(total_index);
        }
      }
      
    }
  }
}

void CNSSolver::BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  unsigned short iVar, jVar, iDim, jDim;
  unsigned long iVertex, iPoint, Point_Normal, total_index;
  
  double *Normal, *Coord_i, *Coord_j, Area, dist_ij, theta2;
  double Twall, Temperature, dTdn, dTdrho, thermal_conductivity;
  double thetax, thetay, thetaz, etax, etay, etaz, pix, piy, piz, factor;
  double ProjGridVel, *GridVel, GridVel2, Pressure, Density, Vel2, Energy;
  double total_viscosity, div_vel, turb_ke, tau_vel[3], UnitNormal[3];
  double laminar_viscosity, eddy_viscosity, **grad_primvar, tau[3][3];
  double delta[3][3] = {{1.0, 0.0, 0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}};
  
  double Prandtl_Lam  = config->GetPrandtl_Lam();
  double Prandtl_Turb = config->GetPrandtl_Turb();
  double Gas_Constant = config->GetGas_ConstantND();
  double cp = (Gamma / Gamma_Minus_One) * Gas_Constant;
  
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool grid_movement  = config->GetGrid_Movement();
  
  Point_Normal = 0;
  
  /*--- Identify the boundary ---*/
  
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  
  /*--- Retrieve the specified wall temperature ---*/
  
  Twall = config->GetIsothermal_Temperature(Marker_Tag);
  
  /*--- Loop over boundary points ---*/
  
  for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Compute dual-grid area and boundary normal ---*/
      
      Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
      Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt (Area);
      for (iDim = 0; iDim < nDim; iDim++)
        UnitNormal[iDim] = -Normal[iDim]/Area;
      
      /*--- Calculate useful quantities ---*/
      
      theta2 = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        theta2 += UnitNormal[iDim]*UnitNormal[iDim];
      
      /*--- Compute closest normal neighbor ---*/
      
      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
      
      /*--- Get coordinates of i & nearest normal and compute distance ---*/
      
      Coord_i = geometry->node[iPoint]->GetCoord();
      Coord_j = geometry->node[Point_Normal]->GetCoord();
      dist_ij = 0;
      for (iDim = 0; iDim < nDim; iDim++)
        dist_ij += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
      dist_ij = sqrt(dist_ij);
      
      /*--- Store the corrected velocity at the wall which will
       be zero (v = 0), unless there is grid motion (v = u_wall)---*/
      
      if (grid_movement) {
        GridVel = geometry->node[iPoint]->GetGridVel();
        for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = GridVel[iDim];
      }
      else {
        for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = 0.0;
      }
      
      /*--- Initialize the convective & viscous residuals to zero ---*/
      
      for (iVar = 0; iVar < nVar; iVar++) {
        Res_Conv[iVar] = 0.0;
        Res_Visc[iVar] = 0.0;
      }
      
      /*--- Set the residual, truncation error and velocity value on the boundary ---*/
      
      node[iPoint]->SetVelocity_Old(Vector);
      
      for (iDim = 0; iDim < nDim; iDim++)
        LinSysRes.SetBlock_Zero(iPoint, iDim+1);
      node[iPoint]->SetVel_ResTruncError_Zero();
      
      /*--- Compute the normal gradient in temperature using Twall ---*/
      
      dTdn = -(node[Point_Normal]->GetPrimVar(0) - Twall)/dist_ij;
      
      /*--- Get transport coefficients ---*/
      
      laminar_viscosity    = node[iPoint]->GetLaminarViscosity();
      eddy_viscosity       = node[iPoint]->GetEddyViscosity();
      thermal_conductivity = cp * ( laminar_viscosity/Prandtl_Lam + eddy_viscosity/Prandtl_Turb);
      
      /*--- Apply a weak boundary condition for the energy equation.
       Compute the residual due to the prescribed heat flux. ---*/
      
      Res_Visc[nDim+1] = thermal_conductivity * dTdn * Area;
      
      /*--- Calculate Jacobian for implicit time stepping ---*/
      
      if (implicit) {
        
        for (iVar = 0; iVar < nVar; iVar ++)
          for (jVar = 0; jVar < nVar; jVar ++)
            Jacobian_i[iVar][jVar] = 0.0;
        
        /*--- Calculate useful quantities ---*/
        
        Density = node[iPoint]->GetPrimVar(nDim+2);
        Energy  = node[iPoint]->GetSolution(nDim+1);
        Temperature = node[iPoint]->GetPrimVar(0);
        Vel2 = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          Vel2 += node[iPoint]->GetPrimVar(iDim+1) * node[iPoint]->GetPrimVar(iDim+1);
        dTdrho = 1.0/Density * ( -Twall + (Gamma-1.0)/Gas_Constant*(Vel2/2.0) );
        
        /*--- Enforce the no-slip boundary condition in a strong way ---*/
        
        for (iVar = 1; iVar <= nDim; iVar++) {
          total_index = iPoint*nVar+iVar;
          Jacobian.DeleteValsRowi(total_index);
        }
        
        /*--- Add contributions to the Jacobian from the weak enforcement of the energy equations ---*/
        
        Jacobian_i[nDim+1][0]      = -thermal_conductivity*theta2/dist_ij * dTdrho * Area;
        Jacobian_i[nDim+1][nDim+1] = -thermal_conductivity*theta2/dist_ij * (Gamma-1.0)/(Gas_Constant*Density) * Area;
        
        /*--- Subtract the block from the Global Jacobian structure ---*/
        
        Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
        
      }
      
      /*--- If the wall is moving, there are additional residual contributions
       due to pressure (p v_wall.n) and shear stress (tau.v_wall.n). ---*/
      
      if (grid_movement) {
        
        /*--- Get the grid velocity at the current boundary node ---*/
        
        GridVel = geometry->node[iPoint]->GetGridVel();
        ProjGridVel = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          ProjGridVel += GridVel[iDim]*UnitNormal[iDim]*Area;
        
        /*--- Retrieve other primitive quantities and viscosities ---*/
        
        Density  = node[iPoint]->GetSolution(0);
        Pressure = node[iPoint]->GetPressure();
        laminar_viscosity = node[iPoint]->GetLaminarViscosity();
        eddy_viscosity    = node[iPoint]->GetEddyViscosity();
        
        total_viscosity   = laminar_viscosity + eddy_viscosity;
        grad_primvar      = node[iPoint]->GetGradient_Primitive();
        
        /*--- Turbulent kinetic energy ---*/
        
        if (config->GetKind_Turb_Model() == SST)
          turb_ke = solver_container[TURB_SOL]->node[iPoint]->GetSolution(0);
        else
          turb_ke = 0.0;
        
        /*--- Divergence of the velocity ---*/
        
        div_vel = 0.0;
        for (iDim = 0 ; iDim < nDim; iDim++)
          div_vel += grad_primvar[iDim+1][iDim];
        
        /*--- Compute the viscous stress tensor ---*/
        
        for (iDim = 0; iDim < nDim; iDim++)
          for (jDim = 0; jDim < nDim; jDim++) {
            tau[iDim][jDim] = total_viscosity*( grad_primvar[jDim+1][iDim]
                                               +grad_primvar[iDim+1][jDim] )
            - TWO3*total_viscosity*div_vel*delta[iDim][jDim]
            - TWO3*Density*turb_ke*delta[iDim][jDim];
          }
        
        /*--- Dot product of the stress tensor with the grid velocity ---*/
        
        for (iDim = 0 ; iDim < nDim; iDim++) {
          tau_vel[iDim] = 0.0;
          for (jDim = 0 ; jDim < nDim; jDim++)
            tau_vel[iDim] += tau[iDim][jDim]*GridVel[jDim];
        }
        
        /*--- Compute the convective and viscous residuals (energy eqn.) ---*/
        
        Res_Conv[nDim+1] = Pressure*ProjGridVel;
        for (iDim = 0 ; iDim < nDim; iDim++)
          Res_Visc[nDim+1] += tau_vel[iDim]*UnitNormal[iDim]*Area;
        
        /*--- Implicit Jacobian contributions due to moving walls ---*/
        
        if (implicit) {
          
          /*--- Jacobian contribution related to the pressure term ---*/
          
          GridVel2 = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            GridVel2 += GridVel[iDim]*GridVel[iDim];
          for (iVar = 0; iVar < nVar; iVar++)
            for (jVar = 0; jVar < nVar; jVar++)
              Jacobian_i[iVar][jVar] = 0.0;
          
          Jacobian_i[nDim+1][0] = 0.5*(Gamma-1.0)*GridVel2*ProjGridVel;
          for (jDim = 0; jDim < nDim; jDim++)
            Jacobian_i[nDim+1][jDim+1] = -(Gamma-1.0)*GridVel[jDim]*ProjGridVel;
          Jacobian_i[nDim+1][nDim+1] = (Gamma-1.0)*ProjGridVel;
          
          /*--- Add the block to the Global Jacobian structure ---*/
          
          Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
          
          /*--- Now the Jacobian contribution related to the shear stress ---*/
          
          for (iVar = 0; iVar < nVar; iVar++)
            for (jVar = 0; jVar < nVar; jVar++)
              Jacobian_i[iVar][jVar] = 0.0;
          
          factor = total_viscosity*Area/(Density*dist_ij);
          
          if (nDim == 2) {
            thetax = theta2 + UnitNormal[0]*UnitNormal[0]/3.0;
            thetay = theta2 + UnitNormal[1]*UnitNormal[1]/3.0;
            
            etaz   = UnitNormal[0]*UnitNormal[1]/3.0;
            
            pix = GridVel[0]*thetax + GridVel[1]*etaz;
            piy = GridVel[0]*etaz   + GridVel[1]*thetay;
            
            Jacobian_i[nDim+1][0] -= factor*(-pix*GridVel[0]+piy*GridVel[1]);
            Jacobian_i[nDim+1][1] -= factor*pix;
            Jacobian_i[nDim+1][2] -= factor*piy;
          }
          else {
            thetax = theta2 + UnitNormal[0]*UnitNormal[0]/3.0;
            thetay = theta2 + UnitNormal[1]*UnitNormal[1]/3.0;
            thetaz = theta2 + UnitNormal[2]*UnitNormal[2]/3.0;
            
            etaz = UnitNormal[0]*UnitNormal[1]/3.0;
            etax = UnitNormal[1]*UnitNormal[2]/3.0;
            etay = UnitNormal[0]*UnitNormal[2]/3.0;
            
            pix = GridVel[0]*thetax + GridVel[1]*etaz   + GridVel[2]*etay;
            piy = GridVel[0]*etaz   + GridVel[1]*thetay + GridVel[2]*etax;
            piz = GridVel[0]*etay   + GridVel[1]*etax   + GridVel[2]*thetaz;
            
            Jacobian_i[nDim+1][0] -= factor*(-pix*GridVel[0]+piy*GridVel[1]+piz*GridVel[2]);
            Jacobian_i[nDim+1][1] -= factor*pix;
            Jacobian_i[nDim+1][2] -= factor*piy;
            Jacobian_i[nDim+1][3] -= factor*piz;
          }
          
          /*--- Subtract the block from the Global Jacobian structure ---*/
          
          Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
        }
        
      }
      
      /*--- Convective contribution to the residual at the wall ---*/
      
      LinSysRes.AddBlock(iPoint, Res_Conv);
      
      /*--- Viscous contribution to the residual at the wall ---*/
      
      LinSysRes.SubtractBlock(iPoint, Res_Visc);
      
      /*--- Enforce the no-slip boundary condition in a strong way by
       modifying the velocity-rows of the Jacobian (1 on the diagonal). ---*/
      
      if (implicit) {
        for (iVar = 1; iVar <= nDim; iVar++) {
          total_index = iPoint*nVar+iVar;
          Jacobian.DeleteValsRowi(total_index);
        }
      }
      
    }
  }
}
