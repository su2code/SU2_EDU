/*!
 * \file solution_direct_turbulent.cpp
 * \brief Main subrotuines for solving direct problems (Euler, Navier-Stokes, etc.).
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

#include "../include/solver_structure.hpp"

CTurbSolver::CTurbSolver(void) : CSolver() {
  
  /*--- Pointer initialization ---*/
  
  FlowPrimVar_i = NULL;
  FlowPrimVar_j = NULL;
  lowerlimit    = NULL;
  upperlimit    = NULL;
  
}

CTurbSolver::CTurbSolver(CConfig *config) : CSolver() {
  
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
  /*--- Pointer initialization ---*/
  
  FlowPrimVar_i = NULL;
  FlowPrimVar_j = NULL;
  lowerlimit    = NULL;
  upperlimit    = NULL;
  
}

CTurbSolver::~CTurbSolver(void) {
  
  if (FlowPrimVar_i != NULL) delete [] FlowPrimVar_i;
  if (FlowPrimVar_j != NULL) delete [] FlowPrimVar_j;
  if (lowerlimit    != NULL) delete [] lowerlimit;
  if (upperlimit    != NULL) delete [] upperlimit;
  
}

void CTurbSolver::Upwind_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, unsigned short iMesh) {
  double *Turb_i, *Turb_j, *Limiter_i = NULL, *Limiter_j = NULL, *V_i, *V_j, **Gradient_i, **Gradient_j, Project_Grad_i, Project_Grad_j;
  unsigned long iEdge, iPoint, jPoint;
  unsigned short iDim, iVar;
  
  bool high_order_diss = (config->GetKind_Upwind_Turb() == SCALAR_UPWIND_2ND);
  bool limiter = (config->GetKind_SlopeLimit() != NONE);
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    /*--- Points in edge and normal vectors ---*/
    
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    numerics->SetNormal(geometry->edge[iEdge]->GetNormal());
    
    /*--- Conservative variables w/o reconstruction ---*/
    
    V_i = solver_container[FLOW_SOL]->node[iPoint]->GetPrimVar();
    V_j = solver_container[FLOW_SOL]->node[jPoint]->GetPrimVar();
    numerics->SetPrimitive(V_i, V_j);
    
    /*--- Turbulent variables w/o reconstruction ---*/
    
    Turb_i = node[iPoint]->GetSolution();
    Turb_j = node[jPoint]->GetSolution();
    numerics->SetTurbVar(Turb_i,Turb_j);
    
    if (high_order_diss) {
      
      for (iDim = 0; iDim < nDim; iDim++) {
        Vector_i[iDim] = 0.5*(geometry->node[jPoint]->GetCoord(iDim) - geometry->node[iPoint]->GetCoord(iDim));
        Vector_j[iDim] = 0.5*(geometry->node[iPoint]->GetCoord(iDim) - geometry->node[jPoint]->GetCoord(iDim));
      }
      
      /*--- Mean flow primitive variables using gradient reconstruction and limiters ---*/
      
      Gradient_i = solver_container[FLOW_SOL]->node[iPoint]->GetGradient_Primitive();
      Gradient_j = solver_container[FLOW_SOL]->node[jPoint]->GetGradient_Primitive();
      if (limiter) {
        Limiter_i = solver_container[FLOW_SOL]->node[iPoint]->GetLimiter_Primitive();
        Limiter_j = solver_container[FLOW_SOL]->node[jPoint]->GetLimiter_Primitive();
      }
      
      for (iVar = 0; iVar < solver_container[FLOW_SOL]->GetnPrimVarGrad(); iVar++) {
        Project_Grad_i = 0.0; Project_Grad_j = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          Project_Grad_i += Vector_i[iDim]*Gradient_i[iVar][iDim];
          Project_Grad_j += Vector_j[iDim]*Gradient_j[iVar][iDim];
        }
        if (limiter) {
          FlowPrimVar_i[iVar] = V_i[iVar] + Limiter_i[iVar]*Project_Grad_i;
          FlowPrimVar_j[iVar] = V_j[iVar] + Limiter_j[iVar]*Project_Grad_j;
        }
        else {
          FlowPrimVar_i[iVar] = V_i[iVar] + Project_Grad_i;
          FlowPrimVar_j[iVar] = V_j[iVar] + Project_Grad_j;
        }
      }
      
      numerics->SetPrimitive(FlowPrimVar_i, FlowPrimVar_j);
      
      /*--- Turbulent variables using gradient reconstruction and limiters ---*/
      
      Gradient_i = node[iPoint]->GetGradient();
      Gradient_j = node[jPoint]->GetGradient();
      //      if (limiter) {
      //        Limiter_i = node[iPoint]->GetLimiter();
      //        Limiter_j = node[jPoint]->GetLimiter();
      //      }
      
      for (iVar = 0; iVar < nVar; iVar++) {
        Project_Grad_i = 0.0; Project_Grad_j = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          Project_Grad_i += Vector_i[iDim]*Gradient_i[iVar][iDim];
          Project_Grad_j += Vector_j[iDim]*Gradient_j[iVar][iDim];
        }
        //        if (limiter) {
        //          Solution_i[iVar] = Turb_i[iVar] + Limiter_i[iVar]*Project_Grad_i;
        //          Solution_j[iVar] = Turb_j[iVar] + Limiter_j[iVar]*Project_Grad_j;
        //        }
        //        else {
        Solution_i[iVar] = Turb_i[iVar] + Project_Grad_i;
        Solution_j[iVar] = Turb_j[iVar] + Project_Grad_j;
        //        }
      }
      
      numerics->SetTurbVar(Solution_i, Solution_j);
      
    }
    
    /*--- Add and subtract residual ---*/
    
    numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
    
    LinSysRes.AddBlock(iPoint, Residual);
    LinSysRes.SubtractBlock(jPoint, Residual);
    
    /*--- Implicit part ---*/
    
    Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
    Jacobian.AddBlock(iPoint, jPoint, Jacobian_j);
    Jacobian.SubtractBlock(jPoint, iPoint, Jacobian_i);
    Jacobian.SubtractBlock(jPoint, jPoint, Jacobian_j);
    
  }
  
}

void CTurbSolver::Viscous_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                   CConfig *config, unsigned short iMesh, unsigned short iRKStep) {
  unsigned long iEdge, iPoint, jPoint;
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    /*--- Points in edge ---*/
    
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    
    /*--- Points coordinates, and normal vector ---*/
    
    numerics->SetCoord(geometry->node[iPoint]->GetCoord(),
                       geometry->node[jPoint]->GetCoord());
    numerics->SetNormal(geometry->edge[iEdge]->GetNormal());
    
    /*--- Conservative variables w/o reconstruction ---*/
    
    numerics->SetPrimitive(solver_container[FLOW_SOL]->node[iPoint]->GetPrimVar(),
                           solver_container[FLOW_SOL]->node[jPoint]->GetPrimVar());
    
    /*--- Turbulent variables w/o reconstruction, and its gradients ---*/
    
    numerics->SetTurbVar(node[iPoint]->GetSolution(), node[jPoint]->GetSolution());
    numerics->SetTurbVarGradient(node[iPoint]->GetGradient(), node[jPoint]->GetGradient());
    
    /*--- Menter's first blending function (only SST)---*/
    if (config->GetKind_Turb_Model() == SST)
      numerics->SetF1blending(node[iPoint]->GetF1blending(),node[jPoint]->GetF1blending());
    
    /*--- Compute residual, and Jacobians ---*/
    
    numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
    
    /*--- Add and subtract residual, and update Jacobians ---*/
    
    LinSysRes.SubtractBlock(iPoint, Residual);
    LinSysRes.AddBlock(jPoint, Residual);
    
    Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
    Jacobian.SubtractBlock(iPoint, jPoint, Jacobian_j);
    Jacobian.AddBlock(jPoint, iPoint, Jacobian_i);
    Jacobian.AddBlock(jPoint, jPoint, Jacobian_j);
    
  }
  
}

void CTurbSolver::BC_Sym_Plane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  /*--- Convective fluxes across symmetry plane are equal to zero. ---*/
}

void CTurbSolver::BC_Euler_Wall(CGeometry *geometry, CSolver **solver_container,
                                CNumerics *numerics, CConfig *config, unsigned short val_marker) {
  /*--- Convective fluxes across euler wall are equal to zero. ---*/
}

void CTurbSolver::ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
  unsigned short iVar;
  unsigned long iPoint, total_index;
  double Delta, Vol, density_old, density;
  
  bool adjoint = config->GetAdjoint();
  
  /*--- Set maximum residual to zero ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }
  
  /*--- Build implicit system ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    
    /*--- Read the volume ---*/
    Vol = geometry->node[iPoint]->GetVolume();
    
    /*--- Modify matrix diagonal to assure diagonal dominance ---*/
    Delta = Vol / (config->GetTurb_CFLRedCoeff()*solver_container[FLOW_SOL]->node[iPoint]->GetDelta_Time());
    Jacobian.AddVal2Diag(iPoint,Delta);
    
    /*--- Right hand side of the system (-Residual) and initial guess (x = 0) ---*/
    for (iVar = 0; iVar < nVar; iVar++) {
      total_index = iPoint*nVar+iVar;
      LinSysRes[total_index] = - LinSysRes[total_index];
      LinSysSol[total_index] = 0.0;
      AddRes_RMS(iVar, LinSysRes[total_index]*LinSysRes[total_index]);
      AddRes_Max(iVar, fabs(LinSysRes[total_index]), geometry->node[iPoint]->GetGlobalIndex());
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
  
  /*--- Solve the linear system (Krylov subspace methods) ---*/
  CMatrixVectorProduct* mat_vec = new CSysMatrixVectorProduct(Jacobian, geometry, config);
  
  CPreconditioner* precond = NULL;
  if (config->GetKind_Linear_Solver_Prec() == JACOBI) {
    Jacobian.BuildJacobiPreconditioner();
    precond = new CJacobiPreconditioner(Jacobian, geometry, config);
  }
  else if (config->GetKind_Linear_Solver_Prec() == LU_SGS) {
    precond = new CLU_SGSPreconditioner(Jacobian, geometry, config);
  }
  else if (config->GetKind_Linear_Solver_Prec() == LINELET) {
    Jacobian.BuildJacobiPreconditioner();
    precond = new CLineletPreconditioner(Jacobian, geometry, config);
  }
  
  CSysSolve system;
  if (config->GetKind_Linear_Solver() == BCGSTAB)
    system.BCGSTAB(LinSysRes, LinSysSol, *mat_vec, *precond, config->GetLinear_Solver_Error(),
                   config->GetLinear_Solver_Iter(), false);
  else if (config->GetKind_Linear_Solver() == FGMRES)
    system.FGMRES(LinSysRes, LinSysSol, *mat_vec, *precond, config->GetLinear_Solver_Error(),
                  config->GetLinear_Solver_Iter(), false);
  
  delete mat_vec;
  delete precond;
  
  /*--- Update solution (system written in terms of increments) ---*/
  
  if (!adjoint) {
    
    /*--- Update and clip trubulent solution ---*/
    
    switch (config->GetKind_Turb_Model()) {
        
      case SA:
        
        for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
          node[iPoint]->AddClippedSolution(0, config->GetLinear_Solver_Relax()*LinSysSol[iPoint],
                                           lowerlimit[0], upperlimit[0]);
        }
        
        break;
        
      case SST:
        
        for (iPoint = 0; iPoint < nPointDomain; iPoint++){
          density_old = solver_container[FLOW_SOL]->node[iPoint]->GetSolution_Old(0);
          density     = solver_container[FLOW_SOL]->node[iPoint]->GetSolution(0);
          
          for (iVar = 0; iVar < nVar; iVar++) {
            node[iPoint]->AddConservativeSolution(iVar, config->GetLinear_Solver_Relax()*LinSysSol[iPoint*nVar+iVar],
                                                  density, density_old, lowerlimit[iVar], upperlimit[iVar]);
          }
          
        }
        
        break;
        
    }
  }
  
  /*--- Compute the root mean square residual ---*/
  
  SetResidual_RMS(geometry, config);
  
}

CTurbSASolver::CTurbSASolver(void) : CTurbSolver() { }

CTurbSASolver::CTurbSASolver(CGeometry *geometry, CConfig *config, unsigned short iMesh) : CTurbSolver() {
  unsigned short iVar, iDim, nLineLets;
  unsigned long iPoint, index;
  double Density_Inf, Viscosity_Inf, Factor_nu_Inf, dull_val;
  
  bool restart = (config->GetRestart() || config->GetRestart_Flow());
  
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
  /*--- Dimension of the problem --> dependent of the turbulent model ---*/
  nVar = 1;
  nPoint = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();
  
  /*--- Define geometry constants in the solver structure ---*/
  nDim = geometry->GetnDim();
  node = new CVariable*[nPoint];
  
  /*--- Single grid simulation ---*/
  if (iMesh == MESH_0) {
    
    /*--- Define some auxiliar vector related with the residual ---*/
    Residual = new double[nVar]; Residual_RMS = new double[nVar];
    Residual_i = new double[nVar]; Residual_j = new double[nVar];
    Residual_Max = new double[nVar]; Point_Max = new unsigned long[nVar];
    
    
    /*--- Define some auxiliar vector related with the solution ---*/
    Solution = new double[nVar];
    Solution_i = new double[nVar]; Solution_j = new double[nVar];
    
    /*--- Define some auxiliar vector related with the geometry ---*/
    Vector_i = new double[nDim]; Vector_j = new double[nDim];
    
    /*--- Define some auxiliar vector related with the flow solution ---*/
    FlowPrimVar_i = new double [nDim+7]; FlowPrimVar_j = new double [nDim+7];
    
    /*--- Jacobians and vector structures for implicit computations ---*/
    Jacobian_i = new double* [nVar];
    Jacobian_j = new double* [nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Jacobian_i[iVar] = new double [nVar];
      Jacobian_j[iVar] = new double [nVar];
    }
    
    /*--- Initialization of the structure of the whole Jacobian ---*/
    cout << "Initialize jacobian structure (SA model)." << endl;
    Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry);
    
    if (config->GetKind_Linear_Solver_Prec() == LINELET) {
      nLineLets = Jacobian.BuildLineletPreconditioner(geometry, config);
      cout << "Compute linelet structure. " << nLineLets << " elements in each line (average)." << endl;
    }
    
    LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
    LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
    
    if (config->GetExtraOutput()) {
      if (nDim == 2){ nOutputVariables = 13; }
      else if (nDim == 3){ nOutputVariables = 19; }
      OutputVariables.Initialize(nPoint, nPointDomain, nOutputVariables, 0.0);
      OutputHeadingNames = new string[nOutputVariables];
    }
    
    /*--- Computation of gradients by least squares ---*/
    if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
      /*--- S matrix := inv(R)*traspose(inv(R)) ---*/
      Smatrix = new double* [nDim];
      for (iDim = 0; iDim < nDim; iDim++)
        Smatrix[iDim] = new double [nDim];
      /*--- c vector := transpose(WA)*(Wb) ---*/
      cvector = new double* [nVar];
      for (iVar = 0; iVar < nVar; iVar++)
        cvector[iVar] = new double [nDim];
    }
    
  }
  
  /* --- Initialize lower and upper limits--- */
  lowerlimit = new double[nVar];
  upperlimit = new double[nVar];
  
  lowerlimit[0] = 1.0e-10;
  upperlimit[0] = 1.0;
  
  /*--- Read farfield conditions from config ---*/
  Density_Inf   = config->GetDensity_FreeStreamND();
  Viscosity_Inf = config->GetViscosity_FreeStreamND();
  
  /*--- Factor_nu_Inf in [3.0, 5.0] ---*/
  Factor_nu_Inf = config->GetNuFactor_FreeStream();
  nu_tilde_Inf  = Factor_nu_Inf*Viscosity_Inf/Density_Inf;
  
  /*--- Eddy viscosity at infinity ---*/
  double Ji, Ji_3, fv1, cv1_3 = 7.1*7.1*7.1;
  double muT_Inf;
  Ji = nu_tilde_Inf/Viscosity_Inf*Density_Inf;
  Ji_3 = Ji*Ji*Ji;
  fv1 = Ji_3/(Ji_3+cv1_3);
  muT_Inf = Density_Inf*fv1*nu_tilde_Inf;
  
  /*--- Restart the solution from file information ---*/
  if (!restart || geometry->GetFinestMGLevel() == false) {
    for (iPoint = 0; iPoint < nPoint; iPoint++)
      node[iPoint] = new CTurbSAVariable(nu_tilde_Inf, muT_Inf, nDim, nVar, config);
  }
  else {
    
    /*--- Restart the solution from file information ---*/
    ifstream restart_file;
    string filename = config->GetSolution_FlowFileName();
    double Density, Laminar_Viscosity, nu, nu_hat, muT, U[5];
    double Temperature, Temperature_Dim, Pressure;
    double Temperature_Ref = config->GetTemperature_Ref();
    double Viscosity_Ref   = config->GetViscosity_Ref();
    double Gas_Constant    = config->GetGas_ConstantND();
    
    /*--- Open the restart file, throw an error if this fails. ---*/
    restart_file.open(filename.data(), ios::in);
    if (restart_file.fail()) {
      cout << "There is no turbulent restart file!!" << endl;
      exit(1);
    }
    
    /*--- In case this is a parallel simulation, we need to perform the
     Global2Local index transformation first. ---*/
    long *Global2Local;
    Global2Local = new long[geometry->GetGlobal_nPointDomain()];
    /*--- First, set all indices to a negative value by default ---*/
    for(iPoint = 0; iPoint < geometry->GetGlobal_nPointDomain(); iPoint++) {
      Global2Local[iPoint] = -1;
    }
    /*--- Now fill array with the transform values only for local points ---*/
    for(iPoint = 0; iPoint < nPointDomain; iPoint++) {
      Global2Local[geometry->node[iPoint]->GetGlobalIndex()] = iPoint;
    }
    
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
      if (iPoint_Local >= 0) {
        
        if (nDim == 2) point_line >> index >> dull_val >> dull_val >> U[0] >> U[1] >> U[2] >> U[3] >> Solution[0];
        if (nDim == 3) point_line >> index >> dull_val >> dull_val >> dull_val >> U[0] >> U[1] >> U[2] >> U[3] >> U[4] >> Solution[0];
        
        Density = U[0];
        if (nDim == 2)
          Pressure = Gamma_Minus_One*(U[3] - (U[1]*U[1] + U[2]*U[2])/(2.0*U[0]));
        else
          Pressure = Gamma_Minus_One*(U[4] - (U[1]*U[1] + U[2]*U[2] + U[3]*U[3])/(2.0*U[0]));
        
        Temperature = Pressure/(Gas_Constant*Density);
        
        /*--- Calculate viscosity from a non-dim. Sutherland's Law ---*/
        Temperature_Dim = Temperature*Temperature_Ref;
        Laminar_Viscosity = 1.853E-5*(pow(Temperature_Dim/300.0,3.0/2.0) * (300.0+110.3)/(Temperature_Dim+110.3));
        Laminar_Viscosity = Laminar_Viscosity/Viscosity_Ref;
        
        nu     = Laminar_Viscosity/Density;
        nu_hat = Solution[0];
        Ji     = nu_hat/nu;
        Ji_3   = Ji*Ji*Ji;
        fv1    = Ji_3/(Ji_3+cv1_3);
        muT    = Density*fv1*nu_hat;
        
        /*--- Instantiate the solution at this node, note that the eddy viscosity should be recomputed ---*/
        node[iPoint_Local] = new CTurbSAVariable(Solution[0], muT, nDim, nVar, config);
      }
      iPoint_Global++;
    }
    
    /*--- Instantiate the variable class with an arbitrary solution
     at any halo/periodic nodes. The initial solution can be arbitrary,
     because a send/recv is performed immediately in the solver. ---*/
    for(iPoint = nPointDomain; iPoint < nPoint; iPoint++) {
      node[iPoint] = new CTurbSAVariable(Solution[0], muT_Inf, nDim, nVar, config);
    }
    
    /*--- Close the restart file ---*/
    restart_file.close();
    
    /*--- Free memory needed for the transformation ---*/
    delete [] Global2Local;
  }
  
}

CTurbSASolver::~CTurbSASolver(void) {
  
}

void CTurbSASolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem) {
  
  unsigned long iPoint;
  
  for (iPoint = 0; iPoint < nPoint; iPoint ++) {
    
    /*--- Initialize the residual vector ---*/
    
    LinSysRes.SetBlock_Zero(iPoint);
    
  }
  
  /*--- Initialize the jacobian matrices ---*/
  
  Jacobian.SetValZero();
  
  /*--- Upwind second order reconstruction ---*/
  
  if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry, config);
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetSolution_Gradient_LS(geometry, config);
  //  if (config->GetKind_SlopeLimit() != NONE) SetSolution_Limiter(geometry, config);
  
}

void CTurbSASolver::Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh) {
  
  double rho, mu, nu, *nu_hat, muT, Ji, Ji_3, fv1;
  double cv1_3 = 7.1*7.1*7.1;
  unsigned long iPoint;
  
  /*--- Compute eddy viscosity ---*/
  
  for (iPoint = 0; iPoint < nPoint; iPoint ++) {
    
    rho = solver_container[FLOW_SOL]->node[iPoint]->GetDensity();
    mu  = solver_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity();
    
    nu  = mu/rho;
    nu_hat = node[iPoint]->GetSolution();
    
    Ji   = nu_hat[0]/nu;
    Ji_3 = Ji*Ji*Ji;
    fv1  = Ji_3/(Ji_3+cv1_3);
    
    muT = rho*fv1*nu_hat[0];
    node[iPoint]->SetmuT(muT);
    
  }
  
}

void CTurbSASolver::Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics,
                                    CConfig *config, unsigned short iMesh) {
  unsigned long iPoint;
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    
    /*--- Conservative variables w/o reconstruction ---*/
    
    numerics->SetPrimitive(solver_container[FLOW_SOL]->node[iPoint]->GetPrimVar(), NULL);
    
    /*--- Gradient of the primitive and conservative variables ---*/
    
    numerics->SetPrimVarGradient(solver_container[FLOW_SOL]->node[iPoint]->GetGradient_Primitive(), NULL);
    
    /*--- Turbulent variables w/o reconstruction, and its gradient ---*/
    
    numerics->SetTurbVar(node[iPoint]->GetSolution(), NULL);
    numerics->SetTurbVarGradient(node[iPoint]->GetGradient(), NULL);
    
    /*--- Set volume ---*/
    
    numerics->SetVolume(geometry->node[iPoint]->GetVolume());
    
    /*--- Set distance to the surface ---*/
    
    numerics->SetDistance(geometry->node[iPoint]->GetWall_Distance(), 0.0);
    
    /*--- Compute the source term ---*/
    
    numerics->ComputeResidual(Residual, Jacobian_i, NULL, config);
    
    /*--- Subtract residual and the jacobian ---*/
    
    LinSysRes.SubtractBlock(iPoint, Residual);
    
    Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
    
  }
  
}

void CTurbSASolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  unsigned long iPoint, iVertex;
  unsigned short iVar;
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Get the velocity vector ---*/
      
      for (iVar = 0; iVar < nVar; iVar++)
        Solution[iVar] = 0.0;
      
      node[iPoint]->SetSolution_Old(Solution);
      LinSysRes.SetBlock_Zero(iPoint);
      
      /*--- includes 1 in the diagonal ---*/
      
      Jacobian.DeleteValsRowi(iPoint);
    }
  }
  
}

void CTurbSASolver::BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                                       unsigned short val_marker) {
  unsigned long iPoint, iVertex;
  unsigned short iVar;
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Get the velocity vector ---*/
      for (iVar = 0; iVar < nVar; iVar++)
        Solution[iVar] = 0.0;
      
      node[iPoint]->SetSolution_Old(Solution);
      LinSysRes.SetBlock_Zero(iPoint);
      
      /*--- Includes 1 in the diagonal ---*/
      
      Jacobian.DeleteValsRowi(iPoint);
    }
  }
  
}

void CTurbSASolver::BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  unsigned long iPoint, iVertex;
  unsigned short iVar, iDim;
  double *Normal, *V_infty, *V_domain;
  
  Normal = new double[nDim];
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Allocate the value at the infinity ---*/
      V_infty = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);
      
      /*--- Retrieve solution at the farfield boundary node ---*/
      V_domain = solver_container[FLOW_SOL]->node[iPoint]->GetPrimVar();
      
      conv_numerics->SetPrimitive(V_domain, V_infty);
      
      /*--- Set turbulent variable at the wall, and at infinity ---*/
      for (iVar = 0; iVar < nVar; iVar++)
        Solution_i[iVar] = node[iPoint]->GetSolution(iVar);
      Solution_j[0] = nu_tilde_Inf;
      conv_numerics->SetTurbVar(Solution_i, Solution_j);
      
      /*--- Set Normal (it is necessary to change the sign) ---*/
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++)
        Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);
      
      /*--- Compute residuals and jacobians ---*/
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      
      /*--- Add residuals and jacobians ---*/
      LinSysRes.AddBlock(iPoint, Residual);
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
    }
  }
  
  delete [] Normal;
  
}

void CTurbSASolver::BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  unsigned short iDim;
  unsigned long iVertex, iPoint, Point_Normal;
  double *V_inlet, *V_domain, *Normal;
  
  Normal = new double[nDim];
  
  string Marker_Tag = config->GetMarker_All_Tag(val_marker);
  
  /*--- Loop over all the vertices on this boundary marker ---*/
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Index of the closest interior node ---*/
      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
      
      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      
      /*--- Allocate the value at the inlet ---*/
      V_inlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);
      
      /*--- Retrieve solution at the farfield boundary node ---*/
      V_domain = solver_container[FLOW_SOL]->node[iPoint]->GetPrimVar();
      
      /*--- Set various quantities in the solver class ---*/
      conv_numerics->SetPrimitive(V_domain, V_inlet);
      
      /*--- Set the turbulent variable states (prescribed for an inflow) ---*/
      Solution_i[0] = node[iPoint]->GetSolution(0);
      Solution_j[0] = nu_tilde_Inf;
      
      conv_numerics->SetTurbVar(Solution_i, Solution_j);
      
      /*--- Set various other quantities in the conv_numerics class ---*/
      conv_numerics->SetNormal(Normal);
      
      /*--- Compute the residual using an upwind scheme ---*/
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Jacobian contribution for implicit integration ---*/
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
      /*--- Viscous contribution ---*/
      visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
      visc_numerics->SetNormal(Normal);
      
      /*--- Conservative variables w/o reconstruction ---*/
      visc_numerics->SetPrimitive(V_domain, V_inlet);
      
      /*--- Turbulent variables w/o reconstruction, and its gradients ---*/
      visc_numerics->SetTurbVar(Solution_i, Solution_j);
      visc_numerics->SetTurbVarGradient(node[iPoint]->GetGradient(), node[iPoint]->GetGradient());
      
      /*--- Compute residual, and Jacobians ---*/
      visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      
      /*--- Subtract residual, and update Jacobians ---*/
      LinSysRes.SubtractBlock(iPoint, Residual);
      Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
      
    }
  }
  
  /*--- Free locally allocated memory ---*/
  delete[] Normal;
  
}

void CTurbSASolver::BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                              CConfig *config, unsigned short val_marker) {
  unsigned long iPoint, iVertex, Point_Normal;
  unsigned short iVar, iDim;
  double *V_outlet, *V_domain, *Normal;
  
  Normal = new double[nDim];
  
  /*--- Loop over all the vertices on this boundary marker ---*/
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Index of the closest interior node ---*/
      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
      
      /*--- Allocate the value at the outlet ---*/
      V_outlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);
      
      /*--- Retrieve solution at the farfield boundary node ---*/
      V_domain = solver_container[FLOW_SOL]->node[iPoint]->GetPrimVar();
      
      /*--- Set various quantities in the solver class ---*/
      conv_numerics->SetPrimitive(V_domain, V_outlet);
      
      /*--- Set the turbulent variables. Here we use a Neumann BC such
       that the turbulent variable is copied from the interior of the
       domain to the outlet before computing the residual.
       Solution_i --> TurbVar_internal,
       Solution_j --> TurbVar_outlet ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        Solution_i[iVar] = node[iPoint]->GetSolution(iVar);
        Solution_j[iVar] = node[iPoint]->GetSolution(iVar);
      }
      conv_numerics->SetTurbVar(Solution_i, Solution_j);
      
      /*--- Set Normal (negate for outward convention) ---*/
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++)
        Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);
      
      /*--- Compute the residual using an upwind scheme ---*/
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Jacobian contribution for implicit integration ---*/
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
      /*--- Viscous contribution ---*/
      visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
      visc_numerics->SetNormal(Normal);
      
      /*--- Conservative variables w/o reconstruction ---*/
      visc_numerics->SetPrimitive(V_domain, V_outlet);
      
      /*--- Turbulent variables w/o reconstruction, and its gradients ---*/
      visc_numerics->SetTurbVar(Solution_i, Solution_j);
      visc_numerics->SetTurbVarGradient(node[iPoint]->GetGradient(), node[iPoint]->GetGradient());
      
      /*--- Compute residual, and Jacobians ---*/
      visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      
      /*--- Subtract residual, and update Jacobians ---*/
      LinSysRes.SubtractBlock(iPoint, Residual);
      Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
      
    }
  }
  
  /*--- Free locally allocated memory ---*/
  delete[] Normal;
  
}

CTurbSSTSolver::CTurbSSTSolver(void) : CTurbSolver() {
  
  /*--- Array initialization ---*/
  constants = NULL;
  
}

CTurbSSTSolver::CTurbSSTSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh) : CTurbSolver() {
  unsigned short iVar, iDim, nLineLets;
  unsigned long iPoint, index;
  double dull_val;
  ifstream restart_file;
  string text_line;
  
  bool restart = (config->GetRestart() || config->GetRestart_Flow());
  
  /*--- Array initialization ---*/
  constants = NULL;
  
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
  /*--- Dimension of the problem --> dependent of the turbulent model ---*/
  nVar = 2;
  nPoint = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();
  
  /*--- Define geometry constants in the solver structure ---*/
  nDim = geometry->GetnDim();
  node = new CVariable*[nPoint];
  
  /*--- Single grid simulation ---*/
  if (iMesh == MESH_0) {
    
    /*--- Define some auxiliary vector related with the residual ---*/
    Residual = new double[nVar]; Residual_RMS = new double[nVar];
    Residual_i = new double[nVar]; Residual_j = new double[nVar];
    Residual_Max = new double[nVar]; Point_Max = new unsigned long[nVar];
    
    
    /*--- Define some auxiliary vector related with the solution ---*/
    Solution = new double[nVar];
    Solution_i = new double[nVar]; Solution_j = new double[nVar];
    
    /*--- Define some auxiliary vector related with the geometry ---*/
    Vector_i = new double[nDim]; Vector_j = new double[nDim];
    
    /*--- Define some auxiliary vector related with the flow solution ---*/
    FlowPrimVar_i = new double [nDim+7]; FlowPrimVar_j = new double [nDim+7];
    
    /*--- Jacobians and vector structures for implicit computations ---*/
    Jacobian_i = new double* [nVar];
    Jacobian_j = new double* [nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Jacobian_i[iVar] = new double [nVar];
      Jacobian_j[iVar] = new double [nVar];
    }
    
    /*--- Initialization of the structure of the whole Jacobian ---*/
    cout << "Initialize jacobian structure (SST model)." << endl;
    Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry);
    
    if (config->GetKind_Linear_Solver_Prec() == LINELET) {
      nLineLets = Jacobian.BuildLineletPreconditioner(geometry, config);
      cout << "Compute linelet structure. " << nLineLets << " elements in each line (average)." << endl;
    }
    
    LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
    LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
  }
  
  /*--- Computation of gradients by least squares ---*/
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
    /*--- S matrix := inv(R)*traspose(inv(R)) ---*/
    Smatrix = new double* [nDim];
    for (iDim = 0; iDim < nDim; iDim++)
      Smatrix[iDim] = new double [nDim];
    /*--- c vector := transpose(WA)*(Wb) ---*/
    cvector = new double* [nVar];
    for (iVar = 0; iVar < nVar; iVar++)
      cvector[iVar] = new double [nDim];
  }
  
  /* --- Initialize value for model constants --- */
  constants = new double[10];
  constants[0] = 0.85;   //sigma_k1
  constants[1] = 1.0;    //sigma_k2
  constants[2] = 0.5;    //sigma_om1
  constants[3] = 0.856;  //sigma_om2
  constants[4] = 0.075;  //beta_1
  constants[5] = 0.0828; //beta_2
  constants[6] = 0.09;   //betaStar
  constants[7] = 0.31;   //a1
  constants[8] = constants[4]/constants[6] - constants[2]*0.41*0.41/sqrt(constants[6]);  //alfa_1
  constants[9] = constants[5]/constants[6] - constants[3]*0.41*0.41/sqrt(constants[6]);  //alfa_2
  
  /* --- Initialize lower and upper limits--- */
  lowerlimit = new double[nVar];
  upperlimit = new double[nVar];
  
  lowerlimit[0] = 1.0e-10;
  upperlimit[0] = 1.0e10;
  
  lowerlimit[1] = 1.0e-4;
  upperlimit[1] = 1.0e15;
  
  /*--- Flow infinity initialization stuff ---*/
  double rhoInf, *VelInf, muLamInf, Intensity, viscRatio, muT_Inf;
  
  rhoInf    = config->GetDensity_FreeStreamND();
  VelInf    = config->GetVelocity_FreeStreamND();
  muLamInf  = config->GetViscosity_FreeStreamND();
  Intensity = config->GetTurbulenceIntensity_FreeStream();
  viscRatio = config->GetTurb2LamViscRatio_FreeStream();
  
  double VelMag = 0;
  for (iDim = 0; iDim < nDim; iDim++)
    VelMag += VelInf[iDim]*VelInf[iDim];
  VelMag = sqrt(VelMag);
  
  kine_Inf  = 3.0/2.0*(VelMag*VelMag*Intensity*Intensity);
  omega_Inf = rhoInf*kine_Inf/(muLamInf*viscRatio);
  
  /*--- Eddy viscosity, initialized without stress limiter at the infinity ---*/
  muT_Inf = rhoInf*kine_Inf/omega_Inf;
  
  /*--- Restart the solution from file information ---*/
  if (!restart || geometry->GetFinestMGLevel() == false) {
    for (iPoint = 0; iPoint < nPoint; iPoint++)
      node[iPoint] = new CTurbSSTVariable(kine_Inf, omega_Inf, muT_Inf, nDim, nVar, constants, config);
  }
  else {
    
    /*--- Restart the solution from file information ---*/
    ifstream restart_file;
    string filename = config->GetSolution_FlowFileName();
    
    /*--- Open the restart file, throw an error if this fails. ---*/
    restart_file.open(filename.data(), ios::in);
    if (restart_file.fail()) {
      cout << "There is no turbulent restart file!!" << endl;
      exit(1);
    }
    
    /*--- In case this is a parallel simulation, we need to perform the
     Global2Local index transformation first. ---*/
    long *Global2Local;
    Global2Local = new long[geometry->GetGlobal_nPointDomain()];
    /*--- First, set all indices to a negative value by default ---*/
    for(iPoint = 0; iPoint < geometry->GetGlobal_nPointDomain(); iPoint++) {
      Global2Local[iPoint] = -1;
    }
    /*--- Now fill array with the transform values only for local points ---*/
    for(iPoint = 0; iPoint < nPointDomain; iPoint++) {
      Global2Local[geometry->node[iPoint]->GetGlobalIndex()] = iPoint;
    }
    
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
      if (iPoint_Local >= 0) {
        
        if (nDim == 2) point_line >> index >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> Solution[0] >> Solution[1];
        if (nDim == 3) point_line >> index >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> Solution[0] >> Solution[1];
        
        /*--- Instantiate the solution at this node, note that the muT_Inf should recomputed ---*/
        node[iPoint_Local] = new CTurbSSTVariable(Solution[0], Solution[1], muT_Inf, nDim, nVar, constants, config);
      }
      iPoint_Global++;
    }
    
    /*--- Instantiate the variable class with an arbitrary solution
     at any halo/periodic nodes. The initial solution can be arbitrary,
     because a send/recv is performed immediately in the solver. ---*/
    for(iPoint = nPointDomain; iPoint < nPoint; iPoint++) {
      node[iPoint] = new CTurbSSTVariable(Solution[0], Solution[1], muT_Inf, nDim, nVar, constants, config);
    }
    
    /*--- Close the restart file ---*/
    restart_file.close();
    
    /*--- Free memory needed for the transformation ---*/
    delete [] Global2Local;
  }
  
}

CTurbSSTSolver::~CTurbSSTSolver(void) {
  
  if (constants != NULL) delete [] constants;
  
}

void CTurbSSTSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem) {
  
  unsigned long iPoint;
  
  for (iPoint = 0; iPoint < nPoint; iPoint ++) {
    
    /*--- Initialize the residual vector ---*/
    
    LinSysRes.SetBlock_Zero(iPoint);
    
  }
  
  /*--- Initialize the jacobian matrices ---*/
  
  Jacobian.SetValZero();
  
  if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry, config);
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetSolution_Gradient_LS(geometry, config);
  //  if (config->GetKind_SlopeLimit() != NONE) SetSolution_Limiter(geometry, config);
  
}

void CTurbSSTSolver::Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh) {
  double rho, mu, dist, omega, kine, vorticity[3], vortMag, strMag, F2, muT, zeta;
  double a1 = constants[7];
  unsigned long iPoint;
  
  /*--- Compute mean flow and turbulence gradients ---*/
  if (config->GetKind_Gradient_Method() == GREEN_GAUSS) {
    solver_container[FLOW_SOL]->SetPrimVar_Gradient_GG(geometry, config);
    SetSolution_Gradient_GG(geometry, config);
  }
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
    solver_container[FLOW_SOL]->SetPrimVar_Gradient_LS(geometry, config);
    SetSolution_Gradient_LS(geometry, config);
  }
  
  for (iPoint = 0; iPoint < nPoint; iPoint ++) {
    /*--- Compute vorticity and rate of strain magnitude ---*/
    solver_container[FLOW_SOL]->node[iPoint]->SetVorticity();
    vorticity[0] = solver_container[FLOW_SOL]->node[iPoint]->GetVorticity(0);
    vorticity[1] = solver_container[FLOW_SOL]->node[iPoint]->GetVorticity(1);
    vorticity[2] = solver_container[FLOW_SOL]->node[iPoint]->GetVorticity(2);
    vortMag = sqrt(vorticity[0]*vorticity[0] + vorticity[1]*vorticity[1] + vorticity[2]*vorticity[2]);
    solver_container[FLOW_SOL]->node[iPoint]->SetStrainMag();
    strMag = solver_container[FLOW_SOL]->node[iPoint]->GetStrainMag();
    
    /*--- Compute blending functions and cross diffusion ---*/
    rho  = solver_container[FLOW_SOL]->node[iPoint]->GetDensity();
    mu   = solver_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity();
    
    dist = geometry->node[iPoint]->GetWall_Distance();
    
    node[iPoint]->SetBlendingFunc(mu, dist, rho);
    F2 = node[iPoint]->GetF2blending();
    
    /*--- Compute the eddy viscosity ---*/
    kine  = node[iPoint]->GetSolution(0);
    omega = node[iPoint]->GetSolution(1);
    zeta = min(1.0/omega,a1/(strMag*F2));
    muT = min(max(rho*kine*zeta,0.0),1.0);
    node[iPoint]->SetmuT(muT);
  }
  
}

void CTurbSSTSolver::Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics,
                                     CConfig *config, unsigned short iMesh) {
  unsigned long iPoint;
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    
    /*--- Conservative variables w/o reconstruction ---*/
    numerics->SetPrimitive(solver_container[FLOW_SOL]->node[iPoint]->GetPrimVar(), NULL);
    
    /*--- Gradient of the primitive and conservative variables ---*/
    numerics->SetPrimVarGradient(solver_container[FLOW_SOL]->node[iPoint]->GetGradient_Primitive(), NULL);
    
    /*--- Turbulent variables w/o reconstruction, and its gradient ---*/
    numerics->SetTurbVar(node[iPoint]->GetSolution(), NULL);
    numerics->SetTurbVarGradient(node[iPoint]->GetGradient(), NULL);
    
    /*--- Set volume ---*/
    numerics->SetVolume(geometry->node[iPoint]->GetVolume());
    
    /*--- Set distance to the surface ---*/
    numerics->SetDistance(geometry->node[iPoint]->GetWall_Distance(), 0.0);
    
    /*--- Menter's first blending function ---*/
    numerics->SetF1blending(node[iPoint]->GetF1blending(),0.0);
    
    /*--- Menter's second blending function ---*/
    numerics->SetF2blending(node[iPoint]->GetF2blending(),0.0);
    
    /*--- Rate of strain magnitude ---*/
    numerics->SetStrainMag(solver_container[FLOW_SOL]->node[iPoint]->GetStrainMag(),0.0);
    
    /*--- Cross diffusion ---*/
    numerics->SetCrossDiff(node[iPoint]->GetCrossDiff(),0.0);
    
    /*--- Compute the source term ---*/
    numerics->ComputeResidual(Residual, Jacobian_i, NULL, config);
    
    /*--- Subtract residual and the jacobian ---*/
    LinSysRes.SubtractBlock(iPoint, Residual);
    Jacobian.SubtractBlock(iPoint,iPoint,Jacobian_i);
    
  }
  
}

void CTurbSSTSolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  unsigned long iPoint, jPoint, iVertex, total_index;
  unsigned short iDim, iVar;
  double distance, density, laminar_viscosity, beta_1;
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- distance to closest neighbor ---*/
      jPoint = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
      distance = 0.0;
      for(iDim = 0; iDim < nDim; iDim++){
        distance += (geometry->node[iPoint]->GetCoord(iDim) - geometry->node[jPoint]->GetCoord(iDim))*
        (geometry->node[iPoint]->GetCoord(iDim) - geometry->node[jPoint]->GetCoord(iDim));
      }
      distance = sqrt(distance);
      
      /*--- Set wall values ---*/
      density = solver_container[FLOW_SOL]->node[jPoint]->GetDensity();
      laminar_viscosity = solver_container[FLOW_SOL]->node[jPoint]->GetLaminarViscosity();
      
      beta_1 = constants[4];
      
      Solution[0] = 0.0;
      Solution[1] = 60.0*laminar_viscosity/(density*beta_1*distance*distance);
      
      /*--- Set the solution values and zero the residual ---*/
      node[iPoint]->SetSolution_Old(Solution);
      node[iPoint]->SetSolution(Solution);
      LinSysRes.SetBlock_Zero(iPoint);
      
      /*--- Change rows of the Jacobian (includes 1 in the diagonal) ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        total_index = iPoint*nVar+iVar;
        Jacobian.DeleteValsRowi(total_index);
      }
      
    }
  }
}

void CTurbSSTSolver::BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                                        unsigned short val_marker) {
  unsigned long iPoint, jPoint, iVertex, total_index;
  unsigned short iDim, iVar;
  double distance, density, laminar_viscosity, beta_1;
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- distance to closest neighbor ---*/
      jPoint = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
      distance = 0.0;
      for(iDim = 0; iDim < nDim; iDim++){
        distance += (geometry->node[iPoint]->GetCoord(iDim) - geometry->node[jPoint]->GetCoord(iDim))*
        (geometry->node[iPoint]->GetCoord(iDim) - geometry->node[jPoint]->GetCoord(iDim));
      }
      distance = sqrt(distance);
      
      /*--- Set wall values ---*/
      density = solver_container[FLOW_SOL]->node[jPoint]->GetDensity();
      laminar_viscosity = solver_container[FLOW_SOL]->node[jPoint]->GetLaminarViscosity();
      
      beta_1 = constants[4];
      
      Solution[0] = 0.0;
      Solution[1] = 60.0*laminar_viscosity/(density*beta_1*distance*distance);
      
      /*--- Set the solution values and zero the residual ---*/
      node[iPoint]->SetSolution_Old(Solution);
      node[iPoint]->SetSolution(Solution);
      LinSysRes.SetBlock_Zero(iPoint);
      
      /*--- Change rows of the Jacobian (includes 1 in the diagonal) ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        total_index = iPoint*nVar+iVar;
        Jacobian.DeleteValsRowi(total_index);
      }
      
    }
  }
}

void CTurbSSTSolver::BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  unsigned long iPoint, iVertex;
  double *Normal, *V_infty, *V_domain;
  unsigned short iVar, iDim;
  
  Normal = new double[nDim];
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Allocate the value at the infinity ---*/
      V_infty = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);
      
      /*--- Retrieve solution at the farfield boundary node ---*/
      V_domain = solver_container[FLOW_SOL]->node[iPoint]->GetPrimVar();
      
      conv_numerics->SetPrimitive(V_domain, V_infty);
      
      /*--- Set turbulent variable at the wall, and at infinity ---*/
      for (iVar = 0; iVar < nVar; iVar++)
        Solution_i[iVar] = node[iPoint]->GetSolution(iVar);
      
      Solution_j[0] = kine_Inf;
      Solution_j[1] = omega_Inf;
      
      conv_numerics->SetTurbVar(Solution_i, Solution_j);
      
      /*--- Set Normal (it is necessary to change the sign) ---*/
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++)
        Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);
      
      /*--- Compute residuals and jacobians ---*/
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      
      /*--- Add residuals and jacobians ---*/
      LinSysRes.AddBlock(iPoint, Residual);
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
    }
  }
  
  delete [] Normal;
}

void CTurbSSTSolver::BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                              unsigned short val_marker) {
  
  unsigned short iVar, iDim;
  unsigned long iVertex, iPoint, Point_Normal;
  double *V_inlet, *V_domain, *Normal;
  
  Normal = new double[nDim];
  
  string Marker_Tag = config->GetMarker_All_Tag(val_marker);
  
  /*--- Loop over all the vertices on this boundary marker ---*/
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Index of the closest interior node ---*/
      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
      
      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      
      /*--- Allocate the value at the inlet ---*/
      V_inlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);
      
      /*--- Retrieve solution at the farfield boundary node ---*/
      V_domain = solver_container[FLOW_SOL]->node[iPoint]->GetPrimVar();
      
      /*--- Set various quantities in the solver class ---*/
      conv_numerics->SetPrimitive(V_domain, V_inlet);
      
      /*--- Set the turbulent variable states. Use free-stream SST
       values for the turbulent state at the inflow. ---*/
      for (iVar = 0; iVar < nVar; iVar++)
        Solution_i[iVar] = node[iPoint]->GetSolution(iVar);
      
      Solution_j[0]= kine_Inf;
      Solution_j[1]= omega_Inf;
      
      conv_numerics->SetTurbVar(Solution_i, Solution_j);
      
      /*--- Set various other quantities in the solver class ---*/
      conv_numerics->SetNormal(Normal);
      
      /*--- Compute the residual using an upwind scheme ---*/
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Jacobian contribution for implicit integration ---*/
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
      /*--- Viscous contribution ---*/
      visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
      visc_numerics->SetNormal(Normal);
      
      /*--- Conservative variables w/o reconstruction ---*/
      visc_numerics->SetPrimitive(V_domain, V_inlet);
      
      /*--- Turbulent variables w/o reconstruction, and its gradients ---*/
      visc_numerics->SetTurbVar(Solution_i, Solution_j);
      visc_numerics->SetTurbVarGradient(node[iPoint]->GetGradient(), node[iPoint]->GetGradient());
      
      /*--- Menter's first blending function ---*/
      visc_numerics->SetF1blending(node[iPoint]->GetF1blending(),node[iPoint]->GetF1blending());
      
      /*--- Compute residual, and Jacobians ---*/
      visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      
      /*--- Subtract residual, and update Jacobians ---*/
      LinSysRes.SubtractBlock(iPoint, Residual);
      Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
      
    }
  }
  
  /*--- Free locally allocated memory ---*/
  delete[] Normal;
}

void CTurbSSTSolver::BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                               CConfig *config, unsigned short val_marker) {
  
  /*--- Local variables and initialization. ---*/
  unsigned long iPoint, iVertex, Point_Normal;
  unsigned short iVar, iDim;
  double *V_outlet, *V_domain, *Normal;
  
  Normal = new double[nDim];
  
  /*--- Loop over all the vertices on this boundary marker ---*/
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Index of the closest interior node ---*/
      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
      
      /*--- Allocate the value at the outlet ---*/
      V_outlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);
      
      /*--- Retrieve solution at the farfield boundary node ---*/
      V_domain = solver_container[FLOW_SOL]->node[iPoint]->GetPrimVar();
      
      /*--- Set various quantities in the solver class ---*/
      conv_numerics->SetPrimitive(V_domain, V_outlet);
      
      /*--- Set the turbulent variables. Here we use a Neumann BC such
       that the turbulent variable is copied from the interior of the
       domain to the outlet before computing the residual.
       Solution_i --> TurbVar_internal,
       Solution_j --> TurbVar_outlet ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        Solution_i[iVar] = node[iPoint]->GetSolution(iVar);
        Solution_j[iVar] = node[iPoint]->GetSolution(iVar);
      }
      conv_numerics->SetTurbVar(Solution_i, Solution_j);
      
      /*--- Set Normal (negate for outward convention) ---*/
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++)
        Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);
      
      /*--- Compute the residual using an upwind scheme ---*/
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Jacobian contribution for implicit integration ---*/
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
      /*--- Viscous contribution ---*/
      visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
      visc_numerics->SetNormal(Normal);
      
      /*--- Conservative variables w/o reconstruction ---*/
      visc_numerics->SetPrimitive(V_domain, V_outlet);
      
      /*--- Turbulent variables w/o reconstruction, and its gradients ---*/
      visc_numerics->SetTurbVar(Solution_i, Solution_j);
      visc_numerics->SetTurbVarGradient(node[iPoint]->GetGradient(), node[iPoint]->GetGradient());
      
      /*--- Menter's first blending function ---*/
      visc_numerics->SetF1blending(node[iPoint]->GetF1blending(),node[iPoint]->GetF1blending());
      
      /*--- Compute residual, and Jacobians ---*/
      visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      
      /*--- Subtract residual, and update Jacobians ---*/
      LinSysRes.SubtractBlock(iPoint, Residual);
      Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
      
    }
  }
  
  /*--- Free locally allocated memory ---*/
  delete[] Normal;
}

double* CTurbSSTSolver::GetConstants() {
  return constants;
}
