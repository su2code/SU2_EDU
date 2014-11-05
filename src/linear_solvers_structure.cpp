 /*!
 * \file linear_solvers_structure.cpp
 * \brief Main classes required for solving linear systems of equations
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

#include "../include/linear_solvers_structure.hpp"

void CSysSolve::ApplyGivens(const double & s, const double & c, double & h1, double & h2) {
  
  double temp = c*h1 + s*h2;
  h2 = c*h2 - s*h1;
  h1 = temp;
}

void CSysSolve::GenerateGivens(double & dx, double & dy, double & s, double & c) {
  
  if ( (dx == 0.0) && (dy == 0.0) ) {
    c = 1.0;
    s = 0.0;
  }
  else if ( fabs(dy) > fabs(dx) ) {
    double tmp = dx/dy;
    dx = sqrt(1.0 + tmp*tmp);
    s = Sign(1.0/dx,dy);
    c = tmp*s;
  }
  else if ( fabs(dy) <= fabs(dx) ) {
    double tmp = dy/dx;
    dy = sqrt(1.0 + tmp*tmp);
    c = Sign(1.0/dy,dx);
    s = tmp*c;
  }
  else {
    // dx and/or dy must be invalid
    dx = 0.0;
    dy = 0.0;
    c = 1.0;
    s = 0.0;
  }
  dx = fabs(dx*dy);
  dy = 0.0;
}

void CSysSolve::SolveReduced(const int & n, const vector<vector<double> > & Hsbg,
                             const vector<double> & rhs, vector<double> & x) {
  // initialize...
  for (int i = 0; i < n; i++)
    x[i] = rhs[i];
  // ... and backsolve
  for (int i = n-1; i >= 0; i--) {
    x[i] /= Hsbg[i][i];
    for (int j = i-1; j >= 0; j--) {
      x[j] -= Hsbg[j][i]*x[i];
    }
  }
}

void CSysSolve::ModGramSchmidt(int i, vector<vector<double> > & Hsbg, vector<CSysVector> & w) {
  
  /*--- Parameter for reorthonormalization ---*/
  static const double reorth = 0.98;
  
  /*--- get the norm of the vector being orthogonalized, and find the
   threshold for re-orthogonalization ---*/
  double nrm = dotProd(w[i+1],w[i+1]);
  double thr = nrm*reorth;
  if (nrm <= 0.0) {
    /*--- The norm of w[i+1] < 0.0 ---*/
    cerr << "CSysSolve::modGramSchmidt: dotProd(w[i+1],w[i+1]) < 0.0" << endl;
    throw(-1);
  }
  else if (nrm != nrm) {
    /*--- This is intended to catch if nrm = NaN, but some optimizations
     may mess it up (according to posts on stackoverflow.com) ---*/
    cerr << "CSysSolve::modGramSchmidt: w[i+1] = NaN" << endl;
    throw(-1);
  }
  
  /*--- Begin main Gram-Schmidt loop ---*/
  for (int k = 0; k < i+1; k++) {
    double prod = dotProd(w[i+1],w[k]);
    Hsbg[k][i] = prod;
    w[i+1].Plus_AX(-prod, w[k]);
    
    /*--- Check if reorthogonalization is necessary ---*/
    if (prod*prod > thr) {
      prod = dotProd(w[i+1],w[k]);
      Hsbg[k][i] += prod;
      w[i+1].Plus_AX(-prod, w[k]);
    }
    
    /*--- Update the norm and check its size ---*/
    nrm -= Hsbg[k][i]*Hsbg[k][i];
    if (nrm < 0.0) nrm = 0.0;
    thr = nrm*reorth;
  }
  
  /*--- Test the resulting vector ---*/
  nrm = w[i+1].norm();
  Hsbg[i+1][i] = nrm;
  if (nrm <= 0.0) {
    /*--- w[i+1] is a linear combination of the w[0:i] ---*/
    cerr << "CSysSolve::modGramSchmidt: w[i+1] linearly dependent on w[0:i]" << endl;
    throw(-1);
  }
  
  /*--- Scale the resulting vector ---*/
  w[i+1] /= nrm;
}

void CSysSolve::WriteHeader(const string & solver, const double & restol, const double & resinit) {
  
  cout << "# " << solver << " residual history" << endl;
  cout << "# Residual tolerance target = " << restol << endl;
  cout << "# Initial residual norm     = " << resinit << endl;
  
}

void CSysSolve::WriteHistory(const int & iter, const double & res, const double & resinit) {
  
  cout << "     " << iter << "     " << res/resinit << endl;
  
}

unsigned long CSysSolve::CG_LinSolver(const CSysVector & b, CSysVector & x, CMatrixVectorProduct & mat_vec,
                                      CPreconditioner & precond, double tol, unsigned long m, bool monitoring) {
	
  int rank = MASTER_NODE;
  
  /*--- Check the subspace size ---*/
  if (m < 1) {
    if (rank == 0) cerr << "CSysSolve::ConjugateGradient: illegal value for subspace size, m = " << m << endl;
    exit(1);
  }
  
  CSysVector r(b);
  CSysVector A_p(b);
  
  /*--- Calculate the initial residual, compute norm, and check if system is already solved ---*/
  mat_vec(x,A_p);
  
  r -= A_p; // recall, r holds b initially
  double norm_r = r.norm();
  double norm0 = b.norm();
  if ( (norm_r < tol*norm0) || (norm_r < eps) ) {
    if (rank == 0) cout << "CSysSolve::ConjugateGradient(): system solved by initial guess." << endl;
    return 0;
  }
  
  double alpha, beta, r_dot_z;
  CSysVector z(r);
  precond(r, z);
  CSysVector p(z);
  
  /*--- Set the norm to the initial initial residual value ---*/
  norm0 = norm_r;
  
  /*--- Output header information including initial residual ---*/
  int i = 0;
  if ((monitoring) && (rank == 0))  {
    WriteHeader("CG", tol, norm_r);
    WriteHistory(i, norm_r, norm0);
  }
  
  /*---  Loop over all search directions ---*/
  for (i = 0; i < m; i++) {
    
    /*--- Apply matrix to p to build Krylov subspace ---*/
    mat_vec(p, A_p);
    
    /*--- Calculate step-length alpha ---*/
    r_dot_z = dotProd(r, z);
    alpha = dotProd(A_p, p);
    alpha = r_dot_z / alpha;
    
    /*--- Update solution and residual: ---*/
    x.Plus_AX(alpha, p);
    r.Plus_AX(-alpha, A_p);
    
    /*--- Check if solution has converged, else output the relative residual if necessary ---*/
    norm_r = r.norm();
    if (norm_r < tol*norm0) break;
    if (((monitoring) && (rank == 0)) && ((i+1) % 5 == 0)) WriteHistory(i+1, norm_r, norm0);
    
    precond(r, z);
    
    /*--- Calculate Gram-Schmidt coefficient beta,
		 beta = dotProd(r_{i+1}, z_{i+1}) / dotProd(r_{i}, z_{i}) ---*/
    beta = 1.0 / r_dot_z;
    r_dot_z = dotProd(r, z);
    beta *= r_dot_z;
    
    /*--- Gram-Schmidt orthogonalization; p = beta *p + z ---*/
    p.Equals_AX_Plus_BY(beta, p, 1.0, z);
  }
  
  
  
  if ((monitoring) && (rank == 0))  {
    cout << "# Conjugate Gradient final (true) residual:" << endl;
    cout << "# Iteration = " << i << ": |res|/|res0| = "  << norm_r/norm0 << endl;
  }
  
  //  /*--- Recalculate final residual (this should be optional) ---*/
  //  mat_vec(x, A_p);
  //  r = b;
  //  r -= A_p;
  //  double true_res = r.norm();
  //
  //  if (fabs(true_res - norm_r) > tol*10.0) {
  //    if (rank == 0) {
  //      cout << "# WARNING in CSysSolve::ConjugateGradient(): " << endl;
  //      cout << "# true residual norm and calculated residual norm do not agree." << endl;
  //      cout << "# true_res - calc_res = " << true_res - norm_r << endl;
  //    }
  //  }
	
	return i;
  
}

unsigned long CSysSolve::FGMRES_LinSolver(const CSysVector & b, CSysVector & x, CMatrixVectorProduct & mat_vec,
                                          CPreconditioner & precond, double tol, unsigned long m, bool monitoring) {
	
  int rank = MASTER_NODE;
  
  /*---  Check the subspace size ---*/
  
  if (m < 1) {
    if (rank == 0) cerr << "CSysSolve::FGMRES: illegal value for subspace size, m = " << m << endl;
    exit(1);
  }
  
  /*---  Check the subspace size ---*/
  
  if (m > 1000) {
    if (rank == 0) cerr << "CSysSolve::FGMRES: illegal value for subspace size (too high), m = " << m << endl;
    exit(1);
  }
  
  /*---  Define various arrays
	 Note: elements in w and z are initialized to x to avoid creating
	 a temporary CSysVector object for the copy constructor ---*/
  
  vector<CSysVector> w(m+1, x);
  vector<CSysVector> z(m+1, x);
  vector<double> g(m+1, 0.0);
  vector<double> sn(m+1, 0.0);
  vector<double> cs(m+1, 0.0);
  vector<double> y(m, 0.0);
  vector<vector<double> > H(m+1, vector<double>(m, 0.0));
  
  /*---  Calculate the norm of the rhs vector ---*/
  
  double norm0 = b.norm();
  
  /*---  Calculate the initial residual (actually the negative residual)
	 and compute its norm ---*/
  
  mat_vec(x,w[0]);
  w[0] -= b;
  
  double beta = w[0].norm();
  
  if ( (beta < tol*norm0) || (beta < eps) ) {
    
    /*---  System is already solved ---*/
    
    if (rank == 0) cout << "CSysSolve::FGMRES(): system solved by initial guess." << endl;
    return 0;
  }
  
  /*---  Normalize residual to get w_{0} (the negative sign is because w[0]
	 holds the negative residual, as mentioned above) ---*/
  
  w[0] /= -beta;
  
  /*---  Initialize the RHS of the reduced system ---*/
  
  g[0] = beta;
  
  /*--- Set the norm to the initial residual value ---*/
  
  norm0 = beta;
  
  /*---  Output header information including initial residual ---*/
  
  int i = 0;
  if ((monitoring) && (rank == 0)) {
    WriteHeader("FGMRES", tol, beta);
    WriteHistory(i, beta, norm0);
  }
  
  /*---  Loop over all search directions ---*/
  
  for (i = 0; i < m; i++) {
    
    /*---  Check if solution has converged ---*/
    
    if (beta < tol*norm0) break;
    
    /*---  Precondition the CSysVector w[i] and store result in z[i] ---*/
    
    precond(w[i], z[i]);
    
    /*---  Add to Krylov subspace ---*/
    
    mat_vec(z[i], w[i+1]);
    
    /*---  Modified Gram-Schmidt orthogonalization ---*/
    
    ModGramSchmidt(i, H, w);
    
    /*---  Apply old Givens rotations to new column of the Hessenberg matrix
		 then generate the new Givens rotation matrix and apply it to
		 the last two elements of H[:][i] and g ---*/
    
    for (int k = 0; k < i; k++)
      ApplyGivens(sn[k], cs[k], H[k][i], H[k+1][i]);
    GenerateGivens(H[i][i], H[i+1][i], sn[i], cs[i]);
    ApplyGivens(sn[i], cs[i], g[i], g[i+1]);
    
    /*---  Set L2 norm of residual and check if solution has converged ---*/
    
    beta = fabs(g[i+1]);
    
    /*---  Output the relative residual if necessary ---*/
    
    if ((((monitoring) && (rank == 0)) && ((i+1) % 100 == 0)) && (rank == 0)) WriteHistory(i+1, beta, norm0);
  }
  
  /*---  Solve the least-squares system and update solution ---*/
  
  SolveReduced(i, H, g, y);
  for (int k = 0; k < i; k++) {
    x.Plus_AX(y[k], z[k]);
  }
  
  if ((monitoring) && (rank == 0)) {
    cout << "# FGMRES final (true) residual:" << endl;
    cout << "# Iteration = " << i << ": |res|/|res0| = " << beta/norm0 << endl;
  }
  
  //  /*---  Recalculate final (neg.) residual (this should be optional) ---*/
  //  mat_vec(x, w[0]);
  //  w[0] -= b;
  //  double res = w[0].norm();
  //
  //  if (fabs(res - beta) > tol*10) {
  //    if (rank == 0) {
  //      cout << "# WARNING in CSysSolve::FGMRES(): " << endl;
  //      cout << "# true residual norm and calculated residual norm do not agree." << endl;
  //      cout << "# res - beta = " << res - beta << endl;
  //    }
  //  }
	
	return i;
  
}

unsigned long CSysSolve::BCGSTAB_LinSolver(const CSysVector & b, CSysVector & x, CMatrixVectorProduct & mat_vec,
                                           CPreconditioner & precond, double tol, unsigned long m, bool monitoring) {
	
  int rank = MASTER_NODE;
  
  /*--- Check the subspace size ---*/
  if (m < 1) {
    if (rank == 0) cerr << "CSysSolve::BCGSTAB: illegal value for subspace size, m = " << m << endl;
    exit(1);
  }
	
  CSysVector r(b);
  CSysVector r_0(b);
  CSysVector p(b);
	CSysVector v(b);
  CSysVector s(b);
	CSysVector t(b);
	CSysVector phat(b);
	CSysVector shat(b);
  CSysVector A_x(b);
  
  /*--- Calculate the initial residual, compute norm, and check if system is already solved ---*/
	mat_vec(x,A_x);
  r -= A_x; r_0 = r; // recall, r holds b initially
  double norm_r = r.norm();
  double norm0 = b.norm();
  if ( (norm_r < tol*norm0) || (norm_r < eps) ) {
    if (rank == 0) cout << "CSysSolve::BCGSTAB(): system solved by initial guess." << endl;
    return 0;
  }
	
	/*--- Initialization ---*/
  double alpha = 1.0, beta = 1.0, omega = 1.0, rho = 1.0, rho_prime = 1.0;
	
  /*--- Set the norm to the initial initial residual value ---*/
  norm0 = norm_r;
  
  /*--- Output header information including initial residual ---*/
  int i = 0;
  if ((monitoring) && (rank == 0)) {
    WriteHeader("BCGSTAB", tol, norm_r);
    WriteHistory(i, norm_r, norm0);
  }
	
  /*---  Loop over all search directions ---*/
  for (i = 0; i < m; i++) {
		
		/*--- Compute rho_prime ---*/
		rho_prime = rho;
		
		/*--- Compute rho_i ---*/
		rho = dotProd(r, r_0);
		
		/*--- Compute beta ---*/
		beta = (rho / rho_prime) * (alpha /omega);
		
		/*--- p_{i} = r_{i-1} + beta * p_{i-1} - beta * omega * v_{i-1} ---*/
		double beta_omega = -beta*omega;
		p.Equals_AX_Plus_BY(beta, p, beta_omega, v);
		p.Plus_AX(1.0, r);
		
		/*--- Preconditioning step ---*/
		precond(p, phat);
		mat_vec(phat, v);
    
		/*--- Calculate step-length alpha ---*/
    double r_0_v = dotProd(r_0, v);
    alpha = rho / r_0_v;
    
		/*--- s_{i} = r_{i-1} - alpha * v_{i} ---*/
		s.Equals_AX_Plus_BY(1.0, r, -alpha, v);
		
		/*--- Preconditioning step ---*/
		precond(s, shat);
		mat_vec(shat, t);
    
		/*--- Calculate step-length omega ---*/
    omega = dotProd(t, s) / dotProd(t, t);
    
		/*--- Update solution and residual: ---*/
    x.Plus_AX(alpha, phat); x.Plus_AX(omega, shat);
		r.Equals_AX_Plus_BY(1.0, s, -omega, t);
    
    /*--- Check if solution has converged, else output the relative residual if necessary ---*/
    norm_r = r.norm();
    if (norm_r < tol*norm0) break;
    if (((monitoring) && (rank == 0)) && ((i+1) % 5 == 0) && (rank == 0)) WriteHistory(i+1, norm_r, norm0);
    
  }
  
  if ((monitoring) && (rank == 0)) {
    cout << "# BCGSTAB final (true) residual:" << endl;
    cout << "# Iteration = " << i << ": |res|/|res0| = "  << norm_r/norm0 << endl;
  }
	
  //  /*--- Recalculate final residual (this should be optional) ---*/
  //	mat_vec(x, A_x);
  //  r = b; r -= A_x;
  //  double true_res = r.norm();
  //
  //  if ((fabs(true_res - norm_r) > tol*10.0) && (rank == 0)) {
  //    cout << "# WARNING in CSysSolve::BCGSTAB(): " << endl;
  //    cout << "# true residual norm and calculated residual norm do not agree." << endl;
  //    cout << "# true_res - calc_res = " << true_res <<" "<< norm_r << endl;
  //  }
	
	return i;
}

void CSysSolve::MultiGrid_LinSolver(CSysMatrix **Jacobian, CSysVector **LinSysRes, CSysVector **LinSysSol, CMatrixVectorProduct & mat_vec, CGeometry **geometry, CConfig *config, unsigned short iMesh, unsigned short mu, double tol, unsigned long m, bool monitoring) {
  
  CSysVector A_x(*LinSysRes[iMesh]);
  CSysVector ResAux(*LinSysRes[iMesh]);
  CSysVector SolAux(*LinSysSol[iMesh]);
  
  unsigned short Smoother = config->GetKind_Linear_Solver_Prec();
  
  /*--- Smooth the solution in the fine grid Jac_h Sol_h = Res_h, because 
   the implementation (assumes that the initial guess is 0)... it is only possible 
   to perform one smoothing ---*/
  
  switch (Smoother) {
    case LU_SGS:
      Jacobian[iMesh]->ComputeLU_SGSPreconditioner(*LinSysRes[iMesh], *LinSysSol[iMesh], geometry[iMesh], config);
      break;
    case JACOBI:
      Jacobian[iMesh]->BuildJacobiPreconditioner();
      Jacobian[iMesh]->ComputeJacobiPreconditioner(*LinSysRes[iMesh], *LinSysSol[iMesh], geometry[iMesh], config);
      break;
    case ILU:
      Jacobian[iMesh]->BuildILUPreconditioner();
      Jacobian[iMesh]->ComputeILUPreconditioner(*LinSysRes[iMesh], *LinSysSol[iMesh], geometry[iMesh], config);
      break;
    case LINELET:
      Jacobian[iMesh]->BuildJacobiPreconditioner();
      Jacobian[iMesh]->ComputeLineletPreconditioner(*LinSysRes[iMesh], *LinSysSol[iMesh], geometry[iMesh], config);
      break;
  }
  
  if (iMesh < config->GetMGLevels()) {
    
    /*--- Compute the residual in the finesh grid ResAux_h = Jac_h.Sol_h - Res_h ---*/
    
    Jacobian[iMesh]->MatrixVectorProduct(*LinSysSol[iMesh], A_x, geometry[iMesh], config);
    
    ResAux -= A_x;
    
    /*--- Restrict the residual to the coarse grid Res_H = I^H_h.ResAux_h ---*/
    
    SetRestricted_Residual(&ResAux, LinSysRes[iMesh+1], geometry[iMesh+1], config);
    
    /*--- Recursive call to MultiGrid_Cycle ---*/
    
    for (unsigned short imu = 0; imu <= mu; imu++) {
      if (iMesh == config->GetMGLevels()-2)
        MultiGrid_LinSolver(Jacobian, LinSysRes, LinSysSol, mat_vec, geometry, config, iMesh+1, 0, tol, m, monitoring);
      else MultiGrid_LinSolver(Jacobian, LinSysRes, LinSysSol, mat_vec, geometry, config, iMesh+1, mu, tol, m, monitoring);
    }
    
    /*--- Prolongate solution to the fine grid solution SolAux_h = I^h_H.Sol_H ---*/
    
    SetProlongated_Solution(&SolAux, LinSysSol[iMesh+1], geometry[iMesh+1], config);

    /*--- Update the fine grid solution Sol_h += SolAux_h ---*/
    
    *LinSysSol[iMesh] += config->GetDamp_Correc_Prolong()*SolAux;
    
  }
 
}

unsigned long CSysSolve::Solve(CSysMatrix **Jacobian, CSysVector **LinSysRes, CSysVector **LinSysSol, CGeometry **geometry,
                               CConfig *config, unsigned short iMesh) {
  
  double SolverTol = config->GetLinear_Solver_Error();
  unsigned long MaxIter = config->GetLinear_Solver_Iter();
  unsigned long IterLinSol = 0;

  /*--- Solve the linear system using a Krylov subspace method ---*/
  
  if (config->GetKind_Linear_Solver() == BCGSTAB || config->GetKind_Linear_Solver() == FGMRES
      || config->GetKind_Linear_Solver() == RFGMRES) {
    
    CMatrixVectorProduct* mat_vec = new CSysMatrixVectorProduct(*Jacobian[iMesh], geometry[iMesh], config);
    
    CPreconditioner* precond = NULL;
    switch (config->GetKind_Linear_Solver_Prec()) {
      case JACOBI:
        Jacobian[iMesh]->BuildJacobiPreconditioner();
        precond = new CJacobiPreconditioner(*Jacobian[iMesh], geometry[iMesh], config);
        break;
      case ILU:
        Jacobian[iMesh]->BuildILUPreconditioner();
        precond = new CILUPreconditioner(*Jacobian[iMesh], geometry[iMesh], config);
        break;
      case LU_SGS:
        precond = new CLU_SGSPreconditioner(*Jacobian[iMesh], geometry[iMesh], config);
        break;
      case LINELET:
        Jacobian[iMesh]->BuildJacobiPreconditioner();
        precond = new CLineletPreconditioner(*Jacobian[iMesh], geometry[iMesh], config);
        break;
    }
    
    switch (config->GetKind_Linear_Solver()) {
      case BCGSTAB:
        IterLinSol = BCGSTAB_LinSolver(*LinSysRes[iMesh], *LinSysSol[iMesh], *mat_vec, *precond, SolverTol, MaxIter, false);
        break;
      case FGMRES:
        IterLinSol = FGMRES_LinSolver(*LinSysRes[iMesh], *LinSysSol[iMesh], *mat_vec, *precond, SolverTol, MaxIter, false);
        break;
      case RFGMRES:
        IterLinSol = 0;
        while (IterLinSol < config->GetLinear_Solver_Iter()) {
          if (IterLinSol + config->GetLinear_Solver_Restart_Frequency() > config->GetLinear_Solver_Iter())
            MaxIter = config->GetLinear_Solver_Iter() - IterLinSol;
          IterLinSol += FGMRES_LinSolver(*LinSysRes[iMesh], *LinSysSol[iMesh], *mat_vec, *precond, SolverTol, MaxIter, false);
          if (LinSysRes[iMesh]->norm() < SolverTol) break;
          SolverTol = SolverTol*(1.0/LinSysRes[iMesh]->norm());
        }
        break;
    }
    
    /*--- Dealocate memory of the Krylov subspace method ---*/
    
    delete mat_vec;
    delete precond;
    
  }
  
  /*--- Solve the linear system using a linear multigrid ---*/

  else if (config->GetKind_Linear_Solver() == MULTIGRID) {
    
    CMatrixVectorProduct* mat_vec = new CSysMatrixVectorProduct(*Jacobian[iMesh], geometry[iMesh], config);

    MultiGrid_LinSolver(Jacobian, LinSysRes, LinSysSol, *mat_vec, geometry, config, iMesh, config->GetMGCycle(), SolverTol, MaxIter, true);
    
    delete mat_vec;
    
  }

  /*--- Smooth the linear system. ---*/
  
  else if (config->GetKind_Linear_Solver() == SMOOTHER_LUSGS || config->GetKind_Linear_Solver() == SMOOTHER_JACOBI
      || config->GetKind_Linear_Solver() == SMOOTHER_ILU || config->GetKind_Linear_Solver() == SMOOTHER_LINELET) {
    switch (config->GetKind_Linear_Solver()) {
      case SMOOTHER_LUSGS:
        Jacobian[iMesh]->ComputeLU_SGSPreconditioner(*LinSysRes[iMesh], *LinSysSol[iMesh], geometry[iMesh], config);
        break;
      case SMOOTHER_JACOBI:
        Jacobian[iMesh]->BuildJacobiPreconditioner();
        Jacobian[iMesh]->ComputeJacobiPreconditioner(*LinSysRes[iMesh], *LinSysSol[iMesh], geometry[iMesh], config);
        break;
      case SMOOTHER_ILU:
        Jacobian[iMesh]->BuildILUPreconditioner();
        Jacobian[iMesh]->ComputeILUPreconditioner(*LinSysRes[iMesh], *LinSysSol[iMesh], geometry[iMesh], config);
        break;
      case SMOOTHER_LINELET:
        Jacobian[iMesh]->BuildJacobiPreconditioner();
        Jacobian[iMesh]->ComputeLineletPreconditioner(*LinSysRes[iMesh], *LinSysSol[iMesh], geometry[iMesh], config);
        break;
        IterLinSol = 1;
    }
  }
  
  return IterLinSol;
  
}

void CSysSolve::SetRestricted_Residual(CSysVector *res_fine, CSysVector *res_coarse, CGeometry *geo_coarse, CConfig *config) {

  /*--- Restric the residual to coarse levels ---*/
  
  unsigned long iVertex, Point_Fine, Point_Coarse;
  unsigned short iMarker, iVar, iChildren, iDim;
  double *Residual_Fine;
  
  unsigned short nVar = res_fine->GetNVar();
  
  double *Residual = new double[nVar];
  
  
  for (Point_Coarse = 0; Point_Coarse < geo_coarse->GetnPointDomain(); Point_Coarse++) {
    
    res_coarse->SetBlock_Zero(Point_Coarse);
    
    for (iVar = 0; iVar < nVar; iVar++) Residual[iVar] = 0.0;
    
    for (iChildren = 0; iChildren < geo_coarse->node[Point_Coarse]->GetnChildren_CV(); iChildren++) {
      Point_Fine = geo_coarse->node[Point_Coarse]->GetChildren_CV(iChildren);
      Residual_Fine = res_fine->GetBlock(Point_Fine);
      for (iVar = 0; iVar < nVar; iVar++) {
        Residual[iVar] += config->GetDamp_Res_Restric()*Residual_Fine[iVar];
      }
    }
    
    res_coarse->SetBlock(Point_Coarse, Residual);
  }
  
  /*--- Set the dirichlet boundary condition (only Navier-Stokes) ---*/
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if ((config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX              ) ||
        (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL             )) {
      for(iVertex = 0; iVertex<geo_coarse->nVertex[iMarker]; iVertex++) {
        Point_Coarse = geo_coarse->vertex[iMarker][iVertex]->GetNode();
        for (iDim = 0; iDim < geo_coarse->GetnDim(); iDim++)
          res_coarse->SetBlock_Zero(Point_Coarse, iDim+1);
      }
    }
  }
  
  delete [] Residual;
  
}

void CSysSolve::SetProlongated_Solution(CSysVector *sol_fine, CSysVector *sol_coarse, CGeometry *geo_coarse, CConfig *config) {
  
  unsigned long Point_Fine, Point_Coarse, iVertex;
  unsigned short iChildren, iMarker, iDim;
  double *Solution_Coarse;

  for (Point_Coarse = 0; Point_Coarse < geo_coarse->GetnPointDomain(); Point_Coarse++) {
    for (iChildren = 0; iChildren < geo_coarse->node[Point_Coarse]->GetnChildren_CV(); iChildren++) {
      Point_Fine = geo_coarse->node[Point_Coarse]->GetChildren_CV(iChildren);
      Solution_Coarse = sol_coarse->GetBlock(Point_Coarse);
      sol_fine->SetBlock(Point_Fine, Solution_Coarse);
    }
  }
  
  /*--- Set the dirichlet boundary condition (only Navier-Stokes) ---*/
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if ((config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX              ) ||
        (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL             )) {
      for(iVertex = 0; iVertex<geo_coarse->nVertex[iMarker]; iVertex++) {
        Point_Coarse = geo_coarse->vertex[iMarker][iVertex]->GetNode();
        for (iDim = 0; iDim < geo_coarse->GetnDim(); iDim++)
          sol_fine->SetBlock_Zero(Point_Coarse, iDim+1);
      }
    }
  }
  
}

