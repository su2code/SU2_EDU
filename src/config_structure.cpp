/*!
 * \file config_structure.cpp
 * \brief Main file for reading the config file.
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

#include "../include/config_structure.hpp"

CConfig::CConfig(char case_filename[200], unsigned short val_software, unsigned short verb_level) {
  
  /*--- Read the config options  ---*/
  
  SetConfig_Options();
  
  /*--- Parse the config file  ---*/
  
  SetParsing(case_filename);
  
  /*--- Configuration file postprocessing ---*/
  
  SetPostprocessing(val_software);
  
  /*--- Set up the boundaries/markers ---*/
  
  SetMarkers(val_software);
  
  /*--- Configuration file console output ---*/
  
  if (verb_level == VERB_HIGH)
    SetOutput(val_software);
  
}


void CConfig::SetConfig_Options() {
  double default_vec_3d[3];
  double default_vec_6d[6];
  nZone = 1;
  iZone = 0;
  
  /*--- Intialize pointers to NULL. If we don't find these values
   in the config file, they will all be set to zero. ---*/
  Kind_GridMovement = NULL;
  Motion_Origin_X = NULL;     Motion_Origin_Y = NULL;     Motion_Origin_Z = NULL;
  Translation_Rate_X = NULL;  Translation_Rate_Y = NULL;  Translation_Rate_Z = NULL;
  Rotation_Rate_X = NULL;     Rotation_Rate_Y = NULL;     Rotation_Rate_Z = NULL;
  Pitching_Omega_X = NULL;    Pitching_Omega_Y = NULL;    Pitching_Omega_Z = NULL;
  Pitching_Ampl_X = NULL;     Pitching_Ampl_Y = NULL;     Pitching_Ampl_Z = NULL;
  Pitching_Phase_X = NULL;    Pitching_Phase_Y = NULL;    Pitching_Phase_Z = NULL;
  Plunging_Omega_X = NULL;    Plunging_Omega_Y = NULL;    Plunging_Omega_Z = NULL;
  Plunging_Ampl_X = NULL;     Plunging_Ampl_Y = NULL;     Plunging_Ampl_Z = NULL;
  RefOriginMoment_X = NULL;   RefOriginMoment_Y = NULL;   RefOriginMoment_Z = NULL;
  MoveMotion_Origin = NULL;
  
  
  /* BEGIN_CONFIG_OPTIONS */
  
  /*--- Options related to problem definition and partitioning ---*/
  /* CONFIG_CATEGORY: Problem Definition */
  
  /* DESCRIPTION: Adjoint type */
  AddEnumOption("REGIME_TYPE", Kind_Regime, Regime_Map, "COMPRESSIBLE");
  
  /* DESCRIPTION: Write extra output */
  AddSpecialOption("EXTRA_OUTPUT", ExtraOutput, SetBoolOption, false);
  
  /* DESCRIPTION: Physical governing equations */
  AddEnumOption("PHYSICAL_PROBLEM", Kind_Solver, Solver_Map, "NONE");
  /* DESCRIPTION: Mathematical problem */
  AddMathProblem("MATH_PROBLEM" , Adjoint, false , OneShot, false, Linearized, false, Restart_Flow, false);
  /* DESCRIPTION: Specify turbulence model */
  AddEnumOption("KIND_TURB_MODEL", Kind_Turb_Model, Turb_Model_Map, "NONE");
  /* DESCRIPTION: Location of the turb model itself */
  AddScalarOption("ML_TURB_MODEL_FILE", ML_Turb_Model_File, string("model.json"));
  /* DESCRIPTION: Location of the check for the proper loading of the turbulence model */
  AddScalarOption("ML_TURB_MODEL_CHECK_FILE", ML_Turb_Model_Check_File, string("check_model.txt"));
  /* DESCRIPTION:  */
  AddScalarOption("MOTION_FILENAME", Motion_Filename, string("mesh_motion.dat"));
  /* DESCRIPTION: Specify transition model */
  AddEnumOption("KIND_TRANS_MODEL", Kind_Trans_Model, Trans_Model_Map, "NONE");
  
  /* DESCRIPTION: Axisymmetric simulation */
  AddSpecialOption("AXISYMMETRIC", Axisymmetric, SetBoolOption, false);
  /* DESCRIPTION: Add the gravity force */
  AddSpecialOption("GRAVITY_FORCE", GravityForce, SetBoolOption, false);
  
  /* DESCRIPTION: Restart solution from native solution file */
  AddSpecialOption("RESTART_SOL", Restart, SetBoolOption, false);
  /* DESCRIPTION: Write a tecplot file for each partition */
  AddSpecialOption("VISUALIZE_PART", Visualize_Partition, SetBoolOption, false);
  
  /*--- Options related to various boundary markers ---*/
  /* CONFIG_CATEGORY: Boundary Markers */
  
  /* DESCRIPTION: Marker(s) of the surface in the surface flow solution file */
  AddMarkerOption("MARKER_PLOTTING", nMarker_Plotting, Marker_Plotting);
  /* DESCRIPTION: Marker(s) of the surface where evaluate the non-dimensional coefficients */
  AddMarkerOption("MARKER_MONITORING", nMarker_Monitoring, Marker_Monitoring);
  /* DESCRIPTION: Marker(s) of the surface where objective function (design problem) will be evaluated */
  AddMarkerOption("MARKER_DESIGNING", nMarker_Designing, Marker_Designing);
  /* DESCRIPTION: Euler wall boundary marker(s) */
  AddMarkerOption("MARKER_EULER", nMarker_Euler, Marker_Euler);
  /* DESCRIPTION: Far-field boundary marker(s) */
  AddMarkerOption("MARKER_FAR", nMarker_FarField, Marker_FarField);
  /* DESCRIPTION: Symmetry boundary condition */
  AddMarkerOption("MARKER_SYM", nMarker_SymWall, Marker_SymWall);
  /* DESCRIPTION: Symmetry boundary condition */
  AddMarkerOption("MARKER_PRESSURE", nMarker_Pressure, Marker_Pressure);
  /* DESCRIPTION: Near-Field boundary condition */
  AddMarkerOption("MARKER_NEARFIELD", nMarker_NearFieldBound, Marker_NearFieldBound);
  /* DESCRIPTION: Zone interface boundary marker(s) */
  AddMarkerOption("MARKER_INTERFACE", nMarker_InterfaceBound, Marker_InterfaceBound);
  /* DESCRIPTION: Dirichlet boundary marker(s) */
  AddMarkerOption("MARKER_DIRICHLET", nMarker_Dirichlet, Marker_Dirichlet);
  /* DESCRIPTION: Neumann boundary marker(s) */
  AddMarkerOption("MARKER_NEUMANN", nMarker_Neumann, Marker_Neumann);
  /* DESCRIPTION: poisson dirichlet boundary marker(s) */
  AddMarkerDirichlet("ELEC_DIRICHLET", nMarker_Dirichlet_Elec, Marker_Dirichlet_Elec, Dirichlet_Value );
  /* DESCRIPTION: poisson neumann boundary marker(s) */
  AddMarkerOption("ELEC_NEUMANN", nMarker_Neumann_Elec, Marker_Neumann_Elec);
  /* DESCRIPTION: Custom boundary marker(s) */
  AddMarkerOption("MARKER_CUSTOM", nMarker_Custom, Marker_Custom);
  /* DESCRIPTION: Periodic boundary marker(s) for use with SU2_PBC
   Format: ( periodic marker, donor marker, rotation_center_x, rotation_center_y,
   rotation_center_z, rotation_angle_x-axis, rotation_angle_y-axis,
   rotation_angle_z-axis, translation_x, translation_y, translation_z, ... ) */
  AddMarkerPeriodic("MARKER_PERIODIC", nMarker_PerBound, Marker_PerBound, Marker_PerDonor,
                    Periodic_RotCenter, Periodic_RotAngles, Periodic_Translation);
  /* DESCRIPTION: Inlet boundary type */
  AddEnumOption("INLET_TYPE", Kind_Inlet, Inlet_Map, "TOTAL_CONDITIONS");
  /* DESCRIPTION: Inlet boundary marker(s) with the following formats,
   Total Conditions: (inlet marker, total temp, total pressure, flow_direction_x,
   flow_direction_y, flow_direction_z, ... ) where flow_direction is
   a unit vector.
   Mass Flow: (inlet marker, density, velocity magnitude, flow_direction_x,
   flow_direction_y, flow_direction_z, ... ) where flow_direction is
   a unit vector. */
  AddMarkerInlet("MARKER_INLET", nMarker_Inlet, Marker_Inlet, Inlet_Ttotal, Inlet_Ptotal, Inlet_FlowDir);
  /* DESCRIPTION: % Supersonic inlet boundary marker(s)
   Format: (inlet marker, temperature, static pressure, velocity_x,
   velocity_y, velocity_z, ... ), i.e. primitive variables specified. */
  AddMarkerInlet("MARKER_SUPERSONIC_INLET", nMarker_Supersonic_Inlet, Marker_Supersonic_Inlet,
                 Inlet_Temperature, Inlet_Pressure, Inlet_Velocity);
  /* DESCRIPTION: Outlet boundary marker(s)
   Format: ( outlet marker, back pressure (static), ... ) */
  AddMarkerOutlet("MARKER_OUTLET", nMarker_Outlet, Marker_Outlet, Outlet_Pressure);
  /* DESCRIPTION: Isothermal wall boundary marker(s)
   Format: ( isothermal marker, wall temperature (static), ... ) */
  AddMarkerOutlet("MARKER_ISOTHERMAL", nMarker_Isothermal, Marker_Isothermal, Isothermal_Temperature);
  /* DESCRIPTION: Specified heat flux wall boundary marker(s)
   Format: ( Heat flux marker, wall heat flux (static), ... ) */
  AddMarkerOutlet("MARKER_HEATFLUX", nMarker_HeatFlux, Marker_HeatFlux, Heat_Flux);
  /* DESCRIPTION: Nacelle inflow boundary marker(s)
   Format: ( nacelle inflow marker, fan face Mach, ... ) */
  AddMarkerOutlet("MARKER_NACELLE_INFLOW", nMarker_NacelleInflow, Marker_NacelleInflow, FanFace_Mach_Target);
  /* DESCRIPTION: Engine subsonic intake region */
  AddSpecialOption("SUBSONIC_NACELLE_INFLOW", Engine_Intake, SetBoolOption, false);
  /* DESCRIPTION: Nacelle exhaust boundary marker(s)
   Format: (nacelle exhaust marker, total nozzle temp, total nozzle pressure, ... )*/
  AddMarkerInlet("MARKER_NACELLE_EXHAUST", nMarker_NacelleExhaust, Marker_NacelleExhaust, Nozzle_Ttotal, Nozzle_Ptotal);
  /* DESCRIPTION: Displacement boundary marker(s) */
  AddMarkerDisplacement("MARKER_NORMAL_DISPL", nMarker_Displacement, Marker_Displacement, Displ_Value);
  /* DESCRIPTION: Load boundary marker(s) */
  AddMarkerLoad("MARKER_NORMAL_LOAD", nMarker_Load, Marker_Load, Load_Value);
  /* DESCRIPTION: Flow load boundary marker(s) */
  AddMarkerFlowLoad("MARKER_FLOWLOAD", nMarker_FlowLoad, Marker_FlowLoad, FlowLoad_Value);
  /* DESCRIPTION: Damping factor for engine inlet condition */
  AddScalarOption("DAMP_NACELLE_INFLOW", Damp_Nacelle_Inflow, 0.1);
  
  /*--- Options related to grid adaptation ---*/
  /* CONFIG_CATEGORY: Grid adaptation */
  
  /* DESCRIPTION: Kind of grid adaptation */
  AddEnumOption("KIND_ADAPT", Kind_Adaptation, Adapt_Map, "NONE");
  /* DESCRIPTION: Percentage of new elements (% of the original number of elements) */
  AddScalarOption("NEW_ELEMS", New_Elem_Adapt, -1.0);
  /* DESCRIPTION: Scale factor for the dual volume */
  AddScalarOption("DUALVOL_POWER", DualVol_Power, 0.5);
  /* DESCRIPTION: Use analytical definition for surfaces */
  AddEnumOption("ANALYTICAL_SURFDEF", Analytical_Surface, Geo_Analytic_Map, "NONE");
  /* DESCRIPTION: Before each computation, implicitly smooth the nodal coordinates */
  AddSpecialOption("SMOOTH_GEOMETRY", SmoothNumGrid, SetBoolOption, false);
  /* DESCRIPTION: Adapt the boundary elements */
  AddSpecialOption("ADAPT_BOUNDARY", AdaptBoundary, SetBoolOption, true);
  /* DESCRIPTION: Divide rectangles into triangles */
  AddSpecialOption("DIVIDE_ELEMENTS", Divide_Element, SetBoolOption, false);
  
  /*--- Options related to time-marching ---*/
  /* CONFIG_CATEGORY: Time-marching */
  
  /* DESCRIPTION: Unsteady simulation  */
  AddEnumOption("UNSTEADY_SIMULATION", Unsteady_Simulation, Unsteady_Map, "NO");
  /* DESCRIPTION:  Courant-Friedrichs-Lewy condition of the finest grid */
  AddScalarOption("CFL_NUMBER", CFLFineGrid, 1.25);
  /* DESCRIPTION: CFL ramp (factor, number of iterations, CFL limit) */
  default_vec_3d[0] = 1.0; default_vec_3d[1] = 100.0; default_vec_3d[2] = 1.0;
  AddArrayOption("CFL_RAMP", 3, CFLRamp, default_vec_3d);
  /* DESCRIPTION: Reduction factor of the CFL coefficient in the adjoint problem */
  AddScalarOption("ADJ_CFL_REDUCTION", Adj_CFLRedCoeff, 0.8);
  /* DESCRIPTION: Reduction factor of the CFL coefficient in the level set problem */
  AddScalarOption("TURB_CFL_REDUCTION", Turb_CFLRedCoeff, 1.0);
  /* DESCRIPTION: Reduction factor of the CFL coefficient in the turbulent adjoint problem */
  AddScalarOption("ADJTURB_CFL_REDUCTION", AdjTurb_CFLRedCoeff, 1.0);
  /* DESCRIPTION: Number of total iterations */
  AddScalarOption("EXT_ITER", nExtIter, 999999);
  // these options share nRKStep as their size, which is not a good idea in general
  /* DESCRIPTION: Runge-Kutta alpha coefficients */
  AddListOption("RK_ALPHA_COEFF", nRKStep, RK_Alpha_Step);
  /* DESCRIPTION: Time Step for dual time stepping simulations (s) */
  AddScalarOption("UNST_TIMESTEP", Delta_UnstTime, 0.0);
  /* DESCRIPTION: Total Physical Time for dual time stepping simulations (s) */
  AddScalarOption("UNST_TIME", Total_UnstTime, 1.0);
  /* DESCRIPTION: Unsteady Courant-Friedrichs-Lewy number of the finest grid */
  AddScalarOption("UNST_CFL_NUMBER", Unst_CFL, 0.0);
  /* DESCRIPTION: Number of internal iterations (dual time method) */
  AddScalarOption("UNST_INT_ITER", Unst_nIntIter, 100);
  /* DESCRIPTION: Integer number of periodic time instances for Time Spectral */
  AddScalarOption("TIME_INSTANCES", nTimeInstances, 1);
  /* DESCRIPTION: Iteration number to begin unsteady restarts (dual time method) */
  AddScalarOption("UNST_RESTART_ITER", Unst_RestartIter, 0);
  /* DESCRIPTION: Starting direct solver iteration for the unsteady adjoint */
  AddScalarOption("UNST_ADJOINT_ITER", Unst_AdjointIter, 0);
  /* DESCRIPTION: Time discretization */
  AddEnumOption("TIME_DISCRE_FLOW", Kind_TimeIntScheme_Flow, Time_Int_Map, "EULER_IMPLICIT");
  /* DESCRIPTION: Time discretization */
  AddEnumOption("TIME_DISCRE_TNE2", Kind_TimeIntScheme_TNE2, Time_Int_Map, "EULER_IMPLICIT");
  /* DESCRIPTION: Time discretization */
  AddEnumOption("TIME_DISCRE_ADJTNE2", Kind_TimeIntScheme_AdjTNE2, Time_Int_Map, "EULER_IMPLICIT");
  /* DESCRIPTION: Time discretization */
  AddEnumOption("TIME_DISCRE_ADJLEVELSET", Kind_TimeIntScheme_AdjLevelSet, Time_Int_Map, "EULER_IMPLICIT");
  /* DESCRIPTION: Time discretization */
  AddEnumOption("TIME_DISCRE_ADJ", Kind_TimeIntScheme_AdjFlow, Time_Int_Map, "EULER_IMPLICIT");
  /* DESCRIPTION: Time discretization */
  AddEnumOption("TIME_DISCRE_LIN", Kind_TimeIntScheme_LinFlow, Time_Int_Map, "EULER_IMPLICIT");
  /* DESCRIPTION: Time discretization */
  AddEnumOption("TIME_DISCRE_TURB", Kind_TimeIntScheme_Turb, Time_Int_Map, "EULER_IMPLICIT");
  /* DESCRIPTION: Time discretization */
  AddEnumOption("TIME_DISCRE_ADJTURB", Kind_TimeIntScheme_AdjTurb, Time_Int_Map, "EULER_IMPLICIT");
  /* DESCRIPTION: Time discretization */
  AddEnumOption("TIME_DISCRE_WAVE", Kind_TimeIntScheme_Wave, Time_Int_Map, "EULER_IMPLICIT");
  /* DESCRIPTION: Time discretization */
  AddEnumOption("TIME_DISCRE_FEA", Kind_TimeIntScheme_FEA, Time_Int_Map, "EULER_IMPLICIT");
  /* DESCRIPTION: Time discretization */
  AddEnumOption("TIME_DISCRE_HEAT", Kind_TimeIntScheme_Heat, Time_Int_Map, "EULER_IMPLICIT");
  /* DESCRIPTION: Time discretization */
  AddEnumOption("TIME_DISCRE_POISSON", Kind_TimeIntScheme_Poisson, Time_Int_Map, "EULER_IMPLICIT");
  
  /*--- Options related to the linear solvers ---*/
  /* CONFIG_CATEGORY: Linear solver definition */
  
  /* DESCRIPTION: Linear solver for the implicit, mesh deformation, or discrete adjoint systems */
  AddEnumOption("LINEAR_SOLVER", Kind_Linear_Solver, Linear_Solver_Map, "FGMRES");
  /* DESCRIPTION: Preconditioner for the Krylov linear solvers */
  AddEnumOption("LINEAR_SOLVER_PREC", Kind_Linear_Solver_Prec, Linear_Solver_Prec_Map, "LU_SGS");
  /* DESCRIPTION: Minimum error threshold for the linear solver for the implicit formulation */
  AddScalarOption("LINEAR_SOLVER_ERROR", Linear_Solver_Error, 1E-5);
  /* DESCRIPTION: Maximum number of iterations of the linear solver for the implicit formulation */
  AddScalarOption("LINEAR_SOLVER_ITER", Linear_Solver_Iter, 2);
  /* DESCRIPTION: Relaxation of the linear solver for the implicit formulation */
  AddScalarOption("LINEAR_SOLVER_RELAX", Linear_Solver_Relax, 1.0);
  /* DESCRIPTION: Roe-Turkel preconditioning for low Mach number flows */
  AddSpecialOption("ROE_TURKEL_PREC", Low_Mach_Precon, SetBoolOption, false);
  /* DESCRIPTION: Time Step for dual time stepping simulations (s) */
  AddScalarOption("MIN_ROE_TURKEL_PREC", Min_Beta_RoeTurkel, 0.01);
  /* DESCRIPTION: Time Step for dual time stepping simulations (s) */
  AddScalarOption("MAX_ROE_TURKEL_PREC", Max_Beta_RoeTurkel, 0.2);
  /* DESCRIPTION: Linear solver for the turbulent adjoint systems */
  AddEnumOption("ADJTURB_LIN_SOLVER", Kind_AdjTurb_Linear_Solver, Linear_Solver_Map, "BCGSTAB");
  /* DESCRIPTION: Preconditioner for the turbulent adjoint Krylov linear solvers */
  AddEnumOption("ADJTURB_LIN_PREC", Kind_AdjTurb_Linear_Prec, Linear_Solver_Prec_Map, "LU_SGS");
  /* DESCRIPTION: Minimum error threshold for the turbulent adjoint linear solver for the implicit formulation */
  AddScalarOption("ADJTURB_LIN_ERROR", AdjTurb_Linear_Error, 1E-5);
  /* DESCRIPTION: Maximum number of iterations of the turbulent adjoint linear solver for the implicit formulation */
  AddScalarOption("ADJTURB_LIN_ITER", AdjTurb_Linear_Iter, 10);
  
  /*--- Options related to dynamic meshes ---*/
  /* CONFIG_CATEGORY: Dynamic mesh definition */
  
  /* DESCRIPTION: Mesh motion for unsteady simulations */
  AddSpecialOption("GRID_MOVEMENT", Grid_Movement, SetBoolOption, false);
  /* DESCRIPTION: Type of mesh motion */
  AddEnumListOption("GRID_MOVEMENT_KIND", nGridMovement, Kind_GridMovement, GridMovement_Map);
  /* DESCRIPTION: Marker(s) of moving surfaces (MOVING_WALL or DEFORMING grid motion). */
  AddMarkerOption("MARKER_MOVING", nMarker_Moving, Marker_Moving);
  /* DESCRIPTION: Mach number (non-dimensional, based on the mesh velocity and freestream vals.) */
  AddScalarOption("MACH_MOTION", Mach_Motion, 0.0);
  /* DESCRIPTION: Coordinates of the rigid motion origin */
  AddListOption("MOTION_ORIGIN_X", nMotion_Origin_X, Motion_Origin_X);
  /* DESCRIPTION: Coordinates of the rigid motion origin */
  AddListOption("MOTION_ORIGIN_Y", nMotion_Origin_Y, Motion_Origin_Y);
  /* DESCRIPTION: Coordinates of the rigid motion origin */
  AddListOption("MOTION_ORIGIN_Z", nMotion_Origin_Z, Motion_Origin_Z);
  /* DESCRIPTION: Translational velocity vector (m/s) in the x, y, & z directions (RIGID_MOTION only) */
  AddListOption("TRANSLATION_RATE_X", nTranslation_Rate_X, Translation_Rate_X);
  /* DESCRIPTION: Translational velocity vector (m/s) in the x, y, & z directions (RIGID_MOTION only) */
  AddListOption("TRANSLATION_RATE_Y", nTranslation_Rate_Y, Translation_Rate_Y);
  /* DESCRIPTION: Translational velocity vector (m/s) in the x, y, & z directions (RIGID_MOTION only) */
  AddListOption("TRANSLATION_RATE_Z", nTranslation_Rate_Z, Translation_Rate_Z);
  /* DESCRIPTION: Angular velocity vector (rad/s) about x, y, & z axes (RIGID_MOTION only) */
  AddListOption("ROTATION_RATE_X", nRotation_Rate_X, Rotation_Rate_X);
  /* DESCRIPTION: Angular velocity vector (rad/s) about x, y, & z axes (RIGID_MOTION only) */
  AddListOption("ROTATION_RATE_Y", nRotation_Rate_Y, Rotation_Rate_Y);
  /* DESCRIPTION: Angular velocity vector (rad/s) about x, y, & z axes (RIGID_MOTION only) */
  AddListOption("ROTATION_RATE_Z", nRotation_Rate_Z, Rotation_Rate_Z);
  /* DESCRIPTION: Pitching angular freq. (rad/s) about x, y, & z axes (RIGID_MOTION only) */
  AddListOption("PITCHING_OMEGA_X", nPitching_Omega_X, Pitching_Omega_X);
  /* DESCRIPTION: Pitching angular freq. (rad/s) about x, y, & z axes (RIGID_MOTION only) */
  AddListOption("PITCHING_OMEGA_Y", nPitching_Omega_Y, Pitching_Omega_Y);
  /* DESCRIPTION: Pitching angular freq. (rad/s) about x, y, & z axes (RIGID_MOTION only) */
  AddListOption("PITCHING_OMEGA_Z", nPitching_Omega_Z, Pitching_Omega_Z);
  /* DESCRIPTION: Pitching amplitude (degrees) about x, y, & z axes (RIGID_MOTION only) */
  AddListOption("PITCHING_AMPL_X", nPitching_Ampl_X, Pitching_Ampl_X);
  /* DESCRIPTION: Pitching amplitude (degrees) about x, y, & z axes (RIGID_MOTION only) */
  AddListOption("PITCHING_AMPL_Y", nPitching_Ampl_Y, Pitching_Ampl_Y);
  /* DESCRIPTION: Pitching amplitude (degrees) about x, y, & z axes (RIGID_MOTION only) */
  AddListOption("PITCHING_AMPL_Z", nPitching_Ampl_Z, Pitching_Ampl_Z);
  /* DESCRIPTION: Pitching phase offset (degrees) about x, y, & z axes (RIGID_MOTION only) */
  AddListOption("PITCHING_PHASE_X", nPitching_Phase_X, Pitching_Phase_X);
  /* DESCRIPTION: Pitching phase offset (degrees) about x, y, & z axes (RIGID_MOTION only) */
  AddListOption("PITCHING_PHASE_Y", nPitching_Phase_Y, Pitching_Phase_Y);
  /* DESCRIPTION: Pitching phase offset (degrees) about x, y, & z axes (RIGID_MOTION only) */
  AddListOption("PITCHING_PHASE_Z", nPitching_Phase_Z, Pitching_Phase_Z);
  /* DESCRIPTION: Plunging angular freq. (rad/s) in x, y, & z directions (RIGID_MOTION only) */
  AddListOption("PLUNGING_OMEGA_X", nPlunging_Omega_X, Plunging_Omega_X);
  /* DESCRIPTION: Plunging angular freq. (rad/s) in x, y, & z directions (RIGID_MOTION only) */
  AddListOption("PLUNGING_OMEGA_Y", nPlunging_Omega_Y, Plunging_Omega_Y);
  /* DESCRIPTION: Plunging angular freq. (rad/s) in x, y, & z directions (RIGID_MOTION only) */
  AddListOption("PLUNGING_OMEGA_Z", nPlunging_Omega_Z, Plunging_Omega_Z);
  /* DESCRIPTION: Plunging amplitude (m) in x, y, & z directions (RIGID_MOTION only) */
  AddListOption("PLUNGING_AMPL_X", nPlunging_Ampl_X, Plunging_Ampl_X);
  /* DESCRIPTION: Plunging amplitude (m) in x, y, & z directions (RIGID_MOTION only) */
  AddListOption("PLUNGING_AMPL_Y", nPlunging_Ampl_Y, Plunging_Ampl_Y);
  /* DESCRIPTION: Plunging amplitude (m) in x, y, & z directions (RIGID_MOTION only) */
  AddListOption("PLUNGING_AMPL_Z", nPlunging_Ampl_Z, Plunging_Ampl_Z);
  /* DESCRIPTION: Value to move motion origins (1 or 0) */
  AddListOption("MOVE_MOTION_ORIGIN", nMoveMotion_Origin, MoveMotion_Origin);
  /* DESCRIPTION:  */
  AddScalarOption("MOTION_FILENAME", Motion_Filename, string("mesh_motion.dat"));
  
  /*--- Options related to convergence ---*/
  /* CONFIG_CATEGORY: Convergence*/
  
  /* DESCRIPTION: Convergence criteria */
  AddEnumOption("CONV_CRITERIA", ConvCriteria, Converge_Crit_Map, "RESIDUAL");
  /* DESCRIPTION: Residual reduction (order of magnitude with respect to the initial value) */
  AddScalarOption("RESIDUAL_REDUCTION", OrderMagResidual, 3.0);
  /* DESCRIPTION: Min value of the residual (log10 of the residual) */
  AddScalarOption("RESIDUAL_MINVAL", MinLogResidual, -8.0);
  /* DESCRIPTION: Iteration number to begin convergence monitoring */
  AddScalarOption("STARTCONV_ITER", StartConv_Iter, 5);
  /* DESCRIPTION: Number of elements to apply the criteria */
  AddScalarOption("CAUCHY_ELEMS", Cauchy_Elems, 100);
  /* DESCRIPTION: Epsilon to control the series convergence */
  AddScalarOption("CAUCHY_EPS", Cauchy_Eps, 1E-10);
  /* DESCRIPTION: Flow functional for the Cauchy criteria */
  AddEnumOption("CAUCHY_FUNC_FLOW", Cauchy_Func_Flow, Objective_Map, "DRAG");
  /* DESCRIPTION: Adjoint functional for the Cauchy criteria */
  AddEnumOption("CAUCHY_FUNC_ADJ", Cauchy_Func_AdjFlow, Sens_Map, "SENS_GEOMETRY");
  /* DESCRIPTION: Linearized functional for the Cauchy criteria */
  AddEnumOption("CAUCHY_FUNC_LIN", Cauchy_Func_LinFlow, Linear_Obj_Map, "DELTA_DRAG");
  /* DESCRIPTION: Epsilon for a full multigrid method evaluation */
  AddScalarOption("FULLMG_CAUCHY_EPS", Cauchy_Eps_FullMG, 1E-4);
  
  /*--- Options related to Multi-grid ---*/
  /* CONFIG_CATEGORY: Multi-grid */
  
  /* DESCRIPTION: Full multi-grid  */
  AddSpecialOption("FULLMG", FullMG, SetBoolOption, false);
  /* DESCRIPTION: Start up iterations using the fine grid only */
  AddScalarOption("START_UP_ITER", nStartUpIter, 0);
  /* DESCRIPTION: Multi-grid Levels */
  AddScalarOption("MGLEVEL", nMultiLevel, 3);
  /* DESCRIPTION: Multi-grid Cycle (0 = V cycle, 1 = W Cycle) */
  AddScalarOption("MGCYCLE", MGCycle, 0);
  /* DESCRIPTION: Multi-grid pre-smoothing level */
  AddListOption("MG_PRE_SMOOTH", nMG_PreSmooth, MG_PreSmooth);
  /* DESCRIPTION: Multi-grid post-smoothing level */
  AddListOption("MG_POST_SMOOTH", nMG_PostSmooth, MG_PostSmooth);
  /* DESCRIPTION: Jacobi implicit smoothing of the correction */
  AddListOption("MG_CORRECTION_SMOOTH", nMG_CorrecSmooth, MG_CorrecSmooth);
  /* DESCRIPTION: Damping factor for the residual restriction */
  AddScalarOption("MG_DAMP_RESTRICTION", Damp_Res_Restric, 0.9);
  /* DESCRIPTION: Damping factor for the correction prolongation */
  AddScalarOption("MG_DAMP_PROLONGATION", Damp_Correc_Prolong, 0.9);
  /* DESCRIPTION: CFL reduction factor on the coarse levels */
  AddScalarOption("MG_CFL_REDUCTION", MG_CFLRedCoeff, 0.9);
  /* DESCRIPTION: Maximum number of children in the agglomeration stage */
  AddScalarOption("MAX_CHILDREN", MaxChildren, 500);
  /* DESCRIPTION: Maximum length of an agglomerated element (relative to the domain) */
  AddScalarOption("MAX_DIMENSION", MaxDimension, 0.1);
  
  /*--- Options related to the spatial discretization ---*/
  /* CONFIG_CATEGORY: Spatial Discretization */
  
  /* DESCRIPTION: Numerical method for spatial gradients */
  AddEnumOption("NUM_METHOD_GRAD", Kind_Gradient_Method, Gradient_Map, "WEIGHTED_LEAST_SQUARES");
  /* DESCRIPTION: Coefficient for the limiter */
  AddScalarOption("LIMITER_COEFF", LimiterCoeff, 0.5);
  
  /* DESCRIPTION: Convective numerical method */
  Kind_ConvNumScheme_Flow = SPACE_CENTERED; Kind_Centered_Flow = JST; Kind_Upwind_Flow = ROE_2ND;
  AddConvectOption("CONV_NUM_METHOD_FLOW", Kind_ConvNumScheme_Flow, Kind_Centered_Flow, Kind_Upwind_Flow);
  /* DESCRIPTION: Viscous numerical method */
  AddEnumOption("VISC_NUM_METHOD_FLOW", Kind_ViscNumScheme_Flow, Viscous_Map, "AVG_GRAD_CORRECTED");
  /* DESCRIPTION: Source term numerical method */
  AddEnumOption("SOUR_NUM_METHOD_FLOW", Kind_SourNumScheme_Flow, Source_Map, "NONE");
  /* DESCRIPTION: Slope limiter */
  AddEnumOption("SLOPE_LIMITER_FLOW", Kind_SlopeLimit_Flow, Limiter_Map, "VENKATAKRISHNAN");
  default_vec_3d[0] = 0.15; default_vec_3d[1] = 0.5; default_vec_3d[2] = 0.02;
  /* DESCRIPTION: 1st, 2nd and 4th order artificial dissipation coefficients */
  AddArrayOption("AD_COEFF_FLOW", 3, Kappa_Flow, default_vec_3d);
  
  /* DESCRIPTION: Slope limiter */
  AddEnumOption("SLOPE_LIMITER_TURB", Kind_SlopeLimit_Turb, Limiter_Map, "VENKATAKRISHNAN");
  /* DESCRIPTION: Convective numerical method */
  AddConvectOption("CONV_NUM_METHOD_TURB", Kind_ConvNumScheme_Turb, Kind_Centered_Turb, Kind_Upwind_Turb);
  /* DESCRIPTION: Viscous numerical method */
  AddEnumOption("VISC_NUM_METHOD_TURB", Kind_ViscNumScheme_Turb, Viscous_Map, "AVG_GRAD_CORRECTED");
  /* DESCRIPTION: Source term numerical method */
  AddEnumOption("SOUR_NUM_METHOD_TURB", Kind_SourNumScheme_Turb, Source_Map, "PIECEWISE_CONSTANT");
  
  /* DESCRIPTION: Source term numerical method */
  AddEnumOption("SOUR_NUM_METHOD_TEMPLATE", Kind_SourNumScheme_Template, Source_Map, "NONE");
  
  
  /*--- Options related to the adjoint and gradient ---*/
  /* CONFIG_CATEGORY: Adjoint and Gradient */
  
  /* DESCRIPTION: Limit value for the adjoint variable */
  AddScalarOption("ADJ_LIMIT", AdjointLimit, 1E6);
  /* DESCRIPTION: Adjoint problem boundary condition */
  AddEnumOption("ADJ_OBJFUNC", Kind_ObjFunc, Objective_Map, "DRAG");
  /* DESCRIPTION: Definition of the airfoil section */
  default_vec_3d[0] = 1E-6; default_vec_3d[1] = 1;
  AddArrayOption("GEO_SECTION_LIMIT", 2, Section_Limit, default_vec_3d);
  /* DESCRIPTION: Mode of the GDC code (analysis, or gradient) */
  AddEnumOption("GEO_MODE", GeometryMode, GeometryMode_Map, "FUNCTION");
  /* DESCRIPTION: Drag weight in sonic boom Objective Function (from 0.0 to 1.0) */
  AddScalarOption("DRAG_IN_SONICBOOM", WeightCd, 0.0);
  /* DESCRIPTION: Sensitivity smoothing  */
  AddEnumOption("SENS_SMOOTHING", Kind_SensSmooth, Sens_Smoothing_Map, "NONE");
  /* DESCRIPTION: Continuous governing equation set  */
  AddEnumOption("CONTINUOUS_EQNS", Continuous_Eqns, ContinuousEqns_Map, "EULER");
  /* DESCRIPTION: Discrete governing equation set */
  AddEnumOption("DISCRETE_EQNS", Discrete_Eqns, DiscreteEqns_Map, "NONE");
  /* DESCRIPTION: Adjoint frozen viscosity */
  AddSpecialOption("FROZEN_VISC", Frozen_Visc, SetBoolOption, true);
  /* DESCRIPTION:  */
  AddScalarOption("CTE_VISCOUS_DRAG", CteViscDrag, 0.0);
  /* DESCRIPTION: Remove sharp edges from the sensitivity evaluation */
  AddSpecialOption("SENS_REMOVE_SHARP", Sens_Remove_Sharp, SetBoolOption, false);
  
  /*--- Options related to input/output files and formats ---*/
  /* CONFIG_CATEGORY: Input/output files and formats */
  
  /* DESCRIPTION: I/O */
  AddEnumOption("OUTPUT_FORMAT", Output_FileFormat, Output_Map, "TECPLOT");
  /* DESCRIPTION: Mesh input file format */
  AddEnumOption("MESH_FORMAT", Mesh_FileFormat, Input_Map, "SU2");
  /* DESCRIPTION: Convert a CGNS mesh to SU2 format */
  AddSpecialOption("CGNS_TO_SU2", CGNS_To_SU2, SetBoolOption, false);
  /* DESCRIPTION:  Mesh input file */
  AddScalarOption("MESH_FILENAME", Mesh_FileName, string("mesh.su2"));
  /* DESCRIPTION: Mesh output file */
  AddScalarOption("MESH_OUT_FILENAME", Mesh_Out_FileName, string("mesh_out.su2"));
  /* DESCRIPTION: Output file convergence history (w/o extension) */
  AddScalarOption("CONV_FILENAME", Conv_FileName, string("history"));
  /* DESCRIPTION: Restart flow input file */
  AddScalarOption("SOLUTION_FLOW_FILENAME", Solution_FlowFileName, string("solution_flow.dat"));
  /* DESCRIPTION: Restart linear flow input file */
  AddScalarOption("SOLUTION_LIN_FILENAME", Solution_LinFileName, string("solution_lin.dat"));
  /* DESCRIPTION: Restart adjoint input file */
  AddScalarOption("SOLUTION_ADJ_FILENAME", Solution_AdjFileName, string("solution_adj.dat"));
  /* DESCRIPTION: Output file restart flow */
  AddScalarOption("RESTART_FLOW_FILENAME", Restart_FlowFileName, string("restart_flow.dat"));
  /* DESCRIPTION: Output file linear flow */
  AddScalarOption("RESTART_LIN_FILENAME",Restart_LinFileName, string("restart_lin.dat"));
  /* DESCRIPTION: Output file restart adjoint */
  AddScalarOption("RESTART_ADJ_FILENAME", Restart_AdjFileName, string("restart_adj.dat"));
  /* DESCRIPTION: Output file restart wave */
  AddScalarOption("RESTART_WAVE_FILENAME", Restart_WaveFileName, string("restart_wave.dat"));
  /* DESCRIPTION: Output file flow (w/o extension) variables */
  AddScalarOption("VOLUME_FLOW_FILENAME", Flow_FileName, string("flow"));
  /* DESCRIPTION: Output file structure (w/o extension) variables */
  AddScalarOption("VOLUME_STRUCTURE_FILENAME", Structure_FileName, string("structure"));
  /* DESCRIPTION: Output file structure (w/o extension) variables */
  AddScalarOption("SURFACE_STRUCTURE_FILENAME", SurfStructure_FileName, string("surface_structure"));
  /* DESCRIPTION: Output file structure (w/o extension) variables */
  AddScalarOption("SURFACE_WAVE_FILENAME", SurfWave_FileName, string("surface_wave"));
  /* DESCRIPTION: Output file structure (w/o extension) variables */
  AddScalarOption("SURFACE_HEAT_FILENAME", SurfHeat_FileName, string("surface_heat"));
  /* DESCRIPTION: Output file wave (w/o extension) variables */
  AddScalarOption("VOLUME_WAVE_FILENAME", Wave_FileName, string("wave"));
  /* DESCRIPTION: Output file wave (w/o extension) variables */
  AddScalarOption("VOLUME_HEAT_FILENAME", Heat_FileName, string("heat"));
  /* DESCRIPTION: Output file adj. wave (w/o extension) variables */
  AddScalarOption("VOLUME_ADJWAVE_FILENAME", AdjWave_FileName, string("adjoint_wave"));
  /* DESCRIPTION: Output file adjoint (w/o extension) variables */
  AddScalarOption("VOLUME_ADJ_FILENAME", Adj_FileName, string("adjoint"));
  /* DESCRIPTION: Output file linear (w/o extension) variables */
  AddScalarOption("VOLUME_LIN_FILENAME", Lin_FileName, string("linearized"));
  /* DESCRIPTION: Output objective function gradient */
  AddScalarOption("GRAD_OBJFUNC_FILENAME", ObjFunc_Grad_FileName, string("of_grad.dat"));
  /* DESCRIPTION: Output objective function */
  AddScalarOption("VALUE_OBJFUNC_FILENAME", ObjFunc_Value_FileName, string("of_func.dat"));
  /* DESCRIPTION: Output file surface flow coefficient (w/o extension) */
  AddScalarOption("SURFACE_FLOW_FILENAME", SurfFlowCoeff_FileName, string("surface_flow"));
  /* DESCRIPTION: Output file surface adjoint coefficient (w/o extension) */
  AddScalarOption("SURFACE_ADJ_FILENAME", SurfAdjCoeff_FileName, string("surface_adjoint"));
  /* DESCRIPTION: Output file surface linear coefficient (w/o extension) */
  AddScalarOption("SURFACE_LIN_FILENAME", SurfLinCoeff_FileName, string("surface_linear"));
  /* DESCRIPTION: Writing solution file frequency */
  AddScalarOption("WRT_SOL_FREQ", Wrt_Sol_Freq, 1000);
  /* DESCRIPTION: Writing solution file frequency */
  AddScalarOption("WRT_SOL_FREQ_DUALTIME", Wrt_Sol_Freq_DualTime, 1);
  /* DESCRIPTION: Writing convergence history frequency */
  AddScalarOption("WRT_CON_FREQ",  Wrt_Con_Freq, 1);
  /* DESCRIPTION: Writing convergence history frequency for the dual time */
  AddScalarOption("WRT_CON_FREQ_DUALTIME",  Wrt_Con_Freq_DualTime, 10);
  /* DESCRIPTION: Write a volume solution file */
  AddSpecialOption("WRT_VOL_SOL", Wrt_Vol_Sol, SetBoolOption, true);
  /* DESCRIPTION: Write a surface solution file */
  AddSpecialOption("WRT_SRF_SOL", Wrt_Srf_Sol, SetBoolOption, true);
  /* DESCRIPTION: Write a surface CSV solution file */
  AddSpecialOption("WRT_CSV_SOL", Wrt_Csv_Sol, SetBoolOption, true);
  /* DESCRIPTION: Write a restart solution file */
  AddSpecialOption("WRT_RESTART", Wrt_Restart, SetBoolOption, true);
  /* DESCRIPTION: Output residual info to solution/restart file */
  AddSpecialOption("WRT_RESIDUALS", Wrt_Residuals, SetBoolOption, false);
  /* DESCRIPTION: Output the rind layers in the solution files */
  AddSpecialOption("WRT_HALO", Wrt_Halo, SetBoolOption, false);
  /* DESCRIPTION: Output sectional forces for specified markers. */
  AddSpecialOption("WRT_SECTIONAL_FORCES", Wrt_Sectional_Forces, SetBoolOption, false);
  
  /*--- Options related to the equivalent area ---*/
  /* CONFIG_CATEGORY: Equivalent Area */
  
  /* DESCRIPTION: Evaluate equivalent area on the Near-Field  */
  AddSpecialOption("EQUIV_AREA", EquivArea, SetBoolOption, false);
  default_vec_3d[0] = 0.0; default_vec_3d[1] = 1.0; default_vec_3d[2] = 1.0;
  /* DESCRIPTION: Integration limits of the equivalent area ( xmin, xmax, Dist_NearField ) */
  AddArrayOption("EA_INT_LIMIT", 3, EA_IntLimit, default_vec_3d);
  
  /*--- Options related to freestream specification ---*/
  /* CONFIG_CATEGORY: Freestream Conditions */
  
  /* DESCRIPTION: Specific gas constant (287.87 J/kg*K (air), only for compressible flows) */
  AddScalarOption("GAS_CONSTANT", Gas_Constant, 287.87);
  /* DESCRIPTION: Ratio of specific heats (1.4 (air), only for compressible flows) */
  AddScalarOption("GAMMA_VALUE", Gamma, 1.4);
  /* DESCRIPTION: Reynolds number (non-dimensional, based on the free-stream values) */
  AddScalarOption("REYNOLDS_NUMBER", Reynolds, 0.0);
  /* DESCRIPTION: Reynolds length (1 m by default) */
  AddScalarOption("REYNOLDS_LENGTH", Length_Reynolds, 1.0);
  /* DESCRIPTION: Laminar Prandtl number (0.72 (air), only for compressible flows) */
  AddScalarOption("PRANDTL_LAM", Prandtl_Lam, 0.72);
  /* DESCRIPTION: Turbulent Prandtl number (0.9 (air), only for compressible flows) */
  AddScalarOption("PRANDTL_TURB", Prandtl_Turb, 0.90);
  /* DESCRIPTION: Value of the Bulk Modulus  */
  AddScalarOption("BULK_MODULUS", Bulk_Modulus, 2.15E9);
  /* DESCRIPTION: Artifical compressibility factor  */
  AddScalarOption("ARTCOMP_FACTOR", ArtComp_Factor, 1.0);
  /* DESCRIPTION:  Mach number (non-dimensional, based on the free-stream values) */
  AddScalarOption("MACH_NUMBER", Mach, 0.0);
  /* DESCRIPTION: Free-stream pressure (101325.0 N/m^2 by default) */
  AddScalarOption("FREESTREAM_PRESSURE", Pressure_FreeStream, 101325.0);
  /* DESCRIPTION: Free-stream density (1.2886 Kg/m^3 (air), 998.2 Kg/m^3 (water)) */
  AddScalarOption("FREESTREAM_DENSITY", Density_FreeStream, -1.0);
  /* DESCRIPTION: Free-stream temperature (273.15 K by default) */
  AddScalarOption("FREESTREAM_TEMPERATURE", Temperature_FreeStream, 273.15);
  /* DESCRIPTION: Free-stream vibrational-electronic temperature (273.15 K by default) */
  AddScalarOption("FREESTREAM_TEMPERATURE_VE", Temperature_ve_FreeStream, 273.15);
  /* DESCRIPTION: Free-stream velocity (m/s) */
  AddArrayOption("FREESTREAM_VELOCITY", 3, Velocity_FreeStream, default_vec_3d);
  /* DESCRIPTION: Free-stream viscosity (1.853E-5 Ns/m^2 (air), 0.798E-3 Ns/m^2 (water)) */
  AddScalarOption("FREESTREAM_VISCOSITY", Viscosity_FreeStream, -1.0);
  /* DESCRIPTION:  */
  AddScalarOption("FREESTREAM_INTERMITTENCY", Intermittency_FreeStream, 1.0);
  /* DESCRIPTION:  */
  AddScalarOption("FREESTREAM_TURBULENCEINTENSITY", TurbulenceIntensity_FreeStream, 0.05);
  /* DESCRIPTION:  */
  AddScalarOption("FREESTREAM_NU_FACTOR", NuFactor_FreeStream, 3.0);
  /* DESCRIPTION:  */
  AddScalarOption("FREESTREAM_TURB2LAMVISCRATIO", Turb2LamViscRatio_FreeStream, 10.0);
  /* DESCRIPTION: Side-slip angle (degrees, only for compressible flows) */
  AddScalarOption("SIDESLIP_ANGLE", AoS, 0.0);
  /* DESCRIPTION: Angle of attack (degrees, only for compressible flows) */
  AddScalarOption("AOA", AoA, 0.0);
  
  /*--- Options related to reference values for nondimensionalization ---*/
  /* CONFIG_CATEGORY: Reference Conditions */
  
  Length_Ref = 1.0; //<---- NOTE: this should be given an option or set as a const
  
  /* DESCRIPTION: X Reference origin for moment computation */
  AddListOption("REF_ORIGIN_MOMENT_X", nRefOriginMoment_X, RefOriginMoment_X);
  /* DESCRIPTION: Y Reference origin for moment computation */
  AddListOption("REF_ORIGIN_MOMENT_Y", nRefOriginMoment_Y, RefOriginMoment_Y);
  /* DESCRIPTION: Z Reference origin for moment computation */
  AddListOption("REF_ORIGIN_MOMENT_Z", nRefOriginMoment_Z, RefOriginMoment_Z);
  /* DESCRIPTION: Reference area for force coefficients (0 implies automatic calculation) */
  AddScalarOption("REF_AREA", RefAreaCoeff, 1.0);
  /* DESCRIPTION: Reference length for pitching, rolling, and yawing non-dimensional moment */
  AddScalarOption("REF_LENGTH_MOMENT", RefLengthMoment, 1.0);
  /* DESCRIPTION: Reference element length for computing the slope limiter epsilon */
  AddScalarOption("REF_ELEM_LENGTH", RefElemLength, 0.1);
  /* DESCRIPTION: Reference coefficient for detecting sharp edges */
  AddScalarOption("REF_SHARP_EDGES", RefSharpEdges, 3.0);
  /* DESCRIPTION: Reference pressure (1.0 N/m^2 by default, only for compressible flows)  */
  AddScalarOption("REF_PRESSURE", Pressure_Ref, 1.0);
  /* DESCRIPTION: Reference temperature (1.0 K by default, only for compressible flows) */
  AddScalarOption("REF_TEMPERATURE", Temperature_Ref, 1.0);
  /* DESCRIPTION: Reference density (1.0 Kg/m^3 by default, only for compressible flows) */
  AddScalarOption("REF_DENSITY", Density_Ref, 1.0);
  /* DESCRIPTION: Reference velocity (incompressible only) */
  AddScalarOption("REF_VELOCITY", Velocity_Ref, -1.0);
  /* DESCRIPTION: Reference viscosity (incompressible only) */
  AddScalarOption("REF_VISCOSITY", Viscosity_Ref, -1.0);
  /* DESCRIPTION: Factor for converting the grid to meters */
  AddScalarOption("CONVERT_TO_METER", Conversion_Factor, 1.0);
  /* DESCRIPTION: Write a new mesh converted to meters */
  AddSpecialOption("WRITE_CONVERTED_MESH", Write_Converted_Mesh, SetBoolOption, false);
  
  /*--- Options related to the reacting gas mixtures ---*/
  /* CONFIG_CATEGORY: Reacting Flow */
  
  /* DESCRIPTION: Specify chemical model for multi-species simulations */
  AddEnumOption("GAS_MODEL", Kind_GasModel, GasModel_Map, "ARGON");
  
  /*--- Options related to the grid deformation ---*/
  // these options share nDV as their size in the option references; not a good idea
  /* CONFIG_CATEGORY: Grid deformation */
  
  /* DESCRIPTION: Kind of deformation */
  AddEnumListOption("DV_KIND", nDV, Design_Variable, Param_Map);
  /* DESCRIPTION: Marker of the surface to which we are going apply the shape deformation */
  AddMarkerOption("DV_MARKER", nMarker_DV, Marker_DV);
  /* DESCRIPTION: New value of the shape deformation */
  AddListOption("DV_VALUE", nDV, DV_Value);
  /* DESCRIPTION: Parameters of the shape deformation
   - HICKS_HENNE ( Lower Surface (0)/Upper Surface (1)/Only one Surface (2), x_Loc )
   - COSINE_BUMP ( Lower Surface (0)/Upper Surface (1)/Only one Surface (2), x_Loc, Thickness )
   - FOURIER ( Lower Surface (0)/Upper Surface (1)/Only one Surface (2), index, cos(0)/sin(1) )
   - NACA_4DIGITS ( 1st digit, 2nd digit, 3rd and 4th digit )
   - PARABOLIC ( Center, Thickness )
   - DISPLACEMENT ( x_Disp, y_Disp, z_Disp )
   - ROTATION ( x_Orig, y_Orig, z_Orig, x_End, y_End, z_End )
   - OBSTACLE ( Center, Bump size )
   - SPHERICAL ( ControlPoint_Index, Theta_Disp, R_Disp )
   - FFD_CONTROL_POINT ( FFDBox ID, i_Ind, j_Ind, k_Ind, x_Disp, y_Disp, z_Disp )
   - FFD_DIHEDRAL_ANGLE ( FFDBox ID, x_Orig, y_Orig, z_Orig, x_End, y_End, z_End )
   - FFD_TWIST_ANGLE ( FFDBox ID, x_Orig, y_Orig, z_Orig, x_End, y_End, z_End )
   - FFD_ROTATION ( FFDBox ID, x_Orig, y_Orig, z_Orig, x_End, y_End, z_End )
   - FFD_CAMBER ( FFDBox ID, i_Ind, j_Ind )
   - FFD_THICKNESS ( FFDBox ID, i_Ind, j_Ind )
   - FFD_VOLUME ( FFDBox ID, i_Ind, j_Ind ) */
  AddDVParamOption("DV_PARAM", nDV, ParamDV, Design_Variable);
  /* DESCRIPTION: Hold the grid fixed in a region */
  AddSpecialOption("HOLD_GRID_FIXED", Hold_GridFixed, SetBoolOption, false);
  default_vec_6d[0] = -1E15; default_vec_6d[1] = -1E15; default_vec_6d[2] = -1E15;
  default_vec_6d[3] =  1E15; default_vec_6d[4] =  1E15; default_vec_6d[5] =  1E15;
  /* DESCRIPTION: Coordinates of the box where the grid will be deformed (Xmin, Ymin, Zmin, Xmax, Ymax, Zmax) */
  AddArrayOption("HOLD_GRID_FIXED_COORD", 6, Hold_GridFixed_Coord, default_vec_6d);
  /* DESCRIPTION: Visualize the deformation */
  AddSpecialOption("VISUALIZE_DEFORMATION", Visualize_Deformation, SetBoolOption, false);
  /* DESCRIPTION: Number of iterations for FEA mesh deformation (surface deformation increments) */
  AddScalarOption("GRID_DEFORM_ITER", GridDef_Iter, 1);
  
  /*--- option related to rotorcraft problems ---*/
  /* CONFIG_CATEGORY: Rotorcraft problem */
  
  AddScalarOption("CYCLIC_PITCH", Cyclic_Pitch, 0.0);
  AddScalarOption("COLLECTIVE_PITCH", Collective_Pitch, 0.0);
  
  
  /*--- Options related to the FEA solver ---*/
  /* CONFIG_CATEGORY: FEA solver */
  
  /* DESCRIPTION: Modulus of elasticity */
  AddScalarOption("ELASTICITY_MODULUS", ElasticyMod, 2E11);
  /* DESCRIPTION: Poisson ratio */
  AddScalarOption("POISSON_RATIO", PoissonRatio, 0.30);
  /* DESCRIPTION: Material density */
  AddScalarOption("MATERIAL_DENSITY", MaterialDensity, 7854);
  
  /*--- options related to the wave solver ---*/
  /* CONFIG_CATEGORY: Wave solver */
  
  /* DESCRIPTION: Constant wave speed */
  AddScalarOption("WAVE_SPEED", Wave_Speed, 331.79);
  
  /*--- options related to the heat solver ---*/
  /* CONFIG_CATEGORY: Heat solver */
  
  /* DESCRIPTION: Thermal diffusivity constant */
  AddScalarOption("THERMAL_DIFFUSIVITY", Thermal_Diffusivity, 1.172E-5);
  
  
  /* END_CONFIG_OPTIONS */
  
}

void CConfig::SetParsing(char case_filename[200]) {
  string text_line, option_name;
  ifstream case_file;
  vector<string> option_value;
    
  /*--- Read the configuration file ---*/
  case_file.open(case_filename, ios::in);
  
  if (case_file.fail()) {
    cout << "There is no configuration file!!" << endl;
    exit(1);
  }
  
  /*--- Parse the configuration file and set the options ---*/
  while (getline (case_file,text_line)) {
    if (TokenizeString(text_line, option_name, option_value)) {
      map<string, CAnyOptionRef*>::iterator it;
      it = param.find(option_name);
      if (it != param.end()) {
        param[option_name]->SetValue(option_value);
      } else {
        if (!GetPython_Option(option_name))
          cout << "WARNING: unrecognized option in the config. file: " << option_name << "." << endl;
      }
    }
  }
  
  case_file.close();
  
}

void CConfig::SetPostprocessing(unsigned short val_software) {
  
  unsigned short iZone;
    
  /*--- Store the SU2 module that we are executing. ---*/
  Kind_SU2 = val_software;
  
  /*--- Only SU2_DDC, and SU2_CFD work with CGNS ---*/
  if ((Kind_SU2 != SU2_DDC) && (Kind_SU2 != SU2_CFD) && (Kind_SU2 != SU2_EDU) && (Kind_SU2 != SU2_SOL)) {
    if (Mesh_FileFormat == CGNS) {
      cout << "This software is not prepared for CGNS, please switch to SU2" << endl;
      exit(1);
    }
  }
  
  /*--- Set default values for the grid based in the Reynolds number for SU2_EDU ---*/
  
    if (Kind_Solver == EULER) Mesh_FileName = "naca0012_inviscid.su2";
    else {
      if (Reynolds < 1E5) Mesh_FileName = "naca0012_re1e5.su2";
      if ((Reynolds >= 1E5) && (Reynolds <= 1E7)) Mesh_FileName = "naca0012_re1e6.su2";
      if (Reynolds > 1E7) Mesh_FileName = "naca0012_re1e7.su2";
    }
  
  /*--- Don't do any deformation if there is no Design variable information ---*/
  if (Design_Variable == NULL) {
    Design_Variable = new unsigned short [1];
    nDV = 1; Design_Variable[0] = NONE;
  }
  
  /*--- If multiple processors the grid should be always in native .su2 format ---*/
  Mesh_FileFormat = SU2;
  
  /*--- Decide whether we should be writing unsteady solution files. ---*/
  Wrt_Unsteady = false;
  
  /*--- Set grid movement kind to NO_MOVEMENT if not specified, which means
   that we also set the Grid_Movement flag to false. We initialize to the
   number of zones here, because we are guaranteed to at least have one. ---*/
  if (Kind_GridMovement == NULL) {
    Kind_GridMovement = new unsigned short[nZone];
    for (unsigned short iZone = 0; iZone < nZone; iZone++ )
      Kind_GridMovement[iZone] = NO_MOVEMENT;
    if (Grid_Movement == true) {
      cout << "GRID_MOVEMENT = YES but no type provided in GRID_MOVEMENT_KIND!!" << endl;
      exit(1);
    }
  }
  
  /*--- If we're solving a purely steady problem with no prescribed grid
   movement (both rotating frame and moving walls can be steady), make sure that
   there is no grid motion ---*/
  if ((Unsteady_Simulation == STEADY) &&
      ((Kind_GridMovement[ZONE_0] != MOVING_WALL) &&
       (Kind_GridMovement[ZONE_0] != ROTATING_FRAME)))
    Grid_Movement = false;
  
  /*--- If it is not specified, set the mesh motion mach number
   equal to the freestream value. ---*/
  if (Grid_Movement && Mach_Motion == 0.0)
    Mach_Motion = Mach;
  
  /*--- Set the boolean flag if we are in a rotating frame (source term). ---*/
  if (Grid_Movement && Kind_GridMovement[ZONE_0] == ROTATING_FRAME)
    Rotating_Frame = true;
  else
    Rotating_Frame = false;
  
  /*--- Check the number of moving markers against the number of grid movement
   types provided (should be equal, except that rigid motion and rotating frame
   do not depend on surface specification). ---*/
  if (Grid_Movement && (Kind_GridMovement[ZONE_0] != RIGID_MOTION) &&
      (Kind_GridMovement[ZONE_0] != ROTATING_FRAME) &&
      (nGridMovement != nMarker_Moving)) {
    cout << "Number of GRID_MOVEMENT_KIND must match number of MARKER_MOVING!!" << endl;
    exit(1);
  }
  
  /*--- Make sure that there aren't more than one rigid motion or
   rotating frame specified in GRID_MOVEMENT_KIND. ---*/
  if (Grid_Movement && (Kind_GridMovement[ZONE_0] == RIGID_MOTION) &&
      (nGridMovement > 1)) {
    cout << "Can not support more than one type of rigid motion in GRID_MOVEMENT_KIND!!" << endl;
    exit(1);
  }
  if (Grid_Movement && (Kind_GridMovement[ZONE_0] == ROTATING_FRAME) &&
      (nGridMovement > 1)) {
    cout << "Can not support more than one rotating frame in GRID_MOVEMENT_KIND!!" << endl;
    exit(1);
  }
  
  /*--- In case the grid movement parameters have not been declared in the
   config file, set them equal to zero for safety. Also check to make sure
   that for each option, a value has been declared for each moving marker. ---*/
  
  unsigned short nMoving;
  if (nGridMovement > nZone) nMoving = nGridMovement;
  else nMoving = nZone;
  
  /*--- Motion Origin: ---*/
  if (Motion_Origin_X == NULL) {
    Motion_Origin_X = new double[nMoving];
    for (iZone = 0; iZone < nMoving; iZone++ )
      Motion_Origin_X[iZone] = 0.0;
  } else {
    if (Grid_Movement && (nMotion_Origin_X != nGridMovement)) {
      cout << "Length of MOTION_ORIGIN_X must match GRID_MOVEMENT_KIND!!" << endl;
      exit(1);
    }
  }
  
  if (Motion_Origin_Y == NULL) {
    Motion_Origin_Y = new double[nMoving];
    for (iZone = 0; iZone < nMoving; iZone++ )
      Motion_Origin_Y[iZone] = 0.0;
  } else {
    if (Grid_Movement && (nMotion_Origin_Y != nGridMovement)) {
      cout << "Length of MOTION_ORIGIN_Y must match GRID_MOVEMENT_KIND!!" << endl;
      exit(1);
    }
  }
  
  if (Motion_Origin_Z == NULL) {
    Motion_Origin_Z = new double[nMoving];
    for (iZone = 0; iZone < nMoving; iZone++ )
      Motion_Origin_Z[iZone] = 0.0;
  } else {
    if (Grid_Movement && (nMotion_Origin_Z != nGridMovement)) {
      cout << "Length of MOTION_ORIGIN_Z must match GRID_MOVEMENT_KIND!!" << endl;
      exit(1);
    }
  }
  
  if (MoveMotion_Origin == NULL) {
    MoveMotion_Origin = new unsigned short[nMoving];
    for (iZone = 0; iZone < nMoving; iZone++ )
      MoveMotion_Origin[iZone] = 0;
  } else {
    if (Grid_Movement && (nMoveMotion_Origin != nGridMovement)) {
      cout << "Length of MOVE_MOTION_ORIGIN must match GRID_MOVEMENT_KIND!!" << endl;
      exit(1);
    }
  }
  
  /*--- Translation: ---*/
  if (Translation_Rate_X == NULL) {
    Translation_Rate_X = new double[nMoving];
    for (iZone = 0; iZone < nMoving; iZone++ )
      Translation_Rate_X[iZone] = 0.0;
  } else {
    if (Grid_Movement && (nTranslation_Rate_X != nGridMovement)) {
      cout << "Length of TRANSLATION_RATE_X must match GRID_MOVEMENT_KIND!!" << endl;
      exit(1);
    }
  }
  
  if (Translation_Rate_Y == NULL) {
    Translation_Rate_Y = new double[nMoving];
    for (iZone = 0; iZone < nMoving; iZone++ )
      Translation_Rate_Y[iZone] = 0.0;
  } else {
    if (Grid_Movement && (nTranslation_Rate_Y != nGridMovement)) {
      cout << "Length of TRANSLATION_RATE_Y must match GRID_MOVEMENT_KIND!!" << endl;
      exit(1);
    }
  }
  
  if (Translation_Rate_Z == NULL) {
    Translation_Rate_Z = new double[nMoving];
    for (iZone = 0; iZone < nMoving; iZone++ )
      Translation_Rate_Z[iZone] = 0.0;
  } else {
    if (Grid_Movement && (nTranslation_Rate_Z != nGridMovement)) {
      cout << "Length of TRANSLATION_RATE_Z must match GRID_MOVEMENT_KIND!!" << endl;
      exit(1);
    }
  }
  
  /*--- Rotation: ---*/
  if (Rotation_Rate_X == NULL) {
    Rotation_Rate_X = new double[nMoving];
    for (iZone = 0; iZone < nMoving; iZone++ )
      Rotation_Rate_X[iZone] = 0.0;
  } else {
    if (Grid_Movement && (nRotation_Rate_X != nGridMovement)) {
      cout << "Length of ROTATION_RATE_X must match GRID_MOVEMENT_KIND!!" << endl;
      exit(1);
    }
  }
  
  if (Rotation_Rate_Y == NULL) {
    Rotation_Rate_Y = new double[nMoving];
    for (iZone = 0; iZone < nMoving; iZone++ )
      Rotation_Rate_Y[iZone] = 0.0;
  } else {
    if (Grid_Movement && (nRotation_Rate_Y != nGridMovement)) {
      cout << "Length of ROTATION_RATE_Y must match GRID_MOVEMENT_KIND!!" << endl;
      exit(1);
    }
  }
  
  if (Rotation_Rate_Z == NULL) {
    Rotation_Rate_Z = new double[nMoving];
    for (iZone = 0; iZone < nMoving; iZone++ )
      Rotation_Rate_Z[iZone] = 0.0;
  } else {
    if (Grid_Movement && (nRotation_Rate_Z != nGridMovement)) {
      cout << "Length of ROTATION_RATE_Z must match GRID_MOVEMENT_KIND!!" << endl;
      exit(1);
    }
  }
  
  /*--- Pitching: ---*/
  if (Pitching_Omega_X == NULL) {
    Pitching_Omega_X = new double[nMoving];
    for (iZone = 0; iZone < nMoving; iZone++ )
      Pitching_Omega_X[iZone] = 0.0;
  } else {
    if (Grid_Movement && (nPitching_Omega_X != nGridMovement)) {
      cout << "Length of PITCHING_OMEGA_X must match GRID_MOVEMENT_KIND!!" << endl;
      exit(1);
    }
  }
  
  if (Pitching_Omega_Y == NULL) {
    Pitching_Omega_Y = new double[nMoving];
    for (iZone = 0; iZone < nMoving; iZone++ )
      Pitching_Omega_Y[iZone] = 0.0;
  } else {
    if (Grid_Movement && (nPitching_Omega_Y != nGridMovement)) {
      cout << "Length of PITCHING_OMEGA_Y must match GRID_MOVEMENT_KIND!!" << endl;
      exit(1);
    }
  }
  
  if (Pitching_Omega_Z == NULL) {
    Pitching_Omega_Z = new double[nMoving];
    for (iZone = 0; iZone < nMoving; iZone++ )
      Pitching_Omega_Z[iZone] = 0.0;
  } else {
    if (Grid_Movement && (nPitching_Omega_Z != nGridMovement)) {
      cout << "Length of PITCHING_OMEGA_Z must match GRID_MOVEMENT_KIND!!" << endl;
      exit(1);
    }
  }
  
  /*--- Pitching Amplitude: ---*/
  if (Pitching_Ampl_X == NULL) {
    Pitching_Ampl_X = new double[nMoving];
    for (iZone = 0; iZone < nMoving; iZone++ )
      Pitching_Ampl_X[iZone] = 0.0;
  } else {
    if (Grid_Movement && (nPitching_Ampl_X != nGridMovement)) {
      cout << "Length of PITCHING_AMPL_X must match GRID_MOVEMENT_KIND!!" << endl;
      exit(1);
    }
  }
  
  if (Pitching_Ampl_Y == NULL) {
    Pitching_Ampl_Y = new double[nMoving];
    for (iZone = 0; iZone < nMoving; iZone++ )
      Pitching_Ampl_Y[iZone] = 0.0;
  } else {
    if (Grid_Movement && (nPitching_Ampl_Y != nGridMovement)) {
      cout << "Length of PITCHING_AMPL_Y must match GRID_MOVEMENT_KIND!!" << endl;
      exit(1);
    }
  }
  
  if (Pitching_Ampl_Z == NULL) {
    Pitching_Ampl_Z = new double[nMoving];
    for (iZone = 0; iZone < nMoving; iZone++ )
      Pitching_Ampl_Z[iZone] = 0.0;
  } else {
    if (Grid_Movement && (nPitching_Ampl_Z != nGridMovement)) {
      cout << "Length of PITCHING_AMPL_Z must match GRID_MOVEMENT_KIND!!" << endl;
      exit(1);
    }
  }
  
  /*--- Pitching Phase: ---*/
  if (Pitching_Phase_X == NULL) {
    Pitching_Phase_X = new double[nMoving];
    for (iZone = 0; iZone < nMoving; iZone++ )
      Pitching_Phase_X[iZone] = 0.0;
  } else {
    if (Grid_Movement && (nPitching_Phase_X != nGridMovement)) {
      cout << "Length of PITCHING_PHASE_X must match GRID_MOVEMENT_KIND!!" << endl;
      exit(1);
    }
  }
  
  if (Pitching_Phase_Y == NULL) {
    Pitching_Phase_Y = new double[nMoving];
    for (iZone = 0; iZone < nMoving; iZone++ )
      Pitching_Phase_Y[iZone] = 0.0;
  } else {
    if (Grid_Movement && (nPitching_Phase_Y != nGridMovement)) {
      cout << "Length of PITCHING_PHASE_Y must match GRID_MOVEMENT_KIND!!" << endl;
      exit(1);
    }
  }
  
  if (Pitching_Phase_Z == NULL) {
    Pitching_Phase_Z = new double[nMoving];
    for (iZone = 0; iZone < nMoving; iZone++ )
      Pitching_Phase_Z[iZone] = 0.0;
  } else {
    if (Grid_Movement && (nPitching_Phase_Z != nGridMovement)) {
      cout << "Length of PITCHING_PHASE_Z must match GRID_MOVEMENT_KIND!!" << endl;
      exit(1);
    }
  }
  
  /*--- Plunging: ---*/
  if (Plunging_Omega_X == NULL) {
    Plunging_Omega_X = new double[nMoving];
    for (iZone = 0; iZone < nMoving; iZone++ )
      Plunging_Omega_X[iZone] = 0.0;
  } else {
    if (Grid_Movement && (nPlunging_Omega_X != nGridMovement)) {
      cout << "Length of PLUNGING_OMEGA_X must match GRID_MOVEMENT_KIND!!" << endl;
      exit(1);
    }
  }
  
  if (Plunging_Omega_Y == NULL) {
    Plunging_Omega_Y = new double[nMoving];
    for (iZone = 0; iZone < nMoving; iZone++ )
      Plunging_Omega_Y[iZone] = 0.0;
  } else {
    if (Grid_Movement && (nPlunging_Omega_Y != nGridMovement)) {
      cout << "Length of PLUNGING_OMEGA_Y must match GRID_MOVEMENT_KIND!!" << endl;
      exit(1);
    }
  }
  
  if (Plunging_Omega_Z == NULL) {
    Plunging_Omega_Z = new double[nMoving];
    for (iZone = 0; iZone < nMoving; iZone++ )
      Plunging_Omega_Z[iZone] = 0.0;
  } else {
    if (Grid_Movement && (nPlunging_Omega_Z != nGridMovement)) {
      cout << "Length of PLUNGING_OMEGA_Z must match GRID_MOVEMENT_KIND!!" << endl;
      exit(1);
    }
  }
  
  /*--- Plunging Amplitude: ---*/
  if (Plunging_Ampl_X == NULL) {
    Plunging_Ampl_X = new double[nMoving];
    for (iZone = 0; iZone < nMoving; iZone++ )
      Plunging_Ampl_X[iZone] = 0.0;
  } else {
    if (Grid_Movement && (nPlunging_Ampl_X != nGridMovement)) {
      cout << "Length of PLUNGING_AMPL_X must match GRID_MOVEMENT_KIND!!" << endl;
      exit(1);
    }
  }
  
  if (Plunging_Ampl_Y == NULL) {
    Plunging_Ampl_Y = new double[nMoving];
    for (iZone = 0; iZone < nMoving; iZone++ )
      Plunging_Ampl_Y[iZone] = 0.0;
  } else {
    if (Grid_Movement && (nPlunging_Ampl_Y != nGridMovement)) {
      cout << "Length of PLUNGING_AMPL_Y must match GRID_MOVEMENT_KIND!!" << endl;
      exit(1);
    }
  }
  
  if (Plunging_Ampl_Z == NULL) {
    Plunging_Ampl_Z = new double[nMoving];
    for (iZone = 0; iZone < nMoving; iZone++ )
      Plunging_Ampl_Z[iZone] = 0.0;
  } else {
    if (Grid_Movement && (nPlunging_Ampl_Z != nGridMovement)) {
      cout << "Length of PLUNGING_AMPL_Z must match GRID_MOVEMENT_KIND!!" << endl;
      exit(1);
    }
  }
  
  
  /*--- Initialize the RefOriginMoment Pointer ---*/
  RefOriginMoment = NULL;
  RefOriginMoment = new double[3];
  RefOriginMoment[0] = 0.0; RefOriginMoment[1] = 0.0; RefOriginMoment[2] = 0.0;
  
  /*--- In case the moment origin coordinates have not been declared in the
   config file, set them equal to zero for safety. Also check to make sure
   that for each marker, a value has been declared for the moment origin.
   Unless only one value was specified, then set this value for all the markers
   being monitored. ---*/
  
  unsigned short iMarker;
  
  
  if ((nRefOriginMoment_X != nRefOriginMoment_Y) || (nRefOriginMoment_X != nRefOriginMoment_Z) ) {
    cout << "ERROR: Length of REF_ORIGIN_MOMENT_X, REF_ORIGIN_MOMENT_Y and REF_ORIGIN_MOMENT_Z must be the same!!" << endl;
    exit(1);
  }
  
  if (RefOriginMoment_X == NULL) {
    RefOriginMoment_X = new double[nMarker_Monitoring];
    for (iMarker = 0; iMarker < nMarker_Monitoring; iMarker++ )
      RefOriginMoment_X[iMarker] = 0.0;
  } else {
    if (nRefOriginMoment_X == 1) {
      
      double aux_RefOriginMoment_X = RefOriginMoment_X[0];
      delete [] RefOriginMoment_X;
      RefOriginMoment_X = new double[nMarker_Monitoring];
      nRefOriginMoment_X = nMarker_Monitoring;
      
      for (iMarker = 0; iMarker < nMarker_Monitoring; iMarker++ )
        RefOriginMoment_X[iMarker] = aux_RefOriginMoment_X;
    }
    else if (nRefOriginMoment_X != nMarker_Monitoring) {
      cout << "ERROR: Length of REF_ORIGIN_MOMENT_X must match number of Monitoring Markers!!" << endl;
      exit(1);
    }
  }
  
  if (RefOriginMoment_Y == NULL) {
    RefOriginMoment_Y = new double[nMarker_Monitoring];
    for (iMarker = 0; iMarker < nMarker_Monitoring; iMarker++ )
      RefOriginMoment_Y[iMarker] = 0.0;
  } else {
    if (nRefOriginMoment_Y == 1) {
      
      double aux_RefOriginMoment_Y = RefOriginMoment_Y[0];
      delete [] RefOriginMoment_Y;
      RefOriginMoment_Y = new double[nMarker_Monitoring];
      nRefOriginMoment_Y = nMarker_Monitoring;
      
      for (iMarker = 0; iMarker < nMarker_Monitoring; iMarker++ )
        RefOriginMoment_Y[iMarker] = aux_RefOriginMoment_Y;
    }
    else if (nRefOriginMoment_Y != nMarker_Monitoring) {
      cout << "ERROR: Length of REF_ORIGIN_MOMENT_Y must match number of Monitoring Markers!!" << endl;
      exit(1);
    }
  }
  
  if (RefOriginMoment_Z == NULL) {
    RefOriginMoment_Z = new double[nMarker_Monitoring];
    for (iMarker = 0; iMarker < nMarker_Monitoring; iMarker++ )
      RefOriginMoment_Z[iMarker] = 0.0;
  } else {
    if (nRefOriginMoment_Z == 1) {
      
      double aux_RefOriginMoment_Z = RefOriginMoment_Z[0];
      delete [] RefOriginMoment_Z;
      RefOriginMoment_Z = new double[nMarker_Monitoring];
      nRefOriginMoment_Z = nMarker_Monitoring;
      
      for (iMarker = 0; iMarker < nMarker_Monitoring; iMarker++ )
        RefOriginMoment_Z[iMarker] = aux_RefOriginMoment_Z;
    }
    else if (nRefOriginMoment_Z != nMarker_Monitoring) {
      cout << "ERROR: Length of REF_ORIGIN_MOMENT_Z must match number of Monitoring Markers!!" << endl;
      exit(1);
    }
  }
  
  if (FullMG) FinestMesh = nMultiLevel;
  else FinestMesh = MESH_0;
  
  if ((Kind_Solver == NAVIER_STOKES) &&
      (Kind_Turb_Model != NONE))
    Kind_Solver = RANS;
  
  Kappa_1st_Flow = Kappa_Flow[0];
  Kappa_2nd_Flow = Kappa_Flow[1];
  Kappa_4th_Flow = Kappa_Flow[2];
  
  // make the MG_PreSmooth, MG_PostSmooth, and MG_CorrecSmooth arrays consistent with nMultiLevel
  unsigned short * tmp_smooth = new unsigned short[nMultiLevel+1];
  
  if ((nMG_PreSmooth != nMultiLevel+1) && (nMG_PreSmooth != 0)) {
    if (nMG_PreSmooth > nMultiLevel+1) {
      
      // truncate by removing unnecessary elements at the end
      for (unsigned int i = 0; i <= nMultiLevel; i++)
        tmp_smooth[i] = MG_PreSmooth[i];
      delete [] MG_PreSmooth;
    } else {
      
      // add additional elements equal to last element
      for (unsigned int i = 0; i < nMG_PreSmooth; i++)
        tmp_smooth[i] = MG_PreSmooth[i];
      for (unsigned int i = nMG_PreSmooth; i <= nMultiLevel; i++)
        tmp_smooth[i] = MG_PreSmooth[nMG_PreSmooth-1];
      delete [] MG_PreSmooth;
    }
    
    nMG_PreSmooth = nMultiLevel+1;
    MG_PreSmooth = new unsigned short[nMG_PreSmooth];
    for (unsigned int i = 0; i < nMG_PreSmooth; i++)
      MG_PreSmooth[i] = tmp_smooth[i];
  }
  if ((nMultiLevel != 0) && (nMG_PreSmooth == 0)) {
    delete [] MG_PreSmooth;
    nMG_PreSmooth = nMultiLevel+1;
    MG_PreSmooth = new unsigned short[nMG_PreSmooth];
    for (unsigned int i = 0; i < nMG_PreSmooth; i++)
      MG_PreSmooth[i] = i+1;
  }
  
  if ((nMG_PostSmooth != nMultiLevel+1) && (nMG_PostSmooth != 0)) {
    if (nMG_PostSmooth > nMultiLevel+1) {
      // truncate by removing unnecessary elements at the end
      for (unsigned int i = 0; i <= nMultiLevel; i++)
        tmp_smooth[i] = MG_PostSmooth[i];
      delete [] MG_PostSmooth;
    } else {
      // add additional elements equal to last element
      for (unsigned int i = 0; i < nMG_PostSmooth; i++)
        tmp_smooth[i] = MG_PostSmooth[i];
      for (unsigned int i = nMG_PostSmooth; i <= nMultiLevel; i++)
        tmp_smooth[i] = MG_PostSmooth[nMG_PostSmooth-1];
      delete [] MG_PostSmooth;
    }
    nMG_PostSmooth = nMultiLevel+1;
    MG_PostSmooth = new unsigned short[nMG_PostSmooth];
    for (unsigned int i = 0; i < nMG_PostSmooth; i++)
      MG_PostSmooth[i] = tmp_smooth[i];
  }
  if ((nMultiLevel != 0) && (nMG_PostSmooth == 0)) {
    delete [] MG_PostSmooth;
    nMG_PostSmooth = nMultiLevel+1;
    MG_PostSmooth = new unsigned short[nMG_PostSmooth];
    for (unsigned int i = 0; i < nMG_PostSmooth; i++)
      MG_PostSmooth[i] = 0;
  }
  
  if ((nMG_CorrecSmooth != nMultiLevel+1) && (nMG_CorrecSmooth != 0)) {
    if (nMG_CorrecSmooth > nMultiLevel+1) {
      // truncate by removing unnecessary elements at the end
      for (unsigned int i = 0; i <= nMultiLevel; i++)
        tmp_smooth[i] = MG_CorrecSmooth[i];
      delete [] MG_CorrecSmooth;
    } else {
      // add additional elements equal to last element
      for (unsigned int i = 0; i < nMG_CorrecSmooth; i++)
        tmp_smooth[i] = MG_CorrecSmooth[i];
      for (unsigned int i = nMG_CorrecSmooth; i <= nMultiLevel; i++)
        tmp_smooth[i] = MG_CorrecSmooth[nMG_CorrecSmooth-1];
      delete [] MG_CorrecSmooth;
    }
    nMG_CorrecSmooth = nMultiLevel+1;
    MG_CorrecSmooth = new unsigned short[nMG_CorrecSmooth];
    for (unsigned int i = 0; i < nMG_CorrecSmooth; i++)
      MG_CorrecSmooth[i] = tmp_smooth[i];
  }
  if ((nMultiLevel != 0) && (nMG_CorrecSmooth == 0)) {
    delete [] MG_CorrecSmooth;
    nMG_CorrecSmooth = nMultiLevel+1;
    MG_CorrecSmooth = new unsigned short[nMG_CorrecSmooth];
    for (unsigned int i = 0; i < nMG_CorrecSmooth; i++)
      MG_CorrecSmooth[i] = 0;
  }
  
  
  delete [] tmp_smooth;
  
  // override MG Smooth parameters
  if (nMG_PreSmooth != 0)
    MG_PreSmooth[MESH_0] = 1;
  if (nMG_PostSmooth != 0) {
    MG_PostSmooth[MESH_0] = 0;
    MG_PostSmooth[nMultiLevel] = 0;
  }
  if (nMG_CorrecSmooth != 0)
    MG_CorrecSmooth[nMultiLevel] = 0;
  
  if (Restart) FullMG = false;
  
  if (Unsteady_Simulation == TIME_STEPPING) {
    nMultiLevel = 0;
    MGCycle = 0;
  }
  
  if (Unsteady_Simulation == TIME_SPECTRAL) {
    Kind_SourNumScheme_Flow = PIECEWISE_CONSTANT;
    Kind_SourNumScheme_AdjFlow  = PIECEWISE_CONSTANT;
  }
  
  if (Rotating_Frame == YES) {
    Kind_SourNumScheme_Flow = PIECEWISE_CONSTANT;
    Kind_SourNumScheme_AdjFlow  = PIECEWISE_CONSTANT;
  }
  
  if (Axisymmetric == YES) {
    Kind_SourNumScheme_Flow = PIECEWISE_CONSTANT;
    Kind_SourNumScheme_AdjFlow  = PIECEWISE_CONSTANT;
  }
  
  if (GravityForce == YES) {
    Kind_SourNumScheme_Flow = PIECEWISE_CONSTANT;
    Kind_SourNumScheme_AdjFlow  = PIECEWISE_CONSTANT;
  }
  
  nCFL = nMultiLevel+1;
  CFL = new double[nCFL];
  CFL[0] = CFLFineGrid;
  if (Adjoint) CFL[0] = CFL[0] * Adj_CFLRedCoeff;
  for (unsigned short iCFL = 1; iCFL < nCFL; iCFL++)
    CFL[iCFL] = CFL[iCFL-1]*MG_CFLRedCoeff;
  
  if (nRKStep == 0) {
    RK_Alpha_Step = new double[1]; RK_Alpha_Step[0] = 1.0;
  }
  
  if (Kind_Solver == NO_SOLVER) {
    cout << "You must define a solver type!!" << endl;
    exit(1);
  }
  
  if (((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS))
      && (Kind_ViscNumScheme_Flow == NONE)) {
    cout << "You must define a viscous numerical method for the flow equations!!" << endl;
    exit(1);
  }
  
  /*--- Set a flag for viscous simulations ---*/
  Viscous = (( Kind_Solver == NAVIER_STOKES          ) ||
             ( Kind_Solver == RANS                   )   );
  
}

void CConfig::SetMarkers(unsigned short val_software) {
  
  nDomain = SINGLE_NODE;
  
  /*--- Boundary (marker) treatment ---*/
  nMarker_All = nMarker_Euler + nMarker_FarField + nMarker_SymWall + nMarker_PerBound + nMarker_NearFieldBound + nMarker_Supersonic_Inlet
  + nMarker_InterfaceBound + nMarker_Dirichlet + nMarker_Neumann + nMarker_Inlet + nMarker_Outlet + nMarker_Isothermal + nMarker_HeatFlux
  + nMarker_NacelleInflow + nMarker_NacelleExhaust + nMarker_Dirichlet_Elec + nMarker_Displacement + nMarker_Load
  + nMarker_FlowLoad + nMarker_Pressure + nMarker_Custom + 2*nDomain;
  
  Marker_All_Tag        = new string[nMarker_All+2];			    // Store the tag that correspond with each marker.
  Marker_All_SendRecv   = new short[nMarker_All+2];						// +#domain (send), -#domain (receive) or 0 (neither send nor receive).
  Marker_All_Boundary   = new unsigned short[nMarker_All+2];	// Store the kind of boundary condition.
  Marker_All_Monitoring = new unsigned short[nMarker_All+2];	// Store whether the boundary should be monitored.
  Marker_All_Designing  = new unsigned short[nMarker_All+2];  // Store whether the boundary should be designed.
  Marker_All_Plotting   = new unsigned short[nMarker_All+2];	// Store whether the boundary should be plotted.
  Marker_All_DV         = new unsigned short[nMarker_All+2];	// Store whether the boundary should be affected by design variables.
  Marker_All_Moving     = new unsigned short[nMarker_All+2];	// Store whether the boundary should be in motion.
  Marker_All_PerBound   = new short[nMarker_All+2];						// Store whether the boundary belongs to a periodic boundary.
  
  unsigned short iMarker_All, iMarker_Config, iMarker_Euler, iMarker_Custom, iMarker_FarField,
  iMarker_SymWall, iMarker_Pressure, iMarker_PerBound, iMarker_NearFieldBound, iMarker_InterfaceBound, iMarker_Dirichlet,
  iMarker_Inlet, iMarker_Outlet, iMarker_Isothermal, iMarker_HeatFlux, iMarker_NacelleInflow, iMarker_NacelleExhaust, iMarker_Displacement, iMarker_Load,
  iMarker_FlowLoad, iMarker_Neumann, iMarker_Monitoring, iMarker_Designing, iMarker_Plotting, iMarker_DV, iMarker_Moving,
  iMarker_Supersonic_Inlet;
  
  for (iMarker_All = 0; iMarker_All < nMarker_All; iMarker_All++) {
    Marker_All_Tag[iMarker_All] = "NONE";
    Marker_All_SendRecv[iMarker_All]   = 0;
    Marker_All_Boundary[iMarker_All]   = 0;
    Marker_All_Monitoring[iMarker_All] = 0;
    Marker_All_Designing[iMarker_All]  = 0;
    Marker_All_Plotting[iMarker_All]   = 0;
    Marker_All_DV[iMarker_All]         = 0;
    Marker_All_Moving[iMarker_All]     = 0;
    Marker_All_PerBound[iMarker_All]   = 0;
  }
  
  nMarker_Config = nMarker_Euler + nMarker_FarField + nMarker_SymWall + nMarker_Pressure + nMarker_PerBound + nMarker_NearFieldBound
  + nMarker_InterfaceBound + nMarker_Dirichlet + nMarker_Neumann + nMarker_Inlet + nMarker_Outlet + nMarker_Isothermal + nMarker_HeatFlux + nMarker_NacelleInflow + nMarker_NacelleExhaust + nMarker_Supersonic_Inlet + nMarker_Displacement + nMarker_Load + nMarker_FlowLoad + nMarker_Custom;
  
  Marker_Config_Tag        = new string[nMarker_Config];
  Marker_Config_Boundary   = new unsigned short[nMarker_Config];
  Marker_Config_Monitoring = new unsigned short[nMarker_Config];
  Marker_Config_Plotting   = new unsigned short[nMarker_Config];
  Marker_Config_DV         = new unsigned short[nMarker_Config];
  Marker_Config_Moving     = new unsigned short[nMarker_Config];
  Marker_Config_Designing  = new unsigned short[nMarker_Config];
  Marker_Config_PerBound   = new unsigned short[nMarker_Config];
  
  for (iMarker_Config = 0; iMarker_Config < nMarker_Config; iMarker_Config++) {
    Marker_Config_Tag[iMarker_Config] = "NONE";
    Marker_Config_Boundary[iMarker_Config]   = 0;
    Marker_Config_Monitoring[iMarker_Config] = 0;
    Marker_Config_Designing[iMarker_Config]  = 0;
    Marker_Config_Plotting[iMarker_Config]   = 0;
    Marker_Config_DV[iMarker_Config]         = 0;
    Marker_Config_Moving[iMarker_Config]     = 0;
    Marker_Config_PerBound[iMarker_Config]   = 0;
  }
  
  iMarker_Config = 0;
  for (iMarker_Euler = 0; iMarker_Euler < nMarker_Euler; iMarker_Euler++) {
    Marker_Config_Tag[iMarker_Config] = Marker_Euler[iMarker_Euler];
    Marker_Config_Boundary[iMarker_Config] = EULER_WALL;
    iMarker_Config++;
  }
  
  for (iMarker_FarField = 0; iMarker_FarField < nMarker_FarField; iMarker_FarField++) {
    Marker_Config_Tag[iMarker_Config] = Marker_FarField[iMarker_FarField];
    Marker_Config_Boundary[iMarker_Config] = FAR_FIELD;
    iMarker_Config++;
  }
  
  for (iMarker_SymWall = 0; iMarker_SymWall < nMarker_SymWall; iMarker_SymWall++) {
    Marker_Config_Tag[iMarker_Config] = Marker_SymWall[iMarker_SymWall];
    Marker_Config_Boundary[iMarker_Config] = SYMMETRY_PLANE;
    iMarker_Config++;
  }
  
  for (iMarker_Pressure = 0; iMarker_Pressure < nMarker_Pressure; iMarker_Pressure++) {
    Marker_Config_Tag[iMarker_Config] = Marker_Pressure[iMarker_Pressure];
    Marker_Config_Boundary[iMarker_Config] = PRESSURE_BOUNDARY;
    iMarker_Config++;
  }
  
  for (iMarker_PerBound = 0; iMarker_PerBound < nMarker_PerBound; iMarker_PerBound++) {
    Marker_Config_Tag[iMarker_Config] = Marker_PerBound[iMarker_PerBound];
    Marker_Config_Boundary[iMarker_Config] = PERIODIC_BOUNDARY;
    Marker_Config_PerBound[iMarker_Config] = iMarker_PerBound + 1;
    iMarker_Config++;
  }
  
  for (iMarker_NearFieldBound = 0; iMarker_NearFieldBound < nMarker_NearFieldBound; iMarker_NearFieldBound++) {
    Marker_Config_Tag[iMarker_Config] = Marker_NearFieldBound[iMarker_NearFieldBound];
    Marker_Config_Boundary[iMarker_Config] = NEARFIELD_BOUNDARY;
    iMarker_Config++;
  }
  
  for (iMarker_InterfaceBound = 0; iMarker_InterfaceBound < nMarker_InterfaceBound; iMarker_InterfaceBound++) {
    Marker_Config_Tag[iMarker_Config] = Marker_InterfaceBound[iMarker_InterfaceBound];
    Marker_Config_Boundary[iMarker_Config] = INTERFACE_BOUNDARY;
    iMarker_Config++;
  }
  
  for (iMarker_Dirichlet = 0; iMarker_Dirichlet < nMarker_Dirichlet; iMarker_Dirichlet++) {
    Marker_Config_Tag[iMarker_Config] = Marker_Dirichlet[iMarker_Dirichlet];
    Marker_Config_Boundary[iMarker_Config] = DIRICHLET;
    iMarker_Config++;
  }
  
  for (iMarker_Inlet = 0; iMarker_Inlet < nMarker_Inlet; iMarker_Inlet++) {
    Marker_Config_Tag[iMarker_Config] = Marker_Inlet[iMarker_Inlet];
    Marker_Config_Boundary[iMarker_Config] = INLET_FLOW;
    iMarker_Config++;
  }
  
  FanFace_Mach = new double[nMarker_NacelleInflow];
  FanFace_Pressure = new double[nMarker_NacelleInflow];
  
  for (iMarker_NacelleInflow = 0; iMarker_NacelleInflow < nMarker_NacelleInflow; iMarker_NacelleInflow++) {
    Marker_Config_Tag[iMarker_Config] = Marker_NacelleInflow[iMarker_NacelleInflow];
    Marker_Config_Boundary[iMarker_Config] = NACELLE_INFLOW;
    FanFace_Mach[iMarker_NacelleInflow] = 0.0;
    FanFace_Pressure[iMarker_NacelleInflow] = 0.0;
    iMarker_Config++;
  }
  
  for (iMarker_NacelleExhaust = 0; iMarker_NacelleExhaust < nMarker_NacelleExhaust; iMarker_NacelleExhaust++) {
    Marker_Config_Tag[iMarker_Config] = Marker_NacelleExhaust[iMarker_NacelleExhaust];
    Marker_Config_Boundary[iMarker_Config] = NACELLE_EXHAUST;
    iMarker_Config++;
  }
  
  for (iMarker_Supersonic_Inlet = 0; iMarker_Supersonic_Inlet < nMarker_Supersonic_Inlet; iMarker_Supersonic_Inlet++) {
    Marker_Config_Tag[iMarker_Config] = Marker_Supersonic_Inlet[iMarker_Supersonic_Inlet];
    Marker_Config_Boundary[iMarker_Config] = SUPERSONIC_INLET;
    iMarker_Config++;
  }
  
  for (iMarker_Neumann = 0; iMarker_Neumann < nMarker_Neumann; iMarker_Neumann++) {
    Marker_Config_Tag[iMarker_Config] = Marker_Neumann[iMarker_Neumann];
    Marker_Config_Boundary[iMarker_Config] = NEUMANN;
    iMarker_Config++;
  }
  
  for (iMarker_Custom = 0; iMarker_Custom < nMarker_Custom; iMarker_Custom++) {
    Marker_Config_Tag[iMarker_Config] = Marker_Custom[iMarker_Custom];
    Marker_Config_Boundary[iMarker_Config] = CUSTOM_BOUNDARY;
    iMarker_Config++;
  }
  
  for (iMarker_Outlet = 0; iMarker_Outlet < nMarker_Outlet; iMarker_Outlet++) {
    Marker_Config_Tag[iMarker_Config] = Marker_Outlet[iMarker_Outlet];
    Marker_Config_Boundary[iMarker_Config] = OUTLET_FLOW;
    iMarker_Config++;
  }
  
  for (iMarker_Isothermal = 0; iMarker_Isothermal < nMarker_Isothermal; iMarker_Isothermal++) {
    Marker_Config_Tag[iMarker_Config] = Marker_Isothermal[iMarker_Isothermal];
    Marker_Config_Boundary[iMarker_Config] = ISOTHERMAL;
    iMarker_Config++;
  }
  
  for (iMarker_HeatFlux = 0; iMarker_HeatFlux < nMarker_HeatFlux; iMarker_HeatFlux++) {
    Marker_Config_Tag[iMarker_Config] = Marker_HeatFlux[iMarker_HeatFlux];
    Marker_Config_Boundary[iMarker_Config] = HEAT_FLUX;
    iMarker_Config++;
  }
  
  for (iMarker_Displacement = 0; iMarker_Displacement < nMarker_Displacement; iMarker_Displacement++) {
    Marker_Config_Tag[iMarker_Config] = Marker_Displacement[iMarker_Displacement];
    Marker_Config_Boundary[iMarker_Config] = DISPLACEMENT_BOUNDARY;
    iMarker_Config++;
  }
  
  for (iMarker_Load = 0; iMarker_Load < nMarker_Load; iMarker_Load++) {
    Marker_Config_Tag[iMarker_Config] = Marker_Load[iMarker_Load];
    Marker_Config_Boundary[iMarker_Config] = LOAD_BOUNDARY;
    iMarker_Config++;
  }
  
  for (iMarker_FlowLoad = 0; iMarker_FlowLoad < nMarker_FlowLoad; iMarker_FlowLoad++) {
    Marker_Config_Tag[iMarker_Config] = Marker_FlowLoad[iMarker_FlowLoad];
    Marker_Config_Boundary[iMarker_Config] = FLOWLOAD_BOUNDARY;
    iMarker_Config++;
  }
  
  for (iMarker_Config = 0; iMarker_Config < nMarker_Config; iMarker_Config++) {
    Marker_Config_Monitoring[iMarker_Config] = NO;
    for (iMarker_Monitoring = 0; iMarker_Monitoring < nMarker_Monitoring; iMarker_Monitoring++)
      if (Marker_Config_Tag[iMarker_Config] == Marker_Monitoring[iMarker_Monitoring])
        Marker_Config_Monitoring[iMarker_Config] = YES;
  }
  
  for (iMarker_Config = 0; iMarker_Config < nMarker_Config; iMarker_Config++) {
    Marker_Config_Designing[iMarker_Config] = NO;
    for (iMarker_Designing = 0; iMarker_Designing < nMarker_Designing; iMarker_Designing++)
      if (Marker_Config_Tag[iMarker_Config] == Marker_Designing[iMarker_Designing])
        Marker_Config_Designing[iMarker_Config] = YES;
  }
  
  for (iMarker_Config = 0; iMarker_Config < nMarker_Config; iMarker_Config++) {
    Marker_Config_Plotting[iMarker_Config] = NO;
    for (iMarker_Plotting = 0; iMarker_Plotting < nMarker_Plotting; iMarker_Plotting++)
      if (Marker_Config_Tag[iMarker_Config] == Marker_Plotting[iMarker_Plotting])
        Marker_Config_Plotting[iMarker_Config] = YES;
  }
  
  for (iMarker_Config = 0; iMarker_Config < nMarker_Config; iMarker_Config++) {
    Marker_Config_DV[iMarker_Config] = NO;
    for (iMarker_DV = 0; iMarker_DV < nMarker_DV; iMarker_DV++)
      if (Marker_Config_Tag[iMarker_Config] == Marker_DV[iMarker_DV])
        Marker_Config_DV[iMarker_Config] = YES;
  }
  
  for (iMarker_Config = 0; iMarker_Config < nMarker_Config; iMarker_Config++) {
    Marker_Config_Moving[iMarker_Config] = NO;
    for (iMarker_Moving = 0; iMarker_Moving < nMarker_Moving; iMarker_Moving++)
      if (Marker_Config_Tag[iMarker_Config] == Marker_Moving[iMarker_Moving])
        Marker_Config_Moving[iMarker_Config] = YES;
  }
  
}

void CConfig::SetOutput(unsigned short val_software) {
  unsigned short iMarker_Euler, iMarker_Custom, iMarker_FarField,
  iMarker_SymWall, iMarker_PerBound, iMarker_Pressure, iMarker_NearFieldBound, iMarker_InterfaceBound, iMarker_Dirichlet,
  iMarker_Inlet, iMarker_Outlet, iMarker_Isothermal, iMarker_HeatFlux, iMarker_NacelleInflow, iMarker_NacelleExhaust, iMarker_Displacement, iMarker_Load, iMarker_FlowLoad,  iMarker_Neumann, iMarker_Monitoring, iMarker_Designing, iMarker_Plotting, iMarker_DV, iMarker_Moving, iMarker_Supersonic_Inlet;
  
  cout << endl <<"------------------------ Physical case definition -----------------------" << endl;
  if ((val_software == SU2_CFD) || (val_software == SU2_EDU)) {
    switch (Kind_Solver) {
      case EULER: cout << "Compressible Euler equations." << endl; break;
      case NAVIER_STOKES:
        cout << "Compressible Laminar Navier-Stokes' equations." << endl; break;
      case RANS:
        cout << "Compressible RANS equations." << endl;
        cout << "Turbulence model: ";
        switch (Kind_Turb_Model) {
          case SA:  cout << "Spalart Allmaras" << endl; break;
          case SST: cout << "Menter's SST"     << endl; break;
        }
        break;
    }
    
    cout << "Mach number: " << Mach <<"."<< endl;
    cout << "Angle of attack (AoA): " << AoA <<" deg, and angle of sideslip (AoS): " << AoS <<" deg."<< endl;
    if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS))
      cout << "Reynolds number: " << Reynolds <<"."<< endl;
    
    if (Restart) { cout << "Read flow solution from: " << Solution_FlowFileName << "." << endl; }
    else { cout << "No restart solution, use the values at infinity (freestream)." << endl; }
    
    if (RefAreaCoeff == 0) cout << "The reference length/area will be computed using y(2D) or z(3D) projection." <<endl;
    else cout << "The reference length/area (force coefficient) is " << RefAreaCoeff << "." <<endl;
    cout << "The reference length (moment computation) is " << RefLengthMoment << "." <<endl;
    
    if ((nRefOriginMoment_X > 1) || (nRefOriginMoment_Y > 1) || (nRefOriginMoment_Z > 1)) {
      cout << "Surface(s) where the force coefficients are evaluated and their reference origin for moment computation: ";
      for (iMarker_Monitoring = 0; iMarker_Monitoring < nMarker_Monitoring; iMarker_Monitoring++) {
        cout << Marker_Monitoring[iMarker_Monitoring] << " (" << RefOriginMoment_X[iMarker_Monitoring] <<", "<<RefOriginMoment_Y[iMarker_Monitoring] <<", "<< RefOriginMoment_Z[iMarker_Monitoring] << ")";
        if (iMarker_Monitoring < nMarker_Monitoring-1) cout << ", ";
        else cout <<"."<<endl;
      }
    }
    else {
      cout << "Reference origin (moment computation) is (" << RefOriginMoment_X[0] << ", " << RefOriginMoment_Y[0] << ", " << RefOriginMoment_Z[0] << ")." << endl;
      cout << "Surface(s) where the force coefficients are evaluated: ";
      for (iMarker_Monitoring = 0; iMarker_Monitoring < nMarker_Monitoring; iMarker_Monitoring++) {
        cout << Marker_Monitoring[iMarker_Monitoring];
        if (iMarker_Monitoring < nMarker_Monitoring-1) cout << ", ";
        else cout <<"."<<endl;
      }
    }
    
    if (nMarker_Designing != 0) {
      cout << "Surface(s) where the objective function is evaluated: ";
      for (iMarker_Designing = 0; iMarker_Designing < nMarker_Designing; iMarker_Designing++) {
        cout << Marker_Designing[iMarker_Designing];
        if (iMarker_Designing < nMarker_Designing-1) cout << ", ";
        else cout <<".";
      }
      cout<<endl;
    }
    
    cout << "Surface(s) plotted in the output file: ";
    for (iMarker_Plotting = 0; iMarker_Plotting < nMarker_Plotting; iMarker_Plotting++) {
      cout << Marker_Plotting[iMarker_Plotting];
      if (iMarker_Plotting < nMarker_Plotting-1) cout << ", ";
      else cout <<".";
    }
    cout<<endl;
    
    if (nMarker_DV != 0) {
      cout << "Surface(s) affected by the design variables: ";
      for (iMarker_DV = 0; iMarker_DV < nMarker_DV; iMarker_DV++) {
        cout << Marker_DV[iMarker_DV];
        if (iMarker_DV < nMarker_DV-1) cout << ", ";
        else cout <<".";
      }
      cout<<endl;
    }
    
    if ((Kind_GridMovement[ZONE_0] == DEFORMING) || (Kind_GridMovement[ZONE_0] == MOVING_WALL)) {
      cout << "Surface(s) in motion: ";
      for (iMarker_Moving = 0; iMarker_Moving < nMarker_Moving; iMarker_Moving++) {
        cout << Marker_Moving[iMarker_Moving];
        if (iMarker_Moving < nMarker_Moving-1) cout << ", ";
        else cout <<".";
      }
      cout<<endl;
    }
    
  }
  
  cout << "Input mesh file name: " << Mesh_FileName << endl;
  
  if ((((val_software == SU2_CFD) || (val_software == SU2_EDU)) && ( Linearized )) || (val_software == SU2_GPC)) {
    cout << endl <<"-------------------- Surface deformation parameters ---------------------" << endl;
    cout << "Geo. design var. definition (markers <-> old def., new def. <-> param):" <<endl;
    for (unsigned short iDV = 0; iDV < nDV; iDV++) {
      switch (Design_Variable[iDV]) {
        case NO_DEFORMATION: cout << "There isn't any deformation." ; break;
        case HICKS_HENNE: cout << "Hicks Henne <-> " ; break;
        case COSINE_BUMP: cout << "Cosine bump <-> " ; break;
        case FOURIER: cout << "Fourier <-> " ; break;
        case SPHERICAL: cout << "Spherical design <-> " ; break;
        case MACH_NUMBER: cout << "Mach number <-> " ; break;
        case DISPLACEMENT: cout << "Displacement design variable."; break;
        case NACA_4DIGITS: cout << "NACA four digits <-> "; break;
        case PARABOLIC: cout << "Parabolic <-> "; break;
        case OBSTACLE: cout << "Obstacle <-> "; break;
        case AIRFOIL: cout << "Airfoil <-> "; break;
        case STRETCH: cout << "Stretch <-> "; break;
        case ROTATION: cout << "Rotation <-> "; break;
        case FFD_CONTROL_POINT: cout << "FFD (control point) <-> "; break;
        case FFD_DIHEDRAL_ANGLE: cout << "FFD (dihedral angle) <-> "; break;
        case FFD_TWIST_ANGLE: cout << "FFD (twist angle) <-> "; break;
        case FFD_ROTATION: cout << "FFD (rotation) <-> "; break;
        case FFD_CAMBER: cout << "FFD (camber) <-> "; break;
        case FFD_THICKNESS: cout << "FFD (thickness) <-> "; break;
        case FFD_VOLUME: cout << "FFD (volume) <-> "; break;
      }
      for (iMarker_DV = 0; iMarker_DV < nMarker_DV; iMarker_DV++) {
        cout << Marker_DV[iMarker_DV];
        if (iMarker_DV < nMarker_DV-1) cout << ", ";
        else cout << " <-> ";
      }
      cout << DV_Value[iDV] << " <-> ";
      
      if (Design_Variable[iDV] == NO_DEFORMATION) nParamDV = 0;
      if (Design_Variable[iDV] == HICKS_HENNE) nParamDV = 2;
      if (Design_Variable[iDV] == COSINE_BUMP) nParamDV = 3;
      if (Design_Variable[iDV] == FOURIER) nParamDV = 3;
      if (Design_Variable[iDV] == SPHERICAL) nParamDV = 3;
      if (Design_Variable[iDV] == DISPLACEMENT) nParamDV = 3;
      if (Design_Variable[iDV] == ROTATION) nParamDV = 6;
      if (Design_Variable[iDV] == NACA_4DIGITS) nParamDV = 3;
      if (Design_Variable[iDV] == PARABOLIC) nParamDV = 2;
      if (Design_Variable[iDV] == OBSTACLE) nParamDV = 2;
      if (Design_Variable[iDV] == AIRFOIL) nParamDV = 2;
      if (Design_Variable[iDV] == STRETCH) nParamDV = 2;
      if (Design_Variable[iDV] == FFD_CONTROL_POINT) nParamDV = 7;
      if (Design_Variable[iDV] == FFD_DIHEDRAL_ANGLE) nParamDV = 7;
      if (Design_Variable[iDV] == FFD_TWIST_ANGLE) nParamDV = 7;
      if (Design_Variable[iDV] == FFD_ROTATION) nParamDV = 7;
      if (Design_Variable[iDV] == FFD_CAMBER) nParamDV = 3;
      if (Design_Variable[iDV] == FFD_THICKNESS) nParamDV = 3;
      if (Design_Variable[iDV] == FFD_VOLUME) nParamDV = 3;
      
      for (unsigned short iParamDV = 0; iParamDV < nParamDV; iParamDV++) {
        if (iParamDV == 0) cout << "( ";
        cout << ParamDV[iDV][iParamDV];
        if (iParamDV < nParamDV-1) cout << ", ";
        else cout <<" )"<<endl;
      }
    }
  }
  
  
  if ((val_software == SU2_CFD) || (val_software == SU2_EDU)) {
    cout << endl <<"---------------------- Space numerical integration ----------------------" << endl;
    
    if (SmoothNumGrid) cout << "There are some smoothing iterations on the grid coordinates." <<endl;
    
    if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
      if ((Kind_ConvNumScheme_Flow == SPACE_CENTERED) && (Kind_Centered_Flow == JST)) {
        cout << "Jameson-Schmidt-Turkel scheme for the flow inviscid terms."<< endl;
        cout << "JST viscous coefficients (1st, 2nd & 4th): " << Kappa_1st_Flow
        << ", " << Kappa_2nd_Flow << ", " << Kappa_4th_Flow <<"."<< endl;
        cout << "The method includes a grid stretching correction (p = 0.3)."<< endl;
      }
      if ((Kind_ConvNumScheme_Flow == SPACE_CENTERED) && (Kind_Centered_Flow == LAX))
        cout << "Lax-Friedrich scheme for the flow inviscid terms."<< endl;
      if (Kind_ConvNumScheme_Flow == SPACE_UPWIND) {
        if (Kind_Upwind_Flow == ROE_1ST) cout << "1st order Roe solver for the flow inviscid terms."<< endl;
        if (Kind_Upwind_Flow == ROE_TURKEL_1ST) cout << "1st order Roe-Turkel solver for the flow inviscid terms."<< endl;
        if (Kind_Upwind_Flow == AUSM_1ST)	cout << "1st order AUSM solver for the flow inviscid terms."<< endl;
        if (Kind_Upwind_Flow == HLLC_1ST)	cout << "1st order HLLC solver for the flow inviscid terms."<< endl;
        if (Kind_Upwind_Flow == SW_1ST)	cout << "1st order Steger-Warming solver for the flow inviscid terms."<< endl;
        if (Kind_Upwind_Flow == MSW_1ST)	cout << "1st order Modified Steger-Warming solver for the flow inviscid terms."<< endl;
      }
      if ((Kind_ConvNumScheme_Flow == SPACE_UPWIND) &&
          ((Kind_Upwind_Flow == ROE_2ND) || (Kind_Upwind_Flow == AUSM_2ND) || (Kind_Upwind_Flow == HLLC_2ND)
           || (Kind_Upwind_Flow == SW_2ND) || (Kind_Upwind_Flow == MSW_2ND) || (Kind_Upwind_Flow == ROE_TURKEL_2ND))) {
            if (Kind_Upwind_Flow == ROE_2ND) cout << "2nd order Roe solver for the flow inviscid terms."<< endl;
            if (Kind_Upwind_Flow == ROE_TURKEL_2ND) cout << "2nd order Roe-Turkel solver for the flow inviscid terms."<< endl;
            if (Kind_Upwind_Flow == AUSM_2ND) cout << "2nd order AUSM solver for the flow inviscid terms."<< endl;
            if (Kind_Upwind_Flow == HLLC_2ND) cout << "2nd order HLLC solver for the flow inviscid terms."<< endl;
            if (Kind_Upwind_Flow == SW_2ND) cout << "2nd order Steger-Warming solver for the flow inviscid terms."<< endl;
            if (Kind_Upwind_Flow == MSW_2ND) cout << "2nd order Modified Steger-Warming solver for the flow inviscid terms."<< endl;
            switch (Kind_SlopeLimit_Flow) {
              case NONE: cout << "Without slope-limiting method." << endl; break;
              case VENKATAKRISHNAN:
                cout << "Venkatakrishnan slope-limiting method, with constant: " << LimiterCoeff <<". "<< endl;
                cout << "The reference element size is: " << RefElemLength <<". "<< endl;
                break;
              case MINMOD:
                cout << "Minmod slope-limiting method." << endl;
                break;
            }
          }
    }
    
    if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
      switch (Kind_ViscNumScheme_Flow) {
        case AVG_GRAD: cout << "Average of gradients (viscous flow terms)." << endl; break;
        case AVG_GRAD_CORRECTED: cout << "Average of gradients with correction (viscous flow terms)." << endl; break;
      }
    }
    
    if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
      if (Kind_SourNumScheme_Flow == PIECEWISE_CONSTANT) cout << "Piecewise constant integration of the flow source terms." << endl;
    }
    
    if (Kind_Solver == RANS) {
      if ((Kind_ConvNumScheme_Turb == SPACE_UPWIND) && (Kind_Upwind_Turb == SCALAR_UPWIND_1ST))
        cout << "Scalar upwind solver (first order) for the turbulence model."<< endl;
      if ((Kind_ConvNumScheme_Turb == SPACE_UPWIND) && (Kind_Upwind_Turb == SCALAR_UPWIND_2ND))
        cout << "Scalar upwind solver (second order) for the turbulence model."<< endl;
    }
    
    if (Kind_Solver == RANS) {
      if (Kind_ViscNumScheme_Turb == AVG_GRAD) cout << "Average of gradients (viscous turbulence terms)." << endl;
      if (Kind_ViscNumScheme_Turb == AVG_GRAD_CORRECTED) cout << "Average of gradients with correction (viscous turbulence terms)." << endl;
      if (Kind_SourNumScheme_Turb == PIECEWISE_CONSTANT) cout << "Piecewise constant integration of the turbulence model source terms." << endl;
    }
    
    switch (Kind_Gradient_Method) {
      case GREEN_GAUSS: cout << "Gradient computation using Green-Gauss theorem." << endl; break;
      case WEIGHTED_LEAST_SQUARES: cout << "Gradient Computation using weighted Least-Squares method." << endl; break;
    }
    
    cout << endl <<"---------------------- Time numerical integration -----------------------" << endl;
    
    if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
      switch (Kind_TimeIntScheme_Flow) {
        case RUNGE_KUTTA_EXPLICIT:
          cout << "Runge-Kutta explicit method for the flow equations." << endl;
          cout << "Number of steps: " << nRKStep << endl;
          cout << "Alpha coefficients: ";
          for (unsigned short iRKStep = 0; iRKStep < nRKStep; iRKStep++) {
            cout << "\t" << RK_Alpha_Step[iRKStep];
          }
          cout << endl;
          break;
        case EULER_EXPLICIT: cout << "Euler explicit method for the flow equations." << endl; break;
        case EULER_IMPLICIT:
          cout << "Euler implicit method for the flow equations." << endl;
          switch (Kind_Linear_Solver) {
            case BCGSTAB:
              cout << "BCGSTAB is used for solving the linear system." << endl;
              cout << "Convergence criteria of the linear solver: "<< Linear_Solver_Error <<"."<<endl;
              cout << "Max number of iterations: "<< Linear_Solver_Iter <<"."<<endl;
              cout << "Relaxation coefficient: "<< Linear_Solver_Relax <<"."<<endl;
              break;
            case FGMRES:
              cout << "FGMRES is used for solving the linear system." << endl;
              cout << "Convergence criteria of the linear solver: "<< Linear_Solver_Error <<"."<<endl;
              cout << "Max number of iterations: "<< Linear_Solver_Iter <<"."<<endl;
              cout << "Relaxation coefficient: "<< Linear_Solver_Relax <<"."<<endl;
              break;
          }
          break;
      }
    }
    
    if (nMultiLevel !=0) {
      if (nStartUpIter != 0) cout << "A total of " << nStartUpIter << " start up iterations on the fine grid."<< endl;
      if (MGCycle == 0) cout << "V Multigrid Cycle, with " << nMultiLevel << " multigrid levels."<< endl;
      if (MGCycle == 1) cout << "W Multigrid Cycle, with " << nMultiLevel << " multigrid levels."<< endl;
      
      cout << "Reduction of the CFL coefficient in the coarse levels: " << MG_CFLRedCoeff <<"."<<endl;
      cout << "Max. number of children in the agglomeration stage: " << MaxChildren <<"."<<endl;
      cout << "Max. length of an agglom. elem. (compared with the domain): " << MaxDimension <<"."<<endl;
      cout << "Damping factor for the residual restriction: " << Damp_Res_Restric <<"."<<endl;
      cout << "Damping factor for the correction prolongation: " << Damp_Correc_Prolong <<"."<<endl;
    }
    
    if (CFLRamp[0] == 1.0) cout << "No CFL ramp." << endl;
    else cout << "CFL ramp definition. factor: "<< CFLRamp[0] <<", every "<< int(CFLRamp[1]) <<" iterations, with a limit of "<< CFLRamp[2] <<"." << endl;
    
    if (nMultiLevel !=0) {
      cout << "Multigrid Level:                  ";
      for (unsigned short iLevel = 0; iLevel < nMultiLevel+1; iLevel++) {
        cout.width(6); cout << iLevel;
      }
      cout << endl;
    }
    
    cout << "Courant-Friedrichs-Lewy number:   ";
    cout.precision(3);
    for (unsigned short iCFL = 0; iCFL < nMultiLevel+1; iCFL++) {
      cout.width(6); cout << CFL[iCFL];
    }
    cout << endl;
    
    if (nMultiLevel !=0) {
      cout.precision(3);
      cout << "MG PreSmooth coefficients:        ";
      for (unsigned short iMG_PreSmooth = 0; iMG_PreSmooth < nMultiLevel+1; iMG_PreSmooth++) {
        cout.width(6); cout << MG_PreSmooth[iMG_PreSmooth];
      }
      cout << endl;
    }
    
    if (nMultiLevel !=0) {
      cout.precision(3);
      cout << "MG PostSmooth coefficients:       ";
      for (unsigned short iMG_PostSmooth = 0; iMG_PostSmooth < nMultiLevel+1; iMG_PostSmooth++) {
        cout.width(6); cout << MG_PostSmooth[iMG_PostSmooth];
      }
      cout << endl;
    }
    
    if (nMultiLevel !=0) {
      cout.precision(3);
      cout << "MG CorrecSmooth coefficients:     ";
      for (unsigned short iMG_CorrecSmooth = 0; iMG_CorrecSmooth < nMultiLevel+1; iMG_CorrecSmooth++) {
        cout.width(6); cout << MG_CorrecSmooth[iMG_CorrecSmooth];
      }
      cout << endl;
    }
    
    if (Kind_Solver == RANS)
      if (Kind_TimeIntScheme_Turb == EULER_IMPLICIT)
        cout << "Euler implicit time integration for the turbulence model." << endl;
  }
  
  if ((val_software == SU2_CFD) || (val_software == SU2_EDU))
    cout << endl <<"------------------------- Convergence criteria --------------------------" << endl;
  
  if ((val_software == SU2_CFD) || (val_software == SU2_EDU)) {
    cout << "Maximum number of iterations: " << nExtIter <<"."<<endl;
    
    if (ConvCriteria == CAUCHY) {
      if (!Adjoint && !Linearized)
        switch (Cauchy_Func_Flow) {
          case LIFT_COEFFICIENT: cout << "Cauchy criteria for Lift using "
            << Cauchy_Elems << " elements and epsilon " <<Cauchy_Eps<< "."<< endl; break;
          case DRAG_COEFFICIENT: cout << "Cauchy criteria for Drag using "
            << Cauchy_Elems << " elements and epsilon " <<Cauchy_Eps<< "."<< endl; break;
        }
      
      if (Adjoint)
        switch (Cauchy_Func_AdjFlow) {
          case SENS_GEOMETRY: cout << "Cauchy criteria for geo. sensitivity using "
            << Cauchy_Elems << " elements and epsilon " <<Cauchy_Eps<< "."<< endl; break;
          case SENS_MACH: cout << "Cauchy criteria for Mach number sensitivity using "
            << Cauchy_Elems << " elements and epsilon " <<Cauchy_Eps<< "."<< endl; break;
        }
      
      if (Linearized)
        switch (Cauchy_Func_LinFlow) {
          case DELTA_LIFT_COEFFICIENT: cout << "Cauchy criteria for linearized Lift using "
            << Cauchy_Elems << " elements and epsilon " <<Cauchy_Eps<< "."<< endl; break;
          case DELTA_DRAG_COEFFICIENT: cout << "Cauchy criteria for linearized Drag using "
            << Cauchy_Elems << " elements and epsilon " <<Cauchy_Eps<< "."<< endl; break;
        }
      
      cout << "Start convergence criteria at iteration " << StartConv_Iter<< "."<< endl;
      if (OneShot) cout << "Cauchy criteria for one shot method " << Cauchy_Eps_OneShot<< "."<< endl;
      if (FullMG) cout << "Cauchy criteria for full multigrid " << Cauchy_Eps_FullMG<< "."<< endl;
    }
    
    
    if (ConvCriteria == RESIDUAL) {
      if (!Adjoint && !Linearized) {
        cout << "Reduce the density residual " << OrderMagResidual << " orders of magnitude."<< endl;
        cout << "The minimum bound for the density residual is 10^(" << MinLogResidual<< ")."<< endl;
        cout << "Start convergence criteria at iteration " << StartConv_Iter<< "."<< endl;
      }
      
      if (Adjoint) {
        cout << "Reduce the adjoint density residual " << OrderMagResidual << " orders of magnitude."<< endl;
        cout << "The minimum value for the adjoint density residual is 10^(" << MinLogResidual<< ")."<< endl;
      }
      
      if (Linearized) {
        cout << "Reduce the linearized density residual " << OrderMagResidual << " orders of magnitude."<< endl;
        cout << "The minimum value for the linearized density residual is 10^(" << MinLogResidual<< ")."<< endl;
      }
      
    }
    
  }
  
  cout << endl <<"-------------------------- Output information ---------------------------" << endl;
  
  if ((Output_FileFormat==STL) && (val_software != SU2_MDC)){
    cerr << "Error: STL output file format only valid for SU2_MDC" << endl; throw(-1);
  }
  
  if ((val_software == SU2_CFD) || (val_software == SU2_EDU)) {
    
    cout << "Writing a flow solution every " << Wrt_Sol_Freq <<" iterations."<<endl;
    cout << "Writing the convergence history every " << Wrt_Con_Freq <<" iterations."<<endl;
    if ((Unsteady_Simulation == DT_STEPPING_1ST) || (Unsteady_Simulation == DT_STEPPING_2ND))  {
      cout << "Writing the dual time flow solution every " << Wrt_Sol_Freq_DualTime <<" iterations."<<endl;
      cout << "Writing the dual time convergence history every " << Wrt_Con_Freq_DualTime <<" iterations."<<endl;
    }
    
    switch (Output_FileFormat) {
      case PARAVIEW: cout << "The output file format is Paraview ASCII (.vtk)." << endl; break;
      case TECPLOT: cout << "The output file format is Tecplot ASCII (.dat)." << endl; break;
      case TECPLOT_BINARY: cout << "The output file format is Tecplot binary (.plt)." << endl; break;
      case CGNS_SOL: cout << "The output file format is CGNS (.cgns)." << endl; break;
    }
    
    cout << "Convergence history file name: " << Conv_FileName << "." << endl;
    
    cout << "Surface flow coefficients file name: " << SurfFlowCoeff_FileName << "." << endl;
    cout << "Flow variables file name: " << Flow_FileName << "." << endl;
    cout << "Restart flow file name: " << Restart_FlowFileName << "." << endl;
    
  }
  
  if (val_software == SU2_SOL) {
    switch (Output_FileFormat) {
      case PARAVIEW: cout << "The output file format is Paraview ASCII (.dat)." << endl; break;
      case TECPLOT: cout << "The output file format is Tecplot ASCII (.dat)." << endl; break;
      case TECPLOT_BINARY: cout << "The output file format is Tecplot binary (.plt)." << endl; break;
      case CGNS_SOL: cout << "The output file format is CGNS (.cgns)." << endl; break;
    }
    cout << "Flow variables file name: " << Flow_FileName << "." << endl;
  }
  
  if (val_software == SU2_MDC) {
    cout << "Output mesh file name: " << Mesh_Out_FileName << ". " << endl;
    if (Visualize_Deformation) cout << "A file will be created to visualize the deformation." << endl;
    else cout << "No file for visualizing the deformation." << endl;
  }
  
  if (val_software == SU2_PBC) {
    cout << "Output mesh file name: " << Mesh_Out_FileName << ". " << endl;
  }
  
  if (val_software == SU2_SMC) {
    cout << "Output mesh file name: " << Mesh_Out_FileName << ". " << endl;
  }
  
  if (val_software == SU2_GPC) {
    cout << "Output gradient file name: " << ObjFunc_Grad_FileName << ". " << endl;
  }
  
  if (val_software == SU2_MAC) {
    cout << "Output mesh file name: " << Mesh_Out_FileName << ". " << endl;
    cout << "Restart flow file name: " << Restart_FlowFileName << "." << endl;
    if ((Kind_Adaptation == FULL_ADJOINT) || (Kind_Adaptation == GRAD_ADJOINT) || (Kind_Adaptation == GRAD_FLOW_ADJ) ||
        (Kind_Adaptation == ROBUST) || (Kind_Adaptation == COMPUTABLE_ROBUST) || (Kind_Adaptation == COMPUTABLE) ||
        (Kind_Adaptation == REMAINING)) {
      if (Kind_ObjFunc == DRAG_COEFFICIENT) cout << "Restart adjoint file name: " << Restart_AdjFileName << "." << endl;
      if (Kind_ObjFunc == EQUIVALENT_AREA) cout << "Restart adjoint file name: " << Restart_AdjFileName << "." << endl;
      if (Kind_ObjFunc == NEARFIELD_PRESSURE) cout << "Restart adjoint file name: " << Restart_AdjFileName << "." << endl;
      if (Kind_ObjFunc == LIFT_COEFFICIENT) cout << "Restart adjoint file name: " << Restart_AdjFileName << "." << endl;
    }
  }
  
  if (val_software == SU2_DDC) {
    if (Visualize_Partition) cout << "Visualize the partitions. " << endl;
    else cout << "Don't visualize the partitions. " << endl;
  }
  
  cout << endl <<"------------------- Config file boundary information --------------------" << endl;
  
  if (nMarker_Euler != 0) {
    cout << "Euler wall boundary marker(s): ";
    for (iMarker_Euler = 0; iMarker_Euler < nMarker_Euler; iMarker_Euler++) {
      cout << Marker_Euler[iMarker_Euler];
      if (iMarker_Euler < nMarker_Euler-1) cout << ", ";
      else cout <<"."<<endl;
    }
  }
  
  if (nMarker_FarField != 0) {
    cout << "Far-field boundary marker(s): ";
    for (iMarker_FarField = 0; iMarker_FarField < nMarker_FarField; iMarker_FarField++) {
      cout << Marker_FarField[iMarker_FarField];
      if (iMarker_FarField < nMarker_FarField-1) cout << ", ";
      else cout <<"."<<endl;
    }
  }
  
  if (nMarker_SymWall != 0) {
    cout << "Symmetry plane boundary marker(s): ";
    for (iMarker_SymWall = 0; iMarker_SymWall < nMarker_SymWall; iMarker_SymWall++) {
      cout << Marker_SymWall[iMarker_SymWall];
      if (iMarker_SymWall < nMarker_SymWall-1) cout << ", ";
      else cout <<"."<<endl;
    }
  }
  
  if (nMarker_Pressure != 0) {
    cout << "Pressure boundary marker(s): ";
    for (iMarker_Pressure = 0; iMarker_Pressure < nMarker_Pressure; iMarker_Pressure++) {
      cout << Marker_Pressure[iMarker_Pressure];
      if (iMarker_Pressure < nMarker_Pressure-1) cout << ", ";
      else cout <<"."<<endl;
    }
  }
  
  if (nMarker_PerBound != 0) {
    cout << "Periodic boundary marker(s): ";
    for (iMarker_PerBound = 0; iMarker_PerBound < nMarker_PerBound; iMarker_PerBound++) {
      cout << Marker_PerBound[iMarker_PerBound];
      if (iMarker_PerBound < nMarker_PerBound-1) cout << ", ";
      else cout <<"."<<endl;
    }
  }
  
  if (nMarker_NearFieldBound != 0) {
    cout << "Near-field boundary marker(s): ";
    for (iMarker_NearFieldBound = 0; iMarker_NearFieldBound < nMarker_NearFieldBound; iMarker_NearFieldBound++) {
      cout << Marker_NearFieldBound[iMarker_NearFieldBound];
      if (iMarker_NearFieldBound < nMarker_NearFieldBound-1) cout << ", ";
      else cout <<"."<<endl;
    }
  }
  
  if (nMarker_InterfaceBound != 0) {
    cout << "Interface boundary marker(s): ";
    for (iMarker_InterfaceBound = 0; iMarker_InterfaceBound < nMarker_InterfaceBound; iMarker_InterfaceBound++) {
      cout << Marker_InterfaceBound[iMarker_InterfaceBound];
      if (iMarker_InterfaceBound < nMarker_InterfaceBound-1) cout << ", ";
      else cout <<"."<<endl;
    }
  }
  
  if (nMarker_Dirichlet != 0) {
    cout << "Dirichlet boundary marker(s): ";
    for (iMarker_Dirichlet = 0; iMarker_Dirichlet < nMarker_Dirichlet; iMarker_Dirichlet++) {
      cout << Marker_Dirichlet[iMarker_Dirichlet];
      if (iMarker_Dirichlet < nMarker_Dirichlet-1) cout << ", ";
      else cout <<"."<<endl;
    }
  }
  
  if (nMarker_FlowLoad != 0) {
    cout << "Flow Load boundary marker(s): ";
    for (iMarker_FlowLoad = 0; iMarker_FlowLoad < nMarker_FlowLoad; iMarker_FlowLoad++) {
      cout << Marker_FlowLoad[iMarker_FlowLoad];
      if (iMarker_FlowLoad < nMarker_FlowLoad-1) cout << ", ";
      else cout <<"."<<endl;
    }
  }
  
  if (nMarker_Neumann != 0) {
    cout << "Neumann boundary marker(s): ";
    for (iMarker_Neumann = 0; iMarker_Neumann < nMarker_Neumann; iMarker_Neumann++) {
      cout << Marker_Neumann[iMarker_Neumann];
      if (iMarker_Neumann < nMarker_Neumann-1) cout << ", ";
      else cout <<"."<<endl;
    }
  }
  
  if (nMarker_Inlet != 0) {
    cout << "Inlet boundary marker(s): ";
    for (iMarker_Inlet = 0; iMarker_Inlet < nMarker_Inlet; iMarker_Inlet++) {
      cout << Marker_Inlet[iMarker_Inlet];
      if (iMarker_Inlet < nMarker_Inlet-1) cout << ", ";
      else cout <<"."<<endl;
    }
  }
  
  if (nMarker_NacelleInflow != 0) {
    cout << "Nacelle inflow boundary marker(s): ";
    for (iMarker_NacelleInflow = 0; iMarker_NacelleInflow < nMarker_NacelleInflow; iMarker_NacelleInflow++) {
      cout << Marker_NacelleInflow[iMarker_NacelleInflow];
      if (iMarker_NacelleInflow < nMarker_NacelleInflow-1) cout << ", ";
      else cout <<"."<<endl;
    }
  }
  
  if (nMarker_NacelleExhaust != 0) {
    cout << "Nacelle exhaust boundary marker(s): ";
    for (iMarker_NacelleExhaust = 0; iMarker_NacelleExhaust < nMarker_NacelleExhaust; iMarker_NacelleExhaust++) {
      cout << Marker_NacelleExhaust[iMarker_NacelleExhaust];
      if (iMarker_NacelleExhaust < nMarker_NacelleExhaust-1) cout << ", ";
      else cout <<"."<<endl;
    }
  }
  
  if (nMarker_Supersonic_Inlet != 0) {
    cout << "Supersonic inlet boundary marker(s): ";
    for (iMarker_Supersonic_Inlet = 0; iMarker_Supersonic_Inlet < nMarker_Supersonic_Inlet; iMarker_Supersonic_Inlet++) {
      cout << Marker_Supersonic_Inlet[iMarker_Supersonic_Inlet];
      if (iMarker_Supersonic_Inlet < nMarker_Supersonic_Inlet-1) cout << ", ";
      else cout <<"."<<endl;
    }
  }
  
  if (nMarker_Outlet != 0) {
    cout << "Outlet boundary marker(s): ";
    for (iMarker_Outlet = 0; iMarker_Outlet < nMarker_Outlet; iMarker_Outlet++) {
      cout << Marker_Outlet[iMarker_Outlet];
      if (iMarker_Outlet < nMarker_Outlet-1) cout << ", ";
      else cout <<"."<<endl;
    }
  }
  
  if (nMarker_Isothermal != 0) {
    cout << "Isothermal wall boundary marker(s): ";
    for (iMarker_Isothermal = 0; iMarker_Isothermal < nMarker_Isothermal; iMarker_Isothermal++) {
      cout << Marker_Isothermal[iMarker_Isothermal];
      if (iMarker_Isothermal < nMarker_Isothermal-1) cout << ", ";
      else cout <<"."<<endl;
    }
  }
  
  if (nMarker_HeatFlux != 0) {
    cout << "Constant heat flux wall boundary marker(s): ";
    for (iMarker_HeatFlux = 0; iMarker_HeatFlux < nMarker_HeatFlux; iMarker_HeatFlux++) {
      cout << Marker_HeatFlux[iMarker_HeatFlux];
      if (iMarker_HeatFlux < nMarker_HeatFlux-1) cout << ", ";
      else cout <<"."<<endl;
    }
  }
  
  if (nMarker_Displacement != 0) {
    cout << "Displacement boundary marker(s): ";
    for (iMarker_Displacement = 0; iMarker_Displacement < nMarker_Displacement; iMarker_Displacement++) {
      cout << Marker_Displacement[iMarker_Displacement];
      if (iMarker_Displacement < nMarker_Displacement-1) cout << ", ";
      else cout <<"."<<endl;
    }
  }
  
  if (nMarker_Load != 0) {
    cout << "Load boundary marker(s): ";
    for (iMarker_Load = 0; iMarker_Load < nMarker_Load; iMarker_Load++) {
      cout << Marker_Load[iMarker_Load];
      if (iMarker_Load < nMarker_Load-1) cout << ", ";
      else cout <<"."<<endl;
    }
  }
  
  if (nMarker_Neumann != 0) {
    cout << "Neumann boundary marker(s): ";
    for (iMarker_Neumann = 0; iMarker_Neumann < nMarker_Neumann; iMarker_Neumann++) {
      cout << Marker_Neumann[iMarker_Neumann];
      if (iMarker_Neumann < nMarker_Neumann-1) cout << ", ";
      else cout <<"."<<endl;
    }
  }
  
  if (nMarker_Custom != 0) {
    cout << "Custom boundary marker(s): ";
    for (iMarker_Custom = 0; iMarker_Custom < nMarker_Custom; iMarker_Custom++) {
      cout << Marker_Custom[iMarker_Custom];
      if (iMarker_Custom < nMarker_Custom-1) cout << ", ";
      else cout <<"."<<endl;
    }
  }
  
}

void CConfig::AddMarkerOption(const string & name, unsigned short & num_marker, string* & marker) {
  //cout << "Adding Marker option " << name << endl;
  num_marker = 0;
  CAnyOptionRef* option_ref = new CMarkerOptionRef(marker, num_marker);
  param.insert( pair<string, CAnyOptionRef*>(name, option_ref) );
}

void CConfig::AddConvectOption(const string & name, unsigned short & space, unsigned short & centered,
                               unsigned short & upwind) {
  //cout << "Adding Convect option " << name << endl;
  centered = NO_CENTERED;
  upwind = NO_UPWIND;
  space = SPACE_CENTERED;
  CAnyOptionRef* option_ref = new CConvOptionRef(space, centered, upwind);
  param.insert( pair<string, CAnyOptionRef*>(name, option_ref) );
}

void CConfig::AddMathProblem(const string & name, bool & Adjoint, const bool & Adjoint_default,
                             bool & OneShot, const bool & OneShot_default,
                             bool & Linearized, const bool & Linearized_default,
                             bool & Restart_Flow, const bool & Restart_Flow_default) {
  //cout << "Adding Math Problem option " << name << endl;
  Adjoint = Adjoint_default;
  OneShot = OneShot_default;
  Linearized = Linearized_default;
  Restart_Flow = Restart_Flow_default;
  CAnyOptionRef* option_ref = new CMathProblemRef(Adjoint, OneShot, Linearized,
                                                  Restart_Flow);
  param.insert( pair<string, CAnyOptionRef*>(name, option_ref) );
}

void CConfig::AddDVParamOption(const string & name, unsigned short & nDV, double** & ParamDV,
                               unsigned short* & Design_Variable) {
  //cout << "Adding DV Param option " << name << endl;
  CAnyOptionRef* option_ref = new CDVParamOptionRef(nDV, ParamDV, Design_Variable);
  param.insert( pair<string, CAnyOptionRef*>(name, option_ref) );
}

void CConfig::AddMarkerPeriodic(const string & name, unsigned short & nMarker_PerBound,
                                string* & Marker_PerBound, string* & Marker_PerDonor,
                                double** & RotCenter, double** & RotAngles, double** & Translation) {
  //cout << "Adding Marker Periodic option " << name << endl;
  nMarker_PerBound = 0;
  CAnyOptionRef* option_ref = new CMarkerPeriodicRef(nMarker_PerBound, Marker_PerBound,
                                                     Marker_PerDonor, RotCenter,
                                                     RotAngles, Translation);
  param.insert( pair<string, CAnyOptionRef*>(name, option_ref) );
}

void CConfig::AddMarkerInlet(const string & name, unsigned short & nMarker_Inlet,
                             string* & Marker_Inlet, double* & Ttotal, double* & Ptotal,
                             double** & FlowDir) {
  nMarker_Inlet = 0;
  CAnyOptionRef* option_ref = new CMarkerInletRef(nMarker_Inlet, Marker_Inlet,
                                                  Ttotal, Ptotal, FlowDir);
  param.insert( pair<string, CAnyOptionRef*>(name, option_ref) );
}

void CConfig::AddMarkerInlet(const string & name, unsigned short & nMarker_Inlet,
                             string* & Marker_Inlet, double* & Ttotal, double* & Ptotal) {
  nMarker_Inlet = 0;
  CAnyOptionRef* option_ref = new CMarkerInletRef_(nMarker_Inlet, Marker_Inlet,
                                                   Ttotal, Ptotal);
  param.insert( pair<string, CAnyOptionRef*>(name, option_ref) );
}

void CConfig::AddMarkerDirichlet(const string & name, unsigned short & nMarker_Dirichlet_Elec,
                                 string* & Marker_Dirichlet_Elec, double* & Dirichlet_Value) {
  nMarker_Dirichlet_Elec = 0;
  CAnyOptionRef* option_ref = new CMarkerDirichletRef(nMarker_Dirichlet_Elec, Marker_Dirichlet_Elec,
                                                      Dirichlet_Value);
  param.insert( pair<string, CAnyOptionRef*>(name, option_ref) );
  
}

void CConfig::AddMarkerOutlet(const string & name, unsigned short & nMarker_Outlet,
                              string* & Marker_Outlet, double* & Pressure) {
  nMarker_Outlet = 0;
  CAnyOptionRef* option_ref = new CMarkerOutletRef(nMarker_Outlet, Marker_Outlet,
                                                   Pressure);
  param.insert( pair<string, CAnyOptionRef*>(name, option_ref) );
}

void CConfig::AddMarkerDisplacement(const string & name, unsigned short & nMarker_Displacement,
                                    string* & Marker_Displacement, double* & Displ) {
  nMarker_Displacement = 0;
  CAnyOptionRef* option_ref = new CMarkerDisplacementRef(nMarker_Displacement, Marker_Displacement, Displ);
  param.insert( pair<string, CAnyOptionRef*>(name, option_ref) );
}

void CConfig::AddMarkerLoad(const string & name, unsigned short & nMarker_Load,
                            string* & Marker_Load, double* & Force) {
  nMarker_Load = 0;
  CAnyOptionRef* option_ref = new CMarkerLoadRef(nMarker_Load, Marker_Load, Force);
  param.insert( pair<string, CAnyOptionRef*>(name, option_ref) );
}

void CConfig::AddMarkerFlowLoad(const string & name, unsigned short & nMarker_FlowLoad,
                                string* & Marker_FlowLoad, double* & FlowForce) {
  nMarker_FlowLoad = 0;
  CAnyOptionRef* option_ref = new CMarkerLoadRef(nMarker_FlowLoad, Marker_FlowLoad, FlowForce);
  param.insert( pair<string, CAnyOptionRef*>(name, option_ref) );
}

void CConfig::SetBoolOption(bool* ref, const vector<string> & value) {
  if ( (value[0] != "YES") && (value[0] != "NO") ) {
    cerr << "Error in CConfig::SetBoolOption(): "
    << "option value provided must be \"YES\" or \"NO\";"
    << "value given is " << value[0] << endl;
    throw(-1);
  }
  if (value[0] == "YES") {
    *ref = true;
  } else {
    *ref = false;
  }
}

bool CConfig::TokenizeString(string & str, string & option_name,
                             vector<string> & option_value) {
  const string delimiters(" ()[]{}:,\t\n\v\f\r");
  // check for comments or empty string
  string::size_type pos, last_pos;
  pos = str.find_first_of("%");
  if ( (str.length() == 0) || (pos == 0) ) {
    // str is empty or a comment line, so no option here
    return false;
  }
  if (pos != string::npos) {
    // remove comment at end if necessary
    str.erase(pos);
  }
  
  // look for line composed on only delimiters (usually whitespace)
  pos = str.find_first_not_of(delimiters);
  if (pos == string::npos) {
    return false;
  }
  
  // find the equals sign and split string
  string name_part, value_part;
  pos = str.find("=");
  if (pos == string::npos) {
    cerr << "Error in TokenizeString(): "
    << "line in the configuration file with no \"=\" sign."
    << endl;
    cout << "Look for: " << str << endl;
    cout << "str.length() = " << str.length() << endl;
    throw(-1);
  }
  name_part = str.substr(0, pos);
  value_part = str.substr(pos+1,string::npos);
  //cout << "name_part  = |" << name_part  << "|" << endl;
  //cout << "value_part = |" << value_part << "|" << endl;
  
  // the first_part should consist of one string with no interior delimiters
  last_pos = name_part.find_first_not_of(delimiters, 0);
  pos = name_part.find_first_of(delimiters, last_pos);
  if ( (name_part.length() == 0) || (last_pos == string::npos) ) {
    cerr << "Error in CConfig::TokenizeString(): "
    << "line in the configuration file with no name before the \"=\" sign."
    << endl;
    throw(-1);
  }
  if (pos == string::npos) pos = name_part.length();
  option_name = name_part.substr(last_pos, pos - last_pos);
  last_pos = name_part.find_first_not_of(delimiters, pos);
  if (last_pos != string::npos) {
    cerr << "Error in TokenizeString(): "
    << "two or more options before an \"=\" sign in the configuration file."
    << endl;
    throw(-1);
  }
  StringToUpperCase(option_name);
  
  //cout << "option_name = |" << option_name << "|" << endl;
  //cout << "pos = " << pos << ": last_pos = " << last_pos << endl;
  
  // now fill the option value vector
  option_value.clear();
  last_pos = value_part.find_first_not_of(delimiters, 0);
  pos = value_part.find_first_of(delimiters, last_pos);
  while (string::npos != pos || string::npos != last_pos) {
    // add token to the vector<string>
    option_value.push_back(value_part.substr(last_pos, pos - last_pos));
    // skip delimiters
    last_pos = value_part.find_first_not_of(delimiters, pos);
    // find next "non-delimiter"
    pos = value_part.find_first_of(delimiters, last_pos);
  }
  if (option_value.size() == 0) {
    cerr << "Error inT okenizeString(): "
    << "option " << option_name << " in configuration file with no value assigned."
    << endl;
    throw(-1);
  }
  
#if 0
  cout << "option value(s) = ";
  for (unsigned int i = 0; i < option_value.size(); i++)
    cout << option_value[i] << " ";
  cout << endl;
#endif
  
  // look for ';' DV delimiters attached to values
  vector<string>::iterator it;
  it = option_value.begin();
  while (it != option_value.end()) {
    if (it->compare(";") == 0) {
      it++;
      continue;
    }
    
    pos = it->find(';');
    if (pos != string::npos) {
      string before_semi = it->substr(0, pos);
      string after_semi= it->substr(pos+1,string::npos);
      if (before_semi.empty()) {
        *it = ";";
        it++;
        option_value.insert(it, after_semi);
      } else {
        *it = before_semi;
        it++;
        vector<string> to_insert;
        to_insert.push_back(";");
        if (!after_semi.empty())
          to_insert.push_back(after_semi);
        option_value.insert(it, to_insert.begin(), to_insert.end());
      }
      it = option_value.begin(); // go back to beginning; not efficient
      continue;
    } else {
      it++;
    }
  }
#if 0
  cout << "option value(s) = ";
  for (unsigned int i = 0; i < option_value.size(); i++)
    cout << option_value[i] << " ";
  cout << endl;
#endif
  // remove any consecutive ";"
  it = option_value.begin();
  bool semi_at_prev = false;
  while (it != option_value.end()) {
    if (semi_at_prev) {
      if (it->compare(";") == 0) {
        option_value.erase(it);
        it = option_value.begin();
        semi_at_prev = false;
        continue;
      }
    }
    if (it->compare(";") == 0) {
      semi_at_prev = true;
    } else {
      semi_at_prev = false;
    }
    it++;
  }
  
#if 0
  cout << "option value(s) = ";
  for (unsigned int i = 0; i < option_value.size(); i++)
    cout << option_value[i] << " ";
  cout << endl;
#endif
  return true;
}

bool CConfig::GetPython_Option(string & option_name) {
  
  bool isPython_Option = false;
  
  /*--- Check option name against all known Python options
   for a match. These are the design options that are
   never read by the SU2 C++ codes, and we would like
   to ignore them while processing the config file. ---*/
  if (option_name == "OBJFUNC")       isPython_Option = true;
  if (option_name == "OBJFUNC_SCALE")    isPython_Option = true;
  if (option_name == "CONST_IEQ")      isPython_Option = true;
  if (option_name == "CONST_IEQ_SCALE")      isPython_Option = true;
  if (option_name == "CONST_IEQ_SIGN")     isPython_Option = true;
  if (option_name == "CONST_IEQ_VALUE") isPython_Option = true;
  if (option_name == "CONST_EQ")      isPython_Option = true;
  if (option_name == "CONST_EQ_SCALE")      isPython_Option = true;
  if (option_name == "CONST_EQ_SIGN")     isPython_Option = true;
  if (option_name == "CONST_EQ_VALUE") isPython_Option = true;
  if (option_name == "DEFINITION_DV") isPython_Option = true;
  if (option_name == "TASKS") isPython_Option = true;
  if (option_name == "OPT_OBJECTIVE") isPython_Option = true;
  if (option_name == "OPT_CONSTRAINT") isPython_Option = true;
  if (option_name == "GRADIENTS") isPython_Option = true;
  if (option_name == "FIN_DIFF_STEP") isPython_Option = true;
  if (option_name == "ADAPT_CYCLES") isPython_Option = true;
  if (option_name == "CONSOLE") isPython_Option = true;
  if (option_name == "DECOMPOSED") isPython_Option = true;
  
  return isPython_Option;
}

unsigned short CConfig::GetMarker_Config_Tag(string val_marker) {
  
  unsigned short iMarker_Config;
  
  for (iMarker_Config = 0; iMarker_Config < nMarker_Config; iMarker_Config++)
    if (Marker_Config_Tag[iMarker_Config] == val_marker)
      return iMarker_Config;
  
  cout <<"The configuration file doesn't have any definition for marker "<< val_marker <<"!!" << endl;
  exit(1);
}

unsigned short CConfig::GetMarker_Config_Boundary(string val_marker) {
  unsigned short iMarker_Config;
  for (iMarker_Config = 0; iMarker_Config < nMarker_Config; iMarker_Config++)
    if (Marker_Config_Tag[iMarker_Config] == val_marker) break;
  return Marker_Config_Boundary[iMarker_Config];
}

unsigned short CConfig::GetMarker_Config_Monitoring(string val_marker) {
  unsigned short iMarker_Config;
  for (iMarker_Config = 0; iMarker_Config < nMarker_Config; iMarker_Config++)
    if (Marker_Config_Tag[iMarker_Config] == val_marker) break;
  return Marker_Config_Monitoring[iMarker_Config];
}

unsigned short CConfig::GetMarker_Config_Designing(string val_marker) {
  unsigned short iMarker_Config;
  for (iMarker_Config = 0; iMarker_Config < nMarker_Config; iMarker_Config++)
    if (Marker_Config_Tag[iMarker_Config] == val_marker) break;
  return Marker_Config_Designing[iMarker_Config];
}

unsigned short CConfig::GetMarker_Config_Plotting(string val_marker) {
  unsigned short iMarker_Config;
  for (iMarker_Config = 0; iMarker_Config < nMarker_Config; iMarker_Config++)
    if (Marker_Config_Tag[iMarker_Config] == val_marker) break;
  return Marker_Config_Plotting[iMarker_Config];
}

unsigned short CConfig::GetMarker_Config_DV(string val_marker) {
  unsigned short iMarker_Config;
  for (iMarker_Config = 0; iMarker_Config < nMarker_Config; iMarker_Config++)
    if (Marker_Config_Tag[iMarker_Config] == val_marker) break;
  return Marker_Config_DV[iMarker_Config];
}

unsigned short CConfig::GetMarker_Config_Moving(string val_marker) {
  unsigned short iMarker_Config;
  for (iMarker_Config = 0; iMarker_Config < nMarker_Config; iMarker_Config++)
    if (Marker_Config_Tag[iMarker_Config] == val_marker) break;
  return Marker_Config_Moving[iMarker_Config];
}

unsigned short CConfig::GetMarker_Config_PerBound(string val_marker) {
  unsigned short iMarker_Config;
  for (iMarker_Config = 0; iMarker_Config < nMarker_Config; iMarker_Config++)
    if (Marker_Config_Tag[iMarker_Config] == val_marker) break;
  return Marker_Config_PerBound[iMarker_Config];
}

CConfig::~CConfig(void)
{
  delete [] RK_Alpha_Step;
  delete [] MG_PreSmooth;
  delete [] MG_PostSmooth;
  delete [] U_FreeStreamND;
  
  /*--- Free memory for unspecified grid motion parameters ---*/
  
  if (Kind_GridMovement != NULL)
    delete [] Kind_GridMovement;
  
  /*--- motion origin: ---*/
  if (Motion_Origin_X != NULL)
    delete [] Motion_Origin_X;
  if (Motion_Origin_Y != NULL)
    delete [] Motion_Origin_Y;
  if (Motion_Origin_Z != NULL)
    delete [] Motion_Origin_Z;
  if (MoveMotion_Origin != NULL)
    delete [] MoveMotion_Origin;
  
  /*--- rotation: ---*/
  if (Rotation_Rate_X != NULL)
    delete [] Rotation_Rate_X;
  if (Rotation_Rate_Y != NULL)
    delete [] Rotation_Rate_Y;
  if (Rotation_Rate_Z != NULL)
    delete [] Rotation_Rate_Z;
  
  /*--- pitching: ---*/
  if (Pitching_Omega_X != NULL)
    delete [] Pitching_Omega_X;
  if (Pitching_Omega_Y != NULL)
    delete [] Pitching_Omega_Y;
  if (Pitching_Omega_Z != NULL)
    delete [] Pitching_Omega_Z;
  
  /*--- pitching amplitude: ---*/
  if (Pitching_Ampl_X != NULL)
    delete [] Pitching_Ampl_X;
  if (Pitching_Ampl_Y != NULL)
    delete [] Pitching_Ampl_Y;
  if (Pitching_Ampl_Z != NULL)
    delete [] Pitching_Ampl_Z;
  
  /*--- pitching phase: ---*/
  if (Pitching_Phase_X != NULL)
    delete [] Pitching_Phase_X;
  if (Pitching_Phase_Y != NULL)
    delete [] Pitching_Phase_Y;
  if (Pitching_Phase_Z != NULL)
    delete [] Pitching_Phase_Z;
  
  /*--- plunging: ---*/
  if (Plunging_Omega_X != NULL)
    delete [] Plunging_Omega_X;
  if (Plunging_Omega_Y != NULL)
    delete [] Plunging_Omega_Y;
  if (Plunging_Omega_Z != NULL)
    delete [] Plunging_Omega_Z;
  
  /*--- plunging amplitude: ---*/
  if (Plunging_Ampl_X != NULL)
    delete [] Plunging_Ampl_X;
  if (Plunging_Ampl_Y != NULL)
    delete [] Plunging_Ampl_Y;
  if (Plunging_Ampl_Z != NULL)
    delete [] Plunging_Ampl_Z;
  
  if (RefOriginMoment != NULL)
    delete [] RefOriginMoment;
  if (RefOriginMoment_X != NULL)
    delete [] RefOriginMoment_X;
  if (RefOriginMoment_Y != NULL)
    delete [] RefOriginMoment_Y;
  if (RefOriginMoment_Z != NULL)
    delete [] RefOriginMoment_Z;
}

void CConfig::SetFileNameDomain(unsigned short val_domain) { }

string CConfig::GetUnsteady_FileName(string val_filename, int val_iter) {
  
  string UnstExt, UnstFilename = val_filename;
  char buffer[50];
  
  /*--- Check that a positive value iteration is requested (for now). ---*/
  if (val_iter < 0) {
    cout << "Requesting a negative iteration number for the restart file!!" << endl;
    exit(1);
  }
  
  /*--- Append iteration number for unsteady cases ---*/
  if ((Wrt_Unsteady) || (Unsteady_Simulation == TIME_SPECTRAL)) {
    unsigned short lastindex = UnstFilename.find_last_of(".");
    UnstFilename = UnstFilename.substr(0, lastindex);
    if ((val_iter >= 0)    && (val_iter < 10))    sprintf (buffer, "_0000%d.dat", val_iter);
    if ((val_iter >= 10)   && (val_iter < 100))   sprintf (buffer, "_000%d.dat",  val_iter);
    if ((val_iter >= 100)  && (val_iter < 1000))  sprintf (buffer, "_00%d.dat",   val_iter);
    if ((val_iter >= 1000) && (val_iter < 10000)) sprintf (buffer, "_0%d.dat",    val_iter);
    if (val_iter >= 10000) sprintf (buffer, "_%d.dat", val_iter);
    string UnstExt = string(buffer);
    UnstFilename.append(UnstExt);
  }
  
  return UnstFilename;
}

string CConfig::GetObjFunc_Extension(string val_filename) {
  
  string AdjExt, Filename = val_filename;
  
  if (Adjoint) {
    
    /*--- Remove filename extension (.dat) ---*/
    unsigned short lastindex = Filename.find_last_of(".");
    Filename = Filename.substr(0, lastindex);
    
    switch (Kind_ObjFunc) {
      case DRAG_COEFFICIENT:      AdjExt = "_cd";   break;
      case LIFT_COEFFICIENT:      AdjExt = "_cl";   break;
      case SIDEFORCE_COEFFICIENT: AdjExt = "_csf";  break;
      case PRESSURE_COEFFICIENT:  AdjExt = "_cp";   break;
      case MOMENT_X_COEFFICIENT:  AdjExt = "_cmx";  break;
      case MOMENT_Y_COEFFICIENT:  AdjExt = "_cmy";  break;
      case MOMENT_Z_COEFFICIENT:  AdjExt = "_cmz";  break;
      case EFFICIENCY:            AdjExt = "_eff";  break;
      case EQUIVALENT_AREA:       AdjExt = "_ea";   break;
      case NEARFIELD_PRESSURE:    AdjExt = "_nfp";  break;
      case FORCE_X_COEFFICIENT:   AdjExt = "_cfx";  break;
      case FORCE_Y_COEFFICIENT:   AdjExt = "_cfy";  break;
      case FORCE_Z_COEFFICIENT:   AdjExt = "_cfz";  break;
      case THRUST_COEFFICIENT:    AdjExt = "_ct";   break;
      case TORQUE_COEFFICIENT:    AdjExt = "_cq";   break;
      case HEAT_LOAD:             AdjExt = "_Q";    break;
      case MAX_HEAT_FLUX:         AdjExt = "_qmax"; break;
      case FIGURE_OF_MERIT:       AdjExt = "_merit";break;
      case FREE_SURFACE:          AdjExt = "_fs";   break;
    }
    Filename.append(AdjExt);
    
    /*--- Lastly, add the .dat extension ---*/
    Filename.append(".dat");
    
  }
  
  return Filename;
}

unsigned short CConfig::GetContainerPosition(unsigned short val_eqsystem) {
  
  switch (val_eqsystem) {
    case RUNTIME_FLOW_SYS:      return FLOW_SOL;
    case RUNTIME_TURB_SYS:      return TURB_SOL;
    case RUNTIME_MULTIGRID_SYS: return 0;
  }
  return 0;
}

void CConfig::SetKind_ConvNumScheme(unsigned short val_kind_convnumscheme,
                                    unsigned short val_kind_centered, unsigned short val_kind_upwind,
                                    unsigned short val_kind_slopelimit) {
  
  Kind_ConvNumScheme = val_kind_convnumscheme;
  Kind_Centered = val_kind_centered;
  Kind_Upwind = val_kind_upwind;
  Kind_SlopeLimit = val_kind_slopelimit;
  
}

void CConfig::UpdateCFL(unsigned long val_iter) {
  double coeff;
  bool change;
  unsigned short iCFL;
  
  if ((val_iter % int(CFLRamp[1]) == 0 ) && (val_iter != 0)) {
    change = false;
    for (iCFL = 0; iCFL <= nMultiLevel; iCFL++) {
      coeff = pow(MG_CFLRedCoeff, double(iCFL));
      if (Adjoint) coeff = coeff * Adj_CFLRedCoeff;
      
      if (CFL[iCFL]*CFLRamp[0] < CFLRamp[2]*coeff) {
        CFL[iCFL] = CFL[iCFL]*CFLRamp[0];
        change = true;
      }
    }
    
    if (change) {
      cout <<"\n New value of the CFL number: ";
      for (iCFL = 0; iCFL < nMultiLevel; iCFL++)
        cout << CFL[iCFL] <<", ";
      cout << CFL[nMultiLevel] <<".\n"<< endl;
    }
  }
}

void CConfig::SetGlobalParam(unsigned short val_solver,
                             unsigned short val_system,
                             unsigned long val_extiter) {
  
  /*--- Set the simulation global time ---*/
  Current_UnstTime = static_cast<double>(val_extiter)*Delta_UnstTime;
  Current_UnstTimeND = static_cast<double>(val_extiter)*Delta_UnstTimeND;
  
  /*--- Set the solver methods ---*/
  switch (val_solver) {
    case EULER:
      if (val_system == RUNTIME_FLOW_SYS) {
        SetKind_ConvNumScheme(GetKind_ConvNumScheme_Flow(), GetKind_Centered_Flow(),
                              GetKind_Upwind_Flow(), GetKind_SlopeLimit_Flow());
        SetKind_SourNumScheme(GetKind_SourNumScheme_Flow());
        SetKind_ViscNumScheme(NONE);
        SetKind_TimeIntScheme(GetKind_TimeIntScheme_Flow());
      }
      break;
    case NAVIER_STOKES:
      if (val_system == RUNTIME_FLOW_SYS) {
        SetKind_ConvNumScheme(GetKind_ConvNumScheme_Flow(), GetKind_Centered_Flow(),
                              GetKind_Upwind_Flow(), GetKind_SlopeLimit_Flow());
        SetKind_SourNumScheme(GetKind_SourNumScheme_Flow());
        SetKind_ViscNumScheme(GetKind_ViscNumScheme_Flow());
        SetKind_TimeIntScheme(GetKind_TimeIntScheme_Flow());
      }
      break;
    case RANS:
      if (val_system == RUNTIME_FLOW_SYS) {
        SetKind_ConvNumScheme(GetKind_ConvNumScheme_Flow(), GetKind_Centered_Flow(),
                              GetKind_Upwind_Flow(), GetKind_SlopeLimit_Flow());
        SetKind_SourNumScheme(GetKind_SourNumScheme_Flow());
        SetKind_ViscNumScheme(GetKind_ViscNumScheme_Flow());
        SetKind_TimeIntScheme(GetKind_TimeIntScheme_Flow());
      }
      if (val_system == RUNTIME_TURB_SYS) {
        SetKind_ConvNumScheme(GetKind_ConvNumScheme_Turb(), GetKind_Centered_Turb(),
                              GetKind_Upwind_Turb(), GetKind_SlopeLimit_Turb());
        SetKind_ViscNumScheme(GetKind_ViscNumScheme_Turb());
        SetKind_SourNumScheme(GetKind_SourNumScheme_Turb());
        SetKind_TimeIntScheme(GetKind_TimeIntScheme_Turb());
      }
      break;
  }
}

double CConfig::GetInlet_Ttotal(string val_marker) {
  unsigned short iMarker_Inlet;
  for (iMarker_Inlet = 0; iMarker_Inlet < nMarker_Inlet; iMarker_Inlet++)
    if (Marker_Inlet[iMarker_Inlet] == val_marker) break;
  return Inlet_Ttotal[iMarker_Inlet];
}

double CConfig::GetInlet_Ptotal(string val_marker) {
  unsigned short iMarker_Inlet;
  for (iMarker_Inlet = 0; iMarker_Inlet < nMarker_Inlet; iMarker_Inlet++)
    if (Marker_Inlet[iMarker_Inlet] == val_marker) break;
  return Inlet_Ptotal[iMarker_Inlet];
}

double* CConfig::GetInlet_FlowDir(string val_marker) {
  unsigned short iMarker_Inlet;
  for (iMarker_Inlet = 0; iMarker_Inlet < nMarker_Inlet; iMarker_Inlet++)
    if (Marker_Inlet[iMarker_Inlet] == val_marker) break;
  return Inlet_FlowDir[iMarker_Inlet];
}

double CConfig::GetOutlet_Pressure(string val_marker) {
  unsigned short iMarker_Outlet;
  for (iMarker_Outlet = 0; iMarker_Outlet < nMarker_Outlet; iMarker_Outlet++)
    if (Marker_Outlet[iMarker_Outlet] == val_marker) break;
  return Outlet_Pressure[iMarker_Outlet];
}

double CConfig::GetIsothermal_Temperature(string val_marker) {
  unsigned short iMarker_Isothermal;
  for (iMarker_Isothermal = 0; iMarker_Isothermal < nMarker_Isothermal; iMarker_Isothermal++)
    if (Marker_Isothermal[iMarker_Isothermal] == val_marker) break;
  return Isothermal_Temperature[iMarker_Isothermal];
}

double CConfig::GetWall_HeatFlux(string val_marker) {
  unsigned short iMarker_HeatFlux;
  for (iMarker_HeatFlux = 0; iMarker_HeatFlux < nMarker_HeatFlux; iMarker_HeatFlux++)
    if (Marker_HeatFlux[iMarker_HeatFlux] == val_marker) break;
  return Heat_Flux[iMarker_HeatFlux];
}

void CConfig::SetNondimensionalization(unsigned short val_nDim) {
  
  double Mach2Vel_FreeStream, ModVel_FreeStream, Energy_FreeStream = 0.0, ModVel_FreeStreamND;
  double Velocity_Reynolds;
  unsigned short iDim;
  
  Velocity_FreeStreamND = new double[val_nDim];
  
  /*--- Local variables and memory allocation ---*/
  double Alpha = AoA*PI_NUMBER/180.0;
  double Beta  = AoS*PI_NUMBER/180.0;
  double Gamma_Minus_One = Gamma - 1.0;
  bool Unsteady = (Unsteady_Simulation != NO);
  bool turbulent = (Kind_Solver == RANS);
  
  Mach2Vel_FreeStream = sqrt(Gamma*Gas_Constant*Temperature_FreeStream);
  
  /*--- Compute the Free Stream velocity, using the Mach number ---*/
  if (val_nDim == 2) {
    Velocity_FreeStream[0] = cos(Alpha)*Mach*Mach2Vel_FreeStream;
    Velocity_FreeStream[1] = sin(Alpha)*Mach*Mach2Vel_FreeStream;
  }
  if (val_nDim == 3) {
    Velocity_FreeStream[0] = cos(Alpha)*cos(Beta)*Mach*Mach2Vel_FreeStream;
    Velocity_FreeStream[1] = sin(Beta)*Mach*Mach2Vel_FreeStream;
    Velocity_FreeStream[2] = sin(Alpha)*cos(Beta)*Mach*Mach2Vel_FreeStream;
  }
  
  /*--- Compute the modulus of the free stream velocity ---*/
  ModVel_FreeStream = 0;
  for (iDim = 0; iDim < val_nDim; iDim++)
    ModVel_FreeStream += Velocity_FreeStream[iDim]*Velocity_FreeStream[iDim];
  ModVel_FreeStream = sqrt(ModVel_FreeStream);
  
  if (Viscous) {
    
    /*--- First, check if there is mesh motion. If yes, use the Mach
     number relative to the body to initialize the flow. ---*/
    if (Grid_Movement)
      Velocity_Reynolds = Mach_Motion*Mach2Vel_FreeStream;
    else
      Velocity_Reynolds = ModVel_FreeStream;
    
    /*--- For viscous flows, pressure will be computed from a density
     that is found from the Reynolds number. The viscosity is computed
     from the dimensional version of Sutherland's law ---*/
    Viscosity_FreeStream = 1.853E-5*(pow(Temperature_FreeStream/300.0,3.0/2.0) * (300.0+110.3)/(Temperature_FreeStream+110.3));
    Density_FreeStream   = Reynolds*Viscosity_FreeStream/(Velocity_Reynolds*Length_Reynolds);
    Pressure_FreeStream  = Density_FreeStream*Gas_Constant*Temperature_FreeStream;
    
  } else {
    /*--- For inviscid flow, density is calculated from the specified
     total temperature and pressure using the gas law. ---*/
    Density_FreeStream  = Pressure_FreeStream/(Gas_Constant*Temperature_FreeStream);
  }
  /*-- Compute the freestream energy. ---*/
  Energy_FreeStream = Pressure_FreeStream/(Density_FreeStream*Gamma_Minus_One)+0.5*ModVel_FreeStream*ModVel_FreeStream;
  
  /*--- Additional reference values defined by Pref, Tref, RHOref. By definition,
   Lref is one because we have converted the grid to meters.---*/
  Length_Ref         = 1.0;
  Velocity_Ref      = sqrt(Pressure_Ref/Density_Ref);
  Time_Ref          = Length_Ref/Velocity_Ref;
  Omega_Ref         = Velocity_Ref/Length_Ref;
  Force_Ref         = Velocity_Ref*Velocity_Ref/Length_Ref;
  Gas_Constant_Ref  = Velocity_Ref*Velocity_Ref/Temperature_Ref;
  Viscosity_Ref     = Density_Ref*Velocity_Ref*Length_Ref;
  Froude            = ModVel_FreeStream/sqrt(STANDART_GRAVITY*Length_Ref);
  
  /*--- Divide by reference values, to compute the non-dimensional free-stream values ---*/
  Pressure_FreeStreamND = Pressure_FreeStream/Pressure_Ref;
  Density_FreeStreamND  = Density_FreeStream/Density_Ref;
  
  for (iDim = 0; iDim < val_nDim; iDim++)
    Velocity_FreeStreamND[iDim] = Velocity_FreeStream[iDim]/Velocity_Ref;
  Temperature_FreeStreamND = Temperature_FreeStream/Temperature_Ref;
  
  Gas_ConstantND = Gas_Constant/Gas_Constant_Ref;
  
  ModVel_FreeStreamND = 0;
  for (iDim = 0; iDim < val_nDim; iDim++)
    ModVel_FreeStreamND += Velocity_FreeStreamND[iDim]*Velocity_FreeStreamND[iDim];
  ModVel_FreeStreamND    = sqrt(ModVel_FreeStreamND);
  Energy_FreeStreamND    = Pressure_FreeStreamND/(Density_FreeStreamND*Gamma_Minus_One)+0.5*ModVel_FreeStreamND*ModVel_FreeStreamND;
  Viscosity_FreeStreamND = Viscosity_FreeStream / Viscosity_Ref;
  Total_UnstTimeND = Total_UnstTime / Time_Ref;
  Delta_UnstTimeND = Delta_UnstTime / Time_Ref;
  
  double kine_Inf  = 3.0/2.0*(ModVel_FreeStreamND*ModVel_FreeStreamND*TurbulenceIntensity_FreeStream*TurbulenceIntensity_FreeStream);
  double omega_Inf = Density_FreeStreamND*kine_Inf/(Viscosity_FreeStreamND*Turb2LamViscRatio_FreeStream);
  
  cout << endl <<"---------------- Flow & Non-dimensionalization information ---------------" << endl;
  
  cout.precision(6);
  
  if (Viscous) {
    cout << "Viscous flow: Computing pressure using the ideal gas law" << endl;
    cout << "based on the freestream temperature and a density computed" << endl;
    cout << "from the Reynolds number." << endl;
  } else {
    cout << "Inviscid flow: Computing density based on freestream" << endl;
    cout << "temperature and pressure using the ideal gas law." << endl;
  }
  
  cout <<"--Input conditions:"<< endl;
  cout << "Grid conversion factor to meters: " << Conversion_Factor << endl;
  
  cout << "Ratio of specific heats: " << Gamma           << endl;
  cout << "Specific gas constant (J/(kg.K)): "   << Gas_Constant  << endl;
  
  cout << "Freestream pressure (N/m^2): "          << Pressure_FreeStream    << endl;
  cout << "Freestream temperature (K): "       << Temperature_FreeStream << endl;
  cout << "Freestream density (kg/m^3): "					 << Density_FreeStream << endl;
  if (val_nDim == 2) {
    cout << "Freestream velocity (m/s): (" << Velocity_FreeStream[0] << ",";
    cout << Velocity_FreeStream[1] << ")" << endl;
  } else if (val_nDim == 3) {
    cout << "Freestream velocity (m/s): (" << Velocity_FreeStream[0] << ",";
    cout << Velocity_FreeStream[1] << "," << Velocity_FreeStream[2] << ")" << endl;
  }
  
  cout << "Freestream velocity magnitude (m/s):"	<< ModVel_FreeStream << endl;
  
  cout << "Freestream energy (kg.m/s^2): "					 << Energy_FreeStream << endl;
  
  if (Viscous)
    cout << "Freestream viscosity (N.s/m^2): "				 << Viscosity_FreeStream << endl;
  
  if (Unsteady) {
    cout << "Total time (s): " << Total_UnstTime << ". Time step (s): " << Delta_UnstTime << endl;
  }
  
  /*--- Print out reference values. ---*/
  cout <<"--Reference values:"<< endl;
  cout << "Reference pressure (N/m^2): "      << Pressure_Ref    << endl;
  
  cout << "Reference temperature (K): "   << Temperature_Ref << endl;
  cout << "Reference energy (kg.m/s^2): "       << Energy_FreeStream/Energy_FreeStreamND     << endl;
  cout << "Reference density (kg/m^3): "       << Density_Ref     << endl;
  cout << "Reference velocity (m/s): "       << Velocity_Ref     << endl;
  
  if (Viscous)
    cout << "Reference viscosity (N.s/m^2): "       << Viscosity_Ref     << endl;
  
  if (Unsteady)
    cout << "Reference time (s): "        << Time_Ref      << endl;
  
  /*--- Print out resulting non-dim values here. ---*/
  cout << "--Resulting non-dimensional state:" << endl;
  cout << "Mach number (non-dimensional): " << Mach << endl;
  if (Viscous) {
    cout << "Reynolds number (non-dimensional): " << Reynolds << endl;
    cout << "Reynolds length (m): "       << Length_Reynolds     << endl;
  }
  if (GravityForce) {
    cout << "Froude number (non-dimensional): " << Froude << endl;
    cout << "Lenght of the baseline wave (non-dimensional): " << 2.0*PI_NUMBER*Froude*Froude << endl;
  }
  
  cout << "Specific gas constant (non-dimensional): "   << Gas_Constant << endl;
  cout << "Freestream temperature (non-dimensional): "  << Temperature_FreeStreamND << endl;
  
  cout << "Freestream pressure (non-dimensional): "     << Pressure_FreeStreamND    << endl;
  cout << "Freestream density (non-dimensional): "      << Density_FreeStreamND     << endl;
  if (val_nDim == 2) {
    cout << "Freestream velocity (non-dimensional): (" << Velocity_FreeStreamND[0] << ",";
    cout << Velocity_FreeStreamND[1] << ")" << endl;
  } else if (val_nDim == 3) {
    cout << "Freestream velocity (non-dimensional): (" << Velocity_FreeStreamND[0] << ",";
    cout << Velocity_FreeStreamND[1] << "," << Velocity_FreeStreamND[2] << ")" << endl;
  }
  cout << "Freestream velocity magnitude (non-dimensional): "	 << ModVel_FreeStreamND << endl;
  
  if (turbulent){
    cout << "Free-stream turb. kinetic energy (non-dimensional): " << kine_Inf << endl;
    cout << "Free-stream specific dissipation (non-dimensional): " << omega_Inf << endl;
  }
  
  cout << "Freestream energy (non-dimensional): "					 << Energy_FreeStreamND << endl;
  
  if (Viscous)
    cout << "Freestream viscosity (non-dimensional): " << Viscosity_FreeStreamND << endl;
  
  if (Unsteady) {
    cout << "Total time (non-dimensional): "				 << Total_UnstTimeND << endl;
    cout << "Time step (non-dimensional): "				 << Delta_UnstTimeND << endl;
  }
  if (Grid_Movement) cout << "Force coefficients computed using MACH_MOTION." << endl;
  else cout << "Force coefficients computed using freestream values." << endl;
  
  cout << "Note: Negative pressure, temperature or density is not allowed!" << endl;
  
  
}

void CConfig::SetSpline(vector<double> &x, vector<double> &y, unsigned long n, double yp1, double ypn, vector<double> &y2) {
  unsigned long i, k;
  double p, qn, sig, un, *u;
  
  u = new double [n];
  
  if (yp1 > 0.99e30)			// The lower boundary condition is set either to be "nat
    y2[0]=u[0]=0.0;			  // -ural"
  else {									// or else to have a specified first derivative.
    y2[0] = -0.5;
    u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
  }
  
  for (i=2; i<=n-1; i++) {									//  This is the decomposition loop of the tridiagonal al-
    sig=(x[i-1]-x[i-2])/(x[i]-x[i-2]);		//	gorithm. y2 and u are used for tem-
    p=sig*y2[i-2]+2.0;										//	porary storage of the decomposed
    y2[i-1]=(sig-1.0)/p;										//	factors.
    u[i-1]=(y[i]-y[i-1])/(x[i]-x[i-1]) - (y[i-1]-y[i-2])/(x[i-1]-x[i-2]);
    u[i-1]=(6.0*u[i-1]/(x[i]-x[i-2])-sig*u[i-2])/p;
  }
  
  if (ypn > 0.99e30)						// The upper boundary condition is set either to be
    qn=un=0.0;									// "natural"
  else {												// or else to have a specified first derivative.
    qn=0.5;
    un=(3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
  }
  y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
  for (k=n-1; k>=1; k--)					// This is the backsubstitution loop of the tridiagonal
    y2[k-1]=y2[k-1]*y2[k]+u[k-1];	  // algorithm.
  
  delete[] u;
  
}

double CConfig::GetSpline(vector<double>&xa, vector<double>&ya, vector<double>&y2a, unsigned long n, double x) {
  unsigned long klo, khi, k;
  double h, b, a, y;
  
  klo=1;										// We will find the right place in the table by means of
  khi=n;										// bisection. This is optimal if sequential calls to this
  while (khi-klo > 1) {			// routine are at random values of x. If sequential calls
    k=(khi+klo) >> 1;				// are in order, and closely spaced, one would do better
    if (xa[k-1] > x) khi=k;		// to store previous values of klo and khi and test if
    else klo=k;							// they remain appropriate on the next call.
  }								// klo and khi now bracket the input value of x
  h=xa[khi-1]-xa[klo-1];
  if (h == 0.0) cout << "Bad xa input to routine splint" << endl;	// The xa’s must be dis-
  a=(xa[khi-1]-x)/h;																					      // tinct.
  b=(x-xa[klo-1])/h;				// Cubic spline polynomial is now evaluated.
  y=a*ya[klo-1]+b*ya[khi-1]+((a*a*a-a)*y2a[klo-1]+(b*b*b-b)*y2a[khi-1])*(h*h)/6.0;
  
  return y;
}
