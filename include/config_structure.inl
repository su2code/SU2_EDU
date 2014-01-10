/*!
 * \file config_structure.inl
 * \brief In-Line subroutines of the <i>config_structure.hpp</i> file.
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

#pragma once

inline unsigned short CConfig::GetnZone(void) { return nZone; }

inline unsigned short CConfig::GetiZone(void) { return iZone; }

inline unsigned short CConfig::GetKind_SU2(void) { return Kind_SU2; }

inline bool CConfig::GetAdjoint(void) { return Adjoint; }

inline bool CConfig::GetViscous(void) { return Viscous; }

inline unsigned long CConfig::GetnExtIter(void) { return nExtIter; }

inline void CConfig::SetExtIter(unsigned long val_iter) { ExtIter = val_iter; }

inline void CConfig::SetIntIter(unsigned long val_iter) { IntIter = val_iter; }

inline unsigned long CConfig::GetExtIter(void) { return ExtIter; }

inline unsigned long CConfig::GetIntIter(void) { return IntIter; }

inline unsigned short CConfig::GetMaxChildren(void) { return MaxChildren; }

inline double *CConfig::GetHold_GridFixed_Coord(void) { return Hold_GridFixed_Coord; }

inline double CConfig::GetMaxDimension(void) { return MaxDimension; }

inline double CConfig::GetRatioDensity(void) { return RatioDensity; }

inline double CConfig::GetRatioViscosity(void) { return RatioViscosity; }

inline unsigned short CConfig::GetAnalytical_Surface(void) { return Analytical_Surface; }

inline double CConfig::GetDualVol_Power(void) { return DualVol_Power; }

inline bool CConfig::GetVisualize_Partition(void) { return Visualize_Partition; }

inline bool CConfig::GetExtraOutput(void) { return ExtraOutput; }

inline bool CConfig::GetVisualize_Deformation(void) { return Visualize_Deformation; }

inline double CConfig::GetRefAreaCoeff(void) { return RefAreaCoeff; }

inline double CConfig::GetRefLengthMoment(void) { return RefLengthMoment; }

inline double CConfig::GetRefElemLength(void) { return RefElemLength; }

inline double CConfig::GetRefSharpEdges(void) { return RefSharpEdges; }

inline double CConfig::GetDomainVolume(void) { return DomainVolume; }

inline void CConfig::SetRefAreaCoeff(double val_area) { RefAreaCoeff = val_area; }

inline void CConfig::SetDomainVolume(double val_volume) { DomainVolume = val_volume; }

inline void CConfig::SetnExtIter(unsigned long val_niter) { nExtIter = val_niter; }

inline double CConfig::GetMach_FreeStreamND(void) { return Mach; }

inline double CConfig::GetGamma(void) { return Gamma; }

inline double CConfig::GetSection_Limit(unsigned short val_var) { return Section_Limit[val_var]; }

inline double CConfig::GetBulk_Modulus(void) { return Bulk_Modulus; }

inline double CConfig::GetArtComp_Factor(void) { return ArtComp_Factor; }

inline double CConfig::GetGas_Constant(void) { return Gas_Constant; }

inline double CConfig::GetGas_ConstantND(void) { return Gas_ConstantND; }

inline double CConfig::GetWallTemperature(void) { return Wall_Temperature; }

inline double CConfig::GetFreeSurface_Zero(void) { return FreeSurface_Zero; }

inline double CConfig::GetFreeSurface_Depth(void) { return FreeSurface_Depth; }

inline double CConfig::GetGas_Constant_Ref(void) { return Gas_Constant_Ref; }

inline double CConfig::GetTemperature_FreeStream(void) { return Temperature_FreeStream; }

inline double CConfig::GetPrandtl_Lam(void) { return Prandtl_Lam; }

inline double CConfig::GetPrandtl_Turb(void) { return Prandtl_Turb; }

inline double CConfig::GetLength_Ref(void) { return Length_Ref; }

inline double CConfig::GetPressure_Ref(void) { return Pressure_Ref; }

inline double CConfig::GetTemperature_Ref(void) { return Temperature_Ref; }

inline double CConfig::GetDensity_Ref(void) { return Density_Ref; }

inline double CConfig::GetVelocity_Ref(void) { return Velocity_Ref; }

inline double CConfig::GetTime_Ref(void) { return Time_Ref; }

inline double CConfig::GetViscosity_Ref(void) { return Viscosity_Ref; }

inline double CConfig::GetOmega_Ref(void) { return Omega_Ref; }

inline double CConfig::GetForce_Ref(void) { return Force_Ref; }

inline double CConfig::GetPressure_FreeStreamND(void) { return Pressure_FreeStreamND; }

inline double CConfig::GetPressure_FreeStream(void) { return Pressure_FreeStream; }

inline double CConfig::GetTemperature_FreeStreamND(void) { return Temperature_FreeStreamND; }

inline double CConfig::GetDensity_FreeStreamND(void) { return Density_FreeStreamND; }

inline double* CConfig::GetVelocity_FreeStreamND(void) { return Velocity_FreeStreamND; }

inline double* CConfig::GetVelocity_FreeStream(void) { return Velocity_FreeStream; }

inline double CConfig::GetEnergy_FreeStreamND(void) { return Energy_FreeStreamND; }

inline double CConfig::GetViscosity_FreeStreamND(void) { return Viscosity_FreeStreamND; }

inline double CConfig::GetNuFactor_FreeStream(void) { return NuFactor_FreeStream; }

inline double CConfig::GetIntermittency_FreeStream(void) { return Intermittency_FreeStream; }

inline double CConfig::GetTurbulenceIntensity_FreeStream(void) { return TurbulenceIntensity_FreeStream; }

inline double CConfig::GetTurb2LamViscRatio_FreeStream(void) {return Turb2LamViscRatio_FreeStream;}

inline double CConfig::GetLength_Reynolds(void) { return Length_Reynolds; }

inline double CConfig::GetConversion_Factor(void) { return Conversion_Factor; }

inline unsigned short CConfig::GetnStartUpIter(void) { return nStartUpIter; }

inline double *CConfig::GetRefOriginMoment(unsigned short val_marker) {
    RefOriginMoment[0] = RefOriginMoment_X[val_marker];
    RefOriginMoment[1] = RefOriginMoment_Y[val_marker];
    RefOriginMoment[2] = RefOriginMoment_Z[val_marker];
    return RefOriginMoment;
}

inline double CConfig::GetRefOriginMoment_X(unsigned short val_marker) { return RefOriginMoment_X[val_marker]; }

inline double CConfig::GetRefOriginMoment_Y(unsigned short val_marker) { return RefOriginMoment_Y[val_marker]; }

inline double CConfig::GetRefOriginMoment_Z(unsigned short val_marker) { return RefOriginMoment_Z[val_marker]; }

inline void CConfig::SetRefOriginMoment_X(unsigned short val_marker, double val_origin) { RefOriginMoment_X[val_marker] = val_origin; }

inline void CConfig::SetRefOriginMoment_Y(unsigned short val_marker, double val_origin) { RefOriginMoment_Y[val_marker] = val_origin; }

inline void CConfig::SetRefOriginMoment_Z(unsigned short val_marker, double val_origin) { RefOriginMoment_Z[val_marker] = val_origin; }

inline double CConfig::GetChargeCoeff(void) { return ChargeCoeff; }

inline double CConfig::GetLimiterCoeff(void) { return LimiterCoeff; }

inline double CConfig::GetSharpEdgesCoeff(void) { return SharpEdgesCoeff; }

inline double CConfig::GetReynolds(void) { return Reynolds; }

inline double CConfig::GetFroude(void) { return Froude; }

inline double CConfig::GetAoA(void) { return AoA; }

inline unsigned short CConfig::GetnDomain(void) { return nDomain; }

inline void CConfig::SetnDomain(unsigned short val_ndomain) { nDomain = val_ndomain; }

inline double CConfig::GetAoS(void) { return AoS; }

inline unsigned short CConfig::GetMGLevels(void) { return nMultiLevel; }

inline void CConfig::SetMGLevels(unsigned short val_nMultiLevel) { nMultiLevel = val_nMultiLevel; }

inline unsigned short CConfig::GetFinestMesh(void) { return FinestMesh; }

inline void CConfig::SetFinestMesh(unsigned short val_finestmesh) { FinestMesh = val_finestmesh; }

inline void CConfig::SubtractFinestMesh(void) { FinestMesh = FinestMesh-1; }

inline unsigned short CConfig::GetDesign_Variable(unsigned short val_dv) { return Design_Variable[val_dv]; }

inline unsigned short CConfig::GetConvCriteria(void) { return ConvCriteria; }

inline unsigned short CConfig::GetMGCycle(void) { return MGCycle; }

inline unsigned short CConfig::GetGeometryMode(void) { return GeometryMode; }

inline double CConfig::GetCFL(unsigned short val_mesh) {	return CFL[val_mesh]; }

inline double CConfig::GetUnst_CFL(void) {	return Unst_CFL; }

inline double CConfig::GetParamDV(unsigned short val_dv, unsigned short val_param) {	return ParamDV[val_dv][val_param]; }

inline unsigned short CConfig::GetnDV(void) {	return nDV; }

inline unsigned short CConfig::GetnRKStep(void) { return nRKStep; }

inline double CConfig::Get_Alpha_RKStep(unsigned short val_step) { return RK_Alpha_Step[val_step]; }

inline unsigned short CConfig::GetMG_PreSmooth(unsigned short val_mesh) {	
	if (nMG_PreSmooth == 0) return 1;
	else return MG_PreSmooth[val_mesh]; 
}

inline unsigned short CConfig::GetMG_PostSmooth(unsigned short val_mesh) { 
	if (nMG_PostSmooth == 0) return 0;
	else return MG_PostSmooth[val_mesh];
}

inline unsigned short CConfig::GetMG_CorrecSmooth(unsigned short val_mesh) { 
	if (nMG_CorrecSmooth == 0) return 0;
	else return MG_CorrecSmooth[val_mesh]; 
}

inline unsigned long CConfig::GetWrt_Sol_Freq(void) { return Wrt_Sol_Freq; }

inline unsigned long CConfig::GetWrt_Con_Freq(void) { return Wrt_Con_Freq; }

inline unsigned short CConfig::GetKind_Solver(void) { return Kind_Solver; }

inline unsigned short CConfig::GetKind_Regime(void) { return Kind_Regime; }

inline double CConfig::GetminTurkelBeta() { return  Min_Beta_RoeTurkel; }

inline double CConfig::GetmaxTurkelBeta() { return  Max_Beta_RoeTurkel; }

inline unsigned short CConfig::GetKind_Gradient_Method(void) { return Kind_Gradient_Method; }

inline unsigned short CConfig::GetKind_Linear_Solver(void) { return Kind_Linear_Solver; }

inline unsigned short CConfig::GetKind_Linear_Solver_Prec(void) { return Kind_Linear_Solver_Prec; }

inline void CConfig::SetKind_Linear_Solver_Prec(unsigned short val_kind_prec) { Kind_Linear_Solver_Prec = val_kind_prec; }

inline double CConfig::GetLinear_Solver_Error(void) { return Linear_Solver_Error; }

inline unsigned long CConfig::GetLinear_Solver_Iter(void) { return Linear_Solver_Iter; }

inline double CConfig::GetLinear_Solver_Relax(void) { return Linear_Solver_Relax; }

inline unsigned long CConfig::GetGridDef_Iter(void) { return GridDef_Iter; }

inline unsigned short CConfig::GetKind_TimeIntScheme(void) { return Kind_TimeNumScheme; }

inline unsigned short CConfig::GetKind_ConvNumScheme(void) { return Kind_ConvNumScheme; }

inline unsigned short CConfig::GetKind_ViscNumScheme(void) { return Kind_ViscNumScheme; }

inline unsigned short CConfig::GetKind_SourNumScheme(void) { return Kind_SourNumScheme; }

inline unsigned short CConfig::GetKind_Centered(void) { return Kind_Centered; }

inline unsigned short CConfig::GetKind_Upwind(void) { return Kind_Upwind; }

inline unsigned short CConfig::GetKind_TimeIntScheme_Flow(void) { return Kind_TimeIntScheme_Flow; }

inline unsigned short CConfig::GetKind_ConvNumScheme_Flow(void) { return Kind_ConvNumScheme_Flow; }

inline unsigned short CConfig::GetKind_ConvNumScheme_Template(void) { return Kind_ConvNumScheme_Template; }

inline unsigned short CConfig::GetKind_ViscNumScheme_Flow(void) { return Kind_ViscNumScheme_Flow; }

inline unsigned short CConfig::GetKind_ViscNumScheme_Template(void) { return Kind_ViscNumScheme_Template; }

inline unsigned short CConfig::GetKind_SourNumScheme_Flow(void) { return Kind_SourNumScheme_Flow; }

inline unsigned short CConfig::GetKind_SourNumScheme_Template(void) { return Kind_SourNumScheme_Template; }

inline unsigned short CConfig::GetKind_Centered_Flow(void) { return Kind_Centered_Flow; }

inline unsigned short CConfig::GetKind_SlopeLimit(void) { return Kind_SlopeLimit; }

inline unsigned short CConfig::GetKind_SlopeLimit_Flow(void) { return Kind_SlopeLimit_Flow; }

inline unsigned short CConfig::GetKind_SlopeLimit_Turb(void) { return Kind_SlopeLimit_Turb; }

inline unsigned short CConfig::GetKind_Upwind_Flow(void) { return Kind_Upwind_Flow; }

inline double CConfig::GetKappa_1st_Flow(void) { return Kappa_1st_Flow; }

inline double CConfig::GetKappa_2nd_Flow(void) { return Kappa_2nd_Flow; }

inline double CConfig::GetKappa_4th_Flow(void) { return Kappa_4th_Flow; }

inline unsigned short CConfig::GetKind_TimeIntScheme_Turb(void) { return Kind_TimeIntScheme_Turb; }

inline unsigned short CConfig::GetKind_ConvNumScheme_Turb(void) { return Kind_ConvNumScheme_Turb; }

inline unsigned short CConfig::GetKind_ViscNumScheme_Turb(void) { return Kind_ViscNumScheme_Turb; }

inline unsigned short CConfig::GetKind_SourNumScheme_Turb(void) { return Kind_SourNumScheme_Turb; }

inline unsigned short CConfig::GetKind_Centered_Turb(void) { return Kind_Centered_Turb; }

inline unsigned short CConfig::GetKind_Upwind_Turb(void) {	return Kind_Upwind_Turb; }

inline unsigned short CConfig::GetKind_Inlet(void) { return Kind_Inlet; }

inline void CConfig::SetKind_TimeIntScheme(unsigned short val_kind_timeintscheme) { Kind_TimeNumScheme = val_kind_timeintscheme; }

inline void CConfig::SetKind_ViscNumScheme(unsigned short val_kind_viscnumscheme) { Kind_ViscNumScheme = val_kind_viscnumscheme; }

inline void CConfig::SetKind_SourNumScheme(unsigned short val_kind_sournumscheme) { Kind_SourNumScheme = val_kind_sournumscheme; }

inline unsigned short CConfig::GetUnsteady_Simulation(void) { return Unsteady_Simulation; }

inline bool CConfig::GetRestart(void) {	return Restart; }

inline bool CConfig::GetRestart_Flow(void) { return Restart_Flow; }

inline bool CConfig::GetFullMG(void) { return FullMG; }

inline bool CConfig::GetEquivArea(void) { return EquivArea; }

inline void CConfig::SetnMarker_All(unsigned short val_nmarker) { nMarker_All = val_nmarker; }

inline string CConfig::GetMarker_All_Tag(unsigned short val_marker) { return Marker_All_Tag[val_marker]; }

inline string CConfig::GetMarker_Monitoring(unsigned short val_marker) { return Marker_Monitoring[val_marker]; }

inline unsigned short CConfig::GetTag_Marker_All(string val_tag) {
	for (unsigned short iMarker = 0; iMarker < nMarker_All; iMarker++) {
		if (val_tag == Marker_All_Tag[iMarker])
		return iMarker; 
	}
	cout <<"Ups, I don't find the boundary: "<< val_tag << endl; return 0;
}

inline unsigned short CConfig::GetMarker_All_Boundary(unsigned short val_marker) { return Marker_All_Boundary[val_marker]; }

inline void CConfig::SetMarker_All_Boundary(unsigned short val_marker, unsigned short val_boundary) { Marker_All_Boundary[val_marker] = val_boundary; }

inline void CConfig::SetMarker_All_Tag(unsigned short val_marker, string val_index) { Marker_All_Tag[val_marker] = val_index; }

inline void CConfig::SetMarker_All_Monitoring(unsigned short val_marker, unsigned short val_monitoring) { Marker_All_Monitoring[val_marker] = val_monitoring; }

inline void CConfig::SetMarker_All_Designing(unsigned short val_marker, unsigned short val_designing) { Marker_All_Designing[val_marker] = val_designing; }

inline void CConfig::SetMarker_All_Plotting(unsigned short val_marker, unsigned short val_plotting) { Marker_All_Plotting[val_marker] = val_plotting; }

inline void CConfig::SetMarker_All_DV(unsigned short val_marker, unsigned short val_DV) { Marker_All_DV[val_marker] = val_DV; }

inline void CConfig::SetMarker_All_Moving(unsigned short val_marker, unsigned short val_moving) { Marker_All_Moving[val_marker] = val_moving; }

inline void CConfig::SetMarker_All_PerBound(unsigned short val_marker, short val_perbound) { Marker_All_PerBound[val_marker] = val_perbound; }

inline short CConfig::GetMarker_All_PerBound(unsigned short val_marker) { return Marker_All_PerBound[val_marker]; }

inline unsigned short CConfig::GetMarker_All_Monitoring(unsigned short val_marker) { return Marker_All_Monitoring[val_marker]; }

inline unsigned short CConfig::GetMarker_All_Designing(unsigned short val_marker) { return Marker_All_Designing[val_marker]; }

inline short CConfig::GetMarker_All_SendRecv(unsigned short val_marker) { return Marker_All_SendRecv[val_marker]; }

inline void CConfig::SetMarker_All_SendRecv(unsigned short val_marker, short val_index) { Marker_All_SendRecv[val_marker] = val_index; }

inline unsigned short CConfig::GetMarker_All_Plotting(unsigned short val_marker) { return Marker_All_Plotting[val_marker]; }

inline unsigned short CConfig::GetMarker_All_DV(unsigned short val_marker) { return Marker_All_DV[val_marker]; }

inline unsigned short CConfig::GetMarker_All_Moving(unsigned short val_marker) { return Marker_All_Moving[val_marker]; }

inline unsigned short CConfig::GetnMarker_All(void) { return nMarker_All; }

inline unsigned short CConfig::GetnMarker_NacelleInflow(void) {	return nMarker_NacelleInflow; }

inline unsigned short CConfig::GetnMarker_NacelleExhaust(void) { return nMarker_NacelleExhaust; }

inline unsigned short CConfig::GetnMarker_InterfaceBound(void) { return nMarker_InterfaceBound; }

inline unsigned short CConfig::GetnMarker_Monitoring(void) { return nMarker_Monitoring; }

inline unsigned short CConfig::GetnMarker_Moving(void) { return nMarker_Moving; }

inline unsigned short CConfig::GetnMarker_NearFieldBound(void) { return nMarker_NearFieldBound; }

inline string CConfig::GetMesh_FileName(void) { return Mesh_FileName; }

inline string CConfig::GetMesh_Out_FileName(void) { return Mesh_Out_FileName; }

inline unsigned short CConfig::GetMesh_FileFormat(void) { return Mesh_FileFormat; }

inline unsigned short CConfig::GetOutput_FileFormat(void) { return Output_FileFormat; }

inline string CConfig::GetConv_FileName(void) { return Conv_FileName; }

inline string CConfig::GetSolution_FlowFileName(void) { return Solution_FlowFileName; }

inline string CConfig::GetSolution_LinFileName(void) { return Solution_LinFileName; }

inline string CConfig::GetSolution_AdjFileName(void) { return Solution_AdjFileName; }

inline string CConfig::GetFlow_FileName(void) { return Flow_FileName; }

inline string CConfig::GetStructure_FileName(void) { return Structure_FileName; }

inline string CConfig::GetSurfStructure_FileName(void) { return SurfStructure_FileName; }

inline string CConfig::GetSurfWave_FileName(void) { return SurfWave_FileName; }

inline string CConfig::GetSurfHeat_FileName(void) { return SurfHeat_FileName; }

inline string CConfig::GetWave_FileName(void) { return Wave_FileName; }

inline string CConfig::GetHeat_FileName(void) { return Heat_FileName; }

inline string CConfig::GetAdjWave_FileName(void) { return AdjWave_FileName; }

inline string CConfig::GetRestart_FlowFileName(void) { return Restart_FlowFileName; }

inline string CConfig::GetRestart_WaveFileName(void) { return Restart_WaveFileName; }

inline string CConfig::GetRestart_HeatFileName(void) { return Restart_HeatFileName; }

inline string CConfig::GetRestart_LinFileName(void) { return Restart_LinFileName; }

inline string CConfig::GetRestart_AdjFileName(void) { return Restart_AdjFileName; }

inline string CConfig::GetAdj_FileName(void) { return Adj_FileName; }

inline string CConfig::GetLin_FileName(void) { return Lin_FileName; }

inline string CConfig::GetObjFunc_Grad_FileName(void) { return ObjFunc_Grad_FileName; }

inline string CConfig::GetObjFunc_Value_FileName(void) { return ObjFunc_Value_FileName; }

inline string CConfig::GetSurfFlowCoeff_FileName(void) { return SurfFlowCoeff_FileName; }

inline string CConfig::GetSurfAdjCoeff_FileName(void) { return SurfAdjCoeff_FileName; }

inline string CConfig::GetSurfLinCoeff_FileName(void) { return SurfLinCoeff_FileName; }

inline unsigned short CConfig::GetCauchy_Func_Flow(void) { return Cauchy_Func_Flow; }

inline unsigned short CConfig::GetCauchy_Func_AdjFlow(void) { return Cauchy_Func_AdjFlow; }

inline unsigned short CConfig::GetCauchy_Func_LinFlow(void) { return Cauchy_Func_LinFlow; }

inline unsigned short CConfig::GetCauchy_Elems(void) { return Cauchy_Elems; }

inline unsigned long CConfig::GetStartConv_Iter(void) { return StartConv_Iter; }

inline double CConfig::GetCauchy_Eps(void) { return Cauchy_Eps; }

inline double CConfig::GetCauchy_Eps_OneShot(void) { return Cauchy_Eps_OneShot; }

inline double CConfig::GetCauchy_Eps_FullMG(void) { return Cauchy_Eps_FullMG; }

inline bool CConfig::GetDivide_Element(void) { return Divide_Element; }

inline double CConfig::GetDV_Value(unsigned short val_dv) { return DV_Value[val_dv]; }

inline double CConfig::GetOrderMagResidual(void) { return OrderMagResidual; }

inline double CConfig::GetMinLogResidual(void) { return MinLogResidual; }

inline double CConfig::GetDamp_Res_Restric(void) { return Damp_Res_Restric; }

inline double CConfig::GetDamp_Correc_Prolong(void) { return Damp_Correc_Prolong; }

inline double CConfig::GetTurb_CFLRedCoeff(void) { return Turb_CFLRedCoeff; }

inline bool CConfig::GetRotating_Frame(void) { return Rotating_Frame; }

inline bool CConfig::GetAxisymmetric(void) { return Axisymmetric; }

inline bool CConfig::GetAdaptBoundary(void) { return AdaptBoundary; }

inline bool CConfig::GetPoissonSolver(void) { return PoissonSolver; }

inline bool CConfig::Low_Mach_Preconditioning(void) { return Low_Mach_Precon; }

inline bool CConfig::GetGravityForce(void) { return GravityForce; }

inline bool CConfig::GetSmoothNumGrid(void) { return SmoothNumGrid; }

inline void CConfig::SetSmoothNumGrid(bool val_smoothnumgrid) { SmoothNumGrid = val_smoothnumgrid; }

inline unsigned short CConfig::GetKind_Turb_Model(void) { return Kind_Turb_Model; }

inline bool CConfig::GetHold_GridFixed(void) { return Hold_GridFixed; }

inline bool CConfig::GetCGNS_To_SU2(void) {return CGNS_To_SU2; }

inline bool CConfig::GetWrite_Converted_Mesh(void) { return Write_Converted_Mesh; }

inline bool CConfig::GetWrt_Vol_Sol(void) { return Wrt_Vol_Sol; }

inline bool CConfig::GetWrt_Srf_Sol(void) { return Wrt_Srf_Sol; }

inline bool CConfig::GetWrt_Csv_Sol(void) { return Wrt_Csv_Sol; }

inline bool CConfig::GetWrt_Restart(void) { return Wrt_Restart; }

inline bool CConfig::GetWrt_Residuals(void) { return Wrt_Residuals; }

inline bool CConfig::GetWrt_Halo(void) { return Wrt_Halo; }

inline bool CConfig::GetWrt_Sectional_Forces(void) { return Wrt_Sectional_Forces; }
