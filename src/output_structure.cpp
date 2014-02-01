/*!
 * \file output_structure.cpp
 * \brief Main subroutines for output solver information.
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

#include "../include/output_structure.hpp"


COutput::COutput(void) {
  
  /*--- Initialize point and connectivity counters to zero. ---*/
  nGlobal_Poin      = 0;
  nSurf_Poin        = 0;
  nGlobal_Elem      = 0;
  nSurf_Elem        = 0;
  nGlobal_Tria      = 0;
  nGlobal_Quad      = 0;
  nGlobal_Tetr      = 0;
  nGlobal_Hexa      = 0;
  nGlobal_Wedg      = 0;
  nGlobal_Pyra      = 0;
  nGlobal_Line      = 0;
  nGlobal_BoundTria = 0;
  nGlobal_BoundQuad = 0;
  
  /*--- Initialize CGNS write flag ---*/
  wrote_base_file = false;
  
  /*--- Initialize CGNS write flag ---*/
  wrote_CGNS_base = false;
  
  /*--- Initialize Tecplot write flag ---*/
  wrote_Tecplot_base = false;
  
  /*--- Initialize Paraview write flag ---*/
  wrote_Paraview_base = false;
  
}

COutput::~COutput(void) { }

void COutput::SetSurfaceCSV_Flow(CConfig *config, CGeometry *geometry,
                                 CSolver *FlowSolver, unsigned long iExtIter) {
  
  unsigned short iMarker;
  unsigned long iPoint, iVertex, Global_Index;
  double PressCoeff = 0.0, SkinFrictionCoeff;
  double xCoord, yCoord, zCoord, Mach, Pressure;
  char cstr[200];
  
  unsigned short solver = config->GetKind_Solver();
  unsigned short nDim = geometry->GetnDim();
  
  char buffer [50];
  ofstream SurfFlow_file;
  
  
  /*--- Write file name with extension if unsteady ---*/
  strcpy (cstr, config->GetSurfFlowCoeff_FileName().c_str());
  sprintf (buffer, ".csv");
  
  strcat (cstr, buffer);
  SurfFlow_file.precision(15);
  SurfFlow_file.open(cstr, ios::out);
  
  SurfFlow_file << "\"Global_Index\", \"x_coord\", \"y_coord\", ";
  if (nDim == 3) SurfFlow_file << "\"z_coord\", ";
  SurfFlow_file << "\"Pressure\", \"Pressure_Coefficient\", ";
  
  switch (solver) {
    case EULER : SurfFlow_file <<  "\"Mach_Number\"" << endl; break;
    case NAVIER_STOKES: case RANS: SurfFlow_file <<  "\"Skin_Friction_Coefficient\"" << endl; break;
  }
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_Plotting(iMarker) == YES) {
      for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        Global_Index = geometry->node[iPoint]->GetGlobalIndex();
        xCoord = geometry->node[iPoint]->GetCoord(0);
        yCoord = geometry->node[iPoint]->GetCoord(1);
        if (nDim == 3) zCoord = geometry->node[iPoint]->GetCoord(2);
        Pressure = FlowSolver->node[iPoint]->GetPressure();
        PressCoeff = FlowSolver->GetCPressure(iMarker,iVertex);
        SurfFlow_file << scientific << Global_Index << ", " << xCoord << ", " << yCoord << ", ";
        if (nDim == 3) SurfFlow_file << scientific << zCoord << ", ";
        SurfFlow_file << scientific << Pressure << ", " << PressCoeff << ", ";
        switch (solver) {
          case EULER :
            Mach = sqrt(FlowSolver->node[iPoint]->GetVelocity2()) / FlowSolver->node[iPoint]->GetSoundSpeed();
            SurfFlow_file << scientific << Mach << endl;
            break;
          case NAVIER_STOKES: case RANS:
            SkinFrictionCoeff = FlowSolver->GetCSkinFriction(iMarker,iVertex);
            SurfFlow_file << scientific << SkinFrictionCoeff << endl;
            break;
        }
      }
    }
  }
  
  SurfFlow_file.close();
  
}

void COutput::MergeConnectivity(CConfig *config, CGeometry *geometry) {
    
  /*--- Merge connectivity for each type of element (excluding halos). Note
   that we only need to merge the connectivity once, as it does not change
   during computation. Check whether the base file has been written. ---*/
  
  if (!wrote_base_file) {
    
    /*--- Merge volumetric grid. ---*/
    
    MergeVolumetricConnectivity(config, geometry, TRIANGLE    );
    
    MergeVolumetricConnectivity(config, geometry, RECTANGLE   );
    
    MergeVolumetricConnectivity(config, geometry, TETRAHEDRON );
    
    MergeVolumetricConnectivity(config, geometry, HEXAHEDRON  );
    
    MergeVolumetricConnectivity(config, geometry, WEDGE       );
    
    MergeVolumetricConnectivity(config, geometry, PYRAMID     );
    
    /*--- Merge surface grid. ---*/
    
    MergeSurfaceConnectivity(config, geometry, LINE);
    
    MergeSurfaceConnectivity(config, geometry, TRIANGLE);
    
    MergeSurfaceConnectivity(config, geometry, RECTANGLE);
    
    /*--- Update total number of volume elements after merge. ---*/
    
    nGlobal_Elem = nGlobal_Tria + nGlobal_Quad + nGlobal_Tetr +
    nGlobal_Hexa + nGlobal_Pyra + nGlobal_Wedg;
    
    /*--- Update total number of surface elements after merge. ---*/
    
    nSurf_Elem = nGlobal_Line + nGlobal_BoundTria + nGlobal_BoundQuad;
        
  }
}

void COutput::MergeCoordinates(CConfig *config, CGeometry *geometry) {
  
  /*--- Local variables needed on all processors ---*/
  
  unsigned short iDim, nDim = geometry->GetnDim();
  unsigned long iPoint, jPoint;
  
  /*--- In serial, the single process has access to all geometry, so simply
   load the coordinates into the data structure. ---*/
  
  /*--- Total number of points in the mesh (excluding halos). ---*/
  nGlobal_Poin = geometry->GetnPointDomain();
  nGlobal_Doma = geometry->GetnPointDomain();
  
  /*--- Allocate the coordinates data structure. ---*/
  
  Coords = new double*[nDim];
  for (iDim = 0; iDim < nDim; iDim++) {
    Coords[iDim] = new double[nGlobal_Poin];
  }
  
  /*--- Loop over the mesh to collect the coords of the local points. ---*/
  
  jPoint = 0;
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Retrieve the current coordinates at this node. ---*/
      for (iDim = 0; iDim < nDim; iDim++) {
        Coords[iDim][jPoint] = geometry->node[iPoint]->GetCoord(iDim);
      }
      
      /*--- Increment a counter since we may be skipping over
       some halo nodes during this loop. ---*/
      jPoint++;
    }
  }
  
}

void COutput::MergeVolumetricConnectivity(CConfig *config, CGeometry *geometry, unsigned short Elem_Type) {
  
  /*--- Local variables needed on all processors ---*/
  
  unsigned short NODES_PER_ELEMENT;
  
  unsigned long iPoint, iNode, jNode;
  unsigned long iElem = 0, jElem = 0;
  unsigned long nLocalElem = 0, nElem_Total = 0;
  
  int *Conn_Elem;
  
  /*--- Store the local number of this element type and the number of nodes
   per this element type. In serial, this will be the total number of this
   element type in the entire mesh. In parallel, it is the number on only
   the current partition. ---*/
  
  switch (Elem_Type) {
    case TRIANGLE:
      nLocalElem = geometry->GetnElemTria();
      NODES_PER_ELEMENT = N_POINTS_TRIANGLE;
      break;
    case RECTANGLE:
      nLocalElem = geometry->GetnElemQuad();
      NODES_PER_ELEMENT = N_POINTS_QUADRILATERAL;
      break;
    case TETRAHEDRON:
      nLocalElem = geometry->GetnElemTetr();
      NODES_PER_ELEMENT = N_POINTS_TETRAHEDRON;
      break;
    case HEXAHEDRON:
      nLocalElem = geometry->GetnElemHexa();
      NODES_PER_ELEMENT = N_POINTS_HEXAHEDRON;
      break;
    case WEDGE:
      nLocalElem = geometry->GetnElemWedg();
      NODES_PER_ELEMENT = N_POINTS_WEDGE;
      break;
    case PYRAMID:
      nLocalElem = geometry->GetnElemPyra();
      NODES_PER_ELEMENT = N_POINTS_PYRAMID;
      break;
    default:
      cout << "Error: Unrecognized element type \n";
      exit(0); break;
  }
  
  /*--- Merge the connectivity in serial or parallel. ---*/
  
  /*--- In serial, the single process has access to all connectivity,
   so simply load it into the data structure. ---*/
  
  /*--- Allocate a temporary array for the connectivity ---*/
  Conn_Elem = new int[nLocalElem*NODES_PER_ELEMENT];
  
  /*--- Load all elements of the current type into the buffer
   to be sent to the master node. ---*/
  jNode = 0; jElem = 0; nElem_Total = 0; bool isHalo;
  for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
    if(geometry->elem[iElem]->GetVTK_Type() == Elem_Type) {
      
      /*--- Check if this is a halo node. ---*/
      isHalo = false;
      for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {
        iPoint = geometry->elem[iElem]->GetNode(iNode);
        if (!geometry->node[iPoint]->GetDomain())
        isHalo = true;
      }
      
      /*--- Loop over all nodes in this element and load the
       connectivity into the temporary array. Do not merge any
       halo cells (periodic BC). Note that we are adding one to
       the index value because CGNS/Tecplot use 1-based indexing. ---*/
      
      if (!isHalo) {
        nElem_Total++;
        for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {
          Conn_Elem[jNode] = (int)geometry->elem[iElem]->GetNode(iNode) + 1;
          
          /*--- Increment jNode as the counter. ---*/
          jNode++;
        }
      }
    }
  }
  
  /*--- Store the particular global element count in the class data,
   and set the class data pointer to the connectivity array. ---*/
  
    switch (Elem_Type) {
      case TRIANGLE:
        nGlobal_Tria = nElem_Total;
        if (nGlobal_Tria > 0) Conn_Tria = Conn_Elem;
        break;
      case RECTANGLE:
        nGlobal_Quad = nElem_Total;
        if (nGlobal_Quad > 0) Conn_Quad = Conn_Elem;
        break;
      case TETRAHEDRON:
        nGlobal_Tetr = nElem_Total;
        if (nGlobal_Tetr > 0) Conn_Tetr = Conn_Elem;
        break;
      case HEXAHEDRON:
        nGlobal_Hexa = nElem_Total;
        if (nGlobal_Hexa > 0) Conn_Hexa = Conn_Elem;
        break;
      case WEDGE:
        nGlobal_Wedg = nElem_Total;
        if (nGlobal_Wedg > 0) Conn_Wedg = Conn_Elem;
        break;
      case PYRAMID:
        nGlobal_Pyra = nElem_Total;
        if (nGlobal_Pyra > 0) Conn_Pyra = Conn_Elem;
        break;
      default:
        cout << "Error: Unrecognized element type \n";
        exit(0); break;
    }
  
}

void COutput::MergeSurfaceConnectivity(CConfig *config, CGeometry *geometry, unsigned short Elem_Type) {
  
  /*--- Local variables needed on all processors ---*/
  
  unsigned short NODES_PER_ELEMENT;
  
  unsigned short iMarker;
  unsigned long iPoint, iNode, jNode;
  unsigned long iElem = 0, jElem = 0;
  unsigned long nLocalElem = 0, nElem_Total = 0;
  
  int *Conn_Elem;
  
  /*--- Store the local number of this element type and the number of nodes
   per this element type. In serial, this will be the total number of this
   element type in the entire mesh. In parallel, it is the number on only
   the current partition. ---*/
  
  nLocalElem = 0;
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_Plotting(iMarker) == YES) {
      for (iElem = 0; iElem < geometry->GetnElem_Bound(iMarker); iElem++) {
        if (geometry->bound[iMarker][iElem]->GetVTK_Type() == Elem_Type) {
          nLocalElem++;
        }
      }
    }
  }
  
  switch (Elem_Type) {
    case LINE:
      NODES_PER_ELEMENT = N_POINTS_LINE;
      break;
    case TRIANGLE:
      NODES_PER_ELEMENT = N_POINTS_TRIANGLE;
      break;
    case RECTANGLE:
      NODES_PER_ELEMENT = N_POINTS_QUADRILATERAL;
      break;
    default:
      cout << "Error: Unrecognized element type \n";
      exit(0); break;
  }
  
  /*--- Merge the connectivity in serial or parallel. ---*/
  
  /*--- In serial, the single process has access to all connectivity,
   so simply load it into the data structure. ---*/
  
  /*--- Allocate a temporary array for the connectivity ---*/
  Conn_Elem = new int[nLocalElem*NODES_PER_ELEMENT];
  
  /*--- Load all elements of the current type into the buffer
   to be sent to the master node. ---*/
  jNode = 0; jElem = 0; nElem_Total = 0; bool isHalo;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
  if (config->GetMarker_All_Plotting(iMarker) == YES)
  for (iElem = 0; iElem < geometry->GetnElem_Bound(iMarker); iElem++) {
    
    if (geometry->bound[iMarker][iElem]->GetVTK_Type() == Elem_Type) {
      
      /*--- Check if this is a halo node. ---*/
      isHalo = false;
      for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {
        iPoint = geometry->bound[iMarker][iElem]->GetNode(iNode);
        if (!geometry->node[iPoint]->GetDomain())
        isHalo = true;
      }
      
      /*--- Loop over all nodes in this element and load the
       connectivity into the temporary array. Do not merge any
       halo cells (periodic BC). Note that we are adding one to
       the index value because CGNS/Tecplot use 1-based indexing. ---*/
      if (!isHalo) {
        nElem_Total++;
        for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {
          Conn_Elem[jNode] = (int)geometry->bound[iMarker][iElem]->GetNode(iNode) + 1;
          
          /*--- Increment jNode as the counter. ---*/
          jNode++;
        }
      }
    }
  }
  
  /*--- Store the particular global element count in the class data,
   and set the class data pointer to the connectivity array. ---*/
  
    switch (Elem_Type) {
      case LINE:
        nGlobal_Line = nElem_Total;
        if (nGlobal_Line > 0) Conn_Line = Conn_Elem;
        break;
      case TRIANGLE:
        nGlobal_BoundTria = nElem_Total;
        if (nGlobal_BoundTria > 0) Conn_BoundTria = Conn_Elem;
        break;
      case RECTANGLE:
        nGlobal_BoundQuad = nElem_Total;
        if (nGlobal_BoundQuad > 0) Conn_BoundQuad = Conn_Elem;
        break;
      default:
        cout << "Error: Unrecognized element type \n";
        exit(0); break;
    }
  
}

void COutput::MergeSolution(CConfig *config, CGeometry *geometry, CSolver **solver) {
  
  /*--- Local variables needed on all processors ---*/
  unsigned short Kind_Solver  = config->GetKind_Solver();
  unsigned short iVar, jVar, FirstIndex = NONE, SecondIndex = NONE, ThirdIndex = NONE;
  unsigned short nVar_First = 0, nVar_Second = 0, nVar_Third = 0, iVar_Eddy = 0, iVar_Sharp = 0;
  unsigned short iVar_PressMach = 0, iVar_TempLam = 0,
  iVar_ViscCoeffs = 0, iVar_Extra = 0;
  
  unsigned long iPoint = 0, jPoint = 0, iVertex = 0, iMarker = 0;
  double RefDensity, RefPressure, factor;
  
  double *Aux_Frict, *Aux_Heat, *Aux_yPlus;
  
  bool flow           = (( config->GetKind_Solver() == EULER             ) ||
                         ( config->GetKind_Solver() == NAVIER_STOKES     ) ||
                         ( config->GetKind_Solver() == RANS              )   );
  
  unsigned short iDim;
  bool nDim               = geometry->GetnDim();
  double RefAreaCoeff     = config->GetRefAreaCoeff();
  double RefVel2;
  
  /*--- Set the non-dimensionalization ---*/
  if (flow) {
      RefVel2 = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
      RefVel2  += solver[FLOW_SOL]->GetVelocity_Inf(iDim)*solver[FLOW_SOL]->GetVelocity_Inf(iDim);
    RefDensity  = solver[FLOW_SOL]->GetDensity_Inf();
    RefPressure = solver[FLOW_SOL]->GetPressure_Inf();
    factor = 1.0 / (0.5*RefDensity*RefAreaCoeff*RefVel2);
  }
  
  /*--- Prepare send buffers for the conservative variables. Need to
   find the total number of conservative variables and also the
   index for their particular solution container. ---*/
  switch (Kind_Solver) {
    case EULER : case NAVIER_STOKES:
      FirstIndex = FLOW_SOL; SecondIndex = NONE; ThirdIndex = NONE;
      break;
    case RANS :
      FirstIndex = FLOW_SOL; SecondIndex = TURB_SOL; ThirdIndex = NONE;
      break;
  }
  nVar_First = solver[FirstIndex]->GetnVar();
  if (SecondIndex != NONE) nVar_Second = solver[SecondIndex]->GetnVar();
  if (ThirdIndex != NONE) nVar_Third = solver[ThirdIndex]->GetnVar();
  nVar_Consv = nVar_First + nVar_Second + nVar_Third;
  
  if (config->GetWrt_Residuals()) nVar_Total = 2*nVar_Consv;
  else nVar_Total = nVar_Consv;
  
  if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
    /*--- Pressure, Cp, Mach ---*/
    iVar_PressMach = nVar_Total;
    nVar_Total += 3;
  }
  
  if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
    /*--- Temperature, Laminar Viscosity ---*/
    iVar_TempLam = nVar_Total;
    nVar_Total += 2;
    /*--- Skin Friction, Heat Flux, & yPlus ---*/
    iVar_ViscCoeffs = nVar_Total;
    nVar_Total += 3;
  }
  
  if (Kind_Solver == RANS) {
    /*--- Eddy Viscosity ---*/
    iVar_Eddy = nVar_Total;
    nVar_Total += 1;
  }
  
  if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
    /*--- Sharp edges ---*/
    iVar_Sharp = nVar_Total;
    nVar_Total += 1;
  }
  
  if (config->GetExtraOutput()) {
    if (Kind_Solver == RANS) {
      iVar_Extra  = nVar_Total;
      nVar_Extra  = solver[TURB_SOL]->GetnOutputVariables();
      nVar_Total += nVar_Extra;
    }
  }
  
  /*--- Merge the solution either in serial or parallel. ---*/
  
  /*--- In serial, the single process has access to all solution data,
   so it is simple to retrieve and store inside Solution_Data. ---*/
  
  nGlobal_Poin = geometry->GetnPointDomain();
  Data = new double*[nVar_Total];
  for (iVar = 0; iVar < nVar_Total; iVar++) {
    Data[iVar] = new double[nGlobal_Poin];
  }
  
  /*--- First, loop through the mesh in order to find and store the
   value of some coefficients at any surface nodes. They
   will be placed in an auxiliary vector and then communicated like
   all other volumetric variables. ---*/
  
  if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
    
    Aux_Frict = new double [geometry->GetnPointDomain()];
    Aux_Heat  = new double [geometry->GetnPointDomain()];
    Aux_yPlus = new double [geometry->GetnPointDomain()];
    
    for(iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
      Aux_Frict[iPoint] = 0.0;
      Aux_Heat[iPoint]  = 0.0;
      Aux_yPlus[iPoint] = 0.0;
    }
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if (config->GetMarker_All_Plotting(iMarker) == YES) {
      for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        Aux_Frict[iPoint] = solver[FLOW_SOL]->GetCSkinFriction(iMarker,iVertex);
        Aux_Heat[iPoint]  = solver[FLOW_SOL]->GetHeatTransferCoeff(iMarker,iVertex);
        Aux_yPlus[iPoint] = solver[FLOW_SOL]->GetYPlus(iMarker,iVertex);
      }
    }
  }
  
  /*--- Loop over all points in the mesh, but only write data
   for nodes in the domain (ignore periodic halo nodes). ---*/
  
  jPoint = 0;
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    
    /*--- Check for halo nodes & only write if requested ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Solution (first, second and third system of equations) ---*/
      jVar = 0;
      for (iVar = 0; iVar < nVar_First; iVar++) {
        Data[jVar][jPoint] = solver[FirstIndex]->node[iPoint]->GetSolution(iVar);
        jVar++;
      }
      
      for (iVar = 0; iVar < nVar_Second; iVar++) {
        Data[jVar][jPoint] = solver[SecondIndex]->node[iPoint]->GetSolution(iVar);
        jVar++;
      }
      
      for (iVar = 0; iVar < nVar_Third; iVar++) {
        Data[jVar][jPoint] = solver[ThirdIndex]->node[iPoint]->GetSolution(iVar);
        jVar++;
      }
      
      /*--- Residual (first, second, and third system of equations) ---*/
      if (config->GetWrt_Residuals()) {
        for (iVar = 0; iVar < nVar_First; iVar++) {
          Data[jVar][jPoint] = solver[FirstIndex]->LinSysRes.GetBlock(iPoint, iVar);
          jVar++;
        }
        
        for (iVar = 0; iVar < nVar_Second; iVar++) {
          Data[jVar][jPoint] = solver[SecondIndex]->LinSysRes.GetBlock(iPoint, iVar);
          jVar++;
        }
        
        for (iVar = 0; iVar < nVar_Third; iVar++) {
          Data[jVar][jPoint] = solver[ThirdIndex]->LinSysRes.GetBlock(iPoint, iVar);
          jVar++;
        }
      }
      
      /*--- Any extra data for output files ---*/
      
      switch (Kind_Solver) {
        case EULER:
          
          /*--- Load buffers with the pressure, Cp, and mach variables. ---*/
            Data[jVar][jPoint] = solver[FLOW_SOL]->node[iPoint]->GetPressure(); jVar++;
            Data[jVar][jPoint] = (solver[FLOW_SOL]->node[iPoint]->GetPressure() - RefPressure)*factor*RefAreaCoeff; jVar++;
            Data[jVar][jPoint] = sqrt(solver[FLOW_SOL]->node[iPoint]->GetVelocity2())/solver[FLOW_SOL]->node[iPoint]->GetSoundSpeed(); jVar++;
          Data[jVar][jPoint] = geometry->node[iPoint]->GetSharpEdge_Distance(); jVar++;
          break;
          /*--- Write pressure, Cp, mach, temperature, laminar viscosity, skin friction, heat transfer, yplus ---*/
        case NAVIER_STOKES:
            Data[jVar][jPoint] = solver[FLOW_SOL]->node[iPoint]->GetPressure(); jVar++;
            Data[jVar][jPoint] = (solver[FLOW_SOL]->node[iPoint]->GetPressure() - RefPressure)*factor*RefAreaCoeff; jVar++;
            Data[jVar][jPoint] = sqrt(solver[FLOW_SOL]->node[iPoint]->GetVelocity2())/
            solver[FLOW_SOL]->node[iPoint]->GetSoundSpeed(); jVar++;
            Data[jVar][jPoint] = solver[FLOW_SOL]->node[iPoint]->GetTemperature(); jVar++;
            Data[jVar][jPoint] = solver[FLOW_SOL]->node[iPoint]->GetLaminarViscosity(); jVar++;
            Data[jVar][jPoint] = Aux_Frict[iPoint]; jVar++;
            Data[jVar][jPoint] = Aux_Heat[iPoint];  jVar++;
            Data[jVar][jPoint] = Aux_yPlus[iPoint]; jVar++;
          Data[jVar][jPoint] = geometry->node[iPoint]->GetSharpEdge_Distance(); jVar++;
          break;
          /*--- Write pressure, Cp, mach, temperature, laminar viscosity, skin friction, heat transfer, yplus, eddy viscosity ---*/
        case RANS:
            Data[jVar][jPoint] = solver[FLOW_SOL]->node[iPoint]->GetPressure(); jVar++;
            Data[jVar][jPoint] = (solver[FLOW_SOL]->node[iPoint]->GetPressure() - RefPressure)*factor*RefAreaCoeff; jVar++;
            Data[jVar][jPoint] = sqrt(solver[FLOW_SOL]->node[iPoint]->GetVelocity2())/
            solver[FLOW_SOL]->node[iPoint]->GetSoundSpeed(); jVar++;
            Data[jVar][jPoint] = solver[FLOW_SOL]->node[iPoint]->GetTemperature(); jVar++;
            Data[jVar][jPoint] = solver[FLOW_SOL]->node[iPoint]->GetLaminarViscosity(); jVar++;
            Data[jVar][jPoint] = Aux_Frict[iPoint]; jVar++;
            Data[jVar][jPoint] = Aux_Heat[iPoint];  jVar++;
            Data[jVar][jPoint] = Aux_yPlus[iPoint]; jVar++;
            Data[jVar][jPoint] = solver[FLOW_SOL]->node[iPoint]->GetEddyViscosity(); jVar++;
          Data[jVar][jPoint] = geometry->node[iPoint]->GetSharpEdge_Distance(); jVar++;
          break;
      }
    }
    
    if (config->GetExtraOutput()) {
      if (Kind_Solver == RANS) {
        for (unsigned short iVar = 0; iVar < nVar_Extra; iVar++) {
          Data[jVar][jPoint] =  solver[TURB_SOL]->OutputVariables[iPoint*nVar_Extra+iVar];
          jVar++;
        }
      }
    }
    
    /*--- Increment jPoint as the counter. We need this because iPoint
     may include halo nodes that we skip over during this loop. ---*/
    
    jPoint++;
    
  }
  
  /*--- Release memory needed for surface coefficients ---*/
  
  if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
    delete [] Aux_Frict; delete [] Aux_Heat; delete [] Aux_yPlus;
  }
  
}

void COutput::SetRestart(CConfig *config, CGeometry *geometry) {
  
  /*--- Local variables ---*/
  
  unsigned short Kind_Solver  = config->GetKind_Solver();
  unsigned short iVar, iDim, nDim = geometry->GetnDim();
  unsigned long iPoint;
  ofstream restart_file;
  string filename;
  
  /*--- Retrieve filename from config ---*/
  if (config->GetAdjoint()) {
    filename = config->GetRestart_AdjFileName();
    filename = config->GetObjFunc_Extension(filename);
  } else {
    filename = config->GetRestart_FlowFileName();
  }
  
  /*--- Open the restart file and write the solution. ---*/
  restart_file.open(filename.c_str(), ios::out);
  restart_file.precision(15);
  
  /*--- Write the header line based on the particular solver ----*/
  restart_file << "\"PointID\"";
  
  /*--- Mesh coordinates are always written to the restart first ---*/
  if (nDim == 2) {
    restart_file << "\t\"x\"\t\"y\"";
  } else {
    restart_file << "\t\"x\"\t\"y\"\t\"z\"";
  }
  
  for (iVar = 0; iVar < nVar_Consv; iVar++) {
    restart_file << "\t\"Conservative_" << iVar+1<<"\"";
  }
  if (config->GetWrt_Residuals()) {
    for (iVar = 0; iVar < nVar_Consv; iVar++) {
      restart_file << "\t\"Residual_" << iVar+1<<"\"";
    }
  }
  
  if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
    restart_file << "\t\"Pressure\"\t\"Pressure_Coefficient\"\t\"Mach\"";
  }
  
  if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
    restart_file << "\t\"Temperature\"\t\"Laminar_Viscosity\"\t\"Skin_Friction_Coefficient\"\t\"Heat_Transfer\"\t\"Y_Plus\"";
  }
  
  if (Kind_Solver == RANS) {
    restart_file << "\t\"Eddy_Viscosity\"";
  }
  
  if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
    restart_file << "\t\"Sharp_Edge_Dist\"";
  }
  
  restart_file << endl;
  
  /*--- Write the restart file ---*/
  
  for (iPoint = 0; iPoint < geometry->GetGlobal_nPointDomain(); iPoint++) {
    
    /*--- Index of the point ---*/
    
    restart_file << iPoint << "\t";
    
    /*--- Write the grid coordinates first ---*/
    
    for (iDim = 0; iDim < nDim; iDim++) {
      restart_file << scientific << Coords[iDim][iPoint] << "\t";
    }
    
    /*--- Loop over the variables and write the values to file ---*/
    for (iVar = 0; iVar < nVar_Total; iVar++) {
      restart_file << scientific << Data[iVar][iPoint] << "\t";
    }
    
    /*--- End the line of the restart file for this point ---*/
    
    restart_file << endl;
    
  }
  
  restart_file.close();
  
}

void COutput::DeallocateCoordinates(CConfig *config, CGeometry *geometry) {
  
  /*--- Local variables and initialization ---*/
  
  unsigned short iDim, nDim = geometry->GetnDim();
  
  /*--- The master node alone owns all data found in this routine. ---*/
  
    /*--- Deallocate memory for coordinate data ---*/
    for (iDim = 0; iDim < nDim; iDim++) {
      delete [] Coords[iDim];
    }
    delete [] Coords;
    
}

void COutput::DeallocateConnectivity(CConfig *config, CGeometry *geometry, bool surf_sol) {
  
    /*--- Deallocate memory for connectivity data ---*/
    if (surf_sol) {
      if (nGlobal_Line > 0) delete [] Conn_Line;
      if (nGlobal_BoundTria > 0) delete [] Conn_BoundTria;
      if (nGlobal_BoundQuad > 0) delete [] Conn_BoundQuad;
    }
    else {
      if (nGlobal_Tria > 0) delete [] Conn_Tria;
      if (nGlobal_Quad > 0) delete [] Conn_Quad;
      if (nGlobal_Tetr > 0) delete [] Conn_Tetr;
      if (nGlobal_Hexa > 0) delete [] Conn_Hexa;
      if (nGlobal_Wedg > 0) delete [] Conn_Wedg;
      if (nGlobal_Pyra > 0) delete [] Conn_Pyra;
    }
    
}

void COutput::DeallocateSolution(CConfig *config, CGeometry *geometry) {
  
    /*--- Deallocate memory for solution data ---*/
    for (unsigned short iVar = 0; iVar < nVar_Total; iVar++) {
      delete [] Data[iVar];
    }
    delete [] Data;
    
}

void COutput::SetHistory_Header(ofstream *ConvHist_file, CConfig *config) {
  
  char cstr[200], buffer[50], turb_resid[1000];
  unsigned short iMarker;
  string Monitoring_Tag, Monitoring_coeff;
  
  bool rotating_frame = config->GetRotating_Frame();
  bool turbulent = config->GetKind_Solver() == RANS;
  
  bool isothermal = false;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
  if (config->GetMarker_All_Boundary(iMarker) == ISOTHERMAL) isothermal = true;
  
  /*--- Write file name with extension ---*/
  string filename = config->GetConv_FileName();
  strcpy (cstr, filename.data());
  
  if ((config->GetOutput_FileFormat() == TECPLOT) ||
      (config->GetOutput_FileFormat() == TECPLOT_BINARY))  sprintf (buffer, ".plt");
  if ((config->GetOutput_FileFormat() == CGNS_SOL) ||
      (config->GetOutput_FileFormat() == PARAVIEW))  sprintf (buffer, ".csv");
  strcat(cstr,buffer);
  
  ConvHist_file->open(cstr, ios::out);
  ConvHist_file->precision(15);
  
  /*--- Begin of the header ---*/
  
  char begin[]= "\"Iteration\"";
  
  /*--- Header for the coefficients ---*/
  
  char flow_coeff[]= ",\"CLift\",\"CDrag\",\"CSideForce\",\"CMx\",\"CMy\",\"CMz\",\"CFx\",\"CFy\",\"CFz\",\"CL/CD\"";
  char heat_coeff[]= ",\"CHeat_Load\",\"CHeat_Max\"";
  char rotating_frame_coeff[]= ",\"CMerit\",\"CT\",\"CQ\"";
  
  /*--- Header for the residuals ---*/
  
  char flow_resid[]= ",\"Res_Flow[0]\",\"Res_Flow[1]\",\"Res_Flow[2]\",\"Res_Flow[3]\",\"Res_Flow[4]\"";
  switch (config->GetKind_Turb_Model()) {
    case SA:	sprintf (turb_resid, ",\"Res_Turb[0]\""); break;
    case SST:	sprintf (turb_resid, ",\"Res_Turb[0]\",\"Res_Turb[1]\""); break;
  }

  
  /*--- End of the header ---*/
  
  char end[]= ",\"Linear_Solver_Iterations\",\"Time(min)\"\n";
  
  if ((config->GetOutput_FileFormat() == TECPLOT) ||
      (config->GetOutput_FileFormat() == TECPLOT_BINARY)) {
    ConvHist_file[0] << "TITLE = \"SU2 Simulation\"" << endl;
    ConvHist_file[0] << "VARIABLES = ";
  }
  
  /*--- Write the header, case depending ---*/
  switch (config->GetKind_Solver()) {
      
    case EULER : case NAVIER_STOKES: case RANS :
      ConvHist_file[0] << begin << flow_coeff;
      if (isothermal) ConvHist_file[0] << heat_coeff;
      if (rotating_frame) ConvHist_file[0] << rotating_frame_coeff;
      ConvHist_file[0] << flow_resid;
      if (turbulent) ConvHist_file[0] << turb_resid;
      ConvHist_file[0] << end;
      break;
      
  }
  
  if (config->GetOutput_FileFormat() == TECPLOT || config->GetOutput_FileFormat() == TECPLOT_BINARY) {
    ConvHist_file[0] << "ZONE T= \"Convergence history\"" << endl;
  }
  
}


void COutput::SetConvergence_History(ofstream *ConvHist_file, CGeometry **geometry, CSolver ***solver_container, CConfig *config, CIntegration **integration, bool DualTime_Iteration, double timeused) {
  
    unsigned long iIntIter = config->GetIntIter();
    unsigned long iExtIter = config->GetExtIter();
    
    /*--- WARNING: These buffers have hard-coded lengths. Note that you
     may have to adjust them to be larger if adding more entries. ---*/
    char begin[1000], direct_coeff[1000], flow_resid[1000], turb_resid[1000], end[1000];
    double dummy = 0.0;
    unsigned short iVar;
    
    unsigned long LinSolvIter = 0;
    double timeiter = timeused/double(iExtIter+1);
    
    unsigned short FinestMesh = config->GetFinestMesh();
    unsigned short nDim = geometry[FinestMesh]->GetnDim();
    
    bool turbulent = ((config->GetKind_Solver() == RANS));
    bool flow = (config->GetKind_Regime() == EULER) || (config->GetKind_Regime() == NAVIER_STOKES) || (config->GetKind_Regime() == RANS);
    
    /*--- Initialize variables to store information from all domains (direct solution) ---*/
    double Total_CLift = 0.0, Total_CDrag = 0.0, Total_CSideForce = 0.0, Total_CMx = 0.0, Total_CMy = 0.0, Total_CMz = 0.0, Total_CEff = 0.0,
    Total_CFx = 0.0, Total_CFy = 0.0, Total_CFz = 0.0;
    
    /*--- Residual arrays ---*/
    double *residual_flow = NULL, *residual_turbulent = NULL, *residual_transition = NULL, *residual_TNE2 = NULL, *residual_levelset = NULL;
    double *residual_adjflow = NULL, *residual_adjturbulent = NULL, *residual_adjTNE2 = NULL, *residual_adjlevelset = NULL;
    double *residual_wave = NULL; double *residual_fea = NULL; double *residual_heat = NULL;
    
    /*--- Coefficients Monitored arrays ---*/
    double *aeroelastic_plunge = NULL, *aeroelastic_pitch = NULL, *Surface_CLift = NULL, *Surface_CDrag = NULL, *Surface_CMx = NULL, *Surface_CMy = NULL, *Surface_CMz = NULL;
    
    /*--- Initialize number of variables ---*/
    unsigned short nVar_Flow = 0, nVar_LevelSet = 0, nVar_Turb = 0, nVar_Trans = 0, nVar_TNE2 = 0, nVar_Wave = 0, nVar_Heat = 0, nVar_FEA = 0,     nVar_AdjFlow = 0, nVar_AdjTNE2 = 0, nVar_AdjLevelSet = 0, nVar_AdjTurb = 0;
    
    /*--- Direct problem variables ---*/
    nVar_Flow = nDim+2;
    if (turbulent) {
      switch (config->GetKind_Turb_Model()){
        case SA:	nVar_Turb = 1; break;
        case SST: nVar_Turb = 2; break;
      }
    }
    
    /*--- Adjoint problem variables ---*/
    nVar_AdjFlow = nDim+2;
    if (turbulent) {
      switch (config->GetKind_Turb_Model()){
        case SA:	nVar_AdjTurb = 1; break;
        case SST: nVar_AdjTurb = 2; break;
      }
    }
    
    /*--- Allocate memory for the residual ---*/
    residual_flow       = new double[nVar_Flow];
    residual_turbulent  = new double[nVar_Turb];
    residual_transition = new double[nVar_Trans];
    residual_TNE2       = new double[nVar_TNE2];
    residual_levelset   = new double[nVar_LevelSet];
    residual_wave       = new double[nVar_Wave];
    residual_fea        = new double[nVar_FEA];
    residual_heat       = new double[nVar_Heat];
    
    residual_adjflow      = new double[nVar_AdjFlow];
    residual_adjturbulent = new double[nVar_AdjTurb];
    residual_adjTNE2      = new double[nVar_AdjTNE2];
    residual_adjlevelset  = new double[nVar_AdjLevelSet];
    
    /*--- Allocate memory for the coefficients being monitored ---*/
    aeroelastic_plunge = new double[config->GetnMarker_Monitoring()];
    aeroelastic_pitch  = new double[config->GetnMarker_Monitoring()];
    Surface_CLift      = new double[config->GetnMarker_Monitoring()];
    Surface_CDrag      = new double[config->GetnMarker_Monitoring()];
    Surface_CMx        = new double[config->GetnMarker_Monitoring()];
    Surface_CMy        = new double[config->GetnMarker_Monitoring()];
    Surface_CMz        = new double[config->GetnMarker_Monitoring()];
    
    /*--- Write information from nodes ---*/
    switch (config->GetKind_Solver()) {
        
      case EULER:                   case NAVIER_STOKES:                   case RANS:
        
        /*--- Flow solution coefficients ---*/
        Total_CLift       = solver_container[FinestMesh][FLOW_SOL]->GetTotal_CLift();
        Total_CDrag       = solver_container[FinestMesh][FLOW_SOL]->GetTotal_CDrag();
        Total_CSideForce  = solver_container[FinestMesh][FLOW_SOL]->GetTotal_CSideForce();
        Total_CEff        = solver_container[FinestMesh][FLOW_SOL]->GetTotal_CEff();
        Total_CMx         = solver_container[FinestMesh][FLOW_SOL]->GetTotal_CMx();
        Total_CMy         = solver_container[FinestMesh][FLOW_SOL]->GetTotal_CMy();
        Total_CMz         = solver_container[FinestMesh][FLOW_SOL]->GetTotal_CMz();
        Total_CFx         = solver_container[FinestMesh][FLOW_SOL]->GetTotal_CFx();
        Total_CFy         = solver_container[FinestMesh][FLOW_SOL]->GetTotal_CFy();
        Total_CFz         = solver_container[FinestMesh][FLOW_SOL]->GetTotal_CFz();
        
        /*--- Flow Residuals ---*/
        
        for (iVar = 0; iVar < nVar_Flow; iVar++)
        residual_flow[iVar] = solver_container[FinestMesh][FLOW_SOL]->GetRes_RMS(iVar);
        
        /*--- Turbulent residual ---*/
        
        if (turbulent) {
          for (iVar = 0; iVar < nVar_Turb; iVar++)
          residual_turbulent[iVar] = solver_container[FinestMesh][TURB_SOL]->GetRes_RMS(iVar);
        }
        
        
        /*--- Iterations of the linear solver ---*/
        
        LinSolvIter = (unsigned long) solver_container[FinestMesh][FLOW_SOL]->GetIterLinSolver();

        
        break;
        
    }
    
    /*--- Header frecuency ---*/
    
    bool Unsteady = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                     (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
    bool In_NoDualTime = (!DualTime_Iteration && (iExtIter % config->GetWrt_Con_Freq() == 0));
    bool In_DualTime_0 = (DualTime_Iteration && (iIntIter % 10 == 0));
    bool In_DualTime_1 = (!DualTime_Iteration && Unsteady);
    bool In_DualTime_2 = (Unsteady && DualTime_Iteration && (iExtIter % config->GetWrt_Con_Freq() == 0));
    bool In_DualTime_3 = (Unsteady && !DualTime_Iteration && (iExtIter % config->GetWrt_Con_Freq() == 0));
    
    bool write_heads;
    if (Unsteady) write_heads = (iIntIter == 0);
    else write_heads = (((iExtIter % (config->GetWrt_Con_Freq()*20)) == 0));
    
    if ((In_NoDualTime || In_DualTime_0 || In_DualTime_1) && (In_NoDualTime || In_DualTime_2 || In_DualTime_3)) {
      
      /*--- Prepare the history file output, note that the dual
       time output don't write to the history file ---*/
      if (!DualTime_Iteration) {
        
        /*--- Write the begining of the history file ---*/
        sprintf (begin, "%12d", int(iExtIter));
        
        /*--- Write the end of the history file ---*/
        sprintf (end, ", %12.10f, %12.10f\n", double(LinSolvIter), timeused/60.0);
        
        /*--- Write the solution and residual of the history file ---*/
        switch (config->GetKind_Solver()) {
            
          case EULER : case NAVIER_STOKES: case RANS:
            
            /*--- Direct coefficients ---*/
            sprintf (direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f",
                     Total_CLift, Total_CDrag, Total_CSideForce, Total_CMx, Total_CMy, Total_CMz, Total_CFx, Total_CFy,
                     Total_CFz, Total_CEff);
            
            /*--- Flow residual ---*/
            if (nDim == 2) {
              sprintf (flow_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", log10 (residual_flow[0]), log10 (residual_flow[1]), log10 (residual_flow[2]), log10 (residual_flow[3]), dummy );
            }
            else {
              sprintf (flow_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", log10 (residual_flow[0]), log10 (residual_flow[1]), log10 (residual_flow[2]), log10 (residual_flow[3]), log10 (residual_flow[4]) );
            }
            
            /*--- Turbulent residual ---*/
            if (turbulent){
              switch(nVar_Turb) {
                case 1: sprintf (turb_resid, ", %12.10f", log10 (residual_turbulent[0])); break;
                case 2: sprintf (turb_resid, ", %12.10f, %12.10f", log10(residual_turbulent[0]), log10(residual_turbulent[1])); break;
              }
            }
            
            break;
            
        }
      }
      
      /*--- Write the screen header---*/
      if ((write_heads) && !(!DualTime_Iteration && Unsteady)) {
        
        if (!Unsteady) {
          switch (config->GetKind_Solver()) {
            case EULER :                  case NAVIER_STOKES:
              cout << endl << " Min Delta Time: " << solver_container[FinestMesh][FLOW_SOL]->GetMin_Delta_Time()<<
              ". Max Delta Time: " << solver_container[FinestMesh][FLOW_SOL]->GetMax_Delta_Time() << ".";
              break;
          }
        }
        else {
          if (flow) {
            cout << endl << " Min DT: " << solver_container[FinestMesh][FLOW_SOL]->GetMin_Delta_Time()<<
            ". Max DT: " << solver_container[FinestMesh][FLOW_SOL]->GetMax_Delta_Time() <<
            ". Dual Time step: " << ".";
          }
          else {
            cout << endl << " Dual Time step: "  << ".";
          }
        }
        
        switch (config->GetKind_Solver()) {
          case EULER :                  case NAVIER_STOKES:
            
            /*--- Visualize the maximum residual ---*/
            cout << endl << " Maximum residual: " << log10(solver_container[FinestMesh][FLOW_SOL]->GetRes_Max(0))
            <<", located at point "<< solver_container[FinestMesh][FLOW_SOL]->GetPoint_Max(0) << "." << endl;
            
            if (!Unsteady) cout << endl << " Iter" << "    Time(s)";
            else cout << endl << " IntIter" << " ExtIter";
            
            cout << "     Res[Rho]" << "     Res[RhoE]" << "   CLift(Total)" << "   CDrag(Total)" << endl;
            
            break;
            
          case RANS :
            
            /*--- Visualize the maximum residual ---*/
            cout << endl << " Maximum residual: " << log10(solver_container[FinestMesh][FLOW_SOL]->GetRes_Max(0))
            <<", located at point "<< solver_container[FinestMesh][FLOW_SOL]->GetPoint_Max(0) << "." << endl;
            
            if (!Unsteady) cout << endl << " Iter" << "    Time(s)";
            else cout << endl << " IntIter" << " ExtIter";
            cout << "      Res[Rho]";
            
            switch (config->GetKind_Turb_Model()){
              case SA:	cout << "       Res[nu]"; break;
              case SST:	cout << "     Res[kine]" << "     Res[omega]"; break;
            }
            
            cout << "   CLift(Total)"   << "   CDrag(Total)"   << endl;
            break;
            
            
        }
        
      }
      
      /*--- Write the solution on the screen and history file ---*/
      cout.precision(6);
      cout.setf(ios::fixed,ios::floatfield);
      
      if (!Unsteady) {
        cout.width(5); cout << iExtIter;
        cout.width(11); cout << timeiter;
        
      } else {
        cout.width(8); cout << iIntIter;
        cout.width(8); cout << iExtIter;
      }
      
      switch (config->GetKind_Solver()) {
        case EULER : case NAVIER_STOKES:
          
          if (!DualTime_Iteration) {
            ConvHist_file[0] << begin << direct_coeff << flow_resid;
            ConvHist_file[0] << end;
            ConvHist_file[0].flush();
          }
          
          cout.precision(6);
          cout.setf(ios::fixed,ios::floatfield);
          cout.width(13); cout << log10(residual_flow[0]);
              if (nDim == 2 ) { cout.width(14); cout << log10(residual_flow[3]); }
              else { cout.width(14); cout << log10(residual_flow[4]); }
          
          cout.width(15); cout << min(1000.0,max(-1000.0, Total_CLift)); cout.width(15); cout << min(1000.0,max(-1000.0, Total_CDrag));
          cout << endl;
          
          break;
          
        case RANS :
          
          if (!DualTime_Iteration) {
            ConvHist_file[0] << begin << direct_coeff << flow_resid << turb_resid << end;
            ConvHist_file[0].flush();
          }
          
          cout.precision(6);
          cout.setf(ios::fixed,ios::floatfield);
          
          cout.width(14);
          cout << log10(residual_flow[0]);
          
          switch(nVar_Turb) {
            case 1: cout.width(14); cout << log10(residual_turbulent[0]); break;
            case 2: cout.width(14); cout << log10(residual_turbulent[0]);
              cout.width(15); cout << log10(residual_turbulent[1]); break;
          }
          
          
          cout.width(15); cout << min(1000.0,max(-1000.0, Total_CLift)); cout.width(15); cout << min(1000.0,max(-1000.0, Total_CDrag));
          cout << endl;
          
          break;
          
          
      }
      cout.unsetf(ios::fixed);
      
      delete [] residual_flow;
      delete [] residual_levelset;
      delete [] residual_TNE2;
      delete [] residual_turbulent;
      delete [] residual_transition;
      delete [] residual_wave;
      delete [] residual_fea;
      delete [] residual_heat;
      
      delete [] residual_adjflow;
      delete [] residual_adjTNE2;
      delete [] residual_adjlevelset;
      delete [] residual_adjturbulent;
      
      delete [] Surface_CLift;
      delete [] Surface_CDrag;
      delete [] Surface_CMx;
      delete [] Surface_CMy;
      delete [] Surface_CMz;
      
    }
}

void COutput::SetResult_Files(CSolver ***solver_container, CGeometry **geometry, CConfig *config, unsigned long iExtIter) {
  
    /*--- Flags identifying the types of files to be written. ---*/
    bool Wrt_Vol = config->GetWrt_Vol_Sol();
    bool Wrt_Srf = config->GetWrt_Srf_Sol();
    bool Wrt_Csv = config->GetWrt_Csv_Sol();
    bool Wrt_Rst = config->GetWrt_Restart();
    
    switch (config->GetKind_Solver()) {
        
      case EULER : case NAVIER_STOKES : case RANS :
        
        if (Wrt_Csv) SetSurfaceCSV_Flow(config, geometry[MESH_0], solver_container[MESH_0][FLOW_SOL], iExtIter);
        break;
        
    }
    
    /*--- Get the file output format ---*/
    
    unsigned short FileFormat = config->GetOutput_FileFormat();
    
    /*--- Merge the node coordinates and connectivity, if necessary. This
     is only performed if a volume solution file is requested, and it
     is active by default. ---*/
    
    if (Wrt_Vol || Wrt_Srf)
    MergeConnectivity(config, geometry[MESH_0]);
    
    /*--- Merge coordinates of all grid nodes (excluding ghost points).
     The grid coordinates are always merged and included first in the
     restart files. ---*/
    
    MergeCoordinates(config, geometry[MESH_0]);
    
    /*--- Merge the solution data needed for volume solutions and restarts ---*/
    
    if (Wrt_Vol || Wrt_Rst)
    MergeSolution(config, geometry[MESH_0], solver_container[MESH_0]);
    
    /*--- Write restart, CGNS, or Tecplot files using the merged data.
     This data lives only on the master, and these routines are currently
     executed by the master proc alone (as if in serial). ---*/
  
      /*--- Write a native restart file ---*/
      if (Wrt_Rst)
        SetRestart(config, geometry[MESH_0]);
      
      if (Wrt_Vol) {
        
        switch (FileFormat) {
            
          case TECPLOT:
            
            /*--- Write a Tecplot ASCII file ---*/
            SetTecplot_ASCII(config, geometry[MESH_0], false);
            DeallocateConnectivity(config, geometry[MESH_0], false);
            break;
            
          case PARAVIEW:
            
            /*--- Write a Paraview ASCII file ---*/
            SetParaview_ASCII(config, geometry[MESH_0], false);
            DeallocateConnectivity(config, geometry[MESH_0], false);
            break;
            
          default:
            break;
        }
        
      }
      
      if (Wrt_Srf) {
        
        switch (FileFormat) {
            
          case TECPLOT:
            
            /*--- Write a Tecplot ASCII file ---*/
            SetTecplot_ASCII(config, geometry[MESH_0], true);
            DeallocateConnectivity(config, geometry[MESH_0], true);
            break;
            
          case PARAVIEW:
            
            /*--- Write a Paraview ASCII file ---*/
            SetParaview_ASCII(config, geometry[MESH_0], true);
            DeallocateConnectivity(config, geometry[MESH_0], true);
            break;
            
          default:
            break;
        }
        
      }
      
      /*--- Release memory needed for merging the solution data. ---*/
      if (((Wrt_Vol) || (Wrt_Srf)) && (FileFormat == TECPLOT ||
                                       FileFormat == TECPLOT_BINARY ||
                                       FileFormat == PARAVIEW))
      DeallocateCoordinates(config, geometry[MESH_0]);
      
      if (Wrt_Vol || Wrt_Rst)
      DeallocateSolution(config, geometry[MESH_0]);
      
  
}

void COutput::SetBaselineResult_Files(CSolver *solver, CGeometry *geometry, CConfig *config,
                                      unsigned long iExtIter) {
    
    /*--- Flags identifying the types of files to be written. ---*/
    
    bool Wrt_Vol = config->GetWrt_Vol_Sol();
    bool Wrt_Srf = config->GetWrt_Srf_Sol();
//    bool Wrt_Rst = config->GetWrt_Restart();
  
    /*--- Get the file output format ---*/
    
    unsigned short FileFormat = config->GetOutput_FileFormat();
    
    /*--- Merge the node coordinates and connectivity if necessary. This
     is only performed if a volume solution file is requested, and it
     is active by default. ---*/
    
    if (Wrt_Vol || Wrt_Srf) {
      cout <<"Merging grid connectivity." << endl;
      MergeConnectivity(config, geometry);
    }
    
    /*--- Merge the solution data needed for volume solutions and restarts ---*/
    
//    if (Wrt_Vol || Wrt_Rst) {
//      if (rank == MASTER_NODE) cout <<"Merging solution." << endl;
//      MergeBaselineSolution(config, geometry, solver, iZone);
//    }
    
    /*--- Write restart, CGNS, Tecplot or Paraview files using the merged data.
     This data lives only on the master, and these routines are currently
     executed by the master proc alone (as if in serial). ---*/
  
      if (Wrt_Vol) {
        
        cout <<"Writing volume solution file." << endl;
        
        switch (FileFormat) {
            
          case TECPLOT:
            
            /*--- Write a Tecplot ASCII file ---*/
            SetTecplot_ASCII(config, geometry, false);
            DeallocateConnectivity(config, geometry, false);
            break;
            
          case PARAVIEW:
            
            /*--- Write a Paraview ASCII file ---*/
            SetParaview_ASCII(config, geometry, false);
            DeallocateConnectivity(config, geometry, false);
            break;
            
          default:
            break;
        }
        
      }
      
      if (Wrt_Srf) {
        
        cout <<"Writing surface solution file." << endl;
        
        switch (FileFormat) {
            
          case TECPLOT:
            
            /*--- Write a Tecplot ASCII file ---*/
            SetTecplot_ASCII(config, geometry, true);
            DeallocateConnectivity(config, geometry, true);
            break;
            
          case PARAVIEW:
            
            /*--- Write a Paraview ASCII file ---*/
            SetParaview_ASCII(config, geometry, true);
            DeallocateConnectivity(config, geometry, true);
            break;
            
          default:
            break;
        }
      }
      
      if (FileFormat == TECPLOT_BINARY) {
        if (!wrote_base_file)
        DeallocateConnectivity(config, geometry, false);
        if (!wrote_surf_file)
        DeallocateConnectivity(config, geometry, wrote_surf_file);
      }
      
      if (Wrt_Vol || Wrt_Srf)
      DeallocateSolution(config, geometry);
  
}

