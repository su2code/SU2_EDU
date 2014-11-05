/*!
 * \file output_structure.cpp
 * \brief Main subroutines for output solver information.
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
                                 CSolver *FlowSolver, unsigned long iExtIter,
                                 unsigned short val_iZone) {
  
  unsigned short iMarker;
  unsigned long iPoint, iVertex, Global_Index;
  double PressCoeff = 0.0, SkinFrictionCoeff, HeatFlux;
  double xCoord = 0.0, yCoord = 0.0, zCoord = 0.0, Mach, Pressure;
  char cstr[200];
  
  unsigned short solver = config->GetKind_Solver();
  unsigned short nDim = geometry->GetnDim();
  
  char buffer [50];
  ofstream SurfFlow_file;
  
  
  /*--- Write file name with extension if unsteady ---*/
  strcpy (cstr, config->GetSurfFlowCoeff_FileName().c_str());
  
  if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) {
    if ((int(iExtIter) >= 0)    && (int(iExtIter) < 10))    sprintf (buffer, "_0000%d.csv", int(iExtIter));
    if ((int(iExtIter) >= 10)   && (int(iExtIter) < 100))   sprintf (buffer, "_000%d.csv",  int(iExtIter));
    if ((int(iExtIter) >= 100)  && (int(iExtIter) < 1000))  sprintf (buffer, "_00%d.csv",   int(iExtIter));
    if ((int(iExtIter) >= 1000) && (int(iExtIter) < 10000)) sprintf (buffer, "_0%d.csv",    int(iExtIter));
    if  (int(iExtIter) >= 10000) sprintf (buffer, "_%d.csv", int(iExtIter));
  }
  else
    sprintf (buffer, ".csv");
  
  strcat (cstr, buffer);
  SurfFlow_file.precision(15);
  SurfFlow_file.open(cstr, ios::out);
  
  SurfFlow_file << "\"Global_Index\", \"x_coord\", \"y_coord\", ";
  if (nDim == 3) SurfFlow_file << "\"z_coord\", ";
  SurfFlow_file << "\"Pressure\", \"Pressure_Coefficient\", ";
  
  switch (solver) {
    case EULER : SurfFlow_file <<  "\"Mach_Number\"" << endl; break;
    case NAVIER_STOKES: case RANS: SurfFlow_file <<  "\"Skin_Friction_Coefficient\", \"Heat_Flux\"" << endl; break;
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
            HeatFlux = FlowSolver->GetHeatFlux(iMarker,iVertex);
            SurfFlow_file << scientific << SkinFrictionCoeff << ", " << HeatFlux << endl;
            break;
        }
      }
    }
  }
  
  SurfFlow_file.close();
  
}

void COutput::SetSurfaceCSV_Adjoint(CConfig *config, CGeometry *geometry, CSolver *AdjSolver, CSolver *FlowSolution, unsigned long iExtIter, unsigned short val_iZone) {
  
  unsigned long iPoint, iVertex;
  double *Solution, xCoord, yCoord, zCoord, *IntBoundary_Jump;
  unsigned short iMarker;
  char cstr[200], buffer[50];
  ofstream SurfAdj_file;
  
  /*--- Write file name with extension if unsteady ---*/
  strcpy (cstr, config->GetSurfAdjCoeff_FileName().c_str());
  
  if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) {
    if ((int(iExtIter) >= 0)    && (int(iExtIter) < 10))    sprintf (buffer, "_0000%d.csv", int(iExtIter));
    if ((int(iExtIter) >= 10)   && (int(iExtIter) < 100))   sprintf (buffer, "_000%d.csv",  int(iExtIter));
    if ((int(iExtIter) >= 100)  && (int(iExtIter) < 1000))  sprintf (buffer, "_00%d.csv",   int(iExtIter));
    if ((int(iExtIter) >= 1000) && (int(iExtIter) < 10000)) sprintf (buffer, "_0%d.csv",    int(iExtIter));
    if  (int(iExtIter) >= 10000) sprintf (buffer, "_%d.csv", int(iExtIter));
  }
  else
    sprintf (buffer, ".csv");
  
  strcat(cstr, buffer);
  SurfAdj_file.precision(15);
  SurfAdj_file.open(cstr, ios::out);
  
  if (geometry->GetnDim() == 2) {
    SurfAdj_file <<  "\"Point\",\"Sensitivity\",\"PsiRho\",\"Phi_x\",\"Phi_y\",\"PsiE\",\"x_coord\",\"y_coord\"" << endl;
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_Plotting(iMarker) == YES)
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          Solution = AdjSolver->node[iPoint]->GetSolution();
          IntBoundary_Jump = AdjSolver->node[iPoint]->GetIntBoundary_Jump();
          xCoord = geometry->node[iPoint]->GetCoord(0);
          yCoord = geometry->node[iPoint]->GetCoord(1);
          SurfAdj_file << scientific << iPoint << ", " << AdjSolver->GetCSensitivity(iMarker,iVertex) << ", " << Solution[0] << ", "
          << Solution[1] << ", " << Solution[2] << ", " << Solution[3] <<", " << xCoord <<", "<< yCoord << endl;
        }
    }
  }
  
  if (geometry->GetnDim() == 3) {
    SurfAdj_file <<  "\"Point\",\"Sensitivity\",\"PsiRho\",\"Phi_x\",\"Phi_y\",\"Phi_z\",\"PsiE\",\"x_coord\",\"y_coord\",\"z_coord\"" << endl;
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_Plotting(iMarker) == YES)
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          Solution = AdjSolver->node[iPoint]->GetSolution();
          xCoord = geometry->node[iPoint]->GetCoord(0);
          yCoord = geometry->node[iPoint]->GetCoord(1);
          zCoord = geometry->node[iPoint]->GetCoord(2);
          
          SurfAdj_file << scientific << iPoint << ", " << AdjSolver->GetCSensitivity(iMarker,iVertex) << ", " << Solution[0] << ", "
          << Solution[1] << ", " << Solution[2] << ", " << Solution[3] << ", " << Solution[4] << ", "<< xCoord <<", "<< yCoord <<", "<< zCoord << endl;
        }
    }
  }
  
  SurfAdj_file.close();
  
}

void COutput::SetSurfaceCSV_Linearized(CConfig *config, CGeometry *geometry, CSolver *LinSolution, string val_filename, unsigned long iExtIter) { }

void COutput::MergeConnectivity(CConfig *config, CGeometry *geometry, unsigned short val_iZone) {
  
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
  
  int rank = MASTER_NODE;
  
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
  
  if (rank == MASTER_NODE) {
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
  
}

void COutput::MergeSurfaceConnectivity(CConfig *config, CGeometry *geometry, unsigned short Elem_Type) {
  
  int rank = MASTER_NODE;
  
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
  
  if (rank == MASTER_NODE) {
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
  
}

void COutput::MergeSolution(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned short val_iZone) {
  
  /*--- Local variables needed on all processors ---*/
  unsigned short Kind_Solver  = config->GetKind_Solver();
  unsigned short iVar = 0, jVar = 0, FirstIndex = NONE, SecondIndex = NONE, ThirdIndex = NONE;
  unsigned short nVar_First = 0, nVar_Second = 0, nVar_Third = 0, iVar_Eddy = 0;
  unsigned short iVar_GridVel = 0, iVar_PressMach = 0, iVar_TempLam = 0,
  iVar_ViscCoeffs = 0;
  
  unsigned long iPoint = 0, jPoint = 0, iVertex = 0, iMarker = 0;
  double Gas_Constant, Mach2Vel, Mach_Motion, RefDensity, RefPressure = 0.0, factor = 0.0;
  
  double *Aux_Frict = NULL, *Aux_Heat = NULL, *Aux_yPlus = NULL;
  
  bool grid_movement  = (config->GetGrid_Movement());
  bool flow           = (( config->GetKind_Solver() == EULER             ) ||
                         ( config->GetKind_Solver() == NAVIER_STOKES     ) ||
                         ( config->GetKind_Solver() == RANS              )   );
  
  unsigned short iDim;
  unsigned short nDim = geometry->GetnDim();
  double RefAreaCoeff = config->GetRefAreaCoeff();
  double Gamma = config->GetGamma();
  double RefVel2;
  
  /*--- Set the non-dimensionalization ---*/
  if (flow) {
    if (grid_movement) {
      Gas_Constant = config->GetGas_ConstantND();
      Mach2Vel = sqrt(Gamma*Gas_Constant*config->GetTemperature_FreeStreamND());
      Mach_Motion = config->GetMach_Motion();
      RefVel2 = (Mach_Motion*Mach2Vel)*(Mach_Motion*Mach2Vel);
    }
    else {
      RefVel2 = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        RefVel2  += solver[FLOW_SOL]->GetVelocity_Inf(iDim)*solver[FLOW_SOL]->GetVelocity_Inf(iDim);
    }
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
  
  /*--- Add the grid velocity to the restart file for the unsteady adjoint ---*/
  if (grid_movement) {
    iVar_GridVel = nVar_Total;
    if (geometry->GetnDim() == 2) nVar_Total += 2;
    else if (geometry->GetnDim() == 3) nVar_Total += 3;
  }
  
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
  
  /*--- Merge the solution either in serial or parallel. ---*/
  
  /*--- In serial, the single process has access to all solution data,
   so it is simple to retrieve and store inside Solution_Data. ---*/
  
  nGlobal_Poin = geometry->GetnPointDomain();
  Data = new double*[nVar_Total];
  for (iVar = 0; iVar < nVar_Total; iVar++) {
    Data[iVar] = new double[nGlobal_Poin];
  }
  
  /*--- In case there is grid movement ---*/
  double *Grid_Vel;
  
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
          Aux_Heat[iPoint]  = solver[FLOW_SOL]->GetHeatFlux(iMarker,iVertex);
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
      
      /*--- Residual (first, second and third system of equations) ---*/
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
      
      /*--- For unsteady problems with grid movement, write the mesh velocities ---*/
      if (grid_movement) {
        Grid_Vel = geometry->node[iPoint]->GetGridVel();
        for (unsigned short iDim = 0; iDim < geometry->GetnDim(); iDim++) {
          Data[jVar][jPoint] = Grid_Vel[iDim];
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
          break;
          
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

void COutput::MergeBaselineSolution(CConfig *config, CGeometry *geometry, CSolver *solver, unsigned short val_iZone) {
  
  /*--- Local variables needed on all processors ---*/
  unsigned short iVar;
  unsigned long iPoint = 0, jPoint = 0;
  
  nVar_Total = config->fields.size() - 1;
  
  /*--- Merge the solution either in serial or parallel. ---*/
  
  /*--- In serial, the single process has access to all solution data,
   so it is simple to retrieve and store inside Solution_Data. ---*/
  nGlobal_Poin = geometry->GetnPointDomain();
  Data = new double*[nVar_Total];
  for (iVar = 0; iVar < nVar_Total; iVar++) {
    Data[iVar] = new double[nGlobal_Poin];
  }
  
  /*--- Loop over all points in the mesh, but only write data
   for nodes in the domain (ignore periodic halo nodes). ---*/
  jPoint = 0;
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    if (geometry->node[iPoint]->GetDomain()) {
      /*--- Solution (first, and second system of equations) ---*/
      unsigned short jVar = 0;
      for (iVar = 0; iVar < nVar_Total; iVar++) {
        Data[jVar][jPoint] = solver->node[iPoint]->GetSolution(iVar);
        jVar++;
      }
    }
    
    /*--- Increment jPoint as the counter. We need this because iPoint
     may include halo nodes that we skip over during this loop. ---*/
    jPoint++;
    
  }
  
}

void COutput::SetRestart(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned short val_iZone) {
  
  /*--- Local variables ---*/
  unsigned short Kind_Solver  = config->GetKind_Solver();
  unsigned short iVar, iDim, nDim = geometry->GetnDim();
  unsigned long iPoint, iExtIter = config->GetExtIter();
  bool grid_movement = config->GetGrid_Movement();
  ofstream restart_file;
  string filename;
  
  /*--- Retrieve filename from config ---*/
  if (config->GetAdjoint()) {
    filename = config->GetRestart_AdjFileName();
    filename = config->GetObjFunc_Extension(filename);
  } else {
    filename = config->GetRestart_FlowFileName();
  }
  
  /*--- Unsteady problems require an iteration number to be appended. ---*/
  if (config->GetWrt_Unsteady()) {
    filename = config->GetUnsteady_FileName(filename, int(iExtIter));
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
  
  /*--- Mesh velocities for dynamic mesh cases ---*/
  if (grid_movement) {
    if (nDim == 2) {
      restart_file << "\t\"Grid_Velx\"\t\"Grid_Vely\"";
    } else {
      restart_file << "\t\"Grid_Velx\"\t\"Grid_Vely\"\t\"Grid_Velz\"";
    }
  }
  
  if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
    restart_file << "\t\"Pressure\"\t\"Pressure_Coefficient\"\t\"Mach\"";
  }
  
  if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
    restart_file << "\t\"Temperature\"\t\"Laminar_Viscosity\"\t\"Skin_Friction_Coefficient\"\t\"Heat_Flux\"\t\"Y_Plus\"";
  }
  
  if (Kind_Solver == RANS) {
    restart_file << "\t\"Eddy_Viscosity\"";
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
    restart_file << endl;
  }
  
  restart_file.close();
  
}

void COutput::DeallocateCoordinates(CConfig *config, CGeometry *geometry) {
  
  int rank = MASTER_NODE;
  
  /*--- Local variables and initialization ---*/
  
  unsigned short iDim, nDim = geometry->GetnDim();
  
  /*--- The master node alone owns all data found in this routine. ---*/
  if (rank == MASTER_NODE) {
    
    /*--- Deallocate memory for coordinate data ---*/
    for (iDim = 0; iDim < nDim; iDim++) {
      delete [] Coords[iDim];
    }
    delete [] Coords;
    
  }
}

void COutput::DeallocateConnectivity(CConfig *config, CGeometry *geometry, bool surf_sol) {
  
  int rank = MASTER_NODE;
  
  /*--- The master node alone owns all data found in this routine. ---*/
  if (rank == MASTER_NODE) {
    
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
}

void COutput::DeallocateSolution(CConfig *config, CGeometry *geometry) {
  
  int rank = MASTER_NODE;
  
  /*--- The master node alone owns all data found in this routine. ---*/
  if (rank == MASTER_NODE) {
    
    /*--- Deallocate memory for solution data ---*/
    for (unsigned short iVar = 0; iVar < nVar_Total; iVar++) {
      delete [] Data[iVar];
    }
    delete [] Data;
    
  }
}

void COutput::SetHistory_Header(ofstream *ConvHist_file, CConfig *config) {
  char cstr[200], buffer[50], turb_resid[1000];
  unsigned short iMarker, iMarker_Monitoring;
  string Monitoring_Tag, monitoring_coeff, aeroelastic_coeff;
  
  bool equiv_area = config->GetEquivArea();
  bool turbulent = ((config->GetKind_Solver() == RANS) || (config->GetKind_Solver() == ADJ_RANS));
  bool inv_design = (config->GetInvDesign_Cp() || config->GetInvDesign_HeatFlux());
  bool output_1d = config->GetWrt_1D_Output();
  bool output_per_surface = false;
  if(config->GetnMarker_Monitoring() > 1) output_per_surface = true;
  
  bool isothermal = false;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if ((config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL             ) ||
        (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL_CATALYTIC   ) ||
        (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL_NONCATALYTIC)   )
      isothermal = true;
  
  /*--- Write file name with extension ---*/
  string filename = config->GetConv_FileName();
  strcpy (cstr, filename.data());
  
  if (config->GetWrt_Unsteady() && config->GetRestart()) {
    long iExtIter = config->GetUnst_RestartIter();
    if (int(iExtIter) < 10) sprintf (buffer, "_0000%d", int(iExtIter));
    if ((int(iExtIter) >= 10) && (int(iExtIter) < 100)) sprintf (buffer, "_000%d", int(iExtIter));
    if ((int(iExtIter) >= 100) && (int(iExtIter) < 1000)) sprintf (buffer, "_00%d", int(iExtIter));
    if ((int(iExtIter) >= 1000) && (int(iExtIter) < 10000)) sprintf (buffer, "_0%d", int(iExtIter));
    if (int(iExtIter) >= 10000) sprintf (buffer, "_%d", int(iExtIter));
    strcat(cstr,buffer);
  }
  
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
  char heat_coeff[]= ",\"HeatFlux_Total\",\"HeatFlux_Maximum\"";
  char equivalent_area_coeff[]= ",\"CEquivArea\",\"CNearFieldOF\"";
  char oneD_outputs[]= ",\"1D_Pt\",\"1D_Mach\",\"1D_Temperature\",\"Flux Avg Pressure\",\"Flux Avg Density\",\"Flux Avg Velocity\",\"Flux Avg Enthalpy\"";
  char Cp_inverse_design[]= ",\"Cp_Diff\"";
  char Heat_inverse_design[]= ",\"HeatFlux_Diff\"";
  
  /* Find the markers being monitored and create a header for them */
  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    Monitoring_Tag = config->GetMarker_Monitoring(iMarker_Monitoring);
    monitoring_coeff += ",\"CLift_"  + Monitoring_Tag + "\"";
    monitoring_coeff += ",\"CDrag_"  + Monitoring_Tag + "\"";
    monitoring_coeff += ",\"CMx_"    + Monitoring_Tag + "\"";
    monitoring_coeff += ",\"CMy_"    + Monitoring_Tag + "\"";
    monitoring_coeff += ",\"CMz_"    + Monitoring_Tag + "\"";
  }
  
  /*--- Header for the residuals ---*/
  
  char flow_resid[]= ",\"Res_Flow[0]\",\"Res_Flow[1]\",\"Res_Flow[2]\",\"Res_Flow[3]\",\"Res_Flow[4]\"";
  switch (config->GetKind_Turb_Model()) {
    case SA:	sprintf (turb_resid, ",\"Res_Turb[0]\""); break;
    case ML:	sprintf (turb_resid, ",\"Res_Turb[0]\""); break;
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
      if (equiv_area) ConvHist_file[0] << equivalent_area_coeff;
      if (inv_design) {
        ConvHist_file[0] << Cp_inverse_design;
        if (isothermal) ConvHist_file[0] << Heat_inverse_design;
      }
      ConvHist_file[0] << flow_resid;
      if (turbulent) ConvHist_file[0] << turb_resid;
      if (output_per_surface) ConvHist_file[0] << monitoring_coeff;
      if (output_1d) ConvHist_file[0] << oneD_outputs;
      ConvHist_file[0] << end;
      break;
      
  }
  
  if (config->GetOutput_FileFormat() == TECPLOT || config->GetOutput_FileFormat() == TECPLOT_BINARY) {
    ConvHist_file[0] << "ZONE T= \"Convergence history\"" << endl;
  }
  
}


void COutput::SetConvergence_History(ofstream *ConvHist_file,
                                     CGeometry **geometry,
                                     CSolver ***solver_container,
                                     CConfig *config,
                                     CIntegration **integration,
                                     bool DualTime_Iteration,
                                     double timeused) {
  
	int rank;
  rank = MASTER_NODE;
  
  /*--- Output using only the master node ---*/
  if (rank == MASTER_NODE) {
    
    unsigned long iIntIter = config->GetIntIter();
    unsigned long iExtIter = config->GetExtIter();
    
    /*--- WARNING: These buffers have hard-coded lengths. Note that you
     may have to adjust them to be larger if adding more entries. ---*/
    char begin[1000], direct_coeff[1000], surface_coeff[1000], monitoring_coeff[1000],
    adjoint_coeff[1000], flow_resid[1000], adj_flow_resid[1000], turb_resid[1000], trans_resid[1000],
    adj_turb_resid[1000], levelset_resid[1000], end[1000];
    char oneD_outputs[1000];
    double dummy = 0.0;
    unsigned short iVar, iMarker, iMarker_Monitoring;
    
    unsigned long LinSolvIter = 0;
    double timeiter = timeused/double(iExtIter+1);
    
    unsigned short FinestMesh = config->GetFinestMesh();
    unsigned short nDim = geometry[FinestMesh]->GetnDim();
    
    bool equiv_area = config->GetEquivArea();
    bool inv_design = (config->GetInvDesign_Cp() || config->GetInvDesign_HeatFlux());
    bool transition = (config->GetKind_Trans_Model() == LM);
    bool isothermal = false;
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
      if ((config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL) ||
          (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL_CATALYTIC) ||
          (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL_NONCATALYTIC))
        isothermal = true;
    bool turbulent = ((config->GetKind_Solver() == RANS) || (config->GetKind_Solver() == ADJ_RANS) ||
                      (config->GetKind_Solver() == FLUID_STRUCTURE_RANS));
    bool adjoint = config->GetAdjoint();
    bool fluid_structure = ((config->GetKind_Solver() == FLUID_STRUCTURE_EULER) || (config->GetKind_Solver() == FLUID_STRUCTURE_NAVIER_STOKES) ||
                            (config->GetKind_Solver() == FLUID_STRUCTURE_RANS));
    bool flow = (config->GetKind_Regime() == EULER) || (config->GetKind_Regime() == NAVIER_STOKES) ||
    (config->GetKind_Regime() == RANS) || (config->GetKind_Regime() == ADJ_EULER) ||
    (config->GetKind_Regime() == ADJ_NAVIER_STOKES) || (config->GetKind_Regime() == ADJ_RANS);
    
    bool output_1d  = (config->GetWrt_1D_Output());
    bool output_per_surface = false;
    if(config->GetnMarker_Monitoring() > 1) output_per_surface = true;
    
    /*--- Initialize variables to store information from all domains (direct solution) ---*/
    double Total_CLift = 0.0, Total_CDrag = 0.0, Total_CSideForce = 0.0, Total_CMx = 0.0, Total_CMy = 0.0, Total_CMz = 0.0, Total_CEff = 0.0,
    Total_CEquivArea = 0.0, Total_CNearFieldOF = 0.0, Total_CFx = 0.0, Total_CFy = 0.0, Total_CFz = 0.0,
    Total_CpDiff = 0.0, Total_HeatFluxDiff = 0.0,
    Total_CFEA = 0.0, Total_Heat = 0.0, Total_MaxHeat = 0.0;
    
    double OneD_Stagnation_Pressure = 0.0, OneD_Mach = 0.0, OneD_Temp = 0.0,
        OneD_FluxAvgP=0.0, OneD_FluxAvgRho=0.0, OneD_FluxAvgU=0.0, OneD_FluxAvgH=0.0;
    
    /*--- Initialize variables to store information from all domains (adjoint solution) ---*/
    double Total_Sens_Geo = 0.0, Total_Sens_Mach = 0.0, Total_Sens_AoA = 0.0;
    double Total_Sens_Press = 0.0, Total_Sens_Temp = 0.0;
    
    /*--- Residual arrays ---*/
    double *residual_flow         = NULL,
           *residual_turbulent    = NULL,
           *residual_transition   = NULL,
           *residual_TNE2         = NULL,
           *residual_levelset     = NULL;
    double *residual_adjflow      = NULL,
           *residual_adjturbulent = NULL,
           *residual_adjTNE2      = NULL,
           *residual_adjlevelset  = NULL;
    double *residual_wave         = NULL;
    double *residual_fea          = NULL;
    double *residual_heat         = NULL;
    
    /*--- Coefficients Monitored arrays ---*/
    double *Surface_CLift      = NULL,
           *Surface_CDrag      = NULL,
           *Surface_CMx        = NULL,
           *Surface_CMy        = NULL,
           *Surface_CMz        = NULL;
    
    /*--- Initialize number of variables ---*/
    unsigned short nVar_Flow = 0, nVar_LevelSet = 0, nVar_Turb = 0,
    nVar_Trans = 0, nVar_TNE2 = 0, nVar_Wave = 0, nVar_Heat = 0, nVar_FEA = 0,
    nVar_AdjFlow = 0, nVar_AdjTNE2 = 0, nVar_AdjLevelSet = 0, nVar_AdjTurb = 0;
    
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
    Surface_CLift      = new double[config->GetnMarker_Monitoring()];
    Surface_CDrag      = new double[config->GetnMarker_Monitoring()];
    Surface_CMx        = new double[config->GetnMarker_Monitoring()];
    Surface_CMy        = new double[config->GetnMarker_Monitoring()];
    Surface_CMz        = new double[config->GetnMarker_Monitoring()];
    
    /*--- Write information from nodes ---*/
    switch (config->GetKind_Solver()) {
        
      case EULER:                   case NAVIER_STOKES:                   case RANS:
      case FLUID_STRUCTURE_EULER:   case FLUID_STRUCTURE_NAVIER_STOKES:   case FLUID_STRUCTURE_RANS:
      case ADJ_EULER:               case ADJ_NAVIER_STOKES:               case ADJ_RANS:
        
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
        
        if (isothermal) {
          Total_Heat     = solver_container[FinestMesh][FLOW_SOL]->GetTotal_HeatFlux();
          Total_MaxHeat  = solver_container[FinestMesh][FLOW_SOL]->GetTotal_MaxHeatFlux();
        }
        
        if (equiv_area) {
          Total_CEquivArea    = solver_container[FinestMesh][FLOW_SOL]->GetTotal_CEquivArea();
          Total_CNearFieldOF  = solver_container[FinestMesh][FLOW_SOL]->GetTotal_CNearFieldOF();
          
          /*--- Note that there is a redefinition of the nearfield based functionals ---*/
          Total_CEquivArea    = config->GetWeightCd()*Total_CDrag + (1.0-config->GetWeightCd())*Total_CEquivArea;
          Total_CNearFieldOF  = config->GetWeightCd()*Total_CDrag + (1.0-config->GetWeightCd())*Total_CNearFieldOF;
        }
        
        if (inv_design) {
          Total_CpDiff  = solver_container[FinestMesh][FLOW_SOL]->GetTotal_CpDiff();
          if (isothermal) {
            Total_HeatFluxDiff = solver_container[FinestMesh][FLOW_SOL]->GetTotal_HeatFluxDiff();
          }
        }
        
        if (output_per_surface) {
          /*--- Look over the markers being monitored and get the desired values ---*/
          for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
            Surface_CLift[iMarker_Monitoring]      = solver_container[FinestMesh][FLOW_SOL]->GetSurface_CLift(iMarker_Monitoring);
            Surface_CDrag[iMarker_Monitoring]      = solver_container[FinestMesh][FLOW_SOL]->GetSurface_CDrag(iMarker_Monitoring);
            Surface_CMx[iMarker_Monitoring]        = solver_container[FinestMesh][FLOW_SOL]->GetSurface_CMx(iMarker_Monitoring);
            Surface_CMy[iMarker_Monitoring]        = solver_container[FinestMesh][FLOW_SOL]->GetSurface_CMy(iMarker_Monitoring);
            Surface_CMz[iMarker_Monitoring]        = solver_container[FinestMesh][FLOW_SOL]->GetSurface_CMz(iMarker_Monitoring);
          }
        }
        
        if (fluid_structure) {
          Total_CFEA  = solver_container[FinestMesh][FEA_SOL]->GetTotal_CFEA();
        }
        
        if (output_1d) {
          /* Get area-averaged and flux-averaged values at the specified surface*/
          OneD_Stagnation_Pressure=solver_container[FinestMesh][FLOW_SOL]->GetOneD_Pt();
          OneD_Mach=solver_container[FinestMesh][FLOW_SOL]->GetOneD_M();
          OneD_Temp=solver_container[FinestMesh][FLOW_SOL]->GetOneD_T();
          OneD_FluxAvgP=solver_container[FinestMesh][FLOW_SOL]->GetOneD_fluxavgP();
          OneD_FluxAvgRho=solver_container[FinestMesh][FLOW_SOL]->GetOneD_fluxavgRho();
          OneD_FluxAvgU=solver_container[FinestMesh][FLOW_SOL]->GetOneD_fluxavgU();
          OneD_FluxAvgH=solver_container[FinestMesh][FLOW_SOL]->GetOneD_fluxavgH();
        }
        
        /*--- Flow Residuals ---*/
        
        for (iVar = 0; iVar < nVar_Flow; iVar++)
          residual_flow[iVar] = solver_container[FinestMesh][FLOW_SOL]->GetRes_RMS(iVar);
        
        /*--- Turbulent residual ---*/
        
        if (turbulent) {
          for (iVar = 0; iVar < nVar_Turb; iVar++)
            residual_turbulent[iVar] = solver_container[FinestMesh][TURB_SOL]->GetRes_RMS(iVar);
        }
        
        /*--- Transition residual ---*/
        
        if (transition) {
          for (iVar = 0; iVar < nVar_Trans; iVar++)
            residual_transition[iVar] = solver_container[FinestMesh][TRANS_SOL]->GetRes_RMS(iVar);
        }
        
        /*--- FEA residual ---*/
        
        if (fluid_structure) {
          for (iVar = 0; iVar < nVar_FEA; iVar++)
            residual_fea[iVar] = solver_container[FinestMesh][FEA_SOL]->GetRes_RMS(iVar);
        }
        
        /*--- Iterations of the linear solver ---*/
        
        LinSolvIter = (unsigned long) solver_container[FinestMesh][FLOW_SOL]->GetIterLinSolver();
        
        break;
        
    }
    
    /*--- Header frecuency ---*/
    
    bool Unsteady = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                     (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
    bool In_NoDualTime = (!DualTime_Iteration && (iExtIter % config->GetWrt_Con_Freq() == 0));
    bool In_DualTime_0 = (DualTime_Iteration && (iIntIter % config->GetWrt_Con_Freq_DualTime() == 0));
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
          case FLUID_STRUCTURE_EULER: case FLUID_STRUCTURE_NAVIER_STOKES: case FLUID_STRUCTURE_RANS:
          case ADJ_EULER: case ADJ_NAVIER_STOKES: case ADJ_RANS:
            
            /*--- Direct coefficients ---*/
            sprintf (direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f",
                     Total_CLift, Total_CDrag, Total_CSideForce, Total_CMx, Total_CMy, Total_CMz, Total_CFx, Total_CFy,
                     Total_CFz, Total_CEff);
            if (isothermal)
              sprintf (direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", Total_CLift, Total_CDrag, Total_CSideForce, Total_CMx, Total_CMy,
                       Total_CMz, Total_CFx, Total_CFy, Total_CFz, Total_CEff, Total_Heat, Total_MaxHeat);
            if (equiv_area)
              sprintf (direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", Total_CLift, Total_CDrag, Total_CSideForce, Total_CMx, Total_CMy, Total_CMz, Total_CFx, Total_CFy, Total_CFz, Total_CEff, Total_CEquivArea, Total_CNearFieldOF);
            if (inv_design) {
              sprintf (direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", Total_CLift, Total_CDrag, Total_CSideForce, Total_CMx, Total_CMy, Total_CMz, Total_CFx, Total_CFy, Total_CFz, Total_CEff, Total_CpDiff);
              Total_CpDiff  = solver_container[FinestMesh][FLOW_SOL]->GetTotal_CpDiff();
              if (isothermal) {
                sprintf (direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", Total_CLift, Total_CDrag, Total_CSideForce, Total_CMx, Total_CMy, Total_CMz, Total_CFx, Total_CFy, Total_CFz, Total_CEff, Total_Heat, Total_MaxHeat, Total_CpDiff, Total_HeatFluxDiff);
              }
            }
            if (fluid_structure)
              sprintf (direct_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", Total_CLift, Total_CDrag, Total_CSideForce, Total_CMx, Total_CMy, Total_CMz,
                       Total_CFx, Total_CFy, Total_CFz, Total_CEff, Total_CFEA);

            if (output_per_surface) {
              for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
                //Append one by one the surface coeff to monitoring coeff. (Think better way do this, maybe use string)
                if (iMarker_Monitoring == 0) {
                  sprintf(monitoring_coeff, ", %12.10f",Surface_CLift[iMarker_Monitoring]);
                }
                else {
                  sprintf(surface_coeff, ", %12.10f",Surface_CLift[iMarker_Monitoring]);
                  strcat(monitoring_coeff, surface_coeff);
                }
                sprintf(surface_coeff, ", %12.10f",Surface_CDrag[iMarker_Monitoring]);
                strcat(monitoring_coeff, surface_coeff);
                sprintf(surface_coeff, ", %12.10f",Surface_CMx[iMarker_Monitoring]);
                strcat(monitoring_coeff, surface_coeff);
                sprintf(surface_coeff, ", %12.10f",Surface_CMy[iMarker_Monitoring]);
                strcat(monitoring_coeff, surface_coeff);
                sprintf(surface_coeff, ", %12.10f",Surface_CMz[iMarker_Monitoring]);
                strcat(monitoring_coeff, surface_coeff);
              }
            }

            
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
            /*---- Averaged stagnation pressure at an exit ---- */
            if (output_1d)
              sprintf( oneD_outputs, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", OneD_Stagnation_Pressure, OneD_Mach, OneD_Temp,OneD_FluxAvgP, OneD_FluxAvgRho, OneD_FluxAvgU, OneD_FluxAvgH);

            /*--- Transition residual ---*/
            if (transition){
              sprintf (trans_resid, ", %12.10f, %12.10f", log10(residual_transition[0]), log10(residual_transition[1]));
            }
            
            /*--- Fluid structure residual ---*/
            if (fluid_structure) {
              if (nDim == 2) sprintf (levelset_resid, ", %12.10f, %12.10f, 0.0", log10 (residual_fea[0]), log10 (residual_fea[1]));
              else sprintf (levelset_resid, ", %12.10f, %12.10f, %12.10f", log10 (residual_fea[0]), log10 (residual_fea[1]), log10 (residual_fea[2]));
            }
            
            if (adjoint) {
              
              /*--- Adjoint coefficients ---*/
              sprintf (adjoint_coeff, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f, 0.0", Total_Sens_Geo, Total_Sens_Mach, Total_Sens_AoA, Total_Sens_Press, Total_Sens_Temp);
              
              /*--- Adjoint flow residuals ---*/
              if (nDim == 2) {
                sprintf (adj_flow_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, 0.0", log10 (residual_adjflow[0]),log10 (residual_adjflow[1]),log10 (residual_adjflow[2]),log10 (residual_adjflow[3]) );
              }
              else {
                sprintf (adj_flow_resid, ", %12.10f, %12.10f, %12.10f, %12.10f, %12.10f", log10 (residual_adjflow[0]),log10 (residual_adjflow[1]),log10 (residual_adjflow[2]),log10 (residual_adjflow[3]), log10 (residual_adjflow[4]) );
              }
              
              /*--- Adjoint turbulent residuals ---*/
              if (turbulent)
                if (!config->GetFrozen_Visc())
                  sprintf (adj_turb_resid, ", %12.10f", log10 (residual_adjturbulent[0]));
              
            }
            
            break;
            
        }
      }
      
      /*--- Write the screen header---*/
      if ((write_heads) && !(!DualTime_Iteration && Unsteady)) {
        
        if (!Unsteady) {
          switch (config->GetKind_Solver()) {
            case EULER :  case NAVIER_STOKES: case RANS:
              cout << endl << " Local time stepping summary:" << endl;
              for (unsigned short iMesh = 0; iMesh <= config->GetMGLevels(); iMesh++)
                cout << " MG level: "<< iMesh << "-> Min. DT: " << solver_container[iMesh][FLOW_SOL]->GetMin_Delta_Time()<<
                ". Max. DT: " << solver_container[iMesh][FLOW_SOL]->GetMax_Delta_Time() << "." << endl;
              
              break;
           }
        }
        else {
          if (flow) {
            cout << endl << " Min DT: " << solver_container[FinestMesh][FLOW_SOL]->GetMin_Delta_Time()<<
            ". Max DT: " << solver_container[FinestMesh][FLOW_SOL]->GetMax_Delta_Time() <<
            ". Dual Time step: " << config->GetDelta_UnstTimeND() << ".";
          }
          else {
            cout << endl << " Dual Time step: " << config->GetDelta_UnstTimeND() << ".";
          }
        }
        
        switch (config->GetKind_Solver()) {
          case EULER : case NAVIER_STOKES:
            
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
              case ML:	cout << "       Res[nu]"; break;
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
        case FLUID_STRUCTURE_EULER: case FLUID_STRUCTURE_NAVIER_STOKES:
          
          if (!DualTime_Iteration) {
            ConvHist_file[0] << begin << direct_coeff << flow_resid;
            /*---- Averaged stagnation pressure at an exit ---- */
            ConvHist_file[0] << end;
            ConvHist_file[0].flush();
          }
          
          cout.precision(6);
          cout.setf(ios::fixed,ios::floatfield);
          cout.width(13); cout << log10(residual_flow[0]);
          if (!fluid_structure && !equiv_area) {
              if (nDim == 2 ) { cout.width(14); cout << log10(residual_flow[3]); }
              else { cout.width(14); cout << log10(residual_flow[4]); }
          }
          else if (fluid_structure) { cout.width(14); cout << log10(residual_fea[0]); }
          
          cout.width(15); cout << min(1000.0,max(-1000.0, Total_CLift)); cout.width(15); cout << min(1000.0,max(-1000.0, Total_CDrag));
          cout << endl;
          
          break;
          
        case RANS :
          
          if (!DualTime_Iteration) {
            ConvHist_file[0] << begin << direct_coeff << flow_resid << turb_resid;
            if (output_per_surface) ConvHist_file[0] << monitoring_coeff;
            ConvHist_file[0] << end;
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
          
          if (transition) { cout.width(14); cout << log10(residual_transition[0]); cout.width(14); cout << log10(residual_transition[1]); }
          
          cout.width(15); cout << min(1000.0,max(-1000.0, Total_CLift)); cout.width(15); cout << min(1000.0,max(-1000.0, Total_CDrag));
          cout << endl;
          
          break;
          
      }
      cout.unsetf(ios::fixed);
      
      delete [] residual_flow;
      delete [] residual_turbulent;
      delete [] residual_transition;
      delete [] residual_TNE2;
      delete [] residual_levelset;
      delete [] residual_wave;
      delete [] residual_fea;
      delete [] residual_heat;
      
      delete [] residual_adjflow;
      delete [] residual_adjturbulent;
      delete [] residual_adjTNE2;
      delete [] residual_adjlevelset;
      
      delete [] Surface_CLift;
      delete [] Surface_CDrag;
      delete [] Surface_CMx;
      delete [] Surface_CMy;
      delete [] Surface_CMz;
    }
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
      
      if (Wrt_Csv) SetSurfaceCSV_Flow(config, geometry[MESH_0], solver_container[MESH_0][FLOW_SOL], iExtIter, ZONE_0);
      break;
      
  }
  
  /*--- Get the file output format ---*/
  
  unsigned short FileFormat = config->GetOutput_FileFormat();
  
  /*--- Merge the node coordinates and connectivity, if necessary. This
   is only performed if a volume solution file is requested, and it
   is active by default. ---*/
  
  if (Wrt_Vol || Wrt_Srf)
    MergeConnectivity(config, geometry[MESH_0], ZONE_0);
  
  /*--- Merge coordinates of all grid nodes (excluding ghost points).
   The grid coordinates are always merged and included first in the
   restart files. ---*/
  
  MergeCoordinates(config, geometry[MESH_0]);
  
  /*--- Merge the solution data needed for volume solutions and restarts ---*/
  
  if (Wrt_Vol || Wrt_Rst)
    MergeSolution(config, geometry[MESH_0], solver_container[MESH_0], ZONE_0);
  
  /*--- Write restart, CGNS, or Tecplot files using the merged data.
   This data lives only on the master, and these routines are currently
   executed by the master proc alone (as if in serial). ---*/
  
  /*--- Write a native restart file ---*/
  if (Wrt_Rst)
    SetRestart(config, geometry[MESH_0], solver_container[MESH_0], ZONE_0);
  
  if (Wrt_Vol) {
    
    switch (FileFormat) {
        
      case TECPLOT:
        
        /*--- Write a Tecplot ASCII file ---*/
        SetTecplot_ASCII(config, geometry[MESH_0], solver_container[MESH_0], ZONE_0, 1, false);
        DeallocateConnectivity(config, geometry[MESH_0], false);
        break;
        
      case PARAVIEW:
        
        /*--- Write a Paraview ASCII file ---*/
        SetParaview_ASCII(config, geometry[MESH_0], ZONE_0, 1, false);
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
        SetTecplot_ASCII(config, geometry[MESH_0], solver_container[MESH_0], ZONE_0, 1, true);
        DeallocateConnectivity(config, geometry[MESH_0], true);
        break;
        
      case PARAVIEW:
        
        /*--- Write a Paraview ASCII file ---*/
        SetParaview_ASCII(config, geometry[MESH_0], ZONE_0, 1, true);
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
