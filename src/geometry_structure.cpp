/*!
 * \file geometry_structure.cpp
 * \brief Main subroutines for creating the primal grid and multigrid structure.
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

#include "../include/geometry_structure.hpp"

CGeometry::CGeometry(void) {
  
  nEdge = 0;
  nPoint = 0;
  nElem = 0;
  
  nElem_Bound = NULL;
  Tag_to_Marker = NULL;
  elem = NULL;
  face = NULL;
  bound = NULL;
  node = NULL;
  edge = NULL;
  vertex = NULL;
  nVertex = NULL;
  newBound = NULL;
  nNewElem_Bound = NULL;
  Marker_All_SendRecv = NULL;
  
  //	PeriodicPoint[MAX_NUMBER_PERIODIC][2].clear();
  //	PeriodicElem[MAX_NUMBER_PERIODIC].clear();
  //	OldBoundaryElems[MAX_NUMBER_MARKER].clear();
  //  SendTransfLocal[MAX_NUMBER_DOMAIN].clear();
  //  ReceivedTransfLocal[MAX_NUMBER_DOMAIN].clear();
  //	SendDomainLocal[MAX_NUMBER_DOMAIN].clear();
  //	ReceivedDomainLocal[MAX_NUMBER_DOMAIN].clear();
  //	XCoordList.clear();
  
  //	Xcoord_plane.clear();
  //	Ycoord_plane.clear();
  //	Zcoord_plane.clear();
  //	FaceArea_plane.clear();
  //	Plane_points.clear();
  
}

CGeometry::~CGeometry(void) {
  unsigned long iElem, iElem_Bound, iPoint, iFace, iVertex, iEdge;
  unsigned short iMarker;
  
  if (elem != NULL) {
    for (iElem = 0; iElem < nElem; iElem++)
      if (elem[iElem] != NULL) delete elem[iElem];
    delete[] elem;
  }
  
  if (bound != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
        if (bound[iMarker][iElem_Bound] != NULL) delete bound[iMarker][iElem_Bound];
      }
    }
    delete[] bound;
  }
  
  if (face != NULL) {
    for (iFace = 0; iFace < nFace; iFace ++)
      if (face[iFace] != NULL) delete face[iFace];
    delete[] face;
  }
  
  if (node != NULL) {
    for (iPoint = 0; iPoint < nPoint; iPoint ++)
      if (node[iPoint] != NULL) delete node[iPoint];
    delete[] node;
  }
  
  if (edge != NULL) {
    for (iEdge = 0; iEdge < nEdge; iEdge ++)
      if (edge[iPoint] != NULL) delete edge[iEdge];
    delete[] edge;
  }
  
  if (vertex != NULL)  {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
        if (vertex[iMarker][iVertex] != NULL) delete vertex[iMarker][iVertex];
      }
    }
    delete[] vertex;
  }
  
  if (newBound != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
        if (newBound[iMarker][iElem_Bound] != NULL) delete newBound[iMarker][iElem_Bound];
      }
    }
    delete[] newBound;
  }
  
  if (nElem_Bound != NULL) delete[] nElem_Bound;
  if (nVertex != NULL) delete[] nVertex;
  if (nNewElem_Bound != NULL) delete[] nNewElem_Bound;
  if (Marker_All_SendRecv != NULL) delete[] Marker_All_SendRecv;
  if (Tag_to_Marker != NULL) delete[] Tag_to_Marker;
  
  //	PeriodicPoint[MAX_NUMBER_PERIODIC][2].~vector();
  //	PeriodicElem[MAX_NUMBER_PERIODIC].~vector();
  //	OldBoundaryElems[MAX_NUMBER_MARKER].~vector();
  //  SendTransfLocal[MAX_NUMBER_DOMAIN].~vector();
  //  ReceivedTransfLocal[MAX_NUMBER_DOMAIN].~vector();
  //	SendDomainLocal[MAX_NUMBER_DOMAIN].~vector();
  //	ReceivedDomainLocal[MAX_NUMBER_DOMAIN].~vector();
  //	XCoordList.~vector();
  
  //	Xcoord_plane.~vector()
  //	Ycoord_plane.~vector()
  //	Zcoord_plane.~vector()
  //	FaceArea_plane.~vector()
  //	Plane_points.~vector()
  
}

double CGeometry::Point2Plane_Distance(double *Coord, double *iCoord, double *jCoord, double *kCoord) {
  double CrossProduct[3], iVector[3], jVector[3], distance, modulus;
  unsigned short iDim;
  
  for (iDim = 0; iDim < 3; iDim ++) {
    iVector[iDim] = jCoord[iDim] - iCoord[iDim];
    jVector[iDim] = kCoord[iDim] - iCoord[iDim];
  }
  
  CrossProduct[0] = iVector[1]*jVector[2] - iVector[2]*jVector[1];
  CrossProduct[1] = iVector[2]*jVector[0] - iVector[0]*jVector[2];
  CrossProduct[2] = iVector[0]*jVector[1] - iVector[1]*jVector[0];
  
  modulus = sqrt(CrossProduct[0]*CrossProduct[0]+CrossProduct[1]*CrossProduct[1]+CrossProduct[2]*CrossProduct[2]);
  
  distance = 0.0;
  for (iDim = 0; iDim < 3; iDim ++)
    distance += CrossProduct[iDim]*(Coord[iDim]-iCoord[iDim]);
  distance /= modulus;
  
  return distance;
  
}

long CGeometry::FindEdge(unsigned long first_point, unsigned long second_point) {
  unsigned long iPoint = 0;
  unsigned short iNode;
  for (iNode = 0; iNode < node[first_point]->GetnPoint(); iNode++) {
    iPoint = node[first_point]->GetPoint(iNode);
    if (iPoint == second_point) break;
  }
  
  if (iPoint == second_point) return node[first_point]->GetEdge(iNode);
  else {
    cout << "\n\n   !!! Error !!!\n" << endl;
    cout <<"Can't find the edge that connects "<< first_point <<" and "<< second_point <<"."<< endl;
    exit(1);
    return -1;
  }
}

bool CGeometry::CheckEdge(unsigned long first_point, unsigned long second_point) {
  unsigned long iPoint = 0;
  unsigned short iNode;
  for (iNode = 0; iNode < node[first_point]->GetnPoint(); iNode++) {
    iPoint = node[first_point]->GetPoint(iNode);
    if (iPoint == second_point) break;
  }
  
  if (iPoint == second_point) return true;
  else return false;
  
}

void CGeometry::SetEdges(void) {
  unsigned long iPoint, jPoint;
  long iEdge;
  unsigned short jNode, iNode;
  long TestEdge = 0;
  
  nEdge = 0;
  for(iPoint = 0; iPoint < nPoint; iPoint++)
    for(iNode = 0; iNode < node[iPoint]->GetnPoint(); iNode++) {
      jPoint = node[iPoint]->GetPoint(iNode);
      for(jNode = 0; jNode < node[jPoint]->GetnPoint(); jNode++)
        if (node[jPoint]->GetPoint(jNode) == iPoint) {
          TestEdge = node[jPoint]->GetEdge(jNode);
          break;
        }
      if (TestEdge == -1) {
        node[iPoint]->SetEdge(nEdge, iNode);
        node[jPoint]->SetEdge(nEdge, jNode);
        nEdge++;
      }
    }
  
  edge = new CEdge*[nEdge];
  
  for(iPoint = 0; iPoint < nPoint; iPoint++)
    for(iNode = 0; iNode < node[iPoint]->GetnPoint(); iNode++) {
      jPoint = node[iPoint]->GetPoint(iNode);
      iEdge = FindEdge(iPoint, jPoint);
      if (iPoint < jPoint) edge[iEdge] = new CEdge(iPoint, jPoint, nDim);
    }
}

void CGeometry::SetFaces(void) {
  //	unsigned long iPoint, jPoint, iFace;
  //	unsigned short jNode, iNode;
  //	long TestFace = 0;
  //
  //	nFace = 0;
  //	for(iPoint = 0; iPoint < nPoint; iPoint++)
  //		for(iNode = 0; iNode < node[iPoint]->GetnPoint(); iNode++) {
  //			jPoint = node[iPoint]->GetPoint(iNode);
  //			for(jNode = 0; jNode < node[jPoint]->GetnPoint(); jNode++)
  //				if (node[jPoint]->GetPoint(jNode) == iPoint) {
  //					TestFace = node[jPoint]->GetFace(jNode);
  //					break;
  //				}
  //			if (TestFace == -1) {
  //				node[iPoint]->SetFace(nFace, iNode);
  //				node[jPoint]->SetFace(nFace, jNode);
  //				nFace++;
  //			}
  //		}
  //
  //	face = new CFace*[nFace];
  //
  //	for(iPoint = 0; iPoint < nPoint; iPoint++)
  //		for(iNode = 0; iNode < node[iPoint]->GetnPoint(); iNode++) {
  //			jPoint = node[iPoint]->GetPoint(iNode);
  //			iFace = FindFace(iPoint, jPoint);
  //			if (iPoint < jPoint) face[iFace] = new CFace(iPoint,jPoint,nDim);
  //		}
}

void CGeometry::TestGeometry(void) {
  
  ofstream para_file;
  
  para_file.open("test_geometry.dat", ios::out);
  
  double *Normal = new double[nDim];
  
  for(unsigned long iEdge = 0; iEdge < nEdge; iEdge++) {
    para_file << "Edge index: " << iEdge << endl;
    para_file << "   Point index: " << edge[iEdge]->GetNode(0) << "\t" << edge[iEdge]->GetNode(1) << endl;
    edge[iEdge]->GetNormal(Normal);
    para_file << "      Face normal : ";
    for(unsigned short iDim = 0; iDim < nDim; iDim++)
      para_file << Normal[iDim] << "\t";
    para_file << endl;
  }
  
  para_file << endl;
  para_file << endl;
  para_file << endl;
  para_file << endl;
  
  for(unsigned short iMarker =0; iMarker < nMarker; iMarker++) {
    para_file << "Marker index: " << iMarker << endl;
    for(unsigned long iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
      para_file << "   Vertex index: " << iVertex << endl;
      para_file << "      Point index: " << vertex[iMarker][iVertex]->GetNode() << endl;
      para_file << "      Point coordinates : ";
      for(unsigned short iDim = 0; iDim < nDim; iDim++) {
        para_file << node[vertex[iMarker][iVertex]->GetNode()]->GetCoord(iDim) << "\t";}
      para_file << endl;
      vertex[iMarker][iVertex]->GetNormal(Normal);
      para_file << "         Face normal : ";
      for(unsigned short iDim = 0; iDim < nDim; iDim++)
        para_file << Normal[iDim] << "\t";
      para_file << endl;
    }
  }
  
}

void CGeometry::SetSpline(vector<double> &x, vector<double> &y, unsigned long n, double yp1, double ypn, vector<double> &y2) {
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

double CGeometry::GetSpline(vector<double>&xa, vector<double>&ya, vector<double>&y2a, unsigned long n, double x) {
  unsigned long klo, khi, k;
  double h, b, a, y;
  
  if (x < xa[0]) x = xa[0];       // Clip max and min values
  if (x > xa[n-1]) x = xa[n-1];
  
  klo = 1;										// We will find the right place in the table by means of
  khi = n;										// bisection. This is optimal if sequential calls to this
  while (khi-klo > 1) {			// routine are at random values of x. If sequential calls
    k = (khi+klo) >> 1;				// are in order, and closely spaced, one would do better
    if (xa[k-1] > x) khi = k;		// to store previous values of klo and khi and test if
    else klo=k;							// they remain appropriate on the next call.
  }								// klo and khi now bracket the input value of x
  h = xa[khi-1] - xa[klo-1];
  if (h == 0.0) cout << "Bad xa input to routine splint" << endl;	// The xa’s must be distinct.
  a = (xa[khi-1]-x)/h;
  b = (x-xa[klo-1])/h;				// Cubic spline polynomial is now evaluated.
  y = a*ya[klo-1]+b*ya[khi-1]+((a*a*a-a)*y2a[klo-1]+(b*b*b-b)*y2a[khi-1])*(h*h)/6.0;
  
  return y;
}

unsigned short CGeometry::ComputeSegmentPlane_Intersection(double *Segment_P0, double *Segment_P1, double *Plane_P0, double *Plane_Normal, double *Intersection) {
  double u[3], v[3], Denominator, Numerator, Aux;
  unsigned short iDim;
  
  for (iDim = 0; iDim < 3; iDim++) {
    u[iDim] = Segment_P1[iDim] - Segment_P0[iDim];
    v[iDim] = Plane_P0[iDim] - Segment_P0[iDim];
  }
  
  Numerator = Plane_Normal[0]*v[0] + Plane_Normal[1]*v[1] + Plane_Normal[2]*v[2];
  Denominator = Plane_Normal[0]*u[0] + Plane_Normal[1]*u[1] + Plane_Normal[2]*u[2];
  
  if (fabs(Denominator) <= 0.0) return 0; // No intersection.
  
  Aux = Numerator / Denominator;
  
  if (Aux < 0.0 || Aux > 1.0) return 0; // No intersection.
  
  for (iDim = 0; iDim < 3; iDim++)
    Intersection[iDim] = Segment_P0[iDim] + Aux * u[iDim];
  
  
  /*--- Check that the intersection is in the segment ---*/
  for (iDim = 0; iDim < 3; iDim++) {
    u[iDim] = Segment_P0[iDim] - Intersection[iDim];
    v[iDim] = Segment_P1[iDim] - Intersection[iDim];
  }
  
  Denominator = Plane_Normal[0]*u[0] + Plane_Normal[1]*u[1] + Plane_Normal[2]*u[2];
  Numerator = Plane_Normal[0]*v[0] + Plane_Normal[1]*v[1] + Plane_Normal[2]*v[2];
  
  Aux = Numerator * Denominator;
  
  if (Aux > 0.0) return 3; // Intersection outside the segment.
  else return 1;
  
}

CPhysicalGeometry::CPhysicalGeometry() : CGeometry() {}

CPhysicalGeometry::CPhysicalGeometry(CConfig *config) : CGeometry() {
  
  /*--- Local variables and initialization ---*/
  string text_line, Marker_Tag;
  ifstream mesh_file;
  unsigned short iNode_Surface, iMarker;
  unsigned long Point_Surface, iElem_Surface;
  double Conversion_Factor = 1.0;
  nZone = 1;
  
  string val_mesh_filename = config->GetMesh_FileName();
  
  /*--- Initialize counters for local/global points & elements ---*/
  
  cout << endl <<"---------------------- Read grid file information -----------------------" << endl;
  
  Read_SU2_Format(config, val_mesh_filename);
  
  if (config->GetKind_SU2() == SU2_CFD) Conversion_Factor = config->GetConversion_Factor();
  else Conversion_Factor = 1.0;
  
  /*--- Loop over the surface element to set the boundaries ---*/
  for (iMarker = 0; iMarker < nMarker; iMarker++)
    for (iElem_Surface = 0; iElem_Surface < nElem_Bound[iMarker]; iElem_Surface++)
      for (iNode_Surface = 0; iNode_Surface < bound[iMarker][iElem_Surface]->GetnNodes(); iNode_Surface++) {
        Point_Surface = bound[iMarker][iElem_Surface]->GetNode(iNode_Surface);
        node[Point_Surface]->SetBoundary(nMarker);
        if (config->GetMarker_All_Boundary(iMarker) != SEND_RECEIVE &&
            config->GetMarker_All_Boundary(iMarker) != INTERFACE_BOUNDARY &&
            config->GetMarker_All_Boundary(iMarker) != NEARFIELD_BOUNDARY &&
            config->GetMarker_All_Boundary(iMarker) != PERIODIC_BOUNDARY)
          node[Point_Surface]->SetPhysicalBoundary(true);
      }
  
  /*--- Write a new copy of the grid in meters if requested ---*/
  if (config->GetKind_SU2() == SU2_CFD)
    if (config->GetWrite_Converted_Mesh()) {
      SetMeshFile(config,config->GetMesh_Out_FileName());
      cout.precision(4);
      cout << "Converted mesh by a factor of " << Conversion_Factor << endl;
      cout << "  and wrote to the output file: " << config->GetMesh_Out_FileName() << endl;
    }
  
}

CPhysicalGeometry::~CPhysicalGeometry(void) {
  
}

void CPhysicalGeometry::Read_SU2_Format(CConfig *config, string val_mesh_filename) {
  
  /*--- Local variables and initialization ---*/
  string text_line, Marker_Tag;
  ifstream mesh_file;
  unsigned short VTK_Type, iMarker, iChar, iCount = 0;
  unsigned long iElem_Bound = 0, iPoint = 0, ielem_div = 0, ielem = 0, *Local2Global = NULL, vnodes_edge[2], vnodes_triangle[3], vnodes_quad[4], vnodes_tetra[4], vnodes_hexa[8], vnodes_wedge[6], vnodes_pyramid[5], dummyLong, GlobalIndex, iElem;
  char cstr[200];
  double Coord_2D[2], Coord_3D[3], Conversion_Factor = 1.0;
  string::size_type position;
  bool domain_flag = false;
  
  /*--- Initialize counters for local/global points & elements ---*/
  FinestMGLevel = true;
  Global_nPoint = 0; Global_nPointDomain = 0; Global_nElem = 0;
  nelem_edge     = 0; Global_nelem_edge     = 0;
  nelem_triangle = 0; Global_nelem_triangle = 0;
  nelem_quad     = 0; Global_nelem_quad     = 0;
  nelem_tetra    = 0; Global_nelem_tetra    = 0;
  nelem_hexa     = 0; Global_nelem_hexa     = 0;
  nelem_wedge    = 0; Global_nelem_wedge    = 0;
  nelem_pyramid  = 0; Global_nelem_pyramid  = 0;
  
  /*--- Open grid file ---*/
  strcpy (cstr, val_mesh_filename.c_str());
  mesh_file.open(cstr, ios::in);
  
  /*--- Check the grid ---*/
  if (mesh_file.fail()) {
    cout << "There is no geometry file (CPhysicalGeometry)!! " << cstr << endl;
    exit(1);
  }
  
  /*--- If divided grid, we need the global index to
   perform the right element division, this is just a hack in the future we
   should read first the coordinates ---*/
  if (config->GetDivide_Element()) {
    
    unsigned long nDim_, nElem_, nPoint_;
    
    while (getline (mesh_file, text_line)) {
      
      position = text_line.find ("NDIME=",0);
      if (position != string::npos) {
        text_line.erase (0,6); nDim_ = atoi(text_line.c_str());
      }
      
      position = text_line.find ("NELEM=",0);
      if (position != string::npos) {
        text_line.erase (0,6); nElem_ = atoi(text_line.c_str());
        for (iElem = 0; iElem < nElem_; iElem ++)
          getline(mesh_file, text_line);
      }
      
      position = text_line.find ("NPOIN=",0);
      if (position != string::npos) {
        text_line.erase (0,6);
        
        /*--- Now read the number of points and ghost points. ---*/
        stringstream  stream_line(text_line);
        stream_line >> nPoint_;
        
        Local2Global = new unsigned long [nPoint_];
        
        for (iPoint = 0; iPoint < nPoint_; iPoint ++) {
          getline(mesh_file, text_line);
          istringstream point_line(text_line);
          Local2Global[iPoint] = iPoint;
        }
        
      }
    }
    
    /*--- Close and open again the grid file ---*/
    mesh_file.close();
    strcpy (cstr, val_mesh_filename.c_str());
    mesh_file.open(cstr, ios::in);
    
  }
  
  /*--- Read grid file with format SU2 ---*/
  while (getline (mesh_file,text_line)) {
    
    /*--- Read the dimension of the problem ---*/
    position = text_line.find ("NDIME=",0);
    if (position != string::npos) {
      if (domain_flag == false) {
        text_line.erase (0,6); nDim = atoi(text_line.c_str());
        if (nDim == 2) cout << "Two dimensional problem." << endl;
        if (nDim == 3) cout << "Three dimensional problem." << endl;
        domain_flag = true;
      } else {
        break;
      }
    }
    
    /*--- Read the information about inner elements ---*/
    position = text_line.find ("NELEM=",0);
    if (position != string::npos) {
      text_line.erase (0,6); nElem = atoi(text_line.c_str());
      
      Global_nElem = nElem;
      
      cout << Global_nElem << " interior elements. " << endl;
      
      /*--- Allocate space for elements ---*/
      if (!config->GetDivide_Element()) elem = new CPrimalGrid*[nElem];
      else {
        if (nDim == 2) elem = new CPrimalGrid*[2*nElem];
        if (nDim == 3) elem = new CPrimalGrid*[6*nElem];
      }
      
      unsigned short IndirectionPrism[6][6], CT_FromVTK_Prism[6], IndirectionHexa[9][9];
      unsigned short temp, zero, one, two, three, four, five, six, seven, eight, iNode, smallestNode = 0, lookupindex;
      unsigned long smallest;
      
      /*--- Indirection matrix for dividing prisms into tets, using vtk format, and conversion table for local numbering of prisms in vtk format ---*/
      CT_FromVTK_Prism[0] = 1;  CT_FromVTK_Prism[1] = 3;  CT_FromVTK_Prism[2] = 2;  CT_FromVTK_Prism[3] = 4;  CT_FromVTK_Prism[4] = 6;  CT_FromVTK_Prism[5] = 5;
      
      IndirectionPrism[0][0] = 0; IndirectionPrism[0][1] = 2; IndirectionPrism[0][2] = 1; IndirectionPrism[0][3] = 3; IndirectionPrism[0][4] = 5; IndirectionPrism[0][5] = 4;
      IndirectionPrism[1][0] = 2; IndirectionPrism[1][1] = 1; IndirectionPrism[1][2] = 0; IndirectionPrism[1][3] = 5; IndirectionPrism[1][4] = 4; IndirectionPrism[1][5] = 3;
      IndirectionPrism[2][0] = 1; IndirectionPrism[2][1] = 0; IndirectionPrism[2][2] = 2; IndirectionPrism[2][3] = 4; IndirectionPrism[2][4] = 3; IndirectionPrism[2][5] = 5;
      IndirectionPrism[3][0] = 3; IndirectionPrism[3][1] = 4; IndirectionPrism[3][2] = 5; IndirectionPrism[3][3] = 0; IndirectionPrism[3][4] = 1; IndirectionPrism[3][5] = 2;
      IndirectionPrism[4][0] = 5; IndirectionPrism[4][1] = 3; IndirectionPrism[4][2] = 4; IndirectionPrism[4][3] = 2; IndirectionPrism[4][4] = 0; IndirectionPrism[4][5] = 1;
      IndirectionPrism[5][0] = 4; IndirectionPrism[5][1] = 5; IndirectionPrism[5][2] = 3; IndirectionPrism[5][3] = 1; IndirectionPrism[5][4] = 2; IndirectionPrism[5][5] = 0;
      
      /*--- Indirection matrix for dividing hexahedron into tets ---*/
      IndirectionHexa[1][1] = 1; IndirectionHexa[1][2] = 2; IndirectionHexa[1][3] = 3; IndirectionHexa[1][4] = 4; IndirectionHexa[1][5] = 5; IndirectionHexa[1][6] = 6; IndirectionHexa[1][7] = 7; IndirectionHexa[1][8] = 8;
      IndirectionHexa[2][1] = 2; IndirectionHexa[2][2] = 1; IndirectionHexa[2][3] = 5; IndirectionHexa[2][4] = 6; IndirectionHexa[2][5] = 3; IndirectionHexa[2][6] = 4; IndirectionHexa[2][7] = 8; IndirectionHexa[2][8] = 7;
      IndirectionHexa[3][1] = 3; IndirectionHexa[3][2] = 2; IndirectionHexa[3][3] = 6; IndirectionHexa[3][4] = 7; IndirectionHexa[3][5] = 4; IndirectionHexa[3][6] = 1; IndirectionHexa[3][7] = 5; IndirectionHexa[3][8] = 8;
      IndirectionHexa[4][1] = 4; IndirectionHexa[4][2] = 1; IndirectionHexa[4][3] = 2; IndirectionHexa[4][4] = 3; IndirectionHexa[4][5] = 8; IndirectionHexa[4][6] = 5; IndirectionHexa[4][7] = 6; IndirectionHexa[4][8] = 7;
      IndirectionHexa[5][1] = 5; IndirectionHexa[5][2] = 1; IndirectionHexa[5][3] = 4; IndirectionHexa[5][4] = 8; IndirectionHexa[5][5] = 6; IndirectionHexa[5][6] = 2; IndirectionHexa[5][7] = 3; IndirectionHexa[5][8] = 7;
      IndirectionHexa[6][1] = 6; IndirectionHexa[6][2] = 2; IndirectionHexa[6][3] = 1; IndirectionHexa[6][4] = 5; IndirectionHexa[6][5] = 7; IndirectionHexa[6][6] = 3; IndirectionHexa[6][7] = 4; IndirectionHexa[6][8] = 8;
      IndirectionHexa[7][1] = 7; IndirectionHexa[7][2] = 3; IndirectionHexa[7][3] = 2; IndirectionHexa[7][4] = 6; IndirectionHexa[7][5] = 8; IndirectionHexa[7][6] = 4; IndirectionHexa[7][7] = 1; IndirectionHexa[7][8] = 5;
      IndirectionHexa[8][1] = 8; IndirectionHexa[8][2] = 4; IndirectionHexa[8][3] = 3; IndirectionHexa[8][4] = 7; IndirectionHexa[8][5] = 5; IndirectionHexa[8][6] = 1; IndirectionHexa[8][7] = 2; IndirectionHexa[8][8] = 6;
      
      
      /*--- Loop over all the volumetric elements ---*/
      while (ielem_div < nElem) {
        getline(mesh_file,text_line);
        istringstream elem_line(text_line);
        
        elem_line >> VTK_Type;
        
        switch(VTK_Type) {
          case TRIANGLE:
            
            elem_line >> vnodes_triangle[0]; elem_line >> vnodes_triangle[1]; elem_line >> vnodes_triangle[2];
            elem[ielem] = new CTriangle(vnodes_triangle[0], vnodes_triangle[1], vnodes_triangle[2], 2);
            ielem_div++; ielem++; nelem_triangle++;
            break;
            
          case RECTANGLE:
            
            elem_line >> vnodes_quad[0]; elem_line >> vnodes_quad[1]; elem_line >> vnodes_quad[2]; elem_line >> vnodes_quad[3];
            if (!config->GetDivide_Element()) {
              elem[ielem] = new CRectangle(vnodes_quad[0], vnodes_quad[1], vnodes_quad[2], vnodes_quad[3], 2);
              ielem++; nelem_quad++; }
            else {
              
              unsigned long index1, index2;
              bool division1 = false;
              index1 = min(Local2Global[vnodes_quad[0]], Local2Global[vnodes_quad[2]]);
              index2 = min(Local2Global[vnodes_quad[1]], Local2Global[vnodes_quad[3]]);
              if (index1 < index2) division1 = true;
              
              if (division1) {
                elem[ielem] = new CTriangle(vnodes_quad[0], vnodes_quad[2], vnodes_quad[1], 2);
                ielem++; nelem_triangle++;
                elem[ielem] = new CTriangle(vnodes_quad[3], vnodes_quad[2], vnodes_quad[0], 2);
                ielem++; nelem_triangle++;
              }
              else {
                elem[ielem] = new CTriangle(vnodes_quad[2], vnodes_quad[1], vnodes_quad[3], 2);
                ielem++; nelem_triangle++;
                elem[ielem] = new CTriangle(vnodes_quad[1], vnodes_quad[0], vnodes_quad[3], 2);
                ielem++; nelem_triangle++;
              }
              
            }
            
            ielem_div++;
            break;
            
          case TETRAHEDRON:
            
            elem_line >> vnodes_tetra[0]; elem_line >> vnodes_tetra[1]; elem_line >> vnodes_tetra[2]; elem_line >> vnodes_tetra[3];
            elem[ielem] = new CTetrahedron(vnodes_tetra[0], vnodes_tetra[1], vnodes_tetra[2], vnodes_tetra[3]);
            ielem_div++; ielem++; nelem_tetra++;
            break;
            
          case HEXAHEDRON:
            
            elem_line >> vnodes_hexa[0]; elem_line >> vnodes_hexa[1]; elem_line >> vnodes_hexa[2];
            elem_line >> vnodes_hexa[3]; elem_line >> vnodes_hexa[4]; elem_line >> vnodes_hexa[5];
            elem_line >> vnodes_hexa[6]; elem_line >> vnodes_hexa[7];
            
            if (!config->GetDivide_Element()) {
              elem[ielem] = new CHexahedron(vnodes_hexa[0], vnodes_hexa[1], vnodes_hexa[2], vnodes_hexa[3],
                                            vnodes_hexa[4], vnodes_hexa[5], vnodes_hexa[6], vnodes_hexa[7]);
              ielem++; nelem_hexa++;
            }
            else {
              
              smallest = Local2Global[vnodes_hexa[0]]; smallestNode = 0;
              for (iNode = 1; iNode < 8; iNode ++) {
                if ( smallest > Local2Global[vnodes_hexa[iNode]]) {
                  smallest = Local2Global[vnodes_hexa[iNode]];
                  smallestNode = iNode;
                }
              }
              
              one  = IndirectionHexa[smallestNode+1][1] - 1;
              two  = IndirectionHexa[smallestNode+1][2] - 1;
              three  = IndirectionHexa[smallestNode+1][3] - 1;
              four = IndirectionHexa[smallestNode+1][4] - 1;
              five = IndirectionHexa[smallestNode+1][5] - 1;
              six = IndirectionHexa[smallestNode+1][6] - 1;
              seven = IndirectionHexa[smallestNode+1][7] - 1;
              eight = IndirectionHexa[smallestNode+1][8] - 1;
              
              unsigned long index1, index2;
              unsigned short code1 = 0, code2 = 0, code3 = 0;
              
              index1 = min(Local2Global[vnodes_hexa[two]], Local2Global[vnodes_hexa[seven]]);
              index2 = min(Local2Global[vnodes_hexa[three]], Local2Global[vnodes_hexa[six]]);
              if (index1 < index2) code1 = 1;
              
              index1 = min(Local2Global[vnodes_hexa[four]], Local2Global[vnodes_hexa[seven]]);
              index2 = min(Local2Global[vnodes_hexa[three]], Local2Global[vnodes_hexa[eight]]);
              if (index1 < index2) code2 = 1;
              
              index1 = min(Local2Global[vnodes_hexa[five]], Local2Global[vnodes_hexa[seven]]);
              index2 = min(Local2Global[vnodes_hexa[six]], Local2Global[vnodes_hexa[eight]]);
              if (index1 < index2) code3 = 1;
              
              /*--- Rotation of 120 degrees ---*/
              if ((!code1 && !code2 && code3) || (code1 && code2 && !code3)) {
                temp = two; two = five; five = four; four = temp;
                temp = six; six = eight; eight = three; three = temp;
              }
              
              /*--- Rotation of 240 degrees ---*/
              if ((!code1 && code2 && !code3) || (code1 && !code2 && code3)) {
                temp = two; two = four; four = five; five = temp;
                temp = six; six = three; three = eight; eight = temp;
              }
              
              if ((code1 + code2 + code3) == 0) {
                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[two], vnodes_hexa[three], vnodes_hexa[six]);
                ielem++; nelem_tetra++;
                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[three], vnodes_hexa[eight], vnodes_hexa[six]);
                ielem++; nelem_tetra++;
                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[three], vnodes_hexa[four], vnodes_hexa[eight]);
                ielem++; nelem_tetra++;
                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[six], vnodes_hexa[eight], vnodes_hexa[five]);
                ielem++; nelem_tetra++;
                elem[ielem] = new CTetrahedron(vnodes_hexa[three], vnodes_hexa[eight], vnodes_hexa[six], vnodes_hexa[seven]);
                ielem++; nelem_tetra++;
              }
              if ((code1 + code2 + code3) == 1) {
                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[six], vnodes_hexa[eight], vnodes_hexa[five]);
                ielem++; nelem_tetra++;
                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[two], vnodes_hexa[eight], vnodes_hexa[six]);
                ielem++; nelem_tetra++;
                elem[ielem] = new CTetrahedron(vnodes_hexa[two], vnodes_hexa[seven], vnodes_hexa[eight], vnodes_hexa[six]);
                ielem++; nelem_tetra++;
                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[eight], vnodes_hexa[three], vnodes_hexa[four]);
                ielem++; nelem_tetra++;
                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[eight], vnodes_hexa[two], vnodes_hexa[three]);
                ielem++; nelem_tetra++;
                elem[ielem] = new CTetrahedron(vnodes_hexa[two], vnodes_hexa[eight], vnodes_hexa[seven], vnodes_hexa[three]);
                ielem++; nelem_tetra++;
                
                //                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[six], vnodes_hexa[eight], vnodes_hexa[five]);
                //                ielem++; nelem_tetra++;
                //                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[two], vnodes_hexa[seven], vnodes_hexa[six]);
                //                ielem++; nelem_tetra++;
                //                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[seven], vnodes_hexa[eight], vnodes_hexa[six]);
                //                ielem++; nelem_tetra++;
                //                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[eight], vnodes_hexa[three], vnodes_hexa[four]);
                //                ielem++; nelem_tetra++;
                //                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[eight], vnodes_hexa[seven], vnodes_hexa[three]);
                //                ielem++; nelem_tetra++;
                //                elem[ielem] = new CTetrahedron(vnodes_hexa[two], vnodes_hexa[one], vnodes_hexa[seven], vnodes_hexa[three]);
                //                ielem++; nelem_tetra++;
                
              }
              if ((code1 + code2 + code3) == 2) {
                //                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[five], vnodes_hexa[six], vnodes_hexa[seven]);
                //                ielem++; nelem_tetra++;
                //                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[four], vnodes_hexa[eight], vnodes_hexa[seven]);
                //                ielem++; nelem_tetra++;
                //                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[eight], vnodes_hexa[five], vnodes_hexa[seven]);
                //                ielem++; nelem_tetra++;
                //                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[two], vnodes_hexa[three], vnodes_hexa[six]);
                //                ielem++; nelem_tetra++;
                //                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[four], vnodes_hexa[seven], vnodes_hexa[three]);
                //                ielem++; nelem_tetra++;
                //                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[seven], vnodes_hexa[six], vnodes_hexa[three]);
                //                ielem++; nelem_tetra++;
                
                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[three], vnodes_hexa[four], vnodes_hexa[seven]);
                ielem++; nelem_tetra++;
                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[five], vnodes_hexa[seven], vnodes_hexa[eight]);
                ielem++; nelem_tetra++;
                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[seven], vnodes_hexa[four], vnodes_hexa[eight]);
                ielem++; nelem_tetra++;
                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[two], vnodes_hexa[three], vnodes_hexa[six]);
                ielem++; nelem_tetra++;
                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[seven], vnodes_hexa[five], vnodes_hexa[six]);
                ielem++; nelem_tetra++;
                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[three], vnodes_hexa[seven], vnodes_hexa[six]);
                ielem++; nelem_tetra++;
              }
              if ((code1 + code2 + code3) == 3) {
                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[three], vnodes_hexa[four], vnodes_hexa[seven]);
                ielem++; nelem_tetra++;
                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[four], vnodes_hexa[eight], vnodes_hexa[seven]);
                ielem++; nelem_tetra++;
                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[eight], vnodes_hexa[five], vnodes_hexa[seven]);
                ielem++; nelem_tetra++;
                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[six], vnodes_hexa[seven], vnodes_hexa[five]);
                ielem++; nelem_tetra++;
                elem[ielem] = new CTetrahedron(vnodes_hexa[two], vnodes_hexa[six], vnodes_hexa[seven], vnodes_hexa[one]);
                ielem++; nelem_tetra++;
                elem[ielem] = new CTetrahedron(vnodes_hexa[two], vnodes_hexa[seven], vnodes_hexa[three], vnodes_hexa[one]);
                ielem++; nelem_tetra++;
                
                //                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[two], vnodes_hexa[seven], vnodes_hexa[six]);
                //                ielem++; nelem_tetra++;
                //                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[seven], vnodes_hexa[eight], vnodes_hexa[five]);
                //                ielem++; nelem_tetra++;
                //                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[six], vnodes_hexa[seven], vnodes_hexa[five]);
                //                ielem++; nelem_tetra++;
                //                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[seven], vnodes_hexa[two], vnodes_hexa[three]);
                //                ielem++; nelem_tetra++;
                //                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[eight], vnodes_hexa[seven], vnodes_hexa[four]);
                //                ielem++; nelem_tetra++;
                //                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[seven], vnodes_hexa[three], vnodes_hexa[four]);
                //                ielem++; nelem_tetra++;
                
                //                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[two], vnodes_hexa[seven], vnodes_hexa[six]);
                //                ielem++; nelem_tetra++;
                //                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[two], vnodes_hexa[three], vnodes_hexa[seven]);
                //                ielem++; nelem_tetra++;
                //                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[three], vnodes_hexa[four], vnodes_hexa[seven]);
                //                ielem++; nelem_tetra++;
                //                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[six], vnodes_hexa[seven], vnodes_hexa[five]);
                //                ielem++; nelem_tetra++;
                //                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[seven], vnodes_hexa[four], vnodes_hexa[eight]);
                //                ielem++; nelem_tetra++;
                //                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[seven], vnodes_hexa[eight], vnodes_hexa[five]);
                //                ielem++; nelem_tetra++;
              }
              
            }
            ielem_div++;
            break;
            
          case WEDGE:
            
            elem_line >> vnodes_wedge[0]; elem_line >> vnodes_wedge[1]; elem_line >> vnodes_wedge[2];
            elem_line >> vnodes_wedge[3]; elem_line >> vnodes_wedge[4]; elem_line >> vnodes_wedge[5];
            
            if (!config->GetDivide_Element()) {
              
              elem[ielem] = new CWedge(vnodes_wedge[0],vnodes_wedge[1],vnodes_wedge[2],vnodes_wedge[3],vnodes_wedge[4],vnodes_wedge[5]);
              ielem++; nelem_wedge++;
              
            }
            else {
              
              smallest = Local2Global[vnodes_wedge[0]]; smallestNode = 0;
              for (iNode = 1; iNode < 6; iNode ++) {
                if ( smallest > Local2Global[vnodes_wedge[iNode]]) {
                  smallest = Local2Global[vnodes_wedge[iNode]];
                  smallestNode = iNode;
                }
              }
              
              lookupindex = (CT_FromVTK_Prism[smallestNode] - 1);
              zero  = IndirectionPrism[lookupindex][0];
              one  = IndirectionPrism[lookupindex][1];
              two  = IndirectionPrism[lookupindex][2];
              three = IndirectionPrism[lookupindex][3];
              four = IndirectionPrism[lookupindex][4];
              five = IndirectionPrism[lookupindex][5];
              
              unsigned long index1, index2;
              bool division = false;
              index1 = min(Local2Global[vnodes_wedge[one]], Local2Global[vnodes_wedge[five]]);
              index2 = min(Local2Global[vnodes_wedge[two]], Local2Global[vnodes_wedge[four]]);
              if (index1 < index2) division = true;
              
              if (division) {
                elem[ielem] = new CTetrahedron(vnodes_wedge[zero], vnodes_wedge[one], vnodes_wedge[two], vnodes_wedge[five]);
                ielem++; nelem_tetra++;
                elem[ielem] = new CTetrahedron(vnodes_wedge[zero], vnodes_wedge[one], vnodes_wedge[five], vnodes_wedge[four]);
                ielem++; nelem_tetra++;
                elem[ielem] = new CTetrahedron(vnodes_wedge[zero], vnodes_wedge[four], vnodes_wedge[five], vnodes_wedge[three]);
                ielem++; nelem_tetra++;
              }
              else {
                elem[ielem] = new CTetrahedron(vnodes_wedge[zero], vnodes_wedge[one], vnodes_wedge[two], vnodes_wedge[four]);
                ielem++; nelem_tetra++;
                elem[ielem] = new CTetrahedron(vnodes_wedge[zero], vnodes_wedge[four], vnodes_wedge[two], vnodes_wedge[five]);
                ielem++; nelem_tetra++;
                elem[ielem] = new CTetrahedron(vnodes_wedge[zero], vnodes_wedge[four], vnodes_wedge[five], vnodes_wedge[three]);
                ielem++; nelem_tetra++;
              }
            }
            
            ielem_div++;
            break;
          case PYRAMID:
            
            elem_line >> vnodes_pyramid[0]; elem_line >> vnodes_pyramid[1]; elem_line >> vnodes_pyramid[2];
            elem_line >> vnodes_pyramid[3]; elem_line >> vnodes_pyramid[4];
            
            if (!config->GetDivide_Element()) {
              
              elem[ielem] = new CPyramid(vnodes_pyramid[0],vnodes_pyramid[1],vnodes_pyramid[2],vnodes_pyramid[3],vnodes_pyramid[4]);
              ielem++; nelem_pyramid++;
              
            }
            else {
              
              unsigned long index1, index2;
              bool division = false;
              index1 = min(Local2Global[vnodes_pyramid[0]], Local2Global[vnodes_pyramid[2]]);
              index2 = min(Local2Global[vnodes_pyramid[1]], Local2Global[vnodes_pyramid[3]]);
              if (index1 < index2) division = true;
              
              if (division) {
                elem[ielem] = new CTetrahedron(vnodes_pyramid[0], vnodes_pyramid[1], vnodes_pyramid[2], vnodes_pyramid[4]);
                ielem++; nelem_tetra++;
                elem[ielem] = new CTetrahedron(vnodes_pyramid[0], vnodes_pyramid[2], vnodes_pyramid[3], vnodes_pyramid[4]);
                ielem++; nelem_tetra++;
              }
              else {
                elem[ielem] = new CTetrahedron(vnodes_pyramid[1], vnodes_pyramid[2], vnodes_pyramid[3], vnodes_pyramid[4]);
                ielem++; nelem_tetra++;
                elem[ielem] = new CTetrahedron(vnodes_pyramid[1], vnodes_pyramid[3], vnodes_pyramid[0], vnodes_pyramid[4]);
                ielem++; nelem_tetra++;
              }
              
            }
            
            ielem_div++;
            break;
        }
      }
      
      if (config->GetDivide_Element()) nElem = nelem_triangle + nelem_quad + nelem_tetra + nelem_hexa + nelem_wedge + nelem_pyramid;
      
      /*--- Communicate the number of each element type to all processors. ---*/
      Global_nelem_triangle = nelem_triangle;
      Global_nelem_quad     = nelem_quad;
      Global_nelem_tetra    = nelem_tetra;
      Global_nelem_hexa     = nelem_hexa;
      Global_nelem_wedge    = nelem_wedge;
      Global_nelem_pyramid  = nelem_pyramid;
      
      /*--- Print information about the elements to the console ---*/
      if (Global_nelem_triangle > 0)
        cout << Global_nelem_triangle << " triangles." << endl;
      if (Global_nelem_quad > 0)
        cout << Global_nelem_quad << " quadrilaterals." << endl;
      if (Global_nelem_tetra > 0)
        cout << Global_nelem_tetra << " tetrahedra." << endl;
      if (Global_nelem_hexa > 0)
        cout << Global_nelem_hexa << " hexahedra." << endl;
      if (Global_nelem_wedge > 0)
        cout << Global_nelem_wedge << " prisms." << endl;
      if (Global_nelem_pyramid > 0)
        cout << Global_nelem_pyramid << " pyramids." << endl;
    }
    
    /*--- Read number of points ---*/
    position = text_line.find ("NPOIN=",0);
    if (position != string::npos) {
      text_line.erase (0,6);
      
      /*--- Check for ghost points. ---*/
      stringstream test_line(text_line);
      while (test_line >> dummyLong)
        iCount++;
      
      /*--- Now read and store the number of points and possible ghost points. ---*/
      stringstream  stream_line(text_line);
      if (iCount == 2) {
        stream_line >> nPoint;
        stream_line >> nPointDomain;
        
        /*--- Set some important point information for parallel simulations. ---*/
        Global_nPoint = nPoint;
        Global_nPointDomain = nPointDomain;
        
        cout << Global_nPointDomain << " points, and " << Global_nPoint-Global_nPointDomain << " ghost points." << endl;
      }
      else if (iCount == 1) {
        stream_line >> nPoint;
        nPointDomain = nPoint;
        Global_nPointDomain = nPoint;
        Global_nPoint = nPoint;
        cout << nPoint << " points." << endl;
      }
      else {
        cout << "NPOIN improperly specified!!" << endl;
        exit(1);
      }
      
      /*--- Retrieve grid conversion factor. The conversion is only
       applied for SU2_CFD. All other SU2 components leave the mesh
       as is. ---*/
      if (config->GetKind_SU2() == SU2_CFD)
        Conversion_Factor = config->GetConversion_Factor();
      else
        Conversion_Factor = 1.0;
      
      node = new CPoint*[nPoint];
      iPoint = 0;
      while (iPoint < nPoint) {
        getline(mesh_file,text_line);
        istringstream point_line(text_line);
        switch(nDim) {
          case 2:
            GlobalIndex = iPoint;
            point_line >> Coord_2D[0]; point_line >> Coord_2D[1];
            node[iPoint] = new CPoint(Conversion_Factor*Coord_2D[0], Conversion_Factor*Coord_2D[1], GlobalIndex, config);
            iPoint++; break;
          case 3:
            GlobalIndex = iPoint;
            point_line >> Coord_3D[0]; point_line >> Coord_3D[1]; point_line >> Coord_3D[2];
            node[iPoint] = new CPoint(Conversion_Factor*Coord_3D[0], Conversion_Factor*Coord_3D[1], Conversion_Factor*Coord_3D[2], GlobalIndex, config);
            iPoint++; break;
        }
      }
    }
    
    
    /*--- Read number of markers ---*/
    position = text_line.find ("NMARK=",0);
    if (position != string::npos) {
      text_line.erase (0,6); nMarker = atoi(text_line.c_str());
      cout << nMarker << " surface markers." << endl;
      config->SetnMarker_All(nMarker);
      bound = new CPrimalGrid**[nMarker];
      nElem_Bound = new unsigned long [nMarker];
      Tag_to_Marker = new string [MAX_INDEX_VALUE];
      
      for (iMarker = 0 ; iMarker < nMarker; iMarker++) {
        getline (mesh_file,text_line);
        text_line.erase (0,11);
        string::size_type position;
        for (iChar = 0; iChar < 20; iChar++) {
          position = text_line.find( " ", 0 );
          if(position != string::npos) text_line.erase (position,1);
          position = text_line.find( "\r", 0 );
          if(position != string::npos) text_line.erase (position,1);
          position = text_line.find( "\n", 0 );
          if(position != string::npos) text_line.erase (position,1);
        }
        Marker_Tag = text_line.c_str();
        
        /*--- Physical boundaries definition ---*/
        if (Marker_Tag != "SEND_RECEIVE") {
          getline (mesh_file,text_line);
          text_line.erase (0,13); nElem_Bound[iMarker] = atoi(text_line.c_str());
          cout << nElem_Bound[iMarker]  << " boundary elements in index "<< iMarker <<" (Marker = " <<Marker_Tag<< ")." << endl;
          
          /*--- Allocate space for elements ---*/
          if (!config->GetDivide_Element()) bound[iMarker] = new CPrimalGrid* [nElem_Bound[iMarker]];
          else {
            if (nDim == 2) bound[iMarker] = new CPrimalGrid* [2*nElem_Bound[iMarker]];;
            if (nDim == 3) bound[iMarker] = new CPrimalGrid* [2*nElem_Bound[iMarker]];;
          }
          
          nelem_edge_bound = 0; nelem_triangle_bound = 0; nelem_quad_bound = 0; ielem = 0;
          for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
            getline(mesh_file,text_line);
            istringstream bound_line(text_line);
            bound_line >> VTK_Type;
            switch(VTK_Type) {
              case LINE:
                
                if (nDim == 3) {
                  cout << "Please remove line boundary conditions from the mesh file!" << endl;
                  exit(1);
                }
                
                bound_line >> vnodes_edge[0]; bound_line >> vnodes_edge[1];
                bound[iMarker][ielem] = new CLine(vnodes_edge[0],vnodes_edge[1],2);
                ielem++; nelem_edge_bound++; break;
                
              case TRIANGLE:
                bound_line >> vnodes_triangle[0]; bound_line >> vnodes_triangle[1]; bound_line >> vnodes_triangle[2];
                bound[iMarker][ielem] = new CTriangle(vnodes_triangle[0],vnodes_triangle[1],vnodes_triangle[2],3);
                ielem++; nelem_triangle_bound++; break;
                
              case RECTANGLE:
                
                bound_line >> vnodes_quad[0]; bound_line >> vnodes_quad[1]; bound_line >> vnodes_quad[2]; bound_line >> vnodes_quad[3];
                
                if (!config->GetDivide_Element()) {
                  
                  bound[iMarker][ielem] = new CRectangle(vnodes_quad[0],vnodes_quad[1],vnodes_quad[2],vnodes_quad[3],3);
                  ielem++; nelem_quad_bound++;
                  
                }
                else {
                  
                  unsigned long index1, index2;
                  bool division1 = false;
                  index1 = min(Local2Global[vnodes_quad[0]], Local2Global[vnodes_quad[2]]);
                  index2 = min(Local2Global[vnodes_quad[1]], Local2Global[vnodes_quad[3]]);
                  if (index1 < index2) division1 = true;
                  
                  if (division1) {
                    bound[iMarker][ielem] = new CTriangle(vnodes_quad[0], vnodes_quad[2], vnodes_quad[1], 3);
                    ielem++; nelem_triangle_bound++;
                    bound[iMarker][ielem] = new CTriangle(vnodes_quad[3], vnodes_quad[2], vnodes_quad[0], 3);
                    ielem++; nelem_triangle_bound++;
                  }
                  else {
                    bound[iMarker][ielem] = new CTriangle(vnodes_quad[2], vnodes_quad[1], vnodes_quad[3], 3);
                    ielem++; nelem_triangle_bound++;
                    bound[iMarker][ielem] = new CTriangle(vnodes_quad[1], vnodes_quad[0], vnodes_quad[3], 3);
                    ielem++; nelem_triangle_bound++;
                  }
                  
                }
                
                break;
                
                
            }
          }
          if (config->GetDivide_Element()) nElem_Bound[iMarker] = nelem_edge_bound + nelem_triangle_bound + nelem_quad_bound;
          
          /*--- Update config information storing the boundary information in the right place ---*/
          Tag_to_Marker[config->GetMarker_Config_Tag(Marker_Tag)] = Marker_Tag;
          config->SetMarker_All_Tag(iMarker, Marker_Tag);
          config->SetMarker_All_Boundary(iMarker, config->GetMarker_Config_Boundary(Marker_Tag));
          config->SetMarker_All_Monitoring(iMarker, config->GetMarker_Config_Monitoring(Marker_Tag));
          config->SetMarker_All_Designing(iMarker, config->GetMarker_Config_Designing(Marker_Tag));
          config->SetMarker_All_Plotting(iMarker, config->GetMarker_Config_Plotting(Marker_Tag));
          config->SetMarker_All_DV(iMarker, config->GetMarker_Config_DV(Marker_Tag));
          config->SetMarker_All_Moving(iMarker, config->GetMarker_Config_Moving(Marker_Tag));
          config->SetMarker_All_PerBound(iMarker, config->GetMarker_Config_PerBound(Marker_Tag));
          config->SetMarker_All_SendRecv(iMarker, NONE);
          
        }
        
        /*--- Send-Receive boundaries definition ---*/
        else {
          unsigned long nelem_vertex = 0, vnodes_vertex;
          unsigned short transform;
          getline (mesh_file,text_line);
          text_line.erase (0,13); nElem_Bound[iMarker] = atoi(text_line.c_str());
          bound[iMarker] = new CPrimalGrid* [nElem_Bound[iMarker]];
          
          nelem_vertex = 0; ielem = 0;
          getline (mesh_file,text_line); text_line.erase (0,8);
          config->SetMarker_All_Boundary(iMarker, SEND_RECEIVE);
          config->SetMarker_All_SendRecv(iMarker, atoi(text_line.c_str()));
          
          for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
            getline(mesh_file,text_line);
            istringstream bound_line(text_line);
            bound_line >> VTK_Type; bound_line >> vnodes_vertex; bound_line >> transform;
            
            bound[iMarker][ielem] = new CVertexMPI(vnodes_vertex, nDim);
            bound[iMarker][ielem]->SetRotation_Type(transform);
            ielem++; nelem_vertex++;
            if (config->GetMarker_All_SendRecv(iMarker) < 0)
              node[vnodes_vertex]->SetDomain(false);
          }
          
        }
        
      }
    }
    
  }
  
  /*--- Close the input file ---*/
  mesh_file.close();
  
  if (config->GetDivide_Element()) {
    if (Local2Global != NULL) delete [] Local2Global;
  }
  
}

void CPhysicalGeometry::Check_Orientation(CConfig *config) {
  unsigned long Point_1, Point_2, Point_3, Point_4, Point_5, Point_6,
  iElem, Point_1_Surface, Point_2_Surface, Point_3_Surface, Point_4_Surface,
  iElem_Domain, Point_Domain = 0, Point_Surface, iElem_Surface;
  double test_1, test_2, test_3, test_4, *Coord_1, *Coord_2, *Coord_3, *Coord_4,
  *Coord_5, *Coord_6, a[3], b[3], c[3], n[3], test;
  unsigned short iDim, iMarker, iNode_Domain, iNode_Surface;
  bool find;
  
  /*--- Loop over all the elements ---*/
  for (iElem = 0; iElem < nElem; iElem++) {
    
    /*--- 2D grid, triangle case ---*/
    if (elem[iElem]->GetVTK_Type() == TRIANGLE) {
      
      Point_1 = elem[iElem]->GetNode(0); Coord_1 = node[Point_1]->GetCoord();
      Point_2 = elem[iElem]->GetNode(1); Coord_2 = node[Point_2]->GetCoord();
      Point_3 = elem[iElem]->GetNode(2); Coord_3 = node[Point_3]->GetCoord();
      
      for(iDim = 0; iDim < nDim; iDim++) {
        a[iDim] = 0.5*(Coord_2[iDim]-Coord_1[iDim]);
        b[iDim] = 0.5*(Coord_3[iDim]-Coord_1[iDim]); }
      test = a[0]*b[1]-b[0]*a[1];
      
      if (test < 0.0) elem[iElem]->Change_Orientation();
    }
    
    /*--- 2D grid, rectangle case ---*/
    if (elem[iElem]->GetVTK_Type() == RECTANGLE) {
      
      Point_1 = elem[iElem]->GetNode(0); Coord_1 = node[Point_1]->GetCoord();
      Point_2 = elem[iElem]->GetNode(1); Coord_2 = node[Point_2]->GetCoord();
      Point_3 = elem[iElem]->GetNode(2); Coord_3 = node[Point_3]->GetCoord();
      Point_4 = elem[iElem]->GetNode(3); Coord_4 = node[Point_4]->GetCoord();
      
      for(iDim = 0; iDim < nDim; iDim++) {
        a[iDim] = 0.5*(Coord_2[iDim]-Coord_1[iDim]);
        b[iDim] = 0.5*(Coord_3[iDim]-Coord_1[iDim]); }
      test_1 = a[0]*b[1]-b[0]*a[1];
      
      for(iDim = 0; iDim < nDim; iDim++) {
        a[iDim] = 0.5*(Coord_3[iDim]-Coord_2[iDim]);
        b[iDim] = 0.5*(Coord_4[iDim]-Coord_2[iDim]); }
      test_2 = a[0]*b[1]-b[0]*a[1];
      
      for(iDim = 0; iDim < nDim; iDim++) {
        a[iDim] = 0.5*(Coord_4[iDim]-Coord_3[iDim]);
        b[iDim] = 0.5*(Coord_1[iDim]-Coord_3[iDim]); }
      test_3 = a[0]*b[1]-b[0]*a[1];
      
      for(iDim = 0; iDim < nDim; iDim++) {
        a[iDim] = 0.5*(Coord_1[iDim]-Coord_4[iDim]);
        b[iDim] = 0.5*(Coord_3[iDim]-Coord_4[iDim]); }
      test_4 = a[0]*b[1]-b[0]*a[1];
      
      if ((test_1 < 0.0) && (test_2 < 0.0) && (test_3 < 0.0) && (test_4 < 0.0))
        elem[iElem]->Change_Orientation();
    }
    
    /*--- 3D grid, tetrahedron case ---*/
    if (elem[iElem]->GetVTK_Type() == TETRAHEDRON) {
      
      Point_1 = elem[iElem]->GetNode(0); Coord_1 = node[Point_1]->GetCoord();
      Point_2 = elem[iElem]->GetNode(1); Coord_2 = node[Point_2]->GetCoord();
      Point_3 = elem[iElem]->GetNode(2); Coord_3 = node[Point_3]->GetCoord();
      Point_4 = elem[iElem]->GetNode(3); Coord_4 = node[Point_4]->GetCoord();
      
      for(iDim = 0; iDim < nDim; iDim++) {
        a[iDim] = 0.5*(Coord_2[iDim]-Coord_1[iDim]);
        b[iDim] = 0.5*(Coord_3[iDim]-Coord_1[iDim]);
        c[iDim] = Coord_4[iDim]-Coord_1[iDim]; }
      n[0] = a[1]*b[2]-b[1]*a[2];
      n[1] = -(a[0]*b[2]-b[0]*a[2]);
      n[2] = a[0]*b[1]-b[0]*a[1];
      
      test = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];
      if (test < 0.0) elem[iElem]->Change_Orientation();
      
    }
    
    /*--- 3D grid, wedge case ---*/
    if (elem[iElem]->GetVTK_Type() == WEDGE) {
      
      Point_1 = elem[iElem]->GetNode(0); Coord_1 = node[Point_1]->GetCoord();
      Point_2 = elem[iElem]->GetNode(1); Coord_2 = node[Point_2]->GetCoord();
      Point_3 = elem[iElem]->GetNode(2); Coord_3 = node[Point_3]->GetCoord();
      Point_4 = elem[iElem]->GetNode(3); Coord_4 = node[Point_4]->GetCoord();
      Point_5 = elem[iElem]->GetNode(4); Coord_5 = node[Point_5]->GetCoord();
      Point_6 = elem[iElem]->GetNode(5); Coord_6 = node[Point_6]->GetCoord();
      
      for(iDim = 0; iDim < nDim; iDim++) {
        a[iDim] = 0.5*(Coord_3[iDim]-Coord_1[iDim]);
        b[iDim] = 0.5*(Coord_2[iDim]-Coord_1[iDim]);
        c[iDim] = (Coord_4[iDim]-Coord_1[iDim])+
        (Coord_5[iDim]-Coord_2[iDim])+
        (Coord_6[iDim]-Coord_3[iDim]); }
      
      /*--- The normal vector should point to the interior of the element ---*/
      n[0] = a[1]*b[2]-b[1]*a[2];
      n[1] = -(a[0]*b[2]-b[0]*a[2]);
      n[2] = a[0]*b[1]-b[0]*a[1];
      
      test_1 = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];
      
      for(iDim = 0; iDim < nDim; iDim++) {
        a[iDim] = 0.5*(Coord_5[iDim]-Coord_4[iDim]);
        b[iDim] = 0.5*(Coord_6[iDim]-Coord_4[iDim]);
        c[iDim] = (Coord_1[iDim]-Coord_4[iDim])+
        (Coord_2[iDim]-Coord_5[iDim])+
        (Coord_3[iDim]-Coord_6[iDim]); }
      
      /*--- The normal vector should point to the interior of the element ---*/
      n[0] = a[1]*b[2]-b[1]*a[2];
      n[1] = -(a[0]*b[2]-b[0]*a[2]);
      n[2] = a[0]*b[1]-b[0]*a[1];
      
      test_2 = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];
      
      if ((test_1 < 0.0) || (test_2 < 0.0))
        elem[iElem]->Change_Orientation();
      
    }
    
    if (elem[iElem]->GetVTK_Type() == HEXAHEDRON) {
      
      Point_1 = elem[iElem]->GetNode(0); Coord_1 = node[Point_1]->GetCoord();
      Point_2 = elem[iElem]->GetNode(1); Coord_2 = node[Point_2]->GetCoord();
      Point_3 = elem[iElem]->GetNode(2); Coord_3 = node[Point_3]->GetCoord();
      Point_4 = elem[iElem]->GetNode(5); Coord_4 = node[Point_4]->GetCoord();
      
      for(iDim = 0; iDim < nDim; iDim++) {
        a[iDim] = 0.5*(Coord_2[iDim]-Coord_1[iDim]);
        b[iDim] = 0.5*(Coord_3[iDim]-Coord_1[iDim]);
        c[iDim] = Coord_4[iDim]-Coord_1[iDim]; }
      n[0] = a[1]*b[2]-b[1]*a[2];
      n[1] = -(a[0]*b[2]-b[0]*a[2]);
      n[2] = a[0]*b[1]-b[0]*a[1];
      
      test_1 = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];
      
      Point_1 = elem[iElem]->GetNode(2); Coord_1 = node[Point_1]->GetCoord();
      Point_2 = elem[iElem]->GetNode(3); Coord_2 = node[Point_2]->GetCoord();
      Point_3 = elem[iElem]->GetNode(0); Coord_3 = node[Point_3]->GetCoord();
      Point_4 = elem[iElem]->GetNode(7); Coord_4 = node[Point_4]->GetCoord();
      
      for(iDim = 0; iDim < nDim; iDim++) {
        a[iDim] = 0.5*(Coord_2[iDim]-Coord_1[iDim]);
        b[iDim] = 0.5*(Coord_3[iDim]-Coord_1[iDim]);
        c[iDim] = Coord_4[iDim]-Coord_1[iDim]; }
      n[0] = a[1]*b[2]-b[1]*a[2];
      n[1] = -(a[0]*b[2]-b[0]*a[2]);
      n[2] = a[0]*b[1]-b[0]*a[1];
      
      test_2 = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];
      
      Point_1 = elem[iElem]->GetNode(1); Coord_1 = node[Point_1]->GetCoord();
      Point_2 = elem[iElem]->GetNode(2); Coord_2 = node[Point_2]->GetCoord();
      Point_3 = elem[iElem]->GetNode(3); Coord_3 = node[Point_3]->GetCoord();
      Point_4 = elem[iElem]->GetNode(6); Coord_4 = node[Point_4]->GetCoord();
      
      for(iDim = 0; iDim < nDim; iDim++) {
        a[iDim] = 0.5*(Coord_2[iDim]-Coord_1[iDim]);
        b[iDim] = 0.5*(Coord_3[iDim]-Coord_1[iDim]);
        c[iDim] = Coord_4[iDim]-Coord_1[iDim]; }
      n[0] = a[1]*b[2]-b[1]*a[2];
      n[1] = -(a[0]*b[2]-b[0]*a[2]);
      n[2] = a[0]*b[1]-b[0]*a[1];
      
      test_3 = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];
      
      Point_1 = elem[iElem]->GetNode(3); Coord_1 = node[Point_1]->GetCoord();
      Point_2 = elem[iElem]->GetNode(0); Coord_2 = node[Point_2]->GetCoord();
      Point_3 = elem[iElem]->GetNode(1); Coord_3 = node[Point_3]->GetCoord();
      Point_4 = elem[iElem]->GetNode(4); Coord_4 = node[Point_4]->GetCoord();
      
      for(iDim = 0; iDim < nDim; iDim++) {
        a[iDim] = 0.5*(Coord_2[iDim]-Coord_1[iDim]);
        b[iDim] = 0.5*(Coord_3[iDim]-Coord_1[iDim]);
        c[iDim] = Coord_4[iDim]-Coord_1[iDim]; }
      n[0] = a[1]*b[2]-b[1]*a[2];
      n[1] = -(a[0]*b[2]-b[0]*a[2]);
      n[2] = a[0]*b[1]-b[0]*a[1];
      
      test_4 = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];
      
      if ((test_1 < 0.0) || (test_2 < 0.0) || (test_3 < 0.0)
          || (test_4 < 0.0)) elem[iElem]->Change_Orientation();
      
    }
    
    if (elem[iElem]->GetVTK_Type() == PYRAMID) {
      
      Point_1 = elem[iElem]->GetNode(0); Coord_1 = node[Point_1]->GetCoord();
      Point_2 = elem[iElem]->GetNode(1); Coord_2 = node[Point_2]->GetCoord();
      Point_3 = elem[iElem]->GetNode(2); Coord_3 = node[Point_3]->GetCoord();
      Point_4 = elem[iElem]->GetNode(4); Coord_4 = node[Point_4]->GetCoord();
      
      for(iDim = 0; iDim < nDim; iDim++) {
        a[iDim] = 0.5*(Coord_2[iDim]-Coord_1[iDim]);
        b[iDim] = 0.5*(Coord_3[iDim]-Coord_1[iDim]);
        c[iDim] = Coord_4[iDim]-Coord_1[iDim]; }
      n[0] = a[1]*b[2]-b[1]*a[2];
      n[1] = -(a[0]*b[2]-b[0]*a[2]);
      n[2] = a[0]*b[1]-b[0]*a[1];
      
      test_1 = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];
      
      Point_1 = elem[iElem]->GetNode(2); Coord_1 = node[Point_1]->GetCoord();
      Point_2 = elem[iElem]->GetNode(3); Coord_2 = node[Point_2]->GetCoord();
      Point_3 = elem[iElem]->GetNode(0); Coord_3 = node[Point_3]->GetCoord();
      Point_4 = elem[iElem]->GetNode(4); Coord_4 = node[Point_4]->GetCoord();
      
      for(iDim = 0; iDim < nDim; iDim++) {
        a[iDim] = 0.5*(Coord_2[iDim]-Coord_1[iDim]);
        b[iDim] = 0.5*(Coord_3[iDim]-Coord_1[iDim]);
        c[iDim] = Coord_4[iDim]-Coord_1[iDim]; }
      n[0] = a[1]*b[2]-b[1]*a[2];
      n[1] = -(a[0]*b[2]-b[0]*a[2]);
      n[2] = a[0]*b[1]-b[0]*a[1];
      
      test_2 = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];
      
      if ((test_1 < 0.0) || (test_2 < 0.0))
        elem[iElem]->Change_Orientation();
      
    }
    
  }
  
  /*--- Surface elements ---*/
  for (iMarker = 0; iMarker < nMarker; iMarker++)
    for (iElem_Surface = 0; iElem_Surface < nElem_Bound[iMarker]; iElem_Surface++) {
      
      iElem_Domain = bound[iMarker][iElem_Surface]->GetDomainElement();
      for (iNode_Domain = 0; iNode_Domain < elem[iElem_Domain]->GetnNodes(); iNode_Domain++) {
        Point_Domain = elem[iElem_Domain]->GetNode(iNode_Domain);
        find = false;
        for (iNode_Surface = 0; iNode_Surface < bound[iMarker][iElem_Surface]->GetnNodes(); iNode_Surface++) {
          Point_Surface = bound[iMarker][iElem_Surface]->GetNode(iNode_Surface);
          if (Point_Surface == Point_Domain) {find = true; break;}
        }
        if (!find) break;
      }
      
      /*--- 2D grid, line case ---*/
      if (bound[iMarker][iElem_Surface]->GetVTK_Type() == LINE) {
        
        Point_1_Surface = bound[iMarker][iElem_Surface]->GetNode(0); Coord_1 = node[Point_1_Surface]->GetCoord();
        Point_2_Surface = bound[iMarker][iElem_Surface]->GetNode(1); Coord_2 = node[Point_2_Surface]->GetCoord();
        Coord_3 = node[Point_Domain]->GetCoord();
        
        for(iDim = 0; iDim < nDim; iDim++) {
          a[iDim] = 0.5*(Coord_2[iDim]-Coord_1[iDim]);
          b[iDim] = 0.5*(Coord_3[iDim]-Coord_1[iDim]); }
        test = a[0]*b[1]-b[0]*a[1];
        
        if (test < 0.0) bound[iMarker][iElem_Surface]->Change_Orientation();
      }
      
      /*--- 3D grid, triangle case ---*/
      if (bound[iMarker][iElem_Surface]->GetVTK_Type() == TRIANGLE) {
        
        Point_1_Surface = bound[iMarker][iElem_Surface]->GetNode(0); Coord_1 = node[Point_1_Surface]->GetCoord();
        Point_2_Surface = bound[iMarker][iElem_Surface]->GetNode(1); Coord_2 = node[Point_2_Surface]->GetCoord();
        Point_3_Surface = bound[iMarker][iElem_Surface]->GetNode(2); Coord_3 = node[Point_3_Surface]->GetCoord();
        Coord_4 = node[Point_Domain]->GetCoord();
        
        for(iDim = 0; iDim < nDim; iDim++) {
          a[iDim] = 0.5*(Coord_2[iDim]-Coord_1[iDim]);
          b[iDim] = 0.5*(Coord_3[iDim]-Coord_1[iDim]);
          c[iDim] = Coord_4[iDim]-Coord_1[iDim]; }
        n[0] = a[1]*b[2]-b[1]*a[2];
        n[1] = -(a[0]*b[2]-b[0]*a[2]);
        n[2] = a[0]*b[1]-b[0]*a[1];
        
        test = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];
        if (test < 0.0) bound[iMarker][iElem_Surface]->Change_Orientation();
      }
      
      if (bound[iMarker][iElem_Surface]->GetVTK_Type() == RECTANGLE) {
        
        Point_1_Surface = bound[iMarker][iElem_Surface]->GetNode(0); Coord_1 = node[Point_1_Surface]->GetCoord();
        Point_2_Surface = bound[iMarker][iElem_Surface]->GetNode(1); Coord_2 = node[Point_2_Surface]->GetCoord();
        Point_3_Surface = bound[iMarker][iElem_Surface]->GetNode(2); Coord_3 = node[Point_3_Surface]->GetCoord();
        Point_4_Surface = bound[iMarker][iElem_Surface]->GetNode(3); Coord_4 = node[Point_4_Surface]->GetCoord();
        Coord_5 = node[Point_Domain]->GetCoord();
        
        for(iDim = 0; iDim < nDim; iDim++) {
          a[iDim] = 0.5*(Coord_2[iDim]-Coord_1[iDim]);
          b[iDim] = 0.5*(Coord_3[iDim]-Coord_1[iDim]);
          c[iDim] = Coord_5[iDim]-Coord_1[iDim]; }
        n[0] = a[1]*b[2]-b[1]*a[2];
        n[1] = -(a[0]*b[2]-b[0]*a[2]);
        n[2] = a[0]*b[1]-b[0]*a[1];
        test_1 = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];
        
        for(iDim = 0; iDim < nDim; iDim++) {
          a[iDim] = 0.5*(Coord_3[iDim]-Coord_2[iDim]);
          b[iDim] = 0.5*(Coord_4[iDim]-Coord_2[iDim]);
          c[iDim] = Coord_5[iDim]-Coord_2[iDim]; }
        n[0] = a[1]*b[2]-b[1]*a[2];
        n[1] = -(a[0]*b[2]-b[0]*a[2]);
        n[2] = a[0]*b[1]-b[0]*a[1];
        test_2 = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];
        
        for(iDim = 0; iDim < nDim; iDim++) {
          a[iDim] = 0.5*(Coord_4[iDim]-Coord_3[iDim]);
          b[iDim] = 0.5*(Coord_1[iDim]-Coord_3[iDim]);
          c[iDim] = Coord_5[iDim]-Coord_3[iDim]; }
        n[0] = a[1]*b[2]-b[1]*a[2];
        n[1] = -(a[0]*b[2]-b[0]*a[2]);
        n[2] = a[0]*b[1]-b[0]*a[1];
        test_3 = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];
        
        for(iDim = 0; iDim < nDim; iDim++) {
          a[iDim] = 0.5*(Coord_1[iDim]-Coord_4[iDim]);
          b[iDim] = 0.5*(Coord_3[iDim]-Coord_4[iDim]);
          c[iDim] = Coord_5[iDim]-Coord_4[iDim]; }
        n[0] = a[1]*b[2]-b[1]*a[2];
        n[1] = -(a[0]*b[2]-b[0]*a[2]);
        n[2] = a[0]*b[1]-b[0]*a[1];
        test_4 = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];
        
        if ((test_1 < 0.0) && (test_2 < 0.0) && (test_3 < 0.0) && (test_4 < 0.0))
          bound[iMarker][iElem_Surface]->Change_Orientation();
      }
    }
}

void CPhysicalGeometry::SetEsuP(void) {
  unsigned long iPoint, iElem;
  unsigned short iNode;
  
  /*--- Loop over all the elements ---*/
  for(iElem = 0; iElem < nElem; iElem++)
  /*--- Loop over all the nodes of an element ---*/
    for(iNode = 0; iNode < elem[iElem]->GetnNodes(); iNode++) {
      iPoint = elem[iElem]->GetNode(iNode);
      /*--- Store the element into the point ---*/
      node[iPoint]->SetElem(iElem);
    }
}

void CPhysicalGeometry::ComputeWall_Distance(CConfig *config) {
  
  double *coord, dist2, dist;
  unsigned short iDim, iMarker;
  unsigned long iPoint, iVertex, nVertex_SolidWall;
  
  cout << "Computing wall distances." << endl;
  
  /*--- Compute the total number of nodes on no-slip boundaries ---*/
  
  nVertex_SolidWall = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if ((config->GetMarker_All_Boundary(iMarker) == HEAT_FLUX) ||
        (config->GetMarker_All_Boundary(iMarker) == ISOTHERMAL))
      nVertex_SolidWall += GetnVertex(iMarker);
  
  /*--- Allocate an array to hold boundary node coordinates ---*/
  
  double **Coord_bound;
  Coord_bound = new double* [nVertex_SolidWall];
  for (iVertex = 0; iVertex < nVertex_SolidWall; iVertex++)
    Coord_bound[iVertex] = new double [nDim];
  
  /*--- Retrieve and store the coordinates of the no-slip boundary nodes ---*/
  
  nVertex_SolidWall = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if ((config->GetMarker_All_Boundary(iMarker) == HEAT_FLUX) ||
        (config->GetMarker_All_Boundary(iMarker) == ISOTHERMAL))
      for (iVertex = 0; iVertex < GetnVertex(iMarker); iVertex++) {
        iPoint = vertex[iMarker][iVertex]->GetNode();
        for (iDim = 0; iDim < nDim; iDim++)
          Coord_bound[nVertex_SolidWall][iDim] = node[iPoint]->GetCoord(iDim);
        nVertex_SolidWall++;
      }
  }
  
  /*--- Loop over all interior mesh nodes and compute the distances to each
   of the no-slip boundary nodes. Store the minimum distance to the wall for
   each interior mesh node. ---*/
  
  for (iPoint = 0; iPoint < GetnPoint(); iPoint++) {
    coord = node[iPoint]->GetCoord();
    dist = 1E20;
    for (iVertex = 0; iVertex < nVertex_SolidWall; iVertex++) {
      dist2 = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        dist2 += (coord[iDim]-Coord_bound[iVertex][iDim])
        *(coord[iDim]-Coord_bound[iVertex][iDim]);
      if (dist2 < dist) dist = dist2;
    }
    node[iPoint]->SetWall_Distance(sqrt(dist));
  }
  
  /*--- Deallocate the vector of boundary coordinates. ---*/
  
  for (iVertex = 0; iVertex < nVertex_SolidWall; iVertex++)
    delete[] Coord_bound[iVertex];
  delete[] Coord_bound;
  
}

void CPhysicalGeometry::ComputeReference_Area(CConfig *config) {
  unsigned short iMarker, Boundary, Monitoring;
  unsigned long iVertex, iPoint;
  double *Normal, PositiveZArea;
  
  PositiveZArea = 0.0;
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    Boundary = config->GetMarker_All_Boundary(iMarker);
    Monitoring = config->GetMarker_All_Monitoring(iMarker);
    
    if (((Boundary == EULER_WALL) || (Boundary == HEAT_FLUX) ||
         (Boundary == ISOTHERMAL)) && (Monitoring == YES))
      for(iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
        iPoint = vertex[iMarker][iVertex]->GetNode();
        if (node[iPoint]->GetDomain()) {
          Normal = vertex[iMarker][iVertex]->GetNormal();
          if (Normal[nDim-1] < 0) PositiveZArea -= Normal[nDim-1];
        }
      }
  }
  
  if (config->GetRefAreaCoeff() == 0.0)
    config->SetRefAreaCoeff(PositiveZArea);
  
  if (nDim == 2) cout << "Area projection in the y-plane = "<< PositiveZArea << "." << endl;
  else cout << "Area projection in the z-plane = "<< PositiveZArea << "." << endl;
  
}

void CPhysicalGeometry::SetPsuP(void) {
  
  unsigned short Node_Neighbor, iNode, iNeighbor;
  unsigned long jElem, Point_Neighbor, iPoint, iElem;
  
  /*--- Loop over all the points ---*/
  
  for(iPoint = 0; iPoint < nPoint; iPoint++)
    
  /*--- Loop over all elements shared by the point ---*/
    
    for(iElem = 0; iElem < node[iPoint]->GetnElem(); iElem++) {
      
      jElem = node[iPoint]->GetElem(iElem);
      
      /*--- If we find the point iPoint in the surronding element ---*/
      
      for(iNode = 0; iNode < elem[jElem]->GetnNodes(); iNode++)
        
        if (elem[jElem]->GetNode(iNode) == iPoint)
          
        /*--- Localize the local index of the neighbor of iPoint in the element ---*/
          
          for(iNeighbor = 0; iNeighbor < elem[jElem]->GetnNeighbor_Nodes(iNode); iNeighbor++) {
            Node_Neighbor = elem[jElem]->GetNeighbor_Nodes(iNode,iNeighbor);
            Point_Neighbor = elem[jElem]->GetNode(Node_Neighbor);
            
            /*--- Store the point into the point ---*/
            
            node[iPoint]->SetPoint(Point_Neighbor);
          }
    }
  
  /*--- Set the number of neighbors variable, this is
   important for JST and multigrid in parallel ---*/
  
  for(iPoint = 0; iPoint < nPoint; iPoint++)
    node[iPoint]->SetnNeighbor(node[iPoint]->GetnPoint());
  
}

void CPhysicalGeometry::SetRCM(CConfig *config) {
  unsigned long iPoint, AdjPoint, AuxPoint, AddPoint, iElem;
  vector<unsigned long> Queue, AuxQueue, Result;
  unsigned short Degree, MinDegree, iNode, jNode, iDim, iMarker;
  bool *AddIndex, inQueue;
  
  AddIndex = new bool [nPoint];
  
  for(iPoint = 0; iPoint < nPoint; iPoint++)
    AddIndex[iPoint] = false;
  
  /*--- Select the node with the lowest degree in the grid. ---*/
  
  MinDegree = node[0]->GetnNeighbor(); AddPoint = 0;
  for(iPoint = 1; iPoint < nPointDomain; iPoint++) {
    Degree = node[iPoint]->GetnPoint();
    if ((Degree < MinDegree) && (AddIndex[iPoint] == false)) {
      MinDegree = Degree;
      AddPoint = iPoint;
    }
  }
  
  /*--- Add the node in the first free position. ---*/
  Result.push_back(AddPoint); AddIndex[AddPoint] = true;
  
  /*--- Loop until reorganize all the nodes ---*/
  
  while (Result.size() < nPointDomain) {
    
    /*--- Add to the queue all the nodes adjacent in the increasing
     order of their degree, checking if the element is already
     in the Queue. ---*/
    
    AuxQueue.clear();
    for (iNode = 0; iNode < node[AddPoint]->GetnPoint(); iNode++) {
      AdjPoint = node[AddPoint]->GetPoint(iNode);
      
      inQueue = false;
      for (jNode = 0; jNode < Queue.size(); jNode++) {
        if (Queue[jNode] == AdjPoint) { inQueue = true; break; }
      }
      
      if ((AddIndex[AdjPoint] == false) && (!inQueue) && (AdjPoint < nPointDomain))
        
        AuxQueue.push_back(AdjPoint);
    }
    
    for (iNode = 0; iNode < AuxQueue.size(); iNode++) {
      for (jNode = 0; jNode < AuxQueue.size() - 1 - iNode; jNode++) {
        if (node[AuxQueue[jNode]]->GetnPoint() > node[AuxQueue[jNode+1]]->GetnPoint()) {
          AuxPoint = AuxQueue[jNode];
          AuxQueue[jNode] = AuxQueue[jNode+1];
          AuxQueue[jNode+1] = AuxPoint;
        }
      }
    }
    
    Queue.insert(Queue.end(), AuxQueue.begin(), AuxQueue.end());
    
    /*--- Extract the first node from the queue and add it in the first free
     position. ---*/
    
    AddPoint = Queue[0];
    Result.push_back(Queue[0]); AddIndex[Queue[0]] = true;
    Queue.erase (Queue.begin(), Queue.begin()+1);
    
    /*--- Add to the queue all the nodes adjacent in the increasing
     order of their degree, checking if the element is already
     in the Queue. ---*/
    
  }
  
  delete[] AddIndex;
  
  reverse(Result.begin(), Result.end());
  
  /*--- Add the MPI points ---*/
  for(iPoint = nPointDomain; iPoint < nPoint; iPoint++) {
    Result.push_back(iPoint);
  }
  
  /*--- Reset old data structures ---*/
  
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    node[iPoint]->ResetElem();
    node[iPoint]->ResetPoint();
    node[iPoint]->ResetBoundary();
    node[iPoint]->SetPhysicalBoundary(false);
  }
  
  /*--- Set the new coordinates ---*/
  
  double **AuxCoord;
  
  AuxCoord = new double* [nPoint];
  for (iPoint = 0; iPoint < nPoint; iPoint++)
    AuxCoord[iPoint] = new double [nDim];
  
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    for (iDim = 0; iDim < nDim; iDim++) {
      AuxCoord[iPoint][iDim] = node[iPoint]->GetCoord(iDim);
    }
  }
  
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    for (iDim = 0; iDim < nDim; iDim++)
      node[iPoint]->SetCoord(iDim, AuxCoord[Result[iPoint]][iDim]);
  }
  
  for (iPoint = 0; iPoint < nPoint; iPoint++)
    delete[] AuxCoord[iPoint];
  delete[] AuxCoord;
  
  /*--- Set the new conectivities ---*/
  
  unsigned long *InvResult;
  InvResult = new unsigned long [nPoint];
  for (iPoint = 0; iPoint < nPoint; iPoint++)
    InvResult[Result[iPoint]] = iPoint;
  
  
  for(iElem = 0; iElem < nElem; iElem++) {
    
    for (iNode = 0; iNode < elem[iElem]->GetnNodes(); iNode++) {
      iPoint = elem[iElem]->GetNode(iNode);
      elem[iElem]->SetNode(iNode, InvResult[iPoint]);
    }
  }
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    for (iElem = 0; iElem < nElem_Bound[iMarker]; iElem++) {
      for (iNode = 0; iNode < bound[iMarker][iElem]->GetnNodes(); iNode++) {
        iPoint = bound[iMarker][iElem]->GetNode(iNode);
        bound[iMarker][iElem]->SetNode(iNode, InvResult[iPoint]);
        node[InvResult[iPoint]]->SetBoundary(nMarker);
        if (config->GetMarker_All_Boundary(iMarker) != SEND_RECEIVE &&
            config->GetMarker_All_Boundary(iMarker) != INTERFACE_BOUNDARY &&
            config->GetMarker_All_Boundary(iMarker) != NEARFIELD_BOUNDARY &&
            config->GetMarker_All_Boundary(iMarker) != PERIODIC_BOUNDARY)
          node[InvResult[iPoint]]->SetPhysicalBoundary(true);
      }
    }
  }
  
  delete[] InvResult;
  
}

void CPhysicalGeometry::SetEsuE(void) {
  unsigned short first_elem_face, second_elem_face, iFace, iNode, jElem;
  unsigned long face_point, Test_Elem, iElem;
  
  /*--- Loop over all the elements, faces and nodes ---*/
  for(iElem = 0; iElem < nElem; iElem++)
    for (iFace = 0; iFace < elem[iElem]->GetnFaces(); iFace++)
      for (iNode = 0; iNode < elem[iElem]->GetnNodesFace(iFace); iNode++) {
        face_point = elem[iElem]->GetNode(elem[iElem]->GetFaces(iFace,iNode));
        /*--- Loop over all elements sharing the face point ---*/
        for(jElem = 0; jElem < node[face_point]->GetnElem(); jElem++) {
          Test_Elem = node[face_point]->GetElem(jElem);
          /*--- If it is a new element in this face ---*/
          if ((elem[iElem]->GetNeighbor_Elements(iFace) == -1) && (iElem < Test_Elem) &&
              (FindFace(iElem, Test_Elem, first_elem_face, second_elem_face))) {
            /*--- Localice which faces are sharing both elements ---*/
            elem[iElem]->SetNeighbor_Elements(Test_Elem,first_elem_face);
            /*--- Store the element for both elements ---*/
            elem[Test_Elem]->SetNeighbor_Elements(iElem,second_elem_face);
            
          }
        }
      }
}

void CPhysicalGeometry::SetBoundVolume(void) {
  unsigned short cont, iMarker, iElem, iNode_Domain, iNode_Surface;
  unsigned long Point_Domain, Point_Surface, Point, iElem_Surface, iElem_Domain;
  bool CheckVol;
  
  for (iMarker = 0; iMarker < nMarker; iMarker++)
    for (iElem_Surface = 0; iElem_Surface < nElem_Bound[iMarker]; iElem_Surface++) {
      
      /*--- Choose and arbitrary point from the surface --*/
      Point = bound[iMarker][iElem_Surface]->GetNode(0);
      CheckVol = false;
      
      for (iElem = 0; iElem < node[Point]->GetnElem(); iElem++) {
        /*--- Look for elements surronding that point --*/
        cont = 0; iElem_Domain = node[Point]->GetElem(iElem);
        for (iNode_Domain = 0; iNode_Domain < elem[iElem_Domain]->GetnNodes(); iNode_Domain++) {
          Point_Domain = elem[iElem_Domain]->GetNode(iNode_Domain);
          for (iNode_Surface = 0; iNode_Surface < bound[iMarker][iElem_Surface]->GetnNodes(); iNode_Surface++) {
            Point_Surface = bound[iMarker][iElem_Surface]->GetNode(iNode_Surface);
            if (Point_Surface == Point_Domain) cont++;
            if (cont == bound[iMarker][iElem_Surface]->GetnNodes()) break;
          }
          if (cont == bound[iMarker][iElem_Surface]->GetnNodes()) break;
        }
        
        if (cont == bound[iMarker][iElem_Surface]->GetnNodes()) {
          bound[iMarker][iElem_Surface]->SetDomainElement(iElem_Domain);
          CheckVol = true;
          break;
        }
      }
      if (!CheckVol) {
        cout << "The surface element ("<< iMarker <<", "<< iElem_Surface << ") doesn't have an associated volume element." << endl;
        exit(1);
      }
    }
}

void CPhysicalGeometry::SetVertex(CConfig *config) {
  unsigned long  iPoint, iVertex, iElem;
  unsigned short iMarker, iNode;
  
  /*--- Initialize the Vertex vector for each node of the grid ---*/
  for (iPoint = 0; iPoint < nPoint; iPoint++)
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      node[iPoint]->SetVertex(-1,iMarker);
  
  /*--- Create and compute the vector with the number of vertex per marker ---*/
  nVertex = new unsigned long [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    /*--- Initialize the number of Bound Vertex for each Marker ---*/
    nVertex[iMarker] = 0;
    for (iElem = 0; iElem < nElem_Bound[iMarker]; iElem++)
      for(iNode = 0; iNode < bound[iMarker][iElem]->GetnNodes(); iNode++) {
        iPoint = bound[iMarker][iElem]->GetNode(iNode);
        /*--- Set the vertex in the node information ---*/
        if ((node[iPoint]->GetVertex(iMarker) == -1) || (config->GetMarker_All_Boundary(iMarker) == SEND_RECEIVE)) {
          iVertex = nVertex[iMarker];
          node[iPoint]->SetVertex(nVertex[iMarker],iMarker);
          nVertex[iMarker]++;
        }
      }
  }
  
  /*--- Initialize the Vertex vector for each node, the previous result is deleted ---*/
  for (iPoint = 0; iPoint < nPoint; iPoint++)
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      node[iPoint]->SetVertex(-1,iMarker);
  
  /*--- Create the bound vertex structure, note that the order
   is the same as in the input file, this is important for Send/Receive part ---*/
  vertex = new CVertex**[nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    vertex[iMarker] = new CVertex* [nVertex[iMarker]];
    nVertex[iMarker] = 0;
    
    /*--- Initialize the number of Bound Vertex for each Marker ---*/
    for (iElem = 0; iElem < nElem_Bound[iMarker]; iElem++)
      for(iNode = 0; iNode < bound[iMarker][iElem]->GetnNodes(); iNode++) {
        iPoint = bound[iMarker][iElem]->GetNode(iNode);
        /*--- Set the vertex in the node information ---*/
        if ((node[iPoint]->GetVertex(iMarker) == -1) || (config->GetMarker_All_Boundary(iMarker) == SEND_RECEIVE)){
          iVertex = nVertex[iMarker];
          vertex[iMarker][iVertex] = new CVertex(iPoint, nDim);
          
          if (config->GetMarker_All_Boundary(iMarker) == SEND_RECEIVE) {
            vertex[iMarker][iVertex]->SetRotation_Type(bound[iMarker][iElem]->GetRotation_Type());
          }
          node[iPoint]->SetVertex(nVertex[iMarker],iMarker);
          nVertex[iMarker]++;
        }
      }
  }
}

void CPhysicalGeometry::SetCG(void) {
  unsigned short nNode, iDim, iMarker, iNode;
  unsigned long elem_poin, edge_poin, iElem, iEdge;
  double **Coord;
  
  /*--- Compute the center of gravity for elements ---*/
  for(iElem = 0; iElem<nElem; iElem++) {
    nNode = elem[iElem]->GetnNodes();
    Coord = new double* [nNode];
    /*--- Store the coordinates for all the element nodes ---*/
    for (iNode = 0; iNode < nNode; iNode++) {
      elem_poin = elem[iElem]->GetNode(iNode);
      Coord[iNode] = new double [nDim];
      for (iDim = 0; iDim < nDim; iDim++)
        Coord[iNode][iDim]=node[elem_poin]->GetCoord(iDim);
    }
    /*--- Compute the element CG coordinates ---*/
    elem[iElem]->SetCG(Coord);
    
    for (iNode = 0; iNode < nNode; iNode++)
      if (Coord[iNode] != NULL) delete[] Coord[iNode];
    if (Coord != NULL) delete[] Coord;
  }
  
  /*--- Center of gravity for face elements ---*/
  for(iMarker = 0; iMarker < nMarker; iMarker++)
    for(iElem = 0; iElem < nElem_Bound[iMarker]; iElem++) {
      nNode = bound[iMarker][iElem]->GetnNodes();
      Coord = new double* [nNode];
      /*--- Store the coordinates for all the element nodes ---*/
      for (iNode = 0; iNode < nNode; iNode++) {
        elem_poin = bound[iMarker][iElem]->GetNode(iNode);
        Coord[iNode] = new double [nDim];
        for (iDim = 0; iDim < nDim; iDim++)
          Coord[iNode][iDim]=node[elem_poin]->GetCoord(iDim);
      }
      /*--- Compute the element CG coordinates ---*/
      bound[iMarker][iElem]->SetCG(Coord);
      for (iNode=0; iNode < nNode; iNode++)
        if (Coord[iNode] != NULL) delete[] Coord[iNode];
      if (Coord != NULL) delete[] Coord;
    }
  
  /*--- Center of gravity for edges ---*/
  for (iEdge = 0; iEdge < nEdge; iEdge++) {
    nNode = edge[iEdge]->GetnNodes();
    Coord = new double* [nNode];
    /*--- Store the coordinates for all the element nodes ---*/
    for (iNode = 0; iNode < nNode; iNode++) {
      edge_poin=edge[iEdge]->GetNode(iNode);
      Coord[iNode] = new double [nDim];
      for (iDim = 0; iDim<nDim; iDim++)
        Coord[iNode][iDim]=node[edge_poin]->GetCoord(iDim);
    }
    /*--- Compute the edge CG coordinates ---*/
    edge[iEdge]->SetCG(Coord);
    
    for (iNode=0; iNode < nNode; iNode++)
      if (Coord[iNode] != NULL) delete[] Coord[iNode];
    if (Coord != NULL) delete[] Coord;
  }
}

void CPhysicalGeometry::SetBoundControlVolume(CConfig *config, unsigned short action) {
  unsigned short Neighbor_Node, iMarker, iNode, iNeighbor_Nodes, iDim;
  unsigned long Neighbor_Point, iVertex, iPoint, iElem;
  long iEdge;
  double Area, *NormalFace = NULL;
  
  /*--- Update values of faces of the edge ---*/
  if (action != ALLOCATE)
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++)
        vertex[iMarker][iVertex]->SetZeroValues();
  
  double *Coord_Edge_CG = new double [nDim];
  double *Coord_Elem_CG = new double [nDim];
  double *Coord_Vertex = new double [nDim];
  
  /*--- Loop over all the markers ---*/
  for (iMarker = 0; iMarker < nMarker; iMarker++)
  /*--- Loop over all the boundary elements ---*/
    for (iElem = 0; iElem < nElem_Bound[iMarker]; iElem++)
    /*--- Loop over all the nodes of the boundary ---*/
      for(iNode = 0; iNode < bound[iMarker][iElem]->GetnNodes(); iNode++) {
        iPoint = bound[iMarker][iElem]->GetNode(iNode);
        iVertex = node[iPoint]->GetVertex(iMarker);
        /*--- Loop over the neighbor nodes, there is a face for each one ---*/
        for(iNeighbor_Nodes = 0; iNeighbor_Nodes < bound[iMarker][iElem]->GetnNeighbor_Nodes(iNode); iNeighbor_Nodes++) {
          Neighbor_Node = bound[iMarker][iElem]->GetNeighbor_Nodes(iNode,iNeighbor_Nodes);
          Neighbor_Point = bound[iMarker][iElem]->GetNode(Neighbor_Node);
          /*--- Shared edge by the Neighbor Point and the point ---*/
          iEdge = FindEdge(iPoint, Neighbor_Point);
          for (iDim = 0; iDim < nDim; iDim++) {
            Coord_Edge_CG[iDim] = edge[iEdge]->GetCG(iDim);
            Coord_Elem_CG[iDim] = bound[iMarker][iElem]->GetCG(iDim);
            Coord_Vertex[iDim] = node[iPoint]->GetCoord(iDim);
          }
          switch (nDim) {
            case 2:
              /*--- Store the 2D face ---*/
              if (iNode == 0) vertex[iMarker][iVertex]->SetNodes_Coord(Coord_Elem_CG, Coord_Vertex, config);
              if (iNode == 1) vertex[iMarker][iVertex]->SetNodes_Coord(Coord_Vertex, Coord_Elem_CG, config);
              break;
            case 3:
              /*--- Store the 3D face ---*/
              if (iNeighbor_Nodes == 0) vertex[iMarker][iVertex]->SetNodes_Coord(Coord_Elem_CG, Coord_Edge_CG, Coord_Vertex, config);
              if (iNeighbor_Nodes == 1) vertex[iMarker][iVertex]->SetNodes_Coord(Coord_Edge_CG, Coord_Elem_CG, Coord_Vertex, config);
              break;
          }
        }
      }
  
  delete[] Coord_Edge_CG;
  delete[] Coord_Elem_CG;
  delete[] Coord_Vertex;
  
  /*--- Check if there is a normal with null area ---*/
  for (iMarker = 0; iMarker < nMarker; iMarker ++)
    for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
      NormalFace = vertex[iMarker][iVertex]->GetNormal();
      Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += NormalFace[iDim]*NormalFace[iDim];
      Area = sqrt(Area);
      if (Area == 0.0) for (iDim = 0; iDim < nDim; iDim++) NormalFace[iDim] = EPS*EPS;
    }
  
}

void CPhysicalGeometry::MatchNearField(CConfig *config) {
  double epsilon = 1e-1;
  
  unsigned short nMarker_NearfieldBound = config->GetnMarker_NearFieldBound();
  
  if (nMarker_NearfieldBound != 0) {
    
    unsigned short iMarker, jMarker;
    unsigned long iVertex, iPoint, jVertex, jPoint = 0, pPoint = 0;
    double *Coord_i, *Coord_j, dist = 0.0, mindist, maxdist;
    
    cout << "Set Near-Field boundary conditions. " <<endl;
    
    maxdist = 0.0;
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
      if (config->GetMarker_All_Boundary(iMarker) == NEARFIELD_BOUNDARY) {
        for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
          iPoint = vertex[iMarker][iVertex]->GetNode();
          Coord_i = node[iPoint]->GetCoord();
          
          mindist = 1e10;
          for (jMarker = 0; jMarker < config->GetnMarker_All(); jMarker++)
            if ((config->GetMarker_All_Boundary(jMarker) == NEARFIELD_BOUNDARY) && (iMarker != jMarker))
              for (jVertex = 0; jVertex < nVertex[jMarker]; jVertex++) {
                jPoint = vertex[jMarker][jVertex]->GetNode();
                Coord_j = node[jPoint]->GetCoord();
                if (nDim == 2) dist = sqrt(pow(Coord_j[0]-Coord_i[0],2.0) + pow(Coord_j[1]-Coord_i[1],2.0));
                if (nDim == 3) dist = sqrt(pow(Coord_j[0]-Coord_i[0],2.0) + pow(Coord_j[1]-Coord_i[1],2.0) + pow(Coord_j[2]-Coord_i[2],2.0));
                if (dist < mindist) { mindist = dist; pPoint = jPoint; }
              }
          maxdist = max(maxdist, mindist);
          vertex[iMarker][iVertex]->SetDonorPoint(pPoint);
          
          if (mindist > epsilon) {
            cout.precision(10);
            cout << endl;
            cout << "   Bad match for point " << iPoint << ".\tNearest";
            cout << " donor distance: " << scientific << mindist << ".";
            vertex[iMarker][iVertex]->SetDonorPoint(iPoint);
            maxdist = min(maxdist, 0.0);
          }
        }
        cout <<"The max distance between points is: " << maxdist <<"."<< endl;
      }
    
  }
  
}

void CPhysicalGeometry::MatchInterface(CConfig *config) {
  double epsilon = 1.5e-1;
  
  unsigned short nMarker_InterfaceBound = config->GetnMarker_InterfaceBound();
  
  if (nMarker_InterfaceBound != 0) {
    
    unsigned short iMarker, jMarker;
    unsigned long iVertex, iPoint, jVertex, jPoint = 0, pPoint = 0;
    double *Coord_i, *Coord_j, dist = 0.0, mindist, maxdist;
    
    cout << "Set Interface boundary conditions." << endl;
    
    maxdist = 0.0;
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
      if (config->GetMarker_All_Boundary(iMarker) == INTERFACE_BOUNDARY) {
        for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
          iPoint = vertex[iMarker][iVertex]->GetNode();
          Coord_i = node[iPoint]->GetCoord();
          
          mindist = 1E6;
          for (jMarker = 0; jMarker < config->GetnMarker_All(); jMarker++)
            if ((config->GetMarker_All_Boundary(jMarker) == INTERFACE_BOUNDARY) && (iMarker != jMarker))
              for (jVertex = 0; jVertex < nVertex[jMarker]; jVertex++) {
                jPoint = vertex[jMarker][jVertex]->GetNode();
                Coord_j = node[jPoint]->GetCoord();
                if (nDim == 2) dist = sqrt(pow(Coord_j[0]-Coord_i[0],2.0) + pow(Coord_j[1]-Coord_i[1],2.0));
                if (nDim == 3) dist = sqrt(pow(Coord_j[0]-Coord_i[0],2.0) + pow(Coord_j[1]-Coord_i[1],2.0) + pow(Coord_j[2]-Coord_i[2],2.0));
                if (dist < mindist) {mindist = dist; pPoint = jPoint;}
              }
          maxdist = max(maxdist, mindist);
          vertex[iMarker][iVertex]->SetDonorPoint(pPoint);
          
          if (mindist > epsilon) {
            cout.precision(10);
            cout << endl;
            cout << "   Bad match for point " << iPoint << ".\tNearest";
            cout << " donor distance: " << scientific << mindist << ".";
            vertex[iMarker][iVertex]->SetDonorPoint(iPoint);
            maxdist = min(maxdist, 0.0);
          }
          
        }
        cout <<"The max distance between points is: " << maxdist <<"."<< endl;
      }
  }
}


void CPhysicalGeometry::SetControlVolume(CConfig *config, unsigned short action) {
  unsigned long face_iPoint = 0, face_jPoint = 0, iPoint, iElem;
  long iEdge;
  unsigned short nEdgesFace = 1, iFace, iEdgesFace, iDim;
  double *Coord_Edge_CG, *Coord_FaceElem_CG, *Coord_Elem_CG, *Coord_FaceiPoint, *Coord_FacejPoint, Area,
  Volume, DomainVolume, my_DomainVolume, *NormalFace = NULL;
  bool change_face_orientation;
  int rank;
  
  rank = MASTER_NODE;
  
  /*--- Update values of faces of the edge ---*/
  if (action != ALLOCATE) {
    for(iEdge = 0; iEdge < nEdge; iEdge++)
      edge[iEdge]->SetZeroValues();
    for(iPoint = 0; iPoint < nPoint; iPoint++)
      node[iPoint]->SetVolume (0.0);
  }
  
  Coord_Edge_CG = new double [nDim];
  Coord_FaceElem_CG = new double [nDim];
  Coord_Elem_CG = new double [nDim];
  Coord_FaceiPoint = new double [nDim];
  Coord_FacejPoint = new double [nDim];
  
  my_DomainVolume = 0.0;
  for(iElem = 0; iElem < nElem; iElem++)
    for (iFace = 0; iFace < elem[iElem]->GetnFaces(); iFace++) {
      
      /*--- In 2D all the faces have only one edge ---*/
      if (nDim == 2) nEdgesFace = 1;
      /*--- In 3D the number of edges per face is the same as the number of point per face ---*/
      if (nDim == 3) nEdgesFace = elem[iElem]->GetnNodesFace(iFace);
      
      /*-- Loop over the edges of a face ---*/
      for (iEdgesFace = 0; iEdgesFace < nEdgesFace; iEdgesFace++) {
        
        /*--- In 2D only one edge (two points) per edge ---*/
        if (nDim == 2) {
          face_iPoint = elem[iElem]->GetNode(elem[iElem]->GetFaces(iFace,0));
          face_jPoint = elem[iElem]->GetNode(elem[iElem]->GetFaces(iFace,1));
        }
        
        /*--- In 3D there are several edges in each face ---*/
        if (nDim == 3) {
          face_iPoint = elem[iElem]->GetNode(elem[iElem]->GetFaces(iFace,iEdgesFace));
          if (iEdgesFace != nEdgesFace-1)
            face_jPoint = elem[iElem]->GetNode(elem[iElem]->GetFaces(iFace,iEdgesFace+1));
          else
            face_jPoint = elem[iElem]->GetNode(elem[iElem]->GetFaces(iFace,0));
        }
        
        /*--- We define a direction (from the smalest index to the greatest) --*/
        change_face_orientation = false;
        if (face_iPoint > face_jPoint) change_face_orientation = true;
        iEdge = FindEdge(face_iPoint, face_jPoint);
        
        for (iDim = 0; iDim < nDim; iDim++) {
          Coord_Edge_CG[iDim] = edge[iEdge]->GetCG(iDim);
          Coord_Elem_CG[iDim] = elem[iElem]->GetCG(iDim);
          Coord_FaceElem_CG[iDim] = elem[iElem]->GetFaceCG(iFace,iDim);
          Coord_FaceiPoint[iDim] = node[face_iPoint]->GetCoord(iDim);
          Coord_FacejPoint[iDim] = node[face_jPoint]->GetCoord(iDim);
        }
        
        switch (nDim) {
          case 2:
            /*--- Two dimensional problem ---*/
            if (change_face_orientation) edge[iEdge]->SetNodes_Coord(Coord_Elem_CG, Coord_Edge_CG, config);
            else edge[iEdge]->SetNodes_Coord(Coord_Edge_CG,Coord_Elem_CG, config);
            Area = edge[iEdge]->GetVolume(Coord_FaceiPoint,Coord_Edge_CG,Coord_Elem_CG);
            node[face_iPoint]->AddVolume(Area); my_DomainVolume +=Area;
            Area = edge[iEdge]->GetVolume(Coord_FacejPoint,Coord_Edge_CG,Coord_Elem_CG);
            node[face_jPoint]->AddVolume(Area); my_DomainVolume +=Area;
            break;
          case 3:
            /*--- Three dimensional problem ---*/
            if (change_face_orientation) edge[iEdge]->SetNodes_Coord(Coord_FaceElem_CG,Coord_Edge_CG,Coord_Elem_CG, config);
            else edge[iEdge]->SetNodes_Coord(Coord_Edge_CG,Coord_FaceElem_CG,Coord_Elem_CG, config);
            Volume = edge[iEdge]->GetVolume(Coord_FaceiPoint,Coord_Edge_CG,Coord_FaceElem_CG, Coord_Elem_CG);
            node[face_iPoint]->AddVolume(Volume); my_DomainVolume +=Volume;
            Volume = edge[iEdge]->GetVolume(Coord_FacejPoint,Coord_Edge_CG,Coord_FaceElem_CG, Coord_Elem_CG);
            node[face_jPoint]->AddVolume(Volume); my_DomainVolume +=Volume;
            break;
        }
      }
    }
  
  /*--- Check if there is a normal with null area ---*/
  for (iEdge = 0; iEdge < nEdge; iEdge++) {
    NormalFace = edge[iEdge]->GetNormal();
    Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += NormalFace[iDim]*NormalFace[iDim];
    Area = sqrt(Area);
    if (Area == 0.0) for (iDim = 0; iDim < nDim; iDim++) NormalFace[iDim] = EPS*EPS;
  }
  
  DomainVolume = my_DomainVolume;
  
  if (action == ALLOCATE) {
    if (nDim == 2) cout <<"Area of the computational grid: "<< DomainVolume <<"."<<endl;
    if (nDim == 3) cout <<"Volume of the computational grid: "<< DomainVolume <<"."<<endl;
  }
  
  config->SetDomainVolume(DomainVolume);
  
  delete[] Coord_Edge_CG;
  delete[] Coord_FaceElem_CG;
  delete[] Coord_Elem_CG;
  delete[] Coord_FaceiPoint;
  delete[] Coord_FacejPoint;
}

void CPhysicalGeometry::SetMeshFile (CConfig *config, string val_mesh_out_filename) {
  unsigned long iElem, iPoint, iElem_Bound;
  unsigned short iMarker, iNodes, iDim;
  ofstream output_file;
  string Grid_Marker;
  char *cstr;
  
  cstr = new char [val_mesh_out_filename.size()+1];
  strcpy (cstr, val_mesh_out_filename.c_str());
  
  /*--- Open .su2 grid file ---*/
  output_file.precision(15);
  output_file.open(cstr, ios::out);
  
  /*--- Write dimension, number of elements and number of points ---*/
  output_file << "NDIME= " << nDim << endl;
  output_file << "NELEM= " << nElem << endl;
  for (iElem = 0; iElem < nElem; iElem++) {
    output_file << elem[iElem]->GetVTK_Type();
    for (iNodes = 0; iNodes < elem[iElem]->GetnNodes(); iNodes++)
      output_file << "\t" << elem[iElem]->GetNode(iNodes);
    output_file << "\t"<<iElem<<endl;
  }
  
  /*--- Write the node coordinates ---*/
  output_file << "NPOIN= " << nPoint << "\t" << nPointDomain << endl;
  output_file.precision(15);
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    for (iDim = 0; iDim < nDim; iDim++)
      output_file << scientific << "\t" << node[iPoint]->GetCoord(iDim) ;
    output_file << "\t" << iPoint << endl;
    
  }
  
  /*--- Loop through and write the boundary info ---*/
  output_file << "NMARK= " << nMarker << endl;
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
    /*--- Ignore SEND_RECEIVE for the moment ---*/
    if (bound[iMarker][0]->GetVTK_Type() != VERTEX) {
      
      Grid_Marker = config->GetMarker_All_Tag(iMarker);
      output_file << "MARKER_TAG= " << Grid_Marker <<endl;
      output_file << "MARKER_ELEMS= " << nElem_Bound[iMarker]<< endl;
      
      if (nDim == 2) {
        for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
          output_file << bound[iMarker][iElem_Bound]->GetVTK_Type() << "\t" ;
          for (iNodes = 0; iNodes < bound[iMarker][iElem_Bound]->GetnNodes(); iNodes++)
            output_file << bound[iMarker][iElem_Bound]->GetNode(iNodes) << "\t" ;
          output_file	<< iElem_Bound << endl;
        }
      }
      
      if (nDim == 3) {
        for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
          output_file << bound[iMarker][iElem_Bound]->GetVTK_Type() << "\t" ;
          for (iNodes = 0; iNodes < bound[iMarker][iElem_Bound]->GetnNodes(); iNodes++)
            output_file << bound[iMarker][iElem_Bound]->GetNode(iNodes) << "\t" ;
          output_file	<< iElem_Bound << endl;
        }
      }
      
    } else if (bound[iMarker][0]->GetVTK_Type() == VERTEX) {
      output_file << "MARKER_TAG= SEND_RECEIVE" << endl;
      output_file << "MARKER_ELEMS= " << nElem_Bound[iMarker]<< endl;
      if (config->GetMarker_All_SendRecv(iMarker) > 0) output_file << "SEND_TO= " << config->GetMarker_All_SendRecv(iMarker) << endl;
      if (config->GetMarker_All_SendRecv(iMarker) < 0) output_file << "SEND_TO= " << config->GetMarker_All_SendRecv(iMarker) << endl;
      
      for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
        output_file << bound[iMarker][iElem_Bound]->GetVTK_Type() << "\t" <<
        bound[iMarker][iElem_Bound]->GetNode(0) << "\t" <<
        bound[iMarker][iElem_Bound]->GetRotation_Type() << endl;
      }
      
    }
  }
  
  output_file.close();
}

void CPhysicalGeometry::SetMeshFile(CConfig *config, string val_mesh_out_filename, string val_mesh_in_filename) {
  unsigned long iElem, iPoint, iElem_Bound, nElem_, nElem_Bound_, vnodes_edge[2], vnodes_triangle[3], vnodes_quad[4], vnodes_tetra[4], vnodes_hexa[8], vnodes_wedge[6], vnodes_pyramid[5], vnodes_vertex;
  unsigned short iMarker, iDim, iChar, VTK_Type, nDim_, nMarker_, transform;
  char *cstr;
  long SendRecv;
  ofstream output_file;
  ifstream input_file;
  string Grid_Marker, text_line, Marker_Tag;
  string::size_type position;
  
  /*--- Open output_file .su2 grid file ---*/
  cstr = new char [val_mesh_out_filename.size()+1];
  strcpy (cstr, val_mesh_out_filename.c_str());
  output_file.precision(15);
  output_file.open(cstr, ios::out);
  
  /*--- Open input_file .su2 grid file ---*/
  cstr = new char [val_mesh_in_filename.size()+1];
  strcpy (cstr, val_mesh_in_filename.c_str());
  input_file.open(cstr, ios::out);
  
  /*--- Read grid file with format SU2 ---*/
  while (getline (input_file, text_line)) {
    
    /*--- Read the dimension of the problem ---*/
    position = text_line.find ("NDIME=",0);
    if (position != string::npos) {
      text_line.erase (0,6); nDim_ = atoi(text_line.c_str());
      output_file << "NDIME= " << nDim_ << endl;
    }
    
    /*--- Read the information about inner elements ---*/
    position = text_line.find ("NELEM=",0);
    if (position != string::npos) {
      text_line.erase (0,6); nElem_ = atoi(text_line.c_str());
      output_file << "NELEM= " << nElem_ << endl;
      
      
      /*--- Loop over all the volumetric elements ---*/
      for (iElem = 0; iElem < nElem_;  iElem++) {
        getline(input_file, text_line);
        istringstream elem_line(text_line);
        
        elem_line >> VTK_Type;
        output_file << VTK_Type;
        
        switch(VTK_Type) {
          case TRIANGLE:
            elem_line >> vnodes_triangle[0]; elem_line >> vnodes_triangle[1]; elem_line >> vnodes_triangle[2];
            output_file << "\t" << vnodes_triangle[0] << "\t" << vnodes_triangle[1] << "\t" << vnodes_triangle[2] << endl;
            break;
          case RECTANGLE:
            elem_line >> vnodes_quad[0]; elem_line >> vnodes_quad[1]; elem_line >> vnodes_quad[2]; elem_line >> vnodes_quad[3];
            output_file << "\t" << vnodes_quad[0] << "\t" << vnodes_quad[1] << "\t" << vnodes_quad[2] << "\t" << vnodes_quad[3] << endl;
            break;
          case TETRAHEDRON:
            elem_line >> vnodes_tetra[0]; elem_line >> vnodes_tetra[1]; elem_line >> vnodes_tetra[2]; elem_line >> vnodes_tetra[3];
            output_file << "\t" << vnodes_tetra[0] << "\t" << vnodes_tetra[1] << "\t" << vnodes_tetra[2] << "\t" << vnodes_tetra[3] << endl;
            break;
          case HEXAHEDRON:
            elem_line >> vnodes_hexa[0]; elem_line >> vnodes_hexa[1]; elem_line >> vnodes_hexa[2];
            elem_line >> vnodes_hexa[3]; elem_line >> vnodes_hexa[4]; elem_line >> vnodes_hexa[5];
            elem_line >> vnodes_hexa[6]; elem_line >> vnodes_hexa[7];
            output_file << "\t" << vnodes_hexa[0] << "\t" << vnodes_hexa[1] << "\t" << vnodes_hexa[2] << "\t" << vnodes_hexa[3] << "\t" << vnodes_hexa[4] << "\t" << vnodes_hexa[5] << "\t" << vnodes_hexa[6] << "\t" << vnodes_hexa[7] << endl;
            break;
          case WEDGE:
            elem_line >> vnodes_wedge[0]; elem_line >> vnodes_wedge[1]; elem_line >> vnodes_wedge[2];
            elem_line >> vnodes_wedge[3]; elem_line >> vnodes_wedge[4]; elem_line >> vnodes_wedge[5];
            output_file << "\t" << vnodes_wedge[0] << "\t" << vnodes_wedge[1] << "\t" << vnodes_wedge[2] << "\t" << vnodes_wedge[3] << "\t" << vnodes_wedge[4] << "\t" << vnodes_wedge[5] << endl;
            break;
          case PYRAMID:
            elem_line >> vnodes_pyramid[0]; elem_line >> vnodes_pyramid[1]; elem_line >> vnodes_pyramid[2];
            elem_line >> vnodes_pyramid[3]; elem_line >> vnodes_pyramid[4];
            output_file << "\t" << vnodes_pyramid[0] << "\t" << vnodes_pyramid[1] << "\t" << vnodes_pyramid[2] << "\t" << vnodes_pyramid[3] << "\t" << vnodes_pyramid[4] << endl;
            break;
        }
      }
    }
    
    /*--- Coordinates ---*/
    position = text_line.find ("NPOIN=",0);
    if (position != string::npos) {
      
      /*--- Skip the lines about the points ---*/
      for (iPoint = 0; iPoint < nPoint;  iPoint++) {
        getline(input_file, text_line);
      }
      
      /*--- Add the new coordinates ---*/
      output_file << "NPOIN= " << nPoint << "\t" << nPointDomain << endl;
      for (iPoint = 0; iPoint < nPoint; iPoint++) {
        for (iDim = 0; iDim < nDim; iDim++)
          output_file << scientific << node[iPoint]->GetCoord(iDim) << "\t";
        output_file << iPoint << endl;
      }
      
    }
    
    /*--- Write the physical boundaries ---*/
    position = text_line.find ("NMARK=",0);
    if (position != string::npos) {
      
      text_line.erase (0,6); nMarker_ = atoi(text_line.c_str());
      output_file << "NMARK= " << nMarker_ << endl;
      
      for (iMarker = 0 ; iMarker < nMarker_; iMarker++) {
        
        getline (input_file,text_line);
        text_line.erase (0,11);
        string::size_type position;
        for (iChar = 0; iChar < 20; iChar++) {
          position = text_line.find( " ", 0 );
          if(position != string::npos) text_line.erase (position,1);
          position = text_line.find( "\r", 0 );
          if(position != string::npos) text_line.erase (position,1);
          position = text_line.find( "\n", 0 );
          if(position != string::npos) text_line.erase (position,1);
        }
        Marker_Tag = text_line.c_str();
        
        /*--- Standart physical boundary ---*/
        if (Marker_Tag != "SEND_RECEIVE") {
          
          getline (input_file, text_line);
          
          text_line.erase (0,13); nElem_Bound_ = atoi(text_line.c_str());
          output_file << "MARKER_TAG= " << Marker_Tag << endl;
          output_file << "MARKER_ELEMS= " << nElem_Bound_<< endl;
          
          for (iElem_Bound = 0; iElem_Bound < nElem_Bound_; iElem_Bound++) {
            
            getline(input_file, text_line);
            istringstream bound_line(text_line);
            
            bound_line >> VTK_Type;
            output_file << VTK_Type;
            
            switch(VTK_Type) {
              case LINE:
                bound_line >> vnodes_edge[0]; bound_line >> vnodes_edge[1];
                output_file << "\t" << vnodes_edge[0] << "\t" << vnodes_edge[1] << endl;
                break;
              case TRIANGLE:
                bound_line >> vnodes_triangle[0]; bound_line >> vnodes_triangle[1]; bound_line >> vnodes_triangle[2];
                output_file << "\t" << vnodes_triangle[0] << "\t" << vnodes_triangle[1] << "\t" << vnodes_triangle[2] << endl;
                break;
              case RECTANGLE:
                bound_line >> vnodes_quad[0]; bound_line >> vnodes_quad[1]; bound_line >> vnodes_quad[2]; bound_line >> vnodes_quad[3];
                output_file << "\t" << vnodes_quad[0] << "\t" << vnodes_quad[1] << "\t" << vnodes_quad[2] << "\t" << vnodes_quad[3] << endl;
                break;
            }
          }
          
        }
        
        /*--- Send-Receive boundaries definition ---*/
        else {
          output_file << "MARKER_TAG= SEND_RECEIVE" << endl;
          getline (input_file,text_line);
          text_line.erase (0,13); nElem_Bound_ = atoi(text_line.c_str());
          output_file << "MARKER_ELEMS= " << nElem_Bound_ << endl;
          getline (input_file, text_line); text_line.erase (0,8);
          SendRecv = atoi(text_line.c_str());
          output_file << "SEND_TO= " << SendRecv << endl;
          
          for (iElem_Bound = 0; iElem_Bound < nElem_Bound_; iElem_Bound++) {
            getline(input_file,text_line);
            istringstream bound_line(text_line);
            bound_line >> VTK_Type; bound_line >> vnodes_vertex; bound_line >> transform;
            output_file << VTK_Type << "\t" << vnodes_vertex << "\t" << transform << endl;
          }
        }
        
      }
    }
  }
  
  input_file.close();
  output_file.close();
  
}

void CPhysicalGeometry::SetTecPlot(char mesh_filename[200]) {
  unsigned long iElem, iPoint;
  unsigned short iDim;
  ofstream Tecplot_File;
  
  Tecplot_File.open(mesh_filename, ios::out);
  Tecplot_File << "TITLE = \"Visualization of the volumetric grid\"" << endl;
  
  if (nDim == 2) {
    Tecplot_File << "VARIABLES = \"x\",\"y\" " << endl;
    Tecplot_File << "ZONE NODES= "<< nPoint <<", ELEMENTS= "<< nElem <<", DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL"<< endl;
  }
  if (nDim == 3) {
    Tecplot_File << "VARIABLES = \"x\",\"y\",\"z\" " << endl;
    Tecplot_File << "ZONE NODES= "<< nPoint <<", ELEMENTS= "<< nElem <<", DATAPACKING=POINT, ZONETYPE=FEBRICK"<< endl;
  }
  
  for(iPoint = 0; iPoint < nPoint; iPoint++) {
    for(iDim = 0; iDim < nDim; iDim++)
      Tecplot_File << scientific << node[iPoint]->GetCoord(iDim) << "\t";
    Tecplot_File << "\n";
  }
  
  for(iElem = 0; iElem < nElem; iElem++) {
    if (elem[iElem]->GetVTK_Type() == TRIANGLE) {
      Tecplot_File <<
      elem[iElem]->GetNode(0)+1 <<" "<< elem[iElem]->GetNode(1)+1 <<" "<<
      elem[iElem]->GetNode(2)+1 <<" "<< elem[iElem]->GetNode(2)+1 << endl;
    }
    if (elem[iElem]->GetVTK_Type() == RECTANGLE) {
      Tecplot_File <<
      elem[iElem]->GetNode(0)+1 <<" "<< elem[iElem]->GetNode(1)+1 <<" "<<
      elem[iElem]->GetNode(2)+1 <<" "<< elem[iElem]->GetNode(3)+1 << endl;
    }
    if (elem[iElem]->GetVTK_Type() == TETRAHEDRON) {
      Tecplot_File <<
      elem[iElem]->GetNode(0)+1 <<" "<< elem[iElem]->GetNode(1)+1 <<" "<<
      elem[iElem]->GetNode(2)+1 <<" "<< elem[iElem]->GetNode(2)+1 <<" "<<
      elem[iElem]->GetNode(3)+1 <<" "<< elem[iElem]->GetNode(3)+1 <<" "<<
      elem[iElem]->GetNode(3)+1 <<" "<< elem[iElem]->GetNode(3)+1 << endl;
    }
    if (elem[iElem]->GetVTK_Type() == HEXAHEDRON) {
      Tecplot_File <<
      elem[iElem]->GetNode(0)+1 <<" "<< elem[iElem]->GetNode(1)+1 <<" "<<
      elem[iElem]->GetNode(2)+1 <<" "<< elem[iElem]->GetNode(3)+1 <<" "<<
      elem[iElem]->GetNode(4)+1 <<" "<< elem[iElem]->GetNode(5)+1 <<" "<<
      elem[iElem]->GetNode(6)+1 <<" "<< elem[iElem]->GetNode(7)+1 << endl;
    }
    if (elem[iElem]->GetVTK_Type() == PYRAMID) {
      Tecplot_File <<
      elem[iElem]->GetNode(0)+1 <<" "<< elem[iElem]->GetNode(1)+1 <<" "<<
      elem[iElem]->GetNode(2)+1 <<" "<< elem[iElem]->GetNode(3)+1 <<" "<<
      elem[iElem]->GetNode(4)+1 <<" "<< elem[iElem]->GetNode(4)+1 <<" "<<
      elem[iElem]->GetNode(4)+1 <<" "<< elem[iElem]->GetNode(4)+1 << endl;
    }
    if (elem[iElem]->GetVTK_Type() == WEDGE) {
      Tecplot_File <<
      elem[iElem]->GetNode(0)+1 <<" "<< elem[iElem]->GetNode(1)+1 <<" "<<
      elem[iElem]->GetNode(1)+1 <<" "<< elem[iElem]->GetNode(2)+1 <<" "<<
      elem[iElem]->GetNode(3)+1 <<" "<< elem[iElem]->GetNode(4)+1 <<" "<<
      elem[iElem]->GetNode(4)+1 <<" "<< elem[iElem]->GetNode(5)+1 << endl;
    }
  }
  
  Tecplot_File.close();
}

void CPhysicalGeometry::SetCoord_Smoothing (unsigned short val_nSmooth, double val_smooth_coeff, CConfig *config) {
  unsigned short iSmooth, nneigh, iMarker;
  double *Coord_Old, *Coord_Sum, *Coord, *Coord_i, *Coord_j, Position_Plane = 0.0;
  unsigned long iEdge, iPoint, jPoint, iVertex;
  double eps = 1E-6;
  bool NearField = false;
  
  Coord = new double [nDim];
  
  for (iPoint = 0; iPoint < GetnPoint(); iPoint++) {
    double *Coord = node[iPoint]->GetCoord();
    node[iPoint]->SetCoord_Old(Coord);
  }
  
  /*--- Jacobi iterations ---*/
  for (iSmooth = 0; iSmooth < val_nSmooth; iSmooth++) {
    
    for (iPoint = 0; iPoint < nPoint; iPoint++)
      node[iPoint]->SetCoord_SumZero();
    
    
    /*--- Loop over Interior edges ---*/
    for(iEdge = 0; iEdge < nEdge; iEdge++) {
      iPoint = edge[iEdge]->GetNode(0);
      Coord_i = node[iPoint]->GetCoord();
      
      jPoint = edge[iEdge]->GetNode(1);
      Coord_j = node[jPoint]->GetCoord();
      
      /*--- Accumulate nearest neighbor Coord to Res_sum for each variable ---*/
      node[iPoint]->AddCoord_Sum(Coord_j);
      node[jPoint]->AddCoord_Sum(Coord_i);
      
    }
    
    /*--- Loop over all mesh points (Update Coords with averaged sum) ---*/
    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      nneigh = node[iPoint]->GetnPoint();
      Coord_Sum = node[iPoint]->GetCoord_Sum();
      Coord_Old = node[iPoint]->GetCoord_Old();
      
      if (nDim == 2) {
        Coord[0] =(Coord_Old[0] + val_smooth_coeff*Coord_Sum[0]) /(1.0 + val_smooth_coeff*double(nneigh));
        Coord[1] =(Coord_Old[1] + val_smooth_coeff*Coord_Sum[1]) /(1.0 + val_smooth_coeff*double(nneigh));
        if ((NearField) && ((Coord_Old[1] > Position_Plane-eps) && (Coord_Old[1] < Position_Plane+eps)))
          Coord[1] = Coord_Old[1];
      }
      
      if (nDim == 3) {
        Coord[0] =(Coord_Old[0] + val_smooth_coeff*Coord_Sum[0]) /(1.0 + val_smooth_coeff*double(nneigh));
        Coord[1] =(Coord_Old[1] + val_smooth_coeff*Coord_Sum[1]) /(1.0 + val_smooth_coeff*double(nneigh));
        Coord[2] =(Coord_Old[2] + val_smooth_coeff*Coord_Sum[2]) /(1.0 + val_smooth_coeff*double(nneigh));
        if ((NearField) && ((Coord_Old[2] > Position_Plane-eps) && (Coord_Old[2] < Position_Plane+eps)))
          Coord[2] = Coord_Old[2];
      }
      
      node[iPoint]->SetCoord(Coord);
    }
    
    /*--- Copy boundary values ---*/
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
        iPoint = vertex[iMarker][iVertex]->GetNode();
        Coord_Old = node[iPoint]->GetCoord_Old();
        node[iPoint]->SetCoord(Coord_Old);
      }
  }
  
  delete[] Coord;
}

bool CPhysicalGeometry::FindFace(unsigned long first_elem, unsigned long second_elem, unsigned short &face_first_elem,
                                 unsigned short &face_second_elem) {
  
  /*--- Find repeated nodes between two elements to identify the common face ---*/
  unsigned long iPoint = 0, jPoint = 0;
  unsigned short face_node, iFace, iNode, jNode, kNode, nNodesFace;
  vector<unsigned long> CommonPoints, PointFaceFirst, PointFaceSecond;
  vector<unsigned long>::iterator IterPoint;
  pair<vector <unsigned long>::iterator, vector <unsigned long>::iterator> mypair;
  bool face_first_found = false, face_second_found =false;
  
  if (first_elem == second_elem) return 0;
  
  kNode = 0;
  for (iNode = 0; iNode < elem[first_elem]->GetnNodes(); iNode++) {
    iPoint = elem[first_elem]->GetNode(iNode);
    for (jNode = 0; jNode < elem[second_elem]->GetnNodes(); jNode++) {
      jPoint = elem[second_elem]->GetNode(jNode);
      if (iPoint == jPoint) {
        CommonPoints.push_back(iPoint);
        break;
      }
    }
  }
  
  /*--- Sort point in face and check that the list is unique ---*/
  sort( CommonPoints.begin(), CommonPoints.end());
  IterPoint = unique( CommonPoints.begin(), CommonPoints.end());
  CommonPoints.resize( distance(CommonPoints.begin(), IterPoint) );
  
  /*--- Search the secuence in the first element ---*/
  for (iFace = 0; iFace < elem[first_elem]->GetnFaces(); iFace++) {
    nNodesFace = elem[first_elem]->GetnNodesFace(iFace);
    for (iNode = 0; iNode < nNodesFace; iNode++) {
      face_node = elem[first_elem]->GetFaces(iFace, iNode);
      PointFaceFirst.push_back(elem[first_elem]->GetNode(face_node));
    }
    
    /*--- Sort face_poin to perform comparison ---*/
    sort( PointFaceFirst.begin(), PointFaceFirst.end());
    
    /*--- List comparison ---*/
    mypair = mismatch (PointFaceFirst.begin(), PointFaceFirst.end(), CommonPoints.begin());
    if (mypair.first == PointFaceFirst.end()) {
      face_first_elem = iFace;
      face_first_found = true;
      break;
    }
    
    PointFaceFirst.erase (PointFaceFirst.begin(),PointFaceFirst.end());
  }
  
  /*--- Search the secuence in the second element ---*/
  for (iFace = 0; iFace < elem[second_elem]->GetnFaces(); iFace++) {
    nNodesFace = elem[second_elem]->GetnNodesFace(iFace);
    for (iNode = 0; iNode < nNodesFace; iNode++) {
      face_node = elem[second_elem]->GetFaces(iFace,iNode);
      PointFaceSecond.push_back(elem[second_elem]->GetNode(face_node));
    }
    
    /*--- Sort face_poin to perform comparison ---*/
    sort( PointFaceSecond.begin(), PointFaceSecond.end());
    
    /*--- List comparison ---*/
    mypair = mismatch (PointFaceSecond.begin(), PointFaceSecond.end(), CommonPoints.begin());
    if (mypair.first == PointFaceSecond.end()) {
      face_second_elem = iFace;
      face_second_found = true;
      break;
    }
    
    PointFaceSecond.erase (PointFaceSecond.begin(),PointFaceSecond.end());
  }
  
  if (face_first_found && face_second_found) return true;
  else return false;
  
}

void CPhysicalGeometry::SetBoundTecPlot (CConfig *config, char mesh_filename[200]) {
  ofstream Tecplot_File;
  unsigned long iPoint, Total_nElem_Bound, iElem, *PointSurface = NULL, nPointSurface = 0;
  unsigned short Coord_i, iMarker;
  
  /*--- It is important to do a renumering to don't add points
   that do not belong to the surfaces ---*/
  PointSurface = new unsigned long[nPoint];
  for (iPoint = 0; iPoint < nPoint; iPoint++)
    if (node[iPoint]->GetBoundary()) {
      PointSurface[iPoint] = nPointSurface;
      nPointSurface++;
    }
  
  /*--- Compute the total number of elements ---*/
  Total_nElem_Bound = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_Plotting(iMarker) == YES) {
      Total_nElem_Bound += nElem_Bound[iMarker];
    }
  }
  
  /*--- Open the tecplot file ---*/
  Tecplot_File.open(mesh_filename, ios::out);
  Tecplot_File << "TITLE = \"Visualization of the surface grid\"" << endl;
  
  if (Total_nElem_Bound != 0) {
    
    /*--- Write the header of the file ---*/
    if (nDim == 2) {
      Tecplot_File << "VARIABLES = \"x\",\"y\" " << endl;
      Tecplot_File << "ZONE NODES= "<< nPointSurface <<", ELEMENTS= "<< Total_nElem_Bound <<", DATAPACKING=POINT, ZONETYPE=FELINESEG"<< endl;
    }
    if (nDim == 3) {
      Tecplot_File << "VARIABLES = \"x\",\"y\",\"z\" " << endl;
      Tecplot_File << "ZONE NODES= "<< nPointSurface <<", ELEMENTS= "<< Total_nElem_Bound <<", DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL"<< endl;
    }
    
    /*--- Only write the coordiantes of the points that are on the surfaces ---*/
    if (nDim == 3) {
      for(iPoint = 0; iPoint < nPoint; iPoint++)
        if (node[iPoint]->GetBoundary()) {
          for(Coord_i = 0; Coord_i < nDim-1; Coord_i++)
            Tecplot_File << node[iPoint]->GetCoord(Coord_i) << " ";
          Tecplot_File << node[iPoint]->GetCoord(nDim-1) << "\n";
        }
    }
    else {
      for(iPoint = 0; iPoint < nPoint; iPoint++)
        if (node[iPoint]->GetBoundary()){
          for(Coord_i = 0; Coord_i < nDim; Coord_i++)
            Tecplot_File << node[iPoint]->GetCoord(Coord_i) << " ";
          Tecplot_File << "\n";
        }
    }
    
    /*--- Write the cells using the new numbering ---*/
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
      if (config->GetMarker_All_Plotting(iMarker) == YES)
        for(iElem = 0; iElem < nElem_Bound[iMarker]; iElem++) {
          if (nDim == 2) {
            Tecplot_File << PointSurface[bound[iMarker][iElem]->GetNode(0)]+1 << " "
            << PointSurface[bound[iMarker][iElem]->GetNode(1)]+1 << endl;
          }
          if (nDim == 3) {
            if (bound[iMarker][iElem]->GetnNodes() == 3) {
              Tecplot_File << PointSurface[bound[iMarker][iElem]->GetNode(0)]+1 << " "
              << PointSurface[bound[iMarker][iElem]->GetNode(1)]+1 << " "
              << PointSurface[bound[iMarker][iElem]->GetNode(2)]+1 << " "
              << PointSurface[bound[iMarker][iElem]->GetNode(2)]+1 << endl;
            }
            if (bound[iMarker][iElem]->GetnNodes() == 4) {
              Tecplot_File << PointSurface[bound[iMarker][iElem]->GetNode(0)]+1 << " "
              << PointSurface[bound[iMarker][iElem]->GetNode(1)]+1 << " "
              << PointSurface[bound[iMarker][iElem]->GetNode(2)]+1 << " "
              << PointSurface[bound[iMarker][iElem]->GetNode(3)]+1 << endl;
            }
          }
        }
  }
  else {
    /*--- No elements in the surface ---*/
    if (nDim == 2) {
      Tecplot_File << "VARIABLES = \"x\",\"y\" " << endl;
      Tecplot_File << "ZONE NODES= 1, ELEMENTS= 1, DATAPACKING=POINT, ZONETYPE=FELINESEG"<< endl;
      Tecplot_File << "0.0 0.0"<< endl;
      Tecplot_File << "1 1"<< endl;
    }
    if (nDim == 3) {
      Tecplot_File << "VARIABLES = \"x\",\"y\",\"z\" " << endl;
      Tecplot_File << "ZONE NODES= 1, ELEMENTS= 1, DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL"<< endl;
      Tecplot_File << "0.0 0.0 0.0"<< endl;
      Tecplot_File << "1 1 1 1"<< endl;
    }
  }
  
  
  /*--- Dealocate memory and close the file ---*/
  delete[] PointSurface;
  Tecplot_File.close();
}

void CPhysicalGeometry::SetBoundSTL (CConfig *config, char mesh_filename[200]) {
  ofstream STL_File;
  unsigned long this_node, iNode, nNode, iElem;
  unsigned short iDim, iMarker;
  double p[3], u[3], v[3], n[3], a;
  
  /*---	STL format:
   solid NAME
   ...
   facet normal 0.00 0.00 1.00
   outer loop
   vertex  2.00  2.00  0.00
   vertex -1.00  1.00  0.00
   vertex  0.00 -1.00  0.00
   endloop
   endfacet
   ...
   end solid
   --- */
  
  /*--- Open the STL file ---*/
  STL_File.open(mesh_filename, ios::out);
  
  /*--- Write the header of the file ---*/
  STL_File << "solid surface_mesh" << endl;
  
  /*--- Write facets of surface markers ---*/
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if (config->GetMarker_All_Plotting(iMarker) == YES)
      for(iElem = 0; iElem < nElem_Bound[iMarker]; iElem++) {
        
        /*--- number of nodes for this elemnt ---*/
        nNode = bound[iMarker][iElem]->GetnNodes();
        
        /*--- Calculate Normal Vector ---*/
        for (iDim=0; iDim<nDim; iDim++){
          p[0] = node[bound[iMarker][iElem]->GetNode(0)]      ->GetCoord(iDim);
          p[1] = node[bound[iMarker][iElem]->GetNode(1)]      ->GetCoord(iDim);
          p[2] = node[bound[iMarker][iElem]->GetNode(nNode-1)]->GetCoord(iDim);
          /*cout << p[0] <<endl;
           cout << p[1] <<endl;
           cout << p[2] <<endl;*/
          u[iDim] = p[1]-p[0];
          v[iDim] = p[2]-p[0];
        }
        
        n[0] = u[1]*v[2]-u[2]*v[1];
        n[1] = u[2]*v[0]-u[0]*v[2];
        n[2] = u[0]*v[1]-u[1]*v[0];
        a = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
        /*cout << n[0] <<endl;
         cout << n[1] <<endl;
         cout << n[2] <<endl;
         cout << a << endl;*/
        
        /*--- Print normal vector ---*/
        STL_File << "  facet normal ";
        for (iDim=0; iDim<nDim; iDim++){
          STL_File << n[iDim]/a << " ";
        }
        STL_File << endl;
        
        /*--- STL Facet Loop --*/
        STL_File << "    outer loop" << endl;
        
        /*--- Print Nodes for Facet ---*/
        for(iNode=0; iNode<nNode; iNode++) {
          this_node = bound[iMarker][iElem]->GetNode(iNode);
          STL_File << "      vertex ";
          for (iDim = 0; iDim < nDim; iDim++)
            STL_File << node[this_node]->GetCoord(iDim) << " ";
          if (nDim==2)
            STL_File << 0.0 << " ";
          STL_File <<  endl;
        }
        STL_File << "    endloop" << endl;
        STL_File << "  endfacet" << endl;
      }
  
  /*--- Done with Surface Mesh ---*/
  STL_File << "endsolid" << endl;
  
  /*--- Close the file ---*/
  STL_File.close();
}

void CPhysicalGeometry::FindNormal_Neighbor(CConfig *config) {
  double cos_max, scalar_prod, norm_vect, norm_Normal, cos_alpha, diff_coord, *Normal;
  unsigned long Point_Normal, jPoint;
  unsigned short iNeigh, iMarker, iDim;
  unsigned long iPoint, iVertex;
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    
    if (config->GetMarker_All_Boundary(iMarker) != SEND_RECEIVE &&
        config->GetMarker_All_Boundary(iMarker) != INTERFACE_BOUNDARY &&
        config->GetMarker_All_Boundary(iMarker) != NEARFIELD_BOUNDARY ) {
      
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
        
        iPoint = vertex[iMarker][iVertex]->GetNode();
        Normal = vertex[iMarker][iVertex]->GetNormal();
        
        /*--- Compute closest normal neighbor, note that the normal are oriented inwards ---*/
        Point_Normal = 0; cos_max = -1.0;
        for (iNeigh = 0; iNeigh < node[iPoint]->GetnPoint(); iNeigh++) {
          jPoint = node[iPoint]->GetPoint(iNeigh);
          scalar_prod = 0.0; norm_vect = 0.0; norm_Normal = 0.0;
          for(iDim = 0; iDim < nDim; iDim++) {
            diff_coord = node[jPoint]->GetCoord(iDim)-node[iPoint]->GetCoord(iDim);
            scalar_prod += diff_coord*Normal[iDim];
            norm_vect += diff_coord*diff_coord;
            norm_Normal += Normal[iDim]*Normal[iDim];
          }
          norm_vect = sqrt(norm_vect);
          norm_Normal = sqrt(norm_Normal);
          cos_alpha = scalar_prod/(norm_vect*norm_Normal);
          
          /*--- Get maximum cosine ---*/
          if (cos_alpha >= cos_max) {
            Point_Normal = jPoint;
            cos_max = cos_alpha;
          }
        }
        vertex[iMarker][iVertex]->SetNormal_Neighbor(Point_Normal);
      }
    }
  }
}

void CPhysicalGeometry::SetGeometryPlanes(CConfig *config) {
  
  bool loop_on;
  unsigned short iMarker = 0;
  double auxXCoord, auxYCoord, auxZCoord,	*Face_Normal = NULL, auxArea, *Xcoord = NULL, *Ycoord = NULL, *Zcoord = NULL, *FaceArea = NULL;
  unsigned long jVertex, iVertex,ixCoord, iPoint, iVertex_Wall, nVertex_Wall = 0;
  
  /*--- Compute the total number of points on the near-field ---*/
  nVertex_Wall = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if ((config->GetMarker_All_Boundary(iMarker) == HEAT_FLUX) ||
        (config->GetMarker_All_Boundary(iMarker) == ISOTHERMAL) ||
        (config->GetMarker_All_Boundary(iMarker) == EULER_WALL))
      nVertex_Wall += nVertex[iMarker];
  
  
  /*--- Create an array with all the coordinates, points, pressures, face area,
   equivalent area, and nearfield weight ---*/
  Xcoord = new double[nVertex_Wall];
  Ycoord = new double[nVertex_Wall];
  if (nDim == 3)	Zcoord = new double[nVertex_Wall];
  FaceArea = new double[nVertex_Wall];
  
  /*--- Copy the boundary information to an array ---*/
  iVertex_Wall = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if ((config->GetMarker_All_Boundary(iMarker) == HEAT_FLUX) ||
        (config->GetMarker_All_Boundary(iMarker) == ISOTHERMAL) ||
        (config->GetMarker_All_Boundary(iMarker) == EULER_WALL))
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
        iPoint = vertex[iMarker][iVertex]->GetNode();
        Xcoord[iVertex_Wall] = node[iPoint]->GetCoord(0);
        Ycoord[iVertex_Wall] = node[iPoint]->GetCoord(1);
        if (nDim==3) Zcoord[iVertex_Wall] = node[iPoint]->GetCoord(2);
        Face_Normal = vertex[iMarker][iVertex]->GetNormal();
        FaceArea[iVertex_Wall] = fabs(Face_Normal[nDim-1]);
        iVertex_Wall ++;
      }
  
  
  //vector<double> XCoordList;
  vector<double>::iterator IterXCoordList;
  
  for (iVertex = 0; iVertex < nVertex_Wall; iVertex++)
    XCoordList.push_back(Xcoord[iVertex]);
  
  sort( XCoordList.begin(), XCoordList.end());
  IterXCoordList = unique( XCoordList.begin(), XCoordList.end());
  XCoordList.resize( IterXCoordList - XCoordList.begin() );
  
  /*--- Create vectors and distribute the values among the different PhiAngle queues ---*/
  Xcoord_plane.resize(XCoordList.size());
  Ycoord_plane.resize(XCoordList.size());
  if (nDim==3) Zcoord_plane.resize(XCoordList.size());
  FaceArea_plane.resize(XCoordList.size());
  Plane_points.resize(XCoordList.size());
  
  
  double dist_ratio;
  unsigned long iCoord;
  
  /*--- Distribute the values among the different PhiAngles ---*/
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    if (node[iPoint]->GetDomain()) {
      loop_on = true;
      for (ixCoord = 0; ixCoord < XCoordList.size()-1 && loop_on; ixCoord++) {
        dist_ratio = (node[iPoint]->GetCoord(0) - XCoordList[ixCoord])/(XCoordList[ixCoord+1]- XCoordList[ixCoord]);
        if (dist_ratio >= 0 && dist_ratio <= 1.0) {
          if (dist_ratio <= 0.5) iCoord = ixCoord;
          else iCoord = ixCoord+1;
          Xcoord_plane[iCoord].push_back(node[iPoint]->GetCoord(0) );
          Ycoord_plane[iCoord].push_back(node[iPoint]->GetCoord(1) );
          if (nDim==3) Zcoord_plane[iCoord].push_back(node[iPoint]->GetCoord(2) );
          FaceArea_plane[iCoord].push_back(node[iPoint]->GetVolume());   ///// CHECK AREA CALCULATION
          Plane_points[iCoord].push_back(iPoint );
          loop_on = false;
        }
      }
    }
  }
  
  unsigned long auxPoint;
  /*--- Order the arrays in ascending values of y ---*/
  for (ixCoord = 0; ixCoord < XCoordList.size(); ixCoord++)
    for (iVertex = 0; iVertex < Xcoord_plane[ixCoord].size(); iVertex++)
      for (jVertex = 0; jVertex < Xcoord_plane[ixCoord].size() - 1 - iVertex; jVertex++)
        if (Ycoord_plane[ixCoord][jVertex] > Ycoord_plane[ixCoord][jVertex+1]) {
          auxXCoord = Xcoord_plane[ixCoord][jVertex]; Xcoord_plane[ixCoord][jVertex] = Xcoord_plane[ixCoord][jVertex+1]; Xcoord_plane[ixCoord][jVertex+1] = auxXCoord;
          auxYCoord = Ycoord_plane[ixCoord][jVertex]; Ycoord_plane[ixCoord][jVertex] = Ycoord_plane[ixCoord][jVertex+1]; Ycoord_plane[ixCoord][jVertex+1] = auxYCoord;
          auxPoint = Plane_points[ixCoord][jVertex]; Plane_points[ixCoord][jVertex] = Plane_points[ixCoord][jVertex+1]; Plane_points[ixCoord][jVertex+1] = auxPoint;
          if (nDim==3) {
            auxZCoord = Zcoord_plane[ixCoord][jVertex]; Zcoord_plane[ixCoord][jVertex] = Zcoord_plane[ixCoord][jVertex+1]; Zcoord_plane[ixCoord][jVertex+1] = auxZCoord;
          }
          auxArea = FaceArea_plane[ixCoord][jVertex]; FaceArea_plane[ixCoord][jVertex] = FaceArea_plane[ixCoord][jVertex+1]; FaceArea_plane[ixCoord][jVertex+1] = auxArea;
        }
  
  /*--- Delete structures ---*/
  delete[] Xcoord; delete[] Ycoord;
  if (nDim==3) delete[] Zcoord;
  delete[] FaceArea;
}

CMultiGridGeometry::CMultiGridGeometry(CGeometry **geometry, CConfig *config_container, unsigned short iMesh) : CGeometry() {
  
  /*--- CGeometry & CConfig pointers to the fine grid level for clarity. We may
   need access to the other zones in the mesh for zone boundaries. ---*/
  
  CGeometry *fine_grid = geometry[iMesh-1];
  CConfig *config = config_container;
  
  /*--- Local variables ---*/
  
  unsigned long iPoint, Index_CoarseCV, CVPoint, iElem, iVertex, jPoint, iteration, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector, iParent, jVertex, *Buffer_Receive_Parent = NULL, *Buffer_Send_Parent = NULL, *Buffer_Receive_Children = NULL, *Buffer_Send_Children = NULL, *Parent_Remote = NULL, *Children_Remote = NULL, *Parent_Local = NULL, *Children_Local = NULL, Local_nPointCoarse, Local_nPointFine, Global_nPointCoarse, Global_nPointFine;;
  short marker_seed;
  int send_to, receive_from, rank;
  bool agglomerate_seed = true, agglomerate_CV = true;
  unsigned short nChildren, iNode, counter, iMarker, jMarker, Marker_Boundary, priority, copy_marker[MAX_NUMBER_MARKER], MarkerS, MarkerR, *nChildren_MPI;
  vector<unsigned long> Suitable_Indirect_Neighbors, Aux_Parent;
  vector<unsigned long>::iterator it;
  
  rank = MASTER_NODE;
  
  unsigned short max_children = config->GetMaxChildren(); // Set the max number of children in the agglomeration
  FinestMGLevel = false; // Set the boolean to indicate that this is a coarse multigrid level.
  nDim = fine_grid->GetnDim(); // Write the number of dimensions of the coarse grid.
  
  /*--- Create a queue system to deo the agglomeration
   1st) More than two markers ---> Vertices (never agglomerate)
   2nd) Two markers ---> Edges (agglomerate if same BC, never agglomerate if different BC)
   3rd) One marker ---> Surface (always agglomarate)
   4th) No marker ---> Internal Volume (always agglomarate) ---*/
  
  /*--- Set a marker to indicate indirect agglomeration ---*/
  
  if (iMesh == MESH_1) {
    for (iPoint = 0; iPoint < fine_grid->GetnPoint(); iPoint ++)
      fine_grid->node[iPoint]->SetAgglomerate_Indirect(false);
    for (iElem = 0; iElem < fine_grid->GetnElem(); iElem++)
      if ((fine_grid->elem[iElem]->GetVTK_Type() == HEXAHEDRON) ||
          (fine_grid->elem[iElem]->GetVTK_Type() == RECTANGLE))
        for (iNode = 0; iNode < fine_grid->elem[iElem]->GetnNodes(); iNode++) {
          iPoint = fine_grid->elem[iElem]->GetNode(iNode);
          fine_grid->node[iPoint]->SetAgglomerate_Indirect(true);
        }
  }
  
  /*--- Create the coarse grid structure using as baseline the fine grid ---*/
  
  CMultiGridQueue MGQueue_InnerCV(fine_grid->GetnPoint());
  
  node = new CPoint*[fine_grid->GetnPoint()];
  for (iPoint = 0; iPoint < fine_grid->GetnPoint(); iPoint ++) {
    
    /*--- Create node structure ---*/
    
    node[iPoint] = new CPoint(nDim, iPoint, config);
    
    /*--- Set the indirect agglomeration to false ---*/
    
    node[iPoint]->SetAgglomerate_Indirect(false);
  }
  
  Index_CoarseCV = 0;
  
  /*--- The first step is the boundary agglomeration. ---*/
  
  for (iMarker = 0; iMarker < fine_grid->GetnMarker(); iMarker++) {
    Marker_Boundary = config->GetMarker_All_Boundary(iMarker);
    
    for (iVertex = 0; iVertex < fine_grid->GetnVertex(iMarker); iVertex++) {
      iPoint = fine_grid->vertex[iMarker][iVertex]->GetNode();
      
      /*--- If the element has not being previously agglomerated and it belongs
       to the physical domain, then the agglomeration is studied ---*/
      
      if ((fine_grid->node[iPoint]->GetAgglomerate() == false) &&
          (fine_grid->node[iPoint]->GetDomain()) &&
          (GeometricalCheck(iPoint, fine_grid, config))) {
        
        nChildren = 1;
        
        /*--- We set an index for the parent control volume ---*/
        
        fine_grid->node[iPoint]->SetParent_CV(Index_CoarseCV);
        
        /*--- We add the seed point (child) to the parent control volume ---*/
        
        node[Index_CoarseCV]->SetChildren_CV(0,iPoint);
        agglomerate_seed = true; counter = 0; marker_seed = iMarker;
        
        /*--- For a particular point in the fine grid we save all the markers
         that are in that point ---*/
        
        for (jMarker = 0; jMarker < fine_grid->GetnMarker(); jMarker ++)
          if (fine_grid->node[iPoint]->GetVertex(jMarker) != -1) {
            copy_marker[counter] = jMarker;
            counter++;
          }
        
        /*--- To aglomerate a vertex it must have only one physical bc!!
         This can be improved. If there is only a marker, it is a good
         candidate for agglomeration ---*/
        
        if (counter == 1) agglomerate_seed = true;
        
        /*--- If there are two markers, we will aglomerate if one of the
         marker is SEND_RECEIVE ---*/
        
        if (counter == 2) {
          if ((config->GetMarker_All_Boundary(copy_marker[0]) == SEND_RECEIVE) ||
              (config->GetMarker_All_Boundary(copy_marker[1]) == SEND_RECEIVE)) agglomerate_seed = true;
          else agglomerate_seed = false;
        }
        
        /*--- If there are more than 2 markers on, the aglomeration will be discarted ---*/
        
        if (counter > 2) agglomerate_seed = false;
        
        /*--- If the seed can be agglomerated, we try to agglomerate more points ---*/
        
        if (agglomerate_seed) {
          
          /*--- Now we do a sweep over all the nodes that surround the seed point ---*/
          
          for (iNode = 0; iNode <	fine_grid->node[iPoint]->GetnPoint(); iNode ++) {
            
            agglomerate_CV = false;
            CVPoint = fine_grid->node[iPoint]->GetPoint(iNode);
            
            /*--- Determine if the CVPoint can be agglomerated ---*/
            
            agglomerate_CV = SetBoundAgglomeration(CVPoint, marker_seed, fine_grid, config);
            
            /*--- The new point can be agglomerated ---*/
            
            if (agglomerate_CV && (nChildren < max_children))  {
              
              /*--- We set the value of the parent ---*/
              
              fine_grid->node[CVPoint]->SetParent_CV(Index_CoarseCV);
              
              /*--- We set the value of the child ---*/
              
              node[Index_CoarseCV]->SetChildren_CV(nChildren,CVPoint);
              nChildren++;
            }
          }
          
          Suitable_Indirect_Neighbors.clear();
          if (fine_grid->node[iPoint]->GetAgglomerate_Indirect())
            SetSuitableNeighbors(&Suitable_Indirect_Neighbors, iPoint, Index_CoarseCV, fine_grid);
          
          /*--- Now we do a sweep over all the indirect nodes that can be added ---*/
          
          for (iNode = 0; iNode <	Suitable_Indirect_Neighbors.size(); iNode ++) {
            agglomerate_CV = false;
            CVPoint = Suitable_Indirect_Neighbors[iNode];
            
            /*--- Determine if the CVPoint can be agglomerated ---*/
            
            agglomerate_CV = SetBoundAgglomeration(CVPoint, marker_seed, fine_grid, config);
            
            /*--- The new point can be agglomerated ---*/
            
            if (agglomerate_CV && (nChildren < max_children))  {
              
              /*--- We set the value of the parent ---*/
              
              fine_grid->node[CVPoint]->SetParent_CV(Index_CoarseCV);
              
              /*--- We set the indirect agglomeration information ---*/
              
              if (fine_grid->node[CVPoint]->GetAgglomerate_Indirect())
                node[Index_CoarseCV]->SetAgglomerate_Indirect(true);
              
              /*--- We set the value of the child ---*/
              
              node[Index_CoarseCV]->SetChildren_CV(nChildren,CVPoint);
              nChildren++;
            }
          }
        }
        
        /*--- Update the number of child of the control volume ---*/
        
        node[Index_CoarseCV]->SetnChildren_CV(nChildren);
        Index_CoarseCV++;
      }
    }
  }
  
  /*--- Agglomerate all the nodes that have more than one physical boundary condition,
   Maybe here we can add the posibility of merging the vertex that have the same number,
   and kind  of markers---*/
  
  for (iMarker = 0; iMarker < fine_grid->GetnMarker(); iMarker++)
    for (iVertex = 0; iVertex < fine_grid->GetnVertex(iMarker); iVertex++) {
      iPoint = fine_grid->vertex[iMarker][iVertex]->GetNode();
      if ((fine_grid->node[iPoint]->GetAgglomerate() == false) &&
          (fine_grid->node[iPoint]->GetDomain())) {
        fine_grid->node[iPoint]->SetParent_CV(Index_CoarseCV);
        node[Index_CoarseCV]->SetChildren_CV(0,iPoint);
        node[Index_CoarseCV]->SetnChildren_CV(1);
        Index_CoarseCV++;
      }
    }
  
  /*--- Update the queue with the results from the boundary agglomeration ---*/
  
  for (iPoint = 0; iPoint < fine_grid->GetnPoint(); iPoint ++) {
    
    /*--- The CV has been agglomerated, remove form the list ---*/
    
    if (fine_grid->node[iPoint]->GetAgglomerate() == true) {
      
      MGQueue_InnerCV.RemoveCV(iPoint);
      
    }
    
    else {
      
      /*--- Count the number of agglomerated neighbors, and modify the queue ---*/
      
      priority = 0;
      for (iNode = 0; iNode <	fine_grid->node[iPoint]->GetnPoint(); iNode ++) {
        jPoint = fine_grid->node[iPoint]->GetPoint(iNode);
        if (fine_grid->node[jPoint]->GetAgglomerate() == true) priority++;
      }
      MGQueue_InnerCV.MoveCV(iPoint, priority);
    }
  }
  
  /*--- Agglomerate the domain nodes ---*/
  
  iteration = 0;
  while (!MGQueue_InnerCV.EmptyQueue() && (iteration < fine_grid->GetnPoint())) {
    
    iPoint = MGQueue_InnerCV.NextCV();
    iteration ++;
    
    /*--- If the element has not being previously agglomerated, belongs to the physical domain,
     and satisfies several geometrical criteria then the seed CV is acepted for agglomeration ---*/
    
    if ((fine_grid->node[iPoint]->GetAgglomerate() == false) &&
        (fine_grid->node[iPoint]->GetDomain()) &&
        (GeometricalCheck(iPoint, fine_grid, config))) {
      
      nChildren = 1;
      
      /*--- We set an index for the parent control volume ---*/
      
      fine_grid->node[iPoint]->SetParent_CV(Index_CoarseCV);
      
      /*--- We add the seed point (child) to the parent control volume ---*/
      
      node[Index_CoarseCV]->SetChildren_CV(0, iPoint);
      
      /*--- Update the queue with the seed point (remove the seed and
       increase the priority of the neighbors) ---*/
      
      MGQueue_InnerCV.Update(iPoint, fine_grid);
      
      /*--- Now we do a sweep over all the nodes that surround the seed point ---*/
      
      for (iNode = 0; iNode <	fine_grid->node[iPoint]->GetnPoint(); iNode ++) {
        
        CVPoint = fine_grid->node[iPoint]->GetPoint(iNode);
        
        /*--- Determine if the CVPoint can be agglomerated ---*/
        
        if ((fine_grid->node[CVPoint]->GetAgglomerate() == false) &&
            (fine_grid->node[CVPoint]->GetDomain()) &&
            (GeometricalCheck(CVPoint, fine_grid, config))) {
          
          /*--- The new point can be agglomerated, note that the applicability of
           max_children depend on the seed ---*/
          
          if (nChildren < max_children)  {
            
            /*--- We set the value of the parent ---*/
            
            fine_grid->node[CVPoint]->SetParent_CV(Index_CoarseCV);
            
            /*--- We set the value of the child ---*/
            
            node[Index_CoarseCV]->SetChildren_CV(nChildren, CVPoint);
            nChildren++;
            
            /*--- Update the queue with the new control volume (remove the CV and
             increase the priority of the neighbors) ---*/
            
            MGQueue_InnerCV.Update(CVPoint, fine_grid);
            
          }
        }
        
      }
      
      /*--- Subrotuine to identify the indirect neighbors ---*/
      
      Suitable_Indirect_Neighbors.clear();
      if (fine_grid->node[iPoint]->GetAgglomerate_Indirect())
        SetSuitableNeighbors(&Suitable_Indirect_Neighbors, iPoint, Index_CoarseCV, fine_grid);
      
      /*--- Now we do a sweep over all the indirect nodes that can be added ---*/
      
      for (iNode = 0; iNode <	Suitable_Indirect_Neighbors.size(); iNode ++) {
        
        agglomerate_CV = false;
        CVPoint = Suitable_Indirect_Neighbors[iNode];
        
        /*--- Determine if the CVPoint can be agglomerated ---*/
        
        if ((fine_grid->node[CVPoint]->GetAgglomerate() == false) && (fine_grid->node[CVPoint]->GetDomain()))
          agglomerate_CV = true;
        
        /*--- The new point can be agglomerated ---*/
        
        if ((agglomerate_CV) && (nChildren < max_children))  {
          
          /*--- We set the value of the parent ---*/
          
          fine_grid->node[CVPoint]->SetParent_CV(Index_CoarseCV);
          
          /*--- We set the indirect agglomeration information ---*/
          
          if (fine_grid->node[CVPoint]->GetAgglomerate_Indirect())
            node[Index_CoarseCV]->SetAgglomerate_Indirect(true);
          
          /*--- We set the value of the child ---*/
          
          node[Index_CoarseCV]->SetChildren_CV(nChildren, CVPoint);
          nChildren++;
          
          /*--- Update the queue with the new control volume (remove the CV and
           increase the priority of the neighbors) ---*/
          
          MGQueue_InnerCV.Update(CVPoint, fine_grid);
          
        }
      }
      
      /*--- Update the number of control of childrens ---*/
      
      node[Index_CoarseCV]->SetnChildren_CV(nChildren);
      Index_CoarseCV++;
    }
    else {
      
      /*--- The seed point can not be agglomerated because of size, domain, streching, etc.
       move the point to the lowest priority ---*/
      
      MGQueue_InnerCV.MoveCV(iPoint, -1);
    }
    
  }
  
  /*--- Add all the elements that have not being agglomerated, in the previous stage ---*/
  
  for (iPoint = 0; iPoint < fine_grid->GetnPoint(); iPoint ++) {
    if ((fine_grid->node[iPoint]->GetAgglomerate() == false) && (fine_grid->node[iPoint]->GetDomain())) {
      
      nChildren = 1;
      fine_grid->node[iPoint]->SetParent_CV(Index_CoarseCV);
      if (fine_grid->node[iPoint]->GetAgglomerate_Indirect())
        node[Index_CoarseCV]->SetAgglomerate_Indirect(true);
      node[Index_CoarseCV]->SetChildren_CV(0,iPoint);
      node[Index_CoarseCV]->SetnChildren_CV(nChildren);
      Index_CoarseCV++;
      
    }
  }
  
  nPointDomain = Index_CoarseCV;
  
  /*--- Dealing with MPI parallelization, the objective is that the received nodes must be agglomerated
   in the same way as the donor nodes. Send the node agglomeration information of the donor
   (parent and children), Sending only occurs with MPI ---*/
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    
    if ((config->GetMarker_All_Boundary(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
      
      MarkerS = iMarker;  MarkerR = iMarker+1;
      
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
      
      nVertexS = fine_grid->nVertex[MarkerS];   nVertexR = fine_grid->nVertex[MarkerR];
      nBufferS_Vector = nVertexS;               nBufferR_Vector = nVertexR;
      
      /*--- Allocate Receive and send buffers  ---*/
      
      Buffer_Receive_Children = new unsigned long [nBufferR_Vector];
      Buffer_Send_Children = new unsigned long [nBufferS_Vector];
      
      Buffer_Receive_Parent = new unsigned long [nBufferR_Vector];
      Buffer_Send_Parent = new unsigned long [nBufferS_Vector];
      
      /*--- Copy the information that should be sended ---*/
      
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = fine_grid->vertex[MarkerS][iVertex]->GetNode();
        Buffer_Send_Children[iVertex] = iPoint;
        Buffer_Send_Parent[iVertex] = fine_grid->node[iPoint]->GetParent_CV();
      }
      
      /*--- Receive information without MPI ---*/
      
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        Buffer_Receive_Children[iVertex] = Buffer_Send_Children[iVertex];
        Buffer_Receive_Parent[iVertex] = Buffer_Send_Parent[iVertex];
      }
      
      /*--- Deallocate send buffer ---*/
      
      delete [] Buffer_Send_Children;
      delete [] Buffer_Send_Parent;
      
      /*--- Create a list of the parent nodes without repeated parents ---*/
      
      Aux_Parent.clear();
      for (iVertex = 0; iVertex < nVertexR; iVertex++)
        Aux_Parent.push_back (Buffer_Receive_Parent[iVertex]);
      
      sort(Aux_Parent.begin(), Aux_Parent.end());
      it = unique(Aux_Parent.begin(), Aux_Parent.end());
      Aux_Parent.resize(it - Aux_Parent.begin());
      
      /*--- Allocate some structures ---*/
      
      Parent_Remote = new unsigned long[nVertexR];
      Children_Remote = new unsigned long[nVertexR];
      Parent_Local = new unsigned long[nVertexR];
      Children_Local = new unsigned long[nVertexR];
      
      /*--- Create the local vector and remote for the parents and the children ---*/
      
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        Parent_Remote[iVertex] = Buffer_Receive_Parent[iVertex];
        
        /*--- We use the same sorting as in the donor domain ---*/
        
        for (jVertex = 0; jVertex < Aux_Parent.size(); jVertex++) {
          if (Parent_Remote[iVertex] == Aux_Parent[jVertex]) {
            Parent_Local[iVertex] = jVertex + Index_CoarseCV;
            break;
          }
        }
        
        Children_Remote[iVertex] = Buffer_Receive_Children[iVertex];
        Children_Local[iVertex] = fine_grid->vertex[MarkerR][iVertex]->GetNode();
        
      }
      
      Index_CoarseCV += Aux_Parent.size();
      
      nChildren_MPI = new unsigned short [Index_CoarseCV];
      for (iParent = 0; iParent < Index_CoarseCV; iParent++)
        nChildren_MPI[iParent] = 0;
      
      /*--- Create the final structure ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Be careful, it is possible that a node change the agglomeration configuration, the priority
         is always, when receive the information ---*/
        
        fine_grid->node[Children_Local[iVertex]]->SetParent_CV(Parent_Local[iVertex]);
        node[Parent_Local[iVertex]]->SetChildren_CV(nChildren_MPI[Parent_Local[iVertex]], Children_Local[iVertex]);
        nChildren_MPI[Parent_Local[iVertex]]++;
        node[Parent_Local[iVertex]]->SetnChildren_CV(nChildren_MPI[Parent_Local[iVertex]]);
        node[Parent_Local[iVertex]]->SetDomain(false);
        
      }
      
      /*--- Deallocate auxiliar structures ---*/
      
      delete[] nChildren_MPI;
      delete[] Parent_Remote;
      delete[] Children_Remote;
      delete[] Parent_Local;
      delete[] Children_Local;
      
      /*--- Deallocate receive buffer ---*/
      
      delete [] Buffer_Receive_Children;
      delete [] Buffer_Receive_Parent;
      
    }
    
  }
  
  /*--- Update the number of points after the MPI agglomeration ---*/
  
  nPoint = Index_CoarseCV;
  
  /*--- Console output with the summary of the agglomeration ---*/
  
  
  
  Local_nPointCoarse = nPoint;
  Local_nPointFine = fine_grid->GetnPoint();
  
  Global_nPointCoarse = Local_nPointCoarse;
  Global_nPointFine = Local_nPointFine;
  
  cout << "CVs of the MG level: " << Global_nPointCoarse << ". Agglom. rate 1/" << double(Global_nPointFine)/double(Global_nPointCoarse) <<". MG level: "<< iMesh <<"."<< endl;
  
}


CMultiGridGeometry::~CMultiGridGeometry(void) {
  
}

bool CMultiGridGeometry::SetBoundAgglomeration(unsigned long CVPoint, short marker_seed, CGeometry *fine_grid, CConfig *config) {
  
  bool agglomerate_CV = false;
  unsigned short counter, jMarker;
  
  /*--- Basic condition, the element has not being previously agglomerated and, it belong to the domain ---*/
  if ((fine_grid->node[CVPoint]->GetAgglomerate() == false) &&
      (fine_grid->node[CVPoint]->GetDomain()) &&
      (GeometricalCheck(CVPoint, fine_grid, config))) {
    
    /*--- If the element belong to the boundary, we must be careful ---*/
    if (fine_grid->node[CVPoint]->GetBoundary()) {
      
      /*--- Identify the markers of the vertex that we whant to agglomerate ---*/
      counter = 0;
      unsigned short copy_marker[MAX_NUMBER_MARKER];
      for (jMarker = 0; jMarker < fine_grid->GetnMarker(); jMarker ++)
        if (fine_grid->node[CVPoint]->GetVertex(jMarker) != -1) {
          copy_marker[counter] = jMarker;
          counter++;
        }
      
      /*--- The basic condition is that the aglomerated vertex must have the same physical marker,
       but eventually a send-receive condition ---*/
      
      /*--- Only one marker in the vertex that is going to be aglomerated ---*/
      if (counter == 1) {
        
        /*--- We agglomerate if there is only a marker and is the same marker as the seed marker ---*/
        if (copy_marker[0] == marker_seed)
          agglomerate_CV = true;
        
        /*--- If there is only a marker, but the marker is the SEND_RECEIVE ---*/
        if (config->GetMarker_All_Boundary(copy_marker[0]) == SEND_RECEIVE)
          agglomerate_CV = true;
      }
      
      /*--- If there are two markers in the vertex that is going to be aglomerated ---*/
      if (counter == 2) {
        
        /*--- First we verify that the seed is a physical boundary ---*/
        if (config->GetMarker_All_Boundary(marker_seed) != SEND_RECEIVE) {
          
          /*--- Then we check that one of the marker is equal to the seed marker, and the other is send/receive ---*/
          if (((copy_marker[0] == marker_seed) && (config->GetMarker_All_Boundary(copy_marker[1]) == SEND_RECEIVE)) ||
              ((config->GetMarker_All_Boundary(copy_marker[0]) == SEND_RECEIVE) && (copy_marker[1] == marker_seed)))
            agglomerate_CV = true;
        }
      }
      
    }
    
    /*--- If the element belong to the domain, it is allways aglomerated ---*/
    else {
      agglomerate_CV = true;
    }
  }
  
  return agglomerate_CV;
}


bool CMultiGridGeometry::GeometricalCheck(unsigned long iPoint, CGeometry *fine_grid, CConfig *config) {
  
  double max_dimension = 1.0/config->GetMaxDimension();
  
  /*--- Evaluate the total size of the element ---*/
  bool Volume = true;
  double ratio = pow(fine_grid->node[iPoint]->GetVolume(), 1.0/double(nDim))*max_dimension;
  double limit = pow(config->GetDomainVolume(), 1.0/double(nDim));
  if ( ratio > limit ) Volume = false;
  
  /*--- Evaluate the stretching of the element ---*/
  bool Stretching = true;
  /*	unsigned short iNode, iDim;
   unsigned long jPoint;
   double *Coord_i = fine_grid->node[iPoint]->GetCoord();
   double max_dist = 0.0 ; double min_dist = 1E20;
   for (iNode = 0; iNode <	fine_grid->node[iPoint]->GetnPoint(); iNode ++) {
   jPoint = fine_grid->node[iPoint]->GetPoint(iNode);
   double *Coord_j = fine_grid->node[jPoint]->GetCoord();
   double distance = 0.0;
   for (iDim = 0; iDim < nDim; iDim++)
   distance += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
   distance = sqrt(distance);
   max_dist = max(distance, max_dist);
   min_dist = min(distance, min_dist);
   }
   if ( max_dist/min_dist > 100.0 ) Stretching = false;*/
  
  return (Stretching && Volume);
  
}

void CMultiGridGeometry::SetSuitableNeighbors(vector<unsigned long> *Suitable_Indirect_Neighbors, unsigned long iPoint,
                                              unsigned long Index_CoarseCV, CGeometry *fine_grid) {
  
  unsigned long jPoint, kPoint,  iOriginNeighbor, lPoint;
  unsigned short iNode, jNode, iNeighbor, jNeighbor, kNode;
  bool SecondNeighborSeed, check_1, ThirdNeighborSeed;
  
  /*--- Create a list with the first neighbors, including the seed ---*/
  vector<unsigned long> First_Neighbor_Points;
  First_Neighbor_Points.push_back(iPoint);
  for (iNode = 0; iNode <	fine_grid->node[iPoint]->GetnPoint(); iNode ++) {
    jPoint = fine_grid->node[iPoint]->GetPoint(iNode);
    First_Neighbor_Points.push_back(jPoint);
  }
  
  /*--- Create a list with the second neighbors, without first, and seed neighbors ---*/
  vector<unsigned long> Second_Neighbor_Points, Second_Origin_Points, Suitable_Second_Neighbors;
  
  for (iNode = 0; iNode <	fine_grid->node[iPoint]->GetnPoint(); iNode ++) {
    jPoint = fine_grid->node[iPoint]->GetPoint(iNode);
    for (jNode = 0; jNode <	fine_grid->node[jPoint]->GetnPoint(); jNode ++) {
      kPoint = fine_grid->node[jPoint]->GetPoint(jNode);
      
      /*--- Check that the second neighbor do not belong to the first neighbor or the seed ---*/
      SecondNeighborSeed = true;
      for (iNeighbor = 0; iNeighbor <	First_Neighbor_Points.size(); iNeighbor ++)
        if (kPoint == First_Neighbor_Points[iNeighbor]) {
          SecondNeighborSeed = false;
          break;
        }
      
      if (SecondNeighborSeed) {
        Second_Neighbor_Points.push_back(kPoint);
        Second_Origin_Points.push_back(jPoint);
      }
      
    }
  }
  
  /*---  Identify those second neighbors that are repeated (candidate to be added) ---*/
  for (iNeighbor = 0; iNeighbor <	Second_Neighbor_Points.size(); iNeighbor ++)
    for (jNeighbor = 0; jNeighbor <	Second_Neighbor_Points.size(); jNeighbor ++)
      
    /*--- Repeated second neighbor with different origin ---*/
      if ((Second_Neighbor_Points[iNeighbor] == Second_Neighbor_Points[jNeighbor]) &&
          (Second_Origin_Points[iNeighbor] != Second_Origin_Points[jNeighbor]) &&
          (iNeighbor < jNeighbor)) {
        
        /*--- Check that the origin nodes are not neighbor ---*/
        check_1 = true;
        for (iNode = 0; iNode <	fine_grid->node[Second_Origin_Points[iNeighbor]]->GetnPoint(); iNode ++) {
          iOriginNeighbor = fine_grid->node[Second_Origin_Points[iNeighbor]]->GetPoint(iNode);
          if (iOriginNeighbor == Second_Origin_Points[jNeighbor]) {
            check_1 = false;
            break;
          }
        }
        
        if (check_1) {
          Suitable_Indirect_Neighbors->push_back(Second_Neighbor_Points[iNeighbor]);
          
          /*--- Create alist with the suitable second neighbor, that we will use
           to compute the third neighbors --*/
          Suitable_Second_Neighbors.push_back(Second_Neighbor_Points[iNeighbor]);
        }
      }
  
  vector<unsigned long>::iterator it;
  
  /*--- Remove repeated from the suitable second neighbors ---*/
  sort(Suitable_Second_Neighbors.begin(), Suitable_Second_Neighbors.end());
  it = unique(Suitable_Second_Neighbors.begin(), Suitable_Second_Neighbors.end());
  Suitable_Second_Neighbors.resize(it - Suitable_Second_Neighbors.begin());
  
  /*--- Remove repeated from first neighbors ---*/
  sort(First_Neighbor_Points.begin(), First_Neighbor_Points.end());
  it = unique(First_Neighbor_Points.begin(), First_Neighbor_Points.end());
  First_Neighbor_Points.resize(it - First_Neighbor_Points.begin());
  
  /*--- Create a list with the third neighbors, without first, second, and seed neighbors ---*/
  vector<unsigned long> Third_Neighbor_Points, Third_Origin_Points;
  
  for (jNode = 0; jNode <	Suitable_Second_Neighbors.size(); jNode ++) {
    kPoint = Suitable_Second_Neighbors[jNode];
    
    for (kNode = 0; kNode <	fine_grid->node[kPoint]->GetnPoint(); kNode ++) {
      lPoint = fine_grid->node[kPoint]->GetPoint(kNode);
      
      /*--- Check that the third neighbor do not belong to the first neighbors or the seed ---*/
      ThirdNeighborSeed = true;
      
      for (iNeighbor = 0; iNeighbor <	First_Neighbor_Points.size(); iNeighbor ++)
        if (lPoint == First_Neighbor_Points[iNeighbor]) {
          ThirdNeighborSeed = false;
          break;
        }
      
      /*--- Check that the third neighbor do not belong to the second neighbors ---*/
      for (iNeighbor = 0; iNeighbor <	Suitable_Second_Neighbors.size(); iNeighbor ++)
        if (lPoint == Suitable_Second_Neighbors[iNeighbor]) {
          ThirdNeighborSeed = false;
          break;
        }
      
      if (ThirdNeighborSeed) {
        Third_Neighbor_Points.push_back(lPoint);
        Third_Origin_Points.push_back(kPoint);
      }
      
    }
  }
  
  /*---  Identify those third neighbors that are repeated (candidate to be added) ---*/
  for (iNeighbor = 0; iNeighbor <	Third_Neighbor_Points.size(); iNeighbor ++)
    for (jNeighbor = 0; jNeighbor <	Third_Neighbor_Points.size(); jNeighbor ++)
      
    /*--- Repeated second neighbor with different origin ---*/
      if ((Third_Neighbor_Points[iNeighbor] == Third_Neighbor_Points[jNeighbor]) &&
          (Third_Origin_Points[iNeighbor] != Third_Origin_Points[jNeighbor]) &&
          (iNeighbor < jNeighbor)) {
        
        /*--- Check that the origin nodes are not neighbor ---*/
        check_1 = true;
        for (iNode = 0; iNode <	fine_grid->node[Third_Origin_Points[iNeighbor]]->GetnPoint(); iNode ++) {
          iOriginNeighbor = fine_grid->node[Third_Origin_Points[iNeighbor]]->GetPoint(iNode);
          if (iOriginNeighbor == Third_Origin_Points[jNeighbor]) {
            check_1 = false;
            break;
          }
        }
        
        if (check_1)
          Suitable_Indirect_Neighbors->push_back(Third_Neighbor_Points[iNeighbor]);
        
      }
  
  /*--- Remove repeated from Suitable Indirect Neighbors List ---*/
  sort(Suitable_Indirect_Neighbors->begin(), Suitable_Indirect_Neighbors->end());
  it = unique(Suitable_Indirect_Neighbors->begin(), Suitable_Indirect_Neighbors->end());
  Suitable_Indirect_Neighbors->resize(it - Suitable_Indirect_Neighbors->begin());
  
}



void CMultiGridGeometry::SetPsuP(CGeometry *fine_grid) {
  unsigned long iFinePoint, iFinePoint_Neighbor, iParent, iCoarsePoint;
  unsigned short iChildren, iNode;
  
  /*--- Set the point suronfding a point ---*/
  for (iCoarsePoint = 0; iCoarsePoint < nPoint; iCoarsePoint ++)
    for (iChildren = 0; iChildren <  node[iCoarsePoint]->GetnChildren_CV(); iChildren ++) {
      iFinePoint = node[iCoarsePoint]->GetChildren_CV(iChildren);
      for (iNode = 0; iNode < fine_grid->node[iFinePoint]->GetnPoint(); iNode ++) {
        iFinePoint_Neighbor = fine_grid->node[iFinePoint]->GetPoint(iNode);
        iParent = fine_grid->node[iFinePoint_Neighbor]->GetParent_CV();
        if (iParent != iCoarsePoint) node[iCoarsePoint]->SetPoint(iParent);
      }
    }
  
  /*--- Set the number of neighbors variable, this is
   important for JST and multigrid in parallel ---*/
  for (iCoarsePoint = 0; iCoarsePoint < nPoint; iCoarsePoint ++)
    node[iCoarsePoint]->SetnNeighbor(node[iCoarsePoint]->GetnPoint());
  
}

void CMultiGridGeometry::SetVertex(CGeometry *fine_grid, CConfig *config) {
  unsigned long  iVertex, iFinePoint, iCoarsePoint;
  unsigned short iMarker, iMarker_Tag, iChildren;
  
  nMarker = fine_grid->GetnMarker();
  
  /*--- If any children node belong to the boundary then the entire control
   volume will belong to the boundary ---*/
  for (iCoarsePoint = 0; iCoarsePoint < nPoint; iCoarsePoint ++)
    for (iChildren = 0; iChildren <	node[iCoarsePoint]->GetnChildren_CV(); iChildren ++) {
      iFinePoint = node[iCoarsePoint]->GetChildren_CV(iChildren);
      if (fine_grid->node[iFinePoint]->GetBoundary()) {
        node[iCoarsePoint]->SetBoundary(nMarker);
        break;
      }
    }
  
  vertex = new CVertex**[nMarker];
  nVertex = new unsigned long [nMarker];
  
  Tag_to_Marker = new string [MAX_INDEX_VALUE];
  for (iMarker_Tag = 0; iMarker_Tag < MAX_INDEX_VALUE; iMarker_Tag++)
    Tag_to_Marker[iMarker_Tag] = fine_grid->GetMarker_Tag(iMarker_Tag);
  
  /*--- Compute the number of vertices to do the dimensionalization ---*/
  for (iMarker = 0; iMarker < nMarker; iMarker++) nVertex[iMarker] = 0;
  
  
  for (iCoarsePoint = 0; iCoarsePoint < nPoint; iCoarsePoint ++) {
    if (node[iCoarsePoint]->GetBoundary()) {
      for (iChildren = 0; iChildren <	node[iCoarsePoint]->GetnChildren_CV(); iChildren ++) {
        iFinePoint = node[iCoarsePoint]->GetChildren_CV(iChildren);
        for (iMarker = 0; iMarker < nMarker; iMarker ++) {
          if ((fine_grid->node[iFinePoint]->GetVertex(iMarker) != -1) && (node[iCoarsePoint]->GetVertex(iMarker) == -1)) {
            iVertex = nVertex[iMarker];
            node[iCoarsePoint]->SetVertex(iVertex,iMarker);
            nVertex[iMarker]++;
          }
        }
      }
    }
  }
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    vertex[iMarker] = new CVertex* [fine_grid->GetnVertex(iMarker)+1];
    nVertex[iMarker] = 0;
  }
  
  for (iCoarsePoint = 0; iCoarsePoint < nPoint; iCoarsePoint ++)
    if (node[iCoarsePoint]->GetBoundary())
      for (iMarker = 0; iMarker < nMarker; iMarker ++)
        node[iCoarsePoint]->SetVertex(-1,iMarker);
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) nVertex[iMarker] = 0;
  
  for (iCoarsePoint = 0; iCoarsePoint < nPoint; iCoarsePoint ++)
    if (node[iCoarsePoint]->GetBoundary())
      for (iChildren = 0; iChildren <	node[iCoarsePoint]->GetnChildren_CV(); iChildren ++) {
        iFinePoint = node[iCoarsePoint]->GetChildren_CV(iChildren);
        for (iMarker = 0; iMarker < fine_grid->GetnMarker(); iMarker ++) {
          if ((fine_grid->node[iFinePoint]->GetVertex(iMarker) != -1) && (node[iCoarsePoint]->GetVertex(iMarker) == -1)) {
            iVertex = nVertex[iMarker];
            vertex[iMarker][iVertex] = new CVertex(iCoarsePoint, nDim);
            node[iCoarsePoint]->SetVertex(iVertex,iMarker);
            
            /*--- Set the transformation to apply ---*/
            unsigned long ChildVertex = fine_grid->node[iFinePoint]->GetVertex(iMarker);
            unsigned short RotationKind = fine_grid->vertex[iMarker][ChildVertex]->GetRotation_Type();
            vertex[iMarker][iVertex]->SetRotation_Type(RotationKind);
            nVertex[iMarker]++;
          }
        }
      }
}

void CMultiGridGeometry::MatchNearField(CConfig *config) {
  
  unsigned short iMarker;
  unsigned long iVertex, iPoint;
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if (config->GetMarker_All_Boundary(iMarker) == NEARFIELD_BOUNDARY)
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
        iPoint = vertex[iMarker][iVertex]->GetNode();
        vertex[iMarker][iVertex]->SetDonorPoint(iPoint);
      }
  
}

void CMultiGridGeometry::MatchInterface(CConfig *config) {
  
  unsigned short iMarker;
  unsigned long iVertex, iPoint;
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if (config->GetMarker_All_Boundary(iMarker) == INTERFACE_BOUNDARY)
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
        iPoint = vertex[iMarker][iVertex]->GetNode();
        vertex[iMarker][iVertex]->SetDonorPoint(iPoint);
      }
  
}


void CMultiGridGeometry::SetControlVolume(CConfig *config, CGeometry *fine_grid, unsigned short action) {
  
  unsigned long iFinePoint,iFinePoint_Neighbor, iCoarsePoint, iEdge, iParent, iPoint, jPoint;
  long FineEdge, CoarseEdge;
  unsigned short iChildren, iNode, iDim;
  bool change_face_orientation;
  double *Normal, Coarse_Volume, Area, *NormalFace = NULL;
  Normal = new double [nDim];
  
  /*--- Compute the area of the coarse volume ---*/
  for (iCoarsePoint = 0; iCoarsePoint < nPoint; iCoarsePoint ++) {
    node[iCoarsePoint]->SetVolume(0.0);
    Coarse_Volume = 0.0;
    for (iChildren = 0; iChildren < node[iCoarsePoint]->GetnChildren_CV(); iChildren ++) {
      iFinePoint = node[iCoarsePoint]->GetChildren_CV(iChildren);
      Coarse_Volume += fine_grid->node[iFinePoint]->GetVolume();
    }
    node[iCoarsePoint]->SetVolume(Coarse_Volume);
  }
  
  /*--- Update or not the values of faces at the edge ---*/
  if (action != ALLOCATE) {
    for(iEdge=0; iEdge < nEdge; iEdge++)
      edge[iEdge]->SetZeroValues();
  }
  
  for (iCoarsePoint = 0; iCoarsePoint < nPoint; iCoarsePoint ++)
    for (iChildren = 0; iChildren < node[iCoarsePoint]->GetnChildren_CV(); iChildren ++) {
      iFinePoint = node[iCoarsePoint]->GetChildren_CV(iChildren);
      
      for (iNode = 0; iNode < fine_grid->node[iFinePoint]->GetnPoint(); iNode ++) {
        iFinePoint_Neighbor = fine_grid->node[iFinePoint]->GetPoint(iNode);
        iParent = fine_grid->node[iFinePoint_Neighbor]->GetParent_CV();
        if ((iParent != iCoarsePoint) && (iParent < iCoarsePoint)) {
          
          FineEdge = fine_grid->FindEdge(iFinePoint, iFinePoint_Neighbor);
          
          change_face_orientation = false;
          if (iFinePoint < iFinePoint_Neighbor) change_face_orientation = true;
          
          CoarseEdge = FindEdge(iParent, iCoarsePoint);
          
          fine_grid->edge[FineEdge]->GetNormal(Normal);
          
          if (change_face_orientation) {
            for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
            edge[CoarseEdge]->AddNormal(Normal);
          }
          else {
            edge[CoarseEdge]->AddNormal(Normal);
          }
        }
      }
    }
  delete[] Normal;
  
  /*--- Check if there isn't any element with only one neighbor...
   a CV that is inside another CV ---*/
  for (iPoint = 0; iPoint < nPoint; iPoint ++) {
    if (node[iPoint]->GetnPoint() == 1) {
      jPoint = node[iPoint]->GetPoint(0);
      node[jPoint]->AddVolume(node[iPoint]->GetVolume());
    }
  }
  
  /*--- Check if there is a normal with null area ---*/
  for (iEdge = 0; iEdge < nEdge; iEdge++) {
    NormalFace = edge[iEdge]->GetNormal();
    Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += NormalFace[iDim]*NormalFace[iDim];
    Area = sqrt(Area);
    if (Area == 0.0) for (iDim = 0; iDim < nDim; iDim++) NormalFace[iDim] = EPS*EPS;
  }
  
}

void CMultiGridGeometry::SetBoundControlVolume(CConfig *config, CGeometry *fine_grid, unsigned short action) {
  unsigned long iCoarsePoint, iFinePoint, FineVertex, iVertex;
  unsigned short iMarker, iChildren, iDim;
  double *Normal, Area, *NormalFace = NULL;
  
  Normal = new double [nDim];
  
  if (action != ALLOCATE) {
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++)
        vertex[iMarker][iVertex]->SetZeroValues();
  }
  
  for (iMarker = 0; iMarker < nMarker; iMarker ++)
    for(iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
      iCoarsePoint = vertex[iMarker][iVertex]->GetNode();
      for (iChildren = 0; iChildren < node[iCoarsePoint]->GetnChildren_CV(); iChildren ++) {
        iFinePoint = node[iCoarsePoint]->GetChildren_CV(iChildren);
        if (fine_grid->node[iFinePoint]->GetVertex(iMarker)!=-1) {
          FineVertex = fine_grid->node[iFinePoint]->GetVertex(iMarker);
          fine_grid->vertex[iMarker][FineVertex]->GetNormal(Normal);
          vertex[iMarker][iVertex]->AddNormal(Normal);
        }
      }
    }
  
  delete[] Normal;
  
  /*--- Check if there is a normal with null area ---*/
  for (iMarker = 0; iMarker < nMarker; iMarker ++)
    for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
      NormalFace = vertex[iMarker][iVertex]->GetNormal();
      Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += NormalFace[iDim]*NormalFace[iDim];
      Area = sqrt(Area);
      if (Area == 0.0) for (iDim = 0; iDim < nDim; iDim++) NormalFace[iDim] = EPS*EPS;
    }
  
}

void CMultiGridGeometry::SetCoord(CGeometry *geometry) {
  unsigned long Point_Fine, Point_Coarse;
  unsigned short iChildren, iDim;
  double Area_Parent, Area_Children;
  double *Coordinates_Fine, *Coordinates;
  Coordinates = new double[nDim];
  
  for (Point_Coarse = 0; Point_Coarse < GetnPoint(); Point_Coarse++) {
    Area_Parent = node[Point_Coarse]->GetVolume();
    for (iDim = 0; iDim < nDim; iDim++) Coordinates[iDim] = 0.0;
    for (iChildren = 0; iChildren < node[Point_Coarse]->GetnChildren_CV(); iChildren++) {
      Point_Fine = node[Point_Coarse]->GetChildren_CV(iChildren);
      Area_Children = geometry->node[Point_Fine]->GetVolume();
      Coordinates_Fine = geometry->node[Point_Fine]->GetCoord();
      for (iDim = 0; iDim < nDim; iDim++)
        Coordinates[iDim] += Coordinates_Fine[iDim]*Area_Children/Area_Parent;
    }
    for (iDim = 0; iDim < nDim; iDim++)
      node[Point_Coarse]->SetCoord(iDim,Coordinates[iDim]);
  }
  delete[] Coordinates;
}


void CMultiGridGeometry::FindNormal_Neighbor(CConfig *config) {
  
  unsigned short iMarker, iDim;
  unsigned long iPoint, iVertex;
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    
    if (config->GetMarker_All_Boundary(iMarker) != SEND_RECEIVE &&
        config->GetMarker_All_Boundary(iMarker) != INTERFACE_BOUNDARY &&
        config->GetMarker_All_Boundary(iMarker) != NEARFIELD_BOUNDARY ) {
      
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
        
        iPoint = vertex[iMarker][iVertex]->GetNode();
        
        /*--- If the node belong to the domain ---*/
        if (node[iPoint]->GetDomain()) {
          
          /*--- Compute closest normal neighbor ---*/
          double cos_max, scalar_prod, norm_vect, norm_Normal, cos_alpha, diff_coord;
          unsigned long Point_Normal = 0, jPoint;
          unsigned short iNeigh;
          double *Normal = vertex[iMarker][iVertex]->GetNormal();
          cos_max = -1.0;
          for (iNeigh = 0; iNeigh < node[iPoint]->GetnPoint(); iNeigh++) {
            jPoint = node[iPoint]->GetPoint(iNeigh);
            scalar_prod = 0.0; norm_vect = 0.0; norm_Normal = 0.0;
            for(iDim = 0; iDim < nDim; iDim++) {
              diff_coord = node[jPoint]->GetCoord(iDim)-node[iPoint]->GetCoord(iDim);
              scalar_prod += diff_coord*Normal[iDim];
              norm_vect += diff_coord*diff_coord;
              norm_Normal += Normal[iDim]*Normal[iDim];
            }
            norm_vect = sqrt(norm_vect);
            norm_Normal = sqrt(norm_Normal);
            cos_alpha = scalar_prod/(norm_vect*norm_Normal);
            
            /*--- Get maximum cosine (not minimum because normals are oriented inwards) ---*/
            if (cos_alpha >= cos_max) {
              Point_Normal = jPoint;
              cos_max = cos_alpha;
            }
          }
          vertex[iMarker][iVertex]->SetNormal_Neighbor(Point_Normal);
        }
      }
    }
  }
}


void CMultiGridGeometry::SetGeometryPlanes(CConfig *config) {
  bool loop_on;
  unsigned short iMarker = 0;
  double auxXCoord, auxYCoord, auxZCoord,	*Face_Normal = NULL, auxArea, *Xcoord = NULL, *Ycoord = NULL, *Zcoord = NULL, *FaceArea = NULL;
  unsigned long jVertex, iVertex,ixCoord, iPoint, iVertex_Wall, nVertex_Wall = 0;
  
  /*--- Compute the total number of points on the near-field ---*/
  nVertex_Wall = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if ((config->GetMarker_All_Boundary(iMarker) == HEAT_FLUX) ||
        (config->GetMarker_All_Boundary(iMarker) == ISOTHERMAL) ||
        (config->GetMarker_All_Boundary(iMarker) == EULER_WALL))
      nVertex_Wall += nVertex[iMarker];
  
  
  /*--- Create an array with all the coordinates, points, pressures, face area,
   equivalent area, and nearfield weight ---*/
  Xcoord = new double[nVertex_Wall];
  Ycoord = new double[nVertex_Wall];
  if (nDim == 3)	Zcoord = new double[nVertex_Wall];
  FaceArea = new double[nVertex_Wall];
  
  /*--- Copy the boundary information to an array ---*/
  iVertex_Wall = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if ((config->GetMarker_All_Boundary(iMarker) == HEAT_FLUX) ||
        (config->GetMarker_All_Boundary(iMarker) == ISOTHERMAL) ||
        (config->GetMarker_All_Boundary(iMarker) == EULER_WALL))
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
        iPoint = vertex[iMarker][iVertex]->GetNode();
        Xcoord[iVertex_Wall] = node[iPoint]->GetCoord(0);
        Ycoord[iVertex_Wall] = node[iPoint]->GetCoord(1);
        if (nDim==3) Zcoord[iVertex_Wall] = node[iPoint]->GetCoord(2);
        Face_Normal = vertex[iMarker][iVertex]->GetNormal();
        FaceArea[iVertex_Wall] = fabs(Face_Normal[nDim-1]);
        iVertex_Wall ++;
      }
  
  
  //vector<double> XCoordList;
  vector<double>::iterator IterXCoordList;
  
  for (iVertex = 0; iVertex < nVertex_Wall; iVertex++)
    XCoordList.push_back(Xcoord[iVertex]);
  
  sort( XCoordList.begin(), XCoordList.end());
  IterXCoordList = unique( XCoordList.begin(), XCoordList.end());
  XCoordList.resize( IterXCoordList - XCoordList.begin() );
  
  /*--- Create vectors and distribute the values among the different PhiAngle queues ---*/
  Xcoord_plane.resize(XCoordList.size());
  Ycoord_plane.resize(XCoordList.size());
  if (nDim==3) Zcoord_plane.resize(XCoordList.size());
  FaceArea_plane.resize(XCoordList.size());
  Plane_points.resize(XCoordList.size());
  
  
  double dist_ratio;
  unsigned long iCoord;
  
  /*--- Distribute the values among the different PhiAngles ---*/
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    if (node[iPoint]->GetDomain()) {
      loop_on = true;
      for (ixCoord = 0; ixCoord < XCoordList.size()-1 && loop_on; ixCoord++) {
        dist_ratio = (node[iPoint]->GetCoord(0) - XCoordList[ixCoord])/(XCoordList[ixCoord+1]- XCoordList[ixCoord]);
        if (dist_ratio >= 0 && dist_ratio <= 1.0) {
          if (dist_ratio <= 0.5) iCoord = ixCoord;
          else iCoord = ixCoord+1;
          Xcoord_plane[iCoord].push_back(node[iPoint]->GetCoord(0) );
          Ycoord_plane[iCoord].push_back(node[iPoint]->GetCoord(1) );
          if (nDim==3) Zcoord_plane[iCoord].push_back(node[iPoint]->GetCoord(2) );
          FaceArea_plane[iCoord].push_back(node[iPoint]->GetVolume());   ///// CHECK AREA CALCULATION
          Plane_points[iCoord].push_back(iPoint );
          loop_on = false;
        }
      }
    }
  }
  
  unsigned long auxPoint;
  /*--- Order the arrays in ascending values of y ---*/
  for (ixCoord = 0; ixCoord < XCoordList.size(); ixCoord++)
    for (iVertex = 0; iVertex < Xcoord_plane[ixCoord].size(); iVertex++)
      for (jVertex = 0; jVertex < Xcoord_plane[ixCoord].size() - 1 - iVertex; jVertex++)
        if (Ycoord_plane[ixCoord][jVertex] > Ycoord_plane[ixCoord][jVertex+1]) {
          auxXCoord = Xcoord_plane[ixCoord][jVertex]; Xcoord_plane[ixCoord][jVertex] = Xcoord_plane[ixCoord][jVertex+1]; Xcoord_plane[ixCoord][jVertex+1] = auxXCoord;
          auxYCoord = Ycoord_plane[ixCoord][jVertex]; Ycoord_plane[ixCoord][jVertex] = Ycoord_plane[ixCoord][jVertex+1]; Ycoord_plane[ixCoord][jVertex+1] = auxYCoord;
          auxPoint = Plane_points[ixCoord][jVertex]; Plane_points[ixCoord][jVertex] = Plane_points[ixCoord][jVertex+1]; Plane_points[ixCoord][jVertex+1] = auxPoint;
          if (nDim==3) {
            auxZCoord = Zcoord_plane[ixCoord][jVertex]; Zcoord_plane[ixCoord][jVertex] = Zcoord_plane[ixCoord][jVertex+1]; Zcoord_plane[ixCoord][jVertex+1] = auxZCoord;
          }
          auxArea = FaceArea_plane[ixCoord][jVertex]; FaceArea_plane[ixCoord][jVertex] = FaceArea_plane[ixCoord][jVertex+1]; FaceArea_plane[ixCoord][jVertex+1] = auxArea;
        }
  
  /*--- Delete structures ---*/
  delete[] Xcoord; delete[] Ycoord;
  if (nDim==3) delete[] Zcoord;
  delete[] FaceArea;
}

CMultiGridQueue::CMultiGridQueue(unsigned long val_npoint) {
  unsigned long iPoint;
  
  nPoint = val_npoint;
  Priority = new short[nPoint];
  RightCV = new bool[nPoint];
  
  QueueCV.resize(1);
  
  /*--- Queue initialization with all the points in the finer grid ---*/
  for (iPoint = 0; iPoint < nPoint; iPoint ++) {
    QueueCV[0].push_back(iPoint);
    Priority[iPoint] = 0;
    RightCV[iPoint] = true;
  }
  
}

CMultiGridQueue::~CMultiGridQueue(void) {
  
  delete[] Priority;
  delete[] RightCV;
  
}

void CMultiGridQueue::AddCV(unsigned long val_new_point, unsigned short val_number_neighbors) {
  
  unsigned short Max_Neighbors = QueueCV.size()-1;
  
  /*--- Basic check ---*/
  if (val_new_point > nPoint) {
    cout << "The index of the CV is greater than the size of the priority list." << endl;
    exit(0);
  }
  
  /*--- Resize the list ---*/
  if (val_number_neighbors > Max_Neighbors)
    QueueCV.resize(val_number_neighbors+1);
  
  /*--- Find the point in the queue ---*/
  bool InQueue = false;
  if (Priority[val_new_point] == val_number_neighbors) InQueue = true;
  
  if (!InQueue) {
    /*--- Add the control volume, and update the priority list ---*/
    QueueCV[val_number_neighbors].push_back(val_new_point);
    Priority[val_new_point] = val_number_neighbors;
  }
  
}

void CMultiGridQueue::RemoveCV(unsigned long val_remove_point) {
  unsigned short iPoint;
  bool check;
  
  /*--- Basic check ---*/
  if (val_remove_point > nPoint) {
    cout << "The index of the CV is greater than the size of the priority list." << endl;
    exit(0);
  }
  
  /*--- Find priority of the Control Volume ---*/
  short Number_Neighbors = Priority[val_remove_point];
  if (Number_Neighbors == -1) {
    cout << "The CV "<< val_remove_point <<" is not in the priority list. (RemoveCV)" << endl;
    exit(0);
  }
  
  /*--- Find the point in the queue ---*/
  vector<unsigned long>::iterator ItQueue = find(QueueCV[Number_Neighbors].begin(),
                                                 QueueCV[Number_Neighbors].end(),
                                                 val_remove_point);
  if( ItQueue != QueueCV[Number_Neighbors].end() ) QueueCV[Number_Neighbors].erase(ItQueue);
  
  Priority[val_remove_point] = -1;
  
  /*--- Check that the size of the queue is the right one ---*/
  unsigned short Size_QueueCV = 0;
  check = false;
  for (iPoint = 0; iPoint < QueueCV.size(); iPoint ++)
    if (QueueCV[iPoint].size() != 0) { Size_QueueCV = iPoint; check = true;}
  
  /*--- Resize the queue, if check = false, the queue is empty, at least
   we need one element in the queue ---*/
  if (check) QueueCV.resize(Size_QueueCV+1);
  else QueueCV.resize(1);
  
}

void CMultiGridQueue::MoveCV(unsigned long val_move_point, short val_number_neighbors) {
  unsigned short Priority;
  
  if (val_number_neighbors < 0) {
    val_number_neighbors = 0;
    RightCV[val_move_point] = false;
  }
  else {
    Priority = val_number_neighbors;
    RightCV[val_move_point] = true;
  }
  
  /*--- Remove the control volume ---*/
  RemoveCV(val_move_point);
  
  /*--- Add a new control volume ---*/
  AddCV(val_move_point, val_number_neighbors);
  
}

void CMultiGridQueue::IncrPriorityCV(unsigned long val_incr_point) {
  
  /*--- Find the priority list ---*/
  short Number_Neighbors = Priority[val_incr_point];
  if (Number_Neighbors == -1) {
    cout << "The CV "<< val_incr_point <<" is not in the priority list. (IncrPriorityCV)" << endl;
    exit(0);
  }
  
  /*--- Remove the control volume ---*/
  RemoveCV(val_incr_point);
  
  /*--- Increase the priority ---*/
  AddCV(val_incr_point, Number_Neighbors+1);
  
}

void CMultiGridQueue::RedPriorityCV(unsigned long val_red_point) {
  
  /*--- Find the priority list ---*/
  short Number_Neighbors = Priority[val_red_point];
  if (Number_Neighbors == -1) {
    cout << "The CV "<< val_red_point <<" is not in the priority list. (RedPriorityCV)" << endl;
    exit(0);
  }
  
  if (Number_Neighbors != 0) {
    
    /*--- Remove the control volume ---*/
    RemoveCV(val_red_point);
    
    /*--- Increase the priority ---*/
    AddCV(val_red_point, Number_Neighbors-1);
    
  }
  
}

void CMultiGridQueue::VisualizeQueue(void) {
  unsigned short iPoint;
  unsigned long jPoint;
  
  cout << endl;
  for (iPoint = 0; iPoint < QueueCV.size(); iPoint ++) {
    cout << "Number of neighbors " << iPoint <<": ";
    for (jPoint = 0; jPoint < QueueCV[iPoint].size(); jPoint ++) {
      cout << QueueCV[iPoint][jPoint] << " ";
    }
    cout << endl;
  }
  
}

void CMultiGridQueue::VisualizePriority(void) {
  unsigned long iPoint;
  
  for (iPoint = 0; iPoint < nPoint; iPoint ++)
    cout << "Control Volume: " << iPoint <<" Priority: " << Priority[iPoint] << endl;
  
}

long CMultiGridQueue::NextCV(void) {
  if (QueueCV.size() != 0) return QueueCV[QueueCV.size()-1][0];
  else return -1;
}

bool CMultiGridQueue::EmptyQueue(void) {
  unsigned short iPoint;
  
  /*--- In case there is only the no agglomerated elements,
   check if they can be agglomerated or we have already finished ---*/
  bool check = true;
  
  if ( QueueCV.size() == 1 ) {
    for (iPoint = 0; iPoint < QueueCV[0].size(); iPoint ++) {
      if (RightCV[QueueCV[0][iPoint]]) { check = false; break; }
    }
  }
  else {
    for (iPoint = 1; iPoint < QueueCV.size(); iPoint ++)
      if (QueueCV[iPoint].size() != 0) { check = false; break;}
  }
  
  return check;
}

unsigned long CMultiGridQueue::TotalCV(void) {
  unsigned short iPoint;
  unsigned long TotalCV;
  
  TotalCV = 0;
  for (iPoint = 0; iPoint < QueueCV.size(); iPoint ++)
    if (QueueCV[iPoint].size() != 0) { TotalCV += QueueCV[iPoint].size(); }
  
  return TotalCV;
}

void CMultiGridQueue::Update(unsigned long iPoint, CGeometry *fine_grid) {
  unsigned short iNode;
  unsigned long jPoint;
  
  RemoveCV(iPoint);
  for (iNode = 0; iNode <	fine_grid->node[iPoint]->GetnPoint(); iNode ++) {
    jPoint = fine_grid->node[iPoint]->GetPoint(iNode);
    if (fine_grid->node[jPoint]->GetAgglomerate() == false)
      IncrPriorityCV(jPoint);
  }
  
}
