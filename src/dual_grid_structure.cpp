/*!
 * \file dual_grid_structure.cpp
 * \brief Main classes for defining the dual grid (points, vertex, and edges).
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

#include "../include/dual_grid_structure.hpp"

unsigned short CDualGrid::nDim = 0;

CDualGrid::CDualGrid(unsigned short val_nDim) { nDim = val_nDim; Color = 0;}

CDualGrid::~CDualGrid() {}

CPoint::CPoint(unsigned short val_nDim, unsigned long val_globalindex,
               CConfig *config) : CDualGrid(val_nDim) {
  
  /*--- Element, point, and edge structure initialization ---*/
  
  Elem.clear();  nElem  = 0;
  Point.clear(); nPoint = 0;
  Edge.clear();
  
  /*--- Pointer initialization ---*/
  
  Volume    = NULL;
  vertex    = NULL;
  coord     = NULL;
  Coord_old = NULL;
  Coord_sum = NULL;
  
  /*--- Volume and coordinates of the control volume ---*/
  
  Volume = new double[1]; Volume[0] = 0.0;
  coord  = new double[nDim];
  
  /*--- Indicator if the control volume has been agglomerated (multigrid) ---*/
  
  Agglomerate = false;
  
  /*--- Indicator if the point is to be moved in a mesh deformation ---*/
  
  Move = true;
  
  /*--- Identify boundaries, physical boundaries (not send-receive
   condition), detect if an element belong to the domain or it must
   be computed with other processor  ---*/
  
  Boundary = false;
  PhysicalBoundary = false;
  Domain = true;
  RepeatedPoint = false;
  
  /*--- Set the global index of this node for parallel simulations ---*/
  
  GlobalIndex = val_globalindex;
  
  /*--- Set the color for mesh partitioning ---*/
  
  Color = 0;
  
  /*--- Arrays for smoothing the numerical grid coordinates ---*/
  
  if (config->GetSmoothNumGrid()) {
    Coord_old = new double[nDim];
    Coord_sum = new double[nDim];
  }
  
}

CPoint::CPoint(double val_coord_0, double val_coord_1,
               unsigned long val_globalindex, CConfig *config) : CDualGrid(2) {
  
  /*--- Element, point and edge structures initialization ---*/
  
  Elem.clear();  nElem = 0;
  Point.clear(); nPoint = 0;
  Edge.clear();
  
  /*--- Array initialization ---*/
  
  Volume    = NULL;
  vertex    = NULL;
  coord     = NULL;
  Coord_old = NULL;
  Coord_sum = NULL;
  
  /*--- Volume and coordinates of the control volume ---*/
  
  Volume = new double[1]; Volume[0] = 0.0;
  coord = new double[nDim]; coord[0] = val_coord_0; coord[1] = val_coord_1;
  
  /*--- Indicator if the control volume has been agglomerated (multigrid) ---*/
  
  Agglomerate = false;
  
  /*--- Indicator if the point is to be moved in a mesh deformation ---*/
  
  Move = true;
  
  /*--- Identify boundaries, physical boundaries (not send-receive
   condition), detect if an element belong to the domain or it must
   be computed with other processor  ---*/
  
  Boundary = false;
  PhysicalBoundary = false;
  Domain = true;
  RepeatedPoint = false;

  /*--- Set the color for mesh partitioning ---*/
  
  Color = 0;
  
  /*--- Set the global index of this node for parallel simulations ---*/
  
  GlobalIndex = val_globalindex;
  
  /*--- Arrays for smoothing the numerical grid coordinates ---*/
  
  if (config->GetSmoothNumGrid()) {
    Coord_old = new double[nDim];
    Coord_sum = new double[nDim];
  }
  
}

CPoint::CPoint(double val_coord_0, double val_coord_1, double val_coord_2,
               unsigned long val_globalindex, CConfig *config) : CDualGrid(3) {
  
  /*--- Element, point and edge structures initialization ---*/
  
  Elem.clear(); nElem = 0;
  Point.clear(); nPoint = 0;
  Edge.clear();
  
  /*--- Array initialization ---*/
  
  Volume    = NULL;
  vertex    = NULL;
  coord     = NULL;
  Coord_old = NULL;
  Coord_sum = NULL;
  
  /*--- Volume and coordinates of the control volume ---*/
  
  Volume = new double[1]; Volume[0] = 0.0;
  coord = new double[nDim];
  coord[0] = val_coord_0; coord[1] = val_coord_1; coord[2] = val_coord_2;
  
  /*--- Indicator if the control volume has been agglomerated (multigrid) ---*/
  
  Agglomerate = false;
  
  /*--- Indicator if the point is to be moved in a mesh deformation ---*/
  
  Move = true;
  
  /*--- Identify boundaries, physical boundaries (not send-receive
   condition), detect if an element belong to the domain or it must
   be computed with other processor  ---*/
  
  Boundary = false;
  PhysicalBoundary = false;
  Domain = true;
  RepeatedPoint = false;

  /*--- Set the color for mesh partitioning ---*/
  
  Color = 0;
  
  /*--- Set the global index of this node for parallel simulations ---*/
  
  GlobalIndex = val_globalindex;
  
  /*--- Arrays for smoothing the numerical grid coordinates ---*/
  
  if (config->GetSmoothNumGrid()) {
    Coord_old = new double[nDim];
    Coord_sum = new double[nDim];
  }
  
}

CPoint::~CPoint() {
  
  Elem.~vector();
  Point.~vector();
  Edge.~vector();
  Children_CV.~vector();
  
  if (Volume    != NULL) delete[] Volume;
  if (vertex    != NULL) delete[] vertex;
  if (coord     != NULL) delete[] coord;
  if (Coord_old != NULL) delete[] Coord_old;
  if (Coord_sum != NULL) delete[] Coord_sum;
  
}

void CPoint::SetPoint(unsigned long val_point) {
  
  unsigned short iPoint;
  bool new_point;
  
  /*--- Check if the point has already been processed to avoid duplicates ---*/
  
  new_point = true;
  for (iPoint = 0; iPoint < GetnPoint(); iPoint++)
    if (Point[iPoint] == val_point) {
      new_point = false;
      break;
    }
  
  /*--- If new, store the index value and update the edge structure ---*/
  
  if (new_point) {
    Point.push_back(val_point);
    Edge.push_back(-1);
    nPoint = Point.size();
  }
  
}

void CPoint::SetBoundary(unsigned short val_nmarker) {
  
  /*--- Check if the vertex has already been processed to avoid duplicates ---*/
  
  if (!Boundary) {
    
    /*--- Allocate & initialization the vertex structure (set to -1) ---*/
    
    vertex = new long[val_nmarker];
    for (unsigned short imarker = 0; imarker < val_nmarker; imarker++)
      vertex[imarker] = -1;
    
  }
  
  /*--- Set the flag to indicate that this is a boundary node ---*/
  
  Boundary = true;
  
}

CEdge::CEdge(unsigned long val_iPoint, unsigned long val_jPoint,
             unsigned short val_ndim) : CDualGrid(val_ndim) {
  
  /*--- Pointer initialization ---*/
  
  Coord_CG = NULL;
  Normal   = NULL;
  Nodes    = NULL;
  
  /*--- Allocate center of gravity coordinates, nodes, and face normal ---*/
  
  Coord_CG = new double [nDim];
  Nodes    = new unsigned long[2];
  Normal   = new double [nDim];
  
  /*--- Initialize the coordinates and normal to zero ---*/
  
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    Coord_CG[iDim] = 0.0;
    Normal[iDim]   = 0.0;
  }
  
  /*--- Initialize the node index values for this edge from the input ---*/
  
  Nodes[0] = val_iPoint;
  Nodes[1] = val_jPoint;
  
  /*--- Initialize the edge color for METIS+OpenMP to zero ---*/
  
  Color = 0;
  
}

CEdge::~CEdge() {
  
  if (Coord_CG != NULL) delete[] Coord_CG;
  if (Normal   != NULL) delete[] Normal;
  if (Nodes    != NULL) delete[] Nodes;
  
}

void CEdge::SetCG(double **val_coord) {
  
  unsigned short iDim, iNode;
  
  /*--- Compute the center of gravity for this edge as the average
   of the coordinates for the two nodes sharing the edge. ---*/
  
  for (iDim = 0; iDim < nDim; iDim++) {
    Coord_CG[iDim] = 0.0;
    for (iNode = 0; iNode < 2; iNode++)
      Coord_CG[iDim] += val_coord[iNode][iDim]/2.0;
  }
  
}

double CEdge::GetVolume(double *val_coord_Edge_CG, double *val_coord_FaceElem_CG,
                        double *val_coord_Elem_CG, double *val_coord_Point) {
  
  unsigned short iDim;
  double vec_a[3], vec_b[3], vec_c[3], vec_d[3], Local_Volume;
  
  /*--- Compute the volume (3-D) of the tetrahedron associated with a node,
   an edge, the CG of a neighboring element, and the CG of an adjacent face
   of that element. These local volumes will be summed to total the volume
   of the dual control volume. ---*/
  
  for (iDim = 0; iDim < nDim; iDim++) {
    vec_a[iDim] = val_coord_Edge_CG[iDim]-val_coord_Point[iDim];
    vec_b[iDim] = val_coord_FaceElem_CG[iDim]-val_coord_Point[iDim];
    vec_c[iDim] = val_coord_Elem_CG[iDim]-val_coord_Point[iDim];
  }
  
  vec_d[0] =   vec_a[1]*vec_b[2]-vec_a[2]*vec_b[1];
  vec_d[1] = -(vec_a[0]*vec_b[2]-vec_a[2]*vec_b[0]);
  vec_d[2] =   vec_a[0]*vec_b[1]-vec_a[1]*vec_b[0];
  
  Local_Volume = fabs(vec_c[0]*vec_d[0] + vec_c[1]*vec_d[1] + vec_c[2]*vec_d[2])/6.0;
  
  return Local_Volume;
}

double CEdge::GetVolume(double *val_coord_Edge_CG, double *val_coord_Elem_CG,
                        double *val_coord_Point) {
  
  unsigned short iDim;
  double vec_a[2], vec_b[2], Local_Volume;
  
  /*--- Compute the area (2-D) of the triangle associated with a node,
   an edge, and the CG of a neighboring element. These local areas
   will be summed to total the area of the dual control volume. ---*/
  
  for (iDim = 0; iDim < nDim; iDim++) {
    vec_a[iDim] = val_coord_Elem_CG[iDim]-val_coord_Point[iDim];
    vec_b[iDim] = val_coord_Edge_CG[iDim]-val_coord_Point[iDim];
  }
  
  Local_Volume = 0.5*fabs(vec_a[0]*vec_b[1]-vec_a[1]*vec_b[0]);
  
  return Local_Volume;
}

void CEdge::SetNodes_Coord(double *val_coord_Edge_CG, double *val_coord_FaceElem_CG,
                           double *val_coord_Elem_CG, CConfig *config) {
  
  unsigned short iDim;
  double vec_a[3], vec_b[3], Dim_Normal[3];
  
  /*--- Compute the normal vector for this dual control volume face.
   This is a face in the interior of the domain (3-D). ---*/
  
  for (iDim = 0; iDim < nDim; iDim++) {
    vec_a[iDim] = val_coord_Elem_CG[iDim]-val_coord_Edge_CG[iDim];
    vec_b[iDim] = val_coord_FaceElem_CG[iDim]-val_coord_Edge_CG[iDim];
  }
  
  Dim_Normal[0] =  0.5*(vec_a[1]*vec_b[2]-vec_a[2]*vec_b[1]);
  Dim_Normal[1] = -0.5*(vec_a[0]*vec_b[2]-vec_a[2]*vec_b[0]);
  Dim_Normal[2] =  0.5*(vec_a[0]*vec_b[1]-vec_a[1]*vec_b[0]);
  
  /*--- Add the normal contribution from this control volume face ---*/
  
  Normal[0] += Dim_Normal[0];
  Normal[1] += Dim_Normal[1];
  Normal[2] += Dim_Normal[2];
  
}

void CEdge::SetNodes_Coord(double *val_coord_Edge_CG, double *val_coord_Elem_CG,
                           CConfig *config) {
  
  double Dim_Normal[2];
  
  /*--- Compute the normal vector for this dual control volume face.
   This is a segment in the interior of the domain (2-D). ---*/
  
  Dim_Normal[0] =   val_coord_Elem_CG[1]-val_coord_Edge_CG[1];
  Dim_Normal[1] = -(val_coord_Elem_CG[0]-val_coord_Edge_CG[0]);
  
  /*--- Add the normal contribution from this control volume face ---*/
  
  Normal[0] += Dim_Normal[0];
  Normal[1] += Dim_Normal[1];
  
}

CVertex::CVertex(unsigned long val_point, unsigned short val_nDim) : CDualGrid(val_nDim) {
  
  /*--- Pointer initialization ---*/
  
  Nodes  = NULL;
  Normal = NULL;
  
  /*--- Allocate the node index and face normal arrays ---*/
  
  Nodes  = new unsigned long[1];
  Normal = new double [nDim];
  
  /*--- Initialize the node index and normal ---*/
  
  Nodes[0] = val_point;
  for (unsigned short iDim = 0; iDim < nDim; iDim ++) Normal[iDim] = 0.0;
  
  /*--- Set to zero the variation of the coordinates (used for tracking
   changes in the shape of the boundaries as part of mesh deformation) ---*/
  
  VarCoord[0] = 0.0; VarCoord[1] = 0.0; VarCoord[2] = 0.0;
  
}

CVertex::~CVertex() {
  
  if (Nodes  != NULL) delete[] Nodes;
  if (Normal != NULL) delete[] Normal;
  
}

void CVertex::SetNodes_Coord(double *val_coord_Edge_CG, double *val_coord_FaceElem_CG,
                             double *val_coord_Elem_CG, CConfig *config) {
  
  double vec_a[3], vec_b[3], Dim_Normal[3];
  unsigned short iDim;
  
  /*--- Compute the normal vector for this dual control volume face.
   This is a face along the boundary of the domain (3-D). ---*/
  
  for (iDim = 0; iDim < nDim; iDim++) {
    vec_a[iDim] = val_coord_Elem_CG[iDim]-val_coord_Edge_CG[iDim];
    vec_b[iDim] = val_coord_FaceElem_CG[iDim]-val_coord_Edge_CG[iDim];
  }
  
  Dim_Normal[0] =  0.5*(vec_a[1]*vec_b[2]-vec_a[2]*vec_b[1]);
  Dim_Normal[1] = -0.5*(vec_a[0]*vec_b[2]-vec_a[2]*vec_b[0]);
  Dim_Normal[2] =  0.5*(vec_a[0]*vec_b[1]-vec_a[1]*vec_b[0]);
  
  /*--- Add the normal contribution from this control volume face ---*/
  
  Normal[0] += Dim_Normal[0]; 
  Normal[1] += Dim_Normal[1];	
  Normal[2] += Dim_Normal[2];
  
}

void CVertex::SetNodes_Coord(double *val_coord_Edge_CG, double *val_coord_Elem_CG,
                             CConfig *config) {
  
  double Dim_Normal[2];
  
  /*--- Compute the normal vector for this dual control volume face. 
   This is a segment along the boundary of the domain (2-D). ---*/
  
  Dim_Normal[0] =   val_coord_Elem_CG[1]-val_coord_Edge_CG[1];
  Dim_Normal[1] = -(val_coord_Elem_CG[0]-val_coord_Edge_CG[0]);
  
  /*--- Add the normal contribution from this control volume face ---*/
  
  Normal[0] += Dim_Normal[0]; 
  Normal[1] += Dim_Normal[1];
  
}

void CVertex::AddNormal(double *val_face_normal) {
  
  /*--- Increment the local normal with the input vector ---*/
  
  Normal[0] += val_face_normal[0]; 
  Normal[1] += val_face_normal[1];
  if (nDim == 3) Normal[2] += val_face_normal[2];
  
}
