/*!
 * \file grid_movement_structure.hpp
 * \brief Headers of the main subroutines for doing the numerical grid 
 *        movement (including volumetric movement, surface movement and Free From 
 *        technique definition). The subroutines and functions are in 
 *        the <i>grid_movement_structure.cpp</i> file.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 1.2.0
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

#pragma once

#include "geometry_structure.hpp"
#include "config_structure.hpp"
#include "matrix_structure.hpp"
#include "vector_structure.hpp"
#include "linear_solvers_structure.hpp"

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <ctime>

using namespace std;

/*!
 * \class CGridMovement
 * \brief Class for moving the surface and volumetric 
 *        numerical grid (2D and 3D problems).
 * \author F. Palacios.
 * \version 1.2.0
 */
class CGridMovement {
public:

	/*! 
	 * \brief Constructor of the class. 
	 */
	CGridMovement(void);

	/*! 
	 * \brief Destructor of the class. 
	 */
	~CGridMovement(void);
  
  
  /*!
	 * \brief A pure virtual member.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	virtual void SetSurface_Deformation(CGeometry *geometry, CConfig *config);
  
};

/*! 
 * \class CVolumetricMovement
 * \brief Class for moving the volumetric numerical grid.
 * \author F. Palacios, A. Bueno, T. Economon, S. Padron.
 * \version 1.2.0
 */
class CVolumetricMovement : public CGridMovement {
protected:

	unsigned short nDim;		/*!< \brief Number of dimensions. */
	unsigned short nVar;		/*!< \brief Number of variables. */
  
	unsigned long nPoint;		/*!< \brief Number of points. */
	unsigned long nPointDomain;		/*!< \brief Number of points in the domain. */

  CSysMatrix StiffMatrix; /*!< \brief Matrix to store the point-to-point stiffness. */
  CSysVector LinSysSol;
  CSysVector LinSysRes;

public:

	/*! 
	 * \brief Constructor of the class.
	 */
	CVolumetricMovement(CGeometry *geometry);
	
	/*! 
	 * \brief Destructor of the class. 
	 */
	~CVolumetricMovement(void);
  
	/*!
	 * \brief Update the value of the coordinates after the grid movement.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void UpdateGridCoord(CGeometry *geometry, CConfig *config);
  
  /*!
	 * \brief Update the dual grid after the grid movement (edges and control volumes).
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void UpdateDualGrid(CGeometry *geometry, CConfig *config);
  
	/*! 
	 * \brief Update the coarse multigrid levels after the grid movement.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void UpdateMultiGrid(CGeometry **geometry, CConfig *config);
  
  /*!
	 * \brief Compute the stiffness matrix for grid deformation using spring analogy.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \return Value of the length of the smallest edge of the grid.
	 */
	double SetFEAMethodContributions_Elem(CGeometry *geometry, CConfig *config);
  
  /*!
	 * \brief Build the stiffness matrix for a 3-D hexahedron element. The result will be placed in StiffMatrix_Elem.
	 * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] StiffMatrix_Elem - Element stiffness matrix to be filled.
	 * \param[in] CoordCorners[8][3] - Index value for Node 1 of the current hexahedron.
	 */
  void SetFEA_StiffMatrix3D(CGeometry *geometry, CConfig *config, double **StiffMatrix_Elem, unsigned long PointCorners[8], double CoordCorners[8][3], unsigned short nNodes, double scale);
	
  /*!
	 * \brief Build the stiffness matrix for a 3-D hexahedron element. The result will be placed in StiffMatrix_Elem.
	 * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] StiffMatrix_Elem - Element stiffness matrix to be filled.
	 * \param[in] CoordCorners[8][3] - Index value for Node 1 of the current hexahedron.
	 */
  void SetFEA_StiffMatrix2D(CGeometry *geometry, CConfig *config, double **StiffMatrix_Elem, unsigned long PointCorners[8], double CoordCorners[8][3], unsigned short nNodes, double scale);
  
  /*!
	 * \brief Shape functions and derivative of the shape functions
   * \param[in] Xi - Local coordinates.
   * \param[in] Eta - Local coordinates.
   * \param[in] Mu - Local coordinates.
	 * \param[in] CoordCorners[8][3] - Coordiantes of the corners.
   * \param[in] shp[8][4] - Shape function information
	 */
  double ShapeFunc_Hexa(double Xi, double Eta, double Mu, double CoordCorners[8][3], double DShapeFunction[8][4]);
  
  /*!
	 * \brief Shape functions and derivative of the shape functions
   * \param[in] Xi - Local coordinates.
   * \param[in] Eta - Local coordinates.
   * \param[in] Mu - Local coordinates.
	 * \param[in] CoordCorners[8][3] - Coordiantes of the corners.
   * \param[in] shp[8][4] - Shape function information
	 */
  double ShapeFunc_Tetra(double Xi, double Eta, double Mu, double CoordCorners[8][3], double DShapeFunction[8][4]);
  
  /*!
	 * \brief Shape functions and derivative of the shape functions
   * \param[in] Xi - Local coordinates.
   * \param[in] Eta - Local coordinates.
   * \param[in] Mu - Local coordinates.
	 * \param[in] CoordCorners[8][3] - Coordiantes of the corners.
   * \param[in] shp[8][4] - Shape function information
	 */
  double ShapeFunc_Pyram(double Xi, double Eta, double Mu, double CoordCorners[8][3], double DShapeFunction[8][4]);
  
  /*!
	 * \brief Shape functions and derivative of the shape functions
   * \param[in] Xi - Local coordinates.
   * \param[in] Eta - Local coordinates.
   * \param[in] Mu - Local coordinates.
	 * \param[in] CoordCorners[8][3] - Coordiantes of the corners.
   * \param[in] shp[8][4] - Shape function information
	 */
  double ShapeFunc_Wedge(double Xi, double Eta, double Mu, double CoordCorners[8][3], double DShapeFunction[8][4]);
  
  /*!
	 * \brief Shape functions and derivative of the shape functions
   * \param[in] Xi - Local coordinates.
   * \param[in] Eta - Local coordinates.
   * \param[in] Mu - Local coordinates.
	 * \param[in] CoordCorners[8][3] - Coordiantes of the corners.
   * \param[in] shp[8][4] - Shape function information
	 */
  double ShapeFunc_Triangle(double Xi, double Eta, double CoordCorners[8][3], double DShapeFunction[8][4]);
  
  /*!
	 * \brief Shape functions and derivative of the shape functions
   * \param[in] Xi - Local coordinates.
   * \param[in] Eta - Local coordinates.
   * \param[in] Mu - Local coordinates.
	 * \param[in] CoordCorners[8][3] - Coordiantes of the corners.
   * \param[in] shp[8][4] - Shape function information
	 */
  double ShapeFunc_Rectangle(double Xi, double Eta, double CoordCorners[8][3], double DShapeFunction[8][4]);
  
  /*!
	 * \brief Compute the shape functions for hexahedron
	 * \param[in] HexaCorners[8][3] - coordinates of the cornes of the hexahedron.
	 */
  double GetHexa_Volume(double CoordCorners[8][3]);
  
  /*!
	 * \brief Compute the shape functions for hexahedron
	 * \param[in] TetCorners[4][3] - coordinates of the cornes of the hexahedron.
	 */
  double GetTetra_Volume(double CoordCorners[8][3]);
  
  /*!
	 * \brief Compute the shape functions for hexahedron
	 * \param[in] TetCorners[4][3] - coordinates of the cornes of the hexahedron.
	 */
  double GetWedge_Volume(double CoordCorners[8][3]);
  
  /*!
	 * \brief Compute the shape functions for hexahedron
	 * \param[in] TetCorners[4][3] - coordinates of the cornes of the hexahedron.
	 */
  double GetPyram_Volume(double CoordCorners[8][3]);
  
  /*!
	 * \brief Compute the shape functions for hexahedron
	 * \param[in] TetCorners[4][3] - coordinates of the cornes of the hexahedron.
	 */
  double GetTriangle_Area(double CoordCorners[8][3]);
  
  /*!
	 * \brief Compute the shape functions for hexahedron
	 * \param[in] TetCorners[4][3] - coordinates of the cornes of the hexahedron.
	 */
  double GetRectangle_Area(double CoordCorners[8][3]);
    
  /*!
	 * \brief Add the stiffness matrix for a 2-D triangular element to the global stiffness matrix for the entire mesh (node-based).
	 * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] StiffMatrix_Elem - Element stiffness matrix to be filled.
	 * \param[in] val_Point_0 - Index value for Node 0 of the current tetrahedron.
   * \param[in] val_Point_1 - Index value for Node 1 of the current tetrahedron.
   * \param[in] val_Point_2 - Index value for Node 2 of the current tetrahedron.
   * \param[in] val_Point_3 - Index value for Node 3 of the current tetrahedron.
	 */
  void AddFEA_StiffMatrix(CGeometry *geometry, double **StiffMatrix_Elem, unsigned long PointCorners[8], unsigned short nNodes);
  
  /*!
	 * \brief Check for negative volumes (all elements) after performing grid deformation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 */
	double Check_Grid(CGeometry *geometry);
  
	/*!
	 * \brief Check the boundary vertex that are going to be moved.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetBoundaryDisplacements(CGeometry *geometry, CConfig *config);
  
  /*!
	 * \brief Grid deformation using the spring analogy method.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 * \param[in] UpdateGeo - Update geometry.
	 */
	void SetVolume_Deformation(CGeometry *geometry, CConfig *config, bool UpdateGeo);
  
  /*!
	 * \brief Compute the determinant of a 3 by 3 matrix.
	 * \param[in] val_matrix 3 by 3 matrix.
	 * \result Determinant of the matrix
	 */
	double Determinant_3x3(double A00, double A01, double A02, double A10, double A11, double A12, double A20, double A21, double A22);

};

/*! 
 * \class CSurfaceMovement
 * \brief Class for moving the surface numerical grid.
 * \author F. Palacios, T. Economon.
 * \version 1.2.0
 */
class CSurfaceMovement : public CGridMovement {
protected:

public:
	
	/*! 
	 * \brief Constructor of the class.
	 */
	CSurfaceMovement(void);
	
	/*! 
	 * \brief Destructor of the class. 
	 */
	~CSurfaceMovement(void);
  
  /*!
	 * \brief Set a obstacle in a channel.
	 * \param[in] boundary - Geometry of the boundary.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetAirfoil(CGeometry *boundary, CConfig *config);
	
	/*! 
	 * \brief Copy the boundary coordinates to each vertex.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void CopyBoundary(CGeometry *geometry, CConfig *config);
  
  /*!
	 * \brief Set the surface/boundary deformation.
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	void SetSurface_Deformation(CGeometry *geometry, CConfig *config);
		
};

#include "grid_movement_structure.inl"
