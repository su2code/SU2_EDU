# SU^2 global Makefile.am
# Stanford University Aerospace Design Lab
#
# Created: 		10/26/12
# Last updated: 	10/30/13
#
# Created by: 		M. Colonno
# Updated by: 		T. Economon
#
##########################

ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS} -I m4

# build third-party optional dependencies first
SUBDIRS = externals

bin_PROGRAMS = bin/SU2_EDU

bin_SU2_EDU_SOURCES = include/definition_structure.hpp \
		include/integration_structure.hpp \
		include/integration_structure.inl \
		include/iteration_structure.hpp \
		include/numerics_structure.hpp \
		include/numerics_structure.inl \
		include/output_structure.hpp \
		include/solver_structure.hpp \
		include/solver_structure.inl \
		include/SU2_EDU.hpp \
		include/variable_structure.hpp \
		include/variable_structure.inl \
		src/definition_structure.cpp \
		src/integration_structure.cpp \
		src/integration_time.cpp \
		src/iteration_structure.cpp \
		src/numerics_direct_mean.cpp \
		src/numerics_direct_turbulent.cpp \
		src/numerics_structure.cpp \
		src/output_structure.cpp \
		src/output_tecplot.cpp \
		src/output_paraview.cpp \
		src/solver_direct_mean.cpp \
		src/solver_direct_turbulent.cpp \
		src/solver_structure.cpp \
		src/SU2_EDU.cpp \
		src/variable_direct_mean.cpp \
		src/variable_direct_turbulent.cpp \
		src/config_structure.cpp \
		src/dual_grid_structure.cpp \
		src/grid_movement_structure.cpp \
		src/geometry_structure.cpp \
		src/matrix_structure.cpp \
		src/primal_grid_structure.cpp \
		src/linear_solvers_structure.cpp \
		src/variable_structure.cpp \
		src/vector_structure.cpp

# always link to built dependencies from ./externals
bin_SU2_EDU_CXXFLAGS = @su2_externals_INCLUDES@
bin_SU2_EDU_LDADD = @su2_externals_LIBS@
