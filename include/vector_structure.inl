/*!
 * \file vector_structure.inl
 * \brief inline subroutines of the <i>vector_structure.hpp</i> file.
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

#pragma once

inline void CSysVector::SetValZero(void) { 
  for (unsigned long i = 0; i < nElm; i++)
		vec_val[i] = 0.0;
}

inline unsigned long CSysVector::GetLocSize() const { return nElm; }

inline unsigned long CSysVector::GetSize() const { return (unsigned long)nElm; }

inline unsigned short CSysVector::GetNVar() const { return nVar; }

inline unsigned long CSysVector::GetNBlk() const { return nBlk; }

inline unsigned long CSysVector::GetNBlkDomain() const { return nBlkDomain; }

inline double & CSysVector::operator[](const unsigned long & i) { return vec_val[i]; }

inline const double & CSysVector::operator[](const unsigned long & i) const { return vec_val[i]; }
