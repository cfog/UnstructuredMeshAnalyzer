/*
 * GMGW_unstr.hxx
 *
 *  Created on: Jan 27, 2017
 *      Author: cfog
 */


/*
    Copyright (C) The University of British Columbia, 2018.

    This file is part of UnstructuredMeshAnalyzer.

    UnstructuredMeshAnalyzer is free software: you can redistribute it
    and/or modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation, either version 3 of
    the License, or (at your option) any later version.

    UnstructuredMeshAnalyzer is distributed in the hope that it will be
    useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with UnstructuredMeshAnalyzer.  If not, see
    <https://www.gnu.org/licenses/>.
*/

#ifndef GMGW_UNSTR_HXX_
#define GMGW_UNSTR_HXX_

#include <string>
#include <GMGW_FileWrapper.hxx>

#define MIN(a,b) ((a) < (b) ? (a) : (b))

// VTK cell type definitions

#define BDRY_TRI 5
#define BDRY_QUAD 9
#define TET 10
#define PYRAMID 14
#define PRISM 13
#define HEX 12

int projectionChecks(const int nBdryVerts,
                     const double bdryCoords[][3], double bdryDist[],
                     int bdrySurf[], std::string nmbFileName);

void buildFaceList(FileWrapper* reader,
                   unsigned int faceToCell[][2],
                   unsigned int& badTris,
                   unsigned int& badQuads);

#endif /* GMGW_UNSTR_HXX_ */
