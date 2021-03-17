/*
 * FileWrapper.cxx
 *
 *  Created on: Jan 30, 2017
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <iostream>

using std::cout;
using std::endl;

#include "config.h"
#include "GMGW_FileWrapper.hxx"
#include "GMGW_UGridFileWrapper.hxx"
#include "GMGW_VTKFileWrapper.hxx"
#ifdef HAVE_TAU
#include "GMGW_TAUFileWrapper.hxx"
#endif

extern "C" {
  void libMeshIO_is_present(void) {
  }    
}

FileWrapper::FileWrapper() :
    input(nullptr), coordsStart(-1), connectStart(-1), nVerts(0), nCells(0), nBdryTris(
	0), nBdryQuads(0), nTets(0), nPyrs(0), nPrisms(0), nHexes(0), nBdryVerts(
	0), m_cellTypes(nullptr), m_isBdryFace(nullptr), m_isBdryVert(nullptr)
{
}

FileWrapper::~FileWrapper()
{
  if (m_cellTypes)
    delete[] m_cellTypes;
  if (m_isBdryFace)
    delete[] m_isBdryFace;
  if (m_isBdryVert)
    delete[] m_isBdryVert;
  if (input)
    fclose(input);
}

void
FileWrapper::seekStartOfCoords() const
{
  fseek(input, coordsStart, SEEK_SET);
}

void
FileWrapper::seekStartOfConnectivity() const
{
  fseek(input, connectStart, SEEK_SET);
}

void
FileWrapper::identifyBdryVerts()
{
  m_isBdryVert = new bool[nVerts];
  for (GMGW_int ii = 0; ii < nVerts; ii++)
    clearBdryVert(ii);
  seekStartOfConnectivity();
  // Now process all the connectivity info.
  for (GMGW_int iC = 0; iC < nCells; iC++) {
    GMGW_int nConn = -1;
    GMGW_int connect[8];
    getNextCellConnectivity(nConn, connect);
    if (isBdryFace(iC)) {
      for (GMGW_int ii = 0; ii < nConn; ii++) {
	setBdryVert(connect[ii]);
      }
    }
  }

  // Finally, count the number of bdry verts
  for (GMGW_int ii = 0; ii < nVerts; ii++) {
    if (isBdryVert(ii)) {
      nBdryVerts++;
    }
  }
}

FileWrapper*
FileWrapper::factory(const char baseName[], const char type[],
		     const char ugridInfix[])
{
  if (strstr(type, "vtk")) {
    return new VTKFileWrapper(baseName);
  }
  else if (strstr(type, "ugrid")) {
    return new UGridFileWrapper(baseName, ugridInfix);
  }
#ifdef HAVE_TAU
  else if (strstr(type, "tau")) {
    return new TAUFileWrapper(baseName);
  }
#endif
  else {
    fprintf(stderr, "Missing or invalid file type: %s\n", type);
    exit(1);
  }
}


void
FileWrapper::writeMeshSizeInfo()
{
  cout << "Scanned mesh and found:" << endl;
  cout.width(16);
  cout << nVerts << " verts" << endl;
  cout.width(16);
  cout << nBdryVerts << " bdry verts" << endl;
  cout.width(16);
  cout << nBdryTris << " bdry tris" << endl;
  cout.width(16);
  cout << nBdryQuads << " bdry quads" << endl;
  cout.width(16);
  cout << nTets << " tets" << endl;
  cout.width(16);
  cout << nPyrs << " pyramids" << endl;
  cout.width(16);
  cout << nPrisms << " prisms" << endl;
  cout.width(16);
  cout << nHexes << " hexes" << endl;
  cout.width(16);
  cout << nCells << " total cells and bdry faces" << endl;
}

