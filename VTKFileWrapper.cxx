/*
 * VTKReader.cxx
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

#include <stdlib.h>
#include <stdio.h>

#include <iostream>

using std::cout;
using std::endl;

#include <GMGW_unstr.hxx>
#include <GMGW_VTKFileWrapper.hxx>

VTKFileWrapper::VTKFileWrapper(const char fileNameBase[]) :
    FileWrapper()
{
  char fileNameIn[1000];
  snprintf(fileNameIn, 1000, "%s.vtk", fileNameBase);
  input = fopen(fileNameIn, "r");
  if (input == nullptr) {
    cout << "Bad file name " << fileNameIn << endl;
    exit(1);
  }
}

void
VTKFileWrapper::consumeLine() const
{
  GMGW_int result;
  result = fscanf(input, "%*[^\n]"); /* Skip to the End of the Line */
  result = fscanf(input, "%*1[\n]"); /* Skip One Newline */
  assert(result == 0);
}

void
VTKFileWrapper::skipHeader() const
{
  for (GMGW_int ii = 0; ii < 4; ii++) {
    consumeLine();
  }
}

void
VTKFileWrapper::scanFile()
{
  rewind(input);
  skipHeader();

  GMGW_int result = fscanf(input, " POINTS %" GMGW_int_format " float \n",
			   &nVerts);
  assert(result == 1);
  coordsStart = ftell(input);
  for (GMGW_int ii = 0; ii < nVerts; ii++) {
    double x, y, z;
    result = fscanf(input, " %lf %lf %lf \n", &x, &y, &z);
    assert(result == 3);
  }

  GMGW_int nInts = 0;
  result = fscanf(input, " CELLS %" GMGW_int_format " %" GMGW_int_format " \n",
		  &nCells, &nInts);
  connectStart = ftell(input);
  assert(result == 2);
  for (GMGW_int ii = 0; ii < nInts; ii++) {
    GMGW_int dummy;
    result = fscanf(input, " %d ", &dummy);
    assert(result == 1);
  }

  m_cellTypes = new char[nCells];
  m_isBdryFace = new bool[nCells];

  GMGW_int data;
  result = fscanf(input, "CELL_TYPES %u\n", &data);
  assert(result == 1);
  assert(data == nCells);
  for (GMGW_int ii = 0; ii < nCells; ii++) {
    char type;
    result = fscanf(input, " %hhd \n", &type);
    m_cellTypes[ii] = type;
    switch (type)
      {
      case BDRY_TRI:
	nBdryTris++;
	m_isBdryFace[ii] = true;
	break;
      case BDRY_QUAD:
	nBdryQuads++;
	m_isBdryFace[ii] = true;
	break;
      case TET:
	nTets++;
	m_isBdryFace[ii] = false;
	break;
      case HEX:
	nHexes++;
	m_isBdryFace[ii] = false;
	break;
      case PYRAMID:
	nPyrs++;
	m_isBdryFace[ii] = false;
	break;
      case PRISM:
	nPrisms++;
	m_isBdryFace[ii] = false;
	break;
      default:
	assert(0);
	break;
      }
  }
  assert(nBdryTris + nBdryQuads + nTets + nHexes + nPyrs + nPrisms == nCells);
  assert(
      nBdryTris * 4 + nBdryQuads * 5 + nTets * 5 + nHexes * 9 + nPyrs * 6
	  + nPrisms * 7 == nInts);

  identifyBdryVerts();

  writeMeshSizeInfo();
}

void
VTKFileWrapper::getNextVertexCoords(double& x, double& y, double& z) const
{
  GMGW_int result = fscanf(input, " %lf %lf %lf ", &x, &y, &z);
  assert(result == 3);
}

void
VTKFileWrapper::getNextCellConnectivity(GMGW_int& nConn,
GMGW_int connect[]) const
{
  GMGW_int result = fscanf(input, " %d ", &nConn);
  assert(result == 1);
  for (GMGW_int ii = 0; ii < nConn; ii++) {
    result = fscanf(input, " %u ", &(connect[ii]));
    assert(result == 1);
  }
}

