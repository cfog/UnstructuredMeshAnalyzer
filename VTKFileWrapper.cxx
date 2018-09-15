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

#include <GMGW_unstr.hxx>
#include <GMGW_VTKFileWrapper.hxx>

VTKFileWrapper::VTKFileWrapper(const char fileNameBase[]) :
    FileWrapper()
{
    char fileNameIn[1000];
    snprintf(fileNameIn, 1000, "%s.vtk", fileNameBase);
    input = fopen(fileNameIn, "r");
    if (input == nullptr) {
        printf("Bad file name %s.\n", fileNameIn);
        exit(1);
    }
}


void VTKFileWrapper::consumeLine() const
{
    int result;
    result = fscanf(input, "%*[^\n]");   /* Skip to the End of the Line */
    result = fscanf(input, "%*1[\n]");   /* Skip One Newline */
    assert(result == 0);
}

void VTKFileWrapper::skipHeader() const
{
    for (int ii = 0; ii < 4; ii++) {
        consumeLine();
    }
}

void VTKFileWrapper::scanFile()
{
    rewind(input);
    skipHeader();

    int result = fscanf(input, " POINTS %lu float \n", &nVerts);
    assert(result == 1);
    coordsStart = ftell(input);
    for (unsigned int ii = 0; ii < nVerts; ii++) {
        double x, y, z;
        result = fscanf(input, " %lf %lf %lf \n", &x, &y, &z);
        assert(result == 3);
    }

    unsigned long nInts = 0;
    result = fscanf(input, " CELLS %lu %lu \n", &nCells, &nInts);
    connectStart = ftell(input);
    assert(result == 2);
    for (unsigned int ii = 0; ii < nInts; ii++) {
        int dummy;
        result = fscanf(input, " %d ", &dummy);
        assert(result == 1);
    }

    m_cellTypes = new char[nCells];
    m_isBdryFace = new bool[nCells];

    unsigned int data;
    result = fscanf(input, "CELL_TYPES %u\n", &data);
    assert(result == 1);
    assert(data == nCells);
    for (unsigned int ii = 0; ii < nCells; ii++) {
        char type;
        result = fscanf(input, " %hhd \n", &type);
        m_cellTypes[ii] = type;
        switch (type) {
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
    assert(nBdryTris*4 + nBdryQuads*5 + nTets*5 + nHexes*9 + nPyrs*6 + nPrisms*7
           == nInts);

    identifyBdryVerts();

    printf("Scanned mesh and found:\n");
    printf("%'13lu verts\n", nVerts);
    printf("%'13lu bdry verts\n", nBdryVerts);
    printf("%'13lu bdry tris\n", nBdryTris);
    printf("%'13lu bdry quads\n", nBdryQuads);
    printf("%'13lu tets\n", nTets);
    printf("%'13lu pyramids\n", nPyrs);
    printf("%'13lu prisms\n", nPrisms);
    printf("%'13lu hexes\n", nHexes);
    printf("%'13lu total cells\n", nCells);
}

void VTKFileWrapper::getNextVertexCoords(double& x, double& y, double& z) const
{
    int result = fscanf(input, " %lf %lf %lf ", &x, &y, &z);
    assert(result == 3);
}

void VTKFileWrapper::getNextCellConnectivity(int& nConn, unsigned int connect[]) const
{
    int result = fscanf(input, " %d ", &nConn);
    assert(result == 1);
    for (int ii = 0; ii < nConn; ii++) {
        result = fscanf(input, " %u ", &(connect[ii]));
        assert(result == 1);
    }
}

