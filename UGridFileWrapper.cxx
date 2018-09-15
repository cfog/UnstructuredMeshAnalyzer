/*
 * UGridFileWrapper.cxx
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
#include <string.h>

#include "GMGW_unstr.hxx"
#include "GMGW_UGridFileWrapper.hxx"

UGridFileWrapper::UGridFileWrapper(const char fileNameBase[],
                                   const char infix[]) :
    FileWrapper(), bigEndian(true), whichCell(0)
{
    char fileNameIn[1000];
    snprintf(fileNameIn, 1000, "%s.%s.ugrid", fileNameBase, infix);
    input = fopen(fileNameIn, "r");
    if (input == nullptr) {
        printf("Bad file name %s.\n", fileNameIn);
        exit(1);
    }
    const char* result = strstr(infix, "l");
    if (result) {
        printf("Little-endian UGrid file\n");
        bigEndian = false; // Looks like little endian.
    }
    else {
        printf("Big-endian UGrid file\n");
        bigEndian = true; // Looks like big endian.
    }
}

unsigned int UGridFileWrapper::convertToInt(const unsigned char raw[4]) const
{
  // This implementation should be okay even with 64 bit ints, whereas a
  // union might be a little tricky.
//   return ( (reinterpret_cast<int>(raw[0]) << 24)
//         + (reinterpret_cast<int>(raw[1]) << 16)
//         + (reinterpret_cast<int>(raw[2]) << 8)
//         + (raw[3]) );
  static const int is = sizeof(int);
  union {
    unsigned char raw2[is];
    unsigned int dummy;
  } a;
  // Linux x86_64 stores data in little-endian order.
  if (!bigEndian) {
      a.raw2[0] = raw[0];
      a.raw2[1] = raw[1];
      a.raw2[2] = raw[2];
      a.raw2[3] = raw[3];
  }
  else {
      a.raw2[is-1] = raw[0];
      a.raw2[is-2] = raw[1];
      a.raw2[is-3] = raw[2];
      a.raw2[is-4] = raw[3];
  }
  return a.dummy;
}

double UGridFileWrapper::convertToDouble(const unsigned char raw[8]) const
{
  union {
    unsigned char raw2[8];
    double dummy;
  } a;
  // Linux x86_64 stores data in little-endian order.
  if (!bigEndian) {
      a.raw2[0] = raw[0];
      a.raw2[1] = raw[1];
      a.raw2[2] = raw[2];
      a.raw2[3] = raw[3];
      a.raw2[4] = raw[4];
      a.raw2[5] = raw[5];
      a.raw2[6] = raw[6];
      a.raw2[7] = raw[7];
  }
  else {
      a.raw2[7] = raw[0];
      a.raw2[6] = raw[1];
      a.raw2[5] = raw[2];
      a.raw2[4] = raw[3];
      a.raw2[3] = raw[4];
      a.raw2[2] = raw[5];
      a.raw2[1] = raw[6];
      a.raw2[0] = raw[7];
  }
  return a.dummy;
}


void UGridFileWrapper::scanFile()
{
    // First, let's read the alleged sizes.
    unsigned char raw[36];
    size_t nRead = fread(raw, 1, 28, input);
    assert(nRead == 28);

    nVerts     = convertToInt(raw);
    nBdryTris  = convertToInt(raw + 4);
    nBdryQuads = convertToInt(raw + 8);
    nTets      = convertToInt(raw + 12);
    nPyrs      = convertToInt(raw + 16);
    nPrisms    = convertToInt(raw + 20);
    nHexes     = convertToInt(raw + 24);

    // Heuristic size check.  Far from foolproof, but should be right nearly
    // always, even if the file name is wrong.

    // Number of equivalent tets is:
    unsigned long numEquivTets = nTets + 2*nPyrs + 3*nPrisms + 6*nHexes;

    // Approximate range of number of tets should be:
    unsigned long minTets = 5*nVerts;
    unsigned long maxTets = 6.5*nVerts; // Some windage added here.

    if (maxTets < minTets ||
            numEquivTets < minTets || numEquivTets > maxTets) {
        printf("Looks like this file has different endianness.  Trying again.\n");
        bigEndian = !bigEndian;
        nVerts     = convertToInt(raw);
        nBdryTris  = convertToInt(raw + 4);
        nBdryQuads = convertToInt(raw + 8);
        nTets      = convertToInt(raw + 12);
        nPyrs      = convertToInt(raw + 16);
        nPrisms    = convertToInt(raw + 20);
        nHexes     = convertToInt(raw + 24);

        numEquivTets = nTets + 2*nPyrs + 3*nPrisms + 6*nHexes;

        minTets = 5*nVerts;
        maxTets = 6.5*nVerts; // Some windage added here.

        if (maxTets < minTets || numEquivTets < minTets ||
                numEquivTets > maxTets) {
            printf("Neither enddianness looks right.  I give up.\n");
        }
    }
    nCells = nTets + nPyrs + nPrisms + nHexes + nBdryTris + nBdryQuads;

    coordsStart = ftell(input);
    connectStart = coordsStart + nVerts*sizeof(double)*3;

    m_cellTypes = new char[nCells];
    m_isBdryFace = new bool[nCells];

    unsigned long int iC = 0, end = nBdryTris;
    for (; iC < end; iC++) {
        m_cellTypes[iC] = BDRY_TRI;
        m_isBdryFace[iC] = true;
    }
    end += nBdryQuads;
    for (; iC < end; iC++) {
        m_cellTypes[iC] = BDRY_QUAD;
        m_isBdryFace[iC] = true;
    }
    end += nTets;
    for (; iC < end; iC++) {
        m_cellTypes[iC] = TET;
        m_isBdryFace[iC] = false;
    }
    end += nPyrs;
    for (; iC < end; iC++) {
        m_cellTypes[iC] = PYRAMID;
        m_isBdryFace[iC] = false;
    }
    end += nPrisms;
    for (; iC < end; iC++) {
        m_cellTypes[iC] = PRISM;
        m_isBdryFace[iC] = false;
    }
    end += nHexes;
    for (; iC < end; iC++) {
        m_cellTypes[iC] = HEX;
        m_isBdryFace[iC] = false;
    }

    m_isBdryVert = new bool[nVerts];
    for (unsigned long int iV = 0; iV < nVerts; iV++) {
        m_isBdryVert[iV] = 0;
    }

    // This is faster than the default implementation, because with UGrid,
    // you know exactly where to find all of the bdry faces, so there's no
    // need to touch the rest.
    seekStartOfConnectivity();
    for (iC = 0; iC < nBdryTris + nBdryQuads; iC++) {
        int nConn = -1;
        unsigned int connect[8];
        getNextCellConnectivity(nConn, connect);
        assert(nConn == 3 || nConn == 4);
        assert(isBdryFace(iC));
        for (int ii = 0; ii < nConn; ii++) {
            setBdryVert(connect[ii]);
        }
    }
    nBdryVerts = 0;
    for (unsigned long int iV = 0; iV < nVerts; iV++) {
        if (m_isBdryVert[iV]) nBdryVerts++;
    }

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

void UGridFileWrapper::seekStartOfConnectivity() const
{
    whichCell = 0;
    FileWrapper::seekStartOfConnectivity();
}

void UGridFileWrapper::getNextCellConnectivity(int& nConn, unsigned int connect[]) const
{
    switch (m_cellTypes[whichCell]) {
        case BDRY_TRI:
            nConn = 3;
            break;
        case BDRY_QUAD:
        case TET:
            nConn = 4;
            break;
        case PYRAMID:
            nConn = 5;
            break;
        case PRISM:
            nConn = 6;
            break;
        case HEX:
            nConn = 8;
            break;
        default:
            nConn = -1;
            assert(0);
            break;
    }

    // Now actually read the data and convert it.
    unsigned char raw[256];
    size_t nRead = fread(raw, 1, sizeof(int)*nConn, input);
    assert(nRead == sizeof(int)*nConn);
    for (int ii = 0; ii < nConn; ii++) {
        connect[ii] = convertToInt(raw + sizeof(int)*ii) - 1;
    }
    if (nConn == 5) {
        // This is a pyramid, so the order of data has to be changed.
        unsigned int temp = connect[1];
        connect[1] = connect[3];
        connect[3] = temp;

        temp = connect[2];
        connect[2] = connect[4];
        connect[4] = temp;
    }

    whichCell++;
    if (whichCell == nBdryTris + nBdryQuads) {
        // Need to seek over the bdry conditions
        long offset = sizeof(int) * (nBdryTris + nBdryQuads);
        fseek(input, offset, SEEK_CUR);
    }
}

void UGridFileWrapper::getNextVertexCoords(double& x, double& y, double& z) const
{
    size_t nRead = 0;
    unsigned char raw[8];
    nRead = fread(raw, 1, sizeof(double), input);
    assert(nRead == 8);
    x = convertToDouble(raw);

    nRead = fread(raw, 1, sizeof(double), input);
    assert(nRead == 8);
    y = convertToDouble(raw);

    nRead = fread(raw, 1, sizeof(double), input);
    assert(nRead == 8);
    z = convertToDouble(raw);
}
