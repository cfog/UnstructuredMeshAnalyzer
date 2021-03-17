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

#include <iostream>

using std::cout;
using std::endl;

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "GMGW_unstr.hxx"
#include "GMGW_UGridFileWrapper.hxx"

UGridFileWrapper::UGridFileWrapper(const char fileNameBase[],
				   const char infix[]) :
    FileWrapper(), bigEndian(true), is64bit(false), whichCell(0)
{
  char fileNameIn[1000];
  snprintf(fileNameIn, 1000, "%s.%s.ugrid", fileNameBase, infix);
  input = fopen(fileNameIn, "r");
  if (input == nullptr) {
    cout << "Bad file name " << fileNameIn << endl;
    exit(1);
  }
  if (infix[0] == 'l') {
    cout << "Little-endian UGrid file" << endl;
    bigEndian = false; // Looks like little endian.
  }
  else {
    cout << "Big-endian UGrid file" << endl;
    bigEndian = true; // Looks like big endian.
  }
  if (infix[2] == 'l') {
    cout << "Reading a UGrid file with 64 bit ints" << endl;
#ifdef GMGW_INT32
    cout << "Reading and storing ints as 32 bit.  This won't end well; aborting." << endl;
    exit(3);
#endif
  }
  else {
    cout << "Reading a UGrid file with 32 bit ints" << endl;
#ifdef GMGW_INT64
    cout << "Reading and storing ints as 64 bit.  This won't end well; aborting." << endl;
    exit(3);
#endif

  }
}

GMGW_int
UGridFileWrapper::convertToInt(const unsigned char raw[]) const
{
#ifdef GMGW_INT32
  return convertToInt32(raw);
#else
  return convertToInt64(raw);
#endif
}

int32_t
UGridFileWrapper::convertToInt32(const unsigned char raw[4]) const
{
  // This implementation should be okay even with 64 bit ints, whereas a
  // union might be a little tricky.
//   return ( (reinterpret_cast<GMGW_int>(raw[0]) << 24)
//         + (reinterpret_cast<GMGW_int>(raw[1]) << 16)
//         + (reinterpret_cast<GMGW_int>(raw[2]) << 8)
//         + (raw[3]) );
  static const size_t is = sizeof(int32_t);
  union {
    unsigned char raw2[is];
    int32_t dummy;
  } a;
  // Linux x86_64 stores data in little-endian order.
  if (!bigEndian) {
    a.raw2[0] = raw[0];
    a.raw2[1] = raw[1];
    a.raw2[2] = raw[2];
    a.raw2[3] = raw[3];
  }
  else {
    a.raw2[is - 1] = raw[0];
    a.raw2[is - 2] = raw[1];
    a.raw2[is - 3] = raw[2];
    a.raw2[is - 4] = raw[3];
  }
  return a.dummy;
}

int64_t
UGridFileWrapper::convertToInt64(const unsigned char raw[8]) const
{
  // This implementation should be okay even with 64 bit ints, whereas a
  // union might be a little tricky.
//   return ( (reinterpret_cast<GMGW_int>(raw[0]) << 24)
//         + (reinterpret_cast<GMGW_int>(raw[1]) << 16)
//         + (reinterpret_cast<GMGW_int>(raw[2]) << 8)
//         + (raw[3]) );
  static const size_t is = sizeof(int64_t);
  union {
    unsigned char raw2[is];
    int64_t dummy;
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
    a.raw2[is - 1] = raw[0];
    a.raw2[is - 2] = raw[1];
    a.raw2[is - 3] = raw[2];
    a.raw2[is - 4] = raw[3];
    a.raw2[is - 5] = raw[4];
    a.raw2[is - 6] = raw[5];
    a.raw2[is - 7] = raw[6];
    a.raw2[is - 8] = raw[7];
  }
  return a.dummy;
}

double
UGridFileWrapper::convertToDouble(const unsigned char raw[8]) const
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

void
UGridFileWrapper::scanFile()
{
  // First, let's read the alleged sizes.
  unsigned char raw[36];
  size_t nRead = fread(raw, 1, 28, input);
  assert(nRead == 28);

  size_t intSize = sizeof(GMGW_int);
  nVerts = convertToInt(raw);
  nBdryTris = convertToInt(raw + intSize);
  nBdryQuads = convertToInt(raw + intSize * 2);
  nTets = convertToInt(raw + intSize * 3);
  nPyrs = convertToInt(raw + intSize * 4);
  nPrisms = convertToInt(raw + intSize * 5);
  nHexes = convertToInt(raw + intSize * 6);

  // Heuristic size check.  Far from foolproof, but should be right nearly
  // always, even if the file name is wrong.

  // Number of equivalent tets is:
  uint64_t numEquivTets = uint64_t(nTets) + 2 * uint64_t(nPyrs)
    + 3 * uint64_t(nPrisms) + 6 * uint64_t(nHexes);

  // Approximate range of number of tets should be:
  uint64_t minTets = 4.5 * uint64_t(nVerts);
  uint64_t maxTets = 6.5 * uint64_t(nVerts); // Some windage added here.

  if (!(nVerts<20000) && (maxTets < minTets || minTets > numEquivTets || numEquivTets > maxTets)) {
    cout << "Looks like this file has different endianness.  Trying again."
	<< endl;
    bigEndian = !bigEndian;
    nVerts = convertToInt(raw);
    nBdryTris = convertToInt(raw + intSize);
    nBdryQuads = convertToInt(raw + intSize * 2);
    nTets = convertToInt(raw + intSize * 3);
    nPyrs = convertToInt(raw + intSize * 4);
    nPrisms = convertToInt(raw + intSize * 5);
    nHexes = convertToInt(raw + intSize * 6);

    numEquivTets = nTets + 2 * nPyrs + 3 * nPrisms + 6 * nHexes;

    minTets = 4.5 * nVerts;
    maxTets = 6.5 * nVerts; // Some windage added here.

    if (!(nVerts<20000) && (maxTets < minTets || minTets>numEquivTets || numEquivTets > maxTets)) {
      cout << "Neither enddianness looks right.  I give up." << endl;
      cout << "Could also be a mismatch in integer size." << endl;
      exit(2);
    }
  }
  nCells = nTets + nPyrs + nPrisms + nHexes + nBdryTris + nBdryQuads;

  coordsStart = ftell(input);
  connectStart = coordsStart + nVerts * sizeof(double) * 3;

  m_cellTypes = new char[nCells];
  m_isBdryFace = new bool[nCells];

  GMGW_int iC = 0, end = nBdryTris;
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
  for (GMGW_int iV = 0; iV < nVerts; iV++) {
    m_isBdryVert[iV] = 0;
  }

  // This is faster than the default implementation, because with UGrid,
  // you know exactly where to find all of the bdry faces, so there's no
  // need to touch the rest.
  seekStartOfConnectivity();
  for (iC = 0; iC < nBdryTris + nBdryQuads; iC++) {
    GMGW_int nConn = -1;
    GMGW_int connect[8];
    getNextCellConnectivity(nConn, connect);
    assert(nConn == 3 || nConn == 4);
    assert(isBdryFace(iC));
    for (GMGW_int ii = 0; ii < nConn; ii++) {
      setBdryVert(connect[ii]);
    }
  }
  nBdryVerts = 0;
  for (GMGW_int iV = 0; iV < nVerts; iV++) {
    if (m_isBdryVert[iV])
      nBdryVerts++;
  }

  writeMeshSizeInfo();
}

void
UGridFileWrapper::seekStartOfConnectivity() const
{
  whichCell = 0;
  FileWrapper::seekStartOfConnectivity();
}

void
UGridFileWrapper::getNextCellConnectivity(GMGW_int& nConn,
					  GMGW_int connect[]) const
{
  switch (m_cellTypes[whichCell])
    {
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
  size_t nRead = fread(raw, 1, sizeof(GMGW_int) * nConn, input);
  assert(nRead == sizeof(GMGW_int) * nConn);
  for (GMGW_int ii = 0; ii < nConn; ii++) {
    connect[ii] = convertToInt(raw + sizeof(GMGW_int) * ii) - 1;
  }
  if (nConn == 5) {
    // This is a pyramid, so the order of data has to be changed.
    GMGW_int temp = connect[1];
    connect[1] = connect[3];
    connect[3] = temp;

    temp = connect[2];
    connect[2] = connect[4];
    connect[4] = temp;
  }

  whichCell++;
  if (whichCell == nBdryTris + nBdryQuads) {
    // Need to seek over the bdry conditions
    size_t offset = sizeof(GMGW_int) * (nBdryTris + nBdryQuads);
    fseek(input, offset, SEEK_CUR);
  }
}

void
UGridFileWrapper::getNextVertexCoords(double& x, double& y, double& z) const
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
