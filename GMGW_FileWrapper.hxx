/*
 * GMGW_FileWrapper.hxx
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

#ifndef GMGW_FILEWRAPPER_HXX_
#define GMGW_FILEWRAPPER_HXX_

#include "config.h"

#include <cstdio>
#include <assert.h>

class FileWrapper {
protected:
  FILE* input;
  size_t coordsStart, connectStart;
  GMGW_int nVerts, nCells, nBdryTris, nBdryQuads, nTets, nPyrs, nPrisms, nHexes,
      nBdryVerts;
  char* m_cellTypes;
  bool *m_isBdryFace, *m_isBdryVert;
public:
  FileWrapper();
  virtual
  ~FileWrapper();

  static FileWrapper*
  factory(const char baseName[], const char type[], const char ugridInfix[]);

  void
  writeMeshSizeInfo();

  virtual void
  scanFile() = 0;

  char
  getCellType(GMGW_int i) const
  {
    assert(i < nCells);
    return m_cellTypes[i];
  }

  bool
  isBdryFace(GMGW_int i) const
  {
    assert(i < nCells);
    return m_isBdryFace[i];
  }

  bool
  isBdryVert(GMGW_int i) const
  {
    assert(i < nVerts);
    return m_isBdryVert[i];
  }

  void
  setBdryVert(GMGW_int i) const
  {
    assert(i < nVerts);
    m_isBdryVert[i] = true;
  }

  void
  clearBdryVert(GMGW_int i) const
  {
    assert(i < nVerts);
    m_isBdryVert[i] = false;
  }

  GMGW_int
  getNumBdryQuads() const
  {
    return nBdryQuads;
  }

  GMGW_int
  getNumBdryTris() const
  {
    return nBdryTris;
  }

  GMGW_int
  getNumBdryVerts() const
  {
    return nBdryVerts;
  }

  GMGW_int
  getNumCells() const
  {
    return nCells;
  }

  GMGW_int
  getNumVerts() const
  {
    return nVerts;
  }

  GMGW_int
  getNumTris() const
  {
    assert(nBdryTris % 2 == 0);
    return (4 * nTets + 4 * nPyrs + 2 * nPrisms + nBdryTris) / 2;
  }

  GMGW_int
  getNumQuads() const
  {
    assert((nBdryQuads + nPyrs + 3 * nPrisms) % 2 == 0);
    return (nBdryQuads + nPyrs + 3 * nPrisms + 6 * nHexes) / 2;
  }

  GMGW_int
  getNumFaces() const
  {
    return getNumTris() + getNumQuads();
  }

  virtual void
  seekStartOfConnectivity() const;
  virtual void
  getNextCellConnectivity(GMGW_int& nConn, GMGW_int connect[]) const = 0;
  virtual void // needs to be virtual to enable reading from streams instead of FILEs (jww)
  seekStartOfCoords() const; 
  virtual void
  getNextVertexCoords(double& x, double& y, double &z) const = 0;
protected:
  void
  identifyBdryVerts();
};

#endif /* GMGW_FILEWRAPPER_HXX_ */
