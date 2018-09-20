/*
 * faceSorting.cxx
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

#include <algorithm>
#include <set>
#include <unordered_set>
#include <iostream>

using std::cout;
using std::endl;

#include <sys/sysinfo.h>

#include "GMGW_sort.hxx"
#include "GMGW_VTKFileWrapper.hxx"

//typedef std::unordered_set<triSort, triHash> triSet;
//typedef std::unordered_set<quadSort, quadHash> quadSet;

typedef std::set<triSort> triSet;
typedef std::set<quadSort> quadSet;

static void
buildFaceListSort(FileWrapper* reader, GMGW_int faceToCell[][2],
		  GMGW_int& badTris, GMGW_int& badQuads)
{
  // The array for the face-to-cell data is already allocated, but
  // the stuff that needs sorting isn't yet.
  const GMGW_int nTris = reader->getNumTris();
  const GMGW_int nQuads = reader->getNumQuads();
  const GMGW_int nCells = reader->getNumCells();
  triSort *triData = new triSort[2 * nTris];
  quadSort *quadData = new quadSort[2 * nQuads];

  GMGW_int tri = 0, quad = 0;
  reader->seekStartOfConnectivity();

  GMGW_int connect[8];
  for (GMGW_int iC = 0; iC < nCells; iC++) {
    GMGW_int nConn;
    reader->getNextCellConnectivity(nConn, connect);

    // Now it's time to do something with that data.
    if (reader->isBdryFace(iC)) {
      switch (nConn)
	{
	case 3:
	  // Definitely a tri.
	  triData[tri].set(iC, connect[0], connect[1], connect[2]);
	  tri++;
	  break;
	case 4:
	  // Definitely a quad.
	  quadData[quad].set(iC, connect[0], connect[1], connect[2],
			     connect[3]);
	  quad++;
	  break;
	default:
	  assert(0);
	  break;
	}
    }
    else {
      // On these, orientation doesn't matter, nor does order of verts
      // within the face, given the way the sort is being done.
      switch (nConn)
	{
	case 4:
	  // A tet
	  triData[tri++].set(iC, connect[0], connect[1], connect[2]);
	  triData[tri++].set(iC, connect[1], connect[2], connect[3]);
	  triData[tri++].set(iC, connect[2], connect[3], connect[0]);
	  triData[tri++].set(iC, connect[3], connect[0], connect[1]);
	  break;
	case 5:
	  // A pyramid
	  quadData[quad++].set(iC, connect[0], connect[1], connect[2],
			       connect[3]);
	  triData[tri++].set(iC, connect[0], connect[1], connect[4]);
	  triData[tri++].set(iC, connect[1], connect[2], connect[4]);
	  triData[tri++].set(iC, connect[2], connect[3], connect[4]);
	  triData[tri++].set(iC, connect[3], connect[0], connect[4]);
	  break;
	case 6:
	  // A prism
	  quadData[quad++].set(iC, connect[0], connect[2], connect[5],
			       connect[3]);
	  quadData[quad++].set(iC, connect[2], connect[1], connect[4],
			       connect[5]);
	  quadData[quad++].set(iC, connect[1], connect[0], connect[3],
			       connect[4]);
	  triData[tri++].set(iC, connect[0], connect[1], connect[2]);
	  triData[tri++].set(iC, connect[3], connect[4], connect[5]);
	  break;
	case 8:
	  // A hex
	  quadData[quad++].set(iC, connect[0], connect[1], connect[2],
			       connect[3]);
	  quadData[quad++].set(iC, connect[0], connect[1], connect[5],
			       connect[4]);
	  quadData[quad++].set(iC, connect[1], connect[2], connect[6],
			       connect[5]);
	  quadData[quad++].set(iC, connect[2], connect[3], connect[7],
			       connect[6]);
	  quadData[quad++].set(iC, connect[3], connect[0], connect[4],
			       connect[7]);
	  quadData[quad++].set(iC, connect[4], connect[5], connect[6],
			       connect[7]);
	  break;
	default:
	  assert(0);
	  break;
	} // end switch
    } // end else
    assert(tri <= 2 * nTris);
    assert(quad <= 2 * nQuads);
    if ((iC + 1) % 5000000 == 0) {
      cout << (iC + 1) / 1000000 << "M ";
      cout.flush();
    }
  } // end loop
  cout << endl;
  assert(tri == 2 * nTris);
  assert(quad == 2 * nQuads);

  cout << "Sorting " << tri << " tri entries..." << endl;
  // Now sort them all
  std::sort(triData, triData + 2 * nTris, triCompare);

  if (nQuads != 0) {
    cout << "Sorting " << quad << " quad entries..." << endl;
    std::sort(quadData, quadData + 2 * nQuads, quadCompare);
  }

  // Identify faces that are single sided (only a cell on one side of them)
  // For the others, transcribe data to the face-to-cell connectivity table
  badTris = badQuads = 0;

  cout << "Transcribing face-to-cell data for tris" << endl;
  GMGW_int face = 0;
  for (tri = 0; tri < 2 * nTris;) {
    if (triData[tri] == triData[tri + 1]) {
      faceToCell[face][0] = triData[tri].m_cell;
      faceToCell[face][1] = triData[tri + 1].m_cell;
      tri += 2;
      face++;
    }
    else {
      badTris++;
      tri++;
    }
  }
  assert(tri == 2 * nTris);
  badTris /= 2;

  cout << "Total tris: " << nTris << "  Bad tris: " << badTris
      << "  Tris transcribed: " << face << endl;

  if (nQuads > 0) {
    badQuads = 0;
    cout << "Transcribing face-to-cell data for quads" << endl;
    for (quad = 0; quad < 2 * nQuads;) {
      if (quadData[quad] == quadData[quad + 1]) {
	faceToCell[face][0] = quadData[quad].m_cell;
	faceToCell[face][1] = quadData[quad + 1].m_cell;
	quad += 2;
	face++;
      }
      else {
	badQuads++;
	quad++;
      }
    }
    assert(quad == 2 * nQuads);
    badQuads /= 2;
    cout << "Total quads: " << nQuads << " Bad quads: " << badQuads
	<< "  Quads transcribed: " << face - nTris << endl;
  }
  if (badTris + badQuads > 0) {
    cout << "Found " << badTris << " hanging tris, " << badQuads
	<< " hanging quads" << endl;
  }

  delete[] triData;
  delete[] quadData;
}

GMGW_int
matchTri(triSet& triData, GMGW_int& tri, size_t& trisInSet,
	 size_t& maxTrisInSet, GMGW_int triToCell[][2], const GMGW_int iC,
	 const GMGW_int v0, const GMGW_int v1, const GMGW_int v2)
{
  triSort triTemp;
  triTemp.set(iC, v0, v1, v2);

  triSet::iterator iter = triData.find(triTemp);
  if (iter == triData.end()) {
    // Didn't find it; add it to the set.
    triData.insert(triTemp);
    trisInSet++;
    if (trisInSet > maxTrisInSet) {
      maxTrisInSet = trisInSet;
    }
    return +1;
  }
  else {
    // Found one!  Make an entry into faceToCell, then delete the set entry.
    triToCell[tri][0] = iC;
    triToCell[tri][1] = iter->m_cell;
    tri++;
    trisInSet--;
    triData.erase(iter);
    return -1;
  }
}

GMGW_int
matchQuad(quadSet& quadData, GMGW_int& quad, size_t& quadsInSet,
	  size_t& maxQuadsInSet, GMGW_int quadToCell[][2], const GMGW_int iC,
	  const GMGW_int v0, const GMGW_int v1, const GMGW_int v2,
	  const GMGW_int v3)
{
  quadSort quadTemp;
  quadTemp.set(iC, v0, v1, v2, v3);

  quadSet::iterator iter = quadData.find(quadTemp);
  if (iter == quadData.end()) {
    // Didn't find it; add it to the set.
    quadData.insert(quadTemp);
    quadsInSet++;
    if (quadsInSet > maxQuadsInSet) {
      maxQuadsInSet = quadsInSet;
    }
    return +1;
  }
  else {
    // Found one!  Make an entry into faceToCell, then delete the set entry.
    quadToCell[quad][0] = iC;
    quadToCell[quad][1] = iter->m_cell;
    quad++;
    quadsInSet--;
    quadData.erase(iter);
    return -1;
  }
}

static void
buildFaceListSTL(FileWrapper* reader, GMGW_int faceToCell[][2],
		 GMGW_int& badTris, GMGW_int& badQuads)
{
  triSet triData;
  quadSet quadData;

  const GMGW_int nTris = reader->getNumTris();
  const GMGW_int nQuads = reader->getNumQuads();
  const GMGW_int nCells = reader->getNumCells();

  GMGW_int (*triToCell)[2] = faceToCell;
  GMGW_int (*quadToCell)[2] = &(faceToCell[nTris]);

  GMGW_int tri = 0, quad = 0;
  size_t trisInSet = 0, quadsInSet = 0, maxTrisInSet = 0, maxQuadsInSet = 0;
  reader->seekStartOfConnectivity();

  GMGW_int connect[8];
  for (GMGW_int iC = 0; iC < nCells; iC++) {
    GMGW_int nConn;
    reader->getNextCellConnectivity(nConn, connect);

    // Now it's time to do something with that data.
    if (reader->isBdryFace(iC)) {
      switch (nConn)
	{
	case 3:
	  // Definitely a tri.
	  matchTri(triData, tri, trisInSet, maxTrisInSet, triToCell, iC,
		   connect[0], connect[1], connect[2]);
	  break;
	case 4:
	  // Definitely a quad.
	  matchQuad(quadData, quad, quadsInSet, maxQuadsInSet, quadToCell, iC,
		    connect[0], connect[1], connect[2], connect[3]);
	  break;
	default:
	  assert(0);
	  break;
	}
    }
    else {
      // On these, orientation doesn't matter, nor does order of verts
      // within the face, given the way the sort is being done.
      switch (nConn)
	{
	case 4:
	  // A tet
	  matchTri(triData, tri, trisInSet, maxTrisInSet, triToCell, iC,
		   connect[0], connect[1], connect[2]);
	  matchTri(triData, tri, trisInSet, maxTrisInSet, triToCell, iC,
		   connect[1], connect[2], connect[3]);
	  matchTri(triData, tri, trisInSet, maxTrisInSet, triToCell, iC,
		   connect[2], connect[3], connect[0]);
	  matchTri(triData, tri, trisInSet, maxTrisInSet, triToCell, iC,
		   connect[3], connect[0], connect[1]);
	  break;
	case 5:
	  // A pyramid
	  matchQuad(quadData, quad, quadsInSet, maxQuadsInSet, quadToCell, iC,
		    connect[0], connect[1], connect[2], connect[3]);
	  matchTri(triData, tri, trisInSet, maxTrisInSet, triToCell, iC,
		   connect[0], connect[1], connect[4]);
	  matchTri(triData, tri, trisInSet, maxTrisInSet, triToCell, iC,
		   connect[1], connect[2], connect[4]);
	  matchTri(triData, tri, trisInSet, maxTrisInSet, triToCell, iC,
		   connect[2], connect[3], connect[4]);
	  matchTri(triData, tri, trisInSet, maxTrisInSet, triToCell, iC,
		   connect[3], connect[0], connect[4]);
	  break;
	case 6:
	  // A prism
	  matchQuad(quadData, quad, quadsInSet, maxQuadsInSet, quadToCell, iC,
		    connect[0], connect[2], connect[5], connect[3]);
	  matchQuad(quadData, quad, quadsInSet, maxQuadsInSet, quadToCell, iC,
		    connect[2], connect[1], connect[4], connect[5]);
	  matchQuad(quadData, quad, quadsInSet, maxQuadsInSet, quadToCell, iC,
		    connect[1], connect[0], connect[3], connect[4]);
	  matchTri(triData, tri, trisInSet, maxTrisInSet, triToCell, iC,
		   connect[0], connect[1], connect[2]);
	  matchTri(triData, tri, trisInSet, maxTrisInSet, triToCell, iC,
		   connect[3], connect[4], connect[5]);
	  break;
	case 8:
	  // A hex
	  matchQuad(quadData, quad, quadsInSet, maxQuadsInSet, quadToCell, iC,
		    connect[0], connect[1], connect[2], connect[3]);
	  matchQuad(quadData, quad, quadsInSet, maxQuadsInSet, quadToCell, iC,
		    connect[0], connect[1], connect[5], connect[4]);
	  matchQuad(quadData, quad, quadsInSet, maxQuadsInSet, quadToCell, iC,
		    connect[1], connect[2], connect[6], connect[5]);
	  matchQuad(quadData, quad, quadsInSet, maxQuadsInSet, quadToCell, iC,
		    connect[2], connect[3], connect[7], connect[6]);
	  matchQuad(quadData, quad, quadsInSet, maxQuadsInSet, quadToCell, iC,
		    connect[3], connect[0], connect[4], connect[7]);
	  matchQuad(quadData, quad, quadsInSet, maxQuadsInSet, quadToCell, iC,
		    connect[4], connect[5], connect[6], connect[7]);
	  break;
	default:
	  assert(0);
	  break;
	} // end switch
    } // end else
    if ((iC + 1) % 1000000 == 0) {
      cout << (iC + 1) / 1000000 << "M cells.  Tri set size: " << trisInSet
	  << " (max " << maxTrisInSet << ")  Quad set size: " << quadsInSet
	  << " (max " << maxQuadsInSet << ")" << endl;
    }
  } // end loop

  // Identify faces that are single sided (only a cell on one side of them)
  // For the others, transcribe data to the face-to-cell connectivity table
  badTris = nTris - tri;
  badQuads = nQuads - quad;

  cout << "Total tris:  " << nTris << "  Bad tris:  " << badTris
      << "  Tris transcribed:  " << tri << endl;
  cout << "Total quads: " << nQuads << "  Bad quads: " << badQuads
      << "  Quads transcribed: " << quad << endl;
}

void
buildFaceList(FileWrapper* reader, GMGW_int faceToCell[][2], GMGW_int& badTris,
	      GMGW_int& badQuads)
{
  // Need decision machinery here.
  struct sysinfo info;
  sysinfo(&info);

  // Available memory; willing to be pushy with other processes, at least a bit.
  size_t avail = info.totalram / 2;
  if (info.totalram > 0x200000000UL)
    avail = info.totalram - 0x100000000UL;

  // Needed for the full lists of face connectivity for sorting.
  const GMGW_int nTris = reader->getNumTris();
  const GMGW_int nQuads = reader->getNumQuads();
  size_t needed = sizeof(triSort) * nTris * 2 + sizeof(quadSort) * nQuads * 2;

  if (needed > avail) {
    cout << "Creating face-cell connectivity the small way" << endl;
    buildFaceListSTL(reader, faceToCell, badTris, badQuads);
  }
  else {
    cout << "Creating face-cell connectivity the fast way" << endl;
    buildFaceListSort(reader, faceToCell, badTris, badQuads);
  }
}
