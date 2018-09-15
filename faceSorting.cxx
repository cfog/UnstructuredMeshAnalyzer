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

#include <sys/sysinfo.h>

#include "GMGW_sort.hxx"
#include "GMGW_VTKFileWrapper.hxx"

//typedef std::unordered_set<triSort, triHash> triSet;
//typedef std::unordered_set<quadSort, quadHash> quadSet;

typedef std::set<triSort> triSet;
typedef std::set<quadSort> quadSet;

static void buildFaceListSort(FileWrapper* reader,
                       unsigned int faceToCell[][2],
                       unsigned int& badTris,
                       unsigned int& badQuads)
{
    // The array for the face-to-cell data is already allocated, but
    // the stuff that needs sorting isn't yet.
    const unsigned long nTris = reader->getNumTris();
    const unsigned long nQuads = reader->getNumQuads();
    const unsigned long nCells = reader->getNumCells();
    triSort *triData = new triSort[2*nTris];
    quadSort *quadData = new quadSort[2*nQuads];

    unsigned long tri = 0, quad = 0;
    reader->seekStartOfConnectivity();

    unsigned int connect[8];
    for (unsigned long iC = 0; iC < nCells; iC++) {
        int nConn;
        reader->getNextCellConnectivity(nConn, connect);

        // Now it's time to do something with that data.
        if (reader->isBdryFace(iC)) {
            switch (nConn) {
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
            switch (nConn) {
                case 4:
                    // A tet
                    triData[tri++].set(iC, connect[0], connect[1], connect[2]);
                    triData[tri++].set(iC, connect[1], connect[2], connect[3]);
                    triData[tri++].set(iC, connect[2], connect[3], connect[0]);
                    triData[tri++].set(iC, connect[3], connect[0], connect[1]);
                    break;
                case 5:
                    // A pyramid
                    quadData[quad++].set(iC, connect[0], connect[1], connect[2], connect[3]);
                    triData[tri++].set(iC, connect[0], connect[1], connect[4]);
                    triData[tri++].set(iC, connect[1], connect[2], connect[4]);
                    triData[tri++].set(iC, connect[2], connect[3], connect[4]);
                    triData[tri++].set(iC, connect[3], connect[0], connect[4]);
                    break;
                case 6:
                    // A prism
                    quadData[quad++].set(iC, connect[0], connect[2], connect[5], connect[3]);
                    quadData[quad++].set(iC, connect[2], connect[1], connect[4], connect[5]);
                    quadData[quad++].set(iC, connect[1], connect[0], connect[3], connect[4]);
                    triData[tri++].set(iC, connect[0], connect[1], connect[2]);
                    triData[tri++].set(iC, connect[3], connect[4], connect[5]);
                    break;
                case 8:
                    // A hex
                    quadData[quad++].set(iC, connect[0], connect[1], connect[2], connect[3]);
                    quadData[quad++].set(iC, connect[0], connect[1], connect[5], connect[4]);
                    quadData[quad++].set(iC, connect[1], connect[2], connect[6], connect[5]);
                    quadData[quad++].set(iC, connect[2], connect[3], connect[7], connect[6]);
                    quadData[quad++].set(iC, connect[3], connect[0], connect[4], connect[7]);
                    quadData[quad++].set(iC, connect[4], connect[5], connect[6], connect[7]);
                    break;
                default:
                    assert(0);
                    break;
            } // end switch
        } // end else
        assert(tri <= 2*nTris);
        assert(quad <= 2*nQuads);
        if ((iC+1)%5000000 == 0) {
            printf("%5luM", (iC+1)/1000000);
            fflush(stdout);
        }
    } // end loop
    printf("\n"); fflush(stdout);
    assert(tri == 2*nTris);
    assert(quad == 2*nQuads);

    printf("Sorting %'ld tri entries...\n", tri); fflush(stdout);
    // Now sort them all
    std::sort(triData, triData+2*nTris, triCompare);
    printf("Sorting %'ld quad entries...\n", quad); fflush(stdout);
    if (nQuads != 0) {
        std::sort(quadData, quadData+2*nQuads, quadCompare);
    }

    // Identify faces that are single sided (only a cell on one side of them)
    // For the others, transcribe data to the face-to-cell connectivity table
    badTris = badQuads = 0;

    printf("Transcribing face-to-cell data for tris\n"); fflush(stdout);
    unsigned long face = 0;
    for (tri = 0; tri < 2*nTris; ) {
        if (triData[tri] == triData[tri+1]) {
            faceToCell[face][0] = triData[tri].m_cell;
            faceToCell[face][1] = triData[tri+1].m_cell;
            tri += 2;
            face++;
        }
        else {
//            printf("%10u %10u %10u %10u\n", triData[tri].m_cell,
//                   triData[tri].m_verts[0], triData[tri].m_verts[1],
//                   triData[tri].m_verts[2]);
            badTris++;
            tri++;
        }
    }
    assert(tri == 2*nTris);
    badTris /= 2;

    printf("Total tris: %'ld  Bad tris:  %'d  Tris transcribed: %'ld\n",
           nTris, badTris, face);

    badQuads = 0;
    printf("Transcribing face-to-cell data for quads\n");
    fflush(stdout);
    for (quad = 0; quad < 2*nQuads; ) {
        if (quadData[quad] == quadData[quad+1]) {
            faceToCell[face][0] = quadData[quad].m_cell;
            faceToCell[face][1] = quadData[quad+1].m_cell;
            quad += 2;
            face++;
        }
        else {
//            printf("%'10u %'10u %'10u %'10u %'10u\n", quadData[quad].m_cell,
//                   quadData[quad].m_verts[0], quadData[quad].m_verts[1],
//                   quadData[quad].m_verts[2], quadData[quad].m_verts[3]);
            badQuads++;
            quad++;
        }
    }
    assert(quad == 2*nQuads);
    badQuads /= 2;
    printf("Total quads: %'ld  Bad quads:  %'d  Quads transcribed: %'ld\n",
           nQuads, badQuads, face - nTris - badQuads);
//    assert(face == nTris + nQuads - badTris - badQuads);

    printf("Found %'u hanging tris, %'u hanging quads\n", badTris, badQuads);
    fflush(stdout);

    delete [] triData;
    delete [] quadData;
}

int matchTri(triSet& triData,
             unsigned long& tri,
             size_t& trisInSet, size_t& maxTrisInSet,
             unsigned int triToCell[][2],
             const unsigned int iC, const unsigned int v0,
             const unsigned int v1, const unsigned int v2)
{
    triSort triTemp;
    triTemp.set(iC, v0, v1, v2);

    triSet::iterator iter = triData.find(triTemp);
    if (iter == triData.end()) {
        // Didn't find it; add it to the set.
        triData.insert(triTemp);
        trisInSet ++;
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

int matchQuad(quadSet& quadData,
             unsigned long& quad,
             size_t& quadsInSet, size_t& maxQuadsInSet,
             unsigned int quadToCell[][2],
             const unsigned int iC, const unsigned int v0,
             const unsigned int v1, const unsigned int v2,
             const unsigned int v3)
{
    quadSort quadTemp;
    quadTemp.set(iC, v0, v1, v2, v3);

    quadSet::iterator iter = quadData.find(quadTemp);
    if (iter == quadData.end()) {
        // Didn't find it; add it to the set.
        quadData.insert(quadTemp);
        quadsInSet ++;
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

static void buildFaceListSTL(FileWrapper* reader,
                   unsigned int faceToCell[][2],
                   unsigned int& badTris,
                   unsigned int& badQuads)
{
    triSet triData;
    quadSet quadData;
//    triSort triTemp, triTwin;
//    quadSort quadTemp, quadTwin;

    const unsigned long nTris = reader->getNumTris();
    const unsigned long nQuads = reader->getNumQuads();
    const unsigned long nCells = reader->getNumCells();

//    unsigned int (*triToCell)[2] = new unsigned int[nTris][2];
//    unsigned int (*quadToCell)[2] = new unsigned int [nQuads][2];
    unsigned int (*triToCell)[2] = faceToCell;
    unsigned int (*quadToCell)[2] = &(faceToCell[nTris]);

    unsigned long tri = 0, quad = 0;
    size_t trisInSet = 0, quadsInSet = 0, maxTrisInSet = 0, maxQuadsInSet = 0;
    reader->seekStartOfConnectivity();

    unsigned int connect[8];
    for (unsigned long iC = 0; iC < nCells; iC++) {
        int nConn;
        reader->getNextCellConnectivity(nConn, connect);

        // Now it's time to do something with that data.
        if (reader->isBdryFace(iC)) {
            switch (nConn) {
                case 3:
                    // Definitely a tri.
                    matchTri(triData, tri, trisInSet, maxTrisInSet,
                             triToCell, iC,
                             connect[0], connect[1], connect[2]);
                    break;
                case 4:
                    // Definitely a quad.
                    matchQuad(quadData, quad, quadsInSet, maxQuadsInSet,
                             quadToCell, iC,
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
            switch (nConn) {
                case 4:
                    // A tet
                    matchTri(triData, tri, trisInSet, maxTrisInSet,
                             triToCell, iC,
                             connect[0], connect[1], connect[2]);
                    matchTri(triData, tri, trisInSet, maxTrisInSet,
                             triToCell, iC,
                             connect[1], connect[2], connect[3]);
                    matchTri(triData, tri, trisInSet, maxTrisInSet,
                             triToCell, iC,
                             connect[2], connect[3], connect[0]);
                    matchTri(triData, tri, trisInSet, maxTrisInSet,
                             triToCell, iC,
                             connect[3], connect[0], connect[1]);
                    break;
                case 5:
                    // A pyramid
                    matchQuad(quadData, quad, quadsInSet, maxQuadsInSet,
                             quadToCell, iC,
                             connect[0], connect[1], connect[2], connect[3]);
                    matchTri(triData, tri, trisInSet, maxTrisInSet,
                             triToCell, iC,
                             connect[0], connect[1], connect[4]);
                    matchTri(triData, tri, trisInSet, maxTrisInSet,
                             triToCell, iC,
                             connect[1], connect[2], connect[4]);
                    matchTri(triData, tri, trisInSet, maxTrisInSet,
                             triToCell, iC,
                             connect[2], connect[3], connect[4]);
                    matchTri(triData, tri, trisInSet, maxTrisInSet,
                             triToCell, iC,
                             connect[3], connect[0], connect[4]);
                    break;
                case 6:
                    // A prism
                    matchQuad(quadData, quad, quadsInSet, maxQuadsInSet,
                             quadToCell, iC,
                             connect[0], connect[2], connect[5], connect[3]);
                    matchQuad(quadData, quad, quadsInSet, maxQuadsInSet,
                             quadToCell, iC,
                             connect[2], connect[1], connect[4], connect[5]);
                    matchQuad(quadData, quad, quadsInSet, maxQuadsInSet,
                             quadToCell, iC,
                             connect[1], connect[0], connect[3], connect[4]);
                    matchTri(triData, tri, trisInSet, maxTrisInSet,
                             triToCell, iC,
                             connect[0], connect[1], connect[2]);
                    matchTri(triData, tri, trisInSet, maxTrisInSet,
                             triToCell, iC,
                             connect[3], connect[4], connect[5]);
                    break;
                case 8:
                    // A hex
                    matchQuad(quadData, quad, quadsInSet, maxQuadsInSet,
                             quadToCell, iC,
                             connect[0], connect[1], connect[2], connect[3]);
                    matchQuad(quadData, quad, quadsInSet, maxQuadsInSet,
                             quadToCell, iC,
                             connect[0], connect[1], connect[5], connect[4]);
                    matchQuad(quadData, quad, quadsInSet, maxQuadsInSet,
                             quadToCell, iC,
                             connect[1], connect[2], connect[6], connect[5]);
                    matchQuad(quadData, quad, quadsInSet, maxQuadsInSet,
                             quadToCell, iC,
                             connect[2], connect[3], connect[7], connect[6]);
                    matchQuad(quadData, quad, quadsInSet, maxQuadsInSet,
                             quadToCell, iC,
                             connect[3], connect[0], connect[4], connect[7]);
                    matchQuad(quadData, quad, quadsInSet, maxQuadsInSet,
                             quadToCell, iC,
                             connect[4], connect[5], connect[6], connect[7]);
                    break;
                default:
                    assert(0);
                    break;
            } // end switch
        } // end else
        if ((iC+1)%1000000 == 0) {
            printf("%5luM Tris: %'13lu (max %'13lu)  Quads: %'13lu (max %'13lu)\n",
                   (iC+1)/1000000, trisInSet, maxTrisInSet, quadsInSet,
                   maxQuadsInSet);
	     fflush(stdout);
//            printf("%5uM", (iC+1)/1000000); fflush(stdout);
        }
    } // end loop
    printf("\n");

    // Identify faces that are single sided (only a cell on one side of them)
    // For the others, transcribe data to the face-to-cell connectivity table
    badTris = nTris -tri;
    badQuads = nQuads - quad;

//    printf("Copying face-to-cell data for tris\n");
//    unsigned int* result = std::copy(&(triToCell[0][0]), &(triToCell[tri][0]),
//                                     &(faceToCell[0][0]));
    printf("Total tris: %'ld  Bad tris:  %'ld  Tris transcribed: %'ld\n",
           nTris, nTris - tri, tri);
    fflush(stdout);
//    printf("Copying face-to-cell data for quads\n");
//    result = std::copy(&(quadToCell[0][0]), &(quadToCell[quad][0]), result);
//
    printf("Total quads: %'ld  Bad quads:  %'ld  Quads transcribed: %'ld\n",
           nQuads, nQuads - quad, quad);

    printf("Found %'u hanging tris, %'u hanging quads\n", badTris, badQuads);
    fflush(stdout);
//    delete [] triToCell;
//    delete [] quadToCell;
}

void buildFaceList(FileWrapper* reader,
                   unsigned int faceToCell[][2],
                   unsigned int& badTris,
                   unsigned int& badQuads)
{
    // Need decision machinery here.
    struct sysinfo info;
    sysinfo(&info);

    // Available memory; willing to be pushy with other processes, at least a bit.
    unsigned long avail = info.totalram/2;
    if (info.totalram > 0x200000000UL) avail = info.totalram - 0x100000000UL;

    // Needed for the full lists of face connectivity for sorting.
    const unsigned int nTris = reader->getNumTris();
    const unsigned int nQuads = reader->getNumQuads();
    unsigned long needed = sizeof(triSort)*nTris*2 + sizeof(quadSort)*nQuads*2;

    if (needed > avail) {
        printf("Creating face-cell connectivity the small way\n");
        buildFaceListSTL(reader, faceToCell, badTris, badQuads);
    }
    else {
        printf("Creating face-cell connectivity the fast way\n");
        buildFaceListSort(reader, faceToCell, badTris, badQuads);
    }
}
