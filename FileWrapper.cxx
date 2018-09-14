/*
 * FileWrapper.cxx
 *
 *  Created on: Jan 30, 2017
 *      Author: cfog
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "GMGW_FileWrapper.hxx"
#include "GMGW_UGridFileWrapper.hxx"
#include "GMGW_VTKFileWrapper.hxx"

FileWrapper::FileWrapper() :
        input(nullptr), coordsStart(-1), connectStart(-1),
        nVerts(0), nCells(0), nBdryTris(0), nBdryQuads(0), nTets(0),
        nPyrs(0), nPrisms(0), nHexes(0), nBdryVerts(0), m_cellTypes(nullptr),
        m_isBdryFace(nullptr), m_isBdryVert(nullptr)
{}

FileWrapper::~FileWrapper() {
     if (m_cellTypes)  delete [] m_cellTypes;
     if (m_isBdryFace) delete [] m_isBdryFace;
     if (m_isBdryVert) delete [] m_isBdryVert;
     if (input) fclose(input);
}

void FileWrapper::seekStartOfCoords() const
{
    fseek(input, coordsStart, SEEK_SET);
}

void FileWrapper::seekStartOfConnectivity() const
{
    fseek(input, connectStart, SEEK_SET);
}

void FileWrapper::identifyBdryVerts()
{
    m_isBdryVert = new bool[nVerts];
    for (unsigned int ii = 0; ii < nVerts; ii++)
        clearBdryVert(ii);
    seekStartOfConnectivity();
    // Now process all the connectivity info.
    for (unsigned int iC = 0; iC < nCells; iC++) {
        int nConn = -1;
        unsigned int connect[8];
        getNextCellConnectivity(nConn, connect);
        if (isBdryFace(iC)) {
            for (int ii = 0; ii < nConn; ii++) {
                setBdryVert(connect[ii]);
            }
        }
    }

    // Finally, count the number of bdry verts
    for (unsigned int ii = 0; ii < nVerts; ii++) {
        if (isBdryVert(ii)) {
            nBdryVerts++;
        }
    }
}

FileWrapper* FileWrapper::factory(const char baseName[],
                                  const char type[],
                                  const char ugridInfix[])
{
    if (strstr(type, "vtk")) {
        return new VTKFileWrapper(baseName);
    }
    else if (strstr(type, "ugrid")) {
        return new UGridFileWrapper(baseName, ugridInfix);
    }
    else {
        fprintf(stderr, "Missing or invalid file type: %s\n",
                type);
        exit(1);
    }
}

