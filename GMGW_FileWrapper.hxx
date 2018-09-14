/*
 * GMGW_FileWrapper.hxx
 *
 *  Created on: Jan 27, 2017
 *      Author: cfog
 */

#ifndef GMGW_FILEWRAPPER_HXX_
#define GMGW_FILEWRAPPER_HXX_

#include "config.h"

#include <cstdio>
#include <assert.h>

class FileWrapper {
protected:
    FILE* input;
    long int coordsStart, connectStart;
    unsigned long nVerts, nCells, nBdryTris, nBdryQuads,
        nTets, nPyrs, nPrisms, nHexes, nBdryVerts;
    char* m_cellTypes;
    bool *m_isBdryFace, *m_isBdryVert;
public:
    FileWrapper();
    virtual ~FileWrapper();

    static FileWrapper* factory(const char baseName[],
                                const char type[],
                                const char ugridInfix[]);

    virtual void scanFile() = 0;

    char getCellType(unsigned int i) const
    {
        assert(i < nCells);
        return m_cellTypes[i];
    }

    bool isBdryFace(unsigned int i) const
    {
        assert(i < nCells);
        return m_isBdryFace[i];
    }

    bool isBdryVert(unsigned int i) const
    {
        assert(i < nVerts);
        return m_isBdryVert[i];
    }

    void setBdryVert(unsigned int i) const
    {
        assert(i < nVerts);
        m_isBdryVert[i] = true;
    }

    void clearBdryVert(unsigned int i) const
    {
        assert(i < nVerts);
        m_isBdryVert[i] = false;
    }

    unsigned long getNumBdryQuads() const
    {
        return nBdryQuads;
    }

    unsigned long getNumBdryTris() const
    {
        return nBdryTris;
    }

    unsigned long getNumBdryVerts() const
    {
        return nBdryVerts;
    }

    unsigned long getNumCells() const
    {
        return nCells;
    }

    unsigned long getNumVerts() const
    {
        return nVerts;
    }

    unsigned long getNumTris() const
    {
        assert(nBdryTris % 2 == 0);
        return (4*nTets + 4*nPyrs + 2*nPrisms + nBdryTris) / 2;
    }

    unsigned long getNumQuads() const
    {
        assert((nBdryQuads + nPyrs + 3*nPrisms) % 2 == 0);
        return (nBdryQuads + nPyrs + 3*nPrisms + 6*nHexes) / 2;
    }

    unsigned long getNumFaces() const
    {
        return getNumTris() + getNumQuads();
    }

    virtual void seekStartOfConnectivity() const;
    virtual void getNextCellConnectivity(int& nConn, unsigned int connect[]) const = 0;
    void seekStartOfCoords() const;
    virtual void getNextVertexCoords(double& x, double& y, double &z) const = 0;
protected:
    void identifyBdryVerts();
};




#endif /* GMGW_FILEWRAPPER_HXX_ */
