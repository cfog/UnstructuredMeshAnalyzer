/*
 * GMGW_unstr.hxx
 *
 *  Created on: Jan 27, 2017
 *      Author: cfog
 */

#ifndef GMGW_UNSTR_HXX_
#define GMGW_UNSTR_HXX_

#include <string>
#include <GMGW_FileWrapper.hxx>

#define MIN(a,b) ((a) < (b) ? (a) : (b))

// VTK cell type definitions

#define BDRY_TRI 5
#define BDRY_QUAD 9
#define TET 10
#define PYRAMID 14
#define PRISM 13
#define HEX 12

int projectionChecks(const int nBdryVerts,
                     const double bdryCoords[][3], double bdryDist[],
                     int bdrySurf[], std::string nmbFileName);

void buildFaceList(FileWrapper* reader,
                   unsigned int faceToCell[][2],
                   unsigned int& badTris,
                   unsigned int& badQuads);

#endif /* GMGW_UNSTR_HXX_ */
