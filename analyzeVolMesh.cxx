/*
 * vol_to_surf_vtk.cxx
 *
 *  Created on: Oct 13, 2016
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
#include <cmath>
#include <string>
#include <vector>

#include <assert.h>
#include <stdlib.h>
#include <cstdio>
#include <string.h>

#include "GMGW_unstr.hxx"
#include "GMGW_geom_utils.hxx"
#include "GMGW_FileWrapper.hxx"

void
outputSizeRatios(FileWrapper* wrapper, const unsigned int validTris,
		 const unsigned int validQuads,
		 const unsigned int faceToCell[][2],
		 const char sizeRatioFileName[], double* allVolumes)
{
  // Now let's tabulate volume ratios for histograms
  double tetTetRatio[50], tetPyrRatio[50], tetPrismRatio[50];
  double pyrPyrViaTriRatio[50], pyrPrismViaTriRatio[50],
      prismPrismViaTriRatio[50];
  double pyrPyrViaQuadRatio[50], pyrPrismViaQuadRatio[50],
      prismPrismViaQuadRatio[50];
  double pyrHexRatio[50], prismHexRatio[50], hexHexRatio[50];
  for (int ii = 0; ii < 50; ii++) {
    tetTetRatio[ii] = 0;
    tetPyrRatio[ii] = 0;
    tetPrismRatio[ii] = 0;
    pyrPyrViaTriRatio[ii] = 0;
    pyrPrismViaTriRatio[ii] = 0;
    prismPrismViaTriRatio[ii] = 0;
    pyrPyrViaQuadRatio[ii] = 0;
    pyrPrismViaQuadRatio[ii] = 0;
    prismPrismViaQuadRatio[ii] = 0;
    pyrHexRatio[ii] = 0;
    prismHexRatio[ii] = 0;
    hexHexRatio[ii] = 0;
  }
  unsigned int tetTetCount = 0, tetPyrCount = 0, tetPrismCount = 0;
  unsigned int pyrPyrViaTriCount = 0, pyrPrismViaTriCount = 0,
      prismPrismViaTriCount = 0;
  unsigned int pyrPyrViaQuadCount = 0, pyrPrismViaQuadCount = 0,
      prismPrismViaQuadCount = 0;
  unsigned int pyrHexCount = 0, prismHexCount = 0, hexHexCount = 0;
  for (unsigned int ui = 0; ui < validTris; ui++) {
    unsigned int cellA = faceToCell[ui][0];
    unsigned int cellB = faceToCell[ui][1];
    double sizeA = allVolumes[cellA];
    double sizeB = allVolumes[cellB];
    if (sizeA <= 0 || sizeB <= 0)
      continue;
    double ratio = std::min(sizeA / sizeB, sizeB / sizeA);
    int bin = floor(ratio * 50);

    char typeA = wrapper->getCellType(cellA);
    char typeB = wrapper->getCellType(cellB);

    if (typeA == TET && typeB == TET) {
      tetTetRatio[bin]++;
      tetTetCount++;
    }
    else if (((typeA == TET) && (typeB == PYRAMID))
	|| ((typeA == PYRAMID) && (typeB == TET))) {
      tetPyrRatio[bin]++;
      tetPyrCount++;
    }
    else if (((typeA == TET) && (typeB == PRISM))
	|| ((typeA == PRISM) && (typeB == TET))) {
      tetPrismRatio[bin]++;
      tetPrismCount++;
    }
    if (typeA == PYRAMID && typeB == PYRAMID) {
      pyrPyrViaTriRatio[bin]++;
      pyrPyrViaTriCount++;
    }
    else if (((typeA == PRISM) && (typeB == PYRAMID))
	|| ((typeA == PYRAMID) && (typeB == PRISM))) {
      pyrPrismViaTriRatio[bin]++;
      pyrPrismViaTriCount++;
    }
    else if ((typeB == PRISM) && (typeA == PRISM)) {
      prismPrismViaTriRatio[bin]++;
      prismPrismViaTriCount++;
    }
  }
  // Now for the quads
  for (unsigned int ui = validTris; ui < validTris + validQuads; ui++) {
    unsigned int cellA = faceToCell[ui][0];
    unsigned int cellB = faceToCell[ui][1];
    double sizeA = allVolumes[cellA];
    double sizeB = allVolumes[cellB];
    if (sizeA <= 0 || sizeB <= 0)
      continue;
    double ratio = std::min(sizeA / sizeB, sizeB / sizeA);
    int bin = floor(ratio * 50);

    char typeA = wrapper->getCellType(cellA);
    char typeB = wrapper->getCellType(cellB);

    if (typeA == HEX && typeB == HEX) {
      hexHexRatio[bin]++;
      hexHexCount++;
    }
    else if (((typeA == HEX) && (typeB == PYRAMID))
	|| ((typeA == PYRAMID) && (typeB == HEX))) {
      pyrHexRatio[bin]++;
      pyrHexCount++;
    }
    else if (((typeA == HEX) && (typeB == PRISM))
	|| ((typeA == PRISM) && (typeB == HEX))) {
      prismHexRatio[bin]++;
      prismHexCount++;
    }
    if (typeA == PYRAMID && typeB == PYRAMID) {
      pyrPyrViaQuadRatio[bin]++;
      pyrPyrViaQuadCount++;
    }
    else if (((typeA == PRISM) && (typeB == PYRAMID))
	|| ((typeA == PYRAMID) && (typeB == PRISM))) {
      pyrPrismViaQuadRatio[bin]++;
      pyrPrismViaQuadCount++;
    }
    else if ((typeB == PRISM) && (typeA == PRISM)) {
      prismPrismViaQuadRatio[bin]++;
      prismPrismViaQuadCount++;
    }
  }
  // Now output all that data, normalized.
  FILE* sizeRatios = fopen(sizeRatioFileName, "w");
  fprintf(
      sizeRatios,
      "# Bin val  tet-tet tet-pyr tet-prism pyr-pyr(t) pyr-prism(t) prism-prism(t) pyr-pyr(q) pyr-prism(q) prism-prism(q) pyr-hex prism-hex hex-hex\n");
  for (int ii = 0; ii < 50; ii++) {
    fprintf(
	sizeRatios,
	"%4f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n",
	(ii + 0.5) / 50,
	tetTetCount ? tetTetRatio[ii] / tetTetCount : 0,
	tetPyrCount ? tetPyrRatio[ii] / tetPyrCount : 0,
	tetPrismCount ? tetPrismRatio[ii] / tetPrismCount : 0,
	pyrPyrViaTriCount ? pyrPyrViaTriRatio[ii] / pyrPyrViaTriCount : 0,
	pyrPrismViaTriCount ? pyrPrismViaTriRatio[ii] / pyrPrismViaTriCount : 0,
	prismPrismViaTriCount ?
	    prismPrismViaTriRatio[ii] / prismPrismViaTriCount : 0,
	pyrPyrViaQuadCount ? pyrPyrViaQuadRatio[ii] / pyrPyrViaQuadCount : 0,
	pyrPrismViaQuadCount ?
	    pyrPrismViaQuadRatio[ii] / pyrPrismViaQuadCount : 0,
	prismPrismViaQuadCount ?
	    prismPrismViaQuadRatio[ii] / prismPrismViaQuadCount : 0,
	pyrHexCount ? pyrHexRatio[ii] / pyrHexCount : 0,
	prismHexCount ? prismHexRatio[ii] / prismHexCount : 0,
	hexHexCount ? hexHexRatio[ii] / hexHexCount : 0);
  }
  fprintf(sizeRatios,
	  "# Counts %8u %8u %8u %8u %8u %8u %8u %8u %8u %8u %8u %8u\n",
	  tetTetCount, tetPyrCount, tetPrismCount, pyrPyrViaTriCount,
	  pyrPrismViaTriCount, prismPrismViaTriCount, pyrPyrViaQuadCount,
	  pyrPrismViaQuadCount, prismPrismViaQuadCount, pyrHexCount,
	  prismHexCount, hexHexCount);
  fclose(sizeRatios);
}



void
printDataStats(const double ignoreValue, unsigned int dataSize, double* data,
	       const char prefix[])
{
  std::sort(data, data + dataSize);
  int ignoreCount = 0;
  while (data[ignoreCount] <= ignoreValue) {
    ignoreCount++;
  }
  int nActualVals = dataSize - ignoreCount;
  int minIndex = ignoreCount;
  int fivePercentileIndex = ignoreCount + (nActualVals * 5) / 100;
  int medianIndex = ignoreCount + nActualVals / 2;
  int ninetyfivePercentileIndex = ignoreCount + (nActualVals * 95) / 100;
  int maxIndex = dataSize - 1;
  printf("%s stats: \n", prefix);
  printf("   min: %12E\n", data[minIndex]);
  printf("    5%%: %12E\n", data[fivePercentileIndex]);
  printf("   50%%: %12E\n", data[medianIndex]);
  printf("   95%%: %12E\n", data[ninetyfivePercentileIndex]);
  printf("   max: %12E\n", data[maxIndex]);
  fflush (stdout);
}

void
writeBdryFaceConnectivity(FILE* output, const int nConn,
			  const unsigned int connect[], const int newIndex[])
{
  fprintf(output, "%d ", nConn);
  for (int ii = 0; ii < nConn; ii++) {
    assert(newIndex[connect[ii]] >= 0);
    fprintf(output, "%u ", newIndex[connect[ii]]);
  }
  fprintf(output, "\n");
}


// Write out a new VTK file with only the surface mesh in it.
bool
writeSurfaceMesh(FileWrapper* wrapper, const char outputFileName[],
		 const char sizeRatioFileName[], const char nmbFileName_arg[],
		 const char angleFileName[], const char distortFileName[],
		 const unsigned int faceToCell[][2],
		 const unsigned int validTris, const unsigned int validQuads)
{
  // Open the output file
  FILE* output = fopen(outputFileName, "w");
  assert(output);
  std::string nmbFilename(nmbFileName_arg);

  // Write the file header
  fprintf(output, "# vtk DataFile Version 2.0\n");
  fprintf(output, "Surface extracted from volume mesh\n");
  fprintf(output, "ASCII\n");
  fprintf(output, "DATASET UNSTRUCTURED_GRID\n");

  // Write bdry vertex data
  wrapper->seekStartOfCoords();
  unsigned int nBdryVerts = wrapper->getNumBdryVerts();
  unsigned int nVerts = wrapper->getNumVerts();
  fprintf(output, "POINTS %u float\n", nBdryVerts);
  int* newIndex = new int[nVerts];

  double (*coords)[3] = new double[nVerts][3];
  double (*bdryCoords)[3] = new double[nBdryVerts][3];
  unsigned int iBV = 0;

  printf("Reading verts... ");
  fflush (stdout);
  double minY = 0;
  for (unsigned int iV = 0; iV < nVerts; iV++) {
    // Read the coords anyway, because we have to move through the file.
    double x, y, z;
    wrapper->getNextVertexCoords(x, y, z);
    if (y < minY)
      minY = y;
    if (wrapper->isBdryVert(iV)) {
      // This one's on the bdry!
      newIndex[iV] = iBV;
      bdryCoords[iBV][0] = x;
      bdryCoords[iBV][1] = y;
      bdryCoords[iBV][2] = z;
      iBV++;
    }
    else {
      newIndex[iV] = -1; // No mistaking this for a valid index!
    }
    coords[iV][0] = x;
    coords[iV][1] = y;
    coords[iV][2] = z;

    if ((iV + 1) % 2000000 == 0) {
      printf("%dM ", (iV + 1) / 1000000);
      fflush(stdout);
    }
  }
  printf("\n");
  double sign = 1;
  if (minY < -10) {
    sign = -1; // Need to flip signs on all the volumes, too.
    printf("Swapping y and z coords...\n");
    fflush(stdout);
    // Flip coords for meshes with y up and z out the wingtip.
    for (unsigned int ii = 0; ii < nBdryVerts; ii++) {
      double tmp = bdryCoords[ii][1];
      bdryCoords[ii][1] = bdryCoords[ii][2];
      bdryCoords[ii][2] = tmp;
    }

    for (unsigned int ii = 0; ii < nVerts; ii++) {
      std::swap(coords[ii][1], coords[ii][2]);
    }
  }
  printf("Writing verts...");
  for (unsigned int ii = 0; ii < nBdryVerts; ii++) {
    double& x = bdryCoords[ii][0];
    double& y = bdryCoords[ii][1];
    double& z = bdryCoords[ii][2];
    fprintf(output, "%15.12f %15.12f %15.12f\n", x, y, z);
    if ((ii + 1) % 100000 == 0) {
      printf("%dK ", (ii + 1) / 1000);
      fflush(stdout);
    }
  }
  printf("\n");

  assert(iBV == nBdryVerts);

  // Write bdry face data (with connectivity indices updated!  Also,
  // accumulate info about nearest point off the surface.
  printf("Transcribing cell data and finding boundary spacing...");
  fflush(stdout);
  int numNeg = 0;

  // Initially, set wall spacing to something ridiculous so that min will
  // always pick the edge length the first time it tries.
  double *volSpacing = new double[nBdryVerts];
  double *skinSpacing = new double[nBdryVerts];
  for (unsigned int ii = 0; ii < nBdryVerts; ii++)
    skinSpacing[ii] = volSpacing[ii] = 1.e100;

  double totalVolume = 0;

  wrapper->seekStartOfConnectivity();

  unsigned int nBdryTris = wrapper->getNumBdryTris();
  unsigned int nBdryQuads = wrapper->getNumBdryQuads();
  unsigned int nCells = wrapper->getNumCells();

  // A spot to store cell volumes for later volume ratio calculations
  double *allVolumes = new double[nCells];
  // Obviously bdry entities have zero volume.

  int triFaceAngles[30], quadFaceAngles[30], quadDistortion[30];
  int dihedralsQuadQuad[30], dihedralsQuadTri[30], dihedralsTriTri[30];
  for (int ii = 0; ii < 30; ii++) {
    triFaceAngles[ii] = quadFaceAngles[ii] = quadDistortion[ii] =
	dihedralsQuadQuad[ii] = dihedralsQuadTri[ii] = dihedralsTriTri[ii] = 0;
  }

  std::vector<unsigned int> badElementVertList;
  std::vector<unsigned int> badElementConnect;
  std::vector<unsigned int> badElementType;
  unsigned int dummyVert = 0;

  static const int VTKtype[] =
    { 0, 0, 0, 0, TET, PYRAMID, PRISM, 0, HEX };

  fprintf(output, "CELLS %u %u\n", nBdryTris + nBdryQuads,
	  4 * nBdryTris + 5 * nBdryQuads);
  unsigned int connect[8];
  for (unsigned int iC = 0; iC < nCells; iC++) {
    int nConn;
    wrapper->getNextCellConnectivity(nConn, connect);
    // Write the connectivity only for bdry tris and quads
    if (wrapper->isBdryFace(iC)) {
      allVolumes[iC] = 0;
      writeBdryFaceConnectivity(output, nConn, connect, newIndex);
      findOnWallSpacing(coords, nConn, connect, newIndex, skinSpacing);

      // For all faces, find their face angles and bin the results.
      // For quad boundary faces, find their non-planarity; bin the result.
      analyzeBdryFace(coords, nConn, connect, triFaceAngles, quadFaceAngles,
		      quadDistortion);
    }
    // Otherwise, check for the closest interior point for any bdry verts.
    else {
      double thisVol = cellVolume(coords, nConn, connect, sign);
      allVolumes[iC] = thisVol;
      if (thisVol < 0) {
	numNeg++;
	badElementConnect.push_back(nConn);
	badElementType.push_back(VTKtype[nConn]);
	for (int ii = 0; ii < nConn; ii++) {
	  badElementVertList.push_back(connect[ii]);
	  badElementConnect.push_back(dummyVert++);
	}
//                printf("Negative volume!  Cell %u, type %d, vol %.10G\n",
//                       iC, nConn, thisVol);
      }
//            printf("Cell %u, type %d, vol %.10G\n", iC, nConn, thisVol);
      totalVolume += thisVol;

      findOffWallSpacing(wrapper, coords, nConn, connect, newIndex, volSpacing);

      // For each cell, this call will do the following:

      // For all faces, find their faces angle and bin the results.
      // For quad boundary faces, find their non-planarity; bin the result.
      // For all edges, find the dihedral angle and bin the results.

      analyzeCellQuality(coords, nConn, connect, triFaceAngles, quadFaceAngles,
			 quadDistortion, dihedralsQuadQuad, dihedralsQuadTri,
			 dihedralsTriTri);
    }
    if ((iC + 1) % 5000000 == 0) {
      printf("%dM ", (iC + 1) / 1000000);
      fflush(stdout);
    }
  }
  printf("\n");
  printf("Total volume: %'.12G.  Number w/ negative volume: %d\n", totalVolume,
	 numNeg);
  delete[] newIndex;

  FILE* angleFile = fopen(angleFileName, "w");
  outputAngleHistograms(angleFile, quadFaceAngles, triFaceAngles,
			dihedralsQuadQuad, dihedralsQuadTri, dihedralsTriTri);
  fclose(angleFile);
  FILE* distortFile = fopen(distortFileName, "w");
  outputDistortionHistogram(distortFile, quadDistortion);
  fclose(distortFile);

  // if (numNeg > 0) {
  //     // Write another VTK file, with all the elements that measured with
  //     // negative volumes.
  //     FILE* badElements = fopen("bad-elements.vtk", "w");
  //     fprintf(badElements, "# vtk DataFile Version 2.0\n");
  //     fprintf(badElements, "Surface extracted from volume mesh\n");
  //     fprintf(badElements, "ASCII\n");
  //     fprintf(badElements, "DATASET UNSTRUCTURED_GRID\n");

  //     fprintf(badElements, "POINTS %lu float\n", badElementVertList.size());
  //     for (unsigned int ii = 0; ii < badElementVertList.size(); ii++) {
  //         unsigned int vert = badElementVertList[ii];
  //         fprintf(badElements, "%.12G %.12G %.12G\n",
  //                 coords[vert][0], coords[vert][1], coords[vert][2]);
  //     }
  //     fprintf(badElements, "CELLS %u %lu\n", numNeg, badElementConnect.size());
  //     for (unsigned int ii = 0; ii < badElementConnect.size(); ii++) {
  //         fprintf(badElements, "%u\n", badElementConnect[ii]);
  //     }
  //     fprintf(badElements, "CELL_TYPES %u\n", numNeg);
  //     for (unsigned int ii = 0; ii < badElementType.size(); ii++) {
  //         fprintf(badElements, "%u\n", badElementType[ii]);
  //     }
  //     fclose(badElements);
  // }
  delete[] coords;

  // Now let's tabulate volume ratios for histograms; then we can get rid of
  // some data.
  outputSizeRatios(wrapper, validTris, validQuads, faceToCell,
		   sizeRatioFileName, allVolumes);
  delete[] allVolumes;
  delete[] faceToCell;

  printf("Writing cell types\n");
  fflush(stdout);
  // Write cell types
  fprintf(output, "CELL_TYPES %u\n", nBdryTris + nBdryQuads);

  for (unsigned int iC = 0; iC < nCells; iC++) {
    if (wrapper->isBdryFace(iC)) {
      char type = wrapper->getCellType(iC);
      fprintf(output, "%d\n", type);
    }
  }

  printf("Writing off-wall spacing\n");
  fflush(stdout);
  // Write off-wall spacing as point data.
  fprintf(output, "POINT_DATA %u\n", nBdryVerts);
  fprintf(output, "SCALARS WallSpacing float\n");
  fprintf(output, "LOOKUP_TABLE default\n");
  for (unsigned int iV = 0; iV < nBdryVerts; iV++) {
    if (volSpacing[iV] > 1.E99)
      volSpacing[iV] = -1;
    fprintf(output, "%.6G\n", volSpacing[iV]);
  }
  printDataStats(-1, nBdryVerts, volSpacing, "Volume spacing");
  delete[] volSpacing;

  printf("Writing on-wall spacing\n");
  fflush(stdout);
  // Write on-wall spacing as point data.
  fprintf(output, "SCALARS SkinSpacing float\n");
  fprintf(output, "LOOKUP_TABLE default\n");
  for (unsigned int iV = 0; iV < nBdryVerts; iV++) {
    if (skinSpacing[iV] > 1.E99)
      skinSpacing[iV] = 0;
    fprintf(output, "%.6G\n", skinSpacing[iV]);
  }
  delete[] skinSpacing;

  if (nmbFilename != "NONE") {
    double* bdryDist = new double[nBdryVerts];
    int* bdrySurf = new int[nBdryVerts];
    for (unsigned int ii = 0; ii < nBdryVerts; ii++) {
      bdrySurf[ii] = -1;
    }
    int retVal = projectionChecks(nBdryVerts, bdryCoords, bdryDist, bdrySurf,
				  nmbFilename);
    if (retVal == 0) {
      // The projection checks will fail if there's no GEODE kernel out
      // there to check against.
      printf("Writing wall projection distances\n");
      fflush(stdout);
      // Write projection distance as point data.
      fprintf(output, "SCALARS ProjDistance float\n");
      fprintf(output, "LOOKUP_TABLE default\n");
      for (unsigned int iV = 0; iV < nBdryVerts; iV++) {
	fprintf(output, "%.6G\n", bdryDist[iV]);
      }
      // Extract some stats (min, 5%-ile, median, 95%-ile, max) from
      // projection distance before nuking the data.
      printDataStats(-1, nBdryVerts, bdryDist, "Boundary projection");
      delete[] bdryDist;

      // Write surface projected to as point data.
      fprintf(output, "SCALARS ProjSurf integer\n");
      fprintf(output, "LOOKUP_TABLE default\n");
      for (unsigned int iV = 0; iV < nBdryVerts; iV++) {
	fprintf(output, "%d\n", bdrySurf[iV]);
      }

      delete[] bdrySurf;
    }
  }
  delete[] bdryCoords;
  fclose(output);

  return true;
}

void
usage()
{
  fprintf(stderr,
	  "analyzeVolMesh ft base_filename [-nmb NMB_geom_filename] [infix]\n");
  fprintf(stderr, "  where ft (filetype) = vtk or ugrid\n");
  fprintf(stderr, "  base_filename has no extension (.vtk, .infix.ugrid)\n");
  fprintf(
      stderr,
      "  If linked using Pointwise's GEODE kernel, the NMB geometry file name must be given.\n");
  fprintf(stderr,
	  "  infix = b8, lb8, etc, is a valid ugrid file type specifier\n");
}

int
main(int argc, char * const argv[])
{
  // Two arguments are required: the file type and base file name (no
  // extension).  For ugrid files, an infix that specifies the particular
  // ugrid flavor is also requires.  Plus the program name is 3 (or 4).
  if (argc < 3) {
    usage();
    exit(1);
  }
  if (sizeof(long) == 4) {
    fprintf(
	stderr,
	"Long int is only four bytes!  Files larger than 4 GB may cause problems!\n");
  }
  char fileNameOut[1024], fileNameSizes[1024], fileNameNMB[1024];
  char fileNameAngles[1024], fileNameDistort[1024];
  snprintf(fileNameOut, 1000, "%s-surf.vtk", argv[2]);
  snprintf(fileNameSizes, 1000, "%s-size.dat", argv[2]);
  snprintf(fileNameAngles, 1000, "%s-angles.dat", argv[2]);
  snprintf(fileNameDistort, 1000, "%s-distort.dat", argv[2]);
  char infixDflt[] = "b8";
  int infixArgNum = 3;
  if (strncmp(argv[3], "-nmb", 4) == 0) {
    snprintf(fileNameNMB, 1000, "%s", argv[4]);
    infixArgNum = 5;
  }
  else {
    snprintf(fileNameNMB, 1000, "NONE");
  }
  char* infixPtr = (argc == infixArgNum) ? infixDflt : argv[infixArgNum];

  FileWrapper* reader = FileWrapper::factory(argv[2], argv[1], infixPtr);

  reader->scanFile();

  unsigned int nTris = reader->getNumTris();
  unsigned int nQuads = reader->getNumQuads();
  unsigned int nFaces = reader->getNumFaces();

  unsigned int (*faceToCell)[2] = new unsigned int[nFaces][2];
  unsigned int nBadTris = 0, nBadQuads = 0;
  buildFaceList(reader, faceToCell, nBadTris, nBadQuads);
  unsigned int nValidTris = nTris - nBadTris;
  unsigned int nValidQuads = nQuads - nBadQuads;

  writeSurfaceMesh(reader, fileNameOut, fileNameSizes, fileNameNMB,
		   fileNameAngles, fileNameDistort, faceToCell,
		   nValidTris, nValidQuads);

  printf("Deleting the last of the data\n");
  return 0;
}
