/*
 * geom_utils.cxx
 *
 *  Created on: Oct 13, 2017
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

#include <cmath>
#include <cassert>
#include <fstream>
#include <iomanip>

#include "GMGW_geom_utils.hxx"
#include "GMGW_unstr.hxx"

static const double angleBinBdrys[] =
  { 0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72, 78, 84, 90, 96, 102, 108,
      114, 120, 126, 132, 138, 144, 150, 156, 162, 168, 174, 180 };
static const GMGW_int nAngleBins = 30;

static const double distortBinBdrys[] =
  { 0, 1. / 1048576, 1. / 524288, 1. / 262144, 1. / 131072, 1. / 65536, 1.
      / 32768, 1. / 16384, 1. / 8192, 1. / 4096, 1. / 2048, 1. / 1024, 1. / 512,
      1. / 256, 1. / 128, 1. / 64, 1. / 32, 1. / 16, 1. / 8, 1. / 4, 1. / 2, 1 };
static const GMGW_int nDistortBins = 21;

void
outputAngleHistograms(FILE* angleFile, const GMGW_int quadFaceAngles[],
		      const GMGW_int triFaceAngles[],
		      const GMGW_int dihedralsQuadQuad[],
		      const GMGW_int dihedralsQuadTri[],
		      const GMGW_int dihedralsTriTri[])
{
  double totalQFA(0), totalTFA(0), totalDQQ(0), totalDQT(0), totalDTT(0);
  for (GMGW_int ii = 0; ii < nAngleBins; ii++) {
    totalQFA += quadFaceAngles[ii];
    totalTFA += triFaceAngles[ii];
    totalDQQ += dihedralsQuadQuad[ii];
    totalDQT += dihedralsQuadTri[ii];
    totalDTT += dihedralsTriTri[ii];
  }
  fprintf(angleFile, "%12s %12s %12s %12s %12s %12s\n", "#Bin-mid", "quad face",
	  "tri face", "q-q dihed", "q-t dihed", "t-t dihed");
  for (GMGW_int ii = 0; ii < nAngleBins; ii++) {
    fprintf(angleFile, "%10.0f %11.3f%% %11.3f%% %11.3f%% %11.3f%% %11.3f%%\n",
	    0.5 * (angleBinBdrys[ii] + angleBinBdrys[ii + 1]),
	    quadFaceAngles[ii] / totalQFA * 100,
	    triFaceAngles[ii] / totalTFA * 100,
	    dihedralsQuadQuad[ii] / totalDQQ * 100,
	    dihedralsQuadTri[ii] / totalDQT * 100,
	    dihedralsTriTri[ii] / totalDTT * 100);
  }
}

void
outputDistortionHistogram(const char* distortFileName,
			  const GMGW_int quadDistortion[])
{
  std::fstream distortFile(distortFileName);

  double total = 0;
  for (GMGW_int ii = 0; ii < nDistortBins; ii++) {
    total += quadDistortion[ii];
  }
  distortFile << "# Bin-val     distort" << std::endl;
  for (GMGW_int ii = 0; ii < nDistortBins; ii++) {
    distortFile << "<" << ii - nDistortBins + 1 << " " << std::setprecision(3)
	<< quadDistortion[ii] / total * 100 << "%" << std::endl;
  }
  distortFile.close();
}

double
distance(const double a[3], const double b[3])
{
  return sqrt(
      (b[0] - a[0]) * (b[0] - a[0]) + (b[1] - a[1]) * (b[1] - a[1])
	  + (b[2] - a[2]) * (b[2] - a[2]));
}

void
findClosestPoint(const GMGW_int nConn, const GMGW_int connect[], bool onBdry[],
		 const double coords[][3], double minDist[])
{
  // First, filter out points on the symmetry plane and on the far field,
  // and mark those as non-bdry.  This will clean up the visualization of
  // these quantities a lot.
  const double radius = 4000; // Outside this is considered farfield.
  const double symmetryCapture = 1.e-4; // y < this is considered on the symmetry plane.

  for (GMGW_int ii = 0; ii < nConn; ii++) {
    const double* myCoords = coords[connect[ii]];
    if (myCoords[1] < symmetryCapture) {
      onBdry[ii] = false;
    }
    else if ((myCoords[0] * myCoords[0] + myCoords[1] * myCoords[1]
	+ myCoords[2] * myCoords[2]) > radius * radius) {
      onBdry[ii] = false;
    }
  }

  bool allOnBdry = true;
  for (GMGW_int ii = 0; ii < nConn; ii++) {
    if (!onBdry[ii]) {
      allOnBdry = false;
      break;
    }
  }

  if (allOnBdry) {
    // Compute all distances between points; each point is assigned the
    // smallest dist.  This will (at least mostly) patch up intersections
    // of the geometry with a symmetry plane.
    for (GMGW_int ii = 0; ii < nConn - 1; ii++) {
      for (GMGW_int jj = ii + 1; jj < nConn; jj++) {
	double dist = distance(coords[connect[ii]], coords[connect[jj]]);
	minDist[ii] = MIN(minDist[ii], dist);
	minDist[jj] = MIN(minDist[jj], dist);
      }
    }
  }
  else {
    // In this case, only compute bdry-interior pair distances
    for (GMGW_int ii = 0; ii < nConn; ii++) {
      if (!onBdry[ii])
	continue;
      for (GMGW_int jj = 0; jj < nConn; jj++) {
	if (onBdry[jj] || (jj == ii))
	  continue;
	double dist = distance(coords[connect[ii]], coords[connect[jj]]);
	minDist[ii] = MIN(minDist[ii], dist);
      }
    }
  }
}

static void
addToBins(const double data, const GMGW_int nBins, const double binBdrys[],
	  GMGW_int binCounts[])
{
  assert(data >= binBdrys[0] && data <= binBdrys[nBins]);
  for (GMGW_int ii = 0; ii < nBins; ii++) {
    if (data < binBdrys[ii + 1]) {
      binCounts[ii]++;
      break;
    }
  }
}

static double
angleBetweenVecs(const double vecA[3], const double vecB[3],
		 const GMGW_int sign = 1)
{
  double cosine = sign * DOT(vecA, vecB);
  double cross[] = CROSS(vecA, vecB);
  double sine = MAG(cross);
  double angle = atan2(sine, cosine);
  assert(angle >= 0 && angle <= M_PI);
  return angle * 180 / M_PI;
}

static double
triArea(const double coords[][3], const GMGW_int v0, const GMGW_int v1,
	const GMGW_int v2)
{
  double e0 = distance(coords[v0], coords[v1]);
  double e1 = distance(coords[v1], coords[v2]);
  double e2 = distance(coords[v2], coords[v0]);

  // This specific ordering is required for numerical stability.
  if (e0 < e1) {
    std::swap(e0, e1);
  }
  // Now e0 > e1
  if (e1 < e2) {
    std::swap(e1, e2);
    // Now e1 > e2, but we don't know about e0 vs the new e1
    if (e0 < e1) {
      std::swap(e0, e1);
    }
  }
  assert(e0 >= e1);
  assert(e1 >= e2);
  if (e2 - (e0 -e1) < 0) {
	  // These three edge lengths can't be a real triangle, because e0
	  // is larger than the sum of e1 and e2.  Treat this triangle as
	  // exactly linear, and return zero area.
	  return 0;
  }
  // The parentheses enforce a numerically stable order of operations.
  double Area = 0.25
      * sqrt(
	  (e0 + (e1 + e2)) * (e2 - (e0 - e1)) * (e2 + (e0 - e1))
	      * (e0 + (e1 - e2)));
  return Area;
}

void
findOnWallSpacing(const double coords[][3], const GMGW_int nConn,
		  const GMGW_int connect[], const GMGW_int newIndex[],
		  double skinSpacing[])
{
  bool onBdry[] =
    { true, true, true, true };
  double minDist[] =
    { 1e100, 1e100, 1e100, 1e100 };
  findClosestPoint(nConn, connect, onBdry, coords, minDist);
  for (GMGW_int ii = 0; ii < nConn; ii++) {
    GMGW_int oldIdx = connect[ii];
    GMGW_int newIdx = newIndex[oldIdx];
    if ((minDist[ii] < 1e99) && (newIdx >= 0)) {
      skinSpacing[newIdx] =
	  (minDist[ii] < skinSpacing[newIdx]) ?
	      minDist[ii] : skinSpacing[newIdx];
    }
  }

}

void
findOffWallSpacing(const FileWrapper* wrapper, const double coords[][3],
		   const GMGW_int nConn, const GMGW_int connect[],
		   const GMGW_int newIndex[], double volSpacing[])
{
  // Which verts are on the bdry?
  bool onBdry[8], anyOnBdry = false;
  for (GMGW_int ii = 0; ii < nConn; ii++) {
    // Is this a bdry vert?
    if (wrapper->isBdryVert(connect[ii])) {
      onBdry[ii] = true;
      anyOnBdry = true;
    }
    else {
      onBdry[ii] = false;
    }
  }
  if (anyOnBdry) {
    double minDist[] =
      { 1e100, 1e100, 1e100, 1e100, 1e100, 1e100, 1e100, 1e100 };
    findClosestPoint(nConn, connect, onBdry, coords, minDist);
    for (GMGW_int ii = 0; ii < nConn; ii++) {
      GMGW_int oldIdx = connect[ii];
      GMGW_int newIdx = newIndex[oldIdx];
      if ((minDist[ii] < 1e99) && (newIdx >= 0)) {
	volSpacing[newIdx] =
	    (minDist[ii] < volSpacing[newIdx]) ?
		minDist[ii] : volSpacing[newIdx];
      }
    }
  }
}

static void
analyzeTri(const double coords[][3], const GMGW_int v0, const GMGW_int v1,
	   const GMGW_int v2, double angles[], double normal[])
{
  double vec01[] = SUB(coords[v0], coords[v1]);
  double vec12[] = SUB(coords[v1], coords[v2]);
  double vec20[] = SUB(coords[v2],
      coords[v0]);
  angles[0] = angleBetweenVecs(vec01, vec12, -1);
  angles[1] = angleBetweenVecs(vec12, vec20, -1);
  angles[2] = angleBetweenVecs(vec20, vec01, -1);
  double junk[] = CROSS(
      vec01, vec12);
  normal[0] = -junk[0];
  normal[1] = -junk[1];
  normal[2] = -junk[2];
}

static void
analyzeQuad(const double coords[][3], const GMGW_int v0, const GMGW_int v1,
	    const GMGW_int v2, const GMGW_int v3, double angles[],
	    double normal[])
{
  double vec01[] = SUB(coords[v0], coords[v1]);
  double vec12[] = SUB(coords[v1], coords[v2]);
  double vec23[] = SUB(coords[v2],
      coords[v3]);
  double vec30[] = SUB(
      coords[v3], coords[v0]);
  angles[0] = angleBetweenVecs(vec01, vec12, -1);
  angles[1] = angleBetweenVecs(vec12, vec23, -1);
  angles[2] = angleBetweenVecs(vec23, vec30, -1);
  angles[3] = angleBetweenVecs(vec30, vec01, -1);

  double vecB[] =
    { 0.25 * (coords[v0][0] - coords[v1][0] - coords[v2][0] + coords[v3][0]),
	0.25 * (coords[v0][1] - coords[v1][1] - coords[v2][1] + coords[v3][1]),
	0.25 * (coords[v0][2] - coords[v1][2] - coords[v2][2] + coords[v3][2]) };
  double vecC[] =
    { 0.25 * (coords[v0][0] + coords[v1][0] - coords[v2][0] - coords[v3][0]),
	0.25 * (coords[v0][1] + coords[v1][1] - coords[v2][1] - coords[v3][1]),
	0.25 * (coords[v0][2] + coords[v1][2] - coords[v2][2] - coords[v3][2]) };
  double junk[] = CROSS(vecB, vecC);
  normal[0] = junk[0];
  normal[1] = junk[1];
  normal[2] = junk[2];
}

static double
tetVolume(const double coordsA[3], const double coordsB[3],
	  const double coordsC[3], const double coordsD[3])
{
  static const double sixth = 1. / 6.;
  double vec1[] = SUB(coordsB, coordsA);
  double vec2[] = SUB(coordsC, coordsA);
  double vec3[] = SUB(coordsD, coordsA);
  double cross[] = CROSS(vec2, vec3);
  double vol = DOT(cross, vec1) * sixth;
  return vol;
}

static double
tetVolume(const double allCoords[][3], const GMGW_int vertA,
	  const GMGW_int vertB, const GMGW_int vertC, const GMGW_int vertD)
{
  return tetVolume(allCoords[vertA], allCoords[vertB], allCoords[vertC],
		   allCoords[vertD]);
}

static double
tetVolume(const double allCoords[][3], const GMGW_int vertA,
	  const GMGW_int vertB, const GMGW_int vertC, const double coordsD[])
{
  return tetVolume(allCoords[vertA], allCoords[vertB], allCoords[vertC],
		   coordsD);
}

static double
pyrVolume(const double allCoords[][3], const GMGW_int vert0,
	  const GMGW_int vert1, const GMGW_int vert2, const GMGW_int vert3,
	  const double coords4[])
{
  // As per VTK file format, verts 0-3 are the base of the pyramid, in
  // right-handed cyclic order.

  // This calculation is exact for a pyramid with a bi-linear base.
  // Straightforward in principle: a coordinate transformation from a canonical
  // pyramid to the real one, then to integrate 1 over the physical pyramid,
  // you integrate det(jacobian) over the canonical pyramid.  In the end,
  // it works out to be exactly the triple product of three vectors.

  double vecA[] =
    { (allCoords[vert0][0] + allCoords[vert1][0] + allCoords[vert2][0]
	+ allCoords[vert3][0]) * 0.25, (allCoords[vert0][1]
	+ allCoords[vert1][1] + allCoords[vert2][1] + allCoords[vert3][1])
	* 0.25, (allCoords[vert0][2] + allCoords[vert1][2] + allCoords[vert2][2]
	+ allCoords[vert3][2]) * 0.25 };
  double vecB[] =
    { (allCoords[vert0][0] - allCoords[vert1][0] - allCoords[vert2][0]
	+ allCoords[vert3][0]) * 0.25, (allCoords[vert0][1]
	- allCoords[vert1][1] - allCoords[vert2][1] + allCoords[vert3][1])
	* 0.25, (allCoords[vert0][2] - allCoords[vert1][2] - allCoords[vert2][2]
	+ allCoords[vert3][2]) * 0.25 };
  double vecC[] =
    { (allCoords[vert0][0] + allCoords[vert1][0] - allCoords[vert2][0]
	- allCoords[vert3][0]) * 0.25, (allCoords[vert0][1]
	+ allCoords[vert1][1] - allCoords[vert2][1] - allCoords[vert3][1])
	* 0.25, (allCoords[vert0][2] + allCoords[vert1][2] - allCoords[vert2][2]
	- allCoords[vert3][2]) * 0.25 };
  double vecE[] =
    { coords4[0] - vecA[0], coords4[1] - vecA[1], coords4[2] - vecA[2] };

  double result = (+vecB[0] * vecC[1] * vecE[2] + vecB[1] * vecC[2] * vecE[0]
      + vecB[2] * vecC[0] * vecE[1] - vecB[0] * vecC[2] * vecE[1]
      - vecB[1] * vecC[0] * vecE[2] - vecB[2] * vecC[1] * vecE[0]) * 4. / 3.;

  return result;
}

static double
pyrVolume(const double allCoords[][3], const GMGW_int vertA,
	  const GMGW_int vertB, const GMGW_int vertC, const GMGW_int vertD,
	  const GMGW_int vertE)
{
  return pyrVolume(allCoords, vertA, vertB, vertC, vertD, allCoords[vertE]);
}

static double
prismVolume(const double allCoords[][3], const GMGW_int vertA,
	    const GMGW_int vertB, const GMGW_int vertC, const GMGW_int vertD,
	    const GMGW_int vertE, const GMGW_int vertF)
{
  static const double sixth = 1. / 6.;
  double result = 0;
  // As per VTK file format, verts ABC form a ring at the bottom,
  // and verts DEF form a ring at the top.  (Orientation in the VTK docs
  // seems backwards compared to all other elements, and both CGNS-to-VTK
  // and my Ugrid-to-VTK converters write properly oriented prisms.  So
  // this routine reflects that, even though it'll puke on "proper" VTK
  // files.

  // Instead of splitting the prism up into three tets (which can be done,
  // but has 8 (=2^3) possibilities, depending on diagonals), I'm going to
  // connect each quad and each tri to the centroid, (technically, the
  // average of the vertex locations) and add their volumes.

  double centroid[] =
    { (allCoords[vertA][0] + allCoords[vertB][0] + allCoords[vertC][0]
	+ allCoords[vertD][0] + allCoords[vertE][0] + allCoords[vertF][0])
	* sixth,

    (allCoords[vertA][1] + allCoords[vertB][1] + allCoords[vertC][1]
	+ allCoords[vertD][1] + allCoords[vertE][1] + allCoords[vertF][1])
	* sixth, (allCoords[vertA][2] + allCoords[vertB][2]
	+ allCoords[vertC][2] + allCoords[vertD][2] + allCoords[vertE][2]
	+ allCoords[vertF][2]) * sixth };

  // First quad face: BADE
  result += pyrVolume(allCoords, vertB, vertA, vertD, vertE, centroid);
  // Second quad face: CBEF
  result += pyrVolume(allCoords, vertC, vertB, vertE, vertF, centroid);
  // Third quad face: ACFD
  result += pyrVolume(allCoords, vertA, vertC, vertF, vertD, centroid);
  // First tri face: ABC
  result += tetVolume(allCoords, vertA, vertB, vertC, centroid);
  // Second tri face: FED
  result += tetVolume(allCoords, vertF, vertE, vertD, centroid);
  return result;
}

static double
hexVolume(const double allCoords[][3], const GMGW_int vertA,
	  const GMGW_int vertB, const GMGW_int vertC, const GMGW_int vertD,
	  const GMGW_int vertE, const GMGW_int vertF, const GMGW_int vertG,
	  const GMGW_int vertH)
{
  double result = 0;
  // Verts ABCD form a ring at the bottom, and verts EFGH form a ring at
  // the top...

  // Instead of splitting the prism up into five or six tets, I'm going to
  // connect each quad to the centroid, (technically, the
  // average of the vertex locations) and add the volumes of six pyramids.

  double centroid[] =
    { (allCoords[vertA][0] + allCoords[vertB][0] + allCoords[vertC][0]
	+ allCoords[vertD][0] + allCoords[vertE][0] + allCoords[vertF][0]
	+ allCoords[vertG][0] + allCoords[vertH][0]) * 0.125,
	(allCoords[vertA][1] + allCoords[vertB][1] + allCoords[vertC][1]
	    + allCoords[vertD][1] + allCoords[vertE][1] + allCoords[vertF][1]
	    + allCoords[vertG][1] + allCoords[vertH][1]) * 0.125,
	(allCoords[vertA][2] + allCoords[vertB][2] + allCoords[vertC][2]
	    + allCoords[vertD][2] + allCoords[vertE][2] + allCoords[vertF][2]
	    + allCoords[vertG][2] + allCoords[vertH][2]) * 0.125 };

  // First quad face: ABCD
  result += pyrVolume(allCoords, vertA, vertB, vertC, vertD, centroid);
  // Second quad face: HGFE
  result += pyrVolume(allCoords, vertH, vertG, vertF, vertE, centroid);
  // Third quad face: BAEF
  result += pyrVolume(allCoords, vertB, vertA, vertE, vertF, centroid);
  // Fourth quad face: CBFG
  result += pyrVolume(allCoords, vertC, vertB, vertF, vertG, centroid);
  // Fifth quad face: DCGH
  result += pyrVolume(allCoords, vertD, vertC, vertG, vertH, centroid);
  // Sixth quad face: ADHE
  result += pyrVolume(allCoords, vertA, vertD, vertH, vertE, centroid);
  return result;
}

double
cellVolume(const double coords[][3], const GMGW_int nConn,
	   const GMGW_int connect[])
{
  double thisVol = -1;
// Volume checks
  switch (nConn)
    {
    case 4:
      // This is a tet.
      thisVol = tetVolume(coords, connect[0], connect[1], connect[2],
			  connect[3]);
      break;
    case 5:
      // This is a pyramid.
      thisVol = pyrVolume(coords, connect[0], connect[1], connect[2],
			  connect[3], connect[4]);
      break;
    case 6:
      // This is a prism.
      thisVol = prismVolume(coords, connect[0], connect[1], connect[2],
			    connect[3], connect[4], connect[5]);
      break;
    case 8:
      // This is a hex.
      thisVol = hexVolume(coords, connect[0], connect[1], connect[2],
			  connect[3], connect[4], connect[5], connect[6],
			  connect[7]);
      break;
    default:
      assert(0);
      break;
    }
  return thisVol;
}

static void
findQuadNonPlanarity(const double coords[][3], const GMGW_int v0,
		     const GMGW_int v1, const GMGW_int v2, const GMGW_int v3,
		     GMGW_int quadDistortion[])
{
  static constexpr double scaling = 2 * M_SQRT2 * pow(3., 1.75); // about 20

// It would be really nice not to have to go through all this;
// importing the calcs for this quality measure might be easier...
  double vol = fabs(tetVolume(coords, v0, v1, v2, v3));
  double area013 = triArea(coords, v0, v1, v3);
  double area123 = triArea(coords, v1, v2, v3);
  double area203 = triArea(coords, v2, v0, v3);
  double area012 = triArea(coords, v0, v1, v2);

  double totalArea = area013 + area123 + area203 + area012;

  double value = scaling * vol / pow(totalArea, 1.5);
  addToBins(value, nDistortBins, distortBinBdrys, quadDistortion);
}

static void
analyzeTetQuality(const double coords[][3], const GMGW_int connect[],
		  GMGW_int triFaceAngles[], GMGW_int dihedralsTriTri[])
{
  double angles[12];
  double norm012[3], norm031[3], norm132[3], norm230[3];
  analyzeTri(coords, connect[0], connect[1], connect[2], angles + 0, norm012);
  analyzeTri(coords, connect[0], connect[3], connect[1], angles + 3, norm031);
  analyzeTri(coords, connect[1], connect[3], connect[2], angles + 6, norm132);
  analyzeTri(coords, connect[2], connect[3], connect[0], angles + 9, norm230);

  for (GMGW_int ii = 0; ii < 12; ii++) {
    addToBins(angles[ii], nAngleBins, angleBinBdrys, triFaceAngles);
  }

  angles[0] = angleBetweenVecs(norm012, norm031, -1);
  angles[1] = angleBetweenVecs(norm012, norm132, -1);
  angles[2] = angleBetweenVecs(norm012, norm230, -1);
  angles[3] = angleBetweenVecs(norm031, norm132, -1);
  angles[4] = angleBetweenVecs(norm031, norm230, -1);
  angles[5] = angleBetweenVecs(norm132, norm230, -1);
  for (GMGW_int ii = 0; ii < 6; ii++) {
    addToBins(angles[ii], nAngleBins, angleBinBdrys, dihedralsTriTri);
  }
}

static void
analyzePyrQuality(const double coords[][3], const GMGW_int connect[],
		  GMGW_int triFaceAngles[], GMGW_int quadFaceAngles[],
		  GMGW_int quadDistortion[], GMGW_int dihedralsQuadTri[],
		  GMGW_int dihedralsTriTri[])
{
  // Verts 0-3 are the base of the pyramid, in RH cyclic order, and vert 4 is
  // the apex.
  const GMGW_int& v0 = connect[0];
  const GMGW_int& v1 = connect[1];
  const GMGW_int& v2 = connect[2];
  const GMGW_int& v3 = connect[3];
  const GMGW_int& v4 = connect[4];

  // There are sixteen face angles here.
  double angles[16];
  // A normal for each face.
  double norm014[3], norm124[3], norm234[3], norm304[3], norm0123[3];
  analyzeTri(coords, v1, v0, v4, angles + 0, norm014);
  analyzeTri(coords, v2, v1, v4, angles + 3, norm124);
  analyzeTri(coords, v3, v2, v4, angles + 6, norm234);
  analyzeTri(coords, v0, v3, v4, angles + 9, norm304);
  analyzeQuad(coords, v0, v1, v2, v3, angles + 12, norm0123);

  for (GMGW_int ii = 0; ii < 12; ii++) {
    addToBins(angles[ii], nAngleBins, angleBinBdrys, triFaceAngles);
  }
  for (GMGW_int ii = 12; ii < 16; ii++) {
    addToBins(angles[ii], nAngleBins, angleBinBdrys, quadFaceAngles);
  }
  findQuadNonPlanarity(coords, v0, v1, v2, v3, quadDistortion);

  // Eight dihedral angles
  angles[0] = angleBetweenVecs(norm014, norm124, -1);
  angles[1] = angleBetweenVecs(norm124, norm234, -1);
  angles[2] = angleBetweenVecs(norm234, norm304, -1);
  angles[3] = angleBetweenVecs(norm304, norm014, -1);
  angles[4] = angleBetweenVecs(norm014, norm0123, -1);
  angles[5] = angleBetweenVecs(norm124, norm0123, -1);
  angles[6] = angleBetweenVecs(norm234, norm0123, -1);
  angles[7] = angleBetweenVecs(norm304, norm0123, -1);
  for (GMGW_int ii = 0; ii < 4; ii++) {
    addToBins(angles[ii], nAngleBins, angleBinBdrys, dihedralsTriTri);
  }
  for (GMGW_int ii = 5; ii < 8; ii++) {
    addToBins(angles[ii], nAngleBins, angleBinBdrys, dihedralsQuadTri);
  }
}

static void
analyzePrismQuality(const double coords[][3], const GMGW_int connect[],
		    GMGW_int triFaceAngles[], GMGW_int quadFaceAngles[],
		    GMGW_int quadDistortion[], GMGW_int dihedralsQuadQuad[],
		    GMGW_int dihedralsQuadTri[])
{
  // Verts 0-3 are the base of the pyramid, in RH cyclic order, and vert 4 is
  // the apex.
  const GMGW_int& vertA = connect[0];
  const GMGW_int& vertB = connect[1];
  const GMGW_int& vertC = connect[2];
  const GMGW_int& vertD = connect[3];
  const GMGW_int& vertE = connect[4];
  const GMGW_int& vertF = connect[5];

  // There are eighteen face angles here.
  double angles[18];
  // A normal for each face.
  double normABC[3], normFED[3], normBADE[3], normCBEF[3], normACFD[3];
  analyzeTri(coords, vertA, vertB, vertC, angles + 0, normABC);
  analyzeTri(coords, vertF, vertE, vertD, angles + 3, normFED);
  analyzeQuad(coords, vertB, vertA, vertD, vertE, angles + 6, normBADE);
  analyzeQuad(coords, vertC, vertB, vertE, vertF, angles + 10, normCBEF);
  analyzeQuad(coords, vertA, vertC, vertF, vertD, angles + 14, normACFD);

  for (GMGW_int ii = 0; ii < 6; ii++) {
    addToBins(angles[ii], nAngleBins, angleBinBdrys, triFaceAngles);
  }
  for (GMGW_int ii = 7; ii < 18; ii++) {
    addToBins(angles[ii], nAngleBins, angleBinBdrys, quadFaceAngles);
  }

  findQuadNonPlanarity(coords, vertB, vertA, vertD, vertE, quadDistortion);
  findQuadNonPlanarity(coords, vertC, vertB, vertE, vertF, quadDistortion);
  findQuadNonPlanarity(coords, vertA, vertC, vertF, vertD, quadDistortion);

  // Nine dihedral angles
  angles[0] = angleBetweenVecs(normABC, normBADE, -1);
  angles[1] = angleBetweenVecs(normABC, normCBEF, -1);
  angles[2] = angleBetweenVecs(normABC, normACFD, -1);
  angles[3] = angleBetweenVecs(normABC, normBADE, -1);
  angles[4] = angleBetweenVecs(normABC, normCBEF, -1);
  angles[5] = angleBetweenVecs(normABC, normACFD, -1);
  angles[6] = angleBetweenVecs(normACFD, normBADE, -1);
  angles[7] = angleBetweenVecs(normBADE, normCBEF, -1);
  angles[8] = angleBetweenVecs(normCBEF, normACFD, -1);
  for (GMGW_int ii = 0; ii < 6; ii++) {
    addToBins(angles[ii], nAngleBins, angleBinBdrys, dihedralsQuadTri);
  }
  for (GMGW_int ii = 7; ii < 9; ii++) {
    addToBins(angles[ii], nAngleBins, angleBinBdrys, dihedralsQuadQuad);
  }
}

static void
analyzeHexQuality(const double coords[][3], const GMGW_int connect[],
		  GMGW_int quadFaceAngles[], GMGW_int quadDistortion[],
		  GMGW_int dihedralsQuadQuad[])
{
  // Verts 0-3 are the base of the pyramid, in RH cyclic order, and vert 4 is
  // the apex.
  const GMGW_int& vertA = connect[0];
  const GMGW_int& vertB = connect[1];
  const GMGW_int& vertC = connect[2];
  const GMGW_int& vertD = connect[3];
  const GMGW_int& vertE = connect[4];
  const GMGW_int& vertF = connect[5];
  const GMGW_int& vertG = connect[6];
  const GMGW_int& vertH = connect[7];

  // There are twenty-four face angles here.
  double angles[24];
  // A normal for each face.
  double normABCD[3], normHGFE[3], normBAEF[3], normCBFG[3], normDCGH[3],
      normADHE[3];
  analyzeQuad(coords, vertA, vertB, vertC, vertD, angles + 0, normABCD);
  analyzeQuad(coords, vertH, vertG, vertF, vertE, angles + 4, normHGFE);
  analyzeQuad(coords, vertB, vertA, vertE, vertF, angles + 8, normBAEF);
  analyzeQuad(coords, vertC, vertB, vertF, vertG, angles + 12, normCBFG);
  analyzeQuad(coords, vertD, vertC, vertG, vertH, angles + 16, normDCGH);
  analyzeQuad(coords, vertA, vertD, vertH, vertE, angles + 20, normADHE);

  for (GMGW_int ii = 0; ii < 24; ii++) {
    addToBins(angles[ii], nAngleBins, angleBinBdrys, quadFaceAngles);
  }

  findQuadNonPlanarity(coords, vertA, vertB, vertC, vertD, quadDistortion);
  findQuadNonPlanarity(coords, vertH, vertG, vertF, vertE, quadDistortion);
  findQuadNonPlanarity(coords, vertB, vertA, vertE, vertF, quadDistortion);
  findQuadNonPlanarity(coords, vertC, vertB, vertF, vertG, quadDistortion);
  findQuadNonPlanarity(coords, vertD, vertC, vertG, vertH, quadDistortion);
  findQuadNonPlanarity(coords, vertA, vertD, vertH, vertE, quadDistortion);

  // Twelve dihedral angles
  angles[0] = angleBetweenVecs(normABCD, normBAEF, -1);
  angles[1] = angleBetweenVecs(normABCD, normCBFG, -1);
  angles[2] = angleBetweenVecs(normABCD, normDCGH, -1);
  angles[3] = angleBetweenVecs(normABCD, normADHE, -1);
  angles[4] = angleBetweenVecs(normHGFE, normBAEF, -1);
  angles[5] = angleBetweenVecs(normHGFE, normCBFG, -1);
  angles[6] = angleBetweenVecs(normHGFE, normDCGH, -1);
  angles[7] = angleBetweenVecs(normHGFE, normADHE, -1);
  angles[8] = angleBetweenVecs(normBAEF, normCBFG, -1);
  angles[9] = angleBetweenVecs(normCBFG, normDCGH, -1);
  angles[10] = angleBetweenVecs(normDCGH, normADHE, -1);
  angles[11] = angleBetweenVecs(normADHE, normBAEF, -1);
  for (GMGW_int ii = 0; ii < 12; ii++) {
    addToBins(angles[ii], nAngleBins, angleBinBdrys, dihedralsQuadQuad);
  }
}

void
analyzeCellQuality(const double coords[][3], const GMGW_int nConn,
		   const GMGW_int connect[], GMGW_int triFaceAngles[],
		   GMGW_int quadFaceAngles[], GMGW_int quadDistortion[],
		   GMGW_int dihedralsQuadQuad[], GMGW_int dihedralsQuadTri[],
		   GMGW_int dihedralsTriTri[])
{
  // For all faces, find their faces angle and bin the results.
  // For quad faces, find their non-planarity; bin the result.
  // For all edges, find the dihedral angle and bin the results.
  switch (nConn)
    {
    case 4:
      analyzeTetQuality(coords, connect, triFaceAngles, dihedralsTriTri);
      break;
    case 5:
      analyzePyrQuality(coords, connect, triFaceAngles, quadFaceAngles,
			quadDistortion, dihedralsQuadTri, dihedralsTriTri);
      break;
    case 6:
      analyzePrismQuality(coords, connect, triFaceAngles, quadFaceAngles,
			  quadDistortion, dihedralsQuadQuad, dihedralsQuadTri);
      break;
    case 8:
      analyzeHexQuality(coords, connect, quadFaceAngles, quadDistortion,
			dihedralsQuadQuad);
      break;
    default:
      assert(0);
break;
}
}

void
analyzeBdryFace(const double coords[][3], const GMGW_int nConn,
		const GMGW_int connect[], GMGW_int triFaceAngles[],
GMGW_int quadFaceAngles[], GMGW_int quadDistortion[])
{
double normal[3], angles[4];
assert(nConn == 3 || nConn == 4);
if (nConn == 3) {
analyzeTri(coords, connect[0], connect[1], connect[2], angles, normal);
addToBins(angles[0], nAngleBins, angleBinBdrys, triFaceAngles);
addToBins(angles[1], nAngleBins, angleBinBdrys, triFaceAngles);
addToBins(angles[2], nAngleBins, angleBinBdrys, triFaceAngles);
}
else {
analyzeQuad(coords, connect[0], connect[1], connect[2], connect[3], angles,
		normal);
addToBins(angles[0], nAngleBins, angleBinBdrys, quadFaceAngles);
addToBins(angles[1], nAngleBins, angleBinBdrys, quadFaceAngles);
addToBins(angles[2], nAngleBins, angleBinBdrys, quadFaceAngles);
addToBins(angles[3], nAngleBins, angleBinBdrys, quadFaceAngles);
findQuadNonPlanarity(coords, connect[0], connect[1], connect[2], connect[3],
			 quadDistortion);
}
}

