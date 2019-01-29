/*
 * TAUFileWrapper.cxx
 *
 *  Created on: Jan 07, 2019
 *      Author: jww
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

#include <GMGW_unstr.hxx>
#include <GMGW_TAUFileWrapper.hxx>
#include "check_malloc.h"
#include "stream_interface.h"
#include "tau_io_names.h"
#include "primgrid_io_pdat.h"
#include <cstring>

static long int m_point_idx=0;
static long int m_cell_idx=0;


TAUFileWrapper::TAUFileWrapper(const char fileNameBase[]) :
    FileWrapper()//, m_pg(nullptr) , m_point_idx(0), m_cell_idx(0)
{
    char fileNameIn[1000];
    char *stream = "primgrid_read";
//    m_pg = init_primgrid();
    snprintf(fileNameIn, 1000, "%s", fileNameBase);
    stream_open(stream, STREAM_BUFFER_ATTS);

    stream_open_file(stream, STREAM_NETCDF, fileNameIn, STREAM_READ_ONLY);
}

TAUFileWrapper::~TAUFileWrapper(void)
{
  char *stream = "primgrid_read";
  stream_close(stream);
//  free_primgrid(&m_pg);
}

void TAUFileWrapper::scanFile()
{
    char *stream = "primgrid_read";

    coordsStart = 0;
    nVerts = stream_get_dim(stream, NO_OF_POINTS);

    connectStart = 0;
    nTets      = stream_get_dim(stream, NO_OF_TETRAS);
    nPrisms    = stream_get_dim(stream, NO_OF_PRISMS);
    nHexes     = stream_get_dim(stream, NO_OF_HEXAS);	
    nPyrs      = stream_get_dim(stream, NO_OF_PYRAS);
    nBdryTris  = stream_get_dim(stream, NO_OF_SURFACETRIANGLES);
    nBdryQuads = stream_get_dim(stream, NO_OF_SURFACEQUADRILATERALS);
    nCells = nTets + nPrisms + nHexes + nPyrs + nBdryTris + nBdryQuads;

    m_cellTypes = new char[nCells];
    m_isBdryFace = new bool[nCells];

    unsigned int i_cell = 0;
    for (unsigned int ii = 0; ii < nTets; ii++,i_cell++) {
        m_isBdryFace[i_cell] = false;
        m_cellTypes[i_cell] = TET;
    }
    for (unsigned int ii = 0; ii < nPrisms; ii++,i_cell++) {
        m_isBdryFace[i_cell] = false;
        m_cellTypes[i_cell] = PRISM;
    }
    for (unsigned int ii = 0; ii < nHexes; ii++,i_cell++) {
        m_isBdryFace[i_cell] = false;
        m_cellTypes[i_cell] = HEX;
    }
    for (unsigned int ii = 0; ii < nPyrs; ii++,i_cell++) {
        m_isBdryFace[i_cell] = false;
        m_cellTypes[i_cell] = PYRAMID;
    }
    for (unsigned int ii = 0; ii < nBdryTris; ii++,i_cell++) {
        m_isBdryFace[i_cell] = true;
        m_cellTypes[i_cell] = BDRY_TRI;
    }
    for (unsigned int ii = 0; ii < nBdryQuads; ii++,i_cell++) {
        m_isBdryFace[i_cell] = true;
        m_cellTypes[i_cell] = BDRY_QUAD;
    }

    identifyBdryVerts();

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

void TAUFileWrapper::getNextVertexCoords(double& x, double& y, double& z) const
{
    char *st = "primgrid_read";
    size_t tmp_offset = m_point_idx;
    size_t tmp_pnts = 1;

    QueVarFormat format = stream_array_format(&tmp_offset, &tmp_pnts);
    TauDouble (*xx)[3] = NULL;
    xx          = stream_get_double3_format_var(st, POINT_COORDINATES,
                                                    &format);
    
    if(xx == NULL) /* if coordinates are stored with old names */
    {
       char *xyzname[3] = { POINT_X_COORDINATES,
                            POINT_Y_COORDINATES,
                            POINT_Z_COORDINATES };
       int n_points = static_cast<int> (nVerts);		    
       read_xyz(st, &xx, &n_points, xyzname, &format);
    }

    x = xx[0][0];    
    y = xx[0][1];    
    z = xx[0][2];    
    check_free(xx);
    m_point_idx++;
}

void TAUFileWrapper::getNextCellConnectivity(GMGW_int& nConn, GMGW_int connect[]) const
{
    char *stream = "primgrid_read";
    assert (m_cell_idx < nCells);
    
    long int cell_offset = m_cell_idx;
    if (cell_offset >=  nTets) {
      cell_offset-= nTets;
      if (cell_offset >=  nPrisms) {
        cell_offset-= nPrisms;
        if (cell_offset >=  nHexes) {
	  cell_offset-= nHexes;
          if (cell_offset >=  nPyrs) {
	    cell_offset-= nPyrs;
            if (cell_offset >=  nBdryTris) {
	      cell_offset-= nBdryTris;
	      assert (cell_offset < nBdryQuads);
	    }
	  }
	}
      }
    }
    
    size_t start[2], tmp_offset[2];
    start[0]	  = cell_offset;
    start[1]	  = 0;
    tmp_offset[0] = 1;
    
    TauIndex act[8];
    QueVarFormat format;
    
    int j;
    switch(m_cellTypes[m_cell_idx])
    {
    case TET:
      nConn = tmp_offset[1] = 4;
      format = stream_array_format(start, tmp_offset);

      TauIndex (*act_tet)[4];
      act_tet = stream_get_index4_format_var(stream, TETRA_POINTS, &format);
      memcpy(act,act_tet[0],tmp_offset[1]*sizeof(TauIndex));
      check_free(act_tet);
	
      break;
    case PRISM:
      nConn = tmp_offset[1] = 6;
      format = stream_array_format(start, tmp_offset);

//      act = stream_get_index6_format_var(stream, PRISM_POINTS, &format);
      TauIndex (*act_prism)[6];
      act_prism = stream_get_index6_format_var(stream, PRISM_POINTS, &format);
      memcpy(act,act_prism[0],tmp_offset[1]*sizeof(TauIndex));
      check_free(act_prism);
	
      break;
    case HEX:
      nConn = tmp_offset[1] = 8;
      format = stream_array_format(start, tmp_offset);

//      act = stream_get_index8_format_var(stream, HEXA_POINTS, &format);
      TauIndex (*act_hexa)[8];
      act_hexa = stream_get_index8_format_var(stream, HEXA_POINTS, &format);
      memcpy(act,act_hexa[0],tmp_offset[1]*sizeof(TauIndex));
      check_free(act_hexa);
	
      break;
    case PYRAMID:
      nConn = tmp_offset[1] = 5;
      format = stream_array_format(start, tmp_offset);

//      act = stream_get_index5_format_var(stream, PYRA_POINTS, &format);
      TauIndex (*act_pyra)[5];
      act_pyra = stream_get_index5_format_var(stream, PYRA_POINTS, &format);
      memcpy(act,act_pyra[0],tmp_offset[1]*sizeof(TauIndex));
      check_free(act_pyra);
	
      break;
    case BDRY_TRI:
      nConn = tmp_offset[1] = 3;
      format = stream_array_format(start, tmp_offset);

//      act = stream_get_index3_format_var(stream, SURFTRI_POINTS, &format);
      TauIndex (*act_stri)[3];
      act_stri = stream_get_index3_format_var(stream, SURFTRI_POINTS, &format);
      memcpy(act,act_stri[0],tmp_offset[1]*sizeof(TauIndex));
      check_free(act_stri);
	
      break;
    case BDRY_QUAD:
      nConn = tmp_offset[1] = 4;
      format = stream_array_format(start, tmp_offset);

//      act = stream_get_index4_format_var(stream, SURFQUAD_POINTS, &format);
      TauIndex (*act_squad)[4];
      act_squad = stream_get_index4_format_var(stream, SURFQUAD_POINTS, &format);
      memcpy(act,act_squad[0],tmp_offset[1]*sizeof(TauIndex));
      check_free(act_squad);
	
      break;
    }
    for (j=0; j<nConn; j++)
      connect[j] = act[j];
//    check_free(act);
    m_cell_idx++;
}

void  TAUFileWrapper::seekStartOfCoords() const {m_point_idx = 0;} 
void  TAUFileWrapper::seekStartOfConnectivity() const {m_cell_idx = 0;}
