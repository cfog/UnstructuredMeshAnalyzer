/*
 * GMGW_TAUFileWrapper.hxx
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

#ifndef GMGW_TAUFILEWRAPPER_HXX_
#define GMGW_TAUFILEWRAPPER_HXX_

#include "GMGW_FileWrapper.hxx"
#include "primgrid.h"


class TAUFileWrapper : public FileWrapper {
public:
    TAUFileWrapper(const char fileNameBase[]);
    virtual ~TAUFileWrapper();

    virtual void scanFile();

    virtual void getNextCellConnectivity(GMGW_int& nConn, GMGW_int connect[]) const;
    virtual void getNextVertexCoords(double& x, double& y, double &z) const ;
  virtual void
  seekStartOfCoords() const;
  virtual void
  seekStartOfConnectivity() const;
private:
//    TauPrimGrid *m_pg;
  mutable long int m_point_idx=0;
  mutable long int m_cell_idx=0;
};

#endif /* GMGW_TAUFILEWRAPPER_HXX_ */
