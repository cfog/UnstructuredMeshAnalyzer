/*
 * GMGW_VTKReader.hxx
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

#ifndef GMGW_VTKFILEWRAPPER_HXX_
#define GMGW_VTKFILEWRAPPER_HXX_

#include "GMGW_FileWrapper.hxx"

class VTKFileWrapper : public FileWrapper {
public:
    VTKFileWrapper(const char fileNameBase[]);
    virtual ~VTKFileWrapper() {}

    virtual void scanFile();

    virtual void getNextCellConnectivity(int& nConn, unsigned int connect[]) const;
    virtual void getNextVertexCoords(double& x, double& y, double &z) const;
private:
    void consumeLine() const;
    void skipHeader() const;
};

#endif /* GMGW_VTKFILEWRAPPER_HXX_ */
