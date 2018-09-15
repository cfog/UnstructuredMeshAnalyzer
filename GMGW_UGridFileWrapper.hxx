/*
 * GMGW_UgridFileWrapper.hxx
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


#ifndef GMGW_UGRIDFILEWRAPPER_HXX_
#define GMGW_UGRIDFILEWRAPPER_HXX_

#include "GMGW_FileWrapper.hxx"

class UGridFileWrapper : public FileWrapper {
    bool bigEndian;
    mutable unsigned long whichCell;
public:
    UGridFileWrapper(const char fileNameBase[], const char infix[]);
    virtual ~UGridFileWrapper() {}

    virtual void scanFile();

    virtual void seekStartOfConnectivity() const;
    virtual void getNextCellConnectivity(int& nConn, unsigned int connect[]) const;
    virtual void getNextVertexCoords(double& x, double& y, double &z) const;
private:
    unsigned int convertToInt(const unsigned char raw[4]) const;
    double convertToDouble(const unsigned char raw[8]) const;
};

#endif /* GMGW_UGRIDFILEWRAPPER_HXX_ */
