/*
 * GMGW_VTKReader.hxx
 *
 *  Created on: Jan 27, 2017
 *      Author: cfog
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
