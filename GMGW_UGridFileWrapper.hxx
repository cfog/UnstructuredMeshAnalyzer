/*
 * GMGW_UgridFileWrapper.hxx
 *
 *  Created on: Jan 27, 2017
 *      Author: cfog
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
