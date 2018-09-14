/*
 * surface_proj.cxx
 *
 *  Created on: Oct 31, 2016
 *      Author: Mike Jefferies, Pointwise
 *
 *      Modified for integration with other parts of the HiLPW/GMGW tool chain
 *      by Carl Ollivier-Gooch.
 */

// GeodeProjectPoints.cpp
//
#include <string>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <vector>
#include <map>
#include <set>
#include <sstream>

#ifdef HAVE_GEODE
#include <geom/Database.h>
#include <geom/DictionaryAttribute.h>
#include <geom/IsectProjPoint.h>
#include <geom/ProjectionBSPTree.h>
#include <geom/Surface.h>
#include <nmb/CurvedCoedge.h>
#include <nmb/CurvedEdge.h>
#include <nmb/CurvedFace.h>
#include <nmb/CurvedModel.h>
#include <nmb/CurvedSheet.h>
#include <nmb/CurvedVertex.h>
#include <nmb/NativeTopologyReader.h>
#include <nmb/TopologyProjectionBSPTreeWrapper.h>

typedef std::map<const GE::Entity *, std::string> EntityToNameMap;
#endif

struct Counter {
    Counter() : count(0) {}
    Counter(unsigned int c) : count(c) {}
    operator unsigned int & () { return count; }
    operator unsigned int const & () const { return count; }
    unsigned int count;
};

typedef std::map<std::string, Counter> NameToCountMap;

#ifdef HAVE_GEODE
struct Args {
    bool showUsage;
    std::string nmbFilename;
    std::string inFilename;
    GE::Int32 inNumFiles;
    GE::Real64 modelSize;
    GE::Int32 bspMaxLevel;
    GE::Int32 bspMaxObjPerCell;
    GE::Real64 bspMaxProjDist;

    Args() :
        showUsage(false),
        nmbFilename(),
        inFilename(),
        inNumFiles(0),
        modelSize(0.0),
        bspMaxLevel(0),
        bspMaxObjPerCell(0),
        bspMaxProjDist(GE::Tolerance::Infinity)
    {
    }

    bool parse(int argc, char **argv)
    {
        bool result = true;
        for (int i = 1; i < argc && result; ++i) {
            std::string arg = argv[i];
            if (arg == "-h") {
                showUsage = true;
            }
            else if (arg == "-nmb") {
                result = (i + 1 < argc);
                if (result) {
                    nmbFilename = argv[++i];
                }
            }
            else if (arg == "-modelSize") {
                result = (i + 1 < argc);
                if (result) {
                    modelSize = atof(argv[++i]);
                    result = (modelSize >= 10.0 && modelSize <= 100000.0);
                }
            }
            else if (arg == "-bspMaxLevel") {
                result = (i + 1 < argc);
                if (result) {
                    bspMaxLevel = atoi(argv[++i]);
                    result = (bspMaxLevel >= 8 && bspMaxLevel <= 30);
                }
            }
            else if (arg == "-bspMaxObjPerCell") {
                result = (i + 1 < argc);
                if (result) {
                    bspMaxObjPerCell = atoi(argv[++i]);
                    result = (bspMaxObjPerCell >= 4 && bspMaxObjPerCell <= 16);
                }
            }
            else if (arg == "-bspMaxProjDist") {
                result = (i + 1 < argc);
                if (result) {
                    bspMaxProjDist = atof(argv[++i]);
                    result = (bspMaxProjDist > 0.0);
                }
            }
            else {
                std::cerr << "Invalid arg: " << arg << std::endl;
                result = false;
            }
        }

        if (result && !showUsage && nmbFilename.length() == 0) {
            std::cerr << "The -nmb parameter is required" << std::endl;
            result = false;
        }

        return result;
    }

    void printUsage()
    {
        std::cout << "GeodeProjectPoints usage:" << std::endl;
        std::cout << "Required parameters" << std::endl;
        std::cout << "   -nmb file : The nmb file to import for projecting onto." << std::endl;
        std::cout << "Optional parameters" << std::endl;
        std::cout << "   -modelSize num : The model size to use in Geode.  In most cases, this" << std::endl;
        std::cout << "         should match the model size used to generate the nmb file.  Valid" << std::endl;
        std::cout << "         values are in the range [10,100000].  If not given, the Geode default" << std::endl;
        std::cout << "         will be used, which is 1000." << std::endl;
        std::cout << "   -bspMaxLevel num : The max number of levels of the BSP tree. Valid values" << std::endl;
        std::cout << "         are in the range [8, 30].  If not given, the parameter will be based" << std::endl;
        std::cout << "         on the projection entities." << std::endl;
        std::cout << "   -bspMaxObjPerCell num : The max number of objects per cell of the BSP tree." << std::endl;
        std::cout << "         Valid values are in the range [4, 16].  If not given, the parameter" << std::endl;
        std::cout << "         will be based on the projection entities." << std::endl;
        std::cout << "   -bspMaxProjDist num : The max projection distance when projecting a point" << std::endl;
        std::cout << "         to the BSP tree.  Valid values are in the range (0, +inf).  If not" << std::endl;
        std::cout << "         given, the parameter defaults to positive infinity." << std::endl;
        std::cout << "   -h : Print command usage and exit" << std::endl;
    }
};


static std::string getEntityName(const GE::Entity *entity)
{
    // find the first dictionary attribute with sub class PW::Common
    // that has a name key and return the value
    std::string result;
    GE::AttributeRegistry *reg = GE::DictionaryAttribute::Singleton_AttributeRegistry();
    if (entity != 0 && reg != 0) {
        GE::EntityList<GE::DictionaryAttribute> dictAttrs;
        entity->Inquire_Attributes(reg, &dictAttrs.Upcast_NonConst<GE::Attribute>());
        for (GE::Int32 i = 0; i < dictAttrs.Size(); i++) {
            if (dictAttrs[i]->Inquire_SubClass() == "PW::Common") {
                result = dictAttrs[i]->Value("name").ConstData();
                if (result.size()) {
                    break;
                }
            }
        }
    }
    return result;
}

static void addEntityNames(const GE::Entity *entity,
    EntityToNameMap &entToName, NameToCountMap &nameToCount,
    NameToCountMap &nameToIndex)
{
    const GE::Surface *surface = GE::Surface::Downcast(entity);
    const GE::CurvedModel *model = GE::CurvedModel::Downcast(entity);

    int surfInd = 0;
    if (surface) {
        // for surfaces, associate the surface with it's name
        // and initialize count to 0, so it will always be reported
        std::string name = getEntityName(surface);
        entToName[surface] = name;
        nameToCount[name] = 0;
        nameToIndex[name] = surfInd++;
    }
    else if (model) {
        // for models,

        // associate the sheets with the sheet name
        // and initialize count to 0, so it will always be reported
        GE::EntityList<GE::CurvedFace> faces;
        model->Inquire_Faces(&faces);
        for (GE::Int32 i = 0; i < faces.Size(); i++) {
            GE::CurvedSheet *sheet = faces[i]->Inquire_Sheet();
            std::string name = getEntityName(sheet);
            entToName[sheet] = name;
            nameToCount[name] = 0;
            nameToIndex[name] = surfInd++;
        }

        // associate edges with the sheet names, don't initialize count
        GE::EntityList<GE::CurvedEdge> edges;
        model->Inquire_Edges(&edges);
        for (GE::Int32 i = 0; i < edges.Size(); i++) {
            std::set<std::string> names;
            edges[i]->Inquire_Faces(&faces);
            for (GE::Int32 f = 0; f < faces.Size(); f++) {
                names.insert(entToName[faces[f]->Inquire_Sheet()]);
            }

            std::ostringstream edgeDescription;
            edgeDescription << "Edge of";
            std::string prefix = " ";
            std::set<std::string>::const_iterator iter;
            for (iter = names.begin(); iter != names.end(); ++iter) {
                edgeDescription << prefix << *iter;
                prefix = ", ";
            }
            entToName[edges[i]] = edgeDescription.str();
            nameToIndex[edgeDescription.str()] = surfInd++;
        }

        // associate vertices with the sheet names, don't initialize count
        GE::EntityList<GE::CurvedVertex> vertices;
        model->Inquire_Vertices(&vertices);
        for (GE::Int32 i = 0; i < vertices.Size(); i++) {
            std::set<std::string> names;
            vertices[i]->Inquire_Faces(&faces);
            for (GE::Int32 f = 0; f < faces.Size(); f++) {
                names.insert(entToName[faces[f]->Inquire_Sheet()]);
            }

            std::ostringstream cornerDescription;
            cornerDescription << "Corner of";
            std::string prefix = " ";
            std::set<std::string>::const_iterator iter;
            for (iter = names.begin(); iter != names.end(); ++iter) {
                cornerDescription << prefix << *iter;
                prefix = ", ";
            }
            entToName[vertices[i]] = cornerDescription.str();
            nameToIndex[cornerDescription.str()] = surfInd++;
        }
    }
}
#endif

int projectionChecks(const int nBdryVerts,
                     const double bdryCoords[][3], double bdryDist[],
                     int bdrySurf[], std::string nmbFilename)
{
#ifdef HAVE_GEODE
    printf("Computing wall projection distances\n");

    Args args;
//    if (!args.parse(argc, argv) || args.showUsage) {
//        args.printUsage();
//        return 0;
//    }
    args.nmbFilename = nmbFilename;

    GE::Error err = GE::Error::No_errors;

    // set the model size, if requested
    if (args.modelSize > 0.0) {
        GE::Tolerance::SetModelSize(args.modelSize);
    }

    // read entities from the file, into the database
    std::cout << "Importing entities from " << args.nmbFilename << std::endl;
    GE::Database database;
    err = GE::NativeTopologyReader::Read(args.nmbFilename.c_str(), &database);
    if (err != GE::Error::No_errors) {
        std::cerr << "Error occurred reading from file: " << args.nmbFilename << std::endl;
        return err.ToInt();
    }

    // get the entities to project onto
    GE::EntityList<GE::Entity> projectEntities;

    // get the curved models from the database and add them to the project entities
    GE::EntityList<GE::Entity> curvedModels;
    database.InquireByClassID_Entities(GE::CurvedModel::Static_ClassID(), curvedModels);
    projectEntities += curvedModels;

    // get the free surfaces from the database and add them to the project entities
    GE::EntityList<GE::Entity> allSurfaces;
    database.InquireByClassID_Entities(GE::Surface::Static_ClassID(), allSurfaces);
    for (GE::Int32 i = 0; i < allSurfaces.Size(); i++) {
        GE::Surface *surf = GE::Surface::Downcast(allSurfaces[i]);
        if (surf) {
            GE::EntityList<GE::CurvedFace> faces;
            GE::CurvedFace::Surface_Inquire_Faces(surf, &faces);
            if (faces.Size() == 0) {
                projectEntities += surf;
            }
        }
    }

    // make sure we have entities to project to
    if (projectEntities.Size() == 0) {
        std::cerr << "The NMB file did not contain any projection entities" << std::endl;
        return 0;
    }

    EntityToNameMap entityToNameMap;
    NameToCountMap hitCounts, surfIndices;

    // add the projection entities to a bsp tree and determine names for the entities
    std::cout << "Adding entities to BSP tree." << std::endl;
    GE::ProjectionBSPTree bspTree;
    for (GE::Int32 i = 0; i < projectEntities.Size(); i++) {
        err = GE::TopologyProjectionBSPTreeWrapper::Add_Entity(&bspTree, projectEntities[i]);
        if (err != GE::Error::No_errors) {
            std::cerr << "Error adding an entity to BSP tree" << std::endl;
            return err.ToInt();
        }
        addEntityNames(projectEntities[i], entityToNameMap, hitCounts,
                       surfIndices);
    }

    // build the bsp tree
    std::cout << "Building BSP tree." << std::endl;
    err = bspTree.Build_BSPTree(args.bspMaxLevel, args.bspMaxObjPerCell);
    if (err != GE::Error::No_errors) {
        std::cerr << "Error occurred building the BSP tree" << std::endl;
        return err.ToInt();
    }

    // now read the points to project from the input files and project the
    // points onto the entities in the BSP tree keeping track of entity hit
    // counts, distance, squared distance and max distance to the BSP tree

    GE::Real64 sum = 0.0, sum2 = 0.0, max = 0.0;

    GE::Array<GE::Vector3D> points;
    for (int ii = 0; ii < nBdryVerts; ii++) {
        points.Append(GE::Vector3D(bdryCoords[ii][0], bdryCoords[ii][1],
                                   bdryCoords[ii][2]));
    }

    GE::Int32 numPoints = points.Size();

    GE::Real64 stepsPerPercent = points.Size() / 100.0;
    GE::Real64 stepsUntilPrint = stepsPerPercent;
    for (GE::Int32 i = 0; i < points.Size(); i++) {
        stepsUntilPrint -= 1.0;
        if (stepsUntilPrint < 0.0) {
            std::cout << "." << std::flush;
            stepsUntilPrint += stepsPerPercent;
        }
        GE::Real64 dist = 0;
        int surf = -1;

        // Don't bother with distant points, or those on the symmetry plane.
        if (points[i].Y() == 0 && hypot(points[i].X(), points[i].Z()) >= 4000) {
        	surf = -1;
        	dist = -1;
        }
        else {
            // project the point
            GE::Interval bounds(0, args.bspMaxProjDist);
            bool boundsExceeded = false;
            GE::IsectProjPoint projPoint;
            err = bspTree.Compute_CoordMinimumDistance(points[i], &bounds,
                                                       &boundsExceeded,
                                                       &projPoint);
            if (err != GE::Error::No_errors) {
                std::cerr << "Error projecting point to BSP tree" << std::endl;
                return err.ToInt();
            }
            // classify the projection result
            if (projPoint.End1.subEntity == 0) {
                if (projPoint.End1.entity != 0) {
                    hitCounts[entityToNameMap[projPoint.End1.entity]]++;
                    surf = surfIndices[entityToNameMap[projPoint.End1.entity]];
                }
            }
            else {
                const GE::Entity *hitEntity = projPoint.End1.subEntity;

                const GE::CurvedVertex *vert = GE::CurvedVertex::Downcast(
                        hitEntity);
                if (vert && 0 != vert->Inquire_OwningTopologyDimension()) {
                    // this is not a "hard" vertex, so consider an edge of the vertex the hit entity
                    hitEntity = vert->Inquire_AnyEdge();
                }

                const GE::CurvedEdge *edge = GE::CurvedEdge::Downcast(
                        hitEntity);
                if (edge && 1 != edge->Inquire_OwningTopologyDimension()) {
                    // this is not a "hard" edge, so consider a face of the edge as the hit entity
                    hitEntity = edge->Inquire_FirstCoedge()->Inquire_Face();
                }

                const GE::CurvedFace *face = GE::CurvedFace::Downcast(
                        hitEntity);
                if (face) {
                    // consider the sheet of the face as the hit entity
                    hitEntity = face->Inquire_Sheet();
                }

                hitCounts[entityToNameMap[hitEntity]]++;
                surf = surfIndices[entityToNameMap[hitEntity]];
            }
            // record the distance and updated the sums and max
            dist = (boundsExceeded ? args.bspMaxProjDist : projPoint.Distance);
        }
        sum += dist;
        sum2 += dist * dist;
        max = (dist > max ? dist : max);
        bdryDist[i] = dist;
        bdrySurf[i] = surf;
    }

    std::cout << std::endl;

    std::ostringstream checkWidth;
    checkWidth << numPoints;
    int width = checkWidth.str().size() + 1;
    std::cout << std::right << std::setfill(' ') << std::endl;
    NameToCountMap::const_iterator hit;
    for (hit = hitCounts.begin(); hit != hitCounts.end(); ++hit) {
    std::cout << std::setw(width) << hit->second << " points on " << hit->first
	<< " (" << surfIndices[hit->first].count << ")" << std::endl;
    }

    std::cout << std::setw(width) << numPoints << " total points." << std::endl;

    std::cout << std::endl;
    std::cout << "Error stats:";
    std::cout.precision(17);
    std::cout << "  L1: " << (sum / numPoints);
    std::cout << "  L2: " << sqrt(sum2 / numPoints);
    std::cout << "  LInf: " << max << std::endl;

    return 0;
#else
    return 1; // Failure; that is, we didn't produce bdry distances.
#endif
}



