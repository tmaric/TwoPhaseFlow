/*---------------------------------------------------------------------------*\

    geom-vof-init

    An OpenFOAM implementation of algorithms for computing volume fractions 
    by intersecting unstructured meshes.
     
    Copyright (C) 2019 Tomislav Maric, TU Darmstadt 

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

Description
    Building geometrical models from the `geophase` library from OpenFOAM 
    types.

Author
    Tomislav Maric maric@mma.tu-darmstadt.de
    Mathematical Modeling and Analysis, TU Darmstadt

\*---------------------------------------------------------------------------*/

#ifndef GeophaseMake_H 
#define GeophaseMake_H 

// STL 
#include <vector>
#include <list>
#include <numeric>
#include <type_traits>

// Geophase
#include "Polygon.H"
#include "Polyhedron.H"
#include "Halfspace.H"
#include "PolyhedronIntersection.H"
#include "Halfspace.H"
#include "Equal.H"
#include "Make.H"
#include "Centroid.H"
#include "Area.H"
#include "Volume.H"
#include "Distance.H"
#include "Intersect.H"
#include "PolyhedronIntersection.H"
#include "WriteVtkPolyData.H"
#include "ReadVtkPolyData.H"
#include "messageStream.H"
#include "vector.H"

// OpenFOAM
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace geophase { 

    // An OpenFOAM vector is a vector so tag it as a vector.
    template<>
    struct tag<Foam::vector>
    {
        using type = vector_tag;
    };

    // An OpenFOAM std::vector polygon is a VectorPolygon.
    using foamVectorPolygon = VectorPolygon<Foam::vector>; 

    // An OpenFOAM std::vector polygon is a VectorPolygon.
    using foamVectorPolyhedron = VectorPolyhedron<foamVectorPolygon>; 

    // An OpenFOAM Halfspace
    using foamHalfspace = Halfspace<Foam::vector>;

    // An OpenFOAM polyhedron intersection. 
    using foamPolyhedronIntersection = PolyhedronIntersection<foamVectorPolyhedron>;
}

using geophase::vectorPolygon;
using geophase::vectorPolyhedron;
using geophase::foamVectorPolygon;
using geophase::foamVectorPolyhedron;
using geophase::make_tetrahedron;
using geophase::foamHalfspace;
using geophase::orient;
using geophase::equal_by_tolerance;
using geophase::signed_distance;
using geophase::intersect_tolerance;
using geophase::foamPolyhedronIntersection;

namespace Foam::GeometricalTransport { 

    template <typename Polygon> 
    geophase::base_tag_t<Polygon, geophase::sequence_polygon_tag>
    make(face const & face, polyMesh const & mesh)
    {
        Polygon polygon;
        const auto& points = mesh.points(); 

        for (auto pointLabel : face)
            polygon.push_back(points[pointLabel]); 

        return polygon;
    }

    //template <typename PolygonSequence>
    //typename std::enable_if
    //<
        //Traits::tag_enabled<PolygonSequence, polygon_sequence_tag>::value, 
        //PolygonSequence
    //>::type 
    //make(const face& meshFace, const polyMesh& mesh, pointField const& displacements)
    //{
        //auto basePolygon = 
            //make<typename PolygonSequence::value_type>(meshFace, mesh); 

        //return make<PolygonSequence>(basePolygon, displacements);
    //}

    //template <typename PointSequence> 
    //typename std::enable_if
    //<
        //Traits::tag_enabled<PointSequence, point_sequence_tag>::value,
        //PointSequence 
    //>::type 
    //make(label const & faceLabel, const polyMesh & mesh)
    //{
        //return make<PointSequence>(mesh.faces()[faceLabel], mesh);
    //}

    //template <typename PolygonSequence> 
    //typename std::enable_if
    //<
        //Traits::tag_enabled<PolygonSequence, polygon_sequence_tag>::value,
        //PolygonSequence
    //>::type 
    //make(label cellLabel, const polyMesh & mesh)
    //{
        //typedef typename PolygonSequence::value_type PolygonType; 
        //PolygonSequence polyhedron; 
        //const auto& own = mesh.faceOwner(); 
        //const auto& cell = mesh.cells()[cellLabel]; 

        //for (auto faceLabel : cell)
        //{
            //auto polygon = make<PolygonType>(faceLabel, mesh);

            //// If the face normal looks into the cell. 
            //if (own[faceLabel] != cellLabel)
            //{
                //std::reverse(polygon.begin(), polygon.end()); 
            //}

            //// Append the poylygon to the PolygonCollection. 
            //polyhedron.push_back(polygon); 
        //}

        //return polyhedron;
    //}

    //template <typename Halfspace>  
    //typename std::enable_if
    //<
        //Traits::tag_enabled<Halfspace, halfspace_tag>::value,
        //Halfspace 
    //>::type 
    //make(const face& meshFace, const pointField& meshPoints)
    //{
        //const point& P0 = meshPoints[meshFace[0]];  
        //const point& P1 = meshPoints[meshFace[1]];  
        //const point& P2 = meshPoints[meshFace[2]];  

        //auto n = (P1 - P0) ^ (P2 - P0);

        //return Halfspace(P0, n); 
    //}
   
    //template <typename Halfspace>  
    //typename std::enable_if
    //<
        //Traits::tag_enabled<Halfspace, halfspace_tag>::value,
        //Halfspace 
    //>::type 
    //make(const label& faceLabel, const fvMesh& mesh)
    //{
        //const auto& Sf = mesh.Sf(); 
        //const auto& Cf = mesh.Cf(); 

        //return Halfspace(Cf[faceLabel], Sf[faceLabel]); 
    //}

    //// Generate a sequence of outward oriented halfspaces from a mesh cell.
    //template <typename HalfspaceCollection> 
    //typename std::enable_if
    //<
        //Traits::tag_enabled<HalfspaceCollection, halfspace_collection_tag>::value,
        //HalfspaceCollection 
    //>::type 
    //make(label cellLabel, const fvMesh& mesh)
    //{
        //HalfspaceCollection halfspaces;

        //typedef typename HalfspaceCollection::value_type Halfspace; 

        //const auto& Sf = mesh.Sf(); 
        //const auto& Cf = mesh.Cf(); 

        //const auto& own = mesh.faceOwner(); 
        //const auto& cell = mesh.cells()[cellLabel]; 
        //const auto& faces = mesh.faces(); 
        //const auto& points = mesh.points(); 

        //forAll(cell, faceI)
        //{
            //auto faceLabel = cell[faceI]; 

            //if (mesh.isInternalFace(faceLabel))
            //{
                //halfspaces.push_back(Halfspace(Cf[cell[faceI]], Sf[cell[faceI]]));            
                //// Flip the face normal if the cell owns the face.
                //// This is to ensure that the set of halfspaces built from a cell
                //// defines a convex set, since it is used as an intersection tool.
                //// Otherwise a flip() of each halfspace is required. TM.
                //if (cellLabel == own[cell[faceI]])
                    //halfspaces.back().flip(); 
            //}
            //else
            //{
                //halfspaces.push_back(Halfspace(-1*faces[faceI].normal(points), points[faces[faceI][0]]));            
            //}
        //}

        //return halfspaces; 
    //}

    //template <typename AABBox> 
    //typename std::enable_if
    //<
        //Traits::tag_enabled<AABBox, aabbox_tag>::value,
        //AABBox 
    //>::type 
    //make(const cell& cell, polyMesh const & mesh)
    //{
        //const auto& meshFaces = mesh.faces(); 
        //const auto& meshPoints = mesh.points(); 
        
        //// Not creating a local copy of cell points for efficiency reasons. TM.
        //std::pair<point,point> extremeCellPoints = std::make_pair(
            //point(GREAT, GREAT, GREAT),
            //point(-GREAT, -GREAT, -GREAT)
        //);

        //forAll(cell, faceI)
        //{
            //const face& meshFace = meshFaces[cell[faceI]]; 
            //forAll(meshFace, pointI)
            //{
                //less_assign(meshPoints[meshFace[pointI]], extremeCellPoints.first); 
                //greater_assign(meshPoints[meshFace[pointI]], extremeCellPoints.second);
            //}
        //}

        //return AABBox(extremeCellPoints.first, extremeCellPoints.second);
    //}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam::GeometricalTransport

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
