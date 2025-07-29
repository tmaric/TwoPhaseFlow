/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2011 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2019 AUTHOR,AFFILIATION
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::TriSurfaceImmersion::adaptiveTriFaceRefinement

Description
    Adaptively refine a set of triangles up to a maximum refinement level.

    The approach decomposes a triangle into 4 smaller triangles using the
    edge centres as newly added vertices. Whether a triangle needs further
    refinement is determined by a refinement criterion (usually evaluating
    if the trinagle is intersected by an interface or not). The maximum
    refinement level cann either be manually prescribed or be computed based
    on the average edge length of the initial triangles and a characteristic
    length of the surface.
    The class uses indirect addressing, meaning that there is a unique set of
    points representing the triangle vertices. A triangle is represented
    as a set of three point indices.
    Due to the adaptivity, the class stores triangles for the different
    refinement levels.

SourceFiles
    adaptiveTriFaceRefinementI.H

\*---------------------------------------------------------------------------*/

#ifndef adaptiveTriFaceRefinement_H
#define adaptiveTriFaceRefinement_H

#include <algorithm>
#include <array>
#include <filesystem>
#include <map>
#include <utility>
#include <vector>

#include "fvCFD.H"
#include "tmp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::TriSurfaceImmersion
{

/*---------------------------------------------------------------------------*\
                   Class AdaptiveTriFaceRefinement Declaration
\*---------------------------------------------------------------------------*/

using indexedTri = std::array<label, 3>;

template<class LevelSet, class RefinementCriterion>
class adaptiveTriFaceRefinement
{
private:

    // Typedefs
    using indexTuple = std::tuple<label, label>;
    using edge = indexTuple;

    // Private data
    //- Reference to the level set function representing the surface/interface.
    const LevelSet& levelSet_;

    //- Criterion determining if a triagnle needs further refinement or not.
    RefinementCriterion criterion_;

    //- Set of unique points used to represent the triangles
    std::vector<point> points_;

    //- Level set values of the corresponding points / vertices
    std::vector<scalar> levelSetValues_;

    //- Set of triangles represented as sets of point indices
    std::vector<indexedTri> tris_;

    //- Store for each triangle whether it requires refinement
    std::vector<bool> refinement_required_;

    //- Maps triangle edges to edge centre point indices
    std::map<edge, label> edge_to_point_id_;

    //- Maps refinement levels to their corresponding point indices
    //  NOTE: the second entry in the tuple denotes the index of the first
    //  element which is NOT part of the level, aka "[begin, end)"
    std::vector<indexTuple> level_to_pointid_range_;

    //- Maps refinement levels to their corresponding triangle indices
    //  NOTE: the second entry in the tuple denotes the index of the first
    //  element which is NOT part of the level, aka "[begin, end)"
    std::vector<indexTuple> level_to_triid_range_;

    //- Maximum number of refinement levels used
    label max_refinement_level_;

    //- Flag indicating whether the decomposition has already been performed
    bool decomposition_performed_ = false;

    //- Number of new triangles obtained from decomposing a single triangle
    const label n_tris_from_decomposition = 4;


    // Private Member Functions
    //- Compute edge lengths for a triangle
    std::array<scalar, 3> edge_lengths(const indexedTri& tet) const;

    //- Compute maximum refinement level if it is not prescribed
    label compute_max_refinement_level();

    //- Invoke the adaptive decomposition of the initial triangles
    void compute_decomposition();

    //- For a given refinement level, determine which triangles need refinement
    label flag_tris_for_refinement(int level);

    //- Increase size of the container storing triangles
    void update_tri_container_sizes(int level, int n_new_tris);

    //- Increase size of edge to point map
    void update_edge_to_point_map(int level);

    //- Add edges of a triangle to the edge-to-point-id map
    void add_to_map(std::array<edge, 3> tri_edges);

    //- Return the edges of a given triangle
    std::array<edge, 3> edges(const indexedTri& tet) const;

    //- Create new triangles through refinement for the given level
    void create_refined_tris(int level);

    //- Decompose the given itriangle  and add it to the set of triangles
    void decompose_and_add_new_tris(const indexedTri& tri, label start_id);

    //- Three points generated by refinement compute their level set values
    void computeLevelSetValues(int level);

    //- Save the resulting decomposition as a VTK legacy file
    void save_decomposition_as_vtk(const std::vector<indexedTri>& tris,
        const std::vector<point>& points,
        const std::vector<scalar>& signed_distance,
        const std::vector<label>& refinement_levels,
        const std::vector<scalar>& area_fractions,
        std::string file_name) const;


public:

    // Constructors
    adaptiveTriFaceRefinement(const LevelSet& levelSet,
        std::vector<point> points,
        std::vector<scalar> levelSetValues,
        std::vector<indexedTri> tris,
        label max_refine_level = -1);


    // Member Functions
    //- Reference to the point set resulting from decompostion
    const std::vector<point>& points();

    //- Reference to the level set values 
    const std::vector<scalar>& signedDistance();

    //- Return the final triangles from decomposition.
    //  The resulting set of triangles is such that there is no volumetric
    //  overlap, meaning that the accumulated area of the initial triangles
    //  is equal to the accumulated area of the triangles after
    //  decomposition.
    std::vector<indexedTri> resultingTris();

    //- Return the maximum refinement level
    label refinementLevel() const;

    //- For each final triangle gives its refinement level 
    std::vector<label> refinementLevels(label n_tris);

    //- Save the final decomposition as VTK for face faceID
    void writeTris(label faceID, const std::vector<scalar>& areaFractions);

    //- Print infos about number of triangles and points for each level
    void printLevelInfos() const;

    //- Print number of triangles and each triangle's indices
    void printTris() const;

    //- Print points and their indices
    void printPoints() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "AdaptiveTriFaceRefinementI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam::TriSurfaceImmersion

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
