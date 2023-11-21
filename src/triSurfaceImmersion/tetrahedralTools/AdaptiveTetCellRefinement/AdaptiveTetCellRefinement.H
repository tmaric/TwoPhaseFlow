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
    Foam::TriSurfaceImmersion::adaptiveTetCellRefinement

Description
    Adaptively refine a set of tetrahedra up to a maximum refinement level.

    The approach decomposes a tetrahedron into 8 smaller tetrahedra using the
    edge centres as newly added vertices. Whether a tetrahedron needs further
    refinement is determined by a refinement criterion (usually evaluating
    if the tetrahedron is intersected by an interface or not). The maximum
    refinement level cann either be manually prescribed or be computed based
    on the average edge length of the initial tetrahedra and a characteristic
    length of the surface.
    The class uses indirect addressing, meaning that there is a unique set of
    points representing the tetrahedra vertices. A tetrahedron is represented
    as a set of four point indices.
    Due to the adaptivity, the class stores tetrahedra for the different
    refinement levels.

SourceFiles
    adaptiveTetCellRefinementI.H

\*---------------------------------------------------------------------------*/

#ifndef adaptiveTetCellRefinement_H
#define adaptiveTetCellRefinement_H

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
                   Class AdaptiveTetCellRefinement Declaration
\*---------------------------------------------------------------------------*/

using indexedTet = std::array<label, 4>;

template<class LevelSet, class RefinementCriterion>
class adaptiveTetCellRefinement
{
private:

    // Typedefs
    using indexTuple = std::tuple<label, label>;
    using edge = indexTuple;

    // Private data
    //- Reference to the level set function representing the surface/interface.
    const LevelSet& levelSet_;

    //- Criterion determining if a tetrahedron needs further refinement or not.
    RefinementCriterion criterion_;

    //- Set of unique points used to represent the tetrahedra
    std::vector<point> points_;

    //- Level set values of the corresponding points / vertices
    std::vector<scalar> levelSetValues_;

    //- Set of tetrahedra represented as sets of point indices
    std::vector<indexedTet> tets_;

    //- Store for each tetrahedron whether it requires refinement
    std::vector<bool> refinement_required_;

    //- Maps tetrahedra edges to edge centre point indices
    std::map<edge, label> edge_to_point_id_;

    //- Maps refinement levels to their corresponding point indices
    //  NOTE: the second entry in the tuple denotes the index of the first
    //  element which is NOT part of the level, aka "[begin, end)"
    std::vector<indexTuple> level_to_pointid_range_;

    //- Maps refinement levels to their corresponding tetrahedra indices
    //  NOTE: the second entry in the tuple denotes the index of the first
    //  element which is NOT part of the level, aka "[begin, end)"
    std::vector<indexTuple> level_to_tetid_range_;

    //- Maximum number of refinement levels used
    label max_refinement_level_;

    //- Flag indicating whether the decomposition has already been performed
    bool decomposition_performed_ = false;

    //- Number of new tetrahedra obtained from decomposing a single tetrahedron
    const label n_tets_from_decomposition = 8;


    // Private Member Functions
    //- Compute edge lengths for a tetrahedron
    std::array<scalar, 6> edge_lengths(const indexedTet& tet) const;

    //- Compute maximum refinement level if it is not prescribed
    label compute_max_refinement_level();

    //- Invoke the adaptive decomposition of the initial tetrahedra
    void compute_decomposition();

    //- For a given refinement level, determine which tetrahedra need refinement
    label flag_tets_for_refinement(int level);

    //- Increase size of the container storing tetrahedra
    void update_tet_container_sizes(int level, int n_new_tets);

    //- Increase size of edge to point map
    void update_edge_to_point_map(int level);

    //- Add edges of a tetrahedron to the edge-to-point-id map
    void add_to_map(std::array<edge, 6> tet_edges);

    //- Return the edges of a given tetrahedron
    std::array<edge, 6> edges(const indexedTet& tet) const;

    //- Create new tetrahedra through refinement for the given level
    void create_refined_tets(int level);

    //- Decompose the given tetrahedron and add it to the set of tertahedra
    void decompose_and_add_new_tets(const indexedTet& tet, label start_id);

    //- For points generated by refinement compute their level set values
    void computeLevelSetValues(int level);

    //- Save the resulting decomposition as a VTL legacy file
    void save_decomposition_as_vtk(const std::vector<indexedTet>& tets,
        const std::vector<point>& points,
        const std::vector<scalar>& signed_distance,
        const std::vector<label>& refinement_levels,
        std::string file_name) const;


public:

    // Constructors
    adaptiveTetCellRefinement(const LevelSet& levelSet,
        std::vector<point> points,
        std::vector<scalar> levelSetValues,
        std::vector<indexedTet> tets,
        label max_refine_level = -1);


    // Member Functions
    //- Reference to the point set resulting from decompostion
    const std::vector<point>& points();

    //- Reference to the level set values 
    const std::vector<scalar>& signedDistance();

    //- Return the final tetrahedra from decomposition.
    //  The resulting set of tetrahedra is such that there is no volumetric
    //  overlap, meaning that the accumulated volume of the initial tetrahedra
    //  is equal to the accumulated volume of the tetrahedra after
    //  decomposition.
    std::vector<indexedTet> resultingTets();

    //- Return the maximum refinement level
    label refinementLevel() const;

    //- For each final tetrahedron gives its refinement level 
    std::vector<label> refinementLevels(label n_tets);

    //- Save the final decomposition as VTK for cell cellID
    void writeTets(label cellID);

    //- Print infos about number of tetrahedra and points for each level
    void printLevelInfos() const;

    //- Print number of tetrahedra and each tetrahedron's indices
    void printTets() const;

    //- Print points and their indices
    void printPoints() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "AdaptiveTetCellRefinementI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam::TriSurfaceImmersion

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
