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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
template<class LevelSet, class RefinementCriterion>
std::array<scalar, 3> adaptiveTriFaceRefinement<LevelSet,
RefinementCriterion>::edge_lengths
(
    const indexedTri& tri
) const
{
    std::array<scalar, 3> lengths{};
    const auto tri_edges = edges(tri);

    for (int idx = 0; idx != 3; ++idx)
    {
        const auto [p_id, q_id] = tri_edges[idx];
        lengths[idx] = mag(points_[p_id] - points_[q_id]);
    }

    return lengths;
}


template<class LevelSet, class RefinementCriterion>
label adaptiveTriFaceRefinement<LevelSet,
RefinementCriterion>::compute_max_refinement_level()
{
    // TODO (TT): use average edge length of tris for now and see, how
    // it works.
    scalar avg_length = 0.0;

    for (const auto& tri : tris_)
    {
        const auto lengths = edge_lengths(tri);
        avg_length += std::accumulate(lengths.begin(), lengths.end(), 0.0);
    }

    avg_length /= 3 * tris_.size();
    auto ref_length = levelSet_.referenceLength();
    label level = 0;

    while (avg_length > ref_length)
    {
        avg_length /= 2.0;
        ++level;
    }

    return level;
}


template<class LevelSet, class RefinementCriterion>
void adaptiveTriFaceRefinement<LevelSet,
RefinementCriterion>::compute_decomposition()
{
    if (decomposition_performed_)
    {
        return;
    }

    for (int level = 0; level != max_refinement_level_; ++level)
    {
        const auto n_refined_tris = flag_tris_for_refinement(level);
        update_tri_container_sizes(
            level, n_refined_tris * n_tris_from_decomposition);
        update_edge_to_point_map(level);
        create_refined_tris(level);
        computeLevelSetValues(level);
    }

    decomposition_performed_ = true;
}


template<class LevelSet, class RefinementCriterion>
label adaptiveTriFaceRefinement<LevelSet,
RefinementCriterion>::flag_tris_for_refinement(const int level)
{
    const auto [start, end] = level_to_triid_range_[level];

    auto n_tris_to_refine = 0;

    for (auto idx = start; idx != end; ++idx)
    {
        refinement_required_[idx] =
            considerIntersected(tris_[idx], points_, levelSetValues_, criterion_);

        if (refinement_required_[idx])
        {
            ++n_tris_to_refine;
        }
    }

    return n_tris_to_refine;
}


template<class LevelSet, class RefinementCriterion>
void adaptiveTriFaceRefinement<LevelSet,
RefinementCriterion>::update_tri_container_sizes
(
     const int level,
    const int n_new_tris
)
{
    // Update vector sizes related to the number of tris
    level_to_triid_range_[level + 1] =
        indexTuple{tris_.size(), tris_.size() + n_new_tris};
    tris_.resize(tris_.size() + n_new_tris);
    refinement_required_.resize(refinement_required_.size() + n_new_tris);

    // NOTE: only tris whose refine flag is false are returned by the public
    // member functions. Thus, default initialize the newly added fields to
    // false.
    const auto [next_start, next_end] = level_to_triid_range_[level + 1];
    for (auto idx = next_start; idx != next_end; ++idx)
    {
        refinement_required_[idx] = false;
    }
}


template<class LevelSet, class RefinementCriterion>
void adaptiveTriFaceRefinement<LevelSet, RefinementCriterion>::add_to_map
(
    const std::array<edge, 3> tri_edges
)
{
    for (const auto tri_edge : tri_edges)
    {
        edge_to_point_id_[tri_edge] = 0;
    }
}


template<class LevelSet, class RefinementCriterion>
void adaptiveTriFaceRefinement<LevelSet,
RefinementCriterion>::update_edge_to_point_map(const int level)
{
    edge_to_point_id_.clear();

    const auto [start, end] = level_to_triid_range_[level];

    for (auto idx = start; idx != end; ++idx)
    {
        if (refinement_required_[idx])
        {
            add_to_map(edges(tris_[idx]));
        }
    }

    // Update container sizes
    const auto n_new_points = edge_to_point_id_.size();
    level_to_pointid_range_[level + 1] =
        indexTuple{points_.size(), points_.size() + n_new_points};
    points_.resize(points_.size() + n_new_points);
    levelSetValues_.resize(levelSetValues_.size() + n_new_points);

    // Add global point ids to the mapping and compute the new points
    label global_point_id = std::get<0>(level_to_pointid_range_[level + 1]);

    for (auto& edge_to_point : edge_to_point_id_)
    {
        edge_to_point.second = global_point_id;
        auto [p1_id, p2_id] = edge_to_point.first;
        points_[global_point_id] = 0.5 * (points_[p1_id] + points_[p2_id]);

        ++global_point_id;
    }
}


template<class LevelSet, class RefinementCriterion>
std::array<std::tuple<label, label>, 3> adaptiveTriFaceRefinement<LevelSet,
RefinementCriterion>::edges(const indexedTri& tri) const
{
    return std::array<edge, 3>
           {
            edge{tri[0], tri[1]},
            edge{tri[0], tri[2]},
            edge{tri[1], tri[2]},
           };
}


template<class LevelSet, class RefinementCriterion>
void adaptiveTriFaceRefinement<LevelSet,
RefinementCriterion>::create_refined_tris(int level)
{
    const auto [start, end] = level_to_triid_range_[level];
    auto refined_tri_id = end;

    for (auto idx = start; idx != end; ++idx)
    {
        if (refinement_required_[idx])
        {
            decompose_and_add_new_tris(tris_[idx], refined_tri_id);
            refined_tri_id += n_tris_from_decomposition;
        }
    }
}


template<class LevelSet, class RefinementCriterion>
void adaptiveTriFaceRefinement<LevelSet,
RefinementCriterion>::decompose_and_add_new_tris
(
    const indexedTri& tri,
    const label start_id
)
{
    auto tri_edges = edges(tri);

    // Translate edges into point ids
    std::array<label, 3> pids{};

    for (int idx = 0; idx != 3; ++idx)
    {
        pids[idx] = edge_to_point_id_[tri_edges[idx]];
    }

    // Define new tris
    // Triangles containing each one of the original triangle vertices
    tris_[start_id] = indexedTri{tri[0], pids[0], pids[1]};
    tris_[start_id + 1] = indexedTri{tri[1], pids[0], pids[2]};
    tris_[start_id + 2] = indexedTri{tri[2], pids[1], pids[2]};

    // Triangles consituted by three new points on the edge midpoints
    tris_[start_id + 3] = indexedTri{pids[0], pids[1], pids[2]};
}


template<class LevelSet, class RefinementCriterion>
void adaptiveTriFaceRefinement<LevelSet,
RefinementCriterion>::computeLevelSetValues(const int level)
{
    const auto [start, end] = level_to_pointid_range_[level + 1];

    for (auto idx = start; idx != end; ++idx)
    {
        levelSetValues_[idx] =
            levelSetValue(levelSet_, points_[idx], criterion_);
    }
}

template<class LevelSet, class RefinementCriterion>
void adaptiveTriFaceRefinement<LevelSet, RefinementCriterion>::
save_decomposition_as_vtk
(
    const std::vector<indexedTri>& tris,
    const std::vector<point>& points,
    const std::vector<scalar>& signed_distance,
    const std::vector<label>& refinement_levels,
    const std::vector<scalar>& area_fractions,
    std::string file_name
) const
{
    // Use VTK legacy format for now for the sake of simplicity (TT)
    std::ofstream out_file;
    // Ensure directory exists
    std::filesystem::create_directory("VTK");
    out_file.open("VTK/" + file_name);

    // Header
    out_file << "# vtk DataFile Version 3.0\n";

    out_file << "Trianglular decomposition of a face.\n";
    out_file << "ASCII\n";
    out_file << "DATASET UNSTRUCTURED_GRID\n";

    // Write points
    out_file << "POINTS " << std::to_string(points.size()) << " double\n";
    for (const auto& p : points)
    {
        out_file << std::to_string(p[0]) << " " << std::to_string(p[1]) << " "
                 << std::to_string(p[2]) << "\n";
    }

    // Write tris
    out_file << "CELLS " << std::to_string(tris.size()) << " "
             << std::to_string(4 * tris.size()) << "\n";
    for (const auto& tri : tris)
    {
        out_file << "3 " << std::to_string(tri[0]) << " "
                         << std::to_string(tri[1]) << " "
                         << std::to_string(tri[2])
                         << "\n";
    }
    out_file << "CELL_TYPES " << std::to_string(tris.size()) << "\n";
    for (uint idx = 0; idx != tris.size(); ++idx)
    {
        out_file << "5\n";
    }

    // Write signed distance
    out_file << "POINT_DATA " << std::to_string(points.size()) << "\n";
    out_file << "Scalars signed_distance double 1\n";
    out_file << "LOOKUP_TABLE default\n";
    for (const auto phi : signed_distance)
    {
        out_file << std::to_string(phi) << "\n";
    }

    // Write refinement level as cell data
    out_file << "CELL_DATA " << std::to_string(tris.size()) << "\n";
    out_file << "Scalars refinement_level double 1\n";
    out_file << "LOOKUP_TABLE default\n";
    for (const label alpha : refinement_levels)
    {
        out_file << std::to_string(alpha) << "\n";
    }

    // Write elementary area fractions as cell data
    //out_file << "CELL_DATA " << std::to_string(tris.size()) << "\n";
    out_file << "Scalars area_fractions double 1\n";
    out_file << "LOOKUP_TABLE default\n";
    for (const label af : area_fractions)
    {
        out_file << std::to_string(af) << "\n";
    }
    
    out_file.close();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
template<class LevelSet, class RefinementCriterion>
adaptiveTriFaceRefinement<LevelSet,
    RefinementCriterion>::adaptiveTriFaceRefinement(const LevelSet& levelSet,
    const std::vector<point> points,
    const std::vector<scalar> levelSetValues,
    const std::vector<indexedTri> tris,
    const label max_refine_level)
    : levelSet_{levelSet}, criterion_{}, points_{points},
      levelSetValues_{levelSetValues}, tris_{tris},
      refinement_required_(tris_.size(), false), edge_to_point_id_{},
      level_to_pointid_range_{}, level_to_triid_range_{},
      max_refinement_level_{max_refine_level}
{
    if (max_refinement_level_ < 0)
    {
        max_refinement_level_ = compute_max_refinement_level();
    }

    level_to_pointid_range_.resize(max_refinement_level_ + 1);
    level_to_pointid_range_[0] = indexTuple{0, points_.size()};

    level_to_triid_range_.resize(max_refinement_level_ + 1);
    level_to_triid_range_[0] = indexTuple{0, tris_.size()};
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
template<class LevelSet, class RefinementCriterion>
const std::vector<point>& adaptiveTriFaceRefinement<LevelSet,
RefinementCriterion>::points()
{
    compute_decomposition();

    return points_;
}


template<class LevelSet, class RefinementCriterion>
const std::vector<scalar>& adaptiveTriFaceRefinement<LevelSet,
RefinementCriterion>::signedDistance()
{
    compute_decomposition();

    return levelSetValues_;
}


template<class LevelSet, class RefinementCriterion>
std::vector<indexedTri> adaptiveTriFaceRefinement<LevelSet,
RefinementCriterion>::resultingTris()
{
    std::vector<indexedTri> final_tris{};

    compute_decomposition();

    const auto n_tris = std::count(
        refinement_required_.begin(), refinement_required_.end(), false);

    final_tris.resize(n_tris);

    int count = 0;
    for (int idx(refinement_required_.size() - 1); idx != -1; --idx)
    {
        if (refinement_required_[idx] == false)
        {
            final_tris[count] = tris_[idx];
            ++count;
        }
    }

    return final_tris;
}


template<class LevelSet, class RefinementCriterion>
label adaptiveTriFaceRefinement<LevelSet,
RefinementCriterion>::refinementLevel() const
{
    return max_refinement_level_;
}


template<class LevelSet, class RefinementCriterion>
void adaptiveTriFaceRefinement<LevelSet, RefinementCriterion>::printLevelInfos()
    const
{
    Info << "Number of tris and points per level\n";

    for (auto level = 0; level <= max_refinement_level_; ++level)
    {
        auto [tstart, tend] = level_to_triid_range_[level];
        auto [pstart, pend] = level_to_pointid_range_[level];
        Info << "Level = " << level << ":\n"
             << "\tn_tris = " << (tend - tstart) << "\n"
             << "\tn_points = " << (pend - pstart) << "\n";
    }
}


template<class LevelSet, class RefinementCriterion>
void adaptiveTriFaceRefinement<LevelSet, RefinementCriterion>::printTris() const
{
    Info << "Indexed tris (n = " << tris_.size() << " in total):\n";

    for (uint idx = 0; idx != tris_.size(); ++idx)
    {
        const auto& tri = tris_[idx];
        Info << "tri_id = " << idx << ": " << tri[0] << "\t" << tri[1] << "\t"
             << tri[2] << "\n";
    }
}


template<class LevelSet, class RefinementCriterion>
void adaptiveTriFaceRefinement<LevelSet, RefinementCriterion>::printPoints() const
{
    Info << "Points:\n";

    for (uint idx = 0; idx != points_.size(); ++idx)
    {
        const auto& p = points_[idx];

        Info << "p_id = " << idx << ": " << p << "\n";
    }
}


template<class LevelSet, class RefinementCriterion>
std::vector<label> adaptiveTriFaceRefinement<LevelSet,
RefinementCriterion>::refinementLevels(const label n_tris)
{
    std::vector<label> tri_ref_levels(n_tris);

    label tri_id = 0;

    compute_decomposition();

    for (int level = max_refinement_level_; level >= 0; --level)
    {
        const auto [idx_end, idx_start] = level_to_triid_range_[level];

        for (int idx = idx_start - 1; idx >= idx_end; --idx)
        {
            if (refinement_required_[idx] == false)
            {
                tri_ref_levels[tri_id] = level;
                ++tri_id;
            }
        }
    }

    return tri_ref_levels;
}


template<class LevelSet, class RefinementCriterion>
void adaptiveTriFaceRefinement<LevelSet, RefinementCriterion>::writeTris
(
    const label faceID,
    const std::vector<scalar>& areaFractions
)
{
    compute_decomposition();

    auto final_tris = resultingTris();

    save_decomposition_as_vtk
    (
        final_tris,
        points_,
        levelSetValues_,
        refinementLevels(final_tris.size()),
        areaFractions,
        "face-" + std::to_string(faceID) + "-decomposition.vtk"
    );
}
// ************************************************************************* //
