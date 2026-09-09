# Data dictionary

All tables are UTF-8 CSV. Units and normalization are stated per column below. Empty cells are permitted only where specified in `data-dictionary.json`; they never mean zero.

| Column | Unit | Meaning |
|---|---|---|
| `N` | 1 | Resolution parameter: N is total longitudinal cells for layered cases, or cells per Cartesian cube edge for sphere cases. cellsPerWavelength is cells per acoustic wavelength. |
| `P_relL2` | 1 | Combined real/imaginary pressure L2 error over 1200 uniformly weighted samples including the PML. |
| `Pim_analytic` | Pa | Real or imaginary component of numerical (sim) or analytical complex pressure. |
| `Pim_relL2` | 1 | Real or imaginary pressure component L2 error, normalized by the corresponding analytical component norm. |
| `Pim_sim` | Pa | Real or imaginary component of numerical (sim) or analytical complex pressure. |
| `Pre_analytic` | Pa | Real or imaginary component of numerical (sim) or analytical complex pressure. |
| `Pre_relL2` | 1 | Real or imaginary pressure component L2 error, normalized by the corresponding analytical component norm. |
| `Pre_sim` | Pa | Real or imaginary component of numerical (sim) or analytical complex pressure. |
| `SPL_analytic_dB` | dB re 20 uPa | 10*log10(0.5*abs(P)^2/(20e-6 Pa)^2), using RMS pressure from the complex peak amplitude. |
| `SPL_sim_dB` | dB re 20 uPa | 10*log10(0.5*abs(P)^2/(20e-6 Pa)^2), using RMS pressure from the complex peak amplitude. |
| `absLinf` | 1 | Maximum absolute error of piston pressure magnitude after normalization by rho*c*u0; despite the column name it is dimensionless. |
| `absolute_relative_error` | 1 | Absolute force difference divided by absolute Gorkov force, for nonzero-reference cases. |
| `avgNonOrthogonality` | degree | Average cell-face nonorthogonality angle in the mesh. |
| `avg_nonorthogonality_deg` | degree | Average cell-face nonorthogonality angle in the mesh. |
| `baselineDifference` | 1 | Relative complex-pressure L2 difference from the original serial baseline field for the same physical configuration. |
| `boundaryMode` | not applicable | Category or identifier. Method and ordering meanings are defined in methods.json. |
| `cell_size_m` | m | Target cell size next to the rigid-sphere surface. |
| `cells` | 1 | Number of acoustic control volumes. |
| `cellsPerWavelength` | 1 | Resolution parameter: N is total longitudinal cells for layered cases, or cells per Cartesian cube edge for sphere cases. cellsPerWavelength is cells per acoustic wavelength. |
| `change_to_next_finer` | 1 | Absolute force change to the next finer grid divided by that finer-grid force magnitude. Empty on the finest grid. |
| `difference_from_gorkov` | 1 | Absolute force difference divided by absolute Gorkov force, for nonzero-reference cases. |
| `difference_normalized_by_peak_force` | 1 | Absolute force difference divided by the peak analytical force over position; finite at analytical zero crossings. |
| `directivity_analytic` | 1 | Far-field pressure magnitude divided by its on-axis value for the same curve. |
| `directivity_sim` | 1 | Far-field pressure magnitude divided by its on-axis value for the same curve. |
| `drop_cell_size_m` | m | Target cell size next to the rigid-sphere surface. |
| `farField_onAxisAmplitude_relError` | 1 | Signed numerical/analytical on-axis far-field magnitude ratio minus one. |
| `farField_pressureMagnitude_relL2` | 1 | Relative L2 error of far-field pressure magnitudes over the sampled polar angles. |
| `force_at_1000Pa_N` | N | Archived integrated axial force for the stated incident pressure amplitude. |
| `force_at_1000Pa_div_1e6_N` | N | 1000-Pa result divided by 1000^2 for comparison with the 1-Pa force. |
| `force_at_1Pa_N` | N | Archived integrated axial force for the stated incident pressure amplitude. |
| `gasRelL2` | 1 | Relative complex-pressure L2 error over the indicated cell subset; see archived sphere/advanced.py metrics definition for the interface band. |
| `gorkov_force_y_N` | N | Analytical rigid-Rayleigh-sphere Gorkov axial force for the documented excitation. |
| `hOverLambda` | 1 | Nominal cell size divided by acoustic wavelength. |
| `h_over_a` | 1 | Target cell size near the sphere divided by sphere radius. |
| `h_over_lambda` | 1 | Nominal cell size divided by acoustic wavelength. |
| `host` | not applicable | Category or identifier. Method and ordering meanings are defined in methods.json. |
| `interfaceRelL2` | 1 | Relative complex-pressure L2 error over the indicated cell subset; see archived sphere/advanced.py metrics definition for the interface band. |
| `job_id` | not applicable | Category or identifier. Method and ordering meanings are defined in methods.json. |
| `ka` | 1 | Acoustic wavenumber times sphere radius. |
| `level` | not applicable | Category or identifier. Method and ordering meanings are defined in methods.json. |
| `liquidRelL2` | 1 | Relative complex-pressure L2 error over the indicated cell subset; see archived sphere/advanced.py metrics definition for the interface band. |
| `matrixNonzeros` | 1 | Stored entries in the PETSc real block matrix, including both pressure components. |
| `maxNonOrthogonality` | degree | Maximum cell-face nonorthogonality angle in the mesh. |
| `max_nonorthogonality_deg` | degree | Maximum cell-face nonorthogonality angle in the mesh. |
| `max_relative_residual` | 1 | Maximum true assembled relative linear-system residual over correction sweeps. Empty only for frozen baseline runs that did not report this diagnostic. |
| `max_skewness` | 1 | Maximum mesh skewness diagnostic. |
| `meshFamily` | not applicable | Category or identifier. Method and ordering meanings are defined in methods.json. |
| `mode` | not applicable | Category or identifier. Method and ordering meanings are defined in methods.json. |
| `mumpsFactorMemorySumMB` | MB as reported by MUMPS | Summed factor-memory diagnostic INFOG(22); not total process peak memory. |
| `n_cells` | 1 | Number of acoustic control volumes. |
| `n_samples` | 1 | Number of included sample points, all with valid probe flags. |
| `n_surface_faces` | 1 | Number of acoustic mesh faces on the sphere boundary. |
| `numerical_force_y_N` | N | Integrated numerical axial radiation force in the archived whole-sphere convention. Multiply by 1e15 for plotted fN. |
| `numerical_to_gorkov` | 1 | Signed ratio of numerical to analytical axial force. |
| `ordering` | not applicable | Category or identifier. Method and ordering meanings are defined in methods.json. |
| `p_abs_analytic` | Pa | Magnitude of numerical or analytical complex pressure at the sample point. |
| `p_abs_sim` | Pa | Magnitude of numerical or analytical complex pressure at the sample point. |
| `p_analytic_over_p0` | 1 | Piston on-axis pressure magnitude divided by rho*c*u0. |
| `p_analytic_pa` | Pa | Piston on-axis numerical or analytical pressure magnitude. |
| `p_far_abs_analytic_pa` | Pa | Numerical or analytical far-field pressure magnitude at r_far_m. |
| `p_far_abs_sim_pa` | Pa | Numerical or analytical far-field pressure magnitude at r_far_m. |
| `p_far_imag` | Pa | Real or imaginary component of pressure reconstructed by the Kirchhoff integral. |
| `p_far_real` | Pa | Real or imaginary component of pressure reconstructed by the Kirchhoff integral. |
| `p_sim_over_p0` | 1 | Piston on-axis pressure magnitude divided by rho*c*u0. |
| `p_sim_pa` | Pa | Piston on-axis numerical or analytical pressure magnitude. |
| `parallelDifference` | 1 | Relative complex-pressure L2 difference between the indicated run and the serial field of the same method, in serial cell ordering. |
| `position_m` | m | Sphere-centre coordinate along the standing-wave axis. |
| `position_mm` | mm | Sphere-centre coordinate along the standing-wave axis in millimetres. |
| `pressureOrder` | 1 | Observed order log(E_coarse/E_fine)/log(h_coarse/h_fine). Empty at the first level, where no preceding level exists. |
| `pressureRelL2` | 1 | Relative combined complex-pressure L2 error. Cell-volume weighting is used for homogeneous/sphere studies; see methods.json. |
| `pressureRelL2_cell_centres` | 1 | Combined complex-pressure relative L2 error at cell centres, distinct from the manuscript 1200-point sampled norm. |
| `pressureRelLinf` | 1 | Maximum complex-pressure error magnitude divided by maximum analytical complex-pressure magnitude. |
| `profile_path` | not applicable | Path relative to the dataset root to the exact 1200-point sampled pressure profile. |
| `r_far_m` | m | Radius of the exterior far-field sampling arc. |
| `radius_m` | m | Rigid-sphere radius. |
| `radius_um` | um | Rigid-sphere radius in micrometres. |
| `ranks` | 1 | Number of MPI processes; 1 denotes a serial run. |
| `relL2` | 1 | Piston near-field relative L2 or Linf error of the pressure magnitude, using the uniformly sampled on-axis profile. |
| `relLinf` | 1 | Piston near-field relative L2 or Linf error of the pressure magnitude, using the uniformly sampled on-axis profile. |
| `relative_scaling_deviation` | 1 | Relative deviation from the expected quadratic pressure-amplitude scaling; original archived statistic. |
| `seconds` | s | Comparison-driver elapsed wall time; a single timing measurement, not a scaling result. |
| `segments` | 1 | Number of generating semicircle surface segments. |
| `sigma_max_s_inv` | s^-1 | Positive maximum PML damping magnitude. Some original reports print its negative implementation coefficient. |
| `theta_deg` | degree | Polar angle from the piston axis, ranging from zero to 90 degrees. |
| `time` | solver pseudo-time | Stationary solver output index; not elapsed physical or wall-clock time. |
| `velocityBoundaryRelL2` | 1 | Velocity L2 error over the indicated cell subset; interior cells are more than 2.5 nominal cell widths from the boundary. |
| `velocityInteriorRelL2` | 1 | Velocity L2 error over the indicated cell subset; interior cells are more than 2.5 nominal cell widths from the boundary. |
| `velocityOrder` | 1 | Observed order log(E_coarse/E_fine)/log(h_coarse/h_fine). Empty at the first level, where no preceding level exists. |
| `velocityRelL2` | 1 | Relative L2 error of the complex velocity vector, with the same cell-volume weights as pressure. |
| `x_m` | m | Sample position along the layered domain, including endpoints. |
| `z_m` | m | Piston on-axis sample distance from the baffle (the solver axial direction is y). |
| `z_over_rayleigh` | 1 | Axial distance divided by Rayleigh distance R0=k*a^2/2. |
