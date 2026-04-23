# cfMeshTest

This case reproduces the Andrade 2019 axisymmetric wedge levitator geometry
from `WedgeLevitatorDropAlpha`, but generates the mesh with `cfMesh`.

Workflow:

1. Render a 2D levitator cross-section and mesh dictionaries from `caseParams.sh`.
2. Generate the axisymmetric cross-section with `cartesian2DMesh`.
3. Extrude the 2D mesh into a 1-degree wedge with `extrudeMesh`.
4. Apply interface-only refinement with a semicircular chain of `line` refinement objects in the 2D meshing step.
5. Initialise `alpha.water` with a spherical drop at the levitator mid-gap.

Run:

The resulting 3D mesh keeps explicit `transducer1`, `transducer2`,
`reflector1`, `reflector2`, `openAir`, `axis`, and wedge patches.
