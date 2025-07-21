#!/usr/bin/env bash
cp -r 0.orig 0
blockMesh
# cartesianMesh
topoSet
createPatch -overwrite
setExprBoundaryFields
decomposePar
