#!/usr/bin/env bash

rm -r constant/polyMesh
gmsh -0 constant/levitatorWedgeHex.geo
gmshToFoam  constant/levitatorWedgeHex.msh
rm constant/levitatorWedgeHex.msh
changeDictionary
