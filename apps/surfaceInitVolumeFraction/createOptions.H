argList::addNote
(
    "Initialize a volume fraction field from a triangulated surface or level set using the SMCI/A algorithm"
);

argList::addOption
(
    "surfaceType",
    "triSurface / levelSet",
    "Surface type, meaning either a triSurface, e.g. from an STL, or a level set. \nDefault: triSurface"
);

argList::addOption
(
    "surfaceFile",
    "Surface mesh file",
    "File containing the interface description. Requires a consistent, inward normal orientation. Use OpenFOAM's 'surfaceOrient' tool for this purpose.\nDefault: surface.stl"
);

argList::addOption
(
    "narrowBandWidth", 
    "scalar",
    "Number of cells in interface normal direction constituting the narrow band. \nDefault: 4.0"
); 

argList::addOption
(
    "fieldName", 
    "alpha.water",
    "Name of the volume fraction field.\nDefault: alpha.water" 
); 

argList::addOption
(
    "surfaceFieldName", 
    "alphaSurface.water",
    "Name of the area fraction field.\nDefault: none" 
); 

argList::addOption
(
    "algorithm", 
    "SMCI/SMCA",
    "Name of the volume fraction calculator.\nDefault: SMCI" 
); 

argList::addOption
(
    "refinementLevel", 
    "label",
    "Maximum refinement level to be used (SMCA only).\nDefault: -1 (auto mode)"
); 

argList::addOption
(
    "relError", 
    "scalar",
    "Only for closed, implicit surfaces: ensure that the global relative volume error is below the given threshold. Note: using this option overrides the '-refinementLevel' option.\nDefault: -1.0 (Disabled)"
); 

argList::addBoolOption
(
    "writeGeometry",
    "Write the intersected geometry (SMCI) or the testrahedral decomposition (SMCA)used to compute the volume fractions"
);

argList::addBoolOption
(
    "invert",
    "Invert the computed volume fraction field with respect to the surface orientation."
);

argList::addBoolOption
(
    "writeAllFields",
    "Write out all fields used by the initialization method."
);

argList::addBoolOption
(
    "checkVolume",
    "Check the volume given by volume fraction compared with the volume of the surface mesh. Works only with closed surfaces."
);
