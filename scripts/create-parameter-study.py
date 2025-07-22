#!/usr/bin/env python3

from argparse import ArgumentParser

import testReportCore as trc
import parameterStudyPreparation as psp

app_description="""Prepare the case directories for a given parameter study. Does not perform any
preprocessing steps like meshing or field initialization.
"""

def main():

    #---- Command line arguments ----------------------------------------------
    parser = ArgumentParser(description=app_description)

    parser.add_argument("parameter_file",
                        help="Name of the parameter file for the study")
    parser.add_argument("-p", "--prefix",
                        help="Prefix added to the name of the study folders. The result is prefix-parameter_file_00... etc.",
                        default="",
                        dest="studyprefix")
    parser.add_argument("-t", "--template-case",
                        help="Name of the template case directory. Default: templateCase",
                        default="templateCase",
                        dest="template_case")
    parser.add_argument("-v", "--variants",
                        help="Only use the specified variations. By default, all variations are used. Argument can either be a single number (e.g. 42), a list of numbers (e.g. '3,5,11') or a range ('3 - 10')",
                        default="all",
                        dest="variants")
    #--------------------------------------------------------------------------

    args = parser.parse_args()

    parameter_file = args.parameter_file
    case_name = args.template_case

    # Create vector of all variants to be set up
    variationNumbers = psp.create_variant_vector(args.variants)

    studyprefix = args.studyprefix
    if studyprefix:
        studyprefix = studyprefix + "-"

    # Assemble the parameter variation command according to the given options
    command =   ["pyFoamRunParameterVariation.py",
                 "--allow-derived-changes",
                 "--every-variant-one-case-execution",
                 "--create-database", 
                 "--no-server-process",
                 "--no-execute-solver",
                 "--no-post-templates",
                 "--no-final-templates",
                 "--no-mesh-create",
                 "--no-case-setup",
                 "--parameter-file=default.parameter",
                 "--cloned-case-prefix="+studyprefix+args.parameter_file.rsplit('.', maxsplit=1)[0],
                 case_name,
                 parameter_file
                 ]

    # Finally, create the variations
    psp.setup_variants(command, variationNumbers)


if __name__ == "__main__":
    main()
