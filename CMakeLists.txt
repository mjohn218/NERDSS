# Use modern CMake for nice things like find_package()
cmake_minimum_required(VERSION 3.0)

# Compilation tests using icc was failing on MARCC, this fixes it.
# Evidently used when cross-compiling.
set(CMAKE_TRY_COMPILE_TARGET_TYPE "STATIC_LIBRARY")

# Name of the executable
project(nerdss)

# Use C++11 and only C++11
# Should automatically set a CXX compiler
# If you want to use a different one, provide the path using the -DCMAKE_CXX_COMPILER
# flag, e.g. cmake .. -DCMAKE_CXX_COMPILER=path/to/compiler
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

# Output to bin/ dir
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Flags for debugging:
#   - `-g` -> prints GDB information, for using GDB/LLDB
#   - `-fsanitize=address` -> uses a sanitizer to check for memory issues, such as out-of-range access
#   - `-fno-omit-frame-pointer` -> makes the program store function frame pointers (needed for sanitizer)
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -fsanitize=address -fno-omit-frame-pointer")

# Flags for release
#   - `-O2` -> second level optimization (`-O3` is not recommended, and `-O2` can often be as fast or faster)
set(CMAKE_CXX_FLAGS_RELEASE "-O2")

set(default_build_type "Release")
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to default, '${default_build_type}'.")
    set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE STRING "Choose the type of build." FORCE)
endif()

# Protects against building in the source/root directories.
if(${PROJECT_SOURCE_DIR} STREQUAL ${PROJECT_BINARY_DIR})
    message(FATAL_ERROR "In-source builds not allowed. Please use a different directory (e.g. 'build').")
endif()

# Get system information
cmake_host_system_information(RESULT NUM_CPUS QUERY NUMBER_OF_PHYSICAL_CORES)
cmake_host_system_information(RESULT NUM_THREADS QUERY NUMBER_OF_LOGICAL_CORES)

# Print system information
message(STATUS "We are on a ${CMAKE_SYSTEM_NAME} system")
message(STATUS "The host processor is ${CMAKE_HOST_SYSTEM_PROCESSOR} with ${NUM_CPUS} physical cores and ${NUM_THREADS} logical cores.")

# Set up the executable
set(SHARED_SOURCES
    ###### INCLUDE ######
    ## CLASSES
    include/classes/class_Coord.hpp
    include/classes/class_MDTimer.hpp
    include/classes/class_MolTemplate.hpp
    include/classes/class_Molecule_Complex.hpp
    include/classes/class_Observable.hpp
    include/classes/class_Parameters.hpp
    include/classes/class_Quat.hpp
    include/classes/class_Rxns.hpp
    include/classes/class_SimulVolume.hpp
    include/classes/class_Vector.hpp
    include/classes/class_bngl_parser.hpp

    ## MATH
    include/math/math_functions.hpp
    include/math/matrix.hpp
    include/math/Faddeeva.hpp
    include/math/rand_gsl.hpp
    include/math/constants.hpp

    ## PARSER
    include/parser/parser_functions.hpp

    ## REACTIONS
    include/reactions/association/association.hpp
    include/reactions/bimolecular/2D_reaction_table_functions.hpp
    include/reactions/bimolecular/bimolecular_reactions.hpp
    include/reactions/shared_reaction_functions.hpp
    include/reactions/unimolecular/unimolecular_reactions.hpp

    ## IO
    include/io/io.hpp

    ## SYSTEM SETUP
    include/system_setup/system_setup.hpp

    ## TRAJECTORY FUNCTIONS
    include/trajectory_functions/trajectory_functions.hpp

    ## BOUNDARY CONDITIONS
    include/boundary_conditions/reflect_functions.hpp

    ## ERROR HANDLING
    include/error_handling.hpp

    ###### SRC ######
    ## CLASSES
    src/classes/class_Coord.cpp
    src/classes/class_MDTimer.cpp
    src/classes/class_Molecule_Complex.cpp
    src/classes/class_MolTemplate.cpp
    src/classes/class_Observable.cpp
    src/classes/class_Parameters.cpp
    src/classes/class_Quat.cpp
    src/classes/class_Rxns.cpp
    src/classes/class_SimulVolume.cpp
    src/classes/class_Vector.cpp
    src/classes/class_bngl_parser_functions.cpp

    ## MATH
    src/math/matrix_functions.cpp
    src/math/rand_gsl.cpp
    src/math/Faddeeva.cpp

    ## PARSER
    src/parser/areSameExceptState.cpp
    src/parser/check_for_state_change.cpp
    src/parser/check_for_valid_states.cpp
    src/parser/create_conjugate_reaction_itrs.cpp
    src/parser/determine_bound_iface_index.cpp
    src/parser/determine_iface_indices.cpp
    src/parser/display_all_MolTemplates.cpp
    src/parser/display_all_reactions.cpp
    src/parser/parse_input.cpp
    src/parser/parse_input_array.cpp
    src/parser/parse_molFile.cpp
    src/parser/parse_molecule_bngl.cpp
    src/parser/parse_observable.cpp
    src/parser/parse_reaction.cpp
    src/parser/parse_states.cpp
    src/parser/populate_reaction_lists.cpp
    src/parser/read_bonds.cpp
    src/parser/read_boolean.cpp
    src/parser/read_internal_coordinates.cpp
    src/parser/remove_comment.cpp
    src/parser/write_mol_iface.cpp

    ## REACTIONS
    ### ASSOCIATION
    src/reactions/calculate_phi.cpp
    src/reactions/phi_rotation.cpp
    src/reactions/calculate_omega.cpp
    src/reactions/associate.cpp
    src/reactions/areParallel.cpp
    src/reactions/determine_normal.cpp
    src/reactions/areSameAngle.cpp
    src/reactions/angleSignIsCorrect.cpp
    src/reactions/check_bases.cpp
    src/reactions/conservedMags.cpp
    src/reactions/rotate.cpp
    src/reactions/transform.cpp
    src/reactions/omega_rotation.cpp
    src/reactions/reverse_rotation.cpp
    src/reactions/create_arbitrary_vector.cpp
    src/reactions/theta_rotation.cpp
    src/reactions/determine_rotation_angles.cpp
    src/reactions/requiresSignFlip.cpp
    src/reactions/orient_crds_to_template.cpp
    src/reactions/check_for_structure_overlap.cpp
    src/reactions/conservedRigid.cpp

    ### BIMOLECULAR
    src/reactions/pirr_pfree_ratio_psF.cpp
    src/reactions/perform_bimolecular_state_change.cpp
    src/reactions/determine_3D_bimolecular_reaction_probability.cpp
    src/reactions/determine_2D_bimolecular_reaction_probability.cpp
    src/reactions/survival_function.cpp
    src/reactions/create_pirMatrix.cpp
    src/reactions/get_distance.cpp
    src/reactions/evaluate_binding_same_complex.cpp
    src/reactions/pir_function.cpp
    src/reactions/passocF.cpp
    src/reactions/calc_pirr.cpp
    src/reactions/norm_function.cpp
    src/reactions/create_DDMatrices.cpp
    src/reactions/integrator.cpp
    src/reactions/get_prevNorm.cpp
    src/reactions/create_normMatrix.cpp
    src/reactions/DDpirr_pfree_ratio_ps.cpp
    src/reactions/create_survMatrix.cpp
    src/reactions/size_lookup.cpp
    src/reactions/get_prevSurv.cpp
    src/reactions/check_bimolecular_reactions.cpp
    src/reactions/check_for_structure_overlap.cpp
    src/reactions/check_for_structure_overlap_system.cpp

    ### UNIMOLECULAR
    src/reactions/check_for_unimolecular_reactions.cpp
    src/reactions/check_dissociation.cpp
    src/reactions/create_molecule_and_complex_from_rxn.cpp
    src/reactions/break_interaction.cpp
    src/reactions/determine_parent_complex.cpp

        ### SHARED
    src/reactions/determine_if_reaction_occurs.cpp
    src/reactions/hasIntangibles.cpp
    src/reactions/find_which_reaction.cpp
    src/reactions/find_reaction_rate_state.cpp
    src/reactions/isReactant.cpp
    src/reactions/initialize_molecule_after_reaction.cpp
    src/reactions/find_which_state_change_reaction.cpp
    src/reactions/update_Nboundpairs.cpp

    ### ZEROTH
    src/reactions/check_for_zeroth_order_creation.cpp

    ## SYSTEM SETUP ##
    src/system_setup/areInVicinity.cpp
    src/system_setup/generate_coordinates.cpp
    src/system_setup/initialize_complex.cpp
    src/system_setup/initialize_molecule.cpp
    src/system_setup/set_rMaxLimit.cpp

    ## IO ##
    src/io/eye_candy.cpp
    src/io/init_NboundPairs.cpp
    src/io/init_counterCopyNums.cpp
    src/io/init_print_dimers.cpp
    src/io/print_complex_hist.cpp
    src/io/print_dimers.cpp
    src/io/print_system_information.cpp
    src/io/read_restart.cpp
    src/io/write_NboundPairs.cpp
    src/io/write_all_species.cpp
    src/io/write_complex_components.cpp
    src/io/write_complex_crds.cpp
    src/io/write_crds.cpp
    src/io/write_observables.cpp
    src/io/write_pdb.cpp
    src/io/write_psf.cpp
    src/io/write_restart.cpp
    src/io/write_timestep_information.cpp
    src/io/write_traj.cpp
    src/io/write_xyz.cpp
    src/io/write_xyz_assoc.cpp

    ## BOUNDARY CONDITIONS ##
    src/boundary_conditions/reflect_traj_complex_rad_rot.cpp
    src/boundary_conditions/reflect_complex_rad_rot.cpp
    src/boundary_conditions/reflect_traj_check_span.cpp
    src/boundary_conditions/reflect_traj_complex_rad_rot_nocheck.cpp
    src/boundary_conditions/check_if_spans_box.cpp

    ## TRAJECTORY FUNCTIONS ##
    src/trajectory_functions/create_complex_propagation_vectors.cpp
    src/trajectory_functions/clear_reweight_vecs.cpp
    src/trajectory_functions/sweep_separation_complex_rot_memtest.cpp
    src/trajectory_functions/sweep_separation_complex_rot.cpp
)

# Default executable (serial)
set(DEFAULT_SOURCES EXEs/nerdss.cpp)
add_executable(nerdss ${SHARED_SOURCES} ${DEFAULT_SOURCES})

# Set version
set(PROJECT_VERSION_MAJOR 0)
set(PROJECT_VERSION_MINOR 90)

# Set up external libraries
find_package(GSL REQUIRED)
target_link_libraries(nerdss GSL::gsl GSL::gslcblas)

# Set up header directories
include_directories(include $(GSL_INCLUDE_DIR))

if (${CMAKE_BUILD_TYPE} MATCHES Debug)
    message(STATUS "Compiled with ${CMAKE_CXX_COMPILER} and flags ${CMAKE_CXX_FLAGS_DEBUG}")
else()
    message(STATUS "Compiled with ${CMAKE_CXX_COMPILER} and flags ${CMAKE_CXX_FLAGS_RELEASE}")
endif()
