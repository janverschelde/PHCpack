void initialize ( void );

/* wrapping functions in phcpack.h starts from here */

static PyObject *py2c_PHCpack_version_string
 ( PyObject *self, PyObject *args );
static PyObject *py2c_set_seed ( PyObject *self, PyObject *args );
static PyObject *py2c_read_target_system ( PyObject *self, PyObject *args );
static PyObject *py2c_read_start_system ( PyObject *self, PyObject *args );
static PyObject *py2c_define_output_file ( PyObject *self, PyObject *args );
static PyObject *py2c_write_target_system ( PyObject *self, PyObject *args );
static PyObject *py2c_write_start_system ( PyObject *self, PyObject *args );
/* moving systems from and to containers */
static PyObject *py2c_copy_target_system_to_container
 ( PyObject *self, PyObject *args );
static PyObject *py2c_copy_dobldobl_target_system_to_container
 ( PyObject *self, PyObject *args );
static PyObject *py2c_copy_quaddobl_target_system_to_container
 ( PyObject *self, PyObject *args );
static PyObject *py2c_copy_multprec_target_system_to_container
 ( PyObject *self, PyObject *args );
static PyObject *py2c_copy_container_to_target_system
 ( PyObject *self, PyObject *args );
static PyObject *py2c_copy_dobldobl_container_to_target_system
 ( PyObject *self, PyObject *args );
static PyObject *py2c_copy_quaddobl_container_to_target_system
 ( PyObject *self, PyObject *args );
static PyObject *py2c_copy_multprec_container_to_target_system
 ( PyObject *self, PyObject *args );
static PyObject *py2c_copy_start_system_to_container
 ( PyObject *self, PyObject *args );
static PyObject *py2c_copy_dobldobl_start_system_to_container
 ( PyObject *self, PyObject *args );
static PyObject *py2c_copy_quaddobl_start_system_to_container
 ( PyObject *self, PyObject *args );
static PyObject *py2c_copy_multprec_start_system_to_container
 ( PyObject *self, PyObject *args );
static PyObject *py2c_copy_container_to_start_system 
 ( PyObject *self, PyObject *args );
static PyObject *py2c_copy_dobldobl_container_to_start_system 
 ( PyObject *self, PyObject *args );
static PyObject *py2c_copy_quaddobl_container_to_start_system 
 ( PyObject *self, PyObject *args );
static PyObject *py2c_copy_multprec_container_to_start_system 
 ( PyObject *self, PyObject *args );
/* creation of homotopy and tracking all paths */
static PyObject *py2c_create_homotopy ( PyObject *self, PyObject *args );
static PyObject *py2c_create_homotopy_with_gamma
 ( PyObject *self, PyObject *args );
static PyObject *py2c_create_dobldobl_homotopy
 ( PyObject *self, PyObject *args );
static PyObject *py2c_create_dobldobl_homotopy_with_gamma
 ( PyObject *self, PyObject *args );
static PyObject *py2c_create_quaddobl_homotopy
 ( PyObject *self, PyObject *args );
static PyObject *py2c_create_quaddobl_homotopy_with_gamma
 ( PyObject *self, PyObject *args );
static PyObject *py2c_create_multprec_homotopy
 ( PyObject *self, PyObject *args );
static PyObject *py2c_create_multprec_homotopy_with_gamma
 ( PyObject *self, PyObject *args );
static PyObject *py2c_clear_homotopy ( PyObject *self, PyObject *args );
static PyObject *py2c_clear_dobldobl_homotopy
 ( PyObject *self, PyObject *args );
static PyObject *py2c_clear_quaddobl_homotopy
 ( PyObject *self, PyObject *args );
static PyObject *py2c_clear_multprec_homotopy
 ( PyObject *self, PyObject *args );
static PyObject *py2c_write_start_solutions ( PyObject *self, PyObject *args );
static PyObject *py2c_tune_continuation_parameters
 ( PyObject *self, PyObject *args );
static PyObject *py2c_show_continuation_parameters
 ( PyObject *self, PyObject *args );
static PyObject *py2c_autotune_continuation_parameters
 ( PyObject *self, PyObject *args );
static PyObject *py2c_determine_output_during_continuation
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solve_by_homotopy_continuation
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solve_by_dobldobl_homotopy_continuation
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solve_by_quaddobl_homotopy_continuation
 ( PyObject *self, PyObject *args );
/* moving solutions from and to containers */
static PyObject *py2c_copy_target_solutions_to_container
 ( PyObject *self, PyObject *args );
static PyObject *py2c_copy_dobldobl_target_solutions_to_container
 ( PyObject *self, PyObject *args );
static PyObject *py2c_copy_quaddobl_target_solutions_to_container
 ( PyObject *self, PyObject *args );
static PyObject *py2c_copy_multprec_target_solutions_to_container
 ( PyObject *self, PyObject *args );
static PyObject *py2c_copy_container_to_target_solutions
 ( PyObject *self, PyObject *args );
static PyObject *py2c_copy_dobldobl_container_to_target_solutions
 ( PyObject *self, PyObject *args );
static PyObject *py2c_copy_quaddobl_container_to_target_solutions
 ( PyObject *self, PyObject *args );
static PyObject *py2c_copy_multprec_container_to_target_solutions
 ( PyObject *self, PyObject *args );
static PyObject *py2c_copy_start_solutions_to_container
 ( PyObject *self, PyObject *args );
static PyObject *py2c_copy_dobldobl_start_solutions_to_container
 ( PyObject *self, PyObject *args );
static PyObject *py2c_copy_quaddobl_start_solutions_to_container
 ( PyObject *self, PyObject *args );
static PyObject *py2c_copy_multprec_start_solutions_to_container
 ( PyObject *self, PyObject *args );
static PyObject *py2c_copy_container_to_start_solutions
 ( PyObject *self, PyObject *args );
static PyObject *py2c_copy_dobldobl_container_to_start_solutions
 ( PyObject *self, PyObject *args );
static PyObject *py2c_copy_quaddobl_container_to_start_solutions
 ( PyObject *self, PyObject *args );
static PyObject *py2c_copy_multprec_container_to_start_solutions
 ( PyObject *self, PyObject *args );
/* black box solver, mixed volume calculator, and Newton step */
static PyObject *py2c_solve_system ( PyObject *self, PyObject *args );
static PyObject *py2c_solve_Laurent_system ( PyObject *self, PyObject *args );
static PyObject *py2c_mixed_volume ( PyObject *self, PyObject *args );
static PyObject *py2c_deflate ( PyObject *self, PyObject *args );
static PyObject *py2c_Newton_step ( PyObject *self, PyObject *args );
static PyObject *py2c_dobldobl_Newton_step ( PyObject *self, PyObject *args );
static PyObject *py2c_quaddobl_Newton_step ( PyObject *self, PyObject *args );
static PyObject *py2c_multprec_Newton_step ( PyObject *self, PyObject *args );

/* wrapping functions in syscon.h starts from here */

static PyObject *py2c_syscon_read_system ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_read_Laurent_system
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_read_dobldobl_system
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_read_quaddobl_system
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_read_multprec_system
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_write_system ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_write_Laurent_system
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_write_dobldobl_system
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_write_quaddobl_system
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_write_multprec_system
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_clear_system ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_clear_Laurent_system
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_clear_dobldobl_system
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_clear_quaddobl_system
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_clear_multprec_system
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_number_of_symbols
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_write_symbols
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_string_of_symbols
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_remove_symbol_name
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_clear_symbol_table
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_number_of_polynomials
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_number_of_dobldobl_polynomials
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_number_of_quaddobl_polynomials
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_number_of_multprec_polynomials
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_number_of_Laurentials
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_initialize_number
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_initialize_number_of_dobldobl_polynomials
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_initialize_number_of_quaddobl_polynomials
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_initialize_number_of_multprec_polynomials
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_initialize_number_of_Laurentials
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_degree_of_polynomial
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_degree_of_dobldobl_polynomial
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_degree_of_quaddobl_polynomial
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_degree_of_multprec_polynomial
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_number_of_terms ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_number_of_Laurent_terms
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_retrieve_term ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_store_polynomial
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_store_dobldobl_polynomial
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_store_quaddobl_polynomial
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_store_multprec_polynomial
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_load_polynomial
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_load_dobldobl_polynomial
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_load_quaddobl_polynomial
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_load_multprec_polynomial
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_store_Laurential
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_total_degree ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_standard_drop_variable_by_index
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syson_standard_drop_variable_by_name
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_dobldobl_drop_variable_by_index
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_dobldobl_drop_variable_by_name
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_quaddobl_drop_variable_by_index
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_quaddobl_drop_variable_by_name
 ( PyObject *self, PyObject *args );

/* wrapping functions in solcon.h starts from here */

static PyObject *py2c_solcon_read_solutions ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_read_dobldobl_solutions
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_read_quaddobl_solutions
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_read_multprec_solutions
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_write_solutions ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_write_dobldobl_solutions
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_write_quaddobl_solutions
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_write_multprec_solutions
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_clear_solutions ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_clear_dobldobl_solutions
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_clear_quaddobl_solutions
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_clear_multprec_solutions
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_open_solution_input_file
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_length_solution_string
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_length_dobldobl_solution_string
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_length_quaddobl_solution_string
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_length_multprec_solution_string
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_write_solution_string
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_write_dobldobl_solution_string
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_write_quaddobl_solution_string
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_write_multprec_solution_string
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_append_solution_string
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_append_dobldobl_solution_string
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_append_quaddobl_solution_string
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_append_multprec_solution_string
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_number_of_solutions
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_number_of_dobldobl_solutions
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_number_of_quaddobl_solutions
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_number_of_multprec_solutions
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_standard_drop_coordinate_by_index
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_standard_drop_coordinate_by_name
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_dobldobl_drop_coordinate_by_index
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_dobldobl_drop_coordinate_by_name
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_quaddobl_drop_coordinate_by_index
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_quaddobl_drop_coordinate_by_name
 ( PyObject *self, PyObject *args );

/* wrapping functions in product.h starts here */

static PyObject *py2c_product_supporting_set_structure 
 ( PyObject *self, PyObject *args );
static PyObject *py2c_product_write_set_structure
 ( PyObject *self, PyObject *args );
static PyObject *py2c_product_linear_product_root_count
 ( PyObject *self, PyObject *args );
static PyObject *py2c_product_random_linear_product_system
 ( PyObject *self, PyObject *args );
static PyObject *py2c_product_solve_linear_product_system
 ( PyObject *self, PyObject *args );
static PyObject *py2c_product_clear_set_structure
 ( PyObject *self, PyObject *args );

/* wrapping functions in celcon.h starts here */

static PyObject *py2c_celcon_number_of_cells ( PyObject *self, PyObject *args );
static PyObject *py2c_celcon_create_random_coefficient_system 
 ( PyObject *self, PyObject *args );
static PyObject *py2c_celcon_copy_into_systems_container
 ( PyObject *self, PyObject *args );
static PyObject *py2c_celcon_create_polyhedral_homotopy
 ( PyObject *self, PyObject *args );
static PyObject *py2c_celcon_solve_start_system
 ( PyObject *self, PyObject *args );
static PyObject *py2c_celcon_track_solution_path
 ( PyObject *self, PyObject *args );
static PyObject *py2c_celcon_copy_target_solution_to_container
 ( PyObject *self, PyObject *args );
static PyObject *py2c_celcon_permute_system
 ( PyObject *self, PyObject *args );
static PyObject *py2c_celcon_clear_container
 ( PyObject *self, PyObject *args );

/* wrapping functions to manipulate algebraic sets */

static PyObject *py2c_embed_system ( PyObject *self, PyObject *args );
static PyObject *py2c_standard_cascade_homotopy
 ( PyObject *self, PyObject *args );
static PyObject *py2c_dobldobl_cascade_homotopy
 ( PyObject *self, PyObject *args );
static PyObject *py2c_quaddobl_cascade_homotopy
 ( PyObject *self, PyObject *args );
static PyObject *py2c_factor_set_to_mute ( PyObject *self, PyObject *args );
static PyObject *py2c_factor_define_output_file
 ( PyObject *self, PyObject *args );
static PyObject *py2c_factor_assign_labels ( PyObject *self, PyObject *args );
static PyObject *py2c_factor_initialize_sampler
 ( PyObject *self, PyObject *args );
static PyObject *py2c_factor_initialize_monodromy
 ( PyObject *self, PyObject *args );
static PyObject *py2c_factor_store_solutions
 ( PyObject *self, PyObject *args );
static PyObject *py2c_factor_restore_solutions
 ( PyObject *self, PyObject *args );
static PyObject *py2c_factor_track_paths ( PyObject *self, PyObject *args );
static PyObject *py2c_factor_swap_slices ( PyObject *self, PyObject *args );
static PyObject *py2c_factor_new_slices ( PyObject *self, PyObject *args );
static PyObject *py2c_factor_set_trace_slice
 ( PyObject *self, PyObject *args );
static PyObject *py2c_factor_store_gammas ( PyObject *self, PyObject *args );
static PyObject *py2c_factor_permutation_after_loop
 ( PyObject *self, PyObject *args );
static PyObject *py2c_factor_update_decomposition
 ( PyObject *self, PyObject *args );
static PyObject *py2c_factor_number_of_components
 ( PyObject *self, PyObject *args );
static PyObject *py2c_factor_witness_points_of_component
 ( PyObject *self, PyObject *args );
static PyObject *py2c_factor_trace_sum_difference
 ( PyObject *self, PyObject *args );
static PyObject *py2c_witness_set_of_hypersurface
 ( PyObject *self, PyObject *args );
static PyObject *py2c_create_diagonal_homotopy
 ( PyObject *self, PyObject *args );
static PyObject *py2c_start_diagonal_cascade_solutions
 ( PyObject *self, PyObject *args );
static PyObject *py2c_extrinsic_top_diagonal_dimension
 ( PyObject *self, PyObject *args );
static PyObject *py2c_collapse_diagonal ( PyObject *self, PyObject *args );

/* wrapping of Pieri homotopies */

static PyObject *py2c_pieri_count ( PyObject *self, PyObject *args );
static PyObject *py2c_localization_poset ( PyObject *self, PyObject *args );
static PyObject *py2c_pieri_homotopies ( PyObject *self, PyObject *args );
static PyObject *py2c_osculating_planes ( PyObject *self, PyObject *args );
static PyObject *py2c_Schubert_system ( PyObject *self, PyObject *args );

/* wrapping functions in mapcon.h starts from here */

static PyObject *py2c_mapcon_solve_system ( PyObject *self, PyObject *args );
static PyObject *py2c_mapcon_write_maps ( PyObject *self, PyObject *args );
static PyObject *py2c_mapcon_clear_maps ( PyObject *self, PyObject *args );
static PyObject *py2c_mapcon_top_dimension ( PyObject *self, PyObject *args );
static PyObject *py2c_mapcon_number_of_maps ( PyObject *self, PyObject *args );
static PyObject *py2c_mapcon_degree_of_map ( PyObject *self, PyObject *args );
static PyObject *py2c_mapcon_coefficients_of_map
 ( PyObject *self, PyObject *args );
static PyObject *py2c_mapcon_exponents_of_map
 ( PyObject *self, PyObject *args );

/* wrapping functions in next_track starts below */

static PyObject *py2c_initialize_standard_homotopy
 ( PyObject *self, PyObject *args );
static PyObject *py2c_initialize_dobldobl_homotopy
 ( PyObject *self, PyObject *args );
static PyObject *py2c_initialize_quaddobl_homotopy
 ( PyObject *self, PyObject *args );
static PyObject *py2c_initialize_multprec_homotopy
 ( PyObject *self, PyObject *args );
static PyObject *py2c_initialize_standard_solution
 ( PyObject *self, PyObject *args );
static PyObject *py2c_initialize_dobldobl_solution
 ( PyObject *self, PyObject *args );
static PyObject *py2c_initialize_quaddobl_solution
 ( PyObject *self, PyObject *args );
static PyObject *py2c_initialize_multprec_solution
 ( PyObject *self, PyObject *args );
static PyObject *py2c_next_standard_solution
 ( PyObject *self, PyObject *args );
static PyObject *py2c_next_dobldobl_solution
 ( PyObject *self, PyObject *args );
static PyObject *py2c_next_quaddobl_solution
 ( PyObject *self, PyObject *args );
static PyObject *py2c_next_multprec_solution
 ( PyObject *self, PyObject *args );
static PyObject *py2c_clear_standard_tracker
 ( PyObject *self, PyObject *args );
static PyObject *py2c_clear_dobldobl_tracker
 ( PyObject *self, PyObject *args );
static PyObject *py2c_clear_multprec_tracker
 ( PyObject *self, PyObject *args );
static PyObject *py2c_clear_quaddobl_tracker
 ( PyObject *self, PyObject *args );
