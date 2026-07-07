!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @page jedi_tlm_forecast_tl_with_jc program

!> @brief Main program for running linear forecast with jedi emulator
!>        objects.
!>
!> @details Setup and run a linear model forecastTL using the JEDI
!>          emulator objects. The linear state trajectory is provided via the
!>          pseudo model forecast. The jedi objects are constructed via an
!>          initialiser call and the forecasts are handled by the model
!>          objects.
!>

! Note: This program file represents generic OOPS code and so it should not be
!       edited. If you need to make changes at the program level then please
!       contact darth@metofice.gov.uk for advice.
program jedi_tlm_forecast_tl_with_jc

  use cli_mod,                      only : parse_command_line
  use config_mod,                   only : config_type
  use constants_mod,                only : PRECISION_REAL, i_def, str_def
  use field_collection_mod,         only : field_collection_type
  use mesh_mod,                     only : mesh_type
  use log_mod,                      only : log_event, log_scratch_space, &
                                           LOG_LEVEL_ALWAYS

  ! Jedi emulator objects
  use jedi_checksum_mod,            only : output_linear_checksum
  use jedi_lfric_duration_mod,      only : jedi_duration_type
  use jedi_run_mod,                 only : jedi_run_type
  use jedi_geometry_mod,            only : jedi_geometry_type
  use jedi_state_mod,               only : jedi_state_type
  use jedi_increment_mod,           only : jedi_increment_type
  use jedi_pseudo_model_mod,        only : jedi_pseudo_model_type
  use jedi_linear_model_mod,        only : jedi_linear_model_type
  use jedi_post_processor_traj_mod, only : jedi_post_processor_traj_type
  use jedi_field_utils_mod,          only : populate_field_collection

  use total_energy_norm_mod,         only : calculate_total_energy_norm

  implicit none

  ! Emulator objects
  type( jedi_geometry_type )            :: jedi_geometry
  type( jedi_state_type )               :: jedi_state
  type( jedi_increment_type )           :: jedi_increment
  type( jedi_pseudo_model_type )        :: jedi_psuedo_model
  type( jedi_linear_model_type )        :: jedi_linear_model
  type( jedi_run_type )                 :: jedi_run
  type( jedi_post_processor_traj_type ) :: pp_traj

  ! Local
  type( config_type ), pointer :: config

  character(:), allocatable  :: filename
  integer( kind=i_def )      :: model_communicator
  type( jedi_duration_type ) :: forecast_length
  character( str_def )       :: forecast_length_str

  type(mesh_type), pointer      :: mesh3d
  type(mesh_type), pointer      :: mesh2d
  type( field_collection_type ) :: jc_increment_fields
  type( field_collection_type ) :: jc_state_fields
  character( len=str_def )      :: jc_increment_term_names(6)
  character( len=str_def )      :: jc_state_term_names(2)


  character(*), parameter :: program_name = "jedi_tlm_forecast_tl_with_jc"

  ! Infrastructure config
  call parse_command_line( filename )

  call log_event( 'Running ' // program_name // ' ...', LOG_LEVEL_ALWAYS )
  write(log_scratch_space,'(A)')                        &
        'Application built with '//trim(PRECISION_REAL)// &
        '-bit real numbers'
  call log_event( log_scratch_space, LOG_LEVEL_ALWAYS )

  ! Run object - handles initialization and finalization of required
  ! infrastructure. Initialize external libraries such as XIOS
  call jedi_run%initialise( program_name, model_communicator )

  ! Ensemble applications would split the communicator here

  ! Initialize LFRic infrastructure
  call jedi_run%initialise_infrastructure( filename, model_communicator )

  ! Get the configuration
  config => jedi_run%get_config()

  ! Get the forecast length
  forecast_length_str = config%jedi_lfric_settings%forecast_length()
  call forecast_length%init(forecast_length_str)

  ! Create geometry
  call jedi_geometry%initialise( model_communicator, config )

  ! Create state
  call jedi_state%initialise( jedi_geometry, config )

  ! Create increment
  call jedi_increment%initialise( jedi_geometry, config )

  ! Create linear model
  call jedi_linear_model%initialise( jedi_geometry, filename )

  ! Initialise trajectory post processor with instance of jedi_linear_model
  call pp_traj%initialise( jedi_linear_model )

  ! Create non-linear model
  call jedi_psuedo_model%initialise( config )

  ! Run non-linear model forecast to populate the trajectory object
  call jedi_psuedo_model%forecast( jedi_state, forecast_length, pp_traj )

  ! Run the linear model TL forecast
  call jedi_linear_model%forecastTL( jedi_increment, forecast_length )

  ! Calculate JC term norm.
  ! = Get required fields into a field collection
  jc_increment_term_names(1) = "theta_factor"
  jc_increment_term_names(2) = "inv_theta_factor"
  jc_increment_term_names(3) = "pressure_factor"
  jc_increment_term_names(4) = "inv_pressure_factor"
  jc_increment_term_names(5) = "wind_factor"
  jc_increment_term_names(6) = "inv_wind_factor"
  call jc_increment_fields%initialise(name = "jc_increment_fields", table_len=100)
  mesh3d => jedi_geometry%get_mesh()
  mesh2d => jedi_geometry%get_twod_mesh()
  call populate_field_collection(mesh3d, mesh2d, jc_increment_term_names, jc_increment_fields)
  
  call jedi_increment%get_to_field_collection(jc_increment_term_names, &
                                              jc_increment_fields)

  ! = Create state fields in field collection.
  jc_state_term_names(1) = "rho"
  jc_state_term_names(2) = "theta"
  call jc_state_fields%initialise(name = "jc_state_fields", table_len=100)
  call populate_field_collection(mesh3d, mesh2d, jc_state_term_names, jc_state_fields)

  call jedi_state%get_to_field_collection(jc_state_term_names, jc_state_fields)

  ! = Run calculation and set the increment
  call calculate_total_energy_norm(config, jc_state_fields, jc_increment_fields)
  call jedi_increment%set_from_field_collection(jc_increment_term_names, &
                                                jc_increment_fields)

  ! Print the final state and increment diagnostics
  call jedi_state%print()
  call jedi_increment%print()

  ! To provide KGO
  call output_linear_checksum( program_name, jedi_linear_model%modeldb )

  call log_event( 'Finalising ' // program_name // ' ...', LOG_LEVEL_ALWAYS )

  call jedi_run%finalise()

end program jedi_tlm_forecast_tl_with_jc
