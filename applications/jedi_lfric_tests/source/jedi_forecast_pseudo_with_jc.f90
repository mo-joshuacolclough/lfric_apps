!-----------------------------------------------------------------------------
! (C) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @page jedi_forecast_pseudo_with_jc program

!> @brief Main program for running pseudo model forecast with jedi emulator
!>        objects.

!> @details Setup and run a pseudo model forecast using the jedi emulator
!>          objects. The jedi objects are constructed via an initialiser call
!>          and the forecast is handled by the model object.
!>

! Note: This program file represents generic OOPS code and so it should not be
!       edited. If you need to make changes at the program level then please
!       contact darth@metofice.gov.uk for advice.
program jedi_forecast_pseudo_with_jc

  use cli_mod,                      only : parse_command_line
  use config_mod,                   only : config_type
  use constants_mod,                only : PRECISION_REAL, i_def, str_def
  use field_collection_mod,         only : field_collection_type
  use log_mod,                      only : log_event, log_scratch_space, &
                                           LOG_LEVEL_ALWAYS

  ! Jedi emulator objects
  use jedi_checksum_mod,             only : output_checksum
  use jedi_lfric_duration_mod,       only : jedi_duration_type
  use jedi_run_mod,                  only : jedi_run_type
  use jedi_geometry_mod,             only : jedi_geometry_type
  use jedi_state_mod,                only : jedi_state_type
  use jedi_increment_mod,            only : jedi_increment_type
  use jedi_pseudo_model_mod,         only : jedi_pseudo_model_type
  use jedi_post_processor_empty_mod, only : jedi_post_processor_empty_type

  use total_energy_norm_mod,         only : calculate_total_energy_norm

  implicit none

  ! Emulator objects
  type( jedi_geometry_type )             :: jedi_geometry
  type( jedi_state_type )                :: jedi_state
  type( jedi_pseudo_model_type )         :: jedi_psuedo_model
  type( jedi_run_type )                  :: jedi_run
  type( jedi_post_processor_empty_type ) :: jedi_pp_empty

  type( jedi_increment_type )            :: jedi_jc_norm_increment

  ! Local
  type( config_type ), pointer :: config

  character(:), allocatable  :: filename
  integer(i_def)             :: model_communicator
  type(jedi_duration_type)   :: forecast_length
  character(str_def)         :: forecast_length_str

  type( field_collection_type ) :: increment_fields
  type( field_collection_type ) :: state_fields
  character( len=str_def )      :: increment_jc_term_names(6)
  character( len=str_def )      :: state_jc_term_names(2)


  character(*), parameter :: program_name = "jedi_forecast_pseudo_with_jc"

  ! Infrastructure configuration
  call parse_command_line( filename )

  ! Run object - handles initialization and finalization of required
  ! infrastructure. Initialize external libraries such as XIOS
  call jedi_run%initialise( program_name, model_communicator )

  ! Ensemble applications would split the communicator here

  ! Initialize LFRic infrastructure
  call jedi_run%initialise_infrastructure( filename, model_communicator )

  call log_event( 'Running ' // program_name // ' ...', LOG_LEVEL_ALWAYS )
  write(log_scratch_space,'(A)')                        &
        'Application built with '//trim(PRECISION_REAL)// &
        '-bit real numbers'
  call log_event( log_scratch_space, LOG_LEVEL_ALWAYS )

  ! Get the configuration
  config => jedi_run%get_config()

  ! Get the forecast length
  forecast_length_str = config%jedi_lfric_settings%forecast_length()
  call forecast_length%init(forecast_length_str)

  ! Create geometry
  call jedi_geometry%initialise( model_communicator, config )

  ! Create state
  call jedi_state%initialise( jedi_geometry, config )

  ! Model
  call jedi_psuedo_model%initialise( config )

  ! Run non-linear model forecast
  call jedi_psuedo_model%forecast( jedi_state, forecast_length, jedi_pp_empty )

  ! Calculate the total energy norm for the Jc term.
  call jedi_jc_norm_increment%initialise( jedi_geometry, config )
  increment_jc_term_names(1) = "theta_factor"
  increment_jc_term_names(2) = "inv_theta_factor"
  increment_jc_term_names(3) = "pressure_factor"
  increment_jc_term_names(4) = "inv_pressure_factor"
  increment_jc_term_names(5) = "wind_factor"
  increment_jc_term_names(6) = "inv_wind_factor"

  call jedi_jc_norm_increment%get_to_field_collection(increment_jc_term_names, &
                                                      increment_fields)

  state_jc_term_names(1) = "rho"
  state_jc_term_names(2) = "theta"
  call jedi_state%get_to_field_collection(state_jc_term_names, state_fields)

  call calculate_total_energy_norm(state_fields, increment_fields)
  call jedi_jc_norm_increment%set_from_field_collection(increment_jc_term_names, &
                                                        increment_fields)

  ! Write a netCDF via XIOS at the last time step.
  ! Passing state%datetime to state to be consistent with the implementation
  ! in JEDI.
  call jedi_state%write_file( jedi_state%valid_time() )

  ! Print the final state diagnostics
  call jedi_state%print()

  ! To provide KGO
  call output_checksum( program_name, jedi_state%io_collection )

  call log_event( 'Finalising ' // program_name // ' ...', LOG_LEVEL_ALWAYS )

  call jedi_run%finalise()

end program jedi_forecast_pseudo_with_jc
