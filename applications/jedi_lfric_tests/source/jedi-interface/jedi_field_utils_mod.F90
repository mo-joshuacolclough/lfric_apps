!-----------------------------------------------------------------------------
! (C) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @details A subroutine is included to setup a field meta data object for a
!>          given list of variable names. The module includes a list of the
!>          presently supported fields that define the required meta-data that
!>          is needed to setup a field:
!>          i)  Function-space enumerator
!>          ii) Logical defining if the field is 2D or 3D
!>
module jedi_field_utils_mod

  use constants_mod,             only : i_def, str_def, l_def
  use fs_continuity_mod,         only : W3, Wtheta
  use jedi_lfric_datetime_mod,   only : jedi_datetime_type
  use jedi_lfric_field_meta_mod, only : jedi_lfric_field_meta_type
  use log_mod,                   only : log_event,          &
                                        log_scratch_space,  &
                                        LOG_LEVEL_ERROR

  implicit none

  private

  public :: setup_field_meta_data
  public :: populate_field_collection

  contains

  !> @brief    Setup field meta data using a list of variables
  !>
  !> @param [out] field_meta_data A field meta data object
  !> @param [in]  variable_names  A list of variable names
  subroutine setup_field_meta_data(field_meta_data, variable_names)

  type( jedi_lfric_field_meta_type ), intent(out) :: field_meta_data
  character( len=str_def ),           intent(in)  :: variable_names(:)

  ! Local
  logical( kind=l_def ), allocatable :: variable_is_2d(:)
  integer( kind=i_def ), allocatable :: variable_function_spaces(:)
  integer( kind=i_def )              :: ivar
  integer( kind=i_def )              :: nvars

  ! Setup arrays
  nvars = size(variable_names)
  allocate( variable_function_spaces(nvars), variable_is_2d(nvars) )

  ! Get the field info
  do ivar = 1, nvars
    call get_field_info( variable_function_spaces(ivar), &
                         variable_is_2d(ivar),           &
                         variable_names(ivar) )
  enddo

  ! Setup field_meta_data
  call field_meta_data%initialise( variable_names,           &
                                   variable_function_spaces, &
                                   variable_is_2d )

  end subroutine setup_field_meta_data

  !> @brief    Get field meta-data for a given variable name
  !>
  !> @param [out] function_space The function space enumerator
  !> @param [out] is_2d          Logical defineing if the field is 2D
  !> @param [in]  variable_name  The variable name
  subroutine get_field_info(function_space, is_2d, variable_name)

    integer( kind=i_def ),    intent(out) :: function_space
    logical( kind=l_def ),    intent(out) :: is_2d
    character( len=str_def ), intent(in)  :: variable_name

    ! Return the function_space and is_2d for a given variable name:
    select case ( variable_name )
      case ( "theta" )
        function_space = Wtheta
        is_2d = .false.
      case ( "theta_factor" )
        function_space = Wtheta
        is_2d = .false.
      case ( "inv_theta_factor" )
        function_space = Wtheta
        is_2d = .false.
      case ( "rho" )
        function_space = W3
        is_2d = .false.
      case ( "exner" )
        function_space = W3
        is_2d = .false.
      case ( "u_in_w3" )
        function_space = W3
        is_2d = .false.
      case ( "v_in_w3" )
        function_space = W3
        is_2d = .false.
      case ( "w_in_wth" )
        function_space = Wtheta
        is_2d = .false.
      case ( "m_v" )
        function_space = Wtheta
        is_2d = .false.
      case ( "m_cl" )
        function_space = Wtheta
        is_2d = .false.
      case ( "m_r" )
        function_space = Wtheta
        is_2d = .false.
      case ( "m_s" )
        function_space = Wtheta
        is_2d = .false.
      case ( "u10m" )
        function_space = W3
        is_2d = .true.
      case ( "land_fraction" )
        function_space = W3
        is_2d = .true.
      case ( "pressure_factor" )
        function_space = W3
        is_2d = .false.
      case ( "inv_pressure_factor" )
        function_space = W3
        is_2d = .false.
      case ( "wind_factor" )
        function_space = W3
        is_2d = .false.
      case ( "inv_wind_factor" )
        function_space = W3
        is_2d = .false.
      case default
        write ( log_scratch_space, '(4A)' )                          &
                "jedi_field_utils_mod::get_field_info:: ", &
                "The variable name: ",                               &
                trim(variable_name),                                 &
                " is not yet supported."
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end select

  end subroutine get_field_info

! ------------------------------------------------------------------------------

  !> @brief    Populate a field_collection with a list of given fields.
  !>
  !> @param [in] mesh3d  3D Mesh object to create fields on.
  !> @param [in] mesh2d  2D Mesh object to create fields on.
  !> @param [in] var_names  A list of variable names
  !> @param [inout] field_collection  The field collection to populate.
  subroutine populate_field_collection(mesh3d, mesh2d, var_names, field_collection)

    use function_space_mod,                only : function_space_type
    use field_collection_mod,              only : field_collection_type
    use field_mod,                         only : field_type
    use finite_element_config_mod,         only : element_order_h, element_order_v
    use mesh_mod,                          only : mesh_type

    use jedi_lfric_utils_mod,              only : add_real_field

    type(mesh_type),     pointer, intent(in) :: mesh3d
    type(mesh_type),     pointer, intent(in) :: mesh2d
    character( len=str_def ),    intent(in)  :: var_names(:)
    type(field_collection_type), intent(inout) :: field_collection

    ! Local
    type(mesh_type), pointer           :: mesh
    integer(kind=i_def)                :: fspace_enum
    logical(kind=l_def)                :: is_2d
    integer(kind=i_def)                :: i

    do i = 1, size(var_names, dim=1)
      call get_field_info(fspace_enum, is_2d, var_names(i))

      if (is_2d) then
        mesh => mesh2d
      else
        mesh => mesh3d
      end if

      call add_real_field(field_collection, mesh, fspace_enum, var_names(i))
    end do

  end subroutine populate_field_collection

end module jedi_field_utils_mod
