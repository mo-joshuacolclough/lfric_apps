!-----------------------------------------------------------------------------
! (C) Crown copyright 2026 Met Office. All rights reserved.
! For further details please refer to the file LICENCE which you should have
! received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief    Calculates pre-factor for potential temperature increment
!!           for the total energy norm of the Jc term.
!> @details  Uses central ls density values + cell volumes on W3 to calculate
!!           the pre-factor for the potential temperature increment on Wtheta.
!!           Currently valid for horizontal and vertical element orders = 0.
module calc_energy_norm_theta_factor_kernel_mod

  use argument_mod,         only : arg_type, func_type,        &
                                   GH_FIELD, GH_SCALAR,        &
                                   GH_REAL, GH_READ, GH_WRITE, &
                                   CELL_COLUMN
  use fs_continuity_mod,    only : W3, Wtheta
  use constants_mod,        only : r_def, i_def, EPS
  use kernel_mod,           only : kernel_type

  implicit none

  private

  !> The type declaration for the kernel. Contains the metadata needed by
  !> the PSy layer.
  !>
  type, public, extends(kernel_type) :: calc_energy_norm_theta_factor_kernel_type
    private
    type(arg_type) :: meta_args(8) = (/                   &
         arg_type(GH_FIELD,   GH_REAL, GH_WRITE, Wtheta), & ! factor
         arg_type(GH_FIELD,   GH_REAL, GH_WRITE, Wtheta), & ! inv_factor
         arg_type(GH_FIELD,   GH_REAL, GH_READ,  Wtheta), & ! ls_theta
         arg_type(GH_FIELD,   GH_REAL, GH_READ,  W3),     & ! ls_rho
         arg_type(GH_FIELD,   GH_REAL, GH_READ,  W3),     & ! cell_vol
         arg_type(GH_SCALAR,  GH_REAL, GH_READ),          & ! area_domain
         arg_type(GH_SCALAR,  GH_REAL, GH_READ),          & ! grav_const
         arg_type(GH_SCALAR,  GH_REAL, GH_READ)           & ! buoyancy_freq
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: calc_energy_norm_theta_factor_kernel_code
  end type

  !-------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-------------------------------------------------------------------------
  public :: calc_energy_norm_theta_factor_kernel_code

contains

  !! @param[in]      nlayers        Number of layers
  !! @param[in,out]  factor         Wtheta field containing the factor of the diagonal norm matrix
  !!                                pertaining to potential temperature increments
  !! @param[in,out]  inv_factor     Wtheta field containing the factor of the inverse diagonal norm matrix
  !!                                pertaining to potential temperature increments
  !! @param[in]      ls_theta       Linearisation state potential temperature
  !! @param[in]      ls_rho         Linearisation state density
  !! @param[in]      cell_vol       Cell volumes on W3
  !! @param[in]      area_domain    Surface area of the model domain
  !! @param[in]      grav_const     Gravitational constant for the planet
  !! @param[in]      buoyancy_freq  Hydrostatic potential temperature buoyancy frequency
  !! @param[in]      ndf_wth        Number of degrees of freedom per cell for Wtheta
  !! @param[in]      undf_wth       Total number of degrees of freedom for Wtheta
  !! @param[in]      map_wth        Dofmap for the cell at the base of the column for Wtheta
  !! @param[in]      ndf_w3         Number of degrees of freedom per cell for W3
  !! @param[in]      undf_w3        Number of degrees of freedom for W3
  !! @param[in]      map_w3         Dofmap for the cell at the base of the column for W3
  subroutine calc_energy_norm_theta_factor_kernel_code(nlayers,       &
                                           factor,        &
                                           inv_factor,    &
                                           ls_theta,      &
                                           ls_rho,        &
                                           cell_vol,      &
                                           area_domain,   &
                                           grav_const,    &
                                           buoyancy_freq, &
                                           ndf_wth,       &
                                           undf_wth,      &
                                           map_wth,       &
                                           ndf_w3,        &
                                           undf_w3,       &
                                           map_w3)

    implicit none

    ! Arguments
    integer(kind=i_def),                             intent(in)    :: nlayers
    integer(kind=i_def),                             intent(in)    :: ndf_wth
    integer(kind=i_def),                             intent(in)    :: ndf_w3
    integer(kind=i_def),                             intent(in)    :: undf_wth
    integer(kind=i_def),                             intent(in)    :: undf_w3
    integer(kind=i_def), dimension(ndf_wth),         intent(in)    :: map_wth
    integer(kind=i_def), dimension(ndf_w3),          intent(in)    :: map_w3
    real(kind=r_def), dimension(undf_wth),           intent(inout) :: factor
    real(kind=r_def), dimension(undf_wth),           intent(inout) :: inv_factor
    real(kind=r_def), dimension(undf_wth),           intent(in)    :: ls_theta
    real(kind=r_def), dimension(undf_w3),            intent(in)    :: ls_rho
    real(kind=r_def), dimension(undf_w3),            intent(in)    :: cell_vol
    real(kind=r_def),                                intent(in)    :: area_domain
    real(kind=r_def),                                intent(in)    :: grav_const
    real(kind=r_def),                                intent(in)    :: buoyancy_freq

    ! Internal variables
    real(kind=r_def)                           :: ls_rho_at_cell
    real(kind=r_def)                           :: volume_at_cell
    integer(kind=i_def)                        :: df_wth
    integer(kind=i_def)                        :: k

    do k = 0, nlayers - 1
      ! Use central values for density and volume
      ls_rho_at_cell = ls_rho(map_w3(1) + k) + sign(EPS, ls_rho(map_w3(1) + k))
      volume_at_cell = cell_vol(map_w3(1) + k) + sign(EPS, cell_vol(map_w3(1) + k))

      do df_wth = 1, ndf_wth
        factor(map_wth(df_wth) + k) = grav_const / ((ls_theta(map_wth(df_wth) + k) &
                                                     + sign(EPS, ls_theta(map_wth(df_wth) + k)))*buoyancy_freq)
        factor(map_wth(df_wth) + k) = ls_rho_at_cell*volume_at_cell*factor(map_wth(df_wth) + k)**2 / &
                                      (2.0_r_def*(area_domain + sign(EPS, area_domain)))
        inv_factor(map_wth(df_wth) + k) = 1.0_r_def/(factor(map_wth(df_wth) + k) + &
                                                     sign(EPS, factor(map_wth(df_wth) + k)))
      end do
    end do

  end subroutine calc_energy_norm_theta_factor_kernel_code

end module calc_energy_norm_theta_factor_kernel_mod
