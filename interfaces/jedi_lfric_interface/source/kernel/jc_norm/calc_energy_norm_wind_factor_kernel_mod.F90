!-----------------------------------------------------------------------------
! (C) Crown copyright 2026 Met Office. All rights reserved.
! For further details please refer to the file LICENCE which you should have
! received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief    Calculates pre-factor for wind increment
!!           for the total energy norm of the Jc term.
!> @details  Uses central ls density values + cell volumes on W3 to calculate
!!           the pre-factor for the wind increment on W3.
!!           Currently valid for horizontal and vertical element orders = 0.
module calc_energy_norm_wind_factor_kernel_mod

  use argument_mod,         only : arg_type, func_type,        &
                                   GH_FIELD, GH_SCALAR,        &
                                   GH_REAL, GH_READ, GH_WRITE, &
                                   CELL_COLUMN
  use fs_continuity_mod,    only : W3
  use constants_mod,        only : r_def, i_def, EPS
  use kernel_mod,           only : kernel_type

  implicit none

  private

  !> The type declaration for the kernel. Contains the metadata needed by
  !> the PSy layer.
  !>
  type, public, extends(kernel_type) :: calc_energy_norm_wind_factor_kernel_type
    private
    type(arg_type) :: meta_args(5) = (/                   &
         arg_type(GH_FIELD,   GH_REAL, GH_WRITE, W3),     & ! factor
         arg_type(GH_FIELD,   GH_REAL, GH_WRITE, W3),     & ! inv_factor
         arg_type(GH_FIELD,   GH_REAL, GH_READ,  W3),     & ! ls_rho
         arg_type(GH_FIELD,   GH_REAL, GH_READ,  W3),     & ! cell_vol
         arg_type(GH_SCALAR,  GH_REAL, GH_READ)           & ! area_domain
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: calc_energy_norm_wind_factor_kernel_code
  end type

  !-------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-------------------------------------------------------------------------
  public :: calc_energy_norm_wind_factor_kernel_code

contains

  !! @param[in]      nlayers        Number of layers
  !! @param[in,out]  factor         W3 field containing the factor of the diagonal norm matrix
  !!                                pertaining to wind increments
  !! @param[in,out]  inv_factor     W3 field containing the factor of the inverse diagonal norm matrix
  !!                                pertaining to wind increments
  !! @param[in]      ls_rho         Linearisation state density
  !! @param[in]      cell_vol       Cell volumes on W3
  !! @param[in]      area_domain    Surface area of the model domain
  !! @param[in]      ndf_w3         Number of degrees of freedom per cell for W3
  !! @param[in]      undf_w3        Number of degrees of freedom for W3
  !! @param[in]      map_w3         Dofmap for the cell at the base of the column for W3
  subroutine calc_energy_norm_wind_factor_kernel_code(nlayers,      &
                                                      factor,       &
                                                      inv_factor,   &
                                                      ls_rho,       &
                                                      cell_vol,     &
                                                      area_domain,  &
                                                      ndf_w3,       &
                                                      undf_w3,      &
                                                      map_w3)

    implicit none

    ! Arguments
    integer(kind=i_def),                             intent(in)    :: nlayers
    integer(kind=i_def),                             intent(in)    :: ndf_w3
    integer(kind=i_def),                             intent(in)    :: undf_w3
    integer(kind=i_def), dimension(ndf_w3),          intent(in)    :: map_w3
    real(kind=r_def), dimension(undf_w3),            intent(inout) :: factor
    real(kind=r_def), dimension(undf_w3),            intent(inout) :: inv_factor
    real(kind=r_def), dimension(undf_w3),            intent(in)    :: ls_rho
    real(kind=r_def), dimension(undf_w3),            intent(in)    :: cell_vol
    real(kind=r_def),                                intent(in)    :: area_domain

    ! Internal variables
    real(kind=r_def)                           :: ls_rho_at_cell
    real(kind=r_def)                           :: volume_at_cell
    integer(kind=i_def)                        :: df_w3
    integer(kind=i_def)                        :: k

    do k = 0, nlayers - 1
      ! Use central values for density and volume
      ls_rho_at_cell = ls_rho(map_w3(1) + k) + sign(EPS, ls_rho(map_w3(1) + k))
      volume_at_cell = cell_vol(map_w3(1) + k) + sign(EPS, cell_vol(map_w3(1) + k))

      do df_w3 = 1, ndf_w3
        factor(map_w3(df_w3) + k) = volume_at_cell*ls_rho_at_cell / &
                                    (2.0_r_def*(area_domain + sign(EPS, area_domain)))
        inv_factor(map_w3(df_w3) + k) = 1.0_r_def / (factor(map_w3(df_w3) + k) + &
                                                     sign(EPS, factor(map_w3(df_w3) + k)))
      end do
    end do

  end subroutine calc_energy_norm_wind_factor_kernel_code

end module calc_energy_norm_wind_factor_kernel_mod
