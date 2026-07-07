!-----------------------------------------------------------------------------
! (C) Crown copyright 2026 Met Office. All rights reserved.
! For further details please refer to the file LICENCE which you should have
! received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief    Calculates pre-factor for pressure increment
!!           for the total energy norm of the Jc term.
!> @details  Uses central ls density values + cell volumes on W3 to calculate
!!           the pre-factor for the pressure increment on T.B.D (W3?).
!!           Currently valid for horizontal and vertical element orders = 0.
module calc_energy_norm_pressure_factor_kernel_mod

  use argument_mod,         only : arg_type, func_type,        &
                                   GH_FIELD, GH_SCALAR,        &
                                   GH_REAL, GH_READ, GH_WRITE, &
                                   ANY_SPACE_1,                &
                                   CELL_COLUMN
  use fs_continuity_mod,    only : W3, Wtheta
  use constants_mod,        only : r_def, i_def, EPS
  use kernel_mod,           only : kernel_type

  implicit none

  private

  !> The type declaration for the kernel. Contains the metadata needed by
  !> the PSy layer.
  !>
  type, public, extends(kernel_type) :: calc_energy_norm_pressure_factor_kernel_type
    private
    type(arg_type) :: meta_args(7) = (/                        &
         arg_type(GH_FIELD,   GH_REAL, GH_WRITE, ANY_SPACE_1), & ! factor
         arg_type(GH_FIELD,   GH_REAL, GH_WRITE, ANY_SPACE_1), & ! inv_factor
         arg_type(GH_FIELD,   GH_REAL, GH_READ,  W3),          & ! ls_rho
         arg_type(GH_FIELD,   GH_REAL, GH_READ,  W3),          & ! cell_vol
         arg_type(GH_FIELD,   GH_REAL, GH_READ,  W3),          & ! avg_radh
         arg_type(GH_SCALAR,  GH_REAL, GH_READ),               & ! area_domain
         arg_type(GH_SCALAR,  GH_REAL, GH_READ)                & ! sound_speed
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: calc_energy_norm_pressure_factor_kernel_code
  end type

  !-------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-------------------------------------------------------------------------
  public :: calc_energy_norm_pressure_factor_kernel_code

contains

  !! @param[in]      nlayers        Number of layers
  !! @param[in,out]  factor         Field containing the factor of the diagonal norm matrix
  !!                                pertaining to pressure increments
  !! @param[in,out]  inv_factor     Field containing the factor of the inverse of the diagonal norm matrix
  !!                                pertaining to pressure increments
  !! @param[in]      ls_rho         Linearisation state density
  !! @param[in]      cell_vol       Cell volumes on W3
  !! @param[in]      avg_radh       Cell radial heights on W3
  !! @param[in]      area_domain    Surface area of the model domain
  !! @param[in]      sound_speed    Speed of sound
  !! @param[in]      ndf_aspc1      Number of degrees of freedom per cell for factor
  !! @param[in]      undf_aspc1     Total number of degrees of freedom for factor
  !! @param[in]      map_aspc1      Dofmap for the cell at the base of the column for factor
  !! @param[in]      ndf_w3         Number of degrees of freedom per cell for W3
  !! @param[in]      undf_w3        Number of degrees of freedom for W3
  !! @param[in]      map_w3         Dofmap for the cell at the base of the column for W3
  subroutine calc_energy_norm_pressure_factor_kernel_code(nlayers,       &
                                                          factor,        &
                                                          inv_factor,    &
                                                          ls_rho,        &
                                                          cell_vol,      &
                                                          avg_radh,      &
                                                          area_domain,   &
                                                          sound_speed,   &
                                                          ndf_aspc1,     &
                                                          undf_aspc1,    &
                                                          map_aspc1,     &
                                                          ndf_w3,        &
                                                          undf_w3,       &
                                                          map_w3)

    implicit none

    ! Arguments
    integer(kind=i_def),                             intent(in)    :: nlayers
    integer(kind=i_def),                             intent(in)    :: ndf_aspc1
    integer(kind=i_def),                             intent(in)    :: ndf_w3
    integer(kind=i_def),                             intent(in)    :: undf_aspc1
    integer(kind=i_def),                             intent(in)    :: undf_w3
    integer(kind=i_def), dimension(ndf_aspc1),       intent(in)    :: map_aspc1
    integer(kind=i_def), dimension(ndf_w3),          intent(in)    :: map_w3
    real(kind=r_def), dimension(undf_aspc1),         intent(inout) :: factor
    real(kind=r_def), dimension(undf_aspc1),         intent(inout) :: inv_factor
    real(kind=r_def), dimension(undf_w3),            intent(in)    :: ls_rho
    real(kind=r_def), dimension(undf_w3),            intent(in)    :: cell_vol
    real(kind=r_def), dimension(undf_w3),            intent(in)    :: avg_radh
    real(kind=r_def),                                intent(in)    :: area_domain
    real(kind=r_def),                                intent(in)    :: sound_speed

    ! Internal variables
    real(kind=r_def)                           :: ls_rho_at_cell
    real(kind=r_def)                           :: volume_at_cell
    real(kind=r_def)                           :: radh_at_cell
    integer(kind=i_def)                        :: df_aspc1
    integer(kind=i_def)                        :: k

    do k = 0, nlayers - 1
      ! Use central values for density, volume and cell radial height
      ls_rho_at_cell = ls_rho(map_w3(1) + k) + sign(EPS, ls_rho(map_w3(1) + k))
      volume_at_cell = cell_vol(map_w3(1) + k) + sign(EPS, cell_vol(map_w3(1) + k))
      radh_at_cell = avg_radh(map_w3(1) + k) + sign(EPS, avg_radh(map_w3(1) + k))

      do df_aspc1 = 1, ndf_aspc1
        factor(map_aspc1(df_aspc1) + k) = (radh_at_cell**2*volume_at_cell) / &
                                          (2.0_r_def*(area_domain + sign(EPS, area_domain))*ls_rho_at_cell*sound_speed**2)
        inv_factor(map_aspc1(df_aspc1) + k) = 1.0_r_def/(factor(map_aspc1(df_aspc1) + k) + &
                                                         sign(EPS, factor(map_aspc1(df_aspc1) + k)))
      end do
    end do

  end subroutine calc_energy_norm_pressure_factor_kernel_code

end module calc_energy_norm_pressure_factor_kernel_mod
