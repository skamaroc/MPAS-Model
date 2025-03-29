! Copyright (c) 2016,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module mpas_basic_state_diagnostics

    use mpas_dmpar
    use mpas_kind_types
    use mpas_derived_types
    use mpas_pool_routines
    use mpas_constants
    use mpas_log, only : mpas_log_write

    type (MPAS_pool_type), pointer :: mesh
    type (MPAS_pool_type), pointer :: state
    type (MPAS_pool_type), pointer :: diag

    type (MPAS_clock_type), pointer :: clock

    public :: basic_state_diagnostics_setup, &
              basic_state_diagnostics_compute

    private

    logical :: need_temperature

    contains

    !-----------------------------------------------------------------------
    !  routine basic_state_diagnostics_setup
    !
    !> \brief Set up the basic state diagnostics module
    !> \author Bill Skamarock
    !> \date   27 March 2025
    !> \details
    !>  This routine sets up the basic state diagnostics module, principally by
    !>  saving pointers to pools that are used in the computation of diagnostics.
    !
    !-----------------------------------------------------------------------
    subroutine basic_state_diagnostics_setup(all_pools, simulation_clock)

        use mpas_derived_types, only : MPAS_pool_type, MPAS_clock_type
        use mpas_pool_routines, only : mpas_pool_get_subpool

        implicit none

        type (MPAS_pool_type), pointer :: all_pools
        type (MPAS_clock_type), pointer :: simulation_clock

        clock => simulation_clock

        call mpas_pool_get_subpool(all_pools, 'mesh', mesh)
        call mpas_pool_get_subpool(all_pools, 'state', state)
        call mpas_pool_get_subpool(all_pools, 'diag', diag)
   
    end subroutine basic_state_diagnostics_setup


    !-----------------------------------------------------------------------
    !  routine basic_state_diagnostics_compute
    !
    !> \brief Compute basic state diagnostic before model output is written
    !> \author Bill Skamarock
    !> \date   27 March 2025
    !> \details
    !>  Compute basic state diagnostic before model output is written. 
    !
    !-----------------------------------------------------------------------
    subroutine basic_state_diagnostics_compute()

        use mpas_atm_diagnostics_utils, only : MPAS_field_will_be_written

        implicit none

        logical :: need_any_diags

        need_any_diags = .false.
        need_temperature = .false.

        ! call mpas_log_write(' checking if temperature diagnostic is needed ')

        need_temperature = MPAS_field_will_be_written('temperature')
        need_any_diags = need_any_diags .or. need_temperature
        if (need_any_diags) then
            call compute_diagnostics(mesh, state, 1, diag)
        end if
   
    end subroutine basic_state_diagnostics_compute


   !==================================================================================================
    subroutine compute_diagnostics(mesh, state, time_lev, diag)
   !==================================================================================================

       !input arguments:
        type (mpas_pool_type), intent(in)  :: mesh
        type (mpas_pool_type), intent(in) :: state
        integer, intent(in) :: time_lev              ! which time level to use from state
       
       !inout arguments:
        type (mpas_pool_type), intent(inout) :: diag
       
        integer :: iCell,k
        integer, pointer :: nCells, nVertLevels
        integer, pointer :: index_qv
       
        real (kind=RKIND), dimension(:,:), pointer :: exner
        real (kind=RKIND), dimension(:,:), pointer :: theta_m
        real (kind=RKIND), dimension(:,:,:), pointer :: scalars
        real (kind=RKIND), dimension(:,:), pointer :: temperature
       
       !--------------------------------------------------------------------------------------------------
       
       ! call mpas_log_write('')
       ! call mpas_log_write('--- enter subroutine interp_diagnostics:')
       
        call mpas_pool_get_dimension(mesh, 'nCells', nCells)
        call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
        call mpas_pool_get_dimension(state, 'index_qv', index_qv)
        call mpas_pool_get_array(state, 'theta_m', theta_m, time_lev)
        call mpas_pool_get_array(state, 'scalars', scalars, time_lev)
        call mpas_pool_get_array(diag, 'exner', exner)
       
        call mpas_pool_get_array(diag, 'temperature', temperature)

        if (need_temperature) then
            ! call mpas_log_write(' computing temperature diagnostic ')
            !calculation of temperature at cell centers:
            do iCell = 1,nCells
            do k = 1,nVertLevels
                temperature(k,iCell) = exner(k,iCell)*(theta_m(k,iCell)/(1._RKIND+rvord*scalars(index_qv,k,iCell)))
            enddo
            enddo
         end if

    end subroutine compute_diagnostics
       
end module mpas_basic_state_diagnostics
