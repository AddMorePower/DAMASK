!--------------------------------------------------------------------------------------------------
!> @author Sharan Roongta, Max-Planck-Institut für Eisenforschung GmbH
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Finite difference solver for chemical transport equation: FDM
!--------------------------------------------------------------------------------------------------
module grid_chemical_FDM
#include <petsc/finclude/petscsnes.h>
#include <petsc/finclude/petscdmda.h>
  use PETScDMDA
  use PETScSNES
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>14) && !defined(PETSC_HAVE_MPI_F90MODULE_VISIBILITY)
  use MPI_f08
#endif

  use prec
  use parallelization
  use IO
  use misc
  use CLI
  use HDF5
  use HDF5_utilities
  use math
  use rotations
  use spectral_utilities
  use grid_utilities
  use config
  use homogenization
  use discretization
  use discretization_grid
  use constants

#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>14) && !defined(PETSC_HAVE_MPI_F90MODULE_VISIBILITY)
  implicit none(type,external)
#else
  implicit none
#endif
  private

  type(tSolutionParams) :: params

  type :: tNumerics
    integer :: &
      itmin, &                                                                                      !< minimum number of iterations
      itmax                                                                                         !< maximum number of iterations
    real(pREAL) :: &
      eps_chemical_atol, &                                                                          !< absolute tolerance for chemical equilibrium
      eps_chemical_rtol                                                                             !< relative tolerance for chemical equilibrium
  end type tNumerics

  type(tNumerics) :: num                                                                            ! numerics parameters. Better name?

!--------------------------------------------------------------------------------------------------
! PETSc data
  SNES     :: chemical_snes
  Vec      :: mu_PETSc, mu_lastInc_PETSc

  type, private :: tScal
    PetscScalar, pointer :: scal(:,:,:,:)
  end type tScal

  real(pReal), private :: delta(3)
  real(pReal), dimension(:,:,:,:), allocatable :: &
    conc_current, &
    conc_lastInc, &
    conc_stagInc

  integer :: totalIter = 0                                                                          !< total iteration in current increment
  integer :: N_components

  public :: &
    grid_chemical_FDM_init, &
    grid_chemical_FDM_solution, &
    grid_chemical_FDM_forward

contains

!--------------------------------------------------------------------------------------------------
!> @brief allocates all neccessary fields and fills them with data
! ToDo: Restart not implemented
!--------------------------------------------------------------------------------------------------
subroutine grid_chemical_FDM_init()

end subroutine grid_chemical_FDM_init

!--------------------------------------------------------------------------------------------------
!> @brief solution for the FDM chemical scheme with internal iterations
!--------------------------------------------------------------------------------------------------
function grid_chemical_FDM_solution(Delta_t) result(solution)

  real(pREAL), intent(in) :: Delta_t                                                                !< increment in time for current solution

  type(tSolutionState) :: solution

end function grid_chemical_FDM_solution

!--------------------------------------------------------------------------------------------------
!> @brief Set DAMASK data to current solver status.
!--------------------------------------------------------------------------------------------------
subroutine grid_chemical_FDM_forward(cutBack)

  logical, intent(in) :: cutBack

end subroutine grid_chemical_FDM_forward

end module grid_chemical_FDM
