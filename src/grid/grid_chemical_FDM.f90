module grid_chemical_FDM
#include <petsc/finclude/petscsnes.h>
#include <petsc/finclude/petscdmda.h>
  use PETScDMDA
  use PETScSNES
#ifndef PETSC_HAVE_MPI_F90MODULE_VISIBILITY
  use MPI_f08
#endif

  use prec
  use parallelization
  use IO
  use spectral_utilities
  use discretization_grid
  use homogenization
  use material
  use types
  use config
  use math

#ifndef PETSC_HAVE_MPI_F90MODULE_VISIBILITY
  implicit none(type,external)
#else
  implicit none
#endif
  private

  type :: tNumerics
    integer :: &
      itmax                                                                                         !< maximum number of iterations
    real(pREAL) :: &
      eps_chemical_atol, &                                                                          !< absolute tolerance for thermal equilibrium
      eps_chemical_rtol                                                                             !< relative tolerance for thermal equilibrium
  end type tNumerics

  type(tNumerics) :: num

!--------------------------------------------------------------------------------------------------
! PETSc data
  SNES :: SNES_chemical
  Vec :: solution_current, solution_lastInc

  type, private :: tScal
    PetscScalar, pointer :: scal(:,:,:,:)
  end type tScal

  real(pREAL), private :: delta(3)
  real(pREAL), dimension(:,:,:,:), allocatable :: &
    conc_current, &
    conc_lastInc, &
    conc_stagInc

  integer :: totalIter = 0                                                                          !< total iteration in current increment
  integer :: N_components
  real(pREAL) :: Delta_t_

  public :: &
    grid_chemical_FDM_init, &
    grid_chemical_FDM_solution, &
    grid_chemical_FDM_forward

contains

!--------------------------------------------------------------------------------------------------
!> @brief Allocate all necessary fields and fill them with data.
! ToDo: Restart not implemented
!--------------------------------------------------------------------------------------------------
subroutine grid_chemical_FDM_init(num_grid)

  type(tDict), pointer, intent(in) :: num_grid

  PetscInt, dimension(0:worldsize-1) :: cells3_global
  integer :: i, j , k, ce, com
  DM :: DM_chemical
  PetscScalar, dimension(:,:,:,:), pointer :: mu
  PetscErrorCode :: err_PETSc
  integer(MPI_INTEGER_KIND) :: err_MPI
  class(tDict), pointer :: &
    num_grid_chemical
  real(pREAL), dimension(:,:,:,:), allocatable :: mu_0
  character(len=:), allocatable :: &
    extmsg, &
    petsc_options


  print'(/,1x,a)', '<<<+-  grid_chemical_FDM init  -+>>>'

  print'(/,1x,a)', 'P. Shanthraj et al., Computer Methods in Applied Mechanics and Engineering, 2020'
  print'(  1x,a)', 'https://doi.org/10.1016/j.cma.2020.113029'

  allocate(mu_0(homogenization_chemical_maxNcomponents-1,cells(1),cells(2),cells3),source=0.0_pREAL)
  do com = 1, homogenization_chemical_maxNcomponents - 1
    mu_0(com,:,:,:) = discretization_grid_getInitialCondition(material_name_species(com))
  end do

  N_components = homogenization_chemical_maxNcomponents - 1
!-------------------------------------------------------------------------------------------------
! read numerical parameters and do sanity checks
  num_grid_chemical => num_grid%get_dict('chemical',defaultVal=emptyDict)

  num%itmax             = num_grid_chemical%get_asInt  ('itmax',           defaultVal=250)
  num%eps_chemical_atol = num_grid_chemical%get_asReal('eps_chemical_atol',defaultVal=1.0e-6_pREAL)
  num%eps_chemical_rtol = num_grid_chemical%get_asReal('eps_chemical_rtol',defaultVal=1.0e-6_pREAL)

  if (num%itmax <= 1)                     call IO_error(301,ext_msg='itmax')
  if (num%eps_chemical_atol <= 0.0_pREAL) call IO_error(301,ext_msg='eps_chemical_atol')
  if (num%eps_chemical_rtol <= 0.0_pREAL) call IO_error(301,ext_msg='eps_chemical_rtol')

!--------------------------------------------------------------------------------------------------
! set default and user defined options for PETSc
  petsc_options = misc_prefixOptions('-snes_type newtonls -snes_ksp_ew -ksp_type fgmres -ksp_max_it 25 ' &
                                        //num_grid_chemical%get_asStr('PETSc_options',defaultVal=''), &
                                     'chemical_')
  call PetscOptionsInsertString(PETSC_NULL_OPTIONS,petsc_options,err_PETSc)
  CHKERRQ(err_PETSc)

!--------------------------------------------------------------------------------------------------
! initialize solver specific parts of PETSc
  call SNESCreate(PETSC_COMM_WORLD,SNES_chemical,err_PETSc); CHKERRQ(err_PETSc)
  call SNESSetOptionsPrefix(SNES_chemical,'chemical_',err_PETSc);CHKERRQ(err_PETSc)
  call MPI_Allgather(int(cells3,pPETSCINT),1_MPI_INTEGER_KIND,MPI_INTEGER,&
                     cells3_global,1_MPI_INTEGER_KIND,MPI_INTEGER,MPI_COMM_WORLD,err_MPI)
  call parallelization_chkerr(err_MPI)
  call DMDACreate3D(PETSC_COMM_WORLD, &
         DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC, &                        ! cut off stencil at boundary
         DMDA_STENCIL_BOX, &                                                                        ! Moore (26) neighborhood around central point
         int(cells(1),pPetscInt),int(cells(2),pPetscInt),int(cells(3),pPetscInt), &                    ! global grid
         1_pPetscInt, 1_pPetscInt, int(worldsize,pPetscInt), &
         int(N_components,pPetscInt), 1_pPetscInt, &                                                ! #dof (mu field), ghost boundary width (domain overlap)
         [int(cells(1),pPetscInt)],[int(cells(2),pPetscInt)],cells3_global, &                                ! local grid
         DM_chemical,err_PETSc)                                                                   ! handle, error
  CHKERRQ(err_PETSc)
  call SNESSetDM(SNES_chemical,DM_chemical,err_PETSc); CHKERRQ(err_PETSc)                         ! connect snes to da
  call DMsetFromOptions(DM_chemical,err_PETSc); CHKERRQ(err_PETSc)
  call DMsetUp(DM_chemical,err_PETSc); CHKERRQ(err_PETSc)
  call DMCreateGlobalVector(DM_chemical,solution_current,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMCreateGlobalVector(DM_chemical,solution_lastInc,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMSNESSetFunctionLocal(DM_chemical,formResidual,PETSC_NULL_SNES,err_PETSc)                 ! residual vector of same shape as solution vector
  CHKERRQ(err_PETSc)
  call DMSNESSetJacobianLocal(DM_chemical,formJacobian,PETSC_NULL_SNES,err_PETSc)                 ! function to evaluate stiffness matrix
  CHKERRQ(err_PETSc)
  call SNESSetMaxLinearSolveFailures(SNES_chemical, huge(1_pPetscInt), err_PETSc)
  CHKERRQ(err_PETSc)
  call SNESSetFromOptions(SNES_chemical,err_PETSc); CHKERRQ(err_PETSc)

!--------------------------------------------------------------------------------------------------
! init fields
  call VecSet(solution_current,0.0_pREAL,err_PETSc); CHKERRQ(err_PETSc)
  call VecSet(solution_lastInc,0.0_pREAL,err_PETSc); CHKERRQ(err_PETSc)

  allocate(conc_current(1:N_components+1,cells(1),cells(2),cells3),source=0.0_pREAL)
  allocate(conc_stagInc(1:N_components+1,cells(1),cells(2),cells3),source=0.0_pREAL)
  allocate(conc_lastInc(1:N_components+1,cells(1),cells(2),cells3),source=0.0_pREAL)

  call DMDAVecGetArray(DM_chemical,solution_current,mu,err_PETSc)
  CHKERRQ(err_PETSc)
  ce = 0
  do k = 1,cells3; do j = 1, cells(2); do i = 1, cells(1)
    ce = ce + 1
    conc_current(:,i,j,k) = homogenization_composition(mu_0(:,i,j,k), 1.0_pREAL, ce)
    conc_lastInc(:,i,j,k) = conc_current(:,i,j,k)
    conc_stagInc(:,i,j,k) = conc_current(:,i,j,k)
    call homogenization_chemical_setField(mu_0(:,i,j,k),conc_current(1:N_components+1,i,j,k), 1.0_pREAL, ce)
  end do; end do; end do
  mu = mu_0
  call DMDAVecRestoreArray(DM_chemical,solution_current,mu,err_PETSc)
  CHKERRQ(err_PETSc)
  call VecCopy(solution_current,solution_lastInc,err_PETSc)
  CHKERRQ(err_PETSc)

  delta = geomSize/real(cells,pREAL)

end subroutine grid_chemical_FDM_init


!--------------------------------------------------------------------------------------------------
!> @brief solution for the FEM chemical scheme with internal iterations
!--------------------------------------------------------------------------------------------------
function grid_chemical_FDM_solution(Delta_t) result(solution)

  real(pREAL), intent(in) :: Delta_t
  type(tSolutionState) :: solution

  DM :: da_local
  PetscScalar, pointer :: x_scal(:,:,:,:)
  integer   :: i, j, k, ce, com
  real(pREAL), dimension(N_components+1) :: conc_min, conc_max, conc_avg, stagNorm
  PetscErrorCode :: err_PETSc
  integer(MPI_INTEGER_KIND) :: err_MPI
  SNESConvergedReason :: reason

  solution%converged = .false.

!--------------------------------------------------------------------------------------------------
! set module wide availabe data
  Delta_t_ = Delta_t

  call SNESSolve(SNES_chemical,PETSC_NULL_VEC,solution_current,err_PETSc)
  CHKERRQ(err_PETSc)
  call SNESGetConvergedReason(SNES_chemical,reason,err_PETSc)
  CHKERRQ(err_PETSc)

#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<23)
  solution%converged = reason > SNES_CONVERGED_ITERATING ! .and. status == STATUS_OK
#else
  solution%converged = reason%v > SNES_CONVERGED_ITERATING%v ! .and. status == STATUS_OK
#endif
  solution%iterationsNeeded = merge(totalIter,num%itmax,solution%converged)

  call SNESGetDM(SNES_chemical,da_local,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMDAVecGetArray(da_local,solution_current,x_scal,err_PETSc)
  CHKERRQ(err_PETSc)

  conc_min(1:N_components+1) =  huge(1.0_pREAL)
  conc_max(1:N_components+1) = -huge(1.0_pREAL)
  conc_avg(1:N_components+1) = 0.0_pREAL
  stagNorm(1:N_components+1) = -huge(1.0_pREAL)

  ce = 0
  do k = 1+cells3offset, cells3+cells3offset; do j = 1, cells(2); do i = 1, cells(1)
    ce = ce + 1
    conc_current(:,i,j,k-cells3offset) = homogenization_composition(x_scal(0:N_components-1,i-1,j-1,k-1), Delta_t_, ce)
    do com = 1, N_components+1
      conc_min(com) = min(conc_min(com),conc_current(com,i,j,k-cells3offset))
      conc_max(com) = max(conc_max(com),conc_current(com,i,j,k-cells3offset))
      conc_avg(com) = conc_avg(com) + conc_current(com,i,j,k-cells3offset)
      stagNorm(com) = max(stagNorm(com),abs(conc_current(com,i,j,k-cells3offset) - conc_stagInc(com,i,j,k-cells3offset)))
    end do
  end do; end do; end do
  call MPI_Allreduce(MPI_IN_PLACE,stagNorm,N_components+1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD,err_MPI)
  call parallelization_chkerr(err_MPI)
  call MPI_Allreduce(MPI_IN_PLACE,conc_avg,N_components+1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,err_MPI)
  call parallelization_chkerr(err_MPI)
  call MPI_Allreduce(MPI_IN_PLACE,conc_max,N_components+1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD,err_MPI)
  call parallelization_chkerr(err_MPI)
  call MPI_Allreduce(MPI_IN_PLACE,conc_min,N_components+1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD,err_MPI)
  call parallelization_chkerr(err_MPI)
  solution%stagConverged = all(stagNorm < num%eps_chemical_atol)
  conc_stagInc = conc_current

!--------------------------------------------------------------------------------------------------
! updating chemical state
  ce = 0
  do k = 1 + cells3offset, cells3+cells3offset; do j = 1, cells(2); do i = 1, cells(1)
    ce = ce + 1
    call homogenization_chemical_setField(x_scal(0:N_components-1,i-1,j-1,k-1),&
                                           conc_current(1:N_components+1,i,j,k-cells3offset), &
                                           Delta_t_, ce)
  end do; end do; end do

  call DMDAVecRestoreArray(da_local,solution_current,x_scal,err_PETSc)
  CHKERRQ(err_PETSc)
  if (solution%converged) then
    print'(/,1x,a)', '... chemical diffusion converged ..................................'
    do com = 1, N_components+1
      print'(/,1x,a,f8.4,2x,f8.4,2x,f8.4,2x,f8.4)', 'Minimum|Maximum|Avg|Delta composition / mole_fraction = ', &
        conc_min(com), conc_max(com), conc_avg(com)*wgt, stagNorm(com)
    end do
    print'(/,1x,a)', '==========================================================================='
    flush(IO_STDOUT)
  end if


end function grid_chemical_FDM_solution


!--------------------------------------------------------------------------------------------------
!> @brief residual function
!--------------------------------------------------------------------------------------------------
subroutine formResidual(da_local,solution_current_local,residual_current_local,dummy,err_PETSc)

  DM  :: da_local
  Vec :: solution_current_local, &
         residual_current_local
  PetscObject :: dummy
  PetscErrorCode :: err_PETSc

  Vec, dimension(N_components) :: mobility_global, &
                                  mobility_local
  type(tScal),dimension(N_components) :: mobility_scal
  PetscScalar, dimension(:,:,:,:), pointer :: solution_current_scal, &
                                              residual_current_scal
  real(pREAL),dimension(3,N_components)    :: solution_grad_forward, &
                                              solution_grad_backward, &
                                              flux_forward, &
                                              flux_backward

  real(pREAL), dimension(3,N_components,N_components) :: mobility_forward, &
                                                         mobility_backward
  integer :: ce, i, j, k, dir, com_i, com_j
  real(pREAL), dimension(N_components,N_components) ::  mobility
  real(pREAL) :: norm

  call VecSet(residual_current_local,0.0_pREAL,err_PETSc); CHKERRQ(err_PETSc)
  call DMDAVecGetArray(da_local,residual_current_local,residual_current_scal,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMDAVecGetArray(da_local,solution_current_local,solution_current_scal,err_PETSc)
  CHKERRQ(err_PETSc)

  do com_i = 1, N_components
    call DMGetGlobalVector(da_local,mobility_global(com_i),err_PETSc)
    CHKERRQ(err_PETSc)
    call DMGetLocalVector(da_local,mobility_local(com_i),err_PETSc)
    CHKERRQ(err_PETSc)
    call VecSet(mobility_global(com_i),0.0_pREAL,err_PETSc)
    CHKERRQ(err_PETSc)
    call VecSet(mobility_local(com_i),0.0_pREAL,err_PETSc)
    CHKERRQ(err_PETSc)
    call DMDAVecGetArray(da_local,mobility_global(com_i),mobility_scal(com_i)%scal,err_PETSc)
    CHKERRQ(err_PETSc)
  end do

  ce = 0
  do k = 1+cells3offset, cells3offset+cells3; do j = 1, cells(2); do i = 1, cells(1)
    ce = ce + 1
    ! Get mobility tensor
    mobility = homogenization_mobility(ce)
    do com_i = 1, N_components
      mobility_scal(com_i)%scal(0:N_components-1,i-1,j-1,k-1) = mobility(com_i,1:N_components)
    end do
  end do; end do; end do

  do com_i = 1, N_components
    call DMDAVecRestoreArray(da_local,mobility_global(com_i),mobility_scal(com_i)%scal,err_PETSc)
    CHKERRQ(err_PETSc)
    call DMGlobalToLocalBegin(da_local,mobility_global(com_i),INSERT_VALUES,mobility_local(com_i),err_PETSc)
    CHKERRQ(err_PETSc)
    call DMGlobalToLocalEnd  (da_local,mobility_global(com_i),INSERT_VALUES,mobility_local(com_i),err_PETSc)
    CHKERRQ(err_PETSc)
    call DMRestoreGlobalVector(da_local,mobility_global(com_i),err_PETSc)
    CHKERRQ(err_PETSc)
    call DMDAVecGetArray(da_local,mobility_local(com_i),mobility_scal(com_i)%scal,err_PETSc)
    CHKERRQ(err_PETSc)
  end do

  ce = 0
  do k = cells3offset, cells3offset+cells3-1; do j = 0, cells(2)-1; do i = 0, cells(1)-1
    ce = ce + 1
    do com_i = 1, N_components
       mobility_forward(1,com_i,:)   =  (mobility_scal(com_i)%scal(:,i  ,j  ,k  ) &
                                     +   mobility_scal(com_i)%scal(:,i+1,j  ,k  ))/2
       mobility_forward(2,com_i,:)   =  (mobility_scal(com_i)%scal(:,i  ,j  ,k  ) &
                                     +   mobility_scal(com_i)%scal(:,i  ,j+1,k  ))/2
       mobility_forward(3,com_i,:)   =  (mobility_scal(com_i)%scal(:,i  ,j  ,k  ) &
                                     +   mobility_scal(com_i)%scal(:,i  ,j  ,k+1))/2
       mobility_backward(1,com_i,:)  =  (mobility_scal(com_i)%scal(:,i  ,j  ,k  ) &
                                     +   mobility_scal(com_i)%scal(:,i-1,j  ,k  ))/2
       mobility_backward(2,com_i,:)  =  (mobility_scal(com_i)%scal(:,i  ,j  ,k  ) &
                                     +   mobility_scal(com_i)%scal(:,i  ,j-1,k  ))/2
       mobility_backward(3,com_i,:)  =  (mobility_scal(com_i)%scal(:,i  ,j  ,k  ) &
                                     +   mobility_scal(com_i)%scal(:,i  ,j  ,k-1))/2
    end do

    ! get potential gradients
    solution_grad_forward (1,:)  = (solution_current_scal(:,i+1,j  ,k  ) &
                                 -  solution_current_scal(:,i  ,j  ,k  ))/delta(1)
    solution_grad_forward (2,:)  = (solution_current_scal(:,i  ,j+1,k  ) &
                                 -  solution_current_scal(:,i  ,j  ,k  ))/delta(2)
    solution_grad_forward (3,:)  = (solution_current_scal(:,i  ,j  ,k+1) &
                                 -  solution_current_scal(:,i  ,j  ,k  ))/delta(3)
    solution_grad_backward(1,:)  = (solution_current_scal(:,i  ,j  ,k  ) &
                                 -  solution_current_scal(:,i-1,j  ,k  ))/delta(1)
    solution_grad_backward(2,:)  = (solution_current_scal(:,i  ,j  ,k  ) &
                                 -  solution_current_scal(:,i  ,j-1,k  ))/delta(2)
    solution_grad_backward(3,:)  = (solution_current_scal(:,i  ,j  ,k  ) &
                                 -  solution_current_scal(:,i  ,j  ,k-1))/delta(3)

    ! get fluxes
    flux_forward  = 0.0_pREAL
    flux_backward = 0.0_pREAL
    do com_i = 1, N_components
      do com_j = 1, N_components
        do dir = 1, 3
          flux_forward (dir,com_i) = flux_forward (dir,com_i) &
                                   + mobility_forward(dir,com_i,com_j) &
                                   * solution_grad_forward(dir,com_j)
          flux_backward(dir,com_i) = flux_backward(dir,com_i) &
                                   + mobility_backward(dir,com_i,com_j) &
                                   * solution_grad_backward(dir,com_j)
        end do
      end do
    end do


    ! form residual (c - c_0 - deltaT * div (M*grad(mu))). source (-fdt not included yet)

    conc_current(:,i+1,j+1,k+1-cells3offset) = homogenization_composition(solution_current_scal(:,i,j,k), Delta_t_, ce)
    residual_current_scal(0:N_components-1,i,j,k) = conc_current(1:N_components,i+1,j+1,k+1-cells3offset) - &
                                      conc_lastInc(1:N_components,i+1,j+1,k+1-cells3offset)
    do com_i = 0, N_components-1
      do dir = 1, 3
        residual_current_scal(com_i,i,j,k) =  residual_current_scal(com_i,i,j,k) &
                                           -  (flux_forward (dir,com_i+1) - flux_backward(dir,com_i+1))*Delta_t_/delta(dir)
      end do
    end do
  end do; end do; end do

  do com_i = 1, N_components
    call DMDAVecRestoreArray(da_local,mobility_local(com_i),mobility_scal(com_i)%scal,err_PETSc)
    CHKERRQ(err_PETSc)
    call DMRestoreLocalVector(da_local,mobility_local(com_i),err_PETSc)
    CHKERRQ(err_PETSc)
  end do

  call DMDAVecRestoreArray(da_local,solution_current_local,solution_current_scal,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMDAVecRestoreArray(da_local,residual_current_local,residual_current_scal,err_PETSc)
  CHKERRQ(err_PETSc)
  call VecNorm(residual_current_local,NORM_INFINITY,norm,err_PETSc)
  CHKERRQ(err_PETSc)

end subroutine formResidual


subroutine formJacobian(da_local,solution_current_local,Jac,Jac_pre,dummy,err_PETSc) ! needed for periodic bc?

  DM  :: da_local
  Vec :: solution_current_local
  Mat :: Jac, Jac_pre
  PetscObject :: dummy
  PetscErrorCode :: err_PETSc

  Vec, dimension(N_components) :: mobility_global, &
                                  mobility_local
  type(tScal),dimension(N_components) :: mobility_scal
  PetscScalar, dimension(:,:,:,:), pointer :: solution_current_scal
  real(pREAL), dimension(3,N_components,N_components) :: mobility_forward, &
                                                         mobility_backward
  integer :: ce, i, j, k, dir, com, com_i
  real(pREAL), dimension(N_components,N_components) ::  mobility
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<23)
  MatStencil,dimension(4,7*N_components) :: col
  MatStencil,dimension(4,N_components) :: row
#else
  MatStencil,dimension(7*N_components) :: col
  MatStencil,dimension(N_components) :: row
#endif
  PetscScalar :: stiffness(N_components,7*N_components)

!--------------------------------------------------------------------------------------------------
! get mobilities
  do com_i = 1, N_components
    call DMGetGlobalVector(da_local,mobility_global(com_i),err_PETSc);
    CHKERRQ(err_PETSc)
    call DMGetLocalVector(da_local,mobility_local(com_i),err_PETSc)
    CHKERRQ(err_PETSc)
    call VecSet(mobility_global(com_i),0.0_pREAL,err_PETSc)
    CHKERRQ(err_PETSc)
    call VecSet(mobility_local(com_i),0.0_pREAL,err_PETSc);
    CHKERRQ(err_PETSc)
    call DMDAVecGetArray(da_local,mobility_global(com_i),mobility_scal(com_i)%scal,err_PETSc)
    CHKERRQ(err_PETSc)
  end do

  ce = 0
  do k = 1+cells3offset, cells3offset+cells3; do j = 1, cells(2); do i = 1, cells(1)
    ce = ce + 1
   ! Get mobility tensor
    mobility = homogenization_mobility(ce)
    do com_i = 1, N_components
      mobility_scal(com_i)%scal(:,i-1,j-1,k-1) = mobility(com_i,:)
    end do
  end do; end do; end do

  do com_i = 1, N_components
    call DMDAVecRestoreArray(da_local,mobility_global(com_i),mobility_scal(com_i)%scal,err_PETSc)
    CHKERRQ(err_PETSc)
    call DMGlobalToLocalBegin(da_local,mobility_global(com_i),INSERT_VALUES,mobility_local(com_i),err_PETSc)
    CHKERRQ(err_PETSc)
    call DMGlobalToLocalEnd  (da_local,mobility_global(com_i),INSERT_VALUES,mobility_local(com_i),err_PETSc)
    CHKERRQ(err_PETSc)
    call DMRestoreGlobalVector(da_local,mobility_global(com_i),err_PETSc)
    CHKERRQ(err_PETSc)
    call DMDAVecGetArray(da_local,mobility_local(com_i),mobility_scal(com_i)%scal,err_PETSc)
    CHKERRQ(err_PETSc)
  end do

!--------------------------------------------------------------------------------------------------
! assemble matrix
  call MatSetOption(Jac,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE,err_PETSc)
  CHKERRQ(err_PETSc)
  call MatSetOption(Jac,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE,err_PETSc)
  CHKERRQ(err_PETSc)
  call MatZeroEntries(Jac,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMDAVecGetArray(da_local,solution_current_local,solution_current_scal,err_PETSc)
  CHKERRQ(err_PETSc)
  ce = 0
  do k = cells3offset, cells3offset+cells3-1; do j = 0, cells(2)-1; do i = 0, cells(1)-1
    ce = ce + 1
    stiffness = 0.0_pREAL
    do com = 0, N_components-1
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<23)
      row(MatStencil_i,  com+1) = i
      row(MatStencil_j,  com+1) = j
      row(MatStencil_k,  com+1) = k
      row(MatStencil_c,  com+1) = com
      col(MatStencil_i,7*com+1) = i
      col(MatStencil_j,7*com+1) = j
      col(MatStencil_k,7*com+1) = k
      col(MatStencil_c,7*com+1) = com
      col(MatStencil_i,7*com+2) = i+1
      col(MatStencil_j,7*com+2) = j
      col(MatStencil_k,7*com+2) = k
      col(MatStencil_c,7*com+2) = com
      col(MatStencil_i,7*com+3) = i
      col(MatStencil_j,7*com+3) = j+1
      col(MatStencil_k,7*com+3) = k
      col(MatStencil_c,7*com+3) = com
      col(MatStencil_i,7*com+4) = i
      col(MatStencil_j,7*com+4) = j
      col(MatStencil_k,7*com+4) = k+1
      col(MatStencil_c,7*com+4) = com
      col(MatStencil_i,7*com+5) = i-1
      col(MatStencil_j,7*com+5) = j
      col(MatStencil_k,7*com+5) = k
      col(MatStencil_c,7*com+5) = com
      col(MatStencil_i,7*com+6) = i
      col(MatStencil_j,7*com+6) = j-1
      col(MatStencil_k,7*com+6) = k
      col(MatStencil_c,7*com+6) = com
      col(MatStencil_i,7*com+7) = i
      col(MatStencil_j,7*com+7) = j
      col(MatStencil_k,7*com+7) = k-1
      col(MatStencil_c,7*com+7) = com
#else
      row(  com+1)%i = i
      row(  com+1)%j = j
      row(  com+1)%k = k
      row(  com+1)%c = com
      col(7*com+1)%i = i
      col(7*com+1)%j = j
      col(7*com+1)%k = k
      col(7*com+1)%c = com
      col(7*com+2)%i = i+1
      col(7*com+2)%j = j
      col(7*com+2)%k = k
      col(7*com+2)%c = com
      col(7*com+3)%i = i
      col(7*com+3)%j = j+1
      col(7*com+3)%k = k
      col(7*com+3)%c = com
      col(7*com+4)%i = i
      col(7*com+4)%j = j
      col(7*com+4)%k = k+1
      col(7*com+4)%c = com
      col(7*com+5)%i = i-1
      col(7*com+5)%j = j
      col(7*com+5)%k = k
      col(7*com+5)%c = com
      col(7*com+6)%i = i
      col(7*com+6)%j = j-1
      col(7*com+6)%k = k
      col(7*com+6)%c = com
      col(7*com+7)%i = i
      col(7*com+7)%j = j
      col(7*com+7)%k = k-1
      col(7*com+7)%c = com
#endif
    end do

    stiffness(1:N_components,1:7*N_components:7) = &
             matmul(math_eye(N_components), &
                    homogenization_compositionTangent(solution_current_scal(:,i,j,k),Delta_t_,ce))

    do com_i = 1, N_components
      mobility_forward(1,com_i,:)   =  (mobility_scal(com_i)%scal(:,i  ,j  ,k  ) &
                                      + mobility_scal(com_i)%scal(:,i+1,j  ,k  ))/2
      mobility_forward(2,com_i,:)   =  (mobility_scal(com_i)%scal(:,i  ,j  ,k  ) &
                                      + mobility_scal(com_i)%scal(:,i  ,j+1,k  ))/2
      mobility_forward(3,com_i,:)   =  (mobility_scal(com_i)%scal(:,i  ,j  ,k  ) &
                                      + mobility_scal(com_i)%scal(:,i  ,j  ,k+1))/2
      mobility_backward(1,com_i,:)  =  (mobility_scal(com_i)%scal(:,i  ,j  ,k  ) &
                                      + mobility_scal(com_i)%scal(:,i-1,j  ,k  ))/2
      mobility_backward(2,com_i,:)  =  (mobility_scal(com_i)%scal(:,i  ,j  ,k  ) &
                                      + mobility_scal(com_i)%scal(:,i  ,j-1,k  ))/2
      mobility_backward(3,com_i,:)  =  (mobility_scal(com_i)%scal(:,i  ,j  ,k  ) &
                                      + mobility_scal(com_i)%scal(:,i  ,j  ,k-1))/2
    end do

    do dir = 1, 3
      stiffness(1    :N_components,1:7*N_components:7) = stiffness(1    :N_components,1:7*N_components:7) &
                                                      + (  &
                                                            mobility_forward (dir,:,:)   &
                                                          + mobility_backward(dir,:,:)  &
                                                         )* Delta_t_/delta(dir)/delta(dir)
      stiffness(1:N_components,1+dir:7*N_components:7) = stiffness(1:N_components,1+dir:7*N_components:7) &
                                                      - (  &
                                                           mobility_forward (dir,:,:)   &
                                                         )*Delta_t_/delta(dir)/delta(dir)
      stiffness(1:N_components,4+dir:7*N_components:7) = stiffness(1:N_components,4+dir:7*N_components:7) &
                                                      - (  &
                                                           mobility_backward(dir,:,:)   &
                                                         )*Delta_t_/delta(dir)/delta(dir)
    end do
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<23)
    call MatSetValuesStencil(Jac,int(N_components,pPetscInt),row,int(7*N_components,pPetscInt),&
                               col,transpose(stiffness),INSERT_VALUES,err_PETSc)
#else
    call MatSetValuesStencil(Jac,int(N_components,pPetscInt),row,int(7*N_components,pPetscInt),&
                             col,reshape(transpose(stiffness),[size(stiffness)]),INSERT_VALUES,err_PETSc)
#endif
    CHKERRQ(err_PETSc)
  end do; end do; end do
  call MatAssemblyBegin(Jac,MAT_FINAL_ASSEMBLY,err_PETSc)
  CHKERRQ(err_PETSc)
  call MatAssemblyEnd(Jac,MAT_FINAL_ASSEMBLY,err_PETSc)
  CHKERRQ(err_PETSc)
  call MatAssemblyBegin(Jac_pre,MAT_FINAL_ASSEMBLY,err_PETSc)
  CHKERRQ(err_PETSc)
  call MatAssemblyEnd(Jac_pre,MAT_FINAL_ASSEMBLY,err_PETSc)
  CHKERRQ(err_PETSc)

  do com_i = 1, N_components
    call DMDAVecRestoreArray(da_local,mobility_local(com_i),mobility_scal(com_i)%scal,err_PETSc)
    CHKERRQ(err_PETSc)
    call DMRestoreLocalVector(da_local,mobility_local(com_i),err_PETSc)
    CHKERRQ(err_PETSc)
  end do
  call DMDAVecRestoreArray(da_local,solution_current_local,solution_current_scal,err_PETSc)
  CHKERRQ(err_PETSc)

end subroutine formJacobian


!--------------------------------------------------------------------------------------------------
!> @brief forwarding routine
!--------------------------------------------------------------------------------------------------
subroutine grid_chemical_FDM_forward(cutBack)

  logical, intent(in) :: cutBack

  integer :: i, j, k, ce
  DM :: da_local
  PetscScalar, dimension(:,:,:,:), pointer :: x_scal
  PetscErrorCode :: err_PETSc

  if (cutBack) then
    conc_current = conc_lastInc
    conc_stagInc = conc_lastInc

    call SNESGetDM(SNES_chemical,da_local,err_PETSc)
    CHKERRQ(err_PETSc)
    call VecCopy(solution_lastInc,solution_current,err_PETSc)
    CHKERRQ(err_PETSc)
    call DMDAVecGetArray(da_local,solution_current,x_scal,err_PETSc)
    CHKERRQ(err_PETSc)
    ce = 0
    do k = 1, cells3offset, cells3+cells3offset; do j = 1, cells(2); do i = 1, cells(1)
      ce = ce + 1
      call homogenization_chemical_setField(x_scal(0:N_components-1,i-1,j-1,k-1),&
                                           conc_current(1:N_components+1,i,j,k-cells3offset), &
                                           Delta_t_, ce)
    end do; end do; end do
    call DMDAVecRestoreArray(da_local,solution_current,x_scal,err_PETSc)
    CHKERRQ(err_PETSc)
  else
    conc_lastInc = conc_current
    call VecCopy(solution_current,solution_lastInc,err_PETSc)
    CHKERRQ(err_PETSc)
  end if


end subroutine grid_chemical_FDM_forward


end module grid_chemical_FDM
