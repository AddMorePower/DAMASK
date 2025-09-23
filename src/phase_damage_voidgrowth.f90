!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, KU Leuven
!> @brief material subroutine incoprorating isotropic void growth
!> @details to be done
!--------------------------------------------------------------------------------------------------
submodule(phase:damage) voidgrowth

  type :: tParameters                                                                               !< container type for internal constitutive parameters
    real(pREAL) :: &
      W_crit                                                                                        !< critical elastic strain energy
    character(len=pSTRLEN), allocatable, dimension(:) :: &
      output
  end type tParameters

  type :: tVoidgrowthState
    real(pREAL), pointer, dimension(:) :: &                                                         !< vectors along Nmembers
      r_W                                                                                           !< ratio between actual and critical strain energy density
  end type tVoidgrowthState

  type(tParameters),      allocatable, dimension(:) :: param                                        !< containers of constitutive parameters (len Ninstances)
  type(tVoidgrowthState), allocatable, dimension(:) :: &
    deltaState, &
    state

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
module function voidgrowth_init() result(mySources)

  logical, dimension(:), allocatable :: mySources

  type(tDict), pointer :: &
    phases, &
    phase, &
    src
  integer :: Nmembers,ph
  character(len=:), allocatable :: &
    refs, &
    extmsg


  mySources = source_active('voidgrowth')
  if (count(mySources) == 0) return

  print'(/,1x,a)', '<<<+-  phase:damage:voidgrowth init  -+>>>'
  print'(/,1x,a,1x,i0)', '# phases:',count(mySources); flush(IO_STDOUT)


  phases => config_material%get_dict('phase')
  allocate(param(size(phases)))
  allocate(state(size(phases)))
  allocate(deltaState(size(phases)))
  extmsg = ''

  do ph = 1, size(phases)
    if (mySources(ph)) then
      phase => phases%get_dict(ph)
      src => phase%get_dict('damage')

      associate(prm => param(ph), dlt => deltaState(ph), stt => state(ph))

        prm%W_crit = src%get_asReal('G_crit')/src%get_asReal('l_c')

        print'(/,1x,a,1x,i0,a)', 'phase',ph,': '//phases%key(ph)
        refs = config_listReferences(src,indent=3)
        if (len(refs) > 0) print'(/,1x,a)', refs

#if defined (__GFORTRAN__)
        prm%output = output_as1dStr(src)
#else
        prm%output = src%get_as1dStr('output',defaultVal=emptyStrArray)
#endif

        ! sanity checks
        if (prm%W_crit <= 0.0_pREAL) extmsg = trim(extmsg)//' W_crit'

        Nmembers = count(material_ID_phase==ph)
        call phase_allocateState(damageState(ph),Nmembers,1,0,1)
        damageState(ph)%atol = src%get_asReal('atol_phi',defaultVal=1.0e-9_pREAL)
        if (any(damageState(ph)%atol < 0.0_pREAL)) extmsg = trim(extmsg)//' atol_phi'

        stt%r_W => damageState(ph)%state(1,:)
        dlt%r_W => damageState(ph)%deltaState(1,:)

      end associate


      if (extmsg /= '') call IO_error(211,ext_msg=trim(extmsg)//'(damage_voidgrowth)')
    end if

  end do

end function voidgrowth_init


!--------------------------------------------------------------------------------------------------
!> @brief
!--------------------------------------------------------------------------------------------------
module subroutine voidgrowth_deltaState(C, Fe, ph,en)

  integer, intent(in) :: ph,en
  real(pREAL),  intent(in), dimension(3,3) :: &
    Fe
  real(pREAL),  intent(in), dimension(6,6) :: &
    C

  real(pREAL), dimension(6) :: &
    epsilon_e
  real(pREAL) :: &
    r_W


  epsilon_e = math_33toVoigt6_strain(0.5_pREAL*(matmul(transpose(Fe),Fe)-math_I3))

  associate(prm => param(ph), stt => state(ph), dlt => deltaState(ph))

    r_W = math_trace33(math_Voigt6to33_stress(matmul(C,epsilon_e)))/3.0_pREAL/prm%W_crit            ! math_spherical33 maps to tensor
    dlt%r_W(en) = merge(r_W - stt%r_W(en), 0.0_pREAL, r_W > stt%r_W(en))

  end associate

end subroutine voidgrowth_deltaState


!--------------------------------------------------------------------------------------------------
!> @brief Write results to HDF5 output file.
!--------------------------------------------------------------------------------------------------
module subroutine voidgrowth_result(phase,group)

  integer,          intent(in) :: phase
  character(len=*), intent(in) :: group

  integer :: o


  associate(prm => param(phase), stt => damageState(phase)%state) ! point to state and output r_W (is scalar, not 1D vector)

    outputsLoop: do o = 1,size(prm%output)
      select case(trim(prm%output(o)))
        case ('r_W')
          call result_writeDataset(stt,group,trim(prm%output(o)),'ratio between actual and critical strain energy density','-')
      end select
    end do outputsLoop

  end associate

end subroutine voidgrowth_result

end submodule voidgrowth
