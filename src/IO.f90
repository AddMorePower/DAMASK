!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Christoph Kords, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief  input/output functions
!--------------------------------------------------------------------------------------------------
module IO
  use, intrinsic :: ISO_fortran_env, only: &
    IO_STDOUT => OUTPUT_UNIT, &
    IO_STDERR => ERROR_UNIT

  use prec
  use constants
  use misc
#ifndef MARC_SOURCE
  use OS
#endif

implicit none(type,external)
  private

  ! For transition period
  interface IO_error
    module procedure IO_error_new
    module procedure IO_error_old
  end interface IO_error
  interface IO_warning
    module procedure IO_warning_new
    module procedure IO_warning_old
  end interface IO_warning

  character, parameter, public :: &
    IO_ESC = achar(27), &                                                                           !< escape character
    IO_EOL = LF                                                                                     !< end of line character
  character(len=*), parameter, public :: &
    IO_FORMATRESET = IO_ESC//'[0m', &                                                               !< reset formatting
    IO_EMPH = IO_ESC//'[3m', &                                                                      !< emphasize (italics)
    IO_QUOTES  = "'"//'"' , &                                                                       !< quotes for strings
    IO_WHITESPACE = achar(44)//achar(32)//achar(9)//achar(10)//achar(13)                            !< whitespace characters

  public :: &
    IO_init, &
    IO_selfTest, &
    IO_read, &
    IO_wrapLines, &
    IO_lc, &
    IO_glueDiffering, &
    IO_intAsStr, &
    IO_realAsStr, &
    IO_strAsInt, &
    IO_strAsReal, &
    IO_strAsBool, &
    IO_color, &
    IO_error, &
    IO_warning, &
    IO_STDOUT, &
    tokenize

contains


!--------------------------------------------------------------------------------------------------
!> @brief Do self test.
!--------------------------------------------------------------------------------------------------
subroutine IO_init()

  print'(/,1x,a)', '<<<+-  IO init  -+>>>'; flush(IO_STDOUT)

  call IO_selfTest()

end subroutine IO_init


!--------------------------------------------------------------------------------------------------
!> @brief Read ASCII file.
!> @details Proper Unix style (LF line endings and LF at EOF) is ensured.
!--------------------------------------------------------------------------------------------------
function IO_read(fileName) result(fileContent)

  character(len=*),  intent(in) :: fileName
  character(len=:), allocatable :: fileContent

  integer ::  &
    fileUnit, &
    myStat
  integer(pI64) ::  &
    fileLength


  inquire(file = fileName, size=fileLength)
  open(newunit=fileUnit, file=fileName, access='stream',&
       status='old', position='rewind', action='read',iostat=myStat)
  if (myStat /= 0) call IO_error(100_pI16,trim(fileName),emph=[1])
  allocate(character(len=fileLength)::fileContent)
  if (fileLength==0) then
    close(fileUnit)
    return
  end if

  read(fileUnit,iostat=myStat) fileContent
  if (myStat /= 0) call IO_error(102_pI16,trim(fileName),emph=[1])
  close(fileUnit)

  if (index(fileContent,CR//LF,kind=pI64) /= 0)     fileContent = CRLF2LF(fileContent)
  if (fileContent(fileLength:fileLength) /= IO_EOL) fileContent = fileContent//IO_EOL               ! ensure EOL@EOF

end function IO_read


!--------------------------------------------------------------------------------------------------
!> @brief Insert EOL at separator trying to keep line length below limit.
!--------------------------------------------------------------------------------------------------
function IO_wrapLines(str,separator,filler,length)

  character(len=*),    intent(in) :: str                                                            !< string to split
  character, optional, intent(in) :: separator                                                      !< line breaks are possible after this character, defaults to ','
  character(len=*), optional, intent(in) :: filler                                                  !< character(s) to insert after line break, defaults to none
  integer,   optional, intent(in) :: length                                                         !< (soft) line limit, defaults to 80
  character(len=:), allocatable :: IO_wrapLines

  integer, dimension(:), allocatable :: pos_sep, pos_split
  integer :: i,s,e


  i = index(str,misc_optional(separator,','))
  if (i == 0) then
    IO_wrapLines = str
  else
    pos_sep = [0]
    s = i
    do while (i /= 0 .and. s < len(str))
      pos_sep = [pos_sep,s]
      i = index(str(s+1:),misc_optional(separator,','))
      s = s + i
    end do
    pos_sep = [pos_sep,len(str)]

    pos_split = emptyIntArray
    s = 1
    e = 2
    IO_wrapLines = ''
    do while (e < size(pos_sep))
      if (pos_sep(e+1) - pos_sep(s) >= misc_optional(length,80)) then
        IO_wrapLines = IO_wrapLines//adjustl(str(pos_sep(s)+1:pos_sep(e)))//IO_EOL//misc_optional(filler,'')
        s = e
      end if
      e = e + 1
    end do
    IO_wrapLines = IO_wrapLines//adjustl(str(pos_sep(s)+1:))
  end if

end function IO_wrapLines


!--------------------------------------------------------------------------------------------------
!> @brief Convert characters in string to lower case.
!--------------------------------------------------------------------------------------------------
pure function IO_lc(str)

  character(len=*), intent(in) :: str                                                               !< string to convert
  character(len=len(str))   :: IO_lc

  integer :: i,n


  do i = 1,len(str)
    n = index(UPPER,str(i:i))
    if (n==0) then
      IO_lc(i:i) = str(i:i)
    else
      IO_lc(i:i) = LOWER(n:n)
    end if
  end do

end function IO_lc


!--------------------------------------------------------------------------------------------------
! @brief Return first (with glued on second if they differ).
!--------------------------------------------------------------------------------------------------
pure function IO_glueDiffering(first,second,glue)

  character(len=*),           intent(in)  :: first
  character(len=*),           intent(in)  :: second
  character(len=*), optional, intent(in)  :: glue
  character(len=:), allocatable :: IO_glueDiffering

  character(len=:), allocatable           :: glue_


  glue_ = misc_optional(glue,'<--')
  IO_glueDiffering = trim(first)
  if (trim(first) /= trim(second)) IO_glueDiffering = IO_glueDiffering//' '//trim(glue_)//' '//trim(second)

end function IO_glueDiffering


!--------------------------------------------------------------------------------------------------
!> @brief Return given int value as string.
!--------------------------------------------------------------------------------------------------
pure function IO_intAsStr(i)

  integer, intent(in)            :: i
  character(len=:), allocatable  :: IO_intAsStr


  allocate(character(len=merge(2,1,i<0) + floor(log10(real(abs(merge(1,i,i==0))))))::IO_intAsStr)
  write(IO_intAsStr,'(i0)') i

end function IO_intAsStr


!--------------------------------------------------------------------------------------------------
!> @brief Return given float value as string.
!--------------------------------------------------------------------------------------------------
pure function IO_realAsStr(f)

  real(pREAL), intent(in)        :: f
  character(len=:), allocatable  :: IO_realAsStr
  character(len=15)              :: tmp


  write(tmp,'(g15.7)') f
  tmp = adjustl(tmp)
  allocate(IO_realAsStr,source=tmp(:len_trim(tmp)))

end function IO_realAsStr


!--------------------------------------------------------------------------------------------------
!> @brief Return integer value from given string.
!--------------------------------------------------------------------------------------------------
integer function IO_strAsInt(str)

  character(len=*), intent(in) :: str                                                               !< string for conversion to int value

  integer :: readStatus


  read(str,*,iostat=readStatus) IO_strAsInt
  if (readStatus /= 0) call IO_error(111_pI16,'cannot represent',str,'as integer',emph=[2])

end function IO_strAsInt


!--------------------------------------------------------------------------------------------------
!> @brief Return real value from given string.
!--------------------------------------------------------------------------------------------------
real(pREAL) function IO_strAsReal(str)

  character(len=*), intent(in) :: str                                                               !< string for conversion to real value

  integer :: readStatus


  read(str,*,iostat=readStatus) IO_strAsReal
  if (readStatus /= 0) call IO_error(111_pI16,'cannot represent',str,'as real',emph=[2])

end function IO_strAsReal


!--------------------------------------------------------------------------------------------------
!> @brief Return logical value from given string.
!> @details: 'True' and 'true' are converted to .true.
!> @details: 'False' and 'false' are converted to .false.
!--------------------------------------------------------------------------------------------------
logical function IO_strAsBool(str)

  character(len=*), intent(in) :: str                                                               !< string for conversion to boolean


  if     (trim(adjustl(str)) == 'True' .or.  trim(adjustl(str)) == 'true') then
    IO_strAsBool = .true.
  elseif (trim(adjustl(str)) == 'False' .or. trim(adjustl(str)) == 'false') then
    IO_strAsBool = .false.
  else
    call IO_error(111_pI16,'cannot represent',str,'as boolean',emph=[2])
  end if

end function IO_strAsBool


!--------------------------------------------------------------------------------------------------
!> @brief Return string to set foreground and/or background color.
!> @details Only active if unit is a TTY. Does nothing for MSC.Marc. No color disables formatting.
!> @details https://stackoverflow.com/questions/4842424
!--------------------------------------------------------------------------------------------------
function IO_color(fg,bg,unit)

  character(len=:), allocatable :: IO_color
  integer, intent(in), dimension(3), optional :: &
    fg, &                                                                                           !< foreground color (8 bit RGB)
    bg                                                                                              !< background color (8 bit RGB)
  integer, intent(in), optional :: unit                                                             !< output unit (default STDOUT)


  IO_color = ''

#ifndef MARC_SOURCE
  if (.not. OS_isaTTY(misc_optional(unit,IO_STDOUT))) return

  if (present(fg)) &
    IO_color = IO_color//IO_ESC//'[38;2;'//IO_intAsStr(fg(1))//';' &
                                         //IO_intAsStr(fg(2))//';' &
                                         //IO_intAsStr(fg(3))//'m'
  if (present(bg)) &
    IO_color = IO_color//IO_ESC//'[48;2;'//IO_intAsStr(bg(1))//';' &
                                         //IO_intAsStr(bg(2))//';' &
                                         //IO_intAsStr(bg(3))//'m'

  if (.not. present(fg) .and. .not. present(bg)) IO_color = IO_FORMATRESET
#endif

end function IO_color


!--------------------------------------------------------------------------------------------------
!> @brief Write error statements and terminate the run with exit #9xxx.
!> @details Should become "IO_error" after completed migration.
!--------------------------------------------------------------------------------------------------
subroutine IO_error_new(error_ID, &
                        info_1,info_2,info_3,info_4,info_5,info_6,info_7,info_8,info_9, &
                        emph)


  integer(pI16),      intent(in) :: error_ID        ! should go back to default integer after completed migration.
  class(*), optional, intent(in) :: info_1,info_2,info_3,info_4,info_5,info_6,info_7,info_8,info_9
  integer, dimension(:), optional, intent(in) :: emph                                               !< which info(s) to emphasize

  character(len=:), allocatable :: msg
  external :: quit


  select case (error_ID)

!--------------------------------------------------------------------------------------------------
! file handling errors
    case (100)
      msg = 'could not open file:'
    case (102)
      msg = 'could not read file:'

!--------------------------------------------------------------------------------------------------
! file parsing errors
    case (110)
      msg = 'invalid chunk selected'
    case (111)
      msg = 'invalid string for conversion'
    case (114)
      msg = 'cannot decode base64 string:'

!--------------------------------------------------------------------------------------------------
! lattice error messages
    case (130)
      msg = 'invalid crystal parameters'
    case (132)
      msg = 'invalid parameters for transformation'
    case (137)
      msg = 'not defined for lattice structure'
    case (138)
      msg = 'not enough interaction parameters given'

!--------------------------------------------------------------------------------------------------
! errors related to the parsing of material.yaml
    case (140)
      msg = 'key not found'
    case (141)
      msg = 'number of chunks in string differs'
    case (142)
      msg = 'empty list'
    case (143)
      msg = 'no value found for key'
    case (144)
      msg = 'negative number systems requested'
    case (145)
      msg = 'too many systems requested'
    case (146)
      msg = 'number of values does not match'
    case (147)
      msg = 'V_e needs to be symmetric'
    case (148)
      msg = 'Nconstituents mismatch between homogenization and material'

!--------------------------------------------------------------------------------------------------
! material error messages and related messages in geometry
    case (150)
      msg = 'index out of bounds'
    case (153)
      msg = 'sum of phase fractions differs from 1'
    case (155)
      msg = 'material index out of bounds'
    case (180)
      msg = 'missing/invalid material definition'
    case (190)
      msg = 'unknown element type:'
    case (191)
      msg = 'mesh contains more than one element type'

!--------------------------------------------------------------------------------------------------
! plasticity error messages
    case (200)
      msg = 'unknown type specified:'

    case (211)
      msg = 'material parameter out of bounds:'
    case (212)
      msg = 'nonlocal model not supported'

!--------------------------------------------------------------------------------------------------
! numerics error messages
    case (301)
      msg = 'numerics parameter out of bounds:'

!--------------------------------------------------------------------------------------------------
! math errors
    case (402)
      msg = 'invalid orientation specified'

!-------------------------------------------------------------------------------------------------
! homogenization errors
    case (500)
      msg = 'unknown homogenization specified'
    case (501)
      msg = 'homogenization description absent'

!--------------------------------------------------------------------------------------------------
! user errors
    case (600)
      msg = 'only one source entry allowed'
    case (603)
      msg = 'invalid data for table'
    case (610)
      msg = 'missing value for command line flag'
    case (611)
      msg = 'invalid value for command line flag'
    case (612)
      msg = 'missing command line flag'
    case (613)
      msg = 'invalid command line flag'
    case (640)
      msg = 'invalid working directory'


!------------------------------------------------------------------------------------------------
! errors related to YAML data
    case (701)
      msg = 'incorrect indent/Null value not allowed'
    case (702)
      msg = 'invalid use of flow YAML'
    case (703)
      msg = 'invalid YAML'
    case (704)
      msg = 'space expected after a colon for <key>: <value> pair'
    case (705)
      msg = 'unsupported feature'
    case (706)
      msg = 'type mismatch in YAML data node'
    case (707)
      msg = 'abrupt end of file'
    case (708)
      msg = '"---" expected after YAML file header'
    case (709)
      msg = 'length mismatch'
    case (710)
      msg = 'closing quotation mark missing in string'

!-------------------------------------------------------------------------------------------------
! errors related to the mesh solver
    case (821)
      msg = 'order not supported'

!-------------------------------------------------------------------------------------------------
! errors related to the grid solver
    case (831)
      msg = 'mask consistency violated in grid load case'
    case (833)
      msg = 'non-positive ratio for geometric progression'
    case (834)
      msg = 'negative time increment in grid load case'
    case (835)
      msg = 'non-positive increments in grid load case'
    case (836)
      msg = 'non-positive result frequency in grid load case'
    case (837)
      msg = 'incomplete loadcase'
    case (838)
      msg = 'mixed boundary conditions allow rotation'
    case (839)
      msg = 'non-positive restart frequency in grid load case'
    case (844)
      msg = 'invalid VTI file'
    case (891)
      msg = 'unknown solver type selected'
    case (892)
      msg = 'unknown filter type selected'
    case (894)
      msg = 'MPI error'

    case (950)
      msg = 'max number of cutbacks exceeded, terminating'

    case default
      error stop 'invalid error number'

  end select

  call panel('error',int(error_ID),msg, &
             info_1,info_2,info_3,info_4,info_5,info_6,info_7,info_8,info_9, &
             emph)
  call quit(9000+int(error_ID))

end subroutine IO_error_new


!--------------------------------------------------------------------------------------------------
!> @brief Write error statements and terminate the run with exit #9xxx.
!> @details Deprecated.
!--------------------------------------------------------------------------------------------------
subroutine IO_error_old(error_ID,ext_msg,label1,ID1,label2,ID2)

  integer,                    intent(in) :: error_ID
  character(len=*), optional, intent(in) :: ext_msg,label1,label2
  integer,          optional, intent(in) :: ID1,ID2

  external                      :: quit
  character(len=:), allocatable :: msg_extra


  if (.not. present(label1) .and. present(ID1)) error stop 'missing label for value 1'
  if (.not. present(label2) .and. present(ID2)) error stop 'missing label for value 2'

  msg_extra = ''
  if (present(ext_msg)) msg_extra = msg_extra//ext_msg//IO_EOL
  if (present(label1)) then
    msg_extra = msg_extra//'at '//label1
    if (present(ID1)) msg_extra = msg_extra//' '//IO_intAsStr(ID1)
    msg_extra = msg_extra//IO_EOL
  end if
  if (present(label2)) then
    msg_extra = msg_extra//'at '//label2
    if (present(ID2)) msg_extra = msg_extra//' '//IO_intAsStr(ID2)
    msg_extra = msg_extra//IO_EOL
  end if

  call IO_error_new(int(error_ID,pI16),msg_extra,IO_EOL)

end subroutine IO_error_old


!--------------------------------------------------------------------------------------------------
!> @brief Write warning statements.
!> @details Should become "IO_warning" after completed migration.
!--------------------------------------------------------------------------------------------------
subroutine IO_warning_new(warning_ID, &
                          info_1,info_2,info_3,info_4,info_5,info_6,info_7,info_8,info_9, &
                          emph)


  integer(pI16),      intent(in) :: warning_ID        ! should go back to default integer after completed migration.
  class(*), optional, intent(in) :: info_1,info_2,info_3,info_4,info_5,info_6,info_7,info_8,info_9
  integer, dimension(:), optional, intent(in) :: emph                                              !< which info(s) to emphasize

  character(len=:), allocatable :: msg

  select case (warning_ID)
    case (10)
      msg = 'deprecated keyword'
    case (47)
      msg = 'invalid parameter for FFTW'
    case (207)
      msg = 'line truncated'
    case (600)
      msg = 'crystallite responds elastically'
    case (601)
      msg = 'stiffness close to zero'
    case (709)
      msg = 'read only the first document'

    case default
      error stop 'invalid warning number'
  end select

  call panel('warning',int(warning_ID),msg, &
             info_1,info_2,info_3,info_4,info_5,info_6,info_7,info_8,info_9, &
             emph)

end subroutine IO_warning_new


!--------------------------------------------------------------------------------------------------
!> @brief Write warning statements.
!--------------------------------------------------------------------------------------------------
subroutine IO_warning_old(warning_ID,ext_msg,label1,ID1,label2,ID2)

  integer,                    intent(in) :: warning_ID
  character(len=*), optional, intent(in) :: ext_msg,label1,label2
  integer,          optional, intent(in) :: ID1,ID2

  character(len=:), allocatable :: msg,msg_extra


  if (.not. present(label1) .and. present(ID1)) error stop 'missing label for value 1'
  if (.not. present(label2) .and. present(ID2)) error stop 'missing label for value 2'

  select case (warning_ID)
    case (47)
      msg = 'invalid parameter for FFTW'
    case (207)
      msg = 'line truncated'
    case (600)
      msg = 'crystallite responds elastically'
    case (601)
      msg = 'stiffness close to zero'
    case (709)
      msg = 'read only the first document'

    case default
      error stop 'invalid warning number'
  end select

  msg_extra = ''
  if (present(ext_msg)) msg_extra = msg_extra//ext_msg//IO_EOL
  if (present(label1)) then
    msg_extra = msg_extra//'at '//label1
    if (present(ID1)) msg_extra = msg_extra//' '//IO_intAsStr(ID1)
    msg_extra = msg_extra//IO_EOL
  end if
  if (present(label2)) then
    msg_extra = msg_extra//'at '//label2
    if (present(ID2)) msg_extra = msg_extra//' '//IO_intAsStr(ID2)
    msg_extra = msg_extra//IO_EOL
  end if
  call panel('warning',warning_ID,msg,msg_extra,IO_EOL)

end subroutine IO_warning_old


!--------------------------------------------------------------------------------------------------
!> @brief Convert Windows (CRLF) to Unix (LF) line endings.
!--------------------------------------------------------------------------------------------------
pure function CRLF2LF(str)

  character(len=*), intent(in)  :: str
  character(len=:), allocatable :: CRLF2LF

  integer(pI64) :: c,n


  allocate(character(len=len_trim(str,pI64))::CRLF2LF)
  if (len(CRLF2LF,pI64) == 0) return

  n = 0_pI64
  do c=1_pI64, len_trim(str,pI64)
    CRLF2LF(c-n:c-n) = str(c:c)
    if (c == len_trim(str,pI64)) exit
    if (str(c:c+1_pI64) == CR//LF) n = n + 1_pI64
  end do

  CRLF2LF = CRLF2LF(:c-n)

end function CRLF2LF

#if ((defined(__INTEL_COMPILER) && __INTEL_COMPILER_BUILD_DATE < 20240000) || !defined(__INTEL_COMPILER))
!--------------------------------------------------------------------------------------------------
!> @brief Fortran 2023 tokenize (first form).
!--------------------------------------------------------------------------------------------------
pure subroutine tokenize(string,set,tokens)

  character(len=*), intent(in) :: string, set
  character(len=:), dimension(:), allocatable, intent(out) :: tokens

  integer, allocatable, dimension(:,:) :: pos
  integer :: i, s, e


  allocate(pos(2,0))
  e = 0
  do while (e < verify(string,set,back=.true.))
    s = e + merge(verify(string(e+1:),set),1,scan(string(e+1:),set)/=0)
    e = s + merge(scan(string(s:),set)-2,len(string(s:))-1,scan(string(s:),set)/=0)
    pos = reshape([pos,[s,e]],[2,size(pos)/2+1])
  end do
  allocate(character(len=merge(maxval(pos(2,:)-pos(1,:))+1,0,size(pos)>0))::tokens(size(pos,2)))
  do i = 1, size(pos,2)
    tokens(i) = string(pos(1,i):pos(2,i))
  end do

end subroutine tokenize
#endif

!--------------------------------------------------------------------------------------------------
!> @brief Write statements to standard error.
!--------------------------------------------------------------------------------------------------
subroutine panel(paneltype,ID,msg, &
                 info_1,info_2,info_3,info_4,info_5,info_6,info_7,info_8,info_9, &
                 emph)

  character(len=*),           intent(in) :: paneltype, msg
  integer,                    intent(in) :: ID
  class(*),         optional, intent(in) :: info_1,info_2,info_3,info_4,info_5,info_6,info_7,info_8,info_9
  integer, dimension(:), optional, intent(in) :: emph                                               !< which info(s) to emphasize

  integer, parameter :: panelwidth = 69
  character(len=*), parameter :: DIVIDER = repeat('─',panelwidth)
  character(len=pSTRLEN) :: formatString
  character(len=:), allocatable :: heading, msg_, info_extra
  character(len=:), dimension(:), allocatable :: info_split
  integer :: len_corrected, &                                                                       !< string length corrected for control characters
             i


  heading = paneltype//' '//IO_intAsStr(ID)

  select case (paneltype)

    case ('error')
       msg_ = IO_color([255,0,0],  unit=IO_STDERR)//trim(msg)//IO_color(unit=IO_STDERR)
    case ('warning')
       msg_ = IO_color([255,192,0],unit=IO_STDERR)//trim(msg)//IO_color(unit=IO_STDERR)
    case default
       error stop 'invalid panel type: '//trim(paneltype)

  end select

  info_extra = as_str(info_1,any(emph==1)) &
            // as_str(info_2,any(emph==2)) &
            // as_str(info_3,any(emph==3)) &
            // as_str(info_4,any(emph==4)) &
            // as_str(info_5,any(emph==5)) &
            // as_str(info_6,any(emph==6)) &
            // as_str(info_7,any(emph==7)) &
            // as_str(info_8,any(emph==8)) &
            // as_str(info_9,any(emph==9))


  !$OMP CRITICAL (write2out)
  write(IO_STDERR,'(/,a)')                ' ┌'       //DIVIDER//        '┐'
  write(formatString,'(a,i2,a)') '(a,24x,a,',max(1,panelwidth-24-len_trim(heading)),'x,a)'
  write(IO_STDERR,formatString)           ' │',    trim(heading),       '│'
  write(IO_STDERR,'(a)')                  ' ├'       //DIVIDER//        '┤'
  write(formatString,'(a,i3.3,a,i3.3,a)') '(a,a',max(1,len_trim(msg_)),',',&
                                                 max(1,panelwidth+3-len_trim(msg)-4),'x,a)'
  write(IO_STDERR,formatString)           ' │ ',     trim(msg_),        '│'
  if (len_trim(info_extra) > 0) then
    call tokenize(info_extra,IO_EOL,info_split)
    do i = 1, size(info_split)
      info_extra = adjustl(info_split(i))
      if (len_trim(info_extra) == 0) then
        write(IO_STDERR,'(a)')            ' │'//repeat(' ',panelwidth)//'│'
      else
        len_corrected = len_trim(info_extra) - count([(info_extra(i:i)==IO_ESC,i=1,len_trim(info_extra))])*4
        write(formatString,'(a,i3.3,a,i3.3,a)') '(a,a',max(1,len_trim(info_extra)),',',&
                                                       max(1,panelwidth+3-len_corrected-4),'x,a)'
        write(IO_STDERR,formatString)     ' │ ',    trim(info_extra),   '│'
      end if
    end do
  endif
  write(IO_STDERR,'(a)')                  ' └'       //DIVIDER//        '┘'
  flush(IO_STDERR)
  !$OMP END CRITICAL (write2out)

  contains

    !-----------------------------------------------------------------------------------------------
    !> @brief Convert to string with white space prefix and optional emphasis.
    !-----------------------------------------------------------------------------------------------
    function as_str(info,emph)

      character(len=:), allocatable :: as_str
      class(*), optional, intent(in) :: info
      logical, intent(in) :: emph


      if (present(info)) then
        select type(info)
          type is (character(*))
            as_str = info
          type is (integer)
            as_str = IO_intAsStr(info)
          type is (real(pREAL))
            as_str = IO_realAsStr(info)
          class default
            error stop 'cannot convert info argument to string'
        end select

        if (emph) then
#ifndef MARC_SOURCE
          if (OS_isaTTY(IO_STDERR)) then
            as_str = IO_EMPH//as_str//IO_FORMATRESET
          else
            as_str = IO_QUOTES(2:2)//as_str//IO_QUOTES(2:2)
          end if
#else
          as_str = IO_QUOTES(2:2)//as_str//IO_QUOTES(2:2)
#endif
        end if
        as_str = ' '//as_str
      else
        as_str = ''
      end if

    end function as_str

end subroutine panel


!--------------------------------------------------------------------------------------------------
!> @brief Check correctness of some IO functions.
!--------------------------------------------------------------------------------------------------
subroutine IO_selfTest()

  character(len=:),      allocatable :: str
  character(len=:), dimension(:), allocatable :: tokens


  if (dNeq(1.0_pREAL, IO_strAsReal('1.0')))          error stop 'IO_strAsReal'
  if (dNeq(1.0_pREAL, IO_strAsReal('1e0')))          error stop 'IO_strAsReal'
  if (dNeq(0.1_pREAL, IO_strAsReal('1e-1')))         error stop 'IO_strAsReal'
  if (dNeq(0.1_pREAL, IO_strAsReal('1.0e-1')))       error stop 'IO_strAsReal'
  if (dNeq(0.1_pREAL, IO_strAsReal('1.00e-1')))      error stop 'IO_strAsReal'
  if (dNeq(10._pREAL, IO_strAsReal(' 1.0e+1 ')))     error stop 'IO_strAsReal'

  if (3112019  /= IO_strAsInt( '3112019'))           error stop 'IO_strAsInt'
  if (3112019  /= IO_strAsInt(' 3112019'))           error stop 'IO_strAsInt'
  if (-3112019 /= IO_strAsInt('-3112019'))           error stop 'IO_strAsInt'
  if (3112019  /= IO_strAsInt('+3112019 '))          error stop 'IO_strAsInt'
  if (3112019  /= IO_strAsInt('03112019 '))          error stop 'IO_strAsInt'
  if (3112019  /= IO_strAsInt('+03112019'))          error stop 'IO_strAsInt'

  if (.not. IO_strAsBool(' true'))                   error stop 'IO_strAsBool'
  if (.not. IO_strAsBool(' True '))                  error stop 'IO_strAsBool'
  if (      IO_strAsBool(' false'))                  error stop 'IO_strAsBool'
  if (      IO_strAsBool('False'))                   error stop 'IO_strAsBool'

  if ('1234' /= IO_intAsStr(1234))                   error stop 'IO_intAsStr'
  if ('-12'  /= IO_intAsStr(-0012))                  error stop 'IO_intAsStr'

  if ('-0.1200000' /= IO_realAsStr(-0.12_pREAL))        error stop 'IO_realAsStr'
  if ('0.1234000E-31' /= IO_realAsStr(123.4e-34_pREAL)) error stop 'IO_realAsStr'

  if (CRLF2LF('') /= '')                             error stop 'CRLF2LF/0'
  if (CRLF2LF(LF)     /= LF)                         error stop 'CRLF2LF/1a'
  if (CRLF2LF(CR//LF) /= LF)                         error stop 'CRLF2LF/1b'
  if (CRLF2LF(' '//LF)     /= ' '//LF)               error stop 'CRLF2LF/2a'
  if (CRLF2LF(' '//CR//LF) /= ' '//LF)               error stop 'CRLF2LF/2b'
  if (CRLF2LF('A'//CR//LF//'B') /= 'A'//LF//'B')     error stop 'CRLF2LF/3'
  if (CRLF2LF('A'//CR//LF//'B'//CR//LF) /= &
              'A'//LF//'B'//LF)                      error stop 'CRLF2LF/4'
  if (CRLF2LF('A'//LF//CR//'B') /= 'A'//LF//CR//'B') error stop 'CRLF2LF/5'

  str='*(HiU!)3';if ('*(hiu!)3' /= IO_lc(str))       error stop 'IO_lc'

  if ('abc, def' /= IO_wrapLines('abc, def')) &
                                                     error stop 'IO_wrapLines/1'
  if ('abc,'//IO_EOL//'def' /= IO_wrapLines('abc,def',length=3)) &
                                                     error stop 'IO_wrapLines/2'
  if ('abc,'//IO_EOL//'def' /= IO_wrapLines('abc,def',length=5)) &
                                                     error stop 'IO_wrapLines/3'
  if ('abc, def' /= IO_wrapLines('abc, def',length=3,separator='.')) &
                                                     error stop 'IO_wrapLines/4'
  if ('abc.'//IO_EOL//'def' /= IO_wrapLines('abc. def',length=3,separator='.')) &
                                                     error stop 'IO_wrapLines/5'
  if ('abc,'//IO_EOL//'defg,'//IO_EOL//'hij' /= IO_wrapLines('abc,defg,hij',length=4)) &
                                                     error stop 'IO_wrapLines/6'
  if ('abc,'//IO_EOL//'xxdefg,'//IO_EOL//'xxhij' /= IO_wrapLines('abc,defg, hij',filler='xx',length=4)) &
                                                     error stop 'IO_wrapLines/7'

#if ((defined(__INTEL_COMPILER) && __INTEL_COMPILER_BUILD_DATE < 20240000) || !defined(__INTEL_COMPILER))
  call tokenize('','$',tokens)
  if (size(tokens) /= 0 .or. len(tokens) /=0) error stop 'tokenize empty'
  call tokenize('abcd','dcba',tokens)
  if (size(tokens) /= 0 .or. len(tokens) /=0) error stop 'tokenize only separators'

  tokens=['a']
  call test_tokenize('a','#',tokens)
  call test_tokenize('#a','#',tokens)
  call test_tokenize('a#','#',tokens)

  tokens=['aa']
  call test_tokenize('aa','#',tokens)
  call test_tokenize('$aa','$',tokens)
  call test_tokenize('aa$','$',tokens)

  tokens=['a','b']
  call test_tokenize('a$b','$',tokens)
  call test_tokenize('@a@$b@','$@',tokens)

  tokens=['aa','bb']
  call test_tokenize('aa$bb','$',tokens)
  call test_tokenize('aa$$bb','$',tokens)
  call test_tokenize('aa$bb$','$',tokens)

  tokens=['aa  ','bbb ','cccc']
  call test_tokenize('aa$bbb$cccc','$',tokens)
  call test_tokenize('$aa$bbb$cccc$','$',tokens)
  call tokenize('#aa@@bbb!!!cccc#','#@!',tokens)


  contains
  pure subroutine test_tokenize(input,delimiter,solution)
    character(len=*), intent(in) :: input, delimiter
    character(len=*), dimension(:), intent(in) :: solution

    character(len=:), dimension(:), allocatable :: tok
    integer :: i


    call tokenize(input,delimiter,tok)
    do i = 1,size(tok)
      !if (solution(i) /= tok(i)) error stop 'tokenize "'//solution(i)//'" vs. "'//tok(i)//'"'      ! requires 2018 standard
      if (solution(i) /= tok(i)) error stop 'tokenize'
    end do

  end subroutine test_tokenize
#endif

end subroutine IO_selfTest

end module IO
