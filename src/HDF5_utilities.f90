!--------------------------------------------------------------------------------------------------
!> @author Vitesh Shah, Max-Planck-Institut für Eisenforschung GmbH
!> @author Yi-Chin Yang, Max-Planck-Institut für Eisenforschung GmbH
!> @author Jennifer Nastola, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!--------------------------------------------------------------------------------------------------
module HDF5_Utilities
  use prec
  use IO
  use HDF5
#ifdef PETSc
  use PETSC
#endif

 implicit none
 private
 integer(HID_T), public, protected :: tempCoordinates, tempResults
 integer(HID_T), private :: resultsFile, currentIncID, plist_id
 integer(pInt),  private :: currentInc

!--------------------------------------------------------------------------------------------------
!> @brief reads pInt or pReal data of defined shape from file
!--------------------------------------------------------------------------------------------------
 interface HDF5_read
   module procedure HDF5_read_pReal_1
   module procedure HDF5_read_pReal_2
   module procedure HDF5_read_pReal_3
   module procedure HDF5_read_pReal_4
   module procedure HDF5_read_pReal_5
   module procedure HDF5_read_pReal_6
   module procedure HDF5_read_pReal_7

   module procedure HDF5_read_pInt_1
   module procedure HDF5_read_pInt_2
   module procedure HDF5_read_pInt_3
   module procedure HDF5_read_pInt_4
   module procedure HDF5_read_pInt_5
   module procedure HDF5_read_pInt_6
   module procedure HDF5_read_pInt_7  !ABOVE 8 DIMENSIONS IT GIVES ERROR: THE CALL TO H5DREAD_F DOESNT WORK

 end interface HDF5_read

!--------------------------------------------------------------------------------------------------
!> @brief writes pInt or pReal data of defined shape to file
!--------------------------------------------------------------------------------------------------
 interface HDF5_write
   module procedure HDF5_write_pReal1
   module procedure HDF5_write_pReal2
   module procedure HDF5_write_pReal3
   module procedure HDF5_write_pReal4
   module procedure HDF5_write_pReal5
   module procedure HDF5_write_pReal6
   module procedure HDF5_write_pReal7

   module procedure HDF5_write_pInt1
   module procedure HDF5_write_pInt2
   module procedure HDF5_write_pInt3
   module procedure HDF5_write_pInt4
   module procedure HDF5_write_pInt5
   module procedure HDF5_write_pInt6
   module procedure HDF5_write_pInt7  !ABOVE 8 DIMENSIONS IT GIVES ERROR: THE CALL TO H5DREAD_F DOESNT WORK

 end interface HDF5_write

 public :: &
   HDF5_Utilities_init, &
   HDF5_mappingPhase, &
   HDF5_mappingHomog, &
   HDF5_mappingCrystallite, &
   HDF5_backwardMappingPhase, &
   HDF5_backwardMappingHomog, &
   HDF5_backwardMappingCrystallite, &
   HDF5_mappingCells, &
   HDF5_addGroup ,&
   HDF5_closeGroup ,&
   HDF5_openGroup, &
   HDF5_openGroup2, &
   HDF5_forwardResults, &
   HDF5_writeVectorDataset, &
   HDF5_writeScalarDataset, &
   HDF5_writeTensorDataset, &
   HDF5_closeJobFile, &
   HDF5_removeLink, &
   HDF5_createFile, &
   HDF5_closeFile, &
   HDF5_addGroup2, &
   HDF5_openFile, &
   HDF5_read, &
   HDF5_write
contains

subroutine HDF5_Utilities_init
 use, intrinsic :: &
   iso_fortran_env  ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)

 implicit none
 integer              :: hdferr
 integer(SIZE_T)      :: typeSize

 write(6,'(/,a)') ' <<<+-  HDF5_Utilities init  -+>>>'
#include "compilation_info.f90"

 !currentInc = -1_pInt ToDo
 !call HDF5_createJobFile ToDo

!--------------------------------------------------------------------------------------------------
!initialize HDF5 library and check if integer and float type size match
 call h5open_f(hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_Utilities_init: h5open_f')
 call h5tget_size_f(H5T_NATIVE_INTEGER,typeSize, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_Utilities_init: h5tget_size_f (int)')
 if (int(pInt,SIZE_T)/=typeSize) call IO_error(0_pInt,ext_msg='pInt does not match H5T_NATIVE_INTEGER')
 call h5tget_size_f(H5T_NATIVE_DOUBLE,typeSize, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_Utilities_init: h5tget_size_f (double)')
 if (int(pReal,SIZE_T)/=typeSize) call IO_error(0_pInt,ext_msg='pReal does not match H5T_NATIVE_DOUBLE')

end subroutine HDF5_Utilities_init

!--------------------------------------------------------------------------------------------------
!> @brief creates and initializes HDF5 output files
!--------------------------------------------------------------------------------------------------
subroutine HDF5_createJobFile
 use hdf5
 use DAMASK_interface, only: &
   getSolverJobName

 implicit none
 integer              :: hdferr
 character(len=1024)  :: path
#ifdef PETSc
#include <petsc/finclude/petscsys.h>

 call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_Utilities_init: h5pcreate_f')
 call h5pset_fapl_mpio_f(plist_id, PETSC_COMM_WORLD, MPI_INFO_NULL, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_Utilities_init: h5pset_fapl_mpio_f')
#endif

!--------------------------------------------------------------------------------------------------
! open file
 path = trim(getSolverJobName())//'.'//'hdf5'
 !call h5fcreate_f(path,H5F_ACC_TRUNC_F,resultsFile,hdferr)
 call h5fcreate_f(path,H5F_ACC_TRUNC_F,resultsFile,hdferr,access_prp = plist_id)
 if (hdferr < 0) call IO_error(100_pInt,ext_msg=path)
 call HDF5_addStringAttribute(resultsFile,'createdBy',DAMASKVERSION)
 call h5pclose_f(plist_id, hdferr)  !neu

end subroutine HDF5_createJobFile


!--------------------------------------------------------------------------------------------------
!> @brief creates and initializes HDF5 output files
!--------------------------------------------------------------------------------------------------
 integer(HID_T) function HDF5_createFile(path)
 use hdf5
 use DAMASK_interface, only: &
   getSolverJobName

 implicit none
 integer                       :: hdferr
 integer(SIZE_T)               :: typeSize
 character(len=*), intent(in)  :: path
#ifdef PETSc
#include <petsc/finclude/petscsys.h>
#endif
  call h5open_f(hdferr) !############################################################ DANGEROUS
#ifdef PETSc
 call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_Utilities_init: h5pcreate_f')
 call h5pset_fapl_mpio_f(plist_id, PETSC_COMM_WORLD, MPI_INFO_NULL, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_Utilities_init: h5pset_fapl_mpio_f')
#endif
!--------------------------------------------------------------------------------------------------
! create a file
 !call h5fcreate_f(path,H5F_ACC_TRUNC_F,resultsFile,hdferr)
 call h5fcreate_f(path,H5F_ACC_TRUNC_F,HDF5_createFile,hdferr,access_prp = plist_id)
 if (hdferr < 0) call IO_error(100_pInt,ext_msg=path)
 !call HDF5_addStringAttribute(HDF5_createFile,'createdBy',DAMASKVERSION)
 call h5pclose_f(plist_id, hdferr)  !neu

end function HDF5_createFile

!--------------------------------------------------------------------------------------------------
!> @brief close the opened HDF5 output file
!--------------------------------------------------------------------------------------------------
subroutine HDF5_closeJobFile()
 use hdf5

 implicit none
 integer  :: hdferr
 call HDF5_removeLink('current')
 call h5fclose_f(resultsFile,hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_closeJobFile: h5fclose_f',el=hdferr)
! call h5close_f(hdferr)

end subroutine HDF5_closeJobFile

!--------------------------------------------------------------------------------------------------
!> @brief open and initializes HDF5 output file
!--------------------------------------------------------------------------------------------------
integer(HID_T) function HDF5_openFile(fileName,mode)

 implicit none
 character(len=*), intent(in)           :: fileName
 character,        intent(in), optional :: mode
 character                              :: m
 integer :: hdferr

 if (present(mode)) then
   m = mode
 else
   m = 'r'
 endif

 if (m == 'w') then
   call h5fcreate_f(fileName,H5F_ACC_TRUNC_F,HDF5_openFile,hdferr)
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_openFile: h5fcreate_f',el=hdferr)
 elseif(m == 'a') then
   call h5fopen_f(fileName,H5F_ACC_RDWR_F,HDF5_openFile,hdferr)
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_openFile: h5fopen_f (a)',el=hdferr)
 elseif(m == 'r') then
   call h5fopen_f(fileName,H5F_ACC_RDONLY_F,HDF5_openFile,hdferr)
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_openFile: h5fopen_f (r)',el=hdferr)
 else
   call IO_error(1_pInt,ext_msg='HDF5_openFile: h5fopen_f unknown access mode',el=hdferr)
 endif

end function HDF5_openFile

!--------------------------------------------------------------------------------------------------
!> @brief close the opened HDF5 output file
!--------------------------------------------------------------------------------------------------
subroutine HDF5_closeFile(fileHandle)

 implicit none
 integer :: hdferr
 integer(HID_T), intent(in) :: fileHandle
 call h5fclose_f(fileHandle,hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_closeFile: h5fclose_f',el=hdferr)

end subroutine HDF5_closeFile

!--------------------------------------------------------------------------------------------------
!> @brief adds a new group to the results file
!--------------------------------------------------------------------------------------------------
integer(HID_T) function HDF5_addGroup(groupName)
 use hdf5

 implicit none
 character(len=*), intent(in) :: groupName
 integer                      :: hdferr

 call h5gcreate_f(resultsFile, trim(groupName), HDF5_addGroup, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg = 'HDF5_addGroup: h5gcreate_f ('//trim(groupName)//')')

end function HDF5_addGroup


!--------------------------------------------------------------------------------------------------
!> @brief adds a new group to the fileHandle (additional to addGroup2)
!--------------------------------------------------------------------------------------------------
integer(HID_T) function HDF5_addGroup2(fileHandle,groupName)
 use hdf5

 implicit none
 character(len=*), intent(in) :: groupName
 integer(HID_T), intent(in)   :: fileHandle
 integer                      :: hdferr

 call h5gcreate_f(fileHandle, trim(groupName), HDF5_addGroup2, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg = 'HDF5_addGroup2: h5gcreate_f ('//trim(groupName)//')')

end function HDF5_addGroup2


!--------------------------------------------------------------------------------------------------
!> @brief open a group from the results file
!--------------------------------------------------------------------------------------------------
integer(HID_T) function HDF5_openGroup(groupName)
 use hdf5

 implicit none
 character(len=*), intent(in) :: groupName
 integer                      :: hdferr

 call h5gopen_f(resultsFile, trim(groupName), HDF5_openGroup, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg = 'HDF5_openGroup: h5gopen_f ('//trim(groupName)//')')

end function HDF5_openGroup

!--------------------------------------------------------------------------------------------------
!> @brief open an existing group of a file
!--------------------------------------------------------------------------------------------------
integer(HID_T) function HDF5_openGroup2(FileReadID,groupName)
 use hdf5

 implicit none
 character(len=*), intent(in) :: groupName
 integer                      :: hdferr
 integer(HID_T), intent(in)   :: FileReadID

 call h5gopen_f(FileReadID, trim(groupName), HDF5_openGroup2, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg = 'HDF5_openGroup2: h5gopen_f ('//trim(groupName)//')')

end function HDF5_openGroup2

!--------------------------------------------------------------------------------------------------
!> @brief set link to object in results file
!--------------------------------------------------------------------------------------------------
subroutine HDF5_setLink(path,link)
 use hdf5

 implicit none
 character(len=*), intent(in) :: path, link
 integer                      :: hdferr
 logical                      :: linkExists

 call h5lexists_f(resultsFile, link,linkExists, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg = 'HDF5_setLink: h5lexists_soft_f ('//trim(link)//')')
 if (linkExists) then
   call h5ldelete_f(resultsFile,link, hdferr)
   if (hdferr < 0) call IO_error(1_pInt,ext_msg = 'HDF5_setLink: h5ldelete_soft_f ('//trim(link)//')')
 endif
 call h5lcreate_soft_f(path, resultsFile, link, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg = 'HDF5_setLink: h5lcreate_soft_f ('//trim(path)//' '//trim(link)//')')

end subroutine HDF5_setLink

!--------------------------------------------------------------------------------------------------
!> @brief remove link to an object
!--------------------------------------------------------------------------------------------------
subroutine HDF5_removeLink(link)
 use hdf5

 implicit none
 character(len=*), intent(in) :: link
 integer                      :: hdferr

 call h5ldelete_f(resultsFile,link, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg = 'HDF5_removeLink: h5ldelete_soft_f ('//trim(link)//')')

end subroutine HDF5_removeLink

!--------------------------------------------------------------------------------------------------
!> @brief close a group
!--------------------------------------------------------------------------------------------------
subroutine HDF5_closeGroup(ID)
 use hdf5

 implicit none
 integer(HID_T), intent(in) :: ID
 integer                    :: hdferr

 call h5gclose_f(ID, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg = 'HDF5_closeGroup: h5gclose_f (el is ID)', el = int(ID,pInt))

end subroutine HDF5_closeGroup

!--------------------------------------------------------------------------------------------------
!> @brief adds a StringAttribute to the results file
!--------------------------------------------------------------------------------------------------
subroutine HDF5_addStringAttribute(entity,attrLabel,attrValue)
 use hdf5

 implicit none
 integer(HID_T),   intent(in)  :: entity
 character(len=*), intent(in)  :: attrLabel, attrValue
 integer                       :: hdferr
 integer(HID_T)                :: attr_id, space_id, type_id

 call h5screate_f(H5S_SCALAR_F,space_id,hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_addStringAttribute: h5screate_f')
 call h5tcopy_f(H5T_NATIVE_CHARACTER, type_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_addStringAttribute: h5tcopy_f')
 call h5tset_size_f(type_id, int(len(trim(attrValue)),HSIZE_T), hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_addStringAttribute: h5tset_size_f')
 call h5acreate_f(entity, trim(attrLabel),type_id,space_id,attr_id,hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_addStringAttribute: h5acreate_f')
 call h5awrite_f(attr_id, type_id, trim(attrValue), int([1],HSIZE_T), hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_addStringAttribute: h5awrite_f')
 call h5aclose_f(attr_id,hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_addStringAttribute: h5aclose_f')
 call h5tclose_f(type_id,hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_addStringAttribute: h5tclose_f')
 call h5sclose_f(space_id,hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_addStringAttribute: h5sclose_f')

end subroutine HDF5_addStringAttribute

!--------------------------------------------------------------------------------------------------
!> @brief adds the unique mapping from spatial position and constituent ID to results
!--------------------------------------------------------------------------------------------------
subroutine HDF5_mappingPhase(mapping,mapping2,Nconstituents,material_phase,phase_name,dataspace_size,mpiOffset,mpiOffset_phase)
 use hdf5

 implicit none
 integer(pInt),    intent(in)                   :: Nconstituents, dataspace_size, mpiOffset
 integer(pInt),    intent(in), dimension(:)     :: mapping, mapping2
 character(len=*), intent(in), dimension(:)     :: phase_name
 integer(pInt),    intent(in), dimension(:)     :: mpiOffset_phase
 integer(pInt),    intent(in), dimension(:,:,:) :: material_phase

 character(len=len(phase_name(1))), dimension(:), allocatable :: namesNA
 character(len=len(phase_name(1)))                            :: a
 character(len=*),                  parameter                 :: n = "NULL"

 integer(pInt)   :: hdferr, NmatPoints, i, j, k
 integer(HID_T)  :: mapping_id, dtype_id, dset_id, space_id, name_id, position_id, plist_id, memspace
 integer(HID_T)  :: dt5_id      ! Memory datatype identifier
 integer(SIZE_T) :: typesize, type_sizec, type_sizei, type_size

 integer(HSIZE_T),  dimension(2)                :: counter
 integer(HSSIZE_T), dimension(2)                :: fileOffset
 integer(pInt),     dimension(:,:), allocatable :: arrOffset

 a = n
 allocate(namesNA(0:size(phase_name)),source=[a,phase_name])
 NmatPoints = size(mapping,1)/Nconstituents
 mapping_ID = HDF5_openGroup("current/mapGeometry")

 allocate(arrOffset(Nconstituents,NmatPoints))
 do i=1_pInt, NmatPoints
   do k=1_pInt, Nconstituents
     do j=1_pInt, size(phase_name)
       if(material_phase(k,1,i) == j) &
         arrOffset(k,i) = mpiOffset_phase(j)
     enddo
   enddo
 enddo

!--------------------------------------------------------------------------------------------------
! create dataspace
 call h5screate_simple_f(2, int([Nconstituents,dataspace_size],HSIZE_T), space_id, hdferr, &
                            int([Nconstituents,dataspace_size],HSIZE_T))
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_writeMapping')

!--------------------------------------------------------------------------------------------------
! compound type
     ! First calculate total size by calculating sizes of each member
     !
     CALL h5tcopy_f(H5T_NATIVE_CHARACTER, dt5_id, hdferr)
     typesize = len(phase_name(1))
     CALL h5tset_size_f(dt5_id, typesize, hdferr)
     CALL h5tget_size_f(dt5_id, type_sizec, hdferr)
     CALL h5tget_size_f(H5T_STD_I32LE,type_sizei, hdferr)
     type_size = type_sizec + type_sizei
 call h5tcreate_f(H5T_COMPOUND_F, type_size, dtype_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_writeMapping: h5tcreate_f dtype_id')

 call h5tinsert_f(dtype_id, "Name",          0_SIZE_T, dt5_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingPhase: h5tinsert_f 0')
 call h5tinsert_f(dtype_id, "Position",  type_sizec, H5T_STD_I32LE, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingPhase: h5tinsert_f 2')

!--------------------------------------------------------------------------------------------------
! create Dataset
 call h5dcreate_f(mapping_id, 'constitutive', dtype_id, space_id, dset_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingPhase')

!--------------------------------------------------------------------------------------------------
! Create memory types (one compound datatype for each member)
 call h5tcreate_f(H5T_COMPOUND_F, int(type_sizec,SIZE_T), name_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingPhase: h5tcreate_f instance_id')
 call h5tinsert_f(name_id, "Name",        0_SIZE_T, dt5_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingPhase: h5tinsert_f instance_id')

 call h5tcreate_f(H5T_COMPOUND_F, int(pInt,SIZE_T), position_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingPhase: h5tcreate_f position_id')
 call h5tinsert_f(position_id, "Position", 0_SIZE_T, H5T_STD_I32LE, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingPhase: h5tinsert_f position_id')

!--------------------------------------------------------------------------------------------------
! Define and select hyperslabs
 counter(1)    = Nconstituents                  ! how big i am
 counter(2)    = NmatPoints
 fileOffset(1) = 0                              ! where i start to write my data
 fileOffset(2) = mpiOffset

 call h5screate_simple_f(2, counter, memspace, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingPhase: h5screate_simple_f')
 call h5dget_space_f(dset_id, space_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingPhase: h5dget_space_f')
 call h5sselect_hyperslab_f(space_id, H5S_SELECT_SET_F, fileOffset, counter, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingPhase: h5sselect_hyperslab_f')

!--------------------------------------------------------------------------------------------------
 ! Create property list for collective dataset write
#ifdef PETSc
 call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingPhase: h5pcreate_f')
 call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingPhase: h5pset_dxpl_mpio_f')
#endif

!--------------------------------------------------------------------------------------------------
! write data by fields in the datatype. Fields order is not important.
 call h5dwrite_f(dset_id, name_id, reshape(namesNA(mapping),[Nconstituents,NmatPoints]), &
                          int([Nconstituents,  dataspace_size],HSIZE_T), hdferr, &
                          file_space_id = space_id, mem_space_id = memspace, xfer_prp = plist_id)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingPhase: h5dwrite_f position_id')

 call h5dwrite_f(dset_id, position_id, reshape(mapping2-1_pInt,[Nconstituents,NmatPoints])+arrOffset, &
                          int([Nconstituents,  dataspace_size],HSIZE_T), hdferr, &
                          file_space_id = space_id, mem_space_id = memspace, xfer_prp = plist_id)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingPhase: h5dwrite_f instance_id')

!--------------------------------------------------------------------------------------------------
! close types, dataspaces
 call h5tclose_f(dtype_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingPhase: h5tclose_f dtype_id')
 call h5tclose_f(position_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingPhase: h5tclose_f position_id')
 call h5tclose_f(name_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingPhase: h5tclose_f name_id  ')
 call h5tclose_f(dt5_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingPhase: h5tclose_f dt5_id')
 call h5dclose_f(dset_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingPhase: h5dclose_f')
 call h5sclose_f(space_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingPhase: h5sclose_f space_id')
 call h5sclose_f(memspace, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingPhase: h5sclose_f memspace')
 call h5pclose_f(plist_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingPhase: h5pclose_f')
 call HDF5_closeGroup(mapping_ID)

end subroutine HDF5_mappingPhase

!--------------------------------------------------------------------------------------------------
!> @brief adds the backward mapping from spatial position and constituent ID to results
!--------------------------------------------------------------------------------------------------
subroutine HDF5_backwardMappingPhase(material_phase,phasememberat,phase_name,dataspace_size,mpiOffset,mpiOffset_phase)
 use hdf5

 implicit none
 integer(pInt),    intent(in), dimension(:,:,:) :: material_phase, phasememberat
 character(len=*), intent(in), dimension(:)     :: phase_name
 integer(pInt),    intent(in), dimension(:)     :: dataspace_size, mpiOffset_phase
 integer(pInt),    intent(in)                   :: mpiOffset

 integer(pInt)   :: hdferr, NmatPoints, Nconstituents, i, j
 integer(HID_T)  :: mapping_id, dtype_id, dset_id, space_id, position_id, plist_id, memspace
 integer(SIZE_T) :: type_size

 integer(pInt), dimension(:,:), allocatable :: arr

 integer(HSIZE_T),  dimension(1) :: counter
 integer(HSSIZE_T), dimension(1) :: fileOffset

 character(len=64) :: phaseID

 Nconstituents = size(phasememberat,1)
 NmatPoints = count(material_phase /=0_pInt)/Nconstituents

 allocate(arr(2,NmatPoints*Nconstituents))

 do i=1_pInt, NmatPoints
   do j=Nconstituents-1_pInt, 0_pInt, -1_pInt
     arr(1,Nconstituents*i-j) = i-1_pInt
   enddo
 enddo
 arr(2,:) = pack(material_phase,material_phase/=0_pInt)

 do i=1_pInt, size(phase_name)
   write(phaseID, '(i0)') i
   mapping_ID = HDF5_openGroup('/current/constitutive/'//trim(phaseID)//'_'//phase_name(i))
   NmatPoints = count(material_phase == i)

!--------------------------------------------------------------------------------------------------
  ! create dataspace
   call h5screate_simple_f(1, int([dataspace_size(i)],HSIZE_T), space_id, hdferr, &
                              int([dataspace_size(i)],HSIZE_T))
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_writeBackwardMapping')

!--------------------------------------------------------------------------------------------------
  ! compound type
   call h5tget_size_f(H5T_STD_I32LE, type_size, hdferr)
   call h5tcreate_f(H5T_COMPOUND_F, type_size, dtype_id, hdferr)
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_writeBackwardMapping: h5tcreate_f dtype_id')

   call h5tinsert_f(dtype_id, "Position",          0_SIZE_T, H5T_STD_I32LE, hdferr)
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_backwardMappingPhase: h5tinsert_f 0')

!--------------------------------------------------------------------------------------------------
  ! create Dataset
   call h5dcreate_f(mapping_id, 'mapGeometry', dtype_id, space_id, dset_id, hdferr)
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_backwardMappingPhase')

!--------------------------------------------------------------------------------------------------
  ! Create memory types (one compound datatype for each member)
   call h5tcreate_f(H5T_COMPOUND_F, int(pInt,SIZE_T), position_id, hdferr)
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_backwardMappingPhase: h5tcreate_f position_id')
   call h5tinsert_f(position_id, "Position", 0_SIZE_T, H5T_STD_I32LE, hdferr)
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_backwardMappingPhase: h5tinsert_f position_id')

!--------------------------------------------------------------------------------------------------
  ! Define and select hyperslabs
   counter = NmatPoints                       ! how big i am
   fileOffset = mpiOffset_phase(i)            ! where i start to write my data

   call h5screate_simple_f(1, counter, memspace, hdferr)
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_backwardMappingPhase: h5screate_simple_f')
   call h5dget_space_f(dset_id, space_id, hdferr)
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_backwardMappingPhase: h5dget_space_f')
   call h5sselect_hyperslab_f(space_id, H5S_SELECT_SET_F, fileOffset, counter, hdferr)
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_backwardMappingPhase: h5sselect_hyperslab_f')

!--------------------------------------------------------------------------------------------------
 ! Create property list for collective dataset write
#ifdef PETSc
   call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdferr)
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_backwardMappingPhase: h5pcreate_f')
   call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdferr)
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_backwardMappingPhase: h5pset_dxpl_mpio_f')
#endif

!--------------------------------------------------------------------------------------------------
  ! write data by fields in the datatype. Fields order is not important.
   call h5dwrite_f(dset_id, position_id, pack(arr(1,:),arr(2,:)==i)+mpiOffset, int([dataspace_size(i)],HSIZE_T),&
                              hdferr, file_space_id = space_id, mem_space_id = memspace, xfer_prp = plist_id)
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_backwardMappingPhase: h5dwrite_f instance_id')

!--------------------------------------------------------------------------------------------------
  !close types, dataspaces
   call h5tclose_f(dtype_id, hdferr)
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_backwardMappingPhase: h5tclose_f dtype_id')
   call h5tclose_f(position_id, hdferr)
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_backwardMappingPhase: h5tclose_f position_id')
   call h5dclose_f(dset_id, hdferr)
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_backwardMappingPhase: h5dclose_f')
   call h5sclose_f(space_id, hdferr)
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_backwardMappingPhase: h5sclose_f space_id')
   call h5sclose_f(memspace, hdferr)
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_backwardMappingPhase: h5sclose_f memspace')
   call h5pclose_f(plist_id, hdferr)
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_backwardMappingPhase: h5pclose_f')
   call HDF5_closeGroup(mapping_ID)

 enddo

end subroutine HDF5_backwardMappingPhase

!--------------------------------------------------------------------------------------------------
!> @brief adds the unique mapping from spatial position and constituent ID to results
!--------------------------------------------------------------------------------------------------
subroutine HDF5_mappingHomog(material_homog,homogmemberat,homogenization_name,dataspace_size,mpiOffset,mpiOffset_homog)
 use hdf5

 implicit none
 integer(pInt),    intent(in), dimension(:,:) :: material_homog, homogmemberat
 character(len=*), intent(in), dimension(:)   :: homogenization_name
 integer(pInt),    intent(in), dimension(:)   :: mpiOffset_homog
 integer(pInt),    intent(in)                 :: dataspace_size, mpiOffset

 integer(pInt)   :: hdferr, NmatPoints, i, j
 integer(HID_T)  :: mapping_id, dtype_id, dset_id, space_id, name_id, position_id, plist_id, memspace

 integer(HID_T)  :: dt5_id      ! Memory datatype identifier
 integer(SIZE_T) :: typesize, type_sizec, type_sizei, type_size

 integer(HSIZE_T),  dimension(1)              :: counter
 integer(HSSIZE_T), dimension(1)              :: fileOffset
 integer(pInt),     dimension(:), allocatable :: arrOffset

 NmatPoints = count(material_homog /=0_pInt)
 mapping_ID = HDF5_openGroup("current/mapGeometry")

 allocate(arrOffset(NmatPoints))
 do i=1_pInt, NmatPoints
   do j=1_pInt, size(homogenization_name)
     if(material_homog(1,i) == j) &
       arrOffset(i) = mpiOffset_homog(j)
   enddo
 enddo

!--------------------------------------------------------------------------------------------------
! create dataspace
 call h5screate_simple_f(1, int([dataspace_size],HSIZE_T), space_id, hdferr, &
                            int([dataspace_size],HSIZE_T))
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_writeMapping')

!--------------------------------------------------------------------------------------------------
! compound type
     ! First calculate total size by calculating sizes of each member
     !
     CALL h5tcopy_f(H5T_NATIVE_CHARACTER, dt5_id, hdferr)
     typesize = len(homogenization_name(1))
     CALL h5tset_size_f(dt5_id, typesize, hdferr)
     CALL h5tget_size_f(dt5_id, type_sizec, hdferr)
     CALL h5tget_size_f(H5T_STD_I32LE,type_sizei, hdferr)
     type_size = type_sizec + type_sizei
 call h5tcreate_f(H5T_COMPOUND_F, type_size, dtype_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_writeMapping: h5tcreate_f dtype_id')

 call h5tinsert_f(dtype_id, "Name",          0_SIZE_T, dt5_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingHomog: h5tinsert_f 0')
 call h5tinsert_f(dtype_id, "Position",  type_sizec, H5T_STD_I32LE, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingHomog: h5tinsert_f 2')

!--------------------------------------------------------------------------------------------------
! create Dataset
 call h5dcreate_f(mapping_id, 'homogenization', dtype_id, space_id, dset_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingHomog')

!--------------------------------------------------------------------------------------------------
! Create memory types (one compound datatype for each member)
 call h5tcreate_f(H5T_COMPOUND_F, int(type_sizec,SIZE_T), name_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingHomog: h5tcreate_f instance_id')
 call h5tinsert_f(name_id, "Name",        0_SIZE_T, dt5_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingHomog: h5tinsert_f instance_id')

 call h5tcreate_f(H5T_COMPOUND_F, int(pInt,SIZE_T), position_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingHomog: h5tcreate_f position_id')
 call h5tinsert_f(position_id, "Position", 0_SIZE_T, H5T_STD_I32LE, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingHomog: h5tinsert_f position_id')

!--------------------------------------------------------------------------------------------------
! Define and select hyperslabs
 counter = NmatPoints                    ! how big i am
 fileOffset = mpiOffset                  ! where i start to write my data

 call h5screate_simple_f(1, counter, memspace, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingHomog: h5screate_simple_f')
 call h5dget_space_f(dset_id, space_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingHomog: h5dget_space_f')
 call h5sselect_hyperslab_f(space_id, H5S_SELECT_SET_F, fileOffset, counter, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingHomog: h5sselect_hyperslab_f')

!--------------------------------------------------------------------------------------------------
! Create property list for collective dataset write
#ifdef PETSc
 call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingHomog: h5pcreate_f')
 call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingHomog: h5pset_dxpl_mpio_f')
#endif

!--------------------------------------------------------------------------------------------------
! write data by fields in the datatype. Fields order is not important.
 call h5dwrite_f(dset_id, name_id, homogenization_name(pack(material_homog,material_homog/=0_pInt)), &
                          int([dataspace_size],HSIZE_T), hdferr, file_space_id = space_id, &
                          mem_space_id = memspace, xfer_prp = plist_id)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingHomog: h5dwrite_f position_id')

 call h5dwrite_f(dset_id, position_id, pack(homogmemberat-1_pInt,homogmemberat/=0_pInt) + arrOffset, &
                          int([dataspace_size],HSIZE_T), hdferr, file_space_id = space_id, &
                          mem_space_id = memspace, xfer_prp = plist_id)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingHomog: h5dwrite_f instance_id')

!--------------------------------------------------------------------------------------------------
!close types, dataspaces
call h5tclose_f(dtype_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingHomog: h5tclose_f dtype_id')
call h5tclose_f(position_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingHomog: h5tclose_f position_id')
call h5tclose_f(name_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingHomog: h5tclose_f name_id  ')
call h5tclose_f(dt5_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingHomog: h5tclose_f dt5_id')
call h5dclose_f(dset_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingHomog: h5dclose_f')
call h5sclose_f(space_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingHomog: h5sclose_f space_id')
call h5sclose_f(memspace, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingHomog: h5sclose_f memspace')
call h5pclose_f(plist_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingHomog: h5pclose_f')
call HDF5_closeGroup(mapping_ID)

end subroutine HDF5_mappingHomog

!--------------------------------------------------------------------------------------------------
!> @brief adds the backward mapping from spatial position and constituent ID to results
!--------------------------------------------------------------------------------------------------
subroutine HDF5_backwardMappingHomog(material_homog,homogmemberat,homogenization_name,dataspace_size,mpiOffset,mpiOffset_homog)
 use hdf5

 implicit none
 integer(pInt),    intent(in), dimension(:,:) :: material_homog, homogmemberat
 character(len=*), intent(in), dimension(:)   :: homogenization_name
 integer(pInt),    intent(in), dimension(:)   :: dataspace_size, mpiOffset_homog
 integer(pInt),    intent(in)                 :: mpiOffset

 integer(pInt)   :: hdferr, NmatPoints, i
 integer(HID_T)  :: mapping_id, dtype_id, dset_id, space_id, position_id, plist_id, memspace
 integer(SIZE_T) :: type_size

 integer(pInt),     dimension(:,:), allocatable :: arr

 integer(HSIZE_T),  dimension(1) :: counter
 integer(HSSIZE_T), dimension(1) :: fileOffset

 character(len=64) :: homogID

 NmatPoints = count(material_homog /=0_pInt)
 allocate(arr(2,NmatPoints))

 arr(1,:) = (/(i, i=0_pint,NmatPoints-1_pInt)/)
 arr(2,:) = pack(material_homog,material_homog/=0_pInt)

 do i=1_pInt, size(homogenization_name)
   write(homogID, '(i0)') i
   mapping_ID = HDF5_openGroup('/current/homogenization/'//trim(homogID)//'_'//homogenization_name(i))

!--------------------------------------------------------------------------------------------------
  ! create dataspace
   call h5screate_simple_f(1, int([dataspace_size(i)],HSIZE_T), space_id, hdferr, &
                              int([dataspace_size(i)],HSIZE_T))
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_writeBackwardMapping')

!--------------------------------------------------------------------------------------------------
  ! compound type
   call h5tget_size_f(H5T_STD_I32LE, type_size, hdferr)
   call h5tcreate_f(H5T_COMPOUND_F, type_size, dtype_id, hdferr)
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_writeBackwardMapping: h5tcreate_f dtype_id')

   call h5tinsert_f(dtype_id, "Position",          0_SIZE_T, H5T_STD_I32LE, hdferr)
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_backwardMappingHomog: h5tinsert_f 0')

!--------------------------------------------------------------------------------------------------
  ! create Dataset
   call h5dcreate_f(mapping_id, 'mapGeometry', dtype_id, space_id, dset_id, hdferr)
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_backwardMappingHomog')

!--------------------------------------------------------------------------------------------------
  ! Create memory types (one compound datatype for each member)
   call h5tcreate_f(H5T_COMPOUND_F, int(pInt,SIZE_T), position_id, hdferr)
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_backwardMappingHomog: h5tcreate_f position_id')
   call h5tinsert_f(position_id, "Position", 0_SIZE_T, H5T_STD_I32LE, hdferr)
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_backwardMappingHomog: h5tinsert_f position_id')

!--------------------------------------------------------------------------------------------------
  ! Define and select hyperslabs
   counter = NmatPoints                           ! how big i am
   fileOffset = mpiOffset_homog(i)                ! where i start to write my data

   call h5screate_simple_f(1, counter, memspace, hdferr)
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_backwardMappingHomog: h5screate_simple_f')
   call h5dget_space_f(dset_id, space_id, hdferr)
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_backwardMappingHomog: h5dget_space_f')
   call h5sselect_hyperslab_f(space_id, H5S_SELECT_SET_F, fileOffset, counter, hdferr)
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_backwardMappingHomog: h5sselect_hyperslab_f')

!--------------------------------------------------------------------------------------------------
 ! Create property list for collective dataset write
#ifdef PETSc
   call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdferr)
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_backwardMappingHomog: h5pcreate_f')
   call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdferr)
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_backwardMappingHomog: h5pset_dxpl_mpio_f')
#endif

!--------------------------------------------------------------------------------------------------
  ! write data by fields in the datatype. Fields order is not important.
   call h5dwrite_f(dset_id, position_id, pack(arr(1,:),arr(2,:)==i)+mpiOffset,int([dataspace_size(i)],HSIZE_T),&
                   hdferr, file_space_id = space_id, mem_space_id = memspace, xfer_prp = plist_id)
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_backwardMappingHomog: h5dwrite_f instance_id')

!--------------------------------------------------------------------------------------------------
  !close types, dataspaces
   call h5tclose_f(dtype_id, hdferr)
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_backwardMappingHomog: h5tclose_f dtype_id')
   call h5tclose_f(position_id, hdferr)
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_backwardMappingHomog: h5tclose_f position_id')
   call h5dclose_f(dset_id, hdferr)
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_backwardMappingHomog: h5dclose_f')
   call h5sclose_f(space_id, hdferr)
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_backwardMappingHomog: h5sclose_f space_id')
   call h5sclose_f(memspace, hdferr)
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_backwardMappingHomog: h5sclose_f memspace')
   call h5pclose_f(plist_id, hdferr)
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_backwardMappingHomog: h5pclose_f')
   call HDF5_closeGroup(mapping_ID)

 enddo

end subroutine HDF5_backwardMappingHomog

!--------------------------------------------------------------------------------------------------
!> @brief adds the unique mapping from spatial position and constituent ID to results
!--------------------------------------------------------------------------------------------------
subroutine HDF5_mappingCrystallite(crystalliteAt,crystmemberAt,crystallite_name,dataspace_size,mpiOffset,mpiOffset_cryst)
 use hdf5

 implicit none
 integer(pInt),    intent(in), dimension(:,:)   :: crystalliteAt
 integer(pInt),    intent(in), dimension(:,:,:) :: crystmemberAt
 character(len=*), intent(in), dimension(:)     :: crystallite_name
 integer(pInt),    intent(in), dimension(:)     :: mpiOffset_cryst
 integer(pInt),    intent(in)                   :: dataspace_size, mpiOffset

 integer         :: hdferr
 integer(pInt)   :: NmatPoints, Nconstituents, i, j
 integer(HID_T)  :: mapping_id, dtype_id, dset_id, space_id, name_id, plist_id, memspace

 integer(HID_T), dimension(:), allocatable :: position_id

 integer(HID_T)  :: dt5_id      ! Memory datatype identifier
 integer(SIZE_T) :: typesize, type_sizec, type_sizei, type_size

 integer(HSIZE_T),  dimension(1)              :: counter
 integer(HSSIZE_T), dimension(1)              :: fileOffset
 integer(pInt),     dimension(:), allocatable :: arrOffset

 character(len=64) :: m

 Nconstituents = size(crystmemberAt,1)
 NmatPoints = count(crystalliteAt /=0_pInt)
 mapping_ID = HDF5_openGroup("current/mapGeometry")

 allocate(position_id(Nconstituents))

 allocate(arrOffset(NmatPoints))
 do i=1_pInt, NmatPoints
   do j=1_pInt, size(crystallite_name)
     if(crystalliteAt(1,i) == j) &
       arrOffset(i) = Nconstituents*mpiOffset_cryst(j)
   enddo
 enddo

!--------------------------------------------------------------------------------------------------
! create dataspace
 call h5screate_simple_f(1, int([dataspace_size],HSIZE_T), space_id, hdferr, &
                            int([dataspace_size],HSIZE_T))
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_writeMapping')

!--------------------------------------------------------------------------------------------------
! compound type
     ! First calculate total size by calculating sizes of each member
     !
     CALL h5tcopy_f(H5T_NATIVE_CHARACTER, dt5_id, hdferr)
     typesize = len(crystallite_name(1))
     CALL h5tset_size_f(dt5_id, typesize, hdferr)
     CALL h5tget_size_f(dt5_id, type_sizec, hdferr)
     CALL h5tget_size_f(H5T_STD_I32LE, type_sizei, hdferr)
     type_size = type_sizec + type_sizei*Nconstituents
 call h5tcreate_f(H5T_COMPOUND_F, type_size, dtype_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_writeMapping: h5tcreate_f dtype_id')

 call h5tinsert_f(dtype_id, "Name",          0_SIZE_T, dt5_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingCrystallite: h5tinsert_f 0')
 do i=1_pInt, Nconstituents
   write(m, '(i0)') i
   call h5tinsert_f(dtype_id, "Position "//trim(m),  type_sizec+(i-1)*type_sizei, H5T_STD_I32LE, hdferr)
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingCrystallite: h5tinsert_f 2 '//trim(m))
 enddo

!--------------------------------------------------------------------------------------------------
! create Dataset
 call h5dcreate_f(mapping_id, 'crystallite', dtype_id, space_id, dset_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingCrystallite')

!--------------------------------------------------------------------------------------------------
! Create memory types (one compound datatype for each member)
 call h5tcreate_f(H5T_COMPOUND_F, int(type_sizec,SIZE_T), name_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingCrystallite: h5tcreate_f instance_id')
 call h5tinsert_f(name_id, "Name",        0_SIZE_T, dt5_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingCrystallite: h5tinsert_f instance_id')

 do i=1_pInt, Nconstituents
   write(m, '(i0)') i
   call h5tcreate_f(H5T_COMPOUND_F, int(pInt,SIZE_T), position_id(i), hdferr)
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingCrystallite: h5tcreate_f position_id')
   call h5tinsert_f(position_id(i), "Position "//trim(m), 0_SIZE_T, H5T_STD_I32LE, hdferr)
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingCrystallite: h5tinsert_f position_id')
 enddo

!--------------------------------------------------------------------------------------------------
! Define and select hyperslabs
 counter = NmatPoints                    ! how big i am
 fileOffset = mpiOffset                  ! where i start to write my data

 call h5screate_simple_f(1, counter, memspace, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingCrystallite: h5screate_simple_f')
 call h5dget_space_f(dset_id, space_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingCrystallite: h5dget_space_f')
 call h5sselect_hyperslab_f(space_id, H5S_SELECT_SET_F, fileOffset, counter, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingCrystallite: h5sselect_hyperslab_f')

!--------------------------------------------------------------------------------------------------
 ! Create property list for collective dataset write
#ifdef PETSc
 call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingCrystallite: h5pcreate_f')
 call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingCrystallite: h5pset_dxpl_mpio_f')
#endif

!--------------------------------------------------------------------------------------------------
! write data by fields in the datatype. Fields order is not important.
 call h5dwrite_f(dset_id, name_id, crystallite_name(pack(crystalliteAt,crystalliteAt/=0_pInt)), &
                          int([dataspace_size],HSIZE_T), hdferr, file_space_id = space_id, &
                          mem_space_id = memspace, xfer_prp = plist_id)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingCrystallite: h5dwrite_f position_id')

 do i=1_pInt, Nconstituents
   call h5dwrite_f(dset_id, position_id(i), pack(crystmemberAt(i,:,:)-1_pInt,crystmemberAt(i,:,:)/=0_pInt)+arrOffset,&
                            int([dataspace_size],HSIZE_T), hdferr, file_space_id = space_id, &
                            mem_space_id = memspace, xfer_prp = plist_id)
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingCrystallite: h5dwrite_f instance_id')
 enddo

!--------------------------------------------------------------------------------------------------
!close types, dataspaces
 call h5tclose_f(dtype_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingCrystallite: h5tclose_f dtype_id')
 do i=1_pInt, Nconstituents
   call h5tclose_f(position_id(i), hdferr)
 enddo
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingCrystallite: h5tclose_f position_id')
 call h5tclose_f(name_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingCrystallite: h5tclose_f name_id')
 call h5tclose_f(dt5_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingCrystallite: h5tclose_f dt5_id')
 call h5dclose_f(dset_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingCrystallite: h5dclose_f')
 call h5sclose_f(space_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingCrystallite: h5sclose_f space_id')
 call h5sclose_f(memspace, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingCrystallite: h5sclose_f memspace')
 call h5pclose_f(plist_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingCrystallite: h5pclose_f')
 call HDF5_closeGroup(mapping_ID)

end subroutine HDF5_mappingCrystallite


!--------------------------------------------------------------------------------------------------
!> @brief adds the backward mapping from spatial position and constituent ID to results
!--------------------------------------------------------------------------------------------------
subroutine HDF5_backwardMappingCrystallite(crystalliteAt,crystmemberAt,crystallite_name,dataspace_size,mpiOffset,mpiOffset_cryst)
 use hdf5

 implicit none
 integer(pInt),    intent(in), dimension(:,:)   :: crystalliteAt
 integer(pInt),    intent(in), dimension(:,:,:) :: crystmemberAt
 character(len=*), intent(in), dimension(:)     :: crystallite_name
 integer(pInt),    intent(in), dimension(:)     :: dataspace_size, mpiOffset_cryst
 integer(pInt),    intent(in)                   :: mpiOffset

 integer         :: hdferr
 integer(pInt)   :: NmatPoints, Nconstituents, i, j
 integer(HID_T)  :: mapping_id, dtype_id, dset_id, space_id, position_id, plist_id, memspace
 integer(SIZE_T) :: type_size

 integer(pInt), dimension(:,:), allocatable     :: h_arr, arr

 integer(HSIZE_T),  dimension(1) :: counter
 integer(HSSIZE_T), dimension(1) :: fileOffset

 character(len=64) :: crystallID

 Nconstituents = size(crystmemberAt,1)
 NmatPoints = count(crystalliteAt /=0_pInt)

 allocate(h_arr(2,NmatPoints))
 allocate(arr(2,Nconstituents*NmatPoints))

 h_arr(1,:) = (/(i, i=0_pInt,NmatPoints-1_pInt)/)
 h_arr(2,:) = pack(crystalliteAt,crystalliteAt/=0_pInt)

 do i=1_pInt, NmatPoints
   do j=Nconstituents-1_pInt, 0_pInt, -1_pInt
     arr(1,Nconstituents*i-j) = h_arr(1,i)
     arr(2,Nconstituents*i-j) = h_arr(2,i)
   enddo
 enddo

 do i=1_pInt, size(crystallite_name)
   if (crystallite_name(i) == 'none') cycle
   write(crystallID, '(i0)') i
   mapping_ID = HDF5_openGroup('/current/crystallite/'//trim(crystallID)//'_'//crystallite_name(i))
   NmatPoints = count(crystalliteAt == i)

!--------------------------------------------------------------------------------------------------
  ! create dataspace
   call h5screate_simple_f(1, int([Nconstituents*dataspace_size(i)],HSIZE_T), space_id, hdferr, &
                              int([Nconstituents*dataspace_size(i)],HSIZE_T))
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_writeBackwardMapping')

!--------------------------------------------------------------------------------------------------
  ! compound type
   call h5tget_size_f(H5T_STD_I32LE, type_size, hdferr)
   call h5tcreate_f(H5T_COMPOUND_F, type_size, dtype_id, hdferr)
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_writeBackwardMapping: h5tcreate_f dtype_id')

   call h5tinsert_f(dtype_id, "Position",          0_SIZE_T, H5T_STD_I32LE, hdferr)
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_backwardMappingCrystallite: h5tinsert_f 0')

!--------------------------------------------------------------------------------------------------
  ! create Dataset
   call h5dcreate_f(mapping_id, 'mapGeometry', dtype_id, space_id, dset_id, hdferr)
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_backwardMappingCrystallite')

!--------------------------------------------------------------------------------------------------
  ! Create memory types
   call h5tcreate_f(H5T_COMPOUND_F, int(pInt,SIZE_T), position_id, hdferr)
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_backwardMappingCrystallite: h5tcreate_f position_id')
   call h5tinsert_f(position_id, "Position", 0_SIZE_T, H5T_STD_I32LE, hdferr)
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_backwardMappingCrystallite: h5tinsert_f position_id')

!--------------------------------------------------------------------------------------------------
  ! Define and select hyperslabs
   counter = Nconstituents*NmatPoints                       ! how big i am
   fileOffset = Nconstituents*mpiOffset_cryst(i)            ! where i start to write my data

   call h5screate_simple_f(1, counter, memspace, hdferr)
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_backwardMappingCrystallite: h5screate_simple_f')
   call h5dget_space_f(dset_id, space_id, hdferr)
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_backwardMappingCrystallite: h5dget_space_f')
   call h5sselect_hyperslab_f(space_id, H5S_SELECT_SET_F, fileOffset, counter, hdferr)
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_backwardMappingCrystallite: h5sselect_hyperslab_f')

!--------------------------------------------------------------------------------------------------
 ! Create property list for collective dataset write
#ifdef PETSc
   call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdferr)
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_backwardMappingCrystallite: h5pcreate_f')
   call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdferr)
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_backwardMappingCrystallite: h5pset_dxpl_mpio_f')
#endif

!--------------------------------------------------------------------------------------------------
  ! write data by fields in the datatype. Fields order is not important.
   call h5dwrite_f(dset_id, position_id, pack(arr(1,:),arr(2,:)==i) + mpiOffset,&
                              int([Nconstituents*dataspace_size(i)],HSIZE_T), hdferr, file_space_id = space_id, &
                              mem_space_id = memspace, xfer_prp = plist_id)
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_backwardMappingCrystallite: h5dwrite_f instance_id')

!--------------------------------------------------------------------------------------------------
  !close types, dataspaces
   call h5tclose_f(dtype_id, hdferr)
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_backwardMappingCrystallite: h5tclose_f dtype_id')
   call h5tclose_f(position_id, hdferr)
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_backwardMappingCrystallite: h5tclose_f position_id')
   call h5dclose_f(dset_id, hdferr)
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_backwardMappingCrystallite: h5dclose_f')
   call h5sclose_f(space_id, hdferr)
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_backwardMappingCrystallite: h5sclose_f space_id')
   call h5sclose_f(memspace, hdferr)
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_backwardMappingCrystallite: h5sclose_f memspace')
   call h5pclose_f(plist_id, hdferr)
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_backwardMappingCrystallite: h5pclose_f')
   call HDF5_closeGroup(mapping_ID)

 enddo

end subroutine HDF5_backwardMappingCrystallite

!--------------------------------------------------------------------------------------------------
!> @brief adds the unique cell to node mapping
!--------------------------------------------------------------------------------------------------
subroutine HDF5_mappingCells(mapping)
 use hdf5

 implicit none
 integer(pInt), intent(in), dimension(:) :: mapping

 integer        :: hdferr, Nnodes
 integer(HID_T) :: mapping_id, dset_id, space_id

 Nnodes=size(mapping)
 mapping_ID = HDF5_openGroup("mapping")

!--------------------------------------------------------------------------------------------------
! create dataspace
 call h5screate_simple_f(1, int([Nnodes],HSIZE_T), space_id, hdferr, &
                            int([Nnodes],HSIZE_T))
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingCells: h5screate_simple_f')

!--------------------------------------------------------------------------------------------------
! create Dataset
 call h5dcreate_f(mapping_id, "Cell",H5T_NATIVE_INTEGER, space_id, dset_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingCells')

!--------------------------------------------------------------------------------------------------
! write data by fields in the datatype. Fields order is not important.
 call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, mapping, int([Nnodes],HSIZE_T), hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingCells: h5dwrite_f instance_id')

!--------------------------------------------------------------------------------------------------
!close types, dataspaces
 call h5dclose_f(dset_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingConstitutive: h5dclose_f')
 call h5sclose_f(space_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='IO_mappingConstitutive: h5sclose_f')
 call HDF5_closeGroup(mapping_ID)

end subroutine HDF5_mappingCells

!--------------------------------------------------------------------------------------------------
!> @brief creates a new 3D Tensor dataset in the given group location   !!!TODO: really necessary?
!--------------------------------------------------------------------------------------------------
subroutine HDF5_addTensor3DDataset(group,Nnodes,tensorSize,label,SIunit)
 use hdf5

 implicit none
 integer(HID_T),   intent(in)   :: group
 integer(pInt),    intent(in)   :: Nnodes, tensorSize
 character(len=*), intent(in)   :: SIunit, label

 integer                        :: hdferr
 integer(HID_T)                 :: space_id, dset_id
 integer(HSIZE_T), dimension(3) :: dataShape

 dataShape = int([tensorSize,tensorSize,Nnodes], HSIZE_T)

!--------------------------------------------------------------------------------------------------
! create dataspace
 call h5screate_simple_f(3, dataShape, space_id, hdferr, dataShape)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_addTensor3DDataset: h5screate_simple_f')

!--------------------------------------------------------------------------------------------------
! create Dataset
 call h5dcreate_f(group, trim(label),H5T_NATIVE_DOUBLE, space_id, dset_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_addTensor3DDataset: h5dcreate_f')
 call HDF5_addStringAttribute(dset_id,'unit',trim(SIunit))

!--------------------------------------------------------------------------------------------------
!close types, dataspaces
 call h5dclose_f(dset_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_addTensor3DDataset: h5dclose_f')
 call h5sclose_f(space_id, hdferr)

end subroutine HDF5_addTensor3DDataset

!--------------------------------------------------------------------------------------------------
!> @brief creates a new vector dataset in the given group location   !!!TODO: really necessary?
!--------------------------------------------------------------------------------------------------
subroutine HDF5_writeVectorDataset(group,dataset,label,SIunit,dataspace_size,mpiOffset)
 use hdf5

 implicit none
 integer(HID_T),   intent(in)                 :: group
 character(len=*), intent(in)                 :: SIunit,label
 integer(pInt),    intent(in)                 :: dataspace_size, mpiOffset
 real(pReal),      intent(in), dimension(:,:) :: dataset

 integer        :: hdferr, vectorSize
 integer(HID_T) :: dset_id, space_id, memspace, plist_id

 integer(HSIZE_T),  dimension(2) :: counter
 integer(HSSIZE_T), dimension(2) :: fileOffset

 if(any(shape(dataset) == 0)) return

 vectorSize = size(dataset,1)

 call HDF5_addVectorDataset(group,dataspace_size,vectorSize,label,SIunit)                   ! here nNodes need to be global
 call h5dopen_f(group, label, dset_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_writeVectorDataset: h5dopen_f')

 ! Define and select hyperslabs
 counter(1) = vectorSize                    ! how big i am
 counter(2) = size(dataset,2)
 fileOffset(1)  = 0                         ! where i start to write my data
 fileOffset(2)  = mpiOffset

 call h5screate_simple_f(2, counter, memspace, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_writeVectorDataset: h5screate_simple_f')
 call h5dget_space_f(dset_id, space_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_writeVectorDataset: h5dget_space_f')
 call h5sselect_hyperslab_f(space_id, H5S_SELECT_SET_F, fileOffset, counter, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_writeVectorDataset: h5sselect_hyperslab_f')

 ! Create property list for collective dataset write
#ifdef PETSc
 call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_writeVectorDataset: h5pcreate_f')
 call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_writeVectorDataset: h5pset_dxpl_mpio_f')
#endif

 ! Write the dataset collectively
 call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, dataset, int([vectorSize, dataspace_size],HSIZE_T), hdferr, &
                 file_space_id = space_id, mem_space_id = memspace, xfer_prp = plist_id)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_writeVectorDataset: h5dwrite_f')

 call h5sclose_f(space_id, hdferr)
 call h5sclose_f(memspace, hdferr)
 call h5dclose_f(dset_id, hdferr)
 call h5pclose_f(plist_id, hdferr)

end subroutine HDF5_writeVectorDataset

!--------------------------------------------------------------------------------------------------
!> @brief creates a new tensor dataset in the given group location
!  by default, a 3x3 tensor is assumed                              !!!TODO: really necessary?
!--------------------------------------------------------------------------------------------------
subroutine HDF5_writeTensorDataset(group,dataset,label,SIunit,dataspace_size,mpiOffset)
 use hdf5

 implicit none
 integer(HID_T),   intent(in)                   :: group
 character(len=*), intent(in)                   :: SIunit,label
 integer(pInt),    intent(in)                   :: dataspace_size, mpiOffset
 real(pReal),      intent(in), dimension(:,:,:) :: dataset

 integer        :: hdferr, tensorSize
 integer(HID_T) :: dset_id, space_id, memspace, plist_id

 integer(HSIZE_T),  dimension(3)                :: counter
 integer(HSSIZE_T), dimension(3)                :: fileOffset

 if(any(shape(dataset) == 0)) return

 tensorSize = size(dataset,1)

 call HDF5_addTensor3DDataset(group,dataspace_size,tensorSize,label,SIunit)
 call h5dopen_f(group, label, dset_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_writeTensorDataset: h5dopen_f')

 ! Define and select hyperslabs
 counter(1) = tensorSize                    ! how big i am
 counter(2) = tensorSize
 counter(3) = size(dataset,3)
 fileOffset(1)  = 0                         ! where i start to write my data
 fileOffset(2)  = 0
 fileOffset(3)  = mpiOffset

 call h5screate_simple_f(3, counter, memspace, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_writeTensorDataset: h5screate_simple_f')
 call h5dget_space_f(dset_id, space_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_writeTensorDataset: h5dget_space_f')
 call h5sselect_hyperslab_f(space_id, H5S_SELECT_SET_F, fileOffset, counter, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_writeTensorDataset: h5sselect_hyperslab_f')

 ! Create property list for collective dataset write
#ifdef PETSc
 call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_writeTensorDataset: h5pcreate_f')
 call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_writeTensorDataset: h5pset_dxpl_mpio_f')
#endif

 ! Write the dataset collectively
 call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, dataset, int([tensorSize, dataspace_size],HSIZE_T), hdferr, &
                 file_space_id = space_id, mem_space_id = memspace, xfer_prp = plist_id)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_writeTensorDataset: h5dwrite_f')

 call h5sclose_f(space_id, hdferr)
 call h5sclose_f(memspace, hdferr)
 call h5dclose_f(dset_id, hdferr)
 call h5pclose_f(plist_id, hdferr)

 end subroutine HDF5_writeTensorDataset

!--------------------------------------------------------------------------------------------------
!> @brief adds a new vector dataset to the given group location
!--------------------------------------------------------------------------------------------------
subroutine HDF5_addVectorDataset(group,nnodes,vectorSize,label,SIunit)
 use hdf5

 implicit none
 integer(HID_T),   intent(in) :: group
 integer(pInt),    intent(in) :: nnodes,vectorSize
 character(len=*), intent(in) :: SIunit,label

 integer        :: hdferr
 integer(HID_T) :: space_id, dset_id

!--------------------------------------------------------------------------------------------------
! create dataspace
 call h5screate_simple_f(2, int([vectorSize,Nnodes],HSIZE_T), space_id, hdferr, &
                            int([vectorSize,Nnodes],HSIZE_T))
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_addVectorDataset: h5screate_simple_f')

!--------------------------------------------------------------------------------------------------
! create Dataset
 call h5dcreate_f(group, trim(label), H5T_NATIVE_DOUBLE, space_id, dset_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_addVectorDataset: h5dcreate_f')
 call HDF5_addStringAttribute(dset_id,'unit',trim(SIunit))

!--------------------------------------------------------------------------------------------------
!close types, dataspaces
 call h5dclose_f(dset_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_addVectorDataset: h5dclose_f')
 call h5sclose_f(space_id, hdferr)

end subroutine HDF5_addVectorDataset

!--------------------------------------------------------------------------------------------------
!> @brief writes to a new scalar dataset in the given group location
!--------------------------------------------------------------------------------------------------
subroutine HDF5_writeScalarDataset(group,dataset,label,SIunit,dataspace_size,mpiOffset)
 use hdf5

 implicit none
 integer(HID_T),   intent(in)               :: group
 character(len=*), intent(in)               :: SIunit,label
 integer(pInt),    intent(in)               :: dataspace_size, mpiOffset
 real(pReal),      intent(in), dimension(:) :: dataset

 integer        :: hdferr, nNodes
 integer(HID_T) :: dset_id, space_id, memspace, plist_id

 integer(HSIZE_T), dimension(1)             :: counter
 integer(HSIZE_T), dimension(1)             :: fileOffset

 nNodes = size(dataset)
 if (nNodes < 1) return

 call HDF5_addScalarDataset(group,dataspace_size,label,SIunit)
 call h5dopen_f(group, label, dset_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_writeScalarDataset: h5dopen_f')

 ! Define and select hyperslabs
 counter = size(dataset)                    ! how big i am
 fileOffset  = mpiOffset                    ! where i start to write my data

 call h5screate_simple_f(1, counter, memspace, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_writeScalarDataset: h5screate_simple_f')
 call h5dget_space_f(dset_id, space_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_writeScalarDataset: h5dget_space_f')
 call h5sselect_hyperslab_f(space_id, H5S_SELECT_SET_F, fileOffset, counter, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_writeScalarDataset: h5sselect_hyperslab_f')

 ! Create property list for collective dataset write
#ifdef PETSc
 call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_writeScalarDataset: h5pcreate_f')
 call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_writeScalarDataset: h5pset_dxpl_mpio_f')
#endif

 ! Write the dataset collectively
 call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, dataset, int([dataspace_size],HSIZE_T), hdferr, &
                 file_space_id = space_id, mem_space_id = memspace, xfer_prp = plist_id)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_writeScalarDataset: h5dwrite_f')

 call h5sclose_f(space_id, hdferr)
 call h5sclose_f(memspace, hdferr)
 call h5dclose_f(dset_id, hdferr)
 call h5pclose_f(plist_id, hdferr)

end subroutine HDF5_writeScalarDataset

!--------------------------------------------------------------------------------------------------
!> @brief subroutine for reading dataset of the type pReal with 1 dimension
!--------------------------------------------------------------------------------------------------
subroutine HDF5_read_pReal_1(dataset,loc_id,datasetName)

 implicit none
 real(pReal),      intent(out), dimension(:) ::    dataset
 integer(HID_T),   intent(in) :: loc_id                                                             !< file or group handle
 character(len=*), intent(in) :: datasetName                                                        !< name of the dataset in the file
 integer(pInt),dimension(:), allocatable  :: myShape
 integer        :: hdferr
 integer(HID_T) :: dset_id
 myShape = shape(dataset)

 call h5dopen_f(loc_id,datasetName,dset_id,hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pReal_shape1: h5dopen_f')
 call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,dataset,int(myShape,HSIZE_T),hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pReal_shape1: h5dread_f')
 call h5dclose_f(dset_id,hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pReal_shape1: h5dclose_f')

end subroutine HDF5_read_pReal_1

!--------------------------------------------------------------------------------------------------
!> @brief subroutine for reading dataset of the type pReal with 2 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_read_pReal_2(dataset,loc_id,datasetName)

 implicit none
 real(pReal),      intent(out), dimension(:,:) ::    dataset
 integer(HID_T),   intent(in) :: loc_id                                                             !< file or group handle
 character(len=*), intent(in) :: datasetName                                                        !< name of the dataset in the file
 integer(pInt),dimension(:), allocatable  :: myShape

 integer        :: hdferr
 integer(HID_T) :: dset_id
 myShape = shape(dataset)

 call h5dopen_f(loc_id,datasetName,dset_id,hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pReal_shape2: h5dopen_f')
 call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,dataset,int(myShape,HSIZE_T),hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pReal_shape2: h5dread_f')
 call h5dclose_f(dset_id,hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pReal_shape2: h5dclose_f')

end subroutine HDF5_read_pReal_2

!--------------------------------------------------------------------------------------------------
!> @brief subroutine for reading dataset of the type pReal with 3 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_read_pReal_3(dataset,loc_id,datasetName)

 implicit none
 real(pReal),      intent(out), dimension(:,:,:) ::    dataset
 integer(HID_T),   intent(in) :: loc_id                                                             !< file or group handle
 character(len=*), intent(in) :: datasetName                                                        !< name of the dataset in the file
 integer(pInt),dimension(:), allocatable  :: myShape
 integer        :: hdferr
 integer(HID_T) :: dset_id
 myShape = shape(dataset)

 call h5dopen_f(loc_id,datasetName,dset_id,hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pReal_shape3: h5dopen_f')
 call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,dataset,int(myShape,HSIZE_T),hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pReal_shape3: h5dread_f')
 call h5dclose_f(dset_id,hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pReal_shape3: h5dclose_f')

end subroutine HDF5_read_pReal_3

!--------------------------------------------------------------------------------------------------
!> @brief subroutine for reading dataset of the type pReal with 4 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_read_pReal_4(dataset,loc_id,datasetName)

 implicit none
 real(pReal),      intent(out), dimension(:,:,:,:) ::    dataset
 integer(HID_T),   intent(in) :: loc_id                                                             !< file or group handle
 character(len=*), intent(in) :: datasetName                                                        !< name of the dataset in the file
 integer(pInt),dimension(:), allocatable  :: myShape

 integer        :: hdferr
 integer(HID_T) :: dset_id
 myShape = shape(dataset)

 call h5dopen_f(loc_id,datasetName,dset_id,hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pReal_shape4: h5dopen_f')
 call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,dataset,int(myShape,HSIZE_T),hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pReal_shape4: h5dread_f')
 call h5dclose_f(dset_id,hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pReal_shape4: h5dclose_f')

end subroutine HDF5_read_pReal_4

!--------------------------------------------------------------------------------------------------
!> @brief subroutine for reading dataset of the type pReal with 5 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_read_pReal_5(dataset,loc_id,datasetName)

 implicit none
 real(pReal),      intent(out), dimension(:,:,:,:,:) ::    dataset
 integer(HID_T),   intent(in) :: loc_id                                                             !< file or group handle
 character(len=*), intent(in) :: datasetName                                                        !< name of the dataset in the file
 integer(pInt),dimension(:), allocatable  :: myShape

 integer        :: hdferr
 integer(HID_T) :: dset_id
 myShape = shape(dataset)

 call h5dopen_f(loc_id,datasetName,dset_id,hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pReal_shape5: h5dopen_f')
 call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,dataset,int(myShape,HSIZE_T),hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pReal_shape5: h5dread_f')
 call h5dclose_f(dset_id,hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pReal_shape5: h5dclose_f')

end subroutine HDF5_read_pReal_5

!--------------------------------------------------------------------------------------------------
!> @brief subroutine for reading dataset of the type pReal with 6 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_read_pReal_6(dataset,loc_id,datasetName)

 implicit none
 real(pReal),      intent(out), dimension(:,:,:,:,:,:) ::    dataset
 integer(HID_T),   intent(in) :: loc_id                                                             !< file or group handle
 character(len=*), intent(in) :: datasetName                                                        !< name of the dataset in the file
 integer(pInt),dimension(:), allocatable  :: myShape

 integer        :: hdferr
 integer(HID_T) :: dset_id
 myShape = shape(dataset)

 call h5dopen_f(loc_id,datasetName,dset_id,hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pReal_shape6: h5dopen_f')
 call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,dataset,int(myShape,HSIZE_T),hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pReal_shape6: h5dread_f')
 call h5dclose_f(dset_id,hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pReal_shape6: h5dclose_f')

end subroutine HDF5_read_pReal_6

!--------------------------------------------------------------------------------------------------
!> @brief subroutine for reading dataset of the type pReal with 7 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_read_pReal_7(dataset,loc_id,datasetName)

 implicit none
 real(pReal),      intent(out), dimension(:,:,:,:,:,:,:) ::    dataset
 integer(HID_T),   intent(in) :: loc_id                                                             !< file or group handle
 character(len=*), intent(in) :: datasetName                                                        !< name of the dataset in the file
 integer(pInt),dimension(:), allocatable  :: myShape

 integer        :: hdferr
 integer(HID_T) :: dset_id
 myShape = shape(dataset)

 call h5dopen_f(loc_id,datasetName,dset_id,hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pReal_shape7: h5dopen_f')
 call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,dataset,int(myShape,HSIZE_T),hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pReal_shape7: h5dread_f')
 call h5dclose_f(dset_id,hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pReal_shape7: h5dclose_f')

end subroutine HDF5_read_pReal_7

!--------------------------------------------------------------------------------------------------
!> @brief subroutine for reading dataset of the type pInt with 1 dimension
!--------------------------------------------------------------------------------------------------
subroutine HDF5_read_pInt_1(dataset,loc_id,datasetName)

 implicit none
 integer(pInt),      intent(out), dimension(:) ::    dataset
 integer(HID_T),   intent(in) :: loc_id                                                             !< file or group handle
 character(len=*), intent(in) :: datasetName                                                        !< name of the dataset in the file
 integer(pInt),dimension(:), allocatable  :: myShape

 integer        :: hdferr
 integer(HID_T) :: dset_id
 myShape = shape(dataset)

 call h5dopen_f(loc_id,datasetName,dset_id,hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pInt__shape1: h5dopen_f')
 call h5dread_f(dset_id,H5T_NATIVE_INTEGER,dataset,int(myShape,HSIZE_T),hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pInt__shape1: h5dread_f')
 call h5dclose_f(dset_id,hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pInt__shape1: h5dclose_f')

end subroutine HDF5_read_pInt_1

!--------------------------------------------------------------------------------------------------
!> @brief subroutine for reading dataset of the type pInt with 2 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_read_pInt_2(dataset,loc_id,datasetName)

 implicit none
 integer(pInt),      intent(out), dimension(:,:) ::    dataset
 integer(HID_T),   intent(in) :: loc_id                                                             !< file or group handle
 character(len=*), intent(in) :: datasetName                                                        !< name of the dataset in the file
 integer(pInt),dimension(:), allocatable  :: myShape

 integer        :: hdferr
 integer(HID_T) :: dset_id
 myShape = shape(dataset)

 call h5dopen_f(loc_id,datasetName,dset_id,hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pInt__shape2: h5dopen_f')
 call h5dread_f(dset_id,H5T_NATIVE_INTEGER,dataset,int(myShape,HSIZE_T),hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pInt__shape2: h5dread_f')
 call h5dclose_f(dset_id,hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pInt__shape2: h5dclose_f')

end subroutine HDF5_read_pInt_2

!--------------------------------------------------------------------------------------------------
!> @brief subroutine for reading dataset of the type pInt with 3 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_read_pInt_3(dataset,loc_id,datasetName)

 implicit none
 integer(pInt),      intent(out), dimension(:,:,:) ::    dataset
 integer(HID_T),   intent(in) :: loc_id                                                             !< file or group handle
 character(len=*), intent(in) :: datasetName                                                        !< name of the dataset in the file
 integer(pInt),dimension(:), allocatable  :: myShape

 integer        :: hdferr
 integer(HID_T) :: dset_id
 myShape = shape(dataset)

 call h5dopen_f(loc_id,datasetName,dset_id,hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pInt__shape3: h5dopen_f')
 call h5dread_f(dset_id,H5T_NATIVE_INTEGER,dataset,int(myShape,HSIZE_T),hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pInt__shape3: h5dread_f')
 call h5dclose_f(dset_id,hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pInt__shape3: h5dclose_f')

end subroutine HDF5_read_pInt_3

!--------------------------------------------------------------------------------------------------
!> @brief subroutine for reading dataset of the type pInt with 4 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_read_pInt_4(dataset,loc_id,datasetName)

 implicit none
 integer(pInt),      intent(out), dimension(:,:,:,:) ::    dataset
 integer(HID_T),   intent(in) :: loc_id                                                             !< file or group handle
 character(len=*), intent(in) :: datasetName                                                        !< name of the dataset in the file
 integer(pInt),dimension(:), allocatable  :: myShape

 integer        :: hdferr
 integer(HID_T) :: dset_id
 myShape = shape(dataset)

 call h5dopen_f(loc_id,datasetName,dset_id,hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pInt__shape4: h5dopen_f')
 call h5dread_f(dset_id,H5T_NATIVE_INTEGER,dataset,int(myShape,HSIZE_T),hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pInt__shape4: h5dread_f')
 call h5dclose_f(dset_id,hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pInt__shape4: h5dclose_f')

end subroutine HDF5_read_pInt_4

!--------------------------------------------------------------------------------------------------
!> @brief subroutine for reading dataset of the type pInt with 5 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_read_pInt_5(dataset,loc_id,datasetName)

 implicit none
 integer(pInt),      intent(out), dimension(:,:,:,:,:) ::    dataset
 integer(HID_T),   intent(in) :: loc_id                                                             !< file or group handle
 character(len=*), intent(in) :: datasetName                                                        !< name of the dataset in the file
 integer(pInt),dimension(:), allocatable  :: myShape

 integer        :: hdferr
 integer(HID_T) :: dset_id
 myShape = shape(dataset)

 call h5dopen_f(loc_id,datasetName,dset_id,hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pInt__shape5: h5dopen_f')
 call h5dread_f(dset_id,H5T_NATIVE_INTEGER,dataset,int(myShape,HSIZE_T),hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pInt__shape5: h5dread_f')
 call h5dclose_f(dset_id,hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pInt__shape5: h5dclose_f')

end subroutine HDF5_read_pInt_5

!--------------------------------------------------------------------------------------------------
!> @brief subroutine for reading dataset of the type pInt with 6 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_read_pInt_6(dataset,loc_id,datasetName)

 implicit none
 integer(pInt),      intent(out), dimension(:,:,:,:,:,:) ::    dataset
 integer(HID_T),   intent(in) :: loc_id                                                             !< file or group handle
 character(len=*), intent(in) :: datasetName                                                        !< name of the dataset in the file
 integer(pInt),dimension(:), allocatable  :: myShape

 integer        :: hdferr
 integer(HID_T) :: dset_id
 myShape = shape(dataset)

 call h5dopen_f(loc_id,datasetName,dset_id,hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pInt__shape6: h5dopen_f')
 call h5dread_f(dset_id,H5T_NATIVE_INTEGER,dataset,int(myShape,HSIZE_T),hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pInt__shape6: h5dread_f')
 call h5dclose_f(dset_id,hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pInt__shape6: h5dclose_f')

end subroutine HDF5_read_pInt_6

!--------------------------------------------------------------------------------------------------
!> @brief subroutine for reading dataset of the type pInt with 7 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_read_pInt_7(dataset,loc_id,datasetName)

 implicit none
 integer(pInt),      intent(out), dimension(:,:,:,:,:,:,:) ::    dataset
 integer(HID_T),   intent(in) :: loc_id                                                             !< file or group handle
 character(len=*), intent(in) :: datasetName                                                        !< name of the dataset in the file
 integer(pInt),dimension(:), allocatable  :: myShape

 integer        :: hdferr
 integer(HID_T) :: dset_id
 myShape = shape(dataset)

 call h5dopen_f(loc_id,datasetName,dset_id,hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pInt__shape7: h5dopen_f')
 call h5dread_f(dset_id,H5T_NATIVE_INTEGER,dataset,int(myShape,HSIZE_T),hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pInt__shape7: h5dread_f')
 call h5dclose_f(dset_id,hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pInt__shape7: h5dclose_f')

end subroutine HDF5_read_pInt_7

!--------------------------------------------------------------------------------------------------
!> @brief subroutine for writing dataset of the type pReal with 1 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_pReal1(dataset,loc_id,datasetName)

 implicit none
 real(pReal),      intent(out), dimension(:) ::    dataset
 integer(HID_T),   intent(in) :: loc_id                                                             !< file or group handle
 character(len=*), intent(in) :: datasetName                                                        !< name of the dataset in the file

 integer(pInt), dimension(:), allocatable :: myShape                                                          !<shape of the dataset
 integer :: hdferr
 integer(HID_T) :: dset_id, space_id

 myShape = shape(dataset)
!--------------------------------------------------------------------------------------------------
! create dataspace
call h5screate_simple_f(size(myShape), int(myShape,HSIZE_T), space_id, hdferr, &
                           int(myShape,HSIZE_T))
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pReal1: h5screate_simple_f')
!--------------------------------------------------------------------------------------------------
! create dataset
call h5dcreate_f(loc_id, trim(datasetName), H5T_NATIVE_DOUBLE, space_id, dset_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pReal1: h5dcreate_f')

CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,dataset,int(myShape,HSIZE_T), hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pReal1: h5dcreate_f')

!--------------------------------------------------------------------------------------------------
!close types, dataspaces
 call h5dclose_f(dset_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_HDF5_write_pReal1: h5dclose_f')
call h5sclose_f(space_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_HDF5_write_pReal1: h5sclose_f')

end subroutine HDF5_write_pReal1

!--------------------------------------------------------------------------------------------------
!> @brief subroutine for writing dataset of the type pReal with 2 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_pReal2(dataset,loc_id,datasetName)

 implicit none
 real(pReal),      intent(out), dimension(:,:) ::    dataset
 integer(HID_T),   intent(in) :: loc_id                                                             !< file or group handle
 character(len=*), intent(in) :: datasetName                                                        !< name of the dataset in the file

 integer(pInt), dimension(:), allocatable :: myShape                                                          !<shape of the dataset
 integer :: hdferr
 integer(HID_T) :: dset_id, space_id

 myShape = shape(dataset)
!--------------------------------------------------------------------------------------------------
! create dataspace
call h5screate_simple_f(size(myShape), int(myShape,HSIZE_T), space_id, hdferr, &
                           int(myShape,HSIZE_T))
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pReal2: h5screate_simple_f')
!--------------------------------------------------------------------------------------------------
! create dataset
call h5dcreate_f(loc_id, trim(datasetName), H5T_NATIVE_DOUBLE, space_id, dset_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pReal2: h5dcreate_f')

CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,dataset,int(myShape,HSIZE_T), hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pReal2: h5dcreate_f')

!--------------------------------------------------------------------------------------------------
!close types, dataspaces
 call h5dclose_f(dset_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_HDF5_write_pReal2: h5dclose_f')
call h5sclose_f(space_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_HDF5_write_pReal2: h5sclose_f')

end subroutine HDF5_write_pReal2

!--------------------------------------------------------------------------------------------------
!> @brief subroutine for writing dataset of the type pReal with 3 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_pReal3(dataset,loc_id,datasetName)

 implicit none
 real(pReal),      intent(out), dimension(:,:,:) ::    dataset
 integer(HID_T),   intent(in) :: loc_id                                                             !< file or group handle
 character(len=*), intent(in) :: datasetName                                                        !< name of the dataset in the file

 integer(pInt), dimension(:), allocatable :: myShape                                                          !<shape of the dataset
 integer :: hdferr
 integer(HID_T) :: dset_id, space_id

 myShape = shape(dataset)
!--------------------------------------------------------------------------------------------------
! create dataspace
call h5screate_simple_f(size(myShape), int(myShape,HSIZE_T), space_id, hdferr, &
                           int(myShape,HSIZE_T))
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pReal3: h5screate_simple_f')
!--------------------------------------------------------------------------------------------------
! create dataset
call h5dcreate_f(loc_id, trim(datasetName), H5T_NATIVE_DOUBLE, space_id, dset_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pReal3: h5dcreate_f')

CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,dataset,int(myShape,HSIZE_T), hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pReal3: h5dcreate_f')

!--------------------------------------------------------------------------------------------------
!close types, dataspaces
 call h5dclose_f(dset_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_HDF5_write_pReal3: h5dclose_f')
call h5sclose_f(space_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_HDF5_write_pReal3: h5sclose_f')

end subroutine HDF5_write_pReal3

!--------------------------------------------------------------------------------------------------
!> @brief subroutine for writing dataset of the type pReal with 4 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_pReal4(dataset,loc_id,datasetName)

 implicit none
 real(pReal),      intent(out), dimension(:,:,:,:) ::    dataset
 integer(HID_T),   intent(in) :: loc_id                                                             !< file or group handle
 character(len=*), intent(in) :: datasetName                                                        !< name of the dataset in the file

 integer(pInt), dimension(:), allocatable :: myShape                                                          !<shape of the dataset
 integer :: hdferr
 integer(HID_T) :: dset_id, space_id

 myShape = shape(dataset)
!--------------------------------------------------------------------------------------------------
! create dataspace
call h5screate_simple_f(size(myShape), int(myShape,HSIZE_T), space_id, hdferr, &
                           int(myShape,HSIZE_T))
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pReal4: h5screate_simple_f')
!--------------------------------------------------------------------------------------------------
! create dataset
call h5dcreate_f(loc_id, trim(datasetName), H5T_NATIVE_DOUBLE, space_id, dset_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pReal4: h5dcreate_f')

CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,dataset,int(myShape,HSIZE_T), hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pReal4: h5dcreate_f')

!--------------------------------------------------------------------------------------------------
!close types, dataspaces
 call h5dclose_f(dset_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_HDF5_write_pReal4: h5dclose_f')
call h5sclose_f(space_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_HDF5_write_pReal4: h5sclose_f')

end subroutine HDF5_write_pReal4

!--------------------------------------------------------------------------------------------------
!> @brief subroutine for writing dataset of the type pReal with 5 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_pReal5(dataset,loc_id,datasetName)

 implicit none
 real(pReal),      intent(out), dimension(:,:,:,:,:) ::    dataset
 integer(HID_T),   intent(in) :: loc_id                                                             !< file or group handle
 character(len=*), intent(in) :: datasetName                                                        !< name of the dataset in the file

 integer(pInt), dimension(:), allocatable :: myShape                                                          !<shape of the dataset
 integer :: hdferr
 integer(HID_T) :: dset_id, space_id

 myShape = shape(dataset)
!--------------------------------------------------------------------------------------------------
! create dataspace
call h5screate_simple_f(size(myShape), int(myShape,HSIZE_T), space_id, hdferr, &
                           int(myShape,HSIZE_T))
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pReal5: h5screate_simple_f')
!--------------------------------------------------------------------------------------------------
! create dataset
call h5dcreate_f(loc_id, trim(datasetName), H5T_NATIVE_DOUBLE, space_id, dset_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pReal5: h5dcreate_f')

CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,dataset,int(myShape,HSIZE_T), hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pReal5: h5dcreate_f')

!--------------------------------------------------------------------------------------------------
!close types, dataspaces
 call h5dclose_f(dset_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_HDF5_write_pReal5: h5dclose_f')
call h5sclose_f(space_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_HDF5_write_pReal5: h5sclose_f')

end subroutine HDF5_write_pReal5

!--------------------------------------------------------------------------------------------------
!> @brief subroutine for writing dataset of the type pReal with 6 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_pReal6(dataset,loc_id,datasetName)

 implicit none
 real(pReal),      intent(out), dimension(:,:,:,:,:,:) ::    dataset
 integer(HID_T),   intent(in) :: loc_id                                                             !< file or group handle
 character(len=*), intent(in) :: datasetName                                                        !< name of the dataset in the file

 integer(pInt), dimension(:), allocatable :: myShape                                                          !<shape of the dataset
 integer :: hdferr
 integer(HID_T) :: dset_id, space_id

 myShape = shape(dataset)
!--------------------------------------------------------------------------------------------------
! create dataspace
call h5screate_simple_f(size(myShape), int(myShape,HSIZE_T), space_id, hdferr, &
                           int(myShape,HSIZE_T))
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pReal6: h5screate_simple_f')
!--------------------------------------------------------------------------------------------------
! create dataset
call h5dcreate_f(loc_id, trim(datasetName), H5T_NATIVE_DOUBLE, space_id, dset_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pReal6: h5dcreate_f')

CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,dataset,int(myShape,HSIZE_T), hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pReal6: h5dcreate_f')

!--------------------------------------------------------------------------------------------------
!close types, dataspaces
 call h5dclose_f(dset_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_HDF5_write_pReal6: h5dclose_f')
call h5sclose_f(space_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_HDF5_write_pReal6: h5sclose_f')

end subroutine HDF5_write_pReal6

!--------------------------------------------------------------------------------------------------
!> @brief subroutine for writing dataset of the type pReal with 7 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_pReal7(dataset,loc_id,datasetName)

 implicit none
 real(pReal),      intent(out), dimension(:,:,:,:,:,:,:) ::    dataset
 integer(HID_T),   intent(in) :: loc_id                                                             !< file or group handle
 character(len=*), intent(in) :: datasetName                                                        !< name of the dataset in the file

 integer(pInt), dimension(:), allocatable :: myShape                                                          !<shape of the dataset
 integer :: hdferr
 integer(HID_T) :: dset_id, space_id

 myShape = shape(dataset)
!--------------------------------------------------------------------------------------------------
! create dataspace
call h5screate_simple_f(size(myShape), int(myShape,HSIZE_T), space_id, hdferr, &
                           int(myShape,HSIZE_T))
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pReal7: h5screate_simple_f')
!--------------------------------------------------------------------------------------------------
! create dataset
call h5dcreate_f(loc_id, trim(datasetName), H5T_NATIVE_DOUBLE, space_id, dset_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pReal7: h5dcreate_f')

CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,dataset,int(myShape,HSIZE_T), hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pReal7: h5dcreate_f')

!--------------------------------------------------------------------------------------------------
!close types, dataspaces
 call h5dclose_f(dset_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_HDF5_write_pReal7: h5dclose_f')
call h5sclose_f(space_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_HDF5_write_pReal7: h5sclose_f')

end subroutine HDF5_write_pReal7

!--------------------------------------------------------------------------------------------------
!> @brief subroutine for writing dataset of the type pInt with 1 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_pInt1(dataset,loc_id,datasetName)

 implicit none
 integer(pInt),      intent(out), dimension(:) ::    dataset
 integer(HID_T),   intent(in) :: loc_id                                                             !< file or group handle
 character(len=*), intent(in) :: datasetName                                                        !< name of the dataset in the file

 integer(pInt), dimension(:), allocatable :: myShape                                                          !<shape of the dataset
 integer :: hdferr
 integer(HID_T) :: dset_id, space_id

 myShape = shape(dataset)
!--------------------------------------------------------------------------------------------------
! create dataspace
call h5screate_simple_f(size(myShape), int(myShape,HSIZE_T), space_id, hdferr, &
                           int(myShape,HSIZE_T))
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pInt1: h5screate_simple_f')
!--------------------------------------------------------------------------------------------------
! create dataset
call h5dcreate_f(loc_id, trim(datasetName), H5T_NATIVE_INTEGER, space_id, dset_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pInt1: h5dcreate_f')

CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER,dataset,int(myShape,HSIZE_T), hdferr)

!--------------------------------------------------------------------------------------------------
!close types, dataspaces
 call h5dclose_f(dset_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_HDF5_write_pInt1: h5dclose_f')
call h5sclose_f(space_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_HDF5_write_pInt1: h5sclose_f')

end subroutine HDF5_write_pInt1

!--------------------------------------------------------------------------------------------------
!> @brief subroutine for writing dataset of the type pInt with 2 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_pInt2(dataset,loc_id,datasetName)

 implicit none
 integer(pInt),      intent(out), dimension(:,:) ::    dataset
 integer(HID_T),   intent(in) :: loc_id                                                             !< file or group handle
 character(len=*), intent(in) :: datasetName                                                        !< name of the dataset in the file

 integer(pInt), dimension(:), allocatable :: myShape                                                          !<shape of the dataset
 integer :: hdferr
 integer(HID_T) :: dset_id, space_id

 myShape = shape(dataset)
!--------------------------------------------------------------------------------------------------
! create dataspace
call h5screate_simple_f(size(myShape), int(myShape,HSIZE_T), space_id, hdferr, &
                           int(myShape,HSIZE_T))
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pInt2: h5screate_simple_f')
!--------------------------------------------------------------------------------------------------
! create dataset
call h5dcreate_f(loc_id, trim(datasetName), H5T_NATIVE_INTEGER, space_id, dset_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pInt2: h5dcreate_f')

CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER,dataset,int(myShape,HSIZE_T), hdferr)

!--------------------------------------------------------------------------------------------------
!close types, dataspaces
 call h5dclose_f(dset_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_HDF5_write_pInt2: h5dclose_f')
call h5sclose_f(space_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_HDF5_write_pInt2: h5sclose_f')

end subroutine HDF5_write_pInt2

!--------------------------------------------------------------------------------------------------
!> @brief subroutine for writing dataset of the type pInt with 3 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_pInt3(dataset,loc_id,datasetName)

 implicit none
 integer(pInt),      intent(out), dimension(:,:,:) ::    dataset
 integer(HID_T),   intent(in) :: loc_id                                                             !< file or group handle
 character(len=*), intent(in) :: datasetName                                                        !< name of the dataset in the file

 integer(pInt), dimension(:), allocatable :: myShape                                                          !<shape of the dataset
 integer :: hdferr
 integer(HID_T) :: dset_id, space_id

 myShape = shape(dataset)
!--------------------------------------------------------------------------------------------------
! create dataspace
call h5screate_simple_f(size(myShape), int(myShape,HSIZE_T), space_id, hdferr, &
                           int(myShape,HSIZE_T))
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pInt3: h5screate_simple_f')
!--------------------------------------------------------------------------------------------------
! create dataset
call h5dcreate_f(loc_id, trim(datasetName), H5T_NATIVE_INTEGER, space_id, dset_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pInt3: h5dcreate_f')

CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER,dataset,int(myShape,HSIZE_T), hdferr)

!--------------------------------------------------------------------------------------------------
!close types, dataspaces
 call h5dclose_f(dset_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_HDF5_write_pInt3: h5dclose_f')
call h5sclose_f(space_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_HDF5_write_pInt3: h5sclose_f')

end subroutine HDF5_write_pInt3

!--------------------------------------------------------------------------------------------------
!> @brief subroutine for writing dataset of the type pInt with 4 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_pInt4(dataset,loc_id,datasetName)

 implicit none
 integer(pInt),      intent(out), dimension(:,:,:,:) ::    dataset
 integer(HID_T),   intent(in) :: loc_id                                                             !< file or group handle
 character(len=*), intent(in) :: datasetName                                                        !< name of the dataset in the file

 integer(pInt), dimension(:), allocatable :: myShape                                                          !<shape of the dataset
 integer :: hdferr
 integer(HID_T) :: dset_id, space_id

 myShape = shape(dataset)
!--------------------------------------------------------------------------------------------------
! create dataspace
call h5screate_simple_f(size(myShape), int(myShape,HSIZE_T), space_id, hdferr, &
                           int(myShape,HSIZE_T))
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pInt4: h5screate_simple_f')
!--------------------------------------------------------------------------------------------------
! create dataset
call h5dcreate_f(loc_id, trim(datasetName), H5T_NATIVE_INTEGER, space_id, dset_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pInt4: h5dcreate_f')

CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER,dataset,int(myShape,HSIZE_T), hdferr)

!--------------------------------------------------------------------------------------------------
!close types, dataspaces
 call h5dclose_f(dset_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_HDF5_write_pInt4: h5dclose_f')
call h5sclose_f(space_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_HDF5_write_pInt4: h5sclose_f')

end subroutine HDF5_write_pInt4

!--------------------------------------------------------------------------------------------------
!> @brief subroutine for writing dataset of the type pInt with 5 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_pInt5(dataset,loc_id,datasetName)

 implicit none
 integer(pInt),      intent(out), dimension(:,:,:,:,:) ::    dataset
 integer(HID_T),   intent(in) :: loc_id                                                             !< file or group handle
 character(len=*), intent(in) :: datasetName                                                        !< name of the dataset in the file

 integer(pInt), dimension(:), allocatable :: myShape                                                          !<shape of the dataset
 integer :: hdferr
 integer(HID_T) :: dset_id, space_id

 myShape = shape(dataset)
!--------------------------------------------------------------------------------------------------
! create dataspace
call h5screate_simple_f(size(myShape), int(myShape,HSIZE_T), space_id, hdferr, &
                           int(myShape,HSIZE_T))
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pInt5: h5screate_simple_f')
!--------------------------------------------------------------------------------------------------
! create dataset
call h5dcreate_f(loc_id, trim(datasetName), H5T_NATIVE_INTEGER, space_id, dset_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pInt5: h5dcreate_f')

CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER,dataset,int(myShape,HSIZE_T), hdferr)

!--------------------------------------------------------------------------------------------------
!close types, dataspaces
 call h5dclose_f(dset_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_HDF5_write_pInt5: h5dclose_f')
call h5sclose_f(space_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_HDF5_write_pInt5: h5sclose_f')

end subroutine HDF5_write_pInt5

!--------------------------------------------------------------------------------------------------
!> @brief subroutine for writing dataset of the type pInt with 6 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_pInt6(dataset,loc_id,datasetName)

 implicit none
 integer(pInt),      intent(out), dimension(:,:,:,:,:,:) ::    dataset
 integer(HID_T),   intent(in) :: loc_id                                                             !< file or group handle
 character(len=*), intent(in) :: datasetName                                                        !< name of the dataset in the file

 integer(pInt), dimension(:), allocatable :: myShape                                                          !<shape of the dataset
 integer :: hdferr
 integer(HID_T) :: dset_id, space_id

 myShape = shape(dataset)
!--------------------------------------------------------------------------------------------------
! create dataspace
call h5screate_simple_f(size(myShape), int(myShape,HSIZE_T), space_id, hdferr, &
                           int(myShape,HSIZE_T))
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pInt6: h5screate_simple_f')
!--------------------------------------------------------------------------------------------------
! create dataset
call h5dcreate_f(loc_id, trim(datasetName), H5T_NATIVE_INTEGER, space_id, dset_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pInt6: h5dcreate_f')

CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER,dataset,int(myShape,HSIZE_T), hdferr)

!--------------------------------------------------------------------------------------------------
!close types, dataspaces
 call h5dclose_f(dset_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_HDF5_write_pInt6: h5dclose_f')
call h5sclose_f(space_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_HDF5_write_pInt6: h5sclose_f')

end subroutine HDF5_write_pInt6

!--------------------------------------------------------------------------------------------------
!> @brief subroutine for writing dataset of the type pInt with 7 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_pInt7(dataset,loc_id,datasetName)

 implicit none
 integer(pInt),      intent(out), dimension(:,:,:,:,:,:,:) ::    dataset
 integer(HID_T),   intent(in) :: loc_id                                                             !< file or group handle
 character(len=*), intent(in) :: datasetName                                                        !< name of the dataset in the file

 integer(pInt), dimension(:), allocatable :: myShape                                                          !<shape of the dataset
 integer :: hdferr
 integer(HID_T) :: dset_id, space_id

 myShape = shape(dataset)
!--------------------------------------------------------------------------------------------------
! create dataspace
call h5screate_simple_f(size(myShape), int(myShape,HSIZE_T), space_id, hdferr, &
                           int(myShape,HSIZE_T))
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pInt7: h5screate_simple_f')
!--------------------------------------------------------------------------------------------------
! create dataset
call h5dcreate_f(loc_id, trim(datasetName), H5T_NATIVE_INTEGER, space_id, dset_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pInt7: h5dcreate_f')

CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER,dataset,int(myShape,HSIZE_T), hdferr)

!--------------------------------------------------------------------------------------------------
!close types, dataspaces
 call h5dclose_f(dset_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_HDF5_write_pInt7: h5dclose_f')
call h5sclose_f(space_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_HDF5_write_pInt7: h5sclose_f')

end subroutine HDF5_write_pInt7

!--------------------------------------------------------------------------------------------------
!> @brief adds a new scalar dataset to the given group location
!--------------------------------------------------------------------------------------------------
subroutine HDF5_addScalarDataset(group,nnodes,label,SIunit)
 use hdf5

 implicit none
 integer(HID_T),   intent(in) :: group
 integer(pInt),    intent(in) :: nnodes
 character(len=*), intent(in) :: SIunit,label

 integer        :: hdferr
 integer(HID_T) :: space_id, dset_id

!--------------------------------------------------------------------------------------------------
! create dataspace
 call h5screate_simple_f(1, int([Nnodes],HSIZE_T), space_id, hdferr, &
                            int([Nnodes],HSIZE_T))
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_addScalarDataset: h5screate_simple_f')

!--------------------------------------------------------------------------------------------------
! create Dataset
 call h5dcreate_f(group, trim(label),H5T_NATIVE_DOUBLE, space_id, dset_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_addScalarDataset: h5dcreate_f')
 call HDF5_addStringAttribute(dset_id,'unit',trim(SIunit))

!--------------------------------------------------------------------------------------------------
!close types, dataspaces
 call h5dclose_f(dset_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_addScalarDataset: h5dclose_f')
 call h5sclose_f(space_id, hdferr)

end subroutine HDF5_addScalarDataset

!--------------------------------------------------------------------------------------------------
!> @brief copies the current temp results to the actual results file
!--------------------------------------------------------------------------------------------------
subroutine HDF5_forwardResults(time)
 use hdf5
 use IO, only: &
   IO_intOut

 implicit none
 integer :: hdferr
 integer(HID_T) :: currentIncID
 real(pReal), intent(in) :: time
 character(len=1024)     :: myName

 currentInc = currentInc +1_pInt
 write(6,*) 'forward results';flush(6)
 write(myName,'(a,'//IO_intOut(currentInc)//')') 'inc',currentInc
 currentIncID = HDF5_addGroup(myName)
 call HDF5_setLink(myName,'current')
! call HDF5_flush(resultsFile)
 call HDF5_closeGroup(currentIncID)

end subroutine HDF5_forwardResults

end module HDF5_Utilities