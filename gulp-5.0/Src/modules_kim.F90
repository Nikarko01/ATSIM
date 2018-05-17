!************************************
!  Module for GULP-KIM interaction  *
!************************************
!
!  10/12 Created from modules.f90
!   4/14 Modified for OpenKIM v1.4
!   7/14 pkim_model now kept as integer pointer
!   7/16 pkim_model changed to an array of type c_ptr
!
!  Conditions of use:
!
!  GULP is available free of charge to academic institutions
!  and non-commerical establishments only. Copies should be
!  obtained from the author only and should not be distributed
!  in any form by the user to a third party without the express
!  permission of the author. This notice applies to all parts
!  of the program, except any library routines which are
!  distributed with the code for completeness. All rights for
!  such routines remain with the original distributor.
!
!  No claim is made that this program is free from errors and
!  no liability will be accepted for any loss or damage that
!  may result. The user is responsible for checking the validity
!  of their results.
!
!  Copyright Curtin University 2016
!
!  Julian Gale, CIC, Curtin University, July 2016
!
!  NB: In this module the precision of some quantities is not set using the GULP defined types
!      since they must remain compatable with OpenKIM
!

!
!  OpenKIM models
!
  module kim_models
    use datatypes
#if KIM_F03
    use KIM_API_F03
    use, intrinsic :: iso_c_binding 
#elif KIM
    use KIM_API
#endif
    integer(i4),                                parameter         :: kim_len = 10000               ! Length of KIM descriptor string
    integer(i4),                                             save :: maxkimmodel = 1               ! Maximum number of KIM models
    character(len=80),        dimension(:),     pointer,     save :: kim_model => null()           ! Names of KIM models
#ifdef KIM_F03
    type(c_ptr),              dimension(:,:),   allocatable, save :: pkim_model                    ! C pointers for each KIM model/configuration
    integer(c_int),                             target,      save :: kim_numat                     ! Integer containing numat
    integer(c_int),                             pointer,     save :: pkim_numat => null()          ! Integer pointer for numat
    integer(c_int),           dimension(:,:),   pointer,     save :: kim_nbc => null()             ! Array containing integer code NBC setting : 0 => cluster, 1=> RVEC_F
    integer(c_int),           dimension(:),     pointer,     save :: kim_ntype => null()           ! Array to hold number of particle types for KIM for each configuration
    integer(c_int),           dimension(:,:),   pointer,     save :: kim_types => null()           ! Array to hold particle types for KIM of each atom in each configuration
    real(c_double),           dimension(:,:),   pointer,     save :: kim_coord => null()           ! Array to hold coordinates for KIM
    real(c_double),           dimension(:),     pointer,     save :: kim_cutoff => null()          ! Array to hold cutoffs for KIM models
    real(c_double),           dimension(:,:),   pointer,     save :: kim_forces => null()          ! Array to hold forces for KIM
    real(c_double),                             target,      save :: kim_energy                    ! Variable to hold KIM energy
    real(c_double),                             pointer,     save :: pkim_energy => null()         ! Pointer to KIM energy
    real(c_double),           dimension(:),     pointer,     save :: kim_virial => null()          ! Variable to hold KIM virial
#elif KIM
    integer(kind=kim_intptr), dimension(:,:),   pointer,     save :: pkim_model => null()          ! Integer references for KIM pointers for each model/configuration
    integer,                  dimension(:,:),   pointer,     save :: kim_nbc => null()             ! Array containing integer code NBC setting : 0 => cluster, 1=> RVEC_F
    integer,                  dimension(:),     pointer,     save :: kim_ntype => null()           ! Array to hold number of particle types for KIM for each configuration
    integer,                  dimension(:,:),   pointer,     save :: kim_types => null()           ! Array to hold particle types for KIM of each atom in each configuration
    real*8,                   dimension(:,:),   pointer,     save :: kim_coord => null()           ! Array to hold coordinates for KIM
    real*8,                   dimension(:),     pointer,     save :: kim_cutoff => null()          ! Array to hold cutoffs for KIM models
    real*8,                   dimension(:,:),   pointer,     save :: kim_forces => null()          ! Array to hold forces for KIM
    real*8,                                                  save :: kim_energy                    ! Variable to hold KIM energy
    real*8,                                                  save :: kim_virial(6)                 ! Variable to hold KIM virial
#endif
    logical,                  dimension(:,:),   pointer,     save :: lkim_model_cfg_OK => null()   ! Stores logical check that configuration is OK for KIM model
!
    character(len=kim_len),   dimension(:,:),   pointer,     save :: kim_test_descriptor => null() ! KIM descriptor for test - used instead of a .kim file
    integer(i4),                                             save :: nkimmodel = 0                 ! Number of KIM models
    integer(i4),                                             save :: nlibnkimmodel = 0             ! Number of KIM models prior to libraries
    integer(i4),                                             save :: kim_types_dummy(1)            ! Dummy variable for KIM particle types
    logical,                                                 save :: lkim_model = .false.          ! If true then OpenKIM model to be used
  end module kim_models
