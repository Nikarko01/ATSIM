  subroutine setkim
!
!  Initialises OpenKIM models for each configuration 
!
!  10/12 Created
!   7/14 Modifications for OpenKIM v1.4 started
!   7/16 Modifications for OpenKIM v1.7.3 added
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
!  Julian Gale, Curtin University, July 2016
!
#if KIM_F03
#include "KIM_API_status.h"
#define TRUEFALSE(TRUTH) merge(1,0,(TRUTH))
  use, intrinsic :: iso_c_binding
  use KIM_API_F03
#elif KIM
#include "KIM_API_status.h"
#define TRUEFALSE(TRUTH) merge(1,0,(TRUTH))
  use KIM_API
#endif
  use configurations, only : ncfg
  use control
  use current
  use iochannels
#if KIM_F03
  use kim_functions,  only : get_kim_neigh, kim_neighobject
  use reallocate
#elif KIM
  use kim_functions,  only : get_kim_neigh, kim_neighobject
#endif
  use kim_models
  use parallel
  use species,        only : maxspec
  implicit none
!
!  Local variables
!
  character(len=5)                                  :: type_string      ! Element symbol for KIM atom typing
  character(len=2), allocatable, dimension(:), save :: kim_ntype_symbol ! Element symbol for KIM types
  integer(i4)                                       :: i                ! Loop counter over number of atoms
  integer(i4)                                       :: j                ! Loop counter over number of atoms
  integer(i4)                                       :: icfg             ! Loop counter over configurations
  integer(i4)                                       :: nkim             ! Loop counter over KIM models
  logical                                           :: lfound           ! Logical to indicate whether a match between types has been found
#ifdef KIM_F03
  integer(i4)                                       :: ierror           ! Error flag
  integer(c_int)                                    :: int_logic        ! Integer used as a logical
  integer(c_int)                                    :: KIMerror         ! KIMerror flag
  integer(c_int)                                    :: KIMerror_d       ! KIMerror flag - dummy
!
!  Integer local variables - use c_int for size
!
  integer(c_int)                                    :: nsize_1          ! Variable that contains the value one in C integer type
  integer(c_int)                                    :: nsize_6          ! Variable that contains the value six in C integer type
  integer(c_int)                                    :: nsize_N          ! Variable that contains the value of numat in C integer type
  integer(c_int)                                    :: nsize_3N         ! Variable that contains the value of 3*numat in C integer type
#elif KIM
  integer(i4)                                       :: int_logic        ! Integer used as a logical
  integer(i4)                                       :: KIMerror         ! KIMerror flag
  integer(i4)                                       :: KIMerror_d       ! KIMerror flag - dummy
!
!  Integer local variables - the size of these values is fixed by KIM
!
  integer(kind=kim_intptr)                          :: nsize_1          ! Variable that contains the value one in KIM integer type
  integer(kind=kim_intptr)                          :: nsize_6          ! Variable that contains the value six in KIM integer type
  integer(kind=kim_intptr)                          :: nsize_N          ! Variable that contains the value of numat in KIM integer type
  integer(kind=kim_intptr)                          :: nsize_3N         ! Variable that contains the value of 3*numat in KIM integer type
  integer(kind=kim_intptr)                          :: loc_numat        ! Location of numat in memory
  integer(kind=kim_intptr)                          :: loc_ntype        ! Location of kim_ntype in memory
  integer(kind=kim_intptr)                          :: loc_types        ! Location of kim_types in memory
  integer(kind=kim_intptr)                          :: loc_coords       ! Location of kim_coord in memory
  integer(kind=kim_intptr)                          :: loc_cutoff       ! Location of kim_cutoff in memory
  integer(kind=kim_intptr)                          :: loc_energy       ! Location of kim_energy in memory
  integer(kind=kim_intptr)                          :: loc_forces       ! Location of kim_forces in memory
  integer(kind=kim_intptr)                          :: loc_virial       ! Location of kim_virial in memory
#endif
!
!  Allocate workspace array for type symbols
!
  allocate(kim_ntype_symbol(maxspec))
#ifdef KIM_F03
!
!  Allocate kim related pointers
!
  call realloc(kim_neighobject,1_i4,ierror)
  if (ierror.ne.0) call outofmemory('setkim','kim_neighobject')
  kim_neighobject = 0
  call realloc(kim_virial,6_i4,ierror)
  if (ierror.ne.0) call outofmemory('setkim','kim_virial')
!
!  Set kim pointers for scalars to targets
!
  pkim_numat => kim_numat
  pkim_energy => kim_energy
#endif
!
!  Loop over configurations
!
  do icfg = 1,ncfg
!
!  Setup for this configuration
!
    ncf = icfg
    call setup(.true.)
!
!  Determine the number of species that will be passed to KIM
!
!  NB: At this stage KIM cannot handle multiple types of the same element
!
    call label(nat(1),0_i4,type_string)
#if KIM_F03
    kim_ntype(icfg) = 1
    kim_ntype_symbol(1) = type_string(1:2)
#elif KIM
    kim_ntype(icfg) = 1
    kim_ntype_symbol(1) = type_string(1:2)
#endif
    do i = 2,numat
      j = 0
      lfound = .false.
      do while (.not.lfound.and.j.lt.i-1)
        j = j + 1
        lfound = (nat(i).eq.nat(j))
      enddo
      if (.not.lfound) then
        call label(nat(i),0_i4,type_string)
#if KIM_F03
        kim_ntype(icfg) = kim_ntype(icfg) + 1
        kim_ntype_symbol(kim_ntype(icfg)) = type_string(1:2)
#elif KIM
        kim_ntype(icfg) = kim_ntype(icfg) + 1
        kim_ntype_symbol(kim_ntype(icfg)) = type_string(1:2)
#endif
      endif
    enddo
!
!  Loop over KIM models
!
    model_loop: do nkim = 1,nkimmodel
#ifdef KIM_F03
!
!  Set KIM test descriptor string for this configuration - model combination
!
      call set_kim_test_descriptor(nkim,icfg,kim_ntype(icfg),kim_ntype_symbol)
!
!  Initialize the KIM objects for this configuration-model coupling
!
      KIMerror = kim_api_string_init(pkim_model(nkim,icfg),trim(kim_test_descriptor(nkim,icfg))//char(0),kim_model(nkim))
      if (KIMerror.lt.KIM_STATUS_OK) then
        lkim_model_cfg_OK(nkim,icfg) = .false.
!
!  Having set flag for this model to false then there is no need to continue here
!
        cycle model_loop
!
!        if (ioproc) then
!          KIMerror_d = kim_api_report_error(__LINE__,"setkim","kim_api_init",KIMerror)
!        endif
!        call stopnow('setkim')
      endif
!
!  Set integer variables containing array sizes using KIM integer type
!
      nsize_1 = 1
      nsize_6 = 6
      nsize_N = numat
      nsize_3N = 3*numat
      kim_numat = numat
!
!  Set KIM memory
!
      if (ndim.gt.0) then
        call kim_api_setm_data(pkim_model(nkim,icfg), KIMerror, &
             "numberOfParticles",   nsize_1,  c_loc(pkim_numat),        TRUEFALSE(.true.), &
             "numberOfSpecies",     nsize_1,  c_loc(kim_ntype(icfg)),   TRUEFALSE(.true.), &
             "particleSpecies",     nsize_N,  c_loc(kim_types(1,icfg)), TRUEFALSE(.true.), &
             "coordinates",         nsize_3N, c_loc(kim_coord(1,1)),    TRUEFALSE(.true.), &
             "cutoff",              nsize_1,  c_loc(kim_cutoff(nkim)),  TRUEFALSE(.true.), &
             "energy",              nsize_1,  c_loc(pkim_energy),       TRUEFALSE(.true.), &
             "forces",              nsize_3N, c_loc(kim_forces(1,1)),   TRUEFALSE(.true.), &
             "virial",              nsize_6,  c_loc(kim_virial(1)),     TRUEFALSE(.true.))
      else
        call kim_api_setm_data(pkim_model(nkim,icfg), KIMerror, &
             "numberOfParticles",   nsize_1,  c_loc(pkim_numat),        TRUEFALSE(.true.), &
             "numberOfSpecies",     nsize_1,  c_loc(kim_ntype(icfg)),   TRUEFALSE(.true.), &
             "particleSpecies",     nsize_N,  c_loc(kim_types(1,icfg)), TRUEFALSE(.true.), &
             "coordinates",         nsize_3N, c_loc(kim_coord(1,1)),    TRUEFALSE(.true.), &
             "cutoff",              nsize_1,  c_loc(kim_cutoff(nkim)),  TRUEFALSE(.true.), &
             "energy",              nsize_1,  c_loc(pkim_energy),       TRUEFALSE(.true.), &
             "forces",              nsize_3N, c_loc(kim_forces(1,1)),   TRUEFALSE(.true.))
      endif
      if (KIMerror.lt.KIM_STATUS_OK) then
        if (ioproc) then
          KIMerror_d = kim_api_report_error(__LINE__,"setkim","kim_api_setm_data",KIMerror)
        endif
        call stopnow('setkim')
      endif
      if (KIMerror.lt.KIM_STATUS_OK) then
        if (ioproc) then
          KIMerror_d = kim_api_report_error(__LINE__,"setkim","kim_api_setm_data",KIMerror)
        endif
        call stopnow('setkim')
      endif
!
!  Initialise the KIM model
!
      KIMerror = kim_api_model_init(pkim_model(nkim,icfg))
      if (KIMerror.lt.KIM_STATUS_OK) then
        if (ioproc) then
          KIMerror_d = kim_api_report_error(__LINE__,"setkim","kim_api_model_init",KIMerror)
        endif
        call stopnow('setkim')
      endif
!
!  The following only needs doing once as it is configuration independent
!
!  Find out whether cluster or rvec_f algorithms are supported - do not test the error flag for these calls!
!
      int_logic = kim_api_get_index(pkim_model(nkim,icfg),"NEIGH_RVEC_F",KIMerror)
      if (int_logic.ge.0) then
        kim_nbc(nkim,icfg) = 1
      else
        int_logic = kim_api_get_index(pkim_model(nkim,icfg),"CLUSTER",KIMerror)
        if (int_logic.ge.0) then
          kim_nbc(nkim,icfg) = 0
        else
          call outerror('KIM model does not support CLUSTER or NEIGH_RVEC_F',0_i4)
          call stopnow('setkim')
        endif
      endif
!
!  Set the function names for get_neigh and associate integer for neighObject if NBC => RVEC_F
!
      if (kim_nbc(nkim,icfg).eq.1) then
        KIMerror = kim_api_set_method(pkim_model(nkim,icfg),"get_neigh",nsize_1,c_funloc(get_kim_neigh))
        if (KIMerror.lt.KIM_STATUS_OK) then
          if (ioproc) then
            KIMerror_d = kim_api_report_error(__LINE__,"setkim","kim_api_set_method", KIMerror)
          endif
          call stopnow('setkim')
        endif
        KIMerror = kim_api_set_data(pkim_model(nkim,icfg),"neighObject",nsize_1,c_loc(kim_neighobject(1)))
        if (KIMerror.lt.KIM_STATUS_OK) then
          if (ioproc) then
            KIMerror_d = kim_api_report_error(__LINE__,"setkim","kim_api_set_data", KIMerror)
          endif
          call stopnow('setkim')
        endif
      endif
!
!  Determine KIM type numbers of species for this configuration
!
      do i = 1,numat
!
!  Create element symbol but without any type number
!
        call label(nat(i),0_i4,type_string)
!
!  Check whether symbol is in configuration
!
        kim_types(i,icfg) = kim_api_get_species_code(pkim_model(nkim,icfg),type_string(1:2),KIMerror)
        if (kim_types(i,icfg).eq.KIM_STATUS_PARTICLE_INVALID_SPECIES) then
          call outerror('species type not supported by OpenKIM model',0_i4)
          call stopnow('setkim')
        endif
        if (KIMerror.lt.KIM_STATUS_OK) then
          if (ioproc) then
            KIMerror_d = kim_api_report_error(__LINE__,"setkim","kim_api_get_species_code",KIMerror)
          endif
!
!  Turn off this configuration-model coupling as model does not support the atom type(s)
!
          lkim_model_cfg_OK(nkim,icfg) = .false.
        endif
      enddo
#elif KIM
!
!  Set KIM test descriptor string for this configuration - model combination
!
      call set_kim_test_descriptor(nkim,icfg,kim_ntype(icfg),kim_ntype_symbol)
!
!  Initialize the KIM objects for this configuration-model coupling
!
      KIMerror = kim_api_string_init_f(pkim_model(nkim,icfg),trim(kim_test_descriptor(nkim,icfg))//char(0),kim_model(nkim))
      if (KIMerror.lt.KIM_STATUS_OK) then
        lkim_model_cfg_OK(nkim,icfg) = .false.
!
!  Having set flag for this model to false then there is no need to continue here
!
        cycle model_loop
!
!        if (ioproc) then
!          KIMerror_d = kim_api_report_error_f(__LINE__,"setkim","kim_api_init_f",KIMerror)
!        endif
!        call stopnow('setkim')
      endif
!
!  Set integer variables containing array sizes using KIM integer type
!
      nsize_1 = 1
      nsize_6 = 6
      nsize_N = numat
      nsize_3N = 3*numat
!
      loc_numat = loc(numat)
      loc_ntype = loc(kim_ntype(icfg))
      loc_types = loc(kim_types(1,icfg))
      loc_coords = loc(kim_coord)
      loc_cutoff = loc(kim_cutoff(nkim))
      loc_energy = loc(kim_energy)
      loc_forces = loc(kim_forces)
      loc_virial = loc(kim_virial)
!
!  Set KIM memory 
!
      if (ndim.gt.0) then
        call kim_api_setm_data_f(pkim_model(nkim,icfg), KIMerror, &
             "numberOfParticles",   nsize_1,  loc_numat,   TRUEFALSE(.true.), &
             "numberParticleTypes", nsize_1,  loc_ntype,   TRUEFALSE(.true.), &
             "particleTypes",       nsize_N,  loc_types,   TRUEFALSE(.true.), &
             "coordinates",         nsize_3N, loc_coords,  TRUEFALSE(.true.), &
             "cutoff",              nsize_1,  loc_cutoff,  TRUEFALSE(.true.), &
             "energy",              nsize_1,  loc_energy,  TRUEFALSE(.true.), &
             "forces",              nsize_3N, loc_forces,  TRUEFALSE(.true.), &
             "virial",              nsize_6,  loc_virial,  TRUEFALSE(.true.))
      else
        call kim_api_setm_data_f(pkim_model(nkim,icfg), KIMerror, &
             "numberOfParticles",   nsize_1,  loc_numat,   TRUEFALSE(.true.), &
             "numberParticleTypes", nsize_1,  loc_ntype,   TRUEFALSE(.true.), &
             "particleTypes",       nsize_N,  loc_types,   TRUEFALSE(.true.), &
             "coordinates",         nsize_3N, loc_coords,  TRUEFALSE(.true.), &
             "cutoff",              nsize_1,  loc_cutoff,  TRUEFALSE(.true.), &
             "energy",              nsize_1,  loc_energy,  TRUEFALSE(.true.), &
             "forces",              nsize_3N, loc_forces,  TRUEFALSE(.true.))
      endif
      if (KIMerror.lt.KIM_STATUS_OK) then
        if (ioproc) then
          KIMerror_d = kim_api_report_error_f(__LINE__,"setkim","kim_api_setm_data_f",KIMerror)
        endif
        call stopnow('setkim')
      endif
!
!  Initialise the KIM model
!
      KIMerror = kim_api_model_init_f(pkim_model(nkim,icfg))
      if (KIMerror.lt.KIM_STATUS_OK) then
        if (ioproc) then
          KIMerror_d = kim_api_report_error_f(__LINE__,"setkim","kim_api_model_init",KIMerror)
        endif
        call stopnow('setkim')
      endif
!
!  The following only needs doing once as it is configuration independent
!
!  Find out whether cluster or rvec_f algorithms are supported - do not test the error flag for these calls!
!
      int_logic = kim_api_get_index_f(pkim_model(nkim,icfg),"NEIGH_RVEC_F",KIMerror)
      if (int_logic.ge.0) then
        kim_nbc(nkim,icfg) = 1
      else
        int_logic = kim_api_get_index_f(pkim_model(nkim,icfg),"CLUSTER",KIMerror)
        if (int_logic.ge.0) then
          kim_nbc(nkim,icfg) = 0
        else
          call outerror('KIM model does not support CLUSTER or NEIGH_RVEC_F',0_i4)
          call stopnow('setkim')
        endif
      endif
!
!  Set the function names for get_neigh and associate integer for neighObject if NBC => RVEC_F
!
      if (kim_nbc(nkim,icfg).eq.1) then
        KIMerror = kim_api_set_data_f(pkim_model(nkim,icfg),"get_neigh",nsize_1,loc(get_kim_neigh))
        if (KIMerror.lt.KIM_STATUS_OK) then
          if (ioproc) then
            KIMerror_d = kim_api_report_error_f(__LINE__,"setkim","kim_api_set_data_f", KIMerror)
          endif
          call stopnow('setkim')
        endif
        KIMerror = kim_api_set_data_f(pkim_model(nkim,icfg),"neighObject",nsize_1,loc(kim_neighobject))
        if (KIMerror.lt.KIM_STATUS_OK) then
          if (ioproc) then
            KIMerror_d = kim_api_report_error_f(__LINE__,"setkim","kim_api_set_data_f", KIMerror)
          endif
          call stopnow('setkim')
        endif
      endif
!
!  Determine KIM type numbers of species for this configuration
!
      do i = 1,numat
!
!  Create element symbol but without any type number
!
        call label(nat(i),0_i4,type_string)
!
!  Check whether symbol is in configuration
!
        kim_types(i,icfg) = kim_api_get_partcl_type_code_f(pkim_model(nkim,icfg),type_string(1:2),KIMerror)
        if (KIMerror.lt.KIM_STATUS_OK) then
          if (ioproc) then
            KIMerror_d = kim_api_report_error_f(__LINE__,"setkim","kim_api_get_partcl_type_code_f",KIMerror)
          endif
!
!  Turn off this configuration-model coupling as model does not support the atom type(s)
!
          lkim_model_cfg_OK(nkim,icfg) = .false.
        endif
      enddo
#else
      lkim_model_cfg_OK(nkim,icfg) = .false.
#endif
!
!  End of loop over KIM models
!
    enddo model_loop
!
!  End of loop over configurations
!
  enddo
!
!  Deallocate workspace array for type symbols
!
  deallocate(kim_ntype_symbol)
!
  return
  end
!
#if KIM_F03
  subroutine set_kim_test_descriptor(nkim,icfg,kim_ntype,kim_ntype_symbol)

  use, intrinsic :: iso_c_binding
  use datatypes
  use KIM_API_F03
  use current,      only : ndim
  use kim_models,   only : kim_test_descriptor, kim_model
  use parallel,     only : ioproc
  implicit none
!
!  Passed variables
!
  integer(i4),                          intent(in)   :: nkim             ! Current KIM model number
  integer(i4),                          intent(in)   :: icfg             ! Current configuration number
  integer(i4),                          intent(in)   :: kim_ntype
  character(len=2),                     intent(in)   :: kim_ntype_symbol(kim_ntype)
!
!  Local variables
!
  type(c_ptr)                                        :: pkim_tmp         ! Local C pointer used for model string
  type(c_ptr)                                        :: pdata            ! Local C pointer used for data testing
  integer(i4)                                        :: i                ! Looping index over KIM types
  integer(c_int)                                     :: KIMerror         ! KIM error flag
  integer(c_int)                                     :: KIMerror_d       ! KIM error flag - dummy
  character(len=103)                                 :: divider          ! Line that divides sections
  character(len=1)                                   :: cr               ! Carriage return
  character(len=56)                                  :: type_line        ! String used in symbol writing
  character(len=24)                                  :: type24           ! String used in symbol writing
  character(len=28)                                  :: type28           ! String used in write NBC line
  character(len=32)                                  :: nbcline          ! String used to write NBC line
  logical                                            :: lcluster
  logical                                            :: ldestroy
  logical                                            :: lBothAccess
  logical                                            :: lIterAccess
  logical                                            :: lLocaAccess
  logical                                            :: lrvecf
!
!  Define frequently used variables
!
  cr = char(10)
  divider = '#######################################################################################################'
!
!  Get model information so that we can find out what we need to write
!
  KIMerror = kim_api_model_info(pkim_tmp,kim_model(nkim))
  if (KIMerror.lt.KIM_STATUS_OK) then
    if (ioproc) then
      KIMerror_d = kim_api_report_error(__LINE__,"set_kim_test_descriptor","kim_api_init",KIMerror)
    endif
    call stopnow('set_kim_test_descriptor')
  endif
!
!  Test whether various options are present or not
!
  pdata = kim_api_get_data(pkim_tmp,"destroy",KIMerror)
  ldestroy = (KIMerror.eq.KIM_STATUS_OK)
!
  pdata = kim_api_get_data(pkim_tmp,"CLUSTER",KIMerror)
  lcluster = (KIMerror.eq.KIM_STATUS_OK)
!
  pdata = kim_api_get_data(pkim_tmp,"NEIGH_RVEC_F",KIMerror)
  lrvecf = (KIMerror.eq.KIM_STATUS_OK)
!
  pdata = kim_api_get_data(pkim_tmp,"Neigh_IterAccess",KIMerror)
  lIterAccess = (KIMerror.eq.KIM_STATUS_OK)
!
  pdata = kim_api_get_data(pkim_tmp,"Neigh_LocaAccess",KIMerror)
  lLocaAccess = (KIMerror.eq.KIM_STATUS_OK)
!
  pdata = kim_api_get_data(pkim_tmp,"Neigh_BothAccess",KIMerror)
  lBothAccess = (KIMerror.eq.KIM_STATUS_OK)
!
!  Initialise string
!
  kim_test_descriptor(nkim,icfg) = ' '
!
!  Write GULP descriptor into string kim_test_descriptor
!
  kim_test_descriptor(nkim,icfg) = &
    divider                                                                     // cr // &
    'TEST_NAME := gulp_kim_test'                                                // cr // &
    'KIM_API_Version  := 1.7.3'                                                 // cr // &
    divider                                                                     // cr // &
    'Unit_length      := A'                                                     // cr // &
    'Unit_energy      := eV'                                                    // cr // &
    'Unit_charge      := e'                                                     // cr // &
    'Unit_temperature := K'                                                     // cr // &
    'Unit_time        := ps'                                                    // cr // &
                                                                                   cr // &
                                                                                   cr // &
    divider                                                                     // cr // &
    'PARTICLE_SPECIES:'                                                         // cr // &
    '# Symbol/name               Type                    code'                  // cr // cr

  do i = 1,kim_ntype
    type24 = kim_ntype_symbol(i)
    write(type_line,'(a24,4x,''spec'',17x,i4)') type24,0
    kim_test_descriptor(nkim,icfg) = trim(kim_test_descriptor(nkim,icfg)) // type_line // cr
  enddo
!
!  Set NBC string
!
  if (ndim.eq.3) then
    if (lrvecf) then
      type28 = 'NEIGH_RVEC_F'
    elseif (lcluster) then
      type28 = 'CLUSTER'
    else
      type28 = ' '
    endif
  else
    if (lcluster) then
      type28 = 'CLUSTER'
    elseif (lrvecf) then
      type28 = 'NEIGH_RVEC_F'
    else
      type28 = ' '
    endif
  endif
  nbcline = type28 // 'flag'
!
!  Write out conventions block
!
  kim_test_descriptor(nkim,icfg) = trim(kim_test_descriptor(nkim,icfg)) // &
                                                                                   cr // &
                                                                                   cr // &
    divider                                                                     // cr // &
    'CONVENTIONS:'                                                              // cr // &
    '# Name                      Type'                                          // cr // &
                                                                                   cr // &
    'OneBasedLists               flag'                                          // cr // &
                                                                                   cr // &
    nbcline                                                                     // cr // &
                                                                                   cr

  if (lBothAccess) then
    kim_test_descriptor(nkim,icfg) = trim(kim_test_descriptor(nkim,icfg)) // &
      'Neigh_BothAccess            flag'                                          // cr // cr
  elseif (lLocaAccess) then
    kim_test_descriptor(nkim,icfg) = trim(kim_test_descriptor(nkim,icfg)) // &
      'Neigh_LocaAccess            flag'                                          // cr // cr
  elseif (lIterAccess) then
    kim_test_descriptor(nkim,icfg) = trim(kim_test_descriptor(nkim,icfg)) // &
      'Neigh_IterAccess            flag'                                          // cr // cr
  endif
!
!  Divider to close conventions block
!
  kim_test_descriptor(nkim,icfg) = trim(kim_test_descriptor(nkim,icfg)) // divider // cr
!
!  Write out model input block
!
  kim_test_descriptor(nkim,icfg) = trim(kim_test_descriptor(nkim,icfg)) // &
    'MODEL_INPUT:'                                                              // cr // &
    '# Name                      Type         Unit       Shape              requirements' // cr // cr // &
    'numberOfParticles           integer      none       []'                    // cr // &
                                                                                   cr // &
    'numberOfSpecies             integer      none       []'                    // cr // &
                                                                                   cr // &
    'particleSpecies             integer      none       [numberOfParticles]'   // cr // &
                                                                                   cr // &
    'coordinates                 double       length     [numberOfParticles,3]' // cr // &
                                                                                   cr
  if (ndim.gt.0) then
    kim_test_descriptor(nkim,icfg) = trim(kim_test_descriptor(nkim,icfg)) // &
      'boxSideLengths              double       length     [3]'                   // cr // cr
  endif
  if (lrvecf) then
    kim_test_descriptor(nkim,icfg) = trim(kim_test_descriptor(nkim,icfg)) // &
      'get_neigh                   method       none       []'                    // cr // cr
    kim_test_descriptor(nkim,icfg) = trim(kim_test_descriptor(nkim,icfg)) // &
      'neighObject                 pointer      none       []'                    // cr // cr
  endif
!
!  Divider to close model input block
!
  kim_test_descriptor(nkim,icfg) = trim(kim_test_descriptor(nkim,icfg)) // divider // cr
!
!  Write out model output block
!
  kim_test_descriptor(nkim,icfg) = trim(kim_test_descriptor(nkim,icfg)) // &
    'MODEL_OUTPUT:'                                                             // cr // &
    '# Name                      Type         Unit       Shape              requirements' // cr // &
                                                                                   cr // &
    'compute                     method       none       []'                    // cr // &
                                                                                   cr // &
    'cutoff                      double       length     []'                    // cr // &
                                                                                   cr // &
    'energy                      double       energy     []'                    // cr // &
                                                                                   cr // &
    'forces                      double       force      [numberOfParticles,3]' // cr // &
                                                                                   cr
!
!  Optional model output strings
!
  if (ndim.gt.0) then
    kim_test_descriptor(nkim,icfg) = trim(kim_test_descriptor(nkim,icfg)) // &
      'virial                      double       energy     [6]'                   // cr // cr
  endif
  if (ldestroy) then
    kim_test_descriptor(nkim,icfg) = trim(kim_test_descriptor(nkim,icfg)) // &
      'destroy                     method       none       []'                    // cr // cr
  endif

  return
  end 
#elif KIM
  subroutine set_kim_test_descriptor(nkim,icfg,kim_ntype,kim_ntype_symbol)

  use datatypes
  use KIM_API
  use current,      only : ndim
  use kim_models,   only : kim_test_descriptor, kim_model
  use parallel,     only : ioproc
  implicit none
!
!  Passed variables
!
  integer(i4),                          intent(in)   :: nkim             ! Current KIM model number
  integer(i4),                          intent(in)   :: icfg             ! Current configuration number
  integer(i4),                          intent(in)   :: kim_ntype
  character(len=2),                     intent(in)   :: kim_ntype_symbol(kim_ntype)
!
!  Local variables
!
  integer(kind=kim_intptr)                           :: pkim_tmp         ! Local pkim for model string
  integer(kind=kim_intptr)                           :: pdata            ! Local pointer used for data testing
  integer(i4)                                        :: i                ! Looping index over KIM types
  integer(i4)                                        :: KIMerror         ! KIM error flag
  integer(i4)                                        :: KIMerror_d       ! KIM error flag - dummy
  character(len=103)                                 :: divider          ! Line that divides sections
  character(len=1)                                   :: cr               ! Carriage return
  character(len=56)                                  :: type_line        ! String used in symbol writing
  character(len=24)                                  :: type24           ! String used in symbol writing
  character(len=28)                                  :: type28           ! String used in write NBC line
  character(len=32)                                  :: nbcline          ! String used to write NBC line
  logical                                            :: lcluster
  logical                                            :: ldestroy
  logical                                            :: lIterAccess
  logical                                            :: lLocaAccess
  logical                                            :: lrvecf
!
!  Define frequently used variables
!
  cr = char(10)
  divider = '#######################################################################################################'
!
!  Get model information so that we can find out what we need to write
!
  KIMerror = kim_api_model_info_f(pkim_tmp,kim_model(nkim))
  if (KIMerror.lt.KIM_STATUS_OK) then
    if (ioproc) then
      KIMerror_d = kim_api_report_error_f(__LINE__,"set_kim_test_descriptor","kim_api_init_f",KIMerror)
    endif
    call stopnow('set_kim_test_descriptor')
  endif
!
!  Test whether various options are present or not
!
  pdata = kim_api_get_data_f(pkim_tmp,"destroy",KIMerror)
  ldestroy = (KIMerror.eq.KIM_STATUS_OK)
!
  pdata = kim_api_get_data_f(pkim_tmp,"CLUSTER",KIMerror)
  lcluster = (KIMerror.eq.KIM_STATUS_OK)
!
  pdata = kim_api_get_data_f(pkim_tmp,"NEIGH_RVEC_F",KIMerror)
  lrvecf = (KIMerror.eq.KIM_STATUS_OK)
!
  pdata = kim_api_get_data_f(pkim_tmp,"Neigh_IterAccess",KIMerror)
  lIterAccess = (KIMerror.eq.KIM_STATUS_OK)
!
  pdata = kim_api_get_data_f(pkim_tmp,"Neigh_LocaAccess",KIMerror)
  lLocaAccess = (KIMerror.eq.KIM_STATUS_OK)
!
!  Initialise string
!
  kim_test_descriptor(nkim,icfg) = ' '
!
!  Write GULP descriptor into string kim_test_descriptor
!
  kim_test_descriptor(nkim,icfg) = &
    divider                                                                     // cr // &
    'TEST_NAME := gulp_kim_test'                                                // cr // &
    'Unit_length      := A'                                                     // cr // &
    'Unit_energy      := eV'                                                    // cr // &
    'Unit_charge      := e'                                                     // cr // &
    'Unit_temperature := K'                                                     // cr // &
    'Unit_time        := ps'                                                    // cr // &
                                                                                   cr // &
                                                                                   cr // &
    divider                                                                     // cr // &
    'SUPPORTED_ATOM/PARTICLES_TYPES:'                                           // cr // &
    '# Symbol/name               Type                    code'                  // cr // cr

  do i = 1,kim_ntype
    type24 = kim_ntype_symbol(i)
    write(type_line,'(a24,4x,''spec'',17x,i4)') type24,0
    kim_test_descriptor(nkim,icfg) = trim(kim_test_descriptor(nkim,icfg)) // type_line // cr
  enddo
!
!  Set NBC string
!
  if (ndim.eq.3) then
    if (lrvecf) then
      type28 = 'NEIGH_RVEC_F'
    elseif (lcluster) then
      type28 = 'CLUSTER'
    else
      type28 = ' '
    endif
  else
    if (lcluster) then
      type28 = 'CLUSTER' 
    elseif (lrvecf) then
      type28 = 'NEIGH_RVEC_F'
    else
      type28 = ' '
    endif                                                                          
  endif
  nbcline = type28 // 'flag'
!
!  Write out conventions block
!
  kim_test_descriptor(nkim,icfg) = trim(kim_test_descriptor(nkim,icfg)) // &
                                                                                   cr // &
                                                                                   cr // &
    divider                                                                     // cr // &
    'CONVENTIONS:'                                                              // cr // &
    '# Name                      Type'                                          // cr // &
                                                                                   cr // &
    'OneBasedLists               flag'                                          // cr // &
                                                                                   cr // &
    nbcline                                                                     // cr // &
                                                                                   cr

  if (lLocaAccess) then
    kim_test_descriptor(nkim,icfg) = trim(kim_test_descriptor(nkim,icfg)) // &
      'Neigh_LocaAccess            flag'                                          // cr // cr
  endif
  if (lIterAccess) then
    kim_test_descriptor(nkim,icfg) = trim(kim_test_descriptor(nkim,icfg)) // &
      'Neigh_IterAccess            flag'                                          // cr // cr
  endif
!
!  Divider to close conventions block
!
  kim_test_descriptor(nkim,icfg) = trim(kim_test_descriptor(nkim,icfg)) // divider // cr
!
!  Write out model input block
!
  kim_test_descriptor(nkim,icfg) = trim(kim_test_descriptor(nkim,icfg)) // &
    'MODEL_INPUT:'                                                              // cr // &
    '# Name                      Type         Unit       Shape              requirements' // cr // cr // &
    'numberOfParticles           integer      none       []'                    // cr // &
                                                                                   cr // &
    'numberParticleTypes         integer      none       []'                    // cr // &
                                                                                   cr // &
    'particleTypes               integer      none       [numberOfParticles]'   // cr // &
                                                                                   cr // &
    'coordinates                 real*8       length     [numberOfParticles,3]' // cr // &
                                                                                   cr
  if (ndim.gt.0) then
    kim_test_descriptor(nkim,icfg) = trim(kim_test_descriptor(nkim,icfg)) // &
      'boxSideLengths              real*8       length     [3]'                   // cr // cr
  endif
  if (lrvecf) then
    kim_test_descriptor(nkim,icfg) = trim(kim_test_descriptor(nkim,icfg)) // &
      'get_neigh                   method       none       []'                    // cr // cr
    kim_test_descriptor(nkim,icfg) = trim(kim_test_descriptor(nkim,icfg)) // &
      'neighObject                 pointer      none       []'                    // cr // cr
  endif
!
!  Divider to close model input block
!
  kim_test_descriptor(nkim,icfg) = trim(kim_test_descriptor(nkim,icfg)) // divider // cr
!
!  Write out model output block
!
  kim_test_descriptor(nkim,icfg) = trim(kim_test_descriptor(nkim,icfg)) // &
    'MODEL_OUTPUT:'                                                             // cr // &
    '# Name                      Type         Unit       Shape              requirements' // cr // &
                                                                                   cr // &
    'compute                     method       none       []'                    // cr // &
                                                                                   cr // &
    'cutoff                      real*8       length     []'                    // cr // &
                                                                                   cr // &
    'energy                      real*8       energy     []'                    // cr // &
                                                                                   cr // &
    'forces                      real*8       force      [numberOfParticles,3]' // cr // &
                                                                                   cr
!
!  Optional model output strings
!
  if (ndim.gt.0) then
    kim_test_descriptor(nkim,icfg) = trim(kim_test_descriptor(nkim,icfg)) // &
      'virial                      real*8       energy     [6]'                   // cr // cr
  endif
  if (ldestroy) then
    kim_test_descriptor(nkim,icfg) = trim(kim_test_descriptor(nkim,icfg)) // &
      'destroy                     method       none       []'                    // cr // cr
  endif

  return
  end 
#endif
