  subroutine changemaxfreqat(nfqat)
!
!  Alters the size of the arrays associated with maxfqat
!
!   4/11 Created from changemaxat
!   8/14 eigv and groupvelocity added
!   8/16 IR intensity array added
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
!  Julian Gale, CIC, Curtin University, August 2016
!
  use datatypes
  use control,        only : lgroupvelocity
  use frequencies,    only : freq, eigv, maxfkpt, maxfqat, lStoreEig, groupvelocity, IRintensity
  use reallocate
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in) :: nfqat
!
!  Local variables
!
  integer(i4)             :: ierror
!
  if (nfqat.gt.maxfqat) then
    maxfqat = nfqat
    call realloc(freq,3_i4*maxfqat,maxfkpt,ierror)
    if (ierror.ne.0) call outofmemory('changemaxfreqat','freq')
    call realloc(IRintensity,3_i4*maxfqat,maxfkpt,ierror)
    if (ierror.ne.0) call outofmemory('changemaxfreqat','IRintensity')
    if (lStoreEig) then
      call realloc(eigv,3_i4*maxfqat,3_i4*maxfqat,maxfkpt,ierror)
      if (ierror.ne.0) call outofmemory('changemaxfreqat','eigv')
    endif
    if (lgroupvelocity) then
      call realloc(groupvelocity,3_i4,3_i4*maxfqat,maxfkpt,ierror)
      if (ierror.ne.0) call outofmemory('changemaxfreqat','groupvelocity')
    endif
  endif
!
  return
  end
