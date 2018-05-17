  module reallocate_kim
!
!  Set of routines and interface for handling memory
!  allocation and reallocation while preserving the
!  contents. This version is specifically for the 
!  kind kim_intptr used by OpenKIM.
!
!  10/12 Created from reallocate
!   4/14 Modified for KIM version 1.4
!   7/14 pkim changed from c_ptr to c_int
!   7/16 Checked for KIM version 1.7.3
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
    use datatypes
    use iochannels
#if KIM_F03
    use kim_models, only : kim_len
    use KIM_API_F03
!
!  Local variables for storing memory usage information
!
    interface realloc_kim

      module procedure &
        realloc_kim_1, &
        realloc_kim_2, &
        realloc_kim_3, &
        realloc_kim_4

    end interface 

    interface realloc_kim_len

      module procedure &
        realloc_kim_len_1, realloc_kim_len_2, &
        realloc_kim_len_3, realloc_kim_len_4

    end interface 

  contains

    subroutine realloc_kim_1(p,n,ierror)

      implicit none
!
!  Passed arguments
!
      integer(c_int),   dimension(:), pointer :: p
      integer(i4), intent(in)                 :: n
      integer(i4), intent(out)                :: ierror
!
!  Local arguments
!
      integer(c_int),   dimension(:), pointer :: plocal
      integer(i4)                             :: nold
      integer(i4)                             :: i
!
!  Allocate new size of memory
!
      nullify(plocal)
      allocate(plocal(n),stat=ierror)
!
!  If there is an error then return
!
      if (ierror.ne.0) then
        write(ioout,'(''!! Allocation error !!'')')
        return
      endif
      plocal = 0_c_int
!
      if (associated(p)) then
!
!  Find size of old array
!
        nold = size(p)
!
!  Preserve data
!
        do i = 1,min(nold,n)
          plocal(i) = p(i)
        enddo
!
!  Free old memory
!
        deallocate(p)
      endif
!
!  Re-assign pointers so that p is returned as the new array
!
      p=>plocal

    end subroutine realloc_kim_1

    subroutine realloc_kim_2(p,m,n,ierror)

      implicit none
!
!  Passed arguments
!
      integer(c_int),           dimension(:,:), pointer :: p
      integer(i4), intent(in)                           :: m,n
      integer(i4), intent(out)                          :: ierror
!
!  Local arguments
!
      integer(c_int),           dimension(:,:), pointer :: plocal
      integer(i4)                                       :: mold,nold
      integer(i4)                                       :: i, j
!
!  Allocate new size of memory
!
      nullify(plocal)
      allocate(plocal(m,n),stat=ierror)
!
!  If there is an error then return
!
      if (ierror.ne.0) then
        write(ioout,'(''!! Allocation error !!'')')
        return
      endif
      plocal = 0_c_int
!
      if (associated(p)) then
!
!  Find size of old array
!
        mold = size(p,1)
        nold = size(p,2)
!
!  Preserve data
!
        do i = 1,min(nold,n)
          do j = 1,min(mold,m)
            plocal(j,i) = p(j,i)
          enddo
        enddo
!
!  Free old memory
!
        deallocate(p)
      endif
!
!  Re-assign pointers so that p is returned as the new array
!
      p=>plocal

    end subroutine realloc_kim_2

    subroutine realloc_kim_3(p,k,m,n,ierror)

      implicit none
!
!  Passed arguments
!
      integer(c_int),           dimension(:,:,:), pointer :: p
      integer(i4), intent(in)                             :: k,m,n
      integer(i4), intent(out)                            :: ierror
!
!  Local arguments
!
      integer(c_int),           dimension(:,:,:), pointer :: plocal
      integer(i4)                                         :: kold,mold,nold
      integer(i4)                                         :: i, j, l
!
!  Allocate new size of memory
!
      nullify(plocal)
      allocate(plocal(k,m,n),stat=ierror)
!
!  If there is an error then return
!
      if (ierror.ne.0) then
        write(ioout,'(''!! Allocation error !!'')')
        return
      endif
      plocal = 0_c_int
!
      if (associated(p)) then
!
!  Find size of old array
!
        kold = size(p,1)
        mold = size(p,2)
        nold = size(p,3)
!
!  Preserve data
!
        do i = 1,min(nold,n)
          do j = 1,min(mold,m)
            do l = 1,min(kold,k)
              plocal(l,j,i) = p(l,j,i)
            enddo
          enddo
        enddo
!
!  Free old memory
!
        deallocate(p)
      endif
!
!  Re-assign pointers so that p is returned as the new array
!
      p=>plocal

    end subroutine realloc_kim_3

    subroutine realloc_kim_4(p,k,l,m,n,ierror)

      implicit none
!
!  Passed arguments
!
      integer(c_int),           dimension(:,:,:,:), pointer :: p
      integer(i4), intent(in)                               :: k,l,m,n
      integer(i4), intent(out)                              :: ierror
!
!  Local arguments
!
      integer(c_int),           dimension(:,:,:,:), pointer :: plocal
      integer(i4)                                           :: kold,lold
      integer(i4)                                           :: mold,nold
!
!  Allocate new size of memory
!
      nullify(plocal)
      allocate(plocal(k,l,m,n),stat=ierror)
!
!  If there is an error then return
!
      if (ierror.ne.0) then
        write(ioout,'(''!! Allocation error !!'')')
        return
      endif
      plocal = 0_c_int
!
      if (associated(p)) then
!
!  Find size of old array
!
        kold = size(p,1)
        lold = size(p,2)
        mold = size(p,3)
        nold = size(p,4)
!
!  Preserve data
!
        plocal(1:min(kold,k),1:min(lold,l),1:min(mold,m), &
          1:min(nold,n)) = p(1:min(kold,k),1:min(lold,l), &
          1:min(mold,m),1:min(nold,n))
!
!  Free old memory
!
        deallocate(p)
      endif
!
!  Re-assign pointers so that p is returned as the new array
!
      p=>plocal

    end subroutine realloc_kim_4

    subroutine realloc_kim_len_1(p,n,ierror)

      implicit none
!
!  Passed arguments
!
      integer(i4), intent(in)                       :: n
      integer(i4), intent(out)                      :: ierror
      character(len=kim_len), dimension(:), pointer :: p
!
!  Local arguments
!
      character(len=kim_len), dimension(:), pointer :: plocal
      integer(i4)                                   :: nold
!
!  Allocate new size of memory
!
      nullify(plocal)
      allocate(plocal(n),stat=ierror)
!
!  If there is an error then return
!
      if (ierror.ne.0) then
        write(ioout,'(''!! Allocation error !!'')')
        return
      endif
      plocal = ' '
!
      if (associated(p)) then
!
!  Find size of old array
!
        nold = size(p)
!
!  Preserve data
!
        plocal(1:min(nold,n)) = p(1:min(nold,n))
!
!  Free old memory
!
        deallocate(p)
      endif
!
!  Re-assign pointers so that p is returned as the new array
!
      p=>plocal

    end subroutine realloc_kim_len_1

    subroutine realloc_kim_len_2(p,m,n,ierror)

      implicit none
!
!  Passed arguments
!
      integer(i4), intent(in)                         :: m,n
      integer(i4), intent(out)                        :: ierror
      character(len=kim_len), dimension(:,:), pointer :: p
!
!  Local arguments
!
      character(len=kim_len), dimension(:,:), pointer :: plocal
      integer(i4)                                     :: mold,nold
!
!  Allocate new size of memory
!
      nullify(plocal)
      allocate(plocal(m,n),stat=ierror)
!
!  If there is an error then return
!
      if (ierror.ne.0) then
        write(ioout,'(''!! Allocation error !!'')')
        return
      endif
      plocal = ' '
!
      if (associated(p)) then
!
!  Find size of old array
!
        mold = size(p,1)
        nold = size(p,2)
!
!  Preserve data
!
        plocal(1:min(mold,m),1:min(nold,n)) = &
          p(1:min(mold,m),1:min(nold,n))
!
!  Free old memory
!
        deallocate(p)
      endif
!
!  Re-assign pointers so that p is returned as the new array
!
      p=>plocal

    end subroutine realloc_kim_len_2

    subroutine realloc_kim_len_3(p,k,m,n,ierror)

      implicit none
!
!  Passed arguments
!
      integer(i4), intent(in)                           :: k,m,n
      integer(i4), intent(out)                          :: ierror
      character(len=kim_len), dimension(:,:,:), pointer :: p
!
!  Local arguments
!
      character(len=kim_len), dimension(:,:,:), pointer :: plocal
      integer(i4)                                       :: mold,nold
      integer(i4)                                       :: kold
!
!  Allocate new size of memory
!
      nullify(plocal)
      allocate(plocal(k,m,n),stat=ierror)
!
!  If there is an error then return
!
      if (ierror.ne.0) then
        write(ioout,'(''!! Allocation error !!'')')
        return
      endif
      plocal = ' '
!
      if (associated(p)) then
!
!  Find size of old array
!
        kold = size(p,1)
        mold = size(p,2)
        nold = size(p,3)
!
!  Preserve data
!
        plocal(1:min(kold,k),1:min(mold,m),1:min(nold,n)) = &
          p(1:min(kold,k),1:min(mold,m),1:min(nold,n))
!
!  Free old memory
!
        deallocate(p)
      endif
!
!  Re-assign pointers so that p is returned as the new array
!
      p=>plocal

    end subroutine realloc_kim_len_3

    subroutine realloc_kim_len_4(p,k,l,m,n,ierror)

      implicit none
!
!  Passed arguments
!
      integer(i4), intent(in)                             :: k,l,m,n
      integer(i4), intent(out)                            :: ierror
      character(len=kim_len), dimension(:,:,:,:), pointer :: p
!
!  Local arguments
!
      character(len=kim_len), dimension(:,:,:,:), pointer :: plocal
      integer(i4)                                         :: mold,nold
      integer(i4)                                         :: kold,lold
!
!  Allocate new size of memory
!
      nullify(plocal)
      allocate(plocal(k,l,m,n),stat=ierror)
!
!  If there is an error then return
!
      if (ierror.ne.0) then
        write(ioout,'(''!! Allocation error !!'')')
        return
      endif
      plocal = ' '
!
      if (associated(p)) then
!
!  Find size of old array
!
        kold = size(p,1)
        lold = size(p,2)
        mold = size(p,3)
        nold = size(p,4)
!
!  Preserve data
!
        plocal(1:min(kold,k),1:min(lold,l),1:min(mold,m), &
          1:min(nold,n)) = p(1:min(kold,k),1:min(lold,l), &
          1:min(mold,m),1:min(nold,n))
!
!  Free old memory
!
        deallocate(p)
      endif
!
!  Re-assign pointers so that p is returned as the new array
!
      p=>plocal

    end subroutine realloc_kim_len_4
#elif KIM
    use kim_models, only : kim_len
    use KIM_API
!
!  Local variables for storing memory usage information
!
    interface realloc_kim

      module procedure &
        realloc_kim_1, &
        realloc_kim_2, &
        realloc_kim_3, &
        realloc_kim_4

    end interface 

    interface realloc_kim_len

      module procedure &
        realloc_kim_len_1, realloc_kim_len_2, &
        realloc_kim_len_3, realloc_kim_len_4

    end interface 

  contains

    subroutine realloc_kim_1(p,n,ierror)

      implicit none
!
!  Passed arguments
!
      integer(kind=kim_intptr), dimension(:), pointer :: p
      integer(i4), intent(in)                         :: n
      integer(i4), intent(out)                        :: ierror
!
!  Local arguments
!
      integer(kind=kim_intptr), dimension(:), pointer :: plocal
      integer(i4)                                     :: nold
      integer(i4)                                     :: i
!
!  Allocate new size of memory
!
      nullify(plocal)
      allocate(plocal(n),stat=ierror)
!
!  If there is an error then return
!
      if (ierror.ne.0) then
        write(ioout,'(''!! Allocation error !!'')')
        return
      endif
      plocal = 0_kim_intptr
!
      if (associated(p)) then
!
!  Find size of old array
!
        nold = size(p)
!
!  Preserve data
!
        do i = 1,min(nold,n)
          plocal(i) = p(i)
        enddo
!
!  Free old memory
!
        deallocate(p)
      endif
!
!  Re-assign pointers so that p is returned as the new array
!
      p=>plocal

    end subroutine realloc_kim_1

    subroutine realloc_kim_2(p,m,n,ierror)

      implicit none
!
!  Passed arguments
!
      integer(kind=kim_intptr), dimension(:,:), pointer :: p
      integer(i4), intent(in)                           :: m,n
      integer(i4), intent(out)                          :: ierror
!
!  Local arguments
!
      integer(kind=kim_intptr), dimension(:,:), pointer :: plocal
      integer(i4)                                       :: mold,nold
      integer(i4)                                       :: i, j
!
!  Allocate new size of memory
!
      nullify(plocal)
      allocate(plocal(m,n),stat=ierror)
!
!  If there is an error then return
!
      if (ierror.ne.0) then
        write(ioout,'(''!! Allocation error !!'')')
        return
      endif
      plocal = 0_kim_intptr
!
      if (associated(p)) then
!
!  Find size of old array
!
        mold = size(p,1)
        nold = size(p,2)
!
!  Preserve data
!
        do i = 1,min(nold,n)
          do j = 1,min(mold,m)
            plocal(j,i) = p(j,i)
          enddo
        enddo
!
!  Free old memory
!
        deallocate(p)
      endif
!
!  Re-assign pointers so that p is returned as the new array
!
      p=>plocal

    end subroutine realloc_kim_2

    subroutine realloc_kim_3(p,k,m,n,ierror)

      implicit none
!
!  Passed arguments
!
      integer(kind=kim_intptr), dimension(:,:,:), pointer :: p
      integer(i4), intent(in)                             :: k,m,n
      integer(i4), intent(out)                            :: ierror
!
!  Local arguments
!
      integer(kind=kim_intptr), dimension(:,:,:), pointer :: plocal
      integer(i4)                                         :: kold,mold,nold
      integer(i4)                                         :: i, j, l
!
!  Allocate new size of memory
!
      nullify(plocal)
      allocate(plocal(k,m,n),stat=ierror)
!
!  If there is an error then return
!
      if (ierror.ne.0) then
        write(ioout,'(''!! Allocation error !!'')')
        return
      endif
      plocal = 0_kim_intptr
!
      if (associated(p)) then
!
!  Find size of old array
!
        kold = size(p,1)
        mold = size(p,2)
        nold = size(p,3)
!
!  Preserve data
!
        do i = 1,min(nold,n)
          do j = 1,min(mold,m)
            do l = 1,min(kold,k)
              plocal(l,j,i) = p(l,j,i)
            enddo
          enddo
        enddo
!
!  Free old memory
!
        deallocate(p)
      endif
!
!  Re-assign pointers so that p is returned as the new array
!
      p=>plocal

    end subroutine realloc_kim_3

    subroutine realloc_kim_4(p,k,l,m,n,ierror)

      implicit none
!
!  Passed arguments
!
      integer(kind=kim_intptr), dimension(:,:,:,:), pointer :: p
      integer(i4), intent(in)                               :: k,l,m,n
      integer(i4), intent(out)                              :: ierror
!
!  Local arguments
!
      integer(kind=kim_intptr), dimension(:,:,:,:), pointer :: plocal
      integer(i4)                                           :: kold,lold
      integer(i4)                                           :: mold,nold
!
!  Allocate new size of memory
!
      nullify(plocal)
      allocate(plocal(k,l,m,n),stat=ierror)
!
!  If there is an error then return
!
      if (ierror.ne.0) then
        write(ioout,'(''!! Allocation error !!'')')
        return
      endif
      plocal = 0_kim_intptr
!
      if (associated(p)) then
!
!  Find size of old array
!
        kold = size(p,1)
        lold = size(p,2)
        mold = size(p,3)
        nold = size(p,4)
!
!  Preserve data
!
        plocal(1:min(kold,k),1:min(lold,l),1:min(mold,m), &
          1:min(nold,n)) = p(1:min(kold,k),1:min(lold,l), &
          1:min(mold,m),1:min(nold,n))
!
!  Free old memory
!
        deallocate(p)
      endif
!
!  Re-assign pointers so that p is returned as the new array
!
      p=>plocal

    end subroutine realloc_kim_4

    subroutine realloc_kim_len_1(p,n,ierror)

      implicit none
!
!  Passed arguments
!
      integer(i4), intent(in)                       :: n
      integer(i4), intent(out)                      :: ierror
      character(len=kim_len), dimension(:), pointer :: p
!
!  Local arguments
!
      character(len=kim_len), dimension(:), pointer :: plocal
      integer(i4)                                   :: nold
!
!  Allocate new size of memory
!
      nullify(plocal)
      allocate(plocal(n),stat=ierror)
!
!  If there is an error then return
!
      if (ierror.ne.0) then
        write(ioout,'(''!! Allocation error !!'')')
        return
      endif
      plocal = ' '
!
      if (associated(p)) then
!
!  Find size of old array
!
        nold = size(p)
!
!  Preserve data
!
        plocal(1:min(nold,n)) = p(1:min(nold,n))
!
!  Free old memory
!
        deallocate(p)
      endif
!
!  Re-assign pointers so that p is returned as the new array
!
      p=>plocal

    end subroutine realloc_kim_len_1

    subroutine realloc_kim_len_2(p,m,n,ierror)

      implicit none
!
!  Passed arguments
!
      integer(i4), intent(in)                         :: m,n
      integer(i4), intent(out)                        :: ierror
      character(len=kim_len), dimension(:,:), pointer :: p
!
!  Local arguments
!
      character(len=kim_len), dimension(:,:), pointer :: plocal
      integer(i4)                                     :: mold,nold
!
!  Allocate new size of memory
!
      nullify(plocal)
      allocate(plocal(m,n),stat=ierror)
!
!  If there is an error then return
!
      if (ierror.ne.0) then
        write(ioout,'(''!! Allocation error !!'')')
        return
      endif
      plocal = ' '
!
      if (associated(p)) then
!
!  Find size of old array
!
        mold = size(p,1)
        nold = size(p,2)
!
!  Preserve data
!
        plocal(1:min(mold,m),1:min(nold,n)) = &
          p(1:min(mold,m),1:min(nold,n))
!
!  Free old memory
!
        deallocate(p)
      endif
!
!  Re-assign pointers so that p is returned as the new array
!
      p=>plocal

    end subroutine realloc_kim_len_2

    subroutine realloc_kim_len_3(p,k,m,n,ierror)

      implicit none
!
!  Passed arguments
!
      integer(i4), intent(in)                           :: k,m,n
      integer(i4), intent(out)                          :: ierror
      character(len=kim_len), dimension(:,:,:), pointer :: p
!
!  Local arguments
!
      character(len=kim_len), dimension(:,:,:), pointer :: plocal
      integer(i4)                                       :: mold,nold
      integer(i4)                                       :: kold
!
!  Allocate new size of memory
!
      nullify(plocal)
      allocate(plocal(k,m,n),stat=ierror)
!
!  If there is an error then return
!
      if (ierror.ne.0) then
        write(ioout,'(''!! Allocation error !!'')')
        return
      endif
      plocal = ' '
!
      if (associated(p)) then
!
!  Find size of old array
!
        kold = size(p,1)
        mold = size(p,2)
        nold = size(p,3)
!
!  Preserve data
!
        plocal(1:min(kold,k),1:min(mold,m),1:min(nold,n)) = &
          p(1:min(kold,k),1:min(mold,m),1:min(nold,n))
!
!  Free old memory
!
        deallocate(p)
      endif
!
!  Re-assign pointers so that p is returned as the new array
!
      p=>plocal

    end subroutine realloc_kim_len_3

    subroutine realloc_kim_len_4(p,k,l,m,n,ierror)

      implicit none
!
!  Passed arguments
!
      integer(i4), intent(in)                             :: k,l,m,n
      integer(i4), intent(out)                            :: ierror
      character(len=kim_len), dimension(:,:,:,:), pointer :: p
!
!  Local arguments
!
      character(len=kim_len), dimension(:,:,:,:), pointer :: plocal
      integer(i4)                                         :: mold,nold
      integer(i4)                                         :: kold,lold
!
!  Allocate new size of memory
!
      nullify(plocal)
      allocate(plocal(k,l,m,n),stat=ierror)
!
!  If there is an error then return
!
      if (ierror.ne.0) then
        write(ioout,'(''!! Allocation error !!'')')
        return
      endif
      plocal = ' '
!
      if (associated(p)) then
!
!  Find size of old array
!
        kold = size(p,1)
        lold = size(p,2)
        mold = size(p,3)
        nold = size(p,4)
!
!  Preserve data
!
        plocal(1:min(kold,k),1:min(lold,l),1:min(mold,m), &
          1:min(nold,n)) = p(1:min(kold,k),1:min(lold,l), &
          1:min(mold,m),1:min(nold,n))
!
!  Free old memory
!
        deallocate(p)
      endif
!
!  Re-assign pointers so that p is returned as the new array
!
      p=>plocal

    end subroutine realloc_kim_len_4
#endif

  end module reallocate_kim
