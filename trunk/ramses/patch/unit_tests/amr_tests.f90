subroutine amr_tests(all_ok)
  implicit none
  logical::all_ok
  logical::amr_ok=.true.

  call hilbert_tests(amr_ok)
  call hash_tests(amr_ok, .false.)
!  call other_test2(amr_ok)

  if(.not. amr_ok)then
     write(*,*)'AMR_TESTS FAILED'
     all_ok=.false.
  end if
end subroutine amr_tests


! =====================================================================================
! =====================================================================================
! ADD UNIT TESTS HERE
! =====================================================================================
! =====================================================================================

subroutine hilbert_tests(all_ok)
  use amr_parameters, only : ndim, nhilbert, int_pre, nvector
  use hilbert,       only: hilbert_nd_reverse, hilbert_nd, bits_per_int
  implicit none

  logical, intent(inout) :: all_ok
  ! A simple test which transforms integer coordinates into hilbert keys
  ! and the hilbert key back into integer coordinate. Results must be equal to input.  

  integer :: ilevel, maxlevel, i, idim
  logical :: ok, ok_test
  real, dimension(1:nvector, 1:ndim) :: xx
  integer(int_pre),dimension(1:nvector, 1:ndim) :: ix, ix_store
  integer(kind=8),dimension(1:nvector, 1:nhilbert) :: hkey
  integer(kind=4),dimension(1:nvector) :: cstate

  maxlevel = nhilbert * bits_per_int(ndim) / ndim 

  ok_test = .true.

  do ilevel = 1, maxlevel
     ok=.true.

     call random_number(xx)

     do idim = 1, ndim
        do i = 1, nvector
           ix(i, idim) = int(xx(i, idim) * 2.0**ilevel, kind=8)
        end do
     end do
     ix_store = ix
     call hilbert_nd(ix, hkey, cstate, 0, ilevel, nvector)
     call hilbert_nd_reverse(ix, hkey, ilevel, nvector)
     
     do idim = 1, ndim
        do i=1,nvector
           if( ix_store(i, idim) .ne. ix(i, idim) )then
              ok = .false.
              ok_test = .false.
              write(*,*)ix_store(i, idim), ix(i, idim)
           end if
        end do
     end do
     
     if (.not. ok)then
        write(*,*)'hilbert test FAILED for ilevel ', ilevel
        all_ok = .false.
        ok_test = .false.
     end if
  end do

  if (ok_test) write(*,*)'hilbert test passed.'

end subroutine hilbert_tests
! =====================================================================================
! =====================================================================================

subroutine hash_tests(all_ok, verbose)
  use hash
  implicit none

  logical, intent(inout) :: all_ok
  logical, intent(in) :: verbose

  type(hash_table)::htable
  integer::i,nfree_store,nfree_chain_store
  logical::ok
  real,dimension(1:3000)::val_float
  real,dimension(0:ndim,1:3000)::key_float
  integer,dimension(1:3000)::val
  integer(int_pre),dimension(0:ndim,1:3000)::key
  character(6) :: hash_type = 'simple'
  
  ok=.true.
 
  call random_number(key_float)
  call random_number(val_float)
  
  do i=1,3000
     key(0:ndim,i)=int(key_float(0:ndim,i)*2.0**12,kind=8)
     val(i)=int(val_float(i)*2000,kind=4) + 1
  end do

  call init_empty_hash(htable,3000, hash_type)

  if(verbose) call hash_stats(htable)

  nfree_store=htable%nfree
  nfree_chain_store=htable%nfree_chain
  
  do i=1,2000
     call hash_set(htable,key(0:ndim,i),val(i))     
  end do

  if(verbose) call hash_stats(htable)


  do i=2000,1,-1
     ok=ok .and. (val(i)==hash_get(htable,key(0:ndim ,i) ))
  end do


  do i=1001,2000
     call hash_free(htable,key(0:ndim ,i) )
  end do

  if(verbose) call hash_stats(htable)

  do i=2001,3000
     call hash_set(htable,key(0:ndim ,i) ,val(i))
  end do
  
  if(verbose) call hash_stats(htable)

  do i=1,1000
     ok=ok .and. (val(i)==hash_get(htable,key(0:ndim ,i) ))
  end do

  do i=2001,3000
     ok=ok .and. (val(i)==hash_get(htable,key(0:ndim ,i) ))
     call hash_free(htable,key(0:ndim ,i) )
  end do

  if(verbose) call hash_stats(htable)

  do i=1000,1,-1
     ok=ok .and. (val(i)==hash_get(htable,key(0:ndim ,i) ))
     call hash_free(htable,key(0:ndim ,i) )
  end do
  
  if(verbose) call hash_stats(htable)
  
  ok=ok .and. (nfree_store==htable%nfree)
  ok=ok .and. (nfree_chain_store==htable%nfree_chain)

  
  if (.not. ok)then
     write(*,*)'hash test FAILED '
     all_ok=.false.
  else
     write(*,*)'hash test passed '
  end if
  
end subroutine hash_tests
