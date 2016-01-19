subroutine pm_tests(all_ok)
  use pm_commons, only: xp, vp, mp, idp, levelp, part_hkey, current_state, part_ind_permutation, part_ind_permutation2, npart
  use amr_parameters, only: boxlen, ndim, nhilbert
  implicit none
  logical::all_ok
  logical::pm_ok=.true.
  integer, parameter :: size = 331



  ! a little bit of init_part
  allocate(xp    (size,ndim))
  allocate(vp    (size,ndim))
  vp=0.
  allocate(mp    (size))
  mp=1
  allocate(levelp(size))
  levelp=1
  allocate(idp   (size))
  allocate(part_hkey(size,1:nhilbert))
  allocate(current_state(size))
  allocate(part_ind_permutation(size))
  allocate(part_ind_permutation2(size))
  npart=size
  boxlen=1





  call sort_particle_tests(pm_ok)
!  call other_test(pm_ok)

  if(.not. pm_ok)then
     write(*,*)'PM_TESTS FAILED'
     all_ok=.false.
  end if

contains

! =====================================================================================
! =====================================================================================
! ADD UNIT TESTS HERE
! =====================================================================================
! =====================================================================================
  subroutine sort_particle_tests(all_ok)
    use hilbert,       only: hilbert_for_particle, bits_per_int
    use sort,          only: lsd_radix_sort_particles, gt_keys, apply_particle_permutation
    implicit none


    logical, intent(inout) :: all_ok

    logical :: ok_test
    integer, parameter :: offs=117
    integer :: ilevel, maxlevel, i

    maxlevel = nhilbert * bits_per_int(ndim) / ndim 

    ok_test = .true.
    
    do ilevel = 1, maxlevel
       call random_number(xp)
       call hilbert_for_particle(offs, size-offs, 0, ilevel)
       call lsd_radix_sort_particles(offs, size-offs, ilevel, ilevel, .true.)
       call apply_particle_permutation(offs, size-offs, ilevel)
       do i = offs + 1, size - 1
          if (gt_keys(part_hkey(i, 1:nhilbert),part_hkey(i + 1, 1:nhilbert)))then
             write(*,*)'particle sort test FAILED for level', ilevel, ilevel
             write(*,*)part_hkey(i, 1:nhilbert)
             write(*,*)part_hkey(i + 1, 1:nhilbert)
             all_ok=.false.
             ok_test = .false.
          end if
       end do
    end do

    if (ok_test) write(*,*)'sort particle test passed' 
    
  end subroutine sort_particle_tests
  
end subroutine pm_tests
            
