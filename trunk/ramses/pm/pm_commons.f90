module pm_commons
  use amr_parameters
  use pm_parameters

  ! Particles related arrays
  real(dp),allocatable,dimension(:,:)       ::xp       ! Positions
  real(dp),allocatable,dimension(:,:)       ::vp       ! Velocities
  real(dp),allocatable,dimension(:)         ::mp       ! Masses

  integer(kind=8),allocatable,dimension(:,:), target :: part_hkey
  integer(kind=4),allocatable,dimension(:  ), target :: current_state
  integer(kind=4),allocatable,dimension(:)  ::part_ind_permutation, part_ind_permutation2
  ! Particle histogram related variables and arrays
  integer(kind=8),allocatable,dimension(:,:)::bin_keys,particle_histogram_keys
  real(dp),allocatable,dimension(:)         ::bin_mass,particle_histogram_mass, bin_count
  integer                                   ::nbins

#ifdef OUTPUT_PARTICLE_POTENTIAL
  real(dp),allocatable,dimension(:)  ::ptcl_phi ! Potential of particle added by AP for output purposes 
#endif
  integer ,allocatable,dimension(:)  ::nextp    ! Next particle in list
  integer ,allocatable,dimension(:)  ::prevp    ! Previous particle in list
  integer ,allocatable,dimension(:)  ::levelp   ! Current level of particle
  integer(id_pre),allocatable,dimension(:)::idp    ! Identity of particle
  integer ,allocatable,dimension(:)  :: part_level_offset    
  integer ,allocatable,dimension(:)  :: bin_start_offset


  ! Tree related arrays
  integer ,allocatable,dimension(:)  ::headp    ! Head particle in grid
  integer ,allocatable,dimension(:)  ::tailp    ! Tail particle in grid
  integer ,allocatable,dimension(:)  ::numbp    ! Number of particles in grid

  ! Global particle linked lists
  integer::headp_free,tailp_free,numbp_free=0,numbp_free_tot=0

  real(dp),allocatable,dimension(:)::dp_part_send_buf,dp_part_recv_buf


  ! Operator overloading for arry_pop(arr, ipop)
  ! Remove element ipop from array and move elements > ipop one position to the left.
  interface array_pop
     module procedure array_pop_i8, array_pop_i4, array_pop_r8, array_pop_r4, array_pop_l
  end interface array_pop
  
contains

! #############################################################
  
  subroutine kill_one_particle(ipart)
    integer, intent(in) :: ipart
    
    ! Kill particle ipart by sliding all particles right
    ! of ipart one position to the left.
    ! Put here in the hope that poeple adding new particle
    ! attributes won't forget to adapt this routine...
    
    integer :: ihilbert, idim

    npart = npart - 1
    do idim = 1, ndim
       call array_pop(xp(1:npart, idim), ipart)
       call array_pop(vp(1:npart, idim), ipart)
    end do
    call array_pop(mp(1:npart), ipart)
    call array_pop(idp(1:npart), ipart)
    call array_pop(levelp(1:npart), ipart)
    do ihilbert = 1, nhilbert
       call array_pop(part_hkey(1:npart, ihilbert), ipart)       
    end do
    call array_pop(current_state(1:npart), ipart)
#ifdef OUTPUT_PARTICLE_POTENTIAL
    call array_pop(ptcl_phi(1:npart), ipart)
#endif
  end subroutine kill_one_particle

! #############################################################

! #############################################################
  
  ! TODO: Can this stupid copying be avoided?
  subroutine array_pop_i4(arr, ipop)
    implicit none
    integer(4), dimension(:), intent(inout) :: arr
    integer, intent(in) :: ipop
    integer :: i
    do i = ipop, size(arr) - 1
       arr(i) = arr(i + 1)
    end do
  end subroutine array_pop_i4
  
  subroutine array_pop_i8(arr, ipop)
    implicit none
    integer(8), dimension(:), intent(inout) :: arr
    integer, intent(in) :: ipop
    integer :: i
    do i = ipop, size(arr) - 1
       arr(i) = arr(i + 1)
    end do
  end subroutine array_pop_i8
  
  subroutine array_pop_r4(arr, ipop)
    implicit none
    real(4), dimension(:), intent(inout) :: arr
    integer, intent(in) :: ipop
    integer :: i
    do i = ipop, size(arr) - 1
       arr(i) = arr(i + 1)
    end do
  end subroutine array_pop_r4
  
  subroutine array_pop_r8(arr, ipop)
    implicit none
    real(8), dimension(:), intent(inout) :: arr
    integer, intent(in) :: ipop
    integer :: i
    do i = ipop, size(arr) - 1
       arr(i) = arr(i + 1)
    end do
  end subroutine array_pop_r8
  
  subroutine array_pop_l(arr, ipop)
    implicit none
    logical, dimension(:), intent(inout) :: arr
    integer, intent(in) :: ipop
    integer :: i
    do i = ipop, size(arr) - 1
       arr(i) = arr(i + 1)
    end do
  end subroutine array_pop_l

! #############################################################

end module pm_commons
