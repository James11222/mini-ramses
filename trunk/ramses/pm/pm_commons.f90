module pm_commons
  use amr_parameters
  use pm_parameters

  real(dp),allocatable,dimension(:,:)       ::xp       ! Positions
  real(dp),allocatable,dimension(:,:)       ::vp       ! Velocities
  real(dp),allocatable,dimension(:)         ::mp       ! Masses
#ifdef OUTPUT_PARTICLE_POTENTIAL
  real(dp),allocatable,dimension(:)         ::phip     ! Potential
#endif
  integer ,allocatable,dimension(:)         ::levelp   ! Current level of particle
  integer(i8b),allocatable,dimension(:)     ::idp      ! Particle unique identifier
  integer ,allocatable,dimension(:)         ::sortp    ! Sorted indices
  integer ,allocatable,dimension(:)         ::workp    ! Work space

  integer ,allocatable,dimension(:)         ::headp    ! Particle levels head
  integer ,allocatable,dimension(:)         ::tailp    ! Particle levels tail

  real(dp)::mp_min=-1.0   ! Minimum particle mass

  integer :: action_kick_only = 1
  integer :: action_kick_drift = 2  

  integer ,dimension(1:MAXLEVEL):: part_patch_acc=3 !
  integer ,dimension(1:MAXLEVEL):: part_patch_rho=3 !
  integer ,dimension(1:MAXLEVEL):: part_patch_ref=3 !
  

  contains
    subroutine swap_parts(istart, iend, swap_table)
      use amr_parameters, only: ndim
      implicit none
      integer, intent(in) :: istart, iend
      integer, dimension(istart: ) :: swap_table

      ! Swap particles using new index table
      real(dp),dimension(1:ndim) :: xp_tmp, vp_tmp
      real(dp) :: mp_tmp
      integer :: levelp_tmp
      integer(i8b) :: idp_tmp
      integer :: ipart, jpart
      
      do ipart = istart, iend
         do while(swap_table(ipart) .NE. ipart)
            ! Swap new index
            jpart=swap_table(ipart)
            swap_table(ipart)=swap_table(jpart)
            swap_table(jpart)=jpart
            ! Swap positions
            xp_tmp(1:ndim)=xp(ipart,1:ndim)
            xp(ipart,1:ndim)=xp(jpart,1:ndim)
            xp(jpart,1:ndim)=xp_tmp(1:ndim)
            ! Swap velocities
            vp_tmp(1:ndim)=vp(ipart,1:ndim)
            vp(ipart,1:ndim)=vp(jpart,1:ndim)
            vp(jpart,1:ndim)=vp_tmp(1:ndim)
            ! Swap masses
            mp_tmp=mp(ipart)
            mp(ipart)=mp(jpart)
            mp(jpart)=mp_tmp
            ! Swap ids
            idp_tmp=idp(ipart)
            idp(ipart)=idp(jpart)
            idp(jpart)=idp_tmp
            ! Swap levels
            levelp_tmp=levelp(ipart)
            levelp(ipart)=levelp(jpart)
            levelp(jpart)=levelp_tmp
         end do
      end do
    end subroutine swap_parts
    subroutine apply_particle_permutation(istart, iend, build_workp)
      use amr_commons,    only: dp, ndim
      implicit none
      integer, intent(in)                          :: istart, iend
      logical, intent(in), value                   :: build_workp

      integer                                      :: ipart, idim
      
      real(dp),        allocatable, dimension(:)   :: extra_storage_dp
      integer(kind=8), allocatable, dimension(:)   :: extra_storage_i8
      integer(kind=4), allocatable, dimension(:)   :: extra_storage_i4
      
      
      ! The permutation sigma built during the radix sort sorts the index array I such that
      ! key(sigma(i)) is sorted. We thus have to invert sigma and apply
      ! this to the key array to get (sigma^-1(key))(i) sorted in memory

      if (build_workp)then
         do ipart = istart, iend
            workp( sortp(ipart) ) = ipart
         end do
      end if
      
      ! Apply  particle_permutation2
      
      ! Rearrange float arrays
      allocate(extra_storage_dp(istart: iend))
      
      do idim = 1, ndim       
         do ipart = istart, iend
            extra_storage_dp(workp(ipart)) = xp(ipart, idim) 
         end do
         xp(istart: iend, idim) = extra_storage_dp(istart: iend)
      end do
      
      do idim = 1, ndim       
         do ipart = istart, iend
            extra_storage_dp(workp(ipart)) = vp(ipart, idim) 
         end do
         vp(istart: iend, idim) = extra_storage_dp(istart: iend)
      end do
      
      do ipart = istart, iend
         extra_storage_dp(workp(ipart)) = mp(ipart) 
      end do
      mp(istart: iend) = extra_storage_dp(istart: iend)

      deallocate(extra_storage_dp)

      ! Rearrange long integer arrays

#if ID_PRECISION == 8
      allocate(extra_storage_i8(istart: iend))
      do ipart = istart, iend
         extra_storage_i8(workp(ipart)) = idp(ipart)
      end do
      idp(istart: iend) = extra_storage_i8(istart: iend)
      deallocate(extra_storage_i8)
#endif



      ! Rearrange short integer arrays
      allocate(extra_storage_i4(istart: iend))
#if ID_PRECISION == 4
      do ipart = istart, iend
         extra_storage_i4(workp(ipart)) = idp(ipart)
      end do
      idp(istart: iend) = extra_storage_i4(istart: iend)
#endif

      do ipart = istart, iend
         extra_storage_i4(workp(ipart)) = levelp(ipart)
      end do
      levelp(istart: iend) = extra_storage_i4(istart: iend)

      deallocate(extra_storage_i4)

      !    ! Reset applied part of permutation
      !    do ipart = offset + 1, offset + np
      !       sortp(ipart) = ipart
      !    end do

    end subroutine apply_particle_permutation




end module pm_commons
