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




end module pm_commons
