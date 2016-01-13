!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine kick(ilevel, previous_timestep)
  use pm_commons,      only: part_level_offset, xp, vp, levelp
  use amr_parameters,  only: dp, ndim, tracer, hydro, static, nvector
  use amr_commons,     only: dtnew, dtold
  implicit none

  integer, intent(in) :: ilevel
  logical, intent(in), value :: previous_timestep


  integer :: offset, nparts, idim, ioft, np, ip
  real(dp), allocatable, dimension(:,:) :: ap
  real(dp), dimension(1:nvector) :: dteff
  
  offset = part_level_offset(ilevel)
  nparts = part_level_offset(ilevel + 1) - part_level_offset(ilevel)

  allocate(ap(offset + 1: offset + nparts, 1:ndim))
  
  call compute_particle_acceleration(ap, offset, nparts, ilevel, tracer .and. hydro)

  dteff = 0.5d0 * dtnew(ilevel)
  
  do ioft = offset, offset + nparts - 1, nvector
     np = min(nvector, offset + nparts - ioft)
     
     ! Compute individual time steps
     if (previous_timestep)then
        do ip = 1, np
           if(levelp(ioft + ip) >= ilevel)then
              dteff(ip) = 0.5d0 * dtnew(levelp(ioft + ip))
           else
              dteff(ip) = 0.5d0 * dtold(levelp(ioft + ip))
           endif
        end do
     end if
     
     ! Update velocity     
     ! TODO: fix static/tracer cases
     do idim = 1, ndim     
        do ip = 1, np
           vp(ioft + ip, idim) = vp(ioft + ip, idim) &
                + ap(ioft + ip, idim) * dteff(ip)
        end do
     end do
  end do
  deallocate(ap)
  
end subroutine kick
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine drift(ilevel)
  use pm_commons, only: xp, vp, part_level_offset
  use amr_parameters, only: dp, ndim, nx, ny, nz, boxlen
  use amr_commons, only: dtnew, period
  implicit none

  integer, intent(in) :: ilevel

  integer :: idim, ipart, offset, nparts

  offset = part_level_offset(ilevel)
  nparts = part_level_offset(ilevel + 1) - part_level_offset(ilevel)
  
  ! Update position
  do idim = 1, ndim
     do ipart = offset + 1, offset + nparts 
        xp(ipart, idim) = xp(ipart, idim) &
             + vp(ipart, idim) * dtnew(ilevel)
     end do
     
     ! Take care of boundary conditions
     ! TODO: non-periodic boundaries!!
     do ipart = offset + 1, offset + nparts
        if (xp(ipart, idim) > boxlen)then
           xp(ipart, idim) = xp(ipart, idim) - boxlen
        end if
        if(xp(ipart, idim) < 0.d0)then
           xp(ipart, idim) = xp(ipart, idim) + boxlen
        end if
     end do
  end do
end subroutine drift
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine update_levelp(ilevel)
  use pm_commons, only: levelp, part_level_offset
  implicit none

  integer, intent(in) :: ilevel
  integer :: ipart, nparts, offset

  offset = part_level_offset(ilevel)
  nparts = part_level_offset(ilevel + 1) - part_level_offset(ilevel)
  
  do ipart = offset + 1, offset + nparts 
     levelp(ipart) = ilevel
  end do

end subroutine update_levelp
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine compute_particle_acceleration(ap, offset, nparts, ilevel, read_gas_velocity)
  use pm_commons,      only: part_level_offset, xp, &
                             part_hkey, npart
  use amr_parameters,  only: dp, nvector, ndim, twotondim, poisson, verbose
  use hydro_commons,   only: uold
  use poisson_commons, only: f
  use amr_commons,     only: dtnew, ncpu, myid, t, son
  use pm_parameters,   only: npartmax
#ifndef WITHOUTMPI
  use particle_communication, only: build_communicator, part_data_to_domain_dp, domain_data_to_part_dp
#endif
  use hilbert,     only: hilbert_for_particle 
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h' 
  integer :: info
#endif



  integer, intent(in) :: ilevel, offset, nparts
  logical, intent(in) :: read_gas_velocity
  real(dp), intent(inout), dimension(offset + 1 : offset + nparts, 1:ndim) :: ap
  
  real(dp), allocatable, dimension(:,:) :: xp_remote, ap_remote
  integer,  dimension(1:ncpu, 1:4)       :: communicator
  integer,  dimension(1:nvector, 1:twotondim), save :: cell_index
  real(dp), dimension(1:nvector, 1:twotondim), save :: vol
  
  integer :: ioft, np, ip, ind, idim, ipart, local_oft, npart_recv, nparts_local

  ! TODO: consistent naming (np, nparts, npart) throughout routines
  ! TODO: try to avoid usage of big ap(1:npartmax) array. For example, sudivide ilevel and call routine
  ! several times and not just once for all parts from offset + 1 to offset_nparts

  if(verbose)write(*,'("Entering compute_particle_acceleration, level " I2)')ilevel 

#ifndef WITHOUTMPI
  call build_communicator(communicator, npart_recv, &
       nparts, nparts_local, local_oft, &
       part_hkey(offset + 1 : offset + nparts, 1), & 
#if NHILBERT > 1
       part_hkey(offset + 1 : offset + nparts, 2), &
#endif
#if NHILBERT > 2
       part_hkey(offset + 1 : offset + nparts, 3), &
#endif
       ilevel)

  allocate(xp_remote(1:npart_recv, 1:ndim), ap_remote(1:npart_recv, 1:ndim))
  call part_data_to_domain_dp(communicator, xp(offset + 1 : offset + nparts, 1), xp_remote(:, 1))
  if (ndim > 1) call part_data_to_domain_dp(communicator, xp(offset + 1 : offset + nparts, 2), xp_remote(:, 2))
  if (ndim > 2) call part_data_to_domain_dp(communicator, xp(offset + 1 : offset + nparts, 3), xp_remote(:, 3))

  
  ! Deal with remote particles
  do ioft = 0, npart_recv - 1, nvector
     np = min(nvector, npart_recv - ioft)
     
     call cic(xp_remote, npart_recv, cell_index, vol, ioft, np, ilevel, 2)

     ap_remote(ioft + 1: ioft + np, 1: ndim) = 0.0D0
     if(read_gas_velocity)then
        do idim = 1, ndim
           do ind = 1, twotondim              
              do ip = 1, np
                 ap_remote(ioft + ip, idim) = ap_remote(ioft + ip, idim) + uold(cell_index(ip,ind),idim+1) * vol(ip,ind)
              end do
           end do
        end do
     endif
     
     if(poisson)then
        do idim = 1,ndim
           do ind = 1,twotondim
              do ip = 1,np
                 ap_remote(ioft + ip, idim) = ap_remote(ioft + ip, idim) + f(cell_index(ip,ind),idim) * vol(ip,ind)
              end do
           end do
        end do
     endif
  end do
  call domain_data_to_part_dp(communicator, ap_remote(:,1), ap(offset + 1 : offset + nparts, 1))
  if (ndim > 1) call domain_data_to_part_dp(communicator, ap_remote(:,2), ap(offset + 1 : offset + nparts, 2))
  if (ndim > 2) call domain_data_to_part_dp(communicator, ap_remote(:,3), ap(offset + 1 : offset + nparts, 3))
  deallocate(xp_remote, ap_remote)

#endif
  
  ! Deal with local particles
  do ioft = offset + local_oft, offset + local_oft + nparts_local - 1, nvector
     np = min(nvector, offset + local_oft + nparts_local - ioft)

     call cic(xp, npartmax, cell_index, vol, ioft, np, ilevel, 2)
     
     ap(ioft + 1: ioft + np, 1: ndim) = 0.0D0
     if(read_gas_velocity)then
        do idim = 1, ndim
           do ind = 1, twotondim              
              do ip = 1, np
                 ap(ioft + ip, idim) = ap(ioft + ip, idim) + uold(cell_index(ip,ind),idim+1) * vol(ip,ind)
              end do
           end do
        end do
     endif
     
     if(poisson)then
        do idim = 1,ndim
           do ind = 1,twotondim
              do ip = 1,np
                  ap(ioft + ip, idim) =  ap(ioft + ip, idim) + f(cell_index(ip,ind),idim) * vol(ip,ind)
               end do
           end do
        end do
     endif
  end do
  
#ifdef OUTPUT_PARTICLE_POTENTIAL
  ! Just a reminder that this option is not built in yet
  print*,"stopping because of OUTPUT_PARTICLE_POTENTIAL"
  call clean_stop
#endif
end subroutine compute_particle_acceleration
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
