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
  use pm_commons,     only: xp, vp, part_level_offset
  use amr_parameters, only: dp, ndim, nx, ny, nz, boxlen
  use amr_commons,    only: dtnew, icoarse_max, icoarse_min, jcoarse_min, kcoarse_min, &
                            nx, ny, nz
  implicit none

  integer, intent(in) :: ilevel

  integer :: idim, ipart, offset, nparts, nx_loc
  real(dp), dimension(1:3) :: skip_loc, xbound
  real(dp) :: scale
  
  xbound(1:3) = (/dble(nx), dble(ny), dble(nz)/)
  nx_loc = (icoarse_max - icoarse_min + 1)
  if(ndim > 0) skip_loc(1) = dble(icoarse_min)
  if(ndim > 1) skip_loc(2) = dble(jcoarse_min)
  if(ndim > 2) skip_loc(3) = dble(kcoarse_min)
  scale = boxlen / dble(nx_loc)
  
  offset = part_level_offset(ilevel)
  nparts = part_level_offset(ilevel + 1) - part_level_offset(ilevel)
  
  ! Update position
  do idim = 1, ndim
     do ipart = offset + 1, offset + nparts 
        xp(ipart, idim) = xp(ipart, idim) &
             + vp(ipart, idim) * dtnew(ilevel)
     end do
  end do

  ! Fix periodic boundaries
  do idim = 1, ndim
     do ipart = offset + 1, offset + nparts
        if (xp(ipart, idim) / scale + skip_loc(idim) < 0.0d0) &
             & xp(ipart, idim) = xp(ipart, idim) + (xbound(idim) - skip_loc(idim)) * scale
        if (xp(ipart, idim) / scale + skip_loc(idim) >= xbound(idim)) &
             & xp(ipart,idim) = xp(ipart, idim) - (xbound(idim) - skip_loc(idim)) * scale
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
  integer :: np, offset
  
  offset = part_level_offset(ilevel)
  np = part_level_offset(ilevel + 1) - part_level_offset(ilevel)
  
  levelp(offset + 1:offset + np) = ilevel
  
end subroutine update_levelp
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine compute_particle_acceleration(ap, offset, nparts, ilevel, read_gas_velocity)
  use pm_commons,      only: part_level_offset, xp, &
                             part_hkey, npart
  use amr_parameters,  only: dp, nvector, ndim, twotondim, poisson, verbose, nhilbert
  use hydro_commons,   only: uold
  use poisson_commons, only: f
  use amr_commons,     only: dtnew, ncpu, myid, t, son
  use pm_parameters,   only: npartmax
#ifndef WITHOUTMPI
  use particle_communication, only: build_communicator, part_data_to_domain, domain_data_to_part
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
  ! several times and not just once for all parts from offset + 1 to offset_nparts

  if(verbose)write(*,'("Entering compute_particle_acceleration, level " I2)')ilevel 

#ifndef WITHOUTMPI
  call build_communicator(communicator, npart_recv, &
       nparts, nparts_local, local_oft, &
       part_hkey(offset + 1 : offset + nparts, 1:nhilbert), & 
       ilevel)

  allocate(xp_remote(1:npart_recv, 1:ndim), ap_remote(1:npart_recv, 1:ndim))
  do idim = 1, ndim
     call part_data_to_domain(communicator, xp(offset + 1 : offset + nparts, idim), xp_remote(:, idim))
  end do
  
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
  do idim = 1, ndim
     call domain_data_to_part(communicator, ap_remote(:,idim), ap(offset + 1 : offset + nparts, idim))
  end do

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
