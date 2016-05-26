!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine rho_fine(ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer :: ilevel
  !------------------------------------------------------------------
  ! ADD NEW DESCRIPTION HERE!
  !------------------------------------------------------------------
  integer :: particle_level, offset, grid_level, part_level, np
  integer::iskip,icpu,ind,i,info,nx_loc,ibound,idim,icell, ncell, ilev
  real(dp)::dx,d_scale,scale,dx_loc,scalar
  real(dp)::d0,m_refine_loc,dx_min,vol_min,mstar,msnk,nISM,nCOM
  real(kind=8)::total,total_all,total2,total2_all,tms
  real(kind=8),dimension(2)::totals_in,totals_out
  logical::multigrid=.false., ok, first
  real(kind=8),dimension(1:ndim+1)::multipole_in,multipole_out

  interface
     subroutine mass_deposit(xpart, mpart, nparts, grid_level, nbits_patch)
       use amr_parameters, only: ndim, dp, int_pre
       use amr_commons,    only: ncpu, ind_table2, boxlen
       implicit none
       integer, intent(in) :: grid_level, nparts, nbits_patch
       real(dp), dimension(:, :), intent(inout) :: xpart
       real(dp), dimension(:), intent(in) :: mpart
     end subroutine mass_deposit
  end interface

  if(.not. poisson)return
  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Mesh spacing in that level
  dx=0.5D0**ilevel 
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  if(ilevel==levelmin)multipole=0d0

  !-------------------------------------------------------
  ! Initialize rho to analytical and baryon density field
  !-------------------------------------------------------
  do i=nlevelmax,ilevel,-1
     ! Compute mass multipole
     if(hydro)call multipole_fine(i)
     ! Perform CIC using pseudo-particle
     call cic_from_multipole(i)
     ! Update boundaries
     call make_virtual_reverse_dp(rho(1),i)
     call make_virtual_fine_dp   (rho(1),i)
  end do

  do ilev = ilevel, nlevelmax
     !--------------------------
     ! Initialize fields to zero
     !--------------------------
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,active(ilev)%ngrid
           phi(active(ilev)%igrid(i)+iskip)=0.0D0
        end do
     end do

     !-------------------------------------------------------------------------
     ! Initialize "number density" field to baryon number density in array phi.
     !-------------------------------------------------------------------------
     if(m_refine(ilev)>-1.0d0)then
        d_scale=max(mass_sph/dx_loc**ndim,smallr)
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           if(hydro)then
              if(ivar_refine>0)then
                 do i=1,active(ilev)%ngrid
                    scalar=uold(active(ilev)%igrid(i)+iskip,ivar_refine) &
                         & /uold(active(ilev)%igrid(i)+iskip,1)
                    if(scalar>var_cut_refine)then
                       phi(active(ilev)%igrid(i)+iskip)= &
                            & rho(active(ilev)%igrid(i)+iskip)/d_scale
                    endif
                 end do
              else
                 do i=1,active(ilev)%ngrid
                    phi(active(ilev)%igrid(i)+iskip)= &
                         & rho(active(ilev)%igrid(i)+iskip)/d_scale
                 end do
              endif
           endif
        end do
     endif

     !-------------------------------------------------------
     ! Initialize rho and phi to zero in virtual boundaries
     !-------------------------------------------------------
     do icpu=1,ncpu
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,reception(icpu,ilev)%ngrid
              rho(reception(icpu,ilev)%igrid(i)+iskip)=0.0D0
              phi(reception(icpu,ilev)%igrid(i)+iskip)=0.0D0
           end do
        end do
     end do
  end do
  !---------------------------------------------------------
  ! Compute particle contribution to density field
  !---------------------------------------------------------
  ! Compute density due to current level particles

  if(pic)then
     do grid_level = ilevel, nlevelmax
        offset = part_level_offset(grid_level)
        call mass_deposit(xp(offset + 1: npart, 1:ndim), mp(offset + 1: npart), npart - offset, grid_level, 3)
     end do
  end if

  do particle_level = ilevel, nlevelmax
     call make_virtual_reverse_dp(rho(1),particle_level)
     call make_virtual_fine_dp   (rho(1),particle_level)
     if(m_refine(particle_level)>-1.0d0)then
        call make_virtual_reverse_dp(phi(1),particle_level)
        call make_virtual_fine_dp   (phi(1),particle_level)
     endif
  end do

  if (ilevel==levelmin) then
     call add_particle_multipole
  end if
  !--------------------------------------------------------------
  ! Compute multipole contribution from all cpus and set rho_tot
  !--------------------------------------------------------------
#ifndef WITHOUTMPI
  if(ilevel==levelmin)then
     multipole_in=multipole
     call MPI_ALLREDUCE(multipole_in,multipole_out,ndim+1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     multipole=multipole_out
  endif
#endif
  if(nboundary==0)then
     rho_tot=multipole(1)/scale**ndim
     if(debug)write(*,*)'rho_average=',rho_tot
  else
     rho_tot=0d0
  endif

  do particle_level = ilevel, nlevelmax
     !----------------------------------------------------
     ! Reset rho and phi in physical boundaries
     !----------------------------------------------------
     do ibound=1,nboundary
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,boundary(ibound,particle_level)%ngrid
              phi(boundary(ibound,particle_level)%igrid(i)+iskip)=0.0
              rho(boundary(ibound,particle_level)%igrid(i)+iskip)=0.0
           end do
        end do
     end do

     !-----------------------------------------
     ! Compute quasi Lagrangian refinement map
     !-----------------------------------------
     if(m_refine(particle_level)>-1.0d0)then
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,active(particle_level)%ngrid
              if(phi(active(particle_level)%igrid(i)+iskip)>=m_refine(particle_level))then
                 cpu_map2(active(particle_level)%igrid(i)+iskip)=1
              else
                 cpu_map2(active(particle_level)%igrid(i)+iskip)=0
              end if
           end do
        end do
        ! Update boundaries
        call make_virtual_fine_int(cpu_map2(1),particle_level)
     end if
  end do

111 format('   Entering rho_fine for level ',I2)
  
end subroutine rho_fine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine multipole_fine(ilevel)
  use amr_commons
  use hydro_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  !-------------------------------------------------------------------
  ! This routine compute array rho (source term for Poisson equation)
  ! by first reseting array rho to zero, then 
  ! by affecting the gas density to leaf cells, and finally
  ! by performing a restriction operation for split cells.
  ! For pure particle runs, the restriction is not necessary and the
  ! routine only set rho to zero. On the other hand, for the Multigrid
  ! solver, the restriction is necessary in any case.
  !-------------------------------------------------------------------
  integer ::ind,i,icpu,ncache,igrid,ngrid,iskip,info,ibound,nx_loc
  integer ::idim,nleaf,nsplit,ix,iy,iz,iskip_son,ind_son,ind_grid_son,ind_cell_son
  integer,dimension(1:nvector),save::ind_grid,ind_cell,ind_leaf,ind_split
  real(dp),dimension(1:nvector,1:ndim),save::xx
  real(dp),dimension(1:nvector),save::dd
  real(kind=8)::vol,dx,dx_loc,scale,vol_loc,mm
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Mesh spacing in that level
  dx=0.5D0**ilevel 
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc=dx_loc**ndim
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do

  ! Initialize fields to zero
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,active(ilevel)%ngrid
        unew(active(ilevel)%igrid(i)+iskip,1)=0.0D0
     end do
     do idim=1,ndim
        do i=1,active(ilevel)%ngrid
           unew(active(ilevel)%igrid(i)+iskip,idim+1)=0.0D0
        end do
     end do
  end do

  ! Compute mass multipoles in each cell
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     
     ! Loop over cells
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        ! Gather cell indices
        do i=1,ngrid
           ind_cell(i)=ind_grid(i)+iskip
        end do

        ! Gather leaf cells and compute cell centers
        nleaf=0
        do i=1,ngrid
           if(son(ind_cell(i))==0)then
              nleaf=nleaf+1
              ind_leaf(nleaf)=ind_cell(i)
              do idim=1,ndim
                 xx(nleaf,idim)=(xg(ind_grid(i),idim)+xc(ind,idim)-skip_loc(idim))*scale
              end do
           end if
        end do
        
        ! Compute gas multipole for leaf cells only
        if(hydro)then
           do i=1,nleaf
              mm=max(uold(ind_leaf(i),1),smallr)*vol_loc
              unew(ind_leaf(i),1)=unew(ind_leaf(i),1)+mm
           end do
           do idim=1,ndim
              do i=1,nleaf
                 mm=max(uold(ind_leaf(i),1),smallr)*vol_loc
                 unew(ind_leaf(i),idim+1)=unew(ind_leaf(i),idim+1)+mm*xx(i,idim)
              end do
           end do
        endif

        ! Add analytical density profile for leaf cells only
        if(gravity_type < 0)then           
           ! Call user defined routine rho_ana
           call rho_ana(xx,dd,dx_loc,nleaf)
           ! Scatter results to array phi
           do i=1,nleaf
              unew(ind_leaf(i),1)=unew(ind_leaf(i),1)+dd(i)*vol_loc
           end do
           do idim=1,ndim
              do i=1,nleaf
                 mm=dd(i)*vol_loc
                 unew(ind_leaf(i),idim+1)=unew(ind_leaf(i),idim+1)+mm*xx(i,idim)
              end do
           end do           
        end if

        ! Gather split cells
        nsplit=0
        do i=1,ngrid
           if(son(ind_cell(i))>0)then
              nsplit=nsplit+1
              ind_split(nsplit)=ind_cell(i)
           end if
        end do

        ! Add children multipoles
        do ind_son=1,twotondim
           iskip_son=ncoarse+(ind_son-1)*ngridmax
           do i=1,nsplit
              ind_grid_son=son(ind_split(i))
              ind_cell_son=iskip_son+ind_grid_son
              unew(ind_split(i),1)=unew(ind_split(i),1)+unew(ind_cell_son,1)
           end do
           do idim=1,ndim
              do i=1,nsplit
                 ind_grid_son=son(ind_split(i))
                 ind_cell_son=iskip_son+ind_grid_son
                 unew(ind_split(i),idim+1)=unew(ind_split(i),idim+1)+unew(ind_cell_son,idim+1)
              end do
           end do
        end do

     end do
  enddo

  ! Update boundaries
  do idim=1,ndim+1
     call make_virtual_fine_dp(unew(1,idim),ilevel)
  end do

111 format('   Entering multipole_fine for level',i2)

end subroutine multipole_fine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cic_from_multipole(ilevel)
  use amr_commons
  use hydro_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  logical::multigrid
  !-------------------------------------------------------------------
  ! This routine compute array rho (source term for Poisson equation)
  ! by first reseting array rho to zero, then 
  ! by affecting the gas density to leaf cells, and finally
  ! by performing a restriction operation for split cells.
  ! For pure particle runs, the restriction is not necessary and the
  ! routine only set rho to zero. On the other hand, for the Multigrid
  ! solver, the restriction is necessary in any case.
  !-------------------------------------------------------------------
  integer ::ind,i,j,icpu,ncache,ngrid,iskip,info,ibound,nx_loc
  integer ::idim,nleaf,ix,iy,iz,igrid
  integer,dimension(1:nvector),save::ind_grid

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Initialize density field to zero
  do icpu=1,ncpu
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,reception(icpu,ilevel)%ngrid
           rho(reception(icpu,ilevel)%igrid(i)+iskip)=0.0D0
        end do
     end do
  end do
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,active(ilevel)%ngrid
        rho(active(ilevel)%igrid(i)+iskip)=0.0D0
     end do
  end do
  ! Reset rho in physical boundaries
  do ibound=1,nboundary
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,boundary(ibound,ilevel)%ngrid
           rho(boundary(ibound,ilevel)%igrid(i)+iskip)=0.0
        end do
     end do
  end do
  
  if(hydro)then
     ! Perform a restriction over split cells (ilevel+1)
     ncache=active(ilevel)%ngrid
     do igrid=1,ncache,nvector
        ! Gather nvector grids
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
        end do
        call cic_cell(ind_grid,ngrid,ilevel)
     end do
  end if

111 format('   Entering cic_from_multipole for level',i2)

end subroutine cic_from_multipole
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cic_cell(ind_grid,ngrid,ilevel)
  use amr_commons
  use poisson_commons
  use hydro_commons, ONLY: unew
  implicit none
  integer::ngrid,ilevel
  integer,dimension(1:nvector)::ind_grid
  !
  !
  integer::i,j,idim,ind_cell_son,iskip_son,np,ind_son,nx_loc,ind
  integer ,dimension(1:nvector),save::ind_cell,ind_cell_father
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  ! Particle-based arrays
  logical ,dimension(1:nvector),save::ok
  real(dp),dimension(1:nvector),save::mmm,ttt
  real(dp),dimension(1:nvector),save::vol2
  real(dp),dimension(1:nvector,1:ndim),save::x,dd,dg
  integer ,dimension(1:nvector,1:ndim),save::ig,id,igg,igd,icg,icd
  real(dp),dimension(1:nvector,1:twotondim),save::vol
  integer ,dimension(1:nvector,1:twotondim),save::igrid,icell,indp,kg
  real(dp),dimension(1:3)::skip_loc
  real(kind=8)::dx,dx_loc,scale,vol_loc
  logical::error
  
  ! Mesh spacing in that level
  dx=0.5D0**ilevel 
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc=dx_loc**ndim
  np=ngrid

  ! Compute father cell index
  do i=1,ngrid
     ind_cell(i)=father(ind_grid(i))
  end do

  ! Gather 3x3x3 neighboring parent cells
  call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids,ngrid,ilevel)

  ! Loop over grid cells
  do ind_son=1,twotondim
     iskip_son=ncoarse+(ind_son-1)*ngridmax

     ! Compute pseudo particle (centre of mass) position
     do idim=1,ndim
        do j=1,np
           ind_cell_son=iskip_son+ind_grid(j)
           x(j,idim)=unew(ind_cell_son,idim+1)/unew(ind_cell_son,1)
        end do
     end do
     
     ! Compute total multipole
     if(ilevel==levelmin)then
        do idim=1,ndim+1
           do j=1,np
              ind_cell_son=iskip_son+ind_grid(j)
              multipole(idim)=multipole(idim)+unew(ind_cell_son,idim)
           end do
        end do
     endif

     ! Rescale particle position at level ilevel
     do idim=1,ndim
        do j=1,np
           x(j,idim)=x(j,idim)/scale+skip_loc(idim)
        end do
     end do
     do idim=1,ndim
        do j=1,np
           x(j,idim)=x(j,idim)-(xg(ind_grid(j),idim)-3d0*dx)
        end do
     end do
     do idim=1,ndim
        do j=1,np
           x(j,idim)=x(j,idim)/dx
        end do
     end do
     
     ! Gather particle mass
     do j=1,np
        ind_cell_son=iskip_son+ind_grid(j)
        mmm(j)=unew(ind_cell_son,1)
     end do
     
     ! CIC at level ilevel (dd: right cloud boundary; dg: left cloud boundary)
     do idim=1,ndim
        do j=1,np
           dd(j,idim)=x(j,idim)+0.5D0
           id(j,idim)=dd(j,idim)
           dd(j,idim)=dd(j,idim)-id(j,idim)
           dg(j,idim)=1.0D0-dd(j,idim)
           ig(j,idim)=id(j,idim)-1
        end do
     end do
     
     ! Check for illegal moves
     error=.false.
     do idim=1,ndim
        do j=1,np
           if(x(j,idim)<0.5D0.or.x(j,idim)>5.5D0)error=.true.
        end do
     end do
     if(error)then
        write(*,*)'problem in cic'
        do idim=1,ndim
           do j=1,np
              if(x(j,idim)<0.5D0.or.x(j,idim)>5.5D0)then
                 write(*,*)x(j,1:ndim)
              endif
           end do
        end do
        stop
     end if

     ! Compute cloud volumes
#if NDIM==1
     do j=1,np
        vol(j,1)=dg(j,1)
        vol(j,2)=dd(j,1)
     end do
#endif
#if NDIM==2
     do j=1,np
        vol(j,1)=dg(j,1)*dg(j,2)
        vol(j,2)=dd(j,1)*dg(j,2)
        vol(j,3)=dg(j,1)*dd(j,2)
        vol(j,4)=dd(j,1)*dd(j,2)
     end do
#endif
#if NDIM==3
     do j=1,np
        vol(j,1)=dg(j,1)*dg(j,2)*dg(j,3)
        vol(j,2)=dd(j,1)*dg(j,2)*dg(j,3)
        vol(j,3)=dg(j,1)*dd(j,2)*dg(j,3)
        vol(j,4)=dd(j,1)*dd(j,2)*dg(j,3)
        vol(j,5)=dg(j,1)*dg(j,2)*dd(j,3)
        vol(j,6)=dd(j,1)*dg(j,2)*dd(j,3)
        vol(j,7)=dg(j,1)*dd(j,2)*dd(j,3)
        vol(j,8)=dd(j,1)*dd(j,2)*dd(j,3)
     end do
#endif
     
     ! Compute parent grids
     do idim=1,ndim
        do j=1,np
           igg(j,idim)=ig(j,idim)/2
           igd(j,idim)=id(j,idim)/2
        end do
     end do
#if NDIM==1
     do j=1,np
        kg(j,1)=1+igg(j,1)
        kg(j,2)=1+igd(j,1)
     end do
#endif
#if NDIM==2
     do j=1,np
        kg(j,1)=1+igg(j,1)+3*igg(j,2)
        kg(j,2)=1+igd(j,1)+3*igg(j,2)
        kg(j,3)=1+igg(j,1)+3*igd(j,2)
        kg(j,4)=1+igd(j,1)+3*igd(j,2)
     end do
#endif
#if NDIM==3
     do j=1,np
        kg(j,1)=1+igg(j,1)+3*igg(j,2)+9*igg(j,3)
        kg(j,2)=1+igd(j,1)+3*igg(j,2)+9*igg(j,3)
        kg(j,3)=1+igg(j,1)+3*igd(j,2)+9*igg(j,3)
        kg(j,4)=1+igd(j,1)+3*igd(j,2)+9*igg(j,3)
        kg(j,5)=1+igg(j,1)+3*igg(j,2)+9*igd(j,3)
        kg(j,6)=1+igd(j,1)+3*igg(j,2)+9*igd(j,3)
        kg(j,7)=1+igg(j,1)+3*igd(j,2)+9*igd(j,3)
        kg(j,8)=1+igd(j,1)+3*igd(j,2)+9*igd(j,3)
     end do
#endif
     do ind=1,twotondim
        do j=1,np
           igrid(j,ind)=son(nbors_father_cells(j,kg(j,ind)))
        end do
     end do
     
     ! Compute parent cell position
     do idim=1,ndim
        do j=1,np
           icg(j,idim)=ig(j,idim)-2*igg(j,idim)
           icd(j,idim)=id(j,idim)-2*igd(j,idim)
        end do
     end do
#if NDIM==1
     do j=1,np
        icell(j,1)=1+icg(j,1)
        icell(j,2)=1+icd(j,1)
     end do
#endif
#if NDIM==2
     do j=1,np
        icell(j,1)=1+icg(j,1)+2*icg(j,2)
        icell(j,2)=1+icd(j,1)+2*icg(j,2)
        icell(j,3)=1+icg(j,1)+2*icd(j,2)
        icell(j,4)=1+icd(j,1)+2*icd(j,2)
     end do
#endif
#if NDIM==3
     do j=1,np
        icell(j,1)=1+icg(j,1)+2*icg(j,2)+4*icg(j,3)
        icell(j,2)=1+icd(j,1)+2*icg(j,2)+4*icg(j,3)
        icell(j,3)=1+icg(j,1)+2*icd(j,2)+4*icg(j,3)
        icell(j,4)=1+icd(j,1)+2*icd(j,2)+4*icg(j,3)
        icell(j,5)=1+icg(j,1)+2*icg(j,2)+4*icd(j,3)
        icell(j,6)=1+icd(j,1)+2*icg(j,2)+4*icd(j,3)
        icell(j,7)=1+icg(j,1)+2*icd(j,2)+4*icd(j,3)
        icell(j,8)=1+icd(j,1)+2*icd(j,2)+4*icd(j,3)
     end do
#endif
     
     ! Compute parent cell adress
     do ind=1,twotondim
        do j=1,np
           indp(j,ind)=ncoarse+(icell(j,ind)-1)*ngridmax+igrid(j,ind)
        end do
     end do
     
     ! Update mass density and number density fields
     do ind=1,twotondim
        do j=1,np
           ok(j)=igrid(j,ind)>0
        end do
        do j=1,np
           vol2(j)=mmm(j)*vol(j,ind)/vol_loc
        end do        
        do j=1,np
           if(ok(j))then
              rho(indp(j,ind))=rho(indp(j,ind))+vol2(j)
           end if
        end do
     end do
     
  end do
  ! End loop over grid cells

end subroutine cic_cell
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine mass_deposit(xpart, mpart, nparts, grid_level, nbits_patch)
  use amr_parameters, only: ndim, dp, int_pre
  use amr_commons,    only: ncpu, ind_table2, boxlen
  use pm_utils,       only: patched_particle_loop
  implicit none
  integer, intent(in) :: grid_level, nparts, nbits_patch
  real(dp), dimension(:, :), intent(inout) :: xpart
  real(dp), dimension(:), intent(in) :: mpart
  
  ! This routine deposits the input particles to the AMR grid at
  ! level grid_level. It uses a regular cartesian grid patch as
  ! a "3d-histogram" before accessing the hash table.

  integer :: patch_size
  real(dp), allocatable, dimension(:,:,:,:), target :: rho_tmp_data
  real(dp) :: dx

  
  ! Place a warning sign to make sure the current limitations on this routine are known.
  if (ncpu > 1)then
     print*, 'particle patch deposition is not yet implemented for MPI'
     stop
  end if

  if (nparts == 0) return  
  patch_size = 2 ** nbits_patch
  dx = boxlen * 0.5d0 ** grid_level
  
  ! Allocate two cell-thick boundaries to make the depostion onto the AMR grid
  ! simpler.
  allocate(rho_tmp_data(-2: patch_size + 1, -2: patch_size + 1, -2: patch_size + 1, 1:2))
  call patched_particle_loop(xpart, nparts, grid_level, 3, mass_deposit_callback)  
  deallocate(rho_tmp_data)

contains

  subroutine mass_deposit_callback(oft, np, grid_offset)
    use amr_parameters, only: nvector
    use cic_stuff,      only: dump_rho_tmp, cic_nvector
    implicit none
    integer, intent(in), value :: oft, np
    integer(int_pre), dimension(1: ndim) :: grid_offset

    integer(int_pre), dimension(1:nvector, 1:ndim, 0:7) :: ix
    real(dp),         dimension(1:nvector, 0:7)         :: vol
    real(dp), pointer                                   :: rho_tmp(:, :, :, :)  
    
    integer :: idim, ip, icell, sweep_offset, sweep_nparts

    ! Use a pointer to get the right indexing for the pre-allocated grid patch
    rho_tmp(grid_offset(1) - 2: grid_offset(1) + patch_size + 1, &
            grid_offset(2) - 2: grid_offset(2) + patch_size + 1, &
            grid_offset(3) - 2: grid_offset(3) + patch_size + 1, &
            1:2) => rho_tmp_data

    rho_tmp = 0.d0
    ! Loop particles in nvector sweeps
    do sweep_offset = 0, np - 1, nvector
       sweep_nparts = min(np - sweep_offset, nvector)
       
       ! Get cloud corner integer coordinates and cloud fractions
       call cic_nvector(xpart(oft + sweep_offset + 1: oft + sweep_offset + sweep_nparts, 1:ndim), ix, vol, sweep_nparts, dx)
       
       ! Add mass number density to temporary grid patches
       do icell = 0, 7
          do ip = 1, sweep_nparts             
             rho_tmp(ix(ip, 1, icell), ix(ip, 2, icell), ix(ip, 3, icell), 1) = &
                  rho_tmp(ix(ip, 1, icell), ix(ip, 2, icell), ix(ip, 3, icell), 1) + mpart(oft + sweep_offset + ip) * vol(ip, icell) 
             rho_tmp(ix(ip, 1, icell), ix(ip, 2, icell), ix(ip, 3, icell), 2) = &
                  rho_tmp(ix(ip, 1, icell), ix(ip, 2, icell), ix(ip, 3, icell), 2) + vol(ip, icell)
          end do
       end do
    end do
    rho_tmp(:,:,:,1) = rho_tmp(:,:,:,1) / (dx**3)
    
    ! Deposit density / number density from patches to AMR grids
    call dump_rho_tmp(rho_tmp, grid_offset, patch_size, grid_level)
  end subroutine mass_deposit_callback

end subroutine mass_deposit


!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine add_particle_multipole
  use amr_parameters, only: ndim
  use amr_commons, only: myid
  use pm_commons, only: xp, mp, npart
  use poisson_commons, only: multipole
  implicit none
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Simple routine to compute the multipole contribution from particles
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer :: ipart, idim

  do ipart = 1, npart
     multipole(1) = multipole(1) + mp(ipart)
  end do
  do idim = 1, ndim
     do ipart = 1, npart
        multipole(idim + 1) = multipole(idim + 1) + mp(ipart) * xp(ipart, idim)
     end do
  end do
       
end subroutine add_particle_multipole


!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine kick_part(xpart, vpart, levelp, nparts, grid_level, nbits_patch, previous_timestep)
  use amr_parameters, only: ndim, dp, int_pre, MASK_VALUE
  use amr_commons,    only: ncpu, ind_table2, boxlen, dtnew, dtold
  use cic_stuff,      only: load_f_tmp
  use pm_utils,       only: patched_particle_loop
  implicit none
  integer, intent(in) :: grid_level, nparts
  integer, value, intent(in) :: nbits_patch
  integer, dimension(:), intent(inout) :: levelp
  real(dp), dimension(:, :), intent(inout) :: xpart, vpart
  logical, intent(in), value :: previous_timestep
  
  ! This routine deposits the input particles to the AMR grid at
  ! level grid_level. It uses a regular cartesian grid patch as
  ! a "3d-histogram" before accessing the hash table.

  integer :: idim, patch_size, patch_size_coarse
  real(dp), allocatable, dimension(:,:,:,:), target :: f_tmp_data, f_tmp_coarse_data
  real(dp) :: dx


  if (nparts==0)return
  
  ! Place a warning sign to make sure the current limitations on this routine are known.
  if (ncpu > 1)then
     print*, 'particle patch deposition is not yet implemented for MPI'
     stop
  end if

  patch_size = 2 ** nbits_patch
  patch_size_coarse  = max(patch_size / 2, 2)
  dx = boxlen * 0.5d0 ** grid_level
  
  ! Allocate two cell-thick boundaries to make the depostion onto the AMR grid
  ! simpler.
  allocate(f_tmp_data(-2: patch_size + 1, -2: patch_size + 1, -2: patch_size + 1, 1:ndim))
  allocate(f_tmp_coarse_data(-2: patch_size_coarse + 1, -2: patch_size_coarse + 1, -2: patch_size_coarse + 1, 1:ndim))

  call patched_particle_loop(xpart, nparts, grid_level, 3, kick_part_callback)

  deallocate(f_tmp_data, f_tmp_coarse_data)

contains

  subroutine kick_part_callback(oft, np, grid_offset)
    use amr_parameters, only: nvector
    use cic_stuff,      only: dump_rho_tmp, cic_nvector
    implicit none
    integer, intent(in), value :: oft, np
    integer(int_pre), dimension(1: ndim) :: grid_offset

    integer(int_pre), dimension(1:nvector, 1:ndim, 0:7) :: ix
    real(dp),         dimension(1:nvector, 0:7)         :: vol
    real(dp),         dimension(1:nvector, 1:ndim)      :: ap
    real(dp),         dimension(1:nvector)              :: dteff
    logical,          dimension(1:nvector)              :: repeat_coarser
    integer(int_pre), dimension(1: ndim) :: grid_offset_coarse
    real(dp), pointer :: f_tmp(:, :, :, :), f_tmp_coarse(:, :, :, :)
    
    integer :: idim, ip, ipart, icell, sweep_offset, sweep_nparts
    logical :: all_ok

    ! Use pointers to get the right indexing for the pre-allocated grid patch
    f_tmp(grid_offset(1) - 2: grid_offset(1) + patch_size + 1, &
          grid_offset(2) - 2: grid_offset(2) + patch_size + 1, &
          grid_offset(3) - 2: grid_offset(3) + patch_size + 1, &
          1:ndim) => f_tmp_data

    grid_offset_coarse = grid_offset / 2
     
    f_tmp_coarse(grid_offset_coarse(1) - 2: grid_offset_coarse(1) + patch_size_coarse + 1, &
                 grid_offset_coarse(2) - 2: grid_offset_coarse(2) + patch_size_coarse + 1, &
                 grid_offset_coarse(3) - 2: grid_offset_coarse(3) + patch_size_coarse + 1, &
                 1:ndim) => f_tmp_coarse_data

    call load_f_tmp(f_tmp, grid_offset, patch_size, grid_level)
    call load_f_tmp(f_tmp_coarse, grid_offset_coarse, patch_size_coarse, grid_level - 1)

    
    ! Loop particles in nvector sweeps
    do sweep_offset = 0, np - 1, nvector
       sweep_nparts = min(np - sweep_offset, nvector)
       
       ! Get cloud corner integer coordinates and cloud fractions
       call cic_nvector(xpart(oft + sweep_offset + 1: oft + sweep_offset + sweep_nparts, 1:ndim), ix, vol, sweep_nparts, dx)
       

       repeat_coarser = .false.
       all_ok = .true.
       ap = 0.d0
       do icell = 0, 7
          do ip = 1, sweep_nparts
             if (f_tmp(ix(ip, 1, icell), ix(ip, 2, icell), ix(ip, 3, icell), 1) == MASK_VALUE) then
                repeat_coarser(ip) = .true.
                all_ok = .false.
             else
                ! Unfortunate array stride here!!!
                ap(ip, 1:ndim) =  ap(ip, 1:ndim) + vol(ip, icell) * f_tmp(ix(ip, 1, icell), ix(ip, 2, icell), ix(ip, 3, icell), 1:ndim)
             end if
          end do
       end do

       ! For particles which are partially in a coarser level, repeat at coarse level
       do ip = 1, sweep_nparts
          if (repeat_coarser(ip)) then
             ! Maybe write cic_one subroutine...
             call cic_nvector(xpart(oft + sweep_offset + ip: oft + sweep_offset + ip, 1:ndim), ix(1:1, 1:ndim, 0:7), vol(1:1, 0:7), 1, 2 * dx)
             ap(ip, 1:3) = 0.d0
             do icell = 0, 7
                ! Unfortunate array stride here!!!
                ap(ip, 1:ndim) =  ap(ip, 1:ndim) + vol(1, icell) * f_tmp_coarse(ix(1, 1, icell), ix(1, 2, icell), ix(1, 3, icell), 1:ndim)
             end do
          end if
       end do

       ! Compute individual time step
       if (previous_timestep)then
          do ip = 1, sweep_nparts
             ipart = oft + sweep_offset + ip 
             if(levelp(ipart) >= grid_level)then
                dteff(ip) = 0.5d0 * dtnew(levelp(ipart))
             else
                dteff(ip) = 0.5d0 * dtold(levelp(ipart))
             endif
          end do
       else
          dteff = 0.5d0 * dtnew(grid_level)
       end if

       ! Finally, apply the kick
       do ip = 1, sweep_nparts
          ipart = oft + sweep_offset + ip
          vpart(ipart, 1:ndim) = vpart(ipart, 1:ndim) + ap(ip, 1:ndim) * dteff(ip)
       end do

    end do
  end subroutine kick_part_callback
end subroutine kick_part
