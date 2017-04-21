!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine godunov_fine_vec_orig(ilevel)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel
  !--------------------------------------------------------------------------
  ! This routine is a wrapper to the second order Godunov solver.
  ! Small grids (2x2x2) are gathered from level ilevel and sent to the
  ! hydro solver. On entry, hydro variables are gathered from array uold.
  ! On exit, unew has been updated. 
  !--------------------------------------------------------------------------
  integer::i,igrid,ngrid
  integer,dimension(1:nvector),save::ind_grid

  if(noct_tot(ilevel)==0)return
  if(static)return
  if(verbose)write(*,111)ilevel

  call open_cache(operation_godunov,domain_decompos_amr)

  ! Loop over active grids by vector sweeps
  do igrid=head(ilevel),tail(ilevel),nvector
     ngrid=MIN(nvector,tail(ilevel)-igrid+1)
     do i=1,ngrid
        ind_grid(i)=igrid+i-1
     end do
     call godfine1_vec_origin(ind_grid,ngrid,ilevel)
  end do

  call close_cache(grid_dict)

111 format('   Entering godunov_fine for level ',i2)

end subroutine godunov_fine_vec_orig
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine godfine1_vec_origin(ind_grid,ngrid,ilevel)
  use amr_commons
  use hash
  use hydro_vec_parameters
  use hydro_parameters
  use poisson_commons
  implicit none
  integer::ilevel,ngrid
  integer,dimension(1:nvector)::ind_grid
  !-------------------------------------------------------------------
  ! This routine gathers first hydro variables from neighboring grids
  ! to set initial conditions in a 6x6x6 grid. It interpolate from
  ! coarser level missing grid variables. It then calls the
  ! Godunov solver that computes fluxes. These fluxes are zeroed at 
  ! coarse-fine boundaries, since contribution from finer levels has
  ! already been taken into account. Conservative variables are updated 
  ! and stored in array unew(:), both at the current level and at the 
  ! coarser level if necessary.
  !-------------------------------------------------------------------
  integer ,dimension(1:nvector,1:threetondim,1:ndim),save::cart_nbor
  integer ,dimension(1:nvector,1:threetondim),save::igrid_nbor
  integer ,dimension(1:nvector,1:threetondim),save::igrid_father_nbor
  integer ,dimension(1:nvector,1:threetondim),save::ind_father_nbor
  integer ,dimension(1:nvector,1:threetondim,0:twondim),save::igrid_interpol_nbor
  integer ,dimension(1:nvector,0:twondim),save::igrid_nbor2
  integer ,dimension(1:nvector,0:twondim),save::ind_nbor2

  real(dp),dimension(1:nvector,0:twondim,1:nvar),save::u1
  real(dp),dimension(1:nvector,1:twotondim,1:nvar),save::u2

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar),save::uloc
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndim),save::gloc=0.0d0
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:nvar,1:ndim),save::flux
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:2,1:ndim),save::tmp
  logical ,dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),save::ok

  integer(kind=8),dimension(1:nvector,0:ndim),save::hash_key
  integer,dimension(1:nvector),save::igrid_father,ind_father,ind_exist,igrid_exist,ind_nexist
  integer,dimension(1:nvector),save::grid_address

  integer::i,ind,ivar,idim,ioct,ind_son,inbor,ind_nbor,nbuffer
  integer::i0,j0,k0,ii0,jj0,kk0,i1,j1,k1,i2,j2,k2,i3,j3,k3,nb_noneigh,nexist
  integer::i1min,i1max,j1min,j1max,k1min,k1max
  integer::i2min,i2max,j2min,j2max,k2min,k2max
  integer::i3min,i3max,j3min,j3max,k3min,k3max
  real(dp)::dx,oneontwotondim

  oneontwotondim = 1.d0/dble(twotondim)

  ! Mesh spacing in that level
  dx=boxlen/2**ilevel

  ! Integer constants
  i1min=0; i1max=0; i2min=0; i2max=0; i3min=1; i3max=1
  j1min=0; j1max=0; j2min=0; j2max=0; j3min=1; j3max=1
  k1min=0; k1max=0; k2min=0; k2max=0; k3min=1; k3max=1
  if(ndim>0)then
     i1max=2; i2max=1; i3max=2
  end if
  if(ndim>1)then
     j1max=2; j2max=1; j3max=2
  end if
  if(ndim>2)then
     k1max=2; k2max=1; k3max=2
  end if

  !-------------------------------------------------
  ! Compute 3x3x3 neighboring grids Cartesian keys
  !-------------------------------------------------
  do k1=k1min,k1max
     do j1=j1min,j1max
        do i1=i1min,i1max     
           ind_nbor=1+i1+3*j1+9*k1
           do i=1,ngrid
#if NDIM>0
              cart_nbor(i,ind_nbor,1)=grid(ind_grid(i))%ckey(1)+i1-1
#endif
#if NDIM>1
              cart_nbor(i,ind_nbor,2)=grid(ind_grid(i))%ckey(2)+j1-1
#endif
#if NDIM>2
              cart_nbor(i,ind_nbor,3)=grid(ind_grid(i))%ckey(3)+k1-1
#endif
           end do
           ! Periodic boundary conditons
           do idim=1,ndim
              do i=1,ngrid
                 if(cart_nbor(i,ind_nbor,idim)<0)cart_nbor(i,ind_nbor,idim)=ckey_max(ilevel)-1
                 if(cart_nbor(i,ind_nbor,idim)==ckey_max(ilevel))cart_nbor(i,ind_nbor,idim)=0
              enddo
           end do
        end do
     end do
  end do
  
  !------------------------------------------
  ! Set hash key level
  !------------------------------------------
  do i=1,ngrid
     hash_key(i,0)=ilevel
  end do

  !------------------------------------------
  ! Compute 3x3x3 neighboring grid index
  !------------------------------------------
  do ind=1,threetondim
     do idim=1,ndim
        do i=1,ngrid
           hash_key(i,idim)=cart_nbor(i,ind,idim)
        end do
     end do
     call get_grid_vec(hash_key,grid_dict,grid_address,.false.,.true.,ngrid)
     ! Store neighboring octs address for subsequent use.
     do i=1,ngrid
        igrid_nbor(i,ind)=grid_address(i)
     end do
  end do

  !-----------------------------------------
  ! Reset neighoring father grids to zero
  !-----------------------------------------
  igrid_father_nbor=0
  igrid_interpol_nbor=0

  !---------------------------
  ! Gather 6x6x6 cells stencil
  !---------------------------
  ! Loop over 3x3x3 neighboring grid
  do k1=k1min,k1max
     do j1=j1min,j1max
        do i1=i1min,i1max
           
           ind_nbor=1+i1+3*j1+9*k1
           
           ! Check if neighboring grid exists or not
           nexist=0
           nbuffer=0
           do i=1,ngrid
              if(igrid_nbor(i,ind_nbor)>0) then
                 nexist=nexist+1
                 ind_exist(nexist)=i
                 igrid_exist(nexist)=igrid_nbor(i,ind_nbor)
              else
                 nbuffer=nbuffer+1
                 ind_nexist(nbuffer)=i
                 do idim=1,ndim
                    hash_key(nbuffer,idim)=cart_nbor(i,ind_nbor,idim)
                 end do
              end if
           end do

           ! Collect father cells for non-existing grids
           if(nbuffer>0)then
              call get_parent_cell_vec(hash_key,grid_dict,igrid_father,ind_father,.true.,.true.,nbuffer)           
              ! Store father cells address for subsequent use.
              do i=1,nbuffer
                 igrid_father_nbor(ind_nexist(i),ind_nbor)=igrid_father(i)
                 ind_father_nbor(ind_nexist(i),ind_nbor)=ind_father(i)
              end do
           endif
           
           ! In case one wants to interpolate using high-order schemes
           if(nbuffer>0.and.interpol_type>0)then
              
              ! Get 2ndim neighboring father cells with read-write cache
              call get_twondim_nbor_parent_cell_vec(hash_key,grid_dict,igrid_nbor2,ind_nbor2,.true.,.true.,nbuffer)
              ! Store interpolation octs address for subsequent use.
              do inbor=0,twotondim
                 do i=1,nbuffer
                    igrid_interpol_nbor(ind_nexist(i),ind_nbor,inbor)=igrid_nbor2(i,inbor)
                 end do
              end do
              
              ! Fill interpolation work space
              do inbor=0,twondim
                 do ivar=1,nvar
                    do i=1,nbuffer
                       u1(i,inbor,ivar)=grid(igrid_nbor2(i,inbor))%uold(ind_nbor2(i,inbor),ivar)
                    end do
                 end do
              end do
              
              ! Interpolate
              call interpol_hydro_vec(u1,u2,nbuffer)
              
           endif
           
           ! Loop over 2x2x2 cells
           do k2=k2min,k2max
              do j2=j2min,j2max
                 do i2=i2min,i2max
                    
                    ind_son=1+i2+2*j2+4*k2
                    
                    i3=1; j3=1; k3=1
                    if(ndim>0)i3=1+2*(i1-1)+i2
                    if(ndim>1)j3=1+2*(j1-1)+j2
                    if(ndim>2)k3=1+2*(k1-1)+k2
                    
                    ! Gather hydro variables
                    do ivar=1,nvar
                       do i=1,nexist
                          uloc(ind_exist(i),i3,j3,k3,ivar)=grid(igrid_exist(i))%uold(ind_son,ivar)
                       end do
                       if(nbuffer>0.and.interpol_type>0)then
                          do i=1,nbuffer
                             uloc(ind_nexist(i),i3,j3,k3,ivar)=u2(i,ind_son,ivar)
                          end do
                       else if(nbuffer>0)then
                          do i=1,nbuffer
                             uloc(ind_nexist(i),i3,j3,k3,ivar)=grid(igrid_father(i))%uold(ind_father(i),ivar)
                          end do
                       endif
                    end do             
#ifdef GRAV
                    ! Gather gravitational acceleration
                    if(poisson)then
                       do idim=1,ndim
                          do i=1,nexist
                             gloc(ind_exist(i),i3,j3,k3,idim)=grid(igrid_exist(i))%f(ind_son,idim)
                          end do
                          ! Use straight injection for buffer cells
                          do i=1,nbuffer
                             gloc(ind_nexist(i),i3,j3,k3,idim)=grid(igrid_father(i))%f(ind_father(i),idim)
                          end do
                       end do
                    end if
#endif                    
                    ! Gather refinement flag
                    do i=1,nexist
                       ok(ind_exist(i),i3,j3,k3)=grid(igrid_exist(i))%refined(ind_son)
                    end do
                    do i=1,nbuffer
                       ok(ind_nexist(i),i3,j3,k3)=.false.
                    end do
                    
                 end do
              end do
           end do
           ! End loop over cells
           
        end do
     end do
  end do
  ! End loop over neighboring grids
  
  !-----------------------------------------------
  ! Compute flux using second-order Godunov method
  !-----------------------------------------------
  call unsplit_vec(uloc,gloc,flux,tmp,dx,dx,dx,dtnew(ilevel),ngrid)

  !------------------------------------------------
  ! Reset flux along direction at refined interface    
  !------------------------------------------------
  do idim=1,ndim
     i0=0; j0=0; k0=0
     if(idim==1)i0=1
     if(idim==2)j0=1
     if(idim==3)k0=1
     do k3=k3min,k3max+k0
        do j3=j3min,j3max+j0
           do i3=i3min,i3max+i0
              do ivar=1,nvar
                 do i=1,ngrid
                    if(ok(i,i3-i0,j3-j0,k3-k0) .or. ok(i,i3,j3,k3))then
                       flux(i,i3,j3,k3,ivar,idim)=0.0d0
                    end if
                 end do
              end do
              if(pressure_fix)then
                 do ivar=1,2
                    do i=1,ngrid
                       if(ok(i,i3-i0,j3-j0,k3-k0) .or. ok(i,i3,j3,k3))then
                          tmp(i,i3,j3,k3,ivar,idim)=0.0d0
                       end if
                    end do
                 end do
              end if
           end do
        end do
     end do
  end do
  !--------------------------------------
  ! Conservative update at level ilevel
  !--------------------------------------
  do idim=1,ndim
     i0=0; j0=0; k0=0
     if(idim==1)i0=1
     if(idim==2)j0=1
     if(idim==3)k0=1
     do k2=k2min,k2max
        do j2=j2min,j2max
           do i2=i2min,i2max
              ind_son=1+i2+2*j2+4*k2
              i3=1+i2
              j3=1+j2
              k3=1+k2
              ! Update conservative variables new state vector
              do ivar=1,nvar
                 do i=1,ngrid
                    grid(ind_grid(i))%unew(ind_son,ivar)=grid(ind_grid(i))%unew(ind_son,ivar)+ &
                         & (flux(i,i3   ,j3   ,k3   ,ivar,idim) &
                         & -flux(i,i3+i0,j3+j0,k3+k0,ivar,idim))
                 end do
              end do
#ifdef DUALENER
              ! Update velocity divergence
              do i=1,ngrid
                 grid(ind_grid(i))%divu(ind_son)=grid(ind_grid(i))%divu(ind_son)+ &
                      & (tmp(i,i3   ,j3   ,k3   ,1,idim) &
                      & -tmp(i,i3+i0,j3+j0,k3+k0,1,idim))
              end do
              ! Update internal energy
              do i=1,ngrid
                 grid(ind_grid(i))%enew(ind_son)=grid(ind_grid(i))%enew(ind_son)+ &
                      & (tmp(i,i3   ,j3   ,k3   ,2,idim) &
                      & -tmp(i,i3+i0,j3+j0,k3+k0,2,idim))
              end do
#endif
           end do
        end do
     end do
  end do
  
  !--------------------------------------
  ! Conservative update at level ilevel-1
  !--------------------------------------
  ! Loop over dimensions
  do idim=1,ndim
     i0=0; j0=0; k0=0
     if(idim==1)i0=1
     if(idim==2)j0=1
     if(idim==3)k0=1
     if(ndim>0)ii0=1
     if(ndim>1)jj0=1
     if(ndim>2)kk0=1
     i1=i1min+ii0
     j1=j1min+jj0
     k1=k1min+kk0
     !----------------------
     ! Left flux at boundary
     !----------------------     
     if(idim==1)i1=i1min
     if(idim==2)j1=j1min
     if(idim==3)k1=k1min
     ! Check if grids sits near left boundary
     ! and gather neighbor father cells index
     ind_nbor=1+i1+3*j1+9*k1
     nb_noneigh=0
     do i=1,ngrid
        if (igrid_nbor(i,ind_nbor)==0) then
           nb_noneigh = nb_noneigh + 1
           igrid_father(nb_noneigh) = igrid_father_nbor(i,ind_nbor)
           ind_father(nb_noneigh) = ind_father_nbor(i,ind_nbor)
           ind_nexist(nb_noneigh) = i
        end if
     end do
     ! Conservative update of new state variables
     do ivar=1,nvar
        ! Loop over boundary cells
        do k3=k3min,k3max-k0
           do j3=j3min,j3max-j0
              do i3=i3min,i3max-i0
                 do i=1,nb_noneigh
                    grid(igrid_father(i))%unew(ind_father(i),ivar)=grid(igrid_father(i))%unew(ind_father(i),ivar) &
                         & -flux(ind_nexist(i),i3,j3,k3,ivar,idim)*oneontwotondim
                 end do
              end do
           end do
        end do
     end do
#ifdef DUALENER
     ! Update velocity divergence
     do k3=k3min,k3max-k0
        do j3=j3min,j3max-j0
           do i3=i3min,i3max-i0
              do i=1,nb_noneigh
                 grid(igrid_father(i))%divu(ind_father(i))=grid(igrid_father(i))%divu(ind_father(i)) &
                      & -tmp(ind_nexist(i),i3,j3,k3,1,idim)*oneontwotondim
              end do
           end do
        end do
     end do
     ! Update internal energy
     do k3=k3min,k3max-k0
        do j3=j3min,j3max-j0
           do i3=i3min,i3max-i0
              do i=1,nb_noneigh
                 grid(igrid_father(i))%enew(ind_father(i))=grid(igrid_father(i))%enew(ind_father(i)) &
                      & -tmp(ind_nexist(i),i3,j3,k3,2,idim)*oneontwotondim
              end do
           end do
        end do
     end do
#endif
     
     !-----------------------
     ! Right flux at boundary
     !-----------------------     
     if(idim==1)i1=i1max
     if(idim==2)j1=j1max
     if(idim==3)k1=k1max
     ! Check if grids sits near right boundary
     ! and gather neighbor father cells index
     ind_nbor=1+i1+3*j1+9*k1
     nb_noneigh=0
     do i=1,ngrid
        if (igrid_nbor(i,ind_nbor)==0) then
           nb_noneigh = nb_noneigh + 1
           igrid_father(nb_noneigh) = igrid_father_nbor(i,ind_nbor)
           ind_father(nb_noneigh) = ind_father_nbor(i,ind_nbor)
           ind_nexist(nb_noneigh) = i
        end if
     end do
     ! Conservative update of new state variables
     do ivar=1,nvar
        ! Loop over boundary cells
        do k3=k3min+k0,k3max
           do j3=j3min+j0,j3max
              do i3=i3min+i0,i3max
                 do i=1,nb_noneigh
                    grid(igrid_father(i))%unew(ind_father(i),ivar)=grid(igrid_father(i))%unew(ind_father(i),ivar) &
                         & +flux(ind_nexist(i),i3+i0,j3+j0,k3+k0,ivar,idim)*oneontwotondim
                 end do
              end do
           end do
        end do
     end do
#ifdef DUALENER
     ! Update velocity divergence
     do k3=k3min+k0,k3max
        do j3=j3min+j0,j3max
           do i3=i3min+i0,i3max
              do i=1,nb_noneigh
                 grid(igrid_father(i))%divu(ind_father(i))=grid(igrid_father(i))%divu(ind_father(i)) &
                      & +tmp(ind_nexist(i),i3+i0,j3+j0,k3+k0,1,idim)*oneontwotondim
              end do
           end do
        end do
     end do
     ! Update internal energy
     do k3=k3min+k0,k3max
        do j3=j3min+j0,j3max
           do i3=i3min+i0,i3max
              do i=1,nb_noneigh
                 grid(igrid_father(i))%enew(ind_father(i))=grid(igrid_father(i))%enew(ind_father(i)) &
                      & +tmp(ind_nexist(i),i3+i0,j3+j0,k3+k0,2,idim)*oneontwotondim
              end do
           end do
        end do
     end do
#endif

  end do
  ! End loop over dimensions

  !--------------------------------------
  ! Unlock all octs
  !--------------------------------------
  do ind=1,threetondim
     do i=1,ngrid
        ioct=igrid_nbor(i,ind)
        if(ioct>ngridmax)locked(ioct-ngridmax)=.false.
        ioct=igrid_father_nbor(i,ind)
        if(ioct>ngridmax)locked(ioct-ngridmax)=.false.
     end do
     if(interpol_type>0)then
        do inbor=0,twotondim
           do i=1,ngrid
              ioct=igrid_interpol_nbor(i,ind,inbor)
              if(ioct>ngridmax)locked(ioct-ngridmax)=.false.
           end do
        end do
     end if
  end do

end subroutine godfine1_vec_origin
