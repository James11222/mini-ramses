!################################################################
!################################################################
!################################################################
!################################################################
! subroutine memory_sort_level(ind_com,np,ilevel,icpu)
!   use pm_commons
!   use amr_commons

!   ! sort all the particles in memory according to their level

!   implicit none
!   integer::np,icpu,ilevel
!   integer,dimension(1:nvector)::ind_com
  
!   integer::i,idim,igrid
!   integer,dimension(1:nvector),save::ind_list,ind_part
!   logical,dimension(1:nvector),save::ok=.true.
!   integer::current_property

!   ! Compute parent grid index
!   do i=1,np
!      igrid=emission(icpu,ilevel)%fp(ind_com(i),1)
!      ind_list(i)=emission(icpu,ilevel)%igrid(igrid)
!   end do

!   ! Add particle to parent linked list
!   call remove_free(ind_part,np)
!   call add_list(ind_part,ind_list,ok,np)

!   ! Scatter particle level and identity
!   do i=1,np
!      levelp(ind_part(i))=emission(icpu,ilevel)%fp(ind_com(i),2)
!      idp   (ind_part(i))=emission(icpu,ilevel)%fp(ind_com(i),3)
!   end do

!   ! Scatter particle position and velocity
!   do idim=1,ndim
!   do i=1,np
!      xp(ind_part(i),idim)=emission(icpu,ilevel)%up(ind_com(i),idim     )
!      vp(ind_part(i),idim)=emission(icpu,ilevel)%up(ind_com(i),idim+ndim)
!   end do
!   end do

!   current_property = twondim+1

!   ! Scatter particle mass
!   do i=1,np
!      mp(ind_part(i))=emission(icpu,ilevel)%up(ind_com(i),current_property)
!   end do
!   current_property = current_property+1

! #ifdef OUTPUT_PARTICLE_POTENTIAL
!   ! Scatter particle phi
!   do i=1,np
!      ptcl_phi(ind_part(i))=emission(icpu,ilevel)%up(ind_com(i),current_property)
!   end do
!   current_property = current_property+1
! #endif

! end subroutine memory_sort_level
! !################################################################
! !################################################################
! !################################################################
! !################################################################
!################################################################
!################################################################
!################################################################
!################################################################
! subroutine part_to_cell_i8(part_array,sortind)
!   use pm_parameters, only: npartmax
!   implicit none
! #ifndef WITHOUTMPI
!   include 'mpif.h'
! #endif
  
!   integer(kind=8),dimension(1:npartmax)::part_array, sortind
  
!   ! Communication routine that sends a particle-based quantity
!   ! to the MPI domain that owns the respective cell
  
  
  
  

! #ifndef WITHOUTMPI
!   do j=1,part_recv_tot
!      ipart=part_recv_buf(j)-ipart_start(myid)
!      int_part_recv_buf(j)=xx(ipart)
!   end do
!   call MPI_ALLTOALLV(int8_part_recv_buf,part_recv_cnt,part_recv_oft,MPI_INTEGER, &
!        &             int8_part_send_buf,part_send_cnt,part_send_oft,MPI_INTEGER,MPI_COMM_WORLD,info)
  
  
!   deallocate(int8_part_send_buf,int8_part_recv_buf)
! #endif
  
! end subroutine part_to_cell_i8


! #########################################################################
! #########################################################################
! #########################################################################
! #########################################################################
subroutine write_ascii_parts
  use pm_commons
  use amr_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  !count particles per level for old particle layout
  integer::ilevel, ilun
  integer::igrid,jgrid,i,ngrid,ncache, ipart, jpart, next_part
  integer::ig,ip,npart1,npart2,npart2_tot,icpu,info
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nlevelmax)::npts
  character(len=80)::filename
  character(LEN=5)::nchar
  
  ilun=myid+100
  
  call title(ilun,nchar)
  filename='OLDParts'//nchar
 
  open(unit=ilun,file=filename,form='formatted')
  
  do ilevel=levelmin,levelmin
     ! Loop over cpus
     do icpu=1,ncpu
        igrid=headl(icpu,ilevel)
        ! Loop over grids
        do jgrid=1,numbl(icpu,ilevel)
           npart1=numbp(igrid)  ! Number of particles in the grid
           if(npart1>0)then
              ipart=headp(igrid)
              ! Loop over particles
              do jpart=1,npart1
                 ! Save next particle   <--- Very important !!!
                 next_part=nextp(ipart)
                 write(ilun,"(6(F20.15),2(I10))")xp(ipart,1:3),vp(ipart,1:3),idp(ipart), levelp(ipart) 
                 ipart=next_part  ! Go to next particle
              end do
           endif
           igrid=next(igrid)   ! Go to next grid                                                              
        end do
     end do
  end do
  close(ilun)

  filename='NEWParts'//nchar
  open(unit=ilun,file=filename,form='formatted')
  do ipart=1,npart
     write(ilun,"(6(F20.15),2(I10))")xp(ipart,1:3),vp(ipart,1:3), &
          idp(ipart), levelp(ipart)
  end do

  close(ilun)
end subroutine write_ascii_parts
! #########################################################################
! #########################################################################
! #########################################################################
! #########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine count_parts
  use pm_commons
  use amr_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  ! ugly routine to count all particles in the simulation level by level
  ! used for debugging
  integer::ilevel
  integer::igrid,jgrid,i,ngrid,ncache
  integer::ig,ip,npart1,npart2,npart2_tot,icpu,info
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nlevelmax)::npts

  npts=0  
  do ilevel=1,nlevelmax
     npart2=0
     ! Loop over cpus
     do icpu=1,ncpu
        igrid=headl(icpu,ilevel)
        ig=0
        ip=0
        ! Loop over grids
        do jgrid=1,numbl(icpu,ilevel)
           npart1=numbp(igrid)  ! Number of particles in the grid
           npart2=npart2+npart1
           igrid=next(igrid)   ! Go to next grid                                                              
        
        end do
     end do
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(npart2,npart2_tot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#else
     npart2_tot=npart2
#endif          
     npts(ilevel)=npart2_tot
  end do
  if (myid==1)print*,'total'
  if (myid==1)print*,npts(1:nlevelmax)
  
  

   npts=0  
   do ilevel=1,nlevelmax
      npart2=0        
      
      ncache=active(ilevel)%ngrid
      do igrid=1,ncache,nvector
         ngrid=MIN(nvector,ncache-igrid+1)
         do i=1,ngrid
            ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
         end do
         do i=1,ngrid
            npart2=npart2+numbp(ind_grid(i))
         end do
        
      end do

#ifndef WITHOUTMPI
      call MPI_ALLREDUCE(npart2,npart2_tot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#else
      npart2_tot=npart2
#endif          
      npts(ilevel)=npart2_tot
   end do
   if (myid==1)print*,'active'
   if (myid==1)print*,npts(1:nlevelmax)
   
   
end subroutine count_parts
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
! subroutine hilbert_allparts(ilevel)
!   use pm_commons
!   use amr_commons
!   use sort, only: quick_sort_keys, qsort_parts_in_mem
!   implicit none
! #ifndef WITHOUTMPI
!   include 'mpif.h'
! #endif


!   integer::ilevel
!   integer::igrid,jgrid,i,ngrid,ncache,ipart,jpart
!   integer::ip,npart1,icpu,info
!   integer,dimension(1:nvector)::ind_grid
!   integer,dimension(1:nlevelmax)::npts
! !  real(qdp),dimension(1:nvector)::order
!   real(dp),dimension(1:nvector,1:ndim)::xtest
!   integer(kind=8),dimension(1:nvector)::hkey0,hkey1,hkey2
!   integer::ntot
!   integer,allocatable,dimension(:)::order


!   ntot=0
!   ! Loop over cpus
!   do icpu=1,ncpu
!      igrid=headl(icpu,ilevel)
!      ip=0
!      ! Loop over grids
!      do jgrid=1,numbl(icpu,ilevel)
!         npart1=numbp(igrid)  ! Number of particles in the grid
!         if(npart1>0)then
!            ipart=headp(igrid)
!            ! Loop over particles
!            do jpart=1,npart1
!               ! Save next particle   <--- Very important !!!
!               ip=ip+1
!               xtest(ip,1:ndim)=xp(ipart,1:ndim)
!               if(ip==nvector)then
!                  call cmp_ordering_int(xtest,hkey2,hkey1,hkey0,ip)
!                  part_hkey(ntot+1:ntot+ip,2)=hkey2(1:ip)
!                  part_hkey(ntot+1:ntot+ip,1)=hkey1(1:ip)
!                  part_hkey(ntot+1:ntot+ip,0)=hkey0(1:ip)
!                  ntot=ntot+ip
!                  ip=0
!               end if
!               ipart=nextp(ipart)  ! Go to next particle
!            end do
!         endif
!         igrid=next(igrid)   ! Go to next grid
!      end do
!      if(ip>0)then 
!         call cmp_ordering_int(xtest,hkey2,hkey1,hkey0,ip)
!         part_hkey(ntot+1:ntot+ip,2)=hkey2(1:ip)
!         part_hkey(ntot+1:ntot+ip,1)=hkey1(1:ip)
!         part_hkey(ntot+1:ntot+ip,0)=hkey0(1:ip)
!         ntot=ntot+ip
!      end if
!   end do

!   print*,'ntot:',ntot,npart,myid

! !  call qsort_parts_in_mem(ntot,1)

!   allocate(order(1:ntot))
!   call quick_sort_keys(order, ntot)
  
! !  do i=1,ntot
! !     write(*,'(A8,I5,I10,3(I20))'),"myid: ",myid,i,order(i),part_hkey(i,1),part_hkey(i,0)
! !  end do
!   deallocate(order)

! end subroutine hilbert_allparts



! subroutine hilbert_allparts
!   use pm_commons
!   use amr_commons
!   use sort, only: qsort_parts_in_mem
!   implicit none
! #ifndef WITHOUTMPI
!   include 'mpif.h'
! #endif


!   integer::ilevel
!   integer::igrid,jgrid,i,ngrid,ncache,jpart
!   integer::ip,npart1,icpu,info
!   integer,dimension(1:nvector)::ind_grid
!   integer,dimension(1:nlevelmax)::npts
! !  real(qdp),dimension(1:nvector)::order
!   real(dp),dimension(1:nvector,1:ndim)::xtest
!   integer(kind=8),dimension(1:nvector)::hkey0,hkey1,hkey2
!   integer::ntot
!   integer,allocatable,dimension(:)::order


!   ntot=0
!   ip=0
!   do while(ntot+ip<npart)
!      ip=ip+1
!      xtest(ip,1:ndim)=xp(ntot+ip,1:ndim)
!      if(ip==nvector)then
!         call cmp_ordering_int(xtest,hkey2,hkey1,hkey0,ip)
!         part_hkey(ntot+1:ntot+ip,2)=hkey2(1:ip)
!         part_hkey(ntot+1:ntot+ip,1)=hkey1(1:ip)
!         part_hkey(ntot+1:ntot+ip,0)=hkey0(1:ip)
!         ntot=ntot+ip
!         ip=0
!      end if
!   end do
!   if(ip>0)then 
!      call cmp_ordering_int(xtest,hkey2,hkey1,hkey0,ip)
!      part_hkey(ntot+1:ntot+ip,2)=hkey2(1:ip)
!      part_hkey(ntot+1:ntot+ip,1)=hkey1(1:ip)
!      part_hkey(ntot+1:ntot+ip,0)=hkey0(1:ip)
!      ntot=ntot+ip
!   end if

!   print*,'ntot:',ntot,npart,myid

! !  call qsort_parts_in_mem(ntot,1)

! !   allocate(order(1:ntot))
! !   call quick_sort_keys(order, ntot)
  
! ! !  do i=1,ntot
! ! !     write(*,'(A8,I5,I10,3(I20))'),"myid: ",myid,i,order(i),part_hkey(i,1),part_hkey(i,0)
! ! !  end do
! !   deallocate(order)

! end subroutine hilbert_allparts
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
!   subroutine cic_amr(xpart, mpart, np, grid_level)
!     use amr_parameters,  only: static, mass_cut_refine
!     use amr_commons,     only: boxlen, icoarse_max, & 
!                                icoarse_min, nvector, ndim, nstep_coarse
!     use poisson_commons, only: multipole, rho, phi
!     use hilbert,         only: hilbert3d
!     implicit none
!     integer,  intent(in)                               :: np, grid_level
!     real(dp), intent(in), dimension(1:nvector)         :: mpart
!     real(dp), intent(in), dimension(1:nvector, 1:ndim) :: xpart
!     ! This routine deposits nvector particles (local or remote) onto the grid (local)
!     ! at level grid_level.

!     ! in:           - particle masses
!     !               - particle positions
!     !               - number of particles
!     !               - grid_level 
    
!     ! out:          - "corrupted" particle positions -> do not reuse xpart outside of 
!     !                 this subroutine  
    
!     ! side effect:  - updates rho field on level grid_level

!     integer(kind=8), dimension(1:nvector, 0:2),    save :: cloud_hkey
! !    integer(kind=8), dimension(1:nvector, 1:ndim), save :: id
!     integer(kind=8), dimension(1:nvector, 1:ndim), save :: ix
!     integer(kind=4), dimension(1:nvector),         save :: dummy_state
!     integer(kind=4), dimension(1:nvector),         save :: parent_cell_level, parent_cell_index
!     real(dp),   dimension(1:nvector, 0:1, 1:ndim), save :: cloud_boundary
!     real(dp),        dimension(1:nvector),         save :: vol, delta
!     real(dp),        dimension(1:nvector, 1:3),    save :: xpart_cart
!     logical,         dimension(1:nvector),         save :: ok
!     integer,         dimension(1:ndim),            save :: ind
!     integer,  save :: idim, nx_loc, ind_cloud, ip
!     real(dp), save :: dx, dx_loc, scale, vol_loc, pos_to_cart
!     integer(kind=8) :: grid_size

!     grid_size = 2**grid_level    
!     nx_loc=(icoarse_max-icoarse_min+1)
!     scale=boxlen/dble(nx_loc)
!     dx = 0.5D0**grid_level
!     dx_loc=dx*scale
!     vol_loc=dx_loc**ndim
    

!     ! Convert particle coordinates in code units
!     ! into "cartesian" coordinates at grid_level
!     pos_to_cart = 2.0_dp**grid_level / dble(boxlen)
!     xpart_cart = xpart * pos_to_cart

!     ! compute distances of cloud boundary from nearest "integer coordinate"
!     do idim=1,ndim       

!        ! upper/right/front boundary of the cloud
!        do ip=1,np
!           cloud_boundary(ip,1,idim) = xpart_cart(ip, idim) + 0.5D0
!        end do

!        ! upper/rigt/front boundary rel to nearest integer (type conversion here...)
!        do ip=1,np
!           cloud_boundary(ip,1,idim) = cloud_boundary(ip,1,idim) - floor(cloud_boundary(ip,1,idim), kind=8)
!        end do
       
!        ! lower/left/back boundary rel to nearest integer
!        do ip=1,np
!           cloud_boundary(ip,0,idim) = 1._dp - cloud_boundary(ip,1,idim)
!        end do
!     end do
    
! #if NDIM<3
!     write(*,*)'add non-3D version of this routine'
!     stop
! #endif
! #if NDIM==1
!     ! Loop cloud/cell intersections
!     do ind_cloud = 0, 1
!        ind(1) = ind_cloud 
       
!        ! Compute cloud volume
!        do ip=1,np
!           vol(ip) = cloud_boundary(ip,ind(1),1) * &
!                cloud_boundary(ip,ind(2),2) 
!        end do
! #endif
! #if NDIM==2
!     ! Loop cloud/cell intersections
!     do ind_cloud = 0, 3
!        ind(1) = ind_cloud/2
!        ind(2) = mod(ind_cloud,2)
       
!        ! Compute cloud volume
!        do ip=1,np
!           vol(ip) = cloud_boundary(ip,ind(1),1) * &
!                cloud_boundary(ip,ind(2),2) 
!        end do
! #endif
! #if NDIM==3
!     ! Loop cloud/cell intersections
!     do ind_cloud = 0, 7
!        ind(1) = ind_cloud/4
!        ind(2) = mod(ind_cloud,4)/2
!        ind(3) = mod(mod(ind_cloud,4),2)
       
!        ! Compute cloud volume
!        do ip=1,np
!           vol(ip) = cloud_boundary(ip,ind(1),1) * &
!                cloud_boundary(ip,ind(2),2) * &
!                cloud_boundary(ip,ind(3),3) 
!        end do
! #endif

!        ! Compute cloud corner offset from cloud center
!        delta(1:ndim) = ind(1:ndim) - 0.5D0       

!        ! Get cell indices which are covered by cloud
!        ! (cartesian key -> hilbert key -> cell index)
!        ! TODO: Add support for non-periodic boundaries
!        ! TODO: Check boundary behaviour - currently particles sitting in cell touching the boundary cause
!        ! deviations from the old code.
!        do idim = 1, ndim
!           do ip = 1, np
!              ix(ip,idim) = floor(xpart_cart(ip,idim) + delta(idim), kind = 8)
!           end do
!        end do
!        do idim = 1, ndim
!           do ip = 1, np
!              if (ix(ip, idim) >= grid_size)then
!                 ix(ip, idim) = ix(ip, idim) - grid_size
!              else if (ix(ip, idim) < 0) then
!                 ix(ip, idim) = ix(ip, idim) + grid_size
!              end if
!           end do
!        end do
       
!        ! call hilbert3d(ix(1:np,1), ix(1:np,2), ix(1:np,3), &
!        !      cloud_hkey(1:np, 2), cloud_hkey(1:np, 1), cloud_hkey(1:np, 0), &
!        !      dummy_state, 0, grid_level, np)
      
!        ! call get_cell_index_from_hilbertkey(parent_cell_index(1:np), parent_cell_level(1:np), &
!        !      cloud_hkey(1:np, 2), cloud_hkey(1:np, 1), cloud_hkey(1:np, 0), np, grid_level)


!        call get_cell_index_from_cartesian_hash(parent_cell_index(1:np), parent_cell_level(1:np), &
!             ix(1:np, 1), ix(1:np, 2), ix(1:np, 3), grid_level, np, grid_level)  
       
!        ! Exclude cloud fraction which lies in coarser level
!        do ip = 1, np        
!           ok(ip) = (parent_cell_level(ip) == grid_level)
!        end do

!        ! Add to number density which is stored in phi
!        do ip=1,np
!           if(ok(ip))then
!              phi(parent_cell_index(ip)) = phi(parent_cell_index(ip)) + vol(ip)
!           end if
!        end do

!        ! compute delta rho and add to rho
!        do ip = 1, np
!           vol(ip) = mpart(ip) * vol(ip) / vol_loc
!        end do

!        do ip = 1, np
!           if (ok(ip)) then
!              rho(parent_cell_index(ip)) = rho(parent_cell_index(ip)) + vol(ip)
!           end if
!        end do

!     end do ! end loop over cloud/cell intersections
!   end subroutine cic_amr
