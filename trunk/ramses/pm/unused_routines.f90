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
