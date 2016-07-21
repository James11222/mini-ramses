!###############################################################
!###############################################################
!###############################################################
!###############################################################
subroutine get_threetondim_nbor_parent_cell_vec(hash_key,hash_dict,igrid_nbor,ind_nbor,flush_cache,fetch_cache,ngrid)
  use amr_commons
  use hash
  implicit none
  integer::ngrid
  logical::flush_cache,fetch_cache
  integer(kind=8),dimension(1:nvector,0:ndim)::hash_key
  type(hash_table)::hash_dict
  integer,dimension(1:nvector,1:threetondim)::igrid_nbor,ind_nbor
  !
  ! This routine computes and acquire the 3**ndim neighboring father cells 
  ! for the input hash_key. The output arrays are the father cells
  ! parent oct indices and their associated cell indices within the oct.
  ! The corresponding data can be accessed using: grid(igrid)%data(ind).
  ! If the grid index is zero, it means that this oct does not exist.
  ! Note that the 2**ndim grids are all locked if remote.
  !
  integer,dimension(1:nvector,1:twotondim),save::igrid_twotondim_nbor
  integer(kind=8),dimension(1:nvector,0:ndim),save::hash_nbor,hash_ref,hash_father
  integer(kind=8),dimension(1:nvector,1:ndim),save::ii
  integer,dimension(1:nvector),save::ind,ipos
  integer,dimension(1:3,1:8),save::shift_oct=reshape(&
       & (/-1,-1,-1,+1,-1,-1,-1,+1,-1,+1,+1,-1,&
       &   -1,-1,+1,+1,-1,+1,-1,+1,+1,+1,+1,+1/),(/3,8/))
  integer::i1,j1,k1
  integer,save::i1min=-1
  integer,save::i1max=+1
  integer,save::j1min=0*(1-ndim/2)-1*(ndim/2)
  integer,save::j1max=0*(1-ndim/2)+1*(ndim/2)
  integer,save::k1min=0*(1-ndim/3)-1*(ndim/3)
  integer,save::k1max=0*(1-ndim/3)+1*(ndim/3)
  integer::i,idim,ilevel,inbor

  ilevel=hash_key(1,0)

  hash_father(1:ngrid,0)=hash_key(1,0)-1

  ! Gather twotondim neighboring father grids
  do inbor=1,twotondim
     do idim=1,ndim
        do i=1,ngrid
           hash_nbor(i,idim)=hash_key(i,idim)+shift_oct(idim,inbor)
        end do
     end do
     ! Periodic boundary conditons
     do idim=1,ndim
        do i=1,ngrid
           if(hash_nbor(i,idim)<0)hash_nbor(i,idim)=ckey_max(ilevel)-1
           if(hash_nbor(i,idim)==ckey_max(ilevel))hash_nbor(i,idim)=0
        enddo
     end do
     do i=1,ngrid
        hash_father(i,1:ndim)=hash_nbor(i,1:ndim)/2
     end do

     ! Store lower left neighbor coordinates 
     if(inbor==1)then
        do i=1,ngrid
           hash_ref(i,1:ndim)=hash_father(i,1:ndim)
        end do
     end if

     ! Get grid into memory
     call get_grid_vec(hash_father,hash_dict,ipos,flush_cache,fetch_cache,ngrid)

     ! Store results into array
     do i=1,ngrid
        igrid_twotondim_nbor(i,inbor)=ipos(i)
     end do
  end do
     
  ! Deal with neighboring father cells
  inbor=0
  do k1=k1min,k1max
     do j1=j1min,j1max
        do i1=i1min,i1max           
           inbor=inbor+1
           do i=1,ngrid
#if NDIM>0
              hash_nbor(i,1)=hash_key(i,1)+i1
#endif
#if NDIM>1
              hash_nbor(i,2)=hash_key(i,2)+j1
#endif
#if NDIM>2
              hash_nbor(i,3)=hash_key(i,3)+k1
#endif
           end do
           ! Periodic boundary conditons
           do idim=1,ndim
              do i=1,ngrid
                 if(hash_nbor(i,idim)<0)hash_nbor(i,idim)=ckey_max(ilevel)-1
                 if(hash_nbor(i,idim)==ckey_max(ilevel))hash_nbor(i,idim)=0
              enddo
           end do
           ! Compute neighboring cell index
           do idim=1,ndim
              do i=1,ngrid
                 hash_father(i,idim)=hash_nbor(i,idim)/2
              end do
           end do
           do idim=1,ndim
              do i=1,ngrid
                 ii(i,idim)=hash_nbor(i,idim)-2*hash_father(i,idim)
              end do
           end do
           ind(1:ngrid)=1
           do idim=1,ndim
              do i=1,ngrid
                 ind(i)=ind(i)+2**(idim-1)*ii(i,idim)
              end do
           end do
           do i=1,ngrid
              ind_nbor(i,inbor)=ind(i)
           end do

           ! Compute neighboring grid index
           do idim=1,ndim
              do i=1,ngrid
                 ii(i,idim)=hash_father(i,idim)-hash_ref(i,idim)
              end do
           end do
           ! Periodic boundary conditons
           do idim=1,ndim
              do i=1,ngrid
                 if(ii(i,idim)<0)ii(i,idim)=ii(i,idim)+ckey_max(ilevel-1)
              end do
           enddo
           ind(1:ngrid)=1
           do idim=1,ndim
              do i=1,ngrid
                 ind(i)=ind(i)+2**(idim-1)*ii(i,idim)
              end do
           end do
           do i=1,ngrid
              igrid_nbor(i,inbor)=igrid_twotondim_nbor(i,ind(i))
           end do
        end do
     end do
  end do

end subroutine get_threetondim_nbor_parent_cell_vec
!###############################################################
!###############################################################
!###############################################################
!###############################################################
subroutine get_twondim_nbor_parent_cell_vec(hash_key,hash_dict,igrid_nbor,ind_nbor,flush_cache,fetch_cache,ngrid)
  use amr_commons
  use hash
  implicit none
  integer::ngrid
  logical::flush_cache,fetch_cache
  integer(kind=8),dimension(1:nvector,0:ndim)::hash_key
  type(hash_table)::hash_dict
  integer,dimension(1:nvector,0:twondim)::igrid_nbor,ind_nbor
  !
  ! This routine computes and acquires the 2xndim neighboring father cells 
  ! for the input hash_key. The output arrays are the father cells
  ! parent oct indices and their associated cell indices within the oct.
  ! The corresponding data can be accessed using: grid(igrid)%data(ind).
  ! The first element (0) stands for the central father cell.
  ! If the grid index is zero, it means that this oct does not exist.
  ! Note that the parent grids are all locked if remote.
  !
  integer(kind=8),dimension(1:nvector,0:ndim),save::hash_nbor
  integer(kind=8),dimension(1:nvector,0:ndim),save::hash_father
  integer(kind=8),dimension(1:nvector,1:ndim),save::ii
  integer,dimension(1:nvector),save::ind,ipos
  integer,dimension(1:3,1:6),save::shift=reshape(&
       & (/-1,0,0,1,0,0,0,-1,0,0,1,0,0,0,-1,0,0,1/),(/3,6/))
  integer::i,idim,ilevel,inbor

  ilevel=hash_key(1,0)

  ! Deal with central parent cell first
  do i=1,ngrid
     hash_father(i,0)=hash_key(i,0)-1
  end do
  do idim=1,ndim
     do i=1,ngrid
        hash_father(i,idim)=hash_key(i,idim)/2
     end do
  end do
  do idim=1,ndim
     do i=1,ngrid
        ii(i,idim)=hash_key(i,idim)-2*hash_father(i,idim)
     end do
  end do
  ind(1:ngrid)=1
  do idim=1,ndim
     do i=1,ngrid
        ind(i)=ind(i)+2**(idim-1)*ii(i,idim)
     end do
  end do

  ! Get grid address into memory
  call get_grid_vec(hash_father,hash_dict,ipos,flush_cache,fetch_cache,ngrid)

  ! Store results in output array
  do i=1,ngrid
     igrid_nbor(i,0)=ipos(i)
     ind_nbor(i,0)=ind(i)
  end do
  
  ! Deal with neighboring father cells
  do inbor=1,twondim
     do idim=1,ndim
        do i=1,ngrid
           hash_nbor(i,idim)=hash_key(i,idim)+shift(idim,inbor)
        end do
     end do
     ! Periodic boundary conditons
     do idim=1,ndim
        do i=1,ngrid
           if(hash_nbor(i,idim)<0)hash_nbor(i,idim)=ckey_max(ilevel)-1
           if(hash_nbor(i,idim)==ckey_max(ilevel))hash_nbor(i,idim)=0
        end do
     enddo
     do idim=1,ndim
        do i=1,ngrid
           hash_father(i,idim)=hash_nbor(i,idim)/2
        end do
     end do
     do idim=1,ndim
        do i=1,ngrid
           ii(i,idim)=hash_nbor(i,idim)-2*hash_father(i,idim)
        end do
     end do
     ind(1:ngrid)=1
     do idim=1,ndim
        do i=1,ngrid
           ind(i)=ind(i)+2**(idim-1)*ii(i,idim)
        end do
     end do

     ! Get grid into memory
     call get_grid_vec(hash_father,hash_dict,ipos,flush_cache,fetch_cache,ngrid)

     ! Store results in output arrays
     do i=1,ngrid
        igrid_nbor(i,inbor)=ipos(i)
        ind_nbor(i,inbor)=ind(i)
     end do
  end do

end subroutine get_twondim_nbor_parent_cell_vec
!###############################################################
!###############################################################
!###############################################################
!###############################################################
subroutine get_parent_cell_vec(hash_key,hash_dict,igrid_prnt,ind_prnt,flush_cache,fetch_cache,ngrid)
  use amr_commons
  use hash
  implicit none
  integer::ngrid
  logical::flush_cache,fetch_cache
  integer(kind=8),dimension(1:nvector,0:ndim)::hash_key
  integer,dimension(1:nvector)::igrid_prnt,ind_prnt
  type(hash_table)::hash_dict
  !
  ! This routine acquires the parent cell of the grid
  ! corresponding to the input hash key. On output, the parent cell 
  ! is defined by the father grid address in memory stored in igrid_prnt,
  ! and the father cell location in the father grid stored in ind_prnt. 
  !
  integer(kind=8),dimension(1:nvector,0:ndim),save::hash_father
  integer(kind=8),dimension(1:nvector,1:ndim),save::ii
  integer::i,idim

  ! Get father level
  do i=1,ngrid
     hash_father(i,0)=hash_key(i,0)-1
  end do

  ! Get father grid Cartesian index
  do idim=1,ndim
     do i=1,ngrid
        hash_father(i,idim)=hash_key(i,idim)/2
     end do
  end do

  ! Get father cell coordinates in father grid
  do idim=1,ndim
     do i=1,ngrid
        ii(i,idim)=hash_key(i,idim)-2*hash_father(i,idim)
     end do
  end do

  ! Translate coordinates into position index
  ind_prnt(1:ngrid)=1
  do idim=1,ndim
     do i=1,ngrid
        ind_prnt(i)=ind_prnt(i)+2**(idim-1)*ii(i,idim)
     end do
  end do

  ! Get father grid memory location
  call get_grid_vec(hash_father,hash_dict,igrid_prnt,flush_cache,fetch_cache,ngrid)

end subroutine get_parent_cell_vec
!###############################################################
!###############################################################
!###############################################################
!###############################################################
subroutine get_grid_vec(hash_key,hash_dict,igrid,flush_cache,fetch_cache,ngrid)
  use amr_commons
  use hilbert
  use hash
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ngrid
  logical::flush_cache,fetch_cache
  integer,dimension(1:nvector)::igrid
  integer(kind=8),dimension(1:nvector,0:ndim)::hash_key
  type(hash_table)::hash_dict
  !
  ! This routine acquires the grid 
  ! corresponding to the input hash key.
  !
  integer::i,iz,n_zero,get_grid
  integer,dimension(1:nvector),save::ind_zero
 
#ifndef WITHOUTMPI
  if(ncpu>1)then
     ! If counter is good, check on incoming messages and perform actions
     if(mail_counter==32)then
        call check_mail(MPI_REQUEST_NULL,hash_dict)
        mail_counter=0
     endif
     mail_counter=mail_counter+1
  endif
#endif

  ! Access hash table
  call hash_get_vec(hash_dict,hash_key,igrid,ngrid)

#ifndef WITHOUTMPI
  if(ncpu>1)then

     ! Count cells with zero address
     n_zero=0
     do i=1,ngrid
        if(igrid(i).EQ.0)then
           n_zero=n_zero+1
           ind_zero(n_zero)=i
        endif
     end do

     ! If grid index is -1, then set it to 0
     ! This means we already know the remote grid does not exist
     do i=1,ngrid
        if(igrid(i).EQ.-1)then
           igrid(i)=0
        endif
     end do

     ! Now we know child_grid=0 locally
     ! Acquire grids one by one and lock them all
     do iz=1,n_zero
        i=ind_zero(iz)
        igrid(i)=get_grid(hash_key(i,0:ndim),hash_dict,flush_cache,fetch_cache)
        if(igrid(i)>ngridmax)locked(igrid(i)-ngridmax)=.true.
     end do

  endif
#endif

end subroutine get_grid_vec
!##############################################################
!##############################################################
!##############################################################
!##############################################################
