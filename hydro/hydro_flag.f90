module hydro_flag_module
contains
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
subroutine hydro_flag(s,ilevel)
  use amr_parameters, only: ndim,twotondim,twondim,dp
  use amr_commons, only: oct,nbor
  use ramses_commons, only: ramses_t
  use hydro_parameters, only: nvar
  use cache_commons
  use cache
  use nbors_utils
  use boundaries, only: init_bound_refine
  implicit none
  type(ramses_t)::s
  integer::ilevel
  ! -------------------------------------------------------------------
  ! This routine flag for refinement cells that satisfies
  ! some user-defined physical criteria at the level ilevel. 
  ! -------------------------------------------------------------------
  integer,dimension(1:3,1:8),save::iii=reshape(&
       & (/0,0,0,1,0,0,0,1,0,1,1,0,0,0,1,1,0,1,0,1,1,1,1,1/),(/3,8/))
  integer,dimension(1:3,1:6),save::shift=reshape(&
       & (/-1,0,0,1,0,0,0,-1,0,0,1,0,0,0,-1,0,0,1/),(/3,6/))
  integer::igrid,ind,idim,ivar,i_nbor
  integer::igridd,igridg,icelld,icellg,igridp,icellp
  integer,dimension(1:twondim)::igridn,icelln
  integer(kind=8),dimension(0:ndim)::hash_key,hash_nbor
  real(dp),dimension(1:nvar)::uug,uum,uud
  logical::ok
  type(nbor),dimension(1:twondim)::gridn
  type(oct),pointer::gridp
  type(msg_realdp)::dummy_realdp

#ifdef HYDRO

  associate(r=>s%r,g=>s%g,m=>s%m)

  if(    r%err_grad_d==-1.0.and.&
       & r%err_grad_p==-1.0.and.&
       & r%err_grad_u==-1.0)return

  hash_key(0)=ilevel+1

  call open_cache(s,table=m%grid_dict,data_size=storage_size(m%grid(1))/32,&
                     hilbert=m%domain, pack_size=storage_size(dummy_realdp)/32,&
                     pack=pack_fetch_hydro,unpack=unpack_fetch_hydro,&
                     bound=init_bound_refine)

  ! Loop over active grids
  do igrid=m%head(ilevel),m%tail(ilevel)

     ! Loop over cells
     do ind=1,twotondim

        ! Compute cell hash key
        hash_key(1:ndim)=2*m%grid(igrid)%ckey(1:ndim)+iii(1:ndim,ind)

        ! Initialize refinement to false
        ok=.false.

        ! If a neighbor cell does not exist,
        ! replace it by its father cell
        do i_nbor=1,twondim
           hash_nbor(0)=hash_key(0)
           ! Periodic boundary conditions
           do idim=1,ndim
              hash_nbor(idim)=hash_key(idim)+shift(idim,i_nbor)
              if(r%periodic(idim))then
                 if(hash_nbor(idim)<m%box_ckey_min(idim,ilevel+1))hash_nbor(idim)=m%box_ckey_max(idim,ilevel+1)-1
                 if(hash_nbor(idim)>=m%box_ckey_max(idim,ilevel+1))hash_nbor(idim)=m%box_ckey_min(idim,ilevel+1)
              endif
           enddo
           call get_parent_cell(s,hash_nbor,m%grid_dict,gridp,icellp,flush_cache=.false.,fetch_cache=.true.,lock=.true.)
           if(associated(gridp))then
              gridn(i_nbor)%p=>gridp
              icelln(i_nbor)=icellp
           else
              hash_nbor(0)=hash_nbor(0)-1
              hash_nbor(1:ndim)=hash_nbor(1:ndim)/2
              call get_parent_cell(s,hash_nbor,m%grid_dict,gridp,icellp,flush_cache=.false.,fetch_cache=.true.,lock=.true.)
              gridn(i_nbor)%p=>gridp
              icelln(i_nbor)=icellp
           endif
        end do

#if NDOF>1
        ! First dimension
        idim=1; j=1; k=1
#if NDIM>2
        do k=1,ndof
#endif
#if NDIM>1
        do j=1,ndof
#endif
        do i=1,ndof
           if(i==1)then
              il=ndof
              ir=2
              idof=i+(j-1)*ndof+(k-1)*ndof*ndof
              idofl=il+(j-1)*ndof+(k-1)*ndof*ndof
              idofr=ir+(j-1)*ndof+(k-1)*ndof*ndof
              icellg=icelln(2*idim-1)
              do ivar=1,nvar
                 uug(ivar)=gridn(2*idim-1)%p%uold(idofl,icellg,ivar)
                 uum(ivar)=m%grid(igrid)%uold(idof,ind,ivar)
                 uud(ivar)=m%grid(igrid)%uold(idofr,ind,ivar)
              end do
           else if(i==ndof)then
              il=ndof-1
              ir=1
              idof=i+(j-1)*ndof+(k-1)*ndof*ndof
              idofl=il+(j-1)*ndof+(k-1)*ndof*ndof
              idofr=ir+(j-1)*ndof+(k-1)*ndof*ndof
              icelld=icelln(2*idim)
              do ivar=1,nvar
                 uug(ivar)=m%grid(igrid)%uold(idofl,ind,ivar)
                 uum(ivar)=m%grid(igrid)%uold(idof,ind,ivar)
                 uud(ivar)=uud(ivar)=gridn(2*idim)%p%uold(idofr,icelld,ivar)
              end do
           else
              il=i-1
              ir=i+1
              idof=i+(j-1)*ndof+(k-1)*ndof*ndof
              idofl=il+(j-1)*ndof+(k-1)*ndof*ndof
              idofr=ir+(j-1)*ndof+(k-1)*ndof*ndof
              do ivar=1,nvar
                 uug(ivar)=m%grid(igrid)%uold(idofl,ind,ivar)
                 uum(ivar)=m%grid(igrid)%uold(idof,ind,ivar)
                 uud(ivar)=m%grid(igrid)%uold(idofr,ind,ivar)
              end do
           endif
           call hydro_refine(r,uug,uum,uud,ok)
        end do
#if NDIM>1
        end do
#endif
#if NDIM>2
        end do
#endif
#if NDIM>1
        ! Second dimension
        idim=2; k=1
#if NDIM>2
        do k=1,ndof
#endif
        do j=1,ndof
        do i=1,ndof
           if(j==1)then
              jl=ndof
              jr=2
              idof=i+(j-1)*ndof+(k-1)*ndof*ndof
              idofl=i+(jl-1)*ndof+(k-1)*ndof*ndof
              idofr=i+(jr-1)*ndof+(k-1)*ndof*ndof
              icellg=icelln(2*idim-1)
              do ivar=1,nvar
                 uug(ivar)=gridn(2*idim-1)%p%uold(idofl,icellg,ivar)
                 uum(ivar)=m%grid(igrid)%uold(idof,ind,ivar)
                 uud(ivar)=m%grid(igrid)%uold(idofr,ind,ivar)
              end do
           else if(j==ndof)then
              jl=ndof-1
              jr=1
              idof=i+(j-1)*ndof+(k-1)*ndof*ndof
              idofl=i+(jl-1)*ndof+(k-1)*ndof*ndof
              idofr=i+(jr-1)*ndof+(k-1)*ndof*ndof
              icelld=icelln(2*idim)
              do ivar=1,nvar
                 uug(ivar)=m%grid(igrid)%uold(idofl,ind,ivar)
                 uum(ivar)=m%grid(igrid)%uold(idof,ind,ivar)
                 uud(ivar)=gridn(2*idim)%p%uold(idofr,icelld,ivar)
              end do
           else
              jl=j-1
              jr=j+1
              idof=i+(j-1)*ndof+(k-1)*ndof*ndof
              idofl=i+(jl-1)*ndof+(k-1)*ndof*ndof
              idofr=i+(jr-1)*ndof+(k-1)*ndof*ndof
              do ivar=1,nvar
                 uug(ivar)=m%grid(igrid)%uold(idofl,ind,ivar)
                 uum(ivar)=m%grid(igrid)%uold(idof,ind,ivar)
                 uud(ivar)=m%grid(igrid)%uold(idofr,ind,ivar)
              end do
           endif
           call hydro_refine(r,uug,uum,uud,ok)
        end do
        end do
#if NDIM>2
        end do
#endif
#endif
#if NDIM>2
        ! Third dimension
        idim=3
        do k=1,ndof
        do j=1,ndof
        do i=1,ndof
           if(k==1)then
              kl=ndof
              kr=2
              idof=i+(j-1)*ndof+(k-1)*ndof*ndof
              idofl=i+(j-1)*ndof+(kl-1)*ndof*ndof
              idofr=i+(j-1)*ndof+(kr-1)*ndof*ndof
              icellg=icelln(2*idim-1)
              do ivar=1,nvar
                 uug(ivar)=gridn(2*idim-1)%p%uold(idofl,icellg,ivar)
                 uum(ivar)=m%grid(igrid)%uold(idof,ind,ivar)
                 uud(ivar)=m%grid(igrid)%uold(idofr,ind,ivar)
              end do
           else if(k==ndof)then
              kl=ndof-1
              kr=1
              idof=i+(j-1)*ndof+(k-1)*ndof*ndof
              idofl=i+(j-1)*ndof+(kl-1)*ndof*ndof
              idofr=i+(j-1)*ndof+(kr-1)*ndof*ndof
              icelld=icelln(2*idim)
              do ivar=1,nvar
                 uug(ivar)=m%grid(igrid)%uold(idofl,ind,ivar)
                 uum(ivar)=m%grid(igrid)%uold(idof,ind,ivar)
                 uud(ivar)=gridn(2*idim)%p%uold(idofr,icelld,ivar)
              end do
           else
              kl=k-1
              kr=k+1
              idof=i+(j-1)*ndof+(k-1)*ndof*ndof
              idofl=i+(j-1)*ndof+(kl-1)*ndof*ndof
              idofr=i+(j-1)*ndof+(kr-1)*ndof*ndof
              do ivar=1,nvar
                 uug(ivar)=m%grid(igrid)%uold(idofl,ind,ivar)
                 uum(ivar)=m%grid(igrid)%uold(idof,ind,ivar)
                 uud(ivar)=m%grid(igrid)%uold(idofr,ind,ivar)
              end do
           endif
           call hydro_refine(r,uug,uum,uud,ok)
        end do
        end do
        end do
#endif
#else
        ! Loop over dimensions
        do idim=1,ndim
           ! Gather hydro variables
           do ivar=1,nvar
              icellg=icelln(2*idim-1)
              icelld=icelln(2*idim  )
              uug(ivar)=gridn(2*idim-1)%p%uold(icellg,ivar)
              uum(ivar)=m%grid(igrid)%uold(ind,ivar)
              uud(ivar)=gridn(2*idim)%p%uold(icelld,ivar)
           end do
           call hydro_refine(r,uug,uum,uud,ok)
        end do
#endif
        
        do i_nbor=1,twondim
           call unlock_cache(s,gridn(i_nbor)%p)
        end do

        ! Count only newly flagged cells
        if(m%grid(igrid)%flag1(ind)==0.and.ok)g%nflag=g%nflag+1
        if(ok)m%grid(igrid)%flag1(ind)=1

     end do
     ! End loop over cells
  end do
  ! End loop over grids

  call close_cache(s,m%grid_dict)

  end associate

#endif

end subroutine hydro_flag
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
subroutine pack_fetch_hydro(grid,msg_size,msg_array)
  use amr_parameters, only: ndim,twotondim
  use hydro_parameters, only: nvar
  use amr_commons, only: oct
  use cache_commons, only: msg_realdp
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  integer::ind,ivar
  type(msg_realdp)::msg

  do ind=1,twotondim
     if(grid%refined(ind))then
        msg%int4(ind)=1
     else
        msg%int4(ind)=0
     endif
  end do
  
#ifdef HYDRO
  do ivar=1,nvar
     do ind=1,twotondim
        msg%realdp(ind,ivar)=grid%uold(ind,ivar)
     end do
  end do
#endif

  msg_array=transfer(msg,msg_array)

end subroutine pack_fetch_hydro
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
subroutine unpack_fetch_hydro(grid,msg_size,msg_array,hash_key)
  use amr_parameters, only: ndim,twotondim
  use hydro_parameters, only: nvar
  use amr_commons, only: oct
  use cache_commons, only: msg_realdp
  type(oct)::grid
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array
  integer(kind=8),dimension(0:ndim)::hash_key

  integer::ind,ivar
  type(msg_realdp)::msg

  grid%lev=hash_key(0)
  grid%ckey(1:ndim)=hash_key(1:ndim)
  msg=transfer(msg_array,msg)

  do ind=1,twotondim
     if(msg%int4(ind)==1)then
        grid%refined(ind)=.true.
     else
        grid%refined(ind)=.false.
     endif
  end do
  
#ifdef HYDRO
  do ivar=1,nvar
     do ind=1,twotondim
        grid%uold(ind,ivar)=msg%realdp(ind,ivar)
     end do
  end do
#endif

end subroutine unpack_fetch_hydro
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
end module hydro_flag_module
