subroutine output_poisson(filename)
  use amr_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  character(LEN=80)::filename

  integer::ilevel,igrid,idim,ilun,ngrid,i
  character(LEN=5)::nchar
  character(LEN=80)::fileloc
  integer,parameter::nio_buffer=1024
  real(dp),dimension(1:twotondim,1:ndim+1,1:nio_buffer)::io_buf

  if(verbose)write(*,*)'Entering output_poisson'
  ilun=ncpu+myid+10
  call title(myid,nchar)
  fileloc=TRIM(filename)//TRIM(nchar)
  open(unit=ilun,file=fileloc,access="stream"&
       & ,status="replace",action="write",form='unformatted')
  write(ilun)ndim
  write(ilun)ndim+1
  write(ilun)levelmin
  write(ilun)nlevelmax
  do ilevel=levelmin,nlevelmax
     write(ilun)noct(ilevel)
  enddo
#ifdef GRAV
  do ilevel=levelmin,nlevelmax
     do igrid=head(ilevel),tail(ilevel),nio_buffer
        ngrid=MIN(nio_buffer,tail(ilevel)-igrid+1)
        if(ngrid==nio_buffer)then
           do i=1,ngrid
              io_buf(1:twotondim,1,i)=grid(igrid+i-1)%phi
              do idim=1,ndim
                 io_buf(1:twotondim,idim+1,i)=grid(igrid+i-1)%f(1:twotondim,idim)
              end do
           end do
           write(ilun)io_buf
        else
           do i=1,ngrid
              write(ilun)grid(igrid+i-1)%phi
              write(ilun)grid(igrid+i-1)%f
           end do
        end if
     end do
  enddo
#endif
  close(ilun)
     
end subroutine output_poisson





