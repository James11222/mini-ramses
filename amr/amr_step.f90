recursive subroutine amr_step(ilevel,icount)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel,icount,ilev,icnt
  !-------------------------------------------------------------------!
  ! This routine is the adaptive-mesh/adaptive-time-step main driver. !
  ! Each routine is called using a specific order, don't change it,   !
  ! unless you check all consequences first                           !
  !-------------------------------------------------------------------!
  logical,save::first_step=.true.

  if(noct_tot(ilevel)==0)return
  if(verbose)write(*,999)icount,ilevel

  if(ilevel==levelmin.or.icount>1)then

  !---------------------
  ! Make new refinements
  !---------------------
  if(ilevel==levelmin.or.icount>1)then
                               call timer('refine','start')
     call refine_fine(ilevel)
                               call timer('load balance','start')
     call load_balance(ilevel)
  endif

  !------------------------------
  ! Balance particles across cpus
  ! Careful rho must have been called once !
  !------------------------------
  if(first_step)then
     first_step=.false.
  else
     if(ilevel==levelmin)then
                               call timer('particles','start')
        if(pic)call balance_part(ilevel)
     endif
  endif

  !------------------------
  ! Output results to files
  !------------------------
                               call timer('output','start')
  if(ilevel==levelmin)then
     if(mod(nstep_coarse,foutput)==0.or.aexp>=aout(iout).or.t>=tout(iout))then
        call dump_all
     endif
  endif
  
  !----------------------------
  ! Output frame to movie dump
  !----------------------------
  if(movie .and. ilevel==levelmin) then
     if(imov.le.imovout)then 
        if(aexp>=amovout(imov).or.t>=tmovout(imov))then
           call output_frame()
        endif
     endif
  end if

  !--------------------
  ! Poisson source term
  !--------------------
  if(poisson)then
                               call timer('rho','start')
     call rho_fine(ilevel)
  endif

  !---------------
  ! Gravity solver
  !---------------
#ifdef GRAV
  if(poisson)then
                               call timer('poisson','start')
     ! Remove gravity source term with half time step and old force
     if(hydro)then
        call synchro_hydro_fine(ilevel,nlevelmax,-0.5_dp)
     endif

     ! Save old potential for time-extrapolation at level boundaries
     call save_phi_old(ilevel)

     ! Gravity solver depends on level
     ! Multigrid [levelmin:cglevelmin-1], Conjugate gradient [cg_levelmin:nlevelmax]
     call multigrid(ilevel, cg_levelmin-1, icount)
     call phi_fine_cg(max(cg_levelmin,ilevel), merge(icount,1,cg_levelmin<=ilevel))

     ! Initial old potential
     if (nstep==0)call save_phi_old(ilevel)

     ! Compute gravitational acceleration
     call force_fine(ilevel,icount)

     ! Perform second kick for particles
                               call timer('particles','start')
     if(pic)call kick_drift_part(ilevel,nlevelmax,action_kick_only)

     ! Add gravity source term with half time step and new force
     if(hydro)then
                               call timer('poisson','start')
        call synchro_hydro_fine(ilevel,nlevelmax,+0.5_dp)
     end if

  end if
#endif

  !----------------------
  ! Compute new time step
  !----------------------
                               call timer('courant','start')
  call newdt_fine(ilevel)
  do ilev=max(ilevel,levelmin+1), nlevelmax
     dtnew(ilev)=MIN(dtnew(ilev-1)/real(nsubcycle(ilev-1)),dtnew(ilev))
  end do
  
  !-----------------------
  ! Set unew equal to uold
  !-----------------------
                               call timer('hydro - set unew','start')
  if(hydro)call set_unew(ilevel)

  end if
  !---------------------------
  ! Recursive call to amr_step
  !---------------------------
                               call timer('recursive call','start')
  if(ilevel<nlevelmax)then
     if(noct_tot(ilevel+1)>0)then
        if(nsubcycle(ilevel)==2)then
           call amr_step(ilevel+1,1)
           call amr_step(ilevel+1,2)
        else
           call amr_step(ilevel+1,1)
        endif
     else 
        ! Otherwise, update time and finer level time-step
        dtold(ilevel+1)=dtnew(ilevel)/dble(nsubcycle(ilevel))
        dtnew(ilevel+1)=dtnew(ilevel)/dble(nsubcycle(ilevel))
        call update_time(ilevel)
     end if
  else
     call update_time(ilevel)
  end if

  !-----------
  ! Hydro step
  !-----------

  if (ilevel==levelmin .or. (icount==1 .and. nsubcycle(ilevel-1)==2)) then
  !-------------------------------
  ! Update coarser level time-step
  !-------------------------------
  do ilev=nlevelmax,ilevel,-1
     if (noct_tot(ilev)==0) cycle
     icnt = merge(2,1,ilev > ilevel)
     if(ilev>levelmin)then
        if(nsubcycle(ilev-1)==1)dtnew(ilev-1)=dtnew(ilev)
        if(icnt==2)dtnew(ilev-1)=dtold(ilev)+dtnew(ilev)
     end if
  end do

  if(hydro)then

     ! Hyperbolic solver
                               call timer('hydro - godunov','start')
     call godunov_fine(ilevel)

     ! Add gravity source terms to unew with half time step
                               call timer('poisson - synchro','start')
     if(poisson)call add_gravity_source_terms(ilevel)

     ! Set uold equal to unew
                               call timer('hydro - set uold','start')
     call set_uold(ilevel)

     ! Add gravity source terms to uold with half time step
     ! to complete the time step (will be removed later)
                               call timer('poisson - synchro','start')
     if(poisson)call synchro_hydro_fine(ilevel,nlevelmax,+0.5_dp)

     ! Restriction operator
                               call timer('hydro - upload','start')
     call upload_fine(ilevel)
  endif

  !----------------------------
  ! Compute cooling/heating
  !----------------------------
                               call timer('cooling','start')
  if(cooling)call cooling_fine(ilev)

  !-------------------------------------------
  ! Perform first kick and drift for particles
  !-------------------------------------------
                               call timer('particles','start')
  if(pic)call kick_drift_part(ilevel,nlevelmax,action_kick_drift)

  do ilev=nlevelmax,ilevel,-1
    if (noct_tot(ilev)==0) cycle
    icnt = merge(2,1,ilev > ilevel)
  !-----------------------
  ! Compute refinement map
  !-----------------------
                               call timer('flag','start')
  if(.not.static) call flag_fine(ilev,icnt)

  enddo
                               call timer('recursive call','start')

  end if

999 format(' Entering amr_step',i1,' for level',i2)

end subroutine amr_step
