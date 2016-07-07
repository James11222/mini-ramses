subroutine amr_step
  USE amr_parameters, only :: tg_max, max_active_time_group, dt_old, dt_new
  implicit none
  !=============================================================================
  ! Basic time step driver routine for looping through time-step groups in one explicit loop
  !=============================================================================
  ! Time is measured as integer steps, with each time group having 2^tgroup sub-steps inside the full step.
  ! To allow the number of time-groups to change during a step but use exact integer arithmetics,
  ! one full step is equal to a very large integer (base=2^63) allowing effectively for 63 time-group levels.
  integer(kind=8), parameter         :: base=2**tg_max ! integer unit used for one full time-step
  integer(kind=8)                    :: t_int          ! current time inside timestep
  integer(kind=8), dimension(tg_max) :: now, future    ! integer current and future time for each time-group
  integer(kind=8), dimension(tg_max) :: dt_int         ! integer dt for each time group
  integer(kind=8) :: tgroup, tgroup_fut
  !
  ! Set up integer time-steps
  t_int = 0
  now   = 0
  do ig=1_8,tg_max
     future(ig) = base / 2_8**(ig-1)
  enddo
  dt_int = future - now

  ! Main loop over one full time-step
  do while (t_int .ne. base)
     !-------------------------------
     ! Determine which groups have reached current time, and have to prepare things
     ! for updating the time-step. This is done at time "now" (e.g. calculating the gravitational potential).
     !-------------------------------
     do tgroup=1,tg_max
        if (now(tgroup)==t_int) exit
     end do
     call amr_prepare(tgroup)

     !----------------------
     ! Compute new real time step and update time-group distribution
     !----------------------
                                  call timer('courant','start')
     call newdt_fine(tgroup)
     call refine_time_groups(tgroup)
     do ig=max(tgroup,2), max_active_time_group
        dtnew(ig)=MIN(dtnew(ig-1)/2.,dtnew(ig))
     end do

     !-------------------------------
     ! Update integer time-step with dt of the highest active time-group
     !-------------------------------
     t_int = t_int + dt(max_active_time_group)

     !-------------------------------
     ! Update real time of simulation with real dt of higest active time-group
     !-------------------------------
     !dtold(ilevel+1)=dtnew(ilevel)/dble(nsubcycle(ilevel))
     !dtnew(ilevel+1)=dtnew(ilevel)/dble(nsubcycle(ilevel))
     call update_time(max_active_time_group)

     !-------------------------------
     ! Determine which groups will evolve to current time
     !-------------------------------
     do tgroup_fut=1,tg_max
        if (future(tgroup_fut)==t_int) exit
     end do

     !-------------------------------
     ! Update coarser time-groups time-step
     !-------------------------------
     do ig=max_active_time_group,tgroup_fut+1,-1
        dtnew(ig-1)=dtold(ig)+dtnew(ig)
     end do

     !-------------------------------
     ! Load balance cells across timegroups
     !-------------------------------
     if (do_timegroup) call load_balance(tgroup_fut,tgroup)

     !-------------------------------
     ! Evolve time groups
     !-------------------------------
     call amr_evolve(tgroup_fut)
                               call timer('recursive call','start')

     !-------------------------------
     ! Update integer timestep counters to account for updated time
     !-------------------------------
     now(tgroup:tg_max)    = t_int
     future(tgroup:tg_max) = future(tgroup:tg_max) + dt(tgroup:tg_max)
  end do
end subroutine amr_step

subroutine amr_prepare(tgroup)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer, intent(in)::tgroup
  !-----------------------------------------------------------------!
  ! This is the first adaptive-mesh/adaptive-time-step main driver. !
  ! It prepares everything for the actual time-step.                !
  ! Each routine is called using a specific order, don't change it, !
  ! unless you check all consequences first.                        !
  !-----------------------------------------------------------------!
  integer::ilevel,ilev,ig,icount
  logical,save::first_step=.true.

  if(verbose)write(*,999)tgroup

  !---------------------
  ! check if timegroups are bound to levels through nsubcycle
  !---------------------
  if (.not. do_timegroups) then
     ig=1; ilevel=levelmin
     do while (ig < tgroup)
        ilevel = ilevel + 1
        ig = ig + nsubcycle(ilevel-1)-1
     end do
  endif

  icount = merge(1,2,tgroup==1)

  !---------------------
  ! Make new refinements
  !---------------------
                               call timer('refine','start')
  call refine_fine(ilevel)
                               call timer('load balance','start')
  call load_balance(ilevel)
  !call load_balance(tgroup,tgroup)

  !------------------------------
  ! Balance particles across cpus
  ! Careful rho must have been called once !
  !------------------------------
  if(first_step)then
     first_step=.false.
  else
     if(tgroup==1)then
                               call timer('particles','start')
        if(pic)call balance_part(ilevel)
     endif
  endif

  !------------------------
  ! Output results to files
  !------------------------
                               call timer('output','start')
  if(tgroup==1)then
     if(mod(nstep_coarse,foutput)==0.or.aexp>=aout(iout).or.t>=tout(iout))then
        call dump_all
     endif
  endif
  
  !----------------------------
  ! Output frame to movie dump
  !----------------------------
  if(movie .and. tgroup==1) then
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
        call synchro_hydro_fine(ilevel,-0.5_dp)
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
     if(pic)call kick_drift_part(ilevel,action_kick_only)

     ! Add gravity source term with half time step and new force
     if(hydro)then
                               call timer('poisson','start')
        call synchro_hydro_fine(ilevel,+0.5_dp)
     end if

  end if
#endif
  
  !-----------------------
  ! Set unew equal to uold
  !-----------------------
                               call timer('hydro - set unew','start')
  if(hydro)call set_unew(ilevel)

end subroutine amr_prepare

subroutine amr_evolve(tgroup)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::tgroup
  !------------------------------------------------------------------!
  ! This is the second adaptive-mesh/adaptive-time-step main driver. !
  ! It evolves all the cells in time-groups >= tgroup                !
  ! Each routine is called using a specific order, don't change it,  !
  ! unless you check all consequences first.                         !
  !------------------------------------------------------------------!
  integer::ilevel,ilev,ig

  !---------------------
  ! check if timegroups are bound to levels through nsubcycle
  !---------------------
  if (.not. do_timegroups) then
     ig=1; ilevel=levelmin
     do while (ig < tgroup)
        ilevel = ilevel + 1
        ig = ig + nsubcycle(ilevel-1)-1
     end do
  endif

  !-----------
  ! Hydro step
  !-----------
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
     if(poisson)call synchro_hydro_fine(ilevel,+0.5_dp)

     ! Restriction operator
                               call timer('hydro - upload','start')
     call upload_fine(ilevel)
  endif

  !----------------------------
  ! Compute cooling/heating
  !----------------------------
                               call timer('cooling','start')
  if(cooling)call cooling_fine(ilevel)

  !-------------------------------------------
  ! Perform first kick and drift for particles
  !-------------------------------------------
                               call timer('particles','start')
  if(pic)call kick_drift_part(ilevel,action_kick_drift)

  do ilev=nlevelmax,ilevel,-1
    if (noct_tot(ilev)==0) cycle
  !-----------------------
  ! Compute refinement map
  !-----------------------
                               call timer('flag','start')
  if(.not.static) call flag_fine(ilev,ilev==ilevel)

  enddo
end subroutine amr_evolve

