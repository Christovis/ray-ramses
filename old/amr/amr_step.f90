recursive subroutine amr_step(ilevel,icount)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  use ray_commons
  use ray_parameters
#ifdef RT
  use rt_hydro_commons
  use SED_module
#endif
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel,icount
  !-------------------------------------------------------------------!
  ! This routine is the adaptive-mesh/adaptive-time-step main driver. !
  ! Each routine is called using a specific order, don't change it,   !
  ! unless you check all consequences first                           !
  !-------------------------------------------------------------------!
  integer::i,idim,ivar,ilevel2
  integer::ii ! Alex: for_exact_output_times
  logical::ok_defrag
  logical,save::first_step=.true.
  real(dp) :: ray_a_s
  real(dp) :: ray_aout                                             ! BAOJIU-17-02-2016

  ! LENSING STUFF !
  real(dp) :: dcom_ray ! this is a function
  real(dp) :: zexp,chi_wave,coverH0_ray,chi_start  ! Sownak: 02/09/2015
  real(dp) :: temp_dist_x,temp_dist_y,temp_dist_z
  real(dp) :: pi

  if(numbtot(1,ilevel)==0)return
    
  if(ray) then

     ray_a_s = 1.0d0/(1.0d0+ray_z_s) 
  
     ! Start: Sownak - ISW (28/10/2015)
     if(.not.allocated(aexp_old_ray)) then
        allocate(aexp_old_ray(1:nlevelmax))
        !aexp_old_ray = aexp
        aexp_old_ray = -1.0D0 ! abs_for_sync

        if(ray_multi_out.and.ray_nout.gt.0) then                   ! BAOJIU-17-02-2016
           ray_iout = 1                                            ! BAOJIU-17-02-2016
        end if                                                     ! BAOJIU-17-02-2016

        if(ray_do_isw) then
           allocate(aexp_new_ray(1:nlevelmax))
           aexp_new_ray = aexp
        endif
     end if
     ! End: Sownak

  end if

  if(verbose)write(*,999)icount,ilevel

  !-------------------------------------------
  ! Alex: for_exact_output_times
  ! Start of modification 
  !-------------------------------------------
  if(cosmo.and.dabs(t_next).lt.1.0D-8) then
     ! Find neighbouring scale factors
     ii=1
     do while(aexp_frw(ii)>aout(iout).and.ii<n_frw)
        ii=ii+1
     end do
     ! Interpolate expansion factor for the next step
     t_next=tau_frw(ii)*(aout(iout)-aexp_frw(ii-1))/(aexp_frw(ii)-aexp_frw(ii-1)) + tau_frw(ii-1)*(aout(iout)-aexp_frw(ii))/(aexp_frw(ii-1)-aexp_frw(ii))
  end if
  ! End of modification 

  !-------------------------------------------
  ! Make new refinements and update boundaries
  !-------------------------------------------
  if(levelmin.lt.nlevelmax .and..not. static)then
     if(ilevel==levelmin.or.icount>1)then
        do i=ilevel,nlevelmax
           if(i>levelmin)then

              !--------------------------
              ! Build communicators
              !--------------------------
              call build_comm(i)

              !--------------------------
              ! Update boundaries
              !--------------------------
              call make_virtual_fine_int(cpu_map(1),i)
              if(hydro)then
#ifdef SOLVERmhd
                 do ivar=1,nvar+3
#else
                 do ivar=1,nvar
#endif
                    call make_virtual_fine_dp(uold(1,ivar),i)
#ifdef SOLVERmhd
                 end do
#else
                 end do
#endif
                 if(simple_boundary)call make_boundary_hydro(i)
              end if
#ifdef RT
              if(rt)then
                 do ivar=1,nrtvar
                    call make_virtual_fine_dp(rtuold(1,ivar),i)
                 end do
                 if(simple_boundary)call rt_make_boundary_hydro(i)
              end if
#endif
              if(poisson)then
                 call make_virtual_fine_dp(phi(1),i)
                 do idim=1,ndim
                    call make_virtual_fine_dp(f(1,idim),i)
                 end do
                 if(simple_boundary)call make_boundary_force(i)
              end if
           end if

           !--------------------------
           ! Refine grids
           !--------------------------
           call refine_fine(i)
           
        end do
     end if
  end if

  !--------------------------
  ! Load balance
  !--------------------------
  ok_defrag=.false.
  if(levelmin.lt.nlevelmax)then
     if(ilevel==levelmin)then
        if(nremap>0)then
           ! Skip first load balance because it has been performed before file dump
           if(nrestart>0.and.first_step)then
              first_step=.false.
           else
              if(MOD(nstep_coarse,nremap)==0)then
                 call load_balance
                 call defrag
                 ok_defrag=.true.
              endif
           end if
        end if
     endif
  end if

  !-----------------
  ! Update sink cloud particle properties
  !-----------------
  if(sink)call update_cloud(ilevel,.false.)

  !-----------------
  ! Particle leakage
  !-----------------
  if(pic)call make_tree_fine(ilevel)
  
  !------------------------
  ! Output results to files
  !------------------------
  if(ilevel==levelmin)then
!     if(mod(nstep_coarse,foutput)==0.or.aexp>=aout(iout).or.t>=tout(iout))then        ! Alex: for_exact_output_times
     if(mod(nstep_coarse,foutput)==0.or.aexp>=aout(iout)-1.0D-8.or.t>=tout(iout))then ! Alex: for_exact_output_times
        if(.not.ok_defrag)then
           call defrag
        endif

        call dump_all

        !-------------------------------------------
        ! Alex: for_exact_output_times
        ! Start of modification 
        !-------------------------------------------
        if(cosmo) then
           ! Find neighbouring scale factors
           ii = 1
           do while(aexp_frw(ii) > aout(iout).and.ii<n_frw)
              ii = ii+1
           end do
           ! Interpolate expansion factor for the next step
           t_next = tau_frw(ii)*(aout(iout)-aexp_frw(ii-1))/(aexp_frw(ii)-aexp_frw(ii-1))+tau_frw(ii-1)*(aout(iout)-aexp_frw(ii))/(aexp_frw(ii-1)-aexp_frw(ii))
        end if
        ! End of modification

        if(gas_analytics) call gas_ana

        ! Run the clumpfinder
        if(clumpfind .and. ndim==3) call clump_finder(.true.)

        ! Dump lightcone
        if(lightcone) call output_cone()

     endif

     ! Important can't be done in sink routines because it must be done after dump all
     if(sink)acc_rate=0.

  endif

  !----------------------------
  ! Output frame to movie dump (without synced levels)
  !----------------------------
  if(movie) then
     if(aexp>=amovout(imov).or.t>=tmovout(imov))then
        call output_frame()
     endif
  end if


  !-----------------------------------------------------------
  ! Put here all stuffs that are done only at coarse time step
  !-----------------------------------------------------------
  if(ilevel==levelmin)then
     !----------------------------------------------------
     ! Kinetic feedback from giant molecular clouds
     !----------------------------------------------------
     if(hydro.and.star.and.eta_sn>0.and.f_w>0)call kinetic_feedback

  endif

  !--------------------
  ! Poisson source term
  !--------------------
  if(poisson)then
     !save old potential for time-extrapolation at level boundaries
     call save_phi_old(ilevel)

     ! Start: Sownak - ISW (28/10/2015) abs_for_sync: this if checking the sign
     ! of aexp_old_ray
     if(ray.and.ray_do_isw) then
        if(aexp_old_ray(ilevel) .lt. 0.0d0 .and. ilevel .ne. levelmin) then
           aexp_old_ray(ilevel) = (aexp_old_ray(ilevel-1) + aexp_new_ray(ilevel-1))/2.0d0
!          write(*,*) 'line1', ilevel, aexp, aexp_new_ray(ilevel), aexp_old_ray(ilevel)
        else
           aexp_old_ray(ilevel) = aexp_new_ray(ilevel)
!          write(*,*) 'line2', ilevel, aexp, aexp_new_ray(ilevel), aexp_old_ray(ilevel)
        end if
     end if
     ! End: Sownak

     call rho_fine(ilevel,icount)
  endif

  

  !-------------------------------------------
  ! Sort particles between ilevel and ilevel+1
  !-------------------------------------------
  if(pic)then
     ! Remove particles to finer levels
     call kill_tree_fine(ilevel)
     ! Update boundary conditions for remaining particles
     call virtual_tree_fine(ilevel)
  end if

  !---------------
  ! Gravity update
  !---------------
  if(poisson)then
 
     ! Remove gravity source term with half time step and old force
     if(hydro)then
        call synchro_hydro_fine(ilevel,-0.5*dtnew(ilevel))
     endif
     
     ! Compute gravitational potential
     if(ilevel>levelmin)then
        if(ilevel .ge. cg_levelmin) then
           call phi_fine_cg(ilevel,icount)
        else
           call multigrid_fine(ilevel,icount)
        end if
     else
        call multigrid_fine(levelmin,icount)
     end if
     
     ! Start: Sownak - ISW (28/10/2015)
     if(ray.and.ray_do_isw) then
        aexp_new_ray(ilevel) = aexp
     end if
     ! End: Sownak
     
     ! For testing rays  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !call ray_analytic_rho(ilevel)  # ALEX-21-02-2016 (otherwise the code does tests and don't want that.)
     
     !when there is no old potential...
     if (nstep==0)call save_phi_old(ilevel)

     ! Compute gravitational acceleration
     call force_fine(ilevel,icount)

     ! Synchronize remaining particles for gravity
     if(pic)then
        call synchro_fine(ilevel)
     end if

     if(hydro)then

        ! Compute Bondi-Hoyle accretion parameters
        if(sink.and.bondi)call bondi_hoyle(ilevel)

        ! Add gravity source term with half time step and new force
        call synchro_hydro_fine(ilevel,+0.5*dtnew(ilevel))

        ! Update boundaries
#ifdef SOLVERmhd
        do ivar=1,nvar+3
#else
        do ivar=1,nvar
#endif
           call make_virtual_fine_dp(uold(1,ivar),ilevel)
#ifdef SOLVERmhd
        end do
#else
        end do
#endif
        if(simple_boundary)call make_boundary_hydro(ilevel)
     end if
  end if

#ifdef RT
  ! Turn on RT in case of rt_stars and first stars just created:
  ! Update photon packages according to star particles
  if(rt .and. rt_star) call update_star_RT_feedback(ilevel)
#endif

  !----------------------
  ! Compute new time step
  !----------------------
  call newdt_fine(ilevel)
  !-----------------------------------
  ! Alex: for_exact_output_times
  ! Start of modification
  !-----------------------------------
  if(ilevel.eq.levelmin .and. cosmo) then
     if(t+dtnew(ilevel).gt.t_next .and. t.le.t_next) then
        dtnew(ilevel) = t_next -t
     end if
  end if
  ! End of modification
  if(ilevel>levelmin)then
     dtnew(ilevel)=MIN(dtnew(ilevel-1)/real(nsubcycle(ilevel-1)),dtnew(ilevel))
  end if

  ! Set unew equal to uold
  if(hydro)call set_unew(ilevel)

#ifdef RT
  ! Set rtunew equal to rtuold
  if(rt)call rt_set_unew(ilevel)
#endif

  !---------------------------
  ! Recursive call to amr_step
  !---------------------------
  if(ilevel<nlevelmax)then
     if(numbtot(1,ilevel+1)>0)then
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

        ! Start: Sownak - ISW (28/10/2015)
!        if(ray.and.ray_do_isw) then
!           aexp_new_ray(ilevel) = aexp
!        end if
        ! End: Sownak
!       write(*,*) 'line3', ilevel, aexp, aexp_new_ray(ilevel), aexp_old_ray(ilevel) !abs_for_sync
        !---------------------------------
        ! Do lensing
        !---------------------------------
        if(ray.and..not.ray_initialized) then

           pi          = 4.0d0*datan(1.0d0)

           coverH0_ray = 299792.458d0/100.d0  !c/100/h [Mpc/h]
           zexp        = 1.0d0/aexp-1.0d0
           chi_wave    = dcom_ray(0.0d0,zexp,omega_m,omega_l,coverH0_ray,1.0d-6)/ray_Lbox! Sownak: 02/09/2015
            
           temp_dist_z = dabs(ray_z_obs)/ray_Lbox+1.0d0
           temp_dist_x = dtan(ray_opening_x/360.0d0*pi)*temp_dist_z
           temp_dist_y = dtan(ray_opening_y/360.0d0*pi)*temp_dist_z
            
           chi_start   = dsqrt(temp_dist_x**2+temp_dist_y**2+temp_dist_z**2) 
           ray_chi_s   = dcom_ray(0.0d0,ray_z_s,omega_m,omega_l,coverH0_ray,1.0d-6)/ray_Lbox

           if(chi_wave.lt.dmin1(chi_start,ray_chi_s)*1.01d0) then
              ilevel2= levelmin
              call init_ray(ilevel2) ! Initialise the rays
           end if
        end if
        if(ray.and.ray_initialized) then                           ! BAOJIU-17-02-2016
           if(ray_multi_out.and.ray_nout.gt.0) then                ! BAOJIU-17-02-2016
              ray_aout = 1.0D0/(1.0D0+ray_zout(ray_iout))          ! BAOJIU-17-02-2016 
              if(ray_aout.lt.aexp.and.ray_aout.gt.aexp_old_ray(ilevel)) then   ! BAOJIU-17-02-2016
                 call ray_step(ilevel,ray_aout)                    ! BAOJIU-17-02-2016
                 ray_iout = ray_iout+1                             ! BAOJIU-17-02-2016
              end if                                               ! BAOJIU-17-02-2016
           end if                                                  ! BAOJIU-17-02-2016
           call ray_step(ilevel,-1.0D0)                            ! BAOJIU-17-02-2016
        end if                                                     ! BAOJIU-17-02-2016
        !---------------------------------
        if(sink)call update_sink(ilevel)
     end if
  else
     call update_time(ilevel)

!    write(*,*) 'line4', ilevel, aexp, aexp_new_ray(ilevel), aexp_old_ray(ilevel) !abs_for_sync
     !---------------------------------
     ! Do lensing
     !---------------------------------
     if(ray.and..not.ray_initialized) then

        pi          = 4.0d0 *datan(1.0d0)

        coverH0_ray = 299792.458d0/100.d0  !c/100/h [Mpc/h]
        zexp        = 1.0d0/aexp-1.0d0
        chi_wave    = dcom_ray(0.0d0,zexp,omega_m,omega_l,coverH0_ray,1.0d-6)/ray_Lbox! Sownak: 02/09/2015

        temp_dist_z = dabs(ray_z_obs)/ray_Lbox+1.0d0
        temp_dist_x = dtan(ray_opening_x/360.0d0*pi)*temp_dist_z
        temp_dist_y = dtan(ray_opening_y/360.0d0*pi)*temp_dist_z
        
        chi_start   = dsqrt(temp_dist_x**2+temp_dist_y**2+temp_dist_z**2)
        ray_chi_s   = dcom_ray(0.0d0,ray_z_s,omega_m,omega_l,coverH0_ray,1.0d-6)/ray_Lbox

        if(chi_wave.lt.dmin1(chi_start,ray_chi_s)*1.01d0) then
           ilevel2= levelmin
           call init_ray(ilevel2) ! Initialise the rays
        end if
     end if
     if(ray.and.ray_initialized) then                              ! BAOJIU-17-02-2016
        if(ray_multi_out.and.ray_nout.gt.0) then                   ! BAOJIU-17-02-2016 
           ray_aout = 1.0D0/(1.0D0+ray_zout(ray_iout))             ! BAOJIU-17-02-2016
           write(*,*) 'xxxxx',aexp,aexp_old_ray(ilevel),ray_aout
           if(ray_aout.lt.aexp.and.ray_aout.gt.aexp_old_ray(ilevel)) then      ! BAOJIU-17-02-2016
              call ray_step(ilevel,ray_aout)                       ! BAOJIU-17-02-2016
              ray_iout = ray_iout+1                                ! BAOJIU-17-02-2016
           end if                                                  ! BAOJIU-17-02-2016
        end if                                                     ! BAOJIU-17-02-2016
        call ray_step(ilevel,-1.0D0)                               ! BAOJIU-17-02-2016
     end if
     !---------------------------------
     if(sink)call update_sink(ilevel)
  end if

  ! Thermal feedback from stars
  if(hydro.and.star.and.eta_sn>0)call thermal_feedback(ilevel)

#ifdef RT
  ! Add stellar radiation sources
  if(rt.and.rt_star) call star_RT_feedback(ilevel,dtnew(ilevel))
#endif

  !---------------
  ! Move particles
  !---------------
  if(pic)then
     call move_fine(ilevel) ! Only remaining particles
  end if

  !-----------
  ! Hydro step
  !-----------
  if(hydro)then

     ! Hyperbolic solver
     call godunov_fine(ilevel)

     ! Reverse update boundaries
#ifdef SOLVERmhd
     do ivar=1,nvar+3
#else
     do ivar=1,nvar
#endif
        call make_virtual_reverse_dp(unew(1,ivar),ilevel)
#ifdef SOLVERmhd
     end do
#else
     end do
#endif
     if(pressure_fix)then
        call make_virtual_reverse_dp(enew(1),ilevel)
        call make_virtual_reverse_dp(divu(1),ilevel)
     endif

     ! Set uold equal to unew
     call set_uold(ilevel)

     ! Density threshold or Bondi accretion onto sink particle
     if(sink)call grow_sink(ilevel)

     ! Add gravity source term with half time step and old force
     ! in order to complete the time step 
     if(poisson)call synchro_hydro_fine(ilevel,+0.5*dtnew(ilevel))

     ! Restriction operator
     call upload_fine(ilevel)

  endif

#ifdef RT
  !---------------
  ! Radiation step
  !---------------
  if(rt)then
     ! Hyperbolic solver
     if(rt_advect) call rt_godunov_fine(ilevel,dtnew(ilevel))

     call add_rt_sources(ilevel,dtnew(ilevel))

     ! Reverse update boundaries
     do ivar=1,nrtvar
        call make_virtual_reverse_dp(rtunew(1,ivar),ilevel)
     end do

     ! Set rtuold equal to rtunew
     call rt_set_uold(ilevel)

     ! Restriction operator
     call rt_upload_fine(ilevel)
  endif
#endif
  
  !-------------------------------
  ! Source term in leaf cells only
  !-------------------------------
  if(neq_chem.or.cooling.or.T2_star>0.0)call cooling_fine(ilevel)

  !----------------------------------
  ! Star formation in leaf cells only
  !----------------------------------
  if(hydro.and.star)call star_formation(ilevel)

  !---------------------------------------
  ! Update physical and virtual boundaries
  !---------------------------------------
  if(hydro)then
#ifdef SOLVERmhd
     do ivar=1,nvar+3
#else
     do ivar=1,nvar
#endif
        call make_virtual_fine_dp(uold(1,ivar),ilevel)
#ifdef SOLVERmhd
     end do
#else
     end do
#endif
     if(simple_boundary)call make_boundary_hydro(ilevel)
  endif
#ifdef RT
  if(rt)then
     do ivar=1,nrtvar
        call make_virtual_fine_dp(rtuold(1,ivar),ilevel)
     end do
     if(simple_boundary)call rt_make_boundary_hydro(ilevel)
  end if
#endif

#ifdef SOLVERmhd
  ! Magnetic diffusion step
 if(hydro)then
     if(eta_mag>0d0.and.ilevel==levelmin)then
        call diffusion
     endif
  end if
#endif

  !-----------------------
  ! Compute refinement map
  !-----------------------
  if(.not.static) call flag_fine(ilevel,icount)

  !----------------------------
  ! Merge finer level particles
  !----------------------------
  if(pic)call merge_tree_fine(ilevel)

  !---------------
  ! Radiation step
  !---------------
#ifdef ATON
  if(aton.and.ilevel==levelmin)then
     call rad_step(dtnew(ilevel))
  endif
#endif

  if(sink)then
     !-------------------------------
     ! Update coarser level sink velocity
     !-------------------------------
     if(ilevel>levelmin)then
        vsold(1:nsink,1:ndim,ilevel-1)=vsnew(1:nsink,1:ndim,ilevel-1)
        if(nsubcycle(ilevel-1)==1)vsnew(1:nsink,1:ndim,ilevel-1)=vsnew(1:nsink,1:ndim,ilevel)
        if(icount==2)vsnew(1:nsink,1:ndim,ilevel-1)= &
             (vsold(1:nsink,1:ndim,ilevel)*dtold(ilevel)+vsnew(1:nsink,1:ndim,ilevel)*dtnew(ilevel))/ &
             (dtold(ilevel)+dtnew(ilevel))
     end if
     !---------------
     ! Sink production
     !---------------
     if(ilevel==levelmin)call create_sink
  end if

  !-------------------------------
  ! Update coarser level time-step
  !-------------------------------
  if(ilevel>levelmin)then
     if(nsubcycle(ilevel-1)==1)dtnew(ilevel-1)=dtnew(ilevel)
     if(icount==2)dtnew(ilevel-1)=dtold(ilevel)+dtnew(ilevel)
  end if

999 format(' Entering amr_step',i1,' for level',i2)

end subroutine amr_step




