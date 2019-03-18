! ------------------------------------------------------------------------
! ray_step
! ------------------------------------------------------------------------
subroutine ray_step(ilevel_arg,aexp_arg)                           ! BAOJIU-17-02-2016: added optional aexp_arg

  use amr_commons
  use pm_commons
  use poisson_commons
  use ray_commons
  use ray_parameters 
  implicit none

#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  integer, intent(in) :: ilevel_arg
  real(dp),intent(in) :: aexp_arg                                  ! BAOJIU-17-02-2016

  integer, dimension(1:3,1:2,1:8) :: iii,jjj

  integer  :: cpu,ind,idim,inbor,icpu
  integer  :: igrid, icell
  integer  :: igshift, igrid_nbor

  integer :: iray, istep

  logical :: ray_myid_finished,ray_allid_finished

  integer :: i

  integer :: ilevel

  real(dp):: tcomm = 0.0D0 

  real(dp) :: chi_target_timestep_myid, chi_target_timestep, chi_target_timestep_allid
  real(dp) :: dray_max, chi_start_timestep, coverH0_ray 
  real(dp) :: E_friedmann_ray, dcom_ray, dcom_integrand_ray

  real(dp) :: zexp, zexp_old, chi_A

  real(dp) :: chi_wave ! SOWNAK-02-09

  ! For taking into account rays moving through corners
  integer ,dimension(1:nvector) :: father_cell
  integer ,dimension(1:nvector,1:threetondim)::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim)::nbors_father_grids

  logical                  :: get_lost 
! logical                  :: additional_comm                      ! BAOJIU-26-09

#ifndef WITHOUTMPI
  integer                  :: comm_counter
  integer                  :: info 
  integer                  :: new_cpu 
  integer, dimension(ncpu) :: nrays_comm
  integer                  :: sum_nrays_comm,nrays_comm_tot 
#endif

  ! For testing
  integer,dimension(1:40000) :: test1
  integer :: jray

  ! For making interpolation (we assume at most 20 fields!!)
  real(dp), dimension(1:RAY_NFIELDS,1:8) :: ray_fields
  real(dp), dimension(1:RAY_NFIELDS) :: tidal_at_center         ! ALEX-21-02-2016
  real(dp) :: phidot_at_center, aexp_mean, zexp_mean            ! ALEX-21-02-2016
  integer :: chunk_size                                         ! BAOJIU-04-10

  !---------------------------------------------------------------------
  ! This routine moves all the rays during one particle time step
  !---------------------------------------------------------------------
  !  1 -> rho
  !  2 -> d2phi_dx2
  !  3 -> d2phi_dy2
  !  4 -> d2phi_dz2
  !  5 -> d2phi_dxdy
  !  6 -> d2phi_dxdz
  !  7 -> d2phi_dydz
  !  8 -> -dphi_dx
  !  9 -> -dphi_dy
  ! 10 -> -dphi_dz
  ! 11 ->  phi
  ! 12 ->  phi_dot (not implemented)

  if(myid==1) write(*,*) 'Entering ray_step'

  iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)
  
  ! ---------------------------- !
  ! ------- BAOJIU-04-10 ------- !
  ! ---------------------------- !

  ray_in_cell = 0
  ray_ncells  = 0

  chunk_size = (ncoarse+twotondim*ngridmax)/50             
  chunk_size = min(chunk_size,ray_nrays)
  chunk_size = max(chunk_size,100)

  allocate(ray_stored_dens(1:chunk_size,1:8))
  allocate(ray_stored_pots(1:chunk_size,1:8))

  if(ray_do_kappa.and.ray_kappa_method.ne.2) then
     allocate(ray_stored_force(1:chunk_size,1:3,1:8))
  end if

  if(ray_do_shear.or.(ray_do_kappa.and.ray_kappa_method.ne.1)) then
     allocate(ray_stored_tidal(1:chunk_size,1:6,1:8))
  end if

  ! +++ add new stuff here to calculate other observables +++

  ! Start: Sownak - ISW (26/10/2015)
  if(ray_do_isw) then
     allocate(ray_stored_phidot(1:chunk_size,1:8))
  end if
  ! End: Sownak

  ! ---------------------------- !
  ! ------- BAOJIU-04-10 ------- !
  ! ---------------------------- !

  coverH0_ray = 299792.458D0/100.D0  !c/100/h [Mpc/h]
  zexp        = 1.0D0/aexp     - 1.0D0
  if(aexp_arg.gt.0.0D0) then                                       ! BAOJIU-17-02-2016
     zexp     = 1.0D0/aexp_arg - 1.0D0                             ! BAOJIU-17-02-2016
  end if                                                           ! BAOJIU-17-02-2016
! zexp_old    = 1.0D0/aexp_old_ray - 1.0D0
! if(ray_z_s.lt.zexp_old.and.ray_z_s.gt.zexp) zexp_old = ray_z_s 
! dray_max    = dcom_ray(zexp ,zexp_old,omega_m,omega_l,coverH0_ray,1.0D-6)/ray_Lbox ! Maximum distance a photon can travel between [aexp_old, aexp]

  chi_wave    = dcom_ray(0.0D0,zexp    ,omega_m,omega_l,coverH0_ray,1.0D-6)/ray_Lbox ! SOWNAK-02-09

  if(ray_nrays.gt.0) then
!    chi_target_timestep_myid = ray_coord(1,1)-dray_max
     chi_target_timestep_myid = dmax1(chi_wave,0.0D0)
  else 
     chi_target_timestep_myid = 0.0D0
  end if

#ifndef WITHOUTMPI
   call MPI_BARRIER(MPI_COMM_WORLD,info)
   call MPI_ALLREDUCE(chi_target_timestep_myid,chi_target_timestep_allid,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
   chi_target_timestep = chi_target_timestep_allid
#else 
   chi_target_timestep = chi_target_timestep_myid
#endif

! additional_comm = .false.                                        ! BAOJIU-26-09

#ifndef WITHOUTMPI

  comm_counter  = 0  

  do while(.true.) ! Communication loop

!    if(.not.additional_comm) nrays_comm(1:ncpu) = 0               ! BAOJIU-26-09
     nrays_comm(1:ncpu) = 0 

#endif

     ray_myid_finished = .true.
     ! Loop over all the rays
!    if(.not.additional_comm) then                                 ! BAOJIU-26-09
     do iray=1,ray_nrays  

        if(ray_coord(iray,1).lt.chi_wave) then
          ray_myid_finished = .false.
          cycle
        end if

        ! Check if the ray already reached the observer
        if(ray_coord(iray,1)*dcos(ray_coord(iray,2))+1.0D0.gt.dabs(ray_z_obs)+1.0D0 .and. ray_coord(iray,1).gt.chi_target_timestep) then 
           get_lost          = .false.
           ray_myid_finished = .false.
        else
           get_lost = .true.
        endif

        ! Identify ray coordinates:
        cpu    = ray_grid(iray,2)/1000
        ilevel = mod(ray_grid(iray,2),100)
        igrid  = ray_grid(iray,1)
        ind    = (ray_grid(iray,2)-cpu*1000)/100
        icell  = ncoarse+(ind-1)*ngridmax+igrid

#ifndef WITHOUTMPI
        if(igrid.lt.0) then                                        ! BAOJIU-06-08
           nrays_comm(cpu) = nrays_comm(cpu)+1                     ! BAOJIU-06-08

           if(comm_counter.eq.0) then
              write(*,*) 'ray_step: grid index negative when ray_step starts; please check.'
              write(*,*) 'debug info:',xg(abs(igrid),1),xg(abs(igrid),2),xg(abs(igrid),3)
              write(*,*) 'debug info:',myid,cpu,ray_id(iray),ray_coord(iray,1)*dcos(ray_coord(iray,2))
              stop
           end if
 
           cycle                                                   ! BAOJIU-06-08
        end if                                                     ! BAOJIU-06-08
#endif
        if(cpu.ne.cpu_map(father(igrid))) then
           write(*,*) 'ray_step: CPU info incorrect; please check CPU-reshuffle-related stuff.'
           write(*,*) 'debug info:',cpu,myid,cpu_map(father(igrid)),ray_id(iray),comm_counter
           stop
        endif

        ! Main time loop:
        !----------------
        istep = 0

        do while(.not.get_lost)

           istep = istep+1 

           ! save ray_coord(iray, 1) at the start of the ray step, to be used in ray_integrate
           chi_A = ray_coord(iray,1) ! save ray_coord(iray, 1) at the start of the ray step

           ! Find new cell and update ray_coord(iray, 1)
           call find_direction(iray,ilevel,ind,idim,inbor)

           ! Check whether inbor == 0 (flagged in find_direction)
           if(inbor.eq.0) then
              get_lost          = .true.
              ray_myid_finished = ray_myid_finished.and..true.
              if (ray_id(iray).eq.ray_Nx*ray_Ny) then
                 write(*,*) 'HERE',myid,istep,comm_counter
                 write(*,*) '123123111',ray_coord(iray,1)*dsin(ray_coord(iray,2))*dcos(ray_coord(iray,3))
                 write(*,*) '123123222',ray_coord(iray,1)*dsin(ray_coord(iray,2))*dsin(ray_coord(iray,3))
                 write(*,*) '123123333',ray_coord(iray,1)*dcos(ray_coord(iray,2)),dabs(ray_z_obs)
                 write(*,*) '123123444',aexp,comm_counter
                 write(*,*) '.........',aexp_old,chi_target_timestep
              end if
           endif

           ! Check whether ray has moved by more than is allowed in timestep -- reduce excess step if so
           if (ray_no_bending) then
!             if (chi_start_timestep-ray_coord(iray,1) .ge. dray_max) then
              if (chi_target_timestep.ge.ray_coord(iray,1)) then
!                ray_coord(iray, 1) = chi_start_timestep - dray_max
                 ray_coord(iray, 1) = chi_target_timestep
                 get_lost          = .true.
                 ray_myid_finished = .false.
                 if (ray_id(iray).eq.ray_Nx*ray_Ny) then
                    write(*,*) '123123111',ray_coord(iray,1)*dsin(ray_coord(iray,2))*dcos(ray_coord(iray,3))
                    write(*,*) '123123222',ray_coord(iray,1)*dsin(ray_coord(iray,2))*dsin(ray_coord(iray,3))
                    write(*,*) '123123333',ray_coord(iray,1)*dcos(ray_coord(iray,2))
                    write(*,*) '123123444',aexp,comm_counter
                    write(*,*) '.........',aexp_old,chi_target_timestep
                 end if
              end if
           else
           ! Ray not moving towards the observer: need to check distance travelled in non-radial directions
           end if

           ! ALEX-21-02-2016 Start of modification (added a series of if statements to bifurcate corners and NGP integration)
           if(ray_do_ngp) then
              ! Use quantities at cell centers
              call tidal_tensor(ilevel, igrid, icell, ind, tidal_at_center)
              phidot_at_center = 0.0d0
              if(ray_do_isw) then
                 if(ilevel.gt.levelmin .and. dabs(phi_old(icell)/phi_old(father(igrid))-1.0D0).lt.1.0D-10) then
                    ! If the cell had just been refined, then compute time derivative using coarse level for stability
                    aexp_mean = (aexp_new_ray(ilevel-1)+aexp_old_ray(ilevel-1))/2.0D0
                    zexp_mean = 1.0D0/aexp_mean - 1.0D0
                    phidot_at_center= (phi(father(igrid))/aexp_new_ray(ilevel-1)**2 - phi_old(father(igrid))/aexp_old_ray(ilevel-1)**2) & 
                         / (aexp_new_ray(ilevel-1)-aexp_old_ray(ilevel-1)) * E_friedmann_ray(zexp_mean,omega_m,omega_l) * aexp_mean
                 else
                    ! If the cell already existed in previous steps, use current cells data
                    aexp_mean = (aexp_new_ray(ilevel)+aexp_old_ray(ilevel))/2.0D0
                    zexp_mean = 1.0D0/aexp_mean - 1.0D0
                    phidot_at_center = (phi(icell)/aexp_new_ray(ilevel)**2&
                                       - phi_old(icell)/aexp_old_ray(ilevel)**2)&
                                       / (aexp_new_ray(ilevel) &
                                          - aexp_old_ray(ilevel)) &
                                       * E_friedmann_ray(zexp_mean,omega_m,omega_l) &
                                       * aexp_mean
                 end if
              end if
           else
              ! Interpolate given quantity from cell centers to cell corners
              call ray_fields_in_corners(ilevel, igrid, ind, ray_fields)
           end if

           ! Integrate given quantity in cell
           if(ray_do_kappa) then                                   ! BAOJIU-29-09
              if(ray_kappa_method.ne.2) then                       ! BAOJIU-29-09
                 if(ray_do_ngp) then
                    write(*,*) 'NGP integration method not available for &
                                &ray_kappa_method = 1 (Method A for kappa)'
                    write(*,*) 'Switch to ray_kappa_method = 2'
                    stop 
                 else
                    ! BAOJIU-29-09; this adds to ray_kappa(iray,1)
                    call ray_integrate(iray, ilevel, ind, chi_A, ray_fields, 1)
                 end if
              end if                                               ! BAOJIU-29-09
              if(ray_kappa_method.ne.1) then                       ! BAOJIU-29-09
                 if(ray_do_ngp) then
                    call ray_integrate_ngp(iray, ilevel, ind, chi_A, &
                                           tidal_at_center, phidot_at_center, 2)
                 else
                    ! BAOJIU-29-09; this adds to ray_kappa(iray,2)
                    call ray_integrate(iray, ilevel, ind, chi_A, &
                                       ray_fields, 2)
                 end if
              end if                                               ! BAOJIU-29-09
           end if                                                  ! BAOJIU-29-09

           if(ray_do_shear) then                                   ! BAOJIU-29-09
              if(ray_do_ngp) then
                 call ray_integrate_ngp(iray, ilevel, ind, chi_A, &
                                        tidal_at_center, phidot_at_center, 3)
                 call ray_integrate_ngp(iray, ilevel, ind, chi_A, &
                                        tidal_at_center, phidot_at_center, 4)
              else
                 ! BAOJIU-29-09; this adds to ray_shear(iray,1)
                 call ray_integrate(iray, ilevel, ind, chi_A, &
                                    ray_fields, 3)
                 ! BAOJIU-29-09; this adds to ray_shear(iray,2)
                 call ray_integrate(iray, ilevel, ind, chi_A, &
                                    ray_fields, 4)
              end if
           end if                                                  ! BAOJIU-29-09
           
           ! Start: Sownak - ISW (26/10/2015)
           if (ray_do_isw) then
              if(ray_do_ngp) then
                 call ray_integrate_ngp(iray, ilevel, ind, chi_A, &
                                        tidal_at_center, phidot_at_center, 5)
              else
                 call ray_integrate(iray, ilevel, ind, chi_A, &
                                    ray_fields, 5)
              end if
           end if

           ! +++ add new stuff here to calculate other observables +++
           
           if(ray_do_deflection) then                                 ! Christoph
              if(ray_do_ngp) then
                  ! FIXME needs to be implemented
                 call ray_integrate_ngp(iray, ilevel, ind, chi_A, &
                                        tidal_at_center, phidot_at_center, 6)
                 call ray_integrate_ngp(iray, ilevel, ind, chi_A, &
                                        tidal_at_center, phidot_at_center, 7)
              else
                 call ray_integrate(iray, ilevel, ind, chi_A, &
                                    ray_fields, 6)
                 call ray_integrate(iray, ilevel, ind, chi_A, &
                                    ray_fields, 7)
              end if
           end if

           ! End: Sownak
           ! ALEX-21-02-2016 End of modification
           ! (added a series of if statements to 
           !  bifurcate corners and NGP integration)

           ! If get_lost set to true above then cycle this ray 
           ! This can either mean end of total integration or end of particle time step
           if(get_lost) exit 

           ! Get neighbour grid
           if(inbor.gt.0) then                                     ! normal case
              igshift = iii(idim,inbor,ind)
              if(igshift==0) then
                 igrid_nbor = igrid
              else
                 igrid_nbor = son(nbor(igrid,igshift))
              end if
           else                                                    ! abnormal case (ray going out from border)
              father_cell(1) = father(igrid)
              call get3cubefather(father_cell,nbors_father_cells,nbors_father_grids,1,ilevel)
              igshift    = ray_kkk(-inbor,ind)
              igrid_nbor = son(nbors_father_cells(1,igshift))
              ! CLAUDIO: CAN WE TAKE THIS AWAY???
              if(igrid_nbor.ne.0) ray_grid(iray,1) = igrid_nbor    ! BAOJIU-06-08
           endif
 
           ! Check if neighbour does NOT exist:
           !-----------------------------------
           if(igrid_nbor.eq.0) then
              ! Find coarse cell
              icell  = father(ray_grid(iray,1))
              ! Find coarse grid
              ind    = (icell-ncoarse-1)/ngridmax+1
              igrid  = icell-ncoarse-(ind-1)*ngridmax

              ! Update ilevel
              ilevel = ilevel-1

              ! Find neighbour grid in the new level
              if(inbor.gt.0) then                                  ! straight neighbour case
                 igshift = iii(idim,inbor,ind)
                 if(igshift==0) then
                    igrid_nbor = igrid
                 else
                    igrid_nbor = son(nbor(igrid,igshift))
                 end if
              else                                                 ! diagonal neighbour case
                 father_cell(1) = father(igrid)
                 call get3cubefather(father_cell,nbors_father_cells,nbors_father_grids,1,ilevel)
                 igshift    = ray_kkk(-inbor,ind)
                 igrid_nbor = son(nbors_father_cells(1,igshift))
                 ! CLAUDIO: CAN WE TAKE THIS AWAY???
                 ray_grid(iray,1) = igrid_nbor                     ! BAOJIU-06-08
              endif

           endif

           ! Store new grid in ray_grid and find new cell:
           !-----------------------------------------------
           igrid = igrid_nbor
           ray_grid(iray,1) = igrid_nbor

           if(inbor.gt.0) then                                     ! straight neighbour case
              icell = igrid+(ncoarse+(    jjj(idim,inbor,ind)-1)*ngridmax)
              ind   = jjj(idim,inbor,ind)
           else                                                    ! diagonal neighbour case
              icell = igrid+(ncoarse+(ray_lll(    -inbor,ind)-1)*ngridmax)
              ind   = ray_lll(-inbor,ind)
           endif

           ! Check to see if we need to communicate (if not then find ref):
           !---------------------------------------------------------------
           new_cpu = cpu_map(father(igrid)) 

           if(new_cpu.eq.myid) then
              ! CLAUDIO: CAN WE TAKE THIS AWAY???
              ray_grid(iray,1) = igrid                             ! BAOJIU-06-08

              ! Find refinements in new cell:
              !------------------------------
              if(son(icell).ne.0) then                             ! BAOJIU-06-08 

                 ! Find finer grid
                 igrid = son(icell)
                 ray_grid(iray,1) = igrid

                 ! Find finer cell (this requires the old ilevel) - requires WHAT?
                 call get_ind_from_ray_position(iray,ilevel+1,ray_grid(iray,1),ind) ! ALEX-21-09  
                 icell = ncoarse+(ind-1)*ngridmax+ray_grid(iray,1)

                 ilevel = ilevel+1

                 ! Check if finer grid belongs to myid             ! BAOJIU-06-08 
                 new_cpu = cpu_map(father(igrid))                  ! BAOJIU-06-08
                 if(new_cpu.ne.myid) then                          ! BAOJIU-06-08
                    ray_grid(iray,1) = -igrid                      ! BAOJIU-06-08
                    nrays_comm(new_cpu) = nrays_comm(new_cpu)+1    ! BAOJIU-06-08
                    get_lost = .true.                              ! BAOJIU-06-08
                 end if                                            ! BAOJIU-06-08
 
              endif

           else ! When the ray is NOT in myid
              ray_grid(iray,1) = -igrid 
              nrays_comm(new_cpu) = nrays_comm(new_cpu)+1 
              get_lost = .true. 
           end if 

           ! Store new cell information in ray_grid
           ray_grid(iray,2) = 1000*new_cpu+100*ind+ilevel 

        end do  ! Loop for time steps (do while)

     end do  ! Loop of the rays (iray)
!    end if                                                        ! BAOJIU-26-09

     ! Take care of parallelisation:
     !------------------------------
#ifndef WITHOUTMPI
     ! Wait for all the rays to reach the end point.
     call MPI_BARRIER(MPI_COMM_WORLD,info)

     ! Calculate total (global) number of rays to be communicated
     sum_nrays_comm = SUM(nrays_comm)
     call MPI_ALLREDUCE(sum_nrays_comm,nrays_comm_tot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)

     ! If we need to communicate ways...
     if (nrays_comm_tot.gt.0) then 
!       write(*,666) nrays_comm_tot  
        do icpu=1,ncpu
           if(icpu.eq.myid) cycle
           if(nrays_comm(icpu).ne.0) then
              write(*,777) myid,icpu,nrays_comm(icpu) 
           end if
        end do
        call ray_communicate(ilevel_arg,1,nrays_comm)
        comm_counter = comm_counter+1                              ! BAOJIU-06-08
!       do iray=1,ray_nrays
!          if(myid.ne.ray_grid(iray,2)/1000) &
!               & write(*,888) comm_counter,myid,ray_id(iray),ray_grid(iray,2)/1000,ray_coord(iray,1),ray_coord(iray,2),ray_coord(iray,3) 
!       end do
     else
        exit ! Nothing else to do; exit communication loop
     endif

     ! Check for refinements for all the rays:
!    nrays_comm(1:ncpu) = 0                                        ! BAOJIU-26-09
!    additional_comm = .false.                                     ! BAOJIU-26-09
     do iray=1,ray_nrays
        cpu = ray_grid(iray,2)/1000
        if(cpu.ne.myid) then
           write(*,*) 'ray_step: communication not properly done; please check.'
           stop
        end if
        ilevel = mod(ray_grid(iray,2),100)
        igrid  = ray_grid(iray,1)
        ind    = (ray_grid(iray,2)-cpu*1000)/100
        icell  = ncoarse+(ind-1)*ngridmax+igrid
 
        ! Check if cell is refined
        if(son(icell).ne.0) then 
 
           ! Find finer grid
           ray_grid(iray,1) = son(icell)
 
           ! Find finer cell (this requires the old ilevel)
           call get_ind_from_ray_position(iray,ilevel+1,ray_grid(iray,1),ind) ! ALEX-21-09   

           ! Check if the finer cell belongs to another CPU 
           if(cpu_map(father(ray_grid(iray,1))).ne.cpu) then       ! BAOJIU-06-08
              cpu = cpu_map(father(ray_grid(iray,1)))              ! BAOJIU-06-08
              ray_grid(iray,1) = - ray_grid(iray,1)                ! BAOJIU-06-08
!             nrays_comm(cpu)  = nrays_comm(cpu)+1                 ! BAOJIU-26-09
!             additional_comm = .true.                             ! BAOJIU-26-09
           end if
           ! Store updated ray_grid
           ray_grid(iray,2) = 1000*cpu+100*ind+(ilevel+1)
 
        endif
     enddo

  end do ! Loop over communications

#endif

  ! ---------------------------- !
  ! ------- BAOJIU-04-10 ------- !
  ! ---------------------------- !

  if(allocated(ray_stored_dens) ) deallocate(ray_stored_dens)
  if(allocated(ray_stored_pots) ) deallocate(ray_stored_pots)
  if(allocated(ray_stored_force)) deallocate(ray_stored_force)
  if(allocated(ray_stored_tidal)) deallocate(ray_stored_tidal)

  ! +++ add new stuff here to calculate other observables +++

  ! Start: Sownak - ISW (26/10/2015)
  if(allocated(ray_stored_phidot)) deallocate(ray_stored_phidot)
  ! End: Sownak
  
  ! ---------------------------- !
  ! ------- BAOJIU-04-10 ------- !
  ! ---------------------------- !

#ifndef WITHOUTMPI
  call MPI_BARRIER(MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(ray_myid_finished,ray_allid_finished,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,info)
  ray_all_finished = ray_allid_finished
#else
  ray_all_finished = ray_myid_finished
#endif

! if(ray_nrays.gt.0) write(*,*) myid,ray_nrays,ray_myid_finished,ray_allid_finished,aexp

  ! ------------- START ------------- !
  ! ------- BAOJIU-17-02-2016 ------- !
  ! ------------- START ------------- !

  if(aexp_arg.gt.0.0D0.and.ray_nrays.gt.0) then
     write(*,*) 'xxxxx',ray_kappa(ray_nrays,1),ray_shear(ray_nrays,1),ray_phidot(ray_nrays)
     call ray_write_final_result(ray_iout)                     
     if(ray_afresh) then                                      
        if(ray_do_kappa) ray_kappa  = 0.0D0                 
        if(ray_do_shear) ray_shear  = 0.0D0                
        if(ray_do_isw)   ray_phidot = 0.0D0        
     end if
  end if

  ! -------------  END  ------------- !
  ! ------- BAOJIU-17-02-2016 ------- !
  ! -------------- END -------------- !

  if(ray_all_finished.and.ray_nrays.gt.0) then

!    write(*,*) "All rays finished."
!    do jray=1,ray_nrays
!       do iray=1,ray_nrays
!          if(ray_id(iray).eq.jray) then
!             test1(jray)=iray
!             exit 
!          end if
!       end do
!    end do
!    open(unit = 1231232, file = 'coord_test.dat')
!    open(unit = 1231234, file = 'kappa_test.dat')
!    open(unit = 1231238, file = 'id_test.dat')
!    do iray=1,ray_nrays
!       write(1231232, 555) ray_id(test1(iray)), ray_coord(test1(iray),1),ray_coord(test1(iray),2),ray_coord(test1(iray),3),ray_kappa(test1(iray))
!       write(1231234, *) ray_kappa(test1(iray))
!       write(1231238, *) ray_id(test1(iray))
!    end do
!    close(unit = 1231232)
!    close(unit = 1231234)
!    close(unit = 1231238)
     
     call ray_write_final_result(ray_nout+1)                       ! BAOJIU-17-02-2016
  end if

  if((ray_end_of_sim.or.ray_stop_after_raytracing).and.ray_all_finished) then
     write(*,*) "Ray integration finished. Run stops here."
     call clean_stop
  end if

555 format(I6,F10.6,F18.6,F18.6,F18.6)
666 format('All CPUs communicating',I10,'  rays in total.')
777 format('CPU',I5,' is sending CPU',I5,' a number of',I10,' rays')
888 format('During communication round',I3,', myid',I5,' will send ray ID',I10,' to CPU',I5,' and coord (',3(F10.6,' '),')')
999 format('After communication round',I3,', myid',I5,' contains ray ID',I10,' with coord (',3(F10.6,' '),')')

end subroutine ray_step

!===============================
! ray_write_final_result
!===============================
subroutine ray_write_final_result(iout_arg)                        ! BAOJIU-17-02-2016: added aexp_arg
  use amr_commons
  use poisson_commons
  use ray_commons
  use pm_commons
  implicit none

  integer,intent(in) :: iout_arg                                   ! BAOJIU-17-02-2016

  ! Writes final result
  
  integer :: iray
  integer :: n_quants                                              ! BAOJIU-29-09

  character*1000 file_name
  character(LEN=5)::nchar,nchar2                                   ! BAOJIU-17-02-2016
  
  !if(myid==1) then
     write(*,*) "The output format for lensing is:",myid
     write(*,*) "ray_id"
     write(*,*) "ray_coord 1"
     write(*,*) "ray_coord 2"
     write(*,*) "ray_coord 3"
     if(ray_do_kappa) then
        write(*,*) "ray_kappa 1"
     endif
     
     if(ray_do_kappa .and. ray_kappa_method.eq.3) then
        write(*,*) "ray_kappa 2"
     endif
     
     if(ray_do_shear) then
        write(*,*) "ray_shear 1"
        write(*,*) "ray_shear 2"
     endif
     
     if(ray_do_isw) then
        write(*,*) "ray_phidot"
     endif
  !endif

  call title(myid,nchar)
  call title(iout_arg,nchar2)                                      ! BAOJIU-17-02-2016
  file_name = 'Ray_maps_o'//trim(nchar2)//'_'//trim(nchar)//'.dat' ! BAOJIU-17-02-2016
  open(unit=100, file=file_name, form='formatted')

  do iray=1,ray_nrays

     write(100, "(I10,E20.8,E20.8,E20.8)", advance='no') &
          ray_id(iray), &
          ray_coord(iray,1), &
          ray_coord(iray,2), &
          ray_coord(iray,3)
     
     if(ray_do_kappa) then
        write(100, "(E20.10)", advance='no') ray_kappa(iray,1)
     endif
     
     if(ray_do_kappa .and. ray_kappa_method.eq.3) then
        write(100, "(E20.10)", advance='no') ray_kappa(iray,2)
     endif
     
     if(ray_do_shear) then
        write(100, "(E20.10,E20.10)", advance='no') ray_shear(iray,1), ray_shear(iray,2)
     endif
     
     if(ray_do_isw) then
        write(100, "(E20.10)", advance='no') ray_phidot(iray)
     endif
     
     write(100, "(A1)") ""  ! this is to get a carriage return
     
  enddo
  
  close(100)


end subroutine ray_write_final_result


!===============================
! ray_write_final_result_old
!===============================
subroutine ray_write_final_result_old()
  use amr_commons
  use poisson_commons
  use ray_commons
  use pm_commons
  implicit none

  ! Writes final result
  
  integer :: iray
  integer :: n_quants                                              ! BAOJIU-29-09

  character*1000 file_name
  character(LEN=5)::nchar
  
  call title(myid,nchar)
  file_name = 'ray_maps_'//trim(nchar)//'.dat'
  open(unit=100, file=file_name, form='formatted')

  
  n_quants = 0                                                     ! BAOJIU-29-09
  if(ray_do_kappa) then                                            ! BAOJIU-29-09
     n_quants = n_quants+1                                         ! BAOJIU-29-09
     if(ray_kappa_method.eq.3) n_quants = n_quants+1               ! BAOJIU-29-09
  end if                                                           ! BAOJIU-29-09
  if(ray_do_shear) then                                            ! BAOJIU-29-09
     n_quants = n_quants+2                                         ! BAOJIU-29-09
  end if                                                           ! BAOJIU-29-09

  ! +++ add new stuff here to calculate other observables +++

  ! Start: Sownak - ISW (26/10/2015)
  if(ray_do_isw) then
     n_quants = n_quants+1
  end if
  ! End: Sownak

  if(n_quants.eq.1) then                                           ! BAOJIU-29-09
     do iray=1,ray_nrays                                           ! BAOJIU-29-09
        write(100,1111) &                                          ! BAOJIU-29-09
             ray_id(iray), &                                       ! BAOJIU-29-09
             ray_coord(iray,1), &                                  ! BAOJIU-29-09
             ray_coord(iray,2), &                                  ! BAOJIU-29-09
             ray_coord(iray,3), &                                  ! BAOJIU-29-09
             ray_kappa(iray,1)                                     ! BAOJIU-29-09
     enddo                                                         ! BAOJIU-29-09
  end if                                                           ! BAOJIU-29-09

  if(n_quants.eq.2.and.ray_do_kappa) then                          ! BAOJIU-29-09
     do iray=1,ray_nrays                                           ! BAOJIU-29-09
        write(100,2222) &                                          ! BAOJIU-29-09
             ray_id(iray), &                                       ! BAOJIU-29-09
             ray_coord(iray,1), &                                  ! BAOJIU-29-09
             ray_coord(iray,2), &                                  ! BAOJIU-29-09
             ray_coord(iray,3), &                                  ! BAOJIU-29-09 
             ray_kappa(iray,1), &                                  ! BAOJIU-29-09
             ray_kappa(iray,2)                                     ! BAOJIU-29-09
     enddo                                                         ! BAOJIU-29-09
  end if                                                           ! BAOJIU-29-09

  if(n_quants.eq.2.and.ray_do_shear) then                          ! BAOJIU-29-09
     do iray=1,ray_nrays                                           ! BAOJIU-29-09
        write(100,2222) &                                          ! BAOJIU-29-09 
             ray_id(iray), &                                       ! BAOJIU-29-09
             ray_coord(iray,1), &                                  ! BAOJIU-29-09
             ray_coord(iray,2), &                                  ! BAOJIU-29-09
             ray_coord(iray,3), &                                  ! BAOJIU-29-09
             ray_shear(iray,1), &                                  ! BAOJIU-29-09
             ray_shear(iray,2)                                     ! BAOJIU-29-09
     enddo                                                         ! BAOJIU-29-09
  end if                                                           ! BAOJIU-29-09

  ! Start: Sownak - ISW (26/10/2015)
  if (n_quants.eq.2.and.ray_do_isw) then
     do iray=1,ray_nrays
        write(100,2222) & 
             ray_id(iray), &
             ray_coord(iray,1), &
             ray_coord(iray,2), &
             ray_coord(iray,3), &
             ray_kappa(iray,1), &
             ray_phidot(iray)
     end do
  end if
  ! End: Sownak

  if(n_quants.eq.3) then                                           ! BAOJIU-29-09
     do iray=1,ray_nrays                                           ! BAOJIU-29-09 
        write(100,3333) &                                          ! BAOJIU-29-09 
             ray_id(iray), &                                       ! BAOJIU-29-09 
             ray_coord(iray,1), &                                  ! BAOJIU-29-09 
             ray_coord(iray,2), &                                  ! BAOJIU-29-09 
             ray_coord(iray,3), &                                  ! BAOJIU-29-09
             ray_kappa(iray,1), &                                  ! BAOJIU-29-09
             ray_shear(iray,1), &                                  ! BAOJIU-29-09
             ray_shear(iray,2)                                     ! BAOJIU-29-09
     enddo                                                         ! BAOJIU-29-09
  end if                                                           ! BAOJIU-29-09

  ! Start: Sownak - ISW (26/10/2015)
  if (n_quants.eq.3.and.ray_do_isw) then
     do iray=1,ray_nrays
        write(100,3333) &
             ray_id(iray), &
             ray_coord(iray,1), &
             ray_coord(iray,2), &
             ray_coord(iray,3), &
             ray_kappa(iray,1), &
             ray_kappa(iray,2), &
             ray_phidot(iray)
     end do
  end if
  ! End: Sownak

  if(n_quants.eq.4) then                                           ! BAOJIU-29-09
     do iray=1,ray_nrays                                           ! BAOJIU-29-09
        write(100,4444) &                                          ! BAOJIU-29-09
             ray_id(iray), &                                       ! BAOJIU-29-09
             ray_coord(iray,1), &                                  ! BAOJIU-29-09
             ray_coord(iray,2), &                                  ! BAOJIU-29-09
             ray_coord(iray,3), &                                  ! BAOJIU-29-09 
             ray_kappa(iray,1), &                                  ! BAOJIU-29-09
             ray_kappa(iray,2), &                                  ! BAOJIU-29-09
             ray_shear(iray,1), &                                  ! BAOJIU-29-09
             ray_shear(iray,2)                                     ! BAOJIU-29-09
     enddo                                                         ! BAOJIU-29-09 
  end if                                                           ! BAOJIU-29-09

  ! Start: Sownak - ISW (26/10/2015)
  if (n_quants.eq.4.and.ray_do_isw) then
     do iray=1,ray_nrays
        write(100,4444) &
             ray_id(iray), &
             ray_coord(iray,1), &
             ray_coord(iray,2), &
             ray_coord(iray,3), &
             ray_kappa(iray,1), &
             ray_shear(iray,1), &
             ray_shear(iray,2), &
             ray_phidot(iray)
     end do
  end if
  ! End: Sownak

  ! Start: Sownak - ISW (26/10/2015)
  if (n_quants.eq.5.and.ray_do_isw) then
     do iray=1,ray_nrays
        write(100,5555) &
             ray_id(iray), &
             ray_coord(iray,1), &
             ray_coord(iray,2), &
             ray_coord(iray,3), &
             ray_kappa(iray,1), &
             ray_kappa(iray,2), &
             ray_shear(iray,1), &
             ray_shear(iray,2), &
             ray_phidot(iray)
     end do
  end if
  ! End: Sownak


  ! +++ add new stuff here to calculate other observables +++

  close(100)

1111 format(I10,E20.12,E20.12,E20.12,E20.12)                        ! BAOJIU-29-09
2222 format(I10,E20.12,E20.12,E20.12,E20.12,E20.12)                 ! BAOJIU-29-09
3333 format(I10,E20.12,E20.12,E20.12,E20.12,E20.12,E20.12)          ! BAOJIU-29-09
4444 format(I10,E20.12,E20.12,E20.12,E20.12,E20.12,E20.12,E20.12)   ! BAOJIU-29-09
! Start: Sownak - ISW (26/10/2015)
5555 format(I10,E20.12,E20.12,E20.12,E20.12,E20.12,E20.12,E20.12,E20.12)
! End: Sownak

end subroutine ray_write_final_result_old

!===============================
! ray_write_cell
!===============================
subroutine ray_write_cell(ilevel,iray,istep)
  use amr_commons
  use poisson_commons
  use ray_commons
  use pm_commons
  implicit none

  ! For testing porpoises
  ! Arguments:
  !            iray: to identify grid
  !            ind:  to identify cell
  !            istep: number of step

  integer :: ilevel, iray, istep

  integer :: flag

  integer :: ix, iy, iz
  real(dp),dimension(1:3)::xc  ! Center of the cells
  real(dp) :: dx

  integer :: i_ind
  integer :: cpu, igrid, ind, icell

  cpu    = ray_grid(iray,2)/1000  
  igrid  = ray_grid(iray,1)
  ind    = (ray_grid(iray,2)-cpu*1000)/100
  icell  = ncoarse+(ind-1)*ngridmax+igrid
        
  dx  = 0.5D0**ilevel

  !do i_ind=1,twotondim
  iz=(ind-1)/4
  iy=(ind-1-4*iz)/2
  ix=(ind-1-2*iy-4*iz)
  xc(1)=(dble(ix)-0.5D0)*dx
  xc(2)=(dble(iy)-0.5D0)*dx
  xc(3)=(dble(iz)-0.5D0)*dx
  write(1000,111) &
       ray_id(iray), &
       xg(igrid,1) + xc(1), &
       xg(igrid,2) + xc(2), &
       xg(igrid,3) + xc(3), &
       ray_coord(iray,1)*dsin(ray_coord(iray,2))*dcos(ray_coord(iray,3)), &
       ray_coord(iray,1)*dsin(ray_coord(iray,2))*dsin(ray_coord(iray,3)), &
       ray_coord(iray,1)*dcos(ray_coord(iray,2)), &
       ilevel, &
       rho(icell)
  
  !enddo
  
  !dsin(ray_coord(iray,2))*dcos(ray_coord(iray,3)), &
  !     dsin(ray_coord(iray,2))*dsin(ray_coord(iray,3)), &
  !     dcos(ray_coord(iray,2)), &
       
  !ray_coord(iray,1)*dsin(ray_coord(iray,2))*dcos(ray_coord(iray,3)), &
  !     ray_coord(iray,1)*dsin(ray_coord(iray,2))*dsin(ray_coord(iray,3)), &
  !     ray_coord(iray,1)*dcos(ray_coord(iray,2)), &

  call flush(1000)

111 format(I4, E20.8, E20.8, E20.8, E20.8, E20.8, E20.8, I4, E20.8)


end subroutine ray_write_cell

!===============================
! ray_write_cell
!===============================
subroutine ray_write_corners(ilevel, iray, icell, ind, istep, ray_fields)
  use amr_commons
  use poisson_commons
  use ray_commons
  use pm_commons
  implicit none

  ! For testing porpoises
  ! Arguments:
  !            iray: to identify grid
  !            ind:  to identify cell
  !            istep: number of step

  integer :: ilevel, iray, icell, ind, istep
  real (dp), dimension(1:RAY_NFIELDS,1:8) :: ray_fields

  integer :: flag

  integer :: ix, iy, iz
  real(dp),dimension(1:3)::xc  ! Center of the cells
  real(dp) :: dx

  integer :: cpu

  !if(iray .ne. 214) return

  cpu = ray_grid(iray,2)/1000  

  dx  = 0.5D0**ilevel
  iz=(ind-1)/4
  iy=(ind-1-4*iz)/2
  ix=(ind-1-2*iy-4*iz)
  xc(1)=(dble(ix)-0.5D0)*dx
  xc(2)=(dble(iy)-0.5D0)*dx
  xc(3)=(dble(iz)-0.5D0)*dx
  xc(1) = xg(abs(ray_grid(iray,1)),1) + xc(1)
  xc(2) = xg(abs(ray_grid(iray,1)),2) + xc(2)
  xc(3) = xg(abs(ray_grid(iray,1)),3) + xc(3)
  !rho(icell) = dble(icell)
  write(1000,222) &
       iray, xc(1), xc(2), xc(3), &
       xc(1)-dx/2.0D0, xc(2)-dx/2.0D0, xc(3)-dx/2.0D0, &
       rho(icell), ray_fields(1,1), 1

  write(1000,222) &
       iray, xc(1), xc(2), xc(3), &
       xc(1)+dx/2.0D0, xc(2)-dx/2.0D0, xc(3)-dx/2.0D0, &
       rho(icell), ray_fields(1,2), 2

  write(1000,222) &
       iray, xc(1), xc(2), xc(3), &
       xc(1)-dx/2.0D0, xc(2)+dx/2.0D0, xc(3)-dx/2.0D0, &
       rho(icell), ray_fields(1,3), 3

  write(1000,222) &
       iray, xc(1), xc(2), xc(3), &
       xc(1)+dx/2.0D0, xc(2)+dx/2.0D0, xc(3)-dx/2.0D0, &
       rho(icell), ray_fields(1,4), 4

  write(1000,222) &
       iray, xc(1), xc(2), xc(3), &
       xc(1)-dx/2.0D0, xc(2)-dx/2.0D0, xc(3)+dx/2.0D0, &
       rho(icell), ray_fields(1,5), 5

  write(1000,222) &
       iray, xc(1), xc(2), xc(3), &
       xc(1)+dx/2.0D0, xc(2)-dx/2.0D0, xc(3)+dx/2.0D0, &
       rho(icell), ray_fields(1,6), 6

  write(1000,222) &
       iray, xc(1), xc(2), xc(3), &
       xc(1)-dx/2.0D0, xc(2)+dx/2.0D0, xc(3)+dx/2.0D0, &
       rho(icell), ray_fields(1,7), 7

  write(1000,222) &
       iray, xc(1), xc(2), xc(3), &
       xc(1)+dx/2.0D0, xc(2)+dx/2.0D0, xc(3)+dx/2.0D0, &
       rho(icell), ray_fields(1,8), 8



  call flush(1000)

222 format(I4, E20.8, E20.8, E20.8, E20.8, E20.8, E20.8, E20.8, E20.8, I4)


end subroutine ray_write_corners


!===============================
! ray_write_grid
!===============================
subroutine ray_write_grid(ilevel, iray, ind, istep, flag)
  use amr_commons
  use poisson_commons
  use ray_commons
  use pm_commons
  implicit none

  ! For testing porpoises
  ! Arguments:
  !            iray: to identify grid
  !            ind:  to identify cell
  !            istep: number of step

  integer :: ilevel, iray, ind, istep
  integer :: flag

  integer :: ix, iy, iz
  real(dp),dimension(1:3)::xc  ! Center of the cells
  real(dp) :: dx

  integer :: i_ind

  !if(iray .ne. 214) return

  write(1000,888) &
       iray, &
       istep, &
       xg(abs(ray_grid(iray,1)),1), &
       xg(abs(ray_grid(iray,1)),2), &
       xg(abs(ray_grid(iray,1)),3), &
       ray_coord(iray,1)*dsin(ray_coord(iray,2))*dcos(ray_coord(iray,3)), &
       ray_coord(iray,1)*dsin(ray_coord(iray,2))*dsin(ray_coord(iray,3)), &
       ray_coord(iray,1)*dcos(ray_coord(iray,2)), &
       ilevel, &
       xp(1,3), &
       ray_coord(iray, 1), &
       flag
  !enddo

  call flush(1000)

888 format(I4, I4, E20.8, E20.8, E20.8, E20.8, E20.8, E20.8, I4, E20.8, E20.8, I4)


end subroutine ray_write_grid


!#############################################################
! Copied from the codes Alex used in the Cluster lensing paper
! Check the latter for generalizing to non LCDM backgrounds
!#############################################################

!---------------------------
! This is H/H_0 (physical H)
!---------------------------
function E_friedmann_ray(z,Om0,Ode0)
  implicit none
  integer,parameter :: dl = KIND(1.D0)
  real(dl) :: E_friedmann_ray
  real(dl), intent(in) :: z,Om0,Ode0
  E_friedmann_ray = sqrt(Om0*(1.0D0+z)**3+Ode0)
end function E_friedmann_ray

!------------------------------------
! Integrand for the comoving distance 
!------------------------------------
function dcom_integrand_ray(z,Om0,Ode0)
  implicit none
  integer,parameter :: dl = KIND(1.D0)
  real(dl) :: dcom_integrand_ray,E_friedmann_ray
  real(dl),intent(in) :: z,Om0,Ode0
  dcom_integrand_ray = 1.0D0/E_friedmann_ray(z,Om0,Ode0)
end function dcom_integrand_ray

!-----------------------------------------
! Comoving distance from z1 to z2 in Mpc/h
!-----------------------------------------
function dcom_ray(z1,z2,Om0,Ode0,coverH0,tol_dist)
  implicit none
  integer,parameter :: dl = KIND(1.D0)
  real(dl) :: dcom_ray, rombint_dcom_integrand_ray,dcom_integrand_ray
  real(dl), intent(in) :: z1,z2,Om0,Ode0,coverH0,tol_dist
  dcom_ray = rombint_dcom_integrand_ray(z1,z2,tol_dist,Om0,Ode0)
  dcom_ray = coverH0*dcom_ray
end function dcom_ray

!----------------------------------------------------------------------
! rombint integrator (based on that in subroutines.f90 in the CAMB code)
!----------------------------------------------------------------------
function rombint_dcom_integrand_ray(a,b,tol,Om0,Ode0)
  !        use Precision
  !  Rombint returns the integral from a to b of using Romberg integration.
  !  The method converges provided that f(x) is continuous in (a,b).
  !  f must be real(dl) and must be declared external in the calling
  !  routine.  tol indicates the desired relative accuracy in the integral.
  !
  implicit none
  integer,parameter :: dl = KIND(1.D0)
  integer,parameter :: MAXITER=20
  integer,parameter :: MAXJ=5
  dimension g(MAXJ+1)
  real(dl) :: rombint_dcom_integrand_ray,dcom_integrand_ray
  real(dl), intent(in) :: a,b,tol,Om0,Ode0
  integer :: nint,i,k,jmax,j
  real(dl) :: h, gmax,error,g,g0,g1,fourj

  h=0.5D0*(b-a)
  gmax=h*(dcom_integrand_ray(a,Om0,Ode0)+dcom_integrand_ray(b,Om0,Ode0))
  g(1)=gmax
  nint=1
  error=1.0d20
  i=0
10 i=i+1
  if (i.gt.MAXITER.or.(i.gt.5.and.abs(error).lt.tol)) &
       go to 40
  !  Calculate next trapezoidal rule approximation to integral.
  g0=0.0D0
  do 20 k=1,nint
     g0=g0+dcom_integrand_ray(a+(k+k-1)*h,Om0,Ode0)
20   continue
     g0=0.5D0*g(1)+h*g0
     h=0.5D0*h
     nint=nint+nint
     jmax=min(i,MAXJ)
     fourj=1.0D0
     do 30 j=1,jmax
        !  Use Richardson extrapolation.
        fourj=4.0D0*fourj
        g1=g0+(g0-g(j))/(fourj-1.0D0)
        g(j)=g0
        g0=g1
30      continue
        if (abs(g0).gt.tol) then
           error=1.0D0-gmax/g0
        else
           error=gmax
        end if
        gmax=g0
        g(jmax+1)=g0
        go to 10
40      rombint_dcom_integrand_ray=g0
        if (i.gt.MAXITER.and.abs(error).gt.tol)  then
           write(*,*) 'Warning: Rombint failed to converge; '
           write (*,*)'integral, error, tol:', rombint_dcom_integrand_ray,error, tol
        end if

end function rombint_dcom_integrand_ray

