
!################################################################
!################################################################
!################################################################
!################################################################
subroutine check_ray_grid_refinement(i_arg)
  use amr_commons
  use amr_parameters
  use ray_commons
  use ray_parameters
  use pm_commons ! for debuging
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer,intent(in) :: i_arg

  integer  :: old_ind,old_grid,old_cell,old_level
  integer  :: new_ind,new_grid,new_cell,new_level
  integer  ::iray
  integer  :: ix,iy,iz
  integer  :: new_cpu,old_cpu
  real(dp) :: dx,x_g,y_g,z_g,x_r,y_r,z_r
  logical  :: out_of_bound                                         ! for test use only

#ifndef WITHOUTMPI
  integer,dimension(ncpu)                 :: nrays_comm
#endif

  !============================================
  ! determine the size of communication commons
  !============================================

  nrays_comm(1:ncpu) = 0
  
  do iray=1,ray_nrays

     old_grid  =     ray_grid(iray,1)                       ! index of old grid (local)
     old_cpu   =     ray_grid(iray,2)/1000
     old_ind   =    (ray_grid(iray,2)-old_cpu*1000)/100     ! index of old cell relative to old grid
     old_level = mod(ray_grid(iray,2)-old_cpu*1000 ,100)    ! old level #
     old_cell  = ncoarse+(old_ind-1)*ngridmax+old_grid      ! index of old cell (local)

     if(old_grid.eq.0) then
        write(*,*) 'check_ray: old_grid index incorrect; please check.'
        stop
     end if

     ! Grid has been destroyed
     if(father(old_grid).eq.0) then

!       if(ray_step_test) write(*,*) "check_ray: grid has been destroyed"
        
        new_level = old_level-1                             ! coarse level #
        new_cell  = ex_father(old_grid)                     ! index of carse cell (local);
                                                            ! ex_father stores ex father cell index of the killed grid
        new_ind   = (new_cell-ncoarse-1)/ngridmax+1         ! index of coarse cell relative to coarse grid (local)
        new_grid  = new_cell-(ncoarse+(new_ind-1)*ngridmax) ! index of coarse grid (local)

        new_cpu   = cpu_map(new_cell)                       ! coarse grid CPU #

!       write(*,*) "check_ray1",ray_id(iray),myid,new_cpu,old_level,new_level,old_grid,new_grid
        
        if(new_grid.eq.0) then
           write(*,*) 'check_ray: new_grid index out of bound; please check (after destroying grid).'
           stop
        end if
        if(new_level.le.0 .or. new_cpu.le.0) then
           write(*,*) 'check_ray: new_level out of bound (1); please check.',new_level,new_cpu
           stop
        end if

     ! Grid has been further refined
     else if(son(old_cell).gt.0) then
        
!       if(ray_step_test) write(*,*) "check_ray: grid has been refined"
        
        new_grid  = son(old_cell)       ! index of find grid (local)
        new_level = old_level+1         ! fine level #
        new_cpu   = cpu_map(old_cell)   ! fine grid CPU #

        call get_ind_from_ray_position(iray,new_level,new_grid,new_ind)

        if(new_grid.eq.0) then
           write(*,*) 'check_ray: new_grid index out of bound (2); please check.'
           stop
        end if

     ! Grid remains unchanged
     else
        new_grid  = old_grid
        new_level = old_level
        new_cpu   = old_cpu
        if(old_cpu.ne.cpu_map(father(old_grid))) new_cpu = cpu_map(father(old_grid))
        new_ind   = old_ind

!       write(*,*) "check_ray3",ray_id(iray),myid,new_cpu,old_level,new_level,old_grid,new_grid
     end if

     if(new_level.le.0 .or. new_cpu.le.0) then
        write(*,*) 'check_ray: new_level out of bound(3); please check.',new_level,new_cpu
        stop
     end if

     ! Treat myid's grids and others' grids differently
     if(new_cpu.eq.myid) then
        ray_grid(iray,1) =  new_grid 
     else
        ray_grid(iray,1) = -new_grid                  ! mark that the new grid belongs to a different CPU
        nrays_comm(new_cpu) = nrays_comm(new_cpu)+1   ! number of rays to be sent to that CPU
     end if
     ray_grid(iray,2) = 1000*new_cpu+100*new_ind+new_level

     ! Test if new grid actually contains the ray
     if(ray_step_test .and. .false.) then
        x_r = ray_coord(iray,1)*dsin(ray_coord(iray,2))*dcos(ray_coord(iray,3))
        y_r = ray_coord(iray,1)*dsin(ray_coord(iray,2))*dsin(ray_coord(iray,3))
        z_r = ray_coord(iray,1)*dcos(ray_coord(iray,2))
        dx  = 0.5D0**new_level
        out_of_bound = .false.
        if(x_r.lt.xg(new_grid,1)-ray_x_obs-dx .or. x_r.gt.xg(new_grid,1)-ray_x_obs+dx) then
           write(*,*) xg(new_grid,1)-ray_x_obs-dx, x_r, xg(new_grid,1)-ray_x_obs+dx
           out_of_bound = .true.
        endif
        if(y_r.lt.xg(new_grid,2)-ray_y_obs-dx .or. y_r.gt.xg(new_grid,2)-ray_y_obs+dx) then
           write(*,*) xg(new_grid,2)-ray_y_obs-dx, y_r, xg(new_grid,2)-ray_y_obs+dx
           out_of_bound = .true.
        endif
        if(z_r.lt.xg(new_grid,3)-ray_z_obs-dx .or. z_r.gt.xg(new_grid,3)-ray_z_obs+dx) then
           write(*,*) xg(new_grid,3)-ray_z_obs-dx, z_r, xg(new_grid,3)-ray_z_obs+dx
           out_of_bound = .true.
        endif
        if(out_of_bound) then
           write(*,*) 'check_ray: ray outside boundary of new grid; please check.', new_level, new_grid, "(", xg(new_grid,1), xg(new_grid,2), xg(new_grid,3), ")  (", x_r, y_r, z_r, ") (", &
                xg(new_grid,1)-ray_x_obs-dx, xg(new_grid,2)-ray_y_obs-dx, xg(new_grid,3)-ray_z_obs-dx, ") (", &
                xg(new_grid,1)-ray_x_obs+dx, xg(new_grid,2)-ray_y_obs+dx, xg(new_grid,3)-ray_z_obs+dx, ")"
           write(*,*) ray_coord(iray,1), dsin(ray_coord(iray,2)), dcos(ray_coord(iray,3))
           !stop
        end if
     endif

  end do ! loop over rays

  !-----------------
  ! Communicate rays
  !-----------------
  call ray_communicate(i_arg,0,nrays_comm)
  

end subroutine check_ray_grid_refinement

!===========================================
! ray_communicate
!===========================================
subroutine ray_communicate(i_arg,j_arg,nrays_comm)

  use amr_commons
  use amr_parameters
  use ray_commons
  use ray_parameters
  use pm_commons ! for debuging
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  integer,intent(in) :: i_arg,j_arg
  integer,intent(in),dimension(ncpu)                 :: nrays_comm
  
  integer  :: new_grid,new_level
  integer  :: iray,jray,kray,nrays1,nrays2,nrays,nrays_tot
  integer  :: icpu,icpu2,ig
  integer  :: test_count

  integer::info,tag=101
  integer,dimension(ncpu) :: reqsend,reqrecv
#ifndef WITHOUTMPI
  integer,dimension(ncpu)                 :: ig_loc
  integer,dimension(ncpu)                 :: sendbuf,recvbuf
  integer,dimension(MPI_STATUS_SIZE,ncpu) :: statuses
  integer                                 :: countsend,countrecv,ncache
  integer, allocatable,dimension(:,:)     :: tmparr1
  real(dp),allocatable,dimension(:,:)     :: tmparr2
  integer                                 :: cur_size,new_size,size_kappa       ! BAOJIU-29-09
  integer                                 :: chunk_size                         ! BAOJIU-02-10
#endif
  
  !-----------------------------------------------------------------------------
  ! - This routine communicate rays among processors.
  ! - The rays to be communicated must be flagged with negative ray_grid(iray,1)
  ! - The number of rays to be communicated to each cpu should be given in
  !    nrays_comm(new_cpu)
  ! - This routine updates nrays_tot.
  !-----------------------------------------------------------------------------
  
#ifndef WITHOUTMPI

  !============================
  ! build communication commons
  !============================

  sendbuf = 0
  recvbuf = 0  
  ! Allocate emission commons to contain grid indices myid is going to send to other CPUs
  do icpu=1,ncpu
     ncache = nrays_comm(icpu)
     if(ncache>0) then
        allocate(ray_emission(icpu)%igrid (1:(nrays_comm(icpu))))
        allocate(ray_emission(icpu)%igrid2(1:(nrays_comm(icpu))))
        allocate(ray_emission(icpu)%ivar  (1:(nrays_comm(icpu))))
        allocate(ray_emission(icpu)%fvar  (1:(nrays_comm(icpu))))
        ray_emission(icpu)%ngrid = nrays_comm(icpu)                   ! number of rays to be sent to icpu
        sendbuf(icpu) = ray_emission(icpu)%ngrid
        ig_loc(icpu)  = 1                                            
     end if
  end do

  call MPI_ALLTOALL(sendbuf,1,MPI_INTEGER,recvbuf,1,MPI_INTEGER,MPI_COMM_WORLD,info)

  ! Allocate grid index
  do icpu=1,ncpu
     ray_reception(icpu)%ngrid=recvbuf(icpu)
     ncache = ray_reception(icpu)%ngrid
     if(ncache>0) then
        allocate(ray_reception(icpu)%igrid (1:ncache))
        allocate(ray_reception(icpu)%igrid2(1:ncache))
        allocate(ray_reception(icpu)%ivar  (1:ncache))
        allocate(ray_reception(icpu)%fvar  (1:ncache))
     end if
  end do

  ! Receive grid list    
  countrecv=0
  do icpu=1,ncpu
     ncache = ray_reception(icpu)%ngrid
     if(ncache>0) then
        countrecv = countrecv+1
        call MPI_IRECV(ray_reception(icpu)%igrid,ncache, &
                       & MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
     end if
  end do

  ! Gather emission arrays
  nrays2 = ray_nrays                                               ! temporarily store number of rays on myid
  do iray=1,ray_nrays
     ! Consider only rays no longer belonging to myid grids
     if(ray_grid(iray,1).gt.0) cycle

     ! Find the new CPU to which the ray belongs
     icpu = ray_grid(iray,2)/1000
     ! Find the new level to which the ray belongs
     new_level = mod(ray_grid(iray,2)-icpu*1000,100)

     if(reception(icpu,new_level)%ngrid.eq.0) cycle

     do ig=1,reception(icpu,new_level)%ngrid
        if(reception(icpu,new_level)%igrid(ig).eq.(abs(ray_grid(iray,1)))) then
           ray_emission(icpu)%igrid(ig_loc(icpu)) = ig
!          if(myid.eq.1) write(*,*) 'emission:',myid,icpu,ig_loc(icpu),new_level,ig,ray_grid(iray,2),ray_id(iray),ray_emission(2)%ngrid,ray_grid(iray,1)
           ig_loc(icpu) = ig_loc(icpu)+1
           exit
        end if
     end do

     nrays2 = nrays2-1
  end do

  ! Send grid list   
  countsend=0
  do icpu=1,ncpu
     ncache = ray_emission(icpu)%ngrid
     if(ncache>0) then
        countsend = countsend+1
        call MPI_ISEND(ray_emission(icpu)%igrid,ncache,  &
                       & MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,reqsend(countsend),info)
     end if
  end do

  ! Wait for full completion of sends
  call MPI_WAITALL(countsend,reqsend,statuses,info)

  ! Wait for full completion of receives
  call MPI_WAITALL(countrecv,reqrecv,statuses,info)

  ! Update ray number on myid. At the end of this: 
  ! ray_nrays remains as the old number of rays;
  ! nrays1 is the old number of rays plus rays to be received from other CPUs;
  ! nrays2 is the old number of rays minus rays to be sent to other CPUs.
  nrays1 = ray_nrays                                               ! temporarily store number of rays on myid
  do icpu=1,ncpu
     ncache = ray_reception(icpu)%ngrid                            ! number of rays received from icpu
     if(ncache.le.0) cycle

     if(icpu.ne.myid) then 
        nrays1 = nrays1+ncache                                     ! update number of rays on myid
     end if
  end do

  !==================
  ! do communications
  !==================

  cur_size = size(ray_id)                                          ! stores the current size of ray-related arrays
  if(ray_do_kappa) size_kappa = size(ray_kappa,2)                  ! BAOJIU-29-09

  kray = ray_nrays
  call make_virtual_ray_grid_int(ray_grid(:,2))

  chunk_size = 0                                                   ! BAOJIU-02-10
  do icpu=1,ncpu                                                   ! BAOJIU-02-10
     chunk_size = chunk_size+ray_reception(icpu)%ngrid             ! BAOJIU-02-10
  end do                                                           ! BAOJIU-02-10
  chunk_size = chunk_size-(cur_size-ray_nrays)+5                   ! BAOJIU-02-10
  if(chunk_size.lt.100) chunk_size = 100                           ! BAOJIU-02-10

  do icpu=1,ncpu

     ncache = ray_reception(icpu)%ngrid
     if(ncache.le.0)  cycle
     if(icpu.eq.myid) cycle
   
     do ig=1,ncache
        kray  = kray+1
        new_level = mod(ray_reception(icpu)%ivar(ig),100)
!       if(myid.eq.2) write(*,*) 'reception:',myid,icpu,ig,new_level,ray_reception(icpu)%igrid(ig),ray_reception(icpu)%ivar(ig),ray_reception(1)%ngrid
        new_grid  = emission(icpu,new_level)%igrid(ray_reception(icpu)%igrid(ig))
        if(new_grid.le.0) then
           write(*,*) 'ray_communicate: ray grid incorrect.',myid,icpu,new_level
           stop
        end if
        ! Extend array size if necessary
        if(kray.gt.cur_size) then
           new_size = cur_size+chunk_size

           allocate(tmparr1(1:cur_size,1:2))
           allocate(tmparr2(1:cur_size,1:3))

           ! Array for grid index info
           tmparr1(1:cur_size,1:2) = ray_grid(1:cur_size,1:2)
           deallocate(ray_grid)
           allocate(ray_grid(1:new_size,1:2))
           ray_grid = 0
           ray_grid(1:cur_size,1:2) = tmparr1(1:cur_size,1:2)

           ! Array for ray global id
           tmparr1(1:cur_size,1) = ray_id(1:cur_size)
           deallocate(ray_id)
           allocate(ray_id(1:new_size))
           ray_id = 0
           ray_id(1:cur_size) = tmparr1(1:cur_size,1)

           ! Array for ray WL convergence 
           if(ray_do_kappa) then                                                        ! BAOJIU-29-09
              tmparr2(1:cur_size,1:size_kappa) = ray_kappa(1:cur_size,1:size_kappa)     ! BAOJIU-29-09
              deallocate(ray_kappa)                                                     ! BAOJIU-29-09
              allocate(ray_kappa(1:new_size,1:size_kappa))                              ! BAOJIU-29-09
              ray_kappa = 0.0D0                                                         ! BAOJIU-29-09
              ray_kappa(1:cur_size,1:size_kappa) = tmparr2(1:cur_size,1:size_kappa)     ! BAOJIU-29-09
           end if                                                                       ! BAOJIU-29-09

           ! Array for ray WL shear
           if(ray_do_shear) then                                                        ! BAOJIU-29-09
              tmparr2(1:cur_size,1:2) = ray_shear(1:cur_size,1:2)                       ! BAOJIU-29-09
              deallocate(ray_shear)                                                     ! BAOJIU-29-09
              allocate(ray_shear(1:new_size,1:2))                                       ! BAOJIU-29-09
              ray_shear = 0.0D0                                                         ! BAOJIU-29-09
              ray_shear(1:cur_size,1:2) = tmparr2(1:cur_size,1:2)                       ! BAOJIU-29-09
           end if                                                                       ! BAOJIU-29-09

           ! +++ add new stuff here to calculate other observables +++

           ! Start: Sownak - ISW (15/01/2016)
           if(ray_do_isw) then
              tmparr2(1:cur_size,1) = ray_phidot(1:cur_size)
              deallocate(ray_phidot)
              allocate  (ray_phidot(1:new_size))
              ray_phidot = 0.0D0
              ray_phidot(1:cur_size) = tmparr2(1:cur_size,1)
           end if
           ! End: Sownak

           ! Array for ray coordinates
           tmparr2(1:cur_size,1:3) = ray_coord(1:cur_size,1:3)
           deallocate(ray_coord)
           allocate(ray_coord(1:new_size,1:3))
           ray_coord = 0.0D0
           ray_coord(1:cur_size,1:3) = tmparr2(1:cur_size,1:3)
 
           deallocate(tmparr1)
           deallocate(tmparr2)

           cur_size = new_size
        end if
        ray_grid(kray,1) = new_grid
        ray_grid(kray,2) = ray_reception(icpu)%ivar(ig)
     end do 
  end do


  kray = ray_nrays
  call make_virtual_ray_grid_int(ray_id(:))
  do icpu=1,ncpu
     ncache = ray_reception(icpu)%ngrid
     if(ncache.le.0)  cycle
     if(icpu.eq.myid) cycle
   
     do ig=1,ncache
        kray = kray+1
        ray_id(kray) = ray_reception(icpu)%ivar(ig)
     end do 
  end do

  kray = ray_nrays
  call make_virtual_ray_grid_dp(ray_coord(:,1))
  do icpu=1,ncpu
     ncache = ray_reception(icpu)%ngrid
     if(ncache.le.0)  cycle
     if(icpu.eq.myid) cycle
   
     do ig=1,ncache
        kray = kray+1
        ray_coord(kray,1) = ray_reception(icpu)%fvar(ig)
     end do 
  end do

  kray = ray_nrays
  call make_virtual_ray_grid_dp(ray_coord(:,2))
  do icpu=1,ncpu
     ncache = ray_reception(icpu)%ngrid
     if(ncache.le.0)  cycle
     if(icpu.eq.myid) cycle
   
     do ig=1,ncache
        kray = kray+1
        ray_coord(kray,2) = ray_reception(icpu)%fvar(ig)
     end do 
  end do

  kray = ray_nrays
  call make_virtual_ray_grid_dp(ray_coord(:,3))
  do icpu=1,ncpu
     ncache = ray_reception(icpu)%ngrid
     if(ncache.le.0)  cycle
     if(icpu.eq.myid) cycle
   
     do ig=1,ncache
        kray = kray+1
        ray_coord(kray,3) = ray_reception(icpu)%fvar(ig)
     end do 
  end do

  if(ray_do_kappa) then                                            ! BAOJIU-29-09
     kray = ray_nrays
     call make_virtual_ray_grid_dp(ray_kappa(:,1))                 ! BAOJIU-29-09 
     do icpu=1,ncpu                       
        ncache = ray_reception(icpu)%ngrid       
        if(ncache.le.0)  cycle                 
        if(icpu.eq.myid) cycle           
   
        do ig=1,ncache            
           kray = kray+1
           ray_kappa(kray,1) = ray_reception(icpu)%fvar(ig)        ! BAOJIU-29-09
        end do 
     end do
  end if                                                           ! BAOJIU-29-09

  if(ray_do_kappa.and.ray_kappa_method.eq.3) then                  ! BAOJIU-29-09
     kray = ray_nrays                                              ! BAOJIU-29-09
     call make_virtual_ray_grid_dp(ray_kappa(:,2))                 ! BAOJIU-29-09 
     do icpu=1,ncpu                                                ! BAOJIU-29-09
        ncache = ray_reception(icpu)%ngrid                         ! BAOJIU-29-09
        if(ncache.le.0)  cycle                                     ! BAOJIU-29-09
        if(icpu.eq.myid) cycle                                     ! BAOJIU-29-09
                                                                   ! BAOJIU-29-09
        do ig=1,ncache                                             ! BAOJIU-29-09
           kray = kray+1                                           ! BAOJIU-29-09
           ray_kappa(kray,2) = ray_reception(icpu)%fvar(ig)        ! BAOJIU-29-09
        end do                                                     ! BAOJIU-29-09
     end do                                                        ! BAOJIU-29-09
  end if                                                           ! BAOJIU-29-09

  if(ray_do_shear) then                                            ! BAOJIU-29-09
     kray = ray_nrays
     call make_virtual_ray_grid_dp(ray_shear(:,1))
     do icpu=1,ncpu
        ncache = ray_reception(icpu)%ngrid
        if(ncache.le.0)  cycle
        if(icpu.eq.myid) cycle
   
        do ig=1,ncache
           kray = kray+1
           ray_shear(kray,1) = ray_reception(icpu)%fvar(ig)
        end do 
     end do

     kray = ray_nrays
     call make_virtual_ray_grid_dp(ray_shear(:,2))
     do icpu=1,ncpu
        ncache = ray_reception(icpu)%ngrid
        if(ncache.le.0)  cycle
        if(icpu.eq.myid) cycle
   
        do ig=1,ncache
           kray = kray+1
           ray_shear(kray,2) = ray_reception(icpu)%fvar(ig)
        end do 
     end do
  end if                                                           ! BAOJIU-29-09

  ! +++ add new stuff here to calculate other observables +++

  ! Start: Sownak - ISW (15/01/2016)
  if(ray_do_isw) then
     kray = ray_nrays
     call make_virtual_ray_grid_dp(ray_phidot(:))
     do icpu=1,ncpu
        ncache = ray_reception(icpu)%ngrid
        if(ncache.le.0) cycle
        if(icpu.eq.myid) cycle
        
        do ig=1,ncache
           kray = kray+1
           ray_phidot(kray) = ray_reception(icpu)%fvar(ig)
        end do
     end do
  end if
  ! End: Sownak

  ! Defrag arrays so that old myid rays are listed first
  jray = 0
  do iray=1,nrays1
     if(ray_grid (iray,1).gt.0) then
        jray                = jray+1
        ray_grid (jray,1:2) = ray_grid (iray,1:2)
        ray_id   (jray    ) = ray_id   (iray    )
        ray_coord(jray,1:3) = ray_coord(iray,1:3)
        if(ray_do_kappa) ray_kappa(jray,1:size_kappa) = ray_kappa(iray,1:size_kappa)! BAOJIU-29-09
        if(ray_do_shear) ray_shear(jray,1:2) = ray_shear(iray,1:2) ! BAOJIU-29-09
        ! +++ add new stuff here to calculate other observables +++

        ! Start: Sownak - ISW (19/01/2016)
        if(ray_do_isw) ray_phidot(jray) = ray_phidot(iray)
        ! End: Sownak

     endif
  end do

  ray_nrays = jray

#endif

  nrays = ray_nrays

#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(nrays,nrays_tot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
  nrays_tot=nrays
#endif
  if(nrays_tot.ne.ray_Nx*ray_Ny) then
     write(*,*) 'ray_communicate: incorrect ray number.'
     write(*,*) 'after communication myid ',myid,' has ', nrays1,' rays.'
     if(myid.eq.1) write(*,*) 'there are ',nrays_tot,' rays in total, while it should be ',ray_Nx*ray_Ny,'!'
     stop
  end if
  if(j_arg.eq.0) then
!    write(*,888) i_arg,nrays_tot,myid,nrays
  else
!    write(*,999) nrays_tot,myid,nrays
  endif

  ! Deallocate the ray communication commons
  do icpu=1,ncpu
     if(ray_emission (icpu)%ngrid>0) then
        ray_emission (icpu)%ngrid = 0
        deallocate(ray_emission (icpu)%igrid )
        deallocate(ray_emission (icpu)%igrid2)
        deallocate(ray_emission (icpu)%ivar  )
        deallocate(ray_emission (icpu)%fvar  )
     end if
     if(ray_reception(icpu)%ngrid>0) then
        ray_reception(icpu)%ngrid = 0
        deallocate(ray_reception(icpu)%igrid )
        deallocate(ray_reception(icpu)%igrid2)
        deallocate(ray_reception(icpu)%ivar  )
        deallocate(ray_reception(icpu)%fvar  )
     end if
  end do

888 format('Level',I3,' refined, nrays_tot =',I8,', nrays (myid =',I4,') =',I6)
999 format('Rays have changed CPUs, nrays_tot =',I8,', nrays (myid =',I4,') =',I6)

end subroutine ray_communicate
