module ray_parameters
  use amr_parameters
 
  ! Number for fuzzy comparison
  real(dp) :: ray_epsilon = 1.0d-13

  ! Number of rays to follow in each direction
  integer  :: ray_Nx = 256 
  integer  :: ray_Ny = 256

  ! Opening of the field of view in x & y directions, in degrees
  real(dp) :: ray_opening_x = 5.0D0 
  real(dp) :: ray_opening_y = 5.0D0

  ! Angles of the centre of field of view, in degress
  real(dp) :: ray_cof_x = 0.0D0
  real(dp) :: ray_cof_y = 0.0D0

  ! Coordinate of the observer (in cartesian coordinates in units of boz size)
  ! By defaul it is the center of the z = 0 face
  real(dp) :: ray_x_obs = 0.5D0
  real(dp) :: ray_y_obs = 0.5D0
  real(dp) :: ray_z_obs = 0.0D0 

  ! Box size in Mpc/h
  ! This is the same as boxlen_ini but with double precision 
  real(dp) :: ray_Lbox = 250.0D0

  ! Flag indicating whether bending of light rays is required
  logical  :: ray_no_bending = .true.

  ! Flag indicating whether want whole simulation to finish after all rays are finished 
  logical  :: ray_stop_after_raytracing = .true.

  ! source redshift
  real(dp) :: ray_z_s = 1.0D0

  ! Tests:
  logical  :: ray_step_test = .false.

  ! Number of fields
  integer, parameter  :: RAY_NFIELDS = 12
 
  ! Set which quantities to be calculated
  logical  :: ray_do_kappa = .true.
  logical  :: ray_do_shear = .true.

  ! Method to computer kappa
  ! 1: method A
  ! 2: method B
  ! 3: method A & B
  integer  :: ray_kappa_method = 3

  ! +++ add new stuff here to calculate other observables +++
  ! +++ don't forget to update amr/read_params.f90 and namelist/cosmo.nml if needed +++

  ! Start: Sownak - ISW (26/10/2015)
  logical  :: ray_do_isw = .true.
  ! End: Sownak

  ! Set if want to integrate fields at cell centers ALEX-21-02-2016
  logical  :: ray_do_ngp = .false.

  ! Ray output related stuff                                       ! BAOJIU-17-02-2016
  integer,parameter :: RAY_MAXOUT    = 1000                        ! BAOJIU-17-02-2016
  logical           :: ray_multi_out = .false.                     ! BAOJIU-17-02-2016
  integer           :: ray_nout      = 0                           ! BAOJIU-17-02-2016               
  real(dp),dimension(1:RAY_MAXOUT) :: ray_zout = -9.9D0            ! BAOJIU-17-02-2016
  logical           :: ray_afresh    = .true.                      ! BAOJIU-17-02-2016

end module ray_parameters
