module ray_commons 
  use amr_commons
  use ray_parameters

  real(dp), allocatable, dimension(:,:) :: ray_coord               ! Ray_coordinates (chi,theta,phi)
  real(dp), allocatable, dimension(:,:) :: ray_kappa               ! WL convergence
  real(dp), allocatable, dimension(:,:) :: ray_shear               ! WL shear (gamma_1, gamma_2)
  ! +++ add new stuff here to calculate other observables +++

  ! Start: Sownak - ISW (26/10/2015)
  real(dp), allocatable, dimension(:)   :: ray_phidot              ! Phi_dot for ISW calculation
  ! End: Sownak

  integer                               :: ray_nrays               ! number of rays on my CPU
  integer,  allocatable, dimension(:)   :: ray_id                  ! ray ID
  integer,  allocatable, dimension(:,:) :: ray_grid                ! (grid ind, info of cell_ind,icpu & ilevel)
  integer,  allocatable, dimension(:)   :: ex_father               ! father cell befoire kill_grid is called

  ! ---------------------------- !
  ! ------- BAOJIU-04-10 ------- !
  ! ---------------------------- !
  integer                                 :: ray_ncells            ! number of cells containing rays
  integer,  allocatable, dimension(:)     :: ray_in_cell           ! number of rays in a given cell 
  real(dp), allocatable, dimension(:,:)   :: ray_stored_dens       ! store corner density values
  real(dp), allocatable, dimension(:,:)   :: ray_stored_pots       ! store corner potential values
  real(dp), allocatable, dimension(:,:,:) :: ray_stored_force      ! store corner force values
  real(dp), allocatable, dimension(:,:,:) :: ray_stored_tidal      ! store corner tidal tensors
  ! +++ add new stuff here to calculate other observables +++

  ! Start: Sownak - ISW (26/10/2015)
  real(dp), allocatable, dimension(:,:)   :: ray_stored_phidot     ! store corner phi_dot values
  ! End: Sownak

  ! ---------------------------- !
  ! ------- BAOJIU-04-10 ------- !
  ! ---------------------------- !

  real(dp) :: ray_chi_s                                            ! source comoving distance

  ! Communication structure for rays
  type ray_communicator
    integer                             :: ngrid                   ! number of rays to be communicated
    integer ,dimension(:),pointer       :: igrid                   ! local index of the grid
    integer ,dimension(:),pointer       :: igrid2      
    integer ,dimension(:),pointer       :: ivar                    ! container for integer  variables to be communicated
    real(dp),dimension(:),pointer       :: fvar                    ! container for floating variables to be communicated
  end type ray_communicator


  ! Ray emission and reception communicators
  type(ray_communicator),allocatable,dimension(:) :: ray_emission
  type(ray_communicator),allocatable,dimension(:) :: ray_reception

  ! To store the value of aexp at the beggining of a particle time step
  ! (used to determine the maximum distance a ray can travel in a time step)

  ! Start: Sownak - ISW (28/10/2015)
  real(dp),allocatable,dimension(:) :: aexp_old_ray,aexp_new_ray
  ! End: Sownak

  ! Arrays for identification of neighbouring nodes
  integer, dimension(1:27,    1:8) :: ray_kkk,ray_lll
  integer, dimension(1:8, 1:8,1:8) :: ray_kkk_mmm,ray_lll_mmm
  integer, dimension(1:8,     1:8) :: ray_ccc
  real(dp),dimension(1:8         ) :: ray_bbbb
  integer, dimension(1:3, 1:2,1:8) :: ray_iii,ray_jjj
  integer, dimension(1:8,     1:8) :: ray_ooo
  
  !---BAOJIU---!
  ! Ray integration status flags
  logical :: ray_all_started                                       ! flag all rays have started  integration
  logical :: ray_all_finished                                      ! flag all rays have finished integration
  logical :: ray_end_of_sim                                        ! flag that we are ready to end the simulation

  ! Box ID in the tiling
  logical :: ray_first_box

  ! Flag to be used to decide ray status
  logical :: ray_initialized = .false.

  ! Ray info output flags                                          ! BAOJIU-17-02-2016
  integer :: ray_iout                                              ! BAOJIU-17-02-2016               
  !---BAOJIU---!

end module ray_commons
