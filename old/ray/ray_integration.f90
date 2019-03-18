!------------------------------------------
! ray_integrate
!------------------------------------------
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>< 
subroutine ray_integrate(iray, ilevel, ind, chi_A, ray_fields, &
                         which_ray_quantity)
  !------------------------------------------------------------------
  ! Integrates [g(chi, chi_S) * quantity] along the ray path in a cell 
  ! 
  ! g(chi, chi_S) = chi(chi_S-chi)/chi_S is the lensing kernel 
  ! quantity is any quantity we want to integrate (density, nabla^2Phi,
  !      nabla_inabla_i Phi, etc.
  ! 
  ! The notation used here tries to follow that of the paper
  ! Point A --> point in cell where the integration ends
  ! Point B --> point in cell where the integration starts
  ! Point 1 --> vertice 1 of cell, ie, that with lowest x,y,z coords
  ! 
  ! a_cell  --> diff. in x-dir between Point A and Point 1
  ! b_cell  --> diff. in y-dir between Point A and Point 1
  ! c_cell  --> diff. in z-dir between Point A and Point 1
  ! 
  ! alpha_* --> coefficients of the trilinear interpolation of int_qua 
  ! d_1234  --> d_N coefficients in the integral
  !
  ! equations are references to paper: https://arxiv.org/pdf/1601.02012.pdf
  !
  ! Parameters:
  ! -----------
  ! which_ray_quantity : int
  !     1 - \delta/a 
  !     2 - \nabla^1\nabla_1_\Phi + \nabla^2\nabla_2_\Phi = \nabla^2_\xi\Phi
  !     3 - \nabla^1\nabla_1_\Phi - \nabla^2\nabla_2_\Phi 
  !     4 - 2*\nabla^1\nabla_2_\Phi
  !     5 - phidot
  !     6 - \nabla_1_\Phi
  !     7 - \nabla_2_\Phi
  !
  !------------------------------------------------------------------
  use pm_commons
  use poisson_commons
  use amr_commons
  use amr_parameters
  use ray_commons
  use ray_parameters
  implicit none

  integer, intent(in )  :: iray,ilevel,ind
  integer, intent(in )  :: which_ray_quantity
  real(dp), intent(in ) :: chi_A
  !integer, intent(in ) :: icell
  real (dp), dimension(1:RAY_NFIELDS,1:8), intent(in) :: ray_fields
 
  real(dp) :: dx
  integer :: i
  real(dp) :: a_cell, b_cell, c_cell      ! Difference in x,y,z of point A and
                                          ! vertice 1 of cell 
  real(dp) :: R_cell                      ! chi_B - chi_A, ending point of the
                                          ! integral (recall paper notation)
  real(dp), dimension(1:3) :: xo          ! Cartesian coordinates of the observer
  real(dp), dimension(1:3) :: chiA_pos    ! Cartesian coordinates of point A
                                          ! wrt observer
  real(dp), dimension(1:3) :: xgo         ! Cartesian coordinates of grid centre
                                          ! wrt observer
  real(dp), dimension(1:3) :: xco         ! Cartesian coordinates of cell centre
                                          ! wrt observer
  real(dp), dimension(1:3) :: delta_xc    ! Cartesian coordinates of cell centre
                                          ! wrt grid centre
  real(dp), dimension(1:4) :: d_1234      ! Coefficients d_N (N=1..4) in the
                                          ! integral formula 
  real(dp), dimension(1:4) :: d_1234_corr ! Coefficients d_N (N=1..4) in the
                                          ! integral formula of the correction
  real(dp) :: alpha1, alpha2, alpha3, alpha4 ! trilinear interpolation coeff.
  real(dp) :: alpha5, alpha6, alpha7, alpha8 ! trilinear interpolation coeff.
  real(dp) :: alpha_nablachiphi1, alpha_nablachiphi2 ! trilinear interpol. coeff.
  real(dp) :: alpha_nablachiphi3, alpha_nablachiphi4 ! trilinear interpol. coeff.
  real(dp) :: alpha_nablachiphi5, alpha_nablachiphi6 ! trilinear interpol. coeff.
  real(dp) :: alpha_nablachiphi7, alpha_nablachiphi8 ! trilinear interpol. coeff.
  real(dp) :: sinthe, costhe, sinphi, cosphi, sin2phi, sin2the, cos2phi
  !real(dp) :: dint1, dint2, dint3, dint4
  real(dp) :: integral_cell, integral_cell_corr
  real(dp) :: Dx_A, Dy_A, Dz_A, Dx_B, Dy_B, Dz_B
  real(dp) :: nablachiphi_at_A, nablachiphi_at_B, g_at_A, g_at_B, correction_term
  real(dp) :: aexp_mean ! ALEX-21-02-2016

  ! Integrand (which will be a linear combination of ray_fields)
  real (dp), dimension(1:8) :: int_qua
  real (dp), dimension(1:8) :: nablachiphi !correction force term in Method A
                                           ! (which_ray_quantity = 1)
  
  ! Useful trigonometric quantities
  sinthe  = dsin(      ray_coord(iray, 2))
  costhe  = dcos(      ray_coord(iray, 2))
  sinphi  = dsin(      ray_coord(iray, 3))
  cosphi  = dcos(      ray_coord(iray, 3))
  sin2the = dsin(2.0D0*ray_coord(iray, 2))
  sin2phi = dsin(2.0D0*ray_coord(iray, 3))
  cos2phi = dcos(2.0D0*ray_coord(iray, 3))

  ! Define quantity to be integrated:
  if(which_ray_quantity==1) then 
     ! Loop over corners of this cell.
     do i=1,8
        ! Quantity is \delta/a eq. 43
        int_qua(i) = (ray_fields(1,i)-1.0D0)/aexp
        ! --------------------------------------------------------------------- !
        ! alexandre_block_alexandre_block_alexandre_block  ==
        ! This is to compute the corrections to the density integral for
        ! Method A , nablachiphi is \nabla_\chi\Phi; text above eq. 44
        nablachiphi(i) = -(sinthe*cosphi)*ray_fields( 8,i) - & !x term
                          (sinthe*sinphi)*ray_fields( 9,i) - & !y term
                          (costhe       )*ray_fields(10,i)     !z term
        ! --------------------------------------------------------------------- !
     enddo
  endif
  if(which_ray_quantity==2) then
     ! Loop over corners of this cell.
     do i=1,8
        ! Quantity is \nabla^2_\xi\Phi eq. 39 + eq. 40
        int_qua(i) = (sinphi**2 + costhe**2*cosphi**2)*ray_fields(2,i) + & !xx term
                     (cosphi**2 + costhe**2*sinphi**2)*ray_fields(3,i) + & !yy term
                     (sinthe**2                      )*ray_fields(4,i) + & !zz term
                     (-sinthe**2*sin2phi             )*ray_fields(5,i) + & !xy term
                     (-cosphi*sin2the                )*ray_fields(6,i) + & !xz term
                     (-sinphi*sin2the                )*ray_fields(7,i)     !yz term
        int_qua(i) = int_qua(i)
     enddo
  endif
  if(which_ray_quantity==3) then
     ! Loop over corners of this cell.
     do i=1,8
        ! Quantity is \nabla^1\nabla_1_\Phi - \nabla^2\nabla_2_\Phi
        ! eq. 39 - eq. 40
        int_qua(i) = ( sinphi**2 - costhe**2*cosphi**2)*ray_fields(2,i) + & !xx term
                     ( cosphi**2 - costhe**2*sinphi**2)*ray_fields(3,i) + & !yy term
                     (-sinthe**2                      )*ray_fields(4,i) + & !zz term
                     (-(costhe**2+1.0D0)*sin2phi      )*ray_fields(5,i) + & !xy term
                     ( cosphi*sin2the                 )*ray_fields(6,i) + & !xz term
                     ( sinphi*sin2the                 )*ray_fields(7,i)     !yz term
        int_qua(i) = int_qua(i)
     enddo
  endif
  if(which_ray_quantity==4) then
     ! Loop over corners of this cell.
     do i=1,8
        ! Quantity is 2*\nabla^1\nabla_2_\Phi eq. 41
        int_qua(i) = 2.0*( &
             (-costhe*cosphi*sinphi        )*ray_fields(2,i) + & !xx term
             ( costhe*cosphi*sinphi        )*ray_fields(3,i) + & !yy term
             ( 0.0D0                       )*ray_fields(4,i) + & !zz term
             ( costhe*(cosphi**2-sinphi**2))*ray_fields(5,i) + & !xy term
             ( sinphi*sinthe               )*ray_fields(6,i) + & !xz term
             (-cosphi*sinthe               )*ray_fields(7,i))    !yz term
        int_qua(i) = int_qua(i)
     enddo
  endif
  ! Start: Sownak - ISW (26/10/2015)
  if (which_ray_quantity==5) then
     ! Loop over corners of this cell.
     do i=1,8
        ! Quantity is phidot
        int_qua(i) = 2.0*ray_fields(12,i)
     end do
  end if
  ! End: Sownak
  ! Start: Christoph - deflection angle (30/01/19)
  if (which_ray_quantity==6) then
     ! Loop over corners of this cell.
     do i=1,8
        ! Quantity is \nable_1_\Phi
        int_qua(i) = (costhe*cosphi)*ray_fields(8,i) + &       !x term
                     (costhe*sinphi)*ray_fields(9,i) + &       !y term
                     (sinthe)*ray_fields(      10,i)           !z term
        int_qua(i) = int_qua(i)
     end do
  end if
  if (which_ray_quantity==7) then
     ! Loop over corners of this cell.
     do i=1,8
        ! Quantity is \nable_2_\Phi
        int_qua(i) = (sinthe*sinphi)*ray_fields(8,i)*(-1) + &  !x term
                     (sinthe*cosphi)*ray_fields(9,i)           !y term
        int_qua(i) = int_qua(i)
     end do
  end if
  ! End: Christoph
  
  ! Ending point of the integral: chi_B - chi_A (recall paper variable)
  R_cell = ray_coord(iray, 1) - chi_A

  ! cell size of ilevel
  dx = 0.5D0**ilevel

  ! Cartesian coordinates of Point A wrt the observer
  chiA_pos(1) = chi_A*dsin(ray_coord(iray,2))*dcos(ray_coord(iray,3))
  chiA_pos(2) = chi_A*dsin(ray_coord(iray,2))*dsin(ray_coord(iray,3))
  chiA_pos(3) = chi_A*dcos(ray_coord(iray,2))

  ! Cartesian coordinates of the observer
  xo(1) = ray_x_obs
  xo(2) = ray_y_obs
  xo(3) = ray_z_obs

  ! Cartesian coordinates of the grid centre wrt the observer
  do i=1,ndim
     xgo(i) = xg(ray_grid(iray,1),i)-xo(i)
  end do
  ! Cartesian coordinates of the cell centre wrt the grid centre
  call get_cell_center(ilevel,iray,ind,delta_xc)
  ! Cartesian coordinates of the cell centre wrt the observer   
  do i=1,ndim
     xco(i) = xgo(i)+delta_xc(i)
  end do

  ! Differences in x,y,z direction between Point A and Point 1  
  a_cell = chiA_pos(1) - (xco(1) - 0.5D0*dx) ! diff in x-dir
  b_cell = chiA_pos(2) - (xco(2) - 0.5D0*dx) ! diff in y-dir
  c_cell = chiA_pos(3) - (xco(3) - 0.5D0*dx) ! diff in z-dir

  ! alpha* coefficients of the trilinear interpolation; eq. 15
  alpha1 = int_qua(1)
  alpha2 = int_qua(2) - int_qua(1)
  alpha3 = int_qua(3) - int_qua(1)
  alpha4 = int_qua(5) - int_qua(1)
  alpha5 = int_qua(4) - int_qua(3) - int_qua(2) + int_qua(1) 
  alpha6 = int_qua(7) - int_qua(5) - int_qua(3) + int_qua(1) 
  alpha7 = int_qua(6) - int_qua(5) - int_qua(2) + int_qua(1) 
  alpha8 = int_qua(8) - int_qua(7) - int_qua(6) - int_qua(4) & 
         + int_qua(2) + int_qua(5) + int_qua(3) - int_qua(1)

  ! --------------------------------------------------------------------------- !
  ! alexandre_block_alexandre_block_alexandre_block  ==
  ! This is to compute the corrections to the density integral
  ! for Method A eq. ??
  alpha_nablachiphi1 = nablachiphi(1)
  alpha_nablachiphi2 = nablachiphi(2) - nablachiphi(1)
  alpha_nablachiphi3 = nablachiphi(3) - nablachiphi(1)
  alpha_nablachiphi4 = nablachiphi(5) - nablachiphi(1)
  alpha_nablachiphi5 = nablachiphi(4) - nablachiphi(3) - &
                       nablachiphi(2) + nablachiphi(1) 
  alpha_nablachiphi6 = nablachiphi(7) - nablachiphi(5) - &
                       nablachiphi(3) + nablachiphi(1) 
  alpha_nablachiphi7 = nablachiphi(6) - nablachiphi(5) - &
                       nablachiphi(2) + nablachiphi(1) 
  alpha_nablachiphi8 = nablachiphi(8) - nablachiphi(7) - &
                       nablachiphi(6) - nablachiphi(4) + & 
                       nablachiphi(2) + nablachiphi(5) + &
                       nablachiphi(3) - nablachiphi(1)
  ! alexandre_block_alexandre_block_alexandre_block
  ! --------------------------------------------------------------------------- !

  ! d_N coefficients in the integral formula; eq. B1 - B4
  d_1234(1) = alpha1 + (alpha2*a_cell + alpha3*b_cell + alpha4*c_cell)/dx &
              + (alpha5*a_cell*b_cell &
               + alpha6*b_cell*c_cell &
               + alpha7*a_cell*c_cell)/dx**2 &
              + alpha8*a_cell*b_cell*c_cell/dx**3
  d_1234(2) = (sinthe*cosphi*alpha2 + sinthe*sinphi*alpha3 + costhe*alpha4)/dx & 
              + (sinthe*cosphi*(alpha7*c_cell + alpha5*b_cell) &
               + sinthe*sinphi*(alpha6*c_cell + alpha5*a_cell) &
               + costhe*(alpha7*a_cell + alpha6*b_cell))/dx**2 &
              + (sinthe*cosphi*alpha8*b_cell*c_cell &
               + sinthe*sinphi*alpha8*a_cell*c_cell &
               + costhe* alpha8*a_cell*b_cell)/dx**3
  d_1234(3) = (sinthe**2*cosphi*sinphi*alpha5 +
               + sinthe*cosphi*costhe*alpha7 + &
               + sinthe*sinphi*costhe*alpha6)/dx**2 & 
              + (sinthe**2*cosphi*sinphi*alpha8*c_cell &
               + sinthe*cosphi*costhe*alpha8*b_cell &
               + sinthe*sinphi*costhe*alpha8*a_cell)/dx**3
  d_1234(4) = alpha8*(sinthe**2*costhe*cosphi*sinphi)/dx**3 
  
  ! Compute the integral in the cell; eq. 17
  integral_cell = 0.0D0
  do i = 1,4
      integral_cell = integral_cell + d_1234(i)*( &
                      R_cell**(i + 2)                          /(dble(i)+2.0D0) + & 
                      R_cell**(i + 1)*(2.0D0*chi_A - ray_chi_s)/(dble(i)+1.0D0) + &
                      R_cell**(i    )*(chi_A - ray_chi_s)*chi_A/(dble(i)))
  end do
  integral_cell = integral_cell/ray_chi_s

  ! Start: Sownak -ISW (29/10/2015)
  if(ray_do_isw.and.which_ray_quantity.eq.5) then
     integral_cell = 0.0D0
     do i = 1,4
        integral_cell = integral_cell + d_1234(i)*(-R_cell**(i) / (dble(i)))
     end do
  end if
  ! End: Sownak


  ! --------------------------------------------------------------------------- !
  ! --------------------------------------------------------------------------- !
  ! alexandre_block_alexandre_block_alexandre_block  ==  
  !     This is to compute the corrections to the density integral for Method A
  !     This modification block computes the term
  !     integral([g'\nabla_\chi\Phi]dchi)
  ! d_N coefficients in the integral formula
  d_1234_corr(1) = alpha_nablachiphi1 &
                   + (alpha_nablachiphi2*a_cell + &
                      alpha_nablachiphi3*b_cell + &
                      alpha_nablachiphi4*c_cell)/dx &
                   + (alpha_nablachiphi5*a_cell*b_cell + &
                      alpha_nablachiphi6*b_cell*c_cell + &
                      alpha_nablachiphi7*a_cell*c_cell )/dx**2 &
                   + alpha_nablachiphi8*a_cell*b_cell*c_cell/dx**3
  d_1234_corr(2) = (sinthe*cosphi*alpha_nablachiphi2 + &
                    sinthe*sinphi*alpha_nablachiphi3 + &
                    costhe*alpha_nablachiphi4)/dx & 
                   + (sinthe*cosphi*(alpha_nablachiphi7*c_cell + &
                                     alpha_nablachiphi5*b_cell) + &
                      sinthe*sinphi*(alpha_nablachiphi6*c_cell + &
                                     alpha_nablachiphi5*a_cell) + &
                      costhe*(alpha_nablachiphi7*a_cell + &
                      alpha_nablachiphi6*b_cell))/dx**2 &
                   + (sinthe*cosphi*alpha_nablachiphi8*b_cell*c_cell + &
                      sinthe*sinphi*alpha_nablachiphi8*a_cell*c_cell + &
                      costhe* alpha_nablachiphi8*a_cell*b_cell)/dx**3
  d_1234_corr(3) = (sinthe**2*cosphi*sinphi*alpha_nablachiphi5 + &
                    sinthe*cosphi*costhe*alpha_nablachiphi7 + &
                    sinthe*sinphi*costhe*alpha_nablachiphi6)/dx**2 & 
                   + (sinthe**2*cosphi*sinphi*alpha_nablachiphi8*c_cell + &
                      sinthe*cosphi*costhe*alpha_nablachiphi8*b_cell + &
                      sinthe*sinphi*costhe*alpha_nablachiphi8*a_cell)/dx**3
  d_1234_corr(4) = alpha_nablachiphi8*(sinthe**2*costhe*cosphi*sinphi)/dx**3 
  
  ! Compute the integral in the cell
  integral_cell_corr = 0.0D0
  do i = 1,4
  integral_cell_corr = integral_cell_corr + d_1234_corr(i) * &
                       (R_cell**(i+1)* 2.0D0                   /(dble(i) + 1.0D0) + &
                        R_cell**(i  )*(2.0D0*chi_A - ray_chi_s)/(dble(i)        ))
  end do
  integral_cell_corr = integral_cell_corr/ray_chi_s

  ! alexandre_block_alexandre_block_alexandre_block
  ! --------------------------------------------------------------------------- !
  
  ! --------------------------------------------------------------------------- !
  ! alexandre_block_alexandre_block_alexandre_block  ==
  !     This is to compute the corrections to the density integral for Method A
  !     This modification block computes the term
  !     integral([g\nabla_\chi\Phi]dchi)
  ! Compute the displacements Dx, Dy, Dz for point A
  Dx_A = a_cell/dx
  Dy_A = b_cell/dx
  Dz_A = c_cell/dx
  ! Compute the displacements Dx, Dy, Dz for point B
  Dx_B = (a_cell + R_cell*sinthe*cosphi)/dx
  Dy_B = (b_cell + R_cell*sinthe*sinphi)/dx
  Dz_B = (c_cell + R_cell*costhe       )/dx
  ! Compute \nabla_\chi\Phi at point A
  nablachiphi_at_A = alpha_nablachiphi1 &
                     + alpha_nablachiphi2*Dx_A &
                     + alpha_nablachiphi3*Dy_A &
                     + alpha_nablachiphi4*Dz_A &
                     + alpha_nablachiphi5*Dx_A*Dy_A &
                     + alpha_nablachiphi6*Dy_A*Dz_A &
                     + alpha_nablachiphi7*Dx_A*Dz_A &
                     + alpha_nablachiphi8*Dx_A*Dy_A*Dz_A
  ! Compute \nabla_\chi\Phi at point B
  nablachiphi_at_B = alpha_nablachiphi1 & 
                     + alpha_nablachiphi2*Dx_B &
                     + alpha_nablachiphi3*Dy_B &
                     + alpha_nablachiphi4*Dz_B &
                     + alpha_nablachiphi5*Dx_B*Dy_B &
                     + alpha_nablachiphi6*Dy_B*Dz_B &
                     + alpha_nablachiphi7*Dx_B*Dz_B &
                     & alpha_nablachiphi8*Dx_B*Dy_B*Dz_B
  ! Compute f(chi, chi_S) at point A and point B
  g_at_A = (ray_chi_s - chi_A             )*chi_A             /ray_chi_s
  g_at_B = (ray_chi_s - ray_coord(iray, 1))*ray_coord(iray, 1)/ray_chi_s

  correction_term = (g_at_B*nablachiphi_at_B - g_at_A*nablachiphi_at_A)

  ! alexandre_block_alexandre_block_alexandre_block
  ! --------------------------------------------------------------------------- !
  ! --------------------------------------------------------------------------- !

  if(which_ray_quantity.eq.1) then
     ! These ray_kappa are being computed in units of (km/s)^2 
     ! density integral + both corrections
     ray_kappa(iray,1) = ray_kappa(iray,1) &
                         + (1.5D0*omega_m*100.D0**2)*integral_cell*ray_Lbox**2 &
                         + (correction_term + integral_cell_corr) &
                         * ray_Lbox**2*100.0D0**2/aexp**2          ! BAOJIU-29-09
  end if
  if(which_ray_quantity.eq.2) then
     ! This should have units of (km/s)^2
     ray_kappa(iray,ray_kappa_method-1) = ray_kappa(iray,ray_kappa_method-1) &
                                          + integral_cell &
                                          * ray_Lbox**2*100.D0**2/aexp**2
  end if
  if(ray_do_shear.and.which_ray_quantity.eq.3) then
     !This should have units of (km/s)^2
     ray_shear(iray,1) = ray_shear(iray,1) &
                         + integral_cell &
                         * ray_Lbox**2*100.D0**2/aexp**2
  end if
  if(ray_do_shear.and.which_ray_quantity.eq.4) then
     !This should have units of (km/s)^2
     ray_shear(iray,2) = ray_shear(iray,2) &
                         + integral_cell &
                         * ray_Lbox**2*100.D0**2/aexp**2
  end if
  if(ray_do_deflection.and.which_ray_quantity.eq.6) then              ! Christoph
     !This should have units of (km/s)^2
     ray_deflection(iray,1) = ray_deflection(iray,1) &
                              + integral_cell &
                              * ray_Lbox**2*100.D0**2/aexp**2
  end if
  if(ray_do_deflection.and.which_ray_quantity.eq.7) then              ! Christoph
     !This should have units of (km/s)^2
     ray_deflection(iray,2) = ray_deflection(iray,2) &
                              + integral_cell &
                              * ray_Lbox**2*100.D0**2/aexp**2
  end if

  ! +++ add new stuff here to calculate other observables +++

  ! Start: Sownak - ISW (26/10/2015)
  if(ray_do_isw.and.which_ray_quantity.eq.5) then
     ! ALEX-21-02-2016 (factor of aexp was missing) 
     aexp_mean = (aexp_new_ray(ilevel)+aexp_old_ray(ilevel))/2.0D0
     ! ALEX-21-02-2016 (factor of aexp was missing)
     ray_phidot(iray) = ray_phidot(iray) &
                        + integral_cell*aexp_mean*ray_Lbox**3*100.0D0**3
  end if
  ! End: Sownak   

end subroutine ray_integrate
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>< 

! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>< 
! ALEX-21-02-2016 (whole subroutine ray_integrate_ngp below)
subroutine ray_integrate_ngp(iray, ilevel, ind, chi_A, ray_tidal_at_center, &
                             ray_phidot_at_center, which_ray_quantity)
  use pm_commons
  use poisson_commons
  use amr_commons
  use amr_parameters
  use ray_commons
  use ray_parameters
  implicit none

  integer, intent(in )  :: iray,ilevel,ind
  integer, intent(in )  :: which_ray_quantity
  real(dp), intent(in ) :: chi_A
  !integer, intent(in ) :: icell
  real (dp), dimension(1:6), intent(in) :: ray_tidal_at_center
  real (dp), intent(in) :: ray_phidot_at_center 
  !------------------------------------------------------------------
  ! Integrates [g(chi, chi_S) * quantity] along the ray path in a cell, but
  ! assuming quantity is constant inside the cell, so in practice in integrates:
  ! g(chi, chi_S) and multiplies by the quantity
  ! 
  ! g(chi, chi_S) = chi(chi_S-chi)/chi_S is the lensing kernel 
  ! quantity is any quantity we want (density, nabla^2Phi,
  !      nabla_inabla_i Phi, etc.
  ! 
  ! The notation used here tries to follow that of the paper
  ! Point A --> point in cell where the integration ends
  ! Point B --> point in cell where the integration starts
  ! 
  ! Note that for ISW, the kernel is a constant
  ! 
  !------------------------------------------------------------------
 
  real(dp) :: dx
  integer :: i
  real(dp) :: chi_B                       ! chi_B
  real(dp), dimension(1:3) :: xo          ! Cartesian coordinates of the observer
  real(dp), dimension(1:3) :: chiA_pos    ! Cartesian coordinates of point A     wrt observer
  real(dp), dimension(1:3) :: xgo         ! Cartesian coordinates of grid centre wrt observer
  real(dp), dimension(1:3) :: xco         ! Cartesian coordinates of cell centre wrt observer
  real(dp), dimension(1:3) :: delta_xc    ! Cartesian coordinates of cell centre wrt grid centre
  real(dp), dimension(1:4) :: d_1234      ! Coefficients d_N (N=1..4) in the integral formula 
  real(dp), dimension(1:4) :: d_1234_corr ! Coefficients d_N (N=1..4) in the integral formula of the correction
  real(dp) :: sinthe, costhe, sinphi, cosphi, sin2phi, sin2the, cos2phi
  real(dp) :: integral_cell, integral_cell_corr
  real(dp) :: int_qua, aexp_mean
  
  ! Useful trigonometric quantities
  sinthe  = dsin(      ray_coord(iray, 2))
  costhe  = dcos(      ray_coord(iray, 2))
  sinphi  = dsin(      ray_coord(iray, 3))
  cosphi  = dcos(      ray_coord(iray, 3))
  sin2the = dsin(2.0D0*ray_coord(iray, 2))
  sin2phi = dsin(2.0D0*ray_coord(iray, 3))
  cos2phi = dcos(2.0D0*ray_coord(iray, 3))

  ! Define quantity to be integrated:
  if(which_ray_quantity==1) then 
     write(*,*) 'In ray_integrate_constantincell : this routine only does Method B for now'
  endif
  if(which_ray_quantity==2) then
     ! Quantity is \nabla^2_\xi\Phi
     int_qua = (sinphi**2 + costhe**2*cosphi**2)*ray_tidal_at_center(1) + & !xx term
               (cosphi**2 + costhe**2*sinphi**2)*ray_tidal_at_center(2) + & !yy term
               (sinthe**2                      )*ray_tidal_at_center(3) + & !zz term
               (-sinthe**2*sin2phi             )*ray_tidal_at_center(4) + & !xy term
               (-cosphi*sin2the                )*ray_tidal_at_center(5) + & !xz term
               (-sinphi*sin2the                )*ray_tidal_at_center(6)     !yz term
  endif
  if(which_ray_quantity==3) then
     ! Quantity is \nabla^1\nabla_1_\Phi - \nabla^2\nabla_2_\Phi
     int_qua = ( sinphi**2 - costhe**2*cosphi**2)*ray_tidal_at_center(1) + & !xx term
               ( cosphi**2 - costhe**2*sinphi**2)*ray_tidal_at_center(2) + & !yy term
               (-sinthe**2                      )*ray_tidal_at_center(3) + & !zz term
               (-(costhe**2+1.0D0)*sin2phi      )*ray_tidal_at_center(4) + & !xy term
               ( cosphi*sin2the                 )*ray_tidal_at_center(5) + & !xz term
               ( sinphi*sin2the                 )*ray_tidal_at_center(6)     !yz term
  endif
  if(which_ray_quantity==4) then
     ! Quantity is \nabla^1\nabla_2_\Phi
     int_qua = (-costhe* cosphi*sinphi       /sinthe)*ray_tidal_at_center(1) + & !xx term
               ( costhe* cosphi*sinphi       /sinthe)*ray_tidal_at_center(2) + & !yy term
               ( 0.0D0                              )*ray_tidal_at_center(3) + & !zz term
               ( costhe*(cosphi**2-sinphi**2)/sinthe)*ray_tidal_at_center(4) + & !xy term
               ( sinphi                             )*ray_tidal_at_center(5) + & !xz term
               (-cosphi                             )*ray_tidal_at_center(6)     !yz term
  endif

  if(which_ray_quantity==5) then                                         ! Sownak
     int_qua = 2.0d0 * ray_phidot_at_center 
  endif
  
  ! Ending point of the integral: chi_B (recall paper variable)
  chi_B = ray_coord(iray, 1)
  
  ! Compute the integral in the cell
  integral_cell = 0.0d0
  if(which_ray_quantity .lt. 5) then
     integral_cell = ray_chi_s*(chi_B**2 - chi_A**2)/2.0d0 - (chi_B**3 - chi_A**3)/3.0d0
     integral_cell = -int_qua*integral_cell/ray_chi_s
  end if
  if(which_ray_quantity==5) then
     integral_cell = int_qua*(chi_A - chi_B)
  endif

  if(which_ray_quantity.eq.1) then
     write(*,*) 'In ray_integrate_constantincell : this routine only does Method B for now'
  end if
  if(which_ray_quantity.eq.2) then
     ! This should have units of (km/s)^2
     ray_kappa(iray,ray_kappa_method-1) = ray_kappa(iray, ray_kappa_method-1) + &
                                          integral_cell*ray_Lbox**2 * &
                                          100.D0**2/aexp**2        ! BAOJIU-29-09
  end if
  if(ray_do_shear.and.which_ray_quantity.eq.3) then
      !This should have units of (km/s)^2
     ray_shear(iray,1) = ray_shear(iray,1) + &
                         integral_cell*ray_Lbox**2*100.D0**2/aexp**2
  end if
  if(ray_do_shear.and.which_ray_quantity.eq.4) then
     !This should have units of (km/s)^2
     ray_shear(iray,2) = ray_shear(iray,2) + &
                         integral_cell*ray_Lbox**2*100.D0**2/aexp**2
  end if
  if(ray_do_isw.and.which_ray_quantity.eq.5) then                        ! Sownak
     aexp_mean        = (aexp_new_ray(ilevel)+aexp_old_ray(ilevel))/2.0D0
     ray_phidot(iray) = ray_phidot(iray) + &
                        integral_cell*aexp_mean*ray_Lbox**3*100.0D0**3
  endif

end subroutine ray_integrate_ngp
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>< 
