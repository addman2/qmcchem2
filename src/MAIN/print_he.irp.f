program print_he

  PROVIDE ezfio_filename

  implicit none
  integer                    :: i, n_theta
  double precision           :: pi, d_theta, theta 

  print *,  'Number of determinants                   : ', det_num
  print *,  'Number of unique alpha/beta determinants : ', det_alpha_num, det_beta_num
  print *,  'Closed-shell MOs                         : ', mo_closed_num
  print *,  'Number of MOs in determinants            : ', num_present_mos

  print *, ' j2e_type = ', j2e_type
  print *, ' mu_erf = ', mu_erf

  print *, ' j1e_type = ', j1e_type
  print *, ' j1e_coef = ', j1e_coef
  print *, ' j1e_expo = ', j1e_expo

  print *, ' env_type', env_type
  print *, ' env_ceof = ', env_coef
  print *, ' env_expo = ', env_expo


  print *, ' jastrow value:'
  print *, jast_value, jast_value_inv
  print *, ' '

  ! ---

  pi      = 3.14d0
  n_theta = 250
  d_theta = 2.d0 * pi / dble(n_theta)

  do i = 1, n_theta

    theta = -pi + dble(i-1) * d_theta

    elec_coord(1,1) = 0.5d0
    elec_coord(1,2) = 0.0d0
    elec_coord(1,3) = 0.0d0
    elec_coord(2,1) = 0.5d0 * dcos(theta)
    elec_coord(2,2) = 0.5d0 * dsin(theta)
    elec_coord(2,3) = 0.0d0
    TOUCH elec_coord

    print *, theta, psidet_right_value, psi_value

  enddo


end

