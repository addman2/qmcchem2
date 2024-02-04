! Mu x envelop + Nu (1 - envelop)
! -------------------------------

! ---

 BEGIN_PROVIDER [double precision, jast_Mu_Nu_value    ]
&BEGIN_PROVIDER [double precision, jast_Mu_Nu_value_inv]

  implicit none
  integer          :: i
  double precision :: argexpo

  argexpo = 0.d0
  !DIR$ VECTOR ALIGNED
  !DIR$ LOOP COUNT (200)
  do i = 1, elec_num
    argexpo += jast_elec_Mu_Nu_value(i)
  enddo

  jast_Mu_Nu_value     = dexp(argexpo)
  jast_Mu_Nu_value_inv = 1.d0 / jast_Mu_Nu_value

END_PROVIDER

! ---

BEGIN_PROVIDER [double precision, jast_elec_Mu_Nu_value, (elec_num_8)]

  BEGIN_DOC  
  !
  ! J(i) = 0.5 \sum_{j!=i} [u_mu(rij) x v(ri) x v(rj) + u_nu(rij) x (1 - v(ri) x v(rj))] + u_{1e}(ri)
  ! 
  !     u_mu(rij) = 0.5 [rij (1-erf(mu rij)) - exp(-(mu rij)**2) / (pi**0.5 mu)]
  !     u_nu(rij) = 0.5 [rij (1-erf(nu rij)) - exp(-(nu rij)**2) / (pi**0.5 nu)]
  !
  END_DOC

  include '../constants.F'

  implicit none
  integer          :: i, j
  double precision :: tmp_ij
  double precision :: jmu_val(elec_num)
  double precision :: jnu_val(elec_num)

  PROVIDE env_type
  PROVIDE j1e_type
  PROVIDE mu_erf nu_erf

  if(env_type .ne. "None") then

    do i = 1, elec_num

      call get_jmu_val_r1(i, dble(mu_erf), elec_num, jmu_val)
      call get_jmu_val_r1(i, dble(nu_erf), elec_num, jnu_val)

      tmp_ij = 0.d0
      do j = 1, elec_num
        if(j==i) cycle
        tmp_ij = tmp_ij + jmu_val(j) * vi_env(i) * vi_env(j) + jnu_val(j) * (1.d0 - vi_env(i) * vi_env(j))
      enddo

      jast_elec_Mu_Nu_value(i) = 0.5d0 * tmp_ij
    enddo

  else

    do i = 1, elec_num

      call get_jmu_val_r1(i, dble(mu_erf), elec_num, jmu_val)

      tmp_ij = 0.d0
      do j = 1, elec_num
        if(j==i) cycle
        tmp_ij = tmp_ij + jmu_val(j)
      enddo

      jast_elec_Mu_Nu_value(i) = 0.5d0 * tmp_ij
    enddo

  endif ! env_type

  if(j1e_type .ne. "None") then
    do i = 1, elec_num
      jast_elec_Mu_Nu_value(i) = jast_elec_Mu_Nu_value(i) + jast_elec_1e_value(i)
    enddo
  endif

END_PROVIDER

! ---

 BEGIN_PROVIDER [double precision, jast_elec_Mu_Nu_grad_x, (elec_num_8)]
&BEGIN_PROVIDER [double precision, jast_elec_Mu_Nu_grad_y, (elec_num_8)]
&BEGIN_PROVIDER [double precision, jast_elec_Mu_Nu_grad_z, (elec_num_8)]
&BEGIN_PROVIDER [double precision, jast_elec_Mu_Nu_lapl  , (elec_num_8)]

  BEGIN_DOC  
  !
  ! \grad_i   J = \sum_{j!=i} \grad_i   u_{sym}(ri, rj)
  ! \grad_i^2 J = \sum_{j!=i} \grad_i^2 u_{sym}(ri, rj)
  !
  END_DOC

  include '../constants.F'

  implicit none
  integer          :: i, j
  double precision :: rij
  double precision :: jmu_val(elec_num), jmu_derx(elec_num), jmu_dery(elec_num), jmu_derz(elec_num), jmu_lapl(elec_num)
  double precision :: jnu_val(elec_num), jnu_derx(elec_num), jnu_dery(elec_num), jnu_derz(elec_num), jnu_lapl(elec_num)

  PROVIDE j1e_type
  PROVIDE env_type
  PROVIDE mu_erf nu_erf

  if(env_type .eq. "None") then

    do i = 1, elec_num

      call get_jmu_der_lap_r1(i, dble(mu_erf), elec_num, jmu_derx, jmu_dery, jmu_derz, jmu_lapl)

      jast_elec_Mu_Nu_grad_x(i) = 0.d0
      jast_elec_Mu_Nu_grad_y(i) = 0.d0
      jast_elec_Mu_Nu_grad_z(i) = 0.d0
      jast_elec_Mu_Nu_lapl  (i) = 0.d0
      do j = 1, elec_num
        if(i==j) cycle

        jast_elec_Mu_Nu_grad_x(i) += jmu_derx(j)
        jast_elec_Mu_Nu_grad_y(i) += jmu_dery(j)
        jast_elec_Mu_Nu_grad_z(i) += jmu_derz(j)
        jast_elec_Mu_Nu_lapl  (i) += jmu_lapl(j)
      enddo
    enddo

  else

    do i = 1, elec_num

      call get_jmu_val_r1(i, dble(mu_erf), elec_num, jmu_val)
      call get_jmu_val_r1(i, dble(nu_erf), elec_num, jnu_val)

      call get_jmu_der_lap_r1(i, dble(mu_erf), elec_num, jmu_derx, jmu_dery, jmu_derz, jmu_lapl)
      call get_jmu_der_lap_r1(i, dble(nu_erf), elec_num, jnu_derx, jnu_dery, jnu_derz, jnu_lapl)

      jast_elec_Mu_Nu_grad_x(i) = 0.d0
      jast_elec_Mu_Nu_grad_y(i) = 0.d0
      jast_elec_Mu_Nu_grad_z(i) = 0.d0
      jast_elec_Mu_Nu_lapl  (i) = 0.d0
      do j = 1, elec_num
        if(i==j) cycle

        jast_elec_Mu_Nu_grad_x(i) += jnu_derx(j) + ((jmu_derx(j)-jnu_derx(j)) * vi_env(i) + (jmu_val(j)-jnu_val(j)) * deriv_env_x(i)) * vi_env(j)
        jast_elec_Mu_Nu_grad_y(i) += jnu_dery(j) + ((jmu_dery(j)-jnu_dery(j)) * vi_env(i) + (jmu_val(j)-jnu_val(j)) * deriv_env_y(i)) * vi_env(j)
        jast_elec_Mu_Nu_grad_z(i) += jnu_derz(j) + ((jmu_derz(j)-jnu_derz(j)) * vi_env(i) + (jmu_val(j)-jnu_val(j)) * deriv_env_z(i)) * vi_env(j)
        jast_elec_Mu_Nu_lapl  (i) += jnu_lapl(j) + ((jmu_lapl(j)-jnu_lapl(j)) * vi_env(i)                &          
                                                   + 2.d0 * ((jmu_derx(j)-jnu_derx(j)) * deriv_env_x(i)  &
                                                           + (jmu_dery(j)-jnu_dery(j)) * deriv_env_y(i)  &
                                                           + (jmu_derz(j)-jnu_derz(j)) * deriv_env_z(i)) &
                                                   + (jmu_val(j)-jnu_val(j)) * lapl_env(i)) * vi_env(j)
      enddo
    enddo

  endif ! env_type

  if(j1e_type .ne. "None") then
    do i = 1, elec_num
      jast_elec_Mu_Nu_grad_x(i) = jast_elec_Mu_Nu_grad_x(i) + jast_elec_1e_grad_x(i)
      jast_elec_Mu_Nu_grad_y(i) = jast_elec_Mu_Nu_grad_y(i) + jast_elec_1e_grad_y(i)
      jast_elec_Mu_Nu_grad_z(i) = jast_elec_Mu_Nu_grad_z(i) + jast_elec_1e_grad_z(i)
      jast_elec_Mu_Nu_lapl  (i) = jast_elec_Mu_Nu_lapl  (i) + jast_elec_1e_lapl  (i)
    enddo
  endif

END_PROVIDER

! ---

