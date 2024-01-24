! Mu Jastrow x envelop
! ---------------------

! ---

 BEGIN_PROVIDER [double precision, jast_Mu_value    ]
&BEGIN_PROVIDER [double precision, jast_Mu_value_inv]

  implicit none
  integer          :: i
  double precision :: argexpo

  argexpo = 0.d0
  !DIR$ VECTOR ALIGNED
  !DIR$ LOOP COUNT (200)
  do i = 1, elec_num
    argexpo += jast_elec_Mu_value(i)
  enddo

  jast_Mu_value     = dexp(argexpo)
  jast_Mu_value_inv = 1.d0 / jast_Mu_value

END_PROVIDER

! ---

BEGIN_PROVIDER [double precision, jast_elec_Mu_value, (elec_num_8)]

  BEGIN_DOC  
  !
  ! J(i) = 0.5 [\sum_{j!=i} u_{2e}(rij) x v(ri) x v(rj)] + u_{1e}(ri)
  ! 
  ! for(j2e_type .eq. "Mu")
  !
  !     u_{2e}(rij) = 0.5 [rij (1-erf(mu rij)) - exp(-(mu rij)**2) / (pi**0.5 mu)]
  !
  END_DOC

  include '../constants.F'

  implicit none
  integer          :: i, j
  double precision :: mu_rij, mu_pi
  double precision :: rij, tmp_ij

  PROVIDE env_type
  PROVIDE j1e_type

  mu_pi = 1.d0 / (dsqpi * mu_erf)

  if(env_type .ne. "None") then

    do i = 1, elec_num

      tmp_ij = 0.d0
      !DIR$ LOOP COUNT(100)
      do j = 1, elec_num

        if(j==i) cycle

        rij    = elec_dist(j,i)
        mu_rij = mu_erf * rij

        tmp_ij += (rij * (1.d0 - derf(mu_rij)) - mu_pi * dexp(-mu_rij*mu_rij)) * vi_env(j)
      enddo

      jast_elec_Mu_value(i) = 0.25d0 * tmp_ij * vi_env(i)
    enddo

  else

    do i = 1, elec_num

      tmp_ij = 0.d0
      !DIR$ LOOP COUNT(100)
      do j = 1, elec_num

        if(j==i) cycle

        rij    = elec_dist(j,i)
        mu_rij = mu_erf * rij

        tmp_ij += (rij * (1.d0 - derf(mu_rij)) - mu_pi * dexp(-mu_rij*mu_rij))
      enddo

      jast_elec_Mu_value(i) = 0.25d0 * tmp_ij
    enddo

  endif ! env_type

  if(j1e_type .ne. "None") then
    do i = 1, elec_num
      jast_elec_Mu_value(i) = jast_elec_Mu_value(i) + jast_elec_1e_value(i)
    enddo
  endif

END_PROVIDER

! ---

 BEGIN_PROVIDER [double precision, jast_elec_Mu_grad_x, (elec_num_8)]
&BEGIN_PROVIDER [double precision, jast_elec_Mu_grad_y, (elec_num_8)]
&BEGIN_PROVIDER [double precision, jast_elec_Mu_grad_z, (elec_num_8)]
&BEGIN_PROVIDER [double precision, jast_elec_Mu_lapl  , (elec_num_8)]

  BEGIN_DOC  
  !
  ! \grad_i   J = \sum_{j!=i} \grad_i   u_{sym}(ri, rj)
  ! \grad_i^2 J = \sum_{j!=i} \grad_i^2 u_{sym}(ri, rj)
  !
  END_DOC

  include '../constants.F'

  implicit none
  integer          :: i, j
  double precision :: mu_div_sqrtpi, mu_sqrtpi_inv, rij, mu_rij
  double precision :: tmp0_ij, tmp1_ij, tmp2_ij, tmp3_ij
  double precision :: vj_lapl_uij, vj_derivx_uij, vj_derivy_uij, vj_derivz_uij, vj_uij

  PROVIDE j1e_type
  PROVIDE env_type

  mu_div_sqrtpi = mu_erf / dsqpi

  if(env_type .eq. "None") then

    do i = 1, elec_num

      jast_elec_Mu_grad_x(i) = 0.d0
      jast_elec_Mu_grad_y(i) = 0.d0
      jast_elec_Mu_grad_z(i) = 0.d0
      jast_elec_Mu_lapl  (i) = 0.d0

      !DIR$ LOOP COUNT (100)
      do j = 1, elec_num
        if(i==j) cycle

        rij     = elec_dist(j,i)
        mu_rij  = mu_erf * rij
        tmp0_ij = dexp(-mu_rij * mu_rij)
        tmp1_ij = 1.d0 - derf(mu_rij)
        tmp2_ij = tmp1_ij * elec_dist_inv(j,i)
        tmp3_ij = -0.5d0 * tmp2_ij

        jast_elec_Mu_grad_x(i) += tmp3_ij * elec_dist_vec_x(j,i)
        jast_elec_Mu_grad_y(i) += tmp3_ij * elec_dist_vec_y(j,i)
        jast_elec_Mu_grad_z(i) += tmp3_ij * elec_dist_vec_z(j,i)
        jast_elec_Mu_lapl  (i) += tmp2_ij - mu_div_sqrtpi * tmp0_ij
      enddo
    enddo

  else

    mu_sqrtpi_inv = 1.d0 / (dsqpi * mu_erf)

    do i = 1, elec_num

      vj_uij        = 0.d0
      vj_derivx_uij = 0.d0 
      vj_derivy_uij = 0.d0 
      vj_derivz_uij = 0.d0 
      vj_lapl_uij   = 0.d0

      !DIR$ LOOP COUNT (100)
      do j = 1, elec_num

        if(i==j) cycle

        rij     = elec_dist(j,i)
        mu_rij  = mu_erf * rij
        tmp0_ij = dexp(-mu_rij * mu_rij)
        tmp1_ij = 1.d0 - derf(mu_rij)
        tmp2_ij = tmp1_ij * elec_dist_inv(j,i)
        tmp3_ij = -0.5d0 * tmp2_ij * vi_env(j)

        vj_uij        += 0.5d0 * (rij * tmp1_ij - mu_sqrtpi_inv * tmp0_ij) * vi_env(j)
        vj_derivx_uij += tmp3_ij * elec_dist_vec_x(j,i)
        vj_derivy_uij += tmp3_ij * elec_dist_vec_y(j,i)
        vj_derivz_uij += tmp3_ij * elec_dist_vec_z(j,i)
        vj_lapl_uij   += (tmp2_ij - mu_div_sqrtpi * tmp0_ij) * vi_env(j)
      enddo

      jast_elec_Mu_grad_x(i) = vj_derivx_uij * vi_env(i) + vj_uij * deriv_env_x(i)
      jast_elec_Mu_grad_y(i) = vj_derivy_uij * vi_env(i) + vj_uij * deriv_env_y(i)
      jast_elec_Mu_grad_z(i) = vj_derivz_uij * vi_env(i) + vj_uij * deriv_env_z(i)
      jast_elec_Mu_lapl  (i) = vj_lapl_uij * vi_env(i)                      &
                                + 2.d0 * ( vj_derivx_uij * deriv_env_x(i)   &
                                         + vj_derivy_uij * deriv_env_y(i)   & 
                                         + vj_derivz_uij * deriv_env_z(i) ) &
                                + vj_uij * lapl_env(i)
    enddo

  endif ! env_type

  if(j1e_type .ne. "None") then
    do i = 1, elec_num
      jast_elec_Mu_grad_x(i) = jast_elec_Mu_grad_x(i) + jast_elec_1e_grad_x(i)
      jast_elec_Mu_grad_y(i) = jast_elec_Mu_grad_y(i) + jast_elec_1e_grad_y(i)
      jast_elec_Mu_grad_z(i) = jast_elec_Mu_grad_z(i) + jast_elec_1e_grad_z(i)
      jast_elec_Mu_lapl  (i) = jast_elec_Mu_lapl  (i) + jast_elec_1e_lapl  (i)
    enddo
  endif

END_PROVIDER

! ---

 BEGIN_PROVIDER [double precision, grad_j_mu_x, (elec_num, elec_num)]
&BEGIN_PROVIDER [double precision, grad_j_mu_y, (elec_num, elec_num)]
&BEGIN_PROVIDER [double precision, grad_j_mu_z, (elec_num, elec_num)]

  BEGIN_DOC  
  !
  ! useful for three_body_mu calculation 
  !
  END_DOC

  include '../constants.F'

  implicit none
  integer          :: i, j
  double precision :: mu_div_sqrtpi, mu_sqrtpi_inv, rij, mu_rij
  double precision :: tmp0_ij, tmp1_ij, tmp2_ij, tmp3_ij, tmp
  double precision :: vj_lapl_uij, vj_derivx_uij, vj_derivy_uij, vj_derivz_uij, vj_uij

  mu_div_sqrtpi = mu_erf / dsqpi

  grad_j_mu_x = 0.d0
  grad_j_mu_y = 0.d0
  grad_j_mu_z = 0.d0

  if(env_type .eq. "None") then

    do i = 1, elec_num
      do j = 1, elec_num
        if(i==j) cycle

        rij     = elec_dist(j,i)
        mu_rij  = mu_erf * rij
        tmp0_ij = dexp(-mu_rij * mu_rij)
        tmp1_ij = 1.d0 - derf(mu_rij)
        tmp2_ij = tmp1_ij * elec_dist_inv(j,i)
        tmp3_ij = -0.5d0 * tmp2_ij

        grad_j_mu_x(j,i) = tmp3_ij * elec_dist_vec_x(j,i)
        grad_j_mu_y(j,i) = tmp3_ij * elec_dist_vec_y(j,i)
        grad_j_mu_z(j,i) = tmp3_ij * elec_dist_vec_z(j,i)
      enddo
    enddo

  else

    mu_sqrtpi_inv = 1.d0 / (dsqpi * mu_erf)

    do i = 1, elec_num
      do j = 1, elec_num

        if(i==j) cycle

        rij     = elec_dist(j,i)
        mu_rij  = mu_erf * rij
        tmp0_ij = dexp(-mu_rij * mu_rij)
        tmp1_ij = 1.d0 - derf(mu_rij)
        tmp2_ij = tmp1_ij * elec_dist_inv(j,i)
        tmp3_ij = -0.5d0 * tmp2_ij * vi_env(j)

        vj_uij        = 0.5d0 * (rij * tmp1_ij - mu_sqrtpi_inv * tmp0_ij) * vi_env(j)
        vj_derivx_uij = tmp3_ij * elec_dist_vec_x(j,i)
        vj_derivy_uij = tmp3_ij * elec_dist_vec_y(j,i)
        vj_derivz_uij = tmp3_ij * elec_dist_vec_z(j,i)

        grad_j_mu_x(j,i) = vj_derivx_uij * vi_env(i) + vj_uij * deriv_env_x(i)
        grad_j_mu_y(j,i) = vj_derivy_uij * vi_env(i) + vj_uij * deriv_env_y(i)
        grad_j_mu_z(j,i) = vj_derivz_uij * vi_env(i) + vj_uij * deriv_env_z(i)
      enddo
    enddo

  endif ! env_type

END_PROVIDER 

! ---

BEGIN_PROVIDER [double precision, deltaE_Jmu_lapl]

  implicit none
  integer          :: i
  double precision :: tmp

  tmp = 0.d0
  do i = 1, elec_num
    tmp -= jast_elec_Mu_lapl(i)
  enddo
  deltaE_Jmu_lapl = 0.5d0 * tmp

END_PROVIDER

! ---

BEGIN_PROVIDER [double precision, deltaE_Jmu_nonh]

  implicit none
  integer          :: i
  double precision :: tmp

  tmp = 0.d0
  do i = 1, elec_num
    tmp -= ( psidet_right_grad_lapl(1,i) * jast_elec_Mu_grad_x(i) &
           + psidet_right_grad_lapl(2,i) * jast_elec_Mu_grad_y(i) &
           + psidet_right_grad_lapl(3,i) * jast_elec_Mu_grad_z(i) ) * psidet_right_inv
  enddo
  deltaE_Jmu_nonh = tmp 

END_PROVIDER

! ---

BEGIN_PROVIDER [double precision, deltaE_Jmu_grad]

  implicit none
  integer          :: i
  double precision :: tmp

  tmp = 0.d0
  do i = 1, elec_num
    tmp -= jast_elec_Mu_grad_x(i) * jast_elec_Mu_grad_x(i) &
         + jast_elec_Mu_grad_y(i) * jast_elec_Mu_grad_y(i) &
         + jast_elec_Mu_grad_z(i) * jast_elec_Mu_grad_z(i)
  enddo
  deltaE_Jmu_grad = 0.5d0 * tmp

END_PROVIDER

! ---

BEGIN_PROVIDER [double precision, deltaE_Jmu_grad_2b]

  implicit none
  integer          :: i, j
  double precision :: mu, tmp, rij, a

  mu  = mu_erf
  tmp = 0.d0
  do i = 1, elec_num
    !DIR$ LOOP COUNT(100)
    do j = 1, elec_num
      if(i==j) cycle
      rij = mu * elec_dist(j,i)
      a   = 1.d0 - derf(rij)
      tmp = tmp - a * a
    enddo
  enddo

  ! x 0.50 for double couting
  ! x 0.25 formula
  deltaE_Jmu_grad_2b = 0.125d0 * tmp

END_PROVIDER

! ---

BEGIN_PROVIDER [double precision, deltaE_Jmu_grad_3b]

  implicit none
  integer :: i, j, k

  deltaE_Jmu_grad_3b = 0.d0
  do i = 1, elec_num
    do j = i+1, elec_num
      do k = j+1, elec_num
        deltaE_Jmu_grad_3b -= grad_j_mu_x(i,j) * grad_j_mu_x(i,k)
        deltaE_Jmu_grad_3b -= grad_j_mu_y(i,j) * grad_j_mu_y(i,k)
        deltaE_Jmu_grad_3b -= grad_j_mu_z(i,j) * grad_j_mu_z(i,k)

        deltaE_Jmu_grad_3b -= grad_j_mu_x(j,i) * grad_j_mu_x(j,k)
        deltaE_Jmu_grad_3b -= grad_j_mu_y(j,i) * grad_j_mu_y(j,k)
        deltaE_Jmu_grad_3b -= grad_j_mu_z(j,i) * grad_j_mu_z(j,k)

        deltaE_Jmu_grad_3b -= grad_j_mu_x(k,i) * grad_j_mu_x(k,j)
        deltaE_Jmu_grad_3b -= grad_j_mu_y(k,i) * grad_j_mu_y(k,j)
        deltaE_Jmu_grad_3b -= grad_j_mu_z(k,i) * grad_j_mu_z(k,j)
      enddo
    enddo
  enddo

END_PROVIDER

! ---


