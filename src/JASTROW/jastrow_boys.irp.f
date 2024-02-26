! Boys's Jastrow 
! ---------------

! ---

 BEGIN_PROVIDER [double precision, jast_Boys_value    ]
&BEGIN_PROVIDER [double precision, jast_Boys_value_inv]

  implicit none
  integer          :: i
  double precision :: argexpo

  argexpo = 0.d0
  !DIR$ VECTOR ALIGNED
  !DIR$ LOOP COUNT (200)
  do i = 1, elec_num
    argexpo += jast_elec_Boys_value(i)
  enddo

  jast_Boys_value     = dexp(argexpo)
  jast_Boys_value_inv = 1.d0 / jast_Boys_value

END_PROVIDER

! ---

BEGIN_PROVIDER [double precision, jast_elec_Boys_value, (elec_num_8)]

  BEGIN_DOC  
  !
  ! J(i) = 0.5 [\sum_{j!=i} u_{2e}(rij) x v(ri) x v(rj)] + u_{1e}(ri)
  ! 
  ! for(j2e_type .eq. "Boys")
  !
  !     u_{2e}(rij) = 0.5 rij / (1 + a x rij)
  !
  END_DOC

  include '../constants.F'

  implicit none
  integer          :: i, j
  double precision :: rij, tmp_ij

  PROVIDE env_type
  PROVIDE j1e_type

  if(env_type .ne. "None") then

    do i = 1, elec_num

      tmp_ij = 0.d0
      !DIR$ LOOP COUNT(100)
      do j = 1, elec_num

        if(j==i) cycle

        rij = elec_dist(j,i)

        tmp_ij += rij / (1.d0 + a_Boys*rij)
      enddo

      jast_elec_Boys_value(i) = 0.25d0 * tmp_ij * vi_env(i)
    enddo

  else

    do i = 1, elec_num

      tmp_ij = 0.d0
      !DIR$ LOOP COUNT(100)
      do j = 1, elec_num

        if(j==i) cycle

        rij = elec_dist(j,i)

        tmp_ij += rij / (1.d0 + a_Boys*rij)
      enddo

      jast_elec_Boys_value(i) = 0.25d0 * tmp_ij
    enddo

  endif ! env_type

  if(j1e_type .ne. "None") then
    do i = 1, elec_num
      jast_elec_Boys_value(i) = jast_elec_Boys_value(i) + jast_elec_1e_value(i)
    enddo
  endif

END_PROVIDER

! ---

 BEGIN_PROVIDER [double precision, jast_elec_Boys_grad_x, (elec_num_8)]
&BEGIN_PROVIDER [double precision, jast_elec_Boys_grad_y, (elec_num_8)]
&BEGIN_PROVIDER [double precision, jast_elec_Boys_grad_z, (elec_num_8)]
&BEGIN_PROVIDER [double precision, jast_elec_Boys_lapl  , (elec_num_8)]

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
  double precision :: tmp0, tmp1, tmp2, tmp3
  double precision :: vj_lapl_uij, vj_derivx_uij, vj_derivy_uij, vj_derivz_uij, vj_uij

  PROVIDE j1e_type
  PROVIDE env_type

  if(env_type .eq. "None") then

    do i = 1, elec_num

      jast_elec_Boys_grad_x(i) = 0.d0
      jast_elec_Boys_grad_y(i) = 0.d0
      jast_elec_Boys_grad_z(i) = 0.d0
      jast_elec_Boys_lapl  (i) = 0.d0

      !DIR$ LOOP COUNT (100)
      do j = 1, elec_num
        if(i==j) cycle

        rij  = elec_dist(j,i)
        tmp0 = 1.d0 / (1.d0 + a_Boys * rij)
        tmp1 = tmp0 * tmp0
        tmp2 = elec_dist_inv(j,i) * tmp1
        tmp3 = -0.5d0 * tmp2

        jast_elec_Boys_grad_x(i) += tmp3 * elec_dist_vec_x(j,i)
        jast_elec_Boys_grad_y(i) += tmp3 * elec_dist_vec_y(j,i)
        jast_elec_Boys_grad_z(i) += tmp3 * elec_dist_vec_z(j,i)
        jast_elec_Boys_lapl  (i) += tmp2 * tmp0
      enddo
    enddo

  else

    do i = 1, elec_num

      vj_uij        = 0.d0
      vj_derivx_uij = 0.d0 
      vj_derivy_uij = 0.d0 
      vj_derivz_uij = 0.d0 
      vj_lapl_uij   = 0.d0

      !DIR$ LOOP COUNT (100)
      do j = 1, elec_num

        if(i==j) cycle

        rij  = elec_dist(j,i)
        tmp0 = 1.d0 / (1.d0 + a_Boys * rij)
        tmp1 = tmp0 * tmp0
        tmp2 = elec_dist_inv(j,i) * tmp1
        tmp3 = -0.5d0 * tmp2 * vi_env(j)

        vj_uij        += 0.5d0 * rij * tmp0 * vi_env(j)
        vj_derivx_uij += tmp3 * elec_dist_vec_x(j,i)
        vj_derivy_uij += tmp3 * elec_dist_vec_y(j,i)
        vj_derivz_uij += tmp3 * elec_dist_vec_z(j,i)
        vj_lapl_uij   += tmp2 * tmp0 * vi_env(j)
      enddo

      jast_elec_Boys_grad_x(i) = vj_derivx_uij * vi_env(i) + vj_uij * deriv_env_x(i)
      jast_elec_Boys_grad_y(i) = vj_derivy_uij * vi_env(i) + vj_uij * deriv_env_y(i)
      jast_elec_Boys_grad_z(i) = vj_derivz_uij * vi_env(i) + vj_uij * deriv_env_z(i)
      jast_elec_Boys_lapl  (i) = vj_lapl_uij * vi_env(i)                   &
                               + 2.d0 * ( vj_derivx_uij * deriv_env_x(i)   &
                                        + vj_derivy_uij * deriv_env_y(i)   &
                                        + vj_derivz_uij * deriv_env_z(i) ) &
                               + vj_uij * lapl_env(i)
    enddo

  endif ! env_type

  if(j1e_type .ne. "None") then
    do i = 1, elec_num
      jast_elec_Boys_grad_x(i) = jast_elec_Boys_grad_x(i) + jast_elec_1e_grad_x(i)
      jast_elec_Boys_grad_y(i) = jast_elec_Boys_grad_y(i) + jast_elec_1e_grad_y(i)
      jast_elec_Boys_grad_z(i) = jast_elec_Boys_grad_z(i) + jast_elec_1e_grad_z(i)
      jast_elec_Boys_lapl  (i) = jast_elec_Boys_lapl  (i) + jast_elec_1e_lapl  (i)
    enddo
  endif

END_PROVIDER

! ---

