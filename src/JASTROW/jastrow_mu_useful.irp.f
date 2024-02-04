! Mu Jastrow x envelop
! ---------------------

! ---

BEGIN_PROVIDER [integer, List_all_comb_b2_size]

  implicit none

  List_all_comb_b2_size = 2**nucl_num

END_PROVIDER

! ---

BEGIN_PROVIDER [integer, List_all_comb_b2, (nucl_num, List_all_comb_b2_size)]

  implicit none
  integer :: i, j

  if(nucl_num .gt. 32) then
    print *, ' nucl_num = ', nucl_num, '> 32'
    stop
  endif

  List_all_comb_b2 = 0

  do i = 0, List_all_comb_b2_size-1
    do j = 0, nucl_num-1
      if (btest(i,j)) then
        List_all_comb_b2(j+1,i+1) = 1
      endif
    enddo
  enddo

END_PROVIDER

! ---

subroutine j_elec_Mu(r1, r2, je)

  BEGIN_DOC  
  !
  ! J(i,j) = 0.5 [ rij (1-erf(mu rij)) - exp(-(mu rij)**2) / (pi**0.5 mu) ] 
  !        x v(riA) x v(rjA)
  !
  END_DOC

  include '../constants.F'

  implicit none
  double precision, intent(in)  :: r1(3), r2(3)
  double precision, intent(out) :: je
  double precision              :: mu_rij
  double precision              :: dx, dy, dz, rij, u_ij, vi, vj

  PROVIDE env_type

  dx     = r1(1) - r2(1)
  dy     = r1(2) - r2(2)
  dz     = r1(3) - r2(3)
  rij    = dsqrt(dx*dx + dy*dy + dz*dz)
  mu_rij = mu_erf * rij
  u_ij   = 0.5d0 * (rij * (1.d0 - derf(mu_rij)) - dexp(-mu_rij*mu_rij)/(dsqpi*mu_erf))

  if(env_type .eq. "None") then
    je = u_ij
  else
    call j_elec_env(r1, vi)
    call j_elec_env(r2, vj)
    je = u_ij * vi * vj
  endif

  return
end

! ---

subroutine j_elec_env(r, jenv)

  implicit none
  double precision, intent(in)  :: r(3)
  double precision, intent(out) :: jenv
  integer                       :: iA
  double precision              :: a, c, riA, r2
  double precision              :: dx, dy, dz

  PROVIDE env_type

  if(env_type .eq. "None") then

    jenv = 0.d0

  elseif(env_type .eq. "Sum_Slat") then

    jenv = 1.d0
    !DIR$ LOOP COUNT (100)
    do iA = 1, nucl_num
      a   = env_expo(iA)
      c   = env_coef(iA)
      dx  = r(1) - nucl_coord(iA,1)
      dy  = r(2) - nucl_coord(iA,2)
      dz  = r(3) - nucl_coord(iA,3)
      riA = dsqrt(dx*dx + dy*dy + dz*dz)
      jenv = jenv - c * dexp(-a*riA)
    enddo

  elseif(env_type .eq. "Prod_Gauss") then

    jenv = 1.d0
    !DIR$ LOOP COUNT (100)
    do iA = 1, nucl_num
      a   = env_expo(iA)
      dx  = r(1) - nucl_coord(iA,1)
      dy  = r(2) - nucl_coord(iA,2)
      dz  = r(3) - nucl_coord(iA,3)
      riA = dsqrt(dx*dx + dy*dy + dz*dz)
      jenv = jenv * (1.d0 - dexp(-a*riA*riA))
    enddo

  elseif(env_type .eq. "Sum_Gauss") then

    jenv = 1.d0
    !DIR$ LOOP COUNT (100)
    do iA = 1, nucl_num
      a   = env_expo(iA)
      c   = env_coef(iA)
      dx  = r(1) - nucl_coord(iA,1)
      dy  = r(2) - nucl_coord(iA,2)
      dz  = r(3) - nucl_coord(iA,3)
      r2  = dx*dx + dy*dy + dz*dz
      jenv = jenv - c * dexp(-a*r2)
    enddo

  elseif(env_type .eq. "Sum_Quartic") then

    jenv = 1.d0
    !DIR$ LOOP COUNT (100)
    do iA = 1, nucl_num
      a   = env_expo(iA)
      c   = env_coef(iA)
      dx  = r(1) - nucl_coord(iA,1)
      dy  = r(2) - nucl_coord(iA,2)
      dz  = r(3) - nucl_coord(iA,3)
      r2  = dx*dx + dy*dy + dz*dz
      jenv = jenv - c * dexp(-a*r2*r2)
    enddo

  else

    print *, ' Error in j_elec_env: Unknown env_type = ', env_type
    stop

  endif

  return
end

! ---

 BEGIN_PROVIDER [double precision, jast_elec_Mu_grad_x_num, (elec_num_8)]
&BEGIN_PROVIDER [double precision, jast_elec_Mu_grad_y_num, (elec_num_8)]
&BEGIN_PROVIDER [double precision, jast_elec_Mu_grad_z_num, (elec_num_8)]
&BEGIN_PROVIDER [double precision, jast_elec_Mu_lapl_num  , (elec_num_8)]

  implicit none
  integer          :: i, j
  double precision :: eps, tmp_der, tmp_lap, je_p, je_m, je_0, tmp
  double precision :: r1(3), r2(3)

  eps     = 1d-3
  tmp_der = 0.5d0 /  eps
  tmp_lap = 1.0d0 / (eps * eps)

  do i = 1, elec_num

    r1(1) = elec_coord(i,1)
    r1(2) = elec_coord(i,2)
    r1(3) = elec_coord(i,3)

    jast_elec_Mu_grad_x_num(i) = 0.d0
    jast_elec_Mu_grad_y_num(i) = 0.d0
    jast_elec_Mu_grad_z_num(i) = 0.d0
    jast_elec_Mu_lapl_num  (i) = 0.d0

    do j = 1, elec_num
      if(j==i) cycle

      r2(1) = elec_coord(j,1)
      r2(2) = elec_coord(j,2)
      r2(3) = elec_coord(j,3)
      call j_elec_Mu(r1, r2, je_0)

      r1(1) += eps
      call j_elec_Mu(r1, r2, je_p)
      r1(1) -= 2.d0 * eps
      call j_elec_Mu(r1, r2, je_m)
      r1(1) += eps
      jast_elec_Mu_grad_x_num(i) += tmp_der * (je_p - je_m)
      jast_elec_Mu_lapl_num  (i) += tmp_lap * (je_p - 2.d0 * je_0 + je_m)

      r1(2) += eps
      call j_elec_Mu(r1, r2, je_p)
      r1(2) -= 2.d0 * eps
      call j_elec_Mu(r1, r2, je_m)
      r1(2) += eps
      jast_elec_Mu_grad_y_num(i) += tmp_der * (je_p - je_m)
      jast_elec_Mu_lapl_num  (i) += tmp_lap * (je_p - 2.d0 * je_0 + je_m)

      r1(3) += eps
      call j_elec_Mu(r1, r2, je_p)
      r1(3) -= 2.d0 * eps
      call j_elec_Mu(r1, r2, je_m)
      r1(3) += eps
      jast_elec_Mu_grad_z_num(i) += tmp_der * (je_p - je_m)
      jast_elec_Mu_lapl_num  (i) += tmp_lap * (je_p - 2.d0 * je_0 + je_m)
    enddo
  enddo

END_PROVIDER

! ---

subroutine get_jmu_val_r1(i, mu, e_num, jval)

  BEGIN_DOC  
  !
  ! u(i,j) = 0.5 [rij (1 - erf(mu rij)) - exp(-(mu rij)**2)/(pi**0.5 mu)] 
  !
  END_DOC

  include '../constants.F'

  implicit none
  integer,          intent(in)  :: i, e_num
  double precision, intent(in)  :: mu
  double precision, intent(out) :: jval(e_num)

  integer                       :: j
  double precision              :: rij, mu_rij, mu_pi

  mu_pi = 1.d0 / (dsqpi * mu)

  do j = 1, e_num

    if(j == i) then
      jval(j) = 0.d0
    endif

    rij    = elec_dist(j,i)
    mu_rij = mu * rij

    jval(j) = 0.5d0 * (rij * (1.d0 - derf(mu_rij)) - mu_pi * dexp(-mu_rij*mu_rij))
  enddo

end

! ---

subroutine get_jmu_der_lap_r1(i, mu, e_num, jderx, jdery, jderz, jlapl)

  BEGIN_DOC  
  !
  ! d/dxi u(i,j) = 0.5 [(1 - erf(mu rij)) / rij] * (xi - xj)
  ! d/dyi u(i,j) = 0.5 [(1 - erf(mu rij)) / rij] * (yi - yj)
  ! d/dzi u(i,j) = 0.5 [(1 - erf(mu rij)) / rij] * (zi - zj)
  !
  ! \Delta_i u(i,j) = (1 - erf(mu rij)) / rij - mu exp[-(mu rij)^2] / sqrt(pi)
  !
  END_DOC

  include '../constants.F'

  implicit none
  integer,          intent(in)  :: i, e_num
  double precision, intent(in)  :: mu
  double precision, intent(out) :: jderx(e_num), jdery(e_num), jderz(e_num), jlapl(e_num)

  integer                       :: j
  double precision              :: rij, mu_rij, mu_div_sqrtpi
  double precision              :: tmp0_ij, tmp1_ij, tmp2_ij, tmp3_ij

  mu_div_sqrtpi = mu / dsqpi

  do j = 1, e_num

    if(j == i) then
      jderx(j) = 0.d0
      jdery(j) = 0.d0 
      jderz(j) = 0.d0 
      jlapl(j) = 0.d0
    endif

    rij    = elec_dist(j,i)
    mu_rij = mu * rij

    tmp0_ij = dexp(-mu_rij * mu_rij)

    tmp1_ij = 1.d0 - derf(mu_rij)
    tmp2_ij = tmp1_ij * elec_dist_inv(j,i)
    tmp3_ij = -0.5d0 * tmp2_ij

    jderx(j) = tmp3_ij * elec_dist_vec_x(j,i)
    jdery(j) = tmp3_ij * elec_dist_vec_y(j,i)
    jderz(j) = tmp3_ij * elec_dist_vec_z(j,i)
    jlapl(j) = tmp2_ij - mu_div_sqrtpi * tmp0_ij
  enddo

end

! ---

