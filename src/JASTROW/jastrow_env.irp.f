
! ---

BEGIN_PROVIDER [double precision, vi_env, (elec_num_8)]

  implicit none
  integer          :: i, iA
  integer          :: ii, phase, b
  double precision :: a, riA, tmp
  double precision :: expo, c

  if(env_type .eq. "None") then

    vi_env = 1.d0

  elseif(env_type .eq. "Sum_Slat") then

    do i = 1, elec_num
      tmp = 1.d0
      !DIR$ LOOP COUNT (100)
      do iA = 1, nucl_num
        a   = env_expo(iA)
        c   = env_coef(iA)
        riA = nucl_elec_dist(iA,i)
        tmp = tmp - c * dexp(-a*riA)
      enddo
      vi_env(i) = tmp
    enddo

  elseif(env_type .eq. "Prod_Gauss") then

    do i = 1, elec_num
      tmp = 1.d0
      !DIR$ LOOP COUNT (100)
      do iA = 1, nucl_num
        a   = env_expo(iA)
        riA = nucl_elec_dist(iA,i)
        tmp = tmp * (1.d0 - dexp(-a*riA*riA))
      enddo
      vi_env(i) = tmp
    enddo

  elseif(env_type .eq. "Sum_Gauss") then

    do i = 1, elec_num
      tmp = 1.d0
      !DIR$ LOOP COUNT (100)
      do iA = 1, nucl_num
        a   = env_expo(iA)
        c   = env_coef(iA)
        riA = nucl_elec_dist(iA,i)
        tmp = tmp - c * dexp(-a*riA*riA)
      enddo
      vi_env(i) = tmp
    enddo

  elseif(env_type .eq. "Sum_Quartic") then

    do i = 1, elec_num
      tmp = 1.d0
      !DIR$ LOOP COUNT (100)
      do iA = 1, nucl_num
        a   = env_expo(iA)
        c   = env_coef(iA)
        riA = nucl_elec_dist(iA,i)
        tmp = tmp - c * dexp(-a*riA*riA*riA*riA)
      enddo
      vi_env(i) = tmp
    enddo

  else

    print *, ' Error in vi_env: Unknown env_type = ', env_type
    stop

  endif

END_PROVIDER

! ---

 BEGIN_PROVIDER [double precision, deriv_env_x, (elec_num_8)]
&BEGIN_PROVIDER [double precision, deriv_env_y, (elec_num_8)]
&BEGIN_PROVIDER [double precision, deriv_env_z, (elec_num_8)]
&BEGIN_PROVIDER [double precision,    lapl_env, (elec_num_8)]

  implicit none
  integer          :: i, ii, iA, phase, b
  double precision :: a, riA, dx, dy, dz, r2, r4
  double precision :: tmp, tmpx, tmpy, tmpz, tmpl
  double precision :: expo, coef, coef_x, coef_y, coef_z, c
  double precision :: arg

  if(env_type .eq. "None") then

    deriv_env_x = 0.d0
    deriv_env_y = 0.d0
    deriv_env_z = 0.d0
    lapl_env    = 0.d0

  elseif(env_type .eq. "Sum_Slat") then

    do i = 1, elec_num
      tmpx = 0.d0
      tmpy = 0.d0
      tmpz = 0.d0
      tmpl = 0.d0
      !DIR$ LOOP COUNT (100)
      do iA = 1, nucl_num
        a = env_expo(iA)
        c = env_coef(iA)
        ! xi - xA = nucl_elec_dist_vec(1,iA,i)
        dx  = nucl_elec_dist_vec(1,iA,i)
        dy  = nucl_elec_dist_vec(2,iA,i)
        dz  = nucl_elec_dist_vec(3,iA,i)
        riA = nucl_elec_dist(iA,i)
        tmp = a * c * dexp(-a*riA) / riA
        tmpx = tmpx + tmp * dx
        tmpy = tmpy + tmp * dy
        tmpz = tmpz + tmp * dz
        tmpl = tmpl + tmp * (2.d0 - a*riA)
      enddo
      deriv_env_x(i) = tmpx
      deriv_env_y(i) = tmpy
      deriv_env_z(i) = tmpz
      lapl_env   (i) = tmpl
    enddo

  elseif(env_type .eq. "Prod_Gauss") then

    do i = 1, elec_num
      deriv_env_x(i) = 0.d0
      deriv_env_y(i) = 0.d0
      deriv_env_z(i) = 0.d0
      lapl_env   (i) = 0.d0
      do ii = 1, List_all_comb_b2_size
        phase  = 0
        expo   = 0.d0
        coef   = 0.d0
        coef_x = 0.d0
        coef_y = 0.d0
        coef_z = 0.d0
        !DIR$ LOOP COUNT (100)
        do iA = 1, nucl_num
          a   = env_expo(iA)
          b   = List_all_comb_b2(iA,ii)
          c   = dble(b) * a
          riA = nucl_elec_dist(iA,i)
          phase  += b
          coef   += c
          expo   += c * riA * riA
          ! xi - xA = nucl_elec_dist_vec(1,iA,i)
          coef_x += c * nucl_elec_dist_vec(1,iA,i)
          coef_y += c * nucl_elec_dist_vec(2,iA,i)
          coef_z += c * nucl_elec_dist_vec(3,iA,i)
        enddo
        tmp = -2.d0 * (-1.d0)**dble(phase) * dexp(-expo)
        deriv_env_x(i) += tmp * coef_x
        deriv_env_y(i) += tmp * coef_y
        deriv_env_z(i) += tmp * coef_z
        lapl_env   (i) += tmp * (3.d0 * coef - 2.d0 * (coef_x*coef_x + coef_y*coef_y + coef_z*coef_z))
      enddo
    enddo

  elseif(env_type .eq. "Sum_Gauss") then

    do i = 1, elec_num
      tmpx = 0.d0
      tmpy = 0.d0
      tmpz = 0.d0
      tmpl = 0.d0
      !DIR$ LOOP COUNT (100)
      do iA = 1, nucl_num
        a = env_expo(iA)
        c = env_coef(iA)
        ! xi - xA = nucl_elec_dist_vec(1,iA,i)
        dx  = nucl_elec_dist_vec(1,iA,i)
        dy  = nucl_elec_dist_vec(2,iA,i)
        dz  = nucl_elec_dist_vec(3,iA,i)
        riA = nucl_elec_dist(iA,i)
        arg = a * riA * riA
        tmp = a * c * dexp(-arg)
        tmpx = tmpx + tmp * dx
        tmpy = tmpy + tmp * dy
        tmpz = tmpz + tmp * dz
        tmpl = tmpl + tmp * (3.d0 - 2.d0 * arg)
      enddo
      deriv_env_x(i) = 2.d0 * tmpx
      deriv_env_y(i) = 2.d0 * tmpy
      deriv_env_z(i) = 2.d0 * tmpz
      lapl_env   (i) = 2.d0 * tmpl
    enddo

  elseif(env_type .eq. "Sum_Quartic") then

    do i = 1, elec_num
      tmpx = 0.d0
      tmpy = 0.d0
      tmpz = 0.d0
      tmpl = 0.d0
      !DIR$ LOOP COUNT (100)
      do iA = 1, nucl_num
        a = env_expo(iA)
        c = env_coef(iA)
        ! xi - xA = nucl_elec_dist_vec(1,iA,i)
        dx  = nucl_elec_dist_vec(1,iA,i)
        dy  = nucl_elec_dist_vec(2,iA,i)
        dz  = nucl_elec_dist_vec(3,iA,i)
        riA = nucl_elec_dist(iA,i)
        r2  = riA * riA
        r4  = r2  * r2
        tmp = a * c * r2 * dexp(-a*r4)
        tmpx = tmpx + tmp * dx
        tmpy = tmpy + tmp * dy
        tmpz = tmpz + tmp * dz
        tmpl = tmpl + tmp * (5.d0 - 4.d0*a*r4)
      enddo
      deriv_env_x(i) = 4.d0 * tmpx
      deriv_env_y(i) = 4.d0 * tmpy
      deriv_env_z(i) = 4.d0 * tmpz
      lapl_env   (i) = 4.d0 * tmpl
    enddo

  else

    print *, ' Error in deriv_env & lapl_env: Unknown env_type = ', env_type
    stop

  endif

END_PROVIDER

! ---

