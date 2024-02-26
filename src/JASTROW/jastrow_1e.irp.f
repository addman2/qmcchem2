! j1e
! ----

 BEGIN_PROVIDER [double precision, jast_elec_1e_value , (elec_num_8)]
&BEGIN_PROVIDER [double precision, jast_elec_1e_grad_x, (elec_num_8)]
&BEGIN_PROVIDER [double precision, jast_elec_1e_grad_y, (elec_num_8)]
&BEGIN_PROVIDER [double precision, jast_elec_1e_grad_z, (elec_num_8)]
&BEGIN_PROVIDER [double precision, jast_elec_1e_lapl  , (elec_num_8)]

  include '../constants.F'

  implicit none
  integer          :: i, j, k, p, iA
  real             :: r(3), ao_val, ao_der(3), ao_lap
  double precision :: c, ckj, a, dx, dy, dz, riA, riA2, tmp1, tmp2
  double precision :: tmp_x, tmp_y, tmp_z, g, f
  double precision :: aok_v, aok_gx, aok_gy, aok_gz, aok_l
  double precision :: aoj_v, aoj_gx, aoj_gy, aoj_gz, aoj_l


  PROVIDE j1e_type

  if(j1e_type .eq. "None") then

    jast_elec_1e_value  = 0.d0
    jast_elec_1e_grad_x = 0.d0 
    jast_elec_1e_grad_y = 0.d0
    jast_elec_1e_grad_z = 0.d0
    jast_elec_1e_lapl   = 0.d0

  elseif(j1e_type .eq. "Gauss") then

    ! \sum_{A} \sum_p c_{p_A} \exp(-\alpha_{p_A} (r - R_A)^2)

    do i = 1, elec_num

      tmp1  = 0.d0
      tmp_x = 0.d0
      tmp_y = 0.d0
      tmp_z = 0.d0
      tmp2  = 0.d0
      do iA = 1, nucl_num

        dx   = nucl_elec_dist_vec(1,iA,i)
        dy   = nucl_elec_dist_vec(2,iA,i)
        dz   = nucl_elec_dist_vec(3,iA,i)
        riA  = nucl_elec_dist(iA,i)
        riA2 = riA * riA

        do p = 1, j1e_size
          a = j1e_expo(p,iA)
          c = j1e_coef(p,iA)
          f = c * dexp(-a*riA2)
          g = a * f

          tmp1  = tmp1  + f
          tmp_x = tmp_x + g * dx
          tmp_y = tmp_y + g * dy
          tmp_z = tmp_z + g * dz
          tmp2  = tmp2  + g * (3.d0 - 2.d0*a*riA2)
        enddo
      enddo

      jast_elec_1e_value (i) = tmp1
      jast_elec_1e_grad_x(i) = -2.d0 * tmp_x 
      jast_elec_1e_grad_y(i) = -2.d0 * tmp_y 
      jast_elec_1e_grad_z(i) = -2.d0 * tmp_z 
      jast_elec_1e_lapl  (i) = -2.d0 * tmp2
    enddo

  elseif(j1e_type .eq. "Charge_Harmonizer_AO") then

    PROVIDE qmckl_ao_vgl
    do i = 1, elec_num

      tmp1  = 0.d0
      tmp_x = 0.d0
      tmp_y = 0.d0
      tmp_z = 0.d0
      tmp2  = 0.d0
      do k = 1, ao_num
        aok_v  = qmckl_ao_vgl(k,1,i)
        aok_gx = qmckl_ao_vgl(k,2,i)
        aok_gy = qmckl_ao_vgl(k,3,i)
        aok_gz = qmckl_ao_vgl(k,4,i)
        aok_l  = qmckl_ao_vgl(k,5,i)

        do j = 1, ao_num
          aoj_v  = qmckl_ao_vgl(j,1,i)
          aoj_gx = qmckl_ao_vgl(j,2,i)
          aoj_gy = qmckl_ao_vgl(j,3,i)
          aoj_gz = qmckl_ao_vgl(j,4,i)
          aoj_l  = qmckl_ao_vgl(j,5,i)

          ckj = j1e_coef_ao2(j,k)

          tmp1  = tmp1  + ckj *  aok_v * aoj_v
          tmp_x = tmp_x + ckj * (aok_v * aoj_gx + aok_gx * aoj_v)
          tmp_y = tmp_y + ckj * (aok_v * aoj_gy + aok_gy * aoj_v)
          tmp_z = tmp_z + ckj * (aok_v * aoj_gz + aok_gz * aoj_v)
          tmp2  = tmp2  + ckj * (aok_v * aoj_l  + aok_l  * aoj_v + 2.d0 * (aok_gx*aoj_gx + aok_gy*aoj_gy + aok_gz*aoj_gz))
        enddo
      enddo

      jast_elec_1e_value (i) = tmp1
      jast_elec_1e_grad_x(i) = tmp_x 
      jast_elec_1e_grad_y(i) = tmp_y 
      jast_elec_1e_grad_z(i) = tmp_z 
      jast_elec_1e_lapl  (i) = tmp2
    enddo


  else

    print *, ' Error in jast_elec_1e: Unknown j1e_type = ', j1e_type
    stop

  endif ! j1e_type

END_PROVIDER

! ---




