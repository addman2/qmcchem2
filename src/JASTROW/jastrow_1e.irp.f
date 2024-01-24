! j1e
! ----

 BEGIN_PROVIDER [double precision, jast_elec_1e_value , (elec_num_8)]
&BEGIN_PROVIDER [double precision, jast_elec_1e_grad_x, (elec_num_8)]
&BEGIN_PROVIDER [double precision, jast_elec_1e_grad_y, (elec_num_8)]
&BEGIN_PROVIDER [double precision, jast_elec_1e_grad_z, (elec_num_8)]
&BEGIN_PROVIDER [double precision, jast_elec_1e_lapl  , (elec_num_8)]

  include '../constants.F'

  implicit none
  integer          :: i, j, p, iA
  real             :: r(3), ao_val, ao_der(3), ao_lap
  double precision :: c, a, dx, dy, dz, riA, riA2, tmp1, tmp2
  double precision :: tmp_x, tmp_y, tmp_z, g, f

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

    do i = 1, elec_num

      r(1) = elec_coord_transp(1,i)
      r(2) = elec_coord_transp(2,i)
      r(3) = elec_coord_transp(3,i)

      tmp1  = 0.d0
      tmp_x = 0.d0
      tmp_y = 0.d0
      tmp_z = 0.d0
      tmp2  = 0.d0
      do j = 1, ao_num

        c = j1e_coef_ao(j)

        call get_ao_val_der_lap(j, r, ao_val, ao_der, ao_lap)

        tmp1  = tmp1  + c * dble(ao_val   )  
        tmp_x = tmp_x + c * dble(ao_der(1))
        tmp_y = tmp_y + c * dble(ao_der(2))
        tmp_z = tmp_z + c * dble(ao_der(3))
        tmp2  = tmp2  + c * dble(ao_lap   )
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




