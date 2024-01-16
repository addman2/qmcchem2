! j1e
! ----

! ---

BEGIN_PROVIDER [double precision, jast_elec_1e_value, (elec_num_8)]

  include '../constants.F'

  implicit none
  integer          :: i, p, iA
  double precision :: c, a, dx, dy, dz, riA, riA2, tmp

  PROVIDE j1e_type

  if(j1e_type .eq. "none") then

    jast_elec_1e_value = 0.d0

  elseif(j1e_type .eq. "gauss") then

    ! \sum_{A} \sum_p c_{p_A} \exp(-\alpha_{p_A} (r - R_A)^2)

    do i = 1, elec_num

      tmp = 0.d0
      do iA = 1, nucl_num

        dx   = nucl_elec_dist_vec(1,iA,i)
        dy   = nucl_elec_dist_vec(2,iA,i)
        dz   = nucl_elec_dist_vec(3,iA,i)
        riA  = nucl_elec_dist(iA,i)
        riA2 = riA * riA

        do p = 1, j1e_size
          a = j1e_expo(p,iA)
          c = j1e_coef(p,iA)
          tmp = tmp + c * dexp(-a*riA2)
        enddo
      enddo

      jast_elec_1e_value(i) = tmp
    enddo

  elseif(j1e_type .eq. "charge-harmonizer") then

    ! TODO
    jast_elec_1e_value = 0.d0

  else

    print *, ' Error in jast_elec_1e_value,: Unknown j1e_type = ', j1e_type
    stop

  endif ! j1e_type

END_PROVIDER

! ---

 BEGIN_PROVIDER [double precision, jast_elec_1e_grad_x, (elec_num_8)]
&BEGIN_PROVIDER [double precision, jast_elec_1e_grad_y, (elec_num_8)]
&BEGIN_PROVIDER [double precision, jast_elec_1e_grad_z, (elec_num_8)]
&BEGIN_PROVIDER [double precision, jast_elec_1e_lapl  , (elec_num_8)]

  include '../constants.F'

  implicit none
  integer          :: i, p, iA
  double precision :: c, a, dx, dy, dz, riA, riA2, tmp
  double precision :: tmp_x, tmp_y, tmp_z, g

  PROVIDE j1e_type

  if(j1e_type .eq. "none") then

    jast_elec_1e_grad_x = 0.d0 
    jast_elec_1e_grad_y = 0.d0
    jast_elec_1e_grad_z = 0.d0
    jast_elec_1e_lapl   = 0.d0

  elseif(j1e_type .eq. "gauss") then

    ! \sum_{A} \sum_p c_{p_A} \exp(-\alpha_{p_A} (r - R_A)^2)

    do i = 1, elec_num

      tmp_x = 0.d0
      tmp_y = 0.d0
      tmp_z = 0.d0
      tmp   = 0.d0
      do iA = 1, nucl_num

        dx   = nucl_elec_dist_vec(1,iA,i)
        dy   = nucl_elec_dist_vec(2,iA,i)
        dz   = nucl_elec_dist_vec(3,iA,i)
        riA  = nucl_elec_dist(iA,i)
        riA2 = riA * riA

        do p = 1, j1e_size
          a = j1e_expo(p,iA)
          c = j1e_coef(p,iA)
          g = c * a * dexp(-a*riA2)

          tmp_x = tmp_x + g * dx
          tmp_y = tmp_y + g * dy
          tmp_z = tmp_z + g * dz
          tmp   = tmp   + g * (-3.d0 + 2.d0*a*riA2)
        enddo
      enddo

      jast_elec_1e_grad_x(i) = -2.d0*tmp_x 
      jast_elec_1e_grad_y(i) = -2.d0*tmp_y 
      jast_elec_1e_grad_z(i) = -2.d0*tmp_z 
      jast_elec_1e_lapl  (i) = tmp 
    enddo

  elseif(j1e_type .eq. "charge-harmonizer") then

    ! TODO
    jast_elec_1e_grad_x = 0.d0 
    jast_elec_1e_grad_y = 0.d0
    jast_elec_1e_grad_z = 0.d0
    jast_elec_1e_lapl   = 0.d0

  else

    print *, ' Error in jast_elec_1e_grad & jast_elec_1e_lapl,: Unknown j1e_type = ', j1e_type
    stop

  endif ! j1e_type

END_PROVIDER

! ---

