
subroutine get_jBH_val_r1(i, e_num, jval)

  BEGIN_DOC  
  !
  ! u(i,j) = \sum_A \sum_{pA} D(m_pA) c_pA (f(r_iA)^m_pA f(r_jA)^n_pA + f(r_iA)^n_pA f(r_jA)^m_pA) g(r_ij)^o_pA
  !
  END_DOC

  include '../constants.F'

  implicit none
  integer,          intent(in)  :: i, e_num
  double precision, intent(out) :: jval(e_num)

  integer                       :: j, i_nucl, p
  integer                       :: mpA, npA, opA
  double precision              :: rn(3), r1(3), r2(3)
  double precision              :: f1A, f2A, g12, tmp

  PROVIDE jBH_en jBH_ee jBH_size jBH_m jBH_n jBH_o jBH_c

  r1(1) = elec_coord_transp(1,i)
  r1(2) = elec_coord_transp(2,i)
  r1(3) = elec_coord_transp(3,i)

  do j = 1, e_num

    jval(j) = 0.d0
    if(j == i) cycle

    r2(1) = elec_coord_transp(1,j)
    r2(2) = elec_coord_transp(2,j)
    r2(3) = elec_coord_transp(3,j)

    do i_nucl = 1, nucl_num

      rn(1) = nucl_coord(i_nucl,1)
      rn(2) = nucl_coord(i_nucl,2)
      rn(3) = nucl_coord(i_nucl,3)

      call jBH_elem_fct(dble(jBH_en(i_nucl)), r1, rn, f1A)
      call jBH_elem_fct(dble(jBH_en(i_nucl)), r2, rn, f2A)
      call jBH_elem_fct(dble(jBH_ee(i_nucl)), r1, r2, g12)

      do p = 1, jBH_size
        mpA = jBH_m(p,i_nucl)
        npA = jBH_n(p,i_nucl)
        opA = jBH_o(p,i_nucl)
        tmp = dble(jBH_c(p,i_nucl))
        if(mpA .eq. npA) then
          tmp = tmp * 0.5d0
        endif

        jval(j) = jval(j) + tmp * g12**dble(opA) * (f1A**dble(mpA) * f2A**dble(npA) + f1A**dble(npA) * f2A**dble(mpA))
      enddo ! p
    enddo ! i_nucl
  enddo ! j

  return
end

! ---

subroutine get_jBH_der_lap_r1(i, e_num, jderx, jdery, jderz, jlapl)

  BEGIN_DOC  
  !
  ! d/dxi u(i,j) = 
  ! d/dyi u(i,j) = 
  ! d/dzi u(i,j) = 
  !
  ! \Delta_i u(i,j) = 
  !
  END_DOC

  include '../constants.F'

  implicit none
  integer,          intent(in)  :: i, e_num
  double precision, intent(out) :: jderx(e_num), jdery(e_num), jderz(e_num), jlapl(e_num)

  integer                       :: j, i_nucl, p
  integer                       :: mpA, npA, opA
  double precision              :: rn(3), r1(3), r2(3)
  double precision              :: f1A, der1_f1A(3), lap1_f1A
  double precision              :: f2A, der2_f2A(3), lap2_f2A
  double precision              :: g12, der1_g12(3), lap1_g12
  double precision              :: der1_f1A_sq, der1_g12_sq, der1_f1A_g12
  double precision              :: tmp, tmp1, tmp2, tmp3


  PROVIDE jBH_en jBH_ee jBH_size jBH_m jBH_n jBH_o jBH_c

  r1(1) = elec_coord_transp(1,i)
  r1(2) = elec_coord_transp(2,i)
  r1(3) = elec_coord_transp(3,i)

  do j = 1, e_num

    jderx(j) = 0.d0
    jdery(j) = 0.d0
    jderz(j) = 0.d0
    jlapl(j) = 0.d0
    if(j == i) cycle

    r2(1) = elec_coord_transp(1,j)
    r2(2) = elec_coord_transp(2,j)
    r2(3) = elec_coord_transp(3,j)

    do i_nucl = 1, nucl_num

      rn(1) = nucl_coord(i_nucl,1)
      rn(2) = nucl_coord(i_nucl,2)
      rn(3) = nucl_coord(i_nucl,3)

      call jBH_elem_fct_der_lap(dble(jBH_en(i_nucl)), r1, rn, f1A, der1_f1A, lap1_f1A)
      call jBH_elem_fct_der_lap(dble(jBH_en(i_nucl)), r2, rn, f2A, der2_f2A, lap2_f2A)
      call jBH_elem_fct_der_lap(dble(jBH_ee(i_nucl)), r1, r2, g12, der1_g12, lap1_g12)

      der1_f1A_sq  = der1_f1A(1) * der1_f1A(1) &
                   + der1_f1A(2) * der1_f1A(2) &
                   + der1_f1A(3) * der1_f1A(3)
      der1_g12_sq  = der1_g12(1) * der1_g12(1) &
                   + der1_g12(2) * der1_g12(2) &
                   + der1_g12(3) * der1_g12(3)
      der1_f1A_g12 = der1_f1A(1) * der1_g12(1) &
                   + der1_f1A(2) * der1_g12(2) &
                   + der1_f1A(3) * der1_g12(3)

      do p = 1, jBH_size
        mpA = jBH_m(p,i_nucl)
        npA = jBH_n(p,i_nucl)
        opA = jBH_o(p,i_nucl)
        tmp = dble(jBH_c(p,i_nucl))
        if(mpA .eq. npA) then
          tmp = tmp * 0.5d0
        endif

        tmp1 = 0.d0
        tmp3 = 0.d0
        if(mpA .gt. 0) then
          tmp1 = tmp1 + dble(mpA) * f1A**dble(mpA-1) * f2A**dble(npA)
          if(mpA .gt. 1) then
            tmp3 = tmp3 + dble(mpA) * dble(mpA-1) * f1A**dble(mpA-2) * f2A**dble(npA)
          endif
        endif
        if(npA .gt. 0) then
          tmp1 = tmp1 + dble(npA) * f1A**dble(npA-1) * f2A**dble(mpA)
          if(npA .gt. 1) then
            tmp3 = tmp3 + dble(npA) * dble(npA-1) * f1A**dble(npA-2) * f2A**dble(mpA)
          endif
        endif
        tmp1 = tmp1 * g12**dble(opA)
        tmp3 = tmp3 * g12**dble(opA) * der1_f1A_sq

        tmp2 = 0.d0
        if(opA .gt. 0) then

          tmp2 = tmp2 + dble(opA) * g12**dble(opA-1) * (f1A**dble(mpA) * f2A**dble(npA) + f1A**dble(npA) * f2A**dble(mpA))

          if(opA .gt. 1) then
            tmp3 = tmp3 + dble(opA) * dble(opA-1) * g12**dble(opA-2) * (f1A**dble(mpA) * f2A**dble(npA) + f1A**dble(npA) * f2A**dble(mpA)) * der1_g12_sq
          endif

          if(mpA .gt. 0) then
            tmp3 = tmp3 + 2.d0 * dble(mpA) * f1A**dble(mpA-1) * f2A**dble(npA) * dble(opA) * g12**dble(opA-1) * der1_f1A_g12
          endif
          if(npA .gt. 0) then
            tmp3 = tmp3 + 2.d0 * dble(npA) * f1A**dble(npA-1) * f2A**dble(mpA) * dble(opA) * g12**dble(opA-1) * der1_f1A_g12
          endif
        endif

        jderx(j) = jderx(j) + tmp * (tmp1 * der1_f1A(1) + tmp2 * der1_g12(1))
        jdery(j) = jdery(j) + tmp * (tmp1 * der1_f1A(2) + tmp2 * der1_g12(2))
        jderz(j) = jderz(j) + tmp * (tmp1 * der1_f1A(3) + tmp2 * der1_g12(3))
        jlapl(j) = jlapl(j) + tmp * (tmp1 * lap1_f1A    + tmp2 * lap1_g12   + tmp3)
      enddo ! p
    enddo ! i_nucl
  enddo ! j

  return
end

! ---

subroutine jBH_elem_fct(alpha, r1, r2, fct)

  implicit none
  double precision, intent(in)  :: alpha, r1(3), r2(3)
  double precision, intent(out) :: fct
  double precision              :: dist, tmp1, tmp2

  dist = dsqrt( (r1(1) - r2(1)) * (r1(1) - r2(1)) &
              + (r1(2) - r2(2)) * (r1(2) - r2(2)) &
              + (r1(3) - r2(3)) * (r1(3) - r2(3)) )

  tmp1 = 1.d0 / (1.d0 + alpha * dist)

  fct = alpha * dist * tmp1

  return
end 

! ---

subroutine jBH_elem_fct_der_lap(alpha, r1, r2, fct, der1_fct, lap1_fct)

  implicit none
  double precision, intent(in)  :: alpha, r1(3), r2(3)
  double precision, intent(out) :: fct, der1_fct(3), lap1_fct
  double precision              :: dist, tmp1, tmp2

  dist = dsqrt( (r1(1) - r2(1)) * (r1(1) - r2(1)) &
              + (r1(2) - r2(2)) * (r1(2) - r2(2)) &
              + (r1(3) - r2(3)) * (r1(3) - r2(3)) )

  tmp1 = 1.d0 / (1.d0 + alpha * dist)

  fct = alpha * dist * tmp1

  lap1_fct = 2.d0 * alpha * tmp1 * tmp1 * tmp1 / dist

  if(dist .lt. 1d-10) then
    der1_fct(1) = 0.d0
    der1_fct(2) = 0.d0
    der1_fct(3) = 0.d0
  else
    tmp2 = alpha * tmp1 * tmp1 / dist
    der1_fct(1) = tmp2 * (r1(1) - r2(1))
    der1_fct(2) = tmp2 * (r1(2) - r2(2))
    der1_fct(3) = tmp2 * (r1(3) - r2(3))
  endif

  return
end 

! ---

