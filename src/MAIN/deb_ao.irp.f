program deb_ao

  implicit none
  integer :: i, ipoint, npoint
  real    :: r(3), ao_val, ao_der(3), ao_lap

  print*, ' testing AOs of QMC=CHEM'

!  r(1) = -0.001483754835917605
!  r(2) = -0.0004945849453058681
!  r(3) = -0.132591278813387     
!
!  print*, ' VAL, GRAD, LAPL of AOs on r ='
!  print*, r

  npoint = 4122
  open(unit=11, name="grid", action="read")
    do ipoint = 1, npoint
      read(11, *) r(1), r(2), r(3)
      do i = 1, ao_num
        call get_ao_val_der_lap(i, r, ao_val, ao_der, ao_lap)
        write(*, '(5(f15.7, 3X))') ao_val, ao_der, ao_lap
      enddo
    enddo
  close(11)

end

