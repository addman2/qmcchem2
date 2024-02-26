! Boys_Handy-Handy's Jastrow 
! ---------------------

! ---

 BEGIN_PROVIDER [double precision, jast_Boys_Handy_value    ]
&BEGIN_PROVIDER [double precision, jast_Boys_Handy_value_inv]

  implicit none
  integer          :: i
  double precision :: argexpo

  argexpo = 0.d0
  !DIR$ VECTOR ALIGNED
  !DIR$ LOOP COUNT (200)
  do i = 1, elec_num
    argexpo += jast_elec_Boys_Handy_value(i)
  enddo

  jast_Boys_Handy_value     = dexp(argexpo)
  jast_Boys_Handy_value_inv = 1.d0 / jast_Boys_Handy_value

END_PROVIDER

! ---

BEGIN_PROVIDER [double precision, jast_elec_Boys_Handy_value, (elec_num_8)]

  BEGIN_DOC  
  !
  ! J(i) = 0.5 \sum_{j!=i} u_{BH}(ri,rj) 
  ! 
  END_DOC

  include '../constants.F'

  implicit none
  integer          :: i, j
  double precision :: tmp_ij, tmp_val(elec_num)

  PROVIDE env_type
  PROVIDE j1e_type

  if(env_type .ne. "None") then

    print*, ' not implemented yet'
    stop

  else

    do i = 1, elec_num

      jast_elec_Boys_Handy_value(i) = 0.d0
      call get_jBH_val_r1(i, elec_num, tmp_val)

      tmp_ij = 0.d0
      !DIR$ LOOP COUNT(100)
      do j = 1, elec_num
        if(j==i) cycle

        jast_elec_Boys_Handy_value(i) = jast_elec_Boys_Handy_value(i) + 0.5d0 *  tmp_val(j)
      enddo
    enddo

  endif ! env_type

  if(j1e_type .ne. "None") then
    do i = 1, elec_num
      jast_elec_Boys_Handy_value(i) = jast_elec_Boys_Handy_value(i) + jast_elec_1e_value(i)
    enddo
  endif

END_PROVIDER

! ---

 BEGIN_PROVIDER [double precision, jast_elec_Boys_Handy_grad_x, (elec_num_8)]
&BEGIN_PROVIDER [double precision, jast_elec_Boys_Handy_grad_y, (elec_num_8)]
&BEGIN_PROVIDER [double precision, jast_elec_Boys_Handy_grad_z, (elec_num_8)]
&BEGIN_PROVIDER [double precision, jast_elec_Boys_Handy_lapl  , (elec_num_8)]

  BEGIN_DOC  
  !
  ! \grad_i   J = \sum_{j!=i} \grad_i   u_{sym}(ri, rj)
  ! \grad_i^2 J = \sum_{j!=i} \grad_i^2 u_{sym}(ri, rj)
  !
  END_DOC

  include '../constants.F'

  implicit none
  integer          :: i, j
  double precision :: jderx(elec_num), jdery(elec_num), jderz(elec_num), jlapl(elec_num)

  PROVIDE j1e_type
  PROVIDE env_type

  if(env_type .eq. "None") then

    do i = 1, elec_num

      jast_elec_Boys_Handy_grad_x(i) = 0.d0
      jast_elec_Boys_Handy_grad_y(i) = 0.d0
      jast_elec_Boys_Handy_grad_z(i) = 0.d0
      jast_elec_Boys_Handy_lapl  (i) = 0.d0

      call get_jBH_der_lap_r1(i, elec_num, jderx, jdery, jderz, jlapl)

      !DIR$ LOOP COUNT (100)
      do j = 1, elec_num
        if(i==j) cycle

        jast_elec_Boys_Handy_grad_x(i) = jast_elec_Boys_Handy_grad_x(i) + jderx(j) 
        jast_elec_Boys_Handy_grad_y(i) = jast_elec_Boys_Handy_grad_y(i) + jdery(j) 
        jast_elec_Boys_Handy_grad_z(i) = jast_elec_Boys_Handy_grad_z(i) + jderz(j) 
        jast_elec_Boys_Handy_lapl  (i) = jast_elec_Boys_Handy_lapl  (i) + jlapl(j) 
      enddo
    enddo

  else

    print*, ' not implemented yet'
    stop

  endif ! env_type

  if(j1e_type .ne. "None") then
    do i = 1, elec_num
      jast_elec_Boys_Handy_grad_x(i) = jast_elec_Boys_Handy_grad_x(i) + jast_elec_1e_grad_x(i)
      jast_elec_Boys_Handy_grad_y(i) = jast_elec_Boys_Handy_grad_y(i) + jast_elec_1e_grad_y(i)
      jast_elec_Boys_Handy_grad_z(i) = jast_elec_Boys_Handy_grad_z(i) + jast_elec_1e_grad_z(i)
      jast_elec_Boys_Handy_lapl  (i) = jast_elec_Boys_Handy_lapl  (i) + jast_elec_1e_lapl  (i)
    enddo
  endif

END_PROVIDER

! ---

