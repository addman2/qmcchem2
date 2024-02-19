! Input data
! ----------

BEGIN_PROVIDER [logical, do_jast]

  BEGIN_DOC  
  ! If true, compute the Jastrow factor
  END_DOC

  implicit none
  include '../types.F'
  do_jast = j2e_type /= t_None
  call linfo(irp_here, 'do_jast', do_jast)

END_PROVIDER

! ---

BEGIN_PROVIDER [logical, do_jpsi]

  BEGIN_DOC  
  ! If true, add Jastrow factor to WF
  END_DOC

  implicit none
  include '../types.F'
  do_jpsi = jpsi_type /= t_None
  call linfo(irp_here, 'do_jpsi', do_jpsi)

END_PROVIDER

! ---

BEGIN_PROVIDER [integer, j2e_type]

  BEGIN_DOC
  ! Type of Jastrow factor : Simple or Core
  END_DOC

  include '../types.F'

  implicit none
  character*(32) :: buffer

  buffer = types(t_Simple)
  j2e_type = t_None
  call get_jastrow_j2e_type(buffer)
  if (buffer == types(t_Simple)) then
    j2e_type = t_Simple
  else if (buffer == types(t_None)) then
    j2e_type = t_None
  else if (buffer == types(t_Core)) then
    j2e_type = t_Core
  else if (buffer == types(t_Mu)) then
    j2e_type = t_Mu
  else if (buffer == types(t_Mu_Nu)) then
    j2e_type = t_Mu_Nu
  else if (buffer == types(t_Mur)) then
    j2e_type = t_Mur
    print*, ' do not forget to increase the block time'
  else if (buffer == types(t_Qmckl)) then
    j2e_type = t_Qmckl
  else if (buffer == types(t_Boys)) then
    j2e_type = t_Boys
  else if (buffer == types(t_Boys_Handy)) then
    j2e_type = t_Boys_Handy
  else
    call abrt(irp_here, 'j2e_type should be (None|Simple|Core|Mu|Mu_r|Qmckl|Boys|Boys_Handy)')
  endif
  call cinfo(irp_here, 'j2e_type', buffer)

END_PROVIDER

! ---

BEGIN_PROVIDER [integer, jpsi_type]

  BEGIN_DOC
  ! Type of Jastrow factor used as function not to sample
  END_DOC

  include '../types.F'

  implicit none
  character*(32) :: buffer

  buffer = types(t_Simple)
  jpsi_type = t_Core
  call get_jastrow_jpsi_type(buffer)
  if (buffer == types(t_Simple)) then
    jpsi_type = t_Simple
  else if (buffer == types(t_None)) then
    jpsi_type = t_None
  else if (buffer == types(t_Core)) then
    jpsi_type = t_Core
  else if (buffer == types(t_Mu)) then
    jpsi_type = t_Mu
  else if (buffer == types(t_Mu_Nu)) then
    jpsi_type = t_Mu_Nu
  else if (buffer == types(t_Mur)) then
    jpsi_type = t_Mur
    print*, ' do not forget to increase the block time'
  else if (buffer == types(t_Qmckl)) then
    jpsi_type = t_Qmckl
  else if (buffer == types(t_Boys)) then
    jpsi_type = t_Boys
  else if (buffer == types(t_Boys_Handy)) then
    jpsi_type = t_Boys_Handy
  else
    call abrt(irp_here, 'jpsi type should be (None|Simple|Core|Mu|Mu_Nu|Mu_r|Qmckl|Boys|Boys_Handy)')
  endif
  call cinfo(irp_here, 'jpsi_type',buffer)

END_PROVIDER

! ---

!BEGIN_PROVIDER [character*(32), j2e_type]
!  implicit none
!  include '../types.F'
!  j2e_type = "None" ! no 2s-Jastrow
!  call get_jastrow_j2e_type(j2e_type)
!END_PROVIDER

BEGIN_PROVIDER [real, mu_erf]
  implicit none
  mu_erf = 0.5
  call get_hamiltonian_mu_erf(mu_erf)
END_PROVIDER

BEGIN_PROVIDER [real, nu_erf]
  implicit none
  nu_erf = 1.0
  call get_jastrow_nu_erf(nu_erf)
END_PROVIDER

BEGIN_PROVIDER [real, a_boys]
  implicit none
  a_boys = 1.0
  call get_jastrow_a_boys(a_boys)
END_PROVIDER

! ---

BEGIN_PROVIDER [integer, mur_type]
  implicit none
  mur_type = 1
  call get_jastrow_mur_type(mur_type)
END_PROVIDER

BEGIN_PROVIDER [real, mu_r_ct]
  implicit none
  mu_r_ct = 0.6203504908994001
  call get_jastrow_mu_r_ct(mu_r_ct)
END_PROVIDER

! ---

BEGIN_PROVIDER [character*(32), env_type]
  implicit none
  include '../types.F'
  env_type = "None" ! no Envelop
  call get_jastrow_env_type(env_type)
END_PROVIDER

BEGIN_PROVIDER [real, env_expo, (nucl_num)]
  implicit none
  include '../types.F'
  env_expo(:) = 1.0
  call get_jastrow_env_expo(env_expo)
END_PROVIDER

BEGIN_PROVIDER [real, env_coef, (nucl_num)]
  implicit none
  include '../types.F'
  env_coef(:) = 1.0
  call get_jastrow_env_coef(env_coef)
END_PROVIDER

! ---

BEGIN_PROVIDER [character*(32), j1e_type]
  implicit none
  include '../types.F'
  j1e_type = "None" ! no 1e-Jastrow
  call get_jastrow_j1e_type(j1e_type)
END_PROVIDER

BEGIN_PROVIDER [integer, j1e_size]
  implicit none
  include '../types.F'
  j1e_size = 1 ! nb of term per nuclei
  call get_jastrow_j1e_size(j1e_size)
END_PROVIDER

BEGIN_PROVIDER [real, j1e_coef, (j1e_size, nucl_num)]
  implicit none
  include '../types.F'
  j1e_coef = 0.0
  call get_jastrow_j1e_coef(j1e_coef)
END_PROVIDER

BEGIN_PROVIDER [real, j1e_expo, (j1e_size, nucl_num)]
  implicit none
  include '../types.F'
  j1e_expo = 1.0
  call get_jastrow_j1e_expo(j1e_expo)
END_PROVIDER

BEGIN_PROVIDER [real, j1e_coef_ao, (ao_num)]
  implicit none
  integer :: i
  include '../types.F'
  j1e_coef_ao = 0.0
  call get_jastrow_j1e_coef_ao(j1e_coef_ao)
END_PROVIDER

BEGIN_PROVIDER [real, j1e_coef_ao2, (ao_num,ao_num)]
  implicit none
  integer :: i
  include '../types.F'
  j1e_coef_ao2 = 0.0
  call get_jastrow_j1e_coef_ao2(j1e_coef_ao2)
END_PROVIDER

! ---

BEGIN_PROVIDER [integer, jbh_size]
  implicit none
  include '../types.F'
  jbh_size = 1 ! nb of term per nuclei
  call get_jastrow_jbh_size(jbh_size)
END_PROVIDER

BEGIN_PROVIDER [real, jbh_ee, (nucl_num)]
  implicit none
  include '../types.F'
  jbh_ee = 1.0
  call get_jastrow_jbh_ee(jbh_ee)
END_PROVIDER

BEGIN_PROVIDER [real, jbh_en, (nucl_num)]
  implicit none
  include '../types.F'
  jbh_en = 1.0
  call get_jastrow_jbh_en(jbh_en)
END_PROVIDER

BEGIN_PROVIDER [real, jbh_c, (jbh_size, nucl_num)]
  implicit none
  include '../types.F'
  jbh_c = 0.0
  call get_jastrow_jbh_c(jbh_c)
END_PROVIDER

BEGIN_PROVIDER [integer, jbh_m, (jbh_size, nucl_num)]
  implicit none
  include '../types.F'
  jbh_m = 0
  call get_jastrow_jbh_m(jbh_m)
END_PROVIDER

BEGIN_PROVIDER [integer, jbh_n, (jbh_size, nucl_num)]
  implicit none
  include '../types.F'
  jbh_n = 0
  call get_jastrow_jbh_n(jbh_n)
END_PROVIDER

BEGIN_PROVIDER [integer, jbh_o, (jbh_size, nucl_num)]
  implicit none
  include '../types.F'
  jbh_o = 0
  call get_jastrow_jbh_o(jbh_o)
END_PROVIDER

! ---

BEGIN_PROVIDER [real, jast_a_up_up]

  BEGIN_DOC
  ! a_{up up} parameters of the Jastrow
  END_DOC

  implicit none

  jast_a_up_up = 0.5
  call get_jastrow_jast_a_up_up(jast_a_up_up)

END_PROVIDER

BEGIN_PROVIDER [real, jast_a_up_dn]

  BEGIN_DOC
  ! a_{up dn} parameters of the Jastrow
  END_DOC

  implicit none

  jast_a_up_dn = 0.5
  call get_jastrow_jast_a_up_dn(jast_a_up_dn)

END_PROVIDER

BEGIN_PROVIDER [real, jast_b_up_up]

  BEGIN_DOC
  ! b_{up up} parameters of the Jastrow
  END_DOC

  implicit none

  jast_b_up_up = 1.
  call get_jastrow_jast_b_up_up(jast_b_up_up)

END_PROVIDER

BEGIN_PROVIDER [real, jast_b_up_dn]

  BEGIN_DOC
  ! b_{up dn} parameters of the Jastrow
  END_DOC

  implicit none

  jast_b_up_dn = 1.
  call get_jastrow_jast_b_up_dn(jast_b_up_dn)

END_PROVIDER

BEGIN_PROVIDER [real, jast_pen, (nucl_num)]

  BEGIN_DOC
  ! penetration parameters of the Jastrow
  END_DOC

  implicit none

  include '../types.F'
  jast_pen(:) = 0.5
  call get_jastrow_jast_pen(jast_pen)

END_PROVIDER

BEGIN_PROVIDER [real, jast_eeN_e_a, (nucl_num)]

  BEGIN_DOC
  ! a parameters of the electron-electron-Nucleus component of the Jastrow
  END_DOC

  implicit none

  include '../types.F'
  jast_eeN_e_a(:) = 0.5
  call get_jastrow_jast_eeN_e_a(jast_eeN_e_a)

END_PROVIDER

BEGIN_PROVIDER [real, jast_eeN_e_b, (nucl_num)]

  BEGIN_DOC
  ! b parameters of the electron-electron-Nucleus component of the Jastrow
  END_DOC

  implicit none

  include '../types.F'
  jast_eeN_e_b(:) = 1.0
  call get_jastrow_jast_eeN_e_b(jast_eeN_e_b)

END_PROVIDER

BEGIN_PROVIDER [real, jast_eeN_N, (nucl_num)]

  BEGIN_DOC
  ! penetration parameters of the electron-electron-nucleus component of the Jastrow
  END_DOC

  implicit none
  include '../types.F'
  integer :: i
  jast_eeN_N(:) = 100.0
  call get_jastrow_jast_eeN_N(jast_eeN_N)

END_PROVIDER

BEGIN_PROVIDER [real, jast_core_a1, (nucl_num)]

  BEGIN_DOC
  ! parameters of the core Jastrow
  END_DOC

  implicit none
  include '../types.F'
  integer :: i
  do i=1,nucl_num
    if (nucl_charge(i) > 0.) then
      jast_core_a1(i) = 0.6/nucl_charge(i)
    else
      jast_core_a1(i) = 0.0
    endif
  enddo
  call get_jastrow_jast_core_a1(jast_core_a1)

END_PROVIDER

BEGIN_PROVIDER [real, jast_core_b1, (nucl_num)]

  BEGIN_DOC
  ! parameters of the core Jastrow
  END_DOC

  implicit none
  include '../types.F'
  jast_core_b1(:) = max(1.,1. - 0.3 * nucl_charge(:))
  call get_jastrow_jast_core_b1(jast_core_b1)

END_PROVIDER

