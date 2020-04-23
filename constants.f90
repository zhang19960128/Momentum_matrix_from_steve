MODULE constants
  !
  IMPLICIT NONE
  !
  ! kind for default double-precision real/complex subtypes
  INTEGER, PARAMETER :: dp = kind(1.0d0)
  INTEGER, PARAMETER :: dpc = kind((1.0_dp,1.0_dp))   !Maybe should not be used
  !
  ! Units and conversions
  !
  REAL(dp), PARAMETER :: bohr_meter = 5.291772d-11
  REAL(dp), PARAMETER :: hartree_ev = 27.21138386d0
  REAL(dp), PARAMETER :: hartree_joule = 4.35975d-18
  !
  ! Physical constants
  !
  REAL(dp), PARAMETER :: e_mass = 9.109383d-31
  REAL(dp), PARAMETER :: e_charge=1.6021766d-19
  REAL(dp), PARAMETER :: h_ev = 4.135667d-15
  REAL(dp), PARAMETER :: h_joule = 6.626068d-34
  REAL(dp), PARAMETER :: h_hartree = 1.51983d-16
  REAL(dp), PARAMETER :: eps_0 = 8.845d-12
  REAL(dp), PARAMETER :: c_e = 2.999d8
  !
  ! Numbers
  !
  REAL(dp), PARAMETER :: pi = 3.141592653589793238462643383279502884197d0
  REAL(dp), PARAMETER :: eq_tol = 1.0d-6
  !
  COMPLEX(dp), PARAMETER :: czero = cmplx(0.d0, 0.d0, kind=dpc)
  COMPLEX(dp), PARAMETER :: cone = cmplx(1.d0, 0.d0, kind=dpc)
  COMPLEX(dp), PARAMETER :: ci = cmplx(0.d0, 1.d0, kind=dpc)
  !
END MODULE constants
