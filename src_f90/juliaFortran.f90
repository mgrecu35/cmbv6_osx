subroutine initGeophys(nmembI,nmfreqI,sysdN)
  use geophysEns
  integer :: nmembI,nmfreqI
  real :: sysdN
  call allocGeophys(9,61,9,nmembI,nmfreqI*nmembI*2)
  call setdNwIcJcL(sysdN,nmembI)
  print*,sysdN,nmembI,nmfreqI

end subroutine initGeophys
