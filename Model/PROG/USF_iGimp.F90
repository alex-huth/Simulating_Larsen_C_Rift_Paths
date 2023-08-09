FUNCTION getpassive(Model,nodenumber,VarIn) RESULT(VarOut)
  USE DefUtils
  implicit none
  !-----------------
  TYPE(Model_t) :: Model
  INTEGER :: nodenumber
  REAL(kind=dp) :: VarIn(1),VarOut

  VarOut=(VarIn(1))*(-1.0_dp)
End FUNCTION getpassive

!icerises,velocity
FUNCTION zerovelaticerises(Model,nodenumber,VarIn) RESULT(VarOut)
  USE DefUtils
  implicit none
  !-----------------
  TYPE(Model_t) :: Model
  INTEGER :: nodenumber
  REAL(kind=dp) :: VarIn(2),VarOut


  IF (VarIn(1) >= 0.0_dp) THEN
     VarOut = 0.0_dp
  Else
     VarOut = VarIn(2)
  END IF

End FUNCTION zerovelaticerises

!icerises,velocity
FUNCTION bctrackoneaticerises(Model,nodenumber,VarIn) RESULT(VarOut)
  USE DefUtils
  implicit none
  !-----------------
  TYPE(Model_t) :: Model
  INTEGER :: nodenumber
  REAL(kind=dp) :: VarIn(2),VarOut


  IF (VarIn(1) >= 0.0_dp) THEN
     VarOut = 1.0_dp
  Else
     VarOut = VarIn(2)
  END IF

End FUNCTION bctrackoneaticerises
