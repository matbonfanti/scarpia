!******************************************************************************
!*                      MODULE FiniteDifference
!******************************************************************************
!
!>  \brief     Compute derivatives with finite difference
!>  \details   Subroutines which compute derivatives by means of \n
!>             finite difference formulas: \n
!>             1) GetGradient compute gradient vector from potential
!>             2) GetHessian  compute hessian matrix from potential
!>             3) GetHessianFromForces computes Hessian matrix as
!>                  first derivatives of the analytic forces
!>             4) TestForces  test implemented analytic forces by
!>                  comparing them with finite difference ones
!
!******************************************************************************
!
!>  \author           Matteo Bonfanti
!>  \version          1.0
!>  \date             20 February 2017
!>
!******************************************************************************
!
!>  \par Updates
!>  \arg DATE : brief description of the change
!>
!>  \todo   introduce a initialization mechanism to define the type of
!>               finite difference formula to use
!>  \todo   include some automatic testing procedures which might be
!>               used to identify a proper interval for the displacement
!>               e.g. compute difference (choosing some appropriate matrix norm)
!>               of the hessian with respect to a reference hessian for
!>               difference values of smalldelta
!>
!******************************************************************************

MODULE FiniteDifference
#include "preprocessoptions.cpp"
!    USE RandomNumberGenerator
   IMPLICIT NONE

      PRIVATE
      PUBLIC :: GetGradient, GetHessian, GetHessianFromForces

      ! Finite difference - 4 points formula for first derivative
      !> Displacements in units of delta for 4pts first derivative formula
!       REAL, DIMENSION(4), PARAMETER :: DeltasI = (/ -2.0,    -1.0,    +1.0,    +2.0    /)
      !> Coefficients for 4pts first derivative formula
!       REAL, DIMENSION(4), PARAMETER :: CoeffsI = (/ +1./12., -8./12., +8./12., -1./12. /)

      ! Finite difference - 6 points formula for first derivative
      !> Displacements in units of delta for 6pts first derivative formula
      REAL, DIMENSION(6), PARAMETER :: DeltasI = (/ -3.0,    -2.0,    -1.0,    +1.0,     +2.0,    +3.0  /)
      !> Coefficients for 6pts first derivative formula
      REAL, DIMENSION(6), PARAMETER :: CoeffsI = (/ -1./60., 3./20.,  -3./4.,  3./4., -3./20.,   1./60. /)

      ! Finite difference - 3 points formula for first derivative at the right
      !> Displacements in units of delta for 3pts first right derivative formula
      REAL, DIMENSION(3), PARAMETER :: ForwardDeltasI = (/  0.0,   +1.0,  +2.0   /)
      !> Displacements in units of delta for 3pts first right derivative formula
      REAL, DIMENSION(3), PARAMETER :: ForwardCoeffsI = (/ -3./2., +2.0,  -1./2. /)

      ! Finite difference - 9 points formula for second derivative
      !> Displacements in units of delta for 9pts second derivative formula
      REAL, DIMENSION(9), PARAMETER :: DeltasII = &
        (/ -4.0    , -3.0   , -2.0  , -1.0 ,  0.0     , +1.0 , +2.0  , +3.0   , +4.0     /)
      !> Coefficients for 9pts second derivative formula
      REAL, DIMENSION(9), PARAMETER :: CoeffsII = &
        (/ -1./560., 8./315., -1./5., 8./5., -205./72., 8./5., -1./5., 8./315., -1./560. /)

      ! Default value for delta for numerical first derivatives
      REAL, PARAMETER :: SmallDeltaI = 0.01
      ! Default value for delta for numerical second derivatives
      REAL, PARAMETER :: SmallDeltaII = 0.005

      !> Number of random generated points for the derivative testing
      INTEGER, PARAMETER :: NPointsTest = 10.**3

#if defined(LOG_FILE)
    CHARACTER(17), SAVE :: LogStr = " FiniteDifference |"
#endif

!       ! use preprocession option to write debug information to log file...
! #if defined(LOG_FILE)
!       __OPEN_LOG_FILE
!       WRITE(__LOG_UNIT,*)  LogStr," Pippo Pluto Paperino"
!       __CLOSE_LOG_FILE
! #endif
!
!       ! ...or to print them as output
! #if defined(VERBOSE_OUTPUT)
!       WRITE(*,*) "Pippo pluto paperino"
! #endif

   CONTAINS

!==============================================================================
!                                 SUBROUTINES
!==============================================================================


!******************************************************************************
!> Compute the gradient of a potential given as external subroutine
!> When DeltaInp is not given, default value of the module is taken.
!>
!> @param AtPoint       Input vector with the coords where to compute gradient
!> @param GetPotential  Function to evaluate the potential
!> @param DeltaInp      Optional, magnitude of the finite coords displacements
!> @returns             Array with the first derivatives of the pot in AtPoint
!******************************************************************************
FUNCTION GetGradient( AtPoint, GetPotential, DeltaInp ) RESULT( Grad )
      IMPLICIT NONE
      REAL, DIMENSION(:), INTENT(IN)   :: AtPoint
      REAL, DIMENSION(size(AtPoint))   :: Grad
      REAL, OPTIONAL                   :: DeltaInp

      INTERFACE
         REAL FUNCTION GetPotential( X )
            REAL, DIMENSION(:), INTENT(IN)  :: X
         END FUNCTION GetPotential
      END INTERFACE

      REAL, DIMENSION(size(AtPoint)) :: Coordinates
      REAL :: Potential, SmallDelta
      INTEGER :: i,k

      ! Define the displacement length
      IF (PRESENT(DeltaInp)) THEN
         SmallDelta = DeltaInp
      ELSE
         SmallDelta = SmallDeltaI
      ENDIF

      ! Initialize gradient to 0
      Grad(:) = 0.0

      DO i = 1, size(AtPoint)       ! Cycle over the number of coordinates
         DO k = 1, size(DeltasI)         !  Cycle over the finite displacements

               ! Define small displacement from the point where compute the derivative
               Coordinates(:) = AtPoint(:)
               Coordinates(i) = Coordinates(i) + DeltasI(k)*SmallDelta

               ! Compute potential and forces in the displaced coordinate
               Potential = GetPotential( Coordinates )
               ! Increment numerical derivative of the analytical derivative
               Grad(i) = Grad(i) + CoeffsI(k)*Potential

         END DO
      END DO

      ! Divide by the finite displament length
      Grad(:) = Grad(:)/SmallDelta

   END FUNCTION GetGradient

!==============================================================================


!******************************************************************************
!> Compute the Hessian of a potential given as external subroutine
!> When DeltaInp is not given, default value of the module is taken.
!>
!> @param AtPoint       Input vector with the coords where to compute Hessian
!> @param GetPotential  Function to evaluate the potential
!> @param DeltaInp      Optional, magnitude of the finite coords displacements
!> @returns             Matrix of 2nd derivatives of the potential in AtPoint
!******************************************************************************
   FUNCTION GetHessian( AtPoint, GetPotential, DeltaInp ) RESULT( Hessian )
      REAL, DIMENSION(:), INTENT(IN)                 :: AtPoint
      REAL, DIMENSION(size(AtPoint),size(AtPoint))   :: Hessian
      REAL, OPTIONAL                                 :: DeltaInp

      INTERFACE
         REAL FUNCTION GetPotential( X )
            REAL, DIMENSION(:), INTENT(IN)  :: X
         END FUNCTION GetPotential
      END INTERFACE

      REAL, DIMENSION(size(AtPoint)) :: Coordinates
      REAL :: Potential, SmallDelta
      INTEGER :: i,j, m,n

      ! Define the displacement length
      IF (PRESENT(DeltaInp)) THEN
         SmallDelta = DeltaInp
      ELSE
         SmallDelta = SmallDeltaII
      ENDIF

      ! Initialize hessian to 0
      Hessian(:,:) = 0.0

      ! Diagonal elements
      DO i = 1, size(AtPoint)
         DO n = 1, size(DeltasII)
            ! Define small displacement from the point where compute the derivative
            Coordinates(:) = AtPoint(:)
            Coordinates(i) = Coordinates(i) + DeltasII(n)*SmallDelta
            ! Compute potential and forces in the displaced coordinate
            Potential = GetPotential( Coordinates )
            ! Increment numerical derivative
            Hessian(i,i) = Hessian(i,i) + CoeffsII(n)*Potential
         END DO
      END DO

      DO i = 2, size(AtPoint)
         DO j = 1, i-1
            ! Off-diagonal elements, lower triangle
            DO m = 1, size(DeltasI)
               DO n = 1, size(DeltasI)
                  ! Define small displacement from the point where compute the derivative
                  Coordinates(:) = AtPoint(:)
                  Coordinates(i) = Coordinates(i) + DeltasI(m)*SmallDelta
                  Coordinates(j) = Coordinates(j) + DeltasI(n)*SmallDelta
                  ! Compute potential and forces in the displaced coordinate
                  Potential = GetPotential( Coordinates )
                  ! Increment numerical derivative
                  Hessian(i,j) = Hessian(i,j) + CoeffsI(m)*CoeffsI(n)*Potential
               END DO
            END DO
            ! Off-diagonal elements, upper triangle
            Hessian(j,i) = Hessian(i,j)
         END DO
      END DO

      ! Divide by the squared finite displament length
      Hessian(:,:) = Hessian(:,:) / SmallDelta**2

   END FUNCTION GetHessian

!==============================================================================


!******************************************************************************
!> Compute the Hessian of the potential from its Forces (Forces = -Gradient)
!>
!> @param AtPoint       Input vector with the coords where to compute Hessian
!> @param GetPotential  Function to evaluate the potential and forces
!> @param DeltaInp      Optional, magnitude of the finite coords displacements
!> @param RightDerivInp Optional, compute right derivative instead of centered (logic mask)
!> @returns             Matrix of 2nd derivatives of the potential in AtPoint
!**************************************************************************************
   FUNCTION GetHessianFromForces( AtPoint, GetPotAndForces, DeltaInp, RightDerivInp ) RESULT( Hessian )
      REAL, DIMENSION(:), INTENT(IN)                 :: AtPoint
      REAL, DIMENSION(size(AtPoint),size(AtPoint))   :: Hessian
      REAL, OPTIONAL                                 :: DeltaInp
      LOGICAL, DIMENSION(size(AtPoint)), OPTIONAL    :: RightDerivInp

      INTERFACE
         REAL FUNCTION GetPotAndForces( X, Force )
            REAL, DIMENSION(:), INTENT(IN)  :: X
            REAL, DIMENSION(:), INTENT(OUT) :: Force
         END FUNCTION GetPotAndForces
      END INTERFACE

      REAL, DIMENSION(size(AtPoint)) :: Coordinates, FirstDerivative
      LOGICAL, DIMENSION(size(AtPoint))    :: RightDeriv
      REAL :: Potential, SmallDelta
      INTEGER :: i,k

      ! Define the displacement length
      IF (PRESENT(DeltaInp)) THEN
         SmallDelta = DeltaInp
      ELSE
         SmallDelta = SmallDeltaI
      ENDIF
      ! Define the right derivative logical mask
      IF (PRESENT(RightDerivInp)) THEN
         RightDeriv = RightDerivInp
      ELSE
         RightDeriv = .FALSE.
      ENDIF

      ! Initialize hessian to 0
      Hessian(:,:) = 0.0

      DO i = 1, size(AtPoint)       ! Cycle over the number of coordinates
         IF (RightDeriv(i)) THEN
            DO k = 1, size(ForwardDeltasI)         !  Cycle over the finite displacements

               ! Define small displacement from the point where compute the derivative
               Coordinates(:) = AtPoint(:)
               Coordinates(i) = Coordinates(i) + ForwardDeltasI(k)*SmallDelta

               ! Compute potential and forces in the displaced coordinate
               Potential = GetPotAndForces( Coordinates, FirstDerivative )
               FirstDerivative = - FirstDerivative

               ! Increment numerical derivative of the analytical derivative
               Hessian(i,:) = Hessian(i,:) + ForwardCoeffsI(k)*FirstDerivative(:)

            END DO
         ELSE IF (.NOT. RightDeriv(i)) THEN
            DO k = 1, size(DeltasI)         !  Cycle over the finite displacements

               ! Define small displacement from the point where compute the derivative
               Coordinates(:) = AtPoint(:)
               Coordinates(i) = Coordinates(i) + DeltasI(k)*SmallDelta

               ! Compute potential and forces in the displaced coordinate
               Potential = GetPotAndForces( Coordinates, FirstDerivative )
               FirstDerivative = - FirstDerivative

               ! Increment numerical derivative of the analytical derivative
               Hessian(i,:) = Hessian(i,:) + CoeffsI(k)*FirstDerivative(:)

            END DO
         ENDIF
      END DO

      ! Divide by the finite displament length
      Hessian(:,:) = Hessian(:,:)/SmallDelta
!       CALL TheOneWithMatrixPrintedLineAfterLine( SecDeriv )

   END FUNCTION GetHessianFromForces

!==============================================================================


!==============================================================================
!                              END OF MODULE
!==============================================================================
END MODULE FiniteDifference
