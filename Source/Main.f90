!***************************************************************************************
!*                              PROGRAM scarpia
!***************************************************************************************
!>  \mainpage      Program scarpia
!>
!>  Optimization of a potential energy surface with gradient computed by \n
!>  finite difference and with a first order optimization algorithm \n
!>
!>  \author        Matteo Bonfanti, Robert Binder
!>  \version       VERSIONTAG
!>  \date          27 November 2017
!>
!***************************************************************************************
PROGRAM scarpia
#include "preprocessoptions.cpp"
   USE FiniteDifference
   USE Optimize

   IMPLICIT NONE

   ! variable definition for preprocessing directives
   __TIME_STRING__

   ! Variable to handle the command line
   INTEGER :: NArgs
   LOGICAL :: Help = .FALSE.

   ! Input file name, set from command line arguments
   CHARACTER(120) :: InputFileName

   ! Number of coordinates of the potential
   INTEGER :: NCoord = 1

   ! Guess geometry of the minimum
   REAL, DIMENSION(:), ALLOCATABLE :: XGuess
   ! Optimized geometry of the minimum
   REAL, DIMENSION(:), ALLOCATABLE :: XOpt

   ! Finite difference step
   REAL, PARAMETER :: DeltaFiniteDiff = 0.01


   !*************************************************************
   !   INITIAL MESSAGES AND MISCELLANOUS STUFF
   !*************************************************************

   PRINT "(/,     '                    ==============================')"
   PRINT "(       '                               scarpia          ')"
   PRINT "(       '                    ==============================',/)"
   PRINT "(       '                       Author: M. Bonfanti, R. Binder'      )"
   PRINT "(       '                       Release: ',A)", VERSIONTAG
   PRINT "(       '                       Compilation: ',A,1X,A,/)", __DATE__, __TIME__

   PRINT "(       '                       << Fu scoperta la fuga! ')"
   PRINT "(       '                  Or Scarpia i suoi sbirri sguinzaglia! >> '/)"
   PRINT "(       '               [G.Giacosa & L.Illica, Tosca (G. Puccini), Atto I]  ',2/)"

#if defined(LOG_FILE)
   __INIT_LOG_FILE
#endif

   !*************************************************************
   !         COMMAND LINE ARGUMENT
   !*************************************************************

   ! Check and read from command line the input file name
   NArgs = COMMAND_ARGUMENT_COUNT()
   IF (NArgs<1) THEN
      Help = .TRUE.
   ELSE
      CALL GET_COMMAND_ARGUMENT( 1, InputFileName )
      IF ( trim(InputFileName) == "help" ) Help = .TRUE.
   ENDIF
   IF (Help) THEN ! Call help
      PRINT*, ' Launch this program as:'
      PRINT*, ' % scarpia "InputFileName" '
      STOP
   ENDIF

   !*************************************************************
   !                 INPUT AND ARRAY ALLOCATION
   !*************************************************************

   ! HERE INPUT DATA

   ! Allocate X vectors
   ALLOCATE( XGuess(NCoord), XOpt(NCoord) )
   ! Define guess geometry
   XGuess = 0.0

   !*************************************************************
   !                 OPTIMIZATION
   !*************************************************************

   XOpt = SteepLocator( GetPotAndForces, XGuess, 1000, 1.E-5, DeltaFiniteDiff )

   !*************************************************************
   !                 OUTPUT AND DEALLOCATION
   !*************************************************************

   ! Print result


   ! Deallocate X vectors
   DEALLOCATE( XGuess, XOpt )


   CONTAINS


      !*************************************************************
      ! Interface to the routine that computes the potential
      !*************************************************************

      REAL FUNCTION GetPotential( X )
         REAL, DIMENSION(:), INTENT(IN)  :: X

         ! DEFINE HERE THE POTENTIAL AS A FUNCTION OF THE COORDINATES X

         GetPotential = 0.0

      END FUNCTION GetPotential


      !*************************************************************
      ! Interface to the finite difference estimate of the gradient of the potential
      !*************************************************************

      REAL FUNCTION GetPotAndForces( X, Force )
         REAL, DIMENSION(:), INTENT(IN)  :: X
         REAL, DIMENSION(:), INTENT(OUT) :: Force

         GetPotAndForces =  GetPotential( X )
         Force = GetGradient( X, GetPotential, DeltaFiniteDiff )

      END FUNCTION GetPotAndForces

END PROGRAM