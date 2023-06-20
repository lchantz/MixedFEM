!>******************************************************************************
!>    The form
!>******************************************************************************
program CahnHilliard
   use, intrinsic :: iso_fortran_env
   USE BiCGStab_mod
   IMPLICIT NONE
   INCLUDE 'CHparameters.inc'

   !> Space coordinates
   REAL(real64) x(NNmax),y(NNmax)

   !> Time coordinates
   REAL(real64) t(Tmax)

   !> Local to global node numbering
   INTEGER nop(NELmax,NNref)

   !> Order parameter - c - of the Cahn Hilliard Equation
   REAL(real64) c(NNmax,Tmax)

   !> Auxiliary variable - w - the chemical potential
   REAL(real64) w(NNmax,Tmax)

   !> Newton solution vector
   REAL(real64) S(2*NNmax)

   !> System matrices
   !> STIFF - Stiffness
   !> MASS - Mass
   !> P - Non linear advective term from double well potential c^3-c
   REAL(real64) STIFF(NNmax,NNmax)
   REAL(real64) MASS(NNmax,NNmax)
   REAL(real64) P(NNmax,Tmax)

   !> Derivative of the non-linear term
   REAL(real64) DP(NNmax,NNmax,Tmax)

   INTEGER i,j,timepoint

   OPEN(unit=2,file='OrderParameter.dat')

   WRITE(*,'(25x,A)') 'Results c(n,t)'
   CALL SpaceDiscretization(x,y)
   CALL GlobalNodeNumbering(nop)
   CALL TimeDiscretization(t)
   CALL WeakFormulation(x,y,STIFF,MASS,P,nop)

   DP=0.
   DO timepoint=1,Tmax
      P(:,timepoint)=0.
      DP(:,:,timepoint)=0.

      CALL NewtonProcedure(x,y,nop,t,timepoint,c,w,S,STIFF,MASS,P,DP)
      DO i=1,NNmax

         WRITE(2,*) t(timepoint), i, c(i,timepoint)

      ENDDO
   ENDDO
   CLOSE(2)
END program CahnHilliard


 !>******************************************************************************
 !> SpaceDiscretization creates the discrete points in space.
 !> -------------
 !> Params:
 !> -------------
 !> x,y:       space coordinates
 !> N:       number of space coordinates
 !> NEL:     number of discrete elements
 !>******************************************************************************
SUBROUTINE SpaceDiscretization(x,y)
   use, intrinsic :: iso_fortran_env
   IMPLICIT NONE
   INCLUDE 'CHparameters.inc'

   !> Space coordinates
   REAL(real64) x(NNmax),y(NNmax)

   REAL(real64) xorigin,yorigin,xlast,ylast
   REAL(real64) deltax,deltay
   INTEGER i,j,NNode

   xorigin=1.
   xlast=5.
   deltax=(xlast-xorigin)/2*NELx

   yorigin=1.
   ylast=5.
   deltay=(ylast-yorigin)/2*NELy

   !> Setting the coordinates of each node
   x(1)=xorigin
   y(1)=yorigin
   DO i=1,NNx
      NNode = (i-1)*NNy + 1 ! Ascending node numbering from bottom left up to right
      x(NNode) = x(1) + (i-1)*deltax
      y(NNode) = y(1)
      DO j=2,NNy
         x(NNode+j-1)=x(NNode)
         y(NNode+j-1)=y(NNode)+(j-1)*deltay

      ENDDO
   ENDDO



END SUBROUTINE SpaceDiscretization

 !>******************************************************************************
 !> SetGlobalSpaceCoordinates creates the mapping between the local coordinates
 !> of a point on the element to the global coordinates in space.
 !> -------------
 !> Params:
 !> -------------
 !> ngl:   global coordinates in space
 !> NEL:                 number of discrete elements
 !>******************************************************************************
SUBROUTINE GlobalNodeNumbering(nop)
   use, intrinsic :: iso_fortran_env
   IMPLICIT NONE
   INCLUDE 'CHparameters.inc'

   !> Local to global space numbering
   INTEGER nop(NELmax,NNref)

   INTEGER i,j,k,l,NEL
   NEL=0
   DO i=1,NELx
      DO j=1,NELy
         NEL=NEL+1
         DO k=1,3
            l=3*k-2
            nop(nel,l)=nny*(2*i+k-3)+2*j-1
            nop(nel,l+1)=nop(nel,l)+1
            nop(nel,l+2)=nop(nel,l)+2
         ENDDO
      ENDDO
   ENDDO
END SUBROUTINE GlobalNodeNumbering

 !>******************************************************************************
 !> TimeDiscretization creates the discrete points of time.
 !> -------------
 !> Params:
 !> -------------
 !> t:                   time coordinates
 !> M:                   number of time coordinates
 !>******************************************************************************
SUBROUTINE TimeDiscretization(t)
   use, intrinsic :: iso_fortran_env
   IMPLICIT NONE
   INCLUDE 'CHparameters.inc'

   !> Time coordinates
   REAL(real64) t(Tmax)

   REAL(real64) :: tfirst
   REAL(real64) :: tlast
   INTEGER j

   tfirst=0.
   tlast=2.
   deltat = (tlast-tfirst)/(Tmax-1)

   t(1)=tfirst
   DO j=2,Tmax
      t(j) = (j-1)*deltat
   ENDDO

END SUBROUTINE TimeDiscretization

!>******************************************************************************
 !> Initial&BoundaryConditions returns the initial solution (w,u) to be used by
 !> the Newton procedure.
 !> -------------
 !> Params:
 !> -------------
 !> x:                   Space coordinates
 !> t:                   Time coordinates
 !>
 !>******************************************************************************
SUBROUTINE InitialBoundaryConditions(x,y,nop, t, c, w)
   use, intrinsic :: iso_fortran_env
   IMPLICIT NONE
   INCLUDE 'CHparameters.inc'

   !> Space coordinates
   REAL(real64), INTENT(IN) :: x(NNmax),y(NNmax)
   !> Local to global space numbering
   INTEGER nop(NELmax,NNref)
   !> Time coordinates
   REAL(real64), INTENT(IN) :: t(Tmax)
   !> Order parameter - c - of the Cahn Hilliard Equation
   REAL(real64) c(NNmax,Tmax)
   !> Auxiliary variable - w - the chemical potential
   REAL(real64) w(NNmax,Tmax)

   LOGICAL:: isBoundaryPoint
   LOGICAL:: isInitialTime
   INTEGER i,timepoint

   !Impose Boundary and initial conditions on u, w
   DO i=1,NNmax
      isBoundaryPoint = ((x(i) == 0.) .or. (x(i) == 1.) .or. (y(i) == 0. .or. y(i) == 1.))
      DO timepoint=1,Tmax
         isInitialTime = (t(timepoint) == 0.)
         IF(isBoundaryPoint) THEN
            c(i,timepoint) = 0.
            w(i,timepoint) = 0.
         ELSEIF (isInitialTime) THEN
            c(i,timepoint) = RAND()
            w(i,timepoint) = 0.
         ELSE
            c(i,timepoint) = 0.
            w(i,timepoint) = 0.
         ENDIF
      ENDDO
   ENDDO

END SUBROUTINE InitialBoundaryConditions


 !>******************************************************************************
 !> WeakFormulation initializes and calculates the system matrices and the nonlinear term
 !> imposes the Boundary conditions
 !>
 !> -------------
 !> Params:
 !> -------------
 !> x:                   Space coordinates
 !> A,K:                 System matrices
 !> ngl:                 global coordinates in space
 !>******************************************************************************








 !>******************************************************************************
 !> Given an initial guess x for a root in n dimensions, take ntrial Newton-Raphson
 !> steps to improve the root. Stop if the root converges in either summed absolute
 !> variable increments tolx or summed absolute function values tolf.
 !> -------------
 !> Params:
 !> -------------
 !> Y:                   Newton solution vector

 !> w:                   Auxiliary vector
 !> A,K:                 System matrices
 !> f:                   Non-linear term of the weak formulation
 !> Df:                  Derivative of the non-linear term of the weak formulation
 !>******************************************************************************
SUBROUTINE NewtonProcedure(x,y,nop,t,timepoint,c,w,S,STIFF,MASS,P,DP)
   use, intrinsic :: iso_fortran_env
   USE BiCGStab_mod
   IMPLICIT NONE
   INCLUDE 'CHparameters.inc'

   !> Space coordinates
   REAL(real64) x(NNmax),y(NNmax)
   !> Local to global space numbering
   INTEGER nop(NEL,NNref)

   !> Time coordinates
   REAL(real64) t(Tmax)
   !> Current time point
   INTEGER timePoint

   !> Order parameter - c - of the Cahn Hilliard Equation
   REAL(real64) c(NNmax,Tmax)
   !> Auxiliary variable - w - the chemical potential
   REAL(real64) w(NNmax,Tmax)

   !> Newton solution vector
   REAL(real64) S(2*NNmax)

   !> System matrices
   !> STIFF - Stiffness
   !> MASS - Mass
   !> P - Non linear advective term from double well potential c^3-c
   REAL(real64) STIFF(NNmax,NNmax)
   REAL(real64) MASS(NNmax,NNmax)
   !> Non-linear term of the weak formulation
   REAL(real64) P(NNmax,Tmax)
   !> Derivative of the non-linear term of the weak formulation
   REAL(real64) DP(NNmax,NNmax,Tmax)



   REAL(real64),PARAMETER:: tolS=1.0E-8_real64, tolf=1.0E-8_real64
   REAL(real64) errorf,errorS,fjac(2*NNmax,2*NNmax),fvec(2*NNmax),deltaS(2*NNmax)
   INTEGER ntrial
   INTEGER trial

   INTEGER i, NEL
   REAL DomainBoundary(NNmax)




   IF (timePoint == 1) THEN
      CALL InitialBoundaryConditions(x,y,nop, t, c, w)
      return
   ENDIF

   !> Update S for the next cycle, update c,w
   S(1:NNmax) = w(:, timePoint)
   S(NNmax+1:2*NNmax) = c(:, timePoint)
   S = S * 0.99_real64


   !> Imposing the boundary conditions on the stationary solution vector S
   DomainBoundary=0.
   DO i=1,NNy
      DomainBoundary(i)=1.
   ENDDO
   DO i=NNy,NNmax,NNy
      DomainBoundary(i)=1.
   ENDDO
   DO i=NNmax,(NNx-1)*NNy+1,-1
      DomainBoundary(i)=1.
   ENDDO
   DO i=(NNx-1)*NNy+1,1,-NNy
      DomainBoundary(i)=1.
   ENDDO

   DO i=1,NNmax
      IF (DomainBoundary(i)==1.) THEN
         S(i) = 0.
      ENDIF
   ENDDO


   fvec=0.
   fjac=0.
   ntrial = 100

   ! Solve
   do trial=1,ntrial

      ! Form non linear term
      DO NEL=1,NELmax
         CALL NonlinearTerm(x,y,NEL,nop,timepoint,c,P,DP)
      ENDDO

      ! Stationary System
      CALL StationarySystem(timePoint,c,w,STIFF,MASS,P,DP,fjac,fvec)
      deltaS = BiCGStab(fjac, -fvec)
      S = S + deltaS
      w(:, timePoint) = S(1:NNmax)
      c(:, timePoint) = S(NNmax+1:2*NNmax)

      ! Update error terms
      errorf = 0.
      DO i=1,2*NNmax
         errorf = errorf + abs(fvec(i))
      ENDDO
      errorS = 0.
      DO i=1,2*NNmax
         errorS = errorS + abs(deltaS(i))
      ENDDO
      IF((errorS.le.tolS) .or. (errorf.le.tolf)) THEN
         w(:, MIN(timePoint+1, Tmax)) = S(1:NNmax)
         c(:, MIN(timePoint+1, Tmax)) = S(NNmax+1:2*NNmax)
         return
      ENDIF
   enddo

END SUBROUTINE


 !>******************************************************************************
 !> TestFunctions initializes the test functions for the given point.
 !> -------------
 !> Params:
 !> -------------
 !> point:               the point at which to calculate test functions
 !> ph:                  array of test functions
 !> phd:                 array of derivatives for test functions
 !>******************************************************************************

SUBROUTINE TestFunctions(x,y,phi,phic,phie)
   use, intrinsic :: iso_fortran_env
   IMPLICIT NONE
   INCLUDE 'CHparameters.inc'

   !> Point at which to form the test functions
   REAL(real64),DIMENSION(NNref):: phi,phic,phie
   REAL(real64), INTENT(IN) :: x,y
   REAL(real64) ::  l1x, l2x, l3x, l1y, l2y, l3y
   REAL(real64) :: dl1x,dl2x,dl3x, dl1y, dl2y, dl3y

   !*** One Statement Functions ***
   l1x=2.*x**2-3.*x+1.
   l2x=-4.*x**2+4.*x
   l3x=2.*x**2-x
   dl1x=4.*x-3.
   dl2x=-8.*x+4.
   dl3x=4.*x-1.

   l1y=2.*y**2-3.*y+1.
   l2y=-4.*y**2+4.*y
   l3y=2.*y**2-y
   dl1y=4.*y-3.
   dl2y=-8.*y+4.
   dl3y=4.*y-1.
   !*******************************

   phi(1)=l1x*l1y
   phi(2)=l1x*l2y
   phi(3)=l1x*l3y
   phi(4)=l2x*l1y
   phi(5)=l2x*l2y
   phi(6)=l2x*l3y
   phi(7)=l3x*l1y
   phi(8)=l3x*l2y
   phi(9)=l3x*l3y
   phic(1)=l1y*dl1x
   phic(2)=l2y*dl1x
   phic(3)=l3y*dl1x
   phic(4)=l1y*dl2x
   phic(5)=l2y*dl2x
   phic(6)=l3y*dl2x
   phic(7)=l1y*dl3x
   phic(8)=l2y*dl3x
   phic(9)=l3y*dl3x
   phie(1)=l1x*dl1y
   phie(2)=l1x*dl2y
   phie(3)=l1x*dl3y
   phie(4)=l2x*dl1y
   phie(5)=l2x*dl2y
   phie(6)=l2x*dl3y
   phie(7)=l3x*dl1y
   phie(8)=l3x*dl2y
   phie(9)=l3x*dl3y

END SUBROUTINE TestFunctions






