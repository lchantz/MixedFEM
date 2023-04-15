!>******************************************************************************
!>    The form of the Burgers equation considered here is
!>
!>      du       du        d^2 u
!>      -- + u * -- = v * -----
!>      dt       dx        dx^2
!>
!>    for 0 < x < 1, and 0 < t < 2
!>
!>    Initial conditions:
!>     u(x,0) =  sin(pi*x).
!>
!>    Boundary conditions:
!>     u(0,t) = u(1,t) = 0.
!>
!>    viscosity parameter:
!>     v
!>******************************************************************************
program MixedFEM
   use, intrinsic :: iso_fortran_env
   USE BiCGStab_mod
   IMPLICIT NONE
   INCLUDE 'parameters.inc'

   !> Space coordinates
   REAL(real64) x(N)
   !> Time coordinates
   REAL(real64) t(M)

   !> Local to global space numbering
   INTEGER ngl(NEL,PPE)

   !> Solution of Burgers Equation
   REAL(real64) u(N,M)
   !> Auxiliary variable
   REAL(real64) w(N,M)
   !> System matrices
   REAL(real64) A(N,N), K(N,N)
   REAL(real64) f(M,N)
   !> Derivative of the non-linear term
   REAL(real64) Df(M,N,N)

   INTEGER i,timepoint

   OPEN(unit=1,file='BurgersSolution.dat')

   WRITE(*,'(25x,A)') 'Results u(x,t)'
   CALL SpaceDiscretization(x)
   CALL SetGlobalNodeNumbering(ngl)
   CALL TimeDiscretization(t)
   CALL WeakFormulation(x,A,K,f,ngl)

   Df=0.
   DO timepoint=1,M
      f(timepoint,:)=0.
      Df(timepoint,:,:)=0.
      CALL NewtonProcedure(timepoint,f,Df,u,w,A,K,ngl,x,t)
      DO i=1,N
         WRITE(1,*) t(timepoint), x(i), u(i,timepoint)
      ENDDO
   ENDDO
   CLOSE(1)
END program MixedFEM


!>******************************************************************************
!> SpaceDiscretization creates the discrete points in space.
!> -------------
!> Params:
!> -------------
!> x:       space coordinates
!> N:       number of space coordinates
!> NEL:     number of discrete elements
!>******************************************************************************
SUBROUTINE SpaceDiscretization(x)
   use, intrinsic :: iso_fortran_env
   IMPLICIT NONE
   INCLUDE 'parameters.inc'

   !> Space coordinates
   REAL(real64) x(N)

   REAL(real64) xfirst
   REAL(real64) xlast
   REAL(real64) deltax
   INTEGER i

   xfirst=0.
   xlast=1.
   deltax=(xlast-xfirst)/NEL
   x(1)=xfirst

   DO i=2,N
      x(i)=x(i-1)+deltax/(PPE-1)
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
SUBROUTINE SetGlobalNodeNumbering(ngl)
   use, intrinsic :: iso_fortran_env
   IMPLICIT NONE
   INCLUDE 'parameters.inc'

   !> Local to global space numbering
   INTEGER ngl(NEL,PPE)

   INTEGER i

   DO i=1,NEL
      ngl(i,1)=1+2*(i-1)
      ngl(i,2)=ngl(i,1)+1
      ngl(i,3)=ngl(i,2)+1
   ENDDO
END SUBROUTINE SetGlobalNodeNumbering

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
   INCLUDE 'parameters.inc'

   !> Time coordinates
   REAL(real64) t(M)

   REAL(real64) :: tfirst
   REAL(real64) :: tlast
   INTEGER j

   tfirst=0.
   tlast=2.
   deltat = (tlast-tfirst)/(M-1)
   t(1)=tfirst

   DO j=2,M
      t(j) = t(j-1) + deltat
   ENDDO

END SUBROUTINE TimeDiscretization

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
SUBROUTINE WeakFormulation(x,A,K,f,ngl)
   use, intrinsic :: iso_fortran_env
   IMPLICIT NONE
   INCLUDE 'parameters.inc'

   !> Space coordinates
   REAL(real64) x(N)
   !> System matrices
   REAL(real64) A(N,N), K(N,N)
   !> Non linear term
   REAL(real64) f(M,N)
   !> Local to global space numbering / ngl --> node global numbering
   INTEGER ngl(NEL,PPE)

   INTEGER element

   A=0.
   K=0.
   f=0.

   DO element=1,NEL
      CALL GaussianQuadrature(element,x,A,K,ngl)
   ENDDO

   !> Imposing Boundary conditions
   CALL BoundaryConditions(A)



END SUBROUTINE WeakFormulation

SUBROUTINE BoundaryConditions(A)

   use, intrinsic :: iso_fortran_env
   IMPLICIT NONE
   INCLUDE 'parameters.inc'


   REAL(real64) A(N,N)


   A(1,1)= A(1,1) + 1.
   A(N,N)= A(N,N) - 1.

END SUBROUTINE



!>******************************************************************************
!> GaussianQuadrature numerical method for evaluation of integrals appearing in
!> system matrices.
!> -------------
!> Params:
!> -------------
!> element:             the current element
!> x:                   space coordinates
!> u:                   Solution of Burgers Equation
!> A,K:                 System matrices
!> ngl:                 global coordinates in space
!>******************************************************************************

SUBROUTINE GaussianQuadrature(element,x,A,K,ngl)
   use, intrinsic :: iso_fortran_env
   IMPLICIT NONE
   INCLUDE 'parameters.inc'

   !> Space coordinates
   REAL(real64) x(N)
   !> Current element
   INTEGER element
   !> System matrices
   REAL(real64) A(N,N), K(N,N)

   !> Local to global space numbering
   INTEGER ngl(NEL,PPE)

   REAL(real64) ph(PPE), phd(PPE), phx(PPE)
   REAL(real64) gp(PPE), gw(PPE)
   REAL(real64) Xcomputational,Xc

   INTEGER gpindex
   INTEGER i,j
   INTEGER m1,n1


   gw  =(/0.27777777777778, 0.444444444444, 0.27777777777778/)
   gp =(/0.1127016654    , 0.5           , 0.8872983346    /)


   ! Loop over qauss points
   do gpindex=1,PPE

      call TestFunctions(gp(gpindex), ph, phd)

      ! Defines the computational domain coordinates and the 1D Jacobian Xc
      Xcomputational=0.
      Xc=0.
      do i=1,PPE
         Xcomputational= Xcomputational + x(ngl(element, i)) * ph(i)
         Xc= Xc + x(ngl(element, i)) * phd(i)
      enddo

      ! df/dx = (1/J) df/dc
      do i=1,PPE
         phx(i) = phd(i)/Xc
      enddo



      DO i=1,PPE
         m1 = ngl(element, i)
         DO j=1,PPE
            n1 = ngl(element, j)
            A(m1,n1) = A(m1,n1) + gw(gpindex) * Xc * ph(j) * phx(i)
            K(m1,n1) = K(m1,n1) + gw(gpindex) * Xc * ph(i) * ph(j)
         ENDDO
      ENDDO
   ENDDO !End of loop over gauss points

END SUBROUTINE GaussianQuadrature



!>******************************************************************************
!> TestFunctions initializes the test functions for the given point.
!> -------------
!> Params:
!> -------------
!> point:               the point at which to calculate test functions
!> ph:                  array of test functions
!> phd:                 array of derivatives for test functions
!>******************************************************************************

SUBROUTINE TestFunctions(r, ph, phd)
   use, intrinsic :: iso_fortran_env
   IMPLICIT NONE
   INCLUDE 'parameters.inc'

   !> Point at which to form the test functions
   REAL(real64):: r
   !> Test functions
   REAL(real64):: ph(PPE)
   !> Derivatives of test functions
   REAL(real64):: phd(PPE)

   ph(1) = 2. * (r**2) - 3. * r + 1.
   ph(2) = -4. * (r**2) + 4. * r
   ph(3) = 2. * (r**2) - r
   phd(1)=4. * r - 3.
   phd(2)=-8. * r + 4.
   phd(3)=4. * r - 1.

END SUBROUTINE TestFunctions

!>******************************************************************************
!> NonlinearTerm evaluates the entries of the f(u) operator and the Df(u) matrix
!> at a fixed time step.
!> -------------
!> Params:
!> -------------
!> f:                   Non-linear term of the weak formulation
!> timePoint:           time point
!> Df:                  Derivative of the non-linear term
!>******************************************************************************


SUBROUTINE NonlinearTerm(f,Df,element,ngl, timePoint,u,x)
   use, intrinsic :: iso_fortran_env
   IMPLICIT NONE
   INCLUDE 'parameters.inc'

   !> Non-linear term of the weak formulation
   REAL(real64) f(M,N)
   !> Derivative of the non-linear term
   REAL(real64) Df(M,N,N)
   !> Current time step to evaluate f
   INTEGER timePoint
   !> Solution of Burgers Equation
   REAL(real64) u(N,M)
   !> Space coordinates
   REAL(real64) x(N)

   !> Local to global space numbering
   INTEGER ngl(NEL,PPE)


   INTEGER i,j,element
   REAL(real64) ui

   REAL(real64) ph(PPE), phd(PPE), phx(PPE)
   REAL(real64) gp(PPE), gw(PPE)
   REAL(real64) Xcomputational,Xc



   INTEGER gpindex


   gw  =(/0.27777777777778, 0.444444444444, 0.27777777777778/)
   gp =(/0.1127016654    , 0.5           , 0.8872983346    /)



   ! Loop over qauss points
   do gpindex=1,PPE

      call TestFunctions(gp(gpindex), ph, phd)

      ! Defines the computational domain coordinates and the 1D Jacobian Xc
      Xcomputational=0.
      Xc=0.
      do i=1,PPE
         Xcomputational= Xcomputational + x(ngl(element, i)) * ph(i)
         Xc= Xc + x(ngl(element, i)) * phd(i)
      enddo

      ! df/dx = (1/J) df/dc
      do i=1,PPE
         phx(i) = phd(i)/Xc
      enddo

      ui = 0.
      DO i=1,PPE
         ui = ui + u(ngl(element,i),timePoint)*ph(i)
      ENDDO

      DO i=1,PPE
         f(timePoint,ngl(element,i)) = f(timePoint,ngl(element,i)) +(0.5_real64)* gw(gpindex) * Xc *((ui)**2) * ph(i)
         DO j=1,PPE
            Df(timePoint,ngl(element,i),ngl(element,j)) =  Df(timePoint,ngl(element,i),ngl(element,j)) &
               + u(ngl(element,j),timePoint) * gw(gpindex) * Xc * (ui)* ph(j)*ph(i)
         ENDDO
      ENDDO

   ENDDO


END SUBROUTINE NonlinearTerm



!>******************************************************************************
!> Given an initial guess x for a root in n dimensions, take ntrial Newton-Raphson
!> steps to improve the root. Stop if the root converges in either summed absolute
!> variable increments tolx or summed absolute function values tolf.
!> -------------
!> Params:
!> -------------
!> Y:                   Newton solution vector
!> u:                   Solution of Burgers Equation
!> w:                   Auxiliary vector
!> A,K:                 System matrices
!> f:                   Non-linear term of the weak formulation
!> Df:                  Derivative of the non-linear term of the weak formulation
!>******************************************************************************
SUBROUTINE NewtonProcedure(timePoint,f,Df,u,w,A,K,ngl,x,t)
   use, intrinsic :: iso_fortran_env
   USE BiCGStab_mod
   IMPLICIT NONE
   INCLUDE 'parameters.inc'

   !> Solution of Burgers Equation
   REAL(real64) u(N,M)
   !> Auxiliary variable
   REAL(real64) w(N,M)
   !> Newton solution vector
   REAL(real64) Y(2*N)
   !> Current time step
   INTEGER timePoint
   !> System matrices
   REAL(real64) A(N,N), K(N,N)
   !> Non-linear term of the weak formulation
   REAL(real64) f(M,N)
   !> Derivative of the non-linear term of the weak formulation
   REAL(real64) Df(M,N,N)
   !> Local to global space numbering
   INTEGER ngl(NEL,PPE)
   !> Space coordinates
   REAL(real64) x(N)
   !> Time coordinates
   REAL(real64) t(M)

   REAL(real64),PARAMETER:: toly=1.0E-8_real64, tolf=1.0E-8_real64
   REAL(real64) errorf,errory,fjac(2*N,2*N),fvec(2*N),deltaY(2*N)
   INTEGER ntrial
   INTEGER trial

   INTEGER i, element



   IF (timePoint == 1) THEN
      CALL InitialSolution(x, t, u, w)
      return
   ENDIF

   ! Update Y, u, w with perturbation
   Y(1:N) = w(:, timePoint)
   Y(N+1:2*N) = u(:, timePoint)
   Y = Y * 0.99_real64
   Y(1) = 0.0_real64
   Y(N) = 0.0_real64
   Y(N+1) = 0.0_real64
   Y(2*N) = 0.0_real64

   fvec=0.
   fjac=0.
   ntrial = 100

   ! Solve
   do trial=1,ntrial
      ! Form non linear term
      DO element=1,NEL
         CALL NonlinearTerm(f,Df,element,ngl,timePoint,u,x)
      ENDDO
      ! Calculate system for timepoint
      CALL FormulateSystem(timePoint,f,Df,fvec,fjac,u,w,A,K)
      deltaY = BiCGStab(fjac, -fvec)
      Y = Y + deltaY
      w(:, timePoint) = Y(1:N)
      u(:, timePoint) = Y(N+1:2*N)
      ! Update error terms
      errorf = 0.
      DO i=1,2*N
         errorf = errorf + abs(fvec(i))
      ENDDO
      errory = 0.
      DO i=1,2*N
         errory = errory + abs(deltaY(i))
      ENDDO
      IF((errory.le.toly) .or. (errorf.le.tolf)) THEN
         w(:, MIN(timePoint+1, M)) = Y(1:N)
         u(:, MIN(timePoint+1, M)) = Y(N+1:2*N)
         return
      ENDIF
   enddo

   return

END SUBROUTINE

!>******************************************************************************
!> InitialSolution returns the initial solution (w,u) to be used by
!> the Newton procedure.
!> -------------
!> Params:
!> -------------
!> x:                   Space coordinates
!> t:                   Time coordinates
!> u:                   Solution of Burgers Equation
!> w:                   Auxiliary vector
!>******************************************************************************
SUBROUTINE InitialSolution(x, t, u, w)
   use, intrinsic :: iso_fortran_env
   IMPLICIT NONE
   INCLUDE 'parameters.inc'

   !> Space coordinates
   REAL(real64), INTENT(IN) :: x(N)
   !> Time coordinates
   REAL(real64), INTENT(IN) :: t(M)
   !> Solution of Burgers Equation
   REAL(real64), INTENT(OUT) :: u(N,M)
   !> Auxiliary variable
   REAL(real64), INTENT(OUT) :: w(N,M)

   LOGICAL:: isBoundaryPoint
   LOGICAL:: isInitialTime
   INTEGER i,j

   !Impose Boundary and initial conditions on u, w
   DO i=1,N
      isBoundaryPoint = ((x(i) == 0.) .or. (x(i) == 1.))
      DO j=1,M
         isInitialTime = (t(j) == 0.)
         IF(isBoundaryPoint) THEN
            u(i,j) = 0.
            w(i,j) = 0.
         ELSEIF (isInitialTime) THEN
            u(i,j) = sin(pi*x(i))
            w(i,j) = 0.
         ELSE
            u(i,j) = 0.
            w(i,j) = 0.
         ENDIF
      ENDDO
   ENDDO

END SUBROUTINE InitialSolution



!>******************************************************************************
!> FormulateSystem constructs the stationary system of each Newton procedure cycle
!> .
!> -------------
!> Params:
!> -------------
!> timepoint:           the current time point
!> u:                   Solution of Burgers Equation
!> w:                   Auxiliary vector
!> A,K:                 System matrices
!> f:                   Non-linear term of the weak formulation
!> Df:                  Derivative of the non-linear term of the weak formulation
!>******************************************************************************
SUBROUTINE FormulateSystem(timePoint,f,Df,fvec,fjac,u,w,A,K)
   use, intrinsic :: iso_fortran_env
   IMPLICIT NONE
   INCLUDE 'parameters.inc'

   !> Current time step
   INTEGER timePoint
   !> Jacobian matix of Newton system
   REAL(real64) fjac(2*N,2*N)
   !> Righ handside vector of Newton system
   REAL(real64) fvec(2*N)
   !> Solution of Burgers Equation
   REAL(real64) u(N,M)
   !> Auxiliary variable
   REAL(real64) w(N,M)
   !> System matrices
   REAL(real64) A(N,N),K(N,N)
   !> Non-linear term of the weak formulation
   REAL(real64) f(M,N)
   !> Derivative of the non-linear term
   REAL(real64) Df(M,N,N)

   INTEGER i,j


   DO i=1,N
      DO j=1,N
         fvec(i)= fvec(i) + K(i,j) * w(j,timePoint) + f(timePoint,i)  + v * A(i,j) * u(j,timePoint)
         fjac(i,j)= fjac(i,j) +  K(i,j)
      ENDDO

      DO j= N+1,2*N
         fjac(i,j)=  fjac(i,j) + Df(timePoint,i,j-N) + v * A(i,j-N)
      ENDDO
   ENDDO

   DO i=N+1,2*N
      DO j=1,N
         fvec(i)= fvec(i) + K(i-N,j) * u(j,timePoint) + deltat * A(i-N,j) *  w(j,timePoint) - K(i-N,j) *  u(j,timePoint-1)
         fjac(i,j)= fjac(i,j) + deltat * A(i-N,j)
      ENDDO
      DO j= N+1,2*N
         fjac(i,j)= fjac(i,j) + K(i-N,j-N)
      ENDDO
   ENDDO

END SUBROUTINE


