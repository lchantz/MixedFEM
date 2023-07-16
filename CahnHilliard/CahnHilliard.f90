!>******************************************************************************
!>    The form of the Cahn Hilliard equation is
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

   INTEGER i,timepoint

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

         WRITE(2,'(2x,E16.9,1x,E16.9,1x,E16.9,1x,E16.9)') t(timepoint), x(i),y(i), c(i,timepoint) 

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

   xorigin=0.
   xlast=1.
   deltax=(xlast-xorigin)/2*NELx

   yorigin=0.
   ylast=1.
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
 !> Initial&BoundaryConditions returns the initial solution (w_0,u_0) to be used by
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

   !Impose Periodic Boundary and initial conditions on c, w
   !can be redefined using the nop correspondence, so it appears on argument list
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


 !>***************************************************************************************
 !> WeakFormulation initializes and calculates the system matrices and the nonlinear term
 !> taking into account the Periodic Boundary conditions
 !>
 !> -------------
 !> Params:
 !> -------------
 !> x,y:                   Space coordinates
 !> STIFF,MASS,P:          System matrices
 !> ngl:                   ngl --> node global numbering
 !>**************************************************************************************

SUBROUTINE WeakFormulation(x,y,STIFF,MASS,P,nop)
   use, intrinsic :: iso_fortran_env
   IMPLICIT NONE
   INCLUDE 'CHparameters.inc'

   !> Space coordinates
   REAL(real64), INTENT(IN) :: x(NNmax),y(NNmax)

   !> System matrices
   !> STIFF - Stiffness
   !> MASS - Mass
   !> P - Non linear term from double well potential c^3-c
   REAL(real64) STIFF(NNmax,NNmax)
   REAL(real64) MASS(NNmax,NNmax)
   REAL(real64) P(NNmax,Tmax)

   !> Local to global space numbering / ngl --> node global numbering
   INTEGER nop(NELmax,NNref)

   INTEGER i,element
   LOGICAL:: isBoundaryPoint

   STIFF=0.
   MASS=0.
   P=0.

   DO element=1,NELmax
      CALL GaussianQuadrature(element,x,y,STIFF,MASS,nop)
   ENDDO

   !> Impose Periodic Boundary conditions on STIFF,MASS,P
   !> A,P set to zero to be eliminated for generating equations on c,w accuratelly
   !> Page 5 of Notes for explanation
   DO i=1,NNmax
      isBoundaryPoint = ((x(i) == 0.) .or. (x(i) == 1.) .or. (y(i) == 0. .or. y(i) == 1.))
      IF (isBoundaryPoint) THEN
         MASS(i,:) = 0.
         MASS(i,i) = 1.
         P(i,:) = 0.
         STIFF(i,:) = 0.
      ENDIF
   ENDDO      
  



END SUBROUTINE WeakFormulation

!>******************************************************************************
!> GaussianQuadrature numerical method for evaluation of integrals appearing in
!> system matrices.
!> -------------
!> Params:
!> -------------
!> element:             the current element
!> x,y:                 spatial coordinates
!> STIFF,MASS:          System matrices
!> ngl:                 Local to global space numbering
!>******************************************************************************

SUBROUTINE GaussianQuadrature(element,x,y,STIFF,MASS,nop)
   use, intrinsic :: iso_fortran_env
   IMPLICIT NONE
   INCLUDE 'CHparameters.inc'


   !> Current element
   INTEGER element

   !> Space coordinates !xpt,ypt
   REAL(real64), INTENT(IN) :: x(NNmax),y(NNmax)

   !> System matrices
   !> STIFF - Stiffness
   !> MASS - Mass
   REAL(real64) STIFF(NNmax,NNmax)
   REAL(real64) MASS(NNmax,NNmax)
  

   !> Local to global space numbering
   INTEGER nop(NELmax,NNref)
   INTEGER ngl(NNref)

   REAL(real64) phi(NNref), tphx(NNref), tphy(NNref),phic(NNref),phie(NNref)
   REAL(real64) gp(3), gw(3)
   REAL(real64) Xcomputational,Xc,Xe,Ycomputational,Yc,Ye,dett

   
   INTEGER i,j,k,l,m,n
   


   gw  =(/0.27777777777778, 0.444444444444, 0.27777777777778/)
   gp =(/0.1127016654    , 0.5           , 0.8872983346    /)

   DO i = 1,NNref
      ngl(i) = nop(element,i)
   ENDDO  

   phi=0.
   phic=0.
   phie=0.
   tphx=0.
   tphy=0.
   ! Loop over qauss points
   DO j = 1,3
      DO k = 1,3

      ! TODO: check test func | tphy,tphx,dett -> nan
      call TestFunctions(gp(j),gp(k),phi,phic,phie)

      ! Defines the computational domain coordinates and the 2-dimensional Jacobian dett
      Xcomputational=0.
      Xc=0.
      Xe=0.
      Ycomputational=0.
      Yc=0.
      Ye=0.
      DO n=1,NNref
         Xcomputational= Xcomputational + x(ngl(n)) * phi(n)
         Xc= Xc + x(ngl(n)) * phic(n)
         Xe= Xe + x(ngl(n)) * phie(n)
         Ycomputational= Ycomputational + y(ngl(n)) * phi(n)
         Yc= Yc + y(ngl(n)) * phic(n)
         Ye= Ye + y(ngl(n)) * phie(n)
      ENDDO
      dett=Xc*Ye-Xe*Yc

      DO i=1,NNref
         tphx(i)=(Ye*phic(i)-Yc*phie(i))/dett
         tphy(i)=(Xc*phie(i)-Xe*phic(i))/dett
      ENDDO

   



      DO l=1,NNref
         DO m=1,NNref
            STIFF(ngl(l),ngl(m)) = STIFF(ngl(l),ngl(m)) + gw(j)*gw(k)*dett*(tphx(l)*tphx(m)+tphy(l)*tphy(m))
            MASS(ngl(l),ngl(m)) = MASS(ngl(l),ngl(m)) + gw(j)*gw(k)*dett*(phi(l)*phi(m))
         ENDDO
      ENDDO
   ENDDO !End of loop over gauss points
   ENDDO !End of loop over gauss points

END SUBROUTINE GaussianQuadrature







 !>******************************************************************************
 !> Given an initial guess x for a root in n dimensions, take ntrial Newton-Raphson
 !> steps to improve the root. Stop if the root converges in either summed absolute
 !> variable increments tolx or summed absolute function values tolf.
 !> -------------
 !> Params:
 !> -------------
 !> S:                   Newton solution vector
 !> c:                   Order parameter   
 !> w:                   Auxiliary vector
 !> STIFF,MASS:                 System matrices
 !> P:                   Non-linear term of the weak formulation
 !> DP:                  Derivative of the non-linear term of the weak formulation
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


   !> Imposing the Periodic (essential)
   !> boundary conditions on the stationary solution vector S
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
   DO trial=1,ntrial

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
   ENDDO

END SUBROUTINE

!>******************************************************************************
!> NonlinearTerm evaluates the entries of the P(c,t) operator and the DP(c,t) matrix
!> at a fixed time step.
!> -------------
!> Params:
!> -------------
!> x,y:                 spatial coordinates
!> NEL:                 current element
!> nop:                 local to global node numbering
!> timePoint:           time point
!> c:                   order parameter
!> P:                   Non-linear term of the weak formulation
!> DP:                  Derivative of the non-linear term
!>******************************************************************************


SUBROUTINE NonlinearTerm(x,y,NEL,nop,timepoint,c,P,DP)
   use, intrinsic :: iso_fortran_env
   IMPLICIT NONE
   INCLUDE 'CHparameters.inc'

  !> Current element
   INTEGER NEL


   !> Space coordinates !xpt,ypt
   REAL(real64), INTENT(IN) :: x(NNmax),y(NNmax)

   !> timepoint
   INTEGER timepoint

   !> Order parameter - c - of the Cahn Hilliard Equation
   REAL(real64) c(NNmax,Tmax)

   !> non linear term
   REAL(real64) P(NNmax,Tmax)
   !> Derivative of the non-linear term
   REAL(real64) DP(NNmax,NNmax,Tmax)


   !> Local to global space numbering
   INTEGER nop(NELmax,NNref)
   INTEGER ngl(NNref)

   REAL(real64) phi(NNref), tphx(NNref), tphy(NNref),phic(NNref),phie(NNref)
   REAL(real64) gp(3), gw(3)
   REAL(real64) Xcomputational,Xc,Xe,Ycomputational,Yc,Ye,dett

   REAL(real64) ci,cx,cy
   INTEGER i,j,k,l,m,n
   


   gw  =(/0.27777777777778, 0.444444444444, 0.27777777777778/)
   gp =(/0.1127016654    , 0.5           , 0.8872983346    /)

   DO i = 1,NNref
      ngl(i) = nop(NEL,i)
   ENDDO

   phi=0.
   phic=0.
   phie=0.
   tphx=0.
   tphy=0.
   ! Loop over qauss points
   DO j = 1,3
      DO k = 1,3

      call TestFunctions(gp(j),gp(k),phi,phic,phie)

      ! Defines the computational domain coordinates and the 2-dimensional Jacobian dett
      Xcomputational=0.
      Xc=0.
      Xe=0.
      Ycomputational=0.
      Yc=0.
      Ye=0.
      DO n=1,NNref
         Xcomputational= Xcomputational + x(ngl(n)) * phi(n)
         Xc= Xc + x(ngl(n)) * phic(n)
         Xe= Xe + x(ngl(n)) * phie(n)
         Ycomputational= Ycomputational + y(ngl(n)) * phi(n)
         Yc= Yc + y(ngl(n)) * phic(n)
         Ye= Ye + y(ngl(n)) * phie(n)
      ENDDO
      dett=Xc*Ye-Xe*Yc

      DO i=1,NNref
         tphx(i)=(Ye*phic(i)-Yc*phie(i))/dett
         tphy(i)=(Xc*phie(i)-Xe*phic(i))/dett
      ENDDO


      ci = 0.
      cx = 0.
      cy = 0.
      DO i=1,NNref
         ci = ci + c(ngl(i),timePoint)*phi(i)
         cx = cx + c(ngl(i),timePoint)*tphx(i)
         cy = cy + c(ngl(i),timePoint)*tphy(i)
      ENDDO

      DO l=1,NNref
         DO m=1,NNref
            P(ngl(l),timepoint) = P(ngl(l),timepoint) + gw(j)*gw(k)*dett*(ci**3 - ci)*phi(l)
            DP(ngl(l),ngl(m),timepoint) =  DP(ngl(l),ngl(m),timepoint) &
               +gw(j)*gw(k)*dett*(3.*(ci**2) - 1.)*phi(l)*phi(m)
         ENDDO
      ENDDO

   ENDDO
ENDDO


END SUBROUTINE NonlinearTerm

!>******************************************************************************
!> Constructs the stationary system of each Newton procedure cycle
!> .
!> -------------
!> Params:
!> -------------
!> timepoint:           the current time point
!> c:                   Order parameter
!> w:                   Chemical Potential
!> STIFF,MASS:          System matrices
!> P:                   Non-linear term of the weak formulation
!> DP:                  Derivative of the non-linear term of the weak formulation
!> fjac:                jacobian matrix of newton procedure
!> fvec:                system of equations of newton procedure
!>*********************************************************************************
SUBROUTINE StationarySystem(timepoint,c,w,STIFF,MASS,P,DP,fjac,fvec)
   use, intrinsic :: iso_fortran_env
   IMPLICIT NONE
   INCLUDE 'CHparameters.inc'

   !> Current time step
   INTEGER timepoint
   !> Order parameter
   REAL(real64) c(NNmax,Tmax)
   !> Auxiliary variable-Chemical potential
   REAL(real64) w(NNmax,Tmax)
   !> System matrices
   !> STIFF - Stiffness
   !> MASS - Mass
   REAL(real64) STIFF(NNmax,NNmax)
   REAL(real64) MASS(NNmax,NNmax)
   !> Non-linear term of the weak formulation
   REAL(real64) P(NNmax,Tmax)
   !> Derivative of the non-linear term of the weak formulation
   REAL(real64) DP(NNmax,NNmax,Tmax)

   !> Jacobian matix of Newton system
   REAL(real64) fjac(2*NNmax,2*NNmax)
   !> Righ handside vector of Newton system
   REAL(real64) fvec(2*NNmax)

   INTEGER i,j


   DO i=1,NNmax
      DO j=1,NNmax
         fvec(i)= fvec(i) + MASS(i,j)*w(j,timepoint) - P(j,timepoint) - (epsilon)**2*STIFF(i,j)*c(j,timepoint)
         fjac(i,j)= fjac(i,j) + MASS(i,j)
      ENDDO

      DO j= NNmax+1,2*NNmax
         fjac(i,j)=  fjac(i,j) - DP(i,j-NNmax,timepoint) - (epsilon)**2*STIFF(i,j-NNmax)
      ENDDO
   ENDDO

   DO i=NNmax+1,2*NNmax
      DO j=1,NNmax
         fvec(i)= fvec(i) + deltat*STIFF(i-NNmax,j)*w(j,timepoint) + MASS(i-NNmax,j)*c(j,timePoint) - MASS(i-NNmax,j)*c(j,timePoint-1)
         fjac(i,j)= fjac(i,j) + deltat*STIFF(i-NNmax,j)
      ENDDO
      DO j= NNmax+1,2*NNmax
         fjac(i,j)= fjac(i,j) + MASS(i-NNmax,j-NNmax)
      ENDDO
   ENDDO

END SUBROUTINE



 !>******************************************************************************
 !> TestFunctions initializes the test functions for the given point.
 !> -------------
 !> Params:
 !> -------------
 !> x,y:               the point at which to calculate test functions
 !> ph:                  array of test functions
 !> phic:                array of derivatives on c for test functions
 !> phie:                array of derivatives on e for test functions
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






