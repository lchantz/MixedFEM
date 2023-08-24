!>******************************************************************************
!>    The form of the Cahn Hilliard equation is...
!>******************************************************************************
program CahnHilliard
   use, intrinsic :: iso_fortran_env
   USE BiCGStab_mod
   IMPLICIT NONE
   INCLUDE 'CHParameters.inc'

   !> Space coordinates
   REAL(real64) x(NNmax),y(NNmax)

   !> Time coordinates
   REAL(real64) t(Tmax)

   !> Local to global node numbering
   INTEGER nop(NELmax,NNref)

   !> Location indices to check for Boundary Elements
   INTEGER BottomElement(NELmax)
   INTEGER RightElement(NELmax)
   INTEGER TopElement(NELmax)
   INTEGER LeftElement(NELmax)

   !> Boundary element control index
   INTEGER BoundaryElement(NELmax)

   !> Order parameter - c - of the Cahn Hilliard Equation
   REAL(real64) c(NNmax,Tmax)


   !> Auxiliary variable - w - the chemical potential
   REAL(real64) w(NNmax,Tmax)
  



   !> System matrices
   !> STIFF - Stiffness
   !> MASS - Mass
   !> P - Non linear advective term from double well potential c^3-c


   REAL(real64) STIFF(NNmax,NNmax)


   REAL(real64) MASS(NNmax,NNmax)
 

   REAL(real64) P(NNmax,Tmax)
 


   !> Derivative of the non-linear term from double well potential
   REAL(real64) DP(NNmax,NNmax,Tmax)
 

   REAL(real64) FLUXBottom(NNmax,NNmax)
   REAL(real64) FLUXRight(NNmax,NNmax)
   REAL(real64) FLUXTop(NNmax,NNmax)
   REAL(real64) FLUXLeft(NNmax,NNmax)

   !> Derivative of the Surface FLUX integral term
   REAL(real64) FLUX(NNmax,NNmax)


   !> Newton tolerance parameters
   REAL(real64),PARAMETER:: tolS=1.0E-6_real64, tolf=1.0E-4_real64
   INTEGER, PARAMETER :: ntrial = 150,m=100
  
   !> Newton solution vector
   REAL(real64) fvec(2*NNmax)
   REAL(real64) fjac(2*NNmax,2*NNmax)


   ! REAL(real64) S(2*NNmax,ntrial)
   ! REAL(real64) dS(2*NNmax,ntrial)
   REAL(real64) dSol(2*NNmax)
   REAL(real64) Sol(2*NNmax)


   !> Norms
   ! REAL(real64) dS_norm
   ! REAL(real64) c_change_norm
   ! REAL(real64) w_change_norm

   INTEGER element,i,timepoint,trial
   LOGICAL CONVERGENT

   OPEN(unit=1,file='InitialCondition.dat')
   OPEN(unit=2,file='OrderParameter.dat')

   WRITE(*,'(25x,A)') 'Results c(n,t)'
   CALL SpaceDiscretization(x,y)
   CALL GlobalNodeNumbering(nop)

   !> Mark boundary elements
   CALL BoundaryElements(BottomElement,RightElement,TopElement,LeftElement,BoundaryElement)

   !> Set time grid
   t=0.0_real64
   CALL TimeDiscretization(t)

   !>****************************************************************************
   !> Calculation of Stiffness,Mass,FLUX for the full dimension (Neumann problem)
   !>****************************************************************************
   !> Initialization of system matrices and vectors
   STIFF=0.0_real64
   MASS=0.0_real64
   DO element=1,NELmax
      CALL StiffnessMassMatrices(element,x,y,nop,STIFF,MASS)
      !> Initialization of the derivative of the FLUX matrix. Independent from c,w
      ! FLUXBottom=0.
      ! FLUXRight=0.
      ! FLUXTop=0.
      ! FLUXLeft=0.
      ! FLUX=0.
      ! FLUXper=0.
      ! CALL BoundaryFluxMatrix(element,x,y,nop,BottomElement,RightElement,TopElement,LeftElement,BoundaryElement,&
      ! FLUXBottom,FLUXRight,FLUXTop,FLUXLeft,FLUX)
   ENDDO

   !>********************************************************
   !> Calculation of the reduced matrices due to periodicity
   !>********************************************************
   !CALL PeriodicBoundaryConditions(STIFF,STIFFper,MASS,MASSper,&
   !FLUX,FLUXper)


   !>****************************************************************
   !> Set the initial conditions for the c,w variables of the problem
   !> For timepoint = 1
   !>*****************************************************************
   ! timepoint = 1
   c=0.0_real64
   w=0.0_real64
   CALL InitialConditions(x,y,c,w)
   DO i=1,NNmax
      WRITE(1,'(2x,E16.9,1x,E16.9,1x,E16.9,1x,E16.9)') x(i),y(i),c(i,1)
   ENDDO
   CLOSE(1)
   
   !>*****************************************************************
   !> Initialization of the stationary system vectors, matrices
   !>*****************************************************************

   
   !> Time loop to attain the solution at every timepoint
   DO timepoint=2,Tmax
      print *, "Time point: ", timepoint
      
      trial=0
      CONVERGENT=.false. 
      Sol=0.0_real64

      !> Choose the initial guess for the first newton iteration
      Sol(1:NNmax) = w(:,timepoint-1) !+ (w(:,timepoint) - w(:,timepoint-1)) * (deltat)
      Sol(NNmax+1:2*NNmax) = c(:,timepoint-1) !+ (c(:,timepoint) - c(:,timepoint-1)) * (deltat)
      
      
      ! !> Start the Newton iteration loop
      ! NewtonTrials: DO trial=1,ntrial
      DO WHILE (.NOT.CONVERGENT .AND. trial < ntrial)
         dSol=0.0_real64
      
            !>*****************************************************************
            !> Calculate the double well potential and the derivative
            !> Depend on c, the order parameter
            !>*****************************************************************
         P=0.0_real64
         DP=0.0_real64
            DO element=1,NELmax
               CALL PotentialTerm(x,y,element,nop,timepoint,c,P,DP)
            ENDDO

            !>******************************************************************
            !> Set the Stationary System to be solved at every Newton iteration
            !> ___________________Fjac.dx = -Fvec_______________________________
            !>******************************************************************
            fvec=0.0_real64
            fjac=0.0_real64
            CALL StationarySystem(timePoint,c,w,STIFF,MASS,P,DP,FLUX,&
               fjac,fvec)

            !> Solve the system - find the correction dSol
               fvec = -fvec
               dSol = BiCGStab(fjac,fvec)
               ! CALL solve_system(NNmax,fjac,fvec,dSol)
               ! CALL gausselim(2*NNmax,fjac,fvec,dSol)
              
               !> Update the solution 
                Sol(:) = Sol(:) + dSol(:)
                 

              !> Check for convergence
              CONVERGENT = (NORM2(dSol) < tolS) 
              print*,"trial", trial
              print*,"CONVERGENCE", CONVERGENT
              
              trial = trial + 1
              w(:,timepoint) = Sol(1:NNmax)
              c(:,timepoint) = Sol(NNmax+1:2*NNmax)

            ENDDO  !trials

         IF (CONVERGENT) THEN
            write(*,*) 'Newton-Raphson converged in', trial, 'iterations'
         ELSE
            write(*,*) 'Newton-Raphson did not converge within', ntrial, 'iterations'
         ENDIF

   ENDDO !> Time loop

   !> Write the results
   DO i=1,NNmax
      WRITE(2,'(2x,E16.9,1x,E16.9,1x,E16.9,1x,E16.9)') x(i),y(i),c(i,7)
   ENDDO

   CLOSE(2)
 
END program CahnHilliard


!>******************************************************************************
!> SpaceDiscretization creates the discrete points in space.
!> -------------
!> Params:
!> -------------
!> x:       space coordinates
!> N:       number of space coordinates
!> NEL:     number of discrete elements
!>******************************************************************************
SUBROUTINE SpaceDiscretization(x,y)
   use, intrinsic :: iso_fortran_env
   IMPLICIT NONE
   INCLUDE 'CHParameters.inc'

   !> Space coordinates
   REAL(real64) x(NNmax),y(NNmax)

   REAL(real64) xorigin,yorigin,xlast,ylast
   REAL(real64) deltax,deltay
   INTEGER i,j,NNode

   x=0.0_real64
   y=0.0_real64

   xorigin=0.0_real64
   xlast=1.0_real64
   deltax=(xlast-xorigin)/real(2*NELx)

   yorigin=0.0_real64
   ylast=1.0_real64
   deltay=(ylast-yorigin)/real(2*NELy)

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
 !> Creates the mapping between the local coordinates
 !> of a point on the element to the global coordinates in space.
 !> -------------
 !> Params:
 !> -------------
 !> nop:   global coordinates in space
 !>******************************************************************************
SUBROUTINE GlobalNodeNumbering(nop)
   use, intrinsic :: iso_fortran_env
   IMPLICIT NONE
   INCLUDE 'CHParameters.inc'

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

SUBROUTINE BoundaryElements(BottomElement,RightElement,TopElement,LeftElement,BoundaryElement)
   use, intrinsic :: iso_fortran_env
   IMPLICIT NONE
   INCLUDE 'CHParameters.inc'

   !> Location indices to check for Boundary Elements
   INTEGER BottomElement(NELmax)
   INTEGER RightElement(NELmax)
   INTEGER TopElement(NELmax)
   INTEGER LeftElement(NELmax)

   !> Boundary element control index
   INTEGER BoundaryElement(NELmax)
   INTEGER i

   !> Initialization of indices
   BottomElement=0
   RightElement=0
   TopElement=0
   LeftElement=0
   BottomElement=0



   !> Check for boundary elements
   DO i=1,(NELx-1)*NELy+1,NELy !bottom side elements
      BottomElement(i)=1
      BoundaryElement(i)=1
   ENDDO
   DO i=(NELx-1)*NELy+1,NELmax,1 !right side elements
      RightElement(i)=1
      BoundaryElement(i)=1
   ENDDO
   DO i=NELmax,NELy,-NELy !top side elements
      TopElement(i)=1
      BoundaryElement(i)=1
   ENDDO
   DO i=NELy,1,-1 !left side elements
      LeftElement(i)=1
      BoundaryElement(i)=1
   ENDDO

END SUBROUTINE



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
   INCLUDE 'CHParameters.inc'

   !> Time coordinates
   REAL(real64) t(Tmax)

   REAL(real64) :: tfirst
   ! REAL(real64) :: tlast
   INTEGER j

   tfirst=0.0_real64
   ! tlast=2.0_real64
   ! deltat = (tlast-tfirst)/real(Tmax-1)

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
SUBROUTINE InitialConditions(x,y,c,w)
   use, intrinsic :: iso_fortran_env
   IMPLICIT NONE
   INCLUDE 'CHParameters.inc'

   !> Space coordinates
   REAL(real64) x(NNmax),y(NNmax)

   !> Order parameter - c - of the Cahn Hilliard Equation
   REAL(real64) c(NNmax,Tmax)
   

   !> Auxiliary variable - w - the chemical potential
   REAL(real64) w(NNmax,Tmax)
   

   INTEGER i

   !> Set to a random distribution for initial time
   DO i=1,NNmax
       call random_number(c(i,1))
         w(i,1) = 0.0_real64
   ENDDO
   


END SUBROUTINE InitialConditions




!>******************************************************************************
!> StiffnessMassMatrices numerical method for evaluation of integrals appearing in
!> system matrices.
!> -------------
!> Params:
!> -------------
!> element:             the current element
!> x,y:                 spatial coordinates
!> STIFF,MASS:          System matrices
!> ngl:                 Local to global space numbering
!>******************************************************************************

SUBROUTINE StiffnessMassMatrices(element,x,y,nop,STIFF,MASS)
   use, intrinsic :: iso_fortran_env
   IMPLICIT NONE
   INCLUDE 'CHParameters.inc'


   !> Current element
   INTEGER element

   !> Space coordinates !xpt,ypt
   REAL(real64), INTENT(IN) :: x(NNmax),y(NNmax)

   !> Local to global space numbering
   INTEGER nop(NELmax,NNref)
   INTEGER ngl(NNref)

   !> System matrices
   !> STIFF - Stiffness
   !> MASS - Mass

   REAL(real64) STIFF(NNmax,NNmax)
   REAL(real64) MASS(NNmax,NNmax)




   REAL(real64) phi(NNref), tphx(NNref), tphy(NNref),phic(NNref),phie(NNref)
   REAL(real64) gp(3), gw(3)
   REAL(real64) Xcomputational,Xc,Xe,Ycomputational,Yc,Ye,dett


   INTEGER i,j,k,l,m,n



   gw  =(/0.27777777777778, 0.444444444444, 0.27777777777778/)
   gp =(/0.1127016654    , 0.5           , 0.8872983346    /)

   DO i = 1,NNref
      ngl(i) = nop(element,i)
   ENDDO

   phi=0.0_real64
   phic=0.0_real64
   phie=0.0_real64
   tphx=0.0_real64
   tphy=0.0_real64
   ! Loop over qauss points
   DO j = 1,3
      DO k = 1,3

         ! TODO: check test func | tphy,tphx,dett -> nan
         call TestFunctions(gp(j),gp(k),phi,phic,phie)

         ! Defines the computational domain coordinates and the 2-dimensional Jacobian dett
         Xcomputational=0.0_real64
         Xc=0.0_real64
         Xe=0.0_real64
         Ycomputational=0.0_real64
         Yc=0.0_real64
         Ye=0.0_real64
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

END SUBROUTINE StiffnessMassMatrices

SUBROUTINE BoundaryFluxMatrix(element,x,y,nop,BottomElement,RightElement,TopElement,LeftElement,BoundaryElement,&
   FLUXBottom,FLUXRight,FLUXTop,FLUXLeft,FLUX)
   use, intrinsic :: iso_fortran_env
   IMPLICIT NONE
   INCLUDE 'CHParameters.inc'


   !> Current element
   INTEGER element

   !> Space coordinates !xpt,ypt
   REAL(real64), INTENT(IN) :: x(NNmax),y(NNmax)

   !> Local to global space numbering
   INTEGER nop(NELmax,NNref)
   INTEGER ngl(NNref)

   !> FLUX contributions
   REAL(real64) FLUXBottom(NNmax,NNmax)
   REAL(real64) FLUXRight(NNmax,NNmax)
   REAL(real64) FLUXTop(NNmax,NNmax)
   REAL(real64) FLUXLeft(NNmax,NNmax)

   !> Derivative of the Surface FLUX integral term
   REAL(real64) FLUX(NNmax,NNmax)

   !> Location indices to check for Boundary Elements
   INTEGER BottomElement(NELmax)
   INTEGER RightElement(NELmax)
   INTEGER TopElement(NELmax)
   INTEGER LeftElement(NELmax)

   !> Boundary element control index
   INTEGER BoundaryElement(NELmax)

   REAL(real64) phi(NNref), tphx(NNref), tphy(NNref),phic(NNref),phie(NNref)
   REAL(real64) gp(3), gw(3)
   REAL(real64) Xcomputational,Xc,Xe,Ycomputational,Yc,Ye,dett


   INTEGER i,j,k,l,m,n



   gw  =(/0.27777777777778, 0.444444444444, 0.27777777777778/)
   gp =(/0.1127016654    , 0.5           , 0.8872983346    /)

   DO i = 1,NNref
      ngl(i) = nop(element,i)
   ENDDO

   phi=0.0_real64
   phic=0.0_real64
   phie=0.0_real64
   tphx=0.0_real64
   tphy=0.0_real64
   !> calculate the contribution of the surface integral term
   IF (BoundaryElement(element)==1) THEN
      !> check which side of the boundary the element belongs to
      IF (BottomElement(element)==1)THEN
         DO j=1,3
            CALL TestFunctions(gp(j),0._real64,phi, phic, phie)

            Xcomputational=0.0_real64
            Xc=0.0_real64
            Xe=0.0_real64
            Ycomputational=0.0_real64
            Yc=0.0_real64
            Ye=0.0_real64
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


            ! ci = 0.
            ! cx = 0.
            ! cy = 0.
            ! DO i=1,NNref
            !    ci = ci + c(ngl(i),timePoint)*phi(i)
            !    cx = cx + c(ngl(i),timePoint)*tphx(i)
            !    cy = cy + c(ngl(i),timePoint)*tphy(i)
            ! ENDDO

            ! wi = 0.
            ! wx = 0.
            ! wy = 0.
            ! DO i=1,NNref
            !    wi = wi + c(ngl(i),timePoint)*phi(i)
            !    wx = wx + c(ngl(i),timePoint)*tphx(i)
            !    wy = wy + c(ngl(i),timePoint)*tphy(i)
            ! ENDDO

            DO m=1,7,3
               DO n=1,NNref
                  FLUXBottom(ngl(m),ngl(n)) = FLUXBottom(ngl(m),ngl(n)) &
                     - gw(j)*Xc*phi(m)*tphy(n)
               ENDDO
            ENDDO
         ENDDO
      ELSEIF (RightElement(element)==1)THEN
         DO j=1,3
            CALL TestFunctions(1._real64,gp(j),phi, phic, phie)

            Xcomputational=0.0_real64
            Xc=0.0_real64
            Xe=0.0_real64
            Ycomputational=0.0_real64
            Yc=0.0_real64
            Ye=0.0_real64
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


            ! ci = 0.
            ! cx = 0.
            ! cy = 0.
            ! DO i=1,NNref
            !    ci = ci + c(ngl(i),timePoint)*phi(i)
            !    cx = cx + c(ngl(i),timePoint)*tphx(i)
            !    cy = cy + c(ngl(i),timePoint)*tphy(i)
            ! ENDDO

            ! wi = 0.
            ! wx = 0.
            ! wy = 0.
            ! DO i=1,NNref
            !    wi = wi + c(ngl(i),timePoint)*phi(i)
            !    wx = wx + c(ngl(i),timePoint)*tphx(i)
            !    wy = wy + c(ngl(i),timePoint)*tphy(i)
            ! ENDDO

            DO m=7,9
               DO n=1,NNref
                  FLUXRight(ngl(m),ngl(n)) = FLUXRight(ngl(m),ngl(n)) &
                     + gw(j)*Ye*phi(m)*tphx(n)
               ENDDO
            ENDDO
         ENDDO
      ELSEIF (TopElement(element)==1)THEN
         DO j=1,3
            CALL TestFunctions(gp(j),1._real64,phi, phic, phie)

            Xcomputational=0.0_real64
            Xc=0.0_real64
            Xe=0.0_real64
            Ycomputational=0.0_real64
            Yc=0.0_real64
            Ye=0.0_real64
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


            ! ci = 0.
            ! cx = 0.
            ! cy = 0.
            ! DO i=1,NNref
            !    ci = ci + c(ngl(i),timePoint)*phi(i)
            !    cx = cx + c(ngl(i),timePoint)*tphx(i)
            !    cy = cy + c(ngl(i),timePoint)*tphy(i)
            ! ENDDO

            ! wi = 0.
            ! wx = 0.
            ! wy = 0.
            ! DO i=1,NNref
            !    wi = wi + c(ngl(i),timePoint)*phi(i)
            !    wx = wx + c(ngl(i),timePoint)*tphx(i)
            !    wy = wy + c(ngl(i),timePoint)*tphy(i)
            ! ENDDO

            DO m=3,9,3
               DO n=1,NNref
                  FLUXTop(ngl(m),ngl(n)) = FLUXTop(ngl(m),ngl(n)) &
                     + gw(j)*Xc*phi(m)*tphy(n)
               ENDDO

            ENDDO
         ENDDO
      ELSEIF (LeftElement(element)==1)THEN
         DO j=1,3
            CALL TestFunctions(0._real64,gp(j),phi, phic, phie)

            Xcomputational=0.0_real64
            Xc=0.0_real64
            Xe=0.0_real64
            Ycomputational=0.0_real64
            Yc=0.0_real64
            Ye=0.0_real64
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


            ! ci = 0.
            ! cx = 0.
            ! cy = 0.
            ! DO i=1,NNref
            !    ci = ci + c(ngl(i),timePoint)*phi(i)
            !    cx = cx + c(ngl(i),timePoint)*tphx(i)
            !    cy = cy + c(ngl(i),timePoint)*tphy(i)
            ! ENDDO

            ! wi = 0.
            ! wx = 0.
            ! wy = 0.
            ! DO i=1,NNref
            !    wi = wi + c(ngl(i),timePoint)*phi(i)
            !    wx = wx + c(ngl(i),timePoint)*tphx(i)
            !    wy = wy + c(ngl(i),timePoint)*tphy(i)
            ! ENDDO

            DO m=1,3
               DO n=1,NNref
                  FLUXLeft(ngl(m),ngl(n)) = FLUXLeft(ngl(m),ngl(n)) &
                     - gw(j)*Ye*phi(m)*tphx(n)
               ENDDO
            ENDDO
         ENDDO


      ENDIF
   ENDIF

   !> Summation of the components to the total FLUX
   DO k=1,NNref
      DO l=1,NNref
         FLUX(ngl(k),ngl(l)) = FLUXBottom(ngl(k),ngl(l)) + FLUXRight(ngl(k),ngl(l)) &
            + FLUXTop(ngl(k),ngl(l)) + FLUXLeft(ngl(k),ngl(l))
      ENDDO

   ENDDO
END SUBROUTINE

SUBROUTINE PeriodicBoundaryConditions(STIFF,STIFFper,MASS,MASSper,&
   FLUX,FLUXper)
   use, intrinsic :: iso_fortran_env
   IMPLICIT NONE
   INCLUDE 'CHParameters.inc'

   !> STIFF - Stiffness
   !> MASS - Mass
   REAL(real64) STIFF(NNmax,NNmax)
   REAL(real64) STIFFper(Nhp,Nhp)

   REAL(real64) MASS(NNmax,NNmax)
   REAL(real64) MASSper(Nhp,Nhp)


   REAL(real64) FLUX(NNmax,NNmax)
   REAL(real64) FLUXper(Nhp,Nhp)


   INTEGER i,j

   !> innner nodes block
   DO i=1,NNy+1,Nhp
      DO j=1,NNy+1,Nhp
         MASSper(i,j) = MASS(i,j)
         STIFFper(i,j) = STIFF(i,j)
         FLUXper(i,j) = FLUX(i,j)
      ENDDO
   ENDDO

   !> 1st block
   DO i=1,NNy
      MASSper(i,i) = MASS(i,i) + MASS(i,Nhp+i) +MASS(Nhp+i,i) + MASS(Nhp+i,Nhp+i)
      STIFFper(i,i) = STIFF(i,i) + STIFF(i,Nhp+i) +STIFF(Nhp+i,i) + STIFF(Nhp+i,Nhp+i)
      FLUXper(i,i) = FLUX(i,i) + FLUX(i,Nhp+i) +FLUX(Nhp+i,i) + FLUX(Nhp+i,Nhp+i)
      DO j=1,NNy
         IF (i>j) THEN
            MASSper(i,j) = MASS(i,j) + MASS(i,Nhp+j)
            STIFFper(i,j) = STIFF(i,j) + STIFF(i,Nhp+j)
            FLUXper(i,j) = FLUX(i,j) + FLUX(i,Nhp+j)
         ELSEIF (i<j) THEN
            MASSper(i,j) = MASS(i,j) + MASS(Nhp+i,j)
            STIFFper(i,j) = STIFF(i,j) + STIFF(Nhp+i,j)
            FLUXper(i,j) = FLUX(i,j) + FLUX(Nhp+i,j)
         ENDIF
      ENDDO
   ENDDO

   !> 2nd block
   DO i=1, NNy
      DO j=NNy+1,Nhp
         MASSper(i,j) = MASS(i,j) + MASS(Nhp+i,j)
         STIFFper(i,j) = STIFF(i,j) + STIFF(Nhp+i,j)
         FLUXper(i,j) = FLUX(i,j) + FLUX(Nhp+i,j)
      ENDDO
   ENDDO

   !>3rd block
   DO i=NNy+1,Nhp
      DO j=1,NNy
         MASSper(i,j) = MASS(i,j) + MASS(i,Nhp+j)
         STIFFper(i,j) = STIFF(i,j) + STIFF(i,Nhp+j)
         FLUXper(i,j) = FLUX(i,j) + FLUX(i,Nhp+j)
      ENDDO
   ENDDO


END SUBROUTINE

! SUBROUTINE SurfaceIntegral(x,y,element,nop,timepoint,c,cper,w,wper,&
!    FLUXBottomc, FLUXRightc,FLUXTopc,FLUXLeftc, &
!    FLUXBottomw, FLUXRightw,FLUXTopw,FLUXLeftw,FLUXc,FLUXcper,FLUXw,FLUXwper)
!    use, intrinsic :: iso_fortran_env
!    IMPLICIT NONE
!    INCLUDE 'CHParameters.inc'


!    !> Space coordinates !xpt,ypt
!    REAL(real64), INTENT(IN) :: x(NNmax),y(NNmax)

!    !> Current element
!    INTEGER element

!    !> Local to global space numbering
!    INTEGER nop(NELmax,NNref)
!    INTEGER ngl(NNref)

!    !> timepoint
!    INTEGER timepoint

!    !> Order parameter - c - of the Cahn Hilliard Equation
!    REAL(real64) c(NNmax,Tmax)
!    REAL(real64) cper(Nhp,Tmax)

!    !> Auxiliary variable - w - the chemical potential
!    REAL(real64) w(NNmax,Tmax)
!    REAL(real64) wper(Nhp,Tmax)

!    REAL(real64) FLUXBottomc(NNmax,Tmax)
!    REAL(real64) FLUXRightc(NNmax,Tmax)
!    REAL(real64) FLUXTopc(NNmax,Tmax)
!    REAL(real64) FLUXLeftc(NNmax,Tmax)

!    REAL(real64) FLUXBottomw(NNmax,Tmax)
!    REAL(real64) FLUXRightw(NNmax,Tmax)
!    REAL(real64) FLUXTopw(NNmax,Tmax)
!    REAL(real64) FLUXLeftw(NNmax,Tmax)

!    REAL(real64) FLUXc(NNmax,Tmax)
!    REAL(real64) FLUXcper(Nhp,Tmax)
!    REAL(real64) FLUXw(NNmax,Tmax)
!    REAL(real64) FLUXwper(Nhp,Tmax)



!    !> Location indices to check for Boundary Elements
!    INTEGER BottomElement(NELx)
!    INTEGER RightElement(NELy)
!    INTEGER TopElement(NELx)
!    INTEGER LeftElement(NELy)

!    !> Boundary element control index
!    INTEGER BoundaryElement(NELmax)



!    REAL(real64) phi(NNref), tphx(NNref), tphy(NNref),phic(NNref),phie(NNref)
!    REAL(real64) gp(3), gw(3)
!    REAL(real64) Xcomputational,Xc,Xe,Ycomputational,Yc,Ye,dett

!    REAL(real64) ci,cx,cy,wi,wx,wy
!    INTEGER i,j,k,m,n



!    gw  =(/0.27777777777778, 0.444444444444, 0.27777777777778/)
!    gp =(/0.1127016654    , 0.5           , 0.8872983346    /)









!    DO i = 1,NNref
!       ngl(i) = nop(element,i)
!    ENDDO

!    !> initialization of the test functions
!    phi=0.
!    phic=0.
!    phie=0.
!    tphx=0.
!    tphy=0.





! !> calculate the contribution of the surface integral term
! IF (BoundaryElement(element)==1) THEN
!    !> check which side of the boundary the element belongs to
!    IF (BottomElement(element)==1)THEN
!       DO j=1,3
!          CALL TestFunctions(gp(j),0._real64,phi, phic, phie)

!          Xcomputational=0.
!          Xc=0.
!          Xe=0.
!          Ycomputational=0.
!          Yc=0.
!          Ye=0.
!          DO n=1,NNref
!             Xcomputational= Xcomputational + x(ngl(n)) * phi(n)
!             Xc= Xc + x(ngl(n)) * phic(n)
!             Xe= Xe + x(ngl(n)) * phie(n)
!             Ycomputational= Ycomputational + y(ngl(n)) * phi(n)
!             Yc= Yc + y(ngl(n)) * phic(n)
!             Ye= Ye + y(ngl(n)) * phie(n)
!          ENDDO
!          dett=Xc*Ye-Xe*Yc

!          DO i=1,NNref
!             tphx(i)=(Ye*phic(i)-Yc*phie(i))/dett
!             tphy(i)=(Xc*phie(i)-Xe*phic(i))/dett
!          ENDDO


!          ci = 0.
!          cx = 0.
!          cy = 0.
!          DO i=1,NNref
!             ci = ci + cper(ngl(i),timePoint)*phi(i)
!             cx = cx + cper(ngl(i),timePoint)*tphx(i)
!             cy = cy + cper(ngl(i),timePoint)*tphy(i)
!          ENDDO

!          wi = 0.
!          wx = 0.
!          wy = 0.
!          DO i=1,NNref
!             wi = wi + wper(ngl(i),timePoint)*phi(i)
!             wx = wx + wper(ngl(i),timePoint)*tphx(i)
!             wy = wy + wper(ngl(i),timePoint)*tphy(i)
!          ENDDO

!             DO m=1,7,3
!                FLUXBottomc(ngl(m),timepoint) = FLUXBottomc(ngl(m),timepoint) &
!                - gw(j)*Xc*phi(m)*cy
!                FLUXBottomw(ngl(m),timepoint) = FLUXBottomw(ngl(m),timepoint) &
!                - gw(j)*Xc*phi(m)*wy
!             ENDDO
!       ENDDO
!    ELSEIF (RightElement(element)==1)THEN
!       DO j=1,3
!          CALL TestFunctions(1._real64,gp(j),phi, phic, phie)

!          Xcomputational=0.
!          Xc=0.
!          Xe=0.
!          Ycomputational=0.
!          Yc=0.
!          Ye=0.
!          DO n=1,NNref
!             Xcomputational= Xcomputational + x(ngl(n)) * phi(n)
!             Xc= Xc + x(ngl(n)) * phic(n)
!             Xe= Xe + x(ngl(n)) * phie(n)
!             Ycomputational= Ycomputational + y(ngl(n)) * phi(n)
!             Yc= Yc + y(ngl(n)) * phic(n)
!             Ye= Ye + y(ngl(n)) * phie(n)
!          ENDDO
!          dett=Xc*Ye-Xe*Yc

!          DO i=1,NNref
!             tphx(i)=(Ye*phic(i)-Yc*phie(i))/dett
!             tphy(i)=(Xc*phie(i)-Xe*phic(i))/dett
!          ENDDO


!          ci = 0.
!          cx = 0.
!          cy = 0.
!          DO i=1,NNref
!             ci = ci + cper(ngl(i),timePoint)*phi(i)
!             cx = cx + cper(ngl(i),timePoint)*tphx(i)
!             cy = cy + cper(ngl(i),timePoint)*tphy(i)
!          ENDDO

!          wi = 0.
!          wx = 0.
!          wy = 0.
!          DO i=1,NNref
!             wi = wi + wper(ngl(i),timePoint)*phi(i)
!             wx = wx + wper(ngl(i),timePoint)*tphx(i)
!             wy = wy + wper(ngl(i),timePoint)*tphy(i)
!          ENDDO

!             DO m=7,9
!                FLUXRightc(ngl(m),timepoint) = FLUXRightc(ngl(m),timepoint) &
!                + gw(j)*Ye*phi(m)*cx
!                FLUXRightw(ngl(m),timepoint) = FLUXRightw(ngl(m),timepoint) &
!                + gw(j)*Ye*phi(m)*wx
!             ENDDO
!       ENDDO
!    ELSEIF (TopElement(element)==1)THEN
!          DO j=1,3
!             CALL TestFunctions(gp(j),1._real64,phi, phic, phie)

!             Xcomputational=0.
!             Xc=0.
!             Xe=0.
!             Ycomputational=0.
!             Yc=0.
!             Ye=0.
!             DO n=1,NNref
!                Xcomputational= Xcomputational + x(ngl(n)) * phi(n)
!                Xc= Xc + x(ngl(n)) * phic(n)
!                Xe= Xe + x(ngl(n)) * phie(n)
!                Ycomputational= Ycomputational + y(ngl(n)) * phi(n)
!                Yc= Yc + y(ngl(n)) * phic(n)
!                Ye= Ye + y(ngl(n)) * phie(n)
!             ENDDO
!             dett=Xc*Ye-Xe*Yc

!             DO i=1,NNref
!                tphx(i)=(Ye*phic(i)-Yc*phie(i))/dett
!                tphy(i)=(Xc*phie(i)-Xe*phic(i))/dett
!             ENDDO


!             ci = 0.
!             cx = 0.
!             cy = 0.
!             DO i=1,NNref
!                ci = ci + cper(ngl(i),timePoint)*phi(i)
!                cx = cx + cper(ngl(i),timePoint)*tphx(i)
!                cy = cy + cper(ngl(i),timePoint)*tphy(i)
!             ENDDO

!             wi = 0.
!             wx = 0.
!             wy = 0.
!             DO i=1,NNref
!                wi = wi + wper(ngl(i),timePoint)*phi(i)
!                wx = wx + wper(ngl(i),timePoint)*tphx(i)
!                wy = wy + wper(ngl(i),timePoint)*tphy(i)
!             ENDDO

!                DO m=3,9,3
!                   FLUXTopc(ngl(m),timepoint) = FLUXTopc(ngl(m),timepoint) &
!                   + gw(j)*Xc*phi(m)*cy
!                   FLUXTopw(ngl(m),timepoint) = FLUXTopw(ngl(m),timepoint) &
!                   + gw(j)*Xc*phi(m)*wy
!                ENDDO
!          ENDDO
!    ELSEIF (LeftElement(element)==1)THEN
!             DO j=1,3
!                CALL TestFunctions(0._real64,gp(j),phi, phic, phie)

!                Xcomputational=0.
!                Xc=0.
!                Xe=0.
!                Ycomputational=0.
!                Yc=0.
!                Ye=0.
!                DO n=1,NNref
!                   Xcomputational= Xcomputational + x(ngl(n)) * phi(n)
!                   Xc= Xc + x(ngl(n)) * phic(n)
!                   Xe= Xe + x(ngl(n)) * phie(n)
!                   Ycomputational= Ycomputational + y(ngl(n)) * phi(n)
!                   Yc= Yc + y(ngl(n)) * phic(n)
!                   Ye= Ye + y(ngl(n)) * phie(n)
!                ENDDO
!                dett=Xc*Ye-Xe*Yc

!                DO i=1,NNref
!                   tphx(i)=(Ye*phic(i)-Yc*phie(i))/dett
!                   tphy(i)=(Xc*phie(i)-Xe*phic(i))/dett
!                ENDDO


!                ci = 0.
!                cx = 0.
!                cy = 0.
!                DO i=1,NNref
!                   ci = ci + cper(ngl(i),timePoint)*phi(i)
!                   cx = cx + cper(ngl(i),timePoint)*tphx(i)
!                   cy = cy + cper(ngl(i),timePoint)*tphy(i)
!                ENDDO

!                wi = 0.
!                wx = 0.
!                wy = 0.
!                DO i=1,NNref
!                   wi = wi + wper(ngl(i),timePoint)*phi(i)
!                   wx = wx + wper(ngl(i),timePoint)*tphx(i)
!                   wy = wy + wper(ngl(i),timePoint)*tphy(i)
!                ENDDO

!                   DO m=1,3
!                      FLUXLeftc(ngl(m),timepoint) = FLUXLeftc(ngl(m),timepoint) &
!                      - gw(j)*Ye*phi(m)*cx
!                      FLUXLeftw(ngl(m),timepoint) = FLUXLeftw(ngl(m),timepoint) &
!                      - gw(j)*Ye*phi(m)*wx
!                   ENDDO
!             ENDDO


!       ENDIF
!    ENDIF

!     !> Summation of the components to the total FLUX
!    DO k=1,NNref
!       FLUXcper(ngl(k),timepoint) = FLUXBottomc(ngl(k),timepoint) + FLUXRightc(ngl(k),timepoint) &
!                                 + FLUXTopc(ngl(k),timepoint) + FLUXLeftc(ngl(k),timepoint)

!       FLUXwper(ngl(k),timepoint) = FLUXBottomw(ngl(k),timepoint) + FLUXRightw(ngl(k),timepoint) &
!                                 + FLUXTopw(ngl(k),timepoint) + FLUXLeftw(ngl(k),timepoint)
!    ENDDO







! END SUBROUTINE SurfaceIntegral









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


!>******************************************************************************
!> PotentialTerm evaluates the entries of the P(c,t) operator and the DP(c,t) matrix
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


SUBROUTINE PotentialTerm(x,y,element,nop,timepoint,c,P,DP)
   use, intrinsic :: iso_fortran_env
   IMPLICIT NONE
   INCLUDE 'CHParameters.inc'

   
   
   !> Space coordinates !xpt,ypt
   REAL(real64), INTENT(IN) :: x(NNmax),y(NNmax)

   !> Current element
   INTEGER element

   !> Local to global space numbering
   INTEGER nop(NELmax,NNref)
   INTEGER ngl(NNref)

   !> timepoint
   INTEGER timepoint

   !> Order parameter - c - of the Cahn Hilliard Equation
   REAL(real64) c(NNmax,Tmax)


   !> non linear term
   REAL(real64) P(NNmax,Tmax)

   !> Derivative of the non-linear term
   REAL(real64) DP(NNmax,NNmax,Tmax)




   REAL(real64) phi(NNref), tphx(NNref), tphy(NNref),phic(NNref),phie(NNref)
   REAL(real64) gp(3), gw(3)
   REAL(real64) Xcomputational,Xc,Xe,Ycomputational,Yc,Ye,dett

   REAL(real64) ci,cx,cy
   INTEGER i,j,k,l,m,n



   gw  =(/0.27777777777778, 0.444444444444, 0.27777777777778/)
   gp =(/0.1127016654    , 0.5           , 0.8872983346    /)

   DO i = 1,NNref
      ngl(i) = nop(element,i)
   ENDDO

   phi=0.0_real64
   phic=0.0_real64
   phie=0.0_real64
   tphx=0.0_real64
   tphy=0.0_real64
   ! Loop over qauss points
   DO j = 1,3
      DO k = 1,3

         call TestFunctions(gp(j),gp(k),phi,phic,phie)

         ! Defines the computational domain coordinates and the 2-dimensional Jacobian dett
         Xcomputational=0.0_real64
         Xc=0.0_real64
         Xe=0.0_real64
         Ycomputational=0.0_real64
         Yc=0.0_real64
         Ye=0.0_real64
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


         ci = 0.0_real64
         cx = 0.0_real64
         cy = 0.0_real64
         DO i=1,NNref
            ci = ci + c(ngl(i),timePoint)*phi(i)
            cx = cx + c(ngl(i),timePoint)*tphx(i)
            cy = cy + c(ngl(i),timePoint)*tphy(i)
         ENDDO

         DO l=1,NNref
            DO m=1,NNref
                DP(ngl(l),ngl(m),timepoint) =  DP(ngl(l),ngl(m),timepoint) + gw(j)*gw(k)*dett*((3.0_real64)*(ci**2) - 1.0_real64)*phi(l)*phi(m)
            ENDDO
            P(ngl(l),timepoint) = P(ngl(l),timepoint) + gw(j)*gw(k)*dett*(ci**3 - ci)*phi(l)
         ENDDO

      ENDDO
   ENDDO


END SUBROUTINE PotentialTerm


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
SUBROUTINE StationarySystem(timePoint,c,w,STIFF,MASS,P,DP,FLUX,&
   fjac,fvec)
   use, intrinsic :: iso_fortran_env
   IMPLICIT NONE
   INCLUDE 'CHParameters.inc'

   !> Current time step
   INTEGER timepoint

   !> Order parameter - c - of the Cahn Hilliard Equation
   REAL(real64) c(NNmax,Tmax)
   
   !> Auxiliary variable - w - the chemical potential
   REAL(real64) w(NNmax,Tmax)
  

   !> System matrices
   !> STIFF - Stiffness
   !> MASS - Mass
   REAL(real64) STIFF(NNmax,NNmax)
   
   REAL(real64) MASS(NNmax,NNmax)
   

   !> non linear term
   REAL(real64) P(NNmax,Tmax)
   
   !> Derivative of the non-linear term
   REAL(real64) DP(NNmax,NNmax,Tmax)
   


   !> The Surface FLUX integral term
   REAL(real64) FLUX(NNmax,NNmax)
   


   !> Jacobian matix of Newton system
   REAL(real64) fjac(2*NNmax,2*NNmax)
   !> Righ handside vector of Newton system
   REAL(real64) fvec(2*NNmax)


   INTEGER i,j

   REAL(real64) e2

   e2 = (epsilon)**2


   !*************************************
   !> CAHN HILLIARD
   !*************************************
   DO i=1,NNmax
      DO j=1,NNmax
         fvec(i)= fvec(i) + MASS(i,j)*w(j,timepoint) - (e2)*(STIFF(i,j))*c(j,timepoint)
         fjac(i,j)= fjac(i,j) + MASS(i,j)
         
      ENDDO
      fvec(i) = fvec(i)  - P(i,timepoint)

      DO j= NNmax+1,2*NNmax
         fjac(i,j)=  fjac(i,j) - (e2)*(STIFF(i,j-NNmax)) - DP(i,j-NNmax,timepoint)
        
      ENDDO
   ENDDO

   DO i=NNmax+1,2*NNmax
      DO j=1,NNmax
         fvec(i)= fvec(i) - MASS(i-NNmax,j)*c(j,timePoint-1)+ (deltat)*(STIFF(i-NNmax,j))*w(j,timepoint) + MASS(i-NNmax,j)*c(j,timePoint) 
        
         fjac(i,j)= fjac(i,j) + (deltat)*(STIFF(i-NNmax,j))
          
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
   INCLUDE 'CHParameters.inc'

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


subroutine solve_system(n, A, b, x)
   use, intrinsic :: iso_fortran_env
   IMPLICIT NONE
   EXTERNAL           DGETRF, DGETRS
   INCLUDE 'CHParameters.inc'

   integer, intent(in) :: n
   real(real64), intent(inout) :: A(n,n)
   real(real64), intent(inout) :: b(n)
   real(real64), intent(out) :: x(n)
   integer :: pivot(n)
   integer :: info

   ! LU decomposition with pivoting
   call DGETRF(n, n, A, n, pivot, info)
   if (info /= 0) then
      write(*,*) "LU decomposition failed.",info
      return
   end if

   ! Solve using LU decomposition
   call DGETRS('N', n, 1, A, n, pivot, b, n, info)
   if (info /= 0) then
      write(*,*) "Solution failed."
      return
   end if

   x = b
end subroutine solve_system



subroutine gausselim(n, A, b, x)
   use, intrinsic :: iso_fortran_env
   IMPLICIT NONE
   INCLUDE 'CHParameters.inc'

   integer, intent(in) :: n
   real(real64), intent(inout) :: A(n,n)
   real(real64), intent(inout) :: b(n)
   real(real64), intent(out) :: x(n)
   integer :: i, j, k, max_row
   real(real64) :: max_val, temp, factor

   ! Forward elimination with partial pivoting
   do k = 1, n-1
      max_row = k
      max_val = abs(A(k,k))
      do i = k+1, n
         if (abs(A(i,k)) > max_val) then
            max_row = i
            max_val = abs(A(i,k))
         end if
      end do

      ! Swap rows to ensure non-zero pivot element
      if (max_row /= k) then
         do j = k, n
            temp = A(k,j)
            A(k,j) = A(max_row,j)
            A(max_row,j) = temp
         end do
         temp = b(k)
         b(k) = b(max_row)
         b(max_row) = temp
      end if

      do i = k+1, n
         factor = A(i,k) / A(k,k)
         do j = k+1, n
            A(i,j) = A(i,j) - factor * A(k,j)
         end do
         b(i) = b(i) - factor * b(k)
      end do
   end do

   ! Back substitution
   x(n) = b(n) / A(n,n)
   do i = n-1, 1, -1
      x(i) = b(i)
      do j = i+1, n
         x(i) = x(i) - A(i,j) * x(j)
      end do
      x(i) = x(i) / A(i,i)
   end do
end subroutine gausselim

