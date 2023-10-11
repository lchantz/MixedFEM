!>******************************************************************************
!>    The form of the Cahn Hilliard equation is...
!>******************************************************************************
program NEUMANN
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
 
    !> Location indices to check for Boundary nodes
    INTEGER BottomNode(NNmax)
    INTEGER RightNode(NNmax)
    INTEGER TopNode(NNmax)
    INTEGER LeftNode(NNmax)
 
    !> Boundary element control index
    INTEGER BoundaryNode(NNmax)
 
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
 
 
    !> ***********************************************************
    !> The discrete approximation matrices of the weak formulation
    !> ************************************************************
 
    !> STIFF - Stiffness
    REAL(real64) STIFF(NNmax,NNmax)
  
    !> MASS - Mass
    REAL(real64) MASS(NNmax,NNmax)
    REAL(real64) PREC(2*NNmax,2*NNmax)
    REAL(real64) PRECinv(2*NNmax,2*NNmax)
 
 
    !> f - Non linear term from inner product of the double well potential c^3-c
    REAL(real64) f(NNmax,Tmax)
 
 
    !> Jacobian of the non-linear term from double well potential
    REAL(real64) J(NNmax,NNmax,Tmax)
 

 
 
 
    !> Newton tolerance parameters
    REAL(real64),PARAMETER:: tolS=1.0E-6_real64, tol = 1.0E-6_real64
    INTEGER, PARAMETER :: ntrial = 100
    INTEGER, PARAMETER :: max_iter = 1000
 
    !> Stationary system
    REAL(real64) fvec(2*NNmax)
    REAL(real64) fjac(2*NNmax,2*NNmax)
 
 
    !> Newton solution vector
    REAL(real64) dSol(2*NNmax)
    REAL(real64) Sol(2*NNmax)
 
 
    INTEGER element,i,timepoint,trial,q,p
    LOGICAL CONVERGENT

    OPEN(unit=1,file='T=0_NEUMANN.dat')

    OPEN(unit=2,file='T=20_NEUMANN.dat')
    OPEN(unit=3,file='T=30_NEUMANN.dat')
    OPEN(unit=4,file='T=40_NEUMANN.dat')
    OPEN(unit=5,file='T=50_NEUMANN.dat')
   !  OPEN(unit=8,file='T=80_NEUMANN.dat')
   !  OPEN(unit=9,file='T=90_NEUMANN.dat')
 
   
   !  OPEN(unit=20,file='T=200_NEUMANN.dat')
   !  OPEN(unit=30,file='T=300_NEUMANN.dat')
   !  OPEN(unit=40,file='T=400_NEUMANN.dat')
   !  OPEN(unit=50,file='T=500_NEUMANN.dat')
   !  OPEN(unit=80,file='T=800_NEUMANN.dat')
   !  OPEN(unit=90,file='T=900_NEUMANN.dat')

   !  OPEN(unit=200,file='T=2000_NEUMANN.dat')
   !  OPEN(unit=300,file='T=3000_NEUMANN.dat')
   !  OPEN(unit=400,file='T=4000_NEUMANN.dat')
   !  OPEN(unit=500,file='T=5000_NEUMANN.dat')
    ! OPEN(unit=800,file='T=8000_NEUMANN.dat')
    ! OPEN(unit=900,file='T=9000_NEUMANN.dat')


 
 
    WRITE(*,'(25x,A)') 'Results c(n,t)'
    CALL SpaceDiscretization(x,y)
    CALL GlobalNodeNumbering(nop)
 
    !> Initialization of indices
    BottomNode=0
    RightNode=0
    TopNode=0
    LeftNode=0
    BoundaryNode=0
    !> Mark boundary nodes
    CALL BoundaryNodes(x,y,BottomNode,RightNode,TopNode,LeftNode,BoundaryNode)
 
    !> Initialization of indices
    BottomElement=0
    RightElement=0
    TopElement=0
    LeftElement=0
    BottomElement=0
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
    ENDDO

    !> PRECONDITIONER SET UP
!     PREC=0.0
!     DO c1=1,NNmax
!       DO c2=1,NNmax
!          PREC(c1,c2) = PREC(c1,c2) + MASS(c1,c2)
!       ENDDO 
!       DO c2=NNmax+1,2*NNmax 
!          PREC(c1,c2) = PREC(c1,c2) -(e)**2 *STIFF(c1,c2-NNmax)
!       ENDDO
!    ENDDO
   
!    DO c1=NNmax+1,2*NNmax 
!       DO c2=1,NNmax
!          PREC(c1,c2) = PREC(c1,c2) + (dt)*STIFF(c1-NNmax,c2)
!       ENDDO 
!       DO c2=NNmax+1,2*NNmax 
!          PREC(c1,c2) = PREC(c1,c2) + MASS(c1-NNmax,c2-NNmax)
!       ENDDO
!    ENDDO   
!    PRECinv=0.0
!   CALL inv(PREC,PRECinv,2*NNmax)
 
 
    !>****************************************************************
    !> Set the initial conditions for the c,w variables of the problem
    !> For timepoint = 1
    !>*****************************************************************
    ! timepoint = 1
    c=0.0_real64
    w=0.0_real64
    CALL InitialConditions(x,y,c,w)
    DO q=1,NNy
      DO p=1,NNx
       WRITE(1,*) q,p,c((p-1)*NNy+q,1)
      ENDDO  
    ENDDO
    CLOSE(1)
    !> Calculate potential at t=1, to start the forward in time solver
          f=0.0_real64
          J=0.0_real64
          DO element=1,NELmax
             CALL PotentialTerm(x,y,element,nop,1,c,f,J)
          ENDDO
 
 !  !>******************************************************************************
 !    !> Given an initial guess x for a root in n dimensions, take ntrial Newton-Raphson
 !  !> steps to improve the root. Stop if the root converges in either summed absolute
 !    !> variable increments tolx or summed absolute function values tolf.
 !  !>******************************************************************************
 
 
 !> Time loop to attain the solution at every timepoint
   !  DO timepoint=2,Tmax
   !     print *, "Time point: ", timepoint
 
   !     trial=0
   !     CONVERGENT=.false.
   !     Sol=0.0_real64
 
   !     !> Choose the initial guess for the first newton iteration
   !     Sol(1:NNmax) = w(:,timepoint-1) !+ (w(:,timepoint) - w(:,timepoint-1)) * (dt)
   !     Sol(NNmax+1:2*NNmax) = c(:,timepoint-1) !+ (c(:,timepoint) - c(:,timepoint-1)) * (dt)
 
 
   !     ! !> Start the Newton iteration loop
   !     ! NewtonTrials: DO trial=1,ntrial
   !     DO WHILE (.NOT.CONVERGENT .AND. trial < ntrial)
   !        dSol=0.0_real64
 
   !        !>*****************************************************************
   !        !> Calculate the double well potential and the derivative
   !        !> Depend on c, the order parameter
   !        !>*****************************************************************
   !        f=0.0_real64
   !        J=0.0_real64
   !        DO element=1,NELmax
   !           CALL PotentialTerm(x,y,element,nop,timepoint,c,f,J)
   !        ENDDO
 
   !        !>******************************************************************
   !        !> Set the Stationary System to be solved at every Newton iteration
   !        !> ___________________Fjac.dx = -Fvec_______________________________
   !        !>******************************************************************
   !        fvec=0.0_real64
   !        fjac=0.0_real64
   !        CALL StationarySystem(timePoint,c,w,STIFF,MASS,f,J,fjac,fvec)
 
 
   !        !> Solve the system - find the correction dSol
   !        fvec = -fvec
   !       !  DO i =1,2*NNmax
   !       !  fvec(i) = dot_product(PRECinv(i,:),fvec)
   !       !  fjac=matmul(PRECinv,fjac)
   !       !  ENDDO
   !       !  CALL LU_decomposition(fjac,2*NNmax)
   !        CALL SolveFJACDSOL( tolS, max_iter, fvec, fjac, dSol)
          
 
   !        !> Update the solution
   !        Sol(:) = Sol(:) + dSol(:)
 
 
   !        !> Check for convergence
   !        CONVERGENT = (NORM2(dSol) < tolS)
   !        print*,"NORM2(dsol): ", NORM2(dSol)
   !       !  print*,"CONVERGENCE", CONVERGENT
 
   !        trial = trial + 1
   !        w(:,timepoint) = dSol(1:NNmax)
   !        c(:,timepoint) = dSol(NNmax+1:2*NNmax)
 
   !     ENDDO  !trials
 
   !     IF (CONVERGENT) THEN
   !        write(*,*) 'Newton-Raphson converged in', trial, 'iterations'
   !     ELSE
   !        write(*,*) 'Newton-Raphson did not converge within', ntrial, 'iterations'
   !     ENDIF
 
   !  ENDDO !> Time loop


   !>*************************************************************
   !> FORWARD EULER & PROVATAS/ELDER 2010
   !>*************************************************************  
    DO timepoint=2,Tmax
      print*,'timepoint',timepoint
      Sol=0.0_real64
   !>*****************************************************************
   !> Calculate the double well potential and the derivative
   !> Depend on c, the order parameter
   !>*****************************************************************
          f=0.0_real64
          J=0.0_real64
          DO element=1,NELmax
             CALL PotentialTerm(x,y,element,nop,timepoint,c,f,J)
          ENDDO

          fvec=0.0_real64
          fjac=0.0_real64
         CALL StationarySystem(timePoint,c,w,STIFF,MASS,f,J,fjac,fvec)
 
 
          !> Solve the system - find the correction dSol
           
          Sol = BiCGStab(fjac,fvec)

           w(:,timepoint) = Sol(1:NNmax)
           c(:,timepoint) = Sol(NNmax+1:2*NNmax)

         ENDDO !TIME LOOP     

 
    !> Write the results
    DO q=1,NNy
      DO p=1,NNx
      WRITE(2,*) q,p,c((p-1)*NNy+q,200)
       WRITE(3,*) q,p,c((p-1)*NNy+q,300)
       WRITE(4,*) q,p,c((p-1)*NNy+q,400)
       WRITE(5,*) q,p,c((p-1)*NNy+q,490)
      !  WRITE(8,'(2x,E16.9,1x,E16.9,1x,E16.9,1x,E16.9)') x(i),y(i),c(i,80)
      !  WRITE(9,'(2x,E16.9,1x,E16.9,1x,E16.9,1x,E16.9)') x(i),y(i),c(i,90)

      !  WRITE(20,'(2x,E16.9,1x,E16.9,1x,E16.9,1x,E16.9)') x(i),y(i),c(i,200)
      !  WRITE(30,'(2x,E16.9,1x,E16.9,1x,E16.9,1x,E16.9)') x(i),y(i),c(i,300)
      !  WRITE(40,'(2x,E16.9,1x,E16.9,1x,E16.9,1x,E16.9)') x(i),y(i),c(i,400)
      !  WRITE(50,'(2x,E16.9,1x,E16.9,1x,E16.9,1x,E16.9)') x(i),y(i),c(i,500)
      !  WRITE(80,'(2x,E16.9,1x,E16.9,1x,E16.9,1x,E16.9)') x(i),y(i),c(i,800)
      !  WRITE(90,'(2x,E16.9,1x,E16.9,1x,E16.9,1x,E16.9)') x(i),y(i),c(i,900)

      !  WRITE(200,'(2x,E16.9,1x,E16.9,1x,E16.9,1x,E16.9)') x(i),y(i),c(i,2000)
      !  WRITE(300,'(2x,E16.9,1x,E16.9,1x,E16.9,1x,E16.9)') x(i),y(i),c(i,3000)
      !  WRITE(400,'(2x,E16.9,1x,E16.9,1x,E16.9,1x,E16.9)') x(i),y(i),c(i,4000)
      !  WRITE(500,'(2x,E16.9,1x,E16.9,1x,E16.9,1x,E16.9)') x(i),y(i),c(i,5000)
    !    WRITE(800,'(2x,E16.9,1x,E16.9,1x,E16.9,1x,E16.9)') x(i),y(i),c(i,8000)
    !    WRITE(900,'(2x,E16.9,1x,E16.9,1x,E16.9,1x,E16.9)') x(i),y(i),c(i,9000)
      ENDDO
    ENDDO
 
    CLOSE(2)
    CLOSE(3)
    CLOSE(4)
    CLOSE(5)
   !  CLOSE(8)
   !  CLOSE(9)

   !  CLOSE(20)
   !  CLOSE(30)
   !  CLOSE(40)
   !  CLOSE(50)
   !  CLOSE(80)
   !  CLOSE(90)

   !  CLOSE(200)
   !  CLOSE(300)
   !  CLOSE(400)
   !  CLOSE(500)
    ! CLOSE(800)
    ! CLOSE(900)

END program NEUMANN
 
 
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
    REAL(real64) dx,dy
    INTEGER i,j,NNode
 
    x=0.0_real64
    y=0.0_real64
 
    xorigin=0.0_real64
    xlast=1.0_real64
    dx=(xlast-xorigin)/real(2*NELx)
 
    yorigin=0.0_real64
    ylast=1.0_real64
    dy=(ylast-yorigin)/real(2*NELy)
 
    !> Setting the coordinates of each node
    x(1)=xorigin
    y(1)=yorigin
    DO i=1,NNx
       NNode = (i-1)*NNy + 1 ! Ascending node numbering from bottom left up to right
       x(NNode) = x(1) + (i-1)*dx
       y(NNode) = y(1)
       DO j=2,NNy
          x(NNode+j-1)=x(NNode)
          y(NNode+j-1)=y(NNode)+(j-1)*dy
 
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
 
 SUBROUTINE BoundaryNodes(x,y,BottomNode,RightNode,TopNode,LeftNode,BoundaryNode)
    use, intrinsic :: iso_fortran_env
    IMPLICIT NONE
    INCLUDE 'CHParameters.inc'
 
    !> Space coordinates
    REAL(real64) x(NNmax),y(NNmax)
 
    !> Location indices to check for Boundary nodes
    INTEGER BottomNode(NNmax)
    INTEGER RightNode(NNmax)
    INTEGER TopNode(NNmax)
    INTEGER LeftNode(NNmax)
 
    !> Boundary element control index
    INTEGER BoundaryNode(NNmax)
    INTEGER i
 
 
 
    !> Check for boundary nodes
    DO i=1,NNmax
       IF (y(i)==0.0_real64) THEN !bottom side nodes
          BottomNode(i)=1
          BoundaryNode(i)=1
       ELSEIF (x(i)==1.0_real64) THEN !right side nodes
          RightNode(i)=1
          BoundaryNode(i)=1
       ELSEIF (y(i)==1.0_real64) THEN !top side nodes
          TopNode(i)=1
          BoundaryNode(i)=1
       ELSEIF (x(i)==0.0_real64) THEN !left side nodes
          LeftNode(i)=1
          BoundaryNode(i)=1
       ENDIF
    ENDDO
    ! !> Check for boundary nodes
    ! DO i=1,(NNx-1)*NNy+1,NNy !bottom side nodes
    !    BottomNode(i)=1
    !    BoundaryNode(i)=1
    ! ENDDO
    ! DO i=(NNx-1)*NNy+1,NNmax,1 !right side nodes
    !    RightNode(i)=1
    !    BoundaryNode(i)=1
    ! ENDDO
    ! DO i=NNmax,NNy,-NNy !top side nodes
    !    TopNode(i)=1
    !    BoundaryNode(i)=1
    ! ENDDO
    ! DO i=NNy,1,-1 !left side nodes
    !    LeftNode(i)=1
    !    BoundaryNode(i)=1
    ! ENDDO
 
 END SUBROUTINE
 
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
 !>******************************************************************************
 SUBROUTINE TimeDiscretization(t)
    use, intrinsic :: iso_fortran_env
    IMPLICIT NONE
    INCLUDE 'CHParameters.inc'
 
    !> Time coordinates
    REAL(real64) t(Tmax)
 
    REAL(real64)  Tinit
 
    INTEGER j
 
    Tinit=0.0_real64
 
 
    t(1)=Tinit
    DO j=2,Tmax
       t(j) = (j-1)*dt
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
 SUBROUTINE InitialConditions(x, y, c, w)
    use, intrinsic :: iso_fortran_env
    IMPLICIT NONE
    INCLUDE 'CHParameters.inc'
 
    !> Space coordinates
    REAL(real64), DIMENSION(NNmax) :: x, y
 
    !> Order parameter - c - of the Cahn Hilliard Equation
    REAL(real64), DIMENSION(NNmax, Tmax) :: c
 
    !> Auxiliary variable - w - the chemical potential
    REAL(real64), DIMENSION(NNmax, Tmax) :: w
 
    INTEGER :: i

 
 
    DO i = 1, NNmax
      
      ! c(i,1) = 0.5_real64*(1.0_real64-dcos(4.0_real64*pi*x(i)))*(1.0_real64-dcos(2.0_real64*pi*y(i))) - 1.0_real64
      ! w(i,1) = 0.0_real64

      CALL RANDOM_NUMBER(c(i,1))
      c(i,1) = 0.63_real64 + 0.02_real64*(0.5_real64-c(i,1))
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
    REAL(real64) Xdomain,Xc,Xe,Ydomain,Yc,Ye,dett
 
 
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
 
          ! Defines the domain domain coordinates and the 2-dimensional Jacobian dett
          Xdomain=0.0_real64
          Xc=0.0_real64
          Xe=0.0_real64
          Ydomain=0.0_real64
          Yc=0.0_real64
          Ye=0.0_real64
          DO n=1,NNref
             Xdomain= Xdomain + x(ngl(n)) * phi(n)
             Xc= Xc + x(ngl(n)) * phic(n)
             Xe= Xe + x(ngl(n)) * phie(n)
             Ydomain= Ydomain + y(ngl(n)) * phi(n)
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
 
 
 SUBROUTINE PotentialTerm(x,y,element,nop,timepoint,c,f,J)
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
    REAL(real64) f(NNmax,Tmax)
 
    !> Derivative of the non-linear term
    REAL(real64) J(NNmax,NNmax,Tmax)
 
 
 
 
    REAL(real64) phi(NNref), tphx(NNref), tphy(NNref),phic(NNref),phie(NNref)
    REAL(real64) gp(3), gw(3)
    REAL(real64) Xdomain,Xc,Xe,Ydomain,Yc,Ye,dett
 
    REAL(real64) ci,cx,cy
    INTEGER i,r,k,l,m,n
 
 
 
    gw  =(/0.27777777777778, 0.444444444444, 0.27777777777778/)
    gp =(/0.1127016654    , 0.5           , 0.8872983346    /)
 
    DO i = 1,NNref
       ngl(i) = nop(element,i)
    ENDDO
 
    ! Loop over qauss points
    DO r = 1,3
       DO k = 1,3
 
          call TestFunctions(gp(r),gp(k),phi,phic,phie)
 
          ! Defines the domain domain coordinates and the 2-dimensional Jacobian dett
          Xdomain=0.0_real64
          Xc=0.0_real64
          Xe=0.0_real64
          Ydomain=0.0_real64
          Yc=0.0_real64
          Ye=0.0_real64
          DO n=1,NNref
             Xdomain= Xdomain + x(ngl(n)) * phi(n)
             Xc= Xc + x(ngl(n)) * phic(n)
             Xe= Xe + x(ngl(n)) * phie(n)
             Ydomain= Ydomain + y(ngl(n)) * phi(n)
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
             ci = ci + c(ngl(i),timePoint)*phi(i) !calculated at t and not t-1 for forward euler
             cx = cx + c(ngl(i),timePoint)*tphx(i)
             cy = cy + c(ngl(i),timePoint)*tphy(i)
          ENDDO
          
         !*********************************
         !W = 100c**2(1-c)**2
         !*********************************  
         !  DO l=1,NNref
         !    f(ngl(l),timepoint) = f(ngl(l),timepoint) + gw(r)*gw(k)*dett*(2.0E+2_real64)*(ci-ci**2)*( 1.0_real64 - 2.0_real64*ci)*phi(l)
         !     DO m=1,NNref
         !        J(ngl(l),ngl(m),timepoint) =  J(ngl(l),ngl(m),timepoint) + gw(r)*gw(k)*dett*(2.0E+2_real64)*(6.0_real64*ci**2 - 6.0_real64*ci +1.0_real64)*phi(l)*phi(m)
         !     ENDDO
         !  ENDDO

         !***************************************************
         !W = 0.25(c**2-1)**2 ginzburg landau 
         !*************************************************** 
          DO l=1,NNref
            f(ngl(l),timepoint) = f(ngl(l),timepoint) + gw(r)*gw(k)*dett*(ci**3 - ci)*phi(l)
            DO m=1,NNref
               J(ngl(l),ngl(m),timepoint) =  J(ngl(l),ngl(m),timepoint) + gw(r)*gw(k)*dett*((3.0_real64)*(ci**2) - 1.0_real64)*phi(l)*phi(m)
            ENDDO
         ENDDO
         
         !****************************************************
         !W = 0.25(c**2-1)**2 ginzburg landau & convex splitting
         !**************************************************** 
         !  DO l=1,NNref
         !    f(ngl(l),timepoint) = f(ngl(l),timepoint) + gw(r)*gw(k)*dett*(ci**3)*phi(l)
         !    DO m=1,NNref
         !       J(ngl(l),ngl(m),timepoint) =  J(ngl(l),ngl(m),timepoint) + gw(r)*gw(k)*dett*((3.0_real64)*(ci**2))*phi(l)*phi(m)
         !    ENDDO
         ! ENDDO
 
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
 SUBROUTINE StationarySystem(timePoint,c,w,STIFF,MASS,f,J,fjac,fvec)
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
    REAL(real64) f(NNmax,Tmax)
 
    !> Derivative of the non-linear term
    REAL(real64) J(NNmax,NNmax,Tmax)

 
 
    !> Jacobian matix of Newton system
    REAL(real64) fjac(2*NNmax,2*NNmax)
    !> Righ handside vector of Newton system
    REAL(real64) fvec(2*NNmax)
 
 
    INTEGER m,n
 
 
 
    !*******************************************************************
    !> CAHN HILLIARD first order semi implicit  theta = 1
    !*******************************************************************
   !  DO m=1,NNmax
   !    DO n=1,NNmax
         
   !       fvec(m) = fvec(m)  + MASS(m,n)*w(n,timepoint) + MASS(m,n)*c(n,timepoint-1) - (e**2)*STIFF(m,n)*c(n,timepoint)
         
   !       fjac(m,n) = fjac(m,n) + MASS(m,n)

   !     ENDDO

   !       fvec(m) = fvec(m) - f(m,timepoint)

   !       DO n= NNmax+1,2*NNmax

   !        fjac(m,n)=  fjac(m,n) - (e**2)*STIFF(m,n-NNmax) - J(m,n-NNmax,timepoint) 
 
   !     ENDDO

   !  ENDDO
 
   !  DO m=NNmax+1,2*NNmax

   !     DO n=1,NNmax

   !       !  fvec(m)= fvec(m) + (1.0_real64/dt)*MASS(m-NNmax,n)*c(n,timePoint) -(1.0_real64/dt)*MASS(m-NNmax,n)*c(n,timePoint-1) &
   !       !  + STIFF(m-NNmax,n)*w(n,timepoint) 

   !       !  fjac(m,n)= fjac(m,n) + STIFF(m-NNmax,n)

   !         fvec(m)= fvec(m) + MASS(m-NNmax,n)*c(n,timePoint) -MASS(m-NNmax,n)*c(n,timePoint-1) &
   !         + (dt)*STIFF(m-NNmax,n)*w(n,timepoint) 

   !         fjac(m,n)= fjac(m,n) + (dt)*STIFF(m-NNmax,n)
 
   !     ENDDO
 
 
   !     DO n= NNmax+1,2*NNmax

   !       !  fjac(m,n)= fjac(m,n) + (1.0_real64/dt)*MASS(m-NNmax,n-NNmax)
   !        fjac(m,n)= fjac(m,n) + MASS(m-NNmax,n-NNmax)
 
   !     ENDDO

   !  ENDDO


    !************************************************************
    ! forward euler linear scheme - Provatas/Elder 2010
    !***********************************************************

     DO m=1,NNmax
           fvec(m) = fvec(m)  + f(m,timepoint-1) ! subroutine calculates at t, for t=1,...
       DO n=1,NNmax
         
         
          fjac(m,n) = fjac(m,n) + MASS(m,n)

        ENDDO

         DO n= NNmax+1,2*NNmax

           fjac(m,n)=  fjac(m,n) -(e**2)*STIFF(m,n-NNmax)  
 
        ENDDO

     ENDDO
 
     DO m=NNmax+1,2*NNmax

        DO n=1,NNmax

           fvec(m)= fvec(m) +  MASS(m-NNmax,n)*c(n,timePoint-1) 

           fjac(m,n)= fjac(m,n) + (dt)*STIFF(m-NNmax,n)
 
        ENDDO
 
 
       DO n= NNmax+1,2*NNmax

          fjac(m,n)= fjac(m,n) + MASS(m-NNmax,n-NNmax)
 
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
 
 subroutine inv(A, Ainv, n)
   use, intrinsic :: iso_fortran_env
    IMPLICIT NONE
    INCLUDE 'CHParameters.inc'
   integer, intent(in) :: n
   real(real64), intent(in) :: A(n, n)
   real(real64), intent(out) :: Ainv(n, n)
   
   real(real64) :: temp(n, 2 * n)
   real(real64) :: pivot, factor
   integer :: i, j, k, maxrow
   
   ! Initialize Ainv as an identity matrix
   Ainv = 0.0d0
   do i = 1, n
     Ainv(i, i) = 1.0d0
   end do
 
   ! Augment A with the identity matrix
   temp(:, 1:n) = A
   temp(:, n+1:2*n) = Ainv
 
   ! Perform Gaussian elimination with partial pivoting
   do i = 1, n
     ! Find the pivot row (row with maximum absolute value in current column)
     maxrow = i
     do j = i+1, n
       if (abs(temp(j, i)) > abs(temp(maxrow, i))) maxrow = j
     end do
     
     ! Swap rows i and maxrow
     if (i /= maxrow) then
       do k = 1, 2*n
         pivot = temp(i, k)
         temp(i, k) = temp(maxrow, k)
         temp(maxrow, k) = pivot
       end do
     end if
     
     ! Scale the current row to make the pivot 1
     factor = 1.0d0 / temp(i, i)
     temp(i, :) = temp(i, :) * factor
     
     ! Eliminate other rows
     do j = 1, n
       if (i /= j) then
         factor = -temp(j, i)
         temp(j, :) = temp(j, :) + factor * temp(i, :)
       end if
     end do
   end do
 
   ! Extract the inverse from the augmented matrix
   Ainv = temp(:, n+1:2*n)
   
   return
 end subroutine inv
 
 SUBROUTINE LU_decomposition(A, N)
   use, intrinsic :: iso_fortran_env
   IMPLICIT NONE

   INTEGER, INTENT(IN) :: N
   REAL(real64), INTENT(INOUT) :: A(N, N)

   INTEGER :: i, j, k
   REAL(real64) :: temp
   REAL(real64), DIMENSION(N) :: scale
   INTEGER, DIMENSION(N) :: pivot
   REAL(real64) :: max_pivot

   ! Initialize pivot vector and scaling vector
   DO i = 1, N
       pivot(i) = i
       max_pivot = 0.0_real64
       DO j = 1, N
           max_pivot = MAX(max_pivot, ABS(A(i, j)))
       ENDDO
       IF (max_pivot == 0.0_real64) THEN
           ! Handle the case where max_pivot is zero (singular matrix)
           WRITE(*, *) 'Singular matrix encountered'
           RETURN
       ENDIF
       scale(i) = 1.0_real64 / max_pivot
   ENDDO

   DO k = 1, N - 1
       max_pivot = 0.0_real64
       DO i = k, N
           temp = scale(i) * ABS(A(i, k))
           IF (temp >= max_pivot) THEN
               max_pivot = temp
               pivot(k) = i
           ENDIF
       ENDDO

       IF (A(pivot(k), k) == 0.0_real64) THEN
           ! Handle the case where A(pivot(k), k) is zero (singular matrix)
           WRITE(*, *) 'Singular matrix encountered'
           RETURN
       ENDIF

       IF (pivot(k) /= k) THEN
           ! Swap rows in matrix A
           DO j = 1, N
               temp = A(k, j)
               A(k, j) = A(pivot(k), j)
               A(pivot(k), j) = temp
           ENDDO
           scale(pivot(k)) = scale(k)
       ENDIF

       DO i = k + 1, N
           temp = A(i, k) / A(k, k)
           A(i, k) = temp
           DO j = k + 1, N
               A(i, j) = A(i, j) - temp * A(k, j)
           ENDDO
       ENDDO
   ENDDO

END SUBROUTINE LU_decomposition

SUBROUTINE SolveFJACDSOL( tolS, max_iter, fvec, fjac, dSol)
   use, intrinsic :: iso_fortran_env
   IMPLICIT NONE
   INCLUDE 'CHParameters.inc'

   INTEGER, INTENT(IN) ::  max_iter
   REAL(real64), INTENT(IN) :: tolS
   REAL(real64), INTENT(INOUT) :: fvec(2*NNmax)
   REAL(real64), INTENT(INOUT) :: fjac(2*NNmax, 2*NNmax)
   REAL(real64), INTENT(OUT) :: dSol(2*NNmax)

   INTEGER :: i, j, k
   REAL(real64) :: temp, sum

   ! LU decomposition of fjac
   CALL LU_decomposition(fjac, NNmax)

   ! Back-substitution to solve the system
   DO i = 1, 2*NNmax
       sum = fvec(i)
       DO j = 1, i - 1
           sum = sum - fjac(i, j) * dSol(j)
       ENDDO
       dSol(i) = sum / fjac(i, i)
   ENDDO

   ! Check for convergence
   IF (NORM2(dSol) < tolS) THEN
       WRITE(*, *) 'Converged!'
   ELSE
       WRITE(*, *) 'Did not converge within', max_iter, 'iterations'
   ENDIF

   RETURN
END SUBROUTINE SolveFJACDSOL

! You would need to implement LU_decomposition and NORM2 functions or use existing implementations.

 