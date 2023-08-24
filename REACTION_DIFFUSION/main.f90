!>******************************************************************************
!>    REACTION_DIFFUSION
!>******************************************************************************
program REACTION_DIFFUSION
    use, intrinsic :: iso_fortran_env
    USE BiCGStab_mod
    IMPLICIT NONE
    INCLUDE 'PARAMETERS.inc'
 
    !> Space coordinates
    REAL(real64) x(NNmax),y(NNmax)
 
    !> Time coordinates
    REAL(real64) t(Tmax)
 
    !> Local to global node numbering
    INTEGER nop(NELmax,NNref)
 
 
    !> Solution
    REAL(real64) u(NNmax,Tmax)
 

   
 
 
 
    !> System matrices
    !> STIFF - Stiffness
    !> MASS - Mass
    
    
    REAL(real64) STIFF(NNmax,NNmax)
    
    
    REAL(real64) MASS(NNmax,NNmax)
    
    
    !> P - Non linear term from  potential u(1-u)
    REAL(real64) f(NNmax,Tmax)
  
    !> Derivative of the non-linear term from double well potential
    REAL(real64) J(NNmax,NNmax,Tmax)
  
 
 
    !> Newton tolerance parameters
    REAL(real64),PARAMETER:: tol=1.0E-4_real64
    INTEGER, PARAMETER :: iter_max = 1000
   
    !> Newton solution vector
    REAL(real64) fvec(NNmax)
    REAL(real64) fjac(NNmax,NNmax)
 
 
    ! REAL(real64) S(2*NNmax,ntrial)
    ! REAL(real64) dS(2*NNmax,ntrial)
    REAL(real64) dSol(NNmax)
    REAL(real64) Sol(NNmax)
 
 
    INTEGER element,i,timepoint,iter
    LOGICAL CONVERGENT
 
    OPEN(unit=1,file='InitialCondition.dat')
    OPEN(unit=2,file='Solution.dat')
 
    WRITE(*,'(25x,A)') 'Results u(n,t)'
    CALL SpaceDiscretization(x,y)
    CALL GlobalNodeNumbering(nop)
 
  
    !> Set time grid
    t=0.0_real64
    CALL TimeDiscretization(t)
 
    
    !> Initialization of system matrices and vectors
    STIFF=0.0_real64
    MASS=0.0_real64
  
    !>****************************************************************************
    !> Calculation of Stiffness,Mass,FLUX for the full dimension (Neumann problem)
    !>****************************************************************************
    DO element=1,NELmax
       CALL StiffnessMassMatrices(element,x,y,nop,STIFF,MASS)
    ENDDO
 
   
    !>****************************************************************
    !> Set the initial conditions for the c,w variables of the problem
    !> For timepoint = 1
    !>*****************************************************************
    ! timepoint = 1
    u=0.0_real64
    CALL InitialConditions(x,y,u)
    DO i=1,NNmax
       WRITE(1,'(2x,E16.9,1x,E16.9,1x,E16.9,1x,E16.9)') x(i),y(i),u(i,1)
    ENDDO
    CLOSE(1)
    
    !> Time loop to attain the solution at every timepoint
    DO timepoint=2,Tmax
       print *, "Time point: ", timepoint
       
       iter=0
       CONVERGENT=.false. 
       Sol=0.0_real64

       !> Choose the initial guess for the first newton iteration
       Sol(1:NNmax) = u(:,timepoint-1) 
       
       
       !> Start the Newton iteration loop
       
       DO WHILE (.NOT.CONVERGENT .AND. iter < iter_max)
         dSol=0.0_real64
             !>*****************************************************************
             !> Calculate the double well potential and the derivative
             !> Depend on c, the order parameter
             !>*****************************************************************
             f=0.0_real64
             J=0.0_real64
 
             DO element=1,NELmax
                CALL NonLinearTerm(x,y,element,nop,timepoint,u,f,J)
             ENDDO
 
             !>******************************************************************
             !> Set the Stationary System to be solved at every Newton iteration
             !> ___________________Fjac.dx = -Fvec_______________________________
             !>******************************************************************
             fvec=0.0_real64
             fjac=0.0_real64

             CALL StationarySystem(timePoint,u,STIFF,MASS,f,J,fjac,fvec)
 
             !> Solve the system - find the correction dSol
                fvec = -fvec
                dSol = BiCGStab(fjac,fvec)
 
                !> Update the solution 
                 Sol(:) = Sol(:) + dSol(:)
 
               !> Check for convergence
               CONVERGENT = (NORM2(dSol) < tol) 
               
               iter = iter + 1
               u(:,timepoint) = Sol(1:NNmax)
       ENDDO        
          IF (CONVERGENT) THEN
             write(*,*) 'Newton-Raphson converged in', iter, 'iterations'
          ELSE
             write(*,*) 'Newton-Raphson did not converge within', iter_max, 'iterations'
          ENDIF
 
       
       
       
    ENDDO !> Time loop
 
    !> Write the results
    DO i=1,NNmax
       WRITE(2,'(2x,E16.9,1x,E16.9,1x,E16.9,1x,E16.9)') x(i),y(i),u(i,7)
    ENDDO
 
    CLOSE(2)
  
 END program REACTION_DIFFUSION
 
 
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
    INCLUDE 'PARAMETERS.inc'
 
    !> Space coordinates
    REAL(real64) x(NNmax),y(NNmax)
 
    REAL(real64) xorigin,yorigin,xlast,ylast
    REAL(real64) deltax,deltay
    INTEGER i,j,NNode
 
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
    INCLUDE 'PARAMETERS.inc'
 
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
    INCLUDE 'PARAMETERS.inc'
 
    !> Time coordinates
    REAL(real64) t(Tmax)
 
    REAL(real64) :: tfirst
   !  REAL(real64) :: tlast
    INTEGER j
 
    tfirst=0.0_real64
   !  tlast=1.0_real64
   !  dt = (tlast-tfirst)/real(Tmax-1)
 
    t(1)=tfirst
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
 SUBROUTINE InitialConditions(x,y,u)
    use, intrinsic :: iso_fortran_env
    IMPLICIT NONE
    INCLUDE 'PARAMETERS.inc'
 
    !> Space coordinates
    REAL(real64) x(NNmax),y(NNmax)
 
    !> Solution
    REAL(real64) u(NNmax,Tmax)
    
    
 
    INTEGER i

   
    !> Set to a random distribution for initial time
    DO i=1,NNmax
       u(i,1) =20.0_real64 + (80.0_real64)*(y(i) - sin(pi*x(i)/2.0_real64)*sin(pi*y(i)/2.0_real64) )
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
    INCLUDE 'PARAMETERS.inc'
 
 
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
  !> J:                  Derivative of the non-linear term of the weak formulation
  !>******************************************************************************
 
 
 !>******************************************************************************
 !> NonLinearTerm evaluates the entries of the P(c,t) operator and the J(c,t) matrix
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
 !> J:                  Derivative of the non-linear term
 !>******************************************************************************
 
 
 SUBROUTINE NonLinearTerm(x,y,element,nop,timepoint,u,f,J)
    use, intrinsic :: iso_fortran_env
    IMPLICIT NONE
    INCLUDE 'PARAMETERS.inc'
 
    
    
    !> Space coordinates !xpt,ypt
    REAL(real64), INTENT(IN) :: x(NNmax),y(NNmax)
 
    !> Current element
    INTEGER element
 
    !> Local to global space numbering
    INTEGER nop(NELmax,NNref)
    INTEGER ngl(NNref)
 
    !> timepoint
    INTEGER timepoint
 
    !> Solution
    REAL(real64) u(NNmax,Tmax)
 
 
    !> non linear term
    REAL(real64) f(NNmax,Tmax)
 
    !> Derivative of the non-linear term
    REAL(real64) J(NNmax,NNmax,Tmax)
 
 
 
 
    REAL(real64) phi(NNref), tphx(NNref), tphy(NNref),phic(NNref),phie(NNref)
    REAL(real64) gp(3), gw(3)
    REAL(real64) Xcomputational,Xc,Xe,Ycomputational,Yc,Ye,dett
 
    REAL(real64) ui,ux,uy
    INTEGER i,k,l,m,n,r
 
 
 
    gw  =(/0.27777777777778, 0.444444444444, 0.27777777777778/)
    gp =(/0.1127016654    , 0.5           , 0.8872983346    /)
 
    DO i = 1,NNref
       ngl(i) = nop(element,i)
    ENDDO
 
    DO r = 1,3
       DO k = 1,3
 
          call TestFunctions(gp(r),gp(k),phi,phic,phie)
 
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
 
 
          ui = 0.0_real64
          ux = 0.0_real64
          uy = 0.0_real64
          DO i=1,NNref
             ui = ui + u(ngl(i),timePoint)*phi(i)
             ux = ux + u(ngl(i),timePoint)*tphx(i)
             uy = uy + u(ngl(i),timePoint)*tphy(i)
          ENDDO
 
          DO l=1,NNref
             DO m=1,NNref
                 J(ngl(l),ngl(m),timepoint) =  J(ngl(l),ngl(m),timepoint) + gw(r)*gw(k)*dett*(1.0_real64-2.0_real64*ui)*phi(l)*phi(m)
             ENDDO
             f(ngl(l),timepoint) = f(ngl(l),timepoint) + gw(r)*gw(k)*dett*(ui*(1.0_real64 - ui))*phi(l)
          ENDDO
 
       ENDDO
    ENDDO
 
 
 END SUBROUTINE NonLinearTerm
 
 
 !>******************************************************************************
 !> Constructs the stationary system of each Newton procedure uycle
 !> .
 !> -------------
 !> Params:
 !> -------------
 !> timepoint:           the current time point
 !> c:                   Order parameter
 !> w:                   Chemical Potential
 !> STIFF,MASS:          System matrices
 !> P:                   Non-linear term of the weak formulation
 !> J:                  Derivative of the non-linear term of the weak formulation
 !> fjac:                jacobian matrix of newton procedure
 !> fvec:                system of equations of newton procedure
 !>*********************************************************************************
 SUBROUTINE StationarySystem(timePoint,u,STIFF,MASS,f,J,fjac,fvec)
    use, intrinsic :: iso_fortran_env
    IMPLICIT NONE
    INCLUDE 'PARAMETERS.inc'
 
    !> Current time step
    INTEGER timepoint
 
    !> Solution
    REAL(real64) u(NNmax,Tmax)

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
    REAL(real64) fjac(NNmax,NNmax)
    !> Righ handside vector of Newton system
    REAL(real64) fvec(NNmax)
 
 
    INTEGER m,n
 
    REAL(real64) e2
 
    e2 = (e)**2

    !*******************************************
    !> REACTION DIFFUSION RESIDUAL AND JACOBIAN
    !*******************************************
    DO m=1,NNmax
       DO n=1,NNmax
        fjac(m,n) = fjac(m,n) + (e2)*(STIFF(m,n)) - (1.0_real64/dt)*MASS(m,n) + J(m,n,timepoint)
        fvec(m) = fvec(m) + ((e2)*(STIFF(m,n)) - (1.0_real64/dt)*MASS(m,n))*u(n,timepoint) + (1.0_real64/dt)*MASS(m,n)*u(n,timepoint-1)
       ENDDO
        fvec(m) = fvec(m) + f(m,timepoint)
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
    INCLUDE 'PARAMETERS.inc'
 
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
 
 
 