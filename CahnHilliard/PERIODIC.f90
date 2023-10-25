!>******************************************************************************
!>    CAHN-HILLIARD WITH PERIODIC BOUNDARY CONDITION
!>******************************************************************************
program PERIODIC
    use, intrinsic :: iso_fortran_env
    USE BiCGStab_mod
    IMPLICIT NONE
    INCLUDE 'PARAMETERS_PERIODIC.inc'
 
    !> Space coordinates
    REAL(real64) x(NNmax),y(NNmax)
 
    !> Time coordinates
    REAL(real64) t(Tmax)
 
    !> Local to global node numbering
    INTEGER nop(NELmax,NNref)
 
    !> Order parameter - u - of the Cahn Hilliard Equation
    REAL(real64) u_new(NNmax)
    REAL(real64) u_old(NNmax)

 
    !> Auxiliary variable - w - the chemical potential
    REAL(real64) w_new(NNmax)
    REAL(real64) w_old(NNmax)
   
 
    !> Newton tolerance parameters
    REAL(real64),PARAMETER:: tolS=1.0E-5_real64
    INTEGER, PARAMETER :: ntrial = 100
    INTEGER, PARAMETER :: max_iter = 1000
 
    !> Stationary system of Newton cycle
    REAL(real64) fvec(2*NNmax)
    REAL(real64) fjac(2*NNmax,2*NNmax)
    !> Local contributions
    REAL(real64) fjac1(NNmax,NNmax)
    REAL(real64) fjac2(NNmax,NNmax)
    REAL(real64) fjac3(NNmax,NNmax)
    REAL(real64) fjac4(NNmax,NNmax)
    
    REAL(real64) fvec1(NNmax)
    REAL(real64) fvec2(NNmax)
 
    !> Newton solution vector
    REAL(real64) dSol(2*NNmax)
    REAL(real64) Sol(2*NNmax)

    INTEGER element,timepoint,trial,n,m,l
    LOGICAL CONVERGENT
    CHARACTER(10) :: timepoint_str
 
    OPEN(unit=1,file='initial_PERIODIC.dat')

    
 
    WRITE(*,'(25x,A)') 'Results c(n,t)'
    CALL SpaceDiscretization(x,y)
    CALL GlobalNodeNumbering(nop)
 
 
    !> Set time grid
    t=0.0_real64
    CALL TimeDiscretization(t)
 
    !>****************************************************************
    !> Set the initial conditions for the c,w variables of the problem
    !> For timepoint = 1
    !>*****************************************************************
    u_old=0.0_real64
    u_new=0.0_real64
    w_old=0.0_real64
    w_new=0.0_real64
    CALL InitialConditions(x,y,u_old,u_new,w_old,w_new)

    DO n=1,NNx
       WRITE(1,*) (u_old((n-1)*NNy+m),m=1,NNy)
    ENDDO
    CLOSE(1)

    
     !>***********************************************************************************
     !> Given an initial guess x for a root in n dimensions, take ntrial Newton-Raphson
     !> steps to improve the root. Stop if the root converges in either summed absolute
     !> variable increments tolx or summed absolute function values tolf.
     !>***********************************************************************************

 
   !> Time loop to attain the solution at every timepoint
    DO timepoint=2,Tmax
      print *, "Time point: ", timepoint

      trial=0
      CONVERGENT=.false.
      Sol=0.0_real64


      
      !> Choose the initial guess for the first newton iteration
      Sol(1:NNmax) = w_old(:) 
      Sol(NNmax+1:2*NNmax) = u_old(:) 
      
      !> Start the Newton iteration loop
      DO WHILE (.NOT.CONVERGENT .AND. trial < ntrial)
         dSol=0.0_real64

         !>******************************************************************
         !> Set the Stationary System to be solved at every Newton iteration
         !> ___________________Fjac.dx = -Fvec_______________________________
         !>******************************************************************
          fvec=0.0_real64
          fjac=0.0_real64

          fjac1=0.0_real64
          fjac2=0.0_real64
          fjac3=0.0_real64
          fjac4=0.0_real64

          fvec1=0.0_real64
          fvec2=0.0_real64

          !> Local contributions
          DO element=1,NELmax
          CALL AssembleSystem(x,y,element,nop,u_old,u_new,w_new,fjac1,fjac2,fjac3,fjac4,fvec1,fvec2)
          ENDDO
          !> Global system 
          fvec(1:NNmax) = fvec1(:)
          fvec(NNmax+1:2*NNmax) = fvec2(:)

         fjac(1:NNmax,1:NNmax) = fjac1(:,:)
         fjac(1:NNmax,NNmax+1:2*NNmax) = fjac2(:,:)
         fjac(NNmax+1:2*NNmax,1:NNmax) = fjac3(:,:)
         fjac(NNmax+1:2*NNmax,NNmax+1:2*NNmax) = fjac4(:,:)

         !> Imposing the periodic conditions
         DO l=1,NNy 
            fvec(l) = fvec(l) + fvec(NNmax-NNy+l)
            fvec(NNmax-NNy+l) = 0.0_real64 

            fjac(l,:) = fjac(l,:) + fjac(NNmax-NNy+l,:)
            fjac(:,l) = fjac(:,l) + fjac(:,NNmax-NNy+l)
            fjac(NNmax-NNy+l,:) = 0.0_real64
            fjac(:,NNmax-NNy+l) = 0.0_real64
            fjac(NNmax-NNy+l,NNmax-NNy+l) = 1.0_real64
         ENDDO        
         DO l=2,NNx-1 
            fvec((l-1)*NNy+1) = fvec((l-1)*NNy+1) + fvec(l*NNy)
            fvec(l*NNy) = 0.0_real64 

            fjac((l-1)*NNy+1,:) = fjac((l-1)*NNy+1,:) + fjac(l*NNy,:)
            fjac(:,(l-1)*NNy+1) = fjac(:,(l-1)*NNy+1) + fjac(:,l*NNy)
            fjac(l*NNy,:) = 0.0_real64
            fjac(:,l*NNy) = 0.0_real64
            fjac(l*NNy,l*NNy) = 1.0_real64
         ENDDO        

         !> Solve the system - find the correction dSol
         fvec = -fvec
         dSol = BiCGStab(fjac,fvec)
         
         !> Recover the right boundary nodes
         DO l=1,NNy 
            dSol(NNmax-NNy+l) = dSol(l)
            dsol(2*NNmax-NNy+l) = dSol(NNmax+l)
         ENDDO 
         !> Recover the top nodes  
         DO l=2,NNx-1 
            dSol(l*NNy) = dSol((l-1)*NNy+1)
            dsol(NNmax+l*NNy) = dSol(NNmax+(l-1)*NNy+1)
         ENDDO 

         !> Update the solution
         Sol(:) = Sol(:) + dSol(:)

         DO l=1,2*NNmax
            IF (Sol(l) <= - 0.9999) THEN 
               Sol(l) = - 0.9999 
            ELSEIF (Sol(l) >= 0.9999) THEN 
               Sol(l) = 0.9999
            ENDIF   
         ENDDO   
     
         !> Check for convergence
         CONVERGENT = (NORM2(dSol) < tolS)
         print*,"CONVERGENCE", CONVERGENT

         trial = trial + 1
         w_new(:) = Sol(1:NNmax)
         u_new(:) = Sol(NNmax+1:2*NNmax)

      ENDDO  !> Newton trial loop
      
      IF (CONVERGENT) THEN

               u_old(:) = u_new(:)
               w_old(:) = w_new(:)

         write(*,*) 'Newton-Raphson converged in', trial, 'iterations'


         IF ((timepoint<=100 .and. MOD(timepoint,10)==0).or.MOD(timepoint, 500) == 0) THEN
            
            WRITE(timepoint_str, '(I0)') timepoint  ! Convert the integer to a string
          
            OPEN(UNIT=2, FILE=TRIM(ADJUSTL(timepoint_str)) // '_PERIODIC.dat')
            DO n = 1, NNx
              WRITE(2, *) (u_new((n - 1) * NNy + m), m = 1, NNy)
            ENDDO
          
            CLOSE(2)
          ENDIF

      

      ELSE
         write(*,*) 'Newton-Raphson did not converge within', ntrial, 'iterations'
      ENDIF


   ENDDO !> Time loop

 END program PERIODIC
 
 
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
    INCLUDE 'PARAMETERS_PERIODIC.inc'
 
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
    INCLUDE 'PARAMETERS_PERIODIC.inc'
 
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
    INCLUDE 'PARAMETERS_PERIODIC.inc'
 
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
   !  DO i=1,(NNx-1)*NNy+1,NNy !bottom side nodes
   !     BottomNode(i)=1
   !     BoundaryNode(i)=1
   !  ENDDO
   !  DO i=(NNx-1)*NNy+1,NNmax,1 !right side nodes
   !     RightNode(i)=1
   !     BoundaryNode(i)=1
   !  ENDDO
   !  DO i=NNmax,NNy,-NNy !top side nodes
   !     TopNode(i)=1
   !     BoundaryNode(i)=1
   !  ENDDO
   !  DO i=NNy,1,-1 !left side nodes
   !     LeftNode(i)=1
   !     BoundaryNode(i)=1
   !  ENDDO
 
 END SUBROUTINE
 
 SUBROUTINE BoundaryElements(BottomElement,RightElement,TopElement,LeftElement,BoundaryElement)
    use, intrinsic :: iso_fortran_env
    IMPLICIT NONE
    INCLUDE 'PARAMETERS_PERIODIC.inc'
 
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
    INCLUDE 'PARAMETERS_PERIODIC.inc'
 
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
 SUBROUTINE InitialConditions(x, y, u_old,u_new,w_old,w_new)
   use, intrinsic :: iso_fortran_env
   IMPLICIT NONE
   INCLUDE 'PARAMETERS_PERIODIC.inc'

   !> Space coordinates
   REAL(real64) x(NNmax),y(NNmax)

   !> Order parameter - c - of the Cahn Hilliard Equation
   REAL(real64), DIMENSION(NNmax) :: u_old
   REAL(real64), DIMENSION(NNmax) :: u_new

   !> Auxiliary variable - w - the chemical potential
   REAL(real64), DIMENSION(NNmax) :: w_old
   REAL(real64), DIMENSION(NNmax) :: w_new

   INTEGER :: i,n,m



   DO i = 1, NNmax
     CALL RANDOM_NUMBER(u_old(i))
     u_old(i) = 0.4_real64 + 0.02_real64*(0.5_real64-u_old(i))
     w_old(i) = 0.0_real64
   ENDDO

 

END SUBROUTINE InitialConditions


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
    INCLUDE 'PARAMETERS_PERIODIC.inc'
 
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
 
 
 SUBROUTINE AssembleSystem(x,y,element,nop,u_old,u_new,w_new,fjac1,fjac2,fjac3,fjac4,fvec1,fvec2)
   use, intrinsic :: iso_fortran_env
   IMPLICIT NONE
   INCLUDE 'PARAMETERS_PERIODIC.inc'



   !> Space coordinates !xpt,ypt
   REAL(real64), INTENT(IN) :: x(NNmax),y(NNmax)

   !> Current element
   INTEGER element

   !> Local to global space numbering
   INTEGER nop(NELmax,NNref)
   INTEGER ngl(NNref)

   !> Order parameter - c - of the Cahn Hilliard Equation
   REAL(real64) u_new(NNmax)
   REAL(real64) u_old(NNmax)
   REAL(real64) w_new(NNmax)
 

 
   !> Local contributions
   REAL(real64) fjac1(NNmax,NNmax)
   REAL(real64) fjac2(NNmax,NNmax)
   REAL(real64) fjac3(NNmax,NNmax)
   REAL(real64) fjac4(NNmax,NNmax)
   
   REAL(real64) fvec1(NNmax)
   REAL(real64) fvec2(NNmax)



   REAL(real64) phi(NNref), tphx(NNref), tphy(NNref),phic(NNref),phie(NNref)
   REAL(real64) gp(3), gw(3)
   REAL(real64) Xdomain,Xc,Xe,Ydomain,Yc,Ye,dett

   REAL(real64) u0i,ui,ux,uy,wi,wx,wy
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

         u0i = 0.0_real64
         ui = 0.0_real64
         ux = 0.0_real64
         uy = 0.0_real64
         wi = 0.0_real64
         wx = 0.0_real64
         wy = 0.0_real64
         DO i=1,NNref
            u0i = u0i + u_old(ngl(i))*phi(i)
            ui = ui + u_new(ngl(i))*phi(i) 
            ux = ux + u_new(ngl(i))*tphx(i)
            uy = uy + u_new(ngl(i))*tphy(i)
         ENDDO
         DO i=1,NNref
            wi = wi + w_new(ngl(i))*phi(i) 
            wx = wx + w_new(ngl(i))*tphx(i)
            wy = wy + w_new(ngl(i))*tphy(i)
         ENDDO
         
        !*********************************
        !W = 100c**2(1-c)**2 FENICS
        !*********************************  
         ! DO l=1,NNref
         !   f(ngl(l)) = f(ngl(l)) + gw(r)*gw(k)*dett*(100.0)*(ui**2)*(( 1.0 - ui)**2)*phi(l)
         !    DO m=1,NNref
         !       J(ngl(l),ngl(m)) =  J(ngl(l),ngl(m)) + gw(r)*gw(k)*dett*(200.0)*(ui*(1.0-ui)*(1.0-2.0*ui))*phi(l)*phi(m)
         !    ENDDO
         ! ENDDO

        !***************************************************
        !W = 0.25(c**2-1)**2 ginzburg landau 
        !*************************************************** 
        !  DO l=1,NNref
        !    f(ngl(l)) = f(ngl(l)) + gw(r)*gw(k)*dett*(ui**3 - ui)*phi(l)
        !    DO m=1,NNref
        !       J(ngl(l),ngl(m)) =  J(ngl(l),ngl(m)) + gw(r)*gw(k)*dett*((3.0_real64)*(ui**2) - 1.0_real64)*phi(l)*phi(m)
        !    ENDDO
        ! ENDDO
        
        !****************************************************
        !W = 0.25(c**2-1)**2 ginzburg landau & convex splitting
        !**************************************************** 
        ! DO l=1,NNref
        !    f(ngl(l)) = f(ngl(l)) + gw(r)*gw(k)*dett*(ui**3)*phi(l)
        !    DO m=1,NNref
        !       J(ngl(l),ngl(m)) =  J(ngl(l),ngl(m)) + gw(r)*gw(k)*dett*((3.0_real64)*(ui**2))*phi(l)*phi(m)
        !    ENDDO
        ! ENDDO

        
         DO l=1,NNref
            
            fvec1(ngl(l)) = fvec1(ngl(l)) + gw(r)*gw(k)*dett*wi*phi(l) -gw(r)*gw(k)*dett*(2.0_real64*ui - 6.0_real64*ui**2 + 4.0_real64*ui**3 )*phi(l) &
             -(e)**2*gw(r)*gw(k)*dett*(ux*tphx(l)+uy*tphy(l))
            
               fvec2(ngl(l)) = fvec2(ngl(l)) + gw(r)*gw(k)*dett*ui*phi(l) - gw(r)*gw(k)*dett*u0i*phi(l) &
               + 1.0*(dt)*gw(r)*gw(k)*dett*(wx*tphx(l)+wy*tphy(l)) 
               
            DO m=1,NNref
               
                  fjac1(ngl(l),ngl(m)) = fjac1(ngl(l),ngl(m)) + gw(r)*gw(k)*dett*phi(l)*phi(m)
               
                  fjac2(ngl(l),ngl(m)) = fjac2(ngl(l),ngl(m)) -(e)**2*gw(r)*gw(k)*dett*(tphx(l)*tphx(m)+tphy(l)*tphy(m)) &
                  - gw(r)*gw(k)*dett*((2.0_real64 - 12.0_real64*ui+ 12.0_real64*ui**2))*phi(l)*phi(m)
               
                  fjac3(ngl(l),ngl(m)) = fjac3(ngl(l),ngl(m)) + (dt)*gw(r)*gw(k)*dett*(tphx(l)*tphx(m)+tphy(l)*tphy(m))
               
                  fjac4(ngl(l),ngl(m)) = fjac4(ngl(l),ngl(m)) + gw(r)*gw(k)*dett*phi(l)*phi(m)
                  
            ENDDO
         ENDDO   
      ENDDO
   ENDDO


END SUBROUTINE AssembleSystem