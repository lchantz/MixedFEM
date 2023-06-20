
module BiCGStab_mod
   use, intrinsic :: iso_fortran_env
   implicit none

contains

   function BiCGStab(A,b) result(x)
      implicit none

!--------------------------PARAMETER AND VARIABLE-------------------------------!
      real    (real64), intent(in )                   :: A (:,:)
      real    (real64), intent(in )                   :: b ( : )
      real    (real64), dimension(1:size(b, dim=1))   :: x

      real    (real64), dimension(1:size(b, dim=1))   :: r, rs, v, p, s, t
      real    (real64), parameter                     :: e = 1d-33,small = 1.e-10_real64
      real    (real64)                                :: rho      , rho_prev
      real    (real64)                                :: alpha    , omega   , beta
      real    (real64)                                :: norm_r   , norm_b


      integer                                 :: it=0
!------------------------END PARAMETER AND VARIABLE-----------------------------!

      if(size(A, dim=1) /= size(A, dim=2)) stop &
         "Error: Improper dimension of matrix A in BiCGStab."






!-------------------------------------------------------!
      x  = 0.0                                !-------> INITIAL GUESS
!-------------------------------------------------------!
      r  = b - matmul(A,x)                    !-------> LINE 1
      rs = r                                          !
!-------------------------------------------------------!
      rho   = 1.0; alpha = 1.0; omega = 1.0  !-------> LINE 2
!-------------------------------------------------------!
      v  = 0.0; p  = 0.0                     !-------> LINE 3
!                                                       !
      norm_r = sqrt(dot_product(r,r))                 !
      norm_b = sqrt(dot_product(b,b))                 !
!-------------------------------------------------------!





      do while(norm_r .GT. e*norm_b .and. it < 200000)                          !-------> START OF LOOP
         !    print *, norm_r, e*norm_b, (norm_r .GT. e*norm_b)
         !-------------------------------------------------------!
         rho_prev = rho                                      !-------> LINE 5
         rho      = dot_product(rs,r)                        !
         !-------------------------------------------------------!
         beta     = dble((rho/MAX(small,rho_prev))) * dble((alpha/omega))           !-------> LINE 6
         !-------------------------------------------------------!
         p        = r + dble(beta) * dble((p - dble(omega)*dble(v)))                 !-------> LINE 7
         !-------------------------------------------------------!
         v        = matmul(A,p)                              !-------> LINE 8
         !-------------------------------------------------------!
         alpha    = dble(rho) / dble(dot_product(rs,v))                    !-------> LINE 9
         !-------------------------------------------------------!
         s        = r - dble(alpha) * dble(v)                              !-------> LINE 10
         !-------------------------------------------------------!
         t        = matmul(A,s)                              !-------> LINE 11
         !-------------------------------------------------------!
         omega    = dble(dot_product(t,s))/dble(dot_product(t,t))        !-------> LINE 12
         !-------------------------------------------------------!
         x        = x + dble(alpha)*dble(p) + dble(omega)*dble(s)                    !-------> LINE 13
         !-------------------------------------------------------!
         r        = s - dble(omega)*dble(t)                              !-------> LINE 17
         !-------------------------------------------------------!
         norm_r   = dble(sqrt(dot_product(r,r)))                   !
         norm_b   = dble(sqrt(dot_product(b,b)))                   !
         !-------------------------------------------------------!
         it = it + 1                                         !
         !-------------------------------------------------------!

      end do                                                      !-------> END OF LOOP

      print*, "Iteration Required :", it

      return
   end function BiCGStab

end module BiCGStab_mod


! module BiCGStab_mod
!    use, intrinsic :: iso_fortran_env
!    implicit none

! contains

!    subroutine apply_preconditioner(A, M, b)
!       implicit none

!       external dgesv
!       real (real64), intent(in) :: A(:,:)
!       real (real64), intent(out) :: M(:,:)
!       real (real64), intent(in) :: b(:)

!       integer :: n, nrhs, lda, ldb, info
!       integer, dimension(:), allocatable :: ipiv

!       n = size(A, 1)
!       nrhs = size(b, 1)
!       lda = n
!       ldb = n

!       allocate(ipiv(n))

!       M = A ! copy A into M
!       call dgesv(n, nrhs, M, lda, ipiv, b, ldb, info)

!       deallocate(ipiv)

!    end subroutine apply_preconditioner

!    function BiCGStab(A, b) result(x)
!       implicit none

!       real, parameter :: small = 1.E-10
!       real    (real64), intent(in)                   :: A (:,:)
!       real    (real64), intent(in)                   :: b ( : )
!       real    (real64), dimension(size(b))           :: x
!       real    (real64), dimension(size(b))           :: r, rs, v, p, s, t
!       real    (real64), parameter                     :: e = 1d-33
!       real    (real64)                                :: rho      , rho_prev
!       real    (real64)                                :: alpha    , omega   , beta
!       real    (real64)                                :: norm_r   , norm_b
!       integer                                         :: it = 0
!       real    (real64), dimension(size(b), size(b))   :: M

!       if(size(A, dim=1) /= size(A, dim=2)) stop &
!          "Error: Improper dimension of matrix A in BiCGStab."

!       ! Apply preconditioner to matrix A
!       call apply_preconditioner(A, M, b)

!       x  = 0.0
!       ! Compute initial residual
!       r  = b - matmul(A,x)
!       rs = r
!       norm_r = sqrt(dot_product(r,r))
!       norm_b = sqrt(dot_product(b,b))

!       ! Initialize parameters
!       rho   = 1.0
!       alpha = 1.0
!       omega = 1.0
!       v     = 0.0
!       p     = 0.0


!       do while(norm_r .GT. e*norm_b .and. it < 100000)                          !-------> START OF LOOP
!          !    print *, norm_r, e*norm_b, (norm_r .GT. e*norm_b)
!          !-------------------------------------------------------!
!          rho_prev = rho                                      !-------> LINE
!          rho      = dot_product(rs,r)                        !
!          !-------------------------------------------------------!
!          beta     = (rho/MAX(small,rho_prev)) * (alpha/omega)           !-------> LINE 6
!          !-------------------------------------------------------!
!          p        = r + beta * (p - omega*v)                 !-------> LINE 7
!          !-------------------------------------------------------!
!          v        = matmul(A,p)                              !-------> LINE 8
!          !-------------------------------------------------------!
!          alpha    = rho/dot_product(rs,v)                    !-------> LINE 9
!          !-------------------------------------------------------!
!          s        = r - alpha*v                              !-------> LINE 10
!          !-------------------------------------------------------!
!          t        = matmul(A,s)                              !-------> LINE 11
!          !-------------------------------------------------------!
!          omega    = dot_product(t,s)/dot_product(t,t)        !-------> LINE 12
!          !-------------------------------------------------------!
!          x        = x + alpha*p + omega*s                    !-------> LINE 13
!          !-------------------------------------------------------!
!          r        = s - omega*t                              !-------> LINE 17
!          !-------------------------------------------------------!
!          norm_r   = sqrt(dot_product(r,r))                   !
!          norm_b   = sqrt(dot_product(b,b))                   !
!          !-------------------------------------------------------!
!          it = it + 1                                         !
!          !-------------------------------------------------------!

!       end do                                                      !-------> END OF LOOP

!       print*, "Iteration Required :", it

!       return
!    end function BiCGStab

! end module BiCGStab_mod
