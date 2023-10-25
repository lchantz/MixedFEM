module BiCGStab_mod
    use, intrinsic :: iso_fortran_env
    implicit none

    integer (kind=2 ), parameter    :: sp = kind(1.000)
    integer (kind=2 ), parameter    :: doublep = kind(1.0d0)
    integer (kind=2 ), parameter    :: qp = kind(1.0q0)
    integer (kind=2 ), parameter    :: wp = doublep         ! can be changed to "qp" to get quad precision.

contains

    function BiCGStab(A,b) result(x)
        implicit none

!--------------------------PARAMETER AND VARIABLE-------------------------------!
        real    (real64), intent(in )                   :: A (:,:)
        real    (real64), intent(in )                   :: b ( : )
        real    (real64), dimension(1:size(b, dim=1))   :: x

        real    (real64), dimension(1:size(b, dim=1))   :: r, rs, v, p, s, t
        real    (real64), parameter                     :: e = 1.0E-4_real64
        real    (real64)                                :: rho      , rho_prev
        real    (real64)                                :: alpha    , omega   , beta
        real    (real64)                                :: norm_r   , norm_b       
        real    (real64)                                :: summesion, temp

        integer (kind=sp)                                :: it=0,err
!------------------------END PARAMETER AND VARIABLE-----------------------------!  

        if(size(A, dim=1) /= size(A, dim=2)) stop &
        "Error: Improper dimension of matrix A in BiCGStab."





!-------------------------------------------------------!
        x  = 0.0_real64                                     !-------> INITIAL GUESS
!-------------------------------------------------------!
        r  = b - matmul(A,x)                            !-------> LINE 1
        rs = r                                          !
!-------------------------------------------------------!
        rho   = 1.0_real64; alpha = 1.0_real64; omega = 1.0_real64  !-------> LINE 2
!-------------------------------------------------------!
        v  = 0.0_real64; p  = 0.0_real64                        !-------> LINE 3
!                                                       !
        norm_r = sqrt(dot_product(r,r))                 !
        norm_b = sqrt(dot_product(b,b))                 !
!-------------------------------------------------------!




        it = 0
        do while(norm_r .GT. e*norm_b)                          !-------> START OF LOOP
            
        !-------------------------------------------------------!
            rho_prev = rho                                   !-------> LINE 5
            rho      = dot_product(rs,r)                        !
        !-------------------------------------------------------!
            if (rho_prev == 0.0_real64) then
                beta = 0.0_real64
            else
                beta = (rho/rho_prev) * (alpha/omega)
            end if
            
            ! beta     = (rho/rho_prev) * (alpha/omega)           !-------> LINE 6
        !-------------------------------------------------------!
            p        = r + beta * (p - omega*v)                 !-------> LINE 7
        !-------------------------------------------------------!
            v        = matmul(A,p)                              !-------> LINE 8
        !-------------------------------------------------------!
            alpha    = rho/dot_product(rs,v)                    !-------> LINE 9
        !-------------------------------------------------------!
            s        = r - alpha*v                              !-------> LINE 10
        !-------------------------------------------------------!
            t        = matmul(A,s)                              !-------> LINE 11
        !-------------------------------------------------------!
            omega    = dot_product(t,s)/dot_product(t,t)        !-------> LINE 12
        !-------------------------------------------------------!
            x        = x + alpha*p + omega*s                    !-------> LINE 13
        !-------------------------------------------------------!
            r        = s - omega*t                              !-------> LINE 17
        !-------------------------------------------------------!
            norm_r   = sqrt(dot_product(r,r))                   !
            norm_b   = sqrt(dot_product(b,b))                   !
        !-------------------------------------------------------!
            it = it + 1                                         !
        !-------------------------------------------------------!

        end do                                                      !-------> END OF LOOP

        print*, "Iteration Required :", it

return
end function BiCGStab 


end module BiCGStab_mod