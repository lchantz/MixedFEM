   
    !> π = 3.14
    REAL(real64), PARAMETER:: PI=4.D0*DATAN(1.D0)

    !> Viscosity coefficient
    REAL(real64), PARAMETER:: V=0.01_real64

    !> Number of elements on x-axis
    INTEGER, PARAMETER:: NEL=50


    !> Number of space points
    INTEGER, PARAMETER:: N=2*NEL+1

    !> Number of time points
    INTEGER, PARAMETER:: M=100

    !> Number of Points Per Element
    INTEGER, PARAMETER:: PPE = 3

    !> Time step
    REAL(real64) :: deltat

    COMMON/TimeStep/deltat