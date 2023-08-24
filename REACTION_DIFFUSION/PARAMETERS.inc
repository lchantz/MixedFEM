   
    !> π = 3.14
    REAL(real64), PARAMETER:: PI=4.D0*DATAN(1.D0)

    REAL(real64), PARAMETER:: e = 1.0_real64

    !> Number of elements on x-axis and on y-axis
    INTEGER, PARAMETER:: NELx=2
    INTEGER, PARAMETER:: NELy=2
    INTEGER, PARAMETER:: NELmax=NELx*NELy

    !> Number of nodes on x and on y direction
    INTEGER, PARAMETER:: NNx=2*NELx + 1
    INTEGER, PARAMETER:: NNy=2*NELy + 1
    INTEGER, PARAMETER:: NNmax=NNx*NNy
    


    !> Number of time points
    INTEGER, PARAMETER:: Tmax=10

    !> Number of Nodes on Reference Element
    INTEGER, PARAMETER:: NNref = 9

    !> Time step Δt
    REAL(real64),PARAMETER :: dt = 1.0e-5_real64

 

    !COMMON/Timeinfo/dt