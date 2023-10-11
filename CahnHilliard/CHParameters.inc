   
    !> π = 3.14
    REAL(real64), PARAMETER:: PI=4.D0*DATAN(1.D0)
    REAL(real64), PARAMETER:: theta=PI/4.0_real64

    REAL(real64), PARAMETER:: e = 6.25E-2_real64

    !> Number of elements on x-axis and on y-axis
    INTEGER, PARAMETER:: NELx=15
    INTEGER, PARAMETER:: NELy=15
    INTEGER, PARAMETER:: NELmax=NELx*NELy

    !> Number of nodes on x and on y direction
    INTEGER, PARAMETER:: NNx=2*NELx + 1
    INTEGER, PARAMETER:: NNy=2*NELy + 1
    INTEGER, PARAMETER:: NNmax=NNx*NNy
    INTEGER, PARAMETER:: NNp= NNmax - NNy ! # of nodes in the periodic case


    !> Number of time points
    INTEGER, PARAMETER:: Tmax=100

    !> Number of Nodes on Reference Element
    INTEGER, PARAMETER:: NNref = 9

    !> Time step
    REAL(real64),PARAMETER :: dt=5.0E-6_real64

 

    !>COMMON/Timeinfo/dt