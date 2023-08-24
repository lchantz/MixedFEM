   
    !> π = 3.14
    REAL(real64), PARAMETER:: PI=4.D0*DATAN(1.D0)

    !> Mobility coefficient
    REAL(real64), PARAMETER:: mob = 1.0_real64
    !> Mixing Energy Density
    REAL(real64), PARAMETER:: lamda = 1.0_real64
    !> Capillary width
    REAL(real64), PARAMETER:: epsilon = 1.0E-3_real64

    !> Number of elements on x-axis and on y-axis
    INTEGER, PARAMETER:: NELx=5
    INTEGER, PARAMETER:: NELy=5
    INTEGER, PARAMETER:: NELmax=NELx*NELy

    !> Number of nodes on x and on y direction
    INTEGER, PARAMETER:: NNx=2*NELx + 1
    INTEGER, PARAMETER:: NNy=2*NELy + 1
    INTEGER, PARAMETER:: NNmax=NNx*NNy
    INTEGER, PARAMETER:: Nhp=(NNx-1)*NNy !nodes in the horizontal periodicity case


    !> Number of time points
    INTEGER, PARAMETER:: Tmax=10

    !> Number of Nodes on Reference Element
    INTEGER, PARAMETER:: NNref = 9

    !> Time step
    REAL(real64),PARAMETER :: deltat=1.0E-8_real64

 

    !>COMMON/Timeinfo/deltat