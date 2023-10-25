set view map
set size square
set xrange [*:*] noextend
set yrange [*:*] noextend
set key title 'con' at 75, 68
set palette defined ( 0 'blue', 0.5 'grey', 1 'red' )
set pm3d map interpolate 9,9
splot 'initial_NEUMANN.dat' matrix with pm3d notitle

set view map
set size square
set xrange [*:*] noextend
set yrange [*:*] noextend
set key title 'con' at 75, 68
set palette defined ( 0 'blue', 0.5 'grey', 1 'red' )
set pm3d map interpolate 9,9
splot 'initial_PERIODIC.dat' matrix with pm3d notitle

set view map
set size square
set xrange [*:*] noextend
set yrange [*:*] noextend
set key title 'con' at 75, 68
set palette defined ( 0 'blue', 0.5 'grey', 1 'red' )
set pm3d map interpolate 9,9
splot 'initial_CONTACT.dat' matrix with pm3d notitle



