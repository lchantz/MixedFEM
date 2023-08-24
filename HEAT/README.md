# RUNNING using Makefile
```shell
make
```

# PLOTTING 
```shell
gnuplot
set yrange [0:1]
set xrange [0:2]
set dgrid 100,100 qnorm2
splot "BurgersSolution.dat" with lines
```