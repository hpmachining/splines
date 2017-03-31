set title "Parabola test"
set key outside below right
set size ratio -1
set colorsequence podo
plot "parabola.dat" title "Control points of original parabola" with linespoints ls 1, \
"parabola.points" title "Calculated points on parabola", \
"parabola.dat" using 1:2 smooth bezier title "Parabola from control points using smooth bezier" 

pause -1
