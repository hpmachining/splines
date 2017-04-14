set title "Spline test"
set key outside below right
set size ratio -1
set colorsequence podo
plot "quint3d.dat" title "Control points of original spline" with linespoints ls 1, \
"quint3d.dat" using 1:2 smooth bezier title "Spline from control points using smooth bezier" 

pause -1
