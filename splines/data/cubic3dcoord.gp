set title "Spline test"
set key outside below right
set size ratio -1
set colorsequence podo
plot "ctrl.dat" title "Control points of original spline" with linespoints ls 1, \
"cubic3dcoord.points" title "Calculated points on spline", \
"ctrl.dat" using 1:2 smooth bezier title "Spline from control points using smooth bezier" 

pause -1
