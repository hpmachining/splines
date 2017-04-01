set title "Cubic 2d Spline test"
set key outside below right
set size ratio -1
set colorsequence podo
plot "cubic2d.dat" title "Control points of Cubic 2d spline" with linespoints ls 1, \
"cubic2dcoord.points" title "Calculated points on spline", \
"cubic2d.dat" using 1:2 smooth bezier title "Spline from control points using smooth bezier" 

pause -1
