set title "Cubic 3d Spline test"
set key outside below right
set size ratio -1
set colorsequence podo
plot "cubic3d.dat" title "Control points of Cubic 3d spline" with linespoints ls 1, \
"cubic3dcoord.points" title "Calculated points on spline", \
"cubic3d.dat" using 1:2 smooth bezier title "Spline from control points using smooth bezier" 

pause -1
