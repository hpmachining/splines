set title "3d Quadratic Bezier Curve test"
set key outside below right
set size ratio -1
set colorsequence podo
plot "quad3d.dat" title "Control points of quadratic curve" with linespoints ls 1, \
"quad3dcoord.points" title "Calculated points on quadratic curve", \
"quad3d.dat" using 1:2 smooth bezier title "Quadratic curve from control points using smooth bezier" 

pause -1
