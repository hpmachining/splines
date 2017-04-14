set title "B-Spline test"
set key outside below right
#set size ratio -1
set colorsequence podo
plot "ctrl.dat" title "Control points of original spline" with linespoints ls 1, \
"spline.dat" title "Polyline from points on spline" with lines, \
"ctrl.dat" using 1:2 smooth bezier title "Spline from control points using smooth bezier", \
"b-spline.dat" using 1:2 smooth bezier title "Spline from b-spline control points using smooth bezier" 

pause -1
