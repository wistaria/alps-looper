set xlabel 'T'

set ylabel 'Energy Density'
plot \
  "< choose REPRESENTATION diagonalization kagome-1-plot-energy.dat" u 1:2 w lp, \
  "< choose REPRESENTATION 'path integral' kagome-1-plot-energy.dat" u 1:2:3 w e, \
  "< choose REPRESENTATION SSE kagome-1-plot-energy.dat" u 1:2:3 w e
pause -1

set ylabel 'Specific Heat'
plot \
  "< choose REPRESENTATION diagonalization kagome-1-plot-sheat.dat" u 1:2 w lp, \
  "< choose REPRESENTATION 'path integral' kagome-1-plot-sheat.dat" u 1:2:3 w e, \
  "< choose REPRESENTATION SSE kagome-1-plot-sheat.dat" u 1:2:3 w e
pause -1

set ylabel 'Magnetization Density^2'
plot \
  "< choose REPRESENTATION diagonalization kagome-1-plot-m2.dat" u 1:2 w lp, \
  "< choose REPRESENTATION 'path integral' kagome-1-plot-m2.dat" u 1:2:3 w e, \
  "< choose REPRESENTATION SSE kagome-1-plot-m2.dat" u 1:2:3 w e
pause -1

set ylabel 'Susceptibility'
plot \
  "< choose REPRESENTATION diagonalization kagome-1-plot-sus.dat" u 1:2 w lp, \
  "< choose REPRESENTATION 'path integral' kagome-1-plot-sus.dat" u 1:2:3 w e, \
  "< choose REPRESENTATION SSE kagome-1-plot-sus.dat" u 1:2:3 w e
pause -1

set ylabel 'Energy Density'
plot \
  "< choose REPRESENTATION diagonalization kagome-2-plot-energy.dat" u 1:2 w lp, \
  "< choose REPRESENTATION 'path integral' kagome-2-plot-energy.dat" u 1:2:3 w e, \
  "< choose REPRESENTATION SSE kagome-2-plot-energy.dat" u 1:2:3 w e
pause -1

set ylabel 'Specific Heat'
plot \
  "< choose REPRESENTATION diagonalization kagome-2-plot-sheat.dat" u 1:2 w lp, \
  "< choose REPRESENTATION 'path integral' kagome-2-plot-sheat.dat" u 1:2:3 w e, \
  "< choose REPRESENTATION SSE kagome-2-plot-sheat.dat" u 1:2:3 w e
pause -1

set ylabel 'Magnetization Density^2'
plot \
  "< choose REPRESENTATION diagonalization kagome-2-plot-m2.dat" u 1:2 w lp, \
  "< choose REPRESENTATION 'path integral' kagome-2-plot-m2.dat" u 1:2:3 w e, \
  "< choose REPRESENTATION SSE kagome-2-plot-m2.dat" u 1:2:3 w e
pause -1

set ylabel 'Susceptibility'
plot \
  "< choose REPRESENTATION diagonalization kagome-2-plot-sus.dat" u 1:2 w lp, \
  "< choose REPRESENTATION 'path integral' kagome-2-plot-sus.dat" u 1:2:3 w e, \
  "< choose REPRESENTATION SSE kagome-2-plot-sus.dat" u 1:2:3 w e
pause -1

set ylabel 'Energy Density'
plot \
  "< choose REPRESENTATION diagonalization kagome-3-plot-energy.dat" u 1:2 w lp, \
  "< choose REPRESENTATION 'path integral' kagome-3-plot-energy.dat" u 1:2:3 w e, \
  "< choose REPRESENTATION SSE kagome-3-plot-energy.dat" u 1:2:3 w e
pause -1

set ylabel 'Specific Heat'
plot \
  "< choose REPRESENTATION diagonalization kagome-3-plot-sheat.dat" u 1:2 w lp, \
  "< choose REPRESENTATION 'path integral' kagome-3-plot-sheat.dat" u 1:2:3 w e, \
  "< choose REPRESENTATION SSE kagome-3-plot-sheat.dat" u 1:2:3 w e
pause -1

set ylabel 'Magnetization Density^2'
plot \
  "< choose REPRESENTATION diagonalization kagome-3-plot-m2.dat" u 1:2 w lp, \
  "< choose REPRESENTATION 'path integral' kagome-3-plot-m2.dat" u 1:2:3 w e, \
  "< choose REPRESENTATION SSE kagome-3-plot-m2.dat" u 1:2:3 w e
pause -1

set ylabel 'Susceptibility'
plot \
  "< choose REPRESENTATION diagonalization kagome-3-plot-sus.dat" u 1:2 w lp, \
  "< choose REPRESENTATION 'path integral' kagome-3-plot-sus.dat" u 1:2:3 w e, \
  "< choose REPRESENTATION SSE kagome-3-plot-sus.dat" u 1:2:3 w e
pause -1

set ylabel 'Sign'
plot \
  "< choose REPRESENTATION 'path integral' kagome-3-plot-sign.dat" u 1:2:3 w e, \
  "< choose REPRESENTATION SSE kagome-3-plot-sign.dat" u 1:2:3 w e
pause -1

set ylabel 'Energy Density'
plot \
  "< choose REPRESENTATION diagonalization kagome-4-plot-energy.dat" u 1:2 w lp, \
  "< choose REPRESENTATION 'path integral' kagome-4-plot-energy.dat" u 1:2:3 w e, \
  "< choose REPRESENTATION SSE kagome-4-plot-energy.dat" u 1:2:3 w e
pause -1

set ylabel 'Specific Heat'
plot \
  "< choose REPRESENTATION diagonalization kagome-4-plot-sheat.dat" u 1:2 w lp, \
  "< choose REPRESENTATION 'path integral' kagome-4-plot-sheat.dat" u 1:2:3 w e, \
  "< choose REPRESENTATION SSE kagome-4-plot-sheat.dat" u 1:2:3 w e
pause -1

set ylabel 'Magnetization Density^2'
plot \
  "< choose REPRESENTATION diagonalization kagome-4-plot-m2.dat" u 1:2 w lp, \
  "< choose REPRESENTATION 'path integral' kagome-4-plot-m2.dat" u 1:2:3 w e, \
  "< choose REPRESENTATION SSE kagome-4-plot-m2.dat" u 1:2:3 w e
pause -1

set ylabel 'Susceptibility'
plot \
  "< choose REPRESENTATION diagonalization kagome-4-plot-sus.dat" u 1:2 w lp, \
  "< choose REPRESENTATION 'path integral' kagome-4-plot-sus.dat" u 1:2:3 w e, \
  "< choose REPRESENTATION SSE kagome-4-plot-sus.dat" u 1:2:3 w e
pause -1

set ylabel 'Sign'
plot \
  "< choose REPRESENTATION 'path integral' kagome-4-plot-sign.dat" u 1:2:3 w e, \
  "< choose REPRESENTATION SSE kagome-4-plot-sign.dat" u 1:2:3 w e
pause -1
