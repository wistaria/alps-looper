#!/bin/sh

./diag < xy-z.ip | sh extract.sh > xy-z.dat
./diag < xy-x.ip | sh extract.sh > xy-x.dat
cat <<EOF > xy.dat
# LATTICE = "chain lattice"
# MODEL   = "spin"
# L = 8
# Jx = 1
# Jy = 1
# Jz = 0
# column  1: temperature
#         2: energy per site
#         3: specific heat
#         4: uniform magnetization^2 in Z direction
#         5: uniform susceptibility in Z direction
#         6: staggered magnetization^2 in Z direction
#         7: staggered susceptibility in Z direction
#         8: uniform magnetization^2 in X direction
#         9: uniform susceptibility in X direction
#        10: staggered magnetization^2 in X direction
#        11: staggered susceptibility in X direction
EOF
join xy-z.dat xy-x.dat | awk '{print $1,$2,$3,$4,$5,$6,$7,$10,$11,$12,$13}' >> xy.dat
rm -f xy-z.dat xy-x.dat

cat <<EOF > heisenberg.dat
# LATTICE = "chain lattice"
# MODEL   = "spin"
# L = 8
# Jx = -1
# Jy = -1
# Jz = -1
# column  1: temperature
#         2: energy per site
#         3: specific heat
#         4: uniform magnetization^2
#         5: uniform susceptibility
#         6: staggered magnetization^2
#         7: staggered susceptibility
EOF
./diag < heisenberg.ip | sh extract.sh >> heisenberg.dat

./diag < ising-z.ip | sh extract.sh > ising-z.dat
./diag < ising-x.ip | sh extract.sh > ising-x.dat
cat <<EOF > ising.dat
# LATTICE = "chain lattice"
# MODEL   = "spin"
# L = 8
# Jx = -1
# Jy = -1
# Jz = -2
# column  1: temperature
#         2: energy per site
#         3: specific heat
#         4: uniform magnetization^2 in Z direction
#         5: uniform susceptibility in Z direction
#         6: staggered magnetization^2 in Z direction
#         7: staggered susceptibility in Z direction
#         8: uniform magnetization^2 in X direction
#         9: uniform susceptibility in X direction
#        10: staggered magnetization^2 in X direction
#        11: staggered susceptibility in X direction
EOF
join ising-z.dat ising-x.dat | awk '{print $1,$2,$3,$4,$5,$6,$7,$10,$11,$12,$13}' >> ising.dat
rm -f ising-z.dat ising-x.dat
