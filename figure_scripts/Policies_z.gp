set term dumb size 218,55

set xtics in scale 0.00001
set ytics in
set xtics nomirror out

dat1 = "plot_data/Policyfnz.txt"

set multiplot layout 2,3

set xlabel "intial profitability"
set xrange[-2.0:2.0]

set title "Investment"
plot dat1 index 0 u 3:($10>0.0 ? $6 : 1/0) w lp notitle "low leverage", \
     dat1 index 1 u 3:($10>0.0 ? $6 : 1/0) w lp notitle "high leverage", \

set title "Investment rate"
plot dat1 index 0 u 3:($10>0.0 ? $7 : 1/0) w lp notitle "low leverage", \
     dat1 index 1 u 3:($10>0.0 ? $7 : 1/0) w lp notitle "low leverage"

set title "Divident"
plot dat1 index 0 u 3:($10>0.0 ? $9 : 1/0) w lp notitle "low leverage", \
     dat1 index 1 u 3:($10>0.0 ? $9 : 1/0) w lp notitle "low leverage"

set yrange[*:*]
set title "Debt"
plot dat1 index 0 u 3:($10>0.0 ? $5 : 1/0) w lp notitle "low leverage", \
     dat1 index 1 u 3:($10>0.0 ? $5 : 1/0) w lp notitle "low leverage"

set title "Leverage"
plot dat1 index 0 u 3:($10>0.0 ? $8 : 1/0) w lp notitle "low leverage", \
     dat1 index 1 u 3:($10>0.0 ? $8 : 1/0) w lp notitle "low leverage"

set yrange[*:*]
set key b r
set title "Value"
plot dat1 index 0 u 3:($10>0.0 ? $10 : 1/0) w lp title "low leverage", \
     dat1 index 1 u 3:($10>0.0 ? $10 : 1/0) w lp title "high leverage"

unset multiplot
