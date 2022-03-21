#set term dumb size 218,59
set term postscript eps enhanced color round size 12cm,4cm font 10
set output "figure_scripts/Pricedebt_defcf_3.eps"

set xtics in scale 0.3
set ytics in scale 0.3
set xtics nomirror out

dat1 = "plot_data/Pricedebtfnp.txt"

set multiplot layout 1,3

set xlabel "next period leverage"

set title "Probability of repayment"
set xrange[0.35:0.95]
set yrange[0.0:1]
set key b l reverse Left samplen 0.5

plot dat1 index 0 u 4:8 w p pt 71 ps 0.4 lc 'black' title "high profitability", \
    dat1 index 1 u 4:8 w p pt 7  ps 0.5 lc 'black' dt(1,1) title "low profitability"

set title "Debt price"
plot dat1 index 0 u 4:5 w p pt 71 ps 0.4 lc 'black' notitle "high profitability", \
     dat1 index 1 u 4:5 w p pt 7  ps 0.5 lc 'black' dt(1,1) notitle "low profitability"

set yrange[*:*]
set title "Credit gross inflow"
plot dat1 index 0 u 4:6 w p pt 71 ps 0.4 lc 'black' notitle "high profitability", \
     dat1 index 1 u 4:6 w p pt 7  ps 0.5 lc 'black' dt(1,1) notitle "low profitability"


unset multiplot
