#set term dumb size 235,59
set term postscript eps enhanced color round size 8cm,4cm font 10
set output "figure_scripts/Policyfns_211.eps"

set xtics in scale 0.3
set ytics in scale 0.3
set xtics nomirror out

#dat1 = "plot_data/Policyfnk.txt"
#dat2 = "plot_data/Policyfnz.txt"
dat3 = "plot_data/Policyfnz_baseline.txt"
dat4 = "plot_data/Policyfnz_abc.txt"

set multiplot layout 1,1

# set xlabel "intial capital"
# #set yrange [-0.1:0.35]
# set xrange[100:1000]
# set title "Investment rate"
# set xtics 200
# set key b l reverse Left  samplen 1
# plot dat1 index 0 u 1:($10>0.0 ? $7 : 1/0) w lp pi -0.5 pt 71 lw 1 ps 0.3 lc 'black' title "low leverage", \
#      dat1 index 1 u 1:($10>0.0 ? $7 : 1/0) w lp pi -0.5 pt 7  lw 1 ps 0.4 lc 'black' title "high leverage"

# set xtics autofreq
# set xlabel "intial profitability"
# set yrange [*:*]
# set xrange[-2.0:2.0]
# set title "Investment rate"
# set key t l reverse Left  samplen 1
# plot dat2 index 0 u 3:($10>0.0 ? $7 : 1/0) w lp pi -0.5 pt 71 lw 1 ps 0.3 lc 'black' title "low leverage", \
#      dat2 index 1 u 3:($10>0.0 ? $7 : 1/0) w lp pi -0.5 pt 7  lw 1 ps 0.4 lc 'black' title "high leverage"


set xrange[-1.15:2]
set yrange[0.55:0.7]
set key t l samplen 2
set xlabel "initial profitability"
set title "Leverage"
plot dat3 index 0 u ($3):($10>0.0 ? $8 : 1/0) w lp pi -0.5 pt 71 lw 1 ps 0.4 lc 'black' title "baseline", \
     dat4 index 0 u ($3):($10>0.0 ? $8 : 1/0) w lp pi -0.5 pt 76 lw 1 ps 0.4 lc 'black' title "asset-based only"

unset multiplot
