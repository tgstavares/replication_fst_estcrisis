set term dumb size 218,55
#set term postscript eps enhanced color round size 10cm,4cm font 10
#set output "figure_scripts/Lev_z_comp_less_mean15sd_zneg1075566.eps"

set xtics in scale 0.00001
set ytics in
set xtics nomirror out

dat1 = "plot_data/Policyfnz_baseline.txt"
dat2 = "plot_data/Policyfnz_abc.txt"

set multiplot layout 1,1

set xlabel "change in profitability {/Symbol \104}z"


# set title "Investment rate"
# plot dat1 index 0 u 3:($10>0.0 ? $7 : 1/0) w lp notitle "baseline", \
#      dat2 index 1 u 3:($10>0.0 ? $7 : 1/0) w lp notitle "abc"


z0 = -1.075566
#set xrange[-1.0-z0:2.0-z0]
set xrange[-2:2]
set yrange[]
set key t l
#set title "Leverage"
plot dat1 index 0 u ($3):($10>0.0 ? $8 : 1/0) w lp pi -1 pt 71 lw 1 ps 0.5 lc 'black' title "baseline", \
     dat2 index 0 u ($3):($10>0.0 ? $8 : 1/0) w lp pi -1 pt 7  lw 1 ps 0.5 lc 'black' title "collateral only constraint"

unset multiplot
