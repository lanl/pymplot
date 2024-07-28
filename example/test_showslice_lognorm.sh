
x_showslice \
  -in=./data/slicesqr.bin -out=test_showslice_lognorm.pdf \
  -n1=100 -n2=200 -n3=200 -norm=log \
  -slice1=0.55 -slice2=0.11 -slice3=0.3 \
  -x2end=0.6 -x2beg=-0.6 -x1beg=0.0 -cmin=9 -cmax=14 \
  -x3end=3.0 \
  -size1=2 -size2=2.5 -size3=4.0 \
  -color=kgbwr \
  -label1='Depth (km)' \
  -label2='Inline position (km)' \
  -label3='Crossline position (km)' \
  -label1size=16 -label2size=16 -label3size=16 \
  -d1=0.01 -d2=0.008 -d3=0.008 \
  -tick1d=0.25 -tick2d=0.2 -tick3d=0.5 \
  -o1=0 -o2=-0.6 -o3=0 \
  -tick2beg=0 -tick2beg=-0.4 -tick3beg=0 \
  -mtick1=4 -mtick2=4 -mtick3=4 \
  -legend=1 -unit='Test Log (m/s)' -lmtick=10 \
  -dpi=400 -sliceline1=on -sliceline2=on -sliceline3=on \
  -sliceline1color=yellow -sliceline2color=k -sliceline3color=b \
  -lloc=bottom -sliceline3style=: -color=gist_ncar

evince test_showslice_lognorm.pdf &
