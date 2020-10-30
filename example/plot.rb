
#-------------------------------------------------------------------------------
# 2D matrix plot -- Mamrousi model
opts = "-in=./data/mar_vp.bin -n1=651 -d1=0.005 -d2=0.005 -label1='Depth (km)' \
-label2='Horizontal Position (km)' -mtick1=9 -mtick2=9 -tick1d=1 -tick2d=1 \
-size1=2 -size2=4.5 -legend=1 -legendloc=bottom -unit='P-wave Velocity (m/s)' -lmtick=9 \
-ltickbeg=500 -ld=1000"
system "x_showmatrix #{opts} -color=jet -font=arial -unit='-color=jet -font=arial' -out=test_mar_1.pdf & "
system "x_showmatrix #{opts} -color=gist_ncar -font=consolas -unit='-color=gist_ncar -font=consolas' -out=test_mar_2.pdf & "
system "x_showmatrix #{opts} -color=bwr -font=times -unit='-color=bwr -font=times' -out=test_mar_3.pdf & "
system "x_showmatrix #{opts} -color=binary -font=courier -unit='-color=binary -font=courier' -out=test_mar_4.pdf & "
system "x_showmatrix #{opts} -color=viridis -font=plex -unit='-color=viridis -font=plex' -out=test_mar_5.pdf & "
system "x_showmatrix #{opts} -color=rainbow256 -font=helvetica -unit='-color=rainbow256 -font=helvetica' -out=test_mar_6.pdf & "

#-------------------------------------------------------------------------------
# 2D contour plot with background image
system "x_showcontour -in=./data/tt_random.bin -n1=320 -d1=0.01 -d2=0.01 -label1='Axis 1' \
-label2='Flexible Tick Labels' -ticks2=0,'Start':1,'$\\alpha$':2.7,'$G(\\omega)$' \
-contourfill=1 -legend=1 -color=rainbow -overlay=1 -interp=Gaussian -ncolor=12 -contourbeg=0 -contourend=0.6 \
-tick1d=1 -tick2d=0.5 -mtick1=9 -mtick2=4 -cmax=0.6 -unit='Unit with Super/Subscripts (N$\\mathregular{_{\\theta}}$/m$\\mathregular{^2})$' \
-size1=4 -contourlevel=0.1 -contourwidth=2 -mcontour=1 -out=./test_contour_time.pdf & "

#-------------------------------------------------------------------------------
# 2D contour plot with irregular topography
system "x_showcontour -in=./data/water_time_surface.bin -n1=610 -color=Blues -size1=4 -overlay=1 \
-legend=1 -contours=0.1,0.3,0.4,1.0,1.2,1.5,1.7,2.0 -mcontour=4 -mcontourwidth=0.5 -contourwidth=1.5 \
-unit='First Arrial Traveltime (s)' -mtick1=9 -mtick2=9 -f1=1000 -tick1d=-1000 -tick2d=1000 -lmtick=9 \
-label1='Elevation (m)' -label2='Horizontal Position (m)' -tick1rot=90 -tick2rot=30 -d1=-8 -d2=8 -clabelcolor=r \
-color=Blues -out=test_contour_surface.pdf & "

#-------------------------------------------------------------------------------
# Wiggle plots
system "x_showwiggle -in=./data/seisdata.bin,./data/seisdata_smooth.bin -n1=530 \
-d2=0.000995 -d1=0.02 -mtick1=1 -label2='Time (s)' -label1='Horizontal Position (km)' \
-every=25 -fill=0 -along=2 -wigglecolor=b,r -wigglewidth=1,2 \
-reverse1=1 -x2end=3.2 -clip=2 -size1=4.5 -size2=5.5 -plotlabel='Original','Filtered' \
-plotlabelloc=upper_right -x1end=5.5 -reverse1=0 -tick2d=0.5 -mtick2=4 -x2end=3 \
-polygonalpha=0.7 -polygon=0,0.1,0,1.2,4.7,0.1 -polygonfacecolor=lime -polygonedgecolor=k \
-out=test_wiggle_horizontal.pdf &"

system "x_showwiggle -in=./data/seisdata.bin -out=test_wiggle_single_1.pdf -n1=4001 \
-d1=0.000995 -d2=0.02 -label1='Time (s)' -label2='Horizontal Position (km)' \
-size1=4 -size2=6 -fill=1 -transpose=1 -along=1 -every=5 -x1end=3 \
-x2end=6 -clip=2 -wigglecolor=k -tick1d=0.5 -mtick1=4 &"

system "x_showwiggle -in=./data/seisdata.bin -out=test_wiggle_single_2.pdf -n1=4001 \
-d1=0.000995 -d2=0.02 -label1='Time (s)' -label2='Horizontal Position (km)' \
-size1=4 -size2=6 -fill=0 -transpose=1 -along=1 -every=20 -x1end=2 \
-x2end=6 -clip=4 -wigglewidth=2 -wigglecolor=r -tick1d=0.5 -backcolor=binary \
-backcmin=-2 -interp=Gaussian -backcmax=2 -background=./data/seisdata.bin -mtick1=4 &"

#-------------------------------------------------------------------------------
# 2D matrix plot -- demonstration of fonts and colors
opts = "-in=./data/mar_vp.bin -n1=651 -d1=0.005 -d2=0.005 -label1='Depth (km)' \
-label2='Horizontal Position (km)' -mtick1=9 -mtick2=9 -tick1d=1 -tick2d=1 \
-size1=2 -size2=4.5 -legend=1 -legendloc=bottom -unit='P-wave Velocity (m/s)' -lmtick=9 \
-ltickbeg=500 -ld=1000"
system "x_showmatrix #{opts} -color=jet -font=arial -unit='-color=jet -font=arial' -out=test_mar_1.pdf & "
system "x_showmatrix #{opts} -color=gist_ncar -font=consolas -unit='-color=gist_ncar -font=consolas' -out=test_mar_2.pdf & "
system "x_showmatrix #{opts} -color=bwr -font=times -unit='-color=bwr -font=times' -out=test_mar_3.pdf & "
system "x_showmatrix #{opts} -color=binary -font=courier -unit='-color=binary -font=courier' -out=test_mar_4.pdf & "

#-------------------------------------------------------------------------------
# 1D data
system "x_showgraph -in=./data/gv.bin -out=test_gvgraph.pdf -transpose=1 -n1=360,360 -n2=4 \
-select=2,3 -ptype=2 -marker=v,o -markersize=4,4 -x1beg=-6500 -x1end=6500 -x2beg=-6500 \
-x2end=6500 -tick1beg=-8000 -tick1d=2000 -tick2beg=-8000 -tick2d=2000 -mtick1=4 -mtick2=4 \
-size1=4.5 -size2=4.5  -label1='Horizontal Velocity (m/s)' -label2='Vertical Velocity (m/s)' \
-close=1,1 -linestyle=-,- -plotlabel='qP':'qSV' &"

system "x_showgraph -in=./data/randdata.bin -n1=100 -ptype=3 -out=test_randgraph.pdf \
 -color=jet -legend=1 -x1beg=10 -x2beg=10 -x1end=1000 -x2end=1000 \
-grid1=on -grid2=on -mgrid2=on -mtick1=5 -mtick2=5 -size1=4 \
-size2=4 -grid1color=gray -grid2color=gray -mgrid1color=gray -color=magma \
-mgrid2color=gray -markersizemin=40 -markersizemax=240 -mgrid1style=: \
-mgrid2style=: -arrow=210,410,810,450 -arrowstyle='<|-|>' \
-arrowconnect='arc3,rad=-0.3' -arrowfacecolor=r -arrowedgecolor=r \
-text='$\\Delta \\omega=80$%' -textloc=210,30 \
-textrotation=-45 -textsize=18 -arroworder=9 -textorder=10 \
-polygon=20.0,20.0,100.0,30.0,100.0,100.0,45.0,90.0 -unit='Scatter Values in magma Colormap' \
-ltickbeg=0 -ld=200 -lmtick=4 -label1='Natural Axis 1' -label2='Logarithmic Axis 2' \
-curve=./data/rand1.txt -curvestyle=line-- -curvecolor=b,r -norm2=log "

#-------------------------------------------------------------------------------
# 3D volume contours
system "x_showvolcon -in=./data/gauss_scaled.bin -label1='$\\omega_1$' \
-label3='$\\omega_3$' -label2='$\\omega_2$' \
-n1=100 -n2=200 -contourfill=1 -angle=40,15 -color=Spectral \
-slice1=-70 -slice2=210 -slice3=-400 -mtick1=2 -mtick2=9 -mtick3=9 \
-d1=-1.221 -d2=2.33 -d3=-4.669 -clabelcolor=k \
-clabelsize=16 -contourwidth=2 -mcontourstyle=: -mcontour=2 \
-octant=--+ -out=test_volcon_2.pdf -norm=log &"

system "x_showvolcon -in=./data/gauss.bin \
-n1=100 -n2=200 -contourfill=1 -angle=20,10 -color=inferno \
-contourlevel=0.1 -contourbeg=0.0 \
-slice1=-70 -slice2=210 -slice3=-400 \
-d1=-1.221 -d2=2.33 -d3=-4.669 -clabelcolor=blue \
-clabelsize=14 -mcontourstyle=: -mcontour=2 \
-octant=--+ -out=test_volcon_1.pdf &"

system "x_showslicon -in=./data/time_200x200x200.bin -n1=200 -n2=200 \
-slicegap=0.15 -d1=0.02 -d2=0.02 -d3=0.02 -label1='Depth (km)' \
-label2='Crossline Position (km)' -label3='Inline Position (km)' \
-tick1d=1 -mtick1=9 -tick2d=3 -mtick2=5 -tick3d=3 -mtick3=5 \
-slice1=1 -slice2=2 -slice3=3 -lmtick=9 -vol3d=over3d.png \
-size1=4 -size2=4 -size3=4 -label1size=18 -label2size=18 -label3size=18 \
-unitsize=18 -contourlevel=0.05 -contourbeg=0 -contourfill=1 \
-vol3d=contour3d.png -clabelcolor=w -color=magma -out=test_slice_2.pdf &"

#-------------------------------------------------------------------------------
# 2D matrix plot with annotations -- BP velocity model
system "x_showmatrix \
-in=./data/bp_vp.bin -out=test_bp.pdf \
-n1=1911 -label1='Depth (km)' -label2='Horizontal Position (km)' \
-size1=3 -size2=7 \
-d1=0.00625 -d2=0.0125 \
-tick1d=2.5 -tick2d=5.0 \
-mtick1=4 -mtick2=4 \
-legend=1 -lmtick=4 -polygon=0.6,58,4,58,4,63,0.6,63 -polygonfacecolor=none \
-polygonedgecolor=w -polygonlinewidth=1.5 -polygonalpha=1 -polygonlinestyle=dashed \
-circle=2.75,27,5,5,0 -circlefacecolor=none -circlealpha=1 -circleedgecolor=k -circlelinewidth=3 \
-text='Object 1':'Object 2' -textloc=2.5,20:0.75,52 -textcolor=w,yellow \
-arrow=0,47,7.5,39 -arrowconnect='arc3,rad=-0.3' -arrowwidth=2 -arrowedgecolor=w \
-unit='P-wave velocity (m/s)' -lwidth=0.15 & "

#-------------------------------------------------------------------------------
# The following datasets are too large to attach in the repository. 
# You can use your own copies of these models -- they are both publicly available. 
# Or you can use your own favorite models -- remember to adjust n1, n2 and n3. 

abort

# SEG/EAGE 3D overthrust model
system "x_showslice -in=./data/overthrust_187x801x801.bin -n1=187 -n2=801 \
-slicegap=0.15 -d1=0.02 -d2=0.02 -d3=0.02 -label1='Depth (km)' \
-label2='Crossline Position (km)' -label3='Inline Position (km)' \
-tick1d=1 -mtick1=9 -tick2d=3 -mtick2=5 -tick3d=3 -mtick3=5 \
-legend=1 -unit='P-wave Velocity (m/s)' -ctruncbeg=0.05 -ctruncend=0.95 \
-slice1=1 -slice2=9 -slice3=6 -lmtick=9 \
-size1=2 -size2=5 -size3=5 -label1size=18 -label2size=18 -label3size=18 \
-unitsize=18 -out=test_slice_1.pdf -vol3d=over3d.png &"

# SEG/EAGE 3D salt model
system "x_showvolume -in=./data/salt_201x676x676.bin -n1=201 -n2=676 \
-d1=0.02 -d2=0.02 -d3=0.02 -label1='Z (km)' \
-label2='Y (km)' -label3='X (km)' -octant=-++ \
-tick1d=1 -mtick1=9 -tick2d=3 -mtick2=5 -tick3d=3 -mtick3=5 \
-legend=1 -unit='P-wave Velocity (m/s)' -ctruncbeg=0.05 -ctruncend=0.95 \
-slice1=3 -slice2=3 -slice3=6 -lmtick=9 -label1size=18 -label2size=18 -label3size=18 \
-unitsize=18 -legendloc=right -x1beg=0.5 -tick1beg=0 -size1=2.5 -size2=5 -size3=5 \
-angle=20,30 -out=test_vol_2.pdf &"

system "x_showvolume -in=./data/salt_201x676x676.bin -n1=201 -n2=676 \
-d1=0.02 -d2=0.02 -d3=0.02 -label1='Z (km)' \
-label2='Y (km)' -label3='X (km)' \
-tick1d=1 -mtick1=9 -tick2d=3 -mtick2=5 -tick3d=3 -mtick3=5 \
-legend=1 -unit='P-wave Velocity (m/s)' -ctruncbeg=0.05 -ctruncend=0.95 \
-slice1=3 -slice2=9 -slice3=6 -lmtick=9 -label1size=18 -label2size=18 -label3size=18 \
-unitsize=18 -legendloc=right -size1=2.5 -size2=5 -size3=5 -x1beg=0.5 \
-tick1beg=0 -out=test_vol_1.pdf &"


