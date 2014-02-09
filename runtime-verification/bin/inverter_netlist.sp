Cmos inverter
.inc PMOS_VTL.inc
.inc NMOS_VTL.inc
.inc INV_X1_lpe_vsub.spi

xinv	vdd 	vss 	vsubp 	vsubn 	vin 	vout	INV_X1
cl		vout 	vss		10f
vsource		vdd		0		0.948845v
vground		vss		0		0.117307v
vsub1		vsubp   0 		1.01868v
vsub2		vsubn 	0 		0.16599v
vinput		vin 	0		1.95474v


.tran 1fs 1e-011 uic
.save type=ic file=ic_999.ic level=all time=1e-011
.load file=ic_0.ic0

.option post 
.print vout
.end