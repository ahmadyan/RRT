Cmos inverter
.inc PMOS_VTL.inc
.inc NMOS_VTL.inc
.inc INV_X1_lpe_vsub.spi

xinv	vdd 	vss 	vsubp 	vsubn 	vin 	vout	INV_X1
cl		vout 	vss		10f
vsource		vdd		0		1.13019v
vground		vss		0		0.114219v
vsub1		vsubp   0 		0.904874v
vsub2		vsubn 	0 		-0.0221625v
vinput		vin 	0		-0.00433973v


.tran 1fs 1e-012 uic
.save type=ic file=ic_9999.ic level=all time=1e-012
.load file=ic_2533.ic0

.option post 
.print vout
.end