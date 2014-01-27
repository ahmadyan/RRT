Cmos inverter
.inc PMOS_VTL.inc
.inc NMOS_VTL.inc
.inc INV_X1_lpe.spi

xinv	vdd 	vss 	vin 	vout	INV_X1
cl		vout 	vss		10f
vsource		vdd		0		-0.0244301v
vground		vss		0		0.0453475v
vinput		vin 	0		-0.0617298v

.tran 1fs 5e-012 uic
.save type=ic file=ic_0 level=all time=5e-012
.load file=ic_1
.option post 
.print vout
.end