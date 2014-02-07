Cmos inverter
.inc PMOS_VTL.inc
.inc NMOS_VTL.inc
.inc INV_X1_lpe.spi

xinv	vdd 	vss 	vin 	vout	INV_X1
cl		vout 	vss		10f
vsource		vdd		0		0.978353v
vground		vss		0		-0.0678091v
vinput		vin 	0		0.0752007v

.tran 1fs 1e-012 uic
.save type=ic file=ic_1999.ic level=all time=1e-012
.load file=ic_1954.ic0

.option post 
.print vout
.end