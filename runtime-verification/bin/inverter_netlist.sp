Cmos inverter
.inc PMOS_VTL.inc
.inc NMOS_VTL.inc
.inc INV_X1_lpe.spi

xinv	vdd 	vss 	vin 	vout	INV_X1
cl		vout 	vss		10f
vsource		vdd		0		1.0593v
vground		vss		0		0.0649648v
vinput		vin 	0		0.932194v

.tran 1fs 1e-012 uic
.save type=ic file=ic_1009.ic level=all time=1e-012
.load file=ic_1008.ic0

.option post 
.print vout
.end