Cmos inverter
.inc PMOS_VTL.inc
.inc NMOS_VTL.inc
.inc INV_X1_lpe.spi

xinv	vdd 	vss 	vin 	vout	INV_X1
cl		vout 	vss		10f
vsource		vdd		0		-0.00804498v
vground		vss		0		0.00449995v
vinput		vin 	0		0.514653v

.tran 1fs 5e-009 uic
.save type=ic file=ic_19.ic level=all time=5e-009
.load file=ic_18.ic0

.option post 
.print vout
.end