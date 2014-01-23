Cmos inverter
.inc PMOS_VTL.inc
.inc NMOS_VTL.inc
.inc INV_X1_lpe.spi

xinv	vdd 	vss 	vin 	vout	INV_X1
cl		vout 	vss		10f
* Rl 		vout	vss		1k
vsource		vdd		0		1
vground		vss		0		0
vinput		vin 	0		pulse( 0 1 5ps 5ps 5ps 40ps 100ps )

.dc vIN start=0 stop=2.5 step=0.01
.tran 1fs 1ns uic
*.save type=nodeset file=test.ic0 level=all time=420ps
.load file=test.ic0
.option post 
.end