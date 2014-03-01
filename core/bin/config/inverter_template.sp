Cmos inverter
.inc PMOS_VTL.inc
.inc NMOS_VTL.inc
.inc INV_X1_lpe.spi

xinv	vdd 	vss 	vin 	vout	INV_X1
cl		vout 	vss		10f
vsource		vdd		0		$PARAM_0v
vground		vss		0		$PARAM_1v
vinput		vin 	0		$PARAM_2v

$tran
$save
$load

.option post 
.print vout
.end