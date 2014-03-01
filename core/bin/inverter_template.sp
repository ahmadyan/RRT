Cmos inverter
.inc PMOS_VTL.inc
.inc NMOS_VTL.inc
.inc INV_X1_lpe_vsub.spi

xinv	vdd 	vss 	vsubp 	vsubn 	vin 	vout	INV_X1
cl		vout 	vss		10f
vsource		vdd		0		$PARAM_0v
vground		vss		0		$PARAM_1v
vsub1		vsubp   0 		$PARAM_2v
vsub2		vsubn 	0 		$PARAM_3v
vinput		vin 	0		$PARAM_4v


$tran
$save
$load

.option post 
.print vout
.end