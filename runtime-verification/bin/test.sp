Cmos inverter
.inc PMOS_VTL.inc
.inc NMOS_VTL.inc
.inc INV_X1_lpe.spi

xinv	vdd 	vss 	vin 	vout	INV_X1
cl		vout 	vss		10f
vsource		vdd		0		$PARAM_0
vground		vss		0		$PARAM_1
vinput		vin 	0		$PARAM_2 

.tran 1fs $DT uic
.save type=nodeset file=$SAVE_FILE level=all time=$DT
.load file=$LOAD_FILE
.option post 
.end