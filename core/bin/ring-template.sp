.title ringosc

.inc PMOS_VTL.inc
.inc NMOS_VTL.inc
.inc INV_X1_lpe_vsub.spi


x1 vdd1 gnd1 vsubp1 vsubn1 1 8   INV_X1
x2 vdd2 gnd2 vsubp2 vsubn2 2 9   INV_X1
x3 vdd3 gnd3 vsubp3 vsubn3 3 10  INV_X1
x4 vdd4 gnd4 vsubp4 vsubn4 4 11  INV_X1
x5 vdd5 gnd5 vsubp5 vsubn5 5 12  INV_X1
x6 vdd6 gnd6 vsubp6 vsubn6 6 13  INV_X1
x7 vdd7 gnd7 vsubp7 vsubn7 7 14  INV_X1

c1 1 0 1fF
c2 2 0 1fF
c3 3 0 1fF
c4 4 0 1fF
c5 5 0 1fF
c6 6 0 1fF
c7 7 0 1fF

vsrc1		vdd1	0       $PARAM_0v
vsrc2		vdd2	0       $PARAM_1v
vsrc3		vdd3	0       $PARAM_2v
vsrc4		vdd4	0       $PARAM_3v
vsrc5		vdd5	0       $PARAM_4v
vsrc6		vdd6	0       $PARAM_5v
vsrc7		vdd7	0       $PARAM_6v

vgnd1        gnd1    0       $PARAM_7v
vgnd2        gnd2    0       $PARAM_8v
vgnd3        gnd3    0       $PARAM_9v
vgnd4        gnd4    0       $PARAM_10v
vgnd5        gnd5    0       $PARAM_11v
vgnd6        gnd6    0       $PARAM_12v
vgnd7        gnd7    0       $PARAM_13v

v1noise  	8 	2  		$PARAM_14v
v2noise 	9 	3  		$PARAM_15v
v3noise 	10 	4  		$PARAM_16v
v4noise 	11 	5  		$PARAM_17v
v5noise 	12 	6  		$PARAM_18v
v6noise 	13 	7  		$PARAM_19v
v7noise 	14 	1  		$PARAM_20v

vsub1		vsubp1   0 		$PARAM_21v
vsub2		vsubp2   0 		$PARAM_22v
vsub3		vsubp3   0 		$PARAM_23v
vsub4		vsubp4   0 		$PARAM_24v
vsub5		vsubp5   0 		$PARAM_25v
vsub6		vsubp6   0 		$PARAM_26v
vsub7		vsubp7   0 		$PARAM_27v
vsub8		vsubn1 	0 		$PARAM_28v
vsub9		vsubn2	0 		$PARAM_29v
vsub10		vsubn3 	0 		$PARAM_30v
vsub11		vsubn4 	0 		$PARAM_31v
vsub12		vsubn5 	0 		$PARAM_32v
vsub13		vsubn6 	0 		$PARAM_33v
vsub14		vsubn7 	0 		$PARAM_34v


.options post
.options snaccuracy=50

$tran
$save
$load

.end
