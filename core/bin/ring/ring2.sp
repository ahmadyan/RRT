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

vsrc1		vdd1	0       0.9
vsrc2		vdd2	0       0.9
vsrc3		vdd3	0       0.9
vsrc4		vdd4	0       0.9
vsrc5		vdd5	0       0.9
vsrc6		vdd6	0       0.9
vsrc7		vdd7	0       0.9

vgnd1        gnd1    0       0
vgnd2        gnd2    0       0
vgnd3        gnd3    0       0
vgnd4        gnd4    0       0
vgnd5        gnd5    0       0
vgnd6        gnd6    0       0
vgnd7        gnd7    0       0


v1noise  	8 	2  		0
v2noise 	9 	3  		0
v3noise 	10 	4  		0
v4noise 	11 	5  		0
v5noise 	12 	6  		0
v6noise 	13 	7  		0
v7noise 	14 	1  		0

vsub1		vsubp1   0 		0.9
vsub2		vsubp2   0 		0.9
vsub3		vsubp3   0 		0.9
vsub4		vsubp4   0 		0.9
vsub5		vsubp5   0 		0.9
vsub6		vsubp6   0 		0.9
vsub7		vsubp7   0 		0.9

vsub8		vsubn1 	0 		0
vsub9		vsubn2	0 		0
vsub10		vsubn3 	0 		0
vsub11		vsubn4 	0 		0
vsub12		vsubn5 	0 		0
vsub13		vsubn6 	0 		0
vsub14		vsubn7 	0 		0


.options post
.options snaccuracy=50

.tran 1fs 1e-09
.save type=ic file=ring.ic level=all time=0
*.load file=ring1.ic0
.ic v(1)=0.9
.end
