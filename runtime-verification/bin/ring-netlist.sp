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

vsrc1		vdd1	0       0.965484v
vsrc2		vdd2	0       0.97047v
vsrc3		vdd3	0       0.868526v
vsrc4		vdd4	0       0.972698v
vsrc5		vdd5	0       0.952947v
vsrc6		vdd6	0       0.970519v
vsrc7		vdd7	0       0.832991v

vgnd1        gnd1    0       -0.0704031v
vgnd2        gnd2    0       0.0989074v
vgnd3        gnd3    0       0.0799188v
vgnd4        gnd4    0       0.0711966v
vgnd5        gnd5    0       -0.0088168v
vgnd6        gnd6    0       -0.0405744v
vgnd7        gnd7    0       -0.058916v

v1noise  	8 	2  		0.0119739v
v2noise 	9 	3  		-0.0346034v
v3noise 	10 	4  		-0.0364803v
v4noise 	11 	5  		-0.0306787v
v5noise 	12 	6  		-0.0246208v
v6noise 	13 	7  		-0.035107v
v7noise 	14 	1  		-0.0109424v

vsub1		vsubp1   0 		0.982269v
vsub2		vsubp2   0 		0.806305v
vsub3		vsubp3   0 		0.973034v
vsub4		vsubp4   0 		0.911258v
vsub5		vsubp5   0 		0.816449v
vsub6		vsubp6   0 		0.9361v
vsub7		vsubp7   0 		0.987036v
vsub8		vsubn1 	0 		0.0579943v
vsub9		vsubn2	0 		0.0135289v
vsub10		vsubn3 	0 		-0.0775567v
vsub11		vsubn4 	0 		0.0268715v
vsub12		vsubn5 	0 		0.0469405v
vsub13		vsubn6 	0 		0.00610675v
vsub14		vsubn7 	0 		-0.0878353v


.options post
.options snaccuracy=50

.tran 1fs 1e-012 uic
.save type=ic file=ic_24999.ic level=all time=1e-012
.load file=ic_2058.ic0

.end
