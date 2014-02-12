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

vsrc1		vdd1	0       0.921444v
vsrc2		vdd2	0       0.936279v
vsrc3		vdd3	0       0.934146v
vsrc4		vdd4	0       0.864173v
vsrc5		vdd5	0       0.901012v
vsrc6		vdd6	0       0.866309v
vsrc7		vdd7	0       0.923208v

vgnd1        gnd1    0       0.00347453v
vgnd2        gnd2    0       0.0318201v
vgnd3        gnd3    0       0.0232627v
vgnd4        gnd4    0       0.0176992v
vgnd5        gnd5    0       -0.0164876v
vgnd6        gnd6    0       0.0224448v
vgnd7        gnd7    0       -0.00467696v

v1noise  	8 	2  		-0.00386868v
v2noise 	9 	3  		-0.000153966v
v3noise 	10 	4  		-0.00353847v
v4noise 	11 	5  		0.00344996v
v5noise 	12 	6  		-0.00335505v
v6noise 	13 	7  		-0.00440886v
v7noise 	14 	1  		-0.00487487v

vsub1		vsubp1   0 		0.865015v
vsub2		vsubp2   0 		0.901473v
vsub3		vsubp3   0 		0.877549v
vsub4		vsubp4   0 		0.886863v
vsub5		vsubp5   0 		0.914599v
vsub6		vsubp6   0 		0.90945v
vsub7		vsubp7   0 		0.905531v
vsub8		vsubn1 	0 		0.05v
vsub9		vsubn2	0 		0.05v
vsub10		vsubn3 	0 		0.05v
vsub11		vsubn4 	0 		0.05v
vsub12		vsubn5 	0 		0.05v
vsub13		vsubn6 	0 		0.05v
vsub14		vsubn7 	0 		0.05v


.options post
.options snaccuracy=50

.tran 1fs 1e-012 uic
.save type=ic file=ic_509.ic level=all time=1e-012
.load file=ic_508.ic0

.end
