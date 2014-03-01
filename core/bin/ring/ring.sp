.title ringosc

.inc PMOS_VTL.inc
.inc NMOS_VTL.inc
.inc INV_X1_lpe_vsub.spi


x1 vdd vss vsubp vsubn 1 2  INV_X1
v1noise  2 3  0.0
x2 vdd vss vsubp vsubn 3 4  INV_X1
v2noise 4 5  0.0
x3 vdd vss vsubp vsubn 5 6  INV_X1
v3noise 6 7  0.0
x4 vdd vss vsubp vsubn 7 8  INV_X1
v4noise 8 9  0.0
x5 vdd vss vsubp vsubn 9 10  INV_X1
v5noise 10 11  0.0
x6 vdd vss vsubp vsubn 11 12  INV_X1
v6noise 12 13  0.0
x7 vdd vss vsubp vsubn 13 14  INV_X1
v7noise 14 1  0.0

c1 1 0 1fF
c2 2 0 1fF
c3 3 0 1fF
c4 4 0 1fF
c5 5 0 1fF
c6 6 0 1fF
c7 7 0 1fF

vsource		vdd		0		0.9
vground		vss		0		0
vsub1		vsubp   0 		0.9
vsub2		vsubn 	0 		0


.options post
.options snaccuracy=50

.tran 1fs 1e-09 uic
.save type=ic file=ring2.ic level=all time=13ps
.load file=ring1.ic0
.end
