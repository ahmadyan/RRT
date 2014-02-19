* passive half-wave limiter with low pass filter

.model D1N4148	D 	(Is=0.1pA Rs=16 CJO=2p Tt=12n Bv=100 Ibv=0.1p)

*the main signal without any noise
v1	1	0	$PARAM_0
inoise  3   0   $PARAM_1

*passive low pass filter
r1	1	2	1k
c1  2   0   6nF

* passive half-wave limiter
r2	2	3	1k
d2	0 	3 	D1N4148

*voltage divier load
r3 3  4  5k
rl 4  0  5k


$tran
$load
$save

.option post delmax=0.01u interp
.option  probe acct
.option relv=1e-5
.end