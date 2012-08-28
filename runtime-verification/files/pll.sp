Phase-Locked Loop

.include param.sp
.option  probe acct
.option relv=1e-5
*
* wideband FM example, Grebene gives:
*  f0=1meg kf=250kHz/V
*  kd=0.1 V/rad
*  R=10K C=1000p
*  f_lock = kf*kd*pi/2 = 39kHz, v_lock = kd*pi/2 = 0.157
*  f_capture/f_lock ~= 1/sqrt(2*pi*R*C*f_lock)
*                    = 0.63,    v_capture ~= 0.100

*******************
* Input Waveforms *
*******************

* These are the input sinusoids at different frequencies
* Hint: You can use the .alter statement to do three versions
* of the transient analysis in one run.

*vin in inb sin(0 0.6 1.00x)
vin in inb sin(0 0.6 1.02x)
*vin in inb sin(0 0.6 0.95x)
.tran 0.1u 5u uic
.SAVE TYPE=.IC LEVEL=ALL TIME=5u
.load

* FM Signal for Demodulation
*vin in inb sffm(0 0.6 1.0x 1 10k)
*.tran 0.1u 1000u

******************************************
*************
* VCO Model *
*************
xvco e eb osc oscb vco f0=1x kf=125k phi=0 out_off=-1 out_amp=0.3

************************
* Phase Detector Model *
************************
xpd in inb osc oscb mout moutb pd kd=0.1 out_off=-3.5

**********
* Filter *
**********

rf mout e 10k
cf  e   0 1000p
rfb moutb eb 10k
cfb eb   0 1000p

****************
* Final Output *
****************
rout out e 100k
cout out 0 100p
routb outb eb 100k
coutb outb 0 100p


.macro vco in inb out outb f0=100k kf=50k phi=0.0 out_off=0.0 out_amp=1.0
gs 0 s poly(2) c 0 in inb    0 '6.2832e-9*f0' 0 0 '6.2832e-9*kf'
gc c 0 poly(2) s 0 in inb    0 '6.2832e-9*f0' 0 0 '6.2832e-9*kf'
cs s 0 1e-9
cc c 0 1e-9
e1 s_clip 0 pwl(1) s 0  -0.1,-0.1 0.1,0.1
e out 0 s_clip 0 out_off '10*out_amp'
eb outb 0 s_clip 0 out_off '-10*out_amp'2
.eom

.macro pd in inb in2 in2b out outb kd=0.1 out_off=0
e1 clip1 0 pwl(1) in inb  -0.1,-0.1 0.1,0.1
e2 clip2 0 pwl(1) in2 in2b  -0.1,-0.1 0.1,0.1
*e3 n1 0 poly(2) clip1 0 clip2 0    0 0 0 0 '78.6*kd'
e3 n1 0 poly(2) clip1 0 clip2 0    var 0 0 0 '78.6*kd'
e4 outb 0 n1 0  out_off 1
e5 out 0 n1 0  out_off -1
.eom

.option post delmax=0.01u interp
.probe v_in=v(in,inb) v_control=v(e,eb) v_vco=v(osc,oscb)

.end
