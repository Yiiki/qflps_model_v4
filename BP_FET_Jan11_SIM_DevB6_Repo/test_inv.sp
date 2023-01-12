.title test_basic
.hdl "amos_qflps3.0.va"
.hdl "amos_qflps3.0Dev12.va"
.hdl "amos_qflps3.0Dev23.va"
.hdl "amos_qflps3.0Dev34.va"
.hdl "amos_qflps3.0Dev45.va"
.hdl "amos_qflps3.0Dev56.va"
.hdl "amos_qflps3.0Dev13.va"
.hdl "amos_qflps3.0Dev24.va"
.hdl "amos_qflps3.0Dev35.va"
.hdl "amos_qflps3.0Dev46.va"
.hdl "amos_qflps3.0Dev14.va"
.hdl "amos_qflps3.0Dev25.va"
.hdl "amos_qflps3.0Dev36.va"
.hdl "amos_qflps3.0Dev15.va"
.hdl "amos_qflps3.0Dev26.va"
.hdl "amos_qflps3.0Dev16.va"

.global gnd vdd

.param wx=1u lx=1u

//could change w or l of the ambpfet
x12 g p2 p1 gnd ambpfet12 $ w=wx l=lx 
x23 g p3 p2 gnd ambpfet23 $ w=wx l=lx 
x34 g p4 p3 gnd ambpfet34 $ w=wx l=lx 
x45 g p5 p4 gnd ambpfet45 $ w=wx l=lx 
x56 g p6 p5 gnd ambpfet56 $ w=wx l=lx 
* x13 g d s gnd ambpfet13 $ w=wx l=lx 
* x24 g d s gnd ambpfet24 $ w=wx l=lx 
* x35 g d s gnd ambpfet35 $ w=wx l=lx 
* x46 g d s gnd ambpfet46 $ w=wx l=lx 
* x14 g d s gnd ambpfet14 $ w=wx l=lx 
* x25 g d s gnd ambpfet25 $ w=wx l=lx 
* x36 g d s gnd ambpfet36 $ w=wx l=lx 
* x15 g d s gnd ambpfet15 $ w=wx l=lx 
* x26 g d s gnd ambpfet26 $ w=wx l=lx 
* x16 g d s gnd ambpfet16 $ w=wx l=lx 

vd p6 gnd 3V
vg g gnd 0V
vs p1 gnd 0V

* .dc vg -3V 3V 0.03V vd 1V 3V 1V
.dc vg 0V 3V 0.1V
//.dc vd 0V 5V 0.05V sweep wx 1u 5u 1u $ use this line to sweep wx
.option post=2
.print v(p5)
.print v(p4)
.print v(p3)
.print v(p2)

.end