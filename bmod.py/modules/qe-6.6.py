import os
from sandbox import SandBox


s = SandBox()
home = os.environ["HOME"]

s.set_mod(mod_name="qe-6.6")
s.set_mod("pkg", home+"/soft/libxc-4.3.4")
s.set_mod("pkg", home+"/soft/libbeef")
s.set_mod("pkg", home+"/soft/elpa-2018.11.001", cmd="rm")
s.set_mod("pkg", home+"/soft/elpa-2019.11.001")
s.set_mod("bin", home+"/soft/qe-6.6/bin")
s.reset("OMP_NUM_THREADS", 1)

s.echo_commands()
