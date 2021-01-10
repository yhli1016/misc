import os
from sandbox import SandBox


s = SandBox()
home = os.environ["HOME"]

s.set_mod(mod_name="qe-6.6")
s.set_mod("bin", home+"/soft/dft/qe-6.6/bin")
s.reset("OMP_NUM_THREADS", 1)

s.echo_commands()
