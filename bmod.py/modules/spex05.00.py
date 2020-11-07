import os
from sandbox import SandBox


s = SandBox()
home = os.environ["HOME"]

s.set_mod(mod_name="spex05.00")
s.set_mod("bin", home+"/soft/spex05.00/bin")
s.reset("OMP_NUM_THREADS", 1)

s.echo_commands()
