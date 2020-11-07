import os
from sandbox import SandBox


s = SandBox()
home = os.environ["HOME"]

s.set_mod(mod_name="fleurMaXR3.1")
s.set_mod("bin", home+"/soft/fleurMaXR3.1-nowan/build.GNU")
s.reset("OMP_NUM_THREADS", 1)

s.echo_commands()
