import os
from sandbox import SandBox


s = SandBox()
home = os.environ["HOME"]

s.set_mod(mod_name="wannier90-2.1.0")
s.set_mod("bin", home+"/soft/fleurMaXR3.1-nowan/build.GNU")
s.set_mod("bin,lib", home+"/soft/wannier90-2.1.0")
s.set_mod("bin", home+"/soft/spex05.00/bin")
s.reset("OMP_NUM_THREADS", 1)

s.echo_commands()
