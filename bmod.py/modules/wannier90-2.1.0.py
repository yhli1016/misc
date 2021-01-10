import os
from sandbox import SandBox


s = SandBox()
home = os.environ["HOME"]

s.set_mod(mod_name="wannier90-2.1.0")
s.set_mod("bin,lib", home+"/soft/dft/wannier90-2.1.0")

s.echo_commands()
