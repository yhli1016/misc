import os
from sandbox import SandBox


s = SandBox()
home = os.environ["HOME"]

s.set_mod(mod_name="vesta")
s.set_mod("bin,lib", home+"/soft/dft/VESTA-gtk3")
s.set_alias("vesta", "VESTA")

s.echo_commands()
