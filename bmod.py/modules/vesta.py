import os
from sandbox import SandBox, get_mod_name


s = SandBox()
home = os.environ["HOME"]

s.set_mod(mod_name=get_mod_name())
s.set_mod("bin,lib", home+"/soft/dft/VESTA-gtk3")
s.set_alias("vesta", "VESTA")

s.echo_commands()
