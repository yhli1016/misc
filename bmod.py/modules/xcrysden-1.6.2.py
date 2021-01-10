import os
from sandbox import SandBox, get_mod_name


s = SandBox()
home = os.environ["HOME"]

s.set_mod(mod_name=get_mod_name())
s.set_mod("bin,lib", home+"/soft/dft/xcrysden-1.6.2-bin-shared")

s.echo_commands()
