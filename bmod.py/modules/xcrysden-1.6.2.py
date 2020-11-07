import os
from sandbox import SandBox


s = SandBox()
home = os.environ["HOME"]

s.set_mod(mod_name="xcrysden-1.6.2")
s.set_mod("bin,lib", home+"/soft/xcrysden-1.6.2-bin-shared")

s.echo_commands()
