import os
from sandbox import SandBox, get_mod_name


s = SandBox()
home = os.environ["HOME"]

s.set_mod(mod_name=get_mod_name())
s.set_mod("pkg", home+"/soft/lib/boost-1.75.0")

s.echo_commands()