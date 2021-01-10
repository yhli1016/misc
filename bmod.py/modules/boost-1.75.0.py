import os
from sandbox import SandBox


s = SandBox()
home = os.environ["HOME"]

s.set_mod(mod_name="boost-1.75.0")
s.set_mod("pkg", home+"/soft/lib/boost-1.75.0")

s.echo_commands()