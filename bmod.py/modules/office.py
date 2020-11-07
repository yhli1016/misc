import os
from sandbox import SandBox


s = SandBox()
home = os.environ["HOME"]

s.set_mod(mod_name="office")
s.set_mod("bin,lib", home+"/soft/doublecmd")
s.set_mod("bin", home+"/soft/Typora-linux-x64")
s.set_mod("bin", home+"/soft/TeXmacs/bin")
s.reset("TEXMACS_PATH", home+"/soft/TeXmacs")

s.echo_commands()