import os
from sandbox import SandBox


s = SandBox()
home = os.environ["HOME"]

s.set_mod(mod_name="office")
s.set_mod("bin,lib", home+"/soft/office/doublecmd")
s.set_mod("bin", home+"/soft/office/Typora-linux-x64")
s.set_mod("bin", home+"/soft/office/pycharm-community-2020.3/bin")
s.set_mod("pkg", home+"/soft/office/codeblocks-20.03")

s.echo_commands()