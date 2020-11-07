import os
from sandbox import SandBox


s = SandBox()
home = os.environ["HOME"]

s.set_mod(mod_name="openmx3.9")
s.set_mod("pkg", home+"/soft/fftw-3.3.8")
s.set_mod("bin", home+"/soft/openmx3.9/source")
s.reset("OMP_NUM_THREADS", 1)

s.echo_commands()