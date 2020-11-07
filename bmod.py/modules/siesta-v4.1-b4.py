import os
from sandbox import SandBox


s = SandBox()
home = os.environ["HOME"]
siesta = home + "/soft/siesta-v4.1-b4"

s.set_mod(mod_name="siesta-v4.1-b4")
s.set_mod("pkg", siesta+"/Docs/build/flook/0.8.1")
s.set_mod("pkg", siesta+"/Docs/build/hdf5/1.8.21")
s.set_mod("pkg", siesta+"/Docs/build/netcdf/4.7.4")
s.set_mod("pkg", siesta+"/Docs/build/zlib/1.2.11")
s.set_mod("bin", siesta+"/Obj")
s.set_mod("pkg", home+"/soft/elpa-2018.11.001", cmd="rm")
s.set_mod("pkg", home+"/soft/elpa-2019.11.001")
s.reset("OMP_NUM_THREADS", 1)

s.echo_commands()