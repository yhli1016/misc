About
-----
Pmod is a tiny module system written in Python. Initially it came out as the
successor of a complicated shell script for managing user-dependent
environmental variables. As shell script was not well suited for such tasks,
and we did not have enough time to learn and deploy advanced module systems
such as lmod, we decided to write our own module system. Installation of pmod
is very easy. Root privileges are not essential. Configuration of pmod requires
some basic python knowledge. Basic capabilities of a common module system are
supported, such as listing available modules, recursively loading and unloading
modules, and printing the status of modules.


Install
-------
Pmod does not have other requisties besides Python. Both Python 2 and Python 3
are supported. The common installation procedure is as below:

1. Put init.sh, modcmd.py and modconfig.py somewhere, e.g. ~/soft/pmod. If you
   plan to install pmod for all the users, we suggest /opt/pmod as the
   installation destination. Note that in this case root privileges are
   essential.

2. Edit init.sh to update PMODPATH to the installation destination of pmod in
   step 1.

3. Add the command to source init.sh in your ~/.bashrc. If pmod is to be
   installed for all the users, copy init.sh to /etc/profile.d.


Configuration
-------------
Configurations of available modules are stored in modconfig.py as a global
dictionary named avail\_mods, with the names of modules as the keys. Each module
is further presented as a dictionary with the keys being the environmental
variables to be modified. The auxiliary function add\_mod() is provided for the
manipulation of avail\_mods.

Modules can have dependencies, which are stored in a list with the key
"DEPENDENCY". The auxiliary function check\_integrity() checks if there are any
unresolved dependencies. It is strongly recommended that you call this
functional at the end of modconfig.py.

Modules can also have additional initialization commands. These commands are
stored in a list with the key "CMD" and will be executed BEFORE setting the
environmental variables.

See modconfig.py for more details and examples.


Usage
-----
The common usage of pmod is "pmod [-a] [-r] operation [mod\_name]". "-a" tells
pmod to load all the modules. "-r" tells pmod to recursively load or unload the
specified module as well as the dependencies, while omitting "-r" loads or
unloads the module only. Operations are:

- avail: listing available modules
- load, add: loading the specified module
- unload, remove, rm: unloading the specified module
- list, ls, status, stat: printing the status of modules, with [L] indicating
  that the module has been loaded, [U] for unloaded modules and [B] for loaded
  but broken modules

See modcmd.py for more details.
