import os
import sys


def print_stdout(command):
    """
    Print commands to stdout, which are then interpreted by shell.

    :param command: string, command to be interpreted by shell
    :return: None
    """
    sys.stdout.write("%s\n" % command)
    sys.stdout.flush()


def print_stderr(message):
    """
    Print message to stderr, which WILL NOT be interpreted by shell.
    This function is also utilized to print banners and tables.

    :param message: string, message to write to stderr
    :return: None
    """
    sys.stderr.write("%s\n" % message)
    sys.stderr.flush()


def get_mod_name():
    """Get module name from sys.argv"""
    return sys.argv[0].split("/")[-1].split(".py")[0]


class SandBox(object):
    """
    Class that records the settings of each module and outputs shell commands
    to stdout.

    self.environ is a copy of os.environ, but has all the elements split by ":".

    self.env_name_changed records the environmental variables modified by the
    modules to load or unload.

    self.command, self.alias and self.unalias records the command and alias of
    each module.
    """
    def __init__(self):
        self.environ = dict()
        for env_name, env_value in os.environ.items():
            self.environ[env_name] = env_value.split(":")
        self.env_name_changed = []
        self.command = []
        self.alias = []
        self.unalias = []

    def match_pattern(self, env_name, pattern):
        """
        Check if pattern equals self.environ[env_name].

        :param env_name: string, name of the environmental variable
        :param pattern: string, pattern to check
        :return: True if self.environ[env_name] equals pattern, False otherwise.
        """
        if env_name not in self.environ.keys():
            return False
        else:
            return pattern == self.environ[env_name][0]

    def has_pattern(self, env_name, pattern):
        """
        Check if pattern is already included in self.environ[env_name].

        :param env_name: string, name of the environmental variable
        :param pattern: string, pattern to check
        :return: True if pattern already included, False otherwise.
        """
        if env_name not in self.environ.keys():
            return False
        else:
            return pattern in self.environ[env_name]

    def reset_env(self, env_name, pattern):
        """
        Reset the value of environmental variable.

        :param env_name: string, name of the environmental variable
        :param pattern: string, new value of environmental variable
        :return: None
        """
        if not self.match_pattern(env_name, pattern):
            self.environ[env_name] = [pattern]
            if env_name not in self.env_name_changed:
                self.env_name_changed.append(env_name)

    def append_env(self, env_name, pattern):
        """
        Append pattern to environmental variable.

        :param env_name: string, name of the environmental variable
        :param pattern: string, with which the environmental variable will be
                        modified
        :return: None
        """
        if not self.has_pattern(env_name, pattern):
            if env_name not in self.environ.keys():
                self.environ[env_name] = [pattern]
            else:
                self.environ[env_name].append(pattern)
            if env_name not in self.env_name_changed:
                self.env_name_changed.append(env_name)

    def prepend_env(self, env_name, pattern):
        """
        Prepend pattern to environmental variable.

        :param env_name: string, name of the environmental variable
        :param pattern: string, with which the environmental variable will be
                        modified
        :return: None
        """
        if not self.has_pattern(env_name, pattern):
            if env_name not in self.environ.keys():
                self.environ[env_name] = [pattern]
            else:
                self.environ[env_name].insert(0, pattern)
            if env_name not in self.env_name_changed:
                self.env_name_changed.append(env_name)

    def remove_env(self, env_name, pattern):
        """
        Remove pattern from environmental variable.

        :param env_name: string, name of the environmental variable
        :param pattern: string, pattern to be removed from the environmental
                        variable
        :return: None
        """
        if self.has_pattern(env_name, pattern):
            while pattern in self.environ[env_name]:
                self.environ[env_name].remove(pattern)
            if env_name not in self.env_name_changed:
                self.env_name_changed.append(env_name)

    def add_alias(self, alias):
        """
        Add aliases to be set to self.alias.

        :param alias: tuple of (alias, full_name)
        :return: None.
        """
        self.alias.append(alias)

    def add_unalias(self, unalias):
        """
        Add unaliases to self.unalias.

        :param unalias: tuple of (alias, full_name)
        :return: None
        """
        self.unalias.append(unalias)

    def add_command(self, command):
        """
        Add command to self.command.

        :param command: list of commands, similar to that of Module class.
        :return: None
        """
        self.command.extend(command)

    def echo_commands(self, shell="bash"):
        """
        Print shell commands to stdout to be evaluated by shell.

        :param shell: string, type of the shell that evaluates the output.
        :return: None
        """
        for env_name in self.env_name_changed:
            env_value = self.environ[env_name]
            env_string = "".join(["%s:" % pattern for pattern in env_value])
            while env_string != "" and env_string[-1] == ":":
                env_string = env_string[:-1]
            if shell == "bash":
                print_stdout("export %s=%s;" % (env_name, env_string))
            else:
                raise NotImplementedError("Shell type %s not supported" % shell)
        if shell == "bash":
            for alias in self.unalias:
                print_stdout("unalias %s;" % alias[0])
            for alias in self.alias:
                print_stdout("alias %s=\"%s\";" % (alias[0], alias[1]))
        for command in self.command:
            print_stdout("%s;" % command)

    def reset(self, env_name=None, pattern=None, cmd="auto"):
        """Wrapper for self.reset_env."""
        if cmd == "auto":
            cmd = sys.argv[1]
        if cmd == "add":
            self.reset_env(env_name, pattern)
        elif cmd == "rm":
            self.reset_env(env_name, "")
        else:
            print_stderr("ERROR: unknown command '%s'" % cmd)
            sys.exit(-1)

    def set_env(self, env_name=None, pattern=None, cmd="auto"):
        "Wrapper for self.prepend_env and self.remove_env."
        if cmd == "auto":
            cmd = sys.argv[1]
        if cmd == "add":
            self.prepend_env(env_name, pattern)
        elif cmd == "rm":
            self.remove_env(env_name, pattern)
        else:
            print_stderr("ERROR: unknown command '%s'" % cmd)
            sys.exit(-1)

    def set_mod(self, presets=None, dest=None, mod_name=None, cmd="auto"):
        """
        Wrapper for self.prepend_env and self.remove_env with pre-defined presets.
        """
        if cmd == "auto":
            cmd = sys.argv[1]
        if cmd in ("add", "rm"):
            func = self.prepend_env if cmd == "add" else self.remove_env
            if presets is not None:
                for preset in presets.split(","):
                    if preset == "pkg":
                        if os.path.exists("%s/bin" % dest):
                            func("PATH", "%s/bin" % dest)
                        for path in ("lib", "lib64"):
                            lib_path = "%s/%s" % (dest, path)
                            if os.path.exists(lib_path):
                                func("LIBRARY_PATH", lib_path)
                                func("LD_RUN_PATH", lib_path)
                                func("LD_LIBRARY_PATH", lib_path)
                                if os.path.exists("%s/pkgconfig" % lib_path):
                                    func("PKG_CONFIG_PATH", "%s/pkgconfig" % lib_path)
                        if os.path.exists("%s/include" % dest):
                            func("C_INCLUDE_PATH", "%s/include" % dest)
                            func("CPLUS_INCLUDE_PATH", "%s/include" % dest)
                    elif preset == "bin":
                        func("PATH", dest)
                    elif preset == "lib":
                        func("LIBRARY_PATH", dest)
                        func("LD_RUN_PATH", dest)
                        func("LD_LIBRARY_PATH", dest)
                        if os.path.exists("%s/pkgconfig" % dest):
                            func("PKG_CONFIG_PATH", "%s/pkgconfig" % dest)
                    elif preset == "inc":
                        func("C_INCLUDE_PATH", dest)
                        func("CPLUS_INCLUDE_PATH", dest)
                    elif preset == "py":
                        func("PYTHONPATH", dest)
                    else:
                        print_stderr("ERROR: undefined preset type %s" % preset)
                        sys.exit(-1)
        else:
            print_stderr("ERROR: unknown command '%s'" % cmd)
            sys.exit(-1)
        if mod_name is not None:
            func("BMOD_LOADED_MODS", mod_name)

    def set_alias(self, alias, full_name, cmd="auto"):
        "Wrapper for self.add_alias and self.add_unalias."
        if cmd == "auto":
            cmd = sys.argv[1]
        if cmd == "add":
            self.add_alias((alias, full_name))
        elif cmd == "rm":
            self.add_unalias((alias, full_name))
        else:
            print_stderr("ERROR: unknown command '%s'" % cmd)
            sys.exit(-1)


if __name__ == "__main__":
    sandbox = SandBox()
    if sys.argv[1] == "reset":
        sandbox.reset(cmd=sys.argv[2], env_name=sys.argv[3], pattern=sys.argv[4])
    elif sys.argv[1] == "set_env":
        sandbox.set_env(cmd=sys.argv[2], env_name=sys.argv[3], pattern=sys.argv[4])
    elif sys.argv[1] == "set_mod":
        sandbox.set_mod(cmd=sys.argv[2], presets=sys.argv[3], dest=sys.argv[4])
    else:
        print_stderr("ERROR: unknown function '%s'" % sys.argv[1])
        sys.exit(-1)
    sandbox.echo_commands()
