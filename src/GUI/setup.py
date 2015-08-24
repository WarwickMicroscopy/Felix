import sys
from cx_Freeze import setup, Executable

# GUI applications require a different base on Windows (the default is for a
# console application).
base = None
if sys.platform == "win32":
    base = "Win32GUI"

setup(  name = "FelixGUI",
        version = "0.1",
        description = "Felixsim GUI Application",
        options = {"build_exe": build_exe_options},
        executables = [Executable("Felix_gui.py", base=base)])
