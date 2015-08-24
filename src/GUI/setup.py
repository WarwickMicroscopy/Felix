import sys
from cx_Freeze import setup, Executable

# Dependencies are automatically detected, but it might need fine tuning.
build_exe_options = {}

# GUI applications require a different base on Windows (the default is for a
# console application).
base = None
if sys.platform == "win32":
    base = "Win32GUI"

setup(  name = "FelixGUI",
        version = "0.1",
        description = "GUI for Felixsim suite of electron diffraction codes",
        options = {"build_exe": build_exe_options},
        executables = [Executable("Felix_gui.py", base=base)])
