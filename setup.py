"""
setup.py -- Builds the C++ WaveQTL engine via `make all` automatically during `pip install .`.
"""

from pathlib import Path
import subprocess
import sys

from setuptools import setup
from setuptools.command.build_py import build_py

try:
    from setuptools.command.editable_wheel import editable_wheel
except ImportError:
    editable_wheel = None


def _build_c_engine() -> None:
    repo_root = Path(__file__).resolve().parent
    try:
        subprocess.check_call(["make", "clean"], cwd=repo_root)
        subprocess.check_call(["make", "all"], cwd=repo_root)
    except (subprocess.CalledProcessError, FileNotFoundError) as exc:
        sys.stderr.write(
            "\n"
            "========================================================================\n"
            "ERROR: Compilation of the C++ WaveQTL engine failed.\n"
            "Please ensure 'g++', 'make', and the GSL library headers are installed:\n"
            "    sudo apt-get update && sudo apt-get install -y build-essential libgsl-dev\n"
            "========================================================================\n\n"
        )
        raise RuntimeError("WaveQTL C++ engine build failed.") from exc


class BuildWithMake(build_py):
    def run(self):
        _build_c_engine()
        super().run()


cmdclass = {"build_py": BuildWithMake}

if editable_wheel is not None:
    class EditableWithMake(editable_wheel):
        def run(self):
            _build_c_engine()
            super().run()

    cmdclass["editable_wheel"] = EditableWithMake

setup(cmdclass=cmdclass)
