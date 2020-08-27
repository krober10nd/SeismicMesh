import os
import re
import sys
import platform
import subprocess

from setuptools.command.build_ext import build_ext
from distutils.version import LooseVersion

try:
    from setuptools import setup, Extension, find_packages
except ImportError:
    print("Setuptools is required to build!")

if sys.version_info < (3, 0):
    print("Python 3.0 or higher required, please upgrade.")
    sys.exit(1)


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=""):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(["cmake", "--version"])
        except OSError:
            raise RuntimeError(
                "CMake must be installed to build the following extensions: "
                + ", ".join(e.name for e in self.extensions)
            )

        if platform.system() == "Windows":
            cmake_version = LooseVersion(
                re.search(r"version\s*([\d.]+)", out.decode()).group(1)
            )
            if cmake_version < "3.1.0":
                raise RuntimeError("CMake >= 3.1.0 is required on Windows")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        cmake_args = [
            "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=" + extdir,
            "-DPYTHON_EXECUTABLE=" + sys.executable,
        ]

        cfg = "Debug" if self.debug else "Release"
        build_args = ["--config", cfg]

        if platform.system() == "Windows":
            cmake_args += [
                "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}".format(cfg.upper(), extdir)
            ]
            if sys.maxsize > 2 ** 32:
                cmake_args += ["-A", "x64"]
            build_args += ["--", "/m"]
        else:
            cmake_args += ["-DCMAKE_BUILD_TYPE=" + cfg]
            build_args += ["--", "-j2"]

        env = os.environ.copy()
        env["CXXFLAGS"] = '{} -DVERSION_INFO=\\"{}\\"'.format(
            env.get("CXXFLAGS", ""), self.distribution.get_version()
        )
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        subprocess.check_call(
            ["cmake", ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env
        )
        subprocess.check_call(
            ["cmake", "--build", "."] + build_args, cwd=self.build_temp
        )


with open("README.rst", "r") as fh:
    long_description = fh.read()

setup(
    name="SeismicMesh",
    url="https://github.com/krober10nd/SeismicMesh",
    version="1.2.41",
    author="Keith Roberts",
    author_email="keithrbt0@gmail.com",
    description="2D/3D serial and parallel triangular mesh generation and mesh improvement for seismology",
    long_description=long_description,
    long_description_content_type="text/x-rst",
    setup_requires=["pybind11"],
    install_requires=[
        "numpy",
        "segyio",
        "scipy",
        "meshio",
        "h5py",
        "pybind11",
        "matplotlib",
        "mpi4py",
        "scikit-fmm",
    ],
    packages=find_packages(),
    package_data={"SeismicMesh": ["*.so", "*.txt"]},
    include_package_data=True,
    ext_modules=[
        CMakeExtension("SeismicMesh/sizing/cpp/FastHJ"),
        CMakeExtension("SeismicMesh/generation/cpp/delaunay"),
        CMakeExtension("SeismicMesh/generation/cpp/delaunay_class"),
        CMakeExtension("SeismicMesh/generation/cpp/delaunay_class3"),
        CMakeExtension("SeismicMesh/migration/cpp/cpputils"),
        CMakeExtension("SeismicMesh/geometry/cpp/fast_geometry"),
    ],
    cmdclass=dict(build_ext=CMakeBuild),
    zip_safe=False,
    classifiers=[
        "Development Status :: 4 - Beta",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Mathematics",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Visualization"],
    python_requires='>=3.0',
)
