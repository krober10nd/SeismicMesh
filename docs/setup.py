import os
import re
import sys
import platform
import subprocess

from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
from distutils.version import LooseVersion

setup(
    name="SeismicMesh",
    version="0.5.0",
    author="Keith Roberts",
    author_email="keithrbt0gmail.com",
    description="2D/3D serial and parallel triangular mesh generation for seismology",
    long_description="",
    zip_safe=False,
)
