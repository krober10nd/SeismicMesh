import os

from setuptools import find_packages, setup

# https://packaging.python.org/single_source_version/
base_dir = os.path.abspath(os.path.dirname(__file__))
about = {}
with open(os.path.join(base_dir, "SeismicMesh", "__about__.py"), "rb") as f:
    exec(f.read(), about)

setup(
    name="SeismicMesh",
    version=about["__version__"],
    author=about["__author__"],
    author_email=about["__author_email__"],
    packages=find_packages(),
    description="A mesh generator for seismology",
    long_description=open("README.rst").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/krober10nd/SeismicMesh",
    download_url="https://github.com/nschloe/SeismicMesh/releases",
    license=about["__license__"],
    platforms="any",
    install_requires=["numpy", "segyio", "scipy", "meshio", "pybind11"],
    extras_require={"all": ["matplotlib"], "plot": ["matplotlib"]},
    classifiers=[
        about["__status__"],
        about["__license__"],
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Mathematics",
    ],
)
