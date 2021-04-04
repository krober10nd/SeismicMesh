from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup

# https://github.com/pybind/python_example/
is_loc = [
    "SeismicMesh/sizing/cpp/FastHJ",
    "SeismicMesh/generation/cpp/delaunay",
    "SeismicMesh/generation/cpp/delaunay_class",
    "SeismicMesh/generation/cpp/delaunay_class3",
    "SeismicMesh/migration/cpp/cpputils",
    "SeismicMesh/geometry/cpp/fast_geometry",
]

files = [
    "SeismicMesh/sizing/cpp/FastHJ.cpp",
    "SeismicMesh/generation/cpp/delaunay.cpp",
    "SeismicMesh/generation/cpp/delaunay_class.cpp",
    "SeismicMesh/generation/cpp/delaunay_class3.cpp",
    "SeismicMesh/migration/cpp/cpputils.cpp",
    "SeismicMesh/geometry/cpp/fast_geometry.cpp",
]

# no CGAL libraries necessary from CGAL 5.0 onwards
ext_modules = [
    Pybind11Extension(loc, [fi], libraries=["gmp", "mpfr"])
    for fi, loc in zip(files, is_loc)
]

if __name__ == "__main__":
    setup(
        cmdclass={"build_ext": build_ext},
        ext_modules=ext_modules,
        zip_safe=False,
    )
