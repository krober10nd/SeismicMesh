from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup

# https://github.com/pybind/python_example/
ext_modules = [
    Pybind11Extension(
        "_SeismicMesh",
        # Sort input source files to ensure bit-for-bit reproducible builds
        # (https://github.com/pybind/python_example/pull/53)
        sorted(
            [
                "SeismicMesh/sizing/cpp/FastHJ.cpp",
                "SeismicMesh/generation/cpp/delaunay.cpp",
                "SeismicMesh/generation/cpp/delaunay_class.cpp",
                "SeismicMesh/generation/cpp/delaunay_class3.cpp",
                "SeismicMesh/migration/cpp/cpputils.cpp",
                "SeismicMesh/geometry/cpp/fast_geometry.cpp",
            ]
        ),
        # no CGAL libraries necessary from CGAL 5.0 onwards
        libraries=["gmp", "mpfr"],
    )
]

if __name__ == "__main__":
    setup(
        cmdclass={"build_ext": build_ext},
        ext_modules=ext_modules,
        zip_safe=False,
    )
