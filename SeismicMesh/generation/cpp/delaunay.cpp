#include <boost/lexical_cast.hpp>
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <assert.h>
#include <vector>

#include <CGAL/Kernel/global_functions.h>

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned int, Kernel> Vb;
typedef CGAL::Triangulation_data_structure_2<Vb> Tds;
typedef CGAL::Delaunay_triangulation_2<Kernel, Tds> Delaunay;
typedef Kernel::Point_2 Point;

typedef CGAL::Triangulation_vertex_base_with_info_3<unsigned int, Kernel> Vb3;
typedef CGAL::Triangulation_data_structure_3<Vb3> Tds3;
typedef CGAL::Delaunay_triangulation_3<Kernel, Tds3> Delaunay3;
typedef Kernel::Point_3 Point3;

std::vector<double> c_circumballs2(std::vector<double> &vertices) {
  int num_faces = vertices.size() / 6;
  std::vector<double> circumcenters;
  for (std::size_t i = 0; i < num_faces; ++i) {
    Point tmp_cc =
        CGAL::circumcenter(Point(vertices[i * 6], vertices[i * 6 + 1]),
                           Point(vertices[i * 6 + 2], vertices[i * 6 + 3]),
                           Point(vertices[i * 6 + 4], vertices[i * 6 + 5]));
    circumcenters.push_back(tmp_cc.x());
    circumcenters.push_back(tmp_cc.y());
    circumcenters.push_back(
        CGAL::squared_radius(Point(vertices[i * 6], vertices[i * 6 + 1]),
                             Point(vertices[i * 6 + 2], vertices[i * 6 + 3]),
                             Point(vertices[i * 6 + 4], vertices[i * 6 + 5])));
  }
  return circumcenters;
}

std::vector<double> c_circumballs3(std::vector<double> &vertices) {
  int num_cells = vertices.size() / 12;
  std::vector<double> circumcenters;
  for (std::size_t i = 0; i < num_cells; ++i) {
    Point3 tmp_cc = CGAL::circumcenter(
        Point3(vertices[i * 12], vertices[i * 12 + 1], vertices[i * 12 + 2]),
        Point3(vertices[i * 12 + 3], vertices[i * 12 + 4],
               vertices[i * 12 + 5]),
        Point3(vertices[i * 12 + 6], vertices[i * 12 + 7],
               vertices[i * 12 + 8]),
        Point3(vertices[i * 12 + 9], vertices[i * 12 + 10],
               vertices[i * 12 + 11]));
    circumcenters.push_back(tmp_cc.x());
    circumcenters.push_back(tmp_cc.y());
    circumcenters.push_back(tmp_cc.z());
    circumcenters.push_back(CGAL::squared_radius(
        Point3(vertices[i * 12], vertices[i * 12 + 1], vertices[i * 12 + 2]),
        Point3(vertices[i * 12 + 3], vertices[i * 12 + 4],
               vertices[i * 12 + 5]),
        Point3(vertices[i * 12 + 6], vertices[i * 12 + 7],
               vertices[i * 12 + 8]),
        Point3(vertices[i * 12 + 9], vertices[i * 12 + 10],
               vertices[i * 12 + 11])));
  }
  return circumcenters;
}

std::vector<int> c_delaunay2(std::vector<double> &x, std::vector<double> &y) {
  int num_points = x.size();
  assert(y.size() != num_points);
  std::vector<std::pair<Point, unsigned>> points;
  // add index information to form face table later
  for (std::size_t i = 0; i < num_points; ++i) {
    points.push_back(std::make_pair(Point(x[i], y[i]), i));
  }

  Delaunay triangulation;
  triangulation.insert(points.begin(), points.end());

  // save the face table
  int num_faces = triangulation.number_of_faces();
  std::vector<int> faces;
  faces.resize(num_faces * 3);

  int i = 0;
  for (Delaunay::Finite_faces_iterator fit = triangulation.finite_faces_begin();
       fit != triangulation.finite_faces_end(); ++fit) {

    Delaunay::Face_handle face = fit;
    faces[i * 3] = face->vertex(0)->info();
    faces[i * 3 + 1] = face->vertex(1)->info();
    faces[i * 3 + 2] = face->vertex(2)->info();
    i += 1;
  }
  return faces;
}

std::vector<int> c_delaunay3(std::vector<double> &x, std::vector<double> &y,
                             std::vector<double> &z) {
  int num_points = x.size();
  assert(y.size() != num_points);
  assert(z.size() != num_points);
  std::vector<std::pair<Point3, unsigned>> points;
  // add index information to form face table later
  for (std::size_t i = 0; i < num_points; ++i) {
    points.push_back(std::make_pair(Point3(x[i], y[i], z[i]), i));
  }
  Delaunay3 triangulation;
  triangulation.insert(points.begin(), points.end());
  // save the indices of all cells
  int num_cells = triangulation.number_of_finite_cells();
  std::vector<int> cells;
  cells.resize(num_cells * 4);

  int i = 0;
  for (Delaunay3::Finite_cells_iterator cit =
           triangulation.finite_cells_begin();
       cit != triangulation.finite_cells_end(); ++cit) {

    Delaunay3::Cell_handle cell = cit;
    cells[i * 4] = cell->vertex(0)->info();
    cells[i * 4 + 1] = cell->vertex(1)->info();
    cells[i * 4 + 2] = cell->vertex(2)->info();
    cells[i * 4 + 3] = cell->vertex(3)->info();
    i += 1;
  }
  return cells;
}

// ----------------
// Python interface
// ----------------
// (from
// https://github.com/tdegeus/pybind11_examples/blob/master/04_numpy-2D_cpp-vector/example.cpp)

namespace py = pybind11;
py::array circumballs2(
    py::array_t<double, py::array::c_style | py::array::forcecast> vertices) {
  // each triangle has 3 vertices with 2 coordinates each
  int sz = vertices.shape()[0];
  std::vector<double> cppvertices(sz);
  std::memcpy(cppvertices.data(), vertices.data(), sz * sizeof(double));
  std::vector<double> circumcenters = c_circumballs2(cppvertices);
  ssize_t soreal = sizeof(double);
  ssize_t num_points = circumcenters.size() / 3;
  ssize_t ndim = 2;
  std::vector<ssize_t> shape = {num_points, 3};
  std::vector<ssize_t> strides = {soreal * 3, soreal};
  // return 2-D NumPy array
  return py::array(
      py::buffer_info(circumcenters.data(), /* data as contiguous array  */
                      sizeof(double),       /* size of one scalar        */
                      py::format_descriptor<double>::format(), /* data type */
                      2,      /* number of dimensions      */
                      shape,  /* shape of the matrix       */
                      strides /* strides for each axis     */
                      ));
}

py::array circumballs3(
    py::array_t<double, py::array::c_style | py::array::forcecast> vertices) {
  // each triangle has 4 vertices with 3 coordinates each
  int sz = vertices.size();
  std::vector<double> cppvertices(sz);
  std::memcpy(cppvertices.data(), vertices.data(), sz * sizeof(double));
  std::vector<double> circumcenters = c_circumballs3(cppvertices);
  ssize_t soreal = sizeof(double);
  ssize_t num_points = circumcenters.size() / 4;
  ssize_t ndim = 2;
  std::vector<ssize_t> shape = {num_points, 4};
  std::vector<ssize_t> strides = {soreal * 4, soreal};
  // return 2-D NumPy array
  return py::array(
      py::buffer_info(circumcenters.data(), /* data as contiguous array  */
                      sizeof(double),       /* size of one scalar        */
                      py::format_descriptor<double>::format(), /* data type */
                      2,      /* number of dimensions      */
                      shape,  /* shape of the matrix       */
                      strides /* strides for each axis     */
                      ));
}

py::array
delaunay2(py::array_t<double, py::array::c_style | py::array::forcecast> x,
          py::array_t<double, py::array::c_style | py::array::forcecast> y) {

  // check input dimensions
  if (x.ndim() != 1)
    throw std::runtime_error("Input should be 2 1D NumPy arrays");
  if (y.ndim() != 1)
    throw std::runtime_error("Input should be 2 1D NumPy arrays");

  int num_points = x.shape()[0];

  // allocate std::vector (to pass to the C++ function)
  std::vector<double> cppx(num_points);
  std::vector<double> cppy(num_points);

  // copy py::array -> std::vector
  std::memcpy(cppx.data(), x.data(), num_points * sizeof(double));
  std::memcpy(cppy.data(), y.data(), num_points * sizeof(double));
  std::vector<int> faces = c_delaunay2(cppx, cppy);

  ssize_t soint = sizeof(int);
  ssize_t num_faces = faces.size() / 3;
  ssize_t ndim = 2;
  std::vector<ssize_t> shape = {num_faces, 3};
  std::vector<ssize_t> strides = {soint * 3, soint};

  // return 2-D NumPy array
  return py::array(
      py::buffer_info(faces.data(), /* data as contiguous array  */
                      sizeof(int),  /* size of one scalar        */
                      py::format_descriptor<int>::format(), /* data type */
                      2,      /* number of dimensions      */
                      shape,  /* shape of the matrix       */
                      strides /* strides for each axis     */
                      ));
}

py::array
delaunay3(py::array_t<double, py::array::c_style | py::array::forcecast> x,
          py::array_t<double, py::array::c_style | py::array::forcecast> y,
          py::array_t<double, py::array::c_style | py::array::forcecast> z) {

  // check input dimensions
  if (x.ndim() != 1)
    throw std::runtime_error("Input should be three 1D NumPy arrays");
  if (y.ndim() != 1)
    throw std::runtime_error("Input should be three 1D NumPy arrays");
  if (z.ndim() != 1)
    throw std::runtime_error("Input should be three 1D NumPy arrays");

  int num_points = x.shape()[0];

  // allocate std::vector (to pass to the C++ function)
  std::vector<double> cppx(num_points);
  std::vector<double> cppy(num_points);
  std::vector<double> cppz(num_points);

  // copy py::array -> std::vector
  std::memcpy(cppx.data(), x.data(), num_points * sizeof(double));
  std::memcpy(cppy.data(), y.data(), num_points * sizeof(double));
  std::memcpy(cppz.data(), z.data(), num_points * sizeof(double));
  std::vector<int> cells = c_delaunay3(cppx, cppy, cppz);

  ssize_t num_cells = cells.size() / 4;
  ssize_t ndim = 2;
  ssize_t soint = sizeof(int);
  std::vector<ssize_t> shape = {num_cells, 4};
  std::vector<ssize_t> strides = {soint * 4, soint};

  // return 2-D NumPy array
  return py::array(
      py::buffer_info(cells.data(), /* data as contiguous array  */
                      sizeof(int),  /* size of one scalar        */
                      py::format_descriptor<int>::format(), /* data type */
                      2,      /* number of dimensions      */
                      shape,  /* shape of the matrix       */
                      strides /* strides for each axis     */
                      ));
}

PYBIND11_MODULE(c_cgal, m) {
  m.def("circumballs3", &circumballs3);
  m.def("circumballs2", &circumballs2);
  m.def("delaunay2", &delaunay2);
  m.def("delaunay3", &delaunay3);
}
