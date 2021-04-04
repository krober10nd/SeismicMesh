#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cstdarg>
#include <cstddef>
#include <cstring>
#include <iterator>
#include <string>
#include <vector>

#include <boost/lexical_cast.hpp>

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>

namespace py = pybind11;

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Vb = CGAL::Triangulation_vertex_base_with_info_3<unsigned int, K>;
using Tds = CGAL::Triangulation_data_structure_3<Vb>;
using DT = CGAL::Delaunay_triangulation_3<K, Tds>;

using Point = K::Point_3;
using Vertex_handle = DT::Vertex_handle;
using Cell_handle = DT::Cell_handle;
using Edge = DT::Edge;
using Vi = DT::Finite_vertices_iterator;
using Cc = DT::Cell_circulator;
using Ci = DT::Finite_cells_iterator;

template <typename T>
      std::vector<T> vectorSortIntArr(std::vector<std::array<T, 2>> v) {
        std::sort(v.begin(), v.end());
        // double t = tmr.elapsed();
        // tmr.reset();
        auto iter = std::unique(v.begin(), v.end());
        // t = tmr.elapsed();
        // std::cout << t << std::endl;

        size_t len = iter - v.begin();
        std::vector<T> outvec;
        outvec.reserve(len * 2);
        for (auto i = v.begin(); i != iter; ++i) {
          outvec.push_back(i->at(0));
          outvec.push_back(i->at(1));
        }
        return outvec;
      }



template <typename T> class TypedInputIterator {
public:
  using iterator_category = std::input_iterator_tag;
  using difference_type = std::ptrdiff_t;
  using value_type = T;
  using pointer = T *;
  using reference = T &;

  explicit TypedInputIterator(py::iterator &py_iter) : py_iter_(py_iter) {}

  explicit TypedInputIterator(py::iterator &&py_iter) : py_iter_(py_iter) {}

  value_type operator*() { return (*py_iter_).template cast<value_type>(); }

  TypedInputIterator operator++(int) {
    auto copy = *this;
    ++py_iter_;
    return copy;
  }

  TypedInputIterator &operator++() {
    ++py_iter_;
    return *this;
  }

  bool operator!=(TypedInputIterator &rhs) { return py_iter_ != rhs.py_iter_; }

  bool operator==(TypedInputIterator &rhs) { return py_iter_ == rhs.py_iter_; }

private:
  py::iterator py_iter_;
};

PYBIND11_MODULE(_delaunay_class3, m) {
  py::class_<Point>(m, "Point")
      .def(py::init<int, int, int>(), py::arg("x"), py::arg("y"), py::arg("z"))
      .def(py::init<double, double, double>(), py::arg("x"), py::arg("y"),
           py::arg("z"))
      .def_property_readonly("x", &Point::x)
      .def_property_readonly("y", &Point::y)
      .def_property_readonly("z", &Point::z)
      .def("__repr__", [](const Point &p) {
        std::string r("Point(");
        r += boost::lexical_cast<std::string>(p.x());
        r += ", ";
        r += boost::lexical_cast<std::string>(p.y());
        r += ", ";
        r += boost::lexical_cast<std::string>(p.z());
        r += ")";
        return r;
      });

  py::class_<Vertex_handle>(m, "VertexHandle")
      .def_property_readonly("point", [](const Vertex_handle &vertex_handle) {
        return vertex_handle->point();
      });

  py::class_<DT>(m, "DelaunayTriangulation3")

      .def(py::init())

      .def("insert",
           [](DT &dt, const std::vector<double> &p) {
             std::vector<std::pair<Point, unsigned>> points;
             int num_points = p.size() / 3;
             // start adding at the end of the current table
             int start = dt.number_of_vertices();
             for (std::size_t i = 0; i < num_points; ++i) {
               // add index information to form face table later
               points.push_back(std::make_pair(
                   Point(p[i * 3 + 0], p[i * 3 + 1], p[i * 3 + 2]), start));
               start += 1;
             }
             return dt.insert(points.begin(), points.end());
           })

      .def("move",
           [](DT &dt, const std::vector<unsigned int> &to_move,
              const std::vector<double> &new_positions) {
             std::vector<Vertex_handle> handles;
             std::vector<Point> new_pos;
             int num_to_move = to_move.size();
             // store all vertex handles
             // Should this be finite_vertices_begin and finite_vertices_end?
             for (Vi vi = dt.finite_vertices_begin();
                  vi != dt.finite_vertices_end(); vi++) {
               handles.push_back(vi);
             }
             // store new positions as a vector of Point
             for (std::size_t i = 0; i < num_to_move; ++i) {
               new_pos.push_back(Point(new_positions[3 * i],
                                       new_positions[3 * i + 1],
                                       new_positions[3 * i + 2]));
             }
             //
             for (std::size_t i = 0; i < num_to_move; ++i) {
               dt.move(handles[to_move[i]], new_pos[i]);
             }
             return dt;
           })

      .def("remove",
           [](DT &dt, const std::vector<unsigned int> &to_remove) {
             int num_to_remove = to_remove.size();
             std::vector<Vertex_handle> handles;
             for (Vi vi = dt.finite_vertices_begin();
                  vi != dt.finite_vertices_end(); vi++) {
               handles.push_back(vi);
             }
             for (std::size_t i = 0; i < num_to_remove; ++i) {
               dt.remove(handles[to_remove[i]]);
             }
             return dt;
           })

      .def("number_of_vertices", [](DT &dt) { return dt.number_of_vertices(); })

      .def("get_edges",[](const DT &dt, const std::vector<int> &cells_to_get){

              // User supplies cells_to_get
              std::unordered_set<Cell_handle> cell_handles;
              cell_handles.reserve(cells_to_get.size());
              int ix = 0;
              for (Ci ci = dt.finite_cells_begin(); ci != dt.finite_cells_end(); ci++) {
                  if(cells_to_get[ix] == 1){
                      cell_handles.insert(ci);
                    }
                  ix += 1;
              }
              std::vector<int> edges;
              edges.reserve(6*cells_to_get.size());
              // this creates each edge of the tetrahedron
              int indices[6][2];
              // first edge
              indices[0][0]=0;
              indices[0][1]=1;
              // second edge
              indices[1][0]=1;
              indices[1][1]=2;
              // third edge
              indices[2][0]=2;
              indices[2][1]=0;
              // fourth edge
              indices[3][0]=0;
              indices[3][1]=3;
              // fifth edge
              indices[4][0]=1;
              indices[4][1]=3;
              // sixth edge
              indices[5][0]=2;
              indices[5][1]=3;
              // If it is in the set report the unique edge by comparing the &* of the two face handle
              int n = 0;
              for (const auto& cell : cell_handles){
                  // for each cell look at its six edges
                  for (std::size_t i = 0; i < 6; ++i){
                     // build edge structure
                     int indx0 = indices[i][0];
                     int indx1 = indices[i][1];
                     auto edge = Edge(cell, indx0, indx1);
                     // circulate around the edge
                     Cc nei_cells = dt.incident_cells(edge, cell);
                     Cc done = nei_cells;
                     bool is_smallest = true;
                     do
                     {
                        Cell_handle nei_cell = nei_cells;
                        if (cell > nei_cell && cell_handles.count(nei_cell)) {
                           ++nei_cells; // advance to next neighoring cell
                           is_smallest = false;
                           break; // report interior edges once
                        }
                        ++nei_cells;
                     }
                     while(nei_cells != done);

                     if (is_smallest)
                     {
                       Vertex_handle vs = cell->vertex(indx0);
                       Vertex_handle vt = cell->vertex(indx1);

                       edges.push_back(vs->info());
                       edges.push_back(vt->info());
                       n += 2;
                     }
                  }
              }
           edges.resize(n);
           ssize_t soint = sizeof(int);
           ssize_t num_edges = edges.size() / 2;
           ssize_t ndim = 2;
           std::vector<ssize_t> shape = {num_edges, 2};
           std::vector<ssize_t> strides = {soint * 2, soint};

           // return 2-D NumPy array
           return py::array(py::buffer_info(
               edges.data(), /* data as contiguous array  */
               sizeof(int),  /* size of one scalar        */
               py::format_descriptor<int>::format(), /* data type */
               2,      /* number of dimensions      */
               shape,  /* shape of the matrix       */
               strides /* strides for each axis     */
               ));
           })

      .def("number_of_cells",
           [](DT &dt) {
             int count = 0;
             for (DT::Finite_cells_iterator fit = dt.finite_cells_begin();
                  fit != dt.finite_cells_end(); ++fit) {
               count += 1;
             }
             return count;
           })

      .def("finite_vertices",
           [](DT &dt) -> py::iterator {
             return py::make_iterator(dt.finite_vertices_begin(),
                                      dt.finite_vertices_end());
           })

      .def("get_finite_cells",
           [](DT &dt) {
             // ouput the cell table
             // YOU MUST CALL get_finite_vertices before if any incremental
             // operations were performed
             std::vector<int> cells;
             cells.resize(dt.number_of_finite_cells() * 4);

             int i = 0;
             for (DT::Finite_cells_iterator fit = dt.finite_cells_begin();
                  fit != dt.finite_cells_end(); ++fit) {

               DT::Cell_handle cell = fit;
               cells[i * 4] = cell->vertex(0)->info();
               cells[i * 4 + 1] = cell->vertex(1)->info();
               cells[i * 4 + 2] = cell->vertex(2)->info();
               cells[i * 4 + 3] = cell->vertex(3)->info();
               i += 1;
             }
             ssize_t soint = sizeof(int);
             ssize_t num_cells = cells.size() / 4;
             ssize_t ndim = 2;
             std::vector<ssize_t> shape = {num_cells, 4};
             std::vector<ssize_t> strides = {soint * 4, soint};

             // return 2-D NumPy array
             return py::array(py::buffer_info(
                 cells.data(), /* data as contiguous array  */
                 sizeof(int),  /* size of one scalar        */
                 py::format_descriptor<int>::format(), /* data type */
                 2,      /* number of dimensions      */
                 shape,  /* shape of the matrix       */
                 strides /* strides for each axis     */
                 ));
           })

      .def("get_finite_vertices", [](DT &dt) {
        // ouput the vertices
        std::vector<double> vertices;
        vertices.resize(dt.number_of_vertices() * 3);

        int i = 0;
        for (DT::Finite_vertices_iterator fit = dt.finite_vertices_begin();
             fit != dt.finite_vertices_end(); ++fit) {

          Vertex_handle vertex = fit;
          // critical! update the point index table so faces comes out correctly
          vertex->info() = i;
          vertices[i * 3] = vertex->point().x();
          vertices[i * 3 + 1] = vertex->point().y();
          vertices[i * 3 + 2] = vertex->point().z();
          i += 1;
        }
        ssize_t sdble = sizeof(double);
        ssize_t num_vertices = vertices.size() / 3;
        ssize_t ndim = 2;
        std::vector<ssize_t> shape = {num_vertices, 3};
        std::vector<ssize_t> strides = {sdble * 3, sdble};

        // return 2-D NumPy array
        return py::array(py::buffer_info(
            vertices.data(), /* data as contiguous array  */
            sizeof(double),  /* size of one scalar        */
            py::format_descriptor<double>::format(), /* data type */
            2,      /* number of dimensions      */
            shape,  /* shape of the matrix       */
            strides /* strides for each axis     */
            ));
      });

  py::class_<DT::Finite_vertices_iterator::value_type>(m, "Vertex")
      .def_property_readonly(
          "point", [](DT::Finite_vertices_iterator::value_type &vertex) {
            return vertex.point();
          });

  py::class_<DT::Finite_cells_iterator::value_type>(m, "Cell").def(
      "vertex_handle",
      [](DT::Finite_cells_iterator::value_type &cell, int index) {
        return cell.vertex(index);
      },
      py::arg("index"));
}
