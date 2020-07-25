#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/numpy.h>

#include <vector>
#include <string>
#include <cstdarg>
#include <cstring>
#include <cstddef>
#include <iterator>

#include <boost/lexical_cast.hpp>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

namespace py = pybind11;

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Vb = CGAL::Triangulation_vertex_base_with_info_2<unsigned int, K>;
using Tds = CGAL::Triangulation_data_structure_2<Vb>;
using DT = CGAL::Delaunay_triangulation_2<K, Tds>;

using Point = K::Point_2;
using Vertex_handle = DT::Vertex_handle;
using Vi = DT::Finite_vertices_iterator;


template <typename T>
class TypedInputIterator
{
public:
    using iterator_category = std::input_iterator_tag;
    using difference_type = std::ptrdiff_t;
    using value_type = T;
    using pointer = T*;
    using reference = T&;

    explicit TypedInputIterator(py::iterator& py_iter) :
        py_iter_(py_iter)
    {
    }

    explicit TypedInputIterator(py::iterator&& py_iter) :
        py_iter_(py_iter)
    {
    }

    value_type operator*()
    {
        return (*py_iter_).template cast<value_type>();
    }

    TypedInputIterator operator++(int)
    {
        auto copy = *this;
        ++py_iter_;
        return copy;
    }

    TypedInputIterator& operator++()
    {
        ++py_iter_;
        return *this;
    }

    bool operator!=(TypedInputIterator &rhs)
    {
        return py_iter_ != rhs.py_iter_;
    }

    bool operator==(TypedInputIterator &rhs)
    {
        return py_iter_ == rhs.py_iter_;
    }

private:
    py::iterator py_iter_;
};


PYBIND11_MODULE(delaunay_class, m)
{
    py::class_<Point>(m, "Point")
            .def(py::init<int, int>(),  py::arg("x"), py::arg("y"))
            .def(py::init<double, double>(), py::arg("x"), py::arg("y"))
            .def_property_readonly("x", &Point::x)
            .def_property_readonly("y", &Point::y)
            .def("__repr__",
            [](const Point &p) {
                std::string r("Point(");
                r += boost::lexical_cast<std::string>(p.x());
                r += ", ";
                r += boost::lexical_cast<std::string>(p.y());
                r += ")";
                return r;
            })
            ;

    py::class_<Vertex_handle>(m, "VertexHandle")
        .def_property_readonly(
        "point",
        [](const Vertex_handle& vertex_handle)
        {
            return vertex_handle->point();
        })
        ;

        py::class_<DT>(m, "DelaunayTriangulation")

            .def(py::init())

            .def("insert", [](DT & dt, const std::vector<double> & p) {
                  std::vector< std::pair<Point,unsigned> > points;
                  int num_points = p.size()/2;
                  // start adding at the end of the current table
                  int start = dt.number_of_vertices();
                  for(std::size_t i = 0; i < num_points; ++i)
                  {
                    // add index information to form face table later
                     points.push_back( std::make_pair( Point(p[i*2+0],p[i*2+1]), start) );
                     start += 1;
                  }
                  return dt.insert(points.begin(),points.end());
                })

            .def("remove", [](DT & dt, const std::vector<unsigned int> & to_remove) {
                    int num_to_remove= to_remove.size();
                    std::vector<Vertex_handle> handles;
                    for (Vi vi = dt.finite_vertices_begin(); vi != dt.finite_vertices_end(); vi++){
                        handles.push_back(vi);
                    }
                    for(std::size_t i=0; i < num_to_remove; ++i){
                        dt.remove(handles[to_remove[i]]);
                    }
                    return dt;
                })

            .def("move", [](DT & dt, const std::vector<unsigned int> & to_move, const std::vector<double> & new_positions){
                    std::vector<Vertex_handle> handles;
                    std::vector<Point> new_pos;
                    int num_to_move= to_move.size();
                    // store all vertex handles
                    for (Vi vi = dt.finite_vertices_begin(); vi != dt.finite_vertices_end(); vi++){
                        handles.push_back(vi);
                    }
                    // store new positions as a vector of Point
                    for(std::size_t i = 0; i < num_to_move; ++i)
                    {
                       new_pos.push_back(Point(new_positions[2*i], new_positions[2*i+1]));
                    }
                    //
                    for(std::size_t i = 0; i < num_to_move; ++i)
                    {
                       dt.move( handles[to_move[i]], new_pos[i] );
                    }
                    return dt;
                    })

            .def("number_of_vertices", &DT::number_of_vertices)

            .def("number_of_faces", [](DT & dt){
                int count=0;
                for(DT::Finite_faces_iterator fit = dt.finite_faces_begin();
                fit != dt.finite_faces_end(); ++fit) {
                    count += 1;
                }
                    return count;
                })

            .def("finite_vertices", [](DT & dt) -> py::iterator
             {
                 return py::make_iterator(dt.finite_vertices_begin(), dt.finite_vertices_end());
             })

            .def("get_finite_cells", [](DT & dt)
            {
              // ouput the face table
              // YOU MUST CALL get_finite_vertices before if any incremental operations
              // were performed
              std::vector<int> faces;
              faces.resize(dt.number_of_faces()*3);

              int i=0;
              for(DT::Finite_faces_iterator fit = dt.finite_faces_begin();
                fit != dt.finite_faces_end(); ++fit) {

                DT::Face_handle face = fit;
                faces[i*3]=face->vertex(0)->info();
                faces[i*3+1]=face->vertex(1)->info();
                faces[i*3+2]=face->vertex(2)->info();
                i+=1;
              }
              ssize_t              soint      = sizeof(int);
              ssize_t              num_faces = faces.size()/3;
              ssize_t              ndim      = 2;
              std::vector<ssize_t> shape     = {num_faces, 3};
              std::vector<ssize_t> strides   = {soint*3, soint};

              // return 2-D NumPy array
              return py::array(py::buffer_info(
                faces.data(),                           /* data as contiguous array  */
                sizeof(int),                          /* size of one scalar        */
                py::format_descriptor<int>::format(), /* data type                 */
                2,                                    /* number of dimensions      */
                shape,                                   /* shape of the matrix       */
                strides                                  /* strides for each axis     */
              ));
            })

            .def("get_finite_vertices", [](DT & dt)
             {
               // ouput the vertices
               std::vector<double> vertices;
               vertices.resize(dt.number_of_vertices()*2);

               int i=0;
               for(DT::Finite_vertices_iterator fit = dt.finite_vertices_begin();
                 fit != dt.finite_vertices_end(); ++fit) {

                 Vertex_handle vertex = fit;
                 // critical! update the point index table so faces comes out correctly
                 vertex->info() = i;
                 vertices[i*2]=vertex->point().x();
                 vertices[i*2+1]=vertex->point().y();
                 i+=1;
               }
               ssize_t              sdble   = sizeof(double);
               ssize_t              num_vertices = vertices.size()/2;
               ssize_t              ndim      = 2;
               std::vector<ssize_t> shape     = {num_vertices, 2};
               std::vector<ssize_t> strides   = {sdble*2, sdble};

               // return 2-D NumPy array
               return py::array(py::buffer_info(
                 vertices.data(),                           /* data as contiguous array  */
                 sizeof(double),                          /* size of one scalar        */
                 py::format_descriptor<double>::format(), /* data type                 */
                 2,                                    /* number of dimensions      */
                 shape,                                   /* shape of the matrix       */
                 strides                                  /* strides for each axis     */
               ));
             })
            ;


    py::class_<DT::Finite_vertices_iterator::value_type>(m, "Vertex")
            .def_property_readonly(
                "point", [](DT::Finite_vertices_iterator::value_type& vertex)
                {
                    return vertex.point();
                }
            )
            ;

    py::class_<DT::Finite_faces_iterator::value_type>(m, "Face")
            .def("vertex_handle",
                [](DT::Finite_faces_iterator::value_type& face, int index)
                {
                    return face.vertex(index);
                },
                py::arg("index")
            )
            ;

}
