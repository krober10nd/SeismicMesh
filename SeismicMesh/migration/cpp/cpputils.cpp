#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/numpy.h>

#include <CGAL/Kernel/global_functions.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/Circle_2.h>
#include <CGAL/intersections.h>

#include <assert.h>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
typedef K::Circle_2  Circle;
typedef K::Iso_rectangle_2 Rectangle;

namespace py = pybind11;

std::vector<double> c_circumballs2(std::vector<double> &vertices)
{
    int num_faces = vertices.size()/6;
    std::vector<double> circumcenters;
    for(std::size_t i=0; i < num_faces; ++i)
    {
        Point tmp_cc =
                CGAL::circumcenter(
                Point(vertices[i*6],vertices[i*6+1]),
                Point(vertices[i*6+2],vertices[i*6+3]),
                Point(vertices[i*6+4],vertices[i*6+5])
                );
        circumcenters.push_back(tmp_cc.x());
        circumcenters.push_back(tmp_cc.y());
        circumcenters.push_back(
                CGAL::squared_radius(
                    Point(vertices[i*6],vertices[i*6+1]),
                    Point(vertices[i*6+2],vertices[i*6+3]),
                    Point(vertices[i*6+4],vertices[i*6+5])
                    )
                );
    }
    return circumcenters;
}

py::array circumballs2(py::array_t<double, py::array::c_style | py::array::forcecast> vertices)
{
    // each triangle has 3 vertices with 2 coordinates each
    int sz = vertices.shape()[0];
    std::vector<double> cppvertices(sz);
    std::memcpy(cppvertices.data(),vertices.data(),sz*sizeof(double));
    std::vector<double> circumcenters = c_circumballs2(cppvertices);
    ssize_t              soreal      = sizeof(double);
    ssize_t              num_points = circumcenters.size()/3;
    ssize_t              ndim      = 2;
    std::vector<ssize_t> shape     = {num_points, 3};
    std::vector<ssize_t> strides   = {soreal*3, soreal};
    // return 2-D NumPy array
    return py::array(py::buffer_info(
        circumcenters.data(),                    /* data as contiguous array  */
        sizeof(double),                          /* size of one scalar        */
        py::format_descriptor<double>::format(), /* data type                 */
        2,                                       /* number of dimensions      */
        shape,                                   /* shape of the matrix       */
        strides                                  /* strides for each axis     */
  ));
}

// determine which rank points need to be exported to
std::vector<double> c_where_to2(std::vector<double> &points, std::vector<int> &faces,
                             std::vector<double> &llc, std::vector<double> &urc, int rank)
{
    int num_faces = faces.size()/3;
    int num_points = points.size()/2;

    // Step 1. Determine all elements connected to each point in points.
    std::vector<int> vtoe;
    std::vector<int> nne;
    // assume each vertex has a max. of 20 elements neigh.
    vtoe.resize(num_points*20, -1);
    nne.resize(num_points, 0);
    for(std::size_t ie = 0; ie < num_faces; ++ie ) {
        for(std::size_t iv =0; iv < 3; ++iv ) {
            int nm1 = faces[ie*3+iv];
            vtoe[nm1*20 + nne[nm1]] = ie;
            nne[nm1] += 1;
            }
        }
    // Step 2. Determine which rank to send the vertex (exports)
    // exports[iv] is either 0 or 1 (0 for block owned by rank-1 and 1 for block owned by rank+1)
    std::vector<int> exports;
    exports.resize(num_points,-1);
    // For each point in points
    for(std::size_t iv=0; iv < num_points; ++iv)
    {
        // For all connected elements to point iv
        for(std::size_t ic=0; ic < nne[iv]; ++ic)
        {
            int nei_ele = vtoe[iv*20 + ic];
            // Indices of element into points
            int nm1 = faces[nei_ele*3];
            int nm2 = faces[nei_ele*3+1];
            int nm3 = faces[nei_ele*3+2];
            // Coordinates of each vertex of element
            Point pnm1 = Point(points[nm1*2], points[nm1*2+1]);
            Point pnm2 = Point(points[nm2*2], points[nm2*2+1]);
            Point pnm3 = Point(points[nm3*2], points[nm3*2+1]);
            // Calculate circumball of element
            Point cc = CGAL::circumcenter(pnm1, pnm2, pnm3);
            double sqr_radius = CGAL::squared_radius(pnm1, pnm2, pnm3);
            Circle circ = Circle(cc, sqr_radius, CGAL::CLOCKWISE);
            // Does this circumball intersect with box above or box below?
            for(std::size_t bx=0; bx< 2; ++bx ){
                Rectangle rect = Rectangle(Point(llc[bx*2], llc[bx*2 +1]),
                        Point(urc[bx*2], urc[bx*2+1]));
                bool intersects = CGAL::do_intersect(circ, rect);
                if(intersects){
                    exports[iv] = bx;
                }
            }
        }
    }

    std::vector<double> pointsToMigrate;
    pointsToMigrate.resize(num_points*2,-1);

    double kount_below = 0.0;
    for(std::size_t iv=0; iv < num_points; ++iv)
    {
        if(exports[iv]==0)
        {
            pointsToMigrate[kount_below*2+2]=points[iv*2];
            pointsToMigrate[kount_below*2+1+2]=points[iv*2+1];
            kount_below += 1.0;
        }
    }

    double kount_above = 0.0;
    for(std::size_t iv=0; iv < num_points; ++iv)
    {
        if(exports[iv]==1)
        {
            pointsToMigrate[kount_below*2 + kount_above*2+2]=points[iv*2];
            pointsToMigrate[kount_below*2 + kount_above*2+1+2]=points[iv*2+1];
            kount_above += 1.0;
        }
    }
    pointsToMigrate[0] = kount_below;
    pointsToMigrate[1] = kount_above;
    return pointsToMigrate;
}

// ----------------
// Python interface for c_where_to2
// ----------------
py::array where_to2(py::array_t<double, py::array::c_style | py::array::forcecast> points,
                    py::array_t<int, py::array::c_style | py::array::forcecast> faces,
                    py::array_t<double, py::array::c_style | py::array::forcecast> llc,
                    py::array_t<double, py::array::c_style | py::array::forcecast> urc,
                    int rank
                    )
{
  int num_faces = faces.size()/3;
  int num_points = points.size()/2;

  // allocate std::vector (to pass to the C++ function)
  std::vector<double> cpppoints(num_points*2);
  std::vector<int> cppfaces(num_faces*3);
  std::vector<double> cppllc(4);
  std::vector<double> cppurc(4);

  // copy py::array -> std::vector
  std::memcpy(cpppoints.data(),points.data(),num_points*2*sizeof(double));
  std::memcpy(cppfaces.data(),faces.data(),num_faces*3*sizeof(int));
  std::memcpy(cppllc.data(), llc.data(),4*sizeof(double));
  std::memcpy(cppurc.data(), urc.data(),4*sizeof(double));

  // call cpp code
  std::vector<double> pointsToMigrate = c_where_to2(cpppoints, cppfaces, cppllc, cppurc, rank);

  ssize_t              sodble    = sizeof(double);
  ssize_t              ndim      = 2;
  std::vector<ssize_t> shape     = {num_points, 2};
  std::vector<ssize_t> strides   = {sodble*2, sodble};

  // return 2-D NumPy array
  return py::array(py::buffer_info(
    pointsToMigrate.data(),                /* data as contiguous array  */
    sizeof(double),                       /* size of one scalar        */
    py::format_descriptor<double>::format(), /* data type                 */
    2,                                    /* number of dimensions      */
    shape,                                   /* shape of the matrix       */
    strides                                  /* strides for each axis     */
  ));
}


// Does a circle intersect with a rectangle?
bool c_do_intersect2(std::vector<double> &center, double radius, std::vector<double> &llc, std::vector<double> &urc )
{
  double sqr_r1 = radius*radius;
  Circle circ1 = Circle(Point(center[0],center[1]), sqr_r1, CGAL::CLOCKWISE);
  Rectangle rect1 = Rectangle(Point(llc[0], llc[1]), Point(urc[0], urc[1]));
  return CGAL::do_intersect(circ1, rect1);
}

// ----------------
// Python interface for do_intersect2
// ----------------
bool do_intersect2(
        py::array_t<double, py::array::c_style | py::array::forcecast> center,
        double radius,
        py::array_t<double, py::array::c_style | py::array::forcecast> llc,
        py::array_t<double, py::array::c_style | py::array::forcecast> urc
        )
{
    std::vector<double> cppCenter(2);
    std::vector<double> cppLLC(2);
    std::vector<double> cppURC(2);
    std::memcpy(cppCenter.data(),center.data(),2*sizeof(double));
    std::memcpy(cppLLC.data(),llc.data(),2*sizeof(double));
    std::memcpy(cppURC.data(),urc.data(),2*sizeof(double));
    bool result = c_do_intersect2(cppCenter,radius,cppLLC,cppURC);
    return result;
}


// A table with the elements connected to a vertex
std::vector<int> c_vertex_to_elements(std::vector<int> &faces, int &num_points, int &num_faces)
{
    std::vector<int> vtoe;
    std::vector<int> nne;
    // assume each vertex has a max. of 20 elements neigh.
    vtoe.resize(num_points*20, -1);
    nne.resize(num_points, 0);
    for(std::size_t ie = 0; ie < num_faces; ++ie ) {
        for(std::size_t iv =0; iv < 3; ++iv ) {
            int nm1 = faces[ie*3+iv];
            vtoe[nm1*20 + nne[nm1]] = ie;
            nne[nm1] += 1;
            }
        }
    return vtoe;
}


// ----------------
// Python interface for vertex_to_elements
// ----------------
py::array vertex_to_elements(py::array_t<int, py::array::c_style | py::array::forcecast> faces,
        int num_points, int num_faces)
{

  // check input dimensions
  if ( faces.ndim() != 2 )
    throw std::runtime_error("Input should be a 2D NumPy array");

  if ( num_points < 3 )
      throw std::runtime_error("Too few points!");

  if (num_faces != faces.size()/3)
      throw std::runtime_error("Number of faces doesn't match!");

  // allocate std::vector (to pass to the C++ function)
  std::vector<int> cppfaces(num_faces*3);

  // copy py::array -> std::vector
  std::memcpy(cppfaces.data(),faces.data(),num_faces*3*sizeof(int));
  std::vector<int> vtoe = c_vertex_to_elements(cppfaces, num_points, num_faces);

  ssize_t              soint      = sizeof(int);
  ssize_t              ndim      = 2;
  std::vector<ssize_t> shape     = {num_points, 20};
  std::vector<ssize_t> strides   = {soint*20, soint};

  // return 2-D NumPy array
  return py::array(py::buffer_info(
    vtoe.data(),                           /* data as contiguous array  */
    sizeof(int),                          /* size of one scalar        */
    py::format_descriptor<int>::format(), /* data type                 */
    2,                                    /* number of dimensions      */
    shape,                                   /* shape of the matrix       */
    strides                                  /* strides for each axis     */
  ));
}



PYBIND11_MODULE(cpputils, m) {
    m.def("do_intersect2", &do_intersect2);
    m.def("where_to2", &where_to2);
    m.def("vertex_to_elements", &vertex_to_elements);
}
