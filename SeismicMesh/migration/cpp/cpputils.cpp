#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/numpy.h>

#include <CGAL/Kernel/global_functions.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/Iso_cuboid_3.h>
#include <CGAL/Circle_2.h>
#include <CGAL/Sphere_3.h>
#include <CGAL/intersections.h>

#include <assert.h>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
typedef K::Point_3 Point3;
typedef K::Circle_2 Circle;
typedef K::Sphere_3 Sphere;
typedef K::Iso_rectangle_2 Rectangle;
typedef K::Iso_cuboid_3 Cuboid;

namespace py = pybind11;


// fixed size calculation of 4x4 determinant
double calc_4x4dete(std::vector<double> &m) {
    return
         m[12] * m[9]  * m[6]  * m[3]   -  m[8] * m[13] * m[6]  * m[3]   -
         m[12] * m[5]  * m[10] * m[3]   +  m[4] * m[13] * m[10] * m[3]   +
         m[8]  * m[5]  * m[14] * m[3]   -  m[4] * m[9]  * m[14] * m[3]   -
         m[12] * m[9]  * m[2]  * m[7]   +  m[8] * m[13] * m[2]  * m[7]   +
         m[12] * m[1]  * m[10] * m[7]   -  m[0] * m[13] * m[10] * m[7]   -
         m[8]  * m[1]  * m[14] * m[7]   +  m[0] * m[9]  * m[14] * m[7]   +
         m[12] * m[5]  * m[2]  * m[11]  -  m[4] * m[13] * m[2]  * m[11]  -
         m[12] * m[1]  * m[6]  * m[11]  +  m[0] * m[13] * m[6]  * m[11]  +
         m[4]  * m[1]  * m[14] * m[11]  -  m[0] * m[5]  * m[14] * m[11]  -
         m[8]  * m[5]  * m[2]  * m[15]  +  m[4] * m[9]  * m[2]  * m[15]  +
         m[8]  * m[1]  * m[6]  * m[15]  -  m[0] * m[9]  * m[6]  * m[15]  -
         m[4]  * m[1]  * m[10] * m[15]  +  m[0] * m[5]  * m[10] * m[15];
}

// determine which rank points need to be exported to in 2D
std::vector<double> c_where_to2(std::vector<double> &points, std::vector<int> &faces,
                             std::vector<int> &vtoe, std::vector<int> &ptr,
                             std::vector<double> &llc, std::vector<double> &urc, int rank)
{
    int num_faces = faces.size()/3;
    int num_points = points.size()/2;

    // Determine which rank to send the vertex (exports)
    // exports[iv] is either 0 or 1 (0 for block owned by rank-1 and 1 for block owned by rank+1)
    std::vector<int> exports;
    exports.resize(num_points,-1);
    // For each point in points
    for(std::size_t iv=0; iv < num_points; ++iv)
    {
        int nneis = ptr[iv+1]-ptr[iv];
        // For all connected elements to point iv
        for(std::size_t ic=0; ic < nneis; ++ic)
        {
            int nei_ele = vtoe[ptr[iv]+ic];
            // Indices of element into points
            int nm1 = faces[nei_ele*3];
            int nm2 = faces[nei_ele*3+1];
            int nm3 = faces[nei_ele*3+2];
            // Coordinates of each vertex of element
            Point pnm1 = Point(points[nm1*2], points[nm1*2+1]);
            Point pnm2 = Point(points[nm2*2], points[nm2*2+1]);
            Point pnm3 = Point(points[nm3*2], points[nm3*2+1]);
            bool isCollinear = CGAL::collinear (pnm1, pnm2, pnm3);
            // if this happens we cannot test accurately
            if(isCollinear){
                //std::cout<<"point is collinear" << std::endl;
                continue;
            }
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
    pointsToMigrate.resize(1+num_points*3,-1);

    double kount_below = 0.0;
    for(std::size_t iv=0; iv < num_points; ++iv)
    {
        if(exports[iv]==0)
        {
            pointsToMigrate[kount_below*3+0+3]=points[iv*2];
            pointsToMigrate[kount_below*3+1+3]=points[iv*2+1];
            pointsToMigrate[kount_below*3+2+3]=(double)iv;
            kount_below += 1.0;
        }
    }

    double kount_above = 0.0;
    for(std::size_t iv=0; iv < num_points; ++iv)
    {
        if(exports[iv]==1)
        {
            pointsToMigrate[kount_below*3 + kount_above*3+0+3]=points[iv*2];
            pointsToMigrate[kount_below*3 + kount_above*3+1+3]=points[iv*2+1];
            pointsToMigrate[kount_below*3 + kount_above*3+2+3]=(double)iv;
            kount_above += 1.0;
        }
    }
    pointsToMigrate[0] = kount_below;
    pointsToMigrate[1] = kount_above;
    pointsToMigrate[2] = 0.0;

    return pointsToMigrate;
}

// ----------------
// Python interface for c_where_to2
// ----------------
py::array where_to2(py::array_t<double, py::array::c_style | py::array::forcecast> points,
                    py::array_t<int, py::array::c_style | py::array::forcecast> faces,
                    py::array_t<int, py::array::c_style | py::array::forcecast> vtoe,
                    py::array_t<int, py::array::c_style | py::array::forcecast> ptr,
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
  std::vector<int> cppvtoe(num_faces*3);
  std::vector<int> cppptr(num_points+1);
  std::vector<double> cppllc(4);
  std::vector<double> cppurc(4);

  // copy py::array -> std::vector
  std::memcpy(cpppoints.data(),points.data(),num_points*2*sizeof(double));
  std::memcpy(cppfaces.data(),faces.data(),num_faces*3*sizeof(int));
  std::memcpy(cppvtoe.data(),vtoe.data(),num_faces*3*sizeof(int));
  std::memcpy(cppptr.data(),ptr.data(),(num_points+1)*sizeof(int));
  std::memcpy(cppllc.data(), llc.data(),4*sizeof(double));
  std::memcpy(cppurc.data(), urc.data(),4*sizeof(double));

  // call cpp code
  std::vector<double> pointsToMigrate = c_where_to2(cpppoints, cppfaces, cppvtoe, cppptr, cppllc, cppurc, rank);

  ssize_t              sodble    = sizeof(double);
  ssize_t              ndim      = 2;
  std::vector<ssize_t> shape     = {1 + num_points, 3};
  std::vector<ssize_t> strides   = {sodble*3, sodble};

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

// determine which rank points need to be exported to in 3D
std::vector<double> c_where_to3(std::vector<double> &points, std::vector<int> &faces,
                             std::vector<int> &vtoe, std::vector<int> &ptr,
                             std::vector<double> &llc, std::vector<double> &urc, int rank)
{
    int num_faces = faces.size()/4;
    int num_points = points.size()/3;

    // Determine which rank to send the vertex (exports)
    // exports[iv] is either 0 or 1 (0 for block owned by rank-1 and 1 for block owned by rank+1)
    std::vector<int> exports;
    exports.resize(num_points,-1);
    // For each point in points
    for(std::size_t iv=0; iv < num_points; ++iv)
    {
        int nneis = ptr[iv+1]-ptr[iv];
        // For all connected elements to point iv
        for(std::size_t ic=0; ic < nneis; ++ic)
        {
            int nei_ele = vtoe[ptr[iv]+ic];
            // Indices of element into points
            int nm1 = faces[nei_ele*4];
            int nm2 = faces[nei_ele*4+1];
            int nm3 = faces[nei_ele*4+2];
            int nm4 = faces[nei_ele*4+3];

            // Coordinates of each vertex of element
            Point3 pnm1 = Point3(points[nm1*3], points[nm1*3+1], points[nm1*3+2]);
            Point3 pnm2 = Point3(points[nm2*3], points[nm2*3+1], points[nm2*3+2]);
            Point3 pnm3 = Point3(points[nm3*3], points[nm3*3+1], points[nm3*3+2]);
            Point3 pnm4 = Point3(points[nm4*3], points[nm4*3+1], points[nm4*3+2]);

            if(volume(pnm1, pnm2, pnm3, pnm4) == 0){
                continue;
            }

            bool isCoplanar = CGAL::coplanar(pnm1, pnm2, pnm3, pnm4);

            if(isCoplanar){
                //std::cout<<"alert"<<std::endl;
                continue;
            }
            //// Calculate circumball of element
            Point3 cc = CGAL::circumcenter(pnm1, pnm2, pnm3, pnm4);
            double sqr_radius = CGAL::squared_radius(pnm1, pnm2, pnm3, pnm4);
            Sphere sphere = Sphere(cc, sqr_radius, CGAL::CLOCKWISE);
            // Does this circumball intersect with box above or box below?
            for(std::size_t bx=0; bx< 2; ++bx ){
                Cuboid cube = Cuboid(Point3(llc[bx*3], llc[bx*3 +1], llc[bx*3 +2]),
                        Point3(urc[bx*3], urc[bx*3+1], urc[bx*3+2]));
                bool intersects = CGAL::do_intersect(sphere, cube);
                if(intersects){
                    exports[iv] = bx;
                }
            }
        }
    }

    std::vector<double> pointsToMigrate;
    pointsToMigrate.resize(1+num_points*4,-1);

    double kount_below = 0.0;
    for(std::size_t iv=0; iv < num_points; ++iv)
    {
        if(exports[iv]==0)
        {
            pointsToMigrate[kount_below*4+0+4]=points[iv*3];
            pointsToMigrate[kount_below*4+1+4]=points[iv*3+1];
            pointsToMigrate[kount_below*4+2+4]=points[iv*3+2];
            pointsToMigrate[kount_below*4+3+4]=(double)iv;
            kount_below += 1.0;
        }
    }

    double kount_above = 0.0;
    for(std::size_t iv=0; iv < num_points; ++iv)
    {
        if(exports[iv]==1)
        {
            pointsToMigrate[kount_below*4 + kount_above*4+0+4]=points[iv*3];
            pointsToMigrate[kount_below*4 + kount_above*4+1+4]=points[iv*3+1];
            pointsToMigrate[kount_below*4 + kount_above*4+2+4]=points[iv*3+2];
            pointsToMigrate[kount_below*4 + kount_above*4+3+4]=(double)iv;
            kount_above += 1.0;
        }
    }
    pointsToMigrate[0] = kount_below;
    pointsToMigrate[1] = kount_above;
    pointsToMigrate[2] = 0.0;
    pointsToMigrate[3] = 0.0;

    return pointsToMigrate;
}

// ----------------
// Python interface for c_where_to3
// ----------------
py::array where_to3(py::array_t<double, py::array::c_style | py::array::forcecast> points,
                    py::array_t<int, py::array::c_style | py::array::forcecast> faces,
                    py::array_t<int, py::array::c_style | py::array::forcecast> vtoe,
                    py::array_t<int, py::array::c_style | py::array::forcecast> ptr,
                    py::array_t<double, py::array::c_style | py::array::forcecast> llc,
                    py::array_t<double, py::array::c_style | py::array::forcecast> urc,
                    int rank
                    )
{
  int num_faces = faces.size()/4;
  int num_points = points.size()/3;

  // allocate std::vector (to pass to the C++ function)
  std::vector<double> cpppoints(num_points*3);
  std::vector<int> cppfaces(num_faces*4);
  std::vector<int> cppvtoe(num_faces*4);
  std::vector<int> cppptr(num_points+1);
  std::vector<double> cppllc(6);
  std::vector<double> cppurc(6);

  // copy py::array -> std::vector
  std::memcpy(cpppoints.data(),points.data(),num_points*3*sizeof(double));
  std::memcpy(cppfaces.data(),faces.data(),num_faces*4*sizeof(int));
  std::memcpy(cppvtoe.data(),vtoe.data(),num_faces*4*sizeof(int));
  std::memcpy(cppptr.data(),ptr.data(),(num_points+1)*sizeof(int));
  std::memcpy(cppllc.data(), llc.data(),6*sizeof(double));
  std::memcpy(cppurc.data(), urc.data(),6*sizeof(double));

  // call cpp code
  std::vector<double> pointsToMigrate = c_where_to3(cpppoints, cppfaces, cppvtoe, cppptr, cppllc, cppurc, rank);

  ssize_t              sodble    = sizeof(double);
  ssize_t              ndim      = 2;
  std::vector<ssize_t> shape     = {1 + num_points, 4};
  std::vector<ssize_t> strides   = {sodble*4, sodble};

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

PYBIND11_MODULE(cpputils, m) {
    m.def("where_to2", &where_to2);
    m.def("where_to3", &where_to3);
}
