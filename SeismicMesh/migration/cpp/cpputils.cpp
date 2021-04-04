#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <CGAL/Circle_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Iso_cuboid_3.h>
#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/Sphere_3.h>
#include <CGAL/intersections.h>

#include <CGAL/Bbox_2.h>
#include <CGAL/Bbox_3.h>

#include <assert.h>
#include <vector>

#include <iterator>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
typedef K::Point_3 Point3;
typedef K::Circle_2 Circle;
typedef K::Sphere_3 Sphere;
typedef K::Iso_rectangle_2 Rectangle;
typedef K::Iso_cuboid_3 Cuboid;

namespace py = pybind11;

// build the vertex-to-element connectivity here for 2d
std::vector<std::vector<int>> build_vtoe2(std::vector<double> &p,
                                          std::vector<int> &t) {
  int sz = p.size() / 2;
  int sz_e = t.size() / 3;
  std::vector<std::vector<int>> vtoe(sz);
  // assume 3 for a reasonable amount of elemental neighbors
  for (auto &v : vtoe) {
    v.reserve(10);
  }
  for (std::size_t e = 0; e < sz_e; ++e) {
    for (std::size_t n = 0; n < 3; n++) {
      vtoe[t[3 * e + n]].push_back(e);
    }
  }
  vtoe.resize(vtoe.capacity());
  return vtoe;
}
// build the vertex-to-element connectivity here for 3d
std::vector<std::vector<int>> build_vtoe3(std::vector<double> &p,
                                          std::vector<int> &t) {
  int sz = p.size() / 3;
  int sz_e = t.size() / 4;
  std::vector<std::vector<int>> vtoe(sz);
  for (auto &v : vtoe) {
    v.reserve(10);
  }
  for (std::size_t e = 0; e < sz_e; ++e) {
    for (std::size_t n = 0; n < 4; n++) {
      vtoe[t[4 * e + n]].push_back(e);
    }
  }
  vtoe.resize(vtoe.capacity());
  return vtoe;
}

// fixed size calculation of 4x4 determinant
double calc_4x4dete(std::vector<double> &m) {
  return m[12] * m[9] * m[6] * m[3] - m[8] * m[13] * m[6] * m[3] -
         m[12] * m[5] * m[10] * m[3] + m[4] * m[13] * m[10] * m[3] +
         m[8] * m[5] * m[14] * m[3] - m[4] * m[9] * m[14] * m[3] -
         m[12] * m[9] * m[2] * m[7] + m[8] * m[13] * m[2] * m[7] +
         m[12] * m[1] * m[10] * m[7] - m[0] * m[13] * m[10] * m[7] -
         m[8] * m[1] * m[14] * m[7] + m[0] * m[9] * m[14] * m[7] +
         m[12] * m[5] * m[2] * m[11] - m[4] * m[13] * m[2] * m[11] -
         m[12] * m[1] * m[6] * m[11] + m[0] * m[13] * m[6] * m[11] +
         m[4] * m[1] * m[14] * m[11] - m[0] * m[5] * m[14] * m[11] -
         m[8] * m[5] * m[2] * m[15] + m[4] * m[9] * m[2] * m[15] +
         m[8] * m[1] * m[6] * m[15] - m[0] * m[9] * m[6] * m[15] -
         m[4] * m[1] * m[10] * m[15] + m[0] * m[5] * m[10] * m[15];
}

// determine which rank points need to be exported to in 2D
std::vector<double> c_where_to2(std::vector<double> &points,
                                std::vector<int> &faces,
                                std::vector<double> &llc,
                                std::vector<double> &urc, int rank) {
  int num_faces = faces.size() / 3;
  int num_points = points.size() / 2;

  // Determine which rank to send the vertex (exports)
  // exports[iv] is either 0 or 1 (0 for block owned by rank-1 and 1 for block
  // owned by rank+1)
  std::vector<int> exports;
  exports.resize(num_points, -1);

  auto vtoe_map = build_vtoe2(points, faces);

  double kount_below = 0.0;
  double kount_above = 0.0;

 Rectangle rect1 = Rectangle(Point(llc[0 * 2], llc[0 * 2 + 1]),Point(urc[0 * 2], urc[0 * 2 + 1]));
 Rectangle rect2 = Rectangle(Point(llc[1 * 2], llc[1 * 2 + 1]),Point(urc[1 * 2], urc[1 * 2 + 1]));

 // create bbox of rectangles to filter tests
 CGAL::Bbox_2 rect1_bbox = CGAL::Bbox_2(llc[0], llc[1], urc[0], urc[1]);
 CGAL::Bbox_2 rect2_bbox = CGAL::Bbox_2(llc[2], llc[3], urc[2], urc[3]);

  // For each point in points
  for (std::size_t iv = 0; iv < num_points; ++iv) {

    for (std::size_t ie = 0; ie < vtoe_map[iv].size(); ++ie) {

      int nei_ele = vtoe_map[iv][ie];

      // Indices of element into points
      int nm1 = faces[nei_ele * 3];
      int nm2 = faces[nei_ele * 3 + 1];
      int nm3 = faces[nei_ele * 3 + 2];

      // Test if this element has a bounding box that intersects with one of the rectangles
      // otherwise go no further
      std::array<double, 3> ele_x = {points[nm1*2], points[nm2*2], points[nm3*2]};
      std::array<double, 3> ele_y = {points[nm1*2+1], points[nm2*2+1], points[nm3*2+1]};

      auto min_x = std::min_element(std::begin(ele_x), std::end(ele_x));
      auto min_y = std::min_element(std::begin(ele_y), std::end(ele_y));

      auto max_x = std::max_element(std::begin(ele_x), std::end(ele_x));
      auto max_y = std::max_element(std::begin(ele_y), std::end(ele_y));

      CGAL::Bbox_2 ball_bbox = CGAL::Bbox_2(*min_x, *min_y, *max_x, *max_y);

      if(CGAL::do_overlap(ball_bbox, rect1_bbox) or CGAL::do_overlap(ball_bbox, rect2_bbox)){

        // Coordinates of each vertex of element
        Point pnm1 = Point(points[nm1 * 2], points[nm1 * 2 + 1]);
        Point pnm2 = Point(points[nm2 * 2], points[nm2 * 2 + 1]);
        Point pnm3 = Point(points[nm3 * 2], points[nm3 * 2 + 1]);
        bool isCollinear = CGAL::collinear(pnm1, pnm2, pnm3);
        // if this happens we cannot test accurately
        if (isCollinear) {
          // std::cout<<"point is collinear" << std::endl;
          continue;
        }
        // Calculate circumball of element
        Point cc = CGAL::circumcenter(pnm1, pnm2, pnm3);
        double sqr_radius = CGAL::squared_radius(pnm1, pnm2, pnm3);
        Circle circ = Circle(cc, sqr_radius, CGAL::CLOCKWISE);
        // Does this circumball intersect with box above or box below?
        bool intersects;
        intersects = CGAL::do_intersect(circ, rect1);
        if (intersects){
            exports[iv] = 0;
            kount_below += 1.0;
            break;
        }
        intersects = CGAL::do_intersect(circ, rect2);
        if (intersects){
            exports[iv] = 1;
            kount_above += 1.0;
            break;
        }
      }
    }
  }

  std::vector<double> pointsToMigrate;
  pointsToMigrate.resize(1 + num_points * 3, -1);

  double below_count = 0.0;
  double above_count = 0.0;

  for (std::size_t iv = 0; iv < num_points; ++iv) {
    if (exports[iv] == -1){
        continue;
    }
    if (exports[iv] == 0) {
      pointsToMigrate[below_count * 3 + 0 + 3] = points[iv * 2];
      pointsToMigrate[below_count * 3 + 1 + 3] = points[iv * 2 + 1];
      pointsToMigrate[below_count * 3 + 2 + 3] = (double)iv;
      below_count += 1.0;
    }
    else if (exports[iv] ==1){
      pointsToMigrate[kount_below * 3 + above_count * 3 + 0 + 3] =
          points[iv * 2];
      pointsToMigrate[kount_below * 3 + above_count * 3 + 1 + 3] =
          points[iv * 2 + 1];
      pointsToMigrate[kount_below * 3 + above_count * 3 + 2 + 3] = (double)iv;
      above_count += 1.0;
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
py::array
where_to2(py::array_t<double, py::array::c_style | py::array::forcecast> points,
          py::array_t<int, py::array::c_style | py::array::forcecast> faces,
          py::array_t<double, py::array::c_style | py::array::forcecast> llc,
          py::array_t<double, py::array::c_style | py::array::forcecast> urc,
          int rank) {
  int num_faces = faces.size() / 3;
  int num_points = points.size() / 2;

  // allocate std::vector (to pass to the C++ function)
  std::vector<double> cpppoints(num_points * 2);
  std::vector<int> cppfaces(num_faces * 3);
  std::vector<double> cppllc(4);
  std::vector<double> cppurc(4);

  // copy py::array -> std::vector
  std::memcpy(cpppoints.data(), points.data(), num_points * 2 * sizeof(double));
  std::memcpy(cppfaces.data(), faces.data(), num_faces * 3 * sizeof(int));
  std::memcpy(cppllc.data(), llc.data(), 4 * sizeof(double));
  std::memcpy(cppurc.data(), urc.data(), 4 * sizeof(double));

  // call cpp code
  std::vector<double> pointsToMigrate =
      c_where_to2(cpppoints, cppfaces, cppllc, cppurc, rank);

  ssize_t sodble = sizeof(double);
  ssize_t ndim = 2;
  std::vector<ssize_t> shape = {1 + num_points, 3};
  std::vector<ssize_t> strides = {sodble * 3, sodble};

  // return 2-D NumPy array
  return py::array(
      py::buffer_info(pointsToMigrate.data(), /* data as contiguous array  */
                      sizeof(double),         /* size of one scalar        */
                      py::format_descriptor<double>::format(), /* data type */
                      2,      /* number of dimensions      */
                      shape,  /* shape of the matrix       */
                      strides /* strides for each axis     */
                      ));
}

// determine which rank points need to be exported to in 3D
std::vector<double> c_where_to3(std::vector<double> &points,
                                std::vector<int> &faces,
                                std::vector<double> &llc,
                                std::vector<double> &urc, int rank) {

  int num_faces = faces.size() / 4;
  int num_points = points.size() / 3;

  // Determine which rank to send the vertex (exports)
  // exports[iv] is either 0 or 1 (0 for block owned by rank-1 and 1 for block
  // owned by rank+1)
  std::vector<int> exports;
  exports.resize(num_points, -1);

  auto vtoe_map = build_vtoe3(points, faces);

  double kount_below = 0.0;
  double kount_above = 0.0;

  Cuboid cube1 = Cuboid(Point3(llc[0 * 3], llc[0 * 3 + 1], llc[0 * 3 + 2]),
                Point3(urc[0 * 3], urc[0 * 3 + 1], urc[0 * 3 + 2]));
  Cuboid cube2 = Cuboid(Point3(llc[1 * 3], llc[1 * 3 + 1], llc[1 * 3 + 2]),
                Point3(urc[1 * 3], urc[1 * 3 + 1], urc[1 * 3 + 2]));

  // create bbox of cuboids to filter tests
  CGAL::Bbox_3 cuboid1_bbox = CGAL::Bbox_3(llc[0], llc[1], llc[2] ,urc[0], urc[1], urc[2]);
  CGAL::Bbox_3 cuboid2_bbox = CGAL::Bbox_3(llc[3], llc[4], llc[5], urc[3], urc[4], urc[5]);

  // For each point in points
  for (std::size_t iv = 0; iv < num_points; ++iv) {

    for (std::size_t ie = 0; ie < vtoe_map[iv].size(); ++ie) {

      int nei_ele = vtoe_map[iv][ie];

      // Indices of element into points
      int nm1 = faces[nei_ele * 4];
      int nm2 = faces[nei_ele * 4 + 1];
      int nm3 = faces[nei_ele * 4 + 2];
      int nm4 = faces[nei_ele * 4 + 3];

      // Test if this element has a bounding box that intersects with one of the rectangles
      // otherwise go no further
      std::array<double, 4> ele_x = {points[nm1*3], points[nm2*3], points[nm3*3], points[nm4*3]};
      std::array<double, 4> ele_y = {points[nm1*3+1], points[nm2*3+1], points[nm3*3+1], points[nm4*3+1]};
      std::array<double, 4> ele_z = {points[nm1*3+2], points[nm2*3+2], points[nm3*3+2], points[nm4*3+2]};

      auto min_x = std::min_element(std::begin(ele_x), std::end(ele_x));
      auto min_y = std::min_element(std::begin(ele_y), std::end(ele_y));
      auto min_z = std::min_element(std::begin(ele_z), std::end(ele_z));

      auto max_x = std::max_element(std::begin(ele_x), std::end(ele_x));
      auto max_y = std::max_element(std::begin(ele_y), std::end(ele_y));
      auto max_z = std::max_element(std::begin(ele_z), std::end(ele_z));

      CGAL::Bbox_3 sphere_bbox = CGAL::Bbox_3(*min_x, *min_y, *min_z, *max_x, *max_y, *max_z);

      if(CGAL::do_overlap(sphere_bbox, cuboid1_bbox) or CGAL::do_overlap(sphere_bbox, cuboid2_bbox)){

         // Coordinates of each vertex of element
         Point3 pnm1 =
             Point3(points[nm1 * 3], points[nm1 * 3 + 1], points[nm1 * 3 + 2]);
         Point3 pnm2 =
             Point3(points[nm2 * 3], points[nm2 * 3 + 1], points[nm2 * 3 + 2]);
         Point3 pnm3 =
             Point3(points[nm3 * 3], points[nm3 * 3 + 1], points[nm3 * 3 + 2]);
         Point3 pnm4 =
             Point3(points[nm4 * 3], points[nm4 * 3 + 1], points[nm4 * 3 + 2]);

         if (volume(pnm1, pnm2, pnm3, pnm4) == 0) {
           continue;
         }

         bool isCoplanar = CGAL::coplanar(pnm1, pnm2, pnm3, pnm4);
         if (isCoplanar) {
           // std::cout<<"alert"<<std::endl;
           continue;
         }

         // Calculate circumball of element
         Point3 cc = CGAL::circumcenter(pnm1, pnm2, pnm3, pnm4);
         double sqr_radius = CGAL::squared_radius(pnm1, pnm2, pnm3, pnm4);
         Sphere sphere = Sphere(cc, sqr_radius, CGAL::CLOCKWISE);
         // Does this circumball intersect with box above or box below?
         bool intersects;
         intersects = CGAL::do_intersect(sphere, cube1);
         if (intersects){
             exports[iv] = 0;
             kount_below += 1.0;
             break;
         }
         intersects = CGAL::do_intersect(sphere, cube2);
         if (intersects){
              exports[iv] = 1;
              kount_above += 1.0;
              break;
          }
        }
     }
  }

  std::vector<double> pointsToMigrate;
  pointsToMigrate.resize(1 + num_points * 4, -1);

  double below_count = 0.0;
  double above_count = 0.0;

  for (std::size_t iv = 0; iv < num_points; ++iv) {
      if (exports[iv] == -1){
          continue;
      }
    if (exports[iv] == 0) {
      pointsToMigrate[below_count * 4 + 0 + 4] = points[iv * 3];
      pointsToMigrate[below_count * 4 + 1 + 4] = points[iv * 3 + 1];
      pointsToMigrate[below_count * 4 + 2 + 4] = points[iv * 3 + 2];
      pointsToMigrate[below_count * 4 + 3 + 4] = (double)iv;
      below_count += 1.0;
    }
    else if (exports[iv] == 1) {
      pointsToMigrate[kount_below * 4 + above_count * 4 + 0 + 4] =
          points[iv * 3];
      pointsToMigrate[kount_below * 4 + above_count * 4 + 1 + 4] =
          points[iv * 3 + 1];
      pointsToMigrate[kount_below * 4 + above_count * 4 + 2 + 4] =
          points[iv * 3 + 2];
      pointsToMigrate[kount_below * 4 + above_count * 4 + 3 + 4] = (double)iv;
      above_count += 1.0;
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
py::array
where_to3(py::array_t<double, py::array::c_style | py::array::forcecast> points,
          py::array_t<int, py::array::c_style | py::array::forcecast> faces,
          py::array_t<double, py::array::c_style | py::array::forcecast> llc,
          py::array_t<double, py::array::c_style | py::array::forcecast> urc,
          int rank) {

  int num_faces = faces.size() / 4;
  int num_points = points.size() / 3;

  // allocate std::vector (to pass to the C++ function)
  std::vector<double> cpppoints(num_points * 3);
  std::vector<int> cppfaces(num_faces * 4);
  std::vector<double> cppllc(6);
  std::vector<double> cppurc(6);

  // copy py::array -> std::vector
  std::memcpy(cpppoints.data(), points.data(), num_points * 3 * sizeof(double));
  std::memcpy(cppfaces.data(), faces.data(), num_faces * 4 * sizeof(int));
  std::memcpy(cppllc.data(), llc.data(), 6 * sizeof(double));
  std::memcpy(cppurc.data(), urc.data(), 6 * sizeof(double));

  // call cpp code
  std::vector<double> pointsToMigrate =
      c_where_to3(cpppoints, cppfaces, cppllc, cppurc, rank);

  ssize_t sodble = sizeof(double);
  ssize_t ndim = 2;
  std::vector<ssize_t> shape = {1 + num_points, 4};
  std::vector<ssize_t> strides = {sodble * 4, sodble};

  // return 2-D NumPy array
  return py::array(
      py::buffer_info(pointsToMigrate.data(), /* data as contiguous array  */
                      sizeof(double),         /* size of one scalar        */
                      py::format_descriptor<double>::format(), /* data type */
                      2,      /* number of dimensions      */
                      shape,  /* shape of the matrix       */
                      strides /* strides for each axis     */
                      ));
}

PYBIND11_MODULE(_cpputils, m) {
  m.def("where_to2", &where_to2);
  m.def("where_to3", &where_to3);
}
