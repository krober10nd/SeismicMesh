#include <vector>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

double l2_norm(std::vector<double> const& u) {
    double accum = 0.;
    for (int i = 0; i < 3; ++i) {
        accum += u[i] * u[i];
    }
    return sqrt(accum);
}

std::vector<double> cross_product(std::vector<double> &vect_A, std::vector<double> &vect_B)
{
    std::vector<double> cross_P;
    cross_P.resize(3);
    cross_P[0] = vect_A[1] * vect_B[2] - vect_A[2] * vect_B[1];
    cross_P[1] = vect_A[2] * vect_B[0] - vect_A[0] * vect_B[2];
    cross_P[2] = vect_A[0] * vect_B[1] - vect_A[1] * vect_B[0];
    return cross_P;
}

double dot_product(std::vector<double> &vect_A, std::vector<double> &vect_B)
{
    double product = 0.0;
    for (unsigned int i = 0; i < 3; i++)
        product = product + vect_A[i] * vect_B[i];
    return product;
}

std::vector<double> c_calc_dihedral_angles(std::vector<double> &points, std::vector<int> &cells)
{
    // compute the 6 dihedral angles of all tetrahedrons in the mesh
  static const std::size_t edges[6][2] = {{2, 3}, {1, 3}, {1, 2}, {0, 3}, {0, 2}, {0, 1}};

  int num_cells = cells.size()/4;

  std::vector<double> dh_angles;
  dh_angles.resize(num_cells*6);

  for (unsigned int c = 0; c < num_cells; ++c)
  {
    for (unsigned int i = 0; i < 6; ++i)
    {
      const std::size_t i0 = cells[4*c + edges[i][0]];
      const std::size_t i1 = cells[4*c + edges[i][1]];
      const std::size_t i2 = cells[4*c + edges[5 - i][0]];
      const std::size_t i3 = cells[4*c + edges[5 - i][1]];

      std::vector<double> p0(points.begin()+i0*3,points.begin()+i0*3+3);
      std::vector<double> v1(points.begin()+i1*3,points.begin()+i1*3+3);
      std::vector<double> v2(points.begin()+i2*3,points.begin()+i2*3+3);
      std::vector<double> v3(points.begin()+i3*3,points.begin()+i3*3+3);

      // subtract p0 vector
      for (unsigned int j = 0; j < 3; ++j){
          v1[j] -= p0[j];
          v2[j] -= p0[j];
          v3[j] -= p0[j];
      }

      // normalize by l2-norm
      const double v1_div = l2_norm(v1);
      const double v2_div = l2_norm(v2);
      const double v3_div = l2_norm(v3);

      for (unsigned int j = 0; j < 3; ++j){
          v1[j] /= v1_div;
          v2[j] /= v2_div;
          v3[j] /= v3_div;
      }

      // compute dot product
      double v2Dotv3 = dot_product(v2,v3);
      double v1Dotv2 = dot_product(v1,v2);
      double v1Dotv3 = dot_product(v1,v3);

      // compute cross product
      std::vector<double> v1Crossv2 = cross_product(v1,v2);
      std::vector<double> v1Crossv3 = cross_product(v1,v3);

      const double norm_v1Crossv2 = l2_norm(v1Crossv2);
      const double norm_v1Crossv3 = l2_norm(v1Crossv3);

      double cphi = (v2Dotv3 - v1Dotv2 * v1Dotv3) / ((norm_v1Crossv2 * norm_v1Crossv3));
      dh_angles[c*6 + i] = std::acos(cphi); // in radians
    }
 }
return dh_angles;
}

py::array calc_dihedral_angles(py::array_t<double, py::array::c_style | py::array::forcecast> pts,
                    py::array_t<int, py::array::c_style | py::array::forcecast> cells)
{

  // check input dimensions
  int num_points = pts.size()/3;
  int num_cells = cells.size()/4;
  // allocate std::vector (to pass to the C++ function)
  std::vector<double> cppPts(num_points*3);
  std::vector<int> cppCells(num_cells*4);

  // copy py::array -> std::vector
  std::memcpy(cppPts.data(),pts.data(),3*num_points*sizeof(double));
  std::memcpy(cppCells.data(),cells.data(),4*num_cells*sizeof(int));

  std::vector<double> dh_angles = c_calc_dihedral_angles(cppPts, cppCells);

  ssize_t              num_angles = num_cells*6;
  ssize_t              ndim      = 2;
  ssize_t              sodble    = sizeof(double);
  std::vector<ssize_t> shape     = {num_angles, 1};
  std::vector<ssize_t> strides   = {sodble, sodble};

  // return 2-D NumPy array
  return py::array(py::buffer_info(
    dh_angles.data(),                           /* data as contiguous array  */
    sizeof(double),                          /* size of one scalar        */
    py::format_descriptor<double>::format(), /* data type                 */
    2,                                    /* number of dimensions      */
    shape,                                   /* shape of the matrix       */
    strides                                  /* strides for each axis     */
  ));
}


// fixed size calculation for 3x3 determinant
double c_calc_3x3determinant(std::vector<double> &m){
    // | (0,0) 0  (0,1) 1 (0,2) 2|
    // | (1,0) 3  (1,1) 4 (1,2) 5|
    // | (2,0) 6  (2,1) 7 (2,2) 8|
    return m[0]*(m[4]*m[8] - m[7]*m[5]) -
           m[1]*(m[3]*m[8] - m[6]*m[5]) +
           m[2]*(m[3]*m[7] - m[6]*m[4]);
}

// ----------------
// Python interface for c_calc_3x3determinant
// ----------------
double calc_3x3determinant(py::array_t<double, py::array::c_style | py::array::forcecast> matrix)
{
  // allocate std::vector (to pass to the C++ function)
  std::vector<double> cppMatrix(9);

  // copy py::array -> std::vector
  std::memcpy(cppMatrix.data(),matrix.data(),9*sizeof(double));

  // call cpp code
  double result= c_calc_3x3determinant(cppMatrix);

  return result;
}

// fixed size calculation for 4x4 determinant
double c_calc_4x4determinant(std::vector<double> &m) {
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

// ----------------
// Python interface for c_calc_4x4determinant
// ----------------
double calc_4x4determinant(py::array_t<double, py::array::c_style | py::array::forcecast> matrix)
{
  // allocate std::vector (to pass to the C++ function)
  std::vector<double> cppMatrix(16);

  // copy py::array -> std::vector
  std::memcpy(cppMatrix.data(),matrix.data(),16*sizeof(double));

  // call cpp code
  double result= c_calc_4x4determinant(cppMatrix);

  return result;
}


// Calcuate the gradient of the circumsphere radius wrt to point0
// used to guide the pertubation to the point
// !!!Assumes p0 is translated to (0,0,0)!!!
std::vector<double> c_calc_circumsphere_grad(std::vector<double> &p1, std::vector<double> &p2, std::vector<double> &p3)
{
    // coordinates of cell (p0 is assumed fixed)
    double x1, x2, x3;
    double y1, y2, y3;
    double z1, z2, z3;
    // the square of the point coordinates
    double p1sq;
    double p2sq;
    double p3sq;
    // some determinants
    double a;
    double Dx;
    double Dy;
    double Dz;
    // gradients of determinants
    std::vector<double> gradient_a;
    std::vector<double> gradient_Dx;
    std::vector<double> gradient_Dy;
    std::vector<double> gradient_Dz;
    // gradient of the cirumsphere radius
    std::vector<double> gradient_cradius;
    // components of the gradient of the circumsphere radius
    double gradient_cradius_x;
    double gradient_cradius_y;
    double gradient_cradius_z;
    // used for getting the determinants
    std::vector<double> tmp;
    // unpack coordinates
    x1 = p1[0];
    y1 = p1[1];
    z1 = p1[2];
    x2 = p2[0];
    y2 = p2[1];
    z2 = p2[2];
    x3 = p3[0];
    y3 = p3[1];
    z3 = p3[2];
    // Calculate squared coordinates for p1,p2,p3
    p1sq = p1[0]*p1[0]+p1[1]*p1[1]+p1[2]*p1[2];
    p2sq = p2[0]*p2[0]+p2[1]*p2[1]+p2[2]*p2[2];
    p3sq = p3[0]*p3[0]+p3[1]*p3[1]+p3[2]*p3[2];
    // Determinant of "a"
    tmp = {x1,y1,z1,x2,y2,z3,x3,y3,z3};
    a = c_calc_3x3determinant(tmp);
    // Determinant of "Dx"
    tmp = {p1sq, y1, z1, p2sq, y2, z2, p3sq, y3, z3};
    Dx = c_calc_3x3determinant(tmp);
    Dx *= -1;
    // Determinant of "Dy"
    tmp = {p1sq, x1, z1, p2sq, x2, z2, p3sq, x3, z3};
    Dy = c_calc_3x3determinant(tmp);
    // Determinant of "Dz"
    tmp = {p1sq, x1, y1, p2sq, x2, y2, p3sq, x3, y3};
    Dz = c_calc_3x3determinant(tmp);
    Dz *= -1;
    // Gradient of "a"
    gradient_a = {y2*z3-y3*z2, -(x2*z3-x3*z2), x2*y3-x3*y2};
    // Gradient of "Dx"
    gradient_Dx = {-2*x1*gradient_a[0], -2*y1*gradient_a[0]+p2sq*z3-p3sq*z2, -2*z1*gradient_a[0]-p2sq*y3+p3sq*y2};
    // Gradient of "Dy"
    gradient_Dy = {-2*x1*gradient_a[1]-p2sq*z3+p3sq*z2, -2*y1*gradient_a[1], -2*z1*gradient_a[1]+p2sq*x3-p3sq*x2};
    // Gradient of "Dz"
    gradient_Dz = {-2*x1*gradient_a[2]+p2sq*y3-p3sq*y2, -2*y1*gradient_a[2]-p2sq*x3+p3sq*x2, -2*z1*gradient_a[2]};
    // The gradient of the circumradius
    gradient_cradius_x = 1/(2*a*a*a)*(a*(Dx*gradient_Dx[0]+Dy*gradient_Dy[0]+Dz*gradient_Dz[0])-gradient_a[0]*(Dx*Dx+Dy*Dy+Dz*Dz));
    gradient_cradius_y = 1/(2*a*a*a)*(a*(Dx*gradient_Dx[1]+Dy*gradient_Dy[1]+Dz*gradient_Dz[1])-gradient_a[1]*(Dx*Dx+Dy*Dy+Dz*Dz));
    gradient_cradius_z = 1/(2*a*a*a)*(a*(Dx*gradient_Dx[2]+Dy*gradient_Dy[2]+Dz*gradient_Dz[2])-gradient_a[2]*(Dx*Dx+Dy*Dy+Dz*Dz));
    gradient_cradius = {gradient_cradius_x, gradient_cradius_y, gradient_cradius_z};
    return gradient_cradius;
}

// Python wrapper accepts the coordinates of the slivers' vertices
py::array calc_circumsphere_grad(py::array_t<double, py::array::c_style | py::array::forcecast> p0,
                                 py::array_t<double, py::array::c_style | py::array::forcecast> p1,
                                 py::array_t<double, py::array::c_style | py::array::forcecast> p2,
                                 py::array_t<double, py::array::c_style | py::array::forcecast> p3)
{

  int num_points = p0.size()/3;
  // allocate std::vector (to pass to the C++ function)
  std::vector<double> cppP0(3*num_points);
  std::vector<double> cppP1(3*num_points);
  std::vector<double> cppP2(3*num_points);
  std::vector<double> cppP3(3*num_points);

  // copy py::array -> std::vector
  std::memcpy(cppP0.data(),p0.data(),3*num_points*sizeof(double));
  std::memcpy(cppP1.data(),p1.data(),3*num_points*sizeof(double));
  std::memcpy(cppP2.data(),p2.data(),3*num_points*sizeof(double));
  std::memcpy(cppP3.data(),p3.data(),3*num_points*sizeof(double));

  std::vector<double> circumsphere_grad;
  circumsphere_grad.resize(num_points*3);
  for(int i=0; i < num_points; ++i){
      // translate points so that P0 is at (0.0,0.0,0.0)
      for(int j=0; j < 3; ++j){
          cppP1[i*3+j]-=cppP0[i*3+j];
          cppP2[i*3+j]-=cppP0[i*3+j];
          cppP3[i*3+j]-=cppP0[i*3+j];
      }
      std::vector<double> tmp = c_calc_circumsphere_grad(cppP1, cppP2, cppP3);
      // unack gradient for each sliver
      for(int j=0; j<3; ++j){
        circumsphere_grad[i*3+j] =tmp[j];
      }
  }

  ssize_t              sodble    = sizeof(double);
  std::vector<ssize_t> shape     = {num_points, 1};
  std::vector<ssize_t> strides   = {sodble, sodble};

  // return 2-D NumPy array
  return py::array(py::buffer_info(
    circumsphere_grad.data(),                           /* data as contiguous array  */
    sizeof(double),                          /* size of one scalar        */
    py::format_descriptor<double>::format(), /* data type                 */
    2,                                    /* number of dimensions      */
    shape,                                   /* shape of the matrix       */
    strides                                  /* strides for each axis     */
  ));
}


PYBIND11_MODULE(fast_geometry, m) {
    m.def("calc_circumsphere_grad", &calc_circumsphere_grad);
    m.def("calc_dihedral_angles", &calc_dihedral_angles);
    m.def("calc_4x4determinant", &calc_4x4determinant);
    m.def("calc_3x3determinant", &calc_3x3determinant);
}
