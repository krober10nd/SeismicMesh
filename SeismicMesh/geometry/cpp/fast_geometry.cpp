#include <vector>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

double l2_norm(const std::array<double, 3> &u) {
    double accum = 0.;
    for (int i = 0; i < 3; ++i) {
        accum += u[i] * u[i];
    }
    return sqrt(accum);
}

std::array<double, 3> cross_product(const std::array<double, 3> &vect_A, const std::array<double, 3> &vect_B)
{
    std::array<double,3> cross_P;
    cross_P[0] = vect_A[1] * vect_B[2] - vect_A[2] * vect_B[1];
    cross_P[1] = vect_A[2] * vect_B[0] - vect_A[0] * vect_B[2];
    cross_P[2] = vect_A[0] * vect_B[1] - vect_A[1] * vect_B[0];
    return cross_P;
}

double dot_product(const std::array<double, 3> &vect_A, const std::array<double, 3> &vect_B)
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

  std::array<double,3> p0;
  std::array<double,3> v1;
  std::array<double,3> v2;
  std::array<double,3> v3;

  for (unsigned int c = 0; c < num_cells; ++c)
  {
    for (unsigned int i = 0; i < 6; ++i)
    {
      const std::size_t i0 = cells[4*c + edges[i][0]];
      const std::size_t i1 = cells[4*c + edges[i][1]];
      const std::size_t i2 = cells[4*c + edges[5 - i][0]];
      const std::size_t i3 = cells[4*c + edges[5 - i][1]];

      // populate point coordinates
      for (unsigned int j = 0; j < 3; ++j){
          p0[j] = points[i0*3+j];
          v1[j] = points[i1*3+j];
          v2[j] = points[i2*3+j];
          v3[j] = points[i3*3+j];
      }

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
      std::array<double, 3> v1Crossv2 = cross_product(v1,v2);
      std::array<double, 3> v1Crossv3 = cross_product(v1,v3);

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

//Calculate the gradient of the tetrahedral's volume wrt to point0
std::vector<double> c_calc_volume_grad(std::vector<double> &p1, std::vector<double> &p2, std::vector<double> &p3){
    int num_points = p1.size()/3;
    double x1, x2, x3;
    double y1, y2, y3;
    double z1, z2, z3;
    std::vector<double> gradient;
    gradient.resize(3*num_points);
    for(int i=0; i < num_points; ++i){
        // unpack coordinates
        x1 = p1[3*i];
        y1 = p1[3*i+1];
        z1 = p1[3*i+2];
        x2 = p2[3*i];
        y2 = p2[3*i+1];
        z2 = p2[3*i+2];
        x3 = p3[3*i];
        y3 = p3[3*i+1];
        z3 = p3[3*i+2];
        //
        gradient[3*i]=(1.0/6.0)*(y2*z3+y1*(z2-z3)-y3*z2-z1*(y2-y3));
        gradient[3*i+1]=(1.0/6.0)*(-x2*z3-x1*(z2-z3)+x3*z2+z1*(x2-x3));
        gradient[3*i+2]=(1.0/6.0)*(x2*y3+x1*(y2-y3)-x3*y2-y1*(x2-x3));
    }
    return gradient;
}

// Python wrapper accepts the coordinates of the slivers' vertices
py::array calc_volume_grad(py::array_t<double, py::array::c_style | py::array::forcecast> p1,
                           py::array_t<double, py::array::c_style | py::array::forcecast> p2,
                           py::array_t<double, py::array::c_style | py::array::forcecast> p3)
{

  int num_points = p1.size()/3;

  std::vector<double> cppP1(3*num_points);
  std::vector<double> cppP2(3*num_points);
  std::vector<double> cppP3(3*num_points);

  std::memcpy(cppP1.data(),p1.data(),3*num_points*sizeof(double));
  std::memcpy(cppP2.data(),p2.data(),3*num_points*sizeof(double));
  std::memcpy(cppP3.data(),p3.data(),3*num_points*sizeof(double));

  std::vector<double> volume_grad = c_calc_volume_grad(cppP1,cppP2,cppP3);

  ssize_t              sodble    = sizeof(double);
  std::vector<ssize_t> shape     = {num_points, 3};
  std::vector<ssize_t> strides   = {sodble*3, sodble};

  // return 2-D NumPy array
  return py::array(py::buffer_info(
    volume_grad.data(),                           /* data as contiguous array  */
    sizeof(double),                          /* size of one scalar        */
    py::format_descriptor<double>::format(), /* data type                 */
    2,                                    /* number of dimensions      */
    shape,                                   /* shape of the matrix       */
    strides                                  /* strides for each axis     */
  ));
}

// used to guide the pertubation to the point
std::vector<double> c_calc_circumsphere_grad(std::vector<double> &p0, std::vector<double> &p1, std::vector<double> &p2, std::vector<double> &p3)
{
    const int num_points = p0.size()/3;

    // coordinates of cell
    double x1;
    double y1;
    double z1;

    double x2;
    double y2;
    double z2;

    double x3;
    double y3;
    double z3;

    std::vector<double>gradient_cradius;
    gradient_cradius.resize(num_points*3);

    for(int i=0; i < num_points; ++i){
        // translate the tet so that fourth point is the origin
        x1 = p0[3*i] - p3[3*i];
        y1 = p0[3*i+1] - p3[3*i+1];
        z1 = p0[3*i+2] - p3[3*i+2];

        x2 = p1[3*i] - p3[3*i];
        y2 = p1[3*i+1] - p3[3*i+1];
        z2 = p1[3*i+2] - p3[3*i+2];

        x3 = p2[3*i] - p3[3*i];
        y3 = p2[3*i+1] - p3[3*i+1];
        z3 = p2[3*i+2] - p3[3*i+2];

        // pre-compute everything
        double sq_p1 = x1*x1 + y1*y1 + z1*z1;
        double sq_p2 = x2*x2 + y2*y2 + z2*z2;
        double sq_p3 = x3*x3 + y3*y3 + z3*z3;

        // every derivative is computed w.r.t p1 (x1, y1, z1)
        double da_dx = y2*z3 - y3*z2;
        double da_dy = z2*x3 - x2*z3;
        double da_dz = x2*y3 - x3*y2;

        double dDx_dx = -2.0*x1*da_dx;
        double dDx_dy = -2.0*y1*da_dx + sq_p2*z3 - sq_p3*z2;
        double dDx_dz = -2.0*z1*da_dx - sq_p2*y3 + sq_p3*y2;

        double dDy_dx = -2.0*x1*da_dy - sq_p2*z3 + sq_p3*z2;
        double dDy_dy = -2.0*y1*da_dy;
        double dDy_dz = -2.0*z1*da_dy + sq_p2*x3 - sq_p3*x2;

        double dDz_dx = -2.0*x1*da_dz + sq_p2*y3 - sq_p3*y2;
        double dDz_dy = -2.0*y1*da_dz - sq_p2*x3 + sq_p3*x2;
        double dDz_dz = -2.0*z1*da_dz;

        double a  = x1*da_dx + y1*da_dy + z1*da_dz;
        //std::cout<<a<<std::endl;
        //if ( CGAL_NTS is_zero(a) )
        //  return CGAL::NULL_VECTOR;

        double Dx = -sq_p1*da_dx + y1*(sq_p2*z3 - sq_p3*z2) - z1*(sq_p2*y3 - sq_p3*y2);
        double Dy = -sq_p1*da_dy - x1*(sq_p2*z3 - sq_p3*z2) + z1*(sq_p2*x3 - sq_p3*x2);
        double Dz = -sq_p1*da_dz + x1*(sq_p2*y3 - sq_p3*y2) - y1*(sq_p2*x3 - sq_p3*x2);

         // compute gradient vector
        double sum_sqD = Dx*Dx + Dy*Dy + Dz*Dz;
        double gx = (Dx*dDx_dx + Dy*dDy_dx + Dz*dDz_dx) / (2.0*a*a) - (da_dx * sum_sqD) / (2.0*a*a*a);
        double gy = (Dx*dDx_dy + Dy*dDy_dy + Dz*dDz_dy) / (2.0*a*a) - (da_dy * sum_sqD) / (2.0*a*a*a);
        double gz = (Dx*dDx_dz + Dy*dDy_dz + Dz*dDz_dz) / (2.0*a*a) - (da_dz * sum_sqD) / (2.0*a*a*a);

        gradient_cradius[3*i]=gx;
        gradient_cradius[3*i+1]=gy;
        gradient_cradius[3*i+2]=gz;
    }

    // the square of the point coordinates
    return gradient_cradius;
}

// Python wrapper accepts the coordinates of the slivers' vertices
py::array calc_circumsphere_grad(py::array_t<double, py::array::c_style | py::array::forcecast> p0,
                                 py::array_t<double, py::array::c_style | py::array::forcecast> p1,
                                 py::array_t<double, py::array::c_style | py::array::forcecast> p2,
                                 py::array_t<double, py::array::c_style | py::array::forcecast> p3)
{

  int num_points = p0.size()/3;

  std::vector<double> cppP0(3*num_points);
  std::vector<double> cppP1(3*num_points);
  std::vector<double> cppP2(3*num_points);
  std::vector<double> cppP3(3*num_points);

  std::memcpy(cppP0.data(),p0.data(),3*num_points*sizeof(double));
  std::memcpy(cppP1.data(),p1.data(),3*num_points*sizeof(double));
  std::memcpy(cppP2.data(),p2.data(),3*num_points*sizeof(double));
  std::memcpy(cppP3.data(),p3.data(),3*num_points*sizeof(double));

  std::vector<double> circumsphere_grad= c_calc_circumsphere_grad(cppP0, cppP1, cppP2, cppP3);

  ssize_t              sodble    = sizeof(double);
  std::vector<ssize_t> shape     = {num_points, 3};
  std::vector<ssize_t> strides   = {sodble*3, sodble};

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
    m.def("calc_volume_grad", &calc_volume_grad);
    m.def("calc_circumsphere_grad", &calc_circumsphere_grad);
    m.def("calc_dihedral_angles", &calc_dihedral_angles);
    m.def("calc_4x4determinant", &calc_4x4determinant);
    m.def("calc_3x3determinant", &calc_3x3determinant);
}
