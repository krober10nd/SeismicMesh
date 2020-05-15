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

// calcuate the gradient of the circumsphere radius wrt to a point
//std::vector<double> c_calc_circumsphere_grad(std::vector<double> &p0, std::vector<double> &p1, std::vector<double> &p2, std::vector<double> &p3)
//{
//
//    return 0;
//}

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

PYBIND11_MODULE(fast_geometry, m) {
    m.def("calc_dihedral_angles", &calc_dihedral_angles);
    m.def("calc_4x4determinant", &calc_4x4determinant);
    m.def("calc_3x3determinant", &calc_3x3determinant);
}
