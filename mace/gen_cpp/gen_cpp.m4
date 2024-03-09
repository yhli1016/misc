changequote([,])dnl

define([FUNC_CLAIM], [
void Allreduce$1__(const MPI_Comm &comm, const Eigen::$1 &source,
                         Eigen::$1 &dest);
])dnl

define([FUNC_DEF], [
void Allreduce$1__(const MPI_Comm &comm, const Eigen::$1 &source,
                         Eigen::$1 &dest) {
  dest.resize(source.rows(), source.cols());
#ifdef WITH_MPI
  int status = MPI_Allreduce(source.data(), dest.data(), source.size(), $2,
                             MPI_SUM, comm);
  if (status != 0) {
    std::cerr << "MPI_Allreduce call failed!\n";
    std::exit(-1);
  }
#else
  dest = source;
#endif
}
])dnl

define([BASE_CLAIM], [
void Allreduce$1(const Eigen::$1 &source, Eigen::$1 &dest) const;
])dnl

define([BASE_DEF], [
void BaseMPIEnv::Allreduce$1(const Eigen::$1 &source,
                                   Eigen::$1 &dest) const {
  Allreduce$1__(m_Comm, source, dest);
}
])dnl

FUNC_CLAIM([MatrixXi], [MPI_INT])
FUNC_CLAIM([MatrixXd], [MPI_DOUBLE])
FUNC_CLAIM([MatrixXcd], [MPI_COMPLEX])
FUNC_DEF([MatrixXi], [MPI_INT])
FUNC_DEF([MatrixXd], [MPI_DOUBLE])
FUNC_DEF([MatrixXcd], [MPI_COMPLEX])
BASE_CLAIM([MatrixXi], [MPI_INT])
BASE_CLAIM([MatrixXd], [MPI_DOUBLE])
BASE_CLAIM([MatrixXcd], [MPI_COMPLEX])
BASE_DEF([MatrixXi], [MPI_INT])
BASE_DEF([MatrixXd], [MPI_DOUBLE])
BASE_DEF([MatrixXcd], [MPI_COMPLEX])
