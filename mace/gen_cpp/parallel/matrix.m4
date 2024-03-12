changequote([,])dnl

define([FUNC_CLAIM], [
void Allreduce$1__(const MPI_Comm &comm, const Eigen::$1 &source, Eigen::$1 &dest);
void InplaceAllreduce$1__(const MPI_Comm &comm, Eigen::$1 &source);
void Bcast$1__(const MPI_Comm &comm, Eigen::$1 &source, int root);
void Send$1__(const MPI_Comm &comm, const Eigen::$1 &source, int destID, int tag);
void Recv$1__(const MPI_Comm &comm, Eigen::$1 &dest, int sourceID, int tag);
])dnl

define([BASE_CLAIM], [
void Allreduce$1(const Eigen::$1 &source, Eigen::$1 &dest) const;
void InplaceAllreduce$1(Eigen::$1 &source) const;
void Bcast$1(Eigen::$1 &source, int root = 0) const;
void Send$1(const Eigen::$1 &source, int destID, int tag) const;
void Recv$1(Eigen::$1 &dest, int sourceID, int tag) const;
])dnl

define([FUNC_DEF], [
void Allreduce$1__(const MPI_Comm &comm, const Eigen::$1 &source,
                         Eigen::$1 &dest) {
  dest.resize(source.rows(), source.cols());
#ifdef WITH_MPI
  int status = MPI_Allreduce(source.data(), dest.data(), source.size(), $2,
                             MPI_SUM, comm);
  if (status != MPI_SUCCESS) {
    std::cerr << "MPI_Allreduce call failed!\n";
    std::exit(-1);
  }
#else
  dest = source;
#endif
}

void InplaceAllreduce$1__(const MPI_Comm &comm, Eigen::$1 &source) {
#ifdef WITH_MPI
  int status = MPI_Allreduce(MPI_IN_PLACE, source.data(), source.size(), $2,
                             MPI_SUM, comm);
  if (status != MPI_SUCCESS) {
    std::cerr << "MPI_Allreduce call failed!\n";
    std::exit(-1);
  }
#endif
}

void Bcast$1__(const MPI_Comm &comm, Eigen::$1 &source, int root) {
#ifdef WITH_MPI
  int status = MPI_Bcast(source.data(), source.size(), $2, root, comm);
  if (status != MPI_SUCCESS) {
    std::cerr << "MPI_Bcast call failed!\n";
    std::exit(-1);
  }
#endif
}

void Send$1__(const MPI_Comm &comm, const Eigen::$1 &source, int destID, int tag) {
#ifdef WITH_MPI
  int status = MPI_Send(source.data(), source.size(), $2, destID, tag, comm);
  if (status != MPI_SUCCESS) {
    std::cerr << "MPI_Send call failed!\n";
    std::exit(-1);
  }
#endif
}

void Recv$1__(const MPI_Comm &comm, Eigen::$1 &dest, int sourceID, int tag) {
#ifdef WITH_MPI
  int status = MPI_Recv(dest.data(), dest.size(), $2, sourceID, tag, comm, MPI_STATUS_IGNORE);
  if (status != MPI_SUCCESS) {
    std::cerr << "MPI_Recv call failed!\n";
    std::exit(-1);
  }
#endif
}
])dnl

define([BASE_DEF], [
void BaseMPIEnv::Allreduce$1(const Eigen::$1 &source,
                                   Eigen::$1 &dest) const {
  Allreduce$1__(m_Comm, source, dest);
}

void BaseMPIEnv::InplaceAllreduce$1(Eigen::$1 &source) const {
  InplaceAllreduce$1__(m_Comm, source);
}

void BaseMPIEnv::Bcast$1(Eigen::$1 &source, int root) const {
  Bcast$1__(m_Comm, source, root);
}

void BaseMPIEnv::Send$1(const Eigen::$1 &source, int destID, int tag) const {
  Send$1__(m_Comm, source, destID, tag);
}

void BaseMPIEnv::Recv$1(Eigen::$1 &dest, int sourceID, int tag) const {
  Recv$1__(m_Comm, dest, sourceID, tag);
}
])dnl

FUNC_CLAIM([MatrixXi], [MPI_INT])
FUNC_CLAIM([MatrixXd], [MPI_DOUBLE])
FUNC_CLAIM([MatrixXcd], [MPI_COMPLEX])
BASE_CLAIM([MatrixXi], [MPI_INT])
BASE_CLAIM([MatrixXd], [MPI_DOUBLE])
BASE_CLAIM([MatrixXcd], [MPI_COMPLEX])
FUNC_DEF([MatrixXi], [MPI_INT])
FUNC_DEF([MatrixXd], [MPI_DOUBLE])
FUNC_DEF([MatrixXcd], [MPI_COMPLEX])
BASE_DEF([MatrixXi], [MPI_INT])
BASE_DEF([MatrixXd], [MPI_DOUBLE])
BASE_DEF([MatrixXcd], [MPI_COMPLEX])
