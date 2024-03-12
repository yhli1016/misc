changequote([,])dnl

define([FUNC_CLAIM], [
void AllreduceScalar__(const MPI_Comm &comm, const $1 &source, $1 &dest);
void InplaceAllreduceScalar__(const MPI_Comm &comm, $1 &source);
void BcastScalar__(const MPI_Comm &comm, $1 &source, int root);
void SendScalar__(const MPI_Comm &comm, const $1 &source, int destID, int tag);
void RecvScalar__(const MPI_Comm &comm, $1 &dest, int sourceID, int tag);
])dnl

define([BASE_CLAIM], [
void AllreduceScalar(const $1 &source, $1 &dest) const;
void InplaceAllreduceScalar($1 &source) const;
void BcastScalar($1 &source, int root = 0) const;
void SendScalar(const $1 &source, int destID, int tag) const;
void RecvScalar($1 &dest, int sourceID, int tag) const;
])dnl

define([FUNC_DEF], [
void AllreduceScalar__(const MPI_Comm &comm, const $1 &source,
                              $1 &dest) {
#ifdef WITH_MPI
  int status = MPI_Allreduce(&source, &dest, 1, $2,
                             MPI_SUM, comm);
  if (status != MPI_SUCCESS) {
    std::cerr << "MPI_Allreduce call failed!\n";
    std::exit(-1);
  }
#else
  dest = source;
#endif
}

void InplaceAllreduceScalar__(const MPI_Comm &comm, $1 &source) {
#ifdef WITH_MPI
  int status = MPI_Allreduce(MPI_IN_PLACE, &source, 1, $2,
                             MPI_SUM, comm);
  if (status != MPI_SUCCESS) {
    std::cerr << "MPI_Allreduce call failed!\n";
    std::exit(-1);
  }
#endif
}

void BcastScalar__(const MPI_Comm &comm, $1 &source, int root) {
#ifdef WITH_MPI
  int status = MPI_Bcast(&source, 1, $2, root, comm);
  if (status != MPI_SUCCESS) {
    std::cerr << "MPI_Bcast call failed!\n";
    std::exit(-1);
  }
#endif
}

void SendScalar__(const MPI_Comm &comm, const $1 &source, int destID, int tag) {
#ifdef WITH_MPI
  int status = MPI_Send(&source, 1, $2, destID, tag, comm);
  if (status != MPI_SUCCESS) {
    std::cerr << "MPI_Send call failed!\n";
    std::exit(-1);
  }
#endif
}

void RecvScalar__(const MPI_Comm &comm, $1 &dest, int sourceID, int tag) {
#ifdef WITH_MPI
  int status = MPI_Recv(&dest, 1, $2, sourceID, tag, comm, MPI_STATUS_IGNORE);
  if (status != MPI_SUCCESS) {
    std::cerr << "MPI_Recv call failed!\n";
    std::exit(-1);
  }
#endif
}
])dnl

define([BASE_DEF], [
void BaseMPIEnv::AllreduceScalar(const $1& source,
                                   $1 &dest) const {
  AllreduceScalar__(m_Comm, source, dest);
}

void BaseMPIEnv::InplaceAllreduceScalar($1 &source) const {
  InplaceAllreduceScalar__(m_Comm, source);
}

void BaseMPIEnv::BcastScalar($1 &source, int root) const {
  BcastScalar__(m_Comm, source, root);
}

void BaseMPIEnv::SendScalar(const $1 &source, int destID, int tag) const {
  SendScalar__(m_Comm, source, destID, tag);
}

void BaseMPIEnv::RecvScalar($1 &dest, int sourceID, int tag) const {
  RecvScalar__(m_Comm, dest, sourceID, tag);
}
])dnl

FUNC_CLAIM([int], [MPI_INT])
FUNC_CLAIM([double], [MPI_DOUBLE])
FUNC_CLAIM([std::complex<double>], [MPI_COMPLEX])
BASE_CLAIM([int], [MPI_INT])
BASE_CLAIM([double], [MPI_DOUBLE])
BASE_CLAIM([std::complex<double>], [MPI_COMPLEX])
FUNC_DEF([int], [MPI_INT])
FUNC_DEF([double], [MPI_DOUBLE])
FUNC_DEF([std::complex<double>], [MPI_COMPLEX])
BASE_DEF([int], [MPI_INT])
BASE_DEF([double], [MPI_DOUBLE])
BASE_DEF([std::complex<double>], [MPI_COMPLEX])
