/* begin func_claim */
void Allreduce<eigen_type>__(const MPI_Comm &comm, const Eigen::<eigen_type> &source,
                         Eigen::<eigen_type> &dest);
void InplaceAllreduce<eigen_type>__(const MPI_Comm &comm, Eigen::<eigen_type> &source);
void Bcast<eigen_type>__(const MPI_Comm &comm, Eigen::<eigen_type> &source, int root);
/* end func_claim */

/* begin class_func_claim */
void Allreduce<eigen_type>(const Eigen::<eigen_type> &source, Eigen::<eigen_type> &dest) const;
void InplaceAllreduce<eigen_type>(Eigen::<eigen_type> &source) const;
void Bcast<eigen_type>(Eigen::<eigen_type> &source, int root = 0) const;
/* end class_func_claim */

/* begin func_def */
void Allreduce<eigen_type>__(const MPI_Comm &comm, const Eigen::<eigen_type> &source,
                         Eigen::<eigen_type> &dest) {
  dest.resize(source.rows(), source.cols());
#ifdef WITH_MPI
  int status = MPI_Allreduce(source.data(), dest.data(), source.size(), <mpi_type>,
                             MPI_SUM, comm);
  if (status != MPI_SUCCESS) {
    std::cerr << "MPI_Allreduce call failed!\n";
    std::exit(-1);
  }
#else
  dest = source;
#endif
}

void InplaceAllreduce<eigen_type>__(const MPI_Comm &comm, Eigen::<eigen_type> &source) {
#ifdef WITH_MPI
  int status = MPI_Allreduce(MPI_IN_PLACE, source.data(), source.size(), <mpi_type>,
                             MPI_SUM, comm);
  if (status != MPI_SUCCESS) {
    std::cerr << "MPI_Allreduce call failed!\n";
    std::exit(-1);
  }
#endif
}

void Bcast<eigen_type>__(const MPI_Comm &comm, Eigen::<eigen_type> &source, int root) {
#ifdef WITH_MPI
  int status = MPI_Bcast(source.data(), source.size(), <mpi_type>, root, comm);
  if (status != MPI_SUCCESS) {
    std::cerr << "MPI_Bcast call failed!\n";
    std::exit(-1);
  }
#endif
}
/* end func_def */

/* begin class_func_def */
void BaseMPIEnv::Allreduce<eigen_type>(const Eigen::<eigen_type> &source,
                                   Eigen::<eigen_type> &dest) const {
  Allreduce<eigen_type>__(m_Comm, source, dest);
}

void BaseMPIEnv::InplaceAllreduce<eigen_type>(Eigen::<eigen_type> &source) const {
  InplaceAllreduce<eigen_type>__(m_Comm, source);
}

void BaseMPIEnv::Bcast<eigen_type>(Eigen::<eigen_type> &source, int root) const {
  Bcast<eigen_type>__(m_Comm, source, root);
}
/* end class_func_def */