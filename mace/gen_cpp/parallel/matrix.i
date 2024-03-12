/* begin func_claim */
void Allreduce<eigen_type>__(const MPI_Comm &comm, const Eigen::<eigen_type> &source, Eigen::<eigen_type> &dest);
void InplaceAllreduce<eigen_type>__(const MPI_Comm &comm, Eigen::<eigen_type> &source);
void Bcast<eigen_type>__(const MPI_Comm &comm, Eigen::<eigen_type> &source, int root);
void Send<eigen_type>__(const MPI_Comm &comm, const Eigen::<eigen_type> &source, int destID, int tag);
void Recv<eigen_type>__(const MPI_Comm &comm, Eigen::<eigen_type> &dest, int sourceID, int tag);
/* end func_claim */

/* begin class_func_claim */
void Allreduce<eigen_type>(const Eigen::<eigen_type> &source, Eigen::<eigen_type> &dest) const;
void InplaceAllreduce<eigen_type>(Eigen::<eigen_type> &source) const;
void Bcast<eigen_type>(Eigen::<eigen_type> &source, int root = 0) const;
void Send<eigen_type>(const Eigen::<eigen_type> &source, int destID, int tag) const;
void Recv<eigen_type>(Eigen::<eigen_type> &dest, int sourceID, int tag) const;
/* end class_func_claim */

/* begin func_def */
/// @brief Reduce source matrix to dest matrix
/// @param[in] comm MPI communicator
/// @param[in] source source matrix
/// @param[out] dest dest matrix
/// @return None
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

/// @brief Reduce source matrix in-place
/// @param[in] comm MPI communicator
/// @param[in,out] source source matrix
/// @return None
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

/// @brief Broadcast source matrix from root process to other processes
/// @param[in] comm MPI communicator
/// @param[in,out] source source matrix
/// @param[in] root rank of root process in the communicator
/// @return None
void Bcast<eigen_type>__(const MPI_Comm &comm, Eigen::<eigen_type> &source, int root) {
#ifdef WITH_MPI
  int status = MPI_Bcast(source.data(), source.size(), <mpi_type>, root, comm);
  if (status != MPI_SUCCESS) {
    std::cerr << "MPI_Bcast call failed!\n";
    std::exit(-1);
  }
#endif
}

/// @brief Send source matrix to destination process
/// @param[in] comm MPI communicator
/// @param[in] source source matrix
/// @param[in] destID rank of destination process in the communicator
/// @param[in] tag tag to distinguish the send/recv process
/// @return None
void Send<eigen_type>__(const MPI_Comm &comm, const Eigen::<eigen_type> &source, int destID, int tag) {
#ifdef WITH_MPI
  int status = MPI_Send(source.data(), source.size(), <mpi_type>, destID, tag, comm);
  if (status != MPI_SUCCESS) {
    std::cerr << "MPI_Send call failed!\n";
    std::exit(-1);
  }
#endif
}

/// @brief Receive matrix from source process and save to destination
/// @param[in] comm MPI communicator
/// @param[out] dest destination matrix
/// @param[in] sourceID rank of source process in the communicator
/// @param[in] tag tag to distinguish the send/recv process
/// @return None
void Recv<eigen_type>__(const MPI_Comm &comm, Eigen::<eigen_type> &dest, int sourceID, int tag) {
#ifdef WITH_MPI
  int status = MPI_Recv(dest.data(), dest.size(), <mpi_type>, sourceID, tag, comm, MPI_STATUS_IGNORE);
  if (status != MPI_SUCCESS) {
    std::cerr << "MPI_Recv call failed!\n";
    std::exit(-1);
  }
#endif
}
/* end func_def */

/* begin class_func_def */
/// @brief Reduce source matrix to dest matrix
/// @param[in] source source matrix
/// @param[out] dest dest matrix
/// @return None
void BaseMPIEnv::Allreduce<eigen_type>(const Eigen::<eigen_type> &source,
                                   Eigen::<eigen_type> &dest) const {
  Allreduce<eigen_type>__(m_Comm, source, dest);
}

/// @brief Reduce source matrix in-place
/// @param[in,out] source source matrix
/// @return None
void BaseMPIEnv::InplaceAllreduce<eigen_type>(Eigen::<eigen_type> &source) const {
  InplaceAllreduce<eigen_type>__(m_Comm, source);
}

/// @brief Broadcast source matrix from root process to other processes
/// @param[in,out] source source matrix
/// @param[in] root rank of root process in the communicator
/// @return None
void BaseMPIEnv::Bcast<eigen_type>(Eigen::<eigen_type> &source, int root) const {
  Bcast<eigen_type>__(m_Comm, source, root);
}

/// @brief Send source matrix to destination process
/// @param[in] source source matrix
/// @param[in] destID rank of destination process in the communicator
/// @param[in] tag tag to distinguish the send/recv process
/// @return None
void BaseMPIEnv::Send<eigen_type>(const Eigen::<eigen_type> &source, int destID, int tag) const {
  Send<eigen_type>__(m_Comm, source, destID, tag);
}

/// @brief Receive matrix from source process and save to destination
/// @param[out] dest destination matrix
/// @param[in] sourceID rank of source process in the communicator
/// @param[in] tag tag to distinguish the send/recv process
/// @return None
void BaseMPIEnv::Recv<eigen_type>(Eigen::<eigen_type> &dest, int sourceID, int tag) const {
  Recv<eigen_type>__(m_Comm, dest, sourceID, tag);
}
/* end class_func_def */