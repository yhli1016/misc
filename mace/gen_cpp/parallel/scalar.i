/* begin func_claim */
void AllreduceScalar__(const MPI_Comm &comm, const <scalar_type> &source, <scalar_type> &dest);
void InplaceAllreduceScalar__(const MPI_Comm &comm, <scalar_type> &source);
void BcastScalar__(const MPI_Comm &comm, <scalar_type> &source, int root);
void SendScalar__(const MPI_Comm &comm, const <scalar_type> &source, int destID, int tag);
void RecvScalar__(const MPI_Comm &comm, <scalar_type> &dest, int sourceID, int tag);
/* end func_claim */

/* begin class_func_claim */
void AllreduceScalar(const <scalar_type> &source, <scalar_type> &dest) const;
void InplaceAllreduceScalar(<scalar_type> &source) const;
void BcastScalar(<scalar_type> &source, int root = 0) const;
void SendScalar(const <scalar_type> &source, int destID, int tag) const;
void RecvScalar(<scalar_type> &dest, int sourceID, int tag) const;
/* end class_func_claim */

/* begin func_def */
/// @brief Reduce source scalar to dest scalar
/// @param[in] comm MPI communicator
/// @param[in] source source scalar
/// @param[out] dest dest scalar
/// @return None
void AllreduceScalar__(const MPI_Comm &comm, const <scalar_type> &source,
                              <scalar_type> &dest) {
#ifdef WITH_MPI
  int status = MPI_Allreduce(&source, &dest, 1, <mpi_type>,
                             MPI_SUM, comm);
  if (status != MPI_SUCCESS) {
    std::cerr << "MPI_Allreduce call failed!\n";
    std::exit(-1);
  }
#else
  dest = source;
#endif
}

/// @brief Reduce source scalar in-place
/// @param[in] comm MPI communicator
/// @param[in,out] source source scalar
/// @return None
void InplaceAllreduceScalar__(const MPI_Comm &comm, <scalar_type> &source) {
#ifdef WITH_MPI
  int status = MPI_Allreduce(MPI_IN_PLACE, &source, 1, <mpi_type>,
                             MPI_SUM, comm);
  if (status != MPI_SUCCESS) {
    std::cerr << "MPI_Allreduce call failed!\n";
    std::exit(-1);
  }
#endif
}

/// @brief Broadcast source scalar from root process to other processes
/// @param[in] comm MPI communicator
/// @param[in,out] source source scalar
/// @param[in] root rank of root process in the communicator
/// @return None
void BcastScalar__(const MPI_Comm &comm, <scalar_type> &source, int root) {
#ifdef WITH_MPI
  int status = MPI_Bcast(&source, 1, <mpi_type>, root, comm);
  if (status != MPI_SUCCESS) {
    std::cerr << "MPI_Bcast call failed!\n";
    std::exit(-1);
  }
#endif
}

/// @brief Send source scalar to destination process
/// @param[in] comm MPI communicator
/// @param[in] source source scalar
/// @param[in] destID rank of destination process in the communicator
/// @param[in] tag tag to distinguish the send/recv process
/// @return None
void SendScalar__(const MPI_Comm &comm, const <scalar_type> &source, int destID, int tag) {
#ifdef WITH_MPI
  int status = MPI_Send(&source, 1, <mpi_type>, destID, tag, comm);
  if (status != MPI_SUCCESS) {
    std::cerr << "MPI_Send call failed!\n";
    std::exit(-1);
  }
#endif
}

/// @brief Receive scalar from source process and save to destination
/// @param[in] comm MPI communicator
/// @param[out] dest destination scalar
/// @param[in] sourceID rank of source process in the communicator
/// @param[in] tag tag to distinguish the send/recv process
/// @return None
void RecvScalar__(const MPI_Comm &comm, <scalar_type> &dest, int sourceID, int tag) {
#ifdef WITH_MPI
  int status = MPI_Recv(&dest, 1, <mpi_type>, sourceID, tag, comm, MPI_STATUS_IGNORE);
  if (status != MPI_SUCCESS) {
    std::cerr << "MPI_Recv call failed!\n";
    std::exit(-1);
  }
#endif
}
/* end func_def */

/* begin class_func_def */
/// @brief Reduce source scalar to dest scalar
/// @param[in] source source scalar
/// @param[out] dest dest scalar
/// @return None
void BaseMPIEnv::AllreduceScalar(const <scalar_type>& source,
                                   <scalar_type> &dest) const {
  AllreduceScalar__(m_Comm, source, dest);
}

/// @brief Reduce source scalar in-place
/// @param[in,out] source source scalar
/// @return None
void BaseMPIEnv::InplaceAllreduceScalar(<scalar_type> &source) const {
  InplaceAllreduceScalar__(m_Comm, source);
}

/// @brief Broadcast source scalar from root process to other processes
/// @param[in,out] source source scalar
/// @param[in] root rank of root process in the communicator
/// @return None
void BaseMPIEnv::BcastScalar(<scalar_type> &source, int root) const {
  BcastScalar__(m_Comm, source, root);
}

/// @brief Send source scalar to destination process
/// @param[in] source source scalar
/// @param[in] destID rank of destination process in the communicator
/// @param[in] tag tag to distinguish the send/recv process
/// @return None
void BaseMPIEnv::SendScalar(const <scalar_type> &source, int destID, int tag) const {
  SendScalar__(m_Comm, source, destID, tag);
}

/// @brief Receive scalar from source process and save to destination
/// @param[out] dest destination scalar
/// @param[in] sourceID rank of source process in the communicator
/// @param[in] tag tag to distinguish the send/recv process
/// @return None
void BaseMPIEnv::RecvScalar(<scalar_type> &dest, int sourceID, int tag) const {
  RecvScalar__(m_Comm, dest, sourceID, tag);
}
/* end class_func_def */