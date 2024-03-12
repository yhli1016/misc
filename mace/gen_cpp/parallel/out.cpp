/* begin func_claim */
void AllreduceScalar__(const MPI_Comm &comm, const int &source, int &dest);
void InplaceAllreduceScalar__(const MPI_Comm &comm, int &source);
void BcastScalar__(const MPI_Comm &comm, int &source, int root);
void SendScalar__(const MPI_Comm &comm, const int &source, int destID, int tag);
void RecvScalar__(const MPI_Comm &comm, int &dest, int sourceID, int tag);
/* end func_claim */

/* begin func_claim */
void AllreduceScalar__(const MPI_Comm &comm, const double &source, double &dest);
void InplaceAllreduceScalar__(const MPI_Comm &comm, double &source);
void BcastScalar__(const MPI_Comm &comm, double &source, int root);
void SendScalar__(const MPI_Comm &comm, const double &source, int destID, int tag);
void RecvScalar__(const MPI_Comm &comm, double &dest, int sourceID, int tag);
/* end func_claim */

/* begin func_claim */
void AllreduceScalar__(const MPI_Comm &comm, const std::complex<double> &source, std::complex<double> &dest);
void InplaceAllreduceScalar__(const MPI_Comm &comm, std::complex<double> &source);
void BcastScalar__(const MPI_Comm &comm, std::complex<double> &source, int root);
void SendScalar__(const MPI_Comm &comm, const std::complex<double> &source, int destID, int tag);
void RecvScalar__(const MPI_Comm &comm, std::complex<double> &dest, int sourceID, int tag);
/* end func_claim */

/* begin class_func_claim */
void AllreduceScalar(const int &source, int &dest) const;
void InplaceAllreduceScalar(int &source) const;
void BcastScalar(int &source, int root = 0) const;
void SendScalar(const int &source, int destID, int tag) const;
void RecvScalar(int &dest, int sourceID, int tag) const;
/* end class_func_claim */

/* begin class_func_claim */
void AllreduceScalar(const double &source, double &dest) const;
void InplaceAllreduceScalar(double &source) const;
void BcastScalar(double &source, int root = 0) const;
void SendScalar(const double &source, int destID, int tag) const;
void RecvScalar(double &dest, int sourceID, int tag) const;
/* end class_func_claim */

/* begin class_func_claim */
void AllreduceScalar(const std::complex<double> &source, std::complex<double> &dest) const;
void InplaceAllreduceScalar(std::complex<double> &source) const;
void BcastScalar(std::complex<double> &source, int root = 0) const;
void SendScalar(const std::complex<double> &source, int destID, int tag) const;
void RecvScalar(std::complex<double> &dest, int sourceID, int tag) const;
/* end class_func_claim */

/* begin func_def */
/// @brief Reduce source scalar to dest scalar
/// @param[in] comm MPI communicator
/// @param[in] source source scalar
/// @param[out] dest dest scalar
/// @return None
void AllreduceScalar__(const MPI_Comm &comm, const int &source,
                              int &dest) {
#ifdef WITH_MPI
  int status = MPI_Allreduce(&source, &dest, 1, MPI_INT,
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
void InplaceAllreduceScalar__(const MPI_Comm &comm, int &source) {
#ifdef WITH_MPI
  int status = MPI_Allreduce(MPI_IN_PLACE, &source, 1, MPI_INT,
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
void BcastScalar__(const MPI_Comm &comm, int &source, int root) {
#ifdef WITH_MPI
  int status = MPI_Bcast(&source, 1, MPI_INT, root, comm);
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
void SendScalar__(const MPI_Comm &comm, const int &source, int destID, int tag) {
#ifdef WITH_MPI
  int status = MPI_Send(&source, 1, MPI_INT, destID, tag, comm);
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
void RecvScalar__(const MPI_Comm &comm, int &dest, int sourceID, int tag) {
#ifdef WITH_MPI
  int status = MPI_Recv(&dest, 1, MPI_INT, sourceID, tag, comm, MPI_STATUS_IGNORE);
  if (status != MPI_SUCCESS) {
    std::cerr << "MPI_Recv call failed!\n";
    std::exit(-1);
  }
#endif
}
/* end func_def */

/* begin func_def */
/// @brief Reduce source scalar to dest scalar
/// @param[in] comm MPI communicator
/// @param[in] source source scalar
/// @param[out] dest dest scalar
/// @return None
void AllreduceScalar__(const MPI_Comm &comm, const double &source,
                              double &dest) {
#ifdef WITH_MPI
  int status = MPI_Allreduce(&source, &dest, 1, MPI_DOUBLE,
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
void InplaceAllreduceScalar__(const MPI_Comm &comm, double &source) {
#ifdef WITH_MPI
  int status = MPI_Allreduce(MPI_IN_PLACE, &source, 1, MPI_DOUBLE,
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
void BcastScalar__(const MPI_Comm &comm, double &source, int root) {
#ifdef WITH_MPI
  int status = MPI_Bcast(&source, 1, MPI_DOUBLE, root, comm);
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
void SendScalar__(const MPI_Comm &comm, const double &source, int destID, int tag) {
#ifdef WITH_MPI
  int status = MPI_Send(&source, 1, MPI_DOUBLE, destID, tag, comm);
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
void RecvScalar__(const MPI_Comm &comm, double &dest, int sourceID, int tag) {
#ifdef WITH_MPI
  int status = MPI_Recv(&dest, 1, MPI_DOUBLE, sourceID, tag, comm, MPI_STATUS_IGNORE);
  if (status != MPI_SUCCESS) {
    std::cerr << "MPI_Recv call failed!\n";
    std::exit(-1);
  }
#endif
}
/* end func_def */

/* begin func_def */
/// @brief Reduce source scalar to dest scalar
/// @param[in] comm MPI communicator
/// @param[in] source source scalar
/// @param[out] dest dest scalar
/// @return None
void AllreduceScalar__(const MPI_Comm &comm, const std::complex<double> &source,
                              std::complex<double> &dest) {
#ifdef WITH_MPI
  int status = MPI_Allreduce(&source, &dest, 1, MPI_COMPLEX,
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
void InplaceAllreduceScalar__(const MPI_Comm &comm, std::complex<double> &source) {
#ifdef WITH_MPI
  int status = MPI_Allreduce(MPI_IN_PLACE, &source, 1, MPI_COMPLEX,
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
void BcastScalar__(const MPI_Comm &comm, std::complex<double> &source, int root) {
#ifdef WITH_MPI
  int status = MPI_Bcast(&source, 1, MPI_COMPLEX, root, comm);
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
void SendScalar__(const MPI_Comm &comm, const std::complex<double> &source, int destID, int tag) {
#ifdef WITH_MPI
  int status = MPI_Send(&source, 1, MPI_COMPLEX, destID, tag, comm);
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
void RecvScalar__(const MPI_Comm &comm, std::complex<double> &dest, int sourceID, int tag) {
#ifdef WITH_MPI
  int status = MPI_Recv(&dest, 1, MPI_COMPLEX, sourceID, tag, comm, MPI_STATUS_IGNORE);
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
void BaseMPIEnv::AllreduceScalar(const int& source,
                                   int &dest) const {
  AllreduceScalar__(m_Comm, source, dest);
}

/// @brief Reduce source scalar in-place
/// @param[in,out] source source scalar
/// @return None
void BaseMPIEnv::InplaceAllreduceScalar(int &source) const {
  InplaceAllreduceScalar__(m_Comm, source);
}

/// @brief Broadcast source scalar from root process to other processes
/// @param[in,out] source source scalar
/// @param[in] root rank of root process in the communicator
/// @return None
void BaseMPIEnv::BcastScalar(int &source, int root) const {
  BcastScalar__(m_Comm, source, root);
}

/// @brief Send source scalar to destination process
/// @param[in] source source scalar
/// @param[in] destID rank of destination process in the communicator
/// @param[in] tag tag to distinguish the send/recv process
/// @return None
void BaseMPIEnv::SendScalar(const int &source, int destID, int tag) const {
  SendScalar__(m_Comm, source, destID, tag);
}

/// @brief Receive scalar from source process and save to destination
/// @param[out] dest destination scalar
/// @param[in] sourceID rank of source process in the communicator
/// @param[in] tag tag to distinguish the send/recv process
/// @return None
void BaseMPIEnv::RecvScalar(int &dest, int sourceID, int tag) const {
  RecvScalar__(m_Comm, dest, sourceID, tag);
}
/* end class_func_def */
/* begin class_func_def */
/// @brief Reduce source scalar to dest scalar
/// @param[in] source source scalar
/// @param[out] dest dest scalar
/// @return None
void BaseMPIEnv::AllreduceScalar(const double& source,
                                   double &dest) const {
  AllreduceScalar__(m_Comm, source, dest);
}

/// @brief Reduce source scalar in-place
/// @param[in,out] source source scalar
/// @return None
void BaseMPIEnv::InplaceAllreduceScalar(double &source) const {
  InplaceAllreduceScalar__(m_Comm, source);
}

/// @brief Broadcast source scalar from root process to other processes
/// @param[in,out] source source scalar
/// @param[in] root rank of root process in the communicator
/// @return None
void BaseMPIEnv::BcastScalar(double &source, int root) const {
  BcastScalar__(m_Comm, source, root);
}

/// @brief Send source scalar to destination process
/// @param[in] source source scalar
/// @param[in] destID rank of destination process in the communicator
/// @param[in] tag tag to distinguish the send/recv process
/// @return None
void BaseMPIEnv::SendScalar(const double &source, int destID, int tag) const {
  SendScalar__(m_Comm, source, destID, tag);
}

/// @brief Receive scalar from source process and save to destination
/// @param[out] dest destination scalar
/// @param[in] sourceID rank of source process in the communicator
/// @param[in] tag tag to distinguish the send/recv process
/// @return None
void BaseMPIEnv::RecvScalar(double &dest, int sourceID, int tag) const {
  RecvScalar__(m_Comm, dest, sourceID, tag);
}
/* end class_func_def */
/* begin class_func_def */
/// @brief Reduce source scalar to dest scalar
/// @param[in] source source scalar
/// @param[out] dest dest scalar
/// @return None
void BaseMPIEnv::AllreduceScalar(const std::complex<double>& source,
                                   std::complex<double> &dest) const {
  AllreduceScalar__(m_Comm, source, dest);
}

/// @brief Reduce source scalar in-place
/// @param[in,out] source source scalar
/// @return None
void BaseMPIEnv::InplaceAllreduceScalar(std::complex<double> &source) const {
  InplaceAllreduceScalar__(m_Comm, source);
}

/// @brief Broadcast source scalar from root process to other processes
/// @param[in,out] source source scalar
/// @param[in] root rank of root process in the communicator
/// @return None
void BaseMPIEnv::BcastScalar(std::complex<double> &source, int root) const {
  BcastScalar__(m_Comm, source, root);
}

/// @brief Send source scalar to destination process
/// @param[in] source source scalar
/// @param[in] destID rank of destination process in the communicator
/// @param[in] tag tag to distinguish the send/recv process
/// @return None
void BaseMPIEnv::SendScalar(const std::complex<double> &source, int destID, int tag) const {
  SendScalar__(m_Comm, source, destID, tag);
}

/// @brief Receive scalar from source process and save to destination
/// @param[out] dest destination scalar
/// @param[in] sourceID rank of source process in the communicator
/// @param[in] tag tag to distinguish the send/recv process
/// @return None
void BaseMPIEnv::RecvScalar(std::complex<double> &dest, int sourceID, int tag) const {
  RecvScalar__(m_Comm, dest, sourceID, tag);
}
/* end class_func_def */
/* begin func_claim */
void AllreduceMatrixXi__(const MPI_Comm &comm, const Eigen::MatrixXi &source, Eigen::MatrixXi &dest);
void InplaceAllreduceMatrixXi__(const MPI_Comm &comm, Eigen::MatrixXi &source);
void BcastMatrixXi__(const MPI_Comm &comm, Eigen::MatrixXi &source, int root);
void SendMatrixXi__(const MPI_Comm &comm, const Eigen::MatrixXi &source, int destID, int tag);
void RecvMatrixXi__(const MPI_Comm &comm, Eigen::MatrixXi &dest, int sourceID, int tag);
/* end func_claim */

/* begin func_claim */
void AllreduceMatrixXd__(const MPI_Comm &comm, const Eigen::MatrixXd &source, Eigen::MatrixXd &dest);
void InplaceAllreduceMatrixXd__(const MPI_Comm &comm, Eigen::MatrixXd &source);
void BcastMatrixXd__(const MPI_Comm &comm, Eigen::MatrixXd &source, int root);
void SendMatrixXd__(const MPI_Comm &comm, const Eigen::MatrixXd &source, int destID, int tag);
void RecvMatrixXd__(const MPI_Comm &comm, Eigen::MatrixXd &dest, int sourceID, int tag);
/* end func_claim */

/* begin func_claim */
void AllreduceMatrixXcd__(const MPI_Comm &comm, const Eigen::MatrixXcd &source, Eigen::MatrixXcd &dest);
void InplaceAllreduceMatrixXcd__(const MPI_Comm &comm, Eigen::MatrixXcd &source);
void BcastMatrixXcd__(const MPI_Comm &comm, Eigen::MatrixXcd &source, int root);
void SendMatrixXcd__(const MPI_Comm &comm, const Eigen::MatrixXcd &source, int destID, int tag);
void RecvMatrixXcd__(const MPI_Comm &comm, Eigen::MatrixXcd &dest, int sourceID, int tag);
/* end func_claim */

/* begin class_func_claim */
void AllreduceMatrixXi(const Eigen::MatrixXi &source, Eigen::MatrixXi &dest) const;
void InplaceAllreduceMatrixXi(Eigen::MatrixXi &source) const;
void BcastMatrixXi(Eigen::MatrixXi &source, int root = 0) const;
void SendMatrixXi(const Eigen::MatrixXi &source, int destID, int tag) const;
void RecvMatrixXi(Eigen::MatrixXi &dest, int sourceID, int tag) const;
/* end class_func_claim */

/* begin class_func_claim */
void AllreduceMatrixXd(const Eigen::MatrixXd &source, Eigen::MatrixXd &dest) const;
void InplaceAllreduceMatrixXd(Eigen::MatrixXd &source) const;
void BcastMatrixXd(Eigen::MatrixXd &source, int root = 0) const;
void SendMatrixXd(const Eigen::MatrixXd &source, int destID, int tag) const;
void RecvMatrixXd(Eigen::MatrixXd &dest, int sourceID, int tag) const;
/* end class_func_claim */

/* begin class_func_claim */
void AllreduceMatrixXcd(const Eigen::MatrixXcd &source, Eigen::MatrixXcd &dest) const;
void InplaceAllreduceMatrixXcd(Eigen::MatrixXcd &source) const;
void BcastMatrixXcd(Eigen::MatrixXcd &source, int root = 0) const;
void SendMatrixXcd(const Eigen::MatrixXcd &source, int destID, int tag) const;
void RecvMatrixXcd(Eigen::MatrixXcd &dest, int sourceID, int tag) const;
/* end class_func_claim */

/* begin func_def */
/// @brief Reduce source matrix to dest matrix
/// @param[in] comm MPI communicator
/// @param[in] source source matrix
/// @param[out] dest dest matrix
/// @return None
void AllreduceMatrixXi__(const MPI_Comm &comm, const Eigen::MatrixXi &source,
                         Eigen::MatrixXi &dest) {
  dest.resize(source.rows(), source.cols());
#ifdef WITH_MPI
  int status = MPI_Allreduce(source.data(), dest.data(), source.size(), MPI_INT,
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
void InplaceAllreduceMatrixXi__(const MPI_Comm &comm, Eigen::MatrixXi &source) {
#ifdef WITH_MPI
  int status = MPI_Allreduce(MPI_IN_PLACE, source.data(), source.size(), MPI_INT,
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
void BcastMatrixXi__(const MPI_Comm &comm, Eigen::MatrixXi &source, int root) {
#ifdef WITH_MPI
  int status = MPI_Bcast(source.data(), source.size(), MPI_INT, root, comm);
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
void SendMatrixXi__(const MPI_Comm &comm, const Eigen::MatrixXi &source, int destID, int tag) {
#ifdef WITH_MPI
  int status = MPI_Send(source.data(), source.size(), MPI_INT, destID, tag, comm);
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
void RecvMatrixXi__(const MPI_Comm &comm, Eigen::MatrixXi &dest, int sourceID, int tag) {
#ifdef WITH_MPI
  int status = MPI_Recv(dest.data(), dest.size(), MPI_INT, sourceID, tag, comm, MPI_STATUS_IGNORE);
  if (status != MPI_SUCCESS) {
    std::cerr << "MPI_Recv call failed!\n";
    std::exit(-1);
  }
#endif
}
/* end func_def */

/* begin func_def */
/// @brief Reduce source matrix to dest matrix
/// @param[in] comm MPI communicator
/// @param[in] source source matrix
/// @param[out] dest dest matrix
/// @return None
void AllreduceMatrixXd__(const MPI_Comm &comm, const Eigen::MatrixXd &source,
                         Eigen::MatrixXd &dest) {
  dest.resize(source.rows(), source.cols());
#ifdef WITH_MPI
  int status = MPI_Allreduce(source.data(), dest.data(), source.size(), MPI_DOUBLE,
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
void InplaceAllreduceMatrixXd__(const MPI_Comm &comm, Eigen::MatrixXd &source) {
#ifdef WITH_MPI
  int status = MPI_Allreduce(MPI_IN_PLACE, source.data(), source.size(), MPI_DOUBLE,
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
void BcastMatrixXd__(const MPI_Comm &comm, Eigen::MatrixXd &source, int root) {
#ifdef WITH_MPI
  int status = MPI_Bcast(source.data(), source.size(), MPI_DOUBLE, root, comm);
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
void SendMatrixXd__(const MPI_Comm &comm, const Eigen::MatrixXd &source, int destID, int tag) {
#ifdef WITH_MPI
  int status = MPI_Send(source.data(), source.size(), MPI_DOUBLE, destID, tag, comm);
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
void RecvMatrixXd__(const MPI_Comm &comm, Eigen::MatrixXd &dest, int sourceID, int tag) {
#ifdef WITH_MPI
  int status = MPI_Recv(dest.data(), dest.size(), MPI_DOUBLE, sourceID, tag, comm, MPI_STATUS_IGNORE);
  if (status != MPI_SUCCESS) {
    std::cerr << "MPI_Recv call failed!\n";
    std::exit(-1);
  }
#endif
}
/* end func_def */

/* begin func_def */
/// @brief Reduce source matrix to dest matrix
/// @param[in] comm MPI communicator
/// @param[in] source source matrix
/// @param[out] dest dest matrix
/// @return None
void AllreduceMatrixXcd__(const MPI_Comm &comm, const Eigen::MatrixXcd &source,
                         Eigen::MatrixXcd &dest) {
  dest.resize(source.rows(), source.cols());
#ifdef WITH_MPI
  int status = MPI_Allreduce(source.data(), dest.data(), source.size(), MPI_COMPLEX,
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
void InplaceAllreduceMatrixXcd__(const MPI_Comm &comm, Eigen::MatrixXcd &source) {
#ifdef WITH_MPI
  int status = MPI_Allreduce(MPI_IN_PLACE, source.data(), source.size(), MPI_COMPLEX,
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
void BcastMatrixXcd__(const MPI_Comm &comm, Eigen::MatrixXcd &source, int root) {
#ifdef WITH_MPI
  int status = MPI_Bcast(source.data(), source.size(), MPI_COMPLEX, root, comm);
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
void SendMatrixXcd__(const MPI_Comm &comm, const Eigen::MatrixXcd &source, int destID, int tag) {
#ifdef WITH_MPI
  int status = MPI_Send(source.data(), source.size(), MPI_COMPLEX, destID, tag, comm);
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
void RecvMatrixXcd__(const MPI_Comm &comm, Eigen::MatrixXcd &dest, int sourceID, int tag) {
#ifdef WITH_MPI
  int status = MPI_Recv(dest.data(), dest.size(), MPI_COMPLEX, sourceID, tag, comm, MPI_STATUS_IGNORE);
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
void BaseMPIEnv::AllreduceMatrixXi(const Eigen::MatrixXi &source,
                                   Eigen::MatrixXi &dest) const {
  AllreduceMatrixXi__(m_Comm, source, dest);
}

/// @brief Reduce source matrix in-place
/// @param[in,out] source source matrix
/// @return None
void BaseMPIEnv::InplaceAllreduceMatrixXi(Eigen::MatrixXi &source) const {
  InplaceAllreduceMatrixXi__(m_Comm, source);
}

/// @brief Broadcast source matrix from root process to other processes
/// @param[in,out] source source matrix
/// @param[in] root rank of root process in the communicator
/// @return None
void BaseMPIEnv::BcastMatrixXi(Eigen::MatrixXi &source, int root) const {
  BcastMatrixXi__(m_Comm, source, root);
}

/// @brief Send source matrix to destination process
/// @param[in] source source matrix
/// @param[in] destID rank of destination process in the communicator
/// @param[in] tag tag to distinguish the send/recv process
/// @return None
void BaseMPIEnv::SendMatrixXi(const Eigen::MatrixXi &source, int destID, int tag) const {
  SendMatrixXi__(m_Comm, source, destID, tag);
}

/// @brief Receive matrix from source process and save to destination
/// @param[out] dest destination matrix
/// @param[in] sourceID rank of source process in the communicator
/// @param[in] tag tag to distinguish the send/recv process
/// @return None
void BaseMPIEnv::RecvMatrixXi(Eigen::MatrixXi &dest, int sourceID, int tag) const {
  RecvMatrixXi__(m_Comm, dest, sourceID, tag);
}
/* end class_func_def */
/* begin class_func_def */
/// @brief Reduce source matrix to dest matrix
/// @param[in] source source matrix
/// @param[out] dest dest matrix
/// @return None
void BaseMPIEnv::AllreduceMatrixXd(const Eigen::MatrixXd &source,
                                   Eigen::MatrixXd &dest) const {
  AllreduceMatrixXd__(m_Comm, source, dest);
}

/// @brief Reduce source matrix in-place
/// @param[in,out] source source matrix
/// @return None
void BaseMPIEnv::InplaceAllreduceMatrixXd(Eigen::MatrixXd &source) const {
  InplaceAllreduceMatrixXd__(m_Comm, source);
}

/// @brief Broadcast source matrix from root process to other processes
/// @param[in,out] source source matrix
/// @param[in] root rank of root process in the communicator
/// @return None
void BaseMPIEnv::BcastMatrixXd(Eigen::MatrixXd &source, int root) const {
  BcastMatrixXd__(m_Comm, source, root);
}

/// @brief Send source matrix to destination process
/// @param[in] source source matrix
/// @param[in] destID rank of destination process in the communicator
/// @param[in] tag tag to distinguish the send/recv process
/// @return None
void BaseMPIEnv::SendMatrixXd(const Eigen::MatrixXd &source, int destID, int tag) const {
  SendMatrixXd__(m_Comm, source, destID, tag);
}

/// @brief Receive matrix from source process and save to destination
/// @param[out] dest destination matrix
/// @param[in] sourceID rank of source process in the communicator
/// @param[in] tag tag to distinguish the send/recv process
/// @return None
void BaseMPIEnv::RecvMatrixXd(Eigen::MatrixXd &dest, int sourceID, int tag) const {
  RecvMatrixXd__(m_Comm, dest, sourceID, tag);
}
/* end class_func_def */
/* begin class_func_def */
/// @brief Reduce source matrix to dest matrix
/// @param[in] source source matrix
/// @param[out] dest dest matrix
/// @return None
void BaseMPIEnv::AllreduceMatrixXcd(const Eigen::MatrixXcd &source,
                                   Eigen::MatrixXcd &dest) const {
  AllreduceMatrixXcd__(m_Comm, source, dest);
}

/// @brief Reduce source matrix in-place
/// @param[in,out] source source matrix
/// @return None
void BaseMPIEnv::InplaceAllreduceMatrixXcd(Eigen::MatrixXcd &source) const {
  InplaceAllreduceMatrixXcd__(m_Comm, source);
}

/// @brief Broadcast source matrix from root process to other processes
/// @param[in,out] source source matrix
/// @param[in] root rank of root process in the communicator
/// @return None
void BaseMPIEnv::BcastMatrixXcd(Eigen::MatrixXcd &source, int root) const {
  BcastMatrixXcd__(m_Comm, source, root);
}

/// @brief Send source matrix to destination process
/// @param[in] source source matrix
/// @param[in] destID rank of destination process in the communicator
/// @param[in] tag tag to distinguish the send/recv process
/// @return None
void BaseMPIEnv::SendMatrixXcd(const Eigen::MatrixXcd &source, int destID, int tag) const {
  SendMatrixXcd__(m_Comm, source, destID, tag);
}

/// @brief Receive matrix from source process and save to destination
/// @param[out] dest destination matrix
/// @param[in] sourceID rank of source process in the communicator
/// @param[in] tag tag to distinguish the send/recv process
/// @return None
void BaseMPIEnv::RecvMatrixXcd(Eigen::MatrixXcd &dest, int sourceID, int tag) const {
  RecvMatrixXcd__(m_Comm, dest, sourceID, tag);
}
/* end class_func_def */
