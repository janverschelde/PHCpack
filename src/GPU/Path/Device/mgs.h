// mgs.h contains the prototypes for the function to MGS kernels

template <class ComplexType, class RealType>
int GPU_MGS
 ( int neq, int nvr, int BS_QR=256 ) // default for double precision
/*
 * MGS method for large problems. */
