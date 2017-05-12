
/**
 *\defgroup IUAmpToolsMPI
 *
 * The IUAmpToolsMPI directory contains classes for parallel-processing
 * implementations of IUAmpTools using the Message Passing Interface (MPI).
 * These classes allow fits to executed in parallel on multiple CPUs.  The
 * classes enable a parallel computation of both the likelihood and associated
 * normalization integrals.  In general the classes inherit from their 
 * counterparts in IUAmpTools.  The design goal was to only implement 
 * MPI transactions in these classes -- this avoids replication of key
 * functions of the IUAmpTools framework.  Where possible the constructors
 * to these clases were kept the same as IUAmpTools.  This should facilitate
 * an easy switch between single process and MPI implementations.
 */
