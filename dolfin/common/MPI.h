// Copyright (C) 2007-2013 Magnus Vikstrøm and Garth N. Wells
//
// This file is part of DOLFIN.
//
// DOLFIN is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// DOLFIN is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with DOLFIN. If not, see <http://www.gnu.org/licenses/>.
//
// Modified by Ola Skavhaug 2008-2009
// Modified by Anders Logg 2008-2011
// Modified by Niclas Jansson 2009
// Modified by Joachim B Haga 2012
//
// First added:  2007-11-30
// Last changed: 2013-01-11

#ifndef __MPI_DOLFIN_WRAPPER_H
#define __MPI_DOLFIN_WRAPPER_H

#include <utility>
#include <vector>

#ifdef HAS_MPI
#define MPICH_IGNORE_CXX_SEEK 1
#include <boost/serialization/utility.hpp>
#include <boost/mpi.hpp>
#include <mpi.h>
#endif

#include <dolfin/log/dolfin_log.h>

namespace dolfin
{

  #ifdef HAS_MPI

  class MPICommunicator
  {

  public:

    /// Create communicator (copy of MPI_COMM_WORLD)
    MPICommunicator();

    /// Destructor
    ~MPICommunicator();

    /// Dereference operator
    MPI_Comm& operator*();

  private:

    MPI_Comm communicator;

  };

  class MPIInfo
  {
  public:
    MPIInfo();
    ~MPIInfo();
    MPI_Info& operator*();
  private:

    MPI_Info info;

  };

  #endif

  class MPINonblocking
  {
    /// This class provides stateful (single communicator) non-blocking
    /// MPI functionality.

  public:

    /// Destroy instance (waits for outstanding requests)
    ~MPINonblocking()
    {
      #ifdef HAS_MPI
      wait_all();
      #endif
    }

    /// Non-blocking send and receive
    template<typename T>
    void send_recv(const T& send_value, unsigned int dest, T& recv_value,
                   unsigned int source);

    /// Non-blocking send and receive with tag
    template<typename T>
    void send_recv(const T& send_value, unsigned int dest_tag, unsigned int dest,
                   T& recv_value, unsigned int source_tag, unsigned int source);

    /// Wait for all requests to finish
    void wait_all();

  private:

    #ifdef HAS_MPI
    MPICommunicator mpi_comm;
    std::vector<boost::mpi::request> reqs;
    #endif

  };

  /// This class provides utility functions for easy communcation
  /// with MPI.

  class MPI
  {
  public:

    /// Return proccess number
    static unsigned int process_number();

    /// Return number of processes
    static unsigned int num_processes();

    /// Determine whether we should broadcast (based on current parallel
    /// policy)
    static bool is_broadcaster();

    /// Determine whether we should receive (based on current parallel
    /// policy)
    static bool is_receiver();

    /// Set a barrier (synchronization point)
    static void barrier();

    /// Send in_values[p0] to process p0 and receive values from process
    /// p1 in out_values[p1]
    template<typename T>
    static void all_to_all(const std::vector<std::vector<T> >& in_values,
                           std::vector<std::vector<T> >& out_values)
    {
      #ifdef HAS_MPI
      MPICommunicator mpi_comm;
      boost::mpi::communicator comm(*mpi_comm, boost::mpi::comm_attach);
      boost::mpi::all_to_all(comm, in_values, out_values);
      #else
      dolfin_assert(in_values.size() == 1);
      out_values = in_values;
      #endif

    }

    /// Distribute local arrays on a group of processes (typically
    /// neighbours from GenericDofMap::neighbours()). It is important
    /// that each process' group includes exactly the processes that
    /// has it in their groups, otherwise it will deadlock.
    template<typename T, typename S>
    static void distribute(const std::set<S> group,
                           const std::map<S, T>& in_values_per_dest,
                           std::map<S, T>& out_values_per_src);

    /// Broadcast value from broadcaster process to all processes
    template<typename T>
    static void broadcast(T& value, unsigned int broadcaster=0)
    {
      #ifdef HAS_MPI
      MPICommunicator mpi_comm;
      boost::mpi::communicator comm(*mpi_comm, boost::mpi::comm_attach);
      boost::mpi::broadcast(comm, value, broadcaster);
      #endif
    }

    /// Scatter in_values[i] to process i
    template<typename T>
    static void scatter(const std::vector<T>& in_values,
                        T& out_value, unsigned int sending_process=0)
    {
      #ifdef HAS_MPI
      MPICommunicator mpi_comm;
      boost::mpi::communicator comm(*mpi_comm, boost::mpi::comm_attach);
      boost::mpi::scatter(comm, in_values, out_value, sending_process);
      #else
      dolfin_assert(sending_process == 0);
      dolfin_assert(in_values.size() == 1);
      out_value = in_values[0];
      #endif
    }

    /// Gather values on one process (wrapper for boost::mpi::gather)
    template<typename T>
    static void gather(const T& in_value, std::vector<T>& out_values,
                       unsigned int receiving_process=0)
    {
      #ifdef HAS_MPI
      MPICommunicator mpi_comm;
      boost::mpi::communicator comm(*mpi_comm, boost::mpi::comm_attach);
      boost::mpi::gather(comm, in_value, out_values, receiving_process);
      #else
      out_values.clear();
      out_values.push_back(in_value);
      #endif
    }

    /// Gather values, one from each process (wrapper for boost::mpi::all_gather)
    template<typename T>
    static void all_gather(const T& in_value, std::vector<T>& out_values)
    {
      #ifdef HAS_MPI
      MPICommunicator mpi_comm;
      boost::mpi::communicator comm(*mpi_comm, boost::mpi::comm_attach);
      boost::mpi::all_gather(comm, in_value, out_values);
      #else
      out_values.clear();
      out_values.push_back(in_value);
      #endif
    }

    /// Return global max value
    template<typename T> static T max(const T& value)
    {
     #ifdef HAS_MPI
     return all_reduce(value, boost::mpi::maximum<T>());
     #else
     return value;
     #endif
    }

    /// Return global min value
    template<typename T> static T min(const T& value)
    {
      #ifdef HAS_MPI
      return all_reduce(value, boost::mpi::minimum<T>());
      #else
      return value;
      #endif
    }

    /// Sum values and return sum
    template<typename T> static T sum(const T& value)
    {
      #ifdef HAS_MPI
      return all_reduce(value, std::plus<T>());
      #else
      return value;
      #endif
    }

    /// All reduce
    template<typename T, typename X> static T all_reduce(const T& value, X op)
    {
      #ifdef HAS_MPI
      MPICommunicator mpi_comm;
      boost::mpi::communicator comm(*mpi_comm, boost::mpi::comm_attach);
      T out;
      boost::mpi::all_reduce(comm, value, out, op);
      return out;
      #else
      dolfin_error("MPI.h",
                   "call MPI::all_reduce",
                   "DOLFIN has been configured without MPI support");
      return T(0);
      #endif
    }

    /// Find global offset (index) (wrapper for MPI_(Ex)Scan with MPI_SUM as
    /// reduction op)
    static std::size_t global_offset(std::size_t range, bool exclusive);

    /// Send-receive data. Note that if the number of posted send-receives
    /// may differ between processes, another interface (such as
    /// MPINonblocking::send_recv) must be used since duplicating the
    /// communicator requires participation from all processes.
    template<typename T>
    static void send_recv(const T& send_value, unsigned int dest,
                          T& recv_value, unsigned int source)
    {
      #ifdef HAS_MPI
      MPINonblocking mpi;
      mpi.send_recv(send_value, dest, recv_value, source);
      #else
      dolfin_error("MPI.h",
                   "call MPI::send_recv",
                   "DOLFIN has been configured without MPI support");
      #endif
    }

    /// Return local range for local process, splitting [0, N - 1] into
    /// num_processes() portions of almost equal size
    static std::pair<std::size_t, std::size_t> local_range(std::size_t N);

    /// Return local range for given process, splitting [0, N - 1] into
    /// num_processes() portions of almost equal size
    static std::pair<std::size_t, std::size_t>
        local_range(unsigned int process, std::size_t N);

    /// Return local range for given process, splitting [0, N - 1] into
    /// num_processes portions of almost equal size
    static std::pair<std::size_t, std::size_t>
        local_range(unsigned int process, std::size_t N,
                    unsigned int num_processes);

    /// Return which process owns index (inverse of local_range)
    static unsigned int index_owner(std::size_t index, std::size_t N);

  private:

    #ifndef HAS_MPI
    static void error_no_mpi(const char *where)
    {
      dolfin_error("MPI.h", where, "DOLFIN has been configured without MPI support");
    }
    #endif

    #ifdef HAS_MPI
    // Return MPI data type
    template<typename T> static MPI_Datatype mpi_type()
    {
      dolfin_error("MPI.h",
                   "perform MPI operation",
                   "MPI data type unknown");
      return MPI_CHAR;
    }
    #endif

  };

  #ifdef HAS_MPI
  // Specialisations for MPI_Datatypes
  template<> inline MPI_Datatype MPI::mpi_type<double>() { return MPI_DOUBLE; }
  template<> inline MPI_Datatype MPI::mpi_type<int>() { return MPI_INT; }
  template<> inline MPI_Datatype MPI::mpi_type<long int>() { return MPI_LONG; }
  template<> inline MPI_Datatype MPI::mpi_type<unsigned int>() { return MPI_UNSIGNED; }
  template<> inline MPI_Datatype MPI::mpi_type<unsigned long>() { return MPI_UNSIGNED_LONG; }
  #endif

  //---------------------------------------------------------------------------
  template<typename T, typename S>
  void dolfin::MPI::distribute(const std::set<S> processes_group,
                               const std::map<S, T>& in_values_per_dest,
                               std::map<S, T>& out_values_per_src)
  {
    #ifdef HAS_MPI
    typedef typename std::map<S, T>::const_iterator map_const_iterator;
    typedef typename std::map<S, T>::iterator map_iterator;
    dolfin::MPINonblocking mpi;
    const T no_data;

    // Send and receive values to all processes in groups
    // (non-blocking). If a given process is not found in
    // in_values_per_dest, send empty data.
    out_values_per_src.clear();
    typename std::set<S>::const_iterator dest;
    for (dest = processes_group.begin(); dest != processes_group.end(); ++dest)
    {
      map_const_iterator values = in_values_per_dest.find(*dest);
      if (values != in_values_per_dest.end())
        mpi.send_recv(values->second, *dest, out_values_per_src[*dest], *dest);
      else
        mpi.send_recv(no_data, *dest, out_values_per_src[*dest], *dest);
    }

    // Wait for all MPI calls before modifying out_values_per_src
    mpi.wait_all();

    // Remove received no_data entries.
    map_iterator it = out_values_per_src.begin();
    while (it != out_values_per_src.end())
    {
      map_iterator tmp = it++;
      if (tmp->second.empty())
        out_values_per_src.erase(tmp);
    }
    #else
    error_no_mpi("call MPI::distribute");
    #endif
  }
  //-----------------------------------------------------------------------------
  template<typename T>
  void dolfin::MPINonblocking::send_recv(const T& send_value, unsigned int dest,
                                         T& recv_value, unsigned int source)
  {
    MPINonblocking::send_recv(send_value, 0, dest, recv_value, 0, source);
  }
  //-----------------------------------------------------------------------------
  template<typename T>
  void dolfin::MPINonblocking::send_recv(const T& send_value, unsigned int dest_tag,
                                         unsigned int dest,
                                         T& recv_value, unsigned int source_tag,
                                         unsigned int source)
  {
    #ifdef HAS_MPI
    boost::mpi::communicator comm(*mpi_comm, boost::mpi::comm_attach);
    reqs.push_back(comm.isend(dest, dest_tag, send_value));
    reqs.push_back(comm.irecv(source, source_tag, recv_value));
    #else
    dolfin_error("MPI.h",
                  "call MPINonblocking::send_recv",
                  "DOLFIN has been configured without MPI support");
    #endif
  }
  //-----------------------------------------------------------------------------
}

#endif
