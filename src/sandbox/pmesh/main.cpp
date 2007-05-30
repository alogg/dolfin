// Copyright (C) 2007 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2007-10-30
// Last changed:
//
// This file is used for testing distribution of the mesh using MPI

#include <dolfin.h>
#include <mpi.h>

using namespace dolfin;
//-----------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  cout << "Starting mesh distribution" << endl;    

  // Initialise MPI  
  MPI_Init(&argc, &argv);

  // Get number of processes
  int num_processes_int;
  MPI_Comm_size(MPI_COMM_WORLD, &num_processes_int);

  if( num_processes_int != 2 )
    error("This test program must be run with 2 processes.");

  // Get this process number
  int process_int;
  MPI_Comm_rank(MPI_COMM_WORLD, &process_int);
  unsigned int this_process = process_int;

  //-- MPI test code ----------------------------------------------------------
  int sent_message = -1;
  int received_message = 0;
  int receive_from = -1;
  int send_to      = -1;
  if(this_process == 0)
  {
    receive_from = 1;
    send_to      = 1;
    sent_message = 20;
  }
  if(this_process == 1)
  {
    receive_from = 0;
    send_to      = 0;
    sent_message = 22;
  }
  int tag = 55;
//  int send_tag = 55;

  MPI_Status status;
  MPI_Send(&sent_message,     1, MPI_INT, send_to,     tag, MPI_COMM_WORLD);
  MPI_Recv(&received_message, 1, MPI_INT, receive_from, 55, MPI_COMM_WORLD, &status);

  cout << "Finished send/rec " << this_process << "  " << received_message << endl;
  MPI_Finalize();

  return 0;
}
