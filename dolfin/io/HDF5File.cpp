// Copyright (C) 2012 Chris N Richardson
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
// Modified by Garth N. Wells, 2012
//
// First added:  2012-06-01
// Last changed: 2012-09-18

#ifdef HAS_HDF5

#include <cstdio>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>

// Use version 1.6 API for stability
// Fairly easy to switch to later version
// But requires adding extra fields to several calls
//
#define H5_USE_16_API
#include <hdf5.h>

#include <dolfin/common/types.h>
#include <dolfin/common/constants.h>
#include <dolfin/common/MPI.h>
#include <dolfin/common/NoDeleter.h>
#include <dolfin/function/Function.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/log/log.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/LocalMeshData.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshPartitioning.h>
#include <dolfin/mesh/MeshEntityIterator.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/mesh/Vertex.h>

#include "HDF5File.h"

using namespace dolfin;

//-----------------------------------------------------------------------------
HDF5File::HDF5File(const std::string filename) : GenericFile(filename, "H5")
{
  // Do nothing

  // FIXME: Create file here in constructor?
}
//-----------------------------------------------------------------------------
HDF5File::~HDF5File()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
void HDF5File::operator>>(Mesh& input_mesh)
{
  // FIXME: Figure out how to handle multiple meshes in file

  // Get list all datasets in the /Mesh folder
  std::vector<std::string> _dataset_list = dataset_list("/Mesh");

  // TODO: should do a more comprehensive check
  // Check that there are thres data sets (should be starting with
  // "Topology", "Coordinates" and "GlobalIndex"
  if(_dataset_list.size() != 3)
  {
    dolfin_error("HDF5File.cpp",
                 "read mesh from file",
                 "Invalid number of Mesh datasets in HDF5 file");
  }

  // FIXME: Is this function only for distributed meshes?
  LocalMeshData mesh_data;
  mesh_data.clear();

  // Coordinates
  std::string coords_name("/Mesh/");
  // FIXME:  // hopefully 'Coordinates' - need to make more robust
  coords_name.append(_dataset_list[0]);
  std::pair<uint, uint> coords_dim = dataset_dimensions(coords_name);

  // FIXME: This looks weird
  const uint num_global_vertices = coords_dim.first;
  mesh_data.num_global_vertices  = num_global_vertices;

  // FIXME: Document what's going on
  const std::pair<uint, uint> vertex_range = MPI::local_range(num_global_vertices);
  uint num_local_vertices = vertex_range.second-vertex_range.first;
  mesh_data.vertex_indices.reserve(num_local_vertices);
  std::cout << "Reserved space for " << num_local_vertices << " vertices" << std::endl;
  std::vector<double> data;

  // FIXME: This looks seriously broken HDF5::read does not resize data
  data.reserve(num_local_vertices*3); // Mesh always saved in 3D regardless, so may need to decimate
  read(data, vertex_range, coords_name, H5T_NATIVE_DOUBLE, 3);

  std::string global_index_name("/Mesh/");
  // FIXME: With luck...
  global_index_name.append(_dataset_list[1]);
  std::vector<uint> global_index_data;

  // FIXME: This looks seriously broken HDF5::read does not resize data
  global_index_data.reserve(num_local_vertices*2);
  read(global_index_data, vertex_range, global_index_name, H5T_NATIVE_INT, 2);

  printf("Loading %d vertices\n", num_local_vertices);

  for(uint i = 0; i < num_local_vertices; i++)
  {
    std::vector<double> v(&data[i*3], &data[i*3 + coords_dim.second]); // copy correct width (2D or 3D)
    mesh_data.vertex_coordinates.push_back(v);
    mesh_data.vertex_indices.push_back(global_index_data[i*2]);
  }

  // Topology

  // FIXME: get these from somewhere
  mesh_data.gdim = 2;
  mesh_data.tdim = 2;

  std::string topo_name("/Mesh/");
  // FIXME: Make this more robust
  // FIXME: Use better names, i.e. not 'topo'
  topo_name.append(_dataset_list[2]);
  std::pair<uint, uint> topo_dim = dataset_dimensions(topo_name);

  const uint num_global_cells = topo_dim.first;
  mesh_data.num_global_cells = num_global_cells;

  std::pair<uint,uint> cell_range = MPI::local_range(num_global_cells);
  uint num_local_cells = cell_range.second-cell_range.first;
  mesh_data.global_cell_indices.reserve(num_local_cells);
  mesh_data.cell_vertices.reserve(num_local_cells);
  uint num_vertices_per_cell = topo_dim.second;
  mesh_data.num_vertices_per_cell = num_vertices_per_cell;

  std::vector<uint> topo_data(num_local_cells*num_vertices_per_cell);

  // FIXME: This looks seriously broken HDF5::read does not resize data
  topo_data.reserve(num_local_cells*num_vertices_per_cell);
  read(topo_data, cell_range, topo_name, H5T_NATIVE_INT, num_vertices_per_cell);

  // FIXME: The same number of processes *does not* guarantee the same
  // partitioning. At different partitioning might be used, and
  // partitioners often use a random seed.

  // This only works if the partitioning is the same as when it was saved,
  // i.e. the same number of processes
  const uint vertex_offset = MPI::global_offset(num_local_vertices, true);

  // This only works if the partitioning is the same as when it was saved,
  // i.e. the same number of processes
  const uint cell_offset = MPI::global_offset(num_local_cells, true);

  // FIXME: Do not use i, j, k, etc for iterators. Use meaningful name
  // FIXME: Use meaningful names, i.e. not 'ci'
  uint ci = cell_offset;
  for(std::vector<uint>::iterator i = topo_data.begin();
          i != topo_data.end(); i += num_vertices_per_cell)
  {
    std::vector<uint> cell;
    mesh_data.global_cell_indices.push_back(ci);
    ci++;

    for(uint j = 0; j < num_vertices_per_cell; j++)
    {
      uint idx = *(i + j) - vertex_offset;
      cell.push_back(mesh_data.vertex_indices[idx]);
    }
    mesh_data.cell_vertices.push_back(cell);
  }

  std::stringstream s;
  s << "MPI: " << MPI::process_number() << std::endl;
  s << "Cells" << std::endl;

  for(uint i = 0; i < num_local_cells; i++)
  {
    s << "[" << mesh_data.global_cell_indices[i] << "] ";
    for(uint j = 0; j < num_vertices_per_cell; j++)
      s << mesh_data.cell_vertices[i][j] << ",";
    s << std::endl;

  }

  s << "Vertices" << std::endl;
  for(uint i = 0; i < num_local_vertices; i++)
  {
    s << "[" << mesh_data.vertex_indices[i] << "] ";
    for(uint j = 0; j < mesh_data.tdim; j++)
      s << mesh_data.vertex_coordinates[i][j] << ",";
    s << std::endl;

  }
  std::cout << s.str();

  // Build distributed mesh
  MeshPartitioning::build_distributed_mesh(input_mesh, mesh_data);
}
//-----------------------------------------------------------------------------
void HDF5File::operator<<(const Mesh& mesh)
{
  write_mesh(mesh, true);
}
//-----------------------------------------------------------------------------
void HDF5File::write_mesh(const Mesh& mesh, bool true_topology_indices)
{
  // Clear file when writing to file for the first time
  if(counter == 0)
    create();

  // Get local mesh data
  const uint cell_dim = mesh.topology().dim();
  const uint num_local_cells = mesh.num_cells();
  const uint num_local_vertices = mesh.num_vertices();
  const CellType::Type _cell_type = mesh.type().cell_type();
  const std::string cell_type = CellType::type2string(_cell_type);

  // Get cell offset and local cell range
  const uint cell_offset = MPI::global_offset(num_local_cells, true);
  const std::pair<uint, uint> cell_range(cell_offset, cell_offset + num_local_cells);

  // Get vertex offset and local vertex range
  const uint vertex_offset = MPI::global_offset(num_local_vertices, true);
  const std::pair<uint, uint> vertex_range(vertex_offset, vertex_offset + num_local_vertices);

  // FIXME: This is a bit clumsy because of lack of good support in DOLFIN
  //        for local/global indices. Replace when support in DOLFIN is
  //        improved
  // Get global vertex indices
  MeshFunction<uint> v_indices(mesh, 0);
  if (MPI::num_processes() == 1)
  {
    for (VertexIterator v(mesh); !v.end(); ++v)
      v_indices[*v] = v->index();
  }
  else
    v_indices = mesh.parallel_data().global_entity_indices(0);

  // Get vertex indices
  std::vector<uint> vertex_indices;
  std::vector<double> vertex_coords;
  vertex_indices.reserve(2*num_local_vertices);
  vertex_coords.reserve(3*num_local_vertices);
  const uint process_number = MPI::process_number();
  for (VertexIterator v(mesh); !v.end(); ++v)
  {
    // Vertex gloabl index and process number
    vertex_indices.push_back(v_indices[*v]);
    vertex_indices.push_back(process_number);

    // Vertex coordinates
    const Point p = v->point();
    vertex_coords.push_back(p.x());
    vertex_coords.push_back(p.y());
    vertex_coords.push_back(p.z());
  }

  // Write vertex data to HDF5 file
  const std::string coord_dataset = mesh_coords_dataset_name(mesh);
  if (!dataset_exists(coord_dataset))
  {
    write(vertex_indices, vertex_range, mesh_index_dataset_name(mesh), 2);
    write(vertex_coords, vertex_range, coord_dataset, 3);
  }

  // Get cell connectivity
  // NOTE: For visualisation via XDMF, the vertex indices correspond
  //       to the local vertex poistion, and not the true vertex indices.
  std::vector<uint> topological_data;
  if (true_topology_indices)
  {
    // Build connectivity using true vertex indices
    for (CellIterator cell(mesh); !cell.end(); ++cell)
      for (VertexIterator v(*cell); !v.end(); ++v)
        topological_data.push_back(v_indices[*v]);
  }
  else
  {
    // Build connectivity using contiguous vertex indices
    for (CellIterator cell(mesh); !cell.end(); ++cell)
      for (VertexIterator v(*cell); !v.end(); ++v)
        topological_data.push_back(v->index() + vertex_range.first);
  }

  // Write connectivity to HDF5 file
  const std::string topology_dataset = mesh_topo_dataset_name(mesh);
  if (!dataset_exists(topology_dataset))
  {
    write(topological_data, cell_range, topology_dataset, cell_dim + 1);
    add_attribute(topology_dataset, "celltype", cell_type);
  }
}
//-----------------------------------------------------------------------------
void HDF5File::operator<< (const GenericVector& x)
{
  // Get local range;
  std::pair<uint, uint> range = x.local_range(0);

  // Get all local data
  std::vector<double> data;
  x.get_local(data);

  // Overwrite any existing file
  if (counter == 0)
    create();

  // Write to HDF5 file
  std::stringstream s("");
  const std::string name = "/Vector/" + boost::lexical_cast<std::string>(counter);
  write(data, range, name.c_str(), 1);

  // Increment counter
  counter++;
}
//-----------------------------------------------------------------------------
void HDF5File::operator>> (GenericVector& input)
{
  const std::pair<uint, uint> range = input.local_range(0);
  std::vector<double> data(range.second - range.first);
  read(data, range, "/Vector/0", H5T_NATIVE_DOUBLE, 1);
  input.set_local(data);
}
//-----------------------------------------------------------------------------
void HDF5File::create()
{
  // make empty HDF5 file
  // overwriting any existing file
  // create some default 'folders' for storing different datasets

  hid_t  file_id;     // file and dataset identifiers
  hid_t  plist_id;    // property list identifier
  hid_t  group_id;
  herr_t status;

  MPICommunicator comm;
  MPIInfo info;

  // Set parallel access with communicator
  plist_id = H5Pcreate(H5P_FILE_ACCESS);
  status = H5Pset_fapl_mpio(plist_id, *comm, *info);
  dolfin_assert(status != HDF5_FAIL);


  file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  dolfin_assert(file_id != HDF5_FAIL);

  // create subgroups suitable for storing different types of data.
  // DataVector - values for visualisation
  group_id = H5Gcreate(file_id, "/DataVector", H5P_DEFAULT);
  dolfin_assert(group_id != HDF5_FAIL);
  status = H5Gclose (group_id);
  dolfin_assert(status != HDF5_FAIL);

  // Vector - for checkpointing etc
  group_id = H5Gcreate(file_id, "/Vector", H5P_DEFAULT);
  dolfin_assert(group_id != HDF5_FAIL);
  status = H5Gclose (group_id);
  dolfin_assert(status != HDF5_FAIL);

  // Mesh
  group_id = H5Gcreate(file_id, "/Mesh", H5P_DEFAULT);
  dolfin_assert(group_id != HDF5_FAIL);
  status = H5Gclose (group_id);
  dolfin_assert(status != HDF5_FAIL);

  // Release file-access template 
  status = H5Pclose(plist_id);
  dolfin_assert(status != HDF5_FAIL);
  status = H5Fclose(file_id);
  dolfin_assert(status != HDF5_FAIL);
}
//-----------------------------------------------------------------------------
template <typename T>
void HDF5File::read(std::vector<T>& data,  const std::pair<uint, uint> range,
                     const std::string dataset_name,
                     const int h5type, const uint width) const
{
  // read a generic block of 2D data from a HDF5 dataset in parallel

  dolfin_assert(data);
  if(data.size() != range.second - range.first)
    data.resize(range.second - range.first);
  
  hid_t file_id;      // HDF5 file ID
  hid_t plist_id;     // File access template
  hid_t filespace;    // File dataspace ID
  hid_t memspace;     // memory dataspace ID
  hid_t dset_id;      // Dataset ID
  herr_t status;      // Generic return value

  // MPI
  MPICommunicator comm;
  MPIInfo info;

  // Hyperslab selection
  hsize_t offset[2] = {range.first, 0};
  hsize_t count[2] = {range.second - range.first, width};
  
  // Setup file access template
  plist_id = H5Pcreate(H5P_FILE_ACCESS);
  dolfin_assert(plist_id != HDF5_FAIL);

  // Set parallel access with communicator
  status = H5Pset_fapl_mpio(plist_id, *comm, *info);
  dolfin_assert(status != HDF5_FAIL);

  // Open the file collectively
  file_id = H5Fopen(filename.c_str(), H5F_ACC_RDWR,plist_id);
  dolfin_assert(file_id != HDF5_FAIL);

  // Release file-access template 
  status = H5Pclose(plist_id);
  dolfin_assert(status != HDF5_FAIL);

  // open the dataset collectively
  dset_id = H5Dopen(file_id, dataset_name.c_str());
  dolfin_assert(dset_id != HDF5_FAIL);

  // Create a file dataspace independently
  filespace = H5Dget_space (dset_id);
  dolfin_assert(filespace != HDF5_FAIL);

  status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL,
                               count, NULL);
  dolfin_assert(status != HDF5_FAIL);

  // Create a memory dataspace independently
  memspace = H5Screate_simple (2, count, NULL);
  dolfin_assert (memspace != HDF5_FAIL);

  // read data independently
  status = H5Dread(dset_id, h5type, memspace, filespace,
                   H5P_DEFAULT, data.data());
  dolfin_assert(status != HDF5_FAIL);

  // close dataset collectively
  status = H5Dclose(dset_id);
  dolfin_assert(status != HDF5_FAIL);

  // release all IDs created
  status = H5Sclose(filespace);
  dolfin_assert(status != HDF5_FAIL);

  // close the file collectively
  status = H5Fclose(file_id);
  dolfin_assert(status != HDF5_FAIL);
}
//-----------------------------------------------------------------------------
void HDF5File::write(const std::vector<double>& data,
                     const std::pair<uint, uint> range,
                     const std::string dataset_name, const uint width)
{
  // Write data to existing HDF file as defined by range blocks on each process
  // range: the local range on this processor
  // width: is the width of the dataitem (e.g. 3 for x, y, z data)
  write(data, range, dataset_name, H5T_NATIVE_DOUBLE, width);
}
//-----------------------------------------------------------------------------
void HDF5File::write(const std::vector<uint>& data,
                     const std::pair<uint, uint> range,
                     const std::string dataset_name, const uint width)
{
  write(data, range, dataset_name, H5T_NATIVE_INT, width);
}
//-----------------------------------------------------------------------------
template <typename T>
void HDF5File::write(const std::vector<T>& data,
                     const std::pair<uint, uint> range,
                     const std::string dataset_name,
                     const int h5type, const uint width) const
{
  // File and dataset identifiers
  hid_t file_id, dset_id;

  // File and memory dataspace identifiers
  hid_t filespace, memspace;

  // Hyperslab selection parameters
  hsize_t count[2]  = {range.second - range.first, width};
  hsize_t offset[2] = {range.first, 0};

  // Dataset dimensions
  hsize_t dimsf[2] = {MPI::sum(count[0]), width} ;

  // Property list identifier
  hid_t  plist_id;
  herr_t status;

  // Get MPI objects
  MPICommunicator comm;
  MPIInfo info;

  // Set parallel access with communicator
  plist_id = H5Pcreate(H5P_FILE_ACCESS);
  status = H5Pset_fapl_mpio(plist_id, *comm, *info);
  dolfin_assert(status != HDF5_FAIL);

  // Try to open existing HDF5 file
  file_id = H5Fopen(filename.c_str(), H5F_ACC_RDWR, plist_id);
  dolfin_assert(file_id != HDF5_FAIL);

  // Release file-access template 
  status = H5Pclose(plist_id);
  dolfin_assert(status != HDF5_FAIL);

  // Create a global 2D data space
  filespace = H5Screate_simple(2, dimsf, NULL);
  dolfin_assert(filespace != HDF5_FAIL);

  // Create global dataset (using dataset_name)
  dset_id = H5Dcreate(file_id, dataset_name.c_str(), h5type, filespace,
                      H5P_DEFAULT);
  dolfin_assert(dset_id != HDF5_FAIL);

  // Close global data space
  status = H5Sclose(filespace);
  dolfin_assert(status != HDF5_FAIL);

  // Create a local 2D data space
  memspace = H5Screate_simple(2, count, NULL);

  // Create a file dataspace within the global space - a hyperslab
  filespace = H5Dget_space(dset_id);
  status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
  dolfin_assert(status != HDF5_FAIL);

  // Set parallel access with communicator
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  status = H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  dolfin_assert(status != HDF5_FAIL);

  // Write local dataset into selected hyperslab
  status = H5Dwrite(dset_id, h5type, memspace, filespace, plist_id, data.data());
  dolfin_assert(status != HDF5_FAIL);

  // close dataset collectively
  status = H5Dclose(dset_id);
  dolfin_assert(status != HDF5_FAIL);

  // close hyperslab
  status = H5Sclose(filespace);
  dolfin_assert(status != HDF5_FAIL);

  // close local dataset
  status = H5Sclose(memspace);
  dolfin_assert(status != HDF5_FAIL);

  // Release file-access template 
  status = H5Pclose(plist_id);
  dolfin_assert(status != HDF5_FAIL);

  // Close file
  status = H5Fclose(file_id);
  dolfin_assert(status != HDF5_FAIL);
}
//-----------------------------------------------------------------------------
bool HDF5File::dataset_exists(const std::string dataset_name) const
{
  MPICommunicator comm;
  MPIInfo info;
  herr_t status;

  // Set parallel access with communicator
  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  status = H5Pset_fapl_mpio(plist_id,*comm, *info);
  dolfin_assert(status != HDF5_FAIL);

  // Try to open existing HDF5 file
  hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDWR, plist_id);
  dolfin_assert(file_id != HDF5_FAIL);

  // Release file-access template 
  status = H5Pclose(plist_id);
  dolfin_assert(status != HDF5_FAIL);

  // disable error reporting
  herr_t (*old_func)(void*);
  void *old_client_data;
  H5Eget_auto(&old_func, &old_client_data);

  // redirect error reporting (to none)
  status = H5Eset_auto(NULL, NULL);
  dolfin_assert(status != HDF5_FAIL);

  // try to open dataset - returns HDF5_FAIL if non-existent
  hid_t dset_id = H5Dopen(file_id, dataset_name.c_str());

  if(dset_id != HDF5_FAIL)
    H5Dclose(dset_id);

  // re-enable error reporting
  status = H5Eset_auto(old_func, old_client_data);
  dolfin_assert(status != HDF5_FAIL);

  status = H5Fclose(file_id);
  dolfin_assert(status != HDF5_FAIL);

  return (dset_id != HDF5_FAIL);
}
//-----------------------------------------------------------------------------
std::vector<std::string> HDF5File::dataset_list(const std::string group_name) const
{
  char namebuf[HDF5_MAXSTRLEN];

  MPICommunicator comm;
  MPIInfo info;
  herr_t status;

  // Set parallel access with communicator
  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  status = H5Pset_fapl_mpio(plist_id,*comm, *info);
  dolfin_assert(status != HDF5_FAIL);

  // Try to open existing HDF5 file
  hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDWR, plist_id);
  dolfin_assert(file_id != HDF5_FAIL);

  // FIXME: document
  // Release file-access template 
  status = H5Pclose(plist_id);
  dolfin_assert(status != HDF5_FAIL);

  // FIXME: document
  hid_t group_id = H5Gopen(file_id,group_name.c_str());
  dolfin_assert(group_id != HDF5_FAIL);

  // count how many datasets in the group
  hsize_t num_obj;
  status = H5Gget_num_objs(group_id, &num_obj);
  dolfin_assert(status != HDF5_FAIL);

  // FIXME: document
  std::vector<std::string> lvec;
  std::string str;
  // go through all objects
  for(hsize_t i=0; i<num_obj; i++)
  {
    H5Gget_objname_by_idx(group_id, i, namebuf, HDF5_MAXSTRLEN);
    str=namebuf;
    lvec.push_back(str);
  }

  // FIXME: document
  status = H5Gclose(group_id);
  dolfin_assert(status != HDF5_FAIL);

  // FIXME: document
  status = H5Fclose(file_id);
  dolfin_assert(status != HDF5_FAIL);

  return lvec;
}
//-----------------------------------------------------------------------------
std::pair<uint, uint> HDF5File::dataset_dimensions(const std::string dataset_name) const
{
  // Get dimensions of a 2D dataset

  hsize_t cur_size[2];   // current dataset dimensions
  hsize_t max_size[2];   // maximum dataset dimensions
  hid_t   space;         // data space
  int     ndims;         // dimensionality

  MPICommunicator comm;
  MPIInfo info;
  herr_t status;

  // Set parallel access with communicator
  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  status = H5Pset_fapl_mpio(plist_id,*comm, *info);
  dolfin_assert(status != HDF5_FAIL);

  // Try to open existing HDF5 file
  hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDWR, plist_id);
  dolfin_assert(file_id != HDF5_FAIL);

  // FIXME: document
  // Release file-access template 
  status = H5Pclose(plist_id);
  dolfin_assert(status != HDF5_FAIL);

  // FIXME: document
  hid_t dset_id = H5Dopen(file_id, dataset_name.c_str());
  dolfin_assert(dset_id != HDF5_FAIL);

  // Create a file dataspace independently
  space = H5Dget_space(dset_id);
  ndims = H5Sget_simple_extent_dims(space, cur_size, max_size);
  dolfin_assert(ndims == 2);

  // close dataset collectively
  status = H5Dclose(dset_id);
  dolfin_assert(status != HDF5_FAIL);

  // FIXME: document
  status = H5Fclose(file_id);
  dolfin_assert(status != HDF5_FAIL);

  return std::pair<uint,uint>(cur_size[0],cur_size[1]);
}
//-----------------------------------------------------------------------------
void HDF5File::add_attribute(const std::string dataset_name,
                             const std::string attribute_name,
                             const std::string attribute_value)
{
  MPICommunicator comm;
  MPIInfo info;
  herr_t status;

  // Set parallel access with communicator
  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  status = H5Pset_fapl_mpio(plist_id,*comm, *info);
  dolfin_assert(status != HDF5_FAIL);

  // Try to open existing HDF5 file
  hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDWR, plist_id);
  dolfin_assert(file_id != HDF5_FAIL);

  // FIXME: document
  // Release file-access template 
  status = H5Pclose(plist_id);
  dolfin_assert(status != HDF5_FAIL);

  // FIXME: document
  hid_t dset_id = H5Dopen(file_id, dataset_name.c_str());
  dolfin_assert(dset_id != HDF5_FAIL);

  // add string attribute
  // FIXME: document better
  hid_t datatype_id = H5Tcopy (H5T_C_S1);
  status = H5Tset_size (datatype_id, attribute_value.size());
  hid_t dataspaces_id = H5Screate (H5S_SCALAR);
  hid_t attribute_id = H5Acreate (dset_id, attribute_name.c_str(), datatype_id,
                                  dataspaces_id, H5P_DEFAULT);
  status = H5Awrite(attribute_id, datatype_id, attribute_value.c_str());
  dolfin_assert(status != HDF5_FAIL);

  // FIXME: document
  status = H5Aclose(attribute_id);
  dolfin_assert(status != HDF5_FAIL);

  // close dataset collectively
  status = H5Dclose(dset_id);
  dolfin_assert(status != HDF5_FAIL);

  // FIXME: document
  status = H5Fclose(file_id);
  dolfin_assert(status != HDF5_FAIL);
}
//-----------------------------------------------------------------------------
std::string HDF5File::get_attribute(const std::string dataset_name,
                                    const std::string attribute_name) const
{
  MPICommunicator comm;
  MPIInfo info;
  herr_t status;

  // Set parallel access with communicator
  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  status = H5Pset_fapl_mpio(plist_id,*comm, *info);
  dolfin_assert(status != HDF5_FAIL);

  // Try to open existing HDF5 file
  hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDWR, plist_id);
  dolfin_assert(file_id != HDF5_FAIL);

  // FIXME: document
  // Release file-access template 
  status = H5Pclose(plist_id);
  dolfin_assert(status != HDF5_FAIL);

  // FIXME: document
  hid_t dset_id = H5Dopen(file_id, dataset_name.c_str());
  dolfin_assert(dset_id != HDF5_FAIL);

  // FIXME: document
  hid_t attr_id = H5Aopen(dset_id, attribute_name.c_str(), H5P_DEFAULT);
  hid_t filetype = H5Aget_type(attr_id);
  int slen = H5Tget_size(filetype);
  slen++;

  //  hid_t space_id = H5Aget_space(attr_id);
  hid_t memtype = H5Tcopy (H5T_C_S1);

  // FIXME: document
  status = H5Tset_size(memtype,slen);
  dolfin_assert(status != HDF5_FAIL);

  // FIXME: document
  std::vector<char> str(slen);
  status = H5Aread(attr_id, memtype, &str[0]);

  // FIXME: document
  status = H5Aclose(attr_id);
  dolfin_assert(status != HDF5_FAIL);

  // close dataset collectively
  status = H5Dclose(dset_id);
  dolfin_assert(status != HDF5_FAIL);

  // FIXME: document
  status = H5Fclose(file_id);
  dolfin_assert(status != HDF5_FAIL);

  return std::string(&str[0]);
}
//-----------------------------------------------------------------------------
std::string HDF5File::mesh_coords_dataset_name(const Mesh& mesh) const
{
  std::stringstream dataset_name;
  dataset_name << "/Mesh/Coordinates_" << std::setfill('0')
          << std::hex << std::setw(8) << mesh.coordinates_hash();
  return dataset_name.str();
}
//-----------------------------------------------------------------------------
std::string HDF5File::mesh_index_dataset_name(const Mesh& mesh) const
{
  std::stringstream dataset_name;
  dataset_name << "/Mesh/GlobalIndex_" << std::setfill('0')
          << std::hex << std::setw(8) << mesh.coordinates_hash();
  return dataset_name.str();
}
//-----------------------------------------------------------------------------
std::string HDF5File::mesh_topo_dataset_name(const Mesh& mesh) const
{
  std::stringstream dataset_name;
  dataset_name << "/Mesh/Topology_" << std::setfill('0')
          << std::hex << std::setw(8) << mesh.topology_hash();
  return dataset_name.str();
}
//-----------------------------------------------------------------------------

#endif
