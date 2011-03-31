// Copyright (C) 2009-20111 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2009-11-11
// Last changed: 2011-03-31

#include <algorithm>
#include <sstream>
#include <boost/scoped_ptr.hpp>

#include <dolfin/io/File.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/la/LinearAlgebraFactory.h>
#include "TimeSeries.h"

using namespace dolfin;

// Template function for storing objects
template <class T>
void store_object(const T& object, double t,
                  std::vector<double>& times,
                  std::string series_name,
                  std::string type_name)
{
  // Write object
  File file_data(TimeSeries::filename_data(series_name, type_name, times.size()));
  file_data << object;

  // Add time
  times.push_back(t);

  // Store times
  File file_times(TimeSeries::filename_times(series_name, type_name));
  file_times << times;
}

//-----------------------------------------------------------------------------
TimeSeries::TimeSeries(std::string name)
  : _name(name), _cleared(false)
{
  not_working_in_parallel("Storing of data to time series");

  // Set default parameters
  parameters = default_parameters();

  // Read vector times
  std::string filename = TimeSeries::filename_times(_name, "vector");
  if (File::exists(filename))
  {
    // Read from file
    File file(filename);
    file >> _vector_times;
    log(PROGRESS, "Found %d vector sample(s) in time series.", _vector_times.size());
  }
  else
    log(PROGRESS, "No vector samples found in time series.");

  // Read mesh times
  filename = TimeSeries::filename_times(_name, "mesh");
  if (File::exists(filename))
  {
    // Read from file
    File file(filename);
    file >> _mesh_times;
    log(PROGRESS, "Found %d mesh sample(s) in time series.", _mesh_times.size());
  }
  else
    log(PROGRESS, "No mesh samples found in time series.");
}
//-----------------------------------------------------------------------------
TimeSeries::~TimeSeries()
{
  // Do nothing (keep files)
}
//-----------------------------------------------------------------------------
void TimeSeries::store(const GenericVector& vector, double t)
{
  // Clear earlier history first time we store a value
  const bool clear_on_write = this->parameters["clear_on_write"];
  if (!_cleared && clear_on_write)
    clear();

  // Store object
  store_object(vector, t, _vector_times, _name, "vector");
}
//-----------------------------------------------------------------------------
void TimeSeries::store(const Mesh& mesh, double t)
{
  // Clear earlier history first time we store a value
  const bool clear_on_write = this->parameters["clear_on_write"];
  if (!_cleared && clear_on_write)
    clear();

  // Store object
  store_object(mesh, t, _mesh_times, _name, "mesh");
}
//-----------------------------------------------------------------------------
void TimeSeries::retrieve(GenericVector& vector, double t, bool interpolate) const
{
  // Interpolate value
  if (interpolate)
  {
    // Find closest pair
    const std::pair<uint, uint> index_pair = find_closest_pair(t, _vector_times, _name, "vector");
    const uint i0 = index_pair.first;
    const uint i1 = index_pair.second;

    // Special case: same index
    if (i0 == i1)
    {
      File f(filename_data(_name, "vector", i0));
      f >> vector;
      log(PROGRESS, "Reading vector value at t = %g.", _vector_times[0]);
      return;
    }

    log(PROGRESS, "Interpolating vector value at t = %g in interval [%g, %g].",
        t, _vector_times[i0], _vector_times[i1]);

    // Read vectors
    GenericVector& x0(vector);
    boost::scoped_ptr<GenericVector> x1(x0.factory().create_vector());
    File f0(filename_data(_name, "vector", i0));
    File f1(filename_data(_name, "vector", i1));
    f0 >> x0;
    f1 >> *x1;

    // Check that the vectors have the same size
    if (x0.size() != x1->size())
      error("Unable to interpolate vector value; vector sizes don't match (%d and %d).",
            x0.size(), x1->size());

    // Compute weights for linear interpolation
    const double dt = _vector_times[i1] - _vector_times[i0];
    assert(std::abs(dt) > DOLFIN_EPS);
    const double w0 = (_vector_times[i1] - t) / dt;
    const double w1 = 1.0 - w0;

    // Interpolate
    x0 *= w0;
    x0.axpy(w1, *x1);
  }

  // Read closest value
  else
  {
    // Find closest index
    const uint index = find_closest_index(t, _vector_times, _name, "vector");

    log(PROGRESS, "Reading vector at t = %g (close to t = %g).",
        _vector_times[index], t);

    // Read vector
    File file(filename_data(_name, "vector", index));
    file >> vector;
  }
}
//-----------------------------------------------------------------------------
void TimeSeries::retrieve(Mesh& mesh, double t) const
{
  // Get index closest to given time
  const uint index = find_closest_index(t, _mesh_times, _name, "mesh");

  log(PROGRESS, "Reading mesh at t = %g (close to t = %g).",
      _mesh_times[index], t);

  // Read vector
  File file(filename_data(_name, "mesh", index));
  file >> mesh;
}
//-----------------------------------------------------------------------------
Array<double> TimeSeries::vector_times() const
{
  Array<double> times(_vector_times.size());
  for (uint i = 0; i < _vector_times.size(); i++)
    times[i] = _vector_times[i];
  return times;
}
//-----------------------------------------------------------------------------
Array<double> TimeSeries::mesh_times() const
{
  Array<double> times(_mesh_times.size());
  for (uint i = 0; i < _mesh_times.size(); i++)
    times[i] = _mesh_times[i];
  return times;
}
//-----------------------------------------------------------------------------
void TimeSeries::clear()
{
  _vector_times.clear();
  _mesh_times.clear();
  _cleared = true;
}
//-----------------------------------------------------------------------------
std::string TimeSeries::filename_data(std::string series_name,
                                      std::string type_name,
                                      uint index)
{
  std::stringstream s;
  s << series_name << "_" << type_name << "_" << index << ".bin";
  return s.str();
}
//-----------------------------------------------------------------------------
std::string TimeSeries::filename_times(std::string series_name,
                                       std::string type_name)
{
  std::stringstream s;
  s << series_name << "_" << type_name << "_times" << ".bin";
  return s.str();
}
//-----------------------------------------------------------------------------
std::string TimeSeries::str(bool verbose) const
{
  std::stringstream s;

  if (verbose)
  {
    s << str(false) << std::endl << std::endl;

    s << "Vectors:";
    for (uint i = 0; i < _vector_times.size(); ++i)
      s << "  " << i << ": " << _vector_times[i];
    s << std::endl;

    s << "Meshes:";
    for (uint i = 0; i < _mesh_times.size(); ++i)
      s << "  " << i << ": " << _mesh_times[i];
    s << std::endl;
  }
  else
  {
    s << "<Time series with "
      << _vector_times.size()
      << " vector(s) and "
      << _mesh_times.size()
      << " mesh(es)>";
  }

  return s.str();
}
//-----------------------------------------------------------------------------
bool TimeSeries::increasing(const std::vector<double>& times)
{
  for (uint i = 0; i < times.size() - 1; i++)
    if (!(times[i + 1] > times[i]))
      return false;
  return true;
}
//-----------------------------------------------------------------------------
dolfin::uint TimeSeries::find_closest_index(double t,
                                            const std::vector<double>& times,
                                            std::string series_name,
                                            std::string type_name)
{
  // Get closest pair
  const std::pair<uint, uint> index_pair = find_closest_pair(t, times, series_name, type_name);
  const uint i0 = index_pair.first;
  const uint i1 = index_pair.second;

  // Check which is closer
  const uint i = (std::abs(t - times[i0]) < std::abs(t - times[i1]) ? i0 : i1);
  dolfin_debug2("Using closest value t[%d] = %g", i, times[i]);

  return i;
}
//-----------------------------------------------------------------------------
std::pair<dolfin::uint, dolfin::uint>
TimeSeries::find_closest_pair(double t,
                              const std::vector<double>& times,
                              std::string series_name,
                              std::string type_name)
{
  //for (uint i = 0; i < times.size(); i++) cout << " " << times[i]; cout << endl;

  // Must have at least one value stored
  if (times.size() == 0)
    error("Unable to retrieve %s, no %s stored in time series.",
          type_name.c_str(), type_name.c_str());

  // Special case: just one value stored
  if (times.size() == 1)
  {
    dolfin_debug("Series has just one value, returning index 0.");
    return std::make_pair(0, 0);
  }

  // Check whether series is reversed
  const bool reversed = times[0] > times[1];

  // Find lower bound. Note that lower_bound() returns first item
  // larger than t or end of vector if no such item exists.
  std::vector<double>::const_iterator lower;
  if (reversed)
    lower = std::lower_bound(times.begin(), times.end(), t, std::greater<double>());
  else
    lower = std::lower_bound(times.begin(), times.end(), t, std::less<double>());

  // Set indexlower and upper bound
  uint i0 = 0;
  uint i1 = 0;
  if (lower == times.begin())
    i0 = i1 = lower - times.begin();
  else if (lower == times.end())
    i0 = i1 = lower - times.begin() - 1;
  else
  {
    i0 = lower - times.begin() - 1;
    i1 = i0 + 1;
  }

  dolfin_debug1("Looking for value at time t = %g", t);
  dolfin_debug4("Neighboring values are t[%d] = %g and t[%d] = %g",
                i0, times[i0], i1, times[i1]);

  return std::make_pair(i0, i1);
}
//-----------------------------------------------------------------------------
