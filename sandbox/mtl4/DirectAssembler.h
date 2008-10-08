#ifndef __DOLFIN_D_ASSEMBLER_H
#define __DOLFIN_D_ASSEMBLER_H

#include <ufc.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/common/Array.h>
#include <dolfin/common/Timer.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/MeshData.h>
#include <dolfin/mesh/BoundaryMesh.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/mesh/SubDomain.h>
#include <dolfin/function/Function.h>
#include "dolfin/fem/Form.h"
#include "dolfin/fem/UFC.h"
#include "dolfin/fem/DofMapSet.h"

#include <dolfin/common/timing.h>

namespace dolfin{

  class DirectAssembler
  {
  public:


    /// Constructor
    DirectAssembler(Mesh& mesh) : mesh(mesh)
    {

    }

    /// Destructor
    ~DirectAssembler()
    {

    }

    /// Assemble tensor from given variational form
    template<typename Tensor_t>
    void assemble(Tensor_t& A, Form& form, bool reset_tensor = true)
    {
      form.updateDofMaps(mesh);
      assemble(A, form.form(), form.coefficients(), 
	       form.dofMaps(), 0, 0, 0, reset_tensor);

    }

    /// Assemble tensor from given variational form over a sub domain
    // template<typename Tensor_t>
    // void assemble(Tensor_t& A, Form& form,
    //               const SubDomain& sub_domain, bool reset_tensor = true)
    // {
    //   // Extract cell domains
    //   MeshFunction<uint>* cell_domains = 0;
    //   if (form.form().num_cell_integrals() > 0)
    // 	{
    // 	  cell_domains = new MeshFunction<uint>(mesh, mesh.topology().dim());
    // 	  (*cell_domains) = 1;
    // 	  sub_domain.mark(*cell_domains, 0);
    // 	}

    //   // Extract facet domains
    //   MeshFunction<uint>* facet_domains = 0;
    //   if (form.form().num_exterior_facet_integrals() > 0 ||
    // 	  form.form().num_interior_facet_integrals() > 0)
    // 	{
    // 	  facet_domains = new MeshFunction<uint>(mesh, mesh.topology().dim()-1);
    // 	  (*facet_domains) = 1;
    // 	  sub_domain.mark(*facet_domains, 0);
    // 	}

    //   // Assemble
    //   form.updateDofMaps(mesh);
    //   assemble(A, form.form(), form.coefficients(), form.dofMaps(),
    // 	       cell_domains, facet_domains, facet_domains, reset_tensor);

    //   // Delete domains
    //   if (cell_domains)
    // 	delete cell_domains;
    //   if (facet_domains)
    // 	delete facet_domains;

    // }

    /// Assemble tensor from given variational form over sub domains
    void assemble(GenericTensor& A, Form& form,
                  const MeshFunction<uint>& cell_domains,
                  const MeshFunction<uint>& exterior_facet_domains,
                  const MeshFunction<uint>& interior_facet_domains, 
		  bool reset_tensor = true);
    
    /// Assemble scalar from given variational form
    double assemble(Form& form);
    
    /// Assemble scalar from given variational form over a sub domain
    double assemble(Form& form, const SubDomain& sub_domain);
    
    /// Assemble scalar from given variational form over sub domains
    double assemble(Form& form,
                  const MeshFunction<uint>& cell_domains,
                  const MeshFunction<uint>& exterior_facet_domains,
                  const MeshFunction<uint>& interior_facet_domains);
    
    /// Assemble tensor from given (UFC) form, coefficients and sub domains.
    /// This is the main assembly function in DOLFIN. All other assembly 
    /// functions end up calling this function.
    ///
    /// The MeshFunction arguments can be used to specify assembly over 
    /// subdomains of the mesh cells, exterior facets and interior facets. 
    /// Either a null pointer or an empty MeshFunction may be used to specify 
    /// that the tensor should  be assembled over the entire set of cells 
    /// or facets.
    template<typename Tensor_t>
    void assemble(Tensor_t& A, const ufc::form& form,
                  const Array<Function*>& coefficients,
                  const DofMapSet& dof_map_set,
                  const MeshFunction<uint>* cell_domains,
                  const MeshFunction<uint>* exterior_facet_domains,
                  const MeshFunction<uint>* interior_facet_domains, 
		  bool reset_tensor = true)
    {
      // Note the importance of treating empty mesh functions as null pointers
      // for the PyDOLFIN interface.

      if(reset_tensor)
	error("Experimental assembler: no initialization of tensor");
  
      // Check arguments
      // check(form, coefficients, mesh);

      // Create data structure for local assembly data
      UFC ufc(form, mesh, dof_map_set);

      // Initialize global tensor
      //initGlobalTensor(A, dof_map_set, ufc, reset_tensor);
      A.zero();

      // Assemble over cells
      assembleCells(A, coefficients, dof_map_set, ufc, cell_domains);

      // Assemble over exterior facets 
      assembleExteriorFacets(A, coefficients, dof_map_set, ufc, 
       			     exterior_facet_domains);

      // Assemble over interior facets
      assembleInteriorFacets(A, coefficients, dof_map_set, ufc, 
       			     interior_facet_domains);

      // Finalise assembly of global tensor
      A.apply();
    }

  private:
     
    // Assemble over cells
    template<typename Tensor_t>
    void assembleCells(Tensor_t& A,
                       const Array<Function*>& coefficients,
                       const DofMapSet& dof_map_set,
                       UFC& ufc,
                       const MeshFunction<uint>* domains) const
    {
      Timer timer("Assembly over cells");

      // Skip assembly if there are no cell integrals
      if (ufc.form.num_cell_integrals() == 0)
	return;

      // Cell integral
      ufc::cell_integral* integral = ufc.cell_integrals[0];

      // Assemble over cells
      Progress p("Assembling over cells", mesh.numCells());
      
      //const uint *cells = mesh.cells();
      //uint i;
   
      //#pragma omp parallel for private(i)
      //for(i = 0; i < mesh.numCells(); i++)
      for (CellIterator cell(mesh); !cell.end(); ++cell)
	{
	  //Cell cell(mesh, cells[i]); 
	  
	  // Get integral for sub domain (if any)
	  if (domains && domains->size() > 0)
	    {
	      const uint domain = (*domains)(*cell);
	      if (domain < ufc.form.num_cell_integrals())
		integral = ufc.cell_integrals[domain];
	      else
		continue;
	    }

	  // Update to current cell
	  ufc.update(*cell);

	  // Interpolate coefficients on cell
	  for (uint i = 0; i < coefficients.size(); i++)
	    coefficients[i]->interpolate(ufc.w[i], ufc.cell, 
					 *ufc.coefficient_elements[i], *cell);
    
	  // Tabulate dofs for each dimension
	  for (uint i = 0; i < ufc.form.rank(); i++)
	    dof_map_set[i].tabulate_dofs(ufc.dofs[i], ufc.cell, cell->index());

	  // Tabulate cell tensor
	  integral->tabulate_tensor(ufc.A, ufc.w, ufc.cell);

	  // Add entries to global tensor
	  A.add(ufc.A, ufc.local_dimensions, ufc.dofs);
    
	  p++;
	}
    }

    // Assemble over exterior facets
    template<typename Tensor_t>
    void assembleExteriorFacets(Tensor_t& A,
                                const Array<Function*>& coefficients,
                                const DofMapSet& dof_map_set,
                                UFC& ufc,
                                const MeshFunction<uint>* domains) const
    {
      // Skip assembly if there are no exterior facet integrals
      if (ufc.form.num_exterior_facet_integrals() == 0)
	return;
  
      // Exterior facet integral
      ufc::exterior_facet_integral* integral = ufc.exterior_facet_integrals[0];

      // Create boundary mesh
      BoundaryMesh boundary(mesh);
      MeshFunction<uint>* cell_map = boundary.data().meshFunction("cell map");
      dolfin_assert(cell_map);

      // Assemble over exterior facets (the cells of the boundary)
      Progress p("Assembling over exterior facets", boundary.numCells());
      for (CellIterator boundary_cell(boundary); !boundary_cell.end(); 
	   ++boundary_cell)
	{
	  // Get mesh facet corresponding to boundary cell
	  Facet mesh_facet(mesh, (*cell_map)(*boundary_cell));

	  // Get integral for sub domain (if any)
	  if (domains && domains->size() > 0)
	    {
	      const uint domain = (*domains)(mesh_facet);
	      if (domain < ufc.form.num_exterior_facet_integrals())
		integral = ufc.exterior_facet_integrals[domain];
	      else
		continue;
	    }

	  // Get mesh cell to which mesh facet belongs 
	  // (pick first, there is only one)
	  dolfin_assert(mesh_facet.numEntities(mesh.topology().dim()) == 1);
	  Cell mesh_cell(mesh, mesh_facet.entities(mesh.topology().dim())[0]);

	  // Get local index of facet with respect to the cell
	  const uint local_facet = mesh_cell.index(mesh_facet);
      
	  // Update to current cell
	  ufc.update(mesh_cell);

	  // Interpolate coefficients on cell
	  for (uint i = 0; i < coefficients.size(); i++)
	    coefficients[i]->interpolate(ufc.w[i], ufc.cell, 
					 *ufc.coefficient_elements[i], 
					 mesh_cell, local_facet);

	  // Tabulate dofs for each dimension
	  for (uint i = 0; i < ufc.form.rank(); i++)
	    dof_map_set[i].tabulate_dofs(ufc.dofs[i], ufc.cell, 
					 mesh_cell.index());

	  // Tabulate exterior facet tensor
	  integral->tabulate_tensor(ufc.A, ufc.w, ufc.cell, local_facet);
    
	  // Add entries to global tensor
	  A.add(ufc.A, ufc.local_dimensions, ufc.dofs);

	  p++;  
	}
    }

    // Assemble over interior facets
    template<typename Tensor_t>
    void assembleInteriorFacets(Tensor_t& A,
                                const Array<Function*>& coefficients,
                                const DofMapSet& dof_map_set,
                                UFC& ufc,
                                const MeshFunction<uint>* domains) const
    {
      // Skip assembly if there are no interior facet integrals
      if (ufc.form.num_interior_facet_integrals() == 0)
	return;
  
      // Interior facet integral
      ufc::interior_facet_integral* integral = ufc.interior_facet_integrals[0];

      // Compute facets and facet - cell connectivity if not already computed
      mesh.init(mesh.topology().dim() - 1);
      mesh.init(mesh.topology().dim() - 1, mesh.topology().dim());
      mesh.order();
  
      // Assemble over interior facets (the facets of the mesh)
      Progress p("Assembling over interior facets", mesh.numFacets());
      for (FacetIterator facet(mesh); !facet.end(); ++facet)
	{
	  // Check if we have an interior facet
	  if ( facet->numEntities(mesh.topology().dim()) != 2 )
	    {
	      p++;
	      continue;
	    }

	  // Get integral for sub domain (if any)
	  if (domains && domains->size() > 0)
	    {
	      const uint domain = (*domains)(*facet);
	      if (domain < ufc.form.num_interior_facet_integrals())
		integral = ufc.interior_facet_integrals[domain];
	      else
		continue;
	    }

	  // Get cells incident with facet
	  Cell cell0(mesh, facet->entities(mesh.topology().dim())[0]);
	  Cell cell1(mesh, facet->entities(mesh.topology().dim())[1]);
      
	  // Get local index of facet with respect to each cell
	  uint facet0 = cell0.index(*facet);
	  uint facet1 = cell1.index(*facet);

	  // Update to current pair of cells
	  ufc.update(cell0, cell1);
    
	  // Interpolate coefficients on cell
	  for (uint i = 0; i < coefficients.size(); i++)
	    {
	      const uint offset=ufc.coefficient_elements[i]->space_dimension();
	      coefficients[i]->interpolate(ufc.macro_w[i], ufc.cell0, 
					   *ufc.coefficient_elements[i], 
					   cell0, facet0);
	      coefficients[i]->interpolate(ufc.macro_w[i] + offset, ufc.cell1, 
					   *ufc.coefficient_elements[i], 
					   cell1, facet1);
	    }

	  // Tabulate dofs for each dimension on macro element
	  for (uint i = 0; i < ufc.form.rank(); i++)
	    {
	      const uint offset = ufc.local_dimensions[i];
	      dof_map_set[i].tabulate_dofs(ufc.macro_dofs[i],
					   ufc.cell0, cell0.index());
	      dof_map_set[i].tabulate_dofs(ufc.macro_dofs[i] + offset, 
					   ufc.cell1, cell1.index());
	    }

	  // Tabulate exterior interior facet tensor on macro element
	  integral->tabulate_tensor(ufc.macro_A, ufc.macro_w, ufc.cell0, 
				    ufc.cell1, facet0, facet1);

	  // Add entries to global tensor
	  A.add(ufc.macro_A, ufc.macro_local_dimensions, ufc.macro_dofs);

	  p++;
	}
    }

    // Check arguments
    void check(const ufc::form& form,
               const Array<Function*>& coefficients,
               const Mesh& mesh) const;

    // Initialize global tensor
    void initGlobalTensor(GenericTensor& A, const DofMapSet& dof_map_set, 
			  UFC& ufc, bool reset_tensor) const;

    // Pretty-printing for progress bar
    std::string progressMessage(uint rank, std::string integral_type) const;

    // The mesh
    Mesh& mesh;


  };
}

#endif
