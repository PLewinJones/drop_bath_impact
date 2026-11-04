// Driver for axisymmetric drop impact onto a liquid bath.

// Flag for debug output from the assign_element function
//#define FULL_ELASSIGN_OUTPUT

// Generic oomph-lib header
#include "generic.h"

// Navier-Stokes headers
#include "axisym_navier_stokes.h"

// Interface headers
#include "fluid_interface.h"

// Constitutive law headers
#include "constitutive.h"

// Solid headers
#include "solid.h"

// The mesh headers
#include "meshes/triangle_mesh.h"

// My element headers
#include "drop_bath_curved_elements.h"

using namespace std;

using namespace oomph;

//==start_of_namespace====================================================
/// Namespace for physical parameters
//========================================================================
namespace Global_Physical_Variables {
/// Reynolds number
double Re = 10.0;

/// Product of Reynolds number and inverse of Froude number (Gravity Number)
double ReInvFr = 0.0;

/// Capillary number
double Ca = 0.1;

/// Viscosity Ratio
double Vis_R = 0.0016;

///  Dim-less Hamaker number
double A_h = 0.0;

// Gas atmospheric pressure
// For incompressible flow, this can be set to 0.
double vapour_pressure_0 = 0.0;

/// Initial position for drop centre
double Drop_start_height = 1.1;

/// Initial Drop Radius
double Drop_radius = 1.0;

/// Initial speed in z direction (toward wall/sym plane)
double Initial_speed = 1.0;

/// Bath depth
double Bath_depth = 1.0;

/// Bath radius
double Bath_radius = 10.0;

/// Maximum time
double t_max = 20.0;

/// Use Gas Kinetic Effects?
bool Use_GKE = true;

/// Non-dim Mean free path of gas
double Mean_Free_Path = 8.625e-5;

/// Never need to change the rest of these

/// Strouhal number
double St = 1.0;

/// Womersley number (Reynolds x Strouhal, computed automatically)
double ReSt;

/// Direction of gravity
Vector<double> G(3);

/// Pseudo-solid Poisson ratio
double Nu = 0.2;

// GKE fuctions from Chubynsky et al 2020 for alpha=1
void GKE_Tau_P_fct(const double &Kn, double &result) { result = 1.0; }

void GKE_Tau_C_fct(const double &Kn, double &result) {
  double delta = sqrt(MathematicalConstants::Pi) / (2.0 * Kn);
  result = 1.0 + 2.0 * Kn *
                     (1.1466 - 0.095 * exp(-0.747 * delta) -
                      0.0516 * exp(-3.724 * delta));
}

void GKE_Phi_P_fct(const double &Kn, double &result) {
  if (Kn > 0.0) {
    result = 1.0 + 6.88 * Kn +
             (6.0 * Kn / MathematicalConstants::Pi) *
                 log1p(2.76 * Kn + 0.127 * Kn * Kn);
  } else {
    // Need to allow Kn<0 to stop a crash during adaptive timestepping when h is
    // very small
    result = 1.0 - 6.88 * Kn -
             (6.0 * Kn / MathematicalConstants::Pi) *
                 log1p(-2.76 * Kn + 0.127 * Kn * Kn);
  }
}

void GKE_Phi_C_fct(const double &Kn, double &result) { result = 1.0; }

// Alternative GKE function from Zhang and Law 2011
void GKE_Phi_P_ZL_fct(const double &Kn, double &result) {
  if (Kn < 1.0) {
    result = 1.0 + 6.0966 * Kn + 0.9650 * Kn * Kn + 0.6967 * Kn * Kn * Kn;
  } else {
    result = 8.7583 * pow(Kn, 1.1551);
  }
}

} // namespace Global_Physical_Variables

//==start_of_namespace====================================================
/// Namespace for Simulation Settings
//========================================================================
namespace Global_Sim_Settings {
/// Result Folder
string result_folder = "RESLT";

/// Restart file

string Restart_file = "";

// Minimum allowed height, simulation stopped if reached.
double Min_height_cutoff = 1.0e-5;

// Maximum timestep
double max_dt = 0.0005;

// Maximum timestep during impact
double max_dt_impact = 0.0001;

// Minimum timestep
double min_dt = 1.0e-12;

// target error for adaptive timestepping
double epsilon_t = 1.0e-5;

// Max element area for initial mesh
double initial_max_element_area = 0.01;

// Number of points on free surface in initial mesh
unsigned initial_n_surface_points = 100;

// Number of timesteps per adaption
unsigned steps_per_adapt = 100;

// Number of steps per full output - set very large to supress full outputs by
// default
unsigned steps_per_full_output = 10000000;

// Number of timesteps per restart file saved
unsigned steps_per_restart = 100;

// Start of Impact (where max_dt_impact is used)
double impact_start_time = 0.09;

// End of impact
double impact_end_time = 0.17;

// Targets for spatial adaptivity
// Max error (elements refined if error is above this)
double max_permitted_error = 5.0e-5;
// Min error (elements unrefined if error is below this)
double min_permitted_error = 1.0e-5;
// Max element size
double max_element_size = 0.1;
// Min element size
double min_element_size = 1.0e-7;

// Tolerance for refining boundary based on curvature
double polyline_refinement_tolerance = 0.002; // Default 0.08
// Tolerance for unrefining boundary based on curvature
double polyline_unrefinement_tolerance = 0.001; // Default 0.04

double drop_surface_max_length = 0.001;

double bath_surface_max_length = 0.01;


// Maximum film thickness
double maximum_film_thickness = 0.3;

// Tolerance for refining boundary based on curvature
double farfield_polyline_refinement_tolerance = 0.08; // Default 0.08
// Tolerance for unrefining boundary based on curvature
double farfield_polyline_unrefinement_tolerance = 0.04; // Default 0.04
// Target area for the far field region
double farfield_target_area = 0.05;

// radial position of the far field area
double farfield_radius = 2.5; //

} // End of namespace Global_Sim_Settings

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

//==start_of_problem_class================================================
/// Impact of a drop onto a bath
//========================================================================
template <class ELEMENT, class DROP_INTERFACE_ELEMENT,
          class BATH_INTERFACE_ELEMENT>
class ImpactProblem : public Problem {
public:
  /// Constructor
  ImpactProblem();

  /// Destructor
  ~ImpactProblem() {
    delete drop_outer_boundary_circle_pt;

    for (unsigned i = 0; i < 2; i++) {
      delete drop_boundary_polyline_pt[i];
    }
    drop_boundary_polyline_pt.clear();
    delete drop_boundary_polygon_pt;
    delete Drop_bulk_mesh_pt;
    delete Drop_surface_mesh_pt;
    delete Bath_surface_mesh_pt;
    delete drop_error_estimator_pt;
    delete Constitutive_law_pt;
    delete External_pressure_data_pt;
    delete this->time_stepper_pt();

    // Close trace file
    Trace_file.close();
  }

  /// Set initial conditions
  void set_initial_condition();

  /// Document the solution
  void doc_solution(bool full_output, bool post_adapt);

  // Function to find minimum height
  double min_height();

  /// Global error norm for adaptive time-stepping
  double global_temporal_error_norm();

  /// Create free drop surface elements
  void assign_external_elements();

  /// Suggestion for next timestep (stored to allow it to be written
  /// to (or read from) restart file)
  double &next_dt() { return Next_dt; }

  /// Return DocInfo object
  DocInfo &doc_info() { return Doc_info; }

  unsigned &adaptcount() { return Adaptcount; }

  unsigned &full_out_count() { return Full_out_count; }

  unsigned &restart_count() { return Restart_count; }

  /// Dump problem data to allow for later restart
  void dump_it(ofstream &dump_file);

  /// Write header for trace file
  void write_trace_file_header();

  /// Restart
  void restart(ifstream &restart_file);

private:
  /// Actions required before solve step (empty)
  void actions_before_newton_solve() {}

  /// No actions required after solve step (empty)
  void actions_after_newton_solve() {}

  /// Actions before adapt: Wipe the mesh of boundary elements
  void actions_before_adapt() {

    // Kill the  elements and wipe surface mesh
    delete_drop_boundary_elements();

    // Kill the  elements and wipe surface mesh
    delete_bath_boundary_elements();

    // Rebuild the Problem's global mesh from its various sub-meshes
    this->rebuild_global_mesh();

  } // end of actions_before_adapt

  /// Actions after adapt: Rebuild the mesh of Lagrange multiplier elements
  void actions_after_adapt() {
    // Reassign the lagrangian coordinates (an updated lagrangian approach)
    Drop_bulk_mesh_pt->set_lagrangian_nodal_coordinates();

    // Reassign the lagrangian coordinates (an updated lagrangian approach)
    Bath_bulk_mesh_pt->set_lagrangian_nodal_coordinates();

    // Create the boundary elements
    create_drop_boundary_elements();

    // Create the boundary elements
    create_bath_boundary_elements();

    // Rebuild the Problem's global mesh from its various sub-meshes
    this->rebuild_global_mesh();

    // Re-assign external elements
    // Need to do this before equation numbers are assigned
    oomph_info << "Assigning External Elements" << std::endl;
    assign_external_elements();

    // Setup the problem again -- remember that fluid mesh has been
    // completely rebuilt and its element's don't have any
    // pointers to Re etc. yet
    complete_problem_setup();

  } // end of actions_after_adapt

  /// \short Actions before the timestep: For maximum stability, reset
  /// the current nodal positions to be the "stress-free" ones.
  void actions_before_implicit_timestep() {
    Drop_bulk_mesh_pt->set_lagrangian_nodal_coordinates();
    Bath_bulk_mesh_pt->set_lagrangian_nodal_coordinates();
  }

  /// Create free drop surface elements
  void create_drop_boundary_elements();

  /// Create free drop surface elements
  void create_bath_boundary_elements();

  /// Delete boundary elements and wipe the associated mesh
  void delete_drop_boundary_elements() {
    // How many surface elements are in the surface mesh
    unsigned n_element = Drop_surface_mesh_pt->nelement();

    // Loop over the surface elements
    for (unsigned e = 0; e < n_element; e++) {
      // Kill surface element
      delete Drop_surface_mesh_pt->element_pt(e);
    }

    // Wipe the mesh
    Drop_surface_mesh_pt->flush_element_and_node_storage();

  } // end of delete_drop_boundary_elements

  /// Delete boundary elements and wipe the associated mesh
  void delete_bath_boundary_elements() {
    // How many surface elements are in the surface mesh
    unsigned n_element = Bath_surface_mesh_pt->nelement();

    // Loop over the surface elements
    for (unsigned e = 0; e < n_element; e++) {
      // Kill surface element
      delete Bath_surface_mesh_pt->element_pt(e);
    }

    // Wipe the mesh
    Bath_surface_mesh_pt->flush_element_and_node_storage();

  } // end of delete_bath_boundary_elements

  void complete_problem_setup();

  /// ----Drop Mesh pointers -----
  /// Pointer to the drop "bulk" mesh
  RefineableSolidTriangleMesh<ELEMENT> *Drop_bulk_mesh_pt;

  // Pointers for drop mesh construction
  Circle *drop_outer_boundary_circle_pt;
  Vector<TriangleMeshCurveSection *> drop_boundary_polyline_pt;
  TriangleMeshPolygon *drop_boundary_polygon_pt;
  TriangleMeshClosedCurve *drop_closed_curve_pt;

  // Pointer to drop error estimator
  Z2ErrorEstimator *drop_error_estimator_pt;

  // ids of the parts of the drop boundary
  // Free boundary at surface of drop
  const unsigned drop_free_boundary_id = 0;
  // Flat symmetry boundary at r=0
  const unsigned drop_fixed_boundary_id = 1;

  /// ----Bath Mesh pointers -----
  /// Pointer to the drop "bulk" mesh
  RefineableSolidTriangleMesh<ELEMENT> *Bath_bulk_mesh_pt;

  // Pointers for bath mesh construction
  Vector<TriangleMeshPolyLine *> bath_boundary_polyline_pt;
  Vector<TriangleMeshCurveSection *> bath_boundary_curvesection_pt;
  TriangleMeshPolygon *bath_boundary_polygon_pt;
  TriangleMeshClosedCurve *bath_closed_curve_pt;

  // Internal bath boundary pointers
  Vector<TriangleMeshOpenCurve *> inner_open_boundaries_pt;

  // Pointer to drop error estimator
  Z2ErrorEstimator *bath_error_estimator_pt;

  // ids of the parts of the drop boundary
  // Free boundary at surface of bath
  const unsigned bath_free_boundary_id = 0;
  // Flat symmetry boundary at r=0
  const unsigned bath_sym_boundary_id = 1;
  // Flat base boundary
  const unsigned bath_base_boundary_id = 2;
  // Flat far base boundary
  const unsigned bath_farbase_boundary_id = 3;
  // Flat wall boundary
  const unsigned bath_wall_boundary_id = 4;
  // Far free boundary
  const unsigned bath_farfree_boundary_id = 5;
  // Internal boundary
  const unsigned bath_internal_boundary_id = 6;

  // -----Other pointers
  /// Pointer to the  drop "surface" mesh
  Mesh *Drop_surface_mesh_pt;
  /// Storage for the back up the surface mesh
  BackupMeshForProjection<TElement<1, 3>> *Backed_up_drop_surface_mesh_pt;

  /// Pointer to the  bath "surface" mesh
  Mesh *Bath_surface_mesh_pt;

  // Pointer to the external pressure
  Data *External_pressure_data_pt;

  // Pointer to the constitutive law used to determine the mesh deformation
  ConstitutiveLaw *Constitutive_law_pt;

  /// Suggested next dt (store in problem for restarting)
  double Next_dt;

  /// Doc info object
  DocInfo Doc_info;

  /// Trace file
  ofstream Trace_file;

  unsigned Adaptcount;
  unsigned Full_out_count;
  unsigned Restart_count;

}; // End of problem class

//==start_of_constructor==================================================
/// Constructor for single fluid free surface problem
//========================================================================
template <class ELEMENT, class DROP_INTERFACE_ELEMENT,
          class BATH_INTERFACE_ELEMENT>
ImpactProblem<ELEMENT, DROP_INTERFACE_ELEMENT,
              BATH_INTERFACE_ELEMENT>::ImpactProblem() {
  // Allocate the timestepper (this constructs the time object as well)
  add_time_stepper_pt(new BDF<2>(true));
#ifdef OOMPH_HAS_MUMPS
  // Use mumps if available
  linear_solver_pt() = new MumpsSolver;
#endif

  // Set output directory
  Doc_info.set_directory(Global_Sim_Settings::result_folder);

  // Initialise counter for solutions
  Doc_info.number() = 0;

  Adaptcount = 1;
  Full_out_count = 1;
  Restart_count = 2; // Take restart snapshot before adapting!

  // Don't reject timesteps above tolerance
  // Keep_temporal_error_below_tolerance = false;

  // Timestep safety factor
  target_error_safety_factor() = 0.9;

  // Allow higher max residuals
  max_residuals() = 100000;

  // Allow more newton solves
  //max_newton_iterations() = 20;

  // Temp storage for coordinate along surface
  Vector<double> zeta(1);

  // temp storage for position vector on GeomObject
  Vector<double> posn(2);

  //  Build the drop mesh
  //===================

  // Circle object used to define the drop free outer boundary
  drop_outer_boundary_circle_pt =
      new Circle(0.0, Global_Physical_Variables::Drop_start_height,
                 Global_Physical_Variables::Drop_radius);
  // Pointer to the closed curve that defines the outer boundary
  drop_closed_curve_pt = 0;

  // Provide storage for pointers to the two parts of the curvilinear boundary
  drop_boundary_polyline_pt.resize(2);

  // Free boundary
  // Define start and end angle, anticlockwise
  double zeta_start = -0.5 * MathematicalConstants::Pi;
  double zeta_end = 0.5 * MathematicalConstants::Pi;
  // Contruct the points on the free boundary
  unsigned npoints = Global_Sim_Settings::initial_n_surface_points;
  double unit_zeta = (zeta_end - zeta_start) / double(npoints - 1);
  Vector<Vector<double>> free_boundary_vertex(npoints);
  for (unsigned ipoint = 0; ipoint < npoints; ipoint++) {
    // Resize the vector
    free_boundary_vertex[ipoint].resize(2);

    // Get the coordinates
    zeta[0] = zeta_start + unit_zeta * double(ipoint);
    drop_outer_boundary_circle_pt->position(zeta, posn);
    free_boundary_vertex[ipoint] = posn;
  }
  drop_boundary_polyline_pt[0] =
      new TriangleMeshPolyLine(free_boundary_vertex, drop_free_boundary_id);

  drop_boundary_polyline_pt[0]->set_maximum_length(
      Global_Sim_Settings::drop_surface_max_length);

  // Flat boundary at r=0
  unsigned n_poly_points = 2;
  Vector<Vector<double>> vertex_coord(n_poly_points);
  for (unsigned i = 0; i < n_poly_points; i++) {
    vertex_coord[i].resize(2);
  }
  // Start at end of free surface
  zeta[0] = zeta_end;
  drop_outer_boundary_circle_pt->position(zeta, posn);
  vertex_coord[0][0] = posn[0];
  vertex_coord[0][1] = posn[1];
  // End at start of free surface
  zeta[0] = zeta_start;
  drop_outer_boundary_circle_pt->position(zeta, posn);
  vertex_coord[1][0] = posn[0];
  vertex_coord[1][1] = posn[1];

  // Build the 2nd boundary polyline
  drop_boundary_polyline_pt[1] =
      new TriangleMeshPolyLine(vertex_coord, drop_fixed_boundary_id);

  // Combine to curvilinear boundary
  //--------------------------------
  drop_boundary_polygon_pt = new TriangleMeshPolygon(drop_boundary_polyline_pt);
  // Set the refinement tolerances
  drop_boundary_polygon_pt->set_polyline_refinement_tolerance(
      Global_Sim_Settings::polyline_refinement_tolerance);

  drop_boundary_polygon_pt->set_polyline_unrefinement_tolerance(
      Global_Sim_Settings::polyline_unrefinement_tolerance);

  // Set the pointer
  drop_closed_curve_pt = drop_boundary_polygon_pt;

  // Use the TriangleMeshParameters object for helping on the manage of the
  // TriangleMesh parameters
  TriangleMeshParameters drop_triangle_mesh_parameters(drop_closed_curve_pt);

  // Specify the maximum area element
  drop_triangle_mesh_parameters.element_area() =
      Global_Sim_Settings::initial_max_element_area;

  // Store as the problem's main "Bulk" mesh
  Drop_bulk_mesh_pt = new RefineableSolidTriangleMesh<ELEMENT>(
      drop_triangle_mesh_parameters, this->time_stepper_pt());

  // Set error estimator for bulk mesh
  drop_error_estimator_pt = new Z2ErrorEstimator;
  Drop_bulk_mesh_pt->spatial_error_estimator_pt() = drop_error_estimator_pt;

  // Set targets for spatial adaptivity
  Drop_bulk_mesh_pt->max_permitted_error() =
      Global_Sim_Settings::max_permitted_error;
  Drop_bulk_mesh_pt->min_permitted_error() =
      Global_Sim_Settings::min_permitted_error;
  Drop_bulk_mesh_pt->max_element_size() = Global_Sim_Settings::max_element_size;
  Drop_bulk_mesh_pt->min_element_size() = Global_Sim_Settings::min_element_size;
  Drop_bulk_mesh_pt->set_print_level_timings_adaptation(3);

  //  Build the bath mesh
  //===================

  // Pointer to the closed curve that defines the outer boundary
  bath_closed_curve_pt = 0;

  // Provide storage for pointers to the parts of the curvilinear boundary
  bath_boundary_curvesection_pt.resize(6);
  bath_boundary_polyline_pt.resize(6);
  if (Global_Physical_Variables::Bath_radius <=
      Global_Sim_Settings::farfield_radius) {
    oomph_info << "Error, Far field radius must be less than Bath radius!"
               << std::endl;
  }
  // Flat boundary
  n_poly_points = 2;
  vertex_coord.resize(n_poly_points);
  for (unsigned i = 0; i < n_poly_points; i++) {
    vertex_coord[i].resize(2);
  }
  // start
  vertex_coord[0][0] = Global_Sim_Settings::farfield_radius; // x
  vertex_coord[0][1] = 0.0;                                  // y
  // end
  vertex_coord[1][0] = 0.0; // x
  vertex_coord[1][1] = 0.0; // y

  // Build the boundary polyline
  bath_boundary_polyline_pt[bath_free_boundary_id] =
      new TriangleMeshPolyLine(vertex_coord, bath_free_boundary_id);
  bath_boundary_curvesection_pt[bath_free_boundary_id] =
      bath_boundary_polyline_pt[bath_free_boundary_id];

  // Set the refinement tolerances
  bath_boundary_curvesection_pt[bath_free_boundary_id]
      ->set_refinement_tolerance(
          Global_Sim_Settings::polyline_refinement_tolerance);

  bath_boundary_curvesection_pt[bath_free_boundary_id]
      ->set_unrefinement_tolerance(
          Global_Sim_Settings::polyline_unrefinement_tolerance);

  bath_boundary_curvesection_pt[bath_free_boundary_id]->set_maximum_length(
      Global_Sim_Settings::bath_surface_max_length);
  // start
  vertex_coord[0][0] = 0.0; // x
  vertex_coord[0][1] = 0.0; // y
  // end
  vertex_coord[1][0] = 0.0;                                    // x
  vertex_coord[1][1] = -Global_Physical_Variables::Bath_depth; // y

  // Build the boundary polyline
  bath_boundary_polyline_pt[bath_sym_boundary_id] =
      new TriangleMeshPolyLine(vertex_coord, bath_sym_boundary_id);
  bath_boundary_curvesection_pt[bath_sym_boundary_id] =
      bath_boundary_polyline_pt[bath_sym_boundary_id];

  // start
  vertex_coord[0][0] = 0.0;                                    // x
  vertex_coord[0][1] = -Global_Physical_Variables::Bath_depth; // y
  // end
  vertex_coord[1][0] = Global_Sim_Settings::farfield_radius;   // x
  vertex_coord[1][1] = -Global_Physical_Variables::Bath_depth; // y

  // Build the boundary polyline
  bath_boundary_polyline_pt[bath_base_boundary_id] =
      new TriangleMeshPolyLine(vertex_coord, bath_base_boundary_id);
  bath_boundary_curvesection_pt[bath_base_boundary_id] =
      bath_boundary_polyline_pt[bath_base_boundary_id];

  // start
  vertex_coord[0][0] = Global_Sim_Settings::farfield_radius;   // x
  vertex_coord[0][1] = -Global_Physical_Variables::Bath_depth; // y
  // end
  vertex_coord[1][0] = Global_Physical_Variables::Bath_radius; // x
  vertex_coord[1][1] = -Global_Physical_Variables::Bath_depth; // y

  // Build the boundary polyline
  bath_boundary_polyline_pt[bath_farbase_boundary_id] =
      new TriangleMeshPolyLine(vertex_coord, bath_farbase_boundary_id);
  bath_boundary_curvesection_pt[bath_farbase_boundary_id] =
      bath_boundary_polyline_pt[bath_farbase_boundary_id];

  // start
  vertex_coord[0][0] = Global_Physical_Variables::Bath_radius; // x
  vertex_coord[0][1] = -Global_Physical_Variables::Bath_depth; // y
  // end
  vertex_coord[1][0] = Global_Physical_Variables::Bath_radius; // x
  vertex_coord[1][1] = 0.0;                                    // y

  // Build the boundary polyline
  bath_boundary_polyline_pt[bath_wall_boundary_id] =
      new TriangleMeshPolyLine(vertex_coord, bath_wall_boundary_id);
  bath_boundary_curvesection_pt[bath_wall_boundary_id] =
      bath_boundary_polyline_pt[bath_wall_boundary_id];

  // start
  vertex_coord[0][0] = Global_Physical_Variables::Bath_radius; // x
  vertex_coord[0][1] = 0.0;                                    // y
  // end
  vertex_coord[1][0] = Global_Sim_Settings::farfield_radius; // x
  vertex_coord[1][1] = 0.0;                                  // y

  // Build the boundary polyline
  bath_boundary_polyline_pt[bath_farfree_boundary_id] =
      new TriangleMeshPolyLine(vertex_coord, bath_farfree_boundary_id);
  bath_boundary_curvesection_pt[bath_farfree_boundary_id] =
      bath_boundary_polyline_pt[bath_farfree_boundary_id];

  // Set the refinement tolerances
  bath_boundary_curvesection_pt[bath_farfree_boundary_id]
      ->set_refinement_tolerance(
          Global_Sim_Settings::farfield_polyline_refinement_tolerance);

  bath_boundary_curvesection_pt[bath_farfree_boundary_id]
      ->set_unrefinement_tolerance(
          Global_Sim_Settings::farfield_polyline_unrefinement_tolerance);

  // Combine to curvilinear boundary
  //--------------------------------
  bath_boundary_polygon_pt =
      new TriangleMeshPolygon(bath_boundary_curvesection_pt);

  // Create internal bath boundary
  // Only one internal boundary, but needs to be a vector.
  inner_open_boundaries_pt.resize(1);

  n_poly_points = 2;
  vertex_coord.resize(n_poly_points);
  for (unsigned i = 0; i < n_poly_points; i++) {
    vertex_coord[i].resize(2);
  }

  // start
  vertex_coord[0][0] = Global_Sim_Settings::farfield_radius;   // x
  vertex_coord[0][1] = -Global_Physical_Variables::Bath_depth; // y
  // end
  vertex_coord[1][0] = Global_Sim_Settings::farfield_radius; // x // x
  vertex_coord[1][1] = 0.0;                                  // y

  TriangleMeshPolyLine *bath_internal_boundary_pt =
      new TriangleMeshPolyLine(vertex_coord, bath_internal_boundary_id);

  // Do the connection with the destination boundary
  bath_internal_boundary_pt->connect_initial_vertex_to_polyline(
      bath_boundary_polyline_pt[bath_farbase_boundary_id], 0);

  // Do the connection with the destination boundary
  bath_internal_boundary_pt->connect_final_vertex_to_polyline(
      bath_boundary_polyline_pt[bath_farfree_boundary_id], 1);

  Vector<TriangleMeshCurveSection *> internal_curve_section_pt(1);
  internal_curve_section_pt[0] = bath_internal_boundary_pt;
  inner_open_boundaries_pt[0] =
      new TriangleMeshOpenCurve(internal_curve_section_pt);

  // Set the pointer
  bath_closed_curve_pt = bath_boundary_polygon_pt;

  // Use the TriangleMeshParameters object for helping on the manage of the
  // TriangleMesh parameters
  TriangleMeshParameters bath_triangle_mesh_parameters(bath_closed_curve_pt);

  // Specify the internal open boundaries
  bath_triangle_mesh_parameters.internal_open_curves_pt() =
      inner_open_boundaries_pt;

  // Specify the maximum area element
  bath_triangle_mesh_parameters.element_area() =
      Global_Sim_Settings::initial_max_element_area;

  // Add a region on the domain
  Vector<double> region(2);

  // Define the coordinates of the regions on the domain
  region[0] = 0.5 * (Global_Physical_Variables::Bath_radius +
                     Global_Sim_Settings::farfield_radius);
  region[1] = -0.5 * Global_Physical_Variables::Bath_depth;

  // Pass information about the defined regions
  bath_triangle_mesh_parameters.add_region_coordinates(1, region);

  // Specify different target area of region 1:
  bath_triangle_mesh_parameters.set_target_area_for_region(
      1, Global_Sim_Settings::farfield_target_area);

  // Store as the problem's main "Bulk" mesh
  Bath_bulk_mesh_pt = new RefineableSolidTriangleMesh<ELEMENT>(
      bath_triangle_mesh_parameters, this->time_stepper_pt());

  // Set error estimator for bulk mesh
  bath_error_estimator_pt = new Z2ErrorEstimator;
  Bath_bulk_mesh_pt->spatial_error_estimator_pt() = bath_error_estimator_pt;

  // Set targets for spatial adaptivity
  Bath_bulk_mesh_pt->max_permitted_error() =
      Global_Sim_Settings::max_permitted_error;
  Bath_bulk_mesh_pt->min_permitted_error() =
      Global_Sim_Settings::min_permitted_error;
  Bath_bulk_mesh_pt->max_element_size() = Global_Sim_Settings::max_element_size;
  Bath_bulk_mesh_pt->min_element_size() = Global_Sim_Settings::min_element_size;
  Bath_bulk_mesh_pt->set_print_level_timings_adaptation(3);

  // -----------------------------
  // Create the drop surface
  // -----------------------------
  // Create the "surface mesh" that will contain only the interface
  // elements. The constructor just creates the mesh without giving
  // it any elements, nodes, etc.
  Drop_surface_mesh_pt = new Mesh;

  create_drop_boundary_elements();

  // -----------------------------
  // Create the bath surface
  // -----------------------------
  // Create the "surface mesh" that will contain only the interface
  // elements. The constructor just creates the mesh without giving
  // it any elements, nodes, etc.
  Bath_surface_mesh_pt = new Mesh;
  oomph_info << "Creating Boundary Elements" << std::endl;
  create_bath_boundary_elements();

  // Add the sub-meshes to the problem
  add_sub_mesh(Drop_bulk_mesh_pt);
  add_sub_mesh(Drop_surface_mesh_pt);
  add_sub_mesh(Bath_bulk_mesh_pt);
  add_sub_mesh(Bath_surface_mesh_pt);

  // Combine all sub-meshes into a single mesh
  build_global_mesh();

  oomph_info << "Assigning External Elements" << std::endl;
  assign_external_elements();

  // Define a constitutive law for the solid equations: generalised Hookean
  Constitutive_law_pt = new GeneralisedHookean(&Global_Physical_Variables::Nu);

  // Create a Data object whose single value stores the dummy external pressure
  External_pressure_data_pt = new Data(1);
  External_pressure_data_pt->pin(0);
  External_pressure_data_pt->set_value(0, 0.0);

  maximum_dt() = Global_Sim_Settings::max_dt;
  minimum_dt() = Global_Sim_Settings::min_dt;

  complete_problem_setup();

  // Setup equation numbering scheme
  oomph_info << "Number of equations: " << assign_eqn_numbers() << std::endl;

} // End of constructor

//========start_of_write_trace_file_header==============================
/// Write the header of the trace file
//=========================================================================
template <class ELEMENT, class DROP_INTERFACE_ELEMENT,
          class BATH_INTERFACE_ELEMENT>
void ImpactProblem<ELEMENT, DROP_INTERFACE_ELEMENT,
                   BATH_INTERFACE_ELEMENT>::write_trace_file_header() {

  char filename[100];
  sprintf(filename, "%s/trace.dat", doc_info().directory().c_str());
  Trace_file.open(filename);
  Trace_file.precision(16);

  // Print Trace file header
  Trace_file << "---------------------------------------------" << std::endl;
  Trace_file << "---------- Leidenfrost Drop Impact ----------" << std::endl;
  Trace_file << "---------------------------------------------" << std::endl;
  Trace_file << "Output Folder:       " << Global_Sim_Settings::result_folder
             << std::endl;
  Trace_file << "Reynolds Number:     " << Global_Physical_Variables::Re
             << std::endl;
  Trace_file << "Gravity Number:      " << Global_Physical_Variables::ReInvFr
             << std::endl;
  Trace_file << "Capilliary Number:   " << Global_Physical_Variables::Ca
             << std::endl;
  Trace_file << "Viscosity Ratio:     " << Global_Physical_Variables::Vis_R
             << std::endl;
  Trace_file << "NDim Hamaker Number: " << Global_Physical_Variables::A_h
             << std::endl;
  Trace_file << "Ndim Mean Free Path: "
             << Global_Physical_Variables::Mean_Free_Path << std::endl;
  Trace_file << "Using GKE:           " << Global_Physical_Variables::Use_GKE
             << std::endl;
  Trace_file << "Maximum time:        " << Global_Physical_Variables::t_max
             << std::endl;
  Trace_file << "Centre Start Height: "
             << Global_Physical_Variables::Drop_start_height << std::endl;
  Trace_file << "Initial Drop Radius: "
             << Global_Physical_Variables::Drop_radius << std::endl;
  Trace_file << "Initial Speed:       "
             << Global_Physical_Variables::Initial_speed << std::endl;
  Trace_file << "Bath Depth:          " << Global_Physical_Variables::Bath_depth
             << std::endl;
  Trace_file << "Bath Radius:         "
             << Global_Physical_Variables::Bath_radius << std::endl;
  Trace_file << "Min allowed height:  "
             << Global_Sim_Settings::Min_height_cutoff << std::endl;
  Trace_file << "Max Timestep:        " << Global_Sim_Settings::max_dt
             << std::endl;
  Trace_file << "Max Timestep Impact: " << Global_Sim_Settings::max_dt_impact
             << std::endl;
  Trace_file << "Min Timestep:        " << Global_Sim_Settings::min_dt
             << std::endl;
  Trace_file << "Target Time Error:   " << Global_Sim_Settings::epsilon_t
             << std::endl;
  Trace_file << "Start Element Area:  "
             << Global_Sim_Settings::initial_max_element_area << std::endl;
  Trace_file << "Init Surface points: "
             << Global_Sim_Settings::initial_n_surface_points << std::endl;
  Trace_file << "Steps per adapt:     " << Global_Sim_Settings::steps_per_adapt
             << std::endl;
  Trace_file << "Steps per full out:  "
             << Global_Sim_Settings::steps_per_full_output << std::endl;
  Trace_file << "Max mesh error:      "
             << Global_Sim_Settings::max_permitted_error << std::endl;
  Trace_file << "Min mesh error:      "
             << Global_Sim_Settings::min_permitted_error << std::endl;
  Trace_file << "max element size:    " << Global_Sim_Settings::max_element_size
             << std::endl;
  Trace_file << "min element size:    " << Global_Sim_Settings::min_element_size
             << std::endl;
  Trace_file << "bound refine tol:    "
             << Global_Sim_Settings::polyline_refinement_tolerance << std::endl;
  Trace_file << "bound unrefine tol:  "
             << Global_Sim_Settings::polyline_unrefinement_tolerance
             << std::endl;
  Trace_file << "Impact Start Time:   "
             << Global_Sim_Settings::impact_start_time << std::endl;
  Trace_file << "Impact End Time:     " << Global_Sim_Settings::impact_end_time
             << std::endl;
  Trace_file << "D Surf Max Length:   "
             << Global_Sim_Settings::drop_surface_max_length << std::endl;
  Trace_file << "B Surf Max Length:   "
             << Global_Sim_Settings::bath_surface_max_length << std::endl;
  Trace_file << "Max Film Thickness:  "
             << Global_Sim_Settings::maximum_film_thickness << std::endl;
  Trace_file << "Far refine tol:      "
             << Global_Sim_Settings::farfield_polyline_refinement_tolerance
             << std::endl;
  Trace_file << "Far unefine tol:     "
             << Global_Sim_Settings::farfield_polyline_unrefinement_tolerance
             << std::endl;
  Trace_file << "Far Target Area:     "
             << Global_Sim_Settings::farfield_target_area << std::endl;
  Trace_file << "Far Radius:          " << Global_Sim_Settings::farfield_radius
             << std::endl;
  Trace_file << "---------------------------------------------" << std::endl;

  // Headers for trace file
  Trace_file << "VARIABLES=\"timestep\",\"time\",\"dt\",\"Global Error\",\"min "
                "z\",\"N Elements\""
             << std::endl;

} // End write_trace_file_header

//========start_of_global_temporal_error_norm==============================
/// Global error norm for adaptive timestepping: RMS error, based on
/// difference between predicted and actual value at all nodes.
//=========================================================================
template <class ELEMENT, class DROP_INTERFACE_ELEMENT,
          class BATH_INTERFACE_ELEMENT>
double ImpactProblem<ELEMENT, DROP_INTERFACE_ELEMENT,
                     BATH_INTERFACE_ELEMENT>::global_temporal_error_norm() {
  double global_position_error = 0.0;

  double global_velocity_error = 0.0;

  // Find out how many nodes there are in the problem
  unsigned n_node = mesh_pt()->nnode();

  // Loop over the nodes and calculate the estimated error in the values
  for (unsigned i = 0; i < n_node; i++) {
    // Get error in solution: Difference between predicted and actual
    // value for position value 0
    double error = mesh_pt()
                       ->node_pt(i)
                       ->position_time_stepper_pt()
                       ->temporal_error_in_position(mesh_pt()->node_pt(i), 0);

    // Add the square of the individual error to the global error
    global_position_error += error * error;

    // Get error in solution: Difference between predicted and actual
    // value for position value 1
    error = mesh_pt()
                ->node_pt(i)
                ->position_time_stepper_pt()
                ->temporal_error_in_position(mesh_pt()->node_pt(i), 1);

    // Add the square of the individual error to the global error
    global_position_error += error * error;

    error = mesh_pt()->node_pt(i)->time_stepper_pt()->temporal_error_in_value(
        mesh_pt()->node_pt(i), 0);

    global_velocity_error += error * error;

    error = mesh_pt()->node_pt(i)->time_stepper_pt()->temporal_error_in_value(
        mesh_pt()->node_pt(i), 1);

    global_velocity_error += error * error;
  }

  // Divide by the number of nodes
  global_position_error /= 2.0 * double(n_node);

  global_velocity_error /= 2.0 * double(n_node);

  // Return sum of square roots of errors
  return sqrt(global_position_error) + sqrt(global_velocity_error);
} // end of global_temporal_error_norm

//============start_of_create_drop_boundary_elements===============
/// Create  drop boundary elements
//=======================================================================
template <class ELEMENT, class DROP_INTERFACE_ELEMENT,
          class BATH_INTERFACE_ELEMENT>
void ImpactProblem<ELEMENT, DROP_INTERFACE_ELEMENT,
                   BATH_INTERFACE_ELEMENT>::create_drop_boundary_elements() {
  unsigned n_element = Drop_bulk_mesh_pt->nboundary_element(0);
  // Loop over those elements adjacent to the free surface
  for (unsigned e = 0; e < n_element; e++) {
    // Set a pointer to the bulk element we wish to our interface
    // element to
    ELEMENT *bulk_element_pt =
        dynamic_cast<ELEMENT *>(Drop_bulk_mesh_pt->boundary_element_pt(0, e));

    // Find the index of the face of element e along boundary b
    int face_index = Drop_bulk_mesh_pt->face_index_at_boundary(0, e);

    // Create the interface element
    DROP_INTERFACE_ELEMENT *interface_element_pt =
        new DROP_INTERFACE_ELEMENT(bulk_element_pt, face_index);

    // Add the interface element to the surface mesh
    this->Drop_surface_mesh_pt->add_element_pt(interface_element_pt);
    interface_element_pt->set_boundary_number_in_bulk_mesh(0);
  }

} // end of create_drop_boundary_elements

//============start_of_create_bath_boundary_elements===============
/// Create  bath boundary elements
//=======================================================================
template <class ELEMENT, class DROP_INTERFACE_ELEMENT,
          class BATH_INTERFACE_ELEMENT>
void ImpactProblem<ELEMENT, DROP_INTERFACE_ELEMENT,
                   BATH_INTERFACE_ELEMENT>::create_bath_boundary_elements() {
  unsigned n_element =
      Bath_bulk_mesh_pt->nboundary_element(bath_free_boundary_id);
  // Loop over those elements adjacent to the free surface
  for (unsigned e = 0; e < n_element; e++) {
    // Set a pointer to the bulk element we wish to our interface
    // element to
    ELEMENT *bulk_element_pt = dynamic_cast<ELEMENT *>(
        Bath_bulk_mesh_pt->boundary_element_pt(bath_free_boundary_id, e));

    // Find the index of the face of element e along boundary b
    int face_index =
        Bath_bulk_mesh_pt->face_index_at_boundary(bath_free_boundary_id, e);

    // Create the interface element
    BATH_INTERFACE_ELEMENT *interface_element_pt =
        new BATH_INTERFACE_ELEMENT(bulk_element_pt, face_index);

    // Add the interface element to the surface mesh
    this->Bath_surface_mesh_pt->add_element_pt(interface_element_pt);
    interface_element_pt->set_boundary_number_in_bulk_mesh(
        bath_free_boundary_id);
  }

  n_element = Bath_bulk_mesh_pt->nboundary_element(bath_farfree_boundary_id);
  // Loop over those elements adjacent to the free surface
  for (unsigned e = 0; e < n_element; e++) {
    // Set a pointer to the bulk element we wish to our interface
    // element to
    ELEMENT *bulk_element_pt = dynamic_cast<ELEMENT *>(
        Bath_bulk_mesh_pt->boundary_element_pt(bath_farfree_boundary_id, e));

    // Find the index of the face of element e along boundary b
    int face_index =
        Bath_bulk_mesh_pt->face_index_at_boundary(bath_farfree_boundary_id, e);

    // Create the interface element
    BATH_INTERFACE_ELEMENT *interface_element_pt =
        new BATH_INTERFACE_ELEMENT(bulk_element_pt, face_index);

    // Add the interface element to the surface mesh
    this->Bath_surface_mesh_pt->add_element_pt(interface_element_pt);
    interface_element_pt->set_boundary_number_in_bulk_mesh(
        bath_farfree_boundary_id);
  }

} // end of create_bath_boundary_elements

//============start_of_assign_external_elements===============
/// assign the external elements
//=======================================================================
template <class ELEMENT, class DROP_INTERFACE_ELEMENT,
          class BATH_INTERFACE_ELEMENT>
void ImpactProblem<ELEMENT, DROP_INTERFACE_ELEMENT,
                   BATH_INTERFACE_ELEMENT>::assign_external_elements() {

  oomph_info << "Sorting Drop Surface" << std::endl;
  // Create Sorted Drop Boundary
  // We want a list of the elements on the drop surface, starting from the
  // bottom. This is in the anticlockwise direction around the drop surface.
  // Adapted from
  // RefineableTriangleMesh<ELEMENT>::create_sorted_face_mesh_representation
  // Create a list of pointers to the element
  std::list<DROP_INTERFACE_ELEMENT *> Sorted_drop_surface_el_pt_list;
  // Upcast from GeneralisedElement to the present element
  DROP_INTERFACE_ELEMENT *drop_el_pt = dynamic_cast<DROP_INTERFACE_ELEMENT *>(
      Drop_surface_mesh_pt->finite_element_pt(0));
  Sorted_drop_surface_el_pt_list.push_back(drop_el_pt);

  // Number of nodes
  unsigned nnod = drop_el_pt->nnode();

  // Count elements that have been done
  unsigned drop_count_done = 0;
  // Keep track of who's done
  std::map<DROP_INTERFACE_ELEMENT *, bool> drop_done_el;
  // Keep track of inverted elements
  std::map<DROP_INTERFACE_ELEMENT *, bool> drop_el_is_inverted;

  // Assign drop elements to bath elements
  unsigned n_drop_interface_element = Drop_surface_mesh_pt->nelement();

  // Fit in the other elements in at most nel^2 loops
  for (unsigned ee = 1; ee < n_drop_interface_element; ee++) {
    // Loop over all elements to check if they fit to the right
    // or the left of the current one
    for (unsigned e = 1; e < n_drop_interface_element; e++) {
      // Candidate element
      drop_el_pt = dynamic_cast<DROP_INTERFACE_ELEMENT *>(
          Drop_surface_mesh_pt->finite_element_pt(e));

      // Is it done yet?
      if (!drop_done_el[drop_el_pt]) {
        // Left and rightmost elements
        DROP_INTERFACE_ELEMENT *first_el_pt =
            (*Sorted_drop_surface_el_pt_list.begin());
        typename std::list<DROP_INTERFACE_ELEMENT *>::iterator it =
            Sorted_drop_surface_el_pt_list.end();
        it--;
        DROP_INTERFACE_ELEMENT *last_el_pt = *it;

        // Left and rightmost nodes
        Node *left_node_pt = first_el_pt->node_pt(0);
        if (drop_el_is_inverted[first_el_pt]) {
          left_node_pt = first_el_pt->node_pt(nnod - 1);
        }
        Node *right_node_pt = last_el_pt->node_pt(nnod - 1);
        if (drop_el_is_inverted[last_el_pt]) {
          right_node_pt = last_el_pt->node_pt(0);
        }

        // New element fits at the left of first element and is not inverted
        if (left_node_pt == drop_el_pt->node_pt(nnod - 1)) {
          Sorted_drop_surface_el_pt_list.push_front(drop_el_pt);
          drop_done_el[drop_el_pt] = true;
          drop_count_done++;
          drop_el_is_inverted[drop_el_pt] = false;
        }
        // New element fits at the left of first element and is inverted

        else if (left_node_pt == drop_el_pt->node_pt(0)) {
          Sorted_drop_surface_el_pt_list.push_front(drop_el_pt);
          drop_done_el[drop_el_pt] = true;
          drop_count_done++;
          drop_el_is_inverted[drop_el_pt] = true;
        }
        // New element fits on the right of last element and is not inverted

        else if (right_node_pt == drop_el_pt->node_pt(0)) {
          Sorted_drop_surface_el_pt_list.push_back(drop_el_pt);
          drop_done_el[drop_el_pt] = true;
          drop_count_done++;
          drop_el_is_inverted[drop_el_pt] = false;
        }
        // New element fits on the right of last element and is inverted

        else if (right_node_pt == drop_el_pt->node_pt(nnod - 1)) {
          Sorted_drop_surface_el_pt_list.push_back(drop_el_pt);
          drop_done_el[drop_el_pt] = true;
          drop_count_done++;
          drop_el_is_inverted[drop_el_pt] = true;
        }

        if (drop_done_el[drop_el_pt]) {
          break;
        }
      }
    }
  }

  // Are we done?
  if (drop_count_done != (n_drop_interface_element - 1)) {
    std::ostringstream error_message;
    error_message
        << "When ordering Elements on Drop Surface, only managed to order \n"
        << drop_count_done << " of " << n_drop_interface_element
        << " face elements.\n"
        << std::endl;
    throw OomphLibError(error_message.str(), OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
  }

  // Now convert the list to a vector of points
  std::vector<DROP_INTERFACE_ELEMENT *> Sorted_drop_surface_el_pt;
  Sorted_drop_surface_el_pt.reserve(Sorted_drop_surface_el_pt_list.size());
  std::copy(std::begin(Sorted_drop_surface_el_pt_list),
            std::end(Sorted_drop_surface_el_pt_list),
            std::back_inserter(Sorted_drop_surface_el_pt));

  // Is this an anticlockwise list?
  DROP_INTERFACE_ELEMENT *first_drop_el_pt = Sorted_drop_surface_el_pt[0];
  Vector<double> s_temp(1, 0.5);
  double z_drop_first = first_drop_el_pt->interpolated_x(s_temp, 1);

  DROP_INTERFACE_ELEMENT *last_drop_el_pt =
      Sorted_drop_surface_el_pt[n_drop_interface_element - 1];
  double z_drop_last = last_drop_el_pt->interpolated_x(s_temp, 1);
  bool drop_list_anticlock = true;

  if (z_drop_first > z_drop_last) {
    drop_list_anticlock = false;
  }
#ifdef FULL_ELASSIGN_OUTPUT
  oomph_info << "Drop is anticlockwise:" << drop_list_anticlock << std::endl;
#endif
  unsigned drop_e_top = 0;
  unsigned drop_e_bottom = n_drop_interface_element - 1;

  if (drop_list_anticlock) {
    drop_e_top = n_drop_interface_element - 1;
    drop_e_bottom = 0;
  }
  // End Sort Drop Surface

  oomph_info << "Sorting Bath Surface" << std::endl;
  // Create Sorted Bath Boundary
  // We want a list of the elements on the bath surface, starting from the
  // outside, This is in the anticlockwise direction around the drop surface.
  // Adapted from
  // RefineableTriangleMesh<ELEMENT>::create_sorted_face_mesh_representation
  // Create a list of pointers to the element
  std::list<BATH_INTERFACE_ELEMENT *> Sorted_bath_surface_el_pt_list;
  BATH_INTERFACE_ELEMENT *bath_el_pt = dynamic_cast<BATH_INTERFACE_ELEMENT *>(
      Bath_surface_mesh_pt->finite_element_pt(0));
  Sorted_bath_surface_el_pt_list.push_back(bath_el_pt);

  // Count elements that have been done
  unsigned bath_count_done = 0;
  // Keep track of who's done
  std::map<BATH_INTERFACE_ELEMENT *, bool> bath_done_el;
  // Keep track of inverted elements
  std::map<BATH_INTERFACE_ELEMENT *, bool> bath_el_is_inverted;

  unsigned n_bath_interface_element = Bath_surface_mesh_pt->nelement();

  // Fit in the other elements in at most nel^2 loops
  for (unsigned ee = 1; ee < n_bath_interface_element; ee++) {
    // Loop over all elements to check if they fit to the right
    // or the left of the current one
    for (unsigned e = 1; e < n_bath_interface_element; e++) {
      // Candidate element
      bath_el_pt = dynamic_cast<BATH_INTERFACE_ELEMENT *>(
          Bath_surface_mesh_pt->finite_element_pt(e));

      // Is it done yet?
      if (!bath_done_el[bath_el_pt]) {
        // Left and rightmost elements
        BATH_INTERFACE_ELEMENT *first_el_pt =
            (*Sorted_bath_surface_el_pt_list.begin());
        typename std::list<BATH_INTERFACE_ELEMENT *>::iterator it =
            Sorted_bath_surface_el_pt_list.end();
        it--;
        BATH_INTERFACE_ELEMENT *last_el_pt = *it;

        // Left and rightmost nodes
        Node *left_node_pt = first_el_pt->node_pt(0);
        if (bath_el_is_inverted[first_el_pt]) {
          left_node_pt = first_el_pt->node_pt(nnod - 1);
        }
        Node *right_node_pt = last_el_pt->node_pt(nnod - 1);
        if (bath_el_is_inverted[last_el_pt]) {
          right_node_pt = last_el_pt->node_pt(0);
        }

        // New element fits at the left of first element and is not inverted
        if (left_node_pt == bath_el_pt->node_pt(nnod - 1)) {
          Sorted_bath_surface_el_pt_list.push_front(bath_el_pt);
          bath_done_el[bath_el_pt] = true;
          bath_count_done++;
          bath_el_is_inverted[bath_el_pt] = false;
        }
        // New element fits at the left of first element and is inverted

        else if (left_node_pt == bath_el_pt->node_pt(0)) {
          Sorted_bath_surface_el_pt_list.push_front(bath_el_pt);
          bath_done_el[bath_el_pt] = true;
          bath_count_done++;
          bath_el_is_inverted[bath_el_pt] = true;
        }
        // New element fits on the right of last element and is not inverted

        else if (right_node_pt == bath_el_pt->node_pt(0)) {
          Sorted_bath_surface_el_pt_list.push_back(bath_el_pt);
          bath_done_el[bath_el_pt] = true;
          bath_count_done++;
          bath_el_is_inverted[bath_el_pt] = false;
        }
        // New element fits on the right of last element and is inverted

        else if (right_node_pt == bath_el_pt->node_pt(nnod - 1)) {
          Sorted_bath_surface_el_pt_list.push_back(bath_el_pt);
          bath_done_el[bath_el_pt] = true;
          bath_count_done++;
          bath_el_is_inverted[bath_el_pt] = true;
        }

        if (bath_done_el[bath_el_pt]) {
          break;
        }
      }
    }
  }

  // Are we done?
  if (bath_count_done != (n_bath_interface_element - 1)) {
    std::ostringstream error_message;
    error_message
        << "When ordering Elements on Bath Surface, only managed to order \n"
        << bath_count_done << " of " << n_bath_interface_element
        << " face elements.\n"
        << std::endl;
    throw OomphLibError(error_message.str(), OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
  }

  // Now convert the list to a vector of points
  std::vector<BATH_INTERFACE_ELEMENT *> Sorted_bath_surface_el_pt;
  Sorted_bath_surface_el_pt.reserve(Sorted_bath_surface_el_pt_list.size());
  std::copy(std::begin(Sorted_bath_surface_el_pt_list),
            std::end(Sorted_bath_surface_el_pt_list),
            std::back_inserter(Sorted_bath_surface_el_pt));

  // Is this an anticlockwise list?
  BATH_INTERFACE_ELEMENT *first_bath_el_pt = Sorted_bath_surface_el_pt[0];
  double r_bath_first = first_bath_el_pt->interpolated_x(s_temp, 0);

  BATH_INTERFACE_ELEMENT *last_bath_el_pt =
      Sorted_bath_surface_el_pt[n_bath_interface_element - 1];
  double r_bath_last = last_bath_el_pt->interpolated_x(s_temp, 0);
  bool bath_list_anticlock = true;

  if (r_bath_first < r_bath_last) {
    bath_list_anticlock = false;
  }
  oomph_info << "Bath is anticlockwise:" << bath_list_anticlock << std::endl;

  unsigned bath_e_end = n_bath_interface_element - 1;
  unsigned bath_e_r0 = 0;

  if (bath_list_anticlock) {
    bath_e_end = 0;
    bath_e_r0 = n_bath_interface_element - 1;
  }

  // End Sort Bath Surface

  // Now want a loop over the drop surface elements, starting from the bottom.
  // During this,keep count of the bath elements that have been checked.
  // There is then a subloop over integration points
  // Start with the final bath element (at r=0, as they are listed
  // anticlockwise)
  oomph_info << "Assigning Drop Elements' external bath elements" << std::endl;
  unsigned test_e_bath = 0;
  if (bath_list_anticlock) {
    test_e_bath = n_bath_interface_element - 1;
  }
  bool film_edge_found = false;
  for (unsigned e_drop_loop = 0; e_drop_loop < n_drop_interface_element;
       e_drop_loop++) {
    unsigned e_drop = n_drop_interface_element - 1 - e_drop_loop;
    if (drop_list_anticlock) {
      e_drop = e_drop_loop;
    }
    DROP_INTERFACE_ELEMENT *drop_el_pt = Sorted_drop_surface_el_pt[e_drop];
    // Clean up previous external elements
    const unsigned interaction_index = 0;
    drop_el_pt->initialise_external_element_storage();
    unsigned n_intpt = drop_el_pt->integral_pt()->nweight();
    for (unsigned ipt = 0; ipt < n_intpt; ipt++) {
      drop_el_pt->external_element_pt(interaction_index, ipt) = 0;
    }
#ifdef FULL_ELASSIGN_OUTPUT
    oomph_info << "Drop El: " << e_drop
               << " Is Inverted: " << drop_el_is_inverted[drop_el_pt]
               << std::endl;
#endif
    // Assume we are above the bath to start.
    drop_el_pt->set_above_bath(true);
    // Loop over integration points`
    for (unsigned ipt_loop = 0; ipt_loop < n_intpt; ipt_loop++) {
      unsigned ipt = ipt_loop;
      if (drop_list_anticlock == drop_el_is_inverted[drop_el_pt]) {
        ipt = n_intpt - 1 - ipt_loop;
      }
      double s_ipt = drop_el_pt->integral_pt()->knot(ipt, 0);
      Vector<double> s_vec(1, s_ipt);
      // Coordinates of the integration point
      double r_ipt = drop_el_pt->interpolated_x(s_vec, 0);
      double z_ipt = drop_el_pt->interpolated_x(s_vec, 1);
      // outer normal
      unsigned n_dim = drop_el_pt->nodal_dimension();
      Vector<double> n_drop(n_dim);
      drop_el_pt->outer_unit_normal(ipt, n_drop);

      double n_r = n_drop[0];
      double n_z = n_drop[1];
#ifdef FULL_ELASSIGN_OUTPUT
      oomph_info << "Drop El: " << e_drop << " ipt:" << ipt << " r=" << r_ipt
                 << " z=" << z_ipt << std::endl;
      oomph_info << " n_r=" << n_r << " n_z=" << n_z << std::endl;
#endif
      // Search for the matching bath element
      bool found_bath_element = false;
      bool looped_frop_start = false;

      bool found_backup_element = false;
      double s_bath_backup = 0.0;
      unsigned bath_element_backup = 0;

      unsigned start_test_e_bath = test_e_bath;
      while (found_bath_element == false) {
        // If we've found the outer edge, then just assign the final bath
        // element
        // and set the "no film" flag
        if (film_edge_found == true) {
#ifdef FULL_ELASSIGN_OUTPUT
          oomph_info << "Assigning end of bath element " << bath_e_end
                     << ", above_bath = false " << std::endl;
#endif
          BATH_INTERFACE_ELEMENT *bath_el_pt =
              Sorted_bath_surface_el_pt[bath_e_end];
          drop_el_pt->external_element_pt(interaction_index, ipt) = bath_el_pt;
          Vector<double> s_bath_vec(1, 0.5);
          drop_el_pt->external_element_local_coord(interaction_index, ipt) =
              s_bath_vec;
          drop_el_pt->set_above_bath(false);
          break;
        } // End  if(film_edge_found==true)
#ifdef FULL_ELASSIGN_OUTPUT
        oomph_info << "Test bath_element:" << test_e_bath << std::endl;
#endif
        BATH_INTERFACE_ELEMENT *bath_el_pt =
            Sorted_bath_surface_el_pt[test_e_bath];
        Vector<double> s_0(1, 0.0);
        Vector<double> s_1(1, 0.5);
        Vector<double> s_2(1, 1.0);
        double r0 = bath_el_pt->interpolated_x(s_0, 0);
        double r2 = bath_el_pt->interpolated_x(s_2, 0);
        double z0 = bath_el_pt->interpolated_x(s_0, 1);
        double z2 = bath_el_pt->interpolated_x(s_2, 1);

        
        double b_r_alt = (r2 - r0);
        double c_r = r0 - r_ipt;
       
        double b_z_alt = (z2 - z0);
        double c_z = z0 - z_ipt;
        
        double s_bath;
        double t_bath;
        
        if (std::abs(n_r) < 0.2) {
          s_bath = -(c_r - c_z * n_r / n_z) / (b_r_alt - b_z_alt * n_r / n_z);
          t_bath = (b_z_alt * s_bath + c_z) / n_z;
        } else {
          s_bath = -(c_z - c_r * n_z / n_r) / (b_z_alt - b_r_alt * n_z / n_r);
          t_bath = (b_r_alt * s_bath + c_r) / n_r;
        }
#ifdef FULL_ELASSIGN_OUTPUT
       
        oomph_info << "s_bath=" << s_bath << " t_bath= " << t_bath << std::endl;
#endif
        

        if ((s_bath <= 1.0) && (s_bath >= 0.0) &&
            (std::abs(t_bath) <= Global_Sim_Settings::maximum_film_thickness)) {
#ifdef FULL_ELASSIGN_OUTPUT
          oomph_info << "Assigning bath element " << test_e_bath
                     << " s = " << s_bath << std::endl;
#endif
          drop_el_pt->external_element_pt(interaction_index, ipt) = bath_el_pt;
          Vector<double> s_bath_vec(1, s_bath);
          drop_el_pt->external_element_local_coord(interaction_index, ipt) =
              s_bath_vec;
          break;
        } else {

          // Check if this is nearly intersecting, and save to assign it if we
          // have to
          double backup_eps = 0.1;
          if (std::abs(t_bath) <= Global_Sim_Settings::maximum_film_thickness) {
            if ((s_bath <= (1.0 + backup_eps)) && (s_bath >= 0.0)) {
              found_backup_element = true;
              s_bath_backup = 1.0;
              bath_element_backup = test_e_bath;
            } else if ((s_bath >= (-backup_eps)) && (s_bath <= 1.0)) {
              found_backup_element = true;
              s_bath_backup = 0.0;
              bath_element_backup = test_e_bath;
            }
          }

          // If we've done a whole loop and found a backup element, assign it
          if ((test_e_bath == start_test_e_bath) &&
              (looped_frop_start == true) && (found_backup_element == true)) {
#ifdef FULL_ELASSIGN_OUTPUT
            oomph_info << "No Element Found, but found a backup element"
                       << bath_element_backup << " s=" << s_bath_backup
                       << std::endl;
#endif
            bath_el_pt = Sorted_bath_surface_el_pt[bath_element_backup];
            drop_el_pt->external_element_pt(interaction_index, ipt) =
                bath_el_pt;
            Vector<double> s_bath_vec(1, s_bath_backup);
            drop_el_pt->external_element_local_coord(interaction_index, ipt) =
                s_bath_vec;
            break;
          }

          // If we've done a whole loop and near r=0, then the normal is going
          // into negative r, so assign r=0
          if ((test_e_bath == start_test_e_bath) &&
              (looped_frop_start == true) && (r_ipt <= 0.1)) {
#ifdef FULL_ELASSIGN_OUTPUT
            oomph_info
                << "No Element Found, but near r=0 so assigning bath element"
                << bath_e_r0 << std::endl;
#endif
            bath_el_pt = Sorted_bath_surface_el_pt[bath_e_r0];
            drop_el_pt->external_element_pt(interaction_index, ipt) =
                bath_el_pt;
            Vector<double> s_bath_vec(1, 0.5);
            drop_el_pt->external_element_local_coord(interaction_index, ipt) =
                s_bath_vec;
            break;
          }
          // If we've done a whole loop and not near r=0, we've found the edge
          // of the film
          if ((test_e_bath == start_test_e_bath) &&
              (looped_frop_start == true)) {
            film_edge_found = true;
          }
          if (test_e_bath == bath_e_end) {
// Hit r=r_max, something's gone wrong, reset and note this has
// happened
#ifdef FULL_ELASSIGN_OUTPUT
            oomph_info << "Looped From Start" << std::endl;
#endif
            looped_frop_start = true;
            test_e_bath = bath_e_r0;
          } else {
            if (bath_list_anticlock) {
              test_e_bath = test_e_bath - 1;
            } else {
              test_e_bath = test_e_bath + 1;
            }
          } // End if (test_e_bath == 0)

        } // End s not in [0,1]
      }   // End while (found_bath_element == false)
    }     // End loop over int points

    // Pin the gas pressure if not under drop
    const unsigned nnodes = drop_el_pt->nnode();
    for (unsigned i = 0; i < nnodes; i++) {
      drop_el_pt->unpin_vapour_pressure(i);
      if (drop_el_pt->above_bath() == false) {
#ifdef FULL_ELASSIGN_OUTPUT
        oomph_info << "Pinning Gas Pressure" << std::endl;
#endif
        drop_el_pt->fix_vapour_pressure(
            i, Global_Physical_Variables::vapour_pressure_0);
      }
    } // End loop over nodes
  }   // End loop over elements

  // Now want a loop over the bath surface elements, starting from the end.
  // During this,keep count of the bath elements that have been checked.
  // There is then a subloop over integration points
  // Start with the first drop element (at r=0, as they are listed
  // anticlockwise)
  oomph_info << "Assigning Bath Elements' external drop elements" << std::endl;
  unsigned test_e_drop = n_drop_interface_element - 1;
  if (drop_list_anticlock) {
    test_e_drop = 0;
  }
  film_edge_found = false;
  for (unsigned e_bath_loop = 0; e_bath_loop < n_bath_interface_element;
       e_bath_loop++) {
    unsigned e_bath = e_bath_loop;
    if (bath_list_anticlock) {
      e_bath = n_bath_interface_element - 1 - e_bath_loop;
    }
    BATH_INTERFACE_ELEMENT *bath_el_pt = Sorted_bath_surface_el_pt[e_bath];
    // Clean up previous external elements
    const unsigned interaction_index = 0;
    bath_el_pt->initialise_external_element_storage();
    unsigned n_intpt = bath_el_pt->integral_pt()->nweight();
    for (unsigned ipt = 0; ipt < n_intpt; ipt++) {
      bath_el_pt->external_element_pt(interaction_index, ipt) = 0;
    }
    // Assume we are under the drop to start.
    bath_el_pt->set_under_drop(true);
#ifdef FULL_ELASSIGN_OUTPUT
    oomph_info << "Bath El: " << e_bath
               << " Is Inverted: " << bath_el_is_inverted[bath_el_pt]
               << std::endl;
#endif
    // Loop over integration points`
    for (unsigned ipt_loop = 0; ipt_loop < n_intpt; ipt_loop++) {
      unsigned ipt = n_intpt - 1 - ipt_loop;
      if (bath_list_anticlock == bath_el_is_inverted[bath_el_pt]) {
        ipt = ipt_loop;
      }
      double s_ipt = bath_el_pt->integral_pt()->knot(ipt, 0);
      Vector<double> s_vec(1, s_ipt);
      // Coordinates of the integration point
      double r_ipt = bath_el_pt->interpolated_x(s_vec, 0);
      double z_ipt = bath_el_pt->interpolated_x(s_vec, 1);
      // outer normal
      unsigned n_dim = bath_el_pt->nodal_dimension();
      Vector<double> n_bath(n_dim);
      bath_el_pt->outer_unit_normal(ipt, n_bath);

      double n_r = n_bath[0];
      double n_z = n_bath[1];
#ifdef FULL_ELASSIGN_OUTPUT
      oomph_info << "Bath El: " << e_bath << " ipt:" << ipt << " r=" << r_ipt
                 << " z=" << z_ipt << std::endl;
      oomph_info << " n_r=" << n_r << " n_z=" << n_z << std::endl;
#endif
      // Search for the matching bath element
      bool found_drop_element = false;
      bool looped_frop_start = false;

      bool found_backup_element = false;
      double s_drop_backup = 0.0;
      unsigned drop_element_backup = 0;

      unsigned start_test_e_drop = test_e_drop;
      while (found_drop_element == false) {
        // If we've found the outer edge, then just assign the final drop
        // element
        // and set the "no film" flag
        if (film_edge_found == true) {
#ifdef FULL_ELASSIGN_OUTPUT
          oomph_info << "Assigning top of drop element: " << drop_e_top
                     << " n_drop_el-1, under_drop=false" << std::endl;
#endif
          DROP_INTERFACE_ELEMENT *drop_el_pt =
              Sorted_drop_surface_el_pt[drop_e_top];
          bath_el_pt->external_element_pt(interaction_index, ipt) = drop_el_pt;
          Vector<double> s_drop_vec(1, 0.5);
          bath_el_pt->external_element_local_coord(interaction_index, ipt) =
              s_drop_vec;
          bath_el_pt->set_under_drop(false);
          break;
        } // End  if(film_edge_found==true)
#ifdef FULL_ELASSIGN_OUTPUT
        oomph_info << "Test drop_element:" << test_e_drop << std::endl;
#endif
        DROP_INTERFACE_ELEMENT *drop_el_pt =
            Sorted_drop_surface_el_pt[test_e_drop];
        Vector<double> s_0(1, 0.0);
        Vector<double> s_1(1, 0.5);
        Vector<double> s_2(1, 1.0);
        double r0 = drop_el_pt->interpolated_x(s_0, 0);
        double r2 = drop_el_pt->interpolated_x(s_2, 0);
        double z0 = drop_el_pt->interpolated_x(s_0, 1);
        double z2 = drop_el_pt->interpolated_x(s_2, 1);

      
        double b_r_alt = (r2 - r0);
        double c_r = r0 - r_ipt;
        
        double b_z_alt = (z2 - z0);
        double c_z = z0 - z_ipt;
       
        double s_drop;
        double t_drop;
       

        if (std::abs(n_r) < 0.2) {
          s_drop = -(c_r - c_z * n_r / n_z) / (b_r_alt - b_z_alt * n_r / n_z);
          t_drop = (b_z_alt * s_drop + c_z) / n_z;
        } else {
          s_drop = -(c_z - c_r * n_z / n_r) / (b_z_alt - b_r_alt * n_z / n_r);
          t_drop = (b_r_alt * s_drop + c_r) / n_r;
        }
#ifdef FULL_ELASSIGN_OUTPUT
        //          oomph_info << "abs(a)<" << a_cutoff << " so using linear
        //          approx"
        //                     << std::endl;
        oomph_info << "s_drop=" << s_drop << " t_drop= " << t_drop << std::endl;
#endif
       
        if ((s_drop <= 1.0) && (s_drop >= 0.0) &&
            (std::abs(t_drop) <= Global_Sim_Settings::maximum_film_thickness)) {
#ifdef FULL_ELASSIGN_OUTPUT
          oomph_info << "Assigning drop element " << test_e_drop
                     << " s = " << s_drop << std::endl;
#endif
          bath_el_pt->external_element_pt(interaction_index, ipt) = drop_el_pt;
          Vector<double> s_drop_vec(1, s_drop);
          bath_el_pt->external_element_local_coord(interaction_index, ipt) =
              s_drop_vec;

          break;
        } else {

          // Check if this is nearly intersecting, and save to assign it if we
          // have to
          double backup_eps = 0.1;
          if (std::abs(t_drop) <= Global_Sim_Settings::maximum_film_thickness) {
            if ((s_drop <= (1.0 + backup_eps)) && (s_drop >= 0.0)) {
              found_backup_element = true;
              s_drop_backup = 1.0;
              drop_element_backup = test_e_drop;
            } else if ((s_drop >= (-backup_eps)) && (s_drop <= 1.0)) {
              found_backup_element = true;
              s_drop_backup = 0.0;
              drop_element_backup = test_e_drop;
            }
          }

          // If we've done a whole loop and found a backup element, assign it
          if ((test_e_drop == start_test_e_drop) &&
              (looped_frop_start == true) && (found_backup_element == true)) {
#ifdef FULL_ELASSIGN_OUTPUT
            oomph_info << "No Element Found, but found a backup element"
                       << drop_element_backup << " s=" << s_drop_backup
                       << std::endl;
#endif
            drop_el_pt = Sorted_drop_surface_el_pt[drop_element_backup];
            bath_el_pt->external_element_pt(interaction_index, ipt) =
                drop_el_pt;
            Vector<double> s_drop_vec(1, s_drop_backup);
            bath_el_pt->external_element_local_coord(interaction_index, ipt) =
                s_drop_vec;
            break;
          }

          // If we've done a whole loop,and near r=0, n is pointing into
          // negative r, so assign r=0
          if ((test_e_drop == start_test_e_drop) &&
              (looped_frop_start == true) && (r_ipt <= 0.1)) {
#ifdef FULL_ELASSIGN_OUTPUT
            oomph_info
                << "No Element Found, but near r=0 so assigning drop element: "
                << drop_e_bottom << std::endl;
#endif
            drop_el_pt = Sorted_drop_surface_el_pt[drop_e_bottom];
            bath_el_pt->external_element_pt(interaction_index, ipt) =
                drop_el_pt;
            Vector<double> s_drop_vec(1, 0.5);
            bath_el_pt->external_element_local_coord(interaction_index, ipt) =
                s_drop_vec;
            break;
          }
          // If we've done a whole loop, and not near r=0 we've found the edge
          // of the film
          if ((test_e_drop == start_test_e_drop) &&
              (looped_frop_start == true)) {
            film_edge_found = true;
          }
          if (test_e_drop == drop_e_top) {
// Hit r=0, something's gone wrong, reset and note this has
// happened
#ifdef FULL_ELASSIGN_OUTPUT
            oomph_info << "Looped From Start" << std::endl;
#endif
            looped_frop_start = true;
            test_e_drop = drop_e_bottom;
          } else {
            if (drop_list_anticlock) {
              test_e_drop = test_e_drop + 1;
            } else {
              test_e_drop = test_e_drop - 1;
            }
          } // End if (test_e_drop == 0)

        } // End s not in [0,1]
      }   // End while (found_drop_element == false)
    }     // End loop over int points

  } // End loop over bath elements
  oomph_info << "Finished Assigning External Elements" << std::endl;
  // Reassign equation numbers due to the re-pinning
  oomph_info << "Number of equations: " << assign_eqn_numbers() << std::endl;
} // end of assign_external_elements

//============start_of_complete problem setup===============
/// Complete the problem setup to make the elements fully functional
//=======================================================================
template <class ELEMENT, class DROP_INTERFACE_ELEMENT,
          class BATH_INTERFACE_ELEMENT>
void ImpactProblem<ELEMENT, DROP_INTERFACE_ELEMENT,
                   BATH_INTERFACE_ELEMENT>::complete_problem_setup() {

  //-----Apply Drop parameters
  // Determine number of bulk elements in mesh
  unsigned n_element_bulk = Drop_bulk_mesh_pt->nelement();

  // Loop over the bulk elements
  for (unsigned e = 0; e < n_element_bulk; e++) {
    // Upcast from GeneralisedElement to the present element
    ELEMENT *el_pt = dynamic_cast<ELEMENT *>(Drop_bulk_mesh_pt->element_pt(e));

    // Set the Reynolds number
    el_pt->re_pt() = &Global_Physical_Variables::Re;

    // Set the Womersley number
    el_pt->re_st_pt() = &Global_Physical_Variables::ReSt;

    // Set the product of the Reynolds number and the inverse of the
    // Froude number
    el_pt->re_invfr_pt() = &Global_Physical_Variables::ReInvFr;

    // Set the direction of gravity
    el_pt->g_pt() = &Global_Physical_Variables::G;

    // Set the constitutive law
    el_pt->constitutive_law_pt() = Constitutive_law_pt;

  } // End of loop over drop bulk elements

  //-----Apply Bath parameters
  // Determine number of bulk elements in mesh
  n_element_bulk = Bath_bulk_mesh_pt->nelement();

  // Loop over the bulk elements
  for (unsigned e = 0; e < n_element_bulk; e++) {
    // Upcast from GeneralisedElement to the present element
    ELEMENT *el_pt = dynamic_cast<ELEMENT *>(Bath_bulk_mesh_pt->element_pt(e));

    // Set the Reynolds number
    el_pt->re_pt() = &Global_Physical_Variables::Re;

    // Set the Womersley number
    el_pt->re_st_pt() = &Global_Physical_Variables::ReSt;

    // Set the product of the Reynolds number and the inverse of the
    // Froude number
    el_pt->re_invfr_pt() = &Global_Physical_Variables::ReInvFr;

    // Set the direction of gravity
    el_pt->g_pt() = &Global_Physical_Variables::G;

    // Set the constitutive law
    el_pt->constitutive_law_pt() = Constitutive_law_pt;

  } // End of loop over bath bulk elements

  // Determine number of 1D interface elements in mesh
  unsigned n_interface_element = Drop_surface_mesh_pt->nelement();
  // Counters to locate the node at the top of the drop
  unsigned top_element = 0;
  unsigned top_node = 0;
  double max_z_found = -1000.0;
  unsigned bottom_element = 0;
  unsigned bottom_node = 0;
  double min_z_found = 1000.0;
  // Loop over the interface elements
  for (unsigned e = 0; e < n_interface_element; e++) {
    // Upcast from GeneralisedElement to the present element
    DROP_INTERFACE_ELEMENT *el_pt = dynamic_cast<DROP_INTERFACE_ELEMENT *>(
        Drop_surface_mesh_pt->element_pt(e));

    // Set the Strouhal number
    el_pt->st_pt() = &Global_Physical_Variables::St;

    // Set the Capillary number
    el_pt->ca_pt() = &Global_Physical_Variables::Ca;

    // Pass the Data item that contains the single external pressure value
    el_pt->set_external_pressure_data(External_pressure_data_pt);

    // Set the viscosity ratio
    el_pt->vapour_viscosity_ratio_pt() = &Global_Physical_Variables::Vis_R;

    // Set the Hamaker number
    el_pt->hamaker_pt() = &Global_Physical_Variables::A_h;

    // Set the MFP
    el_pt->mean_free_path_pt() = &Global_Physical_Variables::Mean_Free_Path;

    // Set the GKE factors
    if (Global_Physical_Variables::Use_GKE) {

      el_pt->gke_Phi_P_fct_pt() = &Global_Physical_Variables::GKE_Phi_P_fct;
      el_pt->gke_Phi_C_fct_pt() = &Global_Physical_Variables::GKE_Phi_C_fct;
      el_pt->gke_Tau_P_fct_pt() = &Global_Physical_Variables::GKE_Tau_P_fct;
      el_pt->gke_Tau_C_fct_pt() = &Global_Physical_Variables::GKE_Tau_C_fct;
    }

    // Loop over nodes
    const unsigned nnodes = el_pt->nnode();
    for (unsigned i = 0; i < nnodes; i++) {
      Node *nod_pt = el_pt->node_pt(i);
      double r = nod_pt->x(0);
      double z = nod_pt->x(1);
      // Identify if node is at the top of the drop
      if (r < 1.0e-15) {
        if (z > max_z_found) {
          // Save number of identified node
          max_z_found = z;
          top_element = e;
          top_node = i;
        }

        if (z < min_z_found) {
          // Save number of identified node
          min_z_found = z;
          bottom_element = e;
          bottom_node = i;
        }
      }
    }

  } // End of loop over interface elements

  // Pin identified Vapour pressure in identified node
  // Upcast from GeneralisedElement to the present element
  DROP_INTERFACE_ELEMENT *el_pt = dynamic_cast<DROP_INTERFACE_ELEMENT *>(
      Drop_surface_mesh_pt->element_pt(top_element));
  Node *nod_pt = el_pt->node_pt(top_node);
  double r = nod_pt->x(0);
  double z = nod_pt->x(1);
  // Fix the arc length
  el_pt->fix_arc_length(top_node, 1.0);
  // Document to check this has happened
  oomph_info << "Drop: Pinned  AL=1 at r=" << r << " and z=" << z << std::endl;

  // Pin Al at bottom in identified node
  // Upcast from GeneralisedElement to the present element
  el_pt = dynamic_cast<DROP_INTERFACE_ELEMENT *>(
      Drop_surface_mesh_pt->element_pt(bottom_element));
  nod_pt = el_pt->node_pt(bottom_node);
  r = nod_pt->x(0);
  z = nod_pt->x(1);
  // Fix the arc length
  el_pt->fix_arc_length(bottom_node, 0.0);
  // Document to check this has happened
  oomph_info << "Drop: Pinned AL=0 at r=" << r << " and z=" << z << std::endl;

  //----------Bath Surface
  unsigned right_element = 0;
  unsigned right_node = 0;
  double max_r_found = -1000.0;
  unsigned left_element = 0;
  unsigned left_node = 0;
  double min_r_found = 1000.0;
  // Determine number of 1D interface elements in mesh
  n_interface_element = Bath_surface_mesh_pt->nelement();
  oomph_info << "Bath Interface Elements:" << n_interface_element << std::endl;

  for (unsigned e = 0; e < n_interface_element; e++) {
    // Upcast from GeneralisedElement to the present element
    BATH_INTERFACE_ELEMENT *el_pt = dynamic_cast<BATH_INTERFACE_ELEMENT *>(
        Bath_surface_mesh_pt->element_pt(e));

    // Set the Strouhal number
    el_pt->st_pt() = &Global_Physical_Variables::St;

    // Set the Capillary number
    el_pt->ca_pt() = &Global_Physical_Variables::Ca;

    // Set the viscosity ratio
    el_pt->vapour_viscosity_ratio_pt() = &Global_Physical_Variables::Vis_R;

    // Pass the Data item that contains the single external pressure value
    el_pt->set_external_pressure_data(External_pressure_data_pt);

    // Set the Hamaker number
    el_pt->hamaker_pt() = &Global_Physical_Variables::A_h;

    // Set the MFP
    el_pt->mean_free_path_pt() = &Global_Physical_Variables::Mean_Free_Path;

    // Set the GKE factors
    if (Global_Physical_Variables::Use_GKE) {
      el_pt->gke_Tau_P_fct_pt() = &Global_Physical_Variables::GKE_Tau_P_fct;
      el_pt->gke_Tau_C_fct_pt() = &Global_Physical_Variables::GKE_Tau_C_fct;
    }

    // Loop over nodes
    const unsigned nnodes = el_pt->nnode();
    for (unsigned i = 0; i < nnodes; i++) {
      Node *nod_pt = el_pt->node_pt(i);
      double r = nod_pt->x(0);
      // Identify if node is at the top of the drop
      if (r > max_r_found) {
        // Save number of identified node
        max_r_found = r;
        right_element = e;
        right_node = i;
      }

      if (r < min_r_found) {
        // Save number of identified node
        min_r_found = r;
        left_element = e;
        left_node = i;
      }
    }

  } // End of loop over interface elements

  // Assign AL

  BATH_INTERFACE_ELEMENT *bath_el_pt = dynamic_cast<BATH_INTERFACE_ELEMENT *>(
      Bath_surface_mesh_pt->element_pt(left_element));
  nod_pt = bath_el_pt->node_pt(left_node);
  r = nod_pt->x(0);
  z = nod_pt->x(1);
  // Fix the arc length
  bath_el_pt->fix_arc_length(left_node, 0.0);
  // Document to check this has happened
  oomph_info << "Bath: Pinned AL=0 at r=" << r << " and z=" << z << std::endl;

  bath_el_pt = dynamic_cast<BATH_INTERFACE_ELEMENT *>(
      Bath_surface_mesh_pt->element_pt(right_element));
  nod_pt = bath_el_pt->node_pt(right_node);
  r = nod_pt->x(0);
  z = nod_pt->x(1);
  // Fix the arc length
  bath_el_pt->fix_arc_length(right_node, 1.0);
  // Document to check this has happened
  oomph_info << "Bath: Pinned AL=1 at r=" << r << " and z=" << z << std::endl;

  // --------------------------------------------
  // Set the boundary conditions for this problem
  // --------------------------------------------

  // All nodes are free by default -- just pin the ones that have
  // Dirichlet conditions here

  // ------------Drop BCs
  // Fix Symmetry boundary
  // Find number of nodes on boundary
  unsigned n_boundary_node =
      Drop_bulk_mesh_pt->nboundary_node(drop_fixed_boundary_id);
  // Loop over nodes on fixed boundary
  for (unsigned n = 0; n < n_boundary_node; n++) {
    // Pin Value of radial velocity to 0
    Drop_bulk_mesh_pt->boundary_node_pt(drop_fixed_boundary_id, n)->pin(0);

    // Set value to 0
    Drop_bulk_mesh_pt->boundary_node_pt(drop_fixed_boundary_id, n)
        ->set_value(0, 0.0);
    // Pin radial position
    Drop_bulk_mesh_pt->boundary_node_pt(drop_fixed_boundary_id, n)
        ->pin_position(0);
  }

  // Pin all azimuthal velocities
  unsigned n_node = Drop_bulk_mesh_pt->nnode();
  for (unsigned n = 0; n < n_node; n++) {
    Drop_bulk_mesh_pt->node_pt(n)->pin(2);
  }

  // ------------Bath BCs
  // Fix Symmetry boundary
  // Find number of nodes on boundary
  n_boundary_node = Bath_bulk_mesh_pt->nboundary_node(bath_sym_boundary_id);
  // Loop over nodes on fixed boundary
  for (unsigned n = 0; n < n_boundary_node; n++) {
    // Pin Value of radial velocity to 0
    Bath_bulk_mesh_pt->boundary_node_pt(bath_sym_boundary_id, n)->pin(0);

    // Set value of radial velocity to 0
    Bath_bulk_mesh_pt->boundary_node_pt(bath_sym_boundary_id, n)
        ->set_value(0, 0.0);
    // Pin radial position
    Bath_bulk_mesh_pt->boundary_node_pt(bath_sym_boundary_id, n)
        ->pin_position(0);
  }

  // Fix Base boundary
  // Find number of nodes on boundary
  n_boundary_node = Bath_bulk_mesh_pt->nboundary_node(bath_base_boundary_id);
  // Loop over nodes on fixed boundary
  for (unsigned n = 0; n < n_boundary_node; n++) {
    // Pin Value of radial velocity to 0
    Bath_bulk_mesh_pt->boundary_node_pt(bath_base_boundary_id, n)->pin(0);
    // Set value of radial velocity to 0
    Bath_bulk_mesh_pt->boundary_node_pt(bath_base_boundary_id, n)
        ->set_value(0, 0.0);
    // Pin radial position
    Bath_bulk_mesh_pt->boundary_node_pt(bath_base_boundary_id, n)
        ->pin_position(0);

    // Pin Value of vert velocity to 0
    Bath_bulk_mesh_pt->boundary_node_pt(bath_base_boundary_id, n)->pin(1);
    // Set value of vert velocity to 0
    Bath_bulk_mesh_pt->boundary_node_pt(bath_base_boundary_id, n)
        ->set_value(1, 0.0);
    // Pin vert position
    Bath_bulk_mesh_pt->boundary_node_pt(bath_base_boundary_id, n)
        ->pin_position(1);
  }

  n_boundary_node = Bath_bulk_mesh_pt->nboundary_node(bath_farbase_boundary_id);
  // Loop over nodes on fixed boundary
  for (unsigned n = 0; n < n_boundary_node; n++) {
    // Pin Value of radial velocity to 0
    Bath_bulk_mesh_pt->boundary_node_pt(bath_farbase_boundary_id, n)->pin(0);
    // Set value of radial velocity to 0
    Bath_bulk_mesh_pt->boundary_node_pt(bath_farbase_boundary_id, n)
        ->set_value(0, 0.0);
    // Pin radial position
    Bath_bulk_mesh_pt->boundary_node_pt(bath_farbase_boundary_id, n)
        ->pin_position(0);

    // Pin Value of vert velocity to 0
    Bath_bulk_mesh_pt->boundary_node_pt(bath_farbase_boundary_id, n)->pin(1);
    // Set value of vert velocity to 0
    Bath_bulk_mesh_pt->boundary_node_pt(bath_farbase_boundary_id, n)
        ->set_value(1, 0.0);
    // Pin vert position
    Bath_bulk_mesh_pt->boundary_node_pt(bath_farbase_boundary_id, n)
        ->pin_position(1);
  }

  // Fix Wall boundary
  // Find number of nodes on boundary
  n_boundary_node = Bath_bulk_mesh_pt->nboundary_node(bath_wall_boundary_id);
  // Loop over nodes on fixed boundary
  for (unsigned n = 0; n < n_boundary_node; n++) {
    // Pin Value of radial velocity to 0
    Bath_bulk_mesh_pt->boundary_node_pt(bath_wall_boundary_id, n)->pin(0);
    // Set value of radial velocity to 0
    Bath_bulk_mesh_pt->boundary_node_pt(bath_wall_boundary_id, n)
        ->set_value(0, 0.0);
    // Pin radial position
    Bath_bulk_mesh_pt->boundary_node_pt(bath_wall_boundary_id, n)
        ->pin_position(0);
    // Pin Value of vert velocity to 0
    // Bath_bulk_mesh_pt->boundary_node_pt(bath_wall_boundary_id, n)->pin(1);
    // Set value of vert velocity to 0
    // Bath_bulk_mesh_pt->boundary_node_pt(bath_wall_boundary_id, n)
    //     ->set_value(1, 0.0);
    // Pin vert position
    //  Bath_bulk_mesh_pt->boundary_node_pt(bath_wall_boundary_id, n)
    //     ->pin_position(1);
  }

  // Pin all azimuthal velocities
  n_node = Bath_bulk_mesh_pt->nnode();
  for (unsigned n = 0; n < n_node; n++) {
    Bath_bulk_mesh_pt->node_pt(n)->pin(2);
  }
} // end of complete_problem_setup

//==start_of_set_initial_condition========================================
/// \short Set initial conditions: Set all nodal velocities to zero and
/// initialise the previous velocities and nodal positions to correspond
/// to an impulsive start
//========================================================================
template <class ELEMENT, class DROP_INTERFACE_ELEMENT,
          class BATH_INTERFACE_ELEMENT>
void ImpactProblem<ELEMENT, DROP_INTERFACE_ELEMENT,
                   BATH_INTERFACE_ELEMENT>::set_initial_condition() {
  if (CommandLineArgs::command_line_flag_has_been_set("--restart_file")) {
    // Pointer to restart file
    ifstream *restart_file_pt = 0;

    // Open restart file
    restart_file_pt =
        new ifstream(Global_Sim_Settings::Restart_file, ios_base::in);
    if (restart_file_pt != 0) {
      cout << "Have opened " << Global_Sim_Settings::Restart_file
           << " for restart. " << std::endl;
    } else {
      std::ostringstream error_stream;
      error_stream << "ERROR while trying to open "
                   << Global_Sim_Settings::Restart_file << " for restart."
                   << std::endl;

      throw OomphLibError(error_stream.str(), OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
    // Read the problem data from the restart file
    restart(*restart_file_pt);

  } else {
    // Initialise timestep -- also sets the weights for all timesteppers
    // in the problem.
    initialise_dt(Global_Sim_Settings::max_dt);
    Next_dt = Global_Sim_Settings::max_dt;

    // Determine number of nodes in mesh
    unsigned n_node = Drop_bulk_mesh_pt->nnode();

    // Loop over all nodes in mesh
    for (unsigned n = 0; n < n_node; n++) {
      // Loop over the velocity components
      // Assign radial and angular to 0, set vertical to Initial velocity
      Drop_bulk_mesh_pt->node_pt(n)->set_value(0, 0.0);
      Drop_bulk_mesh_pt->node_pt(n)->set_value(
          1, -1.0 * Global_Physical_Variables::Initial_speed);
      Drop_bulk_mesh_pt->node_pt(n)->set_value(2, 0.0);
    }

    n_node = Bath_bulk_mesh_pt->nnode();

    // Loop over all nodes in mesh
    for (unsigned n = 0; n < n_node; n++) {
      // Loop over the velocity components
      // Assign all to 0
      Bath_bulk_mesh_pt->node_pt(n)->set_value(0, 0.0);
      Bath_bulk_mesh_pt->node_pt(n)->set_value(1, 0.0);
      Bath_bulk_mesh_pt->node_pt(n)->set_value(2, 0.0);
    }

    // Loop over drop surface mesh

    const unsigned n_interface_element = Drop_surface_mesh_pt->nelement();
    for (unsigned e = 0; e < n_interface_element; e++) {
      DROP_INTERFACE_ELEMENT *el_pt = dynamic_cast<DROP_INTERFACE_ELEMENT *>(
          Drop_surface_mesh_pt->element_pt(e));

      const unsigned nnodes = el_pt->nnode();
      for (unsigned i = 0; i < nnodes; i++) {

        el_pt->set_vapour_pressure(
            i, Global_Physical_Variables::vapour_pressure_0);
      }
    }
    // Initialise the previous velocity values and nodal positions
    // for timestepping corresponding to an impulsive start
    assign_initial_values_impulsive();
  } // End if no restart
} // End of set_initial_condition

// Compute min height
template <class ELEMENT, class DROP_INTERFACE_ELEMENT,
          class BATH_INTERFACE_ELEMENT>
double ImpactProblem<ELEMENT, DROP_INTERFACE_ELEMENT,
                     BATH_INTERFACE_ELEMENT>::min_height() {
  double min_value = 1000.0;
  // Loop over surface elemnets
  const unsigned n_interface_element = Drop_surface_mesh_pt->nelement();
  for (unsigned e = 0; e < n_interface_element; e++) {
    DROP_INTERFACE_ELEMENT *el_pt = dynamic_cast<DROP_INTERFACE_ELEMENT *>(
        Drop_surface_mesh_pt->element_pt(e));
    // Only check the elements bounding the film
    if (el_pt->above_bath() == false) {
      continue;
    }
    // Check if height is less than lowest previously found
    double z = el_pt->min_thickness();
    if (z < min_value) {
      min_value = z;
    }
  }

  return min_value;
} // End of min_height

//==start_of_doc_solution=================================================
/// Document the solution
//========================================================================
template <class ELEMENT, class DROP_INTERFACE_ELEMENT,
          class BATH_INTERFACE_ELEMENT>
void ImpactProblem<ELEMENT, DROP_INTERFACE_ELEMENT,
                   BATH_INTERFACE_ELEMENT>::doc_solution(bool full_output,
                                                         bool post_adapt) {
  // Output the time
  oomph_info << "Outputting at time " << time_pt()->time() << std::endl;

  ofstream some_file;
  char filename_surface[100];
  char filename_bulk[100];
  char filename_bath_bulk[100];
  char filename_bath_surface[100];
  char filename_size[100];
  if (post_adapt) {
    sprintf(filename_bulk, "%s/elements_%i_post_adapt.dat",
            Doc_info.directory().c_str(), Doc_info.number());

    sprintf(filename_size, "%s/elements_size_%i_post_adapt.dat",
            Doc_info.directory().c_str(), Doc_info.number());

    sprintf(filename_bath_bulk, "%s/elements_bath_%i_post_adapt.dat",
            Doc_info.directory().c_str(), Doc_info.number());

    sprintf(filename_surface, "%s/surface_elements_%i_post_adapt.dat",
            Doc_info.directory().c_str(), Doc_info.number());

    sprintf(filename_bath_surface, "%s/surface_elements_bath_%i_post_adapt.dat",
            Doc_info.directory().c_str(), Doc_info.number());
  } else {
    sprintf(filename_bulk, "%s/elements_%i.dat", Doc_info.directory().c_str(),
            Doc_info.number());

    sprintf(filename_bath_bulk, "%s/elements_bath_%i.dat",
            Doc_info.directory().c_str(), Doc_info.number());

    sprintf(filename_size, "%s/elements_size_%i.dat",
            Doc_info.directory().c_str(), Doc_info.number());

    sprintf(filename_surface, "%s/surface_elements_%i.dat",
            Doc_info.directory().c_str(), Doc_info.number());

    sprintf(filename_bath_surface, "%s/surface_elements_bath_%i.dat",
            Doc_info.directory().c_str(), Doc_info.number());
  }
  if (full_output) {
    // Output bulk elements
    some_file.open(filename_bulk);
    Drop_bulk_mesh_pt->output(some_file, 3);
    some_file.close();

    // Output bath bulk elements
    some_file.open(filename_bath_bulk);
    Bath_bulk_mesh_pt->output(some_file, 3);
    some_file.close();

    // // Output element areas
    // some_file.open(filename_size);
    // int n_elements = this->Drop_bulk_mesh_pt->nelement();
    // for (int j = 0; j < n_elements; j++) {
    //   ELEMENT *el_pt =
    //       dynamic_cast<ELEMENT *>(this->Drop_bulk_mesh_pt->element_pt(j));
    //   some_file << j << " " << el_pt->size() << std::endl;
    // }
    // some_file.close();
  }

  // Save the drop surface boundary
  some_file.open(filename_surface);
  Drop_surface_mesh_pt->output(some_file, 3);
  some_file.close();

  // Save the bath surface boundary
  some_file.open(filename_bath_surface);
  Bath_surface_mesh_pt->output(some_file, 3);
  some_file.close();

  if (post_adapt == false) {

    Trace_file << doc_info().number() << " " << time_pt()->time() << " "
               << time_pt()->dt() << " " << global_temporal_error_norm() << " "
               << min_height() << " " << mesh_pt()->nelement() << std::endl;
  }
} // End of doc_solution

//==start_of_dump_it=================================================
/// Dump contents
//========================================================================
template <class ELEMENT, class DROP_INTERFACE_ELEMENT,
          class BATH_INTERFACE_ELEMENT>
void ImpactProblem<ELEMENT, DROP_INTERFACE_ELEMENT,
                   BATH_INTERFACE_ELEMENT>::dump_it(ofstream &dump_file) {

  // Write step number
  dump_file << Doc_info.number() << " # step number" << std::endl;

  // Write suggested next timestep
  dump_file << Next_dt << " # suggested next timestep" << std::endl;

  dump_file << Adaptcount << " # counter for adaption" << std::endl;
  dump_file << Full_out_count << " # counter for full output" << std::endl;
  dump_file << Restart_count << " # counter for restarts" << std::endl;

  // Dump the problem
  Problem::dump(dump_file);

} // end of dump_it

//==start_of_restart=================================================
/// Restart
//========================================================================
template <class ELEMENT, class DROP_INTERFACE_ELEMENT,
          class BATH_INTERFACE_ELEMENT>
void ImpactProblem<ELEMENT, DROP_INTERFACE_ELEMENT,
                   BATH_INTERFACE_ELEMENT>::restart(ifstream &restart_file) {

  if (restart_file.is_open()) {
    oomph_info << "restart file exists, start reading" << std::endl;

    // Read line up to termination sign
    string input_string;
    getline(restart_file, input_string, '#');

    // Ignore rest of line
    restart_file.ignore(80, '\n');

    // Read in step number
    Doc_info.number() = unsigned(atof(input_string.c_str()));

    // Read line up to termination sign
    getline(restart_file, input_string, '#');

    // Ignore rest of line
    restart_file.ignore(80, '\n');

    // Read suggested next timestep
    Next_dt = double(atof(input_string.c_str()));

    // Read line up to termination sign
    getline(restart_file, input_string, '#');

    // Ignore rest of line
    restart_file.ignore(80, '\n');

    // Read adapt count
    Adaptcount = unsigned(atof(input_string.c_str()));

    // Read line up to termination sign
    getline(restart_file, input_string, '#');

    // Ignore rest of line
    restart_file.ignore(80, '\n');

    // Read Full_out_count
    Full_out_count = unsigned(atof(input_string.c_str()));

    // Read line up to termination sign
    getline(restart_file, input_string, '#');

    // Ignore rest of line
    restart_file.ignore(80, '\n');

    // Read Restart_count
    Restart_count = unsigned(atof(input_string.c_str()));

    Problem::read(restart_file);
    oomph_info << "End Restart" << std::endl;

  } else {
    oomph_info << "restart file does not exist" << std::endl;
  }

} // End of restart

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

//==start_of_main=========================================================
/// Driver code for axisym drop impact with vapour layer
//========================================================================
int main(int argc, char *argv[]) {
#ifdef OOMPH_HAS_MPI

  // Initialise MPI
  MPI_Helpers::init(argc, argv);

  oomph_mpi_output.restrict_output_to_single_processor();
#endif

  // Store command line arguments
  CommandLineArgs::setup(argc, argv);

  CommandLineArgs::specify_command_line_flag(
      "--folder", &Global_Sim_Settings::result_folder, "Result Folder");

  // Name of restart file
  CommandLineArgs::specify_command_line_flag(
      "--restart_file", &Global_Sim_Settings::Restart_file);

  CommandLineArgs::specify_command_line_flag(
      "--reynolds", &Global_Physical_Variables::Re, "Reynolds Number");

  CommandLineArgs::specify_command_line_flag(
      "--gravity", &Global_Physical_Variables::ReInvFr, "Gravity Number");

  CommandLineArgs::specify_command_line_flag(
      "--capilliary", &Global_Physical_Variables::Ca, "Capilliary Number");

  CommandLineArgs::specify_command_line_flag(
      "--viscosity", &Global_Physical_Variables::Vis_R, "Viscosity Ratio");

  CommandLineArgs::specify_command_line_flag(
      "--hamaker", &Global_Physical_Variables::A_h, "Dim-less Hamaker number");

  CommandLineArgs::specify_command_line_flag(
      "--mfp", &Global_Physical_Variables::Mean_Free_Path,
      "Dim-less Mean-free-path");

  CommandLineArgs::specify_command_line_flag("--usegke", "Use GKE functions");

  CommandLineArgs::specify_command_line_flag(
      "--tmax", &Global_Physical_Variables::t_max, "Max time");

  CommandLineArgs::specify_command_line_flag(
      "--initheight", &Global_Physical_Variables::Drop_start_height,
      "Initial Drop Centre Height");

  CommandLineArgs::specify_command_line_flag(
      "--initradius", &Global_Physical_Variables::Drop_radius,
      "Initial Drop radius");

  CommandLineArgs::specify_command_line_flag(
      "--initspeed", &Global_Physical_Variables::Initial_speed,
      "Initial Drop Speed");

  CommandLineArgs::specify_command_line_flag(
      "--bathdepth", &Global_Physical_Variables::Bath_depth, "Bath Depth");

  CommandLineArgs::specify_command_line_flag(
      "--bathradius", &Global_Physical_Variables::Bath_radius, "Bath Radius");

  CommandLineArgs::specify_command_line_flag(
      "--maxdt", &Global_Sim_Settings::max_dt, "Max Timestep");

  CommandLineArgs::parse_and_assign();
  CommandLineArgs::doc_specified_flags();

  // Finalise parameters
  // Compute the Womersley number
  Global_Physical_Variables::ReSt =
      Global_Physical_Variables::Re * Global_Physical_Variables::St;
  // Set direction of gravity (vertically downwards)
  Global_Physical_Variables::G[0] = 0.0;
  Global_Physical_Variables::G[1] = -1.0;
  Global_Physical_Variables::G[2] = 0.0;

  if (CommandLineArgs::command_line_flag_has_been_set("--usegke")) {
    Global_Physical_Variables::Use_GKE = true;
  }

  oomph_info << "Driver: Constructing Problem" << std::endl;
  // Create problem with Axisym Taylor hood element, with pseudosolid mesh
  // update and projectability for remeshing
  ImpactProblem<
      // Main NS Element
      ProjectableAxisymmetricTaylorHoodElement<PseudoSolidNodeUpdateElement<
          CleanOutputAxisymmetricTTaylorHoodElement, TPVDElement<2, 3>>>,
      // Drop Surface Lubrication Element
      UpperInterfaceCoupledLubricationElement<
          // Base Lubrication Element
          ElasticAxisymmetricUpperInterfaceLubricationElement<
              ProjectableAxisymmetricTaylorHoodElement<
                  PseudoSolidNodeUpdateElement<
                      CleanOutputAxisymmetricTTaylorHoodElement,
                      TPVDElement<2, 3>>>>,
          // External Element (Bath Surface Element)
          ElasticAxisymmetricLowerInterfaceLubricationElement<
              ProjectableAxisymmetricTaylorHoodElement<
                  PseudoSolidNodeUpdateElement<
                      CleanOutputAxisymmetricTTaylorHoodElement,
                      TPVDElement<2, 3>>>>>,
      // Bath Surface element
      LowerInterfaceCoupledLubricationElement<
          // Base Bath Element
          ElasticAxisymmetricLowerInterfaceLubricationElement<
              ProjectableAxisymmetricTaylorHoodElement<
                  PseudoSolidNodeUpdateElement<
                      CleanOutputAxisymmetricTTaylorHoodElement,
                      TPVDElement<2, 3>>>>,
          // External Element (Drop Surface Lubrication Element)
          ElasticAxisymmetricUpperInterfaceLubricationElement<
              ProjectableAxisymmetricTaylorHoodElement<
                  PseudoSolidNodeUpdateElement<
                      CleanOutputAxisymmetricTTaylorHoodElement,
                      TPVDElement<2, 3>>>>>>
      problem;

  // Set initial conditions
  oomph_info << "Driver: Setting Initial Conditions" << std::endl;
  problem.set_initial_condition();
  problem.doc_info().number()++;
  oomph_info << "Driver:Writing trace header" << std::endl;

  problem.write_trace_file_header();

  // Initial timestep: Use the one used when setting up the initial
  // condition
  double dt = problem.next_dt();

  // Run the unsteady simulation

  // Counters for mesh adaption and output
  bool before_impact = true;
  bool before_end_impact = true;

  // Loop until final time
  while (problem.time_pt()->time() < Global_Physical_Variables::t_max) {
    oomph_info << "Driver: Computing Step " << problem.doc_info().number()
               << " at time " << problem.time_pt()->time() << std::endl;

    if (before_impact) {
      if (problem.time_pt()->time() > Global_Sim_Settings::impact_start_time) {
        before_impact = false;
        problem.maximum_dt() = Global_Sim_Settings::max_dt_impact;
      }
    }
    if (before_end_impact) {
      if (problem.time_pt()->time() > Global_Sim_Settings::impact_end_time) {
        before_end_impact = false;
        problem.maximum_dt() = Global_Sim_Settings::max_dt;
      }
    }
    // Adapt periodically
    if (problem.adaptcount() == Global_Sim_Settings::steps_per_adapt) {
      problem.adaptcount() = 1;
      oomph_info << "Driver: Adapting:" << std::endl;
      problem.adapt();
      oomph_info << "Driver: Finished Adaption" << std::endl;
      // Uncomment to output full solution after adaption (for debugging)
      //      problem.doc_solution(doc_info, true, true);
    } else {
      problem.adaptcount() = problem.adaptcount() + 1;
    }
    problem.assign_external_elements();
    // Take one adaptive timestep
    oomph_info << "Taking Step" << std::endl;
    double dt_next = problem.adaptive_unsteady_newton_solve(
        dt, Global_Sim_Settings::epsilon_t);

    oomph_info << "Driver: Finished Step" << std::endl;

    // Doc solution
    oomph_info << "Driver: Docing Solution" << std::endl;
    if (problem.full_out_count() ==
        Global_Sim_Settings::steps_per_full_output) {
      problem.doc_solution(true, false);
      problem.full_out_count() = 0;
    } else {
      problem.doc_solution(false, false);
    }
    problem.full_out_count() = problem.full_out_count() + 1;
    dt = dt_next;
    problem.next_dt() = dt;

    if (problem.restart_count() == Global_Sim_Settings::steps_per_restart) {
      problem.restart_count() = 1;
      // Write restart file
      ofstream some_file;
      char filename[100];
      sprintf(filename, "%s/restart%i.dat",
              problem.doc_info().directory().c_str(),
              problem.doc_info().number());
      some_file.open(filename);
      some_file.precision(20);
      problem.dump_it(some_file);
      some_file.close();
    } else {
      problem.restart_count() = problem.restart_count() + 1;
    }
    // Increment counter for solutions
    problem.doc_info().number()++;

    double min_height = problem.min_height();
    // Exit loop if contact found
    if (min_height < Global_Sim_Settings::Min_height_cutoff) {
      oomph_info << "Driver: Contact detected, ending simulation. Min height="
                 << min_height << std::endl;
      break;
    }
    // Exit if bounce detected
    if (min_height > 1.1 * (Global_Physical_Variables::Drop_start_height -
                            Global_Physical_Variables::Drop_radius)) {
      oomph_info << "Driver: Rebound detected, ending simulation. Min height="
                 << min_height << std::endl;
      break;
    }
    oomph_info << "Driver: Finished Timestep Loop" << std::endl;
  } // End of timestepping loop

  // Finalise MPI
#ifdef OOMPH_HAS_MPI

  MPI_Helpers::finalize();

#endif
} // End of main
