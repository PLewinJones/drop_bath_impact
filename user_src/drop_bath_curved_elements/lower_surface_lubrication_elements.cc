// Interface lubrication element. Unified problems
#include "lower_surface_lubrication_elements.h"

namespace oomph {
// Define the default physical value to be one
double LowerInterfaceLubricationElement::Default_Physical_Constant_Value = 1.0;

//=====================================================================
/// Get the Arc Length
//=====================================================================
double
LowerInterfaceLubricationElement::interpolated_Al(const Vector<double> &s) {
  // Find number of nodes
  unsigned n_node = this->nnode();

  // Local shape function
  Shape psi(n_node);

  // Find values of shape function
  this->shape(s, psi);

  // Initialise value of Pv
  double Al = 0.0;

  // Loop over the local nodes and sum
  for (unsigned l = 0; l < n_node; l++) {
    Al += this->nodal_value(l, this->Al_index[l]) * psi(l);
  }

  return (Al);
}
//=======================================================================
/// Overload the Helper function to calculate the residuals and
/// jacobian entries. This particular function ensures that the
/// additional entries are calculated inside the integration loop
void LowerInterfaceLubricationElement::
    add_additional_residual_contributions_interface(
        Vector<double> &residuals, DenseMatrix<double> &jacobian,
        const unsigned &flag, const Shape &psif, const DShape &dpsifds,
        const DShape &dpsifdS, const DShape &dpsifdS_div,
        const Vector<double> &s, const Vector<double> &interpolated_x,
        const Vector<double> &interpolated_n, const double &W,
        const double &J) {
  // Find out how many nodes there are
  unsigned n_node = this->nnode();

  // Storage for the local equation numbers and unknowns
  int local_eqn = 0, local_unknown = 0;

  // Surface Thin-film equation

  // Find the index at which the velocity is stored
  Vector<unsigned> u_index = this->U_index_interface;

  // Read out the physical constants
  const double Vis_r = this->vapour_viscosity_ratio();
  const double Ha = this->hamaker();
  const double Mfp = this->mean_free_path();

  // Very Hacky way to get the intergration point. Should rewrite the surface
  // element to avoid this, but for my purposes this works. Will need to do this
  // properly if ever do 2D or 3D impact.
  unsigned ipt_out;
  if (s[0] < 0.4) {
    ipt_out = 0;
  } else if (s[0] > 0.6) {
    ipt_out = 2;
  } else {
    ipt_out = 1;
  };
  // oomph_info << "IN LUBRICATION ELEMENT: ipt_out = "<< ipt_out<< " , s=
  // "<<s[0]<< std::endl;
  // Get base height from external element if exists
  const unsigned n_dim = this->node_pt(0)->ndim();

  Vector<double> x_upper(n_dim, 0.0);
  Vector<double> u_upper(n_dim, 0.0);
  double pv_upper = 0.0;
  Vector<double> grad_pv_upper(n_dim, 0.0);

  get_upper_surface(ipt_out, x_upper, u_upper, pv_upper, grad_pv_upper);

  double h_film = (x_upper[0] - interpolated_x[0])*interpolated_n[0]
  + (x_upper[1] - interpolated_x[1])*interpolated_n[1];

  
  // Get the GKE factors

  double GKE_Tau_P = 1.0;
  double GKE_Tau_C = 1.0;
  double Kn = Mfp / h_film;

  get_GKE_Tau_P(Kn, GKE_Tau_P);
  get_GKE_Tau_C(Kn, GKE_Tau_C);

  // Now calculate the vapour pressure at this point
  // Assuming the same shape functions are used (which they are)

  // The mesh and fluid velocity, and grad of vapour pressure
  Vector<double> interpolated_u(n_dim, 0.0);
  Vector<double> interpolated_grad_Al(n_dim, 0.0);

  // Loop over the shape functions
  for (unsigned l = 0; l < n_node; l++) {
    const double psi = psif(l);
    const double Al_ = this->nodal_value(l, this->Al_index[l]);
    // Velocity and Mesh Velocity
    for (unsigned j = 0; j < n_dim; j++) {
      const double u_ = this->nodal_value(l, u_index[j]);
      interpolated_u[j] += u_ * psi;
      interpolated_grad_Al[j] += Al_ * dpsifdS(l, j);
    }
  }

  Vector<double> interpolated_t(n_dim, 0.0);
  interpolated_t[0] = -interpolated_n[1]; // See elements.cc line 5766
  interpolated_t[1] = interpolated_n[0];

  // Now we add the various terms to the appropriate residuals and jacobians
  // The main loop is over the nodes
  for (unsigned l = 0; l < n_node; l++) {
    if (under_drop_flag == true) {
      // Loop over velocity components
      for (unsigned i = 0; i < n_dim; i++) {
        // Read out the apprporiate local equation for the u residual and
        // Jacobian terms
        local_eqn = this->nodal_local_eqn(l, u_index[i]);
        // If not a boundary condition
        if (local_eqn >= 0) {

          // Normal component
          residuals[local_eqn] +=
              (-Vis_r * pv_upper + Ha / (h_film * h_film * h_film)) *
              interpolated_n[i] * psif(l) * W * J;

          double tangent_term = 0.0;
          for (unsigned i2 = 0; i2 < n_dim; i2++) {
            tangent_term +=
               (u_upper[i2] - interpolated_u[i2]) *
                interpolated_t[i2];
          }

          // and the tangential term:
          residuals[local_eqn] +=
              -Vis_r *((h_film / (2.0 * GKE_Tau_P)) * grad_pv_upper[i]-(tangent_term) * interpolated_t[i]/ (h_film * GKE_Tau_C)) * psif(l) * W * J;

          // We also need to worry about the jacobian terms
          if (flag) {
            // Loop over the nodes again, this time for the unknowns
            for (unsigned l2 = 0; l2 < n_node; l2++) {
              // Loop over the velocity components
              for (unsigned i2 = 0; i2 < n_dim; i2++) {
                // Get the unknown u index
                local_unknown = this->nodal_local_eqn(l2, u_index[i2]);

                // If not a boundary condition
                if (local_unknown >= 0) {
                 jacobian(local_eqn, local_unknown) +=
                      -Vis_r * (psif(l2) / (h_film * GKE_Tau_C))*interpolated_t[i2]*interpolated_t[i] * psif(l) * W *J;
                }
              }
            }
          }
        }
      }
    }

    // Read out the apprporiate local equation for the arc length residual and
    // Jacobian terms
    local_eqn = this->nodal_local_eqn(l, this->Al_index[l]);

    if (local_eqn >= 0) {
      for (unsigned i = 0; i < n_dim; i++) {
        residuals[local_eqn] += interpolated_grad_Al[i] * dpsifdS(l, i) * W * J;
      }
      if (flag) {
        // Loop over the nodes again, this time for the unknowns
        for (unsigned l2 = 0; l2 < n_node; l2++) {
          // Get the unknown Al_index
          local_unknown = this->nodal_local_eqn(l2, this->Al_index[l2]);
          if (local_unknown >= 0) {
            for (unsigned i = 0; i < n_dim; i++) {
              jacobian(local_eqn, local_unknown) +=
                  dpsifdS(l2, i) * dpsifdS(l, i) * W * J;
            }
          }
        }
      }
    }
  }
}

//=======================================================
/// Overload the output function
//=======================================================
void LowerInterfaceLubricationElement::output(std::ostream &outfile,
                                              const unsigned &n_plot) {
  outfile.precision(16);

  const unsigned el_dim = this->dim();
  const unsigned n_dim = this->nodal_dimension();
  const unsigned n_velocity = this->U_index_interface.size();

  // Set output Vector
  Vector<double> s(el_dim);
  Vector<double> n(n_dim);
  Vector<double> u(n_velocity);

  // Loop over plot points
  unsigned num_plot_points = this->nplot_points(n_plot);
  for (unsigned iplot = 0; iplot < num_plot_points; iplot++) {
    // Get local coordinates of plot point
    this->get_s_plot(iplot, n_plot, s);
    // Get the outer unit normal
    this->outer_unit_normal(s, n);

    // double u_n = 0.0;
    for (unsigned i = 0; i < n_velocity; i++) {
      u[i] = this->interpolated_u(s, i);
    }
    // Output the x,y,u,v
    for (unsigned i = 0; i < n_dim; i++)
      outfile << this->interpolated_x(s, i) << " ";
    for (unsigned i = 0; i < n_dim; i++)
      outfile << u[i] << " ";

    // Output the vapour pressure
    outfile << interpolated_Al(s) << " ";

    outfile << std::endl;
  }
}
} // namespace oomph
