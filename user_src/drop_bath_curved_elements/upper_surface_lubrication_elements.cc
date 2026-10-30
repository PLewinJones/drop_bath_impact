// Interface lubrication element. Unified problems
#include "upper_surface_lubrication_elements.h"

namespace oomph {
// Define the default physical value to be one
double UpperInterfaceLubricationElement::Default_Physical_Constant_Value = 1.0;

//=====================================================================
/// Get the Vapour pressure
//=====================================================================
double UpperInterfaceLubricationElement::interpolated_Pv(
    const Vector<double> &s) const {
  // Find number of nodes
  unsigned n_node = this->nnode();

  // Local shape function
  Shape psi(n_node);

  // Find values of shape function
  this->shape(s, psi);

  // Initialise value of Pv
  double Pv = 0.0;

  // Loop over the local nodes and sum
  for (unsigned l = 0; l < n_node; l++) {
    Pv += this->nodal_value(l, this->Pv_index[l]) * psi(l);
  }

  return (Pv);
}

//=====================================================================
/// Get the Arc Length
//=====================================================================
double
UpperInterfaceLubricationElement::interpolated_Al(const Vector<double> &s) {
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
void UpperInterfaceLubricationElement::
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
  Vector<double> x_lower(n_dim, 0.0);
  Vector<double> n_base(n_dim, 0.0);
  Vector<double> u_base(n_dim, 0.0);
   Vector<double> dx_dt_base(n_dim, 0.0);

  get_lower_surface(ipt_out, x_lower, n_base, u_base,dx_dt_base);

  double h_film = (x_lower[0] - interpolated_x[0]) * interpolated_n[0] +
                  (x_lower[1] - interpolated_x[1]) * interpolated_n[1];
  // Get the GKE factors

  double GKE_Tau_P = 1.0;
  double GKE_Tau_C = 1.0;
  double GKE_Phi_P = 1.0;
  double GKE_Phi_C = 1.0;

  // Now calculate the vapour pressure at this point
  // Assuming the same shape functions are used (which they are)
  double interpolated_Pv = 0.0;

  // The mesh and fluid velocity, and grad of vapour pressure
  Vector<double> interpolated_u(n_dim, 0.0);
  Vector<double> interpolated_grad_Pv(n_dim, 0.0);
  Vector<double> interpolated_grad_Al(n_dim, 0.0);
  Vector<double> interpolated_w(n_dim, 0.0);
  double interpolated_div_w = 0.0;

  // Loop over the shape functions
  for (unsigned l = 0; l < n_node; l++) {
    const double psi = psif(l);
    const double Pv_ = this->nodal_value(l, this->Pv_index[l]);
    const double Al_ = this->nodal_value(l, this->Al_index[l]);
    interpolated_Pv += Pv_ * psi;
    // Velocity and Mesh Velocity
    for (unsigned j = 0; j < n_dim; j++) {
      const double u_ = this->nodal_value(l, u_index[j]);
      const double w_ = this->dnodal_position_dt(l, j);
      interpolated_u[j] += u_ * psi;
      interpolated_w[j] += w_ * psi;
      interpolated_div_w += w_ * dpsifdS_div(l, j);
      interpolated_grad_Pv[j] += Pv_ * dpsifdS(l, j);
      interpolated_grad_Al[j] += Al_ * dpsifdS(l, j);
    }
  }
     
  // oomph_info << "IN drop element: r= " << interpolated_x[0]
  //            << " , z= " << interpolated_x[1] << std::endl;

  Vector<double> interpolated_t(n_dim, 0.0);
  interpolated_t[0] = -interpolated_n[1]; // See elements.cc line 5766
  interpolated_t[1] = interpolated_n[0];

  Vector<double> t_base(n_dim, 0.0);
  t_base[0] = -n_base[1]; // See elements.cc line 5766
  t_base[1] = n_base[0];

  double td_dot_db = 0.0;
  double tb_dot_nd = 0.0;
  for (unsigned i = 0; i < n_dim; i++) {
    td_dot_db += interpolated_t[i] * t_base[i];
    tb_dot_nd += t_base[i] * interpolated_n[i];
  }

  double Kn = Mfp / h_film;

  double dhdt = 0.0;
  for (unsigned i = 0; i < n_dim; i++) {
    dhdt +=
        (dx_dt_base[i] - interpolated_w[i]) * interpolated_n[i] -
        (tb_dot_nd / td_dot_db) *
            ((dx_dt_base[i] - interpolated_w[i]) * interpolated_t[i] );
  }

  get_GKE_Tau_P(Kn, GKE_Tau_P);
  get_GKE_Tau_C(Kn, GKE_Tau_C);
  get_GKE_Phi_P(Kn, GKE_Phi_P);
  get_GKE_Phi_C(Kn, GKE_Phi_C);
  // Now we add the various terms to the appropriate residuals and jacobians
  // The main loop is over the nodes
  for (unsigned l = 0; l < n_node; l++) {
    // Read out the apprporiate local equation for the Pv residual and
    // Jacobian terms
    local_eqn = this->nodal_local_eqn(l, this->Pv_index[l]);

    // If not a boundary condition
    if (local_eqn >= 0) {
      // Add Residual terms
      // First part
      residuals[local_eqn] +=
          (dhdt + h_film * interpolated_div_w) * psif(l) * W * J;
      // Second Part
      double couette_term = 0.0;
      for (unsigned i = 0; i < n_dim; i++) {
        couette_term +=
            GKE_Phi_C *
            ((interpolated_u[i] + u_base[i] - 2 * interpolated_w[i]) *
                 interpolated_t[i] ) /
            2.0;
      }
      double velocity_term = 0.0;
      for (unsigned i = 0; i < n_dim; i++) {
        velocity_term +=
            h_film *
            (interpolated_grad_Pv[i] * GKE_Phi_P * h_film * h_film / 12.0 -
             couette_term*interpolated_t[i]) *
            dpsifdS(l, i);
      }
      residuals[local_eqn] += velocity_term * W * J;

      // We also need to compute the jacobian terms
      if (flag) {
        // Loop over the nodes again, this time for the unknowns
        for (unsigned l2 = 0; l2 < n_node; l2++) {
          // Get the unknown pv_index
          local_unknown = this->nodal_local_eqn(l2, this->Pv_index[l2]);
          // If not a boundary condition
          if (local_unknown >= 0) {
            for (unsigned i = 0; i < n_dim; i++) {
              jacobian(local_eqn, local_unknown) +=
                  (GKE_Phi_P * h_film * h_film * h_film / 12.0) *
                  (dpsifdS(l2, i)) * (dpsifdS(l, i)) * W * J;
            }
          }

          // Loop over the velocity components
          for (unsigned i2 = 0; i2 < n_dim; i2++) {
            // Get the unknown u index
            local_unknown = this->nodal_local_eqn(l2, u_index[i2]);

            // If not a boundary condition
            if (local_unknown >= 0) {
              for (unsigned i = 0; i < n_dim; i++) {
              jacobian(local_eqn, local_unknown) +=
                  -(GKE_Phi_C * h_film / 2.0) * psif(l2)*interpolated_t[i2]*interpolated_t[i] * dpsifdS(l, i) *
                  W * J;
              }
            }
          }
        }
      }
    }

    if (above_bath() == true) {
      // Loop over velocity components
      for (unsigned i = 0; i < n_dim; i++) {
        // Read out the apprporiate local equation for the u residual and
        // Jacobian terms
        local_eqn = this->nodal_local_eqn(l, u_index[i]);
        // If not a boundary condition
        if (local_eqn >= 0) {
          // Add Residual terms
          // Normal component
          residuals[local_eqn] +=
              (-Vis_r * interpolated_Pv + Ha / (h_film * h_film * h_film)) *
              interpolated_n[i] * psif(l) * W * J;

          double tangent_term = 0.0;
          for (unsigned i2 = 0; i2 < n_dim; i2++) {
            tangent_term +=
            (interpolated_u[i2] - u_base[i2] ) *interpolated_t[i2];
          }

          // and the tangential term:
          residuals[local_eqn] +=
             -Vis_r *((h_film / (2.0 * GKE_Tau_P)) * interpolated_grad_Pv[i]+tangent_term * interpolated_t[i]/ (h_film * GKE_Tau_C)) * psif(l) * W * J;

          // We also need to worry about the jacobian terms
          if (flag) {
            // Loop over the nodes again, this time for the unknowns
            for (unsigned l2 = 0; l2 < n_node; l2++) {
              // Get the unknown pv_index
              local_unknown = this->nodal_local_eqn(l2, this->Pv_index[l2]);
              if (local_unknown >= 0) {
                jacobian(local_eqn, local_unknown) +=
                    Vis_r * (-psif(l2) * interpolated_n[i]) * psif(l) * W * J;
                jacobian(local_eqn, local_unknown) +=
                    -Vis_r *
                      ((h_film / (2.0 * GKE_Tau_P)) * dpsifdS(l2, i)) * psif(l) * W * J;
              }

              // Loop over the velocity components
              for (unsigned i2 = 0; i2 < n_dim; i2++) {
                // Get the unknown u index
                local_unknown = this->nodal_local_eqn(l2, u_index[i2]);

                // If not a boundary condition
                if (local_unknown >= 0) {
                  jacobian(local_eqn, local_unknown) +=
                      -Vis_r * (psif(l2) / (h_film * GKE_Tau_C)) *
                      interpolated_t[i2] * interpolated_t[i] * psif(l) * W * J;
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
void UpperInterfaceLubricationElement::output(std::ostream &outfile,
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
    outfile << interpolated_Pv(s) << " ";
    outfile << interpolated_Al(s) << " ";

    outfile << std::endl;
  }
}
} // namespace oomph
