#include <octave/oct.h>

#define HELP "\nfdcollisionlb(f, tau, cxs, cys, weights)\n\
Computes rho, velocities, Feq, and relaxes towards equilibrium.\n"

DEFUN_DLD(fdcollisionlb, args, , HELP) {
    octave_value_list retval;

    // Input
    NDArray f(args(0).array_value());  // f(y, x, k)
    double tau = args(1).double_value();
    RowVector cxs(args(2).row_vector_value());
    RowVector cys(args(3).row_vector_value());
    RowVector weights(args(4).row_vector_value());

    const int ny = f.rows();
    const int nx = f.columns();
    const int nflows = f.dim3();

    // Allocate outputs
    Matrix rho(ny, nx, 0.0);
    Matrix ux(ny, nx, 0.0);
    Matrix uy(ny, nx, 0.0);

    NDArray feq(f.dims());

    // Main loop over grid
    for (int i = 0; i < ny; ++i) {
        for (int j = 0; j < nx; ++j) {

            double rho_local = 0.0;
            double ux_local = 0.0;
            double uy_local = 0.0;

            // Compute rho, ux, uy
            for (int k = 0; k < nflows; ++k) {
                double f_val = f(i, j, k);
                rho_local += f_val;
                ux_local += f_val * cxs(k);
                uy_local += f_val * cys(k);
            }

            if (rho_local > 0.0) {  // Avoid division by zero
                ux_local /= rho_local;
                uy_local /= rho_local;
            }

            rho(i, j) = rho_local;
            ux(i, j) = ux_local;
            uy(i, j) = uy_local;

            // Compute Feq and relax f
            for (int k = 0; k < nflows; ++k) {
                double c_dot_u = cxs(k) * ux_local + cys(k) * uy_local;
                double u_sq = ux_local * ux_local + uy_local * uy_local;

                double feq_val = weights(k) * rho_local *
                    (1.0 + 3.0 * c_dot_u + 4.5 * c_dot_u * c_dot_u - 1.5 * u_sq);

                f(i, j, k) -= (1.0 / tau) * (f(i, j, k) - feq_val);
            }
        }
    }

    retval.append(f);
    retval.append(rho);
    retval.append(ux);
    retval.append(uy);
    return retval;
}