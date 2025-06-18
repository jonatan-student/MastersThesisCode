#include <octave/oct.h>
#include <omp.h>

#define HELP ("\n")

DEFUN_DLD(fdcplb2, args, , HELP) {
	octave_value_list retval;

	NDArray f(args(0).array_value());           // f(y,x,k)
	NDArray b(args(1).array_value());           // bounce-back values
	Matrix ux(args(2).matrix_value());          // velocity x
	Matrix uy(args(3).matrix_value());          // velocity y
	RowVector cidy(args(4).row_vector_value()); // y indices (rows)
	RowVector cidx(args(5).row_vector_value()); // x indices (cols)
	NDArray slip(args(6).array_value());        // slip length matrix

	const int lenc = cidy.numel();

	const int cxs[9] = {0, 0, 1, 1, 1, 0, -1, -1, -1};  // x
	const int cys[9] = {0, 1, 1, 0, -1, -1, -1, 0, 1};  // y
	const int opp[9] = {0, 5, 6, 7, 8, 1, 2, 3, 4};

	const int ny = f.rows();
	const int nx = f.columns();

	const double tol = 1e-8;
  // Equilibrium weights for D2Q9
  	double w[9] = {4./9., 1./9., 1./36., 1./9., 1./36., 1./9., 1./36., 1./9., 1./36.};

	for (int idx = 0; idx < lenc; idx++) {
		int i = static_cast<int>(cidy(idx)) - 1; // y
		int j = static_cast<int>(cidx(idx)) - 1; // x

		double lambda = slip(i, j); // slip length at this boundary node
		double alpha = lambda / (lambda + 0.5);

		double ux_local = 0.0;
		double uy_local = 0.0;
		double rho_local = 0.0;

		for (int k = 0; k < 9; k++) {
			int k_opp = opp[k];

			int in = i + cys[k]; // neighbor row (y)
			int jn = j + cxs[k]; // neighbor col (x)

			double f_in = 0.0;
			if (in >= 0 && in < ny && jn >= 0 && jn < nx) {
				f_in = f(in, jn, k);
			}

			double bounced = (1.0 - alpha) * b(i, j, k) + alpha * f_in;
			f(i, j, k) = bounced;
			rho_local += bounced;
			ux_local += bounced * cxs[k];
			uy_local += bounced * cys[k];
		}

		// If no-slip is intended, then alpha == 0.
		// In that case, override the wall node’s distributions with equilibrium
		// values corresponding to zero velocity.
		if (fabs(alpha) <= tol) {
			// Use the computed (or expected) density at the wall.
			// If the bounce-back yields an extremely low density, then set a default.
			if (fabs(rho_local) <= tol)
				rho_local = 1.0;

			// Set f(i,j,k) = w[k] * rho_local so that the momentum is zero.
			for (int k = 0; k < 9; k++) {
				f(i, j, k) = w[k] * rho_local;
			}

			ux_local = 0.0;
			uy_local = 0.0;
		}

		ux(i, j) = ux_local / rho_local;
		uy(i, j) = uy_local / rho_local;
	}

	retval.append(f);
	retval.append(ux);
	retval.append(uy);
	return retval;
}
