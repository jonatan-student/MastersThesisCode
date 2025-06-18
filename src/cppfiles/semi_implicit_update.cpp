#include <octave/oct.h>
#include <octave/parse.h>
#include <cmath>

// Semi-implicit Gauss-Seidel solver for diffusion: (I - dt*D*laplace) * a = rhs
// Note: only diffusion is treated implicitly

#define HELP "Usage: a_next = semi_implicit_update(a, rhs, dt, D, dx, num_iters)"

DEFUN_DLD(semi_implicit_update, args, , HELP)
{
    octave_value_list retval;

    if (args.length() < 6) {
        error("semi_implicit_update: Expected 6 arguments");
        return retval;
    }

    Matrix a(args(0).matrix_value());  // initial guess / previous value
    Matrix rhs(args(1).matrix_value());
    double dt = args(2).double_value();
    double D = args(3).double_value();
    double dx = args(4).double_value();
    int iters = args(5).int_value();

    size_t ny = a.rows();
    size_t nx = a.columns();

    double alpha = dt * D / (dx * dx);

    for (int it = 0; it < iters; ++it) {
        Matrix a_new = a;  // to hold updated values

        for (size_t j = 1; j < ny - 1; ++j) {
            for (size_t i = 1; i < nx - 1; ++i) {
                double lap = a(j,i+1) + a(j,i-1) + a(j+1,i) + a(j-1,i) - 4*a(j,i);
                a_new(j,i) = rhs(j,i) + alpha * lap;
            }
        }

        a = a_new;  // update for next iteration
    }

    retval.append(a);
    return retval;
}
