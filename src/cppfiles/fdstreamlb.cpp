#include <octave/oct.h>

#define HELP ("\nfdstreamlb(f, cxs, cys)\n\
Performs streaming (propagation) step of Lattice Boltzmann method.\n")

DEFUN_DLD(fdstreamlb, args, , HELP) {
    octave_value_list retval;

    // Get inputs
    NDArray f(args(0).array_value());      // f(y, x, k)
    RowVector cxs(args(1).row_vector_value());
    RowVector cys(args(2).row_vector_value());

    const int ny = f.rows();
    const int nx = f.columns();
    const int nflows = f.dim3();

    NDArray f_out(f.dims());  // initialize output

    for (int k = 0; k < nflows; ++k) {
        int shift_x = static_cast<int>(cxs(k));
        int shift_y = static_cast<int>(cys(k));

        for (int y = 0; y < ny; ++y) {
            int src_y = (y - shift_y + ny) % ny;

            for (int x = 0; x < nx; ++x) {
                int src_x = (x - shift_x + nx) % nx;

                f_out(y, x, k) = f(src_y, src_x, k);
            }
        }
    }

    retval.append(f_out);
    return retval;
}
