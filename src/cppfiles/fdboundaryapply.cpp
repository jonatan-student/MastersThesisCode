#include <octave/oct.h>

#define HELP "\n fdboundaryapply(f, ux, uy, bndryF, cidy, cidx, slipLength) \n\
Applies slip-aware bounce-back distributions and zeroes velocities at boundaries.\n"

DEFUN_DLD(fdboundaryapply, args, , HELP) {
    octave_value_list retval;

    NDArray f(args(0).array_value());        // f(y, x, k)
    NDArray ux(args(1).array_value());       // ux(y, x)
    NDArray uy(args(2).array_value());       // uy(y, x)
    NDArray bndryF(args(3).array_value());   // extracted incoming populations
    RowVector cidy(args(4).row_vector_value()); // row indices (y)
    RowVector cidx(args(5).row_vector_value()); // col indices (x)
    NDArray slip(args(6).array_value());        // slip length field

    const int ny = f.rows();
    const int nx = f.columns();
    const int nflows = f.dim3();
    const int nb = cidy.numel();

    const int cxs[9] = {0, 0, 1, 1, 1, 0, -1, -1, -1};
    const int cys[9] = {0, 1, 1, 0, -1, -1, -1, 0, 1};
    const int opp[9] = {0, 5, 6, 7, 8, 1, 2, 3, 4};  // Opposite direction

    for (int idx = 0; idx < nb; ++idx) {
        int i = static_cast<int>(cidy(idx)) - 1;
        int j = static_cast<int>(cidx(idx)) - 1;

        double lambda = slip(i, j); // local slip length
        double alpha = (lambda > 0.0) ? lambda / (lambda + 0.5) : 0.0;

        for (int k = 0; k < nflows; ++k) {
            int k_opp = opp[k];

            int i_n = i + cys[k];
            int j_n = j + cxs[k];

            double f_bb = bndryF(idx, 0, k_opp);  // bounce-back direction
            double f_in = 0.0;
            if (i_n >= 0 && i_n < ny && j_n >= 0 && j_n < nx) {
                f_in = f(i_n, j_n, k);
            }
            double bounced = (1.0 - alpha) * f_bb + alpha * f_in;
            f(i, j, k) = bounced;
        }

        ux(i, j) = 0.0;
        uy(i, j) = 0.0;
    }

    retval.append(f);
    retval.append(ux);
    retval.append(uy);
    return retval;
}
