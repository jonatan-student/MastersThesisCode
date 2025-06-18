#include <octave/oct.h>

#define HELP "\n fdboundaryextract(f, cidy, cidx) \n\
Extracts incoming distributions from neighboring fluid nodes for bounce-back.\n"

DEFUN_DLD(fdboundaryextract, args, , HELP) {
    octave_value_list retval;

    NDArray f(args(0).array_value());  // f(y, x, k)
    RowVector cidy(args(1).row_vector_value());  // row indices (y)
    RowVector cidx(args(2).row_vector_value());  // col indices (x)

    const int ny = f.rows();
    const int nx = f.columns();
    const int nflows = f.dim3();
    const int nb = cidy.numel();

    dim_vector out_dims(nb, 1, nflows);
    NDArray bndryF(out_dims);

    const int cxs[9] = {0, 0, 1, 1, 1, 0, -1, -1, -1};
    const int cys[9] = {0, 1, 1, 0, -1, -1, -1, 0, 1};

    for (int idx = 0; idx < nb; ++idx) {
        int i = static_cast<int>(cidy(idx)) - 1;
        int j = static_cast<int>(cidx(idx)) - 1;

        for (int k = 0; k < nflows; ++k) {
            int i_n = i + cys[k]; // neighbor y
            int j_n = j + cxs[k]; // neighbor x

            double f_in = 0.0;
            if (i_n >= 0 && i_n < ny && j_n >= 0 && j_n < nx) {
                f_in = f(i_n, j_n, k);
            }

            bndryF(idx, 0, k) = f(i, j, k);
        }
    }

    retval.append(bndryF);
    return retval;
}