/*--------------------------------------------------------------------
  steadyStateRDA_mex.cpp – serial replacement for SteadyStateRDA.m

  Compile with:
    mkoctfile -v -Wall -O3 -march=native -ffast-math \
              steadyStateRDA_mex.cpp
---------------------------------------------------------------------*/

#include <octave/oct.h>
#include <octave/Sparse.h>
#include <vector>

#define IDX(i,j,ngrdy) ( static_cast<octave_idx_type>((j)*(ngrdy) + (i)) )

DEFUN_DLD (steadyStateRDA_mex, args, ,
           "steadyStateRDA_mex(u_x, u_y, D, K, obstacle, reaction_mask, "
           "a_in, xmin, xmax)")
{
    // 0) Check & unpack inputs
    if (args.length() != 9)
        error("steadyStateRDA_mex: need exactly 9 input arguments");

    Matrix       u_x   = args(0).matrix_value();
    Matrix       u_y   = args(1).matrix_value();
    double       Dloc  = args(2).double_value();
    double       K     = args(3).double_value();
    boolNDArray  obst  = args(4).bool_array_value();
    boolNDArray  react = args(5).bool_array_value();
    double       a_in  = args(6).double_value();
    auto         xmin  = args(7).idx_type_value();
    auto         xmax  = args(8).idx_type_value();

    octave_idx_type ngrdy = u_x.rows();
    octave_idx_type ngrdx = u_x.columns();
    octave_idx_type N     = ngrdy * ngrdx;

    // 1) Build triplets for A and B
    octave_idx_type nnz_est = 15 * N;
    std::vector<octave_idx_type> TiA; TiA.reserve(nnz_est);
    std::vector<octave_idx_type> TjA; TjA.reserve(nnz_est);
    std::vector<double>          TvA; TvA.reserve(nnz_est);
    std::vector<octave_idx_type> TiB; TiB.reserve(nnz_est);
    std::vector<octave_idx_type> TjB; TjB.reserve(nnz_est);
    std::vector<double>          TvB; TvB.reserve(nnz_est);

    ColumnVector bA(N, 0.0), bB(N, 0.0);
    const int neigh[4][2] = {{1,0},{-1,0},{0,1},{0,-1}};

    for (octave_idx_type j = 0; j < ngrdx; ++j)
    {
        // progress
        octave_stdout << "\rAssembling column " << (j+1)
                      << " / " << ngrdx << std::flush;

        for (octave_idx_type i = 0; i < ngrdy; ++i)
        {
            octave_idx_type k = IDX(i,j,ngrdy);

            // obstacle: a=0, p=0
            if (obst(i,j))
            {
                TiA.push_back(k); TjA.push_back(k); TvA.push_back(1.0);
                TiB.push_back(k); TjB.push_back(k); TvB.push_back(1.0);
                continue;
            }

            // inlet (j==0) Dirichlet a=a_in, homogeneous p
            if (j == 0)
            {
                TiA.push_back(k); TjA.push_back(k); TvA.push_back(1.0);
                bA(k) = a_in;
                TiB.push_back(k); TjB.push_back(k); TvB.push_back(1.0);
                continue;
            }

            // outlet (j==ngrdx-1) Neumann: copy from j-1
            if (j == ngrdx-1)
            {
                octave_idx_type km1 = IDX(i,j-1,ngrdy);
                TiA.push_back(k); TjA.push_back(km1); TvA.push_back(1.0);
                TiB.push_back(k); TjB.push_back(km1); TvB.push_back(1.0);
                continue;
            }

            // interior
            double uxv = u_x(i,j), uyv = u_y(i,j);

            // advection x (A & B)
            if (uxv >= 0)
            {
                TiA.push_back(k); TjA.push_back(k);                TvA.push_back(-uxv);
                TiA.push_back(k); TjA.push_back(IDX(i,j-1,ngrdy)); TvA.push_back( uxv);
                TiB.push_back(k); TjB.push_back(k);                TvB.push_back(-uxv);
                TiB.push_back(k); TjB.push_back(IDX(i,j-1,ngrdy)); TvB.push_back( uxv);
            }
            else
            {
                TiA.push_back(k); TjA.push_back(k);                TvA.push_back( uxv);
                TiA.push_back(k); TjA.push_back(IDX(i,j+1,ngrdy)); TvA.push_back(-uxv);
                TiB.push_back(k); TjB.push_back(k);                TvB.push_back( uxv);
                TiB.push_back(k); TjB.push_back(IDX(i,j+1,ngrdy)); TvB.push_back(-uxv);
            }

            // advection y (A & B)
            if (uyv >= 0)
            {
                auto iu = std::max<octave_idx_type>(i-1,0);
                TiA.push_back(k); TjA.push_back(k);               TvA.push_back(-uyv);
                TiA.push_back(k); TjA.push_back(IDX(iu,j,ngrdy)); TvA.push_back( uyv);
                TiB.push_back(k); TjB.push_back(k);               TvB.push_back(-uyv);
                TiB.push_back(k); TjB.push_back(IDX(iu,j,ngrdy)); TvB.push_back( uyv);
            }
            else
            {
                auto id = std::min<octave_idx_type>(i+1,ngrdy-1);
                TiA.push_back(k); TjA.push_back(k);               TvA.push_back( uyv);
                TiA.push_back(k); TjA.push_back(IDX(id,j,ngrdy)); TvA.push_back(-uyv);
                TiB.push_back(k); TjB.push_back(k);               TvB.push_back( uyv);
                TiB.push_back(k); TjB.push_back(IDX(id,j,ngrdy)); TvB.push_back(-uyv);
            }

            // diffusion + reaction in A
            double diagA = -4.0*Dloc - (react(i,j) ? K : 0.0);
            TiA.push_back(k); TjA.push_back(k); TvA.push_back(diagA);
            for (int d=0; d<4; ++d)
            {
                auto ii = i + neigh[d][0];
                auto jj = j + neigh[d][1];
                if (ii<0||ii>=ngrdy||jj<0||jj>=ngrdx||obst(ii,jj))
                    TiA.push_back(k), TjA.push_back(k), TvA.push_back(Dloc);
                else
                    TiA.push_back(k), TjA.push_back(IDX(ii,jj,ngrdy)), TvA.push_back(Dloc);
            }

            // diffusion only in B
            TiB.push_back(k); TjB.push_back(k); TvB.push_back(-4.0*Dloc);
            for (int d=0; d<4; ++d)
            {
                auto ii = i + neigh[d][0];
                auto jj = j + neigh[d][1];
                if (ii<0||ii>=ngrdy||jj<0||jj>=ngrdx||obst(ii,jj))
                    TiB.push_back(k), TjB.push_back(k),    TvB.push_back(Dloc);
                else
                    TiB.push_back(k), TjB.push_back(IDX(ii,jj,ngrdy)), TvB.push_back(Dloc);
            }
        }
    }
    octave_stdout << "\n";  // end progress line

    // 2) Build Sparse A, B and add tiny identity
    SparseMatrix A(N, N, nnz_est);
    for (size_t m=0; m<TiA.size(); ++m)
        A( TiA[m], TjA[m] ) = TvA[m];
    for (octave_idx_type i=0; i<N; ++i) A(i,i) += 1e-12;

    SparseMatrix B(N, N, nnz_est);
    for (size_t m=0; m<TiB.size(); ++m)
        B( TiB[m], TjB[m] ) = TvB[m];
    for (octave_idx_type i=0; i<N; ++i) B(i,i) += 1e-12;

    // 3) Solve A * a = bA
    ColumnVector xA = A.solve(bA);
    Matrix a_ss(ngrdy, ngrdx);
    for (octave_idx_type j=0; j<ngrdx; ++j)
      for (octave_idx_type i=0; i<ngrdy; ++i)
        a_ss(i,j) = xA(IDX(i,j,ngrdy));

    // build bB = K * reaction_mask .* a_ss
    for (octave_idx_type j=0; j<ngrdx; ++j)
      for (octave_idx_type i=0; i<ngrdy; ++i)
        if (!obst(i,j) && react(i,j))
          bB(IDX(i,j,ngrdy)) = K * a_ss(i,j);

    // 4) Solve B * p = bB
    ColumnVector xB = B.solve(bB);
    Matrix p_ss(ngrdy, ngrdx);
    for (octave_idx_type j=0; j<ngrdx; ++j)
      for (octave_idx_type i=0; i<ngrdy; ++i)
        p_ss(i,j) = -xB(IDX(i,j,ngrdy));  // sign flip

    // optional diagnostic
    octave_stdout << " p_ss range: min=" << p_ss.min()
                  << ", max=" << p_ss.max() << "\n";

    // return
    octave_value_list out;
    out(0) = a_ss;
    out(1) = p_ss;
    return out;
}
