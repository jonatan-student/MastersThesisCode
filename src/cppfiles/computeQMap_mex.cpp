#include <octave/oct.h>

static const int cxs[9]={0,0,1,1,1,0,-1,-1,-1};
static const int cys[9]={0,1,1,0,-1,-1,-1,0,1};

DEFUN_DLD(computeQMap_mex, args, , "computeQMap_mex(obstacle) -> q_map")
{
    boolNDArray obs(args(0).bool_array_value());
    octave_idx_type ny = obs.rows(), nx = obs.columns();
    dim_vector dv(ny,nx,9); NDArray qmap(dv,0.0);

    for(int k=1;k<9;++k){
        int cx=cxs[k], cy=cys[k];
        for(octave_idx_type y=1;y<ny-1;++y)
            for(octave_idx_type x=1;x<nx-1;++x)
                if(obs(y,x) && !obs(y+cy,x+cx))
                    qmap(y,x,k)=0.5;   // halfway bounce‑back
    }
    octave_value_list out; out(0)=qmap; return out;
}