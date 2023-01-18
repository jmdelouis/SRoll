
#define RGSIZE (27664)

void Proj(double *Ax0,double *x0,double *him,int *hidx,double *sval,int *istart,int ndata,int nspline);

int HPR_Cleaner( int *hidx,      // in: one ring of {pixname}_HPRIDX_ABER_TotalFlag_dx11 TOI
                 float *data,    // in: one ring of TOI to project to Splined HPR (SHPR)
                 float *ohpr,    // out: one produced SHPR ring of length RGSIZE, allocated by the caller
                 int ndata,      // number of samples in <hidx> and <data>
                 char *message); // prefix of message to display at the end of function, set to NULL to display nothing
            
