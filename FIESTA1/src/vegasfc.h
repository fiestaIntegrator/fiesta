#ifndef VEGASF_H_
#define VEGASF_H_

#ifdef __cplusplus
extern "C" {
#endif

/************ FORTRAN wrapper: **************/
  extern struct {
    double s1; // integral from all iterations 
    double s2; // err. from all iterations
    double s3; //chi2
    double s4; // integral from last iteration
    double s5; // err. from last iteration
  } RESULT;

  extern struct {
    double calls;
    double ti;
    double tsi;
  } bveg4_;

#define FXN double (*fxn)(double* x, double* wgt)
  
  void VEGAS(FXN,double *,int *,int*,int*,int*,int*);
  void VEGAS1(FXN,double *,int *,int*,int*,int*,int*);
  void VEGAS2(FXN,double *,int *,int*,int*,int*,int*);
  void VEGAS3(FXN,double *,int *,int*,int*,int*,int*);


#ifdef __cplusplus
}
#endif

#endif /*VEGAS_H_*/
