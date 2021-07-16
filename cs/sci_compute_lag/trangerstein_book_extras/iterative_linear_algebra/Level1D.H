#ifndef LEVEL1D_H
#define LEVEL1D_H "$Header:$"

class Level {
  public:
    enum SMOOTHER{
      GAUSS_SEIDEL_SMOOTHER,GAUSS_SEIDEL_RED_BLACK_SMOOTHER,
      RICHARDSON_SMOOTHER
    };
    enum RESTRICTION_PROLONGATION{
      ALGEBRAIC_MULTIGRID_RESTRICTION_PROLONGATION,
      FINITE_ELEMENT_RESTRICTION_PROLONGATION
    };
    static const char *smoother_name[3];
    static const char *restriction_prolongation_name[2];
  private:
    static int smoother_iterations;
    static SMOOTHER smoother;
    static RESTRICTION_PROLONGATION restriction_prolongation;

    Level *coarser,*finer;
    int ncells;
    int level_number;
    int fi,la,stride;
    int ifirst,ilast;
    double mu;
//  double decay;

    double *matrix_c;
    double *Ax_b_c;
    double *first_residual_c;
    double *solution_increment_c;
    double *prolongation_vector_c;
    double *matrix_copy_c;
    double *richardson_residual_c;

    double *matrix;
    double *Ax_b;
    double *first_residual;
    double *solution_increment;
    double *prolongation_vector;
    double *matrix_copy;
    double *richardson_residual;

//  const double *error;

    Level(Level*);
    void postSmooth() const ;
    void preSmooth() const ;
    void prolongation() const ;
    void restriction() const ;
    void setup(const double*);
    void solve() const ;
    void updateResidual() const ;
  public:
    Level(int,const double*,int,SMOOTHER,RESTRICTION_PROLONGATION);
    ~Level();
    int numberLevels() const;
    Level *coarserLevel() const { return coarser; }
    Level *finerLevel() const { return finer; }
    void multigridStep(double *r,double *d) const;
    void computeResidual(const double *x,const double *b,double *r) 
      const;
    void checkProlongationAndRestrictionRandom() const;
    void checkProlongationAndRestriction() const;
    void checkSmootherRandom() const;
    void checkSmoother() const;
    void checkMatrixSymmetryRandom() const;
    void checkMatrixSymmetry() const;
    void checkVCycleSymmetryRandom() const;
    void checkVCycleSymmetry() const;
    void checkSmootherUpdate() const;
    void checkSmootherUpdateRestrict() const;
    void checkCoarseGridProjectionRandom() const;
    void checkCoarseGridProjection() const;
};
#endif
