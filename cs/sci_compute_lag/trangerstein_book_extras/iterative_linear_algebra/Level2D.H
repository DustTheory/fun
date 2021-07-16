#ifndef LEVEL2D_H
#define LEVEL2D_H "$Header:$"

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
    int n[2]; // number of grid cells in each coordinate direction
    int fc[2],lc[2]; // first & last corner indices; include boundary values
    int ifirst[2],ilast[2]; // exclude Dirichlet boundary values
    int fe[2],le[2]; // first & last cell indices; include boundary values
    int level_number;
    int number_corners;
    double mu;

    double *matrix;
    double *Ax_b;
    double *first_residual;
    double *solution_increment;
    double *prolongation_vector;
    double *prolongcell;
    double *prolongside0;
    double *prolongside1;

    Level(Level*);
    void setArrayBounds();
    void postSmooth() const ;
    void preSmooth() const ;
    void prolongation() const ;
    void restriction() const ;
    void setup(const double *m);
    void solve() const ;
    void updateResidual() const ;
  public:
    Level(const int*,const double*,int,SMOOTHER,RESTRICTION_PROLONGATION);
    ~Level();
    int numberLevels() const;
    Level *coarserLevel() const { return coarser; }
    Level *finerLevel() const { return finer; }
    void multigridStep(double *r,double *d) const;
    double computeResidual(const double *x,const double *b,double *r) 
      const;
    void checkProlongationAndRestrictionRandom() const;
    void checkProlongationAndRestriction() const;
    void checkMatrixSymmetryRandom() const;
    void checkMatrixSymmetry() const;
    void checkVCycleSymmetryRandom() const;
    void checkVCycleSymmetry() const;
    void checkSmoother() const;
    void checkSmootherUpdate() const;
    void checkSmootherUpdateRestrict() const;
    void checkCoarseGridProjection() const;
    void checkPointSource() const;
    void checkPointSourceRandom() const;
};
#endif
