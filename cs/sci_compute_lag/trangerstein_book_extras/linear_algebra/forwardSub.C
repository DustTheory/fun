void forwardSub(int m,float **L,float *by) {
//replaces by with the solution of L_{ik} y_k = b_i  
//assumes the diagonal entries of L are all one 
  for (int i=0;i<m;i++) {
    for (int j=0;j<i;j++) by[i] -= L[i][j] * by[j];
  }
}
