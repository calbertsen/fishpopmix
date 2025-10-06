#include <R.h>
#include <Rinternals.h>
#include <stdio.h>
#include <ctype.h>
#include <cmath>
#include <R_ext/Utils.h>
#include <queue>

#ifdef DEBUG_PRINT
#define D_PRINT(...) Rprintf(__VA_ARGS__)
#else
#define D_PRINT(...) 
#endif


// Input: Matrix of labels
// Output: Image of labels cc, histogram H of cc
// Alg6 https://enpc.hal.science/hal-03651336/document 
SEXP connectedComponentLabeling(SEXP l){

}

// Alg1


// Alg4


// Alg5
extern "C" {
  SEXP kmeans_SLIC(SEXP pic, SEXP k, SEXP m, SEXP g, SEXP connect){
    // Initialize centers
    int ki = Rf_asInteger(k);
    double md = Rf_asReal(m); 
    int gi = Rf_asInteger(g);
    int nc = Rf_ncols(pic);
    int nr = Rf_nrows(pic);
    int nn = nc * nr;
    D_PRINT("%d %d %d\n",nc,nr,nn);
#define imR(i,j) pic_p[i + j * nr + 0 * (nr*nc)]
#define imG(i,j) pic_p[i + j * nr + 1 * (nr*nc)]
#define imB(i,j) pic_p[i + j * nr + 2 * (nr*nc)]
#define max(a,b) ((a<b)?b:a)
#define SQUARE(x) ((x)*(x))
    int S = max(1, (int)sqrt(nn / (double)ki));
    int np_c = std::floor(nc / (double)S);
    int np_r = std::floor(nr / (double)S);
    int np = np_c * np_r;
    int pad_c = nc - S * np_c;
    int pad_r = nr - S * np_r;
    D_PRINT("%d %d %d %d %d\n",np_c,np_r,np,pad_c,pad_r);
    SEXP X = PROTECT(Rf_allocMatrix(REALSXP,np,5)); // x, y, R, G, B
    double* X_p = REAL(X);
    double* pic_p = REAL(pic);
    for(int i = 0; i < np_r; ++i)
      for(int j = 0; j < np_c; ++j){
	D_PRINT("%d %d\n",i,j);
	int indx = i + j * np_r;
	int spX = i * S + S/2.0 + pad_r/2.0;
	int spY = j * S + S/2.0 + pad_c/2.0;
	// Gradient correction
	if(gi > 0){
	  double mv = R_PosInf;
	  for(int ki = spX-gi; ki <= spX+gi; ++ki)
	    for(int kj = spY-gi; kj <= spY+gi; ++kj){
	      if(ki >= nr || ki < 0 || kj >= nc || kj < 0) // Outside image
		continue;
	      int ki2 = (ki+1<nr)?(ki + 1):(ki-1);
	      int kj2 = (kj+1<nc)?(kj + 1):(kj-1);
	      double G = SQUARE(imR(ki,kj)-imR(ki2,kj)) + SQUARE(imR(ki,kj)-imR(ki2,kj)) + SQUARE(imR(ki,kj)-imR(ki2,kj)) +
		SQUARE(imR(ki,kj)-imR(ki,kj2)) + SQUARE(imR(ki,kj)-imR(ki,kj2)) + SQUARE(imR(ki,kj)-imR(ki,kj2));
	      if(G < mv){
		mv = G;
		spX = ki;
		spY = kj;
	      }
	    }
	}
	// SP X (row)
	X_p[indx] = spX;
	// SP Y (col)
	X_p[indx + 1 * np] = spY;
	// SP R
	X_p[indx + 2 * np] = imR(spX,spY);
	// SP G
	X_p[indx + 3 * np] = imG(spX,spY);
	// SP B
	X_p[indx + 4 * np] = imB(spX,spY);	
      }
    D_PRINT("A\n");
    double E = R_PosInf;
    SEXP l = PROTECT(Rf_allocMatrix(INTSXP,nr,nc));
    int* l_p = INTEGER(l);
    D_PRINT("B\n");
    for(int i = 0; i < nr*nc; ++i)
      l_p[i] = -1;
    D_PRINT("C\n");   
    D_PRINT("D\n");
    // Initialize distance matrix
    SEXP d = PROTECT(Rf_allocMatrix(REALSXP,nr,nc));
    double* d_p = REAL(d);
    D_PRINT("E\n");
    for(int i = 0; i < nr*nc; ++i)
      d_p[i] = R_PosInf;
    SEXP spX2 = PROTECT(Rf_allocVector(REALSXP, np));
    double* spX2_p = REAL(spX2);    
    SEXP spY2 = PROTECT(Rf_allocVector(REALSXP, np));
    double* spY2_p = REAL(spY2);
    SEXP spR2 = PROTECT(Rf_allocVector(REALSXP, np));
    double* spR2_p = REAL(spR2);    
    SEXP spG2 = PROTECT(Rf_allocVector(REALSXP, np));
    double* spG2_p = REAL(spG2);
    SEXP spB2 = PROTECT(Rf_allocVector(REALSXP, np));
    double* spB2_p = REAL(spB2);
    SEXP spN = PROTECT(Rf_allocVector(INTSXP, np));
    int* spN_p = INTEGER(spN);
    while(E > 0.5){
      D_PRINT("F\n");
      // Loop over superpixels
      for(int k = 0; k < np; ++k){
	int spX = X_p[k];
	int spY = X_p[k + 1 * np];
	double spR = X_p[k + 2 * np];
	double spG = X_p[k + 3 * np];
	double spB = X_p[k + 4 * np];
	// Loop over nearby pixels
	for(int pi = spX-S; pi <= spX+S; ++pi)
	  for(int pj = spY-S; pj <= spY+S; ++pj){
	    if(pi >= nr || pi < 0 || pj >= nc || pj < 0) // Outside image
	      continue;
	    double e = md*md * (SQUARE(spX-pi)+SQUARE(spY-pj)) / (S*S) + (SQUARE(spR-imR(pi,pj))+SQUARE(spG-imG(pi,pj))+SQUARE(spB-imB(pi,pj)));
	    if(d_p[pi+pj*nr] > e){
	      d_p[pi+pj*nr] = e;
	      l_p[pi+pj*nr] = k;
	    }
	  }
      }
      D_PRINT("G\n");
      // Update step
      D_PRINT("H\n");
      for(int i = 0; i < np; ++i){
	spX2_p[i] = 0;
	spY2_p[i] = 0;
	spR2_p[i] = 0;
	spG2_p[i] = 0;
	spB2_p[i] = 0;
	spN_p[i] = 0;
      }
      D_PRINT("I\n");
      for(int p = 0; p < nr*nc; ++p){
	if(l_p[p] >= 0 && l_p[p] < np){
	  int i = p % nr;
	  int j = p / nr;
	  spX2_p[l_p[p]] += i;
	  spY2_p[l_p[p]] += j;
	  spR2_p[l_p[p]] += imR(i,j);
	  spG2_p[l_p[p]] += imG(i,j);
	  spB2_p[l_p[p]] += imB(i,j);
	  spN_p[l_p[p]] += 1;
	}
      }
      D_PRINT("J\n");
      // Compute shift + overwrite
      double Ds = 0;
      for(int i = 0; i < np; ++i){
	D_PRINT("\t%d:\n",i);
	double xO = X_p[i];
	D_PRINT("\t1\n");
	double yO = X_p[i + 1 * np];
	double rO = X_p[i + 2 * np];
	double gO = X_p[i + 3 * np];
	double bO = X_p[i + 4 * np];
	double xN = xO;
	double yN = yO;
	double rN = rO;
	double gN = gO;
	double bN = bO;
	if(spN_p[i] > 0){
	  xN = spX2_p[i] / (double)spN_p[i];
	  yN = spY2_p[i] / (double)spN_p[i];
	  rN = spR2_p[i] / (double)spN_p[i];
	  gN = spG2_p[i] / (double)spN_p[i];
	  bN = spB2_p[i] / (double)spN_p[i];
	}
	D_PRINT("\t2 New: %f, %f\n",xN,yN);
	D_PRINT("\t2 Old: %f, %f\n",xO,yO);
	Ds += SQUARE(xO-xN) + SQUARE(yO-yN);
	D_PRINT("\tS1: %f; S2: %f; Dsk: %f\n",SQUARE(xO-xN), SQUARE(yO-yN),SQUARE(xO-xN) + SQUARE(yO-yN));
	// SP X (row)
	X_p[i] = xN;
	// SP Y (col)
	X_p[i + 1 * np] = yN;
	// SP R
	X_p[i + 2 * np] = rN;
	// SP G
	X_p[i + 3 * np] = gN;
	// SP B
	X_p[i + 4 * np] = bN;
	D_PRINT("\t3\n");
	    
      }
      D_PRINT("K\n");
      D_PRINT("Ds: %f, np: %f\n",Ds,(double)np);
      E = sqrt(Ds / (double)np);
      D_PRINT("E: %f\n",E);
      //UNPROTECT(6);
    }
    // Conectivity
    SEXP cc = PROTECT(Rf_allocMatrix(INTSXP,nr,nc));
    int* cc_p = INTEGER(cc);
    SEXP H = PROTECT(Rf_allocVector(INTSXP,nr*nc));
    int* H_p = INTEGER(H);
    SEXP M = PROTECT(Rf_allocVector(INTSXP,np));
    int* M_p = INTEGER(M);
    if(Rf_asLogical(connect)){
      D_PRINT("L\n");
      std::queue<int> C;
      // Alg 6: Connected component labeling
      D_PRINT("M\n");
      for(int i = 0; i < nr*nc; ++i)
	cc_p[i] = -1;
      int n = 0;
      for(int p = 0; p < nr*nc; ++p){
	if(cc_p[p] == -1){
	  cc_p[p] = n;
	  H_p[n] = 0;
	  C.push(p);
	  while (!C.empty()) {
	    int q = C.front();
	    C.pop();
	    H_p[n]++;
	    int q_i = q % nr;
	    int q_j = q / nr;
	    int r_i, r_j, r;
	    //W
	    r_i = q_i + 0;
	    r_j = q_j - 1;
	    r = r_i + r_j * nr;
	    if(r_i >= 0 && r_i < nr && r_j >= 0 && r_j < nc){
	      if(l_p[r] == l_p[q] && cc_p[r] == -1){
		cc_p[r] = n;
		C.push(r);
	      }
	    }
	    //N
	    r_i = q_i + 1;
	    r_j = q_j + 0;
	    r = r_i + r_j * nr;
	    if(r_i >= 0 && r_i < nr && r_j >= 0 && r_j < nc){
	      if(l_p[r] == l_p[q] && cc_p[r] == -1){
		cc_p[r] = n;
		C.push(r);
	      }
	    }
	    //E
	    r_i = q_i + 0;
	    r_j = q_j + 1;
	    r = r_i + r_j * nr;
	    if(r_i >= 0 && r_i < nr && r_j >= 0 && r_j < nc){
	      if(l_p[r] == l_p[q] && cc_p[r] == -1){
		cc_p[r] = n;
		C.push(r);
	      }
	    }
	    //S
	    r_i = q_i - 1;
	    r_j = q_j + 0;
	    r = r_i + r_j * nr;
	    if(r_i >= 0 && r_i < nr && r_j >= 0 && r_j < nc){
	      if(l_p[r] == l_p[q] && cc_p[r] == -1){
		cc_p[r] = n;
		C.push(r);
	      }
	    }
	  }
	  n++;
	}    
      }
      D_PRINT("N\n");
      // Alg7: Orphans
      for(int k = 0; k < np; ++k){
	M_p[k] = -1;
      }
      for(int p = 0; p < nr*nc; ++p){
	if(M_p[l_p[p]] == -1 || (M_p[l_p[p]] >= 0 && H_p[M_p[l_p[p]]] > H_p[cc_p[p]])) // 
	  M_p[l_p[p]] = cc_p[p];    
      }
      for(int p = 0; p < nr*nc; ++p){
	if(cc_p[p] != M_p[l_p[p]])
	  l_p[p] = R_NegInf;
      }
      D_PRINT("O\n");
      // Alg8: Distance of orphans to labeled pixels
      std::queue<int> Q;
      bool hasChanges = false;
      int dd = -1;
      do{
	hasChanges = false;
	for(int p = 0; p < nr*nc; ++p){
	  if(l_p[p] > dd){
	    int p_i = p % nr;
	    int p_j = p / nr;
	    int q_i, q_j, q;
	    //W
	    q_i = p_i + 0;
	    q_j = p_j - 1;
	    q = q_i + q_j * nr;
	    if(q_i >= 0 && q_i < nr && q_j >= 0 && q_j < nc){
	      if(l_p[q] < dd){
		l_p[q] = dd;
		Q.push(q);
		hasChanges = true;
	      }
	    }
	    //N
	    q_i = p_i + 1;
	    q_j = p_j + 0;
	    q = q_i + q_j * nr;
	    if(q_i >= 0 && q_i < nr && q_j >= 0 && q_j < nc){
	      if(l_p[q] < dd){
		l_p[q] = dd;
		Q.push(q);
		hasChanges = true;
	      }
	    }
	    //E
	    q_i = p_i + 0;
	    q_j = p_j + 1;
	    q = q_i + q_j * nr;
	    if(q_i >= 0 && q_i < nr && q_j >= 0 && q_j < nc){
	      if(l_p[q] < dd){
		l_p[q] = dd;
		Q.push(q);
		hasChanges = true;
	      }
	    }
	    //S
	    q_i = p_i - 1;
	    q_j = p_j + 0;
	    q = q_i + q_j * nr;
	    if(q_i >= 0 && q_i < nr && q_j >= 0 && q_j < nc){
	      if(l_p[q] < dd){
		l_p[q] = dd;
		Q.push(q);
		hasChanges = true;
	      }
	    }
	  }
	}
	dd--;
      }while(hasChanges);
      D_PRINT("P\n");
      // Alg 9: Adoption
      while(!Q.empty()){
	int p = Q.front();
	Q.pop();
	int p_i = p % nr;
	int p_j = p / nr;
	int q_i, q_j, q;
	double Dc = R_PosInf;
	int qkeep;
	//W
	q_i = p_i + 0;
	q_j = p_j - 1;
	q = q_i + q_j * nr;
	if(q_i >= 0 && q_i < nr && q_j >= 0 && q_j < nc && l_p[q] > -1){
	  double rQ = X_p[q + 2 * np];
	  double gQ = X_p[q + 3 * np];
	  double bQ = X_p[q + 4 * np];
	  double rP = pic_p[p + 0 * nr];
	  double gP = pic_p[p + 1 * nr];
	  double bP = pic_p[p + 2 * nr];
	  double tmp = SQUARE(rQ-rP) + SQUARE(gQ-gP) + SQUARE(bQ-bP);
	  if(tmp < Dc){
	    Dc = tmp;
	    qkeep = q;
	  }
	}
	//N
	q_i = p_i + 1;
	q_j = p_j + 0;
	q = q_i + q_j * nr;
	if(q_i >= 0 && q_i < nr && q_j >= 0 && q_j < nc && l_p[q] > -1){
	  double rQ = X_p[q + 2 * np];
	  double gQ = X_p[q + 3 * np];
	  double bQ = X_p[q + 4 * np];
	  double rP = pic_p[p + 0 * nr];
	  double gP = pic_p[p + 1 * nr];
	  double bP = pic_p[p + 2 * nr];
	  double tmp = SQUARE(rQ-rP) + SQUARE(gQ-gP) + SQUARE(bQ-bP);
	  if(tmp < Dc){
	    Dc = tmp;
	    qkeep = q;
	  }
	}
	//E
	q_i = p_i + 0;
	q_j = p_j + 1;
	q = q_i + q_j * nr;
	if(q_i >= 0 && q_i < nr && q_j >= 0 && q_j < nc && l_p[q] > -1){
	  double rQ = X_p[q + 2 * np];
	  double gQ = X_p[q + 3 * np];
	  double bQ = X_p[q + 4 * np];
	  double rP = pic_p[p + 0 * nr];
	  double gP = pic_p[p + 1 * nr];
	  double bP = pic_p[p + 2 * nr];
	  double tmp = SQUARE(rQ-rP) + SQUARE(gQ-gP) + SQUARE(bQ-bP);
	  if(tmp < Dc){
	    Dc = tmp;
	    qkeep = q;
	  }
	}
	//S
	q_i = p_i - 1;
	q_j = p_j + 0;
	q = q_i + q_j * nr;
	if(q_i >= 0 && q_i < nr && q_j >= 0 && q_j < nc && l_p[q] > -1){
	  double rQ = X_p[q + 2 * np];
	  double gQ = X_p[q + 3 * np];
	  double bQ = X_p[q + 4 * np];
	  double rP = pic_p[p + 0 * nr];
	  double gP = pic_p[p + 1 * nr];
	  double bP = pic_p[p + 2 * nr];
	  double tmp = SQUARE(rQ-rP) + SQUARE(gQ-gP) + SQUARE(bQ-bP);
	  if(tmp < Dc){
	    Dc = tmp;
	    qkeep = q;
	  }
	}
	l_p[p] = l_p[qkeep];
      }
      D_PRINT("Q\n");
    }
    // Update superpixel values
    for(int i = 0; i < np; ++i){
      spX2_p[i] = 0;
      spY2_p[i] = 0;
      spR2_p[i] = 0;
      spG2_p[i] = 0;
      spB2_p[i] = 0;
      spN_p[i] = 0;
    }
    D_PRINT("R\n");
    for(int p = 0; p < nr*nc; ++p){
      if(l_p[p] >= 0 && l_p[p] < np){
	int i = p % nr;
	int j = p / nr;
	spX2_p[l_p[p]] += i;
	spY2_p[l_p[p]] += j;
	spR2_p[l_p[p]] += imR(i,j);
	spG2_p[l_p[p]] += imG(i,j);
	spB2_p[l_p[p]] += imB(i,j);
	spN_p[l_p[p]] += 1;
      }
    }
    D_PRINT("S\n");
    for(int i = 0; i < np; ++i){
      double xO = X_p[i];
      double yO = X_p[i + 1 * np];
      double rO = X_p[i + 2 * np];
      double gO = X_p[i + 3 * np];
      double bO = X_p[i + 4 * np];
      double xN = xO;
      double yN = yO;
      double rN = rO;
      double gN = gO;
      double bN = bO;
      if(spN_p[i] > 0){
	xN = spX2_p[i] / (double)spN_p[i];
	yN = spY2_p[i] / (double)spN_p[i];
	rN = spR2_p[i] / (double)spN_p[i];
	gN = spG2_p[i] / (double)spN_p[i];
	bN = spB2_p[i] / (double)spN_p[i];
      }
      // SP X (row)
      X_p[i] = xN;
      // SP Y (col)
      X_p[i + 1 * np] = yN;
      // SP R
      X_p[i + 2 * np] = rN;
      // SP G
      X_p[i + 3 * np] = gN;
      // SP B
      X_p[i + 4 * np] = bN;
    }
    D_PRINT("T\n");
    // Make output image
    SEXP picO = PROTECT(Rf_alloc3DArray(REALSXP, nr, nc, 3));
    double* picO_p = REAL(picO);
    for(int p = 0; p < nr*nc; ++p){
      int sp = l_p[p];
      if(sp >= 0 && sp < np){
	// R
	picO_p[p + 0 * (nr*nc)] = X_p[sp + 2 * np];
	// G
	picO_p[p + 1 * (nr*nc)] = X_p[sp + 3 * np];
	// B
	picO_p[p + 2 * (nr*nc)] = X_p[sp + 4 * np];
      }else{
	// R
	picO_p[p + 0 * (nr*nc)] = 0;
	// G
	picO_p[p + 1 * (nr*nc)] = 0;
	// B
	picO_p[p + 2 * (nr*nc)] = 0;
      }
    }
  
    const char *names[] = {"super_pixels", "labels", "pic","cc","H","M", ""};
    SEXP list = PROTECT(Rf_mkNamed(VECSXP, names));
    SET_VECTOR_ELT(list, 0, X);
    SET_VECTOR_ELT(list, 1, l);
    SET_VECTOR_ELT(list, 2, picO);
    if(Rf_asLogical(connect)){
      SET_VECTOR_ELT(list, 3, cc); //cc);
      SET_VECTOR_ELT(list, 4, H); // H
      SET_VECTOR_ELT(list, 5, M); // M
    }else{
      SET_VECTOR_ELT(list, 3, R_NilValue); //cc);
      SET_VECTOR_ELT(list, 4, R_NilValue); // H
      SET_VECTOR_ELT(list, 5, R_NilValue); // M
    }
    UNPROTECT(14);
    return list;
  }

}
