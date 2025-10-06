#include <R.h>
#include <Rinternals.h>
#include <stdio.h>
#include <ctype.h>
#include <R_ext/Utils.h>

struct WS_PriorityQueue {
  SEXP A; // List of lists A[[i]] = list(Key = {Priority}, Value = {Some SEXP, in this case a pixel index})
  int active_size;
  PROTECT_INDEX ipx_A;

  // Should add constructor with data
  WS_PriorityQueue(){
    PROTECT_WITH_INDEX(A = Rf_allocVector(VECSXP, 5000),&ipx_A);
    active_size = 0;
  }
  ~WS_PriorityQueue(){
    UNPROTECT(1);
  }
    
  int Size(){
    return active_size;
  }
  int Value(int i){
    if(i > active_size)
      Rf_error("Index too high");
    return Rf_asInteger(VECTOR_ELT(VECTOR_ELT(A,i-1),1));
  }
  double Key(int i){
    if(i > active_size)
      Rf_error("Index too high");
    return Rf_asReal(VECTOR_ELT(VECTOR_ELT(A,i-1),0));
  }
      
  int Parent(int i){
    return floor(i/2);
  }
  int Left(int i){
    return 2 * i;
  }
  int Right(int i){
    return 2 * i + 1;
  }

  void Swap(int i, int j){
    if(i > active_size || j > active_size)
      Rf_error("Index too high");
    SEXP tmp_i = PROTECT(VECTOR_ELT(A,i-1));
    SEXP tmp_j = PROTECT(VECTOR_ELT(A,j-1));
    SET_VECTOR_ELT(A,j-1,tmp_i);
    SET_VECTOR_ELT(A,i-1,tmp_j);
    UNPROTECT(2);
    return;
  }

  void DecreaseHeap(){
    SET_VECTOR_ELT(A,active_size-1,R_NilValue);
    active_size--;
    //REPROTECT(A = Rf_lengthgets(A,Size()-1),ipx_A);
  };
  void IncreaseHeap(){
    if(active_size + 1 > Rf_length(A)){
      REPROTECT(A = Rf_lengthgets(A,Rf_length(A)+5000),ipx_A);
    }
    active_size++;
  };
    
  void MaxHeapify(int i){
    int l = Left(i);
    int r = Right(i);
    int largest = i;
    if(l <= Size() && Key(l) > Key(i))
      largest = l;
    if(r <= Size() && Key(r) > Key(largest))
      largest = r;
    if(largest != i){
      Swap(i, largest);
      return MaxHeapify(largest);
    }
    return;	
  }

  SEXP PopNext(){
    if(Size() < 1)
      Rf_error("Underflow");
    SEXP maxv = PROTECT(VECTOR_ELT(A,0));
    SET_VECTOR_ELT(A,0,VECTOR_ELT(A,Size()-1));
    DecreaseHeap();
    MaxHeapify(1);
    UNPROTECT(1);
    return maxv;
  }

  int PopNextValue(){
    return Rf_asInteger(VECTOR_ELT(PopNext(),1));
  }

  void IncreaseKey(int i, double new_key){
    if(new_key < Key(i))
      Rf_error("New key is smaller than current key");
    SEXP tmp = PROTECT(VECTOR_ELT(A,i-1));
    SET_VECTOR_ELT(tmp,0,Rf_ScalarReal(new_key));
    SET_VECTOR_ELT(A,i-1,tmp);
    while(i > 1 and Key(Parent(i)) < Key(i)){
      Swap(i,Parent(i));
      i = Parent(i);
    }
    UNPROTECT(1);
    return;
  };

  void Insert(int new_value, double new_key){
    SEXP tmp = PROTECT(Rf_allocVector(VECSXP,2));
    SET_VECTOR_ELT(tmp,0,Rf_ScalarReal(R_NegInf));
    SET_VECTOR_ELT(tmp,1,Rf_ScalarInteger(new_value));
    IncreaseHeap();
    SET_VECTOR_ELT(A,Size()-1,tmp);
    UNPROTECT(1);
    return IncreaseKey(Size(),new_key);
  };

    
};



#ifdef DEBUG_PRINT
#define D_PRINT(...) Rprintf(__VA_ARGS__)
#else
#define D_PRINT(...)
#endif

extern "C" {

  SEXP watershed(SEXP seed, SEXP priority){
    D_PRINT("A\n");
    int nc = Rf_ncols(seed);
    int nr = Rf_nrows(seed);
    int nn = Rf_length(seed);
    D_PRINT("B\n");

    SEXP lab = PROTECT(Rf_allocMatrix(INTSXP,nr,nc));
    D_PRINT("C\n");
    Rf_copyMatrix(lab,seed,FALSE);
    D_PRINT("D\n");
    int* p_lab = INTEGER(lab);
    double* p_priority = REAL(priority);
    D_PRINT("E\n");
    /*
      mask:
      0: unprocessed
      1: in queue
      2: processed
    */
    SEXP mask = PROTECT(Rf_allocMatrix(INTSXP,nr,nc));
    int* p_mask = INTEGER(mask);
    D_PRINT("F\n");

    WS_PriorityQueue queue;
    D_PRINT("G\n");

#define i2r(i) (i % nr)
#define i2c(i) (i / nr)
#define rc2i(r,c) (r + c * nr)
#define init_neighbour(dr,dc)						\
    rx = rt+dr;								\
    cx = ct+dc;								\
    if(rx >= 0 && rx < nr && cx >= 0 && cx < nc){			\
      int ix = rc2i(rx,cx);						\
      if(p_lab[ix] == 0 && p_mask[ix] == 0){				\
	D_PRINT("\t\tInserting pixel %d (%d,%d) with priority %f into queue\n",ix,rx,cx,p_priority[ix]); \
	queue.Insert(ix, p_priority[ix]);				\
	p_mask[ix] = 1;							\
      }									\
    }
#define process_neighbour(dr,dc)					\
    rx = rt+dr;								\
    cx = ct+dc;								\
    if(rx >= 0 && rx < nr && cx >= 0 && cx < nc){			\
    int ix = rc2i(rx,cx);						\
    D_PRINT("\t\tProcessing neighbour to %d, pixel %d (%d,%d)\n",p_i,ix,rx,cx,p_priority[ix]); \
    if(p_lab[ix] == 0){							\
      if(p_mask[ix] == 0){						\
	D_PRINT("\t\t\tHas no label, inserted into queue with priority %f\n", p_priority[ix]); \
	queue.Insert(ix, p_priority[ix]);				\
	p_mask[ix] = 1;							\
      }else{								\
	D_PRINT("\t\t\tHas no label, but is already processed %f\n", p_priority[ix]); \
      }									\
    }else if(p_lab[ix] > 0){						\
      D_PRINT("\t\t\tHas a label\n");					\
      if(allNeighSameLabel){						\
      D_PRINT("\t\t\tAll so far have the same label\n");			\
      if(curlab == 0){							\
	D_PRINT("\t\t\tThis is the first!\n");				\
	curlab = p_lab[ix];						\
      }else{								\
	D_PRINT("\t\t\tThis is NOT the first!\n");			\
	if(curlab != p_lab[ix]){					\
	  D_PRINT("\t\t\tThis has a different label!\n");			\
	  allNeighSameLabel = false;					\
	  curlab = 0;							\
	}else{								\
	  D_PRINT("\t\t\tThis also has the same!\n");			\
	}								\
      }									\
      }else{								\
	D_PRINT("\t\t\tAlready found different labels.. do nothing!\n");	\
      }									\
    }else{								\
      D_PRINT("\t\t\tHow did we get here??...\n");			\
    }									\
  }

// Initialize mask
   for(int i = 0; i < nn; ++i){
      p_mask[i] = 0;
   }

    // Add neighbours of markers to queue
D_PRINT("H\n");
    for(int i = 0; i < nn; ++i){
      D_PRINT("\t%d\n",i);
      if(p_lab[i] > 0){
	D_PRINT("\tHELLO - pixel %d has a label!\n",i);
	p_mask[i] = 2;
	int rt = i2r(i);
	int ct = i2c(i);
	int rx, cx;
	//W
	init_neighbour(0,-1);
	//NW
	init_neighbour(1,-1);
	//N
	init_neighbour(1,0);
	//NE
	init_neighbour(1,1);
	//E
	init_neighbour(0,1);
	//SE
	init_neighbour(-1,1);
	//S
	init_neighbour(-1,0);
	//SW
	init_neighbour(-1,-1);
      }
    }
D_PRINT("I\n");
D_PRINT("Size of queue: %d\n",queue.Size());
    // Do watershed
    while(queue.Size() > 0){
      D_PRINT("\tSize of queue before: %d\n",queue.Size());
      int p_i = queue.PopNextValue();
      D_PRINT("\tSize of queue after: %d\n",queue.Size());
      D_PRINT("\t%d\n",p_i);
      
      int rt = i2r(p_i);
      int ct = i2c(p_i);
      int rx, cx;
      int curlab = 0;
      bool allNeighSameLabel = true;
      D_PRINT("\t\t%d\n",1);
      //W
      process_neighbour(0,-1);
      D_PRINT("\t\t%d\n",2);
      //NW
      process_neighbour(1,-1);
      D_PRINT("\t\t%d\n",3);
      //N
      process_neighbour(1,0);
      D_PRINT("\t\t%d\n",4);
      //NE
      process_neighbour(1,1);
      D_PRINT("\t\t%d\n",5);
      //E
      process_neighbour(0,1);
      D_PRINT("\t\t%d\n",6);
      //SE
      process_neighbour(-1,1);
      D_PRINT("\t\t%d\n",7);
      //S
      process_neighbour(-1,0);
      D_PRINT("\t\t%d\n",8);
      //SW
      process_neighbour(-1,-1);
      D_PRINT("\t\t%d\n",9);
      // Finalize
      D_PRINT("\t\tSet label of %d to %d\n",p_i,curlab);
      p_lab[p_i] = curlab;
      D_PRINT("\t\t%d\n",10);
      p_mask[p_i] = 2;
    }
    UNPROTECT(2);
    return lab;
  }
  
}
