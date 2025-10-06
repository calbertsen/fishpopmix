//#include <TMB.hpp>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

// #include "fromR.hpp"

//     template<class Float>
//     Float pnorm1_1x(Float x, Float lower_tail, Float log_p){
//       return fromR::pnorm5_raw(x,Float(0.0),Float(1.0),lower_tail,log_p);    
//     };


//     TMB_BIND_ATOMIC(pnorm1_2x,100,pnorm1_1x(x[0], x[1], x[2]))  

//     template<class Float>
//     Float log_ipnorm1_1x(Float x, Float y){
//       if(x > y) Rf_error("y should be larger than x!");
//       if(x < 0.5){ // Use lower tail
// 	Float v1 = fromR::pnorm5_raw(x,Float(0.0),Float(1.0),Float(1),Float(1));
// 	Float v2 = fromR::pnorm5_raw(y,Float(0.0),Float(1.0),Float(1),Float(1));
// 	return v2 + atomic::robust_utils::R_Log1_Exp(v1 - v2);
//       }else{ // Use upper tail
// 	Float v2 = fromR::pnorm5_raw(x,Float(0.0),Float(1.0),Float(0),Float(1));
// 	Float v1 = fromR::pnorm5_raw(y,Float(0.0),Float(1.0),Float(0),Float(1));
// 	return v2 + atomic::robust_utils::R_Log1_Exp(v1 - v2);
//       }    
//     };


//     TMB_BIND_ATOMIC(log_ipnorm1_2x,11,log_ipnorm1_1x(x[0], x[1]))  


  


//   template<class Type>
//   Type pnorm5(Type x, Type mu, Type sigma, int lower_tail, int log_p){
//     vector<Type> tx(4);
//     tx[0] = (x-mu) / sigma;
//     // tx[1] = mu;
//     // tx[2] = sigma;
//     tx[1] = (Type)lower_tail;
//     tx[2] = (Type)log_p;
//     tx[3] = 0; // extra argument for derivative order
//     Type res = pnorm1_2x(CppAD::vector<Type>(tx))[0];
//     return res;
//   }


//   template<class Type>
//   Type log_ipnorm(Type x, Type y, Type mu, Type sigma){
//     vector<Type> tx(3);
//     tx[0] = (x-mu) / sigma;
//     tx[1] = (y-mu) / sigma;
//     tx[2] = 0; // extra argument for derivative order
//     Type res = log_ipnorm1_2x(CppAD::vector<Type>(tx))[0];
//     return res;
//   }


//   template<class Float>
//   Float pt_raw(Float x, Float n, Float lower_tail, Float log_p){
//     return fromR::pt(x,n,(int)trunc(lower_tail),(int)trunc(log_p));    
//   };


// TMB_BIND_ATOMIC(pt_at,1100,pt_raw(x[0], x[1], x[2], x[3]))  

//  template<class Type>
//   Type pt(Type x, Type n, int lower_tail, int log_p){
//     vector<Type> tx(5);
//     tx[0] = x;
//     tx[1] = n;
//     // tx[2] = sigma;
//     tx[2] = (Type)lower_tail;
//     tx[3] = (Type)log_p;
//     tx[4] = 0; // extra argument for derivative order
//     Type res = pt_at(CppAD::vector<Type>(tx))[0];
//     return res;
//   }



// namespace rsam_pde_schemes {

//   namespace schemes_atomic {

//     template<class Float>
//     Float ps_raw (Float Pf) {
//       Float af = 0.0;
//       if(Pf > 10.0){
// 	af = (Pf - 1.0) / Pf;
//       }else if(Pf > 1e-6){
// 	af = (Pf-1.0 + pow(1.0-Pf/10.0,5)) / Pf;
//       }else if(Pf > -1e-6){
// 	af = 0.5;
//       }else if(Pf > -10.0){
// 	af = (pow(1.0+Pf/10.0,5) - 1.0) / Pf;
//       }else{
// 	af = -1.0 / Pf;
//       }
//       return af;
//     }
  
//     TMB_BIND_ATOMIC(ps0,
// 		    1,
// 		    ps_raw(x[0]) )
//   }

//   template<class Type>
//   Type ps(Type x) {  
//     vector<Type> tx(2);
//     tx[0] = x;
//     tx[1] = 0; // order
//     return schemes_atomic::ps0(CppAD::vector<Type>(tx))[0];
//   }

// }

// namespace rsam_logspace {

//   namespace logspace_atomic {
  
//     template<class Float>
//     Float logspace_add2_raw (Float logx, Float logy) {
//       // Was:
//       //  fmax2 (logx, logy) + log1p (exp (-fabs (logx - logy)));
//       if(logx == R_NegInf && logy == R_NegInf)
// 	return(R_NegInf);
//       if(logx == R_NegInf)
// 	return(logy);
//       if(logy == R_NegInf)
// 	return(logx);
//       return ( logx < logy ?
// 	       logy + log1p (exp (logx - logy)) :
// 	       logx + log1p (exp (logy - logx)) );
//     }

  
//     TMB_BIND_ATOMIC(logspace_add2,
// 		    11,
// 		    logspace_add2_raw(x[0], x[1]) )

//     template<class Float>
//     Float logspace_sub2_raw (Float logx, Float logy) {
//       if(logx == logy)
// 	return R_NegInf;
//       if(logy == R_NegInf)
// 	return(logx);
//       if(logx < logy)
// 	Rf_error("logx < logy in logspace_sub2");
//       return logx + atomic::robust_utils::R_Log1_Exp(logy - logx);
//     }

  
//     TMB_BIND_ATOMIC(logspace_sub2,
// 		    11,
// 		    logspace_sub2_raw(x[0], x[1]) )

 
//   }


//   template<class Type>
//   Type logspace_add2(Type logx, Type logy) {
//     if ( !CppAD::Variable(logx) && logx == Type(R_NegInf) && !CppAD::Variable(logy) && logy == Type(R_NegInf))
//       return Type(R_NegInf);
//     if ( !CppAD::Variable(logx) && logx == Type(R_NegInf) )
//       return logy;
//     if ( !CppAD::Variable(logy) && logy == Type(R_NegInf) )
//       return logx;
//     vector<Type> tx(3);
//     tx[0] = logx;
//     tx[1] = logy;
//     tx[2] = 0; // order
//     return logspace_atomic::logspace_add2(CppAD::vector<Type>(tx))[0];
//   }

  
//   template<class Type>
//   Type logspace_sub2(Type logx, Type logy) {
//     if ( !CppAD::Variable(logy) && logy == Type(R_NegInf) )
//       return logx;
//     vector<Type> tx(3);
//     tx[0] = logx;
//     tx[1] = logy;
//     tx[2] = 0; // order
//     return logspace_atomic::logspace_sub2(CppAD::vector<Type>(tx))[0];
//   }

// }



// namespace rsam_quantreg {

//   namespace quantreg_atomic {

//     template<class Float>
//     Float quantreg_raw (Float x, Float tau) {
//       if(x < 0)
// 	return x * (tau - 1.0);
//       return x * tau;
//     }
  
//     TMB_BIND_ATOMIC(quantreg,
// 		    11,
// 		    quantreg_raw(x[0], x[1]) )
//   }

//   template<class Type>
//   Type quantreg(Type x, Type tau) {  
//     vector<Type> tx(3);
//     tx[0] = x;
//     tx[1] = tau;
//     tx[2] = 0; // order
//     return quantreg_atomic::quantreg(CppAD::vector<Type>(tx))[0];
//   }

// }

// namespace log_bessel {

//   namespace log_bessel_atomic {

//     template<class Float>
//     Float login_log_besselI_raw(Float logx, Float nu){
//       if(logx == R_NegInf){
// 	if(fabs(nu) < 1e-16) return(0);
// 	return(R_NegInf);
//       }
//       if(logx == R_PosInf) return(R_PosInf);
//       if(nu < 0){
// 	if(fabs(nu - floor(nu)) > 1e-16)
// 	  Rf_error("Not implemented for negative fractions");
// 	nu = -nu;
//       }
//       Float r = R_NegInf;
//       Float dr = R_PosInf;
//       Float m = 0;
//       while(dr > -30){
// 	dr = (m * 2.0 + nu) * (logx - log(2.0));
// 	m += 1.0;
// 	//dr -= atomic::Rmath::D_lgamma(m,0.0) + atomic::Rmath::D_lgamma(m+nu,0.0);
// 	dr -= fromR::lgammafn(m) + fromR::lgammafn(m+nu);
// 	r = rsam_logspace::logspace_atomic::logspace_add2_raw(r,dr);
//       }
//       return r;
//     }
//    TMB_BIND_ATOMIC(login_log_besselI,
// 		    11,
// 		    login_log_besselI_raw(x[0], x[1]) )
   
//     template<class Float>
//     Float log_besselI_raw (Float x, Float nu) {
//       Float v = atomic::bessel_utils::bessel_i(x, nu, 2.0);
//       return log(v) + x;
//     }
    
//    TMB_BIND_ATOMIC(log_besselI,
// 		    11,
// 		    log_besselI_raw(x[0], x[1]) )


//    template<class Float>
//    Float log_MarcumQ_raw(Float a, Float b, Float nu){
//      if(fabs(b) < 1e-16 || a == R_PosInf || nu == R_PosInf) return 0.0;
//      if(b == R_NegInf) return 0.0;

//      Float r = R_NegInf;
//      Float dr = R_PosInf;
//      Float k = 1 - nu;
//      Float ab = a * b;
//      while(dr > -30){
//        dr = k * (log(a) - log(b)) + log_besselI_raw(ab,(Float)-k);
//        k += 1.0;
//        r = rsam_logspace::logspace_atomic::logspace_add2_raw(r,dr);
//      }
//      return -(a*a + b*b) / 2.0 + r;
//    }
    
//     TMB_BIND_ATOMIC(log_MarcumQ,
// 		    110,
// 		    log_MarcumQ_raw(x[0], x[1], x[2]) )


//  template<class Float>
//    Float login_log_MarcumQ_raw(Float loga, Float logb, Float nu){
//      if(logb == R_NegInf || loga == R_PosInf || nu == R_PosInf) return 0.0;

//      Float r = R_NegInf;
//      Float dr = R_PosInf;
//      Float k = 1 - nu;
//      Float logab = loga + logb;
//      while(dr > -30){
//        dr = k * (loga - logb) + login_log_besselI_raw(logab,(Float)-k);
//        k += 1.0;
//        r = rsam_logspace::logspace_atomic::logspace_add2_raw(r,dr);
//      }
//      return -exp(rsam_logspace::logspace_atomic::logspace_add2_raw(2.0 * loga, 2.0 * logb)) / 2.0 + r;
//    }
    
//     TMB_BIND_ATOMIC(login_log_MarcumQ,
// 		    110,
// 		    login_log_MarcumQ_raw(x[0], x[1], x[2]) )


    
   
//     template<class Float>
//     Float log_Marcum1mQ_raw(Float a, Float b, Float nu){
//       if(fabs(b) < 1e-16 || a == R_PosInf || nu == R_PosInf) return R_NegInf;
//       if(b == R_NegInf) return R_PosInf;

//       Float r = R_NegInf;
//       Float dr = R_PosInf;
//       Float alpha = nu;
//       Float ab = a * b;
//       while(dr > -30){
// 	dr = (Float)(alpha) * (log(b) - log(a)) + log_besselI_raw(ab,alpha);
// 	alpha += 1.0;
// 	r = rsam_logspace::logspace_atomic::logspace_add2_raw(r,dr);
//       }
//       return -(a*a + b*b) / 2.0 + r;
//     }
    
//    TMB_BIND_ATOMIC(log_Marcum1mQ,
// 		    110,
// 		   log_Marcum1mQ_raw(x[0], x[1], x[2]) )


//    template<class Float>
//     Float login_log_Marcum1mQ_raw(Float loga, Float logb, Float nu){
//       if(logb == R_NegInf || loga == R_PosInf || nu == R_PosInf) return R_NegInf;

//       Float r = R_NegInf;
//       Float dr = R_PosInf;
//       Float alpha = nu;
//       Float logab = loga + logb;
//       while(dr > -30){
// 	dr = (Float)(alpha) * (logb - loga) + login_log_besselI_raw(logab,alpha);
// 	alpha += 1.0;
// 	r = rsam_logspace::logspace_atomic::logspace_add2_raw(r,dr);
//       }
//       return -exp(rsam_logspace::logspace_atomic::logspace_add2_raw(2.0 * loga, 2.0 * logb)) / 2.0 + r;
//     }
    
//    TMB_BIND_ATOMIC(login_log_Marcum1mQ,
// 		    110,
// 		   login_log_Marcum1mQ_raw(x[0], x[1], x[2]) )

   
 
//   }

//   template<class Type>
//   Type login_log_besselI(Type logx, Type nu) {  
//     vector<Type> tx(3);
//     tx[0] = logx;
//     tx[1] = nu;
//     tx[2] = 0; // order
//     return log_bessel_atomic::login_log_besselI(CppAD::vector<Type>(tx))[0];
//   }


//   template<class Type>
//   Type log_besselI(Type x, Type nu) {  
//     vector<Type> tx(3);
//     tx[0] = x;
//     tx[1] = nu;
//     tx[2] = 0; // order
//     return log_bessel_atomic::log_besselI(CppAD::vector<Type>(tx))[0];
//   }

//     template<class Type>
//     Type log_Marcum1mQ(Type a, Type b, Type nu) {  
//     vector<Type> tx(4);
//     tx[0] = a;
//     tx[1] = b;
//     tx[2] = nu;
//     tx[3] = 0; // order
//     return log_bessel_atomic::log_Marcum1mQ(CppAD::vector<Type>(tx))[0];
//   }

//      template<class Type>
//     Type log_MarcumQ(Type a, Type b, Type nu) {  
//     vector<Type> tx(4);
//     tx[0] = a;
//     tx[1] = b;
//     tx[2] = nu;
//     tx[3] = 0; // order
//     return log_bessel_atomic::log_MarcumQ(CppAD::vector<Type>(tx))[0];
//   }

//    template<class Type>
//     Type login_log_Marcum1mQ(Type loga, Type logb, Type nu) {  
//     vector<Type> tx(4);
//     tx[0] = loga;
//     tx[1] = logb;
//     tx[2] = nu;
//     tx[3] = 0; // order
//     return log_bessel_atomic::login_log_Marcum1mQ(CppAD::vector<Type>(tx))[0];
//   }

//      template<class Type>
//     Type login_log_MarcumQ(Type loga, Type logb, Type nu) {  
//     vector<Type> tx(4);
//     tx[0] = loga;
//     tx[1] = logb;
//     tx[2] = nu;
//     tx[3] = 0; // order
//     return log_bessel_atomic::login_log_MarcumQ(CppAD::vector<Type>(tx))[0];
//   }

  
   
// }




// template<class Type>
// Type objective_function<Type>::operator() ()
// {
//   DATA_INTEGER(fun);
//   /*
//    * 0: PDE schemes
//    * 1: logspace
//    * 2: pnorm
//    * 3: quantreg
//    * 4: besselI
//    */

//   if(fun == 0){ // PDE schemes
//     DATA_INTEGER(Scheme);
//     PARAMETER(Pf);
//     if(Scheme == 0){
//       return rsam_pde_schemes::ps(Pf);
//     }else{
//       Rf_error("Scheme not implemented");
//     }
//     return 0.0;
//   }else if(fun == 1){ 		// logspace functions
//     DATA_INTEGER(add);
//     PARAMETER(logx);
//     PARAMETER(logy);

//     if(add == 1){
//       return rsam_logspace::logspace_add2(logx,logy);
//     }
//     return rsam_logspace::logspace_sub2(logx, logy);
    
//   }else if(fun == 2){
//     DATA_INTEGER(lower_tail);
//     DATA_INTEGER(log_p);
//     PARAMETER(x);
//     PARAMETER(y);
//     PARAMETER(mu);
//     PARAMETER(sigma);
//     if(log_p < 0)
//       return log_ipnorm(x,y,mu,sigma);
//     return pnorm5(x,mu,sigma,lower_tail,log_p);

//   }else if(fun == 3){
//     PARAMETER(tau);
//     PARAMETER(x);
//     return rsam_quantreg::quantreg(x,tau);
  
//   }else if(fun == 4){
//     PARAMETER(x);
//     PARAMETER(nu);
//     return log_bessel::log_besselI(x,nu);

//   }else if(fun == 5){
//     PARAMETER(a);
//     PARAMETER(b);
//     PARAMETER(nu);
//     return log_bessel::log_Marcum1mQ(a,b,nu);

//  }else if(fun == 6){
//     PARAMETER(a);
//     PARAMETER(b);
//     PARAMETER(nu);
//     return log_bessel::log_MarcumQ(a,b,nu);

//   }else if(fun == 7){
//     PARAMETER(logx);
//     PARAMETER(nu);
//     return log_bessel::login_log_besselI(logx,nu);

//  }else if(fun == 8){
//     PARAMETER(loga);
//     PARAMETER(logb);
//     PARAMETER(nu);
//     return log_bessel::login_log_Marcum1mQ(loga,logb,nu);

//  }else if(fun == 9){
//     PARAMETER(loga);
//     PARAMETER(logb);
//     PARAMETER(nu);
//     return log_bessel::login_log_MarcumQ(loga,logb,nu);

//   }else if(fun == 10){
//     DATA_INTEGER(lower_tail);
//     DATA_INTEGER(log_p);  
//     PARAMETER(x);
//     PARAMETER(df);
//     return pt(x,df,lower_tail,log_p);

    
    
//   }

  
  

  
//   return 0.0;
// }

// // SEXP getSetGlobalPtr(SEXP ptr) {
// // #ifdef TMBAD_FRAMEWORK
// //   SEXP global_ptr_tag = Rf_install("global_ptr");
// //   if (!Rf_isNull(ptr)) {
// //     SEXP tag = R_ExternalPtrTag(ptr);
// //     if (tag != global_ptr_tag) Rf_error("Invalid pointer type");
// //     TMBad::global_ptr = (TMBad::global**) R_ExternalPtrAddr(ptr);
// //   }
// //   SEXP res = R_MakeExternalPtr( (void*) TMBad::global_ptr, global_ptr_tag, R_NilValue);
// //   return res;
// // #else
// //   return R_NilValue;
// // #endif
// // }


extern "C" {

SEXP read_ped(SEXP file0, SEXP file1, SEXP nIndi0, SEXP nLoci0,
	      SEXP hasFID, SEXP hasParents, SEXP hasSex, SEXP hasPheno,
	      SEXP progress0, SEXP geno_NA);
  SEXP watershed(SEXP seed, SEXP priority);  
  SEXP kmeans_SLIC(SEXP pic, SEXP k, SEXP m, SEXP g, SEXP connect);  
  
#define CALLDEF(name,n) {#name, (DL_FUNC) &name, n}

  static const R_CallMethodDef CallEntries[] = {
    //TMB_CALLDEFS,
    CALLDEF(read_ped,10),    
    CALLDEF(watershed,2),
    CALLDEF(kmeans_SLIC,5),
    {NULL, NULL, 0}
  };

#define CALLABLE(name) R_RegisterCCallable("fishpopmix", #name, (DL_FUNC) &name)

  void R_init_fishpopmix(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    //TMB_CCALLABLES("fishpopmix");
  }

}
