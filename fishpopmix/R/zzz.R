##' @useDynLib fishpopmix
.onLoad <- function(libname, pkgname){

    .obj_schemes_ps <<- TMB::MakeADFun(list(fun=0,Scheme=0),list(Pf=1),DLL=pkgname)
    .tape_schemes_ps <<- RTMB::GetTape(.obj_schemes_ps,warn=FALSE)$atomic()

    .obj_logspace_add <<- TMB::MakeADFun(list(fun=1,add=1), list(logx=1,logy=0), DLL=pkgname)
    .tape_logspace_add <<- RTMB::GetTape(.obj_logspace_add,warn=FALSE)$atomic()

    .obj_logspace_sub <<- TMB::MakeADFun(list(fun=1,add=0), list(logx=1,logy=0), DLL=pkgname)
    .tape_logspace_sub <<- RTMB::GetTape(.obj_logspace_sub,warn=FALSE)$atomic()


    .obj_log_low_pnorm <<- TMB::MakeADFun(list(fun=2,lower_tail=1,log_p=1), list(x=2,y=0,mu=0,sigma=1),map=list(y=factor(NA)), DLL=pkgname)
    .tape_log_low_pnorm <<- RTMB::GetTape(.obj_log_low_pnorm,warn=FALSE)$atomic()

    .obj_log_up_pnorm <<- TMB::MakeADFun(list(fun=2,lower_tail=0,log_p=1), list(x=2,y=0,mu=0,sigma=1),map=list(y=factor(NA)), DLL=pkgname)
    .tape_log_up_pnorm <<- RTMB::GetTape(.obj_log_up_pnorm,warn=FALSE)$atomic()

    .obj_log_ipnorm <<- TMB::MakeADFun(list(fun=2,lower_tail=0,log_p=-1), list(x=-0.5,y=0.5,mu=0,sigma=1), DLL=pkgname)
    .tape_log_ipnorm <<- RTMB::GetTape(.obj_log_ipnorm,warn=FALSE)$atomic()

    .obj_quantreg <<- TMB::MakeADFun(list(fun=3), list(tau=0.5,x=1), DLL=pkgname)
    .tape_quantreg <<- RTMB::GetTape(.obj_quantreg,warn=FALSE)$atomic()

    
    invisible()
    
}
