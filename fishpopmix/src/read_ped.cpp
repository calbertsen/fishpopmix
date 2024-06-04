#include <R.h>
#include <Rinternals.h>
#include <stdio.h>
#include <ctype.h>
#include <R_ext/Utils.h>


extern "C" {

  SEXP read_ped(SEXP file0, SEXP file1, SEXP nIndi0, SEXP nLoci0,
		SEXP hasFID, SEXP hasParents, SEXP hasSex, SEXP hasPheno,
		SEXP progress0, SEXP geno_NA){
    const char* file = R_CHAR(Rf_asChar(file0));
    const char* fileM = R_CHAR(Rf_asChar(file1));
    FILE *fp;
    FILE *stream;

    //int buffsize = Rf_asInteger(buffer_size);
    fp = fopen(R_ExpandFileName(file), "r");
    if(fp == NULL)
      Rf_error("PED File could not be opened.");
    rewind(fp);
    
    char buff[2048];
    memset(buff,'\0',sizeof(buff));

    // fscanf(fp, "%s", buff);
    // printf("1 : %s\n", buff );

    const char* NA_char = R_CHAR(Rf_asChar(geno_NA));

    int nIndi = Rf_asInteger(nIndi0); // Adapt dynamically
    int nLoci = Rf_asInteger(nLoci0); // Adapt dynamically during first individual
    int i_indi = 0;
    int i_place = 0;
    int p_loci = 0;

    SEXP locusNames;
    PROTECT(locusNames = Rf_allocVector(STRSXP, nLoci));
    FILE *fpM;
    fpM = fopen(R_ExpandFileName(fileM), "r");
    if(fpM == NULL){
      locusNames = R_NilValue;
      if(file1 != R_NilValue)
	Rf_warning("Map file could not be opened.");
    }else{
      rewind(fpM);
      int ll = 0;
      char cBu[255];
      char vID[255];
      char xBu[255];
      char yBu[255];
      while(fscanf(fpM,"%s %s %s %s",cBu,vID,xBu,yBu)>0){
	SET_STRING_ELT(locusNames,ll,Rf_mkChar(vID));
	if(ll >= nLoci){
	  nLoci *= 2;
	  locusNames = Rf_lengthgets(locusNames, nLoci);
	}
	ll++;
      }
      nLoci = ll;
      locusNames = Rf_lengthgets(locusNames, nLoci);
    }


    bool progress = Rf_asLogical(progress0);
    
    bool hFID = Rf_asLogical(hasFID);
    int plFID = 0;
    int plIID = hFID ? 1 : 0;
    bool hPID = Rf_asLogical(hasParents);
    int plPatPID = plIID + (hPID ? 1 : 0);
    int plMatPID = plIID + (hPID ? 2 : 0);
    bool hSex = Rf_asLogical(hasSex);
    int plSex = plMatPID + (hSex ? 1 : 0);
    bool hPheno = Rf_asLogical(hasPheno);
    int plPheno = plSex + (hPheno ? 1 : 0);
    int plGeno = plPheno + 1;

    PROTECT_INDEX ipx_fam;
    SEXP famNames;
    PROTECT_WITH_INDEX(famNames = Rf_allocVector(STRSXP, nIndi),&ipx_fam);
    PROTECT_INDEX ipx_id;
    SEXP idNames;
    PROTECT_WITH_INDEX(idNames = Rf_allocVector(STRSXP, nIndi),&ipx_id);
    PROTECT_INDEX ipx_pat;
    SEXP patNames;
    PROTECT_WITH_INDEX(patNames = Rf_allocVector(STRSXP, nIndi),&ipx_pat);
    PROTECT_INDEX ipx_mat;
    SEXP matNames;
    PROTECT_WITH_INDEX(matNames = Rf_allocVector(STRSXP, nIndi),&ipx_mat);
    PROTECT_INDEX ipx_sex;
    SEXP sexNames;
    PROTECT_WITH_INDEX(sexNames = Rf_allocVector(STRSXP, nIndi),&ipx_sex);
    PROTECT_INDEX ipx_pheno;
    SEXP phenotypes;
    PROTECT_WITH_INDEX(phenotypes = Rf_allocVector(STRSXP, nIndi),&ipx_pheno);
    PROTECT_INDEX ipx_geno;
    SEXP genotypes;
    PROTECT_WITH_INDEX(genotypes = Rf_allocVector(VECSXP, nIndi),&ipx_geno);
    PROTECT_INDEX ipx_an;
    SEXP aNames;
    PROTECT_WITH_INDEX(aNames = Rf_allocVector(VECSXP, nLoci),&ipx_an);
    PROTECT_INDEX ipx_na;
    SEXP numA;
    PROTECT_WITH_INDEX(numA = Rf_allocVector(INTSXP, nLoci),&ipx_na);
    int* p_na = INTEGER(numA);
    for (int i = 0; i < nLoci; ++i) {
      p_na[i] = 0;
    }
    
    PROTECT_INDEX ipx_gt;
    SEXP gt, av, an;
    while (fgets(buff, sizeof buff, fp) != NULL) {
      int bufActiveSize = strcspn(buff, "\n");
      if(i_place == 0){
	if(progress)
	  Rprintf("Parsing individual %d\n",i_indi);
	PROTECT_WITH_INDEX(gt = Rf_allocVector(VECSXP, nLoci),&ipx_gt);  
      }
      // read until next space
      char pBuff[255];
      memset(pBuff,'\0',sizeof(pBuff));
      stream = fmemopen(buff, strlen(buff), "r");
      rewind(stream);
      while(fscanf(stream,"%s",pBuff)>0){
	int pbufActiveSize = strcspn(pBuff, "\n");
	//Rprintf("\tpBuff size: %d, %s\n",strcspn(pBuff, "\n"),pBuff);
	if(hFID && i_place == plFID){ // family name
	  SET_STRING_ELT(famNames,i_indi,Rf_mkChar(pBuff));
	}else if(i_place == plIID){ // individual name
	  SET_STRING_ELT(idNames,i_indi,Rf_mkChar(pBuff));
	}else if(hPID && i_place == plPatPID){ // Paternal name
	  SET_STRING_ELT(patNames,i_indi,Rf_mkChar(pBuff));
	}else if(hPID && i_place == plMatPID){ // Maternal name
	  SET_STRING_ELT(matNames,i_indi,Rf_mkChar(pBuff));
	}else if(hPID && i_place == plSex){ // Sex name
	  if(strcmp(pBuff,"1") == 0){
	    SET_STRING_ELT(sexNames,i_indi,Rf_mkChar("M"));
	  }else if(strcmp(pBuff,"2") == 0){
	    SET_STRING_ELT(sexNames,i_indi,Rf_mkChar("F"));
	  }else{
	    SET_STRING_ELT(sexNames,i_indi,R_NaString);
	  }
	}else if(hPID && i_place == plPheno){ // phenotype name
	  SET_STRING_ELT(phenotypes,i_indi,Rf_mkChar(pBuff));
	}else{ // Genotype
	  p_loci = (i_place - plGeno) / 2;
	  if(p_loci >= nLoci){
	    // if(progress)
	    //   Rprintf("Increasing number of loci from %d to %d\n",nLoci,2*nLoci);
	    // Grow things
	    nLoci *= 2;
	    REPROTECT(gt = Rf_lengthgets(gt, nLoci), ipx_gt);
	    REPROTECT(aNames = Rf_lengthgets(aNames, nLoci), ipx_an);
	    REPROTECT(numA = Rf_lengthgets(numA, nLoci), ipx_na);
	    p_na = INTEGER(numA);
	  }
	  int allele = (i_place - plGeno) % 2;
	  if(allele == 0){
	    PROTECT(av = Rf_allocVector(INTSXP, p_na[p_loci]));
	    for(int ii = 0; ii < Rf_length(av); ++ii)
	      INTEGER(av)[ii] = 0;
	    PROTECT(an = VECTOR_ELT(aNames, p_loci));
	    if(an == R_NilValue){
	      an = Rf_allocVector(STRSXP, p_na[p_loci]);
	    }
	  }
	  if(strcmp(pBuff,NA_char) != 0){ //Not missing
	    int *pg = INTEGER(av);
	    bool didIt = false;
	    if(Rf_length(an) > 0){
	      for(int ii = 0; ii < Rf_length(av); ++ii){	      
		if(strcmp(pBuff,R_CHAR(STRING_ELT(an,ii)))==0){
		  pg[ii]++;
		  didIt = true;
		  break;
		}
	      }
	    }
	    if(!didIt){
	      p_na[p_loci]++;
	      an = Rf_lengthgets(an, p_na[p_loci]);
	      SET_STRING_ELT(an,p_na[p_loci]-1,Rf_mkChar(pBuff));
	      av = Rf_lengthgets(av, p_na[p_loci]);	      
	      INTEGER(av)[p_na[p_loci]-1] = 1;
	    }	      	    
	  }
	  if(allele == 1){
	    Rf_setAttrib(av,Rf_install("names"),an);
	    SET_VECTOR_ELT(gt,p_loci,av);
	    SET_VECTOR_ELT(aNames,p_loci,an);
	    if(i_indi > 0){
	      int l_last = Rf_length(VECTOR_ELT(VECTOR_ELT(genotypes,i_indi-1),p_loci));
	      if(l_last < Rf_length(av)){
		//Rprintf("Need to correct back");
		// Correct for all individuals going back
		for(int j_indi = 0; j_indi < i_indi; ++j_indi){
		  SEXP gtj = PROTECT(VECTOR_ELT(genotypes,j_indi));
		  SEXP av0 = PROTECT(VECTOR_ELT(gtj,p_loci));
		  int nav0 = Rf_length(av0);
		  av0 = Rf_lengthgets(av0, p_na[p_loci]);
		  for(int qq = nav0; qq < p_na[p_loci]; ++qq)
		    INTEGER(av0)[qq] = 0;
		  Rf_setAttrib(av0,Rf_install("names"),an);
		  SET_VECTOR_ELT(gtj,p_loci,av0);
		  SET_VECTOR_ELT(genotypes,j_indi,gtj);
		  UNPROTECT(2);
		}
	      }
	    }
	    UNPROTECT(2);
	  }	 
	}	  
	if(pbufActiveSize < (int)(sizeof pBuff)-1 || (isspace(pBuff[(sizeof pBuff) -2]) && pBuff[(sizeof pBuff) -2]!='\n') || pBuff[(sizeof pBuff) -2] == '\0'){ // Reading next place
	  i_place++;
	}else{
	  // if(progress)
	  //Rprintf("%d %d %s\n",i_indi, i_place, pBuff);	
	  if(i_place < plGeno){
	    Rf_error("Error reading input. IDs are probably too long.");
	  }else{
	    Rf_warning("There may be an issue with very long allele names!");
	  }
	}
	R_CheckUserInterrupt();
      }
      fclose(stream);
      if(bufActiveSize < (int)(sizeof buff)-1 || buff[(sizeof buff) -2]=='\n'){//Starting new individual next time
	// Save genotype
	if(i_indi == 0){
	  //Rprintf("Resizing loci from %d to %d\n",nLoci,p_loci + 1);
	  nLoci = p_loci + 1;
	  REPROTECT(gt = Rf_lengthgets(gt, nLoci), ipx_gt);
	}
	if(locusNames != R_NilValue)
	  Rf_setAttrib(gt,R_NamesSymbol,locusNames);
	Rf_setAttrib(gt,Rf_install("ploidy"),Rf_ScalarInteger(2));
	Rf_setAttrib(gt,R_ClassSymbol,Rf_mkString("Genotype"));
	SET_VECTOR_ELT(genotypes,i_indi,gt);
	UNPROTECT(1);
	// Prepare for next time
	i_indi++;
	i_place = 0;
	if(i_indi >= nIndi){
	  //Rprintf("Increasing number of individuals from %d to %d\n",nIndi,2*nIndi);
	  // Grow things
	  nIndi *= 2;
	  REPROTECT(famNames = Rf_lengthgets(famNames, nIndi), ipx_fam);
	  REPROTECT(idNames = Rf_lengthgets(idNames, nIndi), ipx_id);
	  REPROTECT(patNames = Rf_lengthgets(patNames, nIndi), ipx_pat);
	  REPROTECT(matNames = Rf_lengthgets(matNames, nIndi), ipx_mat);
	  REPROTECT(sexNames = Rf_lengthgets(sexNames, nIndi), ipx_sex);
	  REPROTECT(phenotypes = Rf_lengthgets(phenotypes, nIndi), ipx_pheno);
	  REPROTECT(genotypes = Rf_lengthgets(genotypes, nIndi), ipx_geno);
	}
      }
    }
  
 
    // Resize to actual number of individuals
    //Rprintf("Resizing output from %d to %d\n",nIndi,i_indi);
    nIndi = i_indi;
    REPROTECT(famNames = Rf_lengthgets(famNames, nIndi), ipx_fam);
    REPROTECT(idNames = Rf_lengthgets(idNames, nIndi), ipx_id);
    REPROTECT(patNames = Rf_lengthgets(patNames, nIndi), ipx_pat);
    REPROTECT(matNames = Rf_lengthgets(matNames, nIndi), ipx_mat);
    REPROTECT(sexNames = Rf_lengthgets(sexNames, nIndi), ipx_sex);
    REPROTECT(phenotypes = Rf_lengthgets(phenotypes, nIndi), ipx_pheno);
    REPROTECT(genotypes = Rf_lengthgets(genotypes, nIndi), ipx_geno);      
    
    Rf_setAttrib(genotypes,R_NamesSymbol,idNames);
    Rf_setAttrib(genotypes,Rf_install("population"),famNames);
    Rf_setAttrib(genotypes,Rf_install("paternalID"),patNames);
    Rf_setAttrib(genotypes,Rf_install("maternalID"),matNames);
    Rf_setAttrib(genotypes,Rf_install("sex"),sexNames);
    Rf_setAttrib(genotypes,Rf_install("phenotype"),phenotypes);
    Rf_setAttrib(genotypes,R_ClassSymbol,Rf_mkString("gen"));
    
    fclose(fp);
    UNPROTECT(10);
    return genotypes;
  }

}

