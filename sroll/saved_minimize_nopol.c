void minimize_gain_nopol(double *ix2,double *gaingi){

  MPI_Status statu;
  long i,rrk,j,k,l1,l2,ib;
  int itermax = (NUMBEROFITER);
  int iter;
  double  delta_new, delta_old, beta;
  double  alpha=1.0; // get rid of gcc "maybe-uninitialized" warning depending on optimsation level

  PIOLONG GAINSTEP2;
  int rank;
  int size;
  int mpi_size;
  int rank_size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank_size);
  rank=rank_size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  mpi_size=size;

  GAINSTEP2=GAINSTEP;
    
  nmatres=newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+nfreefree;

  if (rank==0) {
    fprintf(stderr,"==============================\n\nminimize_gain_nopol(): ITERATION NOPOL %ld/%d %ld \n\n==============================\n",itbogo,Param->NITT,(long) GAINSTEP);
    fprintf(stderr,"GAIN ");
    for (i=0;i<nbolo;i++) fprintf(stderr,"%lg ",gaingi[i*GAINSTEP]);
    fprintf(stderr,"\n");
    if (GAINSTEP>1) {
      fprintf(stderr,"...\nGAIN ");
      for (i=0;i<nbolo;i++) fprintf(stderr,"%lg ",gaingi[i*GAINSTEP+(GAINSTEP-1)]);
      fprintf(stderr,"\n");
    }
  }

  if (itbogo==0) delta0=0;
  MPI_Barrier(MPI_COMM_WORLD); //Synchronization MPI process

  struct timeval tp1,tp2;
  gettimeofday(&tp1,NULL);


  iter = 0;
  memset(b2  ,0,nmatres*sizeof (double));
  memset(d2  ,0,nmatres*sizeof (double));
  memset(q2  ,0,nmatres*sizeof (double));
  memset(r2  ,0,nmatres*sizeof (double));
  memset(s2  ,0,nmatres*sizeof (double));
  memset(hit2,0,nmatres*sizeof (double));

  //===============================================================================================================================================================================
  //=  Compute second member
  //=
  //===============================================================================================================================================================================
  long l,m;

  ////// BUILD B2
  //GetProcMem(&vmem,&phymem);
  //if (rank==0) fprintf(stderr,"Rank: %ld Line=%d MEM %.1lf[%.1lf]MB\n",
  //                  (long) rank, __LINE__,
  //                  (double) vmem/1024./1024.,
  //                  (double) phymem/1024./1024.);

  
  int ptest=0;
  for (k=0;k<nnbpix;k++)  {                   // Pour chaque pixel 
    //long imat=the_stat_pix[k];
    long ndata = loc_nhpix[k];
    hpix *htmp = loc_hpix[k];

    II[k]=0;
     
    for (l1=0;l1<ndata;l1++) {
      long ri1=htmp[l1].rg-globalBeginRing;
      if (flg_rg[htmp[l1].ib][ri1]!=0) {
        II[k]+=htmp[l1].w;
      }
    }


    if (II[k]==0&&flgpix[k]>0) fprintf(stderr,"FLAG PROBLEM\n");

    if (ndata>0&&flgpix[k]>0) {

      double SI=0;

      #ifdef USEDII
      memset(dii ,0,sizeof(double)*nbolo*GAINSTEP);
      #endif
      memset(dcoi+k*nbolo,0,sizeof(double)*nbolo);
      memset(dfri+k*nbolo,0,sizeof(double)*nbolo);
      memset(dthetai+k*nbolo,0,sizeof(double)*nbolo);
      memset(ddusti+k*nbolo,0,sizeof(double)*nbolo);
      memset(dpixi+k*nbolo*npixbeam,0,sizeof(double)*nbolo*npixbeam);


      for (l1=0;l1<ndata;l1++) {
        long ri1=htmp[l1].rg-globalBeginRing;
        double g1=gaingi[htmp[l1].gi+htmp[l1].ib*GAINSTEP];

        if (flg_rg[htmp[l1].ib][ri1]!=0) {
          htmp[l1].wp=0;
          //htmp[l1].thsig=0;

          if (REMHDIP==0) SI+=htmp[l1].w*(htmp[l1].sig*g1-htmp[l1].fsl-htmp[l1].dip-htmp[l1].corr_nl-htmp[l1].corr_cnn);
          else SI+=htmp[l1].w*(htmp[l1].sig*g1-htmp[l1].fsl-htmp[l1].freefree-htmp[l1].corr_nl-htmp[l1].corr_cnn);

          #ifdef USEDII
          dii[htmp[l1].gi+htmp[l1].ib*GAINSTEP ]+=htmp[l1].w*htmp[l1].dip;
          #endif

          if (nmatco>0) {
            dcoi[htmp[l1].ib+k*nbolo] += htmp[l1].w*htmp[l1].comap;
          }
          if (nmatdust>0) {
            ddusti[htmp[l1].ib+k*nbolo] += htmp[l1].w*htmp[l1].dustmap;
          }
          if (nfreefree>0) {
            dfri[htmp[l1].ib+k*nbolo] += htmp[l1].w*htmp[l1].freefree;
          }

          if (ittt>0) {
            for (m=0;m<npixbeam;m++)  {
	            dpixi[nbolo*m+htmp[l1].ib+k*nbolo*npixbeam] += htmp[l1].w*htmp[l1].listofpix[m];
	            if (isnan(dpixi[nbolo*m+htmp[l1].ib+k*nbolo*npixbeam])) {
		            fprintf(stderr,"v %d %ld %lf\n",__LINE__,m,htmp[l1].listofpix[m]);
		            exit(0);
	            }
            }
          }
        }
      }

      SI/=II[k];
      SSI[k]=SI;


      for (ib=0;ib<nbolo;ib++) {
        #ifdef USEDII
        for (j=0;j<GAINSTEP;j++) {
          dii[j+ib*GAINSTEP]=dii[j+ib*GAINSTEP]/II[k];
        }
        #endif

        if (nmatco>0) {
          dcoi[ib+k*nbolo]=dcoi[ib+k*nbolo]/II[k];
        }
        if (nmatdust>0) {
          ddusti[ib+k*nbolo]=ddusti[ib+k*nbolo]/II[k];
        }
        if (nfreefree>0) {
          dfri[ib+k*nbolo]=dfri[ib+k*nbolo]/II[k];
        }

        if (ittt>0) {
          for (m=0;m<npixbeam;m++)  {
      	    dpixi[ib+m*nbolo+k*nbolo*npixbeam]=dpixi[ib+m*nbolo+k*nbolo*npixbeam]/II[k];
          }
        }
      }


      for (l1=0;l1<ndata;l1++) {
	      htmp[l1].vi=UNSEENPIX;
        long ri1=htmp[l1].rg-globalBeginRing;
        if (flg_rg[htmp[l1].ib][ri1]!=0) {
          htmp[l1].vi=htmp[l1].w/II[k];
          #ifndef USEDII
          htmp[l1].lvi=htmp[l1].w*htmp[l1].dip/II[k];
          #endif
        }
      }

      for (l1=0;l1<ndata;l1++) {
        long ri1=htmp[l1].rg-globalBeginRing;
        double g1=gaingi[htmp[l1].gi+htmp[l1].ib*GAINSTEP];

        if (flg_rg[htmp[l1].ib][ri1]!=0) {
          double divi=NEP_tab[htmp[l1].ib]*htmp[l1].hit*g1;
          divi=divi*divi;

          if (divi==0) htmp[l1].wp=0;
          else htmp[l1].wp=1/divi;

          if (itbogo==0) normaoff+=NEP_tab[htmp[l1].ib]*htmp[l1].hit;

        }
      }



      for (l1=0;l1<ndata;l1++) {
        long ri1=htmp[l1].rg-globalBeginRing;
        long iri1=rgord[htmp[l1].ib][ri1]+newnr[htmp[l1].ib];

        double g1=gaingi[htmp[l1].gi+htmp[l1].ib*GAINSTEP];

        if (flg_rg[htmp[l1].ib][ri1]!=0) {
          double ww=htmp[l1].wp;
          double tmp=((htmp[l1].sig*g1-htmp[l1].fsl-htmp[l1].corr_nl-htmp[l1].corr_cnn)-(SI));
          if (REMHDIP==0) tmp-=htmp[l1].dip;
          else tmp-=htmp[l1].freefree;

          long l3;
          #ifndef USEDII
          memset(cdip,0,GAINSTEP*nbolo*sizeof(double));
          #endif

          for (l3=0;l3<ndata;l3++) {
            long ri3=htmp[l3].rg-globalBeginRing;
            if (flg_rg[htmp[l3].ib][ri3]!=0) {
              long ir3=rgord[htmp[l3].ib][ri3]+newnr[htmp[l3].ib];
              ctmp[ir3]=-htmp[l3].vi;
              #ifndef USEDII
              cdip[htmp[l3].ib*GAINSTEP+htmp[l3].gi]-=htmp[l3].lvi;
              #endif
            }
          }
          ctmp[iri1]+=1;

          #ifdef USEDII
          for (ib=0;ib<nbolo;ib++){
            for (j=0;j<GAINSTEP;j++) {
              cdip[ib*GAINSTEP+j]=-dii[j+ib*GAINSTEP];
            }
          }
          #endif
    
          cdip[htmp[l1].ib*GAINSTEP+htmp[l1].gi]+=htmp[l1].dip; 

          for (ib=0;ib<nbolo;ib++) cco[ib]=-dcoi[ib+k*nbolo];
          cco[htmp[l1].ib]+=htmp[l1].comap;

          for (ib=0;ib<nbolo;ib++) cdust[ib]=-ddusti[ib+k*nbolo];
          cdust[htmp[l1].ib]+=htmp[l1].dustmap;

          for (ib=0;ib<nbolo;ib++) ccfree[ib]=-dfri[ib+k*nbolo];
          ccfree[htmp[l1].ib]+=htmp[l1].freefree;

          for (m=0;m<npixbeam;m++) {
	          for (ib=0;ib<nbolo;ib++) cpix[ib+m*nbolo]=-dpixi[ib+m*nbolo+k*nbolo*npixbeam];
	          cpix[htmp[l1].ib+m*nbolo]+=htmp[l1].listofpix[m];
          }

          long ir;
          /////////////////  OFFSET

          for (l2=0;l2<ndata;l2++) {
            long ri2=htmp[l2].rg-globalBeginRing;
            if (flg_rg[htmp[l2].ib][ri2]!=0) {
              long ir=rgord[htmp[l2].ib][ri2]+newnr[htmp[l2].ib];
              b2[ir]+=ww*tmp*ctmp[ir];
              hit2[ir]+=ww*ctmp[ir]*ctmp[ir];
#ifndef USEDII
              b2[newnr[nbolo]+htmp[l2].ib*GAINSTEP+htmp[l2].gi]+=ww*tmp*cdip[htmp[l2].ib*GAINSTEP+htmp[l2].gi];
              hit2[newnr[nbolo]+htmp[l2].ib*GAINSTEP+htmp[l2].gi]+=ww*cdip[htmp[l2].ib*GAINSTEP+htmp[l2].gi]
                *cdip[htmp[l2].ib*GAINSTEP+htmp[l2].gi];
#endif
            }
          }


          /////////////////  DIPOLE FIT

          for (ir=0;ir<nbolo;ir++) {

#ifdef USEDII

            for (j=0;j<GAINSTEP;j++) {
              b2[newnr[nbolo]+ir*GAINSTEP+j]+=ww*tmp*cdip[ir*GAINSTEP+j];
              hit2[newnr[nbolo]+ir*GAINSTEP+j]+=ww*cdip[ir*GAINSTEP+j]*cdip[ir*GAINSTEP+j];
            }

#endif

            ///////////// CO
            if (nmatco>0) {
              b2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+ir]+=ww*tmp*cco[ir];
              hit2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+ir]+=ww*cco[ir]*cco[ir];
            }

            ///////////// DUST

            if (nmatdust>0) {
              b2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+ir]+=ww*tmp*cdust[ir];
              hit2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+ir]+=ww*cdust[ir]*cdust[ir];
            }

            ///////////// FREEFREE

            if (nfreefree>0) {
              b2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+ir]+=ww*tmp*ccfree[ir];
              hit2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+ir]+=ww*ccfree[ir]*ccfree[ir];
            }

            ////////// PIXBEAM
            for (j=0;j<npixbeam;j++) {
              if (NOMOREFITTED!=DOCNN[j]||DOCNN[j]==0) {
	              b2[newnr[nbolo]+nbolo*GAINSTEP2+j*nbolo+ir]+=ww*tmp*cpix[ir+j*nbolo];
              }
              hit2[newnr[nbolo]+nbolo*GAINSTEP2+j*nbolo+ir]+=ww*cpix[ir+j*nbolo]*cpix[ir+j*nbolo];
            }
          }
        }
      }

    }
  }

  #ifdef OPTIMPI
  {
    double *lb = (double *) malloc(sizeof(double)*(nmatres));
    MPI_Reduce(b2,lb,nmatres,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

    memcpy(b2,lb,sizeof(double)*(nmatres));
    free(lb);
  }
  #else
  if (rank==0) {
    double *lb = (double *) malloc(sizeof(double)*(nmatres));
    for (rrk=1;rrk<mpi_size;rrk++) {
      MPI_Recv(lb,sizeof(double)*(nmatres), MPI_BYTE, rrk,1031, MPI_COMM_WORLD,&statu);
      for (l=0;l<nmatres;l++) b2[l]+=lb[l];
    }
    free(lb);
  }
  else MPI_Send(b2, sizeof(double)*(nmatres), MPI_BYTE, 0, 1031, MPI_COMM_WORLD);
  #endif
  #ifdef OPTIMPI
  {
    double *lb = (double *) malloc(sizeof(double)*(nmatres));
    MPI_Reduce(hit2,lb,nmatres,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

    memcpy(hit2,lb,sizeof(double)*(nmatres));
    free(lb);
  }
  #else
  if (rank==0) {
    double *lb = (double *) malloc(sizeof(double)*(nmatres));
    for (rrk=1;rrk<mpi_size;rrk++) {
      MPI_Recv(lb,sizeof(double)*(nmatres), MPI_BYTE, rrk,1033, MPI_COMM_WORLD,&statu);

      for (l=0;l<nmatres;l++) hit2[l]+=lb[l];
    }
    free(lb);
  }
  else {
    MPI_Send(hit2, sizeof(double)*(nmatres), MPI_BYTE, 0, 1033, MPI_COMM_WORLD);
  }
  #endif


  //==============================================================================================================================================================================
  // Compute Ax
  //healspline_fill_s_with_PtP_s (hs_rec, Ndata, vec, x, q);
  // +
  // Preconditionnement
  //
  // ==============================================================================================================================================================================

  for (l=0;l<nmatres;l++) q2[l]=0;
  if (rank==0) {

    double soff=0;
    for (i=0;i<newnr[nbolo];i++) soff+=hit2[0]*ix2[i];
    for (i=0;i<newnr[nbolo];i++) q2[i]=soff;

    if (Param->REMHDIP==1) {
      soff=0;
      for (i=newnr[nbolo];i<newnr[nbolo]+GAINSTEP*nbolo;i++) soff+=hit2[newnr[nbolo]]*ix2[i];
      for (i=newnr[nbolo];i<newnr[nbolo]+GAINSTEP*nbolo;i++) q2[i]=soff;
    }

    if (Param->flag_AVGR0==_PAR_TRUE) {
      soff=0;
      for (i=0;i<nbolo;i++) soff+=hit2[newnr[nbolo]+nbolo*(GAINSTEP2)]*
                              ix2[newnr[nbolo]+nbolo*(GAINSTEP2)+i];
      for (i=0;i<nbolo;i++) q2[newnr[nbolo]+nbolo*(GAINSTEP2)+i]=soff;
      for (i=0;i<nbolo;i++) b2[newnr[nbolo]+nbolo*(GAINSTEP2)+i]+=
                              hit2[newnr[nbolo]+nbolo*(GAINSTEP2)]*Param->AVGR0;

    }

    if (nmatco>0) {
      soff=0;
      if ((singleFreq == 0) && (Param->flag_AVG12CO == _PAR_TRUE)) {
        for (i=0;i<nbolo;i++) if (freqs[i] == Param->AVFCO)
            soff += hit2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)]
                    * ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+i];
        for (i=0;i<nbolo;i++) if (freqs[i] == Param->AVFCO)
            q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+i] = soff;
        for (i=0;i<nbolo;i++) if (freqs[i] == Param->AVFCO)
            b2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+i] +=
            hit2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)] * Param->AVG12CO;
      }
      else {
        for (i=0;i<nbolo;i++) soff+=hit2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)]
                                *ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+i];
        for (i=0;i<nbolo;i++) q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+i]=soff;
        for (i=0;i<nbolo;i++) b2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+i]+=
                                hit2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)]*Param->AVG12CO;
      }
    }
    if (nmatdust>0) {
      soff=0;
      if ((singleFreq == 0) && (Param->flag_AVGDUST100 == _PAR_TRUE)) {
        for (i=0;i<nbolo;i++) if (freqs[i] == Param->AVFDUST)
            soff += hit2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco]
                    * ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+i];
        for (i=0;i<nbolo;i++) if (freqs[i] == Param->AVFDUST)
            q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+i] = soff;
        for (i=0;i<nbolo;i++) if (freqs[i] == Param->AVFDUST)
            b2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+i] +=
            hit2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco] * Param->AVGDUST100;
      }
      else {
        for (i=0;i<nbolo;i++)
          soff+=hit2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco]*
            ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+i];
        for (i=0;i<nbolo;i++) q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+i]=soff;
        for (i=0;i<nbolo;i++) b2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+i]+=
                                hit2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco]*Param->AVGDUST;
      }
    }
    if (nfreefree>0) {
      soff=0;
      if ((singleFreq == 0) && (Param->flag_AVGFREEFREE == _PAR_TRUE)) {
	for (i=0;i<nbolo;i++) if (freqs[i] == Param->AVFFF)
	  soff+=hit2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust]*
	    ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+i];
	for (i=0;i<nbolo;i++) if (freqs[i] == Param->AVFFF)
				q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+i]=soff;
	for (i=0;i<nbolo;i++) if (freqs[i] == Param->AVFFF)
				b2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+i]+=
                                hit2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust]*Param->AVGFREEFREE;
      }
      else {
	for (i=0;i<nbolo;i++)
	  soff+=hit2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust]*
	    ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+i];
	for (i=0;i<nbolo;i++) q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+i]=soff;
      }
    }
    long nj=0;
    if (DOFITANGLE==1) nj++;
    if (DOFITPOLEFF==1) nj++;
    if (DOTDUST==1) nj++;
    if (DOCO13==1) nj++;
    if (DOSYNCHRO==1) nj++;


    if (NORMFITPOL==1) {
#if 1
      if (DOFITPOLEFF==1) {
	j=1+DOCO13+DOSYNCHRO+DOTDUST;
	soff=0;
	for (i=0;i<nbolo;i++) soff+=hit2[newnr[nbolo]+nbolo*(GAINSTEP2+npixbeam-j)]*
				ix2[newnr[nbolo]+nbolo*(GAINSTEP2+npixbeam-j)+i];
	for (i=0;i<nbolo;i++) q2[newnr[nbolo]+nbolo*(GAINSTEP2+npixbeam-j)+i]=soff;
      }
#endif
    }
    
    if (DOFITANGLE==1) {
      j=1+DOCO13+DOSYNCHRO+DOTDUST+DOFITPOLEFF;
      soff=0;
      for (i=0;i<nbolo;i++) soff+=hit2[newnr[nbolo]+nbolo*(GAINSTEP2+npixbeam-j)]*
			      ix2[newnr[nbolo]+nbolo*(GAINSTEP2+npixbeam-j)+i];
      for (i=0;i<nbolo;i++) q2[newnr[nbolo]+nbolo*(GAINSTEP2+npixbeam-j)+i]=soff;
    }

    if (DOCO13==1) {
      j=1+DOSYNCHRO;
      soff=0;
      for (i=0;i<nbolo;i++) if (freqs[i] == Param->AVFCO) soff+=hit2[newnr[nbolo]+nbolo*(GAINSTEP2+npixbeam-j)]*
						   ix2[newnr[nbolo]+nbolo*(GAINSTEP2+npixbeam-j)+i];
      for (i=0;i<nbolo;i++) if (freqs[i] == Param->AVFCO) q2[newnr[nbolo]+nbolo*(GAINSTEP2+npixbeam-j)+i]=soff;
      for (i=0;i<nbolo;i++) if (freqs[i] == Param->AVFCO) b2[newnr[nbolo]+nbolo*(GAINSTEP2+npixbeam-j)+i]+=
							    hit2[newnr[nbolo]+nbolo*(GAINSTEP2+npixbeam-j)] * Param->AVG13CO;
    }

    if (DOTDUST==1) {
      j=1+DOCO13+DOSYNCHRO;
      soff=0;
      for (i=0;i<nbolo;i++) if (freqs[i] == MAXFREQ) soff+=hit2[newnr[nbolo]+nbolo*(GAINSTEP2+npixbeam-j)]*
						       ix2[newnr[nbolo]+nbolo*(GAINSTEP2+npixbeam-j)+i];
      for (i=0;i<nbolo;i++) if (freqs[i] == MAXFREQ) q2[newnr[nbolo]+nbolo*(GAINSTEP2+npixbeam-j)+i]=soff;
    }

    if (DOSYNCHRO==1) {
      soff=0;
      j=1;
      for (i=0;i<nbolo;i++) if (freqs[i] == Param->AVFSYNC) soff+=hit2[newnr[nbolo]+nbolo*(GAINSTEP2+npixbeam-j)]*
						 ix2[newnr[nbolo]+nbolo*(GAINSTEP2+npixbeam-j)+i];
      for (i=0;i<nbolo;i++) if (freqs[i] == Param->AVFSYNC) q2[newnr[nbolo]+nbolo*(GAINSTEP2+npixbeam-j)+i]=soff;
      for (i=0;i<nbolo;i++) if (freqs[i] == Param->AVFSYNC)
			      b2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2-j)+i]+=
                                hit2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2-j)]*Param->AVGSYNCHRO;

    }
  }

  for (k=0;k<nnbpix;k++) if (flgpix[k]>0) {

    long ndata = loc_nhpix[k];
    if (ndata>0) {
      #ifdef TIMING
      gettimeofday(&tp1,NULL);
      #endif
      hpix *htmp = loc_hpix[k];

      double vali=0;

      #ifdef USEDII
      memset(dii ,0,sizeof(double)*nbolo*GAINSTEP);

      for (l1=0;l1<ndata;l1++) {
        long ri1=htmp[l1].rg-globalBeginRing;

        if (flg_rg[htmp[l1].ib][ri1]!=0) {

          dii[htmp[l1].gi+htmp[l1].ib*GAINSTEP ]+=htmp[l1].w*htmp[l1].dip;

        }
      }
      for (ib=0;ib<nbolo;ib++) {
        for (j=0;j<GAINSTEP;j++) {
          if(II[k]==0){
            fprintf(stderr,"nbolo = %d II[k] = %f \n",htmp[l1].ib,II[k]);
          }

          dii[j+ib*GAINSTEP]=dii[j+ib*GAINSTEP]/II[k];
        }
      }
      #endif

      long l3;
      for (l3=0;l3<ndata;l3++) {
        long ri3=htmp[l3].rg-globalBeginRing;
        if (flg_rg[htmp[l3].ib][ri3]!=0) {
          long ir3=rgord[htmp[l3].ib][ri3]+newnr[htmp[l3].ib];
          vali+=htmp[l3].vi*ix2[ir3];
	  if (isnan(vali))  {
	    fprintf(stderr,"v %d %lf\n",__LINE__,htmp[l3].vi);
	    exit(0);
	  }

#ifndef USEDII
          vali+=htmp[l3].lvi*ix2[newnr[nbolo]+htmp[l3].ib*GAINSTEP+htmp[l3].gi];
	  if (isnan(vali))  {
	    fprintf(stderr,"v %d %lf\n",__LINE__,htmp[l3].lvi);
	    exit(0);
	  }
#endif
        }
	}
      //fprintf(stderr,"I0 vali %lg %lg %lg\n",vali,valq,valu);

	  for (ib=0;ib<nbolo;ib++) {
#ifdef USEDII
	  for (j=0;j<GAINSTEP;j++) {
          vali+=dii[j+ib*GAINSTEP]  *ix2[newnr[nbolo]+ib*GAINSTEP+j];

	  if (isnan(vali))  {
	  fprintf(stderr,"v %d %lf %lf \n",__LINE__,dii[j+ib*GAINSTEP] ,ix2[newnr[nbolo]+ib*GAINSTEP+j]);
	  exit(0);
	}
	}
        //fprintf(stderr,"I1 vali %lg %lg %lg\n",vali,valq,valu);
#endif

        if (nmatco>0) {
          vali+=dcoi[ib+k*nbolo]  *ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+ib];
        }

	  if (isnan(vali))  {
	  fprintf(stderr,"v %d %lf \n",__LINE__,dcoi[ib+k*nbolo] );
	  exit(0);
	}
        if (nmatdust>0) {
          vali+=ddusti[ib+k*nbolo]  *ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+ib];
        }
	  if (isnan(vali))  {
	  fprintf(stderr,"v %d %lf \n",__LINE__,ddusti[ib+k*nbolo] );
	  exit(0);
	}
        if (nfreefree>0) {
          vali+=dfri[ib+k*nbolo]  *ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+ib];
        }
	  if (isnan(vali))  {
	  fprintf(stderr,"v %d %lf \n",__LINE__,dfri[ib+k*nbolo] );
	  exit(0);
	}
        //fprintf(stderr,"I3 vali %lg %lg %lg : %lg %lg %lg %lg\n",vali,valq,valu,
        //dfri[ib+k*nbolo],dfrq[ib+k*nbolo],dfru[ib+k*nbolo],
        //ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+ib]);

        for (m=0;m<npixbeam;m++) {
	  if (NOMOREFITTED!=DOCNN[m]||DOCNN[m]==0) {
	    vali+=ix2[newnr[nbolo]+nbolo*GAINSTEP2+ib+m*nbolo]*dpixi[ib+m*nbolo+k*nbolo*npixbeam];
	    if (isnan(vali))  {
	      fprintf(stderr,"v %d %lf \n",__LINE__,dpixi[ib+m*nbolo+k*nbolo*npixbeam] );
	      exit(0);
	    }
	  }
        }

        //fprintf(stderr,"I4 vali %lg %lg %lg\n",vali,valq,valu);
	}
      double qri=0;

      for (l1=0;l1<ndata;l1++) {
        long ri1=htmp[l1].rg-globalBeginRing;
        long iri1=rgord[htmp[l1].ib][ri1]+newnr[htmp[l1].ib];

        if (flg_rg[htmp[l1].ib][ri1]!=0) {
          double ww=htmp[l1].wp;
          double val2=ix2[iri1]-(vali);

          val2+=ix2[newnr[nbolo]+htmp[l1].ib*GAINSTEP+htmp[l1].gi]*htmp[l1].dip;
	  if (isnan(val2))  {
	    fprintf(stderr,"v %d %lf %lf\n",__LINE__,ix2[newnr[nbolo]+htmp[l1].ib*GAINSTEP+htmp[l1].gi],htmp[l1].dip);
	    exit(0);
	  }

          //fprintf(stderr,"I2 val2 %lg\n",val2);
          if (nmatco>0)
            val2+=ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+htmp[l1].ib]*htmp[l1].comap;
	  if (isnan(val2))  {
	    fprintf(stderr,"v %d %lf %lf\n",__LINE__,ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+htmp[l1].ib],htmp[l1].comap);
	    exit(0);
	  }
          if (nmatdust>0)
            val2+=ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+htmp[l1].ib]*htmp[l1].dustmap;
	  if (isnan(val2))  {
	    fprintf(stderr,"v %d %lf %lf\n",__LINE__,ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+htmp[l1].ib],htmp[l1].dustmap);
	    exit(0);
	  }
          if (nfreefree>0)
            val2+=ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+htmp[l1].ib]*htmp[l1].freefree;

          //fprintf(stderr,"I3 val2 %lg\n",val2);
          for (m=0;m<npixbeam;m++) {
	    if (NOMOREFITTED!=DOCNN[m]||DOCNN[m]==0) {
	      val2+=ix2[newnr[nbolo]+nbolo*GAINSTEP2+htmp[l1].ib+m*nbolo]*htmp[l1].listofpix[m];
	      if (isnan(val2))  {
		fprintf(stderr,"v %d %lf %lf\n",__LINE__,ix2[newnr[nbolo]+nbolo*GAINSTEP2+htmp[l1].ib+m*nbolo],htmp[l1].listofpix[m]);
		exit(0);
	      }
	    }

          }

          qri-=ww*val2;

          q2[iri1]+=ww*val2;
	  if (isnan(q2[iri1]))  {
	    fprintf(stderr,"Q2NAN %lf %lf %lf\n",ww,val2,q2[iri1]);
	    exit(0);
	  }
          q2[newnr[nbolo]+htmp[l1].ib*GAINSTEP+htmp[l1].gi]+=ww*val2*htmp[l1].dip;

          ///////////// CO
          if (nmatco>0) {
            q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+htmp[l1].ib]+=ww*val2*htmp[l1].comap;
          }

          ///////////// DUST
          if (nmatdust>0) {
            q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+htmp[l1].ib]+=ww*val2*htmp[l1].dustmap;
          }

          ///////////// FREEFREE
          if (nfreefree>0) {
            q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+htmp[l1].ib]+=ww*val2*htmp[l1].freefree;
          }

          ////////// SYSTE
          for (j=0;j<npixbeam;j++) {
	    if (NOMOREFITTED!=DOCNN[j]||DOCNN[j]==0) {
	      q2[newnr[nbolo]+nbolo*GAINSTEP2+j*nbolo+htmp[l1].ib]+=ww*val2*htmp[l1].listofpix[j];
	    }
          }


        }
      }


      for (l2=0;l2<ndata;l2++) {
        long ri2=htmp[l2].rg-globalBeginRing;
        if (flg_rg[htmp[l2].ib][ri2]!=0) {
          long ir=rgord[htmp[l2].ib][ri2]+newnr[htmp[l2].ib];
          q2[ir]+=qri*htmp[l2].vi;
#ifndef USEDII
          q2[newnr[nbolo]+htmp[l2].ib*GAINSTEP+htmp[l2].gi]+=qri*htmp[l2].lvi;
#endif
        }
      }


      for (ib=0;ib<nbolo;ib++) {

#ifdef USEDII
        for (j=0;j<GAINSTEP;j++) {
          q2[newnr[nbolo]+GAINSTEP*ib+j]+=qri*dii[j+ib*GAINSTEP];
        }
#endif

        ///////////// CO
        if (nmatco>0) {
          q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+ib]+=qri*dcoi[ib+k*nbolo];
        }
        ///////////// DUST
        if (nmatdust>0) {
          q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+ib]+=qri*ddusti[ib+k*nbolo];
        }
        ///////////// FREEFREE
        if (nfreefree>0) {
          q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+ib]+=qri*dfri[ib+k*nbolo];
        }

        for (j=0;j<npixbeam;j++) {
	  if (NOMOREFITTED!=DOCNN[j]||DOCNN[j]==0) {
	    q2[newnr[nbolo]+nbolo*GAINSTEP2+j*nbolo+ib]+=qri*dpixi[ib+j*nbolo+k*nbolo*npixbeam];
	  }
        }
      }
    }

#ifdef TIMING
      gettimeofday(&tp2,NULL);
      double dt=(double)(tp2.tv_sec-tp1.tv_sec)+(1E-6)*(tp2.tv_usec-tp1.tv_usec);
      dthit[ndata]+=dt;
      ndthit[ndata]+=1;
#endif
  }

#ifdef OPTIMPI
  {
    double *lb = (double *) malloc(sizeof(double)*(nmatres));
    MPI_Reduce(q2,lb,nmatres,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

    memcpy(q2,lb,sizeof(double)*(nmatres));
    free(lb);
  }
#else

  if (rank==0) {
    double *lb = (double *) malloc(sizeof(double)*(nmatres));
    for (rrk=1;rrk<mpi_size;rrk++) {
      MPI_Recv(lb,sizeof(double)*(nmatres), MPI_BYTE, rrk,1032, MPI_COMM_WORLD,&statu);
      for (l=0;l<nmatres;l++) q2[l]+=lb[l];
    }
    free(lb);
  }
  else {
    MPI_Send(q2, sizeof(double)*(nmatres), MPI_BYTE, 0, 1032, MPI_COMM_WORLD);
  }
#endif

  if (rank==0) fprintf(stderr,"QQ2 %lg\n",q2[0]);

  if (rank==0) {
    for (i=0; i < nmatres; i++)
      {
        r2[i] = b2[i] - q2[i];
        d2[i] = r2[i] / hit2[i];
      }
  }

  double delta_new_tmp = 0.0;
  if (rank==0) {
    for (i = 0; i < nmatres; i++) {
      delta_new_tmp += b2[i] ;
      if (isnan(b2[i])) {
        fprintf(stderr,"NAN B2 PBS %ld %lg\n",(long) i,q2[i]);
      }
      if (isnan(d2[i])) {
        fprintf(stderr,"NAN D2 PBS %ld %lg %lg\n",(long) i,hit2[i],q2[i]);
      }
    }
    //fprintf(stderr,"B2 %lg\n",delta_new_tmp);

    delta_new_tmp = 0.0;
    for (i = 0; i < nmatres; i++) {
      delta_new_tmp += q2[i] ;
      if (isnan(q2[i])) {
        fprintf(stderr,"NAN Q2 PBS %ld\n",(long) i);
      }
    }
    //fprintf(stderr,"Q2 %lg\n",delta_new_tmp);
    delta_new_tmp = 0.0;
    for (i = 0; i < nmatres; i++) {
      delta_new_tmp += q2[i]-b2[i] ;
      //fprintf(stderr,"B2 Q2 B2-Q2 [%ld]: %lg\t%lg\t%lg\n",(long) i,b2[i],q2[i],q2[i]-b2[i]);
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  delta_new_tmp = 0.0;
  if (rank==0) for (i=0; i < nmatres; i++) {
    delta_new_tmp += r2[i] * d2[i];
  }

  delta_new=0;
  for (i=0;i<mpi_size;i++) {
    double tmp=delta_new_tmp;
    MPI_Bcast(&tmp, sizeof(double), MPI_BYTE, i, MPI_COMM_WORLD);
    delta_new+=tmp;
  }

  if (itbogo==0) delta0 = delta_new;
  if (rank==0) fprintf (stderr, "iter = %d - delta0 = %lg - delta_new = %lg\n", iter, delta0, delta_new);
  int testwrit=0;

  if (itbogo==0) MPI_Bcast(&normaoff, sizeof(double), MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&delta0, sizeof(double), MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&delta_new, sizeof(double), MPI_BYTE, 0, MPI_COMM_WORLD);

// ################################################################################################################
// ################################################################################################################// ################################################################################################################// ################################################################################################################// ################################################################################################################// ################################################################################################################// ################################################################################################################
  while ((iter < itermax)  && ((delta_new) > delta0*1E-24) && ((delta_new) > 1E-20)) //Param->XI2STOP))
    {
      // q <= Ad
      //if (rank==0&&mindelta>delta_new) {
      if (rank==0) {
        memcpy(x2old,ix2,nmatres*sizeof(double));
        testwrit=1;
      }
    //  else testwrit=0;

      //healspline_fill_s_with_PtP_s (hs_rec, Ndata, vec, d, q);
      //for (i=0;i<mpi_size;i++) {
      //MPI_Bcast(d+tab_begr[i]*2, sizeof(double)*(tab_edr[i]-tab_begr[i]+1)*2, MPI_BYTE, i, MPI_COMM_WORLD);
      //}
      MPI_Bcast(d2, sizeof(double)*(nmatres), MPI_BYTE, 0, MPI_COMM_WORLD);

      // ===========================================================================
      // PROJECTION DE d dans q
      //
      for (l=0;l<nmatres;l++) q2[l]=0;
      if (rank==0) {
        double soff=0;
        for (i=0;i<newnr[nbolo];i++) soff+=hit2[0]*d2[i];
        for (i=0;i<newnr[nbolo];i++) q2[i]=soff;

        if (Param->REMHDIP==1) {
          soff=0;
          for (i=newnr[nbolo];i<newnr[nbolo]+GAINSTEP*nbolo;i++) soff+=hit2[newnr[nbolo]]*d2[i];
          for (i=newnr[nbolo];i<newnr[nbolo]+GAINSTEP*nbolo;i++) q2[i]=soff;
        }

        if (Param->flag_AVGR0==_PAR_TRUE) {
          soff=0;
          for (i=0;i<nbolo;i++) soff+=hit2[newnr[nbolo]+nbolo*(GAINSTEP2)]*
                                  d2[newnr[nbolo]+nbolo*(GAINSTEP2)+i];
          for (i=0;i<nbolo;i++) q2[newnr[nbolo]+nbolo*(GAINSTEP2)+i]=soff;
        }

        if (nmatco>0) {
          soff=0;
          if ((singleFreq == 0) && (Param->flag_AVG12CO == _PAR_TRUE)) {
            for (i=0;i<nbolo;i++) if (freqs[i] == Param->AVFCO)
                soff += hit2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)]
                        * d2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+i];
            for (i=0;i<nbolo;i++) if (freqs[i] == Param->AVFCO)
                q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+i] = soff;
          }
          else {
            for (i=0;i<nbolo;i++) soff+=hit2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)]
                                    *d2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+i];
            for (i=0;i<nbolo;i++) q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+i]=soff;
          }
        }

        if (nfreefree>0) {
          soff=0;
	  if ((singleFreq == 0) && (Param->flag_AVGFREEFREE == _PAR_TRUE)) {
	    for (i=0;i<nbolo;i++)  if (freqs[i] == Param->AVFFF) soff+=hit2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust]
				                              *d2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+i];
	    for (i=0;i<nbolo;i++) if (freqs[i] == Param->AVFFF)  q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+i+nmatco+nmatdust]=soff;
	  }
	  else {
	    for (i=0;i<nbolo;i++) soff+=hit2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust]
				    *d2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+i];
	    for (i=0;i<nbolo;i++) q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+i+nmatco+nmatdust]=soff;
	  }
	}

        if (nmatdust>0) {
          soff=0;
          if ((singleFreq == 0) && (Param->flag_AVGDUST100 == _PAR_TRUE)) {
            for (i=0;i<nbolo;i++) if (freqs[i] == Param->AVFDUST)
                soff += hit2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco]
                        * d2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+i];
            for (i=0;i<nbolo;i++) if (freqs[i] == Param->AVFDUST)
                q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+i] = soff;
          }
          else {
            for (i=0;i<nbolo;i++)
              soff+=hit2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco]*
                d2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+i];
            for (i=0;i<nbolo;i++) q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+i]=soff;
          }
        }

        long nj=0;
        if (DOFITANGLE==1) nj++;
        if (DOFITPOLEFF==1) nj++;
	if (DOTDUST==1) nj++;
        if (DOCO13==1) nj++;
	if (DOSYNCHRO==1) nj++;

	if (NORMFITPOL==1) {
#if 1
	  if (DOFITPOLEFF==1) {
	    j=1+DOCO13+DOSYNCHRO+DOTDUST;
	    soff=0;
	    for (i=0;i<nbolo;i++) soff+=hit2[newnr[nbolo]+nbolo*(GAINSTEP2+npixbeam-j)]*
				    d2[newnr[nbolo]+nbolo*(GAINSTEP2+npixbeam-j)+i];
	    for (i=0;i<nbolo;i++) q2[newnr[nbolo]+nbolo*(GAINSTEP2+npixbeam-j)+i]=soff;
	  }
#endif
	}
	
	if (DOFITANGLE==1) {
	  j=1+DOCO13+DOSYNCHRO+DOTDUST+DOFITPOLEFF;
	  soff=0;
	  for (i=0;i<nbolo;i++) soff+=hit2[newnr[nbolo]+nbolo*(GAINSTEP2+npixbeam-j)]*
				  d2[newnr[nbolo]+nbolo*(GAINSTEP2+npixbeam-j)+i];
	  for (i=0;i<nbolo;i++) q2[newnr[nbolo]+nbolo*(GAINSTEP2+npixbeam-j)+i]=soff;
	}

	if (DOCO13==1) {
	  j=1+DOSYNCHRO;
	  soff=0;
	  for (i=0;i<nbolo;i++) if (freqs[i] == Param->AVFCO) soff+=hit2[newnr[nbolo]+nbolo*(GAINSTEP2+npixbeam-j)]*
						       d2[newnr[nbolo]+nbolo*(GAINSTEP2+npixbeam-j)+i];
	  for (i=0;i<nbolo;i++) if (freqs[i] == Param->AVFCO) q2[newnr[nbolo]+nbolo*(GAINSTEP2+npixbeam-j)+i]=soff;
	}

	if (DOTDUST==1) {
	  j=1+DOCO13+DOSYNCHRO;
	  soff=0;
	  for (i=0;i<nbolo;i++) if (freqs[i] == MAXFREQ) soff+=hit2[newnr[nbolo]+nbolo*(GAINSTEP2+npixbeam-j)]*
							   d2[newnr[nbolo]+nbolo*(GAINSTEP2+npixbeam-j)+i];
	  for (i=0;i<nbolo;i++) if (freqs[i] == MAXFREQ) q2[newnr[nbolo]+nbolo*(GAINSTEP2+npixbeam-j)+i]=soff;
	}

	if (DOSYNCHRO==1) {
	    soff=0;
	    j=1;
	    for (i=0;i<nbolo;i++) if (freqs[i] == Param->AVFSYNC)  soff+=hit2[newnr[nbolo]+nbolo*(GAINSTEP2+npixbeam-j)]*
						       d2[newnr[nbolo]+nbolo*(GAINSTEP2+npixbeam-j)+i];
	    for (i=0;i<nbolo;i++) if (freqs[i] == Param->AVFSYNC)  q2[newnr[nbolo]+nbolo*(GAINSTEP2+npixbeam-j)+i]=soff;
	}
      }
      for (k=0;k<nnbpix;k++) if (flgpix[k]>0) {

        //long imat=the_stat_pix[k];
        long ndata = loc_nhpix[k];
        if (ndata>0) {
          hpix *htmp = loc_hpix[k];

          double vali=0;


#ifdef USEDII
          memset(dii ,0,sizeof(double)*nbolo*GAINSTEP);

          for (l1=0;l1<ndata;l1++) {
            long ri1=htmp[l1].rg-globalBeginRing;

            if (flg_rg[htmp[l1].ib][ri1]!=0) {

              dii[htmp[l1].gi+htmp[l1].ib*GAINSTEP ]+=htmp[l1].w*htmp[l1].dip;
            }
          }
          for (ib=0;ib<nbolo;ib++) {
            for (j=0;j<GAINSTEP;j++) {
              dii[j+ib*GAINSTEP]=dii[j+ib*GAINSTEP]/II[k];
            }
          }
#endif
          long l3;
          for (l3=0;l3<ndata;l3++) {
            long ri3=htmp[l3].rg-globalBeginRing;
            if (flg_rg[htmp[l3].ib][ri3]!=0) {
              long ir3=rgord[htmp[l3].ib][ri3]+newnr[htmp[l3].ib];
              vali+=htmp[l3].vi*d2[ir3];
#ifndef USEDII
              vali+=htmp[l3].lvi*d2[newnr[nbolo]+htmp[l3].ib*GAINSTEP+htmp[l3].gi];
#endif
            }
          }
          for (ib=0;ib<nbolo;ib++) {


#ifdef USEDII
            for (j=0;j<GAINSTEP;j++) {
              vali+=dii[j+ib*GAINSTEP]  *d2[newnr[nbolo]+ib*GAINSTEP+j];
            }
#endif

            if (nmatco>0) {
              vali+=dcoi[ib+k*nbolo]  *d2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+ib];
            }
            if (nmatdust>0) {
              vali+=ddusti[ib+k*nbolo]  *d2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+ib];
            }
            if (nfreefree>0) {
              vali+=dfri[ib+k*nbolo]  *d2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+ib];
            }
            for (m=0;m<npixbeam;m++) {
	      if (NOMOREFITTED!=DOCNN[m]||DOCNN[m]==0) {
		vali+=d2[newnr[nbolo]+nbolo*GAINSTEP2+ib+m*nbolo]*dpixi[ib+m*nbolo+k*nbolo*npixbeam];
	      }
            }
          }

          double qri=0;

          for (l1=0;l1<ndata;l1++) {
            long ri1=htmp[l1].rg-globalBeginRing;
            long iri1=rgord[htmp[l1].ib][ri1]+newnr[htmp[l1].ib];

            if (flg_rg[htmp[l1].ib][ri1]!=0) {
              double ww=htmp[l1].wp;
              double val2=d2[iri1]-(vali);

              val2+=d2[newnr[nbolo]+htmp[l1].ib*GAINSTEP+htmp[l1].gi]*htmp[l1].dip;
              if (nmatco>0)
                val2+=d2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+htmp[l1].ib]*htmp[l1].comap;
              if (nmatdust>0)
                val2+=d2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+htmp[l1].ib]*htmp[l1].dustmap;
              if (nfreefree>0)
                val2+=d2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+htmp[l1].ib]*htmp[l1].freefree;

              for (m=0;m<npixbeam;m++) {
		if (NOMOREFITTED!=DOCNN[m]||DOCNN[m]==0) {
		  val2+=d2[newnr[nbolo]+nbolo*GAINSTEP2+htmp[l1].ib+m*nbolo]*htmp[l1].listofpix[m];
		}
              }


              qri-=ww*val2;

              q2[iri1]+=ww*val2;

              q2[newnr[nbolo]+htmp[l1].ib*GAINSTEP+htmp[l1].gi]+=ww*val2*htmp[l1].dip;

              ///////////// CO
              if (nmatco>0) {
                q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+htmp[l1].ib]+=ww*val2*htmp[l1].comap;
              }

              ///////////// DUST
              if (nmatdust>0) {
                q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+htmp[l1].ib]+=ww*val2*htmp[l1].dustmap;
              }

              ///////////// FREEFREE
              if (nfreefree>0) {
                q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+htmp[l1].ib]+=ww*val2*htmp[l1].freefree;
              }


              for (j=0;j<npixbeam;j++) {
		if (NOMOREFITTED!=DOCNN[j]||DOCNN[j]==0) {
		  q2[newnr[nbolo]+nbolo*GAINSTEP2+j*nbolo+htmp[l1].ib]+=ww*val2*htmp[l1].listofpix[j];
		}
              }

            }
          }


          for (l2=0;l2<ndata;l2++) {
            long ri2=htmp[l2].rg-globalBeginRing;
            if (flg_rg[htmp[l2].ib][ri2]!=0) {
              long ir=rgord[htmp[l2].ib][ri2]+newnr[htmp[l2].ib];
              q2[ir]+=qri*htmp[l2].vi;
#ifndef USEDII
              q2[newnr[nbolo]+htmp[l2].ib*GAINSTEP+htmp[l2].gi]+=qri*htmp[l2].lvi;
#endif
            }
          }

          for (ib=0;ib<nbolo;ib++) {


#ifdef USEDII
            for (j=0;j<GAINSTEP;j++) {
              q2[newnr[nbolo]+GAINSTEP*ib+j]+=qri*dii[j+ib*GAINSTEP];
            }
#endif

            ///////////// CO
            if (nmatco>0) {
              q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+ib]+=qri*dcoi[ib+k*nbolo];
            }
            ///////////// DUST
            if (nmatdust>0) {
              q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+ib]+=qri*ddusti[ib+k*nbolo];
            }
            ///////////// FREEFREE
            if (nfreefree>0) {
              q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+ib]+=qri*dfri[ib+k*nbolo];
            }

            for (j=0;j<npixbeam;j++) {
	      if (NOMOREFITTED!=DOCNN[j]||DOCNN[j]==0) {
		q2[newnr[nbolo]+nbolo*GAINSTEP2+j*nbolo+ib]+=qri*dpixi[ib+j*nbolo+k*nbolo*npixbeam];
	      }
            }
          }
        }
      }
#ifdef OPTIMPI
      {

	double *lb = (double *) malloc(sizeof(double)*(nmatres));
	MPI_Reduce(q2,lb,nmatres,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

	memcpy(q2,lb,sizeof(double)*(nmatres));
	free(lb);
      }
#else
      if (rank==0) {
        double *lb = (double *) malloc(sizeof(double)*(nmatres));
        for (rrk=1;rrk<mpi_size;rrk++) {
          MPI_Recv(lb,sizeof(double)*(nmatres), MPI_BYTE, rrk,1034, MPI_COMM_WORLD,&statu);
          for (l=0;l<nmatres;l++) q2[l]+=lb[l];
        }
        free(lb);
      }
      else {
        MPI_Send(q2, sizeof(double)*(nmatres), MPI_BYTE, 0, 1034, MPI_COMM_WORLD);
      }
#endif
      double dtq = 0.0;
      if (rank==0) {
        for (i=0; i < nmatres; i++) dtq += d2[i] * q2[i]; //recalcul delta  (p *projX)
        alpha = delta_new / dtq;
        for (i=0; i < nmatres ; i++) ix2[i] += alpha * d2[i];
      }


      gettimeofday(&tp2,NULL);
      if (rank==0&&iter%10==0) fprintf (stderr,"iter = %ld-%d - delta_new = %lg %ld %3lfs\n", itbogo, iter, delta_new,
                          (long) testwrit,(double)(tp2.tv_sec-tp1.tv_sec)+(1E-6)*(tp2.tv_usec-tp1.tv_usec));

      //if (rank==0) fprintf (stderr,".");
      gettimeofday(&tp1,NULL);

      if (iter % 100 == 0 && iter !=0){ //REMOVE IF
          // Use the best case
          memcpy(ix2,x2old,nmatres*sizeof(double));
          //for (i=0;i<mpi_size;i++) {
          //  MPI_Bcast(x+tab_begr[i]*2, sizeof(double)*(tab_edr[i]-tab_begr[i]+1)*2, MPI_BYTE, i, MPI_COMM_WORLD);
          //}
          MPI_Bcast(ix2, sizeof(double)*(nmatres), MPI_BYTE, 0, MPI_COMM_WORLD);


      // ===========================================================================
      // PROJECTION DE x dans q
      // ===========================================================================

          for (l=0;l<nmatres;l++) q2[l]=0;
          if (rank==0) {
            double soff=0;
            for (i=0;i<newnr[nbolo];i++) soff+=hit2[0]*ix2[i];
            for (i=0;i<newnr[nbolo];i++) q2[i]=soff;

            if (Param->REMHDIP==1) {
              soff=0;
              for (i=newnr[nbolo];i<newnr[nbolo]+GAINSTEP*nbolo;i++) soff+=hit2[newnr[nbolo]]*ix2[i];
              for (i=newnr[nbolo];i<newnr[nbolo]+GAINSTEP*nbolo;i++) q2[i]=soff;
            }

            if (Param->flag_AVGR0==_PAR_TRUE) {
              soff=0;
              for (i=0;i<nbolo;i++) soff+=hit2[newnr[nbolo]+nbolo*(GAINSTEP2)]*
                                      ix2[newnr[nbolo]+nbolo*(GAINSTEP2)+i];
              for (i=0;i<nbolo;i++) q2[newnr[nbolo]+nbolo*(GAINSTEP2)+i]=soff;
            }

            if (nmatco>0) {
              soff=0;
              if ((singleFreq == 0) && (Param->flag_AVG12CO == _PAR_TRUE)) {
                for (i=0;i<nbolo;i++) if (freqs[i] == Param->AVFCO)
                    soff += hit2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)]
                            * ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+i];
                for (i=0;i<nbolo;i++) if (freqs[i] == Param->AVFCO)
                    q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+i] = soff;
              }
              else {
                for (i=0;i<nbolo;i++) soff+=hit2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)]
                                        *ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+i];
                for (i=0;i<nbolo;i++) q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+i]=soff;
              }
            }

            if (nmatdust>0) {
              soff=0;
              if ((singleFreq == 0) && (Param->flag_AVGDUST100 == _PAR_TRUE)) {
                for (i=0;i<nbolo;i++) if (freqs[i] == Param->AVFDUST)
                  soff += hit2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco]
                          * ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+i];
                for (i=0;i<nbolo;i++) if (freqs[i] == Param->AVFDUST)
                  q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+i] = soff;
              }
              else {
                for (i=0;i<nbolo;i++)
                  soff+=hit2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco]*
                    ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+i];
                for (i=0;i<nbolo;i++) q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+i]=soff;
              }
            }

            if (nfreefree>0) {
              soff=0;
	      if ((singleFreq == 0) && (Param->flag_AVGFREEFREE == _PAR_TRUE)) {
		for (i=0;i<nbolo;i++) if (freqs[i] == Param->AVFFF)
		  soff+=hit2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust]*
		    ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+i];
		for (i=0;i<nbolo;i++) if (freqs[i] == Param->AVFFF) q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+i]=soff;
	      }
	      else  {
		for (i=0;i<nbolo;i++)
		  soff+=hit2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust]*
		    ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+i];
		for (i=0;i<nbolo;i++) q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+i]=soff;
	      }
	    }

            long nj=0;
            if (DOFITANGLE==1) nj++;
            if (DOFITPOLEFF==1) nj++;
	    if (DOTDUST==1) nj++;
            if (DOCO13==1) nj++;
	    if (DOSYNCHRO==1) nj++;

	    if (NORMFITPOL==1) {
#if 1
	      if (DOFITPOLEFF==1) {
		j=1+DOCO13+DOSYNCHRO+DOTDUST;
		soff=0;
		for (i=0;i<nbolo;i++) soff+=hit2[newnr[nbolo]+nbolo*(GAINSTEP2+npixbeam-j)]*
					ix2[newnr[nbolo]+nbolo*(GAINSTEP2+npixbeam-j)+i];
		for (i=0;i<nbolo;i++) q2[newnr[nbolo]+nbolo*(GAINSTEP2+npixbeam-j)+i]=soff;
	      }
#endif
	    }
	    
	    if (DOFITANGLE==1) {
	      j=1+DOCO13+DOSYNCHRO+DOTDUST+DOFITPOLEFF;
	      soff=0;
	      for (i=0;i<nbolo;i++) soff+=hit2[newnr[nbolo]+nbolo*(GAINSTEP2+npixbeam-j)]*
				      ix2[newnr[nbolo]+nbolo*(GAINSTEP2+npixbeam-j)+i];
	      for (i=0;i<nbolo;i++) q2[newnr[nbolo]+nbolo*(GAINSTEP2+npixbeam-j)+i]=soff;
	    }

	    if (DOCO13==1) {
	      j=1+DOSYNCHRO;
	      soff=0;
	      for (i=0;i<nbolo;i++) if (freqs[i] == Param->AVFCO) soff+=hit2[newnr[nbolo]+nbolo*(GAINSTEP2+npixbeam-j)]*
							   ix2[newnr[nbolo]+nbolo*(GAINSTEP2+npixbeam-j)+i];
	      for (i=0;i<nbolo;i++) if (freqs[i] == Param->AVFCO) q2[newnr[nbolo]+nbolo*(GAINSTEP2+npixbeam-j)+i]=soff;
	    }

	    if (DOTDUST==1) {
	      j=1+DOCO13+DOSYNCHRO;
	      soff=0;
	      for (i=0;i<nbolo;i++) if (freqs[i] == MAXFREQ) soff+=hit2[newnr[nbolo]+nbolo*(GAINSTEP2+npixbeam-j)]*
							       ix2[newnr[nbolo]+nbolo*(GAINSTEP2+npixbeam-j)+i];
	      for (i=0;i<nbolo;i++) if (freqs[i] == MAXFREQ) q2[newnr[nbolo]+nbolo*(GAINSTEP2+npixbeam-j)+i]=soff;
	    }

	    if (DOSYNCHRO==1) {
	      soff=0;
	      j=1;
	      for (i=0;i<nbolo;i++) if (freqs[i] == Param->AVFSYNC)  soff+=hit2[newnr[nbolo]+nbolo*(GAINSTEP2+npixbeam-j)]*
							 ix2[newnr[nbolo]+nbolo*(GAINSTEP2+npixbeam-j)+i];
	      for (i=0;i<nbolo;i++) if (freqs[i] == Param->AVFSYNC)  q2[newnr[nbolo]+nbolo*(GAINSTEP2+npixbeam-j)+i]=soff;

	    }
          }


          for (k=0;k<nnbpix;k++) if (flgpix[k]>0) {

            //long imat=the_stat_pix[k];
            long ndata = loc_nhpix[k];
            if (ndata>0) {
              hpix *htmp = loc_hpix[k];

              double vali=0;

#ifdef USEDII
              memset(dii ,0,sizeof(double)*nbolo*GAINSTEP);


              for (l1=0;l1<ndata;l1++) {
                long ri1=htmp[l1].rg-globalBeginRing;

                if (flg_rg[htmp[l1].ib][ri1]!=0) {

                  dii[htmp[l1].gi+htmp[l1].ib*GAINSTEP ]+=htmp[l1].w*htmp[l1].dip;
                }
              }
              for (ib=0;ib<nbolo;ib++) {
                for (j=0;j<GAINSTEP;j++) {
                  dii[j+ib*GAINSTEP]=dii[j+ib*GAINSTEP]/II[k];
                }
              }
#endif
              long l3;
              for (l3=0;l3<ndata;l3++) {
                long ri3=htmp[l3].rg-globalBeginRing;
                if (flg_rg[htmp[l3].ib][ri3]!=0) {
                  long ir3=rgord[htmp[l3].ib][ri3]+newnr[htmp[l3].ib];
                  vali+=htmp[l3].vi*ix2[ir3];
#ifndef USEDII
                  vali+=htmp[l3].lvi*ix2[newnr[nbolo]+htmp[l3].ib*GAINSTEP+htmp[l3].gi];
#endif
                }
              }

              for (ib=0;ib<nbolo;ib++) {


#ifdef USEDII
                for (j=0;j<GAINSTEP;j++) {
                  vali+=dii[j+ib*GAINSTEP]  *ix2[newnr[nbolo]+ib*GAINSTEP+j];
                }
#endif

                if (nmatco>0) {
                  vali+=dcoi[ib+k*nbolo]  *ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+ib];
                }
                if (nmatdust>0) {
                  vali+=ddusti[ib+k*nbolo]  *ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+ib];
                }
                if (nfreefree>0) {
                  vali+=dfri[ib+k*nbolo]  *ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+ib];
                }
                for (m=0;m<npixbeam;m++) {
		  if (NOMOREFITTED!=DOCNN[m]||DOCNN[m]==0) {
		    vali+=ix2[newnr[nbolo]+nbolo*GAINSTEP2+ib+m*nbolo]*dpixi[ib+m*nbolo+k*nbolo*npixbeam];
		  }
                }
              }

              double qri=0;

              for (l1=0;l1<ndata;l1++) {
                long ri1=htmp[l1].rg-globalBeginRing;
                long iri1=rgord[htmp[l1].ib][ri1]+newnr[htmp[l1].ib];

                if (flg_rg[htmp[l1].ib][ri1]!=0) {
                  double ww=htmp[l1].wp;
                  double val2=ix2[iri1]-(vali);

                  val2+=ix2[newnr[nbolo]+htmp[l1].ib*GAINSTEP+htmp[l1].gi]*htmp[l1].dip;

                  if (nmatco>0)
                    val2+=ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+htmp[l1].ib]*htmp[l1].comap;
                  if (nmatdust>0)
                    val2+=ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+htmp[l1].ib]*htmp[l1].dustmap;
                  if (nfreefree>0)
                    val2+=ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nfreefree+htmp[l1].ib]*htmp[l1].freefree;

                  for (m=0;m<npixbeam;m++) {
		    if (NOMOREFITTED!=DOCNN[m]||DOCNN[m]==0) {
		      val2+=ix2[newnr[nbolo]+nbolo*GAINSTEP2+htmp[l1].ib+m*nbolo]*htmp[l1].listofpix[m];
		    }
                  }

                  qri-=ww*val2;

                  q2[iri1]+=ww*val2;

                  q2[newnr[nbolo]+htmp[l1].ib*GAINSTEP+htmp[l1].gi]+=ww*val2*htmp[l1].dip;

                  ///////////// CO
                  if (nmatco>0) {
                    q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+htmp[l1].ib]+=ww*val2*htmp[l1].comap;
                  }

                  ///////////// DUST
                  if (nmatdust>0) {
                    q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+htmp[l1].ib]+=ww*val2*htmp[l1].dustmap;
                  }

                  ///////////// FREEFREE
                  if (nfreefree>0) {
                    q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nfreefree+htmp[l1].ib]+=ww*val2*htmp[l1].freefree;
                  }

                  for (j=0;j<npixbeam;j++) {
		    if (NOMOREFITTED!=DOCNN[j]||DOCNN[j]==0) {
		      q2[newnr[nbolo]+nbolo*GAINSTEP2+j*nbolo+htmp[l1].ib]+=ww*val2*htmp[l1].listofpix[j];
		    }
                  }

                }
              }

              for (l2=0;l2<ndata;l2++) {
                long ri2=htmp[l2].rg-globalBeginRing;
                if (flg_rg[htmp[l2].ib][ri2]!=0) {
                  long ir=rgord[htmp[l2].ib][ri2]+newnr[htmp[l2].ib];
                  q2[ir]+=qri*htmp[l2].vi;
#ifndef USEDII
                  q2[newnr[nbolo]+htmp[l2].ib*GAINSTEP+htmp[l2].gi]+=qri*htmp[l2].lvi;
#endif
                }
              }

              for (ib=0;ib<nbolo;ib++) {


#ifdef USEDII
                for (j=0;j<GAINSTEP;j++) {
                  q2[newnr[nbolo]+GAINSTEP*ib+j]+=qri*dii[j+ib*GAINSTEP];
                }
#endif

                ///////////// CO
                if (nmatco>0) {
                  q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+ib]+=qri*dcoi[ib+k*nbolo];
                }
                ///////////// DUST
                if (nmatdust>0) {
                  q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+ib]+=qri*ddusti[ib+k*nbolo];
                }

                ///////////// FREEFREE
                if (nfreefree>0) {
                  q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+ib]+=qri*dfri[ib+k*nbolo];
                }

                for (j=0;j<npixbeam;j++) {
		  if (NOMOREFITTED!=DOCNN[j]||DOCNN[j]==0) {
		    q2[newnr[nbolo]+nbolo*GAINSTEP2+j*nbolo+ib]+=qri*dpixi[ib+j*nbolo+k*nbolo*npixbeam];
		  }
                }
              }
            }
          }
#ifdef OPTIMPI
	  {
	    double *lb = (double *) malloc(sizeof(double)*(nmatres));
	    MPI_Reduce(q2,lb,nmatres,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

	    memcpy(q2,lb,sizeof(double)*(nmatres));
	    free(lb);
	  }
#else
          if (rank==0) {
            double *lb = (double *) malloc(sizeof(double)*(nmatres));
            for (rrk=1;rrk<mpi_size;rrk++) {
              MPI_Recv(lb,sizeof(double)*(nmatres), MPI_BYTE, rrk,1034, MPI_COMM_WORLD,&statu);
              for (l=0;l<nmatres;l++) q2[l]+=lb[l];
            }

            free(lb);
          }
          else {
            MPI_Send(q2, sizeof(double)*(nmatres), MPI_BYTE, 0, 1034, MPI_COMM_WORLD);
          }
#endif

          if (rank==0) for (i=0; i < nmatres; i++) {
            r2[i] = b2[i] - q2[i];
          }

        }
      else
        {
          if (rank==0) for (i=0; i < nmatres; i++) r2[i] -= alpha * q2[i];
        }

      if (rank==0) for (i=0; i < nmatres; i++) s2[i] = r2[i] / hit2[i];

      delta_old = delta_new;
      if (rank==0) {
        delta_new=0;
        for (i=0; i < nmatres ; i++) delta_new += r2[i] * s2[i];
        beta = delta_new / delta_old;
        for (i=0; i < nmatres ; i++) d2[i] = s2[i] + beta * d2[i]; // new_p = new_r + beta*p
      }
      iter ++;
      MPI_Bcast(&delta_new, sizeof(double), MPI_BYTE, 0, MPI_COMM_WORLD);

    }
   // ===========================================================================
  // ===========================================================================    
  if (rank==0) fprintf (stderr,"\niter = %d - delta0 = %lg - delta_new = %lg\n",
                        iter, delta0, delta_new);
  if (rank==0) fprintf (stderr,"CG in iter = %d (max=%d)\n", iter, itermax);

  if (rank==0) memcpy(ix2  ,x2old,nmatres*sizeof (double));
  MPI_Bcast(ix2, sizeof(double)*(nmatres), MPI_BYTE, 0, MPI_COMM_WORLD);
  itbogo++;
}
