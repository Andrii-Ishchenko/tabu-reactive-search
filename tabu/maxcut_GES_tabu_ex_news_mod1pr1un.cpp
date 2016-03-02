#include <malloc.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <time.h>
#include <string.h>

long  aurd,nmy,n_point;
long nfail,total_fmax,n_rep,max_n_iter,n_att,n_attr,n_attempt,maxn_attempt=243;//2163;//81;//27;//27;//81;
long numb_th,kmax;
long fix_f,natak,m_c,coef,fz,fzl;
long *total_zz,*best_solution,*tabu_ar,*tabuf_ar,tabumin=21,tabumax=1000;
//long start,end,istop;
double start,beg,end;
long double epsl=1.e-3,epsu=0.9999,arec_t;
char ngdir[100],dirm[100],indirm[100],text;		   
long neqn,nVar,record_f,rr,n_edges;
long *record_sol;
long elite_size,max_elite_size,worst_member;
long t_tabu=7,tabu_moves,tabu_att;//20
long maxn_point=27,*best_freq;
//int npoint;
void xmain()
{
	void sort(long,long,long *,long double *);
    long double urand(void);
	void sort_simple_ins(long *,long *,long double *);
	void quick_sort(long *,long *,long double *);
	void CalcProb(int nmyu,long *iv,long double *y,
				long double *amyu,long double *sf,long double *nf,
				long double *sf1,long double *nf1,
				long *f0max,long *f1max,long double *y0);
	void ReCalc_mod(long f,long *iv,long double *amyu,int npoint,
			long *fmax,long *yz,long *zz,long *f0max,long *f1max,
			long double *sf,long double *nf,long double *sf1,long double *nf1,
			long *key_max,long key,long *gains,long *gains_max);
	void InitSetting(long *iv,int npoint,
					long *fmax,long *yz,long *zz,
					long *f0max,long *f1max,
					long double *sf,long double *nf,
					long double *sf1,long double *nf1);
	long mem_solution(long key,long f,long *numbf,long *keyar,long *vfar,
						long *left,long *right );
	long check_solution(long key,long f,long *numbf,long *keyar,long *vfar,
					long *left,long *right );
	void Random_1Opt_mod1(long *gg,long *moves,long *x,long *JA,long *JB,
				 long *BN,long *diag_a,long *upmoves,long *key,long *chash,
				 long numbf1,long *keyar1,long *vfar1,long *left1,long *right1,
				 long *iv,long *gains,long *xprev,long *xbest,long *best_gains);
	int EliteHandling_mode(long f,long *x,long *iv,
					long *val_f_elite,long *key_elite,long *numf,
					long *el_keyar,long *el_vfar,long *el_left,long *el_right,
					long *chash,_int8 *elite_sol,long keyf,
					long *gains,long *elite_gains,
					long *sum_sat,long *elite_sum_sat);
	void CalcMaxCutGoalf(long *g,long *x,long *edge_wt,long *iv,
						long *beg_list_edges_var,long *list_edges_var,long *list_nodes_var);
	void CalcMaxCutGains(long j,long *g,long *x,long *edge_wt,
						long *beg_list_edges_var,long *list_edges_var,long *list_nodes_var);
	void MaxCutRandom_KOpt_mod(long *gg,long nmoves,long *moves,long *x,
			long *upmoves,long *key,long *chash,
			long ex_numbf,long *ex_keyar,long *ex_vfar,long *ex_left,long *ex_right,
			long *iv,long *gains,long *xprev,long *xbest,long *best_gains,
			long *sum_sat,long *best_sum_sat,
			long *edge_wt,
			long *beg_list_edges_var,
			long *list_edges_var,
			long *beg_list_minus_edge_var,
			long *list_minus_edge_var,
			long *beg_list_plus_var_edge,
			long *list_plus_var_edge,
			long *beg_list_minus_var_edge,
			long *list_minus_var_edge);
	void ReCalcMaxCutGains(long j,long *gains,long *x,long *edge_wt,long *iv,
						long *beg_list_edges_var,long *list_edges_var,long *list_nodes_var);
	void MaxCutRandom_1_Opt(long *gg,long nmoves,long *moves,long *x,
			long *key,long *chash,long *iv,long *gains,long *xbest,long *best_gains,
			long *edge_wt,long *beg_list_edges_var,long *list_edges_var,long *list_nodes_var);
	void MaxCutRandom_KOpt_mod2(long *gg,long nmoves,long *moves,long *x,
			long *key,long *chash,
			long *iv,long *gains,long *xbest,long *best_gains,
			long *edge_wt,long *beg_list_edges_var,long *list_edges_var,long *list_nodes_var);
	void MaxCutRandom_Tabu(long *gg,long nmoves,long *moves,long *x,long *last_used,
			long *key,long *chash,long *numbf,long *keyar,long *vfar,long *left,long *right,
			long *iv,long *gains,long *xbest,long *best_gains,
			long *edge_wt,long *beg_list_edges_var,long *list_edges_var,long *list_nodes_var,
			long double *amyu,long *fmaxx,long *zz,long *f0max,long *f1max,
			long double *sf,long double *nf,long double *sf1,long double *nf1,
			long *key_maxx,long *gains_max);
	void MaxCutRandom_Tabu1(long *gg,long nmoves,long *moves,long *x,long *last_used,
			long *key,long *chash,long *numbf,long *keyar,long *vfar,long *left,long *right,
			long *iv,long *gains,long *xbest,long *best_gains,
			long *edge_wt,long *beg_list_edges_var,long *list_edges_var,long *list_nodes_var,
			long double *amyu,long *fmaxx,long *zz,long *f0max,long *f1max,
			long double *sf,long double *nf,long double *sf1,long double *nf1,
			long *key_maxx,long *gains_max);
	double cpu_time( void );
	void CalcProb1(int nmyu,long *iv,long double *y,
				long double *amyu,long double *sf,long double *nf,
				long double *sf1,long double *nf1,
				long *f0max,long *f1max,long double *y0);
	void MaxCutRandom_Tabu1Mod(long *gg,long nmoves,long *moves,long *x,long *last_used,
			long *key,long *chash,long *numbf,long *keyar,long *vfar,long *left,long *right,
			long *iv,long *gains,long *xbest,long *best_gains,
			long *edge_wt,long *beg_list_edges_var,long *list_edges_var,long *list_nodes_var,
			long double *amyu,long *fmaxx,long *zz,long *f0max,long *f1max,
			long double *sf,long double *nf,long double *sf1,long double *nf1,
			long *key_maxx,long *gains_max);

     extern long aurd;
     FILE *fi,*fi_MaxCut,*fi7;
     long double *y,*y0;
	 long double epsl=1.e-12,epsu=0.999999999999;
     long *chash,*yz,*iv,*zz,*team_zz,*key_ex_ar,*f_ex_ar,
			*iv1,*iv2,*iv3,*iv4,
			*gains,*xbest,*last_used;
     
     long i,ii,j,jj,jv,numit,maxnumit,nump,key,key_ex,n_ex,ex_flg=0,nmyu0,
				numthri,inumthri,
			numbf,numbf1,oldnumbf,memb;
	 
	 long double *sf, *nf, *sf1, *nf1;
	 long nfeval,*keyar, *left, *right,*vfar;
	 long cfmax,fmax,f,*f0max,*f1max,team_fmax;
     long double s,v,w,*amyu;
     long nmyu,v_nmyu,nvar_c;         
     long k,kk,*edge_wt;
     unsigned long rsize,isize,dim,ij;
	 char cc,dir1[100],dir[50],indir[50],dir3[50],dir_file[100],flg_rep;
	int npoint;	
	epsl=1.e-3,epsu=1.-epsl;
	strcpy_s(dir_file,"d:\\data_maxcut\\problem_name.txt");
	fopen_s(&fi_MaxCut,dir_file,"r");
	for (i=0; i<50; i++)
	{
		fscanf(fi_MaxCut,"%c",&cc);
		if (cc==10)
		{
			ngdir[i]=0;
			break;
		}
		ngdir[i]=cc;
	}
	fclose(fi_MaxCut);

	strcpy_s(indirm,"_bks.txt");
	//strcpy_s(indirm,"_bks_Resende.txt");
	strcpy_s(dirm,"d:\\data_MaxCut\\benchmarks\\");
	strcat(dirm,ngdir);
	strcat(dirm,indirm);
	fopen_s(&fi,dirm,"r");
	fscanf(fi,"%ld\n",&fz);
	fclose(fi);

	/*strcpy_s(indirm,"_bksl.txt");
	//strcpy_s(indirm,"_bks_Resende.txt");
	strcpy_s(dirm,"d:\\data_MaxCut\\benchmarks\\");
	strcat(dirm,ngdir);
	strcat(dirm,indirm);
	fopen_s(&fi,dirm,"r");
	fscanf(fi,"%ld\n",&fzl);
	fclose(fi);*/

	fzl=fz+100;
	fzl=fz;
	
	//fz=fz++;
	//fz=5507000;
	//fz=3923752;
	

	strcpy_s(indirm,".txt");
	strcpy_s(dirm,"d:\\data_MaxCut\\benchmarks\\");
	strcat(dirm,ngdir);
	strcat(dirm,indirm);
	fi=fopen(dirm,"r");
	fscanf(fi,"%ld%ld",&neqn,&n_edges);	
	fclose(fi);

	strcpy_s(dirm,"d:\\data_MaxCut\\res\\protokol_maxcut_GES_Tabu_ex_newsmod1pr1un_");
	strcat(dirm,ngdir);
	strcat(dirm,indirm);
	fopen_s(&fi7,dirm,"w");
	fclose(fi7);
	fopen_s(&fi7,dirm,"a+");

	nVar=neqn;
    
	npoint=200;//100;
	nfeval=10000000;
	//nfeval1=1250000;
	long cn_alg,cm_c,c_a_fmax;
	rsize=sizeof(long double);
    isize=sizeof(long);
	long n_moves;
	m_c=1;
	 
    rsize=sizeof(long double);
    isize=sizeof(long);
	dim=nfeval;
	keyar=(long *)  calloc(dim,isize);
    if (keyar==NULL)
    {
		ij=dim*isize;
        printf("It is not enough free memory for array keyar\n");
        goto LEND;
	}
	vfar=(long *)  calloc(dim,isize);
	if (vfar==NULL)goto LEND;
	left=(long *)  calloc(dim,isize);
    if (left==NULL)
    {
		ij=dim*isize;
        printf("It is not enough free memory for array left\n");
        goto LEND;
	}
	right=(long *)  calloc(dim,isize);
    if (right==NULL)
    {
		ij=dim*isize;
        printf("It is not enough free memory for array right\n");
        goto LEND;
	}

	dim=1000;
    best_freq=(long *)  calloc(dim,isize);
    if (best_freq==NULL)goto LEND;

	key_ex_ar=(long *)  calloc(dim,isize);
    if (key_ex_ar==NULL)goto LEND;
	f_ex_ar=(long *)  calloc(dim,isize);
    if (f_ex_ar==NULL)goto LEND;
	 
	long n_wt,*z,*n_con_edges;
	long oldrecord_f,*xz; 
	long double *wt_ar;

	dim=neqn;
    y=(long double *)  calloc(dim,rsize);
    if (y==NULL)
    {
       ij=dim*rsize;
       printf("It is not enough free memory for array y\n");
       goto LEND;
    }
	y0=(long double *)  calloc(dim,rsize);
    if (y0==NULL)
    {
       ij=dim*rsize;
       printf("It is not enough free memory for array y0\n");
       goto LEND;
    }
	wt_ar=(long double *)  calloc(dim,rsize);
	if (wt_ar==NULL)
    {
       ij=dim*rsize;
       printf("It is not enough free memory for array wt_ar\n");
       goto LEND;
    }
    f0max=(long *)  calloc(dim,isize);
    if (f0max==NULL)
    {
       ij=dim*isize;
       printf("It is not enough free memory for array f0max\n");
       goto LEND;
    }
	f1max=(long *)  calloc(dim,isize);
    if (f1max==NULL)
    {
       ij=dim*isize;
       printf("It is not enough free memory for array f1max\n");
       goto LEND;
    }
	chash=(long *)  calloc(dim,isize);
    if (chash==NULL)
    {
       ij=dim*isize;
       printf("It is not enough free memory for array chash\n");
       goto LEND;
    }
	n_con_edges=(long *)  calloc(dim,isize);
	if (n_con_edges==NULL)goto LEND;

	xz=(long *)  calloc(dim,isize);
    if (xz==NULL)goto LEND;
	best_solution=(long *)  calloc(dim,isize);
    if (best_solution==NULL)goto LEND;

    yz=(long *)  calloc(dim,isize);
    if (yz==NULL)
    {
		ij=dim*isize;
        printf("It is not enough free memory for array yz\n");
        goto LEND;
	}
    zz=(long *)  calloc(dim,isize);
    if (zz==NULL)
    {
       ij=dim*isize;
       printf("It is not enough free memory for array zz\n");
       goto LEND;
    }
	z=(long *)  calloc(dim,isize);
    if (z==NULL)
    {
		ij=dim*isize;
        printf("It is not enough free memory for array z\n");
        goto LEND;
       }
	team_zz=(long *)  calloc(dim,isize);
    if (team_zz==NULL)
    {
        ij=dim*isize;
        printf("It is not enough free memory for array team_zz\n");
        goto LEND;
    }
    iv=(long *)  calloc(dim,isize);
    if (iv==NULL)
    {
		ij=dim*isize;
		printf("It is not enough free memory for array iv\n");
        goto LEND;
    }
	iv1=(long *)  calloc(dim,isize);
    if (iv1==NULL)
       {
        ij=dim*isize;
        printf("It is not enough free memory for array iv\n");
        goto LEND;
       }
	 iv2=(long *)  calloc(dim,isize);
     if (iv2==NULL)
       {
        ij=dim*isize;
        printf("It is not enough free memory for array iv\n");
        goto LEND;
       }
	 iv3=(long *)  calloc(dim,isize);
     if (iv3==NULL)
       {
        ij=dim*isize;
        printf("It is not enough free memory for array iv\n");
        goto LEND;
       }
	iv4=(long *)  calloc(dim,isize);
    if (iv4==NULL)
       {
        ij=dim*isize;
        printf("It is not enough free memory for array iv4\n");
        goto LEND;
    }
	gains=(long *)  calloc(dim,isize);
    if (gains==NULL)goto LEND;
	xbest=(long *)  calloc(dim,isize);
    if (xbest==NULL)goto LEND;
	last_used=(long *)  calloc(dim,isize);
    if (last_used==NULL)goto LEND;
	tabu_ar=(long *)  calloc(dim,isize);
    if (tabu_ar==NULL)goto LEND;
	tabuf_ar=(long *)  calloc(dim,isize);
    if (tabuf_ar==NULL)goto LEND;

	long double *radius;
	dim=npoint+1;
	amyu=(long double *)  calloc(dim+1,sizeof(long double));
	if (amyu==NULL)goto LEND;
	radius=(long double *)  calloc(dim+1,sizeof(long double));
	if (radius==NULL)goto LEND;
	
	sf=(long double *)  calloc(dim,sizeof(long double));
	if (sf==NULL)
	{
		printf("It is not enough free memory for array sf\n");
		goto LEND;
	}
	nf=(long double *)  calloc(dim,sizeof(long double));
	if (nf==NULL)
	{
		printf("It is not enough free memory for array nf\n");
		goto LEND;
	}
	dim=neqn*dim;
	sf1=(long double *)  calloc(dim,sizeof(long double));
	if (sf1==NULL)
	{
		printf("It is not enough free memory for array sf1\n");
		goto LEND;
	}
	nf1=(long double *)  calloc(dim,sizeof(long double));
	if (nf1==NULL)
	{
		printf("It is not enough free memory for array nf1\n");
		goto LEND;
	}
	/*npoint=21;//26;//31;
	radius[0]=0.5*neqn;	
	radius[npoint-1]=20.0;
	s=(10.0-radius[0])/(npoint-1);
	for (i=1; i<npoint-1; i++)
	{
		radius[i]=i*s+radius[0];
		i=i;
	}*/
	
	aurd=90853; 
	v=1000.0;
	for (j=0; j<neqn; j++)
	{
		iv[j]=j;
		chash[j]=(long)(urand()*v);
	}


	long old_team_fmax,l,val_best_solution,gmax;
	
	long *moves,*best_gains,
		*gains_max,key_max,*sum_sat,*gains_mem;
	long *sum_plus_var,*sum_plus_edge,*sum_minus_var,
		*sum_minus_edge;

	dim=neqn;
	moves=(long *)  calloc(dim,isize);
    if (moves==NULL)goto LEND;
	record_sol=(long *)  calloc(dim,isize);
    if (record_sol==NULL)goto LEND;
	total_zz=(long *)  calloc(dim,isize);
    if (total_zz==NULL)goto LEND;
	best_gains=(long *)  calloc(dim,isize);
    if (best_gains==NULL)goto LEND;
	gains_max=(long *)  calloc(dim,isize);
    if (gains_max==NULL)goto LEND;
	gains_mem=(long *)  calloc(dim,isize);
    if (gains_mem==NULL)goto LEND;
	long kv,*beg_list_edges_var,*list_edges_var,*list_nodes_var;

	dim=n_edges+1;
	edge_wt=(long *)  calloc(dim,isize);
	if (edge_wt==NULL)goto LEND;
	sum_sat=(long *)  calloc(dim,isize);
	if (sum_sat==NULL)goto LEND;
	
	dim=neqn+1;
	beg_list_edges_var=(long *)  calloc(dim,isize);
	if (beg_list_edges_var==NULL)goto LEND;

/*===================================================================*/
	coef=(long)180000000/fz;
	//coef=1;
	fz=fz*coef;
	fzl=fzl*coef;

	strcpy_s(indirm,".txt");
	strcpy_s(dirm,"d:\\data_MaxCut\\benchmarks\\");
	strcat(dirm,ngdir);
	strcat(dirm,indirm);

	for (j=0; j<neqn; j++)n_con_edges[j]=0;
	fi=fopen(dirm,"r");
	fscanf(fi,"%ld%ld",&neqn,&n_edges);
	for (j=0; j<n_edges; j++)
	{
		fscanf(fi,"%ld%ld%ld",&k,&l,&i); k--; l--;
		edge_wt[j]=i;
		n_con_edges[k]++; n_con_edges[l]++;
	}
	fclose(fi);

	beg_list_edges_var[0]=0;
	for (j=0; j<neqn; j++)
	{
		beg_list_edges_var[j+1]=beg_list_edges_var[j]+n_con_edges[j];
	}

	dim=beg_list_edges_var[neqn];
	list_edges_var=(long *)  calloc(dim,isize);
	if (list_edges_var==NULL)goto LEND;
	list_nodes_var=(long *)  calloc(dim,isize);
	if (list_nodes_var==NULL)goto LEND;

	kv=99999999;
	for (j=0; j<n_edges; j++)
	{
		if (kv>edge_wt[j])kv=edge_wt[j];
	}
	kv=0;
	for (j=0; j<neqn; j++)n_con_edges[j]=0;
	fi=fopen(dirm,"r");
	fscanf(fi,"%ld%ld",&neqn,&n_edges);
	for (j=0; j<n_edges; j++)
	{
		fscanf(fi,"%ld%ld%ld",&k,&l,&i); k--; l--;
		edge_wt[j]=(i-kv)*coef;
		list_edges_var[beg_list_edges_var[k]+n_con_edges[k]]=j;
		list_edges_var[beg_list_edges_var[l]+n_con_edges[l]]=j;
		list_nodes_var[beg_list_edges_var[k]+n_con_edges[k]]=l;
		list_nodes_var[beg_list_edges_var[l]+n_con_edges[l]]=k;
		n_con_edges[k]++; n_con_edges[l]++;
	}
	fclose(fi);

	long gain,n_thr,keyf,mem,dist,df,df15[10],df20[10],f15[10],f20[10],n_cyc,
		old_fmax,max_dist,max_natak;
	long double e_dist,starttime,fintime,t_time,mean_t,mean_f;
/*===================================================================*/
	ij=0; max_natak=20; t_time=0.0;
	for (j=0; j<neqn; j++)moves[j]=j;
	for (j=0; j<1000; j++)best_freq[j]=0;
	mean_f=0.0;
	mean_t=0.0;
	//fz=13350*coef;
	for (natak=1; natak<=max_natak; natak++)
	{
		printf("natak=%ld\n",natak);
		//time(&start);
		start=cpu_time();
		numb_th=0; i=i;
//================================================================
		npoint=50;//21;//35;//30;//25;//32;//50;//50;i=i;
		n_point=npoint;
		amyu[0]=0.0;
		amyu[1]=7.e-7/coef;//7.e-7;//1//-7//-6

		s=exp( (log(3.0/coef)-log(amyu[1]))/(npoint-2.0) );///0.3/1.0//10//0.00001//0.01//0.1
		//s=exp( (log(0.07/coef)-log(amyu[1]))/19.0 );

		//s=exp( (log(10.0/coef)-log(amyu[1]))/(npoint-2.0) );///0.3/1.0//10//0.00001//0.01//0.1

		for (i=2; i<npoint; i++) amyu[i]=amyu[i-1]*s;
		amyu[npoint]=10.0;

		//npoint=21;//26;//31;
		radius[0]=0.5*neqn;	
		radius[npoint-1]=20.0;
		s=(10.0-radius[0])/(npoint-1);
		for (i=1; i<npoint-1; i++)
		{
			radius[i]=i*s+radius[0];
			i=i;
		}
//================================================================			
		nfail=0;
		total_fmax=0;
//===================================================================
		n_rep=0;
/*===================================================================*/
		aurd=natak*333;
/*===================================================================*/
		maxnumit=0;
		numit=0;

		//npoint=21;//21;//31;//30;

		total_fmax=0;
		record_f=team_fmax=fmax=0;
		for (j=0; j<neqn; j++)
		{
			zz[j]=record_sol[j]=0;
			iv[j]=j;
			y[j]=0.5; y0[j]=0.5; 
		}
		//max_dist=50;
		max_dist=neqn/2;//3;
		
		nVar=neqn;
		//n_point=0;
		starttime=cpu_time();
		numbf1=numbf=0;	
		n_attempt=1;
		for (j=0; j<9; j++)df15[j]=df20[j]=f15[j]=f20[j]=0;
		//df15=0;//10; 
		//df20=0;
		n_ex=0;
//==================================================================================
BEG:;	aurd=natak*777+777777*n_attempt+7;
		InitSetting(iv,npoint,&fmax,yz,zz,f0max,f1max,sf,nf,sf1,nf1);
		numbf1=numbf=0;	
		for (jv=0; jv<n_ex; jv++)
		{
			mem_solution(key_ex_ar[jv],f_ex_ar[jv],&numbf,keyar,vfar,left,right);
		}
		for (j=0; j<nVar; j++)
		{
			iv[j]=j; moves[j]=j;
			if (urand()<0.5)xz[j]=team_zz[j]=yz[j]=0;
			else xz[j]=team_zz[j]=yz[j]=1;
			tabu_ar[j]=tabuf_ar[j]=tabumin;
		}
		max_dist=neqn/2;///3;
		f=team_fmax=key=0;
		n_moves=nVar;
		CalcMaxCutGoalf(&i,xz,edge_wt,iv,beg_list_edges_var,list_edges_var,list_nodes_var);
		f=team_fmax=i;
		i=0;
		for (jv=0; jv<nVar; jv++)
		{
			j=iv[jv];
			if (xz[j])i=i+chash[j];
			CalcMaxCutGains(j,&l,xz,edge_wt,beg_list_edges_var,list_edges_var,list_nodes_var);
			gains[j]=l;
//			xz[j]=1-xz[j];
//			CalcMaxCutGoalf(&i,xz,edge_wt,iv,beg_list_edges_var,list_edges_var,list_nodes_var);
//			xz[j]=1-xz[j];
//			if (i-f!=l)
//			{
//				l=l;
//			}
		}
		key=i;

		nmy=nmyu=0;//Mod
		MaxCutRandom_Tabu1(&f,n_moves,moves,xz,last_used,
							&key,chash,&numbf,keyar,vfar,left,right,
							iv,gains,xbest,best_gains,edge_wt,
							beg_list_edges_var,list_edges_var,list_nodes_var,
							amyu,&fmax,zz,f0max,f1max,sf,nf,sf1,nf1,
							&key_max,gains_max);
		/*CalcMaxCutGoalf(&i,xz,edge_wt,iv,
						beg_list_edges_var,list_edges_var,list_nodes_var);
				if (i!=f)
				{
					printf("ERROR in value f\n");
				}*/

		/*MaxCutRandom_Tabu(&f,n_moves,moves,xz,last_used,
							&key,chash,&numbf,keyar,vfar,left,right,
							iv,gains,xbest,best_gains,edge_wt,
							beg_list_edges_var,list_edges_var,list_nodes_var,
							amyu,&fmax,zz,f0max,f1max,sf,nf,sf1,nf1,
							&key_max,gains_max);*/
		/*MaxCutRandom_KOpt_mod2(&f,n_moves,moves,xz,&key,chash,
							iv,gains,xbest,best_gains,
							edge_wt,beg_list_edges_var,list_edges_var,list_nodes_var);*/
		/*MaxCutRandom_1_Opt(&f,n_moves,moves,xz,
							&key,chash,
							iv,gains,xbest,best_gains,
							edge_wt,beg_list_edges_var,list_edges_var,list_nodes_var);*/

//		CalcMaxCutGoalf(&i,xz,edge_wt,iv,beg_list_edges_var,list_edges_var,list_nodes_var);
//		if (i!=f)
//		{
//			i=i;
//		}
		radius[0]=0.5*neqn;	
		radius[npoint-1]=20.0;
		s=(10.0-radius[0])/(npoint-1);
		for (i=1; i<npoint-1; i++)
		{
			radius[i]=i*s+radius[0];
			i=i;
		}
		
//==================================================================
		mem_solution(key,f,&numbf,keyar,vfar,left,right);

		ReCalc_mod(f,iv,amyu,npoint,&fmax,xz,zz,f0max,f1max,
							sf,nf,sf1,nf1,&key_max,key,gains,gains_max);
		
		if (fmax>=fzl)
		{
			for (j=0; j<neqn; j++)
			{
				yz[j]=xz[j];
			}
			goto SOL;	
		}
		flg_rep=0; n_cyc=0;
CCC:;
		n_rep++;
		
		aurd=natak*333+3333333*n_attempt+3;
		nmyu=nmyu0=0;
		oldrecord_f=fmax;
		for (jv=0; jv<nVar; jv++)
		{
			j=iv[jv];
			y[j]=y0[j];
		}

		while(nmyu<npoint)
		//while(nmyu<23)
		{
			if (nmyu>0)
			{
				v_nmyu=nmyu-1;
				CalcProb(v_nmyu,iv,y,amyu,sf,nf,sf1,nf1,f0max,f1max,y0);
				//CalcProb1(nmyu,iv,y,amyu,sf,nf,sf1,nf1,f0max,f1max,y0);
				/*for (jv=0; jv<nVar; jv++)
				{
					j=iv[jv];
					y[j]=y0[j];
				}
				if (nmyu==20 || nmyu==40)
				{
					nmyu=nmyu;
				}*/

			}
			else
			{
				for (jv=0; jv<nVar; jv++)
				{
					j=iv[jv];
					y[j]=y0[j];
					tabu_ar[j]=tabuf_ar[j]=tabumin;
				}
			}
			nmy=nmyu;
			
			//numthri=neqn;
			//if (nmyu<=18)numthri=50;
			//else numthri=100;
			//numthri=neqn/2; //super
			//numthri=numthri*n_rep;

			//numthri=500;//300
			numthri=81;//50;//81;//81;//243;//81;//300

			//numthri=243;
			//if (nmyu>=15)numthri=700;//500
			//numthri=100;

			//if (n_cyc>0)// && nmyu>15)
			//{
			//	numthri=243;
			//}

			//numthri=243;//81;//243;
			//if (nfail<1)numthri=100;
			//else numthri=500;
			//if (nmyu>=15)numthri=729;//243;//81;//243;//729;
			//if (nmyu>=15 && fmax==record_f)numthri=729;
			//if (nmy>=20 && fmax==record_f && nfail==1)numthri=729;
			//if (nmy>=20 && fmax==record_f)numthri=729;
			//if (nmy>=15 && fmax==record_f)numthri=729;
			oldnumbf=numbf; 

			//time(&beg);
			beg=cpu_time();
			cfmax=0;
			e_dist=0.0;
			n_thr=0;//-1
			old_fmax=fmax;
			for (jv=0; jv<nVar; jv++)
			{
				j=iv[jv];
				wt_ar[j]=-fabs(double(y[j]-0.5));
				//wt_ar[j]=fabs(double(y[j]-0.5));
			}
			for (jv=nVar-1; jv>0; jv--)
			{
				k=(long)(urand()*jv);
				i=iv[k]; iv[k]=iv[jv]; iv[jv]=i;
			}
			
			//if (nmyu%2!=0)
			{
				quick_sort(iv,iv+nVar-1,wt_ar);
			}

			while(n_thr<numthri)
			{
				n_thr++;
				rr=n_thr;
				key=key_max; f=fmax;
				for (jv=0; jv<neqn; jv++)
				{
					j=iv[jv];
					yz[j]=zz[j];
					gains[j]=gains_max[j];
				}
				/*CalcMaxCutGoalf(&i,yz,edge_wt,iv,beg_list_edges_var,list_edges_var,list_nodes_var);
				if (i!=f)
				{
					printf("ERROR in value f\n");
				}*/
				dist=0;
				for (jv=0; jv<nVar; jv++)
				{
					j=iv[jv];
					if (zz[j]==0)
					{
						if (urand()<=y[j])
						//if (urand()<=0.5)
						{
							ReCalcMaxCutGains(j,gains,yz,edge_wt,iv,
								beg_list_edges_var,list_edges_var,list_nodes_var);
							yz[j]=1; dist++;
							key=key+chash[j];
							f=f+gains[j];
//							CalcMaxCutGoalf(&i,yz,edge_wt,iv,
//											beg_list_edges_var,list_edges_var,list_nodes_var);
//							if (i!=f)
//							{
//								i=i;
//							}
							gains[j]=-gains[j];
							/*for (kv=0; kv<neqn; kv++)
							{
								k=iv[kv];
								if (k==j)continue;
								CalcMaxCutGains(k,&l,yz,sum_sat,edge_wt,iv,
												beg_list_edges_var,
												list_edges_var,
												beg_list_minus_edge_var,
												list_minus_edge_var,
												beg_list_plus_var_edge,
												list_plus_var_edge,
												beg_list_minus_var_edge,
												list_minus_var_edge);
								gains[k]=l;
							}*/
						}
					}
					else
					{
						if (urand()>y[j])
						//if (urand()>0.5)
						{
							ReCalcMaxCutGains(j,gains,yz,edge_wt,iv,
								beg_list_edges_var,list_edges_var,list_nodes_var);
							yz[j]=0; dist++;
							key=key-chash[j];
							f=f+gains[j];
//							CalcMaxCutGoalf(&i,yz,edge_wt,iv,
//											beg_list_edges_var,list_edges_var,list_nodes_var);
//							if (i!=f)
//							{
//								i=i;
//							}
							gains[j]=-gains[j];
							/*for (kv=0; kv<neqn; kv++)
							{
								k=iv[kv];
								if (k==j)continue;
								CalcMaxCutGains(k,&l,yz,sum_sat,edge_wt,iv,
												beg_list_edges_var,
												list_edges_var,
												beg_list_minus_edge_var,
												list_minus_edge_var,
												beg_list_plus_var_edge,
												list_plus_var_edge,
												beg_list_minus_var_edge,
												list_minus_var_edge);
								gains[k]=l;
							}*/
						}
					}
					
					//if (dist>max_dist)break;

					if (dist>radius[nmyu])break;

				}
				e_dist=e_dist+dist;
				//mem=check_solution(f,key,&numbf,keyar,vfar,left,right);
				//if (mem>0)
				//{
				//	continue;
				//}
				/*CalcMaxCutGoalf(&i,yz,edge_wt,iv,
						beg_list_edges_var,list_edges_var,list_nodes_var);
				if (i!=f)
				{
					printf("ERROR in value f\n");
				}*/
				MaxCutRandom_Tabu1Mod(&f,n_moves,moves,yz,last_used,
							&key,chash,&numbf,keyar,vfar,left,right,
							iv,gains,xbest,best_gains,edge_wt,
							beg_list_edges_var,list_edges_var,list_nodes_var,
							amyu,&fmax,zz,f0max,f1max,sf,nf,sf1,nf1,
							&key_max,gains_max);

				/*CalcMaxCutGoalf(&i,yz,edge_wt,iv,
						beg_list_edges_var,list_edges_var,list_nodes_var);
				if (i!=f)
				{
					printf("ERROR in value f\n");
				}*/

				if (cfmax<f)
				{
					cfmax=f;
				}
				
				mem=mem_solution(key,f,&numbf,keyar,vfar,left,right);
				if (mem>0)
				{
					continue;
				}
				ReCalc_mod(f,iv,amyu,npoint,&fmax,yz,zz,f0max,f1max,
							sf,nf,sf1,nf1,&key_max,key,gains,gains_max);
				if (fmax>=fzl)
				{
					goto SOL;	
				}
			}
			gmax=fmax;
			end=cpu_time();

			e_dist=e_dist/numthri;
			if (fmax>=fzl)
			{
				team_fmax=fmax;
				goto SOL;	
			}
			//if (end-start>=2250)
			//{
			//	goto SOL;
			//}

			//printf("time=%9.3f fmax=%8ld cf=%8ld eps=%10.6Lf nmyu=%2ld %ld\n",
			//	end-start,fmax/coef,cfmax/coef,100.0*(long double)(fz-cfmax)/(long double)fz,nmyu,nfail);
			//printf("numbf=%8ld numbf1=%8ld e_dist=%10.3Lf radius=%10.2Lf\n",
			//			numbf,numbf1,e_dist,radius[nmyu]);

			if (e_dist<20.0)break;
			/*if (nmyu==10)
			{
				f20[n_cyc]=fmax;
				k=df20[n_cyc];
				for (i=n_cyc+1; i<9; i++)k=k+df20[i];
				//if (fmax/coef<record_f/coef-k && n_attempt>9)break;// && n_cyc>0)
			}
			if (nmyu==15)
			{
				f15[n_cyc]=fmax;
				k=df15[n_cyc];
				for (i=n_cyc+1; i<9; i++)k=k+df20[i];
				//if (fmax/coef<record_f/coef-k && n_attempt>9)break;// && n_cyc>0)
			}*/
			nmyu++;
		}
		//if (df15[n_cyc]<fmax/coef-f15[n_cyc]/coef)df15[n_cyc]=fmax/coef-f15[n_cyc]/coef;
		//if (df20[n_cyc]<fmax/coef-f20[n_cyc]/coef)df20[n_cyc]=fmax/coef-f20[n_cyc]/coef;

		if (oldrecord_f>=fmax)nfail++;
		else nfail=0;
		printf("nfail=%ld fmax=%ld n_rep=%ld df15=%ld df10=%ld n_cyc=%ld nmyu=%ld\n",
			nfail,fmax/coef,n_rep,df15[n_cyc],df20[n_cyc],n_cyc,nmyu);
		//if ( nfail==1 && (record_f/coef-fmax/coef<=3) )
		//{
		//	flg_rep=0;//1;
		//}
	//	else flg_rep=0;		
		k=0;
		//for (i=n_cyc+1; i<9; i++)k=k+df20[i];
		
		//if (n_attempt<=9 && k<3)k=3;
		if (k<3)k=3;
		k=4;//20;//4;
		kmax=k;
		/*if (urand()<=0.33)
		{
			k=4;
		}
		else
		{
			if (urand()<=0.66)k=6;
			else k=8;
		}
		k=4;*/
		if ( nfail<1 && ( (record_f/coef-fmax/coef<=k) ))//|| n_attempt<=12) )//flg_rep==1)
		{
			f=fmax; key=key_max;
			for (jv=0; jv<neqn; jv++)
			{
				j=iv[jv];
				xz[j]=zz[j]; gains[j]=gains_max[j];
			}
			InitSetting(iv,npoint,&fmax,yz,zz,f0max,f1max,sf,nf,sf1,nf1);
			numbf=0;
			mem_solution(key,f,&numbf,keyar,vfar,left,right);
			ReCalc_mod(f,iv,amyu,npoint,&fmax,xz,zz,f0max,f1max,
					sf,nf,sf1,nf1,&key_max,key,gains,gains_max);
		
			i=i; n_cyc++;
			for (jv=0; jv<n_ex; jv++)
			{
				mem_solution(key_ex_ar[jv],f_ex_ar[jv],&numbf,keyar,vfar,left,right);
			}
			if (neqn/2>500)radius[0]=500;	
			else radius[0]=neqn/2;	
			radius[npoint-1]=20.0;
			s=(10.0-radius[0])/(npoint-1);
			for (i=1; i<npoint-1; i++)
			{
				radius[i]=i*s+radius[0];
				i=i;
			}
			goto CCC;
		}
		nfail=0;
		key_ex_ar[n_ex]=key_max; f_ex_ar[n_ex]=fmax; n_ex++;
		//max_dist=50;
		max_dist=neqn/2;
		//if (end-start<=2000)goto BEG;//2276)
		fintime=cpu_time();
		
		printf("time:%9.1Lf record_f=%ld n_attempt=%ld numbf=%ld n_ex=%ld\n",
			(fintime-starttime),record_f/coef,n_attempt,numbf,n_ex);
		n_attempt++;
		//if ( n_attempt<=maxn_attempt && (fintime-starttime)<=360.0 )goto BEG;
		if ( (fintime-starttime)<=3600.0 )goto BEG;
		//if ( n_attempt<=maxn_attempt )goto BEG;


SOL:;
		end=cpu_time();
		fintime=cpu_time();
		t_time=t_time+(fintime-starttime);
		mean_t=mean_t+arec_t;
		mean_f=mean_f+record_f/coef;
		printf("time = %8.3lf fmax=%ld\n",(fintime-starttime),fmax/coef);
		fprintf(fi7,"time = %8.3lf fmax=%ld\n",(fintime-starttime),fmax/coef);
		fflush(fi7);
		for (j=0; j<neqn; j++)yz[j]=record_sol[j];

		strcpy_s(indir,".txt");
		strcpy_s(dir1,"d:\\data_MaxCut\\results\\");
		strcat(dir1,ngdir);
		strcat(dir1,"_results_table_maxcut_GES_Tabu_ex_newsmod1pr1un.txt");
		fi=fopen(dir1,"a+");
		CalcMaxCutGoalf(&i,yz,edge_wt,iv,
						beg_list_edges_var,list_edges_var,list_nodes_var);
		if (i!=record_f)
		{
			printf("ERROR in value f\n");
			fprintf(fi,"ERROR in value f\n");
		}
		fprintf(fi," %3ld %9.2Lf %8.3f %8.3lf %9ld %9ld %9.6Lf %3ld %3ld\n",
		natak,arec_t,(end-start),(fintime-starttime),record_f/coef,fz/coef,100.0*(double)(fz-record_f)/(double)fz,n_attr,kmax);
		if (natak==max_natak)
		{
			fprintf(fi,"%9.2Lf %12.2Lf\n",mean_t/max_natak,mean_f/natak);
		}
		fclose(fi);

		strcpy_s(indir,".txt");
		strcpy_s(dir1,"d:\\data_MaxCut\\results\\");
		strcat(dir1,ngdir);
		strcat(dir1,"_results_solutions_maxcut_GES_Tabu_ex_newsmod1pr1un.txt");
		fi=fopen(dir1,"a+");
		for (j=0; j<neqn; j++)fprintf(fi," %ld",record_sol[j]);
		fprintf(fi,"\n");
		fclose(fi);

	  }
		strcpy_s(indir,".txt");
		strcpy_s(dir1,"d:\\data_MaxCut\\results\\");
		strcat(dir1,ngdir);
		strcat(dir1,"_results_solutions_maxcut_GES_Tabu_ex_newsmod1pr1un.txt");
		fi=fopen(dir1,"a+");
		for (j=0; j<numthri; j++)fprintf(fi," %ld %ld\n",j,best_freq[j]);
		fclose(fi);

LEND:;
     
}
/*-------------------------------------------------------------------*/
      #define ia 843314861
      #define ib 453816693
      #define mic 1693666955
      #define mmu 1073741824
      #define psu 4.65661287307739258e-10
      long double urand()
       {
        aurd=aurd*ia;
        if ( aurd > mic )aurd=(aurd-mmu)-mmu;
        aurd=aurd+ib;
        if ( aurd < 0 )aurd=(aurd+mmu)+mmu;
        return(aurd*psu);
       }
	long double ferf(long double xv)
    {
		long double  y,z,p,x;
        x=fabs(xv);
        y=1.0/(1.0+0.3275911*x);
        z=( ( 1.0614054*y-1.453152 )*y+1.4214137 )*y;
        z=( ( ( z-0.28449673 )*y+0.2548295 )*y )/exp( x*x );
        y=1.0/( 1.0+0.2316419*x );
        p=exp( -x*x/2.0 )/2.50662827;
        z=( (1.3302744*y-1.821256 )*y+1.7814779 )*y;
        z=1.0-( ( ( z-0.35656378 )*y+0.31938153 )*y )*p;
        if ( xv>=0.0 ) return(1.0-z);
        else return(z);
	}
	void integry(long neqn,long nmyu,long *iv,long double *y,
				long double *amyu,long double *sf,long double *nf,
				long double *sf1,long double *nf1,
				long *f0max,long *f1max,long double *y0
				)
	{
		long i,j,jv,ii,jj;
		long double s,v,w,u;
		for (jv=0; jv<neqn; jv++)
		{
			j=iv[jv];
			s=0.0; ii=0;
			for (i=0; i<=nmyu; i++)
			{
				jj=ii+j;
				if (nf1[jj]>1.e-12 && nf[i]-nf1[jj]>1.e-12)
				{
					v=sf1[jj]/nf1[jj];
					if (v>f1max[j])v=f1max[j];
					w=(sf[i]-sf1[jj])/(nf[i]-nf1[jj]);
					if (w>f0max[j])w=f0max[j];
				}
				else 
				{
					v=f1max[j];
					w=f0max[j];
				}
						
				u=w-v;

				jj=jj+neqn;
				if (nf1[jj]>1.e-12 && nf[i+1]-nf1[jj]>1.e-12)
				{
					v=sf1[jj]/nf1[jj];
					if (v>f1max[j])v=f1max[j];
					w=(sf[i+1]-sf1[jj])/(nf[i+1]-nf1[jj]);
					if (w>f0max[j])w=f0max[j];
				}
				else
				{
					v=f1max[j];
					w=f0max[j];
				}
						
				u=0.5*(u+w-v);
				s=s+(amyu[i+1]-amyu[i])*u;
				ii=ii+neqn;
			}
			y[j]=1.0/(1.0+(1.0-y0[j])*exp(s)/y0[j] );
			if (y[j]>epsu)y[j]=epsu;
			if (y[j]<epsl)y[j]=epsl;
		}
	}
	long mem_solution(long key,long *numbf,long *keyar,long *left,long *right )
	{
		long node,item;
		long uk=0;
		if (*numbf==0)
		{
			*numbf=*numbf+1;
			keyar[*numbf]=key;
			left[*numbf]=0;
			right[*numbf]=0;
		}
		else
		{
			node=1;
			while(node)
			{
				if (keyar[node]==key)
				{
					return node;
				}
				else
				{
					item=node;
					if (keyar[node]>key)node=left[node];
					else node=right[node];
				}
			}
			*numbf=*numbf+1;
			keyar[*numbf]=key;
			left[*numbf]=right[*numbf]=0;
			if (keyar[item]>key)left[item]=*numbf;
			else right[item]=*numbf;
		}
		return uk;
	}
	_int8 compare_solution(long key,long *numbf,long *keyar,long *left,long *right )
	{
		long node,item;
		_int8 uk=1;
		if (*numbf==0)return uk;
		else
		{
			node=1;
			while(node && uk)
			{
				if (keyar[node]==key)
				{
					uk=0;
				}
				else
				{
					item=node;
					if (keyar[node]>key)node=left[node];
					else node=right[node];
				}
			}
			return uk;	
		}
	}
	void quick_sort(long *low_ptr,long *high_ptr,long double *veight)
	{
		long *pivot_ptr;
		long *partition(long *low_ptr, long *high_ptr,long double *veight);
		void sort_simple_ins(long *low_ptr,long *high_ptr,long double *veight);
		if (low_ptr<high_ptr && high_ptr-low_ptr>9)
		{
			pivot_ptr=partition(low_ptr,high_ptr,veight);
			if (pivot_ptr-1-low_ptr>9)
				quick_sort(low_ptr,pivot_ptr-1,veight);
			else
				sort_simple_ins(low_ptr,pivot_ptr-1,veight);
			if (high_ptr-pivot_ptr>9)
				quick_sort(pivot_ptr,high_ptr,veight);
			else
				sort_simple_ins(pivot_ptr,high_ptr,veight);
		}
		else sort_simple_ins(low_ptr,high_ptr,veight);
	}
	long *partition(long *low_ptr, long *high_ptr,long double *veight)
	{
		void swap(long *ptr_1,long *ptr_2);
		long double pivot=*(veight+*(low_ptr+(high_ptr-low_ptr)/2) );
		while (low_ptr<=high_ptr)
		{
			while(*(veight+*low_ptr)>pivot)low_ptr++;
			while(*(veight+*high_ptr)<pivot)high_ptr--;
			if (low_ptr<=high_ptr)
				swap(low_ptr++,high_ptr--);
		}
		return low_ptr;
	}
	void swap(long *ptr_1,long *ptr_2)
	{
		long temp=*ptr_1;
		*ptr_1=*ptr_2;  *ptr_2=temp; 
	}

void sort_simple_insert(long low,long high,long *iv,long double *veight)
	{
		long double pivot;
		long i,j,temp;
		j=low;
		while(j<high)
		{
			i=j++; pivot=veight[iv[j]];
			temp=iv[j];
			while( i>=low)
			{
				if (veight[iv[i]]>=pivot)break;
				iv[i+1]=iv[i]; i--;
				
			}
			iv[i+1]=temp;
		}
	}
void sort_simple_ins(long *low_ptr,long *high_ptr,long double *veight)
	{
		long double pivot;
		long *ptr_1, *ptr_2,*temp,ivp;
		ptr_2=low_ptr;
		while(ptr_2<high_ptr)
		{
			ptr_1=ptr_2++; pivot=*(veight+*ptr_2);
			temp=ptr_2; ivp=*ptr_2;  
			while( ptr_1>=low_ptr)
			{
				if ( *(veight+*ptr_1) >= pivot)break;
				*temp=*ptr_1;
				temp=ptr_1--;
				
			}
			*temp=ivp;
		}
	}
	void set_res(long fmax,long *afmax,long numit)
	{
		long i,j,v,w,k,it;
		w=100; k=-1; it=numit;//it=numit/n;
		while(w<10000000)//1000000
		{
			for (j=1; j<10; j++)
			{
				v=j*w; k=k+1;
				if (it>v)continue;
				for (i=k; i<47; i++)//47
				{
					afmax[i]=fmax;
				}
				goto FIN;
			}
			w=10*w;
		}
FIN:;
	}
	void transp(long n,long *iv)
    {
     long double urand(void);
     long i,k,l;
     long double s;
     k=n-1;
     while(k>0)
     {
		s=urand();
		s=1.0/(1.0-s);
		s=log(s)/log(2.0);
		l=k-(long)(s);
		if (l<k)
        {
			if (l<0)break;
			i=iv[k]; iv[k]=iv[l]; iv[l]=i;
		}
		k=k-1;
      }
    }
	void transp1(long n,long *iv)
    {
     long double urand(void);
     long i,k,l;
     long double s;
     k=n-1;
     while(k>0)
     {
		s=urand(); 	s=1.0-s*s*s;
		l=(long)(k*s);
		if (l<k)
        {
			i=iv[k]; iv[k]=iv[l]; iv[l]=i;
		}
		k=k-1;
      }
    }
void CalcProb(int nmyu,long *iv,long double *y,
				long double *amyu,long double *sf,long double *nf,
				long double *sf1,long double *nf1,
				long *f0max,long *f1max,long double *y0)
{
	long i,j,jv,ii,jj,a,b;;
	long double s,v,w,u;
	long double epsl=1.e-7,epsu=0.9999999;
	//long double epsl=1.e-3,epsu=0.999999;
	a=4*(tabumax-tabumin);
	b=-a;
	for (jv=0; jv<nVar; jv++)
	{
		j=iv[jv];
		s=0.0; ii=0;
		for (i=0; i<=nmyu; i++)
		{
			jj=ii+j;
			if (nf1[jj]>1.e-12 && nf[i]-nf1[jj]>1.e-12)
			{
				v=sf1[jj]/nf1[jj];
				if (v>f1max[j])v=f1max[j];
				w=(sf[i]-sf1[jj])/(nf[i]-nf1[jj]);
				if (w>f0max[j])w=f0max[j];
			}
			else 
			{
				v=f1max[j];
				w=f0max[j];
			}
			u=w-v;
			jj=jj+neqn;
			if (nf1[jj]>1.e-12 && nf[i+1]-nf1[jj]>1.e-12)
			{
				v=sf1[jj]/nf1[jj];
				if (v>f1max[j])v=f1max[j];
				w=(sf[i+1]-sf1[jj])/(nf[i+1]-nf1[jj]);
				if (w>f0max[j])w=f0max[j];
			}
			else
			{
				v=f1max[j];
				w=f0max[j];
			}
						
			u=0.5*(u+w-v);
			s=s+(amyu[i+1]-amyu[i])*u;
			ii=ii+neqn;
		}
		y[j]=1.0/(1.0+(1.0-y0[j])*exp(s)/y0[j] );
		if (y[j]>epsu)y[j]=epsu;
		if (y[j]<epsl)y[j]=epsl;

		/*if (y[j]<=0.6 && y[j]>=0.4)tabu_ar[j]=tabumin;
		else
		{
			tabu_ar[j]=(long)(a*y[j]*y[j]+b*y[j]+tabumax);
		}
		if (y[j]<=0.5)
		{
			tabu_ar[j]=(long)(tabumin*(1-y[j])/y[j]);
		}
		else
		{
			tabu_ar[j]=(long)(tabumin*y[j]/(1-y[j]));
		}*/
		if (y[j]<=0.6 && y[j]>=0.4)tabu_ar[j]=tabumin;
		else
		{
			if (y[j]<=0.5)
			{
				//s=(1-y[j])/y[j];
				tabu_ar[j]=(long)(tabumin*(1-y[j])/y[j]);
				//if (tabu_ar[j]>27*tabumin)tabu_ar[j]=27*tabumin;
				tabuf_ar[j]=(long)(tabumin*y[j]/(1-y[j]));//tabumin;//
				//if (tabuf_ar[j]<7)tabu_ar[j]=7;
			}
			else
			{
				tabu_ar[j]=(long)(tabumin*y[j]/(1-y[j]));
				//if (tabu_ar[j]>27*tabumin)tabu_ar[j]=27*tabumin;
				tabuf_ar[j]=(long)(tabumin*(1-y[j])/y[j]);//tabumin;//
				//if (tabuf_ar[j]<7)tabu_ar[j]=7;
			}
		}
		/*if (nmy>20)
		{
			if (y[j]>0.5 && f0max[j]>f1max[j])
			{
				j=j;
			}
			if (y[j]<0.5 && f0max[j]<f1max[j])
			{
				j=j;
			}
		}*/
	}
}
void CalcProb1(int nmyu,long *iv,long double *y,
				long double *amyu,long double *sf,long double *nf,
				long double *sf1,long double *nf1,
				long *f0max,long *f1max,long double *y0)
{
	long i,j,jv,ii,jj;
	long double s,v,w,u;
	long double epsl=1.e-7,epsu=0.9999999;
	//long double epsl=1.e-3,epsu=0.999999;

	for (jv=0; jv<nVar; jv++)
	{
		j=iv[jv];
		s=0.0; ii=nmyu*neqn;
		y[j]=nf1[ii+j]/nf[nmyu];
		if (y[j]>epsu)y[j]=epsu;
		if (y[j]<epsl)y[j]=epsl;
		/*if (nmy>20)
		{
			if (y[j]>0.5 && f0max[j]>f1max[j])
			{
				j=j;
			}
			if (y[j]<0.5 && f0max[j]<f1max[j])
			{
				j=j;
			}
		}*/
	}
}
void ReCalc_mod(long f,long *iv,long double *amyu,int npoint,
			long *fmax,long *yz,long *zz,long *f0max,long *f1max,
			long double *sf,long double *nf,long double *sf1,long double *nf1,
			long *key_max,long key,long *gains,long *gains_max)
{
	double cpu_time( void );
	FILE*fi;
	char dirm[100],dir[100];		   
	long i,jv,j,ii,jj;
	long double s,v;

	if (f==*fmax)f=f-coef+1;
	npoint=0;
	for (jv=0; jv<neqn; jv++)
	{
		j=iv[jv];
		if (yz[j])
		{
			if (f>f1max[j])f1max[j]=f;
		}
		else
		{
			if (f>f0max[j])f0max[j]=f;
		}
	}
								
	if (f>*fmax)
	{
		if (*fmax>0.0)
		{
			ii=neqn;
			for (i=1; i<=npoint; i++)
			{
				s=exp( amyu[i]*(long double)(*fmax-f) );
				sf[i]=sf[i]*s;
				nf[i]=nf[i]*s;
				for (jv=0; jv<nVar; jv++)
				{
					j=iv[jv];
					jj=ii+j;
					sf1[jj]=sf1[jj]*s;
					nf1[jj]=nf1[jj]*s;
				}
				ii=ii+neqn;
			}
		}
		*fmax=f; 
		*key_max=key;
		if (f>total_fmax)
		{
			total_fmax=f; //n_sprob=0;
			for (jv=0; jv<nVar; jv++)
			{
				j=iv[jv];
				total_zz[j]=yz[j];
				//if (yz[j]==1)i=i+c[j];
			}
		}
		if (f>record_f)
		{
			//time(&end);
			end=cpu_time();
			strcpy_s(dirm,"d:\\data_MaxCut\\res\\MaxCutrecords_GES_tabu_ex_newsmod1pr1un_");
			strcat(dirm,ngdir);
			strcat(dirm,"_");
			sprintf(dir,"%d",natak);
			strcat(dirm,dir);
			strcat(dirm,".txt");
			fopen_s(&fi,dirm,"a+");
			fprintf(fi,"%8.3f %ld  ",(end-start),f/coef);
			record_f=f; s=v=0;
			for (j=0; j<neqn; j++)
			{
				s=s+fabs((double)( best_solution[j]-yz[j]));
				v=v+fabs((double)(record_sol[j]-yz[j]));
				record_sol[j]=yz[j];
				//fprintf(fi,"%2ld",record_sol[j]);
			}
			//fprintf(fi,"\n");
			//fprintf(fi,"time(sec)=%6ld distg=%g distp=%g f=%ld nmy=%ld n_thr=%ld %ld\n",
			//	(end-start),s,v,f/coef,nmy,rr,key);
			fprintf(fi,"%8.3f %ld %ld %ld %ld %ld %ld\n",
				(end-start),f/coef,nmy,rr,key,n_att,kmax);
			printf("time(sec)=%8.3f distg=%g distp=%g f=%ld nmy=%ld n_thr=%ld %ld %ld kmax=%ld\n",
				(end-start),s,v,f/coef,nmy,rr,key,n_att,kmax);
			fclose(fi);
			arec_t=end-start; n_attr=n_attempt;
		}
		for (jv=0; jv<neqn; jv++)
		{
			j=iv[jv];
			zz[j]=yz[j];
			gains_max[j]=gains[j];
		}
		best_freq[rr]++;
	}
	v=(long double)f;
	sf[0]=sf[0]+v;
	nf[0]=nf[0]+1.0;
	for (jv=0; jv<nVar; jv++)
	{
		j=iv[jv];
		if (yz[j])
		{
			sf1[j]=sf1[j]+v;
			nf1[j]=nf1[j]+1.0;
		}
	}
	ii=neqn;
	for (i=1; i<=npoint; i++)
	{
		s=exp( amyu[i]*(long double)(f-*fmax) );
		sf[i]=sf[i]+v*s;
		nf[i]=nf[i]+s;
		for (jv=0; jv<nVar; jv++)
		{
			j=iv[jv];
			jj=ii+j;
			if (yz[j])
			{
				sf1[jj]=sf1[jj]+v*s;
				nf1[jj]=nf1[jj]+s;
			}
		}
		ii=ii+neqn;
	}
}
void ReCalc_mod9(long f,long *iv,long double *amyu,
			long *fmax,long *yz,long *zz,long *f0max,long *f1max,
			long double *sf,long double *nf,long double *sf1,long double *nf1,
			long *key_max,long key,long *gains,long *gains_max)
{
	double cpu_time( void );
	FILE*fi;
	char dirm[100],dir[100];		   
	long i,jv,j,ii,jj;
	long double s,v;

	/*if (f==*fmax)f=f-coef+1;

	for (jv=0; jv<neqn; jv++)
	{
		j=iv[jv];
		if (yz[j])
		{
			if (f>f1max[j])f1max[j]=f;
		}
		else
		{
			if (f>f0max[j])f0max[j]=f;
		}
	}*/
								
	if (f>*fmax)
	{
		/*if (*fmax>0.0)
		{
			ii=neqn;
			for (i=1; i<=npoint; i++)
			{
				s=exp( amyu[i]*(long double)(*fmax-f) );
				sf[i]=sf[i]*s;
				nf[i]=nf[i]*s;
				for (jv=0; jv<nVar; jv++)
				{
					j=iv[jv];
					jj=ii+j;
					sf1[jj]=sf1[jj]*s;
					nf1[jj]=nf1[jj]*s;
				}
				ii=ii+neqn;
			}
		}*/
		*fmax=f; 
		*key_max=key;
		if (f>total_fmax)
		{
			total_fmax=f; //n_sprob=0;
			for (jv=0; jv<nVar; jv++)
			{
				j=iv[jv];
				total_zz[j]=yz[j];
				//if (yz[j]==1)i=i+c[j];
			}
		}
		if (f>record_f)
		{
			//time(&end);
			end=cpu_time();
			strcpy_s(dirm,"d:\\data_MaxCut\\res\\MaxCutrecords_Glover_tabu_");
			strcat(dirm,ngdir);
			strcat(dirm,"_");
			sprintf(dir,"%d",natak);
			strcat(dirm,dir);
			strcat(dirm,".txt");
			fopen_s(&fi,dirm,"a+");
			fprintf(fi,"%8.3f %ld  ",(end-start),f/coef);
			record_f=f; s=v=0;
			for (j=0; j<neqn; j++)
			{
				s=s+fabs((double)( best_solution[j]-yz[j]));
				v=v+fabs((double)(record_sol[j]-yz[j]));
				record_sol[j]=yz[j];
				//fprintf(fi,"%2ld",record_sol[j]);
			}
			//fprintf(fi,"\n");
			//fprintf(fi,"time(sec)=%6ld distg=%g distp=%g f=%ld nmy=%ld n_thr=%ld %ld\n",
			//	(end-start),s,v,f/coef,nmy,rr,key);
			fprintf(fi,"%8.3f %ld %ld %ld %ld %ld %ld\n",
				(end-start),f/coef,n_rep,nmy,rr,key,n_att);
			printf("time(sec)=%8.3f distg=%g distp=%g f=%ld nmy=%ld n_thr=%ld %ld %ld\n",
				(end-start),s,v,f/coef,nmy,rr,key,n_att);
			fclose(fi);
			arec_t=end-start;
		}
		for (jv=0; jv<neqn; jv++)
		{
			j=iv[jv];
			zz[j]=yz[j];
			gains_max[j]=gains[j];
		}
	}
	/*v=(long double)f;
	sf[0]=sf[0]+v;
	nf[0]=nf[0]+1.0;
	for (jv=0; jv<nVar; jv++)
	{
		j=iv[jv];
		if (yz[j])
		{
			sf1[j]=sf1[j]+v;
			nf1[j]=nf1[j]+1.0;
		}
	}
	ii=neqn;
	for (i=1; i<=npoint; i++)
	{
		s=exp( amyu[i]*(long double)(f-*fmax) );
		sf[i]=sf[i]+v*s;
		nf[i]=nf[i]+s;
		for (jv=0; jv<nVar; jv++)
		{
			j=iv[jv];
			jj=ii+j;
			if (yz[j])
			{
				sf1[jj]=sf1[jj]+v*s;
				nf1[jj]=nf1[jj]+s;
			}
		}
		ii=ii+neqn;
	}*/
}
void ReCalc_sim(long f,long *iv,int npoint,long double *amyu,
			long *fmax,long *yz,long *zz,long *f0max,long *f1max,
			long double *sf,long double *nf,long double *sf1,long double *nf1,
			long *key_max,long key,long *gains,long *gains_max)
{
	double cpu_time( void );
	FILE*fi;
	char dirm[100],dir[100];		   
	long jv,j;
	long double s;

	//if (f==*fmax)f=f-1;

	for (jv=0; jv<nVar; jv++)
	{
		j=iv[jv];
		if (yz[j])
		{
			if (f>f1max[j])f1max[j]=f;
		}
		else
		{
			if (f>f0max[j])f0max[j]=f;
		}
	}
								
	if (f>*fmax)
	{
		*fmax=f; 
		*key_max=key;
		if (f>total_fmax)
		{
			total_fmax=f; //n_sprob=0;
			for (jv=0; jv<neqn; jv++)
			{
				j=iv[jv];
				total_zz[j]=yz[j];
				//if (yz[j]==1)i=i+c[j];
			}
		}
		if (f>record_f)
		{
			//time(&end);
			end=cpu_time();
			strcpy_s(dirm,"d:\\data_MaxCut\\res\\records_upd_new3_");
			strcat(dirm,ngdir);
			strcat(dirm,"_");
			sprintf(dir,"%d",natak);
			strcat(dirm,dir);
			strcat(dirm,".txt");
			fi=fopen(dirm,"a+");
			fprintf(fi,"%8.3f %ld  ",(end-start),f/coef);
			record_f=f; s=0;
			for (j=0; j<neqn; j++)
			{
				s=s+fabs((double)(best_solution[j]-yz[j]));
				record_sol[j]=yz[j];
				//fprintf(fi,"%2ld",record_sol[j]);
			}
			fprintf(fi,"\n");
			fprintf(fi,"time(sec)=%8.3f dist=%g f=%ld nmy=%ld n_thr=%ld %ld\n",
				(end-start),s,f/coef,nmy,rr,key);
			printf("time(sec)=%8.3f dist=%g f=%ld nmy=%ld n_thr=%ld %ld\n",
				(end-start),s,f/coef,nmy,rr,key);
			fclose(fi);
			arec_t=end-start;
		}
		for (jv=0; jv<neqn; jv++)
		{
			j=iv[jv];
			zz[j]=yz[j];
			gains_max[j]=gains[j];
		}
	}
}
void ReCalc(long f,long *iv,int npoint,long double *amyu,
	long *fmax,long *yz,long *zz,long *f0max,long *f1max,
	long double *sf,long double *nf,long double *sf1,long double *nf1)
{
	double cpu_time( void );
	FILE*fi;
	char dirm[100],dir[100];		   
	long i,jv,j,ii,jj;
	long double s,v;

	if (f==*fmax)f=f-1;

	for (jv=0; jv<nVar; jv++)
	{
		j=iv[jv];
		if (yz[j])
		{
			if (f>f1max[j])f1max[j]=f;
		}
		else
		{
			if (f>f0max[j])f0max[j]=f;
		}
	}
								
	if (f>*fmax)
	{
		if (*fmax>0.0)
		{
			ii=neqn;
			for (i=1; i<=npoint; i++)
			{
				s=exp( amyu[i]*(long double)(*fmax-f) );
				sf[i]=sf[i]*s;
				nf[i]=nf[i]*s;
				for (jv=0; jv<nVar; jv++)
				{
					j=iv[jv];
					jj=ii+j;
					sf1[jj]=sf1[jj]*s;
					nf1[jj]=nf1[jj]*s;
				}
				ii=ii+neqn;
			}
		}
		*fmax=f; 
		if (f>total_fmax)
		{
			total_fmax=f; //n_sprob=0;
			for (jv=0; jv<neqn; jv++)
			{
				j=iv[jv];
				total_zz[j]=yz[j];
				//if (yz[j]==1)i=i+c[j];
			}
		}
		if (f>record_f)
		{
			//time(&end);
			end=cpu_time();
			strcpy_s(dirm,"d:\\data_MaxCut\\res\\records_upd2_");
			strcat(dirm,ngdir);
			strcat(dirm,"_");
			sprintf(dir,"%d",natak);
			strcat(dirm,dir);
			strcat(dirm,".txt");
			fopen_s(&fi,dirm,"a+");
			fprintf(fi,"%8.3f %ld  ",(end-start),f/coef);
			record_f=f; s=0;
			for (j=0; j<neqn; j++)
			{
				s=s+fabs((double)(record_sol[j]-yz[j]));
				record_sol[j]=yz[j];
				fprintf(fi,"%2ld",record_sol[j]);
			}
			fprintf(fi,"\n");
			fprintf(fi,"time(sec)=%8.3f dist=%g f=%ld nmy=%ld n_thr=%ld\n",
				(end-start),s,f/coef,nmy,rr);
			printf("time(sec)=%8.3f dist=%g f=%ld nmy=%ld n_thr=%ld\n",
				(end-start),s,f/coef,nmy,rr);
			fclose(fi);
			arec_t=end-start;
		}
		for (jv=0; jv<neqn; jv++)
		{
			j=iv[jv];
			zz[j]=yz[j];
		}
	}
	v=(long double)f;
	sf[0]=sf[0]+v;
	nf[0]=nf[0]+1.0;
	for (jv=0; jv<nVar; jv++)
	{
		j=iv[jv];
		if (yz[j])
		{
			sf1[j]=sf1[j]+v;
			nf1[j]=nf1[j]+1.0;
		}
	}
	ii=neqn;
	for (i=1; i<=npoint; i++)
	{
		s=exp( amyu[i]*(long double)(f-*fmax) );
		sf[i]=sf[i]+v*s;
		nf[i]=nf[i]+s;
		for (jv=0; jv<nVar; jv++)
		{
			j=iv[jv];
			jj=ii+j;
			if (yz[j])
			{
				sf1[jj]=sf1[jj]+v*s;
				nf1[jj]=nf1[jj]+s;
			}
		}
		ii=ii+neqn;
	}
}
long MemSolution(long key,long *numbf,long *keyar,long *left,long *right )
{
	long node,item;
	if (*numbf==0)
	{
		*numbf=*numbf+1;
		keyar[*numbf]=key;
		left[*numbf]=0;
		right[*numbf]=0;
	}
	else
	{
		node=1;
		while(node)
		{
			if (keyar[node]==key)return (node);
			else
			{
				item=node;
				if (keyar[node]>key)
				{
					if (left[node]>1000000)
					{
						node=node;
					}
					node=left[node];
				}
				else 
				{
					if (right[node]>1000000)
					{
						node=node;
					}
					node=right[node];
				}
			}
		}
		*numbf=*numbf+1;
		keyar[*numbf]=key;
		left[*numbf]=right[*numbf]=0;
		if (keyar[item]>key)left[item]=*numbf;
		else right[item]=*numbf;
	}
	return (0);
}
void InitSetting(long *iv,int npoint,
					long *fmax,long *yz,long *zz,
					long *f0max,long *f1max,
					long double *sf,long double *nf,
					long double *sf1,long double *nf1)
{
	long i,j,jv,ii,jj;
	sf[0]=0.0;
	nf[0]=0.0;
	*fmax=0;
	for (jv=0; jv<nVar; jv++)
	{
		j=iv[jv]; yz[j]=zz[j]=0; 
		sf1[j]=0.0;
		nf1[j]=0.0;
		f1max[j]=f0max[j]=0;
	}
	ii=neqn;
	for (i=1; i<=npoint; i++)
	{
		sf[i]=0.0;
		nf[i]=0.0;
		for (jv=0; jv<nVar; jv++)
		{
			j=iv[jv];
			jj=ii+j;
			sf1[jj]=0.0;
			nf1[jj]=0.0;
		}
		ii=ii+neqn;
	}
}
int EliteHandling_mode(long f,long *x,long *iv,
					long *val_f_elite,long *key_elite,long *numf,
					long *el_keyar,long *el_vfar,long *el_left,long *el_right,
					long *chash,_int8 *elite_sol,long keyf,
					long *gains,long *elite_gains,
					long *sum_sat,long *elite_sum_sat)
{
	long mem_solution(long key,long f,long *numbf,long *keyar,long *vfar,
						long *left,long *right );
	long check_solution(long key,long f,long *numbf,long *keyar,long *vfar,
					long *left,long *right );
	long key,jv,ii,j,k,elite_numbf,ij;
	
	elite_numbf=*numf;
	key=keyf;

//начала работы с элитой
	//dist=0;
	//for (jv=0; jv<neqn; jv++)
	//{
	//	j=iv[jv];
	//	dist=dist+(long)fabs(x[j]-record_sol[j]);
	//}
	//if (dist>nVar)goto LEND;

	if (elite_size<max_elite_size)
	{
		if (mem_solution(key,f,&elite_numbf,el_keyar,el_vfar,el_left,el_right)==0)
		{
			val_f_elite[elite_size]=f;
			key_elite[elite_size]=key;
			ij=neqn*elite_size;
			for (jv=0; jv<nVar; jv++)
			{
				j=iv[jv];
				elite_sol[ij+j]=(_int8)x[j];
				elite_gains[ij+j]=gains[j];
			}
			for (j=0; j<n_edges; j++)
			{
				elite_sum_sat[ij+j]=sum_sat[j];
			}
			elite_size++;
			if (elite_size==max_elite_size)
			{
				ii=999999999;
				for (jv=0; jv<elite_size; jv++)
				{
					if (ii>val_f_elite[jv])
					{
						ii=val_f_elite[jv];
						worst_member=jv;
					}
				}
			}
		}
	}
	else
	{
		if (val_f_elite[worst_member]<f)
		{
			if (check_solution(key,f,&elite_numbf,el_keyar,el_vfar,el_left,el_right)==0)
			{
				k=val_f_elite[worst_member];
				val_f_elite[worst_member]=f;
				key_elite[worst_member]=key;
				ij=worst_member*neqn;
				for (jv=0; jv<nVar; jv++)
				{
					j=iv[jv];
					elite_sol[ij+j]=(_int8)x[j];
					elite_gains[ij+j]=gains[j];
				}
				for (j=0; j<n_edges; j++)
				{
					elite_sum_sat[ij+j]=sum_sat[j];
				}
				ii=999999999; elite_numbf=0;
				for (jv=0; jv<elite_size; jv++)
				{
					mem_solution(key_elite[jv],val_f_elite[jv],&elite_numbf,
						el_keyar,el_vfar,el_left,el_right);
					if (ii>val_f_elite[jv])
					{
						ii=val_f_elite[jv];
						worst_member=jv;
					}
				}
			}
		}
	}
//конец работы с элитой
	*numf=elite_numbf;
	

	return(0);
}
long mem_solution(long key,long f,long *numbf,long *keyar,long *vfar,
						long *left,long *right )
{
	long node,item;
	long uk=0;
	if (*numbf==0)
	{
		*numbf=*numbf+1;
		keyar[*numbf]=key;
		vfar[*numbf]=f;
		left[*numbf]=0;
		right[*numbf]=0;
	}
	else
	{
		node=1;
		while(node)
		{
			if (keyar[node]==key && vfar[node]==f)
			{
				return node;
			}
			else
			{
				item=node;
				if (keyar[node]>key)node=left[node];
				else node=right[node];
			}
		}
		*numbf=*numbf+1;
		keyar[*numbf]=key;
		vfar[*numbf]=f;
		left[*numbf]=right[*numbf]=0;
		if (keyar[item]>key)left[item]=*numbf;
		else right[item]=*numbf;
	}
	return uk;
}
long check_solution(long key,long f,long *numbf,long *keyar,long *vfar,
					long *left,long *right )
{
	long node,item;
	long uk=0;
	if (*numbf>0)
	{
		node=1;
		while(node)
		{
			if (keyar[node]==key && vfar[node]==f)
			{
				return node;
			}
			else
			{
				item=node;
				if (keyar[node]>key)node=left[node];
				else node=right[node];
			}
		}
	}
	return uk;
}
void Random_2Opt_mod(long *gg,long *moves,long *x,long *JA,long *JB,
				 long *BN,long *diag_a,long *upmoves,long *key,long *chash,
				 long numbf1,long *keyar1,long *vfar1,long *left1,long *right1,
				 long *iv,long *gains,long *xprev,long *xbest,long *best_gains)
{
	long double urand(void);
	void calc_goalf(long *g,long *x,long *JA,long *JB,long *BM);
	long check_solution(long key,long f,long *numbf,long *keyar,long *vfar,
					long *left,long *right );
	void calc_df(long k,long *g,long *x,long *JA,long *JB,
			 long *BM,long *diag_a);
	void Random_1Opt_mod1(long *gg,long *moves,long *x,long *JA,long *JB,
			 long *BN,long *diag_a,long *upmoves,long *key,long *chash,
			 long numbf1,long *keyar1,long *vfar1,long *left1,long *right1,
			 long *iv,long *gains,long *xprev,long *xbest,long *best_gains);


	long i,ii,j,f,iaa,iab,kk,k,kv,up_moves,mem,
		keyf,jv,g,keybest,
		flg_stop=0,s;

	*upmoves=0;
	f=*gg;
	keyf=*key;
	up_moves=0;
	goto BBB;

AAA:;
	Random_1Opt_mod1(&f,moves,x,JA,JB,BN,diag_a,&up_moves,
					&keyf,chash,numbf1,keyar1,vfar1,left1,right1,
					iv,gains,xprev,xbest,best_gains);
	*upmoves=*upmoves+up_moves;

BBB:;
	for (jv=nVar-1; jv>0; jv--)
	{
		k=(long)(urand()*jv);
		i=moves[k]; moves[k]=moves[jv]; moves[jv]=i;
	}
	for (kv=0; kv<nVar; kv++)
	{
		k=moves[kv];
		g=gains[k];
		x[k]=1-x[k];
		for (jv=0; jv<nVar; jv++)
		{
			best_gains[jv]=gains[jv];
		}
		keybest=keyf;
		s=0;
		gains[k]=-gains[k];
		if (x[k])
		{
			keyf=keyf+chash[k];
			iaa=JA[k]; iab=JA[k+1];
			for (ii=iaa; ii<iab; ii++)
			{
				i=JB[ii];
				if (i!=k)
				{
					gains[i]=gains[i]+2*BN[ii]*(1-2*x[i]);
					if (gains[i]>s)
					{
						s=gains[i];
						kk=i;
					}
				}
			}
		}
		else
		{
			keyf=keyf-chash[k];
			iaa=JA[k]; iab=JA[k+1];
			for (ii=iaa; ii<iab; ii++)
			{
				i=JB[ii];
				if (i!=k)
				{
					gains[i]=gains[i]-2*BN[ii]*(1-2*x[i]);
					if (gains[i]>s)
					{
						s=gains[i];
						kk=i;
					}
				}
			}
		}
		if (s+g>0)
		{
			g=g+s;
			j=kk;
			x[j]=1-x[j];
			gains[j]=-gains[j];
			if (x[j])
			{
				keyf=keyf+chash[j];
				iaa=JA[j]; iab=JA[j+1];
				for (ii=iaa; ii<iab; ii++)
				{
					i=JB[ii];
					if (i!=j)
					{
						gains[i]=gains[i]+2*BN[ii]*(1-2*x[i]);
					}
				}
			}
			else
			{
				keyf=keyf-chash[j];
				iaa=JA[j]; iab=JA[j+1];
				for (ii=iaa; ii<iab; ii++)
				{
					i=JB[ii];
					if (i!=j)
					{
						gains[i]=gains[i]-2*BN[ii]*(1-2*x[i]);
					}
				}
			}
			f=f+g;
			mem=check_solution(keyf,f,&numbf1,keyar1,vfar1,left1,right1);
			if (mem>0)
			{
				goto LEND;
			}
			else
			{
				goto AAA;
			}
		}
		else
		{
			x[k]=1-x[k];
			for (jv=0; jv<nVar; jv++)
			{
				gains[jv]=best_gains[jv];
			}
			keyf=keybest;
		}
	}
	
LEND:;

	*gg=f;
	*key=keyf;
}
void Random_1Opt_mod1(long *gg,long *moves,long *x,long *JA,long *JB,
				 long *BN,long *diag_a,long *upmoves,long *key,long *chash,
				 long numbf1,long *keyar1,long *vfar1,long *left1,long *right1,
				 long *iv,long *gains,long *xprev,long *xbest,long *best_gains)
{
	long double urand(void);
	void calc_goalf(long *g,long *x,long *JA,long *JB,long *BM);
	long check_solution(long key,long f,long *numbf,long *keyar,long *vfar,
					long *left,long *right );
	void calc_df(long k,long *g,long *x,long *JA,long *JB,
			 long *BM,long *diag_a);

	long i,ii,j,f,iaa,iab,k,keyf,jv,g,flg_stop=0;

	*upmoves=0;
	f=*gg;
	keyf=*key;
	for (jv=nVar-1; jv>0; jv--)
	{
		k=(long)(urand()*jv);
		i=moves[k]; moves[k]=moves[jv]; moves[jv]=i;
	}
	while(1)
	{
		g=0;
		for (jv=0; jv<nVar; jv++)
		{
			j=moves[jv];
			if (gains[j]>=0)
			{
				g=g+gains[j];
				(*upmoves)++;
				x[j]=1-x[j];
				if (x[j]==1)keyf=keyf+chash[j];
				else keyf=keyf-chash[j];
				gains[j]=-gains[j];
				iaa=JA[j]; iab=JA[j+1];
				for (ii=iaa; ii<iab; ii++)
				{
					i=JB[ii];
					if (i!=j)
					{
						gains[i]=gains[i]+2*BN[ii]*(1-2*x[i])*(2*x[j]-1);
					}
				}
			}
		}
		if (g>0)
		{
			f=f+g;
		}
		else break;
		//mem=check_solution(keyf,f,&numbf1,keyar1,vfar1,left1,right1);
		//if (mem>0)
		//{
		//	flg_stop=1; break;
		//}
	}
	*gg=f;
	*key=keyf;
}
void Random_1Opt_mod2(long *gg,long nmoves,long *moves,long *x,long *JA,long *JB,
				 long *BN,long *diag_a,long *upmoves,long *key,long *chash,
				 long ex_numbf,long *ex_keyar,long *ex_vfar,long *ex_left,long *ex_right,
				 long *iv,long *gains,long *xprev,long *xbest,long *best_gains,
				 long *last_used,long *nstep)
{
	long double urand(void);
	void calc_goalf(long *g,long *x,long *JA,long *JB,long *BM);
	long check_solution(long key,long f,long *numbf,long *keyar,long *vfar,
					long *left,long *right );
	void calc_df(long k,long *g,long *x,long *JA,long *JB,
			 long *BM,long *diag_a);

	long i,ii,j,f,iaa,iab,k,
		keyf,jv,g,
		flg_stop=0;

	*upmoves=0;
	f=*gg;
	keyf=*key;

	//keybest=keyf;
	//fbest=f;
	for (jv=nmoves-1; jv>0; jv--)
	{
		k=(long)(urand()*jv);
		i=moves[k]; moves[k]=moves[jv]; moves[jv]=i;
	}
	while(1)
	{
		g=0;
		for (jv=0; jv<nVar; jv++)
		{
			j=moves[jv];
			if ((*nstep)-last_used[j]<t_tabu)continue;
			if (gains[j]>=0)
			{
				g=g+gains[j];
				f=f+gains[j];
				(*upmoves)++;
				(*nstep)++;
				last_used[j]=*nstep;
				x[j]=1-x[j];
				gains[j]=-gains[j];
				iaa=JA[j]; iab=JA[j+1];
				if (x[j]==1)
				{
					keyf=keyf+chash[j];
					for (ii=iaa; ii<iab; ii++)
					{
						i=JB[ii];
						if (i!=j)
						{
							gains[i]=gains[i]+2*BN[ii]*(1-2*x[i]);
						}
					}
				}
				else 
				{
					for (ii=iaa; ii<iab; ii++)
					{
						i=JB[ii];
						if (i!=j)
						{
							gains[i]=gains[i]-2*BN[ii]*(1-2*x[i]);
						}
					}
					keyf=keyf-chash[j];
				}
			}
		}
		if (check_solution(keyf,f,&ex_numbf,ex_keyar,ex_vfar,ex_left,ex_right)==0)
		{
			if (g==0)break;
		}
		else break;
	}
	*gg=f;
	*key=keyf;
}
//improvement
void CalcMaxCutGoalf(long *g,long *x,long *edge_wt,long *iv,
						long *beg_list_edges_var,long *list_edges_var,long *list_nodes_var)
{
	long j,i,jv,k,kv,f;

	f=0;
	for (jv=0; jv<neqn; jv++)
	{
		j=iv[jv];
		if (x[j]==1)
		{
			for (kv=beg_list_edges_var[j]; kv<beg_list_edges_var[j+1]; kv++)
			{
				k=list_nodes_var[kv];
				if (x[k]==0)
				{
					f=f+edge_wt[list_edges_var[kv]];
				}
			}
		}
	}
	*g=f;
}
void CalcMaxCutGains(long j,long *g,long *x,long *edge_wt,
						long *beg_list_edges_var,long *list_edges_var,long *list_nodes_var)
{
	long k,kv,df;

	df=0;
	if (x[j]==1)
	{
		for (kv=beg_list_edges_var[j]; kv<beg_list_edges_var[j+1]; kv++)
		{
			k=list_nodes_var[kv];
			if (x[k]==0)
			{
				df=df-edge_wt[list_edges_var[kv]];
			}
			else
			{
				df=df+edge_wt[list_edges_var[kv]];
			}
		}
	}
	else
	{
		for (kv=beg_list_edges_var[j]; kv<beg_list_edges_var[j+1]; kv++)
		{
			k=list_nodes_var[kv];
			if (x[k]==0)
			{
				df=df+edge_wt[list_edges_var[kv]];
			}
			else
			{
				df=df-edge_wt[list_edges_var[kv]];
			}
		}
	}
	*g=df;
}
void ReCalcMaxCutGains(long j,long *gains,long *x,long *edge_wt,long *iv,
						long *beg_list_edges_var,long *list_edges_var,long *list_nodes_var)
{
	long k,kv,l,lv;
	
	if (x[j]==1)
	{
		for (kv=beg_list_edges_var[j]; kv<beg_list_edges_var[j+1]; kv++)
		{
			k=list_nodes_var[kv];
			if (x[k]==0)
			{
				gains[k]=gains[k]+2*edge_wt[list_edges_var[kv]];
			}
			else
			{
				gains[k]=gains[k]-2*edge_wt[list_edges_var[kv]];
			}
		}
	}
	else
	{
		for (kv=beg_list_edges_var[j]; kv<beg_list_edges_var[j+1]; kv++)
		{
			k=list_nodes_var[kv];
			if (x[k]==0)
			{
				gains[k]=gains[k]-2*edge_wt[list_edges_var[kv]];
			}
			else
			{
				gains[k]=gains[k]+2*edge_wt[list_edges_var[kv]];
			}
		}
	}

}
void MaxCutRandom_1_Opt(long *gg,long nmoves,long *moves,long *x,
			long *key,long *chash,long *iv,long *gains,long *xbest,long *best_gains,
			long *edge_wt,long *beg_list_edges_var,long *list_edges_var,long *list_nodes_var)
{
	long double urand(void);

	void CalcMaxCutGains(long j,long *g,long *x,long *edge_wt,
						long *beg_list_edges_var,long *list_edges_var,long *list_nodes_var);
	void ReCalcMaxCutGains(long j,long *gains,long *x,long *edge_wt,long *iv,
						long *beg_list_edges_var,long *list_edges_var,long *list_nodes_var);

	long i,j,f,k,keyf,jv,g,flg_stop=0;

	
	f=*gg;
	keyf=*key;
	
	while(1)
	{
		for (jv=nVar-1; jv>0; jv--)
		{
			k=(long)(urand()*jv);
			i=moves[k]; moves[k]=moves[jv]; moves[jv]=i;
		}
		g=0;
		for (jv=0; jv<nVar; jv++)
		{
			j=moves[jv];
			if (gains[j]>=0)
			{
				g=g+gains[j];
				ReCalcMaxCutGains(j,gains,x,edge_wt,iv,
								beg_list_edges_var,list_edges_var,list_nodes_var);
				x[j]=1-x[j];
				if (x[j]==1)keyf=keyf+chash[j];
				else keyf=keyf-chash[j];
				gains[j]=-gains[j];
			}
		}
		if (g>0)
		{
			f=f+g;
		}
		else break;
	}
	*gg=f;
	*key=keyf;
}
void MaxCutRandom_KOpt_mod2(long *gg,long nmoves,long *moves,long *x,
			long *key,long *chash,
			long *iv,long *gains,long *xbest,long *best_gains,
			long *edge_wt,long *beg_list_edges_var,long *list_edges_var,long *list_nodes_var)
{
	long double urand(void);

	void CalcMaxCutGains(long j,long *g,long *x,long *edge_wt,
						long *beg_list_edges_var,long *list_edges_var,long *list_nodes_var);
	void ReCalcMaxCutGains(long j,long *gains,long *x,long *edge_wt,long *iv,
						long *beg_list_edges_var,long *list_edges_var,long *list_nodes_var);

	long i,j,f,k,kv,n_moves,keyf,jv,gmax,g,gain,step,last_best_move,keybest,
		flg_stop=0,fbest,end_move,cg;

	/*CalcMaxCutGoalf(&i,x,sum_sat,edge_wt,iv,
					beg_list_edges_var,
					list_edges_var,
					beg_list_minus_edge_var,
					list_minus_edge_var,
					beg_list_plus_var_edge,
					list_plus_var_edge,
					beg_list_minus_var_edge,
					list_minus_var_edge);
	if (i!=(*gg))
	{
		i=i;
	}*/
	f=*gg;
	keyf=*key;
	step=last_best_move=0; 
	for (jv=0; jv<nVar; jv++)
	{
		i=iv[jv];
		best_gains[i]=gains[i];
		xbest[i]=x[i];
	}
	keybest=keyf;
	fbest=f;

	
	
AAA:;
	for (jv=nmoves-1; jv>0; jv--)
	{
		k=(long)(urand()*jv);
		i=moves[k]; moves[k]=moves[jv]; moves[jv]=i;
	}

	
	g=gmax=0;

	n_moves=nVar;
	//beg_move=0;
	end_move=nVar;
	last_best_move=step;

	while(end_move>0)
	{
		while(1)
		{
			for (jv=end_move-1; jv>0; jv--)
			{
				k=(long)(urand()*jv);
				i=moves[k]; moves[k]=moves[jv]; moves[jv]=i;
			}
			//kk=end_move+tabu;
			//for (jv=nmoves-1-kk; jv>0; jv--)
			//{
			//	k=(long)(urand()*jv)+kk;
			//	i=moves[k]; moves[k]=moves[jv]; moves[jv]=i;
			//}
			cg=0;
			for (kv=0; kv<end_move; kv++)
			{
				j=moves[kv];
				if (gains[j]>0)
				{
					cg=cg+gains[j];
					//g=g+gains[j];
					//f=f+gains[j];
					ReCalcMaxCutGains(j,gains,x,edge_wt,iv,
								beg_list_edges_var,list_edges_var,list_nodes_var);
					x[j]=1-x[j];
					/*CalcMaxCutGoalf(&i,x,sum_sat,edge_wt,iv,
							beg_list_edges_var,
							list_edges_var,
							beg_list_minus_edge_var,
							list_minus_edge_var,
							beg_list_plus_var_edge,
							list_plus_var_edge,
							beg_list_minus_var_edge,
							list_minus_var_edge);
					if (i!=f)
					{
						i=i;
					}*/
					if (x[j]==1)
					{
						keyf=keyf+chash[j];
					}
					else keyf=keyf-chash[j];
	
					gains[j]=-gains[j];
					if (f>=fzl)goto LEND;
					step++;
				}
			}
			/*for (kv=end_move+tabu; kv<nmoves; kv++)
			{
				move=moves[kv];
				if (gains[move]<0)
				{
					cg=cg+gains[move];
					f=f+gains[move];
					i=f_ind[move]; j=s_ind[move];
					k=x[i]; x[i]=x[j]; x[j]=k;
					CalcQAPGains(gains,x,a,b);
					if (f<=fz)goto LEND;
					step++;
				}
			}*/
			if (cg>0)
			{
				g=g+cg;
				f=f+cg;
			}
			else break;
		}
		if (g>=gmax)
		{
			for (jv=0; jv<nVar; jv++)
			{
				i=iv[jv];
				best_gains[i]=gains[i];
				xbest[i]=x[i];
			}
			keybest=keyf;
			fbest=f;
			if (f>=fzl)goto LEND;
			if (g>gmax)
			{
				last_best_move=step;
				goto AAA;
			}
			gmax=g;
		}
		gain=-999999999;
		for (jv=0; jv<end_move; jv++)
		{
			i=moves[jv];
			if (gain<gains[i])
			{
				gain=gains[i];
				j=i; kv=jv;
			}
		}
		g=g+gain;
		f=f+gain;
		ReCalcMaxCutGains(j,gains,x,edge_wt,iv,
								beg_list_edges_var,list_edges_var,list_nodes_var);
		x[j]=1-x[j];
		/*CalcMaxCutGoalf(&i,x,sum_sat,edge_wt,iv,
				beg_list_edges_var,
				list_edges_var,
				beg_list_minus_edge_var,
				list_minus_edge_var,
				beg_list_plus_var_edge,
				list_plus_var_edge,
				beg_list_minus_var_edge,
				list_minus_var_edge);
		if (i!=f)
		{
			i=i;
		}*/
		if (x[j]==1)
		{
			keyf=keyf+chash[j];
		}
		else keyf=keyf-chash[j];
		gains[j]=-gains[j];
		/*for (jv=0; jv<neqn; jv++)
		{
			k=iv[jv];
			if (k==j)continue;
			CalcMaxCutGains(k,&l,x,sum_sat,edge_wt,iv,
				beg_list_edges_var,
				list_edges_var,
						beg_list_minus_edge_var,
						list_minus_edge_var,
						beg_list_plus_var_edge,
						list_plus_var_edge,
						beg_list_minus_var_edge,
						list_minus_var_edge);
			gains[k]=l;
		}*/
		moves[kv]=moves[end_move]; moves[end_move]=j;
		end_move--;
		step++;
		if (g>=gmax)
		{
			if (g>gmax)last_best_move=step;
			gmax=g;
			fbest=f;
			for (jv=0; jv<nVar; jv++)
			{
				i=iv[jv];
				best_gains[i]=gains[i];
				xbest[i]=x[i];
			}
			keybest=keyf;
			if (f>=fzl)goto LEND;
		}
		if (step-last_best_move>50)//100
		{
			break;
		}
		
	}
	for (jv=0; jv<nVar; jv++)
	{
		i=iv[jv];
		x[i]=xbest[i];
		gains[i]=best_gains[i];
	}
	keyf=keybest;
	f=fbest;
	if (gmax>0)
	{
		goto AAA;
	}

LEND:;

	*gg=f;
	*key=keyf;
}
void MaxCutRandom_Tabu(long *gg,long nmoves,long *moves,long *x,long *last_used,
			long *key,long *chash,long *numbf,long *keyar,long *vfar,long *left,long *right,
			long *iv,long *gains,long *xbest,long *best_gains,
			long *edge_wt,long *beg_list_edges_var,long *list_edges_var,long *list_nodes_var,
			long double *amyu,long *fmaxx,long *zz,long *f0max,long *f1max,
			long double *sf,long double *nf,long double *sf1,long double *nf1,
			long *key_maxx,long *gains_max)
{
	long double urand(void);

	void CalcMaxCutGains(long j,long *g,long *x,long *edge_wt,
						long *beg_list_edges_var,long *list_edges_var,long *list_nodes_var);
	void ReCalcMaxCutGains(long j,long *gains,long *x,long *edge_wt,long *iv,
						long *beg_list_edges_var,long *list_edges_var,long *list_nodes_var);
	void ReCalc_mod(long f,long *iv,long double *amyu,int npoint,
			long *fmax,long *yz,long *zz,long *f0max,long *f1max,
			long double *sf,long double *nf,long double *sf1,long double *nf1,
			long *key_max,long key,long *gains,long *gains_max);
	void ReCalc_mod9(long f,long *iv,long double *amyu,
					long *fmax,long *yz,long *zz,long *f0max,long *f1max,
					long double *sf,long double *nf,long double *sf1,long double *nf1,
					long *key_max,long key,long *gains,long *gains_max);
	long mem_solution(long key,long f,long *numbf,long *keyar,long *vfar,
						long *left,long *right );
	_int8 flg_best;	
	long i,j,f,k,mem,keyf,jv,gmax,g,gain,step,last_best_move,keybest,
		flg_stop=0,fbest,cg,tabu,n_m,fmax,key_max,numf,lv,l;

	/*CalcMaxCutGoalf(&i,x,sum_sat,edge_wt,iv,
					beg_list_edges_var,
					list_edges_var,
					beg_list_minus_edge_var,
					list_minus_edge_var,
					beg_list_plus_var_edge,
					list_plus_var_edge,
					beg_list_minus_var_edge,
					list_minus_var_edge);
	if (i!=(*gg))
	{
		i=i;
	}*/
	//max_n_iter=neqn/2;//super
	//max_n_iter=neqn;
	//max_n_iter=3*neqn;
	//max_n_iter=5*neqn;
	//max_n_iter=7*neqn;
	max_n_iter=9*neqn;
	//max_n_iter=12*neqn;
	//max_n_iter=27*neqn;

	f=*gg;
	keyf=*key;
	fmax=*fmaxx;
	key_max=*key_maxx;
	numf=*numbf;
	step=last_best_move=0; 
	for (jv=0; jv<nVar; jv++)
	{
		i=iv[jv];
		best_gains[i]=gains[i];
		xbest[i]=x[i];
	}
	keybest=keyf;
	fbest=f;
	tabu=20;
	//tabu=100;
	nmoves=neqn;

	
	
AAA:;
	last_best_move=step; 
	last_used[0]=-9999999;
	for (jv=nmoves-1; jv>0; jv--)
	{
		last_used[jv]=-9999999;
		k=(long)(urand()*jv);
		i=moves[k]; moves[k]=moves[jv]; moves[jv]=i;
	}
	g=gmax=0;

	while(1)
	{
		gain=-999999999; flg_best=0; n_m=0;
		for (jv=0; jv<nmoves; jv++)
		{
			i=moves[jv];
			if (step-last_used[i]<tabu)
			{
				if (f+gains[i]<=total_fmax)continue;
				flg_best=1;
			}
			//iv[n_m]=i; n_m++; weight[i]=-gains[i];
			if (gain<gains[i])
			{
				gain=gains[i]; j=i; 
			}
		}
		if (gain<=0)
		{
			//quick_sort_l(iv,iv+n_m-1,weight);
			//if (step-last_best_move>3)j=0;
			//else j=(long)(urand()*7);
			//j=(long)(urand()*7);
			//j=0;
			//move=iv[j]; 
			g=g+gain;
			f=f+gain;
			ReCalcMaxCutGains(j,gains,x,edge_wt,iv,
								beg_list_edges_var,list_edges_var,list_nodes_var);
			x[j]=1-x[j];
			if (x[j]==1)
			{
				keyf=keyf+chash[j];
			}
			else keyf=keyf-chash[j];
			gains[j]=-gains[j];
//			for (lv=0; lv<neqn; lv++)
//			{
//				CalcMaxCutGains(lv,&l,x,edge_wt,beg_list_edges_var,list_edges_var,list_nodes_var);
//				if (l!=gains[lv])
//				{
//					l=l;
//				}
//			}
			last_used[j]=step; step++;
			if (step-last_best_move>max_n_iter)break;
			continue;
		}
		if (flg_best==1)
		{
			g=g+gain;
			f=f+gain;
			ReCalcMaxCutGains(j,gains,x,edge_wt,iv,
								beg_list_edges_var,list_edges_var,list_nodes_var);
			x[j]=1-x[j];
			if (x[j]==1)
			{
				keyf=keyf+chash[j];
			}
			else keyf=keyf-chash[j];
			gains[j]=-gains[j];
			last_used[j]=step; step++;
		}

		while(1)
		{
			for (jv=nmoves-1; jv>0; jv--)
			{
				k=(long)(urand()*jv);
				i=moves[k]; moves[k]=moves[jv]; moves[jv]=i;
			}
			cg=0;
			for (jv=0; jv<nmoves; jv++)
			{
				j=moves[jv];
				if (step-last_used[j]<tabu)
				{
					if (f+gains[j]<=record_f)continue;
				}
				if (gains[j]>=0)
				{
					cg=cg+gains[j];
					ReCalcMaxCutGains(j,gains,x,edge_wt,iv,
										beg_list_edges_var,list_edges_var,list_nodes_var);

					x[j]=1-x[j];
					if (x[j]==1)
					{
						keyf=keyf+chash[j];
					}
					else keyf=keyf-chash[j];
					gains[j]=-gains[j];
//					for (lv=0; lv<neqn; lv++)
//					{
//						if (lv==721)
//						{
//							lv=lv;
//						}
//						CalcMaxCutGains(lv,&l,x,edge_wt,beg_list_edges_var,list_edges_var,list_nodes_var);
//						if (l!=gains[lv])
//						{
//							l=l;
//						}
//					}
					last_used[j]=step; step++;
				}
			}
			if (cg>0)
			{
				f=f+cg;
				g=g+cg;
			}
			else break;
		}
		/*if (f>(long)(0.9*total_fmax))
		{
			mem=mem_solution(keyf,f,&numf,keyar,vfar,left,right);
			if (mem==0)
			{
				ReCalc_mod(f,iv,amyu,&fmax,x,zz,f0max,f1max,
						sf,nf,sf1,nf1,&key_max,keyf,gains,gains_max,
						sum_sat,sum_sat_max);
			}
		}*/
		if (g>=gmax)
		{
			for (jv=0; jv<nVar; jv++)
			{
				i=iv[jv];
				best_gains[i]=gains[i];
				xbest[i]=x[i];
			}
			keybest=keyf;
			fbest=f;
			//if (f>(long)(0.7*total_fmax))
			{
				//mem=mem_solution(keyf,f,&numf,keyar,vfar,left,right);
				//if (mem==0)
				{
					ReCalc_mod9(f,iv,amyu,&fmax,x,zz,f0max,f1max,
							sf,nf,sf1,nf1,&key_max,keyf,gains,gains_max);
				}
			}
			if (f>=fzl)goto LEND;
			if (g>gmax)
			{
				last_best_move=step;
				goto AAA;
			}
			gmax=g;
		}
		if (step-last_best_move>max_n_iter)//neqn*30)//729//777//1000
		{
			break;
		}
	}
	for (jv=0; jv<nVar; jv++)
	{
		i=iv[jv];
		x[i]=xbest[i];
		gains[i]=best_gains[i];
	}
	keyf=keybest;
	f=fbest;
	if (gmax>0)
	{
		goto AAA;
	}

LEND:;

	*gg=f;
	*key=keyf;
	*fmaxx=fmax;
	*key_maxx=key_max;
	*numbf=numf;
}
void MaxCutRandom_Tabu1Mod(long *gg,long nmoves,long *moves,long *x,long *last_used,
			long *key,long *chash,long *numbf,long *keyar,long *vfar,long *left,long *right,
			long *iv,long *gains,long *xbest,long *best_gains,
			long *edge_wt,long *beg_list_edges_var,long *list_edges_var,long *list_nodes_var,
			long double *amyu,long *fmaxx,long *zz,long *f0max,long *f1max,
			long double *sf,long double *nf,long double *sf1,long double *nf1,
			long *key_maxx,long *gains_max)
{
	long double urand(void);

	void CalcMaxCutGains(long j,long *g,long *x,long *edge_wt,
						long *beg_list_edges_var,long *list_edges_var,long *list_nodes_var);
	void ReCalcMaxCutGains(long j,long *gains,long *x,long *edge_wt,long *iv,
						long *beg_list_edges_var,long *list_edges_var,long *list_nodes_var);
	void ReCalc_mod(long f,long *iv,long double *amyu,int npoint,
			long *fmax,long *yz,long *zz,long *f0max,long *f1max,
			long double *sf,long double *nf,long double *sf1,long double *nf1,
			long *key_max,long key,long *gains,long *gains_max);
	void ReCalc_mod9(long f,long *iv,long double *amyu,
					long *fmax,long *yz,long *zz,long *f0max,long *f1max,
					long double *sf,long double *nf,long double *sf1,long double *nf1,
					long *key_max,long key,long *gains,long *gains_max);
	long mem_solution(long key,long f,long *numbf,long *keyar,long *vfar,
						long *left,long *right );
	_int8 flg_best;	
	long i,j,f,k,mem,keyf,jv,gmax,g,gain,step,last_best_move,keybest,
		flg_stop=0,fbest,cg,tabu,n_m,fmax,key_max,numf,lv,l,npoint;

	npoint=n_point;

	/*CalcMaxCutGoalf(&i,x,sum_sat,edge_wt,iv,
					beg_list_edges_var,
					list_edges_var,
					beg_list_minus_edge_var,
					list_minus_edge_var,
					beg_list_plus_var_edge,
					list_plus_var_edge,
					beg_list_minus_var_edge,
					list_minus_var_edge);
	if (i!=(*gg))
	{
		i=i;
	}*/
	//max_n_iter=neqn/7;//super
	//max_n_iter=neqn;

	//max_n_iter=neqn/3;//100;//3*neqn;
	t_tabu=tabu=tabumin;//9;//21;//21;
	max_n_iter=neqn/10;//10;//
	//if (nmy>=20 && *fmaxx==record_f && nfail==1)
	{
		//max_n_iter=neqn/7; t_tabu=20;
	}
	//max_n_iter=7*neqn;
	//max_n_iter=9*neqn;
	//max_n_iter=12*neqn;
	//max_n_iter=27*neqn; 

	f=*gg;
	keyf=*key;
	fmax=*fmaxx;
	key_max=*key_maxx;
	numf=*numbf;
	step=last_best_move=0; 
	for (jv=0; jv<nVar; jv++)
	{
		i=iv[jv];
		best_gains[i]=gains[i];
		xbest[i]=x[i];
	}
	keybest=keyf;
	fbest=f;
	nmoves=neqn; n_att=0;

	
	
AAA:;
	if (f>=(long)(0.99*total_fmax))
	{
		mem=mem_solution(keyf,f,&numf,keyar,vfar,left,right);
		if (mem==0)
		{
			ReCalc_mod(f,iv,amyu,npoint,&fmax,x,zz,f0max,f1max,
						sf,nf,sf1,nf1,&key_max,keyf,gains,gains_max);
		}
	}
	last_best_move=step; 
	last_used[0]=-999999999;
	for (jv=nmoves-1; jv>0; jv--)
	{
		last_used[jv]=-999999999;
		k=(long)(urand()*jv);
		i=moves[k]; moves[k]=moves[jv]; moves[jv]=i;
	}
	g=gmax=0; 

	while(1)
	{
		while(1)
		{
			for (jv=nmoves-1; jv>0; jv--)
			{
				k=(long)(urand()*jv);
				i=moves[k]; moves[k]=moves[jv]; moves[jv]=i;
			}
			cg=0;
			for (jv=0; jv<nmoves; jv++)
			{
				j=moves[jv];

				if (x[j]==zz[j])tabu=tabu_ar[j];
				else tabu=tabuf_ar[j];

				if (step-last_used[j]<tabu)
				{
					if (f+gains[j]<=record_f)continue;
				}
				if (gains[j]>=0)
				{
					cg=cg+gains[j];
					ReCalcMaxCutGains(j,gains,x,edge_wt,iv,
										beg_list_edges_var,list_edges_var,list_nodes_var);

					x[j]=1-x[j];
					if (x[j]==1)
					{
						keyf=keyf+chash[j];
					}
					else keyf=keyf-chash[j];
					if (keyf<0.0)
					{
						keyf=keyf;
					}
					gains[j]=-gains[j];
					last_used[j]=step; step++;
				}
			}
			if (cg>0)
			{
				f=f+cg;
				g=g+cg;
			}
			else break;
		}
		/*CalcMaxCutGoalf(&i,x,edge_wt,iv,beg_list_edges_var,list_edges_var,list_nodes_var);
		if (i!=f)
		{
			printf("ERROR in value f\n");
		}*/
		if (g>=gmax)
		{
			for (jv=0; jv<nVar; jv++)
			{
				i=iv[jv];
				best_gains[i]=gains[i];
				xbest[i]=x[i];
			}
			keybest=keyf;
			fbest=f;
			
			if (f>=fzl)goto LEND;
			if (g>gmax)
			{
				goto AAA;
			}
			gmax=g;
		}

		gain=-999999999; flg_best=0; n_m=0; j=-1;
		for (jv=0; jv<nmoves; jv++)
		{
			i=moves[jv];

			//if (x[i]==zz[i])tabu=tabu_ar[i];
			//else tabu=tabumin;
			if (x[j]==zz[j])tabu=tabu_ar[j];
			else tabu=tabuf_ar[j];


			if (step-last_used[i]<tabu)
			{
				if (f+gains[i]<=total_fmax)continue;
				flg_best=1;
			}
			if (gain<gains[i])
			{
				gain=gains[i]; j=i; 
			}
		}
		if (j<0)
		{
			j=j;
			goto TE;
		}
		gain=gain;
		{
			g=g+gain;
			f=f+gain;

			ReCalcMaxCutGains(j,gains,x,edge_wt,iv,
								beg_list_edges_var,list_edges_var,list_nodes_var);
			x[j]=1-x[j];
			if (x[j]==1)
			{
				keyf=keyf+chash[j];
			}
			else keyf=keyf-chash[j];
			if (keyf<0.0)
			{
				keyf=keyf;
			}
			gains[j]=-gains[j];

			last_used[j]=step; step++;
			if (step-last_best_move>max_n_iter)break;
			continue;
		}
		
		if (step-last_best_move>max_n_iter)//neqn*30)//729//777//1000
		{
			break;
		}
	}
TE:;
	for (jv=0; jv<nVar; jv++)
	{
		i=iv[jv];
		x[i]=xbest[i];
		gains[i]=best_gains[i];
	}
	keyf=keybest;
	f=fbest;
	if (gmax>0)
	{
		n_att=0;
		goto AAA;
	}
	else
	{
		n_att++;
		if ((n_att<9 && f>=record_f) || n_att<3)goto AAA;//9
	}

LEND:;

	*gg=f;
	*key=keyf;
	*fmaxx=fmax;
	*key_maxx=key_max;
	*numbf=numf;
}
void MaxCutRandom_Tabu1(long *gg,long nmoves,long *moves,long *x,long *last_used,
			long *key,long *chash,long *numbf,long *keyar,long *vfar,long *left,long *right,
			long *iv,long *gains,long *xbest,long *best_gains,
			long *edge_wt,long *beg_list_edges_var,long *list_edges_var,long *list_nodes_var,
			long double *amyu,long *fmaxx,long *zz,long *f0max,long *f1max,
			long double *sf,long double *nf,long double *sf1,long double *nf1,
			long *key_maxx,long *gains_max)
{
	long double urand(void);

	void CalcMaxCutGains(long j,long *g,long *x,long *edge_wt,
						long *beg_list_edges_var,long *list_edges_var,long *list_nodes_var);
	void ReCalcMaxCutGains(long j,long *gains,long *x,long *edge_wt,long *iv,
						long *beg_list_edges_var,long *list_edges_var,long *list_nodes_var);
	void ReCalc_mod(long f,long *iv,long double *amyu,int npoint,
			long *fmax,long *yz,long *zz,long *f0max,long *f1max,
			long double *sf,long double *nf,long double *sf1,long double *nf1,
			long *key_max,long key,long *gains,long *gains_max);
	void ReCalc_mod9(long f,long *iv,long double *amyu,
					long *fmax,long *yz,long *zz,long *f0max,long *f1max,
					long double *sf,long double *nf,long double *sf1,long double *nf1,
					long *key_max,long key,long *gains,long *gains_max);
	long mem_solution(long key,long f,long *numbf,long *keyar,long *vfar,
						long *left,long *right );
	_int8 flg_best;	
	long i,j,f,k,mem,keyf,jv,gmax,g,gain,step,last_best_move,keybest,
		flg_stop=0,fbest,cg,tabu,n_m,fmax,key_max,numf,lv,l,npoint;

	npoint=n_point;

	/*CalcMaxCutGoalf(&i,x,sum_sat,edge_wt,iv,
					beg_list_edges_var,
					list_edges_var,
					beg_list_minus_edge_var,
					list_minus_edge_var,
					beg_list_plus_var_edge,
					list_plus_var_edge,
					beg_list_minus_var_edge,
					list_minus_var_edge);
	if (i!=(*gg))
	{
		i=i;
	}*/
	//max_n_iter=neqn/7;//super
	//max_n_iter=neqn;

	//max_n_iter=neqn/3;//100;//3*neqn;
	t_tabu=tabu=21;//9;//21;//21;
	max_n_iter=neqn/10;//100;//
	//if (nmy>=20 && *fmaxx==record_f && nfail==1)
	{
		//max_n_iter=neqn/7; t_tabu=20;
	}
	//max_n_iter=7*neqn;
	//max_n_iter=9*neqn;
	//max_n_iter=12*neqn;
	//max_n_iter=27*neqn; 

	f=*gg;
	keyf=*key;
	fmax=*fmaxx;
	key_max=*key_maxx;
	numf=*numbf;
	step=last_best_move=0; 
	for (jv=0; jv<nVar; jv++)
	{
		i=iv[jv];
		best_gains[i]=gains[i];
		xbest[i]=x[i];
	}
	keybest=keyf;
	fbest=f;
	nmoves=neqn; n_att=0;

	
	
AAA:;
	if (f>=(long)(0.99*total_fmax))
	{
		mem=mem_solution(keyf,f,&numf,keyar,vfar,left,right);
		if (mem==0)
		{
			ReCalc_mod(f,iv,amyu,npoint,&fmax,x,zz,f0max,f1max,
						sf,nf,sf1,nf1,&key_max,keyf,gains,gains_max);
		}
	}
	last_best_move=step; 
	last_used[0]=-9999999;
	for (jv=nmoves-1; jv>0; jv--)
	{
		last_used[jv]=-999999999;
		k=(long)(urand()*jv);
		i=moves[k]; moves[k]=moves[jv]; moves[jv]=i;
	}
	g=gmax=0; 

	while(1)
	{
		while(1)
		{
			for (jv=nmoves-1; jv>0; jv--)
			{
				k=(long)(urand()*jv);
				i=moves[k]; moves[k]=moves[jv]; moves[jv]=i;
			}
			cg=0;
			for (jv=0; jv<nmoves; jv++)
			{
				j=moves[jv];
				if (step-last_used[j]<tabu)
				{
					if (f+gains[j]<=record_f)continue;
				}
				if (gains[j]>=0)
				{
					cg=cg+gains[j];
					ReCalcMaxCutGains(j,gains,x,edge_wt,iv,
										beg_list_edges_var,list_edges_var,list_nodes_var);

					x[j]=1-x[j];
					if (x[j]==1)
					{
						keyf=keyf+chash[j];
					}
					else keyf=keyf-chash[j];
					if (keyf<0.0)
					{
						keyf=keyf;
					}
					gains[j]=-gains[j];
					last_used[j]=step; step++;
				}
			}
			if (cg>0)
			{
				f=f+cg;
				g=g+cg;
			}
			else break;
		}
		/*CalcMaxCutGoalf(&i,x,edge_wt,iv,beg_list_edges_var,list_edges_var,list_nodes_var);
		if (i!=f)
		{
			printf("ERROR in value f\n");
		}*/
		if (g>=gmax)
		{
			for (jv=0; jv<nVar; jv++)
			{
				i=iv[jv];
				best_gains[i]=gains[i];
				xbest[i]=x[i];
			}
			keybest=keyf;
			fbest=f;
			
			if (f>=fzl)goto LEND;
			if (g>gmax)
			{
				goto AAA;
			}
			gmax=g;
		}

		gain=-999999999; flg_best=0; n_m=0;
		for (jv=0; jv<nmoves; jv++)
		{
			i=moves[jv];
			if (step-last_used[i]<tabu)
			{
				if (f+gains[i]<=total_fmax)continue;
				flg_best=1;
			}
			if (gain<gains[i])
			{
				gain=gains[i]; j=i; 
			}
		}
		//if (gain<=0)
		gain=gain;
		{
			g=g+gain;
			f=f+gain;
			ReCalcMaxCutGains(j,gains,x,edge_wt,iv,
								beg_list_edges_var,list_edges_var,list_nodes_var);
			x[j]=1-x[j];
			if (x[j]==1)
			{
				keyf=keyf+chash[j];
			}
			else keyf=keyf-chash[j];
			if (keyf<0.0)
			{
				keyf=keyf;
			}
			gains[j]=-gains[j];

			last_used[j]=step; step++;
			if (step-last_best_move>max_n_iter)break;
			continue;
		}
		
		if (step-last_best_move>max_n_iter)//neqn*30)//729//777//1000
		{
			break;
		}
	}
	for (jv=0; jv<nVar; jv++)
	{
		i=iv[jv];
		x[i]=xbest[i];
		gains[i]=best_gains[i];
	}
	keyf=keybest;
	f=fbest;
	if (gmax>0)
	{
		n_att=0;
		goto AAA;
	}
	else
	{
		n_att++;
		if ((n_att<9 && f>=record_f) || n_att<3)goto AAA;//9
	}

LEND:;

	*gg=f;
	*key=keyf;
	*fmaxx=fmax;
	*key_maxx=key_max;
	*numbf=numf;
}
#if !defined(CLOCKS_PER_SEC) && defined(CLK_TCK)
#define CLOCKS_PER_SEC CLK_TCK
#endif
double cpu_time( void ) 
{
    clock_t t; static clock_t last;
    //TEST_MESSAGE("clock()");
    t = clock();
    if ( last == 0 ) last = t;
    return (double)(t-last)/CLOCKS_PER_SEC;
}
