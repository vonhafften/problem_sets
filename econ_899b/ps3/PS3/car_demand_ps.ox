#include<oxstd.h>
#import <maximize>
#import <solvenle>
#include <oxdraw.h>
#include <oxfloat.h>
#include <oxprob.h>
#import <maxsqp>
#import<blp_func_ps>

/* Global variables */
decl n,aProductID,vYear,T,mX,vPrice,vShare,mZ,aZ,mIV;
decl Sim,vDelta0,vDelta_iia,A,mG,iprint,eps1;

main()
{
  ranseed(-1);
  /***************************************************************************************/
  /* Main dataset + Panel structure */
  /***************************************************************************************/
  decl i,j,k,l;
  decl spec=1;
  decl aCharactName;
  decl mPanelCharact=loadmat(sprint("Car_demand_characteristics_spec",spec,".dta"),&aCharactName);

  /* Panel structure */
  n=rows(mPanelCharact);
  vYear=unique(mPanelCharact[][find(aCharactName,"Year")]);
  T=columns(vYear);
  aProductID=new array[T];
  
  /***************************************************************************************/
  /* Outcome variables */
  /***************************************************************************************/
  vShare=mPanelCharact[][find(aCharactName,"share")];
  vDelta_iia=mPanelCharact[][find(aCharactName,"delta_iia")];
  vDelta0=vDelta_iia;
  vPrice=mPanelCharact[][find(aCharactName,"price")];
  
  /***************************************************************************************/
  /* Characteristics and instrumentsal variables */
  /***************************************************************************************/

  /* Load characteristics */
  decl varlist={"price","dpm","hp2wt","size","turbo","trans","Year_1986","Year_1987","Year_1988","Year_1989","Year_1990","Year_1991","Year_1992","Year_1993","Year_1994",
		"Year_1995","Year_1996","Year_1997","Year_1998","Year_1999","Year_2000","Year_2001","Year_2002","Year_2003","Year_2004","Year_2005",
		"Year_2006","Year_2007","Year_2008","Year_2009","Year_2010","Year_2011","Year_2012","Year_2013","Year_2014","Year_2015","model_class_2",
		"model_class_3","model_class_4","model_class_5","cyl_2","cyl_4","cyl_6","cyl_8","drive_2","drive_3","Intercept"};
  decl exo_varlist={"dpm","hp2wt","size","turbo","trans","Year_1986","Year_1987","Year_1988","Year_1989","Year_1990","Year_1991","Year_1992","Year_1993","Year_1994",
		    "Year_1995","Year_1996","Year_1997","Year_1998","Year_1999","Year_2000","Year_2001","Year_2002","Year_2003","Year_2004","Year_2005",
		    "Year_2006","Year_2007","Year_2008","Year_2009","Year_2010","Year_2011","Year_2012","Year_2013","Year_2014","Year_2015","model_class_2",
		    "model_class_3","model_class_4","model_class_5","cyl_2","cyl_4","cyl_6","cyl_8","drive_2","drive_3","Intercept"};
  mX=mPanelCharact[][find(aCharactName,varlist)];
  println("/* Mean product characteristics */");
  println("%r",aCharactName[find(aCharactName,varlist)],meanc(mX)');

  /* Load price and differentation IVs */
  decl aIVname;
  decl mDemandIV=loadmat(sprint("Car_demand_iv_spec",spec,".dta"),&aIVname);
  decl aIVlist={"i_import","diffiv_local_0","diffiv_local_1","diffiv_local_2","diffiv_local_3","diffiv_ed_0"};
  decl mExclIV=mDemandIV[][find(aIVname,aIVlist)];
  println("/* Mean cost IV (import) and differentiation measures */");
  println("%r",aIVname[find(aIVname,aIVlist)],meanc(mExclIV)');
  mIV=mPanelCharact[][find(aCharactName,exo_varlist)]~mExclIV;

  /* Non-linear attributes */
  mZ=mPanelCharact[][find(aCharactName,"price")];

  /* Pre-compute the row IDs for each market */
  aProductID=new array[T];
  for(i=0;i<T;i++) aProductID[i]=vecindex(mPanelCharact[][find(aCharactName,"Year")],vYear[i]);
  
  /***************************************************************************************/
  /* Random coefficients */
  /***************************************************************************************/
  decl mEta=loadmat("Simulated_type_distribution.dta");
  Sim=rows(mEta); 
  println("distribution of eta: ",meanc(mEta)~sqrt(varc(mEta))~quantilec(mEta,<0,1/4,1/2,3/4,1>)');

  /* Pre-compute interaction between price and random-coefficient */
  /* Two dimenional arrays of JxSim matrices: T x Nb of variables */
  aZ=new array[T];
  for(i=0;i<T;i++)
    {
      aZ[i]=new array[columns(mZ)];
      for(j=0;j<columns(mZ);j++) (aZ[i])[j]=mZ[aProductID[i]][j].*mEta[][j]';
    }

  
  /***************************************************************************************/
  /* GMM Estimator */
  /***************************************************************************************/
  decl Q,vLParam,vXi;
  decl vParam=new matrix[columns(mZ)][1];
  decl vParam0=vParam; 
  decl step=0;
  
  /* 2SLS weighting matrix */
  A=invert(mIV'mIV);  

  println("/* Plot the iteration process */");
  /* Inversion algorithm */
  iprint=1;
  vParam[0]=0.6;
  /* Contraction mapping */
  inverse(&vDelta0, vParam,0,10^(-12));
   /* Newton */
  vDelta0=vDelta_iia;
  inverse(&vDelta0, vParam,1,10^(-12));  
  iprint=0;
  
  println("/* Grid Search */");
  decl vGrid=range(0,1,0.1);
  decl mQgrid=new matrix[rows(vParam)][columns(vGrid)];
  for(i=0;i<rows(vParam);i++)
    {
      for(j=0;j<columns(vGrid);j++)
	{
	  vParam[i]=vGrid[j];
	  gmm_obj(vParam,&Q, 0,0);
	  println("grid: ",Q~vParam');
	  if(Q!=.NaN) {
	    mQgrid[i][j]=-Q;
	  }
	  else {
	    Q=100;
	    vDelta0=vDelta_iia;
	  }
	  vParam[i]=vGrid[mincindex(mQgrid[i][]')];
	}
      DrawXMatrix(i,mQgrid[i][],"Obj",vGrid,"$\\lambda_p$");	  
    }
  ShowDrawWindow();
  SaveDrawWindow(sprint("Car_demand_grid_spec",spec,".pdf"));

  println("/* Two-step GMM */");
  do{
    vParam0=vParam;
    if(step>0) {
      mG=(vXi.*mIV); mG-=meanc(mG);
      A=invert(mG'mG);
    }
    MaxControl(100,1);

    //MaxSimplex(gmm_obj,&vParam,&Q,constant(1/10,vParam));
    MaxControl(1000,1);
    MaxBFGS(gmm_obj,&vParam,&Q,0,1);
    inverse(&vDelta0, vParam,1,10^(-12));  
    vLParam=ivreg(vDelta0,mX,mIV,A);  
    vXi=vDelta0-mX*vLParam;    
    step+=1;
    println("norm: ",norm(vParam0-vParam));
  }while(step<2);

  println("Parameter estimates: ");
  println("%r",{"price random-coefficient paramter"},vParam);
  println("%r",aCharactName[find(aCharactName,varlist)],vLParam);

}