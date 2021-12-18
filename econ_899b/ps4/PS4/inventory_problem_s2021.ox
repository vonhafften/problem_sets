#include <oxstd.h>
#include <oxdraw.h>
#include <oxfloat.h>
#include <quadpack.h>
#import <maximize>
#include <oxprob.h>


/**********************************************************************************************/
/* Global variables */
decl alpha=2;
decl lambda=-4;
decl beta=0.99;
decl ibar=8;
decl ps=1;
decl pr=4;
decl c=1;
decl gamma=1/2;
decl mF1,mF0,mS;
/* Simulated choices and states */
decl vPhat; /* Estimted choice-probalities Sx1 */
decl vY; /* Observed choices Nx1 */
decl vSid; /* Observed state IDs Nx1 */
/* State space indices */
enum{iI,iC,iP}; /* Columns identifiers for each state variable */
/* Parameter names */
enum{ilambda}; /* Row identifiers for parameter */

/**********************************************************************************************/
/* Choice-specific value-function */
value(amV,const vEV)
{
  decl vU1=alpha*mS[][iC]-mS[][iP];
  decl vU0=alpha*mS[][iC].*(mS[][iI].>0)+lambda*(mS[][iC].>0).*(mS[][iI].==0);
  decl mV=(vU0+beta*mF0*vEV)~(vU1+beta*mF1*vEV);
  amV[0]=mV;
  return 1;
}
/* Expected value function */
emax(aEV,const vEV0)
{
  decl mV;
  value(&mV,vEV0);
  aEV[0]=log(sumr(exp(mV)))+M_EULER;
  return 1;
}
/* CCP mapping */
ccp(aEV,aP,const vP)
{
  decl vU1=alpha*mS[][iC]-mS[][iP];
  decl vU0=alpha*mS[][iC].*(mS[][iI].>0)+lambda*(mS[][iC].>0).*(mS[][iI].==0);
  decl vE1=M_EULER-log(vP);
  decl vE0=M_EULER-log(1-vP);
  decl mF=mF0.*(1-vP)+mF1.*vP;
  decl vEU=(1-vP).*(vU0+vE0)+vP.*(vU1+vE1);
  decl vEVp=invert(unit(rows(vP))-beta*mF)*vEU;
  decl mV;
  value(&mV,vEVp);
  aP[0]=exp(mV[][1])./sumr(exp(mV));
  aEV[0]=vEVp;
  return 1;
}
/**********************************************************************************************/

/**********************************************************************************************/
/* Likelihood function */
lfunc_ccp(const vP, const adFunc, const avScore, const amHessian)
{
  lambda=vP[ilambda];
  decl vCCP,vEV;
  ccp(&vEV,&vCCP,vPhat);
  vCCP=vCCP[vSid];
  decl vL=vCCP.*vY+(1-vCCP).*(1-vY);
  decl LLF;
  LLF=double(sumc(log(vL)));
  adFunc[0]=LLF;
  return 1;
}
lfunc_nfxp(const vP, const adFunc, const avScore, const amHessian)
{
  lambda=vP[ilambda];
  decl vCCP,vCCP0;
  decl it=0;
  decl eps=10^(-10);
  vCCP=vPhat;
  decl vEV=constant(0,rows(mS),1),vEV0;
  /* CCP algorithm */
  do{
    vEV0=vEV;
    vCCP0=vCCP;
    ccp(&vEV,&vCCP,vCCP0);
    it+=1;
  }while(norm(vEV-vEV0)>eps);
  vCCP=vCCP[vSid];  
  decl vL=vCCP.*vY+(1-vCCP).*(1-vY);
  decl LLF;
  LLF=double(sumc(log(vL)));
  adFunc[0]=LLF;
  return 1;
}
/**********************************************************************************************/


main()
{
  format(1000);
  decl mPi=<0.75,0.25;0.9,0.1>;
  decl mGamma=(gamma~(1-gamma))|(gamma~(1-gamma));
  println("Gamma: ",mGamma);
  decl vIgrid=range(0,ibar,c)';
  decl vPgrid=pr|ps;
  decl vCgrid=0|c;
  decl S=rows(vIgrid)*rows(vCgrid)*rows(vPgrid);
  decl mSid=(vIgrid|vIgrid|vIgrid|vIgrid)
    ~(zeros(rows(vIgrid),1)|ones(rows(vIgrid),1)|zeros(rows(vIgrid),1)|ones(rows(vIgrid),1))
    ~(zeros(2*rows(vIgrid),1)|ones(2*rows(vIgrid),1));
  mS=mSid;
  mS[][iC]=vCgrid[mS[][iC]];
  mS[][iP]=vPgrid[mS[][iP]];
  println("State space: ","%c",{"I","C","P"},mS);
  savemat("PS4_state_space.csv",range(0,S-1)'~mS,{"id","I","C","P"});
  println("State points: ",S);

  decl vU1=alpha*mS[][iC]-mS[][iP];
  decl vU0=alpha*mS[][iC].*(mS[][iI].>0)+lambda*(mS[][iC].>0).*(mS[][iI].==0);
  println("%c",{"I","C","P","U1","U0"},mS~vU1~vU0);

  decl i,j;
  mF0=new matrix[S][S];
  mF1=new matrix[S][S];
  decl vI1,vI0;
  vI1=setbounds(mS[][iI]-mS[][iC]+1,0,ibar);
  vI0=setbounds(mS[][iI]-mS[][iC],0,ibar);
  decl aStateID=new array[S];
  for(j=0;j<S;j++)
    {
      mF0[][j]=mPi[mSid[][iP]][mSid[j][iP]].*mGamma[mSid[][iC]][mSid[j][iC]].*(mS[j][iI].==vI0);
      mF1[][j]=mPi[mSid[][iP]][mSid[j][iP]].*mGamma[mSid[][iC]][mSid[j][iC]].*(mS[j][iI].==vI1);
      aStateID[j]=sprint("S",j);
    }
  savemat("PS4_transition_a0.csv",range(0,S-1)'~mF0,{"id"}~aStateID);
  savemat("PS4_transition_a1.csv",range(0,S-1)'~mF1,{"id"}~aStateID);  

  /* CCP iteration */
  decl vP=constant(1/2,S,1),vP0,vEV;
  decl it=0;
  decl eps=10^(-10);
  do{
    vP0=vP;
    ccp(&vEV,&vP,vP0);
    it+=1;
    println("norm : ",norm(vP0-vP));
  }while(norm(vP0-vP)>eps);
  
  decl mF=mF0.*(1-vP)+mF1.*vP;
  decl vF0,vF=constant(1/S,S,1);
  do{
    vF0=vF;
    vF=mF'vF0;
  }while(norm(vF0-vF)>10^(-10));
  decl vCF=cumulate(vF);
  decl mCF=cumulate(mF')';
  decl mCF0=cumulate(mF0')';
  decl mCF1=cumulate(mF1')';  
  println("%c",{"I & ","C & ","P & ","EV & ","Pr(x=1|s) & ","f(s) \\\ "},"%cf",
	  {"%9.6g & ","%9.6g & ","%9.6g & ","%9.6g & ","%9.6g & ","%9.6g \\\ "},mS~vEV~vP~vF);

  /* Simulation */
  decl Sim=5000;
  decl mSim=new matrix[Sim][columns(mS)];
  decl vT=range(0,Sim-1);
  decl s,u=ranu(1,1),sid;
  sid=sumc(u.>vCF);
  mSim[0][]=mS[sid][];
  vY=new matrix[Sim][1];
  vSid=new matrix[Sim][1];
  vSid[0]=sid;
  for(s=0;s<Sim;s++)
    {
      u=ranu(1,1);
      vY[s]=u.<vP[sid];
      u=ranu(1,1);
      sid=sumr(u.>mCF0[sid][])*(1-vY[s])+sumr(u.>mCF1[sid][])*vY[s];
      if(s<Sim-1)  {
	mSim[s+1][]=mS[sid][];      
	vSid[s+1]=sid;
      }
    }
  savemat("PS4_simdata.csv",vY~vSid~mS[vSid][],{"choice","state_id","I","C","P"});
  decl eta=10^(-3);
  decl vFhat=meanc(vSid.==range(0,S-1));
  decl vNhat=sumc(vSid.==range(0,S-1));
  vPhat=((vY'(vSid.==range(0,S-1)))./setbounds(vNhat,1,M_INF_POS))';
  vPhat=setbounds(vPhat,eta,1-eta);
  println("Stationry distribution = ","%c",{"Monte-Carlo & ","Exact \\\ "},"%cf",{"%9.3g & ","%9.3g \\\ "},vFhat'~vF);
  println("Frequency of purchase: ",meanc(vY)[0]);
  println("Prob of purchasing on sale:",meanc(selectifr(vY,mSim[][iP].==ps))[0]);

  /* Inverted value function */
  decl vP2,vEV2;
  ccp(&vEV2,&vP2,vPhat);
  println("%c",{"ID & ","I & ","C & ","P & ","$P(x)^\\ast$ & ","$\\hat{P}(s)$ & ","$EV$ & ","\\hat{EV}$ \\\ "},
	  "%cf",{"%9.3g & ","%9.3g & ","%9.3g & ","%9.3g & ","%9.3g & ","%9.3g & ","%9.3g \\\ "},range(1,S)'~mS~vP~vPhat~vEV~vEV2);
  

  /* Pseudo-likelihood */
  decl L,vParam=new matrix[1][1];
  vParam[ilambda]=lambda;
  decl vParam_true=vParam;


  /* MLE with estimated CCP */
  MaxControl(1000,1);
  MaxBFGS(lfunc_ccp,&vParam,&L,0,1);
  decl vParam_est_ccp=vParam;
   
  /* MLE with true CCP */
  vPhat=vP;
  MaxControl(1000,1);
  MaxBFGS(lfunc_ccp,&vParam,&L,0,1);
  decl vParam_true_ccp=vParam;
  
  /* Nested Fixed Point */
  decl vParam_est_nfxp=vParam_true;
  vParam_est_nfxp[ilambda]=lambda;
  MaxBFGS(lfunc_nfxp,&vParam_est_nfxp,&L,0,1);
  //MaxSimplex(lfunc_nfxp,&vParam_est_nfxp,&L,constant(.25,vParam_est_nfxp));

  println("////////////////////////////////////////////////////");
  println("Estimation results: ");
  println("%r",{"lambda"},"%c",{"true","2-step","True CCP","NFXP"},vParam_true~vParam_est_ccp~vParam_true_ccp~vParam_est_nfxp);

  
  decl vGrid_lambda=range(-10,0,1/10);
  decl vGrid_llf=<>;  
  for(i=0;i<columns(vGrid_lambda);i++)
    {
      vParam_est_nfxp[ilambda]=vGrid_lambda[i];
      lfunc_nfxp(vParam_est_nfxp,&L,0,0);
      vGrid_llf~=L/Sim;
    }
  DrawXMatrix(0,vGrid_llf,"llf",vGrid_lambda,"lambda",0,5);  
  ShowDrawWindow();
  SaveDrawWindow("grid_search_lambda.pdf");
  

}