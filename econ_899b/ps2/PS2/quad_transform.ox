transform(aFunc,aDFunc,const vX, const bound, const a, const b)
{
  decl mRho,mDRho;
  if(bound==0) {
    mRho=-log(vX.^(-1)-1);
    mDRho=ones(rows(mRho),columns(mRho)).*((vX.^(-1)-1).^(-1)).*(vX.^(-2));
  }
  if(bound==1) {
    mRho=-log(1-vX)+a;
    mDRho=ones(rows(mRho),columns(mRho)).*((1-vX).^(-1));
  }
  if(bound==-1) {
    mRho=log(vX)+b;
    mDRho=ones(rows(mRho),columns(mRho)).*(vX.^(-1));
  }
  if(bound==2) {
    mRho=a+(b-a).*vX;
    mDRho=ones(rows(mRho),columns(mRho)).*(b-a);
  }
  aFunc[0]=mRho;
  aDFunc[0]=mDRho;
  return 1;
}
