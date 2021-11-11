halton(const aU,const dim,const grid)
{
  decl vPrime=<2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199>;
  vPrime=vPrime[:dim-1];
  decl max_prime=maxr(vPrime);
  decl vU,vU0,k;
  decl mU=new matrix[grid][dim];
  decl d,it;
  for(d=0;d<dim;d++)
    {
      vU=<0>;
      it=1;
      do{
	vU0=vU;
	for(k=1;k<vPrime[d];k++)
	  {
	    vU|=(vU0+k/(vPrime[d]^it));
	  }
	it+=1;
      }while((rows(vU)-max_prime)<grid);
      mU[][d]=vU[max_prime:(max_prime+grid)-1];
    }
  aU[0]=mU;
  return 1;  
}