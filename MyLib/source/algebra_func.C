#include <TF1.h>            // questo serve per la funzione cout
#include <fstream>

using namespace std;

void swap_double(double *a,double *b)
{
  double c=*a;
  *a=*b;
  *b=c;
}

void swap_int(int *a,int *b)
{
  int c=*a;
  *a=*b;
  *b=c;
}

//perform LU decomposition according to Doolittle algorithm
void LU_Doolittle_decomp(double *lu,int *p,int n,double *m_ext)
{
  //copy m
  double *m=(double*)malloc(sizeof(double)*n*n);
  for(int i=0;i<n*n;i++) m[i]=m_ext[i];
  
  //set permutation to 0
  for(int r=0;r<n;r++) p[r]=r;
  
  //loop on row
  for(int r=0;r<n;r++)
    {
      //find the first non empty line in column i
      int rne=r;
      while(m[rne*n+rne]==0 && rne<n) rne++;
      if(rne==n)
	{
	  fprintf(stderr,"singular matrix!\n");
	  exit(0);
	}
      
      //if first non empty line is different from current one, swap it
      if(rne!=r)
	{
	  for(int c=r;c<n;c++) swap_double(&(m[rne*n+c]),&(m[r*n+c]));
	  swap_int(&(p[r]),&(p[rne]));
	}
      
      //mark L and U
      int cl=r; //column in L is current row
      for(int rl=r+1;rl<n;rl++) lu[rl*n+cl]=m[rl*n+cl]/m[r*n+r];
      for(int cr=r;cr<n;cr++) lu[r*n+cr]=m[r*n+cr];
      
      //subtract from m
      for(int r1=r+1;r1<n;r1++)
	{
	  double l=-m[r1*n+r]/m[r*n+r];
	  for(int c=r;c<n;c++) m[r1*n+c]+=m[r*n+c]*l;
	}
    }
  
  free(m);
}

//solve using lu decomposition
void LU_solve(double *x,double *lu,int *p,int n,double *b)
{
  //solve for L
  double *y=(double*)malloc(sizeof(double)*n);
  for(int i=0;i<n;i++)
    {
      y[i]=b[p[i]];
      for(int j=0;j<i;j++) y[i]-=y[j]*lu[i*n+j];
    }
  
  //solve for U
  for(int i=n-1;i>=0;i--)
    {
      x[i]=y[i];
      for(int j=i+1;j<n;j++) x[i]-=x[j]*lu[i*n+j];
      x[i]/=lu[i*n+i];
    }
  
  free(y);
}

//invert the matrix with LU
void LU_invert(double *minv,double *m,int n)
{
  double *lu=(double*)malloc(sizeof(double)*n*n);
  int *p=(int*)malloc(sizeof(int)*n);
  double *b=(double*)malloc(sizeof(double)*n);
  double *x=(double*)malloc(sizeof(double)*n);
  
  LU_Doolittle_decomp(lu,p,n,m);
  
  for(int i=0;i<n;i++)
    {
      for(int j=0;j<n;j++) b[j]=(i==j);
      LU_solve(x,lu,p,n,b);
      for(int j=0;j<n;j++) minv[j*n+i]=x[j];
    }
  
  free(lu);
  free(p);
  free(b);
  free(x);
}

double determinant(double* f, int x){
  
  int pr, i, j, k, p, q, t;

  double d = 0;
  
  double *b=(double*)malloc(sizeof(double)*x*x);
  double *b1=(double*)malloc(sizeof(double)*(x-1)*(x-1));
  double *c=(double*)malloc(sizeof(double)*x);


  if( x == 2 ){
    
    d = 0;
    d = (f[x*0+0]*f[x*1+1]) - (f[x*0+1]*f[x*1+0]);
    
    free(b);
    free(c);

    return(d);

  }else{
    
    for( j = 0; j < x; j++){

      int r=0,s=0;

      for( p = 0; p < x; p++){	
	for( q = 0; q < x; q++){
	  
	  if( p != 0 && q != j ){
	    
	    b[x*r + s]=f[x*p + q];

	    s++;
	    
	    if( s > x-2 ){
	      
	      r++;
	      s = 0;

	    }
	  }
	  
	}// q
      }// p
      
      for(t = 0, pr = -1; t < (1+j); t++){
       
	pr = (-1)*pr;

      }


      for( i = 0; i < x-1; i++){
	for( k = 0; k < x-1; k++){
	  
	  b1[(x-1)*i+k] = b[x*i+k];
	    
	}// k
      }// i

      c[j]=pr*determinant( b1, (x-1));
     
    }// j
    
    for( j = 0, d = 0; j < x; j++){
       
      d = d + (f[x*0 + j]*c[j]);

    }
     
    free(b);
    free(b1);
    free(c);

     return(d);
     
  }// else
  
}


void print_matrix(int Nrows, int Ncol, double* cov_matrix, int iout){
  
  if(iout == 0){

  	for(int i = 0; i < Nrows; i++){
    	for(int j = 0; j < Ncol; j++){  
      	printf("%+.13f  ", cov_matrix[i*Ncol+j]);
    	}// j
    	printf("\n");
  	}// i
  
  }else if(iout == 1){

  	for(int i = 0; i < Nrows; i++){
    	for(int j = 0; j < Ncol; j++){  
      	printf("%+.5e  ", cov_matrix[i*Ncol+j]);
    	}// j
    	printf("\n");
  	}// i

  }// iout

}// print_matrix


void print_matrix_for_mathematica(int Nrows, int Ncol, double* cov_matrix, int iout){
  
  if(iout == 0){

  	printf("{");
  	for(int i = 0; i < Nrows; i++){
    	printf("{");
    	for(int j = 0; j < Ncol; j++){
      	printf("%+.13f", cov_matrix[i*Ncol+j]);
      	if( j != Ncol-1){
					printf(",");
      	}
    	}// j
    	if( i != Nrows-1){
      	printf("},\n");
    	}else if( i == Nrows-1){
      	printf("}");
    	}
  	}// i
  	printf("};\n\n");

  }else if(iout == 1){

  	printf("{");
  	for(int i = 0; i < Nrows; i++){
    	printf("{");
    	for(int j = 0; j < Ncol; j++){
      	printf("%+.13f", cov_matrix[i*Ncol+j]);
      	if( j != Ncol-1){
					printf(",");
      	}
    	}// j
    	if( i != Nrows-1){
      	printf("},\n");
    	}else if( i == Nrows-1){
      	printf("}");
    	}
  	}// i
  	printf("};\n\n");

  }// iout

}// print_matrix_for_mathematica


void read_matrix(int Nrows, int Ncol, double* matrix, std::string* file_name){

  char *f_name=(char*)malloc(sizeof(char)*(1024));

  sprintf( f_name, "%s", file_name[0].c_str());

  ifstream fr_cov (f_name);
  if (!fr_cov) cout<< "Cannot open "<< f_name <<endl;

  for(int i = 0; i < Nrows; i++){
    for(int j = 0; j < Ncol; j++){
      fr_cov >> matrix[i*Ncol+j];
    }// j
  }// i
  
  fr_cov.close();  
  free(f_name);
}// print_matrix


void two_blocks_matrix( double* matrix, int N, int Size_blk_1){

  for(int i = 0; i < N; i++){
    for(int j = 0; j < N; j++){

      if( (i >= 0 && i < Size_blk_1) && (j < 0 || j >= Size_blk_1) ){
				matrix[i*N+j] = 0;
      }

      if( (i >= Size_blk_1 && i < N) && (j < Size_blk_1 || j >= N) ){
				matrix[i*N+j] = 0;
      }

    }// j
  }// i

}// two_blocks_matrix


void block_matrix( double* matrix, int N, int Nblks, int* Size_blks){

  double *temp=(double*)calloc( 1, sizeof(double));

  for(int i = 0; i < Nblks; i++){
		two_blocks_matrix( matrix, N, *temp + Size_blks[i]);

    *temp += Size_blks[i];
  }// i

  free(temp);

}// block_matrix


