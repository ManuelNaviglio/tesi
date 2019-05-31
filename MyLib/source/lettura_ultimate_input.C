#include <TF1.h>            // questo serve per la funzione cout
#include <fstream>
#include "definitions.h"

using namespace std;

void lettura_ultimate_input( double* mlight, double* mstrange, double* mcharm, double a[Nbeta][Nev+1], double ainv[Nbeta][Nev+1], double* r0, double Zev[Nbeta][Nev+1], int iboot[Nbeta][Nmusea][Nev+1], double* f0, double* B0, double* fkfpi, double ZTev[Nbeta][Nev+1]){
  
  double Z1[Nbeta], Z2[Nbeta];
  double ZT1[Nbeta], ZT2[Nbeta];
  
  ifstream fi ("Input_corr_3pts/ultimate_input.out");
  if (!fi) cout<< "Cannot open "<< "Input_corr_3pts/ultimate_input.out" <<endl;
  
  string line;
  
  ////////////////////////
  //                    //
  //      m_light       //
  //                    //
  ////////////////////////

  getline(fi, line);  // ml (GeV)
  cout<<line<<endl;  
  
  for(int iev = 1; iev <= Nev+1; iev++ ){
    fi >> mlight[iev-1];
  }
  

  
  ////////////////////////
  //                    //
  //      m_strange     //
  //                    //
  ////////////////////////
  
  getline(fi, line); // newline
  cout<<line<<endl;
  
  getline(fi, line); //  ms (GeV)
  cout<<line<<endl;
  
  for(int iev = 1; iev <= Nev+1; iev++ ){
    fi >> mstrange[iev-1];
  }
  

  
  ////////////////////////
  //                    //
  //      m_charm       //
  //                    //
  ////////////////////////
  
  getline(fi, line); // newline
  cout<<line<<endl;
  
  getline(fi, line); // mc(2 GeV) (GeV)
  cout<<line<<endl;
  
  for(int iev = 1; iev <= Nev+1; iev++ ){
    fi >> mcharm[iev-1];
  }
  
  
  ////////////////////////
  //                    //
  //      a^-1 (GeV)    //
  //                    //
  ////////////////////////
  
  getline(fi, line);  // newline
  cout<<line<<endl;
  
  getline(fi, line); // a^-1 (GeV) (1.90     1.95       2.10)
  cout<<line<<endl;
  
  
  for(int iev = 1; iev <= Nev+1; iev++ ){
    fi >> ainv[0][iev-1];
    fi >> ainv[2][iev-1];
    fi >> ainv[4][iev-1];
    ainv[1][iev-1] = ainv[0][iev-1];
    ainv[3][iev-1] = ainv[2][iev-1];
  }
  
  for(int ibeta = 0; ibeta < Nbeta; ibeta++){
    for(int iev = 1; iev <= Nev+1; iev++ ){
      
      a[ibeta][iev-1] = 1/ainv[ibeta][iev-1];
      
    }
  }
  

  
  ////////////////////////
  //                    //
  //      r0 (GeV^-1)   //
  //                    //
  ////////////////////////
  
  getline(fi, line);  // newline
  cout<<line<<endl;
  
  getline(fi, line); //  r0 (GeV^-1)
  cout<<line<<endl;
  
  
  for(int iev = 1; iev <= Nev+1; iev++ ){
    fi >> r0[iev-1];
  }
  
  
  ////////////////////////
  //                    //
  //         Zp         //
  //                    //
  ////////////////////////
  
  getline(fi, line);  // newline
  cout<<line<<endl;
  
  getline(fi, line); //  Zp (1.90     1.95       2.10)
  cout<<line<<endl;
  
  fi >> Z1[0];
  fi >> Z1[2];
  fi >> Z1[4];
  Z1[1] = Z1[0];
  Z1[3] = Z1[2];
  Zev[0][0] = Z1[0];
  Zev[1][0] = Z1[1];
  Zev[2][0] = Z1[2];
  Zev[3][0] = Z1[3];
  Zev[4][0] = Z1[4];
  
  /*
    printf("%f %f %f\n", Z1[0], Z1[2], Z1[4]);
  */
  
  for(int iev = 2; iev <= 401; iev++ ){
    fi >> Zev[0][iev-1];
    fi >> Zev[2][iev-1];
    fi >> Zev[4][iev-1];
    Zev[1][iev-1] = Zev[0][iev-1];
    Zev[3][iev-1] = Zev[2][iev-1];
  }
  
  /*
    for(int iev = 2; iev <= 401; iev++ ){
    printf("%f %f %f\n", Zev[0][iev-1], Zev[2][iev-1], Zev[4][iev-1]);
    }
  */
  
  fi >> Z2[0];
  fi >> Z2[2];
  fi >> Z2[4];
  Z2[1] = Z2[0];
  Z2[3] = Z2[2];
  
  /*
    printf("%f %f %f\n", Z2[0], Z2[2], Z2[4]);
  */
  
  for(int iev = 402; iev <= Nev+1; iev++ ){
    fi >> Zev[0][iev-1];
    fi >> Zev[2][iev-1];
    fi >> Zev[4][iev-1];
    Zev[1][iev-1] = Zev[0][iev-1];
    Zev[3][iev-1] = Zev[2][iev-1];
  }
  

  
  ////////////////////////
  //                    //
  //         iboot      //
  //                    //
  ////////////////////////
  
  getline(fi, line);  // newline
  cout<<line<<endl;
  
  getline(fi, line); // Jackknife numbers ( 0.0030(32),0.0040(32), 0.0050(32), 0.0040(24) ...)
  cout<<line<<endl;
  
  for(int iev = 2; iev <= Nev+1; iev++ ){
    fi >> iboot[0][0][iev-1]; 
    fi >> iboot[0][1][iev-1]; 
    fi >> iboot[0][2][iev-1];
    fi >> iboot[1][0][iev-1]; 
    fi >> iboot[1][1][iev-1]; 
    fi >> iboot[1][2][iev-1]; 
    fi >> iboot[1][3][iev-1];
    fi >> iboot[2][0][iev-1]; 
    fi >> iboot[2][1][iev-1]; 
    fi >> iboot[2][2][iev-1]; 
    fi >> iboot[2][3][iev-1];
    fi >> iboot[3][0][iev-1];
    fi >> iboot[4][0][iev-1]; 
    fi >> iboot[4][1][iev-1]; 
    fi >> iboot[4][2][iev-1];
  }
  
  
  ////////////////////////
  //                    //
  //      f0 (GeV)      //
  //                    //
  ////////////////////////
  
  getline(fi, line);  // newline
  cout<<line<<endl;
  
  getline(fi, line); // f0 (GeV)
  cout<<line<<endl;
  
  for(int iev = 1; iev <= Nev+1; iev++ ){
    fi >> f0[iev-1];
  }
  
  
  ////////////////////////
  //                    //
  //      2*B0 (GeV)    //
  //                    //
  ////////////////////////
  
  getline(fi, line);  // newline
  cout<<line<<endl;
  
  getline(fi, line); // 2*B0 (GeV)
  cout<<line<<endl;
  
  for(int iev = 1; iev <= Nev+1; iev++ ){
    fi >> B0[iev-1];
  }
  
  
  ////////////////////////
  //                    //
  //       fk/fpi       //
  //                    //
  ////////////////////////
  
  getline(fi, line);  // newline
  cout<<line<<endl;
  
  getline(fi, line); // fk/fpi
  cout<<line<<endl;
  
  for(int iev = 1; iev <= Nev+1; iev++ ){
    fi >> fkfpi[iev-1];
  }
  
  ////////////////////////
  //                    //
  //         ZT         //
  //                    //
  ////////////////////////
  
  getline(fi, line);  // newline
  cout<<line<<endl;
  
  getline(fi, line); //  ZT(2GeV) (1.90     1.95       2.10)
  cout<<line<<endl;
  
  fi >> ZT1[0];
  fi >> ZT1[2];
  fi >> ZT1[4];
  ZT1[1] = ZT1[0];
  ZT1[3] = ZT1[2];
  ZTev[0][0] = ZT1[0];
  ZTev[1][0] = ZT1[1];
  ZTev[2][0] = ZT1[2];
  ZTev[3][0] = ZT1[3];
  ZTev[4][0] = ZT1[4];
  
  for(int iev = 2; iev <= 401; iev++ ){
    fi >> ZTev[0][iev-1];
    fi >> ZTev[2][iev-1];
    fi >> ZTev[4][iev-1];
    ZTev[1][iev-1] = ZTev[0][iev-1];
    ZTev[3][iev-1] = ZTev[2][iev-1];
  }
  
  fi >> ZT2[0];
  fi >> ZT2[2];
  fi >> ZT2[4];
  ZT2[1] = ZT2[0];
  ZT2[3] = ZT2[2];
  
  for(int iev = 402; iev <= Nev+1; iev++ ){
    fi >> ZTev[0][iev-1];
    fi >> ZTev[2][iev-1];
    fi >> ZTev[4][iev-1];
    ZTev[1][iev-1] = ZTev[0][iev-1];
    ZTev[3][iev-1] = ZTev[2][iev-1];
  }
  
  
   fi.close();
   
   
   ////////  FINE LETTURA Ultimate Input
   
}
