#include <stdio.h>
#include <iostream>
#include <random>
#include <cmath>

int main()
{
  const int N = 10000;
  const int N_file = 1000;
  const double noise_level = 0.5;
  int tpc_data[2][1024][256];
  double wx, wy;
  double cx, cy;
  double noise;
  std::random_device seed_gen;
  std::mt19937 engine(seed_gen());
  std::normal_distribution<double> dist_w(4.,1.);
  std::normal_distribution<double> dist_cx(128.,10.);
  std::uniform_real_distribution<> dist_cy(30.,1000);
  std::uniform_real_distribution<> dist_noise(0.,2.);

  FILE *file[N/N_file];
  char fname[256];
  for(int n=0;n<N/N_file;n++){
    sprintf(fname,"BG_track-%d.dat",n);
    file[n] = fopen(fname,"w");
  }
  
  for(int n=0;n<N;n++){
    std::cout << n << std::endl;
    
    for(int i=0;i<2;i++){
      for(int j=0;j<1024;j++){
	for(int k=0;k<256;k++){
	  tpc_data[i][j][k] = 0;
	}
      }
    }
    
    wx = dist_w(engine);
    wy = dist_w(engine);
    cx = dist_cx(engine);
    cy = dist_cy(engine);
    
    for(int k=0;k<256;k++){
      noise = dist_noise(engine);
      noise = (noise>noise_level) ? noise : 0;
      for(int j=0;j<1024;j++){
	if(fabs((j-cy)/wy)<noise){
	  tpc_data[0][j][k] = 1;
	}
	if(sqrt(((j-cy)*(j-cy))/(wy*wy)+((k-cx)*(k-cx))/(wx*wx))<noise){
	  tpc_data[1][j][k] = 1;
	}
      }
    }

    for(int i=0;i<2;i++){
      for(int j=0;j<1024;j++){
	for(int k=0;k<256;k++){
	  fprintf(file[n%(N/N_file)],"%d ",tpc_data[i][j][k]);
	}
      }
    }
    fprintf(file[n%(N/N_file)],"\n");
  }

  for(int n=0;n<N%N_file;n++){
    fclose(file[n]);
  }
    
  return 0;
}
