#include "mkdata.h"
#include "para.h"

void clear_raw_wave(double*** raw_wave){
  int i,j,k;
  for(i=0; i<N_AC; i++){
    for(j=0; j<N_TCLK; j++){
      for(k=0; k<N_STRP; k++){
	raw_wave[i][j][k]=0;
      }
    }
  }
}

int get_upic_strp_num(double ele_end_pos){
  return (int)(N_STRP*ele_end_pos/TPC_SIZE);
}

int add_raw_wave(double ele_end_pos[4],   // input data 
		 double drift_time,       // input data drift time in ns
		 int ne,                  // input data
		 TGraph *wave,             // input wave template
		 double gain,             // input data
		 double*** raw_wave){     // output data
  
  int i,j;
  int end_strp[N_AC];
  int clk;
  double ns;
  double pulse_height;
  
  end_strp[0] = get_upic_strp_num(ele_end_pos[3]*cmTomm);
  end_strp[1] = get_upic_strp_num(ele_end_pos[1]*cmTomm);  

  for(i=0; i<N_AC; i++){
    if(end_strp[i]>=0 && end_strp[i]<N_STRP){
      for(j=0; j<N_TCLK; j++){
	ns = j*10.0;
	pulse_height = wave->Eval(ns-drift_time)*ne*gain;
	if(pulse_height<0) pulse_height=0;
	raw_wave[i][j][end_strp[i]] += pulse_height;
      }
    }
  }

  return 0;
}

int add_raw_wave2(double ele_end_pos[4],   // input data 
		  double drift_time,       // input data
		  int ne,                  // input data
		  TSpline5 *wave_spline,   // input wave template
		  double gain,             // input data
		  double*** raw_wave){     // output data
  int i,j;
  int end_strp[N_AC];
  int clk;
  double ns;
  double pulse_height;
  double tmp_height;
  int max_flag=0;
  double max_height;

  double min_height;
  int break_flag;
  min_height = wave_spline->Eval(0)*ne*gain;
  
  end_strp[0] = get_upic_strp_num(ele_end_pos[3]*cmTomm);
  end_strp[1] = get_upic_strp_num(ele_end_pos[1]*cmTomm);  

  for(i=0; i<N_AC; i++){
    break_flag=0;
    max_flag=0;
    if(end_strp[i]>=0 && end_strp[i]<N_STRP){
      for(j=0; j<N_TCLK; j++){
	ns = j*10.0;
	pulse_height = wave_spline->Eval(ns-drift_time)*ne*gain;

	if(max_flag==0 && pulse_height<tmp_height && pulse_height>min_height*2){
	  max_flag = 1;
	  max_height = tmp_height;
	}

	if(pulse_height<0) pulse_height=0;
	if(max_flag==1 && pulse_height<max_height*0.05){
	  pulse_height=0;	
	  break_flag=1;
	}
	raw_wave[i][j][end_strp[i]] += pulse_height;
	tmp_height = pulse_height;
	if(break_flag==1) break;
      }
    }
  }

  return 0;

}
void clear_fadc_data(int*** fadc_data){
  int i,j,k;
  for(i=0; i<N_AC; i++){
    for(j=0; j<N_FADC; j++){
      for(k=0; k<N_ACLK; k++){
	fadc_data[i][j][k]=0;
      }
    }
  }
}

int make_fadc_data(double*** raw_wave, int*** fadc_data){
  int i,j,k;
  int fadc_ch;
  double tmp_fadc_data[N_AC][N_FADC][N_TCLK]={};
  int clk;
  
  for(i=0; i<N_AC; i++){
    for(j=0; j<N_TCLK; j++){
      for(k=0; k<N_STRP; k++){
	fadc_ch  = (int)((k*N_FADC)/N_STRP);
	tmp_fadc_data[i][fadc_ch][j] += raw_wave[i][j][k];
      }
    }
  }

  /* sampling with FADC clock */
  for(i=0; i<N_AC; i++){
    for(j=0; j<N_FADC; j++){
      for(k=0; k<N_ACLK; k++){
	clk = (int)(k*(N_TCLK/N_ACLK));
	fadc_data[i][j][k] = tmp_fadc_data[i][j][clk]*mVToFADC;
      }
    }
  }

  return 0;
}

void clear_tpc_data(unsigned int*** tpc_data){
  int i,j,k;
  for(i=0; i<N_AC; i++){
    for(j=0; j<N_TCLK; j++){
      for(k=0; k<N_STRP; k++){
	tpc_data[i][j][k]=0;
      }
    }
  }
}


int make_tpc_data(double*** raw_wave, unsigned int*** tpc_data, double threshold){
  int i,j,k;
  int cnt=0;
  for(i=0; i<N_AC; i++){
    for(j=0; j<N_TCLK; j++){
      for(k=0; k<N_STRP; k++){
	if(raw_wave[i][j][k]>threshold){
	  tpc_data[i][j][k]=1;
	  cnt++;
	}
      }
    }
  }
  return cnt;
}

double Uniform( void ){
    return ((double)rand()+1.0)/((double)RAND_MAX+2.0);
}

double rand_normal( double mu, double sigma ){
    double z=sqrt( -2.0*log(Uniform()) ) * sin( 2.0*M_PI*Uniform() );
    return mu + sigma*z;
}

int drift_electron(double start[4],
		   double driftv, double diff_tra, double diff_long,
		   double stop[4]){

  double umTocm = 0.0001;
  
  double drift_len = start[2];
  double drift_time = drift_len/driftv*10.0;
  double sigma_tra  = diff_tra*sqrt(drift_len)*umTocm;
  double sigma_long = diff_long*sqrt(drift_len)*umTocm;  

  //  stop[0] = (drift_len+gen_long.Gaus(0, sigma_long))/driftv*100.0;
  stop[0] = (drift_len+rand_normal(0, sigma_long))/driftv*100.0;  
  if(stop[0]<0) stop[0]=0;

  //  stop[1] = start[1] + gen_tra.Gaus(0, sigma_tra);
  stop[1] = start[1] + rand_normal(0, sigma_tra);
  stop[2] = 0.0;
  //  stop[3] = start[3] + gen_tra.Gaus(0, sigma_tra);
  stop[3] = start[3] + rand_normal(0, sigma_tra);

  //  printf("%f\n", rand_normal(0, sigma_tra));
  
  return 0;
}

int drift_electron2(double start[4],
		    double driftv, double diff_tra, double diff_long,
		    TRandom3 *gen_tra, TRandom3 *gen_long,
		    double stop[4]){

  double umTocm = 0.0001;
  
  double drift_len = start[2];
  double drift_time = drift_len/driftv*10.0;
  double sigma_tra  = diff_tra*sqrt(drift_len)*umTocm;
  double sigma_long = diff_long*sqrt(drift_len)*umTocm;  

  stop[0] = (drift_len+gen_long->Gaus(0, sigma_long))/driftv*100.0;
  //  stop[0] = (drift_len+rand_normal(0, sigma_long))/driftv*100.0;  
  if(stop[0]<0) stop[0]=0;

  stop[1] = start[1] + gen_tra->Gaus(0, sigma_tra);
  //  stop[1] = start[1] + rand_normal(0, sigma_tra);
  stop[2] = 0.0;
  stop[3] = start[3] + gen_tra->Gaus(0, sigma_tra);
  //  stop[3] = start[3] + rand_normal(0, sigma_tra);
  
  //  printf("%f\n", gen_tra->Gaus(0, sigma_tra));
  return 0;
}

int add_fadc_noise(int*** fadc_data, TRandom3 *gen_noise, double noise){
  int i,j,k;
  double random_noise;
  for(i=0; i<N_AC; i++){
    for(j=0; j<N_FADC; j++){
      for(k=0; k<N_ACLK; k++){
	random_noise = gen_noise->Uniform(noise*2)-noise;
	fadc_data[i][j][k] += random_noise;
      }
    }
  }
  return 0;
}
