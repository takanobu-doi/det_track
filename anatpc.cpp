#include "anatpc.h"
#include "para.h"

void clear_hough_cnt(unsigned int*** hough_cnt){
  int i,j,k;
  for(i=0; i<N_AC; i++){
    for(j=0; j<DIV_HOUGH_X; j++){
      for(k=0; k<DIV_HOUGH_Y; k++){
	hough_cnt[i][j][k]=0;
      }
    }
  }
}


int hough_tra(unsigned int*** tpc_data, unsigned int*** hough_cnt){
  int hit_cnt=0;
  double  y_ratio = 1024.0/DIV_HOUGH_Y;
  
  int i,j,k,n;
  double track_x, track_y;
  double hough_x, hough_y;
  int hough_y_i;
  
  for(i=0; i<N_AC; i++){
    for(j=0; j<N_TCLK; j++){
      for(k=0; k<N_STRP; k++){
	if(tpc_data[i][j][k]==1){
	  track_x = k;
	  track_y = j;
	  hit_cnt++;
	  for(n=0; n<DIV_HOUGH_X; n++){
	    hough_x = n*2*PI/DIV_HOUGH_X;
	    hough_y = track_x * cos(hough_x) + track_y * sin(hough_x);
	    hough_y_i = (int)(hough_y/y_ratio);
	    if(hough_y_i>=0 && hough_y_i<DIV_HOUGH_Y){
	      hough_cnt[i][n][hough_y_i]++;
	    }
	  }
	}
      }
    }
  }
  
  return hit_cnt;
}

void find_hough_max(unsigned int*** hough_cnt,
		    unsigned int hough_max_pos[N_AC][2],
		    unsigned int hough_max_cnt[N_AC]){
  int i,j,k;
  unsigned int tmp_max_cnt[N_AC]={0,0};
  unsigned int tmp_max_pos[N_AC][2];
  
  for(i=0; i<N_AC; i++){
    for(j=1; j<DIV_HOUGH_X; j++){
      for(k=1; k<DIV_HOUGH_Y; k++){
	if(hough_cnt[i][j][k] > tmp_max_cnt[i]){
	  tmp_max_cnt[i] = hough_cnt[i][j][k];
	  tmp_max_pos[i][0] = j;
	  tmp_max_pos[i][1] = k;	  
	}
      }
    }
  }

  for(i=0; i<N_AC; i++){
    hough_max_pos[i][0] = tmp_max_pos[i][0];
    hough_max_pos[i][1] = tmp_max_pos[i][1];    
    hough_max_cnt[i] = tmp_max_cnt[i];
  }
  
}

int calc_pulse_integ_fadc(int*** fadc_data,
			  int pulse_integ[N_AC][N_FADC][FADC_MAX_PULSE],
			  int pulse_lead[N_AC][N_FADC][FADC_MAX_PULSE],
			  int pulse_width[N_AC][N_FADC][FADC_MAX_PULSE],
			  int integ_th){

  int i,j,k;
  int pulse_cnt;
  int integ_flag;
  int integ_cnt;
  int null_cnt;

  int tmp_integ[N_AC][N_FADC][FADC_MAX_PULSE];
  int tmp_lead[N_AC][N_FADC][FADC_MAX_PULSE];
  int tmp_width[N_AC][N_FADC][FADC_MAX_PULSE];

  /* initialization */
  for(i=0; i<N_AC; i++){
    for(j=0; j<N_FADC; j++){
      for(k=0; k<FADC_MAX_PULSE; k++){
        tmp_integ[i][j][k]=0;
        tmp_lead[i][j][k]=-1;
        tmp_width[i][j][k]=0;
      }
    }
  }
  
  /* temporary integration */
  for(i=0; i<N_AC; i++){
    for(j=0; j<N_FADC; j++){
      pulse_cnt=0;
      integ_flag=0;
      integ_cnt=0;
      null_cnt=0;
      for(k=0; k<N_ACLK; k++){
        if(fadc_data[i][j][k]>=integ_th){
          if(integ_cnt==0){
            integ_flag=1;
            tmp_lead[i][j][pulse_cnt]=k;
          }
          integ_cnt++;
          null_cnt=0;
          tmp_integ[i][j][pulse_cnt]+=fadc_data[i][j][k];
          tmp_width[i][j][pulse_cnt]++;
        }
        if(fadc_data[i][j][k]< integ_th-1){
          null_cnt++;
          if(integ_flag==1 && null_cnt>=2){
            pulse_cnt++;
            integ_flag=0;
          }
          integ_cnt=0;
        }
        if(pulse_cnt>=FADC_MAX_PULSE) break;
      }
    }
  }

  /*delete noise like (integ with only 1clk) pulse */
  for(i=0; i<N_AC; i++){
    for(j=0; j<N_FADC; j++){
      pulse_cnt=0;
      for(k=0; k<FADC_MAX_PULSE; k++){
        if(tmp_width[i][j][k]>1){
          pulse_integ[i][j][pulse_cnt]=tmp_integ[i][j][k];
          pulse_lead[i][j][pulse_cnt]=tmp_lead[i][j][k];
          pulse_width[i][j][pulse_cnt]=tmp_width[i][j][k];
          pulse_cnt++;
        }
        if(tmp_integ[i][j][k]==0) break;
      }
    }
  }
  return 0;
}

