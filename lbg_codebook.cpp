// lbg_codebook.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <map>
#include <numeric>

using namespace std;

// Global Parameters.
long double delta = 0.0001;   // Delta value for distortion calculations.
long double epsilon = 0.03;
int codebook = 8; // Size of the codebook calculated.
int p = 12; // number of values for each Capstral coeff vectors
long double tw[12] = {1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0};   // Tokhura Weights.
int universe_size = 0; //stroes the size of universe after reading the data.
long double universe[8000][12]; // will stores the universe in the universe 2-d array.
long double bucket[8][12]; //bucket will store the centroid and finally the codebook.
int assigned_bucket[8][8000]; // this will keep a mapping of samples of universe to the corresponding bucket.
int assigned_bucket_size[8] = {0,0,0,0,0,0,0,0}; //will keep track of no of samples in each bucket.


/*Reading Input data from the file.. */
void read_universe(char file_name[]){
	long double temp_var;
	int counter=0;
	FILE* file = fopen(file_name, "r");
	while (fscanf(file,"%Lf,",&temp_var)!= EOF)
    {	
		universe[counter][0] = temp_var;
		for(int i=1;i<p;i++){      // One row contains 12 C'i values.
			if(i==11){
				fscanf(file,"%Lf\n",&temp_var);
				universe[counter][i] = temp_var;
			}else{
				fscanf(file,"%Lf,",&temp_var);
				universe[counter][i] = temp_var;
			}
		}
		counter++;
    }
    fclose(file);
	universe_size = counter;
	return;
}

/*Tokhura Distance.*/
long double calc_tokhura_distance(long double cb[], long double cu[]){
	long double td = 0;
	for(int i=0;i<p;i++){
		td += (tw[i] * (cb[i]-cu[i])*(cb[i]-cu[i]));
	}
	return td;
}

/*mapping each universe vector to regions based on min tokhura distance*/
void mapping_to_bucket(int codebook_size){
	for(int i=0; i<codebook_size; i++){
		assigned_bucket_size[i]=0;
	}
	int size = universe_size;
	for(int i=0;i<size;i++){
		long double cu[12];
		//code to copy vector
		for(int ic=0; ic<p; ic++){
			cu[ic] = universe[i][ic];
		}  
		long double min_distance = LDBL_MAX;
		int index = 0;
		for(int j=0;j<codebook_size;j++){
			long double cb[12];
			//code to copy vector
			for(int ic=0; ic<p; ic++){
				cb[ic] = bucket[j][ic];
			}
			long double distance = calc_tokhura_distance(cb, cu);  // Calculate tokhura distance
			if(distance < min_distance){       // update min-distance.
				min_distance = distance;
				index = j;
			}
		}
		assigned_bucket[index][assigned_bucket_size[index]]= i; 
		assigned_bucket_size[index]++;
	}
	/*for(int i=0; i<8; i++){
		printf("%d- %d\n", i,assigned_bucket_size[i]);
	}*/
	return;
}

/*Calculates: Avg total distortion for universe & codebook */
long double calc_total_distortion(int codebook_size){
	long double tDistortion = 0;
	long double total = 0;
	for(int i=0;i<codebook_size;i++){
		int size = assigned_bucket_size[i];
		total += (long double)size;
		for(int j=0;j<size;j++){
			int index = assigned_bucket[i][j];
			tDistortion += calc_tokhura_distance(bucket[i], universe[index]);       
		}
	}
	long double average = tDistortion/total;    // Average distortion.
	return average;
}

/* Method to update regions for the next iterations in the K-Means algorithm. */
void update_bucket(int codebook_size){
	for(int i=0;i<codebook_size;i++){
		long double temp[12];
		for(int ic=0; ic<p; ic++){
			temp[ic] = 0.0;
		}
		int size = assigned_bucket_size[i];
		for(int j=0;j<size;j++){
			int index = assigned_bucket[i][j];
			for(int k=0;k<12;k++){
				temp[k] += universe[index][k];
				//printf("%Lf\n", temp[k]);
			}
		}
		for(int j=0;j<12;j++){
			bucket[i][j] = temp[j]/((long double)size);  
		}
	}
	return;
}

/* method to calculate k-means for different k value -2,4,8 */
void perform_kmeans(int k){
	long double old_dist = 0;   // Previous iteration distortion.

	//handling the centroid split starts--------------------------
	for(int i=0; i<(k/2); i++){
		for(int j=0; j<p; j++){
			bucket[(k/2)+i][j] = bucket[i][j]*(1+epsilon);
		}
	}
	//handling the centroid split ends-----------------------------

	mapping_to_bucket(k);
	long double tot_dist = calc_total_distortion(k);
	int count = 1;
	printf("For k: %d , Iteration no: %d , Old Dist: %Lf, New Dist: %Lf, Change in Distortion: %Lf\n",k, count, old_dist, tot_dist, abs(tot_dist - old_dist));

	while(abs(tot_dist - old_dist) > delta){         
		count++;
		update_bucket(k);
		mapping_to_bucket(k);
		old_dist = tot_dist;
		tot_dist = calc_total_distortion(k);
		printf("For k: %d , Iteration no: %d , Old Dist: %Lf, New Dist: %Lf, Change in Distortion: %Lf\n",k, count, old_dist, tot_dist, abs(tot_dist - old_dist));
	}
	printf("\n---------------------------*************************************************---------------------------------------\n");
	return;
}

/* main method for calling the function*/
int _tmain(int argc, _TCHAR* argv[])
{
	printf("Pls hold on.. Reading the universe. \n");
	char file_name[] = "Universe.csv"; 
	read_universe(file_name); 

	//handling for k=1----------------------------------------------------
	for(int i=0; i<p; i++){
			bucket[0][i] = 0.0;
	}
	for(int i=0;i<universe_size;i++){
		for(int k=0;k<12;k++){
			bucket[0][k] += universe[i][k];
		}
	}
	for(int k=0;k<12;k++){
		bucket[0][k] = bucket[0][k]/((long double)universe_size);
	}

	printf("Next.. Performing K-Means Algorithm ..... \n\n");
	for(int k=2; k<=8; k*=2){
		perform_kmeans(k);
	}
	
	printf("Data in each bucket from universe are ..... \n\n");
	for(int i=0; i<codebook; i++){
		printf("%d - %d\n", i, assigned_bucket_size[i] );
	}
	printf("The Final codebook looks like this ..... \n\n");
	for(int i=0; i<codebook; i++){
		for(int j=0; j<p; j++){
			printf("%Lf ", bucket[i][j]);
		}
		printf("\n");
	}

	getchar();
	return 0;
}

