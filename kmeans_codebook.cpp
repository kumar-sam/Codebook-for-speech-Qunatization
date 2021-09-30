// kmeans_codebook.cpp : Defines the entry point for the console application.
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
int codebook = 8; // Size of the codebook calculated.
int p = 12; // number of values for each Capstral coeff vectors
long double tw[12] = {1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0};   // Tokhura Weights.
int universe_size = 0; //stroes the size of universe after reading the data.
long double universe[8000][12]; // will stores the universe in the universe 2-d array.
long double bucket[8][12]; //bucket will store the centroid and finally the codebook.
int assigned_bucket[8][8000]; // this will keep a mapping of samples of universe to the corresponding bucket.
int assigned_bucket_size[8] = {0,0,0,0,0,0,0,0}; //will keep track of no of samples in each bucket.

/*Reading Input data from the CSV file.. */
void read_universe(char file_name[]){
	long double temp_var;
	int counter=0;
	FILE* file = fopen(file_name, "r");
	while (fscanf(file,"%Lf,",&temp_var)!= EOF)
    {	
		universe[counter][0] = temp_var;
		for(int i=1;i<p;i++){      
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
void mapping_to_bucket(){
	for(int i=0; i<codebook; i++){
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
		for(int j=0;j<codebook;j++){
			long double cb[12];
			//code to copy vector
			for(int ic=0; ic<p; ic++){
				cb[ic] = bucket[j][ic];
			}
			long double distance = calc_tokhura_distance(cb, cu);  // Calculate tokhura distance.
			if(distance < min_distance){       // update min-distance.
				min_distance = distance;
				index = j;
			}
		}
		assigned_bucket[index][assigned_bucket_size[index]]= i;  // Fill the assigned region for each entry.
		assigned_bucket_size[index]++;
	}
	/*for(int i=0; i<8; i++){
		printf("%d- %d\n", i,assigned_bucket_size[i]);
	}*/
	return;
}


/*Calculates: Avg total distortion for universe & codebook */
long double calc_total_distortion(){
	long double tDistortion = 0;
	long double total = 0;
	for(int i=0;i<codebook;i++){
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

/* method to update regions for the next iterations in the K-Means algorithm. */
void update_bucket(){
	for(int i=0;i<codebook;i++){
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

/*main method to perform the kmeans - can say the driver code */
void perform_kmeans(){
	long double old_dist = 0;   // Previous iteration distortion.
	mapping_to_bucket();
	long double tot_dist = calc_total_distortion();
	int count = 1;
	printf("Iteration no: %d , Old Dist: %Lf, New Dist: %Lf, Change in Distortion: %Lf\n",count, old_dist, tot_dist, abs(tot_dist - old_dist));

	while(abs(tot_dist - old_dist) > delta){
		count++;
		update_bucket();
		mapping_to_bucket();
		old_dist = tot_dist;
		tot_dist = calc_total_distortion();
		printf("Iteration no: %d , Old Dist: %Lf, New Dist: %Lf, Change in Distortion: %Lf\n",count, old_dist, tot_dist, abs(tot_dist - old_dist));
	}
	return;
}

/* main method for the code */
int _tmain(int argc, _TCHAR* argv[])
{
	printf("Pls hold on.. Reading the universe. \n");
	char file_name[] = "Universe.csv"; 
	read_universe(file_name); 

	int partition = universe_size/codebook;
	for(int i=0;i<codebook;i++){
		int index = (i*partition) +1;
		//code to copy vector
		for(int ic=0; ic<p; ic++){
			bucket[i][ic] = universe[index][ic];
		}
	}
	/*
	for(int i=0; i<universe_size; i++){
		for(int j=0; j<p; j++){
			printf("%Lf ", bucket[i][j]);
		}
		printf("\n");
	}
	*/
	printf("Next.. Performing K-Means Algorithm ..... \n\n");
	perform_kmeans();
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

