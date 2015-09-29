#include <malloc.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <time.h>
#include <string.h>

FILE *file;
char data_file_path[100],problem_name[100],file_suffix[100], path[100];
long best_known_f, n_edges, n_vertices;
long **edges;										//adjacency matrix for the problem
long *keys, *f_values, *left, *right;				//arrays for storing search history
long *last_used, *neighbourhood_values;		//the actual tabu list itself
int *x, *best_x;									//configuration
long f, best_f;										//f(x) ,f(best_x) 
long rand_seed, step, max_steps, tabu_size;
bool *neighbourhood;								// determines whether neighbourhood[i] can be moved.

void main(){
#pragma region Variables
	long  number_f;
	long temp_vert1, temp_vert2, temp_weight;
#pragma endregion

#pragma region References
	void initRandom();
	bool runCondition();
	void tabu();
	double getStandardRandom();
#pragma endregion

#pragma region Fetch_Graph_Properties

	strcpy_s(path, "d:\\data_maxcut\\G37_bks.txt");
	fopen_s(&file, path, "r");
	fscanf(file, "%ld\n", &best_known_f);
	fclose(file);

	strcpy_s(path, "d:\\data_maxcut\\G37.txt");
	fopen_s(&file, path, "r");
	fscanf(file, "%ld%ld", &n_vertices, &n_edges);
	fclose(file);

#pragma endregion

#pragma region Memory_Allocation
	max_steps = 1000000;

	x = (int *)calloc(n_vertices, sizeof(int));
	if (x == NULL){
		printf("It is not enough free memory for x\n");
		goto NOT_ENOUGH_FREE_MEMORY;
	}

	best_x = (int *)calloc(n_vertices, sizeof(int));
	if (x == NULL){
		printf("It is not enough free memory for best_x\n");
		goto NOT_ENOUGH_FREE_MEMORY;
	}

	last_used = (long *)calloc(n_vertices, sizeof(long));
	if (last_used == NULL){
		printf("It is not enough free memory for array last_used\n");
		goto NOT_ENOUGH_FREE_MEMORY;
	}

	neighbourhood = (bool *)calloc(n_vertices, sizeof(bool));
	if (neighbourhood == NULL){
		printf("It is not enough free memory for array neighbourhood\n");
		goto NOT_ENOUGH_FREE_MEMORY;
	}

	neighbourhood_values = (long *)calloc(n_vertices, sizeof(long));
	if (neighbourhood_values == NULL){
		printf("It is not enough free memory for array neighbourhood_values\n");
		goto NOT_ENOUGH_FREE_MEMORY;
	}

	keys = (long *)calloc(max_steps , sizeof(long));
	if (keys == NULL){
		printf("It is not enough free memory for array keyar\n");
		goto NOT_ENOUGH_FREE_MEMORY;
	}

	f_values = (long *)calloc(max_steps, sizeof(long));
	if (f_values == NULL){
		printf("It is not enough free memory for array f_values\n");
		goto NOT_ENOUGH_FREE_MEMORY;
	}

	left = (long *)calloc(max_steps, sizeof(long));
	if (left == NULL){
		printf("It is not enough free memory for array left\n");
		goto NOT_ENOUGH_FREE_MEMORY;
	}

	right = (long *)calloc(max_steps, sizeof(long));
	if (right == NULL){
		printf("It is not enough free memory for array right\n");
		goto NOT_ENOUGH_FREE_MEMORY;
	}

	edges =  (long **)calloc(n_vertices, sizeof(long *));
	if (edges == NULL){
		printf("It is not enough free memory for array edges\n");
		goto NOT_ENOUGH_FREE_MEMORY;
	}

	for (int i = 0; i < n_vertices; i++){
		edges[i] =(long *) calloc(n_vertices , sizeof(long));
		if (edges[i] == NULL){
			printf("It is not enough free memory for array edges[i]\n");
			goto NOT_ENOUGH_FREE_MEMORY;
		}
	}

#pragma endregion 

#pragma region Fetching_Data
	fopen_s(&file,"d:\\data_maxcut\\G37.txt", "r");

	fscanf(file, "%ld%ld", &n_vertices, &n_edges);

	for (int i = 0; i < n_edges; i++){
		fscanf(file, "%ld%ld%ld", &temp_vert1, &temp_vert2, &temp_weight);
		temp_vert1--;
		temp_vert2--;
		
		edges[temp_vert2][temp_vert1] = temp_weight;
		if (temp_vert1 < temp_vert2){
			edges[temp_vert1][temp_vert2] = temp_weight;
		}
		else {
			edges[temp_vert2][temp_vert1] = temp_weight;
		}		
	}

#pragma endregion

#pragma region Initialization
	initRandom();

	step = max_steps-5; // 0!

	tabu_size = (long) floor((double)n_vertices/2); // should be tuned

	for (int i = 0; i < n_vertices; i++){	
		last_used[i] = LONG_MIN;
	}

	for (int i = 0; i< n_vertices; i++){
		if (getStandardRandom() > 0.5){
			x[i] = 1;
			best_x[i] = 1;
		}
		else{
			x[i] = 0;
			best_x[i] = 0;
		}
	}
#pragma endregion
	
#pragma region Main_Loop

	while (runCondition()){
		tabu();
	}

#pragma endregion

#pragma region Results
	printf("Best F:%ld \n", best_f);
	
	if (best_f > best_known_f){
		if (file != NULL)
			fclose(file);

		strcpy_s(path, "d:\\data_maxcut\\G37_bks.txt");
		fopen_s(&file, path, "w");
		fprintf_s(file, "%ld", best_f);
		fclose(file);
	}
	
	if (file != NULL)
		fclose(file);

	strcpy_s(path, "d:\\data_maxcut\\results\\result.txt");
	fopen_s(&file, path, "a");
	//fprintf_s(file, "\n DateTime: ", ...);
	fprintf_s(file, "Best f :\t %ld\n", best_f);
	fclose(file);
#pragma endregion

#pragma region NOT_ENOUGH_FREE_MEMORY
	return;								// exit the program to prevent executing
NOT_ENOUGH_FREE_MEMORY: ;
#pragma endregion		
}



void tabu(	)
{
	printf("%ld\n", n_edges);
	step++;
}

#pragma region Helper_Functions

bool runCondition(){
	return (step <= max_steps);
}

void copyX(char* source, char* destination, int length){
	for (int i = 0; i < length; i++)
		destination[i] = source[i];	
}

void initRandom(){
	#define ia 843314861
	#define ib 453816693
	#define mic 1693666955
	#define mmu 1073741824
	#define psu 4.65661287307739258e-10
	rand_seed = 90853;
}

double getRandom(){
	rand_seed = rand_seed*ia;
	if (rand_seed > mic)
		rand_seed = (rand_seed - mmu) - mmu;
	rand_seed = rand_seed + ib;
	if (rand_seed < 0)
		rand_seed = (rand_seed + mmu) + mmu;
	return(rand_seed*psu);
}

double getStandardRandom(){
	return (double)rand() / (double)RAND_MAX;
}


long save_solution(long key, long f, long *numbf, long *keyar, long *vfar,
	long *left, long *right)
{
	long node, item;
	if (*numbf == 0)
	{
		*numbf = *numbf + 1;
		keyar[*numbf] = key;
		vfar[*numbf] = f;
		left[*numbf] = 0;
		right[*numbf] = 0;
	}
	else
	{
		node = 1;
		while (node)
		{
			if (keyar[node] == key && vfar[node] == f)
			{
				return node;
			}
			else
			{
				item = node;
				if (keyar[node]>key)node = left[node];
				else node = right[node];
			}
		}
		*numbf = *numbf + 1;
		keyar[*numbf] = key;
		vfar[*numbf] = f;
		left[*numbf] = right[*numbf] = 0;

		if (keyar[item]>key)
			left[item] = *numbf;
		else
			right[item] = *numbf;
	}
	return 0;
}

long check_solution(long key, long f, long *numbf, long *keyar, long *vfar,
	long *left, long *right)
{
	long node, item;
	if (*numbf>0)
	{
		node = 1;
		while (node)
		{
			if (keyar[node] == key && vfar[node] == f)
			{
				return node;
			}
			else
			{
				item = node;
				if (keyar[node]>key)node = left[node];
				else node = right[node];
			}
		}
	}
	return 0;
}

#pragma endregion