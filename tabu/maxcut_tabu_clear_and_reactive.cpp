#include <malloc.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <time.h>
#include <string.h>

FILE *file;
char path[100];
long best_known_f, n_edges, n_vertices;
int **edges;										//adjacency matrix for the problem
long *keys, *f_values, *left, *right;				//arrays for storing search history
long *last_used, *neighbourhood_values;				//the actual tabu list itself, list of f(x) for all neighbours
long *x, *best_x;									//configuration
long f, best_f, stored_f_count;						//f(x) ,f(best_x) , number of stored f in tree.
long rand_seed, step,best_f_step,n_random_steps, max_steps, tabu_size;
bool *neighbourhood;								// determines whether neighbourhood[i] can be moved.
double *hashes;

void main(){
#pragma region Variables
	long  number_f;
	long temp_vert1, temp_vert2, temp_weight;
#pragma endregion

#pragma region References
	void initRandom();
	bool runCondition();
	void tabu(long *f, long *x, long n_vertices, int** edges, long *step, long tabu_size, long *stored_f_count, bool *neighbourhood, long *neighbourhood_values, long *last_used, long *best_f, long *best_x);
	double getStandardRandom();
	long F(long *x, int** matrix);
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

	x = (long *)calloc(n_vertices, sizeof(long));
	if (x == NULL){
		printf("It is not enough free memory for x\n");
		goto NOT_ENOUGH_FREE_MEMORY;
	}

	best_x = (long *)calloc(n_vertices, sizeof(long));
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

keys = (long *)calloc(max_steps, sizeof(long));
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

hashes = (double *)calloc(n_vertices, sizeof(double));
if (hashes == NULL){
	printf("It is not enough free memory for array hashes\n");
	goto NOT_ENOUGH_FREE_MEMORY;
}

edges = (int **)calloc(n_vertices, sizeof(int *));
if (edges == NULL){
	printf("It is not enough free memory for array edges\n");
	goto NOT_ENOUGH_FREE_MEMORY;
}

for (int i = 0; i < n_vertices; i++){
	edges[i] = (int *)calloc(n_vertices, sizeof(int));
	if (edges[i] == NULL){
		printf("It is not enough free memory for array edges[i]\n");
		goto NOT_ENOUGH_FREE_MEMORY;
	}
}

#pragma endregion 

#pragma region Fetching_Data
fopen_s(&file, "d:\\data_maxcut\\G37.txt", "r");

fscanf(file, "%ld%ld", &n_vertices, &n_edges);

for (int i = 0; i < n_edges; i++){
	fscanf(file, "%ld%ld%ld", &temp_vert1, &temp_vert2, &temp_weight);
	temp_vert1--;
	temp_vert2--;
	edges[temp_vert1][temp_vert2] = temp_weight;
	edges[temp_vert2][temp_vert1] = temp_weight;
	
}

#pragma endregion

#pragma region Initialization
initRandom();

tabu_size = 21; // i dont know how big it should be.

for (int i = 0; i < n_vertices; i++){
	last_used[i] = -3*n_vertices;
}


for(int i = 0; i < n_vertices; i++){
	hashes[i] = getStandardRandom();
}

	//Generating random starting vector x;
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

	//calculate f(x) and save as best known.
	f = best_f = F(x, edges);
	step = 0;
	best_f_step = 0;
	stored_f_count = 0;
	n_random_steps = static_cast<int>(n_vertices / 10);
#pragma endregion
	
#pragma region Main_Loop

	while (runCondition()){
		tabu(&f,x,n_vertices,edges,&step,tabu_size, &stored_f_count,neighbourhood, neighbourhood_values, last_used, &best_f, best_x);
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
NOT_ENOUGH_FREE_MEMORY: 
	
	;
#pragma endregion		
}



void tabu(long *f, long *x,long n_vertices,int** edges, long *step,long tabu_size,long *stored_f_count, bool *neighbourhood,long *neighbourhood_values, long *last_used,long *best_f, long *best_x)
{
	long best_neighbour_index = -1 , best_neighbour_value = LONG_MIN;

#pragma region References
		void fetchNeighbourhood(bool* neighbours, long* last_used, long current_step, long tabu_size);
		long F(long *x, long f, long index, int** matrix);
		long F(long *x, int** matrix);
		void copyX(long* source, long* destination, int length);
#pragma endregion

	*step=*step+1;

	fetchNeighbourhood(neighbourhood, last_used, *step, tabu_size);

	for (int i = 0; i < n_vertices; i++){
		if (neighbourhood[i]){
			neighbourhood_values[i] = F(x, *f, i, edges);

			if (neighbourhood_values[i] >= best_neighbour_value)
			{
				best_neighbour_index = i;
				best_neighbour_value = neighbourhood_values[i];
			}
		}
	}
	
	if (best_neighbour_value >= *f)
	{		
		//If there is a step that can increase f, perform it, save value, and add to tabu list;
		*f = best_neighbour_value;
		x[best_neighbour_index] = !x[best_neighbour_index];
		last_used[best_neighbour_index] = *step;
		copyX(x, best_x, n_vertices);
		
		if (best_neighbour_value >= *best_f)
		{
			best_f_step = *step;
			*best_f = best_neighbour_value;
			printf("Step: %ld \t Best F: %ld \t F: %ld\n", *step, *best_f, F(x, edges));
			
			
		}
	}
	else 
	{
		//local optimum
		//need to diversificate search (make random steps)

	}


}

#pragma region Helper_Functions

bool runCondition(){
	return (step <= max_steps) ;
}

void copyX(long* source, long* destination, int length){
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
	srand(time(NULL));
	
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

void clearArray(long* arr){
	for (int i = 0; i, n_vertices; i++)
		arr[i] = 0;
}

long getKey(long *x){
	double key = 0;
	for (int i = 0; i < n_vertices; i++){
		if (x[i] == 1)
			key += hashes[i];
	}
	return key;
}

void fetchNeighbourhood(bool* neighbours, long* last_used ,long current_step, long tabu_size){
	for (int i = 0; i < n_vertices; i++){
		if (current_step - last_used[i] < tabu_size)
			neighbours[i] = false;
		else
			neighbours[i] = true;
	}
}

long F(long *x, int** matrix){
	long sum = 0;
	for (int i = 0; i < n_vertices; i++){
			for (int j = i + 1; j < n_vertices; j++){
				sum += matrix[i][j]*(x[i] ^ x[j]);
			}

			
	}
	return sum;
}

long F(long *x, long f, long index, int** matrix){
	long minus = 0, plus = 0;
	long xi = x[index];
	long nxi = !x[index];
	for (int j = 0; j < n_vertices; j++){
		minus += matrix[index][j] * (xi ^ x[j]);
		plus += matrix[index][j] * (nxi ^ x[j]);
	}

	return f - minus + plus;

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