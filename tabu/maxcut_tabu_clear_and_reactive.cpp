#include <malloc.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <time.h>
#include <string.h>

FILE *data_file, *best_known_solution_file;
char data_file_path[100],problem_name[100],file_suffix[100], path[100];
long best_known_solution, n_edges, n_vertices;
long **edges;
long *keys,*f_values,*left,*right,*edge_w, *n_con_edges;
long f_memory_count;

void main(){
	unsigned long int_size, long_double_size, dimension;
	long number_of_function_eval,number_f;
	long temp_vert1, temp_vert2, temp_weight;

	strcpy_s(path, "d:\\data_maxcut\\G37_bks.txt");
	fopen_s(&best_known_solution_file, path, "r");
	fscanf(best_known_solution_file, "%ld\n", &best_known_solution);
	fclose(best_known_solution_file);

	strcpy_s(path, "d:\\data_maxcut\\G37.txt");
	fopen_s(&data_file, path, "r");
	fscanf(data_file, "%ld%ld", &n_vertices, &n_edges);
	fclose(data_file);

#pragma region Memory_Allocation

	int_size = sizeof(int);
	long_double_size = sizeof(long double);
	number_of_function_eval = 10e7;
	dimension = number_of_function_eval;

	keys = (long *)calloc(dimension, int_size);
	if (keys == NULL){
		printf("It is not enough free memory for array keyar\n");
		goto NOT_ENOUGH_FREE_MEMORY;
	}

	f_values = (long *)calloc(dimension, int_size);
	if (f_values == NULL){
		printf("It is not enough free memory for array f_values\n");
		goto NOT_ENOUGH_FREE_MEMORY;
	}

	left = (long *)calloc(dimension, int_size);
	if (left == NULL){
		printf("It is not enough free memory for array left\n");
		goto NOT_ENOUGH_FREE_MEMORY;
	}

	right = (long *)calloc(dimension, int_size);
	if (right == NULL){
		printf("It is not enough free memory for array right\n");
		goto NOT_ENOUGH_FREE_MEMORY;
	}

	edges =  (long **)malloc(n_vertices * sizeof(long *));
	if (edges == NULL){
		printf("It is not enough free memory for array edges\n");
		goto NOT_ENOUGH_FREE_MEMORY;
	}

	for (int i = 0; i < n_vertices; i++){
		edges[i] =(long *) malloc(n_vertices * sizeof(long));
		if (edges[i] == NULL){
			printf("It is not enough free memory for array edges[i]\n");
			goto NOT_ENOUGH_FREE_MEMORY;
		}
	}



	dimension = n_edges + 1;
	edge_w = (long *)calloc(dimension, int_size);
	if (edge_w == NULL)goto NOT_ENOUGH_FREE_MEMORY;


	dimension = n_vertices; +1;
	n_con_edges = (long *)calloc(dimension, int_size);
	if (n_con_edges == NULL)goto NOT_ENOUGH_FREE_MEMORY;

#pragma endregion 

	for (int i = 0; i < n_vertices; i++){
		n_con_edges[i] = 0;
	}

	data_file = fopen("d:\\data_maxcut\\G37.txt", "r");
	fscanf(data_file, "%ld%ld", &n_vertices, &n_edges);
	for (int i = 0; i < n_edges; i++){
		fscanf(data_file, "%ld%ld%ld", &temp_vert1, &temp_vert2, &temp_weight);
		temp_vert1--;
		temp_vert2--;
		edges[temp_vert1][temp_vert2] = temp_weight;
		edges[temp_vert2][temp_vert1] = temp_weight;
		edge_w[i] = temp_weight; // 
		n_con_edges[temp_vert1] ++;
		n_con_edges[temp_vert2] ++;

	}


#pragma region NOT_ENOUGH_FREE_MEMORY
NOT_ENOUGH_FREE_MEMORY: ;

#pragma endregion		
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