#include <malloc.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <time.h>
#include <windows.h>
#include <string>


#pragma region TYPEDEFS

struct graph_s {
	long n_edges;
	long n_vertices;
	int **edges;

	int * n_con_edges;
	int * beg_list_edges_var;
	int * list_edges_var;
	int * list_nodes_var;
	int * edge_weight;
} graph;

struct history_s {
	long * f_values;
	long * left;
	long * right;
	long * occurence_time;
	double * keys;
	double * hashes;
	long stored_f_count = 0;
	long max_steps = 0;
} history;

struct statistics {

};
#pragma endregion

const double eps = 0.00001;
const int test_time_seconds = 2;
const int UR_multiplier = 150000;

#pragma region GLOBAL_VARS
FILE *file, *file2;
char path[100];
char datetimeprefix[80];
long best_known_f;

long *last_used, *neighbourhood_values;				//the actual tabu list itself, list of f(x) for all neighbours
long *x, *best_x;									//configuration
long f, best_f;					                    //f(x) ,f(best_x) 
long rand_seed;
long step, best_f_step, max_steps, tabu_size, tabu_change_t;
bool *neighbourhood;								// determines whether neighbourhood[i] can be moved.
double x_key;

long * reactive_tabu_size_statistics;
long * index_count;
long * ultra_reactive_index_count;//array[i] is count of i vertices moved
long * ultra_reactive_tabu_sizes;
long previous_solution_found_count = 0;

double * elapsed_time, *elapsed_time_reactive, *elapsed_time_ultra_reactive;
clock_t start, current, finish;
#pragma endregion


void main() 
{
#pragma region Variables
	long temp_vert1, temp_vert2, temp_weight;
	long iterations_count;
	long *iteration_bestf, *iteration_bestf_reactive;
	long *iteration_steps, *iteration_steps_reactive;

	bool isTimeStopCondition = true;
#pragma endregion

#pragma region References
	void initRandom();
	bool runCondition();
	double getKey(long *x);
	void tabu(graph_s graph, long *f, long *x, long n_vertices, int** edges, long *step, long max_steps, long *tabu_size,
		long *stored_f_count, bool *neighbourhood, long *neighbourhood_values, long *last_used, long *best_f, long *best_x, long *best_f_step);
	void reactive_tabu(graph_s graph, long *f, long *x, double *x_key, long *step, long max_steps, long *tabu_size,
		bool *neighbourhood, long *neighbourhood_values, long *last_used, long *best_f, long *best_x, long *best_f_step, long* tabu_change_t, long* reactive_tabu_size_statistics);
	void ultra_reactive_tabu(graph_s graph, long *f, long *x, double *x_key, long *step, long max_steps, long *tabu_size,
		bool *neighbourhood, long *neighbourhood_values, long *last_used, long *best_f, long *best_x, long *best_f_step, long* tabu_change_t, long* ultra_reactive_tabu_sizes, long* ultra_reactive_index_count);
	void ultra_reactive_tabu_noac(graph_s graph, long *f, long *x, double *x_key, long *step, long max_steps, long *tabu_size,
		bool *neighbourhood, long *neighbourhood_values, long *last_used, long *best_f, long *best_x, long *best_f_step, long* tabu_change_t, long* ultra_reactive_tabu_sizes, long* ultra_reactive_index_count);
	double getStandardRandom();
	long F(long *x, int** matrix);
	void tabu_size_log_reactive(long iteration_test_count, long max_steps);
	void graph_move_index_log();
	void graph_move_index_log_reactive();
	void graph_move_index_log_ultra_reactive();
	long save_solution(double key, long f, struct history_s *history, long step);
	long check_solution(double key, long f, struct history_s *history);
	char* getDateTimePrefix();
	void addDateTimePrefixToFileName(char* dest, const char* first, const char* second, char* dtprefix);
	void printResultToFile(int index, const char* alg_name);
	void appendIterationResult(int index, long F, long steps, float elapsed_time, char* filename);
	bool allocateMemory(long** x, long** best_x, long** last_used, bool** neighbourhood, long** neighbourhood_values, history_s* history, graph_s* graph,
		long** iteration_bestf, long** iteration_bestf_reactive, long iterations_count, long** iteration_steps, long** iteration_steps_reactive, long** reactive_tabu_size_statistics,
		long** index_count, long** ultra_reactive_tabu_sizes, double** elapsed_time, double** elapsed_time_reactive, double** elapsed_time_ultra_reactive);
#pragma endregion

#pragma region Fetch_Graph_Properties

	strcpy_s(path, "d:\\data_maxcut\\G37_bks.txt");
	fopen_s(&file, path, "r");
	fscanf(file, "%ld\n", &best_known_f);
	fclose(file);

	strcpy_s(path, "d:\\data_maxcut\\G37.txt");
	fopen_s(&file, path, "r");
	fscanf(file, "%ld%ld", &(graph.n_vertices), &(graph.n_edges));
	fclose(file);

#pragma endregion

	strcpy_s(datetimeprefix, getDateTimePrefix());

#pragma region CONSTANTS

	max_steps = 15000000;
	iterations_count = 16;
	tabu_size = 75;

#pragma endregion

#pragma region Memory_Allocation
	
	if (!allocateMemory(&x, &best_x, &last_used,&neighbourhood, &neighbourhood_values, &history, &graph,&iteration_bestf,&iteration_bestf_reactive, 
		iterations_count,&iteration_steps,&iteration_steps_reactive,&reactive_tabu_size_statistics,&index_count,&ultra_reactive_tabu_sizes,&elapsed_time, &elapsed_time_reactive, &elapsed_time_ultra_reactive))
		goto NOT_ENOUGH_FREE_MEMORY;

//	x = (long *)calloc(graph.n_vertices, sizeof(long));
//	if (x == NULL) {
//		printf("It is not enough free memory for x\n");
//		goto NOT_ENOUGH_FREE_MEMORY;
//	}
//
//	best_x = (long *)calloc(graph.n_vertices, sizeof(long));
//	if (x == NULL) {
//		printf("It is not enough free memory for best_x\n");
//		goto NOT_ENOUGH_FREE_MEMORY;
//	}
//
//	last_used = (long *)calloc(graph.n_vertices, sizeof(long));
//	if (last_used == NULL) {
//		printf("It is not enough free memory for array last_used\n");
//		goto NOT_ENOUGH_FREE_MEMORY;
//	}
//
//	neighbourhood = (bool *)calloc(graph.n_vertices, sizeof(bool));
//	if (neighbourhood == NULL) {
//		printf("It is not enough free memory for array neighbourhood\n");
//		goto NOT_ENOUGH_FREE_MEMORY;
//	}
//
//	neighbourhood_values = (long *)calloc(graph.n_vertices, sizeof(long));
//	if (neighbourhood_values == NULL) {
//		printf("It is not enough free memory for array neighbourhood_values\n");
//		goto NOT_ENOUGH_FREE_MEMORY;
//	}
//
//#pragma region history
//	history.keys = (double *)calloc(max_steps, sizeof(double));
//	if (history.keys == NULL) {
//		printf("It is not enough free memory for array keyar\n");
//		goto NOT_ENOUGH_FREE_MEMORY;
//	}
//
//	history.f_values = (long *)calloc(max_steps, sizeof(long));
//	if (history.f_values == NULL) {
//		printf("It is not enough free memory for array f_values\n");
//		goto NOT_ENOUGH_FREE_MEMORY;
//	}
//
//	history.left = (long *)calloc(max_steps, sizeof(long));
//	if (history.left == NULL) {
//		printf("It is not enough free memory for array left\n");
//		goto NOT_ENOUGH_FREE_MEMORY;
//	}
//
//	history.right = (long *)calloc(max_steps, sizeof(long));
//	if (history.right == NULL) {
//		printf("It is not enough free memory for array right\n");
//		goto NOT_ENOUGH_FREE_MEMORY;
//	}
//
//	history.occurence_time = (long *)calloc(max_steps, sizeof(long));
//	if (history.occurence_time == NULL) {
//		printf("It is not enough free memory for array occurence_count\n");
//		goto NOT_ENOUGH_FREE_MEMORY;
//	}
//
//	history.hashes = (double *)calloc(graph.n_vertices, sizeof(double));
//	if (history.hashes == NULL) {
//		printf("It is not enough free memory for array hashes\n");
//		goto NOT_ENOUGH_FREE_MEMORY;
//	}
//
//#pragma endregion
//
//#pragma region graph
//	graph.edges = (int **)calloc(graph.n_vertices, sizeof(int *));
//	if (graph.edges == NULL) {
//		printf("It is not enough free memory for array edges\n");
//		goto NOT_ENOUGH_FREE_MEMORY;
//	}
//
//	for (int i = 0; i < graph.n_vertices; i++) {
//		graph.edges[i] = (int *)calloc(graph.n_edges, sizeof(int));
//		if (graph.edges[i] == NULL) {
//			printf("It is not enough free memory for array edges[i]\n");
//			goto NOT_ENOUGH_FREE_MEMORY;
//		}
//	}
//
//	graph.n_con_edges = (int *)calloc(graph.n_vertices, sizeof(int));
//	if (graph.n_con_edges == NULL) {
//		printf("It is not enough free memory for array n_con_edges\n");
//		goto NOT_ENOUGH_FREE_MEMORY;
//	}
//
//	graph.beg_list_edges_var = (int *)calloc(graph.n_vertices + 1, sizeof(int));
//	if (graph.beg_list_edges_var == NULL) {
//		printf("It is not enough free memory for array beg_list_edges_var\n");
//		goto NOT_ENOUGH_FREE_MEMORY;
//	}
//
//	graph.list_edges_var = (int *)calloc(2 * graph.n_edges, sizeof(int));
//	if (graph.list_edges_var == NULL) {
//		printf("It is not enough free memory for array list_edges_var\n");
//		goto NOT_ENOUGH_FREE_MEMORY;
//	}
//
//	graph.list_nodes_var = (int *)calloc(2 * graph.n_edges, sizeof(int));
//	if (graph.list_nodes_var == NULL) {
//		printf("It is not enough free memory for array list_nodes_var\n");
//		goto NOT_ENOUGH_FREE_MEMORY;
//	}
//
//	graph.edge_weight = (int *)calloc(graph.n_edges, sizeof(int));
//	if (graph.edge_weight == NULL) {
//		printf("It is not enough free memory for array edge_weight\n");
//		goto NOT_ENOUGH_FREE_MEMORY;
//	}
//
//#pragma endregion
//
//	iteration_bestf = (long *)calloc(iterations_count, sizeof(long));
//	if (iteration_bestf == NULL) {
//		printf("It is not enough free memory for array iteration_bestf\n");
//		goto NOT_ENOUGH_FREE_MEMORY;
//	}
//
//	iteration_bestf_reactive = (long *)calloc(iterations_count, sizeof(long));
//	if (iteration_bestf_reactive == NULL) {
//		printf("It is not enough free memory for array iteration_bestf_reactive\n");
//		goto NOT_ENOUGH_FREE_MEMORY;
//	}
//
//	iteration_steps = (long *)calloc(iterations_count, sizeof(long));
//	if (iteration_steps == NULL) {
//		printf("It is not enough free memory for array iteration_steps\n");
//		goto NOT_ENOUGH_FREE_MEMORY;
//	}
//
//	iteration_steps_reactive = (long *)calloc(iterations_count, sizeof(long));
//	if (iteration_steps_reactive == NULL) {
//		printf("It is not enough free memory for array iteration_steps\n");
//		goto NOT_ENOUGH_FREE_MEMORY;
//	}
//
//	reactive_tabu_size_statistics = (long *)calloc(graph.n_vertices, sizeof(long));
//	if (reactive_tabu_size_statistics == NULL) {
//		printf("It is not enough free memory for array reactive_tabu_size_statistics\n");
//		goto NOT_ENOUGH_FREE_MEMORY;
//	}
//
//	index_count = (long *)calloc(graph.n_vertices, sizeof(long));
//	if (index_count == NULL) {
//		printf("It is not enough free memory for array reactive_tabu_size_statistics\n");
//		goto NOT_ENOUGH_FREE_MEMORY;
//	}
//
//
//	ultra_reactive_tabu_sizes = (long *)calloc(graph.n_vertices, sizeof(long));
//	if (ultra_reactive_tabu_sizes == NULL) {
//		printf("It is not enough free memory for array ultra_reactive_tabu_size_statistics\n");
//		goto NOT_ENOUGH_FREE_MEMORY;
//	}
//
//	elapsed_time = (double *)calloc(iterations_count, sizeof(long));
//	if (elapsed_time == NULL) {
//		printf("It is not enough free memory for array elapsed_time\n");
//		goto NOT_ENOUGH_FREE_MEMORY;
//	}
//
//	elapsed_time_reactive = (double *)calloc(iterations_count, sizeof(long));
//	if (elapsed_time_reactive == NULL) {
//		printf("It is not enough free memory for array elapsed_time\n");
//		goto NOT_ENOUGH_FREE_MEMORY;
//	}
//
//	elapsed_time_ultra_reactive = (double *)calloc(iterations_count, sizeof(long));
//	if (elapsed_time_ultra_reactive == NULL) {
//		printf("It is not enough free memory for array elapsed_time\n");
//		goto NOT_ENOUGH_FREE_MEMORY;
//	}
#pragma endregion 

#pragma region Fetching_Data
	fopen_s(&file, "d:\\data_maxcut\\G37.txt", "r");
	fscanf(file, "%ld%ld", &(graph.n_vertices), &(graph.n_edges));
	/*preparing graph by calculating count of adjacent vertices for every vertice.\*/
	for (int j = 0; j < graph.n_edges; j++)
	{
		fscanf(file, "%ld%ld%ld", &temp_vert1, &temp_vert2, &temp_weight);
		temp_vert1--;
		temp_vert2--;
		graph.edge_weight[j] = temp_weight;
		graph.n_con_edges[temp_vert1]++;
		graph.n_con_edges[temp_vert2]++;
		graph.edges[temp_vert1][temp_vert2] = temp_weight;
		graph.edges[temp_vert2][temp_vert1] = temp_weight;
	}

	graph.beg_list_edges_var[0] = 0;

	for (int j = 0; j < graph.n_vertices; j++)
	{
		graph.beg_list_edges_var[j + 1] = graph.beg_list_edges_var[j] + graph.n_con_edges[j];
	}

	for (int j = 0; j < graph.n_vertices; j++)
		graph.n_con_edges[j] = 0;

	if (file != NULL)
		fclose(file);

	/*Fetching nodes.*/
	fopen_s(&file, "d:\\data_maxcut\\G37.txt", "r");
	fscanf(file, "%ld%ld", &(graph.n_vertices), &(graph.n_edges));

	for (int j = 0; j < graph.n_edges; j++)
	{
		fscanf(file, "%ld%ld%ld", &temp_vert1, &temp_vert2, &temp_weight);
		temp_vert1--;
		temp_vert2--;
		graph.list_edges_var[graph.beg_list_edges_var[temp_vert1] + graph.n_con_edges[temp_vert1]] = j;
		graph.list_edges_var[graph.beg_list_edges_var[temp_vert2] + graph.n_con_edges[temp_vert2]] = j;
		graph.list_nodes_var[graph.beg_list_edges_var[temp_vert1] + graph.n_con_edges[temp_vert1]] = temp_vert2;
		graph.list_nodes_var[graph.beg_list_edges_var[temp_vert2] + graph.n_con_edges[temp_vert2]] = temp_vert1;
		graph.n_con_edges[temp_vert1]++;
		graph.n_con_edges[temp_vert2]++;
	}
	fclose(file);
#pragma endregion

/*
#pragma region testing_tabu

	printf("====================TABU====================\n");
	printf("Start tabu size: %ld\n", tabu_size);
	printf("Restarts:%d\n", iterations_count);
	if (isTimeStopCondition)
		printf("Time per iteration, s. : %ld\n", test_time_seconds);
	else
		printf("Steps per iteration : %ld\n", max_steps);

	strcpy_s(path, "d:\\data_maxcut\\results\\tabu\\log-");
	strcat_s(path, datetimeprefix);
	strcat(path, ".txt");
	fopen_s(&file, path, "a");

	fprintf_s(file, "====================TABU====================\n");
	fprintf_s(file, "Start tabu size: %ld\n", tabu_size);
	fprintf_s(file, "Restarts:%d\n", iterations_count);
	if (isTimeStopCondition)
		fprintf_s(file, "Time per iteration, s. : %ld\n", test_time_seconds);
	else
		fprintf_s(file, "Steps till restart : %ld\n", max_steps);
	fprintf_s(file, "\nIndex\tF\tSteps\tTime\n");
	fclose(file);

	for (int i = 0; i < iterations_count; i++) {
		printf("Restart %d : \t", i + 1);

#pragma region Initialization
		initRandom();

		for (int i = 0; i < graph.n_vertices; i++) {
			last_used[i] = -3 * graph.n_vertices;
		}

		for (int i = 0; i < graph.n_vertices; i++) {
			history.hashes[i] = getStandardRandom() * 100;
		}

		//Generating random starting vector x;
		for (int i = 0; i < graph.n_vertices; i++) {
			if (getStandardRandom() > 0.5) {
				x[i] = 1;
				best_x[i] = 1;
			}
			else {
				x[i] = 0;
				best_x[i] = 0;
			}
		}

		//calculate f(x) and save as best known.
		f = best_f = F(x, graph.edges);
		//TODO: best_f will be written to file so spprove with common tabu testing
		x_key = getKey(x);
		step = 0;
		best_f_step = 0;
		history.stored_f_count = 0;
		tabu_change_t = 0;
		start = clock();
		finish = start + CLOCKS_PER_SEC * (test_time_seconds);
#pragma endregion

#pragma region Main_Loop

		while (runCondition()) {
			tabu(graph, &f, x, graph.n_vertices, graph.edges, &step, max_steps, &tabu_size, &history.stored_f_count, neighbourhood, neighbourhood_values, last_used, &best_f, best_x, &best_f_step);
			//reactive_tabu(&f, x, &x_key, n_vertices, edges, &step, max_steps, &tabu_size, &stored_f_count, neighbourhood, neighbourhood_values, last_used, &best_f, best_x, &best_f_step, &tabu_change_t, keys, f_values, left, right,occurence_time );
		}
		elapsed_time[i] = (clock() - start) / CLOCKS_PER_SEC;
		iteration_steps[i] = step;
#pragma endregion

#pragma region Results
		printf("Best F:%ld\tIterations: %ld\tTime: %f\n", best_f, iteration_steps[i], elapsed_time[i]);
		iteration_bestf[i] = best_f;

		if (best_f > best_known_f) {
			best_known_f = best_f;
			if (file != NULL)
				fclose(file);

			strcpy_s(path, "d:\\data_maxcut\\G37_bks.txt");
			fopen_s(&file, path, "w");
			fprintf_s(file, "%ld", best_f);
			fclose(file);
		}

		if (file != NULL)
			fclose(file);

		printResultToFile(i, "tabu");
		addDateTimePrefixToFileName(path, "d:\\data_maxcut\\results\\tabu\\log-", ".txt", datetimeprefix);
		appendIterationResult(i + 1, best_f, iteration_steps[i], elapsed_time[i], path);
		graph_move_index_log();
#pragma endregion

#pragma region test_reallocate_memory
		//free(history.keys);
		//free(history.f_values);
		//free(history.left);
		//free(history.right);
		//free(history.occurence_time);

		history.keys = (double *)calloc(max_steps, sizeof(double));
		if (history.keys == NULL) {
			printf("It is not enough free memory for array keyar\n");
			goto NOT_ENOUGH_FREE_MEMORY;
		}

		history.f_values = (long *)calloc(max_steps, sizeof(long));
		if (history.f_values == NULL) {
			printf("It is not enough free memory for array f_values\n");
			goto NOT_ENOUGH_FREE_MEMORY;
		}

		history.left = (long *)calloc(max_steps, sizeof(long));
		if (history.left == NULL) {
			printf("It is not enough free memory for array left\n");
			goto NOT_ENOUGH_FREE_MEMORY;
		}

		history.right = (long *)calloc(max_steps, sizeof(long));
		if (history.right == NULL) {
			printf("It is not enough free memory for array right\n");
			goto NOT_ENOUGH_FREE_MEMORY;
		}

		history.occurence_time = (long *)calloc(max_steps, sizeof(long));
		if (history.occurence_time == NULL) {
			printf("It is not enough free memory for array occurence_count\n");
			goto NOT_ENOUGH_FREE_MEMORY;
		}
#pragma endregion

}

#pragma endregion
*/
/*
#pragma region testing_reactive_tabu

	printf("====================REACTIVE====================\n");
	printf("Start tabu size: %ld\n", tabu_size);
	printf("Restarts:%d\n", iterations_count);
	if (isTimeStopCondition)
		printf("Time per iteration, s. : %ld\n", test_time_seconds);
	else
		printf("Steps till restart : %ld\n", max_steps);

	strcpy_s(path, "d:\\data_maxcut\\results\\reactive_tabu\\log-");
	strcat_s(path, datetimeprefix);
	strcat(path, ".txt");

	if (file != NULL)
		fclose(file);
	fopen_s(&file, path, "a");

	fprintf_s(file, "====================REACTIVE====================\n");
	fprintf_s(file, "Start tabu size: %ld\n", tabu_size);
	fprintf_s(file, "Restarts:%d\n", iterations_count);
	if (isTimeStopCondition)
		fprintf_s(file, "Time per iteration, s. : %ld\n", test_time_seconds);
	else
		fprintf_s(file, "Steps till restart : %ld\n", max_steps);
	fprintf_s(file, "\nIndex\tF\tSteps\tTime\n");
	fclose(file);

	for (int i = 0; i < iterations_count; i++) {
		printf("Restart %d : \t", i + 1);

#pragma region Initialization
		initRandom();

		for (int i = 0; i < graph.n_vertices; i++)
			last_used[i] = -3 * graph.n_vertices;

		for (int i = 0; i < graph.n_vertices; i++)
			history.hashes[i] = getStandardRandom() * 100;

		//Generating random starting vector x;
		for (int i = 0; i < graph.n_vertices; i++) {
			if (getStandardRandom() > 0.5) {
				x[i] = 1;
				best_x[i] = 1;
			}
			else {
				x[i] = 0;
				best_x[i] = 0;
			}
		}

		//calculate f(x) and save as best known.
		f = best_f = F(x, graph.edges);
		x_key = getKey(x);
		step = 0;
		best_f_step = 0;
		history.stored_f_count = 0;
		tabu_change_t = 0;
		
		for (int i = 0; i < graph.n_vertices; i++)
			index_count[i] = 0;

		start = clock();
		finish = start + CLOCKS_PER_SEC * (test_time_seconds);
#pragma endregion

#pragma region Main_Loop

		while (runCondition()) {
			//tabu(graph, &f, x, n_vertices, edges, &step, max_steps, &tabu_size, &stored_f_count, neighbourhood, neighbourhood_values, last_used, &best_f, best_x, &best_f_step);
			reactive_tabu(&f, x, &x_key, &step, max_steps, &tabu_size, neighbourhood, neighbourhood_values, last_used, &best_f, best_x, &best_f_step, &tabu_change_t, reactive_tabu_size_statistics);
		}

#pragma endregion

#pragma region Results
		elapsed_time_reactive[i] = (clock() - start) / CLOCKS_PER_SEC;
		iteration_steps_reactive[i] = step;

		printf("Best F:%ld\tIterations: %ld\tTime: %f\n", best_f, iteration_steps_reactive[i], elapsed_time_reactive[i]);
		iteration_bestf_reactive[i] = best_f;

		if (best_f > best_known_f) {
			best_known_f = best_f;
			if (file != NULL)
				fclose(file);

			strcpy_s(path, "d:\\data_maxcut\\G37_bks.txt");
			fopen_s(&file, path, "w");
			fprintf_s(file, "%ld", best_f);
			fclose(file);
		}

		printResultToFile(i, "reactive_tabu");

		strcpy_s(path, "d:\\data_maxcut\\results\\reactive_tabu\\log-");
		strcat_s(path, datetimeprefix);
		strcat(path, ".txt");

		appendIterationResult(i + 1, best_f, iteration_steps_reactive[i], elapsed_time_reactive[i],path);

		graph_move_index_log_reactive();
		tabu_size_log_reactive(iterations_count, step);

#pragma endregion



#pragma region test_reallocate_memory
		free(history.keys);
		free(history.f_values);
		free(history.left);
		free(history.right);
		free(history.occurence_time);


		history.keys = (double *)calloc(max_steps, sizeof(double));
		if (history.keys == NULL) {
			printf("It is not enough free memory for array keyar\n");
			goto NOT_ENOUGH_FREE_MEMORY;
		}

		history.f_values = (long *)calloc(max_steps, sizeof(long));
		if (history.f_values == NULL) {
			printf("It is not enough free memory for array f_values\n");
			goto NOT_ENOUGH_FREE_MEMORY;
		}

		history.left = (long *)calloc(max_steps, sizeof(long));
		if (history.left == NULL) {
			printf("It is not enough free memory for array left\n");
			goto NOT_ENOUGH_FREE_MEMORY;
		}

		history.right = (long *)calloc(max_steps, sizeof(long));
		if (history.right == NULL) {
			printf("It is not enough free memory for array right\n");
			goto NOT_ENOUGH_FREE_MEMORY;
		}

		history.occurence_time = (long *)calloc(max_steps, sizeof(long));
		if (history.occurence_time == NULL) {
			printf("It is not enough free memory for array occurence_count\n");
			goto NOT_ENOUGH_FREE_MEMORY;
		}
#pragma endregion
	}

#pragma endregion
	*/



#pragma region testing_ultra_reactive_tabu


printf("====================ULTRA_REACTIVE_TABU====================\n");
printf("Start tabu size: %ld\n", tabu_size);
printf("Restarts:%d\n", iterations_count);
printf("UR_Constant = %ld\n", UR_multiplier);

if (isTimeStopCondition)
	printf("Time per iteration, s. : %ld\n", test_time_seconds);
else
	printf("Steps per iteration : %ld\n", max_steps);

strcpy_s(path, "d:\\data_maxcut\\results\\ultra_reactive_tabu\\log-");
strcat_s(path, datetimeprefix);
strcat(path, ".txt");
if (file != NULL)
	fclose(file);
fopen_s(&file, path, "a");

//PRINT TO FILE

fprintf_s(file, "====================ULTRA_REACTIVE_TABU====================\n");
fprintf_s(file, "Start tabu size: %ld\n", tabu_size);
fprintf_s(file, "Restarts:%d\n", iterations_count);
fprintf_s(file, "UR_Constant = %ld\n", UR_multiplier);
if (isTimeStopCondition)
	fprintf_s(file, "Time per iteration, s. : %ld\n", test_time_seconds);
else
	fprintf_s(file, "Steps till restart : %ld\n", max_steps);

fprintf_s(file, "\nIndex\tF\tSteps\tTime\n");
fclose(file);

for (int i = 0; i < iterations_count; i++) {
	printf("Restart %d : \t", i + 1);

#pragma region Initialization
	initRandom();

	for (int i = 0; i < graph.n_vertices; i++) {
		last_used[i] = -3 * graph.n_vertices;
	}

	for (int i = 0; i < graph.n_vertices; i++) {
		history.hashes[i] = getStandardRandom() * 100;
	}

	index_count = (long *)calloc(graph.n_vertices, sizeof(long));
	if (index_count == NULL) {
		printf("It is not enough free memory for array reactive_tabu_size_statistics\n");
		goto NOT_ENOUGH_FREE_MEMORY;
	}

	//Generating random starting vector x;
	for (int i = 0; i < graph.n_vertices; i++) {
		if (getStandardRandom() > 0.5) {
			x[i] = 1;
			best_x[i] = 1;
		}
		else {
			x[i] = 0;
			best_x[i] = 0;
		}
	}

	//calculate f(x) and save as best known.
	f = best_f = F(x, graph.edges);
	//TODO: best_f will be written to file so spprove with common tabu testing
	x_key = getKey(x);
	step = 0;
	best_f_step = 0;
	history.stored_f_count = 0;
	tabu_change_t = 0;
	start = clock();
	finish = start + CLOCKS_PER_SEC * (test_time_seconds);
#pragma endregion

#pragma region Main_Loop

	while (runCondition()) {
		ultra_reactive_tabu_noac(graph,&f, x, &x_key, &step, max_steps, &tabu_size, neighbourhood, neighbourhood_values, last_used, &best_f, best_x, &best_f_step, &tabu_change_t, ultra_reactive_tabu_sizes,index_count);
		
	}
	elapsed_time[i] = (clock() - start) / CLOCKS_PER_SEC;
	iteration_steps[i] = step;
#pragma endregion

#pragma region Results
	printf("Best F:%ld\tIterations: %ld\tTime: %f\n", best_f, iteration_steps[i], elapsed_time[i]);
	iteration_bestf[i] = best_f;

	if (best_f > best_known_f) {
		best_known_f = best_f;
		if (file != NULL)
			fclose(file);

		strcpy_s(path, "d:\\data_maxcut\\G37_bks.txt");
		fopen_s(&file, path, "w");
		fprintf_s(file, "%ld", best_f);
		fclose(file);
	}

	if (file != NULL)
		fclose(file);

	printResultToFile(i, "ultra_reactive_tabu");
	addDateTimePrefixToFileName(path, "d:\\data_maxcut\\results\\ultra_reactive_tabu\\log-", ".txt", datetimeprefix);
	appendIterationResult(i + 1, best_f, iteration_steps[i], elapsed_time[i], path);
	graph_move_index_log_ultra_reactive();
#pragma endregion

}

#pragma endregion


/*
	#pragma region COMPARISON_LOG
	strcpy_s(path, "d:\\data_maxcut\\results\\comparison\\tabu-");
	strcat_s(path, datetimeprefix);
	strcat(path, ".txt");
	fopen_s(&file, path, "w");

	fprintf_s(file, "ITERATION\tBEST F\tSTEPS\tTIME\n");
	for (int i = 0; i < iterations_count; i++){
		fprintf_s(file, "%d,\t\t%ld,\t%ld\t%f;\n", i, iteration_bestf[i],iteration_steps[i],elapsed_time[i]);
	}

	fclose(file);


	strcpy_s(path, "d:\\data_maxcut\\results\\comparison\\reactive-");
	strcat_s(path, datetimeprefix);
	strcat(path, ".txt");
	fopen_s(&file, path, "w");

	fprintf_s(file, "ITERATION\tBEST F\tSTEPS\tTIME\n");
	for (int i = 0; i < iterations_count; i++){
		fprintf_s(file, "%d,\t\t%ld,\t%ld\t%f;\n", i, iteration_bestf_reactive[i],iteration_steps_reactive[i],elapsed_time_reactive[i]);
	}

	fclose(file);


#pragma endregion
*/

#pragma region NOT_ENOUGH_FREE_MEMORY
	
NOT_ENOUGH_FREE_MEMORY :
		system("pause");// exit the program to prevent executing				   
		return;
#pragma endregion		
}

#pragma region ALGORITHMS

void tabu(graph_s graph,long *f, long *x, long n_vertices, int** edges, long *step, long max_steps, long *tabu_size, long *stored_f_count, bool *neighbourhood, long *neighbourhood_values, long *last_used, long *best_f, long *best_x, long *best_f_step){
	long best_neighbour_index = -1, best_neighbour_delta_value = LONG_MIN;
	long best_neighbour_ac_index = -1;//aspiration criteria: if move can increase best solution but inside tabu : do it.
	long best_neighbour_ac_delta_value = LONG_MIN;

#pragma region References
	void getNeighbourhood(bool* neighbours, long* last_used, long current_step, long tabu_size);
	long F(long *x, long f, long index, int** matrix);
	long F(long *x, int** matrix);
	long Fdelta(long *x, long index, int** matrix);
	long Fdelta(long *x, long index, struct graph_s graph);
	double getKey(long *x);
	void copyX(long* source, long* destination, int length);
#pragma endregion

	*step = *step + 1;

	getNeighbourhood(neighbourhood, last_used, *step, *tabu_size);

	for (int i = 0; i < graph.n_vertices; i++) {
		neighbourhood_values[i] = Fdelta(x, i, graph);

		if (neighbourhood[i]) {
			if (neighbourhood_values[i] > best_neighbour_delta_value)
			{
				best_neighbour_index = i;
				best_neighbour_delta_value = neighbourhood_values[i];
			}
		}
		else { //check aspiration criteria
			if (neighbourhood_values[i] > best_neighbour_ac_delta_value)
			{
				best_neighbour_ac_index = i;
				best_neighbour_ac_delta_value = neighbourhood_values[i];
			}
		}
	}


	if (*f + best_neighbour_ac_delta_value > *best_f)
	{
		if (best_neighbour_delta_value > best_neighbour_ac_delta_value) 
		{
			*f = *f + best_neighbour_delta_value;
			x[best_neighbour_index] = !x[best_neighbour_index];
			last_used[best_neighbour_index] = *step;
			index_count[best_neighbour_index]++;
		}
		else 
		{
			*f = *f + best_neighbour_ac_delta_value;
			x[best_neighbour_ac_index] = !x[best_neighbour_ac_index];
			last_used[best_neighbour_ac_index] = *step;
			index_count[best_neighbour_ac_index]++;
		}
		
	}
	else
	{
		*f = *f + best_neighbour_delta_value;
		x[best_neighbour_index] = !x[best_neighbour_index];
		last_used[best_neighbour_index] = *step;
		index_count[best_neighbour_index]++;
	}

	
	if (*f + best_neighbour_delta_value > *best_f)
	{
		/**best_f_step = *step;
		*best_f = *f + best_neighbour_delta_value;
		strcpy_s(path, "d:\\data_maxcut\\results\\tabu\\best_x.txt");
		fopen_s(&file, path, "w");

		fprintf_s(file, "=====F=====\n%ld\n=====X=====\n", *best_f);
		for (int i = 0; i < n_vertices; i++)
			fprintf_s(file, "%d", x[i]);
		fclose(file);
		*/
		*best_f_step = *step;
		*best_f = *f + best_neighbour_delta_value;

		for (int ii = 0; ii < graph.n_vertices; ii++)
			best_x[ii] = x[ii];
	}

	best_neighbour_index = -1;
	best_neighbour_ac_index = -1;
	best_neighbour_delta_value = LONG_MIN;
	best_neighbour_ac_delta_value = LONG_MIN;
}

void reactive_tabu(graph_s graph,  long *f, long *x, double *x_key, long *step, long max_steps, long *tabu_size,
	bool *neighbourhood, long *neighbourhood_values, long *last_used, long *best_f, long *best_x, long *best_f_step, long* tabu_change_t, long* reactive_tabu_size_statistics)
{
	long best_neighbour_index = -1, best_neighbour_delta_f = LONG_MIN, node_index = 0, deltaVisitTime = 0;
	long best_neighbour_ac_index = -1;//aspiration criteria: if move can increase best solution but inside tabu : do it.
	long best_neighbour_ac_delta_value = LONG_MIN;

#pragma region References
	void getNeighbourhood(bool* neighbours, long* last_used, long current_step, long tabu_size);
	long F(long *x, long f, long index, int** matrix);
	long F(long *x, int** matrix);
	long Fdelta(long *x, long index, struct graph_s graph);

	double getKey(long *x);

	long save_solution(double key, long f, struct history_s *history, long step);
	long check_solution(double key, long f, struct history_s *history);
	void Increase(long* tabu_size, long n_vertices);
	void Decrease(long* tabu_size, long n_vertices);
	void copyX(long* source, long* destination, int length);
#pragma endregion

	*x_key = getKey(x);

#pragma region Reaction
	if (node_index = check_solution(*x_key, *f, &history) > 0)
	{
		previous_solution_found_count++;
		deltaVisitTime = *step - history.occurence_time[node_index];

		history.occurence_time[node_index] = *step;
		if (deltaVisitTime < graph.n_vertices)
		{
			*tabu_change_t = *step;
			Increase(tabu_size, graph.n_vertices);
		}
	}
	else
	{
		save_solution(*x_key, *f, &history, *step);

		if (*step - *tabu_change_t > 125)
		{
			*tabu_change_t = *step;
			Decrease(tabu_size, graph.n_vertices);
		}
	}
	reactive_tabu_size_statistics[*tabu_size]++;
#pragma endregion

	getNeighbourhood(neighbourhood, last_used, *step, *tabu_size);

	*step = *step + 1;

	for (int i = 0; i < graph.n_vertices; i++) {
		neighbourhood_values[i] = Fdelta(x, i, graph);

		if (neighbourhood[i]) {
			if (neighbourhood_values[i] > best_neighbour_delta_f)
			{
				best_neighbour_index = i;
				best_neighbour_delta_f = neighbourhood_values[i];
			}
		}
		else { //check aspiration criteria
			if (neighbourhood_values[i] > best_neighbour_ac_delta_value)
			{
				best_neighbour_ac_index = i;
				best_neighbour_ac_delta_value = neighbourhood_values[i];
			}
		}
	}

	if (*f + best_neighbour_ac_delta_value > *best_f)
	{
		if (best_neighbour_delta_f > best_neighbour_ac_delta_value)
		{
			*f = *f + best_neighbour_delta_f;
			x[best_neighbour_index] = !x[best_neighbour_index];
			last_used[best_neighbour_index] = *step;
			index_count[best_neighbour_index]++;
		}
		else
		{
			*f = *f + best_neighbour_ac_delta_value;
			x[best_neighbour_ac_index] = !x[best_neighbour_ac_index];
			last_used[best_neighbour_ac_index] = *step;
			index_count[best_neighbour_ac_index]++;
		}

	}
	else
	{
		*f = *f + best_neighbour_delta_f;
		x[best_neighbour_index] = !x[best_neighbour_index];
		last_used[best_neighbour_index] = *step;
		index_count[best_neighbour_index]++;
	}

	/*if (*step % 50 == 0){
		printf("Step: %ld \t F: %ld Best F: %ld  Key: %f  Occurence: %ld \t Tabu: %ld\n", *step, best_neighbour_value, *best_f, *x_key,occurence_time[node_index], *tabu_size);
	}*/

	//copyX(x, best_x, n_vertices);

	if (*f + best_neighbour_delta_f > *best_f)
	{
		*best_f_step = *step;
		*best_f = *f + best_neighbour_delta_f;

		for (int ii = 0; ii < graph.n_vertices; ii++)
			best_x[ii] = x[ii];

	}

	best_neighbour_index = -1;
	best_neighbour_ac_index = -1;
	best_neighbour_delta_f = LONG_MIN;
	best_neighbour_ac_delta_value = LONG_MIN;

}

void ultra_reactive_tabu(graph_s graph,  long *f, long *x, double *x_key, long *step, long max_steps, long *tabu_size,
	bool *neighbourhood, long *neighbourhood_values, long *last_used, long *best_f, long *best_x, long *best_f_step, long* tabu_change_t, long* ultra_reactive_tabu_sizes, long* ultra_reactive_index_count)
{
	long best_neighbour_index = -1, best_neighbour_delta_f = LONG_MIN;
	long best_neighbour_ac_index = -1;//aspiration criteria: if move can increase best solution but inside tabu : do it.
	long best_neighbour_ac_delta_f = LONG_MIN;

#pragma region References

	void getNeighbourhoodUR(bool* neighbours, long* last_used, long current_step, long * tabu_sizes);
	void update_ur_tabu_size(long index, long step, long* tabu_sizes, long* index_count);
	long F(long *x, long f, long index, int** matrix);
	long F(long *x, int** matrix);
	long Fdelta(long *x, long index, struct graph_s graph);

	double getKey(long *x);

	void copyX(long* source, long* destination, int length);
#pragma endregion

	*x_key = getKey(x);


	getNeighbourhoodUR(neighbourhood, last_used, *step, ultra_reactive_tabu_sizes);

	*step = *step + 1;

	for (int i = 0; i < graph.n_vertices; i++) {
		neighbourhood_values[i] = Fdelta(x, i, graph);

		if (neighbourhood[i]) {
			if (neighbourhood_values[i] > best_neighbour_delta_f)
			{
				best_neighbour_index = i;
				best_neighbour_delta_f = neighbourhood_values[i];
			}
		}
		else { //check aspiration criteria
			if (neighbourhood_values[i] > best_neighbour_ac_delta_f)
			{
				best_neighbour_ac_index = i;
				best_neighbour_ac_delta_f = neighbourhood_values[i];
			}
		}
	}

	if (*f + best_neighbour_ac_delta_f > *best_f)
	{
		if (best_neighbour_delta_f > best_neighbour_ac_delta_f)
		{
			*f = *f + best_neighbour_delta_f;
			x[best_neighbour_index] = !x[best_neighbour_index];
			last_used[best_neighbour_index] = *step;
			ultra_reactive_index_count[best_neighbour_index]++;
			update_ur_tabu_size(best_neighbour_index, *step, ultra_reactive_tabu_sizes, ultra_reactive_index_count);

			if (*f + best_neighbour_delta_f > *best_f)
			{
				*best_f_step = *step;
				*best_f = *f + best_neighbour_delta_f;

				for (int ii = 0; ii < graph.n_vertices; ii++)
					best_x[ii] = x[ii];
			}

		}
		else
		{
			*f = *f + best_neighbour_ac_delta_f;
			x[best_neighbour_ac_index] = !x[best_neighbour_ac_index];
			last_used[best_neighbour_ac_index] = *step;
			ultra_reactive_index_count[best_neighbour_ac_index]++;
			update_ur_tabu_size(best_neighbour_ac_index, *step, ultra_reactive_tabu_sizes, ultra_reactive_index_count);

			if (*f + best_neighbour_ac_delta_f > *best_f)
			{
				*best_f_step = *step;
				*best_f = *f + best_neighbour_delta_f;

				for (int ii = 0; ii < graph.n_vertices; ii++)
					best_x[ii] = x[ii];
			}

		}

	}
	else
	{
		*f = *f + best_neighbour_delta_f;
		x[best_neighbour_index] = !x[best_neighbour_index];
		last_used[best_neighbour_index] = *step;
		ultra_reactive_index_count[best_neighbour_index]++;
		update_ur_tabu_size(best_neighbour_index, *step, ultra_reactive_tabu_sizes, ultra_reactive_index_count);

		if (*f + best_neighbour_delta_f > *best_f)
		{
			*best_f_step = *step;
			*best_f = *f + best_neighbour_delta_f;

			for (int ii = 0; ii < graph.n_vertices; ii++)
				best_x[ii] = x[ii];
		}

	}

	
}


void ultra_reactive_tabu_noac(graph_s graph, long *f, long *x, double *x_key, long *step, long max_steps, long *tabu_size,
	bool *neighbourhood, long *neighbourhood_values, long *last_used, long *best_f, long *best_x, long *best_f_step, long* tabu_change_t, long* ultra_reactive_tabu_sizes, long* ultra_reactive_index_count)
{
	long best_neighbour_index = -1, best_neighbour_delta_f = LONG_MIN;
	long best_neighbour_ac_index = -1;//aspiration criteria: if move can increase best solution but inside tabu : do it.
	long best_neighbour_ac_delta_f = LONG_MIN;

#pragma region References

	void getNeighbourhoodUR(bool* neighbours, long* last_used, long current_step, long * tabu_sizes);
	void update_ur_tabu_size(long index, long step, long* tabu_sizes, long* index_count);
	long F(long *x, long f, long index, int** matrix);
	long F(long *x, int** matrix);
	long Fdelta(long *x, long index, struct graph_s graph);

	double getKey(long *x);

	void copyX(long* source, long* destination, int length);
#pragma endregion

	*x_key = getKey(x);


	getNeighbourhoodUR(neighbourhood, last_used, *step, ultra_reactive_tabu_sizes);

	*step = *step + 1;

	for (int i = 0; i < graph.n_vertices; i++) {
		neighbourhood_values[i] = Fdelta(x, i, graph);

		if (neighbourhood[i]) {
			if (neighbourhood_values[i] > best_neighbour_delta_f)
			{
				best_neighbour_index = i;
				best_neighbour_delta_f = neighbourhood_values[i];
			}
		}	
	}


	*f = *f + best_neighbour_delta_f;
	x[best_neighbour_index] = !x[best_neighbour_index];
	last_used[best_neighbour_index] = *step;
	ultra_reactive_index_count[best_neighbour_index]++;
	update_ur_tabu_size(best_neighbour_index, *step, ultra_reactive_tabu_sizes, ultra_reactive_index_count);


	if (*f + best_neighbour_delta_f > *best_f)
	{
		*best_f_step = *step;
		*best_f = *f + best_neighbour_delta_f;

		for (int ii = 0; ii < graph.n_vertices; ii++)
			best_x[ii] = x[ii];
	}
}

#pragma endregion

#pragma region Helper_Functions

bool runCondition() {
	return (current = clock() < finish);
}

bool runCondition_steps() {
	return (step < max_steps);
}

void copyX(long* source, long* destination, int length) {
	for (int i = 0; i < length; i++)
		destination[i] = source[i];
}

void Increase(long* tabu_size, long n_vertices) {
	*tabu_size = (long)fmin(fmax(1.15*(*tabu_size) + 1, *tabu_size + 1), n_vertices - 2);
}

void Decrease(long* tabu_size, long n_vertices) {
	*tabu_size = (long)fmax(fmin(0.75*(*tabu_size), *tabu_size - 1), 1);
}

void initRandom() {
#define ia 843314861
#define ib 453816693
#define mic 1693666955
#define mmu 1073741824
#define psu 4.65661287307739258e-10
	rand_seed = 90853;
	srand(time(NULL));

}

double getRandom() {
	rand_seed = rand_seed*ia;
	if (rand_seed > mic)
		rand_seed = (rand_seed - mmu) - mmu;
	rand_seed = rand_seed + ib;
	if (rand_seed < 0)
		rand_seed = (rand_seed + mmu) + mmu;
	return(rand_seed*psu);
}

double getStandardRandom() {
	return (double)rand() / (double)RAND_MAX;
}

void clearArray(long* arr) {
	for (int i = 0; i < graph.n_vertices; i++)
		arr[i] = 0;
}

double getKey(long *x) {
	double key = 0;
	for (int i = 0; i < graph.n_vertices; i++) {
		if (x[i] == 1)
			key += history.hashes[i];
	}
	return key;
}

void getNeighbourhood(bool* neighbours, long* last_used, long current_step, long tabu_size) {
	for (int i = 0; i < graph.n_vertices; i++) {
		if (current_step - last_used[i] < tabu_size)
			neighbours[i] = false;
		else
			neighbours[i] = true;
	}
}

//ULTRA REACTIVE
void getNeighbourhoodUR(bool* neighbours, long* last_used, long current_step, long * tabu_sizes) {
	for (int i = 0; i < graph.n_vertices; i++) {
		if (current_step - last_used[i] < tabu_sizes[i])
			neighbours[i] = false;
		else
			neighbours[i] = true;
	}
}

void update_ur_tabu_size(long index, long step, long* tabu_sizes, long* index_count) {
	long UR_size(long count, long total_count);

	tabu_sizes[index] = UR_size(index_count[index], step);
}
long UR_size(long count, long total_count) {

	double p = (double)count / (total_count+0.01f);

	return (long)(min(UR_multiplier*(p / (1 - p)), graph.n_vertices - 2));
}

long F(long *x, int** matrix) {
	long sum = 0;
	for (int i = 0; i < graph.n_vertices; i++) {
		for (int j = i + 1; j < graph.n_vertices; j++) {
			sum += matrix[i][j] * (x[i] ^ x[j]);
		}
	}
	return sum;
}

long F(long *x, long f, long index, int** matrix) {
	long minus = 0, plus = 0;
	long xi = x[index];
	long nxi = !x[index];
	for (int j = 0; j < graph.n_vertices; j++) {
		minus += matrix[index][j] * (xi ^ x[j]);
		plus += matrix[index][j] * (nxi ^ x[j]);
	}
	return f - minus + plus;
}

long Fdelta(long *x, long index, int** matrix) {
	long minus = 0, plus = 0;
	long xi = x[index];
	long nxi = !x[index];
	for (int j = 0; j < graph.n_vertices; j++) {
		minus += matrix[index][j] * (xi ^ x[j]);
		plus += matrix[index][j] * (nxi ^ x[j]);
	}
	return plus - minus;
}

//what delta f we will get if swap i-th bit.
long Fdelta(long *x, long i, struct graph_s graph) {
	long k = 0, delta = 0;
	if (x[i] == 1)
	{
		for (int j = graph.beg_list_edges_var[i]; j < graph.beg_list_edges_var[i + 1]; j++)
		{
			k = graph.list_nodes_var[j];

			if (x[k] == 0)
				delta -= graph.edge_weight[graph.list_edges_var[j]];
			else
				delta += graph.edge_weight[graph.list_edges_var[j]];
		}
	}
	else
	{
		for (int j = graph.beg_list_edges_var[i]; j < graph.beg_list_edges_var[i + 1]; j++)
		{
			k = graph.list_nodes_var[j];

			if (x[k] == 0)
				delta += graph.edge_weight[graph.list_edges_var[j]];
			else
				delta -= graph.edge_weight[graph.list_edges_var[j]];
		}
	}
	return delta;
}

long save_solution(double key, long f, struct history_s *h, long step)
{
	long node, item;
	if (h->stored_f_count == 0)
	{
		h->stored_f_count++;
		h->keys[h->stored_f_count] = key;
		h->f_values[h->stored_f_count] = f;
		h->left[h->stored_f_count] = 0;
		h->right[h->stored_f_count] = 0;
		h->occurence_time[h->stored_f_count] = 0;
	}
	else
	{
		node = 1;
		while (node)
		{
			if (h->keys[node] == key && h->f_values[node] == f)
			{
				h->occurence_time[node] = step;
				return node;
			}
			else
			{
				item = node;
				if (h->keys[node] > key)node = h->left[node];
				else node = h->right[node];
			}
		}
		h->stored_f_count++;
		h->keys[h->stored_f_count] = key;
		h->f_values[h->stored_f_count] = f;
		h->left[h->stored_f_count] = h->right[h->stored_f_count] = h->occurence_time[h->stored_f_count] = 0;

		if (h->keys[item] > key)
			h->left[item] = h->stored_f_count;
		else
			h->right[item] = h->stored_f_count;
	}
	return 0;
}

long check_solution(double key, long f, struct history_s *h)
{
	long node, item;
	if (h->stored_f_count > 0)
	{
		node = 1;
		while (node)
		{
			if (fabs(h->keys[node] - key) < eps && fabs(h->f_values[node] - f) < eps)//&& h.f_values[node] == f
			{
				return node;
			}
			else
			{
				item = node;
				if (h->keys[node] > key)node = h->left[node];
				else node = h->right[node];
			}
		}
	}
	return 0;
}

void tabu_size_log_reactive(long iteration_test_count, long steps) {
	long sum = 0;
	for (int i = 0; i < graph.n_vertices; i++) {
		sum += i*reactive_tabu_size_statistics[i];
	}

	double mean = (double)sum / (iteration_test_count*steps);

	strcpy_s(path, "d:\\data_maxcut\\results\\reactive_tabu\\size-");
	strcat_s(path, datetimeprefix);
	strcat(path, ".txt");
	fopen_s(&file, path, "w");

	fprintf_s(file, "steps: %ld\n", steps);
	fprintf_s(file, "mean: \t\t%f\n", mean);
	fprintf_s(file, "previous_sol_found_count :%ld\n", previous_solution_found_count);
	fprintf_s(file, "data:\n");
	fprintf_s(file, "SIZE\tCOUNT\n");
	for (int i = 0; i < graph.n_vertices; i++) {
		fprintf_s(file, "%d,\t%ld;\n", i, reactive_tabu_size_statistics[i]);
	}
	fclose(file);
}

void graph_move_index_log() {

	FILE *f;

	strcpy_s(path, "d:\\data_maxcut\\results\\tabu\\index_stats-");
	strcat_s(path, datetimeprefix);
	strcat(path, ".txt");
	fopen_s(&f, path, "w");

	fprintf_s(f, "Total steps: %ld\n", step);
	for (int i = 0; i < graph.n_vertices; i++) {
		fprintf(f, "%ld,\t%ld;\n", i, index_count[i]);
	}
	fclose(f);
}

void graph_move_index_log_reactive() {

	strcpy_s(path, "d:\\data_maxcut\\results\\reactive_tabu\\index_stats-");
	strcat_s(path, datetimeprefix);
	strcat(path, ".txt");
	fopen_s(&file, path, "w");

	fprintf_s(file, "Total steps: %ld\n", step);
	for (int i = 0; i < graph.n_vertices; i++) {
		fprintf(file, "%ld,\t%ld;\n", i, index_count[i]);
	}
	fclose(file);
}

void graph_move_index_log_ultra_reactive() {

	strcpy_s(path, "d:\\data_maxcut\\results\\ultra_reactive_tabu\\index_stats-");
	strcat_s(path, datetimeprefix);
	strcat(path, ".txt");
	fopen_s(&file, path, "w");

	fprintf_s(file, "Total steps: %ld\n", step);
	for (int i = 0; i < graph.n_vertices; i++) {
		fprintf(file, "%ld,\t%ld;\n", i, index_count[i]);
	}
	fclose(file);
}


char* getDateTimePrefix() {
	time_t now = time(0);
	struct tm * timeinfo;

	timeinfo = localtime(&now);
	char* buffer = new char[80];
	strftime(buffer, 80, "%Y%m%d-%H%M%S", timeinfo);
	return buffer;
}

#pragma region print functions

void appendIterationResult(int index, long F, long steps, float elapsed_time, char* filename) {
	FILE *f;
	char path[200];

	strcpy_s(path, filename);

	fopen_s(&file, path, "a");

	fprintf_s(file, "%d,\t%ld,\t%ld\t%f;\n", index, F, steps, elapsed_time);

	fclose(file);
}

void addDateTimePrefixToFileName(char* dest, const char* first, const char* second, char* dtprefix) {
	strcpy(dest, first);
	strcat(dest, dtprefix);
	strcat(dest, second);
}

void printResultToFile(int ind, const char* alg_name) {
	char intstring[4];
	char path[200];
	FILE *f;

	strcpy(path, "d:\\data_maxcut\\results\\");
	strcat(path, alg_name);
	strcat(path, "\\best_x-");
	strcat(path, datetimeprefix);
	strcat(path, "-");
	sprintf(intstring, "%d", ind);
	strcat(path, intstring);
	strcat(path, ".txt");
	fopen_s(&f, path, "w");

	fprintf_s(f, "=====F=====\n%ld\n=====X=====\n", best_f);
	for (int i = 0; i < graph.n_vertices; i++)
		fprintf_s(f, "%d", best_x[i]);
	fclose(f);

}

bool allocateMemory(long** x,long** best_x, long** last_used,bool** neighbourhood,long** neighbourhood_values,history_s* history, graph_s* graph,
	long** iteration_bestf,long** iteration_bestf_reactive,long iterations_count, long** iteration_steps, long** iteration_steps_reactive,long** reactive_tabu_size_statistics,
	long** index_count,long** ultra_reactive_tabu_sizes, double** elapsed_time, double** elapsed_time_reactive, double** elapsed_time_ultra_reactive) {
	*x = (long *)calloc(graph->n_vertices, sizeof(long));
	if (x == NULL) {
		printf("It is not enough free memory for x\n");
		return false;
	}

	*best_x = (long *)calloc(graph->n_vertices, sizeof(long));
	if (x == NULL) {
		printf("It is not enough free memory for best_x\n");
		return false;
	}

	*last_used = (long *)calloc(graph->n_vertices, sizeof(long));
	if (last_used == NULL) {
		printf("It is not enough free memory for array last_used\n");
		return false;
	}

	*neighbourhood = (bool *)calloc(graph->n_vertices, sizeof(bool));
	if (neighbourhood == NULL) {
		printf("It is not enough free memory for array neighbourhood\n");
		return false;
	}

	*neighbourhood_values = (long *)calloc(graph->n_vertices, sizeof(long));
	if (neighbourhood_values == NULL) {
		printf("It is not enough free memory for array neighbourhood_values\n");
		return false;
	}

#pragma region history
	history->keys = (double *)calloc(max_steps, sizeof(double));
	if (history->keys == NULL) {
		printf("It is not enough free memory for array keyar\n");
		return false;
	}

	history->f_values = (long *)calloc(max_steps, sizeof(long));
	if (history->f_values == NULL) {
		printf("It is not enough free memory for array f_values\n");
		return false;
	}

	history->left = (long *)calloc(max_steps, sizeof(long));
	if (history->left == NULL) {
		printf("It is not enough free memory for array left\n");
		return false;
	}

	history->right = (long *)calloc(max_steps, sizeof(long));
	if (history->right == NULL) {
		printf("It is not enough free memory for array right\n");
		return false;
	}

	history->occurence_time = (long *)calloc(max_steps, sizeof(long));
	if (history->occurence_time == NULL) {
		printf("It is not enough free memory for array occurence_count\n");
		return false;
	}

	history->hashes = (double *)calloc(graph->n_vertices, sizeof(double));
	if (history->hashes == NULL) {
		printf("It is not enough free memory for array hashes\n");
		return false;
	}

#pragma endregion

#pragma region graph
	graph->edges = (int **)calloc(graph->n_vertices, sizeof(int *));
	if (graph->edges == NULL) {
		printf("It is not enough free memory for array edges\n");
		return false;
	}

	for (int i = 0; i < graph->n_vertices; i++) {
		graph->edges[i] = (int *)calloc(graph->n_edges, sizeof(int));
		if (graph->edges[i] == NULL) {
			printf("It is not enough free memory for array edges[i]\n");
			return false;
		}
	}

	graph->n_con_edges = (int *)calloc(graph->n_vertices, sizeof(int));
	if (graph->n_con_edges == NULL) {
		printf("It is not enough free memory for array n_con_edges\n");
		return false;
	}

	graph->beg_list_edges_var = (int *)calloc(graph->n_vertices + 1, sizeof(int));
	if (graph->beg_list_edges_var == NULL) {
		printf("It is not enough free memory for array beg_list_edges_var\n");
		return false;
	}

	graph->list_edges_var = (int *)calloc(2 * graph->n_edges, sizeof(int));
	if (graph->list_edges_var == NULL) {
		printf("It is not enough free memory for array list_edges_var\n");
		return false;
	}

	graph->list_nodes_var = (int *)calloc(2 * graph->n_edges, sizeof(int));
	if (graph->list_nodes_var == NULL) {
		printf("It is not enough free memory for array list_nodes_var\n");
		return false;
	}

	graph->edge_weight = (int *)calloc(graph->n_edges, sizeof(int));
	if (graph->edge_weight == NULL) {
		printf("It is not enough free memory for array edge_weight\n");
		return false;
	}

#pragma endregion

	*iteration_bestf = (long *)calloc(iterations_count, sizeof(long));
	if (iteration_bestf == NULL) {
		printf("It is not enough free memory for array iteration_bestf\n");
		return false;
	}

	*iteration_bestf_reactive = (long *)calloc(iterations_count, sizeof(long));
	if (iteration_bestf_reactive == NULL) {
		printf("It is not enough free memory for array iteration_bestf_reactive\n");
		return false;
	}

	*iteration_steps = (long *)calloc(iterations_count, sizeof(long));
	if (iteration_steps == NULL) {
		printf("It is not enough free memory for array iteration_steps\n");
		return false;
	}

	*iteration_steps_reactive = (long *)calloc(iterations_count, sizeof(long));
	if (iteration_steps_reactive == NULL) {
		printf("It is not enough free memory for array iteration_steps\n");
		return false;
	}

	*reactive_tabu_size_statistics = (long *)calloc(graph->n_vertices, sizeof(long));
	if (reactive_tabu_size_statistics == NULL) {
		printf("It is not enough free memory for array reactive_tabu_size_statistics\n");
		return false;
	}

	*index_count = (long *)calloc(graph->n_vertices, sizeof(long));
	if (index_count == NULL) {
		printf("It is not enough free memory for array reactive_tabu_size_statistics\n");
		return false;
	}


	*ultra_reactive_tabu_sizes = (long *)calloc(graph->n_vertices, sizeof(long));
	if (ultra_reactive_tabu_sizes == NULL) {
		printf("It is not enough free memory for array ultra_reactive_tabu_size_statistics\n");
		return false;
	}

	*elapsed_time = (double *)calloc(iterations_count, sizeof(long));
	if (elapsed_time == NULL) {
		printf("It is not enough free memory for array elapsed_time\n");
		return false;
	}

	*elapsed_time_reactive = (double *)calloc(iterations_count, sizeof(long));
	if (elapsed_time_reactive == NULL) {
		printf("It is not enough free memory for array elapsed_time\n");
		return false;
	}

	*elapsed_time_ultra_reactive = (double *)calloc(iterations_count, sizeof(long));
	if (elapsed_time_ultra_reactive == NULL) {
		printf("It is not enough free memory for array elapsed_time\n");
		return false;
	}
}
#pragma endregion

#pragma endregion
