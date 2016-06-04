#include <malloc.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <time.h>
#include <windows.h>
#include <string>
#include<omp.h>
using namespace std;

#pragma region TYPEDEFS

struct graph_s {
	int g_number;
	long n_edges;
	long n_vertices;
	int **edges;
	int * n_con_edges;
	int * beg_list_edges_var;
	int * list_edges_var;
	int * list_nodes_var;
	int * edge_weight;
};

struct history_s {
	long * f_values;
	long * left;
	long * right;
	long * occurence_time;
	double * keys;
	double * hashes;
	long stored_f_count = 0;

	long previous_solution_found_count = 0;
	long max_steps = 15000000;//size of history
	long* best_x;
	long best_f = MINLONG;
	long best_f_step = 0;
	long tabu_change_t = 0;
	long* index_count;		//a[i] = number of i index changes.
};

struct test_params_s {
	char* dateTimePrefix;
	char* alg_name;
	clock_t start,current, finish;
	bool is_time_stop_condition = true;
	long iterations_count = 5;
	long test_time_seconds = 15;
	long max_steps = 15000000;
	char* f_file_path;
};

struct test_history_s {
	char* start_time;
	char* test_folder_path ;
	long* iteration_steps;
	long* iteration_bestf;
	double* elapsed_time;
	long best_known_f;
};

struct config_s {
	long* x;
	double x_key;
	long f;
	long tabu_size = 50;
	long step;
	long* last_used;
	bool* neighbourhood;
	long* neighbourhood_values;
};

struct reactive_s {
	long* tabu_sizes;
	
	long* size_statistics;	//a[i] = number of steps with tabu size == i

	long steps_to_decrease = 100;

};

#pragma endregion

#pragma region CONST
const double eps = 0.00001;
const int UR_multiplier = 150000;

const char* GRAPHS = "D:\\data_maxcut\\graphs\\";

const char* TABU_RESULTS =		"D:\\data_maxcut\\results\\tabu\\";
const char* REACTIVE_RESULTS =	"D:\\data_maxcut\\results\\reactive_tabu\\";
const char* ULTRA_RESULTS =		"D:\\data_maxcut\\results\\ultra_reactive_tabu\\";

const char* BEST_KNOWN =		"D:\\data_maxcut\\records\\";

const char* TESTS = "D:\\data_maxcut\\tests\\";
#pragma endregion

#pragma region GLOBAL_VARS
long rand_seed;
#pragma endregion

#pragma region Helper_Functions

bool runCondition(test_params_s* params) {
	return ( params->current < params->finish);
}

bool runCondition_steps(long step, long max_steps) {
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

void clearArray(graph_s graph, long* arr) {
	for (int i = 0; i < graph.n_vertices; i++)
		arr[i] = 0;
}

double getKey(history_s* history, graph_s* graph, config_s* config) {
	double key = 0;
	for (int i = 0; i < graph->n_vertices; i++) {
		if (config->x[i] == 1)
			key += history->hashes[i];
	}
	return key;
}

void getNeighbourhood(graph_s* graph, bool* neighbours, long* last_used, long current_step, long tabu_size) {
	for (int i = 0; i < graph->n_vertices; i++) {
		if (current_step - last_used[i] < tabu_size)
			neighbours[i] = false;
		else
			neighbours[i] = true;
	}
}

//ULTRA REACTIVE
void getNeighbourhoodUR(graph_s* graph,bool* neighbours, long* last_used, long current_step, long * tabu_sizes) {
	for (int i = 0; i < graph->n_vertices; i++) {
		if (current_step - last_used[i] < tabu_sizes[i])
			neighbours[i] = false;
		else
			neighbours[i] = true;
	}
}

long UR_size(graph_s* graph,long count, long total_count) {

	double p = (double)count / (total_count + 0.01f);

	return (long)(min(UR_multiplier*(p / (1 - p)), graph->n_vertices - 2));
}

void update_ur_tabu_size(graph_s* graph, long index, long step, long* tabu_sizes, long* index_count) {

	tabu_sizes[index] = UR_size(graph,index_count[index], step);
}

long F(graph_s* graph, long *x, int** matrix) {
	long sum = 0;
	for (int i = 0; i < graph->n_vertices; i++) {
		for (int j = i + 1; j < graph->n_vertices; j++) {
			sum += matrix[i][j] * (x[i] ^ x[j]);
		}
	}
	return sum;
}

long F(graph_s* graph, long *x, long f, long index, int** matrix) {
	long minus = 0, plus = 0;
	long xi = x[index];
	long nxi = !x[index];
	for (int j = 0; j < graph->n_vertices; j++) {
		minus += matrix[index][j] * (xi ^ x[j]);
		plus += matrix[index][j] * (nxi ^ x[j]);
	}
	return f - minus + plus;
}

long Fdelta(graph_s* graph,long *x, long index, int** matrix) {
	long minus = 0, plus = 0;
	long xi = x[index];
	long nxi = !x[index];
	for (int j = 0; j < graph->n_vertices; j++) {
		minus += matrix[index][j] * (xi ^ x[j]);
		plus += matrix[index][j] * (nxi ^ x[j]);
	}
	return plus - minus;
}

//what delta f we will get if swap i-th bit.
long Fdelta(long *x, long i, struct graph_s* graph) {
	long k = 0, delta = 0;
	if (x[i] == 1)
	{
		for (int j = graph->beg_list_edges_var[i]; j < graph->beg_list_edges_var[i + 1]; j++)
		{
			k = graph->list_nodes_var[j];

			if (x[k] == 0)
				delta -= graph->edge_weight[graph->list_edges_var[j]];
			else
				delta += graph->edge_weight[graph->list_edges_var[j]];
		}
	}
	else
	{
		for (int j = graph->beg_list_edges_var[i]; j < graph->beg_list_edges_var[i + 1]; j++)
		{
			k = graph->list_nodes_var[j];

			if (x[k] == 0)
				delta += graph->edge_weight[graph->list_edges_var[j]];
			else
				delta -= graph->edge_weight[graph->list_edges_var[j]];
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

void tabu_size_log_reactive( int iteration_index, long iteration_test_count, long steps, int n_vertices, long* size_statistics, char* test_folder_path,  long prev_solutions_found) {
	FILE *file;
	char path[100];
	char intstring[4];
	long sum = 0;

	for (int i = 0; i < n_vertices; i++) {
		sum += i*size_statistics[i];
	}

	double mean = (double)sum / (steps);//double mean = (double)sum / (iteration_test_count*steps);

	//strcpy_s(path, "d:\\data_maxcut\\results\\reactive_tabu\\size-");
	strcpy_s(path, test_folder_path);
	strcat(path, "\\");
	strcat(path, "size");
	//strcat_s(path, params->dateTimePrefix);
	strcat(path, "-");
	sprintf(intstring, "%d", iteration_index);
	strcat(path, intstring);
	strcat(path, ".txt");
	fopen_s(&file, path, "w");

	fprintf_s(file, "steps: %ld\n", steps);
	fprintf_s(file, "mean: \t\t%f\n", mean);
	fprintf_s(file, "previous_sol_found_count :%ld\n", prev_solutions_found);
	fprintf_s(file, "data:\n");
	fprintf_s(file, "SIZE\tCOUNT\n");
	for (int i = 0; i < n_vertices; i++) {
		fprintf_s(file, "%d,\t%ld;\n", i, size_statistics[i]);
	}
	fclose(file);
}


void graph_move_index_log(long n_vertices, long iteration, long* index_count, char* test_folder_path) {
	FILE *file;
	char path[100];
	char intstr[4];

	strcpy_s(path, test_folder_path);
	strcat(path, "\\");
	strcat(path, "index_stats-");
	sprintf(intstr, "%d", iteration);
	strcat(path, intstr);
	strcat(path, ".txt");

	fopen_s(&file, path, "w");

	for (int i = 0; i < n_vertices; i++) {
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

	fopen_s(&f, path, "a");

	fprintf_s(f, "%d,\t%ld,\t%ld\t%f;\n", index, F, steps, elapsed_time);

	fclose(f);
}

void addDateTimePrefixToFileName(char* dest, const char* first, const char* second, char* dtprefix) {
	strcpy(dest, first);
	strcat(dest, dtprefix);
	strcat(dest, second);
}

void printResultToFile(long n_vertices,long best_f, long* best_x, int ind,char* test_folder_path) {
	char intstring[4];
	char path[200];
	FILE *f;

	strcpy(path, test_folder_path);
	strcat(path, "\\best_x-");
	sprintf(intstring, "%d", ind);
	strcat(path, intstring);
	strcat(path, ".txt");

	fopen_s(&f, path, "w");

	fprintf_s(f, "=====F=====\n%ld\n=====X=====\n", best_f);
	for (int i = 0; i < n_vertices; i++)
		fprintf_s(f, "%d", best_x[i]);
	fclose(f);

}

void log_f_value(long f, long step, test_params_s* params) {

	FILE* file;

	fopen_s(&file, params->f_file_path, "a+");
	fprintf_s(file, "%ld,\t%ld,\t%Lf\n", f, step, (params->current - params->start)/ (double)CLOCKS_PER_SEC);
	fclose(file);
}
#pragma endregion

#pragma region MEMORY
bool allocate_graph(graph_s* graph) {
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
	return true;
}

bool allocate_history(history_s* history, int n_vertices) {
	history->keys = (double *)calloc(history->max_steps, sizeof(double));
	if (history->keys == NULL) {
		printf("It is not enough free memory for array keyar\n");
		return false;
	}

	history->f_values = (long *)calloc(history->max_steps, sizeof(long));
	if (history->f_values == NULL) {
		printf("It is not enough free memory for array f_values\n");
		return false;
	}

	history->left = (long *)calloc(history->max_steps, sizeof(long));
	if (history->left == NULL) {
		printf("It is not enough free memory for array left\n");
		return false;
	}

	history->right = (long *)calloc(history->max_steps, sizeof(long));
	if (history->right == NULL) {
		printf("It is not enough free memory for array right\n");
		return false;
	}

	history->occurence_time = (long *)calloc(history->max_steps, sizeof(long));
	if (history->occurence_time == NULL) {
		printf("It is not enough free memory for array occurence_count\n");
		return false;
	}

	history->hashes = (double *)calloc(n_vertices, sizeof(double));
	if (history->hashes == NULL) {
		printf("It is not enough free memory for array hashes\n");
		return false;
	}

	history->index_count = (long *)calloc(n_vertices, sizeof(long));
	if (history->index_count == NULL) {
		printf("It is not enough free memory for array reactive_tabu_size_statistics\n");
		return false;
	}

	history->best_x = (long *)calloc(n_vertices, sizeof(long));
	if (history->best_x == NULL) {
		printf("It is not enough free memory for best_x\n");
		return false;
	}

	return true;
}

bool allocate_config(config_s* config, int n_vertices) {
	config->x = (long *)calloc(n_vertices, sizeof(long));
	if (config->x == NULL) {
		printf("It is not enough free memory for x\n");
		return false;
	}

	config->last_used = (long *)calloc(n_vertices, sizeof(long));
	if (config->last_used == NULL) {
		printf("It is not enough free memory for array last_used\n");
		return false;
	}

	config->neighbourhood = (bool *)calloc(n_vertices, sizeof(bool));
	if (config->neighbourhood == NULL) {
		printf("It is not enough free memory for array neighbourhood\n");
		return false;
	}

	config->neighbourhood_values = (long *)calloc(n_vertices, sizeof(long));
	if (config->neighbourhood_values == NULL) {
		printf("It is not enough free memory for array neighbourhood_values\n");
		return false;
	}
	return true;
}

bool allocate_test_history(test_history_s* test_history, long iterations_count) {
	test_history->iteration_bestf = (long *)calloc(iterations_count, sizeof(long));
	if (test_history->iteration_bestf == NULL) {
		printf("It is not enough free memory for array iteration_bestf\n");
		return false;
	}

	test_history->iteration_steps = (long *)calloc(iterations_count, sizeof(long));
	if (test_history->iteration_steps == NULL) {
		printf("It is not enough free memory for array iteration_steps\n");
		return false;
	}

	test_history->elapsed_time = (double *)calloc(iterations_count, sizeof(double));
	if (test_history->elapsed_time == NULL) {
		printf("It is not enough free memory for array elapsed_time\n");
		return false;
	}

	return true;
}

bool allocate_reactive(reactive_s* reactive, int n_vertices) {
	reactive->tabu_sizes = (long *)calloc(n_vertices, sizeof(long));
	if (reactive->tabu_sizes == NULL) {
		printf("It is not enough free memory for array reactive_tabu_sizes\n");
		return false;
	}

	reactive->size_statistics = (long*)calloc(n_vertices, sizeof(long));
	if (reactive->size_statistics == NULL) {
		printf("It is not enough free memory for array reactive_tabu_size_statistics\n");
		return false;
	}
	return true;
}

bool allocateMemory(graph_s* graph, history_s* history, config_s* config, test_params_s* params, test_history_s* test_history,reactive_s* reactive) {

	if (!allocate_graph(graph))
		return false;

	if (!allocate_history(history, graph->n_vertices))
		return false;

	if (!allocate_config(config, graph->n_vertices))
		return false;

	if (!allocate_test_history(test_history, params->iterations_count))
		return false;

	if (!allocate_reactive(reactive, graph->n_vertices))
		return false;
	
	return true;
}

void free_graph(graph_s* graph) {
	for (int i = 0; i < graph->n_vertices; i++) {
		free(graph->edges[i]);
	}
	free(graph->edges);
	free(graph->edge_weight);
	free(graph->beg_list_edges_var);
	free(graph->list_edges_var);
	free(graph->list_nodes_var);
	free(graph->n_con_edges);
}

void free_history(history_s* history) {
	free(history->best_x);
	free(history->f_values);
	free(history->hashes);
	free(history->index_count);
	free(history->keys);
	free(history->left);
	free(history->right);
	free(history->occurence_time);
}

void free_config(config_s* config) {
	free(config->x);
	free(config->last_used);
	free(config->neighbourhood);
	free(config->neighbourhood_values);
}

void free_test_history(test_history_s* test_history) {

	free(test_history->iteration_bestf);
	free(test_history->iteration_steps);
	free(test_history->elapsed_time);
}

void free_reactive(reactive_s* reactive) {
	free(reactive->tabu_sizes);
	free(reactive->size_statistics);
}

void cleanup(graph_s* graph, history_s *history, config_s *config, test_params_s *params, test_history_s *test_history, reactive_s *reactive) {
	free_reactive(reactive);
	free_test_history(test_history);
	free_config(config);
	free_history(history);
	free_graph(graph);
}

void iteration_cleanup(graph_s* graph, history_s *history, config_s *config, test_params_s *params, test_history_s *test_history) {
	free_history(history);
	allocate_history(history, graph->n_vertices);
	history->stored_f_count = 0;
	history->best_f = MINLONG;
	history->best_f_step = 0;
	history->previous_solution_found_count = 0;
	history->tabu_change_t = 0;
}
#pragma endregion

char* get_graph_path(int graph_number) {
	char* buffer = new char[100];
	char intstr[4];

	strcpy(buffer, GRAPHS);
	strcat(buffer, "G");
	sprintf(intstr, "%d", graph_number);
	strcat(buffer, intstr);
	strcat(buffer, ".mtx");

	return buffer;
}

void fetch_properties(graph_s *graph, long *best_known_f,int graph_number) {

	char path[100];
	char intstr[4];
	char* graph_path;
	FILE *f;

	graph->g_number = graph_number;

	//strcpy_s(path, "d:\\data_maxcut\\G37_bks.txt");
	strcpy_s(path, BEST_KNOWN);
	strcat(path, "G");
	sprintf(intstr, "%d", graph->g_number);
	strcat(path, intstr);
	strcat(path, ".txt");

	fopen_s(&f, path, "r");
	fscanf(f, "%ld\n", best_known_f);
	fclose(f);

	graph_path = get_graph_path(graph_number);
	fopen_s(&f, graph_path, "r");
	fscanf(f, "%ld%ld%ld", &(graph->n_vertices), &(graph->n_vertices), &(graph->n_edges));
	fclose(f);
}

bool is_graph_has_bool_weights(int i) {
	return (i >= 6 && i <= 13) || (i >= 18 && i <= 21) || (i >= 27 && i <= 34) || (i >= 39 && i <= 42);
}

char* graph_name(int g_number) {
	char* intstr = new char[4];
	char* buff = new char[40];
	sprintf(intstr, "%d", g_number);
	strcpy(buff, "G");
	strcat(buff, intstr);
	return buff;
}

void create_history_folder_paths(graph_s* graph, test_history_s *test_history, const char* alg_name) {
	char* path = new char[200];

	strcpy(path, TESTS);
	strcat(path, alg_name);
	strcat(path, "-");
	strcat(path, graph_name(graph->g_number));
	strcat(path, "-");
	strcat(path, test_history->start_time);
	
	test_history->test_folder_path = path;

	CreateDirectory(test_history->test_folder_path, NULL);
}

void fetch(graph_s* graph) {
	FILE *file;
	char* filename;
	long temp_vert1, temp_vert2, temp_weight;

	filename = get_graph_path(graph->g_number);
	fopen_s(&file, filename, "r");
	fscanf(file, "%ld%ld%ld", &(graph->n_vertices),&(graph->n_vertices), &(graph->n_edges));
	/*preparing graph by calculating count of adjacent vertices for every vertice.\*/
	if (is_graph_has_bool_weights(graph->g_number)) {
		for (int j = 0; j < graph->n_edges; j++)
		{
			fscanf(file, "%ld%ld%ld", &temp_vert1, &temp_vert2, &temp_weight);
			temp_vert1--;
			temp_vert2--;
			graph->edge_weight[j] = temp_weight;
			graph->n_con_edges[temp_vert1]++;
			graph->n_con_edges[temp_vert2]++;
			graph->edges[temp_vert1][temp_vert2] = temp_weight;
			graph->edges[temp_vert2][temp_vert1] = temp_weight;
		}
	}
	else {
		for (int j = 0; j < graph->n_edges; j++)
		{
			fscanf(file, "%ld%ld", &temp_vert1, &temp_vert2);
			temp_vert1--;
			temp_vert2--;
			graph->edge_weight[j] = 1;
			graph->n_con_edges[temp_vert1]++;
			graph->n_con_edges[temp_vert2]++;
			graph->edges[temp_vert1][temp_vert2] = 1;
			graph->edges[temp_vert2][temp_vert1] = 1;
		}
	}
	

	graph->beg_list_edges_var[0] = 0;

	for (int j = 0; j < graph->n_vertices; j++)
	{
		graph->beg_list_edges_var[j + 1] = graph->beg_list_edges_var[j] + graph->n_con_edges[j];
	}

	for (int j = 0; j < graph->n_vertices; j++)
		graph->n_con_edges[j] = 0;

	if (file != NULL)
		fclose(file);

	/*Fetching nodes*/
	fopen_s(&file, filename, "r");
	fscanf(file, "%ld%ld%ld", &(graph->n_vertices), &(graph->n_vertices), &(graph->n_edges));

	if (is_graph_has_bool_weights(graph->g_number)) {
		for (int j = 0; j < graph->n_edges; j++)
		{
			fscanf(file, "%ld%ld%ld", &temp_vert1, &temp_vert2, &temp_weight);
			temp_vert1--;
			temp_vert2--;
			graph->list_edges_var[graph->beg_list_edges_var[temp_vert1] + graph->n_con_edges[temp_vert1]] = j;
			graph->list_edges_var[graph->beg_list_edges_var[temp_vert2] + graph->n_con_edges[temp_vert2]] = j;
			graph->list_nodes_var[graph->beg_list_edges_var[temp_vert1] + graph->n_con_edges[temp_vert1]] = temp_vert2;
			graph->list_nodes_var[graph->beg_list_edges_var[temp_vert2] + graph->n_con_edges[temp_vert2]] = temp_vert1;
			graph->n_con_edges[temp_vert1]++;
			graph->n_con_edges[temp_vert2]++;
		}
	}
	else {
		for (int j = 0; j < graph->n_edges; j++)
		{
			fscanf(file, "%ld%ld", &temp_vert1, &temp_vert2);
			temp_vert1--;
			temp_vert2--;
			graph->list_edges_var[graph->beg_list_edges_var[temp_vert1] + graph->n_con_edges[temp_vert1]] = j;
			graph->list_edges_var[graph->beg_list_edges_var[temp_vert2] + graph->n_con_edges[temp_vert2]] = j;
			graph->list_nodes_var[graph->beg_list_edges_var[temp_vert1] + graph->n_con_edges[temp_vert1]] = temp_vert2;
			graph->list_nodes_var[graph->beg_list_edges_var[temp_vert2] + graph->n_con_edges[temp_vert2]] = temp_vert1;
			graph->n_con_edges[temp_vert1]++;
			graph->n_con_edges[temp_vert2]++;
		}
	}

	
	fclose(file);
}

void init(graph_s* graph, history_s *history, config_s *config, test_params_s *params, test_history_s *test_history) {
	initRandom();

	for (int i = 0; i < graph->n_vertices; i++) {
		(config->last_used)[i] = -3 * graph->n_vertices;
	}

	for (int i = 0; i < graph->n_vertices; i++) {
		history->hashes[i] = getStandardRandom() * 100;
	}

	//Generating random starting vector x;
	for (int i = 0; i < graph->n_vertices; i++) {
		if (getStandardRandom() > 0.5) {
			config->x[i] = 1;
			history->best_x[i] = 1;
		}
		else {
			config->x[i] = 0;
			history->best_x[i] = 0;
		}
	}

	config->f = F(graph, config->x, graph->edges);
	history->best_f = config->f;
	config->x_key = getKey(history, graph, config);

	config->step = 0;
	history->best_f_step = 0;
	history->stored_f_count = 0;
	history->tabu_change_t = 0;
	params->dateTimePrefix = getDateTimePrefix();
	params->start = clock();
	params->current = params->start;
	params->finish = params->start + CLOCKS_PER_SEC * (params->test_time_seconds);

}

//TODO: rename and implement function to clear history values between iterations

#pragma endregion

#pragma region ALGORITHMS

void tabu(graph_s* graph, history_s* history, config_s* config, test_params_s* params, test_history_s* test_history){
	long best_neighbour_index = -1, best_neighbour_delta_value = LONG_MIN;
	long best_neighbour_ac_index = -1;//aspiration criteria: if move can increase best solution but inside tabu : do it.
	long best_neighbour_ac_delta_value = LONG_MIN;

	config->step++;

	getNeighbourhood(graph,config->neighbourhood, config->last_used, config->step, config->tabu_size);

#pragma omp parallel for num_threads(16)
	for (int i = 0; i < graph->n_vertices; i++) {
		config->neighbourhood_values[i] = Fdelta(config->x, i, graph);

		if ( (config->neighbourhood)[i]) {
			if ((config->neighbourhood_values)[i] > best_neighbour_delta_value)
			{
				best_neighbour_index = i;
				best_neighbour_delta_value = (config->neighbourhood_values)[i];
			}
		}
		else { //check aspiration criteria
			if ((config->neighbourhood_values)[i] > best_neighbour_ac_delta_value)
			{
				best_neighbour_ac_index = i;
				best_neighbour_ac_delta_value = (config->neighbourhood_values)[i];
			}
		}
	}


	if (config->f + best_neighbour_ac_delta_value > history->best_f)
	{
		if (best_neighbour_delta_value > best_neighbour_ac_delta_value) 
		{
			config->f = config->f + best_neighbour_delta_value;
			(config->x)[best_neighbour_index] = !( (config->x)[best_neighbour_index] );
			config->last_used[best_neighbour_index] = config->step;
			history->index_count[best_neighbour_index]++;
		}
		else 
		{
			config->f = config->f + best_neighbour_ac_delta_value;
			config->x[best_neighbour_ac_index] = !config->x[best_neighbour_ac_index];
			config->last_used[best_neighbour_ac_index] = config->step;
			history->index_count[best_neighbour_ac_index]++;
		}
		
	}
	else
	{
		config->f = config->f + best_neighbour_delta_value;
		config->x[best_neighbour_index] = !((config->x)[best_neighbour_index]);
		config->last_used[best_neighbour_index] = config->step;
		history->index_count[best_neighbour_index]++;
	}

	
	if (config->f + best_neighbour_delta_value > history->best_f)
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
		history->best_f_step = config->step;
		history->best_f = config->f + best_neighbour_delta_value;

		log_f_value(history->best_f, config->step, params);

		for (int ii = 0; ii < graph->n_vertices; ii++)
			history->best_x[ii] = config->x[ii];
	}

	params->current = clock();
}

void reactive_tabu(graph_s* graph, history_s* history, config_s* config, test_params_s* params, test_history_s* test_history, reactive_s* reactive)
{
	long best_neighbour_index = -1, best_neighbour_delta_f = LONG_MIN, node_index = 0, deltaVisitTime = 0;
	long best_neighbour_ac_index = -1;//aspiration criteria: if move can increase best solution but inside tabu : do it.
	long best_neighbour_ac_delta_value = LONG_MIN;
	double x_key;

	x_key = getKey(history, graph, config);

#pragma region Reaction
	if (node_index = check_solution(x_key, config->f, history) > 0)
	{
		history->previous_solution_found_count++;
		deltaVisitTime = config->step - history->occurence_time[node_index];

		history->occurence_time[node_index] = config->step;
		if (deltaVisitTime < graph->n_vertices)
		{
			history->tabu_change_t = config->step;
			Increase(&(config->tabu_size), graph->n_vertices);
		}
	}
	else
	{
		save_solution(x_key, config->f, history, config->step);

		if (config->step - history->tabu_change_t > 100)
		{
			history->tabu_change_t = config->step;
			Decrease(&(config->tabu_size), graph->n_vertices);
		}
	}
	reactive->size_statistics[config->tabu_size]++;
#pragma endregion

	getNeighbourhood(graph, config->neighbourhood, config->last_used, config->step, config->tabu_size);

	(config->step)++;

#pragma omp parallel for num_threads(16)
	for (int i = 0; i < graph->n_vertices; i++) {
		config->neighbourhood_values[i] = Fdelta(config->x, i, graph);

		if (config->neighbourhood[i]) {
			if (config->neighbourhood_values[i] > best_neighbour_delta_f)
			{
				best_neighbour_index = i;
				best_neighbour_delta_f = config->neighbourhood_values[i];
			}
		}
		else { //check aspiration criteria
			if (config->neighbourhood_values[i] > best_neighbour_ac_delta_value)
			{
				best_neighbour_ac_index = i;
				best_neighbour_ac_delta_value = config->neighbourhood_values[i];
			}
		}
	}

	if (config->f + best_neighbour_ac_delta_value > history->best_f)
	{
		if (best_neighbour_delta_f > best_neighbour_ac_delta_value)
		{
			config->f = config->f + best_neighbour_delta_f;
			config->x[best_neighbour_index] = !((config->x)[best_neighbour_index]);
			config->last_used[best_neighbour_index] = config->step;
			history->index_count[best_neighbour_index]++;
		}
		else
		{
			config->f = config->f + best_neighbour_ac_delta_value;
			config->x[best_neighbour_ac_index] = !((config->x)[best_neighbour_ac_index]);
			config->last_used[best_neighbour_ac_index] = config->step;
			history->index_count[best_neighbour_ac_index]++;
		}

	}
	else
	{
		config->f += best_neighbour_delta_f;
		config->x[best_neighbour_index] = !((config->x)[best_neighbour_index]);
		config->last_used[best_neighbour_index] = config->step;
		history->index_count[best_neighbour_index]++;
	}

	/*if (*step % 50 == 0){
		printf("Step: %ld \t F: %ld Best F: %ld  Key: %f  Occurence: %ld \t Tabu: %ld\n", *step, best_neighbour_value, *best_f, *x_key,occurence_time[node_index], *tabu_size);
	}*/

	//copyX(x, best_x, n_vertices);

	if (config->f + best_neighbour_delta_f > history->best_f)
	{
		history->best_f_step = config->step;
		history->best_f = config->f + best_neighbour_delta_f;

		log_f_value(history->best_f, config->step, params);

		for (int ii = 0; ii < graph->n_vertices; ii++)
			history->best_x[ii] = config->x[ii];

	}

	params->current = clock();
}

void ultra_reactive_tabu(graph_s* graph, history_s* history, config_s* config, test_params_s* params, test_history_s* test_history, reactive_s* reactive)
{
	long best_neighbour_index = -1, best_neighbour_delta_f = LONG_MIN;
	long best_neighbour_ac_index = -1;//aspiration criteria: if move can increase best solution but inside tabu : do it.
	long best_neighbour_ac_delta_f = LONG_MIN;
	long index, value;
	double x_key;

	x_key = getKey(history,graph,config);


	getNeighbourhoodUR(graph,config->neighbourhood, config->last_used, config->step, reactive->tabu_sizes);

	(config->step)++;

	#pragma omp parallel for num_threads(16)
	for (int i = 0; i < graph->n_vertices; i++) {
		config->neighbourhood_values[i] = Fdelta(config->x, i, graph);

		if (config->neighbourhood[i]) {
			if (config->neighbourhood_values[i] > best_neighbour_delta_f)
			{
				best_neighbour_index = i;
				best_neighbour_delta_f = config->neighbourhood_values[i];
			}
		}
		else { //check aspiration criteria
			if (config->neighbourhood_values[i] > best_neighbour_ac_delta_f)
			{
				best_neighbour_ac_index = i;
				best_neighbour_ac_delta_f = config->neighbourhood_values[i];
			}
		}
	}

	if (config->f + best_neighbour_ac_delta_f > history->best_f)
	{
		if (best_neighbour_delta_f > best_neighbour_ac_delta_f)
		{
			config->f += best_neighbour_delta_f;
			(config->x)[best_neighbour_index] = !((config->x)[best_neighbour_index]);
			(config->last_used)[best_neighbour_index] = config->step;
			(history->index_count)[best_neighbour_index]++;
			update_ur_tabu_size(graph, best_neighbour_index, config->step, reactive->tabu_sizes, history->index_count);

			if ( config->f + best_neighbour_delta_f > history->best_f)
			{
				history->best_f_step = config->step;
				history->best_f = config->f + best_neighbour_delta_f;

				log_f_value(history->best_f,config->step,params);

				for (int ii = 0; ii < graph->n_vertices; ii++)
					history->best_x[ii] = config->x[ii];
			}

		}
		else
		{
			config->f += best_neighbour_ac_delta_f;
			config->x[best_neighbour_ac_index] = !(config->x[best_neighbour_ac_index]);
			config->last_used[best_neighbour_ac_index] = config->step;
			history->index_count[best_neighbour_ac_index]++;
			update_ur_tabu_size(graph,best_neighbour_ac_index, config->step, reactive->tabu_sizes, history->index_count);

			if (config->f + best_neighbour_ac_delta_f > history->best_f)
			{
				history->best_f_step = config->step;
				history->best_f = config->f + best_neighbour_delta_f;

				log_f_value(history->best_f, config->step, params);

				for (int ii = 0; ii < graph->n_vertices; ii++)
					history->best_x[ii] = config->x[ii];
			}

		}

	}
	else
	{
		config->f += best_neighbour_delta_f;
		config->x[best_neighbour_index] = !((config->x)[best_neighbour_index]);
		config->last_used[best_neighbour_index] = config->step;
		history->index_count[best_neighbour_index]++;
		update_ur_tabu_size(graph,best_neighbour_index, config->step, reactive->tabu_sizes, history->index_count);

		if (config->f + best_neighbour_delta_f > history->best_f)
		{
			history->best_f_step = config->step;
			history->best_f = config->f + best_neighbour_delta_f;

			log_f_value(history->best_f, config->step, params);

			for (int ii = 0; ii < graph->n_vertices; ii++)
				history->best_x[ii] = config->x[ii];
		}
	}	

	params->current = clock();
}

//void ultra_reactive_tabu_noac(graph_s graph, long *f, long *x, double *x_key, long *step, long max_steps, long *tabu_size,
//	bool *neighbourhood, long *neighbourhood_values, long *last_used, long *best_f, long *best_x, long *best_f_step, long* tabu_change_t, long* ultra_reactive_tabu_sizes, long* ultra_reactive_index_count)
//{
//	long best_neighbour_index = -1, best_neighbour_delta_f = LONG_MIN;
//	long best_neighbour_ac_index = -1;//aspiration criteria: if move can increase best solution but inside tabu : do it.
//	long best_neighbour_ac_delta_f = LONG_MIN;
//
//	*x_key = getKey(x);
//
//
//	getNeighbourhoodUR(neighbourhood, last_used, *step, ultra_reactive_tabu_sizes);
//
//	*step = *step + 1;
//
//	for (int i = 0; i < graph.n_vertices; i++) {
//		neighbourhood_values[i] = Fdelta(x, i, graph);
//
//		if (neighbourhood[i]) {
//			if (neighbourhood_values[i] > best_neighbour_delta_f)
//			{
//				best_neighbour_index = i;
//				best_neighbour_delta_f = neighbourhood_values[i];
//			}
//		}	
//	}
//
//
//	*f = *f + best_neighbour_delta_f;
//	x[best_neighbour_index] = !x[best_neighbour_index];
//	last_used[best_neighbour_index] = *step;
//	ultra_reactive_index_count[best_neighbour_index]++;
//	update_ur_tabu_size(best_neighbour_index, *step, ultra_reactive_tabu_sizes, ultra_reactive_index_count);
//
//
//	if (*f + best_neighbour_delta_f > *best_f)
//	{
//		*best_f_step = *step;
//		*best_f = *f + best_neighbour_delta_f;
//
//		for (int ii = 0; ii < graph.n_vertices; ii++)
//			best_x[ii] = x[ii];
//	}
//}

#pragma endregion

#pragma region TESTS
void test_tabu(graph_s* graph, history_s *history, config_s *config, test_params_s *params, test_history_s *test_history) {

	FILE *file;
	char path[100];

	printf("====================TABU====================\n");
	printf("Tabu size: %ld\n", config->tabu_size);
	printf("Restarts:%d\n", params->iterations_count);

	if (params->is_time_stop_condition)
		printf("Time per iteration, s. : %ld\n", params->test_time_seconds);
	else
		printf("Steps per iteration : %ld\n", params->max_steps);

	test_history->start_time = getDateTimePrefix();

	create_history_folder_paths(graph, test_history, "tabu");

	//indicative file
	strcpy_s(path, test_history->test_folder_path);
	strcat(path, "\\");
	strcat(path, graph_name(graph->g_number));
	strcat(path, ".name");
	fopen_s(&file, path, "w");
	fclose(file);

	strcpy_s(path, test_history->test_folder_path);
	strcat(path, "\\");
	strcat(path, "log.txt");

	fopen_s(&file, path, "a");

	fprintf_s(file, "====================TABU====================\n");
	fprintf_s(file, "Tabu size: %ld\n", config->tabu_size);
	fprintf_s(file, "Restarts:%d\n", params->iterations_count);

	if (params->is_time_stop_condition)
		fprintf_s(file, "Time per iteration, s. : %ld\n", params->test_time_seconds);
	else
		fprintf_s(file, "Steps till restart : %ld\n", params->max_steps);

	fprintf_s(file, "\nIndex\tF\tSteps\tTime\n");
	fclose(file);

	for (int i = 0; i < params->iterations_count; i++) {
		printf("Restart %d : \t", i + 1);

		init(graph,history,config,params,test_history);

#pragma region Main_Loop

		if (params->is_time_stop_condition) {
			while (runCondition(params))
				tabu(graph, history, config, params, test_history);
		}
		else {
			while (runCondition_steps(config->step, history->max_steps))
				tabu(graph, history, config, params, test_history);
		}

		(test_history->elapsed_time)[i] = (clock() - params->start) / CLOCKS_PER_SEC;
		test_history->iteration_steps[i] = config->step;
		test_history->iteration_bestf[i] = history->best_f;
#pragma endregion

#pragma region Results
		printf("Best F:%ld\tIterations: %ld\tTime: %f\n", history->best_f, test_history->iteration_steps[i], test_history->elapsed_time[i]);


		if (history->best_f > test_history->best_known_f) {
			test_history->best_known_f = history->best_f;
			if (file != NULL)
				fclose(file);

			strcpy_s(path, BEST_KNOWN);
			strcat(path, graph_name(graph->g_number));
			strcat(path, ".txt");

			fopen_s(&file, path, "w");
			fprintf_s(file, "%ld", history->best_f);
			fclose(file);
		}

		printResultToFile(graph->n_vertices,history->best_f,history->best_x,i, test_history->test_folder_path);

		strcpy_s(path, test_history->test_folder_path);
		strcat(path, "\\");
		strcat(path, "log.txt");
		appendIterationResult(i + 1, history->best_f, test_history->iteration_steps[i], test_history->elapsed_time[i], path);
		
		graph_move_index_log(graph->n_vertices, i, history->index_count, test_history->test_folder_path);
		iteration_cleanup(graph, history, config, params, test_history);
#pragma endregion


	}

}

void test_tabu_reactive(graph_s* graph, history_s *history, config_s *config, test_params_s *params, test_history_s *test_history, reactive_s* reactive) {
	FILE* file;
	char path[100];

	printf("====================REACTIVE====================\n");
	printf("Graph: G%d\n", graph->g_number);
	printf("Restarts:%d\n", params->iterations_count);

	if (params->is_time_stop_condition)
		printf("Time per iteration, s. : %ld\n", params->test_time_seconds);
	else
		printf("Steps till restart : %ld\n", params->max_steps);

	test_history->start_time= getDateTimePrefix();
	
	create_history_folder_paths(graph, test_history, "reactive");

	//indicative file
	strcpy_s(path, test_history->test_folder_path);
	strcat(path, "\\");
	strcat(path, graph_name(graph->g_number));
	strcat(path, ".name");
	fopen_s(&file, path, "w");
	fclose(file);

	strcpy_s(path, test_history->test_folder_path);
	strcat(path, "\\");
	strcat(path, "log.txt");

	fopen_s(&file, path, "a");

	fprintf_s(file, "====================REACTIVE====================\n");
	fprintf_s(file,"Graph: G%d\n", graph->g_number);
	fprintf_s(file, "Restarts:%d\n", params->iterations_count);

	if (params->is_time_stop_condition)
		fprintf_s(file, "Time per iteration, s. : %ld\n", params->test_time_seconds);
	else
		fprintf_s(file, "Steps till restart : %ld\n", params->max_steps);

	fprintf_s(file, "\nIndex\tF\tSteps\tTime\n");
	fclose(file);

	for (int i = 0; i < params->iterations_count; i++) {
		printf("Restart %d : \t", i + 1);

		init(graph, history, config, params, test_history);

#pragma region Main_Loop
		if (params->is_time_stop_condition) {
			while (runCondition(params))
				reactive_tabu(graph, history, config, params, test_history, reactive);
		}
		else {
			while (runCondition_steps(config->step, history->max_steps))
				reactive_tabu(graph, history, config, params, test_history, reactive);
		}

		test_history->elapsed_time[i] = (clock() - params->start) / (double)CLOCKS_PER_SEC;
		test_history->iteration_steps[i] = config->step;
		test_history->iteration_bestf[i] = history->best_f;
#pragma endregion

#pragma region Results
		
		printf("Best F:%ld\tIterations: %ld\tTime: %f\n", history->best_f, test_history->iteration_steps[i], test_history->elapsed_time[i]);

		if (history->best_f > test_history->best_known_f) {
			test_history->best_known_f = history->best_f;
			if (file != NULL)
				fclose(file);

			strcpy_s(path, BEST_KNOWN);
			strcat(path, graph_name(graph->g_number));
			strcat(path, ".txt");

			fopen_s(&file, path, "w");
			fprintf_s(file, "%ld", history->best_f);
			fclose(file);
		}

		printResultToFile(graph->n_vertices, history->best_f,history->best_x,i, test_history->test_folder_path);

		//strcpy_s(path, "d:\\data_maxcut\\results\\reactive_tabu\\log-");
		//strcat_s(path, params->dateTimePrefix);
		//strcat(path, ".txt");

		strcpy_s(path, test_history->test_folder_path);
		strcat(path, "\\");
		strcat(path, "log.txt");

		appendIterationResult(i + 1, history->best_f,test_history->iteration_steps[i],test_history->elapsed_time[i], path);
		graph_move_index_log(graph->n_vertices,i,history->index_count,test_history->test_folder_path);
		tabu_size_log_reactive(i,params->iterations_count,config->step, graph->n_vertices,reactive->size_statistics, test_history->test_folder_path, history->previous_solution_found_count);
		iteration_cleanup(graph, history, config, params, test_history);
#pragma endregion

	}

#pragma endregion

}

void test_tabu_ultra_reactive(graph_s* graph, history_s *history, config_s *config, test_params_s *params, test_history_s *test_history, reactive_s* reactive) {
	FILE* file;
	char path[100];
	char intstr[4];

	printf("====================ULTRA_REACTIVE====================\n");
	printf("Graph: G%d\n", graph->g_number);
	printf("Restarts:%d\n", params->iterations_count);
	printf("UR_Constant = %ld\n", UR_multiplier);

	if (params->is_time_stop_condition)
		printf("Time per iteration, s. : %ld\n", params->test_time_seconds);
	else
		printf("Steps till restart : %ld\n", params->max_steps);

	test_history->start_time = getDateTimePrefix();

	create_history_folder_paths(graph, test_history, "ultra");

	strcpy_s(path, test_history->test_folder_path);
	strcat(path, "\\");
	strcat(path, graph_name(graph->g_number));
	strcat(path, ".name");
	fopen_s(&file, path, "w");
	fclose(file);

	
	
	strcpy_s(path, test_history->test_folder_path);
	strcat(path, "\\");
	strcat(path, "log.txt");

	fopen_s(&file, path, "a");

	fprintf_s(file, "====================ULTRA_REACTIVE====================\n");
	fprintf_s(file, "Graph: G%d\n", graph->g_number);
	fprintf_s(file, "Restarts:%d\n", params->iterations_count);
	fprintf_s(file, "UR_Constant = %ld\n", UR_multiplier);

	if (params->is_time_stop_condition)
		fprintf_s(file, "Time per iteration, s. : %ld\n", params->test_time_seconds);
	else
		fprintf_s(file, "Steps till restart : %ld\n", params->max_steps);

	fprintf_s(file, "\nIndex\tF\tSteps\tTime\n");
	fclose(file);

	for (int i = 0; i < params->iterations_count; i++) {
		printf("Restart %d : \t", i + 1);

		init(graph, history, config, params, test_history);

		strcpy_s(path, test_history->test_folder_path);
		strcat(path, "\\");
		strcat(path, "f-");
		sprintf(intstr, "%ld", i);
		strcat(path, intstr);
		strcat(path, ".txt");
		fopen_s(&file, path, "w");
		fclose(file);
		params->f_file_path = new char[200];
		strcpy(params->f_file_path, path);

#pragma region Main_Loop

		if (params->is_time_stop_condition) {
			while (runCondition(params))
				ultra_reactive_tabu(graph, history, config, params, test_history, reactive);
		}
		else {
			while (runCondition_steps(config->step, history->max_steps))
				ultra_reactive_tabu(graph, history, config, params, test_history, reactive);
		}

		test_history->elapsed_time[i] = (clock() - params->start) / CLOCKS_PER_SEC;
		test_history->iteration_steps[i] = config->step;
		test_history->iteration_bestf[i] = history->best_f;
#pragma endregion

#pragma region Results

		printf("Best F:%ld\tIterations: %ld\tTime: %f\n", history->best_f, test_history->iteration_steps[i], test_history->elapsed_time[i]);

		if (history->best_f > test_history->best_known_f) {
			test_history->best_known_f = history->best_f;
			if (file != NULL)
				fclose(file);

			strcpy_s(path, BEST_KNOWN);
			strcat(path, graph_name(graph->g_number));
			strcat(path, ".txt");

			fopen_s(&file, path, "w");
			fprintf_s(file, "%ld", history->best_f);
			fclose(file);
		}

		printResultToFile(graph->n_vertices, history->best_f, history->best_x, i, test_history->test_folder_path);

		strcpy_s(path, test_history->test_folder_path);
		strcat(path, "\\");
		strcat(path, "log.txt");

		appendIterationResult(i + 1, history->best_f, test_history->iteration_steps[i], test_history->elapsed_time[i], path);

		//graph_move_index_log_ultra_reactive(graph->n_vertices, config->step, history->index_count); ??? IMPLEMENT  OR REMOVE
		iteration_cleanup(graph, history, config, params, test_history, reactive);
#pragma endregion

	}

#pragma endregion

}
#pragma endregion

void main()
{
	FILE *file;

	for (int g = 1; g <= 54; g++) {

		graph_s* graph = &graph_s();
		history_s* history = &history_s();
		config_s* config = &config_s();
		test_params_s* params = &test_params_s();
		test_history_s* test_history = &test_history_s();
		reactive_s* reactive = &reactive_s();

		fetch_properties(graph, &(test_history->best_known_f), g);

		if (!allocateMemory(graph, history, config, params, test_history, reactive))
			goto NOT_ENOUGH_FREE_MEMORY;

		fetch(graph);

		test_tabu_ultra_reactive(graph, history, config, params, test_history, reactive);

		cleanup(graph, history, config, params, test_history, reactive);
	}

	
NOT_ENOUGH_FREE_MEMORY :
	system("pause");			
}
