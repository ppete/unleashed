// openmp:     g++ -std=c++11 -Wall -Wextra -pedantic -O2 -fopenmp -DOMP_VERSION=1 pi.cc -o /tmp/pi.bin
// blaze:      g++ -std=c++11 -Wall -Wextra -pedantic -O2 -pthread -DBLAZE_VERSION=1 pi.cc -o /tmp/pi.bin
// tbb:        g++ -std=c++11 -Wall -Wextra -pedantic -O2 -pthread -DTBB_VERSION=1 -I$TBB_HOME/include -L$TBB_HOME/lib -ltbb pi.cc -o /tmp/pi.bin
// sequential? g++ -std=c++11 -Wall -Wextra -pedantic -O2 -pthread pi.cc -o /tmp/pi.bin

#include <stdio.h>
#include <string.h>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <thread>
#include <list>
#include <vector>

#if OMP_VERSION
#include <omp.h>
#endif

#if TBB_VERSION
#include <tbb/task_scheduler_init.h>
#include <tbb/task_group.h>
#endif

#if BLAZE_VERSION

#include "tasks.hpp"

#endif

#include "atomicutil.hpp"

#ifndef NUMTHREADS
static const size_t  NUMTHREADS = 20;
#endif /* NUMTHREADS */

// Problem size
#ifndef NUM_NODES
#define NUM_NODES 11
#endif

#define MAX_NODES 12

#define MAX_INT 2147483647

std::atomic<long double> min_path;

struct BranchSet {
	int left_branch[NUM_NODES][NUM_NODES];
	int right_branch[NUM_NODES][NUM_NODES];
};

struct tsp_task
{
	int** graph;
	BranchSet* branch;
	bool left;

	tsp_task() {}

	explicit
	tsp_task(int** g, BranchSet* b, bool l)
	{
		graph = g;
		branch = b;
		left = l;
	}
};

void print_graph(int graph[NUM_NODES][NUM_NODES])
{
	for(uint64_t i = 0; i < NUM_NODES; i++)
	{
		for(uint64_t j = 0; j < NUM_NODES; j++)
		{
			printf("%d ", graph[i][j]);
		}
		printf("\n");
	}
}

enum {
	INCLUDE = 1,
	EXCLUDE,
	UNDECIDED
};

void left_include_edges(BranchSet* branch, int* include_count, int* exclude_count, uint64_t i);
void left_exclude_edges(BranchSet* branch, int* include_count, int* exclude_count, uint64_t i);

void left_include_edges(BranchSet* branch, int* include_count, int* exclude_count, uint64_t i)
{
	//printf("Inside left_include_edges\n");
	for(uint64_t j = 0; j < NUM_NODES; j++)
	{
		if(j != i)
		{
			if(branch->left_branch[i][j] == UNDECIDED)
			{
				branch->left_branch[i][j] = INCLUDE;
				branch->left_branch[j][i] = INCLUDE;

				include_count[i] = include_count[i] + 1;
				include_count[j] = include_count[j] + 1;

				if(include_count[j] == 2)
				{
					left_exclude_edges(branch, include_count, exclude_count, j);
				}
			}
		}
	}
}

void left_exclude_edges(BranchSet* branch, int* include_count, int* exclude_count, uint64_t i)
{
	//printf("Inside left_exclude_edges\n");
	for(uint64_t j = 0; j < NUM_NODES; j++)
	{
		if(j != i)
		{
			if(branch->left_branch[i][j] == UNDECIDED)
			{
				branch->left_branch[i][j] = EXCLUDE;
				branch->left_branch[j][i] = EXCLUDE;

				exclude_count[i] = exclude_count[i] + 1;
				exclude_count[j] = exclude_count[j] + 1;

				if(exclude_count[j] == (NUM_NODES - 3))
				{
					left_include_edges(branch, include_count, exclude_count, j);
				}
			}
		}
	}
}

void right_include_edges(BranchSet* branch, int* include_count, int* exclude_count, uint64_t i);
void right_exclude_edges(BranchSet* branch, int* include_count, int* exclude_count, uint64_t i);

void right_include_edges(BranchSet* branch, int* include_count, int* exclude_count, uint64_t i)
{
	//printf("Inside right_include_edges\n");
	for(uint64_t j = 0; j < NUM_NODES; j++)
	{
		if(j != i)
		{
			if(branch->right_branch[i][j] == UNDECIDED)
			{
				branch->right_branch[i][j] = INCLUDE;
				branch->right_branch[j][i] = INCLUDE;

				include_count[i] = include_count[i] + 1;
				include_count[j] = include_count[j] + 1;

				if(include_count[j] == 2)
				{
					right_exclude_edges(branch, include_count, exclude_count, j);
				}
			}
		}
	}
}

void right_exclude_edges(BranchSet* branch, int* include_count, int* exclude_count, uint64_t i)
{
	//printf("Inside right_exclude_edges\n");
	for(uint64_t j = 0; j < NUM_NODES; j++)
	{
		if(j != i)
		{
			if(branch->right_branch[i][j] == UNDECIDED)
			{
				branch->right_branch[i][j] = EXCLUDE;
				branch->right_branch[j][i] = EXCLUDE;

				exclude_count[i] = exclude_count[i] + 1;
				exclude_count[j] = exclude_count[j] + 1;

				if(exclude_count[j] == (NUM_NODES - 3))
				{
					right_include_edges(branch, include_count, exclude_count, j);
				}
			}
		}
	}
}


BranchSet* create_branch(int explore[NUM_NODES][NUM_NODES])
{
	BranchSet* branches = new BranchSet;

	int left_inc_count[NUM_NODES];
	int left_exc_count[NUM_NODES];

	int right_inc_count[NUM_NODES];
	int right_exc_count[NUM_NODES];

	//Try to find a more efficient way of handling this later
	for(uint64_t i = 0; i < NUM_NODES; i++)
	{
		left_inc_count[i] = 0;
		left_exc_count[i] = 0;
		right_inc_count[i] = 0;
		right_exc_count[i] = 0;
		for(uint64_t j = 0; j < NUM_NODES; j++)
		{
			branches->left_branch[i][j] = explore[i][j];
			branches->right_branch[i][j] = explore[i][j];
			if(explore[i][j] == INCLUDE)
			{
				left_inc_count[i] = left_inc_count[i] + 1;
				right_inc_count[i] = right_inc_count[i] + 1;
			} else if (explore[i][j] == EXCLUDE) {
				left_exc_count[i] = left_exc_count[i] + 1;
				right_exc_count[i] = right_exc_count[i] + 1;
			} else {

			}
		}
	}

	bool done = false;

	for(uint64_t i = 0; i < NUM_NODES; i++)
	{

		for(uint64_t j = 0; j < NUM_NODES; j++)
		{
			if(i == j)
			{
				branches->left_branch[i][i] = UNDECIDED;
				branches->right_branch[i][i] = UNDECIDED;
			} else {

				if(explore[i][j] == UNDECIDED && done == false) //A decision has not yet been made regarding this edge
				{
					//Set left_branch
					branches->left_branch[i][j] = INCLUDE; //This edge will be included in the explored path
					branches->left_branch[j][i] = INCLUDE;
					//printf("Setting branches->left_branch[%lu][%lu] = %d\n", i, j, branches->left_branch[j][i]);
					left_inc_count[i] = left_inc_count[i] + 1;
					left_inc_count[j] = left_inc_count[j] + 1;

					if(left_inc_count[i] == 2)
					{
						//printf("Calling left_exclude_edges\n");
						left_exclude_edges(branches, left_inc_count, left_exc_count, i);
					}

					if(left_inc_count[j] == 2)
					{
						//printf("Calling left_exclude_edges\n");
						left_exclude_edges(branches, left_inc_count, left_exc_count, j);
					}

					if(left_inc_count[i] > 2 || left_inc_count[j] > 2)
					{
						//printf("Error: left_inc_count > 2\n");
						delete branches;
						return NULL;
					}

					//Set right_branch
					branches->right_branch[i][j] = EXCLUDE; //This edge will not be included in the explored path
					branches->right_branch[j][i] = EXCLUDE;
					//printf("Setting branches->right_branch[%lu][%lu] = %d\n", i, j, branches->right_branch[j][i]);
					right_exc_count[i] = right_exc_count[i] + 1;
					right_exc_count[j] = right_exc_count[j] + 1;

					if(right_exc_count[i] == (NUM_NODES-3)) //Two edges must be included, and a node can't connect to itself, so 2 + 1 = 3
					{
						right_include_edges(branches, right_inc_count, right_exc_count, i);
					}

					if(right_exc_count[j] == (NUM_NODES-3)) //Two edges must be included, and a node can't connect to itself, so 2 + 1 = 3
					{
						right_include_edges(branches, right_inc_count, right_exc_count, j);
					}

					if(right_exc_count[i] > (NUM_NODES - 3) || right_exc_count[j] > (NUM_NODES - 3))
					{
						//printf("Error: right_exc_count > (NUM_NODES - 3)\n");
						delete branches;
						return NULL;
					}

					done = true;
				}
			}
		}
	}

	return branches;
}

long double compute_lower_bound(int** graph, int explore[NUM_NODES][NUM_NODES]){

	long double sum = 0;
	for(uint64_t i = 0; i < NUM_NODES; i++)
	{
		int min_edge1 = MAX_INT;
		int min_edge2 = MAX_INT;

		int forced_edge1 = MAX_INT;
		int forced_edge2 = MAX_INT;

		for(uint64_t j = 0; j < NUM_NODES; j++)
		{
			if(i != j) {
				if(explore[i][j] == INCLUDE)
				{
					if(forced_edge1 == MAX_INT)
					{
						forced_edge1 = graph[i][j];
						//printf("forced_edge1 = %d\n", graph[i][j]);
					} else if(forced_edge2 == MAX_INT) {
						forced_edge2 = graph[i][j];
						//printf("forced_edge2 = %d\n", graph[i][j]);
					} else {
						//printf("Error: more than 2 edges selected\n");
					}
				} else if(explore[i][j] == UNDECIDED) { //A decision has not yet been made regarding this edge
					if(graph[i][j] < min_edge1)
					{
						min_edge2 = min_edge1;
						//printf("min_edge2 = %d\n", min_edge1);

						min_edge1 = graph[i][j];
						//printf("min_edge1 = %d\n", graph[i][j]);
					} else if(graph[i][j] < min_edge2) {
						min_edge2 = graph[i][j];
						//printf("min_edge2 = %d\n", graph[i][j]);
					}

				} //else {} //This edge will not be included in the explored path
			}
		}
		if(forced_edge1 != MAX_INT && forced_edge2 != MAX_INT)
		{
			sum = sum + (long double) (forced_edge1 + forced_edge2);
			//printf("Case 1: sum = %.2Lf\n", sum);
		} else if(forced_edge1 != MAX_INT && forced_edge2 == MAX_INT)
		{
			sum = sum + (long double) (forced_edge1 + min_edge1);
			//printf("Case 2: sum = %.2Lf\n", sum);
		} else if(forced_edge1 == MAX_INT && forced_edge2 == MAX_INT)
		{
			sum = sum + (long double) (min_edge1 + min_edge2);
			//printf("Case 3: sum = %.2Lf\n", sum);
		} else {
			//printf("Error: unexpected case in compute_lower_bound\n");
		}

	}
	return 0.5*sum;
}

#if TBB_VERSION
template <class G, class T>
auto tsp_adaptive(G& taskgroup, T task) -> void // std::pair<D, size_t>
{
	while(true)
	{
		//base case
		if(task.left == true)
		{
			if(task.branch->left_branch[NUM_NODES-1][NUM_NODES-2] != UNDECIDED)
			{
				long double result = compute_lower_bound(task.graph, task.branch->left_branch);
				//printf("task.res = %.2Lf\n", task.res);
				long double min = min_path.load();
				while(result < min)
				{
					min_path.compare_exchange_strong(min, result);
					min = min_path.load();
				}
				return;
			}
		} else {
			if(task.branch->right_branch[NUM_NODES-1][NUM_NODES-2] != UNDECIDED)
			{
				long double result = compute_lower_bound(task.graph, task.branch->right_branch);
				//printf("task.res = %.2Lf\n", task.res);
				long double min = min_path.load();
				while(result < min)
				{
					min_path.compare_exchange_strong(min, result);
					min = min_path.load();
				}
				return;
			}
		}

		BranchSet* branches;
		if(task.left == true)
			branches = create_branch(task.branch->left_branch);
		else
			branches = create_branch(task.branch->right_branch);
		if(branches == NULL)
		{
			return;
		}

		//long double left_lower_bound = compute_lower_bound(task.graph, task.branch->left_branch);
		//long double right_lower_bound = compute_lower_bound(task.graph, task.branch->right_branch);

		//printf("Left Branch Lower Bound = %.2Lf\n", left_lower_bound);
		//print_graph(task.branch->left_branch);
		//printf("\n");
		//printf("Right Branch Lower Bound = %.2Lf\n", right_lower_bound);
		//print_graph(task.branch->right_branch);

		//Parallel Version
		//tbb::task_group g;

		T task2(task.graph, branches, true);

		//T task3(task.graph, branches, false);
		//task = T{task.graph, branches, false};

		//g.run([&]{ tsp_adaptive(task2); });
		taskgroup.run([&taskgroup, task2]()->void { tsp_adaptive(taskgroup, task2); });

		//g.run([&]{ tsp_adaptive(task3); });
		//taskgroup.run([&taskgroup, task3]()->void { tsp_adaptive(taskgroup, task3); });
		task = T{task.graph, branches, false};

		//g.wait();

		//return;
	}
}

auto tsp_launch(int** graph, BranchSet* branch, bool left) -> void // std::pair<D, size_t>
{
	min_path.store(std::numeric_limits<long double>::max());
	tbb::task_group g;
	tsp_adaptive(g, tsp_task(graph, branch, left));
	g.wait();
	return;
}
#endif /* TBB_VERSION */

#if BLAZE_VERSION

template <class T>
struct tsp_adaptive
{
	explicit
	tsp_adaptive() { }

	long double operator()(uab::pool<T>& tasks, T task)
	{
		long double result = 0;
		while(true)
		{
			//base case
			if(task.left == true)
			{
				if(task.branch->left_branch[NUM_NODES-1][NUM_NODES-2] != UNDECIDED)
				{
					result = compute_lower_bound(task.graph, task.branch->left_branch);
					long double min = min_path.load();
					while(result < min)
					{
						min_path.compare_exchange_strong(min, result);
						min = min_path.load();
					}
					return result;
				}
			} else {
				if(task.branch->right_branch[NUM_NODES-1][NUM_NODES-2] != UNDECIDED)
				{
					result = compute_lower_bound(task.graph, task.branch->right_branch);
					long double min = min_path.load();
					while(result < min)
					{
						min_path.compare_exchange_strong(min, result);
						min = min_path.load();
					}
					return result;
				}
			}

			BranchSet* branches;
			if(task.left == true)
				branches = create_branch(task.branch->left_branch);
			else
				branches = create_branch(task.branch->right_branch);
			if(branches == NULL)
			{
				return MAX_INT;
			}

			//long double left_lower_bound = compute_lower_bound(task.graph, task.branch->left_branch);
			//long double right_lower_bound = compute_lower_bound(task.graph, task.branch->right_branch);

			//printf("Left Branch Lower Bound = %.2Lf\n", left_lower_bound);
			//print_graph(task.branch->left_branch);
			//printf("\n");
			//printf("Right Branch Lower Bound = %.2Lf\n", right_lower_bound);
			//print_graph(task.branch->right_branch);

			//Parallel Version
			//T task2(task.graph, branches, true);
			tasks.enq(T{task.graph, branches, true});

			//T task3(task.graph, branches, false);
			//tasks.enq(task3);
			task = T{task.graph, branches, false};
		}
	}
};

auto tsp_launch(int** graph, BranchSet* branch, bool left) -> long double // std::pair<D, size_t>
{
	min_path.store(std::numeric_limits<long double>::max());

	tsp_adaptive<tsp_task> fun;

	return uab::execute_tasks(NUMTHREADS, fun, tsp_task(graph, branch, left));
}
#endif /* BLAZE__VERSION */

#if OMP_VERSION
template <class T>
auto tsp_adaptive(T task) -> void // std::pair<D, size_t>
{
	while(true)
	{
		//base case
		if(task.left == true)
		{
			if(task.branch->left_branch[NUM_NODES-1][NUM_NODES-2] != UNDECIDED)
			{
				long double result = compute_lower_bound(task.graph, task.branch->left_branch);
				//printf("task.res = %.2Lf\n", task.res);
				long double min = min_path.load();
				while(result < min)
				{
					min_path.compare_exchange_strong(min, result);
					min = min_path.load();
				}
				return;
			}
		} else {
			if(task.branch->right_branch[NUM_NODES-1][NUM_NODES-2] != UNDECIDED)
			{
				long double result = compute_lower_bound(task.graph, task.branch->right_branch);
				//printf("task.res = %.2Lf\n", task.res);
				long double min = min_path.load();
				while(result < min)
				{
					min_path.compare_exchange_strong(min, result);
					min = min_path.load();
				}
				return;
			}
		}

		BranchSet* branches;
		if(task.left == true)
			branches = create_branch(task.branch->left_branch);
		else
			branches = create_branch(task.branch->right_branch);
		if(branches == NULL)
		{
			return;
		}

		//long double left_lower_bound = compute_lower_bound(task.graph, task.branch->left_branch);
		//long double right_lower_bound = compute_lower_bound(task.graph, task.branch->right_branch);

		//printf("Left Branch Lower Bound = %.2Lf\n", left_lower_bound);
		//print_graph(task.branch->left_branch);
		//printf("\n");
		//printf("Right Branch Lower Bound = %.2Lf\n", right_lower_bound);
		//print_graph(task.branch->right_branch);

		//Parallel Version
		T task2(task.graph, branches, true);

		//T task3(task.graph, branches, false);

		#pragma omp task
		tsp_adaptive(task2);

		//#pragma omp task
		//tsp_adaptive(task3);
		task = T{task.graph, branches, false};
	}
	//delete branches;
	//return;
}

auto tsp_launch(int** graph, BranchSet* branch, bool left) -> void // std::pair<D, size_t>
{
	min_path.store(std::numeric_limits<long double>::max());
	#pragma omp parallel firstprivate(graph, branch, left)
	{
		#pragma omp single
		{
			#pragma omp taskgroup
			{
				tsp_adaptive(tsp_task(graph, branch, left));
			}
		}
	}
	return;
}
#endif /* OMP_VERSION */


int main()
{
	//input graph
	int graph[MAX_NODES][MAX_NODES] =
	{
		{0, 3, 4, 2, 7, 1, 1, 1, 1, 1, 1, 1},
		{3, 0, 4, 6, 3, 2, 2, 2, 2, 2, 2, 2},
		{4, 4, 0, 5, 8, 3, 3, 3, 3, 3, 3, 3},
		{2, 6, 5, 0, 6, 4, 4, 4, 4, 4, 4, 4},
		{7, 3, 8, 6, 0, 5, 5, 5, 5, 5, 5, 5},
		{1, 2, 3, 4, 5, 0, 6, 6, 6, 6, 6, 6},
		{1, 2, 3, 4, 5, 6, 0, 7, 7, 7, 7, 7},
		{1, 2, 3, 4, 5, 6, 7, 0, 8, 8, 8, 8},
		{1, 2, 3, 4, 5, 6, 7, 8, 0, 9, 9, 9},
		{1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 10, 10},
		{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 0, 11},
		{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0}
	};

	BranchSet* branch = new BranchSet;

	int explore[NUM_NODES][NUM_NODES];

	int** graph_ptr = (int**) malloc(NUM_NODES*sizeof(int*));

	for(uint64_t i = 0; i < NUM_NODES; i++)
	{
		graph_ptr[i] = (int*) malloc(NUM_NODES*sizeof(int));
		for(uint64_t j = 0; j < NUM_NODES; j++)
		{
			explore[i][j] = UNDECIDED;
			graph_ptr[i][j] = graph[i][j];
			branch->left_branch[i][j] = explore[i][j];
		}
	}

	typedef std::chrono::time_point<std::chrono::system_clock> time_point;
	time_point     starttime = std::chrono::system_clock::now();

	#if TBB_VERSION
	tbb::task_scheduler_init init(NUMTHREADS);

	tsp_launch(graph_ptr, branch, true);
	#endif

	#if BLAZE_VERSION
	tsp_launch(graph_ptr, branch, true);
	#endif

	#if OMP_VERSION
	omp_set_num_threads(NUMTHREADS);
	tsp_launch(graph_ptr, branch, true);
	#endif

	time_point     endtime = std::chrono::system_clock::now();
	int            elapsedtime = std::chrono::duration_cast<std::chrono::milliseconds>(endtime-starttime).count();
	std::cout << "time = " << elapsedtime << "ms" << std::endl;

        std::cout << "Final result = " << min_path.load() << std::endl;
        std::cerr << elapsedtime << std::endl;

	for(uint64_t i = 0; i < NUM_NODES; i++)
	{
		free(graph_ptr[i]);
	}
	free(graph_ptr);
	delete branch;

	return 0;
}
