/*

*/

#include <stdio.h>
#include <cuda_runtime_api.h>
#include <time.h>
//#include <math.h>
#include <float.h>
//#include <cuda_runtime.h>
#include "unimem.h"
#include "first.h"

#define DEBUG_CUDA_INFO
#define MAX_RESULTS 100000
#define PACK_RESULTS

struct boundaries {
	coord x_min;
	coord y_min;
	coord x_max;
	coord y_max;
	bool intersection;
};

struct node* m_dev_node = NULL;
__constant__ boundaries dev_bonds[1];
__constant__ unsigned dev_threads_count[1];
unsigned m_threads_count;
indexer m_count_branches;
int m_length_of_tree = 0;

//__constant__ struct branch *m_ttt_cuda_first_branch = NULL;

#ifdef PACK_RESULTS
/// compare
int cmp(const void* a, const void* b)
{
	return (int)(*(indexer*)a - *(indexer*)b);
}

#endif

extern "C"
bool init_cuda_device(int deviceID, struct node* node)
{
	if (!node)
		return false;
	//return false;

	struct node *nd = node;
	unsigned count1[64], i = 0;
	m_count_branches = 0;
	for (int j = 0; j < 64; j++) count1[j] = 0;
	//count1[0] = 1;

	struct node *stack_node[64];
	int stack_pos = 0;
	indexer stack_idx[64];
	alignas(16) struct branch *first_branch = NULL;
	alignas(16) struct node* stack_first_node[64];
	for (unsigned i = 0; i < 64; ++i) {
		stack_first_node[i] = NULL;
	}
	while (i < nd->count_child_nodes) {
		if (!stack_first_node[stack_pos] || nd < stack_first_node[stack_pos])
			stack_first_node[stack_pos] = nd;
		if (nd->is_last_node) {
			for (unsigned j = 0; j < nd->count_child_nodes; ++j) {
				struct branch *br = (struct branch*)(nd->child_node[j]);
				if (!first_branch || br < first_branch)
					first_branch = br;
			}
			/*if (!count_br)
				count_br = nd->count_child_nodes;
			if (count_br != nd->count_child_nodes)
				printf("Branches %u vs %u\n", count_br, nd->count_child_nodes);*/
			m_count_branches += nd->count_child_nodes;
				// return from stack
			while (stack_pos > 0) {
				stack_pos--;
				nd = stack_node[stack_pos];
				i = stack_idx[stack_pos] + 1;

				if (i < nd->count_child_nodes) {
					stack_idx[stack_pos] = i;
					stack_node[stack_pos] = nd;
					stack_pos++;
					nd = (struct node*)nd->child_node[i];
					i = 0;
					break;
				}
				else {
					//if (count1[stack_pos]) {
						//if (count1[stack_pos] != nd->count_child_nodes) {
						//	printf("Nodes %u vs %u\n", count1[stack_pos], nd->count_child_nodes);
						//}
						count1[stack_pos] += nd->count_child_nodes;
					//}
					//else {
					//	count1[stack_pos] = nd->count_child_nodes;
					//}
				}
			}
		}
		else {
			stack_idx[stack_pos] = i;
			stack_node[stack_pos] = nd;
			stack_pos++;
			i = 0;
			nd = (struct node*)nd->child_node[i];
			/*if (!count1[stack_pos])
				count1[stack_pos] = nd->count_child_nodes;
			else
				count1[stack_pos] += nd->count_child_nodes;*/
		}
		/*} else if (i < nd->count_child_nodes) {
		i++;*/
		/*if (!count1[stack_pos])
			count1[stack_pos] = nd->count_child_nodes;
		else
			count1[stack_pos] += nd->count_child_nodes;
			*/
	}

	//return false;

	int deviceCount;
	cudaError_t er1 = cudaGetDeviceCount(&deviceCount);
	printf("DevicecheckCudaErrors Count: %d\n", deviceCount);

	if (deviceID == -1)
		deviceID = 0;

	cudaDeviceProp prop;
	for (int ii = 0; ii < deviceCount; ++ii) {
		er1 = cudaGetDeviceProperties(&prop, ii);
		if (prop.major < 2 || prop.canMapHostMemory != 1)
		{
			printf("ERROR: calculation requires GPU devices with compute SM 2.0 or higher, or can not using MapHostMemory.\n");
			printf("Current GPU device has compute SM%d.%d, Exiting...", prop.major, prop.minor);
			//exit(EXIT_WAIVED);
			return false;
		}

		printf("GPU device name is %s\n", prop.name);
		printf("GPU total memory = %.0f Mb\n", prop.totalGlobalMem / 1024.0 / 1024.0);
		printf("Number of multiprocessors on the device = %u\n", prop.multiProcessorCount);
	}

	er1 = cudaSetDevice(deviceID);
	cudaSetDeviceFlags(cudaDeviceMapHost);
	er1 = cudaGetDeviceProperties(&prop, deviceID);
	m_threads_count = prop.multiProcessorCount * prop.warpSize;
	//er1 = cudaMalloc((void**)&dev_threads_count, sizeof(unsigned));
	er1 = cudaMemcpyToSymbol(dev_threads_count, &m_threads_count, sizeof(unsigned));

	// copy rtree
	int pos = 63;

	for (; pos >= 0; --pos) {
		if (count1[pos])
			break;
	}
	m_length_of_tree = pos + 1;

	// allocationg memory for branches
	alignas(16) struct branch* tbr = (struct branch*)aligned_alloc(16, sizeof(struct branch) * m_count_branches);
	// struct branch* first_branch = NULL;
	/*nd = node;
	i = 0;
	unsigned k = 0;
	while (i < nd->count_child_nodes) {
		if (nd->is_last_node) {
			for (unsigned j = 0; j < nd->count_child_nodes; ++j) {
				struct branch *br = (struct branch*)(nd->child_node[j]);
				//if (!first_branch || br < first_branch)
				//	first_branch = br;
				memcpy(tbr + k, br, sizeof(struct branch));
				k++;
			}
			while (stack_pos > 0) {
				stack_pos--;
				nd = stack_node[stack_pos];
				i = stack_idx[stack_pos] + 1;

				if (i < nd->count_child_nodes) {
					// insert to stack
					stack_idx[stack_pos] = i;
					stack_node[stack_pos] = nd;
					stack_pos++;
					nd = (struct node*)nd->child_node[i];
					i = 0;
					break;
				}
			}
		}
		else {
			// insert to stack
			stack_idx[stack_pos] = i;
			stack_node[stack_pos] = nd;
			stack_pos++;
			nd = (struct node*)nd->child_node[i];
		}
	}*/
	/*nd = node;
	i = 0;
	for (int j = 0; j <= pos; ++j)
		nd = (struct node*)(nd->child_node[0]);
	for (indexer j = 0; j < count1[pos]; ++j) {
		for (indexer k = 0; k < nd[j].count_child_nodes; ++k) {
			memcpy(tbr + i, nd[j].child_node[k], sizeof(struct branch));
			i++;
		}
	}
	*/
	memcpy(tbr, first_branch, sizeof(struct branch) * m_count_branches);

	// for debug
	/*printf("\n\n\n======================================================================\n");
	for (indexer i = 0; i < count1[pos]; ++i) {
		if ((struct node*)(stack_first_node[pos + 1])[i].is_last_node) {
			unsigned tt = ((struct node*)(stack_first_node[pos + 1]))[i].count_child_nodes;
			for (indexer ii = 0; ii < tt; ++ii) {
				unsigned idx = (struct branch*)((struct node*)(stack_first_node[pos + 1])[i].child_node[ii]) - first_branch;
				if (!idx)
					printf("0\n");
				else
					printf("%u\n", idx);
			}
		}
		else {
			printf("Error last node %u\n", i);
		}
	}*/

	// copy data of branches to device
	clock_t t1 = clock();
	cudaStream_t stream;
	cudaStreamCreate(&stream);
	for (indexer i = 0; i < m_count_branches; ++i)
	{
		void *data_ptr = tbr[i].leaf_x;
		er1 = cudaMalloc((void**)&(tbr[i].leaf_x), sizeof(coord) * tbr[i].count_leafs);
		er1 = cudaMemcpyAsync(tbr[i].leaf_x, data_ptr, sizeof(coord) * tbr[i].count_leafs, cudaMemcpyHostToDevice, stream);
		data_ptr = tbr[i].leaf_y;
		er1 = cudaMalloc((void**)&(tbr[i].leaf_y), sizeof(coord) * tbr[i].count_leafs);
		er1 = cudaMemcpyAsync(tbr[i].leaf_y, data_ptr, sizeof(coord) * tbr[i].count_leafs, cudaMemcpyHostToDevice, stream);
		data_ptr = tbr[i].leaf_number;
		er1 = cudaMalloc((void**)&(tbr[i].leaf_number), sizeof(indexer) * tbr[i].count_leafs);
		er1 = cudaMemcpyAsync(tbr[i].leaf_number, data_ptr, sizeof(indexer) * tbr[i].count_leafs, cudaMemcpyHostToDevice, stream);
		data_ptr = tbr[i].merge_next_leaf;
		er1 = cudaMalloc((void**)&(tbr[i].merge_next_leaf), sizeof(bool) * tbr[i].count_leafs);
		er1 = cudaMemcpyAsync(tbr[i].merge_next_leaf, data_ptr, sizeof(bool) * tbr[i].count_leafs, cudaMemcpyHostToDevice, stream);
		/*data_ptr = tbr[i].xsh_min;
		er1 = cudaMalloc((void**)&(tbr[i].xsh_min), sizeof(coord) * tbr[i].count_shapes);
		er1 = cudaMemcpyAsync(tbr[i].xsh_min, data_ptr, sizeof(coord) * tbr[i].count_shapes, cudaMemcpyHostToDevice, stream);
		data_ptr = tbr[i].xsh_max;
		er1 = cudaMalloc((void**)&(tbr[i].xsh_max), sizeof(coord) * tbr[i].count_shapes);
		er1 = cudaMemcpyAsync(tbr[i].xsh_max, data_ptr, sizeof(coord) * tbr[i].count_shapes, cudaMemcpyHostToDevice, stream);
		data_ptr = tbr[i].ysh_min;
		er1 = cudaMalloc((void**)&(tbr[i].ysh_min), sizeof(coord) * tbr[i].count_shapes);
		er1 = cudaMemcpyAsync(tbr[i].ysh_min, data_ptr, sizeof(coord) * tbr[i].count_shapes, cudaMemcpyHostToDevice, stream);
		data_ptr = tbr[i].ysh_max;
		er1 = cudaMalloc((void**)&(tbr[i].ysh_max), sizeof(coord) * tbr[i].count_shapes);
		er1 = cudaMemcpyAsync(tbr[i].ysh_max, data_ptr, sizeof(coord) * tbr[i].count_shapes, cudaMemcpyHostToDevice, stream); */
		data_ptr = tbr[i].offset;
		er1 = cudaMalloc((void**)&(tbr[i].offset), sizeof(indexer) * tbr[i].count_shapes);
		er1 = cudaMemcpyAsync(tbr[i].offset, data_ptr, sizeof(indexer) * tbr[i].count_shapes, cudaMemcpyHostToDevice, stream);
	}
	er1 = cudaStreamSynchronize(stream);
	er1 = cudaStreamDestroy(stream);
	clock_t t2 = clock();
	printf("Time copying data to device = %u ms\n", t2 - t1);

	// copy branches to device
	struct branch *dev_br = NULL;
	er1 = cudaMalloc((void**)&dev_br, sizeof(struct branch) * m_count_branches);
	er1 = cudaMemcpy(dev_br, tbr, sizeof(struct branch) * m_count_branches, cudaMemcpyHostToDevice);
	//cudaMemcpyToSymbol(m_ttt_cuda_first_branch, &dev_br, sizeof(struct branch*));

	//return false;
	alignas(16) struct node *to_dev_nd[65];
	//void **to_dev_child[64];
	struct node *dev_nd = NULL, *dev_nd_prev = NULL, *dev_ptr = NULL;
	// to_dev_nd[0] = (struct node*)aligned_alloc(16, sizeof(struct node));
	// memcpy(to_dev_nd[0], nd, sizeof(struct node));
	struct node* tnd = node;
	//for (unsigned j = 0; j <= pos; ++j)
	//	tnd = (struct node*)(tnd->child_node[0]);
	//unsigned j = 0;
	//void* tmp1 = NULL;
	unsigned count = tnd->count_child_nodes, prev_count = 1;
	//printf("\n\n\n======================================================================\n");
	for (int k1 = pos; k1 >= 0; --k1) {
		tnd = node;
		//for (unsigned j = 0; j <= k1; ++j)
		//	tnd = (struct node*)(tnd->child_node[0]);
		// data child node
		to_dev_nd[k1] = (struct node*)aligned_alloc(16, sizeof(struct node) * count1[k1]);
		//memcpy(to_dev_nd[k1], tnd/*->child_node[0]*/, sizeof(struct node) * count1[k1]);
		memcpy(to_dev_nd[k1], stack_first_node[k1 + 1], sizeof(struct node) * count1[k1]);
		// pointer to child_node on host
		for (indexer j = 0; j < count1[k1]; ++j) {
			//(to_dev_nd[k1])[j]->child_node = (void**)aligned_alloc(16, sizeof(void*) * MAX_NODES); // tnd->count_child_nodes);
			//(to_dev_child[k1])[j] = to_dev_nd[j]->child_node;
			dev_ptr = NULL;
			er1 = cudaMalloc((void**)&dev_ptr, sizeof(void*) * MAX_NODES);
			(to_dev_nd[k1])[j].child_node = (void**)dev_ptr;
			for (indexer k2 = 0; k2 < MAX_NODES; ++k2) {
				if (k1 == pos) {
					// copy pointer of branches
					//struct branch *ptr = &(dev_br[k2 + j * MAX_NODES]);
					unsigned idx = (struct branch*)((struct node*)(stack_first_node[k1 + 1])[j].child_node[k2]) - first_branch;
					//if (idx == 4899)
					//	printf("%u\n", idx);
					struct branch *ptr = &(dev_br[idx]);
					er1 = cudaMemcpy((void*)((to_dev_nd[k1])[j].child_node + k2), &ptr, sizeof(struct branch*), cudaMemcpyHostToDevice);
				}
				else {
					// copy pointer of nodes
					//struct node* ptr = &(dev_nd_prev[k2 + j * MAX_NODES]);
					unsigned idx = (struct node*)(stack_first_node[k1 + 1])[j].child_node[k2] - (struct node*)(stack_first_node[k1 + 2]);
					//printf("%u\n", idx);
					struct node *ptr = &(dev_nd_prev[idx]);
					er1 = cudaMemcpy((void*)((to_dev_nd[k1])[j].child_node + k2), &ptr, sizeof(struct node*), cudaMemcpyHostToDevice);
				}
			}
		}
		//printf("==========================================\n\n\n");
		// pointers of child nodes
		er1 = cudaMalloc((void**)&dev_nd, sizeof(struct node) * count1[k1]); // tnd->count_child_nodes);
		cudaMemcpy(dev_nd, to_dev_nd[k1], sizeof(struct node) * count1[k1], cudaMemcpyHostToDevice);
		dev_nd_prev = dev_nd;
	}
	// copy top node (root)
	to_dev_nd[64] = (struct node*)aligned_alloc(16, sizeof(struct node));
	memcpy(to_dev_nd[64], node/*->child_node[0]*/, sizeof(struct node));
	dev_ptr = NULL;
	er1 = cudaMalloc((void**)&dev_ptr, sizeof(void*) * node->count_child_nodes);
	(to_dev_nd[64])[0].child_node = (void**)dev_ptr;
	for (indexer k2 = 0; k2 < node->count_child_nodes; ++k2) {
		// copy pointer of nodes
		//struct node* ptr = &(dev_nd_prev[k2]);
		unsigned idx = (struct node*)(stack_first_node[0])[0].child_node[k2] - (struct node*)(stack_first_node[1]);
		struct node* ptr = &(dev_nd_prev[idx]);
		er1 = cudaMemcpy((void*)((to_dev_nd[64])[0].child_node + k2), &ptr, sizeof(struct node*), cudaMemcpyHostToDevice);
	}
	// pointers of child nodes
	er1 = cudaMalloc((void**)&dev_nd, sizeof(struct node)); // tnd->count_child_nodes);
	er1 = cudaMemcpy(dev_nd, to_dev_nd[64], sizeof(struct node), cudaMemcpyHostToDevice);
	m_dev_node = dev_nd;
	printf("============== 0x%llx, 0x%llx, prev = 0x%llx\n", m_dev_node, dev_nd, dev_nd_prev);

	// free memory
	for (int k1 = pos; k1 >= 0; --k1) {
		_aligned_free(to_dev_nd[k1]);
	}
	_aligned_free(to_dev_nd[64]);

	// allocating memory for root
	//er1 = cudaMalloc((void**)&m_dev_node, sizeof(struct node));
	// copy to device root of tree
	//cudaMemcpy(m_dev_node, to_dev_nd[0], sizeof(struct node), cudaMemcpyHostToDevice);

	return true;
}

extern "C"
bool destroy_cuda_device()
{
	//cudaFree(dev_threads_count);
	cudaError_t er1 = cudaDeviceReset();
	return er1 == cudaSuccess ? true : false;
}

extern "C"
#if defined(CALC_CIRCLE) || defined(CALC_POINT)
/* searchin items in selected rectangle on cuda device */
indexer* cuda_search_rect2(struct node *nd, coord x_min, coord y_min, coord x_max, coord y_max, bool intersection, /*out*/indexer *count_items, ret_callback2_circle callback = NULL, void *data = NULL);
/* searchin items in selected rectangle on cuda device imlementation */
__global__ void cuda_search_rect2_impl1(void **nd, indexer *iter_count, indexer *atomic_iter, /*out*/ void **next_nd, /*out*/ indexer *idxs, /*out*/indexer *count_items, ret_callback2_circle callback = NULL, void *data = NULL);
__global__ void cuda_search_rect2_impl2(void **nd, int *iter_count, indexer *atomic_iter, /*out*/ void **next_nd, /*out*/ indexer *idxs, /*out*/indexer *count_items);
__global__ void cuda_search_rect2_impl3(indexer *unpack_idxs, indexer *iter_count, indexer *atomic_iter, /*out*/ void **next_nd, /*out*/ indexer *idxs, /*out*/indexer *count_items);
#else
indexer* search_rect2(struct node *nd, coord x_min, coord y_min, coord x_max, coord y_max, bool intersection, /*out*/indexer *count_items)
__global__ indexer* search_rect2_impl(void *nd_ptr, indexer iter_count, /*out*/indexer *count_items)
#endif // CALC_POINT

#if defined(CALC_CIRCLE) || defined(CALC_POINT)
/* searchin items in selected rectangle on cuda device */
indexer* cuda_search_rect2(node * nd, coord x_min, coord y_min, coord x_max, coord y_max, bool intersection, indexer * count_items, ret_callback2_circle callback, void * data)
#else
indexer* search_rect2(struct node *nd, coord x_min, coord y_min, coord x_max, coord y_max, bool intersection, /*out*/indexer *count_items)
#endif // CALC_POINT
{

	// memory for result
	/*size_t mem_size = 256;
	size_t count_mem = 1;
	alignas(16) indexer* idxs = (indexer*)aligned_alloc(16, sizeof(indexer) * mem_size * count_mem); */
	cudaError_t er1;
	cudaStream_t stream;
	cudaStreamCreate(&stream);

	indexer *host_idxs = NULL, *dev_idxs = NULL, *dev_tmp_idxs = NULL;;
	cudaHostAlloc((void**)&host_idxs, sizeof(indexer) * MAX_RESULTS, cudaHostAllocMapped);
	cudaHostGetDevicePointer((void**)&dev_idxs, host_idxs, 0);
	cudaMalloc((void**)&dev_tmp_idxs, sizeof(indexer) * MAX_RESULTS);

	// searching
	cudaEvent_t start, stop;
	float gtime = 0.0;
	int device_id;
	cudaDeviceProp prop;
	er1 = cudaGetDevice(&device_id);
	er1 = cudaGetDeviceProperties(&prop, device_id);
	dim3 grid_size = dim3(prop.multiProcessorCount, 1, 1), block_size = dim3(prop.warpSize, 1, 1);
	// store boundaries
	boundaries b1;
	b1.intersection = intersection; b1.x_max = x_max; b1.x_min = x_min; b1.y_max = y_max; b1.y_min = y_min;
	//cudaMalloc((void**)dev_bonds, sizeof(struct boundaries));
	cudaMemcpyToSymbolAsync(dev_bonds, &b1, sizeof(struct boundaries), 0, cudaMemcpyHostToDevice, stream);
	// for store count of iterations to next step
	indexer *dev_atomic_iter = NULL;
	cudaMalloc((void**)&dev_atomic_iter, sizeof(indexer));
	cudaMemsetAsync(dev_atomic_iter, 0, 1, stream);
	// store pointers for next step
	void **dev_ptr = NULL, **dev_ptr2 = NULL;
	cudaMalloc((void**)&dev_ptr, sizeof(void*) * m_count_branches);
	//printf("======================= 0x%llx; 0x%llx, count_br = %u\n", &m_dev_node, m_dev_node, m_count_branches);
	void **tptr = (void**)(&m_dev_node);
	cudaMemcpyAsync(dev_ptr, tptr, sizeof(void*), cudaMemcpyHostToDevice, stream);
	cudaMalloc((void**)&dev_ptr2, sizeof(void*) * m_count_branches);
	//printf("======================= 0x%llx; 0x%llx; dev_ptr = 0x%llx\n", &m_dev_node, m_dev_node, dev_ptr);
	// count items
	indexer *dev_count_items = NULL;
	cudaMalloc((void**)&dev_count_items, sizeof(indexer));
	// count of iterations
	indexer *dev_iter_count = NULL;
	cudaMalloc((void**)&dev_iter_count, sizeof(indexer));
	cudaMemsetAsync(dev_iter_count, 1, 1, stream);
	cudaStreamSynchronize(stream);

	indexer atomic_iter = 0;

	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, stream);
	
	// calculating nodes
	for (int i = 0; i < m_length_of_tree + 1; ++i) {
		er1 = cudaMemsetAsync(dev_atomic_iter, 0, sizeof(indexer), stream);
		cuda_search_rect2_impl1 << <grid_size, block_size, 0, stream >> > ((void**)dev_ptr, dev_iter_count, dev_atomic_iter, dev_ptr2, dev_idxs, count_items, callback, data);

		er1 = cudaMemcpyAsync(&atomic_iter, dev_atomic_iter, sizeof(indexer), cudaMemcpyDeviceToHost, stream);
		er1 = cudaMemcpyAsync(dev_ptr, dev_ptr2, sizeof(void*) * atomic_iter, cudaMemcpyDeviceToDevice, stream);
		er1 = cudaMemcpyAsync(dev_iter_count, dev_atomic_iter, sizeof(indexer), cudaMemcpyDeviceToDevice, stream);
		cudaStreamSynchronize(stream);
		printf("===== Iter %i: next = %u (%s)\n", i, atomic_iter, er1 == cudaSuccess ? "true" : "false");
		//cudaThreadSynchronize();
	}

	cudaEventRecord(stop, stream);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&gtime, start, stop);
	printf("Kernel 1 time = %f ms\n", gtime);

	// calculating branches
	grid_size = dim3(atomic_iter, 1, 1);
	cudaEventRecord(start, stream);
	er1 = cudaMemsetAsync(dev_atomic_iter, 0, sizeof(indexer), stream);
	cuda_search_rect2_impl2 << <grid_size, block_size, 0, stream >> > ((void**)dev_ptr, NULL, dev_atomic_iter, dev_ptr2, dev_idxs, count_items);
	er1 = cudaMemcpyAsync(&atomic_iter, dev_atomic_iter, sizeof(indexer), cudaMemcpyDeviceToHost, stream);
	cudaStreamSynchronize(stream);
	cudaThreadSynchronize();
	cudaEventRecord(stop, stream);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&gtime, start, stop);
	printf("Kernel 2 time = %f ms\n", gtime);

	/*grid_size = dim3(prop.multiProcessorCount, 1, 1);
	cudaEventRecord(start, stream);
	er1 = cudaMemcpyAsync(dev_iter_count, dev_atomic_iter, sizeof(indexer), cudaMemcpyDeviceToDevice, stream);
	er1 = cudaMemsetAsync(dev_atomic_iter, 0, sizeof(indexer), stream);
	cuda_search_rect2_impl3 << <grid_size, block_size, 0, stream >> > (dev_tmp_idxs, dev_iter_count, dev_atomic_iter, dev_ptr2, dev_idxs, count_items);
	cudaStreamSynchronize(stream);
	cudaThreadSynchronize();
	cudaEventRecord(stop, stream);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&gtime, start, stop);
	printf("Kernel 3 time = %f ms\n", gtime);*/

	cudaMemcpyAsync(count_items, dev_atomic_iter, sizeof(indexer), cudaMemcpyDeviceToHost, stream);

	indexer *idxs = (indexer*)aligned_alloc(16, sizeof(indexer) * *count_items);
	er1 = cudaMemcpyAsync(idxs, host_idxs, sizeof(indexer) * *count_items, cudaMemcpyHostToHost, stream);
	cudaStreamSynchronize(stream);

	// freeing and destroying
	cudaStreamDestroy(stream);

	er1 = cudaFree(dev_iter_count);
	er1 = cudaFree(dev_ptr);
	er1 = cudaFree(dev_ptr2);
	er1 = cudaFree(dev_tmp_idxs);
	er1 = cudaFree(dev_count_items);
	er1 = cudaFree(dev_atomic_iter);
	cudaEventDestroy(stop);
	cudaEventDestroy(start);
	cudaFreeHost(host_idxs);

#ifdef PACK_RESULTS
	if (*count_items) {
		qsort(idxs, *count_items, sizeof(indexer), cmp);
		indexer j = 1;
		indexer offset = 0;
		for (indexer i = 0; i < *count_items - 1 - offset; ++i) {
			if (idxs[i] == idxs[i + 1 + offset]) {
				offset++;
				idxs[i + 1] = idxs[i + 1 + offset];
			}
			if (offset)
				idxs[i + 1] = idxs[i + 1 + offset];
			j++;
		}
		*count_items = j;
		idxs = (indexer*)_aligned_realloc(idxs, sizeof(indexer) * j, 16);
	}
#endif
	return idxs;
}

/* searchin items in selected rectangle on cuda device imlementation (step 1) */
__global__ void cuda_search_rect2_impl1(void **nd_ptr, indexer *iter_count, indexer *atomic_iter, /*out*/ void** next_nd, /*out*/ indexer *idxs, /*out*/indexer *count_items, ret_callback2_circle callback, void *data)
{
	int idxx = threadIdx.x + blockIdx.x * blockDim.x;
	//printf("Thread = %u\n", idxx);
	//int idxx_t = idxx % (*dev_threads_count);
	indexer iter_count_t = *iter_count % (*dev_threads_count);
	if (!iter_count_t)
		iter_count_t = *iter_count;

	//if (!idxx)
	//	printf("Calc NODE iter_count = %u (0x%llx)\n", *iter_count, iter_count);

	struct node** nd = (struct node**)nd_ptr;
	//indexer idx = 0;

#ifdef CALC_POINT
	//coord tmp_dist = FLT_MAX;
	//indexer tmp_idx = -1;
#endif // CALC_POINT
	int t = (int)ceilf((float)*iter_count / (float)*dev_threads_count);
	int t1 = t - 1;
	for (int j = 0; j < t; ++j) {
		//printf("Thread = %i, j = %u, t1 = %u, idxx = %u, >= iter_count_t = %u (%s => %s)\n", idxx, j, t1, idxx, iter_count_t, idxx >= iter_count_t ? "true" : "false", j == t1 && idxx >= iter_count_t ? "true" : "false");
		if (j == t1 && idxx >= iter_count_t) {// idxx_t >= *iter_count) {
			return;
		}
		//printf("Thread %i (%i): x1 = %f, x2 = %f, y1 = %f, y2 = %f\n", idxx + j * (*dev_threads_count), j, nd[idxx + j * (*dev_threads_count)]->x1, nd[idxx + j * (*dev_threads_count)]->x2, nd[idxx + j * (*dev_threads_count)]->y1, nd[idxx + j * (*dev_threads_count)]->y2);
		// node in bounrary or bounrary in node
		if (nd[idxx + j * (*dev_threads_count)]->x1 <= dev_bonds->x_max && nd[idxx + j * (*dev_threads_count)]->x2 >= dev_bonds->x_min && nd[idxx + j * (*dev_threads_count)]->y1 <= dev_bonds->y_max && nd[idxx + j * (*dev_threads_count)]->y2 >= dev_bonds->y_min) {
			//printf("Thread %i (%i) ================================ ====\n", idxx + j * (*dev_threads_count), j);
			// check node fully in the boundary
			//if (nd[idxx].x1 >= dev_bonds->x_min && nd[idxx].y1 >= dev_bonds->y_min && nd[idxx].x2 <= dev_bonds->x_max && nd[idxx].y2 <= dev_bonds->y_max) {
				// node is fully in the boundary
			//}
			//else {
				// node isn't fully in the boundary, than add to calculation to next iteration
				indexer t1 = atomicAdd(atomic_iter, nd[idxx + j * (*dev_threads_count)]->count_child_nodes);
				//printf("Increase %i: %u to %u (%u)\n", idxx, t1, *atomic_iter, nd[idxx]->count_child_nodes);
				/*if (t1 + nd->count_child_nodes >= 10000)
					return; */
				for (unsigned k = t1, t2 = 0; k < t1 + nd[idxx + j * (*dev_threads_count)]->count_child_nodes; ++k, ++t2) {
					//void **ptr = &next_nd;
					next_nd[k] = nd[idxx + j * (*dev_threads_count)]->child_node[t2];
					//printf("Next index = %u\n", (struct branch*)(nd[idxx + j * (*dev_threads_count)]->child_node[t2]) - m_ttt_cuda_first_branch);
				}
			//}
		}
		else {
			// node and boundary isn't intersection
		}
	}
}

/* searchin items in selected rectangle on cuda device imlementation (step 1) */
__global__ void cuda_search_rect2_impl2(void **br_ptr, int *iter_count, indexer *atomic_iter, /*out*/ void** next_nd, /*out*/ indexer *idxs, /*out*/indexer *count_items)
{
	int idxx = threadIdx.x;
	int idx_gr_br = blockIdx.x;
	//printf("Thread = %u\n", idxx);
	//if (!idxx /*&& !idx_gr_br*/)
	//	printf("============================================== block = %d, %u\n", idx_gr_br, *atomic_iter);

	// for store temporary results
	__shared__ indexer temp_res[33]; // must be as blockDim.x size + 1 (for rpevious result)
	//__shared__ indexer temp_res2[32];
	__shared__ char temp_res_flag[32];
	//__shared__ int atom_index[1];
	temp_res[idxx] = (indexer)-1;
	//temp_res2[idxx] = (indexer)-1;
	if (!idxx)
		temp_res[32] = (indexer)-1;
	temp_res_flag[idxx] = -1;
	//__shared__ int index = 0;

	//if (!idxx)
	//	printf("Calc BR iter_count = %u (0x%llx), grid.x = %i\n", *iter_count, iter_count, idx_gr_br);

	//struct branch* br = (struct branch*)br_ptr[0];
	struct branch** br = (struct branch**)br_ptr;
	//if (!idxx)
	//printf("Thread branch %i: x1 = %f, x2 = %f, y1 = %f, y2 = %f (u)\n", idx_gr_br, br[idx_gr_br]->x_min, br[idx_gr_br]->x_max, br[idx_gr_br]->y_min, br[idx_gr_br]->y_max); // , br[idx_gr_br] - m_ttt_cuda_first_branch);

	if (br[idx_gr_br]->x_min <= dev_bonds->x_max && br[idx_gr_br]->x_max >= dev_bonds->x_min && br[idx_gr_br]->y_min <= dev_bonds->y_max && br[idx_gr_br]->y_max >= dev_bonds->y_min) {
		//if (!idxx)
		//	printf("------------------ %i, %u\n", idx_gr_br, br[idx_gr_br]->count_shapes);
		//int idxx_t = idxx % blockDim.x;
		int t = (int)ceilf((float)br[idx_gr_br]->count_leafs / (float)blockDim.x);
		int t1 = t - 1;
		for (int j = 0; j < t; ++j) {
			int curr_offset = j * blockDim.x;
			if (j == t1 && idxx + curr_offset >= br[idx_gr_br]->count_leafs) {
				break;
			}

			// check points to enter to boundary
			if (br[idx_gr_br]->leaf_x[idxx + curr_offset] >= dev_bonds->x_min && br[idx_gr_br]->leaf_x[idxx + curr_offset] <= dev_bonds->x_max && br[idx_gr_br]->leaf_y[idxx + curr_offset] >= dev_bonds->y_min && br[idx_gr_br]->leaf_y[idxx + curr_offset] <= dev_bonds->y_max) {
				temp_res[idxx] = br[idx_gr_br]->leaf_number[idxx + curr_offset];
				/*int t2 = atomicAdd(atomic_iter, 1);
				if (t2 >= MAX_RESULTS - 1) {
					// can not store result
					atomicSub(atomic_iter, 1);
				}
				else {
					// can store result
					//idxs[t2] = 
				}*/
			}
			//if (!idxx)
			//	atom_index[0] = 0;
			__syncthreads();

			// fill empty places
			/*if (temp_res2[idxx] != (indexer)-1) {
				int t2 = atomicAdd(atom_index, 1);
				temp_res[t2] = temp_res2[idxx];
				printf("Temp (%i): %i => %i, res = %u (%u)\n", idxx, t2, atom_index[0], temp_res[t2], temp_res2[idxx]);
				if (t2 >= 32)
					printf("Vai vai vai %i => %i\n", t2, atom_index[0]);
			}
			//__threadfence();
			__syncthreads();
			if (!idxx)
			printf("1 === %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i\n",
				temp_res[0], temp_res[1], temp_res[2], temp_res[3], temp_res[4], temp_res[5], temp_res[6], temp_res[7], temp_res[8], temp_res[9],
				temp_res[10], temp_res[11], temp_res[12], temp_res[13], temp_res[14], temp_res[15], temp_res[16], temp_res[17], temp_res[18], temp_res[19],
				temp_res[20], temp_res[21], temp_res[22], temp_res[23], temp_res[24], temp_res[25], temp_res[26], temp_res[27], temp_res[28], temp_res[29],
				temp_res[30], temp_res[31]);*/

			// packing temporary results
			if (temp_res[idxx] == temp_res[idxx + 1]) {
				__threadfence();
				temp_res[idxx + 1] = -1;
				//__threadfence();
				//if (temp_res[idxx] != -1)
				//	printf("Index of equals items = %u + %u => %u\n", idxx, idxx + 1, temp_res[idxx + 1]);
			}
			else {
				//__threadfence();
			}
			__syncthreads();
			/*if (!idxx)
			printf("2 === %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i\n",
				temp_res[0], temp_res[1], temp_res[2], temp_res[3], temp_res[4], temp_res[5], temp_res[6], temp_res[7], temp_res[8], temp_res[9],
				temp_res[10], temp_res[11], temp_res[12], temp_res[13], temp_res[14], temp_res[15], temp_res[16], temp_res[17], temp_res[18], temp_res[19],
				temp_res[20], temp_res[21], temp_res[22], temp_res[23], temp_res[24], temp_res[25], temp_res[26], temp_res[27], temp_res[28], temp_res[29],
				temp_res[30], temp_res[31]);*/
			// ckeck for previous result
			if (!idxx) {
				//if (temp_res[0] != (unsigned)-1) printf("%i:%i Prev step = %u (%u)\n", idx_gr_br, j, temp_res[32], temp_res[0]);
				if (temp_res[32] == temp_res[0]) {
					//if (temp_res[0] != (unsigned)-1) printf("From prev step = %u\n", temp_res[0]);
					temp_res[0] = -1;
				}
			}
			//__syncthreads();

			// store temporary results to global array2
			if (temp_res[idxx] != -1) {
				int t2 = atomicAdd(atomic_iter, 1);
				if (t2 >= MAX_RESULTS - 1) {
					// can not store result
					atomicSub(atomic_iter, 1);
				}
				else {
					// can store result (idxs2 - temporary)
					idxs[t2] = temp_res[idxx];
					//printf("Number = %u, value = %u (%u)\n", t2, temp_res[idxx], idxs[t2]);
					temp_res_flag[idxx] = idxx;
					//printf("Increase %i: %u to %u\n", idxx, t2, *iter_count);
				}
			}

			__syncthreads();

			// store previous result
			for (int t2 = blockDim.x / 2; t2 > 0; t2 >>= 1)
			{
				if (idxx < t2) {
					if (temp_res_flag[idxx] < temp_res_flag[idxx + t2])
						temp_res_flag[idxx] = temp_res_flag[idxx + t2];
				}
				__syncthreads();
			}
			if (!idxx) {
				if (temp_res_flag[idxx] != -1) {
					//printf("INDEX === %u\n", temp_res_flag[idxx]);
					temp_res[32] = temp_res[temp_res_flag[idxx]]; // idxs[idxx];
				}
				//if (temp_res[32] != (unsigned)-1) printf("%i:%i Next step = %u\n", idx_gr_br, j, temp_res[32]);
			}
			
			// reset temporary resulats
			temp_res[idxx] = -1;
			//temp_res2[idxx] = -1;
			temp_res_flag[idxx] = -1;
			__syncthreads();
		}
	}
}

__global__ void cuda_search_rect2_impl3(indexer *unpack_idxs, indexer *iter_count, indexer *atomic_iter, /*out*/ void **next_nd, /*out*/ indexer *idxs, /*out*/indexer *count_items)
{
	int idxx = threadIdx.x + blockIdx.x * blockDim.x;

	// for store temporary results
	__shared__ indexer temp_res[33]; // must be as blockDim.x size + 1 (for rpevious result)
	__shared__ char temp_res_flag[32];
	temp_res[idxx] = (indexer)-1;
	if (!idxx)
		temp_res[32] = (indexer)-1;
	temp_res_flag[idxx] = -1;

	int t = (int)ceilf((float)*iter_count / (float)blockDim.x);
	int t1 = t - 1;
	temp_res[idxx] = -1;
	temp_res_flag[idxx] = -1;
	if (!idxx) {
		temp_res[32] = -1;
		printf("---------------------- Count temp results = %u (t = %i) ------------------\n", *iter_count, t);
	}

	for (int j = 0; j < t; ++j) {
		int curr_offset = j * blockDim.x;
		if (j == t1 && idxx + curr_offset >= *iter_count) {
			break;
		}
		temp_res[idxx] = unpack_idxs[idxx + curr_offset];

		if (temp_res[idxx] == temp_res[idxx + 1]) {
			__threadfence();
			temp_res[idxx + 1] = -1;
		}
		__syncthreads();

		if (!idxx) {
			if (temp_res[32] == temp_res[0]) {
				temp_res[0] = -1;
			}
		}

		//__threadfence();
		if (temp_res[idxx] == 2501949) {
			printf("1 =========================================================== 250 ======================================== %i, j = %i\n", idxx, j);
			printf("1 === %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i\n",
				temp_res[0], temp_res[1], temp_res[2], temp_res[3], temp_res[4], temp_res[5], temp_res[6], temp_res[7], temp_res[8], temp_res[9],
				temp_res[10], temp_res[11], temp_res[12], temp_res[13], temp_res[14], temp_res[15], temp_res[16], temp_res[17], temp_res[18], temp_res[19],
				temp_res[20], temp_res[21], temp_res[22], temp_res[23], temp_res[24], temp_res[25], temp_res[26], temp_res[27], temp_res[28], temp_res[29],
				temp_res[30], temp_res[31]);
		}

		// store results to global array
		if (temp_res[idxx] != -1) {
			int t2 = atomicAdd(atomic_iter, 1);
			if (t2 >= MAX_RESULTS - 1) {
				// can not store result
				atomicSub(atomic_iter, 1);
			}
			else {
				// can store result (idxs2 - temporary)
				idxs[t2] = temp_res[idxx];
				//printf("Number = %u, value = %u (%u)\n", t2, temp_res[idxx], idxs[t2]);
				temp_res_flag[idxx] = idxx;
				//printf("Increase %i: %u to %u\n", idxx, t2, *iter_count);
			}
		}
		// store previous result
		for (int t2 = blockDim.x / 2; t2 > 0; t2 >>= 1)
		{
			if (idxx < t2) {
				if (temp_res_flag[idxx] < temp_res_flag[idxx + t2])
					temp_res_flag[idxx] = temp_res_flag[idxx + t2];
			}
			__syncthreads();
		}
		if (!idxx) {
			if (temp_res_flag[idxx] != -1) {
				//printf("INDEX === %u\n", temp_res_flag[idxx]);
				temp_res[32] = temp_res[temp_res_flag[idxx]]; // idxs[idxx];
			}
			//if (temp_res[32] != (unsigned)-1) printf("%i:%i Next step = %u\n", idx_gr_br, j, temp_res[32]);
		}
	}
}
