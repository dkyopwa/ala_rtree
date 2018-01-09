/*

*/

#include <stdio.h>
#include <cuda_runtime_api.h>
//#include <cuda_runtime.h>
#include "unimem.h"
#include "first.h"

#define DEBUG_CUDA_INFO

struct node* m_dev_node = NULL;

extern "C"
bool init_cuda_device(int deviceID, struct node* nd)
{
	if (!nd)
		return false;
	//return false;

	//struct node *tnd = nd;
	unsigned count1[64], i = 0, count_br = 0;
	for (int j = 0; j < 64; j++) count1[j] = 0;
	//count1[0] = 1;

	struct node *stack_node[64];
	int stack_pos = 0;
	indexer stack_idx[64];
	while (i < nd->count_child_nodes) {
		// node in bounrary or bounrary in node
				// node not fully in the boundaty
		if (nd->is_last_node) {
			/*for (unsigned j = 0; j < nd->count_child_nodes; ++j) {
				struct branch *br = (struct branch*)(nd->child_node)[j];
			}*/
			/*if (!count_br)
				count_br = nd->count_child_nodes;
			if (count_br != nd->count_child_nodes)
				printf("Branches %u vs %u\n", count_br, nd->count_child_nodes);*/
			count_br += nd->count_child_nodes;
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
	for (int i = 0; i < deviceCount; ++i) {
		er1 = cudaGetDeviceProperties(&prop, i);
		if (prop.major < 2)
		{
			printf("ERROR: calculation requires GPU devices with compute SM 2.0 or higher.\n");
			printf("Current GPU device has compute SM%d.%d, Exiting...", prop.major, prop.minor);
			//exit(EXIT_WAIVED);
			return false;
		}

		printf("GPU device name is %s\n", prop.name);
		printf("GPU total memory = %.0f Mb\n", prop.totalGlobalMem / 1024.0 / 1024.0);
		printf("Number of multiprocessors on the device = %u\n", prop.multiProcessorCount);
	}

	er1 = cudaSetDevice(deviceID);

	// copy rtree
	int pos = 63;

	for (; pos >= 0; --pos) {
		if (count1[pos])
			break;
	}

	// allocationg memory for branches
	alignas(16) struct branch* tbr = (struct branch*)aligned_alloc(16, sizeof(struct branch) * count_br);
	// struct branch* first_branch = NULL;
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
	}
	struct branch *dev_br = NULL;
	er1 = cudaMalloc((void**)&dev_br, sizeof(struct branch) * count_br);
	er1 = cudaMemcpy(dev_br, tbr, sizeof(struct branch) * count_br, cudaMemcpyHostToDevice);

	//return false;
	alignas(16) struct node *to_dev_nd[64];
	void **to_dev_child[64];
	struct node *dev_nd = NULL, *dev_nd_prev = NULL, *dev_ptr = NULL;
	// to_dev_nd[0] = (struct node*)aligned_alloc(16, sizeof(struct node));
	// memcpy(to_dev_nd[0], nd, sizeof(struct node));
	struct node* tnd = nd;
	for (unsigned j = 0; j <= pos; ++j)
		tnd = (struct node*)(tnd->child_node[0]);
	//unsigned j = 0;
	//void* tmp1 = NULL;
	unsigned count = tnd->count_child_nodes, prev_count = 1;
	for (int k1 = pos; k1 >= 0; --k1) {
		// data child node
		to_dev_nd[k1] = (struct node*)aligned_alloc(16, sizeof(struct node) * count1[k1]);
		memcpy(to_dev_nd[k1], tnd/*->child_node[0]*/, sizeof(struct node) * count1[k1]);
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
					struct branch *ptr = &(dev_br[k2 + j * MAX_NODES]);
					er1 = cudaMemcpy((void*)((to_dev_nd[k1])[j].child_node + k2), &ptr, sizeof(struct branch*), cudaMemcpyHostToDevice);
				}
				else {
					// copy pointer of nodes
					struct node* ptr = &(dev_nd_prev[k2 + j * MAX_NODES]);
					er1 = cudaMemcpy((void*)((to_dev_nd[k1])[j].child_node + k2), &ptr, sizeof(struct node*), cudaMemcpyHostToDevice);
				}
			}
		}
		// pointers of child nodes
		er1 = cudaMalloc((void**)&dev_nd, sizeof(struct node*) * count1[k1]); // tnd->count_child_nodes);
		cudaMemcpy(dev_nd, to_dev_nd[k1], sizeof(struct node) * count1[k1], cudaMemcpyHostToDevice);
		dev_nd_prev = dev_nd;
		continue;

		indexer k3 = 0;
		for (unsigned k2 = 0; k2 < count1[k1]; ++k2) {
			for (indexer j = 0; j < MAX_NODES; ++j) {
				((to_dev_nd[k1])[k2]).child_node[j] = &(dev_nd[k3]);
				k3++;
			}
		}


		break;
		// copy to device
		//cudaMemcpy(dev_nd, to_dev_nd[j]->child_node[0], sizeof(struct node*) * count, cudaMemcpyHostToDevice);
		// pointer to child_node on device
	/*	er1 = cudaMalloc((void**)&dev_nd, sizeof(void*) * count);
		to_dev_nd[j]->child_node = (void**)dev_nd;
		// copy to device child_node
		//cudaMemcpy(dev_nd, to_dev_nd[j]->child_node, sizeof(void*) * count, cudaMemcpyHostToDevice);

		// prepare device to copy nodes
		//cudaMalloc((void**)&dev_nd, sizeof(struct node) * prev_count);
		// copy to device nodes
		//cudaMemcpy(dev_nd, tnd, sizeof(struct node) * prev_count, cudaMemcpyHostToDevice);
		if (j == 1) {
			//to_dev_nd[0]->
		}

		prev_count = count;

		// free local memory
		//_aligned_free(tmp1); // to_dev_nd[j]->child_node);
		//_aligned_free(to_dev_nd[j]);

		j++;
		count = 0;
		if (!((struct node*)tnd->child_node[0])->is_last_node)
			for (unsigned k = 0; k < tnd->count_child_nodes; ++k) {
				count += ((struct node*)(tnd->child_node[k]))->count_child_nodes;
			}
		tnd = (struct node*)(tnd->child_node[0]);*/
	}

	// allocating memory for root
	//er1 = cudaMalloc((void**)&m_dev_node, sizeof(struct node));
	// copy to device root of tree
	//cudaMemcpy(m_dev_node, to_dev_nd[0], sizeof(struct node), cudaMemcpyHostToDevice);

	return true;
}

extern "C"
bool destroy_cuda_device()
{
	cudaError_t er1 = cudaDeviceReset();
	return er1 == cudaSuccess ? true : false;
}