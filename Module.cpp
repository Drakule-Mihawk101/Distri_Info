/*
 *	Author:	Seung-Hee Bae (shbae@cs.washington.edu)
 *	Date:	Mar. 2014
 *	Copyright (C) since 2013,  Seung-Hee Bae, Bill Howe, Database Group at the University of Washington
 */

#include <iostream>
#include <sstream>
#include <fstream>
#include <omp.h>
#include <cstdlib>
#if defined __GNUCC__ || defined __APPLE__
#include <ext/hash_map>
#else
#include <hash_map>
#endif
#include <ctime>
#include <sys/time.h>
#include "Module.h"
#include "Node.h"
#include "timing.h"
#include <mpi.h>
#include <cstdio>
#include <algorithm>
#include <iterator>

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"

//typedef __gnu_cxx::hash_map<int, double> flowmap;
typedef map<int, double> flowmap;
typedef map<int, pair<double, double> > modInfo;	// <modID, <exitPr, sumPr> >

void findAssignedPart(int* start, int* end, int numNodes, int numTh, int myID);

const double InitLength = 10000.0;

using namespace std;

struct MoveSummary {
	double diffCodeLen;			// delta(L(M))
	int newModule;
	double sumPr1, sumPr2;// updated sumPr1 and sumPr2.  1 --> oldM, 2 --> newM.
	double exitPr1, exitPr2;	// updated exitPr1 and exitPr2.
	double newSumExitPr;// SUM(q_i) = this.sumAllExitPr + q_1_new + q_2_new - q_1_old - q_2_old.
};

// Constructors of Module class.
Module::Module() :
		index(0), exitPr(0.0), stayPr(0.0), sumPr(0.0), sumTPWeight(0.0), sumDangling(
				0.0), numMembers(1) {
}

Module::Module(int idx, double exitPr, double sumPr) :
		index(idx), exitPr(exitPr), stayPr(exitPr + sumPr), sumPr(sumPr), sumTPWeight(
				0.0), sumDangling(0.0), numMembers(1) {
}

//Module::Module(int idx, int nNode, Node nd)
Module::Module(int idx, Node * nd) :
		index(idx), exitPr(0.0), stayPr(0.0), sumPr(nd->Size()), sumTPWeight(
				nd->TeleportWeight()), sumDangling(0.0), numMembers(1) {
	double sumExitFlow = 0.0;

	int nOutLinks = nd->outLinks.size();

	for (int i = 0; i < nOutLinks; i++)
		sumExitFlow += nd->outLinks[i].second;// after call initiate(), w_ab is updated as p_a * w_ab.

	if (nd->IsSuper()) {// If SuperNode, we already calculated exitPr in the corresponding module.
		exitPr = nd->ExitPr();
		sumDangling = nd->DanglingSize();
	}
	// exitPr = tau * (1 - sum(tau_a))*sum(p_a) + (1-tau) * sumExitFlow.
	else if (!nd->IsDangling()) {
		exitPr = Network::alpha * (1.0 - nd->TeleportWeight()) * sumPr
				+ Network::beta * sumExitFlow;
		nd->setExitPr(exitPr);
	} else {
		exitPr = (1.0 - nd->TeleportWeight()) * sumPr;
		nd->setExitPr(exitPr);
		sumDangling = nd->Size();
	}

	stayPr = exitPr + sumPr;	// exitPr + sum (P_a).
	members.push_back(nd);
}

SubModule::SubModule() :
		modIdx(0), numMembers(0), sumPr(0.0), sumTPWeight(0.0), sumDangling(0.0) {
}

/**
 *	Module mod : module of newNetwork (network based on a module in original network.)
 */
SubModule::SubModule(Module& mod, map<int, int>& origNodeID, int modIndex) {
	numMembers = mod.numMembers;
	sumPr = mod.sumPr;
	sumTPWeight = mod.sumTPWeight;
	sumDangling = mod.sumDangling;
	modIdx = modIndex;

	for (vector<Node*>::iterator it = mod.members.begin();
			it != mod.members.end(); it++) {
		members.push_back(origNodeID[(*it)->ID()]);
	}

	// numMembers should be equal to members.size() ...
	if (numMembers != members.size())
		cout
				<< "SOMETHING WRONG!!! -- numMembers != members.size() in SubModule()..."
				<< endl;
}

SubModule::SubModule(Module& mod) {
	modIdx = mod.index;
	numMembers = 1;
	sumPr = mod.sumPr;
	sumTPWeight = mod.sumTPWeight;
	sumDangling = mod.sumDangling;

	members.push_back(mod.members[0]->ID());
}

// Constructors of Network class.
Network::Network() :
		level(0), codeLength(InitLength), nNode(0), nEdge(0), nEmptyMod(0), nModule(
				0), totNodeWeights(0.0), nDanglings(0), allNodes_log_allNodes(
				0.0), sumAllExitPr(0.0) {
}

Network::Network(int numNode, int numEdge, int level) :
		level(level), codeLength(InitLength), nNode(numNode), nEdge(numEdge), nEmptyMod(
				0), nModule(numNode), totNodeWeights(0.0), nDanglings(0), allNodes_log_allNodes(
				0.0), sumAllExitPr(0.0) {
	emptyModules.reserve(numNode);
}

Network::Network(int numNode, int numEdge, int level, double codeLen) :
		level(level), codeLength(codeLen), nNode(numNode), nEdge(numEdge), nEmptyMod(
				0), nModule(numNode), totNodeWeights(0.0), nDanglings(0), allNodes_log_allNodes(
				0.0), sumAllExitPr(0.0) {
	emptyModules.reserve(numNode);
}

/*
 * Other Member functions of Network class.
 */

void Network::findDanglingNodes() {

	nDanglings = 0;		// reset the number of dangling nodes.

	for (int i = 0; i < nNode; i++) {
		if (nodes[i].outLinks.empty()) {
			danglings.push_back(i);
			nDanglings++;
			nodes[i].setIsDangling(true);
		}
	}

}

//void Network::initiate() {
void Network::initiate(int numTh) {

	int size, rank;

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int nDangNodes = 0;

	struct timeval startT, endT;

	omp_set_num_threads(numTh);

	gettimeofday(&startT, NULL);

#pragma omp parallel reduction (+:nDangNodes)
	{
		int myID = omp_get_thread_num();
		int nTh = omp_get_num_threads();

		int start, end;
		findAssignedPart(&start, &end, nNode, nTh, myID);

		for (int i = start; i < end; i++) {
			if (nodes[i].outLinks.empty()) {
#pragma omp critical (updateDang)
				{
					danglings.push_back(i);
				}
				nDangNodes++;
				nodes[i].setIsDangling(true);
			} else {	// normal nodes. --> Normalize edge weights.
				int nOutLinks = nodes[i].outLinks.size();
				//double sum = nodes[i].selfLink;	// don't support selfLink yet.
				double sum = 0.0;
				for (int j = 0; j < nOutLinks; j++)
					sum += nodes[i].outLinks[j].second;

				//nodes[i].selfLink /= sum;
				for (int j = 0; j < nOutLinks; j++)
					nodes[i].outLinks[j].second /= sum;
			}
		}
	}

	gettimeofday(&endT, NULL);

	nDanglings = nDangNodes;

	cout << "Level " << level << ": the number of dangling nodes = "
			<< nDanglings << endl;
	cout << "Time for finding dangling nodes : "
			<< elapsedTimeInSec(startT, endT) << " (sec)" << endl;

	gettimeofday(&startT, NULL);
	calculateSteadyState(numTh);
	gettimeofday(&endT, NULL);
	cout << "Time for calculating steady state of nodes (eigenvector): "
			<< elapsedTimeInSec(startT, endT) << " (sec)" << endl;

	gettimeofday(&startT, NULL);

	// Update edges to represent flow based on calculated steady state (aka size).
#pragma omp parallel
	{
		int myID = omp_get_thread_num();
		int nTh = omp_get_num_threads();

		int start, end;
		findAssignedPart(&start, &end, nNode, nTh, myID);

		for (int i = start; i < end; i++) {

			if (!nodes[i].IsDangling()) {
				int nOutLinks = nodes[i].outLinks.size();
				for (int j = 0; j < nOutLinks; j++) {
					nodes[i].outLinks[j].second = nodes[i].Size()
							* nodes[i].outLinks[j].second;
				}
			} else {
				nodes[i].setDanglingSize(nodes[i].Size());
			}
		}
	}

	gettimeofday(&endT, NULL);
	cout << "Time for updating edge weights based on flow information: "
			<< elapsedTimeInSec(startT, endT) << " (sec)" << endl;

	// Calculate SUM of p_log_p over all nodes.
	double allNds_log_allNds = 0.0;
	double allNds_log_allNds_s = 0.0;

	/*#pragma omp parallel reduction (+:allNds_log_allNds)
	 {
	 int myID = omp_get_thread_num();
	 int nTh = omp_get_num_threads();

	 int start, end;
	 findAssignedPart(&start, &end, nNode, nTh, myID);

	 for (int i = start; i < end; i++)
	 allNds_log_allNds += pLogP(nodes[i].Size());
	 }*/

	int start, end;
	findAssignedPart(&start, &end, nNode, size, rank);

	for (int i = start; i < end; i++) {
		allNds_log_allNds_s += pLogP(nodes[i].Size());
	}
	if (rank != 0) {
		MPI_Send(&allNds_log_allNds_s, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}

	if (rank == 0) {
		allNds_log_allNds = allNds_log_allNds_s;
		for (int prId = 1; prId < size; prId++) {
			MPI_Recv(&allNds_log_allNds_s, 1, MPI_DOUBLE, prId, 0,
			MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
			allNds_log_allNds += allNds_log_allNds_s;
		}
		for (int prId = 1; prId < size; prId++) {
			MPI_Send(&allNds_log_allNds, 1, MPI_DOUBLE, prId, 0,
			MPI_COMM_WORLD);

		}
	}

	if (rank != 0) {
		MPI_Recv(&allNds_log_allNds, 1, MPI_DOUBLE, 0, 0,
		MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
	}

	allNodes_log_allNodes = allNds_log_allNds;

	/////////////////////////////
	// Make modules from nodes //
	/////////////////////////////
	/*#pragma omp parallel
	 {*/
	int myID = omp_get_thread_num();
	int nTh = omp_get_num_threads();

	findAssignedPart(&start, &end, nNode, nTh, myID);

	for (int i = start; i < end; i++) {
		modules[i] = Module(i, &nodes[i]);// Assign each Node to a corresponding Module object.
										  // Initially, each node is in its own module.
		nodes[i].setModIdx(i);
	}
	//}

	nModule = nNode;

	gettimeofday(&startT, NULL);
	calibrate(numTh);
	gettimeofday(&endT, NULL);
	cout << "Time for calculating of initial code length: "
			<< elapsedTimeInSec(startT, endT) << " (sec)" << endl;

}

// calculating steady state of nodes by Power Iteration Method.
// same as Greedy::eigenvector() in Infomap implementation.
// Modify for the design of this implementation.

void Network::calculateSteadyState(int numTh) {
	// initial probability distribution = 1 / N.

	MPI_Request request;
	//cout << "inside calculate steady state" << endl;
	int size;
	int rank;
	int iteration = 0;

	vector<double> size_tmp = vector<double>(nNode, 1.0 / nNode);

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	//cout << "sizeof dangling:" << danglings.size() << endl;

	/*for (int i = 0; i < danglings.size(); i++)
	 cout << "rank: " << rank << " danglings:" << danglings[i] << endl;*/

	//cout << "inside calculate steady state 1" << endl;
	/*struct data {
	 int start;
	 int end;
	 double danglingsz;
	 };

	 data chunkobject;
	 MPI_Datatype chunktype;
	 int structlen = 3;
	 int blocklengths[structlen];
	 MPI_Datatype types[structlen];
	 MPI_Aint displacements[structlen];
	 // where are the components relative to the structure?
	 blocklengths[0] = 1;
	 types[0] = MPI_INT;
	 displacements[0] = (size_t) &(chunkobject.start) - (size_t) &chunkobject;
	 blocklengths[1] = 1;
	 types[1] = MPI_INT;
	 displacements[1] = (size_t) &(chunkobject.end) - (size_t) &chunkobject;
	 blocklengths[2] = 1;
	 types[2] = MPI_DOUBLE;
	 displacements[2] = (size_t) &(chunkobject.danglingsz)
	 - (size_t) &chunkobject;

	 MPI_Type_create_struct(structlen, blocklengths, displacements, types,
	 &chunktype);
	 MPI_Type_commit(&chunktype);
	 {
	 MPI_Aint typesize;
	 MPI_Type_extent(chunktype, &typesize);
	 if (rank == 0)
	 printf("Type extent: %d bytes\n", typesize);
	 }*/

	int iter = 0;
	double danglingSize = 0.0;
	double sqdiff = 1.0;
	double sum = 0.0;
	double sum_s = 0.0;
	double danglingsz = 0.0;
	int value = 0;
	// Generate addedSize array per each thread, so that we don't need to use lock for each nodeSize addition.

	/*
	 double** addedSize = new double*[numTh];
	 for (int i = 0; i < numTh; i++)
	 addedSize[i] = new double[nNode];
	 */

	do {

		// calculate sum of the size of dangling nodes.
		iteration++;
		danglingSize = 0.0;
		danglingsz = 0.0;
		int start, end, id;
		findAssignedPart(&start, &end, nDanglings, size, rank);
		// assign dangling array computation to itself rank=0 process
		for (int i = start; i < end; i++) {
			/*cout << "i:" << i << " danglings[i]: " << danglings[i] << " rank"
			 << rank << " size_tmp[danglings[i]]: "
			 << size_tmp[danglings[i]] << " iteration number: "
			 << iteration << endl;*/
			danglingsz += size_tmp[danglings[i]];
		}
		if (rank != 0) {
			MPI_Send(&danglingsz, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		}
		if (rank == 0) {
			danglingSize += danglingsz;
			for (int prId = 1; prId < size; prId++) {
				MPI_Recv(&danglingsz, 1, MPI_DOUBLE, prId, 0, MPI_COMM_WORLD,
				MPI_STATUSES_IGNORE);
				danglingSize += danglingsz;
			}
			//cout << "final summed up dangling value:" << danglingSize << endl;
			//we need to send the value of danglingSize to all of the process back again so that they have the updated value
			for (int prId = 1; prId < size; prId++) {
				/*cout << "before sending from:" << prId << " value of dangling:"
				 << danglingSize << endl;*/
				MPI_Send(&danglingSize, 1, MPI_DOUBLE, prId, 0, MPI_COMM_WORLD);
			}
		}

		if (rank != 0) {
			danglingSize = 0.0;
			MPI_Recv(&danglingSize, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,
			MPI_STATUSES_IGNORE);
			/*cout << "danglingSize value: " << danglingSize << " rank:" << rank
			 << endl;*/
		}

		/*		// calculate sum of the size of dangling nodes.
		 #pragma omp parallel reduction (+:danglingSize)
		 {
		 int myID = omp_get_thread_num();
		 int nTh = omp_get_num_threads();

		 int start, end;
		 findAssignedPart(&start, &end, nDanglings, nTh, myID);

		 for (int i = start; i < end; i++) {
		 cout << "i:" << i << " danglings[i]: " << danglings[i]
		 << " rank" << rank << " size_tmp[danglings[i]]: "
		 << size_tmp[danglings[i]] << " iteration number: "
		 << iteration << endl;
		 danglingSize += size_tmp[danglings[i]];
		 }
		 }
		 cout << "final summed up dangling value:" << danglingSize << endl;*/

		// flow via teleportation.
		// faysal: I think in the following block I do not need to call findAssignedPart again
		//faysal: I need to send update of the node vector back to all of the process, but I did not do that, what I did, I just made all the process to update the nodes vector
//#pragma omp parallel
		{
			//int myID = omp_get_thread_num();
			//int nTh = omp_get_num_threads();

			//int start, end;
			//findAssignedPart(&start, &end, nNode, nTh, myID);

			for (int i = 0; i < nodes.size(); i++) {
				nodes[i].setSize(
						(alpha + beta * danglingSize)
								* nodes[i].TeleportWeight());//alpha is 0.15, beta is 1-0.15 or 0.85, teleportation weight is individual nodeweight/totalNodeweight,
								//size is p_alpha, hence (0.15+0.85*danglingsize)*nodes[i].TeleportWeight()
			}
		}

		/*		for (int i = start; i < end; i++) {
		 nodes[i].setSize(
		 (alpha + beta * danglingSize) * nodes[i].TeleportWeight());	//alpha is 0.15, beta is 1-0.15 or 0.85, teleportation weight is individual nodeweight/totalNodeweight,
		 //size is p_alpha, hence (0.15+0.85*danglingsize)*nodes[i].TeleportWeight()
		 }*/

		/*int realNumTh = 0;

		 // flow from network steps via following edges.
		 #pragma omp parallel
		 {
		 int myID = omp_get_thread_num();
		 int nTh = omp_get_num_threads();

		 #pragma omp master
		 {
		 realNumTh = nTh;
		 }

		 //here instead of myID I am using rank

		 double *myAddedSize = addedSize[rank];
		 for (int i = 0; i < nNode; i++)
		 myAddedSize[i] = 0.0;// initialize for the temporary addedSize array.

		 int start, end;
		 findAssignedPart(&start, &end, nNode, nTh, myID);

		 // here instead of a part, I am putting the size of nNode here
		 for (int i = start; i < end; i++) {
		 int nOutLinks = nodes[i].outLinks.size();
		 for (int j = 0; j < nOutLinks; j++) {
		 myAddedSize[nodes[i].outLinks[j].first] += beta
		 * nodes[i].outLinks[j].second * size_tmp[i];//0.85*nodes[i].outlinks[j].weight*(1/nNodes)

		 }
		 }
		 }
		 //set nodes[i].size by added addedSize[] values.
		 #pragma omp parallel
		 {
		 int myID = omp_get_thread_num();
		 int nTh = omp_get_num_threads();

		 int start, end;
		 findAssignedPart(&start, &end, nNode, nTh, myID);

		 for (int i = start; i < end; i++) {
		 for (int j = 0; j < realNumTh; j++)
		 nodes[i].addSize(addedSize[j][i]);
		 }
		 }
		 */

		//I am using myAddedSize array where each and every process will have its' own copy
		//double myAddedSize = new double[nNode];
		for (int i = 0; i < nNode; i++) {
			//myAddedSize[i] = 0.0;	//initialize the temporary addedSize array
			int nOutLinks = nodes[i].outLinks.size();
			for (int j = 0; j < nOutLinks; j++) {
				nodes[nodes[i].outLinks[j].first].addSize(
						beta * nodes[i].outLinks[j].second * size_tmp[i]);
			}
		}

		// Normalize of node size.
		sum = 0.0;
		sum_s = 0.0;
		/*
		 #pragma omp parallel reduction (+:sum)
		 {
		 int myID = omp_get_thread_num();
		 int nTh = omp_get_num_threads();

		 int start, end;
		 findAssignedPart(&start, &end, nNode, nTh, myID);

		 for (int i = start; i < end; i++)
		 sum += nodes[i].Size();
		 }
		 */

		findAssignedPart(&start, &end, nNode, size, rank);
		for (int i = start; i < end; i++) {
			sum_s += nodes[i].Size();
		}
		if (rank != 0) {
			MPI_Send(&sum_s, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		}

		if (rank == 0) {
			sum = sum_s;
			for (int prId = 1; prId < size; prId++) {
				MPI_Recv(&sum_s, 1, MPI_DOUBLE, prId, 0, MPI_COMM_WORLD,
				MPI_STATUSES_IGNORE);
				sum += sum_s;
			}
			for (int prId = 1; prId < size; prId++) {
				MPI_Send(&sum, 1, MPI_DOUBLE, prId, 0, MPI_COMM_WORLD);
			}
		}

		if (rank != 0) {
			MPI_Recv(&sum, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,
			MPI_STATUSES_IGNORE);
		}

		/*cout << "final value of sum:" << sum << " iteration:" << iteration
		 << endl;*/
		sqdiff = 0.0;
		/*#pragma omp parallel reduction (+:sqdiff)
		 {
		 int myID = omp_get_thread_num();
		 int nTh = omp_get_num_threads();

		 int start, end;
		 findAssignedPart(&start, &end, nNode, nTh, myID);

		 for (int i = start; i < end; i++) {
		 nodes[i].setSize(nodes[i].Size() / sum);
		 sqdiff += fabs(nodes[i].Size() - size_tmp[i]);
		 size_tmp[i] = nodes[i].Size();
		 }
		 }*/

		for (int i = 0; i < nNode; i++) {
			nodes[i].setSize(nodes[i].Size() / sum);
			sqdiff += fabs(nodes[i].Size() - size_tmp[i]);
			size_tmp[i] = nodes[i].Size();
		}

		iter++;

	} while ((iter < 200) && (sqdiff > 1.0e-15 || iter < 50));

	/*	// deallocate 2D array.
	 for (int i = 0; i < numTh; i++)
	 delete[] addedSize[i];
	 delete[] addedSize;*/

	cout << "Calculating flow done in " << iter << " iterations!" << endl;

}

// This function calculate current codeLength.
// This implementation is modified version of infomap implementation.
void Network::calibrate(int numTh) {
	//This is the calculation of Equation (4) in the paper.

	/*	int size;
	 int rank;

	 MPI_Comm_size(MPI_COMM_WORLD, &size);
	 MPI_Comm_rank(MPI_COMM_WORLD, &rank);*/

	/*double *data = new double[3];
	 double *data_s = new double[3];
	 */

	int size;
	int rank;

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	double sum_exit_log_exit = 0.0;
	double sum_stay_log_stay = 0.0;
	double sumExit = 0.0;

	/*
	 double sum_exit_log_exit_s = 0.0;
	 double sum_stay_log_stay_s = 0.0;
	 double sumExit_s = 0.0;
	 */

	omp_set_num_threads(numTh);
//#pragma omp parallel reduction (+:sum_exit_log_exit, sum_stay_log_stay, sumExit)
	//{
	int myID = omp_get_thread_num();
	int nTh = omp_get_num_threads();

	int start, end;
	findAssignedPart(&start, &end, nModule, nTh, myID);

	for (int i = start; i < end; i++) {
		sum_exit_log_exit += pLogP(modules[i].exitPr);
		sum_stay_log_stay += pLogP(modules[i].stayPr);
		sumExit += modules[i].exitPr;
		//printf("the exit probability for rank:%d and module:%d is:%f\n",
		//rank, i, modules[i].exitPr);
	}
	//}

	/*cout << "ping" << endl;
	 sum_exit_log_exit = 0.0; //sum_exit_log_exit
	 sum_stay_log_stay = 0.0;	//sum_stay_log_stay
	 sumExit = 0.0;  //sumExit



	 sum_exit_log_exit_s = 0.0; //sum_exit_log_exit_s
	 sum_stay_log_stay_s = 0.0; //sum_stay_log_stay_s
	 sumExit_s = 0.0; //sumExit_s


	 int start, end;
	 findAssignedPart(&start, &end, nModule, numTh, rank);

	 for (int i = start; i < end; i++) {
	 sum_exit_log_exit_s += pLogP(modules[i].ExitPr());
	 sum_stay_log_stay_s += pLogP(modules[i].StayPr());
	 sumExit_s += modules[i].ExitPr();
	 }
	 if (rank != 0) {
	 MPI_Send(&sum_exit_log_exit_s, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	 MPI_Send(&sum_stay_log_stay_s, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	 MPI_Send(&sumExit_s, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	 }
	 cout << "pong2" << endl;

	 if (rank == 0) {
	 sum_exit_log_exit = sum_exit_log_exit_s;
	 sum_stay_log_stay = sum_stay_log_stay_s;
	 sumExit = sumExit_s;
	 for (int prId = 1; prId < size; prId++) {
	 MPI_Recv(&sum_exit_log_exit_s, 1, MPI_DOUBLE, prId, 0,
	 MPI_COMM_WORLD,
	 MPI_STATUSES_IGNORE);
	 MPI_Recv(&sum_stay_log_stay_s, 1, MPI_DOUBLE, prId, 0,
	 MPI_COMM_WORLD,
	 MPI_STATUSES_IGNORE);
	 MPI_Recv(&sumExit_s, 1, MPI_DOUBLE, prId, 0,
	 MPI_COMM_WORLD,
	 MPI_STATUSES_IGNORE);
	 sum_exit_log_exit += sum_exit_log_exit_s;
	 sum_stay_log_stay += sum_stay_log_stay_s;
	 sumExit += sumExit_s;
	 }

	 for (int prId = 1; prId < size; prId++) {
	 MPI_Send(&sum_exit_log_exit, 1, MPI_DOUBLE, prId, 0,
	 MPI_COMM_WORLD);
	 MPI_Send(&sum_stay_log_stay, 1, MPI_DOUBLE, prId, 0,
	 MPI_COMM_WORLD);
	 MPI_Send(&sumExit, 1, MPI_DOUBLE, prId, 0, MPI_COMM_WORLD);
	 }
	 }
	 cout << "pong1" << endl;

	 if (rank != 0) {
	 MPI_Recv(&sum_exit_log_exit, 1, MPI_DOUBLE, 0, 0,
	 MPI_COMM_WORLD,
	 MPI_STATUSES_IGNORE);
	 MPI_Recv(&sum_stay_log_stay, 1, MPI_DOUBLE, 0, 0,
	 MPI_COMM_WORLD,
	 MPI_STATUSES_IGNORE);
	 MPI_Recv(&sumExit, 1, MPI_DOUBLE, 0, 0,
	 MPI_COMM_WORLD,
	 MPI_STATUSES_IGNORE);
	 }

	 cout << "pong" << endl;*/
	sumAllExitPr = sumExit;
	double sumExit_log_sumExit = pLogP(sumExit);

	codeLength = sumExit_log_sumExit - 2.0 * sum_exit_log_exit
			+ sum_stay_log_stay - allNodes_log_allNodes;

	/*	sumAllExitPr = data[2];
	 double sumExit_log_sumExit = pLogP(data[2]);
	 codeLength = sumExit_log_sumExit - 2.0 * data[0] + data[1]
	 - allNodes_log_allNodes;*/

	//cout << "abcd" << endl;
	/*	delete[] data;
	 delete[] data_s;*/

}

/*
 * This function implements the 1) procedure of stochastic_greedy_partition(Network &network) function in SeqInfomap.cpp.
 *
 *	1) in random sequential order, each node is moved to its neighbor module that results in the largest gain of the map eq.
 *	   If no move results in a gain of the map equation, the node stays in its original module.
 *	
 *	return: the number of moves.
 */
int Network::move() {

// Generate random sequential order of nodes.
	vector<int> randomOrder(nNode);
	for (int i = 0; i < nNode; i++)
		randomOrder[i] = i;

	for (int i = 0; i < nNode; i++) {
		int target = R->randInt(nNode - 1);

		// swap numbers between i and target.
		int tmp = randomOrder[i];
		randomOrder[i] = randomOrder[target];
		randomOrder[target] = tmp;
	}

	int numMoved = 0;	// the counter for the number of movements.

// Move each node to one of its neighbor modules in random sequential order.
	for (int i = 0; i < nNode; i++) {
		Node& nd = nodes[randomOrder[i]];// look at i_th Node of the random sequential order.
		int oldMod = nd.ModIdx();

		int nModLinks = 0;// The number of links to/from between the current node and other modules.

		flowmap outFlowToMod;	// <modID, flow> for outFlow...
		flowmap inFlowFromMod;	// <modID, flow> for inFlow...

		// count other modules that are connected to the current node.
		// During this counting, aggregate outFlowNotToOldMod and inFlowFromOldMod values.
		for (link_iterator linkIt = nd.outLinks.begin();
				linkIt != nd.outLinks.end(); linkIt++) {
			int newMod = nodes[linkIt->first].ModIdx();

			if (outFlowToMod.count(newMod) > 0) {
				outFlowToMod[newMod] += beta * linkIt->second;
			} else {
				outFlowToMod[newMod] = beta * linkIt->second;// initialization of the outFlow of the current newMod.
				inFlowFromMod[newMod] = 0.0;
				nModLinks++;
			}
		}

		for (link_iterator linkIt = nd.inLinks.begin();
				linkIt != nd.inLinks.end(); linkIt++) {
			int newMod = nodes[linkIt->first].ModIdx();

			if (inFlowFromMod.count(newMod) > 0) {
				inFlowFromMod[newMod] += beta * linkIt->second;
			} else {
				outFlowToMod[newMod] = 0.0;
				inFlowFromMod[newMod] = beta * linkIt->second;
				nModLinks++;
			}
		}

		if (nModLinks != outFlowToMod.size())
			cout
					<< "ALERT: nModLinks != outFlowToMod.size() in Network::move()."
					<< endl;

		// copy node specific values for easy use and efficiency.
		double ndSize = nd.Size();		// p_nd.
		double ndTPWeight = nd.TeleportWeight();		// tau_nd.
		double ndDanglingSize = nd.DanglingSize();

		double oldExitPr1 = modules[oldMod].exitPr;
		double oldSumPr1 = modules[oldMod].sumPr;
		double oldSumDangling1 = modules[oldMod].sumDangling;
		double oldModTPWeight = modules[oldMod].sumTPWeight;

		double additionalTeleportOutFlow = (alpha * ndSize
				+ beta * ndDanglingSize) * (oldModTPWeight - ndTPWeight);
		double additionalTeleportInFlow = (alpha * (oldSumPr1 - ndSize)
				+ beta * (oldSumDangling1 - ndDanglingSize)) * ndTPWeight;

		// For teleportation and danling nodes.
		for (flowmap::iterator it = outFlowToMod.begin();
				it != outFlowToMod.end(); it++) {
			int newMod = it->first;
			if (newMod == oldMod) {
				outFlowToMod[newMod] += additionalTeleportOutFlow;
				inFlowFromMod[newMod] += additionalTeleportInFlow;
			} else {
				outFlowToMod[newMod] += (alpha * ndSize + beta * ndDanglingSize)
						* modules[newMod].sumTPWeight;
				inFlowFromMod[newMod] += (alpha * modules[newMod].sumPr
						+ beta * modules[newMod].sumDangling) * ndTPWeight;
			}
		}

		// Calculate flow to/from own module (default value if no links to own module).
		double outFlowToOldMod = additionalTeleportOutFlow;
		double inFlowFromOldMod = additionalTeleportInFlow;
		if (outFlowToMod.count(oldMod) > 0) {
			outFlowToOldMod = outFlowToMod[oldMod];
			inFlowFromOldMod = inFlowFromMod[oldMod];
		}

		//////////////////// THE OPTION TO MOVE TO EMPTY MODULE ////////////////
		if (modules[oldMod].members.size() > 1 && emptyModules.size() > 0) {
			int emptyTarget = emptyModules.back();
			outFlowToMod[emptyTarget] = 0.0;
			inFlowFromMod[emptyTarget] = 0.0;
			nModLinks++;
		}

		//vector<MoveSummary> moveResults(nModLinks);
		MoveSummary currentResult;
		MoveSummary bestResult;

		double newExitPr1 = oldExitPr1 - nd.ExitPr() + outFlowToOldMod
				+ inFlowFromOldMod;

		bestResult.diffCodeLen = 0.0;// This is the default value, if we can't find diffCodeLen < 0, then don't move the node.

		for (flowmap::iterator it = outFlowToMod.begin();
				it != outFlowToMod.end(); it++) {
			int newMod = it->first;
			double outFlowToNewMod = it->second;
			double inFlowFromNewMod = inFlowFromMod[newMod];

			if (newMod != oldMod) {
				// copy module specific values...
				double oldExitPr2 = modules[newMod].exitPr;
				double oldSumPr2 = modules[newMod].sumPr;

				// Calculate status of current investigated movement of the node nd.
				currentResult.newModule = newMod;
				currentResult.sumPr1 = oldSumPr1 - ndSize;// This should be 0.0, because oldModule will be empty module.
				currentResult.sumPr2 = oldSumPr2 + ndSize;
				currentResult.exitPr1 = newExitPr1;
				currentResult.exitPr2 = oldExitPr2 + nd.ExitPr()
						- outFlowToNewMod - inFlowFromNewMod;

				currentResult.newSumExitPr = sumAllExitPr + newExitPr1
						+ currentResult.exitPr2 - oldExitPr1 - oldExitPr2;

				// Calculate delta_L(M) = L(M)_new - L(M)_old
				double delta_allExit_log_allExit = pLogP(
						currentResult.newSumExitPr) - pLogP(sumAllExitPr);
				double delta_exit_log_exit = pLogP(currentResult.exitPr1)
						+ pLogP(currentResult.exitPr2) - pLogP(oldExitPr1)
						- pLogP(oldExitPr2);
				double delta_stay_log_stay = pLogP(
						currentResult.exitPr1 + currentResult.sumPr1)
						+ pLogP(currentResult.exitPr2 + currentResult.sumPr2)
						- pLogP(oldExitPr1 + oldSumPr1)
						- pLogP(oldExitPr2 + oldSumPr2);

				// delta_L(M) = delta_allExit - 2.0 * delta_exit_log_exit + delta_stay_log_stay.
				currentResult.diffCodeLen = delta_allExit_log_allExit
						- 2.0 * delta_exit_log_exit + delta_stay_log_stay;

				if (currentResult.diffCodeLen < bestResult.diffCodeLen) {
					// we need to update bestResult with currentResult.
					bestResult.diffCodeLen = currentResult.diffCodeLen;
					bestResult.newModule = currentResult.newModule;
					bestResult.sumPr1 = currentResult.sumPr1;
					bestResult.sumPr2 = currentResult.sumPr2;
					bestResult.exitPr1 = currentResult.exitPr1;
					bestResult.exitPr2 = currentResult.exitPr2;
					bestResult.newSumExitPr = currentResult.newSumExitPr;
				}
			}
		}

		// Make best possible move for the current node nd.
		if (bestResult.diffCodeLen < 0.0) {
			// update related to newMod...
			int newMod = bestResult.newModule;

			if (modules[newMod].numMembers == 0) {
				newMod = emptyModules.back();
				emptyModules.pop_back();
				nEmptyMod--;
				nModule++;
			}
			nd.setModIdx(newMod);

			modules[newMod].numMembers++;
			modules[newMod].exitPr = bestResult.exitPr2;
			modules[newMod].sumPr = bestResult.sumPr2;
			modules[newMod].stayPr = bestResult.exitPr2 + bestResult.sumPr2;
			modules[newMod].sumTPWeight += ndTPWeight;

			if (nd.IsDangling()) {
				modules[newMod].sumDangling += ndSize;
				modules[oldMod].sumDangling -= ndSize;
			}

			// update related to the oldMod...
			modules[oldMod].numMembers--;
			modules[oldMod].exitPr = bestResult.exitPr1;
			modules[oldMod].sumPr = bestResult.sumPr1;
			modules[oldMod].stayPr = bestResult.exitPr1 + bestResult.sumPr1;
			modules[oldMod].sumTPWeight -= ndTPWeight;

			if (modules[oldMod].numMembers == 0) {
				nEmptyMod++;
				nModule--;
				emptyModules.push_back(oldMod);
			}

			sumAllExitPr = bestResult.newSumExitPr;

			codeLength += bestResult.diffCodeLen;
			numMoved++;
		}
	}

	if (nodes.size() != nModule + nEmptyMod) {
		cout << "Something wrong!! nodes.size() != nModule + nEmptyMod."
				<< endl;
	}

	return numMoved;
}

/*
 * This function implements a prioritized version of move() above.
 *
 *	return: the number of movements.
 */
int Network::prioritize_move(double vThresh, int iteration, bool inWhile) {

	int size;
	int rank;
	int nActive = 0;
	int elementCount = 0;
	const unsigned int intPackSize = 3 * nNode + 3; // 3*nActive+1 index will save emptyModCount and 3*nActive+2 index will save nModuleCount
	const unsigned int doublePackSize = 7 * nNode;
	set<int> emptyModuleSet;

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	printf("*************entering prioritize move:rank:%d*************\n",
			rank);

	if (rank == 0) {
		nActive = activeNodes.size();
	}
	MPI_Bcast(&nActive, 1, MPI_INT, 0, MPI_COMM_WORLD);

	int* randomGlobalArray = new int[nActive];

	for (int i = 0; i < nActive; i++) {
		randomGlobalArray[i] = activeNodes[i];
	}

	for (int i = 0; i < nActive; i++) {
		int target = R->randInt(nActive - 1);
		// swap numbers between i and target.
		int tmp = randomGlobalArray[i];
		randomGlobalArray[i] = randomGlobalArray[target];
		randomGlobalArray[target] = tmp;
	}
	int start, end;
	findAssignedPart(&start, &end, nActive, size, rank);

	int numberOfElements = end - start;	//number of elements assigned to this processor

// Generate random sequential order of active nodes.
	int intSendPack[intPackSize] = { };	// multiplier 3 for 3 different elements of index, newModules, oldModules and +1 for elementCount
	int intReceivePack[intPackSize] = { };
	double doubleSendPack[doublePackSize] = { };// multiplier 6 for 6 different elements of diffCodeLen, sumPr1, sumPr2, exitPr1, exitPr2, newSumExitPr
	double doubleReceivePack[doublePackSize] = { };

	int l = nActive / size;

	double codeLengthReduction = 0.0;
	double currentSumAllExitPr = sumAllExitPr;

	int numMoved = 0;

	elementCount = 0;
	int emptyModCount = 0;
	int nModuleCount = 0;
	//emptyMods[rank][0] = emptyModCount;	//this is to keep the number of emptyMod for each of the process

	if (iteration == 0) {
		printf("start:%d, end:%d, rank:%d\n", start, end, rank);
	}

	for (int i = start; i < end; i++) {
		int counter = 0;
		bool moduleUpdate = false;
		Node& nd = nodes[randomGlobalArray[i]];	// look at i_th Node of the random sequential order.
		int oldMod = nd.ModIdx();
		//resultObj.elementCount = 0;
		if (iteration == 3 && randomGlobalArray[i] == 3 && inWhile == false) {
			printf(
					"unique oldMod:%d, rank:%d, nodes[%d].modIdx:%d, ModIdx():%d, randomGlobalArray[%d]:%d, nodes[randomGlobalArray[%d]].ModIdx():%d\n",
					oldMod, rank, randomGlobalArray[i],
					nodes[randomGlobalArray[i]].modIdx,
					nodes[randomGlobalArray[i]].ModIdx(), i,
					randomGlobalArray[i], i,
					nodes[randomGlobalArray[i]].ModIdx());
		}

		int nModLinks = 0;// The number of links to/from between the current node and other modules.

		flowmap outFlowToMod;	// <modID, flow> for outFlow...
		flowmap inFlowFromMod;		// <modID, flow> for inFlow...

		// count other modules that are connected to the current node.
		// During this counting, aggregate outFlowNotToOldMod and inFlowFromOldMod values.
		for (link_iterator linkIt = nd.outLinks.begin();
				linkIt != nd.outLinks.end(); linkIt++) {
			int newMod = nodes[linkIt->first].ModIdx();
			//printf("analyzing outnode:%d, outlinked node:%d, rank:%d\n",
			//nodes[randomGlobalArray[i]].ID() + 1, linkIt->first + 1,
			//rank);

			if (outFlowToMod.count(newMod) > 0) {
				outFlowToMod[newMod] += beta * linkIt->second;
			} else {
				outFlowToMod[newMod] = beta * linkIt->second;// initialization of the outFlow of the current newMod.
				inFlowFromMod[newMod] = 0.0;
				nModLinks++;
			}
		}

		for (link_iterator linkIt = nd.inLinks.begin();
				linkIt != nd.inLinks.end(); linkIt++) {
			int newMod = nodes[linkIt->first].ModIdx();
			if (inFlowFromMod.count(newMod) > 0) {
				inFlowFromMod[newMod] += beta * linkIt->second;
			} else {
				outFlowToMod[newMod] = 0.0;
				inFlowFromMod[newMod] = beta * linkIt->second;
				nModLinks++;
			}
		}

		if (nModLinks != outFlowToMod.size())
			cout
					<< "ALERT: nModLinks != outFlowToMod.size() in Network::prioritize_move()."
					<< endl;

		// copy node specific values for easy use and efficiency.
		double ndSize = nd.Size();					// p_nd.
		double ndTPWeight = nd.TeleportWeight();		// tau_nd.
		double ndDanglingSize = nd.DanglingSize();

		double oldExitPr1 = modules[oldMod].exitPr;

		double oldSumPr1 = modules[oldMod].sumPr;
		double oldSumDangling1 = modules[oldMod].sumDangling;
		double oldModTPWeight = modules[oldMod].sumTPWeight;

		double additionalTeleportOutFlow = (alpha * ndSize
				+ beta * ndDanglingSize) * (oldModTPWeight - ndTPWeight);
		double additionalTeleportInFlow = (alpha * (oldSumPr1 - ndSize)
				+ beta * (oldSumDangling1 - ndDanglingSize)) * ndTPWeight;

		int newMod;
		// For teleportation and danling nodes.
		for (flowmap::iterator it = outFlowToMod.begin();
				it != outFlowToMod.end(); it++) {
			newMod = it->first;
			if (newMod == oldMod) {
				outFlowToMod[newMod] += additionalTeleportOutFlow;
				inFlowFromMod[newMod] += additionalTeleportInFlow;
			} else {
				outFlowToMod[newMod] += (alpha * ndSize + beta * ndDanglingSize)
						* modules[newMod].sumTPWeight;
				inFlowFromMod[newMod] += (alpha * modules[newMod].sumPr
						+ beta * modules[newMod].sumDangling) * ndTPWeight;
			}
		}

		// Calculate flow to/from own module (default value if no links to own module).
		double outFlowToOldMod = additionalTeleportOutFlow;
		double inFlowFromOldMod = additionalTeleportInFlow;
		if (outFlowToMod.count(oldMod) > 0) {
			outFlowToOldMod = outFlowToMod[oldMod];
			inFlowFromOldMod = inFlowFromMod[oldMod];
		}

		MoveSummary currentResult;
		MoveSummary bestResult;

		double newExitPr1 = oldExitPr1 - nd.ExitPr() + outFlowToOldMod
				+ inFlowFromOldMod;

		bestResult.diffCodeLen = 0.0; // This is the default value, if we can't find diffCodeLen < 0, then don't move the node.
		bestResult.newModule = oldMod;

		counter++;
		for (flowmap::iterator it = outFlowToMod.begin();
				it != outFlowToMod.end(); it++) {
			int newMod = it->first;
			double outFlowToNewMod = it->second;
			double inFlowFromNewMod = inFlowFromMod[newMod];

			if (newMod != oldMod) {
				// copy module specific values...
				double oldExitPr2 = modules[newMod].exitPr;
				double oldSumPr2 = modules[newMod].sumPr;

				// Calculate status of current investigated movement of the node nd.
				currentResult.newModule = newMod;
				currentResult.sumPr1 = oldSumPr1 - ndSize; // This should be 0.0, because oldModule will be empty module.
				currentResult.sumPr2 = oldSumPr2 + ndSize;
				currentResult.exitPr1 = newExitPr1;
				currentResult.exitPr2 = oldExitPr2 + nd.ExitPr()
						- outFlowToNewMod - inFlowFromNewMod;

				currentResult.newSumExitPr = sumAllExitPr + newExitPr1
						+ currentResult.exitPr2 - oldExitPr1 - oldExitPr2;

				// Calculate delta_L(M) = L(M)_new - L(M)_old
				double delta_allExit_log_allExit = pLogP(
						currentResult.newSumExitPr) - pLogP(sumAllExitPr);
				double delta_exit_log_exit = pLogP(currentResult.exitPr1)
						+ pLogP(currentResult.exitPr2) - pLogP(oldExitPr1)
						- pLogP(oldExitPr2);
				double delta_stay_log_stay = pLogP(
						currentResult.exitPr1 + currentResult.sumPr1)
						+ pLogP(currentResult.exitPr2 + currentResult.sumPr2)
						- pLogP(oldExitPr1 + oldSumPr1)
						- pLogP(oldExitPr2 + oldSumPr2);

				// delta_L(M) = delta_allExit - 2.0 * delta_exit_log_exit + delta_stay_log_stay.
				currentResult.diffCodeLen = delta_allExit_log_allExit
						- 2.0 * delta_exit_log_exit + delta_stay_log_stay;

				if (iteration == 3 && randomGlobalArray[i] == 3
						&& inWhile == false) {
					printf(
							"specific oldMod:%d, rank:%d, nodes[%d].modIdx:%d, ModIdx():%d, randomGlobalArray[%d]:%d, nodes[randomGlobalArray[%d]].ModIdx():%d, newMod:%d, currentResult.diffCodeLen:%f, bestResult.diffCodeLen:%f\n",
							oldMod, rank, randomGlobalArray[i],
							nodes[randomGlobalArray[i]].modIdx,
							nodes[randomGlobalArray[i]].ModIdx(), i,
							randomGlobalArray[i], i,
							nodes[randomGlobalArray[i]].ModIdx(), newMod,
							currentResult.diffCodeLen, bestResult.diffCodeLen);
				}

				if (currentResult.diffCodeLen < bestResult.diffCodeLen) {
					// we need to update bestResult with currentResult.

					// update the moduleUpdate flag to indicate that an update has been made for the current vertex i
					moduleUpdate = true;
					//save index
					intSendPack[elementCount] = i;

					//save old module id
					intSendPack[1 * nNode + elementCount] = oldMod;
					//save new module id
					bestResult.newModule =
							intSendPack[2 * nNode + elementCount] =
									currentResult.newModule;

					if (iteration == 3 && randomGlobalArray[i] == 3
							&& inWhile == false) {
						printf(
								"special oldMod:%d, rank:%d, nodes[%d].modIdx:%d, ModIdx():%d, randomGlobalArray[%d]:%d, nodes[randomGlobalArray[%d]].ModIdx():%d, newMod:%d\n",
								oldMod, rank, randomGlobalArray[i],
								nodes[randomGlobalArray[i]].modIdx,
								nodes[randomGlobalArray[i]].ModIdx(), i,
								randomGlobalArray[i], i,
								nodes[randomGlobalArray[i]].ModIdx(),
								bestResult.newModule);
					}

					bestResult.diffCodeLen = doubleSendPack[elementCount] =
							currentResult.diffCodeLen;

					bestResult.sumPr1 = doubleSendPack[1 * nActive
							+ elementCount] = currentResult.sumPr1;

					bestResult.sumPr2 = doubleSendPack[2 * nActive
							+ elementCount] = currentResult.sumPr2;

					bestResult.exitPr1 = doubleSendPack[3 * nActive
							+ elementCount] = currentResult.exitPr1;

					bestResult.exitPr2 = doubleSendPack[4 * nActive
							+ elementCount] = currentResult.exitPr2;

					bestResult.newSumExitPr = currentResult.newSumExitPr;

				}
			} else { //this element/vertex is not moving to a new module, so both the newModule and the oldModule are same for this element
				intSendPack[elementCount] = i;

				// update the moduleUpdate flag to indicate that an update has been made for the current vertex i
				moduleUpdate = true;

				intSendPack[1 * nNode + elementCount] = oldMod;
				intSendPack[2 * nNode + elementCount] = bestResult.newModule;
				if (iteration == 3 && randomGlobalArray[i] == 3
						&& inWhile == false) {
					printf(
							"speci oldMod:%d, rank:%d, i:%d, nodes[%d].modIdx:%d, ModIdx():%d, randomGlobalArray[%d]:%d, nodes[randomGlobalArray[%d]].ModIdx():%d, newMod:%d, currentResult.diffCodeLen:%f, bestResult.diffCodeLen:%f, bestResult.newModule:%d\n",
							oldMod, rank, i, randomGlobalArray[i],
							nodes[randomGlobalArray[i]].modIdx,
							nodes[randomGlobalArray[i]].ModIdx(), i,
							randomGlobalArray[i], i,
							nodes[randomGlobalArray[i]].ModIdx(), newMod,
							currentResult.diffCodeLen, bestResult.diffCodeLen,
							bestResult.newModule);
				}
			}
		}

		// Make best possible move for the current node nd.
		//if (bestResult.diffCodeLen < 0.0) {
		if (bestResult.diffCodeLen < vThresh) {
			// update related to newMod...
			int newMod = bestResult.newModule;

			nd.setModIdx(newMod);

			modules[newMod].exitPr = bestResult.exitPr2;
			modules[newMod].sumPr = bestResult.sumPr2;
			modules[newMod].stayPr = bestResult.exitPr2 + bestResult.sumPr2;

			modules[newMod].sumTPWeight += ndTPWeight;
			//doubleSendPack[5 * nNode + newMod] += ndTPWeight;

			if (nd.IsDangling()) {
				modules[newMod].sumDangling += ndSize;
				modules[oldMod].sumDangling -= ndSize;
			}

			modules[oldMod].exitPr = bestResult.exitPr1;
			modules[oldMod].sumPr = bestResult.sumPr1;
			modules[oldMod].stayPr = bestResult.exitPr1 + bestResult.sumPr1;

			modules[oldMod].sumTPWeight -= ndTPWeight;

			sumAllExitPr = bestResult.newSumExitPr;

			/*			if (iteration == 0) {
			 printf("best result sum all exit pr:%f rank:%d\n",
			 bestResult.newSumExitPr, rank);
			 }*/

			codeLengthReduction += bestResult.diffCodeLen;
			//codeLength += bestResult.diffCodeLen;

			numMoved++;

			// update activeNodes and isActives vectors.
			// We have to add the following nodes in activeNodes: neighbors, members in oldMod & newMod.
			for (link_iterator linkIt = nd.outLinks.begin();
					linkIt != nd.outLinks.end(); linkIt++)
				isActives[linkIt->first] = 1;	// set as an active nodes.

			for (link_iterator linkIt = nd.inLinks.begin();
					linkIt != nd.inLinks.end(); linkIt++)
				isActives[linkIt->first] = 1;	// set as an active nodes.
		}

		if (moduleUpdate == false) {
			intSendPack[elementCount] = i;
			intSendPack[1 * nNode + elementCount] = oldMod;
			intSendPack[2 * nNode + elementCount] = oldMod;
		}
		elementCount++;
	}

	if (numberOfElements != elementCount) {
		if (iteration == 0) {
			printf(
					"something is not right, for rank:%d, numberOfElements:%d does not match elementCount:%d\n",
					rank, numberOfElements, elementCount);
		}
	}

	intSendPack[intPackSize - 1] = elementCount;
	doubleSendPack[doublePackSize - 1] = sumAllExitPr - currentSumAllExitPr;
	codeLength += codeLengthReduction;
	doubleSendPack[doublePackSize - 2] = codeLengthReduction;

	for (int prId = 0; prId < size; prId++) {
		if (prId != rank) {
			MPI_Send(intSendPack, intPackSize, MPI_INT, prId, iteration,
			MPI_COMM_WORLD);
			MPI_Send(doubleSendPack, doublePackSize, MPI_DOUBLE, prId,
					iteration,
					MPI_COMM_WORLD);
		}
	}

	for (int sId = 0; sId < size; sId++) {
		if (sId != rank) {
			MPI_Recv(intReceivePack, intPackSize, MPI_INT, sId, iteration,
			MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
			MPI_Recv(doubleReceivePack, doublePackSize, MPI_DOUBLE, sId,
					iteration,
					MPI_COMM_WORLD, MPI_STATUSES_IGNORE);

			int totalElementSentFromSender = intReceivePack[intPackSize - 1];

			sumAllExitPr += doubleReceivePack[doublePackSize - 1];
			codeLength += doubleReceivePack[doublePackSize - 2];

			for (int z = 0; z < totalElementSentFromSender; z++) {

				int nodeIndex = intReceivePack[z];

				Node& nod = nodes[randomGlobalArray[nodeIndex]];
				int oldModule = intReceivePack[1 * nNode + z];
				int newModule = intReceivePack[2 * nNode + z];

				//we need to do some tweaking here to get rid of the extra reduction of codelength because of the circular moves in distributed informap
				//int nodeIndex = nod.ID();
				nod.setModIdx(newModule);

				if (oldModule != newModule) {
					numMoved++;

					modules[newModule].exitPr = doubleReceivePack[4 * nActive
							+ z];
					modules[newModule].sumPr =
							doubleReceivePack[2 * nActive + z];
					modules[newModule].stayPr = modules[newModule].exitPr
							+ modules[newModule].sumPr;

					modules[newModule].sumTPWeight += nod.TeleportWeight();

					if (nod.IsDangling()) {
						modules[newModule].sumDangling += nod.Size();
						modules[oldModule].sumDangling -= nod.Size();
					}

					modules[oldModule].exitPr = doubleReceivePack[3 * nActive
							+ z];
					modules[oldModule].sumPr =
							doubleReceivePack[1 * nActive + z];
					modules[oldModule].stayPr = modules[oldModule].exitPr
							+ modules[oldModule].sumPr;

					modules[oldModule].sumTPWeight -= nod.TeleportWeight();

					// update activeNodes and isActives vectors.
					// We have to add the following nodes in activeNodes: neighbors, members in oldMod & newMod.
					for (link_iterator linkIt = nod.outLinks.begin();
							linkIt != nod.outLinks.end(); linkIt++)
						isActives[linkIt->first] = 1;// set as an active nodes.

					for (link_iterator linkIt = nod.inLinks.begin();
							linkIt != nod.inLinks.end(); linkIt++)
						isActives[linkIt->first] = 1;// set as an active nodes.
				}
			}

		}
	}

	vector<int>().swap(activeNodes);
	for (int i = 0; i < isActives.size(); i++) {
		if (isActives[i] == 1) {
			activeNodes.push_back(i);
			isActives[i] = 0;	// reset the flag of isActives[i].
		}
	}

	/*	if (iteration == 3 && inWhile == false) {
	 for (int i = 0; i < nNode; i++) {
	 printf(
	 "Inside prioritize_move distributed output for rank:%d, nodes[%d], module Id:%d\n",
	 rank, i, nodes[i].modIdx);
	 }
	 }*/

	printf("*************leaving prioritize move:rank:%d*************\n", rank);

	delete[] randomGlobalArray;

	return numMoved;
}

/*
 *	+ This is a parallel implementation of Network::move() function via OpenMP.
 *		We will loose the criteria for simple implementation and it may help to avoid local optima problem.
 */
int Network::parallelMove(int numTh, double& tSequential) {

	struct timeval tStart, tEnd;

	gettimeofday(&tStart, NULL);

// Generate random sequential order of nodes.
	vector<int> randomOrder(nNode);
	for (int i = 0; i < nNode; i++)
		randomOrder[i] = i;

	for (int i = 0; i < nNode; i++) {
		int target = R->randInt(nNode - 1);

		// swap numbers between i and target.
		int tmp = randomOrder[i];
		randomOrder[i] = randomOrder[target];
		randomOrder[target] = tmp;
	}

	int numMoved = 0;

	omp_set_num_threads(numTh);

	const int emptyTarget = nNode + 1;// This will be an indicator of moving to emptyModule.

	gettimeofday(&tEnd, NULL);
	tSequential += elapsedTimeInSec(tStart, tEnd);

#pragma omp parallel for 
// Move each node to one of its neighbor modules in random sequential order.
	for (int i = 0; i < nNode; i++) {
		Node& nd = nodes[randomOrder[i]];// look at i_th Node of the random sequential order.
		int oldMod = nd.ModIdx();

		int nModLinks = 0;// The number of links to/from between the current node and other modules.

		flowmap outFlowToMod;	// <modID, flow> for outFlow...
		flowmap inFlowFromMod;	// <modID, flow> for inFlow...

		// count other modules that are connected to the current node.
		for (link_iterator linkIt = nd.outLinks.begin();
				linkIt != nd.outLinks.end(); linkIt++) {
			int newMod = nodes[linkIt->first].ModIdx();

			if (outFlowToMod.count(newMod) > 0) {
				outFlowToMod[newMod] += beta * linkIt->second;
			} else {
				outFlowToMod[newMod] = beta * linkIt->second;// initialization of the outFlow of the current newMod.
				inFlowFromMod[newMod] = 0.0;
				nModLinks++;
			}
		}

		for (link_iterator linkIt = nd.inLinks.begin();
				linkIt != nd.inLinks.end(); linkIt++) {
			int newMod = nodes[linkIt->first].ModIdx();

			if (inFlowFromMod.count(newMod) > 0) {
				inFlowFromMod[newMod] += beta * linkIt->second;
			} else {
				outFlowToMod[newMod] = 0.0;
				inFlowFromMod[newMod] = beta * linkIt->second;
				nModLinks++;
			}
		}

		if (nModLinks != outFlowToMod.size())
			cout
					<< "ALERT: nModLinks != outFlowToMod.size() in Network::move()."
					<< endl;

		// copy node specific values for easy use and efficiency.
		double ndSize = nd.Size();	// p_nd.
		double ndTPWeight = nd.TeleportWeight();	// tau_nd.
		double ndDanglingSize = nd.DanglingSize();

		// Below values are related current module information, but it could be changed in the middle of this decision process.
		// However, we used this snapshot for finding next module of current node.
		// These will be a guideline of decision process and the correct values will be calculated at the end of this iteration.
		double oldExitPr1 = modules[oldMod].exitPr;
		double oldSumPr1 = modules[oldMod].sumPr;
		double oldSumDangling1 = modules[oldMod].sumDangling;
		double oldModTPWeight = modules[oldMod].sumTPWeight;

		double additionalTeleportOutFlow = (alpha * ndSize
				+ beta * ndDanglingSize) * (oldModTPWeight - ndTPWeight);
		double additionalTeleportInFlow = (alpha * (oldSumPr1 - ndSize)
				+ beta * (oldSumDangling1 - ndDanglingSize)) * ndTPWeight;

		// For teleportation and danling nodes.
		for (flowmap::iterator it = outFlowToMod.begin();
				it != outFlowToMod.end(); it++) {
			int newMod = it->first;
			if (newMod == oldMod) {
				outFlowToMod[newMod] += additionalTeleportOutFlow;
				inFlowFromMod[newMod] += additionalTeleportInFlow;
			} else {
				outFlowToMod[newMod] += (alpha * ndSize + beta * ndDanglingSize)
						* modules[newMod].sumTPWeight;
				inFlowFromMod[newMod] += (alpha * modules[newMod].sumPr
						+ beta * modules[newMod].sumDangling) * ndTPWeight;
			}
		}

		// Calculate flow to/from own module (default value if no links to own module).
		double outFlowToOldMod = additionalTeleportOutFlow;
		double inFlowFromOldMod = additionalTeleportInFlow;
		if (outFlowToMod.count(oldMod) > 0) {
			outFlowToOldMod = outFlowToMod[oldMod];
			inFlowFromOldMod = inFlowFromMod[oldMod];
		}

		//////////////////// THE OPTION TO MOVE TO EMPTY MODULE ////////////////
		if (modules[oldMod].members.size() > 1 && emptyModules.size() > 0) {
			outFlowToMod[emptyTarget] = 0.0;
			inFlowFromMod[emptyTarget] = 0.0;
			nModLinks++;
		}

		MoveSummary currentResult;
		MoveSummary bestResult;

		double newExitPr1 = oldExitPr1 - nd.ExitPr() + outFlowToOldMod
				+ inFlowFromOldMod;

		bestResult.diffCodeLen = 0.0;// This is the default value, if we can't find diffCodeLen < 0, then don't move the node.

		for (flowmap::iterator it = outFlowToMod.begin();
				it != outFlowToMod.end(); it++) {
			int newMod = it->first;
			double outFlowToNewMod = it->second;
			double inFlowFromNewMod = inFlowFromMod[newMod];

			if (newMod != oldMod) {

				// copy module specific values...
				double oldExitPr2 = modules[newMod].exitPr;
				double oldSumPr2 = modules[newMod].sumPr;

				// Calculate status of current investigated movement of the node nd.
				currentResult.newModule = newMod;
				currentResult.sumPr1 = oldSumPr1 - ndSize;// This should be 0.0, because oldModule will be empty module.
				currentResult.sumPr2 = oldSumPr2 + ndSize;
				currentResult.exitPr1 = newExitPr1;
				currentResult.exitPr2 = oldExitPr2 + nd.ExitPr()
						- outFlowToNewMod - inFlowFromNewMod;

				currentResult.newSumExitPr = sumAllExitPr + newExitPr1
						+ currentResult.exitPr2 - oldExitPr1 - oldExitPr2;

				// Calculate delta_L(M) = L(M)_new - L(M)_old
				double delta_allExit_log_allExit = pLogP(
						currentResult.newSumExitPr) - pLogP(sumAllExitPr);
				double delta_exit_log_exit = pLogP(currentResult.exitPr1)
						+ pLogP(currentResult.exitPr2) - pLogP(oldExitPr1)
						- pLogP(oldExitPr2);
				double delta_stay_log_stay = pLogP(
						currentResult.exitPr1 + currentResult.sumPr1)
						+ pLogP(currentResult.exitPr2 + currentResult.sumPr2)
						- pLogP(oldExitPr1 + oldSumPr1)
						- pLogP(oldExitPr2 + oldSumPr2);

				// delta_L(M) = delta_allExit - 2.0 * delta_exit_log_exit + delta_stay_log_stay.
				currentResult.diffCodeLen = delta_allExit_log_allExit
						- 2.0 * delta_exit_log_exit + delta_stay_log_stay;

				if (currentResult.diffCodeLen < bestResult.diffCodeLen) {
					// we need to update bestResult with currentResult.
					bestResult.diffCodeLen = currentResult.diffCodeLen;
					bestResult.newModule = currentResult.newModule;
				}
			}
		}

		// store the best possilbe movement information if necessary.
		// In this version, we are trying to move as soon as possible, i.e. right after decision making...
		if (bestResult.diffCodeLen < 0.0) {

			bool isEmptyTarget = false;
			bool validMove = true;// This will indicate the validity of the decided movement.
			int newMod = bestResult.newModule;

#pragma omp critical (moveUpdate)
			{
				gettimeofday(&tStart, NULL);

				// if newMod == emptyTarget, it indicates moves to empty module.
				if ((nEmptyMod > 0) && (newMod == emptyTarget)
						&& (modules[oldMod].numMembers > 1)) {
					newMod = emptyModules.back();
					isEmptyTarget = true;
				} else if (newMod == emptyTarget) {
					validMove = false;
				} else if (modules[newMod].numMembers == 0) {
					// This is the case that the algorithm thought there are some nodes in the new module since newMod != emptyTarget,
					// but the nodes are all moved to other modules so there are no nodes in there.
					// Thus, we don't know whether generating a new module will be better or not.
					// Therefore, don't move this option.
					//continue;
					validMove = false;
				}

				MoveSummary moveResult;

				if (validMove) {

					////////////////////////////////////////////////////////////////////
					/// THIS IS THE PART FOR EXAMINING QUALITY IMPROVEMENT OR NOT... ///
					////////////////////////////////////////////////////////////////////

					outFlowToOldMod = 0.0;
					double outFlowToNewMod = 0.0;
					inFlowFromOldMod = 0.0;
					double inFlowFromNewMod = 0.0;

					if (!isEmptyTarget) {
						for (link_iterator linkIt = nd.outLinks.begin();
								linkIt != nd.outLinks.end(); linkIt++) {
							int toMod = nodes[linkIt->first].ModIdx();

							if (toMod == oldMod) {
								outFlowToOldMod += beta * linkIt->second;
							} else if (toMod == newMod) {
								outFlowToNewMod += beta * linkIt->second;
							}
						}

						for (link_iterator linkIt = nd.inLinks.begin();
								linkIt != nd.inLinks.end(); linkIt++) {
							int fromMod = nodes[linkIt->first].ModIdx();

							if (fromMod == oldMod) {
								inFlowFromOldMod += beta * linkIt->second;
							} else if (fromMod == newMod) {
								inFlowFromNewMod += beta * linkIt->second;
							}
						}
					} else {
						for (link_iterator linkIt = nd.outLinks.begin();
								linkIt != nd.outLinks.end(); linkIt++) {
							int toMod = nodes[linkIt->first].ModIdx();
							if (toMod == oldMod) {
								outFlowToOldMod += beta * linkIt->second;
							}
						}

						for (link_iterator linkIt = nd.inLinks.begin();
								linkIt != nd.inLinks.end(); linkIt++) {
							int fromMod = nodes[linkIt->first].ModIdx();
							if (fromMod == oldMod) {
								inFlowFromOldMod += beta * linkIt->second;
							}
						}
					}

					oldExitPr1 = modules[oldMod].exitPr;
					oldSumPr1 = modules[oldMod].sumPr;
					oldSumDangling1 = modules[oldMod].sumDangling;
					oldModTPWeight = modules[oldMod].sumTPWeight;

					// For teleportation and danling nodes.
					outFlowToOldMod += (alpha * ndSize + beta * ndDanglingSize)
							* (oldModTPWeight - ndTPWeight);
					inFlowFromOldMod += (alpha * (oldSumPr1 - ndSize)
							+ beta * (oldSumDangling1 - ndDanglingSize))
							* ndTPWeight;
					outFlowToNewMod += (alpha * ndSize + beta * ndDanglingSize)
							* modules[newMod].sumTPWeight;
					inFlowFromNewMod += (alpha * modules[newMod].sumPr
							+ beta * modules[newMod].sumDangling) * ndTPWeight;

					if (isEmptyTarget) {
						outFlowToNewMod = 0.0;
						inFlowFromNewMod = 0.0;
					}

					moveResult.exitPr1 = oldExitPr1 - nd.ExitPr()
							+ outFlowToOldMod + inFlowFromOldMod;

					// copy module specific values...
					double oldExitPr2 = modules[newMod].exitPr;
					double oldSumPr2 = modules[newMod].sumPr;

					// Calculate status of current investigated movement of the node nd.
					moveResult.newModule = newMod;
					moveResult.sumPr1 = oldSumPr1 - ndSize;	// This should be 0.0, because oldModule will be empty module.
					moveResult.sumPr2 = oldSumPr2 + ndSize;
					moveResult.exitPr2 = oldExitPr2 + nd.ExitPr()
							- outFlowToNewMod - inFlowFromNewMod;

					moveResult.newSumExitPr = sumAllExitPr + moveResult.exitPr1
							+ moveResult.exitPr2 - oldExitPr1 - oldExitPr2;

					// Calculate delta_L(M) = L(M)_new - L(M)_old
					double delta_allExit_log_allExit = pLogP(
							moveResult.newSumExitPr) - pLogP(sumAllExitPr);
					double delta_exit_log_exit = pLogP(moveResult.exitPr1)
							+ pLogP(moveResult.exitPr2) - pLogP(oldExitPr1)
							- pLogP(oldExitPr2);
					double delta_stay_log_stay = pLogP(
							moveResult.exitPr1 + moveResult.sumPr1)
							+ pLogP(moveResult.exitPr2 + moveResult.sumPr2)
							- pLogP(oldExitPr1 + oldSumPr1)
							- pLogP(oldExitPr2 + oldSumPr2);

					moveResult.diffCodeLen = delta_allExit_log_allExit
							- 2.0 * delta_exit_log_exit + delta_stay_log_stay;

					////////////////////////////////////
					////// THE END OF EXAMINATION //////
					////////////////////////////////////

					if (isEmptyTarget) {
						emptyModules.pop_back();

						nEmptyMod--;
						nModule++;
					}

					nd.setModIdx(newMod);

					modules[newMod].numMembers++;
					modules[newMod].exitPr = moveResult.exitPr2;
					modules[newMod].sumPr = moveResult.sumPr2;
					modules[newMod].stayPr = moveResult.exitPr2
							+ moveResult.sumPr2;
					modules[newMod].sumTPWeight += ndTPWeight;

					if (nd.IsDangling()) {
						modules[newMod].sumDangling += ndDanglingSize;
						modules[oldMod].sumDangling -= ndDanglingSize;
					}

					// update related to the oldMod...
					modules[oldMod].numMembers--;
					modules[oldMod].exitPr = moveResult.exitPr1;
					modules[oldMod].sumPr = moveResult.sumPr1;
					modules[oldMod].stayPr = moveResult.exitPr1
							+ moveResult.sumPr1;
					modules[oldMod].sumTPWeight -= ndTPWeight;

					if (modules[oldMod].numMembers == 0) {
						nEmptyMod++;
						nModule--;
						emptyModules.push_back(oldMod);
					}

					sumAllExitPr = moveResult.newSumExitPr;
					codeLength += moveResult.diffCodeLen;

					numMoved++;	// This is inside of critical section, so no race condition happens.
				}

				gettimeofday(&tEnd, NULL);

				tSequential += elapsedTimeInSec(tStart, tEnd);
			}	// END critical

		}

	}	// END parallel for

	if (nodes.size() != nModule + nEmptyMod) {
		cout << "Something wrong!! nodes.size() != nModule + nEmptyMod."
				<< endl;
	}

	return numMoved;
}

/*
 *	+ This is a parallel implementation of Network::prioritize_move() function via OpenMP.
 *		We will loose the criteria for simple implementation and it may help to avoid local optima problem.
 */
int Network::prioritize_parallelMove(int numTh, double& tSequential,
		double vThresh) {

	struct timeval tStart, tEnd;

	gettimeofday(&tStart, NULL);

	int nActive = activeNodes.size();
	int nNextActive = 0;

// Generate random sequential order of nodes.
	vector<int> randomOrder(nActive);
	for (int i = 0; i < nActive; i++)
		randomOrder[i] = activeNodes[i];

	for (int i = 0; i < nActive; i++) {
		int target = R->randInt(nActive - 1);

		// swap numbers between i and target.
		int tmp = randomOrder[i];
		randomOrder[i] = randomOrder[target];
		randomOrder[target] = tmp;
	}

// Now randomOrder vector already had all activeNodes info,
// so we can reset activeNodes vector for adding for new active nodes.
	vector<int>().swap(activeNodes);

	int numMoved = 0;

	omp_set_num_threads(numTh);

	const int emptyTarget = nNode + 1;// This will be an indicator of moving to emptyModule.

	gettimeofday(&tEnd, NULL);
	tSequential += elapsedTimeInSec(tStart, tEnd);

#pragma omp parallel for 
// Move each node to one of its neighbor modules in random sequential order.
	for (int i = 0; i < nActive; i++) {
		Node& nd = nodes[randomOrder[i]];// look at i_th Node of the random sequential order.
		int oldMod = nd.ModIdx();

		int nModLinks = 0;// The number of links to/from between the current node and other modules.

		flowmap outFlowToMod;	// <modID, flow> for outFlow...
		flowmap inFlowFromMod;	// <modID, flow> for inFlow...

		// count other modules that are connected to the current node.
		for (link_iterator linkIt = nd.outLinks.begin();
				linkIt != nd.outLinks.end(); linkIt++) {
			int newMod = nodes[linkIt->first].ModIdx();

			if (outFlowToMod.count(newMod) > 0) {
				outFlowToMod[newMod] += beta * linkIt->second;
			} else {
				outFlowToMod[newMod] = beta * linkIt->second;// initialization of the outFlow of the current newMod.
				inFlowFromMod[newMod] = 0.0;
				nModLinks++;
			}
		}

		for (link_iterator linkIt = nd.inLinks.begin();
				linkIt != nd.inLinks.end(); linkIt++) {
			int newMod = nodes[linkIt->first].ModIdx();

			if (inFlowFromMod.count(newMod) > 0) {
				inFlowFromMod[newMod] += beta * linkIt->second;
			} else {
				outFlowToMod[newMod] = 0.0;
				inFlowFromMod[newMod] = beta * linkIt->second;
				nModLinks++;
			}
		}

		if (nModLinks != outFlowToMod.size())
			cout
					<< "ALERT: nModLinks != outFlowToMod.size() in Network::move()."
					<< endl;

		// copy node specific values for easy use and efficiency.
		double ndSize = nd.Size();	// p_nd.
		double ndTPWeight = nd.TeleportWeight();	// tau_nd.
		double ndDanglingSize = nd.DanglingSize();

		// Below values are related current module information, but it could be changed in the middle of this decision process.
		// However, we used this snapshot for finding next module of current node.
		// These will be a guideline of decision process and the correct values will be calculated at the end of this iteration.
		double oldExitPr1 = modules[oldMod].exitPr;
		double oldSumPr1 = modules[oldMod].sumPr;
		double oldSumDangling1 = modules[oldMod].sumDangling;
		double oldModTPWeight = modules[oldMod].sumTPWeight;

		double additionalTeleportOutFlow = (alpha * ndSize
				+ beta * ndDanglingSize) * (oldModTPWeight - ndTPWeight);
		double additionalTeleportInFlow = (alpha * (oldSumPr1 - ndSize)
				+ beta * (oldSumDangling1 - ndDanglingSize)) * ndTPWeight;

		// For teleportation and danling nodes.
		for (flowmap::iterator it = outFlowToMod.begin();
				it != outFlowToMod.end(); it++) {
			int newMod = it->first;
			if (newMod == oldMod) {
				outFlowToMod[newMod] += additionalTeleportOutFlow;
				inFlowFromMod[newMod] += additionalTeleportInFlow;
			} else {
				outFlowToMod[newMod] += (alpha * ndSize + beta * ndDanglingSize)
						* modules[newMod].sumTPWeight;
				inFlowFromMod[newMod] += (alpha * modules[newMod].sumPr
						+ beta * modules[newMod].sumDangling) * ndTPWeight;
			}
		}

		// Calculate flow to/from own module (default value if no links to own module).
		double outFlowToOldMod = additionalTeleportOutFlow;
		double inFlowFromOldMod = additionalTeleportInFlow;
		if (outFlowToMod.count(oldMod) > 0) {
			outFlowToOldMod = outFlowToMod[oldMod];
			inFlowFromOldMod = inFlowFromMod[oldMod];
		}

		//////////////////// THE OPTION TO MOVE TO EMPTY MODULE ////////////////
		if (modules[oldMod].members.size() > 1 && emptyModules.size() > 0) {
			outFlowToMod[emptyTarget] = 0.0;
			inFlowFromMod[emptyTarget] = 0.0;
			nModLinks++;
		}

		MoveSummary currentResult;
		MoveSummary bestResult;

		double newExitPr1 = oldExitPr1 - nd.ExitPr() + outFlowToOldMod
				+ inFlowFromOldMod;

		bestResult.diffCodeLen = 0.0;// This is the default value, if we can't find diffCodeLen < 0, then don't move the node.

		for (flowmap::iterator it = outFlowToMod.begin();
				it != outFlowToMod.end(); it++) {
			int newMod = it->first;
			double outFlowToNewMod = it->second;
			double inFlowFromNewMod = inFlowFromMod[newMod];

			if (newMod != oldMod) {

				// copy module specific values...
				double oldExitPr2 = modules[newMod].exitPr;
				double oldSumPr2 = modules[newMod].sumPr;

				// Calculate status of current investigated movement of the node nd.
				currentResult.newModule = newMod;
				currentResult.sumPr1 = oldSumPr1 - ndSize;// This should be 0.0, because oldModule will be empty module.
				currentResult.sumPr2 = oldSumPr2 + ndSize;
				currentResult.exitPr1 = newExitPr1;
				currentResult.exitPr2 = oldExitPr2 + nd.ExitPr()
						- outFlowToNewMod - inFlowFromNewMod;

				currentResult.newSumExitPr = sumAllExitPr + newExitPr1
						+ currentResult.exitPr2 - oldExitPr1 - oldExitPr2;

				// Calculate delta_L(M) = L(M)_new - L(M)_old
				double delta_allExit_log_allExit = pLogP(
						currentResult.newSumExitPr) - pLogP(sumAllExitPr);
				double delta_exit_log_exit = pLogP(currentResult.exitPr1)
						+ pLogP(currentResult.exitPr2) - pLogP(oldExitPr1)
						- pLogP(oldExitPr2);
				double delta_stay_log_stay = pLogP(
						currentResult.exitPr1 + currentResult.sumPr1)
						+ pLogP(currentResult.exitPr2 + currentResult.sumPr2)
						- pLogP(oldExitPr1 + oldSumPr1)
						- pLogP(oldExitPr2 + oldSumPr2);

				// delta_L(M) = delta_allExit - 2.0 * delta_exit_log_exit + delta_stay_log_stay.
				currentResult.diffCodeLen = delta_allExit_log_allExit
						- 2.0 * delta_exit_log_exit + delta_stay_log_stay;

				if (currentResult.diffCodeLen < bestResult.diffCodeLen) {
					// we need to update bestResult with currentResult.
					bestResult.diffCodeLen = currentResult.diffCodeLen;
					bestResult.newModule = currentResult.newModule;
				}
			}
		}

		// store the best possilbe movement information if necessary.
		// In this version, we are trying to move as soon as possible, i.e. right after decision making...
		//if (bestResult.diffCodeLen < 0.0) {
		if (bestResult.diffCodeLen < vThresh) {

			bool isEmptyTarget = false;
			bool validMove = true;// This will indicate the validity of the decided movement.
			int newMod = bestResult.newModule;

#pragma omp critical (moveUpdate)
			{
				gettimeofday(&tStart, NULL);

				// if newMod == emptyTarget, it indicates moves to empty module.
				if ((nEmptyMod > 0) && (newMod == emptyTarget)
						&& (modules[oldMod].numMembers > 1)) {
					newMod = emptyModules.back();
					isEmptyTarget = true;
				} else if (newMod == emptyTarget) {
					validMove = false;
				} else if (modules[newMod].numMembers == 0) {
					// This is the case that the algorithm thought there are some nodes in the new module since newMod != emptyTarget,
					// but the nodes are all moved to other modules so there are no nodes in there.
					// Thus, we don't know whether generating a new module will be better or not.
					// Therefore, don't move this option.
					//continue;
					validMove = false;
				}

				MoveSummary moveResult;

				if (validMove) {

					////////////////////////////////////////////////////////////////////
					/// THIS IS THE PART FOR EXAMINING QUALITY IMPROVEMENT OR NOT... ///
					////////////////////////////////////////////////////////////////////

					outFlowToOldMod = 0.0;
					double outFlowToNewMod = 0.0;
					inFlowFromOldMod = 0.0;
					double inFlowFromNewMod = 0.0;

					if (!isEmptyTarget) {
						for (link_iterator linkIt = nd.outLinks.begin();
								linkIt != nd.outLinks.end(); linkIt++) {
							int toMod = nodes[linkIt->first].ModIdx();

							if (toMod == oldMod) {
								outFlowToOldMod += beta * linkIt->second;
							} else if (toMod == newMod) {
								outFlowToNewMod += beta * linkIt->second;
							}
						}

						for (link_iterator linkIt = nd.inLinks.begin();
								linkIt != nd.inLinks.end(); linkIt++) {
							int fromMod = nodes[linkIt->first].ModIdx();

							if (fromMod == oldMod) {
								inFlowFromOldMod += beta * linkIt->second;
							} else if (fromMod == newMod) {
								inFlowFromNewMod += beta * linkIt->second;
							}
						}
					} else {
						for (link_iterator linkIt = nd.outLinks.begin();
								linkIt != nd.outLinks.end(); linkIt++) {
							int toMod = nodes[linkIt->first].ModIdx();
							if (toMod == oldMod) {
								outFlowToOldMod += beta * linkIt->second;
							}
						}

						for (link_iterator linkIt = nd.inLinks.begin();
								linkIt != nd.inLinks.end(); linkIt++) {
							int fromMod = nodes[linkIt->first].ModIdx();
							if (fromMod == oldMod) {
								inFlowFromOldMod += beta * linkIt->second;
							}
						}
					}

					oldExitPr1 = modules[oldMod].exitPr;
					oldSumPr1 = modules[oldMod].sumPr;
					oldSumDangling1 = modules[oldMod].sumDangling;
					oldModTPWeight = modules[oldMod].sumTPWeight;

					// For teleportation and danling nodes.
					outFlowToOldMod += (alpha * ndSize + beta * ndDanglingSize)
							* (oldModTPWeight - ndTPWeight);
					inFlowFromOldMod += (alpha * (oldSumPr1 - ndSize)
							+ beta * (oldSumDangling1 - ndDanglingSize))
							* ndTPWeight;
					outFlowToNewMod += (alpha * ndSize + beta * ndDanglingSize)
							* modules[newMod].sumTPWeight;
					inFlowFromNewMod += (alpha * modules[newMod].sumPr
							+ beta * modules[newMod].sumDangling) * ndTPWeight;

					if (isEmptyTarget) {
						outFlowToNewMod = 0.0;
						inFlowFromNewMod = 0.0;
					}

					moveResult.exitPr1 = oldExitPr1 - nd.ExitPr()
							+ outFlowToOldMod + inFlowFromOldMod;

					// copy module specific values...
					double oldExitPr2 = modules[newMod].exitPr;
					double oldSumPr2 = modules[newMod].sumPr;

					// Calculate status of current investigated movement of the node nd.
					moveResult.newModule = newMod;
					moveResult.sumPr1 = oldSumPr1 - ndSize;	// This should be 0.0, because oldModule will be empty module.
					moveResult.sumPr2 = oldSumPr2 + ndSize;
					moveResult.exitPr2 = oldExitPr2 + nd.ExitPr()
							- outFlowToNewMod - inFlowFromNewMod;

					moveResult.newSumExitPr = sumAllExitPr + moveResult.exitPr1
							+ moveResult.exitPr2 - oldExitPr1 - oldExitPr2;

					// Calculate delta_L(M) = L(M)_new - L(M)_old
					double delta_allExit_log_allExit = pLogP(
							moveResult.newSumExitPr) - pLogP(sumAllExitPr);
					double delta_exit_log_exit = pLogP(moveResult.exitPr1)
							+ pLogP(moveResult.exitPr2) - pLogP(oldExitPr1)
							- pLogP(oldExitPr2);
					double delta_stay_log_stay = pLogP(
							moveResult.exitPr1 + moveResult.sumPr1)
							+ pLogP(moveResult.exitPr2 + moveResult.sumPr2)
							- pLogP(oldExitPr1 + oldSumPr1)
							- pLogP(oldExitPr2 + oldSumPr2);

					// delta_L(M) = delta_allExit - 2.0 * delta_exit_log_exit + delta_stay_log_stay.
					moveResult.diffCodeLen = delta_allExit_log_allExit
							- 2.0 * delta_exit_log_exit + delta_stay_log_stay;

					////////////////////////////////////
					////// THE END OF EXAMINATION //////
					////////////////////////////////////

					if (isEmptyTarget) {
						emptyModules.pop_back();

						nEmptyMod--;
						nModule++;
					}

					nd.setModIdx(newMod);

					modules[newMod].numMembers++;
					modules[newMod].exitPr = moveResult.exitPr2;
					modules[newMod].sumPr = moveResult.sumPr2;
					modules[newMod].stayPr = moveResult.exitPr2
							+ moveResult.sumPr2;
					modules[newMod].sumTPWeight += ndTPWeight;

					if (nd.IsDangling()) {
						modules[newMod].sumDangling += ndDanglingSize;
						modules[oldMod].sumDangling -= ndDanglingSize;
					}

					// update related to the oldMod...
					modules[oldMod].numMembers--;
					modules[oldMod].exitPr = moveResult.exitPr1;
					modules[oldMod].sumPr = moveResult.sumPr1;
					modules[oldMod].stayPr = moveResult.exitPr1
							+ moveResult.sumPr1;
					modules[oldMod].sumTPWeight -= ndTPWeight;

					if (modules[oldMod].numMembers == 0) {
						nEmptyMod++;
						nModule--;
						emptyModules.push_back(oldMod);
					}

					sumAllExitPr = moveResult.newSumExitPr;
					codeLength += moveResult.diffCodeLen;

					numMoved++;
				}

				gettimeofday(&tEnd, NULL);

				tSequential += elapsedTimeInSec(tStart, tEnd);
			}	// END critical

			// update activeNodes and isActives vectors.
			// We have to add the following nodes in activeNodes: neighbors, members in oldMod & newMod.
			// This can be done in parallel without any locking, since the written value is always same.
			for (link_iterator linkIt = nd.outLinks.begin();
					linkIt != nd.outLinks.end(); linkIt++)
				isActives[linkIt->first] = 1;		// set as an active nodes.

			for (link_iterator linkIt = nd.inLinks.begin();
					linkIt != nd.inLinks.end(); linkIt++)
				isActives[linkIt->first] = 1;		// set as an active nodes.
		}

	}	// END parallel for

	for (int i = 0; i < isActives.size(); i++) {
		if (isActives[i] == 1) {
			activeNodes.push_back(i);
			isActives[i] = 0;	// reset the flag of isActives[i].
		}
	}

	if (nodes.size() != nModule + nEmptyMod) {
		cout << "Something wrong!! nodes.size() != nModule + nEmptyMod."
				<< endl;
	}

	return numMoved;
}

/*
 *	same as Network::move() except that the moving unit is SuperNode.
 */

int Network::moveSuperNodes() {

	int nSuperNodes = superNodes.size();

// Generate random sequential order of nodes.
	vector<int> randomOrder(nSuperNodes);
	for (int i = 0; i < nSuperNodes; i++)
		randomOrder[i] = i;

	for (int i = 0; i < nSuperNodes; i++) {
		int target = R->randInt(nSuperNodes - 1);

		// swap numbers between i and target.
		int tmp = randomOrder[i];
		randomOrder[i] = randomOrder[target];
		randomOrder[target] = tmp;
	}

	int numMoved = 0;

// Move each node to one of its neighbor modules in random sequential order.
	for (int i = 0; i < nSuperNodes; i++) {
		SuperNode& nd = superNodes[randomOrder[i]];	// look at i_th Node of the random sequential order.
		int oldMod = nd.ModIdx();

		unsigned int nModLinks = 0;	// The number of links to/from between the current node and other modules.

		flowmap outFlowToMod;	// <modID, flow> for outFlow...
		flowmap inFlowFromMod;	// <modID, flow> for inFlow...

		// count other modules that are connected to the current node.
		for (link_iterator linkIt = nd.outLinks.begin();
				linkIt != nd.outLinks.end(); linkIt++) {
			int newMod = superNodes[linkIt->first].ModIdx();

			if (outFlowToMod.count(newMod) > 0) {
				outFlowToMod[newMod] += beta * linkIt->second;
			} else {
				outFlowToMod[newMod] = beta * linkIt->second;// initialization of the outFlow of the current newMod.
				inFlowFromMod[newMod] = 0.0;
				nModLinks++;
			}
		}

		for (link_iterator linkIt = nd.inLinks.begin();
				linkIt != nd.inLinks.end(); linkIt++) {
			int newMod = superNodes[linkIt->first].ModIdx();

			if (inFlowFromMod.count(newMod) > 0) {
				inFlowFromMod[newMod] += beta * linkIt->second;
			} else {
				outFlowToMod[newMod] = 0.0;
				inFlowFromMod[newMod] = beta * linkIt->second;
				nModLinks++;
			}
		}

		if (nModLinks != outFlowToMod.size())
			cout << "ALERT: nModLinks != outFlowToMod.size()." << endl;

		// copy node specific values for easy use and efficiency.
		double ndSize = nd.Size();	// p_nd.
		double ndTPWeight = nd.TeleportWeight();	// tau_nd.
		double ndDanglingSize = nd.DanglingSize();

		double oldExitPr1 = modules[oldMod].exitPr;
		double oldSumPr1 = modules[oldMod].sumPr;
		double oldSumDangling1 = modules[oldMod].sumDangling;
		double oldModTPWeight = modules[oldMod].sumTPWeight;

		double additionalTeleportOutFlow = (alpha * ndSize
				+ beta * ndDanglingSize) * (oldModTPWeight - ndTPWeight);
		double additionalTeleportInFlow = (alpha * (oldSumPr1 - ndSize)
				+ beta * (oldSumDangling1 - ndDanglingSize)) * ndTPWeight;

		// For teleportation and danling nodes.
		for (flowmap::iterator it = outFlowToMod.begin();
				it != outFlowToMod.end(); it++) {
			int newMod = it->first;
			if (newMod == oldMod) {
				outFlowToMod[newMod] += additionalTeleportOutFlow;
				inFlowFromMod[newMod] += additionalTeleportInFlow;
			} else {
				outFlowToMod[newMod] += (alpha * ndSize + beta * ndDanglingSize)
						* modules[newMod].sumTPWeight;
				inFlowFromMod[newMod] += (alpha * modules[newMod].sumPr
						+ beta * modules[newMod].sumDangling) * ndTPWeight;
			}
		}

		// Calculate flow to/from own module (default value if no links to own module).
		double outFlowToOldMod = additionalTeleportOutFlow;
		double inFlowFromOldMod = additionalTeleportInFlow;
		if (outFlowToMod.count(oldMod) > 0) {
			outFlowToOldMod = outFlowToMod[oldMod];
			inFlowFromOldMod = inFlowFromMod[oldMod];
		}

		//////////////////// TODO: NEED TO IMPLEMENT THE OPTION TO MOVE TO EMPTY MODULE ////////////////
		if (modules[oldMod].members.size() > ndSize
				&& emptyModules.size() > 0) {
			int emptyTarget = emptyModules.back();
			outFlowToMod[emptyTarget] = 0.0;
			inFlowFromMod[emptyTarget] = 0.0;
			nModLinks++;
		}

		//vector<MoveSummary> moveResults(nModLinks);
		MoveSummary currentResult;
		MoveSummary bestResult;

		double newExitPr1 = oldExitPr1 - nd.ExitPr() + outFlowToOldMod
				+ inFlowFromOldMod;

		bestResult.diffCodeLen = 0.0;// This is the default value, if we can't find diffCodeLen < 0, then don't move the node.

		for (flowmap::iterator it = outFlowToMod.begin();
				it != outFlowToMod.end(); it++) {
			int newMod = it->first;
			double outFlowToNewMod = it->second;
			double inFlowFromNewMod = inFlowFromMod[newMod];

			if (newMod != oldMod) {
				// copy module specific values...
				double oldExitPr2 = modules[newMod].exitPr;
				double oldSumPr2 = modules[newMod].sumPr;

				// Calculate status of current investigated movement of the node nd.
				currentResult.newModule = newMod;
				currentResult.sumPr1 = oldSumPr1 - ndSize;// This should be 0.0, because oldModule will be empty module.
				currentResult.sumPr2 = oldSumPr2 + ndSize;
				currentResult.exitPr1 = newExitPr1;
				currentResult.exitPr2 = oldExitPr2 + nd.ExitPr()
						- outFlowToNewMod - inFlowFromNewMod;

				currentResult.newSumExitPr = sumAllExitPr + newExitPr1
						+ currentResult.exitPr2 - oldExitPr1 - oldExitPr2;

				// Calculate delta_L(M) = L(M)_new - L(M)_old
				double delta_allExit_log_allExit = pLogP(
						currentResult.newSumExitPr) - pLogP(sumAllExitPr);
				double delta_exit_log_exit = pLogP(currentResult.exitPr1)
						+ pLogP(currentResult.exitPr2) - pLogP(oldExitPr1)
						- pLogP(oldExitPr2);
				double delta_stay_log_stay = pLogP(
						currentResult.exitPr1 + currentResult.sumPr1)
						+ pLogP(currentResult.exitPr2 + currentResult.sumPr2)
						- pLogP(oldExitPr1 + oldSumPr1)
						- pLogP(oldExitPr2 + oldSumPr2);

				// delta_L(M) = delta_allExit - 2.0 * delta_exit_log_exit + delta_stay_log_stay.
				currentResult.diffCodeLen = delta_allExit_log_allExit
						- 2.0 * delta_exit_log_exit + delta_stay_log_stay;

				if (currentResult.diffCodeLen < bestResult.diffCodeLen) {
					// we need to update bestResult with currentResult.
					bestResult.diffCodeLen = currentResult.diffCodeLen;
					bestResult.newModule = currentResult.newModule;
					bestResult.sumPr1 = currentResult.sumPr1;
					bestResult.sumPr2 = currentResult.sumPr2;
					bestResult.exitPr1 = currentResult.exitPr1;
					bestResult.exitPr2 = currentResult.exitPr2;
					bestResult.newSumExitPr = currentResult.newSumExitPr;
				}
			}
		}

		// Make best possible move for the current node nd.
		if (bestResult.diffCodeLen < 0.0) {
			// update related to newMod...
			int newMod = bestResult.newModule;
			int spMembers = nd.members.size();

			if (modules[newMod].numMembers == 0) {
				newMod = emptyModules.back();
				emptyModules.pop_back();
				nEmptyMod--;
				nModule++;
			}

			nd.setModIdx(newMod);

			for (int j = 0; j < spMembers; j++) {
				nd.members[j]->setModIdx(newMod);
			}

			modules[newMod].numMembers += spMembers;
			modules[newMod].exitPr = bestResult.exitPr2;
			modules[newMod].sumPr = bestResult.sumPr2;
			modules[newMod].stayPr = bestResult.exitPr2 + bestResult.sumPr2;
			modules[newMod].sumTPWeight += ndTPWeight;

			double spDanglingSize = nd.DanglingSize();
			if (spDanglingSize > 0.0) {
				modules[newMod].sumDangling += spDanglingSize;
				modules[oldMod].sumDangling -= spDanglingSize;
			}

			// update related to the oldMod...
			modules[oldMod].numMembers -= spMembers;
			modules[oldMod].exitPr = bestResult.exitPr1;
			modules[oldMod].sumPr = bestResult.sumPr1;
			modules[oldMod].stayPr = bestResult.exitPr1 + bestResult.sumPr1;
			modules[oldMod].sumTPWeight -= ndTPWeight;

			if (modules[oldMod].numMembers == 0) {
				nEmptyMod++;
				nModule--;
				emptyModules.push_back(oldMod);
			}

			sumAllExitPr = bestResult.newSumExitPr;

			codeLength += bestResult.diffCodeLen;
			numMoved += spMembers;// although we moved a superNode, a superNode is actually a set of nodes.
		}
	}

// the following should be true: modules.size() == nModule.
	if (nodes.size() != nModule + nEmptyMod) {
		cout << "Something wrong!! nodes.size() != nModule + nEmptyMod."
				<< endl;
	}

	return numMoved;
}

void Network::showOutput(int iteration, int prioritizeSPMove, bool inWhile) {
	int rank, size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	for (int i = 0; i < nNode; i++) {
		printf(
				"before distributed output for rank:%d, iteration:%d, nodes:%d, module Id:%d, number of active nodes:%d, nNode:%d,  prioritizeSPMove:%d, insideWhile:%d\n",
				rank, iteration, i, nodes[i].modIdx, activeNodes.size(), nNode,
				prioritizeSPMove, inWhile);
	}
}

/*
 * This function implements a prioritized version of Network::moveSuperNodes() above.
 */

int Network::prioritize_moveSPnodes(double vThresh, int iteration, int first,
		int tag) {

	int nNextActive = 0;
	int SPNodeCountforSending = 0;

	bool incrementElementCount = false;
	set<int> emptyModuleSet;

	int size;
	int rank;

	int numSPNode = superNodes.size();

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (iteration == 0) {
		printf("this is called for rank:%d, where 5\n", rank);
	}

	if (iteration == 0 && tag == 1) {
		for (int it = 0; it < numSPNode; it++) {
			int noOutLinks = superNodes[it].outLinks.size();
			printf(
					"iltering by the process:%d, spNode:%d, number of outlinks:%d, iteration:%d, total number of active nodes:%d\n",
					rank, it, noOutLinks, iteration, activeNodes.size());
			for (int ite = 0; ite < noOutLinks; ite++) {
				printf(
						"ilter by the process:%d, superNodes [%d].outLinks[%d]->otherEnd:%d, weight:%f\n",
						rank, it, ite, superNodes[it].outLinks[ite].first,
						superNodes[it].outLinks[ite].second);
			}
		}
	}

	const int nActive = activeNodes.size();

	int randomGlobalArray[nActive];

	for (int i = 0; i < nActive; i++) {
		randomGlobalArray[i] = activeNodes[i];
	}

	for (int i = 0; i < nActive; i++) {
		int target = R->randInt(nActive - 1);
		// swap numbers between i and target.
		int tmp = randomGlobalArray[i];
		randomGlobalArray[i] = randomGlobalArray[target];
		randomGlobalArray[target] = tmp;
	}

// Generate random sequential order of active nodes.
	const unsigned int intPackSize = 3 * nNode + 3; // 3*nActive+1 index will save emptyModCount and 3*nActive+2 index will save nModuleCount
	const unsigned int doublePackSize = 7 * nNode;
	int intSendPack[intPackSize] = { }; // multiplier 3 for 3 different elements of index, newModules, oldModules and +1 for elementCount
	int intReceivePack[intPackSize] = { };
	double doubleSendPack[doublePackSize] = { }; // multiplier 6 for 6 different elements of diffCodeLen, sumPr1, sumPr2, exitPr1, exitPr2, newSumExitPr
	double doubleReceivePack[doublePackSize] = { };

	double codeLengthReduction = 0.0;
	double currentSumAllExitPr = sumAllExitPr;

	int numMoved = 0;

	SPNodeCountforSending = 0;
	if (iteration == 0) {
		printf("this is called for rank:%d, where 6\n", rank);
	}

	int start, end;
	findAssignedPart(&start, &end, nActive, size, rank);

	int numberOfSPNode = end - start; //number of elements assigned to this processor

	for (int i = 0; i < nActive; i++) {
		printf("r1nk:%d, activeNodes[%d]:%d, randomGlobalArray[%d]:%d\n", rank,
				i, activeNodes[i], i, randomGlobalArray[i]);
	}

	if (iteration == 0) {
		printf("this is called for rank:%d, where 7 number of elements:%d\n",
				rank, numberOfSPNode);
	}

// Move each node to one of its neighbor modules in random sequential order.
	for (int i = start; i < end; i++) {
		SuperNode& nd = superNodes[randomGlobalArray[i]]; // look at i_th Node of the random sequential order.

		bool moduleUpdate = false;

		int oldMod = nd.ModIdx();

		if (iteration == 0) {
			printf(
					"this is called for rank:%d, where 8 number of elements:%d, oldMod:%d\n",
					rank, numberOfSPNode, oldMod);
		}

		unsigned int nModLinks = 0;	// The number of links to/from between the current node and other modules.

		flowmap outFlowToMod;	// <modID, flow> for outFlow...
		flowmap inFlowFromMod;	// <modID, flow> for inFlow...

		// count other modules that are connected to the current node.
		for (link_iterator linkIt = nd.outLinks.begin();
				linkIt != nd.outLinks.end(); linkIt++) {
			int newMod = superNodes[linkIt->first].ModIdx();

			if (outFlowToMod.count(newMod) > 0) {
				outFlowToMod[newMod] += beta * linkIt->second;
			} else {
				outFlowToMod[newMod] = beta * linkIt->second;// initialization of the outFlow of the current newMod.
				inFlowFromMod[newMod] = 0.0;
				nModLinks++;
			}
		}

		for (link_iterator linkIt = nd.inLinks.begin();
				linkIt != nd.inLinks.end(); linkIt++) {
			int newMod = superNodes[linkIt->first].ModIdx();

			if (inFlowFromMod.count(newMod) > 0) {
				inFlowFromMod[newMod] += beta * linkIt->second;
			} else {
				outFlowToMod[newMod] = 0.0;
				inFlowFromMod[newMod] = beta * linkIt->second;
				nModLinks++;
			}
		}

		if (nModLinks != outFlowToMod.size())
			cout << "ALERT: nModLinks != outFlowToMod.size()." << endl;

		// copy node specific values for easy use and efficiency.
		double ndSize = nd.Size();	// p_nd.
		double ndTPWeight = nd.TeleportWeight();	// tau_nd.
		double ndDanglingSize = nd.DanglingSize();

		double oldExitPr1 = modules[oldMod].exitPr;
		double oldSumPr1 = modules[oldMod].sumPr;
		double oldSumDangling1 = modules[oldMod].sumDangling;
		double oldModTPWeight = modules[oldMod].sumTPWeight;

		double additionalTeleportOutFlow = (alpha * ndSize
				+ beta * ndDanglingSize) * (oldModTPWeight - ndTPWeight);
		double additionalTeleportInFlow = (alpha * (oldSumPr1 - ndSize)
				+ beta * (oldSumDangling1 - ndDanglingSize)) * ndTPWeight;

		if (iteration == 0) {
			printf(
					"this is called for rank:%d, where 9 number of elements:%d, oldMod:%d\n",
					rank, numberOfSPNode, oldMod);
		}

		// For teleportation and danling nodes.
		for (flowmap::iterator it = outFlowToMod.begin();
				it != outFlowToMod.end(); it++) {
			int newMod = it->first;
			if (newMod == oldMod) {
				outFlowToMod[newMod] += additionalTeleportOutFlow;
				inFlowFromMod[newMod] += additionalTeleportInFlow;
			} else {
				outFlowToMod[newMod] += (alpha * ndSize + beta * ndDanglingSize)
						* modules[newMod].sumTPWeight;
				inFlowFromMod[newMod] += (alpha * modules[newMod].sumPr
						+ beta * modules[newMod].sumDangling) * ndTPWeight;
			}
		}

		// Calculate flow to/from own module (default value if no links to own module).
		double outFlowToOldMod = additionalTeleportOutFlow;
		double inFlowFromOldMod = additionalTeleportInFlow;
		if (outFlowToMod.count(oldMod) > 0) {
			outFlowToOldMod = outFlowToMod[oldMod];
			inFlowFromOldMod = inFlowFromMod[oldMod];
		}

		MoveSummary currentResult;
		MoveSummary bestResult;

		double newExitPr1 = oldExitPr1 - nd.ExitPr() + outFlowToOldMod
				+ inFlowFromOldMod;

		bestResult.diffCodeLen = 0.0; // This is the default value, if we can't find diffCodeLen < 0, then don't move the node.
		bestResult.newModule = oldMod; // I have to do that for the sole purpose of sending the information of moving to a new module when sending through MPI although may be there is no module change at all

		if (iteration == 0) {
			printf(
					"this is called for rank:%d, where 10 number of elements:%d, oldMod:%d\n",
					rank, numberOfSPNode, oldMod);
		}

		for (flowmap::iterator it = outFlowToMod.begin();
				it != outFlowToMod.end(); it++) {
			int newMod = it->first;
			double outFlowToNewMod = it->second;
			double inFlowFromNewMod = inFlowFromMod[newMod];

			if (iteration == 0) {
				printf("there is the process 3:%d, oldMod:%d, newMod:%d\n",
						rank, oldMod, newMod);
			}

			if (newMod != oldMod) {

				if (iteration == 0) {
					printf("there is the process 1:%d, oldMod:%d, newMod:%d\n",
							rank, oldMod, newMod);
				}
				// copy module specific values...
				double oldExitPr2 = modules[newMod].exitPr;
				double oldSumPr2 = modules[newMod].sumPr;

				// Calculate status of current investigated movement of the node nd.
				currentResult.newModule = newMod;
				currentResult.sumPr1 = oldSumPr1 - ndSize; // This should be 0.0, because oldModule will be empty module.
				currentResult.sumPr2 = oldSumPr2 + ndSize;
				currentResult.exitPr1 = newExitPr1;
				currentResult.exitPr2 = oldExitPr2 + nd.ExitPr()
						- outFlowToNewMod - inFlowFromNewMod;

				currentResult.newSumExitPr = sumAllExitPr + newExitPr1
						+ currentResult.exitPr2 - oldExitPr1 - oldExitPr2;

				// Calculate delta_L(M) = L(M)_new - L(M)_old
				double delta_allExit_log_allExit = pLogP(
						currentResult.newSumExitPr) - pLogP(sumAllExitPr);
				double delta_exit_log_exit = pLogP(currentResult.exitPr1)
						+ pLogP(currentResult.exitPr2) - pLogP(oldExitPr1)
						- pLogP(oldExitPr2);
				double delta_stay_log_stay = pLogP(
						currentResult.exitPr1 + currentResult.sumPr1)
						+ pLogP(currentResult.exitPr2 + currentResult.sumPr2)
						- pLogP(oldExitPr1 + oldSumPr1)
						- pLogP(oldExitPr2 + oldSumPr2);

				// delta_L(M) = delta_allExit - 2.0 * delta_exit_log_exit + delta_stay_log_stay.
				currentResult.diffCodeLen = delta_allExit_log_allExit
						- 2.0 * delta_exit_log_exit + delta_stay_log_stay;

				if (iteration == 0) {
					printf(
							"this is called for rank:%d, where 11 number of elements:%d, oldMod:%d\n",
							rank, numberOfSPNode, oldMod);
				}

				if (currentResult.diffCodeLen < bestResult.diffCodeLen) {

					moduleUpdate = true;

					intSendPack[SPNodeCountforSending] = i;
					// we need to update bestResult with currentResult.

					//save old module id
					intSendPack[1 * nNode + SPNodeCountforSending] = oldMod;
					//save new module id
					bestResult.newModule = intSendPack[2 * nNode
							+ SPNodeCountforSending] = currentResult.newModule;

					bestResult.diffCodeLen =
							doubleSendPack[SPNodeCountforSending] =
									currentResult.diffCodeLen;
					//bestResult.diffCodeLen = currentResult.diffCodeLen;
					//bestResult.newModule = currentResult.newModule;
					bestResult.sumPr1 = doubleSendPack[1 * nActive
							+ SPNodeCountforSending] = currentResult.sumPr1;
					//bestResult.sumPr1 = currentResult.sumPr1;
					bestResult.sumPr2 = doubleSendPack[2 * nActive
							+ SPNodeCountforSending] = currentResult.sumPr2;
					//bestResult.sumPr2 = currentResult.sumPr2;
					bestResult.exitPr1 = doubleSendPack[3 * nActive
							+ SPNodeCountforSending] = currentResult.exitPr1;
					//bestResult.exitPr1 = currentResult.exitPr1;
					bestResult.exitPr2 = doubleSendPack[4 * nActive
							+ SPNodeCountforSending] = currentResult.exitPr2;
					//bestResult.exitPr2 = currentResult.exitPr2;
					bestResult.newSumExitPr = currentResult.newSumExitPr;
					//bestResult.newSumExitPr = currentResult.newSumExitPr;

				}
			} else { //this element/vertex is not moving to a new module, so both the newModule and the oldModule are same for this element

				moduleUpdate = true;

				intSendPack[SPNodeCountforSending] = i;
				intSendPack[1 * nNode + SPNodeCountforSending] = oldMod;
				intSendPack[2 * nNode + SPNodeCountforSending] =
						bestResult.newModule;
				if (iteration == 0 && tag == 1) {
					printf("there is the process:%d, oldMod:%d, newMod:%d\n",
							rank, oldMod, bestResult.newModule);
				}
			}
		}

		// Make best possible move for the current node nd.
		//if (bestResult.diffCodeLen < 0.0) {

		if (bestResult.diffCodeLen < vThresh) {
			// update related to newMod...

			int newMod = bestResult.newModule;
			int spMembers = nd.members.size();

			if (oldMod != newMod) {

				nd.setModIdx(newMod);

				for (int j = 0; j < spMembers; j++) {
					nd.members[j]->setModIdx(newMod);
				}

				//modules[newMod].numMembers += spMembers;
				modules[newMod].exitPr = bestResult.exitPr2;
				modules[newMod].sumPr = bestResult.sumPr2;
				modules[newMod].stayPr = bestResult.exitPr2 + bestResult.sumPr2;
				modules[newMod].sumTPWeight += ndTPWeight;

				double spDanglingSize = nd.DanglingSize();
				if (spDanglingSize > 0.0) {
					modules[newMod].sumDangling += spDanglingSize;
					modules[oldMod].sumDangling -= spDanglingSize;
				}

				// update related to the oldMod...
				//modules[oldMod].numMembers -= spMembers;
				modules[oldMod].exitPr = bestResult.exitPr1;
				modules[oldMod].sumPr = bestResult.sumPr1;
				modules[oldMod].stayPr = bestResult.exitPr1 + bestResult.sumPr1;
				modules[oldMod].sumTPWeight -= ndTPWeight;

				if (iteration == 0) {
					printf("this is called for rank:%d\n", rank);
				}
				sumAllExitPr = bestResult.newSumExitPr;

				codeLengthReduction += bestResult.diffCodeLen;
				//codeLength += bestResult.diffCodeLen;
				numMoved += spMembers;

				// update activeNodes and isActives vectors.
				// We have to add the following nodes in activeNodes: neighbors, members in oldMod & newMod.
				for (link_iterator linkIt = nd.outLinks.begin();
						linkIt != nd.outLinks.end(); linkIt++)
					isActives[linkIt->first] = 1;	// set as an active nodes.

				for (link_iterator linkIt = nd.inLinks.begin();
						linkIt != nd.inLinks.end(); linkIt++)
					isActives[linkIt->first] = 1;	// set as an active nodes.*/
			}
		}

		if (moduleUpdate == false) {
			intSendPack[SPNodeCountforSending] = i;
			intSendPack[1 * nNode + SPNodeCountforSending] = oldMod;
			intSendPack[2 * nNode + SPNodeCountforSending] = oldMod;
		}

		SPNodeCountforSending++;
	}

	if (numberOfSPNode != SPNodeCountforSending) {
		if (iteration == 0) {
			printf(
					"something is not right in prioritize_movespnodes, for rank:%d, numberOfSPNode:%d does not match SPNodeCountforSending:%d\n",
					rank, numberOfSPNode, SPNodeCountforSending);
		}
	}

	intSendPack[intPackSize - 1] = SPNodeCountforSending;
	doubleSendPack[doublePackSize - 1] = sumAllExitPr - currentSumAllExitPr;
	codeLength += codeLengthReduction;
	doubleSendPack[doublePackSize - 2] = codeLengthReduction;

	if (iteration == 0) {
		printf("this is called for rank:%d, where 1\n", rank);
	}

//good point: if someone has no data to send, find some way to indicate that
	for (int prId = 0; prId < size; prId++) {
		if (prId != rank) {
			MPI_Send(intSendPack, intPackSize, MPI_INT, prId, iteration,
			MPI_COMM_WORLD);
			MPI_Send(doubleSendPack, doublePackSize, MPI_DOUBLE, prId,
					iteration,
					MPI_COMM_WORLD);
		}
	}

	for (int sId = 0; sId < size; sId++) {//my going through loop from sId 0->size I am imposing an ordering which is not good for MPI Performance
		if (sId != rank) {
			MPI_Recv(intReceivePack, intPackSize, MPI_INT, sId, iteration,
			MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
			MPI_Recv(doubleReceivePack, doublePackSize, MPI_DOUBLE, sId,
					iteration,
					MPI_COMM_WORLD, MPI_STATUSES_IGNORE);

			int totalSPNodefromSender = intReceivePack[intPackSize - 1];

			if (iteration == 0 && tag == 1) {
				printf(
						"woh, rank:%d, sender Process:%d, number of elements sent:%d\n",
						rank, sId, totalSPNodefromSender);
			}
			sumAllExitPr += doubleReceivePack[doublePackSize - 1];
			codeLength += doubleReceivePack[doublePackSize - 2];

			if (iteration == 0) {
				printf("this is called for rank:%d, where 2\n", rank);
			}

			for (int z = 0; z < totalSPNodefromSender; z++) {

				int nodeIndex = intReceivePack[z];
				SuperNode& nod = superNodes[randomGlobalArray[nodeIndex]];// I didn't need to send the indices from the processes because every process has already this information from counts and displs
				int spMembers = nod.members.size();

				int oldModule = intReceivePack[1 * nNode + z];
				int newModule = intReceivePack[2 * nNode + z];

				if (iteration == 0) {
					printf("this is called for rank:%d, where 3\n", rank);
				}

				if (oldModule != newModule) {

					numMoved += spMembers;

					if (iteration == 0) {
						printf("for process:%d, 2 number Moved:%d\n", rank,
								numMoved);

					}

					nod.setModIdx(newModule);

					for (int j = 0; j < spMembers; j++) {
						nod.members[j]->setModIdx(newModule);
						if (iteration == 0 && tag == 1) {
							printf(
									"exposing rank:%d, spMembers:%d, nd.members[%d]:%d, nd.members[%d].module:%d, oldMod:%d, newMod:%d\n",
									rank, spMembers, j, nod.members.at(j)->id,
									j, nod.members.at(j)->modIdx, oldModule,
									newModule);
						}
					}

					modules[newModule].exitPr = doubleReceivePack[4 * nActive
							+ z];
					modules[newModule].sumPr =
							doubleReceivePack[2 * nActive + z];
					modules[newModule].stayPr = modules[newModule].exitPr
							+ modules[newModule].sumPr;
					modules[newModule].sumTPWeight += nod.TeleportWeight();

					double spDanglingSize = nod.DanglingSize();

					if (spDanglingSize > 0.0) {
						modules[newModule].sumDangling += spDanglingSize;
						modules[oldModule].sumDangling -= spDanglingSize;
					}

					//modules[oldModule].numMembers -= spMembers;

					modules[oldModule].exitPr = doubleReceivePack[3 * nActive
							+ z];
					modules[oldModule].sumPr =
							doubleReceivePack[1 * nActive + z];
					modules[oldModule].stayPr = modules[oldModule].exitPr
							+ modules[oldModule].sumPr;
					modules[oldModule].sumTPWeight -= nod.TeleportWeight();

					// update activeNodes and isActives vectors.
					// We have to add the following nodes in activeNodes: neighbors, members in oldMod & newMod.
					for (link_iterator linkIt = nod.outLinks.begin();
							linkIt != nod.outLinks.end(); linkIt++)
						isActives[linkIt->first] = 1;// set as an active nodes.

					for (link_iterator linkIt = nod.inLinks.begin();
							linkIt != nod.inLinks.end(); linkIt++)
						isActives[linkIt->first] = 1;// set as an active nodes.
				}
			}
		}
	}

	vector<int>().swap(activeNodes);
	for (int i = 0; i < isActives.size(); i++) {
		if (isActives[i] == 1) {
			activeNodes.push_back(i);
			isActives[i] = 0;	// reset the flag of isActives[i].
		}
	}

	if (iteration == 2 && tag == 1) {
		for (int i = 0; i < nNode; i++) {
			printf("distributed output for rank:%d, nodes[%d], module Id:%d\n",
					rank, i, nodes[i].modIdx);
		}
		double q = calculateModularityScore();
		printf("distributed output for rank:%d, modularity:%f\n", rank, q);
	}

//delete[] randomGlobalArray;
//printf("numberMoved:%d\n", numMoved);
	return numMoved;
}

double Network::calculateModularityScore() {
	const int D = nNode;
	double eij[D][D] = { 0.0 };
	double a[D] = { 0.0 };
	double Q = 0.0;

	for (int i = 0; i < nNode; i++) {
		for (int j = 0; j < nNode; j++) {
			eij[i][j] = 0.0;
		}
		a[i] = 0.0;
	}

	printf("before calculating modularity score\n");
	/*	for (int i = 0; i < nNode; i++) {
	 for (int j = 0; j < nNode; j++) {
	 printf("%.1f ", eij[i][j]);
	 }
	 printf("\n");
	 }*/

	printf("calculating modularity score\n");
	for (int i = 0; i < nNode; i++) {
		Node& nd = nodes[i];
		for (link_iterator linkIt = nd.outLinks.begin();
				linkIt != nd.outLinks.end(); linkIt++) {
			Node& otherNd = nodes[linkIt->first];
			eij[nd.modIdx][otherNd.modIdx]++;
		}
		for (link_iterator linkIt = nd.inLinks.begin();
				linkIt != nd.inLinks.end(); linkIt++) {
			Node& otherNd = nodes[linkIt->first];
			eij[otherNd.modIdx][nd.modIdx]++; // for the inLinks it is not eij rather it is eji i.e. eji represents the fraction of edges from module j to module i
		}
	}

	for (int i = 0; i < nNode; i++) {
		for (int j = 0; j < nNode; j++) {
			eij[i][j] /= (2 * nEdge);
		}
	}

	for (int i = 0; i < nNode; i++) {
		for (int j = 0; j < nNode; j++) {
			a[i] += eij[i][j];
		}
	}

	for (int i = 0; i < nNode; i++) {
		Q += eij[i][i] - pow(a[i], 2);
	}

	/*	for (int i = 0; i < nNode; i++) {
	 for (int j = 0; j < nNode; j++) {
	 printf("%.1f ", eij[i][j]);
	 }
	 printf("\n");
	 }*/

	return Q;
}

/*
 *	+ This is a parallel implementation of Network::moveSuperNodes() function via OpenMP.
 *		We will loose the criteria for simple implementation and it may help to avoid local optima problem.
 */

int Network::parallelMoveSuperNodes(int numTh, double& tSequential) {

	struct timeval tStart, tEnd;

	gettimeofday(&tStart, NULL);

	int nSuperNodes = superNodes.size();

// Generate random sequential order of nodes.
	vector<int> randomOrder(nSuperNodes);
	for (int i = 0; i < nSuperNodes; i++)
		randomOrder[i] = i;

	for (int i = 0; i < nSuperNodes; i++) {
		int target = R->randInt(nSuperNodes - 1);

		// swap numbers between i and target.
		int tmp = randomOrder[i];
		randomOrder[i] = randomOrder[target];
		randomOrder[target] = tmp;
	}

	int numMoved = 0;

	omp_set_num_threads(numTh);

	const int emptyTarget = nNode + 1; // This will be an indicator of moving to emptyModule.

	gettimeofday(&tEnd, NULL);

	tSequential += elapsedTimeInSec(tStart, tEnd);

#pragma omp parallel for
// Move each node to one of its neighbor modules in random sequential order.
	for (int i = 0; i < nSuperNodes; i++) {
		SuperNode& nd = superNodes[randomOrder[i]];	// look at i_th Node of the random sequential order.
		int oldMod = nd.ModIdx();

		unsigned int nModLinks = 0;	// The number of links to/from between the current node and other modules.

		flowmap outFlowToMod;	// <modID, flow> for outFlow...
		flowmap inFlowFromMod;	// <modID, flow> for inFlow...

		// count other modules that are connected to the current node.
		for (link_iterator linkIt = nd.outLinks.begin();
				linkIt != nd.outLinks.end(); linkIt++) {
			int newMod = superNodes[linkIt->first].ModIdx();

			if (outFlowToMod.count(newMod) > 0) {
				outFlowToMod[newMod] += beta * linkIt->second;
			} else {
				outFlowToMod[newMod] = beta * linkIt->second;// initialization of the outFlow of the current newMod.
				inFlowFromMod[newMod] = 0.0;
				nModLinks++;
			}
		}

		for (link_iterator linkIt = nd.inLinks.begin();
				linkIt != nd.inLinks.end(); linkIt++) {
			int newMod = superNodes[linkIt->first].ModIdx();

			if (inFlowFromMod.count(newMod) > 0) {
				inFlowFromMod[newMod] += beta * linkIt->second;
			} else {
				outFlowToMod[newMod] = 0.0;
				inFlowFromMod[newMod] = beta * linkIt->second;
				nModLinks++;
			}
		}

		if (nModLinks != outFlowToMod.size())
			cout << "ALERT: nModLinks != outFlowToMod.size()." << endl;

		// copy node specific values for easy use and efficiency.
		double ndSize = nd.Size();	// p_nd.
		double ndTPWeight = nd.TeleportWeight();	// tau_nd.
		double ndDanglingSize = nd.DanglingSize();

		double oldExitPr1 = modules[oldMod].exitPr;
		double oldSumPr1 = modules[oldMod].sumPr;
		double oldSumDangling1 = modules[oldMod].sumDangling;
		double oldModTPWeight = modules[oldMod].sumTPWeight;

		double additionalTeleportOutFlow = (alpha * ndSize
				+ beta * ndDanglingSize) * (oldModTPWeight - ndTPWeight);
		double additionalTeleportInFlow = (alpha * (oldSumPr1 - ndSize)
				+ beta * (oldSumDangling1 - ndDanglingSize)) * ndTPWeight;

		// For teleportation and danling nodes.
		for (flowmap::iterator it = outFlowToMod.begin();
				it != outFlowToMod.end(); it++) {
			int newMod = it->first;
			if (newMod == oldMod) {
				outFlowToMod[newMod] += additionalTeleportOutFlow;
				inFlowFromMod[newMod] += additionalTeleportInFlow;
			} else {
				outFlowToMod[newMod] += (alpha * ndSize + beta * ndDanglingSize)
						* modules[newMod].sumTPWeight;
				inFlowFromMod[newMod] += (alpha * modules[newMod].sumPr
						+ beta * modules[newMod].sumDangling) * ndTPWeight;
			}
		}

		// Calculate flow to/from own module (default value if no links to own module).
		double outFlowToOldMod = additionalTeleportOutFlow;
		double inFlowFromOldMod = additionalTeleportInFlow;
		if (outFlowToMod.count(oldMod) > 0) {
			outFlowToOldMod = outFlowToMod[oldMod];
			inFlowFromOldMod = inFlowFromMod[oldMod];
		}

		//////////////////// TODO: NEED TO IMPLEMENT THE OPTION TO MOVE TO EMPTY MODULE ////////////////
		if (modules[oldMod].sumPr > ndSize && emptyModules.size() > 0) {
			outFlowToMod[emptyTarget] = 0.0;
			inFlowFromMod[emptyTarget] = 0.0;
			nModLinks++;
		}

		MoveSummary currentResult;
		MoveSummary bestResult;

		double newExitPr1 = oldExitPr1 - nd.ExitPr() + outFlowToOldMod
				+ inFlowFromOldMod;

		bestResult.diffCodeLen = 0.0;// This is the default value, if we can't find diffCodeLen < 0, then don't move the node.

		for (flowmap::iterator it = outFlowToMod.begin();
				it != outFlowToMod.end(); it++) {
			int newMod = it->first;
			double outFlowToNewMod = it->second;
			double inFlowFromNewMod = inFlowFromMod[newMod];

			if (newMod != oldMod) {

				// copy module specific values...
				double oldExitPr2 = modules[newMod].exitPr;
				double oldSumPr2 = modules[newMod].sumPr;

				// Calculate status of current investigated movement of the node nd.
				currentResult.newModule = newMod;
				currentResult.sumPr1 = oldSumPr1 - ndSize;// This should be 0.0, because oldModule will be empty module.
				currentResult.sumPr2 = oldSumPr2 + ndSize;
				currentResult.exitPr1 = newExitPr1;
				currentResult.exitPr2 = oldExitPr2 + nd.ExitPr()
						- outFlowToNewMod - inFlowFromNewMod;

				currentResult.newSumExitPr = sumAllExitPr + newExitPr1
						+ currentResult.exitPr2 - oldExitPr1 - oldExitPr2;

				// Calculate delta_L(M) = L(M)_new - L(M)_old
				double delta_allExit_log_allExit = pLogP(
						currentResult.newSumExitPr) - pLogP(sumAllExitPr);
				double delta_exit_log_exit = pLogP(currentResult.exitPr1)
						+ pLogP(currentResult.exitPr2) - pLogP(oldExitPr1)
						- pLogP(oldExitPr2);
				double delta_stay_log_stay = pLogP(
						currentResult.exitPr1 + currentResult.sumPr1)
						+ pLogP(currentResult.exitPr2 + currentResult.sumPr2)
						- pLogP(oldExitPr1 + oldSumPr1)
						- pLogP(oldExitPr2 + oldSumPr2);

				// delta_L(M) = delta_allExit - 2.0 * delta_exit_log_exit + delta_stay_log_stay.
				currentResult.diffCodeLen = delta_allExit_log_allExit
						- 2.0 * delta_exit_log_exit + delta_stay_log_stay;

				if (currentResult.diffCodeLen < bestResult.diffCodeLen) {
					// we need to update bestResult with currentResult.
					bestResult.diffCodeLen = currentResult.diffCodeLen;
					bestResult.newModule = currentResult.newModule;
				}
			}

		}

		// store the best possilbe movement information if necessary.
		if (bestResult.diffCodeLen < 0.0) {

			bool isEmptyTarget = false;
			bool validMove = true;// This will indicate the validity of the decided movement.
			int newMod = bestResult.newModule;
			int spMembers = nd.members.size();

#pragma omp critical (moveUpdate)
			{
				gettimeofday(&tStart, NULL);

				// if newMod == emptyTarget, it indicates moves to empty module.
				if ((nEmptyMod > 0) && (newMod == emptyTarget)
						&& (modules[oldMod].numMembers > 1)) {
					newMod = emptyModules.back();
					isEmptyTarget = true;
				} else if (newMod == emptyTarget) {
					validMove = false;
				} else if (modules[newMod].numMembers == 0) {
					// This is the case that the algorithm thought there are some nodes in the new module since newMod != emptyTarget,
					// but the nodes are all moved to other modules so there are no nodes in there.
					// Thus, we don't know whether generating a new module will be better or not.
					// Therefore, don't move this option.
					//continue;
					validMove = false;
				}

				MoveSummary moveResult;

				if (validMove) {

					////////////////////////////////////////////////////////////////////
					/// THIS IS THE PART FOR EXAMINING QUALITY IMPROVEMENT OR NOT... ///
					////////////////////////////////////////////////////////////////////

					outFlowToOldMod = 0.0;
					double outFlowToNewMod = 0.0;
					inFlowFromOldMod = 0.0;
					double inFlowFromNewMod = 0.0;

					if (!isEmptyTarget) {
						for (link_iterator linkIt = nd.outLinks.begin();
								linkIt != nd.outLinks.end(); linkIt++) {
							int toMod = superNodes[linkIt->first].ModIdx();

							if (toMod == oldMod) {
								outFlowToOldMod += beta * linkIt->second;
							} else if (toMod == newMod) {
								outFlowToNewMod += beta * linkIt->second;
							}
						}

						for (link_iterator linkIt = nd.inLinks.begin();
								linkIt != nd.inLinks.end(); linkIt++) {
							int fromMod = superNodes[linkIt->first].ModIdx();

							if (fromMod == oldMod) {
								inFlowFromOldMod += beta * linkIt->second;
							} else if (fromMod == newMod) {
								inFlowFromNewMod += beta * linkIt->second;
							}
						}
					} else {
						for (link_iterator linkIt = nd.outLinks.begin();
								linkIt != nd.outLinks.end(); linkIt++) {
							int toMod = superNodes[linkIt->first].ModIdx();
							if (toMod == oldMod) {
								outFlowToOldMod += beta * linkIt->second;
							}
						}

						for (link_iterator linkIt = nd.inLinks.begin();
								linkIt != nd.inLinks.end(); linkIt++) {
							int fromMod = superNodes[linkIt->first].ModIdx();
							if (fromMod == oldMod) {
								inFlowFromOldMod += beta * linkIt->second;
							}
						}
					}

					oldExitPr1 = modules[oldMod].exitPr;
					oldSumPr1 = modules[oldMod].sumPr;
					oldSumDangling1 = modules[oldMod].sumDangling;
					oldModTPWeight = modules[oldMod].sumTPWeight;

					// For teleportation and danling nodes.
					outFlowToOldMod += (alpha * ndSize + beta * ndDanglingSize)
							* (oldModTPWeight - ndTPWeight);
					inFlowFromOldMod += (alpha * (oldSumPr1 - ndSize)
							+ beta * (oldSumDangling1 - ndDanglingSize))
							* ndTPWeight;
					outFlowToNewMod += (alpha * ndSize + beta * ndDanglingSize)
							* modules[newMod].sumTPWeight;
					inFlowFromNewMod += (alpha * modules[newMod].sumPr
							+ beta * modules[newMod].sumDangling) * ndTPWeight;

					if (isEmptyTarget) {
						outFlowToNewMod = 0.0;
						inFlowFromNewMod = 0.0;
					}

					moveResult.exitPr1 = oldExitPr1 - nd.ExitPr()
							+ outFlowToOldMod + inFlowFromOldMod;

					// copy module specific values...
					double oldExitPr2 = modules[newMod].exitPr;
					double oldSumPr2 = modules[newMod].sumPr;

					// Calculate status of current investigated movement of the node nd.
					moveResult.newModule = newMod;
					moveResult.sumPr1 = oldSumPr1 - ndSize;	// This should be 0.0, because oldModule will be empty module.
					moveResult.sumPr2 = oldSumPr2 + ndSize;
					moveResult.exitPr2 = oldExitPr2 + nd.ExitPr()
							- outFlowToNewMod - inFlowFromNewMod;

					moveResult.newSumExitPr = sumAllExitPr + moveResult.exitPr1
							+ moveResult.exitPr2 - oldExitPr1 - oldExitPr2;

					// Calculate delta_L(M) = L(M)_new - L(M)_old
					double delta_allExit_log_allExit = pLogP(
							moveResult.newSumExitPr) - pLogP(sumAllExitPr);
					double delta_exit_log_exit = pLogP(moveResult.exitPr1)
							+ pLogP(moveResult.exitPr2) - pLogP(oldExitPr1)
							- pLogP(oldExitPr2);
					double delta_stay_log_stay = pLogP(
							moveResult.exitPr1 + moveResult.sumPr1)
							+ pLogP(moveResult.exitPr2 + moveResult.sumPr2)
							- pLogP(oldExitPr1 + oldSumPr1)
							- pLogP(oldExitPr2 + oldSumPr2);

					// delta_L(M) = delta_allExit - 2.0 * delta_exit_log_exit + delta_stay_log_stay.
					moveResult.diffCodeLen = delta_allExit_log_allExit
							- 2.0 * delta_exit_log_exit + delta_stay_log_stay;

					////////////////////////////////////
					////// THE END OF EXAMINATION //////
					////////////////////////////////////

					if (isEmptyTarget) {
						emptyModules.pop_back();

						nEmptyMod--;
						nModule++;
					}

					nd.setModIdx(newMod);
					for (int k = 0; k < spMembers; k++) {
						nd.members[k]->setModIdx(newMod);
					}

					modules[newMod].numMembers += spMembers;
					modules[newMod].exitPr = moveResult.exitPr2;
					modules[newMod].sumPr = moveResult.sumPr2;
					modules[newMod].stayPr = moveResult.exitPr2
							+ moveResult.sumPr2;
					modules[newMod].sumTPWeight += ndTPWeight;

					if (ndDanglingSize > 0.0) {
						modules[newMod].sumDangling += ndDanglingSize;
						modules[oldMod].sumDangling -= ndDanglingSize;
					}

					// update related to the oldMod...
					modules[oldMod].numMembers -= spMembers;
					modules[oldMod].exitPr = moveResult.exitPr1;
					modules[oldMod].sumPr = moveResult.sumPr1;
					modules[oldMod].stayPr = moveResult.exitPr1
							+ moveResult.sumPr1;
					modules[oldMod].sumTPWeight -= ndTPWeight;

					if (modules[oldMod].numMembers == 0) {
						nEmptyMod++;
						nModule--;
						emptyModules.push_back(oldMod);
					}

					sumAllExitPr = moveResult.newSumExitPr;
					codeLength += moveResult.diffCodeLen;

					numMoved += spMembers;
				}

				gettimeofday(&tEnd, NULL);

				tSequential += elapsedTimeInSec(tStart, tEnd);
			}	// END critical

		}

	}	// END parallel for

	return numMoved;
}

/*
 *	+ This is a parallel implementation of Network::prioritize_moveSPnodes() function via OpenMP.
 *		We will loose the criteria for simple implementation and it may help to avoid local optima problem.
 */

int Network::prioritize_parallelMoveSPnodes(int numTh, double& tSequential,
		double vThresh) {

	struct timeval tStart, tEnd;

	gettimeofday(&tStart, NULL);

	int nActive = activeNodes.size();
	int nNextActive = 0;

	cout << "inside movesuper node" << endl;
// Generate random sequential order of nodes.
	vector<int> randomOrder(nActive);
	for (int i = 0; i < nActive; i++)
		randomOrder[i] = activeNodes[i];

	for (int i = 0; i < nActive; i++) {
		int target = R->randInt(nActive - 1);

		// swap numbers between i and target.
		int tmp = randomOrder[i];
		randomOrder[i] = randomOrder[target];
		randomOrder[target] = tmp;
	}

// Now randomOrder vector already had all activeNodes info,
// so we can reset activeNodes vector for adding for new active nodes.
	vector<int>().swap(activeNodes);

	int numMoved = 0;

	omp_set_num_threads(numTh);

	const int emptyTarget = nNode + 1;// This will be an indicator of moving to emptyModule.

	gettimeofday(&tEnd, NULL);

	tSequential += elapsedTimeInSec(tStart, tEnd);

#pragma omp parallel for
// Move each node to one of its neighbor modules in random sequential order.
	for (int i = 0; i < nActive; i++) {
		SuperNode& nd = superNodes[randomOrder[i]];	// look at i_th Node of the random sequential order.
		int oldMod = nd.ModIdx();

		unsigned int nModLinks = 0;	// The number of links to/from between the current node and other modules.

		flowmap outFlowToMod;	// <modID, flow> for outFlow...
		flowmap inFlowFromMod;	// <modID, flow> for inFlow...

		// count other modules that are connected to the current node.
		for (link_iterator linkIt = nd.outLinks.begin();
				linkIt != nd.outLinks.end(); linkIt++) {
			int newMod = superNodes[linkIt->first].ModIdx();

			if (outFlowToMod.count(newMod) > 0) {
				outFlowToMod[newMod] += beta * linkIt->second;
			} else {
				outFlowToMod[newMod] = beta * linkIt->second;// initialization of the outFlow of the current newMod.
				inFlowFromMod[newMod] = 0.0;
				nModLinks++;
			}
		}

		for (link_iterator linkIt = nd.inLinks.begin();
				linkIt != nd.inLinks.end(); linkIt++) {
			int newMod = superNodes[linkIt->first].ModIdx();

			if (inFlowFromMod.count(newMod) > 0) {
				inFlowFromMod[newMod] += beta * linkIt->second;
			} else {
				outFlowToMod[newMod] = 0.0;
				inFlowFromMod[newMod] = beta * linkIt->second;
				nModLinks++;
			}
		}

		if (nModLinks != outFlowToMod.size())
			cout << "ALERT: nModLinks != outFlowToMod.size()." << endl;

		// copy node specific values for easy use and efficiency.
		double ndSize = nd.Size();	// p_nd.
		double ndTPWeight = nd.TeleportWeight();	// tau_nd.
		double ndDanglingSize = nd.DanglingSize();

		double oldExitPr1 = modules[oldMod].exitPr;
		double oldSumPr1 = modules[oldMod].sumPr;
		double oldSumDangling1 = modules[oldMod].sumDangling;
		double oldModTPWeight = modules[oldMod].sumTPWeight;

		double additionalTeleportOutFlow = (alpha * ndSize
				+ beta * ndDanglingSize) * (oldModTPWeight - ndTPWeight);
		double additionalTeleportInFlow = (alpha * (oldSumPr1 - ndSize)
				+ beta * (oldSumDangling1 - ndDanglingSize)) * ndTPWeight;

		// For teleportation and danling nodes.
		for (flowmap::iterator it = outFlowToMod.begin();
				it != outFlowToMod.end(); it++) {
			int newMod = it->first;
			if (newMod == oldMod) {
				outFlowToMod[newMod] += additionalTeleportOutFlow;
				inFlowFromMod[newMod] += additionalTeleportInFlow;
			} else {
				outFlowToMod[newMod] += (alpha * ndSize + beta * ndDanglingSize)
						* modules[newMod].sumTPWeight;
				inFlowFromMod[newMod] += (alpha * modules[newMod].sumPr
						+ beta * modules[newMod].sumDangling) * ndTPWeight;
			}
		}

		// Calculate flow to/from own module (default value if no links to own module).
		double outFlowToOldMod = additionalTeleportOutFlow;
		double inFlowFromOldMod = additionalTeleportInFlow;
		if (outFlowToMod.count(oldMod) > 0) {
			outFlowToOldMod = outFlowToMod[oldMod];
			inFlowFromOldMod = inFlowFromMod[oldMod];
		}

		//////////////////// TODO: NEED TO IMPLEMENT THE OPTION TO MOVE TO EMPTY MODULE ////////////////
		if (modules[oldMod].sumPr > ndSize && emptyModules.size() > 0) {
			outFlowToMod[emptyTarget] = 0.0;
			inFlowFromMod[emptyTarget] = 0.0;
			nModLinks++;
		}

		MoveSummary currentResult;
		MoveSummary bestResult;

		double newExitPr1 = oldExitPr1 - nd.ExitPr() + outFlowToOldMod
				+ inFlowFromOldMod;

		bestResult.diffCodeLen = 0.0;// This is the default value, if we can't find diffCodeLen < 0, then don't move the node.

		for (flowmap::iterator it = outFlowToMod.begin();
				it != outFlowToMod.end(); it++) {
			int newMod = it->first;
			double outFlowToNewMod = it->second;
			double inFlowFromNewMod = inFlowFromMod[newMod];

			if (newMod != oldMod) {

				// copy module specific values...
				double oldExitPr2 = modules[newMod].exitPr;
				double oldSumPr2 = modules[newMod].sumPr;

				// Calculate status of current investigated movement of the node nd.
				currentResult.newModule = newMod;
				currentResult.sumPr1 = oldSumPr1 - ndSize;// This should be 0.0, because oldModule will be empty module.
				currentResult.sumPr2 = oldSumPr2 + ndSize;
				currentResult.exitPr1 = newExitPr1;
				currentResult.exitPr2 = oldExitPr2 + nd.ExitPr()
						- outFlowToNewMod - inFlowFromNewMod;

				currentResult.newSumExitPr = sumAllExitPr + newExitPr1
						+ currentResult.exitPr2 - oldExitPr1 - oldExitPr2;

				// Calculate delta_L(M) = L(M)_new - L(M)_old
				double delta_allExit_log_allExit = pLogP(
						currentResult.newSumExitPr) - pLogP(sumAllExitPr);
				double delta_exit_log_exit = pLogP(currentResult.exitPr1)
						+ pLogP(currentResult.exitPr2) - pLogP(oldExitPr1)
						- pLogP(oldExitPr2);
				double delta_stay_log_stay = pLogP(
						currentResult.exitPr1 + currentResult.sumPr1)
						+ pLogP(currentResult.exitPr2 + currentResult.sumPr2)
						- pLogP(oldExitPr1 + oldSumPr1)
						- pLogP(oldExitPr2 + oldSumPr2);

				// delta_L(M) = delta_allExit - 2.0 * delta_exit_log_exit + delta_stay_log_stay.
				currentResult.diffCodeLen = delta_allExit_log_allExit
						- 2.0 * delta_exit_log_exit + delta_stay_log_stay;

				if (currentResult.diffCodeLen < bestResult.diffCodeLen) {
					// we need to update bestResult with currentResult.
					bestResult.diffCodeLen = currentResult.diffCodeLen;
					bestResult.newModule = currentResult.newModule;
				}
			}

		}

		// store the best possilbe movement information if necessary.
		//if (bestResult.diffCodeLen < 0.0) {
		if (bestResult.diffCodeLen < vThresh) {

			bool isEmptyTarget = false;
			bool validMove = true;// This will indicate the validity of the decided movement.
			int newMod = bestResult.newModule;
			int spMembers = nd.members.size();

#pragma omp critical (moveUpdate)
			{
				gettimeofday(&tStart, NULL);

				// if newMod == emptyTarget, it indicates moves to empty module.
				if ((nEmptyMod > 0) && (newMod == emptyTarget)
						&& (modules[oldMod].numMembers > 1)) {
					newMod = emptyModules.back();
					isEmptyTarget = true;
				} else if (newMod == emptyTarget) {
					validMove = false;
				} else if (modules[newMod].numMembers == 0) {
					// This is the case that the algorithm thought there are some nodes in the new module since newMod != emptyTarget,
					// but the nodes are all moved to other modules so there are no nodes in there.
					// Thus, we don't know whether generating a new module will be better or not.
					// Therefore, don't move this option.
					//continue;
					validMove = false;
				}

				MoveSummary moveResult;

				if (validMove) {

					////////////////////////////////////////////////////////////////////
					/// THIS IS THE PART FOR EXAMINING QUALITY IMPROVEMENT OR NOT... ///
					////////////////////////////////////////////////////////////////////

					outFlowToOldMod = 0.0;
					double outFlowToNewMod = 0.0;
					inFlowFromOldMod = 0.0;
					double inFlowFromNewMod = 0.0;

					if (!isEmptyTarget) {
						for (link_iterator linkIt = nd.outLinks.begin();
								linkIt != nd.outLinks.end(); linkIt++) {
							int toMod = superNodes[linkIt->first].ModIdx();

							if (toMod == oldMod) {
								outFlowToOldMod += beta * linkIt->second;
							} else if (toMod == newMod) {
								outFlowToNewMod += beta * linkIt->second;
							}
						}

						for (link_iterator linkIt = nd.inLinks.begin();
								linkIt != nd.inLinks.end(); linkIt++) {
							int fromMod = superNodes[linkIt->first].ModIdx();

							if (fromMod == oldMod) {
								inFlowFromOldMod += beta * linkIt->second;
							} else if (fromMod == newMod) {
								inFlowFromNewMod += beta * linkIt->second;
							}
						}
					} else {
						for (link_iterator linkIt = nd.outLinks.begin();
								linkIt != nd.outLinks.end(); linkIt++) {
							int toMod = superNodes[linkIt->first].ModIdx();
							if (toMod == oldMod) {
								outFlowToOldMod += beta * linkIt->second;
							}
						}

						for (link_iterator linkIt = nd.inLinks.begin();
								linkIt != nd.inLinks.end(); linkIt++) {
							int fromMod = superNodes[linkIt->first].ModIdx();
							if (fromMod == oldMod) {
								inFlowFromOldMod += beta * linkIt->second;
							}
						}
					}

					oldExitPr1 = modules[oldMod].exitPr;
					oldSumPr1 = modules[oldMod].sumPr;
					oldSumDangling1 = modules[oldMod].sumDangling;
					oldModTPWeight = modules[oldMod].sumTPWeight;

					// For teleportation and danling nodes.
					outFlowToOldMod += (alpha * ndSize + beta * ndDanglingSize)
							* (oldModTPWeight - ndTPWeight);
					inFlowFromOldMod += (alpha * (oldSumPr1 - ndSize)
							+ beta * (oldSumDangling1 - ndDanglingSize))
							* ndTPWeight;
					outFlowToNewMod += (alpha * ndSize + beta * ndDanglingSize)
							* modules[newMod].sumTPWeight;
					inFlowFromNewMod += (alpha * modules[newMod].sumPr
							+ beta * modules[newMod].sumDangling) * ndTPWeight;

					if (isEmptyTarget) {
						outFlowToNewMod = 0.0;
						inFlowFromNewMod = 0.0;
					}

					moveResult.exitPr1 = oldExitPr1 - nd.ExitPr()
							+ outFlowToOldMod + inFlowFromOldMod;

					// copy module specific values...
					double oldExitPr2 = modules[newMod].exitPr;
					double oldSumPr2 = modules[newMod].sumPr;

					// Calculate status of current investigated movement of the node nd.
					moveResult.newModule = newMod;
					moveResult.sumPr1 = oldSumPr1 - ndSize;	// This should be 0.0, because oldModule will be empty module.
					moveResult.sumPr2 = oldSumPr2 + ndSize;
					moveResult.exitPr2 = oldExitPr2 + nd.ExitPr()
							- outFlowToNewMod - inFlowFromNewMod;

					moveResult.newSumExitPr = sumAllExitPr + moveResult.exitPr1
							+ moveResult.exitPr2 - oldExitPr1 - oldExitPr2;

					// Calculate delta_L(M) = L(M)_new - L(M)_old
					double delta_allExit_log_allExit = pLogP(
							moveResult.newSumExitPr) - pLogP(sumAllExitPr);
					double delta_exit_log_exit = pLogP(moveResult.exitPr1)
							+ pLogP(moveResult.exitPr2) - pLogP(oldExitPr1)
							- pLogP(oldExitPr2);
					double delta_stay_log_stay = pLogP(
							moveResult.exitPr1 + moveResult.sumPr1)
							+ pLogP(moveResult.exitPr2 + moveResult.sumPr2)
							- pLogP(oldExitPr1 + oldSumPr1)
							- pLogP(oldExitPr2 + oldSumPr2);

					// delta_L(M) = delta_allExit - 2.0 * delta_exit_log_exit + delta_stay_log_stay.
					moveResult.diffCodeLen = delta_allExit_log_allExit
							- 2.0 * delta_exit_log_exit + delta_stay_log_stay;

					////////////////////////////////////
					////// THE END OF EXAMINATION //////
					////////////////////////////////////

					if (isEmptyTarget) {
						emptyModules.pop_back();

						nEmptyMod--;
						nModule++;
					}

					nd.setModIdx(newMod);
					for (int k = 0; k < spMembers; k++) {
						nd.members[k]->setModIdx(newMod);
					}

					modules[newMod].numMembers += spMembers;
					modules[newMod].exitPr = moveResult.exitPr2;
					modules[newMod].sumPr = moveResult.sumPr2;
					modules[newMod].stayPr = moveResult.exitPr2
							+ moveResult.sumPr2;
					modules[newMod].sumTPWeight += ndTPWeight;

					if (ndDanglingSize > 0.0) {
						modules[newMod].sumDangling += ndDanglingSize;
						modules[oldMod].sumDangling -= ndDanglingSize;
					}

					// update related to the oldMod...
					modules[oldMod].numMembers -= spMembers;
					modules[oldMod].exitPr = moveResult.exitPr1;
					modules[oldMod].sumPr = moveResult.sumPr1;
					modules[oldMod].stayPr = moveResult.exitPr1
							+ moveResult.sumPr1;
					modules[oldMod].sumTPWeight -= ndTPWeight;

					if (modules[oldMod].numMembers == 0) {
						nEmptyMod++;
						nModule--;
						emptyModules.push_back(oldMod);
					}

					sumAllExitPr = moveResult.newSumExitPr;
					codeLength += moveResult.diffCodeLen;

					numMoved += spMembers;
				}

				gettimeofday(&tEnd, NULL);

				tSequential += elapsedTimeInSec(tStart, tEnd);
			}	// END critical

			// update activeNodes and isActives vectors.
			// We have to add the following nodes in activeNodes: neighbors, members in oldMod & newMod.
			// This can be done in parallel without any locking, since the written value is always same.
			for (link_iterator linkIt = nd.outLinks.begin();
					linkIt != nd.outLinks.end(); linkIt++)
				isActives[linkIt->first] = 1;		// set as an active nodes.

			for (link_iterator linkIt = nd.inLinks.begin();
					linkIt != nd.inLinks.end(); linkIt++)
				isActives[linkIt->first] = 1;		// set as an active nodes.
		}

	}	// END parallel for

	for (int i = 0; i < isActives.size(); i++) {
		if (isActives[i] == 1) {
			activeNodes.push_back(i);
			isActives[i] = 0;	// reset the flag of isActives[i].
		}
	}

	return numMoved;
}

/**
 *	This function will update members vector in modules correspondingly.
 */
void Network::updateMembersInModule() {

	int size, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int numModules = modules.size();

	vector<int>().swap(smActiveMods);
	smActiveMods.reserve(nModule);
	vector<int>().swap(lgActiveMods);
	lgActiveMods.reserve(nModule);

	for (int i = 0; i < numModules; i++) {
		vector<Node*>().swap(modules[i].members);
		if (modules[i].numMembers > 10000)
			lgActiveMods.push_back(i);
		else if (modules[i].numMembers > 0)
			smActiveMods.push_back(i);
	}

	for (int i = 0; i < nNode; i++) {
		modules[nodes[i].ModIdx()].members.push_back(&nodes[i]);
	}
}

void Network::updateSPMembersInModule() {
	int numModules = modules.size();

	for (int i = 0; i < numModules; i++)
		vector<Node*>().swap(modules[i].members);

	int numSPNodes = superNodes.size();
	for (int i = 0; i < numSPNodes; i++)
		modules[superNodes[i].ModIdx()].members.push_back(&superNodes[i]);
}

// calculate exit-probability based on node information and current module assignment.
void Network::updateCodeLength(int numTh, bool isSPNode) {

// calculate exit-probability for each module.
	double tempSumAllExit = 0.0;
	double exit_log_exit = 0.0;
	double stay_log_stay = 0.0;

	omp_set_num_threads(numTh);

#pragma omp parallel for reduction (+:tempSumAllExit, exit_log_exit, stay_log_stay)
	for (unsigned int i = 0; i < nNode; i++) {
		if (modules[i].numMembers > 0) {
			// exitPr w.r.t. teleportation.
			double exitPr = Network::alpha * (1.0 - modules[i].sumTPWeight)
					* modules[i].sumPr;
			int curMod = modules[i].index;

			double sumOutFlow = 0.0;
			int nMembers = modules[i].members.size();
			for (int j = 0; j < nMembers; j++) {
				Node * nd = modules[i].members[j];
				int nOutLinks = nd->outLinks.size();

				for (int k = 0; k < nOutLinks; k++) {
					if (!isSPNode
							&& nodes[nd->outLinks[k].first].ModIdx() != curMod)
						sumOutFlow += nd->outLinks[k].second;
					else if (isSPNode
							&& superNodes[nd->outLinks[k].first].ModIdx()
									!= curMod)
						sumOutFlow += nd->outLinks[k].second;
				}
			}

			exitPr += Network::beta
					* (sumOutFlow
							+ (1.0 - modules[i].sumTPWeight)
									* modules[i].sumDangling);

			modules[i].exitPr = exitPr;
			modules[i].stayPr = exitPr + modules[i].sumPr;

			tempSumAllExit += exitPr;
			exit_log_exit += pLogP(exitPr);
			stay_log_stay += pLogP(modules[i].stayPr);
		}
	}

	sumAllExitPr = tempSumAllExit;

	codeLength = pLogP(sumAllExitPr) - 2.0 * exit_log_exit + stay_log_stay
			- allNodes_log_allNodes;
}

// calculate exit-probability based on node information and current module assignment.
double Network::calculateCodeLength() {

// calculate exit-probability for each module.
	double tempSumAllExit = 0.0;
	double exit_log_exit = 0.0;
	double stay_log_stay = 0.0;

#pragma omp parallel for reduction (+:tempSumAllExit, exit_log_exit, stay_log_stay)
	for (unsigned int i = 0; i < nNode; i++) {
		if (modules[i].numMembers > 0) {
			// exitPr w.r.t. teleportation.
			double exitPr = Network::alpha * (1.0 - modules[i].sumTPWeight)
					* modules[i].sumPr;
			int curMod = modules[i].index;

			double sumOutFlow = 0.0;
			int nMembers = modules[i].members.size();
			for (int j = 0; j < nMembers; j++) {
				int nOutLinks = modules[i].members[j]->outLinks.size();
				for (int k = 0; k < nOutLinks; k++) {
					if (nodes[modules[i].members[j]->outLinks[k].first].ModIdx()
							!= curMod)
						sumOutFlow += modules[i].members[j]->outLinks[k].second;
				}
			}

			exitPr += Network::beta
					* (sumOutFlow
							+ (1.0 - modules[i].sumTPWeight)
									* modules[i].sumDangling);

			tempSumAllExit += exitPr;
			exit_log_exit += pLogP(exitPr);
			stay_log_stay += pLogP(exitPr + modules[i].sumPr);
		}
	}

	double computedCodeLength = pLogP(tempSumAllExit) - 2.0 * exit_log_exit
			+ stay_log_stay - allNodes_log_allNodes;

	return computedCodeLength;
}

// similar function of Greedy::level() in the original infomap_dir implementation.
// make a group of nodes as a Node (here, SuperNode)
// Granularity of the network is same but a group of nodes move together instead of moving individually.
void Network::convertModulesToSuperNodes(int tag) {

//initialize superNodes vector for updating w.r.t. the current module status.

	vector<SuperNode>().swap(superNodes);
	superNodes.reserve(nModule);
// parallelization starting from here for now

	int size, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int numSPNodesforCurrentRank;
	const int toSpNodesPackSize = nEdge + 1;
	int toSpNodesSend[toSpNodesPackSize] = { };
	double linkWeightstoSpNodesSend[toSpNodesPackSize] = { };
	int toSpNodesReceive[toSpNodesPackSize] = { };
	const int totalEdges = nEdge;
	double linkWeightstoSpNodesReceive[totalEdges] = { };

	vector<unsigned int> modToSPNode(nNode);// indicate corresponding superNode ID from each module.

// I am not going to parallelize the following code snippet for now, but will parallelize later
	int idx = 0;
	for (unsigned int i = 0; i < nNode; i++) {
		modules[i].numMembers = modules[i].members.size(); //a memory error is happening because of this statement, I don't know why, maybe this is the source of the error or maybe the error is being thrown from another source and this code is revealing the error of another source
		if (modules[i].numMembers > 0) {
			superNodes.push_back(SuperNode(modules[i], idx));
			modToSPNode[i] = idx;
			idx++;
		}
	}

	/*
	 *	Calculate outLinks and inLinks between superNodes...
	 */
	int numSPNode = superNodes.size();
//nModule = superNodes.size();

	int start, end;
	findAssignedPart(&start, &end, numSPNode, size, rank);
	numSPNodesforCurrentRank = end - start;
	toSpNodesSend[nEdge] = numSPNodesforCurrentRank;

	const int PackSize = 2 * nNode;
	int spNodeCountTrackerSend[PackSize] = { };
	int spNodeCountTrackerReceive[PackSize] = { };
	int countSPNode = 0;

	int totalLinksCount = 0;

	for (int i = start; i < end; i++) {

		spNodeCountTrackerSend[countSPNode] = i;
		int numNodesInSPNode = superNodes[i].members.size();

		typedef map<int, double> EdgeMap;
		EdgeMap newOutLinks;
		map<int, int> NodeIndexMap;
		int noOutLinks = 0;
		int toSPNode;
		int currentSPNodeTotalLinks = 0;

		for (int j = 0; j < numNodesInSPNode; j++) {
			// Calculate newOutLinks from a superNode to other superNodes.
			Node* nd = superNodes[i].members[j];
			int nOutEdge = nd->outLinks.size();
			//int totalLinksCount = 0;

			for (int k = 0; k < nOutEdge; k++) {
				toSPNode = modToSPNode[nodes[nd->outLinks[k].first].ModIdx()];
				if (toSPNode != i) {	// outLinks to another superNode...
					pair<map<int, int>::iterator, bool> ret =
							NodeIndexMap.insert(
									make_pair(toSPNode, totalLinksCount));
					if (ret.second) {
						toSpNodesSend[totalLinksCount] = toSPNode;
						linkWeightstoSpNodesSend[totalLinksCount] =
								nd->outLinks[k].second;
						totalLinksCount++;
						currentSPNodeTotalLinks++;
					} else {
						int spNdIndex = ret.first->second;
						linkWeightstoSpNodesSend[spNdIndex] +=
								nd->outLinks[k].second;
					}
				}
			}
		}
		spNodeCountTrackerSend[nNode + countSPNode] = currentSPNodeTotalLinks;
		countSPNode++;

	}

	int totalLinks = 0;
	for (int m = 0; m < toSpNodesSend[nEdge]; m++) {
		int spNodeIndex = spNodeCountTrackerSend[m];
		int numLinks = spNodeCountTrackerSend[nNode + m];
		superNodes[m].outLinks.reserve(numLinks);

		for (int it = totalLinks; it < totalLinks + numLinks; it++) {
			superNodes[spNodeIndex].outLinks.push_back(
					make_pair(toSpNodesSend[it], linkWeightstoSpNodesSend[it]));
		}
		totalLinks += numLinks;
	}

	for (int prId = 0; prId < size; prId++) {
		if (prId != rank) {
			MPI_Send(spNodeCountTrackerSend, PackSize, MPI_INT, prId, tag,
			MPI_COMM_WORLD);
			MPI_Send(toSpNodesSend, toSpNodesPackSize, MPI_INT, prId, tag,
			MPI_COMM_WORLD);
			MPI_Send(linkWeightstoSpNodesSend, toSpNodesPackSize, MPI_DOUBLE,
					prId, tag,
					MPI_COMM_WORLD);
		}
	}

	for (int sId = 0; sId < size; sId++) {
		if (sId != rank) {
			MPI_Recv(spNodeCountTrackerReceive, PackSize, MPI_INT, sId, tag,
			MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
			MPI_Recv(toSpNodesReceive, toSpNodesPackSize, MPI_INT, sId, tag,
			MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
			MPI_Recv(linkWeightstoSpNodesReceive, toSpNodesPackSize, MPI_DOUBLE,
					sId, tag,
					MPI_COMM_WORLD, MPI_STATUSES_IGNORE);

			int spNodeReceived = toSpNodesReceive[nEdge];

			int totLinks = 0;

			for (int it = 0; it < spNodeReceived; it++) {
				int spNodeIndex = spNodeCountTrackerReceive[it];
				int numLinks = spNodeCountTrackerReceive[nNode + it];
				superNodes[spNodeIndex].outLinks.reserve(numLinks);

				for (int iterator = totLinks; iterator < totLinks + numLinks;
						iterator++) {
					superNodes[spNodeIndex].outLinks.push_back(
							make_pair(toSpNodesReceive[iterator],
									linkWeightstoSpNodesReceive[iterator]));
				}
				totLinks += numLinks;
			}
		}
	}

// update inLinks in SEQUENTIAL..
	for (int i = 0; i < numSPNode; i++) {
		int nOutLinks = superNodes[i].outLinks.size();
		{
			for (int j = 0; j < nOutLinks; j++)
				superNodes[superNodes[i].outLinks[j].first].inLinks.push_back(
						make_pair(i, superNodes[i].outLinks[j].second));
		}
	}

}

void Network::displayOutlinksforSuperNodes() {

	int size, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int numSPNode = superNodes.size();
	for (int it = 0; it < numSPNode; it++) {
		int noOutLinks = superNodes[it].outLinks.size();
		printf(
				"filtering by the process:%d, spNode:%d, number of outlinks:%d\n",
				rank, it, noOutLinks);
		for (int ite = 0; ite < noOutLinks; ite++) {
			printf(
					"filter by the process:%d, superNodes [%d].outLinks[%d]->otherEnd:%d, weight:%f\n",
					rank, it, ite, superNodes[it].outLinks[ite].first,
					superNodes[it].outLinks[ite].second);
		}
	}
}

void Network::displayInlinksforSuperNodes() {

	int size, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int numSPNode = superNodes.size();
	for (int it = 0; it < numSPNode; it++) {
		int noinLinks = superNodes[it].inLinks.size();
		printf("refining by the process:%d, spNode:%d, number of inlinks:%d\n",
				rank, it, noinLinks);
		for (int ite = 0; ite < noinLinks; ite++) {
			printf(
					"refine by the process:%d, superNodes [%d].inLinks[%d]->otherEnd:%d, weight:%f\n",
					rank, it, ite, superNodes[it].inLinks[ite].first,
					superNodes[it].inLinks[ite].second);
		}
	}
}

//void Network::generateSuperNodesFromSubModules() {
void Network::generateSuperNodesFromSubModules(int numTh) {
//initialize superNodes vector for updating w.r.t. the current module status.
	vector<SuperNode>().swap(superNodes);

	int nSubMods = subModules.size();
	superNodes.reserve(nSubMods);

	for (int i = 0; i < nSubMods; i++) {
		superNodes.push_back(SuperNode(subModules[i], i, *this));
	}

	/*
	 *	Calculate outLinks and inLinks between superNodes...
	 */
	int numSPNode = superNodes.size();

	typedef map<int, double> EdgeMap;

	omp_set_num_threads(numTh);

#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < numSPNode; i++) {
		SuperNode* spNode = &(superNodes[i]);

		int numNodesInSPNode = spNode->members.size();

		EdgeMap newOutLinks;

		for (int j = 0; j < numNodesInSPNode; j++) {
			// Calculate newOutLinks from a superNode to other superNodes.
			Node* nd = spNode->members[j];
			int nOutEdge = nd->outLinks.size();

			for (int k = 0; k < nOutEdge; k++) {
				int toSPNode = ndToSubMod[nd->outLinks[k].first];

				if (toSPNode != i) {	// outLinks to another superNode...
					pair<EdgeMap::iterator, bool> ret = newOutLinks.insert(
							make_pair(toSPNode, nd->outLinks[k].second));
					if (!ret.second)
						ret.first->second += nd->outLinks[k].second;
				}
			}
		}

		spNode->outLinks.reserve(newOutLinks.size());

		double sumExitFlow = 0.0;

		for (EdgeMap::iterator it = newOutLinks.begin();
				it != newOutLinks.end(); it++) {
			spNode->outLinks.push_back(make_pair(it->first, it->second));
			sumExitFlow += it->second;
		}

		double teleportExitFlow = 0.0;
		for (int j = 0; j < numNodesInSPNode; j++) {
			Node* nd = spNode->members[j];
			teleportExitFlow += (alpha * nd->Size() + beta * nd->DanglingSize())
					* (1.0 - nd->TeleportWeight());
		}

		if (sumExitFlow == 0.0) {
			spNode->setExitPr(teleportExitFlow);
		} else {
			double exitProb = teleportExitFlow + beta * sumExitFlow;
			spNode->setExitPr(exitProb);
		}

	}

// update inLinks in SEQUENTIAL..
	for (int i = 0; i < numSPNode; i++) {
		int nOutLinks = superNodes[i].outLinks.size();
		for (int j = 0; j < nOutLinks; j++)
			superNodes[superNodes[i].outLinks[j].first].inLinks.push_back(
					make_pair(i, superNodes[i].outLinks[j].second));
	}
}

void Network::copyModule(Module * newM, Module * oldM) {
	newM->index = oldM->index;
	newM->exitPr = oldM->exitPr;
	newM->stayPr = oldM->stayPr;
	newM->sumPr = oldM->sumPr;
	newM->sumTPWeight = oldM->sumTPWeight;
	newM->sumDangling = oldM->sumDangling;
	newM->numMembers = oldM->numMembers;

	newM->members.clear();
	int nMembers = oldM->members.size();
	for (int i = 0; i < nMembers; i++)
		newM->members.push_back(oldM->members[i]);
}

void Network::copyModule(int newM, int oldM) {
	modules[newM].index = modules[oldM].index;
	modules[newM].exitPr = modules[oldM].exitPr;
	modules[newM].stayPr = modules[oldM].stayPr;
	modules[newM].sumPr = modules[oldM].sumPr;
	modules[newM].sumTPWeight = modules[oldM].sumTPWeight;
	modules[newM].sumDangling = modules[oldM].sumDangling;
	modules[newM].numMembers = modules[oldM].numMembers;

	modules[newM].members.clear();
	int nMembers = modules[oldM].members.size();
	for (int i = 0; i < nMembers; i++)
		modules[newM].members.push_back(modules[oldM].members[i]);
}

// Miscellaneous function..
void findAssignedPart(int* start, int* end, int numNodes, int numTh, int myID) {
	if (numNodes % numTh == 0) {
		*start = (numNodes / numTh) * myID;
		*end = (numNodes / numTh) * (myID + 1);
	} else {
		int block = numNodes / numTh;
		int modular = numNodes % numTh;

		if (myID < modular) {
			*start = (block + 1) * myID;
			*end = (block + 1) * (myID + 1);
		} else {
			*start = block * myID + modular;
			*end = block * (myID + 1) + modular;
		}
	}
}

