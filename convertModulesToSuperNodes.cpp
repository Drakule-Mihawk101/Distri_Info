int size, rank;

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

//initialize superNodes vector for updating w.r.t. the current module status.
	vector<SuperNode>().swap(superNodes);
	superNodes.reserve(nModule);

	vector<unsigned int> modToSPNode(nNode);// indicate corresponding superNode ID from each module.

	int idx = 0;
	for (unsigned int i = 0; i < nNode; i++) {
		if (modules[i].numMembers > 0) {
			superNodes.push_back(SuperNode(modules[i], idx));
			modToSPNode[i] = idx;
			idx++;
		}
	}

	//cout<<"size: "<<size<<" ,rank: "<<rank<<endl;

	/*
	 *	Calculate outLinks and inLinks between superNodes...
	 */
	int numSPNode = superNodes.size();
	int maxSPNodesPerProcess = (numSPNode / size) + (size - 1);
/*	if (rank == 0) {
		maxSPNodesPerProcess = (numSPNode / size) + (size - 1);
	}

	MPI_Bcast(&maxSPNodesPerProcess, 1, MPI_INT, 0, MPI_COMM_WORLD);*/

	struct edgeList {
		int totalLinksCount;
		int spNodeCount;
		int* listSpNodes;		//list all the superNodes under a process
		int* listLinkSizes;		//number of outlinks conennected to superNode i
		int* toSpNodes;		//this is the superNode array for the Outlinks
		double* linkWeightstoSpNodes;//this is the edgeWeight array for the outlinks
	};

	edgeList list1, list2;

	MPI_Datatype listType1;
	int sctlen = 6;
	int blklengths[sctlen];
	MPI_Datatype tp[sctlen];
	MPI_Aint displs[sctlen];

	list1.totalLinksCount = 0;
	list1.spNodeCount = 0;
	list1.listSpNodes = new int[maxSPNodesPerProcess]();
	list1.listLinkSizes = new int[maxSPNodesPerProcess]();
	list1.toSpNodes = new int[nEdge]();
	list1.linkWeightstoSpNodes = new double[nEdge]();

	list2.totalLinksCount = 0;
	list2.spNodeCount = 0;
	list2.listSpNodes = new int[maxSPNodesPerProcess]();
	list2.listLinkSizes = new int[maxSPNodesPerProcess]();
	list2.toSpNodes = new int[nEdge]();
	list2.linkWeightstoSpNodes = new double[nEdge]();

	cout << "maxSPNodesPerProcess:" << maxSPNodesPerProcess << endl;

	blklengths[0] = 1;
	tp[0] = MPI_INT;
	displs[0] = (size_t) &(list1.totalLinksCount) - (size_t) &list1;

	blklengths[1] = 1;
	tp[1] = MPI_INT;
	displs[1] = (size_t) &(list1.spNodeCount) - (size_t) &list1;

	blklengths[2] = maxSPNodesPerProcess;
	tp[2] = MPI_INT;
	displs[2] = (size_t) &(list1.listSpNodes[0]) - (size_t) &list1;

	blklengths[3] = maxSPNodesPerProcess;
	tp[3] = MPI_INT;
	displs[3] = (size_t) &(list1.listLinkSizes[0]) - (size_t) &list1;

	blklengths[4] = nEdge;
	tp[4] = MPI_INT;
	displs[4] = (size_t) &(list1.toSpNodes[0]) - (size_t) &list1;

	blklengths[5] = nEdge;
	tp[5] = MPI_DOUBLE;
	displs[5] = (size_t) &(list1.linkWeightstoSpNodes[0]) - (size_t) &list1;

	MPI_Type_create_struct(sctlen, blklengths, displs, tp, &listType1);

	MPI_Type_commit(&listType1);
	{
		MPI_Aint typesize1;
		MPI_Type_extent(listType1, &typesize1);
	}

	MPI_Datatype listType2;
	int sctlen2 = 6;
	int blklengths2[sctlen2];
	MPI_Datatype tp2[sctlen2];
	MPI_Aint displs2[sctlen2];

	blklengths2[0] = 1;
	tp2[0] = MPI_INT;
	displs2[0] = (size_t) &(list2.totalLinksCount) - (size_t) &list2;

	blklengths2[1] = 1;
	tp2[1] = MPI_INT;
	displs2[1] = (size_t) &(list2.spNodeCount) - (size_t) &list2;

	blklengths2[2] = maxSPNodesPerProcess;
	tp2[2] = MPI_INT;
	displs2[2] = (size_t) &(list2.listSpNodes[0]) - (size_t) &list2;

	blklengths2[3] = maxSPNodesPerProcess;
	tp2[3] = MPI_INT;
	displs2[3] = (size_t) &(list2.listLinkSizes[0]) - (size_t) &list2;

	blklengths2[4] = nEdge;
	tp2[4] = MPI_INT;
	displs2[4] = (size_t) &(list2.toSpNodes[0]) - (size_t) &list2;

	blklengths2[5] = nEdge;
	tp2[5] = MPI_DOUBLE;
	displs2[5] = (size_t) &(list2.linkWeightstoSpNodes[0]) - (size_t) &list2;

	MPI_Type_create_struct(sctlen2, blklengths2, displs2, tp2, &listType2);

	MPI_Type_commit(&listType2);
	{
		MPI_Aint typesize2;
		MPI_Type_extent(listType2, &typesize2);
	}
	//omp_set_num_threads(numTh);

	/*#pragma omp parallel
	 {
	 int myID = omp_get_thread_num();
	 int nTh = omp_get_num_threads();*/

	int start, end;
	findAssignedPart(&start, &end, numSPNode, size, rank);

	//cout<<"size: "<<size<<" ,rank: "<<rank<<endl;

	int countSPNode = 0;
	for (int i = start; i < end; i++) {

		int numNodesInSPNode = superNodes[i].members.size();

		/*edgeList1.listSuperNodes[countSPNode] = i;
		 cout << "rank:" << rank << "i:" << i << endl;*/

		list1.listSpNodes[countSPNode] = i;
		printf("$$$$rank:%d,countSPNode:%d,superNode:%d\n", rank, countSPNode,
				list1.listSpNodes[countSPNode]);
		typedef map<int, double> EdgeMap;
		EdgeMap newOutLinks;
		map<int, int> NodeIndexMap;
		int noOutLinks = 0;

		for (int j = 0; j < numNodesInSPNode; j++) {
			// Calculate newOutLinks from a superNode to other superNodes.
			Node* nd = superNodes[i].members[j];
			int nOutEdge = nd->outLinks.size();

			for (int k = 0; k < nOutEdge; k++) {
				int toSPNode =
						modToSPNode[nodes[nd->outLinks[k].first].ModIdx()];

				if (toSPNode != i) {	// outLinks to another superNode...
					/*						pair<EdgeMap::iterator, bool> ret = newOutLinks.insert(
					 make_pair(toSPNode, nd->outLinks[k].second));*/
					pair<map<int, int>::iterator, bool> re =
							NodeIndexMap.insert(
									make_pair(toSPNode, list1.totalLinksCount));

					if (re.second) { //it means the superNode is already inserted for an existing link
						list1.toSpNodes[list1.totalLinksCount] = toSPNode;
						list1.linkWeightstoSpNodes[list1.totalLinksCount] =
								nd->outLinks[k].second;
						list1.totalLinksCount++;
						noOutLinks++;
					} else {
						int spIndex = re.first->second; //get the index value from the superNode <--> index pair
						//NodeIndexMap.find(ret.first->first)->second; //get the index value from the superNode <--> index pair
						//ret.first->second += nd->outLinks[k].second;
						list1.linkWeightstoSpNodes[spIndex] +=
								nd->outLinks[k].second;
					}
				}
			}
		}
		list1.listLinkSizes[countSPNode] = noOutLinks;
		countSPNode++;

		/*		superNodes[i].outLinks.reserve(newOutLinks.size());

		 for (EdgeMap::iterator it = newOutLinks.begin();
		 it != newOutLinks.end(); it++) {
		 superNodes[i].outLinks.push_back(make_pair(it->first, it->second));
		 }*/
		//cout << "aaaaaaaaa" << endl;
	}
	list1.spNodeCount = countSPNode;

	printf("rank:%d,start:%d,end:%d\n", rank, start, end);
	for (int n = 0; n < countSPNode; n++) {
		printf("##########rank:%d,i:%d###########\n", rank,
				list1.listSpNodes[n]);
	}
	for (int prId = 0; prId < size; prId++) {
		if (prId != rank) {
			MPI_Send(&list1, 1, listType1, prId, 0, MPI_COMM_WORLD);
			cout << "ccccccccc:rank:prId:" << rank << prId << endl;
		}
	}

	for (int sId = 0; sId < size; sId++) {
		//if (rank == 0) {
		if (sId != rank) {
			MPI_Recv(&list2, 1, listType2, MPI_ANY_SOURCE, 0,
			MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
			int totLinks = 0;
			printf("************************************list2.spNodeCount:%d\n",
					list2.spNodeCount);
			for (int j = 0; j < list2.spNodeCount; j++) {

				printf("list2.listSpNodes[%d], value:%d\n", j,
						list2.listSpNodes[j]);
				printf("list2.toSpNodes[%d], value:%d\n", j,
						list2.toSpNodes[j]);
				int spIdx = list2.listSpNodes[j];
				int numLinks = list2.listLinkSizes[j];
				printf("spIdx:%d,numLinks:%d\n", spIdx, numLinks);
				//cout << "spIdx:" << spIdx << "numLinks:" << numLinks << endl;
				superNodes[spIdx].outLinks.reserve(numLinks);
				for (int k = totLinks; k < totLinks + numLinks; k++) {
					superNodes[spIdx].outLinks.push_back(
							make_pair(list2.toSpNodes[k],
									list2.linkWeightstoSpNodes[k]));
				}
				totLinks += numLinks;
			}
			//cout << "ddddddddd:rank:" << rank << endl;
		}
	}
	/*for (int sId = 0; sId < size; sId++) {
	 if (sId != rank) {
	 MPI_Recv(&edgeList2, 1, edgeListType, MPI_ANY_SOURCE, 0,
	 MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
	 cout << "ddddddddd:rank:sId:" << rank << sId << endl;
	 }
	 }*/
//}
// update inLinks in SEQUENTIAL..
	for (int i = 0; i < numSPNode; i++) {
		int nOutLinks = superNodes[i].outLinks.size();
		for (int j = 0; j < nOutLinks; j++)
			superNodes[superNodes[i].outLinks[j].first].inLinks.push_back(
					make_pair(i, superNodes[i].outLinks[j].second));
	}

	MPI_Type_free(&listType1);
	MPI_Type_free(&listType2);


