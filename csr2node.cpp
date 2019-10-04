#include<iostream>
#include<fstream>
#include<vector>
#include<cstdio>
#include<cstdlib>
#include<cstring>
#include<set>
#include<utility>

using namespace std;


int main(int argc, char* argv[])
{
	if(argc<3)
	{
		cout<<"usage: program inputfile outputfile argument"<<endl;
		return 0;
	}
	
	long vertices = 0, edges=0;
	char* token;
	FILE* fp = fopen(argv[1], "r");
	ofstream output(argv[2]);
	long line_counter = 0;
	
	if(fp==NULL)
	{
		exit(EXIT_FAILURE);
	}
	
	cout<<"File open successful\n";
	
	char* line = NULL;
	size_t len = 0;
	set<long> Nodes;
	
	if((getline(&line, &len, fp)!= -1))
	{
		token = strtok(line," \n\t");
		if(token!= NULL)
		{
			vertices = atoi(token);
			printf("number of vertices:%s\n", token);
		}
		token = strtok(NULL," \n\t");
		if(token!= NULL)
		{
			edges = atoi(token);
			printf("number of edges:%s\n", token);
		}
	}
	
	while((getline(&line, &len, fp)!= -1))
	{

		printf("%s", line);
		
		if(strlen(line)<=2)
		{
			printf("why not reached?\n");
			break;
		}
		token = strtok(line," \n\t");
		while(token != NULL){
			long edge = atoi(token);
			printf("line_counter:%ld, edge:%ld\n",line_counter, edge);
			Nodes.insert(edge);
			token = strtok(NULL," \n\t");
		}
		line_counter++;
	}
	fclose(fp);
	if(line)
	{
		free(line);
	}
	
	for (set<long>::iterator it = Nodes.begin(); it != Nodes.end(); it++)
	{
	    output<<*it<<endl; // Note the "*" here
	}
	
	output.close();
	
 return 0;	
}
