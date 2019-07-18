/*****************************************************

	Use parallel genetic algorithm to find
	the minimum of Griewank function

******************************************************/

#include <bits/stdc++.h>
#include <mpi.h>

using namespace std;

//number of individuals in each generation
#define POPULATION_SIZE 150

//dimensions
int N = pow(2,9);

//precision
double prs = 1.0;

//genes
vector<double> GENES;

//target [0, ..., 0]
vector<double> TARGET(N, 0);

//function to generate random numbers in given range 
int random_num(int start, int end){
	int range = (end-start)+1;
	int random_int = start+(rand()%range);
	return random_int;
}

//create random genes for mutation 
double mutated_genes(){
	int r = random_num(0, GENES.size()-1);
	return GENES[r];
}

//create chromosome or string of genes 
vector<double> create_gnome(){
	vector<double> gnome;
	for(int i=0; i<TARGET.size(); i++){
		gnome.push_back(mutated_genes());
	}
	return gnome;
}

//perform mating and produce new offspring 
vector<double> mate(vector<double> par1, vector<double> par2){
	//chromosome for offspring
	vector<double> child_chromosome;

	int len = par1.size();
	for(int i = 0;i<len;i++){
		//random probability
		double p = random_num(0, 100)/100;

		//if prob is less than 0.45, insert gene
		//from parent 1
		if(p < 0.45){
			child_chromosome.push_back(par1[i]);
		}
		//if prob is between 0.45 and 0.90, insert 
		//gene from parent 2 
		else if(p < 0.90){
			child_chromosome.push_back(par2[i]);
		}
		//otherwise insert random gene(mutate), 
		//for maintaining diversity 
		else{
			child_chromosome.push_back(mutated_genes());
		}
	}
	//create new Individual(offspring) using
	//generated chromosome for offspring
	return child_chromosome;
};

//calculate fittness score, it is the number of
//characters in string which differ from target
//string.
int cal_fitness(vector<double> chromosome){
	int len = TARGET.size();
	int fitness = 0;
	for(int i = 0;i<len;i++){
		if(chromosome[i] != TARGET[i])
			fitness++;
	}
	return fitness;
};

//compare the individuals of population
bool compare(const vector<double> chromo1, const vector<double> chromo2){
	return cal_fitness(chromo1) < cal_fitness(chromo2);
}

//overloading << operator
ostream &operator<<(std::ostream &out, const vector<double> &vec){
	for(int i=0; i<vec.size(); i++){
		out << vec[i] << " ";
	}
	out << "\t";
	return out;
}

//driver code
int main(){
	srand((unsigned)(time(NULL)));
	
	int rank, core;
	MPI_Init(NULL,NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &core);

	double T = 0.0;
	T -= MPI_Wtime();

	//initialize GENES [-1, 1]
	for(int i=0; i<1/prs*2; i++){
		GENES.push_back(-1+i*1/prs);
	}

	//current generation
	int generation = 0;

	vector<vector<double> > population;
	double mpi_population[POPULATION_SIZE][N];
	bool found = false;

	if(rank == 0){
		//create initial population (random)
		for(int i = 0;i<POPULATION_SIZE;i++){
			vector<double> gnome = create_gnome();
			for(int j=0; j<N; j++){
				mpi_population[i][j] = gnome[j];
			}
		}
	}
	MPI_Bcast(mpi_population, POPULATION_SIZE*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	for(int i=0; i<POPULATION_SIZE; i++){
		vector<double> vec_population;
		for(int j=0; j<N; j++){
			vec_population.push_back(mpi_population[i][j]);
		}
		population.push_back(vec_population);
	}

	while(! found){
		//sort the population in increasing order of fitness score
		sort(population.begin(), population.end(), compare);

		//if the individual having lowest fitness score ie.
		//0 then we know that we have reached to the target
		//and break the loop
		if(cal_fitness(population[0]) <= 0){
			found = true;
			break;
		}

		//otherwise generate new offsprings for new generation
		vector<vector<double> > new_generation;

		//perform Elitism, that mean 10-20% of fittest population
		//goes to the next generation
		int s = 38;
		for(int i=0; i<s; i++){
			new_generation.push_back(population[i]);	
		}

		//from 50% of fittest population, Individuals
		//will mate to produce offspring
		int count = (POPULATION_SIZE-s)/core;
		double temp_generation[count][N];

		for(int i=0; i<count; i++){
			int len = population.size();
			int r = random_num(0, POPULATION_SIZE/2);
			vector<double> parent1 = population[r];
			r = random_num(0, POPULATION_SIZE/2);
			vector<double> parent2 = population[r];
			vector<double> offspring = mate(parent1, parent2);
			for(int j=0; j<N; j++){
				temp_generation[i][j] = offspring[j];
			}
		}
		
		//gather temp_generation
		double gather_generation[s][N];
		MPI_Gather(temp_generation, count*N, MPI_DOUBLE, gather_generation, count*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(gather_generation, POPULATION_SIZE*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		for(int i=0; i<count*core; i++){
			vector<double> vec_generation;
			for(int j=0; j<N; j++){
				vec_generation.push_back(gather_generation[i][j]);
			}
			new_generation.push_back(vec_generation);
		}
		population = new_generation;

		if(rank == 0){
			cout<< "Generation: " << generation << "\t";
			cout<< "Fitness: "<< cal_fitness(population[0]) << "\n";
			generation++;
		}
		
	}

	T += MPI_Wtime();
	if(rank == 0){
	cout<< "Generation: " << generation << "\t";
	cout<< "Solution: "<< population[0] <<"\t";
	cout<< "Fitness: "<< cal_fitness(population[0]) << "\n";
	cout << "Time = " << T << endl;
	}
	MPI_Finalize();
}
