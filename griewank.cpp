/************************************************

	The minimum of Griewank function

*************************************************/

#include <bits/stdc++.h>

using namespace std;

//number of individuals in each generation
int POPULATION_SIZE = 100;

//dimensions
int N = 512;

//genes
vector<double> GENES;

//max iteration
int Max_iter = 1000;

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
	for(int i=0; i<N; i++){
		gnome.push_back(mutated_genes());
	}
	return gnome;
}

//perform mating and produce new offspring 
vector<double> mate(vector<double> par1, vector<double> par2){
	//chromosome for offspring
	vector<double> child_chromosome;
	for(int i = 0;i<N;i++){
		//random probability
		float p = random_num(0, 100)/100;

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

//calculate fittness score (griewank function), it is the number
//of characters in string which differ from target string.
double cal_fitness(vector<double> chromosome){
	double part1 = 0.0, part2 = 1.0, fitness = 0.0;
	for(int i=0; i<N; i++){
		part1 += pow(chromosome[i],2);
		part2 *= cos(chromosome[i]/sqrt(i+1));
	}
	fitness = 1+1/4000*part1-part2;
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
	
	//initialize GENES [-1, 1], precision = 0.1
	for(int i=0; i<20; i++){
		GENES.push_back(-1+i*0.1);
	}

	//current generation
	int generation = 0;

	//create initial population
	vector<vector<double> > population;
	for(int i = 0;i<POPULATION_SIZE;i++){
		vector<double> gnome = create_gnome();
		population.push_back(gnome);
	}

	while(generation < Max_iter){
		//sort the population in increasing order of fitness score
		sort(population.begin(), population.end(), compare);

		//if the individual having lowest fitness score ie.
		//0 then we know that we have reached to the target
		//and break the loop

		//otherwise generate new offsprings for new generation
		vector<vector<double> > new_generation;

		//perform Elitism, that mean 10% of fittest population
		//goes to the next generation
		int s = (10*POPULATION_SIZE)/100; 
		for(int i = 0;i<s;i++){
			new_generation.push_back(population[i]);	
		}

		//from 50% of fittest population, Individuals
		//will mate to produce offspring
		s = (90*POPULATION_SIZE)/100;
		for(int i = 0;i<s;i++){
			int len = population.size();
			int r = random_num(0, 50);
			vector<double> parent1 = population[r];
			r = random_num(0, 50);
			vector<double> parent2 = population[r];
			vector<double> offspring = mate(parent1, parent2);
			new_generation.push_back(offspring);
		}
		
		population = new_generation;
		cout<< "Generation: " << generation << "\t";
		cout<< "Fitness: "<< cal_fitness(population[0]) << "\n";

		generation++;
	}
	cout<< "Generation: " << generation << "\t";
	cout<< "Fitness: "<< cal_fitness(population[0]) << "\n";
	
	return 0;
}
