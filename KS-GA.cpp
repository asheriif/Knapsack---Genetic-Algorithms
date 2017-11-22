#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>
#include<iostream>
#include<algorithm>

using namespace std;

#define rD (rand() / (double) RAND_MAX)		// random decimal

#define popsize 200
#define ngenerations 400
#define elite 1
const double pCross = 0.85;			// probability of crossover
const double pMutate = 0.1;			// probability of mutation

struct chromosome
{
	unsigned int genes;
	int fitness;
	int weight;
	chromosome() { }
};

int ngenes,s;		// number of items, size of the knapsack
int w[55];		// weights
int v[55];		// values
double cdf[popsize];	// cumulative probability distribution

chromosome cPopulation[popsize];	// current population
chromosome nPopulation[popsize];	// next population

void init()
{
	memset(cPopulation, 0 , sizeof(cPopulation));
	srand(time(0));
	for(int i = 0; i < popsize; i++)
	{	
		int cW = 0;						// current weight
		for(int j = 0; j < ngenes && cW < s; j++)
		{
			if(rand() % 2 && cW + w[j] <= s)
			{
				cPopulation[i].genes |= (1<<j);
				cW += w[j];
			}
		}
	}
}

void calcFitness()
{
	int totalFitness = 0;
	memset(cdf,0,sizeof(cdf));
	for(int i = 0; i < popsize; i++)
	{
		cPopulation[i].fitness = 0;
		cPopulation[i].weight = 0;
		for(int j = 0; j < ngenes; j++)
		{
			cPopulation[i].fitness += ((cPopulation[i].genes >> j) & 1) * v[j];
			cPopulation[i].weight  += ((cPopulation[i].genes >> j) & 1) * w[j];
		}
		if(cPopulation[i].weight > s) cPopulation[i].fitness = 0;					// shouldn't happen since all cases are handled, but just in case
		totalFitness += cPopulation[i].fitness;
	}
	for(int i = 0; i < popsize; i++)									// create CDF
	{
		double probability_of_i = cPopulation[i].fitness / (double) totalFitness;
		cdf[i] = (i > 0 ? cdf[i-1] + probability_of_i : probability_of_i);
	}
}

bool valid(chromosome* pC)
{
	int totalW = 0;
	for(int i=0; i < ngenes; i++)
		totalW += ((pC->genes >> i) & 1) * w[i];
	return totalW <= s;
}

void evolve()
{
	int cInd = 0;
	if(elite)
	{
		pair<int, int> best (0,0), secondBest (0,0);								// elitism, let the two fittest chromosomes survive with every generation	
		for(int i = 0; i < popsize; i++)
		{
			if(cPopulation[i].fitness >= best.second)
			{
				secondBest = best;
				best.first = i;
				best.second = cPopulation[i].fitness;
			}	
			else if(cPopulation[i].fitness > secondBest.second)
			{
				secondBest.first = i;
				secondBest.second = cPopulation[i].fitness;
			}
		}
		nPopulation[cInd++] = cPopulation[best.first];
		nPopulation[cInd++] = cPopulation[secondBest.first];
	}

	// begin natural selection
	
	while(cInd < popsize)
	{
		chromosome p1 = cPopulation[upper_bound(cdf, cdf + popsize, rD) - cdf]; //first parent
		chromosome p2 = cPopulation[upper_bound(cdf, cdf + popsize, rD) - cdf]; //second parent
		chromosome c1 = p1;
		chromosome c2 = p2;
		nPopulation[cInd++] = p1;
		nPopulation[cInd++] = p2;
		bool willCross = (rD <= pCross);
		if(willCross)
		{
			int crossover_point = rand() % ngenes;
			if(crossover_point == 0) crossover_point++;
			else if(crossover_point == ngenes) crossover_point--;
			chromosome c1 = p1;
			chromosome c2 = p2;
			for(int i = crossover_point + 1; i < ngenes; i++)
			{
				if(p1.genes & (1<<(i-1))) 		// if the i'th bit is on in parent 1, it should be on in child 2
					c2.genes |= 1<<(i-1);
				else			 		// if the i'th bit is off in parent 1, it should be off in child 2
					c2.genes &= ~(1<<(i-1));

				if(p2.genes & (1<<(i-1))) 		// if the i'th bit is on in parent 2, it should be on in child 1
					c1.genes |= 1<<(i-1);
				else			 		// if the i'th bit is off in parent 2, it should be off in child 1
					c1.genes &= ~(1<<(i-1));
			}
		}
		if(valid(&c1) && cInd < popsize) nPopulation[cInd++] = c1;
		if(valid(&c2) && cInd < popsize) nPopulation[cInd++] = c2;
	}

	// end natural selection
	// begin random mutation
	
	for(int i = 0; i < popsize; i++)
	{
		if(rD <= pMutate)
		{
			int flip = rand() % ngenes;			// random index to mutate (flip)
			chromosome mutated = nPopulation[i];
			mutated.genes ^= (1<<flip);
			if(valid(&mutated)) nPopulation[i] = mutated;
		}
	}

	// end random mutation
}

int main()
{
	int c;
	cin >> c;;
	while(c--)
	{
		cin >> ngenes >> s;
		for(int i = 0; i < ngenes; i++)
			cin >> w[i] >> v[i];
		init();
		for(int i = 0; i < ngenerations; i++)
		{ 
			calcFitness();							// calculates fitness for the current generation and creates the CDF table
			evolve();							// natural selection, crossover, and random mutation
			memcpy(cPopulation, nPopulation, sizeof(nPopulation));		// replace current population with the new population
		}
		int globalOptimum = 0, gOInd = 0;
		calcFitness();
		for(int i = 0; i < popsize; i++)
		{
			if(cPopulation[i].fitness > globalOptimum)
			{
				globalOptimum = cPopulation[i].fitness;
				gOInd = i;
			}
		}
		for(int i = 0; i < ngenes; i++)
			cout << ((cPopulation[gOInd].genes >> i) & 1);
		cout << " = " << globalOptimum << endl;
	
	}
	return 0;
}
