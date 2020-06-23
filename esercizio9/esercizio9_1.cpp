//LSN_Exercise_09
//Exercise 09.1

#include <iostream>
#include <fstream>
#include <cmath>
#include <numeric>
#include <algorithm>
#include "esercizio9_1.h"

using namespace std;


int main() {
    const int max_steps = 10000;
    Random rnd;
    init_random_gen(rnd);

//Initialize population
    population p(make_population_circ(rnd, 100, 32, 1.));
    p.set_generator(&rnd);  //Set random generator

//Run simulation
    const int istop = 10000;    //Here, run until the end
    int iprint = 0, stop = 0;
    double old_fitness = 0., new_fitness = 0.;

    p.prepare_upgrade();
    p.best_fitness("fitness_circle.dat");
    p.best_half_fitness("ave_L1_circle.dat");

    for (int i = 0; i < max_steps; i++) {
        p.upgrade_generation();
        p.prepare_upgrade();    //prepare for next iteration
        p.best_half_fitness("ave_L1_circle.dat");  //Write ave fitness for best half population on file
        p.best_fitness("fitness_circle.dat");

        new_fitness = p.get_chromosome(0).get_fitness();

        if (iprint % 500 == 0)
            cout << new_fitness << endl;

        if (old_fitness == new_fitness) {
            stop++;
        }
        else {
            old_fitness = new_fitness;
            stop = 0;
        }

        if (stop == istop) {    //exit after istops times new_fitness == old_fitness
            cout << "Stop at: " << i << endl;
            break;
        }

        iprint++;
    }

    cout << endl << "Best fitness: " << p.get_chromosome(0).get_fitness() << endl;
    p.print_best_path("path_circle.dat");

    return 0;
}



//CLASS **city**
city::city (double x = 0., double y = 0.) : _x(x), _y(y) {}

city::city (const city& c) : _x(c._x), _y(c._y) {}

//distance between 2 cities
double city::distance(const city& c) const {

    return sqrt(distance2(c));
}

//distance squared
double city::distance2(const city& c) const {

    return pow(_x - c.get_x(), 2) + pow(_y - c.get_y(), 2);
}

void city::print_city() const {

    cout << _x << " " << _y << endl;
}

//Operator overload
bool city::operator==(const city& c) const {

    return ( (_x == c._x) && (_y == c._y) );
}

bool city::operator!=(const city& c) const {

    return ( (_x != c._x) || (_y != c._y) );
}


//CLASS **path**
path::path(const std::vector<unsigned int>& cs) : _cities_index(cs), _fitness(0.) {}

path::path(std::vector<unsigned int>&& cs) {

    _cities_index = std::move(cs);
     _fitness = 0.;
}

path::path(const path& p) : _cities_index(p._cities_index), _fitness(0.) {}


void path::set_cities(std::vector<unsigned int>&& cs) {

    _cities_index = std::move(cs);
}


//Path length using L(1) norm -> loss function
double path::path_length_L1(const std::vector<city>& pos_cities) {

    if (pos_cities.size() != _cities_index.size()) {
        cerr << "ERROR: Number of cities is different form Number of indexes!" << endl;
        return 0.;
    }

    double length = 0.;
    unsigned int index1 = 0, index2 = 0;

    for (unsigned int i = 0; i < _cities_index.size() - 1; i++) {
        index1 = _cities_index[i];
        index2 = _cities_index[i+1];
        length += pos_cities[index1].distance(pos_cities[index2]);
    }

    index1 = _cities_index.back();  //last index
    length += pos_cities[index1].distance(pos_cities[0]);    //come back to 1st city

//set fitness and return it
    _fitness = length;

    return length;
}

//Path length using L(2) norm -> loss function
double path::path_length_L2(const std::vector<city>& pos_cities) {

    if (pos_cities.size() != _cities_index.size()) {
        cerr << "ERROR: Number of cities is different form Number of indexes!" << endl;
        return 0.;
    }

    double length = 0.;
    unsigned int index1 = 0, index2 = 0;

    for (unsigned int i = 0; i < _cities_index.size() - 1; i++) {
        index1 = _cities_index[i];
        index2 = _cities_index[i+1];
        length += pos_cities[index1].distance2(pos_cities[index2]);
    }

    index1 = _cities_index.back();  //last index
    length += pos_cities[index1].distance2(pos_cities[0]);    //come back to 1st city

//set fitness and return it
    _fitness = length;

    return length;
}


//Mutation: swap 2 cities (except for the 1st)
void path::pair_permutation(Random& rnd) {
    int index1 = int(rnd.Rannyu(1., _cities_index.size()));
    int index2;

    do {
        index2 = int(rnd.Rannyu(1., _cities_index.size()));
    } while(index2 == index1);

//permutate selected cities
    std::swap(_cities_index[index1], _cities_index[index2]);
}

//Mutation: Shift m contiguos cities (except for the 1st) of +n pos
void path::shift(Random& rnd) {
    unsigned int size = _cities_index.size();
    int index1 = int(rnd.Rannyu(1., size - 1.));
    int index2, index_min;

    do {
        index2 = int(rnd.Rannyu(1., size - 1.));
    } while(index2 == index1);

//find min and max and shift
    index_min = std::min(index1, index2);
    index2 = std::max(index1, index2);
    int n = ((size - 1 - index2) == 1) ? 1 : int(rnd.Rannyu(1., size - index2));

//execute the shift
    std::vector<unsigned int> appo_path; //vector with n deleted elements
    for (int i = index2 + 1; i <= index2 + n ; i++)
        appo_path.push_back(_cities_index[i]);

//shifted values
    std::move_backward(_cities_index.begin() + index_min, _cities_index.begin() + index2 + 1,
                       _cities_index.begin() + index2 + n + 1);

//add overwritten values
    auto appo_iter = appo_path.begin();
    for (int i = index_min; i < index_min + n ; i++) {
        _cities_index[i] = *appo_iter;
        ++appo_iter;
    }
}

//Mutation: Swap m contiguos cities with other m different cities
void path::swap_range(Random& rnd) {
    unsigned int size = _cities_index.size();    //need _cities_index.size() > 5
    double even = (size % 2 == 0) ? 0.5 : 0.;    //needed for bilancing probability

    int m = int(rnd.Rannyu(1., (size - 1.)/2. + even)); //m <= (N-1)/2
    int index2 = int(rnd.Rannyu(m, size - m));

    int index_min = index2 - m + 1;

//swap the selected range (swap elem in range [first1, last1) with [first2,...)
    std::swap_ranges(_cities_index.begin() + index_min, _cities_index.begin() + index2 + 1,
                     _cities_index.begin() + index2 + 1);
}

//Mutation: reverse the order in which m cities appears(except for the 1st)
void path::inversion(Random& rnd) {
    int index1 = int(rnd.Rannyu(1., _cities_index.size()));  //m <= N-1
    int index2, index_min;

    do {
        index2 = int(rnd.Rannyu(1., _cities_index.size()));
    } while(index2 == index1);

//find min and max
    index_min = std::min(index1, index2);
    index2 = std::max(index1, index2);

//reverse interval (+1 'cos need to include element[index2])
    std::reverse(_cities_index.begin() + index_min, _cities_index.begin() + index2 + 1);
}

void path::print_path() const {

    for(unsigned int elem : _cities_index)
        cout << elem << " ";

    cout << endl;
}

//Operator overload
bool path::operator<(const path& p) const {

    return _fitness < p._fitness;
}

bool path::operator>(const path& p) const {

    return _fitness > p._fitness;
}



//CLASS **POPULATION**
population::population(population&& p) {

    _ordered_cities = std::move(p._ordered_cities);
    _chroms = std::move(p._chroms);
    _rnd = p._rnd;
}

population::population(const std::vector<city>& v_cs, const std::vector<path>& v_ph, Random* rnd = nullptr)
                      : _ordered_cities(v_cs), _chroms(v_ph), _rnd(rnd) {}

population::population(std::vector<city>&& v_cs, std::vector<path>&& v_ph, Random* rnd = nullptr) {

     _ordered_cities = std::move(v_cs);
     _chroms = std::move(v_ph);
     _rnd = rnd;
}

void population::set_population(const population& p) {

    _ordered_cities = p._ordered_cities;
    _chroms = p._chroms;
    _rnd = p._rnd;
}

void population::set_population(population&& p) {

    _ordered_cities = std::move(p._ordered_cities);
    _chroms = std::move(p._chroms);
    _rnd = p._rnd;
}

//Set a chromosome at position pos (if inside the array)
void population::set_chromosome(path&& chrom, unsigned int pos) {

    if (pos >= _chroms.size()) {
        cerr << "Tryng to insert a CHROMOSOME over the population size!" << endl;
        cerr << "Population size: " << _chroms.size() << endl;
        cerr << "Desired position: " << pos << endl;
    }
    else if (chrom.get_size() != (_chroms[0]).get_size()) {
        cerr << endl << "Trying to insert a wrong CHROMOSOME! Two sizes are different." << endl;
    }
    else {  //all good
        _chroms[pos] = std::move(chrom);
    }
}

void population::set_chromosome(const path& chrom, unsigned int pos) {

    if (pos >= _chroms.size()) {
        cerr << endl << "Tryng to insert a CHROMOSOME over the population size!" << endl;
        cerr << "Population size: " << _chroms.size() << endl;
        cerr << "Desired position: " << pos << endl;
    }
    else if (chrom.get_size() != (_chroms[0]).get_size()) {
        cerr << endl << "Trying to insert a wrong CHROMOSOME! Two sizes are different." << endl;
    }
    else {  //all good
        _chroms[pos] = chrom;
    }
}

void population::set_all_chromosomes(const std::vector<path>& chroms) {
    unsigned int right_size = _chroms[0].get_size();

    for (unsigned int i = 0; i < chroms.size(); i++) {
        if (chroms[i].get_size() != right_size) {
            cerr << endl << "ERROR: Unable to substitute population, because sizes are different!" << endl;
            return;
        }
    }

    _chroms = chroms;    //if all good
}

void population::set_all_chromosomes(std::vector<path>&& chroms) {
    unsigned int right_size = _chroms[0].get_size();

    for (unsigned int i = 0; i < chroms.size(); i++) {
        if (chroms[i].get_size() != right_size) {
            cerr << endl << "ERROR: Unable to substitute population, because sizes are different!" << endl;
            return;
        }
    }

    _chroms = std::move(chroms);    //if all good
}

path& population::get_chromosome(unsigned int pos) {

    if (pos >= _chroms.size())
        cerr << endl << "ERROR: chromosome at position " << pos << " does not exist." << endl;

    return _chroms.at(pos);
}


//Sort vector of path on fitness L(1) or L(2)
void population::sort_on_fitness() {

//If sort_on_fitness is public decomment one of the following lines
//    compute_fitness_L1();
//    compute_fitness_L2();
    selectionSort();
}


//Method of selection for individuals
void population::selection() {

    _selected1 = int(_chroms.size()*pow(_rnd->Rannyu(), _p_exp));
    _selected2 = int(_chroms.size()*pow(_rnd->Rannyu(), _p_exp));
}


//Crossover operator
void population::crossover() {
    unsigned int size = _ordered_cities.size();
    int icut = int(_rnd->Rannyu(1., size - 1)); // cut index

//generate 2 vectors with new paths, copying the 1st part and crossing the second
    std::vector<std::vector<unsigned int>> new_chroms(2, std::vector<unsigned int>(icut));
    for (int i = 0; i < icut; i++) {
        new_chroms[0][i] = _chroms[_selected1].get_city(i);
        new_chroms[1][i] = _chroms[_selected2].get_city(i);
    }

//add remaining elements in the correct order
    int ichr = 1;
    for (int index : {_selected1, _selected2}) {
        for (unsigned int elem : _chroms[index].get_cities()) {
            if (new_chroms[ichr].size() == size)  //if v is full, end for loop
                break;

            if ( std::find(new_chroms[ichr].begin(), new_chroms[ichr].end(), elem) == new_chroms[ichr].end() )
                new_chroms[ichr].push_back(elem);
        }

        ichr--;
    }

//add new chromosomes to appo population
    _new_chroms.push_back(std::move(new_chroms[0]));
    _new_chroms.push_back(std::move(new_chroms[1]));
}

//Method that tries sorts the paths and then to execute all mutation and the crossover
void population::mutations() {

    if (_rnd->Rannyu() <= _p_cross) {  //Try the crossover and check bonds
        crossover();
        fulfill_bonds(_ordered_cities, _new_chroms[_size_new]);
        fulfill_bonds(_ordered_cities, _new_chroms[_size_new+1]);
    }
    else {
        _new_chroms.push_back(_chroms[_selected1]);     //if no crossover happens, add parents to new pop and try other mutations
        _new_chroms.push_back(_chroms[_selected2]);
    }

//Try all other mutations
    for (int i : {_size_new, _size_new + 1}) {
        if (_rnd->Rannyu() <= _p_pair) {     //pair permutation
            _new_chroms[i].pair_permutation(*_rnd);
            fulfill_bonds(_ordered_cities, _new_chroms[i]);
        }

        if (_rnd->Rannyu() <= _p_shift) {
            _new_chroms[i].shift(*_rnd);
            fulfill_bonds(_ordered_cities, _new_chroms[i]);
        }

        if (_rnd->Rannyu() <= _p_srange) {     //swap range
            _new_chroms[i].swap_range(*_rnd);
            fulfill_bonds(_ordered_cities, _new_chroms[i]);
        }

        if (_rnd->Rannyu() <= _p_inv) {         //inversion
            _new_chroms[i].inversion(*_rnd);
            fulfill_bonds(_ordered_cities, _new_chroms[i]);
        }
    }

    _size_new += 2;     //update new vector size
}

//Initialize fitness and sort _chorms
void population::prepare_upgrade() {

    if (_chroms[0].get_fitness() == 0)
        compute_fitness_L1();
//        compute_fitness_L2();

    sort_on_fitness();  //Sorting the actual generation on fitness
}


//Method for upgrading the actual generation to next one
void population::upgrade_generation() {
    unsigned int old_chroms_size = _chroms.size();

//Create new generation
    while (_size_new != old_chroms_size) {   //repeat upgrade until new generatio is created
            selection();
            mutations();
    }

//when the new population is ready, substitute the older one
    if (old_chroms_size != _size_new) {
        cerr << "WARNING: new population size is different." << endl;
        cerr << "Old size: " << old_chroms_size << " New size: " << _size_new << endl;
    }

    _chroms = std::move(_new_chroms);
    _size_new = 0;

    if (!_new_chroms.empty())
        cout << "Not empty!!" << "  " << _new_chroms.size() << endl << endl;
}

//Function for setting the fitness of all paths
void population::compute_fitness_L1() {

    for (auto& P : _chroms)
        P.path_length_L1(_ordered_cities);
}

//Function for setting the fitness of all paths
void population::compute_fitness_L2() {

    for (auto& P : _chroms)
        P.path_length_L2(_ordered_cities);
}

//Print 1st path
void population::print_best_path(const std::string& fname) const {
    ofstream fileout(fname);

    if (fileout.is_open()) {
        unsigned int index = 0;

        for (unsigned int i = 0; i < _ordered_cities.size(); i++) {
            index = _chroms[0].get_city(i);
            fileout << _ordered_cities[index].get_x() << " " << _ordered_cities[index].get_y() << endl;
        }
    }
    else {
        cerr << endl << "ERROR: Unable to open file: " << fname << endl;
    }

    fileout.close();
}


//Find the position of the path with maxlength
unsigned int population::posMax(unsigned int N) const {
  	path appo = _chroms[0];
  	int pos = 0;

	for (unsigned int i = 0; i < N; i++){
   		if (appo.get_fitness() < _chroms[i].get_fitness()) {
    		appo = _chroms[i];
    		pos = i;
   		}
 	}

    return pos;
}


void population::selectionSort() {
	path appo;
   	int pm = 0;

	for (int i = _chroms.size() - 1; i > 0; i--) {
        pm = posMax(i + 1);
     	appo = _chroms[i];
    	_chroms[i] = _chroms[pm];
		_chroms[pm] = appo;
   	}
}

//Compute best fitness and print to file (append)
void population::best_fitness(const std::string& fname) const {
    ofstream fileout(fname, ios::app);

    if (fileout.is_open())
        fileout << _chroms[0].get_fitness() << endl;
     else
        cerr << endl << "ERROR: Unable to open file " << fname << endl;

    fileout.close();
}

//Compute average value of fitness for the best half of the population and print to file (append)
void population::best_half_fitness(const std::string& fname) const {
    ofstream fileout(fname, ios::app);
    int half_size = _chroms.size()/2;
    double sum = 0.;

    for (int i = 0; i < half_size; i++)
        sum += _chroms[i].get_fitness();

    if (fileout.is_open())
        fileout << sum/half_size << endl;
    else
        cerr << endl << "ERROR: Unable to open file " << fname << endl;


    fileout.close();
}
/*-------------------------------------------------------------------------------------------------------*/
//GENERIC FUNCTIONS

//Check if cities are in different positions
bool check_cities(std::vector<city> v_cs) {

    if (std::unique(v_cs.begin(), v_cs.end()) != v_cs.end()) {
        cout << "Some cities are in the same position ..." << endl;
        cout << "Generate again." << endl;
        return false;
    }

    return true;
}


//Make a new population and returns it
//Here #chromosomes = #path and #N_genes = #cities for each path
population make_population(Random& rnd, unsigned int N_chrom, unsigned int N_genes) {
    std::vector<path> chromosomes;
    std::vector<city> genes(N_genes);
    std::vector<unsigned int> permutations(N_genes);

//fill vector of permutations with increasing numbers
    std::iota(permutations.begin(), permutations.end(), 0);
    chromosomes.push_back(path(permutations));

//make 1st path via random generating cicies
    do {
        for (unsigned int j = 0; j < N_genes; j++) {
            genes[j].set_x(rnd.Rannyu(-10,10));
            genes[j].set_y(rnd.Rannyu(-10,10));
        }
    } while (!check_cities(genes));

//make permutations and fill chromosomes
    for (unsigned int i = 1; i < N_chrom; i++) {
//        std::next_permutation(permutations.begin() + 1, permutations.end());
        std::random_shuffle(permutations.begin() + 1, permutations.end());
        chromosomes.push_back(path(permutations));
    }

//check if the new population fulfills the bonds
    for (unsigned int i = 0; i < N_chrom; i++) {
        if (!fulfill_bonds(genes, chromosomes[i])) {
            cout << endl << "WARNING: initial population does not fulfill bonds." << endl;
            return population();
        }
    }


    return population(genes, chromosomes);
}


//Make a new poulation (with cities distribute on a circumference of radius r) and returns it
//Here #chromosomes = #path and #N_genes = #cities for each path
population make_population_circ(Random& rnd, unsigned int N_chrom, unsigned int N_genes, double r) {
    std::vector<path> chromosomes;
    std::vector<city> genes(N_genes);
    std::vector<unsigned int> permutations(N_genes);


//fill vector of permutations with increasing numbers
    std::iota(permutations.begin(), permutations.end(), 0);
    chromosomes.push_back(path(permutations));

//make 1st path via random generating cities
    double theta = 0.;
    do {
        for (unsigned int j = 0; j < N_genes; j++) {
            theta = rnd.Rannyu(0., 2*M_PI);
            genes[j].set_x(r*cos(theta));
            genes[j].set_y(r*sin(theta));
        }
    } while (!check_cities(genes));

//make permutations and fill chromosomes
    for (unsigned int i = 1; i < N_chrom; i++) {
//        std::next_permutation(permutations.begin() + 1, permutations.end());
        std::random_shuffle(permutations.begin() + 1, permutations.end());
        chromosomes.push_back(path(permutations));
    }

//check if the new population fulfills the bonds
    for (unsigned int i = 0; i < N_chrom; i++) {
        if (!fulfill_bonds(genes, chromosomes[i])) {
            cout << endl << "WARNING: initial population does not fulfill bonds." << endl;
            return population();
        }
    }

    return population(genes, chromosomes);
}


//Make a new poulation (with cities distribute on in a square of side 2l) and returns it
//Here #chromosomes = #path and #N_genes = #cities for each path
population make_population_square(Random& rnd, unsigned int N_chrom, unsigned int N_genes, double l) {
    std::vector<path> chromosomes;
    std::vector<city> genes(N_genes);
    std::vector<unsigned int> permutations(N_genes);

//fill vector of permutations with increasing numbers
    std::iota(permutations.begin(), permutations.end(), 0);
    chromosomes.push_back(path(permutations));

//make 1st path via random generating cicies
    do {
        for (unsigned int j = 0; j < N_genes; j++) {
            genes[j].set_x(rnd.Rannyu(-l, l));
            genes[j].set_y(rnd.Rannyu(-l, l));
        }
    } while (!check_cities(genes));

//make permutations and fill chromosomes
    for (unsigned int i = 1; i < N_chrom; i++) {
//        std::next_permutation(permutations.begin() + 1, permutations.end());
        std::random_shuffle(permutations.begin() + 1, permutations.end());
        chromosomes.push_back(path(permutations));
    }

//check if the new population fulfills the bonds
    for (unsigned int i = 0; i < N_chrom; i++) {
        if (!fulfill_bonds(genes, chromosomes[i])) {
            cout << endl << "WARNING: initial population does not fulfill bonds." << endl;
            return population();
        }
    }

    return population(genes, chromosomes);
}


//check if path fulfills the bonds -> forse fatto così non è il modo più efficiente
bool fulfill_bonds (const std::vector<city>& v_cs, const path& Path) {
    unsigned int psize = Path.get_size();
    std::vector<unsigned int> general(psize);

//check if v_cs and Path have the same dimention
    if (v_cs.size() != Path.get_size()) {
        cout << "size" << endl;
        return false;
    }
//check if 0 is in the 1st position
    if (Path.get_cities()[0] != 0) {
        cout << "no 0" << endl;
        return false;
    }
                                            
//check if others are OK
    std::iota(general.begin(), general.end(), 0);

    return std::is_permutation(general.begin(), general.end(), (Path.get_cities()).begin());
}

//Function for initializing random generator - if it fails then program is terminated
void init_random_gen (Random& gen) {
    int seed[4];
    int p1, p2;
    bool success = true;

    ifstream Primes("Parallel_Generator/Primes");
    if (Primes.is_open())
       Primes >> p1 >> p2;
    else {
        cerr << "PROBLEM: Unable to open Primes" << endl;
        success = false;
    }
    Primes.close();

    ifstream input("Parallel_Generator/seed.in");
    string property;
    if (input.is_open()) {
       while ( !input.eof() ) {
          input >> property;
          if( property == "RANDOMSEED" ) {
             input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
             gen.SetRandom(seed, p1, p2);
          }
       }
       input.close();
    }
    else {
        cerr << "PROBLEM: Unable to open seed.in" << endl;
        success = false;
    }

    if (!success) {
        cerr << "ERROR: Unable to initialize random generator." << endl;
        cerr << "---STOP---" << endl;
        exit (EXIT_FAILURE);
    }
}
