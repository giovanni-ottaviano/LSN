//LSN_Exercise_10
//Exercise 10.1 - Simulated anneling for TSP

#include <iostream>
#include <fstream>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <iomanip>
#include "esercizio10_1.h"

using namespace std;


int main() {
    const int N_schedule = 500;  //Number of different beta
    const double T_min = 0.0025, T_max = 10.;   //Fictious temperatures
    std::vector<city> cities;
    path Path;

    std::vector<double> v_beta(N_schedule);     //Inverse temperature (anneling schedule)
    std::vector<int> v_N(N_schedule);           //Steps per beta (anneling schedule)

    std::vector<double> lengths;                //Vector for saving paths lenghts

    Random rnd;
    init_random_gen(rnd);

//Initialize vector of cities and path
    initialize_circ(rnd, cities, Path, 32);
//    initialize_square(rnd, cities, Path, 32);

//inizializza i vettori di scheduling
    for (int i = 0; i < N_schedule; i++) {
        v_beta[i] = 1./T_max - (1./T_max - 1./T_min)*i/N_schedule;
        v_N[i] = (v_beta[i] <= 50) ? 100 : 1000;
    }

//run simulation (on anneling schedule)
    for (int isched = 0; isched < N_schedule; isched++) {
        double beta = v_beta[isched];
        int N = v_N[isched];

        for (int step = 0; step < N; step++) {
            metropolis(rnd, cities, Path, beta);
            lengths.push_back(Path.get_fitness());
        }

        cout << "Ended step " << isched + 1 << ". Path length: " << lengths.back() << " beta: " << std::setprecision(8) << beta << endl;
    }

//Print results to file (lengths and best path)
    print_vector("lengths_circle_SA.dat", lengths);
    print_best_path("path_circle_SA.dat", Path, cities);
//    print_vector("lengths_square_SA.dat", lengths, 8);
//    print_best_path("path_square_SA.dat", Path, cities);


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


/******************************************************************************/
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

//Fill given vectors of cities and path with the desire number of elements on a circumference
void initialize_circ(Random& rnd, std::vector<city>& cities, path& Path, unsigned int N_cities, double r) {
    std::vector<city> v_cs(N_cities);
    std::vector<unsigned int> permutations(N_cities);


//make new path
    std::iota(permutations.begin(), permutations.end(), 0);
    std::random_shuffle(permutations.begin() + 1, permutations.end());

//make 1st path via random generating cities
    double theta = 0.;
    do {
        for (unsigned int j = 0; j < N_cities; j++) {
            theta = rnd.Rannyu(0., 2*M_PI);
            v_cs[j].set_x(r*cos(theta));
            v_cs[j].set_y(r*sin(theta));
        }
    } while (!check_cities(v_cs));


//check if the new population fulfills the bonds
    if (!fulfill_bonds(v_cs, path(permutations))) {
        cerr << endl << "ERROR: initial components does not fulfill bonds." << endl;
        return;
    }

//Move new path and vector of cities
    Path = std::move(path(permutations));
    cities = std::move(v_cs);
}

//Fill given vectors of cities and path with the desire number of elements in a square
void initialize_square(Random& rnd, std::vector<city>& cities, path& Path, unsigned int N_cities, double l) {
    std::vector<city> v_cs(N_cities);
    std::vector<unsigned int> permutations(N_cities);

//fill vector of permutations with increasing numbers
    std::iota(permutations.begin(), permutations.end(), 0);
    std::random_shuffle(permutations.begin() + 1, permutations.end());

//make 1st path via random generating cicies
    do {
        for (unsigned int j = 0; j < N_cities; j++) {
            v_cs[j].set_x(rnd.Rannyu(-l, l));
            v_cs[j].set_y(rnd.Rannyu(-l, l));
        }
    } while (!check_cities(v_cs));


//check if the new population fulfills the bonds
    if (!fulfill_bonds(v_cs, path(permutations))) {
        cerr << endl << "ERROR: initial components does not fulfill bonds." << endl;
        return;
    }

//Move new path and vector of cities
    Path = path(std::move(permutations));
    cities = std::move(v_cs);
}

//check if path fulfills the bonds -> forse fatto così non è il modo più efficiente
bool fulfill_bonds(const std::vector<city>& v_cs, const path& Path) {
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
void init_random_gen(Random& gen) {
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


//Metropolis algorithm
void metropolis(Random& rnd, const std::vector<city>& v_cs, path& Path, const double& beta) {
    path new_Path = Path;
    double r = rnd.Rannyu();    //select a mutation

//Trial move on new_Path (using mutations)
    if (r < 0.25)
        new_Path.pair_permutation(rnd);
    else if (r >= 0.25 && r < 0.5)
        new_Path.shift(rnd);
    else if (r >= 0.5 && r < 0.75)
        new_Path.swap_range(rnd);
    else
        new_Path.inversion(rnd);

//Compute both fitness (new/old)
    double fitness_delta = new_Path.path_length_L1(v_cs) - Path.path_length_L1(v_cs);

//Verify the move
    if (fitness_delta <= 0.)
        Path = new_Path;
    else if (rnd.Rannyu() <= exp(-beta*fitness_delta))
        Path = new_Path;
}

//Print single vector
template <typename T>
void print_vector(const std::string& name, const std::vector<T>& v1, int prec) {
    ofstream fileout(name);

    if (fileout.is_open()) {
        for (T elem : v1)
            fileout << std::setprecision(prec) << std::fixed << elem << endl;
    }
    else {
        cerr << "ERROR: Unable to write file " << name << endl;
    }

    fileout.close();
}


//Print 1st path
void print_best_path(const std::string& fname, const path& Path, const std::vector<city>& cities) {
    ofstream fileout(fname);

    if (fileout.is_open()) {
        unsigned int index = 0;

        for (unsigned int i = 0; i < cities.size(); i++) {
            index = Path.get_city(i);
            fileout << cities[index].get_x() << " " << cities[index].get_y() << endl;
        }
    }
    else {
        cerr << endl << "ERROR: Unable to open file: " << fname << endl;
    }

    fileout.close();
}
