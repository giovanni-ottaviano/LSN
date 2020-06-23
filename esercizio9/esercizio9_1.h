//LSN_Exercise_09
//Exercise 09.1

#ifndef _TSP_h_
#define _TSP_h_

#include <vector>
#include <string>
#include "Parallel_Generator/random.h"


class city {

    public:
        city(double, double);   //CTOR
        city(const city&);      //Copy CTOR
        ~city() {};             //DTOR

        void set_x(const double& x) {_x = x;}
        void set_y(const double& y) {_y = y;}
        double get_x() const {return _x;}
        double get_y() const {return _y;}

        double distance(const city&) const;
        double distance2(const city&) const;

        void print_city() const;

        bool operator==(const city&) const;
        bool operator!=(const city&) const;

    private:
        double _x, _y;
};


//Individual for genetic algorithm
//index of first city is always fixed to 0 (decrease degeneration of the best path)
class path {

    public:
        path() {}           //CTOR
        path(const path&);  //copy CTOR
        path(const std::vector<unsigned int>&);
        path(std::vector<unsigned int>&&);
        ~path() {}

        void set_cities(const std::vector<unsigned int>& cs) {_cities_index = cs;}
        void set_cities(std::vector<unsigned int>&&);
        unsigned int get_size() const {return _cities_index.size();}
        std::vector<unsigned int>& get_cities() {return _cities_index;}
        const std::vector<unsigned int>& get_cities() const {return _cities_index;}
        unsigned int get_city(int i) const {return _cities_index.at(i);}
        double get_fitness() const {return _fitness;}

        double path_length_L1(const std::vector<city>&);
        double path_length_L2(const std::vector<city>&);

        void print_path() const;
        bool operator<(const path&) const;
        bool operator>(const path&) const;

        friend class population;    //population can access private members of path


    private:
        std::vector<unsigned int> _cities_index;
        double _fitness;

        void pair_permutation(Random&);
        void shift(Random&);
        void swap_range(Random&);
        void inversion(Random&);
};


class population {

    public:
        population() : _rnd(nullptr) {}
        population(const population& p) : _ordered_cities(p._ordered_cities), _chroms(p._chroms), _rnd(p._rnd) {}
        population(const std::vector<city>&, const std::vector<path>&, Random*);
        population(std::vector<city>&&, std::vector<path>&&, Random*);
        population(population&&);
        ~population() {}

        void set_population(const population&);
        void set_population(population&&);
        void set_ordered_cities(const std::vector<city>& cities) {_ordered_cities = cities;}
        void set_chromosome(const path&, unsigned int);
        void set_chromosome(path&&, unsigned int);
        void set_all_chromosomes(const std::vector<path>&);
        void set_all_chromosomes(std::vector<path>&&);
        void set_generator(Random* rnd) {_rnd = rnd;}
        std::vector<city>& get_ordered_cities() {return _ordered_cities;};
        path& get_chromosome(unsigned int);
        std::vector<path>& get_all_chromosomes() {return _chroms;}
        Random* get_generator_pointer() {return _rnd;}

        void print_best_path(const std::string&) const;

        void sort_on_fitness();
        void prepare_upgrade();
        void upgrade_generation();
        void compute_fitness_L1();
        void compute_fitness_L2();
        void best_fitness(const std::string&) const;
        void best_half_fitness(const std::string&) const;

    private:
        std::vector<city> _ordered_cities;
        std::vector<path> _chroms, _new_chroms;
        Random* _rnd;

        double _p_cross = 0.55;
        double _p_pair = 0.15;
        double _p_shift = 0.15;
        double _p_srange = 0.15;
        double _p_inv = 0.15;
        double _p_exp = 4.;     //convenient exponent for rigged wheel (selection)
        int _selected1 = 0, _selected2 = 0; //individuals selected by method selection()
        unsigned int _size_new = 0;

        bool compare_paths_L1(path, path);
        bool compare_paths_L2(path, path);

        void selection();
        void crossover();
        void mutations();
        unsigned int posMax(unsigned int) const;
        void selectionSort();
};

//Generic functions
bool checke_cities(std::vector<city>);
population make_population(Random&, unsigned int, unsigned int);
population make_population_circ(Random&, unsigned int, unsigned int, double);
population make_population_square(Random&, unsigned int, unsigned int, double);
bool fulfill_bonds(const std::vector<city>&, const path&);
void init_random_gen(Random&);

#endif
