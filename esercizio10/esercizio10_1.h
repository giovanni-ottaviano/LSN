//LSN_Exercise_10
//Exercise 10.1

#ifndef _TSP_SA_h_
#define _TSP_SA_h_

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
//index of first city is always fixed to 0
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

        void pair_permutation(Random&);
        void shift(Random&);
        void swap_range(Random&);
        void inversion(Random&);

    private:
        std::vector<unsigned int> _cities_index;
        double _fitness;
};


//Generic functions
bool check_cities(std::vector<city>);
void initialize_circ(Random&, std::vector<city>&, path&, unsigned int, double r = 1.);
void initialize_square(Random&, std::vector<city>&, path&, unsigned int, double l = 1.);
bool fulfill_bonds(const std::vector<city>&, const path&);
void metropolis(Random&, const std::vector<city>&, path&, const double&);
void init_random_gen(Random&);
void print_best_path(const std::string&, const path&, const std::vector<city>&);
template <typename T> void print_vector(const std::string&, const std::vector<T>&, int prec = 5);

#endif
