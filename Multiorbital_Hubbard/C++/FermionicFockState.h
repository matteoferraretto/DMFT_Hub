#pragma once

class FermionicFockState{
private:
    bool** state; // state

public:
    int L; // number of sites in the impurity model
    int Norb; // number of orbitals 
    int f; // number of flavors

// constructors
    FermionicFockState() { };
    FermionicFockState(int L, int Norb, int f);

// methods
    int* get_size();

    int* get_size(bool printQ);

    bool get_element(int i, int orb, int sigma);

    void print_state();

    void get_memory();

    void bit_set(int i, int orb, int sigma);

    int count_fermions(int site, int orb, int sigma);

// destructor
    ~FermionicFockState();
};