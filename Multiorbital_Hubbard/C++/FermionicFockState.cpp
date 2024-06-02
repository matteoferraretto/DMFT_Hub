#include <FermionicFockState.h>
#include <iostream>
#include <cmath>

// constructor
FermionicFockState::FermionicFockState(int L, int Norb, int f){
    this->L = L; this->Norb = Norb; this->f = f;
    // state is a rank 2 tensor of dimension L x (f*Norb)
    state = new bool* [f*Norb];
    for(int i=0; i<f*Norb; i++){
        state[i] = new bool [L];
    }
    for(int i=0; i<f*Norb; i++){
        for(int j=0; j<L; j++){
            state[i][j] = 0;
        }
    }
};

// methods
int* FermionicFockState::get_size(){
    static int dim[2];
    dim[0] = f*Norb; dim[1] = L;
    return dim;
}
int* FermionicFockState::get_size(bool printQ){
    static int dim[2];
    dim[0] = f*Norb; dim[1] = L;
    if(printQ) std::cout << "size of state: (" << dim[0] << ", " << dim[1] << ")\n";
    return dim;
}

bool FermionicFockState::get_element(int i, int orb, int sigma){
    return state[f*orb+sigma][i];
}

void FermionicFockState::print_state(){
    for(int i=0; i<f*Norb; i++){
        for(int j=0; j<L; j++){
            std::cout << state[i][j] << " ";
        }
        std::cout << "\n";
    }
}

void FermionicFockState::get_memory(){
    std::cout << "Memory required: " << sizeof(state) * L * Norb * f << " bytes.\n";
}

void FermionicFockState::bit_set(int i, int orb, int sigma){
    state[f*orb+sigma][i] = 1;
}

int FermionicFockState::count_fermions(int site, int orb, int sigma){
    int count = 0;
    int index = 0;
    int max_index = (f*orb + sigma)*L + site;
    if(max_index == 0){ return 0; }
    for(int i=0; i<f*Norb; i++){
        for(int j=0; j<L; j++){
            count += state[i][j];
            index++;
            if(index > max_index) break;
        }
        if(index > max_index) break;
    }
    return count;
}

// destructor
FermionicFockState::~FermionicFockState(){
    for(int i=0; i<f*Norb; i++){
        delete [] this->state[i];
    }
    delete [] this->state; this->state = NULL;
};