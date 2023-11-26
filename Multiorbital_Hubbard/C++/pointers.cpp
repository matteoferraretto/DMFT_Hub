#include<iostream>
#include<cmath>

class FermionicFockState{
private:
    bool** state; // state

public:
    int L; // number of sites in the impurity model
    int Norb; // number of orbitals 
    int f; // number of flavors

// constructor
    FermionicFockState(int L, int Norb, int f){
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
    int* get_size(){
        static int dim[2];
        dim[0] = f*Norb; dim[1] = L;
        return dim;
    }
    int* get_size(bool printQ){
        static int dim[2];
        dim[0] = f*Norb; dim[1] = L;
        if(printQ) std::cout << "size of state: (" << dim[0] << ", " << dim[1] << ")\n";
        return dim;
    }

    bool get_element(int i, int orb, int sigma){
        return state[f*orb+sigma][i];
    }

    void print_state(){
        for(int i=0; i<f*Norb; i++){
            for(int j=0; j<L; j++){
                std::cout << state[i][j] << " ";
            }
            std::cout << "\n";
        }
    }

    void get_memory(){
        std::cout << "Memory required: " << sizeof(state) * L * Norb * f << " bytes.\n";
    }

    void bit_set(int i, int orb, int sigma){
        state[f*orb+sigma][i] = 1;
    }

    int count_fermions(int site, int orb, int sigma){
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
    ~FermionicFockState(){
        for(int i=0; i<f*Norb; i++){
            delete [] this->state[i];
        }
        delete [] this->state; this->state = NULL;
    };
};


// Application 1: function with multiple outputs
void Min_and_Max(double array[], int size, double* min, double* max){
    for(int i=0; i < size; i++){
        if(array[i] < *min){ *min = array[i]; }
        if(array[i] > *max){ *max = array[i]; }
    }
}

int main(){

    std::cout << "state:\n";
    FermionicFockState psi(3, 2, 1);
    //show content of psi
    psi.get_size(true);
    std::cout << "\n";
    /* set a bit */
    psi.bit_set(2,0,0);
    std::cout << psi.get_element(2,0,0) << "\n";
    psi.print_state();
    psi.get_memory();
    std::cout << psi.count_fermions(2,0,0) << "\n";
    std::cout << psi.count_fermions(0,1,0) << "\n";

    return 0;
}