#include <iostream>
#include <fstream>
#include <random>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <list>
#include <vector>
#include <map>
using namespace std;

typedef pair<int, int> Pair;
typedef vector<int> List;
typedef List::iterator iter;
typedef map<Pair, bool> Dictionary;

const double TL = 1.0;
const double TR = 2.0;

enum Direction {
    horizontal,
    vertical
 };

struct interaction
{
    double time;
    int location;
    Pair coordinate;
    List adjacent_clocks;
    Direction direction;
    interaction* left;
    interaction* right;
};

List find_adjacent_clocks(int i, int j, Direction direction, int M, int N, int current_index) {
    List lst;
    int total = -1;
    switch (direction) {
      case vertical:
          if(i - 1 >= 0) { //Check existence of the previous line (i - 1)
              if(j - 1 >= 0) { //Check existence of clock at coordinate (i - 1, j - 1)
                  total = current_index - (M + 1);
                  lst.push_back(total);
                  //Compute the index in 1D array
              }

              if(j < M) { //Check existence of clock at coordinate (i - 1, j)
                  total = current_index - M;
                  lst.push_back(total);
                  //Compute the index in 1D array
              }

          }

          if(j - 1 >= 0) { //Check existence of clock at coordinate (i, j - 1)
              //Compute the index in 1D array
              total = current_index - 1;
              lst.push_back(total);
          }

          if(j + 1 <= M) { //Check existence of clock at coordinate (i, j + 1)
              // Compute the index in 1D array
              total = current_index + 1;
              lst.push_back(total);
          }

          if(i + 1 < N) { //Check existence of the next line (i + 1)
            if(j < M) { //Check existence of clock at coordinate (i + 1, j)
                //Compute the index in 1D array
                total = current_index + (M + 1);
                lst.push_back(total);
            }

            if(j - 1 >= 0) { //Check existence of clock at coordinate (i + 1, j - 1)
                //Compute the index in 1D array
                total = current_index + M;
                lst.push_back(total);
            }
          }
          break;
      case horizontal:
          if(i - 2 >= 0) { //Check existence of the previous line (i - 2) Previous horizontal line
              //Compute the index of (i - 2, j) in 1D array
              total = current_index - (2 * M + 1);
              lst.push_back(total);
          }

          if(i - 1 >= 0) { //Check existence of the previous line (i - 1)
              if(j >= 0) { //Check existence of clock at coordinate (i - 1, j) // It might always exist
                  //Compute the index in 1D array
                  total = current_index - (M + 1);
                  lst.push_back(total);
              }

              if(j + 1 < M) { //Check existence of clock at coordinate (i - 1, j + 1) It might always exist
                  //Compute the index in 1D array
                  total = current_index - M;
                  lst.push_back(total);
              }
          }

          if(i + 1 >= 0) { //Check existence of the next line (i + 1)
              if(j >= 0) { //Check existence of clock at coordinate (i + 1, j) // It might always exist
                  //Compute the index in 1D array
                  total = current_index + M;
                  lst.push_back(total);
              }

              if(j + 1 < M) { //Check existence of clock at coordinate (i + 1, j + 1) It might always exist
                  //Compute the index in 1D array
                  total = current_index + (M + 1);
                  lst.push_back(total);
              }
          }

          if(i + 2 >= 0) { //Check existence of the previous line (i + 2) Previous horizontal line
              //Compute the index of (i + 2, j) in 1D array
              total = current_index + (2 * M + 1);
              lst.push_back(total);
          }
          break;
    }
    //lst.push_back(10);
    return lst;
}

double energy_summation(double **energy_array, int x, int y, Direction direction) {
    double ret = 10.0;
    switch (direction) {
      case vertical:
        ret = sqrt(energy_array[x / 2][y] + energy_array[x / 2][y+1]);
        break;
      case horizontal:
        ret = sqrt(energy_array[int(x / 2)][y] + energy_array[int(x / 2) + 1][y]);
        break;
    }
    return 10;
}

void print_v(double* Array, int size)
{
    for(int i = 0; i < size; i++)
        cout<< Array[i] << "  ";
    cout<<endl;
}

int find_min(interaction* Link)
{
        double tmp = Link->time;
        interaction* pt = Link;
        interaction* tmp_pt = Link;
        while(tmp_pt != NULL)
        {
            if(tmp_pt->time < tmp)
            {
                tmp = tmp_pt->time;
                pt = tmp_pt;
            }
            tmp_pt = tmp_pt->right;
        }

        return pt->location;

}

void remove(interaction** Link, interaction* pt)
//remove pt from link but keep pt for future use
{
    if( *Link == NULL || pt == NULL)
        return;
    else if( pt->left == NULL && pt->right == NULL)
    {
        *Link = NULL;
        return;
    }
    else
    {
        if((*Link) == pt )
            *Link = pt->right;
        if( pt->right != NULL)
            pt->right->left = pt->left;
        if( pt->left != NULL)
            pt->left->right = pt->right;
    }
    pt->left = NULL;
    pt->right = NULL;
}

void push_front(interaction** Link, interaction* pt)
//push the interaction pointed by pt into the front of the list
{
    pt->left = NULL;
    if(*Link == NULL)
    {
        pt->right = NULL;
        *Link = pt;
    }
    else
    {
        //cout<<"ok";
        pt->right = *Link;
        //if((*Link)->left == NULL) cout<<"NULL"<<endl;
        (*Link)->left = pt;
         *Link = pt;
    }
}

void print_list(interaction* Link)
{
    interaction* tmp = Link;
    while(tmp!= NULL)
    {
        //cout<<" location: "<< tmp->location <<" time: " << tmp->time << "  " ;
        tmp = tmp->right;
    }
    cout<<endl;
}


void big_step_distribute(interaction** &clock_time_in_step, interaction* time_array, const int N, const double small_tau, const int ratio, const int Step)
//distribute clock times of a big step into vectors that represent small steps. If the clock time is bigger than a big tau, then it is arranged in the right location
{
    for(int i = 0; i < N; i++)
    {
        int tmp;
        if(time_array[i].time > (Step + 1)*ratio*small_tau)
        {
            tmp = ratio;
        }
        else
        {
            tmp = int((time_array[i].time - ratio*small_tau*Step)/small_tau);
        }

        if(tmp < 0 || tmp > ratio) {
            tmp = ratio;
        }

        push_front(&clock_time_in_step[tmp], &time_array[i] );
        //cout<<tmp<<endl;


    }
}



void move_interaction(interaction** &clock_time_in_step, interaction* pt, const double small_tau, const int ratio, const int Step, const double new_time)
//move the interaction pointed by *pt from old bucket to new bucket
//n_move1: move without relinking pointers
//n_move2: move with relinking pointers
{
    double old_time = pt->time;
    int old_level, new_level;
    pt->time = new_time;
    if (old_time > (Step + 1)*ratio*small_tau )
    {
        old_level = ratio;
    }
    else
    {
        old_level = int( (old_time - Step*ratio*small_tau)/small_tau );
    }

    if ( new_time > (Step + 1)*ratio*small_tau )
    {
        new_level = ratio;
    }
    else
    {
        new_level = int( (new_time - Step*ratio*small_tau)/small_tau );
    }
    //    cout<<"start to move "<< pt->location << " from " << old_level << " to " << new_level<<endl;
    if(old_level == new_level )
    {
        pt->time = new_time;
    }
    else
    {
        remove(&(clock_time_in_step[old_level]), pt);
        push_front(&(clock_time_in_step[new_level]), pt);
    }
}



void update(interaction** &clock_time_in_step, const int level, const int N, const int M, const double small_tau, const int ratio, const int Step, interaction* time_array, double **energy_array, uniform_real_distribution<double> &u, mt19937 &mt, int *count)
//update clock_time_in_step[level]
{
    cout<<"Start!"<<endl;
    double next_time = (Step*ratio + level + 1)*small_tau;
    int min_loc = find_min(clock_time_in_step[level]); //find_min should return index of clock in time_array
    cout<<"OK, Min_Loc: "<<min_loc<<endl;
    interaction* pt = &time_array[min_loc];
    double current_time = pt->time;
//    cout<<"at level " << level << endl;
    Direction direction_min_loc = pt -> direction;
    Pair coordinate_min_loc = pt -> coordinate; // Find the coordinate of min_loc in 2D energy_array;
    int x = coordinate_min_loc.first / 2; //Compute the x coordinate of cell that the min_loc is corresponding to.
    int y = coordinate_min_loc.second; //Compute the y coordinate of cell......
    List adjacent_clocks_min_loc = pt -> adjacent_clocks;
    Dictionary old_e_value_classifier; // To determine old_e_value I should use to compute the tmp_double 
    Pair left = make_pair(x, y);
    Pair right = (direction_min_loc == vertical) ? make_pair(x, y + 1) : make_pair(x + 1, y);
    old_e_value_classifier[left] = true;
    old_e_value_classifier[right] = true;

    Pair tmp_loc;
    Pair tmp_left, tmp_right;
    int tmp_x = -1;
    int tmp_y = -1;
    Direction dir = vertical;
    int a = 0;
    cout<<"I'm In the Function"<<endl;
    while(current_time < next_time)
    {
        //cout<<"OK"<<a<<endl;
        (*count)++;
        // a++;

        //Step 1: update min interaction and energy
        double total_energy = (direction_min_loc == vertical) ? energy_array[x][y] + energy_array[x][y + 1] : energy_array[x][y] + energy_array[x + 1][y]; //If the direction is horizontal, then we are going to sum the energy of cell at index (x, y) and (x, y + 1). Otherwise, we have to compute the total energy of cell at index (x , y) and (x + 1, y);
        double tmp_double = -log(1 - u(mt))/sqrt(total_energy);
//        cout<<"added time = " << tmp_double << endl;
        double old_e_left = energy_array[x][y];
        double old_e_right = (direction_min_loc == vertical) ? energy_array[x][y + 1] : energy_array[x + 1][y];
        Pair old_e_left_coordinate, old_e_right_coordinate;
        old_e_left_coordinate = make_pair(x, y);
        old_e_right_coordinate = (direction_min_loc == vertical) ? make_pair(x, y + 1) : make_pair(x + 1, y);
        double tmp_rnd_uni = u(mt);
        if(y == 0 && direction_min_loc == vertical)
        {
            total_energy = old_e_right  -log(u(mt))*TL;
        }
        if(y == M && direction_min_loc == vertical)
        {
            total_energy = old_e_left  -log(u(mt))*TR;
        }

        if(direction_min_loc == vertical) {
            if(y != 0) {
                energy_array[x][y] = tmp_rnd_uni*total_energy;
            }

            if(y != M) {
                energy_array[x][y + 1] = (1 - tmp_rnd_uni)*total_energy;//update energy
            }
        }
        else {
            energy_array[x][y] = tmp_rnd_uni*total_energy;
            energy_array[x + 1][y] = (1 - tmp_rnd_uni)*total_energy;//update energy
        }
        move_interaction(clock_time_in_step, pt,small_tau,ratio, Step,  current_time + tmp_double);

        //Step 2: update other interactions
        for(iter i = adjacent_clocks_min_loc.begin() ; i != adjacent_clocks_min_loc.end() ; i++) {
            pt = &time_array[(*i)];
            dir = pt -> direction;
            tmp_loc = pt -> coordinate;
            tmp_x = tmp_loc.first / 2;
            tmp_y = tmp_loc.second;
            tmp_left = make_pair(tmp_x, tmp_y);
            if(dir == vertical) {
                //cout<<"Vertical"<<endl;
                tmp_right = make_pair(tmp_x, tmp_y + 1);
                if(old_e_value_classifier[tmp_left] && direction_min_loc == vertical) {
                    tmp_double = (pt->time - current_time)*sqrt(energy_array[tmp_right.first][tmp_right.second] + old_e_right)/energy_summation(energy_array, tmp_loc.first, tmp_loc.second, dir) + current_time;
                    move_interaction(clock_time_in_step, pt,small_tau,ratio, Step, tmp_double);
                }
                
                if(old_e_value_classifier[tmp_right] && direction_min_loc == vertical) {
                    tmp_double = (pt->time - current_time)*sqrt(energy_array[tmp_left.first][tmp_left.second] + old_e_left)/energy_summation(energy_array, tmp_loc.first, tmp_loc.second, dir) + current_time;
                    move_interaction(clock_time_in_step, pt,small_tau,ratio, Step, tmp_double);
                }
                
                if(old_e_value_classifier[tmp_left] && direction_min_loc == horizontal) {
                    if(tmp_left == left) {
                        tmp_double = (pt->time - current_time)*sqrt(energy_array[tmp_right.first][tmp_right.second] + old_e_left)/energy_summation(energy_array, tmp_loc.first, tmp_loc.second, dir) + current_time;
                        move_interaction(clock_time_in_step, pt,small_tau,ratio, Step, tmp_double);
                    }
                    else if(tmp_left == right){
                        tmp_double = (pt->time - current_time)*sqrt(energy_array[tmp_right.first][tmp_right.second] + old_e_right)/energy_summation(energy_array, tmp_loc.first, tmp_loc.second, dir) + current_time;
                        move_interaction(clock_time_in_step, pt,small_tau,ratio, Step, tmp_double);
                    }
                }
                
                if(old_e_value_classifier[tmp_right] && direction_min_loc == horizontal) {
                    if(tmp_right == left) {
                      tmp_double = (pt->time - current_time)*sqrt(energy_array[tmp_left.first][tmp_left.second] + old_e_left)/energy_summation(energy_array, tmp_loc.first, tmp_loc.second, dir) + current_time;
                      move_interaction(clock_time_in_step, pt,small_tau,ratio, Step, tmp_double);
                    }
                    else if(tmp_left == right){
                      tmp_double = (pt->time - current_time)*sqrt(energy_array[tmp_left.first][tmp_left.second] + old_e_right)/energy_summation(energy_array, tmp_loc.first, tmp_loc.second, dir) + current_time;
                      move_interaction(clock_time_in_step, pt,small_tau,ratio, Step, tmp_double);
                    }
                }
            }
            else {
                if(old_e_value_classifier[tmp_right] && direction_min_loc == vertical) {
                    if(tmp_right == left) {
                        
                    }
                    else if(tmp_left == right) {
                        
                    }
                }
                //cout<<"Horizotal"<<endl;
                //tmp_double = (pt->time - current_time)*sqrt(energy_array[0][0] + old_e_right)/sqrt(energy_array[0][0] + energy_array[0][0]) + current_time;
                  //move_interaction(clock_time_in_step, pt,small_tau,ratio, Step, tmp_double);
            }
            /*if(dir == vertical) {
                tmp_right = make_pair(tmp_x, tmp_y + 1);
                if(tmp_y != 0) {
                    if(tmp_right.first == old_e_left_coordinate.first && tmp_right.second == old_e_left_coordinate.second) {
                        tmp_double = (pt -> time - current_time) * sqrt(energy_array[tmp_left.first][tmp_left.second] + old_e_left)/sqrt(energy_array[tmp_left.first][tmp_left.second]);
                    }
                    else if(tmp_right.first == old_e_right_coordinate.first && tmp_right.second == old_e_right_coordinate.second) {
                        tmp_double = (pt -> time - current_time) * sqrt(energy_array[tmp_left.first][tmp_left.second] + old_e_right)/sqrt(energy_array[tmp_left.first][tmp_left.second]);
                    }
                }

                if(tmp_y != M) {
                    if(tmp_left.first == old_e_right_coordinate.first && tmp_left.second == old_e_right_coordinate.second) {
                        tmp_double = (pt -> time - current_time) * sqrt(energy_array[tmp_right.first][tmp_right.second] + old_e_right)/sqrt(energy_array[tmp_left.first][tmp_left.second]);
                    }
                    else if(tmp_left.first == old_e_left_coordinate.first && tmp_left.second == old_e_left_coordinate.second){
                        tmp_double = (pt -> time - current_time) * sqrt(energy_array[tmp_right.first][tmp_right.second] + old_e_left)/sqrt(energy_array[tmp_right.first][tmp_right.second]);
                    }

                }
            }
            else {
                tmp_right = make_pair(tmp_x + 1, tmp_y);
                if(tmp_right.first == old_e_left_coordinate.first && tmp_right.second == old_e_left_coordinate.second) {
                    tmp_double = (pt -> time - current_time) * sqrt(energy_array[tmp_left.first][tmp_left.second] + old_e_left)/sqrt(energy_array[tmp_left.first][tmp_left.second]);
                }
                else if(tmp_right.first == old_e_right_coordinate.first && tmp_right.second == old_e_right_coordinate.second) {
                    tmp_double = (pt -> time - current_time) * sqrt(energy_array[tmp_left.first][tmp_left.second] + old_e_right)/sqrt(energy_array[tmp_left.first][tmp_left.second]);
                }
                else if(tmp_left.first == old_e_left_coordinate.first && tmp_left.second == old_e_left_coordinate.second) {
                    tmp_double = (pt -> time - current_time) * sqrt(energy_array[tmp_right.first][tmp_right.second] + old_e_left)/sqrt(energy_array[tmp_right.first][tmp_right.second]);
                }
                else if(tmp_left.first == old_e_right_coordinate.first && tmp_left.second == old_e_right_coordinate.second) {
                    tmp_double = (pt -> time - current_time) * sqrt(energy_array[tmp_right.first][tmp_right.second] + old_e_right)/sqrt(energy_array[tmp_left.first][tmp_left.second]);
                }
            }
            move_interaction(clock_time_in_step, pt,small_tau,ratio, Step, tmp_double);*/
        }
        //cout<<"\n";

        //Step 3: update current time
        if(clock_time_in_step[level] != NULL)
        {
            min_loc = find_min(clock_time_in_step[level]);
            pt = &time_array[min_loc];
            current_time = pt->time;
        }
        else
        {
            current_time = next_time + 1;
        }


    }

}

int main(int argc, char** argv)
{
    struct timeval t1, t2;
    ofstream myfile;
    myfile.open("HL_KMP.txt", ios_base::app);
    int N = 5, M = 5;
    if(argc > 1)
    {
        N = strtol(argv[1], NULL, 10);
        M = strtol(argv[2], NULL, 10);
    }
    int N_times_M = N * M; // N is row, M is column
    double big_tau = 0.2;//big time step of tau leaping
    const int ratio = int(N/10);//ratio of big step and small step
    double small_tau = big_tau/double(ratio);//small time step

    int size_of_energy_array = N * (M + 2); //Make the energy array in 2D
    double** energy_array = new double*[N];
    for(int i = 0 ; i < N ; i ++) {
        energy_array[i] = new double[M+2];
    }
    double *E_avg = new double[N_times_M];
    double* last_update = new double[N_times_M];
    for(int i = 0; i <N_times_M; i++)
    {
        E_avg[i] = 0;
        last_update[i] = 0;
    }
    random_device rd;
    mt19937 mt(rd());
    uniform_real_distribution<double> u(0,1);

    //Initialize energy_array
    for(int i = 0 ; i < N ; i++) {
        for(int j = 0 ; j < M + 2 ; j++) {
            if(j == 0) {
                energy_array[i][j] = TL;
            }
            else if(j == M + 1) {
                energy_array[i][j] = TR;
            }
            else {
                energy_array[i][j] = 1;
            }
        }
    }

    int size_of_time_array = N * (M + 1) + (N - 1) * M;
    interaction* time_array = new interaction[size_of_time_array];
    int index = 0;
    int total_number_of_row = 2 * N - 2;
    // int test_index = 0;
    // cout<<total_number_of_row<<" "<<
    for(int i = 0 ; i <= total_number_of_row ; i++) {
        if(i % 2 == 0) {
            for(int j = 0 ; j <= M ; j++) {
                time_array[index].time = -log(1 - u(mt))/energy_summation(energy_array, i, j, vertical); //sqrt(energy_array[index] + energy_array[index+1]);
                time_array[index].location = index;
                time_array[index].coordinate = make_pair(i, j);
                time_array[index].direction = vertical;
                time_array[index].adjacent_clocks = find_adjacent_clocks(i, j, vertical, N, M, index);
                time_array[index].left = NULL;
                time_array[index].right = NULL;
                index++;
            }
        }
        else {
            for(int j = 0 ; j < M ; j++) {
                time_array[index].time = -log(1 - u(mt))/energy_summation(energy_array, i, j, horizontal);
                time_array[index].location = index;
                time_array[index].coordinate = make_pair(i, j);
                time_array[index].direction = horizontal;
                time_array[index].adjacent_clocks = find_adjacent_clocks(i, j, horizontal, N, M, index);
                time_array[index].left = NULL;
                time_array[index].right = NULL;
                index++;
            }
        }
        //cout<<index;
    }

    // for(int i = 0; i < size_of_time_array; i++) {
    //     cout<<"Current Index: "<<time_array[i].location<<endl;
    //     cout<<"Adjacent List: ";
    //     for(iter j = time_array[i].adjacent_clocks.begin() ; j != time_array[i].adjacent_clocks.end() ; j++) {
    //         cout<<(*j)<<" ";
    //     }
    //     cout<<endl;
    //     //<<time_array[i].adjacent_clocks
    // }
    //cout<<N_times_M<<endl;
    //cout<<index;


    /*for(int n = 0; n < size_time_array; n++)
    {
        time_array[n].time = -log(1 - u(mt))/sqrt(energy_array[n] + energy_array[n+1]);
        time_array[n].location = n;
        time_array[n].left = NULL;
        time_array[n].right = NULL;
    }*/
    int count = 0;

    interaction** clock_time_in_step = new interaction*[ratio + 1];//each element in the array is the head of a list
    for(int i = 0; i < ratio + 1; i++)
    {
        clock_time_in_step[i] = NULL;
    }

    gettimeofday(&t1,NULL);



    int Step = 50;
    //double fixed_number = 0.0;
    for(int out_n = 0; out_n < Step; out_n++)
    {
        cout<<out_n<<endl;
        //cin.ignore();
        big_step_distribute(clock_time_in_step,time_array, size_of_time_array ,small_tau,ratio,out_n);
        cout<<"-----------------"<<endl;
        // for(int i = 0 ; i < size_of_time_array ; i++) {
        //     //cout<<-log(1 - u(mt))<<endl;
        //     time_array[i].time = -log(1 - u(mt));
        // }
        //cout<<"Number Of Buckets: "<<ratio<<endl;
        for(int in_n = 0; in_n < ratio; in_n++)
        {
            //cout<<"OK"<<out_n<<endl;
            if(clock_time_in_step[in_n]!= NULL)
            {
                //cout<<"OKOK"<<out_n<<endl;
                /*fixed_number = 0.1 * double(out_n);
                for(interaction *i = clock_time_in_step[in_n] ; i != NULL ; i = i -> right) {
                    //cout<< i -> time <<endl;
                    i -> time = fixed_number;//-log(1 - u(mt));
                    fixed_number += 0.01;
                }*/
                update(clock_time_in_step, in_n, N, M, small_tau, ratio, out_n, time_array, energy_array, u, mt, &count);
            }
        }

        clock_time_in_step[ratio] = NULL;





    }

    gettimeofday(&t2, NULL);
    delete[] energy_array;
    delete[] E_avg;
    delete[] time_array;
    delete[] clock_time_in_step;

    double delta = ((t2.tv_sec  - t1.tv_sec) * 1000000u +
                    t2.tv_usec - t1.tv_usec) / 1.e6;

    cout<<" N = "<<N <<endl;
    myfile<< 1000000*delta/double(count)<<"  ";
    myfile.close();

}
