#ifndef __HAMERLY_H__
#define __HAMERLY_H__

#include <iostream>
#include <string>
#include <vector>
#include <list>
#include "../lloyd/lloyd.h"
using namespace std;

class Point;

class hamerly : public lloyd
{
    int C; // Total number of criticals
    float maxdist; // last distance taken from the centroid which moved the most
    float secmaxdist; // last distance taken from the second centroid which moved the most
    int maxdist_index; // index the centroid which moved the most

    Point ** criticals;

    float * cdistances; // this array contains the distances moved by the centroids during their last update

    list<Point *> * critical_in_threads; // this array contains a list for each thread that in turn contains all the criticals found from that thread

    int * nC; // this array contains per each thread the number of criticals that it found 

    void first_assignation( int nThreads );
    void refresh_assignation( int nThreads );
    void update_criticals( int nThreads );
    public:
        hamerly( int _d, int _N, int _K );
        ~hamerly();

        void update_centroids( int nThreads );
        void assign_points( int nThreads );

        void print_criticals();
        int get_criticalsNumber();
};

#endif
