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
    int nThreads;
    int C; // Total number of criticals
    double maxdist; // last distance taken from the centroid which moved the most
    double secmaxdist; // last distance taken from the second centroid which moved the most
    int maxdist_index; // index the centroid which moved the most

    Point ** criticals;

    double * cdistances; // this array contains the distances moved by the centroids during their last update

    list<Point *> * critical_in_threads; // this array contains a list for each thread that in turn contains all the criticals found from that thread

    int * nC; // this array contains per each thread the number of criticals that it found 

    bool first_ass;
    
    void first_assignation();
    void refresh_assignation();
    void update_criticals();
    public:
        hamerly( int _d, int _N, int _K, int _nThreads );
        ~hamerly();

        void update_centroids();
        void assign_points();

        void print_criticals();
        int get_criticalsNumber();
};

#endif
