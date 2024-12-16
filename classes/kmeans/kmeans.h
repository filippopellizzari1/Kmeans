#ifndef __KMEANS_H__
#define __KMEANS_H__

#include <iostream>
#include <string>
#include <vector>
#include <list>
using namespace std;

class Point;

class kMeans
{
    int d;  // Data dimensionality
    int N;  // Total number of points
    int K;  // Total number of centroids
    int C; // Total number of criticals
    float maxdist; // last distance taken from the centroid which moved the most
    float secmaxdist; // last distance taken from the second centroid which moved the most
    int maxdist_index; // index the centroid which moved the most

    Point * points;
    Point * centroids;
    Point ** criticals;

    float ** average_per_class; // this array contains the averaged coordinate of all points in a class
    int * points_per_class;
    float * cdistances; // this array contains the distances moved by the centroids during their last update

    list<Point *> * critical_in_threads; // this array contains a list for each thread that in turn contains all the criticals found from that thread

    int * nC; // this array contains per each thread the number of criticals that it found 

    void first_assignation( int nThreads );
    void refresh_assignation( int nThreads );
    void update_criticals( int nThreads );
    public:
        kMeans( int _d, int _N, int _K );
        ~kMeans();

        float random_float( int min, int max );
        void generate_points( int min, int max ); // generate random points
        void generate_centroids( int min, int max ); // generate random centroids
        void insert_point( int i, const Point * p );
        void insert_centroid( int i, const Point * c );
        Point * get_point( int i );
        Point * get_centroid( int i );

        void update_centroid( int nThreads );
        void assign_points( int nThreads );

        void export_data( string filename );
        void load_data( string filename );
        void print_points();
        void print_centroids();
        void print_criticals();
        int get_c();
};

#endif
