#ifndef __LLOYD_H__
#define __LLOYD_H__

#include <string>
using namespace std;

class Point;

class lloyd
{

    int d;  // Data dimensionality
    int N;  // Total number of points
    int K;  // Total number of centroids

    Point * points;
    Point * centroids;

    float ** average_per_class; // this array contains the averaged coordinate of all points in a class
    int * points_per_class;

    float random_float( int min, int max );
    public:
        lloyd();
        lloyd( int _d, int _N, int _K );
        ~lloyd();

        void generate_points( int min, int max ); // generate random points
        void generate_centroids( int min, int max ); // generate random centroids
        void insert_point( int i, const Point * p );
        void insert_centroid( int i, const Point * c );
        Point * get_point( int i );
        Point * get_centroid( int i );

        void update_centroids();
        void assign_points();
        
        void export_data( string filename );
        void load_data( string filename );

        void print_points();
        void print_centroids();
};

#endif
