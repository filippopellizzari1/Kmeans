#ifndef __DATASET_H__
#define __DATASET_H__

#include <iostream>
using namespace std;

class dataset
{
    private:
        int d; // data dimensionality
        int N; // nnumber of points
        int K; // number of centroids

        double ** points; // matrix of points
        double ** centroids; // matrix of centorids

        static bool double_eq( double f1, double f2 );
    public:
        dataset( int _d, int _N, int _K );
        ~dataset();

        void rnd_dataset( double min, double max );
        void export_dataset( string points_file, string centroids_file );
        void load_dataset( string points_file, string centroids_file );
        double * get_point( int i );
        double * get_centroid( int i );
        double random_double( double min, double max );

        friend bool operator == ( const dataset &d1, const dataset &d2 );
};

#endif
