#include <iostream>
#include "classes/kmeans/kmeans.h"

using namespace std;

#define N 100 // Number of points
#define K 3 // Number of centroids
#define D 2 // Data dimensionality
#define T 3 // number of threads

int main()
{
    kMeans model(D, N, K);
    model.generate_points( 0, 50 );
    model.generate_centroids( 0, 50 );

    cout << "Points" << endl;
    model.print_points();
    cout << "Centroids" << endl;
    model.print_centroids();

    model.assign_points( 1 );

    cout << "Points after assignation" << endl;
    model.print_points();
    cout << "Centroids" << endl;
    model.print_centroids();

    model.update_centroid( 1 );

    cout << "Points after centriod moved" << endl;
    model.print_points();
    cout << "Centroids" << endl;
    model.print_centroids();

    cout << "Criticals : " << model.get_c() << endl;
}