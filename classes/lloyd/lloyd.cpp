#include "lloyd.h"
#include "../point/point.h"

#include <string>
#include <cmath>
#include <fstream>

lloyd::lloyd( int _d, int _N, int _K )
{
    d = _d;
    N = _N;
    K = _K;

    if ( K < 2 )
        throw invalid_argument( "Cannot instantiate less than 2 centroids" );

    points = new Point[N];
    centroids = new Point[K];

    average_per_class = new double[K*d];
    for ( int i = 0; i < K*d; i++ ) 
        average_per_class[i] = 0;
    
    points_per_class = new int[K];
    for ( int i = 0; i < K; i++ )
        points_per_class[i] = 0;

    srand( time( NULL ) );
}

lloyd::~lloyd()
{
    delete[] points;
    delete[] centroids; 

    delete[] average_per_class;

    delete[] points_per_class;
}

void lloyd::generate_points( int min, int max )
{
    double * coords = new double[d]; // Coordinates for each points
    for ( int i = 0; i < N; i++ )
    {
        // Generate coord
        for ( int j = 0; j < d; j++ )
            coords[j] = random_double( min, max );

        // Insert the point in the list
        Point temp( d, coords );
        points[i] = temp;
    }

    delete[] coords;
}

void lloyd::generate_centroids( int min, int max )
{
    double * coords = new double[d]; // Coordinates for each points
    for ( int i = 0; i < K; i++ )
    {
        // Generate coord
        for ( int j = 0; j < d; j++ )
            coords[j] = random_double( min, max );

        // Insert the point in the list
        Point temp( d, coords );
        centroids[i] = temp;
    }

    delete[] coords;
}

void lloyd::insert_point( int i, const Point * p )
{
    points[i] = *p;
}

void lloyd::insert_centroid( int i, const Point * c )
{
    centroids[i] = *c;
}

Point * lloyd::get_point( int i )
{
    return &points[i];
}

Point * lloyd::get_centroid( int i )
{
    return &centroids[i];
}

void lloyd::assign_points()
{
    for ( int i = 0; i < N; i++ )
    {
        Point * p = &points[i];

        // initiate mindist
        double mindist = p -> distance( centroids[0] );
        int mindist_index = 0;

        for ( int j = 1; j < K; j++ )
        {
            // calculate actual distance
            double actdist = p -> distance( centroids[j] );

            // if distance is lower than mindist
            if ( actdist < mindist )
            {
                // update the mindist
                mindist = actdist;
                mindist_index = j;
            }
        }

        // the point must be account for the average of the class that has been assigned to it
        for ( int j = 0; j < d; j++ )
        {
            double coord = p -> get_coord( j );
            average_per_class[mindist_index * d + j] += coord; 
                    
            //update average of all points in that class
            if ( p -> centroid_index != -1 ) // if the centroid_index is different from -1 it means that the point has been previously assigned to a class, thereforse its contribuition must be removed from that class and added to the new one
                average_per_class[p -> centroid_index * d + j] -= coord;
        }
        
        points_per_class[mindist_index] += 1;
        if ( p -> centroid_index != -1 )
        {
            points_per_class[ p -> centroid_index ] -= 1;
        }

        // assign the point to the right centroid
        p -> set_centroid( &centroids[mindist_index], mindist_index );
    }

    // for ( int i = 0; i < K; i++ )
    //     cout << "Lloyd class" << i << ": " << points_per_class[i] << endl;
}

void lloyd::update_centroids()
{
    for ( int i = 0; i < K; i++ )
    {
        Point * c = &centroids[i];

        int point_in_class = points_per_class[i];
        
        if (point_in_class > 0)
            for ( int j = 0; j < d; j++ )
            {
                double new_coord = average_per_class[i * d + j] / point_in_class;
                c -> set_coord( j, new_coord );
            }
    }
}

void lloyd::print_points()
{
    for ( int i = 0; i < N; i++ )
        cout << points[i] << endl;
}

void lloyd::print_centroids()
{
    for ( int i = 0; i < K; i++ )
        cout << centroids[i] << endl;
}

double lloyd::random_double( int min, int max )
{
    max -= 1;
    return min + (rand() % (max - min)) + (double)(rand()) / (double)(RAND_MAX);
}

void lloyd::export_data( string filename )
{
    // open the file
    fstream file;
    file.open( filename, ios::out );

    if ( !file.is_open() )
        throw runtime_error( "Couldn't open the file " + filename );

    file << "number_points: " << N ;
    file << "; number_centroids: " << K << endl;

    // Write all the centroids
    file << "Centroids" << endl;
    for ( int i = 0; i < K; i++ )
        file << centroids[i] << endl;

    // white all the points
    file << "Points" << endl;
    for ( int i = 0; i < N; i++ )
        file << points[i] << endl;

    file.close();
}

bool operator == ( const lloyd &l1, const lloyd &l2 )
{
    if ( l1.N != l2.N || l1.K != l2.K )
        return false;

    bool equal = true;

    for ( int i = 0; i < l1.N && equal; i++ )
        if ( l1.points[i] != l2.points[i] )
            equal = false;

    for ( int i = 0; i < l1.K && equal; i++ )
        if ( !(l1.centroids[i] == l2.centroids[i]) )
            equal = false;

    return equal;
}