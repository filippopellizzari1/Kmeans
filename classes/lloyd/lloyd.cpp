#include "lloyd.h"
#include "../point/point.h"

#include <string>
#include <cmath>
#include <fstream>

lloyd::lloyd()
{
    d = 0;
    N = 0;
    K = 0;
    
    points = NULL;
    centroids = NULL;

    average_per_class = NULL;
    points_per_class = NULL;
}

lloyd::lloyd( int _d, int _N, int _K )
{
    d = _d;
    N = _N;
    K = _K;

    points = new Point[N];
    centroids = new Point[K];

    average_per_class = new float*[K];
    for ( int i = 0; i < K; i++ )
    {
        average_per_class[i] = new float[d];
        for ( int j = 0; j < d; j++ )
            average_per_class[i][j] = 0;
    }
    
    points_per_class = new int[K];
    for ( int i = 0; i < K; i++ )
        points_per_class[i] = 0;

    srand( time( NULL ) );
}

lloyd::~lloyd()
{
    delete[] points;
    delete[] centroids; 

    for ( int i = 0; i < K; i++ )
        delete[] average_per_class[i];
    delete[] average_per_class;

    delete[] points_per_class;
}

void lloyd::generate_points( int min, int max )
{
    float coords[ d ]; // Coordinates for each points
    for ( int i = 0; i < N; i++ )
    {
        // Generate coord
        for ( int j = 0; j < d; j++ )
            coords[j] = random_float( min, max );

        // Insert the point in the list
        Point temp( d, coords );
        points[i] = temp;
    }
}

void lloyd::generate_centroids( int min, int max )
{
    float coords[ d ]; // Coordinates for each points
    for ( int i = 0; i < K; i++ )
    {
        // Generate coord
        for ( int j = 0; j < d; j++ )
            coords[j] = random_float( min, max );

        // Insert the point in the list
        Point temp( d, coords );
        centroids[i] = temp;
    }
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
        float mindist = p -> distance( centroids[0] );
        int mindist_index = 0;

        for ( int j = 1; j < K; j++ )
        {
            // calculate actual distance
            float actdist = p -> distance( centroids[j] );

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
            float coord = p -> get_coord(j);
            average_per_class[mindist_index][j] += coord; 
            points_per_class[mindist_index] += 1;
                    
            // update average of all points in that class
            if ( p -> centroid_index != -1 ) // if the centroid_index is minus different from minus one it means that the point has been previously assigned to a class, thereforse its contribuition must be removed from that class and added to the new one
            {
                average_per_class[ p -> centroid_index ][j] -= coord;
                points_per_class[ p -> centroid_index ] -= 1;
            }
        }

        // assign the point to the right centroid
        p -> set_centroid( &centroids[mindist_index], mindist_index );
    }
}

void lloyd::update_centroids()
{
    for ( int i = 0; i < K; i++ )
    {
        Point * c = &centroids[i];
        for ( int j = 0; j < d; j++ )
        {
            float new_coord = average_per_class[i][j] / points_per_class[i];
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

float lloyd::random_float( int min, int max )
{
    max -= 1;
    return min + (rand() % (max - min)) + (float)(rand()) / (float)(RAND_MAX);
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

void lloyd::load_data( string filename )
{
    fstream file;
    file.open( filename, ios::in );

    if ( !file.is_open() )
        throw runtime_error( "Couldn't open the file " + filename );

    // read the first line to find the number of centroids and points
    string line;
    getline( file, line );

    // check if the first line syntax is correct
    int ncent_index = line.find( "number_centroids: " );
    int npoint_index = line.find( "number_points: " );
    if ( ncent_index < 0 || npoint_index < 0 )
        throw runtime_error( "number_centroids or number_points not found in " + filename );
    
    // extrapolate number of centroids ( +18 = "number_centroids: ".length )
    K = stoi( line.substr( ncent_index + 18, -1 ) );
    // extrapolate number of points ( +15 = "number_points: ".length )
    N = stoi( line.substr( npoint_index + 15, line.find( ";" ) ) );
    
    // remove old points and centroids
    delete[] points;
    delete[] centroids;

    // new array for points and centroids
    points = new Point[N];
    centroids = new Point[K];

    // iterate through all the remaining lines
    bool cent_tag_found = false;
    bool points_tag_found = false;

    // supporting variables
    int delimiter; 
    string dim_s, coords_s, cent_index_s, lb_s, ub_s, dist_s;
    int dim, cent_index, lb, ub;
    float dist;
    float * coords;
    int next_pointindex = 0;
    int next_centindex = 0;
    while ( getline( file, line ) )
    {
        if ( cent_tag_found && line != "Points" ) // load a generic point
        { 
            // find the dimensionality
            delimiter = line.find( "D: " ) + 3;
            dim_s = line.substr( delimiter, line.find( "; LB" ) - delimiter );
            dim = stoi( dim_s );

            // find the coordinates
            delimiter = line.find( "Coords: " ) + 8;
            coords_s = line.substr( delimiter, line.find( "; C" ) - delimiter );

            coords = new float[dim];
            string coord_s; // i-th coords_s
            int last_pos = 0; // index of the last found ',' in coords_s
            for ( int i = 0; i < dim; i++ )
            {
                int left_bound = last_pos;
                int right_bound = coords_s.find( ", ", last_pos + 1 );
                coord_s = coords_s.substr( left_bound, right_bound );

                coords[i] = stof( coord_s );

                // update index of last found ',' ( +2 = ','.length )
                last_pos = right_bound + 2;
            }

            // find centroid_index
            delimiter = line.find( "C: " ) + 3;
            cent_index_s = line.substr( delimiter, line.find( "; D" ) - delimiter );
            cent_index = stoi( cent_index_s );

            // find lower bound
            delimiter = line.find( "LB: " ) + 4;
            lb_s = line.substr( delimiter, line.find( "; UB" ) - delimiter );
            lb = stoi( lb_s );

            // find upper bound
            delimiter = line.find( "UB: " ) + 4;
            ub_s = line.substr( delimiter, line.find( "; DIST" )  - delimiter );
            ub = stoi( ub_s );

            // find the Point.dist parameter
            delimiter = line.find( "DIST: " ) + 6;
            dist_s = line.substr( delimiter, -1 );
            dist = stof( dist_s );

            cout << dim << "," << coords_s << "," << cent_index << "," << lb << "," << ub << "," << dist << endl;
 
            Point temp( dim, coords, dist, ub, lb, cent_index );
            if ( cent_tag_found && !points_tag_found ) // the point is a centroid
                centroids[next_centindex++] = temp;
            else if ( cent_tag_found && points_tag_found ) // the point is an effective point
                points[next_pointindex++] = temp;

            // Bound check in case the file is not formatted in the right way
            if ( next_centindex >= K )
                throw runtime_error( "Declared number of centroids is lower than the effective number" );

            if ( next_pointindex >= N )
                throw runtime_error( "Declared number of centroids is lower than the effective number" );

            delete[] coords;
        }
        else if ( cent_tag_found && line == "Points" )
            points_tag_found = true;

        if ( !cent_tag_found && line != "Centroids" )
            throw runtime_error( "Couldn't find Centroids tag in " + filename );
        else
            cent_tag_found = true;
    }

    cout << "HERE" << endl;

    cout << points_tag_found << "," << cent_tag_found << endl;
    if ( !points_tag_found )
        throw runtime_error( "Couldn't find Points tag in " + filename );

    file.close();
}