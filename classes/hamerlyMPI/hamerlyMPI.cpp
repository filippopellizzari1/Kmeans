#include <iostream>
#include <cmath>
using namespace std;

#include "../lloyd/lloyd.h"
#include "../point/point.h"
#include "hamerlyMPI.h"
#include <list>

void print( Dot * d )
{
    cout << "(";
    for ( int i = 0; i < DIMENSION; i++ )
    {
        cout << d->coords[i];

        if ( i < DIMENSION - 1 )
            cout << ", ";
    }
    cout << ") " << endl;
}

MPI_Datatype createType()
{
    // create a new type where we only send x, y and n
    int B[] = {
        DIMENSION, // D int's
        1, // size
        1, // lb
        1, // ub
        1, // centroid_index
        1 // point_index
    };
    MPI_Aint disp[] = {
        offsetof(struct Dot, coords),
        offsetof(struct Dot, size),
        offsetof(struct Dot, lb),
        offsetof(struct Dot, ub),
        offsetof(struct Dot, centroid_index),
        offsetof(struct Dot, point_index)
    };
    MPI_Datatype T[] = {
        MPI_DOUBLE,
        MPI_INT,
        MPI_DOUBLE,
        MPI_DOUBLE,
        MPI_INT,
        MPI_INT
    };

    MPI_Datatype mpi_dt_mystruct;
    MPI_Type_create_struct(6, B, disp, T, &mpi_dt_mystruct);
    MPI_Type_commit(&mpi_dt_mystruct);

    return mpi_dt_mystruct;
}

hamerlyMPI::hamerlyMPI( int _n, int _k,  int * argc, char *** argv )
{
    n = _n;
    k = _k;
    MPI_Init( argc, argv );
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    dot_type = createType();

    points = NULL;
    centroids = NULL;
    items_per_process = NULL;
    items_disp_per_process = NULL;
    cent_per_process = NULL;
    cent_disp_per_process = NULL;
    local_points = NULL;
    local_criticals = NULL;
    average_per_class = NULL;
    relative_average_per_class = NULL;
    points_per_class = NULL;
    relative_points_per_class = NULL;

    // initialization of variables used only by process 0
    if ( rank == 0 )
    {
        points = new Dot[n];

        average_per_class = new double[k * DIMENSION];
        relative_average_per_class = new double[k * DIMENSION];
        points_per_class = new int[k];
        relative_points_per_class = new int[k];

        for ( int i = 0; i < k; i++ )
        {
            points_per_class[i] = 0;

            for ( int j = 0; j < DIMENSION; j++ )
                average_per_class[i * DIMENSION + j] = 0;
        }

        // calculate number of points per each process for the first assignation
        items_per_process = new int[size];
        items_disp_per_process = new int[size];

        for ( int i = 0; i < size; i++ )
        {
            items_per_process[i] = n / size;
            int rem = n % size;

            if ( i < rem )
            {
                items_per_process[i]++;
                rem = 0;
            }

            items_disp_per_process[i] = i * items_per_process[i]  + rem;
        }

        // calculate the number of centorids for every process for the centroid assignation
        cent_per_process = new int[size];
        cent_disp_per_process = new int[size];

        for ( int i = 0; i < size; i++ )
        {
            cent_per_process[i] = k / size;
            int rem = k % size;

            if ( i < rem )
            {
                cent_per_process[i]++;
                rem = 0;
            }

            cent_disp_per_process[i] = i * cent_per_process[i]  + rem;
        }

        n_found_crit_disp = new int[size];
    }

    // initialization of the variables in every process
    int cent_in_process = k / size;
    if ( rank < k % size )
        cent_in_process++;

    local_centroids = new Dot[cent_in_process];
    local_average_per_class = new double[k * DIMENSION];
    local_points_per_class = new int[k];
    local_cdistances = new double[cent_in_process];
    cdistances = new double[k];
    centroids = new Dot[k];
    n_criticals = new int[size];
    first_ass = false;
}

hamerlyMPI::~hamerlyMPI()
{
    delete[] centroids;
    if ( rank == 0 )
    {
        delete[] points;
        delete[] items_per_process;
        delete[] items_disp_per_process;
        delete[] cent_per_process;
        delete[] cent_disp_per_process;
        delete[] average_per_class;
        delete[] relative_average_per_class;
        delete[] points_per_class;
        delete[] relative_points_per_class;

        if ( criticals != NULL )
            delete[] criticals;
 
        if ( n_found_crit_disp != NULL )
            delete[] n_found_crit_disp;
    }

    if ( local_points != NULL )
        delete[] local_points;

    if ( local_centroids != NULL )
        delete[] local_centroids;

    if ( cdistances != NULL )
        delete[] cdistances;

    if ( local_average_per_class != NULL )
        delete[] local_average_per_class;

    if ( local_points_per_class != NULL )
        delete[] local_points_per_class;

    if ( local_cdistances != NULL )
        delete[] local_cdistances;

    if ( n_criticals != NULL )
        delete[] n_criticals;

    MPI_Type_free( &dot_type );
    MPI_Finalize();
}

void hamerlyMPI::init_point( int i, double * coords, int dim )
{
    if ( rank == 0 )
    {
        points[i].size = dim;
        points[i].lb = 0;
        points[i].ub = 0;
        points[i].centroid_index = -1;
        points[i].point_index = i;

        for ( int j = 0; j < dim; j++ )
            points[i].coords[j] = coords[j];

        cout << "adding point" << endl;
    }
}

void hamerlyMPI::init_centroid( int i, double * coords, int dim )
{
    if ( rank == 0 )
    {
        centroids[i].size = dim;
        centroids[i].lb = 0;
        centroids[i].ub = 0;
        centroids[i].centroid_index = -1;
        centroids[i].point_index = i;

        for ( int j = 0; j < dim; j++ )
            centroids[i].coords[j] = coords[j];

        cout << "adding centroid" << endl;
    }

    MPI_Bcast( centroids, k, dot_type, 0, MPI_COMM_WORLD );
}

double hamerlyMPI::distanceA( Dot * p1, Dot * p2 )
{
    if ( p1->size != p2->size )
        throw runtime_error( "Cannot calculate distance between two Dots with differen size A" );

    double distance = 0;

    for ( int i = 0; i < p1->size; i++ )
        distance += pow( p1->coords[i] - p2->coords[i], 2 );

    return sqrt(distance);
}

double hamerlyMPI::distanceB( Dot * p1, Dot * p2 )
{
    if ( p1->size != p2->size )
        throw runtime_error( "Cannot calculate distance between two Dots with differen size B" );

    double distance = 0;

    for ( int i = 0; i < p1->size; i++ )
        distance += pow( p1->coords[i] - p2->coords[i], 2 );

    return sqrt(distance);
}

void hamerlyMPI::first_assignation()
{
    // calculate the number of points assigned to the proccess
    int item_in_process = n / size;
    if ( rank < n % size ) 
        item_in_process++;

    local_points = new Dot[ item_in_process ];
    
    for ( int i = 0; i < k; i++ )
    {
        local_points_per_class[i] = 0;

        for ( int j = 0; j < DIMENSION; j++ )
            local_average_per_class[i * DIMENSION + j] = 0;
    }

    // Distribuite the points accross all processes
    MPI_Scatterv( points, items_per_process, items_disp_per_process, dot_type, local_points, item_in_process, dot_type, 0, MPI_COMM_WORLD );

    // loop through the assigned points, and for every one find its closest centroid
    for ( int i = 0; i < item_in_process; i++ )
    {
        Dot * p = &local_points[i];
        double mindist = distanceA( p, &centroids[0] );
        double secmindist = distanceA( p, &centroids[1] );
        int cent_index = 0;

        if (mindist > secmindist)
        {
            double temp = secmindist;
            secmindist = mindist;
            mindist = temp;
            cent_index = 1;
        }

        for ( int j = 2; j < k; j++ )
        {
            double actdist = distanceA( p, &centroids[j] );

            if ( actdist < mindist )
            {
                secmindist = mindist;
                mindist = actdist;
                cent_index = j;
            }
            else if ( actdist < secmindist )
                secmindist = actdist;
        }

        // account point's coordinates inside the new class
        // update average of all points in that class
        if ( cent_index != p -> centroid_index )
        {
            // the point is changing its class so the number of point per each class must be updated
            if ( p -> centroid_index != -1 )
                local_points_per_class[ p -> centroid_index ] -= 1; 
            
            local_points_per_class[ cent_index ] += 1;

            // since the point is in another class it no longer affect the average of the previous class
            for ( int j = 0; j < DIMENSION; j++ )
            {
                double coord = p->coords[j];

                if ( p->centroid_index != -1 )
                    local_average_per_class[ p -> centroid_index * DIMENSION + j] -= coord;

                local_average_per_class[cent_index * DIMENSION + j] += coord;
            }
        }

        // assign the point to the closest centroid
        p -> centroid_index = cent_index;
        p -> ub = mindist;
        p -> lb = secmindist;
    }

    MPI_Gatherv( local_points, item_in_process, dot_type, points, items_per_process, items_disp_per_process, dot_type, 0, MPI_COMM_WORLD );
    MPI_Reduce( local_points_per_class, relative_points_per_class, k, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );
    MPI_Reduce( local_average_per_class, relative_average_per_class, k * DIMENSION, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );

    if ( rank == 0 )
    {
        for ( int i = 0; i < k; i++ )
        {
            points_per_class[i] += relative_points_per_class[i];

            for ( int j = 0; j < k; j++ )
                average_per_class[i * DIMENSION + j] += relative_average_per_class[i * DIMENSION + j];
        }
    }
    
    first_ass = true;
}   

void hamerlyMPI::update_centroids()
{
    // calculate the number of centroids to assign to each 
    int centroids_in_process = k / size;
    if ( rank < k % size )
        centroids_in_process++;

    // calculate cent_per_process * d to scatter the average_per_class array
    int cent_per_process_dim[size];
    int cent_disp_dim[size];
    int points_per_class_scattered[ centroids_in_process ];
    double average_per_class_scattered[ centroids_in_process * DIMENSION ];

    if ( rank == 0 )
    {
        for ( int i = 0; i < size; i++ )
        {
            cent_per_process_dim[i] = cent_per_process[i] * DIMENSION;
            cent_disp_dim[i] = cent_disp_per_process[i] * DIMENSION;
        }
    }

    MPI_Scatterv( centroids, cent_per_process, cent_disp_per_process, dot_type, local_centroids, centroids_in_process, dot_type, 0, MPI_COMM_WORLD  );
    MPI_Scatterv( points_per_class, cent_per_process, cent_disp_per_process, MPI_INT, points_per_class_scattered, centroids_in_process, MPI_INT, 0, MPI_COMM_WORLD );
    MPI_Scatterv( average_per_class, cent_per_process_dim, cent_disp_dim, MPI_DOUBLE, average_per_class_scattered, centroids_in_process * DIMENSION, MPI_DOUBLE, 0, MPI_COMM_WORLD );

    for ( int i = 0; i < centroids_in_process; i++ )
    {
        Dot * c = &local_centroids[i];
        double sumofsquare = 0;

        if ( points_per_class_scattered[i] != 0 )
            for ( int j = 0; j < DIMENSION; j++ )
            {
                double old_coord = c->coords[j];
                double new_coord = average_per_class_scattered[i * DIMENSION + j] / points_per_class_scattered[i];
                
                sumofsquare += pow( old_coord - new_coord, 2 );
                c->coords[j] = new_coord;
            }

        local_cdistances[i] = sqrt(sumofsquare);
    }

    // Send local centroids, and cdistances back to root process 0
    MPI_Gatherv( local_centroids, centroids_in_process, dot_type, centroids, cent_per_process, cent_disp_per_process, dot_type, 0, MPI_COMM_WORLD );
    MPI_Gatherv( local_cdistances, centroids_in_process, MPI_DOUBLE, cdistances, cent_per_process, cent_disp_per_process, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( cdistances, k, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( centroids, k, dot_type, 0, MPI_COMM_WORLD );

    maxdist = cdistances[0];
    secmaxdist = cdistances[1];
    maxdist_index = 0;

    if (maxdist < secmaxdist)
    {
        double temp = maxdist;
        maxdist = secmaxdist;
        secmaxdist = temp;
        maxdist_index = 1;
    }

    for ( int i = 2; i < k; i++ )
    {
        double act_dist = cdistances[i];

        if ( act_dist > maxdist )
        {
            secmaxdist = maxdist;
            maxdist = act_dist;
            maxdist_index = i;
        }
        else if ( act_dist > secmaxdist )
            secmaxdist = act_dist;
    }

    update_criticals();
}   

void hamerlyMPI::update_criticals()
{
    list<Dot *> found_criticals;
    int n_found = 0;

    int item_in_process = n / size;
    if ( rank < n % size )
        item_in_process++;

    MPI_Scatterv( points, items_per_process, items_disp_per_process, dot_type, local_points, item_in_process, dot_type, 0, MPI_COMM_WORLD );

    for ( int i = 0; i < item_in_process; i++ )
    {
        Dot * p = &local_points[i];
        // update ub
        p -> ub += cdistances[ p -> centroid_index ];

        // update lb
        if ( p -> centroid_index == maxdist_index )
            p -> lb -= secmaxdist;
        else
            p -> lb -= maxdist;

        if ( p -> ub > p -> lb )
        {
            // add critical to the list
            found_criticals.push_back( p );

            // increase the number of 
            n_found += 1;
        }
    }

    MPI_Allgather( &n_found, 1, MPI_INT, n_criticals, 1, MPI_INT, MPI_COMM_WORLD );

    if ( rank == 0 )
    {
        tot_criticals = 0;
        for ( int i = 0; i < size; i++ )
        {
            tot_criticals += n_criticals[i];
        }

        criticals = new Dot[tot_criticals];
    }

    int start_index = 0;
    for ( int i = 0; i < rank; i++ )
        start_index += n_criticals[i];

    if ( local_criticals != NULL )
        delete[] local_criticals;

    local_criticals = new Dot[n_found];
    for ( int i = 0; i < n_found; i++ )
    {
        Dot * p = found_criticals.back();
        local_criticals[i].size = p->size;
        local_criticals[i].lb = p->lb;
        local_criticals[i].ub = p->ub;
        local_criticals[i].centroid_index = p->centroid_index;
        local_criticals[i].point_index = p->point_index;

        for ( int j = 0; j < DIMENSION; j++ )
            local_criticals[i].coords[j] = p -> coords[j];

        found_criticals.pop_back();
    }

    MPI_Gather( &start_index, 1, MPI_INT, n_found_crit_disp, 1, MPI_INT, 0, MPI_COMM_WORLD );
    MPI_Bcast( &tot_criticals, 1, MPI_INT, 0, MPI_COMM_WORLD );
    MPI_Gatherv( local_criticals, n_found, dot_type, criticals, n_criticals, n_found_crit_disp, dot_type, 0, MPI_COMM_WORLD );
}

void hamerlyMPI::refresh_assignation()
{
    if ( rank == 0 )
    {
        for ( int i = 0; i < size; i++ )
        {
            items_per_process[i] = tot_criticals / size;
            int rem = tot_criticals % size;

            if ( i < rem )
            {
                items_per_process[i]++;
                rem = 0;
            }

            items_disp_per_process[i] = i * items_per_process[i] + rem;
        }
    }

    int item_in_process = tot_criticals / size;
    if ( rank < tot_criticals % size )
        item_in_process++;

    if ( local_points != NULL )
        delete[] local_points;

    local_points = new Dot[item_in_process];

    for ( int i = 0; i < k; i++ )
    {
        local_points_per_class[i] = 0;

        for ( int j = 0; j < DIMENSION; j++ )
            local_average_per_class[i * DIMENSION + j] = 0;
    }

    // Distribuite the points accross all processes
    MPI_Scatterv( criticals, items_per_process, items_disp_per_process, dot_type, local_points, item_in_process, dot_type, 0, MPI_COMM_WORLD );

    // loop through the assigned points, and for every one find its closest centroid
    for ( int i = 0; i < item_in_process; i++ )
    {
        Dot * p = &local_points[i];

        cout << "rank" << rank << " " << i << " p: "; print(p);

        double mindist = distanceA( p, &centroids[0] );
        double secmindist = distanceB( p, &centroids[1] );
        int cent_index = 0;

        if (mindist > secmindist)
        {
            double temp = secmindist;
            secmindist = mindist;
            mindist = temp;
            cent_index = 1;
        }

        for ( int j = 2; j < k; j++ )
        {
            double actdist = distanceB( p, &centroids[j] );

            if ( actdist < mindist )
            {
                secmindist = mindist;
                mindist = actdist;
                cent_index = j;
            }
            else if ( actdist < secmindist )
                secmindist = actdist;
        }

        // account point's coordinates inside the new class
        // update average of all points in that class
        if ( cent_index != p -> centroid_index )
        {
            // the point is changing its class so the number of point per each class must be updated
            if ( p -> centroid_index != -1 )
                local_points_per_class[ p -> centroid_index ] -= 1; 
            
            local_points_per_class[ cent_index ] += 1;

            // since the point is in another class it no longer affect the average of the previous class
            for ( int j = 0; j < DIMENSION; j++ )
            {
                double coord = p->coords[j];

                if ( p->centroid_index != -1 )
                    local_average_per_class[ p -> centroid_index * DIMENSION + j] -= coord;

                local_average_per_class[cent_index * DIMENSION + j] += coord;
            }
        }

        // assign the point to the closest centroid
        p -> centroid_index = cent_index;
        p -> ub = mindist;
        p -> lb = secmindist;
    }

    MPI_Gatherv( local_points, item_in_process, dot_type, criticals, items_per_process, items_disp_per_process, dot_type, 0, MPI_COMM_WORLD );
    MPI_Reduce( local_points_per_class, relative_points_per_class, k, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );
    MPI_Reduce( local_average_per_class, relative_average_per_class, k * DIMENSION, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );

    if ( rank == 0 )
    {
        for ( int i = 0; i < k; i++ )
        {
            points_per_class[i] += relative_points_per_class[i];

            for ( int j = 0; j < k; j++ )
                average_per_class[i * DIMENSION + j] += relative_average_per_class[i * DIMENSION + j];
        }

        for ( int i = 0; i < tot_criticals; i++ )
        {
            Dot * p = &criticals[i];
            points[p->point_index].lb = p->lb;
            points[p->point_index].ub = p->ub;
            points[p->point_index].centroid_index = p->centroid_index;
        }
    }
}

bool hamerlyMPI::double_eq( double f1, double f2 )
{
    bool upper_cond = f1 < f2 + EPS;
    bool lower_cond = f1 > f2 - EPS;

    return upper_cond && lower_cond;
}

void hamerlyMPI::checkequal( lloyd * l )
{
    if ( rank == 0 )
    {
        bool equal = true;

        for ( int i = 0; i < n && equal; i++ )
        {
            Point * lloyd_coord = l->get_point(i);
            for ( int j = 0; j < DIMENSION && equal; j++ )
            {
                if (lloyd_coord->centroid_index != points[i].centroid_index )
                {
                    cout << "not same cent" << endl;
                    equal = false;
                }
            }
        }

        for ( int i = 0; i < k && equal; i++ )
        {
            Point * lloyd_coord = l->get_centroid(i);
            for ( int j = 0; j < DIMENSION && equal; j++ )
            {
                if ( !hamerlyMPI::double_eq( lloyd_coord -> get_coord(j), centroids[i].coords[j] ) )
                {
                    cout << "centroid not same position" << endl;
                    equal = false;
                }
            }
        }

        if (equal)
            cout << "equal" << endl;
        else
            cout << "not equal" << endl;
    }
}

void hamerlyMPI::assign_points()
{
    if ( first_ass )
        refresh_assignation();
    else
        first_assignation();
}

    // for ( int i = 0; i < item_in_process; i++ )
    // {
    //     cout << "rank: " << rank << " (";
    //     for ( int j = 0; j < DIMENSION; j++ )
    //     {
    //         cout << local_points[i].coords[j];
    //         if ( j < DIMENSION - 1 )
    //             cout << ", ";
    //     }
    //     cout << ") " << endl;
    // }

    // for ( int i = 0; i < k; i++ )
    // {
    //     cout << "CENT rank: " << rank << " (";
    //     for ( int j = 0; j < DIMENSION; j++ )
    //     {
    //         cout << centroids[i].coords[j];
    //         if ( j < DIMENSION - 1 )
    //             cout << ", ";
    //     }
    //     cout << ") " << endl;
    // }  
