#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;

#include <chrono>
#include "../lloyd/lloyd.h"
#include "../point/point.h"
#include "MPIHamerly.h"
#include <list>

void MPIHamerly::splitItems( int nItems )
{
    if ( rank == 0 )
    {
        for ( int i = 0; i < size; i++ )
        {
            int items = nItems / size;
            int rem = nItems % size;

            if ( i < rem )
            {
                items++;
                rem = 0;
            }

            items_per_proc[i] = items;
            disp_proc[i] = i * items + rem;
        }
    }

    proc_items = nItems / size;
    if ( rank < nItems % size )
        proc_items++;

    if ( local_items != NULL )
        delete[] local_items;

    local_items = new Dot[proc_items];
}

void printDot( Dot * d )
{
    cout << "(";
    for ( int i = 0; i < D; i++ )
    {
        cout << d->coords[i];

        if ( i < D - 1 )
            cout << ", ";
    }
    cout << ") cent: " << d -> centroid_index << " lb: " << d->lb << " ub: " << d->ub << endl;
}

void MPIHamerly::print()
{
    if ( rank == 0 )
    {
        for ( int i = 0; i < n; i++ )
            printDot(&points[i]);
    }
}

MPI_Datatype createType()
{
    // create a new type where we only send x, y and n
    int B[] = {
        D, // D int's
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
    MPI_Datatype types[] = {
        MPI_DOUBLE,
        MPI_INT,
        MPI_DOUBLE,
        MPI_DOUBLE,
        MPI_INT,
        MPI_INT
    };

    MPI_Datatype mpi_dt_mystruct;
    MPI_Type_create_struct(6, B, disp, types, &mpi_dt_mystruct);
    MPI_Type_commit(&mpi_dt_mystruct);

    return mpi_dt_mystruct;
}

MPIHamerly::MPIHamerly(int * argc, char *** argv)
{
    n = N;
    k = K;
    MPI_Init( argc, argv );
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    dot_type = createType();

    centroids = new Dot[k];
    first_ass = false;
    local_items = NULL;
    local_avg_class = new double[k * D];
    local_points_class = new int[k];
    cdistances = new double[k];
    crits_proc = new int[size];

    int cent_per_proc = k / size;
    if ( rank < k % size )
        cent_per_proc++;

    local_cdistances = new double[cent_per_proc];
    if ( rank == 0 )
    {
        points = new Dot[n];

        items_per_proc = new int[size];
        disp_proc = new int[size];

        rel_avg_class = new double[k * D];
        rel_points_class = new int[k];
        avg_class = new double[k * D];
        points_class = new int[k];
        for ( int i = 0; i < k; i++ )
        {
            points_class[i] = 0;

            for ( int j = 0; j < D; j++ )
                avg_class[ i * D + j ] = 0;
        }

        criticals = NULL;
    }
}

MPIHamerly::~MPIHamerly()
{
    delete[] centroids;
    delete[] local_avg_class;
    delete[] local_points_class;
    delete[] local_cdistances;
    delete[] cdistances;
    delete[] crits_proc;

    if ( local_items != NULL )
        delete[] local_items;

    if ( rank == 0 )
    {
        delete[] points;
        delete[] items_per_proc;
        delete[] disp_proc;
        delete[] avg_class;
        delete[] points_class;
        delete[] rel_avg_class;
        delete[] rel_points_class;

        if ( criticals != NULL )
            delete[] criticals;
    }

    MPI_Type_free( &dot_type );
    MPI_Finalize();
}

void MPIHamerly::init_point( int i, double * coords )
{
    if ( rank == 0 )
    {
        points[i].size = D;
        points[i].lb = 0;
        points[i].ub = 0;
        points[i].centroid_index = -1;
        points[i].point_index = i;

        for ( int j = 0; j < D; j++ )
            points[i].coords[j] = coords[j];
    }
}

void MPIHamerly::init_centroid( int i, double * coords )
{
    if ( rank == 0 )
    {
        centroids[i].size = D;
        centroids[i].lb = 0;
        centroids[i].ub = 0;
        centroids[i].centroid_index = -1;
        centroids[i].point_index = i;

        for ( int j = 0; j < D; j++ )
            centroids[i].coords[j] = coords[j];
    }

    MPI_Bcast( centroids, k, dot_type, 0, MPI_COMM_WORLD );
}

double MPIHamerly::distance( Dot * p1, Dot * p2 )
{
    if ( p1->size != p2->size )
        throw runtime_error( "Cannot calculate distance between two Dots with differen size A" );

    double distance = 0;

    for ( int i = 0; i < p1->size; i++ )
        distance += pow( p1->coords[i] - p2->coords[i], 2 );

    return sqrt(distance);
}

int MPIHamerly::closest_centroid( Dot * p, double * distmin, double * distsecmin )
{
    double mindist = distance( p, &centroids[0] );
    double secmindist = distance( p, &centroids[1] );
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
        double actdist = distance( p, &centroids[j] );

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
            local_points_class[ p -> centroid_index ] -= 1; 
        
        local_points_class[ cent_index ] += 1;

        // since the point is in another class it no longer affect the average of the previous class
        for ( int j = 0; j < D; j++ )
        {
            double coord = p->coords[j];

            if ( p->centroid_index != -1 )
                local_avg_class[ p -> centroid_index * D + j] -= coord;

            local_avg_class[cent_index * D + j] += coord;
        }
    }

    *distmin = mindist;
    *distsecmin = secmindist;
    return cent_index;
}

void MPIHamerly::refresh_assignation()
{
    if ( n_crit > 0 )
    {
        splitItems( n_crit );

        MPI_Scatterv( criticals, items_per_proc, disp_proc, dot_type, local_items, proc_items, dot_type, 0, MPI_COMM_WORLD );

        for ( int i = 0; i < proc_items; i++ )
        {
            Dot * p = &local_items[i];
            double mindist, secmindist;
            int cent_index = closest_centroid( p, &mindist, &secmindist );

            // assign to the point the closest centroid
            p -> centroid_index = cent_index;
            p -> ub = mindist;
            p -> lb = secmindist;
        }

        MPI_Gatherv( local_items, proc_items, dot_type, criticals, items_per_proc, disp_proc, dot_type, 0, MPI_COMM_WORLD );
        MPI_Reduce( local_points_class, rel_points_class, k, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );
        MPI_Reduce( local_avg_class, rel_avg_class, k * D, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );

        // process zero updates avg_class and point_class by summing the respecctive relative variable which contains the total changes done by all processes in this iteration
        if ( rank == 0 )
        {
            for ( int i = 0; i < k; i++ )
            {
                points_class[i] += rel_points_class[i];

                for ( int j = 0; j < D; j++ )
                    avg_class[i * D + j] += rel_avg_class[i * D + j];
            }

            // process zero has also to update the points by changing the values of the critical points with the one in the criticals array (which has been modified by the processes)
            for ( int i = 0; i < n_crit; i++ )
            {
                Dot * p = &criticals[i];
                points[p->point_index].lb = p->lb;
                points[p->point_index].ub = p->ub;
                points[p->point_index].centroid_index = p->centroid_index;
            }
        }
    }
}

void MPIHamerly::first_assignation()
{
    splitItems( n ); // calculate how many points to assign per each process

    MPI_Scatterv( points, items_per_proc, disp_proc, dot_type, local_items, proc_items, dot_type, 0, MPI_COMM_WORLD );

    for ( int i = 0; i < proc_items; i++ )
    {
        Dot * p = &local_items[i];
        double mindist, secmindist;
        int cent_index = closest_centroid( p, &mindist, &secmindist );

        p -> centroid_index = cent_index;
        p -> ub = mindist;
        p -> lb = secmindist;
    }

    MPI_Gatherv( local_items, proc_items, dot_type, points, items_per_proc, disp_proc, dot_type, 0, MPI_COMM_WORLD );
    MPI_Reduce( local_points_class, rel_points_class, k, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );
    MPI_Reduce( local_avg_class, rel_avg_class, k * D, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
    
    // process zero updates avg_class and point_class by summing the respecctive relative variable which contains the total changes done by all processes in this iteration
    if ( rank == 0 )
    {
        for ( int i = 0; i < k; i++ )
        {
            points_class[i] += rel_points_class[i];

            for ( int j = 0; j < D; j++ )
                avg_class[i * D + j] += rel_avg_class[i * D + j];
        }
    }

    first_ass = true;
}

void MPIHamerly::preAssignation_init()
{
    for ( int i = 0; i < k; i++ )
    {
        local_points_class[i] = 0;

        for ( int j = 0; j < D; j++ )
            local_avg_class[ i * D + j ] = 0;
    }
}

void MPIHamerly::assign_points()
{
    preAssignation_init();

    if ( first_ass )
        refresh_assignation();
    else
        first_assignation();
}

int findTwoMax( double * arr, int dim, double * max, double * secmax )
{
    if ( dim < 2 )
        throw runtime_error( "cannot find max and second max in an array with size:" + dim );

    int maxdist = arr[0];
    int secmaxdist = arr[1];
    int maxdist_index = 0;

    if (maxdist < secmaxdist)
    {
        double temp = maxdist;
        maxdist = secmaxdist;
        secmaxdist = temp;
        maxdist_index = 1;
    }

    for ( int i = 2; i < dim; i++ )
    {
        double act_dist = arr[i];

        if ( act_dist > maxdist )
        {
            secmaxdist = maxdist;
            maxdist = act_dist;
            maxdist_index = i;
        }
        else if ( act_dist > secmaxdist )
            secmaxdist = act_dist;
    }

    *max = maxdist;
    *secmax = secmaxdist;
    
    return maxdist_index;
}

void MPIHamerly::update_centroids()
{
    splitItems( k );

    int items_per_procD[size]; // calculate items_per_proc * D and disp_proc * D because the array avg_class has a lenght of k * D
    int disp_procD[size];

    if ( rank == 0 )
    {
        for ( int i = 0; i < size; i++ )
        {
            items_per_procD[i] = items_per_proc[i] * D;
            disp_procD[i] = disp_proc[i] * D;
        }
    }

    MPI_Scatterv( centroids, items_per_proc, disp_proc, dot_type, local_items, proc_items, dot_type, 0, MPI_COMM_WORLD );
    MPI_Scatterv( points_class, items_per_proc, disp_proc, MPI_INT, local_points_class, proc_items, MPI_INT, 0, MPI_COMM_WORLD );
    MPI_Scatterv( avg_class, items_per_procD, disp_procD, MPI_DOUBLE, local_avg_class, proc_items * D, MPI_DOUBLE, 0, MPI_COMM_WORLD );

    // update the centroids
    for ( int i = 0; i < proc_items; i++ )
    {
        Dot * c = &local_items[i];
        double sumofsquare = 0;

        if ( local_points_class[i] != 0 )
        {
            for ( int j = 0; j < D; j++ )
            {
                double old_coord = c->coords[j];
                double new_coord = local_avg_class[i * D + j] / local_points_class[i];
                
                sumofsquare += pow( old_coord - new_coord, 2 );
                c->coords[j] = new_coord;
            }
        }

        local_cdistances[i] = sqrt(sumofsquare);
    }

    MPI_Gatherv( local_items, proc_items, dot_type, centroids, items_per_proc, disp_proc, dot_type, 0, MPI_COMM_WORLD );
    MPI_Gatherv( local_cdistances, proc_items, MPI_DOUBLE, cdistances, items_per_proc, disp_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( cdistances, k, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( centroids, k, dot_type, 0, MPI_COMM_WORLD );

    maxdist_index = findTwoMax( cdistances, k, &maxdist, &secmaxdist );


    update_criticals();
}

void MPIHamerly::update_criticals()
{
    splitItems( n );


    MPI_Scatterv( points, items_per_proc, disp_proc, dot_type, local_items, proc_items, dot_type, 0, MPI_COMM_WORLD );

    list<Dot *> found_criticals; // criticals found by each process
    int n_found = 0;    // number of criticals found by each process


    for ( int i = 0; i < proc_items; i++ )
    {
        Dot * p = &local_items[i];
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

    MPI_Allgather( &n_found, 1, MPI_INT, crits_proc, 1, MPI_INT, MPI_COMM_WORLD );
    MPI_Allreduce( &n_found, &n_crit, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );

    // process zero create the array that will contains all the criticals
    if ( rank == 0 )
    {
        if ( criticals != NULL )
            delete[] criticals;

        if ( n_crit > 0 ) criticals = new Dot[n_crit];
        else criticals = NULL;
    }

    // every proccess sums the number of crits found by the processes with a lower rank. This way it calculates the index from where it can start to insert its crit points in the array "criticals"
    int start = 0;
    for ( int i = 0; i < rank; i++ )
        start += crits_proc[i];

    Dot * local_crits = NULL;
    if ( n_found > 0 )
    {
        local_crits = new Dot[n_found];

        for ( int i = 0; i < n_found; i++ )
        {
            Dot * p = found_criticals.back();
            local_crits[i].size = p->size;
            local_crits[i].lb = p->lb;
            local_crits[i].ub = p->ub;
            local_crits[i].centroid_index = p->centroid_index;
            local_crits[i].point_index = p->point_index;

            for ( int j = 0; j < D; j++ )
                local_crits[i].coords[j] = p -> coords[j];

            found_criticals.pop_back();
        }
    }

    if ( n_crit > 0 )
    {
        // process 0 create the displacement array for all the other processes to gather all the criticals
        int crit_disp[size];
        MPI_Gather( &start, 1, MPI_INT, crit_disp, 1, MPI_INT, 0, MPI_COMM_WORLD );

        // gather all the criticals in a single array ( "criticals" )
        MPI_Gatherv( local_crits, n_found, dot_type, criticals, crits_proc, crit_disp, dot_type, 0, MPI_COMM_WORLD );
    }

    if ( local_crits != NULL )
        delete[] local_crits;
}

bool MPIHamerly::double_eq( double f1, double f2 )
{
    bool upper_cond = f1 < f2 + EPS;
    bool lower_cond = f1 > f2 - EPS;

    return upper_cond && lower_cond;
}

void MPIHamerly::check_equal( lloyd * l )
{   
    if ( rank == 0 )
    {
        bool eq = true;

        for ( int i = 0; i < n && eq; i++ )
        {
            Point * lp = l -> get_point( i );
            Dot * hp = &points[i];

            if ( lp->centroid_index != hp->centroid_index ) eq = false;
        }

        for ( int i = 0; i < k && eq; i++ )
        {
            Point * lc = l -> get_centroid( i );
            Dot * hc = &centroids[i];

            for ( int j = 0; j < D && eq; j++ )
            {
                if ( !double_eq(lc->get_coord( j ), hc->coords[j]) )
                    eq = false;
            }
        }

        if ( eq ) cout << "equal";
        else cout << "not equal";
    }
}

void MPIHamerly::start_time()
{
    if ( rank == 0 )
        start = chrono::high_resolution_clock::now();
}

void MPIHamerly::elapsed()
{
    if ( rank == 0 )
    {
        chrono::high_resolution_clock::time_point end = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed = chrono::duration_cast<chrono::duration<double>>(end - start);
        double time = elapsed.count();

        cout << time << endl;
    }
}

void MPIHamerly::export_data( string filename )
{
    if ( rank == 0 )
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
        {
            for ( int j = 0; j < D; j++ )
            {
                file << centroids[i].coords[j] << ", ";
            }
            
            file << centroids[i].ub << ", ";
            file << centroids[i].lb << ", ";
            file << centroids[i].centroid_index << endl;
        }

        // white all the points
        file << "Points" << endl;
        for ( int i = 0; i < N; i++ )
        {
            for ( int j = 0; j < D; j++ )
            {
                file << points[i].coords[j] << ", ";
            }
            
            file << points[i].ub << ", ";
            file << points[i].lb << ", ";
            file << points[i].centroid_index << endl;
        }

        file.close();
    }
}
