#ifndef __MPIHAMERLY_H__
#define __MPIHAMERLY_H__

#include <mpi.h>
#define DIMENSION 4
#define EPS 0.0001

class lloyd;
struct Dot
{
    double coords[DIMENSION];
    int size;
    double lb;
    double ub;
    int centroid_index;
    int point_index;
};

class MPIHamerly
{
    int size, n, k, rank;
    MPI_Datatype dot_type;

    // variable in common between all processes
    Dot * centroids;
    int n_crit; // number of criticals found in the last iteration
    bool first_ass;
    int proc_items; // number of items (centroids or points) to iterate
    Dot * local_items; // local array containing the items ( points or centroids ) that have to be iterated
    double * local_avg_class; // sum of coordinates of points in the same calss ( length = k x DIMENSION )
    int * local_points_class; // number of points the same class ( length = k )
    double * local_cdistances;
    double * cdistances; // distances moved by each centroid during the last iteration
    double maxdist, secmaxdist; // max and second max distances moved the centroids during the last iterations
    int maxdist_index; // index of the centroid which moved the most during last iteration
    int * crits_proc; // number of cirticals found by each process during the last iteration

    // Variables used only from process 0
    Dot * points;
    Dot * criticals;
    int * items_per_proc; // this and the following variable are used to store the number of item ( points or centorids ) to assign to every process and the displacement of where they have to start to read them (displacement)
    int  * disp_proc;
    int * points_class; // number of points the same class ( length = k )
    double * avg_class; // sum of coordinates of points in the same calss ( length = k x DIMENSION )
    int * rel_points_class; // number of points the same class ( length = k )
    double * rel_avg_class; // sum of coordinates of points in the same calss ( length = k x DIMENSION )

    void refresh_assignation();
    void first_assignation();
    void update_criticals();
    void splitItems( int nItems );
    void preAssignation_init(); // initiate the variables local_avg_class and local_points_class before a new point assignation is made
    double distance( Dot * p1, Dot * p2 );
    int closest_centroid( Dot * p, double * mindist, double * secmindist );
    bool double_eq( double f1, double f2 );
    public:
        MPIHamerly( int _n, int _k, int * argc, char *** argv );
        ~MPIHamerly();

        void init_point(int i, double * coords);
        void init_centroid(int i, double * coords);
        
        void assign_points();
        void update_centroids();

        void print();
        void check_equal( lloyd * l );
};

#endif
