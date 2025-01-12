#ifndef __HAMERLYMPI_H__
#define __HAMERLYMPI_H__

#include <mpi.h>
#define DIMENSION 3
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

class hamerlyMPI
{
    int size, n, k;

    // attributes used only by process 0
    bool first_ass;
    Dot * points;
    Dot * criticals;
    int * items_per_process;
    int * items_disp_per_process;
    int * cent_per_process;
    int * cent_disp_per_process;
    double * average_per_class;
    double * relative_average_per_class;
    int * points_per_class;
    int * relative_points_per_class;
    int * n_found_crit_disp;
    int tot_criticals;

    // attributes used by all processes
    Dot * centroids;
    int * n_criticals;
    double * cdistances;
    double maxdist, secmaxdist;
    double maxdist_index;
    Dot * local_points;
    Dot * local_centroids;
    Dot * local_criticals;
    double * local_cdistances;
    double * local_average_per_class;
    int * local_points_per_class;

    MPI_Datatype dot_type;

    void refresh_assignation();
    void first_assignation();
    void update_criticals();
    double distanceA( Dot * p1, Dot * p2 );
    double distanceB( Dot * p1, Dot * p2 );
    public:
        int rank;
        
        static bool double_eq( double f1, double f2 );

        hamerlyMPI( int _n, int _k, int * argc, char *** argv );
        ~hamerlyMPI();
        
        void init_point( int i, double * coords, int dim );
        void init_centroid( int i, double * coords, int dim );

        void assign_points();
        void update_centroids();

        void checkequal( lloyd * l );
};

#endif
