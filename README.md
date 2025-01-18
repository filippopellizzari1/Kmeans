
# Parallel Hamerly's clustering algorithm
## C++ Parallel K-Means clustering using both OpenMP and MPI framework

This project is intended to provide 3 main classes that can be used to perform a K-Means clustering in 3 different ways.

- Lloyd class: perform the standard k-means clusterring algorithm
- Hamerly class: perform the Hamerly's algorithm in parallel using OpenMP
- MPIHamerly class: perform Hamerly's algorithm in parallel using MPI

In addition, there are some utility classes

- Point class: used by the other classes to represent a point in the dataset
- Dataset class: can be used to import a dataset from a file or generate a random one


## The settings.h file

The settings file is not required when working with the classes expect for MPIHamerly. The latter needs to access it to initiate the MPI_Datatype of the struct Dot which, in this class, is used instead of the class Point to represent the points and centroids. 

Therefore know that to initiate the MPIHamerly class, first modify the settings file to fit the dataset. For all the other classes, using the settings file is not required as the dataset dimensions can be given as the constructor parameter.

## Usage
After downloading the project, copy the classes folder in your project, then in you c++ file include the classes you want to use.

### Dataset class
To use the dataset class first initiate it by passing as parameters the number of points (N), the number of centroids (K) and the data dimensionality (D), then use one of the method to handle the dataset.

```cpp
void rnd_dataset( double min, double max )
void export_dataset( string points_file, string centroids_file )
void load_dataset( string points_file, string centroids_file )
double * get_point( int i )
double get_point_coord( int i , int j )
double * get_centroid( int i )
double get_centroid_coord( int i, int j )
```

The methods `get_point()` and `get_ceentroid()` return an array containing the coordinates of the selected point or centroids.

To reproduce the results shown in the report, the used dataset can be load from the datasets folder.

### Lloyd class
To use the lloyd (sequential) algorithm, initiate the class passing N, K, D. Then initiate the points and centroids using the methods 

```cpp
void insert_point( int i, const Point * p )
void insert_centroid( int i, const Point * c )
```

To initaite a Point pass as parameter the dimensionality and the double array containing the coordinates of the point.

Note: the dataset class can be used to load a dataset (or generate one) and then initiate the other classes (lloyd, hamerly, MPIHamerly) by extracting the points and passing them as parameter for the points and centroids init methods.

At this point, use the methods below to perform the classification for an arbitrary number of time.

```cpp
void update_centroids();
void assign_points();
```

### Hamerly (OMP) class
To use Hamerly's OMP version initiate the class hamerly passing as parameter N, K, D and the number of threads to use T.

Then, initiate the points and run the classification as in the lloyd class.

### MPIHamerly class
To use Hamerly's MPI version initiate the class MPIHamerly passing the address of the argc and argv variable of your main. 

Then, initiate the points using the methods 
```cpp
void init_point(int i, double * coords);
void init_centroid(int i, double * coords);
```

Finally perform the classification using the same lloyd methods.

## Compile and run the code
### OMP version
To compile the OMP version, use the g++ command with the -fopenmp option and the path to the required .cpp files
```console
g++ yourfile.cpp ./classes/hamerly/hamerly.cpp ./classes/lloyd/lloyd.cpp ./classes/point/point.cpp ./classes/dataset/dataset.cpp -fopenmp
```
To run it, execute the generated file
```console
./outfile
```
### MPI version
To compile an MPI code there are two possible ways. The first one is to use g++ with the option -lmpich, the path to the mpi library and the required .cpp files

```console
g++ yourfile.cpp ./classes/MPIHamerly/MPIHamerly.cpp ./classes/dataset/dataset.cpp -I /path/to/mpi.h -L /path/to/mpi -lmpich 
```

The second one is to use the command mpicxx 

```console
mpicxx yourfile.cpp ./classes/MPIHamerly/MPIHamerly.cpp
```

To run it use the command 
```
mpirun -np <num_processes> ./outfile
```

#### Windows 
To run MPI on windows it is possible to install the Microsoft MPI library (msmpi). To compile using msmpi you can use the g++ command above changing -lmpich to -lmsmpi

To run the code use
```
mpiexec –n <num_processes> ./<executable file>
```
