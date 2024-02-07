/*****************************************************************
* Filename: ModellingHeatFlowinAFusionReactor
* Description: This program solves the partial differential equation used to
* model the heat flow around the first wall of a tokamak reactor. This is
* achieved by transforming the 2D grid of coordinates into a 1D grid with
* periodic boundary conditions taken into account. A coefficient matrix for
* this 1D grid is then generated and then it's corresponding band matrix is 
* produced. This is inverted to get the temperatures at every point in the 
* grid at time t0. The process is then repeated to produce a solution for 
* a prespecified. To run the file have an input.txt and a coefficient.txt file in the
* same directory as this file.
* DATE 31/01/2023
* Author: Kaylen Smith Darnbrook
* Prior to compilation execute following lines on nenneke:
* module purge
* module load intel impi imkl
* Then:
* Compile:  gcc -g -std=c99 -Wall -o HFTAT HeatFlowThroughATorus.c -lm -lmkl -liomp5
* Run: ./Assignment5
******************************************************************/

//Import the relevant libraries
#include <stdlib.h>
#include <stdio.h> 
#include <math.h>
#include <mkl_lapacke.h>
#include <string.h>

//Define the pi constant
#define M_PI 3.14159265358979323846

/* The 'two-dimensional banded matrix example code that is used as a template is sourced from
https://moodle.warwick.ac.uk/course/view.php?id=51816
*/

//Define the banded matrix datatype using C structures
struct band_mat{
  long ncol;        /* Number of columns in band matrix            */
  long nbrows;      /* Number of rows (bands in original matrix)   */
  long nbands_up;   /* Number of bands above diagonal           */
  long nbands_low;  /* Number of bands below diagonal           */
  double *array;    /* Storage for the matrix in banded format  */
  /* Internal temporary storage for solving inverse problem */
  long nbrows_inv;  /* Number of rows of inverse matrix   */
  double *array_inv;/* Store the inverse if this is generated */
  int *ipiv;        /* Additional inverse information         */
};

//Declare the banded matrix as a datastructure.
typedef struct band_mat band_mat;

//Function Prototypes.
long read_input(char*, long*, long*, double*, long*); // This function will read in the input.txt
long read_coefficients(char*, long, long, double*, double*, double*, double*, double* , double*, double*, double, double); // This function will read in the coefficients.txt
long output_file(char*, long, double*, double*, double*, double*); //This function will output the displacements to an output text.
int init_banded_matrix(band_mat*, long, long, long); //This function will initialise the band matrix. 
void finalise_banded_matrix(band_mat*); //Frees the heap memory to prevent memory leakage.
double* getp(band_mat*, long, long); //Returns a pointer that points to a location in the banded matrix using the row and column indices of the full matrix.
double getv(band_mat*, long, long); //Returns the corresponding value in the banded matrix to the pointer given by getp.
void setv(band_mat*, long, long, double); //Sets the corresponding value in the banded matrix to the pointer given by getp.
int solve_Ax_eq_b(band_mat*, double*, double*); //Solves the matrix equation Ax=b for a matrix stored in banded format with x and b being real arrays.
int printmat(band_mat*); //Outputs the b_matrix.

/* This function reads in the number of grid points for the theta and Zeta axes, the final time for the time evolution mode and the minimum number of timesteps to be 
*  to be used to model the heat flow in the simulation.
*/
long read_input(char* filename, long* Number_Theta_Grid_Points, long* Number_Zeta_Grid_Points, double* final_time, long* Num_timesteps) {

    FILE* fptr = fopen(filename,"r"); //Opens the file in read mode and returns its address in memory.
    printf("he");
    if (fptr == NULL) return 1; //If the file could not be read or does not exist the function returns 1.
    printf("Hello");
    if (4!=fscanf_s(fptr,"%ld %ld %lf %ld", Number_Theta_Grid_Points, Number_Zeta_Grid_Points, final_time, Num_timesteps)) { //If the number of elements are not consistent with the given scheme the file returns 1.
        return 1;
    }
    fclose(fptr); //Closes the file.
    return 0;

}

/*This function reads in the three (2 are equal) conductivity coefficients, the source term and the reservoir term used
* in modelling the heat flow in the first wall of the tokamak reactor. 
*/
long read_coefficients(char* filename_2, long Number_Theta_Grid_Points, long Number_Zeta_Grid_Points, double* Theta_Positions, double* Zeta_Positions, double* Q_11, double* Q_22, double* Q_12, double* Source, double* Resvoir, double dTheta, double dZeta) {
    
    //Creates dummies to read in the three conductivity coefficients, the source term and reservoir term
    double a;
    double b; 
    double c;
    double d;
    double e;
    FILE* fptr2 = fopen(filename_2, "r"); //Opens the file in read mode and returns its address in memory.
    if (fptr2 == NULL) return 1; //If the file could not be read or does not exist the function returns 1.
    for(int increment_i = 0; increment_i < Number_Theta_Grid_Points; increment_i++){// Iterates over the arrays writing data to the files.
        for (long increment_j = 0; increment_j < Number_Zeta_Grid_Points; increment_j++){
            long L_increment = increment_i * Number_Zeta_Grid_Points + increment_j;
            if (5!=fscanf_s(fptr2,"%lf %lf %lf %lf %lf\n", &a, &b, &c, &d, &e)) { //If the number of elements are not consistent with the given scheme the file returns 1.
                return 1;
            }
            
            //Writes the Theta_Positions, Zeta_Positions, conductivities, sources and resvoir terms to the arrays.
            Theta_Positions[L_increment] = (2.0 * M_PI * increment_i)/(Number_Theta_Grid_Points);
            Zeta_Positions[L_increment] = (2.0 * M_PI * increment_j)/(Number_Zeta_Grid_Points);
            Q_11[L_increment] = a;
            Q_22[L_increment] = b;
            Q_12[L_increment] = c;
            Source[L_increment] = d;
            Resvoir[L_increment] = e;

        }
    }
    
   fclose(fptr2); //Closes the file.
   return 0;
} 

/*Writes the Temperatures T over time t to the output text file */
long output_file(char* filename, long Num_Of_Grid_Points, double* Time, double* Theta_Positions, double* Zeta_Positions, double* Temperatures) {
    
    FILE* fptr = fopen(filename, "w");
    if (fptr == NULL) return 1;
      for(long increment = 0; increment < Num_Of_Grid_Points; increment++) {
        fprintf(fptr, "%lf %lf %lf %lf\n", Time[increment], Theta_Positions[increment], Zeta_Positions[increment], Temperatures[increment]);
      }
    fclose(fptr);
    return 0;
}

/*Initialises a null banded matrix with a certain size, allocates the necessary memory and sets the parameters.
    The equivalent in OOP is the class constructor function*/
int init_band_mat(band_mat *bmat, long nbands_lower, long nbands_upper, long n_columns) {
  bmat->nbrows = nbands_lower + nbands_upper + 1;// Sets the number of rows the banded matrix will have. 
  bmat->ncol   = n_columns; // Sets the number of columns the banded matrix will have.
  bmat->nbands_up = nbands_upper; // Sets the number of rows above the diagonal.
  bmat->nbands_low= nbands_lower; // Sets the number of rows below the diagonal.
  bmat->array      = (double *) malloc(sizeof(double)*bmat->nbrows*bmat->ncol); // Initilises the banded matrix using heap memory.
  bmat->nbrows_inv = bmat->nbands_up*2 + bmat->nbands_low + 1; // Initialises the number of rows the inverse matrix will have.
  bmat->array_inv  = (double *) malloc(sizeof(double)*(bmat->nbrows+bmat->nbands_low)*bmat->ncol); // Initialises the inverse banded matrix.
  bmat->ipiv       = (int *) malloc(sizeof(int)*bmat->ncol); // Initialises the store of additional information
  
  //Checks whether the array and inverse array have been properly initialised and if terminates the program.
  if (bmat->array==NULL||bmat->array_inv==NULL) {
    return 0;
  }

  /* Initialise array to zero */
  long i;
  for (i=0;i<bmat->nbrows*bmat->ncol;i++) {
    bmat->array[i] = 0.0; //Sets each element in the array to zero.
  }
  return 1;
};

/* Finalise function: should free memory as required */
void finalise_band_mat(band_mat *bmat) {
  free(bmat->array);
  free(bmat->array_inv);
  free(bmat->ipiv);
}

/*Returns a pointer that points to a location in the band matrix using
    * using the corresponding row and column indexes of the full matrix.
    */
double *getp(band_mat *bmat, long row, long column) {
  int bandno = bmat->nbands_up + row - column;
  
  //Checks whether the element being queried exists in the banded_matrix.
  if(row<0 || column<0 || row>=bmat->ncol || column>=bmat->ncol ) {
    printf("Indexes out of bounds in getp: %ld %ld %ld \n",row,column,bmat->ncol);
    exit(1);
  }
  return &bmat->array[bmat->nbrows*column + bandno]; //Accesses the element in the banded matrix using the banded matrix indexing rule.
}


/*Returns the corresponding band matrix value from the pointer acquired from getp using the 
    rows and column indexes of the full matrix.*/
double getv(band_mat *bmat, long row, long column) {
  return *getp(bmat,row,column); //Gets the pointer using getp and deferences it.
}

/*Sets an element of a band matrix to a desired value.
    using a getp to obtain the location in band matrix and then dereferencing.*/
void setv(band_mat *bmat, long row, long column, double val) {
  *getp(bmat,row,column) = val+*getp(bmat,row,column);
}

/* This function outputs the elements of the banded matrix */
int printmat(band_mat *bmat) {
  long i,j;
  for(i=0; i<bmat->ncol;i++) {
    for(j=0; j<bmat->nbrows; j++) {
      printf("%ld %ld %g \n",i,j,bmat->array[bmat->nbrows*i + j]);
    }
  }
 return 0;
}

/*Check that a grid point has valid coordinates */
int is_valid(long i, long j, long Num_i_Grid_Points, long Num_j_Grid_Points) {
  return (i>=0)&&(i<Num_i_Grid_Points)&&(j>=0)&&(j<Num_j_Grid_Points);
}

/*Solve the equation Ax=b for a matrix stored in banded format 
    with x and b being real arrays*/
int solve_Ax_eq_b(band_mat *bmat, double *x, double *b) {
   //Copies the bmat array into the temporary store
  int i,bandno;
  for(i=0;i<bmat->ncol;i++) { //Loop over the rows 
    for (bandno=0;bandno<bmat->nbrows;bandno++) { // loop over the columns
      //Transfer contents of the bmat array to the inverse bmat array
      bmat->array_inv[bmat->nbrows_inv*i+(bandno+bmat->nbands_low)] = bmat->array[bmat->nbrows*i+bandno];
    }
    x[i] = b[i]; //Assigns the source to the x values.
  }

  //Below uses LAPACKE to invert and solve the banded matrix.
  long nrhs = 1;
  long ldab = bmat->nbands_low*2 + bmat->nbands_up + 1;
  long info = 0;
  info = LAPACKE_dgbsv( LAPACK_COL_MAJOR, bmat->ncol, bmat->nbands_low, bmat->nbands_up, nrhs, bmat->array_inv, ldab, bmat->ipiv, x, bmat->ncol);
  printf("info");
  return info;
}

/* Maps the index of a grid point to it's periodic index*/
long periodic_boundary(long index, long Num_Grid_Points) {
    if ( ((double) index) < ((double) Num_Grid_Points/2)) {
        return 2*index;
    }
    else{
        return 2*(Num_Grid_Points-index)-1;
    }
}

//Maps the index of a 2D grid point to 1D periodic index
long index_periodic_boundary_forward_converter(long index_i,long index_j, long Num_i_Grid_Points, long Num_j_Grid_Points) {
    return periodic_boundary(index_i, Num_i_Grid_Points) * Num_j_Grid_Points + periodic_boundary(index_j,Num_j_Grid_Points);
}

//Maps the index of 1D periodic index back to a 2D index
long index_periodic_boundary_backward_converter(long Periodic_L_index, long Num_i_Grid_Points, long Num_j_Grid_Points){
    long index_i;
    long index_j;
    long L_index;

    long Periodic_i_index = Periodic_L_index / Num_j_Grid_Points;
    long Periodic_j_index = Periodic_L_index % Num_j_Grid_Points;

    long i_index_1 = Periodic_i_index/2;
    long i_index_2 = Num_i_Grid_Points - (Periodic_i_index+1)/2;

    long j_index_1 = Periodic_j_index/2;
    long j_index_2 = Num_j_Grid_Points - (Periodic_j_index+1)/2;

    if (index_periodic_boundary_forward_converter(i_index_1, j_index_1, Num_i_Grid_Points, Num_j_Grid_Points) == Periodic_L_index) {
        index_i = i_index_1;
        index_j = j_index_1;
    }
    else if (index_periodic_boundary_forward_converter(i_index_1, j_index_2, Num_i_Grid_Points, Num_j_Grid_Points) == Periodic_L_index){
        index_i = i_index_1;
        index_j = j_index_2;
    }
    else if (index_periodic_boundary_forward_converter(i_index_2, j_index_1, Num_i_Grid_Points, Num_j_Grid_Points) == Periodic_L_index){
        index_i = i_index_2;
        index_j = j_index_1;
    }
    else if (index_periodic_boundary_forward_converter(i_index_2, j_index_2, Num_i_Grid_Points, Num_j_Grid_Points) == Periodic_L_index){
        index_i = i_index_2;
        index_j = j_index_2;
    }
    else {
        index_i = 0;
        index_j = 0;
    } 
    L_index = index_i * Num_j_Grid_Points + index_j;
    return L_index;
}

//Creates the banded matrix to invert using the heat equation derived for the first wall of the tomak reactor.
void MatrixCreator(long Num_i_Points, long Num_j_Points, double dTime, double di, double dj, double* Q_11, double* Q_22, double* Q_12, double* Resvoir, band_mat banded_coefficient_matrix) {
    
  //Defines variables used to enforce periodicity.
  long index_Plus_i = 0;
  long index_Minus_i = 0;
  long index_Plus_j = 0;
  long index_Minus_j = 0;

  for (long index_i = 0; index_i < Num_i_Points; index_i++) {
      if (index_i == 0) {
        index_Minus_i = Num_i_Points - 1;
      }
      else {
        index_Minus_i = index_i - 1;
      }
      if (index_i == Num_i_Points - 1){
        index_Plus_i = 0;
      }
      else {
        index_Plus_i = index_i + 1;
      }

      for (long index_j = 0; index_j < Num_j_Points; index_j++) {
    
        if (index_j == 0) {
        index_Minus_j = Num_j_Points - 1;
        }
        else {
        index_Minus_j = index_j - 1;
        }
        if (index_j == Num_j_Points - 1){
        index_Plus_j = 0;
         }
        else {
        index_Plus_j = index_j + 1;
        }

        //Entering the coefficients into the null matrix
        long index_current = index_periodic_boundary_forward_converter(index_i, index_j, Num_i_Points, Num_j_Points);
    
        if (is_valid(index_Plus_i, index_Plus_j, Num_i_Points, Num_j_Points)) {
          setv(&banded_coefficient_matrix, index_current, index_periodic_boundary_forward_converter(index_Plus_i, index_Plus_j, Num_i_Points, Num_j_Points), (-(1.0)/((2.0)*(di)*(dj))) * (Q_12[index_periodic_boundary_forward_converter(index_i, index_j, Num_i_Points, Num_j_Points)]));
        }
        if (is_valid(index_Plus_i, index_j, Num_i_Points, Num_j_Points)) {
          setv(&banded_coefficient_matrix, index_current, index_periodic_boundary_forward_converter(index_Plus_i, index_j, Num_i_Points, Num_j_Points), ((1.0)/(di)) * ((Q_11[index_periodic_boundary_forward_converter(index_Minus_i, index_j, Num_i_Points, Num_j_Points)] - Q_11[index_periodic_boundary_forward_converter(index_Plus_i, index_j, Num_i_Points, Num_j_Points)])/(4.0*(di)) - Q_11[index_periodic_boundary_forward_converter(index_i,index_j,Num_i_Points,Num_j_Points)]/(di) + (Q_12[index_periodic_boundary_forward_converter(index_i, index_Minus_j,Num_i_Points, Num_j_Points)] - Q_12[index_periodic_boundary_forward_converter(index_i, index_Plus_j,Num_i_Points, Num_j_Points)])/(4.0*(dj))));
        }
        if (is_valid(index_i, index_Plus_j, Num_i_Points, Num_j_Points)) {
          setv(&banded_coefficient_matrix, index_current, index_periodic_boundary_forward_converter(index_i, index_Plus_j, Num_i_Points, Num_j_Points), ((1.0)/(dj)) * ((Q_22[index_periodic_boundary_forward_converter(index_i, index_Minus_j, Num_i_Points, Num_j_Points)] - Q_22[index_periodic_boundary_forward_converter(index_i, index_Plus_j, Num_i_Points, Num_j_Points )])/(4.0*(dj)) - Q_22[index_periodic_boundary_forward_converter(index_i,index_j,Num_i_Points,Num_j_Points)]/(dj) + (Q_12[index_periodic_boundary_forward_converter(index_Minus_i, index_j, Num_i_Points, Num_j_Points)]- Q_12[index_periodic_boundary_forward_converter(index_Plus_i, index_j, Num_i_Points,Num_j_Points)])/(4.0*(di))));
        }
        if (is_valid(index_i, index_j, Num_i_Points, Num_j_Points)) {
          setv(&banded_coefficient_matrix, index_current, index_periodic_boundary_forward_converter(index_i, index_j, Num_i_Points, Num_j_Points), ((2.0)) * ((Q_11[index_periodic_boundary_forward_converter(index_i,index_j,Num_i_Points,Num_j_Points)])/(pow(di,2.0)) + (1.0 + 0.5*Resvoir[index_periodic_boundary_forward_converter(index_i,index_j,Num_i_Points,Num_j_Points)]*dTime)/(dTime) + Q_22[index_periodic_boundary_forward_converter(index_i,index_j,Num_i_Points,Num_j_Points)]/(pow(dj,2.0))));
        }
        if (is_valid(index_Minus_i, index_j, Num_i_Points, Num_j_Points)) {
          setv(&banded_coefficient_matrix, index_current, index_periodic_boundary_forward_converter(index_Minus_i, index_j, Num_i_Points, Num_j_Points), ((1.0)/(di)) * ((Q_11[index_periodic_boundary_forward_converter(index_Plus_i, index_j, Num_i_Points, Num_j_Points)] - Q_11[index_periodic_boundary_forward_converter(index_Minus_i, index_j, Num_i_Points, Num_j_Points)])/(4.0*di) - Q_11[index_periodic_boundary_forward_converter(index_i, index_j, Num_i_Points, Num_j_Points)]/(di) + (Q_12[index_periodic_boundary_forward_converter(index_i, index_Plus_j, Num_i_Points, Num_j_Points)] - Q_12[index_periodic_boundary_forward_converter(index_i, index_Minus_j, Num_i_Points, Num_j_Points)])/(4.0*(dj))));
        }
        if (is_valid(index_i, index_Minus_j, Num_i_Points, Num_j_Points)) {
          setv(&banded_coefficient_matrix, index_current, index_periodic_boundary_forward_converter(index_i, index_Minus_j, Num_i_Points, Num_j_Points), ((1.0)/(dj)) * ((Q_22[index_periodic_boundary_forward_converter(index_i, index_Plus_j, Num_i_Points, Num_j_Points)] - Q_22[index_periodic_boundary_forward_converter(index_i, index_Minus_j, Num_i_Points, Num_j_Points)])/(4.0*dj) - Q_22[index_periodic_boundary_forward_converter(index_i, index_j, Num_i_Points, Num_j_Points)]/(dj) + (Q_12[index_periodic_boundary_forward_converter(index_Plus_i, index_j, Num_i_Points, Num_j_Points)] - Q_12[index_periodic_boundary_forward_converter(index_Minus_i, index_j, Num_i_Points, Num_j_Points)])/(4.0*(di))));
        }
        if (is_valid(index_Minus_i, index_Plus_j, Num_i_Points, Num_j_Points)) {
          setv(&banded_coefficient_matrix, index_current, index_periodic_boundary_forward_converter(index_Minus_i, index_Plus_j, Num_i_Points, Num_j_Points), ((1.0)/(2.0*(di)*(dj)))*(Q_12[index_periodic_boundary_forward_converter(index_i, index_j, Num_i_Points, Num_j_Points)]));
        }
        if (is_valid(index_Plus_i, index_Minus_j, Num_i_Points, Num_j_Points)) {
          setv(&banded_coefficient_matrix, index_current, index_periodic_boundary_forward_converter(index_Plus_i, index_Minus_j, Num_i_Points, Num_j_Points), ((1.0)/(2.0*(di)*(dj)))*(Q_12[index_periodic_boundary_forward_converter(index_i, index_j, Num_i_Points, Num_j_Points)]));
        }
        if (is_valid(index_Minus_i, index_Minus_j, Num_i_Points, Num_j_Points)) {
          setv(&banded_coefficient_matrix, index_current, index_periodic_boundary_forward_converter(index_Minus_i, index_Minus_j, Num_i_Points, Num_j_Points), ((-1.0)/(2.0*(di)*(dj)))*(Q_12[index_periodic_boundary_forward_converter(index_i, index_j, Num_i_Points, Num_j_Points)]));
        }
    }
  }
}

//Generates the values of the product coefficient matrix and current temperature solution at the previous current time step for use in predicting the next timestep
void phi_generator(long Num_i_Points, long Num_j_Points, double dTime, double di, double dj, double* Q_11, double* Q_22, double* Q_12, double* Resvoir, double* phi, double* Temperatures_Solution) {
  
  //Defines variables used to enforce periodicity.
  long index_Plus_i;
  long index_Minus_i;
  long index_Plus_j;
  long index_Minus_j;
  
  //Enforce periodicity.
  for (long index_i = 0; index_i < Num_i_Points; index_i++) {
    if (index_i == 0) {
        index_Minus_i = Num_i_Points - 1;
      }
    else {
        index_Minus_i = index_i - 1;
      }
      if (index_i == Num_i_Points - 1){
        index_Plus_i = 0;
      }
    else {
        index_Plus_i = index_i + 1;
      }
    for (long index_j = 0; index_j < Num_j_Points; index_j++) {
        
        if (index_j == 0) {
        index_Minus_j = Num_j_Points - 1;
        }
        else {
        index_Minus_j = index_j - 1;
        }
        if (index_j == Num_j_Points - 1){
        index_Plus_j = 0;
         }
        else {
        index_Plus_j = index_j + 1;
        }
        //Generate the solution to the matrix equation at current time step.
        double contribution_0 = ((1.0)/(2.0*(di)*(dj)))*(Q_12[index_periodic_boundary_forward_converter(index_i,index_j, Num_i_Points, Num_j_Points)]) * Temperatures_Solution[index_periodic_boundary_forward_converter(index_Plus_i, index_Plus_j, Num_i_Points, Num_j_Points)];
        double contribution_1 = ((1.0)/(di)) * ((Q_11[index_periodic_boundary_forward_converter(index_Plus_i, index_j, Num_i_Points, Num_j_Points)] - Q_11[index_periodic_boundary_forward_converter(index_Minus_i, index_j, Num_i_Points, Num_j_Points)])/(4.0*di) + (Q_11[index_periodic_boundary_forward_converter(index_i, index_j, Num_i_Points, Num_j_Points)])/(di) + (Q_12[index_periodic_boundary_forward_converter(index_i, index_Plus_j, Num_i_Points, Num_j_Points)] - Q_12[index_periodic_boundary_forward_converter(index_i, index_Minus_j, Num_i_Points, Num_j_Points)])/(4.0*dj)) * Temperatures_Solution[index_periodic_boundary_forward_converter(index_Plus_i, index_j, Num_i_Points, Num_j_Points)];
        double contribution_2 = ((1.0)/(dj)) * ((Q_22[index_periodic_boundary_forward_converter(index_i, index_Plus_j, Num_i_Points, Num_j_Points)] - Q_22[index_periodic_boundary_forward_converter(index_i, index_Minus_j, Num_i_Points, Num_j_Points)])/(4.0*dj) + (Q_22[index_periodic_boundary_forward_converter(index_i, index_j, Num_i_Points, Num_j_Points)])/(dj) + (Q_12[index_periodic_boundary_forward_converter(index_Plus_i, index_j, Num_i_Points, Num_j_Points)] - Q_12[index_periodic_boundary_forward_converter(index_Minus_i, index_j, Num_i_Points, Num_j_Points)])/(4.0*di)) * Temperatures_Solution[index_periodic_boundary_forward_converter(index_i, index_Plus_j, Num_i_Points, Num_j_Points)];
        double contribution_3 = (-2.0)*(((Q_11[index_periodic_boundary_forward_converter(index_i, index_j, Num_i_Points, Num_j_Points)]))/(pow(di,2.0)) + (-1.0 + 0.5*Resvoir[index_periodic_boundary_forward_converter(index_i,index_j, Num_i_Points, Num_j_Points)]*dTime)/(dTime) + (Q_22[index_periodic_boundary_forward_converter(index_i,index_j, Num_i_Points, Num_j_Points)])/(pow(dj,2.0))) * Temperatures_Solution[index_periodic_boundary_forward_converter(index_i, index_j, Num_i_Points, Num_j_Points)];
        double contribution_4 = ((1.0)/(di)) * ((Q_11[index_periodic_boundary_forward_converter(index_Minus_i, index_j, Num_i_Points, Num_j_Points)] - Q_11[index_periodic_boundary_forward_converter(index_Plus_i, index_j, Num_i_Points, Num_j_Points)])/(4.0*di) + (Q_11[index_periodic_boundary_forward_converter(index_i, index_j, Num_i_Points, Num_j_Points)])/(di) + (Q_12[index_periodic_boundary_forward_converter(index_i, index_Minus_j, Num_i_Points, Num_j_Points)] - Q_12[index_periodic_boundary_forward_converter(index_i, index_Plus_j, Num_i_Points, Num_j_Points)])/(4.0*dj)) * Temperatures_Solution[index_periodic_boundary_forward_converter(index_Minus_i, index_j, Num_i_Points, Num_j_Points)];
        double contribution_5 = ((1.0)/(dj)) * ((Q_22[index_periodic_boundary_forward_converter(index_i, index_Minus_j, Num_i_Points, Num_j_Points)] - Q_22[index_periodic_boundary_forward_converter(index_i, index_Plus_j, Num_i_Points, Num_j_Points)])/(4.0*dj) + (Q_22[index_periodic_boundary_forward_converter(index_i, index_j, Num_i_Points, Num_j_Points)])/(dj) + (Q_12[index_periodic_boundary_forward_converter(index_Minus_i, index_j, Num_i_Points, Num_j_Points)] - Q_12[index_periodic_boundary_forward_converter(index_Plus_i, index_j, Num_i_Points, Num_j_Points)])/(4.0*di)) * Temperatures_Solution[index_periodic_boundary_forward_converter(index_i, index_Minus_j, Num_i_Points, Num_j_Points)];
        double contribution_6 = (-(1.0)/(2.0*(di)*(dj))) * (Q_12[index_periodic_boundary_forward_converter(index_i,index_j, Num_i_Points, Num_j_Points)]) * Temperatures_Solution[index_periodic_boundary_forward_converter(index_Minus_i, index_Plus_j, Num_i_Points, Num_j_Points)];
        double contribution_7 = (-(1.0)/(2.0*(di)*(dj))) * (Q_12[index_periodic_boundary_forward_converter(index_i,index_j, Num_i_Points, Num_j_Points)]) * Temperatures_Solution[index_periodic_boundary_forward_converter(index_Plus_i, index_Minus_j, Num_i_Points, Num_j_Points)];
        double contribution_8 = (((1.0)/(2.0*(di)*(dj)))) * (Q_12[index_periodic_boundary_forward_converter(index_i,index_j, Num_i_Points, Num_j_Points)]) * Temperatures_Solution[index_periodic_boundary_forward_converter(index_Minus_i, index_Minus_j, Num_i_Points, Num_j_Points)];
        
        phi[index_periodic_boundary_forward_converter(index_i, index_j, Num_i_Points, Num_j_Points)] = contribution_0 + contribution_1 + contribution_2 + contribution_3 + contribution_4 + contribution_5 + contribution_6 + contribution_7 + contribution_8;
    }
  }
}

int main(void) {
    /*
    The main function runs three tasks
    1. Input the file and the coefficient file. Store the inputs and variables 
    and the coefficients as arrays.
    2. Use indexing functions to map the two dimensional grid of theta and zeta to a 
    1D grid of temperature values.
    3. Create a coefficient matrix, convert it into a banded matrix and invert it to 
    solve the matrix equation and return a vector of results for time 1.
    4. Loop through steps 2 and 3 for the range of required timesteps to show how the 
    temperature evolves with time.
    5. Return the solution to a text file.
    */
    
    //Define the variables storing the input values.
    long Number_Theta_Grid_Points; //Stores the number of points in the Theta axis.
    long Number_Zeta_Grid_Points; //Stores the number of points in the Zeta Axis.
    double final_time; //Stores the final time of the simulation.
    long Num_Timesteps; // Stores the number of timesteps in the simulation.

    char* Input_Text = "C:\\Users\\kdarn\\OneDrive\\Documents\\Life's Portfolio\\Projects\\NeuralODEProject\\input.txt"; //Stores the filename of the inputs.
    if (read_input(Input_Text,&Number_Theta_Grid_Points, &Number_Zeta_Grid_Points, &final_time, &Num_Timesteps)){ //Reads in the grid points, length and angular frequency
        printf("Input text file read error"); // Outputs whether or not the input text could be read for the purposes of the purpose of debugging 
        return 1;
    }

    //Defines the infinitesimals
    double dTheta = (2.0*M_PI)/((double) Number_Theta_Grid_Points); 
    double dZeta = (2.0*M_PI)/((double) Number_Zeta_Grid_Points);
    double dTime = (final_time)/((double) Num_Timesteps);

    //Define the number of grid points
    long Num_Of_Grid_Points;
    Num_Of_Grid_Points = Number_Theta_Grid_Points * Number_Zeta_Grid_Points;

    //Outputs the Number of Theta Grid_Points, Number of Zeta Grid_Points,
    // the final time to evolve the solution to and the number of timesteps.
    //printf("The number of Theta grid points is %ld, The number of Zeta grid points is %ld, The final timestep is %lf, The number of timesteps is %ld\n", Number_Theta_Grid_Points, Number_Zeta_Grid_Points,final_time,Num_Timesteps);

    //Outputs the Theta and Zeta Infinitesmial
    //printf("The dTheta is %lf. The dZeta is %lf. The dTime is %lf\n", dTheta, dZeta, dTime);

    //This section defines the coefficients, source and resvoir terms used to build the coeffient matrix
    double* Theta_Positions = (double*) malloc(sizeof(double) * Num_Of_Grid_Points); // Stores the Theta positions of the wall.
    double* Zeta_Positions = (double*) malloc(sizeof(double) * Num_Of_Grid_Points); // Stores the Zeta positions of the wall.
    double* Q_11 = (double*) malloc(sizeof(double) * Num_Of_Grid_Points); //Stores the Q_11 Conductivity Coefficient.
    double* Q_22 = (double*) malloc(sizeof(double) * Num_Of_Grid_Points); //Stores the Q_22 Conductivity Coefficient.
    double* Q_12 = (double*) malloc(sizeof(double) * Num_Of_Grid_Points); //Stores the Q_12 Coductivitty Coefficient.
    double* Source = (double*) malloc(sizeof(double) * Num_Of_Grid_Points); //Stores the Source of the heat into the wall.
    double* Resvoir = (double*) malloc(sizeof(double) * Num_Of_Grid_Points); //Stores the Resvoir coefficients by which heat leaves.
    
    if (Theta_Positions == NULL|| Zeta_Positions == NULL || Q_11 == NULL || Q_22 == NULL || Source == NULL || Resvoir == NULL) {
        printf("Allocation error, unable to obtain sufficient memory\n"); //Checks to see whether the computer is able to allocate sufficient dynamic memory.
        return 1;
        }

    //Cleans the arrays to ensure that locations in memory that the pointers point to contain only zeros.
    for (long L_increment = 0; L_increment < Num_Of_Grid_Points; L_increment++) {
        Theta_Positions[L_increment] = 0.0;
        Zeta_Positions[L_increment] = 0.0;
        Q_11[L_increment] = 0.0;
        Q_22[L_increment] = 0.0;
        Q_12[L_increment] = 0.0;
        Source[L_increment] = 0.0;
        Resvoir[L_increment] = 0.0;
    }

    char* Coefficients_Text = "C:\\Users\\kdarn\\OneDrive\\Documents\\Life's Portfolio\\Projects\\NeuralODEProject\\coefficients.txt";
    if (read_coefficients(Coefficients_Text, Number_Theta_Grid_Points, Number_Zeta_Grid_Points, Theta_Positions, Zeta_Positions, Q_11, Q_22, Q_12, Source, Resvoir, dTheta, dZeta)){ //Reads in the conductivity coefficients, the source and the resvoir coefficient. 
        printf("Coefficients text file read error\n"); //Outputs whether or not the coefficient text could be read for the purpose of debugging.
        return 1;
    }
    //In this submission file I have for debugging purposes left the print commands for the conductivity coefficients, 
    //Source and Resvoir coefficients at theta and zeta points used during run time.

    /*
    for (long increment_theta = 0; increment_theta < Number_Theta_Grid_Points; increment_theta++) {
        for (long increment_zeta = 0; increment_zeta < Number_Zeta_Grid_Points; increment_zeta++) {
            long L_increment = increment_theta * Number_Zeta_Grid_Points + increment_zeta;
            printf("At L_increment %ld, the Theta position is %lf, the Zeta position is %lf, the Q_11 conductivity is %lf, the Q_22 conductivity is %lf, the Q_12 conductivity is %lf ", L_increment, Theta_Positions[L_increment], Zeta_Positions[L_increment], Q_11[L_increment], Q_22[L_increment], Q_12[L_increment]);
            printf("At L_increment %ld, the Source is %lf, the Resvoir is %lf\n",L_increment,Source[L_increment], Resvoir[L_increment]);
        }
    }*/

    //Define the number of bands above and below the diagonal
    long nbands_below = 4*Number_Zeta_Grid_Points -1;
    long nbands_above = nbands_below;
    long Num_Iterations = Num_Timesteps;
    

    //Initialises the vectors used to solve the differential equation.
    double* Temperature_Solution = (double*) malloc(sizeof(double) * Num_Of_Grid_Points);
    double* phi = (double*) malloc(sizeof(double) * Num_Of_Grid_Points);
    double* b = (double*) malloc(sizeof(double) * Num_Of_Grid_Points);
    double* Time = (double*) malloc(sizeof(double) * Num_Of_Grid_Points);

    //Pointer validation for the vectors used to solve the differential equation.
    if (Temperature_Solution == NULL || phi == NULL || b == NULL || Time == NULL) {
        printf("Allocation error, unable to obtain sufficient memory\n"); //Checks to see whether the computer is able to allocate sufficient dynamic memory.
        return 1;
        }

    //Initializes phi vector to have zero in all entries.
    for (long L_increment = 0; L_increment < Num_Of_Grid_Points; L_increment++) {
          Temperature_Solution[L_increment] = 0.0;
          phi[L_increment] = 0.0;
          b[L_increment] = 0.0;
          Time[L_increment] = final_time;
          }


    //Creates the for loop that solves the coefficient matrix for each time step.
    for (long index = 0; index <= Num_Iterations; index++) {
      //Initialise the banded matrix
      band_mat bmat;
      init_band_mat(&bmat, nbands_below, nbands_above, Num_Of_Grid_Points);

      //Zeros the Temperature Solution
      for (long L_increment = 0; L_increment < Num_Of_Grid_Points; L_increment++) {
          Temperature_Solution[L_increment] = 0.0;
      
       }

       //Creates the current timestep vector b for the matrix equation
      for (long L_increment = 0; L_increment < Num_Of_Grid_Points; L_increment++) {
          b[L_increment] = phi[L_increment] + 2.0*Source[L_increment];
       }

      //Creates the banded matrix to invert.
      MatrixCreator(Number_Theta_Grid_Points, Number_Zeta_Grid_Points, dTime, dTheta, dZeta, Q_11, Q_22, Q_12, Resvoir, bmat);
      
      printf("Start");
      //Solves the banded matrix.
      //printmat(&bmat);
      int info = 0;
      info = solve_Ax_eq_b(&bmat, Temperature_Solution, b);
      printf("\n info = %d\n", info);
      printf("End");
      //Generates the next phi vector
      phi_generator(Number_Theta_Grid_Points, Number_Zeta_Grid_Points, dTime, dTheta, dZeta, Q_11, Q_22, Q_12, Resvoir, phi, Temperature_Solution);

      //Writes to the array storing the temperatures.
      for (long Periodic_L_increment = 0; Periodic_L_increment < Num_Of_Grid_Points; Periodic_L_increment++) {
          Temperature_Solution[Periodic_L_increment] = Temperature_Solution[Periodic_L_increment];
      }
      
      //Free the memory associated with the banded matrix
      finalise_band_mat(&bmat);
    }
    //prints the temperature at the position theta and zeta at final time t.
    
    for(long increment = 0; increment < Num_Of_Grid_Points; increment++) {
        printf("%lf %lf %lf %lf\n", Time[increment], Theta_Positions[increment], Zeta_Positions[increment], Temperature_Solution[increment]);
      }
    

    //Outputs the temperatures for every theta and zeta on the grid at the final time t.
    char* outputs = "C:\\Users\\kdarn\\OneDrive\\Documents\\Life's Portfolio\\Projects\\NeuralODEProject\\output.txt";
    if (output_file(outputs, Num_Of_Grid_Points, Time, Theta_Positions, Zeta_Positions, Temperature_Solution)) {
        printf("Output text file cannot be written to");
        return 1;
    } 

    printf("end");

    //Free the memory from the stack
    free(Theta_Positions);
    free(Zeta_Positions);
    free(Temperature_Solution);
    free(b);
    free(Time);
    free(phi);
    free(Q_11);
    free(Q_22);
    free(Q_12);
    free(Source);
    free(Resvoir);

    return 0;
}







