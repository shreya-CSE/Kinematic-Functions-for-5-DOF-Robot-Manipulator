#include <stdio.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.1415927
#endif

void rotation_x(double rot_x[4][4], double angle);
void rotation_y(double rot_y[4][4], double angle);
void rotation_z(double rot_z[4][4], double angle);
void translation_x(double Dx[4][4], double q);
void translation_y(double Dy[4][4], double q);
void translation_z(double Dz[4][4], double q);
void matrix_multiplication(double product[4][4], double matrix_1[4][4], double matrix_2[4][4]);


// Rx 
void rotation_x(double rot_x[4][4], double angle){
    double matrix[][4] = { {1, 0, 0, 0}, {0, cos(angle), -sin(angle), 0}, {0, sin(angle), cos(angle), 0}, {0, 0, 0, 1} };
    
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            rot_x[i][j] = matrix[i][j];
        }
    }
}

// Ry
void rotation_y(double rot_y[4][4], double angle) {
    double matrix[][4] = { {cos(angle), 0, sin(angle), 0}, {0, 1, 0, 0}, {-sin(angle), 0, cos(angle), 0}, {0, 0, 0, 1} };
 
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            rot_y[i][j] = matrix[i][j];
        }
    }
}

// Rz 
void rotation_z(double rot_z[4][4], double angle) {
    double matrix[][4] = { {cos(angle), -sin(angle), 0, 0}, {sin(angle), cos(angle), 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1} };
    
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            rot_z[i][j] = matrix[i][j];
        }
    }
} 

// Dx 
void translation_x(double Dx[4][4], double q) { // for the question q can either be displacement or offset
    double matrix[][4] = { {1, 0, 0, q}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1} };
  
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            Dx[i][j] = matrix[i][j];
        }
    }
}

// Dy
void translation_y(double Dy[4][4], double q) { 
    double matrix[][4]  = { {1, 0, 0, 0}, {0, 1, 0, q}, {0, 0, 1, 0}, {0, 0, 0, 1} };

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            Dy[i][j] = matrix[i][j];
        }
    }
}

// Dz 
void translation_z(double Dz[4][4], double q) {
    double matrix[][4] = { {1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, q}, {0, 0, 0, 1} };

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            Dz[i][j] = matrix[i][j];
        }
    }
}

// Matrix Multiplication Function
void matrix_multiplication(double product[4][4], double matrix_1[4][4], double matrix_2[4][4]){
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            product[i][j] = 0;
            for (int k = 0; k < 4; k++) {
                product[i][j] += matrix_1[i][k] * matrix_2[k][j];
            }
        } 
    } 
} 

 
fwd_kin(theta, x)
double *theta;
double x[3];
{
    double Rx_12[4][4];         // Rx(theta4) 
    rotation_x(Rx_12, theta[4]);
    double Dx_13[4][4];         // Dx(l3)
    translation_x(Dx_13, 0.15); 
    double product_1[4][4];    
    matrix_multiplication(product_1, Rx_12, Dx_13);
   
    double Dz_11[4][4];
    translation_z(Dz_11, -0.04);  // Dz(d4)
    double product_2[4][4];
    matrix_multiplication(product_2, Dz_11, product_1); 

    double Dy_10[4][4];
    translation_y(Dy_10, -0.04);    // Dy(d3)
    double product_3[4][4];
    matrix_multiplication(product_3, Dy_10, product_2);

    double Ry_9[4][4];
    rotation_y(Ry_9, theta[3]);   // Ry(theta3)
    double product_4[4][4];
    matrix_multiplication(product_4, Ry_9, product_3);

    double Dx_8[4][4];
    translation_x(Dx_8, 0.2);   // Dx(l2)
    double product_5[4][4];
    matrix_multiplication(product_5, Dx_8, product_4);

    double Dy_7[4][4];
    translation_y(Dy_7, 0.04);   // Dy(d2)
    double product_6[4][4];
    matrix_multiplication(product_6, Dy_7, product_5);

    double Ry_6[4][4];
    rotation_y(Ry_6, theta[2]);   // Ry(theta2)
    double product_7[4][4];
    matrix_multiplication(product_7, Ry_6, product_6);

    double Dx_5[4][4];
    translation_x(Dx_5, 0.2);    // Dx(l1) 
    double product_8[4][4];
    matrix_multiplication(product_8, Dx_5, product_7);
  
    double Ry_4[4][4];
    rotation_y(Ry_4, theta[1]);   // Ry(theta 4)
    double product_9[4][4];
    matrix_multiplication(product_9, Ry_4, product_8);

    double Dy_3[4][4];
    translation_y(Dy_3, -0.04);   // Dy(d1)
    double product_10[4][4];
    matrix_multiplication(product_10, Dy_3, product_9);

    double Dz_2[4][4];
    translation_z(Dz_2, 0.25);    // Dz(lo)
    double product_11[4][4];
    matrix_multiplication(product_11, Dz_2, product_10);

    double Rz_1[4][4];
    rotation_z(Rz_1, theta[0]);    //Rz(theta 0)
    double product_12[4][4];
    matrix_multiplication(product_12, Rz_1, product_11);

    x[0] = product_12[0][3];   // first row last column - x
    x[1] = product_12[1][3];   // second row last column - y
    x[2] = product_12[2][3];   // third row last column - z
}


inv_kin(x, theta)
double *x;
double theta[6];
{
    double gamma;
    double alpha;
    double gamma_1;          // to evaluate if it is +ve or -ve

    gamma = acos(x[0]/sqrt((x[0]*x[0]) + (x[1] * x[1])));
    gamma_1 = atan(x[1]/x[0]);
    alpha = asin((-0.04)/sqrt((x[0]*x[0])+(x[1]*x[1])));

    if (gamma_1 < 0){
        if (gamma > M_PI/2) {
	    theta[0] = (gamma - alpha);
        }
        else {
            theta[0] = -(gamma + alpha);
        }
    }
    else {
        if (gamma > M_PI/2) {
            theta[0] = -(gamma + alpha);
        }
        else {
            theta[0] = (gamma - alpha);
        }
    }
    theta[4] = 0;         // Theta 4 is always 0

    //  calculating (x', y' z')

    double l3 = 0.15;

//    double x_prime = x[0] + 0.04; 
//    double y_prime = x[1] + 0.04;
//    double z_prime = x[2] + l3;
   
    double d4 = -0.04;
    double d3 = -0.04;
    double l2 = 0.2;
    double l1 = 0.2;
    double l0 = 0.25;


    double x_prime = x[0] + ( - d4 * cos(theta[0])) + (d3 * sin(theta[0])); 
    double y_prime = x[1] + (- d4 * sin(theta[0])) - (d3 * cos(theta[0])); 
    double z_prime = x[2] + l3;

    double x_DoublePrime = sqrt((x_prime * x_prime) + (y_prime * y_prime)); 
    double y_DoublePrime = 0; 
    double z_DoublePrime = z_prime - l0; 

     
    double gamma_2 = atan(z_DoublePrime/x_DoublePrime);
    theta[2] = acos(((x_DoublePrime * x_DoublePrime) + (z_DoublePrime * z_DoublePrime) - (l1*l1) - (l2 * l2))/ (2 * l1 * l2)); 


    double alpha_2 = acos(((l1 * l1) + (x_DoublePrime * x_DoublePrime) + (z_DoublePrime * z_DoublePrime) - (l2 * l2)) / (2 * l1 * sqrt((x_DoublePrime * x_DoublePrime) + (z_DoublePrime * z_DoublePrime)))); 
    theta[1] = - alpha_2 - gamma_2;

    double phi = M_PI/2;
    theta[3] = phi - theta[1] - theta[2]; 
     
}
        









