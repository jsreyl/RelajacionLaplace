//Quiero resolver ecución de Laplace 2d: d2U/dx2 + d2U/dy2 = 0

//HACER CON UNA FUENTE A.K.A. POISSON
//HACER CON UNA BARRA A POTENCIAL = 50V EN LA MITAD DEL PLANO
#include<iostream>
#include<vector>

const double Delta = 0.1; //tamaño de la variación
const double L = 1.479; //tamaño del dominio
const int N = int(L/Delta)+1; //número de divisiones de la malla
const int STEPS = 100;

typedef std::vector<double> Matrix;

void initial_conditions(Matrix & m);
void boundary_conditions(Matrix & m);
void evolve(Matrix & m);
void print(const Matrix & m);

int main(int argc, char **argv){
  double xmax;
  double ymax;
  
  Matrix data(N*N);

  //Condiciones iniciales: llenamos la matriz de ceros
  initial_conditions(data);

  //Condiciones de frontera: llenamos solo los datos que necesitan condición
  boundary_conditions(data);
  
  //Hacemos loop de la función de evolución
  for (int ii=0; ii<=STEPS; ++ii){
    evolve(data);
  }
  print(data); 
  return 0;
}

void initial_conditions(Matrix & m){
  for (int ii= 0; ii<N; ++ii)
    for(int jj=0; jj<N; ++jj)
      m[ii*N+jj] = 1.0;
}

void print(const Matrix & m){
  for (int ii= 0; ii<N; ++ii){
    for(int jj=0; jj<N; ++jj){
      //std::cout<<m[ii*N+jj]<<" "; //Para imprimir la matriz
      std::cout<<ii*Delta<<" "<<jj*Delta<<" "<<m[ii*N+jj]<<std::endl;
    }
    std::cout<<"\n";
  }

}

void boundary_conditions(Matrix & m){
  int ii=0, jj=0;  
  //A la pared de arriba le asigno 100
  ii = 0;
  for(int jj=0; jj<N; ++jj){
    m[ii*N+jj]= 100.0;
  }
  //A la pared de abajo le asigno 0
  ii = N-1;
  for(int jj=0; jj<N; ++jj){
    m[ii*N+jj]= 0.0;
  }
  //A la pared de la izquierda le asigno 0
  jj = 0;
  for(int ii = 1; ii<N-1; ++ii){ //Puesto que ya llené ii=0 e ii= N
    m[ii*N+jj]= 0.0;
  }
  //A la pared de la derecha le asigno 0
  jj = N-1;
  for(int ii = 1; ii<N-1; ++ii){ //Puesto que ya llené ii=0 e ii= N
    m[ii*N+jj]= 0.0;
  }
  
}

void evolve(Matrix & m){
  for(int ii=0; ii<N; ++ii){
    for(int jj = 0; jj<N; ++jj){
      if(ii == 0) continue;
      if(ii == N-1) continue;
      if(jj == 0) continue;
      if(jj == N-1) continue;
      m[ii*N+jj]=(m[(ii+1)*N+jj] +
		  m[(ii-1)*N+jj] +
		  m[ii*N+(jj+1)] +
		  m[ii*N + (jj-1)])/4.0;
	}
  }
}
