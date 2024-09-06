//---------------------------------------------------------------------------

/***************************************************/
/* Functions prototypes by Prof. Anand Subramanian */
/***************************************************/

#ifndef Separation_H
#define Separation_H

#include <ilcplex/ilocplex.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <list>
#include <vector>
#include <algorithm>
#include <iterator>
#include <pthread.h>
#include <set>
#include <list>
#include <iterator>
#include <tuple>
#include "dataFunction.h"

// #define DEBUGGING

#define EPSILON 0.00001

using namespace std;

typedef struct{
    vector<int> nodes;
    int id;
    bool is_vertexA;
} vertex_type;

/* ---------------------------- helper functions ---------------------------- */
inline bool isg(double a, double b) { return a > b + EPSILON; }
inline bool isl(double a, double b) { return a < b - EPSILON; }
inline bool iseq(double a, double b) { return fabs(a - b) < EPSILON; }

/* --------------------------- max-back algorithm --------------------------- */
vector <vector<int> > MaxBack(double** x, int n);

/* ------------------- max-back algorithm helper functions ------------------ */
/**
 * @brief Calcula o valor de b(v) com relação a S.
 */
void init_b(double** x, int n, std::vector<double>& b, const std::set<int> S);

/**
 * @brief
 */
double compMinCut(double** x, int n, const std::set<int>& S);

/**
 * @brief Retorna o index do vértice com maior valor de max-back.
*/
int getMaxBack(const std::vector<double>& b);

/**
 * @brief Atualiza b(t) para todo t fora de S.
*/
void update_b(double** x, int n, std::vector<double>& b, 
              const std::set<int> S, int v);

vector <vector<int> > MaxBack(double** x, int n);

/* ------------------- min-cut algorithm helper functions ------------------- */
/**
 * @brief Calcula os pesos.
*/
void initWeights(double** x, 
                 int n, 
                 std::vector<double>& w,
                 const std::vector<std::vector<int>>& G, 
                 const std::set<int>& A);

/**
 * @brief Retorna o índice do vértice cujo arco tem maior peso.
 */
int getMaxWeight(const std::vector<double>& w);

/**
 * @brief Para todos os vértices de G ainda não adicionados a A, atualiza o peso
 * dos arcos conecatos ao nó de G recém adicionado a A.
 */
void updateWeights(double** x, 
                   const std::vector<std::vector<int>>& G,
                   int i,
                   const std::set<int>& A,
                   std::vector<double>& w);

/**
 * @brief Calcula o peso do corte que separa o vértice adicionado por último t 
 * do restante do grafo G.
 */
double compCutOfThePhase(double** x, 
                         const std::vector<std::vector<int>>& G,
                         int t);

/* ---------------------------- min-cut algorithm --------------------------- */               
/**
 * @brief Calcula o cut-of-the-phase e retorna iterators para os últimos 
 * vértices adicionados.
 */
std::tuple<int, int, double> MinCutPhase(double** x, int n, int a,
                                        const std::vector<std::vector<int>>& G);
                 
vector <vector<int> > MinCut(double** x, int n);

#endif

//---------------------------------------------------------------------------
