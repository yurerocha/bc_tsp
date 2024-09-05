#include "separation.h"

bool isg(double a, double b) {
    return a > b + EPSILON;
}

bool isl(double a, double b) {
    return a < b - EPSILON;
}

bool iseq(double a, double b) {
    return fabs(a - b) < EPSILON;
}

/**
 * @brief Calcula o valor de b(v) com relação a S.
*/
void init_b(double** x, int n, std::vector<double>& b, const std::set<int> S) {
    for(int v = 0; v < b.size(); ++v) {
        // if v is not in set S
        if(S.find(v) == S.end()) {
            for(int u : S) {
                if(v > u) {
                    b[v] += x[u][v];
                }
            }
        }
    }
}

double compCutMin(double** x, int n, const std::set<int>& S) {
    double cutval = 0.0;
    for(int u = 0; u < n; ++u) {
        // if u is a node outside S
        if(S.find(u) == S.end()) {
            for(int v : S) {
                if(v > u) {
                    cutval += x[u][v];
                }
            }
        }
    }
    return cutval;
}

/**
 * @brief Retorna o index do vértice com maior valor de max-back.
*/
int getMaxBack(const std::vector<double>& b) {
    double max_b = b[0];
    int v = 0;
    for(int i = 1; i < b.size(); ++i) {
        if(isg(b[i], max_b)) {
            max_b = b[i];
            v = i;
        }
    }
    return v;
}

/**
 * @brief Atualiza b(t) para todo t fora de S.
*/
void update_b(double** x, int n, std::vector<double>& b, const std::set<int> S, int v) {
    for(int t = 0; t < n; ++t) {
        if(S.find(t) == S.end()) {
            if(t > v) {
                b[t] += x[v][t];
            } else {
                b[t] += x[t][v];
            }
        }
    }
}

vector <vector<int> > MaxBack(double** x, int n) {
    /*
        Notação:
            -delta(S) representa o conjunto de nós com exatamente 1 endnode em S
            -escrevemos delta(v) em vez de delta({v}) para v in V
            -a função de peso x* é 1 p arcos utilizados e 0 c.c.?
            -e representa uma aresta no grafo
            -u, v representam vértices no grafo
    */
    /*
        Começar com um S0 arbitrário ex: S0 = {0} ou {}
        Cutmin = somatório dos pesos de todas as arestas com um endnode em S0
            será necessário passar por todas as linhas dos elementos em S0
        Para todo v fora de S0, calcular b(v) com relação a S0, ou seja, 
        o valor de max-back de v com realação a S0
            também será necessário passar por todos os elementos de S0 e de 
            V - S0
        S = Smin = S0, Cutval = Cutmin
        Enquanto S != V, selecionar o v fora de S com maior valor de max-back
        Atualizar:
            S = S U {v}
                essa operação pode ser mais custosa
            Cutval = Cutval + 2 - 2 * b(v)
            Considerando que v já está em S:
                Para todo t fora de S, atualizar b(t) para b(t) + x*vt 
                (b(t) = b(t) + x*vt)
            Se Cutval < Cutmin
                Cutmin = Cutval
                Smin = S
        
        Ao final do algoritmo, os vértices em Smin não contêm subtours
    */
    std::vector<std::vector<int>> subtours;
    // run the max-back algorithm for multiple starting nodes
    // due to symmetry, we only need to run it for n / 2 starting nodes
    for(int v0 = 0; v0 < n / 2; ++v0) {
        std::set<int> S({v0});
        auto Smin = S;
        auto b = std::vector<double>(n, 0.0);
        double cutmin = compCutMin(x, n, S);
        double cutval = cutmin;
        init_b(x, n, b, S);
        b[0] = -1;
        #ifdef DEBUGGING
            std::cout << "\n\nmax-back" << std::endl;
            std::cout << "sol:" << std::endl;
            for(int i = 0; i < n; ++i) {
                for(int j = i + 1; j < n; ++j) {
                    if(x[i][j] > EPSILON) {
                        std::cout << i << " " << j << ":" << x[i][j] << std::endl;
                    }
                }
            }
            std::cout << std::endl;
        #endif
        while(S.size() < n) {
            // choose v in S of maximum max-back value
            #ifdef DEBUGGING
                std:cout << "b" << std::endl;
                for(int i = 0; i < b.size(); ++i) {
                    std::cout << b[i] << " ";
                }
                std::cout << std::endl;
            #endif
            int v = getMaxBack(b);
            // S = S U {v}
            S.insert(v);
            // cutval = cutval + 2 - 2 * b(v)
            cutval = cutval + 2.0 - 2.0 * b[v];
            b[v] = -1.0;
            #ifdef DEBUGGING
                std::cout << S.size() << " " << v << " " << cutval << " " << cutmin << std::endl;
                getchar();
            #endif
            // update b(t) for all t outside S
            update_b(x, n, b, S, v);

            if(isl(cutval, cutmin)) {
                cutmin = cutval;
                Smin = S;
                #ifdef DEBUGGING
                    std::cout << "cut found" << std::endl;
                    std::vector<int> tmp(Smin.begin(), Smin.end());
                    for(int i = 0; i < tmp.size(); ++i) {
                        cout << tmp[i] << " ";
                    }
                    std::cout << std::endl;
                #endif
            }
        }
        // if there is more than one node in Smin, then it is a subtour
        // also, discard the tour composed of all vertices (i.e., a feas sol)
        if(Smin.size() > 1 && Smin.size() < n) {
            subtours.push_back(std::vector<int>(Smin.begin(), Smin.end()));
        #ifdef DEBUGGING
            std::vector<int> SminVec(Smin.begin(), Smin.end());
            std::cout << "Cutmin: " << cutmin << std::endl;
            for(int i = 0; i < SminVec.size(); ++i) {
                cout << SminVec[i] << " ";
            }
        #endif
        }
    }

   // retornar o Smin como um Vec2D. No MaxCut, só seria necessário retornar
   // um vector com os arcos que devem formar a restrição de corte
   // palpite: o MinCut retornará multiplos vectors, i.e., cortes
   return subtours;
}

vector <vector<int> > MinCut(double** x, int n) {
    /*
    MinCut:
        Enquanto |V| > 1
            MinCutPhase
            Se o novo cut-of-the-phase for menor do que o atual, atualizar o
            cut-of-the-phase
        Retornar o menor cut-of-the-phase
    */
    // initialize G, with one vertex for each node  
    std::list<std::list<int>> G;
    for(int v = 0; v < n; ++v) {
        G.push_back(std::list<int>({v}));
    }

    // min cut algorithm
    double mincut = INFINITE;
    std::vector<int> Smin;
    while(G.size() > 1) {
        double cut = 0.0;
        auto [it_s, it_t] = MinCutPhase(x, n, 0, G, cut);
        // merge the last two nodes added to A
        auto it = it_t;
        if(it_s->size() < it_t->size()) {
            it_t->insert(it_t->end(), it_s->begin(), it_s->end());
            G.erase(it_s);
        } else {
            it = it_s;
            it_s->insert(it_s->end(), it_t->begin(), it_t->end());
            G.erase(it_t);
        }
        if(isl(cut, mincut)) {
            mincut = cut;
            Smin = std::vector<int>(it->begin(), it->end());
        }
    }
}

/**
 * @brief Calcula os pesos.
*/

void initWeights(double** x, 
                 int n, 
                 std::vector<double>& w,
                 const std::list<std::list<int>>& G, 
                 const std::set<int>& A) {
    // iterate over the nodes in G
    for(int i = 0; i < G.size(); ++i) {
        for(int v : G[i]) {
            // if one element of G[i] is in A, then all elements of G[i] are
            if(A.find(v) != A.end()) {
                break;
            }

            for(int u : A) {
                if(v > u) {
                    w[i] += x[u][v];
                } else if(u > v) {
                    w[i] += x[v][u];
                }
            }
        }
    }
}

/**
 * @brief Retorna o índice do vértice cujo arco tem maior peso.
*/
auto getMaxWeight(const std::vector<double>& w, 
                  const std::list<std::list<int>>& G) {
    double max_w = w[0];
    int max_i = 0;
    auto it_i = G.begin();
    auto it_max_i = it_i;
    for(int i = 1; i < w.size(); ++i) {
        if(isg(w[i], max_w)) {
            max_w = w[i];
            max_i = i;
            it_max_i = it_i;
        }
        ++it_i;
    }
    return std::make_pair(max_i, it_max_i);
}

/**
 * @brief Atualiza os pesos dos arcos conectados aos vértices do nó G recém
 * adicionado a A.
*/
void updateWeights(double** x, 
                   int n, 
                   const std::list<std::list<int>>& G,
                   std::iterator<std::list<std::list<int>>> it_i,
                   const std::set<int>& A,
                   std::vector<double>& w) {
    for(auto it_j = G.begin(); it_j != G.end(); ++it_j) {
        if(it_j == it_i) {
            continue;
        })
        for(int v : *it_j) {
            // if one element of G[j] is in A, then all elements of G[j] are
            if(A.find(v) != A.end()) {
                break;
            }
            for(int u : *it_i) {
                if(v > u) {
                    w[j] += x[u][v];
                } else if(u > v) {
                    w[j] += x[v][u];
                }
            }
        }
    }
}

/**
 * @brief Calcula o cut-of-the-phase e retorna os últimos vértices adicionados.
*/
auto MinCutPhase(double** x, 
                 int n, 
                 int a,
                 const std::list<std::list<int>>& G, 
                 double& cut) {
    /*
    MinCutPhase:
        A = {a}
        Para todo v em V, calcular a soma dos pesos de todos os arcos com um 
        endnode em A
        Enquanto A != V, adicionar ao grafo A o vértice com maior valor da soma 
        calculada no passo anterior
        Ao finalizar o loop anterior, armazenar o "cut-of-the-phase", i.e., a 
        soma dos pesos dos arcos que conectam o último vértice t adicionado a A
        Fazer o merge, no grafo original, dos últimos nós s e t adicionados a A
            -a aresta conectando s e t "some"
            -o peso da aresta desse novo nó é igual à soma dos pesos
            das arestas que conectavam s e t ao grafo A
        Retornar o último cut-of-the-phase
    */
    std::set<int> A;
    A.insert(a);
    std::vector<double> weights(G.size(), 0.0);
    initWeights(x, n, weights, G, A);
    auto it_s = G.begin();
    auto it_t = G.begin();
    while(A.size() < n) {
        auto [i, it] = getMaxWeight(weights);
        // store the cut-of-the-phase
        cut = weights[i];
        A.merge(*it);
        weights[i] = -1.0;
        updateWeights(x, n, G, it, A, weights);
        // update last vertices s and t inserted in A
        it_s = it_t;
        it_t = it;
    }

    return std::make_pair(it_s, it_t);
}
