#include "separation.h"

/* ------------------- max-back algorithm helper functions ------------------ */
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

double compMinCut(double** x, int n, const std::set<int>& S) {
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

void update_b(double** x, int n, std::vector<double>& b, 
              const std::set<int> S, int v) {
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

/* ------------------- min-cut algorithm helper functions ------------------- */
void initWeights(double** x, 
                 int n, 
                 std::vector<double>& w,
                 const std::vector<std::vector<int>>& G, 
                 const std::set<int>& A) {
    // iterate over the nodes in G
    for(auto i = 0; i < G.size(); ++i) {
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

int getMaxWeight(const std::vector<double>& w) {
    double max_w = w[0];
    int max_i = 0;
    for(int i = 1; i < w.size(); ++i) {
        if(isg(w[i], max_w)) {
            max_w = w[i];
            max_i = i;
        }
    }
    return max_i;
}

void updateWeights(double** x, 
                   const std::vector<std::vector<int>>& G,
                   int i,
                   const std::set<int>& A,
                   std::vector<double>& w) {
    for(int j = 0; j < G.size(); ++j) {
        if(j == i) {
            continue;
        }
        for(int v : G[j]) {
            // if one element of G[j] is in A, then all elements of G[j] are
            if(A.find(v) != A.end()) {
                break;
            }
            for(int u : G[i]) {
                if(v > u) {
                    w[j] += x[u][v];
                } else if(u > v) {
                    w[j] += x[v][u];
                }
            }
        }
    }
}

double compCutOfThePhase(double** x,
                         const std::vector<std::vector<int>>& G,
                         int t) {
    double cut = 0.0;
    for(int v : G[t]) {
        for(int i = 0; i < G.size(); ++i) {
            if(i == t) {
                continue;
            }
            for(int u : G[i]) {
                if(v > u) {
                    cut += x[u][v];
                } else if(u > v) {
                    cut += x[v][u];
                }
            }
        }
    }
    return cut;
}

/* --------------------------- max-back algorithm --------------------------- */
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
        double cutmin = compMinCut(x, n, S);
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

/* ---------------------------- min-cut algorithm --------------------------- */
std::tuple<int, int, double> MinCutPhase(double** x, int n, int a,
                                       const std::vector<std::vector<int>>& G) {
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
    #ifdef DEBUGGING
        auto printSetLambda = [G](int i, char label) {
            std::cout << "{";
            for(auto g : G[i]) {
                std::cout << g + 1 << " ";
            }
            std::cout << "}:" << label << " ";
        };
        char k = 'a' + 1;
    #endif

    std::set<int> A;
    A.insert(a);
    std::vector<double> weights(G.size(), 0.0);
    initWeights(x, n, weights, G, A);
    int s = 0;
    int t = 0;
    while(A.size() < n) {
        int i = getMaxWeight(weights);
        #ifdef DEBUGGING
            std::cout << "\nweights: ";
            for(int j = 0; j < G.size(); ++j) {
                if(weights[j] > 0.5) {
                    printSetLambda(j, int(weights[j]) + '0');
                }
            }
            std::cout << "selected:";
            printSetLambda(i, k);
            ++k;
        #endif
        A.insert(G[i].begin(), G[i].end());
        weights[i] = -1.0;
        updateWeights(x, G, i, A, weights);
        // update last vertices s and t inserted in A
        s = t;
        t = i;
    }
    #ifdef DEBUGGING
        std::cout << std::endl;
        printSetLambda(s, 's');
        printSetLambda(t, 't');
    #endif

    double w = compCutOfThePhase(x, G, t);

    return std::make_tuple(s, t, w);
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
    std::vector<std::vector<int>> G;
    for(int v = 0; v < n; ++v) {
        G.push_back(std::vector<int>({v}));
    }

    // min cut algorithm
    double mincut = INFINITE;
    std::vector<int> Smin;
    while(G.size() > 1) {
        #ifdef DEBUGGING
            std::cout << "---\nG:" << std::endl;
            for(int i = 0; i < G.size(); ++i) {
                std::cout << "node " << i + 1 << ":";
                for(int j = 0; j < G[i].size(); ++j) {
                    std::cout << G[i][j] + 1 << " ";
                }
                std::cout << std::endl;
            }
        #endif

        auto [s, t, w] = MinCutPhase(x, n, 1, G);

        #ifdef DEBUGGING
            std::cout << "\nCut-of-the-phase: " << w << std::endl;
            getchar();
        #endif

        // if a new min cut is found, update mincut and Smin
        if(isl(w, mincut)) {
            mincut = w;
            Smin = std::vector<int>(G[t].begin(), G[t].end());
        }

        // merge the last two nodes added to A
        // the node with less vertices is merged into the node with more 
        // vertices 
        if(G[s].size() < G[t].size()) {
            G[t].insert(G[t].end(), G[s].begin(), G[s].end());
            G.erase(G.begin() + s);
        } else {
            G[s].insert(G[s].end(), G[t].begin(), G[t].end());
            G.erase(G.begin() + t);
        }
    }
    #ifdef DEBUGGING
        std::cout << "\nMin-cut value:" << mincut 
                  << " size:" << Smin.size() 
                  << " cut:{";
        for(int i = 0; i < Smin.size(); ++i) {
            std::cout << Smin[i] + 1 << " ";
        }
        std::cout << "}\n" << std::endl;
    #endif

    if(Smin.size() > 2) {
        return std::vector<std::vector<int>>({Smin});
    } else {
        return std::vector<std::vector<int>>();
    }
}