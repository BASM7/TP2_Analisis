#pragma once

#include <unordered_set>
#include <set>
#include <algorithm>
#include <stack>
#include <queue>
#include <map>
#include <iterator>
#include <vector>

#include "GraphB.h"

typedef std::unordered_set<Vertex*> dic_vertices;
typedef std::vector<std::pair<Vertex*, Vertex*>> set_of_edges;
typedef std::vector<double> vector_weights;

// Funciones auxiliares.
bool isVertexInDic(dic_vertices diccionary, Vertex* vertex) {
    return (diccionary.find(vertex) != diccionary.end());
}

int getIndex(set_of_edges edges, std::pair<Vertex*, Vertex*> edge){
    auto it = find(edges.begin(), edges.end(), edge);
    if (it != edges.end()){
        return it - edges.begin();
    }
    return -1;
}

void printGraphWidth(Graph* graph) {
    if (!graph->isEmpty()) {
        Vertex* currentVertex = graph->getHead();
        std::queue<Vertex*> queue;
        std::unordered_set<Vertex*> visitedVertices;
        std::unordered_set<Vertex*> alreadyPushedVertices;

        while (currentVertex != nullptr) {
            if (visitedVertices.find(currentVertex) == visitedVertices.end()) {
                visitedVertices.insert(currentVertex);
                queue.push(currentVertex);

                while (!queue.empty()) {
                    Vertex* tempVertex = queue.front();
                    queue.pop();

                    Vertex* adjVertex = graph->getFirstAdj(tempVertex);

                    while (adjVertex != nullptr) {
                        if (visitedVertices.find(adjVertex) == visitedVertices.end()) {
                            if (alreadyPushedVertices.find(adjVertex) == alreadyPushedVertices.end()) {
                                queue.push(adjVertex);
                                alreadyPushedVertices.insert(adjVertex);
                            }
                            std::cout << graph->getLabel(tempVertex) << " -> ";
                            std::cout << graph->getWeight(tempVertex, adjVertex);
                            std::cout << " -> " << graph->getLabel(adjVertex) << std::endl;
                        }


                        adjVertex = graph->getNextAdj(tempVertex, adjVertex);
                    }
                    visitedVertices.insert(tempVertex);
                }
            }

            currentVertex = graph->getNextVert(currentVertex);
        }
    }
}

// Algoritmos.

// Algoritmo A.
bool hasCycles(Graph* graph) {
    bool result = false;
    if (!graph->isEmpty()) {
        Vertex* currentVertex = graph->getHead();
        std::queue<Vertex*> queue;
        std::unordered_set<Vertex*> visitedVertices;
        std::unordered_set<Vertex*> alreadyPushedVertices;

        while (currentVertex != nullptr) {
            if (visitedVertices.find(currentVertex) == visitedVertices.end()) {
                visitedVertices.insert(currentVertex);
                queue.push(currentVertex);

                while (!queue.empty()) {
                    Vertex* tempVertex = queue.front();
                    queue.pop();

                    Vertex* adjVertex = graph->getFirstAdj(tempVertex);

                    while (adjVertex != nullptr) {
                        if (visitedVertices.find(adjVertex) == visitedVertices.end()) {
                            if (alreadyPushedVertices.find(adjVertex) == alreadyPushedVertices.end()) {
                                queue.push(adjVertex);
                                alreadyPushedVertices.insert(adjVertex);
                            }
                            else {
                                result = true;
                            }
                        }
                        adjVertex = graph->getNextAdj(tempVertex, adjVertex);
                    }
                    visitedVertices.insert(tempVertex);
                }
            }
            currentVertex = graph->getNextVert(currentVertex);
        }
    }
    return result;
}


// Algoritmo B.


// Algoritmo C.


// Algoritmo D.


// Algoritmo E.


// Algoritmo F.


// Algoritmo G.


// Algoritmo H.


// Algoritmo I.
Vertex* getMinEdgesPrim(dic_vertices pivotsUsed, set_of_edges edges, vector_weights weights) {
    Vertex* newPivot = nullptr;
    double min_weight = -1.0;
    for (int i = 0; i < edges.size(); i++) {
        Vertex* tempPivot = edges[i].second;
        double tempWeight = weights[i];
        if (!isVertexInDic(pivotsUsed, tempPivot)) {
            if (min_weight == -1.0 || tempWeight < min_weight) {
                min_weight = tempWeight;
                newPivot = tempPivot;
            }
        }
    }
    return newPivot;
}

void prim(Graph* graph) {
    if (!graph->isEmpty()) {
        Vertex* basePivot = graph->getHead();

        dic_vertices pivotsUsed;
        set_of_edges edges;
        vector_weights weights;

        // Inicialización.
        Vertex* tempVertex = graph->getNextVert(basePivot);
        while (tempVertex != nullptr) {
            edges.push_back(std::make_pair(basePivot, tempVertex));
            weights.push_back(graph->getWeight(basePivot, tempVertex));
            tempVertex = graph->getNextVert(tempVertex);
        }

        while (pivotsUsed.size() < graph->getNumVertices()) {
            basePivot = getMinEdgesPrim(pivotsUsed, edges, weights);
            for (int i = 0; i < edges.size(); i++) {
                Vertex* tempVertex = edges[i].second;
                if (!isVertexInDic(pivotsUsed, tempVertex)) {
                    double tempWeight = weights[i];
                    double edgeWeight = graph->getWeight(basePivot, tempVertex);
                    if (tempWeight == -1.0 || (tempWeight > edgeWeight && edgeWeight != -1.0)) {
                        edges[i].first = basePivot;
                        weights[i] = edgeWeight;
                    }
                }
            }
            pivotsUsed.insert(basePivot);
        }

        std::cout << "Arbol de minimo costo: " << std::endl;
        for (int i = 0; i < edges.size(); i++) {
            std::cout
                << edges[i].first->getValue() << " -> "
                << edges[i].second->getValue() << " con peso: "
                << weights[i] << std::endl;
        }
    }
}

// Algoritmo J.

int aux_kruskalGetSetIndex(std::vector<std::set<Vertex*>> setOfSet, Vertex* vertex) {    
    bool found = false;
    int i = 0;
    int setIndex1 = i;
    while (i < setOfSet.size() && !found) {
        std::set<Vertex*> tempSet = setOfSet[i];
        if (tempSet.find(vertex) != tempSet.end()) {
            setIndex1 = i;
            found = true;
        }            
        i++;
    }
    return setIndex1;
}

void kruskal(Graph* graph) {
    // 1. Meter todas las aristas.
    // 2. Recorrer aristas, seleccionar la arista con menor peso.
    // 3. Chequear si al agregar esa arista se forma un ciclo. Si hay, rechaza arista.
    // 4. Si no hay ciclo, añade arista.
    // 5. Cuando cant aristas añadidas = nunVertices-1 termina, sino continua.

    if (!graph->isEmpty()) {

        std::priority_queue<std::pair<double, std::pair<Vertex*, Vertex*>>> edges;
        std::set<std::pair<Vertex*, Vertex*>> diccionaryEdges;
        std::vector<std::set<Vertex*>> setOfSet;
        std::vector<std::pair<Vertex*, Vertex*>> outputEdges;

        // Paso 1. Inicialización.
        Vertex* tempVertex = graph->getHead();
        while( tempVertex != nullptr ){
            //  Al conjunto de conjuntos le asigna un unico elemento.
            std::set<Vertex*> newSet = { tempVertex };
            setOfSet.push_back(newSet);
            Vertex* adjVertex = graph->getFirstAdj(tempVertex);

            while (adjVertex != nullptr) {

                std::pair<Vertex*, Vertex*> newEdgeToFrom(tempVertex, adjVertex);
                std::pair<Vertex*, Vertex*> newEdgeFromTo(adjVertex, tempVertex);
                if (diccionaryEdges.find(newEdgeToFrom) == diccionaryEdges.end() &&
                    diccionaryEdges.find(newEdgeFromTo) == diccionaryEdges.end()) {
                    diccionaryEdges.insert(newEdgeToFrom);

                    // Guarda la arista en la cola de prioridad.
                    // El peso es multiplicado por -1, para que que la salida sea de menor a mayor.
                    std::pair<double, std::pair<Vertex*, Vertex*>> edge(-1 * graph->getWeight(tempVertex, adjVertex), newEdgeToFrom);
                    edges.push(edge);
                }
                adjVertex = graph->getNextAdj(tempVertex, adjVertex);
            }
            tempVertex = graph->getNextVert(tempVertex);
        }
        
        // Paso 2. Escoger las aristas.
        int chosenEdges = 0;
        while (chosenEdges < graph->getNumVertices() - 1 && !edges.empty()) {

            // Sacamos la arista de la cola.
            std::pair<Vertex*, Vertex*> edge = edges.top().second;
            edges.pop();

            int index_set1 = aux_kruskalGetSetIndex(setOfSet, edge.first);
            int index_set2 = aux_kruskalGetSetIndex(setOfSet, edge.second);

            if (index_set1 != index_set2) {
                outputEdges.push_back(edge);

                std::set<Vertex*> set1 = setOfSet[index_set1];
                std::set<Vertex*> set2 = setOfSet[index_set2];

                std::set<Vertex*> newSet;
                std::set_union(set1.begin(), set1.end(), set2.begin(), set2.end(), std::inserter(newSet, newSet.begin()));

                std::vector<std::set<Vertex*>>::iterator iter;
                iter = std::find(setOfSet.begin(), setOfSet.end(), set1);
                setOfSet.erase(iter);
                iter = std::find(setOfSet.begin(), setOfSet.end(), set2);
                setOfSet.erase(iter);

                setOfSet.push_back(newSet);

                chosenEdges++;
            }
        }

        std::cout << "Arbol de minimo costo: " << std::endl;
        for (int i = 0; i < outputEdges.size(); i++) {
            std::cout
                << outputEdges[i].first->getValue() << " -> "
                << outputEdges[i].second->getValue() << " con peso: "
                << graph->getWeight(outputEdges[i].first, outputEdges[i].second) << std::endl;
        }
    }
}

// Algoritmo K.


// Algoritmo L.
double totalWeight(Graph* graph) {
    double total = 0.0;
    if (!graph->isEmpty()) {
        Vertex* currentVertex = graph->getHead();
        std::queue<Vertex*> queue;
        std::unordered_set<Vertex*> visitedVertices;
        std::unordered_set<Vertex*> alreadyPushedVertices;

        while (currentVertex != nullptr) {
            if (visitedVertices.find(currentVertex) == visitedVertices.end()) {
                visitedVertices.insert(currentVertex);
                queue.push(currentVertex);

                while (!queue.empty()) {
                    Vertex* tempVertex = queue.front();
                    queue.pop();

                    Vertex* adjVertex = graph->getFirstAdj(tempVertex);

                    while (adjVertex != nullptr) {
                        if (visitedVertices.find(adjVertex) == visitedVertices.end()) {
                            if (alreadyPushedVertices.find(adjVertex) == alreadyPushedVertices.end()) {
                                queue.push(adjVertex);
                                alreadyPushedVertices.insert(adjVertex);
                            }
                            total += graph->getWeight(tempVertex, adjVertex);
                        }
                        adjVertex = graph->getNextAdj(tempVertex, adjVertex);

                    }
                    visitedVertices.insert(tempVertex);
                }
            }
            currentVertex = graph->getNextVert(currentVertex);
        }
    }
    return total;
}

// Algoritmo M.
Vertex* getVertex(Graph* graph, char labelToSearch) {
    Vertex* foundVertex = graph->getHead();
    bool found = false;
    while (foundVertex != nullptr && !found) {
        if (graph->getLabel(foundVertex) == labelToSearch) {
            return foundVertex;
        }
        foundVertex = graph->getNextVert(foundVertex);
    }
    return nullptr;
}



int main()
{
    Graph* test_graph = new Graph();
    //Test case 1.
    Vertex* vA = test_graph->addVert('A');
    Vertex* vB = test_graph->addVert('B');
    Vertex* vC = test_graph->addVert('C');
    Vertex* vD = test_graph->addVert('D');
    Vertex* vE = test_graph->addVert('E');
    Vertex* vF = test_graph->addVert('F');

    test_graph->createEdge(vA, vB, 2.0);
    test_graph->createEdge(vA, vC, 8.0);
    test_graph->createEdge(vA, vD, 6.0);
    test_graph->createEdge(vC, vB, 3.0);
    test_graph->createEdge(vA, vF, 3.0);
    test_graph->createEdge(vA, vE, 7.0);
    test_graph->createEdge(vB, vD, 9.0);
    test_graph->createEdge(vB, vF, 5.0);
    test_graph->createEdge(vC, vF, 6.0);
    test_graph->createEdge(vC, vE, 1.0);
    test_graph->createEdge(vE, vF, 4.0);
    test_graph->createEdge(vF, vD, 9.0);

    //printGraphWidth(test_graph);

    std::cout << hasCycles(test_graph) << std::endl;

    //prim(test_graph);
    //kruskal(test_graph);

    return 0;
}