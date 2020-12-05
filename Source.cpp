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
void recursivityDepthFirst(Graph* graph, Vertex* vertex, std::unordered_set<Vertex*>& visitedVertices, int& visited) {
    Vertex* toDo = graph->getFirstAdj(vertex);
    while (toDo != nullptr) {
        if (visitedVertices.find(toDo) == visitedVertices.end()) {
            visited++;
            visitedVertices.insert(toDo);
            recursivityDepthFirst(graph, toDo, visitedVertices, visited);
        }
        toDo = graph->getNextAdj(vertex, toDo);
    }
}

bool isRelatedDepthFirst(Graph* graph) {
    bool related = true;
    if (!graph->isEmpty()) {
        Vertex* vertex = graph->getHead();
        std::unordered_set<Vertex*> visitedVertices;
        int finalVisited = 1;
        visitedVertices.insert(vertex);
        recursivityDepthFirst(graph, vertex, visitedVertices, finalVisited);
        if (finalVisited != graph->getNumVertices()) {
            related = false;
        }
    }
    return related;
}

// Algoritmo D. Warshall.
void doMatrix(std::vector <std::vector<bool>>& matrix, int v) {
    int row = 0;
    int colum;
    while (row < matrix.size()) {
        if (row != v && matrix[row][v] == 1) {
            colum = 0;
            while (colum < matrix.size()) {
                if (matrix[row][colum] == 0 && matrix[v][colum] == 1) {
                    matrix[row][colum] = 1;
                }
                colum++;
            }
        }
        row++;
    }
}

bool isRelatedWarshall(Graph* graph) {
    bool related = true;
    if (!graph->isEmpty()) {
        int sizeGrafo = graph->getNumVertices();
        std::vector<std::vector<bool>> matrix;
        std::vector<bool> toPush(sizeGrafo, 0);

        for (int i = 0; i < sizeGrafo; i++) {
            matrix.push_back(toPush);
        }


        Vertex* vertex = graph->getHead();;
        Vertex* secondVertex;
        int row = 0;
        int colum;
        while (vertex != nullptr) {
            colum = 0;
            secondVertex = graph->getHead();
            while (secondVertex != nullptr) {
                if (graph->getWeight(vertex, secondVertex) > 0) {
                    matrix[row][colum] = 1;
                }
                else {
                    matrix[row][colum] = 0;
                }
                colum++;
                secondVertex = graph->getNextVert(secondVertex);
            }
            row++;
            vertex = graph->getNextVert(vertex);
        }

        int haciendoVertice = 0;

        while (haciendoVertice < sizeGrafo) {
            doMatrix(matrix, haciendoVertice);
            haciendoVertice++;
        }

        row = 0;
        while (row < sizeGrafo && related) {
            colum = 0;
            while (colum < sizeGrafo && related) {
                if (matrix[row][colum] == 0) {
                    related = false;
                }
                else {
                    colum++;
                }
            }
            row++;
        }
    }
    return related;
}

// Algoritmo E. Dijkstra
void dijkstra(Graph* graph, Vertex* vertex) {// a b d e f
    std::vector<double> weight_vector;//peso
    std::vector<Vertex*> path;//donde voy
    std::vector<Vertex*> pivots;//donde llego
    std::vector<bool> done;
    Vertex* temp_vertex;
    if (!graph->isEmpty()) {
        temp_vertex = graph->getHead();
        while (temp_vertex != nullptr) {
            if (temp_vertex != vertex) {
                done.push_back(0);
                path.push_back(vertex);
                pivots.push_back(temp_vertex);
                weight_vector.push_back(graph->getWeight(vertex, temp_vertex));
            }
            temp_vertex = graph->getNextVert(temp_vertex);
        }

        for (int i = 0; i < path.size(); i++) {
            std::cout << graph->getLabel(path[i]) << " to " << graph->getLabel(pivots[i]) << " peso:" << weight_vector[i] << std::endl;
        }

        int do_pivot = 0;
        int index;
        while (do_pivot < pivots.size() - 1) {
            index = 0;
            for (int i = 0; i < pivots.size(); i++) {
                if (done[index] && weight_vector[i] != -1) {
                    index = i;
                }
                if (weight_vector[index] > weight_vector[i] && done[i] != 1 && weight_vector[i] != -1) {
                    index = i;
                }
            }

            done[index] = true;
            temp_vertex = pivots[index];

            for (int i = 0; i < pivots.size(); i++) {
                if (index != i) {
                    if (weight_vector[i] > graph->getWeight(temp_vertex, pivots[i]) && !done[i] && weight_vector[i] != -1) {
                        weight_vector[i] = graph->getWeight(temp_vertex, pivots[i]);
                        path[i] = temp_vertex;
                    }
                    else if (weight_vector[i] == -1 && graph->getWeight(temp_vertex, pivots[i]) != -1) {
                        weight_vector[i] = graph->getWeight(temp_vertex, pivots[i]);
                        path[i] = temp_vertex;
                    }
                }
            }

            do_pivot++;
        }
        for (int i = 0; i < pivots.size(); i++) {
            std::cout << "Vertice [" << graph->getLabel(pivots[i]) << "]" << " = " << graph->getLabel(path[i]) << std::endl;
        }
    }
}

// Algoritmo F. Floyd.
void getMinimunCost(Graph* graph) {
    std::vector<double> weight_vector;
    std::vector<Vertex*> pivots;
    if (!graph->isEmpty()) {
        int num_vert = graph->getNumVertices(); //cantidad de nodos 
        Vertex* temp_vertex = graph->getHead();
        Vertex* help_vertex;

        while (temp_vertex != nullptr) {
            help_vertex = graph->getHead();
            while (help_vertex != nullptr) {
                weight_vector.push_back(graph->getWeight(temp_vertex, help_vertex));
                help_vertex = graph->getNextVert(help_vertex);
            }
            pivots.push_back(temp_vertex);
            temp_vertex = graph->getNextVert(temp_vertex);
        }

        std::cout << std::endl;
        double weight = 0.0;
        for (int pivot = 0; pivot < num_vert; pivot++) {
            for (int second_pivot = 0; second_pivot < num_vert; second_pivot++) {
                if (pivot != second_pivot) {
                    for (int third_pivot = 0; third_pivot < num_vert; third_pivot++) {
                        if (second_pivot != third_pivot && graph->isEdge(pivots[pivot], pivots[second_pivot]) && graph->isEdge(pivots[second_pivot], pivots[third_pivot])) {
                            weight = weight_vector[(pivot * num_vert) + second_pivot] + weight_vector[(second_pivot * num_vert) + third_pivot];
                            if (weight_vector[(pivot * num_vert) + third_pivot] > weight && weight > 0) {
                                std::cout << weight << " ";
                                weight_vector[(pivot * num_vert) + third_pivot] = weight;
                            }
                        }
                    }
                }

            }
        }
        std::cout << std::endl;
        for (int i = 0; i < weight_vector.size(); i++) {
            if (i % num_vert == 0 && i != 0) {
                std::cout << std::endl;
                std::cout << weight_vector[i] << " ";
            }
            else {
                std::cout << weight_vector[i] << "  ";
            }
        }
    }
}

// Algoritmo G. Hamilton.
bool find(std::vector<Vertex*> path, Vertex* v) {
    int row = 0;
    bool finded = false;
    while (row < path.size() && !finded) {
        if (path[row] == v) {
            finded = true;
        }
        row++;
    }
    return finded;
}
void doHamiltonCycle(Graph* graph, std::vector<Vertex*>& bestPath, std::vector<Vertex*>& path, int pos, int& bestWeight, int& pathWeight) {
    if (pos == path.size() && graph->isEdge(path[pos - 1], path[0])) {
        if (bestWeight > pathWeight || bestWeight == 0) {
            bestPath = path;
            bestWeight = pathWeight;
        }
    }
    else {
        Vertex* adyancent = graph->getFirstAdj(path[pos - 1]);
        while (adyancent != nullptr) {
            if (!find(path, adyancent)) {
                path[pos] = adyancent;
                pathWeight += graph->getWeight(path[pos - 1], path[pos]);
                doHamiltonCycle(graph, bestPath, path, pos + 1, bestWeight, pathWeight);
                pathWeight -= graph->getWeight(path[pos - 1], path[pos]);
                path[pos] = nullptr;
            }
            adyancent = graph->getNextAdj(path[pos - 1], adyancent);
        }
    }
}

void hamilton(Graph* graph) {
    if (!graph->isEmpty()) {
        std::vector<Vertex*> bestPath(graph->getNumVertices());
        std::vector<Vertex*> path(graph->getNumVertices());
        int bestWeight = 0;
        int pathWeight = 0;

        Vertex* vertex = graph->getHead();
        while (vertex != nullptr) {
            path[0] = vertex;
            doHamiltonCycle(graph, bestPath, path, 1, bestWeight, pathWeight);
            vertex = graph->getNextVert(vertex);
        }

        for (int i = 0; i < path.size(); i++) {
            if (i + 1 == path.size()) {
                std::cout << graph->getLabel(bestPath[i]) << " size:" << bestWeight << std::endl;
            }
            else {
                std::cout << graph->getLabel(bestPath[i]) << "->";
            }
        }
    }
}

// Algoritmo H.
void colorGraph(Graph* graph) {
    if (!graph->isEmpty()){
        int numV = graph->getNumVertices();
        bool* matrixAdj = new bool[numV*numV];

        Vertex* vertex = graph->getHead();
        Vertex* secondVertex;
        int row = 0;
        int colum;
        while (vertex != nullptr) {
            colum = 0;
            secondVertex = graph->getHead();
            while (secondVertex != nullptr) {
                if (graph->getWeight(vertex, secondVertex) > 0) {
                    matrixAdj[row*numV+colum] = 1;
                }
                else {
                    matrixAdj[row * numV + colum] = 0;
                }
                colum++;
                secondVertex = graph->getNextVert(secondVertex);
            }
            row++;
            vertex = graph->getNextVert(vertex);
        }

    }

}

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
void isolateVertex(Graph* graph, Vertex* vertex) {
    Vertex* adjVertex = graph->getFirstAdj(vertex);
    while (adjVertex != nullptr) {
        graph->deleteEdge(vertex, adjVertex);
        adjVertex = graph->getNextAdj(vertex, adjVertex);
    }
}

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
    // Grafo de prueba.
    Graph* test_graph = new Graph();    
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

    //std::cout << hasCycles(test_graph) << std::endl;

    //prim(test_graph);
    kruskal(test_graph);

    return 0;
}