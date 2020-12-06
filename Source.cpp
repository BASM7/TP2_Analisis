#pragma once

#include <windows.h>
#include <unordered_set>
#include <set>
#include <algorithm>
#include <stack>
#include <queue>
#include <map>
#include <iterator>
#include <vector>
#include <locale>
#include <conio.h>
#include <iomanip>

//#include "GraphA.h"
#include "GraphB.h"

typedef std::unordered_set<Vertex*> dic_vertices;
typedef std::vector<std::pair<Vertex*, Vertex*>> set_of_edges;
typedef std::vector<double> vector_weights;

std::map<std::string, Graph*> graphRegistry;

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
bool hasCycleWidth(Graph* graph) {
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
void getCicle(Graph* graph, std::unordered_set<Vertex*> visitedVertices, Vertex* vertex, Vertex* adyacent, bool& cicle) {
    Vertex* adyacentAdyacent = graph->getFirstAdj(adyacent);
    while (adyacentAdyacent != nullptr && !cicle) {
        if (visitedVertices.find(adyacentAdyacent) != visitedVertices.end() && graph->isEdge(vertex, adyacentAdyacent) && adyacentAdyacent != vertex) {
            cicle = true;
        }
        else {
            adyacentAdyacent = graph->getNextAdj(adyacent, adyacentAdyacent);
        }
    }
}

void depthRecur(Graph* graph, Vertex* vertex, std::unordered_set<Vertex*>& visitedVertices, bool& cicle) { // a encola b y b encola a por eso 

    visitedVertices.insert(vertex);
    Vertex* adyacent = graph->getFirstAdj(vertex);
    while (adyacent != nullptr && !cicle) {
        if (visitedVertices.find(adyacent) == visitedVertices.end()) {
            getCicle(graph, visitedVertices, vertex, adyacent, cicle);
            depthRecur(graph, adyacent, visitedVertices, cicle);
        }
        adyacent = graph->getNextAdj(vertex, adyacent);
    }
}

bool cicleGraphDepth(Graph* graph) {
    bool cicle = false;
    if (!graph->isEmpty()) {
        Vertex* vertex = graph->getHead();
        std::unordered_set<Vertex*> visitedVertices;
        while (vertex != nullptr && !cicle) {
            if (visitedVertices.find(vertex) == visitedVertices.end()) {
                depthRecur(graph, vertex, visitedVertices, cicle);
            }
            vertex = graph->getNextVert(vertex);
        }
    }
    return cicle;
}

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
        done[0] = true;
        int do_pivot = 0;
        int index_less;
        while (do_pivot < pivots.size() - 1) {
            index_less = 0;
            for (int i = 0; i < pivots.size(); i++) {
                if (weight_vector[i] != -1 && weight_vector[index_less] == -1 || weight_vector[index_less] > weight_vector[i]) {
                    index_less = i;
                }
            }

            done[index_less] = true;
            temp_vertex = pivots[index_less];

            for (int i = 0; i < pivots.size(); i++) {
                if (index_less != i) {
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
            std::cout << "Ir " << graph->getLabel(pivots[i])
                << " por " << graph->getLabel(path[i]) << std::endl;
        }
    }
}

// Algoritmo F. Floyd.
void floyd(Graph* graph) {
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

        double weight = 0.0;
        for (int pivot = 0; pivot < num_vert; pivot++) {
            for (int second_pivot = 0; second_pivot < num_vert; second_pivot++) {
                if (pivot != second_pivot && graph->isEdge(pivots[pivot], pivots[second_pivot])) {
                    for (int third_pivot = 0; third_pivot < num_vert; third_pivot++) {
                        if (second_pivot != third_pivot && pivot != third_pivot && graph->isEdge(pivots[second_pivot], pivots[third_pivot])) {
                            weight = weight_vector[(pivot * num_vert) + second_pivot] + weight_vector[(second_pivot * num_vert) + third_pivot];
                            if (weight_vector[(pivot * num_vert) + third_pivot] > weight && weight > 0 || weight_vector[(pivot * num_vert) + third_pivot] == -1 && weight != -1) {
                                weight_vector[(pivot * num_vert) + third_pivot] = weight;
                            }
                        }
                    }
                }

            }
        }
        for (int i = 0; i < weight_vector.size(); i++) {
            if (i % num_vert == 0 && i != 0) {
                std::cout << std::endl;
                std::cout << std::setw(4) << weight_vector[i] << " ";
            }
            else {
                std::cout << std::setw(4)<< weight_vector[i] << "  ";
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
bool bestSolutionColoring(Graph* graph, bool* matrixAdj, int numPosibleColors, int indexV, int* &colors) {
    if (indexV == graph->getNumVertices()) {
        return true;
    }

    for (int i = 0; i < numPosibleColors + 1; i++) {
        bool isSafeToPaint = true;
        int j = 0;
        while (j < graph->getNumVertices() && isSafeToPaint) {
            if (matrixAdj[indexV*graph->getNumVertices()+j] == 1 && colors[j] == i) {
                isSafeToPaint = false;
            }
            j++;
        }

        if (isSafeToPaint) {
            colors[indexV] = i;
            if (bestSolutionColoring(graph, matrixAdj, numPosibleColors, indexV + 1, colors)) {
                return true;
            }
            else {
                colors[indexV] = 0;
            }
        }
    }
    return false;
}

void colorGraph(Graph* graph) {
    if (!graph->isEmpty()){
        int numV = graph->getNumVertices();
        int sizeMatrix = numV * numV;
        bool* matrixAdj = new bool[sizeMatrix];

        Vertex* vertex = graph->getHead();
        Vertex* secondVertex;
        int row = 0;
        while (vertex != nullptr) {
            int colum = 0;
            secondVertex = graph->getHead();
            while (secondVertex != nullptr) {
                int indexMatrix = row * numV + colum;
                if (graph->getWeight(vertex, secondVertex) > 0) {
                    matrixAdj[indexMatrix] = 1;
                }
                else {
                    matrixAdj[indexMatrix] = 0;
                }
                colum++;
                secondVertex = graph->getNextVert(secondVertex);
            }
            row++;
            vertex = graph->getNextVert(vertex);
        }

        int* colors = new int[numV];
        for (int i = 1; i < numV; i++) {
            Vertex* tempV = graph->getHead();
            if (bestSolutionColoring(graph, matrixAdj, i, 0, colors)) {
                for (int j = 0; j < numV; j++) {
                    std::cout << graph->getLabel(tempV) << " con color : " << colors[j] << std::endl;
                    tempV = graph->getNextVert(tempV);
                }
                break;
            }
        }

        delete[] matrixAdj;
        delete[] colors;
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
                << graph->getLabel(edges[i].first) << " -> "
                << graph->getLabel(edges[i].second) << " con peso: "
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
                << graph->getLabel(outputEdges[i].first) << " -> "
                << graph->getLabel(outputEdges[i].second) << " con peso: "
                << graph->getWeight(outputEdges[i].first, outputEdges[i].second) << std::endl;
        }
    }
}

// Algoritmo K.
void isolateVertex(Graph* graph, Vertex* vertex) {
    Vertex* adjVertex = graph->getFirstAdj(vertex);
    while (adjVertex != nullptr) {
        graph->deleteEdge(vertex, adjVertex);
        adjVertex = graph->getFirstAdj(vertex);
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

void clear(){
#if defined _WIN32
    system("cls");
#elif defined (__LINUX__) || defined(__gnu_linux__) || defined(__linux__)
    system("clear");
#elif defined (__APPLE__)
    system("clear");
#endif
}

void printGraphNames() {
    for (auto it = graphRegistry.begin(); it != graphRegistry.end(); ++it) {
        std::cout << "\t" << it->first << std::endl;
    }
}

void pressEnter() {
    std::cout << "\tIngrese Enter para continuar..." << std::endl;
    _getch();
}

Vertex* getValidVertex(Graph* graph, std::string message) {
    Vertex* vertex = nullptr;
    bool valid = false;
    char character;
    while (!valid) {
        std::cout << message << std::endl;
        std::cin >> character;
        vertex = getVertex(graph, character);
        if (vertex == nullptr) {
            std::cout << "\tNo se pudo encontrar ese valor..." << std::endl;
        }
        else {
            valid = true;
        }
    }
    return vertex;
}

std::string getGraphName() {
    bool validName = false;
    std::string name = "";
    while (!validName) {
        std::cout << "\tIngrese el nombre del Grafo: " << std::endl;
        std::cin >> name;
        if (!(graphRegistry.find(name) == graphRegistry.end())) {
            validName = true;
        }
    }
    return name;
}

void showMenu() {

    setlocale(LC_ALL, "spanish");
    SetConsoleCP(1252);
    SetConsoleOutputCP(1252);

    bool finished = false;

    int option;

    do {
        clear();

        std::cout << "| Menú de prueba para operadores y algoritmos básicos del modelo Grafo." << std::endl;
        std::cout << "| Grafos actuales:" << std::endl;
        printGraphNames();
        std::cout << "| =====================================================================" << std::endl;
        std::cout << "| 1) Iniciar.\t\t16) TieneCiclosAncho." << std::endl;
        std::cout << "| 2) Destruir.\t\t17) TieneCiclosProfundidad." << std::endl;
        std::cout << "| 3) Vaciar.\t\t18) esConexoProfundidad." << std::endl;
        std::cout << "| 4) Vacio.\t\t19) esConexoWarshall." << std::endl;
        std::cout << "| 5) AgregarVertice.\t20) Dijkstra." << std::endl;
        std::cout << "| 6) EliminarVertice.\t21) Floyd." << std::endl;
        std::cout << "| 7) ModificarEtiqueta.\t22) Hamilton." << std::endl;
        std::cout << "| 8) AgregarArista.\t23) Colorear." << std::endl;
        std::cout << "| 9) EliminarArista.\t24) Prim." << std::endl;
        std::cout << "| 10) ModificarPeso.\t25) Kruskal." << std::endl;
        std::cout << "| 11) ObtenerPeso.\t26) AislarVertice." << std::endl;
        std::cout << "| 12) ExisteArista.\t27) SumaPesos." << std::endl;
        std::cout << "| 13) NumAristas. \t28) ImprimirGrafo." << std::endl;
        std::cout << "| 14) NumVertices." << std::endl;
        std::cout << "| 15) NumVerticesAdyacentes." << std::endl;

        std::cout << std::endl;
        std::cout << "| 0) Salir." << std::endl;
        std::cout << "| =====================================================================" << std::endl;
        std::cin >> option;

        Graph* temp_graph;
        Vertex* temp_vertex1 = nullptr;
        Vertex* temp_vertex2 = nullptr;
        char character;
        double weigth;
        std::string name = "";

        switch (option)
        {
        case 0:
            finished = true;

            //delete temp_graph;

            break;
        case 1:
            std::cout << "\tNombre del Grafo: " << std::endl;
            std::cin >> name;
            temp_graph = new Graph();
            graphRegistry.insert({ name, temp_graph });
            pressEnter();
            break;
        case 2:
            // Destruir.
            std::cout << "\tNombre del Grafo: " << std::endl;
            std::cin >> name;
            temp_graph = graphRegistry.find(name)->second;
            graphRegistry.erase(graphRegistry.find(name));
            break;
        case 3:
            // Vaciar
            temp_graph = graphRegistry.find(getGraphName())->second;
            temp_graph->clear();
            pressEnter();
            break;
        case 4:
            // Vacio.
            temp_graph = graphRegistry.find(getGraphName())->second;
            if (temp_graph->isEmpty())
                std::cout << "\t El grafo esta vacio." << std::endl;
            else
                std::cout << "\t El grafo no esta vacio." << std::endl;
            pressEnter();
            break;
        case 5:
            // AgregarVertice.
            temp_graph = graphRegistry.find(getGraphName())->second;
            std::cout << "\tIngrese el valor: " << std::endl;
            std::cin >> character;
            temp_graph->addVert(character);
            pressEnter();
            break;
        case 6:
            // EliminarVertice.
            temp_graph = graphRegistry.find(getGraphName())->second;
            temp_vertex1 = getValidVertex(temp_graph, "\tIngrese el vertice a borrar: ");
            temp_graph->deleteVert(temp_vertex1);
            break;
        case 7:
            // ModificarEtiqueta.
            temp_graph = graphRegistry.find(getGraphName())->second;
            temp_vertex1 = getValidVertex(temp_graph, "\tIngrese el valor a cambiar: ");
            std::cout << "\tIngrese el nuevo valor: " << std::endl;
            std::cin >> character;
            temp_graph->changeLabel(temp_vertex1, character);
            break;
        case 8:
            // AgregarArista.
            temp_graph = graphRegistry.find(getGraphName())->second;
            temp_vertex1 = getValidVertex(temp_graph, "\tIngrese el primer vertice: ");
            temp_vertex2 = getValidVertex(temp_graph, "\tIngrese el segundo vertice: ");
            std::cout << "\tIngrese el peso: " << std::endl;
            std::cin >> weigth;
            temp_graph->createEdge(temp_vertex1, temp_vertex2, weigth);
            break;
        case 9:
            // EliminarArista.
            temp_graph = graphRegistry.find(getGraphName())->second;
            temp_vertex1 = getValidVertex(temp_graph, "\tIngrese el primer vertice: ");
            temp_vertex2 = getValidVertex(temp_graph, "\tIngrese el segundo vertice: ");

            temp_graph->deleteEdge(temp_vertex1, temp_vertex2);
            break;
        case 10:
            // ModificarPeso
            temp_graph = graphRegistry.find(getGraphName())->second;
            temp_vertex1 = getValidVertex(temp_graph, "\tIngrese el primer vertice: ");
            temp_vertex2 = getValidVertex(temp_graph, "\tIngrese el segundo vertice: ");
            std::cout << "\tIngrese el nuevo peso: " << std::endl;
            std::cin >> weigth;

            temp_graph->changeWeight(temp_vertex1, temp_vertex2, weigth);
            break;
        case 11:
            // ObtenerPeso			
            std::cout << "\tObtener peso..." << std::endl;
            temp_graph = graphRegistry.find(getGraphName())->second;
            temp_vertex1 = getValidVertex(temp_graph, "\tIngrese el primer vertice: ");
            temp_vertex2 = getValidVertex(temp_graph, "\tIngrese el segundo vertice: ");

            std::cout << "\tEl peso es: " << temp_graph->getWeight(temp_vertex1, temp_vertex2);
            pressEnter();
            break;
        case 12:
            // ExisteArista.
            temp_graph = graphRegistry.find(getGraphName())->second;
            temp_vertex1 = getValidVertex(temp_graph, "\tIngrese el primer vertice: ");
            temp_vertex2 = getValidVertex(temp_graph, "\tIngrese el segundo vertice: ");

            if (temp_graph->isEdge(temp_vertex1, temp_vertex2))
                std::cout << "\t Si existe esa arista." << std::endl;
            else
                std::cout << "\t No existe esa arista." << std::endl;
            pressEnter();
            break;
        case 13:
            // NumAristas.
            std::cout << "\tObtener numero de aristas..." << std::endl;
            temp_graph = graphRegistry.find(getGraphName())->second;
            std::cout << "\t La cantidad de aristas es: " << temp_graph->getNumEdges() << std::endl;
            pressEnter();
            break;
        case 14:
            // NumVertices.
            std::cout << "\tObtener numero de vertices..." << std::endl;
            temp_graph = graphRegistry.find(getGraphName())->second;
            std::cout << "\t La cantidad de vertices es: " << temp_graph->getNumVertices() << std::endl;
            pressEnter();
            break;
        case 15:
            // NumVerticesAdj.
            temp_graph = graphRegistry.find(getGraphName())->second;
            temp_vertex1 = getValidVertex(temp_graph, "\tIngrese el valor: ");
            std::cout << "\t La cantidad de vertices adyacentes es: " << temp_graph->getNumAdjVertices(temp_vertex1) << std::endl;
            pressEnter();
            break;
        case 16:
            // TieneCiclosAncho.
            temp_graph = graphRegistry.find(getGraphName())->second;
            if (hasCycleWidth(temp_graph))
                std::cout << "\t Si tiene ciclos." << std::endl;
            else
                std::cout << "\t No tiene ciclos." << std::endl;
            pressEnter();
            break;
        case 17:
            // TieneCiclosProfundidad
            temp_graph = graphRegistry.find(getGraphName())->second;
            if (cicleGraphDepth(temp_graph))
                std::cout << "\t Si tiene ciclos." << std::endl;
            else
                std::cout << "\t No tiene ciclos." << std::endl;
            pressEnter();
            break;
        case 18:
            // esConexoProfundidad 
            temp_graph = graphRegistry.find(getGraphName())->second;
            if (isRelatedDepthFirst(temp_graph))
                std::cout << "\t Si es conexo." << std::endl;
            else
                std::cout << "\t No es conexo." << std::endl;
            pressEnter();
            break;
        case 19:
            // esConexoWarshall
            temp_graph = graphRegistry.find(getGraphName())->second;
            if (isRelatedWarshall(temp_graph))
                std::cout << "\t Si es conexo." << std::endl;
            else
                std::cout << "\t No es conexo." << std::endl;
            pressEnter();
            break;
        case 20:
            // Dijkstra.
            temp_graph = graphRegistry.find(getGraphName())->second;
            temp_vertex1 = getValidVertex(temp_graph, "\tIngrese el valor: ");
            dijkstra(temp_graph, temp_vertex1);
            pressEnter();
            break;
        case 21:
            // Floyd.
            temp_graph = graphRegistry.find(getGraphName())->second;
            floyd(temp_graph);
            pressEnter();
            break;
        case 22:
            // Hamilton.
            temp_graph = graphRegistry.find(getGraphName())->second;
            hamilton(temp_graph);
            pressEnter();
            break;
        case 23:
            // Colorear.
            temp_graph = graphRegistry.find(getGraphName())->second;
            colorGraph(temp_graph);
            pressEnter();
            break;
        case 24:
            // Prim.
            temp_graph = graphRegistry.find(getGraphName())->second;
            prim(temp_graph);
            pressEnter();
            break;
        case 25:
            // Kruskal.
            temp_graph = graphRegistry.find(getGraphName())->second;
            kruskal(temp_graph);
            pressEnter();
            break;
        case 26:
            // ArislarVertice.
            temp_graph = graphRegistry.find(getGraphName())->second;
            temp_vertex1 = getValidVertex(temp_graph, "\tIngrese el valor: ");
            isolateVertex(temp_graph, temp_vertex1);
            
            break;
        case 27:
            // SumaPesos.
            temp_graph = graphRegistry.find(getGraphName())->second;
            std::cout << "\t El peso total es: " << totalWeight(temp_graph) << std::endl;
            pressEnter();
            break;

        case 28:
            // Imprimir.
            temp_graph = graphRegistry.find(getGraphName())->second;
            printGraphWidth(temp_graph);
            pressEnter();
            break;

        default:
            break;
        }

    } while (!finished);

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

    //Vertex* vA = test_graph->addVert('A');
    //Vertex* vB = test_graph->addVert('B');
    //Vertex* vC = test_graph->addVert('C');
    //Vertex* vD = test_graph->addVert('D');

    //test_graph->createEdge(vA, vC, 2.0);
    //test_graph->createEdge(vA, vB, 20.0);
    //test_graph->createEdge(vA, vD, 4.0);
    //test_graph->createEdge(vB, vD, 3.0);
    //test_graph->createEdge(vB, vC, 10.0);
    //test_graph->createEdge(vD, vC, 5.0);

    graphRegistry.insert({ "gA", test_graph });

    //printGraphWidth(test_graph);

    ////isolateVertex(test_graph, vF);
    ////test_graph->deleteEdge(vD, vB);
    //test_graph->deleteVert(vA);
    //std::cout << std::endl;

    //printGraphWidth(test_graph);

    showMenu();

    return 0;
}