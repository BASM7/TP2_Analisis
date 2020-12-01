#pragma once

#include <unordered_set>
#include <algorithm>
#include <stack>
#include <queue>

#include "GraphB.h"

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
    Vertex* vA = test_graph->addVert('A');
    Vertex* vB = test_graph->addVert('B');
    Vertex* vC = test_graph->addVert('C');
    Vertex* vD = test_graph->addVert('D');

    test_graph->createEdge(vA, vB, 5.0);
    test_graph->createEdge(vB, vC, 4.0);
    test_graph->createEdge(vC, vA, 9.0);

    printGraphWidth(test_graph);
    
    std::cout << totalWeight(test_graph) << std::endl;
    std::cout << getVertex(test_graph, 'B')->getValue() << std::endl;
    std::cout << getVertex(test_graph, 'B')->nextAdj->to->getValue() << std::endl;

    return 0;
}