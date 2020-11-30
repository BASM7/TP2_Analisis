#pragma once

#include <map>
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
                            queue.push(adjVertex);
                            visitedVertices.insert(adjVertex);
                        }

                        std::cout << graph->getLabel(tempVertex) << " -> ";
                        std::cout << graph->getWeight(tempVertex, adjVertex);
                        std::cout << " -> " << graph->getLabel(adjVertex) << std::endl;

                        adjVertex = graph->getNextAdj(tempVertex, adjVertex);
                    }
                }
            }
            currentVertex = graph->getNextVert(currentVertex);
        }
    }
}

int main()
{
    Graph* test_graph = new Graph();


    return 0;
}