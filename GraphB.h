#pragma once

#include <iostream>
#include <string>
#include <sstream>

class Edge;

class Vertex {
    char value;
public:    
    Vertex* nextVertex = nullptr;
    Edge* nextAdj = nullptr;

    Vertex(char value, Vertex* nextVertex) {
        this->value = value;
        this->nextVertex = nextVertex;
        this->nextAdj = nullptr;
    }

    char getValue() { 
        return this->value; 
    };
    void setValue(char newValue) {
        this->value = newValue;
    }

    ~Vertex() {
        delete nextAdj;
        delete nextVertex;
    };
};

class Edge {
public:
    Vertex* to = nullptr;
    double weight = 0.0;
    Edge* nextEdge = nullptr;

    Edge(Vertex* to, double weight, Edge* nextEdge) {
        this->to = to;
        this->weight = weight;
        this->nextEdge = nextEdge;
    }

    ~Edge() {};
};

class Graph{
    int cantVertex = 0;
    int cantEdges = 0;
    Vertex* head = nullptr;
    void isolate(Vertex* vertex);

public:    

    Graph() {}

    bool isEmpty() { return cantVertex == 0; };
    void clear() {
        this->head = nullptr;
        this->cantEdges = 0;
        this->cantVertex = 0;
    };

    Vertex* addVert(char newLabel);
    void deleteVert(Vertex* vertexToDelete);
    void changeLabel(Vertex* vertexToChange, char newLabel);
    void createEdge(Vertex* vertex1, Vertex* vertex2, double newWeight);
    void deleteEdge(Vertex* vertex1, Vertex* vertex2);
    void changeWeight(Vertex* vertex1, Vertex* vertex2, double newWeight);
    double getWeight(Vertex* vertex1, Vertex* vertex2);
    char getLabel(Vertex* vertex) { return vertex->getValue(); };
    Vertex* getHead() { return this->head; };
    Vertex* getNextVert(Vertex* vertex);
    Vertex* getFirstAdj(Vertex* vertex);
    Vertex* getNextAdj(Vertex* vertex1, Vertex* vertex2);
    bool isEdge(Vertex* vertex1, Vertex* vertex2);
    int getNumEdges() { return this->cantEdges; };
    int getNumVertices() { return this->cantVertex; };
    int getNumAdjVertices(Vertex* vertex);

    ~Graph() { }
};

void Graph::isolate(Vertex* vertex) {
    Vertex* adjVertex = this->getFirstAdj(vertex);
    while (adjVertex != nullptr) {
        this->deleteEdge(vertex, adjVertex);
        adjVertex = this->getFirstAdj(vertex);
    }
}

Vertex* Graph::addVert(char newLabel) {
    Vertex* newVertex = new Vertex(newLabel, this->head);
    this->head = newVertex;
    this->cantVertex++;
    return newVertex;
}

void Graph::deleteVert(Vertex* vertexToDelete) {
    if (vertexToDelete != nullptr){
        Vertex* tempV = this->getHead();
        Vertex* previousV = nullptr;
        bool found = false;
        while (tempV != nullptr && !found) {
            
            if (tempV == vertexToDelete) {
                found = true;
               

                if (previousV != nullptr) {
                    previousV->nextVertex = tempV->nextVertex;
                }
                else {
                    this->head = tempV->nextVertex;
                }
                isolate(tempV);
                cantVertex--;
            }
            previousV = tempV;
            tempV = this->getNextVert(tempV);
        }
        cantVertex--;
    }
}


void Graph::changeLabel(Vertex* vertexToChange, char newLabel) {
    if (vertexToChange)
        vertexToChange->setValue(newLabel);
}

void Graph::createEdge(Vertex* vertex1, Vertex* vertex2, double newWeight) {
    if (vertex1 && vertex2) {
        Edge* newEdge1 = new Edge(vertex1, newWeight, vertex2->nextAdj);
        vertex2->nextAdj = newEdge1;
        Edge* newEdge2 = new Edge(vertex2, newWeight, vertex1->nextAdj);
        vertex1->nextAdj = newEdge2;
        this->cantEdges++;
    }
}

void Graph::deleteEdge(Vertex* vertex1, Vertex* vertex2) {
    Edge* previousEdge = nullptr;
    Edge* tempEdge = vertex1->nextAdj;
    bool done = false;
    while (tempEdge != nullptr && !done) {
        
        if (tempEdge->to == vertex2) {
            if (previousEdge != nullptr) {
                previousEdge->nextEdge = tempEdge->nextEdge;
            }
            else {
                vertex1->nextAdj = tempEdge->nextEdge;
            }
            //tempEdge->nextEdge = nullptr;
            delete tempEdge;

            done = true;
        }
        previousEdge = tempEdge;
        tempEdge = tempEdge->nextEdge;
    }

    previousEdge = nullptr;
    Edge* tempEdge2 = vertex2->nextAdj;
    done = false;
    while (tempEdge2 != nullptr && !done) {
        if (tempEdge2->to == vertex1) {

            if (previousEdge != nullptr) {
                previousEdge->nextEdge = tempEdge2->nextEdge;
            }
            else {
                vertex2->nextAdj = tempEdge2->nextEdge;
            }
            //tempEdge2->nextEdge = nullptr;
            delete tempEdge2;
            done = true;
        }
        previousEdge = tempEdge2;
        tempEdge2 = tempEdge2->nextEdge;
    }

    this->cantEdges--;
}

void Graph::changeWeight(Vertex* vertex1, Vertex* vertex2, double newWeight) {
    double weight = 0.0;
    bool done = false;
    if (vertex1 && vertex2) {
        Edge* tempEdge = vertex1->nextAdj;
        while (tempEdge != nullptr && !done) {
            if (tempEdge->to == vertex2) {
                weight = tempEdge->weight;
                done = true;
            }
            tempEdge = tempEdge->nextEdge;
        }
    }
}

double Graph::getWeight(Vertex* vertex1, Vertex* vertex2) {
    double weight = -1.0;
    bool done = false;
    if (vertex1 && vertex2) {
        Edge* tempEdge = vertex1->nextAdj;
        while (tempEdge != nullptr && !done) {
            if (tempEdge->to == vertex2) {
                weight = tempEdge->weight;
                done = true;
            }
            tempEdge = tempEdge->nextEdge;
        }
    }
    return weight;
}

Vertex* Graph::getNextVert(Vertex* vertex) {
    if (vertex != nullptr)
        return vertex->nextVertex;
    return nullptr;
}

Vertex* Graph::getFirstAdj(Vertex* vertex) {
    if (vertex != nullptr)
        if (vertex->nextAdj != nullptr)
            return vertex->nextAdj->to;
    return nullptr;
}

Vertex* Graph::getNextAdj(Vertex* vertex1, Vertex* vertex2) {
    bool found = false;
    if (vertex1 && vertex2) {
        Edge* tempEdge = vertex1->nextAdj;
        while (tempEdge != nullptr && !found) {
            if (tempEdge->to == vertex2) {
                found = true;
                if (tempEdge->nextEdge != nullptr) {
                    return tempEdge->nextEdge->to;
                }
            }
            tempEdge = tempEdge->nextEdge;
        }
    }
    return nullptr;
}

bool Graph::isEdge(Vertex* vertex1, Vertex* vertex2) {
    bool found = false;
    if (vertex1 && vertex2) {
        Edge* tempEdge = vertex1->nextAdj;
        while (tempEdge != nullptr && !found) {
            if (tempEdge->to == vertex2) {
                found = true;
            }
            tempEdge = tempEdge->nextEdge;
        }
    }
    return found;
}

int Graph::getNumAdjVertices(Vertex* vertex) {
    int cant = 0;

    return cant;
}



