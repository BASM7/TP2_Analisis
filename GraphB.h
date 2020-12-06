#pragma once

#include <iostream>
#include <string>
#include <sstream>

class Edge;

class Vertex {
    char value;
public:    
    Vertex* nextVertex;
    Edge* nextAdj;

    Vertex(char value, Vertex* nextVertex) {
        this->value = value;
        this->nextVertex = nextVertex;
        this->nextAdj = nullptr;
    }

    char getValue() { return this->value; };
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

public:    

    Graph() {}

    bool isEmpty() { return cantVertex == 0; };
    void clear();
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

    ~Graph() {}
};

Vertex* Graph::addVert(char newLabel) {
    Vertex* newVertex = new Vertex(newLabel, this->head);
    this->head = newVertex;
    this->cantVertex++;
    return newVertex;
}

void Graph::deleteVert(Vertex* vertexToDelete) {

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
    if (vertex)
        if (vertex->nextAdj)
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
                if (tempEdge->nextEdge) {
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



