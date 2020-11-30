#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <sstream>

class Vertex {
public:
    char value;
    Vertex* nextVertex;

    Vertex(char value, Vertex* nextVertex) {
        this->value = value;
        this->nextVertex = nextVertex;
    }

    ~Vertex() {};
};

class Edge {
public:
    Vertex* to;
    double weight;
    Edge* nextEdge;

    Edge(Vertex* to, double weight, Edge* nextEdge) {
        this->to = to;
        this->weight = weight;
        this->nextEdge = nextEdge;
    }

    ~Edge() {};
};

class Graph{
    int cantVertex;
    int cantEdges;
    Vertex* head;

public:    

    Graph() {
        head = nullptr;
        cantVertex = 0;
        cantEdges = 0;
    }

    bool isEmpty() { return cantVertex == 0; };
    Vertex* addVert(char newLabel);
    void deleteVert(Vertex* vertexToDelete);
    void changeLabel(Vertex* vertexToChange, char newLabel);
    void createEdge(Vertex* vertex1, Vertex* vertex2);
    void deleteEdge(Vertex* vertex1, Vertex* vertex2);
    void changeWeight(Vertex* vertex1, Vertex* vertex2, double newWeight);
    double getWeight(Vertex* vertex1, Vertex* vertex2);
    char getLabel(Vertex* vertex) { return vertex->value; };
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

}

void Graph::deleteVert(Vertex* vertexToDelete) {

}


void Graph::changeLabel(Vertex* vertexToChange, char newLabel) {

}

void Graph::createEdge(Vertex* vertex1, Vertex* vertex2) {

}

void Graph::deleteEdge(Vertex* vertex1, Vertex* vertex2) {

}

void Graph::changeWeight(Vertex* vertex1, Vertex* vertex2, double newWeight) {

}

double Graph::getWeight(Vertex* vertex1, Vertex* vertex2) {

}

Vertex* Graph::getNextVert(Vertex* vertex) {

}

Vertex* Graph::getFirstAdj(Vertex* vertex) {

}

Vertex* Graph::getNextAdj(Vertex* vertex, Vertex* vertex2) {

}

bool Graph::isEdge(Vertex* vertex1, Vertex* vertex2) {

}

int Graph::getNumAdjVertices(Vertex* vertex) {

}



