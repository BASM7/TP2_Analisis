#pragma once

#include <iostream>
/*
* TP1 - Analisis de Algoritmos.
* @author B93986 Luis Alfonso Jim�nez
* @author B95346 Jes�s Alonso Moreno Montero
* @author B95092 V�ctor Jes�s Mora Abarca
*/

#pragma once

int const SIZE = 10; //cambiar para probar grafos mas grandes.

class Vertex {
	char label;
	int index;
public:

	Vertex() {
	}

	~Vertex() {}


	char getLabel() { return label; };
	void setLabel(char new_label) { label = new_label; };
	int getIndex() { return index; };
	void setIndex(int new_index) { index = new_index; };;
};

class Graph {
public:
	Graph() {

	};

	~Graph() {};

	void clear();
	bool isEmpty();
	Vertex* addVert(char label);
	void deleteVert(int vertex);
	void changeLabel(Vertex* vertex, char new_label);
	char getLabel(Vertex* vertex);
	void createEdge(Vertex* vertex1, Vertex* vertex2, double weight);
	void deleteEdge(Vertex* vertex1, Vertex* vertex2);
	void changeWeight(Vertex* vertex1, Vertex* vertex2, double new_weight);
	double getWeight(Vertex* vertex1, Vertex* vertex2);
	Vertex* getHead();
	Vertex* getNextVert(Vertex* vertex);
	Vertex* getFirstAdj(Vertex* vertex);
	Vertex* getNextAdj(Vertex* vertex, Vertex* vertex2);
	bool isEdge(Vertex* vertex1, Vertex* vertex2);
	int getNumEdges() { return cant_edges; };
	int getNumVertices() { return cant_vert; };
	int getNumAdjVertices(Vertex* vertex);

	void imprimir();
private:
	double matrix[SIZE * SIZE];
	Vertex* names[SIZE];
	int cant_vert = 0;
	int cant_edges = 0;
};


bool Graph::isEmpty() {
	return cant_vert == 0;
}

void Graph::clear() {
	cant_vert = 0;
	cant_edges = 0;
}

Vertex* Graph::addVert(char label) {
	if (cant_vert < SIZE) {
		Vertex* vertex = new Vertex();
		vertex->setIndex(cant_vert);
		vertex->setLabel(label);
		names[cant_vert] = vertex;
		for (int index = cant_vert * SIZE; index < (cant_vert + 1) * SIZE; index++) {
			matrix[index] = -1;
		}
		cant_vert++;

		return vertex;
	}
	return nullptr;
}

void Graph::deleteVert(int vertex) { // Not rdy
	if (vertex - 1 < cant_vert) {
		std::cout << cant_vert << std::endl;
		//Corrimiento de las filas
		for (int row = vertex - 1; row < cant_vert; row++) {
			names[row] = names[row + 1];
			//matrix[row] = matrix[row + SIZE];
		}


	}

}

void Graph::changeLabel(Vertex* vertex, char new_label) {
	names[vertex->getIndex()]->setLabel(new_label);
}

char Graph::getLabel(Vertex* vertex) {
	return names[vertex->getIndex()]->getLabel();
}

void Graph::createEdge(Vertex* vertex1, Vertex* vertex2, double weight) {
	if (vertex1->getLabel() != vertex2->getLabel()) {
		matrix[(vertex1->getIndex() * SIZE) + vertex2->getIndex()] = weight;
		matrix[(vertex2->getIndex() * SIZE) + vertex1->getIndex()] = weight;
		cant_edges++;
	}
}

void Graph::deleteEdge(Vertex* vertex1, Vertex* vertex2) {
	if (vertex1->getIndex() != vertex2->getIndex() && vertex1->getIndex() < cant_vert && vertex2->getIndex() < cant_vert) {
		matrix[(vertex1->getIndex()* SIZE)  + vertex2->getIndex()] = -1.0;
		matrix[(vertex2->getIndex()* SIZE)  + vertex1->getIndex()] = -1.0;
		cant_edges--;
	}
}

void Graph::changeWeight(Vertex* vertex1, Vertex* vertex2, double new_weight) {
	if (vertex1->getIndex() != vertex2->getIndex() && vertex1->getIndex() < cant_vert && vertex2->getIndex() < cant_vert) {
		matrix[(vertex1->getIndex()* SIZE) + vertex2->getIndex()] = new_weight;
		matrix[(vertex2->getIndex()* SIZE) + vertex1->getIndex()] = new_weight;
	}
}

double Graph::getWeight(Vertex* vertex1, Vertex* vertex2) {
	if (vertex1->getIndex() != vertex2->getIndex() && vertex1->getIndex() < cant_vert && vertex2->getIndex() < cant_vert) {
		return matrix[(vertex1->getIndex()* SIZE)  + vertex2->getIndex()];
	}
	return -1.0;
}

Vertex* Graph::getHead() {
	if (cant_vert > 0) {
		return names[0];
	}
	return nullptr;
}

Vertex* Graph::getNextVert(Vertex* vertex) {
	if (vertex->getIndex() + 1 < cant_vert) {
		return names[vertex->getIndex() + 1];
	}
	return nullptr;
}

Vertex* Graph::getFirstAdj(Vertex* vertex) {
	Vertex* vertex_toReturn = nullptr;
	int row = vertex->getIndex() * SIZE;
	int colum = 0;
	while (colum < SIZE) {
		if (matrix[row + colum] != -1.0) {
			vertex_toReturn = names[colum];
			break;
		}
		colum++;
	}
	return vertex_toReturn;
}


Vertex* Graph::getNextAdj(Vertex* vertex, Vertex* vertex2) {
	Vertex* vertex_toReturn = nullptr;
	int row = vertex->getIndex() * SIZE;
	int colum = vertex2->getIndex() + 1;
	while (colum < SIZE) {
		if (matrix[row + colum] != -1.0) {
			vertex_toReturn = names[colum];
			break;
		}
		colum++;
	}
	return vertex_toReturn;
}

bool Graph::isEdge(Vertex* vertex1, Vertex* vertex2) {
	return matrix[(vertex1->getIndex()) * cant_vert + (vertex2->getIndex())] != 1.0;
}

int Graph::getNumAdjVertices(Vertex* vertex) {
	int num_adj = 0;
	int row = vertex->getIndex() * SIZE;
	int colum = 0;
	while (colum < SIZE) {
		if (matrix[row + colum] != -1.0) {
			num_adj++;
		}
		colum++;
	}
	return num_adj;
}

void Graph::imprimir() {
	for (int index = 0; index < cant_vert * SIZE; index++) {
		if (index % SIZE == 0) {
			std::cout << std::endl;
			std::cout << names[index / SIZE]->getLabel() << "     ";
		}
		if (matrix[index] != -1) {
			std::cout << matrix[index];
		}
		else {
			std::cout << "0";
		}
		std::cout << " ";
	}
}
