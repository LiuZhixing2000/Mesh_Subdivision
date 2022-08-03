#include "object.h"

std::vector<std::string> seperate_string(std::string origin) {
	std::vector<std::string> result;
	std::stringstream ss(origin);
	while (ss >> origin) result.push_back(origin);
	return result;
}

/* judge if s is on the left of edge p to q */
int toLeftTest(Vertex p, Vertex q, Vertex s) {
	double determinant = p.x() * q.y() - p.y() * q.x()
		+ q.x() * s.y() - q.y() * s.x()
		+ s.x() * p.y() - s.y() * p.x();
	if (abs(determinant - 0) <= ER)
		return 0;
	else if (determinant > 0)
		return 1;
	else
		return -1;
}

int Object::readObjFile(const char* filename) {
	std::ifstream objFile;
	objFile.open(filename);
	if (!objFile.is_open()) {
		std::cout << "No such file. - " << filename << std::endl;
		return -1;
	}
	else
		std::cout << "Successfully open file - " << filename << std::endl;
	char buffer[BUFFER_LENGTH];
	while (!objFile.eof()) {
		objFile.getline(buffer, BUFFER_LENGTH);
		if (buffer[0] == 'v') {
			std::vector<std::string> words = seperate_string(std::string(buffer));
			if (words.size() < 4) {
				continue;
			}
			for (int ii = 1; ii < 4; ii++) {
				int wordLen = 0;
				while (wordLen < words[ii].length()) {
					if (words[ii][wordLen] == '/') {
						break;
					}
					else
						wordLen++;
				}
				vertices.push_back(stod(words[ii].substr(0, wordLen)));
			}
		}
		if (buffer[0] == 'f') {
			std::vector<std::string> words = seperate_string(std::string(buffer));
			if (words.size() < 4) {
				continue;
			}
			polygons.push_back(std::vector<int>(0));
			for (int ii = 1; ii < words.size(); ii++) {
				int wordLen = 0;
				while (wordLen < words[ii].length()) {
					if (words[ii][wordLen] == '/') {
						break;
					}
					else
						wordLen++;
				}
				polygons[polygons.size() -1].push_back(stoi(words[ii].substr(0, wordLen)) - 1);
			}
		}
	}
	std::cout << "Read file end. Vertices: " << vertices.size() << " Polygons: " << polygons.size() << std::endl;
	return 1;
}

int Object::DooSabin(int times) {
	std::cout << "Doo-Sabin Surface Subdivision Start..." << std::endl;
	for (int n = 0; n < times; n++) {
		clock_t start_t, end_t;
		double totalTime = -1;
		start_t = clock();
		newVertices.clear();
		newPolygons.clear();
		newVertices.resize(0);
		newPolygons.resize(0);
		edges.clear();
		vertex2Polygons.clear();
		edge2Polygons.clear();
		polygonOffset.clear();
		establishPolygonOffset();
		establishEdges();
		establishVertex2Polygons();
		establishEdge2Polygons();
		// ��ÿһ���������Ƭ
		// �����ڲ��γ��µĵ�
		for (int k = 0; k < polygons.size(); k++) {
			newPolygons.push_back(std::vector<int>(0));
			// face point
			double fpxx = 0, fpyy = 0, fpzz = 0;
			for (int i = 0; i < polygons[k].size(); i++) {
				fpxx += vertices[3 * polygons[k][i]];
				fpyy += vertices[3 * polygons[k][i] + 1];
				fpzz += vertices[3 * polygons[k][i] + 2];
			}
			fpxx /= polygons[k].size();
			fpyy /= polygons[k].size();
			fpzz /= polygons[k].size();
			// ÿ��һ���㣬��������һ����
			for (int i = 0; i < polygons[k].size(); i++) {
				// �����ɵĵ�vi
				double vixx = 0, viyy = 0, vizz = 0;
				// edge point
				int vp, vs;
				if (i == 0) {
					vp = polygons[k][polygons[k].size() - 1];
					vs = polygons[k][1];
				}
				else if (i == polygons[k].size() - 1) {
					vp = polygons[k][polygons[k].size() - 2];
					vs = polygons[k][0];
				}
				else {
					vp = polygons[k][i - 1];
					vs = polygons[k][i + 1];
				}
				double epxx = 0, epyy = 0, epzz = 0;
				epxx = (vertices[3 * polygons[k][i]] + vertices[3 * vp]) / 2 + (vertices[3 * polygons[k][i]] + vertices[3 * vs]) / 2;
				epyy = (vertices[3 * polygons[k][i] + 1] + vertices[3 * vp + 1]) / 2 + (vertices[3 * polygons[k][i] + 1] + vertices[3 * vs + 1]) / 2;
				epzz = (vertices[3 * polygons[k][i] + 2] + vertices[3 * vp + 2]) / 2 + (vertices[3 * polygons[k][i] + 2] + vertices[3 * vs + 2]) / 2;
				vixx = (epxx + fpxx + vertices[3 * polygons[k][i]]) / 4;
				viyy = (epyy + fpyy + vertices[3 * polygons[k][i] + 1]) / 4;
				vizz = (epzz + fpzz + vertices[3 * polygons[k][i] + 2]) / 4;
				newVertices.push_back(vixx);
				newVertices.push_back(viyy);
				newVertices.push_back(vizz);
				newPolygons[newPolygons.size() - 1].push_back(newVertices.size() / 3 - 1);
			}
		}

		// ��ÿ���߽��е��߲���
		int edgeId = 0;
		for (auto i = edges.begin(); i != edges.end(); i++) {
			newPolygons.push_back(std::vector<int>(0));
			auto k = i->begin();
			int v0 = *k; k++;
			int v1 = *k;
			int plg0 = edge2Polygons[edgeId][0];
			int plg1 = edge2Polygons[edgeId][1];
			edgeId++;
			int a = -1, b = -1, c = -1, d = -1;
			for (int j = 0; j < polygons[plg0].size(); j++) {
				if (polygons[plg0][j] == v0)
					a = j;
				if (polygons[plg0][j] == v1)
					b = j;
				if (a >= 0 && b >= 0)
					break;
			}
			for (int j = 0; j < polygons[plg1].size(); j++) {
				if (polygons[plg1][j] == v1)
					c = j;
				if (polygons[plg1][j] == v0)
					d = j;
				if (c >= 0 && d >= 0)
					break;
			}
			newPolygons[newPolygons.size() - 1].push_back(polygonOffset[plg0] + a);
			newPolygons[newPolygons.size() - 1].push_back(polygonOffset[plg0] + b);
			newPolygons[newPolygons.size() - 1].push_back(polygonOffset[plg1] + c);
			newPolygons[newPolygons.size() - 1].push_back(polygonOffset[plg1] + d);
		}
		// ��ÿ������е��ǲ���
		for (int i = 0; i < vertex2Polygons.size(); i++) {
			newPolygons.push_back(std::vector<int>(0));
			for (int j = 0; j < vertex2Polygons[i].size(); j++) {
				for (int k = 0; k < polygons[vertex2Polygons[i][j]].size(); k++) {
					if (polygons[vertex2Polygons[i][j]][k] == i)
						newPolygons[newPolygons.size() - 1].push_back(polygonOffset[vertex2Polygons[i][j]] + k);
				}
			}
		}
		vertices.resize(newVertices.size());
		for (int i = 0; i < vertices.size(); i++)
			vertices[i] = newVertices[i];
		polygons.clear();
		for (int i = 0; i < newPolygons.size(); i++) {
			polygons.push_back(std::vector<int>(0));
			for (int j = 0; j < newPolygons[i].size(); j++) {
				polygons[polygons.size() - 1].push_back(newPolygons[i][j]);
			}
			/*
			if (newPolygons[i].size() == 3) {
				polygons.push_back(std::vector<int>{newPolygons[i][0], newPolygons[i][1], newPolygons[i][2]});
			}
			else {
				for (int j = 2; j < newPolygons[i].size(); j++) {
					polygons.push_back(std::vector<int>({ newPolygons[i][0], newPolygons[i][j - 1], newPolygons[i][j] }));
				}
			}*/
		}
		newVertices.clear();
		newPolygons.clear();
		end_t = clock();
		totalTime = (double)(end_t - start_t) / CLOCKS_PER_SEC;
		std::cout << '\t' << "Subdivision times: " << n + 1 << '\t' << "Time cost: " << totalTime << 's' << std::endl;
	}
	std::cout << "Doo-Sabin Surface Subdivision Success." << std::endl;
	return 1;
}

int Object::CatmullClark(int times) {
	for (int time = 0; time < times; time++) {
		clock_t start_t, end_t;
		double totalTime = -1;
		start_t = clock();
		newVertices.clear();
		newPolygons.clear();
		newVertices.resize(0);
		newPolygons.resize(0);
		edges.clear();
		vertex2Polygons.clear();
		edge2Polygons.clear();
		polygonOffset.clear();
		vertex2Edges.clear();
		edgesV.clear();
		establishPolygonOffset();
		establishEdges();
		establishVertex2Polygons();
		establishEdge2Polygons();
		establishVertex2Edges();
		establishEdgesV();
		// face point
		std::vector<double> facePoints;
		for (int i = 0; i < polygons.size(); i++) {
			double fpxx = 0, fpyy = 0, fpzz = 0;
			for (int j = 0; j < polygons[i].size(); j++) {
				fpxx += vertices[3 * polygons[i][j]];
				fpyy += vertices[3 * polygons[i][j] + 1];
				fpzz += vertices[3 * polygons[i][j] + 2];
			}
			fpxx /= polygons[i].size();
			fpyy /= polygons[i].size();
			fpzz /= polygons[i].size();
			facePoints.push_back(fpxx);
			facePoints.push_back(fpyy);
			facePoints.push_back(fpzz);
		}
		// edge point
		std::vector<double> edgePoints;
		int edgeId = 0;
		for (auto i = edges.begin(); i != edges.end(); i++) {
			double epxx = 0, epyy = 0, epzz = 0;
			auto temp = i->begin();
			int v0 = *temp; temp++;
			int v1 = *temp;
			int p0 = edge2Polygons[edgeId][0];
			int p1 = edge2Polygons[edgeId][1];
			epxx = (facePoints[3 * p0] + facePoints[3 * p1] + (vertices[3 * v0] + vertices[3 * v1]) / 2) / 3;
			epyy = (facePoints[3 * p0 + 1] + facePoints[3 * p1 + 1] + (vertices[3 * v0 + 1] + vertices[3 * v1 + 1]) / 2) / 3;
			epzz = (facePoints[3 * p0 + 2] + facePoints[3 * p1 + 2] + (vertices[3 * v0 + 2] + vertices[3 * v1 + 2]) / 2) / 3;
			edgePoints.push_back(epxx);
			edgePoints.push_back(epyy);
			edgePoints.push_back(epzz);
			edgeId++;
		}
		// vertex point
		std::vector<double> vertexPoints;
		for (int i = 0; i < vertex2Polygons.size(); i++) {
			// Q: the average of the new face points of all faces adjacent to the original face point
			std::vector<int>& plgs = vertex2Polygons[i];
			double qxx = 0, qyy = 0, qzz = 0;
			for (int j = 0; j < plgs.size(); j++) {
				qxx += facePoints[3 * plgs[j]];
				qyy += facePoints[3 * plgs[j] + 1];
				qzz += facePoints[3 * plgs[j] + 2];
			}
			qxx /= plgs.size();
			qyy /= plgs.size();
			qzz /= plgs.size();
			// R: is the average of the midpoints of all original edges incident on the original vertex point
			std::vector<int>& es = vertex2Edges[i];
			double rxx = 0, ryy = 0, rzz = 0;
			for (int j = 0; j < es.size(); j++) {
				rxx += (vertices[3 * edgesV[es[j]][0]] + vertices[3 * edgesV[es[j]][1]]) / 2;
				ryy += (vertices[3 * edgesV[es[j]][0] + 1] + vertices[3 * edgesV[es[j]][1] + 1]) / 2;
				rzz += (vertices[3 * edgesV[es[j]][0] + 2] + vertices[3 * edgesV[es[j]][1] + 2]) / 2;
			}
			rxx /= es.size();
			ryy /= es.size();
			rzz /= es.size();
			// S: is the original vertex point
			double sxx = vertices[3 * i];
			double syy = vertices[3 * i + 1];
			double szz = vertices[3 * i + 2];
			// vertex point
			double vpxx, vpyy, vpzz;
			int n = vertex2Polygons[i].size();
			vpxx = (qxx + 2 * rxx + (n - 3)*sxx) / n;
			vpyy = (qyy + 2 * ryy + (n - 3)*syy) / n;
			vpzz = (qzz + 2 * rzz + (n - 3)*szz) / n;
			vertexPoints.push_back(vpxx);
			vertexPoints.push_back(vpyy);
			vertexPoints.push_back(vpzz);
		}
		std::vector< std::vector<int> > newEdges;
		for (int i = 0; i < facePoints.size(); i++) 
			newVertices.push_back(facePoints[i]);
		for (int i = 0; i < edgePoints.size(); i++)
			newVertices.push_back(edgePoints[i]);
		for (int i = 0; i < vertexPoints.size(); i++)
			newVertices.push_back(vertexPoints[i]);
		int fpNum = facePoints.size() / 3;
		int epNum = edgePoints.size() / 3;
		// Each new face point is connected to the new edge points of the edges defining the original face.
		// Each new vertex point is connected to the new edge points of all original edges incident on the original vertex point.
		for (int i = 0; i < vertex2Polygons.size(); i++) {
			for (int j = 0; j < vertex2Polygons[i].size(); j++) {
				int pId = vertex2Polygons[i][j];
				std::vector<int>& ves = vertex2Edges[i];
				int e0 = -1, e1 = -1;
				for (int k = 0; k < ves.size(); k++) {
					if (edge2Polygons[ves[k]][0] == pId || edge2Polygons[ves[k]][1] == pId) {
						if (e0 == -1)
							e0 = ves[k];
						else
							e1 = ves[k];
					}
					if (e0 >= 0 && e1 >= 0)
						break;
				}
				newPolygons.push_back(std::vector<int>({ pId, e0 + fpNum, i + fpNum + epNum, e1 + fpNum }));
			}
		}
		vertices.resize(newVertices.size());
		for (int i = 0; i < vertices.size(); i++)
			vertices[i] = newVertices[i];
		polygons.clear();
		for (int i = 0; i < newPolygons.size(); i++) {
			polygons.push_back(std::vector<int>(0));
			for (int j = 0; j < newPolygons[i].size(); j++) {
				polygons[polygons.size() - 1].push_back(newPolygons[i][j]);
			}
			/*
			if (newPolygons[i].size() == 3) {
				polygons.push_back(std::vector<int>{newPolygons[i][0], newPolygons[i][1], newPolygons[i][2]});
			}
			else {
				for (int j = 2; j < newPolygons[i].size(); j++) {
					polygons.push_back(std::vector<int>({ newPolygons[i][0], newPolygons[i][j - 1], newPolygons[i][j] }));
				}
			}*/
		}
		newVertices.clear();
		newPolygons.clear();
		end_t = clock();
		totalTime = (double)(end_t - start_t) / CLOCKS_PER_SEC;
		std::cout << '\t' << "Subdivision times: " << time + 1 << '\t' << "Time cost: " << totalTime << 's' << std::endl;
	}
	std::cout << "Catmull-Clark Surface Subdivision Success." << std::endl;
	return 1;
}

// ��������������Ƭ������ϸ��
int Object::Loop(int times) {
	for (int time = 0; time < times; time++) {
		clock_t start_t, end_t;
		double totalTime = -1;
		start_t = clock();
		newVertices.clear();
		newPolygons.clear();
		newVertices.resize(0);
		newPolygons.resize(0);
		edges.clear();
		vertex2Polygons.clear();
		edge2Polygons.clear();
		polygonOffset.clear();
		vertex2Edges.clear();
		edgesV.clear();
		polygon2Edges.clear();
		establishPolygonOffset();
		establishEdges();
		establishVertex2Polygons();
		establishEdge2Polygons();
		establishVertex2Edges();
		establishEdgesV();
		establishPolygon2Edges();
		// edge point
		std::vector<double> edgePoints;
		for (int i = 0; i < edgesV.size(); i++) {
			int plg0 = edge2Polygons[i][0];
			int plg1 = edge2Polygons[i][1];
			int vp0 = -1, vp1 = -1;
			for (int j = 0; j < 3; j++) {
				if (polygons[plg0][j] != edgesV[i][0] && polygons[plg0][j] != edgesV[i][1]) 
					vp0 = polygons[plg0][j];
				if (polygons[plg1][j] != edgesV[i][0] && polygons[plg1][j] != edgesV[i][1])
					vp1 = polygons[plg1][j];
				if (vp0 >= 0 && vp1 >= 0)
					break;
			}
			double epxx, epyy, epzz;
			epxx = vertices[3 * vp0] / 8 + vertices[3 * vp1] / 8 + vertices[3 * edgesV[i][0]] * 3 / 8 + vertices[3 * edgesV[i][1]] * 3 / 8;
			epyy = vertices[3 * vp0 + 1] / 8 + vertices[3 * vp1 + 1] / 8 + vertices[3 * edgesV[i][0] + 1] * 3 / 8 + vertices[3 * edgesV[i][1] + 1] * 3 / 8;
			epzz = vertices[3 * vp0 + 2] / 8 + vertices[3 * vp1 + 2] / 8 + vertices[3 * edgesV[i][0] + 2] * 3 / 8 + vertices[3 * edgesV[i][1] + 2] * 3 / 8;
			edgePoints.push_back(epxx);
			edgePoints.push_back(epyy);
			edgePoints.push_back(epzz);
		}
		// extraordinary point
		std::vector<double> extraordinaryPoints;
		for (int i = 0; i < vertices.size() / 3; i++) {
			// Q
			double qxx = 0, qyy = 0, qzz = 0;
			for (int j = 0; j < vertex2Edges[i].size(); j++) {
				int vId;
				if (edgesV[vertex2Edges[i][j]][0] != i)
					vId = edgesV[vertex2Edges[i][j]][0];
				else
					vId = edgesV[vertex2Edges[i][j]][1];
				qxx += vertices[3 * vId];
				qyy += vertices[3 * vId + 1];
				qzz += vertices[3 * vId + 2];
			}
			qxx /= vertex2Edges[i].size();
			qyy /= vertex2Edges[i].size();
			qzz /= vertex2Edges[i].size();
			// alpha
			double alpha = 0.625 - (0.375 + 0.25*cos(6.28 / vertex2Edges[i].size()))*(0.375 + 0.25*cos(6.28 / vertex2Edges[i].size()));
			double epxx, epyy, epzz;
			epxx = (1 - alpha)*vertices[3 * i] + alpha * qxx;
			epyy = (1 - alpha)*vertices[3 * i + 1] + alpha * qyy;
			epzz = (1 - alpha)*vertices[3 * i + 2] + alpha * qzz;
			extraordinaryPoints.push_back(epxx);
			extraordinaryPoints.push_back(epyy);
			extraordinaryPoints.push_back(epzz);
		}
		for (int i = 0; i < edgePoints.size(); i++)
			newVertices.push_back(edgePoints[i]);
		for (int i = 0; i < extraordinaryPoints.size(); i++)
			newVertices.push_back(extraordinaryPoints[i]);
		for (int i = 0; i < polygons.size(); i++) {
			newPolygons.push_back(std::vector<int>(0));
			auto j = polygon2Edges[i].begin();
			for (int k = 0; k < 3; k++) {
				newPolygons[newPolygons.size() - 1].push_back(*j);
				j++;
			}
			for (int k = 0; k < 3; k++) {
				newPolygons.push_back(std::vector<int>(0));
				int vId = polygons[i][k];
				int v0Id, v1Id;
				if (k == 0) {
					v0Id = polygons[i][1]; v1Id = polygons[i][2];
				}
				else if (k == 1) {
					v0Id = polygons[i][2]; v1Id = polygons[i][0];
				}
				else {
					v0Id = polygons[i][0]; v1Id = polygons[i][1];
				}
				std::vector<int>& ves = vertex2Edges[vId];
				int e0 = -1, e1 = -1;
				for (int kk = 0; kk < ves.size(); kk++) {
					if (edgesV[ves[kk]][0] == v0Id || edgesV[ves[kk]][1] == v0Id)
						e0 = ves[kk];
					if (edgesV[ves[kk]][0] == v1Id || edgesV[ves[kk]][1] == v1Id)
						e1 = ves[kk];
					if (e0 >= 0 && e1 >= 0)
						break;
				}
				newPolygons[newPolygons.size() - 1].push_back(vId + edges.size());
				newPolygons[newPolygons.size() - 1].push_back(e0);
				newPolygons[newPolygons.size() - 1].push_back(e1);
			}
		}
		vertices.resize(newVertices.size());
		for (int i = 0; i < vertices.size(); i++)
			vertices[i] = newVertices[i];
		polygons.clear();
		for (int i = 0; i < newPolygons.size(); i++) {
			polygons.push_back(std::vector<int>(0));
			for (int j = 0; j < newPolygons[i].size(); j++) {
				polygons[polygons.size() - 1].push_back(newPolygons[i][j]);
			}
		}
		newVertices.clear();
		newPolygons.clear();
		end_t = clock();
		totalTime = (double)(end_t - start_t) / CLOCKS_PER_SEC;
		std::cout << '\t' << "Subdivision times: " << time + 1 << '\t' << "Time cost: " << totalTime << 's' << std::endl;
	}
	std::cout << "Loop Surface Subdivision Success." << std::endl;
	return 1;
}

int Object::establishVertex2Polygons() {
	vertex2Polygons.clear();
	vertex2Polygons.resize(vertices.size() / 3);
	for (int i = 0; i < polygons.size(); i++) {
		for (int j = 0; j < polygons[i].size(); j++) {
			vertex2Polygons[polygons[i][j]].push_back(i);
		}
	}
	for (int i = 0; i < vertex2Polygons.size(); i++) {
		std::vector<int> sortedPolygons;
		std::vector<int>& ps = vertex2Polygons[i];
		int entry;
		for (int j = 0; j < polygons[ps[0]].size(); j++) {
			if (polygons[ps[0]][j] == i) {
				entry = j;
				break;
			}
		}
		if (entry == polygons[ps[0]].size() - 1)
			entry = polygons[ps[0]][0];
		else
			entry = polygons[ps[0]][entry + 1];
		sortedPolygons.push_back(ps[0]);
		while (sortedPolygons.size() < ps.size()) {
			for (int j = 0; j < ps.size(); j++) {
				if (ps[j] == sortedPolygons[sortedPolygons.size() - 1])
					continue;
				int flag = 0;
				for (int k = 0; k < polygons[ps[j]].size(); k++) {
					if (polygons[ps[j]][k] == entry) {
						sortedPolygons.push_back(ps[j]);
						int a, b, c, d;
						if (k == 0) {
							a = polygons[ps[j]][polygons[ps[j]].size() - 1];
							b = polygons[ps[j]][polygons[ps[j]].size() - 2];
						}
						else if (k == 1) {
							a = polygons[ps[j]][0];
							b = polygons[ps[j]][polygons[ps[j]].size() - 1];
						}
						else {
							a = polygons[ps[j]][k - 1];
							b = polygons[ps[j]][k - 2];
						}
						if (k == polygons[ps[j]].size() - 1) {
							c = polygons[ps[j]][0];
							d = polygons[ps[j]][1];
						}
						else if (k == polygons[ps[j]].size() - 2) {
							c = polygons[ps[j]][polygons[ps[j]].size() - 1];
							d = polygons[ps[j]][0];
						}
						else {
							c = polygons[ps[j]][k + 1];
							d = polygons[ps[j]][k + 2];
						}
						if (a == i)
							entry = b;
						else if (c == i)
							entry = d;
						flag++;
						break;
					}
				}
				if (flag)
					break;
			}
		}
		for (int j = 0; j < sortedPolygons.size(); j++) {
			ps[j] = sortedPolygons[j];
		}
	}
	return 1;
}

int Object::establishEdge2Polygons() {
	edge2Polygons.clear();
	edge2Polygons.resize(0);
	for (auto i = edges.begin(); i != edges.end(); i++) {
		edge2Polygons.push_back(std::vector<int>(0));
		auto j = i->begin();
		int v0 = *j; j++;
		int v1 = *j;
		std::vector<int>& polygonsOfV0 = vertex2Polygons[v0];
		std::vector<int>& polygonsOfV1 = vertex2Polygons[v1];
		for (int k = 0; k < polygonsOfV0.size(); k++) {
			auto position = std::find(polygonsOfV1.begin(), polygonsOfV1.end(), polygonsOfV0[k]);
			if (position != polygonsOfV1.end()) {
				edge2Polygons[edge2Polygons.size() - 1].push_back(*position);
			}
			if (edge2Polygons[edge2Polygons.size() - 1].size() >= 2)
				break;
		}
	}
	return 1;
}

int Object::establishEdges() {
	edges.clear();
	for (int i = 0; i < polygons.size(); i++) {
		for (int j = 0; j < polygons[i].size() - 1; j++) {
			edges.insert(std::set<int>({ polygons[i][j], polygons[i][j + 1] }));
		}
		edges.insert(std::set<int>({ polygons[i][polygons[i].size() - 1], polygons[i][0] }));
	}
	return 1;
}

int Object::establishEdgesV() {
	edgesV.clear();
	edgesV.resize(edges.size());
	int edgeId = 0;
	for (auto i = edges.begin(); i != edges.end(); i++) {
		auto temp = i->begin();
		int v0 = *temp; temp++;
		int v1 = *temp;
		edgesV[edgeId].push_back(v0);
		edgesV[edgeId].push_back(v1);
		edgeId++;
	}
	return 1;
}

int Object::establishVertex2Edges() {
	vertex2Edges.clear();
	vertex2Edges.resize(vertices.size() / 3);
	int edgeId = 0;
	for (auto i = edges.begin(); i != edges.end(); i++) {
		auto temp = i->begin();
		int v0 = *temp; temp++;
		int v1 = *temp;
		vertex2Edges[v0].push_back(edgeId);
		vertex2Edges[v1].push_back(edgeId);
		edgeId++;
	}
	return 1;
}

int Object::establishPolygon2Edges() {
	polygon2Edges.clear();
	polygon2Edges.resize(polygons.size());
	for (int i = 0; i < polygons.size(); i++) {
		for (int j = 0; j < polygons[i].size(); j++) {
			std::vector<int>& ves = vertex2Edges[polygons[i][j]];
			for (int k = 0; k < ves.size(); k++) {
				if ((edgesV[ves[k]][0] == polygons[i][0] && edgesV[ves[k]][1] == polygons[i][1]) ||
					(edgesV[ves[k]][0] == polygons[i][1] && edgesV[ves[k]][1] == polygons[i][0]) ||
					(edgesV[ves[k]][0] == polygons[i][1] && edgesV[ves[k]][1] == polygons[i][2]) ||
					(edgesV[ves[k]][0] == polygons[i][2] && edgesV[ves[k]][1] == polygons[i][1]) ||
					(edgesV[ves[k]][0] == polygons[i][2] && edgesV[ves[k]][1] == polygons[i][0]) ||
					(edgesV[ves[k]][0] == polygons[i][0] && edgesV[ves[k]][1] == polygons[i][2]))
					polygon2Edges[i].insert(ves[k]);
			}
		}
	}
	return 1;
}

int Object::writeObjFile(const char* filename) {
	std::ofstream fout(filename);
	for (int i = 0; i < vertices.size() / 3; i++) {
		fout << 'v' << ' ' << vertices[3 * i] << ' ' << vertices[3 * i + 1] << ' ' << vertices[3 * i + 2] << '\n';
	}
	for (int i = 0; i < polygons.size(); i++) {
		fout << 'f';
		for (int j = 0; j < polygons[i].size(); j++) {
			fout << ' ' << polygons[i][j] + 1;
		}
		fout << '\n';
	}
	return 1;
}

int Object::establishPolygonOffset() {
	polygonOffset.clear();
	polygonOffset.resize(polygons.size());
	int sum = 0;
	for (int i = 0; i < polygonOffset.size(); i++) {
		polygonOffset[i] = sum;
		sum += polygons[i].size();
	}
	return 1;
}
int Object::setColors(){
	colors.clear();
	colors.resize(polygons.size() * 3);
	for (int i = 0; i < polygons.size(); i++) {
		colors[3 * i] = rand() % 256;
		colors[3 * i + 1] = rand() % 256;
		colors[3 * i + 2] = rand() % 256;
	}
	return 1;
}
int Object::tri() {
	triangles.clear();
	for (int i = 0; i < polygons.size(); i++) {
		if (polygons[i].size() == 3) {
			triangles.push_back(polygons[i][0]);
			triangles.push_back(polygons[i][1]);
			triangles.push_back(polygons[i][2]);
			triColors.push_back(colors[3 * i]);
			triColors.push_back(colors[3 * i + 1]);
			triColors.push_back(colors[3 * i + 2]);
		}
		else {
			triangles.push_back(polygons[i][0]);
			triangles.push_back(polygons[i][1]);
			triangles.push_back(polygons[i][2]);
			triangles.push_back(polygons[i][0]);
			triangles.push_back(polygons[i][2]);
			triangles.push_back(polygons[i][3]);
			triColors.push_back(colors[3 * i]);
			triColors.push_back(colors[3 * i + 1]);
			triColors.push_back(colors[3 * i + 2]);
			triColors.push_back(colors[3 * i]);
			triColors.push_back(colors[3 * i + 1]);
			triColors.push_back(colors[3 * i + 2]);
		}
	}
	return 1;
}

int Object::establishVerticesUVNs(Vertex viewPoint, Vector viewDirection) {
	Vector up(0, 1, 0);
	if (abs(abs(up * viewDirection) - up.magnitude()*viewDirection.magnitude()) <= 1.0e-6) {
		up.setY(0);
		up.setZ(1);
	}
	Vector N = viewDirection; N.normalize();
	Vector V; V = viewDirection ^ up; V.normalize();
	Vector U; U = V ^ viewDirection; U.normalize();
	verticesUVNs.resize(vertices.size());
	for (int i = 0; i < vertices.size() / 3; i++) {
		verticesUVNs[3 * i] = U.x * vertices[3 * i] + U.y * vertices[3 * i + 1] + U.z *vertices[3 * i + 2]
			- viewPoint.x() * U.x - viewPoint.y() * U.y - viewPoint.z() * U.z;
		verticesUVNs[3 * i + 1] = V.x * vertices[3 * i] + V.y * vertices[3 * i + 1] + V.z*vertices[3 * i + 2]
			- viewPoint.x() * V.x - viewPoint.y() * V.y - viewPoint.z() * V.z;
		verticesUVNs[3 * i + 2] = N.x * vertices[3 * i] + N.y * vertices[3 * i + 1] + N.z*vertices[3 * i + 2]
			- viewPoint.x() * N.x - viewPoint.y() * N.y - viewPoint.z() * N.z;
	}
	return 1;
}

int Object::clearFrameBuffer() {
	for (int i = 0; i < LENGTH*HEIGHT * 3; i++)
		frameBuffer[i] = 0;
	return 1;
}
int Object::clearZBuffer() {
	for (int i = 0; i < LENGTH*HEIGHT; i++)
		zBuffer[i] = LARGEST_DISTANCE;
	return 1;
}

int Object::establishVerticesUPVPs(double distance) {
	verticesUPVPs.resize(verticesUVNs.size() * 2 / 3);
	double min = LENGTH / 2;
	for (int i = 0; i < verticesUVNs.size() / 3; i++) {
		verticesUPVPs[2 * i] = verticesUVNs[3 * i] * verticesUVNs[3 * i + 2] / distance;
		verticesUPVPs[2 * i + 1] = verticesUVNs[3 * i + 1] * verticesUVNs[3 * i + 2] / distance;
	}
	return 1;
}

int Object::ZBuffer(Vertex viewPoint, Vector viewDirection) {
	std::cout << "ZBuffer visible surface determinating, please wait..." << std::endl;
	clearFrameBuffer();
	clearZBuffer();
	clock_t start_t, finish_t;
	double totalTime = -1;
	start_t = clock();
	establishVerticesUVNs(viewPoint, viewDirection);
	establishVerticesUPVPs(VIEWPORT_DISTANCE);
	for (int i = 0; i < triangles.size() / 3; i++) {
		// ��ÿһ������Σ��ȼ��������ӵ������е�����;
		int vIds[3];
		vIds[0] = triangles[3 * i];
		vIds[1] = triangles[3 * i + 1];
		vIds[2] = triangles[3 * i + 2];
		std::vector<double> coordUVN(FACET * 3);
		for (int j = 0; j < FACET; j++) {
			coordUVN[3 * j] = verticesUVNs[3 * vIds[j]];
			coordUVN[3 * j + 1] = verticesUVNs[3 * vIds[j] + 1];
			coordUVN[3 * j + 2] = verticesUVNs[3 * vIds[j] + 2];
		}
		// Ȼ����ͶӰ���Ӵ�ƽ�棬�ж��Ƿ��������Ӵ��ڲ�
		std::vector<double> coordUPVP(FACET * 2);
		for (int j = 0; j < FACET; j++) {
			coordUPVP[2 * j] = verticesUPVPs[2 * vIds[j]];
			coordUPVP[2 * j + 1] = verticesUPVPs[2 * vIds[j] + 1];
		}
		// �����������Ӵ��ڲ������������ֵ���ж��Ƿ�д��zbuffer
		int flag = 0;
		for (int j = 0; j < FACET; j++) {
			if (coordUPVP[2 * j] >= -LENGTH / 2 && coordUPVP[2 * j] < LENGTH / 2 &&
				coordUPVP[2 * j + 1] >= -HEIGHT / 2 && coordUPVP[2 * j + 1] < HEIGHT / 2) {
				flag++;
				break;
			}
		}
		if (flag) {// ��������Ӵ��ڲ����������ֵ
			double maxY = std::max({ coordUPVP[1], coordUPVP[3], coordUPVP[5] });
			double minY = std::min({ coordUPVP[1], coordUPVP[3], coordUPVP[5] });
			double maxX = std::max({ coordUPVP[0], coordUPVP[2], coordUPVP[4] });
			double minX = std::min({ coordUPVP[0], coordUPVP[2], coordUPVP[4] });
			if (maxX > LENGTH / 2 - 1)
				maxX = LENGTH / 2 - 1;
			if (maxY > HEIGHT / 2 - 1)
				maxY = HEIGHT / 2 - 1;
			if (minY < -HEIGHT / 2)
				minY = -HEIGHT / 2;
			if (minX < -LENGTH / 2)
				minX = -LENGTH / 2;
			// ���������ƽ�淽�̣���������������Ƭ
			double A = coordUVN[1] * (coordUVN[1 * 3 + 2] - coordUVN[2 * 3 + 2])
				+ coordUVN[1 * 3 + 1] * (coordUVN[2 * 3 + 2] - coordUVN[2])
				+ coordUVN[2 * 3 + 1] * (coordUVN[2] - coordUVN[1 * 3 + 2]);
			double B = coordUVN[2] * (coordUVN[3 * 1] - coordUVN[3 * 2])
				+ coordUVN[1 * 3 + 2] * (coordUVN[3 * 2] - coordUVN[0])
				+ coordUVN[2 * 3 + 2] * (coordUVN[0] - coordUVN[3 * 1]);
			double C = coordUVN[0] * (coordUVN[3 * 1 + 1] - coordUVN[3 * 2 + 1])
				+ coordUVN[1 * 3] * (coordUVN[3 * 2 + 1] - coordUVN[1])
				+ coordUVN[2 * 3] * (coordUVN[1] - coordUVN[3 * 1 + 1]);
			double D = -coordUVN[0] * (coordUVN[3 * 1 + 1] * coordUVN[3 * 2 + 2] - coordUVN[3 * 2 + 1] * coordUVN[3 * 1 + 2])
				- coordUVN[3 * 1] * (coordUVN[3 * 2 + 1] * coordUVN[2] - coordUVN[1] * coordUVN[3 * 2 + 2])
				- coordUVN[3 * 2] * (coordUVN[1] * coordUVN[3 * 1 + 2] - coordUVN[3 * 1 + 1] * coordUVN[2]);
			for (int j = maxY; j > minY; j--) {
				for (int k = minX; k < maxX; k++) {
					int flagss = 0;
					// ���ж����ص��Ƿ�����Ƭ�ڲ�
					// ���ж�������������������Ƭ
					Vertex s(k, j);
					for (int ii = 0; ii < FACET; ii++) {
						Vertex p(coordUPVP[2 * ii], coordUPVP[2 * ii + 1]);
						Vertex q;
						Vertex ss;
						if (ii == 0) {
							q.coordinate = std::vector<double>({ coordUPVP[2 * 1], coordUPVP[2 * 1 + 1], 0 });
							ss.coordinate = std::vector<double>({ coordUPVP[2 * 2], coordUPVP[2 * 2 + 1], 0 });
						}
						else if (ii == 1) {
							q.coordinate = std::vector<double>({ coordUPVP[2 * 2], coordUPVP[2 * 2 + 1], 0 });
							ss.coordinate = std::vector<double>({ coordUPVP[0], coordUPVP[1], 0 });
						}
						else if (ii == 2) {
							q.coordinate = std::vector<double>({ coordUPVP[0], coordUPVP[1], 0 });
							ss.coordinate = std::vector<double>({ coordUPVP[2 * 1], coordUPVP[2 * 1 + 1], 0 });
						}
						int flags = toLeftTest(p, q, s);
						if (flags == 0) {
							flagss = -1;
							break;
						}
						else if (flags == toLeftTest(p, q, ss))
							continue;
						else {
							flagss++;
							break;
						}
					}
					// ������������Ƭ�ڲ�
					// debug
					//if(1){
					if (!flagss) {
						double nn = -D / (A*((double)k) / VIEWPORT_DISTANCE + B * ((double)j) / VIEWPORT_DISTANCE + C);
						int pixelId = LENGTH * (HEIGHT / 2 - 1 - j) + k + LENGTH / 2;
						if (nn < zBuffer[pixelId]) {
							zBuffer[pixelId] = nn;
							// ��������صĹ���ֵ���Բ�д��frameBuffer
							if (flagss == -1) {
								frameBuffer[pixelId * 3] = 0;
								frameBuffer[pixelId * 3 + 1] = 255;
								frameBuffer[pixelId * 3 + 2] = 0;
							}
							else {
								frameBuffer[pixelId * 3] = triColors[i * 3];
								frameBuffer[pixelId * 3 + 1] = triColors[i * 3 + 1];
								frameBuffer[pixelId * 3 + 2] = triColors[i * 3 + 2];
							}
						}
					}
				}
			}
		}
	}
	finish_t = clock();
	totalTime = (double)(finish_t - start_t) / CLOCKS_PER_SEC;
	std::cout << "ZBuffer visible surface determination Success!" << std::endl;
	std::cout << "ZBuffer visible surface determination Time cost: " << totalTime << 's' << std::endl;
	return 0;
}

int Object::writePpmFile(const char* filename) {
	std::cout << "Vsible surface determination result writing to picture..." << std::endl;
	int w = LENGTH; int h = HEIGHT;
	std::ofstream fout(filename);
	fout << "P3" << '\n';
	fout << "# " << filename << '\n';
	fout << w << ' ' << h << '\n';
	fout << "255" << '\n';
	for (int i = 0; i < h; i++) {
		for (int j = 0; j < w; j++) {
			fout << frameBuffer[3 * (i*w + j)] << ' ' << frameBuffer[3 * (i*w + j) + 1] << ' ' << frameBuffer[3 * (i*w + j) + 2] << ' ';
		}
		fout << '\n';
	}
	std::cout << "Picture file writing Success!" << std::endl;
	return 1;
}