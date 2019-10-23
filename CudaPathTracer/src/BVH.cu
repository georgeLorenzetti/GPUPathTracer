#pragma once
#include <precomp.h>

using namespace glm;
/////////////////////
/////BVH NODE //////
////////////////////

void BVHNode::CalculateBounds(std::vector<AABB> t_AABBs, int first, int count, int* indices, int tri_count) {
	if (tri_count <= 0) {
		return;
	}
	AABB aabb;
	aabb = t_AABBs[indices[first]];

	for (int i = first + 1; i < (first + count); i++) {
		if (aabb.bounds[1].x < t_AABBs[indices[i]].bounds[1].x) {
			aabb.bounds[1].x = t_AABBs[indices[i]].bounds[1].x;
		}

		if (aabb.bounds[0].x > t_AABBs[indices[i]].bounds[0].x) {
			aabb.bounds[0].x = t_AABBs[indices[i]].bounds[0].x;
		}

		if (aabb.bounds[1].y < t_AABBs[indices[i]].bounds[1].y) {
			aabb.bounds[1].y = t_AABBs[indices[i]].bounds[1].y;
		}

		if (aabb.bounds[0].y > t_AABBs[indices[i]].bounds[0].y) {
			aabb.bounds[0].y = t_AABBs[indices[i]].bounds[0].y;
		}

		if (aabb.bounds[1].z < t_AABBs[indices[i]].bounds[1].z) {
			aabb.bounds[1].z = t_AABBs[indices[i]].bounds[1].z;
		}

		if (aabb.bounds[0].z > t_AABBs[indices[i]].bounds[0].z) {
			aabb.bounds[0].z = t_AABBs[indices[i]].bounds[0].z;
		}
	}

	aabb.CalculateCentre();
	this->aabb = aabb;

	//printf("%f %f %f | %f %f %f \n", aabb.bounds[0].x, aabb.bounds[0].y, aabb.bounds[0].z, aabb.bounds[1].x, aabb.bounds[1].y, aabb.bounds[1].z);
}

/////////////
////AABB/////
/////////////
void AABB::CalculateCentre() {
	this->centre_point = vec3((bounds[0].x + bounds[1].x) / 2.0f, (bounds[0].y + bounds[1].y) / 2.0f, (bounds[0].z + bounds[1].z) / 2.0f);
}

///////////////
/////BVH //////
//////////////
BVH::BVH() {}

BVH::BVH(vec3* t_vertices_gpu, vec3* t_normals_gpu, int* t_indices_gpu) {
	this->t_vertices_gpu = t_vertices_gpu;
	this->t_normals_gpu = t_normals_gpu;
	this->t_indices_gpu = t_indices_gpu;
}

void BVH::ConstructBVH(std::vector<vec3> t_vertices, std::vector<int> t_indices, int tri_count) {
	this->t_AABBs = std::vector<AABB>();
	indices = new int[tri_count];


	for (int i = 0; i < tri_count; i++) {
		indices[i] = i;

		vec3 v1 = t_vertices[t_indices[i * 3 + 1]];
		vec3 v2 = t_vertices[t_indices[i * 3 + 2]];
		vec3 v3 = t_vertices[t_indices[i * 3]];

		AABB bounds = AABB();
		bounds.bounds[1].x = max(v1.x, max(v2.x, v3.x));
		bounds.bounds[0].x = min(v1.x, min(v2.x, v3.x));
		bounds.bounds[1].y = max(v1.y, max(v2.y, v3.y));
		bounds.bounds[0].y = min(v1.y, min(v2.y, v3.y));
		bounds.bounds[1].z = max(v1.z, max(v2.z, v3.z));
		bounds.bounds[0].z = min(v1.z, min(v2.z, v3.z));
		bounds.CalculateCentre();

		t_AABBs.push_back(bounds);
	}

	//create the root node of the tree
	this->root_node = new BVHNode();
	this->root_node->first = 0;
	this->root_node->count = tri_count;
	this->root_node->CalculateBounds(t_AABBs, this->root_node->first, this->root_node->count, indices, tri_count);

	//subdivide to create the child nodes recursively
	Subdivide(root_node, t_vertices, t_indices, tri_count);
}

void BVH::Subdivide(BVHNode* node, std::vector<vec3> t_vertices, std::vector<int> t_indices, int tri_count) {
	//base case, this is a leaf node if there are less than 3 objects inside
	if (node->count < 3) {
		node->is_leaf = true;
		for (int i = node->first; i < node->first + node->count; i++) {
			std::cout << i << " ";
		}
		std::cout << std::endl;
		return;
	}

	node->left = new BVHNode();
	node->right = new BVHNode();

	//calculate partition
	int split = Partition(node);
	//if no split was made because child nodes do not improve cost then make this a leaf
	if (split == -1) {
		node->is_leaf = true;
		for (int i = node->first; i < node->first + node->count; i++) {
			std::cout << i << " ";
		}
		std::cout << std::endl;
		return;
	}
	if (split == node->first) {
		return;
	}
	//else subdivide
	node->left->first = node->first;
	node->left->count = split - node->first;
	node->left->CalculateBounds(t_AABBs, node->left->first, node->left->count, indices, tri_count);

	node->right->first = split;
	node->right->count = node->count - node->left->count;
	node->right->CalculateBounds(t_AABBs, node->right->first, node->right->count, indices, tri_count);

	//recursive calls to partition child nodes
	if (node->left->count > 0) {
		Subdivide(node->left, t_vertices, t_indices, tri_count);
	}
	else {
		node->left->is_leaf = true;
	}
	if (node->right->count > 0) {
		Subdivide(node->right, t_vertices, t_indices, tri_count);
	}
	else {
		node->right->is_leaf = true;
	}

	node->is_leaf = false;

}

int BVH::Partition(BVHNode* node) {
	/*--FINDING SPLIT AXIS--*/
	//Using the Arbitrary Acyclic Algorithm from: http://graphicsinterface.org/wp-content/uploads/gi1989-22.pdf
	//spatial medians
	float xSpatialMedian = (node->aabb.bounds[0].x + node->aabb.bounds[1].x) / 2.0f;
	float ySpatialMedian = (node->aabb.bounds[0].y + node->aabb.bounds[1].y) / 2.0f;
	float zSpatialMedian = (node->aabb.bounds[0].z + node->aabb.bounds[1].z) / 2.0f;

	//object medians
	float xObjectMedian;
	float yObjectMedian;
	float zObjectMedian;

	float* xVals = new float[node->count];
	float* yVals = new float[node->count];
	float* zVals = new float[node->count];

	for (int i = node->first; i < (node->first + node->count); i++) {
		xVals[i - node->first] = t_AABBs[indices[i]].centre_point.x;
		yVals[i - node->first] = t_AABBs[indices[i]].centre_point.y;
		zVals[i - node->first] = t_AABBs[indices[i]].centre_point.z;
	}

	std::sort(xVals, xVals + node->count);
	std::sort(yVals, yVals + node->count);
	std::sort(zVals, zVals + node->count);

	xObjectMedian = xVals[node->count / 2];
	yObjectMedian = yVals[node->count / 2];
	zObjectMedian = zVals[node->count / 2];

	//choose 9 equally spaced possibilities on each axis and check cost.
	float height = node->aabb.bounds[1].y - node->aabb.bounds[0].y;
	float width = node->aabb.bounds[1].x - node->aabb.bounds[0].x;
	float depth = node->aabb.bounds[1].z - node->aabb.bounds[0].z;

	float SA = 2 * (height * width) + 2 * (height * depth) + 2 * (width * depth);
	float parentCost = SA * node->count;
	float splitCost = parentCost;
	char splitAxis = 'x';
	float splitPosition = xSpatialMedian;
	//x axis
	float min = fminf(xObjectMedian, xSpatialMedian);
	float max = fmaxf(xObjectMedian, xSpatialMedian);
	float step = (max - min) / 9.0f;
	for (int i = 0; i < 9; i++) {
		float split = (min + (step * i));
		float tempCost;

		float widthL = split - node->aabb.bounds[0].x;
		float widthR = node->aabb.bounds[1].x - split;

		float SALeft = 2 * (height * widthL) + 2 * (height * depth) + 2 * (widthL * depth);
		float SARight = 2 * (height * widthR) + 2 * (height * depth) + 2 * (widthR * depth);

		int n = 0;
		while (xVals[n] < split) {
			n++;
		}

		tempCost = (SALeft * (float)n) + (SARight * ((float)node->count - (float)n));
		if (tempCost < splitCost) {
			splitCost = tempCost;
			splitAxis = 'x';
			splitPosition = split;
		}
	}

	//y axis
	min = fminf(yObjectMedian, ySpatialMedian);
	max = fmaxf(yObjectMedian, ySpatialMedian);
	step = (max - min) / 9.0f;
	for (int i = 0; i < 9; i++) {
		float split = (min + (step * i));
		float tempCost;

		float heightL = split - node->aabb.bounds[0].y;
		float heightR = node->aabb.bounds[1].y - split;

		float SALeft = 2 * (heightL * width) + 2 * (heightL * depth) + 2 * (width * depth);
		float SARight = 2 * (heightR * width) + 2 * (heightR * depth) + 2 * (width * depth);

		int n = 0;
		while (yVals[n] < split) {
			n++;
		}

		tempCost = (SALeft * (float)n) + (SARight * ((float)node->count - (float)n));
		if (tempCost < splitCost) {
			splitCost = tempCost;
			splitAxis = 'y';
			splitPosition = split;
		}
	}

	//z axis
	min = fminf(zObjectMedian, zSpatialMedian);
	max = fmaxf(zObjectMedian, zSpatialMedian);
	step = (max - min) / 9.0f;
	for (int i = 0; i < 9; i++) {
		float split = (min + (step * i));
		float tempCost;

		float depthL = split - node->aabb.bounds[0].z;
		float depthR = node->aabb.bounds[1].z - split;

		float SALeft = 2 * (height * width) + 2 * (height * depthL) + 2 * (width * depthL);
		float SARight = 2 * (height * width) + 2 * (height * depthR) + 2 * (width * depthR);

		int n = 0;
		while (zVals[n] < split) {
			n++;
		}

		tempCost = (SALeft * (float)n) + (SARight * ((float)node->count - (float)n));
		if (tempCost < splitCost) {
			splitCost = tempCost;
			splitAxis = 'z';
			splitPosition = split;
		}
	}

	//if best cost is the same as the parent cost then don't split and just have a leaf node
	if (splitCost == parentCost){
		return -1;
	}
	std::vector<int> lessEqual;
	std::vector<int> more;
	//sort into 2 arrays of object indices on the left and right of partition
	for (int i = node->first; i < (node->first + node->count); i++) {
		switch (splitAxis) {
		case 'x':
			if (t_AABBs[indices[i]].centre_point.x <= splitPosition) {
				lessEqual.push_back(indices[i]);
			}
			else {
				more.push_back(indices[i]);
			}
			break;
		case 'y':
			if (t_AABBs[indices[i]].centre_point.y <= splitPosition) {
				lessEqual.push_back(indices[i]);
			}
			else {
				more.push_back(indices[i]);
			}
			break;
		case 'z':
			if (t_AABBs[indices[i]].centre_point.z <= splitPosition) {
				lessEqual.push_back(indices[i]);
			}
			else {
				more.push_back(indices[i]);
			}
			break;
		}
	}

	//if next partition is same as current partition
	if (lessEqual.size() == 0 || more.size() == 0) {
		return -1;
	}

	//rearrange indices according to partition
	int partitionPoint = node->first + node->count - 1;
	for (int i = 0; i < lessEqual.size(); i++) {
		indices[node->first + i] = lessEqual[i];
	}
	for (int i = 0; i < more.size(); i++) {
		indices[node->first + lessEqual.size() + i] = more[i];
	}
	//return the first index for the right child as the point of partition
	return node->first + lessEqual.size();
}