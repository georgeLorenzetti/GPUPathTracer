#pragma once
#include <precomp.h>

using namespace glm;
/////////////////////
/////BVH NODE //////
////////////////////

void BVHNode::CalculateBounds(std::vector<AABB> AABBs) {
	if (AABBs.size() <= 0) {
		return;
	}
	AABB aabb;
	aabb = AABBs[0];

	for (int i = 0; i < AABBs.size(); i++) {
		if (aabb.bounds[1].x < AABBs[i].bounds[1].x) {
			aabb.bounds[1].x = AABBs[i].bounds[1].x;
		}

		if (aabb.bounds[0].x > AABBs[i].bounds[0].x) {
			aabb.bounds[0].x = AABBs[i].bounds[0].x;
		}

		if (aabb.bounds[1].y < AABBs[i].bounds[1].y) {
			aabb.bounds[1].y = AABBs[i].bounds[1].y;
		}

		if (aabb.bounds[0].y > AABBs[i].bounds[0].y) {
			aabb.bounds[0].y = AABBs[i].bounds[0].y;
		}

		if (aabb.bounds[1].z < AABBs[i].bounds[1].z) {
			aabb.bounds[1].z = AABBs[i].bounds[1].z;
		}

		if (aabb.bounds[0].z > AABBs[i].bounds[0].z) {
			aabb.bounds[0].z = AABBs[i].bounds[0].z;
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

void BVH::ConstructBVH(std::vector<vec3> t_vertices, std::vector<int> t_indices, int tri_count) {

	std::vector<AABB> AABBs = std::vector<AABB>();
	//indices = new int[tri_count];

	for (int i = 0; i < tri_count; i++) {
		//indices[i] = i;

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
		bounds.index = i;
		AABBs.push_back(bounds);
	}

	//create the root node of the tree
	this->root_node = new BVHNode();
	this->root_node->CalculateBounds(AABBs);

	//subdivide to create the child nodes recursively
	Subdivide(root_node, AABBs, 1);
}

void BVH::Subdivide(BVHNode* node, std::vector<AABB> AABBs, bool side) {
	//base case, this is a leaf node if there are less than 3 objects inside
	if (AABBs.size() < 4) {
		node->is_leaf = true;
		for (int i = 0; i < AABBs.size(); i++) {
			node->leaf_triangles.push_back(AABBs[i].index);
		}
		return;
	}

	node->left = new BVHNode();
	node->right = new BVHNode();

	//calculate partition
	float split_position = -1;
	int split = Partition(node, AABBs, split_position);

	//if no split was made because child nodes do not improve cost then make this a leaf
	if (split == -1) {
		node->is_leaf = true;
		for (int i = 0; i < AABBs.size(); i++) {
			node->leaf_triangles.push_back(AABBs[i].index);
		}
		return;
	}
	
	//else subdivide
	std::vector<AABB> left_AABBs = std::vector<AABB>();
	std::vector<AABB> right_AABBs = std::vector<AABB>();

	float centre_value;
	if (split == 0) {
		centre_value = node->aabb.centre_point.x;
	}
	else if (split == 1) {
		centre_value = node->aabb.centre_point.y;
	}
	else {
		centre_value = node->aabb.centre_point.z;
	}

	for (int i = 0; i < AABBs.size(); i++) {
		float centre_check;
		if (split == 0) {
			centre_check = AABBs[i].centre_point.x;
		}
		else if (split == 1) {
			centre_check = AABBs[i].centre_point.y;
		}
		else {
			centre_check = AABBs[i].centre_point.z;
		}

		if (centre_check <= centre_value) {
			left_AABBs.push_back(AABBs[i]);
		}
		else {
			right_AABBs.push_back(AABBs[i]);
		}
	}


	node->left->CalculateBounds(left_AABBs);
	node->right->CalculateBounds(right_AABBs);

	if (left_AABBs.size() == 0 || right_AABBs.size() == 0) {
		node->is_leaf = true;
		for (int i = 0; i < AABBs.size(); i++) {
			node->leaf_triangles.push_back(AABBs[i].index);
		}
		return;
	}
	//recursive calls to partition child nodes
	Subdivide(node->left, left_AABBs, 1);
	Subdivide(node->right, right_AABBs, 0);

	node->is_leaf = false;

}

int BVH::Partition(BVHNode* node, std::vector<AABB> aabbs, float& split_position) {
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

	std::vector<float> xVals = std::vector<float>();
	std::vector<float> yVals = std::vector<float>();
	std::vector<float> zVals = std::vector<float>();
	xVals.reserve(aabbs.size());
	yVals.reserve(aabbs.size());
	zVals.reserve(aabbs.size());

	for (int i = 0; i < aabbs.size(); i++) {
		xVals[i] = aabbs[i].centre_point.x;
		yVals[i] = aabbs[i].centre_point.y;
		zVals[i] = aabbs[i].centre_point.z;
	}

	std::sort(xVals.begin(), xVals.end());
	std::sort(yVals.begin(), yVals.end());
	std::sort(zVals.begin(), zVals.end());


	xObjectMedian = xVals[aabbs.size() / 2];
	yObjectMedian = yVals[aabbs.size() / 2];
	zObjectMedian = zVals[aabbs.size() / 2];

	//choose 9 equally spaced possibilities on each axis and check cost.
	float height = node->aabb.bounds[1].y - node->aabb.bounds[0].y;
	float width = node->aabb.bounds[1].x - node->aabb.bounds[0].x;
	float depth = node->aabb.bounds[1].z - node->aabb.bounds[0].z;

	float SA = 2 * (height * width) + 2 * (height * depth) + 2 * (width * depth);
	float parentCost = SA * aabbs.size();
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

		tempCost = (SALeft * (float)n) + (SARight * ((float)aabbs.size() - (float)n));
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

		tempCost = (SALeft * (float)n) + (SARight * ((float)aabbs.size() - (float)n));
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

		tempCost = (SALeft * (float)n) + (SARight * ((float)aabbs.size() - (float)n));
		if (tempCost < splitCost) {
			splitCost = tempCost;
			splitAxis = 'z';
			splitPosition = split;
		}
	}

	//if best cost is the same as the parent cost then don't split and just have a leaf node
	if (splitCost == parentCost) {
		return -1;
	}

	//return the first index for the right child as the point of partition
	split_position = splitPosition;

	switch (splitAxis) {
	case 'x':
		return 0;
		break;
	case 'y':
		return 1;
		break;
	case 'z':
		return 2;
		break;
	}
	return 0;
}

int count_nodes(BVHNode* root) {
	if (!root->is_leaf) {
		return 1 + count_nodes(root->left) + count_nodes(root->right);
	}

	return 1;
}

void BVH::PopulateCFBVH(unsigned int& cumulative_index, unsigned int& triangle_index, BVHNode* node) {
	int current_index = cumulative_index;
	this->cf_bvh[current_index].bounds[0] = node->aabb.bounds[0];
	this->cf_bvh[current_index].bounds[1] = node->aabb.bounds[1];

	if (!(node->is_leaf)) {
		//std::cout << "INNER " << cumulative_index << std::endl;
		cumulative_index++;
		unsigned int left_index = cumulative_index;
		PopulateCFBVH(cumulative_index, triangle_index, node->left);

		cumulative_index++;
		unsigned int right_index = cumulative_index;
		PopulateCFBVH(cumulative_index, triangle_index, node->right);

		cf_bvh[current_index].u.inner.index_left = left_index;
		cf_bvh[current_index].u.inner.index_right = right_index;
	}
	else {
		unsigned int count = (unsigned int)node->leaf_triangles.size();
		cf_bvh[current_index].u.leaf.count = 0x80000000 | count; //basically sets the first bit of count a
		cf_bvh[current_index].u.leaf.index_first_tri = triangle_index;
		//std::cout << "LEAF " << cumulative_index << std::endl;
		for (unsigned int i = 0; i < count; i++) {
			this->triangle_indices[triangle_index] = node->leaf_triangles[i];
			triangle_index++;
		}
	}
}

void BVH::ConstructCacheFriendly(int tri_count) {
	int node_count = count_nodes(this->root_node);

	this->cf_bvh = new BVHNode_CacheFriendly[node_count];
	this->triangle_indices = new int[tri_count];

	unsigned int cumulative_index = 0;
	unsigned int triangle_index = 0;
	PopulateCFBVH(cumulative_index, triangle_index, this->root_node);

	cudaAssert(Malloc(&(cf_bvh_gpu), node_count * sizeof(BVHNode_CacheFriendly)));
	cudaAssert(Malloc(&(triangle_indices_gpu), tri_count * sizeof(int)));
	cudaAssert(Memcpy(cf_bvh_gpu, cf_bvh, node_count * sizeof(BVHNode_CacheFriendly), cudaMemcpyHostToDevice));
	cudaAssert(Memcpy(triangle_indices_gpu, triangle_indices, tri_count * sizeof(int), cudaMemcpyHostToDevice));
}