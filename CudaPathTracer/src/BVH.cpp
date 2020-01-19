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

	//add an epsilon because my intersection code is trash and doesnt count intersections when absolute values are the same
	if (abs(aabb.bounds[0].x) == abs(aabb.bounds[1].x)) {
		aabb.bounds[0].x -= 0.000001f;
	}
	if (abs(aabb.bounds[0].y) == abs(aabb.bounds[1].y)) {
		aabb.bounds[0].y -= 0.000001f;
	}
	if (abs(aabb.bounds[0].z) == abs(aabb.bounds[1].z)) {
		aabb.bounds[0].z -= 0.000001f;
	}
	aabb.CalculateCentre();
	this->aabb = aabb;
}

/////////////
//MBVHNode//
////////////
void MBVHNode::CalculateBounds() {
	if (this->children.size() <= 0) {
		return;
	}
	AABB aabb = AABB();
	for (int i = 0; i < this->children.size(); i++) {

		if (aabb.bounds[1].x < this->children[i]->aabb.bounds[1].x) {
			aabb.bounds[1].x = this->children[i]->aabb.bounds[1].x;
		}

		if (aabb.bounds[0].x > this->children[i]->aabb.bounds[0].x) {
			aabb.bounds[0].x = this->children[i]->aabb.bounds[0].x;
		}

		if (aabb.bounds[1].y < this->children[i]->aabb.bounds[1].y) {
			aabb.bounds[1].y = this->children[i]->aabb.bounds[1].y;
		}

		if (aabb.bounds[0].y > this->children[i]->aabb.bounds[0].y) {
			aabb.bounds[0].y = this->children[i]->aabb.bounds[0].y;
		}

		if (aabb.bounds[1].z < this->children[i]->aabb.bounds[1].z) {
			aabb.bounds[1].z = this->children[i]->aabb.bounds[1].z;
		}

		if (aabb.bounds[0].z > this->children[i]->aabb.bounds[0].z) {
			aabb.bounds[0].z = this->children[i]->aabb.bounds[0].z;
		}
	}

	aabb.CalculateCentre();
	this->aabb = aabb;
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

void BVH::ConstructBVH(std::vector<vec3>& t_vertices, std::vector<int>& t_indices, int tri_count) {

	std::vector<AABB> AABBs = std::vector<AABB>();
	//indices = new int[tri_count];

	for (int i = 0; i < tri_count; i++) {
		//indices[i] = i;

		vec3 v1 = t_vertices[t_indices[i * 3 + 1]];
		vec3 v2 = t_vertices[t_indices[i * 3 + 2]];
		vec3 v3 = t_vertices[t_indices[i * 3]];

		AABB bounds = AABB();
		bounds.bounds[1].x = glm::max(v1.x, glm::max(v2.x, v3.x));
		bounds.bounds[0].x = glm::min(v1.x, glm::min(v2.x, v3.x));
		bounds.bounds[1].y = glm::max(v1.y, glm::max(v2.y, v3.y));
		bounds.bounds[0].y = glm::min(v1.y, glm::min(v2.y, v3.y));
		bounds.bounds[1].z = glm::max(v1.z, glm::max(v2.z, v3.z));
		bounds.bounds[0].z = glm::min(v1.z, glm::min(v2.z, v3.z));
		bounds.CalculateCentre();
		bounds.index = i;
		AABBs.push_back(bounds);
	}

	//create the root node of the tree
	this->root_node = new BVHNode();
	this->root_node->CalculateBounds(AABBs);

	//subdivide to create the child nodes recursively
	Subdivide(root_node, AABBs, 1, t_vertices, t_indices);
}

void BVH::Subdivide(BVHNode* node, std::vector<AABB> AABBs, bool side, std::vector<glm::vec3>& t_vertices, std::vector<int>& t_indices) {
	//base case, this is a leaf node if there are less than 3 objects inside
	if (AABBs.size() < 4) {
		node->is_leaf = true;
		for (int i = 0; i < AABBs.size(); i++) {
			node->leaf_triangles.push_back(AABBs[i].index);
		}
		node->CalculateBounds(AABBs);
		return;
	}

	node->left = new BVHNode();
	node->right = new BVHNode();

	//calculate partition
	float split_cost_object = FLT_MAX;
	float split_position_object = -1.0f;
	int split_axis_object = ObjectPartition(node, AABBs, split_position_object, split_cost_object);

	float split_cost_spatial = FLT_MAX;
	float split_position_spatial = -1.0f;
	//int split_axis_spatial = SpatialPartition(node, AABBs, split_position_spatial, split_cost_spatial, t_vertices, t_indices);

	//if no split was made because child nodes do not improve cost then make this a leaf
	if (split_axis_object == -1) {
		node->is_leaf = true;
		for (int i = 0; i < AABBs.size(); i++) {
			node->leaf_triangles.push_back(AABBs[i].index);
		}
		node->CalculateBounds(AABBs);
		return;
	}
	//else subdivide
	std::vector<AABB> left_AABBs = std::vector<AABB>();
	std::vector<AABB> right_AABBs = std::vector<AABB>();

	MakeObjectSplit(left_AABBs, right_AABBs, split_axis_object, split_position_object, node, AABBs);


	if (left_AABBs.size() == 0 || right_AABBs.size() == 0) {
		node->is_leaf = true;
		for (int i = 0; i < AABBs.size(); i++) {
			node->leaf_triangles.push_back(AABBs[i].index);
		}
		node->CalculateBounds(AABBs);
		return;
	}
	//recursive calls to partition child nodes
	Subdivide(node->left, left_AABBs, 1, t_vertices, t_indices);
	Subdivide(node->right, right_AABBs, 0, t_vertices, t_indices);

	node->is_leaf = false;
}

int BVH::ObjectPartition(BVHNode* node, std::vector<AABB> AABBs, float& split_position, float& split_cost) {
	/*--FINDING SPLIT AXIS--*/
	//Using the Arbitrary Acyclic Algorithm from: http://graphicsinterface.org/wp-content/uploads/gi1989-22.pdf
	//spatial medians
	float x_spatial_median = (node->aabb.bounds[0].x + node->aabb.bounds[1].x) / 2.0f;
	float y_spatial_median = (node->aabb.bounds[0].y + node->aabb.bounds[1].y) / 2.0f;
	float z_spatial_median = (node->aabb.bounds[0].z + node->aabb.bounds[1].z) / 2.0f;

	//object medians
	float x_object_median;
	float y_object_median;
	float z_object_median;

	std::vector<float> x_vals = std::vector<float>();
	std::vector<float> y_vals = std::vector<float>();
	std::vector<float> z_vals = std::vector<float>();
	x_vals.reserve(AABBs.size());
	y_vals.reserve(AABBs.size());
	z_vals.reserve(AABBs.size());

	for (int i = 0; i < AABBs.size(); i++) {
		x_vals[i] = AABBs[i].centre_point.x;
		y_vals[i] = AABBs[i].centre_point.y;
		z_vals[i] = AABBs[i].centre_point.z;
	}

	std::sort(x_vals.begin(), x_vals.end());
	std::sort(y_vals.begin(), y_vals.end());
	std::sort(z_vals.begin(), z_vals.end());


	x_object_median = x_vals[AABBs.size() / 2];
	y_object_median = y_vals[AABBs.size() / 2];
	z_object_median = z_vals[AABBs.size() / 2];

	//choose 9 equally spaced possibilities on each axis and check cost.
	float height = node->aabb.bounds[1].y - node->aabb.bounds[0].y;
	float width = node->aabb.bounds[1].x - node->aabb.bounds[0].x;
	float depth = node->aabb.bounds[1].z - node->aabb.bounds[0].z;

	float SA = 2 * (height * width) + 2 * (height * depth) + 2 * (width * depth);
	float parent_cost = SA * AABBs.size();
	split_cost = parent_cost;
	char split_axis = 'x';
	split_position = x_spatial_median;
	//x axis
	float min = fminf(x_object_median, x_spatial_median);
	float max = fmaxf(x_object_median, x_spatial_median);
	float step = (max - min) / 9.0f;
	for (int i = 0; i < 9; i++) {
		float split = (min + (step * i));
		float tempCost;

		float widthL = split - node->aabb.bounds[0].x;
		float widthR = node->aabb.bounds[1].x - split;

		float SALeft = 2 * (height * widthL) + 2 * (height * depth) + 2 * (widthL * depth);
		float SARight = 2 * (height * widthR) + 2 * (height * depth) + 2 * (widthR * depth);

		int n = 0;
		while (x_vals[n] < split) {
			n++;
		}

		tempCost = (SALeft * (float)n) + (SARight * ((float)AABBs.size() - (float)n));
		if (tempCost < split_cost) {
			split_cost = tempCost;
			split_axis = 'x';
			split_position = split;
		}
	}

	//y axis
	min = fminf(y_object_median, y_spatial_median);
	max = fmaxf(y_object_median, y_spatial_median);
	step = (max - min) / 9.0f;
	for (int i = 0; i < 9; i++) {
		float split = (min + (step * i));
		float tempCost;

		float heightL = split - node->aabb.bounds[0].y;
		float heightR = node->aabb.bounds[1].y - split;

		float SALeft = 2 * (heightL * width) + 2 * (heightL * depth) + 2 * (width * depth);
		float SARight = 2 * (heightR * width) + 2 * (heightR * depth) + 2 * (width * depth);

		int n = 0;
		while (y_vals[n] < split) {
			n++;
		}

		tempCost = (SALeft * (float)n) + (SARight * ((float)AABBs.size() - (float)n));
		if (tempCost < split_cost) {
			split_cost = tempCost;
			split_axis = 'y';
			split_position = split;
		}
	}

	//z axis
	min = fminf(z_object_median, z_spatial_median);
	max = fmaxf(z_object_median, z_spatial_median);
	step = (max - min) / 9.0f;
	for (int i = 0; i < 9; i++) {
		float split = (min + (step * i));
		float tempCost;

		float depthL = split - node->aabb.bounds[0].z;
		float depthR = node->aabb.bounds[1].z - split;

		float SALeft = 2 * (height * width) + 2 * (height * depthL) + 2 * (width * depthL);
		float SARight = 2 * (height * width) + 2 * (height * depthR) + 2 * (width * depthR);

		int n = 0;
		while (z_vals[n] < split) {
			n++;
		}

		tempCost = (SALeft * (float)n) + (SARight * ((float)AABBs.size() - (float)n));
		if (tempCost < split_cost) {
			split_cost = tempCost;
			split_axis = 'z';
			split_position = split;
		}
	}

	//if best cost is the same as the parent cost then don't split and just have a leaf node
	if (split_cost >= parent_cost) {
		return -1;
	}

	//return the first index for the right child as the point of partition
	split_position = split_position;

	switch (split_axis) {
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

void split_reference(AABB& left, AABB& right, AABB& current, int axis, float split_position, std::vector<vec3> t_vertices, std::vector<int> t_indices) {
	vec3 vertices[3] = { t_vertices[t_indices[current.index * 3]], t_vertices[t_indices[current.index * 3 + 1]], t_vertices[t_indices[current.index * 3 + 2]] };

	vec3 v1 = vertices[2];
	for (int i = 0; i < 3; i++) {
		vec3 v0 = v1;
		v1 = vertices[i];
		float v0_position = v0[axis];
		float v1_position = v1[axis];

		if (v0_position <= split_position) {
			left.bounds[0] = glm::min(left.bounds[0], v0);
			left.bounds[1] = glm::max(left.bounds[1], v0);
		}
		if (v1_position >= split_position) {
			right.bounds[0] = glm::min(right.bounds[0], v0);
			right.bounds[1] = glm::max(right.bounds[1], v0);
		}

		//if the edge instersects the plane then it is in both boxes so grow both
		if ((v0_position < split_position && v1_position > split_position) || (v0_position > split_position && v1_position < split_position))
		{
			float t = (split_position - v0_position) / (v1_position - split_position);
			vec3 intersection = v0 * (1.f - t) + v1 * t;

			left.bounds[0] = glm::min(left.bounds[0], intersection);
			left.bounds[1] = glm::max(left.bounds[1], intersection);

			right.bounds[0] = glm::min(right.bounds[0], intersection);
			right.bounds[1] = glm::max(right.bounds[1], intersection);
		}
	}

	left.bounds[1][axis] = split_position;
	left.bounds[0] = max(left.bounds[0], current.bounds[0]);
	left.bounds[1] = min(left.bounds[1], current.bounds[1]);

	right.bounds[0][axis] = split_position;
	right.bounds[0] = max(right.bounds[0], right.bounds[0]);
	right.bounds[1] = min(right.bounds[1], right.bounds[1]);
}

int BVH::SpatialPartition(BVHNode* node, std::vector<AABB> AABBs, float& split_position, float& split_cost, std::vector<vec3>& t_vertices, std::vector<int>& t_indices) {
	
	float height = node->aabb.bounds[1].y - node->aabb.bounds[0].y;
	float width = node->aabb.bounds[1].x - node->aabb.bounds[0].x;
	float depth = node->aabb.bounds[1].z - node->aabb.bounds[0].z;
	float SA = 2 * (height * width) + 2 * (height * depth) + 2 * (width * depth);
	float parent_cost = SA * AABBs.size();
	split_cost = parent_cost;

	//setup bins
	vec3 bin_sizes = (node->aabb.bounds[1] - node->aabb.bounds[0]) / 9.0f;
	int bin_index[3][9] = { {0,1,2,3,4,5,6,7,8}, {9,10,11,12,13,14,15,16,17}, {18,19,20,21,22,23,24,25,26} };
	AABB bin_bounds[27];
	float bin_enter[27];
	float bin_exit[27];

	float bin_num = 9.0f;

	for (int i = 0; i < 27; i++) {
		bin_bounds[i] = AABB();
		bin_enter[i] = 0;
		bin_exit[i] = 0;
	}

	//chopped binning
	for (int i = 0; i < AABBs.size(); i++) {
		AABB current_reference = AABBs[i];

		vec3 first_bin = (current_reference.bounds[0] - node->aabb.bounds[0]) / bin_sizes;
		first_bin = clamp(first_bin, vec3(0, 0, 0), vec3(bin_num - 1, bin_num - 1, bin_num - 1));

		vec3 last_bin = (current_reference.bounds[1] - node->aabb.bounds[0]) / bin_sizes;
		last_bin = clamp(last_bin, first_bin, vec3(bin_num - 1, bin_num - 1, bin_num - 1));

		//in case all AABBs are aligned on an axis
		for (int j = 0; j < 3; j++) {
			if (bin_sizes[j] == 0) {
				first_bin[j] = 0.0f;
				last_bin[j] = 0.0f;
			}
		}

		int first_bin_index[3] = { first_bin.x, first_bin.y, first_bin.z };
		int last_bin_index[3] = { last_bin.x, last_bin.y, last_bin.z };

		for (int j = 0; j < 3; j++) {
			for (int k = first_bin_index[j]; k < last_bin_index[j]; k++) {
				AABB left = AABB();
				AABB right = AABB();

				split_reference(left, right, current_reference, j, node->aabb.bounds[0][j] + bin_sizes[j] * (i + 1.0f), t_vertices, t_indices);
			}

			bin_bounds[bin_index[j][last_bin_index[j]]].bounds[0] = vec3(0);

			bin_bounds[bin_index[j][last_bin_index[j]]].bounds[0] = glm::min(bin_bounds[bin_index[j][last_bin_index[j]]].bounds[0], current_reference.bounds[0]);
			bin_bounds[bin_index[j][last_bin_index[j]]].bounds[1] = glm::max(bin_bounds[bin_index[j][last_bin_index[j]]].bounds[1], current_reference.bounds[0]);
			
			bin_bounds[bin_index[j][last_bin_index[j]]].bounds[0] = glm::min(bin_bounds[bin_index[j][last_bin_index[j]]].bounds[0], current_reference.bounds[1]);
			bin_bounds[bin_index[j][last_bin_index[j]]].bounds[1] = glm::max(bin_bounds[bin_index[j][last_bin_index[j]]].bounds[1], current_reference.bounds[1]);
			
			bin_enter[bin_index[j][first_bin_index[j]]]++;
			bin_exit[bin_index[j][last_bin_index[j]]]++;
		}
	}

	AABB right_bins[8];
	AABB left_bins[9];
	int split_axis = 0;
	//find best split plane
	for (int i = 0; i < 3; i++) {

		//prepopulate all the right bins
		AABB right_temp;
		for (int k = 8; k > 0; k--) {
			right_temp.bounds[0] = glm::min(right_temp.bounds[0], bin_bounds[bin_index[i][k]].bounds[0]);
			right_temp.bounds[1] = glm::max(right_temp.bounds[1], bin_bounds[bin_index[i][k]].bounds[0]);
				 								  
			right_temp.bounds[0] = glm::min(right_temp.bounds[0], bin_bounds[bin_index[i][k]].bounds[1]);
			right_temp.bounds[1] = glm::max(right_temp.bounds[1], bin_bounds[bin_index[i][k]].bounds[1]);

			right_bins[k - 1] = right_temp;
		}

		
		AABB left_temp;
		int right_count = AABBs.size();
		int left_count = 0;
		for (int k = 1; k < 9; k++) {
			left_temp.bounds[0] = glm::min(left_temp.bounds[0], bin_bounds[bin_index[i][k]].bounds[0]);
			left_temp.bounds[1] = glm::max(left_temp.bounds[1], bin_bounds[bin_index[i][k]].bounds[0]);
										   
			left_temp.bounds[0] = glm::min(left_temp.bounds[0], bin_bounds[bin_index[i][k]].bounds[1]);
			left_temp.bounds[1] = glm::max(left_temp.bounds[1], bin_bounds[bin_index[i][k]].bounds[1]);

			left_count += bin_enter[bin_index[i][k]];
			right_count -= bin_exit[bin_index[i][k]];

			float heightL = left_temp.bounds[1][1] - left_temp.bounds[0][1];
			float widthL = left_temp.bounds[1][0] - left_temp.bounds[0][0];
			float depthL = left_temp.bounds[1][2] - left_temp.bounds[0][2];

			float heightR= right_bins[k - 1].bounds[1][1] - right_bins[k - 1].bounds[0][1];
			float widthR = right_bins[k - 1].bounds[1][0] - right_bins[k - 1].bounds[0][0];
			float depthR = right_bins[k - 1].bounds[1][2] - right_bins[k - 1].bounds[0][2];

			float left_SA = (2 * heightL * widthL) + (2 * heightL * depthL) + (2 * widthL * depthL);
			float right_SA = (2 * heightR * widthR) + (2 * heightR * depthR) + (2 * widthR * depthR);

			float SAH = left_SA * left_count + right_SA * right_count;

			if (SAH < split_cost) {
				split_cost = SAH;
				split_axis = i;
				split_position = node->aabb.bounds[0][i] + (bin_sizes[i] * i);
			}
		}
	}

	return split_axis;
}

void BVH::MakeObjectSplit(std::vector<AABB>& left_AABBs, std::vector<AABB>& right_AABBs, int split_axis_object, float split_position, BVHNode* node, std::vector<AABB>& AABBs) {

	for (int i = 0; i < AABBs.size(); i++) {
		float centre_check;
		if (split_axis_object == 0) {
			centre_check = AABBs[i].centre_point.x;
		}
		else if (split_axis_object == 1) {
			centre_check = AABBs[i].centre_point.y;
		}
		else {
			centre_check = AABBs[i].centre_point.z;
		}

		if (centre_check <= split_position) {
			left_AABBs.push_back(AABBs[i]);
		}
		else {
			right_AABBs.push_back(AABBs[i]);
		}
	}


	node->left->CalculateBounds(left_AABBs);
	node->right->CalculateBounds(right_AABBs);
}

void BVH::MakeSpatialSplit(std::vector<AABB>& left_AABBs, std::vector<AABB>& right_AABBs, int split_axis_spatial, float split_position, BVHNode* node, std::vector<AABB>& AABBs){
	for (int i = 0; i < AABBs.size(); i++) {
		if (AABBs[i].bounds[1][split_axis_spatial] < split_position) {
			//printf("in left\n");
			left_AABBs.push_back(AABBs[i]);
		}else if(AABBs[i].bounds[0][split_axis_spatial] >= split_position) {
			//printf("in right\n");
			right_AABBs.push_back(AABBs[i]);
		}
		else {
			left_AABBs.push_back(AABBs[i]);
			right_AABBs.push_back(AABBs[i]);
		}
	}

	//printf("%d %d %d \n", left_AABBs.size(), right_AABBs.size(), AABBs.size());
	node->left->CalculateBounds(left_AABBs);
	node->right->CalculateBounds(right_AABBs);
	node->left->aabb.bounds[1][split_axis_spatial] = split_position;
	node->right->aabb.bounds[0][split_axis_spatial] = split_position;
}

void printBTr(const std::string& prefix, const MBVHNode* node, bool isLeft)
{
	if (node != nullptr)
	{
		std::cout << prefix;

		std::cout << (isLeft ? "---" : "---");

		// print the value of the node
		if (!node->is_leaf) {
			std::cout << "I" << std::endl;
		}
		else {
			std::cout << node->leaf_triangles.size() << std::endl;
		}
		// enter the next tree level - left and right branch
		for (int i = 0; i < node->children.size(); i++) {
			printBTr(prefix + (isLeft ? "   " : "    "), node->children[i], true);
		}
	}
}

void printBT(const MBVHNode* node)
{
	printBTr("", node, false);
}

// pass the root node of your binary tree

MBVHNode* mergeNodes(BVHNode* node) {
	MBVHNode* m_node = new MBVHNode();
	if (!node->is_leaf) {
		if (!node->left->is_leaf) {
			m_node->children.push_back(mergeNodes(node->left->left));
			m_node->children.push_back(mergeNodes(node->left->right));
		}
		else {
			m_node->children.push_back(mergeNodes(node->left));
		}

		if (!node->right->is_leaf) {
			m_node->children.push_back(mergeNodes(node->right->left));
			m_node->children.push_back(mergeNodes(node->right->right));
		}
		else {
			m_node->children.push_back(mergeNodes(node->right));
		}

		m_node->is_leaf = false;
		m_node->CalculateBounds();
	}
	else {
		m_node->is_leaf = true;
		m_node->leaf_triangles = node->leaf_triangles;
		m_node->aabb = node->aabb;
	}

	return m_node;
}

void BVH::Collapse() {
	MBVHNode* root = mergeNodes(this->root_node);
	this->root_node_mbvh = root;
}

/**CACHE FRIENDLY BVH FUNCTIONS**/
int BVH::count_nodes(BVHNode* root) {
	if (!root->is_leaf) {
		//printf("BVH Leaf: %f %f %f | %f %f %f \n", root->aabb.bounds[0][0], root->aabb.bounds[0][1], root->aabb.bounds[0][2], root->aabb.bounds[1][0], root->aabb.bounds[1][1], root->aabb.bounds[1][2]);
		return 1 + count_nodes(root->left) + count_nodes(root->right);
	}
	else {
		this->test_b++;
	}

	return 1;
}

int BVH::count_nodes(MBVHNode* root) {
	if (!root->is_leaf) {
		int count = 0;
		//printf("MBVH Leaf: %f %f %f | %f %f %f \n", root->aabb.bounds[0][0], root->aabb.bounds[0][1], root->aabb.bounds[0][2], root->aabb.bounds[1][0], root->aabb.bounds[1][1], root->aabb.bounds[1][2]);
		for (int i = 0; i < root->children.size(); i++) {
			count += count_nodes(root->children[i]);
		}
		return 1 + count;
	}
	else {
		this->test_m++;
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

void BVH::PopulateCFBVH(unsigned int& current_index, unsigned int& cumulative_index, unsigned int& triangle_index, MBVHNode* node) {
	this->cf_mbvh[current_index].bounds[0] = node->aabb.bounds[0];
	this->cf_mbvh[current_index].bounds[1] = node->aabb.bounds[1];

	if (!(node->is_leaf)) {
		//std::cout << "INNER " << cumulative_index << std::endl;
		std::vector<unsigned int> child_indices = std::vector<unsigned int>();
		cumulative_index += node->children.size();
		cf_mbvh[current_index].u.inner.child_index = cumulative_index + 1 - node->children.size();
		cf_mbvh[current_index].u.inner.child_count = node->children.size();

		for (int i = 0; i < node->children.size(); i++) {
			unsigned int u = cf_mbvh[current_index].u.inner.child_index + i;
			PopulateCFBVH(u, cumulative_index, triangle_index, node->children[i]);
		}

	}
	else {
		unsigned int count = (unsigned int)node->leaf_triangles.size();
		cf_mbvh[current_index].u.leaf.count = 0x80000000 | count; //basically sets the first bit of count a
		cf_mbvh[current_index].u.leaf.index_first_tri = triangle_index;
		//std::cout << "LEAF " << cumulative_index << std::endl;
		for (unsigned int i = 0; i < count; i++) {
			this->mbvh_triangle_indices[triangle_index] = node->leaf_triangles[i];
			triangle_index++;
		}
	}
}

void printBTmr(const std::string& prefix, int current, bool isLeft, MBVHNode_CacheFriendly* cf_mbvh)
{
	printf("Node Sizes MBVH = %d BHV = %d \n", sizeof(MBVHNode_CacheFriendly), sizeof(BVHNode_CacheFriendly));
}

void BVH::ConstructCacheFriendly(int tri_count) {
	int bvh_count = count_nodes(this->root_node);
	unsigned int cumulative_index = 0;
	unsigned int current_index = 0;
	unsigned int triangle_index = 0;
	this->cf_bvh = new BVHNode_CacheFriendly[bvh_count];
	this->triangle_indices = new int[tri_count];
	PopulateCFBVH(cumulative_index, triangle_index, this->root_node);

	int mbvh_count = count_nodes(this->root_node_mbvh);
	cumulative_index = 0;
	current_index = 0;
	triangle_index = 0;
	this->cf_mbvh = new MBVHNode_CacheFriendly[mbvh_count];
	this->mbvh_triangle_indices = new int[tri_count];
	PopulateCFBVH(current_index, cumulative_index, triangle_index, this->root_node_mbvh);
	printBTmr("", mbvh_count, false, this->cf_mbvh);
	//allocate mbvh in cuda memory
	cudaAssert(Malloc(&(cf_bvh_gpu), bvh_count * sizeof(BVHNode_CacheFriendly)));
	cudaAssert(Malloc(&(triangle_indices_gpu), tri_count * sizeof(int)));
	cudaAssert(Memcpy(cf_bvh_gpu, cf_bvh, bvh_count * sizeof(BVHNode_CacheFriendly), cudaMemcpyHostToDevice));
	cudaAssert(Memcpy(triangle_indices_gpu, triangle_indices, tri_count * sizeof(int), cudaMemcpyHostToDevice));

	cudaAssert(Malloc(&(cf_mbvh_gpu), mbvh_count * sizeof(MBVHNode_CacheFriendly)));
	cudaAssert(Malloc(&(mbvh_triangle_indices_gpu), tri_count * sizeof(int)));
	cudaAssert(Memcpy(cf_mbvh_gpu, cf_mbvh, mbvh_count * sizeof(MBVHNode_CacheFriendly), cudaMemcpyHostToDevice));
	cudaAssert(Memcpy(mbvh_triangle_indices_gpu, mbvh_triangle_indices, tri_count * sizeof(int), cudaMemcpyHostToDevice));
}