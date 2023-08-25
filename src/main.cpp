
#include "hungarian_myself.h"

int main(void)
{
    // please use "-std=c++11" for this initialization of vector.
	vector<vector<int>> costMatrix = { { 10, 19, 8, 15, 0 }, 
									   { 10, 18, 7, 17, 0 }, 
									   { 13, 16, 9, 14, 0 }, 
									   { 12, 19, 8, 18, 0 },
									   {14, 32, 41, 34, 0} };

	Hungarian HungAlgo;
	vector<int> assignment;
	auto start_time = std::chrono::high_resolution_clock::now();
	HungAlgo.solve(costMatrix, assignment);

	auto end_time = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = end_time - start_time;

	std::cout << "Consume Time: " << elapsed.count() << " seconds" << std::endl;

	for (unsigned int x = 0; x < costMatrix.size(); x++)
		std::cout << x << "," << assignment[x] << "\t";
    cout<<endl;
	//std::cout << "\ncost: " << cost << std::endl;

	return 0;
}