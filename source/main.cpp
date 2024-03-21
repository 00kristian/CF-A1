//Computational Fabrication Assignment #1
// By David Levin 2014
#include <iostream>
#include <vector>
#include "../include/CompFab.h"
#include "../include/Mesh.h"





//Ray-Triangle Intersection
//Returns 1 if triangle and ray intersect, 0 otherwise
int rayTriangleIntersection(CompFab::Ray &ray, CompFab::Triangle &triangle)
{
    /********* ASSIGNMENT *********/
    /* Ray-Triangle intersection test: Return 1 if ray intersects triangle, 
     * 0 otherwise */

     //Code inspired from https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm#:~:text=The%20M%C3%B6ller%E2%80%93Trumbore%20ray%2Dtriangle,the%20plane%20containing%20the%20triangle 

    //Find the vectors for the two edges sharing the vertex v1: V2V1 and V3V1 
    CompFab::Vec3 edge1 = (triangle.m_v2 - triangle.m_v1);
    CompFab::Vec3 edge2 = (triangle.m_v3 - triangle.m_v1);


    //find the vector normal vector of the plane containing the triangle
    CompFab::Vec3 n = ray.m_direction % edge2;

    //define the determinant 
    float det = edge1 * n;

    //Check if the ray is parallel to the triangle -> if the determinant is close to 0 / 0 then the ray is parallel to the triangle
    if (det > -EPSILON && det < EPSILON)
    {
        return 0;
    }

    //This is the inside-outside test

    //Find the inverse of the determinant
    float inv_det = 1.0 / det;

    //Find the vector from the origin of the ray to the v1 of the triangle
    CompFab::Vec3 s = ray.m_origin - triangle.m_v1;

    // Find the u parameter
    float u = inv_det * (s * n);

    //Check if the u parameter is outside the triangle
    if (u < 0.0 || u > 1.0)
    {
        return 0;
    }

    //find the vector from the origin of the ray to the v1 of the triangle
    CompFab::Vec3 q = s % edge1;

    //Find the v parameter
    float v = inv_det * (ray.m_direction * q);

    //Check if the v parameter is outside the triangle
    if (v < 0.0 || u + v > 1.0)
    {
        return 0;
    }

    //Find the t parameter which is the distance from the origin of the ray to the intersection point

    float t = inv_det * (edge2 * q);

    //Check if the t parameter is close to 0, if it is bigger than 0 then we have an intersection
    if (t > EPSILON)
    {
        return 1;
    }
    else
    {
        return 0;
    }

}

//Triangle list (global)
typedef std::vector<CompFab::Triangle> TriangleList;

TriangleList g_triangleList;
CompFab::VoxelGrid *g_voxelGrid;

//Number of intersections with surface made by a ray originating at voxel and cast in direction.
int numSurfaceIntersections(CompFab::Vec3 &voxelPos, CompFab::Vec3 &dir)
{ 
    //Initialize the number of hits to 0
    unsigned int numHits = 0;
    
    /********* ASSIGNMENT *********/
    /* Check and return the number of times a ray cast in direction dir, 
     * from voxel center voxelPos intersects the surface */

    //Create ray with voxelposistion and direction
    CompFab::Ray ray(voxelPos, dir);

    //Iterate over all the triangles in the triangle list and see if the ray intersects with them 

    for (int i = 0; i < g_triangleList.size(); i++)
    {
        //If the ray intersects with the triangle then we increment the number of hits
        if (rayTriangleIntersection(ray, g_triangleList[i]))
        {
			numHits+=1;
		}
	}

    return numHits;
}

bool loadMesh(char *filename, unsigned int dim)
{
    g_triangleList.clear();
    
    Mesh *tempMesh = new Mesh(filename, true);
    
    CompFab::Vec3 v1, v2, v3;

    //copy triangles to global list
    for (unsigned int tri = 0; tri < tempMesh->t.size(); ++tri)
    {
        v1 = tempMesh->v[tempMesh->t[tri][0]];
        v2 = tempMesh->v[tempMesh->t[tri][1]];
        v3 = tempMesh->v[tempMesh->t[tri][2]];
        g_triangleList.push_back(CompFab::Triangle(v1, v2, v3));
    }

    //Create Voxel Grid
    CompFab::Vec3 bbMax, bbMin;
    BBox(*tempMesh, bbMin, bbMax);

    //Build Voxel Grid
    double bbX = bbMax[0] - bbMin[0];
    double bbY = bbMax[1] - bbMin[1];
    double bbZ = bbMax[2] - bbMin[2];
    double spacing;

    if (bbX > bbY && bbX > bbZ)
    {
        spacing = bbX / (double)(dim - 2);
    }
    else if (bbY > bbX && bbY > bbZ) {
        spacing = bbY / (double)(dim - 2);
    }
    else {
        spacing = bbZ / (double)(dim - 2);
    }

    CompFab::Vec3 hspacing(0.5 * spacing, 0.5 * spacing, 0.5 * spacing);

    g_voxelGrid = new CompFab::VoxelGrid(bbMin - hspacing, dim, dim, dim, spacing);

    delete tempMesh;

    return true;

}

void saveVoxelsToObj(const char* outfile)
{

    Mesh box;
    Mesh mout;
    int nx = g_voxelGrid->m_dimX;
    int ny = g_voxelGrid->m_dimY;
    int nz = g_voxelGrid->m_dimZ;
    double spacing = g_voxelGrid->m_spacing;

    CompFab::Vec3 hspacing(0.5 * spacing, 0.5 * spacing, 0.5 * spacing);

    for (int ii = 0; ii < nx; ii++) {
        for (int jj = 0; jj < ny; jj++) {
            for (int kk = 0; kk < nz; kk++) {
                if (!g_voxelGrid->isInside(ii, jj, kk)) {
                    continue;
                }
                CompFab::Vec3 coord(0.5f + ((double)ii) * spacing, 0.5f + ((double)jj) * spacing, 0.5f + ((double)kk) * spacing);
                CompFab::Vec3 box0 = coord - hspacing;
                CompFab::Vec3 box1 = coord + hspacing;
                makeCube(box, box0, box1);
                mout.append(box);
            }
        }
    }

    mout.save_obj(outfile);
}


int main(int argc, char** argv)
{

    unsigned int dim = 32; //dimension of voxel grid (e.g. 32x32x32)

    //Load OBJ
    if (argc < 3)
    {
        std::cout << "Usage: Voxelizer InputMeshFilename OutputMeshFilename \n";
        exit(0);
    }
     
    std::cout << "Load Mesh : " << argv[1] << "\n";
    loadMesh(argv[1], dim);




    //Cast ray, check if voxel is inside or outside
    //even number of surface intersections = outside (OUT then IN then OUT)
    // odd number = inside (IN then OUT)
    CompFab::Vec3 voxelPos;
    CompFab::Vec3 direction(1.0, 0.0, 0.0);


    /********* ASSIGNMENT *********/
    /* Iterate over all voxels in g_voxelGrid and test whether they are inside our outside of the
     * surface defined by the triangles in g_triangleList */

    //iterate over the x dimension
    for (int i = 0; i < g_voxelGrid->m_dimX; i++)
    {
        //iterate over the y dimension
        for (int j = 0; j < g_voxelGrid->m_dimY; j++)
        {
            //iterate over the z dimension
            for (int k = 0; k < g_voxelGrid->m_dimZ; k++)
            {
                //Create a voxel position with the current i,j,k and the spacing of the voxel grid
				voxelPos = g_voxelGrid->m_lowerLeft + CompFab::Vec3(i * g_voxelGrid->m_spacing, j * g_voxelGrid->m_spacing, k * g_voxelGrid->m_spacing);

                //Check if the number of surface intersections is even or odd
                if (numSurfaceIntersections(voxelPos, direction) % 2 == 0)
                {
                    //even number = outside
					g_voxelGrid->isInside(i,j,k) = false;
				}
                else
                {
                    //odd number = inside
					g_voxelGrid->isInside(i,j,k) = true;
				}             
			}
		}
        
    }
    
    //Write out voxel data as obj
    saveVoxelsToObj(argv[2]);
    
    delete g_voxelGrid;
}