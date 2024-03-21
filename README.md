# CF-A1 by Kristian Berlin Jensen

I had some issues with the CMake command and couldn't create the voxeilizer.exe file. I tried to fix it but ended up doing this command when I build the project:
cmake .. -G "WinGW Makefiles". Github also didn't allow me to push the .obj files to the cloud. 


This allowed me to create a Makefile such that I could use make all. This created the voxeilizer.exe file which could be run with the additional OBJ files. 

To run the program with the sphere.obj file
1. Open the terminal in the folder
2. Go to the build folder (where the a1.exe is located)
3. Run the command: make all
4. Copy the sphere.obj file to the data/sphere/ folder and create a new file called: SphereVoxels32.obj  
5. Run the command: .\a1.exe ..\data\sphere\sphere.obj ..\data\sphere\SphereVoxels32.obj 


To run everything from the start with the spjere object:
1. Open the terminal in the folder
2. If the build folder is empty use the command: cmake .. -G "WinGW Makefiles"
3. Go to the build folder (where the voxeilizer.exe is located)
4. Run the command: make all
5. Copy the sphere.obj file to the data/sphere/ folder and create a new file called: SphereVoxels32.obj  
6. Run the command: .\a1.exe ..\data\sphere\sphere.obj ..\data\sphere\SphereVoxels32.obj

