# SPAAM_BaseCode
Base code for performing the basic mathematical operations of a Single Point Active Alignment Method (SPAAM) calibration.

/*************************************************************************
This code is taken almost directly from the Ubitrack library:

Ubitrack - Library for Ubiquitous Tracking
 * Copyright 2006, Technische Universitaet Muenchen,
 *
 * This is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2.1 of
 * the License, or (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this software; if not, write to the Free
 * Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
 * 02110-1301 USA, or see the FSF site: http://www.fsf.org.

The source code is open source under the GLPL license. Please feel free
to use and modfiy the source code with proper citation where applicable
*************************************************************************/

# Project Setup
This project relies on the Boost extension library for C++ and the LAPACK (Linear Algebra Package) add-on library for Boost. The LAPACK static libraries (.lib) and dynamic link libraries (.dll) are included with this project. The Boost header files, however, are not included and must be downloaded seperately. Please visit the Boost official website to download the Boost library and header files ([Boost Library](http://www.boost.org)).

Once the Boost library is downloaded, extract the Boost header files either into the 'include' folder for this project or into another directory of your choosing (perhaps C:/Boost).

You must then set the include directory for the Boost headers within the Visual Studio 2010 project file.

![Include Directory](images/project_properties.PNG?raw=true)

This should be the only action required to setup the Visual Studio project (for Release build). The LAPACK libraries should already be present in the lib folder within the project as well as the dll's.

# Usage
This project is intended to provide a simple example of how the SPAAM_SVD and Correspondence_Pair objects defined in the SPAAM_SVD.h header file can be used to perform the calculations of a SPAAM calibration. Once the project is built, and the executable created, running the program should result in a 3x4 projection matrix being printed to the console window.

![Console Window](images/console_output.PNG?raw=true)

The exact numbers shown in the above image should always appear in the order shown.

## Correspondence_Pair Object
The Correspondence_Pair object is a simple contained for a 2D/3D alignment set. Each instance of Correspondence_Pair contains 2 vectors: worldPoint & screenPoint. The worldPoint stores the x, y, and z position data of a real world point and the screenPoint member stoes the x and y screen coordinate locations of the pixel aligned with the worldPoint during the calibration. The premise behind this object is to provide a stragihtforward container method for saving alignment data.

## SPAAM_SVD Object
The SPAAM_SVD object provides the interface by which a list of Correspondence_Pair objects can be stored and the Singular Value Decomposition (SVD) operation performed on them.

The list (really an std::vector of Correspondence_Pair objects) which stores the alignment pair data is the member 'corr_points'. By using the push_back method of this vector (i.e. corr_points.push_back...) a new alignment pair can be added to the list. Other operations such as reordering, deleting, inserting, etc. can also be performed on this list.

Once the list of Correspondnece_Pairs contains at least 6 pairs, the SVD operation can be performed to produce the homography, or mapping, from the 3D world points to the 2D screen points. The function 'projectionDLTImpl()' will perform the SVD operation on the 'corr_points' values. The resulting 3x4 matrix will be returned by the function but will also be stored in the SPAAM_SVD data member 'Proj3x4'.

The 'BuildGLMatrix3x4' function can be called was the SVD operation produces a result. Even though the name of the functions contains 3x4, the result of the operation is actually a 4x4 projection matrix usable by OpenGL. The data member 'projMat3x4', which stores the result, is a 16 element array of doubles which can be dierctly utilized by OpenGL for rendering. The 6 parameters that must be passed to the 'BuildGLMatrix3x4' function describe an orthographic viewing projection which is combined with the 3x4 matrix result of the SVD operation. The parameters to be passed in (in order): near clipping plane (float); far clipping plane (float); right most pixel (int); left most pixel (int); top most pixel (int); bottom most pixel (int).
