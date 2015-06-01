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

#ifndef __SPAAM_SVD__
#define __SPAAM_SVD__

#include <vector>

////////BOOST and LAPACK includes////////
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/lapack/gesvd.hpp>

////////use the BOOST NUMERIC namespace///////
using namespace boost::numeric;

/********************************************************************************************
Master Class which stores Vector and Matrix classes, as well as SPAAM methods to simplify the
use of the SVD functions
*********************************************************************************************/
template< typename W, typename S >
class Correspondence_Pair
{
	public:
		//Default Constructor//
		Correspondence_Pair( )
		{ worldPoint.clear(); screenPoint.clear(); }

		//Parameter Constructor//
		Correspondence_Pair( W x1, W y1, W z1, S x2, S y2 )
		{ worldPoint(0) = x1; worldPoint(1) = y1; worldPoint(2) = z1;
			screenPoint(0) = x2; screenPoint(1) = y2; }

		//Default Destructor//
		~Correspondence_Pair()
		{		}

		//Correspondence Points//
		ublas::c_vector< W, 3 > worldPoint;
		ublas::c_vector< S, 2 > screenPoint;
};

/******************************************************************************************
SPAAM_SVD

This class provides the interface for storing a list of Correspondence_Pair objects
and performing a Singular Value Decomposition (SVD) on them. The result of the SVD is
a 3x4 projection matrix that describes the homography (mapping) of 3D world points to
2D screen points (or vice versa).
******************************************************************************************/
template< typename W, typename S >
class SPAAM_SVD
{
	public:

		//Default Constructor//
		SPAAM_SVD( )
		{ modMatrixScreen.resize(3, 3); modMatrixWorld.resize(4, 4); Proj3x4.resize(3, 4); }

		//Default Destructor//
		~SPAAM_SVD()
		{		}
			
		////////Vector containing all of the correspondence pair point values////////
		std::vector< Correspondence_Pair<W, S> > corr_points;

	private:
		////////Data Members for storing intermediate values required for performing the SVD operation////////
		
		////Normalization Components for World Points////
		ublas::c_vector< W, 3 > fromShift;
		ublas::c_vector< W, 3 > fromScale;
		////Normalization Components for Screen Points////
		ublas::c_vector< S, 2 > toShift;
		ublas::c_vector< S, 2 > toScale;
		////Normalization Matrix for World Points////
		ublas::matrix< W, boost::numeric::ublas::column_major > modMatrixWorld;
		////Normalization Matrix for Screen Points////
		ublas::matrix< S, boost::numeric::ublas::column_major > modMatrixScreen;

		////Final 3 x 4 Projection Matrox////
		ublas::matrix< double, boost::numeric::ublas::column_major> Proj3x4;

		////////methods for calculating the SVD decomposition of the SPAAM linear equation matrices///////

		/*****************************************************************************
		This function is used to normalize the screen and world points (the correspondence
		point pairs) so that their values span the same range. This is required in order
		to use values from 2 different spaces (world and screen space) in the same
		calculations. The function should be called before the SVD operation begins.
		The actual values of the 2D and 3D points are not changed, but the normalization value
		for each is calculated for use later in the SVD operation.

		There is no return value for this function
		*****************************************************************************/
		void estimateNormalizationParameters( )
		{
			////determine the number of points to be normalized////
			const std::size_t n_pts = std::distance( corr_points.begin(), corr_points.end() );

			////start with all 0's////
			fromShift.clear(); 
			fromScale.clear(); 
			toShift.clear();
			toScale.clear();
			
			////compute mean and mean of square////
			for ( vector< Correspondence_Pair<W, S> >::const_iterator it( corr_points.begin() ); it < corr_points.end(); ++it )
			{	
				fromShift = fromShift + (*it).worldPoint;
				fromScale = fromScale + boost::numeric::ublas::element_prod( (*it).worldPoint, (*it).worldPoint );

				toShift = toShift + (*it).screenPoint;
				toScale = toScale + boost::numeric::ublas::element_prod( (*it).screenPoint, (*it).screenPoint );
			}
			fromShift *= static_cast< W >( 1 ) / n_pts;
			fromScale *= static_cast< W >( 1 ) / n_pts;
			toShift *= static_cast< S >( 1 ) / n_pts;
			toScale *= static_cast< S >( 1 ) / n_pts;
	
			////compute standard deviation////
			for ( std::size_t i = 0; i < 3; i++ )
				fromScale( i ) = std::sqrt( fromScale( i ) - ( fromShift( i ) * fromShift( i ) ) );
			for ( std::size_t i = 0; i < 2; i++ )
				toScale( i ) = std::sqrt( static_cast<double>(toScale( i ) - ( toShift( i ) * toShift( i )) ) );

			////end of function////
		}

		/*****************************************************************************
		This function is used to create the normalization matrix based upon the normalized
		values from the correspondence point pairs. This is required in order to transform
		the resulting matrix back into real space. The SVD operation is performed on the
		normalized values of the 2D and 3D points. Therefore, the result is also in the
		normalized space. This operation produces the reverse value so that the SVD
		result can be returned back into the proper scale range.

		There is no return value for this function
		*****************************************************************************/
		void generateNormalizationMatrix( )
		{
			////compute correction matrix////
			////Start with All 0's////
			modMatrixWorld.clear();
			modMatrixScreen.clear();

			////create homogeneous matrix////
			modMatrixWorld( 3, 3 ) = static_cast< W >( 1 );
			modMatrixScreen( 2, 2 ) = static_cast< S >( 1 );
	
			////honestly I'm not sure what this is for////
			//if ( true )
			{
				for ( std::size_t i = 0; i < 2; i++ )
				{
					modMatrixScreen( i, i ) = toScale( i );
					modMatrixScreen( i, 2 ) = toShift( i );
				}
			}//else
			{
				for ( std::size_t i = 0; i < 3; i++ )
				{
					modMatrixWorld( i, i ) = static_cast< W >( 1 ) / fromScale( i );
					modMatrixWorld( i, 3 ) = -modMatrixWorld( i, i ) * fromShift( i );
			
				}
			}
		}

		/*****************************************************************************
		This function is used to perform the Singular Value Decomposition on the list
		of Correspondence_Pair objects. The normalization parameters for each set of
		values is determined using the 'estimateNormalizationParameters()' function.

		A system of linear equations (constructed to solve for the unknown parameters
		of the 3x4 projection matrix for the SPAAM calibration) is also constructed and
		solved using the SVD calculation.

		The resulting 3x4 matrix is then converted from the normalized range back into
		the proper scale using the values obtained from the 'generateNormalizationMatrix( )'
		function call. The final 3x4 projection matrix is stored in the class member Proj3x4
		and is used for the return value.

		The return value for this function is a 3x4 Matrix
		*****************************************************************************/
		public:
		ublas::matrix< double, boost::numeric::ublas::column_major> projectionDLTImpl( )
		{
			////minimum of 6 correspondence points required to solve////
			assert( corr_points.size() >= 6 );

			// normalize input points
			estimateNormalizationParameters( );//fromPoints.begin(), fromPoints.end(), fromShift, fromScale );

			// construct equation system
			ublas::matrix< double, boost::numeric::ublas::column_major > A(2*corr_points.size(), 12);
			for ( unsigned i = 0; i < corr_points.size(); i++ )
			{
				ublas::c_vector< S, 2 > to = ublas::element_div( corr_points[i].screenPoint - toShift, toScale );
				ublas::c_vector< W, 3 > from = ublas::element_div( corr_points[i].worldPoint - fromShift, fromScale );

				A( i * 2,  0 ) = A( i * 2, 1 ) = A( i * 2, 2 ) = A( i * 2, 3 ) = 0;
				A( i * 2,  4 ) = -from( 0 );
				A( i * 2,  5 ) = -from( 1 );
				A( i * 2,  6 ) = -from( 2 );
				A( i * 2,  7 ) = -1;
				A( i * 2,  8 ) = to( 1 ) * from( 0 );
				A( i * 2,  9 ) = to( 1 ) * from( 1 );
				A( i * 2, 10 ) = to( 1 ) * from( 2 );
				A( i * 2, 11 ) = to( 1 );
				A( i * 2 + 1,  0 ) = from( 0 );
				A( i * 2 + 1,  1 ) = from( 1 );
				A( i * 2 + 1,  2 ) = from( 2 );
				A( i * 2 + 1,  3 ) = 1;
				A( i * 2 + 1,  4 ) = A( i * 2 + 1, 5 ) = A( i * 2 + 1, 6 ) = A( i * 2 + 1, 7 ) = 0;
				A( i * 2 + 1,  8 ) = -to( 0 ) * from( 0 );
				A( i * 2 + 1,  9 ) = -to( 0 ) * from( 1 );
				A( i * 2 + 1, 10 ) = -to( 0 ) * from( 2 );
				A( i * 2 + 1, 11 ) = -to( 0 );
			}

			// solve using SVD
			ublas::c_vector< double, 12 > s;
			ublas::matrix< double, boost::numeric::ublas::column_major > Vt(12, 12);
			ublas::matrix< double, boost::numeric::ublas::column_major > U( 2 * corr_points.size(), 2 * corr_points.size() );
			boost::numeric::bindings::lapack::gesvd( 'N', 'A', A, s, U, Vt );
			
			// copy result to 3x4 matrix
			Proj3x4( 0, 0 ) = Vt( 11, 0 ); Proj3x4( 0, 1 ) = Vt( 11, 1 ); Proj3x4( 0, 2 ) = Vt( 11,  2 ); Proj3x4( 0, 3 ) = Vt( 11,  3 );
			Proj3x4( 1, 0 ) = Vt( 11, 4 ); Proj3x4( 1, 1 ) = Vt( 11, 5 ); Proj3x4( 1, 2 ) = Vt( 11,  6 ); Proj3x4( 1, 3 ) = Vt( 11,  7 );
			Proj3x4( 2, 0 ) = Vt( 11, 8 ); Proj3x4( 2, 1 ) = Vt( 11, 9 ); Proj3x4( 2, 2 ) = Vt( 11, 10 ); Proj3x4( 2, 3 ) = Vt( 11, 11 );

			// reverse normalization
			generateNormalizationMatrix( );
			const ublas::matrix< double, boost::numeric::ublas::column_major > toCorrect(( modMatrixScreen ));
			ublas::matrix< W, boost::numeric::ublas::column_major > Ptemp(3, 4); Ptemp = ( ublas::prod( toCorrect, Proj3x4 ) );
			const ublas::matrix< double, boost::numeric::ublas::column_major > fromCorrect(( modMatrixWorld ));
			ublas::noalias( Proj3x4 ) = ublas::prod( Ptemp, fromCorrect );

			// normalize result to have a viewing direction of length 1 (optional)
			double fViewDirLen = sqrt( Proj3x4( 2, 0 ) * Proj3x4( 2, 0 ) + Proj3x4( 2, 1 ) * Proj3x4( 2, 1 ) + Proj3x4( 2, 2 ) * Proj3x4( 2, 2 ) );
	
			// if first point is projected onto a negative z value, negate matrix
			const ublas::vector< double > p1st(corr_points[0].worldPoint);
			if ( Proj3x4( 2, 0 ) * p1st( 0 ) + Proj3x4( 2, 1 ) * p1st( 1 ) + Proj3x4( 2, 2 ) * p1st( 2 ) + Proj3x4( 2, 3 ) < 0 )
				fViewDirLen = -fViewDirLen;
	
			Proj3x4 *= double( 1 ) / fViewDirLen;

			return Proj3x4;
		}
};

#endif //__SPAAM_SVD__