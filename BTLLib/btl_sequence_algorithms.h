//
// btl_sequence_algorithms.h
//
// This file contains the header file Sequence Analysis algorithms 
//  
//
// Copyright (C) 1997, 1998 Birkbeck College, Malet Street, London, U.K.
// Copyright (C) 2004, 2005 University College, Gower Street, London, U.K.
//
// This library is free software; you can redistribute it and/or modify it 
// under the terms of the GNU Library General Public License as published by 
// the Free Software Foundation; either version 2 of the License, or (at your
// option) any later version.  This library is distributed in the hope
// that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
// PURPOSE.  See the GNU Library General Public License for more details.
// You should have received a copy of the GNU Library General Public
// License along with this library; if not, write to the Free Software
// Foundation, 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
///////////////////////////////////////////////////////////////////////////
//
// Authors: B.Sweeney, W.Pitt, M.A.Williams                
//         
///////////////////////////////////////////////////////////////////////////
// Modification History:
// Original version: B.S. (assistance from W.P)             Date: 30/08/98
// Last modified by: M.A.W                                  Date: 16/06/05
// Last modified by: R.Alcraft								Date: 24/02/20
// Comments:   
///////////////////////////////////////////////////////////////////////////

#if !defined (BTL_SEQUENCEALGORITHMS_H)
#define BTL_SEQUENCEALGORITHMS_H 1

#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <iterator>

#include "BTL.h"
#include "btl_numeric_vector.h"
#include "btl_matrix.h"


_BTL_BEGIN_NAMESPACE

using namespace std;

/**#: [Description ="
	 This template class contains generic algorithms for manipulating sequences
	 of data elements.<P>"]

	[Summary = "Sequence alignment and related algorithms"]
	[Authors = "B.Sweeney, W.R.Pitt, M.A.Williams."]
	[Files = "<A HREF=./btl_sequence_algorithms.h>btl_sequence_algorithms.h</A>"]
	[Dependencies="<A HREF=#numeric_vector>btl_numeric_vector.h</A>,
				   <A HREF=#matrix>btl_matrix.h</A>"]
*/

//.....................................................................................
// NEEDLEMAN & WUNSCH SIMILARITY SCORE
//.....................................................................................
/**#: [Description="Pairwise similarity score for two sequences."]
	  [Restrictions="The containers holding the sequences must have at least
					 forward iterators and an equality operator (==) should be defined
					 for pairs of elements. The score for any pair of elements is
			match_score (or mismatch_score) - total gap penalty.
			The total gap penalty = gap_penalty + (length_of_gap * gap_length_multiplier)"] */

template <class InputIterator, class T>
T needleman_wunsch_similarity(InputIterator first1, InputIterator last1,
	InputIterator first2, InputIterator last2,
	float match_score, float mismatch_score,
	float gap_penalty, float gap_length_multiplier,
	T init)
{

	// get the lengths of the sequences + 1
	int length_s = (last1 - first1) + 1;
	int length_t = (last2 - first2) + 1;

	// Create a score matrix
	matrix<float> matrix_s(length_s, length_t);

	// Calculate matrix values for the last row and column.

	// The scoring matrix is padded at the end by a blank

	matrix_s[length_s - 1][length_t - 1] = 0.0;

	int i, j;
	for (i = (length_s - 2); i > -1; i--)
	{
		matrix_s[i][length_t - 1] = (i - length_s + 1) * gap_length_multiplier - gap_penalty;
	}
	for (j = (length_t - 2); j > -1; j--)
	{
		matrix_s[length_s - 1][j] = (j - length_t + 1) * gap_length_multiplier - gap_penalty;
	}


	// Calculate the other values in the matrix.

	float cval, rval, dval, tmp_max;
	int ii, jj;
	for (i = (length_s - 2); i > -1; i--)
	{
		for (j = (length_t - 2); j > -1; j--)
		{
			// We now calculate the maximum of the row cross-ways and to the right, 
			// and the column cross-ways and down.

			// Get the cross-value.
			cval = matrix_s[i + 1][j + 1];

			// Get the down-value.
			dval = matrix_s[i + 1][j + 1] - gap_length_multiplier;
			for (ii = i + 2; ii < length_s; ii++)
			{
				if ((matrix_s[ii][j + 1] - (ii - i - 1) * gap_length_multiplier) > dval)
					dval = matrix_s[ii][j + 1] - (ii - i - 1) * gap_length_multiplier;
			}


			// Get the right-value       
			rval = matrix_s[i + 1][j + 1] - gap_length_multiplier;
			for (jj = j + 2; jj < length_t; jj++)
			{
				if ((matrix_s[i + 1][jj] - (jj - j - 1) * gap_length_multiplier) > rval)
					rval = matrix_s[i + 1][jj] - (jj - j - 1) * gap_length_multiplier;
			}

			tmp_max = max(dval, rval);
			tmp_max = max((tmp_max - gap_penalty), cval);
			if (*(first1 + i) == *(first2 + j))
				matrix_s[i][j] = tmp_max + match_score;
			else
				matrix_s[i][j] = tmp_max + mismatch_score;

		}
	}

	//       for (i=0;i<length_s ; i++) 
	//       {
	//           for (j=0; j<length_t; j++) 
	//               cout << matrix_s[i][j] << " ";
	//           cout << "\n";
	//       }


		   // Find the maximum score in the first row or column

	int i_index = 0, j_index = 0;
	float max_score = matrix_s[0][0];
	for (j = 0; j < length_t; j++)
	{
		if (matrix_s[0][j] >= max_score)
		{
			max_score = matrix_s[0][j];
			j_index = j;
		}
	}
	for (i = 0; i < length_s; i++)
	{
		if (matrix_s[i][0] >= max_score)
		{
			max_score = matrix_s[i][0];
			i_index = i;
			j_index = 0;
		}
	}

	return init += max_score;

}

//.....................................................................................
// NEEDLEMAN & WUNSCH ALIGNMENT
//.....................................................................................

/**#: [Description="Pairwise alignment for two complete sequences.
					 The score for any pair of elements is
			match_score (or mismatch_score) - total gap penalty.
			The total gap penalty = gap_penalty + (length_of_gap * gap_length_multiplier).
					 The aligned sequence elements are placed alternately in
					 the output container i.e with the odd elements representing
					 the first sequence."]
	  [Restrictions="The containers holding the sequences must have at least
					 forward iterators and an equality operator (==) should be defined
					 for pairs of elements. The total gap penalty =
					 gap_penalty + (length_of_gap * gap_length_multiplier).
					 The result container should either have no fixed size or
					 be created at least twice the size of the longer sequence."] */


template <class InputIterator, class OutputIterator>
OutputIterator needleman_wunsch_alignment(InputIterator first1, InputIterator last1,
	InputIterator first2, InputIterator last2,
	float match_score, float mismatch_score,
	float gap_penalty, float gap_length_multiplier,
	OutputIterator result)
{

	// get the lengths of the sequences + 1
	int length_s = (last1 - first1) + 1;
	int length_t = (last2 - first2) + 1;

	// Create a score matrix
	matrix<float> matrix_s(length_s, length_t);

	// Calculate matrix values for the last row and column.

	// The scoring matrix is padded at the end by a blank

	matrix_s[length_s - 1][length_t - 1] = 0.0;

	int i, j;
	for (i = (length_s - 2); i > -1; i--)
	{
		matrix_s[i][length_t - 1] = (i - length_s + 1) * gap_length_multiplier - gap_penalty;
	}
	for (j = (length_t - 2); j > -1; j--)
	{
		matrix_s[length_s - 1][j] = (j - length_t + 1) * gap_length_multiplier - gap_penalty;
	}


	// Calculate the other values in the matrix.

	float cval, rval, dval, tmp_max;
	int ii, jj;
	for (i = (length_s - 2); i > -1; i--)
	{
		for (j = (length_t - 2); j > -1; j--)
		{
			// We now calculate the maximum of the row cross-ways and to the right, 
			// and the column cross-ways and down.

			// Get the cross-value.
			cval = matrix_s[i + 1][j + 1];

			// Get the down-value.
			dval = matrix_s[i + 1][j + 1] - gap_length_multiplier;
			for (ii = i + 2; ii < length_s; ii++)
			{
				if ((matrix_s[ii][j + 1] - (ii - i - 1) * gap_length_multiplier) > dval)
					dval = matrix_s[ii][j + 1] - (ii - i - 1) * gap_length_multiplier;
			}


			// Get the right-value       
			rval = matrix_s[i + 1][j + 1] - gap_length_multiplier;
			for (jj = j + 2; jj < length_t; jj++)
			{
				if ((matrix_s[i + 1][jj] - (jj - j - 1) * gap_length_multiplier) > rval)
					rval = matrix_s[i + 1][jj] - (jj - j - 1) * gap_length_multiplier;
			}

			tmp_max = max(dval, rval);
			tmp_max = max((tmp_max - gap_penalty), cval);
			if (*(first1 + i) == *(first2 + j))
				matrix_s[i][j] = tmp_max + match_score;
			else
				matrix_s[i][j] = tmp_max + mismatch_score;

		}
	}


	// Trace path back through the matrix: choosing successive maxima of the 
	// right, down and diagonal elements
	// In event of tie choose diagonal in preference to right or down
	// and down in preference to right i.e. choose shortest path through matrix 
	// and add gaps to the sequence in preference to the target. 
	// This choice mechanism is somewhat arbitrary from a theoretical standpoint 
	// yet seems sensible from a biological one. There may be alternative start points
	// where scores of two alignments match here we choose such that the target
	// has minimum gaps


	// Find the maximum score in the first row or column

	int i_index = 0, j_index = 0;
	float max_score = matrix_s[0][0];
	for (j = 0; j < length_t; j++)
	{
		if (matrix_s[0][j] >= max_score)
		{
			max_score = matrix_s[0][j];
			j_index = j;
		}
	}
	for (i = 0; i < length_s; i++)
	{
		if (matrix_s[i][0] >= max_score)
		{
			max_score = matrix_s[i][0];
			i_index = i;
			j_index = 0;
		}
	}

	// the starting sequence element is (i_index,j_index)

	i = 0; j = 0;
	typename iterator_traits<OutputIterator>::value_type load_element_s;
	typename iterator_traits<OutputIterator>::value_type load_element_t;
	const typename iterator_traits<OutputIterator>::value_type blank = ' ';

	if (i_index == 0)
	{
		for (j = 0; j < j_index; j++)
		{
			*result = blank; ++result;
			*result = *first2; ++result; first2++;
		}
		load_element_s = *first1; ++first1;
		load_element_t = *first2; ++first2;
		j = j_index;
	}
	else
	{
		for (i = 0; i < i_index; i++)
		{
			*result = *first1; ++result; first1++;
			*result = blank; ++result;
		}
		load_element_s = *first1; ++first1;
		load_element_t = *first2; ++first2;
		i = i_index;
	}


	while (i < (length_s - 1) && j < (length_t - 1))
	{
		if ((matrix_s[i + 1][j + 1] >= matrix_s[i][j + 1])
			&& (matrix_s[i + 1][j + 1] >= matrix_s[i + 1][j]))
		{
			*result = load_element_s; ++result;       // diagonal   
			*result = load_element_t; ++result;
			
			if (first1 != last1)//R.Alcraft end of iterator error fixed (24/2/20)
			{
				load_element_s = *first1; ++first1;				
			}
			else
			{
				load_element_s = blank;
			}
			if (first2 != last2)//R.Alcraft end of iterator error fixed (24/2/20)
			{
				load_element_t = *first2; ++first2;				
				
			}
			else
			{
				load_element_t = blank;
			}
			++i;
			++j;
		}
		else
		{
			if (matrix_s[i][j + 1] > matrix_s[i + 1][j])       // right > down
			{
				*result = load_element_s; ++result;
				*result = load_element_t; ++result;
				load_element_s = blank;                 // right	
				load_element_t = *first2; ++first2;
				++j;

			}
			else
			{
				*result = load_element_s; ++result;
				*result = load_element_t; ++result;
				load_element_s = *first1; ++first1;       // down	
				load_element_t = blank;
				++i;
			}
		}
	}

	// pad out the ends of sequences where necessary     

	for (i; i < (length_s - 1); ++i)
	{
		*result = load_element_s; ++result;
		*result = load_element_t; ++result;
		if (first1 != last1)//R.Alcraft end of iterator error fixed (24/2/20)
		{
			load_element_s = *first1; ++first1;
		}
		else
		{
			load_element_s = blank;
		}						
		load_element_t = blank;
	}

	for (j; j < (length_t - 1); ++j)
	{
		*result = load_element_s; ++result;
		*result = load_element_t; ++result;
		load_element_s = blank;
		if (first2 != last2)//R.Alcraft end of iterator error fixed (24/2/20)
		{
			load_element_t = *first2; ++first2;
		}
		else
		{
			load_element_t = blank;
		}		
	}

	return result;
}

_BTL_END_NAMESPACE

#endif   // End of btl_sequence_algorithms.h 

// GENERALIZE FURTHER AND SWITCH TO GOTOH ALGORITHM FOR CONSTRUCTING THE 
// SCORING MATRIX FOR CONCAVE PENALTIES ASAP

