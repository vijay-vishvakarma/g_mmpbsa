/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015, by the GROMACS development team, led by Mark
 * Abraham, David van der Spoel, Berk Hess, and Erik Lindahl, and
 * including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */

/*! \internal \file
 *  \brief
 *  Kokkos versions of create/grow/destroy multi-dimensional arrays
 *
 *  \author Sikandar Y. Mashayak <symashayak@gmail.com>
 *  \ingroup module_mdlib
 */

#include <string>

/* ----------------------------------------------------------------------
   create a 1d array
------------------------------------------------------------------------- */

template <typename TYPE>
TYPE create_kokkos(TYPE &data, typename TYPE::value_type *&array, 
                   int n1, const char *name)
{
  data = TYPE(name,n1);
  array = data.h_view.ptr_on_device();
  return data;
}

// unmanaged 1-d view
template <typename TYPE>
TYPE create_kokkos(TYPE &data, typename TYPE::value_type *&array, 
                   int n1)
{
  data = TYPE(array,n1);
  return data;
}



template <typename TYPE, typename HTYPE>
  TYPE create_kokkos(TYPE &data, HTYPE &h_data, 
                     typename TYPE::value_type *&array, int n1, 
                     const char *name)
{
  data = TYPE(std::string(name),n1);
#ifndef KOKKOS_USE_CUDA_UVM
  h_data = Kokkos::create_mirror_view(data);
#else
  h_data = data;
#endif
  array = h_data.ptr_on_device();
  return data;
}


template <typename TYPE, typename HTYPE>
  TYPE create_kokkos(TYPE &data, HTYPE &h_data,
                     int n1, const char *name)
{
  data = TYPE(std::string(name),n1);
#ifndef KOKKOS_USE_CUDA_UVM
  h_data = Kokkos::create_mirror_view(data);
#else
  h_data = data;
#endif
  return data;
}

/* ----------------------------------------------------------------------
   grow or shrink 1st dim of a 1d array
   last dim must stay the same
------------------------------------------------------------------------- */

template <typename TYPE>
TYPE grow_kokkos(TYPE &data, typename TYPE::value_type *&array, 
                 int n1, const char *name)
{
  if (array == NULL) return create_kokkos(data,array,n1,name);
  
  data.resize(n1);
  array = data.h_view.ptr_on_device();
  return data;
}

template <typename TYPE>
void destroy_kokkos(TYPE data, typename TYPE::value_type* &array)
{
  if (array == NULL) return;
  data = TYPE();
  array = NULL;
}

/* ----------------------------------------------------------------------
   create a 2d array
------------------------------------------------------------------------- */

template <typename TYPE>
TYPE destroy_kokkos(TYPE &data)
{
  /*if(data.ptr_on_device()!=NULL)
    free(data.ptr_on_device());*/
  data = TYPE();
  return data;
}

template <typename TYPE>
TYPE create_kokkos(TYPE &data, int n1, const char *name)
{
  /*typename TYPE::non_const_value_type* ptr = (typename TYPE::non_const_value_type*)
    malloc(n1*sizeof(typename TYPE::non_const_value_type)*4);*/
  data = TYPE(std::string(name),n1);
  return data;
}

template <typename TYPE>
TYPE create_kokkos(TYPE &data, int n1, int n2, const char *name)
{
  /*typename TYPE::non_const_value_type* ptr = (typename TYPE::non_const_value_type*)
    malloc(n1*n2*sizeof(typename TYPE::non_const_value_type)*4);*/
  data = TYPE(std::string(name),n1,n2);
  return data;
}

template <typename TYPE>
TYPE create_kokkos(TYPE &data, int n1, int n2, int n3 ,const char *name)
{
  /*typename TYPE::non_const_value_type* ptr = (typename TYPE::non_const_value_type*)
    malloc(n1*n2*n3*sizeof(typename TYPE::non_const_value_type)*4);*/
  data = TYPE(std::string(name),n1,n2,n3);
  return data;
}

template <typename TYPE>
TYPE create_kokkos(TYPE &data, int n1, int n2, int n3, int n4 ,const char *name)
{
  /*typename TYPE::non_const_value_type* ptr = (typename TYPE::non_const_value_type*)
    malloc(n1*n2*n3*n4*sizeof(typename TYPE::non_const_value_type)*4);*/
  data = TYPE(std::string(name),n1,n2,n3,n4);
  return data;
}

template <typename TYPE>
TYPE create_kokkos(TYPE &data, int n1, int n2, int n3, int n4, int n5 ,const char *name)
{
  /*typename TYPE::non_const_value_type* ptr = (typename TYPE::non_const_value_type*)
    malloc(n1*n2*n3*n4*n5*sizeof(typename TYPE::non_const_value_type)*4);*/
  data = TYPE(std::string(name),n1,n2,n3,n4,n5);
  return data;
}

template <typename TYPE>
TYPE create_kokkos(TYPE &data, int n1, int n2, int n3, int n4, int n5 , int n6 ,const char *name)
{
  /*typename TYPE::non_const_value_type* ptr = (typename TYPE::non_const_value_type*)
    malloc(n1*n2*n3*n4*n5*n6*sizeof(typename TYPE::non_const_value_type)*4);*/
  data = TYPE(std::string(name) ,n1,n2,n3,n4,n5,n6);
  return data;
}



template <typename TYPE, typename HTYPE>
  TYPE create_kokkos(TYPE &data, HTYPE &h_data, int n1, int n2, 
                     const char *name)
{
  data = TYPE(std::string(name),n1,n2);
#ifndef KOKKOS_USE_CUDA_UVM
  h_data = Kokkos::create_mirror_view(data);
#else
  h_data = data;
#endif
  return data;
}


template <typename TYPE>
TYPE create_kokkos(TYPE &data, typename TYPE::value_type **&array, 
                   int n1, int n2, const char *name)
{
  data = TYPE(std::string(name),n1,n2);
  bigint nbytes = ((bigint) sizeof(typename TYPE::value_type *)) * n1;
  array = (typename TYPE::value_type **) smalloc(nbytes,name);
  
  bigint n = 0;
  for (int i = 0; i < n1; i++) {
    array[i] = &data.h_view(i,0);
    n += n2;
  }
  return data;
}

// unmanaged 2-d view
template <typename TYPE>
TYPE create_kokkos(TYPE &data, typename TYPE::value_type **&array, 
                   int n1, int n2)
{
  data = TYPE(&array[0][0],n1,n2);
  return data;
}

template <typename TYPE, typename HTYPE>
  TYPE create_kokkos(TYPE &data, HTYPE &h_data, 
                     typename TYPE::value_type **&array, int n1, int n2, 
                     const char *name)
{
  data = TYPE(std::string(name),n1,n2);
#ifndef KOKKOS_USE_CUDA_UVM
  h_data = Kokkos::create_mirror_view(data);
#else
  h_data = data;
#endif
  bigint nbytes = ((bigint) sizeof(typename TYPE::value_type *)) * n1;
  array = (typename TYPE::value_type **) smalloc(nbytes,name);
  
  bigint n = 0;
  for (int i = 0; i < n1; i++) {
    array[i] = &h_data(i,0);
    n += n2;
  }
  return data;
}

/* ----------------------------------------------------------------------
   grow or shrink 1st dim of a 2d array
   last dim must stay the same
------------------------------------------------------------------------- */

template <typename TYPE>
TYPE grow_kokkos(TYPE &data, typename TYPE::value_type **&array, 
                 int n1, int n2, const char *name)
{
  if (array == NULL) return create_kokkos(data,array,n1,n2,name);
  data.resize(n1,n2);
  bigint nbytes = ((bigint) sizeof(typename TYPE::value_type *)) * n1;
  array = (typename TYPE::value_type**) srealloc(array,nbytes,name);
  
  for (int i = 0; i < n1; i++)
    array[i] = &data.h_view(i,0);
  
  return data;
}

template <typename TYPE>
TYPE create_kokkos(TYPE &data, typename TYPE::value_type **&array, 
                   int n1, const char *name)
{
  data = TYPE(std::string(name),n1);
  bigint nbytes = ((bigint) sizeof(typename TYPE::value_type *)) * n1;
  array = (typename TYPE::value_type **) smalloc(nbytes,name);
  
  for (int i = 0; i < n1; i++)
    array[i] = &data.h_view(i,0);
  
  return data;
}

template <typename TYPE>
TYPE grow_kokkos(TYPE &data, typename TYPE::value_type **&array, 
                 int n1, const char *name)
{
  if (array == NULL) return create_kokkos(data,array,n1,name);
  
  data.resize(n1);
  
  bigint nbytes = ((bigint) sizeof(typename TYPE::value_type *)) * n1;
  array = (typename TYPE::value_type **) smalloc(nbytes,name);
  
  for (int i = 0; i < n1; i++)
    array[i] = &data.h_view(i,0);
  
  return data;
}

/* ----------------------------------------------------------------------
   destroy a 2d array
------------------------------------------------------------------------- */

template <typename TYPE>
void destroy_kokkos(TYPE data, typename TYPE::value_type** &array)
{
  if (array == NULL) return;
  data = TYPE();
  sfree(array);
  array = NULL;
}
