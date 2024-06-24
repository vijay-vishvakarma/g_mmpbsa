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
 *  Data types used internally in the nbnxn_kokkos module.
 *
 *  \author Sikandar Y. Mashayak <symashayak@gmail.com>
 *  \ingroup module_mdlib
 */

#ifndef GMX_KOKKOS_TYPE_H
#define GMX_KOKKOS_TYPE_H

#include "config.h"

#include <impl/Kokkos_Timer.hpp>
#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>
#include <Kokkos_Vectorization.hpp>

#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/nbnxn_pairlist.h"

// set GMXHostype and GMXDeviceType from Kokkos Default Types
typedef Kokkos::DefaultExecutionSpace GMXDeviceType;
typedef Kokkos::HostSpace::execution_space GMXHostType;

// set ExecutionSpace stuct with variable "space"

// template<class Device>
// struct ExecutionSpaceFromDevice;

// template<>
// struct ExecutionSpaceFromDevice<LMPHostType> {
//   static const LAMMPS_NS::ExecutionSpace space = LAMMPS_NS::Host;
// };
// #ifdef KOKKOS_HAVE_CUDA
// template<>
// struct ExecutionSpaceFromDevice<Kokkos::Cuda> {
//   static const LAMMPS_NS::ExecutionSpace space = LAMMPS_NS::Device;
// };
// #endif

// GROMACS types


template <class DeviceType>
struct ArrayTypes;

template <>
struct ArrayTypes<GMXDeviceType> {

    // scalar types

    typedef Kokkos::
        DualView<int, GMXDeviceType::array_layout, GMXDeviceType> tdual_int_scalar;
    typedef tdual_int_scalar::t_dev t_int_scalar;
    typedef tdual_int_scalar::t_dev_const t_int_scalar_const;
    typedef tdual_int_scalar::t_dev_um t_int_scalar_um;
    typedef tdual_int_scalar::t_dev_const_um t_int_scalar_const_um;

    // generic array types

    // 1d unmanaged view for int array n with right layout
    typedef Kokkos::View<int*, GMXDeviceType::array_layout, GMXDeviceType, Kokkos::MemoryTraits<Kokkos::Unmanaged>> t_un_int_1d;

    // 1d unmanaged view for real array n with right layout
    typedef Kokkos::View<real*, GMXDeviceType::array_layout, GMXDeviceType, Kokkos::MemoryTraits<Kokkos::Unmanaged>> t_un_real_1d;

    // 2d (n*3) unmanaged view for real array n with right layout
    typedef Kokkos::View<real*[3], GMXDeviceType::array_layout, GMXDeviceType, Kokkos::MemoryTraits<Kokkos::Unmanaged>> t_un_real_1d3;

    // 1d unmanaged and atomic view for real array n with right layout
    typedef Kokkos::View<real*, GMXDeviceType::array_layout, GMXDeviceType, Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::Atomic>> t_un_at_real_1d;

    // pairlist related views

    // i-cluster list: unmanaged view on host
    typedef Kokkos::View<nbnxn_ci_t*, GMXDeviceType::array_layout, GMXDeviceType, Kokkos::MemoryTraits<Kokkos::Unmanaged>> t_un_ci_1d;

    // j-cluster list: unmanaged view on host
    typedef Kokkos::View<nbnxn_cj_t*, GMXDeviceType::array_layout, GMXDeviceType, Kokkos::MemoryTraits<Kokkos::Unmanaged>> t_un_cj_1d;

};

//default Gromacs Types
typedef struct ArrayTypes<GMXDeviceType> DAT;
typedef struct ArrayTypes<GMXHostType> HAT;

template<class DeviceType>
struct MemsetZeroFunctor {
    typedef DeviceType  device_type ;
    void* ptr;
    KOKKOS_INLINE_FUNCTION void operator()(const int i) const {
        ((int*)ptr)[i] = 0;
    }
};

#endif /* GMX_KOKKOS_TYPE_H */
