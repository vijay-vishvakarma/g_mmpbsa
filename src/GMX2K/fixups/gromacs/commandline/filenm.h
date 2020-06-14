/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017,2018,2019, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
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
/*! \file
 * \brief
 * Declares t_filenm for old-style command-line parsing of file name options.
 *
 * \inpublicapi
 * \ingroup module_commandline
 */
#ifndef GMX_COMMANDLINE_FILENM_H
#define GMX_COMMANDLINE_FILENM_H

#include <string>
#include <vector>

//#include "gromacs/compat/string_view.h"
#include "gromacs/fileio/filetypes.h"
//#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"


//! \addtogroup module_commandline
//! \{

/*! \brief
 * File name option definition for C code.
 *
 * \inpublicapi
 */
struct t_filenm
{
    int                      ftp; //!< File type, see enum in filetypes.h
    const char*              opt; //!< Command line option, can be nullptr in which case the commandline module, including all opt2??? functions below, will use the default option for the file type
    const char*              fn; //!< File name (as set in source code), can be nullptr in which case the commandline module will use the default file name for the file type
    unsigned long            flag;      //!< Flag for all kinds of info (see defs)
    std::vector<std::string> filenames; //!< File names
};

//! Whether a file name option is set.
#define ffSET 1 << 0
//! Whether a file name option specifies an input file.
#define ffREAD 1 << 1
//! Whether a file name option specifies an output file.
#define ffWRITE 1 << 2
//! Whether a file name option specifies an optional file.
#define ffOPT 1 << 3
//! Whether a file name option specifies a library file.
#define ffLIB 1 << 4
//! Whether a file name option accepts multiple file names.
#define ffMULT 1 << 5
//! Whether an input file name option accepts non-existent files.
#define ffALLOW_MISSING 1 << 6
//! Convenience flag for an input/output file.
#define ffRW (ffREAD | ffWRITE)
//! Convenience flag for an optional input file.
#define ffOPTRD (ffREAD | ffOPT)
//! Convenience flag for an optional output file.
#define ffOPTWR (ffWRITE | ffOPT)
//! Convenience flag for an optional input/output file.
#define ffOPTRW (ffRW | ffOPT)
//! Convenience flag for a library input file.
#define ffLIBRD (ffREAD | ffLIB)
//! Convenience flag for an optional library input file.
#define ffLIBOPTRD (ffOPTRD | ffLIB)
//! Convenience flag for an input file that accepts multiple files.
#define ffRDMULT (ffREAD | ffMULT)
//! Convenience flag for an optional input file that accepts multiple files.
#define ffOPTRDMULT (ffRDMULT | ffOPT)
//! Convenience flag for an output file that accepts multiple files.
#define ffWRMULT (ffWRITE | ffMULT)
//! Convenience flag for an optional output file that accepts multiple files.
#define ffOPTWRMULT (ffWRMULT | ffOPT)

/*! \brief
 * Returns the filename belonging to cmd-line option opt, or NULL when
 * no such option.
 */
const char* opt2fn(const char* opt, int nfile, const t_filenm fnm[]);

//! Returns the first file name with type ftp, or NULL when none found.
const char* ftp2fn(int ftp, int nfile, const t_filenm fnm[]);

/*! \brief
 * Returns the file name belonging top cmd-line option opt, or NULL when
 * no such option.
 *
 * Also return NULL when opt is optional and option is not set.
 */
const char* opt2fn_null(const char* opt, int nfile, const t_filenm fnm[]);

/*! \brief
 * Returns the first file name with type ftp, or NULL when none found.
 *
 * Also return NULL when ftp is optional and option is not set.
 */
const char* ftp2fn_null(int ftp, int nfile, const t_filenm fnm[]);

#endif
