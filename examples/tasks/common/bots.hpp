/*
 * BLAZE Task testing framework
 *
 * This file is part of the Blaze testing and sample codes. It defines
 * inlineable alternatives to functions and macros defined in the Barcelona
 * OpenMP Tasks Suite's bots.h file.
 *
 * (c) Peter Pirkelbauer, 2018
 * UAB - University of Alabama at Birmingham
 */

/****************************************************************************/
/*  This program is free software; you can redistribute it and/or modify    */
/*  it under the terms of the GNU General Public License as published by    */
/*  the Free Software Foundation; either version 2 of the License, or       */
/*  (at your option) any later version.                                     */
/*                                                                          */
/*  This program is distributed in the hope that it will be useful,         */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/*  GNU General Public License for more details.                            */
/*                                                                          */
/*  You should have received a copy of the GNU General Public License       */
/*  along with this program; if not, write to the Free Software Foundation, */
/*  Inc., 51 Franklin Street, Fifth Floor, Boston, MA, 02110-1301,USA       */
/****************************************************************************/

#ifndef BOTS_HPP

#define BOTS_HPP 1

#include <iostream>

#define FALSE 0
#define TRUE 1

#define BOTS_RESULT_NA 0
#define BOTS_RESULT_SUCCESSFUL 1
#define BOTS_RESULT_UNSUCCESSFUL 2
#define BOTS_RESULT_NOT_REQUESTED 3

inline
void bots_error(int error, const char* message)
{
  std::cerr << "error: " << error << " " << message << std::endl;
  exit(1);
}

inline
void bots_warning(int warning, const char* message)
{
  std::cerr << "warn: " << warning << " " << message << std::endl;
}

inline
void bots_message(const char*, ...) {}

inline
void bots_debug(const char*, ...) {}

#endif /* BOTS_HPP 1 */
