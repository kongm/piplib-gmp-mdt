/******************************************************************************
 *                     PIP : Parametric Integer Programming                   *
 ******************************************************************************
 *                                  type.h                                    *
 ******************************************************************************
 *                                                                            *
 * Copyright Paul Feautrier, 1988, 1993, 1994, 1996, 2002                     *
 *                                                                            *
 * This is free software; you can redistribute it and/or modify it under the  *
 * terms of the GNU General Public License as published by the Free Software  *
 * Foundation; either version 2 of the License, or (at your option) any later *
 * version.							              *
 *                                                                            *
 * This software is distributed in the hope that it will be useful, but       *
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY *
 * or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License   *
 * for more details.							      *
 *                                                                            *
 * You should have received a copy of the GNU General Public License along    *
 * with software; if not, write to the Free Software Foundation, Inc.,        *
 * 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA                     *
 *                                                                            *
 * Written by Paul Feautrier                                                  *
 *                                                                            *
 ******************************************************************************/

#ifndef TYPEMP_H
#define TYPEMP_H
#if defined(__cplusplus)
extern "C" 
  {
#endif 

/* Modified by Serge Torres to handle very big problems (since 1.3.4 we can put
 * any value we want: sol_space is allocated dynamically), but it is left by
 * default to 4096 because of time/space reasons for most people.
 * #define SOL_SIZE 67108864
 */
#ifndef SOL_SIZE    
#define SOL_SIZE (4096*4)
    //#define SOL_SIZE 67108864
#endif
    
extern EntierMP pipmp_UN;
extern EntierMP pipmp_ZERO;
    
#ifndef Pip_True
#define Pip_True 1
#define Pip_False 0
#endif

#ifdef TC
#define DEBUG 8
#endif
    
#ifndef Q    
#define Q if(cross_product>=limit)
#endif
    
#ifndef  MAXCOL
#define MAXCOL 512
#define MAXPARM 150
    
#undef MAXCOL
#define MAXCOL 4096
    
    
#endif
    
#if defined(__cplusplus)
  }
#endif 
#endif /* define _H */
