/******************************************************************************
 *                     PIP : Parametric Integer Programming                   *
 ******************************************************************************
 *                                 funcall.h                                  *
 ******************************************************************************
 *                                                                            *
 * Copyright Paul Feautrier, 1988-2005                                        *
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
 *****************************************************************************/

#ifndef FUNCALLMP_H
#define FUNCALLMP_H
#if defined(__cplusplus)
extern "C" 
  {
#endif 

int pipmp_traiter(TableauMP *, TableauMP *, int, int, int, int, int, int);
int pipmp_integrer(TableauMP **, TableauMP **, int *, int *, int *, int *);

int solmp_hwm(void);
void solmp_simplify(int);
int is_not_Nil(int);
int solmp_edit(FILE *, int);
void tabmp_reset(struct pipmp_high_water_mark);
void solmp_reset(int);
struct pipmp_high_water_mark tabmp_hwm(void);
TableauMP *tabmp_get(FILE *, int,int,int);
void solmp_init(void);
void solmp_close(void);
void tabmp_init(void);
void tabmp_close(void);
void solmp_if(void);
void solmp_forme(int);
void solmp_val(EntierMP, EntierMP);
void solmp_nil(void);
void solmp_error(int);
TableauMP * tabmp_alloc(int, int, int);
void solmp_list(int);
void tabmp_display(TableauMP *, FILE *);
TableauMP * pipmp_expanser(TableauMP *, int, int, int, int, int, int);
void solmp_new(int);
void solmp_div(void);

#if defined(__cplusplus)
  }
#endif 
#endif /* define _H */
