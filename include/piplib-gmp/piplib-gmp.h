/******************************************************************************
 *                     PIP : Parametric Integer Programming                   *
 ******************************************************************************
 *                                 piplib.h                                   *
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
 * Written by Cedric Bastoul                                                  *
 *                                                                            *
 ******************************************************************************/

/* Premiere version du 18 septembre 2002. */

#ifndef PIPLIBMP_H
# define PIPLIBMP_H


#ifndef PIPMP_LINEAR_VALUE_IS_MP
# define PIPMP_LINEAR_VALUE_IS_MP
#endif

# include <gmp.h>
# define EntierMP   mpz_t
# define PIPMP_FORMAT   "%d"
# define GMP_INPUT_FORMAT   "%lZd"

# if ! defined(MPZVALUE_TO_INT)
#  define MPZVALUE_TO_INT(val) ((int)mpz_get_si(val))
# endif


#if ! defined(mpzvalue_addto)
#define mpzvalue_addto(ref,val1,val2)     	(mpz_add((ref),(val1),(val2)))
#endif
#if ! defined(mpzvalue_assign)
#define mpzvalue_assign(v1,v2)	    	(mpz_set((v1),(v2)))
#endif
#if ! defined(mpzvalue_clear)
#define mpzvalue_clear(val)       		(mpz_clear((val)))
#endif
#if ! defined(mpzvalue_divexact)
#define mpzvalue_divexact(d,v1,v2)	    	(mpz_divexact((d),(v1),(v2)))
#endif
#if ! defined(mpzvalue_gcd)
#define mpzvalue_gcd(g,v1,v2)	    	(mpz_gcd((g),(v1),(v2)))
#endif
#if ! defined(mpzvalue_inti)
#define mpzvalue_init(val)        	    	(mpz_init((val)))
#endif
#if ! defined(mpzvalue_init_set)
#define mpzvalue_init_set(v1,v2)	    	(mpz_init_set((v1),(v2)))
#endif
#if ! defined(mpzvalue_oppose)
#define mpzvalue_oppose(ref,val)       	(mpz_neg((ref),(val)))
#endif
#if ! defined(mpzvalue_set_si)
#define mpzvalue_set_si(val,i)    		(mpz_set_si((val),(i)))
#endif
#if ! defined(mpzvalue_substract)
#define mpzvalue_subtract(ref,val1,val2) 	(mpz_sub((ref),(val1),(val2)))
#endif
#if ! defined(mpzvalue_eq)
#define mpzvalue_eq(v1,v2) 	    	(mpz_cmp((v1),(v2)) == 0)
#endif
#if ! defined(mpzvalue_ne)
#define mpzvalue_ne(v1,v2) 	    	(mpz_cmp((v1),(v2)) != 0)
#endif
#if ! defined(mpzvalue_notzero_p)
#define mpzvalue_notzero_p(val)        	(mpz_sgn(val) != 0)
#endif


#if defined(__cplusplus)
extern "C"
  {
#endif

# include <piplib-gmp/type.h>
# include <piplib-gmp/sol.h>
# include <piplib-gmp/tab.h>
# include <piplib-gmp/funcall.h>


/* Structure PipMPMatrix :
 * Structure de matrice au format PolyLib. Le premier element d'une ligne
 * indique quand il vaut 1 que la ligne decrit une inequation de la forme
 * p(x)>=0 et quand il vaut 0, que la ligne decrit une egalite de la forme
 * p(x)=0. Le dernier element de chaque ligne correspond au coefficient
 * constant.
 */
struct pipmatrixmp
{ unsigned NbRows, NbColumns ;
  EntierMP **p ;
  EntierMP *p_Init ;
  int p_Init_size;	        /* Only for PolyLib compatibility under MP
                                 * version: PolyLib makes sometimes
				 * overestimates on the size of the matrices,
				 * in order to go faster. Thus
				 * NbRows*NbColumns is not the number of
				 * allocated elements. With MP version, we
				 * have to think to mpz_clear() all the
				 * initialized elements before freing, then
				 * we need to know the number of allocated
				 * elements: p_Init_size.
				 */
} ;
typedef struct pipmatrixmp PipMPMatrix ;


/* Structure PipMPVector :
 * Cette structure contient un Vector de 'nb_elements' la ieme composante de
 * ce vecteur vaut the_vector[i]/the_deno[i].
 */
struct pipvectormp
{ int nb_elements ;             /* Nombre d'elements du vecteur. */
  EntierMP * the_vector ;         /* Numerateurs du vecteur. */
  EntierMP * the_deno ;           /* Denominateurs du vecteur. */
} ;
typedef struct pipvectormp PipMPVector ;


/* Structure PipMPNewparm :
 * Liste chainee de Newparm, les informations d'un newparm etant son rang, un
 * vecteur de coefficients et un denominateur. Le newparm est egal a la division
 * du vecteur par le denominateur.
 */
struct pipnewparmmp
{ int rank ;                    /* Rang du 'newparm'. */
  PipMPVector * vector ;          /* Le vector decrivant le newparm. */
  EntierMP deno ;                 /* Denominateur du 'newparm'. */
  struct pipnewparmmp * next ;    /* Pointeur vers le newparm suivant. */
} ;
typedef struct pipnewparmmp PipMPNewparm ;


/* Structure PipMPList :
 * Liste chainee de Vector.
 */
struct piplistmp
{ PipMPVector * vector ;          /* Le vector contenant la partie de solution. */
  struct piplistmp * next ;       /* Pointeur vers l'element suivant. */
} ;
typedef struct piplistmp PipMPList ;


/* Structure pipquast :
 * Arbre binaire. Conformement a la grammaire de sortie (voir mode d'emploi), un
 * noeud de l'arbre des solutions debute par une liste de 'newparm'. Il continue
 * ensuite soit par une 'list' (alors condition vaut null), soit par un 'if'
 * (alors le champ condition contient la condition).
 */
struct pipquastmp
{ PipMPNewparm * newparm ;        /* Les 'newparm'. */
  PipMPList * list ;              /* La 'list' si pas de 'if'. */
  PipMPVector * condition ;       /* La condition si 'if'. */
  struct pipquastmp * next_then ; /* Noeud si condition et si verifiee. */
  struct pipquastmp * next_else ; /* Noeud si condition et si non verifiee. */
  struct pipquastmp * father ;    /* Pointeur vers le quast pere. */
} ;
typedef struct pipquastmp PipMPQuast ;


/* Structure pipoptions:
 * This structure contains each option that can be set to change the PIP
 * behaviour.
 */
struct pipoptionsmp
{ int Nq ;                      /* 1 if an integer solution is needed,
                                 * 0 otherwise.
				 */
  int Verbose ;                 /* -1 -> absolute silence,
                                 *  0 -> relative silence,
                                 *  1 -> information on cuts when an integer
				 *       solution is needed,
                                 *  2 -> information sur les pivots et les
				 *       déterminants,
                                 *  3 -> information on arrays,
                                 * Each option include the preceding.
				 */
  int Simplify ;                /* Set to 1 to eliminate some trivial
                                 * solutions, 0 otherwise.
				 */
  int Deepest_cut ;             /* Set to 1 to include deepest cut
                                 * algorithm.
				 */
  int Maximize;                 /* Set to 1 if maximum is needed. */
  int Urs_parms;             	/* -1 -> all parameters may be negative
				 *  0 -> all parameters are non-negative
				 */
  int Urs_unknowns;             /* -1 -> all unknowns may be negative
				 *  0 -> all unknowns are non-negative
				 */
} ;
typedef struct pipoptionsmp PipMPOptions ;


/* Prototypes des fonctions d'affichages des structures de la PipMPLib. */
void pipmp_matrix_print(FILE *, PipMPMatrix *) ;
void pipmp_vector_print(FILE *, PipMPVector *) ;
void pipmp_newparm_print(FILE * foo, PipMPNewparm *, int indent) ;
void pipmp_list_print(FILE * foo, PipMPList *, int indent) ;
void pipmp_quast_print(FILE *, PipMPQuast *, int) ;


/* Prototypes des fonctions de liberation memoire des structures de la PipMPLib.*/
void pipmp_matrix_free(PipMPMatrix *) ;
void pipmp_vector_free(PipMPVector *) ;
void pipmp_newparm_free(PipMPNewparm *) ;
void pipmp_list_free(PipMPList *) ;
void pipmp_quast_free(PipMPQuast *) ;
void pipmp_options_free(PipMPOptions *) ;


/* Prototypes des fonctions d'acquisition de matrices de contraintes et
 * options.
 */
PipMPMatrix * pipmp_matrix_alloc(unsigned, unsigned) ;
PipMPMatrix * pipmp_matrix_read(FILE *) ;
PipMPOptions * pipmp_options_init(void) ;


/* initialization of pip library */
void pipmp_init();
void pipmp_close();


/* Prototype de la fonction de resolution :
 * pipmp_solve resoud le probleme qu'on lui passe en parametre, suivant les
 * options elles aussi en parametre. Elle renvoie la solution sous forme
 * d'un arbre de PipMPQuast. Parametres :
 * - probleme :
 * 1 PipMPMatrix  : systeme des inequations definissant le domaine des inconnues,
 * 2 PipMPMatrix  : systeme des inequations satisfaites par les parametres,
 * 3 int        : column rank of the bignum, or negative value if there
 *                is no big parameter.
 * 4 PipMPOptions : options for PIP.
 */
PipMPQuast * pipmp_solve(PipMPMatrix *, PipMPMatrix *, int, PipMPOptions *) ;

/* Ced : ajouts specifiques a la PipMPLib pour funcall. */
TableauMP * tabmp_Matrix2TableauMP(PipMPMatrix *, int, int, int, int, int, int);
    //TableauMP * tabmp_Matrix2TableauMPMax(PipMPMatrix *, int, int, int, int) ;
#ifndef SOL_SHIFT 
# define SOL_SHIFT		(1 << 0)    /* Shift solution over -bigparam */
# define SOL_NEGATE		(1 << 1)    /* Negate solution */
# define SOL_REMOVE		(1 << 2)    /* Remove big parameter */
# define SOL_MAX			(SOL_SHIFT | SOL_NEGATE)
    /* Maximum was computed */
#endif     
PipMPQuast *solmp_quast_edit(int *i, PipMPQuast *father, int Bg, int Urs_p, int flags);

#if defined(__cplusplus)
  }
#endif
#endif /* define PIPLIBMP_H */
