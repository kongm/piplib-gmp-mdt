/******************************************************************************
 *                     PIP : Parametric Integer Programming                   *
 ******************************************************************************
 *                                   tab.h                                    *
 ******************************************************************************
 *                                                                            *
 * Copyright Paul Feautrier, 1988, 1993, 1994, 1996, 2002                     *
 *                                                                            *
 * This is free software; you can redistribute it and/or modify it under the  *
 * terms of the GNU General Public License as published by the Free Software  *
 * Foundation; either version 2 of the License, or (at your option) any later *
 * version.                                                                   *
 *                                                                            *
 * This software is distributed in the hope that it will be useful, but       *
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY *
 * or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License   *
 * for more details.                                                          *
 *                                                                            *
 * You should have received a copy of the GNU General Public License along    *
 * with software; if not, write to the Free Software Foundation, Inc.,        *
 * 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA                     *
 *                                                                            *
 * Written by Paul Feautrier and Cedric Bastoul                               *
 *                                                                            *
 *****************************************************************************/


#include <stdio.h>
#include <stdlib.h>

#include <piplib-gmp/piplib-gmp.h>

#ifndef PIPMP_LINEAR_VALUE_IS_MP
# define PIPMP_LINEAR_VALUE_IS_MP
#endif


#define TAB_CHUNK (4096*8)*sizeof(EntierMP)

static char *tabmp_free, *tabmp_top;
static struct pipmp_A *tabmp_base;

extern int pipmp_allocation;
extern long int pipmp_cross_product, pipmp_limit;
static int pipmp_chunk_count;

int pipmp_dgetc(FILE *);
#if defined(PIPMP_LINEAR_VALUE_IS_MP)
int pipmp_dscanf(FILE *, EntierMP);
#else
int pipmp_dscanf(FILE *, EntierMP *);
#endif

extern FILE * pipmp_dump;

void tabmp_init(void)
{
 tabmp_free = malloc(sizeof (struct pipmp_A));
 if(tabmp_free == NULL)
     {fprintf(stderr, "Your computer doesn't have enough memory\n");
      exit(1);
     }
 pipmp_allocation = 1;
 tabmp_top = tabmp_free + sizeof (struct pipmp_A);
 tabmp_base = (struct pipmp_A *)tabmp_free;
 tabmp_free += sizeof(struct pipmp_A);
 tabmp_base->precedent = NULL;
 tabmp_base->bout = tabmp_top;
 tabmp_base->free = tabmp_free;
 pipmp_chunk_count = 1;
}
 
 
void tabmp_close(void)
{
  if (tabmp_base) free(tabmp_base);
}


struct pipmp_high_water_mark tabmp_hwm(void)
{struct pipmp_high_water_mark p;
 p.chunk = pipmp_chunk_count;
 p.top = tabmp_free;
 return p;
}


#if defined(PIPMP_LINEAR_VALUE_IS_MP)
/* the clear_tab routine clears the GMP objects which may be referenced
   in the given TableauMP.
*/
void tabmp_clear(TableauMP *tp)
{
  int i, j;
  /* clear the determinant */
  mpz_clear(tp->determinant);

  for(i=0; i<tp->height; i++){
    /* clear the denominator */
    mpz_clear(Denom(tp, i));
    if((Flag(tp, i) & Unit) == 0)
      for(j=0; j<tp->width; j++)
        mpz_clear(Index(tp,i,j));
  }
}
#endif

void tabmp_reset(struct pipmp_high_water_mark by_the_mark)

{struct pipmp_A *g;
 char *p;
 while(pipmp_chunk_count > by_the_mark.chunk)
     {
      g = tabmp_base->precedent;
      
      #if defined(PIPMP_LINEAR_VALUE_IS_MP)
      /* Before actually freeing the memory, one has to clear the
       * included TableauMPx. If this is not done, the GMP objects
       * referenced in the TableauMPx will be orphaned.
       */

      /* Enumerate the included tableaux. */
      p = (char *)tabmp_base + sizeof(struct pipmp_A);
      while(p < tabmp_base->free){
        TableauMP *pt;
        pt = (TableauMP *) p;
	tabmp_clear(pt);
        p += pt->taille;
      } 
      #endif
      
      free(tabmp_base);
      tabmp_base = g;
      tabmp_top = tabmp_base->bout;
      pipmp_chunk_count--;
     }
 if(pipmp_chunk_count > 0) {
     #if defined(PIPMP_LINEAR_VALUE_IS_MP)
     /* Do not forget to clear the tables in the current chunk above the
        high water mark */
     p = (char *)by_the_mark.top;
     while(p < tabmp_base->free) {
        TableauMP *pt;
        pt = (TableauMP *) p;
        tabmp_clear(pt);
        p += pt->taille;
        } 
     #endif   
     tabmp_free = by_the_mark.top;
     tabmp_base->free = tabmp_free;
     }
 else {
     fprintf(stderr, "Syserr: tabmp_reset : error in memory allocation\n");
     exit(1);
     }
}

TableauMP * tabmp_alloc(int h, int w, int n)

/* h : le nombre de ligne reelles;
   n : le nombre de lignes virtuelles
*/
{
 char *p; TableauMP *tp;
 EntierMP *q;
 unsigned long taille;
 int i, j;
 taille = sizeof(struct pipmp_T) + (h+n-1) * sizeof (struct pipmp_L)
	  + h * w * sizeof (EntierMP);
 if(tabmp_free + taille >= tabmp_top)
     {struct pipmp_A * g;
      unsigned long d;
      d = taille + sizeof(struct pipmp_A);
      if(d < TAB_CHUNK) d = TAB_CHUNK;
      tabmp_free = malloc(d);
      if(tabmp_free == NULL)
	  {printf("Memory overflow\n");
	   exit(23);
	  }
      pipmp_chunk_count++;
      g = (struct pipmp_A *)tabmp_free;
      g->precedent = tabmp_base;
      tabmp_top = tabmp_free + d;
      tabmp_free += sizeof(struct pipmp_A);
      tabmp_base = g;
      g->bout = tabmp_top;
     }
 p = tabmp_free;
 tabmp_free += taille;
 tabmp_base->free = tabmp_free;
 tp = (TableauMP *)p;
 q = (EntierMP *)(p +  sizeof(struct pipmp_T) + (h+n-1) * sizeof (struct pipmp_L));
 #if defined(PIPMP_LINEAR_VALUE_IS_MP)
 mpz_init_set_ui(tp->determinant,1);
 #else
 tp->determinant[0] = (EntierMP) 1;
 tp->l_determinant = 1;
 #endif
 for(i = 0; i<n ; i++){
   tp->row[i].flags = Unit;
   tp->row[i].objet.unit = i;
   #if defined(PIPMP_LINEAR_VALUE_IS_MP)
   mpz_init_set_ui(Denom(tp, i), 1);
   #else
   Denom(tp, i) = pipmp_UN ;
   #endif
 }
 for(i = n; i < (h+n); i++){
   tp->row[i].flags = 0;
   tp->row[i].objet.val = q;
   for(j = 0; j < w; j++)
   #if defined(PIPMP_LINEAR_VALUE_IS_MP)
   mpz_init_set_ui(*q++, 0); /* loop body. */
   mpz_init_set_ui(Denom(tp, i), 0);
   #else
   *q++ = 0;                 /* loop body. */
   Denom(tp, i) = pipmp_ZERO ;
   #endif
 }
 tp->height = h + n; tp->width = w;
 #if defined(PIPMP_LINEAR_VALUE_IS_MP)
 tp->taille = taille ;
 #endif
 
 return(tp);
}

TableauMP * tabmp_get(foo, h, w, n)
FILE * foo;
int h, w, n;
{
 TableauMP *p;
 int i, j, c;
 EntierMP x;
 #if defined(PIPMP_LINEAR_VALUE_IS_MP)
 mpz_init(x);
 #endif
 
 p = tabmp_alloc(h, w, n);
 while((c = pipmp_dgetc(foo)) != EOF)
      if(c == '(')break;
 for(i = n; i<h+n; i++)
     {p->row[i].flags = Unknown;
      #if defined(PIPMP_LINEAR_VALUE_IS_MP)
      mpz_set_ui(Denom(p, i), 1);
      #else
      Denom(p, i) = pipmp_UN;
      #endif
      while((c = pipmp_dgetc(foo)) != EOF)if(c == '[')break;
      for(j = 0; j<w; j++){
        #if defined(PIPMP_LINEAR_VALUE_IS_MP)
	if(pipmp_dscanf(foo, x) < 0)
          return NULL;
        else
	  mpz_set(p->row[i].objet.val[j], x);
        #else
	if(pipmp_dscanf(foo, &x) < 0)
          return NULL;
        else
	  p->row[i].objet.val[j] = x;
        #endif
        }
      } 
      while((c = pipmp_dgetc(foo)) != EOF)if(c == ']')break;
 
 #if defined(PIPMP_LINEAR_VALUE_IS_MP)
 mpz_clear(x);
 #endif
     
 return(p);
}


/* Fonction tabmp_Matrix2TableauMP :
 * Cette fonction effectue la conversion du format de matrice de la polylib
 * vers le format de traitement de Pip. matrix est la matrice a convertir.
 * Nineq est le nombre d'inequations necessaires (dans le format de la
 * polylib, le premier element d'une ligne indique si l'equation decrite
 * est une inequation ou une egalite. Pip ne gere que les inequations. On
 * compte donc le nombre d'inequations total pour reserver la place
 * necessaire, et on scinde toute egalite p(x)=0 en p(x)>=0 et -p(x)>=0).
 * Nv est le nombre de variables dans la premiere serie de variables (c'est
 * a dire que si les premiers coefficients dans les lignes de la matrice
 * sont ceux des inconnues, Nv est le nombre d'inconnues, resp. parametres).
 * n est le nombre de lignes 'virtuelles' contenues dans la matrice (c'est
 * a dire en fait le nombre d'inconnues). Si Shift vaut 0, on va rechercher
 * le minimum lexicographique non-negatif, sinon on recherche le maximum 
 * (Shift = 1) ou bien le minimum tout court (Shift = -1). La fonction
 * met alors en place le bignum s'il n'y est pas deja et prepare les
 * contraintes au calcul du maximum lexicographique.
 * 27 juillet 2001 : Premiere version, Ced.
 * 30 juillet 2001 : Nombreuses modifications. Le calcul du nombre total
 *                   d'inequations (Nineq) se fait a present a l'exterieur.
 *  3 octobre 2001 : Pas mal d'ameliorations.
 * 18 octobre 2003 : Mise en place de la possibilite de calculer le
 *                   maximum lexicographique (parties 'if (Max)').
 */
TableauMP * tabmp_Matrix2TableauMP(matrix, Nineq, Nv, n, Shift, Bg, Urs_parms)
PipMPMatrix * matrix ;
int Nineq, Nv, n, Shift, Bg, Urs_parms;
{ TableauMP * p ;
  unsigned i, j, k, current, new, nb_columns, decal=0, bignum_is_new ;
  int inequality;
  EntierMP bignum;
  
  mpzvalue_init(bignum) ;
  nb_columns = matrix->NbColumns - 1 ;
  /* S'il faut un BigNum et qu'il n'existe pas, on lui reserve sa place. */
  bignum_is_new = Shift && (Bg > (matrix->NbColumns - 2));
  if (bignum_is_new)
    nb_columns++;
  /* Ce sont juste des parametres. */
  if (Bg <= Nv)
    Shift = 0;

  p = tabmp_alloc(Nineq,nb_columns+Urs_parms,n) ;
    
  /* La variable decal sert a prendre en compte les lignes supplementaires
   * issues des egalites.
   */
  for (i = 0; i < matrix->NbRows; i++) {
    current = i + n + decal;
    Flag(p,current) = Unknown ;
    mpzvalue_set_si(Denom(p,current), 1);
    if (Shift)
      mpzvalue_set_si(bignum, 0);
    /* Pour passer l'indicateur d'egalite/inegalite. */
    inequality = mpzvalue_notzero_p(matrix->p[i][0]);
         
    /* Dans le format de la polylib, l'element constant est place en
     * dernier. Dans le format de Pip, il se trouve apres la premiere
     * serie de variables (inconnues ou parametres). On remet donc les
     * choses dans l'ordre de Pip. Ici pour p(x) >= 0.
     */
    for (j=0;j<Nv;j++) {
      if (bignum_is_new && 1+j == Bg)
	continue;
      if (Shift)
	mpzvalue_addto(bignum, bignum, matrix->p[i][1+j]);
      if (Shift > 0)
	mpzvalue_oppose(p->row[current].objet.val[j], matrix->p[i][1+j]);
      else
	mpzvalue_assign(p->row[current].objet.val[j], matrix->p[i][1+j]);
    }
    for (k=j=Nv+1;j<nb_columns;j++) {
	if (bignum_is_new && j == Bg)
	  continue;
	mpzvalue_assign(p->row[current].objet.val[j], matrix->p[i][k++]);
    }
    for (j=0; j < Urs_parms; ++j) {
	int pos = nb_columns - Urs_parms + j;
	if (pos <= Nv)
	    --pos;
	if (pos <= Bg)
	    --pos;
	mpzvalue_oppose(p->row[current].objet.val[nb_columns+j],
		     p->row[current].objet.val[pos]);
    }
    mpzvalue_assign(p->row[current].objet.val[Nv], 
		 matrix->p[i][matrix->NbColumns-1]);
    if (Shift) {
      if (Shift < 0)
	mpzvalue_oppose(bignum, bignum);

      if (bignum_is_new)
	mpzvalue_assign(p->row[current].objet.val[Bg], bignum);
      else
	mpzvalue_addto(p->row[current].objet.val[Bg], 
		    p->row[current].objet.val[Bg], bignum);
    }
    
    /* Et ici lors de l'ajout de -p(x) >= 0 quand on traite une egalite. */
    if (!inequality) {
      decal ++ ;
      new = current + 1 ;
      Flag(p,new)= Unknown ;
      mpzvalue_set_si(Denom(p,new), 1);
      
      for (j=0;j<nb_columns+Urs_parms;j++)
	mpzvalue_oppose(p->row[new].objet.val[j], p->row[current].objet.val[j]);
    }
  }
  mpzvalue_clear(bignum);

  return(p);
}


char *Attr[] = {"Unit", "+", "-", "0", "*", "?"};

void tabmp_display(p, foo)
FILE *foo;
TableauMP *p;
{

 int i, j, ff, fff, n;
 EntierMP x, d;
 #if defined(PIPMP_LINEAR_VALUE_IS_MP)
 mpz_init(d);
 #endif

 fprintf(foo, "%ld/[%d * %d]\n", pipmp_cross_product, p->height, p->width);
 for(i = 0; i<p->height; i++){
   fff = ff = p->row[i].flags;
   /* if(fff ==0) continue; */
   #if defined(PIPMP_LINEAR_VALUE_IS_MP)
   mpz_set(d, Denom(p, i));
   #else
   d = Denom(p, i);
   #endif
   n = 0;
   while(fff){
     if(fff & 1) fprintf(foo, "%s ",Attr[n]);
     n++; fff >>= 1;
   }
   fprintf(foo, "%f #[", p->row[i].size);
   if(ff & Unit)
     for(j = 0; j<p->width; j++)
       fprintf(foo, " /%d/",(j == p->row[i].objet.unit)? 1: 0);
   else
     for(j = 0; j<p->width; j++){
       #if defined(PIPMP_LINEAR_VALUE_IS_MP)
       mpz_out_str(foo, 10, Index(p, i, j));
       putc(' ', foo);
       #else
       x = Index(p,i,j);
       fprintf(foo, FORMAT, x);
       fprintf(foo, " ");
       #endif
     }
   fprintf(foo, "]/");
   #if defined(PIPMP_LINEAR_VALUE_IS_MP)
   mpz_out_str(foo, 10, d);
   #else
   fprintf(foo, "%d", (int)d);
   #endif
   putc('\n', foo);
 }
 #if defined(PIPMP_LINEAR_VALUE_IS_MP)
 mpz_clear(d);
 #endif
}
