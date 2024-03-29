/******************************************************************************
 *                     PIP : Parametric Integer Programming                   *
 ******************************************************************************
 *                                 traiter.c                                  *
 ******************************************************************************
 *                                                                            *
 * Copyright Paul Feautrier, 1988-2005                                        *
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
 * Written by Paul Feautrier                                                  *
 *                                                                            *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include <piplib-gmp/piplib-gmp.h>


#ifndef PIPMP_LINEAR_VALUE_IS_MP
# define PIPMP_LINEAR_VALUE_IS_MP
#endif

#define max(x,y) ((x) > (y)? (x) : (y))


extern long int pipmp_cross_product, pipmp_limit;
extern int pipmp_verbose;
extern FILE *pipmp_dump;
extern int pipmp_profondeur;
extern int pipmp_compa_count;

#if !defined(PIPMP_LINEAR_VALUE_IS_MP)
int llog(EntierMP x)
{int n = 0;
/* x must be positive, you dummy */
 if(x<0) x=-x;
 while(x) x >>= 1, n++;
 return(n);
}
#endif

static
int chercher(TableauMP *p, int masque, int n)
{int i;
 for(i = 0; i<n; i++)
     if(p->row[i].flags & masque) break;
 return(i);
}

/* il est convenu que pipmp_traiter ne doit modifier ni le tableau, ni le contexte;
   le tableau peut grandir en cas de coupure (+1 en hauteur et +1 en largeur
   si nparm != 0) et en cas de partage (+1 en hauteur)(seulement si nparm != 0).
   le contexte peut grandir en cas de coupure (+2 en hauteur et +1 en largeur)
   (seulement si nparm !=0) et en cas de partage (+1 en hauteur)(nparm !=0).
   On estime le nombre de coupures a llog(D) et le nombre de partages a
   ni.
*/

TableauMP *pipmp_expanser(TableauMP *tp, int virt, int reel, int ncol, 
			int off, int dh, int dw)
{
 int i, j, ff;
 char *q; EntierMP *pq;
 EntierMP *pp, *qq;
 TableauMP *rp;
 if(tp == NULL) return(NULL);
 rp = tabmp_alloc(reel+dh, ncol+dw, virt);

 #if defined(PIPMP_LINEAR_VALUE_IS_MP)
 mpz_set(rp->determinant, tp->determinant);
 #else
 rp->l_determinant = tp->l_determinant;
 for(i=0; i<tp->l_determinant; i++)
     rp->determinant[i] = tp->determinant[i];
 #endif
 pq = (EntierMP *) & (rp->row[virt+reel+dh]);
 for(i = off; i<virt + reel; i++)
     {ff = Flag(rp, i) = Flag(tp, i-off);
      #if defined(PIPMP_LINEAR_VALUE_IS_MP)
      mpz_set(Denom(rp, i), Denom(tp, i-off));
      #else
      Denom(rp, i) = Denom(tp, i-off);
      #endif
      if(ff & Unit) rp->row[i].objet.unit = tp->row[i-off].objet.unit;
      else {
	  rp->row[i].objet.val = pq;
	  pq +=(ncol + dw);
	  pp = tp->row[i-off].objet.val;
	  qq = rp->row[i].objet.val;
	  for(j = 0; j<ncol; j++)
             #if defined(PIPMP_LINEAR_VALUE_IS_MP)
	     mpz_set(*qq++, *pp++);
             #else
	     *qq++ = *pp++;
             #endif
	  }
      }
 return(rp);
}

static
int exam_coef(TableauMP *tp, int nvar, int ncol, int bigparm)
{int i, j ;
 int ff, fff;
 #if defined(PIPMP_LINEAR_VALUE_IS_MP)
 int x;
 #else
 EntierMP x;
 #endif
 EntierMP *p;
 
 for(i = 0; i<tp->height; i++)
     {ff = Flag(tp,i);
      if(ff == 0) break;
      if(ff == Unknown) {
           #if defined(PIPMP_LINEAR_VALUE_IS_MP)
           if(bigparm >= 0){
	      x = mpz_sgn(Index(tp,i, bigparm));
	   #else
	   if(bigparm >= 0 && (x = Index(tp,i, bigparm))) {
           #endif
	      if(x<0)    {
	                  Flag(tp, i) = Minus;
			  return(i);
		         }
	       else      
                         #if defined(PIPMP_LINEAR_VALUE_IS_MP)
	                 if(x>0)
                         #endif
	                 Flag(tp, i) = Plus;
	       continue;
	   }
	   ff = Zero;
	   p = &(tp->row[i].objet.val[nvar+1]);
	   for(j = nvar+1; j<ncol; j++) {
                #if defined(PIPMP_LINEAR_VALUE_IS_MP)
	        x = mpz_sgn(*p); p++ ;
	        #else
	        x = *p++;
                #endif
		if(x<0) fff = Minus;
		else if (x>0) fff = Plus;
		else fff = Zero;
		if(fff != Zero && fff != ff)
		    if(ff == Zero) ff = fff;
		    else {ff = Unknown;
			  break;
			 }
	       }
/* bug de'tecte' par [paf], 16/2/93 !
   Si tous les coefficients des parame`tres sont ne'gatifs
   et si le terme constant est nul, le signe est inconnu!!
   On traite donc spe'cialement le terme constant. */
           #if defined(PIPMP_LINEAR_VALUE_IS_MP)
	   x = mpz_sgn(Index(tp, i, nvar));
	   #else
	   x = Index(tp, i, nvar);
           #endif
	   if(x<0) fff = Minus;
	   else if(x>0) fff = Plus;
	   else fff = Zero;
/* ici on a le signe du terme constant */
	   switch(ff){
/* le signe est inconnu si les coefficients sont positifs et
   le terme constant ne'gatif */
	   case Plus: if(fff == Minus) ff = Unknown; break;
/* si les coefficients sont tous nuls, le signe est celui
   du terme constant */
	   case Zero: ff = fff; break;
/* le signe est inconnu si les coefficients sont ne'gatifs,
   sauf si le terme constant est egalement negatif. */
	   case Minus: if(fff != Minus) ff = Unknown; break;
/* enfin, il n'y a rien a` dire si le signe des coefficients est inconnu */
	   }
	   Flag(tp, i) = ff;
	   if(ff == Minus) return(i);
	  }
      }
 return(i);
}
 
static
void compa_test(TableauMP *tp, TableauMP *context,
		int ni, int nvar, int nparm, int nc)
{
 EntierMP discr[MAXPARM];
 int i, j;
 int ff;
 int cPlus, cMinus, isCritic;
 int verbold;
 TableauMP *tPlus, *tMinus;
 int p;
 struct pipmp_high_water_mark q;

 if(nparm == 0) return;
 if(nparm >= MAXPARM) {
     fprintf(stderr, "Too much parameters : %d\n", nparm);
     exit(1);
     }
 q = tabmp_hwm();
 #if defined(PIPMP_LINEAR_VALUE_IS_MP)
 for(i=0; i<=nparm; i++)
   mpz_init(discr[i]);
 #endif

 for(i = 0; i<ni + nvar; i++)
     {ff = Flag(tp,i);
      if(ff & (Critic | Unknown))
	  {isCritic = Pip_True;
	   for(j = 0; j<nvar; j++)
                 #if defined(PIPMP_LINEAR_VALUE_IS_MP)
		 if(mpz_sgn(Index(tp, i, j)) > 0)
                 #else
	         if(Index(tp, i, j) > 0)
                 #endif
		 {isCritic = Pip_False;
		  break;
		 }
           pipmp_compa_count++;
	   for(j = 0; j < nparm; j++)
           #if defined(PIPMP_LINEAR_VALUE_IS_MP)
	   mpz_set(discr[j], Index(tp, i, j+nvar+1));   /* loop body. */
	   mpz_set(discr[nparm], Index(tp, i, nvar));
           mpz_sub_ui(discr[nparm], discr[nparm], (isCritic ? 0 : 1));
           #else
	   discr[j] = Index(tp, i, j+nvar+1);           /* loop body. */
	   discr[nparm] = Index(tp, i, nvar)- (isCritic ? 0 : 1);
           #endif
           /* NdCed : Attention au contexte == NULL ! */
	   tPlus = pipmp_expanser(context, nparm, nc, nparm+1, nparm, 1, 0);
	   Flag(tPlus, nparm+nc) = Unknown;
	   for(j = 0; j<=nparm; j++)
           #if defined(PIPMP_LINEAR_VALUE_IS_MP)
	   mpz_set(Index(tPlus, nparm+nc, j),discr[j]); /* loop body. */
	   mpz_set(Denom(tPlus, nparm+nc), pipmp_UN);
           #else
	   Index(tPlus, nparm+nc, j) = discr[j];        /* loop body. */
	   Denom(tPlus, nparm+nc) = pipmp_UN;
           #endif
	   
	   p = solmp_hwm();
	   pipmp_traiter(tPlus, NULL, Pip_True, nparm, 0, nc+1, 0, -1);
	   cPlus = pipmp_is_not_Nil(p);
	   if(pipmp_verbose>0){
	     fprintf(pipmp_dump, "\nThe positive case has been found ");
	     fprintf(pipmp_dump, cPlus? "possible\n": "impossible\n");
	     fflush(pipmp_dump);
	   }

	   solmp_reset(p);
	   for(j = 0; j<nparm+1; j++)
           #if defined(PIPMP_LINEAR_VALUE_IS_MP)
	   mpz_neg(discr[j], discr[j]);                 /* loop body. */
	   mpz_sub_ui(discr[nparm], discr[nparm], (isCritic ? 1 : 2));
           #else
	   discr[j] = -discr[j];                        /* loop body. */
	   discr[nparm] = discr[nparm] - (isCritic ? 1 : 2);
           #endif
	   tMinus = pipmp_expanser(context, nparm, nc, nparm+1, nparm, 1, 0);
	   Flag(tMinus, nparm+nc) = Unknown;
	   for(j = 0; j<= nparm; j++)
           #if defined(PIPMP_LINEAR_VALUE_IS_MP)
	   mpz_set(Index(tMinus, nparm+nc, j),discr[j]);/* loop body. */
	   mpz_set(Denom(tMinus, nparm+nc), pipmp_UN);
           #else
	   Index(tMinus, nparm+nc, j) = discr[j];       /* loop body. */
	   Denom(tMinus, nparm+nc) = pipmp_UN;
           #endif
	   pipmp_traiter(tMinus, NULL, Pip_True, nparm, 0, nc+1, 0, -1);
	   cMinus = pipmp_is_not_Nil(p);
	   if(pipmp_verbose>0){
	     fprintf(pipmp_dump, "\nThe negative case has been found ");
	     fprintf(pipmp_dump, cMinus? "possible\n": "impossible\n");
	     fflush(pipmp_dump);
	   }

	   solmp_reset(p);
	   if (cPlus && cMinus) {
	       Flag(tp,i) = isCritic ? Critic : Unknown;
	     }
	   else if (cMinus)
	      {Flag(tp,i) = Minus;
	       break;
	      }
	   else {
	     Flag(tp,i) = cPlus ? Plus : Zero;
	   }
	  }
     }
 tabmp_reset(q);

 #if defined(PIPMP_LINEAR_VALUE_IS_MP)
 for(i=0; i<=nparm; i++)
   mpz_clear(discr[i]);
 #endif
 
 return;
}
 
 static
EntierMP *valeur(TableauMP *tp, int i, int j)
{
 if(Flag(tp, i) & Unit)
     return(tp->row[i].objet.unit == j ? &Denom(tp,i) : &pipmp_ZERO);
 else return(&Index(tp, i, j));
}
 
 static
void solution(TableauMP *tp, int nvar, int nparm)
{int i, j;
 int ncol = nvar + nparm + 1;

 solmp_list(nvar);
 for(i = 0; i<nvar; i++)
     {solmp_forme(nparm+1);
      for(j = nvar+1; j<ncol; j++)
	 solmp_val(*valeur(tp, i, j), Denom(tp,i));
      solmp_val(*valeur(tp, i, nvar), Denom(tp,i));
     }
}
 
 static
int choisir_piv(TableauMP *tp, int pivi, int nvar, int nligne)
{
 int j, k;
 EntierMP pivot, foo, x, y;
 int sgn_x, pivj = -1;

 #if defined(PIPMP_LINEAR_VALUE_IS_MP)
 mpz_init(pivot); mpz_init(foo); mpz_init(x); mpz_init(y);
 #endif
 
 for(j = 0; j<nvar; j++) {
    #if defined(PIPMP_LINEAR_VALUE_IS_MP)
    mpz_set(foo, Index(tp, pivi, j));
    if(mpz_sgn(foo) <= 0) continue;
    if(pivj < 0)
	{pivj = j;
         mpz_set(pivot, foo);
	 continue;
	}
    for(k = 0; k<nligne; k++)
        {mpz_mul(x, pivot, *valeur(tp, k, j)); 
         mpz_mul(y, *valeur(tp, k, pivj), foo);
         mpz_sub(x, x, y);
         pipmp_cross_product++;
         sgn_x = mpz_sgn(x);
         if(sgn_x) break;
	}
    if(sgn_x < 0)
        {pivj = j;
         mpz_set(pivot, foo);
        }
    #else
    if((foo = Index(tp, pivi, j)) <= 0) continue;
    if(pivj < 0)
	{pivj = j;
	 pivot = foo;
	 continue;
	}
    for(k = 0; k<nligne; k++)
	{x = pivot * (*valeur(tp, k, j)) - (*valeur(tp, k, pivj)) * foo;
	 pipmp_cross_product++;
	 if(x) break;
	}
    if(x < 0)
	{pivj = j;
	 pivot = foo;
	}
    #endif
 }
 
 #if defined(PIPMP_LINEAR_VALUE_IS_MP)
 mpz_clear(pivot); mpz_clear(foo); mpz_clear(x); mpz_clear(y);
 #endif

 return(pivj);
}

 
 static
int pivoter(TableauMP *tp, int pivi, int nvar, int nparm, int ni, int iq)

{int pivj;
 int ncol = nvar + nparm + 1;
 int nligne = nvar + ni;
 int i, j, k;
 EntierMP x, y, d, gcd, u, dpiv;
 int ff, fff;
 EntierMP pivot, foo, z;
 EntierMP ppivot, dppiv;
 EntierMP new[MAXCOL], *p, *q;
 EntierMP lpiv;
 int sgn_x;
 #if !defined(PIPMP_LINEAR_VALUE_IS_MP)
 char format_format[32];

 sprintf(format_format, "\nPivot %s/%s\n", FORMAT, FORMAT);
 #endif

 if(ncol >= MAXCOL) {
   fprintf(stdout, "Too much variables\n");
   exit(1);
 }
 if(0 > pivi || pivi >= nligne || Flag(tp, pivi) == Unit) {
   fprintf(stdout, "Syserr : pivoter : wrong pivot row\n");
   exit(1);
 }

 pivj = choisir_piv(tp, pivi, nvar, nligne);
 if(pivj < 0) return(-1);
 if(pivj >= nvar) {
   fprintf(stdout, "Syserr : pivoter : wrong pivot\n");
   exit(1);
 }

 #if defined(PIPMP_LINEAR_VALUE_IS_MP)
 mpz_init(x); mpz_init(y); mpz_init(d); 
 mpz_init(gcd); mpz_init(u); mpz_init(dpiv);
 mpz_init(lpiv); mpz_init(pivot); mpz_init(foo);
 mpz_init(z); mpz_init(ppivot); mpz_init(dppiv);

 for(i=0; i<ncol; i++)
   mpz_init(new[i]);

 mpz_set(pivot, Index(tp, pivi, pivj));
 mpz_set(dpiv, Denom(tp, pivi));
 mpz_gcd(d, pivot, dpiv);
 mpz_divexact(ppivot, pivot, d);
 mpz_divexact(dppiv, dpiv, d);
 #else
 pivot = Index(tp, pivi, pivj);
 dpiv = Denom(tp, pivi);
 d = pgcd(pivot, dpiv);
 ppivot = pivot/d;
 dppiv = dpiv/d;
 #endif
 
 if(pipmp_verbose>1){
   #if defined(PIPMP_LINEAR_VALUE_IS_MP)
   fprintf(pipmp_dump, "Pivot ");
   mpz_out_str(pipmp_dump, 10, ppivot);
   putc('/', pipmp_dump);
   mpz_out_str(pipmp_dump, 10, dppiv);
   putc('\n', pipmp_dump);
   #else
   fprintf(pipmp_dump, format_format, ppivot, dppiv);
   #endif
   fprintf(pipmp_dump, "%d x %d\n", pivi, pivj);
 }

 #if defined(PIPMP_LINEAR_VALUE_IS_MP)
 mpz_fdiv_qr(x, y, tp->determinant, dppiv); 
 #else
 for(i=0; i< tp->l_determinant; i++){
     d=pgcd(tp->determinant[i], dppiv);
     tp->determinant[i] /= d;
     dppiv /= d;
     }
 #endif

 #if defined(PIPMP_LINEAR_VALUE_IS_MP)
 if(mpz_sgn(y) != 0){ 
 #else
 if(dppiv != 1) {
 #endif
   fprintf(stderr, "Integer overflow\n");
   if(pipmp_verbose>0) fflush(pipmp_dump);
   //exit(1);
   return 42;
 }
 
 #if defined(PIPMP_LINEAR_VALUE_IS_MP)
 mpz_mul(tp->determinant, x, ppivot);
 #else
 for(i=0; i<tp->l_determinant; i++)
     if(llog(tp->determinant[i]) + llog(ppivot) < 8*sizeof(EntierMP)){
	 tp->determinant[i] *= ppivot;
	 break;
	 }
 if(i >= tp->l_determinant){
     tp->l_determinant++;
     if(tp->l_determinant >= MAX_DETERMINANT){
	 fprintf(stderr, "Integer overflow : %d\n", tp->l_determinant);
	 //exit(1);
	 return 42;
	 }
     tp->determinant[i] = ppivot;
     }
 #endif

 if(pipmp_verbose>1){
   fprintf(pipmp_dump, "determinant ");
   #if defined(PIPMP_LINEAR_VALUE_IS_MP)
   mpz_out_str(pipmp_dump, 10, tp->determinant);
   #else
   for(i=0; i<tp->l_determinant; i++)
	fprintf(pipmp_dump, FORMAT, tp->determinant[i]);
   #endif
   fprintf(pipmp_dump, "\n");
 }

 
 for(j = 0; j<ncol; j++)
   #if defined(PIPMP_LINEAR_VALUE_IS_MP)
   if(j==pivj)
     mpz_set(new[j], dpiv);
   else 
     mpz_neg(new[j], Index(tp, pivi, j));
   #else
   new[j] = (j == pivj ? dpiv : -Index(tp, pivi, j));
   #endif

 for(k = 0; k<nligne; k++){
   if(Flag(tp,k) & Unit)continue;
   if(k == pivi)continue;
   #if defined(PIPMP_LINEAR_VALUE_IS_MP)
   mpz_set(foo, Index(tp, k, pivj));
   mpz_gcd(d, pivot, foo);
   mpz_divexact(lpiv, pivot, d);
   mpz_divexact(foo, foo, d);
   mpz_set(d, Denom(tp,k));
   mpz_mul(gcd, lpiv, d);
   mpz_set(Denom(tp, k), gcd);
   #else
   foo = Index(tp, k, pivj);
   d = pgcd(pivot, foo);
   lpiv = pivot/d;
   foo /= d;
   d = Denom(tp,k);
   gcd = lpiv * d;
   Denom(tp, k) = gcd;
   #endif
   p = tp->row[k].objet.val;
   q = tp->row[pivi].objet.val;
   for(j = 0; j<ncol; j++){
     if(j == pivj)
     #if defined(PIPMP_LINEAR_VALUE_IS_MP)
       mpz_mul(z, dpiv, foo);
     #else
       z = dpiv * foo;
     #endif
     else {
     #if defined(PIPMP_LINEAR_VALUE_IS_MP)
       mpz_mul(z, *p, lpiv);
       mpz_mul(y, *q, foo);
       mpz_sub(z, z, y);
     #else
       z = (*p) * lpiv - (*q) * foo;
     #endif
     }
     q++;
     pipmp_cross_product++;
     #if defined(PIPMP_LINEAR_VALUE_IS_MP)
     mpz_set(*p, z);
     p++;
     if(mpz_cmp_ui(gcd, 1) != 0)
       mpz_gcd(gcd, gcd, z);
     #else
     *p++ = z;
     if(gcd != 1)
       gcd = pgcd(gcd, z);
     #endif
   }
   #if defined(PIPMP_LINEAR_VALUE_IS_MP)
   if(mpz_cmp_ui(gcd, 1) != 0){
     p = tp->row[k].objet.val;
     for(j = 0; j<ncol; j++){
       mpz_divexact(*p, *p, gcd);
       p++;
     }
   }
   mpz_divexact(Denom(tp,k), Denom(tp,k), gcd);
   #else
   if(gcd != 1) {
    p = tp->row[k].objet.val;
    for(j = 0; j<ncol; j++)
      *p++ /= gcd;
      Denom(tp,k) = Denom(tp,k)/gcd;
   }
   #endif
 }
 p = tp->row[pivi].objet.val;
 for(k = 0; k<nligne; k++)
   if((Flag(tp, k) & Unit) && tp->row[k].objet.unit == pivj) break;
 Flag(tp, k) = Plus;
 tp->row[k].objet.val = p;
 for(j = 0; j<ncol; j++)
   #if defined(PIPMP_LINEAR_VALUE_IS_MP)
   mpz_set(*p++, new[j]);
   #else
   *p++ = new[j];
   #endif

 #if defined(PIPMP_LINEAR_VALUE_IS_MP)
 mpz_set(Denom(tp, k), pivot);
 Flag(tp, pivi) = Unit | Zero;
 mpz_set(Denom(tp, pivi), pipmp_UN);
 #else
 Denom(tp, k) = pivot; 
 Flag(tp, pivi) = Unit | Zero;
 Denom(tp, pivi) = pipmp_UN;
 #endif
 tp->row[pivi].objet.unit = pivj;

 for(k = 0; k<nligne; k++){
   ff = Flag(tp, k);
   if(ff & Unit) continue;
   #if defined(PIPMP_LINEAR_VALUE_IS_MP)
   sgn_x = mpz_sgn(Index(tp, k, pivj));
   #else
   sgn_x = Index(tp, k, pivj);
   #endif
   if(sgn_x < 0) fff = Minus;
   else if(sgn_x == 0) fff = Zero;
   else fff = Plus;
   if(fff != Zero && fff != ff)
     if(ff == Zero) ff = (fff == Minus ? Unknown : fff);
     else ff = Unknown;
   Flag(tp, k) = ff;
 }

 if(pipmp_verbose>2){
   fprintf(pipmp_dump, "just pivoted\n");
   tabmp_display(tp, pipmp_dump);
 }

 #if defined(PIPMP_LINEAR_VALUE_IS_MP)
 mpz_clear(x); mpz_clear(y); mpz_clear(d); mpz_clear(gcd);
 mpz_clear(u); mpz_clear(dpiv); mpz_clear(lpiv);
 mpz_clear(pivot); mpz_clear(foo); mpz_clear(z);
 mpz_clear(ppivot); mpz_clear(dppiv);

 for(i=0; i<ncol; i++)
   mpz_clear(new[i]);
 #endif

 return(0);
}

/* dans cette version, "pipmp_traiter" modifie ineq; par contre
   le contexte est immediatement recopie' */

int pipmp_traiter(tp, ctxt, iq, nvar, nparm, ni, nc, bigparm)
TableauMP *tp, *ctxt;
int iq, nvar, nparm, ni, nc, bigparm;
{
  int val_ret;
 int j;
 int pivi, nligne, ncol;
 struct pipmp_high_water_mark x;
 TableauMP *context;
 int dch, dcw;
 double s, t, d, smax;
 int i;
 struct pipmp_L temp;
 EntierMP discr[MAXPARM];

 #if !defined(PIPMP_LINEAR_VALUE_IS_MP)
 EntierMP D = pipmp_UN;
 #endif

 #if defined(PIPMP_LINEAR_VALUE_IS_MP)
 for(i=0; i<MAXPARM; i++)
   mpz_init(discr[i]);
 dcw = mpz_sizeinbase(tp->determinant, 2);
 #else
 dcw = 0;
 for(i=0; i<tp->l_determinant; i++)
   dcw += llog(tp->determinant[i]);
 #endif
 dch = 2 * dcw + 1;
 x = tabmp_hwm();
 nligne = nvar+ni;


 context = pipmp_expanser(ctxt, 0, nc, nparm+1, 0, dch, dcw);
 
 /*
 sort the rows in increasing order of the largest coefficient
*/
   
 smax = 0.;

 for(i=nvar; i<nligne; i++){
   if(Flag(tp,i) & Unit) continue;
   s = 0.;
   #if defined(PIPMP_LINEAR_VALUE_IS_MP)
   d = mpz_get_d(Denom(tp,i));     
   for(j=0; j<nvar; j++){
     t = mpz_get_d(Index(tp,i,j))/d;
     s = max(s, abs(t));
   }
   #else
   d = (float) Denom(tp,i);     
   for(j=0; j<nvar; j++){
       t = Index(tp,i,j)/d;
       s = max(s, abs(t));
       }
   #endif
   tp->row[i].size = s;
   smax = max(s, smax);
 }
     
 for(i=nvar; i<nligne; i++){
   if(Flag(tp,i) & Unit) continue;
   s = smax;
   pivi = i;
   for(j=i; j<nligne; j++){
     if(Flag(tp,j) & Unit) continue;
     if(tp->row[j].size < s){
       s = tp->row[i].size;
       pivi = j;
     }
   }
   if(pivi != i) {
     temp = tp->row[pivi];
     tp->row[pivi] = tp->row[i];
     tp->row[i]=temp;
   }
 }   


 for(;;) {
   if(pipmp_verbose>2){
     fprintf(pipmp_dump, "debut for\n");
     tabmp_display(tp, pipmp_dump);
     fflush(pipmp_dump);
   }
   nligne = nvar+ni; ncol = nvar+nparm+1;
   if(nligne > tp->height || ncol > tp->width) {
     fprintf(stdout, "Syserr : pipmp_traiter : tableau too small\n");
     exit(1);
   }
   pivi = chercher(tp, Minus, nligne);
   if(pivi < nligne) goto pirouette;	       /* There is a negative row   */
   
   pivi = exam_coef(tp, nvar, ncol, bigparm);

   if(pipmp_verbose>2){
     fprintf(pipmp_dump, "coefs examined\n");
     tabmp_display(tp, pipmp_dump);
     fflush(pipmp_dump);
   }

   if(pivi < nligne) goto pirouette;
   /* There is a row whose coefficients are negative */
   compa_test(tp, context, ni, nvar, nparm, nc);
   if(pipmp_verbose>2){
     fprintf(pipmp_dump, "compatibility tested\n");
     tabmp_display(tp, pipmp_dump);
     fflush(pipmp_dump);
   }

   
   pivi = chercher(tp, Minus, nligne);
   if(pivi < nligne) goto pirouette;
   /* The compatibility test has found a negative row */
   pivi = chercher(tp, Critic, nligne);
   if(pivi >= nligne)pivi = chercher(tp, Unknown, nligne);
   /* Here, the problem tree splits        */
   if(pivi < nligne) {
     TableauMP * ntp;
     EntierMP com_dem;
     struct pipmp_high_water_mark q;
     if(nc >= context->height) {
       #if defined(PIPMP_LINEAR_VALUE_IS_MP)
       dcw = mpz_sizeinbase(context->determinant,2);
       #else
       dcw = 0;
       for(i=0; i<tp->l_determinant; i++)
       dcw += llog(tp->determinant[i]);
       #endif
       dch = 2 * dcw + 1;
       context = pipmp_expanser(context, 0, nc, nparm+1, 0, dch, dcw);
     }
     if(nparm >= MAXPARM) {
       fprintf(stdout, "Too much parameters : %d\n", nparm);
       exit(2);
     }
     q = tabmp_hwm();
     if(pipmp_verbose>1)
       fprintf(stdout,"pipmp_profondeur %d %lx\n", pipmp_profondeur, q.top);
     ntp = pipmp_expanser(tp, nvar, ni, ncol, 0, 0, 0);
     fflush(stdout);
     solmp_if();
     solmp_forme(nparm+1);
     #if defined(PIPMP_LINEAR_VALUE_IS_MP)
     mpz_init_set_ui(com_dem, 0);
     for(j = 0; j<nparm; j++) {
       mpz_set(discr[j], Index(tp, pivi, j + nvar +1));
       mpz_gcd(com_dem, com_dem, discr[j]);
     }
     mpz_set(discr[nparm], Index(tp, pivi, nvar));
     mpz_gcd(com_dem, com_dem, discr[nparm]);
     for(j = 0; j<=nparm; j++) {
       mpz_divexact(discr[j], discr[j], com_dem);
       mpz_set(Index(context, nc, j), discr[j]);
       solmp_val(discr[j], pipmp_UN);
     }
     mpz_clear(com_dem);
     Flag(context, nc) = Unknown;
     mpz_set(Denom(context, nc), pipmp_UN);
     #else
     com_dem = 0;
     for(j = 0; j<nparm; j++) {
       discr[j] = Index(tp, pivi, j + nvar +1);
       com_dem = pgcd(com_dem, discr[j]);
     }
     discr[nparm] = Index(tp, pivi, nvar);
     com_dem = pgcd(com_dem, discr[nparm]);
     for(j = 0; j<=nparm; j++) {
       discr[j] /= com_dem;
       Index(context, nc, j) = discr[j];
       solmp_val(discr[j], pipmp_UN);
     }
     Flag(context, nc) = Unknown;
     Denom(context, nc) = pipmp_UN;
     #endif
     Flag(ntp, pivi) = Plus;
     pipmp_profondeur++;
     fflush(stdout);
     if(pipmp_verbose > 0) fflush(pipmp_dump);
     #if defined(PIPMP_LINEAR_VALUE_IS_MP)
     pipmp_traiter(ntp, context, iq, nvar, nparm, ni, nc+1, bigparm);
     pipmp_profondeur--;
     tabmp_reset(q);
     if(pipmp_verbose>1)
       fprintf(stdout, "descente %d %lx\n", pipmp_profondeur, tabmp_hwm().top);
     for(j = 0; j<nparm; j++)
       mpz_neg(Index(context, nc, j), Index(context, nc, j));
     mpz_add_ui(Index(context, nc, nparm), Index(context, nc, nparm), 1);
     mpz_neg(Index(context, nc, nparm), Index(context, nc, nparm));
     Flag(tp, pivi) = Minus;
     mpz_set(Denom(context, nc), pipmp_UN);
     #else
     pipmp_traiter(ntp, context, iq, nvar, nparm, ni, nc+1, bigparm);
     pipmp_profondeur--;
     tabmp_reset(q);
     if(pipmp_verbose>1)
       fprintf(stderr, "descente %d %lx\n", pipmp_profondeur, tabmp_hwm().top);
     for(j = 0; j<nparm; j++)
       Index(context, nc, j) = - Index(context, nc, j);
     Index(context, nc, nparm) = - Index(context, nc, nparm) -1;
     Flag(tp, pivi) = Minus;
     Denom(context, nc) = pipmp_UN;
     #endif
     nc++;
     goto pirouette;
   }

   /* Here, all rows are positive. Do we need an integral solution?      */
   if(!iq) {
     solution(tp, nvar, nparm);
     break;
   }
/* Yes we do! */

   pivi = pipmp_integrer(&tp, &context, &nvar, &nparm, &ni, &nc);

   if(pivi > 0) goto pirouette;
		    /* A cut has been inserted and is always negative */
/* Here, either there is an integral solution, */
   if(pivi == 0) solution(tp, nvar, nparm);
/* or no solution exists */
   else solmp_nil();
   break;

/* Here, a negative row has been found. The call to <<pivoter>> executes
      a pivoting step                                                 */

pirouette :
  val_ret = pivoter(tp, pivi, nvar, nparm, ni, iq);
     if(val_ret < 0) {
       solmp_nil();
       break;
     }
     else if (val_ret == 42)
       return 42;
 }

 /* Danger : a premature return would induce memory leaks   */
 tabmp_reset(x);
 #if defined(PIPMP_LINEAR_VALUE_IS_MP)
 for(i=0; i<MAXPARM; i++)
   mpz_clear(discr[i]);
 #endif

 return 0;
}
