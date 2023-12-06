#ifndef CMAPPY_H
#define CMAPPY_H

#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include "minimap.h"
#include "kseq.h"
#include "ksw2.h"
#include "kthread.h"
#include "kvec.h"
#include "kalloc.h"
#include "sdust.h"
#include "mmpriv.h"
#include "bseq.h"
#include "khash.h"

KSEQ_DECLARE(gzFile)

typedef struct {
	const char *ctg;
	int32_t ctg_start, ctg_end;
	int32_t qry_start, qry_end;
	int32_t blen, mlen, NM, ctg_len;
	uint8_t mapq, is_primary;
	int8_t strand, trans_strand;
	int32_t seg_id;
	int32_t n_cigar32;
	uint32_t *cigar32;
} mm_hitpy_t;



static inline void mm_reg2hitpy(const mm_idx_t *mi, mm_reg1_t *r, mm_hitpy_t *h)
{
	h->ctg = mi->seq[r->rid].name;
	h->ctg_len = mi->seq[r->rid].len;
	h->ctg_start = r->rs, h->ctg_end = r->re;
	h->qry_start = r->qs, h->qry_end = r->qe;
	h->strand = r->rev? -1 : 1;
	h->mapq = r->mapq;
	h->mlen = r->mlen;
	h->blen = r->blen;
	h->NM = r->blen - r->mlen + r->p->n_ambi;
	h->trans_strand = r->p->trans_strand == 1? 1 : r->p->trans_strand == 2? -1 : 0;
	h->is_primary = (r->id == r->parent);
	h->seg_id = r->seg_id;
	h->n_cigar32 = r->p->n_cigar;
	h->cigar32 = r->p->cigar;
}


static inline void mm_kfree(void *km, void *ptr)
{
	kfree(km, ptr);
}


static inline void mm_free_reg1(mm_reg1_t *r)
{
	free(r->p);
}

static inline kseq_t *mm_fastx_open(const char *fn)
{
	gzFile fp;
	fp = fn && strcmp(fn, "-") != 0? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	return kseq_init(fp);
}

static inline void mm_fastx_close(kseq_t *ks)
{
	gzFile fp;
	fp = ks->f->f;
	kseq_destroy(ks);
	gzclose(fp);
}

static inline int mm_verbose_level(int v)
{
	if (v >= 0) mm_verbose = v;
	return mm_verbose;
}

static inline void mm_reset_timer(void)
{
	extern double realtime(void);
	mm_realtime0 = realtime();
}

extern unsigned char seq_comp_table[256];
static inline s_mm128_t *mm_map_aux(const mm_idx_t *mi, const char *seq1, const char *seq2, int *n_regs, mm_tbuf_t *b, const mm_mapopt_t *opt, int *n_aa)
{
	mm_reg1_t *r;
	s_mm128_t *aa;

	Py_BEGIN_ALLOW_THREADS

	int qlen=strlen(seq1);   
	aa = mm_map_frag_(mi, 1, &qlen, &seq1, n_regs, &r, b, opt, NULL, n_aa);


	Py_END_ALLOW_THREADS
	return aa;
}




static inline int *mappy_k_cigar(const uint8_t *target, const uint8_t *query, int *n_cigarcode, int *zdropped, int *max_q, int *max_t, int match, int mismatch, int gap_open_1, int gap_extend_1, int gap_open_2, int gap_extend_2, int bw, int zdropvalue, int *delcount, int *inscount)
{
	int i, a = match, b = mismatch, *cigarcode; // a>0 and b<0
	extern unsigned char seq_nt4_table[256];
	int8_t mat[25] = { a,b,b,b,0, b,a,b,b,0, b,b,a,b,0, b,b,b,a,0, 0,0,0,0,0 };
	int tl = strlen(target), ql = strlen(query);
	uint8_t *ts, *qs;
	ksw_extz_t ez;

	memset(&ez, 0, sizeof(ksw_extz_t));

	ts = (uint8_t*)malloc(tl);
	qs = (uint8_t*)malloc(ql);
    
	for (i = 0; i < tl; ++i) ts[i] = seq_nt4_table[(uint8_t)target[i]]; // encode to 0/1/2/3
	for (i = 0; i < ql; ++i) qs[i] = seq_nt4_table[(uint8_t)query[i]];
	ksw_extd2_sse(0, ql, qs, tl, ts, 5, mat, gap_open_1, gap_extend_1, gap_open_2, gap_extend_2, bw, zdropvalue, 0, 0, &ez);
	if(ez.n_cigar > 0){
		cigarcode = (int*)malloc(2 * ez.n_cigar * sizeof(int));
		if(cigarcode != NULL){
			*n_cigarcode = 2 * ez.n_cigar;
			for (i = 0; i < ez.n_cigar; ++i){ // print CIGAR
				cigarcode[i * 2] = (int)ez.cigar[i]>>4;
				cigarcode[i * 2 + 1] = -((int)ez.cigar[i]&0xf) - 1;
				if( ( ( (int)ez.cigar[i]&0xf ) != 0 ) && ( ( (int)ez.cigar[i]>>4 ) > 9 ) )
					if( ( (int)ez.cigar[i]&0xf ) == 1 )
						*inscount += 1;
					else
						*delcount += 1;
			}
		}
		else{
			*n_cigarcode = 0;
			cigarcode = NULL;
		}
	}
	else{
		*n_cigarcode = 0;   
		cigarcode = NULL;
	}  
	*zdropped = ez.zdropped; 
	*max_q = ez.max_q;
	*max_t = ez.max_t;    

	free(ez.cigar); free(ts); free(qs);
	return cigarcode;

}

static char *mappy_fetch_seq(const mm_idx_t *mi, const char *name, int st, int en, int *len)
{
	int i, rid;
	char *s;
	*len = 0;
	rid = mm_idx_name2id(mi, name);
	if (rid < 0) return 0;
	if ((uint32_t)st >= mi->seq[rid].len || st >= en) return 0;
	if (en < 0 || (uint32_t)en > mi->seq[rid].len)
		en = mi->seq[rid].len;
	s = (char*)malloc(en - st + 1);
	*len = mm_idx_getseq(mi, rid, st, en, (uint8_t*)s);
	for (i = 0; i < *len; ++i)
		s[i] = "ACGTN"[(uint8_t)s[i]];
	s[*len] = 0;
	return s;
}

static mm_idx_t *mappy_idx_seq(int w, int k, int is_hpc, int bucket_bits, const char *seq, int len)
{
	const char *fake_name = "N/A";
	char *s;
	mm_idx_t *mi;
	s = (char*)calloc(len + 1, 1);
	memcpy(s, seq, len);
	mi = mm_idx_str(w, k, is_hpc, bucket_bits, 1, (const char**)&s, (const char**)&fake_name);
	free(s);
	return mi;
}

#endif
