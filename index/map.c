#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <errno.h>
#include "kthread.h"
#include "kvec.h"
#include "kalloc.h"
#include "sdust.h"
#include "mmpriv.h"
#include "bseq.h"
#include "khash.h"

struct mm_tbuf_s {
	void *km;
	int rep_len, frag_gap;
};

mm_tbuf_t *mm_tbuf_init(void)
{
	mm_tbuf_t *b;
	b = (mm_tbuf_t*)calloc(1, sizeof(mm_tbuf_t));
	if (!(mm_dbg_flag & 1)) b->km = km_init();
	return b;
}

void mm_tbuf_destroy(mm_tbuf_t *b)
{
	if (b == 0) return;
	km_destroy(b->km);
	free(b);
}

void *mm_tbuf_get_km(mm_tbuf_t *b)
{
	return b->km;
}

static int mm_dust_minier(void *km, int n, mm128_t *a, int l_seq, const char *seq, int sdust_thres)
{
	int n_dreg, j, k, u = 0;
	const uint64_t *dreg;
	sdust_buf_t *sdb;
	if (sdust_thres <= 0) return n;
	sdb = sdust_buf_init(km);
	dreg = sdust_core((const uint8_t*)seq, l_seq, sdust_thres, 64, &n_dreg, sdb);
	for (j = k = 0; j < n; ++j) { // squeeze out minimizers that significantly overlap with LCRs
		int32_t qpos = (uint32_t)a[j].y>>1, span = a[j].x&0xff;
		int32_t s = qpos - (span - 1), e = s + span;
		while (u < n_dreg && (int32_t)dreg[u] <= s) ++u;
		if (u < n_dreg && (int32_t)(dreg[u]>>32) < e) {
			int v, l = 0;
			for (v = u; v < n_dreg && (int32_t)(dreg[v]>>32) < e; ++v) { // iterate over LCRs overlapping this minimizer
				int ss = s > (int32_t)(dreg[v]>>32)? s : dreg[v]>>32;
				int ee = e < (int32_t)dreg[v]? e : (uint32_t)dreg[v];
				l += ee - ss;
			}
			if (l <= span>>1) a[k++] = a[j]; // keep the minimizer if less than half of it falls in masked region
		} else a[k++] = a[j];
	}
	sdust_buf_destroy(sdb);
	return k; // the new size
}

static void collect_minimizers(void *km, const mm_mapopt_t *opt, const mm_idx_t *mi, int n_segs, const int *qlens, const char **seqs, mm128_v *mv)
{
	int i, n, sum = 0;
	mv->n = 0;
	for (i = n = 0; i < n_segs; ++i) {
		size_t j;
		mm_sketch(km, seqs[i], qlens[i], mi->w, mi->k, i, mi->flag&MM_I_HPC, mv);
		for (j = n; j < mv->n; ++j)
			mv->a[j].y += sum << 1;
		if (opt->sdust_thres > 0) // mask low-complexity minimizers
			mv->n = n + mm_dust_minier(km, mv->n - n, mv->a + n, qlens[i], seqs[i], opt->sdust_thres);
		sum += qlens[i], n = mv->n;
	}
}

#include "ksort.h"
#define heap_lt(a, b) ((a).x > (b).x)
KSORT_INIT(heap, mm128_t, heap_lt)

static inline int skip_seed(int flag, uint64_t r, const mm_seed_t *q, const char *qname, int qlen, const mm_idx_t *mi, int *is_self)
{
	*is_self = 0;
	if (qname && (flag & (MM_F_NO_DIAG|MM_F_NO_DUAL))) {
		const mm_idx_seq_t *s = &mi->seq[r>>32];
		int cmp;
		cmp = strcmp(qname, s->name);
		if ((flag&MM_F_NO_DIAG) && cmp == 0 && (int)s->len == qlen) {
			if ((uint32_t)r>>1 == (q->q_pos>>1)) return 1; // avoid the diagnonal anchors
			if ((r&1) == (q->q_pos&1)) *is_self = 1; // this flag is used to avoid spurious extension on self chain
		}
		if ((flag&MM_F_NO_DUAL) && cmp > 0) // all-vs-all mode: map once
			return 1;
	}
	if (flag & (MM_F_FOR_ONLY|MM_F_REV_ONLY)) {
		if ((r&1) == (q->q_pos&1)) { // forward strand
			if (flag & MM_F_REV_ONLY) return 1;
		} else {
			if (flag & MM_F_FOR_ONLY) return 1;
		}
	}
	return 0;
}

int compare_func(const s_mm128_t *a, const s_mm128_t *b) {   
	return (llabs(a->x) > llabs(b->x)) - (llabs(a->x) < llabs(b->x));
}

int compare_func_d_length(const s_mm128_t *b, const s_mm128_t *a) {   
	return ((a->y - a->x) > (b->y - b->x)) - ((a->y - a->x) < (b->y - b->x));
}

static s_mm128_t *collect_seed_hits(void *km, const mm_mapopt_t *opt, int max_occ, const mm_idx_t *mi, const char *qname, const mm128_v *mv, int qlen, int64_t *n_a, int *rep_len,
								  int *n_mini_pos, uint64_t **mini_pos)
{
	int i, n_m;
	int64_t ref_pos, read_pos, pos_bias = mi->k - 1;
	mm_seed_t *m;
	mm128_t p;
	s_mm128_t *s_a;
	m = mm_collect_matches(km, &n_m, qlen, max_occ, opt->max_max_occ, opt->occ_dist, mi, mv, n_a, rep_len, n_mini_pos, mini_pos);
	s_a = (s_mm128_t*)kmalloc(km, *n_a * sizeof(s_mm128_t));
	for (i = 0, *n_a = 0; i < n_m; ++i) {
		mm_seed_t *q = &m[i];
		const uint64_t *r = q->cr;
		uint32_t k;
		for (k = 0; k < q->n; ++k) {
			int32_t is_self, rpos = (uint32_t)r[k] >> 1;
			s_mm128_t *s_p;
			if (skip_seed(opt->flag, r[k], q, qname, qlen, mi, &is_self)) continue;
			//p = &a[(*n_a)];
			s_p = &s_a[(*n_a)++];
			if ((r[k]&1) == (q->q_pos&1)) { // forward strand
				p.x = (r[k]&0xffffffff00000000ULL) | (rpos);
				p.y = (uint64_t)q->q_span << 32 | (q->q_pos) >> 1;
				ref_pos = (int64_t)(int32_t)p.x + (int64_t)mi->seq[p.x<<1>>33].offset - pos_bias;
			} else{ 
				p.x = 1ULL<<63 | (r[k]&0xffffffff00000000ULL) | (rpos);
				p.y = (uint64_t)q->q_span << 32 | (q->q_pos) >> 1;
				ref_pos = -((int64_t)(int32_t)p.x + (int64_t)mi->seq[p.x<<1>>33].offset - pos_bias);
			}
			p.y |= (uint64_t)q->seg_id << MM_SEED_SEG_SHIFT;
			if (q->is_tandem) p.y |= MM_SEED_TANDEM;
			if (is_self) p.y |= MM_SEED_SELF;           
			read_pos = (int64_t)(int32_t)p.y - pos_bias;
			s_p->x = ref_pos;
			s_p->y = read_pos;
		}
	}
	kfree(km, m);//mi->seq[a[i].x<<1>>33].name, (int32_t)a[i].x, "+-"[a[i].x>>63], (int32_t)a[i].y
	qsort(s_a, *n_a, sizeof(s_mm128_t), compare_func);
	return s_a;
}


s_mm128_t *mm_map_frag_(const mm_idx_t *mi, int n_segs, const int *qlens, const char **seqs, int *n_regs, mm_reg1_t **regs, mm_tbuf_t *b, const mm_mapopt_t *opt, const char *qname, int *n_aa)
{
	int i, j, rep_len, qlen_sum, n_mini_pos;
	int64_t n_a;
	int64_t pre_rpos;
	int64_t cluster_count, n_return_aa = 0, return_cluster_count = *n_aa;
	int64_t c_bias = 5000;
	int64_t pcluster_count, min_cluster_size = 2;
	uint64_t *mini_pos;
	s_mm128_t *aa, *cluster_info, *s_p, *return_aa;
	mm128_v mv = {0,0,0};
	km_stat_t kmst;
	*n_aa = 1;
	aa = (s_mm128_t*)kmalloc(b->km, sizeof(s_mm128_t));
	for (i = 0, qlen_sum = 0; i < n_segs; ++i)
		qlen_sum += qlens[i], n_regs[i] = 0, regs[i] = 0;
	if (qlen_sum == 0 || n_segs <= 0 || n_segs > MM_MAX_SEG) return aa;
	if (opt->max_qlen > 0 && qlen_sum > opt->max_qlen) return aa;
	kfree(b->km, aa);
	collect_minimizers(b->km, opt, mi, n_segs, qlens, seqs, &mv);
	if (opt->q_occ_frac > 0.0f) mm_seed_mz_flt(b->km, &mv, opt->mid_occ, opt->q_occ_frac);
	aa = collect_seed_hits(b->km, opt, opt->mid_occ, mi, qname, &mv, qlen_sum, &n_a, &rep_len, &n_mini_pos, &mini_pos); 
	kfree(b->km, mv.a);
	kfree(b->km, mini_pos);
	if(return_cluster_count == -1){ 
		*n_aa = n_a;    
		if (b->km) {
			km_stat(b->km, &kmst);
			assert(kmst.n_blocks == kmst.n_cores); // otherwise, there is a memory leak
			if (kmst.largest > 1U<<28 || (opt->cap_kalloc > 0 && kmst.capacity > opt->cap_kalloc)) {
				km_destroy(b->km);
				b->km = km_init();
			}
		}  
		return aa;
	}
	pcluster_count = 0;
	pre_rpos = 0;
	for (i = 1; i < n_a; i++){
		if((llabs(aa[i].x) - llabs(aa[i - 1].x)) > c_bias) {
			if((i - pre_rpos) > min_cluster_size) pcluster_count++;  
			pre_rpos = i;
		}              
	} 
	if((i - pre_rpos) > min_cluster_size) pcluster_count++;
	if(pcluster_count == 0){ 
		kfree(b->km, aa);
		*n_aa = 1;        
		aa = (s_mm128_t*)kmalloc(b->km, sizeof(s_mm128_t));
		if (b->km) {
			km_stat(b->km, &kmst);
			assert(kmst.n_blocks == kmst.n_cores); // otherwise, there is a memory leak
			if (kmst.largest > 1U<<28 || (opt->cap_kalloc > 0 && kmst.capacity > opt->cap_kalloc)) {
				km_destroy(b->km);
				b->km = km_init();
			}
		}  
		return aa;
	}
	cluster_info = (s_mm128_t*)kmalloc(b->km, pcluster_count * sizeof(s_mm128_t));
	pcluster_count = 0;
	pre_rpos = 0;
	for (i = 1; i < n_a; i++){
		if((llabs(aa[i].x) - llabs(aa[i - 1].x)) > c_bias) {
			if((i - pre_rpos) > min_cluster_size){ 
				s_p = &cluster_info[pcluster_count++];
				s_p->x = pre_rpos;
				s_p->y = i;
			}            
			pre_rpos = i;            
		}              
	}
	if((i - pre_rpos) > min_cluster_size) {
		s_p = &cluster_info[pcluster_count++];
		s_p->x = pre_rpos;
		s_p->y = i;
	}

	qsort(cluster_info, pcluster_count, sizeof(s_mm128_t), compare_func_d_length);
	if(pcluster_count < return_cluster_count) return_cluster_count = pcluster_count;
	for (i = 0; i < return_cluster_count; i++) {    
		n_return_aa += cluster_info[i].y - cluster_info[i].x;
	} 
	*n_aa = 0;
	return_aa = (s_mm128_t*)kmalloc(b->km, n_return_aa * sizeof(s_mm128_t));
	for (i = 0; i < return_cluster_count; i++) {
		for (j = cluster_info[i].x; j < cluster_info[i].y; j++) {
			s_p = &return_aa[*n_aa];
			*n_aa += 1;
			s_p->x = aa[j].x;
			s_p->y = aa[j].y;   
		}
	}
	kfree(b->km, cluster_info);
	kfree(b->km, aa);
	if (b->km) {
		km_stat(b->km, &kmst);
		assert(kmst.n_blocks == kmst.n_cores); // otherwise, there is a memory leak
		if (kmst.largest > 1U<<28 || (opt->cap_kalloc > 0 && kmst.capacity > opt->cap_kalloc)) {
			km_destroy(b->km);
			b->km = km_init();
		}
	}   
	return return_aa;
}

















