#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bwtgap.h"
#include "bwtaln.h"

#define STATE_M 0
#define STATE_I 1
#define STATE_D 2

#define aln_score(m,o,e,p) ((m)*(p)->s_mm + (o)*(p)->s_gapo + (e)*(p)->s_gape)

gap_stack_t *gap_init_stack(int max_mm, int max_gapo, int max_gape, const gap_opt_t *opt)
{
	int i;
	gap_stack_t *stack;
	stack = (gap_stack_t*)calloc(1, sizeof(gap_stack_t));
	stack->n_stacks = aln_score(max_mm+1, max_gapo+1, max_gape+1, opt);
	stack->stacks = (gap_stack1_t*)calloc(stack->n_stacks, sizeof(gap_stack1_t));
	for (i = 0; i != stack->n_stacks; ++i) {
		gap_stack1_t *p = stack->stacks + i;
		p->m_entries = 4;
		p->stack = (gap_entry_t*)calloc(p->m_entries, sizeof(gap_entry_t));
	}
	return stack;
}

void gap_destroy_stack(gap_stack_t *stack)
{
	int i;
	for (i = 0; i != stack->n_stacks; ++i)
        free(stack->stacks[i].stack);
	free(stack->stacks);
	free(stack);
}

static void gap_reset_stack(gap_stack_t *stack)
{
	int i;
	for (i = 0; i != stack->n_stacks; ++i)
		stack->stacks[i].n_entries = 0;
	stack->best = stack->n_stacks;
	stack->n_entries = 0;
}

static inline void gap_push(gap_stack_t *stack, int i, bwtint_t k, bwtint_t l,
        bwtint_t rev_k, bwtint_t rev_l, int n_mm, int n_gapo, int n_gape,
		int state, int is_diff, const gap_opt_t *opt)
{
	int score;
	gap_entry_t *p;
	gap_stack1_t *q;
	score = aln_score(n_mm, n_gapo, n_gape, opt);
	q = stack->stacks + score;
	if (q->n_entries == q->m_entries) {
		q->m_entries <<= 1; // q->m_entries扩大一倍
		q->stack = (gap_entry_t*)realloc(q->stack, sizeof(gap_entry_t) * q->m_entries);
	}
	p = q->stack + q->n_entries;
	p->info = (u_int32_t)score<<21 | i;
    p->k = k;
    p->l = l;
    p->rev_k = rev_k;
    p->rev_l = rev_l;
	p->n_mm = n_mm; 
    p->n_gapo = n_gapo;
    p->n_gape = n_gape;
    p->state = state;
	p->last_diff_pos = is_diff? i : 0; // can it saving more info ?
	++(q->n_entries);
	++(stack->n_entries);
	if (stack->best > score)
        //best score is the lowest score
        stack->best = score;
}

static inline void gap_pop(gap_stack_t *stack, gap_entry_t *e)
{
	gap_stack1_t *q;
	q = stack->stacks + stack->best;
	*e = q->stack[q->n_entries - 1];
	--(q->n_entries);
	--(stack->n_entries);
	if (q->n_entries == 0 && stack->n_entries) { // reset best
		int i;
		for (i = stack->best + 1; i < stack->n_stacks; ++i)
			if (stack->stacks[i].n_entries != 0)
                break;
		stack->best = i;
	} else if (stack->n_entries == 0)
        stack->best = stack->n_stacks;
}

static inline void gap_shadow(int x, int len, bwtint_t max, int last_diff_pos, bwt_width_t *w)
{
	int i, j;
	for (i = j = 0; i < last_diff_pos; ++i) {
		if (w[i].w > x)
            w[i].w -= x;
		else if (w[i].w == x) {
			w[i].bid = 1;
			w[i].w = max - (++j);
		} // else should not happen
	}
}

static inline int int_log2(uint32_t v)
{
	int c = 0;
	if (v & 0xffff0000u) { v >>= 16; c |= 16; }
	if (v & 0xff00) { v >>= 8; c |= 8; }
	if (v & 0xf0) { v >>= 4; c |= 4; }
	if (v & 0xc) { v >>= 2; c |= 2; }
	if (v & 0x2) c |= 1;
	return c;
}

bwt_aln1_t *bwt_match_gap(bwt_aux_t *aux, int *_n_aln)
{
    // init additional struct from input 
    gap_opt_t *opt = aux->opt;
    gap_stack_t *stack = aux->stack;
    ubyte_t *seq = aux->strand == 1? aux->rc_seq: aux->seq;
    Idx2BWT *bi_bwt = aux->bi_bwt;
    bwt_width_t *width = aux->width_back;
    bwt_width_t *width_seed = aux->width_seed;

	int best_score = aln_score(opt->max_diff+1, opt->max_gapo+1, opt->max_gape+1, opt);
    //max_diff is cal by error rate and length
	int best_diff = opt->max_diff + 1, max_diff = opt->max_diff;
	int best_cnt = 0;
	int max_entries = 0, j, n_aln, m_aln;
	bwt_aln1_t *aln, *aln_new;
    int len = aux->len;
	gap_entry_t e;

	m_aln = 10; n_aln = 0;
	aln = (bwt_aln1_t *)calloc(m_aln, sizeof(bwt_aln1_t));

	//for (j = 0; j != len; ++j) printf("#0 %d: [%d,%u]\t[%d,%u]\n", j, w[0][j].bid, w[0][j].w, w[1][j].bid, w[1][j].w);
	gap_reset_stack(stack); // reset stack
	gap_push(stack, len, 0, bi_bwt->bwt->textLength, 0, bi_bwt->bwt->textLength, 0, 0, 0, 0, 0, opt);

	while (stack->n_entries) {
		int i, m, m_seed = 0, hit_found, allow_diff, allow_M, tmp;
		bwtint_t k, l, rev_k, rev_l, sa_k[4], sa_l[4], rev_sa_k[4], rev_sa_l[4], occ;

		if (max_entries < stack->n_entries)
            max_entries = stack->n_entries;
		if (stack->n_entries > opt->max_entries)
            break;
		gap_pop(stack, &e); // get the best entry
		k = e.k; 
        l = e.l; // SA interval
        rev_k = e.rev_k;
        rev_l = e.rev_l;
		i = e.info&0xffff; // leaving length not mapped
		if (!(opt->mode & BWA_MODE_NONSTOP) && e.info>>21 > best_score + opt->s_mm)
            break; // no need to proceed

		m = max_diff - (e.n_mm + e.n_gapo);
        // m: available diff number
		if (opt->mode & BWA_MODE_GAPE)
            m -= e.n_gape;
		if (m < 0)
            continue;
		if (width_seed) { // apply seeding
			m_seed = opt->max_seed_diff - (e.n_mm + e.n_gapo);
			if (opt->mode & BWA_MODE_GAPE)
                m_seed -= e.n_gape;
		}
		if (i > 0 && m < width[i-1].bid)
            continue;

		// check whether a hit is found
		hit_found = 0;
		if (i == 0)
            //all chars have been mapped
            hit_found = 1;
		else if (m == 0 && (e.state == STATE_M || (opt->mode&BWA_MODE_GAPE) || e.n_gape == opt->max_gape)) {
            // no diff allowed
			if (bwt_match_exact(bi_bwt, seq, i, &k, &l, &rev_k, &rev_l))
                hit_found = 1;
			else
                continue; // no hit, skip
		}

		if (hit_found) { // action for found hits
			int score = aln_score(e.n_mm, e.n_gapo, e.n_gape, opt);
			int do_add = 1;
			//printf("#2 hits found: %d:(%u,%u)\n", e.n_mm+e.n_gapo, k, l);
			if (n_aln == 0) {
                //如果是第一个hit。。。
				best_score = score;
				best_diff = e.n_mm + e.n_gapo;
				if (opt->mode & BWA_MODE_GAPE)
                    best_diff += e.n_gape;
				if (!(opt->mode & BWA_MODE_NONSTOP))
					max_diff = (best_diff + 1 > opt->max_diff)? opt->max_diff : best_diff + 1; // top2 behaviour
			}
			if (score == best_score)
                best_cnt += l - k + 1;
			else if (best_cnt > opt->max_top2)
                break; // top2b behaviour
			if (e.n_gapo) { 
                // check whether the hit has been found. this may happen when a gap occurs in a tandem repeat
                // check whether the last mapping pos is gap extension
				for (j = 0; j != n_aln; ++j)
					if (aln[j].k == k && aln[j].l == l)
                        break;
				if (j < n_aln)
                    do_add = 0;
			}
			if (do_add) { // append
                //能mapping上所有，并且不与之前的hit有相同的k，l
				bwt_aln1_t *p;
				gap_shadow(l - k + 1, len, bi_bwt->bwt->textLength, e.last_diff_pos, width);
				if (n_aln == m_aln) {
					m_aln <<= 1;
					aln_new = (bwt_aln1_t *)realloc(aln, m_aln * sizeof(bwt_aln1_t)); /* mem leak */
                    if(aln_new != NULL){
                        aln = aln_new;
                    }
                    else 
                        break;
					memset(aln + m_aln/2, 0, m_aln/2*sizeof(bwt_aln1_t));
				}
				p = aln + n_aln;
                //aln保存所有满足diff限制的匹配情况
				p->n_mm = e.n_mm;
                p->n_gapo = e.n_gapo;
                p->n_gape = e.n_gape;
				p->k = k;
                p->l = l;
                p->strand = aux->strand;
                p->rev_k = rev_k;
                p->rev_l = rev_l;
				p->score = score;
				++n_aln;
                //n_aln保存read在满足diff限制下，能够匹配上ref的个数
			}
			continue;
		}

		--i;    //i：剩余reads长度
        //bwt_2occ4(bwt, k - 1, l, cnt_k, cnt_l); 
        BWTAllSARangesBackward_Bidirection(bi_bwt, k, l, rev_k, rev_l, sa_k, sa_l,rev_sa_k, rev_sa_l);
        // retrieve Occ values
		occ = l - k + 1;
		// test whether diff is allowed
		allow_diff = allow_M = 1;
		if (i > 0) {
			int ii = i - (len - opt->seed_len);
            //ii: 剩余的seed长度
			if (width[i-1].bid > m-1)
                allow_diff = 0;
			else if (width[i-1].bid == m-1 && width[i].bid == m-1 && width[i-1].w == width[i].w)
                allow_M = 0;
			if (width_seed && ii > 0) {
				if (width_seed[ii-1].bid > m_seed-1) allow_diff = 0;
				else if (width_seed[ii-1].bid == m_seed-1 && width_seed[ii].bid == m_seed-1
						 && width_seed[ii-1].w == width_seed[ii].w)
                    allow_M = 0;
			}
		}
		// indels
        tmp = (opt->mode & BWA_MODE_LOGGAP)? int_log2(e.n_gape + e.n_gapo)/2+1 : e.n_gapo + e.n_gape;
		if (allow_diff && i >= opt->indel_end_skip + tmp && len - i >= opt->indel_end_skip + tmp) {
			if (e.state == STATE_M) { // gap open
				if (e.n_gapo < opt->max_gapo) { // gap open is allowed
					// insertion
					gap_push(stack, i, k, l, rev_k, rev_l, e.n_mm, e.n_gapo + 1, e.n_gape, STATE_I, 1, opt);
					// deletion
					for (j = 0; j != ALPHABET_SIZE; ++j) {
						k = sa_k[j];
						l = sa_l[j];
                        rev_k = rev_sa_k[j];
                        rev_l = rev_sa_l[j];
						if (k <= l)
                            gap_push(stack, i + 1, k, l,rev_k,rev_l, e.n_mm, e.n_gapo + 1, e.n_gape, STATE_D, 1, opt);
					}
				}
			} else if (e.state == STATE_I) { // extention of an insertion
				if (e.n_gape < opt->max_gape) // gap extention is allowed
					gap_push(stack, i, k, l, rev_k, rev_l, e.n_mm, e.n_gapo, e.n_gape + 1, STATE_I, 1, opt);
			} else if (e.state == STATE_D) { // extention of a deletion
				if (e.n_gape < opt->max_gape) { // gap extention is allowed
					if (e.n_gape + e.n_gapo < max_diff || occ < opt->max_del_occ) {
						for (j = 0; j != 4; ++j) {
						    k = sa_k[j];
						    l = sa_l[j];
                            rev_k = rev_sa_k[j];
                            rev_l = rev_sa_l[j];
							if (k <= l)
                                gap_push(stack, i + 1, k, l, rev_k, rev_l, e.n_mm, e.n_gapo, e.n_gape + 1, STATE_D, 1, opt);
						}
					}
				}
			}
		}
		// mismatches
		if (allow_diff && allow_M) { // mismatch is allowed
			for (j = 1; j <= 4; ++j) {
                int c = (seq[i] + j) & 3;
				int is_mm = (j != 4 || seq[i] > 3);
                // when j = 4, c = seq[i]
                //TODO c is ??, and the result of BWTAllOccValue ??
				k = sa_k[c];
				l = sa_l[c];
                rev_k = rev_sa_k[c];
                rev_l = rev_sa_l[c];
				if (k <= l)
                    gap_push(stack, i, k, l, rev_k, rev_l, e.n_mm + is_mm, e.n_gapo, e.n_gape, STATE_M, is_mm, opt);
			}
		} else if (seq[i] < 4) {
            // seq[i] is not N
            // try exact match only
            int c = seq[i] & 3;
			k = sa_k[c];
			l = sa_l[c];
            rev_k = rev_sa_k[c];
            rev_l = rev_sa_l[c];
			if (k <= l)
                gap_push(stack, i, k, l, rev_k, rev_l, e.n_mm, e.n_gapo, e.n_gape, STATE_M, 0, opt);
		}
	}

	*_n_aln = n_aln;
    //aln保存着n_aln个aln结果
	return aln;
}

/*
 * function is used for backtracing base on opt limition
 * is_backward: mark the backtracing is forward or backward
 * len: the length we need to extend
 * width: is from the integerated one 
 * seq: is also the integerated one
 * max_pos: 0 - base, foreward or backward extended pos, relative to
 *          origin seq
 * !!! only return longest result
 * result: -1 - if it can not extend to _max_pos 
 *          1 - leaving length equal 0
 *          2 - leaving length can be not 0
 */
static int bwt_backtracing_search(bwt_aux_t *aux, bwt_aln1_t *aln, int is_backward,
        int *_max_pos)
{
    gap_opt_t *opt = aux->opt;
	int best_score = aln_score(opt->max_diff+1, opt->max_gapo+1, opt->max_gape+1, opt);
	int best_diff = opt->max_diff + 1, max_diff = opt->max_diff;
	int max_entries = 0, j, max_pos;
	gap_entry_t e;
	int i, m, tmp, start, end;
	bwtint_t k, l, rev_k, rev_l, sa_k[4], sa_l[4], rev_sa_k[4], rev_sa_l[4];
    bwtint_t occ;
    int allow_diff;
    int len = aux->len, real_pos;
    gap_stack_t *stack = aux->stack;

    // init
    ubyte_t *seq = aux->strand == 0? aux->seq: aux->rc_seq;
    bwt_width_t *width = (is_backward == 1? aux->width_back : aux->width_fore);

    start = aln->start, end = aln->end;
    max_pos = *_max_pos;
    
    max_entries = 0;
    //real_max = (is_backward == 0? (max_pos-start+1):(max+len-aln->end));
    // DFS
	while (stack->n_entries != 0) {
		if (max_entries < stack->n_entries)
            max_entries = stack->n_entries;
		if (stack->n_entries > opt->max_entries)
            break;
		gap_pop(stack, &e); // get the best entry
		k = e.k; 
        l = e.l; // SA interval
        rev_k = e.rev_k;
        rev_l = e.rev_l; //rev SA interval
		i = e.info&0xffff; // leaving length have not been extended
		if (!(opt->mode & BWA_MODE_NONSTOP) && e.info>>21 > best_score + opt->s_mm)
            break; // no need to proceed

		m = max_diff - (e.n_mm + e.n_gapo);
        //m: available diff number
		if (opt->mode & BWA_MODE_GAPE)
            m -= e.n_gape;
		if (m <= 0 || i == 0){
            if(m ==0 && i != 0){
                // extend the aln
                bwt_extend_exact(aux->bi_bwt, seq,(is_backward==0?(end+len-i+1):(start-len+i-1)) 
                        , &i, is_backward, &k, &l, &rev_k, &rev_l);
            }
            if(is_backward == 1 && max_pos >= start + i -len &&
                    aln->start > start+i-len) {
                //aln->start = i;
                // length changed, but start and end is origin coordinate
                aln->start = start+i-len;
                max_pos = aln->start;
            } else if(is_backward == 0 && max_pos <= end+len-i &&
                    aln->end < end+len-i) {
                aln->end = end+len-i;
                max_pos = aln->end;
            } else 
                continue;

            aln->k = k;
            aln->l = l;
            aln->type = BWA_TYPE_SPLICING;
            aln->rev_k = rev_k;
            aln->rev_l = rev_l;
            aln->n_mm = e.n_mm;
            aln->n_gapo = e.n_gapo;
            aln->n_gape = e.n_gape;
            aln->score = e.info>>21;
            if(i == 0){
                *_max_pos = max_pos;
                return 1;
            }
            continue;
        }

		--i;    //i：剩余reads长度(0 based)
        real_pos = is_backward == 1? start-len+i: len+end-i;
        if(is_backward == 1) {
            BWTAllSARangesBackward_Bidirection(aux->bi_bwt, k, l, rev_k, rev_l,
                    sa_k, sa_l, rev_sa_k, rev_sa_l);
        } else {
            BWTAllSARangesForward_Bidirection(aux->bi_bwt, k, l, rev_k, rev_l,
                    sa_k, sa_l, rev_sa_k, rev_sa_l);
        }
        // retrieve Occ values
		occ = l - k + 1;

        //TODO use width to reduce backtracing
        allow_diff = 1;
        {
            if(is_backward == 1 && max_pos  < real_pos && (width[real_pos].bid - width[max_pos].bid > m ||
                        (width[real_pos].bid - width[max_pos].bid == m &&
                         width[max_pos].bid != width[max_pos +1].bid)) )
                allow_diff = 0;
            if(is_backward == 0 && max_pos > real_pos && (width[real_pos].bid - width[max_pos].bid > m ||
                        (width[real_pos].bid - width[max_pos].bid== m &&
                         width[max_pos].bid != width[max_pos - 1].bid))){
                allow_diff = 0;
            }
        }

		// indels
        // TODO !!! 添加exact_extionsion, 当没有可用的max_diff时，用之，减少迭代次数
        tmp = (opt->mode & BWA_MODE_LOGGAP)? int_log2(e.n_gape + e.n_gapo)/2+1 : e.n_gapo + e.n_gape;
		if (allow_diff && i >= opt->indel_end_skip + tmp &&
                len - i >= opt->indel_end_skip + tmp) {
			if (e.state == STATE_M) { // gap open
				if (e.n_gapo < opt->max_gapo) { // gap open is allowed
					// insertion
					gap_push(stack, i, k, l, rev_k, rev_l,
                            e.n_mm, e.n_gapo + 1, e.n_gape, STATE_I, 1, opt);
					// deletion
					for (j = 0; j != 4; ++j) {
						k = sa_k[j];
						l = sa_l[j];
                        rev_k = rev_sa_k[j];
                        rev_l = rev_sa_l[j];
						if ((is_backward == 1 &&k <= l) || (is_backward == 0 && rev_k <= rev_l))
                            gap_push(stack, i + 1, k, l,rev_k, rev_l,
                                    e.n_mm, e.n_gapo + 1, e.n_gape, STATE_D, 1, opt);
					}
				}
			} else if (e.state == STATE_I) { // extention of an insertion
				if (e.n_gape < opt->max_gape) // gap extention is allowed
					gap_push(stack, i, k, l,rev_k, rev_l,
                            e.n_mm, e.n_gapo, e.n_gape + 1, STATE_I, 1, opt);
			} else if (e.state == STATE_D) { // extention of a deletion
				if (e.n_gape < opt->max_gape) { // gap extention is allowed
					if (e.n_gape + e.n_gapo < max_diff || occ < opt->max_del_occ) {
						for (j = 0; j != 4; ++j) {
						    k = sa_k[j];
						    l = sa_l[j];
                            rev_k = rev_sa_k[j];
                            rev_l = rev_sa_l[j];
							if (k <= l)
                                gap_push(stack, i + 1, k, l, rev_k, rev_l,
                                        e.n_mm, e.n_gapo, e.n_gape + 1, STATE_D, 1, opt);
						}
					}
				}
			}
		}
		// mismatches and exact
        if(allow_diff == 1){
  		    for (j = 1; j <= 4; ++j) {
                int c =(seq[real_pos] + j) & 3;
  		    	int is_mm = (j != 4 || seq[real_pos] > 3 );
  		    	k = sa_k[c];
  		    	l = sa_l[c];
                rev_k = rev_sa_k[c];
                rev_l = rev_sa_l[c];
  		    	if ((is_backward == 1 &&k <= l) || (is_backward == 0 && rev_k <= rev_l))
                    gap_push(stack, i, k, l, rev_k, rev_l, e.n_mm + is_mm, e.n_gapo, e.n_gape, STATE_M, is_mm, opt);
  		    }
        }
	}
    // when do not have an exact extension position
    if(*_max_pos != max_pos){
        *_max_pos = max_pos;
        return 2;
    } else
        return -1;
}

/*
 * find the potential splice site from have been mapped position
 * !! base on canonical or non-canonical motif
 * _left|......|_right
 * _left - _right : the range of motif may take place and !!!
 * _left or _right : 0 based and have been mapped, and relative to 
 *              origin seq
 * ext: n_gapo + n_gape
 * result: possible splice site
 */
static bwtint_t *splice_site_search_from_pos(const Idx2BWT *bi_bwt, int is_backward, int strand,
        bwtint_t _pos, ubyte_t *seq, int ext, int _left, int _right, int *site_num)
{
    int l = 0, ref_len = _right - _left - 1, i;
    bwtint_t k;
    ubyte_t *ref_seq;
    HSP *hsp = bi_bwt->hsp;
    int max_site_num = 10, j, ref_pos, seq_splice_pos;

    // init canonial and non-canonical motif form
    ubyte_t motif_posv[12] = {2,3,0,2,2,1,0,2,0,3,0,1};
    ubyte_t motif_neg[12] = {1,3,0,1,1,3,2,1,2,3,0,3};
    ubyte_t *motif;

    // the highest two bit save intron form type
    // positive motif: type 1: GTAG 2:GCAG 3: ATAC
    // negetive motif: type 1: CTAC 2:CTGC 3: GTAT

    bwtint_t *site_pos = (bwtint_t *)calloc(max_site_num, sizeof(bwtint_t));    
    bwtint_t *site_pos_new;

    ref_seq  = (ubyte_t *)calloc(ref_len, 1);
    if(is_backward == 1)
        _pos -= ref_len;
    else 
        _pos += _left + 1 + ext;

    for(k = _pos; k < _pos + ref_len && k < hsp->dnaLength; ++k)
        ref_seq[l++] = hsp->packedDNA[k>>4] >>((~k & 15) <<  1) & 3;

    *site_num = 0;
    for(i = 0; i < 3; ++i){
        motif = (strand == 0? motif_posv : motif_neg);
        motif = (is_backward == 0? (motif+4*i) : (motif+4*i+2));
        for(j = 2; j < ref_len-1; ++j)
        {
            ref_pos = (is_backward == 0? j : ref_len-2-j);
            seq_splice_pos = (is_backward == 0? _left+j:_right-j);
            if(ref_seq[ref_pos] == motif[0] && ref_seq[ref_pos+1] == motif[1])
            {
                if((is_backward == 0 && //ref_seq[ref_pos-2]== seq[seq_splice_pos-1] &&
                    ref_seq[ref_pos-1] == seq[seq_splice_pos]) ||
                    (is_backward == 1 && ref_seq[ref_pos+2]== seq[seq_splice_pos]))
                    //&&ref_seq[ref_pos+3] == seq[seq_splice_pos+1]))
                {
                    site_pos[*site_num] = (i << 30) | (seq_splice_pos);
                    *site_num += 1;
                    if(*site_num == max_site_num){
                        site_pos_new = (bwtint_t *)realloc(site_pos, sizeof(bwtint_t) * (max_site_num<<1));
                        if(site_pos_new != NULL){
                            site_pos = site_pos_new;
                            max_site_num <<= 1;
                        } else {
                            free(ref_seq);
                            return site_pos;
                        }
                    }
                }
            }
        }
        //if(i == 1 && (*site_num) != 0)
        //    break;

        /* the first two motif have the same fore part */
        if((is_backward == 0 && strand == 1 && i == 0)||
                (is_backward == 1 && strand == 0 && i == 0))
            i += 1;

    }
    free(ref_seq);
    return site_pos;
}

/*
 * check segment extension result
 * check by the intron end type if it is a pair of GT-AG GC-AG and AT-AC
 * type : intron type
 * result: motif type, 3: do not have match motif
 */
static inline int check_site_by_intron_end(const Idx2BWT *bi_bwt, bwt_aln1_t *aln,
        int is_backward, int type, int strand)
{
    ubyte_t ref[2];
    bwtint_t k, m;
    HSP *hsp = bi_bwt->hsp;
    unsigned int seq_id;
    bwtint_t ori_pos, occ_pos;
    ubyte_t motif_posv[12] = {2,3,0,2,2,1,0,2,0,3,0,1};
    ubyte_t motif_neg[12] = {1,3,0,1,1,3,2,1,2,3,0,3};
    ubyte_t *motif, *tmp_motif;
    motif = (strand == 0 ? motif_posv : motif_neg);
    for(m = aln->k; m <= aln->l; ++m){
        BWTRetrievePositionFromSAIndex(bi_bwt, m, &seq_id, &ori_pos, &occ_pos);
        k = (is_backward == 1? (occ_pos - 2) : (occ_pos+aln->end+1));
        ref[0] = hsp->packedDNA[k>>4] >>((~k & 15) <<  1) & 3;
        k += 1;
        ref[1] = hsp->packedDNA[k>>4] >>((~k & 15) <<  1) & 3;
        tmp_motif = (is_backward == 0 ? (motif+type*4) : (motif+type*4+2));
        if(((ref[0]^tmp_motif[0]) | (ref[1]^tmp_motif[1])) == 0)
            return type; 
        else if((type == 0)&&
                ((is_backward == 1 && strand == 1) || (is_backward == 0 && strand == 0)))
        {
            type = 1;
            tmp_motif = (is_backward == 0 ? (motif+type*4) : (motif+type*4+2));
            if(((ref[0]^tmp_motif[0]) | (ref[1]^tmp_motif[1])) == 0)
                return type;
        }
        else
            continue;
    }
    return 3;
}

/*
 * extend aln backward, extend at least to pos *_left
 */
int bwt_extend_backward(bwt_aux_t *aux, bwt_aln1_t *aln, int *_left)
{
    int b_aln;
    gap_reset_stack(aux->stack);
    gap_push(aux->stack, aux->len, aln->k,aln->l,aln->rev_k, aln->rev_l,
            aln->n_mm, aln->n_gapo, aln->n_gape,
            0, 0, aux->opt);
    b_aln = bwt_backtracing_search(aux, aln, 1, _left);
    return b_aln;
}

/*
 * extend aln foreward, extend at least to pos *_right
 */
int bwt_extend_foreward(bwt_aux_t *aux, bwt_aln1_t *aln, int *_right)
{
    int b_aln;
    gap_reset_stack(aux->stack);
	gap_push(aux->stack, aux->len, aln->k, aln->l, aln->rev_k, aln->rev_l,
            aln->n_mm, aln->n_gapo, aln->n_gape,
            0, 0, aux->opt);
    b_aln = bwt_backtracing_search(aux, aln, 0, _right);
    return b_aln;
}

/*
 * seed result check, to find the most correlated aln
 * result: bool
 */
static int bwt_aln_corelate_check(Idx2BWT *bi_bwt, bwt_seed_aln_t *p, bwt_seed_aln_t* q)
{
    int i, k, tot_cnt, max_cnt, cur_cnt, seq_id;
    bwtint_t j, ori_pos, min_dist, res_pos, occ_pos, sa_fore, sa_back;
    bwt_pos_t *pos_info, *tmp_p, *pos_info_new;
    bwt_aln1_t *aln, *tmp_res1, *tmp_res2;
    tot_cnt = 0; cur_cnt = 0; max_cnt = 20;
    pos_info = (bwt_pos_t *)calloc(max_cnt, sizeof(bwt_pos_t));
    tmp_p = pos_info;
    min_dist = 0xffffffff;
    res_pos = 0xffffffff;
    tmp_res1 = (bwt_aln1_t *)malloc(sizeof(bwt_aln1_t));
    tmp_res2 = (bwt_aln1_t *)malloc(sizeof(bwt_aln1_t));

    // only the first 10 aln is validation
    for(i = 0; i < p->n_aln && i < 10; ++i) {
        aln = p->aln + i;
        tot_cnt += aln->l - aln->k + 1 >10? 10: aln->l - aln->k + 1;
        if(tot_cnt > max_cnt){
            max_cnt <<= 1;
            pos_info_new = (bwt_pos_t *)realloc(pos_info, max_cnt * sizeof(bwt_pos_t));
            if(pos_info_new != NULL){
                pos_info = pos_info_new;
            }
            else
                break;
            tmp_p = pos_info + cur_cnt;
        }
        //TODO  only 10 cnt will retrieved
        for(j = aln->k; j <= aln->l && j < aln->k + 10; ++j){
            BWTRetrievePositionFromSAIndex(bi_bwt,j, &seq_id, &tmp_p->ori_pos, &tmp_p->occ_pos);
            tmp_p->seq_id = seq_id;
            //tmp_p->sa = j;
            tmp_p->r_aln = i;
            tmp_p += 1;
            cur_cnt += 1;
        }
    }
    int t_r1= 0, t_r2 = 0;
    int dist = min_dist;
    for(i = 0; i < q->n_aln; ++i){
        aln = q->aln+i;
        for(j = aln->k; j <= aln->l && j < aln->k+50; ++j){
            BWTRetrievePositionFromSAIndex(bi_bwt, j, &seq_id, &ori_pos, &occ_pos);
            for(k = 0; k < tot_cnt; ++k){
                dist = occ_pos- (pos_info+k)->occ_pos;
                if((pos_info+k)->seq_id == seq_id && dist > 50 &&
                        dist < 50000 && min_dist > dist)
                {
                    min_dist = dist;
                    //sa_back = j;
                    //sa_fore = (pos_info+k)->sa;
                    res_pos = (pos_info + k)->occ_pos;
                    t_r1 = (pos_info + k)->r_aln;
                    t_r2 = i;
                }
            }
        }
    }
    if(min_dist != 0xffffffff){
        memcpy(tmp_res1, p->aln + t_r1, sizeof(bwt_aln1_t));
        memcpy(tmp_res2, q->aln + t_r2, sizeof(bwt_aln1_t));
        //tmp_res1->k = tmp_res1->l = sa_fore;
        //tmp_res2->k = tmp_res2->l = sa_back;
        p->n_aln = 1; free(p->aln); p->aln = tmp_res1;
        q->n_aln = 1; free(q->aln); q->aln = tmp_res2;
    } else {
        free(tmp_res1);
        free(tmp_res2);
    }
    free(pos_info);
    // TODO free do not used aln
    return res_pos;
}

/*
 * the function use seed-extend strategy mapping splicing reads to ref
 * input: all things needs are put in bwt_aux_t
 */
bwt_aln1_t *bwt_splice_match(bwt_aux_t *aux, int *_n_aln)
{
    int i, j, _j, seed_num, n_aln, max_aln_n;
    bwtint_t seq_id=0, offset=0;
    bwt_seed_aln_t *multi_seed_aln, *q;
    bwt_aux_t *aux_seed, *aux_ext;

    // init aux struct for seed and extend
    aux_seed = (bwt_aux_t *)malloc(sizeof(bwt_aux_t));
    memcpy(aux_seed, aux, sizeof(bwt_aux_t));
    aux_ext = (bwt_aux_t *)malloc(sizeof(bwt_aux_t));
    memcpy(aux_ext, aux, sizeof(bwt_aux_t));

    // TODO make seed
    int seed_len, len_align, len = aux->len;   // origin seq length
    seed_num = 3;
    seed_len = len/seed_num;
    //aux_seed->len = aux_seed->opt->seed_len = seed_len;
    max_aln_n = seed_num * 2;
 
    // initiate seed opt, do not allow gap
    aux_seed->opt = (gap_opt_t *)malloc(sizeof(gap_opt_t));
    memcpy(aux_seed->opt, aux->opt, sizeof(gap_opt_t));
    aux_seed->opt->mode &= ~BWA_MODE_GAPE;   //do not allow gap
    aux_seed->opt->max_gapo = 0;
    aux_seed->opt->max_gape = 0;
    aux_seed->opt->max_diff = aux->opt->max_seed_diff;

    // init extension opt
    aux_ext->opt = (gap_opt_t *)malloc(sizeof(gap_opt_t));
    memcpy(aux_ext->opt, aux->opt, sizeof(gap_opt_t));
    //aux_ext->opt->mode &= ~BWA_MODE_GAPE;   //do not allow gap
    //aux_ext->opt->max_gapo = 0;
    //aux_ext->opt->max_gape = 0;
    aux_ext->opt->max_gape = 3;

    bwt_array_t *arr = aux->arr;

    // initiate multi seed aln, which save aln result
    multi_seed_aln = (bwt_seed_aln_t *)calloc(seed_num * 2, sizeof(bwt_seed_aln_t));
    int strand = 0; // 1 : reverse complement
    // mapping seed to reference
    _j = 0; // saving total num of seed which can be mapped to ref
    uint8_t seg_mtype = 0;

    // because width is already generated for common seq
    // mapping seed to ref, for both positive and reverse complement seq
    bwtint_t seq_pos;
    // seq_pos is an occ type position
    for(i = 0; i < seed_num * 2; ++i){
        q = multi_seed_aln + i;
        aux_seed->strand = q->strand = i/3;
        len_align = seed_len;
        len_align += (i%3 == 2? (len%3) : 0);
        aux_seed->len = aux_seed->opt->seed_len = len_align;

        memset(aux_seed->width_seed, 0, sizeof(bwt_width_t)*(aux->max_len+1));
        (i < 3)? (aux_seed->seq = aux->seq + (i%3) * seed_len) : 
            (aux_seed->rc_seq = aux->rc_seq + (i%3)*seed_len);
        bwt_cal_width(aux->bi_bwt, len_align, (q->strand == 0)? aux->seq : aux->rc_seq,
                aux_seed->width_seed, 1);
        aux_seed->width_back = aux_seed->width_seed;

        gap_reset_stack(aux->stack);
        q->aln = bwt_match_gap(aux_seed, &q->n_aln);
        if(q->n_aln != 0) {
            ++_j;
            seg_mtype += (1 <<(i%3));
            for(j = 0; j < q->n_aln; ++j){
                (q->aln+j)->start = (i%3) * seed_len;
                (q->aln+j)->end = (q->aln+j)->start+len_align-1;
            }
        }
        if(i == 1 && _j == 0){
            i += 1; _j = 0; seg_mtype = 0;
            continue;
        }
        if((i == 2 || i == 5) && _j > 1){
            strand = i < 3? 0:1;
            q = multi_seed_aln+ strand *3;
            switch (seg_mtype)
            {
                case 5:
                case 7:
                    seq_pos = bwt_aln_corelate_check(aux->bi_bwt, q, q+2);
                    break;
                case 3:
                    seq_pos = bwt_aln_corelate_check(aux->bi_bwt, q , q + 1);
                    break;
                case 6:
                    seq_pos = bwt_aln_corelate_check(aux->bi_bwt, q+1, q+2);
                    break;
            }
            if(seq_pos != 0xffffffff)
                break;
        }
        if(i ==2)
        {
            _j = 0; seg_mtype = 0;
        }
    }

    // check segment mapping result, if both positive and reverse complement
    // cannot meet a condition, dicard the read.
    *_n_aln = 0;
    bwt_aln1_t *res_aln, *tmp_aln_motif; //res_aln: before the final result, it save tmp aln results;
    res_aln = (bwt_aln1_t *)calloc(2, sizeof(bwt_aln1_t));
    bwtint_t *site_pos, occ_of_first;
    site_pos = (bwtint_t *)calloc(1, sizeof(bwtint_t));
    if(_j < 2 || seq_pos == 0xffffffff)
        goto end;

/*
 * seed alignment have four consideration, retrieve sa to pos to check correlation, 
 * then extend seed, report splice site
 */
    ubyte_t *seq;
    if(strand == 0){
        seq = aux->seq;
        bwt_cal_width(aux->bi_bwt, aux->len, seq, aux->width_back, 1);
        bwt_cal_width(aux->bi_bwt, aux->len, seq, aux->width_fore, 0);
    } else {
        seq = aux->rc_seq;
        bwt_cal_width(aux->bi_bwt, aux->len, seq, aux->width_back, 1);
        bwt_cal_width(aux->bi_bwt, aux->len, seq, aux->width_fore, 0);
    }

    // base on the hit count
    aux->strand = aux_ext->strand = strand;
    int max_pos, m, motif_type_checked, split_pos, b_ext, b_mapping_stat, motif_type;
    int site_num = 0;
    b_mapping_stat = 0;
    aux_ext->seq = aux->seq; aux_ext->rc_seq = aux->rc_seq;
    if(seg_mtype == 3) {
        // first two seed mapping to ref, extend first seed,
        // do not find suitable aln pair
        //if(seq_pos != 0xffffffff) {
        // look up possible splicing site by mapping record
        split_pos = bwt_find_split_pos_by_record(arr, seq_pos + q->aln->end,
                aux->len - q->aln->end -1, 0);
        if(split_pos != -1){
            site_num = 1;
            split_pos = q->aln->end + split_pos;
            *site_pos = ((arr->data[arr->last_mid * 2 + 1]<<1) & 0xd0000000) | split_pos;
        } else {
            // look up possible splicing site by intron type
            free(site_pos);
            site_pos = splice_site_search_from_pos(aux->bi_bwt, 0, strand, seq_pos, seq,
                    q->aln->n_gapo+ q->aln->n_gape,
                    (q+1)->aln->end, len, &site_num);
        }
        
        // extend first part to the end of second part
        max_pos = (q+1)->aln->end;
        aux_ext->len = max_pos - (q->aln->end);
        b_ext = bwt_extend_foreward(aux_ext, q->aln, &max_pos);
        if(b_ext != 1)
            goto end;

        // extend first seed, extend length at least accross second seed
        motif_type = site_pos[0]>>30;
        memcpy(res_aln, q->aln, sizeof(bwt_aln1_t));
        {                                   /* get the alignment of last 8 chars */
            aux_ext->len = 12; 
            strand == 0? (aux_ext->seq = aux->seq+len-12):(aux_ext->rc_seq = aux->rc_seq+len-12);
            bwt_width_t *width_tmp = (bwt_width_t *)calloc(aux_ext->len+1, sizeof(bwt_width_t));
            aux_ext->width_back = width_tmp;
            bwt_cal_width(aux->bi_bwt, aux_ext->len, strand==0? aux_ext->seq:aux_ext->rc_seq,
                    aux_ext->width_back, 1);
            aux_ext->width_seed = NULL;
            free((q+2)->aln);
            (q+2)->aln = bwt_match_gap(aux_ext, &(q+2)->n_aln);
            if((q+2)->n_aln == 0)
                goto end;
            else
            {
                for(j = 0; j < (q+2)->n_aln; ++j){
                    (q+2)->aln->start = len-12;
                    (q+2)->aln->end = aux->len -1;
                }
            }
            free(width_tmp);
            aux_ext->width_back = aux->width_back;
            strand == 0?(aux_ext->seq = aux->seq) : (aux_ext->rc_seq = aux->rc_seq);
        }
        occ_of_first = bwt_aln_corelate_check(aux->bi_bwt, q , q + 2);
        if(occ_of_first == 0x3fffffff)
            goto end;
        memcpy(res_aln+1, (q+2)->aln, sizeof(bwt_aln1_t));
        for(m = 0; m < site_num; ++m) {
            max_pos = site_pos[m] & 0x3fffffff;
            if(motif_type != site_pos[m]>>30){
                memcpy(res_aln, q->aln, sizeof(bwt_aln1_t));
                motif_type = site_pos[m]>>30;
            }
            aux_ext->len = max_pos - q->aln->end;
            b_ext = bwt_extend_foreward(aux_ext, q->aln, &max_pos);
            // if it can not extend to max_pos, break
            if(b_ext != 1) {
                while((m+1<site_num) && ((site_pos[m+1]>>30) == motif_type))
                    m++;
                continue;
            }
            // extend leaving part to ref
            max_pos += 1;
            aux_ext->len = (res_aln+1)->start - max_pos;
            // only leaving part seq length great than 8 bp
            // TODO could set the parameter by user
            if(aux_ext->len > 0) {
                memcpy((q+2)->aln, res_aln+1, sizeof(bwt_aln1_t));
                b_ext = bwt_extend_backward(aux_ext, (q+2)->aln, &max_pos);

                if(b_ext != 1)
                    continue;
                else
                {
                    // filter result by distance 
                    seq_pos = bwt_aln_corelate_check(aux->bi_bwt, q , q + 2);
                    if(seq_pos == 0xffffffff)
                        continue;
                    // check the aln can meet the common intron type
                    motif_type_checked = check_site_by_intron_end(aux->bi_bwt, (q+2)->aln,
                            1, motif_type, strand);
                    if(motif_type_checked != 3){
                        *_n_aln = 2;
                        b_mapping_stat = 1;
                        // if only one cnt, insert it to array
                        // TODO record two side of intron positon info
                        if(res_aln->k == res_aln->l)
                            bwt_array_insert(arr, seq_pos + max_pos, motif_type_checked);
                        break;
                    } else
                        continue;
                }
            } else {
                // if leaving length less than 8, do not map leaving part,
                // and save the first part;
                *_n_aln = 1;
                b_mapping_stat = 1;
                continue;
            }
        }
        //} else {
        //    b_mapping_stat = 1;
        //}
        // if it do not have common intron end type
        // or mapping to ref use certain introns failed
        if(b_mapping_stat == 0) {
            //strand == 0? (aux_ext->seq = aux->seq): (aux_ext->rc_seq = aux->rc_seq);
            memcpy((q+2)->aln, res_aln+1, sizeof(bwt_aln1_t));
            aux_ext->len = (q+2)->aln->start - q->aln->end - 1;
            max_pos = (q+1)->aln->end + 1;
            b_ext = bwt_extend_foreward(aux_ext, q->aln, &max_pos);
            *_n_aln = 1;
            if(b_ext == 2){
                // TODO trim extension result
                // mapping leaving part to ref
                max_pos += 1;
                aux_ext->len = (q+2)->aln->start - max_pos;
                if(aux_ext->len > 0)
                {
                    b_ext = bwt_extend_backward(aux_ext, (q+2)->aln, &max_pos);
                    if(b_ext == 1) {
                        (q+2)->aln->start = max_pos;
                        // filter result by distance 
                        seq_pos = bwt_aln_corelate_check(aux->bi_bwt, q , q + 2);
                        if(seq_pos != 0xffffffff){
                            *_n_aln = 2;
                            //insert seq_pos to array
                        }
                    }
                }
            }
        }
        // write result to res_aln
        if(*_n_aln != 0) {
            memcpy(res_aln, q->aln, sizeof(bwt_aln1_t));
            res_aln->type = BWA_TYPE_SPLICING;
        }
        if(*_n_aln == 2) {
            memcpy(res_aln+1, (q+2)->aln, sizeof(bwt_aln1_t));
            (res_aln+1)->type = BWA_TYPE_SPLICING;
            bwt_array_insert(arr, seq_pos + max_pos, 3);
        }
    } else if(seg_mtype == 5 || seg_mtype == 7){
        // first and last seed can mapping to ref, check aln, and extend
        //if(seq_pos != 0xffffffff) {
        // foreward, split_pos is relative to q->aln->end;
        split_pos = bwt_find_split_pos_by_record(arr, seq_pos + q->aln->end,
                (q+2)->aln->start - q->aln->end, 0);
        if(split_pos != -1){
            site_num = 1;
            split_pos += q->aln->end;
            *site_pos = ((arr->data[arr->last_mid*2 + 1]<<1)& 0xd0000000) | split_pos;
        } else {
            // look up possible splicing site
            free(site_pos);
            site_pos = splice_site_search_from_pos(aux->bi_bwt, 0,strand, seq_pos, seq,
                    q->aln->n_gapo+ q->aln->n_gape,
                    q->aln->end, (q+2)->aln->start, &site_num);
        }

        // if the last part of seq have not been mapped
        //if((q+2)->aln->end < aux->len-1){
        //    max_pos = aux->len - 1;
        //    aux_ext->len = aux->len;
        //    b_ext = bwt_extend_foreward(aux_ext, (q+2)->aln, &max_pos);
        //}

        // tmp_aln: for forward extension
        //bwt_aln1_t *tmp_aln;
        //tmp_aln = (bwt_aln1_t *)malloc(sizeof(bwt_aln1_t));
        //memset(tmp_aln, 0, sizeof(bwt_aln1_t));
        // TODO if site_num == 0 ? extend and trim end
        aux_ext->opt->mode &= ~BWA_MODE_GAPE;   //do not allow gap
        motif_type = site_pos[0]>>30;
        memcpy(res_aln, q->aln, sizeof(bwt_aln1_t));
        memcpy(res_aln+1, (q+2)->aln, sizeof(bwt_aln1_t));
        for(m = 0; m < site_num; ++m){
            max_pos = site_pos[m] & 0x3fffffff;
            if(motif_type != site_pos[m]>>30){
                memcpy(q->aln, res_aln, sizeof(bwt_aln1_t));
                motif_type = site_pos[m]>>30;
            }
            memcpy((q+2)->aln, res_aln+1, sizeof(bwt_aln1_t));

            aux_ext->len = max_pos - q->aln->end;
            //strand == 0? (aux_ext->seq = aux->seq): (aux_ext->rc_seq = aux->rc_seq);
            b_ext = bwt_extend_foreward(aux_ext, q->aln, &max_pos);
            if(b_ext != 1) {
                while((m+1<site_num) && ((site_pos[m+1]>>30) == motif_type))
                    m++;
                continue;
            }
            max_pos = max_pos + 1;
            //(strand == 0) ? (aux_ext->seq = aux->seq + max_pos) : 
            //    (aux_ext->rc_seq = aux->rc_seq + max_pos);
            // len is leaving length for mapping
            aux_ext->len = (q+2)->aln->start - max_pos;
            //memcpy(tmp_aln, (q+2)->aln, sizeof(bwt_aln1_t));
            b_ext = bwt_extend_backward(aux_ext, (q+2)->aln, &max_pos);
            // TODO check whether the aln is changed when extension failed 
            if(b_ext != 1)
                continue;
            // check extend result
            motif_type_checked = check_site_by_intron_end(aux->bi_bwt, (q+2)->aln,
                    1, motif_type, strand);
            if(motif_type_checked != 3){
                *_n_aln = 2;
                b_mapping_stat = 1;
                //// if only one cnt, insert it to array
                //// TODO record two side of intron positon info
                //if(res_aln->k == res_aln->l)
                //    bwt_array_insert(arr, seq_pos + max_pos, motif_type_checked);
                bwt_array_insert(arr, seq_pos + max_pos, motif_type_checked);
                break;
            } else
                continue;
            //*_n_aln = 2;
            //memcpy((q+2)->aln, tmp_aln, sizeof(bwt_aln1_t));
            //b_mapping_stat = 1;
            //break;
        }
        //free(tmp_aln);
        //} else {
        //    *_n_aln = 0;
        //    b_mapping_stat = 1;
        //}
        // if it do not have common intron end type
        // or mapping to ref use certain introns failed
        if(b_mapping_stat == 0){
            aux_ext->opt->mode &= ~BWA_MODE_GAPE;   //do not allow gap
            // extens the first part of the seq
            max_pos = q->aln->end + 1;
            aux_ext->len = (q+2)->aln->start - q->aln->end - 1;
            //strand == 0? (aux_ext->seq = aux->seq): (aux_ext->rc_seq = aux->rc_seq);
            b_ext = bwt_extend_foreward(aux_ext, q->aln, &max_pos);
            if(b_ext == -1)
                *_n_aln = 0;
            else if (b_ext == 1)
                *_n_aln = 2;
            else {
                // extend leaving part to q->aln->end
                // allow 3 indels
                max_pos = q->aln->end + 1;
                aux_ext->len = (q+2)->aln->start - max_pos;
                // allow gap
                aux_ext->opt->mode |= BWA_MODE_GAPE;
                //(strand == 0) ? (aux_ext->seq = aux->seq+q->aln->end+1) : 
                //    (aux_ext->rc_seq = aux->rc_seq + max_pos);
                b_ext = bwt_extend_backward(aux_ext, (q+2)->aln, &max_pos);
                seq_pos = bwt_aln_corelate_check(aux->bi_bwt, q , q + 2);
                if(b_ext == -1 || seq_pos == 0xffffffff)
                    *_n_aln = 0;
                else
                    *_n_aln = 2;
            }
        }
        if(*_n_aln != 0) {
            memcpy(res_aln, q->aln, sizeof(bwt_aln1_t));
            memcpy(res_aln+1, (q+2)->aln, sizeof(bwt_aln1_t));
            bwt_array_insert(arr, seq_pos + max_pos, 3);
            res_aln->type = BWA_TYPE_SPLICING;
            (res_aln+1)->type = BWA_TYPE_SPLICING;
        }
    } else if(seg_mtype == 6){
        //last two seed can mapping to ref
        //if(seq_pos != 0xffffffff){
        split_pos = bwt_find_split_pos_by_record(arr, seq_pos, (q+2)->aln->start, 1);
        if(split_pos != -1){
            site_num = 1;
            split_pos = (q+2)->aln->start - split_pos;
            *site_pos = ((arr->data[arr->last_mid*2 + 1]<<1) & 0xd0000000) | split_pos;
        } else {
            free(site_pos);
            site_pos = splice_site_search_from_pos(aux->bi_bwt, 1, strand,
                    seq_pos - seed_len, seq, 0, 0, (q+1)->aln->start, &site_num);
        }

        // extend the third part to the beginning of second part
        {
            aux_ext->len = (q+2)->aln->start - (q+1)->aln->start;
            max_pos = (q+1)->aln->start;
            //strand == 0? aux_ext->seq = aux->seq+max_pos : 
            //    aux_ext->rc_seq = aux->rc_seq+max_pos;
            b_ext = bwt_extend_backward(aux_ext, (q+2)->aln, &max_pos);
            if(b_ext != 1)
                goto end;
        }

        // if the last part of seq have not been mapped
        //if((q+2)->aln->end < aux->len-1){
        //    max_pos = aux->len - 1;
        //    aux_ext->len = aux->len;
        //    b_ext = bwt_extend_foreward(aux_ext, (q+2)->aln, &max_pos);
        //}

        // extend first seed, extend length at least accross second seed
        motif_type = site_pos[0]>>30;
        aux_ext->len = 12; 
        aux_ext->width_back = aux->width_back;
        aux_ext->width_seed = NULL;
        strand == 0? (aux_ext->seq = aux->seq):(aux_ext->rc_seq = aux->rc_seq);
        free(q->aln);
        q->aln = bwt_match_gap(aux_ext, &q->n_aln);
        if(q->n_aln == 0)
            goto end;
        else
        {
            for(j = 0; j < q->n_aln; ++j){
                q->aln->start = 0;
                q->aln->end = 11;
            }
        }
        occ_of_first = bwt_aln_corelate_check(aux->bi_bwt, q , q + 2);
        if(occ_of_first == 0x3fffffff)
            goto end;
        memcpy(res_aln, q->aln, sizeof(bwt_aln1_t));
        memcpy(res_aln+1, (q+2)->aln, sizeof(bwt_aln1_t));
        motif_type = site_pos[0]>>30;
        for(m = 0; m < site_num; ++m){
            max_pos = site_pos[m] & 0x3fffffff;
            aux_ext->len = (q+2)->aln->start - max_pos;
            if(motif_type != site_pos[m]>>30){
                memcpy((q+2)->aln, res_aln+1, sizeof(bwt_aln1_t));
                motif_type = site_pos[m]>>30;
            }
            //(strand == 0) ? (aux_ext->seq = aux->seq + max_pos) : 
            //    (aux_ext->rc_seq = aux->rc_seq + max_pos);
            //aux_ext->width_back = aux->width_back+max_pos;
            //aux_ext->width_seed = NULL;
            b_ext = bwt_extend_backward(aux_ext, (q+2)->aln, &max_pos);
            if(b_ext != 1){
                while((m+1<site_num) && ((site_pos[m+1]>>30) == motif_type))
                    m++;
                continue;
            }
            // mapping leaving part to ref
            max_pos -= 1;
            aux_ext->len = max_pos - q->aln->end;
            // only leaving part seq length gt 12 bp
            // TODO could set the parameter by user
            if(aux_ext->len > 0){
                //strand == 0?(aux_ext->seq = aux->seq+ max_pos) : 
                //    (aux_ext->rc_seq = aux->rc_seq+ max_pos);
                //aux_ext->width_seed = NULL;
                memcpy(q->aln, res_aln, sizeof(bwt_aln1_t));
                b_ext = bwt_extend_foreward(aux_ext, q->aln, &max_pos);
                if(b_ext != 1)
                    continue;
                else
                {
                    // filter result by distance 
                    seq_pos = bwt_aln_corelate_check(aux->bi_bwt, q , q + 2);
                    if(seq_pos == 0xffffffff)
                        continue;
                    // check the aln can meet the common intron type
                    motif_type_checked = check_site_by_intron_end(aux->bi_bwt, q->aln,
                            0, motif_type, strand);
                    if(motif_type_checked != 3){
                        b_mapping_stat = 1;
                        *_n_aln = 2;

                        // if only one cnt, insert it to array
                        if((q+2)->aln->k == (q+2)->aln->l)
                            bwt_array_insert(arr, seq_pos + max_pos, site_pos[m]>>30);
                        break;
                    } else
                        continue;
                }
            } else {
                // if leaving length less than 12, do not map leaving part,
                *_n_aln = 1;
                b_mapping_stat = 1;
                break;
            }
        }
        //} else {
        //    b_mapping_stat = 1;
        //}
        // if it do not have common intron end type
        // or mapping to ref use certain introns failed
        if(b_mapping_stat == 0){
            aux_ext->len = (q+2)->aln->start - q->aln->end - 1;
            max_pos = (q+2)->aln->start - 1;
            //strand == 0? (aux_ext->seq = aux->seq): (aux_ext->rc_seq = aux->rc_seq);
            //(strand == 0) ? (aux_ext->seq = aux->seq + ) : 
            //    (aux_ext->rc_seq = aux->rc_seq + max_pos);
            b_ext = bwt_extend_backward(aux_ext, (q+2)->aln, &max_pos);
            if(b_ext == 2)
                *_n_aln = 1;
            // TODO trim extension result
            // mapping leaving part to ref
            max_pos -= 1;
            aux_ext->len = max_pos - q->aln->end;

            // leaving part should be gt 12 bp
            if(aux_ext->len > 0){
                //strand == 0? aux_ext->seq = aux->seq: aux_ext->rc_seq = aux->rc_seq;
                b_ext = bwt_extend_foreward(aux_ext, q->aln, &max_pos);
                if(b_ext == 1){
                    // filter result by distance 
                    seq_pos = bwt_aln_corelate_check(aux->bi_bwt, q , q + 2);
                    if(seq_pos != 0xffffffff) {
                        //insert seq_pos to array
                        *_n_aln = 2;
                        bwt_array_insert(arr, seq_pos + max_pos, 3);// intron type 3
                    }
                }
            }
        }
        if(*_n_aln == 2) {
            memcpy(res_aln, q->aln, sizeof(bwt_aln1_t));
            memcpy(res_aln+1, (q+2)->aln, sizeof(bwt_aln1_t));
            res_aln->type = BWA_TYPE_SPLICING;
            (res_aln+1)->type = BWA_TYPE_SPLICING;
        }
        else if(*_n_aln > 0) {
            memcpy(res_aln, (q+2)->aln, sizeof(bwt_aln1_t));
            res_aln->type = BWA_TYPE_SPLICING;
        }
    } else {

    }

	// clean up the unused data in the record
end:
    free(site_pos);
    // free aln which saved in multi_seed_aln
    q = multi_seed_aln;
    for(i = 0; i < seed_num * 2; ++i){
        if(q->aln != NULL){
            free(q->aln);
            q->aln = NULL;
        }
        q += 1;
    }
    free(multi_seed_aln);
    free(aux_seed->opt);
    free(aux_ext->opt);
    free(aux_seed);
    free(aux_ext);

    return res_aln;
}
