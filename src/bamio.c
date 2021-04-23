#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <stdio.h>
#include <errno.h>
#include <unistd.h>
#include <limits.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "bamio.h"

/**********************
 * beging of khashl.h *
 **********************/
#ifndef __AC_KHASHL_H
#define __AC_KHASHL_H

#define AC_VERSION_KHASHL_H "0.1"

/************************************
 * Compiler specific configurations *
 ************************************/

#if UINT_MAX == 0xffffffffu
typedef unsigned int khint32_t;
#elif ULONG_MAX == 0xffffffffu
typedef unsigned long khint32_t;
#endif

#if ULONG_MAX == ULLONG_MAX
typedef unsigned long khint64_t;
#else
typedef unsigned long long khint64_t;
#endif

#ifndef kh_inline
#ifdef _MSC_VER
#define kh_inline __inline
#else
#define kh_inline inline
#endif
#endif /* kh_inline */

#ifndef klib_unused
#if (defined __clang__ && __clang_major__ >= 3) || (defined __GNUC__ && __GNUC__ >= 3)
#define klib_unused __attribute__ ((__unused__))
#else
#define klib_unused
#endif
#endif /* klib_unused */

#define KH_LOCAL static kh_inline klib_unused

typedef khint32_t khint_t;

/******************
 * malloc aliases *
 ******************/

#ifndef kcalloc
#define kcalloc(N,Z) calloc(N,Z)
#endif
#ifndef kmalloc
#define kmalloc(Z) malloc(Z)
#endif
#ifndef krealloc
#define krealloc(P,Z) realloc(P,Z)
#endif
#ifndef kfree
#define kfree(P) free(P)
#endif

/****************************
 * Simple private functions *
 ****************************/

#define __kh_used(flag, i)       (flag[i>>5] >> (i&0x1fU) & 1U)
#define __kh_set_used(flag, i)   (flag[i>>5] |= 1U<<(i&0x1fU))
#define __kh_set_unused(flag, i) (flag[i>>5] &= ~(1U<<(i&0x1fU)))

#define __kh_fsize(m) ((m) < 32? 1 : (m)>>5)

static kh_inline khint_t __kh_h2b(khint_t hash, khint_t bits) { return hash * 2654435769U >> (32 - bits); }

/*******************
 * Hash table base *
 *******************/

#define __KHASHL_TYPE(HType, khkey_t) \
    typedef struct HType { \
        khint_t bits, count; \
        khint32_t *used; \
        khkey_t *keys; \
    } HType;

#define __KHASHL_PROTOTYPES(HType, prefix, khkey_t) \
    extern HType *prefix##_init(void); \
    extern void prefix##_destroy(HType *h); \
    extern void prefix##_clear(HType *h); \
    extern khint_t prefix##_getp(const HType *h, const khkey_t *key); \
    extern int prefix##_resize(HType *h, khint_t new_n_buckets); \
    extern khint_t prefix##_putp(HType *h, const khkey_t *key, int *absent); \
    extern void prefix##_del(HType *h, khint_t k);

#define __KHASHL_IMPL_BASIC(SCOPE, HType, prefix) \
    SCOPE HType *prefix##_init(void) { \
        return (HType*)kcalloc(1, sizeof(HType)); \
    } \
    SCOPE void prefix##_destroy(HType *h) { \
        if (!h) return; \
        kfree((void *)h->keys); kfree(h->used); \
        kfree(h); \
    } \
    SCOPE void prefix##_clear(HType *h) { \
        if (h && h->used) { \
            uint32_t n_buckets = 1U << h->bits; \
            memset(h->used, 0, __kh_fsize(n_buckets) * sizeof(khint32_t)); \
            h->count = 0; \
        } \
    }

#define __KHASHL_IMPL_GET(SCOPE, HType, prefix, khkey_t, __hash_fn, __hash_eq) \
    SCOPE khint_t prefix##_getp(const HType *h, const khkey_t *key) { \
        khint_t i, last, n_buckets, mask; \
        if (h->keys == 0) return 0; \
        n_buckets = 1U << h->bits; \
        mask = n_buckets - 1U; \
        i = last = __kh_h2b(__hash_fn(*key), h->bits); \
        while (__kh_used(h->used, i) && !__hash_eq(h->keys[i], *key)) { \
            i = (i + 1U) & mask; \
            if (i == last) return n_buckets; \
        } \
        return !__kh_used(h->used, i)? n_buckets : i; \
    } \
    SCOPE khint_t prefix##_get(const HType *h, khkey_t key) { return prefix##_getp(h, &key); }

#define __KHASHL_IMPL_RESIZE(SCOPE, HType, prefix, khkey_t, __hash_fn, __hash_eq) \
    SCOPE int prefix##_resize(HType *h, khint_t new_n_buckets) { \
        khint32_t *new_used = 0; \
        khint_t j = 0, x = new_n_buckets, n_buckets, new_bits, new_mask; \
        while ((x >>= 1) != 0) ++j; \
        if (new_n_buckets & (new_n_buckets - 1)) ++j; \
        new_bits = j > 2? j : 2; \
        new_n_buckets = 1U << new_bits; \
        if (h->count > (new_n_buckets>>1) + (new_n_buckets>>2)) return 0; /* requested size is too small */ \
        new_used = (khint32_t*)kmalloc(__kh_fsize(new_n_buckets) * sizeof(khint32_t)); \
        memset(new_used, 0, __kh_fsize(new_n_buckets) * sizeof(khint32_t)); \
        if (!new_used) return -1; /* not enough memory */ \
        n_buckets = h->keys? 1U<<h->bits : 0U; \
        if (n_buckets < new_n_buckets) { /* expand */ \
            khkey_t *new_keys = (khkey_t*)krealloc((void*)h->keys, new_n_buckets * sizeof(khkey_t)); \
            if (!new_keys) { kfree(new_used); return -1; } \
            h->keys = new_keys; \
        } /* otherwise shrink */ \
        new_mask = new_n_buckets - 1; \
        for (j = 0; j != n_buckets; ++j) { \
            khkey_t key; \
            if (!__kh_used(h->used, j)) continue; \
            key = h->keys[j]; \
            __kh_set_unused(h->used, j); \
            while (1) { /* kick-out process; sort of like in Cuckoo hashing */ \
                khint_t i; \
                i = __kh_h2b(__hash_fn(key), new_bits); \
                while (__kh_used(new_used, i)) i = (i + 1) & new_mask; \
                __kh_set_used(new_used, i); \
                if (i < n_buckets && __kh_used(h->used, i)) { /* kick out the existing element */ \
                    { khkey_t tmp = h->keys[i]; h->keys[i] = key; key = tmp; } \
                    __kh_set_unused(h->used, i); /* mark it as deleted in the old hash table */ \
                } else { /* write the element and jump out of the loop */ \
                    h->keys[i] = key; \
                    break; \
                } \
            } \
        } \
        if (n_buckets > new_n_buckets) /* shrink the hash table */ \
            h->keys = (khkey_t*)krealloc((void *)h->keys, new_n_buckets * sizeof(khkey_t)); \
        kfree(h->used); /* free the working space */ \
        h->used = new_used, h->bits = new_bits; \
        return 0; \
    }

#define __KHASHL_IMPL_PUT(SCOPE, HType, prefix, khkey_t, __hash_fn, __hash_eq) \
    SCOPE khint_t prefix##_putp(HType *h, const khkey_t *key, int *absent) { \
        khint_t n_buckets, i, last, mask; \
        n_buckets = h->keys? 1U<<h->bits : 0U; \
        *absent = -1; \
        if (h->count >= (n_buckets>>1) + (n_buckets>>2)) { /* rehashing */ \
            if (prefix##_resize(h, n_buckets + 1U) < 0) \
                return n_buckets; \
            n_buckets = 1U<<h->bits; \
        } /* TODO: to implement automatically shrinking; resize() already support shrinking */ \
        mask = n_buckets - 1; \
        i = last = __kh_h2b(__hash_fn(*key), h->bits); \
        while (__kh_used(h->used, i) && !__hash_eq(h->keys[i], *key)) { \
            i = (i + 1U) & mask; \
            if (i == last) break; \
        } \
        if (!__kh_used(h->used, i)) { /* not present at all */ \
            h->keys[i] = *key; \
            __kh_set_used(h->used, i); \
            ++h->count; \
            *absent = 1; \
        } else *absent = 0; /* Don't touch h->keys[i] if present */ \
        return i; \
    } \
    SCOPE khint_t prefix##_put(HType *h, khkey_t key, int *absent) { return prefix##_putp(h, &key, absent); }

#define __KHASHL_IMPL_DEL(SCOPE, HType, prefix, khkey_t, __hash_fn) \
    SCOPE int prefix##_del(HType *h, khint_t i) { \
        khint_t j = i, k, mask, n_buckets; \
        if (h->keys == 0) return 0; \
        n_buckets = 1U<<h->bits; \
        mask = n_buckets - 1U; \
        while (1) { \
            j = (j + 1U) & mask; \
            if (j == i || !__kh_used(h->used, j)) break; /* j==i only when the table is completely full */ \
            k = __kh_h2b(__hash_fn(h->keys[j]), h->bits); \
            if ((j > i && (k <= i || k > j)) || (j < i && (k <= i && k > j))) \
                h->keys[i] = h->keys[j], i = j; \
        } \
        __kh_set_unused(h->used, i); \
        --h->count; \
        return 1; \
    }

#define KHASHL_DECLARE(HType, prefix, khkey_t) \
    __KHASHL_TYPE(HType, khkey_t) \
    __KHASHL_PROTOTYPES(HType, prefix, khkey_t)

#define KHASHL_INIT(SCOPE, HType, prefix, khkey_t, __hash_fn, __hash_eq) \
    __KHASHL_TYPE(HType, khkey_t) \
    __KHASHL_IMPL_BASIC(SCOPE, HType, prefix) \
    __KHASHL_IMPL_GET(SCOPE, HType, prefix, khkey_t, __hash_fn, __hash_eq) \
    __KHASHL_IMPL_RESIZE(SCOPE, HType, prefix, khkey_t, __hash_fn, __hash_eq) \
    __KHASHL_IMPL_PUT(SCOPE, HType, prefix, khkey_t, __hash_fn, __hash_eq) \
    __KHASHL_IMPL_DEL(SCOPE, HType, prefix, khkey_t, __hash_fn)

/*****************************
 * More convenient interface *
 *****************************/

#define __kh_packed __attribute__ ((__packed__))
#define __kh_cached_hash(x) ((x).hash)

#define KHASHL_SET_INIT(SCOPE, HType, prefix, khkey_t, __hash_fn, __hash_eq) \
    typedef struct { khkey_t key; } __kh_packed HType##_s_bucket_t; \
    static kh_inline khint_t prefix##_s_hash(HType##_s_bucket_t x) { return __hash_fn(x.key); } \
    static kh_inline int prefix##_s_eq(HType##_s_bucket_t x, HType##_s_bucket_t y) { return __hash_eq(x.key, y.key); } \
    KHASHL_INIT(KH_LOCAL, HType, prefix##_s, HType##_s_bucket_t, prefix##_s_hash, prefix##_s_eq) \
    SCOPE HType *prefix##_init(void) { return prefix##_s_init(); } \
    SCOPE void prefix##_destroy(HType *h) { prefix##_s_destroy(h); } \
    SCOPE void prefix##_resize(HType *h, khint_t new_n_buckets) { prefix##_s_resize(h, new_n_buckets); } \
    SCOPE khint_t prefix##_get(const HType *h, khkey_t key) { HType##_s_bucket_t t; t.key = key; return prefix##_s_getp(h, &t); } \
    SCOPE int prefix##_del(HType *h, khint_t k) { return prefix##_s_del(h, k); } \
    SCOPE khint_t prefix##_put(HType *h, khkey_t key, int *absent) { HType##_s_bucket_t t; t.key = key; return prefix##_s_putp(h, &t, absent); }

#define KHASHL_MAP_INIT(SCOPE, HType, prefix, khkey_t, kh_val_t, __hash_fn, __hash_eq) \
    typedef struct { khkey_t key; kh_val_t val; } __kh_packed HType##_m_bucket_t; \
    static kh_inline khint_t prefix##_m_hash(HType##_m_bucket_t x) { return __hash_fn(x.key); } \
    static kh_inline int prefix##_m_eq(HType##_m_bucket_t x, HType##_m_bucket_t y) { return __hash_eq(x.key, y.key); } \
    KHASHL_INIT(KH_LOCAL, HType, prefix##_m, HType##_m_bucket_t, prefix##_m_hash, prefix##_m_eq) \
    SCOPE HType *prefix##_init(void) { return prefix##_m_init(); } \
    SCOPE void prefix##_destroy(HType *h) { prefix##_m_destroy(h); } \
    SCOPE khint_t prefix##_get(const HType *h, khkey_t key) { HType##_m_bucket_t t; t.key = key; return prefix##_m_getp(h, &t); } \
    SCOPE int prefix##_del(HType *h, khint_t k) { return prefix##_m_del(h, k); } \
    SCOPE khint_t prefix##_put(HType *h, khkey_t key, int *absent) { HType##_m_bucket_t t; t.key = key; return prefix##_m_putp(h, &t, absent); }

#define KHASHL_CSET_INIT(SCOPE, HType, prefix, khkey_t, __hash_fn, __hash_eq) \
    typedef struct { khkey_t key; khint_t hash; } __kh_packed HType##_cs_bucket_t; \
    static kh_inline int prefix##_cs_eq(HType##_cs_bucket_t x, HType##_cs_bucket_t y) { return x.hash == y.hash && __hash_eq(x.key, y.key); } \
    KHASHL_INIT(KH_LOCAL, HType, prefix##_cs, HType##_cs_bucket_t, __kh_cached_hash, prefix##_cs_eq) \
    SCOPE HType *prefix##_init(void) { return prefix##_cs_init(); } \
    SCOPE void prefix##_destroy(HType *h) { prefix##_cs_destroy(h); } \
    SCOPE khint_t prefix##_get(const HType *h, khkey_t key) { HType##_cs_bucket_t t; t.key = key; t.hash = __hash_fn(key); return prefix##_cs_getp(h, &t); } \
    SCOPE int prefix##_del(HType *h, khint_t k) { return prefix##_cs_del(h, k); } \
    SCOPE khint_t prefix##_put(HType *h, khkey_t key, int *absent) { HType##_cs_bucket_t t; t.key = key, t.hash = __hash_fn(key); return prefix##_cs_putp(h, &t, absent); }

#define KHASHL_CMAP_INIT(SCOPE, HType, prefix, khkey_t, kh_val_t, __hash_fn, __hash_eq) \
    typedef struct { khkey_t key; kh_val_t val; khint_t hash; } __kh_packed HType##_cm_bucket_t; \
    static kh_inline int prefix##_cm_eq(HType##_cm_bucket_t x, HType##_cm_bucket_t y) { return x.hash == y.hash && __hash_eq(x.key, y.key); } \
    KHASHL_INIT(KH_LOCAL, HType, prefix##_cm, HType##_cm_bucket_t, __kh_cached_hash, prefix##_cm_eq) \
    SCOPE HType *prefix##_init(void) { return prefix##_cm_init(); } \
    SCOPE void prefix##_destroy(HType *h) { prefix##_cm_destroy(h); } \
    SCOPE khint_t prefix##_get(const HType *h, khkey_t key) { HType##_cm_bucket_t t; t.key = key; t.hash = __hash_fn(key); return prefix##_cm_getp(h, &t); } \
    SCOPE int prefix##_del(HType *h, khint_t k) { return prefix##_cm_del(h, k); } \
    SCOPE khint_t prefix##_put(HType *h, khkey_t key, int *absent) { HType##_cm_bucket_t t; t.key = key, t.hash = __hash_fn(key); return prefix##_cm_putp(h, &t, absent); }

/**************************
 * Public macro functions *
 **************************/

#define kh_bucket(h, x) ((h)->keys[x])
#define kh_size(h) ((h)->count)
#define kh_capacity(h) ((h)->keys? 1U<<(h)->bits : 0U)
#define kh_end(h) kh_capacity(h)

#define kh_key(h, x) ((h)->keys[x].key)
#define kh_val(h, x) ((h)->keys[x].val)
#define kh_exist(h, x) __kh_used((h)->used, (x))

/**************************************
 * Common hash and equality functions *
 **************************************/

#define kh_eq_generic(a, b) ((a) == (b))
#define kh_eq_str(a, b) (strcmp((a), (b)) == 0)
#define kh_hash_dummy(x) ((khint_t)(x))

static kh_inline khint_t kh_hash_uint32(khint_t key) {
    key += ~(key << 15);
    key ^=  (key >> 10);
    key +=  (key << 3);
    key ^=  (key >> 6);
    key += ~(key << 11);
    key ^=  (key >> 16);
    return key;
}

static kh_inline khint_t kh_hash_uint64(khint64_t key) {
    key = ~key + (key << 21);
    key = key ^ key >> 24;
    key = (key + (key << 3)) + (key << 8);
    key = key ^ key >> 14;
    key = (key + (key << 2)) + (key << 4);
    key = key ^ key >> 28;
    key = key + (key << 31);
    return (khint_t)key;
}

static kh_inline khint_t kh_hash_str(const char *s) {
    khint_t h = (khint_t)*s;
    if (h) for (++s ; *s; ++s) h = (h << 5) - h + (khint_t)*s;
    return h;
}

#endif /* __AC_KHASHL_H */

/*******************
 * end of khashl.h *
 *******************/

typedef struct {
    int size;
    uint8_t *block;
    int64_t end_offset;
} cache_t;
KHASHL_MAP_INIT(, kh_t, kh, uint64_t, cache_t, kh_hash_uint64, kh_eq_generic)

typedef uint8_t bgzf_byte_t;
static const int DEFAULT_BLOCK_SIZE = 64 * 1024;
static const int MAX_BLOCK_SIZE = 64 * 1024;
static const int BLOCK_HEADER_LENGTH = 18;
static const int BLOCK_FOOTER_LENGTH = 8;

static const int GZIP_ID1 = 31;
static const int GZIP_ID2 = 139;
static const int CM_DEFLATE = 8;
static const int FLG_FEXTRA = 4;
static const int OS_UNKNOWN = 255;
static const int BGZF_ID1 = 66; // 'B'
static const int BGZF_ID2 = 67; // 'C'
static const int BGZF_LEN = 2;
static const int BGZF_XLEN = 6; // BGZF_LEN+4
static const int GZIP_WINDOW_BITS = -15; // no zlib header
static const int Z_DEFAULT_MEM_LEVEL = 8;

/*********************
 * from bam_endian.c *
 *********************/
static inline int bam_is_big_endian()
{
    long one= 1;
    return !(*((char *)(&one)));
}
static inline uint16_t bam_swap_endian_2(uint16_t v)
{
    return (uint16_t)(((v & 0x00FF00FFU) << 8) | ((v & 0xFF00FF00U) >> 8));
}
static inline void *bam_swap_endian_2p(void *x)
{
    *(uint16_t*)x = bam_swap_endian_2(*(uint16_t*)x);
    return x;
}
static inline uint32_t bam_swap_endian_4(uint32_t v)
{
    v = ((v & 0x0000FFFFU) << 16) | (v >> 16);
    return ((v & 0x00FF00FFU) << 8) | ((v & 0xFF00FF00U) >> 8);
}
static inline void *bam_swap_endian_4p(void *x)
{
    *(uint32_t*)x = bam_swap_endian_4(*(uint32_t*)x);
    return x;
}
static inline uint64_t bam_swap_endian_8(uint64_t v)
{
    v = ((v & 0x00000000FFFFFFFFLLU) << 32) | (v >> 32);
    v = ((v & 0x0000FFFF0000FFFFLLU) << 16) | ((v & 0xFFFF0000FFFF0000LLU) >> 16);
    return ((v & 0x00FF00FF00FF00FFLLU) << 8) | ((v & 0xFF00FF00FF00FF00LLU) >> 8);
}
static inline void *bam_swap_endian_8p(void *x)
{
    *(uint64_t*)x = bam_swap_endian_8(*(uint64_t*)x);
    return x;
}

/**************
 * from bam.c *
 **************/
int bam_is_be;

bam_hdr_t *bam_hdr_init()
{
    bam_is_be = bam_is_big_endian();
    return (bam_hdr_t*)calloc(1, sizeof(bam_hdr_t));
}

void bam_hdr_destroy(bam_hdr_t *header)
{
    int32_t i;
    if (header == 0) return;
    if (header->target_name) {
        for (i = 0; i < header->n_targets; ++i)
            if (header->target_name[i]) free(header->target_name[i]);
        if (header->target_len) free(header->target_len);
        free(header->target_name);
    }
    if (header->text) free(header->text);
    free(header);
}

bam_hdr_t *bam_hdr_read(bamFile fp)
{
    bam_hdr_t *header;
    char buf[4];
    int magic_len;
    int32_t i = 1, name_len;
    // read "BAM1"
    magic_len = bam_read(fp, buf, 4);
    if (magic_len != 4 || strncmp(buf, "BAM\001", 4) != 0) {
        fprintf(stderr, "[bam_hdr_read] invalid BAM binary header (this is not a BAM file).\n");
        return NULL;
    }
    header = bam_hdr_init();
    // read plain text and the number of reference sequences
    if (bam_read(fp, &header->l_text, 4) != 4) goto fail; 
    if (bam_is_be) bam_swap_endian_4p(&header->l_text);
    header->text = (char*)calloc(header->l_text + 1, 1);
    if (bam_read(fp, header->text, header->l_text) != header->l_text) goto fail;
    if (bam_read(fp, &header->n_targets, 4) != 4) goto fail;
    if (bam_is_be) bam_swap_endian_4p(&header->n_targets);
    // read reference sequence names and lengths
    header->target_name = (char**)calloc(header->n_targets, sizeof(char*));
    header->target_len = (uint32_t*)calloc(header->n_targets, 4);
    for (i = 0; i != header->n_targets; ++i) {
        if (bam_read(fp, &name_len, 4) != 4) goto fail;
        if (bam_is_be) bam_swap_endian_4p(&name_len);
        header->target_name[i] = (char*)calloc(name_len, 1);
        if (bam_read(fp, header->target_name[i], name_len) != name_len) {
            goto fail;
        }
        if (bam_read(fp, &header->target_len[i], 4) != 4) goto fail;
        if (bam_is_be) bam_swap_endian_4p(&header->target_len[i]);
    }
    return header;
 fail:
    bam_hdr_destroy(header);
    return NULL;
}

static void swap_endian_data(const bam1_core_t *c, int l_data, uint8_t *data)
{
    uint8_t *s;
    uint32_t i, *cigar = (uint32_t*)(data + c->l_qname);
    s = data + c->n_cigar*4 + c->l_qname + c->l_qseq + (c->l_qseq + 1)/2;
    for (i = 0; i < c->n_cigar; ++i) bam_swap_endian_4p(&cigar[i]);
    while (s < data + l_data) {
        uint8_t type;
        s += 2; // skip key
        type = toupper(*s); ++s; // skip type
        if (type == 'C' || type == 'A') ++s;
        else if (type == 'S') { bam_swap_endian_2p(s); s += 2; }
        else if (type == 'I' || type == 'F') { bam_swap_endian_4p(s); s += 4; }
        else if (type == 'D') { bam_swap_endian_8p(s); s += 8; }
        else if (type == 'Z' || type == 'H') { while (*s) ++s; ++s; }
    }
}

int bam_hdr_write(bamFile fp, const bam_hdr_t *header)
{
    char buf[4];
    int32_t i, name_len, x;
    // write "BAM1"
    strncpy(buf, "BAM\001", 4);
    bam_write(fp, buf, 4);
    // write plain text and the number of reference sequences
    if (bam_is_be) {
        x = bam_swap_endian_4(header->l_text);
        bam_write(fp, &x, 4);
        if (header->l_text) bam_write(fp, header->text, header->l_text);
        x = bam_swap_endian_4(header->n_targets);
        bam_write(fp, &x, 4);
    } else {
        bam_write(fp, &header->l_text, 4);
        if (header->l_text) bam_write(fp, header->text, header->l_text);
        bam_write(fp, &header->n_targets, 4);
    }
    // write sequence names and lengths
    for (i = 0; i != header->n_targets; ++i) {
        char *p = header->target_name[i];
        name_len = strlen(p) + 1;
        if (bam_is_be) {
            x = bam_swap_endian_4(name_len);
            bam_write(fp, &x, 4);
        } else bam_write(fp, &name_len, 4);
        bam_write(fp, p, name_len);
        if (bam_is_be) {
            x = bam_swap_endian_4(header->target_len[i]);
            bam_write(fp, &x, 4);
        } else bam_write(fp, &header->target_len[i], 4);
    }
    bgzf_flush(fp);
    return 0;
}

int bam_read1(bamFile fp, bam1_t *b)
{
    bam1_core_t *c = &b->core;
    int32_t block_len, ret, i;
    uint32_t x[8];

    if ((ret = bam_read(fp, &block_len, 4)) != 4) {
        if (ret == 0) return -1; // normal end-of-file
        else return -2; // truncated
    }
    if (bam_read(fp, x, sizeof(bam1_core_t)) != sizeof(bam1_core_t)) return -3;
    if (bam_is_be) {
        bam_swap_endian_4p(&block_len);
        for (i = 0; i < 8; ++i) bam_swap_endian_4p(x + i);
    }
    c->tid = x[0]; c->pos = x[1];
    c->bin = x[2]>>16; c->qual = x[2]>>8&0xff; c->l_qname = x[2]&0xff;
    c->flag = x[3]>>16; c->n_cigar = x[3]&0xffff;
    c->l_qseq = x[4];
    c->mtid = x[5]; c->mpos = x[6]; c->isize = x[7];
    b->l_data = block_len - sizeof(bam1_core_t);
    if (b->m_data < b->l_data) {
        b->m_data = b->l_data;
        kroundup32(b->m_data);
        b->data = (uint8_t*)realloc(b->data, b->m_data);
    }
    if (bam_read(fp, b->data, b->l_data) != b->l_data) return -4;
    b->l_aux = b->l_data - c->n_cigar * 4 - c->l_qname - c->l_qseq - (c->l_qseq+1)/2;
    if (bam_is_be) swap_endian_data(c, b->l_data, b->data);
    return 4 + block_len;
}

static BGZF *bgzf_read_init()
{
    BGZF *fp;
    fp = calloc(1, sizeof(BGZF));
    fp->uncompressed_block_size = MAX_BLOCK_SIZE;
    fp->uncompressed_block = malloc(MAX_BLOCK_SIZE);
    fp->compressed_block_size = MAX_BLOCK_SIZE;
    fp->compressed_block = malloc(MAX_BLOCK_SIZE);
    fp->cache_size = 0;
    fp->cache = kh_init();
    return fp;
}

static BGZF *open_read(int fd)
{
    FILE* file = fdopen(fd, "r");
    BGZF* fp;
    if (file == 0) return 0;
    fp = bgzf_read_init();
    fp->file_descriptor = fd;
    fp->open_mode = 'r';
    fp->file = file;
    return fp;
}

static BGZF *open_write(int fd, int compress_level) // compress_level==-1 for the default level
{
    FILE* file = fdopen(fd, "w");
    BGZF* fp;
    if (file == 0) return 0;
    fp = malloc(sizeof(BGZF));
    fp->file_descriptor = fd;
    fp->open_mode = 'w';
    fp->owned_file = 0;
    fp->compress_level = compress_level < 0? Z_DEFAULT_COMPRESSION : compress_level; // Z_DEFAULT_COMPRESSION==-1
    if (fp->compress_level > 9) fp->compress_level = Z_DEFAULT_COMPRESSION;
    fp->file = file;
    fp->uncompressed_block_size = DEFAULT_BLOCK_SIZE;
    fp->uncompressed_block = NULL;
    fp->compressed_block_size = MAX_BLOCK_SIZE;
    fp->compressed_block = malloc(MAX_BLOCK_SIZE);
    fp->block_address = 0;
    fp->block_offset = 0;
    fp->block_length = 0;
    fp->error = NULL;
    return fp;
}

BGZF *bgzf_open(const char* __restrict path, const char* __restrict mode)
{
    BGZF* fp = NULL;
    if (strchr(mode, 'r') || strchr(mode, 'R')) { /* The reading mode is preferred. */
        int fd, oflag = O_RDONLY;
        fd = open(path, oflag);
        if (fd == -1) return 0;
        fp = open_read(fd);
    } else if (strchr(mode, 'w') || strchr(mode, 'W')) {
        int fd, compress_level = -1, oflag = O_WRONLY | O_CREAT | O_TRUNC;
        fd = open(path, oflag, 0666);
        if (fd == -1) return 0;
        { // set compress_level
            int i;
            for (i = 0; mode[i]; ++i)
                if (mode[i] >= '0' && mode[i] <= '9') break;
            if (mode[i]) compress_level = (int)mode[i] - '0';
            if (strchr(mode, 'u')) compress_level = 0;
        }
        fp = open_write(fd, compress_level);
    }
    if (fp != NULL) fp->owned_file = 1;
    return fp;
}

BGZF *bgzf_fdopen(int fd, const char * __restrict mode)
{
    if (fd == -1) return 0;
    if (mode[0] == 'r' || mode[0] == 'R') {
        return open_read(fd);
    } else if (mode[0] == 'w' || mode[0] == 'W') {
        int i, compress_level = -1;
        for (i = 0; mode[i]; ++i)
            if (mode[i] >= '0' && mode[i] <= '9') break;
        if (mode[i]) compress_level = (int)mode[i] - '0';
        if (strchr(mode, 'u')) compress_level = 0;
        return open_write(fd, compress_level);
    } else {
        return NULL;
    }
}

static inline void packInt16(uint8_t* buffer, uint16_t value)
{
    buffer[0] = value;
    buffer[1] = value >> 8;
}

static inline int unpackInt16(const uint8_t* buffer)
{
    return (buffer[0] | (buffer[1] << 8));
}

static inline void packInt32(uint8_t* buffer, uint32_t value)
{
    buffer[0] = value;
    buffer[1] = value >> 8;
    buffer[2] = value >> 16;
    buffer[3] = value >> 24;
}

static inline int bgzf_min(int x, int y)
{
    return (x < y) ? x : y;
}

static void report_error(BGZF* fp, const char* message) {
    fp->error = message;
}

static int deflate_block(BGZF* fp, int block_length)
{
    // Deflate the block in fp->uncompressed_block into fp->compressed_block.
    // Also adds an extra field that stores the compressed block length.

    bgzf_byte_t* buffer = fp->compressed_block;
    int buffer_size = fp->compressed_block_size;

    // Init gzip header
    buffer[0] = GZIP_ID1;
    buffer[1] = GZIP_ID2;
    buffer[2] = CM_DEFLATE;
    buffer[3] = FLG_FEXTRA;
    buffer[4] = 0; // mtime
    buffer[5] = 0;
    buffer[6] = 0;
    buffer[7] = 0;
    buffer[8] = 0;
    buffer[9] = OS_UNKNOWN;
    buffer[10] = BGZF_XLEN;
    buffer[11] = 0;
    buffer[12] = BGZF_ID1;
    buffer[13] = BGZF_ID2;
    buffer[14] = BGZF_LEN;
    buffer[15] = 0;
    buffer[16] = 0; // placeholder for block length
    buffer[17] = 0;

    // loop to retry for blocks that do not compress enough
    int input_length = block_length;
    int compressed_length = 0;
    while (1) {
        z_stream zs;
        zs.zalloc = NULL;
        zs.zfree = NULL;
        zs.next_in = fp->uncompressed_block;
        zs.avail_in = input_length;
        zs.next_out = (void*)&buffer[BLOCK_HEADER_LENGTH];
        zs.avail_out = buffer_size - BLOCK_HEADER_LENGTH - BLOCK_FOOTER_LENGTH;

        int status = deflateInit2(&zs, fp->compress_level, Z_DEFLATED,
                                  GZIP_WINDOW_BITS, Z_DEFAULT_MEM_LEVEL, Z_DEFAULT_STRATEGY);
        if (status != Z_OK) {
            report_error(fp, "deflate init failed");
            return -1;
        }
        status = deflate(&zs, Z_FINISH);
        if (status != Z_STREAM_END) {
            deflateEnd(&zs);
            if (status == Z_OK) {
                // Not enough space in buffer.
                // Can happen in the rare case the input doesn't compress enough.
                // Reduce the amount of input until it fits.
                input_length -= 1024;
                if (input_length <= 0) {
                    // should never happen
                    report_error(fp, "input reduction failed");
                    return -1;
                }
                continue;
            }
            report_error(fp, "deflate failed");
            return -1;
        }
        status = deflateEnd(&zs);
        if (status != Z_OK) {
            report_error(fp, "deflate end failed");
            return -1;
        }
        compressed_length = zs.total_out;
        compressed_length += BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH;
        if (compressed_length > MAX_BLOCK_SIZE) {
            // should never happen
            report_error(fp, "deflate overflow");
            return -1;
        }
        break;
    }

    packInt16((uint8_t*)&buffer[16], compressed_length-1);
    uint32_t crc = crc32(0L, NULL, 0L);
    crc = crc32(crc, fp->uncompressed_block, input_length);
    packInt32((uint8_t*)&buffer[compressed_length-8], crc);
    packInt32((uint8_t*)&buffer[compressed_length-4], input_length);

    int remaining = block_length - input_length;
    if (remaining > 0) {
        if (remaining > input_length) {
            // should never happen (check so we can use memcpy)
            report_error(fp, "remainder too large");
            return -1;
        }
        memcpy(fp->uncompressed_block,
               fp->uncompressed_block + input_length,
               remaining);
    }
    fp->block_offset = remaining;
    return compressed_length;
}

static int inflate_block(BGZF* fp, int block_length)
{
    // Inflate the block in fp->compressed_block into fp->uncompressed_block

    z_stream zs;
    int status;
    zs.zalloc = NULL;
    zs.zfree = NULL;
    zs.next_in = fp->compressed_block + 18;
    zs.avail_in = block_length - 16;
    zs.next_out = fp->uncompressed_block;
    zs.avail_out = fp->uncompressed_block_size;

    status = inflateInit2(&zs, GZIP_WINDOW_BITS);
    if (status != Z_OK) {
        report_error(fp, "inflate init failed");
        return -1;
    }
    status = inflate(&zs, Z_FINISH);
    if (status != Z_STREAM_END) {
        inflateEnd(&zs);
        report_error(fp, "inflate failed");
        return -1;
    }
    status = inflateEnd(&zs);
    if (status != Z_OK) {
        report_error(fp, "inflate failed");
        return -1;
    }
    return zs.total_out;
}

static int check_header(const bgzf_byte_t* header)
{
    return (header[0] == GZIP_ID1 &&
            header[1] == (bgzf_byte_t) GZIP_ID2 &&
            header[2] == Z_DEFLATED &&
            (header[3] & FLG_FEXTRA) != 0 &&
            unpackInt16((uint8_t*)&header[10]) == BGZF_XLEN &&
            header[12] == BGZF_ID1 &&
            header[13] == BGZF_ID2 &&
            unpackInt16((uint8_t*)&header[14]) == BGZF_LEN);
}

static void free_cache(BGZF *fp)
{
    khint_t k;
    kh_t *h = (kh_t *)fp->cache;
    if (fp->open_mode != 'r') return;
    for (k = 0; k < kh_end(h); ++k)
        if (kh_exist(h, k)) free(kh_val(h, k).block);
    kh_destroy(h);
}

static int load_block_from_cache(BGZF *fp, int64_t block_address)
{
    khint_t k;
    cache_t *p;
    kh_t *h = (kh_t *)fp->cache;
    k = kh_get(h, block_address);
    if (k == kh_end(h)) return 0;
    p = &kh_val(h, k);
    if (fp->block_length != 0) fp->block_offset = 0;
    fp->block_address = block_address;
    fp->block_length = p->size;
    memcpy(fp->uncompressed_block, p->block, MAX_BLOCK_SIZE);
    fseeko(fp->file, p->end_offset, SEEK_SET);
    return p->size;
}

static void cache_block(BGZF *fp, int size)
{
    int ret;
    khint_t k;
    cache_t *p;
    kh_t *h = (kh_t *)fp->cache;
    if (MAX_BLOCK_SIZE >= fp->cache_size) return;
    if ((kh_size(h) + 1) * MAX_BLOCK_SIZE > fp->cache_size) {
        /* A better way would be to remove the oldest block in the
         * cache, but here we remove a random one for simplicity. This
         * should not have a big impact on performance. */
        for (k = 0; k < kh_end(h); ++k)
            if (kh_exist(h, k)) break;
        if (k < kh_end(h)) {
            free(kh_val(h, k).block);
            kh_del(h, k);
        }
    }
    k = kh_put(h, fp->block_address, &ret);
    if (ret == 0) return; // if this happens, a bug!
    p = &kh_val(h, k);
    p->size = fp->block_length;
    p->end_offset = fp->block_address + size;
    p->block = malloc(MAX_BLOCK_SIZE);
    memcpy(kh_val(h, k).block, fp->uncompressed_block, MAX_BLOCK_SIZE);
}

int bgzf_read_block(BGZF* fp)
{
    bgzf_byte_t header[BLOCK_HEADER_LENGTH];
    int count, size = 0, block_length, remaining;
    int64_t block_address = ftello(fp->file);
    if (load_block_from_cache(fp, block_address)) return 0;
    count = fread(header, 1, sizeof(header), fp->file);
    if (count == 0) {
        fp->block_length = 0;
        return 0;
    }
    size = count;
    if (count != sizeof(header)) {
        report_error(fp, "read failed");
        return -1;
    }
    if (!check_header(header)) {
        report_error(fp, "invalid block header");
        return -1;
    }
    block_length = unpackInt16((uint8_t*)&header[16]) + 1;
    bgzf_byte_t* compressed_block = (bgzf_byte_t*) fp->compressed_block;
    memcpy(compressed_block, header, BLOCK_HEADER_LENGTH);
    remaining = block_length - BLOCK_HEADER_LENGTH;
    count = fread(&compressed_block[BLOCK_HEADER_LENGTH], 1, remaining, fp->file);
    if (count != remaining) {
        report_error(fp, "read failed");
        return -1;
    }
    size += count;
    count = inflate_block(fp, block_length);
    if (count < 0) return -1;
    if (fp->block_length != 0) {
        // Do not reset offset if this read follows a seek.
        fp->block_offset = 0;
    }
    fp->block_address = block_address;
    fp->block_length = count;
    cache_block(fp, size);
    return 0;
}

int bgzf_read(BGZF* fp, void* data, int length)
{
    if (length <= 0) {
        return 0;
    }
    if (fp->open_mode != 'r') {
        report_error(fp, "file not open for reading");
        return -1;
    }

    int bytes_read = 0;
    bgzf_byte_t* output = data;
    while (bytes_read < length) {
        int copy_length, available = fp->block_length - fp->block_offset;
        bgzf_byte_t *buffer;
        if (available <= 0) {
            if (bgzf_read_block(fp) != 0) {
                return -1;
            }
            available = fp->block_length - fp->block_offset;
            if (available <= 0) {
                break;
            }
        }
        copy_length = bgzf_min(length-bytes_read, available);
        buffer = fp->uncompressed_block;
        memcpy(output, buffer + fp->block_offset, copy_length);
        fp->block_offset += copy_length;
        output += copy_length;
        bytes_read += copy_length;
    }
    if (fp->block_offset == fp->block_length) {
        fp->block_address = ftello(fp->file);
        fp->block_offset = 0;
        fp->block_length = 0;
    }
    return bytes_read;
}

int bgzf_flush(BGZF* fp)
{
    while (fp->block_offset > 0) {
        int count, block_length;
        block_length = deflate_block(fp, fp->block_offset);
        if (block_length < 0) return -1;
        count = fwrite(fp->compressed_block, 1, block_length, fp->file);
        if (count != block_length) {
            report_error(fp, "write failed");
            return -1;
        }
        fp->block_address += block_length;
    }
    return 0;
}

int bgzf_flush_try(BGZF *fp, int size)
{
    if (fp->block_offset + size > fp->uncompressed_block_size)
        return bgzf_flush(fp);
    return -1;
}

int bgzf_write(BGZF* fp, const void* data, int length)
{
    const bgzf_byte_t *input = data;
    int block_length, bytes_written;
    if (fp->open_mode != 'w') {
        report_error(fp, "file not open for writing");
        return -1;
    }

    if (fp->uncompressed_block == NULL)
        fp->uncompressed_block = malloc(fp->uncompressed_block_size);

    input = data;
    block_length = fp->uncompressed_block_size;
    bytes_written = 0;
    while (bytes_written < length) {
        int copy_length = bgzf_min(block_length - fp->block_offset, length - bytes_written);
        bgzf_byte_t* buffer = fp->uncompressed_block;
        memcpy(buffer + fp->block_offset, input, copy_length);
        fp->block_offset += copy_length;
        input += copy_length;
        bytes_written += copy_length;
        if (fp->block_offset == block_length) {
            if (bgzf_flush(fp) != 0) {
                break;
            }
        }
    }
    return bytes_written;
}

int bgzf_close(BGZF* fp)
{
    if (fp->open_mode == 'w') {
        if (bgzf_flush(fp) != 0) return -1;
        { // add an empty block
            int count, block_length = deflate_block(fp, 0);
            count = fwrite(fp->compressed_block, 1, block_length, fp->file);
        }
        if (fflush(fp->file) != 0) {
            report_error(fp, "flush failed");
            return -1;
        }
    }
    if (fp->owned_file)
        if (fclose(fp->file) != 0) return -1;
    free(fp->uncompressed_block);
    free(fp->compressed_block);
    free_cache(fp);
    free(fp);
    return 0;
}


inline int bam_write1_core(bamFile fp, const bam1_core_t *c, int data_len, uint8_t *data)
{
    uint32_t x[8], block_len = data_len + BAM_CORE_SIZE, y;
    int i;
    assert(BAM_CORE_SIZE == 32);
    x[0] = c->tid;
    x[1] = c->pos;
    x[2] = (uint32_t)c->bin<<16 | c->qual<<8 | c->l_qname;
    x[3] = (uint32_t)c->flag<<16 | c->n_cigar;
    x[4] = c->l_qseq;
    x[5] = c->mtid;
    x[6] = c->mpos;
    x[7] = c->isize;
    bgzf_flush_try(fp, 4 + block_len);
    if (bam_is_be) {
        for (i = 0; i < 8; ++i) bam_swap_endian_4p(x + i);
        y = block_len;
        bam_write(fp, bam_swap_endian_4p(&y), 4);
        swap_endian_data(c, data_len, data);
    } else bam_write(fp, &block_len, 4);
    bam_write(fp, x, BAM_CORE_SIZE);
    bam_write(fp, data, data_len);
    if (bam_is_be) swap_endian_data(c, data_len, data);
    return 4 + block_len;
}

int bam_write1(bamFile fp, const bam1_t *b)
{
    return bam_write1_core(fp, &b->core, b->l_data, b->data);
}

//***** bam aux ******//
int get_aux_type2size(uint8_t t)
{
    if (t == 'C' || t == 'c' || t == 'A') return 1;
    else if (t == 's' || t == 'S') return 2;
    else return 4; // i I f F    
}

#define _skip_aux(s) do { \
    int c = toupper(*(s)); \
    ++(s); \
    if (c == 'H' || c == 'Z') { while (*(s)) ++s; ++s;} \
    else if (c == 'B') (s) += 5 + get_aux_type2size(*((s)+4)) * (*(int32_t*)((s)+1)); \
    else s += get_aux_type2size(c); \
} while(0)

uint8_t *bam_aux_get(const bam1_t *b, const char tag[2])
{
    uint8_t *s;
    uint16_t y = (uint16_t) tag[0]<<8 | tag[1];
    s = bam1_aux(b);
    while (s < b->data + b->l_data) {
        uint16_t x = (uint16_t) s[0]<<8 | s[1];
        s += 2;
        if (x == y) return s; 
        else _skip_aux(s);
    }
    return 0;
}
