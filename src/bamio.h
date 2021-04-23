#ifndef BAMIO_H_
#define BAMIO_H_

#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include <zlib.h>

#define bam_open(fn, mode)      bgzf_open(fn, mode)
#define bam_dopen(fd, mode)     bgzf_fdopen(fd, mode)
#define bam_close(fp)           bgzf_close(fp)
#define bam_read(fp, buf, size) bgzf_read(fp, buf, size)

#define BAM_CORE_SIZE   sizeof(bam1_core_t)
#define BAM_FPAIRED        1
#define BAM_FPROPER_PAIR   2
#define BAM_FUNMAP         4
#define BAM_FMUNMAP        8
#define BAM_FREVERSE      16
#define BAM_FMREVERSE     32
#define BAM_FREAD1        64
#define BAM_FREAD2       128
#define BAM_FSECONDARY   256
#define BAM_FQCFAIL      512
#define BAM_FDUP        1024
#define BAM_FSUPPLEMENTARY 2048

#define BAM_CIGAR_STR   "MIDNSHP=XB"
#define BAM_CIGAR_SHIFT 4
#define BAM_CIGAR_MASK  ((1 << BAM_CIGAR_SHIFT) - 1)
#define BAM_CIGAR_TYPE  0x3C1A7

#define bam_cigar_op(c) ((c)&BAM_CIGAR_MASK)
#define bam_cigar_oplen(c) ((c)>>BAM_CIGAR_SHIFT)
#define bam_cigar_opchr(c) (BAM_CIGAR_STR[bam_cigar_op(c)])
#define bam_cigar_gen(l, o) ((l)<<BAM_CIGAR_SHIFT|(o))
#define bam_cigar_type(o) (BAM_CIGAR_TYPE>>((o)<<1)&3)

#define BAM_CMATCH      0
#define BAM_CINS        1
#define BAM_CDEL        2
#define BAM_CREF_SKIP   3
#define BAM_CSOFT_CLIP  4
#define BAM_CHARD_CLIP  5
#define BAM_CPAD        6

typedef struct {
    int file_descriptor;
    char open_mode;  // 'r' or 'w'
    int16_t owned_file, compress_level;
    FILE* file;
    int uncompressed_block_size;
    int compressed_block_size;
    void* uncompressed_block;
    void* compressed_block;
    int64_t block_address;
    int block_length;
    int block_offset;
    int cache_size;
    const char* error;
    void *cache; // a pointer to a hash table
} BGZF;

typedef BGZF *bamFile;
#define bam_write(fp, buf, size) bgzf_write(fp, buf, size)

typedef struct {
    int32_t n_targets;
    char **target_name;
    uint32_t *target_len;
    size_t l_text, n_text;
    char *text;
} bam_hdr_t;

typedef struct {
    int32_t tid;
    int32_t pos;
    uint32_t bin:16, qual:8, l_qname:8;
    uint32_t flag:16, n_cigar:16;
    int32_t l_qseq;
    int32_t mtid;
    int32_t mpos;
    int32_t isize;
} bam1_core_t;

typedef struct {
    bam1_core_t core;
    int l_aux, l_data, m_data;
    uint8_t *data;
} bam1_t;

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

#define bam1_strand(b) (((b)->core.flag&BAM_FREVERSE) != 0)
#define bam1_mstrand(b) (((b)->core.flag&BAM_FMREVERSE) != 0)
#define bam1_cigar(b) ((uint32_t*)((b)->data + (b)->core.l_qname))
#define bam1_qname(b) ((char*)((b)->data))
#define bam1_seq(b) ((b)->data + (b)->core.n_cigar*4 + (b)->core.l_qname)
#define bam1_qual(b) ((b)->data + (b)->core.n_cigar*4 + (b)->core.l_qname + (((b)->core.l_qseq + 1)>>1))
#define bam1_seqi(s, i) ((s)[(i)/2] >> 4*(1-(i)%2) & 0xf)
#define bam1_aux(b) ((b)->data + (b)->core.n_cigar*4 + (b)->core.l_qname + (b)->core.l_qseq + ((b)->core.l_qseq + 1)/2)
#define bam1_qual(b) ((b)->data + (b)->core.n_cigar*4 + (b)->core.l_qname + (((b)->core.l_qseq + 1)>>1))
#define bam_init1() ((bam1_t*)calloc(1, sizeof(bam1_t)))
#define bam_destroy1(b) do { if (b) { free((b)->data); free(b); }} while (0)

extern int bam_is_be;

#ifdef __cplusplus
extern "C" {
#endif
    BGZF* bgzf_fdopen(int fd, const char* __restrict mode);
    BGZF* bgzf_open(const char* path, const char* __restrict mode);
    int bgzf_close(BGZF* fp);
    int bgzf_read(BGZF* fp, void* data, int length);
    int bgzf_write(BGZF* fp, const void* data, int length);
    int bgzf_flush(BGZF* fp);
    int bgzf_flush_try(BGZF *fp, int size);
    bam_hdr_t *bam_hdr_init(void);
    bam_hdr_t *bam_hdr_read(bamFile fp);
    int bam_hdr_write(bamFile fp, const bam_hdr_t *header);
    void bam_hdr_destroy(bam_hdr_t *header);
    int bam_read1(bamFile fp, bam1_t *b);
    uint8_t *bam_aux_get(const bam1_t *b, const char tag[2]);
    int bam_write1_core(bamFile fp, const bam1_core_t *c, int data_len, uint8_t *data);
    int bam_write1(bamFile fp, const bam1_t *b);
#ifdef __cplusplus
}
#endif

#endif
