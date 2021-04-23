#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdarg.h>
#include <ctype.h>
#include <limits.h>
#include <unistd.h>
#include <inttypes.h>

#include "ketopt.h"
#include "bamio.h"

#define VERSION "0.1.0"
#define basename(str) (strrchr(str, '/') ? strrchr(str, '/') + 1 : str)

#define MMLEN 45

static void usage(char *str);
static int get_nm(bam1_t *b);

void error(const char *format, ...)
{
	va_list ap;
	va_start(ap, format);
	vfprintf(stderr, format, ap);
	va_end(ap);
	exit(EXIT_FAILURE);
}

int main(int argc, char *argv[])
{
	int c, m = MMLEN;
	char *o = 0;
	ketopt_t opt = KETOPT_INIT;
	while ((c = ketopt(&opt, argc, argv, 1, "o:m:hv", 0)) >= 0)
	{
		if (c == 'h') usage(argv[0]);
		else if (c == 'v') return(puts(VERSION));
		else if (c == 'o') o = opt.arg;
		else if (c == 'm') m = atoi(opt.arg);
		else return(puts("Unknown option"));
	}
	bamFile fp = 0;
	if (argc - opt.ind < 1)
	{
		if (!isatty(fileno(stdin)))
			fp = bam_dopen(fileno(stdin), "rb");
		else
			usage(argv[0]);
	}
	else if (argc - opt.ind == 1)
	{
		if (strcmp(argv[opt.ind], "-") == 0)
		{
			if (!isatty(fileno(stdin)))
				fp = bam_dopen(fileno(stdin), "rb");
			else
				usage(argv[0]);
		}
		else if (strncmp(argv[opt.ind], "-h", 2) == 0 || strncmp(argv[opt.ind], "--h", 2) == 0)
			usage(argv[0]);
		else
			fp = bam_open(argv[opt.ind], "rb");
	}
	else
		usage(argv[0]);
	BGZF *fo = o ? bgzf_open(o, "w") : bgzf_fdopen(fileno(stdout), "w");
	if (!fp || !fo) return(puts("I/O error"));
	// check io
	if (argc - opt.ind == 1 && o && strcmp(argv[opt.ind], o) == 0)
		error("Error: conflicting input and output bams!\n");
	// write header
	bam_hdr_t *h = bam_hdr_read(fp);
	if (bam_hdr_write(fo, h) < 0)
		error("Error writing bam header!\n");
	bam_hdr_destroy(h);
	bam1_t *b = bam_init1();
	while (bam_read1(fp, b) >= 0)
	{
		if (b->core.flag & (BAM_FQCFAIL | BAM_FSECONDARY | BAM_FSUPPLEMENTARY))
			continue;
		if (b->core.flag & BAM_FUNMAP || get_nm(b) < m)
			if (bam_write1(fo, b) < 0)
				error("Error writing bam\n");
	}
	bam_destroy1(b);
	bam_close(fp);
	bam_close(fo);
	return 0;
}

static int get_nm(bam1_t *b)
{
	int i = 0;
	uint32_t *c = bam1_cigar(b), m = 0;
	int o = b->core.n_cigar;
	for(i = 0; i < o; i++)
	{
		int32_t op = bam_cigar_op(c[i]);
		int32_t ol = bam_cigar_oplen(c[i]);
		switch(op)
		{
			case BAM_CMATCH: // M is match or mismatch!
				m += ol;
				break;
			default:
				break;
		}
	}
	return m;
}

static void usage(char *str)
{
	putchar('\n');
	puts("Filter bam by aligned length");
	putchar('\n');
	fprintf(stdout, "Usage: \e[1;31m%s\e[0;0m [options] <bam>\n", basename(str));
	putchar('\n');
	puts("  -o STR output bam file [stdout]");
	printf("  -m INT max aligned length [%d]", MMLEN);
	putchar('\n');
	puts("  -h     print help");
	puts("  -v     print version");
	putchar('\n');
	exit(1);
}
