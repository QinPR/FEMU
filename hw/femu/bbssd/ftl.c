#include "ftl.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

//#define FEMU_DEBUG_FTL

static void *ftl_thread(void *arg);


static int get_blk_index(struct ssd *ssd, struct ppa *ppa)
{
    /* get the global index of a block in whole ssd */
    struct ssdparams *spp = &ssd->sp;
    return ppa->g.blk + ppa->g.pl * spp->blks_per_pl + ppa->g.lun * spp->blks_per_lun + ppa->g.ch * spp->blks_per_ch;
}


static uint8_t check_block_status(struct ssd *ssd, int blk_index)
{
    /* Given the global index of block, check whether the block worn out */
    return ssd->blk_status_list[blk_index];
}


static inline bool should_gc(struct ssd *ssd)
{
    return (ssd->lm.free_line_cnt <= ssd->sp.gc_thres_lines);
}

static inline bool should_gc_high(struct ssd *ssd)
{
    return (ssd->lm.free_line_cnt <= ssd->sp.gc_thres_lines_high);
}

static inline struct ppa get_maptbl_ent(struct ssd *ssd, uint64_t lpn)
{
    return ssd->maptbl[lpn];
}

static inline void set_maptbl_ent(struct ssd *ssd, uint64_t lpn, struct ppa *ppa)
{
    ftl_assert(lpn < ssd->sp.tt_pgs);
    ssd->maptbl[lpn] = *ppa;
}

static uint64_t ppa2pgidx(struct ssd *ssd, struct ppa *ppa)
{
    struct ssdparams *spp = &ssd->sp;
    uint64_t pgidx;

    pgidx = ppa->g.ch  * spp->pgs_per_ch  + \
            ppa->g.lun * spp->pgs_per_lun + \
            ppa->g.pl  * spp->pgs_per_pl  + \
            ppa->g.blk * spp->pgs_per_blk + \
            ppa->g.pg;

    ftl_assert(pgidx < spp->tt_pgs);

    return pgidx;
}

static inline uint64_t get_rmap_ent(struct ssd *ssd, struct ppa *ppa)
{
    uint64_t pgidx = ppa2pgidx(ssd, ppa);

    return ssd->rmap[pgidx];
}

/* set rmap[page_no(ppa)] -> lpn */
static inline void set_rmap_ent(struct ssd *ssd, uint64_t lpn, struct ppa *ppa)
{
    uint64_t pgidx = ppa2pgidx(ssd, ppa);

    ssd->rmap[pgidx] = lpn;
}

static inline int victim_line_cmp_pri(pqueue_pri_t next, pqueue_pri_t curr)
{
    return (next > curr);
}

static inline pqueue_pri_t victim_line_get_pri(void *a)
{
    return ((struct line *)a)->vpc;
}

static inline void victim_line_set_pri(void *a, pqueue_pri_t pri)
{
    ((struct line *)a)->vpc = pri;
}

static inline size_t victim_line_get_pos(void *a)
{
    return ((struct line *)a)->pos;
}

static inline void victim_line_set_pos(void *a, size_t pos)
{
    ((struct line *)a)->pos = pos;
}

static void ssd_init_lines(struct ssd *ssd)
{
    struct ssdparams *spp = &ssd->sp;
    struct line_mgmt *lm = &ssd->lm;
    struct line *line;

    lm->tt_lines = spp->blks_per_pl;
    ftl_assert(lm->tt_lines == spp->tt_lines);
    lm->lines = g_malloc0(sizeof(struct line) * lm->tt_lines);

    QTAILQ_INIT(&lm->free_line_list);
    lm->victim_line_pq = pqueue_init(spp->tt_lines, victim_line_cmp_pri,
            victim_line_get_pri, victim_line_set_pri,
            victim_line_get_pos, victim_line_set_pos);
    QTAILQ_INIT(&lm->full_line_list);

    lm->free_line_cnt = 0;
    for (int i = 0; i < lm->tt_lines; i++) {
        line = &lm->lines[i];
        line->id = i;
        line->ipc = 0;
        line->vpc = 0;
        line->pos = 0;
        /* initialize all the lines as free lines */
        QTAILQ_INSERT_TAIL(&lm->free_line_list, line, entry);
        lm->free_line_cnt++;
    }

    ftl_assert(lm->free_line_cnt == lm->tt_lines);
    lm->victim_line_cnt = 0;
    lm->full_line_cnt = 0;
}

static void ssd_init_write_pointer(struct ssd *ssd)
{
    struct write_pointer *wpp = &ssd->wp;
    struct line_mgmt *lm = &ssd->lm;
    struct line *curline = NULL;

    curline = QTAILQ_FIRST(&lm->free_line_list);
    QTAILQ_REMOVE(&lm->free_line_list, curline, entry);
    lm->free_line_cnt--;
    printf("[FEMU Dbg] <ssd_init_write_pointer>: the free line cnt is: %d and gc threshold is %d.\n", lm->free_line_cnt,  ssd->sp.gc_thres_lines_high);

    /* wpp->curline is always our next-to-write super-block */
    wpp->curline = curline;
    wpp->ch = 0;
    wpp->lun = 0;
    wpp->pg = 0;
    wpp->blk = 0;
    wpp->pl = 0;
}


static void ssd_init_blks_status(struct ssd *ssd)
{
    struct ssdparams *spp = &ssd->sp;
    ssd->blk_status_list = g_malloc0(sizeof(uint8_t) * spp->total_blk_num); 
    
    uint8_t init_block_status = BLOCK_ALIVE;
    for (int i = 0; i < spp->total_blk_num; i++) {
        ssd->blk_status_list[i] = init_block_status;
    }
}


static inline void check_addr(int a, int max)
{
    ftl_assert(a >= 0 && a < max);
}

static struct line *get_next_free_line(struct ssd *ssd)
{
    struct line_mgmt *lm = &ssd->lm;
    struct line *curline = NULL;

    curline = QTAILQ_FIRST(&lm->free_line_list);    // QTAILQ_FIRST is to get the first element in queue
    if (!curline) {
        ftl_err("No free lines left in [%s] !!!!\n", ssd->ssdname);
        return NULL;
    }

    QTAILQ_REMOVE(&lm->free_line_list, curline, entry);
    lm->free_line_cnt--;
    return curline;
}


static void check_params(struct ssdparams *spp)
{
    /*
     * we are using a general write pointer increment method now, no need to
     * force luns_per_ch and nchs to be power of 2
     */

    //ftl_assert(is_power_of_2(spp->luns_per_ch));
    //ftl_assert(is_power_of_2(spp->nchs));
}

static void ssd_init_params(struct ssdparams *spp)
{
    spp->secsz = 512;
    spp->secs_per_pg = 8;
    spp->pgs_per_blk = 256;
    spp->blks_per_pl = 256; /* 16GB */
    spp->pls_per_lun = 1;
    spp->luns_per_ch = 8;
    spp->nchs = 8;

    spp->pg_rd_lat = NAND_READ_LATENCY;
    spp->pg_wr_lat = NAND_PROG_LATENCY;
    spp->blk_er_lat = NAND_ERASE_LATENCY;
    spp->ch_xfer_lat = 0;



    /* calculated values */
    spp->secs_per_blk = spp->secs_per_pg * spp->pgs_per_blk;
    spp->secs_per_pl = spp->secs_per_blk * spp->blks_per_pl;
    spp->secs_per_lun = spp->secs_per_pl * spp->pls_per_lun;
    spp->secs_per_ch = spp->secs_per_lun * spp->luns_per_ch;
    spp->tt_secs = spp->secs_per_ch * spp->nchs;

    spp->pgs_per_pl = spp->pgs_per_blk * spp->blks_per_pl;
    spp->pgs_per_lun = spp->pgs_per_pl * spp->pls_per_lun;
    spp->pgs_per_ch = spp->pgs_per_lun * spp->luns_per_ch;
    spp->tt_pgs = spp->pgs_per_ch * spp->nchs;

    spp->blks_per_lun = spp->blks_per_pl * spp->pls_per_lun;
    spp->blks_per_ch = spp->blks_per_lun * spp->luns_per_ch;
    spp->tt_blks = spp->blks_per_ch * spp->nchs;

    spp->pls_per_ch =  spp->pls_per_lun * spp->luns_per_ch;
    spp->tt_pls = spp->pls_per_ch * spp->nchs;

    spp->tt_luns = spp->luns_per_ch * spp->nchs;

    spp->total_blk_num = spp->blks_per_pl * spp->pls_per_lun * spp->luns_per_ch * spp->nchs;

    /* line is special, put it at the end */
    spp->blks_per_line = spp->tt_luns; /* TODO: to fix under multiplanes */
    spp->pgs_per_line = spp->blks_per_line * spp->pgs_per_blk;
    spp->secs_per_line = spp->pgs_per_line * spp->secs_per_pg;
    spp->tt_lines = spp->blks_per_lun; /* TODO: to fix under multiplanes */

    // spp->gc_thres_pcent = 0.75;
    spp->gc_thres_pcent = 0.02;
    spp->gc_thres_lines = (int)((1 - spp->gc_thres_pcent) * spp->tt_lines);    // `tt_lines` means number of lines = blks_per_pl
    // spp->gc_thres_pcent_high = 0.95;
    spp->gc_thres_pcent_high = 0.02;
    spp->gc_thres_lines_high = (int)((1 - spp->gc_thres_pcent_high) * spp->tt_lines);
    spp->enable_gc_delay = true;


    check_params(spp);
}

static void ssd_init_nand_page(struct nand_page *pg, struct ssdparams *spp)
{
    pg->nsecs = spp->secs_per_pg;
    pg->sec = g_malloc0(sizeof(nand_sec_status_t) * pg->nsecs);
    for (int i = 0; i < pg->nsecs; i++) {
        pg->sec[i] = SEC_FREE;
    }
    pg->status = PG_FREE;
}

static void ssd_init_nand_blk(struct nand_block *blk, struct ssdparams *spp)
{
    blk->npgs = spp->pgs_per_blk;
    blk->pg = g_malloc0(sizeof(struct nand_page) * blk->npgs);
    for (int i = 0; i < blk->npgs; i++) {
        ssd_init_nand_page(&blk->pg[i], spp);
    }
    blk->ipc = 0;
    blk->vpc = 0;
    blk->erase_cnt = 0;
    blk->wp = 0;

    // Belowing initialization is about wear leveling
    blk -> fail_possibility = BLOCK_FAILURE_POSSIBILTY;
    blk -> max_wr_count = BLOCK_MAX_WR_COUNT;
    blk -> cur_wr_count = 0;
    blk -> block_wear_status = BLOCK_ALIVE;
}

static void ssd_init_nand_plane(struct nand_plane *pl, struct ssdparams *spp)
{
    pl->nblks = spp->blks_per_pl;
    pl->blk = g_malloc0(sizeof(struct nand_block) * pl->nblks);
    for (int i = 0; i < pl->nblks; i++) {
        ssd_init_nand_blk(&pl->blk[i], spp);
    }
}

static void ssd_init_nand_lun(struct nand_lun *lun, struct ssdparams *spp)
{
    lun->npls = spp->pls_per_lun;
    lun->pl = g_malloc0(sizeof(struct nand_plane) * lun->npls);
    for (int i = 0; i < lun->npls; i++) {
        ssd_init_nand_plane(&lun->pl[i], spp);
    }
    lun->next_lun_avail_time = 0;
    lun->busy = false;
}

static void ssd_init_ch(struct ssd_channel *ch, struct ssdparams *spp)
{
    ch->nluns = spp->luns_per_ch;
    ch->lun = g_malloc0(sizeof(struct nand_lun) * ch->nluns);
    for (int i = 0; i < ch->nluns; i++) {
        ssd_init_nand_lun(&ch->lun[i], spp);
    }
    ch->next_ch_avail_time = 0;
    ch->busy = 0;
}

static void ssd_init_maptbl(struct ssd *ssd)
{
    struct ssdparams *spp = &ssd->sp;

    ssd->maptbl = g_malloc0(sizeof(struct ppa) * spp->tt_pgs);
    for (int i = 0; i < spp->tt_pgs; i++) {
        ssd->maptbl[i].ppa = UNMAPPED_PPA;
    }
}

static void ssd_init_rmap(struct ssd *ssd)
{
    struct ssdparams *spp = &ssd->sp;

    ssd->rmap = g_malloc0(sizeof(uint64_t) * spp->tt_pgs);
    for (int i = 0; i < spp->tt_pgs; i++) {
        ssd->rmap[i] = INVALID_LPN;
    }
}

void ssd_init(FemuCtrl *n)
{
    /*
        Initialize the SSD (to mimic the behavior of ssd.). The things to be initialized are:
        1. FTL(Flash Translation Layer): The FTL maps the higher-level disk block addresses to the lower-level flash block addresses, 
                                        so that the operating system can access data on the flash memory device as if it were accessing data on an HDD.
        2. Internal layout architecture.
        3. Maptbl (Page table mapping table).
        4. Rmap
        5. Lines
        6. Write pointer
        Finally it will create a qemu thread.
    */

    femu_debug("[ssd init] Initialize things of SSD to be emulated, such as FTL, Internal layout architecture, Maptbl... \n");

    struct ssd *ssd = n->ssd;
    struct ssdparams *spp = &ssd->sp;

    ftl_assert(ssd);    // Same as the `assert` function, exit when ssd is none.

    femu_debug("[ssd init] Initialize the parameters of SSD. \n");
    ssd_init_params(spp);

    /* initialize ssd internal layout architecture */
    femu_debug("[ssd init] Initialize the internal layout architecture of SSD. \n");
    ssd->ch = g_malloc0(sizeof(struct ssd_channel) * spp->nchs);     // ssd->ch stands for ssd channel.
    for (int i = 0; i < spp->nchs; i++) {
        ssd_init_ch(&ssd->ch[i], spp);
    }

    /* initialize maptbl */
    /*
        P.S. Page-level mapping tables are used in virtual memory systems to map virtual addresses to physical addresses. 
        The page-level mapping table is a data structure that is used by the operating system to keep track of the mapping 
        between virtual addresses and physical addresses. The page-level mapping table is typically implemented as a tree structure, 
        with each node in the tree representing a page of memory.
    */
    femu_debug("[ssd init] Initialize the Page-level mapping tables of SSD. \n");
    ssd_init_maptbl(ssd);

    /* initialize rmap */
    /*
        P.S. rmap stands for reserved page-level mapping table.
        It is used in some virtual memory systems to reserve a portion of the virtual address space for kernel use. 
        The reserved page-level mapping table is typically used to map the kernel’s code and data into the virtual address space. 
        This allows the kernel to access its own code and data without having to switch to a different address space1.
    */
    femu_debug("[ssd init] Initialize the reversed Page-level mapping tables of SSD. \n");
    ssd_init_rmap(ssd);

    /* initialize all the lines */
    /*
        P.S. In the context of SSDs, a line refers to a group of contiguous memory locations that are read or written together. 
        The size of a line is typically 512 bytes or 4 kilobytes1.
    */
    femu_debug("[ssd init] Initialize the line of SSD. \n");
    ssd_init_lines(ssd);

    /* initialize write pointer, this is how we allocate new pages for writes */
    /*
        P.S. The write pointer in the context of SSDs refers to the location in the flash memory where the next write operation will occur. 
        The write pointer is used by the SSD controller to manage the allocation of free space on the drive,
        and to ensure that writes are distributed evenly across the available memory.
    */
    femu_debug("[ssd init] Initialize the write pointer of SSD. \n");
    ssd_init_write_pointer(ssd);

    /* initialize the bitmap for blks status */
    ssd_init_blks_status(ssd);

    qemu_thread_create(&ssd->ftl_thread, "FEMU-FTL-Thread", ftl_thread, n,
                       QEMU_THREAD_JOINABLE);
}


static inline bool valid_ppa(struct ssd *ssd, struct ppa *ppa)
{
    struct ssdparams *spp = &ssd->sp;
    int ch = ppa->g.ch;
    int lun = ppa->g.lun;
    int pl = ppa->g.pl;
    int blk = ppa->g.blk;
    int pg = ppa->g.pg;
    int sec = ppa->g.sec;

    if (ch >= 0 && ch < spp->nchs && lun >= 0 && lun < spp->luns_per_ch && pl >=
        0 && pl < spp->pls_per_lun && blk >= 0 && blk < spp->blks_per_pl && pg
        >= 0 && pg < spp->pgs_per_blk && sec >= 0 && sec < spp->secs_per_pg)
        return true;

    return false;
}

static inline bool valid_lpn(struct ssd *ssd, uint64_t lpn)
{
    return (lpn < ssd->sp.tt_pgs);
}

static inline bool mapped_ppa(struct ppa *ppa)
{
    return !(ppa->ppa == UNMAPPED_PPA);
}

static inline struct ssd_channel *get_ch(struct ssd *ssd, struct ppa *ppa)
{
    return &(ssd->ch[ppa->g.ch]);
}

static inline struct nand_lun *get_lun(struct ssd *ssd, struct ppa *ppa)
{
    // is this kind of fetching data?
    struct ssd_channel *ch = get_ch(ssd, ppa);
    return &(ch->lun[ppa->g.lun]);
}

static inline struct nand_plane *get_pl(struct ssd *ssd, struct ppa *ppa)
{
    struct nand_lun *lun = get_lun(ssd, ppa);
    return &(lun->pl[ppa->g.pl]);
}

static inline struct nand_block *get_blk(struct ssd *ssd, struct ppa *ppa)
{
    /* This function retrieve a block according to physical page address */
    struct nand_plane *pl = get_pl(ssd, ppa);
    return &(pl->blk[ppa->g.blk]);
}

static inline struct line *get_line(struct ssd *ssd, struct ppa *ppa)
{
    return &(ssd->lm.lines[ppa->g.blk]);
}

static inline struct nand_page *get_pg(struct ssd *ssd, struct ppa *ppa)
{
    struct nand_block *blk = get_blk(ssd, ppa);
    return &(blk->pg[ppa->g.pg]);
}


static void update_block_status(struct ssd *ssd, struct ppa *ppa, uint8_t new_status)
{
    // update the bit map, so that it can be indexed in O(1) complexity.
    int blk_index = get_blk_index(ssd, ppa);
    ssd->blk_status_list[blk_index] = new_status;
}


static bool is_block_wear_out(struct ssd *ssd, struct ppa *ppa)
{
    struct nand_block *target_blk = get_blk(ssd, ppa);
    /* 1. check whether the block has already worn out */
    if (target_blk->block_wear_status == BLOCK_FAILURE) {
        printf("[FEMU Dbg] <wear_out_simulator> block has already worn before. \n");
        return true;
    }
    /* 2. check whether the block wears out by failure possibility */
    srand(qemu_clock_get_us(QEMU_CLOCK_REALTIME));        // initialize the random number generator with current time, so that it will generate different random number
    if (rand() % 1000000 < target_blk->fail_possibility) {
        printf("[FEMU Dbg] <wear_out_simulator> block wears out because of possibility. \n");
        return true;
    }
    /* 3. check whether the block wears out by max write count */
    if (target_blk->cur_wr_count >= target_blk->max_wr_count) {
        printf("[FEMU Dbg] <wear_out_simulator> block wears out because of reaching max_wr_count \n");
        return true;
    }
    return false;
}


static void wear_out_simulator(struct ssd *ssd, struct ppa *ppa)
{
    /*
        Simulate the wearing out of cell in a block.
        The block may wear out because of two reasons:
        1. Wear out by accident (the BLOCK_WEAR_OUT_POSSIBILITY in `ftl.h`).
        2. Reach the max # of write (the BLOCK_MAX_WR_COUNT in `ftl.h`).
    */
    struct nand_block *target_block = get_blk(ssd, ppa);
    if (is_block_wear_out(ssd, ppa)) {
        target_block -> block_wear_status = BLOCK_FAILURE;
        update_block_status(ssd, ppa, BLOCK_FAILURE);
        printf("[FEMU Dbg] <wear_out_simulator> The block (index = %d) wears out. \n", get_blk_index(ssd, ppa));
    }
}


static struct ppa get_new_page(struct ssd *ssd)
{
    struct write_pointer *wpp = &ssd->wp;
    struct ppa ppa;

    ppa.ppa = 0;
    ppa.g.ch = wpp->ch;
    ppa.g.lun = wpp->lun;
    ppa.g.pg = wpp->pg;
    ppa.g.blk = wpp->blk;
    ppa.g.pl = wpp->pl;
    ftl_assert(ppa.g.pl == 0);
    return ppa;
}


static void find_alive_block(struct ssd *ssd)
{
    /* 
        Will be called when the block that write_pointer wears out.
        This function will find an alive block in the same line, and make the write pointer point to it.
        If the same line is full, or all the blocks in the same line worn out, move to the new line.
    */
    struct ssdparams *spp = &ssd->sp;
    struct write_pointer *wpp = &ssd->wp;
    struct line_mgmt *lm = &ssd->lm;
    int original_ch = wpp->ch;
    int original_lun = wpp->lun;
    int original_pl = wpp->pl;

    do {
        wpp->ch++;
        if (wpp->ch == spp->nchs) {
            wpp->ch = 0;
            wpp->lun++;
            if (wpp->lun == spp->luns_per_ch) {
                wpp->lun = 0;
                int blk_index = wpp->curline->id + wpp->pl * spp->blks_per_pl + wpp->lun * spp->blks_per_lun + wpp->ch * spp->blks_per_ch;
                if (ssd->blk_status_list[blk_index] == BLOCK_ALIVE && wpp->pg != spp->pgs_per_blk) {
                    return;
                }
            }
        }
    } while (wpp->ch != original_ch && wpp->lun != original_lun && wpp->pl != original_pl);

    /* special situation that all the blocks in the same line wear out. make the write pointer point to next line. */
    wpp->pg = 0;
    if (wpp->curline->vpc == spp->pgs_per_line) {
        QTAILQ_INSERT_TAIL(&lm->full_line_list, wpp->curline, entry); 
        lm->full_line_cnt++;
    } else {
        pqueue_insert(lm->victim_line_pq, wpp->curline);
        lm->victim_line_cnt++;
    }
    wpp->curline = NULL;
    wpp->curline = get_next_free_line(ssd);   
    wpp->blk = wpp->curline->id;
}


static void ssd_wear_leveling_handler(struct ssd *ssd, struct ppa *cur_ppa, int Operation)
{
    /*
        handler the block wear out. return the latency of handling.
    */
    struct ssdparams *spp = &ssd->sp;
    struct write_pointer *wpp = &ssd->wp;
    struct line_mgmt *lm = &ssd->lm;
    struct ppa ppa_temp;
    int cur_blk_index;

    if (Operation == NAND_WRITE) {
        /* if Block wearing out during WRITE operation, choose to write to another block. */
        check_addr(wpp->ch, spp->nchs);
        wpp->ch++;
        if (wpp->ch == spp->nchs) {
            wpp->ch = 0;
            check_addr(wpp->lun, spp->luns_per_ch);
            wpp->lun++;
            /* in this case, we should go to next lun */
            if (wpp->lun == spp->luns_per_ch) {
                wpp->lun = 0;
                /* go to next page in the block */
                check_addr(wpp->pg, spp->pgs_per_blk);
                wpp->pg++;
                if (wpp->pg == spp->pgs_per_blk) {
                    wpp->pg = 0;
                    /* move current line to {victim,full} line list */
                    if (wpp->curline->vpc == spp->pgs_per_line) {
                        /* all pgs are still valid, move to full line list */
                        ftl_assert(wpp->curline->ipc == 0);
                        QTAILQ_INSERT_TAIL(&lm->full_line_list, wpp->curline, entry);
                        lm->full_line_cnt++;
                    } else {
                        ftl_assert(wpp->curline->vpc >= 0 && wpp->curline->vpc < spp->pgs_per_line);
                        /* there must be some invalid pages in this line */
                        ftl_assert(wpp->curline->ipc > 0);
                        pqueue_insert(lm->victim_line_pq, wpp->curline);
                        lm->victim_line_cnt++;
                    }
                    /* current line is used up, pick another empty line */
                    check_addr(wpp->blk, spp->blks_per_pl);
                    wpp->curline = NULL;
                    wpp->curline = get_next_free_line(ssd);   
                    if (!wpp->curline) {
                        /* TODO */
                        abort();
                    }
                    wpp->blk = wpp->curline->id;
                    check_addr(wpp->blk, spp->blks_per_pl);
                    /* make sure we are starting from page 0 in the super block */
                    ftl_assert(wpp->pg == 0);
                    ftl_assert(wpp->lun == 0);
                    ftl_assert(wpp->ch == 0);
                    /* TODO: assume # of pl_per_lun is 1, fix later */
                    ftl_assert(wpp->pl == 0);
                } 
            }
        }
        if (ENABLE_WEAR_OUT_SIMULATION) {
            while (true) {
                ppa_temp.ppa = 0;
                ppa_temp.g.ch = wpp->ch;
                ppa_temp.g.lun = wpp->lun;
                ppa_temp.g.pg = wpp->pg;
                ppa_temp.g.blk = wpp->blk;
                ppa_temp.g.pl = wpp->pl;
                cur_blk_index = get_blk_index(ssd, &ppa_temp);
                if (check_block_status(ssd, cur_blk_index) != BLOCK_ALIVE){
                    /* Find a block which is ALIVE in the same line. If the same line is full or all the blocks in the line worn out, move to new free line. */
                    find_alive_block(ssd);
                } else {
                    break; 
                }
            }
        }
    } else if (Operation == NAND_READ) {
        /* If Block wear out during READ Operation, choose to repair the current block to be read. */
        /* Repairation of the block to be read. */
        struct nand_block *target_block = get_blk(ssd, cur_ppa);
        target_block -> block_wear_status = BLOCK_ALIVE;
        target_block -> cur_wr_count = 0;
        update_block_status(ssd, cur_ppa, BLOCK_ALIVE);
    } else {
        printf("[FEMU Dbg] <ssd_wear_leveling_handler>: Operation not supported. \n");
    }
}


static uint64_t ssd_advance_status(struct ssd *ssd, struct ppa *ppa, struct
        nand_cmd *ncmd)
{
    int c = ncmd->cmd;
    uint64_t cmd_stime = (ncmd->stime == 0) ? \
        qemu_clock_get_ns(QEMU_CLOCK_REALTIME) : ncmd->stime;        // ncmd->stime is the request arrive time.
    uint64_t nand_stime;
    struct ssdparams *spp = &ssd->sp;
    struct nand_lun *lun = get_lun(ssd, ppa);     // get the specific lun (logical unit) (p.s. one channel may have multiple luns, one lun may have multiple planes, one plane may have multiple blocks)
    uint64_t lat = 0;

    switch (c) {
    case NAND_READ:
        /* read: perform NAND cmd first */
        nand_stime = (lun->next_lun_avail_time < cmd_stime) ? cmd_stime : \
                     lun->next_lun_avail_time;
        lun->next_lun_avail_time = nand_stime + spp->pg_rd_lat;
        lat = lun->next_lun_avail_time - cmd_stime;
#if 0
        lun->next_lun_avail_time = nand_stime + spp->pg_rd_lat;

        /* read: then data transfer through channel */
        chnl_stime = (ch->next_ch_avail_time < lun->next_lun_avail_time) ? \
            lun->next_lun_avail_time : ch->next_ch_avail_time;
        ch->next_ch_avail_time = chnl_stime + spp->ch_xfer_lat;

        lat = ch->next_ch_avail_time - cmd_stime;
#endif
        break;

    case NAND_WRITE:
        /* write: transfer data through channel first */
        nand_stime = (lun->next_lun_avail_time < cmd_stime) ? cmd_stime : \
                     lun->next_lun_avail_time;
        if (ncmd->type == USER_IO) {      // what is the difference between these two lines?
            lun->next_lun_avail_time = nand_stime + spp->pg_wr_lat;
        } else {
            lun->next_lun_avail_time = nand_stime + spp->pg_wr_lat;
        }
        lat = lun->next_lun_avail_time - cmd_stime;

#if 0
        chnl_stime = (ch->next_ch_avail_time < cmd_stime) ? cmd_stime : \
                     ch->next_ch_avail_time;
        ch->next_ch_avail_time = chnl_stime + spp->ch_xfer_lat;

        /* write: then do NAND program */
        nand_stime = (lun->next_lun_avail_time < ch->next_ch_avail_time) ? \
            ch->next_ch_avail_time : lun->next_lun_avail_time;
        lun->next_lun_avail_time = nand_stime + spp->pg_wr_lat;

        lat = lun->next_lun_avail_time - cmd_stime;
#endif
        break;

    case NAND_ERASE:
        /* erase: only need to advance NAND status */
        nand_stime = (lun->next_lun_avail_time < cmd_stime) ? cmd_stime : \
                     lun->next_lun_avail_time;
        lun->next_lun_avail_time = nand_stime + spp->blk_er_lat;

        lat = lun->next_lun_avail_time - cmd_stime;
        break;

    default:
        ftl_err("Unsupported NAND command: 0x%x\n", c);
    }

    return lat;
}


/* update SSD status about one page from PG_VALID -> PG_VALID */
static void mark_page_invalid(struct ssd *ssd, struct ppa *ppa)
{
    struct line_mgmt *lm = &ssd->lm;
    struct ssdparams *spp = &ssd->sp;
    struct nand_block *blk = NULL;
    struct nand_page *pg = NULL;
    bool was_full_line = false;
    struct line *line;

    /* update corresponding page status */
    pg = get_pg(ssd, ppa);
    ftl_assert(pg->status == PG_VALID);
    pg->status = PG_INVALID;

    /* update corresponding block status */
    blk = get_blk(ssd, ppa);    
    ftl_assert(blk->ipc >= 0 && blk->ipc < spp->pgs_per_blk);
    blk->ipc++;
    ftl_assert(blk->vpc > 0 && blk->vpc <= spp->pgs_per_blk);
    blk->vpc--;

    /* update corresponding line status */
    line = get_line(ssd, ppa);
    ftl_assert(line->ipc >= 0 && line->ipc < spp->pgs_per_line);
    if (line->vpc == spp->pgs_per_line) {
        ftl_assert(line->ipc == 0);
        was_full_line = true;
    }
    line->ipc++;
    ftl_assert(line->vpc > 0 && line->vpc <= spp->pgs_per_line);
    /* Adjust the position of the victime line in the pq under over-writes */
    if (line->pos) {
        /* Note that line->vpc will be updated by this call */
        pqueue_change_priority(lm->victim_line_pq, line->vpc - 1, line);
    } else {
        line->vpc--;
    }

    if (was_full_line) {
        /* move line: "full" -> "victim" */
        // QTAILQ_REMOVE(&lm->full_line_list, line, entry);   // bug here, TODO: to fix 
        lm->full_line_cnt--;
        pqueue_insert(lm->victim_line_pq, line);
        lm->victim_line_cnt++;
    }
}

static void mark_page_valid(struct ssd *ssd, struct ppa *ppa)
{
    struct nand_block *blk = NULL;
    struct nand_page *pg = NULL;
    struct line *line;

    /* update page status */
    pg = get_pg(ssd, ppa);
    ftl_assert(pg->status == PG_FREE);
    pg->status = PG_VALID;

    /* update corresponding block status */
    blk = get_blk(ssd, ppa);     // change here
    ftl_assert(blk->vpc >= 0 && blk->vpc < ssd->sp.pgs_per_blk);
    blk->vpc++;

    /* update corresponding line status */
    line = get_line(ssd, ppa);
    ftl_assert(line->vpc >= 0 && line->vpc < ssd->sp.pgs_per_line);
    line->vpc++;
}

static void mark_block_free(struct ssd *ssd, struct ppa *ppa)
{
    struct ssdparams *spp = &ssd->sp;
    struct nand_block *blk = get_blk(ssd, ppa);
    struct nand_page *pg = NULL;

    for (int i = 0; i < spp->pgs_per_blk; i++) {
        /* reset page status */
        pg = &blk->pg[i];
        ftl_assert(pg->nsecs == spp->secs_per_pg);
        pg->status = PG_FREE;
    }

    /* reset block status */
    ftl_assert(blk->npgs == spp->pgs_per_blk);
    blk->ipc = 0;
    blk->vpc = 0;
    blk->erase_cnt++;
}

static void gc_read_page(struct ssd *ssd, struct ppa *ppa)
{
    if (ENABLE_WEAR_OUT_SIMULATION) {
        /* Add simulation of block wear out during reading */
        while (true) {
            /* Add simulation of block wearing */
            wear_out_simulator(ssd, ppa);
            int blk_index = get_blk_index(ssd, ppa);
            if (check_block_status(ssd, blk_index) == BLOCK_ALIVE){
                /* Keep while looping until the block to be read is not worn out (ALIVE) */
                break;
            } else {
                /* handle the block wearing out: repair the current block */
                ssd_wear_leveling_handler(ssd, ppa, NAND_READ);
            }
        }
    }
    /* advance ssd status, we don't care about how long it takes */
    if (ssd->sp.enable_gc_delay) {
        struct nand_cmd gcr;
        gcr.type = GC_IO;
        gcr.cmd = NAND_READ;
        gcr.stime = 0;
        ssd_advance_status(ssd, ppa, &gcr);
    }
}

/* move valid page data (already in DRAM) from victim line to a new page */
static uint64_t gc_write_page(struct ssd *ssd, struct ppa *old_ppa)
{
    struct ppa new_ppa;
    struct nand_lun *new_lun;
    uint64_t lpn = get_rmap_ent(ssd, old_ppa);

    ftl_assert(valid_lpn(ssd, lpn));

    if (ENABLE_WEAR_OUT_SIMULATION) {
        // judge whether the block worn out once there is a WRITE Operation.
        while (true) {
            new_ppa = get_new_page(ssd);
            /* Add simulation of block wearing */
            wear_out_simulator(ssd, &new_ppa);
            /* handle the block wearing out: update the wr_pointer to an alive block */
            ssd_wear_leveling_handler(ssd, &new_ppa, NAND_WRITE);
            int blk_index = get_blk_index(ssd, &new_ppa);
            if (check_block_status(ssd, blk_index) == BLOCK_ALIVE){
                /* Keep while looping until the block to be written is not worn out (ALIVE) */
                break;
            } 
        } 
    }

    /* update maptbl */
    set_maptbl_ent(ssd, lpn, &new_ppa);
    /* update rmap */
    set_rmap_ent(ssd, lpn, &new_ppa);

    mark_page_valid(ssd, &new_ppa);

    if (ENABLE_WEAR_OUT_SIMULATION) {
        /* make the wr count += 1 */
        struct nand_block *target_block = get_blk(ssd, &new_ppa);
        target_block -> cur_wr_count += 1;
    }

    if (ssd->sp.enable_gc_delay) {
        struct nand_cmd gcw;
        gcw.type = GC_IO;
        gcw.cmd = NAND_WRITE;
        gcw.stime = 0;
        ssd_advance_status(ssd, &new_ppa, &gcw);
    }

    /* advance per-ch gc_endtime as well */
#if 0
    new_ch = get_ch(ssd, &new_ppa);
    new_ch->gc_endtime = new_ch->next_ch_avail_time;
#endif

    new_lun = get_lun(ssd, &new_ppa);
    new_lun->gc_endtime = new_lun->next_lun_avail_time;

    return 0;
}

static struct line *select_victim_line(struct ssd *ssd, bool force)
{
    struct line_mgmt *lm = &ssd->lm;
    struct line *victim_line = NULL;


    victim_line = pqueue_peek(lm->victim_line_pq);
    if (!victim_line) {
        return NULL;
    }

    if (!force && victim_line->ipc < ssd->sp.pgs_per_line / 8) {
        return NULL;
    }

    pqueue_pop(lm->victim_line_pq);
    victim_line->pos = 0;
    lm->victim_line_cnt--;

    /* victim_line is a danggling node now */
    return victim_line;
}

/* here ppa identifies the block we want to clean */
static void clean_one_block(struct ssd *ssd, struct ppa *ppa)
{
    struct ssdparams *spp = &ssd->sp;
    struct nand_page *pg_iter = NULL;
    int cnt = 0;

    for (int pg = 0; pg < spp->pgs_per_blk; pg++) {
        ppa->g.pg = pg;
        pg_iter = get_pg(ssd, ppa);
        /* there shouldn't be any free page in victim blocks */
        ftl_assert(pg_iter->status != PG_FREE);
        if (pg_iter->status == PG_VALID) {
            gc_read_page(ssd, ppa);
            /* delay the maptbl update until "write" happens */
            gc_write_page(ssd, ppa);
            cnt++;
        }
    }

    ftl_assert(get_blk(ssd, ppa)->vpc == cnt);
}

static void mark_line_free(struct ssd *ssd, struct ppa *ppa)
{
    struct line_mgmt *lm = &ssd->lm;
    struct line *line = get_line(ssd, ppa);
    line->ipc = 0;
    line->vpc = 0;
    /* move this line to free line list */
    QTAILQ_INSERT_TAIL(&lm->free_line_list, line, entry);
    lm->free_line_cnt++;
    
}

static int do_gc(struct ssd *ssd, bool force)
{
    struct line *victim_line = NULL;
    struct ssdparams *spp = &ssd->sp;
    struct nand_lun *lunp;
    struct ppa ppa;
    int ch, lun;
    victim_line = select_victim_line(ssd, force);
    if (!victim_line) {
        return -1;
    }

    ppa.g.blk = victim_line->id;
    printf("[FEMU Dbg] <do_gc>: GC-ing line = %d,victim line's invalid page count = %d,victim line count = %d,full line count = %d,free line count = %d\n", ppa.g.blk,
              victim_line->ipc, ssd->lm.victim_line_cnt, ssd->lm.full_line_cnt,
              ssd->lm.free_line_cnt);

    /* copy back valid data */
    for (ch = 0; ch < spp->nchs; ch++) {
        for (lun = 0; lun < spp->luns_per_ch; lun++) {
            ppa.g.ch = ch;
            ppa.g.lun = lun;
            ppa.g.pl = 0;
            lunp = get_lun(ssd, &ppa);
            clean_one_block(ssd, &ppa);
            mark_block_free(ssd, &ppa);

            if (spp->enable_gc_delay) {
                struct nand_cmd gce;
                gce.type = GC_IO;
                gce.cmd = NAND_ERASE;
                gce.stime = 0;
                ssd_advance_status(ssd, &ppa, &gce);
            }

            lunp->gc_endtime = lunp->next_lun_avail_time;
        }
    }

    /* update line status */
    mark_line_free(ssd, &ppa);

    return 0;
}

static uint64_t ssd_read(struct ssd *ssd, NvmeRequest *req)
{
    struct ssdparams *spp = &ssd->sp;      // 1. get the parameters of SSD
    uint64_t lba = req->slba;              // 2. get “sector logical block addressing”, which is a method used in storage systems to address individual sectors on a disk.
    int nsecs = req->nlb;             // 3. get “number of logical blocks”, which is a term used in storage systems to describe the number of logical blocks that are transferred in a single I/O operation. 
    struct ppa ppa;             
    uint64_t start_lpn = lba / spp->secs_per_pg;     // 4. get the starting logical page and ending logical page.
    uint64_t end_lpn = (lba + nsecs - 1) / spp->secs_per_pg;
    uint64_t lpn;
    uint64_t sublat, maxlat = 0;

    if (end_lpn >= spp->tt_pgs) {
        ftl_err("start_lpn=%"PRIu64",tt_pgs=%d\n", start_lpn, ssd->sp.tt_pgs);
    }

    /* normal IO read path */
    for (lpn = start_lpn; lpn <= end_lpn; lpn++) {
        ppa = get_maptbl_ent(ssd, lpn);       // 5. translate from logical page to physical page.
        if (!mapped_ppa(&ppa) || !valid_ppa(ssd, &ppa)) {
            //printf("%s,lpn(%" PRId64 ") not mapped to valid ppa\n", ssd->ssdname, lpn);
            //printf("Invalid ppa,ch:%d,lun:%d,blk:%d,pl:%d,pg:%d,sec:%d\n",
            //ppa.g.ch, ppa.g.lun, ppa.g.blk, ppa.g.pl, ppa.g.pg, ppa.g.sec);
            continue;
        }
        uint64_t wear_out_lat = 0;
        if (ENABLE_WEAR_OUT_SIMULATION) {
            while (true) {
                /* Add simulation of block wearing */
                wear_out_simulator(ssd, &ppa);
                int blk_index = get_blk_index(ssd, &ppa);
                if (check_block_status(ssd, blk_index) == BLOCK_ALIVE){
                    /* Keep while looping until the block to be written is not worn out (ALIVE) */
                    break;
                } else {
                    /* handle the block wearing out: repair the current block */
                    ssd_wear_leveling_handler(ssd, &ppa, NAND_READ);
                    wear_out_lat += RD_WEAR_OUT_LATENCY;
                }
            }
        }

        struct nand_cmd srd;
        srd.type = USER_IO;
        srd.cmd = NAND_READ;
        srd.stime = req->stime;
        if (ENABLE_WEAR_OUT_SIMULATION) {
            sublat = ssd_advance_status(ssd, &ppa, &srd) + wear_out_lat;
        } else {
            sublat = ssd_advance_status(ssd, &ppa, &srd);
        }
        if (wear_out_lat > 0) {
            // printf("[FEMU Dbg] <ssd_read>: the curlat now is %" PRIu64 " where wear_out_lat is %" PRIu64 " \n", sublat, wear_out_lat);
        }
        maxlat = (sublat > maxlat) ? sublat : maxlat;       // it is calculated to be max of each page, maybe because it is read parallely.
    }
    // printf("[FEMU Dbg] <ssd_read>: the max latency is %" PRIu64 " \n", maxlat);

    return maxlat;
}

static uint64_t ssd_write(struct ssd *ssd, NvmeRequest *req)
{
    uint64_t lba = req->slba;      // 1. get “sector logical block addressing”, which is a method used in storage systems to address individual sectors on a disk.
    struct ssdparams *spp = &ssd->sp;  // 2. retrieve the parameters of SSD.
    int len = req->nlb;  // 3. get “number of logical blocks”, which is a term used in storage systems to describe the number of logical blocks that are transferred in a single I/O operation.
    
    /*
        P.S. LPN stands for “logical page number” and is a term used in storage systems to describe the logical address of a page of data. 
        The LPN is used by the disk controller to translate the logical address into a physical location on the disk.
    */
    uint64_t start_lpn = lba / spp->secs_per_pg;     // 4. get the starting logical page and ending logical page.
    uint64_t end_lpn = (lba + len - 1) / spp->secs_per_pg;
    struct ppa ppa;
    uint64_t lpn;
    uint64_t curlat = 0, maxlat = 0;
    int r;

    if (end_lpn >= spp->tt_pgs) {
        ftl_err("start_lpn=%"PRIu64",tt_pgs=%d\n", start_lpn, ssd->sp.tt_pgs);
    }

    // 5. Do GC if needed.
    while (should_gc_high(ssd)) {
        /* perform GC here until !should_gc(ssd) */
        r = do_gc(ssd, true);
        if (r == -1)
            break;
    }

    for (lpn = start_lpn; lpn <= end_lpn; lpn++) {
        ppa = get_maptbl_ent(ssd, lpn);     // 5. physical page address and update the old message
        if (mapped_ppa(&ppa)) {   
            /* update old page information first */
            mark_page_invalid(ssd, &ppa);
            set_rmap_ent(ssd, INVALID_LPN, &ppa);
        }

        /* new write */  /* block wear out added */
        uint64_t wear_out_lat = 0;
        if (ENABLE_WEAR_OUT_SIMULATION) {
            while (true) {
                ppa = get_new_page(ssd);
                /* Add simulation of block wearing */
                wear_out_simulator(ssd, &ppa);
                /* handle the block wearing out: update the wr_pointer to an alive block */
                ssd_wear_leveling_handler(ssd, &ppa, NAND_WRITE);
                int blk_index = get_blk_index(ssd, &ppa);
                if (check_block_status(ssd, blk_index) == BLOCK_ALIVE){
                    /* Keep while looping until the block to be written is not worn out (ALIVE) */
                    break;
                } else if (check_block_status(ssd, blk_index) == BLOCK_FAILURE){
                    wear_out_lat += WR_WEAR_OUT_LATENCY;
                }
            } 
        }
        else {
            ppa = get_new_page(ssd);
        }

        /* update maptbl */
        set_maptbl_ent(ssd, lpn, &ppa);
        /* update rmap */
        set_rmap_ent(ssd, lpn, &ppa);
        // 6. bind logical page address to the newly allocated page.

        mark_page_valid(ssd, &ppa);

        if (ENABLE_WEAR_OUT_SIMULATION) {
            struct nand_block *target_block = get_blk(ssd, &ppa);
            target_block -> cur_wr_count += 1;
        } else {
            ssd_wear_leveling_handler(ssd, &ppa, NAND_WRITE);
        }

        struct nand_cmd swr;
        swr.type = USER_IO;
        swr.cmd = NAND_WRITE;
        swr.stime = req->stime;
        /* get latency statistics */
        if (ENABLE_WEAR_OUT_SIMULATION) {
            curlat = ssd_advance_status(ssd, &ppa, &swr) + wear_out_lat;    // emulate the WRITE latency + (possible) WEAR OUT latency.
        } else {
            curlat = ssd_advance_status(ssd, &ppa, &swr);
        }
        // if (wear_out_lat > 0) {
        //     // printf("[FEMU Dbg] <ssd_write>: the curlat now is %" PRIu64 " where wear_out_lat is %" PRIu64 " \n", curlat, wear_out_lat);
        // }
        maxlat = (curlat > maxlat) ? curlat : maxlat;
    }
    // printf("[FEMU Dbg] <ssd_write>: the max latency is %" PRIu64 " \n", maxlat);
    return maxlat;
}

static void *ftl_thread(void *arg)
{
    /*
        The key function of understanding the request flow of READ and WRITE Operation.
        contains:
        1. ssd_write
        2. ssd_read
    */
    FemuCtrl *n = (FemuCtrl *)arg;
    struct ssd *ssd = n->ssd;
    NvmeRequest *req = NULL;
    uint64_t lat = 0;
    int rc;
    int i;

    while (!*(ssd->dataplane_started_ptr)) {    // kind of emualating the resouse contend of dataplane?
        usleep(100000);
    }

    /* FIXME: not safe, to handle ->to_ftl and ->to_poller gracefully */
    ssd->to_ftl = n->to_ftl;
    ssd->to_poller = n->to_poller;

    while (1) {
        for (i = 1; i <= n->nr_pollers; i++) {    // is it what mentioned in paper? (using polling to check the status of device)
            if (!ssd->to_ftl[i] || !femu_ring_count(ssd->to_ftl[i]))
                continue;

            rc = femu_ring_dequeue(ssd->to_ftl[i], (void *)&req, 1);
            if (rc != 1) {
                printf("FEMU: FTL to_ftl dequeue failed\n");
            }

            ftl_assert(req);
            switch (req->cmd.opcode) {
            case NVME_CMD_WRITE:
                // femu_debug("[ftl_thread] prepare to do ssd WRITE. \n");
                lat = ssd_write(ssd, req);
                break;
            case NVME_CMD_READ:
                // femu_debug("[ftl_thread] prepare to do ssd READ. \n");
                lat = ssd_read(ssd, req);
                break;
            case NVME_CMD_DSM:
                // femu_debug("[ftl_thread] prepare to do ssd DSM. \n");
                lat = 0;
                break;
            default:
                femu_debug("[ftl_thread] do nothing. \n");
                //ftl_err("FTL received unkown request type, ERROR\n");
                ;
            }

            // attach the request latency (READ/ WRITE/ ...) back to the request.
            req->reqlat = lat;
            req->expire_time += lat;

            // enqueue the request back.
            rc = femu_ring_enqueue(ssd->to_poller[i], (void *)&req, 1);
            if (rc != 1) {
                ftl_err("FTL to_poller enqueue failed\n");
            }

            /* clean one line if needed (in the background) */
            if (should_gc(ssd)) {
                do_gc(ssd, false);
            }
        }
    }

    return NULL;
}

