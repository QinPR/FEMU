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
    lm->user_lines = spp->users_blks_per_pl;
    ftl_assert(lm->tt_lines == spp->tt_lines);
    lm->lines = g_malloc0(sizeof(struct line) * lm->tt_lines);

    QTAILQ_INIT(&lm->free_line_list);
    QTAILQ_INIT(&lm->free_over_provisioning_line_list);
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
        if (i < lm->user_lines) {
            line->block_type = USER_BLOCK;
            QTAILQ_INSERT_TAIL(&lm->free_line_list, line, entry);
        } else {
            line->block_type = OVER_PROVISIONING_BLOCK;
            QTAILQ_INSERT_TAIL(&lm->free_over_provisioning_line_list, line, entry);
        }
        lm->free_line_cnt++;
    }

    ftl_assert(lm->free_line_cnt == lm->tt_lines);
    lm->victim_line_cnt = 0;
    lm->full_line_cnt = 0;
}

static void ssd_init_write_pointer(struct ssd *ssd)
{
    /* 1. For the user blocks' write pointer */
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

    /* 2. For the over-provisioning blocks' write pointer */
    struct write_pointer *over_provisioning_wpp = &ssd->over_provisioning_wp;
    struct ssdparams *spp = &ssd->sp;
    struct line *over_provisioning_curline = NULL;

    over_provisioning_curline = QTAILQ_FIRST(&lm->free_over_provisioning_line_list);
    QTAILQ_REMOVE(&lm->free_over_provisioning_line_list, over_provisioning_curline, entry);
    lm->free_line_cnt--;
    over_provisioning_wpp->curline = over_provisioning_curline;
    over_provisioning_wpp->pg = -1;
    over_provisioning_wpp->blk = spp->users_blks_per_pl;  // the over-provisioning blocks are after the users' block
}


static void ssd_init_blks_status(struct ssd *ssd)
{
    struct ssdparams *spp = &ssd->sp;
    ssd->blk_status_list = g_malloc0(sizeof(uint8_t) * spp->total_blk_num); 
    
    uint8_t init_block_status = HEALTHY_BLOCK;
    for (int i = 0; i < spp->total_blk_num; i++) {
        ssd->blk_status_list[i] = init_block_status;
    }
}


static inline void check_addr(int a, int max)
{
    ftl_assert(a >= 0 && a < max);
}

static struct line *get_next_free_line(struct ssd *ssd, int blk_type)
{
    struct line_mgmt *lm = &ssd->lm;
    struct line *curline = NULL;
    int selected_block = blk_type;

    if (blk_type == USER_BLOCK) {
        /* user cannot write access the over-provisioning blocks */
        curline = QTAILQ_FIRST(&lm->free_line_list);
    } else if (blk_type == OVER_PROVISIONING_BLOCK) {
        curline = QTAILQ_FIRST(&lm->free_over_provisioning_line_list); 
        if (!curline) {   // if there is no free line in over_provisioning areas, use user blocks.
            curline = QTAILQ_FIRST(&lm->free_line_list);
            selected_block = USER_BLOCK;
        }
    } else {
        ftl_err("Unsupportted block type. \n");
    }

    if (!curline) {
        ftl_err("No suitable free lines left in [%s] !!!!\n", ssd->ssdname);
        return NULL;
    }

    if (selected_block == USER_BLOCK) {
        QTAILQ_REMOVE(&lm->free_line_list, curline, entry);
    } else if (selected_block == OVER_PROVISIONING_BLOCK) {
        QTAILQ_REMOVE(&lm->free_over_provisioning_line_list, curline, entry);
    } else {
        ftl_err("Unsupportted block type. \n");
    }
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

    /* A block in a plane are comprised of: 1. blocks for user space, 2. over-provisioning blocks */
    spp->users_blks_per_pl = 256;     /* 16GB */ 
    spp->overprovisioning_blks_per_pl = 32;  /* 2GB for over-provisioning */ 
    spp->blks_per_pl = spp->users_blks_per_pl + spp->overprovisioning_blks_per_pl; /* users blks + over-overprovisioning blks */

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

static void ssd_init_nand_blk(struct nand_block *blk, struct ssdparams *spp, int blk_type)
{
    blk->npgs = spp->pgs_per_blk;
    blk->pg = g_malloc0(sizeof(struct nand_page) * blk->npgs);
    for (int i = 0; i < blk->npgs; i++) {
        ssd_init_nand_page(&blk->pg[i], spp);
    }
    blk->ipc = 0;
    blk->vpc = 0;
    blk->wp = 0;
    blk->block_type = blk_type;

    // Belowing initialization is about wear leveling
    blk -> fail_possibility = BLOCK_WEAR_OUT_POSSIBILITY;
    srand(qemu_clock_get_us(QEMU_CLOCK_REALTIME));
    // Randomly generate the max PE count (within the range specifed in ftl.h) of a block. 
    blk -> max_PE_count = rand() % (BLOCK_MAX_PE_COUNT_HIGH - BLOCK_MAX_PE_COUNT_LOW) + BLOCK_MAX_PE_COUNT_HIGH;
    blk -> cur_PE_count = 0;
    blk -> block_wear_status = HEALTHY_BLOCK;
}

static void ssd_init_nand_plane(struct nand_plane *pl, struct ssdparams *spp)
{
    pl->nblks = spp->blks_per_pl;
    pl->n_users_blks = spp->users_blks_per_pl;
    pl->n_over_provisioning_blks = spp->overprovisioning_blks_per_pl;

    pl->blk = g_malloc0(sizeof(struct nand_block) * pl->nblks);
    for (int i = 0; i < pl->n_users_blks; i++) {
        ssd_init_nand_blk(&pl->blk[i], spp, USER_BLOCK);
    }
    for (int i = pl->n_users_blks; i < pl->nblks; i++) {
        ssd_init_nand_blk(&pl->blk[i], spp, OVER_PROVISIONING_BLOCK);
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
    if (target_blk->block_wear_status == BAD_BLOCK) {
        // printf("[FEMU Dbg] <is_block_wear_out> block has already worn before. \n");
        return true;
    }
    /* 2. check whether the block wears out by failure possibility */
    srand(qemu_clock_get_us(QEMU_CLOCK_REALTIME));        // initialize the random number generator with current time, so that it will generate different random number
    if (rand() % 1000000 < target_blk->fail_possibility) {
        // printf("[FEMU Dbg] <is_block_wear_out> block wears out because of possibility. \n");
        return true;
    }
    /* 3. check whether the block wears out because of reaching its max PE count */
    if (target_blk->cur_PE_count >= target_blk->max_PE_count) {
        printf("[FEMU Dbg] <is_block_wear_out> block wears out because of reaching max_PE_count \n");
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
        2. Reach the max # of write (the BLOCK_MAX_PE_COUNT_LOW in `ftl.h`).
    */
    struct nand_block *target_block = get_blk(ssd, ppa);
    if (is_block_wear_out(ssd, ppa)) {
        target_block -> block_wear_status = BAD_BLOCK;
        // printf("[FEMU Dbg] <wear_out_simulator> The block (index = %d), (blk = %d, lun = %d, pl = %d, ch = %d) wears out. \n", get_blk_index(ssd, ppa), ppa->g.blk, ppa->g.lun, ppa->g.pl, ppa->g.ch);
    }
}


static struct ppa get_new_page(struct ssd *ssd, int Operation)
{
    struct write_pointer *wpp = NULL;
    if (Operation == NAND_WRITE) {
        wpp = &ssd->wp;
    } else if (Operation == OVER_PROVISIONING_WRITE) {
        wpp = &ssd->over_provisioning_wp;
    } else {
        ftl_err("Unsupported Operation. \n");
    }
    struct ppa ppa;

    ppa.ppa = 0;
    ppa.g.ch = wpp->ch;
    ppa.g.lun = wpp->lun;
    ppa.g.pg = wpp->pg;
    ppa.g.blk = wpp->blk;
    ppa.g.pl = wpp->pl;
    return ppa;
}


static void find_healthy_block(struct ssd *ssd)
{
    /* 
        Will be called when the block that write_pointer wears out.
        This function will find an healthy block in the same line, and make the write pointer point to it.
        If the same line is full, or all the blocks in the same line worn out, move to the new line.
    */
    struct ssdparams *spp = &ssd->sp;
    struct write_pointer *wpp = &ssd->wp;
    struct line_mgmt *lm = &ssd->lm;
    int original_ch = wpp->ch;
    int original_lun = wpp->lun;
    int original_pl = wpp->pl;
    int blk_index = 0;

    do {
        wpp->ch++;
        if (wpp->ch == spp->nchs) {
            wpp->ch = 0;
            wpp->lun++;
            if (wpp->lun == spp->luns_per_ch) {
                wpp->lun = 0;
            }
        }
        blk_index = wpp->curline->id + wpp->pl * spp->blks_per_pl + wpp->lun * spp->blks_per_lun + wpp->ch * spp->blks_per_ch;
        // If block is healthy, return.
        if (ssd->blk_status_list[blk_index] == HEALTHY_BLOCK && wpp->pg != spp->pgs_per_blk) {
            return;
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
    wpp->curline = get_next_free_line(ssd, USER_BLOCK);   
    wpp->blk = wpp->curline->id;
}


static void ssd_wear_leveling_handler(struct ssd *ssd, int Operation, struct ppa *old_ppa)
{
    /*
        handler the block wear out. return the latency of handling.
    */
    struct ssdparams *spp = &ssd->sp;
    struct write_pointer *wpp = &ssd->wp;
    struct write_pointer *over_provisioning_wpp = &ssd->over_provisioning_wp;
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
                    wpp->curline = get_next_free_line(ssd, USER_BLOCK);   
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
                if (check_block_status(ssd, cur_blk_index) == BAD_BLOCK){
                    /* Find a block which is HEALTHY_BLOCK in the same line. If the same line is full or all the blocks in the line worn out, move to new free line. */
                    find_healthy_block(ssd);
                } else {
                    break; 
                }
            }
        }
        if (wpp->blk > spp->users_blks_per_pl) {
            ftl_err("user cannot reach the overprovisioning blocks. \n");
        }
    } else if (Operation == OVER_PROVISIONING_WRITE) {
        /*
            Find an over-provisioning block, which is inside the same plane with the bad block.
            Update the write pointer to this place.
        */
        struct nand_page *pg = NULL;
        ppa_temp.ppa = 0;
        ppa_temp.g.ch = old_ppa->g.ch;
        ppa_temp.g.lun = old_ppa->g.lun;
        ppa_temp.g.pl = old_ppa->g.pl;

        ppa_temp.g.blk = over_provisioning_wpp->blk;
        if (over_provisioning_wpp->ch == old_ppa->g.ch && over_provisioning_wpp->lun == old_ppa->g.lun && over_provisioning_wpp->pl == old_ppa->g.pl) {
            // the same over-provisioning block
            ppa_temp.g.pg = over_provisioning_wpp->pg + 1;
        } else {
            // move to a over-provisioning block in new plane
            ppa_temp.g.pg = 0;     
        }
        while (ppa_temp.g.pg < spp->pgs_per_blk) {
            pg = get_pg(ssd, &ppa_temp);
            if (pg->status == PG_FREE) {
                break;
            } 
            ppa_temp.g.pg++;
            if (ppa_temp.g.pg == spp->pgs_per_blk) {
                // If the pages in block are all filled, move to the next free line.
                ppa_temp.g.pg = 0;
                if (over_provisioning_wpp->curline->vpc == spp->pgs_per_line) {
                    QTAILQ_INSERT_TAIL(&lm->full_line_list, over_provisioning_wpp->curline, entry);
                    lm->full_line_cnt++;
                } else {
                    pqueue_insert(lm->victim_line_pq, over_provisioning_wpp->curline);
                    lm->victim_line_cnt++;
                }
                over_provisioning_wpp->curline = get_next_free_line(ssd, OVER_PROVISIONING_BLOCK);   
                over_provisioning_wpp->blk = over_provisioning_wpp->curline->id;
                ppa_temp.g.blk = over_provisioning_wpp->blk;
            }
        }
        over_provisioning_wpp->ch = old_ppa->g.ch;
        over_provisioning_wpp->lun = old_ppa->g.lun;
        over_provisioning_wpp->pl = old_ppa->g.pl;
        over_provisioning_wpp->pg = ppa_temp.g.pg;
        if (ENABLE_WEAR_OUT_SIMULATION) {
            while (true) {
                cur_blk_index = get_blk_index(ssd, &ppa_temp);
                // printf("Relocate to the over-provisioning block, index = %d. (blk = %d, lun = %d, pl = %d, ch = %d) \n", cur_blk_index, ppa_temp.g.blk, ppa_temp.g.lun, ppa_temp.g.pl, ppa_temp.g.ch);
                if (check_block_status(ssd, cur_blk_index) == BAD_BLOCK){
                    /* If the over-provisioning block is also a bad block, move to next over-provisioning block. */
                    // printf("The over-provisioning block (%d) is also a bad block. \n", cur_blk_index);
                    over_provisioning_wpp->pg = 0;
                    if (over_provisioning_wpp->curline->vpc == spp->pgs_per_line) {
                        QTAILQ_INSERT_TAIL(&lm->full_line_list, over_provisioning_wpp->curline, entry);
                        lm->full_line_cnt++;
                    } else {
                        pqueue_insert(lm->victim_line_pq, over_provisioning_wpp->curline);
                        lm->victim_line_cnt++;
                    }
                    over_provisioning_wpp->curline = get_next_free_line(ssd, OVER_PROVISIONING_BLOCK);   
                    over_provisioning_wpp->blk = over_provisioning_wpp->curline->id;
                    ppa_temp.g.blk = over_provisioning_wpp->blk;
                } else {
                    break;
                }
            }
        }
    } else {
        ftl_err("Operation not supported. \n");
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


static uint64_t rellocate_read_page(struct ssd *ssd, struct ppa *ppa, uint64_t stime)
{
    uint64_t lat = 0;
    struct nand_cmd srd;
    srd.type = USER_IO;
    srd.cmd = NAND_READ;
    srd.stime = stime;
    lat = ssd_advance_status(ssd, ppa, &srd);
    return lat;
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


static uint64_t rellocate_write_page(struct ssd *ssd, struct ppa *old_ppa, uint64_t stime)
{
    struct ppa new_ppa;
    uint64_t lat = 0;
    uint64_t lpn = get_rmap_ent(ssd, old_ppa);

    new_ppa = get_new_page(ssd, OVER_PROVISIONING_WRITE);

    ssd_wear_leveling_handler(ssd, OVER_PROVISIONING_WRITE, old_ppa); /* Update the write pointer to a over-provisioning block */

    /* update maptbl */
    set_maptbl_ent(ssd, lpn, &new_ppa);
    /* update rmap */
    set_rmap_ent(ssd, lpn, &new_ppa);

    mark_page_valid(ssd, &new_ppa);

    struct nand_cmd swr;
    swr.type = USER_IO;
    swr.cmd = NAND_WRITE;
    swr.stime = stime;
    lat = ssd_advance_status(ssd, &new_ppa, &swr);
    return lat;
}


static uint64_t bad_block_management(struct ssd *ssd, struct ppa *ppa, uint64_t req_stime)
{
    /* 
        When the block wears out, conduct the following steps of bad block management:
    */
    uint64_t stime = req_stime;
    uint64_t BBM_lat = 0;
    struct ssdparams *spp = &ssd->sp;
    struct nand_lun *lun = get_lun(ssd, ppa);

    update_block_status(ssd, ppa, BAD_BLOCK);    // 1. Update the bad block list.
    struct nand_page *pg_iter = NULL;
    /* for all the valid pages inside the bad block, relocate it to the overprovisioning block within same plane */
    ssd_wear_leveling_handler(ssd, OVER_PROVISIONING_WRITE, ppa); /* Update the write pointer to a over-provisioning block */
    for (int pg = 0; pg < spp->pgs_per_blk; pg++) {
        ppa->g.pg = pg;
        pg_iter = get_pg(ssd, ppa);
        if (pg_iter->status == PG_VALID) {
            BBM_lat += rellocate_read_page(ssd, ppa, stime);   // 2. Read all the valid data of target block out.
            stime = lun->next_lun_avail_time;   // the next req will be executed at the completion time of current read operation
            BBM_lat += DATA_CORRECT_LAT;                     // 3. Correct the error data.
            BBM_lat += rellocate_write_page(ssd, ppa, stime);  // 4. Rellocate the data to a healthy block.
            stime = req_stime + BBM_lat;   // the next req will be executed at the completion time of current write operation
        }
    }
    return BBM_lat;
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


static void mark_block_free(struct ssd *ssd, struct ppa *ppa)
{
    /*
        Simulation of ERASE operation.
    */
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
    blk->cur_PE_count++;
}

static void gc_read_page(struct ssd *ssd, struct ppa *ppa)
{
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

    ssd_wear_leveling_handler(ssd, OVER_PROVISIONING_WRITE, old_ppa);    /* Update the write pointer to a over-provisioning block (intend to enhance GC efficiency). */

    new_ppa = get_new_page(ssd, OVER_PROVISIONING_WRITE);

    /* update maptbl */
    set_maptbl_ent(ssd, lpn, &new_ppa);
    /* update rmap */
    set_rmap_ent(ssd, lpn, &new_ppa);

    mark_page_valid(ssd, &new_ppa);

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
    if (line->block_type == USER_BLOCK) {
        QTAILQ_INSERT_TAIL(&lm->free_line_list, line, entry);
    } else if (line->block_type == OVER_PROVISIONING_BLOCK) {
        QTAILQ_INSERT_TAIL(&lm->free_over_provisioning_line_list, line, entry);
    }
    lm->free_line_cnt++;
}

static int do_gc(struct ssd *ssd, bool force)
{
    struct line *victim_line = NULL;
    struct ssdparams *spp = &ssd->sp;
    struct nand_lun *lunp;
    struct ppa ppa;
    int ch, lun;
    int blk_index;
    victim_line = select_victim_line(ssd, force);
    if (!victim_line) {
        return -1;
    }

    ppa.g.blk = victim_line->id;
    printf("[FEMU Dbg] <do_gc>: GC-ing line = %d,victim line's invalid page count = %d,victim line count = %d,full line count = %d,free line count = %d\n", ppa.g.blk,
              victim_line->ipc, ssd->lm.victim_line_cnt, ssd->lm.full_line_cnt,
              ssd->lm.free_line_cnt);

    /* copy back valid data of all the blocks within the victim line. */
    for (ch = 0; ch < spp->nchs; ch++) {
        for (lun = 0; lun < spp->luns_per_ch; lun++) {
            ppa.g.ch = ch;
            ppa.g.lun = lun;
            ppa.g.pl = 0;
            blk_index = get_blk_index(ssd, &ppa);
            // before conducting GC, check whether the block has worn out before. If so, then skip this block.
            if (check_block_status(ssd, blk_index) == BAD_BLOCK){
                continue;
            }

            lunp = get_lun(ssd, &ppa);
            clean_one_block(ssd, &ppa);   
            /*
                wear out simulation before erase operation.
            */
            if (ENABLE_WEAR_OUT_SIMULATION) {
                wear_out_simulator(ssd, &ppa);
                blk_index = get_blk_index(ssd, &ppa);
                if (check_block_status(ssd, blk_index) == BAD_BLOCK){
                    update_block_status(ssd, &ppa, BAD_BLOCK);
                    break;      // abandon the current erase opertion and move to next block.
                }
            }
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
    struct nand_block *target_blk = NULL;
    struct ssdparams *spp = &ssd->sp;   
    uint64_t lba = req->slba;             
    int nsecs = req->nlb;             
    struct ppa ppa;             
    uint64_t start_lpn = lba / spp->secs_per_pg;     
    uint64_t end_lpn = (lba + nsecs - 1) / spp->secs_per_pg;
    uint64_t lpn;
    uint64_t sublat, maxlat = 0;

    if (end_lpn >= spp->tt_pgs) {
        ftl_err("start_lpn=%"PRIu64",tt_pgs=%d\n", start_lpn, ssd->sp.tt_pgs);
    }

    /* normal IO read path */
    for (lpn = start_lpn; lpn <= end_lpn; lpn++) {
        ppa = get_maptbl_ent(ssd, lpn);       
        if (!mapped_ppa(&ppa) || !valid_ppa(ssd, &ppa)) {
            //printf("%s,lpn(%" PRId64 ") not mapped to valid ppa\n", ssd->ssdname, lpn);
            //printf("Invalid ppa,ch:%d,lun:%d,blk:%d,pl:%d,pg:%d,sec:%d\n",
            //ppa.g.ch, ppa.g.lun, ppa.g.blk, ppa.g.pl, ppa.g.pg, ppa.g.sec);
            continue;
        }
        uint64_t wear_out_lat = 0;
        if (ENABLE_WEAR_OUT_SIMULATION) {
            /* Simulation of block wearing */
            wear_out_simulator(ssd, &ppa);
            target_blk = get_blk(ssd, &ppa);
            if (target_blk->block_wear_status == BAD_BLOCK){
                /* Conduct Bad block management. */
                wear_out_lat = bad_block_management(ssd, &ppa, req->stime);
                ppa = get_maptbl_ent(ssd, lpn);    /* get the ppa after rellocation of bad block's data. */
                req->stime += wear_out_lat;   /* re-execute the write request after the bad-block-management is completed. */
                // printf("[FEMU Dbg] <wear_out_simulator> Wear out happens during read operation. Extra latency is: %"PRIu64". \n", wear_out_lat);
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
        maxlat = (sublat > maxlat) ? sublat : maxlat;       // it is calculated to be max of each page, maybe because it is read parallely.
    }

    return maxlat;
}


static uint64_t ssd_write(struct ssd *ssd, NvmeRequest *req)
{
    struct nand_block *target_blk = NULL;
    uint64_t lba = req->slba;      
    struct ssdparams *spp = &ssd->sp;  
    int len = req->nlb;  
    
    uint64_t start_lpn = lba / spp->secs_per_pg;     // 4. get the starting logical page and ending logical page.
    uint64_t end_lpn = (lba + len - 1) / spp->secs_per_pg;
    struct ppa ppa;
    uint64_t lpn;
    uint64_t curlat = 0, maxlat = 0;
    int r;

    if (end_lpn >= spp->tt_pgs) {
        ftl_err("start_lpn=%"PRIu64",tt_pgs=%d\n", start_lpn, ssd->sp.tt_pgs);
    }

    while (should_gc_high(ssd)) {
        /* perform GC here until !should_gc(ssd) */
        r = do_gc(ssd, true);
        if (r == -1)
            break;
    }

    for (lpn = start_lpn; lpn <= end_lpn; lpn++) {
        ppa = get_maptbl_ent(ssd, lpn);     
        if (mapped_ppa(&ppa)) {   
            /* update old page information first */
            mark_page_invalid(ssd, &ppa);
            set_rmap_ent(ssd, INVALID_LPN, &ppa);
        }

        /* new write */  /* block wear out added */
        ppa = get_new_page(ssd, NAND_WRITE);
        uint64_t wear_out_lat = 0;
        if (ENABLE_WEAR_OUT_SIMULATION) {
            /* Simulation of block wearing */
            wear_out_simulator(ssd, &ppa);
            target_blk = get_blk(ssd, &ppa);
            if (target_blk->block_wear_status == BAD_BLOCK){
                wear_out_lat = bad_block_management(ssd, &ppa, req->stime);
                ppa = get_new_page(ssd, OVER_PROVISIONING_WRITE);    /* After conducting bad block management, update the ppa to a healthy block. */
                req->stime += wear_out_lat;   /* re-execute the write request after the bad-block-management is completed. */
                // printf("[FEMU Dbg] <wear_out_simulator> Wear out happens during write operation. Extra latency is: %"PRIu64". \n", wear_out_lat);
            }
        }

        /* update maptbl */
        set_maptbl_ent(ssd, lpn, &ppa);
        /* update rmap */
        set_rmap_ent(ssd, lpn, &ppa);
        // 6. bind logical page address to the newly allocated page.

        mark_page_valid(ssd, &ppa);

        ssd_wear_leveling_handler(ssd, NAND_WRITE, NULL);   // point the write pointer to next (healthy) block

        struct nand_cmd swr;
        swr.type = USER_IO;
        swr.cmd = NAND_WRITE;
        swr.stime = req->stime;
        if (ENABLE_WEAR_OUT_SIMULATION) {
            curlat = ssd_advance_status(ssd, &ppa, &swr) + wear_out_lat;    // emulate the WRITE latency + (possible) WEAR OUT latency.
        } else {
            curlat = ssd_advance_status(ssd, &ppa, &swr);
        }
        maxlat = (curlat > maxlat) ? curlat : maxlat;
    }
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
                lat = ssd_write(ssd, req);
                break;
            case NVME_CMD_READ:
                lat = ssd_read(ssd, req);
                break;
            case NVME_CMD_DSM:
                lat = 0;
                break;
            default:
                femu_debug("[ftl_thread] do nothing. \n");
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

